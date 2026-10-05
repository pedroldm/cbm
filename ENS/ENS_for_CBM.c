/******************************************************************************

This computer code accompanies the paper:

   " Exponential neighborhood search for Consécutive Block Minimization "

by: Salim Haddadi

Submitted to:  International Transactions in Operational Research

May 15, 2021

Note: This code uses the solver "Linkern" (http://www.math.uwaterloo.ca/tsp/concorde.html)
or, alternatively, LKH.

Experiment-harness changes (the search itself is unchanged):
- every random decision derives from --seed: the column sampling uses a
  portable PRNG and TSP call k gets the seed cbm_derive_seed(seed, k);
- all scratch files live in a private --workDir, so concurrent runs cannot
  collide; solvers are spawned without a shell and their exit status checked;
- TSP tours are rotated to the dummy node 0 (Linkern starts its tour anywhere)
  and ENS tours that break a forced sub-matrix edge are rejected and counted;
- the result is re-validated against the input matrix and written as JSON.

Usage:
  ENS --filePath=<instance> --algorithm=<linkern|lkh> --solverPath=<binary>
      [--outputPath=<json>] [--seed=<n>] [--timeLimit=<s>] [--iterations=<n>]
      [--workDir=<dir>]

******************************************************************************/

#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "../src/common/cbm_lkh_params.h"
#include "../src/common/cbm_seed.h"
#include "../src/common/cbm_spawn.h"

#define RATIO 5
#define M     -1000

typedef enum { SOLVER_LINKERN, SOLVER_LKH } Solver;

/* ── Options ──────────────────────────────────────────────────────────────── */
static Solver solver = SOLVER_LINKERN;
static const char *input_file  = NULL;
static const char *output_file = NULL;
static const char *solver_path = NULL;
static const char *work_dir_opt = NULL;
static uint64_t seed = 1;
static double time_limit = 7200.0;
static int nb_iter = 500;

/* ── Search state (as in the original code) ──────────────────────────────── */
unsigned int m, n, numbcos, numbcos_best;
char **a, **a_best, **b, **orig;
int iter, nbcol, i, j;
int *pi_best;
static cbm_rng rng;

/* ── Instrumentation ──────────────────────────────────────────────────────── */
static struct timespec t_start;
static char work_dir[4096];
static int own_work_dir = 0;
static char tsp_file[4200], par_file[4200], tour_file[4200], log_file[4200];

static long tsp_calls = 0, tsp_limit_hits = 0, invalid_tours = 0;
static double tsp_seconds = 0.0, time_to_best = 0.0;
static uint32_t *call_seeds;
static unsigned int initial_blocks;

typedef struct { int iteration; unsigned int blocks; double elapsed; } HistoryEntry;
static HistoryEntry *history;
static int history_len = 0;

void read_data(void);
void fill_file_in_TSPLIB_format(void);
void ens(void);
void Compute_numbcos_best(void);

static double elapsed_seconds(void)
{
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (double)(now.tv_sec - t_start.tv_sec) + (double)(now.tv_nsec - t_start.tv_nsec) * 1e-9;
}

static void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "ENS error: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(1);
}

static void *xmalloc(size_t size)
{
    void *p = malloc(size);
    if (!p) die("out of memory allocating %zu bytes", size);
    return p;
}

/* ---------------------------------------------------------------------------
 * External TSP solver
 * ------------------------------------------------------------------------- */

/* Rotate a tour so that it starts at node 0 (the dummy column / depot). */
static void rotate_to_zero(int *tour, int len)
{
    int start = 0, *tmp;
    while (start < len && tour[start] != 0) start++;
    if (start == len) die("TSP tour does not visit node 0");
    if (start == 0) return;
    tmp = xmalloc(len * sizeof(int));
    for (int k = 0; k < len; k++) tmp[k] = tour[(start + k) % len];
    memcpy(tour, tmp, len * sizeof(int));
    free(tmp);
}

static void read_tour_linkern(const char *path, int *pi, int expected_nodes)
{
    FILE *fp = fopen(path, "r");
    int n_nodes, n_edges;
    if (!fp) die("cannot open Linkern tour %s", path);
    if (fscanf(fp, "%d %d", &n_nodes, &n_edges) != 2 || n_nodes != expected_nodes)
        die("Linkern tour %s: expected %d nodes", path, expected_nodes);
    for (int k = 0; k < n_nodes; k++) {
        int a_node, b_node, w;
        if (fscanf(fp, "%d %d %d", &a_node, &b_node, &w) != 3) die("Linkern tour %s is truncated", path);
        pi[k] = a_node;
    }
    fclose(fp);
}

static void read_tour_lkh(const char *path, int *pi, int expected_nodes)
{
    FILE *fp = fopen(path, "r");
    char line[256];
    int in_section = 0, count = 0;
    if (!fp) die("cannot open LKH tour %s", path);
    while (fgets(line, sizeof(line), fp)) {
        char *p = line;
        int node;
        while (*p == ' ' || *p == '\t') p++;
        if (!in_section) {
            if (strncmp(p, "TOUR_SECTION", 12) == 0) in_section = 1;
            continue;
        }
        if (sscanf(p, "%d", &node) != 1) continue;
        if (node == -1) break;
        if (count >= expected_nodes) die("LKH tour %s has more than %d nodes", path, expected_nodes);
        pi[count++] = node - 1;
    }
    fclose(fp);
    if (count != expected_nodes) die("LKH tour %s has %d nodes, expected %d", path, count, expected_nodes);
}

/* Did the solver stop on its time bound? LKH says so explicitly (TRACE_LEVEL
 * >= 1, its default); Linkern does not, so we fall back to the elapsed time. */
static int hit_time_limit(double seconds, int limit)
{
    if (solver == SOLVER_LKH) {
        FILE *fp = fopen(log_file, "r");
        char line[512];
        int hit = 0;
        if (!fp) return 0;
        while (!hit && fgets(line, sizeof(line), fp))
            if (strstr(line, "Time limit exceeded")) hit = 1;
        fclose(fp);
        return hit;
    }
    return seconds >= 0.99 * limit;
}

/* Solve TSP_FILE with the configured solver; `call` is the logical call index
 * (0 = initial tour, k = ENS iteration k) and fixes the solver seed. */
static void solve_tsp(int call, int *tour, int nodes)
{
    char seed_str[32], limit_str[32];
    double remaining = time_limit - elapsed_seconds(), start;
    int limit = remaining > 1.0 ? (int)ceil(remaining) : 1;
    uint32_t call_seed = cbm_derive_seed(seed, (uint64_t)call);
    int rc;

    call_seeds[call] = call_seed;
    snprintf(seed_str, sizeof(seed_str), "%u", call_seed);
    snprintf(limit_str, sizeof(limit_str), "%d", limit);
    remove(tour_file);
    remove(log_file);

    if (solver == SOLVER_LKH) {
        FILE *fp = fopen(par_file, "w");
        char *argv[] = {(char *)solver_path, par_file, NULL};
        if (!fp) die("cannot write %s: %s", par_file, strerror(errno));
        fprintf(fp, "PROBLEM_FILE = %s\nTOUR_FILE = %s\nMOVE_TYPE = %d\nPATCHING_C = %d\nPATCHING_A = %d\nRUNS = 1\nTIME_LIMIT = %d\nSEED = %u\n",
                tsp_file, tour_file, CBM_LKH_MOVE_TYPE, CBM_LKH_PATCHING_C, CBM_LKH_PATCHING_A, limit, call_seed);
        fclose(fp);
        start = elapsed_seconds();
        rc = cbm_run(argv, log_file);
    } else {
        char *argv[] = {(char *)solver_path, "-Q", "-s", seed_str, "-t", limit_str, "-o", tour_file, tsp_file, NULL};
        start = elapsed_seconds();
        rc = cbm_run(argv, log_file);
    }

    double seconds = elapsed_seconds() - start;
    tsp_calls++;
    tsp_seconds += seconds;
    if (rc != 0) die("TSP solver %s failed on call %d (status %d); see %s", solver_path, call, rc, log_file);
    if (hit_time_limit(seconds, limit)) tsp_limit_hits++;

    if (solver == SOLVER_LKH)
        read_tour_lkh(tour_file, tour, nodes);
    else
        read_tour_linkern(tour_file, tour, nodes);
    rotate_to_zero(tour, nodes);
}

/* ---------------------------------------------------------------------------
 * Validation and output
 * ------------------------------------------------------------------------- */
static unsigned int count_blocks_original(const int *perm)
{
    unsigned int blocks = 0;
    for (unsigned int r = 0; r < m; r++) {
        int prev = 0;
        for (unsigned int k = 1; k <= n; k++) {
            int cur = orig[r][perm[k]];
            if (cur && !prev) blocks++;
            prev = cur;
        }
    }
    return blocks;
}

static void validate_best(void)
{
    char *seen = calloc(n + 1, 1);
    if (!seen) die("out of memory");
    for (unsigned int k = 1; k <= n; k++) {
        if (pi_best[k] < 1 || pi_best[k] > (int)n || seen[pi_best[k]]) die("best solution is not a permutation");
        seen[pi_best[k]] = 1;
    }
    free(seen);
    unsigned int recount = count_blocks_original(pi_best);
    if (recount != numbcos_best) die("block count mismatch: tracked %u, recomputed %u", numbcos_best, recount);
}

static void write_json(FILE *fp, int iterations_done, const char *stop_reason)
{
    fprintf(fp, "{\n  \"method\": \"ENS\",\n  \"tsp_backend\": \"%s\",\n", solver == SOLVER_LKH ? "lkh" : "linkern");
    fprintf(fp, "  \"instance\": \"%s\",\n  \"rows\": %u,\n  \"cols\": %u,\n", input_file, m, n);
    fprintf(fp, "  \"seed\": %llu,\n", (unsigned long long)seed);
    fprintf(fp, "  \"parameters\": {\"iterations\": %d, \"ratio\": %d, \"forced_edge_weight\": %d, \"time_limit_s\": %.3f},\n",
            nb_iter, RATIO, M, time_limit);
    fprintf(fp, "  \"initial_blocks\": %u,\n  \"best_blocks\": %u,\n  \"validated\": true,\n", initial_blocks, numbcos_best);
    fprintf(fp, "  \"iterations_completed\": %d,\n  \"stop_reason\": \"%s\",\n", iterations_done, stop_reason);
    fprintf(fp, "  \"time_limit_reached\": %s,\n", strcmp(stop_reason, "time_limit") == 0 ? "true" : "false");
    fprintf(fp, "  \"elapsed_s\": %.6f,\n  \"time_to_best_s\": %.6f,\n", elapsed_seconds(), time_to_best);
    fprintf(fp, "  \"tsp\": {\"calls\": %ld, \"time_s\": %.6f, \"time_limit_hits\": %ld, \"invalid_tours\": %ld, \"seeds\": [",
            tsp_calls, tsp_seconds, tsp_limit_hits, invalid_tours);
    for (long k = 0; k < tsp_calls; k++) fprintf(fp, "%s%u", k ? ", " : "", call_seeds[k]);
    fprintf(fp, "]},\n  \"history\": [");
    for (int k = 0; k < history_len; k++)
        fprintf(fp, "%s{\"iteration\": %d, \"blocks\": %u, \"elapsed_s\": %.6f}", k ? ", " : "", history[k].iteration,
                history[k].blocks, history[k].elapsed);
    fprintf(fp, "],\n  \"permutation\": [");
    for (unsigned int k = 1; k <= n; k++) fprintf(fp, "%s%d", k > 1 ? ", " : "", pi_best[k]);
    fprintf(fp, "]\n}\n");
}

/* Atomic: a reader never sees a half-written result file. */
static void write_result(int iterations_done, const char *stop_reason)
{
    if (!output_file) {
        write_json(stdout, iterations_done, stop_reason);
        return;
    }
    char tmp[4200];
    snprintf(tmp, sizeof(tmp), "%s.tmp", output_file);
    FILE *fp = fopen(tmp, "w");
    if (!fp) die("cannot open %s: %s", tmp, strerror(errno));
    write_json(fp, iterations_done, stop_reason);
    if (fflush(fp) != 0 || fsync(fileno(fp)) != 0 || fclose(fp) != 0) die("cannot write %s", tmp);
    if (rename(tmp, output_file) != 0) die("cannot rename %s: %s", tmp, strerror(errno));
}

static void remove_work_dir(void)
{
    remove(tsp_file);
    remove(par_file);
    remove(tour_file);
    remove(log_file);
    if (own_work_dir) rmdir(work_dir);
}

static void parse_args(int argc, char *argv[])
{
    for (int k = 1; k < argc; k++) {
        const char *arg = argv[k];
        if (strncmp(arg, "--filePath=", 11) == 0) input_file = arg + 11;
        else if (strncmp(arg, "--outputPath=", 13) == 0) output_file = arg + 13;
        else if (strncmp(arg, "--solverPath=", 13) == 0) solver_path = arg + 13;
        else if (strncmp(arg, "--workDir=", 10) == 0) work_dir_opt = arg + 10;
        else if (strncmp(arg, "--seed=", 7) == 0) seed = strtoull(arg + 7, NULL, 10);
        else if (strncmp(arg, "--timeLimit=", 12) == 0) time_limit = atof(arg + 12);
        else if (strncmp(arg, "--iterations=", 13) == 0) nb_iter = atoi(arg + 13);
        else if (strncmp(arg, "--algorithm=", 12) == 0) {
            if (strcmp(arg + 12, "lkh") == 0) solver = SOLVER_LKH;
            else if (strcmp(arg + 12, "linkern") == 0) solver = SOLVER_LINKERN;
            else die("unknown algorithm '%s' (use linkern or lkh)", arg + 12);
        } else die("unknown argument: %s", arg);
    }
    if (!input_file) die("missing --filePath=<path>");
    if (!solver_path) solver_path = getenv(solver == SOLVER_LKH ? "LKH_PATH" : "LINKERN_PATH");
    if (!solver_path || access(solver_path, X_OK) != 0)
        die("TSP solver not executable: pass --solverPath=<binary> (or set %s)", solver == SOLVER_LKH ? "LKH_PATH" : "LINKERN_PATH");
    if (time_limit <= 0 || nb_iter < 0) die("--timeLimit must be > 0 and --iterations >= 0");
}

static void setup_work_dir(void)
{
    if (work_dir_opt) {
        snprintf(work_dir, sizeof(work_dir), "%s", work_dir_opt);
        if (mkdir(work_dir, 0755) != 0 && errno != EEXIST) die("cannot create %s: %s", work_dir, strerror(errno));
    } else {
        const char *tmp = getenv("TMPDIR");
        snprintf(work_dir, sizeof(work_dir), "%s/ens-XXXXXX", tmp ? tmp : "/tmp");
        if (!mkdtemp(work_dir)) die("cannot create a work directory: %s", strerror(errno));
        own_work_dir = 1;
    }
    snprintf(tsp_file, sizeof(tsp_file), "%s/tsp", work_dir);
    snprintf(par_file, sizeof(par_file), "%s/tsp.par", work_dir);
    snprintf(tour_file, sizeof(tour_file), "%s/tour", work_dir);
    snprintf(log_file, sizeof(log_file), "%s/solver.log", work_dir);
}

int main(int argc, char *argv[])
{
    const char *stop_reason = "iterations";
    int iterations_done = 0;

    clock_gettime(CLOCK_MONOTONIC, &t_start);
    parse_args(argc, argv);
    setup_work_dir();

    rng.state = cbm_splitmix64(seed);
    call_seeds = xmalloc((nb_iter + 1) * sizeof(uint32_t));
    history = xmalloc((nb_iter + 1) * sizeof(HistoryEntry));

    read_data();

    fill_file_in_TSPLIB_format();
    Compute_numbcos_best();
    initial_blocks = numbcos_best;
    time_to_best = elapsed_seconds();
    history[history_len++] = (HistoryEntry){0, numbcos_best, time_to_best};

    numbcos = numbcos_best;
    for (i = 0; i < m; i++)
        for (j = 0; j <= n; j++) a[i][j] = a_best[i][j];

    for (iter = 1; iter <= nb_iter; iter++) {
        /* Soft time limit, checked between iterations as in the original. */
        if (elapsed_seconds() >= time_limit) {
            stop_reason = "time_limit";
            break;
        }
        ens();
        iterations_done = iter;
    }

    validate_best();
    write_result(iterations_done, stop_reason);
    remove_work_dir();
    return 0;
}

/* =========================================================================
 * ens() — one exponential-neighborhood step (search logic as published)
 * ========================================================================= */
void ens(void)
{
    FILE *ff;
    int i, j, k, e, f = 0, p, kk;
    long *numbers, max, nnbcol, nbsub;
    int *col, **small, *limits, **lenght;
    int *chosen;

    nbcol = n / RATIO;
    numbcos = numbcos_best;
    for (i = 0; i < m; i++)
        for (j = 0; j <= n; j++) a[i][j] = a_best[i][j];
    nnbcol = 2 * nbcol + 1;

    numbers = xmalloc((n + 1) * sizeof(long));
    col = xmalloc((nnbcol + 1) * sizeof(int));

    /* Sample nnbcol distinct columns: those with the largest random keys. */
    for (j = 1; j <= n; j++) numbers[j] = (long)(cbm_rng_next(&rng) >> 33);
    for (i = 1; i <= nnbcol; i++) {
        max = 0;
        for (j = 1; j <= n; j++)
            if (max < numbers[j]) {
                max = numbers[j];
                f = j;
            }
        numbers[f] = 0;
        col[i] = f;
    }

    for (i = 1; i <= nnbcol - 1; i++)
        for (j = 1; j <= nnbcol - i; j++)
            if (col[j] >= col[j + 1]) {
                f = col[j + 1];
                col[j + 1] = col[j];
                col[j] = f;
            }
    col[1] = 1;
    col[nnbcol] = n;

    for (i = 1; i <= nnbcol; i++) {
        f = i / 2;
        if (i - 2 * f == 0) col[i] = 0;
    }
    nbcol = 0;
    for (i = 1; i <= nnbcol; i++) {
        if (col[i] > 0) {
            nbcol++;
            col[nbcol] = col[i];
        }
    }
    nbsub = nbcol - 1;
    p = 2 * nbsub;

    limits = xmalloc((p + 1) * sizeof(int));
    small = xmalloc(m * sizeof(int *));
    for (i = 0; i < m; i++) small[i] = xmalloc((p + 1) * sizeof(int));

    for (j = 1; j <= nbsub; j++) {
        limits[2 * j - 1] = col[j] + 1;
        limits[2 * j] = col[j + 1];
    }
    limits[1] = 1;

    for (j = 1; j <= p; j++)
        for (i = 0; i < m; i++) small[i][j] = a[i][limits[j]];
    for (i = 0; i < m; i++) small[i][0] = 0;

    lenght = xmalloc((p + 1) * sizeof(int *));
    for (i = 0; i <= p; i++) lenght[i] = xmalloc((p + 1) * sizeof(int));

    for (i = 0; i <= p - 1; i++)
        for (j = i + 1; j <= p; j++) {
            e = 0;
            for (k = 0; k < m; k++) e += (1 - small[k][i]) * small[k][j] + (1 - small[k][j]) * small[k][i];
            lenght[i][j] = e;
            lenght[j][i] = e;
        }
    for (i = 0; i <= nbsub; i++) lenght[i][i] = 0;
    /* Sub-matrix endpoints are tied together by a strongly negative edge. */
    for (j = 1; j <= nbsub; j++) {
        lenght[2 * j - 1][2 * j] = M;
        lenght[2 * j][2 * j - 1] = M;
    }

    if ((ff = fopen(tsp_file, "w")) == NULL) die("cannot write %s: %s", tsp_file, strerror(errno));
    fprintf(ff, "NAME: tsp\n");
    fprintf(ff, "TYPE: TSP\n");
    fprintf(ff, "DIMENSION: %d\n", p + 1);
    fprintf(ff, "EDGE_WEIGHT_TYPE: EXPLICIT\n");
    fprintf(ff, "EDGE_WEIGHT_FORMAT: UPPER_ROW\n");
    fprintf(ff, "EDGE_WEIGHT_SECTION\n");
    for (i = 0; i <= p - 1; i++) {
        for (j = i + 1; j <= p; j++) fprintf(ff, "%d ", lenght[i][j]);
        fprintf(ff, "\n");
    }
    fprintf(ff, "EOF");
    fclose(ff);

    chosen = xmalloc((p + 1) * sizeof(int));
    solve_tsp(iter, chosen, p + 1);

    /* After rotating to node 0, positions (2j-1, 2j) must hold the two
     * endpoints of one sub-matrix. A tour that drops a forced edge cannot be
     * decoded into a permutation, so the step is rejected. */
    int valid = 1;
    for (j = 1; j <= nbsub && valid; j++) {
        int lo = chosen[2 * j - 1] < chosen[2 * j] ? chosen[2 * j - 1] : chosen[2 * j];
        int hi = chosen[2 * j - 1] ^ chosen[2 * j] ^ lo;
        if (lo % 2 != 1 || hi != lo + 1) valid = 0;
    }

    if (valid) {
        int *src = xmalloc((n + 1) * sizeof(int));
        e = 0;
        for (j = 1; j <= nbsub; j++) {
            k = chosen[2 * j - 1];
            kk = chosen[2 * j];

            if (k < kk) {
                for (f = limits[k]; f <= limits[kk]; f++) {
                    e++;
                    src[e] = pi_best[f];
                    for (i = 0; i < m; i++) b[i][e] = a[i][f];
                }
            } else if (k > kk) {
                for (f = limits[k]; f >= limits[kk]; f--) {
                    e++;
                    src[e] = pi_best[f];
                    for (i = 0; i < m; i++) b[i][e] = a[i][f];
                }
            }
        }
        for (i = 0; i < m; i++)
            for (j = 1; j <= n; j++) a[i][j] = b[i][j];

        numbcos = 0;
        for (i = 0; i < m; i++) {
            for (j = 1; j <= n - 1; j++)
                if ((a[i][j] == 1) && (a[i][j + 1] == 0)) numbcos++;
            if (a[i][n] == 1) numbcos++;
        }
        if (numbcos_best > numbcos) {
            numbcos_best = numbcos;
            for (i = 0; i < m; i++)
                for (j = 1; j <= n; j++) a_best[i][j] = a[i][j];
            for (j = 1; j <= n; j++) pi_best[j] = src[j];
            time_to_best = elapsed_seconds();
            history[history_len++] = (HistoryEntry){iter, numbcos_best, time_to_best};
        }
        free(src);
    } else {
        invalid_tours++;
    }

    free(numbers);
    free(col);
    free(limits);
    for (i = 0; i < m; i++) free(small[i]);
    free(small);
    for (i = 0; i <= p; i++) free(lenght[i]);
    free(lenght);
    free(chosen);
}

void read_data(void)
{
    FILE *f;
    int e;
    unsigned int k;
    unsigned int i, j;
    int *th, *h;
    long total = 0;

    if ((f = fopen(input_file, "r")) == NULL) die("cannot open instance %s", input_file);
    if (fscanf(f, "%u%u", &m, &n) != 2 || m == 0 || n == 0) die("malformed instance header in %s", input_file);

    a_best = xmalloc(m * sizeof(char *));
    a = xmalloc(m * sizeof(char *));
    b = xmalloc(m * sizeof(char *));
    orig = xmalloc(m * sizeof(char *));
    for (i = 0; i < m; i++) {
        a_best[i] = xmalloc((n + 1) * sizeof(char));
        a[i] = xmalloc((n + 1) * sizeof(char));
        b[i] = xmalloc((n + 1) * sizeof(char));
        orig[i] = calloc(n + 1, sizeof(char));
        if (!orig[i]) die("out of memory");
    }

    th = xmalloc((m + 1) * sizeof(int));
    h = NULL;

    /* Read row by row (the original sized a buffer of m*n ints up front). */
    k = 0;
    th[0] = 0;
    for (i = 1; i <= m; i++) {
        if (fscanf(f, "%d", &e) != 1 || e < 0) die("malformed row %u in %s", i, input_file);
        h = realloc(h, (total + e + 1) * sizeof(int));
        if (!h) die("out of memory");
        for (j = k + 1; j <= k + e; j++) {
            if (fscanf(f, "%d", &h[j]) != 1 || h[j] < 1 || h[j] > (int)n) die("malformed row %u in %s", i, input_file);
        }
        k += e;
        total += e;
        th[i] = k;
    }

    fclose(f);

    for (i = 1; i <= m; i++)
        if (th[i] - th[i - 1] == 0) die("null row %u in %s", i, input_file);

    for (i = 0; i < m; i++)
        for (j = 0; j <= n; j++) a[i][j] = 0;

    for (i = 0; i < m; i++)
        for (j = th[i] + 1; j <= (unsigned int)th[i + 1]; j++) a[i][h[j]] = 1;
    for (i = 0; i < m; i++) memcpy(orig[i], a[i], n + 1);

    free(h);
    free(th);
}

void fill_file_in_TSPLIB_format(void)
{
    unsigned int i, j;
    FILE *f;
    int e, k;
    int **length;

    length = xmalloc((n + 2) * sizeof(int *));
    for (i = 0; i <= n + 1; i++) length[i] = xmalloc((n + 2) * sizeof(int));

    for (i = 0; i <= n; i++)
        for (j = i + 1; j <= n; j++) {
            e = 0;
            for (k = 0; k < m; k++) e += (1 - a[k][i]) * a[k][j] + a[k][i] * (1 - a[k][j]);
            length[i][j] = e;
            length[j][i] = e;
        }
    for (i = 0; i <= n; i++) length[i][i] = 0;

    if ((f = fopen(tsp_file, "w")) == NULL) die("cannot write %s: %s", tsp_file, strerror(errno));
    fprintf(f, "NAME: tsp\n");
    fprintf(f, "TYPE: TSP\n");
    fprintf(f, "COMMENT: \n");
    fprintf(f, "DIMENSION: %d\n", n + 1);
    fprintf(f, "EDGE_WEIGHT_TYPE: EXPLICIT\n");
    fprintf(f, "EDGE_WEIGHT_FORMAT: FULL_MATRIX\n");
    fprintf(f, "EDGE_WEIGHT_SECTION\n");
    for (i = 0; i <= n; i++) {
        for (j = 0; j <= n; j++) fprintf(f, " %5d", length[i][j]);
        fprintf(f, "\n");
    }
    fprintf(f, "EOF\n");
    if (fclose(f) != 0) die("cannot write %s", tsp_file);

    for (i = 0; i <= n + 1; i++) free(length[i]);
    free(length);
}

void Compute_numbcos_best(void)
{
    int *pi = xmalloc((n + 1) * sizeof(int));
    pi_best = xmalloc((n + 1) * sizeof(int));

    solve_tsp(0, pi, n + 1);

    for (i = 0; i < m; i++)
        for (j = 0; j <= n; j++) a_best[i][j] = a[i][pi[j]];

    for (j = 0; j <= n; j++) pi_best[j] = pi[j];

    numbcos_best = 0;
    for (i = 0; i < m; i++) {
        for (j = 1; j <= n - 1; j++)
            if (a_best[i][j] == 1 && a_best[i][j + 1] == 0) numbcos_best++;
        if (a_best[i][n] == 1) numbcos_best++;
    }
    free(pi);
}
