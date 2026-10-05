/******************************************************************************

This computer code accompanies the paper:

   " Exponential neighborhood search for Consécutive Block Minimization "

by: Salim Haddadi

Submitted to:  International Transactions in Operational Research

May 15, 2021

Note: This code uses the solver "Linkern" which is freely downloadable from
the website:

      http://www.math.uwaterloo.ca/tsp/concorde.html

The executable of the solver should be put in a place where it can be accessed

Computing times are provided by the linux command:

$ time -p ./exec

where exec the executable code of this computer code

The running time of Linkern is given by an internal procedure

Results:

- The best configuration is in binary matrix a_best[][]
- The number of 1-blocks in a_best[][] is numbcos_best

******************************************************************************/

#define _POSIX_C_SOURCE 199309L

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <signal.h>       /* SIGALRM                  */
#include <time.h>         /* clock_gettime, CLOCK_MONOTONIC */
#include <unistd.h>       /* alarm()                  */

#define INFINI  999999999
#define NB_ITER 500
#define RATIO   5
#define M       -1000

/* ── Wall-clock budget (seconds) ─────────────────────────────────────────── */
#define MAX_SECONDS 7200   /* 2 hours — change here if needed */

/* ── Solver paths ─────────────────────────────────────────────────────────── */
#define LKH_PATH      "/home/pedro/cbm/src/LKH3/LKH"
#define LINKERN_PATH  "/home/pedro/cbm/linkern"

/* ── TSP instance directory ───────────────────────────────────────────────── */
#define TSP_DIR       "/home/pedro/cbm/instances/tsp"
#define TSP_FILE      TSP_DIR "/tsp"
#define PAR_FILE      TSP_DIR "/tsp.par"
#define SOL_FILE      TSP_DIR "/sol"
#define TSPSOL_FILE   TSP_DIR "/tsp.sol"
#define AUX_FILE      TSP_DIR "/aux"

/* ── Solver selection ─────────────────────────────────────────────────────── */
typedef enum { SOLVER_LINKERN, SOLVER_LKH } Solver;
Solver solver;

unsigned int m, n, numbcos, numbcos_best;
char *input_file  = NULL;
char *output_file = NULL;
char **a, **a_best;
int iter, nbcol, i, j;
char **b;
int *pi_best;

/* ── Timing state ─────────────────────────────────────────────────────────── */
static struct timespec t_start;

/* Return elapsed wall-clock seconds since t_start. */
static double elapsed_seconds(void)
{
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (double)(now.tv_sec  - t_start.tv_sec)
         + (double)(now.tv_nsec - t_start.tv_nsec) * 1e-9;
}

void read_data(void);
void fill_file_in_TSPLIB_format(void);
void ens(void);
void Compute_numbcos_best(void);

/* ── Write the best result to the output file (or stdout if none given) ───── */
static void write_best(void)
{
    FILE *fp = NULL;

    if (output_file != NULL)
    {
        fp = fopen(output_file, "w");
        if (!fp)
        {
            perror("Cannot open output file");
            /* Fall back to stdout so we don't lose the result. */
        }
    }

    if (!fp) fp = stdout;

    fprintf(fp, "Best number of bco's %d\n", numbcos_best);
    fprintf(fp, "Best column permutation (1-indexed):");
    for (int j = 1; j <= (int)n; j++)
        fprintf(fp, " %d", pi_best[j]);
    fprintf(fp, "\n");

    if (fp != stdout) fclose(fp);
}

/* ---------------------------------------------------------------------------
 * LKH parameter-file helpers
 * ------------------------------------------------------------------------- */
static void write_lkh_par(const char *tsp_file, const char *tour_file)
{
    FILE *fp = fopen(PAR_FILE, "w");
    if (!fp) { perror("Cannot write LKH .par file"); exit(1); }
    fprintf(fp, "PROBLEM_FILE = %s\n", tsp_file);
    fprintf(fp, "TOUR_FILE    = %s\n", tour_file);
    fprintf(fp, "RUNS         = 1\n");
    fprintf(fp, "TIME_LIMIT   = 7200\n");
    fclose(fp);
}

static void run_solver_initial(void)
{
    char cmd[512];
    if (solver == SOLVER_LKH)
    {
        write_lkh_par(TSP_FILE, SOL_FILE);
        snprintf(cmd, sizeof(cmd), "%s %s", LKH_PATH, PAR_FILE);
    }
    else
    {
        snprintf(cmd, sizeof(cmd), "%s -Q -t 7200 -o %s %s",
                 LINKERN_PATH, SOL_FILE, TSP_FILE);
    }
    system(cmd);
}

static void run_solver_ens(void)
{
    char cmd[512];
    if (solver == SOLVER_LKH)
    {
        write_lkh_par(TSP_FILE, TSPSOL_FILE);
        snprintf(cmd, sizeof(cmd), "%s %s", LKH_PATH, PAR_FILE);
    }
    else
    {
        snprintf(cmd, sizeof(cmd), "%s -Q -o %s %s >> %s",
                 LINKERN_PATH, TSPSOL_FILE, TSP_FILE, AUX_FILE);
    }
    system(cmd);
}

/* ---------------------------------------------------------------------------
 * Tour readers
 * ------------------------------------------------------------------------- */
static void read_tour_linkern(const char *path, int *pi, int expected_nodes)
{
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "Cannot open Linkern sol: %s\n", path); exit(1); }

    int n_nodes, tour_len;
    fscanf(fp, "%d %d", &n_nodes, &tour_len);

    if (n_nodes != expected_nodes)
    {
        fprintf(stderr,
                "Linkern tour has %d nodes but expected %d\n",
                n_nodes, expected_nodes);
        exit(1);
    }

    for (int i = 0; i < n_nodes; i++)
    {
        int a_node, b_node, ew;
        fscanf(fp, "%d %d %d", &a_node, &b_node, &ew);
        pi[i] = a_node;
    }
    fclose(fp);
}

static void read_tour_lkh(const char *path, int *pi, int expected_nodes)
{
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "Cannot open LKH tour: %s\n", path); exit(1); }

    char line[256];
    int in_section = 0;
    int count = 0;

    while (fgets(line, sizeof(line), fp))
    {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;

        if (!in_section)
        {
            if (strncmp(p, "TOUR_SECTION", 12) == 0) in_section = 1;
            continue;
        }

        int node;
        if (sscanf(p, "%d", &node) != 1) continue;
        if (node == -1) break;

        if (count >= expected_nodes)
        {
            fprintf(stderr, "LKH tour has more nodes than expected (%d)\n",
                    expected_nodes);
            exit(1);
        }
        pi[count++] = node - 1;
    }
    fclose(fp);

    if (count != expected_nodes)
    {
        fprintf(stderr,
                "LKH tour has %d nodes but expected %d\n",
                count, expected_nodes);
        exit(1);
    }
}

static void read_tour(const char *path, int *pi, int expected_nodes)
{
    if (solver == SOLVER_LKH)
        read_tour_lkh(path, pi, expected_nodes);
    else
        read_tour_linkern(path, pi, expected_nodes);
}

int main(int argc, char *argv[])
{
    /* ── Start the wall clock ─────────────────────────────────────────────── */
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    for (int i = 1; i < argc; i++)
    {
        if (strncmp(argv[i], "--filePath=", 11) == 0)
        {
            input_file = argv[i] + 11;
        }
        else if (strncmp(argv[i], "--outputPath=", 13) == 0)
        {
            output_file = argv[i] + 13;
        }
        else if (strncmp(argv[i], "--algorithm=", 12) == 0)
        {
            char *alg = argv[i] + 12;
            if (strcmp(alg, "lkh") == 0)
                solver = SOLVER_LKH;
            else if (strcmp(alg, "linkern") == 0)
                solver = SOLVER_LINKERN;
            else
            {
                fprintf(stderr, "Unknown algorithm '%s'. Use 'linkern' or 'lkh'.\n", alg);
                exit(1);
            }
        }
        else
        {
            fprintf(stderr, "Unknown argument: %s\n", argv[i]);
            exit(1);
        }
    }

    if (input_file == NULL)
    {
        fprintf(stderr, "Missing required argument --filePath=<path>\n");
        exit(1);
    }

    srand(time(0));

    read_data();

    fill_file_in_TSPLIB_format();
    run_solver_initial();
    Compute_numbcos_best();

    numbcos = numbcos_best;
    for (i = 0; i < m; i++)
    for (j = 0; j <= n; j++) a[i][j] = a_best[i][j];

    for (iter = 1; iter <= NB_ITER; iter++)
    {
        /* ── Soft time-limit check: finish this iteration's setup then stop ── */
        double elapsed = elapsed_seconds();
        if (elapsed >= (double)MAX_SECONDS)
        {
            fprintf(stderr,
                    "\n[TIMEOUT] %.1f s elapsed — stopping after %d/%d iterations.\n",
                    elapsed, iter - 1, NB_ITER);
            break;
        }

        ens();
    }

    write_best();
    return 0;
}

/* =========================================================================
 * ens() — unchanged logic, only local variable declarations updated for C99
 * ========================================================================= */
void ens(void)
{
    FILE *ff;
    int i, j, k, e, f, p, kk;
    long *numbers, max, nnbcol, nbsub;
    int *col, *perm, **small, *limits, **lenght;
    int *chosen;

    nbcol = n / RATIO;
    numbcos = numbcos_best;
    for (i = 0; i < m; i++)
    for (j = 0; j <= n; j++) a[i][j] = a_best[i][j];
    nnbcol = 2 * nbcol + 1;

    numbers = (long *) malloc((n + 1) * sizeof(long));
    col = (int *) malloc((nnbcol + 1) * sizeof(int));

    for (j = 1; j <= n; j++) numbers[j] = rand();
    for (i = 1; i <= nnbcol; i++)
    {
        max = 0;
        for (j = 1; j <= n; j++) if (max < numbers[j])
        {
            max = numbers[j];
            f = j;
        }
        numbers[f] = 0;
        col[i] = f;
    }

    perm = (int *) malloc((nnbcol + 1) * sizeof(int));
    for (i = 1; i <= nnbcol; i++) perm[i] = i;
    for (i = 1; i <= nnbcol - 1; i++)
    for (j = 1; j <= nnbcol - i; j++) if (col[j] >= col[j + 1])
    {
        f = col[j + 1];
        e = perm[j + 1];
        col[j + 1] = col[j];
        perm[j + 1] = perm[j];
        col[j] = f;
        perm[j] = e;
    }
    col[1] = 1;
    col[nnbcol] = n;

    for (i = 1; i <= nnbcol; i++)
    {
        f = i / 2;
        if (i - 2 * f == 0) col[i] = 0;
    }
    nbcol = 0;
    for (i = 1; i <= nnbcol; i++)
    {
        if (col[i] > 0)
        {
            nbcol++;
            col[nbcol] = col[i];
        }
    }
    nbsub = nbcol - 1;
    p = 2 * nbsub;

    limits = (int *) malloc((p + 1) * sizeof(int));
    small = (int **) malloc(m * sizeof(int *));
    for (i = 0; i < m; i++) small[i] = (int *) malloc(p * sizeof(int));

    for (j = 1; j <= nbsub; j++)
    {
        limits[2 * j - 1] = col[j] + 1;
        limits[2 * j]     = col[j + 1];
    }
    limits[1] = 1;

    for (j = 1; j <= p; j++)
    for (i = 0; i < m; i++)
        small[i][j] = a[i][limits[j]];
    for (i = 0; i < m; i++) small[i][0] = 0;

    lenght = (int **) malloc((p + 1) * sizeof(int *));
    for (i = 0; i <= p; i++) lenght[i] = (int *) malloc((p + 1) * sizeof(int));

    for (i = 0; i <= p - 1; i++)
    for (j = i + 1; j <= p; j++)
    {
        e = 0;
        for (k = 0; k < m; k++)
            e += (1 - small[k][i]) * small[k][j] + (1 - small[k][j]) * small[k][i];
        lenght[i][j] = e;
        lenght[j][i] = e;
    }
    for (i = 0; i <= nbsub; i++) lenght[i][i] = 0;
    for (j = 1; j <= nbsub; j++)
    {
        lenght[2 * j - 1][2 * j] = M;
        lenght[2 * j][2 * j - 1] = M;
    }

    if ((ff = fopen(TSP_FILE, "w")) == NULL)
    {
        puts("erreur d'ouverture de fichier");
        exit(1);
    }
    fprintf(ff, "NAME: tsp\n");
    fprintf(ff, "TYPE: TSP\n");
    fprintf(ff, "DIMENSION: %d\n", p + 1);
    fprintf(ff, "EDGE_WEIGHT_TYPE: EXPLICIT\n");
    fprintf(ff, "EDGE_WEIGHT_FORMAT: UPPER_ROW\n");
    fprintf(ff, "EDGE_WEIGHT_SECTION\n");
    for (i = 0; i <= p - 1; i++)
    {
        for (j = i + 1; j <= p; j++) fprintf(ff, "%d ", lenght[i][j]);
        fprintf(ff, "\n");
    }
    fprintf(ff, "EOF");
    fclose(ff);

    run_solver_ens();

    chosen = (int *) malloc((p + 1) * sizeof(int));
    read_tour(TSPSOL_FILE, chosen, p + 1);

    int *src = (int *) malloc((n + 1) * sizeof(int));
    e = 0;
    for (j = 1; j <= nbsub; j++)
    {
        k  = chosen[2 * j - 1];
        kk = chosen[2 * j];

        if (k < kk)
        {
            for (f = limits[k]; f <= limits[kk]; f++)
            {
                e++;
                src[e] = pi_best[f];
                for (i = 0; i < m; i++) b[i][e] = a[i][f];
            }
        }
        else if (k > kk)
        {
            for (f = limits[k]; f >= limits[kk]; f--)
            {
                e++;
                src[e] = pi_best[f];
                for (i = 0; i < m; i++) b[i][e] = a[i][f];
            }
        }
    }
    for (i = 0; i < m; i++)
    for (j = 1; j <= n; j++) a[i][j] = b[i][j];

    numbcos = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 1; j <= n - 1; j++)
            if ((a[i][j] == 1) && (a[i][j + 1] == 0)) numbcos++;
        if (a[i][n] == 1) numbcos++;
    }
    if (numbcos_best > numbcos)
    {
        numbcos_best = numbcos;
        for (i = 0; i < m; i++)
        for (j = 1; j <= n; j++) a_best[i][j] = a[i][j];
        for (j = 1; j <= n; j++) pi_best[j] = src[j];
    }

    free(src);
    free(numbers);
    free(col);
    free(perm);
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
    static int e;
    unsigned int k;
    register unsigned int i, j;
    int *th, *h;

    if ((f = fopen(input_file, "r")) == NULL)
    {
        fprintf(stderr, "Cannot open file %s\n", input_file);
        exit(1);
    }

    fscanf(f, "%d%d", &m, &n);

    a_best = (char **) malloc(m * sizeof(char *));
    for (i = 0; i < m; i++) a_best[i] = (char *) malloc((n + 1) * sizeof(char));

    a = (char **) malloc(m * sizeof(char *));
    for (i = 0; i < m; i++) a[i] = (char *) malloc((n + 1) * sizeof(char));

    b = (char **) malloc(m * sizeof(char *));
    for (i = 0; i < m; i++) b[i] = (char *) malloc((n + 1) * sizeof(char));

    th = (int *) malloc((m + 1) * sizeof(int));

    k = m * n * 0.2;
    h = (int *) malloc((m * n + 1) * sizeof(int));

    k = 0;
    for (i = 1; i <= m; i++)
    {
        fscanf(f, "%d", &e);
        for (j = k + 1; j <= k + e; j++) fscanf(f, "%d", &h[j]);
        k += e;
        th[i] = k;
    }

    fclose(f);

    th[0] = 0;
    for (i = 1; i <= m; i++)
    {
        if (th[i] - th[i - 1] == 0)
        {
            puts("null row");
            printf("i= %d ", i);
            exit(1);
        }
    }

    for (i = 0; i < m; i++)
    for (j = 0; j <= n; j++) a[i][j] = 0;

    for (i = 0; i < m; i++)
    for (j = th[i] + 1; j <= th[i + 1]; j++) a[i][h[j]] = 1;

    free(h);
    free(th);
}

void fill_file_in_TSPLIB_format(void)
{
    register unsigned int i, j;
    FILE *f;
    int e, k;
    int **length;

    length = (int **) malloc((n + 2) * sizeof(int *));
    for (i = 0; i <= n + 1; i++) length[i] = (int *) malloc((n + 2) * sizeof(int));

    for (i = 0; i <= n; i++)
    for (j = i + 1; j <= n; j++)
    {
        e = 0;
        for (k = 0; k < m; k++)
            e += (1 - a[k][i]) * a[k][j] + a[k][i] * (1 - a[k][j]);
        length[i][j] = e;
        length[j][i] = e;
    }
    for (i = 0; i <= n; i++) length[i][i] = 0;

    if ((f = fopen(TSP_FILE, "w")) == NULL)
    {
        puts("erreur d'ouverture du fichier");
        exit(1);
    }
    fprintf(f, "NAME: tsp\n");
    fprintf(f, "TYPE: TSP\n");
    fprintf(f, "COMMENT: \n");
    fprintf(f, "DIMENSION: %d\n", n + 1);
    fprintf(f, "EDGE_WEIGHT_TYPE: EXPLICIT\n");
    fprintf(f, "EDGE_WEIGHT_FORMAT: FULL_MATRIX\n");
    fprintf(f, "EDGE_WEIGHT_SECTION\n");
    for (i = 0; i <= n; i++)
    {
        for (j = 0; j <= n; j++) fprintf(f, " %5d", length[i][j]);
        fprintf(f, "\n");
    }
    fprintf(f, "EOF\n");
    fclose(f);

    for (i = 0; i <= n + 1; i++) free(length[i]);
    free(length);
}

void Compute_numbcos_best(void)
{
    int *pi = (int *) malloc((n + 1) * sizeof(int));
    pi_best  = (int *) malloc((n + 1) * sizeof(int));

    read_tour(SOL_FILE, pi, n + 1);

    for (i = 0; i < m; i++)
    for (j = 0; j <= n; j++) a_best[i][j] = a[i][pi[j]];

    for (j = 0; j <= n; j++) pi_best[j] = pi[j];

    numbcos_best = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 1; j <= n - 1; j++)
            if (a_best[i][j] == 1 && a_best[i][j + 1] == 0) numbcos_best++;
        if (a_best[i][n] == 1) numbcos_best++;
    }
    printf("\nInitial number of bcos : %d", numbcos_best);
    free(pi);
}