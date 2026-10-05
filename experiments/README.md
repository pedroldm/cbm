# CBM thesis experiments: runner

`run_experiments.py` runs ENS, ILS, standalone LKH and CBMLKH on every instance in
`instances/`, with reproducible seeds, process isolation, a fixed timeout policy,
crash-safe checkpointing and CSV reports. Python ≥ 3.8, standard library only.

## Quick start

```bash
make                                   # ENS/ENS, ILS/ils, src/CBMLKH/{main_prd,lkh_standalone}, src/LKH3/LKH
export LINKERN_PATH=/path/to/linkern   # or --linkern-path; Linkern is not built here

python3 experiments/run_experiments.py plan            # dry run: jobs, CPUs, largest memory estimate
nohup python3 experiments/run_experiments.py run > outputs/final.nohup 2>&1 &
python3 experiments/run_experiments.py status          # repetitions per method and status
python3 experiments/run_experiments.py report          # regenerate CSVs (also done while running)
```

**Resume** after Ctrl-C, `kill`, a crash, a reboot or a power loss: run the same
`run` command again, or just `run --name final`. Completed repetitions are never
re-run. Interrupted and stale jobs are re-run with the same seed. Failed and
hard-timeout jobs are final outcomes and are kept, unless you pass
`--retry-failed` or `--retry-timeouts`.

**On another machine**, all solver paths are configurable, each with an
environment-variable alternative:

| option | env | default |
|---|---|---|
| `--linkern-path` | `LINKERN_PATH` | `<repo>/bin/linkern` |
| `--lkh-path` | `LKH_PATH` | `<repo>/src/LKH3/LKH` |
| `--ens-bin`, `--ils-bin`, `--cbmlkh-bin`, `--lkh-standalone-bin` | `CBM_ENS_BIN`, … | built locations |

Also useful there: `--cpus 0-15` or `--max-cpus N`, `--mem-budget-gb`, and `--work-root`
(put scratch on a large local disk, not tmpfs). If you move an experiment
directory, pass `--instance-dir` to remap instance paths. Instances are checked by
SHA-256 before use.

## Experimental design (defaults = the full design)

| Method | Repetitions / instance | Process | TSP solver |
|---|---:|---|---|
| ENS (Haddadi 2021; 500 iterations) | 10 | 1 thread | Linkern |
| ILS (NBITER=24, φ=0.2, V=1) | 10 | 1 thread | Linkern |
| Standalone LKH (whole instance as one TSP) | 10 | 1 thread | LKH |
| CBMLKH (irace config #135, `maxIterations` = 1000, `lkhMaxTime` = 20 % of the budget) | 10 | `T` threads per execution | LKH |

Every LKH call (standalone, inside CBMLKH, and ENS/ILS if switched to LKH) uses
`MOVE_TYPE = 5`, `PATCHING_C = 3`, `PATCHING_A = 2`, compiled in from
`src/common/cbm_lkh_params.h`. Linkern runs with its defaults plus a seed (`-s`)
and the remaining time budget (`-t`). ENS and ILS use the parameters of their
original codes.

CBMLKH's parameters are the last completed irace race's best configuration
(id 135, `tunning/output.txt`; checked by `tests/test_tuning.py`), with two
exceptions: each trajectory stops after `maxIterations = 1000` (irace: 1 000 000) or
the time budget, whichever comes first; and each LKH call is limited by
`lkhMaxTime = 0.2 × soft limit` (1440 s for 7200 s) instead of irace's 300 s. They, and the LKH parameters above, are locked:
`--method-params` cannot change them. The testing-only `--allow-param-overrides`
lifts the lock, and its experiment is always marked partial.

A CBMLKH execution with `T` threads produces `T` repetitions (one per trajectory),
so `T = min(10, CPU budget)` and ⌈10 / T⌉ executions per instance. Override with
`--cbmlkh-threads`.

Root seed default: `20261005` (`--seed`). Partial designs (`--methods`,
`--instances`, `--repetitions`, other limits, or `--method-params` overrides) are
flagged `partial: true` and go to `outputs/partial-<hash>` by default. They can never
be named `final`. Resuming with a different design under the same name is refused
(design fingerprint mismatch).

## Timeout policy

* **Soft limit**: 2 h (`--soft-limit`). Each method gets it as its own time budget:
  ENS/ILS `--timeLimit`, CBMLKH `maxTime`, standalone LKH `TIME_LIMIT`. ENS and ILS
  also cap every Linkern/LKH call at the remaining budget. The runner itself only
  records the crossing (`soft_limit_reached` in the job's events and in the record);
  it never kills anything at the soft limit.
* **Hard limit**: 5 h (`--hard-limit`). The runner sends SIGTERM to the job's whole
  process group, then SIGKILL after `--kill-grace` (30 s). The status becomes
  `hard_timeout` and no objective is recorded.
* **Why the soft limit can be exceeded**: budgets are checked between iterations,
  and some phases are not interruptible: ENS/ILS build an O(n²·m) distance matrix,
  LKH's preprocessing is not covered by TIME_LIMIT, and a CBMLKH iteration may run
  one LKH call of up to `lkhMaxTime` (1440 s) plus that call's preprocessing. The same policy applies to standalone
  LKH.

## Reproducibility

**Seeds.** For ENS, ILS and LKH, `run = run_seed(root, method, instance, rep)`
(SHA-256 based, `seeds.py`). CBMLKH uses `run_seed(root, "cbmlkh", instance, -1)`,
and trajectory `r` uses `cbm_derive_seed(run, r)`. Inside the solvers
(`src/common/cbm_seed.h`, mirrored in Python):

* ENS: Linkern call `k` (0 = initial tour, k = iteration k) gets `-s cbm_derive_seed(run, k)`,
  and column sampling uses a portable splitmix64 stream instead of `rand()`/`srand(time)`.
* ILS: TSP call `k` gets `cbm_derive_seed(run, k)`. The perturbation `mt19937` is seeded
  with `cbm_derive_seed(run, 2^64−1)`.
* LKH: `SEED = run`.
* CBMLKH: trajectory `r` gets `cbm_derive_seed(run, r)` whichever OS thread runs it.
  Each LKH sub-problem is solved in **canonical (sorted) column order** with seed
  `cbm_derive_seed(run, hash(columns))`. LKH's answer is therefore a pure function
  of the column set, so the shared cache cannot make one trajectory depend on
  another, and splitting repetitions across executions (`trajectoryOffset`) does
  not change them. The test suite verified that 1×4 threads equals 2×2 threads,
  trajectory for trajectory.

All seeds lie in [1, 2³¹−1]. Linkern treats seed 0 as "use the clock".
Every derived seed is stored in the record (`seeds.derived`). For CBMLKH's LKH
calls the rule is stored instead, since there can be thousands of them.

**What is guaranteed.** Same root seed, instance, method, repetition and
configuration, on the same binaries and platform, give the same solution, **provided
no time bound binds**. Tests check this for all four methods, parallel vs serial.

**Limitations (inherent, documented in each record):**
* A run that stops on its time budget performs a load-dependent amount of work.
  `stopping.method_time_limit_reached` flags it. With a 2 h budget most large
  instances will be time-bound, so their exact replay is not guaranteed.
* An LKH call that hits `TIME_LIMIT` returns a load-dependent tour.
  `tsp_solver.time_limit_hits` counts these: exact for LKH (it prints
  `*** Time limit exceeded ***`), heuristic for Linkern (call time ≥ 99 % of its bound).
* Bit-identical results across machines additionally need the same compiler and
  libstdc++. `std::uniform_*_distribution` and `pow()` are implementation-defined.
  ENS uses its own PRNG and is not affected.
* CBMLKH trajectories of one execution share the LKH cache and the machine. Their
  *results* are seed-determined, but under a time budget a cache hit buys extra
  iterations, so `--cbmlkh-threads` is part of the design fingerprint.

## Concurrency and isolation

Audit (from the sources): ENS and ILS are single-threaded and call one Linkern
process at a time. `lkh_standalone` runs one single-threaded LKH. CBMLKH runs
exactly `threads` OpenMP workers, each running at most one LKH child, so it uses
`threads` CPUs.

* Every job runs in its own session/process group, pinned to its reserved CPUs
  (`sched_setaffinity`, inherited by LKH/Linkern), with `OMP_NUM_THREADS` = its threads.
* Every attempt has a private work directory (also `TMPDIR`), deleted afterwards
  (`*.log` files ≤ 1 MB are kept). Solvers no longer use fixed paths: the old ENS
  wrote `instances/tsp/{tsp,sol,...}`, and CBMLKH emptied a shared `/tmp/LKH/`.
* **CBMLKH exclusion.** At most one CBMLKH job runs at a time:
  * the scheduler never starts a second one;
  * a machine-wide `flock` (`--cbmlkh-lock`, default `$XDG_RUNTIME_DIR/cbm-experiments/cbmlkh.lock`)
    covers other runner processes and other experiments;
  * the lock descriptor is inherited by the CBMLKH process, so it stays held even if
    the runner dies.

  While CBMLKH work remains, `T` CPUs (and its memory estimate) are kept as a lane
  for it, so single-thread jobs cannot starve it.
* **Resource budget.** CPU budget defaults to one logical CPU per physical core. Jobs
  also reserve estimated RAM (budget defaults to 85 % of MemTotal) and scratch disk
  (`methods.estimate_resources`: the explicit distance matrices dominate). A job
  larger than the whole budget runs alone.
* **One runner per experiment** (`runner.lock`). `--allow-concurrent-runners` lifts
  this; per-job `flock`s then still prevent double claims (tested).

Platform: Linux only (flock, `/proc`, `prctl(PR_SET_PDEATHSIG)`, `sched_setaffinity`).
flock is not reliable on some network filesystems, so keep `outputs/` local.

## Fault tolerance

* `runs/<job>/state.json` is the authoritative status. It is rewritten atomically
  (write temp, fsync, rename, fsync dir) on start, soft limit, hard limit and finish.
  Every transition is also appended to its `events` list.
* Records are written before the state turns `completed`. On startup, each
  completed job's records are re-validated against the schema, and a missing or
  damaged record resets the job to `pending`.
* A `running` state left by a dead runner: if no process holds the job lock, the job
  becomes `interrupted` (`stale_running_recovered`). If orphans still hold it, their
  process group is killed first (`orphan_killed`). Solver processes die with the
  runner via `PR_SET_PDEATHSIG`.
* SIGINT/SIGTERM/SIGHUP to the runner terminates running jobs and marks them
  `interrupted`. A second signal escalates to SIGKILL. Exit code 130.

## Output layout and formats

```
outputs/<name>/
  experiment.json        spec, fingerprint, partial flag, immutable job list (seeds, CPUs, estimates)
  manifest.json          snapshot: every job's status/attempts/record paths (refreshed every ≤5 s)
  runs/<job_id>/state.json                     status + events (authoritative)
  runs/<job_id>/<run_id>.json                  record (schema below), one per repetition
  runs/<job_id>/<run_id>.solution.json         best permutation (1-based column ids)
  runs/<job_id>/attempt-NN/                    command.json, stdout.log, stderr.log, native.json,
                                               cbmlkh.cfg, work_logs/, superseded/ (records of a retried attempt)
  results/{detailed,summary,comparison,methods}.csv
  logs/runner.log, logs/sessions.jsonl
```

`job_id` = `<method>__<instance>__rNN` (CBMLKH: `rNN-MM`), and `run_id` = `<method>__<instance>__rNN`.

**Record** (`schema/run_record.schema.json`, version `1.0`). `null` always means
"not available", never zero.
`status` ∈ completed | failed | hard_timeout | interrupted.
`objective.value` = number of 1-blocks (minimize), independently recounted by the
runner from the instance. A mismatch makes the run `failed`. Other fields:
`timing` (wall, method-reported, time-to-best, soft/hard limits, soft-limit crossing),
`stopping`, `search` (iterations, improvements), `tsp_solver` (calls, time,
time-limit hits, rejected tours), `seeds`, `execution` (attempt, exit code/signal,
CPUs, command, logs), `error` (type, message, stderr tail), `software` (runner
version, git commit + dirty flag, Python, platform, SHA-256 of every binary used),
and `method_stats` (compact native statistics; the full native JSON is in
`attempt-NN/native.json`).

**CSVs** are regenerated from the records, idempotently and atomically, after
completions (at most every 10 s) and at exit:
* `detailed.csv`: one row per planned repetition, *including* pending, failed and
  timed-out ones.
* `summary.csv`: per (instance, method). It counts every status. Objective, wall-time
  and time-to-best statistics (mean, median, sample std, min, max) use **only
  `completed` repetitions** (all of which are validated). `hits_best_observed`
  counts repetitions equal to the best observed value.
* `comparison.csv`: per instance, the best observed value, the methods reaching it,
  and each method's min, mean and gaps.
* `methods.csv`: per method across instances: status counts, instances where it
  reached the best observed value, and mean gaps.
* `gap_pct(x) = 100·(x − ref)/ref`, where `ref` is the best value observed for the
  instance by any method. With `--reference best_known.csv` (`instance,value`), the
  `*_to_reference` columns use those values, and the summary gaps use the lower of
  the two.

## Tests

```bash
cd experiments/tests
python3 -m unittest -v                       # fake solvers; ~1 min
LINKERN_PATH=/path/to/linkern python3 -m unittest -v test_integration   # real solvers; ~5 min
```

`test_runner.py` covers statuses and schema conformance, invalid outputs, soft vs
hard limits (including grandchildren and SIGTERM-ignoring solvers), seed
derivation, private work dirs, parallel = serial results, CPU budget and pinning,
CBMLKH exclusion across two runner processes, no double claims, crash (SIGKILL)
and graceful-interrupt resume with orphan cleanup, lost-record recovery, retry
policy, design protection and CSV/record consistency. `test_seeds.py` compiles
`cbm_seed.h` and checks the Python mirror against it. `test_integration.py` runs
all four real solvers twice (4 CPUs vs 2, with an interrupted session) and
requires identical permutations.
