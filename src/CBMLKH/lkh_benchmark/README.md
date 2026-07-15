# lkh_benchmark — compare pure LKH vs. CBMLKH

A small OOP Python package that runs **plain LKH** and the **CBMLKH** C++ solver
over the same CBM instances and reports them side by side.

- **Pure LKH** converts the whole instance to TSP exactly like the C++
  `LKHWrapper::writeTSP` — an `EXPLICIT` / `FULL_MATRIX` problem with a dummy
  depot node 0 (`distance(depot, col) = onesCount(col)`, turning the cycle into
  an open Hamiltonian path) and Hamming distances between columns — solves it,
  parses the tour like `LKHWrapper::getResultTour` (`node - 2`, rotate the depot
  to the front and drop it), and **validates** it like the C++ `Validator` (the
  permutation must contain every column once and the reported cost must equal the
  recomputed 1-block count, computed two independent ways).
- **CBMLKH** shells out to the C++ binary with a generated config and parses its
  JSON report (`global.bestCost`, `runtimeMs`, …). The binary validates every
  improving solution internally, so its `bestCost` is trusted.

Both run through the same resumable benchmark loop, so a long `nohup` run can be
restarted without recomputing or overwriting finished work.

Execution is **interleaved per instance**: every enabled solver runs on the first
instance, then every solver on the next, and so on (a1 LKH, a1 CBMLKH, a2 LKH, a2
CBMLKH, …). Solvers run one at a time, so LKH and CBMLKH are never in flight
simultaneously and never compete for cores.

## Layout

| File | Responsibility |
|------|----------------|
| **`config.py`**      | **Single source of truth for every parameter** — a top-level `PARAMETERS` block holds every tunable value (`LKH_*`, `CBMLKH_CONFIG`, paths, selection, resume flags); `BenchmarkConfig` and `LKHParams` only reference those constants as defaults, never redefining them. |
| `instance.py`        | `CBMInstance` — parse the sparse format; `onesCount`, `onesToZeros`, `hamming_matrix`, `path_cost`, `count_one_blocks`. |
| `tsp_converter.py`   | `TSPConverter` — build & write the depot+Hamming TSP matrix. |
| `lkh_runner.py`      | `LKHRunner` (write par, shell out, parse tour). |
| `validator.py`       | `SolutionValidator` / `ValidationError`. |
| `solvers/`           | `Solver` ABC + `PureLKHSolver` and `CBMLKHSolver`. |
| `result_store.py`    | `ResultStore` — resumable, append-only persistence (one JSON per instance + `index.csv`). |
| `benchmark.py`       | `Benchmark` — solver-agnostic orchestration + resume; `discover_instances`. |
| `report.py`          | `build_comparison` / `write_comparison_csv` / `print_comparison`. |
| `pipeline.py`        | `run_comparison(config)` — runs the configured solvers and compares them. |
| `__main__.py`        | Thin entry point: `run_comparison(BenchmarkConfig())`. |

## Running (for `nohup`)

**All parameters live in the `PARAMETERS` block at the top of `config.py`.**
Edit the values there (nowhere else), then run the module — there are no
command-line arguments:

```bash
cd src/CBMLKH
nohup python -m lkh_benchmark > lkh_run.log 2>&1 &
```

`BenchmarkConfig` fields:

- **Paths/selection**: `instances_dir`, `results_dir`; `run_all` (every file in
  `instances_dir`) or the explicit `instances` list.
- **Which solvers**: `run_lkh`, `run_cbmlkh` (enable both to compare).
- **Pure LKH**: `lkh` — an `LKHParams(lkh_path, time_limit, runs, move_type,
  patching_c, patching_a, seed, extra)`. (`lkh_workers` only applies to a
  standalone `Benchmark.run(..., workers=N)` batch; the interleaved comparison
  runs one solver at a time, so it does not affect `run_comparison`.)
- **CBMLKH**: `cbmlkh_binary` (auto-detected if `None`), `cbmlkh_config` (mirrors
  `parameters.md`; `maxTime` is the per-instance budget), `cbmlkh_timeout`.
- **Resume**: `lkh_run_tag`, `cbmlkh_run_tag`, `skip_done`.

At the end it prints a comparison table and writes `<results_dir>/comparison.csv`.

## Resuming

Each solver writes to `<results-dir>/<run-tag>/`, where `run-tag` defaults to a
signature of that solver's parameters (e.g. `lkh_t300_r1_mt5`,
`cbmlkh_t300_th8_RANDOM`). A finished instance has a `<name>.json` file;
re-running the **same config** skips those and **never overwrites** existing
results. Different parameters → different `run-tag` → isolated directory, so
configurations never collide.

## Programmatic / notebook usage

Build a `BenchmarkConfig` and hand it to `run_comparison`:

```python
from lkh_benchmark import BenchmarkConfig, LKHParams, run_comparison

config = BenchmarkConfig(
    instances=["i3", "scpd1"],
    results_dir="./results",
    lkh=LKHParams(time_limit=120),
    cbmlkh_config={"maxTime": 120, "threads": 8},
)
rows = run_comparison(config)
```

For finer control you can still drive a single solver directly:

```python
from lkh_benchmark import Benchmark, PureLKHSolver, LKHParams
Benchmark(PureLKHSolver(LKHParams(time_limit=120)), instances_dir, results_dir).run(["i3"])
```

## Stored fields (per-instance JSON)

Common: `instance`, `rows`, `cols`, `solver`, `cost`, `runtimeSec`, `runTag`,
`timestamp`.

- Pure LKH adds: `blockCount`, `validated`, `validationError`, `solution` (the
  column permutation), `lkhParams`.
- CBMLKH adds: `wallSec`, `config` (the input config), and `report` — the
  solver's **complete** end-of-run JSON: `global` metrics (best cost, accepted /
  rejected moves, `lkhCache` stats, per-operator `operators` effectiveness) and
  every entry of `trajectories` with its full metrics (`blockSize`, `lkhCalls`,
  `diversifications`, `neighborBiasHistory`, improvement `history`, …). Each
  `history` entry additionally carries a `histogram`: the per-column new-block
  contribution of that improving solution (length `cols`, `histogram[i]` = 1-blocks
  opened at column `i`, summing to the entry's `cost`) — a column-by-column view of
  how dense regions flatten over the trajectory. `validated` is `null` — the C++
  binary validates internally.

## Requirements

`numpy`, a built LKH binary (`cd src/LKH3 && make`; default path
`/home/pedroldm/MSc/cbm/src/LKH3/LKH`), and a built CBMLKH binary (`make prd` —
the package auto-detects `main_prd` / `main` / `main_debug`).
