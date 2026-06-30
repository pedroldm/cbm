"""Run the configured solvers over the chosen instances and compare them.

A thin orchestration layer driven entirely by a :class:`BenchmarkConfig`, so the
same entry point works from the CLI (``python -m lkh_benchmark``) and from a
notebook (``run_comparison(BenchmarkConfig(...))``).
"""

from __future__ import annotations

from pathlib import Path

from .benchmark import Benchmark, discover_instances
from .config import COMPARISON_CSV, CBMLKH_WORK_DIRNAME, LKH_WORK_DIRNAME, LOG_PREFIX, BenchmarkConfig
from .report import build_comparison, print_comparison, write_comparison_csv
from .solvers import CBMLKHSolver, PureLKHSolver


def run_comparison(config: BenchmarkConfig) -> list[dict]:
    """Execute the benchmark described by ``config`` and return comparison rows.

    Execution is **interleaved per instance**: every enabled solver runs on the
    first instance, then every solver on the second, and so on (a1 LKH, a1 CBMLKH,
    a2 LKH, a2 CBMLKH, …). Solvers run one at a time, so LKH and CBMLKH are never
    in flight simultaneously and never compete for cores.
    """
    instance_names = discover_instances(config.instances_dir) if config.run_all else list(config.instances)
    if not instance_names:
        print(f"{LOG_PREFIX} no instances selected; nothing to do.")
        return []

    benches = _build_benchmarks(config)
    if not benches:
        return []

    total = len(instance_names)
    for i, name in enumerate(instance_names, start=1):
        for bench in benches.values():
            bench.process_one(name, i, total, skip_done=config.skip_done)

    labeled_tags = {label: bench.run_tag for label, bench in benches.items()}

    rows = build_comparison(config.results_dir, labeled_tags)
    out_csv = Path(config.results_dir) / COMPARISON_CSV
    write_comparison_csv(rows, out_csv)
    print()
    print_comparison(rows, list(labeled_tags))
    print(f"\n{LOG_PREFIX} comparison written to {out_csv}")
    return rows


def _build_benchmarks(config: BenchmarkConfig) -> dict[str, Benchmark]:
    """Construct a ``Benchmark`` per enabled solver, keyed by label, in run order."""
    benches: dict[str, Benchmark] = {}

    if config.run_lkh:
        solver = PureLKHSolver(params=config.lkh, work_dir=Path(config.results_dir) / LKH_WORK_DIRNAME)
        benches["lkh"] = Benchmark(solver, config.instances_dir, config.results_dir, run_tag=config.lkh_run_tag)

    if config.run_cbmlkh:
        solver = CBMLKHSolver(
            config=config.cbmlkh_config,
            binary_path=config.cbmlkh_binary,
            work_dir=Path(config.results_dir) / CBMLKH_WORK_DIRNAME,
            timeout=config.cbmlkh_timeout,
        )
        benches["cbmlkh"] = Benchmark(solver, config.instances_dir, config.results_dir, run_tag=config.cbmlkh_run_tag)

    return benches
