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

    The two solvers run in **strictly sequential phases** and never overlap: the
    LKH phase runs to completion (its thread pool joins, so every LKH subprocess
    has exited) before the CBMLKH phase starts. This keeps LKH (many
    single-threaded processes) and CBMLKH (one internally multi-threaded process)
    from ever competing for cores.
    """
    instance_names = discover_instances(config.instances_dir) if config.run_all else list(config.instances)
    if not instance_names:
        print(f"{LOG_PREFIX} no instances selected; nothing to do.")
        return []

    labeled_tags: dict[str, str] = {}

    # ---- Phase 1: pure LKH (instances solved concurrently; each LKH is 1 thread).
    if config.run_lkh:
        labeled_tags["lkh"] = _run_lkh_phase(config, instance_names)

    # ---- Phase 2: CBMLKH. Only reached after phase 1 has fully drained above, so
    #      the two solvers are never in flight at the same time.
    if config.run_cbmlkh:
        labeled_tags["cbmlkh"] = _run_cbmlkh_phase(config, instance_names)

    if not labeled_tags:
        return []

    rows = build_comparison(config.results_dir, labeled_tags)
    out_csv = Path(config.results_dir) / COMPARISON_CSV
    write_comparison_csv(rows, out_csv)
    print()
    print_comparison(rows, list(labeled_tags))
    print(f"\n{LOG_PREFIX} comparison written to {out_csv}")
    return rows


def _run_lkh_phase(config: BenchmarkConfig, instance_names: list[str]) -> str:
    """Run the LKH phase to completion and return its run tag.

    ``Benchmark.run`` uses a ``ThreadPoolExecutor`` context manager, so this call
    only returns once every worker thread — and therefore every blocking LKH
    subprocess — has finished. The CBMLKH phase cannot start until then.
    """
    print(f"{LOG_PREFIX} === phase: LKH (workers={config.lkh_workers}) ===")
    solver = PureLKHSolver(params=config.lkh, work_dir=Path(config.results_dir) / LKH_WORK_DIRNAME)
    bench = Benchmark(solver, config.instances_dir, config.results_dir, run_tag=config.lkh_run_tag)
    bench.run(instance_names, skip_done=config.skip_done, workers=config.lkh_workers)
    print(f"{LOG_PREFIX} === phase: LKH complete ===")
    return bench.run_tag


def _run_cbmlkh_phase(config: BenchmarkConfig, instance_names: list[str]) -> str:
    """Run the CBMLKH phase to completion and return its run tag.

    Instances run one at a time (``workers=1``): CBMLKH is already internally
    multi-threaded, so running several at once would oversubscribe cores.
    """
    print(f"{LOG_PREFIX} === phase: CBMLKH ===")
    solver = CBMLKHSolver(
        config=config.cbmlkh_config,
        binary_path=config.cbmlkh_binary,
        work_dir=Path(config.results_dir) / CBMLKH_WORK_DIRNAME,
        timeout=config.cbmlkh_timeout,
    )
    bench = Benchmark(solver, config.instances_dir, config.results_dir, run_tag=config.cbmlkh_run_tag)
    bench.run(instance_names, skip_done=config.skip_done, workers=1)
    print(f"{LOG_PREFIX} === phase: CBMLKH complete ===")
    return bench.run_tag
