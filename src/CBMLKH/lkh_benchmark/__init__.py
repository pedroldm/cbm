"""Benchmark harness to compare pure LKH against the CBMLKH solver.

Pure LKH converts each instance to TSP exactly as the C++ ``LKHWrapper`` does and
solves the whole instance; the CBMLKH solver shells out to the C++ binary. Both
run through the same resumable benchmark loop so their results are directly
comparable.

Every tunable parameter lives in :mod:`lkh_benchmark.config`
(:class:`BenchmarkConfig`). Typical programmatic use (e.g. from a notebook)::

    from lkh_benchmark import BenchmarkConfig, LKHParams, run_comparison

    config = BenchmarkConfig(
        instances=["i3", "scpd1"],
        lkh=LKHParams(time_limit=120),
        cbmlkh_config={"maxTime": 120, "threads": 8},
    )
    rows = run_comparison(config)

For finer control you can still drive a single solver directly::

    from lkh_benchmark import Benchmark, PureLKHSolver
    Benchmark(PureLKHSolver(LKHParams(time_limit=120)), instances_dir, results_dir).run(["i3"])
"""

from .benchmark import Benchmark, discover_instances
from .config import (
    DEFAULT_CBMLKH_CONFIG,
    DEFAULT_LKH_PATH,
    BenchmarkConfig,
    LKHParams,
)
from .instance import CBMInstance
from .lkh_runner import LKHRunner
from .pipeline import run_comparison
from .report import build_comparison, print_comparison, write_comparison_csv
from .result_store import ResultStore
from .solvers import CBMLKHSolver, PureLKHSolver, Solver, find_cbmlkh_binary
from .tsp_converter import TSPConverter
from .validator import SolutionValidator, ValidationError

__all__ = [
    "BenchmarkConfig",
    "run_comparison",
    "Benchmark",
    "discover_instances",
    "Solver",
    "PureLKHSolver",
    "CBMLKHSolver",
    "DEFAULT_CBMLKH_CONFIG",
    "find_cbmlkh_binary",
    "CBMInstance",
    "LKHParams",
    "LKHRunner",
    "DEFAULT_LKH_PATH",
    "ResultStore",
    "TSPConverter",
    "SolutionValidator",
    "ValidationError",
    "build_comparison",
    "print_comparison",
    "write_comparison_csv",
]
