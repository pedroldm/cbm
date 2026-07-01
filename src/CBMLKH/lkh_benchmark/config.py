"""Single source of truth for every benchmark parameter.

**Edit parameters only in the PARAMETERS block below.** Every tunable value the
benchmark uses is a module-level constant in that one block; the config
dataclasses (:class:`BenchmarkConfig`, :class:`LKHParams`) merely reference those
constants as their defaults and never redefine or override them. There are no
placeholder defaults anywhere else in the codebase, so changing a value here is
guaranteed to take effect for every run.

``config.py`` only depends on the standard library, so it never causes import
cycles.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

# Path anchor used to derive the default parameter paths below. src/CBMLKH is one
# level up from this package directory.
CBMLKH_DIR = Path(__file__).resolve().parents[1]


# =========================================================================== #
#  PARAMETERS — the single source of truth. Edit values HERE and nowhere else.
# =========================================================================== #

# ---- Paths ----------------------------------------------------------------- #
INSTANCES_DIR = str(CBMLKH_DIR.parent.parent / "instances")
RESULTS_DIR = str(CBMLKH_DIR / "lkh_benchmark" / "results")

# ---- Instance selection ---------------------------------------------------- #
# RUN_ALL runs every file in INSTANCES_DIR; otherwise the explicit list is used.
RUN_ALL = False
INSTANCES = ["a1", "a2", "i3"]

# ---- Which solvers to run (run both to compare them on the same instances) -- #
RUN_LKH = True
RUN_CBMLKH = True

# ---- Pure-LKH parameters (forwarded to LKH; match ``LKHWrapper::writePar``) - #
LKH_PATH = str(CBMLKH_DIR.parent / "LKH3" / "LKH")  # external LKH3 executable
LKH_TIME_LIMIT = 300  # seconds (TIME_LIMIT)
LKH_RUNS = 1
LKH_MOVE_TYPE = 5
LKH_PATCHING_C = 3
LKH_PATCHING_A = 2
LKH_SEED: int | None = None
# Any additional LKH parameter as {"KEY": value} (e.g. {"MAX_TRIALS": 1000}).
LKH_EXTRA: dict = {}
# Parallelism for a *standalone* LKH batch run via ``Benchmark.run(..., workers=N)``
# (each LKH process is single-threaded). The interleaved comparison pipeline runs
# one solver at a time, so this does NOT affect ``run_comparison``.
LKH_WORKERS = 4

# ---- CBMLKH parameters (the C++ solver) ------------------------------------ #
# cbmlkh_binary=None auto-detects the binary in src/CBMLKH.
CBMLKH_BINARY: str | None = None
# Subprocess hard timeout in seconds (None = no cap).
CBMLKH_TIMEOUT: float | None = None
# Config fed to the C++ binary (keys mirror parameters.md; instancePath is set
# per run). Every value here is used verbatim — nothing overrides it downstream.
CBMLKH_CONFIG: dict = {
    "threads": 8,
    "blockMovement": "RANDOM",
    "maxIterations": 1000,
    "maxTime": 300,  # seconds
    "lkhMaxTime": 1,  # minutes
    "constructionBias": 2.5,
    "neighborBias": 1.0,
    "minNeighborBias": 0.2,
    "minSegmentScore": 10.0,
    "minSegmentScoreLowerBound": 2.0,
    "maxSegmentSize": 0.1,  # fraction of column count
    "maxSegmentSizeUpperBound": 0.2,  # fraction of column count
    "segmentSizeGrowthFactor": 1.25,
    "segmentScoreDecayFactor": 0.85,
    "neighborBiasDecayFactor": 0.90,
    "adaptationInterval": 20,
}

# ---- Output / resume control ----------------------------------------------- #
# run_tag=None derives the results subdir from the solver's parameters.
LKH_RUN_TAG: str | None = None
CBMLKH_RUN_TAG: str | None = None
SKIP_DONE = True  # skip instances that already have a result (never overwrites)


# =========================================================================== #
#  Fixed library constants (internal wiring, not run parameters)
# =========================================================================== #

# CBMLKH C++ binaries to look for (in order) when none is given explicitly.
CBMLKH_BINARY_CANDIDATES = ("main_prd", "main", "main_debug")

# Scratch subdirectory names (.tsp/.par/.tour and generated .cfg files).
LKH_WORK_DIRNAME = "_lkh_work"
CBMLKH_WORK_DIRNAME = "_cbmlkh_cfg"

# Output / logging.
COMPARISON_CSV = "comparison.csv"
LOG_PREFIX = "[lkh_benchmark]"

# Columns written to each solver's index.csv summary.
INDEX_FIELDS = ["instance", "solver", "rows", "cols", "cost", "runtimeSec", "validated", "timestamp"]


# =========================================================================== #
#  Config dataclasses — defaults reference the PARAMETERS block above only.
# =========================================================================== #
@dataclass
class LKHParams:
    """Knobs forwarded to LKH. Defaults come from the PARAMETERS block above."""

    lkh_path: str = LKH_PATH
    time_limit: int = LKH_TIME_LIMIT  # seconds (TIME_LIMIT)
    runs: int = LKH_RUNS
    move_type: int = LKH_MOVE_TYPE
    patching_c: int = LKH_PATCHING_C
    patching_a: int = LKH_PATCHING_A
    seed: int | None = LKH_SEED
    extra: dict = field(default_factory=lambda: dict(LKH_EXTRA))

    def to_par_lines(self, problem_file: str, tour_file: str) -> list[str]:
        lines = [
            f"PROBLEM_FILE = {problem_file}",
            f"TOUR_FILE = {tour_file}",
            f"MOVE_TYPE = {self.move_type}",
            f"PATCHING_C = {self.patching_c}",
            f"PATCHING_A = {self.patching_a}",
            f"RUNS = {self.runs}",
            f"TIME_LIMIT = {self.time_limit}",
        ]
        if self.seed is not None:
            lines.append(f"SEED = {self.seed}")
        for key, value in self.extra.items():
            lines.append(f"{key} = {value}")
        return lines

    def signature(self) -> str:
        """Short, filesystem-safe tag identifying this configuration."""
        parts = [f"t{self.time_limit}", f"r{self.runs}", f"mt{self.move_type}"]
        if self.seed is not None:
            parts.append(f"s{self.seed}")
        for key, value in sorted(self.extra.items()):
            parts.append(f"{key}{value}")
        return "_".join(parts)

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class BenchmarkConfig:
    """Everything a run needs. Defaults come from the PARAMETERS block above;
    construct one explicitly (e.g. from a notebook) only to deviate from them."""

    # Paths.
    instances_dir: str = INSTANCES_DIR
    results_dir: str = RESULTS_DIR

    # Instance selection.
    run_all: bool = RUN_ALL
    instances: list[str] = field(default_factory=lambda: list(INSTANCES))

    # Which solvers to run.
    run_lkh: bool = RUN_LKH
    run_cbmlkh: bool = RUN_CBMLKH

    # Pure-LKH parameters.
    lkh: LKHParams = field(default_factory=LKHParams)
    lkh_workers: int = LKH_WORKERS

    # CBMLKH parameters (the C++ solver).
    cbmlkh_binary: str | None = CBMLKH_BINARY
    cbmlkh_config: dict = field(default_factory=lambda: dict(CBMLKH_CONFIG))
    cbmlkh_timeout: float | None = CBMLKH_TIMEOUT

    # Output / resume control.
    lkh_run_tag: str | None = LKH_RUN_TAG
    cbmlkh_run_tag: str | None = CBMLKH_RUN_TAG
    skip_done: bool = SKIP_DONE
