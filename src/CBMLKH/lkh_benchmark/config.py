"""Single source of truth for every benchmark parameter.

Everything tunable lives here:

- :class:`BenchmarkConfig` — the run specification you edit (which instances,
  which solvers, their parameters, resume/output options).
- :class:`LKHParams` — the pure-LKH knobs (also a field of ``BenchmarkConfig``).
- The fixed library constants (binary locations, scratch dir names, the index
  columns, the log prefix).

Every other module imports its defaults from this file, so there is exactly one
place to change anything. ``config.py`` only depends on the standard library, so
it never causes import cycles.
"""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field
from pathlib import Path

# --------------------------------------------------------------------------- #
# Fixed locations and library constants (rarely changed)
# --------------------------------------------------------------------------- #

# src/CBMLKH (one level up from this package directory).
CBMLKH_DIR = Path(__file__).resolve().parents[1]

# External LKH3 executable used by the pure-LKH solver (src/LKH3/LKH).
DEFAULT_LKH_PATH = str(CBMLKH_DIR.parent / "LKH3" / "LKH")

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

# Default CBMLKH config, mirroring parameters.md (instancePath is set per run).
DEFAULT_CBMLKH_CONFIG: dict = {
    "threads": 8,
    "blockMovement": "RANDOM",
    "maxIterations": 1000,
    "maxTime": 60,  # seconds
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


# --------------------------------------------------------------------------- #
# Pure-LKH parameters
# --------------------------------------------------------------------------- #
@dataclass
class LKHParams:
    """Knobs forwarded to LKH. Defaults match ``LKHWrapper::writePar``.

    ``extra`` lets you pass any additional LKH parameter as ``{"KEY": value}``
    without touching this class (e.g. ``{"MAX_TRIALS": 1000}``).
    """

    lkh_path: str = DEFAULT_LKH_PATH
    time_limit: int = 60  # seconds (TIME_LIMIT)
    runs: int = 1
    move_type: int = 5
    patching_c: int = 3
    patching_a: int = 2
    seed: int | None = None
    extra: dict = field(default_factory=dict)

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


# --------------------------------------------------------------------------- #
# Full run specification (edit the defaults below to configure a run)
# --------------------------------------------------------------------------- #
@dataclass
class BenchmarkConfig:
    """Everything a run needs. Edit the field defaults here, or construct one
    explicitly (e.g. from a notebook) and pass it to ``run_comparison``."""

    # Paths.
    instances_dir: str = str(CBMLKH_DIR.parent.parent / "instances")
    results_dir: str = str(CBMLKH_DIR / "lkh_benchmark" / "results")

    # Instance selection: RUN_ALL runs every file in instances_dir; otherwise the
    # explicit `instances` list is used.
    run_all: bool = False
    instances: list[str] = field(default_factory=lambda: ["a1", "a2", "i3"])

    # Which solvers to run (run both to compare them on the same instances).
    run_lkh: bool = True
    run_cbmlkh: bool = True

    # Pure-LKH parameters.
    lkh: LKHParams = field(default_factory=lambda: LKHParams(time_limit=300))

    # Parallelism for a *standalone* LKH batch run via ``Benchmark.run(..., workers=N)``
    # (each LKH process is single-threaded). The interleaved comparison pipeline
    # runs one solver at a time so LKH and CBMLKH never overlap, so this does NOT
    # affect ``run_comparison``.
    lkh_workers: int = 4

    # CBMLKH parameters (the C++ solver). cbmlkh_binary=None auto-detects.
    cbmlkh_binary: str | None = None
    cbmlkh_config: dict = field(default_factory=lambda: {**DEFAULT_CBMLKH_CONFIG, "maxTime": 300, "threads": 8})
    cbmlkh_timeout: float | None = None  # subprocess hard timeout in seconds (None = no cap)

    # Output / resume control. run_tag=None derives the results subdir from params.
    lkh_run_tag: str | None = None
    cbmlkh_run_tag: str | None = None
    skip_done: bool = True  # skip instances that already have a result (never overwrites)
