"""Per-method adapters: defaults, resource estimates, command lines, output parsing.

Audit of internal parallelism (verified in the sources, see README):
  * ENS, ILS: single-threaded; each TSP call is one single-threaded Linkern/LKH
    child process, run sequentially.
  * standalone LKH: one single-threaded LKH process (LKH-3 has no threads).
  * CBMLKH: OpenMP with exactly ``threads`` workers, each of which runs one
    LKH child at a time, i.e. ``threads`` busy CPUs.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from .seeds import cbm_derive_seed

GIB = 1024**3

METHOD_IDS = ("ens", "ils", "lkh", "cbmlkh")
METHOD_NAMES = {"ens": "ENS", "ils": "ILS", "lkh": "LKH", "cbmlkh": "CBMLKH"}

# MOVE_TYPE / PATCHING_C / PATCHING_A of every LKH call (src/common/cbm_lkh_params.h,
# compiled into CBMLKH, lkh_standalone, ENS and ILS). Recorded here for the
# experiment metadata; they cannot be overridden.
LKH_PARAMS = {"move_type": 5, "patching_c": 3, "patching_a": 2}

# CBMLKH limits each LKH call (lkhMaxTime, LKH's TIME_LIMIT) to this fraction
# of the run's time budget: 0.2 * 7200 s = 1440 s.
CBMLKH_LKH_TIME_FRACTION = 0.2

# Deliberate departures from the irace configuration below, for the experiments:
#  * maxIterations 1000 (irace: 1000000, i.e. effectively time-bounded);
#  * every block handed to LKH spans 10%..25% of the columns: regions (PEAK /
#    INTERVAL windows) start at 10% (minSegmentSizeFraction, maxSegmentSize) and
#    may widen to 12.5% (maxSegmentSizeUpperBound, rounded down; one 1.5x
#    diversification step), and a MERGE of two regions is capped at 2 x 12.5% = 25%
#    (irace: minSegmentSize 5 columns, regions 30%..60%, MERGE up to 100%).
CBMLKH_EXPERIMENT_OVERRIDES = {
    "maxIterations": 1000,
    "minSegmentSizeFraction": 0.10,
    "maxSegmentSize": 0.10,
    "maxSegmentSizeUpperBound": 0.125,
}
CBMLKH_MAX_ITERATIONS = CBMLKH_EXPERIMENT_OVERRIDES["maxIterations"]

# irace 4.2.0 best-so-far configuration 135 (last completed race, 2026-08-20;
# tunning/output.txt), verified by tests/test_tuning.py. Deviations, applied in
# DEFAULT_PARAMS: CBMLKH_EXPERIMENT_OVERRIDES above, and lkhMaxTime (derived from
# the budget by plan.make_spec; 300 in irace).
CBMLKH_IRACE_CONFIG = {
    "blockMovement": "MERGE",
    "maxIterations": 1000000,
    "maxSegmentSize": 0.30,
    "maxSegmentSizeUpperBound": 0.60,
    "segmentSizeGrowthFactor": 1.5,
    "minSegmentSize": 5,
    "neighborBias": 2.0,
    "minNeighborBias": 0.1,
    "neighborBiasDecayFactor": 0.90,
    "minSegmentScore": 2.0,
    "minSegmentScoreLowerBound": 0.15,
    "segmentScoreDecayFactor": 0.85,
    "adaptationInterval": 10,
    "constructionBias": 2.5,
}

# Parameters an experiment may not override.
LOCKED_PARAMS = {"cbmlkh": "all", "lkh": set(LKH_PARAMS) | {"tsp_backend"}}


def cbmlkh_lkh_max_time(soft_limit_s: float) -> int:
    return max(1, int(round(CBMLKH_LKH_TIME_FRACTION * soft_limit_s)))


DEFAULT_PARAMS: Dict[str, Dict[str, Any]] = {
    # Haddadi (2021): NB_ITER = 500, RATIO = 5, Linkern.
    "ens": {"tsp_backend": "linkern", "iterations": 500},
    # ILS defaults of ILS/ils.cpp: NBITER = 24, phi = 0.2, V = 1; Linkern.
    "ils": {"tsp_backend": "linkern", "iterations": 24, "phi": 0.2, "V": 1, "runs": 1},
    # The whole instance as one TSP, with the project-wide LKH parameters.
    "lkh": {"tsp_backend": "lkh", "runs": 1, **LKH_PARAMS},
    # threads, maxTime, seed and paths are set by the runner; lkhMaxTime is
    # added by plan.make_spec from the time budget.
    "cbmlkh": {"tsp_backend": "lkh", **CBMLKH_IRACE_CONFIG, **CBMLKH_EXPERIMENT_OVERRIDES, **{"lkh_" + k: v for k, v in LKH_PARAMS.items()}},
}



@dataclass(frozen=True)
class Estimate:
    mem_bytes: int
    disk_bytes: int


def estimate_resources(method: str, rows: int, cols: int, params: Dict[str, Any], threads: int) -> Estimate:
    """Peak memory and scratch-disk use, from the data structures each code allocates.

    Heuristic but conservative for ENS/ILS/LKH (dominated by explicit distance
    matrices); for CBMLKH it assumes the largest possible sub-problem (a MERGE
    of two maximally widened regions) in every thread at once.
    """
    slack = 64 * 1024**2
    n = cols
    if method == "ens":
        # a, a_best, b, orig (chars) + int (n+2)^2 matrix while writing the TSP;
        # Linkern then holds ~2n^2 bytes. File: FULL_MATRIX, 6 bytes per entry.
        return Estimate(4 * rows * (n + 1) + 4 * (n + 2) ** 2 + slack, 6 * (n + 1) ** 2)
    if method == "ils":
        # W and its perturbed copy Wp, both (c+2)^2 ints. File: FULL_MATRIX.
        return Estimate(8 * (n + 2) ** 2 + slack, 6 * (n + 2) ** 2)
    if method == "lkh":
        # LKH's triangular int cost matrix + per-node structures. File: UPPER_ROW.
        return Estimate(2 * (n + 1) ** 2 + 2048 * n + slack, 3 * (n + 1) ** 2)
    if method == "cbmlkh":
        # MERGE spans at most two regions of maxSegmentSizeUpperBound columns.
        seg = min(n, 2 * int(math.floor(float(params.get("maxSegmentSizeUpperBound", 1.0)) * n + 1e-9)))
        per_thread_mem = 2 * (seg + 1) ** 2 + 2048 * seg
        per_thread_disk = 6 * (seg + 1) ** 2
        return Estimate(threads * per_thread_mem + rows * n // 8 + slack, threads * per_thread_disk)
    raise ValueError(f"unknown method {method}")


def build_command(job: dict, spec: dict, paths: dict, native_out: str, work_dir: str, cfg_file: Optional[str]) -> List[str]:
    """argv for one job attempt. For CBMLKH the config file is written by write_cbmlkh_config()."""
    method, params = job["method"], spec["method_params"][job["method"]]
    inst = job["instance"]["path"]
    soft = spec["soft_limit_s"]
    if method in ("ens", "ils"):
        backend = params["tsp_backend"]
        argv = [
            paths[method],
            f"--filePath={inst}",
            f"--algorithm={backend}",
            f"--solverPath={paths['linkern'] if backend == 'linkern' else paths['lkh']}",
            f"--outputPath={native_out}",
            f"--seed={job['seed']}",
            f"--timeLimit={soft:g}",
            f"--iterations={params['iterations']}",
            f"--workDir={work_dir}",
        ]
        if method == "ils":
            argv += [f"--phi={params['phi']}", f"--V={params['V']}", f"--runs={params['runs']}"]
        return argv
    if method == "lkh":
        return [
            paths["lkh_standalone"],
            f"--filePath={inst}",
            f"--lkhPath={paths['lkh']}",
            f"--seed={job['seed']}",
            f"--timeLimit={soft:g}",
            f"--outputPath={native_out}",
            f"--workDir={work_dir}",
            f"--runs={params['runs']}",
        ]
    if method == "cbmlkh":
        return [paths["cbmlkh"], cfg_file]
    raise ValueError(method)


def cbmlkh_config_text(job: dict, spec: dict, paths: dict, native_out: str, work_dir: str) -> str:
    params = spec["method_params"]["cbmlkh"]
    lines = [
        f"instancePath={job['instance']['path']}",
        "iRace=false",
        f"lkhPath={paths['lkh']}",
        f"lkhTmpDir={work_dir}",
        f"outputPath={native_out}",
        f"seed={job['seed']}",
        f"threads={job['threads']}",
        f"trajectoryOffset={job['repetitions'][0]}",
        # CBMLKH's maxTime is an integer number of seconds.
        f"maxTime={int(math.ceil(spec['soft_limit_s']))}",
    ]
    # tsp_backend and lkh_* are documentation only (LKH parameters are compiled in).
    lines += [f"{k}={v}" for k, v in sorted(params.items()) if k != "tsp_backend" and not k.startswith("lkh_")]
    return "\n".join(lines) + "\n"


def _num(value: Any) -> Optional[float]:
    return float(value) if isinstance(value, (int, float)) and not isinstance(value, bool) else None


def parse_native(job: dict, native: dict) -> Dict[int, dict]:
    """Normalize a solver's own JSON into per-repetition partial records.

    Returns {repetition: {...}} with the fields the record builder expects.
    Raises ValueError if the output is inconsistent with the job.
    """
    method = job["method"]
    if method in ("ens", "ils"):
        tsp = native["tsp"]
        derived = {"tsp_call_seeds": tsp["seeds"], "tsp_call_seed_rule": "cbm_derive_seed(run_seed, call_index); call 0 = initial tour"}
        if method == "ils":
            derived["perturbation_seed"] = native["perturbation_seed"]
        if int(native["seed"]) != job["seed"]:
            raise ValueError(f"solver ran with seed {native['seed']}, expected {job['seed']}")
        rep = job["repetitions"][0]
        return {
            rep: {
                "value": native["best_blocks"],
                "initial_value": native["initial_blocks"],
                "permutation": native["permutation"],
                "method_time_s": _num(native["elapsed_s"]),
                "time_to_best_s": _num(native["time_to_best_s"]),
                "iterations": native["iterations_completed"],
                "improvements": max(0, len(native["history"]) - 1),
                "stop_reason": native["stop_reason"],
                "method_time_limit_reached": bool(native["time_limit_reached"]),
                "tsp": {"calls": tsp["calls"], "time_s": tsp["time_s"], "time_limit_hits": tsp["time_limit_hits"], "invalid_tours": tsp["invalid_tours"]},
                "derived_seeds": derived,
                "method_stats": {k: native[k] for k in ("parameters", "history", "validated")},
            }
        }
    if method == "lkh":
        if int(native["seed"]) != job["seed"]:
            raise ValueError(f"solver ran with seed {native['seed']}, expected {job['seed']}")
        parts = [native.get("tsp_build_time_s"), native.get("lkh_preprocessing_time_s"), native.get("lkh_time_to_best_s")]
        time_to_best = sum(parts) if all(isinstance(p, (int, float)) for p in parts) else None
        hit = bool(native["time_limit_reached"])
        rep = job["repetitions"][0]
        return {
            rep: {
                "value": native["best_blocks"],
                "initial_value": None,
                "permutation": native["permutation"],
                "method_time_s": _num(native["elapsed_s"]),
                "time_to_best_s": time_to_best,
                "iterations": None,
                "improvements": native.get("lkh_improvements"),
                "stop_reason": native["stop_reason"],
                "method_time_limit_reached": hit,
                "tsp": {"calls": 1, "time_s": native["lkh_time_s"], "time_limit_hits": int(hit), "invalid_tours": 0},
                "derived_seeds": {"lkh_seed": job["seed"]},
                "method_stats": {
                    k: native.get(k)
                    for k in ("parameters", "tsp_build_time_s", "lkh_time_s", "lkh_preprocessing_time_s", "lkh_run_time_s", "lkh_time_to_best_s", "lkh_improvements")
                },
            }
        }
    if method == "cbmlkh":
        trajectories = {int(t["index"]): t for t in native["trajectories"]}
        if sorted(trajectories) != sorted(job["repetitions"]):
            raise ValueError(f"CBMLKH reported trajectories {sorted(trajectories)}, expected {job['repetitions']}")
        cfg = native["config"]
        if int(cfg["seed"]) != job["seed"]:
            raise ValueError(f"CBMLKH ran with seed {cfg['seed']}, expected {job['seed']}")
        execution = {
            "threads": cfg["threads"],
            "runtime_ms": native["global"]["runtimeMs"],
            "best_cost": native["global"]["bestCost"],
            "lkh_cache": native["global"]["lkhCache"],
            "resolved_config": cfg,
        }
        out = {}
        for rep, t in trajectories.items():
            expected_seed = cbm_derive_seed(job["seed"], rep)
            if int(t["seed"]) != expected_seed:
                raise ValueError(f"trajectory {rep} seed {t['seed']} != derived {expected_seed}")
            compact = {k: v for k, v in t.items() if k not in ("history", "bestPermutation", "neighborBiasHistory")}
            compact["improvementHistory"] = [{k: h[k] for k in ("iteration", "cost", "elapsedMs", "move")} for h in t["history"]]
            out[rep] = {
                "value": t["bestCost"],
                "initial_value": t["initialCost"],
                "permutation": t["bestPermutation"],
                "method_time_s": t["elapsedMs"] / 1000.0,
                "time_to_best_s": t["timeToBestMs"] / 1000.0,
                "iterations": t["iterations"],
                "improvements": t["acceptedMoves"],
                "stop_reason": t["stopReason"],
                "method_time_limit_reached": t["stopReason"] == "maxTime",
                # Only cache misses start an LKH process; lkhTimeMs also counts
                # cache lookups and file I/O.
                "tsp": {"calls": t["lkhCacheMisses"], "time_s": t["lkhTimeMs"] / 1000.0, "time_limit_hits": t["lkhTimeLimitHits"], "invalid_tours": 0},
                "derived_seeds": {
                    "trajectory_index": rep,
                    "trajectory_seed": t["seed"],
                    "trajectory_seed_rule": "cbm_derive_seed(run_seed, trajectory_index)",
                    "lkh_seed_rule": "cbm_derive_seed(run_seed, hash(sorted sub-problem columns)); see CBMLKH::applyLKH",
                },
                "method_stats": {"trajectory": compact, "execution": execution},
            }
        return out
    raise ValueError(method)
