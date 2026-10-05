"""Experiment design: the immutable list of jobs and its fingerprint.

A *repetition* is one independent run of a method on an instance (one row of
the detailed CSV). A *job* is one process launch: exactly one repetition for
ENS/ILS/LKH, and a batch of consecutive repetitions (one per thread) for
CBMLKH. The plan is fixed when an experiment directory is created; resuming
never re-plans.
"""

from __future__ import annotations

import copy
import hashlib
import json
from typing import Any, Dict, List

from .instances import InstanceInfo
from .methods import DEFAULT_PARAMS, LOCKED_PARAMS, METHOD_IDS, cbmlkh_lkh_max_time, estimate_resources
from .seeds import run_seed

DEFAULT_REPETITIONS = 10
DEFAULT_SOFT_LIMIT_S = 2 * 3600.0
DEFAULT_HARD_LIMIT_S = 5 * 3600.0
DEFAULT_ROOT_SEED = 20261005


def make_spec(
    root_seed: int,
    methods: List[str],
    repetitions: Dict[str, int],
    soft_limit_s: float,
    hard_limit_s: float,
    cbmlkh_threads: int,
    param_overrides: Dict[str, Dict[str, Any]],
    instances: List[InstanceInfo],
    instance_dir: str,
    all_instances_selected: bool,
    allow_locked_overrides: bool = False,
) -> dict:
    if not methods:
        raise ValueError("no methods selected")
    for m in methods:
        if m not in METHOD_IDS:
            raise ValueError(f"unknown method '{m}' (choose from {METHOD_IDS})")
    if not 0 < soft_limit_s <= hard_limit_s:
        raise ValueError("need 0 < soft limit <= hard limit")
    if cbmlkh_threads < 1:
        raise ValueError("cbmlkh threads must be >= 1")

    params = {}
    for m in methods:
        merged = copy.deepcopy(DEFAULT_PARAMS[m])
        if m == "cbmlkh":
            merged["lkhMaxTime"] = cbmlkh_lkh_max_time(soft_limit_s)
        for key, value in param_overrides.get(m, {}).items():
            locked = LOCKED_PARAMS.get(m, set())
            if (locked == "all" or key in locked) and not allow_locked_overrides:
                raise ValueError(f"{m} parameter '{key}' is fixed (tuned / project-wide) and cannot be overridden")
            merged[key] = value
        params[m] = merged

    spec = {
        "root_seed": int(root_seed),
        "methods": list(methods),
        "repetitions": {m: int(repetitions[m]) for m in methods},
        "soft_limit_s": float(soft_limit_s),
        "hard_limit_s": float(hard_limit_s),
        "cbmlkh_threads": int(cbmlkh_threads),
        "method_params": params,
        "instance_dir": instance_dir,
        "instances": [i.as_dict() for i in instances],
    }
    if allow_locked_overrides:
        spec["locked_params_overridden"] = True  # testing only; never a final design
    spec["partial"] = not (
        all_instances_selected
        and sorted(methods) == sorted(METHOD_IDS)
        and all(spec["repetitions"][m] == DEFAULT_REPETITIONS for m in methods)
        and soft_limit_s == DEFAULT_SOFT_LIMIT_S
        and hard_limit_s == DEFAULT_HARD_LIMIT_S
        and not any(param_overrides.get(m) for m in methods)
        and not allow_locked_overrides
    )
    return spec


def fingerprint(spec: dict) -> str:
    """Hash of everything that defines the experiment's results.

    Machine-specific facts (absolute instance paths, binary locations, CPU
    budget) are excluded so an experiment can be moved and resumed elsewhere;
    instance *contents* are included through their SHA-256.
    """
    canonical = {k: v for k, v in spec.items() if k not in ("instance_dir", "instances")}
    canonical["instances"] = [(i["name"], i["sha256"]) for i in spec["instances"]]
    blob = json.dumps(canonical, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode()).hexdigest()


def build_jobs(spec: dict) -> List[dict]:
    """Deterministic job list. CBMLKH repetitions are packed `cbmlkh_threads` per job."""
    jobs = []
    for inst in spec["instances"]:
        for method in spec["methods"]:
            reps = spec["repetitions"][method]
            params = spec["method_params"][method]
            if method == "cbmlkh":
                # One seed for all of an instance's trajectories: trajectory r
                # derives cbm_derive_seed(seed, r), so results do not depend on
                # how repetitions are split into executions.
                seed = run_seed(spec["root_seed"], method, inst["name"], -1)
                width = spec["cbmlkh_threads"]
                for start in range(0, reps, width):
                    batch = list(range(start, min(reps, start + width)))
                    jobs.append(_job(method, inst, batch, seed, len(batch), params))
            else:
                for rep in range(reps):
                    jobs.append(_job(method, inst, [rep], run_seed(spec["root_seed"], method, inst["name"], rep), 1, params))
    return jobs


def _job(method: str, inst: dict, reps: List[int], seed: int, threads: int, params: dict) -> dict:
    est = estimate_resources(method, inst["rows"], inst["cols"], params, threads)
    if len(reps) == 1:
        rep_tag = f"r{reps[0]:02d}"
    else:
        rep_tag = f"r{reps[0]:02d}-{reps[-1]:02d}"
    return {
        "job_id": f"{method}__{inst['name']}__{rep_tag}",
        "method": method,
        "instance": inst,
        "repetitions": reps,
        "seed": seed,
        "threads": threads,
        "cpus": threads,
        "est_mem_bytes": est.mem_bytes,
        "est_disk_bytes": est.disk_bytes,
    }


def run_id(method: str, instance: str, rep: int) -> str:
    return f"{method}__{instance}__r{rep:02d}"
