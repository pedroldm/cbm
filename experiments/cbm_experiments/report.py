"""CSV reports, regenerated from the persisted records (idempotent, read-only on records).

Files written to ``<experiment>/results/``:

detailed.csv
    One row per planned repetition, whatever its status (pending rows included,
    so the file doubles as a progress view).
summary.csv
    One row per (instance, method). Objective and time statistics use only
    repetitions with status ``completed`` (always validated); every other
    status is counted in its own column and never contributes a value.
comparison.csv
    One row per instance: the best observed objective over all methods, and
    each method's best and mean value and their gaps to it.
methods.csv
    One row per method, aggregated over instances.

Objective: number of 1-blocks, **lower is better** (``sense = minimize`` in
every record). gap_pct(x) = 100 * (x - ref) / ref, where ref is the best value
observed for the instance by any method in this experiment (or the reference
value given with --reference, when that is lower).
"""

from __future__ import annotations

import csv
import io
import math
import statistics
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from .plan import run_id
from .state import JOB_STATUSES, PENDING
from .util import atomic_write_text, read_json

DETAILED_FIELDS = [
    "experiment", "instance", "rows", "cols", "method", "tsp_backend", "repetition", "run_id", "job_id", "seed", "status",
    "objective", "initial_objective", "validated", "wall_time_s", "method_time_s", "time_to_best_s", "soft_limit_reached",
    "method_time_limit_reached", "stop_reason", "iterations", "improvements", "tsp_calls", "tsp_time_s", "tsp_time_limit_hits",
    "invalid_tours", "attempt", "exit_code", "signal", "error", "best_observed", "gap_to_best_observed_pct", "reference_value",
    "gap_to_reference_pct", "record_path",
]


def _stats(values: List[float], prefix: str) -> Dict[str, Optional[float]]:
    if not values:
        return {f"{prefix}_{k}": None for k in ("mean", "median", "std", "min", "max")}
    return {
        f"{prefix}_mean": statistics.fmean(values),
        f"{prefix}_median": statistics.median(values),
        f"{prefix}_std": statistics.stdev(values) if len(values) > 1 else None,
        f"{prefix}_min": min(values),
        f"{prefix}_max": max(values),
    }


def _gap(value: Optional[float], ref: Optional[float]) -> Optional[float]:
    if value is None or ref is None or ref <= 0:
        return None
    return 100.0 * (value - ref) / ref


def _fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return "" if math.isnan(value) else f"{value:.6g}"
    return str(value)


def _write_csv(path: Path, fields: List[str], rows: Iterable[dict]) -> None:
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({k: _fmt(row.get(k)) for k in fields})
    atomic_write_text(path, buf.getvalue())


def load_reference(path: Optional[Path]) -> Dict[str, float]:
    """Optional CSV with columns instance,value (e.g. best known solutions)."""
    if not path:
        return {}
    with open(path, newline="") as f:
        return {row["instance"]: float(row["value"]) for row in csv.DictReader(f)}


class RecordCache:
    """Avoids re-reading unchanged record files on every regeneration."""

    def __init__(self):
        self._cache: Dict[str, tuple] = {}

    def load(self, path: Path) -> Optional[dict]:
        try:
            st = path.stat()
        except OSError:
            return None
        key = (st.st_mtime_ns, st.st_size)
        hit = self._cache.get(str(path))
        if hit and hit[0] == key:
            return hit[1]
        data = read_json(path)
        self._cache[str(path)] = (key, data)
        return data


def collect_rows(exp_dir: Path, experiment: dict, cache: Optional[RecordCache] = None) -> List[dict]:
    cache = cache or RecordCache()
    rows = []
    for job in experiment["jobs"]:
        job_dir = exp_dir / "runs" / job["job_id"]
        state = read_json(job_dir / "state.json") or {"status": PENDING}
        for rep in job["repetitions"]:
            rid = run_id(job["method"], job["instance"]["name"], rep)
            rec_path = job_dir / f"{rid}.json"
            rec = cache.load(rec_path) if state.get("status") not in (PENDING, "running") else None
            # job state is authoritative: a record left by an attempt that never
            # finalized (e.g. power loss mid-write) does not count.
            if rec and rec.get("status") != state.get("status"):
                rec = None
            row = {
                "experiment": experiment["name"],
                "instance": job["instance"]["name"],
                "rows": job["instance"]["rows"],
                "cols": job["instance"]["cols"],
                "method": job["method"],
                "tsp_backend": experiment["spec"]["method_params"][job["method"]]["tsp_backend"],
                "repetition": rep,
                "run_id": rid,
                "job_id": job["job_id"],
                "seed": job["seed"],
                "status": state.get("status", PENDING),
            }
            if rec:
                obj, tim, srch, tsp, ex = rec["objective"], rec["timing"], rec["search"], rec["tsp_solver"], rec["execution"]
                row.update(
                    status=rec["status"],
                    objective=obj["value"],
                    initial_objective=obj["initial_value"],
                    validated=obj["validated"],
                    wall_time_s=tim["wall_time_s"],
                    method_time_s=tim["method_time_s"],
                    time_to_best_s=tim["time_to_best_s"],
                    soft_limit_reached=tim["soft_limit_reached"],
                    method_time_limit_reached=rec["stopping"]["method_time_limit_reached"],
                    stop_reason=rec["stopping"]["reason"],
                    iterations=srch["iterations"],
                    improvements=srch["improvements"],
                    tsp_calls=tsp["calls"],
                    tsp_time_s=tsp["time_s"],
                    tsp_time_limit_hits=tsp["time_limit_hits"],
                    invalid_tours=tsp["invalid_tours"],
                    attempt=ex["attempt"],
                    exit_code=ex["exit_code"],
                    signal=ex["signal"],
                    error=(rec["error"] or {}).get("message"),
                    record_path=str(rec_path.relative_to(exp_dir)),
                )
            rows.append(row)
    return rows


def generate(exp_dir: Path, reference: Optional[Dict[str, float]] = None, cache: Optional[RecordCache] = None) -> Dict[str, Path]:
    exp_dir = Path(exp_dir)
    experiment = read_json(exp_dir / "experiment.json")
    if experiment is None:
        raise FileNotFoundError(f"{exp_dir} has no experiment.json")
    reference = reference or {}
    rows = collect_rows(exp_dir, experiment, cache)
    methods = experiment["spec"]["methods"]
    instances = [i["name"] for i in experiment["spec"]["instances"]]

    def ok(r: dict) -> bool:
        return r["status"] == "completed" and r.get("objective") is not None

    best = {}
    for r in rows:
        if ok(r):
            best[r["instance"]] = min(best.get(r["instance"], math.inf), r["objective"])
    for r in rows:
        ref = reference.get(r["instance"])
        r["best_observed"] = best.get(r["instance"])
        r["reference_value"] = ref
        if ok(r):
            r["gap_to_best_observed_pct"] = _gap(r["objective"], best.get(r["instance"]))
            r["gap_to_reference_pct"] = _gap(r["objective"], ref)

    out_dir = exp_dir / "results"
    paths = {"detailed": out_dir / "detailed.csv", "summary": out_dir / "summary.csv", "comparison": out_dir / "comparison.csv", "methods": out_dir / "methods.csv"}
    _write_csv(paths["detailed"], DETAILED_FIELDS, rows)

    # --- summary: instance x method
    by_pair: Dict[tuple, List[dict]] = {}
    for r in rows:
        by_pair.setdefault((r["instance"], r["method"]), []).append(r)
    summary_rows = []
    for inst in instances:
        ref = min([v for v in (best.get(inst), reference.get(inst)) if v is not None], default=None)
        for m in methods:
            group = by_pair.get((inst, m), [])
            good = [r for r in group if ok(r)]
            values = [float(r["objective"]) for r in good]
            row = {"instance": inst, "method": m, "planned": len(group)}
            for s in JOB_STATUSES:
                row[f"n_{s}"] = sum(1 for r in group if r["status"] == s)
            row["n_used"] = len(good)
            row.update(_stats(values, "objective"))
            row.update(_stats([float(r["wall_time_s"]) for r in good if r.get("wall_time_s") is not None], "wall_time_s"))
            row.update(_stats([float(r["time_to_best_s"]) for r in good if r.get("time_to_best_s") is not None], "time_to_best_s"))
            row["n_soft_limit_reached"] = sum(1 for r in group if r.get("soft_limit_reached"))
            row["best_observed"] = best.get(inst)
            row["reference_value"] = reference.get(inst)
            row["gap_min_pct"] = _gap(row["objective_min"], ref)
            row["gap_mean_pct"] = _gap(row["objective_mean"], ref)
            row["hits_best_observed"] = sum(1 for v in values if best.get(inst) is not None and v == best[inst])
            summary_rows.append(row)
    stat_cols = [f"{p}_{k}" for p in ("objective", "wall_time_s", "time_to_best_s") for k in ("mean", "median", "std", "min", "max")]
    summary_fields = (
        ["instance", "method", "planned"] + [f"n_{s}" for s in JOB_STATUSES] + ["n_used"] + stat_cols
        + ["n_soft_limit_reached", "best_observed", "reference_value", "gap_min_pct", "gap_mean_pct", "hits_best_observed"]
    )
    _write_csv(paths["summary"], summary_fields, summary_rows)

    # --- comparison: one row per instance
    summary_index = {(r["instance"], r["method"]): r for r in summary_rows}
    comp_rows, comp_fields = [], ["instance", "best_observed", "best_methods", "reference_value"]
    for m in methods:
        comp_fields += [f"{m}_min", f"{m}_mean", f"{m}_gap_min_pct", f"{m}_gap_mean_pct", f"{m}_n_used"]
    for inst in instances:
        row = {"instance": inst, "best_observed": best.get(inst), "reference_value": reference.get(inst)}
        row["best_methods"] = ";".join(m for m in methods if best.get(inst) is not None and summary_index[(inst, m)]["objective_min"] == best[inst])
        for m in methods:
            s = summary_index[(inst, m)]
            row.update({f"{m}_min": s["objective_min"], f"{m}_mean": s["objective_mean"], f"{m}_gap_min_pct": s["gap_min_pct"],
                        f"{m}_gap_mean_pct": s["gap_mean_pct"], f"{m}_n_used": s["n_used"]})
        comp_rows.append(row)
    _write_csv(paths["comparison"], comp_fields, comp_rows)

    # --- per method over instances
    method_rows = []
    for m in methods:
        srows = [summary_index[(i, m)] for i in instances]
        mean_gaps = [s["gap_mean_pct"] for s in srows if s["gap_mean_pct"] is not None]
        min_gaps = [s["gap_min_pct"] for s in srows if s["gap_min_pct"] is not None]
        row = {
            "method": m,
            "instances": len(instances),
            "planned": sum(s["planned"] for s in srows),
            **{f"n_{st}": sum(s[f"n_{st}"] for s in srows) for st in JOB_STATUSES},
            "instances_with_result": sum(1 for s in srows if s["n_used"] > 0),
            "instances_best_observed": sum(1 for i in instances if m in comp_rows[instances.index(i)]["best_methods"].split(";")),
            "mean_of_gap_mean_pct": statistics.fmean(mean_gaps) if mean_gaps else None,
            "mean_of_gap_min_pct": statistics.fmean(min_gaps) if min_gaps else None,
        }
        method_rows.append(row)
    _write_csv(paths["methods"], ["method", "instances", "planned"] + [f"n_{s}" for s in JOB_STATUSES]
               + ["instances_with_result", "instances_best_observed", "mean_of_gap_mean_pct", "mean_of_gap_min_pct"], method_rows)
    return paths
