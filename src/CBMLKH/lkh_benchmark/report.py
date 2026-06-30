"""Join per-solver result stores into a side-by-side comparison.

Reads the JSON results written by each solver's :class:`ResultStore` and builds
one row per instance with each solver's cost and runtime, the winner (lowest
cost) and the relative gap, so pure LKH and CBMLKH can be compared directly.
"""

from __future__ import annotations

import csv
from pathlib import Path

from .result_store import ResultStore


def _load_by_instance(results_dir: str | Path, run_tag: str) -> dict[str, dict]:
    store = ResultStore(results_dir, run_tag)
    return {rec["instance"]: rec for rec in store.load_all()}


def build_comparison(results_dir: str | Path, labeled_tags: dict[str, str]) -> list[dict]:
    """Build comparison rows across solvers.

    Args:
        results_dir: directory holding the per-solver ``<run_tag>/`` folders.
        labeled_tags: maps a short label (e.g. ``"lkh"``) to its ``run_tag``.

    Returns:
        One dict per instance with ``cost_<label>`` / ``timeSec_<label>`` columns,
        plus ``best`` (winning label), ``bestCost`` and ``gapPct`` (worst vs best).
    """
    per_label = {label: _load_by_instance(results_dir, tag) for label, tag in labeled_tags.items()}

    instances = sorted({name for recs in per_label.values() for name in recs})
    labels = list(labeled_tags)

    rows = []
    for name in instances:
        row: dict = {"instance": name}
        costs: dict[str, float] = {}
        for label in labels:
            rec = per_label[label].get(name)
            cost = rec.get("cost") if rec else None
            row[f"cost_{label}"] = cost
            row[f"timeSec_{label}"] = rec.get("runtimeSec") if rec else None
            if cost is not None:
                costs[label] = cost
            # carry rows/cols from whichever solver has them
            if rec and "rows" in rec:
                row.setdefault("rows", rec["rows"])
                row.setdefault("cols", rec["cols"])

        if costs:
            best_label = min(costs, key=costs.get)
            best_cost = costs[best_label]
            worst_cost = max(costs.values())
            row["best"] = best_label
            row["bestCost"] = best_cost
            row["gapPct"] = round(100.0 * (worst_cost - best_cost) / best_cost, 4) if best_cost else None
        rows.append(row)

    return rows


def write_comparison_csv(rows: list[dict], path: str | Path) -> Path:
    path = Path(path)
    if not rows:
        path.write_text("")
        return path
    # Union of keys, preserving a sensible leading order.
    leading = ["instance", "rows", "cols"]
    keys = leading + [k for k in dict.fromkeys(k for r in rows for k in r) if k not in leading]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


def print_comparison(rows: list[dict], labels: list[str]) -> None:
    """Pretty-print a compact comparison table to stdout."""
    if not rows:
        print("[lkh_benchmark] no comparison rows.")
        return

    headers = ["instance", "cols"]
    for label in labels:
        headers += [f"cost[{label}]", f"t[{label}]s"]
    headers += ["best", "gap%"]

    def fmt(value) -> str:
        return "" if value is None else str(value)

    table = [headers]
    for r in rows:
        line = [fmt(r.get("instance")), fmt(r.get("cols"))]
        for label in labels:
            line += [fmt(r.get(f"cost_{label}")), fmt(r.get(f"timeSec_{label}"))]
        line += [fmt(r.get("best")), fmt(r.get("gapPct"))]
        table.append(line)

    widths = [max(len(row[i]) for row in table) for i in range(len(headers))]
    for row in table:
        print("  ".join(cell.ljust(widths[i]) for i, cell in enumerate(row)))
