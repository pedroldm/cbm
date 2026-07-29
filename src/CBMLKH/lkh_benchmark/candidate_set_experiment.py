"""Investigate the impact of *not* using the 1-tree on LKH's startup.

LKH's default candidate set (``CANDIDATE_SET_TYPE = ALPHA``) builds a minimum
1-tree and runs a subgradient ascent (Held--Karp) to derive alpha-nearness — the
dominant *startup* cost before local search even begins. This script sweeps the
full cross product of ``CANDIDATE_SET_TYPE`` × ``SUBGRADIENT`` ∈ {YES, NO} over a
small, representative subset of instances and reports, per (instance, config):

  * ``preprocess_sec`` — LKH's ``Preprocessing time`` (candidate build + ascent),
    the quantity the 1-tree dominates;
  * ``ascent_sec`` / ``lower_bound`` — the 1-tree subgradient ascent portion
    specifically (absent when no 1-tree is built, e.g. NEAREST-NEIGHBOR, or when
    ``SUBGRADIENT = NO``);
  * ``cbm_cost`` — CBM 1-block count of the final tour (quality the cheaper
    candidate set buys or costs), cross-checked two independent ways (``valid``).

The candidate-set method and the subgradient toggle are the only parameters that
vary; every other LKH knob is left at its default. The TSP is the same
``EXPLICIT`` / ``FULL_MATRIX`` Hamming problem the production ``LKHWrapper``
builds (reused via :class:`TSPConverter`), so coordinate-only candidate sets
(DELAUNAY, QUADRANT) are *illegal* here — the script records LKH's rejection
rather than crashing, which is itself part of the finding.

Run from ``src/CBMLKH``::

    python -m lkh_benchmark.candidate_set_experiment

Everything tunable lives in the PARAMETERS block below.
"""

from __future__ import annotations

import csv
import re
import subprocess
import time
import uuid
from pathlib import Path

from .config import CBMLKH_DIR, INSTANCES_DIR, LKH_PATH
from .instance import CBMInstance
from .lkh_runner import LKHRunner
from .tsp_converter import TSPConverter

# =========================================================================== #
#  PARAMETERS — edit here.
# =========================================================================== #

# Representative subset (restricted to the a--i, scp* and cbm_*.txt families).
# Chosen to span the TSP node count (= column count), which is what drives the
# 1-tree startup cost, from ~200 up to ~11000, across all three families, plus a
# cbm density sweep at fixed size. Comment out the heavy tail for a quick pass.
INSTANCES = [
    # a--i family (100 rows; 200 / 500 / 1000 columns)
    "a1",   # 100 x 200
    "d1",   # 100 x 500
    "g1",   # 100 x 1000
    # scp family (50 -> 10000 columns)
    "scpe1",  # 50 x 500
    "scpa1",  # 300 x 3000
    "scpd1",  # 400 x 4000
    "scpf1",  # 500 x 5000
    "scph1",  # 1000 x 10000  (heavy)
    # cbm density sweep at fixed size (1100 x 11000, density 2/5/10/15) (heavy)
    "cbm_1100_11000_2.txt",
    "cbm_1100_11000_5.txt",
    "cbm_1100_11000_10.txt",
    "cbm_1100_11000_15.txt",
]

# The two parameters under study. Every combination of the two is run; every
# other LKH knob is left at its default.
#   - ALPHA / POPMUSIC build a 1-tree; NEAREST-NEIGHBOR does not.
#   - DELAUNAY / QUADRANT are coordinate-only and illegal for EXPLICIT weights
#     (kept so the CSV documents their inapplicability empirically).
CANDIDATE_SET_TYPES = ["ALPHA", "NEAREST-NEIGHBOR", "POPMUSIC", "DELAUNAY", "QUADRANT"]
SUBGRADIENT_OPTIONS = ["YES", "NO"]  # YES is the LKH default

# Full cross product -> {"alpha/sg": {...}, "alpha/nosg": {...}, ...}.
CANDIDATE_CONFIGS: dict[str, dict[str, str]] = {
    f"{cst.lower()}/{'sg' if sg == 'YES' else 'nosg'}": {"CANDIDATE_SET_TYPE": cst, "SUBGRADIENT": sg}
    for cst in CANDIDATE_SET_TYPES
    for sg in SUBGRADIENT_OPTIONS
}

# Non-default LKH knobs: a per-run time budget and a single run so the sweep
# terminates (one startup + one search per config, as production uses).
TIME_LIMIT = 3600  # seconds per (instance, config)
RUNS = 1

WORK_DIR = Path("/tmp/LKH_candset")
RESULTS_DIR = CBMLKH_DIR / "lkh_benchmark" / "results" / "candidate_set"
CSV_NAME = "candidate_set_experiment.csv"

# =========================================================================== #
#  Internals
# =========================================================================== #

_RE_PREPROCESS = re.compile(r"Preprocessing time\s*=\s*([\d.]+)\s*sec")
_RE_ASCENT = re.compile(r"Lower bound\s*=\s*([\d.]+),\s*Ascent time\s*=\s*([\d.]+)\s*sec")
_RE_COST_MIN = re.compile(r"Cost\.min\s*=\s*(\d+)")
_RE_RUN_COST = re.compile(r"Run \d+:\s*Cost\s*=\s*(\d+),\s*Time\s*=\s*([\d.]+)\s*sec")
_RE_ERROR = re.compile(r"\*\*\* Error \*\*\*\s*\n?\s*(.*)")

CSV_FIELDS = [
    "instance", "rows", "cols", "config", "status",
    "preprocess_sec", "ascent_sec", "lower_bound",
    "cbm_cost", "cbm_blocks", "valid", "lkh_tour_cost", "wall_sec", "note",
]


def _par_lines(problem_file: str, tour_file: str, extra: dict[str, str]) -> list[str]:
    # Only PROBLEM_FILE/TOUR_FILE (mandatory I/O), the time budget, and the
    # candidate-set line under study. Everything else stays at LKH's default.
    lines = [
        f"PROBLEM_FILE = {problem_file}",
        f"TOUR_FILE = {tour_file}",
        f"TIME_LIMIT = {TIME_LIMIT}",
        f"RUNS = {RUNS}",
    ]
    lines += [f"{k} = {v}" for k, v in extra.items()]
    return lines


def _parse_stdout(out: str) -> dict:
    prep = _RE_PREPROCESS.search(out)
    asc = _RE_ASCENT.search(out)
    # Best across all runs (RUNS is left at LKH's default), falling back to the
    # last per-run line if the summary is absent.
    cmin = _RE_COST_MIN.search(out)
    runs = _RE_RUN_COST.findall(out)
    return {
        "preprocess_sec": float(prep.group(1)) if prep else None,
        "lower_bound": float(asc.group(1)) if asc else None,
        "ascent_sec": float(asc.group(2)) if asc else None,  # None => no 1-tree built
        "lkh_tour_cost": int(cmin.group(1)) if cmin else (int(runs[-1][0]) if runs else None),
    }


def _run_one(lkh: str, tsp_file: Path, instance: CBMInstance, label: str, extra: dict[str, str]) -> dict:
    exec_id = uuid.uuid4().hex[:8]
    safe = label.replace("/", "-")  # labels contain '/', not filesystem-safe
    par_file = WORK_DIR / f"{instance.name}_{safe}_{exec_id}.par"
    tour_file = WORK_DIR / f"{instance.name}_{safe}_{exec_id}.tour"
    par_file.write_text("\n".join(_par_lines(str(tsp_file), str(tour_file), extra)) + "\n")

    t0 = time.perf_counter()
    proc = subprocess.run([lkh, str(par_file)], capture_output=True, text=True)
    wall = time.perf_counter() - t0

    rec = {
        "instance": instance.name, "rows": instance.rows, "cols": instance.cols,
        "config": label, "wall_sec": round(wall, 3),
        "preprocess_sec": None, "ascent_sec": None, "lower_bound": None,
        "cbm_cost": None, "cbm_blocks": None, "valid": None, "lkh_tour_cost": None,
        "status": None, "note": "",
    }

    if not tour_file.exists() or proc.returncode != 0:
        err = _RE_ERROR.search(proc.stdout + "\n" + proc.stderr)
        rec["status"] = "failed"
        rec["note"] = (err.group(1).strip() if err else f"rc={proc.returncode}, no tour")[:160]
        _cleanup(par_file, tour_file)
        return rec

    rec.update(_parse_stdout(proc.stdout))
    # LKH's tour cost is a symmetric-Hamming proxy, not the CBM objective. Score
    # quality by the CBM 1-block count of the returned permutation, computed two
    # independent ways (path_cost formula vs. direct recount) as a self-check.
    perm = LKHRunner._parse_tour(tour_file)
    if len(perm) == instance.cols and sorted(perm) == list(range(instance.cols)):
        rec["cbm_cost"] = instance.path_cost(perm)
        rec["cbm_blocks"] = instance.count_one_blocks(perm)
        rec["valid"] = rec["cbm_cost"] == rec["cbm_blocks"]
    else:
        rec["valid"] = False
        rec["note"] = f"tour not a permutation of {instance.cols} columns"
    rec["status"] = "ok"
    _cleanup(par_file, tour_file)
    return rec


def _cleanup(*paths: Path) -> None:
    for p in paths:
        p.unlink(missing_ok=True)


def main() -> None:
    lkh = LKH_PATH
    if not Path(lkh).exists():
        raise FileNotFoundError(f"LKH executable not found: {lkh}")
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    instances_dir = Path(INSTANCES_DIR)

    rows: list[dict] = []
    for name in INSTANCES:
        path = instances_dir / name
        if not path.exists():
            print(f"[candset] SKIP missing instance: {name}")
            continue
        print(f"[candset] loading {name} ...", flush=True)
        instance = CBMInstance.from_file(path)
        tsp_file = WORK_DIR / f"{name}.tsp"
        # Build the (cols+1)^2 Hamming matrix once; reuse across all configs.
        t0 = time.perf_counter()
        TSPConverter(instance).write(tsp_file, name)
        print(f"[candset]   {name}: {instance.rows}x{instance.cols}, tsp built in {time.perf_counter()-t0:.1f}s", flush=True)

        for label, extra in CANDIDATE_CONFIGS.items():
            rec = _run_one(lkh, tsp_file, instance, label, extra)
            rows.append(rec)
            print(
                f"[candset]   {name:<26} {label:<22} {rec['status']:<7} "
                f"prep={rec['preprocess_sec']} ascent={rec['ascent_sec']} "
                f"cbm_cost={rec['cbm_cost']} valid={rec['valid']} {rec['note']}",
                flush=True,
            )
        tsp_file.unlink(missing_ok=True)

    csv_path = RESULTS_DIR / CSV_NAME
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"\n[candset] wrote {len(rows)} rows -> {csv_path}")
    _print_summary(rows)


def _print_summary(rows: list[dict]) -> None:
    """Per-instance table: preprocessing time and cost for each config."""
    labels = list(CANDIDATE_CONFIGS)
    by_inst: dict[str, dict[str, dict]] = {}
    for r in rows:
        by_inst.setdefault(r["instance"], {})[r["config"]] = r

    def cell(r: dict | None) -> str:
        if r is None or r["status"] != "ok":
            return "n/a"
        return f"{(r['preprocess_sec'] or 0):.2f}/{r['cbm_cost']}"

    # Width each column to the widest of its label and its values.
    widths = {l: max(len(l), max((len(cell(cfgs.get(l))) for cfgs in by_inst.values()), default=0)) + 2 for l in labels}

    print("\n=== preprocess_sec / cbm_cost  ('/sg' = SUBGRADIENT on, '/nosg' = off) ===")
    head = f"{'instance':<26}{'cols':>7}  " + "".join(f"{l:<{widths[l]}}" for l in labels)
    print(head)
    print("-" * len(head))
    for inst, cfgs in by_inst.items():
        cols = next(iter(cfgs.values()))["cols"]
        line = f"{inst:<26}{cols:>7}  " + "".join(f"{cell(cfgs.get(l)):<{widths[l]}}" for l in labels)
        print(line)


if __name__ == "__main__":
    main()
