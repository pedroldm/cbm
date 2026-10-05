#!/usr/bin/env python3
"""Stand-in for ENS / ILS / lkh_standalone / CBMLKH in runner tests.

Speaks each binary's command line, writes a valid native JSON (a seed-dependent
permutation with its true block count) and logs what it did to
$FAKE_SOLVER_LOG/<pid>.json. Behaviour per instance comes from the JSON file
$FAKE_SOLVER_CONTROL: {"<instance name>": {"sleep": s, "mode": ..., "grandchild": s}}
with mode one of ok, fail, bad_output, wrong_blocks, ignore_term.
"""

import json
import os
import random
import signal
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cbm_experiments.instances import CBMInstance  # noqa: E402
from cbm_experiments.seeds import cbm_derive_seed  # noqa: E402


def parse() -> dict:
    args = sys.argv[1:]
    if len(args) == 1 and not args[0].startswith("--"):  # CBMLKH: a config file
        cfg = dict(line.split("=", 1) for line in Path(args[0]).read_text().splitlines() if "=" in line)
        return {"kind": "cbmlkh", "instance": cfg["instancePath"], "out": cfg["outputPath"], "seed": int(cfg["seed"]),
                "threads": int(cfg["threads"]), "offset": int(cfg["trajectoryOffset"]), "workdir": cfg["lkhTmpDir"]}
    opts = dict(a[2:].split("=", 1) for a in args)
    kind = "lkh" if "lkhPath" in opts else ("ils" if "phi" in opts else "ens")
    return {"kind": kind, "instance": opts["filePath"], "out": opts["outputPath"], "seed": int(opts["seed"]), "workdir": opts["workDir"]}


def permutation(n: int, seed: int) -> list:
    perm = list(range(1, n + 1))
    random.Random(seed).shuffle(perm)
    return perm


def main() -> int:
    a = parse()
    name = Path(a["instance"]).name
    control = json.loads(Path(os.environ["FAKE_SOLVER_CONTROL"]).read_text()).get(name, {}) if os.environ.get("FAKE_SOLVER_CONTROL") else {}
    mode = control.get("mode", "ok")
    log = {"pid": os.getpid(), "pgid": os.getpgid(0), "kind": a["kind"], "instance": name, "seed": a["seed"], "start": time.time(),
           "affinity": sorted(os.sched_getaffinity(0)), "workdir": a["workdir"], "cwd": os.getcwd(), "omp": os.environ.get("OMP_NUM_THREADS")}
    log_dir = Path(os.environ["FAKE_SOLVER_LOG"]) if os.environ.get("FAKE_SOLVER_LOG") else None

    def write_log(**extra) -> None:
        if log_dir:
            log.update(extra)
            tmp = log_dir / f".{os.getpid()}.tmp"
            tmp.write_text(json.dumps(log))
            os.replace(tmp, log_dir / f"{os.getpid()}.json")

    write_log()
    if mode == "ignore_term":
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
    if control.get("grandchild"):
        child = subprocess.Popen(["sleep", str(control["grandchild"])])
        write_log(grandchild=child.pid)
    # Scratch file in the private work dir, to check isolation.
    Path(a["workdir"]).mkdir(parents=True, exist_ok=True)
    (Path(a["workdir"]) / "scratch.txt").write_text(f"{name} {a['seed']}\n")
    time.sleep(float(control.get("sleep", 0.2)))

    if mode == "fail":
        print("fake solver failure", file=sys.stderr)
        return 3
    if mode == "bad_output":
        Path(a["out"]).write_text("{not json")
        return 0

    inst = CBMInstance(Path(a["instance"]))
    if a["kind"] == "cbmlkh":
        trajectories = []
        for g in range(a["offset"], a["offset"] + a["threads"]):
            tseed = cbm_derive_seed(a["seed"], g)
            perm = permutation(inst.n_cols, tseed)
            cost = inst.count_blocks(perm)
            trajectories.append({"index": g, "seed": tseed, "bestCost": cost, "initialCost": cost + 1, "bestPermutation": perm,
                                 "elapsedMs": 100, "timeToBestMs": 50, "iterations": 3, "acceptedMoves": 1, "stopReason": "maxIterations",
                                 "lkhCacheMisses": 2, "lkhCalls": 3, "lkhTimeMs": 10, "lkhTimeLimitHits": 0,
                                 "history": [{"iteration": 1, "cost": cost, "elapsedMs": 50, "move": "MERGE", "histogram": []}],
                                 "neighborBiasHistory": [1.0]})
        native = {"config": {"seed": a["seed"], "threads": a["threads"], "trajectoryOffset": a["offset"]},
                  "global": {"runtimeMs": 100, "bestCost": min(t["bestCost"] for t in trajectories), "lkhCache": {"hits": 1, "misses": 2}},
                  "trajectories": trajectories}
    else:
        perm = permutation(inst.n_cols, a["seed"])
        blocks = inst.count_blocks(perm) + (1 if mode == "wrong_blocks" else 0)
        native = {"seed": a["seed"], "best_blocks": blocks, "initial_blocks": blocks, "permutation": perm, "elapsed_s": 0.2,
                  "time_to_best_s": 0.1, "iterations_completed": 2, "history": [{"iteration": 0, "blocks": blocks, "elapsed_s": 0.1}],
                  "stop_reason": "iterations", "time_limit_reached": False, "parameters": {}, "validated": True,
                  "tsp": {"calls": 3, "time_s": 0.1, "time_limit_hits": 0, "invalid_tours": 0, "seeds": [1, 2, 3]}, "perturbation_seed": 7,
                  "tsp_build_time_s": 0.01, "lkh_time_s": 0.1, "lkh_preprocessing_time_s": 0.01, "lkh_time_to_best_s": 0.05, "lkh_improvements": 1}
    tmp = a["out"] + ".tmp"
    Path(tmp).write_text(json.dumps(native))
    os.replace(tmp, a["out"])
    write_log(end=time.time())
    return 0


if __name__ == "__main__":
    sys.exit(main())
