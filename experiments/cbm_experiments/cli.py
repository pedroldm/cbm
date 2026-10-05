"""Command line: ``python experiments/run_experiments.py {run,plan,status,report} ...``"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

from . import report
from .instances import discover_instances
from .methods import METHOD_IDS
from .plan import DEFAULT_HARD_LIMIT_S, DEFAULT_REPETITIONS, DEFAULT_ROOT_SEED, DEFAULT_SOFT_LIMIT_S, build_jobs, fingerprint, make_spec
from .runner import REPO_ROOT, ExperimentError, Options, Runner, create_or_load, mem_total_bytes, physical_core_cpus, required_binaries
from .state import JobState
from .util import read_json

log = logging.getLogger("cbm_experiments")

# Options that define the experiment design. When resuming they may be omitted;
# if given, they must reproduce the stored design exactly.
DESIGN_OPTIONS = ("instance_names", "methods", "repetitions", "seed", "soft_limit", "hard_limit", "cbmlkh_threads", "method_params", "allow_param_overrides")


def parse_cpu_list(text: str) -> List[int]:
    cpus = []
    for part in text.split(","):
        if "-" in part:
            lo, hi = part.split("-")
            cpus.extend(range(int(lo), int(hi) + 1))
        elif part.strip():
            cpus.append(int(part))
    return sorted(set(cpus))


def solver_paths(args) -> Dict[str, str]:
    def pick(cli: Optional[str], env: str, default: Path) -> str:
        return str(Path(cli or os.environ.get(env) or default).expanduser().resolve())

    return {
        "ens": pick(args.ens_bin, "CBM_ENS_BIN", REPO_ROOT / "ENS" / "ENS"),
        "ils": pick(args.ils_bin, "CBM_ILS_BIN", REPO_ROOT / "ILS" / "ils"),
        "cbmlkh": pick(args.cbmlkh_bin, "CBM_CBMLKH_BIN", REPO_ROOT / "src" / "CBMLKH" / "main_prd"),
        "lkh_standalone": pick(args.lkh_standalone_bin, "CBM_LKH_STANDALONE_BIN", REPO_ROOT / "src" / "CBMLKH" / "lkh_standalone"),
        "lkh": pick(args.lkh_path, "LKH_PATH", REPO_ROOT / "src" / "LKH3" / "LKH"),
        "linkern": pick(args.linkern_path, "LINKERN_PATH", REPO_ROOT / "bin" / "linkern"),
    }


def add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--outputs", type=Path, default=REPO_ROOT / "outputs", help="root directory of all experiments (default: <repo>/outputs)")
    p.add_argument("--name", help="experiment name (default: 'final' for the full design, 'partial-<hash>' otherwise)")


def add_design(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("experiment design (fixed at creation; omit when resuming)")
    g.add_argument("--instance-dir", type=Path, default=None, help="instance directory (default: <repo>/instances)")
    g.add_argument("--instances", dest="instance_names", default=None, help="comma-separated subset of instance names (partial experiment)")
    g.add_argument("--methods", default=None, help=f"comma-separated subset of {','.join(METHOD_IDS)} (partial experiment)")
    g.add_argument("--repetitions", default=None, help=f"repetitions per method: N or method=N,... (default {DEFAULT_REPETITIONS})")
    g.add_argument("--seed", type=int, default=None, help=f"root random seed (default {DEFAULT_ROOT_SEED})")
    g.add_argument("--soft-limit", type=float, default=None, help=f"seconds; the methods' own time budget, never kills (default {DEFAULT_SOFT_LIMIT_S:g})")
    g.add_argument("--hard-limit", type=float, default=None, help=f"seconds; the process group is killed after this (default {DEFAULT_HARD_LIMIT_S:g})")
    g.add_argument("--cbmlkh-threads", type=int, default=None, help="threads (= repetitions) per CBMLKH execution (default min(10, CPU budget))")
    g.add_argument("--method-params", type=Path, default=None, help="JSON file {method: {param: value}} overriding ENS/ILS/LKH parameters")
    g.add_argument("--allow-param-overrides", action="store_true", default=None,
                   help="TESTING ONLY: let --method-params change the tuned CBMLKH and fixed LKH parameters (always a partial experiment)")


def add_resources(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("resources")
    g.add_argument("--cpus", default=None, help="CPU ids the runner may use, e.g. '0-15' (default: one per physical core)")
    g.add_argument("--max-cpus", type=int, default=None, help="use only the first N CPUs of the default set")
    g.add_argument("--no-pin", action="store_true", help="do not pin jobs to their reserved CPUs")
    g.add_argument("--mem-budget-gb", type=float, default=None, help="RAM the jobs may reserve (default 85%% of MemTotal)")
    g.add_argument("--disk-margin-gb", type=float, default=5.0, help="free disk to always keep on the work filesystem")
    g.add_argument("--work-root", type=Path, default=None, help="scratch directory for solver files (default <experiment>/work)")
    g.add_argument("--cbmlkh-lock", type=Path, default=None, help="machine-wide CBMLKH lock file (default $XDG_RUNTIME_DIR/cbm-experiments/cbmlkh.lock)")


def add_solvers(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("solver binaries (CLI > environment variable > default under the repository)")
    g.add_argument("--linkern-path", help="Linkern executable ($LINKERN_PATH; default <repo>/bin/linkern)")
    g.add_argument("--lkh-path", help="LKH executable ($LKH_PATH; default <repo>/src/LKH3/LKH)")
    g.add_argument("--ens-bin", help="ENS executable ($CBM_ENS_BIN; default <repo>/ENS/ENS)")
    g.add_argument("--ils-bin", help="ILS executable ($CBM_ILS_BIN; default <repo>/ILS/ils)")
    g.add_argument("--cbmlkh-bin", help="CBMLKH executable ($CBM_CBMLKH_BIN; default <repo>/src/CBMLKH/main_prd)")
    g.add_argument("--lkh-standalone-bin", help="standalone LKH driver ($CBM_LKH_STANDALONE_BIN; default <repo>/src/CBMLKH/lkh_standalone)")


def resolve_cpus(args) -> List[int]:
    cpus = parse_cpu_list(args.cpus) if args.cpus else physical_core_cpus()
    if args.max_cpus:
        cpus = cpus[: args.max_cpus]
    if not cpus:
        raise ExperimentError("no CPUs available")
    return cpus


def design_given(args) -> bool:
    return any(getattr(args, k, None) is not None for k in DESIGN_OPTIONS)


def build_spec(args, cpus: List[int]) -> dict:
    instance_dir = (args.instance_dir or REPO_ROOT / "instances").resolve()
    names = [n for n in args.instance_names.split(",") if n] if args.instance_names else None
    instances, skipped = discover_instances(instance_dir, names)
    for s in skipped:
        log.info("skipping %s", s)
    if not instances:
        raise ExperimentError(f"no instances found in {instance_dir}")
    methods = args.methods.split(",") if args.methods else list(METHOD_IDS)
    reps = {m: DEFAULT_REPETITIONS for m in methods}
    if args.repetitions:
        if "=" in args.repetitions:
            for item in args.repetitions.split(","):
                m, n = item.split("=")
                reps[m] = int(n)
        else:
            reps = {m: int(args.repetitions) for m in methods}
    overrides = json.loads(args.method_params.read_text()) if args.method_params else {}
    threads = args.cbmlkh_threads or min(reps.get("cbmlkh", DEFAULT_REPETITIONS), len(cpus))
    try:
        return _make_spec(args, methods, reps, overrides, threads, instances, instance_dir, names)
    except ValueError as exc:  # invalid design (e.g. overriding a locked parameter)
        raise ExperimentError(str(exc)) from exc


def _make_spec(args, methods, reps, overrides, threads, instances, instance_dir, names) -> dict:
    return make_spec(
        root_seed=DEFAULT_ROOT_SEED if args.seed is None else args.seed,
        methods=methods,
        repetitions=reps,
        soft_limit_s=args.soft_limit or DEFAULT_SOFT_LIMIT_S,
        hard_limit_s=args.hard_limit or DEFAULT_HARD_LIMIT_S,
        cbmlkh_threads=threads,
        param_overrides=overrides,
        instances=instances,
        instance_dir=str(instance_dir),
        all_instances_selected=names is None,
        allow_locked_overrides=bool(args.allow_param_overrides),
    )


def experiment_dir(args, spec: Optional[dict]) -> Path:
    if args.name:
        return args.outputs / args.name
    if spec is None:
        return args.outputs / "final"
    return args.outputs / ("partial-" + fingerprint(spec)[:10] if spec["partial"] else "final")


def setup_logging(exp_dir: Optional[Path]) -> None:
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    log.setLevel(logging.INFO)
    if not log.handlers:
        h = logging.StreamHandler(sys.stderr)
        h.setFormatter(fmt)
        log.addHandler(h)
    for h in [h for h in log.handlers if isinstance(h, logging.FileHandler)]:
        log.removeHandler(h)
        h.close()
    if exp_dir is not None:
        (exp_dir / "logs").mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(exp_dir / "logs" / "runner.log")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def cmd_run(args) -> int:
    cpus = resolve_cpus(args)
    exists = args.name and (args.outputs / args.name / "experiment.json").exists()
    spec = build_spec(args, cpus) if (design_given(args) or not exists) else None
    exp_dir = experiment_dir(args, spec)
    exp_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(exp_dir)
    experiment = create_or_load(exp_dir, exp_dir.name, spec)

    paths = solver_paths(args)
    missing = [f"{k} ({paths[k]})" for k in sorted(required_binaries(experiment["spec"])) if not os.access(paths[k], os.X_OK)]
    if missing:
        raise ExperimentError("solver binaries not found or not executable: " + ", ".join(missing) + " (build with `make`, or pass the paths)")

    options = Options(
        cpu_ids=cpus,
        pin_cpus=not args.no_pin,
        mem_budget_bytes=int((args.mem_budget_gb * 1024**3) if args.mem_budget_gb else 0.85 * mem_total_bytes()),
        disk_margin_bytes=int(args.disk_margin_gb * 1024**3),
        cbmlkh_lock_path=args.cbmlkh_lock,
        work_root=args.work_root,
        retry_failed=args.retry_failed,
        retry_timeouts=args.retry_timeouts,
        keep_work=args.keep_work,
        allow_concurrent_runners=args.allow_concurrent_runners,
        poll_interval_s=args.poll_interval,
        kill_grace_s=args.kill_grace,
        reference=report.load_reference(args.reference),
        max_jobs=args.max_jobs,
    )
    runner = Runner(exp_dir, experiment, paths, options, instance_dir=args.instance_dir if (exists and args.instance_dir) else None)
    log.info("experiment %s (%s, %d jobs) with CPUs %s, memory budget %.1f GiB", exp_dir, "partial" if experiment["partial"] else "full design",
             len(experiment["jobs"]), cpus, options.mem_budget_bytes / 1024**3)
    return runner.run()


def cmd_plan(args) -> int:
    cpus = resolve_cpus(args)
    spec = build_spec(args, cpus)
    jobs = build_jobs(spec)
    by_method: Dict[str, int] = {}
    for j in jobs:
        by_method[j["method"]] = by_method.get(j["method"], 0) + 1
    biggest = max(jobs, key=lambda j: j["est_mem_bytes"])
    print(json.dumps({
        "experiment_dir": str(experiment_dir(args, spec)),
        "partial": spec["partial"],
        "fingerprint": fingerprint(spec),
        "instances": len(spec["instances"]),
        "jobs": len(jobs),
        "jobs_by_method": by_method,
        "repetitions_by_method": spec["repetitions"],
        "cbmlkh_threads": spec["cbmlkh_threads"],
        "cpu_ids": cpus,
        "largest_memory_estimate_gib": round(biggest["est_mem_bytes"] / 1024**3, 2),
        "largest_memory_job": biggest["job_id"],
    }, indent=2))
    return 0


def _load_existing(args) -> Path:
    exp_dir = args.outputs / (args.name or "final")
    if read_json(exp_dir / "experiment.json") is None:
        raise ExperimentError(f"no experiment at {exp_dir}")
    return exp_dir


def cmd_status(args) -> int:
    exp_dir = _load_existing(args)
    experiment = read_json(exp_dir / "experiment.json")
    counts: Dict[str, Dict[str, int]] = {}
    for job in experiment["jobs"]:
        st = JobState(exp_dir / "runs" / job["job_id"], job["job_id"]).status
        counts.setdefault(job["method"], {}).setdefault(st, 0)
        counts[job["method"]][st] += len(job["repetitions"])
    print(json.dumps({"experiment": str(exp_dir), "partial": experiment["partial"], "repetitions_by_status": counts}, indent=2))
    return 0


def cmd_report(args) -> int:
    exp_dir = _load_existing(args)
    paths = report.generate(exp_dir, report.load_reference(args.reference))
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="run_experiments", description="Run and resume the CBM thesis experiments.")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="create or resume an experiment and execute its pending jobs")
    add_common(run)
    add_design(run)
    add_resources(run)
    add_solvers(run)
    g = run.add_argument_group("resume behaviour")
    g.add_argument("--retry-failed", action="store_true", help="re-run jobs whose last attempt failed")
    g.add_argument("--retry-timeouts", action="store_true", help="re-run jobs that hit the hard limit")
    g.add_argument("--allow-concurrent-runners", action="store_true", help="let several runners share one experiment (job locks still prevent double claims)")
    g.add_argument("--keep-work", action="store_true", help="keep solver scratch directories")
    g.add_argument("--max-jobs", type=int, default=None, help="launch at most N jobs in this session")
    g.add_argument("--poll-interval", type=float, default=1.0, help=argparse.SUPPRESS)
    g.add_argument("--kill-grace", type=float, default=30.0, help="seconds between SIGTERM and SIGKILL")
    g.add_argument("--reference", type=Path, default=None, help="CSV instance,value of best known solutions for gap columns")

    plan = sub.add_parser("plan", help="print the job plan without running anything")
    add_common(plan)
    add_design(plan)
    add_resources(plan)

    status = sub.add_parser("status", help="repetition counts per method and status")
    add_common(status)

    rep = sub.add_parser("report", help="regenerate the CSV reports from the persisted records")
    add_common(rep)
    rep.add_argument("--reference", type=Path, default=None, help="CSV instance,value of best known solutions")

    args = parser.parse_args(argv)
    if args.command != "run":
        setup_logging(None)
    try:
        return {"run": cmd_run, "plan": cmd_plan, "status": cmd_status, "report": cmd_report}[args.command](args)
    except ExperimentError as exc:
        log.error("%s", exc)
        return 2
