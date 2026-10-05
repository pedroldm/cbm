"""Fault-tolerant experiment runner.

Layout of an experiment directory (``<outputs>/<name>/``)::

    experiment.json          immutable design: spec, fingerprint, job list
    manifest.json            derived snapshot of every job's status (refreshed)
    runner.lock              held by the active runner process
    runs/<job_id>/state.json authoritative job status + event log (atomic)
    runs/<job_id>/lock       held by the processes executing the job
    runs/<job_id>/attempt-NN/{command.json,stdout.log,stderr.log,native.json,...}
    runs/<job_id>/<run_id>.json           shared-schema record per repetition
    runs/<job_id>/<run_id>.solution.json  best permutation of that repetition
    results/*.csv            reports (see report.py)
    logs/runner.log, logs/sessions.jsonl
    work/                    per-attempt scratch dirs (deleted after each job)

Statuses: pending -> running -> completed | failed | hard_timeout | interrupted.
A job is "completed" only after its records are on disk and pass validation;
on restart, records are re-checked and a completed job with a missing or
invalid record is reset to pending.
"""

from __future__ import annotations

import json
import logging
import os
import platform
import shutil
import signal
import socket
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set

from . import methods as M
from . import procs, report
from .instances import InstanceCache
from .plan import build_jobs, fingerprint, run_id
from .schema import SCHEMA_VERSION, load_schema, validate
from .state import COMPLETED, FAILED, HARD_TIMEOUT, INTERRUPTED, PENDING, RUNNING, FileLock, JobState
from .util import atomic_write_json, atomic_write_text, read_json, sha256_file, tail, utc_now

RUNNER_VERSION = "1.0.0"
REPO_ROOT = Path(__file__).resolve().parents[2]
log = logging.getLogger("cbm_experiments")


class ExperimentError(RuntimeError):
    pass


@dataclass
class Options:
    poll_interval_s: float = 1.0
    kill_grace_s: float = 30.0
    retry_failed: bool = False
    retry_timeouts: bool = False
    max_attempts: int = 5
    cpu_ids: List[int] = field(default_factory=list)
    pin_cpus: bool = True
    mem_budget_bytes: int = 0
    disk_margin_bytes: int = 5 * 1024**3
    cbmlkh_lock_path: Optional[Path] = None
    allow_concurrent_runners: bool = False
    keep_work: bool = False
    work_root: Optional[Path] = None
    report_interval_s: float = 10.0
    manifest_interval_s: float = 5.0
    reference: Dict[str, float] = field(default_factory=dict)
    max_jobs: Optional[int] = None  # stop launching after this many (testing / partial sessions)


@dataclass
class Running:
    job: dict
    state: JobState
    proc: subprocess.Popen
    attempt: int
    attempt_dir: Path
    work_dir: Path
    native_out: Path
    command: List[str]
    cpu_ids: List[int]
    started_mono: float
    started_at: str
    cbmlkh_lock: Optional[FileLock] = None
    soft_reached_at: Optional[float] = None
    kill_reason: Optional[str] = None  # "hard_timeout" or "shutdown"
    kill_sent_mono: Optional[float] = None
    kill_escalated: bool = False


# ---------------------------------------------------------------------------
# Environment helpers
# ---------------------------------------------------------------------------
def physical_core_cpus() -> List[int]:
    """One logical CPU per physical core, within this process's affinity mask."""
    allowed = sorted(os.sched_getaffinity(0))
    chosen, seen = [], set()
    for cpu in allowed:
        try:
            core = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list").read_text().strip()
        except OSError:
            return allowed
        if core not in seen:
            seen.add(core)
            chosen.append(cpu)
    return chosen


def mem_total_bytes() -> int:
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemTotal:"):
                return int(line.split()[1]) * 1024
    except OSError:
        pass
    return 0


def git_info() -> dict:
    def git(*args) -> Optional[str]:
        try:
            return subprocess.run(["git", "-C", str(REPO_ROOT), *args], capture_output=True, text=True, timeout=30, check=True).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return None

    status = git("status", "--porcelain", "--untracked-files=no")
    return {"commit": git("rev-parse", "HEAD"), "dirty": None if status is None else bool(status)}


def binary_info(paths: Dict[str, str]) -> Dict[str, dict]:
    info = {}
    for name, path in sorted(paths.items()):
        p = Path(path)
        info[name] = {"path": str(p), "sha256": sha256_file(p) if p.is_file() else None}
    return info


def required_binaries(spec: dict) -> Set[str]:
    needed = set()
    for m in spec["methods"]:
        backend = spec["method_params"][m]["tsp_backend"]
        if m in ("ens", "ils"):
            needed |= {m, backend}
        elif m == "lkh":
            needed |= {"lkh_standalone", "lkh"}
        elif m == "cbmlkh":
            needed |= {"cbmlkh", "lkh"}
    return needed


# ---------------------------------------------------------------------------
# Experiment directory
# ---------------------------------------------------------------------------
def create_or_load(exp_dir: Path, name: str, spec: Optional[dict]) -> dict:
    """Return experiment.json content, creating it from `spec` if absent.

    An existing experiment is only resumed with an identical design: a
    different spec under the same name is refused instead of mixing results.
    """
    existing = read_json(exp_dir / "experiment.json")
    if existing is not None:
        if spec is not None and fingerprint(spec) != existing["fingerprint"]:
            raise ExperimentError(
                f"{exp_dir} holds a different experiment design (fingerprint {existing['fingerprint'][:12]} vs "
                f"{fingerprint(spec)[:12]}). Use another --name, or resume without changing design options."
            )
        return existing
    if spec is None:
        raise ExperimentError(f"{exp_dir} does not contain an experiment to resume")
    if spec["partial"] and name == "final":
        raise ExperimentError("this design deviates from the full experiment; it cannot be named 'final'")
    experiment = {
        "schema_version": SCHEMA_VERSION,
        "name": name,
        "fingerprint": fingerprint(spec),
        "partial": spec["partial"],
        "created_at": utc_now(),
        "runner_version": RUNNER_VERSION,
        "spec": spec,
        "jobs": build_jobs(spec),
    }
    atomic_write_json(exp_dir / "experiment.json", experiment)
    return experiment


class Runner:
    def __init__(self, exp_dir: Path, experiment: dict, paths: Dict[str, str], options: Options, instance_dir: Optional[Path] = None):
        self.exp_dir = Path(exp_dir)
        self.experiment = experiment
        self.spec = experiment["spec"]
        self.paths = paths
        self.opt = options
        self.jobs: List[dict] = experiment["jobs"]
        if instance_dir is not None:  # experiment moved to another machine
            for job in self.jobs:
                job["instance"] = dict(job["instance"], path=str(Path(instance_dir).resolve() / job["instance"]["name"]))
        self.states: Dict[str, JobState] = {j["job_id"]: JobState(self.exp_dir / "runs" / j["job_id"], j["job_id"]) for j in self.jobs}
        self.running: Dict[str, Running] = {}
        self.stop_requested = 0
        self.launched = 0
        self.schema = load_schema()
        self.instances = InstanceCache()
        self.record_cache = report.RecordCache()
        self.work_root = Path(options.work_root or self.exp_dir / "work")
        self.cbmlkh_lock_path = Path(options.cbmlkh_lock_path or Path(os.environ.get("XDG_RUNTIME_DIR", "/tmp")) / "cbm-experiments" / "cbmlkh.lock")
        self.cpu_ids = list(options.cpu_ids)
        self.host = socket.gethostname()
        self.boot = procs.boot_id()
        self._dirty_reports = True
        self._dirty_manifest = True
        self._last_report = 0.0
        self._last_manifest = 0.0
        self._verified_instances: Set[str] = set()
        self.runner_lock = FileLock(self.exp_dir / "runner.lock")

    # -- lifecycle -----------------------------------------------------------
    def run(self) -> int:
        if not self.runner_lock.try_acquire() and not self.opt.allow_concurrent_runners:
            raise ExperimentError(f"another runner is active on {self.exp_dir} (runner.lock is held)")
        self._check_resources_fit()
        self.software = {
            "runner_version": RUNNER_VERSION,
            "git": git_info(),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "binaries": binary_info({k: v for k, v in self.paths.items() if k in required_binaries(self.spec)}),
        }
        self._session_event("start")
        previous = {s: signal.getsignal(s) for s in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)}
        for s in previous:
            signal.signal(s, self._on_signal)
        try:
            self.verify_completed_jobs()
            self.recover_stale_jobs()
            self._loop()
        finally:
            for s, handler in previous.items():
                signal.signal(s, handler)
            self.write_manifest()
            self.write_reports()
            self._session_event("stop", stop_requested=bool(self.stop_requested))
            self.runner_lock.release()
        counts = self.counts()
        log.info("session finished: %s", counts)
        return 0 if not self.stop_requested else 130

    def _on_signal(self, signum, _frame) -> None:
        self.stop_requested += 1
        log.warning("received signal %d: stopping (running jobs are terminated and marked interrupted)", signum)
        for rj in list(self.running.values()):
            if rj.kill_reason is None:
                rj.kill_reason = "shutdown"
                rj.kill_sent_mono = time.monotonic()
                procs.signal_group(rj.proc.pid, signal.SIGTERM)
            elif self.stop_requested > 1:
                procs.signal_group(rj.proc.pid, signal.SIGKILL)

    def _session_event(self, kind: str, **extra) -> None:
        line = {"at": utc_now(), "event": kind, "host": self.host, "pid": os.getpid(), "cpu_ids": self.cpu_ids,
                "mem_budget_bytes": self.opt.mem_budget_bytes, "paths": self.paths, **extra}
        (self.exp_dir / "logs").mkdir(parents=True, exist_ok=True)
        with open(self.exp_dir / "logs" / "sessions.jsonl", "a") as f:
            f.write(json.dumps(line) + "\n")
            f.flush()
            os.fsync(f.fileno())

    def _check_resources_fit(self) -> None:
        widest = max(j["cpus"] for j in self.jobs)
        if widest > len(self.cpu_ids):
            raise ExperimentError(
                f"the plan contains {widest}-thread CBMLKH jobs but only {len(self.cpu_ids)} CPUs are available to the runner; "
                "running it would oversubscribe the machine (create the experiment with a matching --cbmlkh-threads)"
            )

    # -- startup checks ------------------------------------------------------
    def _records_valid(self, job: dict) -> Optional[str]:
        job_dir = self.exp_dir / "runs" / job["job_id"]
        for rep in job["repetitions"]:
            rid = run_id(job["method"], job["instance"]["name"], rep)
            rec = read_json(job_dir / f"{rid}.json")
            if rec is None:
                return f"record {rid}.json missing or unreadable"
            errors = validate(rec, self.schema)
            if errors:
                return f"record {rid}.json violates the schema: {errors[:3]}"
            if rec["status"] != COMPLETED or rec["objective"]["value"] is None or rec["objective"]["validated"] is not True:
                return f"record {rid}.json is not a validated completed run"
            if read_json(job_dir / f"{rid}.solution.json") is None:
                return f"solution {rid}.solution.json missing or unreadable"
        return None

    def verify_completed_jobs(self) -> None:
        for job in self.jobs:
            st = self.states[job["job_id"]]
            if st.status == COMPLETED:
                problem = self._records_valid(job)
                if problem:
                    log.warning("%s: completed but %s; resetting to pending", job["job_id"], problem)
                    st.transition(PENDING, "result_invalid", reason=problem)
                    self._dirty_manifest = True

    def recover_stale_jobs(self) -> None:
        """Resolve 'running' states left by a runner that crashed or lost power."""
        for job in self.jobs:
            st = self.states[job["job_id"]]
            if st.status != RUNNING:
                continue
            attempt = st.data.get("current_attempt") or {}
            if st.lock.is_held_elsewhere():
                if self.opt.allow_concurrent_runners:
                    continue  # presumably executed by another runner
                pgid = attempt.get("pgid")
                members = procs.group_members(pgid) if pgid and attempt.get("boot_id") == self.boot else []
                log.warning("%s: orphaned processes %s still hold the job lock; killing them", job["job_id"], members)
                if pgid:
                    procs.signal_group(pgid, signal.SIGKILL)
                deadline = time.monotonic() + 30
                while st.lock.is_held_elsewhere() and time.monotonic() < deadline:
                    time.sleep(0.2)
                if st.lock.is_held_elsewhere():
                    log.error("%s: lock still held after killing process group %s; leaving it alone", job["job_id"], pgid)
                    continue
                st.transition(INTERRUPTED, "orphan_killed", attempt=attempt.get("number"), pgid=pgid, pids=members)
            else:
                same_boot = attempt.get("boot_id") == self.boot
                st.transition(INTERRUPTED, "stale_running_recovered", attempt=attempt.get("number"),
                              reason="runner stopped without finalizing (crash, kill or power loss)" if same_boot else "machine rebooted")
            self._dirty_manifest = True

    # -- scheduling ----------------------------------------------------------
    def _eligible(self, st: JobState) -> bool:
        if st.data["attempts"] >= self.opt.max_attempts and st.status != PENDING:
            return False
        if st.status in (PENDING, INTERRUPTED):
            return True
        return (st.status == FAILED and self.opt.retry_failed) or (st.status == HARD_TIMEOUT and self.opt.retry_timeouts)

    def _free_cpus(self) -> List[int]:
        used = {c for rj in self.running.values() for c in rj.cpu_ids}
        return [c for c in self.cpu_ids if c not in used]

    def _loop(self) -> None:
        while True:
            for job_id in list(self.running):
                self._monitor(self.running[job_id])
            if self.stop_requested:
                if not self.running:
                    break
            else:
                pending = [j for j in self.jobs if j["job_id"] not in self.running and self._eligible(self.states[j["job_id"]])]
                if not pending and not self.running:
                    break
                limit_reached = self.opt.max_jobs is not None and self.launched >= self.opt.max_jobs
                if not limit_reached:
                    self._schedule(pending)
                elif not self.running:
                    break
            self._maybe_flush()
            time.sleep(self.opt.poll_interval_s)

    def _schedule(self, pending: List[dict]) -> None:
        cbmlkh_pending = [j for j in pending if j["method"] == "cbmlkh"]
        cbmlkh_running = any(rj.job["method"] == "cbmlkh" for rj in self.running.values())
        # CBMLKH gets a dedicated lane while any of it remains, so single-thread
        # jobs cannot starve it of the CPUs it needs at once.
        next_cbmlkh = cbmlkh_pending[0] if (cbmlkh_pending and not cbmlkh_running) else None
        lane_cpus = next_cbmlkh["cpus"] if next_cbmlkh else 0
        lane_mem = next_cbmlkh["est_mem_bytes"] if next_cbmlkh and next_cbmlkh["est_mem_bytes"] <= self.opt.mem_budget_bytes else 0
        for job in cbmlkh_pending[:1] + [j for j in pending if j["method"] != "cbmlkh"]:
            free = self._free_cpus()
            is_cbmlkh = job["method"] == "cbmlkh"
            if is_cbmlkh:
                if cbmlkh_running or len(free) < job["cpus"]:
                    continue
            elif len(free) - lane_cpus < job["cpus"]:
                break
            verdict = self._memory_and_disk_allow(job, 0 if is_cbmlkh else lane_mem)
            if verdict == "drain":
                break  # an oversized job waits for the machine to empty; start nothing else
            if verdict != "ok":
                continue
            if self._launch(job, free[: job["cpus"]]) and is_cbmlkh:
                cbmlkh_running, lane_cpus, lane_mem = True, 0, 0
            if self.opt.max_jobs is not None and self.launched >= self.opt.max_jobs:
                return

    def _memory_and_disk_allow(self, job: dict, reserved_mem: int) -> str:
        """'ok', 'wait' (try other jobs) or 'drain' (job needs the whole machine)."""
        used_mem = sum(rj.job["est_mem_bytes"] for rj in self.running.values())
        used_disk = sum(rj.job["est_disk_bytes"] for rj in self.running.values())
        self.work_root.mkdir(parents=True, exist_ok=True)
        free_disk = shutil.disk_usage(self.work_root).free - self.opt.disk_margin_bytes
        if used_mem + reserved_mem + job["est_mem_bytes"] <= self.opt.mem_budget_bytes and used_disk + job["est_disk_bytes"] <= free_disk:
            return "ok"
        oversized = job["est_mem_bytes"] > self.opt.mem_budget_bytes or job["est_disk_bytes"] > free_disk + used_disk
        if not oversized:
            return "wait"
        if self.running:
            return "drain"
        # Larger than the whole budget: run it alone rather than never.
        log.warning("%s: estimated %.1f GiB RAM / %.1f GiB disk exceeds the budget; running it alone",
                    job["job_id"], job["est_mem_bytes"] / M.GIB, job["est_disk_bytes"] / M.GIB)
        return "ok"

    # -- launch --------------------------------------------------------------
    def _launch(self, job: dict, cpu_ids: List[int]) -> bool:
        st = self.states[job["job_id"]]
        if not st.lock.try_acquire():
            return False  # another process holds this job
        st.reload()
        if not self._eligible(st):
            st.lock.release()
            return False
        cbmlkh_lock = None
        if job["method"] == "cbmlkh":
            cbmlkh_lock = FileLock(self.cbmlkh_lock_path)
            if not cbmlkh_lock.try_acquire():
                st.lock.release()
                return False

        try:
            return self._start(job, st, cpu_ids, cbmlkh_lock)
        except Exception as exc:  # setup failure: record it, never leave the job 'running'
            log.exception("%s: could not start", job["job_id"])
            st.data["attempts"] += 1
            st.transition(FAILED, "launch_failed", error=f"{type(exc).__name__}: {exc}")
            st.lock.release()
            if cbmlkh_lock:
                cbmlkh_lock.release()
            self._dirty_manifest = self._dirty_reports = True
            return False

    def _start(self, job: dict, st: JobState, cpu_ids: List[int], cbmlkh_lock: Optional[FileLock]) -> bool:
        inst = job["instance"]
        if inst["path"] not in self._verified_instances:
            if sha256_file(Path(inst["path"])) != inst["sha256"]:
                raise ExperimentError(f"instance {inst['path']} changed since the experiment was planned (SHA-256 mismatch)")
            self._verified_instances.add(inst["path"])

        attempt = st.data["attempts"] + 1
        job_dir = self.exp_dir / "runs" / job["job_id"]
        attempt_dir = job_dir / f"attempt-{attempt:02d}"
        attempt_dir.mkdir(parents=True, exist_ok=True)
        self._supersede_records(job, job_dir, attempt - 1)

        work_dir = self.work_root / job["job_id"] / f"attempt-{attempt:02d}"
        shutil.rmtree(work_dir, ignore_errors=True)
        work_dir.mkdir(parents=True)
        native_out = attempt_dir / "native.json"
        cfg_file = None
        if job["method"] == "cbmlkh":
            cfg_file = attempt_dir / "cbmlkh.cfg"
            atomic_write_text(cfg_file, M.cbmlkh_config_text(job, self.spec, self.paths, str(native_out), str(work_dir)))
        argv = M.build_command(job, self.spec, self.paths, str(native_out), str(work_dir), str(cfg_file) if cfg_file else None)
        env = dict(os.environ)
        env.update({"OMP_NUM_THREADS": str(job["threads"]), "OMP_DYNAMIC": "false", "TMPDIR": str(work_dir), "LC_ALL": "C"})
        pinned = cpu_ids if self.opt.pin_cpus else None
        atomic_write_json(attempt_dir / "command.json", {"argv": argv, "cpu_ids": pinned, "env": {k: env[k] for k in ("OMP_NUM_THREADS", "OMP_DYNAMIC", "TMPDIR", "LC_ALL")}})

        fds = [st.lock.fd] + ([cbmlkh_lock.fd] if cbmlkh_lock else [])
        proc = procs.launch(argv, work_dir, env, attempt_dir / "stdout.log", attempt_dir / "stderr.log", pinned, fds)
        started_at = utc_now()
        st.data["attempts"] = attempt
        st.data["current_attempt"] = {"number": attempt, "pid": proc.pid, "pgid": proc.pid, "host": self.host, "boot_id": self.boot,
                                      "started_at": started_at, "cpu_ids": cpu_ids, "dir": str(attempt_dir.relative_to(self.exp_dir))}
        st.transition(RUNNING, "started", attempt=attempt, pid=proc.pid, cpu_ids=cpu_ids)
        self.running[job["job_id"]] = Running(job, st, proc, attempt, attempt_dir, work_dir, native_out, argv, cpu_ids,
                                              time.monotonic(), started_at, cbmlkh_lock)
        self.launched += 1
        self._dirty_manifest = True
        log.info("started %s (attempt %d, pid %d, cpus %s)", job["job_id"], attempt, proc.pid, cpu_ids)
        return True

    def _supersede_records(self, job: dict, job_dir: Path, previous_attempt: int) -> None:
        """Move records of an earlier, non-successful attempt out of the way."""
        target = job_dir / f"attempt-{max(previous_attempt, 0):02d}" / "superseded"
        for rep in job["repetitions"]:
            rid = run_id(job["method"], job["instance"]["name"], rep)
            for name in (f"{rid}.json", f"{rid}.solution.json"):
                src = job_dir / name
                if src.exists():
                    target.mkdir(parents=True, exist_ok=True)
                    os.replace(src, target / name)

    # -- monitoring ----------------------------------------------------------
    def _monitor(self, rj: Running) -> None:
        now = time.monotonic()
        elapsed = now - rj.started_mono
        if rj.soft_reached_at is None and elapsed >= self.spec["soft_limit_s"]:
            # The soft limit is the method's own time budget: record it, keep running.
            rj.soft_reached_at = elapsed
            rj.state.event("soft_limit_reached", elapsed_s=round(elapsed, 3))
            rj.state.save()
            log.info("%s reached the soft limit (%.0f s); still running", rj.job["job_id"], elapsed)
        rc = rj.proc.poll()
        if rc is None:
            if rj.kill_reason is None and elapsed >= self.spec["hard_limit_s"]:
                rj.kill_reason = "hard_timeout"
                rj.kill_sent_mono = now
                rj.state.event("hard_limit_reached", elapsed_s=round(elapsed, 3))
                rj.state.save()
                log.warning("%s exceeded the hard limit (%.0f s); terminating its process group", rj.job["job_id"], elapsed)
                procs.signal_group(rj.proc.pid, signal.SIGTERM)
            elif rj.kill_reason and not rj.kill_escalated and now - rj.kill_sent_mono >= self.opt.kill_grace_s:
                rj.kill_escalated = True
                procs.signal_group(rj.proc.pid, signal.SIGKILL)
            return
        # Leader gone: make sure no descendant (LKH/Linkern) outlives it.
        procs.signal_group(rj.proc.pid, signal.SIGKILL)
        self._finalize(rj, rc, elapsed)

    # -- finalization ----------------------------------------------------------
    def _finalize(self, rj: Running, rc: int, elapsed: float) -> None:
        job = rj.job
        finished_at = utc_now()
        parsed: Optional[Dict[int, dict]] = None
        error = None
        status = FAILED
        stderr_tail = tail(rj.attempt_dir / "stderr.log")

        if rc == 0:
            try:
                parsed = self._parse_and_validate(job, rj.native_out)
                status = COMPLETED
            except Exception as exc:
                error = {"type": "invalid_output", "message": f"{type(exc).__name__}: {exc}", "stderr_tail": stderr_tail}
            if rj.kill_reason:
                rj.state.event("finished_despite_termination", reason=rj.kill_reason)
        elif rj.kill_reason == "hard_timeout":
            status = HARD_TIMEOUT
            error = {"type": "hard_timeout", "message": f"terminated after exceeding the hard limit of {self.spec['hard_limit_s']:g} s", "stderr_tail": stderr_tail}
        elif rj.kill_reason == "shutdown":
            status = INTERRUPTED
            error = {"type": "interrupted", "message": "terminated because the runner was stopped", "stderr_tail": stderr_tail}
        else:
            what = f"killed by signal {-rc}" if rc < 0 else f"exit status {rc}"
            error = {"type": "process_failed", "message": f"solver {what}", "stderr_tail": stderr_tail}

        job_dir = self.exp_dir / "runs" / job["job_id"]
        result_files = []
        try:
            for rep in job["repetitions"]:
                rid = run_id(job["method"], job["instance"]["name"], rep)
                part = parsed.get(rep) if parsed else None
                perm_file = None
                if part is not None:
                    perm_file = f"{rid}.solution.json"
                    atomic_write_json(job_dir / perm_file, {"run_id": rid, "indexing": "1-based column ids, left to right", "permutation": part["permutation"]})
                record = self._record(rj, rep, status, part, error, rc, elapsed, finished_at, perm_file)
                problems = validate(record, self.schema)
                if problems:
                    raise ExperimentError(f"internal error: record for {rid} violates the schema: {problems}")
                atomic_write_json(job_dir / f"{rid}.json", record)
                result_files.append(f"{rid}.json")
        except Exception as exc:
            log.exception("%s: could not persist results", job["job_id"])
            # Never leave a partial set of records behind a non-completed job.
            for rep in job["repetitions"]:
                rid = run_id(job["method"], job["instance"]["name"], rep)
                for name in (f"{rid}.json", f"{rid}.solution.json"):
                    try:
                        (job_dir / name).unlink()
                    except OSError:
                        pass
            result_files = []
            status = FAILED
            error = {"type": "persist_failed", "message": f"{type(exc).__name__}: {exc}", "stderr_tail": stderr_tail}
            rj.state.event("persist_failed", error=error["message"])

        rj.state.data["result_files"] = result_files
        rj.state.data["current_attempt"] = None
        rj.state.transition(status, "finished", attempt=rj.attempt, exit_code=rc, wall_time_s=round(elapsed, 3),
                            error=(error or {}).get("message"))
        self._release(rj)
        log.info("finished %s: %s (%.1f s)", job["job_id"], status, elapsed)

    def _release(self, rj: Running) -> None:
        self.running.pop(rj.job["job_id"], None)
        rj.state.lock.release()
        if rj.cbmlkh_lock:
            rj.cbmlkh_lock.release()
        # Keep solver logs, drop bulky scratch files.
        if rj.work_dir.exists():
            for p in rj.work_dir.rglob("*.log"):
                if p.is_file() and p.stat().st_size <= 1024**2:
                    dest = rj.attempt_dir / "work_logs" / p.relative_to(rj.work_dir)
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(p, dest)
            if not self.opt.keep_work:
                shutil.rmtree(rj.work_dir, ignore_errors=True)
                try:
                    rj.work_dir.parent.rmdir()  # the job's directory, once no attempt is left in it
                except OSError:
                    pass
        self._dirty_manifest = self._dirty_reports = True

    def _parse_and_validate(self, job: dict, native_out: Path) -> Dict[int, dict]:
        native = read_json(native_out)
        if native is None:
            raise ValueError(f"solver exited 0 but {native_out.name} is missing or not valid JSON")
        parsed = M.parse_native(job, native)
        inst = self.instances.get(job["instance"]["path"])
        for rep, part in parsed.items():
            defect = inst.check_permutation(part["permutation"])
            if defect:
                raise ValueError(f"repetition {rep}: {defect}")
            recount = inst.count_blocks(part["permutation"])
            if recount != part["value"]:
                raise ValueError(f"repetition {rep}: solver reported {part['value']} blocks, independent recount gives {recount}")
        return parsed

    def _record(self, rj: Running, rep: int, status: str, part: Optional[dict], error: Optional[dict], rc: int, elapsed: float,
                finished_at: str, perm_file: Optional[str]) -> dict:
        job, spec = rj.job, self.spec
        params = spec["method_params"][job["method"]]
        inst = job["instance"]
        tsp = (part or {}).get("tsp") or {}
        criteria = {"soft_limit_s": spec["soft_limit_s"], "hard_limit_s": spec["hard_limit_s"]}
        if job["method"] in ("ens", "ils"):
            criteria["iterations"] = params["iterations"]
        elif job["method"] == "cbmlkh":
            criteria.update(maxIterations=params["maxIterations"], maxTime=int(-(-spec["soft_limit_s"] // 1)), lkhMaxTime=params["lkhMaxTime"])
        elif job["method"] == "lkh":
            criteria["runs"] = params["runs"]
        return {
            "schema_version": SCHEMA_VERSION,
            "run_id": run_id(job["method"], inst["name"], rep),
            "job_id": job["job_id"],
            "experiment": {"name": self.experiment["name"], "fingerprint": self.experiment["fingerprint"], "partial": self.experiment["partial"]},
            "method": {"id": job["method"], "name": M.METHOD_NAMES[job["method"]], "tsp_backend": params["tsp_backend"], "parameters": params},
            "instance": {"name": inst["name"], "path": inst["path"], "sha256": inst["sha256"], "rows": inst["rows"], "cols": inst["cols"]},
            "repetition": rep,
            "seeds": {
                "root": spec["root_seed"],
                "run": job["seed"],
                "derivation": "run = seeds.run_seed(root, method, instance, repetition) (repetition = -1 for CBMLKH, whose trajectories derive from it)",
                "derived": (part or {}).get("derived_seeds") or {},
            },
            "status": status,
            "objective": {
                "name": "one_blocks",
                "sense": "minimize",
                "value": part["value"] if part else None,
                "blocks": part["value"] if part else None,
                "validated": True if part else None,
                "initial_value": part["initial_value"] if part else None,
            },
            "solution": {"permutation_file": perm_file},
            "timing": {
                "started_at": rj.started_at,
                "finished_at": finished_at,
                "wall_time_s": round(elapsed, 6),
                "method_time_s": (part or {}).get("method_time_s"),
                "time_to_best_s": (part or {}).get("time_to_best_s"),
                "soft_limit_s": spec["soft_limit_s"],
                "hard_limit_s": spec["hard_limit_s"],
                "soft_limit_reached": rj.soft_reached_at is not None,
                "soft_limit_reached_at_s": round(rj.soft_reached_at, 3) if rj.soft_reached_at is not None else None,
            },
            "stopping": {"criteria": criteria, "reason": (part or {}).get("stop_reason") or (status if status != COMPLETED else None),
                         "method_time_limit_reached": (part or {}).get("method_time_limit_reached")},
            "search": {"iterations": (part or {}).get("iterations"), "improvements": (part or {}).get("improvements")},
            "tsp_solver": {"name": params["tsp_backend"], "calls": tsp.get("calls"), "time_s": tsp.get("time_s"),
                           "time_limit_hits": tsp.get("time_limit_hits"), "invalid_tours": tsp.get("invalid_tours")},
            "execution": {
                "attempt": rj.attempt,
                "exit_code": rc if rc >= 0 else None,
                "signal": -rc if rc < 0 else None,
                "cpus": job["cpus"],
                "cpu_ids": rj.cpu_ids if self.opt.pin_cpus else None,
                "threads": job["threads"],
                "host": self.host,
                "command": rj.command,
                "stdout_log": str((rj.attempt_dir / "stdout.log").relative_to(self.exp_dir)),
                "stderr_log": str((rj.attempt_dir / "stderr.log").relative_to(self.exp_dir)),
                "native_output": str(rj.native_out.relative_to(self.exp_dir)) if rj.native_out.exists() else None,
            },
            "error": error,
            "software": self.software,
            "method_stats": (part or {}).get("method_stats"),
        }

    # -- derived files -------------------------------------------------------
    def counts(self) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for st in self.states.values():
            out[st.status] = out.get(st.status, 0) + 1
        return out

    def write_manifest(self) -> None:
        jobs = []
        for job in self.jobs:
            st = self.states[job["job_id"]]
            jobs.append({
                "job_id": job["job_id"], "method": job["method"], "instance": job["instance"]["name"], "repetitions": job["repetitions"],
                "seed": job["seed"], "cpus": job["cpus"], "status": st.status, "attempts": st.data["attempts"],
                "state_file": f"runs/{job['job_id']}/state.json",
                "records": [f"runs/{job['job_id']}/{run_id(job['method'], job['instance']['name'], r)}.json" for r in job["repetitions"]],
                "updated_at": st.data.get("updated_at"),
            })
        atomic_write_json(self.exp_dir / "manifest.json", {
            "schema_version": SCHEMA_VERSION, "experiment": self.experiment["name"], "fingerprint": self.experiment["fingerprint"],
            "partial": self.experiment["partial"], "generated_at": utc_now(), "note": "derived snapshot; runs/*/state.json is authoritative",
            "counts": self.counts(), "jobs": jobs,
        })
        self._dirty_manifest = False
        self._last_manifest = time.monotonic()

    def write_reports(self) -> None:
        try:
            report.generate(self.exp_dir, self.opt.reference, self.record_cache)
        except Exception:
            log.exception("report generation failed")
        self._dirty_reports = False
        self._last_report = time.monotonic()

    def _maybe_flush(self) -> None:
        now = time.monotonic()
        if self._dirty_manifest and now - self._last_manifest >= self.opt.manifest_interval_s:
            self.write_manifest()
        if self._dirty_reports and now - self._last_report >= self.opt.report_interval_s:
            self.write_reports()
