"""Persistent per-job state and advisory locks.

``runs/<job_id>/state.json`` is the source of truth for a job's status; it is
rewritten atomically on every transition, together with an append-only event
list. ``runs/<job_id>/lock`` is flock()ed by whoever executes the job; the lock
file descriptor is inherited by the solver process, so the lock is held exactly
as long as any process of the job is alive, even if the runner itself dies.
"""

from __future__ import annotations

import fcntl
import os
from pathlib import Path
from typing import Optional

from .util import atomic_write_json, read_json, utc_now

PENDING, RUNNING = "pending", "running"
COMPLETED, FAILED, HARD_TIMEOUT, INTERRUPTED = "completed", "failed", "hard_timeout", "interrupted"
JOB_STATUSES = (PENDING, RUNNING, COMPLETED, FAILED, HARD_TIMEOUT, INTERRUPTED)


class FileLock:
    """Non-blocking exclusive flock on a file. Released on close or process death."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.fd: Optional[int] = None

    def try_acquire(self) -> bool:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        fd = os.open(str(self.path), os.O_RDWR | os.O_CREAT, 0o644)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            os.close(fd)
            return False
        self.fd = fd
        return True

    def release(self) -> None:
        if self.fd is not None:
            try:
                fcntl.flock(self.fd, fcntl.LOCK_UN)
            finally:
                os.close(self.fd)
                self.fd = None

    def is_held_elsewhere(self) -> bool:
        if self.fd is not None:
            return False
        if self.try_acquire():
            self.release()
            return False
        return True


class JobState:
    def __init__(self, job_dir: Path, job_id: str):
        self.dir = Path(job_dir)
        self.path = self.dir / "state.json"
        self.lock = FileLock(self.dir / "lock")
        self.data = read_json(self.path) or {
            "job_id": job_id,
            "status": PENDING,
            "attempts": 0,
            "current_attempt": None,
            "result_files": [],
            "events": [],
        }

    @property
    def status(self) -> str:
        return self.data["status"]

    def reload(self) -> None:
        data = read_json(self.path)
        if data is not None:
            self.data = data

    def event(self, kind: str, **details) -> None:
        self.data["events"].append({"at": utc_now(), "event": kind, **details})

    def save(self) -> None:
        self.data["updated_at"] = utc_now()
        atomic_write_json(self.path, self.data)

    def transition(self, status: str, kind: str, **details) -> None:
        self.data["status"] = status
        self.event(kind, status=status, **details)
        self.save()
