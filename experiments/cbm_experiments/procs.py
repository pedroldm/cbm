"""Child-process management for solver runs (Linux).

Every job runs as the leader of a new session/process group, so one killpg()
reaches the solver and every LKH/Linkern process it spawned. The child is
pinned to its reserved CPUs (inherited by grandchildren) and gets SIGKILL if
the runner dies (PR_SET_PDEATHSIG), so a crashed runner leaves at most short-
lived grandchildren behind, which the next runner kills via the job's process
group before reusing the job (see runner.recover_stale_jobs).
"""

from __future__ import annotations

import ctypes
import os
import signal
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Sequence

_PR_SET_PDEATHSIG = 1

try:
    _libc = ctypes.CDLL("libc.so.6", use_errno=True)
except OSError:  # non-glibc platform: parent-death signal unavailable
    _libc = None


def boot_id() -> str:
    try:
        return Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    except OSError:
        return "unknown"


def _make_preexec(cpu_ids: Optional[Sequence[int]]):
    parent = os.getpid()

    def preexec() -> None:
        if _libc is not None:
            _libc.prctl(_PR_SET_PDEATHSIG, signal.SIGKILL, 0, 0, 0)
            if os.getppid() != parent:  # the runner died before prctl took effect
                os._exit(1)
        if cpu_ids:
            os.sched_setaffinity(0, set(cpu_ids))

    return preexec


def launch(
    argv: List[str],
    cwd: Path,
    env: Dict[str, str],
    stdout_path: Path,
    stderr_path: Path,
    cpu_ids: Optional[Sequence[int]],
    pass_fds: Sequence[int] = (),
) -> subprocess.Popen:
    with open(stdout_path, "ab") as out, open(stderr_path, "ab") as err:
        return subprocess.Popen(
            argv,
            cwd=str(cwd),
            env=env,
            stdin=subprocess.DEVNULL,
            stdout=out,
            stderr=err,
            start_new_session=True,
            pass_fds=tuple(pass_fds),
            preexec_fn=_make_preexec(cpu_ids),
            close_fds=True,
        )


def signal_group(pgid: int, sig: int) -> bool:
    """Send `sig` to a process group; False if the group no longer exists."""
    try:
        os.killpg(pgid, sig)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return False


def group_members(pgid: int) -> List[int]:
    """PIDs whose process group is `pgid` (scans /proc)."""
    members = []
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            stat = Path(f"/proc/{entry}/stat").read_text()
        except OSError:
            continue
        # Field 5 is pgrp; the command name (field 2) may contain spaces, so
        # split after its closing parenthesis.
        fields = stat[stat.rfind(")") + 2 :].split()
        if len(fields) > 2 and int(fields[2]) == pgid:
            members.append(int(entry))
    return members
