"""Benchmark orchestrator: run a solver over a chosen subset of instances.

Solver-agnostic — it loads each instance, delegates solving to a :class:`Solver`
(pure LKH or the CBMLKH binary), and persists comparable, resumable results so a
long ``nohup`` run can be restarted without recomputing or overwriting finished
instances.
"""

from __future__ import annotations

import datetime as _dt
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .config import LOG_PREFIX
from .instance import CBMInstance
from .result_store import ResultStore
from .solvers.base import Solver


def _now_iso() -> str:
    return _dt.datetime.now().astimezone().isoformat(timespec="seconds")


class Benchmark:
    """Run one :class:`Solver` on CBM instances and store resumable results."""

    def __init__(
        self,
        solver: Solver,
        instances_dir: str | Path,
        results_dir: str | Path,
        run_tag: str | None = None,
    ) -> None:
        self.solver = solver
        self.instances_dir = Path(instances_dir)
        self.run_tag = run_tag or solver.run_tag()
        self.store = ResultStore(results_dir, self.run_tag)

    # ------------------------------------------------------------------ public
    def run(self, instance_names, skip_done: bool = True, workers: int = 1) -> list[dict]:
        """Run the solver over ``instance_names`` and return all stored records.

        With ``workers > 1`` the instances are processed concurrently on a thread
        pool. This suits solvers whose per-instance work is a single-threaded
        external process (pure LKH); the CBMLKH solver is already internally
        parallel, so leave it at ``workers = 1`` to avoid oversubscribing cores.
        """
        names = list(instance_names)
        total = len(names)
        workers = max(1, workers)
        self._log(f"solver='{self.solver.key}' run_tag='{self.run_tag}' over {total} instance(s), workers={workers}")

        if workers == 1:
            for i, name in enumerate(names, start=1):
                self._process_one(name, i, total, skip_done)
        else:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = [pool.submit(self._process_one, name, i, total, skip_done) for i, name in enumerate(names, start=1)]
                for future in futures:
                    future.result()  # surface any unexpected (non-caught) error

        return self.store.load_all()

    def _process_one(self, name: str, i: int, total: int, skip_done: bool) -> None:
        if skip_done and self.store.is_done(name):
            self._log(f"[{i}/{total}] skip (already done): {name}")
            return
        try:
            record = self.run_instance(name)
            self.store.save(record)
            self._log(
                f"[{i}/{total}] done: {name} cost={record.get('cost')} "
                f"time={record.get('runtimeSec')}s validated={record.get('validated')}"
            )
        except Exception as exc:  # noqa: BLE001 - keep the batch going
            self._log(f"[{i}/{total}] ERROR on {name}: {exc}")
            traceback.print_exc()

    def run_instance(self, name: str) -> dict:
        """Solve a single instance and return its result record (does not persist)."""
        path = self.instances_dir / name
        if not path.exists():
            raise FileNotFoundError(f"instance not found: {path}")

        instance = CBMInstance.from_file(path)
        record = self.solver.solve(instance, path)

        return {
            "instance": name,
            "rows": instance.rows,
            "cols": instance.cols,
            "solver": self.solver.key,
            "runTag": self.run_tag,
            **record,
            "timestamp": _now_iso(),
        }

    # ----------------------------------------------------------------- helpers
    @staticmethod
    def _log(message: str) -> None:
        print(f"{_now_iso()} {LOG_PREFIX} {message}", flush=True)


def discover_instances(instances_dir: str | Path) -> list[str]:
    """Return the sorted names of all regular files in ``instances_dir``."""
    instances_dir = Path(instances_dir)
    return sorted(p.name for p in instances_dir.iterdir() if p.is_file())
