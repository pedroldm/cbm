"""Parametrizable LKH invocation, mirroring ``LKHWrapper`` (par file + tour parse).

``LKHRunner`` writes the ``.tsp`` / ``.par``, shells out to the LKH binary and
parses the resulting tour back into a column permutation, replicating
``LKHWrapper::getResultTour`` (``node - 2`` plus depot rotation). The tunable
knobs live in :class:`~lkh_benchmark.config.LKHParams`.
"""

from __future__ import annotations

import subprocess
import time
import uuid
from pathlib import Path

from .config import DEFAULT_LKH_PATH, LKHParams  # re-exported for backwards compatibility
from .instance import CBMInstance
from .tsp_converter import TSPConverter

__all__ = ["LKHRunner", "LKHParams", "DEFAULT_LKH_PATH"]


class LKHRunner:
    """Runs LKH on a full :class:`CBMInstance` and returns ``(permutation, seconds)``."""

    def __init__(self, params: LKHParams, work_dir: str | Path) -> None:
        self.params = params
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

    def solve(self, instance: CBMInstance, keep_files: bool = False):
        exec_id = uuid.uuid4().hex[:12]
        base = f"{instance.name}_{exec_id}"
        tsp_file = self.work_dir / f"{base}.tsp"
        par_file = self.work_dir / f"{base}.par"
        tour_file = self.work_dir / f"{base}.tour"

        TSPConverter(instance).write(tsp_file, base)
        par_file.write_text("\n".join(self.params.to_par_lines(str(tsp_file), str(tour_file))) + "\n")

        start = time.perf_counter()
        self._run_lkh(par_file)
        elapsed = time.perf_counter() - start

        if not tour_file.exists():
            raise RuntimeError(f"LKH produced no tour file for instance '{instance.name}'")
        perm = self._parse_tour(tour_file)

        if not keep_files:
            for p in (tsp_file, par_file, tour_file):
                p.unlink(missing_ok=True)
        return perm, elapsed

    def _run_lkh(self, par_file: Path) -> None:
        lkh = Path(self.params.lkh_path)
        if not lkh.exists():
            raise FileNotFoundError(f"LKH executable not found: {lkh}")
        proc = subprocess.run(
            [str(lkh), str(par_file)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if proc.returncode != 0:
            # LKH can exit non-zero even when a usable tour was written; warn and
            # let tour parsing decide whether the run is salvageable.
            print(f"[lkh_benchmark] WARNING: LKH returned code {proc.returncode} for {par_file.name}")

    @staticmethod
    def _parse_tour(tour_file: Path) -> list[int]:
        """Parse a TSPLIB TOUR_SECTION into a column permutation.

        Replicates ``LKHWrapper::getResultTour``: each node is mapped via
        ``node - 2`` (so the depot becomes ``-1`` and column nodes become
        0-based column ids), then the list is rotated so the depot is first and
        the depot is dropped, yielding the open path's column order.
        """
        tour: list[int] = []
        in_section = False
        with tour_file.open() as f:
            for line in f:
                stripped = line.strip()
                if not in_section:
                    if stripped == "TOUR_SECTION":
                        in_section = True
                    continue
                done = False
                for tok in stripped.split():
                    node = int(tok)
                    if node == -1:
                        done = True
                        break
                    tour.append(node - 2)
                if done:
                    break

        if -1 in tour:
            depot = tour.index(-1)
            tour = tour[depot + 1:] + tour[:depot]
        return tour
