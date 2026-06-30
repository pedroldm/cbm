"""CBMLKH solver: shells out to the C++ binary and parses its JSON report.

The binary takes a single ``key=value`` config file (see ``parameters.md``) and
prints one JSON document on stdout (``main.cpp``). It validates every improving
solution internally via the C++ ``Validator``, so we trust its reported
``bestCost`` rather than re-validating a permutation (the report does not include
one).

The complete report the solver emits is stored verbatim under the record's
``report`` key, so every metric it produces at the end of a run (global
accepted/rejected moves, cache stats, per-operator effectiveness and every
trajectory's metrics + improvement history) is preserved.
"""

from __future__ import annotations

import json
import subprocess
import time
import uuid
from pathlib import Path

from ..config import CBMLKH_BINARY_CANDIDATES, CBMLKH_DIR, CBMLKH_WORK_DIRNAME, DEFAULT_CBMLKH_CONFIG
from ..instance import CBMInstance
from .base import Solver

__all__ = ["CBMLKHSolver", "find_cbmlkh_binary", "DEFAULT_CBMLKH_CONFIG"]


def find_cbmlkh_binary() -> str:
    """Return the first existing CBMLKH binary in src/CBMLKH."""
    for name in CBMLKH_BINARY_CANDIDATES:
        candidate = CBMLKH_DIR / name
        if candidate.exists():
            return str(candidate)
    raise FileNotFoundError(f"No CBMLKH binary found in {CBMLKH_DIR} (build with `make prd`).")


class CBMLKHSolver(Solver):
    key = "cbmlkh"

    def __init__(
        self,
        config: dict | None = None,
        binary_path: str | None = None,
        work_dir: str | Path | None = None,
        tag: str | None = None,
        timeout: float | None = None,
    ) -> None:
        # instancePath is injected per-instance; keep it out of the base config.
        self.config = {k: v for k, v in (config or DEFAULT_CBMLKH_CONFIG).items() if k != "instancePath"}
        self.binary_path = binary_path or find_cbmlkh_binary()
        self.work_dir = Path(work_dir) if work_dir else CBMLKH_DIR / CBMLKH_WORK_DIRNAME
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self._tag = tag
        self.timeout = timeout

    def signature(self) -> str:
        if self._tag:
            return self._tag
        max_time = self.config.get("maxTime", "?")
        threads = self.config.get("threads", "?")
        movement = self.config.get("blockMovement", "?")
        return f"t{max_time}_th{threads}_{movement}"

    def solve(self, instance: CBMInstance, instance_path: Path) -> dict:
        cfg = dict(self.config)
        cfg["instancePath"] = str(instance_path)

        cfg_file = self.work_dir / f"{instance.name}_{uuid.uuid4().hex[:12]}.cfg"
        cfg_file.write_text("".join(f"{k}={v}\n" for k, v in cfg.items()))
        try:
            start = time.perf_counter()
            proc = subprocess.run(
                [self.binary_path, str(cfg_file)],
                cwd=str(CBMLKH_DIR),
                capture_output=True,
                text=True,
                timeout=self.timeout,
            )
            wall = time.perf_counter() - start
        finally:
            cfg_file.unlink(missing_ok=True)

        if proc.returncode != 0:
            raise RuntimeError(f"CBMLKH exited with {proc.returncode} on '{instance.name}'.\nstderr:\n{proc.stderr.strip()}")

        report = self._parse_json(proc.stdout)
        glob = report.get("global", {})
        runtime_ms = glob.get("runtimeMs", round(wall * 1000))
        return {
            # Lightweight summary used by the index and the comparison report.
            "cost": glob.get("bestCost"),
            "runtimeSec": round(runtime_ms / 1000.0, 4),
            "wallSec": round(wall, 4),
            "validated": None,  # validated internally by the C++ Validator at runtime
            "config": cfg,  # the input config we fed the binary
            # The solver's complete end-of-run report: instance echo, resolved
            # config, global metrics (accepted/rejected moves, lkhCache, per-operator
            # stats) and every trajectory's metrics + improvement history.
            "report": report,
        }

    @staticmethod
    def _parse_json(stdout: str) -> dict:
        try:
            return json.loads(stdout)
        except json.JSONDecodeError:
            # Be forgiving if anything leaked before the JSON document.
            start = stdout.find("{")
            if start == -1:
                raise RuntimeError("CBMLKH produced no JSON output.")
            return json.loads(stdout[start:])
