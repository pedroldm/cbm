"""Pure-LKH solver: convert the whole instance to TSP and solve it with LKH."""

from __future__ import annotations

import tempfile
from pathlib import Path

from ..config import LKH_WORK_DIRNAME, LKHParams
from ..instance import CBMInstance
from ..lkh_runner import LKHRunner
from ..validator import SolutionValidator, ValidationError
from .base import Solver


class PureLKHSolver(Solver):
    key = "lkh"

    def __init__(self, params: LKHParams | None = None, work_dir: str | Path | None = None) -> None:
        self.params = params or LKHParams()
        work_dir = work_dir or Path(tempfile.gettempdir()) / LKH_WORK_DIRNAME
        self.runner = LKHRunner(self.params, work_dir)

    def signature(self) -> str:
        return self.params.signature()

    def solve(self, instance: CBMInstance, instance_path: Path) -> dict:
        perm, lkh_time = self.runner.solve(instance)
        cost = instance.path_cost(perm)

        validator = SolutionValidator(instance)
        validated = False
        validation_error = None
        try:
            block_count = validator.validate(perm, cost)
            validated = True
        except ValidationError as exc:
            block_count = instance.count_one_blocks(perm)
            validation_error = str(exc)

        return {
            "cost": cost,
            "runtimeSec": round(lkh_time, 4),
            "blockCount": block_count,
            "validated": validated,
            "validationError": validation_error,
            "solution": [int(c) for c in perm],
            "lkhParams": self.params.as_dict(),
        }
