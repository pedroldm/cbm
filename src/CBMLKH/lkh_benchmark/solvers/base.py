"""Solver abstraction so different methods plug into the same benchmark loop.

A :class:`Solver` knows how to produce a solution record for one instance. The
benchmark handles instance loading, resume and persistence; the solver only
implements :meth:`solve` and a :meth:`signature` used to name its results dir.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

from ..instance import CBMInstance


class Solver(ABC):
    #: short, stable identifier (e.g. "lkh", "cbmlkh"); set by each subclass.
    key: str = "solver"

    @abstractmethod
    def signature(self) -> str:
        """A filesystem-safe tag identifying this solver's parameters."""

    @abstractmethod
    def solve(self, instance: CBMInstance, instance_path: Path) -> dict:
        """Solve one instance and return a result record.

        The record must contain at least ``cost`` and ``runtimeSec``; anything
        else (solution, validation flags, solver-specific metrics) is preserved
        as-is. The benchmark adds ``instance``, ``rows``, ``cols``, ``solver``,
        ``runTag`` and ``timestamp``.
        """

    def run_tag(self) -> str:
        return f"{self.key}_{self.signature()}"
