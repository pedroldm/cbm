"""Pluggable solvers for the benchmark."""

from .base import Solver
from .cbmlkh import CBMLKH_CONFIG, CBMLKHSolver, find_cbmlkh_binary
from .pure_lkh import PureLKHSolver

__all__ = [
    "Solver",
    "PureLKHSolver",
    "CBMLKHSolver",
    "CBMLKH_CONFIG",
    "find_cbmlkh_binary",
]
