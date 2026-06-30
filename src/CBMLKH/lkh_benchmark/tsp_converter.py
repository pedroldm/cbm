"""CBM -> TSP conversion, mirroring ``LKHWrapper::writeTSP``.

The whole instance is converted (pure LKH always solves the full column set).
A dummy depot node 0 is added whose distance to column ``j`` equals
``onesCount(j)``; this turns the Hamiltonian *cycle* LKH solves into an open
Hamiltonian *path* over the columns. Columns are connected by Hamming distance.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .instance import CBMInstance


class TSPConverter:
    """Builds and serializes the ``EXPLICIT`` / ``FULL_MATRIX`` TSP for an instance."""

    def __init__(self, instance: CBMInstance) -> None:
        self.instance = instance

    def build_matrix(self) -> np.ndarray:
        """Return the ``(cols+1) x (cols+1)`` integer distance matrix.

        Row/col 0 is the depot: ``D[0, j] = D[j, 0] = onesCount(col_{j-1})``.
        The remaining block is the column-vs-column Hamming matrix.
        """
        inst = self.instance
        n = inst.cols
        D = np.zeros((n + 1, n + 1), dtype=np.int64)
        D[0, 1:] = inst.ones_count
        D[1:, 0] = inst.ones_count
        D[1:, 1:] = inst.hamming_matrix()
        return D

    def write(self, path: str | Path, name: str) -> Path:
        """Write the TSPLIB problem file and return its path."""
        D = self.build_matrix()
        header = [
            f"NAME : {name}",
            "TYPE : TSP",
            f"DIMENSION : {D.shape[0]}",
            "EDGE_WEIGHT_TYPE : EXPLICIT",
            "EDGE_WEIGHT_FORMAT : FULL_MATRIX",
            "EDGE_WEIGHT_SECTION",
        ]
        path = Path(path)
        with path.open("w") as f:
            f.write("\n".join(header) + "\n")
            np.savetxt(f, D, fmt="%d")
        return path
