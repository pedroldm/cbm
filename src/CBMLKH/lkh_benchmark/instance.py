"""CBM instance representation and the column-pair quantities used everywhere.

Mirrors the C++ ``ColumnStore`` (bitset column statistics) and ``Validator``
(dense binary matrix + 1-block count) so that TSP conversion, cost evaluation
and validation all share a single source of truth.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


class CBMInstance:
    """A Consecutive Block Minimization instance.

    Sparse row format (same as ``instances/``)::

        <rows> <cols>
        <k> <c1> <c2> ... <ck>      # k 1-indexed column positions holding a 1
        ...

    Attributes:
        name:   instance file name.
        matrix: dense ``rows x cols`` uint8 binary matrix.
        rows:   number of matrix rows (l).
        cols:   number of columns / TSP cities (c).
    """

    def __init__(self, name: str, matrix: np.ndarray) -> None:
        if matrix.ndim != 2:
            raise ValueError("matrix must be 2-D (rows x cols)")
        self.name = name
        self.matrix = np.ascontiguousarray(matrix, dtype=np.uint8)
        self.rows, self.cols = self.matrix.shape
        self._ones_count = self.matrix.sum(axis=0, dtype=np.int64)  # per column
        self._ones_to_zeros: np.ndarray | None = None  # lazy cols x cols

    # ------------------------------------------------------------------ parsing
    @classmethod
    def from_file(cls, path: str | Path) -> "CBMInstance":
        path = Path(path)
        with path.open() as f:
            tokens = f.read().split()

        it = iter(tokens)
        try:
            rows = int(next(it))
            cols = int(next(it))
        except StopIteration as exc:
            raise ValueError(f"Malformed instance header: {path}") from exc

        matrix = np.zeros((rows, cols), dtype=np.uint8)
        for r in range(rows):
            try:
                k = int(next(it))
            except StopIteration as exc:
                raise ValueError(f"Malformed instance body (row {r}): {path}") from exc
            for _ in range(k):
                c = int(next(it))
                if c < 1 or c > cols:
                    raise ValueError(f"Column index {c} out of range in: {path}")
                matrix[r, c - 1] = 1
        return cls(path.name, matrix)

    # ------------------------------------------------------- column statistics
    @property
    def ones_count(self) -> np.ndarray:
        """Per-column number of rows holding a 1 (``ColumnStore::onesCount``)."""
        return self._ones_count

    @property
    def ones_to_zeros(self) -> np.ndarray:
        """``C[i, j]`` = #rows where column i == 1 and column j == 0.

        Equivalent to ``ColumnStore::onesToZeros(i, j)``. Built once via a single
        matrix product (``onesToOnes = A.T @ A``; ``onesToZeros = onesCount_i - onesToOnes``).
        """
        if self._ones_to_zeros is None:
            A = self.matrix.astype(np.int64)  # rows x cols
            ones_to_ones = A.T @ A  # cols x cols
            self._ones_to_zeros = self._ones_count[:, None] - ones_to_ones
        return self._ones_to_zeros

    def hamming_matrix(self) -> np.ndarray:
        """``cols x cols`` symmetric Hamming distances (``ColumnStore::hamming``).

        ``hamming(i, j) = onesToZeros(i, j) + zerosToOnes(i, j) = C[i, j] + C[j, i]``.
        """
        C = self.ones_to_zeros
        return C + C.T

    # -------------------------------------------------------------- evaluation
    def path_cost(self, perm) -> int:
        """CBM cost via the ``completeEval`` formula.

        ``cost = onesCount[perm[0]] + Σ zerosToOnes(perm[i-1], perm[i])`` where
        ``zerosToOnes(a, b) = onesToZeros(b, a) = C[b, a]``.
        """
        perm = np.asarray(perm, dtype=np.int64)
        if perm.size == 0:
            return 0
        C = self.ones_to_zeros
        cost = int(self._ones_count[perm[0]])
        if perm.size > 1:
            cost += int(C[perm[1:], perm[:-1]].sum())  # zerosToOnes(prev, cur) = C[cur, prev]
        return cost

    def count_one_blocks(self, perm) -> int:
        """Independent 1-block count (``Validator::countOneBlocks``), vectorized.

        A 1-block starts at a column whose value is 1 and whose predecessor (or
        the array start) is 0. Computed directly on the reordered binary matrix,
        so it is a genuine cross-check of :meth:`path_cost`.
        """
        perm = np.asarray(perm, dtype=np.int64)
        if perm.size == 0:
            return 0
        reordered = self.matrix[:, perm]
        blocks = int(reordered[:, 0].sum())  # blocks touching the left edge
        if reordered.shape[1] > 1:
            starts = (reordered[:, 1:] == 1) & (reordered[:, :-1] == 0)
            blocks += int(starts.sum())
        return blocks
