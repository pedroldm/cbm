"""CBM instance discovery and independent solution validation.

Instance format: first line ``<rows> <cols>``; then one line per row holding a
count ``k`` followed by ``k`` 1-indexed column positions of the row's 1s.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .util import sha256_file


@dataclass(frozen=True)
class InstanceInfo:
    name: str
    path: str
    rows: int
    cols: int
    size_bytes: int
    sha256: str

    def as_dict(self) -> dict:
        return asdict(self)


def read_header(path: Path) -> Optional[Tuple[int, int]]:
    """(rows, cols) if the file starts like a CBM instance, else None."""
    try:
        with open(path, "rb") as f:
            first = f.readline(256).split()
    except OSError:
        return None
    if len(first) != 2:
        return None
    try:
        rows, cols = int(first[0]), int(first[1])
    except ValueError:
        return None
    return (rows, cols) if rows > 0 and cols > 0 else None


def discover_instances(directory: Path, names: Optional[Sequence[str]] = None) -> Tuple[List[InstanceInfo], List[str]]:
    """Regular files directly under `directory` that parse as instances.

    Returns (instances sorted by name, skipped entries with the reason). Only
    the top level is scanned, so scratch subdirectories are never mistaken for
    instances.
    """
    directory = Path(directory)
    wanted = set(names) if names else None
    found, skipped = [], []
    for entry in sorted(directory.iterdir(), key=lambda p: p.name):
        if entry.name.startswith("."):
            continue
        if wanted is not None and entry.name not in wanted:
            continue
        if not entry.is_file():
            skipped.append(f"{entry.name}: not a regular file")
            continue
        header = read_header(entry)
        if header is None:
            skipped.append(f"{entry.name}: no '<rows> <cols>' header")
            continue
        found.append(
            InstanceInfo(
                name=entry.name,
                path=str(entry.resolve()),
                rows=header[0],
                cols=header[1],
                size_bytes=entry.stat().st_size,
                sha256=sha256_file(entry),
            )
        )
    if wanted is not None:
        missing = wanted - {i.name for i in found}
        if missing:
            raise FileNotFoundError(f"instances not found or invalid in {directory}: {sorted(missing)}")
    return found, skipped


class CBMInstance:
    """Row-wise sparse matrix, used to recount 1-blocks independently of the solvers."""

    def __init__(self, path: Path):
        with open(path) as f:
            tokens = f.read().split()
        rows, cols = int(tokens[0]), int(tokens[1])
        pos = 2
        self.rows_ones: List[List[int]] = []
        for r in range(rows):
            k = int(tokens[pos])
            pos += 1
            row = [int(t) for t in tokens[pos : pos + k]]
            pos += k
            if len(row) != k or any(c < 1 or c > cols for c in row):
                raise ValueError(f"{path}: malformed row {r + 1}")
            self.rows_ones.append(row)
        self.n_rows, self.n_cols = rows, cols

    def count_blocks(self, permutation: Sequence[int]) -> int:
        """Number of maximal runs of consecutive 1s, `permutation` 1-indexed."""
        position = [0] * (self.n_cols + 1)
        for idx, col in enumerate(permutation):
            position[col] = idx
        blocks = 0
        for row in self.rows_ones:
            if not row:
                continue
            ordered = sorted({position[c] for c in row})
            blocks += 1 + sum(1 for a, b in zip(ordered, ordered[1:]) if b != a + 1)
        return blocks

    def check_permutation(self, permutation: Sequence[int]) -> Optional[str]:
        """None if valid, else a description of the defect."""
        if len(permutation) != self.n_cols:
            return f"permutation has {len(permutation)} entries, expected {self.n_cols}"
        if sorted(permutation) != list(range(1, self.n_cols + 1)):
            return "permutation is not a bijection on 1..cols"
        return None


class InstanceCache:
    """Keeps the most recently used parsed instances (they can be ~40 MB of text)."""

    def __init__(self, capacity: int = 4):
        self.capacity = capacity
        self._items: Dict[str, CBMInstance] = {}

    def get(self, path: str) -> CBMInstance:
        if path in self._items:
            inst = self._items.pop(path)
        else:
            inst = CBMInstance(Path(path))
            while len(self._items) >= self.capacity:
                self._items.pop(next(iter(self._items)))
        self._items[path] = inst
        return inst
