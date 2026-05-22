from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
import json
import os
from pathlib import Path

import pandas as pd
import tqdm


@dataclass
class Instance:
    """
    Reads a CBM instance file and stores the binary matrix.
    """

    file_path: str | Path
    l: int = 0  # number of rows
    c: int = 0  # number of columns
    binary_matrix: list[list[bool]] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.file_path = Path(self.file_path)
        self._read()

    # -------------------------------------------------------------------------
    # Reading
    # -------------------------------------------------------------------------

    def _read(self) -> None:
        if not self.file_path.exists():
            raise FileNotFoundError(f"File not found: {self.file_path}")

        with self.file_path.open("r", encoding="utf-8") as f:
            tokens = iter(f.read().split())

            self.l = int(next(tokens))
            self.c = int(next(tokens))

            self.binary_matrix = [
                [False] * self.c
                for _ in range(self.l)
            ]

            for row in range(self.l):
                n = int(next(tokens))
                for _ in range(n):
                    col = int(next(tokens)) - 1
                    self.binary_matrix[row][col] = True

    # -------------------------------------------------------------------------
    # Counting methods
    # -------------------------------------------------------------------------

    def count_ones_in_row(self, row: int) -> int:
        return sum(self.binary_matrix[row])

    def count_ones_in_column(self, col: int) -> int:
        return sum(self.binary_matrix[row][col] for row in range(self.l))

    def count_ones_per_row(self) -> list[int]:
        return [sum(row) for row in self.binary_matrix]

    def count_ones_per_column(self) -> list[int]:
        return [
            sum(self.binary_matrix[row][col] for row in range(self.l))
            for col in range(self.c)
        ]

    # -------------------------------------------------------------------------
    # Empty rows/columns
    # -------------------------------------------------------------------------

    def get_empty_rows(self) -> list[int]:
        """Return row indices with no 1s."""
        return [
            row
            for row in range(self.l)
            if self.count_ones_in_row(row) == 0
        ]

    def get_empty_columns(self) -> list[int]:
        """Return column indices with no 1s."""
        return [
            col
            for col in range(self.c)
            if self.count_ones_in_column(col) == 0
        ]

    # -------------------------------------------------------------------------
    # Subset relationships
    # -------------------------------------------------------------------------

    def is_row_subset(self, row_a: int, row_b: int) -> bool:
        """
        Return True if row_a is a subset of row_b.
        That is, every 1 in row_a is also 1 in row_b.
        """
        return all(
            not self.binary_matrix[row_a][col]
            or self.binary_matrix[row_b][col]
            for col in range(self.c)
        )

    def is_column_subset(self, col_a: int, col_b: int) -> bool:
        """
        Return True if column_a is a subset of column_b.
        That is, every 1 in column_a is also 1 in column_b.
        """
        return all(
            not self.binary_matrix[row][col_a]
            or self.binary_matrix[row][col_b]
            for row in range(self.l)
        )

    # -------------------------------------------------------------------------
    # New methods to identify rows/columns with exactly one 1
    # -------------------------------------------------------------------------

    def get_singleton_rows(self) -> list[int]:
        """
        Return row indices that contain exactly one 1.
        """
        row_counts = self.count_ones_per_row()
        return [
            row
            for row, count in enumerate(row_counts)
            if count == 1
        ]


    def get_singleton_columns(self) -> list[int]:
        """
        Return column indices that contain exactly one 1.
        """
        col_counts = self.count_ones_per_column()
        return [
            col
            for col, count in enumerate(col_counts)
            if count == 1
        ]


    # -------------------------------------------------------------------------
    # Updated subset methods
    # Exclude rows/columns with <= 1 ones from subset analysis.
    # -------------------------------------------------------------------------

    def get_row_subsets(self) -> list[tuple[int, int]]:
        """
        Return all pairs (a, b) where row a is a subset of row b.

        Rows with 0 or 1 ones are ignored.
        """
        subsets = []
        row_counts = self.count_ones_per_row()

        valid_rows = [
            row
            for row, count in enumerate(row_counts)
            if count > 1
        ]

        for a in valid_rows:
            for b in valid_rows:
                if a == b:
                    continue

                # Proper subset only
                if (
                    row_counts[a] < row_counts[b]
                    and self.is_row_subset(a, b)
                ):
                    subsets.append((a, b))

        return subsets


    def get_column_subsets(self) -> list[tuple[int, int]]:
        """
        Return all pairs (a, b) where column a is a subset of column b.

        Columns with 0 or 1 ones are ignored.
        """
        subsets = []
        col_counts = self.count_ones_per_column()

        valid_cols = [
            col
            for col, count in enumerate(col_counts)
            if count > 1
        ]

        for a in valid_cols:
            for b in valid_cols:
                if a == b:
                    continue

                # Proper subset only
                if (
                    col_counts[a] < col_counts[b]
                    and self.is_column_subset(a, b)
                ):
                    subsets.append((a, b))

        return subsets

        # -------------------------------------------------------------------------
        # Utility methods
        # -------------------------------------------------------------------------

    def print_matrix(self) -> None:
        for row in self.binary_matrix:
            print(" ".join("1" if cell else "0" for cell in row))

    def total_ones(self) -> int:
        return sum(sum(row) for row in self.binary_matrix)

    def __repr__(self) -> str:
        return (
            f"Instance(file='{self.file_path.name}', "
            f"rows={self.l}, cols={self.c}, total_ones={self.total_ones()})"
        )

def process_instance(instance_path: str) -> dict:
    instance_obj = Instance(instance_path)

    # ------------------------------------------------------------------
    # Compute subset pairs
    # Each pair (a, b) means: column/row a is a subset of column/row b
    # ------------------------------------------------------------------
    row_subset_pairs = instance_obj.get_row_subsets()
    col_subset_pairs = instance_obj.get_column_subsets()

    # ------------------------------------------------------------------
    # Precompute indices of 1s
    # ------------------------------------------------------------------
    row_ones = {
        row: [
            col
            for col, value in enumerate(instance_obj.binary_matrix[row])
            if value
        ]
        for row in range(instance_obj.l)
    }

    col_ones = {
        col: [
            row
            for row in range(instance_obj.l)
            if instance_obj.binary_matrix[row][col]
        ]
        for col in range(instance_obj.c)
    }

    # ------------------------------------------------------------------
    # Build detailed column subset information
    # ------------------------------------------------------------------
    column_subset_details = []

    for subset_col, superset_col in col_subset_pairs:
        column_subset_details.append(
            {
                "subset_column": subset_col,
                "superset_column": superset_col,
                "subset_ones": col_ones[subset_col],
                "superset_ones": col_ones[superset_col],
            }
        )

    # ------------------------------------------------------------------
    # Build detailed row subset information
    # ------------------------------------------------------------------
    row_subset_details = []

    for subset_row, superset_row in row_subset_pairs:
        row_subset_details.append(
            {
                "subset_row": subset_row,
                "superset_row": superset_row,
                "subset_ones": row_ones[subset_row],
                "superset_ones": row_ones[superset_row],
            }
        )

    # ------------------------------------------------------------------
    # Final result
    # ------------------------------------------------------------------
    singleton_row_indices = instance_obj.get_singleton_rows()
    singleton_col_indices = instance_obj.get_singleton_columns()

# ------------------------------------------------------------------
# Final result
# ------------------------------------------------------------------
    return {
        "filename": os.path.basename(instance_path),
        "num_rows": instance_obj.l,
        "num_cols": instance_obj.c,
        "total_ones": instance_obj.total_ones(),

        # Empty rows/columns
        "empty_rows": len(instance_obj.get_empty_rows()),
        "empty_cols": len(instance_obj.get_empty_columns()),

        # Singleton rows/columns
        "singleton_rows": len(singleton_row_indices),
        "singleton_cols": len(singleton_col_indices),
        "singleton_row_indices": singleton_row_indices,
        "singleton_col_indices": singleton_col_indices,

        # Subset counts
        "row_subsets": len(row_subset_pairs),
        "col_subsets": len(col_subset_pairs),

        # Detailed subset information
        "row_subset_details": row_subset_details,
        "column_subset_details": column_subset_details,
    }

def compute_or_load_metrics(
    instance_dir: str,
    cache_file: str = "instance_metrics.json",
    max_workers: int | None = None,
) -> list[dict]:
    cache_path = Path(cache_file)

    # Load cached results if they already exist
    if cache_path.exists():
        print(f"Loading cached metrics from {cache_path}")
        with cache_path.open("r", encoding="utf-8") as f:
            return json.load(f)

    # Collect all instance file paths
    instance_paths = sorted(
        os.path.join(instance_dir, filename)
        for filename in os.listdir(instance_dir)
        if os.path.isfile(os.path.join(instance_dir, filename))
    )

    if max_workers is None:
        max_workers = os.cpu_count()

    results = []

    # Process in parallel using multiple processes
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(process_instance, path)
            for path in instance_paths
        ]

        for future in tqdm.tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Processing instances",
        ):
            results.append(future.result())

    # Keep deterministic ordering
    results.sort(key=lambda x: x["filename"])

    # Save cache
    with cache_path.open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print(f"Metrics saved to {cache_path}")

    return results


