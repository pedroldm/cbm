"""
Convert a Harwell-Boeing (HB) sparse matrix file to a simple row format:

    rows cols
    count idx1 idx2 ...   (one line per row, 1-indexed)

Supports pattern-only matrices (no value array), including symmetric types
(PSA, SSA, etc.) where only the upper/lower triangle is stored — the script
mirrors the off-diagonal entries to produce the full symmetric matrix.

Usage:
    python hb_to_rowformat.py input.hb output.txt

    # explicit output path
    python hb_to_rowformat.py bcsstk29.hb result.txt
"""

import sys
import math


# ─── HB parsing helpers ───────────────────────────────────────────────────────

def parse_fortran_format(fmt: str) -> tuple[int, int]:
    """
    Parse a Fortran integer format string like '(10I8)' → (count=10, width=8).
    Only I (integer) descriptors are needed for pattern-only HB files.
    """
    fmt = fmt.strip().strip("()")
    # e.g. "10I8"
    for i, ch in enumerate(fmt):
        if ch in "IiRrEeDdFf":
            count = int(fmt[:i]) if fmt[:i] else 1
            width = int(fmt[i + 1:])
            return count, width
    raise ValueError(f"Unrecognised Fortran format: '{fmt}'")


def read_fixed_integers(lines: list[str], start_line: int,
                        count: int, width: int) -> tuple[list[int], int]:
    """
    Read `count` fixed-width integers of `width` chars, starting at
    `lines[start_line]`. Returns (values, next_line_index).
    """
    values: list[int] = []
    line_idx = start_line
    per_line = 80 // width  # max integers per line given the field width

    while len(values) < count:
        if line_idx >= len(lines):
            raise EOFError(f"File ended early; expected {count} integers, got {len(values)}")
        line = lines[line_idx]
        line_idx += 1
        # Pad to full width so the last field is always readable
        line = line.rstrip("\n").ljust(per_line * width)
        for start in range(0, per_line * width, width):
            chunk = line[start:start + width].strip()
            if chunk:
                values.append(int(chunk))
            if len(values) == count:
                break

    return values, line_idx


# ─── Main conversion ──────────────────────────────────────────────────────────

def hb_to_row_format(hb_path: str, out_path: str) -> None:
    with open(hb_path, "r") as f:
        lines = f.readlines()

    # ── Line 1: title (ignored) ───────────────────────────────────────────────
    # ── Line 2: totcrd, ptrcrd, indcrd, valcrd, rhscrd ───────────────────────
    line2 = lines[1].split()
    totcrd, ptrcrd, indcrd, valcrd = int(line2[0]), int(line2[1]), int(line2[2]), int(line2[3])

    # ── Line 3: mxtype, nrow, ncol, nnzero, neltvl ───────────────────────────
    line3        = lines[2].split()
    mxtype       = line3[0].upper()   # e.g. "PSA", "RUA", "SSA" …
    nrow         = int(line3[1])
    ncol         = int(line3[2])
    nnzero       = int(line3[3])      # stored non-zeros (may be half for symmetric)

    is_symmetric = mxtype[1] in ("S", "s", "H", "h", "Z", "z")

    # ── Line 4: format strings for ptr, ind, val, rhs ────────────────────────
    line4        = lines[3].split()
    ptr_fmt      = line4[0]
    ind_fmt      = line4[1]
    # val_fmt and rhs_fmt are ignored (pattern-only or not needed here)

    ptr_count, ptr_width = parse_fortran_format(ptr_fmt)
    ind_count, ind_width = parse_fortran_format(ind_fmt)

    # ── Data starts at line 5 (index 4) ──────────────────────────────────────
    data_start = 4

    # Column pointers: ncol + 1 values
    col_ptr, next_line = read_fixed_integers(lines, data_start, ncol + 1, ptr_width)

    # Row indices: nnzero values
    row_ind, _         = read_fixed_integers(lines, next_line,  nnzero,   ind_width)

    # ── Build adjacency sets (0-indexed internally) ───────────────────────────
    # col_ptr and row_ind are 1-indexed in HB format
    adj: list[set[int]] = [set() for _ in range(nrow)]

    for col_0 in range(ncol):
        start = col_ptr[col_0]     - 1   # convert to 0-indexed
        end   = col_ptr[col_0 + 1] - 1
        for ptr in range(start, end):
            row_0 = row_ind[ptr] - 1
            adj[row_0].add(col_0)
            if is_symmetric and row_0 != col_0:
                adj[col_0].add(row_0)   # mirror for full symmetric matrix

    # ── Write output ──────────────────────────────────────────────────────────
    with open(out_path, "w") as f:
        f.write(f"{nrow} {ncol}\n")
        for row in range(nrow):
            sorted_cols = sorted(adj[row])                     # 1-indexed output
            indices_str = " ".join(str(c + 1) for c in sorted_cols)
            f.write(f"{len(sorted_cols)} {indices_str}\n")

    print(f"Done. {nrow}×{ncol} matrix → '{out_path}'")
    print(f"  Stored non-zeros : {nnzero}")
    print(f"  Symmetric mirror : {is_symmetric}")
    total_nnz = sum(len(s) for s in adj)
    print(f"  Output non-zeros : {total_nnz}")


# ─── Entry point ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python hb_to_rowformat.py <input.hb> <output.txt>")
        sys.exit(1)

    hb_to_row_format(sys.argv[1], sys.argv[2])