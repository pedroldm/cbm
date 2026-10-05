"""Seed derivation.

Two layers, both pure functions of logical identifiers (never of the clock,
PIDs, worker ids or completion order):

* :func:`run_seed` maps (root seed, method, instance, repetition) to the seed a
  solver process receives. It hashes a canonical JSON encoding with SHA-256, so
  it is stable across Python versions and machines.
* :func:`cbm_derive_seed` mirrors ``cbm_derive_seed`` in
  ``src/common/cbm_seed.h``, which the solvers use internally to derive the
  seed of TSP call k (ENS, ILS) or of trajectory k (CBMLKH). Keeping a Python
  copy lets the runner record, and the tests check, those inner seeds.

All seeds lie in [1, 2^31 - 1]: valid for LKH's SEED and Linkern's ``-s``
(where 0 would mean "seed from the clock").
"""

from __future__ import annotations

import hashlib
import json

SEED_MODULUS = 2**31 - 1
_MASK64 = (1 << 64) - 1


def run_seed(root_seed: int, method: str, instance: str, repetition: int) -> int:
    payload = json.dumps(["cbm-run-seed/1", int(root_seed), method, instance, int(repetition)], separators=(",", ":"))
    digest = hashlib.sha256(payload.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % SEED_MODULUS + 1


def splitmix64(x: int) -> int:
    x = (x + 0x9E3779B97F4A7C15) & _MASK64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & _MASK64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & _MASK64
    return x ^ (x >> 31)


def cbm_derive_seed(seed: int, index: int) -> int:
    h = splitmix64(splitmix64(seed & _MASK64) ^ (index & _MASK64))
    return h % SEED_MODULUS + 1
