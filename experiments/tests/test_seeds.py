"""Seed derivation: stability, separation and agreement with the C header."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import helpers  # noqa: F401  (sets sys.path)

from cbm_experiments.seeds import SEED_MODULUS, cbm_derive_seed, run_seed

HEADER = Path(__file__).resolve().parents[2] / "src" / "common" / "cbm_seed.h"


class TestRunSeed(unittest.TestCase):
    def test_pinned_values(self):
        # Changing the derivation would silently change every experiment.
        self.assertEqual(run_seed(20261005, "ens", "a1", 0), run_seed(20261005, "ens", "a1", 0))
        self.assertEqual(cbm_derive_seed(99, 0), 594804891)  # observed from the CBMLKH binary
        self.assertEqual(cbm_derive_seed(99, 1), 1404848343)
        self.assertEqual(cbm_derive_seed(7, 0), 1880445815)  # observed from ENS/ILS (TSP call 0, seed 7)

    def test_range_and_separation(self):
        seen = set()
        for method in ("ens", "ils", "lkh", "cbmlkh"):
            for inst in ("a1", "a2", "scpa1"):
                for rep in range(-1, 10):
                    s = run_seed(1, method, inst, rep)
                    self.assertTrue(1 <= s <= SEED_MODULUS)
                    seen.add(s)
        self.assertEqual(len(seen), 4 * 3 * 11)
        self.assertNotEqual(run_seed(1, "ens", "a1", 0), run_seed(2, "ens", "a1", 0))
        calls = [cbm_derive_seed(123, k) for k in range(5000)]
        self.assertEqual(len(set(calls)), len(calls))
        self.assertTrue(all(1 <= c <= SEED_MODULUS for c in calls))


@unittest.skipUnless(shutil.which("gcc"), "gcc not available")
class TestMatchesCHeader(unittest.TestCase):
    def test_python_mirror_equals_c(self):
        cases = [(1, 0), (99, 0), (99, 1), (2**63 + 5, 17), (12345, 2**64 - 1), (0, 0)]
        src = "#include <stdio.h>\n#include \"%s\"\nint main(void){\n" % HEADER
        for seed, idx in cases:
            src += f'printf("%u\\n", cbm_derive_seed({seed}ULL, {idx}ULL));\n'
        src += "return 0;}\n"
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "t.c").write_text(src)
            subprocess.run(["gcc", "-o", f"{d}/t", f"{d}/t.c"], check=True)
            out = subprocess.run([f"{d}/t"], capture_output=True, text=True, check=True).stdout.split()
        self.assertEqual([int(x) for x in out], [cbm_derive_seed(s, i) for s, i in cases])


if __name__ == "__main__":
    unittest.main()
