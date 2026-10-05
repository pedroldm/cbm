"""CBMLKH runs with the irace-tuned configuration and the fixed LKH parameters."""

import re
import unittest
from pathlib import Path

import helpers  # noqa: F401  (sets sys.path)

from cbm_experiments import methods as M
from cbm_experiments.instances import InstanceInfo
from cbm_experiments.plan import build_jobs, make_spec

REPO = Path(__file__).resolve().parents[2]
IRACE_OUTPUT = REPO / "tunning" / "output.txt"
LKH_HEADER = REPO / "src" / "common" / "cbm_lkh_params.h"


def last_best_configuration(text: str) -> dict:
    """The configuration of the last 'Best-so-far configuration' block irace printed."""
    lines = text.splitlines()
    idx = max(i for i, line in enumerate(lines) if line.startswith("Best-so-far configuration:"))
    header = lines[idx + 2].split()
    values = lines[idx + 3].split()
    row = dict(zip(header, values[1:]))  # values[0] is the row name (the id)
    row.pop(".PARENT.", None)
    return row


def spec(soft=7200.0, overrides=None, unlock=False):
    inst = InstanceInfo("x", "/x", 10, 100, 1, "0" * 64)
    return make_spec(1, ["lkh", "cbmlkh"], {"lkh": 1, "cbmlkh": 2}, soft, 18000.0, 2, overrides or {}, [inst], "/", True, unlock)


@unittest.skipUnless(IRACE_OUTPUT.exists(), "tunning/output.txt not present")
class TestIraceConfiguration(unittest.TestCase):
    def test_pinned_config_is_irace_best(self):
        best = last_best_configuration(IRACE_OUTPUT.read_text())
        self.assertEqual(best[".ID."], "135")
        for key, value in M.CBMLKH_IRACE_CONFIG.items():
            expected = best[key]
            self.assertEqual(str(value) if isinstance(value, str) else float(value), expected if isinstance(value, str) else float(expected), key)
        # Everything irace tuned is pinned, except the deliberately derived ones.
        self.assertEqual(set(best) - set(M.CBMLKH_IRACE_CONFIG) - {".ID."}, {"threads", "maxTime", "lkhMaxTime"})


class TestCbmlkhConfig(unittest.TestCase):
    def test_config_file_carries_tuned_values_and_lkh_time_fraction(self):
        s = spec()
        self.assertEqual(s["method_params"]["cbmlkh"]["lkhMaxTime"], 1440)
        job = [j for j in build_jobs(s) if j["method"] == "cbmlkh"][0]
        text = M.cbmlkh_config_text(job, s, {"lkh": "/lkh"}, "/out.json", "/work")
        cfg = dict(line.split("=", 1) for line in text.splitlines())
        self.assertEqual(cfg["lkhMaxTime"], "1440")
        self.assertEqual(cfg["maxTime"], "7200")
        self.assertEqual(cfg["maxIterations"], "1000")
        for key, value in M.CBMLKH_IRACE_CONFIG.items():
            if key != "maxIterations":
                self.assertEqual(cfg[key], str(value), key)
        self.assertFalse([k for k in cfg if k.startswith("lkh_") or k == "tsp_backend"])

    def test_lkh_time_follows_the_budget(self):
        self.assertEqual(spec(soft=3600)["method_params"]["cbmlkh"]["lkhMaxTime"], 720)

    def test_locked_parameters_cannot_be_overridden(self):
        with self.assertRaises(ValueError):
            spec(overrides={"cbmlkh": {"neighborBias": 1.0}})
        with self.assertRaises(ValueError):
            spec(overrides={"lkh": {"move_type": 3}})
        unlocked = spec(overrides={"cbmlkh": {"maxIterations": 5}}, unlock=True)
        self.assertTrue(unlocked["partial"])
        self.assertTrue(unlocked["locked_params_overridden"])

    def test_lkh_parameters_match_the_compiled_header(self):
        header = LKH_HEADER.read_text()
        for key, value in M.LKH_PARAMS.items():
            self.assertRegex(header, rf"#define CBM_LKH_{key.upper()} {value}\b")
        self.assertEqual(M.LKH_PARAMS, {"move_type": 5, "patching_c": 3, "patching_a": 2})


if __name__ == "__main__":
    unittest.main()
