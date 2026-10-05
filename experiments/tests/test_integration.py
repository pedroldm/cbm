"""End-to-end with the real solvers on small instances.

Skipped unless the binaries exist (``make`` at the repository root) and a
Linkern executable is given by $LINKERN_PATH (or <repo>/bin/linkern).
"""

import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

from helpers import EXPERIMENTS, available_cpus

from cbm_experiments import cli
from cbm_experiments.schema import validate
from cbm_experiments.seeds import cbm_derive_seed

REPO = EXPERIMENTS.parent
BINARIES = [REPO / "ENS" / "ENS", REPO / "ILS" / "ils", REPO / "src" / "CBMLKH" / "main_prd", REPO / "src" / "CBMLKH" / "lkh_standalone", REPO / "src" / "LKH3" / "LKH"]
LINKERN = os.environ.get("LINKERN_PATH", str(REPO / "bin" / "linkern"))
READY = all(os.access(b, os.X_OK) for b in BINARIES) and os.access(LINKERN, os.X_OK) and len(available_cpus(4)) == 4

# Short runs; CBMLKH's tuned parameters are locked, so shortening it needs the
# testing-only --allow-param-overrides.
SMALL_PARAMS = {
    "ens": {"iterations": 5},
    "ils": {"iterations": 3},
    "cbmlkh": {"maxIterations": 15},
}


@unittest.skipUnless(READY, "real solver binaries / Linkern / 4 CPUs not available")
class TestRealSolvers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="cbm-integration-"))
        cls.instances = cls.tmp / "instances"
        cls.instances.mkdir()
        for name in ("a1", "d1"):
            shutil.copy(REPO / "instances" / name, cls.instances / name)
        (cls.tmp / "params.json").write_text(json.dumps(SMALL_PARAMS))
        cls.outputs = cls.tmp / "outputs"

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def args(self, name, cpus, extra=()):
        return ["run", "--outputs", str(self.outputs), "--name", name, "--instance-dir", str(self.instances),
                "--cpus", ",".join(map(str, available_cpus(cpus))), "--linkern-path", LINKERN, "--cbmlkh-lock", str(self.tmp / "cbmlkh.lock"),
                "--methods", "ens,ils,lkh,cbmlkh", "--repetitions", "2", "--cbmlkh-threads", "2", "--seed", "11",
                "--soft-limit", "120", "--hard-limit", "600", "--method-params", str(self.tmp / "params.json"), "--allow-param-overrides", "--poll-interval", "0.2", *extra]

    def records(self, name):
        out = {}
        for p in (self.outputs / name / "runs").glob("*/*__r??.json"):
            r = json.loads(p.read_text())
            sol = json.loads((p.parent / r["solution"]["permutation_file"]).read_text())["permutation"]
            out[(r["method"]["id"], r["instance"]["name"], r["repetition"])] = (r, sol)
        return out

    def test_parallel_equals_serial_and_resume(self):
        # Parallel run, interrupted after 3 launches and then resumed.
        self.assertEqual(cli.main(self.args("par", 4, ["--max-jobs", "3"])), 0)
        self.assertEqual(cli.main(self.args("par", 4)), 0)
        self.assertEqual(cli.main(self.args("ser", 2)), 0)
        par, ser = self.records("par"), self.records("ser")
        self.assertEqual(len(par), 2 * 4 * 2)
        self.assertEqual(par.keys(), ser.keys())
        for key, (rec, sol) in par.items():
            self.assertEqual(validate(rec), [], key)
            self.assertEqual(rec["status"], "completed", (key, rec["error"]))
            self.assertTrue(rec["objective"]["validated"])
            self.assertEqual(rec["objective"]["value"], ser[key][0]["objective"]["value"], key)
            self.assertEqual(sol, ser[key][1], f"{key}: same seed, different solution")
            derived = rec["seeds"]["derived"]
            if key[0] in ("ens", "ils"):
                seeds = derived["tsp_call_seeds"]
                self.assertEqual(len(set(seeds)), len(seeds), "two TSP calls shared a seed")
            if key[0] == "ens":
                self.assertEqual(derived["tsp_call_seeds"], [cbm_derive_seed(rec["seeds"]["run"], k) for k in range(len(derived["tsp_call_seeds"]))])
            if key[0] == "cbmlkh":
                self.assertEqual(derived["trajectory_seed"], cbm_derive_seed(rec["seeds"]["run"], key[2]))
                self.assertEqual(rec["tsp_solver"]["time_limit_hits"], 0, "LKH time limit bound: replay not guaranteed")
        # Different repetitions really are different runs.
        self.assertNotEqual(par[("lkh", "d1", 0)][0]["seeds"]["run"], par[("lkh", "d1", 1)][0]["seeds"]["run"])
        states = [json.loads(p.read_text()) for p in (self.outputs / "par" / "runs").glob("*/state.json")]
        self.assertTrue(all(s["status"] == "completed" and s["attempts"] == 1 for s in states), "resume re-ran a job")
        self.assertFalse(any((self.outputs / "par" / "work").iterdir()), "scratch directories were left behind")


if __name__ == "__main__":
    unittest.main()
