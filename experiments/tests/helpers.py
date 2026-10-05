"""Shared fixtures: temporary instances, fake solvers and runner invocations."""

import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

EXPERIMENTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXPERIMENTS))

from cbm_experiments import cli  # noqa: E402

FAKE = str(Path(__file__).resolve().parent / "fake_solver.py")
RUN_SCRIPT = str(EXPERIMENTS / "run_experiments.py")

TINY = {
    # rows, then each row's 1-based columns
    "t1": (8, [[1, 3], [2, 4, 5], [1, 6], [7, 8], [3, 5, 7]]),
    "t2": (6, [[1, 2], [3, 6], [2, 5], [4, 6], [1, 4]]),
    "t3": (7, [[1, 7], [2, 3], [4, 5, 6], [1, 2, 3], [6, 7]]),
}


def available_cpus(n):
    cpus = sorted(os.sched_getaffinity(0))
    return cpus[:n]


class RunnerTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="cbm-runner-test-"))
        self.instances = self.tmp / "instances"
        self.instances.mkdir()
        for name, (cols, rows) in TINY.items():
            lines = [f"{len(rows)} {cols}"] + [" ".join(map(str, [len(r)] + r)) for r in rows]
            (self.instances / name).write_text("\n".join(lines) + "\n")
        (self.instances / "scratch_dir").mkdir()  # must be ignored by discovery
        self.outputs = self.tmp / "outputs"
        self.log_dir = self.tmp / "fake_log"
        self.log_dir.mkdir()
        self.control = self.tmp / "control.json"
        self.control.write_text("{}")
        self.lock = self.tmp / "cbmlkh.lock"
        self.env_backup = dict(os.environ)
        os.environ["FAKE_SOLVER_CONTROL"] = str(self.control)
        os.environ["FAKE_SOLVER_LOG"] = str(self.log_dir)

    def tearDown(self):
        os.environ.clear()
        os.environ.update(self.env_backup)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def set_control(self, data):
        self.control.write_text(json.dumps(data))

    def base_args(self, name="exp", cpus=2, extra=()):
        return [
            "run", "--outputs", str(self.outputs), "--name", name,
            "--instance-dir", str(self.instances),
            "--cpus", ",".join(map(str, available_cpus(cpus))),
            "--ens-bin", FAKE, "--ils-bin", FAKE, "--cbmlkh-bin", FAKE, "--lkh-standalone-bin", FAKE,
            "--lkh-path", "/bin/true", "--linkern-path", "/bin/true",
            "--cbmlkh-lock", str(self.lock),
            "--poll-interval", "0.05", "--kill-grace", "1",
            *extra,
        ]

    def run_cli(self, args):
        return cli.main(args)

    def spawn_cli(self, args):
        return subprocess.Popen([sys.executable, RUN_SCRIPT, *args], env=dict(os.environ), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    def exp_dir(self, name="exp"):
        return self.outputs / name

    def fake_logs(self):
        out = []
        for p in self.log_dir.glob("*.json"):
            try:
                out.append(json.loads(p.read_text()))
            except ValueError:
                pass
        return out

    def states(self, name="exp"):
        exp = json.loads((self.exp_dir(name) / "experiment.json").read_text())
        return {j["job_id"]: json.loads((self.exp_dir(name) / "runs" / j["job_id"] / "state.json").read_text())
                for j in exp["jobs"] if (self.exp_dir(name) / "runs" / j["job_id"] / "state.json").exists()}, exp

    def records(self, name="exp"):
        return [json.loads(p.read_text()) for p in sorted((self.exp_dir(name) / "runs").glob("*/*__r??.json"))]

    def read_csv(self, name, which):
        with open(self.exp_dir(name) / "results" / f"{which}.csv", newline="") as f:
            return list(csv.DictReader(f))

    @staticmethod
    def pid_alive(pid):
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        try:  # a zombie is dead for our purposes
            return Path(f"/proc/{pid}/stat").read_text().split(")")[1].split()[0] != "Z"
        except OSError:
            return False

    def wait_for(self, predicate, timeout=30.0):
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            if predicate():
                return True
            time.sleep(0.05)
        return False
