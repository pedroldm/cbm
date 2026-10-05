"""Runner behaviour with fake solvers: statuses, timeouts, exclusion, resources, resume, CSVs."""

import json
import os
import signal
import statistics
import time
import unittest

from helpers import RunnerTestCase, available_cpus

from cbm_experiments.plan import run_id
from cbm_experiments.schema import validate
from cbm_experiments.seeds import cbm_derive_seed, run_seed


def overlaps(intervals):
    intervals = sorted(intervals)
    return any(b[0] < a[1] for a, b in zip(intervals, intervals[1:]))


class TestStatusesAndSchema(RunnerTestCase):
    def test_every_outcome_yields_a_schema_valid_record(self):
        self.set_control({
            "t1": {"mode": "ok", "sleep": 1.5},          # crosses the soft limit, must complete
            "t2": {"mode": "fail"},
            "t3": {"mode": "ignore_term", "sleep": 60},  # needs SIGKILL after the grace period
        })
        rc = self.run_cli(self.base_args(cpus=3, extra=["--methods", "ens", "--repetitions", "1", "--soft-limit", "1", "--hard-limit", "3"]))
        self.assertEqual(rc, 0)
        recs = {r["instance"]["name"]: r for r in self.records()}
        self.assertEqual(recs["t1"]["status"], "completed")
        self.assertTrue(recs["t1"]["timing"]["soft_limit_reached"], "soft limit crossed but not recorded")
        self.assertGreaterEqual(recs["t1"]["timing"]["wall_time_s"], 1.5)
        self.assertEqual(recs["t2"]["status"], "failed")
        self.assertIn("fake solver failure", recs["t2"]["error"]["stderr_tail"])
        self.assertIsNone(recs["t2"]["objective"]["value"])
        self.assertEqual(recs["t3"]["status"], "hard_timeout")
        self.assertLess(recs["t3"]["timing"]["wall_time_s"], 3 + 1 + 2, "process outlived hard limit + grace")
        for r in recs.values():
            self.assertEqual(validate(r), [], r["run_id"])

    def test_invalid_outputs_are_failures_not_successes(self):
        self.set_control({"t1": {"mode": "bad_output"}, "t2": {"mode": "wrong_blocks"}, "t3": {"mode": "ok"}})
        self.run_cli(self.base_args(cpus=3, extra=["--methods", "ils", "--repetitions", "1"]))
        recs = {r["instance"]["name"]: r for r in self.records()}
        self.assertEqual(recs["t1"]["status"], "failed")
        self.assertEqual(recs["t1"]["error"]["type"], "invalid_output")
        self.assertEqual(recs["t2"]["status"], "failed")
        self.assertIn("independent recount", recs["t2"]["error"]["message"])
        self.assertEqual(recs["t3"]["status"], "completed")
        self.assertTrue(recs["t3"]["objective"]["validated"])

    def test_soft_limit_does_not_kill(self):
        self.set_control({"t1": {"sleep": 2.0}})
        self.run_cli(self.base_args(extra=["--instances", "t1", "--methods", "lkh", "--repetitions", "1", "--soft-limit", "0.5", "--hard-limit", "10"]))
        (rec,) = self.records()
        self.assertEqual(rec["status"], "completed")
        self.assertTrue(rec["timing"]["soft_limit_reached"])
        states, _ = self.states()
        events = [e["event"] for e in next(iter(states.values()))["events"]]
        self.assertIn("soft_limit_reached", events)
        self.assertNotIn("hard_limit_reached", events)

    def test_hard_timeout_kills_descendants(self):
        self.set_control({"t1": {"sleep": 60, "grandchild": 120}})
        self.run_cli(self.base_args(extra=["--instances", "t1", "--methods", "ens", "--repetitions", "1", "--soft-limit", "0.5", "--hard-limit", "1.5"]))
        (rec,) = self.records()
        self.assertEqual(rec["status"], "hard_timeout")
        (entry,) = self.fake_logs()
        self.assertFalse(self.pid_alive(entry["pid"]))
        self.assertTrue(self.wait_for(lambda: not self.pid_alive(entry["grandchild"]), 5), "grandchild survived the group kill")


class TestSeedsAndIsolation(RunnerTestCase):
    def test_seeds_are_logical_and_distinct(self):
        self.run_cli(self.base_args(cpus=3, extra=["--methods", "ens,ils,lkh,cbmlkh", "--repetitions", "3", "--cbmlkh-threads", "2", "--seed", "77"]))
        recs = self.records()
        self.assertEqual(len(recs), 3 * 4 * 3)
        seeds = {}
        for r in recs:
            m, inst, rep = r["method"]["id"], r["instance"]["name"], r["repetition"]
            if m == "cbmlkh":
                base = run_seed(77, m, inst, -1)
                self.assertEqual(r["seeds"]["run"], base)
                self.assertEqual(r["seeds"]["derived"]["trajectory_seed"], cbm_derive_seed(base, rep))
                seeds[(m, inst, rep)] = r["seeds"]["derived"]["trajectory_seed"]
            else:
                self.assertEqual(r["seeds"]["run"], run_seed(77, m, inst, rep))
                seeds[(m, inst, rep)] = r["seeds"]["run"]
        self.assertEqual(len(set(seeds.values())), len(seeds), "two repetitions share a seed")

    def test_concurrent_runs_use_private_work_dirs_and_are_reproducible(self):
        args = ["--methods", "ens,ils,lkh", "--repetitions", "4", "--seed", "5"]
        self.run_cli(self.base_args(name="parallel", cpus=4, extra=args))
        self.run_cli(self.base_args(name="serial", cpus=1, extra=args))
        logs = self.fake_logs()
        par = [e for e in logs if "/parallel/" in e["workdir"]]
        self.assertEqual(len(par), 3 * 3 * 4)
        self.assertEqual(len({e["workdir"] for e in par}), len(par), "two runs shared a work directory")
        self.assertTrue(overlaps([(e["start"], e["end"]) for e in par]), "jobs did not actually run concurrently")

        def key(r):
            return (r["method"]["id"], r["instance"]["name"], r["repetition"])

        a = {key(r): r for r in self.records("parallel")}
        b = {key(r): r for r in self.records("serial")}
        self.assertEqual(a.keys(), b.keys())
        for k in a:
            self.assertEqual(a[k]["objective"]["value"], b[k]["objective"]["value"], k)
            sol_a = json.loads((self.exp_dir("parallel") / "runs" / a[k]["job_id"] / a[k]["solution"]["permutation_file"]).read_text())
            sol_b = json.loads((self.exp_dir("serial") / "runs" / b[k]["job_id"] / b[k]["solution"]["permutation_file"]).read_text())
            self.assertEqual(sol_a["permutation"], sol_b["permutation"], k)


class TestResources(RunnerTestCase):
    def test_cpu_budget_and_pinning(self):
        cpus = available_cpus(3)
        self.set_control({n: {"sleep": 0.4} for n in ("t1", "t2", "t3")})
        self.run_cli(self.base_args(cpus=3, extra=["--methods", "ens,ils,cbmlkh", "--repetitions", "4", "--cbmlkh-threads", "2"]))
        logs = self.fake_logs()
        events = sorted([(e["start"], 1, e) for e in logs] + [(e["end"], -1, e) for e in logs], key=lambda x: (x[0], x[1]))
        busy, peak, active = 0, 0, []
        for _, delta, e in events:
            width = len(e["affinity"])
            busy += delta * width
            peak = max(peak, busy)
            if delta > 0:
                taken = {c for a in active for c in a["affinity"]}
                self.assertFalse(taken & set(e["affinity"]), "two running jobs pinned to the same CPU")
                active.append(e)
            else:
                active.remove(e)
        self.assertLessEqual(peak, len(cpus))
        for e in logs:
            self.assertTrue(set(e["affinity"]) <= set(cpus))
            expected = 2 if e["kind"] == "cbmlkh" else 1
            self.assertEqual(len(e["affinity"]), expected)
            self.assertEqual(e["omp"], str(expected))

    def test_cbmlkh_wider_than_budget_is_refused(self):
        rc = self.run_cli(self.base_args(cpus=2, extra=["--methods", "cbmlkh", "--repetitions", "4", "--cbmlkh-threads", "4"]))
        self.assertEqual(rc, 2)


class TestCbmlkhExclusion(RunnerTestCase):
    def test_never_two_cbmlkh_at_once_even_across_runners(self):
        self.set_control({n: {"sleep": 0.5} for n in ("t1", "t2", "t3")})
        common = ["--methods", "cbmlkh", "--repetitions", "4", "--cbmlkh-threads", "1"]
        procs = [self.spawn_cli(self.base_args(name=f"r{i}", cpus=4, extra=common)) for i in range(2)]
        for p in procs:
            out, _ = p.communicate(timeout=120)
            self.assertEqual(p.returncode, 0, out)
        logs = [e for e in self.fake_logs() if e["kind"] == "cbmlkh"]
        self.assertEqual(len(logs), 2 * 3 * 4)
        self.assertFalse(overlaps([(e["start"], e["end"]) for e in logs]), "CBMLKH executions overlapped")

    def test_concurrent_runners_never_duplicate_a_job(self):
        self.set_control({n: {"sleep": 0.3} for n in ("t1", "t2", "t3")})
        args = self.base_args(name="shared", cpus=2, extra=["--methods", "ens,ils", "--repetitions", "3", "--allow-concurrent-runners"])
        first = self.spawn_cli(args)
        time.sleep(0.5)
        second = self.spawn_cli(args)
        for p in (first, second):
            out, _ = p.communicate(timeout=120)
            self.assertEqual(p.returncode, 0, out)
        starts = {}
        for e in self.fake_logs():
            starts[(e["kind"], e["instance"], e["seed"])] = starts.get((e["kind"], e["instance"], e["seed"]), 0) + 1
        self.assertEqual(len(starts), 2 * 3 * 3)
        self.assertEqual(set(starts.values()), {1}, "a job was executed twice")


class TestResume(RunnerTestCase):
    def _completed_attempts(self, name="exp"):
        states, _ = self.states(name)
        return {j: s["attempts"] for j, s in states.items() if s["status"] == "completed"}

    def test_crash_and_restart(self):
        self.set_control({n: {"sleep": 0.6, "grandchild": 60} for n in ("t1", "t2", "t3")})
        args = self.base_args(cpus=2, extra=["--methods", "ens,lkh", "--repetitions", "3"])
        runner = self.spawn_cli(args)
        self.assertTrue(self.wait_for(lambda: len(self._completed_attempts()) >= 3 if (self.exp_dir() / "experiment.json").exists() else False, 60))
        runner.kill()  # simulated crash: no cleanup at all
        runner.wait()
        before = self._completed_attempts()
        states, exp = self.states()
        crashed = [j for j, s in states.items() if s["status"] == "running"]
        self.assertTrue(crashed, "no job was running at the crash")
        orphan_pids = [e["grandchild"] for e in self.fake_logs() if "grandchild" in e]

        self.assertEqual(self.run_cli(args), 0)
        states, exp = self.states()
        self.assertEqual({s["status"] for s in states.values()}, {"completed"})
        self.assertEqual(len(states), len(exp["jobs"]))
        for job_id, attempts in before.items():
            self.assertEqual(states[job_id]["attempts"], attempts, f"{job_id} was re-run although it had completed")
        for job_id in crashed:
            events = [e["event"] for e in states[job_id]["events"]]
            self.assertTrue({"stale_running_recovered", "orphan_killed"} & set(events), events)
            self.assertEqual(states[job_id]["attempts"], 2)
        for pid in orphan_pids:
            self.assertTrue(self.wait_for(lambda: not self.pid_alive(pid), 5), f"orphan {pid} survived")
        # Exactly one record per repetition, all valid.
        recs = self.records()
        self.assertEqual(len(recs), 2 * 3 * 3)
        self.assertTrue(all(validate(r) == [] and r["status"] == "completed" for r in recs))

    def test_graceful_interrupt_and_resume(self):
        self.set_control({n: {"sleep": 0.6} for n in ("t1", "t2", "t3")})
        args = self.base_args(cpus=2, extra=["--methods", "ils", "--repetitions", "4"])
        runner = self.spawn_cli(args)
        self.assertTrue(self.wait_for(lambda: (self.exp_dir() / "experiment.json").exists() and len(self._completed_attempts()) >= 2, 60))
        runner.send_signal(signal.SIGINT)
        out, _ = runner.communicate(timeout=60)
        self.assertEqual(runner.returncode, 130, out)
        states, _ = self.states()
        statuses = {s["status"] for s in states.values()}
        self.assertNotIn("running", statuses)
        interrupted = [j for j, s in states.items() if s["status"] == "interrupted"]
        before = self._completed_attempts()
        self.assertEqual(self.run_cli(args), 0)
        states, _ = self.states()
        self.assertEqual({s["status"] for s in states.values()}, {"completed"})
        for job_id, attempts in before.items():
            self.assertEqual(states[job_id]["attempts"], attempts)
        for job_id in interrupted:
            self.assertEqual(states[job_id]["attempts"], 2)

    def test_completed_job_with_lost_record_is_rerun(self):
        args = self.base_args(extra=["--instances", "t1", "--methods", "ens", "--repetitions", "2"])
        self.run_cli(args)
        victim = self.exp_dir() / "runs" / "ens__t1__r01" / f"{run_id('ens', 't1', 1)}.json"
        victim.write_text('{"truncated": ')  # e.g. a damaged file
        self.run_cli(args)
        states, _ = self.states()
        self.assertEqual(states["ens__t1__r01"]["attempts"], 2)
        self.assertEqual(states["ens__t1__r00"]["attempts"], 1)
        self.assertIn("result_invalid", [e["event"] for e in states["ens__t1__r01"]["events"]])

    def test_failed_jobs_are_kept_unless_retry_requested(self):
        self.set_control({"t1": {"mode": "fail"}})
        args = self.base_args(extra=["--instances", "t1", "--methods", "ens", "--repetitions", "1"])
        self.run_cli(args)
        self.run_cli(args)
        states, _ = self.states()
        self.assertEqual(states["ens__t1__r00"]["attempts"], 1)
        self.set_control({})
        self.run_cli(args + ["--retry-failed"])
        states, _ = self.states()
        self.assertEqual(states["ens__t1__r00"]["status"], "completed")
        self.assertEqual(states["ens__t1__r00"]["attempts"], 2)
        superseded = list((self.exp_dir() / "runs" / "ens__t1__r00").glob("attempt-01/superseded/*.json"))
        self.assertTrue(superseded, "the failed attempt's record was not preserved")


class TestDesignProtection(RunnerTestCase):
    def test_partial_design_cannot_be_named_final(self):
        rc = self.run_cli(self.base_args(name="final", extra=["--methods", "ens"]))
        self.assertEqual(rc, 2)
        self.assertFalse((self.exp_dir("final") / "experiment.json").exists())

    def test_partial_default_name_and_changed_design_refused(self):
        args = self.base_args(extra=["--instances", "t1", "--methods", "ens", "--repetitions", "1"])
        self.run_cli(args)
        exp = json.loads((self.exp_dir() / "experiment.json").read_text())
        self.assertTrue(exp["partial"])
        changed = self.base_args(extra=["--instances", "t1", "--methods", "ens", "--repetitions", "2"])
        self.assertEqual(self.run_cli(changed), 2)
        resume_without_design = [a for a in self.base_args() if a not in ("--instance-dir", str(self.instances))]
        self.assertEqual(self.run_cli(resume_without_design), 0)


class TestReports(RunnerTestCase):
    def test_csvs_match_records(self):
        self.set_control({"t3": {"mode": "fail"}})
        self.run_cli(self.base_args(cpus=3, extra=["--methods", "ens,lkh,cbmlkh", "--repetitions", "3", "--cbmlkh-threads", "3"]))
        recs = self.records()
        detailed = self.read_csv("exp", "detailed")
        self.assertEqual(len(detailed), 3 * 3 * 3)
        by_id = {r["run_id"]: r for r in recs}
        for row in detailed:
            rec = by_id[row["run_id"]]
            self.assertEqual(row["status"], rec["status"])
            self.assertEqual(row["objective"], "" if rec["objective"]["value"] is None else str(rec["objective"]["value"]))
        summary = self.read_csv("exp", "summary")
        for row in summary:
            values = [r["objective"]["value"] for r in recs if r["instance"]["name"] == row["instance"] and r["method"]["id"] == row["method"] and r["status"] == "completed"]
            self.assertEqual(int(row["n_used"]), len(values))
            if values:
                self.assertAlmostEqual(float(row["objective_mean"]), statistics.fmean(values), places=4)
                self.assertEqual(float(row["objective_min"]), min(values))
            else:
                self.assertEqual(row["objective_mean"], "")
        failed = [r for r in summary if r["instance"] == "t3" and r["method"] in ("ens", "lkh")]
        self.assertTrue(all(r["n_failed"] == "3" for r in failed))
        # Regenerating from disk gives identical files.
        before = {w: (self.exp_dir() / "results" / f"{w}.csv").read_text() for w in ("detailed", "summary", "comparison", "methods")}
        self.run_cli(["report", "--outputs", str(self.outputs), "--name", "exp"])
        after = {w: (self.exp_dir() / "results" / f"{w}.csv").read_text() for w in before}
        self.assertEqual(before, after)
        manifest = json.loads((self.exp_dir() / "manifest.json").read_text())
        self.assertEqual(len(manifest["jobs"]), 3 * 2 * 3 + 3)


if __name__ == "__main__":
    unittest.main()
