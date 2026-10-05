"""CPU budget resolution: one core is left free by default."""

import argparse
import unittest
from unittest import mock

import helpers  # noqa: F401  (sets sys.path)

from cbm_experiments import cli
from cbm_experiments.runner import ExperimentError

DETECTED = [0, 2, 4, 6, 8, 10, 12, 14, 16, 17, 18, 19, 20, 21, 22, 23]  # one CPU per physical core


def resolve(**kw):
    args = argparse.Namespace(**{"cpus": None, "reserve_cores": 1, "max_cpus": None, **kw})
    with mock.patch.object(cli, "physical_core_cpus", return_value=list(DETECTED)):
        return cli.resolve_cpus(args)


class TestCpuBudget(unittest.TestCase):
    def test_one_core_reserved_by_default(self):
        self.assertEqual(resolve(), DETECTED[1:])

    def test_reserve_zero_and_more(self):
        self.assertEqual(resolve(reserve_cores=0), DETECTED)
        self.assertEqual(resolve(reserve_cores=3), DETECTED[3:])

    def test_explicit_list_is_exact(self):
        self.assertEqual(resolve(cpus="0-3"), [0, 1, 2, 3])

    def test_max_cpus_after_reservation(self):
        self.assertEqual(resolve(max_cpus=4), DETECTED[1:5])

    def test_nothing_left_is_an_error(self):
        with self.assertRaises(ExperimentError):
            resolve(reserve_cores=16)

    def test_cbmlkh_threads_follow_the_budget(self):
        # 15 usable CPUs here -> still one 10-thread execution per instance.
        self.assertEqual(min(10, len(resolve())), 10)


if __name__ == "__main__":
    unittest.main()
