"""Entry point: run pure LKH and/or CBMLKH on the same instances and compare.

There are no command-line arguments — all parameters live in ``config.py``
(:class:`BenchmarkConfig` and the constants around it). Edit them there, then run
the module. Designed for long unattended runs, e.g.::

    cd src/CBMLKH
    nohup python -m lkh_benchmark > lkh_run.log 2>&1 &

Re-running resumes: finished instances are skipped per solver and never
overwritten.
"""

from __future__ import annotations

from .config import BenchmarkConfig
from .pipeline import run_comparison


def main() -> int:
    run_comparison(BenchmarkConfig())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
