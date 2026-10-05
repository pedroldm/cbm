#!/usr/bin/env python3
"""Entry point: python3 experiments/run_experiments.py {run,plan,status,report} --help"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cbm_experiments.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
