"""Resumable, append-only persistence for benchmark results.

Results live under ``<base_dir>/<run_tag>/``, one JSON file per instance plus an
``index.csv`` summary. An existing per-instance file marks that instance as done
so a re-run resumes where it left off; :meth:`save` never overwrites an existing
result.
"""

from __future__ import annotations

import csv
import json
import threading
from pathlib import Path

from .config import INDEX_FIELDS


class ResultStore:
    def __init__(self, base_dir: str | Path, run_tag: str) -> None:
        self.run_tag = run_tag
        self.dir = Path(base_dir) / run_tag
        self.dir.mkdir(parents=True, exist_ok=True)
        self.index_path = self.dir / "index.csv"
        # Guards the shared index.csv append when instances are saved concurrently.
        self._lock = threading.Lock()

    def result_path(self, instance_name: str) -> Path:
        return self.dir / f"{instance_name}.json"

    def is_done(self, instance_name: str) -> bool:
        return self.result_path(instance_name).exists()

    def save(self, record: dict) -> Path:
        """Persist a result atomically. Refuses to overwrite an existing one.

        Thread-safe: per-instance JSON files are written to distinct paths, and
        the shared index.csv append is serialized with a lock.
        """
        path = self.result_path(record["instance"])
        if path.exists():
            raise FileExistsError(f"result already exists, refusing to overwrite: {path}")

        tmp = path.parent / f"{path.name}.{threading.get_ident()}.tmp"
        tmp.write_text(json.dumps(record, indent=2))
        tmp.replace(path)  # atomic on POSIX
        with self._lock:
            self._append_index(record)
        return path

    def load_all(self) -> list[dict]:
        return [json.loads(p.read_text()) for p in sorted(self.dir.glob("*.json"))]

    def _append_index(self, record: dict) -> None:
        write_header = not self.index_path.exists()
        with self.index_path.open("a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=INDEX_FIELDS, extrasaction="ignore")
            if write_header:
                writer.writeheader()
            writer.writerow(record)
