"""Shared run-record schema (``experiments/schema/run_record.schema.json``).

The schema file is the reference; :func:`validate` implements the subset of
JSON Schema it uses (type, required, properties, enum, const, items, minimum)
so that validation needs no third-party package.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, List

SCHEMA_VERSION = "1.0"
SCHEMA_PATH = Path(__file__).resolve().parent.parent / "schema" / "run_record.schema.json"

# Terminal statuses of a repetition record. "pending" and "running" exist only
# in job state, never in a record.
RECORD_STATUSES = ("completed", "failed", "hard_timeout", "interrupted")

_TYPES = {
    "object": dict,
    "array": list,
    "string": str,
    "boolean": bool,
    "null": type(None),
}


def load_schema() -> dict:
    with open(SCHEMA_PATH) as f:
        return json.load(f)


def _type_ok(value: Any, name: str) -> bool:
    if name == "integer":
        return isinstance(value, int) and not isinstance(value, bool)
    if name == "number":
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    return isinstance(value, _TYPES[name])


def _check(value: Any, schema: dict, where: str, errors: List[str]) -> None:
    if "const" in schema and value != schema["const"]:
        errors.append(f"{where}: expected {schema['const']!r}, got {value!r}")
    if "enum" in schema and value not in schema["enum"]:
        errors.append(f"{where}: {value!r} not in {schema['enum']}")
    if "type" in schema:
        types = schema["type"] if isinstance(schema["type"], list) else [schema["type"]]
        if not any(_type_ok(value, t) for t in types):
            errors.append(f"{where}: expected type {types}, got {type(value).__name__}")
            return
    if "minimum" in schema and isinstance(value, (int, float)) and value < schema["minimum"]:
        errors.append(f"{where}: {value} < {schema['minimum']}")
    if isinstance(value, dict):
        for key in schema.get("required", []):
            if key not in value:
                errors.append(f"{where}: missing required field '{key}'")
        for key, sub in schema.get("properties", {}).items():
            if key in value:
                _check(value[key], sub, f"{where}.{key}", errors)
    if isinstance(value, list) and "items" in schema:
        for i, item in enumerate(value):
            _check(item, schema["items"], f"{where}[{i}]", errors)


def validate(record: Any, schema: dict = None) -> List[str]:
    """List of schema violations (empty if the record conforms)."""
    errors: List[str] = []
    _check(record, schema or load_schema(), "record", errors)
    return errors
