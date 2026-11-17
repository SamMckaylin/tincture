"""Docs stay in sync with the code."""

from __future__ import annotations

import runpy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_api_reference_is_up_to_date():
    mod = runpy.run_path(str(ROOT / "docs" / "make_api.py"), run_name="make_api")
    assert mod["main"](["--check"]) == 0, "docs/api.md is stale: run python docs/make_api.py"
