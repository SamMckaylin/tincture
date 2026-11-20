"""User-study kit: answer matching and item definitions (no data needed)."""

from __future__ import annotations

import runpy
from pathlib import Path

US = runpy.run_path(str(Path(__file__).resolve().parents[1] / "case_studies" / "user_study.py"))

FACTS = {  # shape of CS1 results.json (values are synthetic, for the test only)
    "retrievable_by_dose": {"10 µM": True},
    "compass_top_reference": "HDAC",
    "compass_unexplained": 0.03,
    "abexinostat_bends": False,
    "panobinostat_bend_above": "1 µM",
    "fingerprint_top_pathways": [{"pathway": "E2F targets", "direction": "down"}],
    "fingerprint_nearest": "X1",
    "fingerprint_nearest_name": "M344",
    "ladder_counts": {"up": 10, "down": 20, "total": 100},
    "ladder_top_gene": "AKR1B10",
}


def test_answer_matching():
    its = {i.id: i for i in US["items"](FACTS)}
    ok = US["is_correct"]
    assert ok(its["Q02"], "HDAC") and ok(its["Q02"], "histone deacetylase inhibitors")
    assert not ok(its["Q02"], "Aurora kinase")
    assert ok(its["Q05"], "yes, above 1 µM") and ok(its["Q05"], "above 1uM")
    assert not ok(its["Q05"], "above 10 µM") and not ok(its["Q05"], "no")
    assert ok(its["Q09"], "more down") and not ok(its["Q09"], "increased")
    assert ok(its["Q08"], "akr1b10") and not ok(its["Q08"], "ABCC2")


def test_items_flag_answerability():
    its = US["items"](FACTS)
    assert len({i.id for i in its}) == len(its)
    primary = [i.id for i in its if i.both_answerable]
    assert {"Q01", "Q02", "Q04", "Q05", "Q08", "Q09"} == set(primary)
