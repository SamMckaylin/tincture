from __future__ import annotations

import math

import numpy as np
import pytest

from tincture.errors import TinctureSpecError
from tincture.report import TEMPLATES, Caption, Methods, fmt_conc, fmt_p
from tincture.report.captions import template_fields
from tincture.report.classic import build_classic, registered


def test_caption_fills_template():
    c = Caption().add(
        TEMPLATES["compass"], compound="X", top_ref="HDACi", w=0.71, pct_unexplained=0.23, n_refs=12
    )
    assert c.render() == (
        "X most resembles HDACi (weight 0.71); 23% of its signature is unexplained by the "
        "12 references."
    )


@pytest.mark.parametrize("missing", [None, math.nan, np.float32("nan"), np.nan])
def test_sentence_dropped_when_stat_unavailable(missing):
    c = Caption().add("A is {a}.", a=1).add("B p {p}.", p=missing)
    assert c.render() == "A is 1."
    assert c.dropped == ["B p {p}."]


def test_template_field_not_supplied_is_a_bug():
    with pytest.raises(TinctureSpecError, match="needs fields"):
        Caption().add("{a} and {b}", a=1)


def test_weak_notes_always_shown_once():
    c = Caption().add("Main {x}.", x=1).weak("fit failed").weak("fit failed")
    assert c.render() == "Main 1. Note: fit failed"


def test_add_if():
    assert Caption().add_if(False, "x {a}", a=1).render() == ""


def test_all_brief_templates_parse():
    assert set(template_fields(TEMPLATES["dose_paths"])) == {
        "compound",
        "ref",
        "dose_ok",
        "dose_bend",
        "stress_axis",
        "p",
    }


@pytest.mark.parametrize("p,out", [(0.0004, "< 0.001"), (0.0312, "= 0.031"), (None, "n/a")])
def test_fmt_p(p, out):
    assert fmt_p(p) == out


@pytest.mark.parametrize(
    "m,out",
    [
        (3.2e-7, "320 nM"),
        (1e-5, "10 µM"),
        (2e-3, "2 mM"),
        (5e-12, "5 pM"),
        (0, "n/a"),
        (999.9e-9, "1 µM"),
        (1e-6, "1 µM"),
        (2.5e-16, "0.25 fM"),
    ],
)
def test_fmt_conc(m, out):
    assert fmt_conc(m) == out


def test_methods_blocks_and_citations():
    m = Methods().use("bootstrap_ci", level=0.95, n_boot=1000, unit="replicates")
    m.use("bh_fdr", alpha=0.05).use("bh_fdr", alpha=0.05)  # second use ignored
    text = m.render()
    assert "95% percentile bootstrap" in text and "Benjamini" in text
    assert text.count("Benjamini–Hochberg procedure") == 1
    assert "Software: tincture" in text and "References: [1]" in text


def test_methods_unknown_or_unfilled_block():
    with pytest.raises(TinctureSpecError, match="Unknown"):
        Methods().use("nope")
    with pytest.raises(TinctureSpecError, match="needs fields"):
        Methods().use("bh_fdr")


def test_classic_registry():
    import tincture.plots.demo  # noqa: F401  (registers the demo fallback)

    assert registered()["demo"] == "scatter of raw values"
    assert build_classic("no-such-plot") is None
