"""Retrievability: mAP statistics (incl. copairs golden), nulls, guards, plot, gate."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import tincture as tc
from tests.fixtures.synthetic_profiles import synthetic_profiles
from tincture.errors import TinctureDataError
from tincture.plots.retrievability import GATE_KEY, compute, gate_mask
from tincture.spec import NoticeMark, PointMark, RefLineMark
from tincture.stats.map_retrieval import (
    ap_from_relevance,
    average_precision,
    cat_curve,
    min_achievable_p,
    random_ap,
    retrievability,
)
from tincture.theme import check_colorblind

GOLDEN = Path(__file__).parent / "fixtures" / "copairs_golden.json"


# --------------------------------------------------------------------------- AP basics
def test_ap_hand_computed():
    # relevance 1,0,1,0 -> precisions at hits 1/1 and 2/3 -> AP = (1 + 2/3) / 2
    assert ap_from_relevance(np.array([[1, 0, 1, 0]]))[0] == pytest.approx(5 / 6)
    assert ap_from_relevance(np.array([[0, 0, 0, 1]]))[0] == pytest.approx(1 / 4)


def test_random_ap_mean_matches_theory():
    rng = np.random.default_rng(0)
    # one positive among n: E[AP] = H_n / n
    n = 10
    null = random_ap(50_000, 1, n, rng)
    assert null.mean() == pytest.approx(sum(1 / k for k in range(1, n + 1)) / n, rel=0.02)


# --------------------------------------------------------------------------- copairs parity
def test_matches_copairs_golden():
    g = json.loads(GOLDEN.read_text())
    X, lab = np.array(g["X"]), np.array(g["labels"])
    r = retrievability(X, lab, lab == "DMSO", null="copied", null_size=g["null_size"], seed=1)
    ours = r.ap.set_index("profile")
    gold = pd.DataFrame(g["ap"]).set_index("profile")
    assert np.allclose(ours.loc[gold.index, "ap"], gold["average_precision"], atol=1e-9)
    assert (ours.loc[gold.index, "n_pos"] == gold["n_pos_pairs"]).all()
    assert (ours.loc[gold.index, "n_total"] == gold["n_total_pairs"]).all()
    t = r.table.set_index("group")
    gm = pd.DataFrame(g["map"]).set_index("compound")
    assert np.allclose(t.loc[gm.index, "map"], gm["mean_average_precision"], atol=1e-9)
    # p-values: independent Monte Carlo nulls, so compare within MC error
    assert np.allclose(t.loc[gm.index, "p"], gm["p_value"], atol=0.01)
    assert (t.loc[gm.index, "q"] < 0.05).tolist() == (gm["corrected_p_value"] < 0.05).tolist()


# --------------------------------------------------------------------------- permutation null
def test_permutation_null_detects_signal_and_is_exact():
    ps = synthetic_profiles()
    r = compute(ps)
    t = r.table.set_index("group")
    assert t.loc["synthetic_strong|1e-06", "retrievable"]
    assert not t.loc["synthetic_none|1e-06", "retrievable"]
    assert t["exact"].all()  # C(19, 3) = 969 <= 10,000 relabelings enumerated


def test_permutation_null_is_calibrated():
    ps = synthetic_profiles(
        effects={f"synthetic_null{i}": 0.0 for i in range(40)}, n_rep=2, single=False, seed=3
    )
    p = compute(ps).table["p"]
    assert 0.35 < p.mean() < 0.65
    assert (p < 0.05).mean() < 0.15


def test_monte_carlo_when_too_many_relabelings():
    ps = synthetic_profiles(n_controls=40, n_rep=4)
    t = compute(ps, null_size=1000).table
    assert not t["exact"].any() and t["p"].min() >= 1 / 1001


def test_min_achievable_p():
    assert min_achievable_p(2, 32, "copied", 10_000) == pytest.approx(1 / 33)
    assert min_achievable_p(2, 32, "permutation", 10_000) == pytest.approx(1 / 561)


# --------------------------------------------------------------------------- guards
def test_design_that_cannot_pass_says_so():
    """Guard: duplicates + copied null over many conditions cannot reach q < 0.05."""
    ps = synthetic_profiles(
        effects={f"synthetic_c{i}": 4.0 for i in range(60)}, n_rep=2, n_controls=10, single=False
    )
    res = tc.retrievability(ps, null="copied", null_size=2000)
    assert res.meta["design_can_pass"] is False
    assert "cannot reach q < 0.05" in res.caption
    assert any("cannot reach" in m.text for m in res.spec.marks(NoticeMark))


def test_single_replicate_excluded_and_reported():
    res = tc.retrievability(synthetic_profiles())
    assert res.meta["excluded"] == ["synthetic_single|1e-06"]
    assert "single replicate could not be scored" in res.caption
    assert any("unscored" in m.text for m in res.spec.marks(NoticeMark))


def test_non_retrievable_greyed_retrievable_blue_focus_terracotta():
    res = tc.retrievability(synthetic_profiles(), focus=["synthetic_medium|1e-06"])
    pts = res.spec.marks(PointMark)
    by_legend = {p.legend: p for p in pts}
    assert by_legend["not retrievable"].color == "ink_soft" and by_legend["not retrievable"].hollow
    assert (
        by_legend[next(k for k in by_legend if k and k.startswith("retrievable"))].color
        == "reference"
    )
    assert by_legend["focus"].color == "focus"
    from tincture.spec import TextMark

    texts = {m.text: m for m in res.spec.panel("strip").marks if isinstance(m, TextMark)}
    assert texts["synthetic_medium 1 µM"].color == "focus"  # direct label, no legend box
    assert "retrievable (q < 0.05)" in texts and "not retrievable" in texts
    assert res.spec.panel("strip").legend is None


def test_q_line_drawn_when_separable():
    res = tc.retrievability(synthetic_profiles())
    lines = [m for m in res.spec.marks(RefLineMark) if m.label and m.label.startswith("q =")]
    assert len(lines) == 1 and res.meta["q_boundary"] is not None


def test_footer_and_methods_name_the_null():
    res = tc.retrievability(synthetic_profiles())
    assert any("permutation null" in f for f in res.spec.footer)
    assert "Benjamini" in res.methods and "relabelling" in res.methods


def test_no_controls_raises():
    ps = synthetic_profiles()
    with pytest.raises(TinctureDataError, match="control"):
        ps2 = ps.subset(obs=(ps.obs["dose"] > 0).to_numpy())
        compute(ps2)


def test_no_replicates_raises():
    ps = synthetic_profiles(n_rep=1, single=False)
    with pytest.raises(TinctureDataError, match="replicate"):
        compute(ps)


def test_nan_and_zero_profiles_raise():
    X = np.ones((4, 3))
    X[0, 0] = np.nan
    with pytest.raises(TinctureDataError, match="NaN"):
        average_precision(X, ["a", "a", "c", "c"], np.array([0, 0, 1, 1], bool))
    with pytest.raises(TinctureDataError, match="all-zero"):
        average_precision(np.zeros((4, 3)), ["a", "a", "c", "c"], np.array([0, 0, 1, 1], bool))


def test_permutation_requires_control_negatives():
    ps = synthetic_profiles()
    with pytest.raises(ValueError, match="permutation"):
        compute(ps, negatives="all")
    assert len(compute(ps, negatives="all", null="copied").table) == 3


def test_unknown_focus_raises():
    with pytest.raises(TinctureDataError, match="no retrievability score"):
        tc.retrievability(synthetic_profiles(), focus=["nope|1"])


# --------------------------------------------------------------------------- gate + classic
def test_gate_mask_and_uns():
    ps = synthetic_profiles()
    r = compute(ps)
    mask = gate_mask(ps, r)
    assert mask[ps.is_vehicle].all()
    assert mask[(ps.obs["compound_id"] == "synthetic_strong").to_numpy()].all()
    assert not mask[(ps.obs["compound_id"] == "synthetic_none").to_numpy()].any()
    assert "synthetic_strong|1e-06" in ps.adata.uns[GATE_KEY]["retrievable_groups"]


def test_classic_fallback_is_replicate_scatter():
    res = tc.retrievability(synthetic_profiles(), focus=["synthetic_strong|1e-06"])
    assert res.classic is not None and "Pearson r" in res.classic.caption
    assert res.classic.meta["classic_of"] == "retrievability"


# --------------------------------------------------------------------------- CAT
def test_cat_identical_and_random():
    rng = np.random.default_rng(0)
    a = rng.normal(size=500)
    same = cat_curve(a, a)
    assert np.allclose(same.observed, 1.0)
    rand = cat_curve(a, rng.normal(size=500))
    inside = (rand.observed >= rand.null_lo - 1e-12) & (rand.observed <= rand.null_hi + 1e-12)
    assert inside.mean() > 0.8
    assert np.allclose(rand.null_mean, rand.k / 500)
    with pytest.raises(TinctureDataError):
        cat_curve(a, a[:10])


# --------------------------------------------------------------------------- visual
@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_retrievability_baseline():
    fig = tc.retrievability(synthetic_profiles(), focus=["synthetic_strong|1e-06"]).figure()
    check_colorblind(fig)
    return fig
