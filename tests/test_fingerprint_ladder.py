"""Fingerprint Strip and Effect Ladder: statistics, guards, encoding, budget, baselines."""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
import pytest

import tincture as tc
from tincture.errors import TinctureDataError
from tincture.model.signature import Signature
from tincture.plots._pathways import curated_order, group_breaks, pathway_label
from tincture.spec import HeatmapMark, IntervalMark, NoticeMark, RectMark
from tincture.theme import check_colorblind

P = 300
FEATS = np.array([f"g{i}" for i in range(P)])


def _sets() -> pd.DataFrame:
    rows = []
    blocks = {
        "HALLMARK_E2F_TARGETS": range(0, 30),
        "HALLMARK_P53_PATHWAY": range(30, 60),
        "HALLMARK_GLYCOLYSIS": range(60, 90),
        "HALLMARK_INTERFERON_GAMMA_RESPONSE": range(120, 150),
        "SYNTHETIC_E2F_COPY": range(0, 26),  # Jaccard 26/30 with E2F → dropped (smaller)
    }
    for s, idx in blocks.items():
        rows += [(s, f"g{i}") for i in idx]
    return pd.DataFrame(rows, columns=["source", "target"])


def _sig(cid: str, eff: np.ndarray, se: float = 0.1, dose: float = 1e-6) -> Signature:
    return Signature(
        features=FEATS,
        effect=eff,
        se=np.full(P, se),
        compound_id=cid,
        dose=dose,
        method="logcpm_eb+shrink",
        meta={"compound_name": cid.replace("synthetic_", "")},
    )


def _profiles(seed: int = 0) -> list[Signature]:
    rng = np.random.default_rng(seed)
    out = []
    for i in range(4):  # E2F down (HDAC-like)
        e = rng.normal(0, 0.1, P)
        e[0:30] -= 1.0
        out.append(_sig(f"synthetic_e2f{i}", e))
    for i in range(3):  # p53 up
        e = rng.normal(0, 0.1, P)
        e[30:60] += 0.8
        out.append(_sig(f"synthetic_p53{i}", e))
    e = rng.normal(0, 0.1, P)
    e[60:90] += 0.6  # glycolysis up
    out.append(_sig("synthetic_gly", e))
    return out


# --------------------------------------------------------------------------- pathways
def test_curated_order_and_groups():
    names = ["GLYCOLYSIS", "HALLMARK_P53_PATHWAY", "ZZZ_CUSTOM", "E2F_TARGETS", "COMPLEMENT"]
    order = curated_order(names)
    assert order == ["E2F_TARGETS", "HALLMARK_P53_PATHWAY", "GLYCOLYSIS", "COMPLEMENT",
                     "ZZZ_CUSTOM"]  # fmt: skip
    assert group_breaks(order) == [0, 1, 2, 3]  # Cell cycle | Stress | Metabolism | Immune | Other
    assert pathway_label("HALLMARK_TNFA_SIGNALING_VIA_NFKB") == "TNFα via NF-κB"


# --------------------------------------------------------------------------- fingerprint
def test_fingerprint_encoding_and_guards():
    sigs = _profiles()
    res = tc.fingerprint(sigs, _sets(), focus="synthetic_e2f3",
                         references=["synthetic_e2f0", "synthetic_p530"])  # fmt: skip
    m = res.meta
    assert m["rows"][:2] == ["synthetic_e2f0", "synthetic_p530"]  # references first, as given
    hm = next(iter(res.spec.marks(HeatmapMark)))
    assert hm.row_gap_after == [1]  # gap after the references
    assert m["columns"][:3] == ["E2F_TARGETS", "P53_PATHWAY", "GLYCOLYSIS"]  # curated order
    assert hm.col_gap_after == list(range(len(m["columns"]) - 1))  # one group per column here
    # non-significant cells are neutral (masked)
    q = m["q"].to_numpy()
    assert np.array_equal(np.asarray(hm.mask), ~(q < 0.05))
    assert m["q"].loc["synthetic_p531", "E2F_TARGETS"] >= 0.05
    # clip shown in legend and footer; ordering stated
    assert res.spec.panels[0].colorbar.clip_note == f"clipped at ±{m['clip']:g}"
    assert "Rows: references first in the given order, then average-linkage" in res.caption
    # drop_empty hides columns without any significant cell, and says how many
    full = tc.fingerprint(sigs, _sets(), drop_empty=False)
    shown = [c for c in full.meta["columns"] if (full.meta["q"][c] < 0.05).any()]
    assert m["columns"] == shown
    hidden = len(full.meta["columns"]) - len(shown)
    assert (f"{hidden} pathway(s) with no significant cell" in res.caption) == (hidden > 0)
    # dedup: the near-copy of E2F was dropped (Jaccard > 0.5, keep the larger)
    assert "1 near-duplicate sets" in res.methods
    # focus tick and nearest row
    tick = [r for r in res.spec.marks(RectMark) if r.color == "focus"]
    assert tick and not tick[0].clip
    assert m["nearest"].startswith("synthetic_e2f")
    assert "strongest effect: E2F targets (down)" in res.caption
    assert res.classic is not None and res.classic.plot == "volcano"


def test_fingerprint_given_order_gate_and_errors():
    sigs = _profiles()
    ids = [s.compound_id for s in sigs]
    res = tc.fingerprint(sigs, _sets(), order="given", gate=set(ids[:-1]))
    assert res.meta["rows"] == ids[:-1]
    assert "Rows: the given order" in res.caption
    assert "failed the retrievability gate" in res.caption
    with pytest.raises(TinctureDataError, match="Focus"):
        tc.fingerprint(sigs, _sets(), focus="nope")
    flat = [_sig(f"synthetic_n{i}", np.random.default_rng(i).normal(0, 0.1, P)) for i in range(3)]
    with pytest.raises(TinctureDataError, match="No pathway is significant"):
        tc.fingerprint(flat, _sets())


def test_fingerprint_budget_1000_by_50():
    rng = np.random.default_rng(7)
    sets = pd.DataFrame(
        [(f"SYNTHETIC_SET_{k:02d}", f"g{i}") for k in range(50) for i in range(k * 6, k * 6 + 6)],
        columns=["source", "target"],
    )
    base = rng.normal(0, 1, (5, P))
    sigs = [_sig(f"synthetic_c{i}", base[i % 5] + rng.normal(0, 0.3, P)) for i in range(1000)]
    t0 = time.perf_counter()
    res = tc.fingerprint(sigs, sets, drop_empty=False)
    res.figure()
    elapsed = time.perf_counter() - t0
    assert res.meta["score"].shape == (1000, 50)
    assert elapsed < 5.0, f"1,000 × 50 took {elapsed:.1f} s (budget 5 s)"
    assert "Only reference and focus rows are labelled" in res.caption


@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_fingerprint_baseline():
    fig = tc.fingerprint(_profiles(), _sets(), focus="synthetic_e2f3",
                         references=["synthetic_e2f0", "synthetic_p530"]).figure()  # fmt: skip
    check_colorblind(fig)
    return fig


# --------------------------------------------------------------------------- effect ladder
def _ladder_sig() -> Signature:
    rng = np.random.default_rng(3)
    e = rng.normal(0, 0.05, P)
    e[0:10] = -np.linspace(1.0, 2.0, 10)  # E2F members down
    e[30:35] = np.linspace(1.0, 1.5, 5)  # p53 members up
    e[200:205] = 1.2  # outside every set
    se = np.full(P, 0.1)
    se[200] = 1.0  # large SE: interval crosses zero
    return Signature(
        features=FEATS, effect=e, se=se, compound_id="synthetic_lad", dose=1e-6,
        method="logcpm_eb+shrink", meta={"compound_name": "lad"},
    )  # fmt: skip


def test_effect_ladder_count_strip_ci_and_printed_n():
    res = tc.effect_ladder(_ladder_sig(), gene_sets=_sets(), top_n=20)
    t, c = res.meta["table"], res.meta["counts"]
    assert len(t) == 20 and c["total"] == P
    assert c["down"] >= 10 and c["up"] >= 9
    assert res.spec.panels[0].key == "counts"  # count strip is always drawn, on top
    assert any(m.text.startswith(f"{c['up']:,} up") for m in res.spec.panels[0].marks
               if hasattr(m, "text"))  # fmt: skip
    assert res.spec.footer[0].startswith("N = 20 features; shrinkage: limma-style")
    # grouped in curated order; the unassigned group last
    groups = list(dict.fromkeys(t["gene_set"]))
    assert groups == ["E2F_TARGETS", "P53_PATHWAY", ""]
    # CI colouring: crossing zero → grey
    iv = next(iter(res.spec.marks(IntervalMark)))
    i200 = int(np.flatnonzero(t["feature"] == "g200")[0])
    assert iv.color[i200] == "mid_grey" and t["call"][i200] == "ns"
    assert "intervals crossing zero" in res.caption
    assert res.classic is not None and res.classic.plot == "volcano"


def test_effect_ladder_drill_down_and_errors():
    res = tc.effect_ladder(_ladder_sig(), gene_sets=_sets(), pathway="P53_PATHWAY", top_n=8)
    assert set(res.meta["table"]["symbol"]) <= {f"g{i}" for i in range(30, 60)}
    assert "in p53 pathway" in res.caption
    with pytest.raises(TinctureDataError, match="not found"):
        tc.effect_ladder(_ladder_sig(), gene_sets=_sets(), pathway="NOPE")
    nose = Signature(features=FEATS, effect=np.zeros(P), compound_id="synthetic_x")
    with pytest.raises(TinctureDataError, match="standard errors"):
        tc.effect_ladder(nose)


def test_effect_ladder_no_significant_features_says_so():
    s = Signature(features=FEATS, effect=np.random.default_rng(0).normal(0, 0.01, P),
                  se=np.full(P, 1.0), compound_id="synthetic_flat")  # fmt: skip
    res = tc.effect_ladder(s)
    assert res.meta["counts"]["up"] + res.meta["counts"]["down"] == 0
    assert any(isinstance(m, NoticeMark) for m in res.spec.panels[0].marks)
    assert "shrinkage: as supplied (method not recorded)" in res.spec.footer[0]


@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_effect_ladder_baseline():
    fig = tc.effect_ladder(_ladder_sig(), gene_sets=_sets(), top_n=20).figure()
    check_colorblind(fig)
    return fig
