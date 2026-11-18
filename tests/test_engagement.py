"""Engagement Ladder: shifts, the three flags, ATP wedge, cell-system and missing-rung
guards (HANDOVER §5.10), baseline."""

from __future__ import annotations

import numpy as np
import pytest

import tincture as tc
from tests.fixtures.synthetic_ladder import ladder, typical
from tincture.errors import TinctureDataError
from tincture.plots.engagement import atp_shift
from tincture.spec import LineMark, PointMark, PolygonMark, TextMark
from tincture.theme import check_colorblind

PAIR = [("synthetic_a", "SYNK1")]


def _flags(res):
    return res.meta["results"][0].flags


def test_shift_and_no_flags_when_consistent():
    res = tc.engagement_ladder(typical(), PAIR)
    r = res.meta["results"][0]
    assert r.shift_bc == pytest.approx(1.5) and not r.flags
    assert "1.5 log shift biochemical→cellular" in res.caption
    assert "phenotype is 3× weaker than cellular target engagement" in res.caption


def test_off_target_flag_when_phenotype_10x_more_potent():
    res = tc.engagement_ladder(typical(pheno_px=8.7, cell_px=7.5), PAIR)  # 16×
    assert "off_target" in _flags(res)
    assert "likely driven by something other than SYNK1" in res.caption
    assert any(
        isinstance(m, TextMark) and m.text.startswith("Flag:") and m.color == "focus"
        for m in res.spec.marks()
    )
    near = tc.engagement_ladder(typical(pheno_px=8.3, cell_px=7.5), PAIR)  # 6×: no flag
    assert "off_target" not in _flags(near)


def test_censored_values_only_used_when_conclusive():
    # phenotype more potent than the lowest dose: still conclusive
    res = tc.engagement_ladder(typical(pheno_px=9.0, cell_px=7.5, pheno_cens="<"), PAIR)
    assert "≥" in _flags(res)["off_target"]
    hollow = [m for m in res.spec.marks(PointMark) if m.hollow is True]
    assert hollow  # censored drawn hollow
    # phenotype only weaker-than bound: inconclusive, says so
    res2 = tc.engagement_ladder(typical(pheno_px=9.0, cell_px=7.5, pheno_cens=">"), PAIR)
    r = res2.meta["results"][0]
    assert "off_target" not in r.flags and "off_target" in r.not_evaluated


def test_atp_wedge_and_shift_flag():
    lo, hi = atp_shift(1e-5)  # Km 10 µM, ATP 1–5 mM
    assert lo == pytest.approx(np.log10(101)) and hi == pytest.approx(np.log10(501))
    within = tc.engagement_ladder(typical(cell_px=7.0), PAIR, km_atp={"SYNK1": 1e-5})  # 2.0 log
    assert "shift" not in _flags(within)
    assert any(isinstance(m, PolygonMark) for m in within.spec.marks())  # grey wedge
    beyond = tc.engagement_ladder(typical(cell_px=5.5), PAIR, km_atp={"SYNK1": 1e-5})  # 3.5 log
    assert "permeability, efflux or binding issue" in _flags(beyond)["shift"]
    no_km = tc.engagement_ladder(typical(cell_px=5.5), PAIR)
    r = no_km.meta["results"][0]
    assert "shift" not in r.flags and "shift" in r.not_evaluated
    assert "no ATP Km" in no_km.caption
    assert not any(isinstance(m, PolygonMark) for m in no_km.spec.marks())


def test_exposure_flag():
    t = typical()
    t.data.loc[len(t.data)] = {
        **t.data.iloc[0].to_dict(),
        "rung": "exposure",
        "px": 7.6,
        "assay_type": "IC50",
        "cell_system": "plasma",
    }
    res = tc.engagement_ladder(t, PAIR)
    assert "unlikely to reach active concentration in vivo" in _flags(res)["exposure"]


def test_dashed_between_cell_systems_and_missing_rungs_skipped():
    res = tc.engagement_ladder(typical(cell_sys="SYN-LINE", pheno_sys="SYN-LINE"), PAIR)
    styles = [m.style for m in res.spec.marks(LineMark)]
    assert styles == ["dashed", "solid"]  # cell-free→cell dashed, same line solid
    assert "Dashed segments join" in res.caption
    two = ladder(
        [
            {
                "compound_id": "c",
                "target": "T",
                "rung": "biochemical",
                "px": 8.0,
                "assay_type": "Kd",
            },
            {"compound_id": "c", "target": "T", "rung": "phenotype", "px": 6.0},
        ]
    )
    r2 = tc.engagement_ladder(two, [("c", "T")])
    assert [g.rung for g in r2.meta["results"][0].rungs] == ["biochemical", "phenotype"]
    assert len(list(r2.spec.marks(LineMark))) == 1  # joined directly, nothing interpolated
    assert r2.spec.panels[0].y.tick_labels == ["Biochemical", "Phenotype"]


def test_errors():
    t = typical()
    with pytest.raises(TinctureDataError, match="No potencies"):
        tc.engagement_ladder(t, [("synthetic_a", "NOPE")])
    bad = t.data.copy()
    bad.loc[0, "rung"] = "biochem"
    from tincture.model.potency import PotencyTable

    with pytest.raises(TinctureDataError, match="Unknown rungs"):
        tc.engagement_ladder(PotencyTable(bad), PAIR)
    dup = t.data.copy()
    dup.loc[1, "rung"] = "biochemical"
    with pytest.raises(TinctureDataError, match="more than one value"):
        tc.engagement_ladder(PotencyTable(dup), PAIR)
    with pytest.raises(TinctureDataError, match="ATP Km must be positive"):
        atp_shift(0)


@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_engagement_baseline():
    fig = tc.engagement_ladder(
        typical(pheno_px=8.7, cell_px=6.5), PAIR, km_atp={"SYNK1": 2e-5}
    ).figure()
    check_colorblind(fig)
    return fig
