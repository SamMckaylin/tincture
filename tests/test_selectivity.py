"""PotencyTable and Selectivity Spectrum: ranking, exposure band, mixed-assay guard."""

from __future__ import annotations

import pandas as pd
import pytest

import tincture as tc
from tests.fixtures.synthetic_potency import synthetic_panel
from tincture.errors import TinctureDataError
from tincture.model.potency import PotencyTable, px_from_molar
from tincture.spec import BandMark, PointMark
from tincture.theme import check_colorblind


def test_potency_table_validation():
    assert px_from_molar(1e-9) == pytest.approx(9.0)
    with pytest.raises(TinctureDataError, match="missing columns"):
        PotencyTable(pd.DataFrame({"compound_id": ["a"]}))
    bad = pd.DataFrame({"compound_id": ["a"], "target": ["t"], "px": [7.0], "assay_type": ["pIC"]})
    with pytest.raises(TinctureDataError, match="Unknown assay"):
        PotencyTable(bad)
    t = PotencyTable.from_molar(
        pd.DataFrame(
            {"compound_id": ["a"], "target": ["t"], "value_molar": [1e-8], "assay_type": ["Kd"]}
        )
    )
    assert t.data["px"].iloc[0] == pytest.approx(8.0) and t.data["censored"].iloc[0] == ""


def test_ranking_coverage_and_fold_rule_without_cmax():
    res = tc.selectivity(
        synthetic_panel(), "synthetic_cpd", "SYN_TARGET", panel_size=518, panel_name="kinases"
    )
    r = res.meta["ranked"]
    assert r.iloc[0]["target"] == "SYN_TARGET"  # strongest first
    assert list(r[r["censored"] == ">"]["rank"]) == list(range(len(r) - 9, len(r) + 1))
    # no exposure: off-targets within 10-fold of the intended target (pX ≥ 7.5)
    assert res.meta["off_targets"] == ["SYN_OFF0", "SYN_OFF1"]
    assert "within 10-fold" in res.caption and "No exposure" in res.caption
    assert res.spec.footer[0] == "Panel coverage: 74 of ~518 kinases tested"
    hollow = [m for m in res.spec.marks(PointMark) if m.hollow is True]
    assert hollow and len(hollow[0].x) == 10  # censored drawn hollow
    assert res.classic is not None


def test_cmax_band_defines_off_targets():
    res = tc.selectivity(synthetic_panel(), "synthetic_cpd", "SYN_TARGET", cmax=(3e-8, 3e-7))
    assert any(
        isinstance(m, BandMark) and "free Cmax" in (m.inner_label or "") for m in res.spec.marks()
    )
    # above the band = potency stronger than 30 nM (pX > 7.52)
    assert res.meta["off_targets"] == ["SYN_OFF0", "SYN_OFF1"]
    res2 = tc.selectivity(synthetic_panel(), "synthetic_cpd", "SYN_TARGET", cmax=(2e-7, 1e-6))
    assert res2.meta["off_targets"] == ["SYN_OFF0", "SYN_OFF1", "SYN_OFF2"]


def test_mixed_assay_types_refused_unless_allowed():
    a = synthetic_panel(assay="Kd").data
    b = synthetic_panel(assay="IC50", seed=1).data.head(5).assign(target=lambda d: d.target + "_c")
    mixed = PotencyTable(pd.concat([a, b]))
    with pytest.raises(TinctureDataError, match="Mixed assay types"):
        tc.selectivity(mixed, "synthetic_cpd", "SYN_TARGET")
    res = tc.selectivity(mixed, "synthetic_cpd", "SYN_TARGET", allow_mixed=True)
    markers = {m.marker for m in res.spec.marks(PointMark)}
    assert len(markers) >= 2 and "mixed" in res.caption


def test_unknown_intended_target():
    with pytest.raises(TinctureDataError, match="not in the panel"):
        tc.selectivity(synthetic_panel(), "synthetic_cpd", "NOPE")


@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_selectivity_baseline():
    fig = tc.selectivity(
        synthetic_panel(),
        "synthetic_cpd",
        "SYN_TARGET",
        cmax=(3e-8, 3e-7),
        panel_size=518,
        panel_name="kinases",
    ).figure()
    check_colorblind(fig)
    return fig
