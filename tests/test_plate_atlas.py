"""Plate Atlas: Z′, edge effect, drift, encoding guards, visual baseline."""

from __future__ import annotations

import numpy as np
import pytest

import tincture as tc
from tests.fixtures.synthetic_plates import synthetic_plates
from tincture.errors import TinctureDataError
from tincture.spec import HeatmapMark, PointMark, RectMark, RefLineMark
from tincture.stats.zprime import Selector, plate_qc, zprime
from tincture.theme import check_colorblind


def test_zprime_formula():
    pos, neg = np.array([10.0, 12.0, 14.0]), np.array([100.0, 102.0, 104.0])
    assert zprime(pos, neg) == pytest.approx(1 - 3 * (2 + 2) / 90)
    assert zprime(pos, neg, robust=True) == pytest.approx(1 - 3 * (1.4826 * 2 * 2) / 90)
    assert np.isnan(zprime(np.array([1.0]), neg))
    assert zprime(neg, neg) == -np.inf


def test_good_and_bad_plates_flagged():
    pt = synthetic_plates(bad_plates=(2,))
    qc = plate_qc(pt)
    t = qc.per_plate.set_index("plate")
    assert t.loc["P02", "fails"] and not t.drop("P02")["fails"].any()


def test_edge_effect_detected_when_injected():
    t = plate_qc(synthetic_plates(edge_factor=0.6)).per_plate
    assert t["edge_ratio"].median() == pytest.approx(0.6, abs=0.05)
    assert (t["edge_depth"] == 0).all()
    t0 = plate_qc(synthetic_plates(edge_factor=1.0)).per_plate
    assert t0["edge_ratio"].median() == pytest.approx(1.0, abs=0.05)


def test_unseeded_outer_ring_uses_next_ring():
    qc = plate_qc(synthetic_plates(edge_factor=0.7, unseeded_outer=True))
    assert (qc.per_plate["edge_depth"] == 1).all() and qc.n_unseeded_outer == 6
    assert qc.per_plate["edge_ratio"].median() == pytest.approx(0.7, abs=0.06)


def test_drift_detected_and_cell_line_confound_removed():
    qc = plate_qc(synthetic_plates(n_plates=12, drift_per_plate=0.03))
    assert qc.drift_rho > 0.9 and qc.drift_change == pytest.approx(0.33, abs=0.08)
    # bright and dim lines alternating, no real drift: normalising within line removes it
    pt = synthetic_plates(n_plates=12, cell_lines=("A", "B"), line_brightness=(1.0, 3.0))
    qc2 = plate_qc(pt)
    assert qc2.drift_basis == "same cell line"
    assert abs(qc2.drift_change) < 0.05
    assert qc2.per_plate["neg_level"].between(0.9, 1.1).all()


def test_missing_controls_raise_with_fix():
    pt = synthetic_plates()
    with pytest.raises(TinctureDataError, match="negative-control"):
        plate_qc(pt, neg=Selector(well_types=("nope",)))
    with pytest.raises(TinctureDataError, match="positive-control"):
        plate_qc(pt, pos=Selector(tags=("nope",)))


def test_failing_plates_terracotta_and_threshold_line():
    res = tc.plate_atlas(synthetic_plates(bad_plates=(1, 4)))
    fails = [
        m for m in res.spec.panel("zprime").marks if isinstance(m, RectMark) and m.color == "focus"
    ]
    assert len(fails) == 1 and np.size(fails[0].x0) == 2
    lines = [m for m in res.spec.panel("zprime").marks if isinstance(m, RefLineMark)]
    assert lines[0].value == 0.5 and "2 (33%) fall below 0.5" in res.caption
    assert res.meta["shown_plate"] in {"P01", "P04"} and res.meta["why"] == "worst Z′"


def test_density_view_for_1536_and_shapes_for_384():
    r1536 = tc.plate_atlas(synthetic_plates(n_plates=2, n_wells=1536))
    assert (
        r1536.meta["density_view"] and r1536.spec.panel("plate").marks[0].__class__ is HeatmapMark
    )
    r384 = tc.plate_atlas(synthetic_plates(n_plates=2, n_wells=384))
    markers = {m.marker for m in r384.spec.panel("plate").marks if isinstance(m, PointMark)}
    assert markers == {"o", "s"}  # controls circles, compounds squares


def test_missing_wells_reported_not_zero():
    res = tc.plate_atlas(synthetic_plates(n_plates=2, n_wells=1536, drop_wells=50))
    assert "50 well positions are not reported" in res.caption
    m = res.spec.panel("plate").marks[0]
    assert np.isnan(np.asarray(m.values)).sum() > 0


def test_too_few_controls_gives_nan_and_note():
    pt = synthetic_plates(n_plates=3)
    d = pt.data
    keep = ~((d["plate"] == "P00") & (d["well_type"] == "pos_ctrl") & (d["row"] > 3))
    from tincture.model import PlateTable

    res = tc.plate_atlas(PlateTable(d[keep].reset_index(drop=True), 384))
    assert "could not be computed for 1 plate" in res.caption


def test_classic_and_methods():
    res = tc.plate_atlas(synthetic_plates())
    assert res.classic is not None and res.classic.meta["classic_of"] == "plate_atlas"
    assert "Zhang" in res.methods and "outermost seeded ring" in res.methods


@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_plate_atlas_baseline():
    fig = tc.plate_atlas(
        synthetic_plates(n_wells=96, n_plates=8, bad_plates=(3,), edge_factor=0.7)
    ).figure()
    check_colorblind(fig)
    return fig


def test_many_plates_switch_to_binned_summary():
    from tincture.plots import plate_atlas as pa

    old = pa.BIN_THRESHOLD
    pa.BIN_THRESHOLD = 5
    try:
        res = tc.plate_atlas(synthetic_plates(n_plates=8, n_wells=96, bad_plates=(2,)))
    finally:
        pa.BIN_THRESHOLD = old
    from tincture.spec import BandMark

    assert any(isinstance(m, BandMark) for m in res.spec.panel("edge").marks)
    fails = [m for m in res.spec.panel("zprime").marks if isinstance(m, PointMark)]
    assert fails and fails[0].color == "focus"  # failing plates still shown individually
    assert "equal-count bins" in res.caption
