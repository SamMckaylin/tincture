"""Combination Surface and Lift: guards from HANDOVER §5.11 and baselines."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import tincture as tc
from tests.fixtures.synthetic_combo import anchored_matrix, synthetic_matrix
from tincture.errors import TinctureDataError
from tincture.spec import ContourMark, LineMark, PointMark, PolygonMark, RectMark
from tincture.theme import check_colorblind


def _outline(res) -> list[LineMark]:
    return [
        m for m in res.spec.panel("surface").marks if isinstance(m, LineMark) and m.color == "focus"
    ]


def test_surface_prints_model_name():
    for model, name in (("bliss", "Bliss independence"), ("zip", "ZIP"), ("hsa", "HSA")):
        res = tc.combination_surface(synthetic_matrix(synergy=0.25), model=model)
        assert name in res.spec.footer[0]
        assert name in res.caption or model == "hsa"
        assert res.spec.panel("excess").title.startswith("excess over")


def test_surface_ceiling_cells_hatched():
    res = tc.combination_surface(synthetic_matrix(synergy=0.25))
    r = res.meta["result"]
    assert r.masked.sum() > 0
    hatched = [m for m in res.spec.panel("surface").marks if isinstance(m, RectMark) and m.hatch]
    assert sum(np.size(m.x0) for m in hatched) == r.masked.sum()
    assert "hatched" in res.caption and "not interpretable" in res.caption


def test_surface_outlines_significant_region_only():
    res = tc.combination_surface(synthetic_matrix(synergy=0.25))
    assert res.meta["result"].n_significant >= 5 and _outline(res)
    null = tc.combination_surface(synthetic_matrix(synergy=0.0))
    assert not _outline(null)
    assert "no cell shows significant excess" in null.caption


def test_surface_single_replicate_draws_no_region_and_says_so():
    res = tc.combination_surface(synthetic_matrix(synergy=0.25, n_rep=1))
    assert not _outline(res)
    assert "without significance" in res.caption


def test_surface_contours_stay_inside_grid():
    cm = synthetic_matrix()
    res = tc.combination_surface(cm)
    (c,) = [m for m in res.spec.marks() if isinstance(m, ContourMark)]
    np.testing.assert_allclose(c.x, cm.levels_a)
    np.testing.assert_allclose(c.y, cm.levels_b)
    assert np.shape(c.z_values) == (len(cm.levels_b), len(cm.levels_a))


def test_surface_cmax_box_and_reachability():
    res = tc.combination_surface(synthetic_matrix(synergy=0.25))
    assert "clinical reachability is not assessed" in res.caption
    box = tc.combination_surface(
        synthetic_matrix(synergy=0.25), cmax_a=(3e-8, 1e-7), cmax_b=(3e-8, 1e-7)
    )
    assert any(isinstance(m, PolygonMark) for m in box.spec.panel("surface").marks)
    assert "inside the clinically reachable box" in box.caption
    with pytest.raises(TinctureDataError, match="molar"):
        tc.combination_surface(synthetic_matrix(), cmax_a=(1.0, 2.0), cmax_b=(1e-7, 1e-6))


def test_surface_refuses_anchored_design():
    with pytest.raises(TinctureDataError, match="anchored") as e:
        tc.combination_surface(anchored_matrix())
    assert "combination_lift" in str(e.value)


def _lift_df(seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(30):
        grp = "pos" if i < 18 else "neg"
        a, b = rng.uniform(0.5, 0.95, 2)
        gain = 0.35 if grp == "pos" else 0.02
        for _rep in range(2):
            rows.append(
                {
                    "cell_line": f"L{i}",
                    "single_a": a,
                    "single_b": b,
                    "combo": min(a, b) - gain + rng.normal(0, 0.02),
                    "grp": grp,
                    "fit_ok": i != 0,
                }
            )
    return pd.DataFrame(rows)


def test_lift_threshold_and_model_printed():
    res = tc.combination_lift(_lift_df(), "A", "B", group="grp", model_note="Bliss for synergy.")
    assert "ΔEmax ≥ 20%" in res.spec.footer[0] and "Bliss" in res.spec.footer[0]
    assert "≥ 20% above the best single agent" in res.caption
    st = res.meta["stats"]["None"]
    assert st["k_pos"] == 17 and st["k_neg"] == 0
    assert st["fisher_p"] < 0.001


def test_lift_failed_fits_hollow_and_not_counted():
    res = tc.combination_lift(_lift_df(), "A", "B", group="grp")
    hollow = [m for m in res.spec.marks(PointMark) if np.any(m.hollow)]
    assert len(hollow) == 1 and np.size(hollow[0].x) == 1
    assert res.meta["stats"]["None"]["n"] == 29
    assert "1 cell line(s) with failed fits are hollow" in res.caption


def test_lift_input_errors():
    with pytest.raises(TinctureDataError, match="missing columns"):
        tc.combination_lift(pd.DataFrame({"cell_line": ["x"]}), "A", "B")
    with pytest.raises(TinctureDataError, match="threshold"):
        tc.combination_lift(_lift_df(), "A", "B", threshold=20)


def test_classics_exist():
    assert tc.combination_surface(synthetic_matrix()).classic is not None
    assert tc.combination_lift(_lift_df(), "A", "B").classic is not None


@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_combination_surface_baseline():
    fig = tc.combination_surface(
        synthetic_matrix(synergy=0.25), cmax_a=(3e-8, 1e-7), cmax_b=(3e-8, 1e-7)
    ).figure()
    check_colorblind(fig)
    return fig


@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_combination_lift_baseline():
    fig = tc.combination_lift(_lift_df(), "A", "B", group="grp").figure()
    check_colorblind(fig)
    return fig
