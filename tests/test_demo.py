"""Demo plot: degenerate inputs, weak-result path and the visual baseline."""

from __future__ import annotations

import numpy as np
import pytest

import tincture as tc
from tincture.errors import TinctureDataError
from tincture.spec import NoticeMark
from tincture.theme import check_colorblind


def test_result_has_caption_methods_classic(dose_table):
    res = tc.demo(dose_table, focus="synthetic_cpd2")
    assert res.caption.startswith("3 compounds at 4 concentrations; synthetic_cpd2")
    assert "Software:" in res.methods
    assert res.classic is not None and res.meta["focus"] == "synthetic_cpd2"
    assert res.warnings == []


def test_missing_columns_raise_with_fix(dose_table):
    with pytest.raises(TinctureDataError, match="How to fix") as e:
        tc.demo(dose_table.drop(columns=["dose"]))
    assert e.value.fix is not None


def test_no_vehicle(dose_table):
    with pytest.raises(TinctureDataError, match="vehicle"):
        tc.demo(dose_table.loc[dose_table["dose"] > 0])


def test_negative_dose(dose_table):
    d = dose_table.copy()
    d.loc[0, "dose"] = -1
    with pytest.raises(TinctureDataError, match="Negative"):
        tc.demo(d)


def test_all_nan_values(dose_table):
    d = dose_table.copy()
    d["value"] = np.nan
    with pytest.raises(TinctureDataError):
        tc.demo(d)


def test_unknown_focus(dose_table):
    with pytest.raises(TinctureDataError, match="Focus compound"):
        tc.demo(dose_table, focus="nope")


def test_single_compound(dose_table):
    d = dose_table.loc[dose_table["compound_id"].isin(["vehicle", "synthetic_cpd1"])]
    res = tc.demo(d)
    res.figure()


def test_single_replicate_says_so(dose_table_single_rep):
    res = tc.demo(dose_table_single_rep)
    assert "Note: only one replicate" in res.caption
    assert res.warnings and res.meta["single_replicate"]
    assert any("no CI" in m.text for m in res.spec.marks(NoticeMark))


@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_demo_baseline(dose_table):
    fig = tc.demo(dose_table, focus="synthetic_cpd1").figure("nature_double")
    check_colorblind(fig)  # accessibility check on every baseline
    return fig


@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_kitchen_sink_baseline():
    from tests.test_render import kitchen_sink

    fig = tc.render_mpl(kitchen_sink(), "nature_double")
    check_colorblind(fig)
    return fig
