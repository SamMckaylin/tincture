from __future__ import annotations

import numpy as np
import pytest

from tincture.errors import TinctureSpecError
from tincture.spec import (
    BandMark,
    HeatmapMark,
    IntervalMark,
    LineMark,
    Panel,
    PlotSpec,
    PointMark,
    Scale,
)


def _spec(*marks, **kw):
    return PlotSpec(plot="t", panels=[Panel(marks=list(marks), **kw)])


def test_point_lengths_must_match():
    with pytest.raises(TinctureSpecError, match="different lengths"):
        _spec(PointMark(x=[1, 2, 3], y=[1, 2])).validate()


def test_colour_list_length_checked():
    with pytest.raises(TinctureSpecError, match="colour list"):
        _spec(PointMark(x=[1, 2], y=[1, 2], color=["focus"])).validate()


def test_interval_lo_gt_hi():
    with pytest.raises(TinctureSpecError, match="lo > hi"):
        _spec(IntervalMark(pos=[0], lo=[2], hi=[1])).validate()


def test_interval_allows_nan():
    _spec(IntervalMark(pos=[0, 1], lo=[np.nan, 0], hi=[np.nan, 1])).validate()


def test_categorical_scale_needs_categories():
    with pytest.raises(TinctureSpecError, match="categories"):
        Scale("categorical")


def test_log_domain_positive():
    with pytest.raises(TinctureSpecError, match="positive"):
        Scale("log", domain=(0, 1))


def test_tick_label_mismatch():
    with pytest.raises(TinctureSpecError):
        Scale(ticks=[1, 2], tick_labels=["a"])


def test_diverging_heatmap_needs_clip():
    with pytest.raises(TinctureSpecError, match="clip"):
        _spec(HeatmapMark(values=np.zeros((2, 2)))).validate()


def test_heatmap_mask_shape():
    with pytest.raises(TinctureSpecError, match="mask"):
        _spec(HeatmapMark(values=np.zeros((2, 2)), mask=np.zeros((3, 2)), clip=1)).validate()


def test_band_span_requires_scalars():
    with pytest.raises(TinctureSpecError, match="scalar"):
        _spec(BandMark(lo=[0, 1], hi=[1, 2], span=True)).validate()


def test_empty_spec_raises():
    with pytest.raises(TinctureSpecError, match="no panels"):
        PlotSpec(plot="x").validate()


def test_overlapping_panels():
    a = Panel(key="a", row=0, col=0, colspan=2)
    b = Panel(key="b", row=0, col=1)
    with pytest.raises(TinctureSpecError, match="overlaps"):
        PlotSpec(plot="x", panels=[a, b], ncols=2).validate()


def test_hollow_length_checked():
    with pytest.raises(TinctureSpecError, match="different lengths"):
        _spec(PointMark(x=[1, 2], y=[1, 2], hollow=[True, False, True])).validate()


def test_panel_outside_grid():
    with pytest.raises(TinctureSpecError, match="outside"):
        PlotSpec(plot="x", panels=[Panel(col=2)], ncols=2).validate()


def test_polar_requires_angle_radius():
    with pytest.raises(TinctureSpecError, match="polar"):
        _spec(projection="polar").validate()


def test_auto_layout_row_major():
    ps = [Panel(key=str(i)) for i in range(5)]
    PlotSpec(plot="x", panels=ps, nrows=2, ncols=3).validate()
    assert [(p.row, p.col) for p in ps] == [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]


def test_marks_filter_and_panel_lookup():
    spec = _spec(PointMark(x=[1], y=[1]), LineMark(x=[1, 2], y=[1, 2]), key="main")
    assert len(spec.marks(PointMark)) == 1
    assert spec.panel("main").key == "main"
    with pytest.raises(KeyError):
        spec.panel("nope")
