"""Onset Timeline: interval-censored onset (HANDOVER §7), guards from §5.12, baseline."""

from __future__ import annotations

import numpy as np
import pytest

import tincture as tc
from tests.fixtures.synthetic_timecourse import onset_long, typical
from tincture.errors import TinctureDataError
from tincture.spec import IntervalMark, PointMark
from tincture.stats.onset import onset_table
from tincture.theme import check_colorblind


def test_interval_censored_onset_2_6_24():
    d = onset_long([2, 6, 24], {"f": [(0.5, 0.6), (3.0, 0.01), (4.0, 0.001)]})
    (o,), _ = onset_table(d)
    assert o.onset == (2.0, 6.0)
    assert o.peak_time == 24.0 and o.sustained


def test_active_at_first_time_point_is_left_censored():
    d = onset_long([3, 6, 12], {"f": [(-3.0, 0.01), (-1.0, 0.5), (0.1, 0.9)]})
    (o,), _ = onset_table(d)
    assert o.onset == (0.0, 3.0) and o.direction == -1
    assert o.resolution == (3.0, 6.0) and not o.sustained


def test_threshold_and_p_to_q():
    d = onset_long([3, 6, 12], {"f": [(1.0, 0.001), (2.5, 0.001), (3.0, 0.001)]})
    (o,), _ = onset_table(d, threshold=2.0)
    assert o.onset == (3.0, 6.0)
    d2 = d.drop(columns="q").assign(p=[0.001, 0.2, 0.3])
    (o2,), x = onset_table(d2)
    assert "q" in x and o2.onset == (0.0, 3.0)


def test_onset_never_drawn_as_a_point():
    res = tc.onset_timeline(typical())
    iv = [m for m in res.spec.marks() if isinstance(m, IntervalMark)]
    assert len(iv) == 3  # one bracket per active pathway
    for m in iv:
        assert np.all(np.asarray(m.hi) > np.asarray(m.lo))
    # dots are peaks only, one per active pathway
    assert sum(np.size(m.x) for m in res.spec.marks(PointMark)) == 3


def test_order_caption_and_never_active():
    res = tc.onset_timeline(typical())
    assert res.meta["order"][:3] == ["direct", "transient", "stress"]
    assert "First to respond (onset in (0, 3] h): direct (down" in res.caption
    assert "never become active: quiet" in res.caption
    assert any("Never active: quiet" in f for f in res.spec.footer)
    assert "true onset can lie anywhere inside the bracket" in res.caption


def test_two_time_points_warn_one_refused():
    d = onset_long([6, 24], {"f": [(0.1, 0.9), (3.0, 0.01)]})
    res = tc.onset_timeline(d)
    assert any("only 2 time points" in w for w in res.warnings)
    with pytest.raises(TinctureDataError, match="at least 2 time points"):
        tc.onset_timeline(onset_long([6], {"f": [(3.0, 0.01)]}))
    with pytest.raises(TinctureDataError, match="needs columns"):
        tc.onset_timeline(typical().drop(columns="q"))


def test_horizon_only_for_dense_series():
    with pytest.raises(TinctureDataError, match="dense time series"):
        tc.onset_timeline(typical(), style="horizon")
    times = [float(t) for t in range(0, 48, 6)]
    rows = {f: [(np.sin(t / 10 + k), 0.01) for t in times] for k, f in enumerate("abc")}
    res = tc.onset_timeline(onset_long(times, rows), style="horizon")
    assert res.spec.plot == "onset_horizon"


@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_onset_baseline():
    fig = tc.onset_timeline(typical()).figure()
    check_colorblind(fig)
    return fig
