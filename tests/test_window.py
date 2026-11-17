"""4PL fitting and the Window plot: statistics, every §5.7 guard, baseline."""

from __future__ import annotations

import numpy as np
import pytest

import tincture as tc
from tests.fixtures.synthetic_curves import DOSES, synthetic_curve, synthetic_units
from tincture.errors import TinctureDataError
from tincture.model.dose_response import FourPL
from tincture.spec import BandMark, NoticeMark, TextMark
from tincture.stats.hill import (
    bootstrap_4pl,
    ec_drop,
    ec_relative,
    fit_4pl,
    median_of_units,
)
from tincture.theme import check_colorblind

POTENT = FourPL(bottom=0.05, top=1.0, log10_ec50=-8.0, hill=1.5)  # on-target, strong
TOXIC = FourPL(bottom=0.1, top=1.0, log10_ec50=-5.6, hill=2.0)  # kills only at > 1 µM
FLAT = FourPL(bottom=1.0, top=1.0, log10_ec50=-7.0, hill=1.0)


# --------------------------------------------------------------------------- statistics
def test_fit_recovers_parameters():
    dr = synthetic_curve(POTENT, noise=0.02)
    h = fit_4pl(dr.dose, dr.response)
    assert h.converged and h.r2 > 0.98
    assert h.fit.log10_ec50 == pytest.approx(-8.0, abs=0.1)
    assert h.fit.hill == pytest.approx(1.5, rel=0.2)


def test_effective_concentrations():
    f = FourPL(0.0, 1.0, -7.0, 1.0)
    assert ec_relative(f, 0.5) == pytest.approx(1e-7)
    assert ec_relative(f, 0.9) == pytest.approx(9e-7)  # EC90 = EC50 · 9^(1/h)
    assert ec_drop(f, 0.1) == pytest.approx(1e-7 / 9)  # viability 0.9 ⇔ 10% of span here
    assert ec_drop(FLAT, 0.1) == float("inf")


def test_bootstrap_interval_covers_truth():
    dr = synthetic_curve(POTENT, noise=0.04, seed=3)
    b = bootstrap_4pl(dr.dose, dr.response, n_boot=150)
    lo, hi = b.interval(lambda f: ec_relative(f, 0.9))
    assert lo < ec_relative(POTENT, 0.9) < hi and b.n_failed < 15


def test_median_of_units_and_unit_bootstrap():
    dr = synthetic_units(TOXIC)
    x, m = median_of_units(dr.dose, dr.response, dr.replicate)
    assert np.allclose(x, DOSES) and m[0] > 0.9 > m[-1]
    b = bootstrap_4pl(dr.dose, dr.response, n_boot=60, cluster=dr.replicate)
    assert b.params.shape[0] >= 50


# --------------------------------------------------------------------------- guards
def _window(on, tox, **kw):
    return tc.window(on, tox, n_boot=60, compound="synthetic_cpd", **kw)


def test_window_exists_and_is_drawn_with_context():
    res = _window(synthetic_curve(POTENT), synthetic_curve(TOXIC, seed=1))
    r = res.meta["result"]
    assert r.exists is True and r.on.ec < r.limiting.ec
    assert any(isinstance(m, BandMark) and m.inner_label == "window" for m in res.spec.marks())
    texts = [m.text for m in res.spec.marks(TextMark)]
    assert any("SYNTH-1" in t and "synthetic 3 d" in t for t in texts)  # assay context printed
    assert "a window of" in res.caption and "No exposure (free Cmax)" in res.caption


def test_no_band_when_on_target_effect_below_half():
    weak = FourPL(bottom=0.7, top=1.0, log10_ec50=-8.0, hill=1.5)  # 30% max effect
    res = _window(synthetic_curve(weak), synthetic_curve(TOXIC, seed=1))
    assert res.meta["result"].exists is None
    assert not [m for m in res.spec.marks(BandMark) if m.inner_label == "window"]
    assert "maximal effect 3" in res.caption and "No window is drawn" in res.caption
    assert any(isinstance(m, NoticeMark) for m in res.spec.marks())


def test_no_band_when_fit_quality_poor():
    res = _window(synthetic_curve(POTENT, noise=0.35, seed=5), synthetic_curve(TOXIC, seed=1))
    r = res.meta["result"]
    assert r.exists is None and any("R²" in p for p in r.on.problems)


def test_flat_toxicity_is_censored_not_failed():
    res = _window(synthetic_curve(POTENT), synthetic_curve(FLAT, seed=2))
    r = res.meta["result"]
    assert r.exists is True and r.limiting.censored and not np.isfinite(r.limiting.ec)
    assert "above the top tested dose" in res.caption and "at least" in res.caption


def test_no_window_when_toxicity_comes_first():
    early_tox = FourPL(bottom=0.1, top=1.0, log10_ec50=-8.6, hill=2.0)
    res = _window(synthetic_curve(POTENT), synthetic_curve(early_tox, seed=4))
    assert res.meta["result"].exists is False and res.caption.startswith("No window")


def test_most_sensitive_toxicity_curve_limits_the_window():
    mild = FourPL(bottom=0.1, top=1.0, log10_ec50=-5.0, hill=2.0)
    res = _window(synthetic_curve(POTENT), [synthetic_curve(mild, seed=1),
                                            synthetic_curve(TOXIC, seed=2)])  # fmt: skip
    r = res.meta["result"]
    assert r.limiting is r.tox[1]


def test_median_of_units_toxicity_and_cmax_band():
    res = _window(synthetic_curve(POTENT), synthetic_units(TOXIC), cmax=(1e-8, 5e-8))
    assert any(isinstance(m, BandMark) and m.inner_label == "free Cmax" for m in res.spec.marks())
    assert "No exposure" not in res.caption
    assert "interquartile" in res.spec.footer[0]
    with pytest.raises(TinctureDataError, match="cmax"):
        _window(synthetic_curve(POTENT), synthetic_units(TOXIC), cmax=(5.0, 1.0))
    assert res.classic is not None and res.classic.plot == "window_classic"


@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_window_baseline():
    fig = _window(synthetic_curve(POTENT), synthetic_curve(TOXIC, seed=1)).figure()
    check_colorblind(fig)
    return fig
