"""Synergy statistics: Bliss / HSA hand-computed on a 3×3 example, ceiling mask, ZIP,
uncertainty modes and design guards (HANDOVER §5.11, §7)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tests.fixtures.synthetic_combo import anchored_matrix, synthetic_matrix
from tincture.errors import TinctureDataError
from tincture.model.combo import ComboMatrix
from tincture.stats.synergy import synergy

VA = [0.9, 0.6, 0.3]  # drug A alone at 1, 10, 100 nM
VB = [0.8, 0.5, 0.2]  # drug B alone at 1, 10, 100 nM
OBS = [  # combination viability, rows = B (low → high), cols = A (low → high)
    [0.70, 0.40, 0.20],
    [0.40, 0.20, 0.10],
    [0.15, 0.10, 0.04],
]
LEVELS = [1e-9, 1e-8, 1e-7]


def hand_matrix() -> ComboMatrix:
    rows = []
    for i, a in enumerate(LEVELS):
        rows.append({"conc_a": a, "conc_b": 0.0, "viability": VA[i]})
    for j, b in enumerate(LEVELS):
        rows.append({"conc_a": 0.0, "conc_b": b, "viability": VB[j]})
        for i, a in enumerate(LEVELS):
            rows.append({"conc_a": a, "conc_b": b, "viability": OBS[j][i]})
    return ComboMatrix(pd.DataFrame(rows), "A", "B")


def test_bliss_matches_hand_computed_3x3():
    r = synergy(hand_matrix(), "bliss", single_agent="observed")
    # expected viability vA·vB, e.g. B 1 nM row: 0.72, 0.48, 0.24
    exp = np.array([[0.72, 0.48, 0.24], [0.45, 0.30, 0.15], [0.18, 0.12, 0.06]])
    np.testing.assert_allclose(r.expected, exp)
    exc = np.array([[0.02, 0.08, 0.04], [0.05, 0.10, 0.05], [0.03, 0.02, 0.02]])
    np.testing.assert_allclose(r.excess, exc, atol=1e-12)


def test_hsa_matches_hand_computed_3x3():
    r = synergy(hand_matrix(), "hsa", single_agent="observed")
    exp = np.array([[0.8, 0.6, 0.3], [0.5, 0.5, 0.3], [0.2, 0.2, 0.2]])
    np.testing.assert_allclose(r.expected, exp)
    np.testing.assert_allclose(r.excess, exp - np.array(OBS), atol=1e-12)


def test_ceiling_mask_applied():
    r = synergy(hand_matrix(), "bliss", single_agent="observed")
    # expected inhibition 1 − 0.06 = 0.94 > 0.9 only at (B 100 nM, A 100 nM)
    assert r.masked.sum() == 1 and r.masked[2, 2]
    assert not np.isfinite(r.q[2, 2])  # never tested
    loose = synergy(hand_matrix(), "bliss", single_agent="observed", ceiling=0.95)
    assert not loose.masked.any()


def test_null_surface_has_no_significant_cells():
    cm = synthetic_matrix(synergy=0.0)
    for model in ("bliss", "zip"):
        r = synergy(cm, model)
        assert r.uncertainty == "replicates"
        assert r.n_significant == 0
        assert abs(r.score) < 0.03


def test_injected_synergy_found_inside_its_region():
    cm = synthetic_matrix(synergy=0.25, region=(3, 6))
    inside = np.zeros((7, 7), bool)
    inside[3:6, 3:6] = True
    for model in ("bliss", "zip"):
        r = synergy(cm, model)
        assert r.n_significant >= 5
        assert not (r.significant & ~inside).any()


def test_hsa_is_benefit_not_synergy():
    # a Bliss-independent pair still beats its best single agent: HSA excess > 0
    r = synergy(synthetic_matrix(synergy=0.0), "hsa")
    assert r.score > 0.02


def test_single_replicate_without_noise_claims_nothing():
    cm = synthetic_matrix(synergy=0.25, n_rep=1)
    r = synergy(cm, "bliss")
    assert r.uncertainty == "none" and r.n_significant == 0
    assert np.isnan(r.lo).all()
    r2 = synergy(cm, "bliss", noise_sd=0.02)
    assert r2.uncertainty == "noise_sd" and r2.n_significant >= 5


def test_min_excess_threshold():
    cm = synthetic_matrix(synergy=0.08, sd=0.005)
    assert synergy(cm, "bliss").n_significant == 0  # clear but below 10%
    assert synergy(cm, "bliss", min_excess=0.05).n_significant >= 5


def test_anchored_design_refused():
    cm = anchored_matrix()
    assert cm.design == "anchored"
    with pytest.raises(TinctureDataError, match="anchored") as e:
        synergy(cm)
    assert "combination_lift" in str(e.value)


def test_partial_grid_refused_and_bad_inputs():
    d = synthetic_matrix().data
    d = d[~((d["conc_a"] == d["conc_a"].max()) & (d["conc_b"] == d["conc_b"].max()))]
    cm = ComboMatrix(d, "A", "B")
    assert cm.design == "partial"
    with pytest.raises(TinctureDataError, match="incomplete"):
        synergy(cm)
    with pytest.raises(TinctureDataError, match="Unknown synergy model"):
        synergy(synthetic_matrix(), "loewe")  # type: ignore[arg-type]
    with pytest.raises(TinctureDataError, match="missing columns"):
        ComboMatrix(pd.DataFrame({"conc_a": [1.0]}), "A", "B")
    with pytest.raises(TinctureDataError, match="≥ 0"):
        ComboMatrix(pd.DataFrame({"conc_a": [-1.0], "conc_b": [0], "viability": [1]}), "A", "B")


def test_zip_needs_four_levels():
    with pytest.raises(TinctureDataError, match="ZIP needs at least 4"):
        synergy(hand_matrix(), "zip")


def test_grid_noise_sd_recovers_known_noise():
    from tincture.stats.synergy import grid_noise_sd, matrix_noise_sd

    ests = []
    for seed in range(20):
        cm = synthetic_matrix(n_rep=1, sd=0.05, seed=seed)
        ests.append(matrix_noise_sd(cm))
    assert all(e is not None for e in ests)
    assert 0.035 < float(np.median(ests)) < 0.07  # true SD 0.05
    # saturated grid: no estimate rather than a falsely small one
    assert grid_noise_sd(np.full((7, 7), 0.02)) is None


def test_matrix_noise_option_single_replicate():
    cm = synthetic_matrix(synergy=0.25, n_rep=1, sd=0.03)
    r = synergy(cm, "bliss", noise_sd="matrix")
    assert r.uncertainty == "noise_sd" and "second differences" in r.noise_method
    assert r.n_significant >= 5
    null = synergy(synthetic_matrix(synergy=0.0, n_rep=1, sd=0.03, seed=3), noise_sd="matrix")
    assert null.n_significant <= 1
