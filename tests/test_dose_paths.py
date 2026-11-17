"""Dose Paths: bend test (straight vs injected bend), axes, guards, encoding, baseline."""

from __future__ import annotations

import numpy as np
import pytest

import tincture as tc
from tincture.errors import TinctureDataError
from tincture.model.reference import library_from_signatures
from tincture.model.signature import Signature
from tincture.spec import EllipseMark, LineMark, NoticeMark, TextMark
from tincture.stats.paths import PathAxes, dose_path, signed_turns
from tincture.stats.signatures import SignatureSet
from tincture.theme import check_colorblind

P = 300
FEATS = np.array([f"g{i}" for i in range(P)])
STRESS = {f"g{i}" for i in range(40)}  # first 40 genes are the "stress" set
DOSES = [1e-8, 1e-7, 1e-6, 1e-5]


def _dirs(seed: int = 0) -> np.ndarray:
    """Orthonormal class directions with no net stress loading."""
    rng = np.random.default_rng(seed)
    D = rng.normal(size=(3, P))
    D[:, :40] -= D[:, :40].mean(axis=1, keepdims=True)
    D -= D.mean(axis=1, keepdims=True)
    q, _ = np.linalg.qr(D.T)
    return q.T


def _stress_dir() -> np.ndarray:
    v = np.zeros(P)
    v[:40] = 1.0
    return v - v.mean()


def _sig(cid: str, dose: float, eff: np.ndarray, se: float = 0.02) -> Signature:
    return Signature(
        features=FEATS,
        effect=eff,
        se=np.full(P, se),
        compound_id=cid,
        dose=dose,
        meta={"compound_name": cid.replace("synthetic_", "")},
    )


def _world(extra: dict[str, list[tuple[float, np.ndarray, float]]] | None = None):
    """Reference classes A/B/C (4 members × 4 doses) plus query compounds."""
    rng = np.random.default_rng(1)
    D = _dirs()
    sigs, top, classes = [], {}, {}
    for name, d in zip("ABC", D, strict=True):
        ids = [f"synthetic_{name}{i}" for i in range(4)]
        classes[name] = ids
        for cid in ids:
            for k, dose in enumerate(DOSES):
                e = (2.0 * 2**k) * d + rng.normal(0, 0.02, P)
                sigs.append(_sig(cid, dose, e))
                top[cid] = e
    for cid, rows in (extra or {}).items():
        for dose, e, se in rows:
            sigs.append(_sig(cid, dose, e + rng.normal(0, se, P), se=se))
    lib = library_from_signatures(top, FEATS, classes, version="synthetic v1", source="synth")
    ss = SignatureSet(
        signatures=sigs,
        stratum={"cell_line": "synthetic"},
        features=FEATS,
        df_prior=4.0,
        df_residual=4.0,
        n_vehicle=8,
        method="synthetic",
    )
    return ss, lib, D


def _straight(D):
    return [(dose, (1.5 * 2**k) * D[0], 0.02) for k, dose in enumerate(DOSES)]


def _bent(D, se: float = 0.02):
    S = _stress_dir()
    pts = [3 * D[0], 6 * D[0], 9 * D[0] + 1.2 * S, 12 * D[0] + 2.4 * S]
    return [(dose, e, se) for dose, e in zip(DOSES, pts, strict=True)]


# --------------------------------------------------------------------------- statistics
def test_signed_turns_geometry():
    xy = np.array([[1.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0]])
    th = np.degrees(signed_turns(xy))
    assert th == pytest.approx([0.0, 90.0, 90.0])  # straight, left turn, left turn


def test_straight_path_has_no_bend():
    ss, _lib, D = _world({"synthetic_straight": _straight(_dirs())})
    sigs = [s for s in ss.signatures if s.compound_id == "synthetic_straight"]
    mask = np.isin(FEATS, list(STRESS))
    dp = dose_path(sigs, PathAxes(D[0], mask, scale=(1.0, 0.05)))
    assert dp.bend_index is None and dp.bend_p > 0.5
    assert np.all(np.diff(dp.x) > 0)  # moves monotonically along the reference


def test_injected_bend_is_detected_at_the_right_dose():
    D = _dirs()
    ss, _lib, _ = _world({"synthetic_bent": _bent(D)})
    sigs = [s for s in ss.signatures if s.compound_id == "synthetic_bent"]
    mask = np.isin(FEATS, list(STRESS))
    dp = dose_path(sigs, PathAxes(D[0], mask, scale=(1.0, 0.05)))
    assert dp.bend_index == 2 and dp.bend_p < 0.05  # turns at the 2nd dose (100 nM)
    assert abs(np.degrees(dp.turn_angle[1])) > 60


def test_bend_within_noise_is_not_called():
    D = _dirs()
    ss, _lib, _ = _world({"synthetic_noisy": _bent(D, se=0.6)})
    sigs = [s for s in ss.signatures if s.compound_id == "synthetic_noisy"]
    dp = dose_path(sigs, PathAxes(D[0], np.isin(FEATS, list(STRESS)), scale=(1.0, 0.05)))
    assert dp.bend_index is None


def test_bend_angle_is_measured_on_scaled_axes():
    """Angles are taken after dividing each axis by its spread (PathAxes.scale): an axis
    with a huge spread flattens the same turn into a straight line."""
    D = _dirs()
    ss, _, _ = _world({"synthetic_bent": _bent(D)})
    sigs = [s for s in ss.signatures if s.compound_id == "synthetic_bent"]
    mask = np.isin(FEATS, list(STRESS))
    flat = dose_path(sigs, PathAxes(D[0], mask, scale=(1.0, 100.0)))
    scaled = dose_path(sigs, PathAxes(D[0], mask, scale=(1.0, 0.05)))
    assert flat.bend_index is None and scaled.bend_index == 2


# --------------------------------------------------------------------------- plot + guards
def test_plot_caption_and_reference_path():
    D = _dirs()
    ss, lib, _ = _world({"synthetic_bent": _bent(D), "synthetic_straight": _straight(D)})
    res = tc.dose_paths(ss, ["synthetic_bent", "synthetic_straight"], lib, STRESS,
                        focus="synthetic_bent")  # fmt: skip
    assert res.meta["refs"]["synthetic_bent"] == "A"
    assert "bent tracks A up to 100 nM; above 100 nM it moves toward stress" in res.caption
    assert any(m.end_label == "A (class mean)" for m in res.spec.marks(LineMark))
    assert any(isinstance(m, EllipseMark) for m in res.spec.marks())
    assert res.classic is not None and res.classic.plot == "dose_paths_classic"
    assert "Bonferroni" in res.methods


def test_non_significant_bend_has_no_annotation():
    D = _dirs()
    ss, lib, _ = _world({"synthetic_straight": _straight(D), "synthetic_noisy": _bent(D, 0.6)})
    res = tc.dose_paths(ss, ["synthetic_straight", "synthetic_noisy"], lib, STRESS)
    assert not [m for m in res.spec.marks(TextMark) if "bend" in m.text]
    assert "without a significant bend" in res.caption


def test_too_few_doses_draws_points_not_a_path():
    D = _dirs()
    two = [(1e-7, 3 * D[0], 0.02), (1e-6, 6 * D[0], 0.02)]
    ss, lib, _ = _world({"synthetic_two": two, "synthetic_straight": _straight(D)})
    res = tc.dose_paths(ss, ["synthetic_two", "synthetic_straight"], lib, STRESS,
                        focus="synthetic_two")  # fmt: skip
    assert "too few doses for a path" in res.caption
    assert any(isinstance(m, NoticeMark) for m in res.spec.marks())
    assert not any(m.end_label == "two" for m in res.spec.marks(LineMark))


def test_more_than_eight_compounds_become_small_multiples():
    D = _dirs()
    q = {f"synthetic_q{i}": _straight(D) for i in range(9)}
    ss, lib, _ = _world(q)
    res = tc.dose_paths(ss, list(q), lib, STRESS, n_boot=100)
    assert res.meta["small_multiples"] and len(res.spec.panels) == 9
    assert "small multiples" in res.caption
    eight = tc.dose_paths(ss, list(q)[:8], lib, STRESS, n_boot=100)
    assert not eight.meta["small_multiples"] and len(eight.spec.panels) == 1


def test_gate_and_errors():
    D = _dirs()
    ss, lib, _ = _world({"synthetic_bent": _bent(D), "synthetic_straight": _straight(D)})
    res = tc.dose_paths(ss, ["synthetic_bent", "synthetic_straight"], lib, STRESS,
                        gate={"synthetic_bent"})  # fmt: skip
    assert res.meta["dropped"] == ["synthetic_straight"]
    assert "failed the retrievability gate" in res.caption
    with pytest.raises(TinctureDataError, match="stress genes"):
        tc.dose_paths(ss, ["synthetic_bent"], lib, {"nope"})
    with pytest.raises(TinctureDataError, match="Focus"):
        tc.dose_paths(ss, ["synthetic_bent"], lib, STRESS, focus="synthetic_straight")
    with pytest.raises(TinctureDataError, match="cell_loss"):
        tc.dose_paths(ss, ["synthetic_bent"], lib, y_axis="cell_loss")


def test_cell_loss_axis():
    D = _dirs()
    ss, lib, _ = _world({"synthetic_straight": _straight(D)})
    loss = {(s.compound_id, s.dose): 0.1 * k for s in ss.signatures
            for k, d in enumerate(DOSES) if s.dose == d}  # fmt: skip
    res = tc.dose_paths(ss, ["synthetic_straight"], lib, y_axis="cell_loss", cell_loss=loss)
    assert res.meta["paths"]["synthetic_straight"].y == pytest.approx([0, 0.1, 0.2, 0.3])
    assert "cell loss" in res.methods


@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_dose_paths_baseline():
    D = _dirs()
    ss, lib, _ = _world({"synthetic_bent": _bent(D), "synthetic_straight": _straight(D)})
    fig = tc.dose_paths(ss, ["synthetic_bent", "synthetic_straight"], lib, STRESS,
                        focus="synthetic_bent").figure()  # fmt: skip
    check_colorblind(fig)
    return fig
