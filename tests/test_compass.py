"""Compass: NNLS statistics, noise correction, leave-one-out, guards, encoding, baseline."""

from __future__ import annotations

import numpy as np
import pytest

import tincture as tc
from tincture.errors import TinctureDataError
from tincture.model.reference import library_from_signatures
from tincture.model.signature import Signature
from tincture.plots.compass import class_coherence, select_spokes, short_name
from tincture.spec import NoticeMark, PointMark, TextMark
from tincture.stats.projection import compass_xy, nnls_project, seriate, spoke_angles
from tincture.theme import check_colorblind

P = 400
FEATS = np.array([f"g{i}" for i in range(P)])


def _dirs(seed: int = 0, k: int = 4) -> np.ndarray:
    rng = np.random.default_rng(seed)
    D = rng.normal(size=(k, P))
    q, _ = np.linalg.qr(D.T)  # orthonormal directions
    return q.T[:k]


def _members(
    direction: np.ndarray, n: int, scale: float, rng, noise: float = 0.05
) -> list[np.ndarray]:
    return [scale * direction + rng.normal(0, noise, P) for _ in range(n)]


def _sig(
    cid: str, effect: np.ndarray, se: float | None = 0.05, name: str | None = None
) -> Signature:
    return Signature(
        features=FEATS,
        effect=effect,
        se=None if se is None else np.full(P, se),
        compound_id=cid,
        meta={"compound_name": name or cid},
    )


def _library(seed: int = 0, correlated: bool = False):
    rng = np.random.default_rng(seed)
    D = _dirs(seed)
    if correlated:  # make class C nearly parallel to class A
        D[2] = 0.95 * D[0] + np.sqrt(1 - 0.95**2) * D[2]
    sigs, classes = {}, {}
    for name, d in zip(["A", "B", "C"], D[:3], strict=True):
        ms = _members(d, 4, 10.0, rng)
        ids = [f"synthetic_{name}{i}" for i in range(4)]
        sigs.update(dict(zip(ids, ms, strict=True)))
        classes[name] = ids
    lib = library_from_signatures(
        sigs, FEATS, classes, version="synthetic v1", source="synthetic refs"
    )
    return lib, D, sigs


# --------------------------------------------------------------------------- statistics
def test_nnls_recovers_mixture_and_residual():
    D = _dirs()
    y = 0.7 * D[0] + 0.3 * D[1]
    p = nnls_project(y, D[:3])
    assert p.shares[:2] == pytest.approx([0.7, 0.3], abs=1e-6)
    assert p.residual == pytest.approx(0.0, abs=1e-9)
    ortho = nnls_project(D[3], D[:3])
    assert ortho.residual == pytest.approx(1.0, abs=1e-9)


def test_noise_correction_and_signal_test():
    rng = np.random.default_rng(1)
    D = _dirs()
    noise = rng.normal(0, 1.0, P)
    pure = nnls_project(noise, D[:3], se=np.ones(P))
    assert not pure.has_signal and pure.magnitude < 0.15 * np.sqrt(P)  # raw norm would be ~20
    sig = 30 * D[0] + noise
    p = nnls_project(sig, D[:3], se=np.ones(P))
    assert p.has_signal and p.magnitude == pytest.approx(30, rel=0.15)
    assert p.residual < 0.1  # noise no longer counted as "unexplained"
    uncorrected = nnls_project(sig, D[:3])
    assert uncorrected.residual > 0.2 > 2 * p.residual  # noise inflates the naive residual


def test_radius_is_linear_in_magnitude_and_direction_from_weights():
    D = _dirs()
    ang = spoke_angles(3)
    for scale in (1.0, 2.0, 5.0):
        x, y, conc = compass_xy(nnls_project(scale * D[1], D[:3]), ang)
        assert np.hypot(x, y) == pytest.approx(scale)
        assert np.arctan2(y, x) == pytest.approx(ang[1])
        assert conc == pytest.approx(1.0)


def test_seriation_puts_similar_references_together():
    sim = np.array([[1, 0.1, 0.9, 0.1], [0.1, 1, 0.1, 0.8], [0.9, 0.1, 1, 0.1], [0.1, 0.8, 0.1, 1]])
    order = list(seriate(sim))
    pos = {k: order.index(k) for k in range(4)}
    assert abs(pos[0] - pos[2]) == 1 and abs(pos[1] - pos[3]) == 1


# --------------------------------------------------------------------------- library
def test_leave_one_out():
    lib, _, sigs = _library()
    loo = lib.without("synthetic_A0")
    assert "synthetic_A0" not in loo.member_ids["A"] and len(loo.member_ids["A"]) == 3
    assert lib.reference_of("synthetic_A0") == ["A"]
    # the compass must project members against the leave-one-out library
    s = {"synthetic_A0": _sig("synthetic_A0", sigs["synthetic_A0"])}
    res = tc.compass(s, lib, focus="synthetic_A0", max_spokes=3)
    p = res.meta["points"][0]
    assert p.projection.shares[res.meta["spokes"].index("A")] > 0.9


def test_incoherent_classes_are_not_spokes():
    lib, _, sigs = _library()
    rng = np.random.default_rng(3)
    junk = {f"synthetic_J{i}": rng.normal(0, 10, P) for i in range(4)}  # members disagree
    sigs2 = {**sigs, **junk}
    classes = {n: lib.member_ids[n] for n in lib.names} | {"Junk": list(junk)}
    lib2 = library_from_signatures(sigs2, FEATS, classes, version="v")
    assert class_coherence(lib2, "Junk")[0] < 0.1
    kept, left = select_spokes(lib2)
    assert "Junk" not in kept.names and left == ["Junk"]


def test_too_few_coherent_classes_raise():
    rng = np.random.default_rng(4)
    sigs = {f"synthetic_{c}{i}": rng.normal(0, 1, P) for c in "XY" for i in range(3)}
    lib = library_from_signatures(
        sigs, FEATS, {c: [f"synthetic_{c}{i}" for i in range(3)] for c in "XY"}, version="v"
    )
    with pytest.raises(TinctureDataError, match="agree"):
        select_spokes(lib)


# --------------------------------------------------------------------------- guards
def _queries(D, rng):
    return {
        "synthetic_q_mix": _sig("synthetic_q_mix", 9 * D[0] + 4 * D[1] + rng.normal(0, 0.05, P)),
        "synthetic_q_novel": _sig("synthetic_q_novel", 12 * D[3] + rng.normal(0, 0.05, P)),
        "synthetic_q_weak": _sig("synthetic_q_weak", rng.normal(0, 0.05, P)),
        "synthetic_q_B": _sig("synthetic_q_B", 8 * D[1] + rng.normal(0, 0.05, P)),
    }


def test_footer_names_reference_set_and_version():
    lib, D, _ = _library()
    res = tc.compass(_queries(D, np.random.default_rng(5)), lib, focus="synthetic_q_mix")
    assert any("synthetic refs synthetic v1" in f for f in res.spec.footer)


def test_poorly_explained_flag_only_with_signal():
    lib, D, _ = _library()
    res = tc.compass(_queries(D, np.random.default_rng(5)), lib, focus="synthetic_q_mix")
    assert res.meta["novel"] == ["synthetic_q_novel"]
    assert "synthetic_q_weak" in res.meta["weak"]  # weak: never called novel
    texts = [m for m in res.spec.marks(TextMark) if "poorly explained" in m.text]
    assert texts and texts[0].color == "novel"
    assert "1 with enough signal are poorly explained" in res.caption


def test_correlated_spokes_guard():
    lib, D, _ = _library(correlated=True)
    res = tc.compass(_queries(D, np.random.default_rng(5)), lib, focus="synthetic_q_B")
    assert res.meta["correlated"], "expected a correlated neighbouring pair"
    assert "not independent" in res.caption
    assert any(isinstance(m, NoticeMark) for m in res.spec.marks())
    assert any(isinstance(m, TextMark) and m.text.startswith("similar (") for m in res.spec.marks())


def test_caption_template_and_residual():
    lib, D, _ = _library()
    res = tc.compass(_queries(D, np.random.default_rng(5)), lib, focus="synthetic_q_mix")
    assert "synthetic_q_mix most resembles A (weight 0.6" in res.caption  # 9/(9+4) = 0.69
    assert "unexplained by the 3 references" in res.caption


def test_outlier_drawn_on_rim_with_magnitude():
    lib, D, _ = _library()
    q = _queries(D, np.random.default_rng(5))
    for i in range(12):  # many ordinary compounds so the 95th percentile ignores the outlier
        q[f"synthetic_n{i}"] = _sig(f"synthetic_n{i}", (6 + i * 0.2) * D[i % 3])
    q["synthetic_big"] = _sig("synthetic_big", 100 * D[2], name="Big")
    res = tc.compass(q, lib, focus="synthetic_q_mix", label_top=20)
    assert "synthetic_big" in res.meta["clipped"]
    assert any(m.text.startswith("Big (100)") for m in res.spec.marks(TextMark))


def test_gate_drop_and_grey():
    lib, D, _ = _library()
    q = _queries(D, np.random.default_rng(5))
    gate = {"synthetic_q_mix", "synthetic_q_B"}
    res = tc.compass(q, lib, focus="synthetic_q_mix", gate=gate)
    assert "2 compound(s) failed the retrievability gate and are not shown" in res.caption
    grey = tc.compass(q, lib, focus="synthetic_q_mix", gate=gate, gate_mode="grey")
    assert any(isinstance(m, PointMark) and np.all(m.hollow) for m in grey.spec.marks(PointMark))


def test_unknown_focus_and_classic():
    lib, D, _ = _library()
    q = _queries(D, np.random.default_rng(5))
    with pytest.raises(TinctureDataError, match="Focus"):
        tc.compass(q, lib, focus="nope")
    res = tc.compass(q, lib, focus="synthetic_q_mix")
    assert res.classic is not None and res.classic.plot == "compass_classic"
    one = tc.compass({"synthetic_q_mix": q["synthetic_q_mix"]}, lib, focus="synthetic_q_mix")
    assert one.classic is None  # PCA needs ≥ 3 signatures


@pytest.mark.parametrize(
    "raw,clean",
    [
        ("Abexinostat (PCI-24781)", "Abexinostat"),
        ("(+)-JQ1", "(+)-JQ1"),
        ("Quisinostat (JNJ-26481585) 2HCl", "Quisinostat"),
        ("Toremifene Citrate", "Toremifene"),
        ("ENMD-2076 L-(+)-Tartaric acid", "ENMD-2076"),
        ("SAHA", "SAHA"),
        ("Lapatinib (GW-572016) Ditosylate", "Lapatinib"),
    ],
)
def test_short_name(raw, clean):
    assert short_name(raw) == clean


@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_compass_baseline():
    lib, D, _ = _library()
    fig = tc.compass(_queries(D, np.random.default_rng(5)), lib, focus="synthetic_q_mix").figure()
    check_colorblind(fig)
    return fig


def test_compass_budget_500_compounds():
    """HANDOVER budget: 500 compounds in under 3 s (fit + spec + figure)."""
    import time

    lib, D, _ = _library()
    rng = np.random.default_rng(11)
    q = {
        f"synthetic_b{i}": _sig(f"synthetic_b{i}", (4 + rng.random() * 6) * D[i % 3]
                                + rng.normal(0, 0.05, P))
        for i in range(500)
    }  # fmt: skip
    t0 = time.perf_counter()
    res = tc.compass(q, lib, focus="synthetic_b0")
    res.figure()
    elapsed = time.perf_counter() - t0
    assert len(res.meta["points"]) == 500
    assert elapsed < 3.0, f"500 compounds took {elapsed:.1f} s (budget 3 s)"
