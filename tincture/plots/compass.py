"""Compass: what known mechanism does this compound resemble, and how strongly? (HANDOVER §5.1)

Replaces PCA. Spokes are reference mechanisms (not principal components), ordered so that
similar references sit next to each other. For each compound, a non-negative least-squares
fit on the references gives the **direction** (weighted towards the spokes it resembles);
the **distance from centre is the signature's magnitude** (linear, not area-scaled).
A dashed halo shows how much of the signature the references do *not* explain (residual);
compounds with residual > 0.5 are purple and labelled "poorly explained" (the novelty flag).

Guards (each tested):
* correlated references: an outer ring shows similarity between neighbouring spokes and any
  pair with cosine > 0.8 is called out as "not independent";
* the reference set and version are printed in the footer;
* the residual is always computed and shown; > 0.5 is labelled;
* leave-one-out: a compound is never explained by a reference built from itself.
Classic fallback: PCA of the same signatures.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np

from tincture.errors import TinctureDataError
from tincture.model.reference import ReferenceLibrary
from tincture.model.signature import Signature
from tincture.report.captions import TEMPLATES, Caption
from tincture.report.classic import build_classic, register_classic
from tincture.report.methods import Methods, register_block
from tincture.result import PlotResult
from tincture.spec import (
    EllipseMark,
    LineMark,
    NoticeMark,
    Panel,
    PlotSpec,
    PointMark,
    Scale,
    TextMark,
)
from tincture.stats.projection import (
    Projection,
    bootstrap_ellipse,
    compass_xy,
    nnls_project,
    seriate,
    spoke_angles,
)

CORRELATED = 0.8
POORLY_EXPLAINED = 0.5
MAX_SPOKES = 8

register_block(
    "compass",
    "Each signature was fitted as a non-negative least-squares combination of {n_refs} "
    "reference signatures (scaled to unit length), giving per-reference weights and the "
    "reconstruction R²; residual = 1 − R². On the compass, references are spokes ordered by "
    "optimal-leaf-ordered clustering of their cosine similarities; a compound's direction is "
    "the weighted resultant of its spokes and its distance from the centre is the norm of its "
    "shrunk signature (linear scale). Compounds that belong to a reference class were fitted "
    "against references rebuilt without them. Ellipses are 95% regions from {n_boot} "
    "posterior draws of the signature.",
)


@dataclass
class CompassPoint:
    compound_id: str
    label: str
    x: float
    y: float
    projection: Projection
    concentration: float
    gated_out: bool = False
    ellipse: tuple[float, float, float, float, float] | None = None

    @property
    def residual(self) -> float:
        return self.projection.residual


def select_spokes(
    lib: ReferenceLibrary, max_spokes: int = MAX_SPOKES
) -> tuple[ReferenceLibrary, list[str]]:
    """Keep the ``max_spokes`` strongest references (by class-mean magnitude), then order
    them by similarity. Returns the ordered library and the names left out."""
    mags = np.linalg.norm(lib.matrix, axis=1)
    keep = np.argsort(-mags)[:max_spokes]
    sub = lib.reorder(sorted(keep))
    order = seriate(sub.similarity())
    left_out = [lib.names[i] for i in range(len(lib)) if i not in set(keep)]
    return sub.reorder([int(i) for i in order]), left_out


def _nice_step(rmax: float) -> float:
    raw = rmax / 3
    mag = 10 ** np.floor(np.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if m * mag >= raw:
            return float(m * mag)
    return float(10 * mag)


def short_name(name: str) -> str:
    """Drop salt forms and code synonyms in parentheses: 'Abexinostat (PCI-24781)' →
    'Abexinostat'. Names that are only a code are kept."""
    import re

    base = re.sub(r"\s*\(.*?\)", "", name).strip()
    base = re.sub(r"\s+(2?HCl|HCl|hydrochloride|mesylate|sodium salt)$", "", base, flags=re.I)
    return base or name


def _raw(sig: Signature, feats: np.ndarray) -> tuple[np.ndarray, np.ndarray | None]:
    """Unshrunk effect and SE when available (noise correction needs unbiased effects),
    else the stored effect and SE."""
    a = sig.align(feats)
    raw, raw_se = sig.meta.get("raw_effect"), sig.meta.get("raw_se")
    if raw is not None and raw_se is not None:
        pos = {f: i for i, f in enumerate(sig.features)}
        idx = np.array([pos.get(f, -1) for f in feats])
        take = np.where(idx >= 0, idx, 0)
        y = np.where(idx >= 0, np.asarray(raw, float)[take], np.nan)
        se = np.where(idx >= 0, np.asarray(raw_se, float)[take], np.nan)
        return y, se
    return a.effect, a.se


def compute_points(
    signatures: Mapping[str, Signature],
    lib: ReferenceLibrary,
    angles: np.ndarray,
    gate: set[str] | None = None,
    ellipses_for: Sequence[str] = (),
    n_boot: int = 200,
) -> list[CompassPoint]:
    pts = []
    feats = lib.features
    for cid, sig in signatures.items():
        y, se = _raw(sig, feats)
        own = lib.reference_of(cid)
        use = lib.without(cid) if own else lib
        # map weights back onto the full spoke set (a dropped reference gets weight 0)
        R = use.matrix
        proj = nnls_project(y, R, se)
        if own:
            w = np.zeros(len(lib))
            for j, n in enumerate(use.names):
                w[lib.names.index(n)] = proj.weights[j]
            proj = Projection(
                weights=w, r2=proj.r2, magnitude=proj.magnitude, signal_z=proj.signal_z
            )
        x, y_, conc = compass_xy(proj, angles)
        ell = None
        if cid in ellipses_for and se is not None:
            sub_angles = np.array([angles[lib.names.index(n)] for n in use.names])
            ell = bootstrap_ellipse(
                np.nan_to_num(y), np.nan_to_num(se), R, sub_angles, n_boot=n_boot
            )
        label = short_name(str(sig.meta.get("compound_name") or cid))
        pts.append(
            CompassPoint(
                cid,
                label,
                x,
                y_,
                proj,
                conc,
                gated_out=gate is not None and cid not in gate,
                ellipse=ell,
            )
        )
    return pts


def build_spec(
    points: list[CompassPoint],
    lib: ReferenceLibrary,
    angles: np.ndarray,
    focus: str | None = None,
    label_top: int = 4,
    spoke_labels: Mapping[str, str] | None = None,
    left_out: Sequence[str] = (),
    title: str = "Compass",
    n_boot: int = 200,
    footer: Sequence[str] = (),
) -> PlotSpec:
    if not points:
        raise TinctureDataError("No signatures to place on the compass.")
    spoke_labels = spoke_labels or {}
    by_id = {p.compound_id: p for p in points}
    if focus is not None and focus not in by_id:
        raise TinctureDataError(
            f"Focus {focus!r} is not among the signatures.", fix=f"Use one of {sorted(by_id)[:5]} …"
        )
    rmax_data = max(np.hypot(p.x, p.y) for p in points)
    step = _nice_step(rmax_data)
    rings = [step * k for k in range(1, int(np.ceil(rmax_data / step)) + 1)]
    R = rings[-1]
    panel = Panel(
        key="compass",
        x=Scale("linear", visible=False, domain=(-1.55 * R, 1.55 * R)),
        y=Scale("linear", visible=False, domain=(-1.32 * R, 1.32 * R)),
        aspect="equal",
        frame=False,
        title="",
    )
    t = np.linspace(0, 2 * np.pi, 241)
    # rings (magnitude grid) + their labels along the upward spoke gap
    for i, r in enumerate(rings):
        panel.add(
            LineMark(
                x=r * np.cos(t),
                y=r * np.sin(t),
                color="grid",
                width=0.9 if i == len(rings) - 1 else 0.6,
                style="solid" if i == len(rings) - 1 else "dashed",
                z=0.5,
            )
        )
        a = np.pi / 2 - np.pi / len(lib)  # between the first two spokes
        panel.add(
            TextMark(
                x=r * np.cos(a),
                y=r * np.sin(a),
                text=f"{r:g}",
                color="ink_soft",
                ha="center",
                va="center",
                size=None,
                z=0.6,
            )
        )
    # outer similarity ring between neighbouring spokes
    sim = lib.similarity()
    n = len(lib)
    correlated = []
    ring_r = 1.1 * R
    for k in range(n):
        j = (k + 1) % n
        if n < 3 and j == 0:
            continue
        a0, a1 = angles[k], angles[j]
        if a1 > a0:
            a1 -= 2 * np.pi
        arc = np.linspace(a0, a1, 40)
        s = float(sim[k, j])
        strong = s > CORRELATED
        if strong:
            correlated.append((lib.names[k], lib.names[j], s))
        panel.add(
            LineMark(
                x=ring_r * np.cos(arc),
                y=ring_r * np.sin(arc),
                color="ink_soft" if strong else "grid",
                width=0.4 + 2.2 * max(s, 0),
                z=0.4,
            )
        )
        if strong:
            am = (a0 + a1) / 2
            panel.add(
                TextMark(
                    x=1.17 * R * np.cos(am),
                    y=1.17 * R * np.sin(am),
                    text=f"similar ({s:.2f})",
                    color="ink_soft",
                    ha="left" if np.cos(am) > 0.2 else ("right" if np.cos(am) < -0.2 else "center"),
                    va="center",
                )
            )
    # spokes
    for k, name in enumerate(lib.names):
        ca, sa = np.cos(angles[k]), np.sin(angles[k])
        panel.add(LineMark(x=[0, R * ca], y=[0, R * sa], color="neutral", width=0.8, z=0.7))
        panel.add(
            PointMark(x=[R * ca], y=[R * sa], color="reference", size=22, categorical=True, z=3)
        )
        ha: Literal["left", "center", "right"] = (
            "left" if ca > 0.25 else ("right" if ca < -0.25 else "center")
        )
        va: Literal["top", "center", "bottom"] = (
            "bottom" if sa > 0.25 else ("top" if sa < -0.25 else "center")
        )
        lab_r = 1.25 * R if ha == "center" else 1.24 * R
        panel.add(
            TextMark(
                x=lab_r * ca,
                y=lab_r * sa,
                text=spoke_labels.get(name, name),
                color="reference",
                ha=ha,
                va=va,
                categorical=True,
            )
        )
    panel.add(PointMark(x=[0], y=[0], color="ink", size=10, z=3))
    # compounds
    gated = [p for p in points if p.gated_out]
    active = [p for p in points if not p.gated_out]
    novel = [
        p
        for p in active
        if p.projection.has_signal and p.residual > POORLY_EXPLAINED and p.compound_id != focus
    ]
    weak = [p for p in active if not p.projection.has_signal and p.compound_id != focus]
    plain = [p for p in active if p not in novel and p not in weak and p.compound_id != focus]
    if gated:
        panel.add(
            PointMark(
                x=[p.x for p in gated],
                y=[p.y for p in gated],
                color="neutral",
                size=8,
                hollow=True,
                z=1.5,
            )
        )
    if weak:
        panel.add(
            PointMark(x=[p.x for p in weak], y=[p.y for p in weak], color="neutral", size=8, z=1.8)
        )
    if plain:
        panel.add(
            PointMark(
                x=[p.x for p in plain], y=[p.y for p in plain], color="mid_grey", size=14, z=2
            )
        )
    halo_scale = 0.22 * R
    for p in novel:
        rr = p.residual * halo_scale
        panel.add(
            EllipseMark(
                cx=[p.x],
                cy=[p.y],
                width=[2 * rr],
                height=[2 * rr],
                color="novel",
                linestyle="dashed",
                linewidth=0.8,
                z=2.2,
                categorical=True,
            )
        )
    if novel:
        panel.add(
            PointMark(
                x=[p.x for p in novel],
                y=[p.y for p in novel],
                color="novel",
                size=16,
                z=2.5,
                categorical=True,
            )
        )
    # labels: focus, novel, and the strongest few others
    to_label = sorted(plain, key=lambda p: -np.hypot(p.x, p.y))[:label_top]
    for p in to_label:
        panel.add(
            TextMark(x=p.x, y=p.y, text=f"  {p.label}", color="ink_soft", ha="left", va="center")
        )
    for p in sorted(novel, key=lambda p: -p.residual)[:label_top]:
        panel.add(
            TextMark(
                x=p.x,
                y=p.y - p.residual * halo_scale,
                text=f"{p.label} · poorly explained",
                color="novel",
                ha="center",
                va="top",
                categorical=True,
            )
        )
    fp = by_id.get(focus) if focus else None
    if fp is not None:
        if fp.ellipse is not None:
            cx, cy, w, h, ang = fp.ellipse
            panel.add(
                EllipseMark(
                    cx=[cx],
                    cy=[cy],
                    width=[w],
                    height=[h],
                    angle=[ang],
                    color="focus",
                    fill=True,
                    alpha=0.25,
                    linewidth=0.6,
                    z=2.6,
                    categorical=True,
                )
            )
        panel.add(
            LineMark(x=[0, fp.x], y=[0, fp.y], color="focus", style="dashed", width=1.0, z=2.7)
        )
        if fp.residual > POORLY_EXPLAINED and fp.projection.has_signal:
            rr = fp.residual * halo_scale
            panel.add(
                EllipseMark(
                    cx=[fp.x],
                    cy=[fp.y],
                    width=[2 * rr],
                    height=[2 * rr],
                    color="novel",
                    linestyle="dashed",
                    linewidth=0.8,
                    z=2.6,
                    categorical=True,
                )
            )
        panel.add(PointMark(x=[fp.x], y=[fp.y], color="focus", size=34, z=4, categorical=True))
        side: Literal["left", "right"] = "left" if fp.x >= 0 else "right"
        panel.add(
            TextMark(
                x=fp.x,
                y=fp.y,
                text=f"  {fp.label}  ",
                color="focus",
                ha=side,
                va="bottom",
                categorical=True,
                z=5,
            )
        )

    # caption
    cap = Caption()
    if fp is not None:
        sh = fp.projection.shares
        top_k = int(np.argmax(sh)) if sh.sum() > 0 else None
        cap.add(
            TEMPLATES["compass"],
            compound=fp.label,
            top_ref=spoke_labels.get(lib.names[top_k], lib.names[top_k])
            if top_k is not None
            else None,
            w=float(sh[top_k]) if top_k is not None else None,
            pct_unexplained=fp.residual,
            n_refs=len(lib),
        )
        if fp.concentration < 0.6 and fp.projection.weights.sum() > 0:
            cap.add(
                "Its weights are spread across distant references (concentration {c:.2f}), so its "
                "direction is a mixture, not a single mechanism.",
                c=fp.concentration,
            )
    if fp is not None and not fp.projection.has_signal:
        cap.weak(
            f"{fp.label} shows no clear signal (z = {fp.projection.signal_z:.1f}); its "
            "position and residual are not interpretable."
        )
    cap.add(
        "{n} compounds shown; distance from centre = noise-corrected signature magnitude "
        "(log-fold changes vs vehicle). {k} with a clear signal are poorly explained by the "
        "references (residual > {thr}); {w} are too weak to judge (light grey).",
        n=len(active),
        k=len(novel)
        + int(fp is not None and fp.projection.has_signal and fp.residual > POORLY_EXPLAINED),
        thr=POORLY_EXPLAINED,
        w=len(weak),
    )
    for ref_a, ref_b, cos_ab in correlated:
        cap.weak(
            f"the {spoke_labels.get(ref_a, ref_a)} and {spoke_labels.get(ref_b, ref_b)} "
            f"references are similar (cosine {cos_ab:.2f}): these spokes are not independent."
        )
    if correlated:
        panel.add(
            NoticeMark(
                text="some neighbouring spokes are not independent (outer ring)",
                loc="lower left",
                level="weak",
            )
        )
    if gated:
        cap.weak(
            f"{len(gated)} compound(s) failed the retrievability gate and are drawn as hollow grey."
        )
    meth = Methods().use("compass", n_refs=len(lib), n_boot=n_boot)
    feet = [f"References: {lib.label()} ({len(lib)} mechanism classes; class means)"]
    if left_out:
        feet.append(
            "Weaker classes not used as spokes: "
            + ", ".join(spoke_labels.get(x, x) for x in left_out)
        )
    feet += list(footer)
    return PlotSpec(
        plot="compass",
        panels=[panel],
        title=title,
        caption=cap.render(),
        methods=meth.render(),
        footer=feet,
        warnings=list(cap.notes),
        aspect=0.8,
        meta={
            "points": points,
            "correlated": correlated,
            "rings": rings,
            "spokes": list(lib.names),
            "angles": angles,
        },
    ).validate()


@register_classic("compass", "PCA")
def build_classic_spec(
    signatures: Mapping[str, Signature] | None = None,
    features: np.ndarray | None = None,
    focus: str | None = None,
    **_: object,
) -> PlotSpec:
    if not signatures or features is None:
        raise TinctureDataError("PCA fallback needs signatures.")
    ids = list(signatures)
    X = np.vstack([np.nan_to_num(signatures[c].align(features).effect) for c in ids])
    Xc = X - X.mean(axis=0)
    U, S, _vt = np.linalg.svd(Xc, full_matrices=False)
    pcs = U[:, :2] * S[:2]
    ve = S**2 / (S**2).sum()
    p = Panel(x=Scale(label=f"PC1 ({ve[0]:.0%})"), y=Scale(label=f"PC2 ({ve[1]:.0%})"))
    mask = np.array([c == focus for c in ids])
    p.add(PointMark(x=pcs[~mask, 0], y=pcs[~mask, 1], color="mid_grey", size=10))
    if mask.any():
        lab = str(signatures[focus].meta.get("compound_name") or focus) if focus else ""
        p.add(
            PointMark(
                x=pcs[mask, 0],
                y=pcs[mask, 1],
                color="focus",
                size=24,
                categorical=True,
                labels=[lab],
            )
        )
    return PlotSpec(
        plot="compass_classic",
        panels=[p],
        title="PCA of signatures",
        caption="PCA of the same signatures (axes carry no mechanistic meaning).",
        aspect=0.8,
    ).validate()


def plot(
    signatures: Mapping[str, Signature],
    library: ReferenceLibrary,
    focus: str | None = None,
    gate: set[str] | None = None,
    spoke_labels: Mapping[str, str] | None = None,
    max_spokes: int = MAX_SPOKES,
    label_top: int = 4,
    n_boot: int = 200,
    title: str = "Compass",
    footer: Sequence[str] = (),
) -> PlotResult:
    """Compass plot.

    ``signatures``: compound id → :class:`Signature` (one per compound, e.g. top dose).
    ``gate``: compound ids that passed retrievability; others are drawn hollow grey.
    """
    lib, left_out = select_spokes(library, max_spokes)
    angles = spoke_angles(len(lib))
    pts = compute_points(
        signatures, lib, angles, gate, ellipses_for=[focus] if focus else [], n_boot=n_boot
    )
    spec = build_spec(
        pts, lib, angles, focus, label_top, spoke_labels, left_out, title, n_boot, footer
    )
    classic = build_classic("compass", signatures=signatures, features=lib.features, focus=focus)
    return PlotResult(spec=spec, classic=classic)
