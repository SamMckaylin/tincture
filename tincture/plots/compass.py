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
    LegendSpec,
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
MIN_COHERENCE = 0.2
#: Two-member classes must agree more strongly to be trusted as a spoke.
MIN_COHERENCE_PAIR = 0.4
#: Novelty is judged only for compounds with at least this share of the rim magnitude.
NOVELTY_MIN_MAGNITUDE = 0.25
POORLY_EXPLAINED = 0.5
MAX_SPOKES = 8

register_block(
    "compass",
    "Each signature was fitted as a non-negative least-squares combination of {n_refs} "
    "reference signatures (class means of member compounds, scaled to unit length), giving "
    "per-reference weights. Spokes are the most coherent classes (mean pairwise cosine of "
    "members × √members), ordered by optimal-leaf-ordered clustering of their cosine "
    "similarities. A compound's direction is the weighted resultant of its spokes; its "
    "distance from the centre is its noise-corrected magnitude √(‖y‖² − Σ se²) on unshrunk "
    "log fold changes (linear scale), and the share of that signal explained by the "
    "references is R²; residual = 1 − R². Residual is only judged for compounds with signal "
    "(χ² test over features, z > 3). Compounds that belong to a reference class were fitted "
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
    """Keep the ``max_spokes`` most *coherent* references, then order them by similarity.

    Coherence = mean pairwise cosine between a class's member signatures (do its members
    agree?) × √(members). Magnitude alone would favour classes with one or two very strong
    compounds and drop well-defined classes with moderate effects (e.g. Aurora kinase
    inhibitors in sci-Plex), which would then wrongly look "novel". Classes whose members do
    not agree are never used: mean cosine ≥ ``MIN_COHERENCE`` with ≥ 3 members, or
    ≥ ``MIN_COHERENCE_PAIR`` for a two-member class.
    Returns the ordered library and the names left out.
    """
    scores = np.array([class_coherence(lib, n) for n in lib.names])
    ok = np.where(
        scores[:, 1] >= 3, scores[:, 0] >= MIN_COHERENCE, scores[:, 0] >= MIN_COHERENCE_PAIR
    )
    rank = np.argsort(-(scores[:, 0] * np.sqrt(scores[:, 1])))
    keep = [int(i) for i in rank if ok[i]][:max_spokes]
    if len(keep) < 2:
        raise TinctureDataError(
            "Fewer than two reference classes have members that agree with each other.",
            fix="Check the class annotations or lower MIN_COHERENCE.",
        )
    sub = lib.reorder(sorted(keep))
    order = seriate(sub.similarity())
    left_out = [lib.names[i] for i in range(len(lib)) if i not in set(keep)]
    out = sub.reorder([int(i) for i in order])
    out.meta["coherence"] = {lib.names[i]: float(scores[i, 0]) for i in range(len(lib))}
    return out, left_out


def class_coherence(lib: ReferenceLibrary, name: str) -> tuple[float, int]:
    """(mean pairwise cosine between members, number of members)."""
    M = lib.member_signatures[name]
    n = M.shape[0]
    if n < 2:
        return 0.0, n
    U = M / np.linalg.norm(M, axis=1, keepdims=True)
    S = U @ U.T
    return float(S[np.triu_indices(n, 1)].mean()), n


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

    base = name.strip()
    salts = (
        r"2?HCl|hydrochloride|dihydrochloride|mesylate|sodium salt|citrate|maleate|"
        r"phosphate|sulfate|tosylate|L-\(\+\)-tartaric acid|tartrate"
    )
    while True:  # strip trailing "(synonym)" and salt words, in any order
        new = re.sub(r"\s+\([^()]*\)$", "", base).strip()
        new = re.sub(rf"\s+(?:di|bi|tri)?({salts})$", "", new, flags=re.I).strip()
        if new == base:
            break
        base = new
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


@dataclass
class _Label:
    x: float
    y: float
    text: str
    color: str
    priority: int  # lower = placed first
    forced: bool = False


def place_labels(
    labels: list[_Label], R: float, occupied: list[tuple[float, float, float, float]] | None = None
) -> list[TextMark]:
    """Greedy label placement: try right / left / above / below each point; skip a label
    that would overlap one already placed (the focus label is always placed)."""
    boxes = list(occupied or [])
    h = 0.075 * R
    out: list[TextMark] = []

    def overlaps(b: tuple[float, float, float, float]) -> bool:
        x0, y0, x1, y1 = b
        return any(not (x1 < c0 or x0 > c1 or y1 < d0 or y0 > d1) for c0, d0, c1, d1 in boxes)

    for lab in sorted(labels, key=lambda item: item.priority):
        w = 0.034 * R * (len(lab.text) + 2)
        gap = 0.03 * R
        cands: list[tuple[tuple[float, float, float, float], float, float, str, str]] = [
            (
                (lab.x + gap, lab.y - h / 2, lab.x + gap + w, lab.y + h / 2),
                lab.x + gap,
                lab.y,
                "left",
                "center",
            ),
            (
                (lab.x - gap - w, lab.y - h / 2, lab.x - gap, lab.y + h / 2),
                lab.x - gap,
                lab.y,
                "right",
                "center",
            ),
            (
                (lab.x - w / 2, lab.y + gap, lab.x + w / 2, lab.y + gap + h),
                lab.x,
                lab.y + gap,
                "center",
                "bottom",
            ),
            (
                (lab.x - w / 2, lab.y - gap - h, lab.x + w / 2, lab.y - gap),
                lab.x,
                lab.y - gap,
                "center",
                "top",
            ),
        ]
        chosen = next((c for c in cands if not overlaps(c[0])), None)
        if chosen is None and lab.forced:
            chosen = cands[0]
        if chosen is None:
            continue
        box, tx, ty, ha, va = chosen
        boxes.append(box)
        out.append(
            TextMark(
                x=tx,
                y=ty,
                text=lab.text,
                color=lab.color,
                ha=ha,  # type: ignore[arg-type]
                va=va,  # type: ignore[arg-type]
                categorical=lab.color in ("focus", "novel"),
                z=5,
            )
        )
    return out


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
    gate_mode: str = "drop",
) -> PlotSpec:
    if not points:
        raise TinctureDataError("No signatures to place on the compass.")
    if gate_mode not in ("drop", "grey"):
        raise ValueError("gate_mode must be 'drop' or 'grey'.")
    spoke_labels = spoke_labels or {}
    by_id = {p.compound_id: p for p in points}
    if focus is not None and focus not in by_id:
        raise TinctureDataError(
            f"Focus {focus!r} is not among the signatures.", fix=f"Use one of {sorted(by_id)[:5]} …"
        )
    gated = [p for p in points if p.gated_out and p.compound_id != focus]
    active = [p for p in points if not p.gated_out or p.compound_id == focus]
    if not active:
        raise TinctureDataError(
            "No compound passed the retrievability gate.",
            fix="Show them anyway with gate=None, or check the gate.",
        )

    # --- radial scale: rim at the 95th percentile so one outlier cannot squash the rest
    mags = np.array([p.projection.magnitude for p in active])
    r_data = float(np.quantile(mags, 0.95)) if len(mags) > 10 else float(mags.max())
    r_data = max(r_data, 1e-9)
    step = _nice_step(r_data)
    rings = [step * k for k in range(1, int(np.ceil(r_data / step)) + 1)]
    R = rings[-1]
    clipped = {p.compound_id for p in active if p.projection.magnitude > R}

    def xy(p: CompassPoint) -> tuple[float, float]:
        if p.compound_id in clipped:
            f = R / p.projection.magnitude
            return p.x * f, p.y * f
        return p.x, p.y

    panel = Panel(
        key="compass",
        x=Scale("linear", visible=False, domain=(-1.6 * R, 1.6 * R)),
        y=Scale("linear", visible=False, domain=(-1.36 * R, 1.36 * R)),
        aspect="equal",
        frame=False,
    )
    t = np.linspace(0, 2 * np.pi, 241)
    # ring labels go in the gap between spokes that holds the fewest compounds
    gaps = angles - np.pi / len(lib)
    pt_ang = np.array([np.arctan2(*xy(p)[::-1]) for p in active])

    def crowd(a: float) -> int:
        d = np.angle(np.exp(1j * (pt_ang - a)))
        return int(np.sum(np.abs(d) < np.pi / (2 * len(lib))))

    gap_angle = float(min(gaps, key=crowd))
    for i, r in enumerate(rings):
        last = i == len(rings) - 1
        panel.add(
            LineMark(
                x=r * np.cos(t),
                y=r * np.sin(t),
                color="grid",
                width=0.9 if last else 0.6,
                style="solid" if last else "dashed",
                z=0.5,
            )
        )
        panel.add(
            TextMark(
                x=r * np.cos(gap_angle),
                y=r * np.sin(gap_angle),
                text=f"{r:g}",
                color="ink_soft",
                ha="center",
                va="center",
                z=0.6,
            )
        )

    # --- outer ring: similarity between neighbouring spokes (guard)
    sim = lib.similarity()
    n = len(lib)
    correlated = []
    ring_r = 1.09 * R
    for k in range(n):
        j = (k + 1) % n
        a0, a1 = angles[k], angles[j]
        if a1 > a0:
            a1 -= 2 * np.pi
        arc = np.linspace(a0, a1, 40)
        s_kj = float(sim[k, j])
        strong = s_kj > CORRELATED
        if strong:
            correlated.append((lib.names[k], lib.names[j], s_kj))
        panel.add(
            LineMark(
                x=ring_r * np.cos(arc),
                y=ring_r * np.sin(arc),
                color="ink_soft" if strong else "grid",
                width=0.4 + 2.4 * max(s_kj, 0),
                style="solid",
                z=0.4,
            )
        )
        if strong:
            am = (a0 + a1) / 2
            ca_m = np.cos(am)
            panel.add(
                TextMark(
                    x=1.15 * R * ca_m,
                    y=1.15 * R * np.sin(am),
                    text=f"similar ({s_kj:.2f})",
                    color="ink_soft",
                    ha="left" if ca_m > 0.2 else ("right" if ca_m < -0.2 else "center"),
                    va="center",
                )
            )

    # --- spokes
    spoke_boxes = []
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
        lx, ly = 1.22 * R * ca, 1.22 * R * sa
        text = spoke_labels.get(name, name)
        panel.add(
            TextMark(x=lx, y=ly, text=text, color="reference", ha=ha, va=va, categorical=True)
        )
        w = 0.034 * R * (len(text) + 1)
        x0 = lx if ha == "left" else (lx - w if ha == "right" else lx - w / 2)
        spoke_boxes.append((x0, ly - 0.05 * R, x0 + w, ly + 0.05 * R))
    panel.add(PointMark(x=[0], y=[0], color="ink", size=10, z=3))

    # --- compounds
    def judged(p: CompassPoint) -> bool:
        return p.projection.has_signal and p.projection.magnitude >= NOVELTY_MIN_MAGNITUDE * R

    others = [p for p in active if p.compound_id != focus]
    weak = [p for p in others if not judged(p)]
    novel = [p for p in others if judged(p) and p.residual > POORLY_EXPLAINED]
    plain = [p for p in others if judged(p) and p.residual <= POORLY_EXPLAINED]
    if gated and gate_mode == "grey":
        g = [xy(p) for p in gated]
        panel.add(
            PointMark(
                x=[a for a, _ in g],
                y=[b for _, b in g],
                color="neutral",
                size=8,
                hollow=True,
                z=1.5,
            )
        )
    for group, colour, size, z in ((weak, "neutral", 9, 1.8), (plain, "mid_grey", 16, 2.0)):
        if group:
            g = [xy(p) for p in group]
            panel.add(
                PointMark(x=[a for a, _ in g], y=[b for _, b in g], color=colour, size=size, z=z)
            )
    halo_scale = 0.12 * R
    for p in novel:
        px, py = xy(p)
        rr = p.residual * halo_scale
        panel.add(
            EllipseMark(
                cx=[px],
                cy=[py],
                width=[2 * rr],
                height=[2 * rr],
                color="novel",
                linestyle="dashed",
                linewidth=0.6,
                z=2.2,
                categorical=True,
            )
        )
    if novel:
        g = [xy(p) for p in novel]
        panel.add(
            PointMark(
                x=[a for a, _ in g],
                y=[b for _, b in g],
                color="novel",
                size=18,
                z=2.5,
                categorical=True,
            )
        )

    fp = by_id.get(focus) if focus else None
    labels: list[_Label] = []

    def name_of(p: CompassPoint) -> str:
        return f"{p.label} ({p.projection.magnitude:.0f})" if p.compound_id in clipped else p.label

    if fp is not None:
        fx, fy = xy(fp)
        tight = fp.ellipse is not None and max(fp.ellipse[2], fp.ellipse[3]) < 0.04 * R
        if fp.ellipse is not None and not tight:
            cx, cy, w_, h_, ang = fp.ellipse
            if fp.compound_id in clipped:
                f = R / fp.projection.magnitude
                cx, cy, w_, h_ = cx * f, cy * f, w_ * f, h_ * f
            panel.add(
                EllipseMark(
                    cx=[cx],
                    cy=[cy],
                    width=[w_],
                    height=[h_],
                    angle=[ang],
                    color="focus",
                    fill=True,
                    alpha=0.22,
                    linewidth=0.6,
                    z=2.6,
                    categorical=True,
                )
            )
        panel.add(LineMark(x=[0, fx], y=[0, fy], color="focus", style="dashed", width=1.0, z=2.7))
        if judged(fp) and fp.residual > POORLY_EXPLAINED:
            rr = fp.residual * halo_scale
            panel.add(
                EllipseMark(
                    cx=[fx],
                    cy=[fy],
                    width=[2 * rr],
                    height=[2 * rr],
                    color="novel",
                    linestyle="dashed",
                    linewidth=0.8,
                    z=2.6,
                    categorical=True,
                )
            )
        panel.add(PointMark(x=[fx], y=[fy], color="focus", size=34, z=4, categorical=True))
        labels.append(_Label(fx, fy, name_of(fp), "focus", 0, forced=True))
    for i, p in enumerate(sorted(novel, key=lambda q: -(q.residual * q.projection.magnitude))[:3]):
        px, py = xy(p)
        labels.append(_Label(px, py, f"{name_of(p)} · poorly explained", "novel", 10 + i))
    for i, p in enumerate(sorted(plain, key=lambda q: -q.projection.magnitude)[:label_top]):
        px, py = xy(p)
        labels.append(_Label(px, py, name_of(p), "ink_soft", 20 + i))
    for tm in place_labels(labels, R, occupied=spoke_boxes):
        panel.add(tm)

    # --- caption
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
            pct_unexplained=fp.residual if judged(fp) else None,
            n_refs=len(lib),
        )
        if fp.ellipse is not None and max(fp.ellipse[2], fp.ellipse[3]) < 0.04 * R:
            cap.add("Its position is precise: the 95% region is smaller than the marker.")
        if not judged(fp):
            cap.weak(
                f"{fp.label}'s signal is too weak (noise-corrected magnitude "
                f"{fp.projection.magnitude:.0f}) to say how much of it the references explain."
            )
        elif fp.concentration < 0.6 and fp.projection.weights.sum() > 0:
            cap.add(
                "Its weights are spread across distant references (concentration {c:.2f}), so its "
                "direction is a mixture, not a single mechanism.",
                c=fp.concentration,
            )
    cap.add(
        "{n} compounds; distance from centre = noise-corrected signature magnitude (log-fold "
        "changes vs vehicle). {k} with enough signal are poorly explained by the references "
        "(residual > {thr}, purple); {w} are too weak to judge (light grey).",
        n=len(active),
        k=len(novel),
        thr=POORLY_EXPLAINED,
        w=len(weak),
    )
    if clipped:
        cap.add(
            "{c} compound(s) beyond the outer ring ({R:g}) are drawn on it, with their magnitude "
            "in brackets.",
            c=len(clipped),
            R=R,
        )
    for ref_a, ref_b, cos_ab in correlated:
        cap.weak(
            f"the {spoke_labels.get(ref_a, ref_a)} and {spoke_labels.get(ref_b, ref_b)} references "
            f"are similar (cosine {cos_ab:.2f}): these spokes are not independent."
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
        verb = "are not shown" if gate_mode == "drop" else "are drawn as hollow grey"
        cap.weak(f"{len(gated)} compound(s) failed the retrievability gate and {verb}.")
    meth = Methods().use("compass", n_refs=len(lib), n_boot=n_boot)
    feet = [f"References: {lib.label()} ({len(lib)} mechanism classes)"]
    if left_out:
        feet.append(
            "Classes not used as spokes (members disagree, or too few): "
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
        aspect=0.82,
        meta={
            "points": points,
            "correlated": correlated,
            "rings": rings,
            "spokes": list(lib.names),
            "angles": angles,
            "clipped": sorted(clipped),
            "novel": [p.compound_id for p in novel],
            "weak": [p.compound_id for p in weak],
            "rim": R,
        },
    ).validate()


@register_classic("compass", "PCA")
def build_classic_spec(
    signatures: Mapping[str, Signature] | None = None,
    features: np.ndarray | None = None,
    focus: str | None = None,
    classes: Mapping[str, Sequence[str]] | None = None,
    class_labels: Mapping[str, str] | None = None,
    **_: object,
) -> PlotSpec:
    """The best conventional comparator: PCA with reference-class members coloured by class
    (legend), the focus compound labelled and everything else grey."""
    if not signatures or features is None:
        raise TinctureDataError("PCA fallback needs signatures.")
    ids = list(signatures)
    X = np.vstack([np.nan_to_num(signatures[c].align(features).effect) for c in ids])
    Xc = X - X.mean(axis=0)
    U, S, _vt = np.linalg.svd(Xc, full_matrices=False)
    pcs = U[:, :2] * S[:2]
    ve = S**2 / (S**2).sum()
    p = Panel(
        x=Scale(label=f"PC1 ({ve[0]:.0%})"),
        y=Scale(label=f"PC2 ({ve[1]:.0%})"),
        legend=LegendSpec(title="mechanism class", loc="best") if classes else None,
    )
    mask = np.array([c == focus for c in ids])
    member = np.zeros(len(ids), bool)
    palette = ["reference", "novel", "good", "focus_light", "reference_light", "novel_light"]
    for k, (name, members) in enumerate((classes or {}).items()):
        if k >= len(palette):
            break
        m = np.array([c in set(members) for c in ids]) & ~mask
        member |= m
        if m.any():
            lab = (class_labels or {}).get(name, name)
            p.add(PointMark(x=pcs[m, 0], y=pcs[m, 1], color=palette[k], size=12,
                            categorical=True, legend=lab))  # fmt: skip
    rest = ~mask & ~member
    p.add(PointMark(x=pcs[rest, 0], y=pcs[rest, 1], color="mid_grey", size=10,
                    legend="other" if classes else None))  # fmt: skip
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
    gate_mode: str = "drop",
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
        pts,
        lib,
        angles,
        focus,
        label_top,
        spoke_labels,
        left_out,
        title,
        n_boot,
        footer,
        gate_mode,
    )
    classic = (
        build_classic(
            "compass",
            signatures=signatures,
            features=lib.features,
            focus=focus,
            classes={n: lib.member_ids[n] for n in lib.names},
            class_labels=spoke_labels,
        )
        if len(signatures) >= 3
        else None
    )
    return PlotResult(spec=spec, classic=classic)
