"""Dose Paths: does a compound's mechanism change with dose? (HANDOVER §5.2)

Replaces PCA coloured by dose. Each compound is a path from vehicle (origin) through its
doses; dots grow with log dose, an arrowhead marks the top dose and a 95% ellipse from
posterior draws shows each dose's uncertainty. x = likeness to the compound's top reference
(from the Compass fit, leave-one-out); y = stress (stress gene-set score, or cell loss). The
reference class's own path (class-mean signature per dose, queries excluded) is drawn in
blue. A bend is annotated only when the bend test is significant
(:mod:`tincture.stats.paths`).

Guards (each tested): > 8 compounds → small multiples; < 3 doses → points only with
"too few doses for a path"; non-significant bend → no annotation.
Classic fallback: PCA of the same signatures, coloured by dose.
"""

from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from tincture.errors import TinctureDataError
from tincture.model.reference import ReferenceLibrary
from tincture.model.signature import Signature
from tincture.report.captions import TEMPLATES, Caption, fmt_conc
from tincture.report.classic import build_classic, register_classic
from tincture.report.methods import Methods, register_block
from tincture.result import PlotResult
from tincture.spec import (
    BandMark,
    ColorbarSpec,
    EllipseMark,
    LineMark,
    NoticeMark,
    Panel,
    PlotSpec,
    PointMark,
    RectMark,
    Scale,
    TextMark,
)
from tincture.stats.paths import MIN_BEND_DEG, DosePath, PathAxes, dose_path
from tincture.stats.projection import cov_ellipse, nnls_project
from tincture.stats.signatures import SignatureSet

MAX_SINGLE_PANEL = 8

register_block(
    "dose_paths",
    "Dose paths place each dose's signature at x = its projection on the unit vector of the "
    "compound's top reference (the reference with the largest NNLS weight at the top dose, "
    "fitted without the compound itself) and y = {y_desc}. The stress zone starts at the "
    "99th percentile of y over all {low} signatures; the on-target zone covers x beyond half "
    "the reference class's own top-dose x. Ellipses are 95% regions from {n_boot} posterior "
    "draws of each signature. Bends were tested on axes scaled by their spread across all "
    "signatures in the stratum: a bend is a turn of more than {min_bend:.0f}° between "
    "successive dose deltas, p = P(|turn| ≤ {min_bend:.0f}°) under the draws, "
    "Bonferroni-corrected over turning points; only bends with corrected p < 0.05 are "
    "annotated.",
)


@dataclass
class _Stratum:
    """All signatures of the stratum, aligned to the library's features."""

    E: np.ndarray  # signatures × features (NaN → 0)
    doses: np.ndarray
    lowest: float


def _aligned(sigset: SignatureSet, cid: str, feats: np.ndarray) -> list[Signature]:
    return [s.align(feats) for s in sigset.signatures if s.compound_id == cid]


def _top_reference(sig: Signature, lib: ReferenceLibrary) -> str:
    """Name of the reference with the largest NNLS weight (leave-one-out)."""
    use = lib.without(sig.compound_id) if lib.reference_of(sig.compound_id) else lib
    p = nnls_project(np.nan_to_num(sig.effect), use.matrix)
    return use.names[int(np.argmax(p.weights))]


def _ref_unit(lib: ReferenceLibrary, ref: str, cid: str) -> np.ndarray:
    use = lib.without(cid) if lib.reference_of(cid) else lib
    v = use.matrix[use.names.index(ref)]
    return v / np.linalg.norm(v)


def _reference_path(
    sigset: SignatureSet, lib: ReferenceLibrary, ref: str, exclude: set[str], axes: PathAxes
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Class-mean path of the reference's members (excluding the queries)."""
    members = set(lib.member_ids[ref]) - exclude
    by_dose: dict[float, list[np.ndarray]] = {}
    for s in sigset.signatures:
        if s.compound_id in members:
            by_dose.setdefault(s.dose, []).append(np.nan_to_num(s.align(lib.features).effect))
    doses = np.array(sorted(by_dose))
    if doses.size == 0:
        return doses, doses, doses
    E = np.vstack([np.mean(by_dose[d], axis=0) for d in doses])
    return doses, np.asarray(axes.x(E)), np.asarray(axes.y(E))


def _sizes(doses: np.ndarray, dmin: float, dmax: float) -> np.ndarray:
    span = np.log10(dmax) - np.log10(dmin) if dmax > dmin else 1.0
    return 10 + 36 * (np.log10(doses) - np.log10(dmin)) / span


def compute(
    sigset: SignatureSet,
    compounds: Sequence[str],
    library: ReferenceLibrary,
    stress_mask: np.ndarray | None,
    stratum: _Stratum,
    cell_loss: Mapping[tuple[str, float], float] | None = None,
    shared_reference: str | None = None,
    n_boot: int = 500,
) -> tuple[dict[str, DosePath], dict[str, str], dict[str, PathAxes]]:
    """Dose path per compound, its top reference and the axes used."""
    feats = library.features
    paths, refs, axes_of = {}, {}, {}
    if stress_mask is None:
        assert cell_loss is not None
        vals = np.array([v for v in cell_loss.values() if np.isfinite(v)])
        sy = float(vals.std()) if vals.size > 1 else 1.0
    else:
        y_all = stratum.E[:, stress_mask].mean(axis=1) - stratum.E.mean(axis=1)
        sy = float(y_all.std())
    for cid in compounds:
        sigs = _aligned(sigset, cid, feats)
        if not sigs:
            raise TinctureDataError(f"No signatures for {cid!r}.", fix="Check the compound id.")
        top = max(sigs, key=lambda s: s.dose)
        ref = shared_reference or _top_reference(top, library)
        unit = _ref_unit(library, ref, cid)
        sx = float((stratum.E @ unit).std())
        axes = PathAxes(ref_unit=unit, stress_mask=stress_mask, scale=(sx or 1.0, sy or 1.0))
        yv = None
        if stress_mask is None:
            assert cell_loss is not None
            yv = np.array([cell_loss.get((cid, s.dose), np.nan) for s in sigs])
            if not np.all(np.isfinite(yv)):
                raise TinctureDataError(f"cell_loss is missing doses of {cid!r}.")
        paths[cid] = dose_path(sigs, axes, yv, n_boot=n_boot)
        refs[cid] = ref
        axes_of[cid] = axes
    return paths, refs, axes_of


def _bend_words(dp: DosePath, b: int, axes: PathAxes, ref_label: str, stress_word: str) -> str:
    """Where the path heads after the bend, judged on scaled axes."""
    dx = (dp.x[b] - dp.x[b - 1]) / axes.scale[0]
    dy = (dp.y[b] - dp.y[b - 1]) / axes.scale[1]
    if abs(dy) >= abs(dx):
        return stress_word if dy > 0 else f"less {stress_word}"
    return ref_label if dx > 0 else f"away from {ref_label}"


def _panel(
    paths: list[DosePath],
    ref: str,
    axes: PathAxes,
    sigset: SignatureSet,
    library: ReferenceLibrary,
    names: Mapping[str, str],
    focus: str | None,
    y_label: str,
    stress_thr: float,
    ref_label: str,
    stress_word: str,
    key: str,
    title: str,
    dmin: float,
    dmax: float,
    show_ref_path: bool,
    label_x: bool = True,
    label_y: bool = True,
    compact: bool = False,
    y_domain: tuple[float, float] | None = None,
) -> tuple[Panel, list[str]]:
    notes: list[str] = []
    exclude = {dp.compound_id for dp in paths}
    if show_ref_path:
        rd, rx, ry = _reference_path(sigset, library, ref, exclude, axes)
    else:
        rd = rx = ry = np.array([])
    xs, ys = [0.0, *rx], [0.0, *ry, stress_thr]
    ells: dict[str, list[tuple[float, float, float, float, float]]] = {}
    for dp in paths:
        xs += list(dp.x)
        ys += list(dp.y)
        if dp.draws is not None and (dp.compound_id == focus or len(paths) == 1):
            ells[dp.compound_id] = [cov_ellipse(dp.draws[:, j, :]) for j in range(dp.n_doses)]
    x_lo, x_hi = min(xs), max(xs)
    y_lo, y_hi = min(ys), max(ys)
    xr, yr = (x_hi - x_lo) or 1.0, (y_hi - y_lo) or 1.0
    dom_x = (x_lo - 0.06 * xr, x_hi + 0.24 * xr)  # room for end labels
    dom_y = y_domain or (y_lo - 0.1 * yr, y_hi + 0.12 * yr)
    p = Panel(
        key=key,
        title=title,
        x=Scale("linear", label=f"{ref_label}-likeness →" if label_x else "", domain=dom_x),
        y=Scale("linear", label=y_label if label_y else "", domain=dom_y),
    )
    p.add(
        BandMark(
            lo=stress_thr,
            hi=dom_y[1],
            span=True,
            orientation="h",
            color="focus_light",
            alpha=0.28,
            inner_label=f"↑ {stress_word}",
            z=0.3,
        )
    )
    if rx.size and rx.max() > 0:
        x_on = 0.5 * float(rx.max())
        p.add(
            RectMark(
                x0=[x_on],
                y0=[dom_y[0]],
                width=[dom_x[1] - x_on],
                height=[stress_thr - dom_y[0]],
                color="window",
                edgecolor="good",
                linestyle="dashed",
                alpha=0.3,
                z=0.3,
            )
        )
        p.add(
            TextMark(
                x=dom_x[1],
                y=stress_thr,
                text="on-target " if not compact else "on-target",
                color="good",
                ha="right",
                va="top",
                categorical=True,
            )
        )
        p.add(
            LineMark(
                x=[0, *rx],
                y=[0, *ry],
                color="reference",
                width=1.6,
                arrow_end=True,
                categorical=True,
                z=2,
                end_label=ref_label if compact else f"{ref_label} (class mean)",
            )
        )
        p.add(
            PointMark(
                x=rx,
                y=ry,
                color="reference_light",
                size=_sizes(rd, dmin, dmax),
                edgecolor="reference",
                z=2.1,
            )
        )
    p.add(PointMark(x=[0], y=[0], color="vehicle", size=24, edgecolor="ink_soft", z=3))
    if not compact:
        p.add(TextMark(x=0, y=0, text="vehicle  ", color="ink_soft", ha="right", va="top"))
    mid = 0.5 * (dom_x[0] + dom_x[1])
    for dp in paths:
        is_focus = dp.compound_id == focus or len(paths) == 1
        col = "focus" if is_focus else "mid_grey"
        name = names.get(dp.compound_id, dp.compound_id)
        for cx, cy, w, h, a in ells.get(dp.compound_id, []):
            p.add(
                EllipseMark(
                    cx=[cx],
                    cy=[cy],
                    width=[w],
                    height=[h],
                    angle=a,
                    color=col,
                    fill=True,
                    alpha=0.16,
                    linewidth=0,
                    z=1.5,
                )
            )
        if dp.n_doses >= 3:
            p.add(
                LineMark(
                    x=[0, *dp.x],
                    y=[0, *dp.y],
                    color=col,
                    width=1.8 if is_focus else 1.0,
                    arrow_end=True,
                    categorical=is_focus,
                    z=2.6 if is_focus else 2.4,
                    end_label=name,
                )
            )
        else:
            notes.append(f"{name}: too few doses for a path")
            p.add(
                TextMark(
                    x=float(dp.x[-1]),
                    y=float(dp.y[-1]),
                    text=f"  {name}",
                    color=col,
                    ha="left",
                    va="center",
                    categorical=is_focus,
                )
            )
        p.add(
            PointMark(
                x=dp.x,
                y=dp.y,
                color=col,
                size=_sizes(dp.doses, dmin, dmax),
                categorical=is_focus,
                z=2.7 if is_focus else 2.5,
            )
        )
        b = dp.bend_index
        if b is not None:
            where = _bend_words(dp, b, axes, ref_label, stress_word)
            d_txt = fmt_conc(float(dp.doses[b - 1]))
            words = f"bend > {d_txt}" if compact else f"bends above {d_txt} → {where}"
            bx, by = float(dp.x[b - 1]), float(dp.y[b - 1])
            right = bx > mid
            p.add(
                TextMark(
                    x=bx,
                    y=by,
                    text=f"{words}  " if right else f"  {words}",
                    color="ink" if is_focus else "ink_soft",
                    ha="right" if right else "left",
                    va="bottom",
                    size=None if is_focus else 6.5,
                )
            )
    if notes:
        p.add(NoticeMark(text="too few doses for a path", loc="upper left", level="weak"))
    return p, notes


def build_spec(
    paths: dict[str, DosePath],
    refs: dict[str, str],
    axes_of: dict[str, PathAxes],
    sigset: SignatureSet,
    library: ReferenceLibrary,
    names: Mapping[str, str],
    focus: str | None,
    y_axis: str,
    stress_thr: float,
    spoke_labels: Mapping[str, str],
    n_boot: int,
    title: str,
    footer: Sequence[str],
    dropped: Sequence[str] = (),
) -> PlotSpec:
    stress_word = "stress" if y_axis == "stress" else "cell loss"
    y_label = (
        "↑ stress (stress genes vs all genes, log2FC)"
        if y_axis == "stress"
        else "↑ cell loss (−log2 cells vs vehicle)"
    )
    all_d = np.concatenate([dp.doses for dp in paths.values()])
    dmin, dmax = float(all_d.min()), float(all_d.max())
    lab = {c: spoke_labels.get(refs[c], refs[c]) for c in paths}
    cap = Caption()
    ids = list(paths)
    panels: list[Panel] = []
    notes: list[str] = []
    small = len(ids) > MAX_SINGLE_PANEL
    show_ref = y_axis == "stress"
    if not small:
        anchor = focus if focus in paths else ids[0]
        pnl, notes = _panel(
            [paths[c] for c in ids], refs[anchor], axes_of[anchor], sigset, library, names,
            focus, y_label, stress_thr, lab[anchor], stress_word, "main", "", dmin, dmax,
            show_ref,
        )  # fmt: skip
        panels, nrows, ncols = [pnl], 1, 1
    else:
        ncols = 4
        nrows = int(np.ceil(len(ids) / ncols))
        short_y = f"↑ {stress_word} (log2FC)" if y_axis == "stress" else f"↑ {stress_word}"
        # shared y so stress is comparable across panels (x differs: one reference each)
        ys = [0.0, stress_thr, *(v for dp in paths.values() for v in dp.y)]
        if show_ref:
            for c in ids:
                ys += list(_reference_path(sigset, library, refs[c], {c}, axes_of[c])[2])
        yr = (max(ys) - min(ys)) or 1.0
        y_dom = (min(ys) - 0.08 * yr, max(ys) + 0.1 * yr)
        for i, c in enumerate(ids):
            r, col = divmod(i, ncols)
            bottom = i + ncols >= len(ids)
            pnl, nt = _panel(
                [paths[c]], refs[c], axes_of[c], sigset, library, names, c, short_y,
                stress_thr, lab[c], stress_word, f"p{i}", names.get(c, c), dmin, dmax,
                show_ref, label_x=bottom, label_y=col == 0, compact=True, y_domain=y_dom,
            )  # fmt: skip
            pnl.row, pnl.col = r, col
            panels.append(pnl)
            notes += nt
        cap.add(
            "{n} compounds, one panel each (small multiples above {k} compounds).",
            n=len(ids),
            k=MAX_SINGLE_PANEL,
        )
    lead = focus if focus in paths else ids[0]
    dp = paths[lead]
    nm = names.get(lead, lead)
    b = dp.bend_index
    if dp.n_doses < 3:
        cap.weak(f"{nm} has {dp.n_doses} dose(s): too few doses for a path.")
    elif b is not None:
        d = fmt_conc(float(dp.doses[b - 1]))
        cap.add(
            TEMPLATES["dose_paths"],
            compound=nm,
            ref=lab[lead],
            dose_ok=d,
            dose_bend=d,
            stress_axis=_bend_words(dp, b, axes_of[lead], lab[lead], stress_word),
            p=dp.bend_p,
        )
    else:
        cap.add(
            "{c} moves from vehicle towards {r} across {lo}–{hi} without a significant bend "
            "(smallest corrected p = {p:.2g}).",
            c=nm,
            r=lab[lead],
            lo=fmt_conc(float(dp.doses[0])),
            hi=fmt_conc(float(dp.doses[-1])),
            p=dp.bend_p,
        )
    others = [c for c in ids if c != lead and paths[c].bend_index is not None]
    if len(ids) > 1:
        cap.add(
            "{k} of {n} other compounds bend significantly{w}.",
            k=len(others),
            n=len(ids) - 1,
            w=(": " + ", ".join(names.get(c, c) for c in others)) if others else "",
        )
    if dropped:
        cap.add("{n} compound(s) failed the retrievability gate and are not shown.", n=len(dropped))
    for note in notes:
        cap.weak(note + ".")
    low = float(min(s.dose for s in sigset.signatures))
    meth = Methods().use(
        "dose_paths",
        y_desc=(
            "the mean log2 fold change of the stress genes (union of MSigDB Hallmark apoptosis, "
            "p53 pathway and unfolded protein response) minus the mean over all genes"
            if y_axis == "stress"
            else "cell loss, −log2(cells per well ÷ vehicle cells per well)"
        ),
        low=fmt_conc(low),
        n_boot=n_boot,
        min_bend=MIN_BEND_DEG,
    )
    return PlotSpec(
        plot="dose_paths",
        panels=panels,
        nrows=nrows,
        ncols=ncols,
        title=title,
        caption=cap.render(),
        methods=meth.render(),
        footer=[
            f"References: {library.label()}. Grey dot = vehicle; dot size grows with log dose; "
            "shaded ellipses = 95% regions (focus compound).",
            *footer,
        ],
        warnings=list(cap.notes),
        aspect=0.62 if not small else 0.27 * nrows,
        meta={
            "paths": paths,
            "refs": refs,
            "stress_threshold": stress_thr,
            "small_multiples": small,
            "dropped": list(dropped),
        },
    ).validate()


@register_classic("dose_paths", "PCA coloured by dose")
def build_classic_spec(
    sigset: SignatureSet | None = None,
    compounds: Sequence[str] = (),
    features: np.ndarray | None = None,
    **_: object,
) -> PlotSpec:
    if sigset is None or features is None:
        raise TinctureDataError("PCA fallback needs signatures.")
    keep = set(compounds)
    sigs = [s for s in sigset.signatures if s.compound_id in keep]
    if len(sigs) < 3:
        raise TinctureDataError("PCA fallback needs at least 3 signatures.")
    X = np.vstack([np.nan_to_num(s.align(features).effect) for s in sigs])
    U, S, _vt = np.linalg.svd(X - X.mean(axis=0), full_matrices=False)
    pcs = U[:, :2] * S[:2]
    ve = S**2 / (S**2).sum()
    doses = np.array([s.dose for s in sigs])
    ud = np.unique(doses)
    p = Panel(
        x=Scale(label=f"PC1 ({ve[0]:.0%})"),
        y=Scale(label=f"PC2 ({ve[1]:.0%})"),
        colorbar=ColorbarSpec(
            label="dose",
            ticks=[float(v) for v in np.log10(ud)],
            tick_labels=[fmt_conc(float(d)) for d in ud],
        ),
    )
    p.add(PointMark(x=pcs[:, 0], y=pcs[:, 1], values=np.log10(doses), cmap="viridis", size=16))
    return PlotSpec(
        plot="dose_paths_classic",
        panels=[p],
        title="PCA coloured by dose",
        caption="PCA of the same signatures, coloured by dose.",
        aspect=0.75,
    ).validate()


def plot(
    sigset: SignatureSet,
    compounds: Sequence[str],
    library: ReferenceLibrary,
    stress_genes: Collection[str] | None = None,
    focus: str | None = None,
    y_axis: str = "stress",
    cell_loss: Mapping[tuple[str, float], float] | None = None,
    gate: Collection[str] | None = None,
    spoke_labels: Mapping[str, str] | None = None,
    reference: str | None = None,
    n_boot: int = 500,
    title: str = "Dose Paths",
    footer: Sequence[str] = (),
) -> PlotResult:
    """Dose Paths for ``compounds`` from one stratum's :class:`SignatureSet`.

    ``stress_genes``: gene symbols for the stress axis (see
    :func:`tincture.stats.paths.stress_gene_set`). ``y_axis="cell_loss"`` uses
    ``cell_loss[(compound, dose)]`` instead. ``gate``: compound ids that passed the
    retrievability gate; others are dropped and the caption says so. ``reference``: use
    this reference for x in every panel (default: the focus compound's top reference for one
    panel, each compound's own top reference in small multiples).
    """
    from tincture.plots.compass import select_spokes, short_name

    if y_axis not in ("stress", "cell_loss"):
        raise ValueError("y_axis must be 'stress' or 'cell_loss'.")
    compounds = list(dict.fromkeys(compounds))
    dropped = [c for c in compounds if gate is not None and c not in gate]
    compounds = [c for c in compounds if c not in dropped]
    if not compounds:
        raise TinctureDataError("No compounds to plot (none given, or all failed the gate).")
    if focus is not None and focus not in compounds:
        raise TinctureDataError(f"Focus {focus!r} is not among the plotted compounds.")
    lib, _ = select_spokes(library)
    feats = lib.features
    E = np.vstack([np.nan_to_num(s.align(feats).effect) for s in sigset.signatures])
    doses = np.array([s.dose for s in sigset.signatures])
    stratum = _Stratum(E=E, doses=doses, lowest=float(doses.min()))
    mask = None
    if y_axis == "stress":
        if not stress_genes:
            raise TinctureDataError(
                "No stress genes given.", fix="Pass tincture.stats.paths.stress_gene_set(...)."
            )
        sg = set(stress_genes)
        sym = sigset.feature_symbols if sigset.feature_symbols is not None else sigset.features
        sym_of = dict(zip(np.asarray(sigset.features).astype(str), sym, strict=True))
        mask = np.array([str(sym_of.get(f, f)) in sg for f in feats])
        if mask.sum() < 10:
            raise TinctureDataError(
                f"Only {int(mask.sum())} stress genes among the features.",
                fix="Check gene identifiers (symbols vs Ensembl).",
            )
        lo = E[doses == stratum.lowest]
        thr = float(np.quantile(lo[:, mask].mean(axis=1) - lo.mean(axis=1), 0.99))
    else:
        if cell_loss is None:
            raise TinctureDataError("y_axis='cell_loss' needs cell_loss[(compound, dose)].")
        vals = np.array([v for (_c, d), v in cell_loss.items() if d == stratum.lowest])
        vals = vals[np.isfinite(vals)]
        thr = float(np.quantile(vals, 0.99)) if vals.size else 0.0
    shared = reference
    if reference is not None and reference not in lib.names:
        raise TinctureDataError(
            f"Reference {reference!r} is not a spoke.", fix=f"Use one of {lib.names}."
        )
    if shared is None and len(compounds) <= MAX_SINGLE_PANEL:  # x must mean the same
        anchor = focus or compounds[0]
        sigs = _aligned(sigset, anchor, feats)
        if not sigs:
            raise TinctureDataError(f"No signatures for {anchor!r}.")
        shared = _top_reference(max(sigs, key=lambda s: s.dose), lib)
    paths, refs, axes_of = compute(sigset, compounds, lib, mask, stratum, cell_loss, shared, n_boot)
    names: dict[str, str] = {}
    for s in sigset.signatures:
        names.setdefault(
            s.compound_id, short_name(str(s.meta.get("compound_name") or s.compound_id))
        )
    spec = build_spec(
        paths, refs, axes_of, sigset, lib, names, focus, y_axis, thr, spoke_labels or {},
        n_boot, title, footer, dropped,
    )  # fmt: skip
    classic = None
    if sum(s.compound_id in set(compounds) for s in sigset.signatures) >= 3:
        classic = build_classic("dose_paths", sigset=sigset, compounds=compounds, features=feats)
    return PlotResult(spec=spec, classic=classic)
