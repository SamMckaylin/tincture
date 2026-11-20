"""Plate Atlas: can we trust these plates? (HANDOVER §5.6; QC / supplementary figure)

Left: one plate as a map, values as % of that plate's negative-control median. 96/384-well
plates draw controls as circles and other wells as squares; 1536-well plates switch to a
density (heatmap) view automatically. Wells absent from the source are drawn grey.
Right, against plate order: Z′ per plate (0.5 line; failing plates terracotta), edge ratio
(1.0 = no edge effect), and negative-control level (drift, Theil–Sen trend).
Classic fallback: a plain heatmap of the same plate.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from tincture.errors import TinctureDataError
from tincture.model.plate import PlateTable, row_label
from tincture.report.captions import Caption, fmt_p
from tincture.report.classic import build_classic, register_classic
from tincture.report.methods import Methods, register_block
from tincture.result import PlotResult
from tincture.spec import (
    BandMark,
    ColorbarSpec,
    HeatmapMark,
    LineMark,
    NoticeMark,
    Panel,
    PlotSpec,
    PointMark,
    RectMark,
    RefLineMark,
    Scale,
)
from tincture.stats.zprime import PlateQC, Selector, plate_qc

register_block(
    "zprime",
    "Plate quality was summarised by Z′ = 1 − 3(σ₊ + σ₋)/|μ₊ − μ₋| with {pos} as the "
    "positive and {neg} as the negative reference; plates with Z′ < {thr} were flagged. "
    "Edge effect is the median signal of the outermost seeded ring over that of wells at "
    "least two rings further in ({edge} wells with signal above 10% of the plate's "
    "negative-control median); drift is the Spearman correlation and Theil–Sen trend of "
    "each plate's negative-control median, relative to comparable plates, against plate "
    "order.",
    (
        "Zhang JH, Chung TDY, Oldenburg KR. A simple statistical parameter for use in "
        "evaluation and validation of high throughput screening assays. J Biomol Screen. "
        "1999;4:67–73.",
    ),
)

DENSITY_THRESHOLD = 1536
#: Above this many plates, per-plate panels show binned median + 5–95% bands.
BIN_THRESHOLD = 1000
N_BINS = 120


def _binned(x: np.ndarray, y: np.ndarray, n_bins: int = N_BINS) -> tuple[np.ndarray, ...]:
    """Median and 5–95% quantiles of ``y`` in equal-count bins of ``x``."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    order = np.argsort(x)
    chunks = [c for c in np.array_split(order, min(n_bins, max(1, len(order)))) if len(c)]
    cx = np.array([x[c].mean() for c in chunks])
    med = np.array([np.median(y[c]) for c in chunks])
    lo = np.array([np.quantile(y[c], 0.05) for c in chunks])
    hi = np.array([np.quantile(y[c], 0.95) for c in chunks])
    return cx, med, lo, hi


def _summary_marks(panel: Panel, x: np.ndarray, y: np.ndarray, colour: str = "reference") -> None:
    cx, med, lo, hi = _binned(x, y)
    panel.add(
        BandMark(
            x=cx, lo=lo, hi=hi, color=f"{colour}_light", alpha=0.75, inner_label="5–95% per bin"
        ),
        LineMark(x=cx, y=med, color=colour, width=1.4, categorical=True, end_label="median"),
    )


def _choose_plate(qc: PlateQC, plate: str | None) -> tuple[str, str]:
    t = qc.per_plate
    if plate is not None:
        if plate not in set(t["plate"]):
            raise TinctureDataError(f"Plate {plate!r} not in data.", fix="See PlateTable.plates.")
        return plate, "selected"
    valid = t.dropna(subset=["z_prime"])
    if valid.empty:
        return str(t["plate"].iloc[0]), "first"
    if valid["fails"].any():
        return str(valid.sort_values("z_prime").iloc[0]["plate"]), "worst Z′"
    mid = valid.iloc[(valid["z_prime"] - valid["z_prime"].median()).abs().argsort().iloc[0]]
    return str(mid["plate"]), "median Z′"


#: Plate-map colour: log2(signal / plate's negative-control median), clipped at ±1
#: (half / double); the neutral grey bin is ±0.2 log2 (≈ ±15%).
PLATE_CLIP = 1.0


def _plate_panel(pt: PlateTable, qc: PlateQC, plate: str, why: str) -> Panel:
    """Plate map in the whiteboard style: rounded tiles with gutters, controls as circles,
    diverging colour vs the plate's negative controls, unreported wells as empty outlines."""
    nr, nc = pt.shape
    d = pt.data[pt.data["plate"].astype(str) == plate]
    neg_med = float(qc.per_plate.set_index("plate").loc[plate, "neg_median"])
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = (
            np.log2(np.maximum(d["value"].to_numpy(float), 1e-9) / neg_med)
            if neg_med > 0
            else np.full(len(d), np.nan)
        )
    m = np.full((nr, nc), np.nan)
    r, c = d["row"].to_numpy(), d["col"].to_numpy()
    m[r, c] = rel
    circles = np.zeros((nr, nc), bool)
    ctrl = d["well_type"].astype(str).isin(["neg_ctrl", "pos_ctrl", "blank"]).to_numpy()
    circles[r[ctrl], c[ctrl]] = True
    step_c = 1 if nc <= 12 else (4 if nc <= 24 else 12)
    step_r = 1 if nr <= 8 else (4 if nr <= 16 else 8)
    cols_t = [0, *range(step_c - 1, nc, step_c)] if step_c > 1 else list(range(nc))
    rows_t = [0, *range(step_r - 1, nr, step_r)] if step_r > 1 else list(range(nr))
    cols_t = sorted(set(cols_t) | {nc - 1})
    rows_t = sorted(set(rows_t) | {nr - 1})
    p = Panel(
        key="plate",
        x=Scale(
            "linear",
            label="Column",
            ticks=[float(i) for i in cols_t],
            tick_labels=[str(i + 1) for i in cols_t],
        ),
        y=Scale(
            "linear",
            label="Row",
            ticks=[float(i) for i in rows_t],
            tick_labels=[row_label(i) for i in rows_t],
        ),
        colorbar=ColorbarSpec(
            label="signal vs negative controls",
            clip_note=f"log2, clipped at ±{PLATE_CLIP:g}",
            orientation="horizontal",
            ticks=[-PLATE_CLIP, 0.0, PLATE_CLIP],
            tick_labels=["÷2", "same", "×2"],
            compact=True,
        ),
        title=f"Plate {plate} ({why})",
        aspect="equal",
        frame=False,
    )
    p.add(
        HeatmapMark(
            values=m,
            cmap="diverging",
            clip=PLATE_CLIP,
            circles=circles,
            style="rounded",
            cell_edge=None,
        )
    )
    return p


def build_spec(
    pt: PlateTable,
    qc: PlateQC,
    plate: str | None = None,
    title: str = "Plate Atlas",
    footer: tuple[str, ...] = (),
) -> PlotSpec:
    t = qc.per_plate
    shown, why = _choose_plate(qc, plate)
    left = _plate_panel(pt, qc, shown, why)
    left.rowspan = 3
    left.row, left.col = 0, 0
    n = len(t)
    x = t["order"].to_numpy(float)
    z = t["z_prime"].to_numpy(float)
    zmin = -1.0
    clipped = int(np.sum(z < zmin))
    zc = np.where(z < zmin, zmin, z)
    zlo = float(min(np.nanmin(zc) - 0.1, qc.threshold - 0.2)) if np.isfinite(zc).any() else zmin
    zlo = max(zlo, zmin - 0.05)
    xs = Scale(
        "linear", label="Plate (scan order)", domain=(x.min() - 0.5, x.max() + 0.5) if n else None
    )
    zp = Panel(
        key="zprime",
        x=xs,
        y=Scale("linear", label="Z′", domain=(zlo, 1.0)),
        row=0,
        col=1,
        title="Z′ per plate",
    )
    fail = t["fails"].to_numpy()
    passed = ~fail & np.isfinite(z)
    binned = n > BIN_THRESHOLD
    if binned:
        fin = np.isfinite(zc)
        _summary_marks(zp, x[fin], zc[fin], colour="good")
        if fail.any():
            zp.add(
                PointMark(
                    x=x[fail],
                    y=zc[fail],
                    color="focus",
                    size=2,
                    categorical=True,
                    edgecolor="focus",
                )
            )
    elif n <= 40:
        ok_w = 0.8
        if passed.any():
            zp.add(
                RectMark(
                    x0=x[passed] - ok_w / 2,
                    y0=np.full(passed.sum(), zlo),
                    width=ok_w,
                    height=zc[passed] - zlo,
                    color="window",
                    edgecolor="good",
                    categorical=True,
                )
            )
        if fail.any():
            zp.add(
                RectMark(
                    x0=x[fail] - ok_w / 2,
                    y0=np.full(fail.sum(), zlo),
                    width=ok_w,
                    height=np.maximum(zc[fail] - zlo, 0),
                    color="focus",
                    edgecolor="focus",
                    categorical=True,
                )
            )
    else:
        if passed.any():
            zp.add(
                PointMark(
                    x=x[passed],
                    y=zc[passed],
                    color="good",
                    size=6,
                    alpha=0.8,
                    categorical=True,
                )
            )
        if fail.any():
            zp.add(
                PointMark(
                    x=x[fail],
                    y=zc[fail],
                    color="focus",
                    size=8,
                    categorical=True,
                )
            )
    zp.add(
        RefLineMark(axis="y", value=qc.threshold, label=f"Z′ = {qc.threshold:g}", style="dashed")
    )
    if clipped:
        zp.add(NoticeMark(text=f"{clipped} plates with Z′ < −1 drawn at −1", loc="lower right"))

    ep = Panel(
        key="edge",
        x=xs,
        y=Scale("linear", label="Edge ÷ interior"),
        row=1,
        col=1,
        title="Edge effect",
    )
    er = t["edge_ratio"].to_numpy(float)
    ok = np.isfinite(er)
    if ok.any() and binned:
        _summary_marks(ep, x[ok], er[ok])
    elif ok.any():
        ep.add(
            PointMark(
                x=x[ok],
                y=er[ok],
                color="reference",
                size=3 if n > 40 else 12,
                alpha=0.6 if n > 40 else 1.0,
                categorical=True,
            )
        )
    else:
        ep.add(NoticeMark(text="edge ratio not estimable", loc="center"))
    ep.add(RefLineMark(axis="y", value=1.0, label="no edge effect", style="dotted"))

    cp = Panel(
        key="drift",
        x=xs,
        y=Scale(
            "log",
            label="Control ÷ median",
            ticks=[0.25, 0.5, 1.0, 2.0, 4.0],
            tick_labels=["0.25", "0.5", "1", "2", "4"],
        ),
        row=2,
        col=1,
        title=f"Control drift (vs {qc.drift_basis})",
    )
    lv = t["neg_level"].to_numpy(float)
    okl = np.isfinite(lv)
    if binned:
        _summary_marks(cp, x[okl], lv[okl])
    else:
        cp.add(
            PointMark(
                x=x[okl],
                y=lv[okl],
                color="reference",
                size=3 if n > 40 else 12,
                alpha=0.6 if n > 40 else 1.0,
                categorical=True,
            )
        )
    cp.add(RefLineMark(axis="y", value=1.0, style="dotted"))
    if np.isfinite(qc.drift_change) and okl.sum() >= 3:
        xx = np.array([np.nanmin(x[okl]), np.nanmax(x[okl])])
        mid = float(np.exp(np.nanmedian(np.log(lv[okl]))))
        yy = mid * np.array([1.0, 1.0 + qc.drift_change]) / np.sqrt(1.0 + qc.drift_change)
        cp.add(LineMark(x=xx, y=yy, color="ink", width=1.0))

    cap = Caption()
    zf = t["z_prime"].dropna()
    n_fail = int(t["fails"].sum())
    cap.add(
        "{n:,} plates ({w}-well); median Z′ {med:.2f} (IQR {q1:.2f}–{q3:.2f}); "
        "{nf:,} ({pct:.0%}) fall below {thr:g}.",
        n=n,
        w=pt.n_wells,
        med=float(zf.median()) if len(zf) else None,
        q1=float(zf.quantile(0.25)) if len(zf) else None,
        q3=float(zf.quantile(0.75)) if len(zf) else None,
        nf=n_fail,
        pct=n_fail / max(len(zf), 1),
        thr=qc.threshold,
    )
    if binned:
        cap.add(
            "Right panels summarise plates in {b} equal-count bins (blue line = median, band = "
            "5–95%); every failing plate is shown as a terracotta point.",
            b=min(N_BINS, n),
        )
    if 0 < n_fail <= 5:
        cap.add("Failing plates: {lst}.", lst=", ".join(t.loc[t["fails"], "plate"].astype(str)))
    med_edge = float(np.nanmedian(er)) if ok.any() else None
    cap.add(
        "The outermost seeded ring reads {e:.0%} of interior wells (median over plates; "
        "{basis} wells with cells).",
        e=med_edge,
        basis=qc.edge.describe(),
    )
    if qc.n_unseeded_outer:
        cap.add(
            "{k:,} of {n:,} plates leave the outer ring unseeded; their edge is the next seeded "
            "ring.",
            k=qc.n_unseeded_outer,
            n=n,
        )
    cap.add(
        "Negative-control level (relative to {basis}) changes {s:+.0%} across the {n:,} plates in "
        "scan order (Theil–Sen; Spearman ρ = {rho:.2f}, p {p}).",
        basis=qc.drift_basis,
        s=qc.drift_change if np.isfinite(qc.drift_change) else None,
        n=n,
        rho=qc.drift_rho if np.isfinite(qc.drift_rho) else None,
        p=fmt_p(qc.drift_p),
    )
    n_nan = int(t["z_prime"].isna().sum())
    if n_nan:
        cap.weak(
            f"Z′ could not be computed for {n_nan} plate(s) with fewer than 3 positive or "
            "negative controls."
        )
    if qc.n_missing_wells:
        cap.weak(
            f"{qc.n_missing_wells:,} well positions are not reported in the source; they are "
            "drawn grey, not as zero."
        )
    meth = Methods().use(
        "zprime",
        pos=qc.pos.describe(),
        neg=qc.neg.describe(),
        thr=qc.threshold,
        edge=qc.edge.describe(),
    )
    return PlotSpec(
        plot="plate_atlas",
        panels=[left, zp, ep, cp],
        nrows=3,
        ncols=2,
        width_ratios=[1.35, 1.0],
        title=title,
        caption=cap.render(),
        methods=meth.render(),
        footer=[
            f"Z′: positive = {qc.pos.describe()}, negative = {qc.neg.describe()}; "
            f"edge ratio from {qc.edge.describe()} wells",
            *footer,
        ],
        warnings=list(cap.notes),
        meta={
            "shown_plate": shown,
            "why": why,
            "n_fail": n_fail,
            "density_view": pt.n_wells >= DENSITY_THRESHOLD,
            "qc": t,
        },
        aspect=0.62,
    ).validate()


@register_classic("plate_atlas", "plain plate heatmap")
def build_classic_spec(
    pt: PlateTable | None = None, plate: str | None = None, **_: object
) -> PlotSpec:
    if pt is None or plate is None:
        raise TinctureDataError("Classic plate heatmap needs a PlateTable and a plate id.")
    m = pt.plate_matrix(plate)
    p = Panel(
        x=Scale("linear", label="Column"),
        y=Scale("linear", label="Row"),
        colorbar=ColorbarSpec(label=pt.value_label),
        aspect="equal",
    )
    finite = m[np.isfinite(m)]
    p.add(
        HeatmapMark(
            values=m,
            cmap="viridis",
            vmin=float(finite.min()) if finite.size else 0,
            vmax=float(finite.max()) if finite.size else 1,
            cell_edge=None,
        )
    )
    return PlotSpec(
        plot="plate_atlas_classic",
        panels=[p],
        title=f"Plate {plate}",
        caption=f"Raw {pt.value_label} for plate {plate}.",
        aspect=0.7,
    ).validate()


def plot(
    pt: PlateTable,
    neg: Selector | None = None,
    pos: Selector | None = None,
    edge: Selector | str = "auto",
    plate: str | None = None,
    title: str = "Plate Atlas",
    footer: tuple[str, ...] = (),
) -> PlotResult:
    """Plate QC figure. See :func:`tincture.stats.zprime.plate_qc` for the selectors."""
    qc = plate_qc(pt, neg=neg, pos=pos, edge=edge)
    spec = build_spec(pt, qc, plate, title, footer)
    classic = build_classic("plate_atlas", pt=pt, plate=spec.meta["shown_plate"])
    spec.meta["qc_object"] = qc
    return PlotResult(spec=spec, classic=classic)


def qc_table(result: PlotResult) -> pd.DataFrame:
    return result.meta["qc"]
