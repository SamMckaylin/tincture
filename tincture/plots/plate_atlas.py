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
    ColorbarSpec,
    HeatmapMark,
    LegendSpec,
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
from tincture.theme.palettes import sample_cmap

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


def _plate_panel(pt: PlateTable, qc: PlateQC, plate: str, why: str) -> Panel:
    nr, nc = pt.shape
    d = pt.data[pt.data["plate"].astype(str) == plate]
    neg_med = float(qc.per_plate.set_index("plate").loc[plate, "neg_median"])
    rel = d["value"].to_numpy(float) / neg_med * 100 if neg_med > 0 else np.full(len(d), np.nan)
    vmax = float(max(120.0, np.nanpercentile(rel, 99))) if np.isfinite(rel).any() else 120.0
    title = f"Plate {plate} ({why})"
    if pt.n_wells >= DENSITY_THRESHOLD:
        m = np.full((nr, nc), np.nan)
        m[d["row"].to_numpy(), d["col"].to_numpy()] = rel
        p = Panel(
            key="plate",
            x=Scale(
                "linear",
                label="Column",
                ticks=[0, 11, 23, 35, nc - 1],
                tick_labels=["1", "12", "24", "36", str(nc)],
            ),
            y=Scale(
                "linear",
                label="Row",
                ticks=[0, 7, 15, 23, nr - 1],
                tick_labels=[row_label(i) for i in (0, 7, 15, 23, nr - 1)],
            ),
            colorbar=ColorbarSpec(
                label="% of negative-control median", clip_note=f"scale 0–{vmax:.0f}%"
            ),
            title=title,
            aspect="equal",
        )
        p.add(HeatmapMark(values=m, cmap="viridis", vmin=0.0, vmax=vmax, cell_edge=None))
    else:
        p = Panel(
            key="plate",
            x=Scale(
                "linear",
                label="Column",
                domain=(-0.7, nc - 0.3),
                ticks=list(range(0, nc, max(1, nc // 6))),
                tick_labels=[str(i + 1) for i in range(0, nc, max(1, nc // 6))],
            ),
            y=Scale(
                "linear",
                label="Row",
                domain=(-0.7, nr - 0.3),
                invert=True,
                ticks=list(range(nr)),
                tick_labels=[row_label(i) for i in range(nr)],
            ),
            legend=LegendSpec(loc="upper left", outside=True),
            title=title,
            aspect="equal",
        )
        ctrl = d["well_type"].isin(["neg_ctrl", "pos_ctrl", "blank"]).to_numpy()
        cols = sample_cmap("viridis", rel, 0.0, vmax)
        for mask, marker, lab in ((ctrl, "o", "controls"), (~ctrl, "s", "other wells")):
            if mask.any():
                p.add(
                    PointMark(
                        x=d["col"].to_numpy()[mask],
                        y=d["row"].to_numpy()[mask],
                        color=[c for c, k in zip(cols, mask, strict=True) if k],
                        marker=marker,
                        size=60 if nr <= 8 else 14,
                        legend=lab,
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
    colour_note = (
        []
        if pt.n_wells >= DENSITY_THRESHOLD
        else ["Well colour: % of the plate's negative-control median (viridis; dark = low)"]
    )
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
        title="Z′ per plate (terracotta: Z′ < 0.5)",
    )
    fail = t["fails"].to_numpy()
    passed = ~fail & np.isfinite(z)
    if n <= 40:
        ok_w = 0.8
        if passed.any():
            zp.add(
                RectMark(
                    x0=x[passed] - ok_w / 2,
                    y0=np.full(passed.sum(), zlo),
                    width=ok_w,
                    height=zc[passed] - zlo,
                    color="ink_soft",
                    legend="Z′ ≥ 0.5",
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
                    categorical=True,
                    legend="failing (Z′ < 0.5)",
                )
            )
    else:
        if passed.any():
            zp.add(
                PointMark(
                    x=x[passed],
                    y=zc[passed],
                    color="ink_soft",
                    size=3,
                    alpha=0.6,
                    legend="Z′ ≥ 0.5",
                )
            )
        if fail.any():
            zp.add(
                PointMark(
                    x=x[fail],
                    y=zc[fail],
                    color="focus",
                    size=5,
                    categorical=True,
                    legend="failing (Z′ < 0.5)",
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
        y=Scale("linear", label="Outer seeded ring / interior"),
        row=1,
        col=1,
        title="Edge effect",
    )
    er = t["edge_ratio"].to_numpy(float)
    ok = np.isfinite(er)
    if ok.any():
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
            label=f"Neg. control / median\n({qc.drift_basis})",
            ticks=[0.25, 0.5, 1.0, 2.0, 4.0],
            tick_labels=["0.25", "0.5", "1", "2", "4"],
        ),
        row=2,
        col=1,
        title="Control drift",
    )
    lv = t["neg_level"].to_numpy(float)
    okl = np.isfinite(lv)
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
        "{n} plates ({w}-well); median Z′ {med:.2f} (IQR {q1:.2f}–{q3:.2f}); "
        "{nf} ({pct:.0%}) fall below {thr:g}.",
        n=n,
        w=pt.n_wells,
        med=float(zf.median()) if len(zf) else None,
        q1=float(zf.quantile(0.25)) if len(zf) else None,
        q3=float(zf.quantile(0.75)) if len(zf) else None,
        nf=n_fail,
        pct=n_fail / max(len(zf), 1),
        thr=qc.threshold,
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
            "{k} of {n} plates leave the outer ring unseeded; their edge is the next seeded ring.",
            k=qc.n_unseeded_outer,
            n=n,
        )
    cap.add(
        "Negative-control level (relative to {basis}) changes {s:+.0%} across the {n} plates in "
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
            *colour_note,
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
