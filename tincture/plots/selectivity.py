"""Selectivity Spectrum: how specific is this compound? (HANDOVER §5.8)

Potency strip across a target panel, strongest first (pX = −log10 molar). The intended
target(s) are terracotta; off-targets engaged at exposure (above the free-Cmax band) are
peach. Without exposure data there is no band, and off-targets are highlighted only if
within 10-fold of the intended target, and the caption says so. Censored values ("not
detected up to X") are hollow markers at their limit, never measurements. Panel coverage
("n of ~518 kinases tested") is printed in the footer.

Guard (tested): mixed assay types (e.g. Kd and IC50) are refused unless
``allow_mixed=True``, in which case each type gets its own marker.
Classic fallback: bar chart of the most potent targets.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from tincture.errors import TinctureDataError
from tincture.model.potency import PotencyTable
from tincture.report.captions import Caption, fmt_conc
from tincture.report.classic import build_classic, register_classic
from tincture.report.methods import Methods, register_block
from tincture.result import PlotResult
from tincture.spec import (
    BandMark,
    LegendSpec,
    Panel,
    PlotSpec,
    PointMark,
    RectMark,
    Scale,
    TextMark,
)

FOLD_RULE = 10.0
MARKERS = ["o", "s", "^", "D", "v", "P"]

register_block(
    "selectivity",
    "Potencies are shown as pX = −log10(molar) for {n} tested targets ({assays}), sorted "
    "from strongest. {rule} Censored values (no binding / no activity up to the highest "
    "concentration) are drawn at their limit as hollow markers and are not ranked as "
    "measurements.",
)


def _molar(px: float) -> float:
    return float(10.0 ** (-px))


def build_spec(
    d: pd.DataFrame,
    compound: str,
    intended: Sequence[str],
    cmax: tuple[float, float] | None,
    allow_mixed: bool,
    panel_size: int | None,
    panel_name: str,
    top_labels: int,
    title: str,
    footer: Sequence[str],
) -> PlotSpec:
    meas = d[d["censored"] == ""].sort_values("px", ascending=False)
    cens = d[d["censored"] == ">"]
    ints = [t for t in intended if t in set(d["target"])]
    if not ints:
        raise TinctureDataError(
            f"Intended target(s) {list(intended)} not in the panel for {compound}.",
            fix="Check target names against the table's `target` column.",
        )
    ranked = pd.concat([meas, cens.sort_values("px", ascending=False)])
    ranked = ranked.assign(rank=np.arange(1, len(ranked) + 1))
    is_int = ranked["target"].isin(ints).to_numpy()
    is_meas = (ranked["censored"] == "").to_numpy()
    px_int = (
        float(ranked.loc[ranked["target"].isin(ints) & (ranked["censored"] == ""), "px"].max())
        if (is_int & is_meas).any()
        else float("nan")
    )
    cap = Caption()
    if cmax is not None:
        engaged_px = -np.log10(cmax[0])  # above the band = potent even at the low exposure
        off = ~is_int & is_meas & (ranked["px"].to_numpy() > engaged_px)
        rule = (
            f"Off-targets in peach have potency stronger than the low end of the free Cmax "
            f"range ({fmt_conc(cmax[0])}–{fmt_conc(cmax[1])})."
        )
    else:
        off = ~is_int & is_meas & (ranked["px"].to_numpy() >= px_int - np.log10(FOLD_RULE))
        rule = (
            "No exposure (free Cmax) was supplied, so off-targets are highlighted if within "
            f"{FOLD_RULE:.0f}-fold of the intended target's potency."
        )
    lo_px = float(ranked["px"].min()) - 0.4
    hi_px = float(ranked["px"].max()) + 0.6
    p = Panel(
        key="spectrum",
        x=Scale("linear", label=f"{panel_name} ranked by potency", domain=(0, len(ranked) + 1)),
        y=Scale("linear", label="potency pX = −log10(M)", domain=(lo_px, hi_px)),
        legend=LegendSpec(loc="upper right", title="assay") if allow_mixed else None,
    )
    if cmax is not None:
        p.add(
            BandMark(
                lo=-np.log10(cmax[1]),
                hi=-np.log10(cmax[0]),
                span=True,
                orientation="h",
                color="neutral",
                alpha=0.45,
                inner_label=f"free Cmax {fmt_conc(cmax[0])}–{fmt_conc(cmax[1])}",
                z=0.2,
            )
        )
    types = list(dict.fromkeys(ranked["assay_type"]))
    for k, at in enumerate(types):
        m_t = (ranked["assay_type"] == at).to_numpy()
        mk = MARKERS[k % len(MARKERS)] if allow_mixed else "o"
        groups = [
            (m_t & is_meas & ~is_int & ~off, "mid_grey", False, False),
            (m_t & off, "focus_light", True, False),
            (m_t & is_int & is_meas, "focus", True, False),
            (m_t & ~is_meas, "mid_grey", False, True),
        ]
        first = True
        for sel, col, cat, hollow in groups:
            if not sel.any():
                continue
            p.add(
                PointMark(
                    x=ranked["rank"].to_numpy()[sel],
                    y=ranked["px"].to_numpy()[sel],
                    color=col,
                    size=26 if cat else 12,
                    marker=mk,
                    hollow=hollow,
                    categorical=cat,
                    legend=at if (allow_mixed and first) else None,
                    z=3 if cat else 2,
                )
            )
            first = False
    # labels: intended + strongest off-targets
    lab_rows = ranked[is_int | off].head(top_labels + len(ints))
    for _, r in lab_rows.iterrows():
        col = "focus" if r["target"] in ints else "ink_soft"
        p.add(
            TextMark(
                x=float(r["rank"]) + 0.6,
                y=float(r["px"]),
                text=str(r["target"]),
                ha="left",
                va="center",
                color=col,
                size=6.5,
                categorical=r["target"] in ints,
            )
        )
    if len(cens):
        lim = float(cens["px"].max())
        p.add(
            TextMark(
                x=float(len(meas)) + 0.5,
                y=lim + 0.25,
                text=f"{len(cens)} not detected (weaker than {fmt_conc(_molar(lim))})",
                ha="left",
                va="bottom",
                color="ink_soft",
                size=6.5,
            )
        )
    if np.isfinite(px_int):
        p.add(
            RectMark(
                x0=[0],
                y0=[px_int - 0.005],
                width=[len(ranked) + 1],
                height=[0.01],
                color="focus_light",
                alpha=0.6,
            )
        )
    n_off = int(off.sum())
    n_meas = int(is_meas.sum())
    int_txt = ", ".join(ints)
    if np.isfinite(px_int):
        rank_int = int(ranked.loc[ranked["target"].isin(ints), "rank"].min())
        cap.add(
            "{c} binds its intended target {t} at {pot} (rank {r} of {n} tested). {k} "
            "off-target(s) are highlighted ({s:.1%} of the panel){strong}.",
            c=compound,
            t=int_txt,
            pot=fmt_conc(_molar(px_int)),
            r=rank_int,
            n=len(ranked),
            k=n_off,
            s=n_off / max(len(ranked), 1),
            strong=(
                "; the strongest is "
                + str(ranked[off].iloc[0]["target"])
                + " at "
                + fmt_conc(_molar(float(ranked[off].iloc[0]["px"])))
            )
            if n_off
            else "",
        )
    else:
        cap.weak(f"The intended target {int_txt} was not detected in this panel.")
    cap.add("{rule}", rule=rule)
    if allow_mixed and len(types) > 1:
        cap.weak(
            "assay types are mixed ({t}); potencies of different types are not "
            "directly comparable.".format(t=", ".join(types))
        )
    cov = (
        f"{len(ranked)} of ~{panel_size} {panel_name} tested"
        if panel_size
        else f"{len(ranked)} {panel_name} tested"
    )
    meth = Methods().use("selectivity", n=len(ranked), assays=", ".join(types), rule=rule)
    _ = n_meas
    return PlotSpec(
        plot="selectivity",
        panels=[p],
        title=title,
        caption=cap.render(),
        methods=meth.render(),
        footer=[f"Panel coverage: {cov}", *footer],
        warnings=list(cap.notes),
        aspect=0.55,
        meta={"ranked": ranked, "off_targets": list(ranked[off]["target"]), "rule": rule},
    ).validate()


@register_classic("selectivity", "bar chart of most potent targets")
def build_classic_spec(ranked: pd.DataFrame | None = None, n: int = 25, **_: object) -> PlotSpec:
    if ranked is None:
        raise TinctureDataError("Classic needs the ranked table.")
    r = ranked[ranked["censored"] == ""].head(n)
    p = Panel(
        x=Scale(
            "linear",
            ticks=list(range(len(r))),
            tick_labels=list(r["target"]),
            tick_rotation=60,
            label="",
        ),
        y=Scale(label="pX"),
    )
    p.add(
        RectMark(
            x0=np.arange(len(r)) - 0.4,
            y0=np.zeros(len(r)),
            width=np.full(len(r), 0.8),
            height=r["px"].to_numpy(),
            color="mid_grey",
        )
    )
    return PlotSpec(
        plot="selectivity_classic",
        panels=[p],
        title="Most potent targets",
        caption="Most potent targets (bar chart).",
        aspect=0.5,
    ).validate()


def plot(
    table: PotencyTable,
    compound: str,
    intended: str | Sequence[str],
    cmax: tuple[float, float] | None = None,
    allow_mixed: bool = False,
    panel_size: int | None = None,
    panel_name: str = "targets",
    top_labels: int = 8,
    title: str = "Selectivity Spectrum",
    footer: Sequence[str] = (),
) -> PlotResult:
    """Selectivity Spectrum for one compound. ``cmax``: free Cmax range (molar); never pass
    a guess. ``panel_size``: size of the full family (e.g. ~518 human kinases) for the
    coverage footer."""
    d = table.compound(compound)
    types = sorted(d["assay_type"].unique())
    if len(types) > 1 and not allow_mixed:
        raise TinctureDataError(
            f"Mixed assay types for {compound}: {types}. Kd, Ki and IC50 are not interchangeable.",
            fix="Filter to one assay type, or pass allow_mixed=True (markers then differ by "
            "type and the caption flags it).",
        )
    if cmax is not None and not (0 < cmax[0] <= cmax[1] < 1):
        raise TinctureDataError("cmax must be a (low, high) range in molar.")
    ints = [intended] if isinstance(intended, str) else list(intended)
    spec = build_spec(
        d, compound, ints, cmax, allow_mixed, panel_size, panel_name, top_labels, title, footer
    )
    return PlotResult(spec=spec, classic=build_classic("selectivity", ranked=spec.meta["ranked"]))
