"""A minimal plot exercising the whole M0 pipeline (spec → render → caption → export).

It is *not* one of the question-first plots; it exists so the scaffold can be tested
end to end. Input is the canonical long table (section 4.1): ``compound_id``,
``dose`` (molar, 0 = vehicle), ``replicate``, ``value``.

Left panel: mean response vs dose per compound with 95% t-intervals, vehicle mean as
a reference line. Right panel: compound × dose heatmap of the change from vehicle on
the diverging scale (clip limit stated on the colourbar).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from tincture.errors import TinctureDataError
from tincture.report.captions import Caption, fmt_conc
from tincture.report.classic import build_classic, register_classic
from tincture.report.methods import Methods
from tincture.result import PlotResult
from tincture.spec import (
    ColorbarSpec,
    HeatmapMark,
    IntervalMark,
    LegendSpec,
    LineMark,
    NoticeMark,
    Panel,
    PlotSpec,
    PointMark,
    RefLineMark,
    Scale,
)

REQUIRED = ("compound_id", "dose", "replicate", "value")


def _check(df: pd.DataFrame) -> pd.DataFrame:
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise TinctureDataError(
            f"Missing required column(s) {missing}.",
            fix="Provide a long table with compound_id, dose (molar; 0 = vehicle), "
            "replicate and value.",
        )
    d = df.loc[:, list(REQUIRED)].copy()
    d["dose"] = pd.to_numeric(d["dose"], errors="coerce")
    d["value"] = pd.to_numeric(d["value"], errors="coerce")
    d = d.dropna(subset=["dose", "value"])
    if (d["dose"] < 0).any():
        raise TinctureDataError("Negative doses found.", fix="Doses must be molar and ≥ 0.")
    if not (d["dose"] == 0).any():
        raise TinctureDataError(
            "No vehicle rows (dose == 0).",
            fix="Include vehicle wells with dose 0; every plot needs a reference.",
        )
    if d.loc[d["dose"] > 0, "compound_id"].nunique() == 0:
        raise TinctureDataError("No treated rows (dose > 0).", fix="Add compound-treated rows.")
    return d


def _summary(d: pd.DataFrame) -> pd.DataFrame:
    g = d.loc[d["dose"] > 0].groupby(["compound_id", "dose"])["value"]
    s = g.agg(["mean", "std", "count"]).reset_index()
    tcrit = stats.t.ppf(0.975, np.maximum(s["count"] - 1, 1))
    half = np.where(s["count"] > 1, tcrit * s["std"] / np.sqrt(s["count"]), np.nan)
    s["lo"] = s["mean"] - half
    s["hi"] = s["mean"] + half
    return s


def build_spec(data: pd.DataFrame, focus: str | None = None, clip: float = 1.0) -> PlotSpec:
    d = _check(data)
    vehicle = float(d.loc[d["dose"] == 0, "value"].mean())
    s = _summary(d)
    compounds = sorted(s["compound_id"].unique())
    if focus is not None and focus not in compounds:
        raise TinctureDataError(
            f"Focus compound {focus!r} not in data.", fix=f"Use one of {compounds}."
        )
    focus = focus or compounds[0]
    doses = np.sort(s["dose"].unique())

    left = Panel(
        key="curves",
        x=Scale("log", label="Concentration (M)"),
        y=Scale("linear", label="Response"),
        legend=LegendSpec(loc="lower left"),
        title="Response vs dose",
    )
    left.add(RefLineMark(axis="y", value=vehicle, label="vehicle", z=1))
    others = [c for c in compounds if c != focus]

    def legend_for(c: str) -> str | None:
        if c == focus:
            return f"{focus} (focus)"
        if others and c == others[0]:
            return "other compounds" if len(others) > 1 else c
        return None

    for c in compounds:
        sc = s.loc[s["compound_id"] == c].sort_values("dose")
        col = "focus" if c == focus else "reference"
        left.add(
            LineMark(
                x=sc["dose"].to_numpy(),
                y=sc["mean"].to_numpy(),
                color=col,
                categorical=True,
                width=1.2 if c == focus else 0.8,
            ),
            IntervalMark(
                pos=sc["dose"].to_numpy(),
                lo=sc["lo"].to_numpy(),
                hi=sc["hi"].to_numpy(),
                color=col,
                categorical=True,
            ),
            PointMark(
                x=sc["dose"].to_numpy(),
                y=sc["mean"].to_numpy(),
                color=col,
                size=14,
                categorical=True,
                legend=legend_for(c),
            ),
        )

    single_rep = bool((s["count"] < 2).all())
    cap = Caption()
    if single_rep:
        left.add(NoticeMark(text="single replicate: no CI", loc="upper right"))
        cap.weak("only one replicate per dose, so no confidence intervals are shown.")

    mat = (
        s.pivot(index="compound_id", columns="dose", values="mean").reindex(
            index=compounds, columns=doses
        )
        - vehicle
    )
    right = Panel(
        key="heatmap",
        x=Scale("categorical", label="Concentration", categories=[fmt_conc(x) for x in doses]),
        y=Scale("categorical", categories=compounds),
        colorbar=ColorbarSpec(label="Δ vs vehicle", clip_note=f"clipped at ±{clip:g}"),
        title="Change from vehicle",
    )
    right.add(HeatmapMark(values=mat.to_numpy(), cmap="diverging", clip=clip))

    fs = s.loc[s["compound_id"] == focus]
    top = fs.loc[fs["dose"].idxmax()]
    cap.add(
        "{n} compounds at {k} concentrations; {focus} changes the response by {delta:+.2f} "
        "at {dose} relative to vehicle.",
        n=len(compounds),
        k=len(doses),
        focus=focus,
        delta=float(top["mean"] - vehicle),
        dose=fmt_conc(float(top["dose"])),
    )
    meth = Methods().text(
        "Points are means over replicates; whiskers are two-sided 95% t-intervals."
    )
    meth.use("colour")

    return PlotSpec(
        plot="demo",
        panels=[left, right],
        ncols=2,
        width_ratios=[1.3, 1.0],
        title="tincture demo",
        caption=cap.render(),
        methods=meth.render(),
        footer=[f"Reference: vehicle (n = {int((d['dose'] == 0).sum())} wells)"],
        warnings=list(cap.notes),
        meta={"vehicle_mean": vehicle, "focus": focus, "single_replicate": single_rep},
        aspect=0.45,
    ).validate()


@register_classic("demo", "scatter of raw values")
def build_classic_spec(data: pd.DataFrame, **_: object) -> PlotSpec:
    d = _check(data)
    d = d.loc[d["dose"] > 0]
    p = Panel(x=Scale("log", label="Concentration (M)"), y=Scale("linear", label="Response"))
    p.add(PointMark(x=d["dose"].to_numpy(), y=d["value"].to_numpy(), color="ink_soft", size=6))
    return PlotSpec(plot="demo_classic", panels=[p], title="Raw values", aspect=0.6).validate()


def demo(data: pd.DataFrame, focus: str | None = None, clip: float = 1.0) -> PlotResult:
    """Build the demo plot. Returns a :class:`~tincture.result.PlotResult`."""
    return PlotResult(
        spec=build_spec(data, focus=focus, clip=clip),
        classic=build_classic("demo", data=data),
    )
