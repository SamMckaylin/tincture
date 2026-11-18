"""Combination Surface and Combination Lift (HANDOVER §5.11).

The question both views answer: *is there a dose pair, reachable in patients, where the
combination clearly beats what the drugs do on their own?* Synergy and benefit are reported
separately, because a combination can be synergistic only where both drugs barely work.

**Combination Surface** (one cell line, full matrix such as 7 × 7):
* left: filled contours of the **observed** viability, interpolated only inside the
  measured grid (no extrapolation); grid points as dots; the significant-excess region
  outlined in terracotta; ceiling-masked cells hatched; the clinically reachable box (free
  Cmax of each drug) dashed when supplied; single-agent curves in the margins.
* right: excess over the reference model per cell (diverging; values printed), significant
  cells outlined, masked cells hatched.

Statistics in :mod:`tincture.stats.synergy`: Bliss (default), ZIP or HSA; the caption also
reports benefit over the best single agent (HSA) whatever the headline model.

Guards (each tested): the model name is printed; masked ceiling cells are hatched and never
tested; an anchored matrix is refused with a pointer to the Lift; a single replicate without
a noise estimate draws no significant region and the caption says so.
Classic fallback: synergy heatmap.

**Combination Lift** (many cell lines; works with anchored 2 × 7 designs):
x = effect of the best single agent, y = effect of the combination (1 − viability at the
highest library concentration, as Emax); the diagonal is "no benefit" and the band above it
is the benefit threshold (default ΔEmax ≥ 0.2 over the best single agent). Points are cell
lines coloured by a group (tissue or biomarker); failed fits are hollow and not counted.
Guards (tested): threshold and model printed; failed fits hollow.
Classic fallback: waterfall of ΔEmax.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

import numpy as np
import pandas as pd
from scipy import stats as sps

from tincture.errors import TinctureDataError
from tincture.model.combo import ComboMatrix
from tincture.report.captions import Caption, fmt_conc, fmt_p
from tincture.report.classic import build_classic, register_classic
from tincture.report.methods import Methods, register_block
from tincture.result import PlotResult
from tincture.spec import (
    ColorbarSpec,
    ContourMark,
    HeatmapMark,
    LegendSpec,
    LineMark,
    Panel,
    PlotSpec,
    PointMark,
    PolygonMark,
    RectMark,
    Scale,
    TextMark,
)
from tincture.stats.synergy import CEILING, MIN_EXCESS, MODEL_NAMES, SynergyResult, synergy

LIFT_THRESHOLD = 0.2
EXCESS_CLIP = 0.3
GROUP_COLOURS = ("focus", "reference", "novel", "good", "mid_grey")

register_block(
    "combination_surface",
    "Viability was normalised to vehicle. Expected viability under {model_name} was computed "
    "from {single} single-agent responses at the same concentrations; excess = expected − "
    "observed viability (positive = more inhibition than expected). Cells with expected "
    "inhibition above {ceiling:.0%} were masked (no room for excess). {uncertainty} Cells "
    "with BH q < {alpha} and excess ≥ {min_excess:.0%} form the significant region. The "
    "observed surface is interpolated only within the measured {na} × {nb} grid.",
)
register_block(
    "combination_lift",
    "For each cell line, the effect of each drug alone and of the combination is taken at "
    "the highest tested library concentration (Emax, as 1 − viability). Benefit = combination "
    "effect − effect of the best single agent; a cell line benefits if this is ≥ "
    "{threshold:.0%}. {model_note} Lines with failed curve fits are shown hollow and excluded "
    "from rates.",
)


# ---------------------------------------------------------------------------------------
# Combination Surface
# ---------------------------------------------------------------------------------------
def _edges(levels: np.ndarray) -> np.ndarray:
    """Cell edges in concentration (geometric midpoints; half a step beyond the ends)."""
    lx = np.log10(levels)
    step = np.diff(lx)
    mid = (lx[:-1] + lx[1:]) / 2
    e = np.r_[lx[0] - step[0] / 2, mid, lx[-1] + step[-1] / 2]
    return np.asarray(10.0**e, float)


def _region_outline(
    mask: np.ndarray, ex: np.ndarray, ey: np.ndarray, color: str, width: float
) -> list[LineMark]:
    """Boundary segments of the union of True cells (exact, no smoothing)."""
    nb, na = mask.shape
    out: list[LineMark] = []
    for j in range(nb):
        for i in range(na):
            if not mask[j, i]:
                continue
            x0, x1, y0, y1 = ex[i], ex[i + 1], ey[j], ey[j + 1]
            if i == 0 or not mask[j, i - 1]:
                out.append(LineMark(x=[x0, x0], y=[y0, y1], color=color, width=width, z=4))
            if i == na - 1 or not mask[j, i + 1]:
                out.append(LineMark(x=[x1, x1], y=[y0, y1], color=color, width=width, z=4))
            if j == 0 or not mask[j - 1, i]:
                out.append(LineMark(x=[x0, x1], y=[y0, y0], color=color, width=width, z=4))
            if j == nb - 1 or not mask[j + 1, i]:
                out.append(LineMark(x=[x0, x1], y=[y1, y1], color=color, width=width, z=4))
    return out


def _uncertainty_text(r: SynergyResult) -> str:
    if r.uncertainty == "replicates":
        return (
            f"Per-cell intervals and p-values come from a t-test across {r.n_replicates} "
            "replicate matrices (excess computed within each replicate)."
        )
    if r.uncertainty == "noise_sd":
        return (
            "Only one replicate matrix: intervals come from a parametric bootstrap with a "
            f"well-level noise SD of {r.noise_sd:.3g} viability ({r.noise_method})."
        )
    return "No replicate matrices and no noise estimate: no significance was assessed."


def _surface_spec(
    cm: ComboMatrix,
    r: SynergyResult,
    hsa: SynergyResult,
    cmax_a: tuple[float, float] | None,
    cmax_b: tuple[float, float] | None,
    title: str,
    footer: Sequence[str],
) -> PlotSpec:
    la, lb = r.levels_a, r.levels_b
    ex, ey = _edges(la), _edges(lb)
    name = MODEL_NAMES[r.model]
    short = {"bliss": "Bliss", "hsa": "HSA", "zip": "ZIP"}[r.model]
    xs = Scale("log", label=f"{cm.drug_a} (M)", domain=(ex[0], ex[-1]))
    ys = Scale("log", label=f"{cm.drug_b} (M)", domain=(ey[0], ey[-1]))
    main = Panel(
        key="surface",
        x=xs,
        y=ys,
        row=1,
        col=0,
        colorbar=ColorbarSpec(label="observed viability", orientation="horizontal"),
    )
    main.add(
        ContourMark(
            x=la, y=lb, z_values=np.clip(r.observed, 0, 1.2), levels=np.linspace(0, 1.2, 13),
            cmap="cividis", z=0.5,
        )
    )  # fmt: skip
    ga, gb = np.meshgrid(la, lb)
    main.add(PointMark(x=ga.ravel(), y=gb.ravel(), color="paper", size=5, edgecolor="ink", z=3))
    # ceiling mask: hatched cells
    if r.masked.any():
        jj, ii = np.nonzero(r.masked)
        gx = np.clip(ex, la[0], la[-1])
        gy = np.clip(ey, lb[0], lb[-1])
        for p in (main,):
            p.add(
                RectMark(
                    x0=gx[ii], y0=gy[jj], width=gx[ii + 1] - gx[ii], height=gy[jj + 1] - gy[jj],
                    color="none", edgecolor="paper", hatch="////", linewidth=0, alpha=0.9, z=2,
                )
            )  # fmt: skip
    for seg in _region_outline(r.significant, ex, ey, "focus", 1.6):
        main.add(seg)
    in_box = None
    if cmax_a is not None and cmax_b is not None:
        main.add(
            PolygonMark(
                x=[cmax_a[0], cmax_a[1], cmax_a[1], cmax_a[0]],
                y=[cmax_b[0], cmax_b[0], cmax_b[1], cmax_b[1]],
                color="none", edgecolor="ink", alpha=1.0, z=4.5,
            )
        )  # fmt: skip
        main.add(
            TextMark(
                x=cmax_a[0], y=cmax_b[1] * 1.08, text="free Cmax", ha="left", va="bottom",
                color="ink", size=6.5,
            )
        )  # fmt: skip
        reach = (ga <= cmax_a[1]) & (gb <= cmax_b[1])
        in_box = int((r.significant & reach).sum())
    # margins: single agents
    top = Panel(
        key="margin_a",
        x=Scale("log", label="", domain=xs.domain, visible=False),
        y=Scale("linear", label="viability", domain=(0, 1.2), ticks=[0, 0.5, 1.0]),
        row=0,
        col=0,
    )
    right = Panel(
        key="margin_b",
        x=Scale("linear", label="viability", domain=(0, 1.2), ticks=[0, 0.5, 1.0]),
        y=Scale("log", label="", domain=ys.domain, visible=False),
        row=1,
        col=1,
    )
    for drug, panel, levels, fit in (("a", top, la, r.fit_a), ("b", right, lb, r.fit_b)):
        s = cm.single(drug)  # type: ignore[arg-type]
        m = s.groupby("conc")["viability"].mean()
        xx, yy = m.index.to_numpy(float), m.to_numpy(float)
        grid = np.geomspace(levels[0], levels[-1], 80)
        if drug == "a":
            panel.add(PointMark(x=xx, y=yy, color="mid_grey", size=8, z=3))
            if fit is not None:
                panel.add(LineMark(x=grid, y=fit(grid), color="ink", width=0.9))
        else:
            panel.add(PointMark(x=yy, y=xx, color="mid_grey", size=8, z=3))
            if fit is not None:
                panel.add(LineMark(x=fit(grid), y=grid, color="ink", width=0.9))
    # excess heatmap (rows: highest drug-B concentration at the top)
    na_, nb_ = len(la), len(lb)
    vals = np.flipud(np.clip(r.excess, -EXCESS_CLIP, EXCESS_CLIP))
    heat = Panel(
        key="excess",
        x=Scale(
            "categorical", label=f"{cm.drug_a}", categories=[fmt_conc(v) for v in la],
            tick_rotation=60,
        ),
        y=Scale("categorical", label=f"{cm.drug_b}", categories=[fmt_conc(v) for v in lb[::-1]]),
        row=1,
        col=2,
        title=f"excess over {short}",
        colorbar=ColorbarSpec(
            label=f"excess over {short} (viability)", clip_note=f"clipped at ±{EXCESS_CLIP:g}",
            compact=True,
        ),
    )  # fmt: skip
    heat.add(
        HeatmapMark(
            values=vals, mask=np.flipud(r.masked | ~np.isfinite(r.excess)), cmap="diverging",
            clip=EXCESS_CLIP, style="flat",
        )
    )  # fmt: skip
    fm = np.flipud(r.masked)
    fs = np.flipud(r.significant)
    fe = np.flipud(r.excess)
    for jj in range(nb_):
        for ii in range(na_):
            if fm[jj, ii]:
                heat.add(
                    RectMark(
                        x0=[ii - 0.5], y0=[jj - 0.5], width=[1], height=[1], color="none",
                        edgecolor="mid_grey", hatch="////", linewidth=0, z=2,
                    )
                )  # fmt: skip
                continue
            if np.isfinite(fe[jj, ii]):
                heat.add(
                    TextMark(
                        x=ii, y=jj, text=f"{round(100 * fe[jj, ii]) + 0:.0f}",
                        ha="center", va="center",
                        size=6, color="ink",
                    )
                )  # fmt: skip
            if fs[jj, ii]:
                heat.add(
                    RectMark(
                        x0=[ii - 0.47], y0=[jj - 0.47], width=[0.94], height=[0.94],
                        color="none", edgecolor="ink", linewidth=1.2, z=4,
                    )
                )  # fmt: skip
    # caption
    cap = Caption()
    testable = int((~r.masked & np.isfinite(r.excess)).sum())
    cap.add(
        "{a} + {b} in {line}: observed viability reaches {vmin:.0%} at the highest doses "
        "({n_rep} replicate {es}).",
        a=cm.drug_a, b=cm.drug_b, line=cm.cell_line or "this cell line",
        vmin=max(0.0, float(np.nanmin(r.observed))), n_rep=r.n_replicates,
        es="matrix" if r.n_replicates == 1 else "matrices",
    )  # fmt: skip
    if r.uncertainty == "none":
        cap.add(
            "Under {name}, the mean excess is {score:+.0%} (largest {mx:+.0%}).",
            name=name, score=r.score, mx=float(np.nanmax(np.where(r.masked, np.nan, r.excess))),
        )  # fmt: skip
        cap.weak(
            "single replicate and no noise estimate: excess is shown without significance, "
            "so no region is outlined."
        )
    else:
        if r.n_significant:
            k = np.nanargmax(np.where(r.significant, r.excess, -np.inf))
            j, i = np.unravel_index(k, r.excess.shape)
            cap.add(
                "Under {name}, {k} of {n} testable cells show significant excess ≥ {mn:.0%} "
                "(BH q < {al}); the largest is {mx:+.0%} at {ca} {a} + {cb} {b}.",
                name=name, k=r.n_significant, n=testable, mn=r.min_excess, al=r.alpha,
                mx=float(r.excess[j, i]), ca=fmt_conc(la[i]), a=cm.drug_a,
                cb=fmt_conc(lb[j]), b=cm.drug_b,
            )  # fmt: skip
        else:
            cap.add(
                "Under {name}, no cell shows significant excess ≥ {mn:.0%} (BH q < {al}; mean "
                "excess {score:+.0%}).",
                name=name, mn=r.min_excess, al=r.alpha, score=r.score,
            )  # fmt: skip
        if r.model != "hsa":
            cap.add(
                "Benefit over the best single agent (HSA) ≥ {mn:.0%} is significant in {k} "
                "cell(s).",
                mn=hsa.min_excess, k=hsa.n_significant,
            )  # fmt: skip
        if r.n_replicates in (2, 3) and r.uncertainty == "replicates":
            cap.weak(f"only {r.n_replicates} replicates: intervals are wide.")
    if r.masked.any():
        cap.add(
            "{m} cell(s) are hatched: expected inhibition exceeds {c:.0%}, so excess there is "
            "not interpretable.",
            m=int(r.masked.sum()), c=CEILING,
        )  # fmt: skip
    if in_box is not None:
        cap.add(
            "{k} significant cell(s) lie inside the clinically reachable box (free Cmax).",
            k=in_box,
        )
    else:
        cap.add("No exposure (free Cmax) was supplied, so clinical reachability is not assessed.")
    for n in r.notes:
        cap.weak(n)
    meth = Methods().use(
        "combination_surface",
        model_name=name,
        single="fitted (4PL)" if r.single_agent == "fitted" else "observed mean",
        ceiling=CEILING,
        uncertainty=_uncertainty_text(r),
        alpha=r.alpha,
        min_excess=r.min_excess,
        na=len(la),
        nb=len(lb),
    )
    return PlotSpec(
        plot="combination_surface",
        panels=[top, main, right, heat],
        nrows=2,
        ncols=3,
        width_ratios=[1.0, 0.24, 0.95],
        height_ratios=[0.24, 1.0],
        title=title,
        caption=cap.render(),
        methods=meth.render(),
        footer=[
            f"Model: {name}; single agents: {r.single_agent}; replicates: {r.n_replicates}",
            *footer,
        ],
        warnings=list(cap.notes),
        aspect=0.5,
        meta={"result": r, "hsa": hsa, "significant_in_box": in_box},
    ).validate()


@register_classic("combination_surface", "synergy heatmap")
def build_surface_classic(result: SynergyResult | None = None, **_: object) -> PlotSpec:
    if result is None:
        raise TinctureDataError("Classic needs the synergy result.")
    p = Panel(
        x=Scale("categorical", label="drug A", categories=[fmt_conc(v) for v in result.levels_a]),
        y=Scale(
            "categorical",
            label="drug B",
            categories=[fmt_conc(v) for v in result.levels_b[::-1]],
        ),
        colorbar=ColorbarSpec(label="synergy score", clip_note=f"clipped at ±{EXCESS_CLIP:g}"),
    )
    p.add(
        HeatmapMark(
            values=np.flipud(np.clip(result.excess, -EXCESS_CLIP, EXCESS_CLIP)),
            cmap="diverging",
            clip=EXCESS_CLIP,
            style="flat",
        )
    )
    return PlotSpec(
        plot="combination_surface_classic",
        panels=[p],
        title=f"{MODEL_NAMES[result.model]} synergy",
        caption="Synergy heatmap.",
        aspect=0.8,
    ).validate()


def surface(
    cm: ComboMatrix,
    model: Literal["bliss", "zip", "hsa"] = "bliss",
    single_agent: Literal["fitted", "observed"] = "fitted",
    noise_sd: float | Literal["matrix"] | None = None,
    cmax_a: tuple[float, float] | None = None,
    cmax_b: tuple[float, float] | None = None,
    min_excess: float = MIN_EXCESS,
    alpha: float = 0.05,
    n_boot: int = 200,
    title: str = "Combination Surface",
    footer: Sequence[str] = (),
) -> PlotResult:
    """Combination Surface for one full concentration matrix. ``cmax_a``/``cmax_b``: free
    Cmax ranges (molar) of each drug; never pass a guess. ``noise_sd``: well-level SD of
    viability, used only when there is a single replicate matrix; ``"matrix"`` estimates it
    from the matrix itself (``stats.synergy.matrix_noise_sd``)."""
    for c in (cmax_a, cmax_b):
        if c is not None and not (0 < c[0] <= c[1] < 1):
            raise TinctureDataError("cmax ranges must be (low, high) in molar.")
    kw = {"single_agent": single_agent, "noise_sd": noise_sd, "alpha": alpha,
          "min_excess": min_excess, "n_boot": n_boot}  # fmt: skip
    r = synergy(cm, model, **kw)  # type: ignore[arg-type]
    hsa = r if model == "hsa" else synergy(cm, "hsa", **kw)  # type: ignore[arg-type]
    spec = _surface_spec(cm, r, hsa, cmax_a, cmax_b, title, footer)
    return PlotResult(spec=spec, classic=build_classic("combination_surface", result=r))


# ---------------------------------------------------------------------------------------
# Combination Lift
# ---------------------------------------------------------------------------------------
LIFT_COLUMNS = ("cell_line", "single_a", "single_b", "combo")


def lift_table(d: pd.DataFrame, group: str | None = None, facet: str | None = None) -> pd.DataFrame:
    """Per cell line (and facet): mean viabilities over replicates, ΔEmax with a t-interval
    when there are ≥ 2 replicates. ``fit_ok`` (optional, bool) marks usable fits."""
    miss = [c for c in LIFT_COLUMNS if c not in d.columns]
    if miss:
        raise TinctureDataError(
            f"Lift table missing columns {miss}.",
            fix="Provide cell_line, single_a, single_b, combo (viabilities at Emax).",
        )
    x = d.copy()
    if "fit_ok" not in x:
        x["fit_ok"] = True
    x["fit_ok"] = x["fit_ok"].astype(bool) & np.isfinite(
        x[["single_a", "single_b", "combo"]].to_numpy(float)
    ).all(axis=1)
    x["delta"] = np.minimum(x["single_a"], x["single_b"]) - x["combo"]
    keys = ["cell_line"] + [k for k in (facet,) if k] + [k for k in (group,) if k]
    ok = x[x["fit_ok"]]

    def agg(g: pd.DataFrame) -> pd.Series:
        n = len(g)
        m = g["delta"].mean()
        se = g["delta"].std(ddof=1) / np.sqrt(n) if n >= 2 else np.nan
        tc = sps.t.ppf(0.975, n - 1) if n >= 2 else np.nan
        return pd.Series(
            {
                "single_a": g["single_a"].mean(),
                "single_b": g["single_b"].mean(),
                "combo": g["combo"].mean(),
                "delta": m,
                "delta_lo": m - tc * se,
                "delta_hi": m + tc * se,
                "n_rep": n,
            }
        )

    good = ok.groupby(keys).apply(agg, include_groups=False).reset_index() if len(ok) else None
    failed = x[~x["cell_line"].isin(ok["cell_line"])].groupby(keys).size().reset_index()
    out = (
        good.assign(fit_ok=True)
        if good is not None
        else pd.DataFrame(columns=[*keys, "single_a", "single_b", "combo", "delta", "fit_ok"])
    )
    if len(failed):
        fa = x[~x["cell_line"].isin(ok["cell_line"])].groupby(keys)[
            ["single_a", "single_b", "combo"]
        ].mean().reset_index()  # fmt: skip
        fa["delta"] = np.minimum(fa["single_a"], fa["single_b"]) - fa["combo"]
        fa["fit_ok"] = False
        fa["n_rep"] = 0
        out = pd.concat([out, fa], ignore_index=True)
    out["best_single_effect"] = 1 - np.minimum(out["single_a"], out["single_b"])
    out["combo_effect"] = 1 - out["combo"]
    return out


def _lift_panel(
    t: pd.DataFrame,
    group: str | None,
    groups: list[str],
    threshold: float,
    key: str,
    ptitle: str,
    show_legend: bool,
) -> Panel:
    p = Panel(
        key=key,
        title=ptitle,
        x=Scale("linear", label="effect of best single agent (1 − viability)", domain=(0, 1)),
        y=Scale("linear", label="effect of combination", domain=(0, 1)),
        aspect="equal",
        legend=LegendSpec(loc="lower right", title=group) if (group and show_legend) else None,
    )
    p.add(
        PolygonMark(
            x=[0, 1 - threshold, 0], y=[threshold, 1, 1], color="focus_light", alpha=0.35, z=0.4
        )
    )
    p.add(LineMark(x=[0, 1], y=[0, 1], color="ink_soft", width=0.7, style="dashed", z=1))
    p.add(
        TextMark(
            x=0.03, y=0.97, text=f"benefit ≥ {threshold:.0%}", va="top", size=6.5, color="ink_soft"
        )
    )
    p.add(
        TextMark(x=0.97, y=0.04, text="diagonal: no benefit", ha="right", size=6.5,
                 color="ink_soft")
    )  # fmt: skip
    for k, g in enumerate(groups if group else [None]):
        sel = t if g is None else t[t[group] == g]
        col = GROUP_COLOURS[k % len(GROUP_COLOURS)] if g is not None else "focus"
        for ok in (True, False):
            s = sel[sel["fit_ok"] == ok]
            if s.empty:
                continue
            p.add(
                PointMark(
                    x=s["best_single_effect"].to_numpy(float),
                    y=s["combo_effect"].to_numpy(float),
                    color=col,
                    size=16,
                    hollow=not ok,
                    alpha=0.85,
                    categorical=True,
                    legend=(str(g) if ok and show_legend and g is not None else None),
                    z=3,
                )
            )
    return p


def _lift_caption(
    t: pd.DataFrame, group: str | None, groups: list[str], threshold: float, label: str
) -> tuple[str, dict]:
    ok = t[t["fit_ok"]]
    ben = ok["delta"] >= threshold
    stats: dict = {"n": len(ok), "k": int(ben.sum()), "failed": int((~t["fit_ok"]).sum())}
    s = (
        f"{label}{stats['k']} of {stats['n']} cell lines ({stats['k'] / max(stats['n'], 1):.0%}) "
        f"show benefit (combination effect ≥ {threshold:.0%} above the best single agent)"
    )
    if group and groups:
        parts = []
        for g in groups:
            gg = ok[ok[group] == g]
            kb = int((gg["delta"] >= threshold).sum())
            parts.append(f"{g} {kb}/{len(gg)} ({kb / max(len(gg), 1):.0%})")
            stats[f"rate_{g}"] = kb / max(len(gg), 1)
            stats[f"n_{g}"] = len(gg)
            stats[f"k_{g}"] = kb
        s += "; by " + group + ": " + ", ".join(parts)
        if len(groups) == 2:
            a, b = (ok[ok[group] == g] for g in groups)
            tab = [
                [int((a["delta"] >= threshold).sum()), int((a["delta"] < threshold).sum())],
                [int((b["delta"] >= threshold).sum()), int((b["delta"] < threshold).sum())],
            ]
            if min(len(a), len(b)) > 0:
                fp = float(sps.fisher_exact(tab)[1])
                mw = float(sps.mannwhitneyu(a["delta"], b["delta"]).pvalue)
                stats.update(fisher_p=fp, mannwhitney_p=mw)
                s += f" (Fisher p {fmt_p(fp)}; ΔEmax Mann–Whitney p {fmt_p(mw)})"
    return s + ".", stats


def lift(
    table: pd.DataFrame,
    drug_a: str,
    drug_b: str,
    group: str | None = None,
    facet: str | None = None,
    threshold: float = LIFT_THRESHOLD,
    model_note: str = "",
    title: str = "Combination Lift",
    footer: Sequence[str] = (),
) -> PlotResult:
    """Combination Lift across cell lines. ``table``: one row per cell line × replicate (×
    facet, e.g. anchor concentration) with viabilities ``single_a``, ``single_b``, ``combo``
    at Emax and optional ``fit_ok``, group and facet columns."""
    if not 0 < threshold < 1:
        raise TinctureDataError("threshold must be a fraction of viability, e.g. 0.2.")
    t = lift_table(table, group, facet)
    groups = sorted(t[group].dropna().astype(str).unique()) if group else []
    if group:
        t[group] = t[group].astype(str)
    facets = list(dict.fromkeys(t[facet])) if facet else [None]
    panels, caps, stats = [], [], {}
    for k, f in enumerate(facets):
        tf = t if f is None else t[t[facet] == f]
        lab = "" if f is None else f"{f}: "
        panels.append(
            _lift_panel(tf, group, groups, threshold, f"lift_{k}", "" if f is None else str(f),
                        show_legend=(k == len(facets) - 1))
        )  # fmt: skip
        c, st = _lift_caption(tf, group, groups, threshold, lab)
        caps.append(c)
        stats[str(f)] = st
    cap = Caption()
    cap.add("{a} + {b}: {s}", a=drug_a, b=drug_b, s=" ".join(caps))
    nfail = int(t.loc[~t["fit_ok"], "cell_line"].nunique())
    if nfail:
        cap.add("{n} cell line(s) with failed fits are hollow and not counted.", n=nfail)
    if model_note:
        cap.add("{m}", m=model_note)
    meth = Methods().use(
        "combination_lift", threshold=threshold, model_note=model_note or "No synergy model."
    )
    spec = PlotSpec(
        plot="combination_lift",
        panels=panels,
        ncols=len(panels),
        title=title,
        caption=cap.render(),
        methods=meth.render(),
        footer=[
            f"Benefit threshold: ΔEmax ≥ {threshold:.0%} over the best single agent"
            + (f"; {model_note}" if model_note else ""),
            *footer,
        ],
        warnings=list(cap.notes),
        aspect=1.0 / len(panels) + 0.08,
        meta={"table": t, "stats": stats, "threshold": threshold},
    ).validate()
    return PlotResult(spec=spec, classic=build_classic("combination_lift", table=t))


@register_classic("combination_lift", "waterfall of ΔEmax")
def build_lift_classic(table: pd.DataFrame | None = None, **_: object) -> PlotSpec:
    if table is None:
        raise TinctureDataError("Classic needs the lift table.")
    t = table[table["fit_ok"]].sort_values("delta", ascending=False)
    p = Panel(
        x=Scale("linear", label="cell lines (sorted)", ticks=[]),
        y=Scale("linear", label="ΔEmax vs best single agent"),
    )
    p.add(
        RectMark(
            x0=np.arange(len(t)) - 0.4,
            y0=np.zeros(len(t)),
            width=np.full(len(t), 0.8),
            height=t["delta"].to_numpy(float),
            color="mid_grey",
        )
    )
    return PlotSpec(
        plot="combination_lift_classic",
        panels=[p],
        title="Combination benefit (waterfall)",
        caption="ΔEmax per cell line, sorted.",
        aspect=0.5,
    ).validate()
