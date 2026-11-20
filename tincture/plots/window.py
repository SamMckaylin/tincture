"""Window: is there a dose that works without the toxicity? (HANDOVER §5.7)

On-target and toxicity dose-response curves on one log-concentration axis, raw points
included. The window is the concentration range between EC90 of the on-target curve and
EC10 of the most sensitive toxicity curve (definitions in :mod:`tincture.stats.hill`); it
exists if EC90(on) < EC10(tox). Optional free Cmax is drawn as a range band, never a line.

Guards (each tested):
* on-target fit failed (no convergence, R² < 0.8, or maximal effect < 50%) → no band, and the
  caption says why;
* a toxicity curve with no decline beyond the EC10 threshold anywhere in the tested range is
  *censored*, not failed: EC10 is reported as "> top dose" and the band is open-ended;
* a toxicity curve that declines but fits badly → no band;
* assay context (cell line, readout, duration) is printed on each curve.
Classic fallback: overlaid curves with IC50s.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from tincture.errors import TinctureDataError
from tincture.model.dose_response import DoseResponse, FourPL
from tincture.report.captions import Caption, fmt_conc
from tincture.report.classic import build_classic, register_classic
from tincture.report.methods import Methods, register_block
from tincture.result import PlotResult
from tincture.spec import (
    BandMark,
    IntervalMark,
    LegendSpec,
    LineMark,
    NoticeMark,
    Panel,
    PlotSpec,
    PointMark,
    RefLineMark,
    Scale,
    TextMark,
)
from tincture.stats.hill import (
    HillFit,
    bootstrap_4pl,
    ec_drop,
    ec_relative,
    fit_4pl,
    median_of_units,
)

MIN_R2 = 0.8
MIN_EFFECT = 0.5

register_block(
    "window",
    "Curves are 4-parameter log-logistic fits (bounded least squares on log10 dose, soft-L1 "
    "loss). EC{on_pct:.0f} of the on-target curve is the concentration giving "
    "{on_pct:.0f}% of its fitted span; EC{tox_pct:.0f} of a toxicity curve is the "
    "concentration at which viability falls {tox_pct:.0f}% below its fitted top. 95% intervals "
    "from {n_boot} bootstrap refits (replicates resampled within dose). A window exists if "
    "EC{on_pct:.0f}(on) < EC{tox_pct:.0f}(tox), using the most sensitive toxicity curve. Fits "
    "are rejected when they do not converge, R² < {min_r2}, or (on-target) the maximal effect "
    "is < {min_eff:.0%}. A toxicity curve that never falls {tox_pct:.0f}% in the tested range is "
    "reported as EC{tox_pct:.0f} > top dose.",
)


@dataclass
class CurveCall:
    """Fit, effective concentration and verdict for one curve."""

    dr: DoseResponse
    hill: HillFit
    ec: float  # EC90 (on) or EC10 (tox); inf = censored above the tested range
    ec_ci: tuple[float, float]
    problems: list[str] = field(default_factory=list)
    censored: bool = False

    @property
    def usable(self) -> bool:
        return not self.problems


def _fit_problems(h: HillFit, need_effect: bool) -> list[str]:
    out = []
    if not h.converged or not np.isfinite(h.fit.hill):
        out.append("fit did not converge")
        return out
    if not np.isfinite(h.r2) or h.r2 < MIN_R2:
        out.append(f"R² {h.r2:.2f} < {MIN_R2}")
    if need_effect and np.isfinite(h.fit.hill) and h.max_effect < MIN_EFFECT:
        out.append(f"maximal effect {h.max_effect:.0%} < {MIN_EFFECT:.0%}")
    return out


def _no_decline(dr: DoseResponse, drop: float, units: np.ndarray | None = None) -> bool:
    """True if the response (per-dose mean, or median of units) never falls ``drop`` below
    its low-dose level."""
    if units is not None:
        u, means = median_of_units(dr.dose, dr.response, units)
    else:
        d, y = dr.dose[dr.dose > 0], dr.response[dr.dose > 0]
        u = np.unique(d)
        means = np.array([np.nanmean(y[d == v]) for v in u])
    if u.size < 2:
        return False
    base = float(np.mean(means[: max(1, u.size // 4)]))
    return bool(np.all(means >= base - drop))


def call_on(dr: DoseResponse, frac: float, n_boot: int, seed: int = 0) -> CurveCall:
    h = fit_4pl(dr.dose, dr.response)
    probs = _fit_problems(h, need_effect=True)
    ec = ec_relative(h.fit, frac) if not probs else float("nan")
    ci = (float("nan"), float("nan"))
    if not probs:
        b = bootstrap_4pl(dr.dose, dr.response, n_boot=n_boot, seed=seed)
        ci = b.interval(lambda f: ec_relative(f, frac))
    return CurveCall(dr, h, ec, ci, probs)


def _units(dr: DoseResponse) -> np.ndarray | None:
    """Unit labels when ``dr`` summarises many units (``context["summary"] == "median"``);
    the replicate column then holds the unit (e.g. cell line)."""
    if (dr.context or {}).get("summary") == "median":
        if dr.replicate is None:
            raise TinctureDataError("A median-of-units curve needs unit ids in `replicate`.")
        return dr.replicate
    return None


def call_tox(dr: DoseResponse, drop: float, n_boot: int, seed: int = 0) -> CurveCall:
    units = _units(dr)
    h = (
        fit_4pl(*median_of_units(dr.dose, dr.response, units))
        if units is not None
        else fit_4pl(dr.dose, dr.response)
    )
    top = dr.dose_range[1]
    if _no_decline(dr, drop, units):
        return CurveCall(dr, h, float("inf"), (top, float("inf")), [], censored=True)
    probs = _fit_problems(h, need_effect=False)
    ec = ec_drop(h.fit, drop) if not probs else float("nan")
    ci = (float("nan"), float("nan"))
    if not probs:
        b = bootstrap_4pl(dr.dose, dr.response, n_boot=n_boot, seed=seed, cluster=units)
        ci = b.interval(lambda f: ec_drop(f, drop))
    return CurveCall(dr, h, ec, ci, probs, censored=not np.isfinite(ec))


def _context_label(dr: DoseResponse, fallback: str) -> str:
    c = dr.context or {}
    parts = [
        str(c.get(k))
        for k in ("cell_line", "ccle_name", "readout", "screen_id", "duration")
        if c.get(k) not in (None, "", "nan")
    ]
    return " · ".join(dict.fromkeys(parts)) or (dr.assay or fallback)


def _curve(f: FourPL, lo: float, hi: float) -> tuple[np.ndarray, np.ndarray]:
    x = np.logspace(np.log10(lo), np.log10(hi), 120)
    return x, f(x)


@dataclass
class WindowResult:
    on: CurveCall
    tox: list[CurveCall]
    limiting: CurveCall | None  # most sensitive usable toxicity curve
    exists: bool | None  # None = cannot be judged (a fit failed)
    reason: str = ""
    overlap: bool = False  # CIs of EC(on) and EC(tox) overlap


def compute(
    on: DoseResponse,
    tox: Sequence[DoseResponse],
    on_frac: float = 0.9,
    tox_drop: float = 0.1,
    n_boot: int = 200,
) -> WindowResult:
    c_on = call_on(on, on_frac, n_boot)
    c_tox = [call_tox(t, tox_drop, n_boot, seed=i + 1) for i, t in enumerate(tox)]
    if not c_on.usable:
        return WindowResult(
            c_on, c_tox, None, None, "on-target fit rejected: " + "; ".join(c_on.problems)
        )
    bad = [t for t in c_tox if not t.usable]
    if bad:
        why = "; ".join(f"{_context_label(t.dr, 'tox')}: {', '.join(t.problems)}" for t in bad)
        return WindowResult(c_on, c_tox, None, None, "toxicity fit rejected: " + why)
    lim = min(c_tox, key=lambda t: t.ec)
    exists = bool(c_on.ec < lim.ec)
    overlap = bool(np.isfinite(lim.ec) and c_on.ec_ci[1] >= lim.ec_ci[0])
    return WindowResult(c_on, c_tox, lim, exists, overlap=overlap)


def build_spec(
    res: WindowResult,
    on_frac: float,
    tox_drop: float,
    cmax: tuple[float, float] | None,
    n_boot: int,
    title: str,
    footer: Sequence[str],
    compound: str,
) -> PlotSpec:
    on_pct, tox_pct = on_frac * 100, tox_drop * 100
    allc = [res.on, *res.tox]
    lo = min(c.dr.dose_range[0] for c in allc) / 3
    hi = max(c.dr.dose_range[1] for c in allc) * 3
    p = Panel(
        key="window",
        x=Scale("log", label="concentration (M)", domain=(lo, hi)),
        y=Scale("linear", label="viability (fraction of vehicle)", domain=(-0.05, 1.3)),
        legend=None,
    )
    cap = Caption()
    if res.exists and res.limiting is not None:
        right = res.limiting.ec if np.isfinite(res.limiting.ec) else hi
        p.add(
            BandMark(
                lo=res.on.ec,
                hi=right,
                span=True,
                orientation="v",
                color="window",
                alpha=0.55,
                inner_label="window",
                z=0.2,
            )
        )
    if cmax is not None:
        p.add(
            BandMark(
                lo=cmax[0],
                hi=cmax[1],
                span=True,
                orientation="v",
                color="focus_light",
                alpha=0.35,
                hatch="//",
                inner_label="free Cmax",
                z=0.3,
            )
        )
    for c, col, role in [
        (res.on, "focus", "on-target"),
        *[(t, "reference", "toxicity") for t in res.tox],
    ]:
        dr = c.dr
        m = dr.dose > 0
        lab = f"{role}: {_context_label(dr, role)}"
        units = _units(dr)
        if units is not None:  # median of many units: show medians with the interquartile range
            import pandas as pd

            df = pd.DataFrame({"d": dr.dose, "y": dr.response, "u": units})
            df = df[df["d"] > 0].groupby(["u", "d"])["y"].mean().groupby(level="d")
            q = df.quantile([0.25, 0.5, 0.75]).unstack()
            xs_ = q.index.to_numpy(float)
            p.add(
                IntervalMark(
                    pos=xs_,
                    lo=q[0.25].to_numpy(),
                    hi=q[0.75].to_numpy(),
                    orientation="v",
                    color=col,
                    width=1.0,
                )
            )
            p.add(PointMark(x=xs_, y=q[0.5].to_numpy(), color=col, size=16, categorical=True))
        else:
            p.add(
                PointMark(
                    x=dr.dose[m], y=dr.response[m], color=col, size=9, alpha=0.6, categorical=True
                )
            )
        if np.isfinite(c.hill.fit.hill):
            xs, ys = _curve(c.hill.fit, *dr.dose_range)
            p.add(
                LineMark(
                    x=xs,
                    y=ys,
                    color=col,
                    width=1.6 if c.usable else 0.8,
                    style="solid" if c.usable else "dashed",
                    categorical=True,
                    z=2,
                )
            )
        else:
            p.add(
                TextMark(
                    x=dr.dose_range[1],
                    y=float(np.nanmean(dr.response[m])),
                    text=f"  {lab}",
                    color=col,
                    ha="left",
                    categorical=True,
                )
            )
    # EC markers with CIs
    y_on = res.on.hill.fit.top - on_frac * (res.on.hill.fit.top - res.on.hill.fit.bottom)
    if res.on.usable and np.isfinite(res.on.ec):
        p.add(
            PointMark(
                x=[res.on.ec], y=[y_on], color="focus", size=34, marker="D", z=3, categorical=True
            )
        )
        if np.all(np.isfinite(res.on.ec_ci)):
            p.add(
                IntervalMark(
                    pos=[y_on],
                    lo=[res.on.ec_ci[0]],
                    hi=[res.on.ec_ci[1]],
                    orientation="h",
                    color="focus",
                    width=1.4,
                )
            )
        p.add(
            TextMark(
                x=res.on.ec,
                y=y_on - 0.07,
                text=f"EC{on_pct:.0f} {fmt_conc(res.on.ec)}",
                color="focus",
                ha="center",
                va="top",
                categorical=True,
            )
        )
    lim = res.limiting
    if lim is not None:
        if np.isfinite(lim.ec):
            y_t = lim.hill.fit.top - tox_drop
            p.add(
                PointMark(
                    x=[lim.ec],
                    y=[y_t],
                    color="reference",
                    size=34,
                    marker="D",
                    z=3,
                    categorical=True,
                )
            )
            if np.all(np.isfinite(lim.ec_ci)):
                p.add(
                    IntervalMark(
                        pos=[y_t],
                        lo=[lim.ec_ci[0]],
                        hi=[lim.ec_ci[1]],
                        orientation="h",
                        color="reference",
                        width=1.4,
                    )
                )
            p.add(
                TextMark(
                    x=lim.ec,
                    y=y_t + 0.06,
                    text=f"EC{tox_pct:.0f} {fmt_conc(lim.ec)}",
                    color="reference",
                    ha="center",
                    va="bottom",
                    categorical=True,
                )
            )
        else:
            p.add(
                TextMark(
                    x=lim.dr.dose_range[1],
                    y=1.12,
                    text=f"tox EC{tox_pct:.0f} > {fmt_conc(lim.dr.dose_range[1])} →",
                    color="reference",
                    ha="right",
                    va="bottom",
                    categorical=True,
                )
            )
    p.add(RefLineMark(axis="y", value=1.0, color="ink_soft", style="dotted", width=0.5))
    # assay context of every curve, as a key in the lower right (long labels never clip)
    key = [(res.on, "focus", "on-target"), *[(t, "reference", "toxicity") for t in res.tox]]
    for k, (c, col, role) in enumerate(reversed(key)):
        p.add(
            TextMark(
                x=0.985,
                y=0.035 + 0.055 * k,
                coords="axes",
                text=f"{role}: {_context_label(c.dr, role)}",
                color=col,
                ha="right",
                va="bottom",
                categorical=True,
            )
        )
    # caption
    if res.exists is None:
        p.add(NoticeMark(text="no window drawn: " + res.reason, loc="upper right", level="warning"))
        cap.weak(f"No window is drawn because the {res.reason}.")
    elif res.exists:
        assert lim is not None
        upper = (
            f"EC{tox_pct:.0f}(tox) {fmt_conc(lim.ec)}"
            if np.isfinite(lim.ec)
            else f"above the top tested dose ({fmt_conc(lim.dr.dose_range[1])})"
        )
        cap.add(
            "{c} reaches {on:.0f}% of its effect at {ec} (95% CI {lo}–{hi}) before the toxicity "
            "readout loses {tox:.0f}% viability ({up}): a window of {fold}.",
            c=compound,
            on=on_pct,
            ec=fmt_conc(res.on.ec),
            lo=fmt_conc(res.on.ec_ci[0]),
            hi=fmt_conc(res.on.ec_ci[1]),
            tox=tox_pct,
            up=upper,
            fold=(
                f"{lim.ec / res.on.ec:.0f}-fold"
                if np.isfinite(lim.ec)
                else f"at least {lim.dr.dose_range[1] / res.on.ec:.0f}-fold"
            ),
        )
        if res.overlap:
            cap.weak(
                "the confidence intervals of the two effective concentrations overlap; "
                "the window is not resolved from noise."
            )
    else:
        assert lim is not None
        cap.add(
            "No window: {c} reaches EC{on:.0f} at {ec}, above the toxicity EC{tox:.0f} ({t}).",
            c=compound,
            on=on_pct,
            ec=fmt_conc(res.on.ec),
            tox=tox_pct,
            t=fmt_conc(lim.ec),
        )
    if cmax is None:
        cap.add(
            "No exposure (free Cmax) was supplied, so the window is not placed against "
            "clinical concentrations."
        )
    meth = Methods().use(
        "window", on_pct=on_pct, tox_pct=tox_pct, n_boot=n_boot, min_r2=MIN_R2, min_eff=MIN_EFFECT
    )
    return PlotSpec(
        plot="window",
        panels=[p],
        title=title,
        caption=cap.render(),
        methods=meth.render(),
        footer=[
            f"Diamonds = EC{on_pct:.0f}(on) / EC{tox_pct:.0f}(tox) with 95% bootstrap "
            "intervals; points = replicate wells"
            + (
                "; median-of-lines curves show per-dose medians with interquartile ranges"
                if any(_units(t.dr) is not None for t in res.tox)
                else ""
            ),
            *footer,
        ],
        warnings=list(cap.notes),
        aspect=0.62,
        meta={"result": res},
    ).validate()


@register_classic("window", "overlaid dose-response curves")
def build_classic_spec(res: WindowResult | None = None, **_: object) -> PlotSpec:
    if res is None:
        raise TinctureDataError("Classic Window needs the fitted curves.")
    allc = [res.on, *res.tox]
    lo = min(c.dr.dose_range[0] for c in allc) / 3
    hi = max(c.dr.dose_range[1] for c in allc) * 3
    p = Panel(
        x=Scale("log", label="concentration (M)", domain=(lo, hi)),
        y=Scale(label="viability"),
        legend=LegendSpec(loc="lower left"),
    )
    for c, col, role in [
        (res.on, "focus", "on-target"),
        *[(t, "reference", "toxicity") for t in res.tox],
    ]:
        units = _units(c.dr)
        if units is not None:
            xs_, ys_ = median_of_units(c.dr.dose, c.dr.response, units)
        else:
            m = c.dr.dose > 0
            xs_, ys_ = c.dr.dose[m], c.dr.response[m]
        p.add(PointMark(x=xs_, y=ys_, color=col, size=8, categorical=True))
        if np.isfinite(c.hill.fit.hill):
            xs, ys = _curve(c.hill.fit, *c.dr.dose_range)
            ic50 = c.hill.fit.ec50
            p.add(
                LineMark(
                    x=xs,
                    y=ys,
                    color=col,
                    categorical=True,
                    legend=f"{role} (IC50 {fmt_conc(ic50)})",
                )
            )
    return PlotSpec(
        plot="window_classic",
        panels=[p],
        title="Dose-response curves",
        caption="Overlaid dose-response curves with IC50s.",
        aspect=0.62,
    ).validate()


def plot(
    on: DoseResponse,
    tox: DoseResponse | Sequence[DoseResponse],
    cmax: tuple[float, float] | None = None,
    on_frac: float = 0.9,
    tox_drop: float = 0.1,
    n_boot: int = 200,
    compound: str | None = None,
    title: str = "Window",
    footer: Sequence[str] = (),
) -> PlotResult:
    """Window plot. ``cmax``: free Cmax range in molar (low, high), drawn as a band; omit if
    no exposure data are available (never pass a guess)."""
    toxs = [tox] if isinstance(tox, DoseResponse) else list(tox)
    if not toxs:
        raise TinctureDataError("Window needs at least one toxicity / counter-screen curve.")
    if cmax is not None and not (0 < cmax[0] <= cmax[1] < 1):
        raise TinctureDataError("cmax must be a (low, high) range in molar.")
    res = compute(on, toxs, on_frac, tox_drop, n_boot)
    spec = build_spec(
        res, on_frac, tox_drop, cmax, n_boot, title, footer, compound or on.compound_id
    )
    return PlotResult(spec=spec, classic=build_classic("window", res=res))
