"""Onset Timeline: what happens first, and what follows? (HANDOVER §5.12)

Gantt-style rows (one per pathway, ordered by onset) on a time axis. The **onset** is an
interval bracket between the last inactive and the first active time point (interval-
censored, :mod:`tincture.stats.onset`), the **peak** a dot at the measured time point with
the largest effect, and the **sustained** effect a bar from the first to the last active time
point (an arrow-free open end when still active at the last time point). Colour = direction
(terracotta up, blue down). Pathways that never become active are listed in grey at the bottom.

Guards (each tested): onset is never drawn as a single point; ≥ 3 time points (2 allowed with
a warning in the caption); ``style="horizon"`` only for dense series (≥ 6 time points).
Classic fallback: score-by-time line chart.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

import numpy as np
import pandas as pd

from tincture.errors import TinctureDataError
from tincture.report.captions import Caption
from tincture.report.classic import build_classic, register_classic
from tincture.report.methods import Methods, register_block
from tincture.result import PlotResult
from tincture.spec import (
    IntervalMark,
    LineMark,
    Panel,
    PlotSpec,
    PointMark,
    PolygonMark,
    RectMark,
    Scale,
)
from tincture.stats.onset import Onset, onset_table

HORIZON_MIN_POINTS = 6

register_block(
    "onset",
    "A pathway is active at a time point when its effect has BH q < {alpha} and |score| ≥ "
    "{threshold:g} ({score_label}). Onset is interval-censored: it lies between the last "
    "inactive and the first active time point (sampled at {times}). The peak is the active "
    "time point with the largest |score|; the effect resolves when it next becomes inactive.",
)


def _fmt_t(t: float, unit: str) -> str:
    return f"{t:g} {unit}"


def _interval_text(iv: tuple[float, float], unit: str) -> str:
    lo, hi = iv
    return f"({lo:g}, {hi:g}] {unit}"


def plot(
    table: pd.DataFrame,
    feature: str = "feature",
    threshold: float = 0.0,
    alpha: float = 0.05,
    time_unit: str = "h",
    score_label: str = "effect score",
    style: Literal["gantt", "horizon"] = "gantt",
    title: str = "Onset Timeline",
    footer: Sequence[str] = (),
) -> PlotResult:
    """Onset Timeline from a long table (``feature``, ``time``, ``score``, ``q`` or ``p``)."""
    onsets, x = onset_table(table, feature, threshold, alpha)
    times = np.sort(x["time"].unique().astype(float))
    if style == "horizon":
        if len(times) < HORIZON_MIN_POINTS:
            raise TinctureDataError(
                f"style='horizon' needs a dense time series (≥ {HORIZON_MIN_POINTS} time "
                f"points); this one has {len(times)}.",
                fix="Use the default Gantt style for sparse time points.",
            )
        spec = _horizon_spec(x, feature, onsets, times, time_unit, score_label, title, footer)
        return PlotResult(spec=spec, classic=build_classic("onset", table=x, feature=feature))
    act = [o for o in onsets if o.onset is not None]
    act.sort(key=lambda o: (o.onset[1], o.onset[0], o.peak_time))  # type: ignore[index]
    inact = [o for o in onsets if o.onset is None]
    rows = act
    n = max(len(rows), 1)
    use_log = times[0] > 0 and times[-1] / times[0] >= 8
    lo_dom = times[0] / 1.6 if use_log else min(0.0, times[0]) - 0.04 * times[-1]
    hi_dom = times[-1] * 1.25 if use_log else times[-1] * 1.04
    left = lo_dom * 1.08 if use_log else lo_dom + 0.01 * times[-1]
    p = Panel(
        key="timeline",
        x=Scale(
            "log" if use_log else "linear",
            label=f"time after treatment ({time_unit})",
            domain=(lo_dom, hi_dom),
            ticks=[float(t) for t in times],
            tick_labels=[f"{t:g}" for t in times],
        ),
        y=Scale(
            "linear",
            label="",
            domain=(n - 0.4, -0.6),
            ticks=list(range(n)),
            tick_labels=[o.feature for o in rows],
        ),
    )
    for t in times:
        p.add(LineMark(x=[t, t], y=[-0.6, n - 0.4], color="grid", width=0.5, z=0.2))
    for i, o in enumerate(rows):
        if o.onset is None:
            continue
        col = "up" if o.direction > 0 else "down"
        light = "up_light" if o.direction > 0 else "down_light"
        lo, hi = o.onset
        a = lo if lo > 0 or not use_log else left
        end = o.last_active if o.last_active is not None else hi
        bar_end = hi_dom / 1.05 if o.sustained else end
        if bar_end > hi:
            p.add(
                RectMark(
                    x0=[hi], y0=[i - 0.16], width=[bar_end - hi], height=[0.32],
                    color=light, z=1,
                )
            )  # fmt: skip
        p.add(
            IntervalMark(
                pos=[i], lo=[a], hi=[hi], orientation="h", color=col, width=1.4, cap=0.32, z=2
            )
        )
        p.add(
            PointMark(
                x=[o.peak_time], y=[i], color=col, size=22, categorical=True, z=3,
            )
        )  # fmt: skip
    # caption
    cap = Caption()
    if len(times) == 2:
        cap.weak("only 2 time points: onsets are bounded very coarsely.")
    if act:
        first = act[0]
        tied = [o for o in act if o.onset == first.onset]
        cap.add(
            "First to respond (onset in {iv}): {fs}.",
            iv=_interval_text(first.onset, time_unit),  # type: ignore[arg-type]
            fs=", ".join(
                f"{o.feature} ({'up' if o.direction > 0 else 'down'}, peak at "
                f"{_fmt_t(o.peak_time, time_unit)})"  # type: ignore[arg-type]
                for o in tied
            ),
        )
        later = [o for o in act[1:] if o.onset[0] >= first.onset[1]]  # type: ignore[index]
        if later:
            cap.add(
                "{k} pathway(s) start strictly later, e.g. {f} (onset in {iv}).",
                k=len(later),
                f=later[0].feature,
                iv=_interval_text(later[0].onset, time_unit),  # type: ignore[arg-type]
            )
        sus = [o.feature for o in act if o.sustained]
        if sus:
            cap.add(
                "{k} remain active at the last time point ({t}).",
                k=len(sus),
                t=_fmt_t(times[-1], time_unit),
            )
    else:
        cap.weak("no pathway is active at any time point.")
    if inact:
        cap.add(
            "{k} pathway(s) never become active: {fs}.",
            k=len(inact),
            fs=", ".join(o.feature for o in inact),
        )
    cap.add(
        "Onsets are intervals between sampled time points ({ts} {u}); the true onset can lie "
        "anywhere inside the bracket.",
        ts=", ".join(f"{t:g}" for t in times),
        u=time_unit,
    )
    meth = Methods().use(
        "onset",
        alpha=alpha,
        threshold=threshold,
        score_label=score_label,
        times=", ".join(_fmt_t(t, time_unit) for t in times),
    )
    spec = PlotSpec(
        plot="onset",
        panels=[p],
        title=title,
        caption=cap.render(),
        methods=meth.render(),
        footer=[
            f"Active: q < {alpha}, |score| ≥ {threshold:g}. Bracket = onset between sampled "
            "time points; dot = peak; bar = still active (open end: active at the last point)",
            *(["Never active: " + ", ".join(o.feature for o in inact)] if inact else []),
            *footer,
        ],
        warnings=list(cap.notes),
        aspect=min(1.1, 0.16 + 0.034 * n),
        meta={"onsets": onsets, "table": x, "order": [o.feature for o in rows + inact]},
    ).validate()
    return PlotResult(spec=spec, classic=build_classic("onset", table=x, feature=feature))


def _horizon_spec(
    x: pd.DataFrame,
    feature: str,
    onsets: list[Onset],
    times: np.ndarray,
    unit: str,
    score_label: str,
    title: str,
    footer: Sequence[str],
) -> PlotSpec:
    """Two-band horizon chart: per row, |score| folded into two bands, coloured by sign."""
    feats = [o.feature for o in onsets]
    vmax = float(np.nanmax(np.abs(x["score"]))) or 1.0
    band = vmax / 2
    p = Panel(
        key="horizon",
        x=Scale("linear", label=f"time after treatment ({unit})",
                domain=(float(times[0]), float(times[-1]))),
        y=Scale("linear", label="", domain=(len(feats), 0), ticks=[i + 0.5 for i in
                range(len(feats))], tick_labels=feats),
    )  # fmt: skip
    for i, f in enumerate(feats):
        g = x[x[feature] == f].set_index("time").reindex(times)["score"].to_numpy(float)
        g = np.nan_to_num(g)
        for sign, cols in ((1, ("up_light", "up")), (-1, ("down_light", "down"))):
            v = np.clip(sign * g, 0, None)
            for b, col in enumerate(cols):
                h = np.clip(v - b * band, 0, band) / band * 0.9
                if not (h > 0).any():
                    continue
                p.add(
                    PolygonMark(
                        x=np.r_[times, times[::-1]],
                        y=np.r_[i + 1 - h, np.full(len(times), i + 1.0)],
                        color=col, alpha=1.0, z=1 + b,
                    )
                )  # fmt: skip
    cap = Caption()
    cap.add(
        "Horizon chart of {s}: each row folds |score| into two bands of {b:.2g} (darker = "
        "stronger); terracotta up, blue down.",
        s=score_label,
        b=band,
    )
    return PlotSpec(
        plot="onset_horizon",
        panels=[p],
        title=title,
        caption=cap.render(),
        footer=list(footer),
        aspect=min(1.2, 0.18 + 0.045 * len(feats)),
        meta={"onsets": onsets, "band": band},
    ).validate()


@register_classic("onset", "score-by-time line chart")
def build_classic_spec(
    table: pd.DataFrame | None = None, feature: str = "feature", **_: object
) -> PlotSpec:
    if table is None:
        raise TinctureDataError("Classic needs the onset table.")
    p = Panel(x=Scale("linear", label="time"), y=Scale("linear", label="score"))
    for _f, g in table.groupby(feature):
        g = g.sort_values("time")
        p.add(LineMark(x=g["time"].to_numpy(float), y=g["score"].to_numpy(float),
                       color="mid_grey", width=0.8))  # fmt: skip
    return PlotSpec(
        plot="onset_classic", panels=[p], title="Scores over time", caption="Scores over time.",
        aspect=0.6,
    ).validate()  # fmt: skip
