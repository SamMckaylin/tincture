"""Onset, peak and resolution of a response over sparse time points (HANDOVER §5.12).

Input: per feature (e.g. pathway) and time point, an effect ``score`` and a ``q``-value (or
``p``, adjusted here with Benjamini–Hochberg across all feature × time tests).
A time point is **active** when q < ``alpha`` and |score| ≥ ``threshold``, in the direction
of the feature's strongest active effect.

* **onset** is *interval-censored*: between the last inactive time point before the first
  active one and that first active one, i.e. ``(t_prev, t_first]``; ``(0, t_1]`` if already
  active at the first time point. Never a single precise time.
* **peak**: the active time point with the largest |score| (a measured time point).
* **resolution**: the response returns within threshold between the last active time point
  after the peak and the next inactive one, ``(t_last, t_next]``; if still active at the
  last time point it is right-censored (``sustained``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from tincture.errors import TinctureDataError


@dataclass
class Onset:
    feature: str
    direction: int  # +1 up, −1 down, 0 never active
    onset: tuple[float, float] | None  # (lo, hi] interval
    peak_time: float | None
    peak_score: float | None
    last_active: float | None
    resolution: tuple[float, float] | None  # (lo, hi]; None if sustained or never active
    sustained: bool
    active: list[float]  # active time points


def onset_table(
    d: pd.DataFrame,
    feature: str = "feature",
    threshold: float = 0.0,
    alpha: float = 0.05,
) -> tuple[list[Onset], pd.DataFrame]:
    """Onsets for every feature. ``d`` has columns ``feature``, ``time``, ``score`` and ``q``
    (or ``p``)."""
    need = {feature, "time", "score"}
    if not need <= set(d.columns) or not ({"q", "p"} & set(d.columns)):
        raise TinctureDataError(
            f"Onset table needs columns {sorted(need)} and q (or p).",
            fix="One row per feature × time point with an effect score and its q-value.",
        )
    x = d.copy()
    if "q" not in x:
        ok = np.isfinite(x["p"])
        x["q"] = np.nan
        x.loc[ok, "q"] = multipletests(x.loc[ok, "p"], method="fdr_bh")[1]
    times = np.sort(x["time"].unique().astype(float))
    if len(times) < 2:
        raise TinctureDataError(
            "An onset needs at least 2 time points (3 or more recommended).",
            fix="Add time points, or show a single-time-point plot instead.",
        )
    if (times < 0).any():
        raise TinctureDataError("Time points must be ≥ 0 (time since treatment).")
    out: list[Onset] = []
    for f, g in x.groupby(feature, sort=False):
        g = g.set_index("time").reindex(times)
        s = g["score"].to_numpy(float)
        sig = (g["q"].to_numpy(float) < alpha) & (np.abs(s) >= threshold) & np.isfinite(s)
        if not sig.any():
            out.append(Onset(str(f), 0, None, None, None, None, None, False, []))
            continue
        k = int(np.nanargmax(np.where(sig, np.abs(s), -np.inf)))
        direction = int(np.sign(s[k]))
        act = sig & (np.sign(s) == direction)
        first = int(np.argmax(act))
        lo = 0.0 if first == 0 else float(times[first - 1])
        later_inactive = [j for j in range(k + 1, len(times)) if not act[j]]
        if later_inactive:
            j = later_inactive[0]
            res: tuple[float, float] | None = (float(times[j - 1]), float(times[j]))
            last = float(times[j - 1])
            sustained = False
        else:
            res = None
            last = float(times[-1])
            sustained = True
        out.append(
            Onset(
                feature=str(f),
                direction=direction,
                onset=(lo, float(times[first])),
                peak_time=float(times[k]),
                peak_score=float(s[k]),
                last_active=last,
                resolution=res,
                sustained=sustained,
                active=[float(t) for t in times[act]],
            )
        )
    return out, x
