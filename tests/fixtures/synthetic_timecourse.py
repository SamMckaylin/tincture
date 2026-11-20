"""Synthetic onset tables for unit tests only (never real data)."""

from __future__ import annotations

import pandas as pd


def onset_long(times: list[float], rows: dict[str, list[tuple[float, float]]]) -> pd.DataFrame:
    """``rows``: feature -> [(score, q) per time point]."""
    out = []
    for f, vals in rows.items():
        for t, (s, q) in zip(times, vals, strict=True):
            out.append({"feature": f, "time": t, "score": s, "q": q})
    return pd.DataFrame(out)


def typical() -> pd.DataFrame:
    """Direct pathway down by 3 h; stress up from 12–24 h; one never active."""
    times = [3.0, 6.0, 12.0, 24.0, 48.0]
    ns = (0.5, 0.6)
    return onset_long(
        times,
        {
            "direct": [(-5, 0.001), (-8, 0.001), (-9, 0.001), (-7, 0.001), (-6, 0.001)],
            "transient": [(1, 0.5), (4, 0.01), (5, 0.001), (1, 0.4), ns],
            "stress": [ns, ns, (1, 0.3), (4, 0.001), (6, 0.001)],
            "quiet": [ns, ns, ns, ns, ns],
        },
    )
