"""4-parameter log-logistic (Hill) fits with bootstrap uncertainty.

Model (:class:`tincture.model.dose_response.FourPL`):
``y = bottom + (top − bottom) / (1 + (x / EC50) ** hill)``, x in molar, fitted on
log10(x) with bounded least squares (scipy ``least_squares``, soft-L1 loss so a single bad
well does not drag the curve).

Bounds: ``top, bottom ∈ [−0.5, 1.5]`` (viability as fraction of vehicle), ``hill ∈ [0.1, 10]``,
``log10 EC50`` within two decades of the tested range.

**Effective concentrations** (documented choices, used by the Window plot):

* ``ec_relative(fit, f)``: concentration giving fraction ``f`` of the curve's *own* span,
  e.g. EC90 = 90% of the way from top to bottom. Meaningful only if the span is large, so
  the Window plot requires a maximal effect ≥ 50% (handover guard).
* ``ec_drop(fit, d)``: concentration at which the response falls by ``d`` (absolute) below
  the curve's top, e.g. EC10 of a toxicity readout = 10% loss of viability. Curves that never
  fall that far inside the tested range return ``inf`` (censored: "> top dose").

**Bootstrap:** replicates are resampled with replacement within each dose and the curve is
refitted; percentile intervals.

**Median-of-units curves** (``cluster=``): when a curve summarises many units (e.g. the
typical response of hundreds of cell lines), raw wells from different units are not
replicates. The curve is then fitted to the per-dose median of unit means, and the bootstrap
resamples whole units, so the interval reflects between-unit variation.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from tincture.model.dose_response import FourPL

BOUNDS_Y = (-0.5, 1.5)
BOUNDS_HILL = (0.1, 10.0)
DECADES = 2.0


@dataclass
class HillFit:
    fit: FourPL
    r2: float
    converged: bool
    n_points: int
    dose_range: tuple[float, float]

    @property
    def max_effect(self) -> float:
        """Top − bottom (fraction of vehicle lost at infinite dose)."""
        return float(self.fit.top - self.fit.bottom)


def _model(p: np.ndarray, lx: np.ndarray) -> np.ndarray:
    bottom, top, le50, hill = p
    with np.errstate(over="ignore", invalid="ignore"):
        return bottom + (top - bottom) / (1.0 + 10.0 ** (hill * (lx - le50)))


def _jac(p: np.ndarray, lx: np.ndarray) -> np.ndarray:
    """Analytic Jacobian of :func:`_model` (columns: bottom, top, log10 EC50, hill)."""
    bottom, top, le50, hill = p
    with np.errstate(over="ignore", invalid="ignore"):
        u = 10.0 ** (hill * (lx - le50))
        d = 1.0 + u
        g = (top - bottom) * u * np.log(10.0) / d**2
        j = np.column_stack([u / d, 1.0 / d, hill * g, -(lx - le50) * g])
    return np.nan_to_num(j)


def fit_4pl(dose: np.ndarray, response: np.ndarray, start: FourPL | None = None) -> HillFit:
    """Fit a 4PL to treated points (dose > 0, molar). ``start`` warm-starts a single
    optimisation (used by the bootstrap); otherwise six starting points are tried."""
    d = np.asarray(dose, float)
    y = np.asarray(response, float)
    ok = (d > 0) & np.isfinite(y)
    d, y = d[ok], y[ok]
    lx = np.log10(d)
    lo_x, hi_x = float(lx.min()), float(lx.max())
    if np.unique(lx).size < 4:
        return HillFit(FourPL(np.nan, np.nan, np.nan, np.nan), np.nan, False, len(y),
                       (10**lo_x, 10**hi_x))  # fmt: skip
    top0 = float(np.clip(np.median(y[lx <= np.quantile(lx, 0.25)]), *BOUNDS_Y))
    bot0 = float(np.clip(np.median(y[lx >= np.quantile(lx, 0.75)]), *BOUNDS_Y))
    lb = [BOUNDS_Y[0], BOUNDS_Y[0], lo_x - DECADES, BOUNDS_HILL[0]]
    ub = [BOUNDS_Y[1], BOUNDS_Y[1], hi_x + DECADES, BOUNDS_HILL[1]]
    best = None
    starts = (
        [(start.bottom, start.top, start.log10_ec50, start.hill)]
        if start is not None
        else [(bot0, top0, le0, h0) for le0 in (np.median(lx), lo_x, hi_x) for h0 in (1.0, 3.0)]
    )
    for st in starts:
        p0 = np.clip(st, np.array(lb) + 1e-9, np.array(ub) - 1e-9)
        r = least_squares(
            lambda p: _model(p, lx) - y,
            p0,
            jac=lambda p: _jac(p, lx),
            bounds=(lb, ub),
            loss="soft_l1",
            f_scale=0.1,
        )
        if best is None or r.cost < best.cost:
            best = r
    assert best is not None
    pred = _model(best.x, lx)
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    b, t, le, h = (float(v) for v in best.x)
    return HillFit(FourPL(bottom=b, top=t, log10_ec50=le, hill=h), r2, bool(best.success),
                   len(y), (10**lo_x, 10**hi_x))  # fmt: skip


def ec_relative(fit: FourPL, frac: float) -> float:
    """Concentration giving ``frac`` of the curve's span (EC90: frac = 0.9)."""
    if not np.isfinite(fit.hill) or fit.hill <= 0 or not 0 < frac < 1:
        return float("nan")
    return float(fit.ec50 * (frac / (1 - frac)) ** (1.0 / fit.hill))


def ec_drop(fit: FourPL, drop: float) -> float:
    """Concentration at which the response is ``drop`` below the top (inf if never)."""
    target = fit.top - drop
    if not np.isfinite(fit.top) or fit.bottom >= target:
        return float("inf")
    return fit.inverse(target)


def median_of_units(
    dose: np.ndarray, response: np.ndarray, cluster: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Per-dose median over units of each unit's mean response (dose > 0)."""
    import pandas as pd

    df = pd.DataFrame({"d": dose, "y": response, "u": np.asarray(cluster).astype(str)})
    df = df[(df["d"] > 0) & np.isfinite(df["y"])]
    m = df.groupby(["u", "d"])["y"].mean().groupby(level="d").median()
    return m.index.to_numpy(float), m.to_numpy(float)


@dataclass
class HillBootstrap:
    params: np.ndarray  # (n_ok, 4): bottom, top, log10_ec50, hill
    n_failed: int

    def fits(self) -> list[FourPL]:
        return [FourPL(*map(float, p)) for p in self.params]

    def interval(self, fn: Callable[[FourPL], float], level: float = 0.95) -> tuple[float, float]:
        """Percentile interval of ``fn(FourPL)`` over the draws (inf allowed)."""
        v = np.array([fn(f) for f in self.fits()], float)
        v = v[~np.isnan(v)]
        if v.size == 0:
            return (float("nan"), float("nan"))
        a = (1 - level) / 2
        return float(np.quantile(v, a)), float(np.quantile(v, 1 - a))


def bootstrap_4pl(
    dose: np.ndarray,
    response: np.ndarray,
    n_boot: int = 200,
    seed: int = 0,
    cluster: np.ndarray | None = None,
) -> HillBootstrap:
    """Refit ``n_boot`` times on points resampled with replacement within each dose, or,
    with ``cluster``, on the median-of-units curve of resampled units."""
    d = np.asarray(dose, float)
    y = np.asarray(response, float)
    if cluster is not None:
        return _bootstrap_units(d, y, np.asarray(cluster).astype(str), n_boot, seed)
    ok = (d > 0) & np.isfinite(y)
    d, y = d[ok], y[ok]
    groups = [np.flatnonzero(d == u) for u in np.unique(d)]
    rng = np.random.default_rng(seed)
    point = fit_4pl(d, y)
    warm = point.fit if np.isfinite(point.fit.hill) else None
    out, failed = [], 0
    for _ in range(n_boot):
        idx = np.concatenate([rng.choice(g, g.size, replace=True) for g in groups])
        f = fit_4pl(d[idx], y[idx], start=warm)
        if f.converged and np.isfinite(f.fit.hill):
            out.append([f.fit.bottom, f.fit.top, f.fit.log10_ec50, f.fit.hill])
        else:
            failed += 1
    return HillBootstrap(np.array(out).reshape(-1, 4), failed)


def _bootstrap_units(
    d: np.ndarray, y: np.ndarray, u: np.ndarray, n_boot: int, seed: int
) -> HillBootstrap:
    units = np.unique(u)
    idx_of = {k: np.flatnonzero(u == k) for k in units}
    point = fit_4pl(*median_of_units(d, y, u))
    warm = point.fit if np.isfinite(point.fit.hill) else None
    rng = np.random.default_rng(seed)
    out, failed = [], 0
    for _ in range(n_boot):
        pick = rng.choice(units, units.size, replace=True)
        idx = np.concatenate([idx_of[k] for k in pick])
        lab = np.concatenate([np.full(idx_of[k].size, f"{j}") for j, k in enumerate(pick)])
        f = fit_4pl(*median_of_units(d[idx], y[idx], lab), start=warm)
        if f.converged and np.isfinite(f.fit.hill):
            out.append([f.fit.bottom, f.fit.top, f.fit.log10_ec50, f.fit.hill])
        else:
            failed += 1
    return HillBootstrap(np.array(out).reshape(-1, 4), failed)
