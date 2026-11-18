"""Combination effect models: Bliss (default), HSA and ZIP, with honest uncertainty.

All responses are viabilities (1 = vehicle). **Excess = expected viability − observed
viability**, i.e. extra inhibition beyond the reference model; positive = better than
expected (synergy for Bliss/ZIP, benefit over the best single agent for HSA).

Reference models (single-agent viabilities vA, vB at the same concentrations):

* **Bliss independence** (default): expected viability = vA · vB (expected inhibition
  = iA + iB − iA·iB). Single agents are taken from 4PL fits of their own wells
  (``single_agent="fitted"``, as GDSC does) or from observed means (``"observed"``).
* **HSA** (highest single agent): expected viability = min(vA, vB).
* **ZIP** (zero interaction potency; Yadav et al., Comput Struct Biotechnol J 2015): the
  expected surface is Bliss on *fitted* single agents. The observed surface is replaced by
  curve fits: for each level of drug B, drug A's dose-response is refitted with drug B's
  single-agent effect as its baseline, ``y = (yB + E·(x/m)^n) / (1 + (x/m)^n)`` (inhibition),
  and vice versa; the two fits are averaged. δ = averaged fit − expected (inhibition), so δ
  is smoothed by construction and needs ≥ 4 levels of each drug.

**Ceiling mask:** cells whose expected inhibition exceeds 90% are masked: there is no room
left to show an excess, so excess there is meaningless (handover §5.11).

**Uncertainty** (per cell, then Benjamini–Hochberg across unmasked cells):

* ≥ 2 replicate matrices: the excess is computed **within each replicate** (single agents
  refitted per replicate when every replicate has its own single-agent wells; otherwise the
  pooled single agents are used and the interval ignores their uncertainty, noted), and a
  one-sample t-test / t-interval across replicates is used. With 2–3 replicates the
  intervals are wide, which is the honest answer; a bootstrap over so few replicates is
  overconfident.
* one replicate, with ``noise_sd`` given (well-level SD of viability, e.g. the residual SD
  of single-agent fits): parametric bootstrap adding N(0, noise_sd) to every well; two-sided
  p = (2·min(#≤0, #≥0) + 1) / (B + 1), so the smallest attainable p is 1/(B + 1)
  (use B ≥ 200 so that BH over ~50 cells can reach q < 0.05);
* one replicate with ``noise_sd="matrix"``: the well-level SD is estimated from the matrix
  itself (:func:`matrix_noise_sd`, second differences along rows and columns of unsaturated
  wells), then the same parametric bootstrap;
* otherwise: no interval and no significance claim (the caption says so).

A cell is in the **significant region** when q < ``alpha`` **and** excess ≥ ``min_excess``
(default 0.1, i.e. 10% of vehicle viability), so that a statistically clear but trivial
excess is not outlined.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import pandas as pd
from scipy import stats as sps
from scipy.optimize import least_squares
from statsmodels.stats.multitest import multipletests

from tincture.errors import TinctureDataError
from tincture.model.combo import ComboMatrix
from tincture.model.dose_response import FourPL
from tincture.stats.hill import fit_4pl

Model = Literal["bliss", "hsa", "zip"]
MODELS: tuple[str, ...] = ("bliss", "hsa", "zip")
MODEL_NAMES = {"bliss": "Bliss independence", "hsa": "highest single agent (HSA)", "zip": "ZIP"}
CEILING = 0.9  # expected inhibition above which excess is masked
MIN_EXCESS = 0.1
NOISE_RANGE = (0.1, 0.9)  # viability range of wells used to estimate matrix noise
NOISE_MIN_TRIPLETS = 10


def grid_noise_sd(
    grid: np.ndarray,
    lo: float = NOISE_RANGE[0],
    hi: float = NOISE_RANGE[1],
    min_triplets: int = NOISE_MIN_TRIPLETS,
) -> float | None:
    """Well-level SD from one concentration grid: second differences y₋₁ − 2y₀ + y₊₁ along
    rows and columns (removes linear dose trends), over triplets whose three viabilities all
    lie in (lo, hi) (saturated wells understate noise); robust scale
    σ = 1.4826 · median|Δ²| / √6. None if fewer than ``min_triplets`` triplets qualify."""
    d: list[float] = []
    for arr in (grid, grid.T):
        for row in arr:
            for k in range(1, len(row) - 1):
                t = row[k - 1 : k + 2]
                if np.isfinite(t).all() and ((t > lo) & (t < hi)).all():
                    d.append(float(t[0] - 2 * t[1] + t[2]))
    if len(d) < min_triplets:
        return None
    return float(1.4826 * np.median(np.abs(d)) / np.sqrt(6))


def matrix_noise_sd(cm: ComboMatrix) -> float | None:
    """Median over replicate matrices of :func:`grid_noise_sd` (single-well SD)."""
    est = [grid_noise_sd(cm.grid(g)) for _, g in cm.data.groupby("replicate")]
    vals = [e for e in est if e is not None]
    return float(np.median(vals)) if vals else None


@dataclass
class SynergyResult:
    model: str
    levels_a: np.ndarray
    levels_b: np.ndarray
    observed: np.ndarray  # (nb, na) mean viability
    expected: np.ndarray  # (nb, na) expected viability under the model
    excess: np.ndarray  # (nb, na) expected − observed viability (ZIP: δ)
    masked: np.ndarray  # bool (nb, na): ceiling
    lo: np.ndarray
    hi: np.ndarray
    p: np.ndarray
    q: np.ndarray
    significant: np.ndarray  # bool: q < alpha and excess ≥ min_excess
    antagonistic: np.ndarray  # bool: q < alpha and excess ≤ −min_excess
    uncertainty: Literal["replicates", "noise_sd", "none"]
    n_replicates: int
    noise_sd: float | None
    single_agent: str
    fit_a: FourPL | None
    fit_b: FourPL | None
    alpha: float = 0.05
    min_excess: float = MIN_EXCESS
    noise_method: str = ""
    notes: list[str] = field(default_factory=list)

    @property
    def score(self) -> float:
        """Mean excess over unmasked measured cells."""
        v = self.excess[~self.masked & np.isfinite(self.excess)]
        return float(v.mean()) if v.size else float("nan")

    @property
    def n_significant(self) -> int:
        return int(self.significant.sum())


# ---------------------------------------------------------------------------------------
# single agents
# ---------------------------------------------------------------------------------------
def _single_values(
    s: pd.DataFrame, levels: np.ndarray, how: str, start: FourPL | None = None
) -> tuple[np.ndarray, FourPL | None, str | None]:
    """Single-agent viability at ``levels`` (fitted or observed means)."""
    if s.empty:
        raise TinctureDataError(
            "No single-agent wells for one of the drugs.",
            fix="Include rows with the other drug at concentration 0.",
        )
    obs = s.groupby("conc")["viability"].mean()
    if how == "fitted":
        h = fit_4pl(s["conc"].to_numpy(), s["viability"].to_numpy(), start=start)
        if h.converged and np.isfinite(h.fit.top):
            return np.asarray(h.fit(levels), float), h.fit, None
        note = "single-agent fit failed; observed means used"
    else:
        note = None
    vals = obs.reindex(levels).to_numpy(float)
    if np.isnan(vals).any():
        raise TinctureDataError(
            "Observed single-agent values are missing at some combination concentrations.",
            fix="Use single_agent='fitted' or measure single agents at every level.",
        )
    return vals, None, note


def _zip_fit_row(
    x: np.ndarray, y: np.ndarray, base: float, start: FourPL | None
) -> np.ndarray | None:
    """Fit y = (base + E·(x/m)^n)/(1 + (x/m)^n) (inhibition) at concentrations x."""
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return None
    lx = np.log10(x[ok])
    yy = y[ok]

    def f(p: np.ndarray, lxx: np.ndarray) -> np.ndarray:
        e, lm, n = p
        r = 10.0 ** (n * (lxx - lm))
        return (base + e * r) / (1 + r)

    lm0 = start.log10_ec50 if start is not None and np.isfinite(start.log10_ec50) else lx.mean()
    n0 = start.hill if start is not None and np.isfinite(start.hill) else 1.0
    lb = [min(base, 0.0) - 1e-6, lx.min() - 2, 0.1]
    ub = [1.5, lx.max() + 2, 10.0]
    p0 = np.clip([max(base, yy.max()), lm0, n0], np.array(lb) + 1e-9, np.array(ub) - 1e-9)
    r = least_squares(lambda p: f(p, lx) - yy, p0, bounds=(lb, ub), loss="soft_l1", f_scale=0.1)
    return np.asarray(f(r.x, np.log10(x)), float)


def _point(
    cm: ComboMatrix,
    data: pd.DataFrame,
    model: str,
    single_agent: str,
    starts: tuple[FourPL | None, FourPL | None] = (None, None),
    singles: pd.DataFrame | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, FourPL | None, FourPL | None, list[str]]:
    """Observed, expected and excess grids for one data set. ``singles`` overrides the
    single-agent wells (pooled singles for per-replicate combination grids)."""
    la, lb = cm.levels_a, cm.levels_b
    src = data if singles is None else singles
    sub = ComboMatrix(src, cm.drug_a, cm.drug_b) if src is not cm.data else cm
    how = "fitted" if model == "zip" else single_agent
    va, fa, na = _single_values(sub.single("a"), la, how, starts[0])
    vb, fb, nb_ = _single_values(sub.single("b"), lb, how, starts[1])
    notes = [n for n in (na, nb_) if n]
    obs = cm.grid(data)
    VA, VB = np.meshgrid(va, vb)  # (nb, na)
    if model == "hsa":
        exp = np.minimum(VA, VB)
        return obs, exp, exp - obs, fa, fb, notes
    exp = VA * VB
    if model == "bliss":
        return obs, exp, exp - obs, fa, fb, notes
    # ZIP on inhibition
    if fa is None or fb is None:
        raise TinctureDataError("ZIP needs converged single-agent fits for both drugs.")
    ia, ib = 1 - va, 1 - vb
    inh = 1 - obs
    row_fit = np.full_like(obs, np.nan)  # vary A at each B level
    col_fit = np.full_like(obs, np.nan)  # vary B at each A level
    for j in range(len(lb)):
        r = _zip_fit_row(la, inh[j, :], ib[j], fa)
        if r is not None:
            row_fit[j, :] = r
    for i in range(len(la)):
        r = _zip_fit_row(lb, inh[:, i], ia[i], fb)
        if r is not None:
            col_fit[:, i] = r
    yc = np.nanmean(np.stack([row_fit, col_fit]), axis=0)
    delta = yc - (1 - exp)  # inhibition: fitted − expected
    return obs, exp, delta, fa, fb, notes


def synergy(
    cm: ComboMatrix,
    model: Model = "bliss",
    single_agent: Literal["fitted", "observed"] = "fitted",
    noise_sd: float | Literal["matrix"] | None = None,
    n_boot: int = 200,
    alpha: float = 0.05,
    min_excess: float = MIN_EXCESS,
    ceiling: float = CEILING,
    seed: int = 0,
) -> SynergyResult:
    """Excess over a reference model on a full concentration matrix."""
    if model not in MODELS:
        raise TinctureDataError(f"Unknown synergy model {model!r}.", fix=f"Use one of {MODELS}.")
    if cm.design == "anchored":
        raise TinctureDataError(
            f"{cm.drug_a} + {cm.drug_b} is an anchored design (one drug at ≤ 2 "
            "concentrations); a combination surface cannot be estimated from it.",
            fix="Use tc.combination_lift for anchored screens.",
        )
    if cm.design != "full":
        raise TinctureDataError(
            "The combination matrix is incomplete (not every concentration pair measured).",
            fix="Provide a full grid (e.g. 7 × 7) with single-agent wells.",
        )
    if model == "zip" and min(len(cm.levels_a), len(cm.levels_b)) < 4:
        raise TinctureDataError("ZIP needs at least 4 concentrations of each drug.")
    obs, exp, exc, fa, fb, notes = _point(cm, cm.data, model, single_agent)
    masked = (1 - exp) > ceiling
    reps = cm.replicates
    rng = np.random.default_rng(seed)
    shape = obs.shape
    lo = np.full(shape, np.nan)
    hi = np.full(shape, np.nan)
    p = np.full(shape, np.nan)
    kind: Literal["replicates", "noise_sd", "none"] = "none"
    if len(reps) >= 2:
        kind = "replicates"
        per = []
        own_singles = (
            all(
                g[(g["conc_a"] > 0) & (g["conc_b"] == 0)]["conc_a"].nunique() >= 4
                and g[(g["conc_b"] > 0) & (g["conc_a"] == 0)]["conc_b"].nunique() >= 4
                for _, g in cm.data.groupby("replicate")
            )
            or single_agent == "observed"
        )
        if not own_singles:
            notes.append(
                "replicates lack their own single-agent wells: pooled single agents used, "
                "so intervals ignore single-agent uncertainty"
            )
        for _, g in cm.data.groupby("replicate"):
            try:
                per.append(
                    _point(
                        cm, g, model, single_agent, (fa, fb),
                        None if own_singles else cm.data,
                    )[2]
                )  # fmt: skip
            except TinctureDataError:
                continue
        E = np.stack(per) if per else np.empty((0, *shape))
        n = np.isfinite(E).sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            m = np.nanmean(E, axis=0) if len(E) else np.full(shape, np.nan)
            se = np.nanstd(E, axis=0, ddof=1) / np.sqrt(n) if len(E) >= 2 else np.nan * m
            tcrit = sps.t.ppf(0.975, np.maximum(n - 1, 1))
            lo = exc - tcrit * se
            hi = exc + tcrit * se
            t = m / se
            p = 2 * sps.t.sf(np.abs(t), np.maximum(n - 1, 1))
        p = np.where((n >= 2) & np.isfinite(p), p, np.nan)
        p = np.where((n >= 2) & (se == 0), 0.0, p)
    noise_method = ""
    if len(reps) < 2 and noise_sd == "matrix":
        est = matrix_noise_sd(cm)
        if est is None:
            notes.append(
                "too few unsaturated wells to estimate the matrix noise; no significance claimed"
            )
            noise_sd = None
        else:
            noise_sd = est
            noise_method = "estimated from the matrix (second differences)"
    elif isinstance(noise_sd, float | int):
        noise_method = "stated"
    if len(reps) < 2 and isinstance(noise_sd, float | int) and noise_sd > 0:
        kind = "noise_sd"
        boots = []
        for _ in range(n_boot):
            d = cm.data.assign(
                viability=cm.data["viability"] + rng.normal(0, noise_sd, len(cm.data))
            )
            try:
                boots.append(_point(cm, d, model, single_agent, (fa, fb))[2])
            except TinctureDataError:
                continue
        if len(boots) >= 20:
            B = np.stack(boots)
            lo, hi = np.nanpercentile(B, [2.5, 97.5], axis=0)
            le = (B <= 0).sum(axis=0)
            ge = (B >= 0).sum(axis=0)
            p = np.minimum(1.0, (2 * np.minimum(le, ge) + 1) / (len(B) + 1))
        else:
            notes.append("too few successful bootstrap refits; no significance claimed")
            kind = "none"
    q = np.full(shape, np.nan)
    sig = np.zeros(shape, bool)
    ant = np.zeros(shape, bool)
    test = ~masked & np.isfinite(exc) & np.isfinite(p)
    if kind != "none" and test.any():
        q[test] = multipletests(p[test], method="fdr_bh")[1]
        sig = test & (q < alpha) & (exc >= min_excess)
        ant = test & (q < alpha) & (exc <= -min_excess)
    return SynergyResult(
        model=model,
        levels_a=cm.levels_a,
        levels_b=cm.levels_b,
        observed=obs,
        expected=exp,
        excess=exc,
        masked=masked,
        lo=lo,
        hi=hi,
        p=p,
        q=q,
        significant=sig,
        antagonistic=ant,
        uncertainty=kind,
        n_replicates=len(reps),
        noise_sd=float(noise_sd) if kind == "noise_sd" else None,  # type: ignore[arg-type]
        single_agent="fitted" if model == "zip" else single_agent,
        fit_a=fa,
        fit_b=fb,
        alpha=alpha,
        min_excess=min_excess,
        noise_method=noise_method if kind == "noise_sd" else "",
        notes=notes,
    )
