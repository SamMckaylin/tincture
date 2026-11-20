"""Empirical-Bayes moderation of gene-wise variances (Smyth 2004, limma ``squeezeVar``).

Given residual variances ``s2`` with ``df`` degrees of freedom, fit a scaled inverse-χ²
prior (``s0²``, ``d0``) and return posterior variances
``(d0·s0² + df·s2) / (d0 + df)``. With ``covariate`` (e.g. average log-expression) the
prior location follows a lowess trend (limma-trend).

Reference: Smyth GK. Linear models and empirical Bayes methods for assessing differential
expression in microarray experiments. Stat Appl Genet Mol Biol 2004;3:Article3.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import digamma, polygamma


def trigamma(x: np.ndarray | float) -> np.ndarray:
    return np.asarray(polygamma(1, x))


def trigamma_inverse(x: float) -> float:
    """Solve trigamma(y) = x for y > 0 (Newton iteration as in limma::trigammaInverse)."""
    if x > 1e7:
        return 1.0 / np.sqrt(x)
    if x < 1e-6:
        return 1.0 / x
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = float(trigamma(y))
        dif = tri * (1 - tri / x) / float(polygamma(2, y))
        y += dif
        if -dif / y < 1e-8:
            break
    return float(y)


@dataclass(frozen=True)
class FDistFit:
    """Prior for the variances: ``s0²`` (scale, per gene if trended) and ``d0`` (df)."""

    s0_sq: np.ndarray
    d0: float


def fit_f_dist(
    s2: np.ndarray, df: float | np.ndarray, covariate: np.ndarray | None = None
) -> FDistFit:
    """Moment estimation of the scaled-F prior (limma ``fitFDist``)."""
    s2 = np.asarray(s2, dtype=float)
    df_arr = np.broadcast_to(np.asarray(df, dtype=float), s2.shape)
    ok = np.isfinite(s2) & (s2 > 0) & (df_arr > 0)
    if ok.sum() < 3:
        raise ValueError("Need at least 3 positive, finite variances to fit the prior.")
    x = np.maximum(s2[ok], 1e-5 * np.median(s2[ok]))
    d1 = df_arr[ok]
    z = np.log(x)
    e = z - digamma(d1 / 2) + np.log(d1 / 2)
    if covariate is None:
        emean = np.full_like(e, e.mean())
    else:
        from statsmodels.nonparametric.smoothers_lowess import lowess

        cov = np.asarray(covariate, dtype=float)[ok]
        frac = min(1.0, max(0.1, 200 / ok.sum()))
        emean = lowess(e, cov, frac=frac, return_sorted=False)
    evar = np.sum((e - emean) ** 2) / (e.size - (1 if covariate is None else 4))
    evar -= float(np.mean(trigamma(d1 / 2)))
    if evar > 0:
        d0 = 2 * trigamma_inverse(evar)
        s0_ok = np.exp(emean + digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0 = np.inf
        s0_ok = np.exp(emean)
    s0 = np.full(s2.shape, np.nan)
    s0[ok] = s0_ok
    if covariate is not None and (~ok).any():
        cov_all = np.asarray(covariate, dtype=float)
        order = np.argsort(cov_all[ok])
        s0[~ok] = np.interp(cov_all[~ok], cov_all[ok][order], s0_ok[order])
    elif (~ok).any():
        s0[~ok] = float(np.median(s0_ok))
    return FDistFit(s0_sq=s0, d0=float(d0))


@dataclass(frozen=True)
class SqueezedVar:
    var_post: np.ndarray
    df_prior: float
    var_prior: np.ndarray
    df_total: np.ndarray


def squeeze_var(
    s2: np.ndarray, df: float | np.ndarray, covariate: np.ndarray | None = None
) -> SqueezedVar:
    """Posterior (moderated) variances. ``df_total = df + d0`` (capped as in limma)."""
    s2 = np.asarray(s2, dtype=float)
    fit = fit_f_dist(s2, df, covariate)
    df_arr = np.broadcast_to(np.asarray(df, dtype=float), s2.shape)
    if np.isinf(fit.d0):
        post = fit.s0_sq.copy()
        df_total = np.full(s2.shape, np.inf)
    else:
        post = (fit.d0 * fit.s0_sq + df_arr * np.nan_to_num(s2)) / (fit.d0 + df_arr)
        df_total = np.minimum(df_arr + fit.d0, df_arr.sum() if df_arr.ndim else np.inf)
    return SqueezedVar(var_post=post, df_prior=fit.d0, var_prior=fit.s0_sq, df_total=df_total)
