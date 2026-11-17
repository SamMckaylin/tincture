"""Projection of signatures onto reference mechanisms (Compass statistics).

* **NNLS fit:** ``y ≈ Σ w_k r̂_k`` with ``w ≥ 0``, references scaled to unit norm so weights
  are comparable.
* **Noise correction.** Over thousands of features, measurement noise alone gives a
  signature a large norm, and noise is by construction "unexplained". With per-feature SEs,
  ``signal² = max(‖y‖² − Σ se², 0)`` is an unbiased estimate of the true signal energy, and
  ``R²_signal = min(1, ‖ŷ‖² / signal²)``: the share of the *signal* the references explain.
  **Residual = 1 − R²_signal** (the novelty signal). Whether there is signal at all is tested
  with ``z = (Σ (y/se)² − p) / √(2p)`` (χ²_p under no effect); novelty is only judged when
  ``z > 3``. Without SEs, the uncorrected R² is used.
* **Compass position:** spoke k sits at angle θ_k. Direction = angle of ``Σ w_k u_k``;
  **radius = noise-corrected magnitude √signal²**, linear, not area-scaled. ``concentration`` =
  ``‖Σ w_k u_k‖ / Σ w_k`` (1 = one spoke; low = a mixture of distant spokes).
* **Uncertainty:** posterior draws of the signature (normal, per-feature SD) are re-fitted;
  the 95% ellipse comes from the draws' covariance.
* **Spoke order:** optimal-leaf-ordered average-linkage clustering of reference cosine
  distances, so correlated references sit next to each other.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import nnls
from scipy.stats import chi2

from tincture.errors import TinctureDataError

SIGNAL_Z = 3.0


@dataclass
class Projection:
    weights: np.ndarray  # one per reference (NNLS, unit-norm references)
    r2: float
    magnitude: float  # noise-corrected when SEs were given
    signal_z: float = float("inf")  # χ² z-score for "any signal at all"

    @property
    def has_signal(self) -> bool:
        return bool(self.signal_z > SIGNAL_Z)

    @property
    def residual(self) -> float:
        return float(1.0 - self.r2)

    @property
    def shares(self) -> np.ndarray:
        s = self.weights.sum()
        return self.weights / s if s > 0 else np.zeros_like(self.weights)


def unit_rows(R: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(R, axis=1, keepdims=True)
    if np.any(n == 0):
        raise TinctureDataError("A reference signature is all zeros.", fix="Drop empty references.")
    return R / n


def nnls_project(y: np.ndarray, R: np.ndarray, se: np.ndarray | None = None) -> Projection:
    """Non-negative least squares of ``y`` on the rows of ``R`` (refs × features).

    With ``se`` (per-feature standard errors of ``y``), magnitude and R² are corrected for
    noise (see module docstring).
    """
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(y) & np.isfinite(R).all(axis=0)
    if se is not None:
        se = np.asarray(se, dtype=float)
        ok &= np.isfinite(se) & (se > 0)
    if ok.sum() < R.shape[0] + 1:
        raise TinctureDataError("Too few shared finite features to project the signature.")
    Ru = unit_rows(R[:, ok])
    yy = y[ok]
    w, rnorm = nnls(Ru.T, yy)
    ss = float(yy @ yy)
    fitted = float(ss - rnorm**2)
    if se is None:
        r2 = 0.0 if ss == 0 else fitted / ss
        return Projection(weights=w, r2=max(0.0, min(1.0, r2)), magnitude=float(np.sqrt(ss)))
    s2 = se[ok] ** 2
    signal2 = max(ss - float(s2.sum()), 0.0)
    p = int(ok.sum())
    z = (float(np.sum(yy**2 / s2)) - p) / np.sqrt(2 * p)
    r2 = min(1.0, fitted / signal2) if signal2 > 0 else 0.0
    return Projection(weights=w, r2=max(0.0, r2), magnitude=float(np.sqrt(signal2)), signal_z=z)


def spoke_angles(n: int, start: float = np.pi / 2) -> np.ndarray:
    """Evenly spaced spoke angles, the first pointing up, going clockwise."""
    return start - np.arange(n) * 2 * np.pi / n


def compass_xy(p: Projection, angles: np.ndarray) -> tuple[float, float, float]:
    """(x, y, concentration) for a projection: direction from weights, radius = magnitude."""
    u = np.column_stack([np.cos(angles), np.sin(angles)])
    v = p.weights @ u
    total = p.weights.sum()
    if total <= 0 or np.hypot(*v) == 0:
        return 0.0, 0.0, 0.0
    d = v / np.hypot(*v)
    return float(d[0] * p.magnitude), float(d[1] * p.magnitude), float(np.hypot(*v) / total)


def bootstrap_ellipse(
    effect: np.ndarray,
    sd: np.ndarray,
    R: np.ndarray,
    angles: np.ndarray,
    n_boot: int = 200,
    level: float = 0.95,
    seed: int = 0,
) -> tuple[float, float, float, float, float]:
    """95% ellipse (cx, cy, width, height, angle_deg) of the compass position under
    posterior draws ``effect + sd · z``."""
    rng = np.random.default_rng(seed)
    pts = np.empty((n_boot, 2))
    for b in range(n_boot):
        draw = effect + np.nan_to_num(sd) * rng.standard_normal(effect.shape)
        x, y, _ = compass_xy(nnls_project(draw, R, sd), angles)
        pts[b] = (x, y)
    return cov_ellipse(pts, level)


def cov_ellipse(pts: np.ndarray, level: float = 0.95) -> tuple[float, float, float, float, float]:
    """Normal-theory ellipse (cx, cy, width, height, angle_deg) covering ``level`` of 2-D
    points (n × 2)."""
    c = pts.mean(axis=0)
    cov = np.cov(pts.T)
    vals, vecs = np.linalg.eigh(cov)
    vals = np.clip(vals, 0, None)
    k = np.sqrt(chi2.ppf(level, 2))
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    ang = float(np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0])))
    return (
        float(c[0]),
        float(c[1]),
        float(2 * k * np.sqrt(vals[0])),
        float(2 * k * np.sqrt(vals[1])),
        ang,
    )


def seriate(sim: np.ndarray) -> np.ndarray:
    """Order references so similar ones are adjacent (optimal leaf ordering)."""
    from scipy.cluster.hierarchy import leaves_list, linkage, optimal_leaf_ordering
    from scipy.spatial.distance import squareform

    n = sim.shape[0]
    if n <= 2:
        return np.arange(n)
    d = np.clip(1 - sim, 0, 2)
    np.fill_diagonal(d, 0)
    cond = squareform(d, checks=False)
    Z = optimal_leaf_ordering(linkage(cond, method="average"), cond)
    return leaves_list(Z)
