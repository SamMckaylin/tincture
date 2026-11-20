"""Dose-path statistics: 2-D positions per dose and the bend test.

Axes (defaults agreed with the user):
* **x — reference likeness:** scalar projection of the dose's signature on the unit vector of
  the compound's top reference (log2 fold-change units along that direction).
* **y — stress:** mean log2 fold change of a stress gene set (default: union of MSigDB
  Hallmark APOPTOSIS, P53_PATHWAY and UNFOLDED_PROTEIN_RESPONSE) minus the mean over all
  genes; or **cell loss** = −log2(cells per well ÷ vehicle cells per well).

**Bend test.** The path starts at vehicle (0, 0) and visits doses in order. Angles depend on
how the two axes are scaled, so they are measured after dividing each axis by its spread
across every signature in the stratum (``PathAxes.scale``; standard deviations), not in raw
units. At each interior point the turning angle θ between successive dose deltas is
recomputed under posterior draws of every dose's position. Both axes are linear in the
signature, so with independent normal per-feature errors each dose's (x, y) is exactly
bivariate normal (variances Σ u² sd², Σ w² sd², covariance Σ u w sd², where u and w are the
x and y weight vectors); draws are taken from that, not feature by feature. A *bend*
is a turn larger than ``MIN_BEND_DEG`` (30°): p = P(|θ| ≤ 30°) under the draws (with the usual
+1 correction), Bonferroni-corrected over the turning points. Only bends with corrected
p < 0.05 are reported, so a small but precisely measured curve is not called a bend and a
large turn within noise is not either.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from tincture.errors import TinctureDataError
from tincture.model.signature import Signature

MIN_BEND_DEG = 30.0

STRESS_SETS = (
    "HALLMARK_APOPTOSIS",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
)


def stress_gene_set(hallmark_table, sets: Sequence[str] = STRESS_SETS) -> set[str]:
    """Union of the stress Hallmark sets (``source``/``target`` table, names with or
    without the ``HALLMARK_`` prefix)."""
    names = {s.replace("HALLMARK_", "") for s in sets}
    src = hallmark_table["source"].astype(str).str.replace("HALLMARK_", "", regex=False)
    genes = set(hallmark_table.loc[src.isin(names), "target"].astype(str))
    if not genes:
        raise TinctureDataError(
            f"None of the stress sets {sorted(names)} found in the gene-set table."
        )
    return genes


@dataclass
class PathAxes:
    """How to turn an effect vector into (x, y)."""

    ref_unit: np.ndarray  # unit vector of the top reference (aligned features)
    stress_mask: np.ndarray | None  # bool over features (stress genes); None → cell-loss axis
    scale: tuple[float, float] = (1.0, 1.0)  # axis spreads used to measure angles

    def x(self, eff: np.ndarray) -> np.ndarray:
        e = np.nan_to_num(eff)
        return e @ self.ref_unit

    @property
    def y_weights(self) -> np.ndarray | None:
        """y as a linear map: y = eff @ w (stress mean minus overall mean)."""
        if self.stress_mask is None:
            return None
        m = np.asarray(self.stress_mask, bool)
        return m / m.sum() - 1.0 / m.size

    def y(self, eff: np.ndarray) -> np.ndarray | float:
        if self.stress_mask is None:
            raise ValueError("stress axis not defined (cell-loss axis in use)")
        e = np.nan_to_num(eff)
        if e.ndim == 1:
            return float(e[self.stress_mask].mean() - e.mean())
        return e[:, self.stress_mask].mean(axis=1) - e.mean(axis=1)


@dataclass
class DosePath:
    compound_id: str
    doses: np.ndarray  # molar, ascending (vehicle excluded)
    x: np.ndarray
    y: np.ndarray
    turn_p: np.ndarray  # Bonferroni-corrected p per interior point (len = n_doses - 1)
    turn_angle: np.ndarray  # observed signed turning angle (radians, scaled axes) per point
    draws: np.ndarray | None = None  # (n_boot, n_doses, 2)

    @property
    def n_doses(self) -> int:
        return int(self.doses.size)

    @property
    def bend_index(self) -> int | None:
        """Dose index where the path bends significantly (smallest p < 0.05), else None."""
        if self.turn_p.size == 0:
            return None
        i = int(np.argmin(self.turn_p))
        return i + 1 if self.turn_p[i] < 0.05 else None

    @property
    def bend_p(self) -> float:
        return float(self.turn_p.min()) if self.turn_p.size else float("nan")


def signed_turns(xy: np.ndarray) -> np.ndarray:
    """Signed turning angles at interior points of a path that starts at the origin.

    ``xy``: (..., n_points, 2) for the doses; the origin (vehicle) is prepended.
    Returns (..., n_points - 1).
    """
    origin = np.zeros((*xy.shape[:-2], 1, 2))
    pts = np.concatenate([origin, xy], axis=-2)
    seg = np.diff(pts, axis=-2)
    a, b = seg[..., :-1, :], seg[..., 1:, :]
    cross = a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]
    dot = (a * b).sum(axis=-1)
    return np.arctan2(cross, dot)


def dose_path(
    sigs: Sequence[Signature],
    axes: PathAxes,
    y_values: np.ndarray | None = None,
    n_boot: int = 500,
    seed: int = 0,
) -> DosePath:
    """Positions and bend test for one compound's signatures (any order; sorted by dose).

    ``y_values`` overrides the y axis (cell-loss axis): one value per signature, no draws on
    y (cell counts carry their own, smaller uncertainty, not modelled here).
    """
    order = np.argsort([s.dose for s in sigs])
    sigs = [sigs[i] for i in order]
    doses = np.array([s.dose for s in sigs], dtype=float)
    if np.any(doses <= 0):
        raise TinctureDataError("Dose paths take treated signatures only (dose > 0).")
    E = np.vstack([np.nan_to_num(s.effect) for s in sigs])
    x = np.array([axes.x(e) for e in E])
    y = (
        np.asarray(y_values, float)[order]
        if y_values is not None
        else np.array([axes.y(e) for e in E])
    )
    k = len(sigs)
    if k < 3:
        return DosePath(sigs[0].compound_id, doses, x, y, np.array([]), np.array([]))
    rng = np.random.default_rng(seed)
    draws = np.empty((n_boot, k, 2))
    u = np.asarray(axes.ref_unit, float)
    w = axes.y_weights if y_values is None else None
    for j, s in enumerate(sigs):
        v = np.nan_to_num(s.se) ** 2 if s.se is not None else np.zeros_like(E[j])
        vx = float(v @ u**2)
        vy = float(v @ w**2) if w is not None else 0.0
        cxy = float(v @ (u * w)) if w is not None else 0.0
        L = np.linalg.cholesky(np.array([[vx, cxy], [cxy, vy]]) + 1e-18 * np.eye(2))
        draws[:, j, :] = np.array([x[j], y[j]]) + rng.standard_normal((n_boot, 2)) @ L.T
    sc = np.asarray(axes.scale, float)
    if np.any(~np.isfinite(sc)) or np.any(sc <= 0):
        raise TinctureDataError("Axis scales must be positive to measure bend angles.")
    obs = signed_turns(np.column_stack([x, y]) / sc)
    th = signed_turns(draws / sc)
    n = th.shape[0]
    small = np.sum(np.abs(th) <= np.radians(MIN_BEND_DEG), axis=0)
    p = np.minimum(1.0, (small + 1) / (n + 1) * th.shape[1])  # Bonferroni over turns
    return DosePath(sigs[0].compound_id, doses, x, y, p, obs, draws)
