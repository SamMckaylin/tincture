"""Dose-response containers.

* :class:`FourPL` — parameters of a 4-parameter log-logistic curve.
* :class:`DoseResponse` — one compound × assay context: raw points + optional fit.
* :class:`DoseResponseSet` — many curves in two tables (points, fits), as produced by
  screen loaders (PRISM, GDSC). ``.get()`` returns a :class:`DoseResponse`.

Fitting itself (``tincture.stats.hill``) arrives with the Window plot; loaders can attach
a provider's fit (e.g. PRISM's drc fits) with ``fit_source`` recording where it came from.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from tincture.errors import TinctureDataError


@dataclass(frozen=True)
class FourPL:
    """y = bottom + (top - bottom) / (1 + (x / ec50) ** hill), x in molar.

    ``top`` is the response at zero dose and ``bottom`` at infinite dose (for an
    inhibitor with ``hill > 0``). This matches drc's LL.4 ``d`` (upper) and ``c`` (lower).
    """

    bottom: float
    top: float
    log10_ec50: float
    hill: float

    @property
    def ec50(self) -> float:
        return float(10.0**self.log10_ec50)

    def __call__(self, x: np.ndarray | float) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            return self.bottom + (self.top - self.bottom) / (1.0 + (x / self.ec50) ** self.hill)

    def inverse(self, y: float) -> float:
        """Concentration giving response ``y``; NaN if the curve never reaches it."""
        lo, hi = sorted((self.bottom, self.top))
        if not lo < y < hi or self.hill == 0:
            return float("nan")
        r = (self.top - self.bottom) / (y - self.bottom) - 1.0
        return float(self.ec50 * r ** (1.0 / self.hill))


#: Fit-quality thresholds used to flag curves (documented in methods text).
FIT_MIN_R2 = 0.8
FIT_ASYMPTOTE_RANGE = (-0.5, 1.5)


@dataclass
class DoseResponse:
    """Raw points and (optionally) a fitted curve for one compound in one assay context."""

    compound_id: str
    dose: np.ndarray
    response: np.ndarray
    replicate: np.ndarray | None = None
    assay: str = ""
    context: dict[str, Any] = field(default_factory=dict)  # cell line, readout, duration ...
    fit: FourPL | None = None
    r2: float = np.nan
    converged: bool | None = None
    ec50_ci: tuple[float, float] | None = None
    fit_source: str = ""
    response_label: str = "viability (fraction of vehicle)"

    def __post_init__(self) -> None:
        self.dose = np.asarray(self.dose, dtype=float)
        self.response = np.asarray(self.response, dtype=float)
        if self.dose.shape != self.response.shape or self.dose.ndim != 1:
            raise TinctureDataError("dose and response must be 1-D arrays of equal length.")
        if np.any(self.dose < 0) or np.any(self.dose > 1):
            raise TinctureDataError(
                "DoseResponse doses must be molar in [0, 1].",
                fix="Convert to molar (µM × 1e-6, nM × 1e-9).",
            )
        if self.replicate is not None:
            self.replicate = np.asarray(self.replicate).astype(str)
            if self.replicate.shape != self.dose.shape:
                raise TinctureDataError("replicate must have one entry per point.")

    @property
    def n_doses(self) -> int:
        return int(np.unique(self.dose[self.dose > 0]).size)

    @property
    def dose_range(self) -> tuple[float, float]:
        d = self.dose[self.dose > 0]
        return (float(d.min()), float(d.max())) if d.size else (np.nan, np.nan)

    def quality_flags(self) -> list[str]:
        """Reasons this curve's fit should not be trusted (empty list = usable)."""
        flags: list[str] = []
        if self.fit is None:
            return ["no fit"]
        if self.converged is False:
            flags.append("did not converge")
        if not np.isfinite(self.r2) or self.r2 < FIT_MIN_R2:
            flags.append(
                f"R² {self.r2:.2f} < {FIT_MIN_R2}" if np.isfinite(self.r2) else "R² missing"
            )
        lo, hi = FIT_ASYMPTOTE_RANGE
        for name in ("top", "bottom"):
            v = getattr(self.fit, name)
            if not (np.isfinite(v) and lo <= v <= hi):
                flags.append(f"{name} asymptote {v:.3g} outside [{lo}, {hi}]")
        dmin, dmax = self.dose_range
        if np.isfinite(dmin) and not (dmin <= self.fit.ec50 <= dmax):
            flags.append("EC50 outside tested range")
        return flags

    @property
    def fit_ok(self) -> bool:
        return not self.quality_flags()


#: Columns of DoseResponseSet.points and .fits.
POINT_COLS = ("compound_id", "context_id", "dose", "replicate", "response")
FIT_COLS = ("compound_id", "context_id", "bottom", "top", "log10_ec50", "hill", "r2")


@dataclass
class DoseResponseSet:
    """Many dose-response curves as two tables.

    ``points``: one row per measurement with :data:`POINT_COLS` (+ any metadata).
    ``fits``: one row per curve with :data:`FIT_COLS` (+ provider extras, e.g. auc, ic50).
    ``contexts``: one row per ``context_id`` (e.g. cell line) with descriptive columns.
    """

    points: pd.DataFrame
    fits: pd.DataFrame | None = None
    contexts: pd.DataFrame | None = None
    compounds: pd.DataFrame | None = None
    assay: str = ""
    fit_source: str = ""
    source: str = ""

    def __post_init__(self) -> None:
        miss = [c for c in POINT_COLS if c not in self.points.columns]
        if miss:
            raise TinctureDataError(f"DoseResponseSet.points missing columns {miss}.")
        if self.fits is not None:
            miss = [c for c in FIT_COLS if c not in self.fits.columns]
            if miss:
                raise TinctureDataError(f"DoseResponseSet.fits missing columns {miss}.")
        if (self.points["dose"] < 0).any() or (self.points["dose"] > 1).any():
            raise TinctureDataError("DoseResponseSet doses must be molar in [0, 1].")

    def __len__(self) -> int:
        return int(self.points[["compound_id", "context_id"]].drop_duplicates().shape[0])

    def keys(self) -> pd.DataFrame:
        return self.points[["compound_id", "context_id"]].drop_duplicates().reset_index(drop=True)

    def get(self, compound_id: str, context_id: str) -> DoseResponse:
        p = self.points
        sel = p[(p["compound_id"] == compound_id) & (p["context_id"] == context_id)]
        if sel.empty:
            raise TinctureDataError(
                f"No points for compound {compound_id!r} in context {context_id!r}.",
                fix="See DoseResponseSet.keys() for available pairs.",
            )
        fit = None
        r2 = np.nan
        extras: dict[str, Any] = {}
        if self.fits is not None:
            f = self.fits
            row = f[(f["compound_id"] == compound_id) & (f["context_id"] == context_id)]
            if not row.empty:
                r = row.iloc[0]
                fit = FourPL(
                    float(r["bottom"]), float(r["top"]), float(r["log10_ec50"]), float(r["hill"])
                )
                r2 = float(r["r2"])
                extras = {k: r[k] for k in row.columns if k not in FIT_COLS}
        context: dict[str, Any] = {"context_id": context_id, **extras}
        if self.contexts is not None and context_id in self.contexts.index:
            context.update({str(k): v for k, v in self.contexts.loc[context_id].to_dict().items()})
        return DoseResponse(
            compound_id=compound_id,
            dose=sel["dose"].to_numpy(),
            response=sel["response"].to_numpy(),
            replicate=sel["replicate"].to_numpy(),
            assay=self.assay,
            context=context,
            fit=fit,
            r2=r2,
            fit_source=self.fit_source,
        )
