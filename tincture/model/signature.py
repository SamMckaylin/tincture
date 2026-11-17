"""Signature: one perturbation's effect vector with uncertainty."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from tincture.errors import TinctureDataError


@dataclass
class Signature:
    """Effect of one perturbation (compound × dose × time × cell line) on each feature.

    Bootstrap draws are stored (not just point estimates) so downstream confidence
    intervals are honest. ``effect`` is usually a (shrunk) log fold change vs vehicle.
    """

    features: np.ndarray
    effect: np.ndarray
    se: np.ndarray | None = None
    boot: np.ndarray | None = None  # shape (n_boot, n_features)
    compound_id: str = ""
    dose: float = np.nan
    time: float = np.nan
    cell_line: str = ""
    n_replicates: int = 0
    method: str = ""
    meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.features = np.asarray(self.features).astype(str)
        self.effect = np.asarray(self.effect, dtype=float)
        n = self.features.size
        if self.effect.shape != (n,):
            raise TinctureDataError(
                f"Signature effect has shape {self.effect.shape}, expected ({n},)."
            )
        if len(set(self.features.tolist())) != n:
            raise TinctureDataError("Signature features must be unique.")
        if self.se is not None:
            self.se = np.asarray(self.se, dtype=float)
            if self.se.shape != (n,):
                raise TinctureDataError("Signature se must have one value per feature.")
            if np.any(self.se[~np.isnan(self.se)] < 0):
                raise TinctureDataError("Signature se must be non-negative.")
        if self.boot is not None:
            self.boot = np.asarray(self.boot, dtype=float)
            if self.boot.ndim != 2 or self.boot.shape[1] != n:
                raise TinctureDataError("Signature boot must have shape (n_boot, n_features).")

    @property
    def n_features(self) -> int:
        return int(self.features.size)

    @property
    def magnitude(self) -> float:
        """Euclidean norm of the effect vector (NaNs ignored)."""
        return float(np.sqrt(np.nansum(self.effect**2)))

    def ci(self, level: float = 0.95) -> tuple[np.ndarray, np.ndarray]:
        """Per-feature CI: bootstrap percentiles if draws exist, else normal from ``se``."""
        if not 0 < level < 1:
            raise ValueError("level must be in (0, 1).")
        a = (1 - level) / 2
        if self.boot is not None and self.boot.shape[0] >= 20:
            lo, hi = np.nanquantile(self.boot, [a, 1 - a], axis=0)
            return lo, hi
        if self.se is None:
            raise TinctureDataError(
                "This signature has no uncertainty (no se and no bootstrap draws).",
                fix="Compute signatures from replicate-level data so tincture can estimate SE.",
            )
        from scipy.stats import norm

        z = norm.ppf(1 - a)
        return self.effect - z * self.se, self.effect + z * self.se

    def align(self, features: np.ndarray | list[str]) -> Signature:
        """Reorder/restrict to ``features``; missing features become NaN."""
        target = np.asarray(features).astype(str)
        if target.shape == self.features.shape and np.array_equal(target, self.features):
            pos = np.arange(target.size)
        else:
            src = pd.Index(self.features.astype(str))
            if not src.is_unique:  # keep the last occurrence, as a dict lookup would
                idx = {f: i for i, f in enumerate(self.features.astype(str))}
                pos = np.array([idx.get(f, -1) for f in target])
            else:
                pos = src.get_indexer(target)
        take = np.where(pos >= 0, pos, 0)

        def pick(a: np.ndarray | None, axis: int = 0) -> np.ndarray | None:
            if a is None:
                return None
            out = np.take(a, take, axis=axis).astype(float)
            if axis == 0:
                out[pos < 0] = np.nan
            else:
                out[:, pos < 0] = np.nan
            return out

        return Signature(
            features=np.asarray(features).astype(str),
            effect=pick(self.effect),  # type: ignore[arg-type]
            se=pick(self.se),
            boot=pick(self.boot, axis=1),
            compound_id=self.compound_id,
            dose=self.dose,
            time=self.time,
            cell_line=self.cell_line,
            n_replicates=self.n_replicates,
            method=self.method,
            meta=dict(self.meta),
        )
