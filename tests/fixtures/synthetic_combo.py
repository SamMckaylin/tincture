"""Synthetic combination matrices for unit tests only (never real data)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from tincture.model.combo import ComboMatrix

LEVELS = 1e-6 * 10.0 ** (-np.arange(7)[::-1] * 0.5)  # 7 levels, 0.5-log steps, top 1 µM


def hill_v(x: np.ndarray, ec50: float, h: float = 1.2, emax: float = 0.9) -> np.ndarray:
    """Viability of a single agent (1 at zero dose, 1 − emax at infinite dose)."""
    x = np.asarray(x, float)
    return 1 - emax * x**h / (x**h + ec50**h)


def synthetic_matrix(
    synergy: float = 0.0,
    n_rep: int = 3,
    sd: float = 0.02,
    ec_a: float = 1e-7,
    ec_b: float = 3e-7,
    region: tuple[int, int] = (3, 6),
    seed: int = 0,
    levels_b: np.ndarray | None = None,
) -> ComboMatrix:
    """Bliss-independent 7×7 surface; ``synergy`` subtracts extra viability in the
    sub-grid ``region`` (index range on both axes, e.g. mid-high doses)."""
    rng = np.random.default_rng(seed)
    la = LEVELS
    lb = LEVELS if levels_b is None else levels_b
    rows = []
    for r in range(n_rep):
        for a in np.r_[0.0, la]:
            for b in np.r_[0.0, lb]:
                if a == 0 and b == 0:
                    continue
                v = float(hill_v(a, ec_a) * hill_v(b, ec_b))
                ia = int(np.searchsorted(la, a)) if a > 0 else -1
                ib = int(np.searchsorted(lb, b)) if b > 0 else -1
                if region[0] <= ia < region[1] and region[0] <= ib < region[1]:
                    v -= synergy
                rows.append(
                    {"conc_a": a, "conc_b": b, "viability": v + rng.normal(0, sd), "replicate": r}
                )
    return ComboMatrix(pd.DataFrame(rows), "synthetic A", "synthetic B", "synthetic line")


def anchored_matrix() -> ComboMatrix:
    """2 × 7 anchored design (drug B at two anchor concentrations)."""
    return synthetic_matrix(levels_b=LEVELS[[3, 5]], n_rep=2)
