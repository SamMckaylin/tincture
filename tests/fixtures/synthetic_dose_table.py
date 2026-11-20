"""SYNTHETIC data for unit tests only. Never used in case studies (ground rule 2)."""

from __future__ import annotations

import numpy as np
import pandas as pd


def synthetic_dose_table(
    n_compounds: int = 3,
    doses: tuple[float, ...] = (1e-8, 1e-7, 1e-6, 1e-5),
    n_rep: int = 3,
    n_vehicle: int = 6,
    seed: int = 0,
) -> pd.DataFrame:
    """Long table: compound_id, dose (M), replicate, value. Hill-like inhibition + noise."""
    rng = np.random.default_rng(seed)
    rows = [("vehicle", 0.0, r, 1.0 + rng.normal(0, 0.03)) for r in range(n_vehicle)]
    for i in range(n_compounds):
        ec50 = 10 ** (-7 + 0.5 * i)
        for d in doses:
            mu = 1.0 - 0.8 * d / (d + ec50)
            for r in range(n_rep):
                rows.append((f"synthetic_cpd{i + 1}", d, r, mu + rng.normal(0, 0.04)))
    return pd.DataFrame(rows, columns=["compound_id", "dose", "replicate", "value"])
