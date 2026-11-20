"""SYNTHETIC long tables for unit tests only (ground rule 2)."""

from __future__ import annotations

import numpy as np
import pandas as pd


def synthetic_long_table(seed: int = 0) -> pd.DataFrame:
    """2 compounds × 3 doses × 3 replicates + 3 vehicle samples, 5 features."""
    rng = np.random.default_rng(seed)
    rows = []
    samples = [("vehicle", 0.0, f"r{i}") for i in range(1, 4)]
    samples += [
        (c, d, f"r{i}")
        for c in ("synthetic_cpdA", "synthetic_cpdB")
        for d in (1e-7, 1e-6, 1e-5)
        for i in range(1, 4)
    ]
    for cpd, dose, rep in samples:
        for j in range(5):
            rows.append((cpd, dose, rep, f"g{j}", float(rng.normal(j, 1)), "A549"))
    return pd.DataFrame(
        rows, columns=["compound_id", "dose", "replicate", "feature", "value", "cell_line"]
    )
