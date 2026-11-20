"""Synthetic potency panel for Selectivity Spectrum tests (not real data)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from tincture.model.potency import PotencyTable


def synthetic_panel(n: int = 60, seed: int = 0, assay: str = "Kd") -> PotencyTable:
    rng = np.random.default_rng(seed)
    rows = [
        {"compound_id": "synthetic_cpd", "target": "SYN_TARGET", "px": 8.5, "assay_type": assay}
    ]
    for k, px in enumerate([8.0, 7.6, 7.0]):  # three strong off-targets
        rows.append(
            {"compound_id": "synthetic_cpd", "target": f"SYN_OFF{k}", "px": px, "assay_type": assay}
        )
    for k in range(n):
        rows.append(
            {
                "compound_id": "synthetic_cpd",
                "target": f"SYN_K{k}",
                "px": float(rng.uniform(5.0, 6.3)),
                "assay_type": assay,
            }
        )
    for k in range(10):  # not detected up to 10 µM
        rows.append(
            {
                "compound_id": "synthetic_cpd",
                "target": f"SYN_ND{k}",
                "px": 5.0,
                "assay_type": assay,
                "censored": ">",
            }
        )
    return PotencyTable(pd.DataFrame(rows), source="synthetic")
