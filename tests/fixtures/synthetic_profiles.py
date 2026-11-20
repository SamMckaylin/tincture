"""SYNTHETIC replicate-level profiles for retrievability tests only (ground rule 2)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from tincture.model import PerturbationSet


def synthetic_profiles(
    n_controls: int = 16,
    n_rep: int = 3,
    n_features: int = 50,
    effects: dict[str, float] | None = None,
    single: bool = True,
    seed: int = 0,
) -> PerturbationSet:
    """Controls (dose 0) + compounds with a fixed direction scaled by ``effects``; plus one
    single-replicate compound when ``single``."""
    rng = np.random.default_rng(seed)
    effects = effects or {"synthetic_strong": 4.0, "synthetic_medium": 1.5, "synthetic_none": 0.0}
    rows, obs = [], []
    for i in range(n_controls):
        rows.append(rng.normal(0, 1, n_features))
        obs.append(("vehicle", 0.0, f"c{i}"))
    for name, size in effects.items():
        d = rng.normal(0, 1, n_features)
        for r in range(n_rep):
            rows.append(size * d + rng.normal(0, 1, n_features))
            obs.append((name, 1e-6, f"r{r}"))
    if single:
        rows.append(rng.normal(0, 1, n_features) + 3)
        obs.append(("synthetic_single", 1e-6, "r0"))
    o = pd.DataFrame(obs, columns=["compound_id", "dose", "replicate"])
    o["compound_name"] = o["compound_id"]
    o.index = [f"p{i}" for i in range(len(o))]
    X = np.vstack(rows)
    return PerturbationSet(
        ad.AnnData(X=X, obs=o, var=pd.DataFrame(index=[f"f{j}" for j in range(n_features)]))
    )
