"""SYNTHETIC count data for unit tests only (ground rule 2)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from tincture.model import PerturbationSet


def synthetic_count_set(
    n_genes: int = 400,
    n_vehicle: int = 8,
    n_rep: int = 3,
    up_genes: int = 20,
    lfc: float = 1.0,
    cell_lines: tuple[str, ...] = ("LINE1",),
    seed: int = 0,
) -> PerturbationSet:
    """Negative-binomial counts. ``synthetic_cpdA`` raises the first ``up_genes`` genes by
    ``lfc`` (log2) at its top dose (half that at the low dose); ``synthetic_cpdB`` does
    nothing; ``synthetic_cpdC`` is a replicate-noisy copy of A (for retrieval tests)."""
    rng = np.random.default_rng(seed)
    base = rng.lognormal(mean=3.0, sigma=1.0, size=n_genes)
    rows, obs = [], []
    for line in cell_lines:
        conds = [("vehicle", 0.0, n_vehicle)]
        conds += [
            (c, d, n_rep)
            for c in ("synthetic_cpdA", "synthetic_cpdB", "synthetic_cpdC")
            for d in (1e-7, 1e-6)
        ]
        for cpd, dose, n in conds:
            mu = base.copy()
            if cpd in ("synthetic_cpdA", "synthetic_cpdC") and dose > 0:
                scale = lfc if dose == 1e-6 else lfc / 2
                mu[:up_genes] *= 2**scale
            for r in range(n):
                lib = rng.uniform(0.8, 1.2)
                disp = 0.05 if cpd != "synthetic_cpdC" else 0.3
                p = 1 / (1 + mu * lib * disp)
                rows.append(rng.negative_binomial(1 / disp, p))
                obs.append(
                    {
                        "compound_id": cpd,
                        "dose": dose,
                        "replicate": f"r{r}",
                        "cell_line": line,
                        "time": 24.0,
                    }
                )
    X = np.vstack(rows).astype(float)
    var = pd.DataFrame(
        {"gene_symbol": [f"G{i}" for i in range(n_genes)]},
        index=[f"ENSG{i:05d}" for i in range(n_genes)],
    )
    o = pd.DataFrame(obs)
    o.index = [f"s{i}" for i in range(len(o))]
    return PerturbationSet(ad.AnnData(X=X, obs=o, var=var), source="synthetic")
