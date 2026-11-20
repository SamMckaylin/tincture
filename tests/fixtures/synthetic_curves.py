"""Synthetic dose-response curves for Window / 4PL tests (not real data)."""

from __future__ import annotations

import numpy as np

from tincture.model.dose_response import DoseResponse, FourPL

DOSES = 10 ** np.linspace(-9, -5, 8)


def synthetic_curve(
    fit: FourPL,
    noise: float = 0.03,
    reps: int = 3,
    seed: int = 0,
    name: str = "synthetic_cpd",
    context: dict | None = None,
) -> DoseResponse:
    rng = np.random.default_rng(seed)
    d = np.repeat(DOSES, reps)
    y = fit(d) + rng.normal(0, noise, d.size)
    return DoseResponse(
        compound_id=name,
        dose=d,
        response=y,
        replicate=np.tile(np.arange(reps).astype(str), DOSES.size),
        assay="synthetic viability",
        context=context or {"cell_line": "SYNTH-1", "readout": "synthetic 3 d"},
    )


def synthetic_units(fit: FourPL, n_units: int = 40, spread: float = 0.15, seed: int = 0):
    """Many units (e.g. cell lines) around one curve; replicate = unit id."""
    rng = np.random.default_rng(seed)
    d, y, u = [], [], []
    for k in range(n_units):
        shift = rng.normal(0, spread)
        f = FourPL(fit.bottom, fit.top, fit.log10_ec50 + shift * 3, fit.hill)
        for x in DOSES:
            for _ in range(2):
                d.append(x)
                y.append(float(f(x)) + rng.normal(0, 0.05))
                u.append(f"unit{k}")
    return DoseResponse(
        compound_id="synthetic_cpd",
        dose=np.array(d),
        response=np.array(y),
        replicate=np.array(u),
        context={"cell_line": f"median of {n_units} synthetic lines", "summary": "median"},
    )
