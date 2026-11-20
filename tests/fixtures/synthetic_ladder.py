"""Synthetic Engagement Ladder tables (not real data)."""

from __future__ import annotations

import pandas as pd

from tincture.model.potency import PotencyTable


def ladder(rows: list[dict]) -> PotencyTable:
    base = {"assay_type": "IC50", "censored": "", "cell_system": ""}
    return PotencyTable(pd.DataFrame([{**base, **r} for r in rows]), source="synthetic")


def typical(
    pheno_px: float = 7.0,
    cell_px: float = 7.5,
    pheno_cens: str = "",
    cell_sys: str = "SYN-HEK",
    pheno_sys: str = "SYN-LINE",
) -> PotencyTable:
    return ladder(
        [
            {
                "compound_id": "synthetic_a",
                "target": "SYNK1",
                "rung": "biochemical",
                "px": 9.0,
                "assay_type": "Kd",
                "cell_system": "cell-free",
            },
            {
                "compound_id": "synthetic_a",
                "target": "SYNK1",
                "rung": "cellular",
                "px": cell_px,
                "px_lo": cell_px - 0.2,
                "px_hi": cell_px + 0.2,
                "assay_type": "IC50_cell",
                "cell_system": cell_sys,
            },
            {
                "compound_id": "synthetic_a",
                "target": "SYNK1",
                "rung": "phenotype",
                "px": pheno_px,
                "censored": pheno_cens,
                "cell_system": pheno_sys,
            },
        ]
    )
