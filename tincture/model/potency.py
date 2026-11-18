"""PotencyTable: compound × target (or rung) potency with assay type and context.

HANDOVER §4.2: compound × target × value (pIC50/pKd/pEC50) × assay type × assay context
(cell line, ATP concentration, readout). Values are stored as **pX = −log10(molar)**;
censored measurements (e.g. "Kd > 10 µM", not detected at the top concentration) keep the
limit as the value with ``censored = ">"`` (true potency weaker than the limit, e.g. not
detected up to the top dose) or ``"<"`` (more potent than the limit, e.g. already fully
active at the lowest dose), so they are never drawn as measurements.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from tincture.errors import TinctureDataError

ASSAY_TYPES = ("Kd", "Ki", "IC50", "EC50", "Kd_app", "IC50_cell")
REQUIRED = ("compound_id", "target", "px", "assay_type")


def px_from_molar(m: np.ndarray | float) -> np.ndarray:
    m = np.asarray(m, float)
    if np.any(m <= 0):
        raise TinctureDataError("Potencies must be positive concentrations (molar).")
    return -np.log10(m)


@dataclass
class PotencyTable:
    """Long table with columns ``compound_id, target, px, assay_type`` (+ optional
    ``censored`` ('', '>' = weaker than the limit, '<' = more potent than the limit),
    ``context`` and any metadata)."""

    data: pd.DataFrame
    source: str = ""

    def __post_init__(self) -> None:
        miss = [c for c in REQUIRED if c not in self.data.columns]
        if miss:
            raise TinctureDataError(f"PotencyTable missing columns {miss}.")
        d = self.data.copy()
        if "censored" not in d:
            d["censored"] = ""
        d["censored"] = d["censored"].fillna("").astype(str)
        bad = ~d["censored"].isin(["", ">", "<"])
        if bad.any():
            raise TinctureDataError(
                "`censored` must be '', '>' (weaker than the limit) or '<' (more potent)."
            )
        unknown = sorted(set(d["assay_type"]) - set(ASSAY_TYPES))
        if unknown:
            raise TinctureDataError(
                f"Unknown assay types {unknown}.", fix=f"Use one of {list(ASSAY_TYPES)}."
            )
        if not np.all(np.isfinite(d["px"].astype(float))):
            raise TinctureDataError("pX values must be finite (use `censored` for limits).")
        d["px"] = d["px"].astype(float)
        self.data = d.reset_index(drop=True)

    @classmethod
    def from_molar(
        cls,
        df: pd.DataFrame,
        value_col: str = "value_molar",
        source: str = "",
    ) -> PotencyTable:
        d = df.copy()
        d["px"] = px_from_molar(d[value_col].to_numpy())
        return cls(d, source=source)

    def compound(self, compound_id: str) -> pd.DataFrame:
        d = self.data[self.data["compound_id"] == compound_id]
        if d.empty:
            raise TinctureDataError(f"No potencies for {compound_id!r}.")
        return d

    @property
    def assay_types(self) -> list[str]:
        return sorted(self.data["assay_type"].unique())
