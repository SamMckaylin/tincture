"""ComboMatrix: drug A × drug B concentration grid of viability (HANDOVER §4.2).

Long table with columns ``conc_a, conc_b`` (molar; ``0`` = drug not added), ``viability``
(fraction of vehicle, 1 = no effect) and ``replicate``. Single-agent wells are rows where the
other drug's concentration is 0. The design is inferred from the grid:

* ``full``: at least 3 non-zero levels of each drug, every pair measured (e.g. 7 × 7);
* ``anchored``: one drug at ≤ 2 levels (e.g. GDSC's 2 × 7 anchored screens);
* ``partial``: anything else (e.g. a ray design), not supported by the Combination Surface.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import pandas as pd

from tincture.errors import TinctureDataError

Design = Literal["full", "anchored", "partial"]
REQUIRED = ("conc_a", "conc_b", "viability")


@dataclass
class ComboMatrix:
    data: pd.DataFrame
    drug_a: str
    drug_b: str
    cell_line: str = ""
    source: str = ""
    meta: dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        miss = [c for c in REQUIRED if c not in self.data.columns]
        if miss:
            raise TinctureDataError(
                f"ComboMatrix missing columns {miss}.",
                fix="Provide conc_a, conc_b (molar, 0 = not added), viability and replicate.",
            )
        d = self.data.copy()
        if "replicate" not in d:
            d["replicate"] = 0
        d["replicate"] = d["replicate"].astype(str)
        for c in ("conc_a", "conc_b", "viability"):
            d[c] = d[c].astype(float)
        if (d["conc_a"] < 0).any() or (d["conc_b"] < 0).any():
            raise TinctureDataError("Concentrations must be ≥ 0 (molar; 0 = drug not added).")
        if not np.isfinite(d[["conc_a", "conc_b"]].to_numpy()).all():
            raise TinctureDataError("Concentrations must be finite.")
        d = d[np.isfinite(d["viability"])]
        if d.empty:
            raise TinctureDataError("ComboMatrix has no finite viability values.")
        self.data = d.reset_index(drop=True)

    # --- structure ------------------------------------------------------------------------
    @property
    def levels_a(self) -> np.ndarray:
        return np.sort(self.data.loc[self.data["conc_a"] > 0, "conc_a"].unique())

    @property
    def levels_b(self) -> np.ndarray:
        return np.sort(self.data.loc[self.data["conc_b"] > 0, "conc_b"].unique())

    @property
    def combo(self) -> pd.DataFrame:
        return self.data[(self.data["conc_a"] > 0) & (self.data["conc_b"] > 0)]

    @property
    def design(self) -> Design:
        ca = np.sort(self.combo["conc_a"].unique())
        cb = np.sort(self.combo["conc_b"].unique())
        if min(len(ca), len(cb)) <= 2 and max(len(ca), len(cb)) >= 3:
            return "anchored"
        pairs = set(map(tuple, self.combo[["conc_a", "conc_b"]].drop_duplicates().to_numpy()))
        if len(ca) >= 3 and len(cb) >= 3 and len(pairs) == len(ca) * len(cb):
            return "full"
        return "partial"

    @property
    def replicates(self) -> list[str]:
        return sorted(self.data["replicate"].unique())

    def single(self, drug: Literal["a", "b"]) -> pd.DataFrame:
        """Single-agent wells of one drug: columns conc, viability, replicate."""
        own, other = ("conc_a", "conc_b") if drug == "a" else ("conc_b", "conc_a")
        s = self.data[(self.data[own] > 0) & (self.data[other] == 0)]
        return s.rename(columns={own: "conc"})[["conc", "viability", "replicate"]]

    def grid(self, data: pd.DataFrame | None = None) -> np.ndarray:
        """Mean combination viability on the (levels_b × levels_a) grid; NaN if unmeasured."""
        d = self.combo if data is None else data[(data["conc_a"] > 0) & (data["conc_b"] > 0)]
        la, lb = self.levels_a, self.levels_b
        m = d.groupby(["conc_b", "conc_a"])["viability"].mean()
        out = np.full((len(lb), len(la)), np.nan)
        ia = {v: i for i, v in enumerate(la)}
        ib = {v: i for i, v in enumerate(lb)}
        for key, v in m.items():
            b, a = key  # type: ignore[misc]
            if a in ia and b in ib:
                out[ib[b], ia[a]] = v
        return out
