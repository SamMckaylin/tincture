"""PlateTable: well-level values with plate geometry."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from tincture.errors import TinctureDataError
from tincture.model.measurements import WELL_TYPES

#: Standard microplate formats: n_wells -> (rows, cols).
GEOMETRY: dict[int, tuple[int, int]] = {96: (8, 12), 384: (16, 24), 1536: (32, 48)}

REQUIRED = ("plate", "row", "col", "well_type", "value")


def row_label(i: int) -> str:
    """0 -> 'A', 25 -> 'Z', 26 -> 'AA', 31 -> 'AF' (1536-well convention)."""
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def position_to_rowcol(
    position: np.ndarray, n_wells: int, one_based: bool = True, order: str = "row-major"
) -> tuple[np.ndarray, np.ndarray]:
    """Convert linear well positions to 0-based (row, col)."""
    if n_wells not in GEOMETRY:
        raise TinctureDataError(
            f"Unknown plate format {n_wells}.", fix=f"Use one of {list(GEOMETRY)}."
        )
    nr, nc = GEOMETRY[n_wells]
    p = np.asarray(position, dtype=int) - (1 if one_based else 0)
    if np.any(p < 0) or np.any(p >= n_wells):
        raise TinctureDataError(f"Well positions outside 1..{n_wells}.")
    if order == "row-major":
        return p // nc, p % nc
    if order == "col-major":
        return p % nr, p // nr
    raise ValueError("order must be 'row-major' or 'col-major'.")


@dataclass
class PlateTable:
    """Well-level measurements on one or more plates of a single format.

    ``data`` columns: ``plate``, ``row``, ``col`` (0-based), ``well_type`` (one of
    :data:`tincture.model.measurements.WELL_TYPES`), ``value``; optional ``compound_id``,
    ``dose`` (molar), ``cell_line``, ``tag`` (source-specific label), ``plate_order``.

    Large loaders (GDSC: 13.5 M wells) store text columns as pandas categoricals to save
    memory; use :meth:`relative_to` rather than ``.map`` on them for per-plate maths.
    """

    data: pd.DataFrame
    n_wells: int
    source: str = ""
    value_label: str = "intensity"

    def __post_init__(self) -> None:
        miss = [c for c in REQUIRED if c not in self.data.columns]
        if miss:
            raise TinctureDataError(
                f"PlateTable is missing columns {miss}.",
                fix="Provide plate, row, col (0-based), well_type and value.",
            )
        if self.n_wells not in GEOMETRY:
            raise TinctureDataError(
                f"Unsupported plate format {self.n_wells}.", fix=f"Use one of {list(GEOMETRY)}."
            )
        nr, nc = GEOMETRY[self.n_wells]
        r, c = self.data["row"].to_numpy(), self.data["col"].to_numpy()
        if (r < 0).any() or (r >= nr).any() or (c < 0).any() or (c >= nc).any():
            raise TinctureDataError(f"Row/col outside a {nr}×{nc} ({self.n_wells}-well) plate.")
        unknown = set(self.data["well_type"].unique()) - WELL_TYPES
        if unknown:
            raise TinctureDataError(
                f"Unknown well_type value(s) {sorted(unknown)}.", fix=f"Use {sorted(WELL_TYPES)}."
            )
        dup = self.data.duplicated(["plate", "row", "col"])
        if dup.any():
            raise TinctureDataError(
                f"{int(dup.sum())} duplicated plate/row/col entries.",
                fix="Each well may appear once per plate.",
            )

    @property
    def shape(self) -> tuple[int, int]:
        return GEOMETRY[self.n_wells]

    @property
    def plates(self) -> list[str]:
        if "plate_order" in self.data:
            o = self.data.drop_duplicates("plate").sort_values("plate_order")
            return o["plate"].astype(str).tolist()
        return sorted(self.data["plate"].astype(str).unique())

    def is_edge(self) -> np.ndarray:
        """Boolean mask of outermost-ring wells."""
        nr, nc = self.shape
        r, c = self.data["row"].to_numpy(), self.data["col"].to_numpy()
        return (r == 0) | (r == nr - 1) | (c == 0) | (c == nc - 1)

    def plate_matrix(self, plate: str, well_type: str | None = None) -> np.ndarray:
        """Values of one plate as a (rows × cols) array; NaN where absent/filtered."""
        nr, nc = self.shape
        d = self.data[self.data["plate"].astype(str) == str(plate)]
        if d.empty:
            raise TinctureDataError(f"No plate {plate!r}.", fix="See PlateTable.plates.")
        if well_type is not None:
            d = d[d["well_type"] == well_type]
        m = np.full((nr, nc), np.nan)
        m[d["row"].to_numpy(), d["col"].to_numpy()] = d["value"].to_numpy(dtype=float)
        return m

    def relative_to(self, mask: np.ndarray | pd.Series) -> pd.Series:
        """Each well's value divided by its plate's median over the wells in ``mask``
        (e.g. negative controls). Works with categorical plate columns."""
        m = np.asarray(mask, dtype=bool)
        plate = self.data["plate"].astype(str)
        ref = self.data.loc[m, "value"].groupby(plate[m]).median()
        return self.data["value"].astype(float) / plate.map(ref).astype(float)

    def well_name(self) -> pd.Series:
        return self.data["row"].map(row_label) + (self.data["col"] + 1).astype(str)

    def subset_plates(self, plates: list[str]) -> PlateTable:
        d = self.data[self.data["plate"].astype(str).isin([str(p) for p in plates])]
        return PlateTable(d.reset_index(drop=True), self.n_wells, self.source, self.value_label)
