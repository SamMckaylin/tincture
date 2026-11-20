"""The canonical long measurement table (HANDOVER §4.1) and its validator."""

from __future__ import annotations

import numpy as np
import pandas as pd

from tincture.errors import TinctureDataError

#: Required columns of the long table.
REQUIRED: tuple[str, ...] = ("compound_id", "dose", "replicate", "feature", "value")
#: Optional columns, with their dtype family.
OPTIONAL: dict[str, str] = {
    "compound_name": "str",
    "smiles": "str",
    "time": "float",
    "plate": "str",
    "well": "str",
    "well_type": "str",
    "cell_line": "str",
    "se": "float",
}
#: Allowed ``well_type`` values. ``blank`` (media only, no cells) extends the handover's
#: list because GDSC plates use blanks as their full-kill reference (docs/data_status.md).
WELL_TYPES: frozenset[str] = frozenset({"compound", "neg_ctrl", "pos_ctrl", "blank", "empty"})


def validate_long(df: pd.DataFrame, require_vehicle: bool = False) -> pd.DataFrame:
    """Validate and normalise a long measurement table.

    Returns a copy with canonical dtypes (``dose``/``value``/``time``/``se`` float, ids as
    strings). Rows whose ``value`` is NaN are kept (missing is information); rows with a
    missing ``dose`` are rejected because dose is a first-class axis.

    Raises
    ------
    TinctureDataError
        With a fix suggestion, for missing columns, negative or non-numeric doses,
        unknown ``well_type`` values, or (if ``require_vehicle``) no vehicle rows.
    """
    if not isinstance(df, pd.DataFrame):
        raise TinctureDataError(
            f"Expected a pandas DataFrame, got {type(df).__name__}.",
            fix="Convert your data to a long DataFrame (one row per measurement).",
        )
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise TinctureDataError(
            f"Long table is missing required column(s): {missing}.",
            fix=f"Provide columns {list(REQUIRED)} (dose in molar, 0 for vehicle). "
            "If you have no dose, use tc.model.add_constant_dose(df) to state that explicitly.",
        )
    out = df.copy()
    for c in ("compound_id", "feature", "replicate"):
        out[c] = out[c].astype(str)
    out["dose"] = pd.to_numeric(out["dose"], errors="coerce")
    bad_dose = out["dose"].isna() & df["dose"].notna()
    if bad_dose.any():
        raise TinctureDataError(
            f"{int(bad_dose.sum())} dose value(s) are not numeric, e.g. "
            f"{df.loc[bad_dose, 'dose'].iloc[0]!r}.",
            fix="Convert doses to molar floats (e.g. 1e-6 for 1 µM).",
        )
    if out["dose"].isna().any():
        raise TinctureDataError(
            f"{int(out['dose'].isna().sum())} row(s) have no dose.",
            fix="Fill dose (molar; 0 for vehicle) or drop those rows explicitly.",
        )
    if (out["dose"] < 0).any():
        raise TinctureDataError("Negative doses found.", fix="Doses must be molar and ≥ 0.")
    if (out["dose"] > 1.0).any():
        raise TinctureDataError(
            f"Doses up to {out['dose'].max():g} M found; this looks like µM or nM, not molar.",
            fix="Convert to molar (µM × 1e-6, nM × 1e-9).",
        )
    out["value"] = pd.to_numeric(out["value"], errors="coerce").astype(float)
    for c, kind in OPTIONAL.items():
        if c not in out.columns:
            continue
        if kind == "float":
            out[c] = pd.to_numeric(out[c], errors="coerce").astype(float)
        else:
            out[c] = out[c].astype("string")
    if "well_type" in out.columns:
        unknown = set(out["well_type"].dropna().unique()) - WELL_TYPES
        if unknown:
            raise TinctureDataError(
                f"Unknown well_type value(s): {sorted(unknown)}.",
                fix=f"Use one of {sorted(WELL_TYPES)}.",
            )
    if require_vehicle and not (out["dose"] == 0).any():
        raise TinctureDataError(
            "No vehicle rows (dose == 0).",
            fix="Include vehicle measurements with dose 0; every tincture plot needs a reference.",
        )
    if np.isinf(out["value"]).any():
        raise TinctureDataError(
            "Infinite values found.", fix="Replace ±inf with NaN or finite values."
        )
    return out


def round_dose(dose: np.ndarray | pd.Series, sig: int = 6) -> np.ndarray:
    """Round molar doses to ``sig`` significant figures.

    Unit conversion (e.g. µM / 1e6) leaves binary noise (9.769999999999999e-09), which
    breaks grouping by dose; rounding restores the intended values (9.77e-09).
    """
    d = np.asarray(dose, dtype=float)
    uniq = np.unique(d[np.isfinite(d)])
    lookup = {v: float(f"{v:.{sig}g}") for v in uniq}
    return np.array([lookup.get(v, v) for v in d.tolist()], dtype=float) if d.size else d


def add_constant_dose(df: pd.DataFrame, dose: float = np.nan) -> pd.DataFrame:
    """Explicitly add a ``dose`` column for data without dose information.

    With the default (NaN) the table will still fail validation: this is deliberate so
    that dose-free data can only enter dose-aware plots after the user states a dose.
    """
    out = df.copy()
    out["dose"] = dose
    return out
