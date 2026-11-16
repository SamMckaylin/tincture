"""Loader for the PRISM Repurposing 20Q2 secondary screen (Corsello et al., Nat Cancer 2020).

figshare 10.6084/m9.figshare.20564034.v1 (CC BY 4.0). Verified facts (docs/data_status.md):

* Doses in the files are **µM** and are the actual plated concentrations (they vary
  slightly per compound around an 8-point, 4-fold series from ~10 µM); converted to M.
* Four screens: HTS002 (main), MTS010 (oncology redo, recommended "when available"),
  MTS006, MTS005. ``screen="auto"`` uses, per compound, MTS010 > HTS002 > MTS006 > MTS005.
* Curve parameters follow drc's LL.4: ``upper_limit`` is the zero-dose asymptote (median
  1.0) and ``lower_limit`` the high-dose asymptote, despite the README's swapped wording.
  There is no ``convergence`` column; many fits are pathological, so fits are flagged
  via :meth:`tincture.model.DoseResponse.quality_flags`, never silently trusted.
* Replicate-level log2 fold changes (vs DMSO, ComBat-corrected) give raw points;
  viability = 2 ** LFC.
* Cell-line rows suffixed ``_FAILED_STR`` failed STR profiling and are excluded unless
  ``include_failed_str=True``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from tincture.datasets.cache import fetch
from tincture.datasets.registry import PRISM_SECONDARY
from tincture.errors import TinctureDataError
from tincture.model.dose_response import DoseResponseSet
from tincture.model.measurements import round_dose

SOURCE = "PRISM Repurposing 20Q2 secondary screen, Corsello et al. Nat Cancer 2020 (CC BY 4.0)"
SCREEN_PREFERENCE = ("MTS010", "HTS002", "MTS006", "MTS005")


def compound_core(broad_id: pd.Series) -> pd.Series:
    """'BRD-K12345678-001-02-3' -> 'BRD-K12345678' (compound, not batch)."""
    return broad_id.astype(str).str.slice(0, 13)


def _path(key: str, root: Path | None, download: bool) -> Path:
    return fetch(PRISM_SECONDARY[key], root=root, download=download)


def read_treatments(root: Path | None = None, download: bool = True) -> pd.DataFrame:
    t = pd.read_csv(_path("treatment_info", root, download))
    t["compound_id"] = compound_core(t["broad_id"])
    t["dose_M"] = round_dose(t["dose"] / 1e6)  # µM -> M
    return t


def read_cell_lines(root: Path | None = None, download: bool = True) -> pd.DataFrame:
    c = pd.read_csv(_path("cell_line_info", root, download))
    # Control barcodes (PR300_CONTROL_BARCODE_*) repeat row_name and have no DepMap id.
    c = c[c["depmap_id"].notna()].drop_duplicates("row_name")
    return c.set_index("row_name")


def choose_screens(t: pd.DataFrame, screen: str) -> pd.Series:
    """Map compound_id -> screen_id to use."""
    e = t[t["perturbation_type"] == "experimental_treatment"]
    if screen != "auto":
        if screen not in SCREEN_PREFERENCE:
            raise TinctureDataError(
                f"Unknown screen {screen!r}.", fix=f"Use 'auto' or one of {SCREEN_PREFERENCE}."
            )
        ids = e.loc[e["screen_id"] == screen, "compound_id"].unique()
        return pd.Series(screen, index=pd.Index(ids, name="compound_id"))
    have = e.groupby("compound_id")["screen_id"].agg(set)
    return have.map(lambda s: next(x for x in SCREEN_PREFERENCE if x in s))


def load_prism_secondary(
    compounds: list[str] | None = None,
    cell_lines: list[str] | None = None,
    screen: str = "auto",
    include_failed_str: bool = False,
    root: Path | None = None,
    download: bool = True,
) -> DoseResponseSet:
    """PRISM secondary viability as a :class:`DoseResponseSet`.

    Parameters
    ----------
    compounds
        Compound names (e.g. ``"dabrafenib"``) or Broad compound ids (``"BRD-K…"``).
        ``None`` loads everything (≈1,500 compounds; several GB in memory).
    cell_lines
        DepMap ids (``"ACH-…"``) or CCLE names to keep.
    screen
        ``"auto"`` (per-compound preference MTS010 > HTS002 > MTS006 > MTS005) or one id.
    """
    t = read_treatments(root, download)
    lines = read_cell_lines(root, download)
    chosen = choose_screens(t, screen)
    e = t[t["perturbation_type"] == "experimental_treatment"].copy()
    e = e[e["screen_id"] == e["compound_id"].map(chosen)]
    if compounds is not None:
        names = e["name"].astype(str).str.lower()
        want = {c.lower() for c in compounds}
        sel = names.isin(want) | e["compound_id"].str.lower().isin(want)
        found = set(e.loc[sel, "name"].str.lower()) | set(e.loc[sel, "compound_id"].str.lower())
        missing = sorted(want - found)
        if missing:
            raise TinctureDataError(
                f"Compound(s) not in PRISM secondary: {missing[:5]}.",
                fix="Use names as in the treatment table (lower-case match) or BRD-K ids.",
            )
        e = e[sel]
    # replicate index within compound × dose × screen
    e = e.sort_values(["compound_id", "dose", "detection_plate", "column_name"])
    e["replicate"] = e.groupby(["compound_id", "dose_M"]).cumcount().astype(str)

    lfc_path = _path("logfold_change", root, download)
    header = pd.read_csv(lfc_path, nrows=0).columns
    cols = [c for c in e["column_name"] if c in header]
    missing_cols = len(e) - len(cols)
    lfc = pd.read_csv(lfc_path, usecols=[header[0], *cols], index_col=0)

    keep_rows = lines.index.intersection(lfc.index)
    if not include_failed_str:
        keep_rows = keep_rows[~keep_rows.str.endswith("_FAILED_STR")]
    lfc = lfc.loc[keep_rows]
    info = lines.loc[keep_rows]
    if cell_lines is not None:
        want = set(cell_lines)
        m = info["depmap_id"].isin(want) | info["ccle_name"].isin(want)
        if not m.any():
            raise TinctureDataError(
                f"None of {sorted(want)[:5]} are in PRISM secondary.",
                fix="Use DepMap ids (ACH-…) or CCLE names.",
            )
        lfc, info = lfc.loc[m.to_numpy()], info.loc[m.to_numpy()]

    long = lfc.rename_axis(index="row_name", columns="column_name").stack().rename("lfc")
    long = long.reset_index()
    long["context_id"] = long["row_name"].map(info["depmap_id"])
    meta = e.set_index("column_name")[["compound_id", "name", "dose_M", "replicate", "screen_id"]]
    long = long.join(meta, on="column_name")
    points = pd.DataFrame(
        {
            "compound_id": long["compound_id"],
            "context_id": long["context_id"],
            "dose": long["dose_M"],
            "replicate": long["replicate"],
            "response": np.power(2.0, long["lfc"]),
            "screen_id": long["screen_id"],
            "compound_name": long["name"],
        }
    ).dropna(subset=["response"])

    curves = pd.read_csv(_path("curve_params", root, download), low_memory=False)
    curves["compound_id"] = compound_core(curves["broad_id"])
    curves = curves[curves["compound_id"].isin(points["compound_id"].unique())]
    curves = curves[curves["screen_id"] == curves["compound_id"].map(chosen)]
    curves = curves[curves["depmap_id"].isin(points["context_id"].unique())]
    with np.errstate(divide="ignore", invalid="ignore"):
        log_ec50 = np.log10(curves["ec50"].astype(float) / 1e6)
    fits = pd.DataFrame(
        {
            "compound_id": curves["compound_id"],
            "context_id": curves["depmap_id"],
            "bottom": curves["lower_limit"].astype(float),
            "top": curves["upper_limit"].astype(float),
            "log10_ec50": log_ec50,
            "hill": curves["slope"].astype(float),
            "r2": curves["r2"].astype(float),
            "auc": curves["auc"].astype(float),
            "ic50": curves["ic50"].astype(float) / 1e6,
            "screen_id": curves["screen_id"],
            "broad_id": curves["broad_id"],
        }
    ).drop_duplicates(["compound_id", "context_id"])

    contexts = info.drop_duplicates("depmap_id").set_index("depmap_id")[
        ["ccle_name", "primary_tissue", "secondary_tissue", "passed_str_profiling"]
    ]
    compounds_tbl = (
        e.drop_duplicates("compound_id")
        .set_index("compound_id")[["name", "moa", "target", "smiles", "phase", "screen_id"]]
        .copy()
    )
    # The smiles field repeats the structure once per batch, comma-separated.
    compounds_tbl["smiles"] = compounds_tbl["smiles"].astype(str).str.split(",").str[0].str.strip()
    out = DoseResponseSet(
        points=points.reset_index(drop=True),
        fits=fits.reset_index(drop=True),
        contexts=contexts,
        compounds=compounds_tbl,
        assay="PRISM pooled viability (5 d, Luminex barcodes)",
        fit_source="PRISM drc LL.4 (published)",
        source=SOURCE,
    )
    if missing_cols:
        out.points.attrs["missing_columns"] = missing_cols
    return out
