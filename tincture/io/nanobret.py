"""Loader for NanoBRET live-cell kinase occupancy (Vasta et al., Cell Chem Biol 2018; D7).

Facts this loader relies on were verified against the file and the article
(docs/data_status.md, D7):

* Table S1, one sheet; three header rows; **180 kinase rows** (the paper says 178 kinases;
  LCK appears twice, once "Not Determined").
* Columns: kinase, protein ID, Nluc orientation, probe, probe concentration, Z′, probe
  EC50, then **% occupancy, two replicates each, for crizotinib and dasatinib**. Values are
  integers, ``"No Binding"`` or ``"Not Determined"``.
* Per the article, occupancy was measured in live **HEK293** cells expressing full-length
  kinase–NanoLuc fusions at **1 µM crizotinib** and **100 nM dasatinib** (single dose).

**Apparent cellular IC50 (model-based).** A single-dose occupancy ``f`` at concentration
``C`` implies, for 1:1 equilibrium binding with a Hill slope of 1, an apparent IC50 of
``C·(1−f)/f``. This is an inference, not a measured IC50, and is labelled
``NanoBRET apparent IC50 (from single-dose occupancy)``. Replicates give an interval
(min–max over replicates). Occupancies outside [5%, 95%] are clipped for the conversion
(they only bound the IC50), and ``"No Binding"`` is censored (IC50 > C·(1−0.05)/0.05).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from tincture.datasets.cache import fetch
from tincture.datasets.registry import NANOBRET_VASTA2018
from tincture.model.potency import PotencyTable

SOURCE = "Vasta et al., Cell Chem Biol 2018, Table S1 (NanoBRET occupancy, HEK293)"
DOSES_M = {"crizotinib": 1e-6, "dasatinib": 1e-7}
CELL_SYSTEM = "HEK293, full-length kinase–NanoLuc"
F_MIN, F_MAX = 0.05, 0.95


def read_table_s1(root: Path | None = None) -> pd.DataFrame:
    """Table S1 as a long table: kinase, compound, replicate, occupancy (%, NaN = not
    determined), no_binding flag."""
    path = fetch(NANOBRET_VASTA2018, root)
    raw = pd.read_excel(path, header=None).iloc[3:]
    raw.columns = [
        "kinase",
        "protein_id",
        "orientation",
        "probe",
        "probe_conc",
        "zprime",
        "probe_ec50",
        "crizotinib_1",
        "crizotinib_2",
        "dasatinib_1",
        "dasatinib_2",
    ]
    raw = raw[raw["kinase"].notna()]
    long = raw.melt(
        id_vars=["kinase", "protein_id", "orientation", "probe", "probe_conc"],
        value_vars=["crizotinib_1", "crizotinib_2", "dasatinib_1", "dasatinib_2"],
        var_name="col",
        value_name="raw",
    )
    long["compound_id"] = long["col"].str.split("_").str[0]
    long["replicate"] = long["col"].str.split("_").str[1]
    s = long["raw"].astype(str).str.strip()
    long["no_binding"] = s.eq("No Binding")
    long["occupancy"] = pd.to_numeric(long["raw"], errors="coerce")
    long.loc[long["no_binding"], "occupancy"] = 0.0
    long.loc[s.eq("Not Determined"), "occupancy"] = np.nan
    return long.drop(columns=["col"])


def occupancy(root: Path | None = None) -> pd.DataFrame:
    """Mean occupancy (fraction) per compound × kinase over determined replicates (kinases
    listed twice keep the row with determined values)."""
    t = read_table_s1(root).dropna(subset=["occupancy"])
    g = t.groupby(["compound_id", "kinase"])
    out = g["occupancy"].agg(["mean", "min", "max", "size"]).reset_index()
    out["no_binding_all"] = g["no_binding"].all().to_numpy()
    for c in ("mean", "min", "max"):
        out[c] = out[c] / 100.0
    out["dose_m"] = out["compound_id"].map(DOSES_M)
    return out.rename(columns={"size": "n_rep"})


def apparent_ic50(f: np.ndarray, dose: float) -> np.ndarray:
    """C·(1−f)/f with f clipped to [F_MIN, F_MAX] (1:1 binding, Hill 1)."""
    ff = np.clip(np.asarray(f, float), F_MIN, F_MAX)
    return dose * (1 - ff) / ff


def load_nanobret(root: Path | None = None, davis_names: bool = True) -> PotencyTable:
    """Apparent cellular IC50s as a :class:`PotencyTable` (assay type ``IC50_cell``).

    ``davis_names``: rename kinases to Davis 2011 gene symbols (curated synonyms in
    :mod:`tincture.io.kinase_names`) so the two panels join; the original name is kept in
    ``nanobret_name``."""
    from tincture.io.kinase_names import to_davis

    occ = occupancy(root)
    ic = apparent_ic50(occ["mean"].to_numpy(), 1.0) * occ["dose_m"]
    lo = apparent_ic50(occ["max"].to_numpy(), 1.0) * occ["dose_m"]  # higher occupancy → lower IC50
    hi = apparent_ic50(occ["min"].to_numpy(), 1.0) * occ["dose_m"]
    cens = occ["no_binding_all"] | (occ["mean"] <= F_MIN)
    d = pd.DataFrame(
        {
            "compound_id": occ["compound_id"],
            "target": occ["kinase"].map(to_davis) if davis_names else occ["kinase"],
            "nanobret_name": occ["kinase"],
            "px": -np.log10(ic),
            "px_lo": -np.log10(hi),
            "px_hi": -np.log10(lo),
            "assay_type": "IC50_cell",
            "censored": np.where(cens, ">", ""),
            "occupancy": occ["mean"],
            "n_rep": occ["n_rep"],
            "context": CELL_SYSTEM,
            "method": "NanoBRET apparent IC50 (from single-dose occupancy, 1:1 binding)",
        }
    )
    return PotencyTable(d.reset_index(drop=True), source=SOURCE)
