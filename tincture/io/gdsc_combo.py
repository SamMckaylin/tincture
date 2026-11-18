"""GDSC combination screens: full 7×7 matrices (D5, Bashi et al. 2024) and anchored 2×7
screens (D4, Jaaks et al. 2022; D6, Vis et al. 2024). File facts in docs/data_status.md.

* :func:`bashi_results` — per combination × cell line summary (GDSC's own Bliss/HSA metrics).
* :func:`load_bashi_matrices` — :class:`ComboMatrix` objects from the per-well file:
  combination wells (``inhibition`` → viability = 1 − inhibition) plus single-agent values
  taken from the matched monotherapy responses (``lib*_y``, measured) at each concentration.
  One replicate = one ``BARCODE`` × ``cmatrix``. ``meta["noise_sd"]`` is the root mean
  square of the two single-agent fit RMSEs (GDSC's own fits), offered as the well-level noise
  for single-replicate matrices.
* :func:`load_jaaks2022` / :func:`load_vis2024` — the anchored fitted tables, and
  :func:`anchored_lift_table` to turn them into Combination Lift input.

Concentrations in the GDSC files are µM; tincture uses molar.
"""

from __future__ import annotations

import io
import zipfile
from collections.abc import Sequence

import numpy as np
import pandas as pd

from tincture.datasets.cache import fetch
from tincture.datasets.registry import BASHI2024, CMP_MODEL_LIST, JAAKS2022, VIS2024_SUPPLEMENT
from tincture.errors import TinctureDataError
from tincture.model.combo import ComboMatrix

UM = 1e-6
_SCREENS = ("gdsc007", "gdsc008", "gdsc009", "gdsc010", "sandpiper01")
_WELL_FILES = {
    "gdsc-007": "combo_well_stats_GDSC_007_05May22.csv",
    "gdsc-008": "combo_well_stats_GDSC_008_12May22.csv",
    "gdsc-009": "combo_well_stats_GDSC_009_12May22.csv",
    "gdsc-010": "combo_well_stats_GDSC_010_11May22.csv",
    "sandpiper-01": "combo_well_stats_Sandpiper-01_04May22.csv",
}
_WELL_COLS = [
    "BARCODE", "CELL_LINE_NAME", "TISSUE", "CANCER_TYPE", "cmatrix", "lib1_drug_id",
    "lib1_conc", "lib1_y", "lib1_RMSE", "lib2_drug_id", "lib2_conc", "lib2_y", "lib2_RMSE",
    "inhibition",
]  # fmt: skip


def bashi_results() -> pd.DataFrame:
    """All per-screen matrix summaries, with a ``screen`` column (e.g. 'gdsc-008')."""
    out = []
    for k in _SCREENS:
        p = fetch(BASHI2024[k])
        out.append(pd.read_csv(p, low_memory=False).assign(screen=p.name.split("_matrix")[0]))
    return pd.concat(out, ignore_index=True)


def _pair(res: pd.DataFrame, drug_1: str, drug_2: str) -> pd.DataFrame:
    m = res[(res["lib1_name"] == drug_1) & (res["lib2_name"] == drug_2)]
    if m.empty:
        m = res[(res["lib1_name"] == drug_2) & (res["lib2_name"] == drug_1)]
    if m.empty:
        raise TinctureDataError(
            f"{drug_1} + {drug_2} is not in the GDSC matrix screens.",
            fix="Check names against bashi_results()['lib1_name'/'lib2_name'].",
        )
    return m


def load_bashi_matrices(
    drug_1: str, drug_2: str, cell_lines: Sequence[str] | None = None
) -> dict[str, ComboMatrix]:
    """7×7 matrices for one drug pair (orientation as screened: drug A = lib1)."""
    res = _pair(bashi_results(), drug_1, drug_2)
    if cell_lines is not None:
        res = res[res["CELL_LINE_NAME"].isin(cell_lines)]
        missing = sorted(set(cell_lines) - set(res["CELL_LINE_NAME"]))
        if missing:
            raise TinctureDataError(f"Not screened with {drug_1} + {drug_2}: {missing}.")
    id1, id2 = int(res["lib1_ID"].iloc[0]), int(res["lib2_ID"].iloc[0])
    name_a, name_b = str(res["lib1_name"].iloc[0]), str(res["lib2_name"].iloc[0])
    lines = set(res["CELL_LINE_NAME"])
    zf = zipfile.ZipFile(fetch(BASHI2024["well_stats"]))
    parts = []
    for screen in sorted(res["screen"].unique()):
        with zf.open(_WELL_FILES[screen]) as fh:
            for ch in pd.read_csv(fh, usecols=_WELL_COLS, chunksize=500_000):
                sel = ch[
                    (ch["lib1_drug_id"] == id1)
                    & (ch["lib2_drug_id"] == id2)
                    & ch["CELL_LINE_NAME"].isin(lines)
                ]
                if len(sel):
                    parts.append(sel)
    if not parts:
        raise TinctureDataError(f"No wells found for {drug_1} + {drug_2}.")
    w = pd.concat(parts, ignore_index=True)
    out: dict[str, ComboMatrix] = {}
    for line, g in w.groupby("CELL_LINE_NAME"):
        g = g.assign(replicate=g["BARCODE"].astype(str) + "_" + g["cmatrix"].astype(str))
        combo = pd.DataFrame(
            {
                "conc_a": g["lib1_conc"] * UM,
                "conc_b": g["lib2_conc"] * UM,
                "viability": 1 - g["inhibition"],
                "replicate": g["replicate"],
            }
        )
        sa = g.drop_duplicates(["replicate", "lib1_conc"])
        sb = g.drop_duplicates(["replicate", "lib2_conc"])
        singles = pd.concat(
            [
                pd.DataFrame(
                    {"conc_a": sa["lib1_conc"] * UM, "conc_b": 0.0,
                     "viability": 1 - sa["lib1_y"], "replicate": sa["replicate"]}
                ),
                pd.DataFrame(
                    {"conc_a": 0.0, "conc_b": sb["lib2_conc"] * UM,
                     "viability": 1 - sb["lib2_y"], "replicate": sb["replicate"]}
                ),
            ]
        )  # fmt: skip
        rmse = np.sqrt(np.nanmean(np.r_[g["lib1_RMSE"].to_numpy(), g["lib2_RMSE"].to_numpy()] ** 2))
        summary = res[res["CELL_LINE_NAME"] == line]
        out[str(line)] = ComboMatrix(
            pd.concat([combo, singles], ignore_index=True),
            name_a,
            name_b,
            cell_line=str(line),
            source="Bashi et al. 2024 (GDSC matrix screens, figshare project 163378)",
            meta={
                "noise_sd": float(rmse),
                "cancer_type": str(g["CANCER_TYPE"].iloc[0]),
                "gdsc_summary": summary,
            },
        )
    return out


# ---------------------------------------------------------------------------------------
# anchored screens
# ---------------------------------------------------------------------------------------
def load_jaaks2022(screen: str = "original") -> pd.DataFrame:
    """Jaaks et al. 2022 fitted anchored data (one row per anchor conc × library × line ×
    replicate plate)."""
    key = {"original": "original_fitted", "validation": "validation_fitted"}.get(screen)
    if key is None:
        raise TinctureDataError("screen must be 'original' or 'validation'.")
    return pd.read_csv(fetch(JAAKS2022[key]), low_memory=False)


def load_vis2024() -> pd.DataFrame:
    """Vis et al. 2024 pan-cancer anchored fitted data (Data S1)."""
    z = zipfile.ZipFile(fetch(VIS2024_SUPPLEMENT))
    inner = zipfile.ZipFile(io.BytesIO(z.read("mmc7.zip")))
    with inner.open("data-code/SupplementalMaterials/pancan_combi_fitted_data_jan21.tsv.gz") as fh:
        return pd.read_csv(fh, sep="\t", compression="gzip", low_memory=False)


def cell_model_annotations() -> pd.DataFrame:
    """Cell Model Passports model list (model_id = SIDM, msi_status, tissue …)."""
    return pd.read_csv(fetch(CMP_MODEL_LIST), low_memory=False)


def anchored_lift_table(
    d: pd.DataFrame, anchor: str, library: str, max_rmse: float | None = None
) -> pd.DataFrame:
    """Lift input from a GDSC anchored table for one anchor × library pair: viabilities at
    the highest library concentration (anchor alone, library Emax, combination Emax) per
    cell line × anchor concentration × replicate. ``fit_ok`` is False when any value is
    missing or (if given) a curve RMSE exceeds ``max_rmse``."""
    m = d[(d["ANCHOR_NAME"] == anchor) & (d["LIBRARY_NAME"] == library)]
    if m.empty:
        raise TinctureDataError(f"No rows for anchor {anchor} × library {library}.")
    ok = np.isfinite(m[["ANCHOR_VIABILITY", "LIBRARY_EMAX", "SYNERGY_OBS_EMAX"]]).all(axis=1)
    if max_rmse is not None:
        ok &= (m["LIBRARY_RMSE"] <= max_rmse) & (m["SYNERGY_RMSE"] <= max_rmse)
    out = pd.DataFrame(
        {
            "cell_line": m["CELL_LINE_NAME"].astype(str),
            "sidm": m.get("SIDM", pd.NA),
            "anchor_conc_um": pd.to_numeric(m["ANCHOR_CONC"], errors="coerce"),
            "single_a": m["ANCHOR_VIABILITY"],
            "single_b": m["LIBRARY_EMAX"],
            "combo": m["SYNERGY_OBS_EMAX"],
            "bliss_expected": m["SYNERGY_EXP_EMAX"],
            "fit_ok": ok.to_numpy(),
        }
    )
    out["anchor"] = [f"{anchor} {c:g} µM" for c in out["anchor_conc_um"]]
    return out.reset_index(drop=True)
