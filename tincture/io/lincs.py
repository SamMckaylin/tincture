"""Loader for LINCS L1000 phase II, Level 5 signatures (GEO GSE70138).

Facts this loader relies on were verified against the files (docs/data_status.md, D9):

* Level 5 ``COMPZ_n118050x12328`` is GCTX (HDF5): ``0/DATA/0/matrix`` is float32 with one
  **row per signature** (118,050) and one column per gene (12,328); ``0/META/COL/id`` holds
  ``sig_id`` and ``0/META/ROW/id`` the Entrez ``pr_gene_id``. Values are moderated z-scores
  (MODZ) of replicate profiles vs plate population.
* Only the 978 **landmark** genes (``gene_info.pr_is_lm == 1``) are measured; the rest are
  inferred and are not used.
* ``sig_info``: ``pert_id`` (Broad ID, ``BRD-K12345678…``), ``pert_iname``, ``pert_type``
  (``trt_cp`` = compound), ``cell_id``, ``pert_idose`` (e.g. ``"10.0 um"``), ``pert_itime``
  (``"24 h"``). ``sig_metrics.tas`` = transcriptional activity score.
* ``pert_info`` carries ``canonical_smiles`` and ``inchi_key``.

The 5 GB Level 5 file is downloaded once (SHA-512 checked), unpacked to a temporary file,
reduced to landmark genes in a ~0.5 GB cache (``data/processed/lincs_gse70138_landmark.h5ad``)
and the temporary file is deleted.
"""

from __future__ import annotations

import gzip
import re
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from tincture.datasets.cache import data_home, fetch
from tincture.datasets.registry import LINCS_GSE70138
from tincture.errors import TinctureDataError

SOURCE = "LINCS L1000 phase II, Level 5 (GEO GSE70138; Subramanian et al., Cell 2017)"
_UNIT = {"um": 1e-6, "nm": 1e-9, "mm": 1e-3, "pm": 1e-12}


def parse_dose(s: object) -> float:
    """``"10.0 um"`` → 1e-5 (molar); unknown or control doses (``-666``) → NaN."""
    m = re.match(r"^\s*([0-9.]+)\s*([a-zA-Zµ]+)\s*$", str(s))
    if not m:
        return float("nan")
    unit = m.group(2).lower().replace("µ", "u")
    return float(m.group(1)) * _UNIT[unit] if unit in _UNIT else float("nan")


def parse_time(s: object) -> float:
    m = re.match(r"^\s*([0-9.]+)\s*h\s*$", str(s))
    return float(m.group(1)) if m else float("nan")


def _meta(root: Path | None, download: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    rd = {k: fetch(LINCS_GSE70138[k], root, download) for k in LINCS_GSE70138 if k != "level5"}
    si = pd.read_csv(rd["sig_info"], sep="\t", low_memory=False)
    sm = pd.read_csv(rd["sig_metrics"], sep="\t", low_memory=False)
    pi = pd.read_csv(rd["pert_info"], sep="\t")
    gi = pd.read_csv(rd["gene_info"], sep="\t")
    obs = si.merge(
        sm[["sig_id", "tas", "distil_cc_q75", "distil_nsample"]], on="sig_id", how="left"
    ).merge(pi[["pert_id", "canonical_smiles", "inchi_key"]], on="pert_id", how="left")
    obs["dose"] = obs["pert_idose"].map(parse_dose)
    obs["time"] = obs["pert_itime"].map(parse_time)
    obs["pert_core"] = obs["pert_id"].astype(str).str[:13]
    lm = gi[gi["pr_is_lm"] == 1].copy()
    return obs, lm


def build_landmark_cache(
    root: Path | None = None, out: Path | None = None, download: bool = True, block: int = 4000
) -> Path:
    """Download Level 5, keep the 978 landmark genes, write the h5ad cache."""
    import h5py

    out = out or (data_home().parent / "processed" / "lincs_gse70138_landmark.h5ad")
    obs, lm = _meta(root, download)
    gz = fetch(LINCS_GSE70138["level5"], root, download)
    tmp = out.parent / "_lincs_level5_tmp.gctx"
    out.parent.mkdir(parents=True, exist_ok=True)
    if not tmp.exists():
        with gzip.open(gz, "rb") as src, open(tmp, "wb") as dst:
            shutil.copyfileobj(src, dst, length=1 << 24)
    with h5py.File(tmp, "r") as f:
        rows = f["0/META/ROW/id"][:].astype(str)
        cols = f["0/META/COL/id"][:].astype(str)
        pos = {g: i for i, g in enumerate(rows)}
        want = [str(g) for g in lm["pr_gene_id"]]
        missing = [g for g in want if g not in pos]
        if missing:
            raise TinctureDataError(f"{len(missing)} landmark genes missing from Level 5.")
        idx = np.array([pos[g] for g in want])
        order = np.argsort(idx)
        mat = f["0/DATA/0/matrix"]
        X = np.empty((mat.shape[0], idx.size), dtype=np.float32)
        for start in range(0, mat.shape[0], block):
            chunk = mat[start : start + block, :]
            X[start : start + block] = chunk[:, idx[order]][:, np.argsort(order)]
    o = obs.set_index("sig_id").reindex(cols)
    if o.index.isna().any() or o["pert_id"].isna().any():
        raise TinctureDataError("Level 5 signature ids do not all match sig_info.")
    var = lm.set_index(lm["pr_gene_id"].astype(str))[["pr_gene_symbol", "pr_gene_title"]]
    a = ad.AnnData(X=X, obs=o.astype({c: str for c in o.columns if o[c].dtype == object}),
                   var=var)  # fmt: skip
    a.uns["source"] = SOURCE
    a.write_h5ad(out)
    tmp.unlink(missing_ok=True)
    return out


def load_lincs(
    cell_lines: list[str] | None = None,
    pert_type: str | None = "trt_cp",
    time: float | None = 24.0,
    min_tas: float | None = None,
    cache: Path | None = None,
    build: bool = False,
) -> ad.AnnData:
    """LINCS GSE70138 Level 5 landmark signatures (signatures × 978 genes) with metadata
    (``pert_id``, ``pert_iname``, ``cell_id``, ``dose`` in molar, ``time`` in hours, ``tas``,
    ``canonical_smiles``, ``inchi_key``)."""
    cache = cache or (data_home().parent / "processed" / "lincs_gse70138_landmark.h5ad")
    if not cache.exists():
        if not build:
            raise TinctureDataError(
                f"LINCS landmark cache not found at {cache}.",
                fix="Pass build=True (downloads 5 GB from GEO once) or run build_landmark_cache().",
            )
        build_landmark_cache(out=cache)
    a = ad.read_h5ad(cache)
    keep = np.ones(a.n_obs, bool)
    o = a.obs
    if pert_type is not None:
        keep &= (o["pert_type"] == pert_type).to_numpy()
    if cell_lines is not None:
        keep &= o["cell_id"].isin(cell_lines).to_numpy()
    if time is not None:
        keep &= np.isclose(o["time"].astype(float).to_numpy(), time)
    if min_tas is not None:
        keep &= (o["tas"].astype(float) >= min_tas).to_numpy()
    return a[keep].copy()
