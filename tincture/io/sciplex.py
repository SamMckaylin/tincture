"""Loader for sci-Plex3 (Srivatsan et al., Science 2020; GEO GSE139944, GSM4150378).

Facts this loader relies on were verified against the files (docs/data_status.md):

* ``UMI.count.matrix.gz`` is a headerless **tab-separated** triplet file
  ``gene_index cell_index count`` (1-based), indexing rows of ``gene.annotations``
  (header line, then 110,983 genes: human ENSG + mouse ENSMUSG) and
  ``cell.annotations`` (799,317 cells, same order as ``pData``).
* ``pData`` (space-separated, quoted) has one row per cell. 36,522 cells have no hash
  assignment (NaN condition) and are dropped. Assigned cells already pass the authors'
  QC (≥ 500 UMIs, hash top-to-second ratio ≥ 2).
* Doses are nM (10, 100, 1,000, 10,000; vehicle 0); converted to molar here.
* Each ``Combo`` (well oligo + plate oligo) is one well with one condition, so the
  natural replicate unit is the well. Treated conditions have two wells (rep1, rep2);
  vehicle has 16 wells per replicate per cell line at 24 h.

The full matrix is ~3.1 GB compressed; it is streamed once into a per-well pseudobulk
(human genes only) cached as ``data/processed/sciplex3_pseudobulk_wells.h5ad``.
"""

from __future__ import annotations

import gzip
import re
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from tincture.datasets.cache import data_home, fetch
from tincture.datasets.registry import SCIPLEX3
from tincture.errors import TinctureDataError
from tincture.model.perturbation import PerturbationSet

SOURCE = "sci-Plex3, Srivatsan et al. Science 2020, GEO GSE139944 (GSM4150378)"
_CXSMILES_SUFFIX = re.compile(r"\s*\|.*\|\s*$")


def _paths(root: Path | None, download: bool) -> dict[str, Path]:
    return {k: fetch(f, root=root, download=download) for k, f in SCIPLEX3.items()}


def clean_smiles(s: object) -> str | None:
    """Strip a ChemAxon CXSMILES extension (``|...|``) to plain SMILES."""
    if not isinstance(s, str) or s in ("", "NA"):
        return None
    return _CXSMILES_SUFFIX.sub("", s).strip() or None


def read_cell_metadata(root: Path | None = None, download: bool = True) -> pd.DataFrame:
    """Per-cell metadata in tincture's vocabulary, indexed by the 0-based matrix column.

    Columns: ``cell``, ``well``, ``plate``, ``cell_line``, ``compound_id``,
    ``compound_name``, ``dose`` (M), ``time`` (h), ``replicate``, ``target``, ``pathway``,
    ``pathway_level_1``, ``pathway_level_2``, ``smiles``, ``n_umi``, ``assigned``.
    """
    paths = _paths(root, download)
    p = pd.read_csv(paths["pdata"], sep=" ", low_memory=False)
    cells = pd.read_csv(paths["cell_annotations"], sep="\t", header=None)
    if len(p) != len(cells) or not (p["cell"].to_numpy() == cells[0].to_numpy()).all():
        raise TinctureDataError(
            "sci-Plex3 pData and cell.annotations are not in the same order.",
            fix="Re-download both files; the loader relies on the published order.",
        )
    hm = pd.read_csv(paths["hash_metadata"], sep="\t")
    smiles = (
        hm.loc[~hm["vehicle"].astype(bool), ["treatment", "SMILES"]]
        .drop_duplicates("treatment")
        .set_index("treatment")["SMILES"]
        .map(clean_smiles)
    )
    vehicle = p["vehicle"].astype("boolean").fillna(False).astype(bool)
    out = pd.DataFrame(
        {
            "cell": p["cell"].astype(str),
            "well": p["Combo"].astype(str),
            "plate": p["culture_plate"].astype(str),
            "cell_line": p["cell_type"],
            "compound_id": np.where(vehicle, "vehicle", p["catalog_number"].astype(str)),
            "compound_name": np.where(vehicle, "vehicle (DMSO)", p["product_name"].astype(str)),
            "dose": np.where(vehicle, 0.0, pd.to_numeric(p["dose"], errors="coerce") / 1e9),
            "time": pd.to_numeric(p["time_point"], errors="coerce"),
            "replicate": p["replicate"],
            "target": p["target"],
            "pathway": p["pathway"],
            "pathway_level_1": p["pathway_level_1"],
            "pathway_level_2": p["pathway_level_2"],
            "n_umi": pd.to_numeric(p["n.umi"], errors="coerce"),
        }
    )
    out["smiles"] = out["compound_id"].map(smiles)
    out["assigned"] = p["cell_type"].notna().to_numpy()
    out.loc[~out["assigned"], ["compound_id", "compound_name"]] = np.nan
    return out


def read_genes(root: Path | None = None, download: bool = True) -> pd.DataFrame:
    """Gene table indexed by 0-based matrix row, with a ``human`` flag."""
    paths = _paths(root, download)
    g = pd.read_csv(paths["gene_annotations"], sep=" ")
    g.columns = ["gene_id", "gene_symbol"]
    g["human"] = g["gene_id"].str.startswith("ENSG")
    return g


def _well_obs(cells: pd.DataFrame) -> pd.DataFrame:
    a = cells[cells["assigned"]]
    cols = [
        "well",
        "plate",
        "cell_line",
        "compound_id",
        "compound_name",
        "dose",
        "time",
        "replicate",
        "target",
        "pathway",
        "pathway_level_1",
        "pathway_level_2",
        "smiles",
    ]
    obs = a[cols].drop_duplicates("well").set_index("well", drop=False)
    obs["n_cells"] = a.groupby("well").size().reindex(obs.index).to_numpy()
    obs["total_umi"] = a.groupby("well")["n_umi"].sum().reindex(obs.index).to_numpy()
    obs.index.name = None
    return obs.sort_values(["cell_line", "time", "compound_id", "dose", "replicate", "well"])


def build_pseudobulk(
    root: Path | None = None,
    out: Path | None = None,
    chunksize: int = 20_000_000,
    progress: Callable[[str], None] | None = print,
) -> Path:
    """Stream the triplet matrix once and sum counts per well (human genes only).

    Writes an ``.h5ad`` with raw summed counts in ``X`` (CSR, float32) and returns its path.
    """
    paths = _paths(root, download=True)
    out = out or (data_home().parent / "processed" / "sciplex3_pseudobulk_wells.h5ad")
    out.parent.mkdir(parents=True, exist_ok=True)
    cells = read_cell_metadata(root, download=False)
    genes = read_genes(root, download=False)
    obs = _well_obs(cells)
    well_index = pd.Series(np.arange(len(obs)), index=obs.index)
    cell_to_well = np.full(len(cells) + 1, -1, dtype=np.int64)  # 1-based lookup
    assigned = cells["assigned"].to_numpy()
    cell_to_well[1:][assigned] = well_index.reindex(cells.loc[assigned, "well"]).to_numpy()
    human = genes["human"].to_numpy()
    gene_to_col = np.full(len(genes) + 1, -1, dtype=np.int64)
    gene_to_col[1:][human] = np.arange(int(human.sum()))
    n_wells, n_genes = len(obs), int(human.sum())

    total = sp.csr_matrix((n_wells, n_genes), dtype=np.float64)
    n_lines = 0
    t0 = time.time()
    with gzip.open(paths["umi_matrix"], "rt") as fh:
        reader = pd.read_csv(
            fh,
            sep="\t",
            header=None,
            names=["g", "c", "n"],
            dtype=np.int64,
            chunksize=chunksize,
            engine="c",
        )
        for chunk in reader:
            g = gene_to_col[chunk["g"].to_numpy()]
            w = cell_to_well[chunk["c"].to_numpy()]
            keep = (g >= 0) & (w >= 0)
            part = sp.coo_matrix(
                (chunk["n"].to_numpy()[keep].astype(np.float64), (w[keep], g[keep])),
                shape=(n_wells, n_genes),
            ).tocsr()
            total = total + part
            n_lines += len(chunk)
            if progress:
                progress(f"  {n_lines / 1e6:,.0f}M triplets, {time.time() - t0:,.0f}s")
    var = genes.loc[human, ["gene_id", "gene_symbol"]].reset_index(drop=True)
    var.index = var["gene_id"].astype(str)
    adata = ad.AnnData(X=total.astype(np.float32), obs=obs, var=var)
    adata.uns["tincture_source"] = SOURCE
    adata.uns["n_triplets"] = n_lines
    adata.uns["dose_unit_note"] = "dose converted from nM to molar"
    for c in adata.obs.columns:
        if adata.obs[c].dtype == object:
            adata.obs[c] = adata.obs[c].astype(str)
    adata.write_h5ad(out, compression="gzip")
    return out


def load_sciplex3(
    root: Path | None = None,
    cache: Path | None = None,
    cell_lines: list[str] | None = None,
    time: float | list[float] | None = 24.0,
    compounds: list[str] | None = None,
    build: bool = False,
) -> PerturbationSet:
    """sci-Plex3 as a per-well pseudobulk :class:`PerturbationSet` (raw summed UMI counts).

    Parameters
    ----------
    cell_lines, time, compounds
        Optional filters. ``time=24`` (default) keeps the 24 h screen; A549 also has 72 h.
        ``compounds`` are catalog ids or names; vehicle is always kept.
    build
        Build the pseudobulk cache if missing (streams the 3.1 GB matrix; tens of minutes).
    """
    cache = cache or (data_home().parent / "processed" / "sciplex3_pseudobulk_wells.h5ad")
    if not cache.exists():
        if not build:
            raise TinctureDataError(
                f"sci-Plex3 pseudobulk cache not found at {cache}.",
                fix="Run tincture.io.sciplex.build_pseudobulk() once (or pass build=True).",
            )
        build_pseudobulk(root, cache)
    adata = ad.read_h5ad(cache)
    # Text columns were stored as strings; restore real missing values (e.g. vehicle SMILES).
    for c in ("smiles", "target", "pathway", "pathway_level_1", "pathway_level_2"):
        if c in adata.obs:
            col = adata.obs[c].astype(object)
            adata.obs[c] = col.where(~col.isin(["nan", "None", ""]), None)
    ps = PerturbationSet(adata, source=SOURCE)
    sel: dict[str, Any] = {}
    if cell_lines is not None:
        sel["cell_line"] = cell_lines
    if time is not None:
        sel["time"] = time if isinstance(time, list) else [time]
    if sel:
        ps = ps.subset(**sel)
    if compounds is not None:
        o = ps.obs
        keep = (
            o["compound_id"].isin(compounds) | o["compound_name"].isin(compounds) | (o["dose"] == 0)
        )
        found = set(o.loc[keep & (o["dose"] > 0), "compound_id"]) | set(
            o.loc[keep & (o["dose"] > 0), "compound_name"]
        )
        unknown = sorted(set(compounds) - found)
        if unknown:
            raise TinctureDataError(
                f"Compounds not in sci-Plex3: {unknown[:5]}.",
                fix="Use catalog ids (e.g. 'S1090') or product names from obs['compound_name'].",
            )
        ps = ps.subset(obs=keep.to_numpy())
    return ps
