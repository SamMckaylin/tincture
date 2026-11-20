"""MIX-seq (McFarland, Paolella et al., Nat Commun 2020; figshare 10298696 v3). Facts in
docs/data_status.md (D10).

:func:`pseudobulk_trametinib_tc` streams the 10x matrix of experiment 5 (trametinib 0.1 µM vs
DMSO, hash-tagged at 3, 6, 12, 24 and 48 h, 24 pooled cell lines) straight from the zip and
sums UMI counts per cell line × hash tag, keeping QC-passing singlets
(``cell_quality == "normal"``) with a single hash tag. The result is cached in
``data/processed/mixseq/`` (parquet), keyed by the zip's checksum.
"""

from __future__ import annotations

import io
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from tincture.datasets.cache import data_home, fetch
from tincture.datasets.registry import MIXSEQ
from tincture.errors import TinctureDataError

_DIR = "trametinib_tc_expt5"
_HASH = re.compile(r"^(DMSO|Tram|Untreated)_(\d+)hr$")


def _cache_dir() -> Path:
    p = data_home().parent / "processed" / "mixseq"
    p.mkdir(parents=True, exist_ok=True)
    return p


def pseudobulk_trametinib_tc(min_cells: int = 10) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns ``(counts, groups)``: counts = genes (symbol index, Ensembl in column
    ``ensembl``… dropped) × groups; groups = one row per column with ``depmap_id``,
    ``condition`` (DMSO / trametinib / untreated), ``time_h`` and ``n_cells``. Groups with
    fewer than ``min_cells`` cells are dropped."""
    zpath = fetch(MIXSEQ["trametinib_tc"])
    tag = MIXSEQ["trametinib_tc"].md5 or "x"
    cfile = _cache_dir() / f"trametinib_tc_pseudobulk_{tag[:8]}.parquet"
    gfile = cfile.with_suffix(".groups.parquet")
    if not (cfile.exists() and gfile.exists()):
        _build(zpath, cfile, gfile)
    counts = pd.read_parquet(cfile)
    groups = pd.read_parquet(gfile)
    keep = groups["n_cells"] >= min_cells
    groups = groups[keep].reset_index(drop=True)
    return counts[groups["group"]], groups


def _build(zpath: Path, cfile: Path, gfile: Path) -> None:
    z = zipfile.ZipFile(zpath)
    bc = pd.read_csv(z.open(f"{_DIR}/barcodes.tsv"), header=None)[0]
    genes = pd.read_csv(z.open(f"{_DIR}/genes.tsv"), header=None, sep="\t", names=["ens", "sym"])
    cl = pd.read_csv(z.open(f"{_DIR}/classifications.csv")).set_index("barcode").reindex(bc)
    ok = (cl["cell_quality"] == "normal") & cl["hash_tag"].astype(str).str.match(_HASH)
    key = (cl["DepMap_ID"].astype(str) + "|" + cl["hash_tag"].astype(str)).where(ok)
    codes, uniq = pd.factorize(key)  # -1 for excluded cells
    n_g, n_grp = len(genes), len(uniq)
    if n_grp == 0:
        raise TinctureDataError("No QC-passing singlet cells found in the MIX-seq file.")
    acc = np.zeros(n_grp * n_g, dtype=np.float64)
    with z.open(f"{_DIR}/matrix.mtx") as fh:
        text = io.TextIOWrapper(fh)
        header = text.readline()
        if "coordinate" not in header:
            raise TinctureDataError("Unexpected MatrixMarket header in MIX-seq matrix.")
        line = text.readline()
        while line.startswith("%"):
            line = text.readline()
        for ch in pd.read_csv(
            text, sep=" ", header=None, names=["g", "c", "v"], dtype=np.int64, chunksize=5_000_000
        ):
            gi = ch["g"].to_numpy() - 1
            grp = codes[ch["c"].to_numpy() - 1]
            m = grp >= 0
            acc += np.bincount(grp[m] * n_g + gi[m], weights=ch["v"].to_numpy()[m],
                               minlength=n_grp * n_g)  # fmt: skip
    mat = acc.reshape(n_grp, n_g).T
    names = [f"g{i}" for i in range(n_grp)]
    counts = pd.DataFrame(mat, index=genes["sym"].astype(str), columns=names)
    counts = counts.groupby(level=0).sum()  # duplicated symbols summed
    rows = []
    n_cells = pd.Series(codes[codes >= 0]).value_counts()
    for i, k in enumerate(uniq):
        dep, tag = str(k).split("|")
        mm = _HASH.match(tag)
        assert mm is not None
        cond = {"DMSO": "DMSO", "Tram": "trametinib", "Untreated": "untreated"}[mm.group(1)]
        rows.append(
            {"group": names[i], "depmap_id": dep, "hash_tag": tag, "condition": cond,
             "time_h": float(mm.group(2)), "n_cells": int(n_cells.get(i, 0))}
        )  # fmt: skip
    counts.to_parquet(cfile)
    pd.DataFrame(rows).to_parquet(gfile)
