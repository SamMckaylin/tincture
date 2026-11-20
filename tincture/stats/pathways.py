"""Pathway activity scores: univariate linear model (ULM), as in decoupler.

For a signature ``y`` (features) and a gene set with weights ``w`` (1 for members, 0
otherwise; PROGENy supplies signed weights), fit ``y = a + b·w`` across all features and
report the t-statistic of ``b`` (score) and its two-sided p-value (df = n − 2). This is
decoupler's ``ulm``; implemented here in matrix form so it runs without the ``[omics]``
extra, and checked against decoupler in the test suite when it is installed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from tincture.errors import TinctureDataError


@dataclass
class PathwayScores:
    """Scores and p-values, signatures × pathways."""

    score: pd.DataFrame
    pvalue: pd.DataFrame
    n_features: pd.Series  # members present per pathway
    source: str = ""


def network_matrix(
    net: pd.DataFrame, features: np.ndarray, min_size: int = 5
) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Features × pathways weight matrix from a long ``source, target[, weight]`` table."""
    if not {"source", "target"} <= set(net.columns):
        raise TinctureDataError("Gene-set table needs 'source' and 'target' columns.")
    w = net["weight"] if "weight" in net.columns else pd.Series(1.0, index=net.index)
    n = net.assign(weight=w.astype(float)).drop_duplicates(["source", "target"])
    pos = {f: i for i, f in enumerate(features)}
    n = n[n["target"].isin(pos)]
    sources = sorted(n["source"].unique())
    sizes = n.groupby("source").size()
    sources = [s for s in sources if sizes[s] >= min_size]
    M = np.zeros((len(features), len(sources)))
    if not sources:
        return M, sources, np.zeros(0, dtype=int)
    col = {s: j for j, s in enumerate(sources)}
    n = n[n["source"].isin(col)]
    rows = n["target"].map(pos).to_numpy(dtype=int)
    cols = n["source"].map(col).to_numpy(dtype=int)
    M[rows, cols] = n["weight"].to_numpy()
    return M, sources, sizes.reindex(sources).to_numpy()


def ulm(Y: pd.DataFrame, net: pd.DataFrame, min_size: int = 5, source: str = "") -> PathwayScores:
    """ULM scores for each row of ``Y`` (signatures × features) and each gene set."""
    feats = Y.columns.to_numpy().astype(str)
    M, sources, _ = network_matrix(net, feats, min_size)
    if not sources:
        raise TinctureDataError(
            f"No gene set has ≥ {min_size} members among the {len(feats)} features.",
            fix="Check feature identifiers match the gene-set targets (e.g. symbols vs Ensembl).",
        )
    y = Y.to_numpy(dtype=float)
    ok = np.isfinite(y).all(axis=0)
    y, M = y[:, ok], M[ok]
    n = y.shape[1]
    xc = M - M.mean(axis=0)
    yc = y - y.mean(axis=1, keepdims=True)
    sxx = (xc**2).sum(axis=0)
    sxy = yc @ xc
    b = sxy / sxx
    syy = (yc**2).sum(axis=1, keepdims=True)
    rss = np.maximum(syy - b * sxy, 1e-300)
    se_b = np.sqrt(rss / (n - 2) / sxx)
    t = b / se_b
    p = 2 * stats.t.sf(np.abs(t), df=n - 2)
    return PathwayScores(
        score=pd.DataFrame(t, index=Y.index, columns=sources),
        pvalue=pd.DataFrame(p, index=Y.index, columns=sources),
        n_features=pd.Series((M != 0).sum(axis=0), index=sources),
        source=source,
    )


def dedup_gene_sets(
    net: pd.DataFrame, max_jaccard: float = 0.5
) -> tuple[pd.DataFrame, list[tuple[str, str, float]]]:
    """Drop the smaller of any pair of sets with Jaccard > ``max_jaccard``.

    Returns the filtered table and the list of (dropped, kept, jaccard) decisions.
    """
    members = net.groupby("source")["target"].agg(set)
    order = members.map(len).sort_values(ascending=False).index.tolist()
    kept: list[str] = []
    dropped: list[tuple[str, str, float]] = []
    for s in order:
        hit = None
        for k in kept:
            a, b = members[s], members[k]
            j = len(a & b) / len(a | b)
            if j > max_jaccard:
                hit = (s, k, j)
                break
        if hit:
            dropped.append(hit)
        else:
            kept.append(s)
    return net[net["source"].isin(kept)].copy(), dropped


def symbols_for(features: np.ndarray, symbols: np.ndarray | None) -> np.ndarray:
    """Map feature ids to gene symbols (first occurrence wins; unmapped keep their id)."""
    if symbols is None:
        return np.asarray(features).astype(str)
    s = pd.Series(np.asarray(symbols).astype(str))
    dup = s.duplicated()
    out = s.where(~dup, pd.Series(np.asarray(features).astype(str)))
    return out.to_numpy()
