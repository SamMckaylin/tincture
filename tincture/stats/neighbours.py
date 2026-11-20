"""Nearest annotated neighbours of a query signature, with an empirical null.

* **Similarity:** cosine between the query and each reference signature.
* **Null:** cosine similarities of distinct reference pairs in the library (the library is
  its own background). ``"all"`` (HANDOVER §5.9, default) uses every pair. ``"unrelated"``
  drops pairs that share an annotation (e.g. the same MoA): related pairs crowd the upper
  tail of the all-pairs null, exactly where neighbours are tested, which makes it
  conservative (in LINCS phase II, 24 of the 39 pairs with cosine ≥ 0.9 share an MoA).
  For a similarity ``s``: percentile = share of null pairs below ``s``; empirical
  p = (#null ≥ s + 1) / (n_null + 1).
* **Multiple testing:** BH over all references the query was compared with.
* **Confidence interval** of each similarity: from resampled query replicates (e.g. the
  query's signatures in several cell lines or doses), or from posterior draws when the
  query has per-feature SEs; otherwise none, and the plot says so.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from statsmodels.stats.multitest import multipletests

from tincture.errors import TinctureDataError
from tincture.model.reference import ReferenceLibrary


def _unit(M: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(M, axis=-1, keepdims=True)
    n[n == 0] = np.nan
    return M / n


@dataclass
class NeighbourResult:
    names: list[str]  # references, sorted by similarity (descending)
    similarity: np.ndarray
    ci: np.ndarray | None  # (n, 2) or None
    p: np.ndarray
    q: np.ndarray
    percentile: np.ndarray
    null: np.ndarray  # sorted null similarities
    n_null: int

    def confident(self, alpha: float = 0.05) -> list[str]:
        return [n for n, qq in zip(self.names, self.q, strict=True) if qq < alpha]


def library_null(
    lib: ReferenceLibrary,
    kind: str = "all",
    key: str = "moa",
    max_pairs: int = 2_000_000,
    seed: int = 0,
) -> np.ndarray:
    """Cosine similarities of distinct library pairs (subsampled above ``max_pairs``).

    ``kind="unrelated"`` keeps only pairs whose ``annotations[key]`` differ (pairs with a
    missing annotation are kept)."""
    if kind not in ("all", "unrelated"):
        raise ValueError("null kind must be 'all' or 'unrelated'.")
    U = _unit(np.nan_to_num(lib.matrix))
    n = U.shape[0]
    if n < 3:
        raise TinctureDataError("The null needs a library of at least 3 references.")
    iu = np.triu_indices(n, 1)
    if kind == "unrelated":
        lab = np.array([lib.annotations.get(r, {}).get(key, "") for r in lib.names], dtype=object)
        diff = (lab[iu[0]] != lab[iu[1]]) | (lab[iu[0]] == "") | (lab[iu[1]] == "")
        iu = (iu[0][diff], iu[1][diff])
    if iu[0].size > max_pairs:
        rng = np.random.default_rng(seed)
        pick = rng.choice(iu[0].size, max_pairs, replace=False)
        iu = (iu[0][pick], iu[1][pick])
    sims = np.einsum("ij,ij->i", U[iu[0]], U[iu[1]])
    return np.sort(sims[np.isfinite(sims)])


def neighbours(
    query: np.ndarray,
    lib: ReferenceLibrary,
    query_replicates: np.ndarray | None = None,
    query_se: np.ndarray | None = None,
    null: np.ndarray | None = None,
    n_boot: int = 200,
    seed: int = 0,
    null_kind: str = "all",
) -> NeighbourResult:
    """Rank library references by cosine similarity to ``query`` (aligned to lib.features)."""
    q = np.nan_to_num(np.asarray(query, float))
    if not np.any(q):
        raise TinctureDataError("The query signature is all zeros.")
    R = _unit(np.nan_to_num(lib.matrix))
    sims = R @ (q / np.linalg.norm(q))
    null = library_null(lib, null_kind) if null is None else np.sort(np.asarray(null, float))
    n_null = null.size
    n_ge = n_null - np.searchsorted(null, sims, side="left")
    p = (n_ge + 1) / (n_null + 1)
    pct = np.searchsorted(null, sims, side="right") / n_null
    qv = multipletests(p, method="fdr_bh")[1]
    ci = None
    rng = np.random.default_rng(seed)
    draws = None
    if query_replicates is not None and np.asarray(query_replicates).shape[0] >= 2:
        Q = np.nan_to_num(np.asarray(query_replicates, float))
        idx = rng.integers(0, Q.shape[0], size=(n_boot, Q.shape[0]))
        draws = Q[idx].mean(axis=1)
    elif query_se is not None:
        se = np.nan_to_num(np.asarray(query_se, float))
        draws = q + se * rng.standard_normal((n_boot, q.size))
    if draws is not None:
        S = _unit(draws) @ R.T  # (n_boot, n_refs)
        ci = np.column_stack([np.nanquantile(S, 0.025, axis=0), np.nanquantile(S, 0.975, axis=0)])
    order = np.argsort(-sims)
    return NeighbourResult(
        names=[lib.names[i] for i in order],
        similarity=sims[order],
        ci=None if ci is None else ci[order],
        p=p[order],
        q=qv[order],
        percentile=pct[order],
        null=null,
        n_null=n_null,
    )
