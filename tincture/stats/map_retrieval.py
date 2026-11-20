"""Replicate retrievability: mean average precision (mAP) with a permutation null.

Definitions follow the JUMP Cell Painting ``copairs`` framework (Kalinin et al.;
github.com/cytomining/copairs), "phenotypic activity" set-up:

* For each non-control profile (query), **positives** are the other replicates of the same
  group (e.g. compound × dose); **negatives** are control profiles (``negatives="controls"``,
  default) or all profiles outside the group (``negatives="all"``).
* Rank positives and negatives by cosine similarity to the query (descending).
  **AP** = mean over positives of precision at that positive's rank.
* **mAP** of a group = mean AP of its members.
* **Null.** Two options:

  - ``"permutation"`` (default): under H0 a group's replicates are exchangeable with the
    control profiles. Relabel which members of the pool (group ∪ controls) form the group
    and recompute mAP from the *observed* similarity matrix. This keeps the dependence
    between replicates' APs (they share sim(a, b)). When the number of relabelings
    C(pool, n_rep) ≤ ``null_size`` all are enumerated (exact p; the observed labelling is
    one of them, so ``p = #{null ≥ mAP} / N``); otherwise Monte Carlo with
    ``p = (#{null ≥ mAP} + 1) / (N + 1)``.
  - ``"copied"`` (copairs ≤ 0.5.5, and copairs 0.6 ``null="copied"``): each member gets
    the AP of a random ranking with the same numbers of positives and candidates; members
    with equal configuration share draws; null mAP = their mean. With 2 replicates and
    n controls, p can never fall below 1/(n + 1), so duplicate designs cannot pass FDR
    control over many conditions. Kept for exact parity with copairs.

  ``negatives="all"`` supports only ``"copied"``.
* **q**: Benjamini–Hochberg across groups.

Agreement with copairs 0.6.0 is tested against a frozen golden file
(``tests/fixtures/copairs_golden.json``); copairs itself is not a dependency because it
pulls in duckdb (~65 MB).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from tincture.errors import TinctureDataError


def _cosine(X: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(X, axis=1, keepdims=True)
    if np.any(norm == 0):
        raise TinctureDataError(
            "Some profiles are all-zero; cosine similarity is undefined.",
            fix="Drop empty profiles or check the feature scaling.",
        )
    Z = X / norm
    return Z @ Z.T


def ap_from_relevance(rel: np.ndarray) -> np.ndarray:
    """AP for each row of a 0/1 relevance matrix already sorted by decreasing score."""
    rel = np.asarray(rel, dtype=float)
    hits = np.cumsum(rel, axis=1)
    ranks = np.arange(1, rel.shape[1] + 1)
    prec = hits / ranks
    n_pos = rel.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (prec * rel).sum(axis=1) / n_pos


def random_ap(n_draws: int, n_pos: int, total: int, rng: np.random.Generator) -> np.ndarray:
    """Null APs: random rankings with ``n_pos`` positives among ``total`` items."""
    keys = rng.random((n_draws, total))
    order = np.argsort(keys, axis=1)
    rel = order < n_pos  # first n_pos item ids are the positives
    return ap_from_relevance(rel)


def average_precision(
    X: np.ndarray,
    groups: Sequence[str] | np.ndarray,
    is_control: np.ndarray,
    negatives: str = "controls",
) -> pd.DataFrame:
    """Per-profile AP (one row per non-control profile that has ≥ 1 positive)."""
    X = np.asarray(X, dtype=float)
    g = np.asarray(groups).astype(str)
    ctrl = np.asarray(is_control, dtype=bool)
    if negatives not in ("controls", "all"):
        raise ValueError("negatives must be 'controls' or 'all'.")
    if negatives == "controls" and ctrl.sum() == 0:
        raise TinctureDataError(
            "No control profiles: retrieval against controls is undefined.",
            fix="Include vehicle/negative-control profiles, or use negatives='all'.",
        )
    if not np.isfinite(X).all():
        raise TinctureDataError(
            "Profiles contain NaN or inf.", fix="Impute or drop those features."
        )
    S = _cosine(X)
    rows = []
    for i in np.flatnonzero(~ctrl):
        pos = np.flatnonzero((g == g[i]) & ~ctrl)
        pos = pos[pos != i]
        if pos.size == 0:
            continue
        neg = (
            np.flatnonzero(ctrl)
            if negatives == "controls"
            else np.flatnonzero((g != g[i]) & ~ctrl | ctrl)
        )
        neg = neg[neg != i]
        cand = np.concatenate([pos, neg])
        sims = S[i, cand]
        order = np.argsort(-sims, kind="stable")
        rel = np.zeros(cand.size)
        rel[: pos.size] = 1
        ap = ap_from_relevance(rel[order][None, :])[0]
        rows.append((i, g[i], float(ap), int(pos.size), int(cand.size)))
    return pd.DataFrame(rows, columns=["profile", "group", "ap", "n_pos", "n_total"])


@dataclass
class MapResult:
    """Per-group retrievability."""

    table: pd.DataFrame  # group, map, p, q, retrievable, n_profiles, null_lo, null_hi
    ap: pd.DataFrame  # per-profile AP
    null_size: int
    alpha: float
    negatives: str
    excluded: list[str]  # groups with < 2 replicates
    null: str = "copied"
    min_p: np.ndarray | None = None  # smallest achievable p per group (design limit)
    design_can_pass: bool = True


def _ap_rows(S_pool: np.ndarray, members: np.ndarray) -> np.ndarray:
    """mAP for each row of ``members`` (draws × n_rep indices into the pool)."""
    n_draw, n_rep = members.shape
    m = S_pool.shape[0]
    in_group = np.zeros((n_draw, m), dtype=bool)
    np.put_along_axis(in_group, members, True, axis=1)
    aps = np.zeros((n_draw, n_rep))
    for j in range(n_rep):
        q = members[:, j]
        sims = S_pool[q].copy()
        sims[np.arange(n_draw), q] = -np.inf  # exclude self
        order = np.argsort(-sims, axis=1, kind="stable")[:, : m - 1]
        rel = np.take_along_axis(in_group, order, axis=1)
        aps[:, j] = ap_from_relevance(rel)
    return aps.mean(axis=1)


def permutation_null(
    S: np.ndarray,
    group_idx: np.ndarray,
    control_idx: np.ndarray,
    null_size: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, bool]:
    """Null mAPs by relabelling group members among group ∪ controls.

    Returns (null values, exact): ``exact`` means all C(pool, n_rep) relabelings were
    enumerated (including the observed labelling).
    """
    from itertools import combinations
    from math import comb

    pool = np.concatenate([group_idx, control_idx])
    S_pool = S[np.ix_(pool, pool)]
    m, k = pool.size, group_idx.size
    if comb(m, k) <= null_size:
        members = np.array(list(combinations(range(m), k)), dtype=int)
        return _ap_rows(S_pool, members), True
    members = np.argsort(rng.random((null_size, m)), axis=1)[:, :k]
    return _ap_rows(S_pool, members), False


def min_achievable_p(n_rep: int, n_controls: int, null: str, null_size: int) -> float:
    """Smallest p-value the design can produce (perfect retrieval)."""
    from math import comb

    if null == "copied":
        return max(1.0 / comb(n_rep - 1 + n_controls, n_rep - 1), 1.0 / (null_size + 1))
    n_comb = comb(n_rep + n_controls, n_rep)
    return 1.0 / n_comb if n_comb <= null_size else 1.0 / (null_size + 1)


def mean_average_precision(
    ap: pd.DataFrame,
    null_size: int = 10_000,
    alpha: float = 0.05,
    seed: int = 0,
    all_groups: Sequence[str] | None = None,
    negatives: str = "controls",
) -> MapResult:
    """mAP per group with the copied permutation null and BH q-values."""
    if null_size < 100:
        raise ValueError("null_size must be ≥ 100 for a usable p-value resolution.")
    rng = np.random.default_rng(seed)
    confs = ap[["n_pos", "n_total"]].drop_duplicates().itertuples(index=False)
    nulls = {(int(a), int(b)): random_ap(null_size, int(a), int(b), rng) for a, b in confs}
    rows = []
    for grp, d in ap.groupby("group", sort=True):
        m = float(d["ap"].mean())
        null = np.mean(
            [nulls[(int(a), int(b))] for a, b in zip(d["n_pos"], d["n_total"], strict=True)], axis=0
        )
        tol = 1e-6
        p = (np.sum(null >= m - tol) + 1) / (null_size + 1)
        lo, hi = np.quantile(null, [0.05, 0.95])
        rows.append((grp, m, p, len(d), float(lo), float(hi)))
    t = pd.DataFrame(rows, columns=["group", "map", "p", "n_profiles", "null_lo", "null_hi"])
    t["q"] = multipletests(t["p"], method="fdr_bh")[1] if len(t) else []
    t["retrievable"] = t["q"] < alpha
    excluded = sorted(set(all_groups or []) - set(t["group"]))
    return MapResult(
        table=t.sort_values("map", ascending=False).reset_index(drop=True),
        ap=ap,
        null_size=null_size,
        alpha=alpha,
        negatives=negatives,
        excluded=excluded,
    )


def retrievability(
    X: np.ndarray,
    groups: Sequence[str] | np.ndarray,
    is_control: np.ndarray,
    negatives: str = "controls",
    null_size: int = 10_000,
    alpha: float = 0.05,
    seed: int = 0,
    null: str = "permutation",
) -> MapResult:
    """AP per profile + mAP, p, q per group in one call (see module docstring for nulls)."""
    if null not in ("permutation", "copied"):
        raise ValueError("null must be 'permutation' or 'copied'.")
    if null == "permutation" and negatives != "controls":
        raise ValueError("The permutation null is defined for negatives='controls' only.")
    g = np.asarray(groups).astype(str)
    ctrl = np.asarray(is_control, dtype=bool)
    ap = average_precision(X, g, ctrl, negatives)
    if ap.empty:
        raise TinctureDataError(
            "No group has ≥ 2 replicate profiles, so nothing can be retrieved.",
            fix="Retrievability needs replicates (e.g. two wells per compound × dose).",
        )
    all_groups = np.unique(g[~ctrl]).tolist()
    if null == "copied":
        res = mean_average_precision(ap, null_size, alpha, seed, all_groups, negatives)
        res.null = "copied"
    else:
        rng = np.random.default_rng(seed)
        S = _cosine(np.asarray(X, dtype=float))
        ctrl_idx = np.flatnonzero(ctrl)
        rows = []
        for grp, d in ap.groupby("group", sort=True):
            m = float(d["ap"].mean())
            members = d["profile"].to_numpy()
            null_vals, exact = permutation_null(S, members, ctrl_idx, null_size, rng)
            n_ge = int(np.sum(null_vals >= m - 1e-9))
            p = n_ge / null_vals.size if exact else (n_ge + 1) / (null_vals.size + 1)
            lo, hi = np.quantile(null_vals, [0.05, 0.95])
            rows.append((grp, m, p, len(d), float(lo), float(hi), exact))
        t = pd.DataFrame(
            rows, columns=["group", "map", "p", "n_profiles", "null_lo", "null_hi", "exact"]
        )
        t["q"] = multipletests(t["p"], method="fdr_bh")[1]
        t["retrievable"] = t["q"] < alpha
        res = MapResult(
            table=t.sort_values("map", ascending=False).reset_index(drop=True),
            ap=ap,
            null_size=null_size,
            alpha=alpha,
            negatives=negatives,
            excluded=sorted(set(all_groups) - set(t["group"])),
            null="permutation",
        )
    n_ctrl = int(ctrl.sum()) if negatives == "controls" else int((~ctrl).sum())
    reps = res.table["n_profiles"].to_numpy()
    res.min_p = np.array([min_achievable_p(int(r), n_ctrl, res.null, null_size) for r in reps])
    # Best case for BH: every condition at its minimum p.
    n = len(res.table)
    res.design_can_pass = bool(n) and _bh_reachable(res.min_p, alpha)
    return res


def _bh_reachable(min_p: np.ndarray, alpha: float) -> bool:
    """Can any condition reach q < alpha if every condition achieved its minimum p?"""
    q = multipletests(min_p, method="fdr_bh")[1]
    return bool(np.any(q < alpha))


# ----------------------------------------------------------------------------- CAT curve


@dataclass
class CatCurve:
    k: np.ndarray
    observed: np.ndarray  # fraction of top-k features shared
    null_lo: np.ndarray
    null_hi: np.ndarray
    null_mean: np.ndarray
    n_features: int


def cat_curve(
    a: np.ndarray, b: np.ndarray, ks: Sequence[int] | None = None, level: float = 0.95
) -> CatCurve:
    """Correspondence at the top: fraction of the top-k features (by |value|) shared
    between two replicate profiles, with an exact hypergeometric null band (overlap of
    two random k-subsets of n features)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.shape != b.shape or a.ndim != 1:
        raise TinctureDataError("CAT curve needs two 1-D profiles over the same features.")
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    n = a.size
    if ks is None:
        k_arr = np.unique(
            np.clip(np.round(np.logspace(0.7, np.log10(max(n // 2, 6)), 25)), 5, n)
        ).astype(int)
    else:
        k_arr = np.asarray(ks, dtype=int)
    ra = np.argsort(-np.abs(a), kind="stable")
    rb = np.argsort(-np.abs(b), kind="stable")
    obs = np.array([len(np.intersect1d(ra[:k], rb[:k], assume_unique=True)) / k for k in k_arr])
    lo = np.array([stats.hypergeom.ppf((1 - level) / 2, n, k, k) / k for k in k_arr])
    hi = np.array([stats.hypergeom.ppf(1 - (1 - level) / 2, n, k, k) / k for k in k_arr])
    mean = k_arr / n
    return CatCurve(k=k_arr, observed=obs, null_lo=lo, null_hi=hi, null_mean=mean, n_features=n)
