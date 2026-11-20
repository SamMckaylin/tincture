"""Perturbation signatures from replicate-level count profiles.

Pipeline (method ``"logcpm_eb"``, the default):

1. **Strata.** Signatures are always computed against vehicle samples of the *same*
   stratum (cell line × time, by default), never across cell lines.
2. **Filter.** Keep features with CPM ≥ ``min_cpm`` in at least ``min_frac`` of the
   stratum's samples.
3. **Transform.** ``log2((count + 0.5) / (lib + 1) · 1e6)`` (logCPM).
4. **Linear model.** One group mean per condition (compound × dose); residual variance
   pooled across all groups in the stratum.
5. **Moderation.** Empirical-Bayes squeezing of residual variances with a mean-variance
   trend (Smyth 2004; limma-trend), essential with 2 replicates per condition.
6. **Effect.** Δ logCPM vs vehicle, SE = √(s̃² (1/n_c + 1/n_v)), moderated t, p.
7. **Shrinkage.** Normal-prior empirical-Bayes shrinkage of each signature's effects
   (prior variance by method of moments); noisy features shrink most. Posterior SD is
   stored as ``se``.

Uncertainty: ``Signature.se`` is the posterior SD. Draws for downstream bootstrap CIs are
generated on demand from this posterior (seeded); with only two replicates per condition
a replicate bootstrap would have three distinct resamples and would understate
uncertainty.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from tincture.errors import TinctureDataError, TinctureWarning
from tincture.model.perturbation import PerturbationSet
from tincture.model.signature import Signature
from tincture.stats.ebayes import squeeze_var

DEFAULT_STRATA: tuple[str, ...] = ("cell_line", "time")
CONDITION: tuple[str, ...] = ("compound_id", "dose")


def logcpm(counts: np.ndarray, prior: float = 0.5) -> np.ndarray:
    """log2 counts per million (samples in rows)."""
    lib = counts.sum(axis=1, keepdims=True)
    if np.any(lib <= 0):
        raise TinctureDataError(
            "Some samples have zero total counts.",
            fix="Remove empty samples before computing signatures.",
        )
    return np.log2((counts + prior) / (lib + 1.0) * 1e6)


def shrink_normal(effect: np.ndarray, se: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Normal-normal empirical-Bayes shrinkage of effects towards 0.

    Prior effect ~ N(0, τ²), τ² = max(mean(effect² − se²), floor). Returns posterior mean,
    posterior SD and τ².
    """
    ok = np.isfinite(effect) & np.isfinite(se)
    if ok.sum() == 0:
        return effect.copy(), se.copy(), float("nan")
    tau2 = float(np.mean(effect[ok] ** 2 - se[ok] ** 2))
    floor = 1e-4 * float(np.median(se[ok] ** 2))
    tau2 = max(tau2, floor)
    w = tau2 / (tau2 + se**2)
    return effect * w, np.sqrt(w * se**2), tau2


@dataclass
class SignatureSet:
    """All signatures of one stratum plus the shared statistics behind them."""

    signatures: list[Signature]
    stratum: dict[str, Any]
    features: np.ndarray
    df_prior: float
    df_residual: float
    n_vehicle: int
    method: str
    pvalues: dict[tuple[str, float], np.ndarray] = field(default_factory=dict)
    feature_symbols: np.ndarray | None = None

    def get(self, compound_id: str, dose: float) -> Signature:
        for s in self.signatures:
            if s.compound_id == compound_id and np.isclose(s.dose, dose):
                return s
        raise KeyError((compound_id, dose))

    def matrix(self, shrunk: bool = True) -> pd.DataFrame:
        """Signatures × features effect matrix (shrunk by default)."""
        key = "effect" if shrunk else "raw_effect"
        rows = {
            f"{s.compound_id}@{s.dose:g}": (s.effect if shrunk else s.meta[key])
            for s in self.signatures
        }
        return pd.DataFrame(rows, index=self.features).T


def _dense(X: Any) -> np.ndarray:
    return np.asarray(X.toarray() if sp.issparse(X) else X, dtype=float)


def compute_signatures(
    ps: PerturbationSet,
    strata: Sequence[str] = DEFAULT_STRATA,
    min_cpm: float = 5.0,
    min_frac: float = 0.25,
    trend: bool = True,
    shrink: bool = True,
    min_replicates: int = 1,
) -> list[SignatureSet]:
    """Compute one :class:`SignatureSet` per stratum (see module docstring).

    Raises
    ------
    TinctureDataError
        No vehicle samples in a stratum, or no feature passes the expression filter.
    """
    import warnings

    obs = ps.obs
    strata = [s for s in strata if s in obs.columns]
    keys = (
        obs[strata].astype(str).agg("|".join, axis=1)
        if strata
        else pd.Series("all", index=obs.index)
    )
    out: list[SignatureSet] = []
    X_all = ps.X
    for key in sorted(keys.unique()):
        idx = np.flatnonzero((keys == key).to_numpy())
        o = obs.iloc[idx]
        stratum = dict(zip(strata, key.split("|"), strict=True)) if strata else {}
        is_veh = (o["dose"] == 0).to_numpy()
        if is_veh.sum() == 0:
            raise TinctureDataError(
                f"No vehicle samples in stratum {stratum}.",
                fix="Signatures are computed against same-stratum vehicle; add vehicle wells "
                "or change `strata`.",
            )
        counts = _dense(X_all[idx]) if not sp.issparse(X_all) else _dense(X_all[idx, :])
        cpm = counts / counts.sum(axis=1, keepdims=True) * 1e6
        keep = (cpm >= min_cpm).mean(axis=0) >= min_frac
        if keep.sum() < 10:
            raise TinctureDataError(
                f"Only {int(keep.sum())} features pass the expression filter in {stratum}.",
                fix="Lower min_cpm/min_frac, or check that X holds raw counts.",
            )
        y = logcpm(counts[:, keep])
        feats = ps.var_names.to_numpy()[keep]
        symbols = (
            ps.var["gene_symbol"].to_numpy()[keep] if "gene_symbol" in ps.var.columns else None
        )
        groups = o[list(CONDITION)].astype(str).agg("|".join, axis=1).to_numpy()
        codes, uniq = pd.factorize(groups)
        n_g = np.bincount(codes)
        sums = np.zeros((len(uniq), y.shape[1]))
        np.add.at(sums, codes, y)
        means = sums / n_g[:, None]
        resid = y - means[codes]
        df_res = float(len(codes) - len(uniq))
        if df_res < 1:
            raise TinctureDataError(
                f"No replicates in stratum {stratum}: residual variance cannot be estimated.",
                fix="Signatures need replicate samples for at least some conditions.",
            )
        s2 = (resid**2).sum(axis=0) / df_res
        amean = y.mean(axis=0)
        sq = squeeze_var(s2, df_res, covariate=amean if trend else None)
        veh_code = pd.factorize(groups)[0][is_veh][0]
        veh_mean = y[is_veh].mean(axis=0)
        n_v = int(is_veh.sum())
        sigs: list[Signature] = []
        pvals: dict[tuple[str, float], np.ndarray] = {}
        for g in range(len(uniq)):
            if g == veh_code:
                continue
            rows = codes == g
            first = o.iloc[np.flatnonzero(rows)[0]]
            n_c = int(rows.sum())
            if n_c < min_replicates:
                continue
            eff = means[g] - veh_mean
            se = np.sqrt(sq.var_post * (1.0 / n_c + 1.0 / n_v))
            t = eff / se
            p = 2 * stats.t.sf(np.abs(t), df=sq.df_total)
            if shrink:
                post, post_sd, tau2 = shrink_normal(eff, se)
            else:
                post, post_sd, tau2 = eff, se, float("nan")
            sig = Signature(
                features=feats,
                effect=post,
                se=post_sd,
                compound_id=str(first["compound_id"]),
                dose=float(first["dose"]),
                time=float(first["time"]) if "time" in o.columns else float("nan"),
                cell_line=str(first["cell_line"]) if "cell_line" in o.columns else "",
                n_replicates=n_c,
                method="logcpm_eb" + ("+shrink" if shrink else ""),
                meta={
                    "raw_effect": eff,
                    "raw_se": se,
                    "moderated_t": t,
                    "tau2": tau2,
                    "compound_name": first.get("compound_name", ""),
                    "n_vehicle": n_v,
                },
            )
            sigs.append(sig)
            pvals[(sig.compound_id, sig.dose)] = p
        if n_v < 2:
            warnings.warn(
                f"Stratum {stratum} has a single vehicle sample; effects are relative to one well.",
                TinctureWarning,
                stacklevel=2,
            )
        out.append(
            SignatureSet(
                signatures=sigs,
                stratum=stratum,
                features=feats,
                df_prior=sq.df_prior,
                df_residual=df_res,
                n_vehicle=n_v,
                method="logcpm_eb" + ("+shrink" if shrink else ""),
                pvalues=pvals,
                feature_symbols=symbols,
            )
        )
    return out


def replicate_profiles(
    ps: PerturbationSet,
    strata: Sequence[str] = DEFAULT_STRATA,
    min_cpm: float = 5.0,
    min_frac: float = 0.25,
    n_features: int | None = 2000,
) -> PerturbationSet:
    """Per-sample profiles for retrieval: logCPM minus same-stratum vehicle mean, scaled
    by the vehicle SD per feature (robust z vs vehicle).

    Feature selection is **label-free** (highest variance across all samples of the
    stratum after the expression filter), so retrieval is not inflated by choosing
    features that separate the very conditions being retrieved.
    """
    obs = ps.obs
    strata = [s for s in strata if s in obs.columns]
    keys = (
        obs[strata].astype(str).agg("|".join, axis=1)
        if strata
        else pd.Series("all", index=obs.index)
    )
    blocks: list[pd.DataFrame] = []
    for key in sorted(keys.unique()):
        idx = np.flatnonzero((keys == key).to_numpy())
        o = obs.iloc[idx]
        is_veh = (o["dose"] == 0).to_numpy()
        if is_veh.sum() < 2:
            raise TinctureDataError(
                f"Stratum {key!r} has {int(is_veh.sum())} vehicle samples; need ≥ 2 to scale "
                "profiles.",
                fix="Add vehicle replicates or pool strata.",
            )
        counts = _dense(ps.X[idx, :]) if sp.issparse(ps.X) else _dense(ps.X[idx])
        cpm = counts / counts.sum(axis=1, keepdims=True) * 1e6
        keep = (cpm >= min_cpm).mean(axis=0) >= min_frac
        y = logcpm(counts[:, keep])
        mu = y[is_veh].mean(axis=0)
        sd = y[is_veh].std(axis=0, ddof=1)
        sd = np.maximum(sd, np.median(sd[sd > 0]) if (sd > 0).any() else 1.0)
        z = (y - mu) / sd
        feats = ps.var_names.to_numpy()[keep]
        if n_features is not None and z.shape[1] > n_features:
            top = np.argsort(-z.var(axis=0))[:n_features]
            z, feats = z[:, top], feats[top]
        blocks.append(pd.DataFrame(z, index=ps.obs_names[idx], columns=feats))
    common = (
        sorted(set.intersection(*(set(b.columns) for b in blocks)))
        if len(blocks) > 1
        else list(blocks[0].columns)
    )
    import anndata as ad

    Z = pd.concat([b[common] for b in blocks]).loc[ps.obs_names]
    out = ad.AnnData(X=Z.to_numpy(), obs=ps.obs.copy(), var=pd.DataFrame(index=pd.Index(common)))
    return PerturbationSet(out, source=ps.source)
