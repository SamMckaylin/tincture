"""PerturbationSet: samples × features with perturbation metadata, wrapping AnnData."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from tincture.errors import TinctureDataError
from tincture.model.measurements import validate_long

#: obs columns every PerturbationSet must have.
OBS_REQUIRED: tuple[str, ...] = ("compound_id", "dose", "replicate")
#: obs columns tincture understands when present.
OBS_KNOWN: tuple[str, ...] = (
    "compound_name",
    "smiles",
    "time",
    "cell_line",
    "plate",
    "well",
    "well_type",
    "n_cells",
    "target",
    "pathway",
)


class PerturbationSet:
    """Samples (cells, wells or pseudobulk groups) × features with perturbation metadata.

    ``obs`` holds one row per sample with at least ``compound_id``, ``dose`` (molar; 0 =
    vehicle) and ``replicate``. ``var`` holds features (genes, pathways, morphology
    features). ``X`` may be dense or sparse. Additional AnnData layers are preserved.

    Signatures (``.signatures()``) are added in Milestone 2 (``tincture.stats``).
    """

    def __init__(self, adata: ad.AnnData, source: str = "") -> None:
        missing = [c for c in OBS_REQUIRED if c not in adata.obs.columns]
        if missing:
            raise TinctureDataError(
                f"AnnData.obs is missing {missing}.",
                fix="Add obs columns compound_id, dose (molar; 0 = vehicle) and replicate.",
            )
        obs = adata.obs
        dose = pd.to_numeric(obs["dose"], errors="coerce")
        if dose.isna().any():
            raise TinctureDataError(
                f"{int(dose.isna().sum())} sample(s) have no numeric dose.",
                fix="Fill obs['dose'] in molar (0 for vehicle) or subset those samples out.",
            )
        if (dose < 0).any() or (dose > 1).any():
            raise TinctureDataError(
                "obs['dose'] must be molar in [0, 1].",
                fix="Convert to molar (µM × 1e-6, nM × 1e-9).",
            )
        adata.obs["dose"] = dose.astype(float)
        for c in ("compound_id", "replicate"):
            adata.obs[c] = adata.obs[c].astype(str)
        adata.obs_names_make_unique()
        self.adata = adata
        self.source = source or str(adata.uns.get("tincture_source", ""))
        self.adata.uns["tincture_source"] = self.source

    # ------------------------------------------------------------------ basics
    @property
    def obs(self) -> pd.DataFrame:
        return cast(pd.DataFrame, self.adata.obs)

    @property
    def var(self) -> pd.DataFrame:
        return cast(pd.DataFrame, self.adata.var)

    @property
    def obs_names(self) -> pd.Index:
        return self.adata.obs_names

    @property
    def var_names(self) -> pd.Index:
        """Feature names (e.g. Ensembl gene ids)."""
        return self.adata.var_names

    @property
    def X(self) -> Any:
        return self.adata.X

    @property
    def shape(self) -> tuple[int, int]:
        return (int(self.adata.n_obs), int(self.adata.n_vars))

    @property
    def compounds(self) -> list[str]:
        """Treated compound ids (vehicle excluded), sorted."""
        o = self.obs
        return sorted(o.loc[o["dose"] > 0, "compound_id"].unique())

    @property
    def is_vehicle(self) -> np.ndarray:
        return (self.obs["dose"] == 0).to_numpy()

    def doses(self, compound_id: str) -> np.ndarray:
        o = self.obs
        return np.sort(o.loc[o["compound_id"] == compound_id, "dose"].unique())

    def __len__(self) -> int:
        return self.adata.n_obs

    def __repr__(self) -> str:
        o = self.obs
        parts = [
            f"{self.adata.n_obs} samples × {self.adata.n_vars} features",
            f"{len(self.compounds)} compounds",
            f"{int(self.is_vehicle.sum())} vehicle samples",
        ]
        if "cell_line" in o:
            parts.append(f"cell lines: {sorted(o['cell_line'].astype(str).unique())[:5]}")
        return f"PerturbationSet({', '.join(parts)}; source={self.source!r})"

    # ------------------------------------------------------------------ selection
    def subset(
        self,
        obs: np.ndarray | pd.Series | None = None,
        features: Sequence[str] | np.ndarray | None = None,
        **equals: Any,
    ) -> PerturbationSet:
        """Subset samples and/or features.

        ``subset(cell_line="A549", time=24)`` keeps samples whose obs columns equal the
        given values (lists mean "any of"). ``obs`` is an optional boolean mask;
        ``features`` a list of feature names or a boolean mask over ``var``.
        """
        mask = np.ones(self.adata.n_obs, dtype=bool) if obs is None else np.asarray(obs, bool)
        for col, val in equals.items():
            if col not in self.obs.columns:
                raise TinctureDataError(
                    f"Unknown obs column {col!r}.", fix=f"Use one of {list(self.obs.columns)}."
                )
            vals = val if isinstance(val, (list, tuple, set, np.ndarray)) else [val]
            mask &= self.obs[col].isin(list(vals)).to_numpy()
        if not mask.any():
            raise TinctureDataError(
                f"Subset {equals or 'mask'} selects no samples.",
                fix="Check the values against PerturbationSet.obs.",
            )
        ad_sub = self.adata[mask]
        if features is not None:
            f = np.asarray(features)
            if f.dtype == bool:
                ad_sub = ad_sub[:, f]
            else:
                unknown = sorted(set(f) - set(self.adata.var_names))
                if unknown:
                    raise TinctureDataError(
                        f"{len(unknown)} unknown feature(s), e.g. {unknown[:3]}.",
                        fix="Check names against PerturbationSet.var_names.",
                    )
                ad_sub = ad_sub[:, list(f)]
        return PerturbationSet(ad_sub.copy(), source=self.source)

    # ------------------------------------------------------------------ statistics
    def signatures(self, **kwargs: Any) -> list[Any]:
        """Signatures vs same-stratum vehicle.

        See :func:`tincture.stats.signatures.compute_signatures`.
        """
        from tincture.stats.signatures import compute_signatures

        return compute_signatures(self, **kwargs)

    # ------------------------------------------------------------------ aggregation
    def pseudobulk(self, by: Sequence[str], min_samples: int = 1) -> PerturbationSet:
        """Sum ``X`` over samples sharing the values of ``by`` (e.g. well or replicate).

        ``by`` must include ``compound_id``, ``dose`` and ``replicate`` (or columns that
        determine them) so the result is still a valid PerturbationSet. The result has an
        ``n_cells`` obs column (summing any existing ``n_cells``).
        """
        by = list(by)
        missing = [c for c in by if c not in self.obs.columns]
        if missing:
            raise TinctureDataError(f"pseudobulk: unknown obs column(s) {missing}.")
        keys = self.obs[by].astype(str).agg("\x1f".join, axis=1)
        codes, uniques = pd.factorize(keys, sort=True)
        n_groups = len(uniques)
        ind = sp.csr_matrix(
            (np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(n_groups, len(codes))
        )
        X = self.adata.X
        summed = ind @ X
        summed = summed if sp.issparse(summed) else np.asarray(summed)
        first = pd.Series(np.arange(len(codes))).groupby(codes).first().to_numpy()
        obs = self.obs.iloc[first].copy()
        base = self.obs["n_cells"].to_numpy() if "n_cells" in self.obs else np.ones(len(codes))
        obs["n_cells"] = np.bincount(codes, weights=base, minlength=n_groups)
        obs.index = pd.Index([str(u).replace("\x1f", "|") for u in uniques])
        keep = obs["n_cells"].to_numpy() >= min_samples
        out = ad.AnnData(X=summed, obs=obs, var=self.var.copy())
        return PerturbationSet(out[keep].copy(), source=self.source)

    # ------------------------------------------------------------------ long format
    @classmethod
    def from_long(cls, df: pd.DataFrame, source: str = "") -> PerturbationSet:
        """Build from the canonical long table (one row per sample × feature)."""
        d = validate_long(df)
        sample_cols = [
            c
            for c in ("compound_id", "dose", "replicate", "time", "cell_line", "plate", "well")
            if c in d.columns
        ]
        d["_sample"] = d[sample_cols].astype(str).agg("|".join, axis=1)
        if d.duplicated(["_sample", "feature"]).any():
            raise TinctureDataError(
                "Several rows share the same sample and feature.",
                fix=f"Make samples unique via {sample_cols} (e.g. add a replicate or well column).",
            )
        wide = d.pivot(index="_sample", columns="feature", values="value")
        meta_cols = [c for c in d.columns if c not in ("feature", "value", "se", "_sample")]
        obs = d.drop_duplicates("_sample").set_index("_sample")[meta_cols].loc[wide.index]
        for c in obs.columns:
            if pd.api.types.is_string_dtype(obs[c]):
                obs[c] = obs[c].astype(object)
        adata = ad.AnnData(
            X=wide.to_numpy(dtype=float),
            obs=obs,
            var=pd.DataFrame(index=wide.columns.astype(str)),
        )
        return cls(adata, source=source)

    def to_long(self) -> pd.DataFrame:
        """Inverse of :meth:`from_long` (dense; intended for small sets and tests)."""
        X: Any = self.adata.X
        X = X.toarray() if sp.issparse(X) else np.asarray(X)
        wide = pd.DataFrame(X, index=self.adata.obs_names, columns=self.adata.var_names)
        long = wide.rename_axis(index="_sample", columns="feature").stack(future_stack=True)
        long = long.rename("value").reset_index()
        meta = self.obs.copy()
        meta.index.name = "_sample"
        out = long.merge(meta.reset_index(), on="_sample").drop(columns="_sample")
        return out
