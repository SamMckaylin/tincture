"""Plate quality: Z′, edge effect and control drift.

* **Z′** (Zhang, Chung & Oldenburg, J Biomol Screen 1999):
  ``Z′ = 1 − 3(σ_pos + σ_neg) / |μ_pos − μ_neg|``. A robust variant replaces mean/SD with
  median/MAD×1.4826. Plates with Z′ < 0.5 are flagged as failing.
* **Edge effect**: median(edge) / median(interior) over drug-free wells that contain cells
  (signal > ``seeded_frac`` × the plate's negative-control median). Many plate layouts leave
  the outermost ring unseeded, so "edge" is the **outermost seeded ring** (ring depth 0 =
  outermost; the first depth with ≥ ``min_controls`` seeded wells) and "interior" is depth
  ≥ edge depth + 2. The depth used and unseeded wells are reported.
* **Drift**: each plate's negative-control median relative to the median of comparable
  plates (same cell line when a ``cell_line`` column exists, otherwise all plates),
  against plate order; Spearman ρ and the Theil–Sen change across the observed plates.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from tincture.errors import TinctureDataError
from tincture.model.plate import PlateTable

Z_PRIME_PASS = 0.5
MAD_TO_SD = 1.4826


def zprime(pos: np.ndarray, neg: np.ndarray, robust: bool = False) -> float:
    """Z′ factor; NaN when either control group has fewer than 2 finite values."""
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if pos.size < 2 or neg.size < 2:
        return float("nan")
    if robust:
        mp, mn = np.median(pos), np.median(neg)
        sp_ = MAD_TO_SD * np.median(np.abs(pos - mp))
        sn = MAD_TO_SD * np.median(np.abs(neg - mn))
    else:
        mp, mn = pos.mean(), neg.mean()
        sp_, sn = pos.std(ddof=1), neg.std(ddof=1)
    sep = abs(mp - mn)
    if sep == 0:
        return float("-inf")
    return float(1 - 3 * (sp_ + sn) / sep)


@dataclass
class Selector:
    """Selects wells by ``well_type`` values and/or source ``tag`` values (OR)."""

    well_types: Sequence[str] = ()
    tags: Sequence[str] = ()

    def mask(self, d: pd.DataFrame) -> np.ndarray:
        m = d["well_type"].isin(list(self.well_types)).to_numpy().copy()
        if self.tags and "tag" in d.columns:
            m |= d["tag"].isin(list(self.tags)).to_numpy()
        return m

    def describe(self) -> str:
        parts = [*self.well_types, *[f"tag {t}" for t in self.tags]]
        return " + ".join(parts) if parts else "nothing"


@dataclass
class PlateQC:
    per_plate: pd.DataFrame  # plate, order, z_prime, z_prime_robust, n_pos, n_neg, edge_ratio, ...
    neg: Selector
    pos: Selector
    edge: Selector
    drift_rho: float
    drift_p: float
    drift_change: float  # Theil–Sen relative change in control level, first → last plate
    n_missing_wells: int
    threshold: float = Z_PRIME_PASS
    drift_basis: str = "all plates"
    n_unseeded_outer: int = 0

    @property
    def failing(self) -> pd.DataFrame:
        return self.per_plate[self.per_plate["fails"]]


def plate_qc(
    pt: PlateTable,
    neg: Selector | None = None,
    pos: Selector | None = None,
    edge: Selector | str = "auto",
    min_controls: int = 3,
    threshold: float = Z_PRIME_PASS,
    seeded_frac: float = 0.1,
) -> PlateQC:
    """Per-plate Z′, edge ratio and control level, plus drift across plate order.

    ``edge="auto"`` uses negative controls when every plate has ≥ ``min_controls`` of them
    both on the edge and inside; otherwise all drug-free cell wells (neg_ctrl + empty).
    """
    neg = neg or Selector(well_types=("neg_ctrl",))
    pos = pos or Selector(well_types=("pos_ctrl",))
    d = pt.data.copy()
    d["_edge"] = pt.is_edge()
    d["_neg"] = neg.mask(d)
    d["_pos"] = pos.mask(d)
    if not d["_neg"].any():
        raise TinctureDataError(
            f"No negative-control wells ({neg.describe()}).",
            fix="Plate QC needs negative controls; pass neg=Selector(...) to choose them.",
        )
    if not d["_pos"].any():
        raise TinctureDataError(
            f"No positive-control wells ({pos.describe()}).",
            fix="Pass pos=Selector(...) (e.g. well_types=('blank',) for a full-kill reference).",
        )
    if isinstance(edge, str):
        if edge != "auto":
            raise ValueError("edge must be 'auto' or a Selector.")
        cnt = d[d["_neg"]].groupby(["plate", "_edge"], observed=True).size().unstack(fill_value=0)
        ok = cnt.reindex(columns=[False, True], fill_value=0).min(axis=1) >= min_controls
        edge_sel = (
            neg
            if (len(ok) == d["plate"].nunique() and ok.all())
            else Selector(well_types=("neg_ctrl", "empty"))
        )
    else:
        edge_sel = edge
    d["_edgewell"] = edge_sel.mask(d)
    nr, nc = pt.shape
    r, c = d["row"].to_numpy(), d["col"].to_numpy()
    d["_depth"] = np.minimum(np.minimum(r, c), np.minimum(nr - 1 - r, nc - 1 - c))
    order = (
        d.drop_duplicates("plate").set_index("plate")["plate_order"]
        if "plate_order" in d.columns
        else pd.Series(np.arange(d["plate"].nunique()), index=sorted(d["plate"].unique()))
    )
    rows = []
    for plate, g in d.groupby("plate", sort=False, observed=True):
        pv, nv = g.loc[g["_pos"], "value"].to_numpy(), g.loc[g["_neg"], "value"].to_numpy()
        neg_med = float(np.median(nv)) if nv.size else np.nan
        ew = g[g["_edgewell"]]
        seeded = ew[ew["value"] > seeded_frac * neg_med] if np.isfinite(neg_med) else ew.iloc[:0]
        outer_unseeded = bool(
            ((ew["_depth"] == 0).sum() >= min_controls)
            and ((seeded["_depth"] == 0).sum() < min_controls)
        )
        depth_counts = seeded.groupby("_depth").size()
        edge_depth = next((int(k) for k, v in depth_counts.items() if v >= min_controls), None)
        ratio = np.nan
        if edge_depth is not None:
            e_out = seeded.loc[seeded["_depth"] == edge_depth, "value"]
            e_in = seeded.loc[seeded["_depth"] >= edge_depth + 2, "value"]
            if len(e_in) >= min_controls and e_in.median() > 0:
                ratio = float(e_out.median() / e_in.median())
        rows.append(
            {
                "plate": str(plate),
                "order": float(order.get(plate, np.nan)),
                "z_prime": zprime(pv, nv) if min(pv.size, nv.size) >= min_controls else np.nan,
                "z_prime_robust": zprime(pv, nv, robust=True)
                if min(pv.size, nv.size) >= min_controls
                else np.nan,
                "n_pos": int(pv.size),
                "n_neg": int(nv.size),
                "neg_median": float(np.median(nv)) if nv.size else np.nan,
                "pos_median": float(np.median(pv)) if pv.size else np.nan,
                "edge_ratio": ratio,
                "edge_depth": edge_depth if edge_depth is not None else np.nan,
                "outer_unseeded": outer_unseeded,
                "cell_line": str(g["cell_line"].iloc[0]) if "cell_line" in g.columns else "",
                "n_wells": len(g),
            }
        )
    t = pd.DataFrame(rows).sort_values("order").reset_index(drop=True)
    t["fails"] = t["z_prime"] < threshold
    if "cell_line" in pt.data.columns and t["cell_line"].nunique() > 1:
        ref = t.groupby("cell_line")["neg_median"].transform("median")
        basis = "same cell line"
    else:
        ref = t["neg_median"].median()
        basis = "all plates"
    lvl = t["neg_median"] / ref
    ok = lvl.notna() & t["order"].notna() & (lvl > 0)
    if ok.sum() >= 3:
        rho, p = stats.spearmanr(t.loc[ok, "order"], lvl[ok])
        xo = t.loc[ok, "order"].to_numpy()
        slope = stats.theilslopes(np.log(lvl[ok].to_numpy()), xo)[0]
        change = float(np.exp(slope * (xo.max() - xo.min())) - 1)
    else:
        rho, p, change = np.nan, np.nan, np.nan
    t["neg_level"] = lvl
    missing = int(len(t) * nr * nc - t["n_wells"].sum())
    return PlateQC(
        per_plate=t,
        neg=neg,
        pos=pos,
        edge=edge_sel,
        drift_rho=float(rho),
        drift_p=float(p),
        drift_change=change,
        n_missing_wells=missing,
        threshold=threshold,
        drift_basis=basis,
        n_unseeded_outer=int(t["outer_unseeded"].sum()),
    )
