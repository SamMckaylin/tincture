"""Retrievability: which compounds have a real, reproducible phenotype? (HANDOVER §5.5)

This is the **gate** before every other plot. For each compound × dose, replicate profiles
should retrieve each other against vehicle profiles (mAP with a permutation null;
:mod:`tincture.stats.map_retrieval`). Compounds that fail are greyed here and can be
greyed or dropped by other plots (``gate=True``).

Panels: (A) strip of conditions sorted by mAP, with the per-condition null 5–95% band,
the q = 0.05 boundary and focus compounds labelled; (B, optional) CAT curve for a focus
compound's two replicates against the hypergeometric null band.
Classic fallback: replicate-vs-replicate scatter with Pearson r.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from tincture.errors import TinctureDataError
from tincture.model.perturbation import PerturbationSet
from tincture.report.captions import Caption, fmt_conc
from tincture.report.classic import build_classic, register_classic
from tincture.report.methods import Methods, register_block
from tincture.result import PlotResult
from tincture.spec import (
    BandMark,
    LegendSpec,
    LineMark,
    NoticeMark,
    Panel,
    PlotSpec,
    PointMark,
    RefLineMark,
    Scale,
    TextMark,
)
from tincture.stats.map_retrieval import CatCurve, MapResult, cat_curve, retrievability

GATE_KEY = "tincture_retrievable"

register_block(
    "map_retrieval",
    "Replicate retrievability was scored as mean average precision (mAP): for each "
    "replicate profile, the other replicates of the same {unit} were ranked among "
    "{negatives} by cosine similarity, average precision was computed per profile and "
    "averaged per {unit}. P values come from {null_desc}, corrected with the "
    "Benjamini–Hochberg procedure; conditions with q < {alpha} were called retrievable.",
    (
        "Kalinin AA, et al. A versatile information retrieval framework for evaluating "
        "profile strength and similarity. bioRxiv 2024 (copairs).",
    ),
)
register_block(
    "cat_curve",
    "Correspondence at the top (CAT): fraction of the top-k features (by |z|) shared by "
    "two replicate profiles, compared with the exact hypergeometric {level:.0%} band for "
    "two random k-subsets.",
)


def _group_keys(obs: pd.DataFrame, group_by: Sequence[str]) -> pd.Series:
    return obs[list(group_by)].astype(str).agg("|".join, axis=1)


def _label(obs: pd.DataFrame, key: str, group_by: Sequence[str]) -> str:
    keys = _group_keys(obs, group_by)
    row = obs[keys == key].iloc[0]
    name = str(row.get("compound_name", row["compound_id"]))
    return f"{name} {fmt_conc(float(row['dose']))}" if "dose" in group_by else name


def compute(
    profiles: PerturbationSet,
    group_by: Sequence[str] = ("compound_id", "dose"),
    negatives: str = "controls",
    null_size: int = 10_000,
    alpha: float = 0.05,
    seed: int = 0,
    null: str = "permutation",
) -> MapResult:
    """Retrievability statistics for replicate-level profiles (controls = dose 0)."""
    o = profiles.obs
    missing = [c for c in group_by if c not in o.columns]
    if missing:
        raise TinctureDataError(f"group_by columns {missing} not in obs.")
    X = profiles.X
    X = np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=float)
    return retrievability(
        X,
        _group_keys(o, group_by).to_numpy(),
        profiles.is_vehicle,
        negatives,
        null_size,
        alpha,
        seed,
        null,
    )


def gate_mask(
    profiles: PerturbationSet, result: MapResult, group_by: Sequence[str] = ("compound_id", "dose")
) -> pd.Series:
    """Boolean per sample: True if its condition is retrievable (vehicle always True).

    Also stored in ``profiles.adata.uns['tincture_retrievable']`` for other plots' ``gate``.
    """
    keys = _group_keys(profiles.obs, group_by)
    ok = set(result.table.loc[result.table["retrievable"], "group"])
    mask = keys.isin(ok) | pd.Series(profiles.is_vehicle, index=keys.index)
    profiles.adata.uns[GATE_KEY] = {
        "retrievable_groups": sorted(ok),
        "group_by": list(group_by),
        "alpha": result.alpha,
        "null_size": result.null_size,
    }
    return mask


def _q_boundary(t: pd.DataFrame) -> float | None:
    """mAP value separating retrievable from not, if the split is monotone in mAP."""
    if t["retrievable"].all() or not t["retrievable"].any():
        return None
    lo_sig = t.loc[t["retrievable"], "map"].min()
    hi_non = t.loc[~t["retrievable"], "map"].max()
    return float((lo_sig + hi_non) / 2) if lo_sig > hi_non else None


def build_spec(
    result: MapResult,
    labels: dict[str, str] | None = None,
    focus: Sequence[str] = (),
    cat: CatCurve | None = None,
    cat_label: str = "",
    title: str = "Retrievability",
    footer: Sequence[str] = (),
) -> PlotSpec:
    t = result.table.copy()
    if t.empty:
        raise TinctureDataError("No conditions to plot.")
    labels = labels or {}
    t["rank"] = np.arange(1, len(t) + 1)
    focus_set = set(focus)
    unknown = focus_set - set(t["group"])
    if unknown:
        raise TinctureDataError(
            f"Focus condition(s) {sorted(unknown)[:3]} have no retrievability score.",
            fix="Use group keys from result.table['group'] (e.g. 'S1090|1e-06').",
        )

    strip = Panel(
        key="strip",
        x=Scale(
            "linear",
            label=f"Conditions ranked by mAP (n = {len(t)})",
            domain=(0, len(t) + 1),
            ticks=[float(i) for i in range(1, len(t) + 1)] if len(t) <= 10 else None,
        ),
        y=Scale("linear", label="mAP (replicates vs vehicle)", domain=(0, 1.04)),
        legend=LegendSpec(loc="lower left"),
        title="Is the phenotype reproducible?",
    )
    strip.add(
        BandMark(
            lo=float(t["null_lo"].median()),
            hi=float(t["null_hi"].median()),
            span=True,
            orientation="h",
            color="neutral",
            alpha=0.6,
            legend="null 5–95% (median over conditions)",
        ),
    )
    ret = t["retrievable"].to_numpy()
    fmask = t["group"].isin(focus_set).to_numpy()
    if (ret & ~fmask).any():
        strip.add(
            PointMark(
                x=t["rank"][ret & ~fmask],
                y=t["map"][ret & ~fmask],
                color="reference",
                size=8,
                categorical=True,
                legend=f"retrievable (q < {result.alpha:g})",
            )
        )
    if (~ret & ~fmask).any():
        strip.add(
            PointMark(
                x=t["rank"][~ret & ~fmask],
                y=t["map"][~ret & ~fmask],
                color="ink_soft",
                size=8,
                hollow=True,
                legend="not retrievable",
            )
        )
    if fmask.any():
        f = t[fmask]
        strip.add(
            PointMark(
                x=f["rank"],
                y=f["map"],
                color="focus",
                size=22,
                categorical=True,
                labels=[labels.get(g, g) for g in f["group"]],
                z=4,
                legend="focus",
            )
        )
    boundary = _q_boundary(t)
    if boundary is not None:
        strip.add(RefLineMark(axis="y", value=boundary, label=f"q = {result.alpha:g}"))

    n_ret, n = int(t["retrievable"].sum()), len(t)
    panels = [strip]
    cap = Caption()
    cap.add(
        "{n_ret} of {n} conditions ({pct:.0%}) have a reproducible phenotype (mAP, q < {alpha:g}; "
        "replicates ranked against {neg}); median mAP {med:.2f}.",
        n_ret=n_ret,
        n=n,
        pct=n_ret / n,
        alpha=result.alpha,
        neg="vehicle profiles" if result.negatives == "controls" else "all other profiles",
        med=float(t["map"].median()),
    )
    for g in focus:
        row = t.set_index("group").loc[g]
        cap.add(
            "{name}: mAP {m:.2f}, q {q} ({verdict}).",
            name=labels.get(g, g),
            m=float(row["map"]),
            q=("< 0.001" if row["q"] < 1e-3 else f"= {row['q']:.2g}"),
            verdict="retrievable"
            if row["retrievable"]
            else "not retrievable — interpret its other plots with caution",
        )
    if not result.design_can_pass:
        best = float(np.min(result.min_p)) if result.min_p is not None else float("nan")
        cap.weak(
            f"this design cannot reach q < {result.alpha:g}: even perfect retrieval gives p ≥ "
            f"{best:.2g} with {n} conditions tested, so 'not retrievable' here means "
            "'not testable', not 'no phenotype'."
        )
        strip.add(
            NoticeMark(text="design cannot reach significance", loc="upper right", level="warning")
        )
    if result.excluded:
        cap.weak(
            f"{len(result.excluded)} condition(s) with a single replicate could not be scored "
            "and are not shown."
        )
        strip.add(
            NoticeMark(text=f"{len(result.excluded)} unscored (1 replicate)", loc="lower right")
        )
    if boundary is None and 0 < n_ret < n:
        cap.weak(
            "the q-value boundary is not a single mAP value here (null depends on replicate "
            "counts); colours show per-condition calls."
        )

    meth = Methods().use(
        "map_retrieval",
        null_desc=(
            "relabelling which members of the replicate ∪ vehicle pool form the condition and "
            "recomputing mAP from the observed similarities (all relabelings enumerated when "
            "feasible)"
            if result.null == "permutation"
            else f"{result.null_size:,} random rankings with the same numbers of positives and "
            "candidates (the 'copied' null of copairs)"
        ),
        unit="condition (compound × dose)",
        negatives="vehicle profiles" if result.negatives == "controls" else "all other profiles",
        null_size=result.null_size,
        alpha=result.alpha,
    )
    ncols, width_ratios = 1, None
    if cat is not None:
        cp = Panel(
            key="cat",
            x=Scale("log", label="Top-k features"),
            y=Scale("linear", label="Fraction shared by replicates", domain=(0, 1.02)),
            legend=LegendSpec(loc="center left"),
            title="Replicate agreement (CAT)",
        )
        if cat_label:
            cp.add(
                TextMark(
                    x=0.03,
                    y=0.97,
                    coords="axes",
                    va="top",
                    text=cat_label,
                    color="focus",
                    categorical=True,
                )
            )
        cp.add(
            BandMark(
                x=cat.k,
                lo=cat.null_lo,
                hi=cat.null_hi,
                color="neutral",
                alpha=0.6,
                legend="random 95%",
            ),
            LineMark(
                x=cat.k,
                y=cat.observed,
                color="focus",
                width=1.4,
                categorical=True,
                legend="replicates",
            ),
        )
        above = float(np.mean(cat.observed > cat.null_hi))
        cap.add(
            "Replicate top-feature overlap{lbl} exceeds the random band at {frac:.0%} of k values.",
            lbl=f" for {cat_label}" if cat_label else "",
            frac=above,
        )
        if above == 0:
            cp.add(NoticeMark(text="no better than random", loc="upper left", level="weak"))
        meth.use("cat_curve", level=0.95)
        panels.append(cp)
        ncols, width_ratios = 2, [1.6, 1.0]
    meth.use("bh_fdr", alpha=result.alpha)
    feet = [
        f"mAP vs {'vehicle' if result.negatives == 'controls' else 'all other'} profiles; "
        f"{result.null} null; BH q < {result.alpha:g}",
        *footer,
    ]
    return PlotSpec(
        plot="retrievability",
        panels=panels,
        ncols=ncols,
        width_ratios=width_ratios,
        title=title,
        caption=cap.render(),
        methods=meth.render(),
        footer=feet,
        warnings=list(cap.notes),
        meta={
            "n_retrievable": n_ret,
            "n_conditions": n,
            "q_boundary": boundary,
            "table": t,
            "excluded": result.excluded,
            "null": result.null,
            "design_can_pass": result.design_can_pass,
        },
        aspect=0.42 if cat is not None else 0.55,
    ).validate()


@register_classic("retrievability", "replicate correlation scatter")
def build_classic_spec(
    rep_a: np.ndarray | None = None, rep_b: np.ndarray | None = None, label: str = "", **_: object
) -> PlotSpec:
    if rep_a is None or rep_b is None:
        raise TinctureDataError("Classic replicate scatter needs two replicate profiles.")
    a, b = np.asarray(rep_a, float), np.asarray(rep_b, float)
    r = float(np.corrcoef(a, b)[0, 1])
    lim = float(np.nanmax(np.abs(np.concatenate([a, b])))) * 1.05
    p = Panel(
        x=Scale(label="Replicate 1 (z vs vehicle)", domain=(-lim, lim)),
        y=Scale(label="Replicate 2 (z vs vehicle)", domain=(-lim, lim)),
        aspect="equal",
    )
    p.add(
        RefLineMark(axis="x", value=0, style="dotted"),
        RefLineMark(axis="y", value=0, style="dotted"),
        PointMark(x=a, y=b, color="ink_soft", size=4, alpha=0.6),
        TextMark(x=0.03, y=0.97, coords="axes", va="top", text=f"Pearson r = {r:.2f}"),
    )
    return PlotSpec(
        plot="retrievability_classic",
        panels=[p],
        title=f"Replicate concordance{': ' + label if label else ''}",
        caption=f"Replicate profiles of {label or 'the focus condition'}: Pearson r = {r:.2f}.",
        aspect=0.9,
    ).validate()


def plot(
    profiles: PerturbationSet,
    group_by: Sequence[str] = ("compound_id", "dose"),
    focus: Sequence[str] = (),
    negatives: str = "controls",
    null_size: int = 10_000,
    alpha: float = 0.05,
    seed: int = 0,
    title: str = "Retrievability",
    footer: Sequence[str] = (),
    null: str = "permutation",
) -> PlotResult:
    """Retrievability gate for replicate-level profiles.

    ``focus`` holds group keys (``"<compound_id>|<dose>"``) to highlight; the first one
    also gets the CAT panel and is used for the classic replicate scatter.
    """
    res = compute(profiles, group_by, negatives, null_size, alpha, seed, null)
    o = profiles.obs
    labels = {g: _label(o, g, group_by) for g in res.table["group"]}
    cat = None
    rep_a = rep_b = None
    cat_label = ""
    if focus:
        keys = _group_keys(o, group_by)
        idx = np.flatnonzero((keys == focus[0]).to_numpy())
        if idx.size >= 2:
            X = profiles.X
            X = np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=float)
            rep_a, rep_b = X[idx[0]], X[idx[1]]
            cat = cat_curve(rep_a, rep_b)
            cat_label = labels.get(focus[0], focus[0])
    spec = build_spec(res, labels, focus, cat, cat_label, title, footer)
    gate_mask(profiles, res, group_by)
    classic = (
        build_classic("retrievability", rep_a=rep_a, rep_b=rep_b, label=cat_label)
        if rep_a is not None
        else None
    )
    return PlotResult(spec=spec, classic=classic)
