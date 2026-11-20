"""Effect Ladder: which features moved, and how confident are we? (HANDOVER §5.4)

Replaces the volcano. A forest plot of the top N features by |shrunk log fold change|,
grouped by gene set (curated order), with 95% intervals: intervals that exclude zero are
coloured (terracotta up, blue down), intervals that cross zero are grey. A **count strip**
on top keeps the volcano's overview: how many features moved up/down significantly across
the whole transcriptome.

Guards (each tested): N and the shrinkage method are printed; the count strip is always
drawn. ``pathway=`` drills down into one gene set (from the Fingerprint Strip).
Classic fallback: volcano.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from tincture.errors import TinctureDataError
from tincture.model.signature import Signature
from tincture.plots._pathways import bare, curated_order, feature_sets, pathway_label
from tincture.report.captions import Caption, fmt_conc
from tincture.report.classic import build_classic, register_classic
from tincture.report.methods import Methods, register_block
from tincture.result import PlotResult
from tincture.spec import (
    IntervalMark,
    NoticeMark,
    Panel,
    PlotSpec,
    PointMark,
    RectMark,
    RefLineMark,
    Scale,
    TextMark,
)
from tincture.stats.pathways import symbols_for
from tincture.stats.signatures import SignatureSet

TOP_N = 30
GROUP_GAP = 1.0

SHRINKAGE = {
    "logcpm_eb+shrink": "limma-style empirical-Bayes variances with normal-prior LFC shrinkage",
    "logcpm_eb": "limma-style empirical-Bayes variances, no LFC shrinkage",
}

register_block(
    "effect_ladder",
    "The Effect Ladder shows the {n} features with the largest |log2 fold change| "
    "({shrink}){scope}, grouped by MSigDB Hallmark gene set (each gene under the first set "
    "containing it in the curated order); bars are 95% intervals (effect ± 1.96 SE). The "
    "count strip counts features with BH-adjusted q < {alpha:g} across all {total} "
    "features.",
)


def _resolve(
    data: Signature | SignatureSet, compound: str | None, dose: float | None
) -> tuple[Signature, np.ndarray | None, np.ndarray | None]:
    """(signature, p-values, symbols) from a Signature or a SignatureSet."""
    if isinstance(data, Signature):
        return data, data.meta.get("pvalues"), data.meta.get("feature_symbols")
    if compound is None:
        raise TinctureDataError("Pass compound= to pick a signature from the SignatureSet.")
    cands = [s for s in data.signatures if s.compound_id == compound]
    if not cands:
        raise TinctureDataError(f"No signature for {compound!r}.", fix="Check the compound id.")
    sig = max(cands, key=lambda s: s.dose) if dose is None else data.get(compound, dose)
    p = data.pvalues.get((sig.compound_id, sig.dose))
    return sig, p, data.feature_symbols


def _pvalues(sig: Signature, p: np.ndarray | None) -> np.ndarray:
    if p is not None:
        return np.asarray(p, float)
    if sig.se is None:
        raise TinctureDataError(
            "The Effect Ladder needs standard errors or p-values.",
            fix="Use signatures from tc.signatures(...) or pass a table with SEs.",
        )
    from scipy import stats

    z = sig.effect / sig.se
    return 2 * stats.norm.sf(np.abs(z))


def compute(
    sig: Signature,
    pvalues: np.ndarray | None,
    symbols: np.ndarray | None,
    gene_sets: pd.DataFrame | None,
    top_n: int = TOP_N,
    pathway: str | None = None,
    alpha: float = 0.05,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Table of the chosen features (in drawing order) and transcriptome-wide counts."""
    sym = symbols_for(sig.features, symbols)
    if sig.se is None:
        raise TinctureDataError("The Effect Ladder needs per-feature standard errors.")
    p = _pvalues(sig, pvalues)
    ok = np.isfinite(sig.effect) & np.isfinite(sig.se) & np.isfinite(p)
    q = np.full(p.shape, np.nan)
    q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    t = pd.DataFrame(
        {"feature": sig.features, "symbol": sym, "lfc": sig.effect, "se": sig.se, "p": p, "q": q}
    )[ok]
    counts = {
        "up": int(((t.q < alpha) & (t.lfc > 0)).sum()),
        "down": int(((t.q < alpha) & (t.lfc < 0)).sum()),
        "total": len(t),
    }
    sets: dict[str, str] = {}
    if gene_sets is not None:
        sets = feature_sets(gene_sets, t["symbol"].to_numpy())
    if pathway is not None:
        if gene_sets is None:
            raise TinctureDataError("pathway= needs gene_sets.")
        src = gene_sets["source"].astype(str).map(bare)
        members = set(gene_sets.loc[src == bare(pathway), "target"].astype(str))
        if not members:
            raise TinctureDataError(
                f"Gene set {pathway!r} not found.", fix="Use a name from the Fingerprint Strip."
            )
        t = t[t["symbol"].isin(members)]
        sets = {g: bare(pathway) for g in t["symbol"]}
        if t.empty:
            raise TinctureDataError(f"No members of {pathway!r} among the features.")
    t = t.assign(abs_lfc=t.lfc.abs()).nlargest(top_n, "abs_lfc")
    t["gene_set"] = [sets.get(g, "") for g in t["symbol"]]
    order = {s: k for k, s in enumerate(curated_order(t["gene_set"][t["gene_set"] != ""].unique()))}
    t["g_rank"] = t["gene_set"].map(lambda s: order.get(s, len(order)))
    t = t.sort_values(["g_rank", "lfc"], ascending=[True, False]).reset_index(drop=True)
    lo, hi = t.lfc - 1.96 * t.se, t.lfc + 1.96 * t.se
    t["call"] = np.where(lo > 0, "up", np.where(hi < 0, "down", "ns"))
    t["lo"], t["hi"] = lo, hi
    return t, counts


def _positions(groups: Sequence[str]) -> np.ndarray:
    pos, y = [], 0.0
    for i, g in enumerate(groups):
        if i and g != groups[i - 1]:
            y += GROUP_GAP
        pos.append(y)
        y += 1
    return np.asarray(pos)


def build_spec(
    t: pd.DataFrame,
    counts: dict[str, int],
    name: str,
    shrink: str,
    pathway: str | None,
    alpha: float,
    title: str,
    footer: Sequence[str],
    dose: float | None = None,
) -> PlotSpec:
    y = _positions(list(t["gene_set"]))
    colors = [{"up": "up", "down": "down", "ns": "mid_grey"}[c] for c in t["call"]]
    span = float(max(np.abs(t["lo"]).max(), np.abs(t["hi"]).max(), 0.1))
    xdom = (-1.12 * span, 1.12 * span)
    # count strip --------------------------------------------------------------------------
    strip = Panel(
        key="counts",
        x=Scale("linear", domain=(-1.0, 1.0), visible=False),
        y=Scale("linear", domain=(-0.6, 1.0), visible=False),
        frame=False,
    )
    tot = max(counts["up"] + counts["down"], 1)
    share_dn, share_up = counts["down"] / tot, counts["up"] / tot
    strip.add(
        RectMark(x0=[-share_dn], y0=[0], width=[share_dn], height=[0.45], color="down"),
        RectMark(x0=[0], y0=[0], width=[share_up], height=[0.45], color="up"),
        RefLineMark(axis="x", value=0, color="ink_soft", style="solid", width=0.6),
        TextMark(x=-1.0, y=0.75, text=f"{counts['down']:,} down", color="down", ha="left",
                 va="center", categorical=True),
        TextMark(x=1.0, y=0.75, text=f"{counts['up']:,} up", color="up", ha="right",
                 va="center", categorical=True),
        TextMark(x=0, y=0.75, text=f"of {counts['total']:,} features, q < {alpha:g}",
                 color="ink_soft", ha="center", va="center"),
    )  # fmt: skip
    if counts["up"] + counts["down"] == 0:
        strip.add(NoticeMark(text="no feature passes q < 0.05", loc="center", level="weak"))
    # forest -------------------------------------------------------------------------------
    f = Panel(
        key="forest",
        x=Scale("linear", label="log2 fold change (shrunk) with 95% interval", domain=xdom),
        y=Scale(
            "linear",
            ticks=list(y),
            tick_labels=list(t["symbol"]),
            domain=(float(y.max()) + 0.7 if len(y) else 1, -1.6),
        ),
    )
    f.add(RefLineMark(axis="x", value=0.0, style="dashed", color="ink_soft"))
    f.add(IntervalMark(pos=y, lo=t["lo"], hi=t["hi"], orientation="h", color=colors, width=1.4))
    f.add(PointMark(x=t["lfc"], y=y, color=colors, size=16, z=3))
    starts = [i for i in range(len(t)) if i == 0 or t["gene_set"][i] != t["gene_set"][i - 1]]
    for i in starts:
        g = t["gene_set"][i]
        f.add(
            TextMark(
                x=xdom[0] + 0.015 * span,
                y=float(y[i]) - 0.75,
                text=pathway_label(g) if g else "No Hallmark set",
                color="ink_soft",
                ha="left",
                va="center",
                weight="bold",
                size=6.5,
            )
        )
    n_ci0 = int((t["call"] == "ns").sum())
    cap = Caption()
    what = name + (f" at {fmt_conc(dose)}" if dose is not None and np.isfinite(dose) else "")
    cap.add(
        "Top {n} features of {what} by |shrunk log2 fold change|{scope}; {k} have 95% "
        "intervals excluding zero. Across all {total:,} features, {up:,} go up and {dn:,} go "
        "down (q < {a:g}).",
        n=len(t),
        what=what,
        scope=f" in {pathway_label(pathway)}" if pathway else "",
        k=len(t) - n_ci0,
        total=counts["total"],
        up=counts["up"],
        dn=counts["down"],
        a=alpha,
    )
    if n_ci0:
        cap.weak(f"{n_ci0} of the top {len(t)} have intervals crossing zero (grey).")
    meth = Methods().use(
        "effect_ladder",
        n=len(t),
        shrink=shrink,
        scope=f" within {pathway_label(pathway)}" if pathway else "",
        alpha=alpha,
        total=f"{counts['total']:,}",
    )
    return PlotSpec(
        plot="effect_ladder",
        panels=[strip, f],
        nrows=2,
        ncols=1,
        height_ratios=[0.09, 1.0],
        title=title,
        caption=cap.render(),
        methods=meth.render(),
        footer=[f"N = {len(t)} features; shrinkage: {shrink}", *footer],
        warnings=list(cap.notes),
        aspect=0.42 + 0.032 * len(t),
        meta={"table": t, "counts": counts, "pathway": pathway},
    ).validate()


def volcano_spec(
    sig: Signature, pvalues: np.ndarray | None, symbols: np.ndarray | None, alpha: float = 0.05
) -> PlotSpec:
    """Classic volcano: shrunk log2 FC vs −log10 p, q < alpha coloured."""
    p = _pvalues(sig, pvalues)
    ok = np.isfinite(sig.effect) & np.isfinite(p)
    q = np.ones_like(p)
    q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    x, yv = sig.effect[ok], -np.log10(np.clip(p[ok], 1e-300, 1))
    sig_up = (q[ok] < alpha) & (x > 0)
    sig_dn = (q[ok] < alpha) & (x < 0)
    rest = ~(sig_up | sig_dn)
    pn = Panel(x=Scale(label="log2 fold change"), y=Scale(label="−log10 p"))
    pn.add(PointMark(x=x[rest], y=yv[rest], color="neutral", size=4))
    pn.add(PointMark(x=x[sig_dn], y=yv[sig_dn], color="down", size=5, categorical=True))
    pn.add(PointMark(x=x[sig_up], y=yv[sig_up], color="up", size=5, categorical=True))
    sym = symbols_for(sig.features, symbols)[ok]
    top = np.argsort(-yv)[:8]
    for i in top:
        pn.add(TextMark(x=float(x[i]), y=float(yv[i]), text=f" {sym[i]}", size=6))
    return PlotSpec(
        plot="volcano",
        panels=[pn],
        title="Volcano",
        caption=f"Volcano: {int(sig_up.sum())} up, {int(sig_dn.sum())} down at q < {alpha:g}.",
        aspect=0.8,
    ).validate()


@register_classic("effect_ladder", "volcano")
def _classic(
    sig: Signature | None = None,
    pvalues: np.ndarray | None = None,
    symbols: np.ndarray | None = None,
    alpha: float = 0.05,
    **_: object,
) -> PlotSpec:
    if sig is None:
        raise TinctureDataError("Volcano needs a signature.")
    return volcano_spec(sig, pvalues, symbols, alpha)


def plot(
    data: Signature | SignatureSet,
    compound: str | None = None,
    dose: float | None = None,
    gene_sets: pd.DataFrame | None = None,
    pathway: str | None = None,
    top_n: int = TOP_N,
    alpha: float = 0.05,
    name: str | None = None,
    title: str = "Effect Ladder",
    footer: Sequence[str] = (),
) -> PlotResult:
    """Effect Ladder for one signature (a :class:`Signature`, or ``compound``/``dose`` of a
    :class:`SignatureSet`). ``gene_sets``: long ``source, target`` table (e.g.
    ``load_geneset("hallmark").table``) for grouping and ``pathway=`` drill-down."""
    from tincture.plots.compass import short_name

    if top_n < 1:
        raise ValueError("top_n must be ≥ 1.")
    sig, p, symbols = _resolve(data, compound, dose)
    t, counts = compute(sig, p, symbols, gene_sets, top_n, pathway, alpha)
    method = sig.method or (data.method if isinstance(data, SignatureSet) else "")
    shrink = SHRINKAGE.get(method, method or "as supplied (method not recorded)")
    nm = name or short_name(str(sig.meta.get("compound_name") or sig.compound_id or "signature"))
    spec = build_spec(t, counts, nm, shrink, pathway, alpha, title, footer, sig.dose)
    classic = build_classic("effect_ladder", sig=sig, pvalues=p, symbols=symbols, alpha=alpha)
    return PlotResult(spec=spec, classic=classic)
