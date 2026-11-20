"""Fingerprint Strip: what did each compound do, at pathway level? (HANDOVER §5.3)

Replaces volcano + GSEA bar charts. Rows = compounds (references first, separated by a
gap), columns = pathways in a fixed curated order grouped Cell cycle / Stress / Metabolism /
Signalling / Immune / Other. Cells are ULM pathway scores (t) on a 5-step diverging scale;
cells that are not significant (BH q ≥ alpha over the whole matrix) are drawn neutral; the
focus compound gets a terracotta tick on the left.

Guards (each tested): colour clip limits shown in the legend; row ordering method stated;
gene sets deduplicated (Jaccard > 0.5, keep the larger); drill-down via
``tc.effect_ladder(..., pathway=...)``. Classic fallback: volcano for one compound.
"""

from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from tincture.errors import TinctureDataError
from tincture.model.signature import Signature
from tincture.plots._pathways import (
    bare,
    curated_order,
    group_breaks,
    group_of,
    pathway_label,
)
from tincture.report.captions import TEMPLATES, Caption, fmt_conc
from tincture.report.classic import build_classic
from tincture.report.methods import Methods, register_block
from tincture.result import PlotResult
from tincture.spec import ColorbarSpec, HeatmapMark, Panel, PlotSpec, RectMark, Scale, TextMark
from tincture.stats.pathways import dedup_gene_sets, symbols_for, ulm
from tincture.stats.signatures import SignatureSet

MAX_ROW_LABELS = 60
NICE_CLIPS = (2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50)

register_block(
    "fingerprint",
    "Pathway activity per signature is the univariate linear model (ULM) t-score of the "
    "signature's log2 fold changes on gene-set membership ({source}; sets with ≥ 5 measured "
    "members; {n_dedup} near-duplicate sets with Jaccard > {jac:g} dropped, keeping the "
    "larger). P-values were BH-adjusted over all {cells:,} cells; cells with q ≥ {alpha:g} "
    "are drawn neutral. Colours use 5 steps clipped at ±{clip:g}. Rows: {ordering}.",
)


def _rows(
    data: SignatureSet | Sequence[Signature] | Mapping[str, Signature],
    dose: float | str | None,
) -> list[Signature]:
    if isinstance(data, SignatureSet):
        sigs = list(data.signatures)
    elif isinstance(data, Mapping):
        sigs = list(data.values())
    else:
        sigs = list(data)
    if dose == "top" or (dose is None and isinstance(data, SignatureSet)):
        top: dict[str, Signature] = {}
        for s in sigs:
            if s.compound_id not in top or s.dose > top[s.compound_id].dose:
                top[s.compound_id] = s
        sigs = list(top.values())
    elif isinstance(dose, int | float):
        sigs = [s for s in sigs if np.isclose(s.dose, float(dose))]
    if not sigs:
        raise TinctureDataError("No signatures to fingerprint.")
    return sigs


def _row_id(s: Signature, multi: bool) -> str:
    return f"{s.compound_id}@{s.dose:g}" if multi else s.compound_id


def nice_clip(scores: np.ndarray) -> float:
    """Symmetric colour limit: the 95th percentile of |score|, rounded up to a nice value."""
    a = np.abs(scores[np.isfinite(scores)])
    if a.size == 0:
        return NICE_CLIPS[0]
    q = float(np.quantile(a, 0.95))
    return float(next((c for c in NICE_CLIPS if c >= q), NICE_CLIPS[-1]))


def _order(score: pd.DataFrame, method: str) -> list[str]:
    if method == "given" or len(score) <= 2:
        return list(score.index)
    from scipy.cluster.hierarchy import leaves_list, linkage, optimal_leaf_ordering
    from scipy.spatial.distance import pdist

    X = np.nan_to_num(score.to_numpy())
    d = pdist(X, metric="correlation") if X.shape[1] > 1 else pdist(X)
    d = np.nan_to_num(d, nan=1.0)
    Z = linkage(d, method="average")
    if len(score) <= 400:
        Z = optimal_leaf_ordering(Z, d)
    return [score.index[i] for i in leaves_list(Z)]


def compute(
    sigs: list[Signature],
    gene_sets: pd.DataFrame,
    symbols: np.ndarray | None,
    max_jaccard: float = 0.5,
    min_size: int = 5,
) -> tuple[pd.DataFrame, pd.DataFrame, list[tuple[str, str, float]]]:
    """ULM scores, BH q-values (rows × pathways) and the dedup decisions."""
    multi = len({s.compound_id for s in sigs}) < len(sigs)
    feats = sigs[0].features
    sym = symbols_for(feats, symbols)
    Y = pd.DataFrame(
        np.vstack([s.align(feats).effect for s in sigs]),
        index=[_row_id(s, multi) for s in sigs],
        columns=sym,
    )
    net = gene_sets.assign(source=gene_sets["source"].astype(str).map(bare))
    net, dropped = dedup_gene_sets(net, max_jaccard)
    ps = ulm(Y, net, min_size=min_size)
    p = ps.pvalue.to_numpy()
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    cols = curated_order(ps.score.columns)
    qd = pd.DataFrame(q, index=ps.score.index, columns=ps.score.columns)
    return ps.score[cols], qd[cols], dropped


def build_spec(
    score: pd.DataFrame,
    q: pd.DataFrame,
    names: Mapping[str, str],
    references: Sequence[str],
    focus: str | None,
    alpha: float,
    clip: float,
    ordering: str,
    source: str,
    n_dedup: int,
    max_jaccard: float,
    hidden: int,
    title: str,
    footer: Sequence[str],
) -> PlotSpec:
    rows = list(score.index)
    n_ref = sum(r in set(references) for r in rows)
    cols = list(score.columns)
    vals = score.to_numpy()
    mask = ~(q.to_numpy() < alpha)
    breaks = group_breaks(cols)
    xpos = np.arange(len(cols), dtype=float)
    for g in breaks:
        xpos[g + 1 :] += 0.6
    ypos = np.arange(len(rows), dtype=float) + (np.arange(len(rows)) >= n_ref) * (
        0.6 if n_ref else 0
    )
    labels = [names.get(r, r) for r in rows]
    if len(rows) > MAX_ROW_LABELS:
        labels = [lab if (i < n_ref or rows[i] == focus) else "" for i, lab in enumerate(labels)]
    p = Panel(
        key="fingerprint",
        x=Scale("categorical", categories=[pathway_label(c) for c in cols], tick_side="top",
                tick_rotation=50, domain=(-1.1, float(xpos[-1]) + 0.5)),
        y=Scale("categorical", categories=labels),
        colorbar=ColorbarSpec(label="pathway score (ULM t)", clip_note=f"clipped at ±{clip:g}",
                              compact=True),
        frame=False,
    )  # fmt: skip
    p.add(
        HeatmapMark(
            values=vals,
            mask=mask,
            cmap="diverging",
            clip=clip,
            row_gap_after=[n_ref - 1] if 0 < n_ref < len(rows) else [],
            col_gap_after=breaks,
        )
    )
    # group names under the matrix
    bottom = float(ypos[-1]) + 0.9
    starts = [0] + [b + 1 for b in breaks]
    ends = [b + 1 for b in breaks] + [len(cols)]
    for a, b in zip(starts, ends, strict=True):
        p.add(
            TextMark(
                x=float(xpos[a] + xpos[b - 1]) / 2,
                y=bottom,
                text=group_of(cols[a]),
                color="ink_soft",
                ha="center",
                va="top",
                size=6.5,
            )
        )
    if focus in rows:
        i = rows.index(focus)
        p.add(
            RectMark(
                x0=[-0.95],
                y0=[float(ypos[i]) - 0.4],
                width=[0.2],
                height=[0.8],
                color="focus",
                clip=False,
                categorical=True,
            )
        )
    # caption ----------------------------------------------------------------------------
    cap = Caption()
    lead: str = focus if focus is not None and focus in rows else rows[min(n_ref, len(rows) - 1)]
    others = [r for r in rows if r != lead]
    v = score.loc[lead].to_numpy(dtype=float)
    qv = q.loc[lead].to_numpy(dtype=float)
    nearest, r_best = "", np.nan
    if others:
        cors = {r: np.corrcoef(v, score.loc[r].to_numpy())[0, 1] for r in others}
        nearest = max(cors, key=lambda r: -np.inf if np.isnan(cors[r]) else cors[r])
        r_best = cors[nearest]
    sig_idx = [j for j in range(len(cols)) if qv[j] < alpha]
    if sig_idx:
        jb = max(sig_idx, key=lambda j: abs(v[j]))
        best = cols[jb]
        cap.add(
            TEMPLATES["fingerprint"],
            n=len(rows),
            m=len(cols),
            focus=names.get(lead, lead),
            nearest_row=f"{names.get(nearest, nearest)} (r = {r_best:.2f})" if nearest else "—",
            pathway=pathway_label(best),
            direction="up" if v[jb] > 0 else "down",
        )
    else:
        cap.add("{n} compounds × {m} pathways.", n=len(rows), m=len(cols))
        cap.weak(f"No pathway is significant for {names.get(lead, lead)} (q < {alpha:g}).")
    cap.add("Rows: {o}.", o=ordering)
    if hidden:
        cap.add("{h} pathway(s) with no significant cell are not shown.", h=hidden)
    if len(rows) > MAX_ROW_LABELS:
        cap.add("Only reference and focus rows are labelled ({n} rows).", n=len(rows))
    meth = Methods().use(
        "fingerprint",
        source=source,
        n_dedup=n_dedup,
        jac=max_jaccard,
        cells=int(np.isfinite(q.to_numpy()).sum()),
        alpha=alpha,
        clip=clip,
        ordering=ordering,
    )
    n_rows = len(rows)
    return PlotSpec(
        plot="fingerprint",
        panels=[p],
        title=title,
        caption=cap.render(),
        methods=meth.render(),
        footer=[
            f"Gene sets: {source}. Neutral = not significant (q ≥ {alpha:g}); "
            f"colour clipped at ±{clip:g}.",
            *footer,
        ],
        warnings=list(cap.notes),
        aspect=min(2.4, 0.3 + 0.022 * n_rows * (50 / max(len(cols), 10))),
        meta={
            "score": score,
            "q": q,
            "rows": rows,
            "columns": cols,
            "nearest": nearest,
            "clip": clip,
            "ordering": ordering,
        },
    ).validate()


def plot(
    data: SignatureSet | Sequence[Signature] | Mapping[str, Signature],
    gene_sets: pd.DataFrame | None = None,
    focus: str | None = None,
    references: Sequence[str] = (),
    dose: float | str | None = None,
    symbols: np.ndarray | None = None,
    alpha: float = 0.05,
    clip: float | None = None,
    order: str = "cluster",
    drop_empty: bool = True,
    max_jaccard: float = 0.5,
    gate: Collection[str] | None = None,
    source: str = "MSigDB Hallmark",
    title: str = "Fingerprint Strip",
    footer: Sequence[str] = (),
) -> PlotResult:
    """Fingerprint Strip for a set of signatures.

    ``data``: a :class:`SignatureSet` (top dose per compound unless ``dose`` is given) or
    signatures. ``gene_sets``: long ``source, target`` table; default MSigDB Hallmark from
    the cache. ``references``: compound ids drawn first, above a gap, in the given order.
    ``order``: ``"cluster"`` (average-linkage on score correlation, optimal leaf order) or
    ``"given"``. ``gate``: compound ids that passed retrievability; others are dropped.
    """
    from tincture.plots.compass import short_name

    if order not in ("cluster", "given"):
        raise ValueError("order must be 'cluster' or 'given'.")
    if gene_sets is None:
        from tincture.datasets.genesets import load_geneset

        gene_sets = load_geneset("hallmark").table
    if symbols is None and isinstance(data, SignatureSet):
        symbols = data.feature_symbols
    sigs = _rows(data, dose)
    dropped = [s.compound_id for s in sigs if gate is not None and s.compound_id not in gate
               and s.compound_id not in set(references)]  # fmt: skip
    sigs = [s for s in sigs if s.compound_id not in set(dropped)]
    if not sigs:
        raise TinctureDataError("All signatures failed the retrievability gate.")
    feats = sigs[0].features
    if any(not np.array_equal(s.features, feats) for s in sigs[1:]):
        sigs = [s.align(feats) for s in sigs]
    score, q, dedup = compute(sigs, gene_sets, symbols, max_jaccard)
    hidden = 0
    if drop_empty:
        keep = [c for c in score.columns if (q[c] < alpha).any()]
        hidden = len(score.columns) - len(keep)
        if not keep:
            raise TinctureDataError(
                f"No pathway is significant in any row (q < {alpha:g}).",
                fix="Pass drop_empty=False to show the empty matrix.",
            )
        score, q = score[keep], q[keep]
    multi = any("@" in r for r in score.index)
    names: dict[str, str] = {}
    for s in sigs:
        rid = _row_id(s, multi)
        nm = short_name(str(s.meta.get("compound_name") or s.compound_id))
        names[rid] = f"{nm} {fmt_conc(s.dose)}" if multi else nm
    refs = [r for r in references if r in score.index]
    queries = [r for r in score.index if r not in set(refs)]
    q_order = _order(score.loc[queries], order)
    ordering = (
        "references first in the given order, then "
        + ("average-linkage clustering of score correlations (optimal leaf order)"
           if order == "cluster" and len(queries) > 2 else "the given order")
    ) if refs else (
        "average-linkage clustering of score correlations (optimal leaf order)"
        if order == "cluster" and len(queries) > 2 else "the given order"
    )  # fmt: skip
    rows = refs + q_order
    score, q = score.loc[rows], q.loc[rows]
    c = clip if clip is not None else nice_clip(score.to_numpy()[q.to_numpy() < alpha])
    foc = focus
    if focus is not None and focus not in rows:
        cand = [r for r in rows if r.split("@")[0] == focus]
        if not cand:
            raise TinctureDataError(f"Focus {focus!r} is not among the rows.")
        foc = cand[-1]
    extra = [f"{len(dropped)} compound(s) failed the retrievability gate"] if dropped else []
    spec = build_spec(score, q, names, refs, foc, alpha, c, ordering, source, len(dedup),
                      max_jaccard, hidden, title, [*extra, *footer])  # fmt: skip
    if dropped:
        spec.caption += (
            f" {len(dropped)} compound(s) failed the retrievability gate and are not shown."
        )
    lead_sig = next(
        (s for s in sigs if _row_id(s, multi) == (foc or rows[len(refs) if refs else 0])),
        sigs[0],
    )
    p = data.pvalues.get((lead_sig.compound_id, lead_sig.dose)) if isinstance(
        data, SignatureSet) else None  # fmt: skip
    classic = build_classic("effect_ladder", sig=lead_sig, pvalues=p, symbols=symbols)
    if classic is not None:
        classic.meta["classic_of"] = "fingerprint"
    return PlotResult(spec=spec, classic=classic)
