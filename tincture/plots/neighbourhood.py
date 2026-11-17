"""Neighbourhood Card: which annotated compounds does this one resemble? (HANDOVER §5.9)

Query structure (RDKit) next to its top-k neighbours in an annotated reference library:
similarity bars with 95% intervals, coloured by mechanism of action, against the library's
own null (similarity of random library pairs). Statistics in
:mod:`tincture.stats.neighbours`. A query that is itself in the library is compared with
the library *without* it.

Guard (tested): if no neighbour passes q < 0.05, the card says **"no confident
neighbour"** in place of a ranked list.
Classic fallback: ranked similarity bars without null, intervals or q values.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from tincture.errors import TinctureDataError
from tincture.model.reference import ReferenceLibrary
from tincture.model.signature import Signature
from tincture.report.captions import Caption, fmt_p
from tincture.report.classic import build_classic, register_classic
from tincture.report.methods import Methods, register_block
from tincture.result import PlotResult
from tincture.spec import (
    ImageMark,
    IntervalMark,
    NoticeMark,
    Panel,
    PlotSpec,
    RectMark,
    RefLineMark,
    Scale,
    TextMark,
)
from tincture.stats.neighbours import NeighbourResult, neighbours

MOA_COLOURS = ["reference", "novel", "good", "focus_light", "reference_light", "good_light"]

register_block(
    "neighbourhood",
    "Neighbours are ranked by cosine similarity between the query signature and each "
    "reference in {library} ({n_refs} references{loo}). The null distribution is the "
    "similarity of {n_null:,} distinct reference pairs in the library ({null_desc}); "
    "p = share of null "
    "pairs at least as similar (+1 correction), BH-adjusted over all references (q). {ci}",
)


def structure_image(smiles: str | None, size: int = 360) -> np.ndarray | None:
    """RGBA image of the molecule, or None if no SMILES / RDKit unavailable / unparsable."""
    if not smiles or smiles in ("nan", "None", "-666"):
        return None
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
    except ImportError:
        return None
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    img = Draw.MolToImage(mol, size=(size, size))
    return np.asarray(img.convert("RGBA"))


def _ann(lib: ReferenceLibrary, ref: str, key: str, default: str = "") -> str:
    v = lib.annotations.get(ref, {}).get(key, default)
    return default if v in ("", "nan", "None") else v


def build_spec(
    res: NeighbourResult,
    lib: ReferenceLibrary,
    query_name: str,
    smiles: str | None,
    k: int,
    alpha: float,
    loo: bool,
    ci_note: str,
    title: str,
    footer: Sequence[str],
    null_kind: str = "all",
) -> PlotSpec:
    top = list(range(min(k, len(res.names))))
    conf = [i for i in top if res.q[i] < alpha]
    # structure panel
    left = Panel(
        key="structure",
        x=Scale("linear", domain=(0, 1), visible=False),
        y=Scale("linear", domain=(0, 1), visible=False),
        frame=False,
        title=query_name,
    )
    img = structure_image(smiles)
    if img is not None:
        left.add(ImageMark(image=img, extent=(0.0, 1.0, 0.0, 1.0)))
    else:
        left.add(
            TextMark(x=0.5, y=0.5, text="structure not available", ha="center", color="ink_soft")
        )
    # neighbour panel
    cap = Caption()
    hi = float(max(0.2, np.nanmax(res.similarity[top]) if top else 0.2))
    if res.ci is not None and top:
        hi = max(hi, float(np.nanmax(res.ci[top, 1])))
    right = Panel(
        key="neighbours",
        x=Scale("linear", label="cosine similarity to the query", domain=(0, hi * 1.25)),
        y=Scale("linear", visible=False, domain=(k - 0.4, -1.2)),
        frame=True,
    )
    p95, p99 = (float(np.quantile(res.null, v)) for v in (0.95, 0.99))
    right.add(RefLineMark(axis="x", value=p95, color="ink_soft", style="dotted", label="null 95%"))
    right.add(RefLineMark(axis="x", value=p99, color="ink_soft", style="dashed", label="99%"))
    if not conf:
        best = 0
        right.add(NoticeMark(text="no confident neighbour", loc="center", level="warning"))
        cap.weak(
            f"No neighbour passes q < {alpha:g}: the closest, "
            f"{_ann(lib, res.names[best], 'name', res.names[best])} "
            f"({_ann(lib, res.names[best], 'moa', 'no MoA')}), has cosine "
            f"{res.similarity[best]:.2f} (q {fmt_p(res.q[best])}), within what random "
            "library pairs reach."
        )
    else:
        moas = list(dict.fromkeys(_ann(lib, res.names[i], "moa", "no MoA") for i in conf))
        col_of = {m: MOA_COLOURS[j] for j, m in enumerate(moas[: len(MOA_COLOURS)])}
        for row, i in enumerate(top):
            ref = res.names[i]
            moa = _ann(lib, ref, "moa", "no MoA")
            ok = res.q[i] < alpha
            col = col_of.get(moa, "mid_grey") if ok else "neutral"
            right.add(
                RectMark(
                    x0=[0],
                    y0=[row - 0.12],
                    width=[res.similarity[i]],
                    height=[0.5],
                    color=col,
                    categorical=ok,
                )
            )
            if res.ci is not None:
                right.add(
                    IntervalMark(
                        pos=[row + 0.13],
                        lo=[res.ci[i, 0]],
                        hi=[res.ci[i, 1]],
                        orientation="h",
                        color="ink_soft",
                        width=1.0,
                    )
                )
            right.add(
                TextMark(
                    x=hi * 1.24,
                    y=row + 0.13,
                    text=f"q {fmt_p(res.q[i])}",
                    ha="right",
                    color="ink" if ok else "ink_soft",
                    size=6.5,
                )
            )
            right.add(
                TextMark(
                    x=0.005,
                    y=row - 0.16,
                    text=f"{_ann(lib, ref, 'name', ref)} · {moa}",
                    ha="left",
                    va="bottom",
                    size=6.5,
                    color="ink" if ok else "ink_soft",
                )
            )
        b = conf[0]
        same = sum(_ann(lib, res.names[i], "moa") == _ann(lib, res.names[b], "moa") for i in conf)
        cap.add(
            "{q} most resembles {n} ({m}; cosine {s:.2f}, q {qq}). {c} of the top {k} "
            "neighbours pass q < {a:g}; {same} of them share this MoA.",
            q=query_name,
            n=_ann(lib, res.names[b], "name", res.names[b]),
            m=_ann(lib, res.names[b], "moa", "no MoA"),
            s=res.similarity[b],
            qq=fmt_p(res.q[b]),
            c=len(conf),
            k=len(top),
            a=alpha,
            same=same,
        )
        if len(conf) < len(top):
            cap.add("Grey bars do not pass q < {a:g}.", a=alpha)
    if res.ci is None:
        cap.weak("no replicate signatures or SEs for the query, so similarities have no interval.")
    meth = Methods().use(
        "neighbourhood",
        library=lib.label(),
        n_refs=len(res.names),
        loo="; the query itself was removed" if loo else "",
        n_null=res.n_null,
        null_desc="all pairs" if null_kind == "all" else "pairs with different MoA annotations",
        ci=ci_note,
    )
    return PlotSpec(
        plot="neighbourhood",
        panels=[left, right],
        ncols=2,
        width_ratios=[1.0, 2.2],
        title=title,
        caption=cap.render(),
        methods=meth.render(),
        footer=[
            f"Library: {lib.label()} ({len(res.names)} references); null = "
            + ("all library pairs" if null_kind == "all" else "library pairs with different MoA"),
            *footer,
        ],
        warnings=list(cap.notes),
        aspect=0.5,
        meta={"result": res, "confident": [res.names[i] for i in conf]},
    ).validate()


@register_classic("neighbourhood", "ranked similarity bars")
def build_classic_spec(
    res: NeighbourResult | None = None,
    lib: ReferenceLibrary | None = None,
    k: int = 10,
    **_: object,
) -> PlotSpec:
    if res is None or lib is None:
        raise TinctureDataError("Classic needs the neighbour result.")
    top = list(range(min(k, len(res.names))))
    p = Panel(
        x=Scale(label="similarity"),
        y=Scale(
            "linear",
            ticks=[float(i) for i in top],
            tick_labels=[_ann(lib, res.names[i], "name", res.names[i]) for i in top],
            domain=(len(top) - 0.4, -0.6),
        ),
    )
    for row, i in enumerate(top):
        p.add(
            RectMark(
                x0=[0], y0=[row - 0.35], width=[res.similarity[i]], height=[0.7], color="mid_grey"
            )
        )
    return PlotSpec(
        plot="neighbourhood_classic",
        panels=[p],
        title="Top similar compounds",
        caption="Top compounds by similarity.",
        aspect=0.6,
    ).validate()


def plot(
    query: Signature,
    library: ReferenceLibrary,
    smiles: str | None = None,
    k: int = 10,
    alpha: float = 0.05,
    query_replicates: np.ndarray | None = None,
    name: str | None = None,
    null: np.ndarray | None = None,
    null_kind: str = "all",
    title: str = "Neighbourhood Card",
    footer: Sequence[str] = (),
) -> PlotResult:
    """Neighbourhood Card for ``query`` against an annotated ``library``.

    ``query_replicates``: (n, features) signatures of the same compound (e.g. other cell
    lines) for bootstrap intervals; otherwise the query's SEs are used if present.
    ``null``: precomputed library-pair similarities (reuse across many queries).
    ``null_kind``: ``"all"`` library pairs (HANDOVER default) or ``"unrelated"`` pairs (pairs
    sharing an MoA annotation removed; see :mod:`tincture.stats.neighbours`).
    """
    if k < 1:
        raise ValueError("k must be ≥ 1.")
    loo = bool(library.reference_of(query.compound_id))
    lib = library.without(query.compound_id) if loo else library
    if null is not None and loo:
        null = None  # the null must come from the library actually searched
    q = query.align(lib.features)
    reps = None
    if query_replicates is not None:
        reps = np.asarray(query_replicates, float)
    res = neighbours(
        q.effect, lib, query_replicates=reps, query_se=q.se, null=null, null_kind=null_kind
    )
    ci_note = (
        "Intervals: 95% bootstrap over the query's replicate signatures."
        if reps is not None and reps.shape[0] >= 2
        else "Intervals: 95% from posterior draws of the query (per-feature SE)."
        if q.se is not None
        else "No intervals (no replicates or SEs for the query)."
    )
    nm = name or str(query.meta.get("compound_name") or query.compound_id)
    smi = smiles if smiles is not None else query.meta.get("smiles")
    spec = build_spec(res, lib, nm, smi, k, alpha, loo, ci_note, title, footer, null_kind)
    return PlotResult(spec=spec, classic=build_classic("neighbourhood", res=res, lib=lib, k=k))
