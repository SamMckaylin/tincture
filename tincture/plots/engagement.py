"""Engagement Ladder: does potency hold up from protein to cell to phenotype, and where does
it break? (HANDOVER §5.10)

Rungs, top to bottom: biochemical (Kd / enzymatic IC50) → cellular target engagement →
pathway (optional) → phenotype (viability IC50) → exposure (free Cmax, optional). x is a
shared log-concentration axis; one polyline per compound, interval whiskers per rung, small
multiples per target.

Statistics: the potency shift between consecutive rungs (Δlog10). The expected ATP shift
for ATP-competitive kinase inhibitors (Cheng–Prusoff, ``log10(1 + [ATP]/Km,ATP)`` for
cellular ATP 1–5 mM) is drawn as a grey wedge, only when the kinase's ATP Km is supplied.

Automatic flags (each tested):
* **phenotype > 10× more potent than target engagement** → "phenotype likely driven by
  something other than {target}" (terracotta);
* **biochemical → cellular shift beyond the ATP wedge** → "permeability, efflux or binding
  issue" (purple); not evaluated without an ATP Km, and the caption says so;
* **exposure below the phenotype rung** → "unlikely to reach active concentration in vivo".

Guards (tested): rungs measured in different cell systems are joined by **dashed**
segments and each point carries its cell system; missing rungs are skipped, never
interpolated; censored values are hollow with ">" and only used where the bound is
conclusive. Classic fallback: grouped bars of potency by assay.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from itertools import pairwise

import numpy as np
import pandas as pd

from tincture.errors import TinctureDataError
from tincture.model.potency import PotencyTable
from tincture.report.captions import Caption, fmt_conc
from tincture.report.classic import build_classic, register_classic
from tincture.report.methods import Methods, register_block
from tincture.result import PlotResult
from tincture.spec import (
    IntervalMark,
    LineMark,
    NoticeMark,
    Panel,
    PlotSpec,
    PointMark,
    PolygonMark,
    RectMark,
    Scale,
    TextMark,
)

RUNGS = ("biochemical", "cellular", "pathway", "phenotype", "exposure")
RUNG_LABELS = {
    "biochemical": "Biochemical",
    "cellular": "Cellular engagement",
    "pathway": "Pathway",
    "phenotype": "Phenotype",
    "exposure": "Exposure (free Cmax)",
}
OFFTARGET_FOLD = 10.0
CELLULAR_ATP = (1e-3, 5e-3)
LINE_COLOURS = ["reference", "novel", "good", "mid_grey"]

register_block(
    "engagement",
    "Each rung's potency is plotted on a shared log-concentration axis ({rungs}); whiskers "
    "are the intervals supplied with each value. The biochemical→cellular shift is "
    "log10(cellular IC50 / biochemical potency). {atp} A phenotype is flagged as likely "
    "off-target when it is more than {fold:.0f}-fold more potent than cellular engagement of "
    "the named target (biochemical potency when no cellular rung exists). Rungs measured in "
    "different cell systems are joined by dashed segments; missing rungs are skipped.",
)


@dataclass
class Rung:
    rung: str
    px: float
    px_lo: float
    px_hi: float
    censored: str
    cell_system: str
    label: str

    @property
    def molar(self) -> float:
        return float(10.0 ** (-self.px))


@dataclass
class LadderResult:
    compound: str
    target: str
    rungs: list[Rung]
    shift_bc: float | None  # log10 shift biochemical → cellular (positive = weaker in cells)
    shift_bc_bound: str = ""  # ">" when the cellular value is censored
    atp_shift: tuple[float, float] | None = None
    flags: dict[str, str] = field(default_factory=dict)
    not_evaluated: dict[str, str] = field(default_factory=dict)

    def get(self, rung: str) -> Rung | None:
        return next((r for r in self.rungs if r.rung == rung), None)


def atp_shift(km_atp: float, atp: tuple[float, float] = CELLULAR_ATP) -> tuple[float, float]:
    """Expected Cheng–Prusoff shift (log10) for ATP-competitive inhibitors."""
    if km_atp <= 0:
        raise TinctureDataError("ATP Km must be positive (molar).")
    return (float(np.log10(1 + atp[0] / km_atp)), float(np.log10(1 + atp[1] / km_atp)))


def compute(
    d: pd.DataFrame,
    compound: str,
    target: str,
    km_atp: float | None = None,
    atp: tuple[float, float] = CELLULAR_ATP,
) -> LadderResult:
    rows = d[(d["compound_id"] == compound) & (d["target"] == target)]
    if rows.empty:
        raise TinctureDataError(f"No potencies for {compound} on {target}.")
    dup = rows["rung"].value_counts()
    if (dup > 1).any():
        raise TinctureDataError(
            f"{compound}/{target}: more than one value for rung(s) {list(dup[dup > 1].index)}.",
            fix="Keep one value per rung (e.g. the primary assay).",
        )
    rungs = []
    for name in RUNGS:
        r = rows[rows["rung"] == name]
        if r.empty:
            continue  # skipped, never interpolated
        x = r.iloc[0]
        px = float(x["px"])
        rungs.append(
            Rung(
                rung=name,
                px=px,
                px_lo=float(x["px_lo"]) if "px_lo" in x and pd.notna(x["px_lo"]) else px,
                px_hi=float(x["px_hi"]) if "px_hi" in x and pd.notna(x["px_hi"]) else px,
                censored=str(x.get("censored", "") or ""),
                cell_system=str(x.get("cell_system", "") or ""),
                label=str(x.get("label", "") or ""),
            )
        )
    res = LadderResult(compound, target, rungs, None)
    b, c, ph, ex = (res.get(k) for k in ("biochemical", "cellular", "phenotype", "exposure"))
    if b is not None and c is not None and not b.censored:
        res.shift_bc = b.px - c.px
        res.shift_bc_bound = {">": ">", "<": "<"}.get(c.censored, "")
    if km_atp is not None:
        res.atp_shift = atp_shift(km_atp, atp)
    # flag 1: phenotype more potent than target engagement
    eng = c or b
    if ph is not None and eng is not None:
        # a censored value is used only when its bound makes the comparison conclusive:
        # phenotype '<' (even more potent) and engagement '>' (even weaker) both strengthen it
        conclusive = ph.censored in ("", "<") and eng.censored in ("", ">")
        gap = ph.px - eng.px
        if gap > np.log10(OFFTARGET_FOLD) and conclusive:
            fold = 10**gap
            more = "≥ " if (ph.censored or eng.censored) else ""
            res.flags["off_target"] = (
                f"phenotype {more}{fold:.0f}× more potent than {target} engagement: likely "
                f"driven by something other than {target}"
            )
        elif gap > np.log10(OFFTARGET_FOLD):  # bound points the wrong way: cannot conclude
            res.not_evaluated["off_target"] = "a censored potency makes it inconclusive"
    # flag 2: biochemical → cellular shift beyond the ATP wedge
    if res.shift_bc is not None:
        if res.atp_shift is None:
            res.not_evaluated["shift"] = f"no ATP Km supplied for {target}"
        elif res.shift_bc > res.atp_shift[1]:
            res.flags["shift"] = (
                f"{res.shift_bc_bound}{res.shift_bc:.1f} log shift into cells, beyond the "
                f"{res.atp_shift[1]:.1f} expected from ATP: permeability, efflux or binding issue"
            )
    # flag 3: exposure below phenotype
    if ex is not None and ph is not None and ph.censored in ("", ">") and ex.px > ph.px:
        res.flags["exposure"] = (
            f"free Cmax {fmt_conc(ex.molar)} is below the phenotype potency "
            f"{fmt_conc(ph.molar)}: unlikely to reach active concentration in vivo"
        )
    return res


def _panel(results: list[LadderResult], key: str, lo: float, hi: float) -> Panel:
    present = [r for r in RUNGS if any(res.get(r) for res in results)]
    ypos = {r: i for i, r in enumerate(present)}
    p = Panel(
        key=key,
        title=results[0].target,
        x=Scale("log", label="concentration (M)", domain=(lo, hi)),
        y=Scale(
            "linear",
            ticks=[float(ypos[r]) for r in present],
            tick_labels=[RUNG_LABELS[r] for r in present],
            domain=(len(present) - 0.4, -0.7),
        ),
    )
    for k, res in enumerate(results):
        col = LINE_COLOURS[k % len(LINE_COLOURS)]
        rs = res.rungs
        b, c = res.get("biochemical"), res.get("cellular")
        if res.atp_shift is not None and b is not None and c is not None and not b.censored:
            x0 = b.molar
            yb, yc = ypos["biochemical"], ypos["cellular"]
            p.add(
                PolygonMark(
                    x=[x0, x0 * 10 ** res.atp_shift[0], x0 * 10 ** res.atp_shift[1]],
                    y=[yb, yc, yc],
                    color="neutral",
                    alpha=0.55,
                    z=0.5,
                )
            )
            p.add(
                TextMark(
                    x=x0 * 10 ** res.atp_shift[1],
                    y=yc + 0.32,
                    text="expected from ATP",
                    ha="center",
                    va="top",
                    color="ink_soft",
                    size=6.0,
                )
            )
        for a, bb in pairwise(rs):
            same = a.cell_system and a.cell_system == bb.cell_system
            p.add(
                LineMark(
                    x=[a.molar, bb.molar],
                    y=[ypos[a.rung], ypos[bb.rung]],
                    color=col,
                    width=1.6,
                    style="solid" if same else "dashed",
                    categorical=True,
                    z=2,
                )
            )
        for r in rs:
            y = ypos[r.rung]
            if r.px_lo != r.px_hi:
                p.add(
                    IntervalMark(
                        pos=[y],
                        lo=[10 ** (-r.px_hi)],
                        hi=[10 ** (-r.px_lo)],
                        orientation="h",
                        color=col,
                        width=1.2,
                    )
                )
            p.add(
                PointMark(
                    x=[r.molar],
                    y=[y],
                    color=col,
                    size=30,
                    hollow=bool(r.censored),
                    categorical=True,
                    z=3,
                )
            )
            tag = f"{r.censored + ' ' if r.censored else ''}{fmt_conc(r.molar)}"
            sys_ = f" · {r.cell_system}" if r.cell_system else ""
            p.add(
                TextMark(
                    x=r.molar * 1.8,
                    y=y - 0.2,
                    text=f"{tag}{sys_}",
                    ha="left",
                    va="bottom",
                    color="ink_soft",
                    size=6.0,
                )
            )
        if len(results) > 1:
            last = rs[-1]
            p.add(
                TextMark(
                    x=last.molar,
                    y=ypos[last.rung] + 0.3,
                    text=res.compound,
                    color=col,
                    ha="left",
                    va="top",
                    categorical=True,
                )
            )
    # flags (the point of the plot)
    lines = []
    for res in results:
        who = f"{res.compound}: " if len(results) > 1 else ""
        for kind, colour in (("off_target", "focus"), ("shift", "novel"), ("exposure", "warn")):
            if kind in res.flags:
                lines.append((f"Flag: {who}{res.flags[kind]}", colour))
    for j, (text, colour) in enumerate(lines):
        p.add(
            TextMark(
                x=0.01,
                y=0.02 + 0.075 * (len(lines) - 1 - j),
                coords="axes",
                text=text,
                color=colour,
                ha="left",
                va="bottom",
                size=6.5,
                categorical=True,
            )
        )
    if not lines:
        p.add(NoticeMark(text="no flags", loc="lower right", level="info"))
    return p


def build_spec(
    results: list[LadderResult], title: str, footer: Sequence[str], atp: tuple[float, float]
) -> PlotSpec:
    by_target: dict[str, list[LadderResult]] = {}
    for r in results:
        by_target.setdefault(r.target, []).append(r)
    allx = [rg.molar for r in results for rg in r.rungs]
    allx += [10 ** (-rg.px_hi) for r in results for rg in r.rungs]
    allx += [10 ** (-rg.px_lo) for r in results for rg in r.rungs]
    for r in results:
        b = r.get("biochemical")
        if r.atp_shift is not None and b is not None:
            allx.append(b.molar * 10 ** r.atp_shift[1])
    lo, hi = min(allx) / 6, max(allx) * 400  # room for value labels on the right
    panels = []
    ncols = 1 if len(by_target) <= 4 else 2  # a ladder needs a wide concentration axis
    for i, rs in enumerate(by_target.values()):
        pnl = _panel(rs, f"t{i}", lo, hi)
        pnl.row, pnl.col = divmod(i, ncols)
        panels.append(pnl)
    nrows = int(np.ceil(len(panels) / ncols))
    cap = Caption()
    for r in results:
        c = r.get("cellular")
        ph = r.get("phenotype")
        if r.shift_bc is not None:
            exp = (
                f"{r.atp_shift[0]:.1f}–{r.atp_shift[1]:.1f} expected from ATP"
                if r.atp_shift
                else "ATP expectation not computed: no ATP Km"
            )
        eng = c or r.get("biochemical")
        if ph is not None and eng is not None and ph.censored in ("", "<") and not eng.censored:
            ratio = 10 ** (ph.px - eng.px)
            rel = (
                f"{ratio:.0f}× more potent than"
                if ratio >= 1.5
                else f"{1 / ratio:.0f}× weaker than"
                if ratio <= 1 / 1.5
                else "in line with"
            )
            rel += " cellular" if c is not None else " biochemical"
        else:
            rel = "not comparable to"
        if r.shift_bc is not None:
            cap.add(
                "{c} on {t}: {b}{s:.1f} log shift biochemical→cellular ({e}); phenotype is "
                "{rel} target engagement.",
                c=r.compound,
                t=r.target,
                b=r.shift_bc_bound,
                s=r.shift_bc,
                e=exp,
                rel=rel,
            )
        else:
            cap.add(
                "{c} on {t}: phenotype is {rel} target engagement.",
                c=r.compound,
                t=r.target,
                rel=rel,
            )
        for v in r.flags.values():
            cap.add("Flag: {v}.", v=v)
        for k, why in r.not_evaluated.items():
            cap.weak(
                f"{r.compound}/{r.target}: {'ATP-shift' if k == 'shift' else 'off-target'} "
                f"flag not evaluated ({why})."
            )
    if any(
        rg.cell_system and nx.cell_system != rg.cell_system
        for r in results
        for rg, nx in pairwise(r.rungs)
    ):
        cap.add("Dashed segments join rungs measured in different cell systems.")
    rungs_used = [r for r in RUNGS if any(x.get(r) for x in results)]
    meth = Methods().use(
        "engagement",
        rungs=", ".join(RUNG_LABELS[r].lower() for r in rungs_used),
        atp=(
            "Expected ATP shift: Cheng–Prusoff log10(1 + [ATP]/Km,ATP) with cellular ATP "
            f"{fmt_conc(atp[0])}–{fmt_conc(atp[1])} (grey wedge)."
            if any(r.atp_shift for r in results)
            else "No ATP Km was supplied, so no expected ATP shift is drawn and the "
            "permeability/efflux flag is not evaluated."
        ),
        fold=OFFTARGET_FOLD,
    )
    return PlotSpec(
        plot="engagement_ladder",
        panels=panels,
        nrows=nrows,
        ncols=ncols,
        title=title,
        caption=cap.render(),
        methods=meth.render(),
        footer=[
            "Hollow = censored (> weaker / < more potent than the value shown); dashed = "
            "different cell systems",
            *footer,
        ],
        warnings=list(cap.notes),
        aspect=(0.24 + 0.3 * nrows) if ncols == 1 else 0.36 * nrows,
        meta={"results": results},
    ).validate()


@register_classic("engagement_ladder", "grouped bar chart of potencies")
def build_classic_spec(results: list[LadderResult] | None = None, **_: object) -> PlotSpec:
    if not results:
        raise TinctureDataError("Classic needs ladder results.")
    labels, vals, cols = [], [], []
    for k, r in enumerate(results):
        for rg in r.rungs:
            labels.append(f"{r.compound} {r.target} · {RUNG_LABELS[rg.rung]}")
            vals.append(rg.px)
            cols.append(LINE_COLOURS[k % len(LINE_COLOURS)])
    x = np.arange(len(vals))
    p = Panel(
        x=Scale("linear", ticks=[float(v) for v in x], tick_labels=labels, tick_rotation=60),
        y=Scale(label="pIC50 / pKd"),
    )
    p.add(
        RectMark(
            x0=x - 0.4,
            y0=np.zeros(len(x)),
            width=np.full(len(x), 0.8),
            height=np.array(vals),
            color=cols,
        )
    )
    return PlotSpec(
        plot="engagement_classic",
        panels=[p],
        title="Potency by assay",
        caption="Potency by assay (bar chart).",
        aspect=0.6,
    ).validate()


def plot(
    table: PotencyTable,
    pairs: Sequence[tuple[str, str]],
    km_atp: Mapping[str, float] | None = None,
    atp: tuple[float, float] = CELLULAR_ATP,
    title: str = "Engagement Ladder",
    footer: Sequence[str] = (),
) -> PlotResult:
    """Engagement Ladder for (compound, target) ``pairs``.

    ``table`` needs a ``rung`` column with values from :data:`RUNGS`; optional columns
    ``px_lo``/``px_hi`` (interval), ``cell_system`` and ``label``. ``km_atp``: ATP Km
    (molar) per target, for the Cheng–Prusoff wedge; never pass a guess.
    """
    d = table.data
    if "rung" not in d.columns:
        raise TinctureDataError("The PotencyTable needs a `rung` column.", fix=f"Use {RUNGS}.")
    bad = sorted(set(d["rung"]) - set(RUNGS))
    if bad:
        raise TinctureDataError(f"Unknown rungs {bad}.", fix=f"Use {list(RUNGS)}.")
    km = dict(km_atp or {})
    results = [compute(d, c, t, km.get(t), atp) for c, t in pairs]
    spec = build_spec(results, title, footer, atp)
    return PlotResult(spec=spec, classic=build_classic("engagement_ladder", results=results))
