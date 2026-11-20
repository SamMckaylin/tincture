"""CS2: therapeutic windows in PRISM (HANDOVER §6.2), dabrafenib in BRAF-driven melanoma.

Pre-specified (fixed before the first run):
* Drug: **dabrafenib** (BRAF V600E inhibitor). Sensitive lineage: **skin** (melanoma;
  PRISM ``primary_tissue``). Not every melanoma line carries BRAF V600E, so the check is a
  comparison of rates, not "all skin lines".
* On-target curve: each cell line's own viability curve. Toxicity proxy ("viability as the
  tox side"): one shared curve, the **median of all non-skin lines**, bootstrapped over lines.
* Window per line: EC90(line) < EC10(median non-skin), with the Window plot's fit guards.
* **Ground truth:** skin lines have a window more often than non-skin lines
  (one-sided Fisher exact p < 0.01).
* Reported, not scored: vemurafenib (same target) as a robustness check; a per-lineage strip
  ("Cell-line Atlas preview").

    python case_studies/cs2_prism_windows.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from _report import Report, md_table
from scipy.stats import fisher_exact

import tincture as tc
from tincture.io.prism import load_prism_secondary
from tincture.model.dose_response import DoseResponse
from tincture.plots.window import call_on, call_tox
from tincture.result import PlotResult
from tincture.spec import Panel, PlotSpec, PointMark, RefLineMark, Scale, TextMark

SENSITIVE = "skin"
FEATURED = {"ACH-000219": "A375 (melanoma, BRAF V600E)", "ACH-000681": "A549 (lung)"}
READOUT = "PRISM secondary viability"
SOURCE = ("PRISM Repurposing 20Q2 secondary (Corsello et al. 2020), CC BY 4.0",)


def windows(drug: str, n_boot: int = 100) -> tuple[pd.DataFrame, DoseResponse, object]:
    drs = load_prism_secondary(compounds=[drug])
    ctx, pts = drs.contexts, drs.points
    assert ctx is not None
    non = pts[~pts["context_id"].map(ctx["primary_tissue"]).eq(SENSITIVE)]
    tox = DoseResponse(
        drug,
        non["dose"].to_numpy(),
        non["response"].to_numpy(),
        replicate=non["context_id"].to_numpy(),
        assay=READOUT,
        context={
            "cell_line": f"median of {non['context_id'].nunique()} non-{SENSITIVE} lines",
            "readout": READOUT,
            "summary": "median",
        },
    )
    t = call_tox(tox, 0.1, n_boot=n_boot)
    cid = drs.keys()["compound_id"].iloc[0]
    rows = []
    for line in drs.keys()["context_id"]:
        dr = drs.get(cid, line)
        c = call_on(dr, 0.9, n_boot=0)  # point estimates for the screen-wide count
        exists = bool(c.usable and np.isfinite(c.ec) and c.ec < t.ec)
        rows.append(
            {
                "line": line,
                "ccle_name": ctx.loc[line, "ccle_name"],
                "lineage": ctx.loc[line, "primary_tissue"],
                "ec90": c.ec,
                "usable": c.usable,
                "why": "; ".join(c.problems),
                "window": exists,
                "fold": t.ec / c.ec if exists else np.nan,
            }
        )
    return pd.DataFrame(rows), tox, (drs, cid, t)


def atlas(df: pd.DataFrame, tox_ec: float, drug: str) -> PlotResult:
    """Per-lineage strip of EC90 (lines with a usable fit); window lines in terracotta. Each
    lineage is annotated with how many lines have a window and how many could not be judged
    (no usable fit), so omitted lines are never hidden."""
    order = df.groupby("lineage")["window"].mean().sort_values(ascending=False).index.tolist()[:14]
    weak = ~df["usable"] & df["why"].str.contains("maximal effect")
    counts = (
        df.assign(weak=weak, poor=~df["usable"] & ~weak)
        .groupby("lineage")
        .agg(n=("line", "size"), w=("window", "sum"), weak=("weak", "sum"), poor=("poor", "sum"))
    )
    d = df[df["usable"]].copy()
    d = d[d["lineage"].isin(order)]
    pos = {k: i for i, k in enumerate(order)}
    rng = np.random.default_rng(0)
    y = d["lineage"].map(pos).to_numpy() + rng.uniform(-0.25, 0.25, len(d))
    p = Panel(
        x=Scale("log", label=f"{drug} EC90 (M)", domain=(1e-9, 2e-3)),
        y=Scale(
            "linear",
            ticks=list(range(len(order))),
            tick_labels=order,
            domain=(len(order) - 0.4, -0.6),
        ),
    )
    w = d["window"].to_numpy()
    p.add(PointMark(x=d["ec90"][~w], y=y[~w], color="mid_grey", size=8, alpha=0.7))
    p.add(PointMark(x=d["ec90"][w], y=y[w], color="focus", size=10, categorical=True))
    p.add(
        RefLineMark(
            axis="x", value=tox_ec, color="reference", style="dashed", label="non-skin EC10"
        )
    )
    for k in order:
        c = counts.loc[k]
        p.add(
            TextMark(
                x=1.8e-3,
                y=pos[k],
                text=f"{int(c.w)}/{int(c.n)} window · {int(c.weak)} too little effect · "
                f"{int(c.poor)} poor fit",
                ha="right",
                va="center",
                color="ink_soft",
                size=6.5,
            )
        )
    spec = PlotSpec(
        plot="cs2_atlas",
        panels=[p],
        title=f"Cell-line Atlas preview: {drug} EC90 by lineage",
        caption=(
            "Each point is a cell line with a usable on-target fit; terracotta = window "
            "(EC90 below the toxicity EC10, dashed). Lines without a usable fit are not plotted "
            "but are counted on the right: 'too little effect' = maximal effect < 50% (a real "
            "negative), 'poor fit' = no convergence or R² < 0.8 (not judged). Lineages sorted "
            "by the share of lines with a window."
        ),
        footer=list(SOURCE),
        aspect=0.75,
    ).validate()
    return PlotResult(spec=spec)


def run() -> tuple[Path, bool]:
    rep = Report("cs2_prism_windows", "CS2: therapeutic windows in PRISM (dabrafenib)")
    rep.text(
        "Pre-specified in the script header: dabrafenib; sensitive lineage = skin; toxicity "
        "proxy = median of all non-skin lines; window = EC90(line) < EC10(proxy); ground "
        "truth = skin lines have a window more often (one-sided Fisher p < 0.01)."
    )
    df, tox, (drs, cid, t) = windows("dabrafenib")
    for line, label in FEATURED.items():
        on = drs.get(cid, line)
        on.context.update({"cell_line": label, "readout": READOUT})
        res = tc.window(
            on,
            tox,
            compound="Dabrafenib",
            n_boot=300,
            title=f"Window: dabrafenib in {label}",
            footer=SOURCE,
        )
        rep.text(f"## {label}")
        rep.figure(f"window_{line}", res)
    skin = df["lineage"] == SENSITIVE
    tab = pd.DataFrame(
        {
            "group": ["skin", "non-skin"],
            "lines": [int(skin.sum()), int((~skin).sum())],
            "usable on-target fit": [int(df[skin]["usable"].sum()), int(df[~skin]["usable"].sum())],
            "window": [int(df[skin]["window"].sum()), int(df[~skin]["window"].sum())],
        }
    )
    tab["share with window"] = (tab["window"] / tab["lines"]).map("{:.0%}".format)
    odds, p = fisher_exact(
        [
            [int(df[skin]["window"].sum()), int((~df[skin]["window"]).sum())],
            [int(df[~skin]["window"].sum()), int((~df[~skin]["window"]).sum())],
        ],
        alternative="greater",
    )
    by_lin = (
        df.groupby("lineage")
        .agg(lines=("line", "size"), window=("window", "sum"))
        .assign(share=lambda x: x["window"] / x["lines"])
        .sort_values("share", ascending=False)
        .head(8)
        .reset_index()
    )
    by_lin["share"] = by_lin["share"].map("{:.0%}".format)
    rep.text(
        "## Screen-wide",
        f"Toxicity proxy EC10 = {tc.report.fmt_conc(t.ec)} "
        f"(95% CI {tc.report.fmt_conc(t.ec_ci[0])}–{tc.report.fmt_conc(t.ec_ci[1])}).",
        md_table(tab),
        "Top lineages by share of lines with a window:",
        md_table(by_lin),
    )
    rep.figure("atlas_preview", atlas(df, t.ec, "dabrafenib"))
    # robustness (reported, not scored)
    df_v, _, (_, _, t_v) = windows("vemurafenib")
    sk = df_v["lineage"] == SENSITIVE
    rob = (
        f"Vemurafenib (same target, not scored): {df_v[sk]['window'].mean():.0%} of skin vs "
        f"{df_v[~sk]['window'].mean():.0%} of non-skin lines have a window "
        f"(proxy EC10 {tc.report.fmt_conc(t_v.ec)})."
    )
    rep.text("## Robustness", rob, "## Ground-truth check")
    ok = rep.check(
        "Sensitive lineage shows a window that other lines lack",
        bool(odds > 1 and p < 0.01),
        f"{int(df[skin]['window'].sum())}/{int(skin.sum())} skin lines vs "
        f"{int(df[~skin]['window'].sum())}/{int((~skin).sum())} non-skin lines have a window "
        f"(one-sided Fisher exact, odds ratio {odds:.1f}, p = {p:.2g}; criterion p < 0.01).",
    )
    rep.text("", f"**Overall: {'PASS' if ok else 'FAIL'}**")
    return rep.write(), ok


if __name__ == "__main__":
    path, ok = run()
    print(path, "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
