"""CS8: Onset Timeline on the MIX-seq trametinib time course (HANDOVER §6.2; D10).

Data: McFarland, Paolella et al. 2020, experiment 5: 24 pooled cancer cell lines, trametinib
0.1 µM or DMSO, hash-tagged at 3, 6, 12, 24 and 48 h (figshare 10298696).

Pre-specified (fixed before the first run):
* Pseudobulk UMI counts per cell line × condition × time (QC-passing singlets); log1p CPM;
  genes with mean CPM ≥ 10. Signature per cell line and time = trametinib − time-matched
  DMSO. Pathway scores per signature with ULM (PROGENy top-500, all 14; MSigDB Hallmark:
  apoptosis, p53 pathway, unfolded protein response, E2F targets, G2M checkpoint, MYC
  targets v1).
* Per pathway and time: mean score across cell lines and a one-sample t-test; BH across
  all pathway × time tests. Active: q < 0.05 and |mean score| ≥ 2.
* **Check:** the direct pathway output, PROGENy MAPK, has an onset that is earlier than that
  of every stress/apoptosis programme (Hallmark apoptosis, Hallmark p53 pathway, PROGENy p53,
  Hallmark unfolded protein response): for each, either it never becomes active or its
  onset interval starts at or after the end of MAPK's onset interval.

    python case_studies/cs8_onset_mixseq.py
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
from scipy import stats as sps

import tincture as tc
from tincture.datasets.genesets import load_geneset
from tincture.io.mixseq import pseudobulk_trametinib_tc
from tincture.stats.pathways import ulm

HALLMARKS = {
    "APOPTOSIS": "apoptosis (Hallmark)",
    "P53_PATHWAY": "p53 pathway (Hallmark)",
    "UNFOLDED_PROTEIN_RESPONSE": "unfolded protein response (Hallmark)",
    "E2F_TARGETS": "E2F targets (Hallmark)",
    "G2M_CHECKPOINT": "G2/M checkpoint (Hallmark)",
    "MYC_TARGETS_V1": "MYC targets (Hallmark)",
}
DIRECT = "MAPK (PROGENy)"
STRESS = [
    "apoptosis (Hallmark)",
    "p53 pathway (Hallmark)",
    "p53 (PROGENy)",
    "unfolded protein response (Hallmark)",
]
THRESHOLD = 2.0
FOOT = ("MIX-seq, McFarland, Paolella et al., Nat Commun 2020 (CC BY 4.0)",)


def signatures() -> tuple[pd.DataFrame, pd.DataFrame]:
    counts, groups = pseudobulk_trametinib_tc()
    cpm = counts / counts.sum(axis=0) * 1e6
    keep = cpm.mean(axis=1) >= 10
    lc = np.log1p(cpm[keep])
    g = groups.set_index("group")
    rows, meta = [], []
    for (dep, t), gg in g[g["condition"] != "untreated"].groupby(["depmap_id", "time_h"]):
        cond = dict(zip(gg["condition"], gg.index, strict=True))
        if {"DMSO", "trametinib"} <= set(cond):
            rows.append((lc[cond["trametinib"]] - lc[cond["DMSO"]]).rename(f"{dep}|{t:g}"))
            meta.append({"sig": f"{dep}|{t:g}", "depmap_id": dep, "time": float(t)})
    return pd.DataFrame(rows), pd.DataFrame(meta).set_index("sig")


def run() -> tuple[Path, bool]:
    rep = Report("cs8_onset_mixseq", "CS8: Onset Timeline, MIX-seq trametinib time course")
    rep.text("Pre-specified in the script header.")
    Y, meta = signatures()
    prog = load_geneset("progeny")
    hall = load_geneset("hallmark")
    sp = ulm(Y, prog.table, source=prog.version).score
    sp.columns = [f"{c} (PROGENy)" for c in sp.columns]
    hset = hall.table[hall.table["source"].isin(HALLMARKS)]
    sh = ulm(Y, hset, source=hall.version).score.rename(columns=HALLMARKS)
    S = pd.concat([sp, sh], axis=1).join(meta)
    long = []
    for (t,), g in S.groupby(["time"]):
        for f in [c for c in S.columns if c not in meta.columns]:
            v = g[f].to_numpy(float)
            v = v[np.isfinite(v)]
            res = sps.ttest_1samp(v, 0.0)
            long.append(
                {"feature": f, "time": t, "score": float(v.mean()), "p": float(res.pvalue),
                 "n_lines": len(v)}
            )  # fmt: skip
    L = pd.DataFrame(long)
    n_lines = int(L["n_lines"].min())
    res = tc.onset_timeline(
        L,
        threshold=THRESHOLD,
        score_label="mean ULM score across cell lines, trametinib − DMSO",
        title="Onset Timeline: trametinib 0.1 µM in 24 pooled cell lines (MIX-seq)",
        footer=(
            *FOOT,
            "Gene sets: PROGENy top-500, MSigDB Hallmark (OmniPath, fetched 2026-11-17)",
        ),
    )
    rep.text(f"Signatures: {len(Y)} (cell line × time), {Y.shape[1]} genes, ≥ {n_lines} lines "
             "per time point.")  # fmt: skip
    rep.figure("onset_trametinib", res)
    ons = {o.feature: o for o in res.meta["onsets"]}
    tab = pd.DataFrame(
        [
            {
                "pathway": o.feature,
                "direction": {1: "up", -1: "down", 0: "—"}[o.direction],
                "onset": f"({o.onset[0]:g}, {o.onset[1]:g}] h" if o.onset else "never",
                "peak": f"{o.peak_time:g} h ({o.peak_score:+.1f})" if o.onset else "—",
                "resolution": (
                    "sustained" if o.sustained
                    else f"({o.resolution[0]:g}, {o.resolution[1]:g}] h" if o.resolution
                    else "—"
                ),
            }
            for o in (ons[f] for f in res.meta["order"])
        ]
    )  # fmt: skip
    rep.text(md_table(tab))
    rep.text("## Ground-truth check")
    d = ons[DIRECT]
    if d.onset is None:
        ok = rep.check("MAPK output responds before stress/apoptosis", False, "MAPK never active.")
    else:
        detail, ok = [], True
        for s in STRESS:
            o = ons[s]
            if o.onset is None:
                detail.append(f"{s}: never active")
            else:
                later = o.onset[0] >= d.onset[1]
                ok &= later
                detail.append(f"{s}: onset ({o.onset[0]:g}, {o.onset[1]:g}] h")
        ok = rep.check(
            "MAPK output responds before stress/apoptosis",
            ok,
            f"MAPK onset ({d.onset[0]:g}, {d.onset[1]:g}] h, "
            f"{'down' if d.direction < 0 else 'up'}; " + "; ".join(detail) + ".",
        )
    rep.text("", f"**Overall: {'PASS' if ok else 'FAIL'}**")
    return rep.write(), ok


if __name__ == "__main__":
    path, ok = run()
    print(path, "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
