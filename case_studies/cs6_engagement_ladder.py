"""CS6: Engagement Ladder on real kinase inhibitors (HANDOVER §6.2).

Data: Davis 2011 Kd (D8, biochemical), NanoBRET occupancy in HEK293 (D7, cellular; apparent
IC50 from single-dose occupancy), PRISM viability (D2, phenotype; our 4PL fit).

Pre-specified (fixed before the first run):
1. **Selectivity (Vasta et al. 2018):** among kinases in both panels, the number engaged
   > 50% in cells (NanoBRET at 1 µM crizotinib / 100 nM dasatinib) vs the number predicted
   from Kd (Kd ≤ dose ⇔ > 50% occupancy in the ATP-free binding assay). Pass: crizotinib
   engages fewer kinases in cells than predicted biochemically (paper: 13 vs 69).
2. **Ladders** (intended target; phenotype line driven by it): dasatinib→ABL1 / K562
   (BCR-ABL), crizotinib→ALK / KARPAS-299 (NPM-ALK), crizotinib→MET / MKN45 (MET-amplified).
   Pass: at least one gets the off-target flag (phenotype > 10× more potent than cellular
   engagement).
3. The ATP-shift (purple) flag needs kinase ATP Km values, which none of the data sources
   provides. **Update (2026-11-18, user decision):** Km values found and verified in three
   open-access papers (docs/data_status.md) are now applied, labelled per source. Rule, fixed
   before rerunning: where two enzyme forms are reported (ALK, MET), the wedge spans both
   (``atp_shift`` with a Km range), so the flag fires only beyond the widest expectation.
   Reported, not pass/fail: the Km values come from different constructs and assays.

Phenotype potency: absolute IC50 (50% viability) from a 4PL fit of PRISM replicate wells
with 200 bootstrap refits; a curve already below 50% at the lowest dose is censored "<",
one never reaching 50% is censored ">", and a fit with R² < 0.8 is not used (rung skipped).

    python case_studies/cs6_engagement_ladder.py
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

import tincture as tc
from tincture.io.davis import load_davis2011
from tincture.io.kinase_names import to_davis
from tincture.io.nanobret import DOSES_M, load_nanobret, occupancy
from tincture.io.prism import load_prism_secondary
from tincture.model.potency import PotencyTable
from tincture.stats.hill import bootstrap_4pl, fit_4pl

PAIRS = [
    ("dasatinib", "ABL1", "ACH-000551", "K562"),
    ("crizotinib", "ALK", "ACH-000053", "KARPAS-299"),
    ("crizotinib", "MET", "ACH-000356", "MKN45"),
]
PAPER = {"crizotinib": (13, 69), "dasatinib": (22, 44)}
#: ATP Km (molar), verified against the papers' full text (docs/data_status.md).
KM_ATP: dict[str, float | tuple[float, float]] = {
    "ABL1": 9.78e-6,  # Grover et al., PLoS One 2015 (ABL core, ADP Quest)
    "ALK": (1.34e-4, 1.59e-4),  # Bresler et al., Sci Transl Med 2011 (unphospho / phospho)
    "MET": (6.0e-5, 2.0e-4),  # Behshad et al., Curr Chem Genomics 2010 (phospho / unphospho)
}
KM_SOURCES = (
    "ATP Km: ABL1 9.8 µM (Grover 2015); ALK 134–159 µM (Bresler 2011); MET 60–200 µM "
    "(Behshad 2010); different constructs and assays"
)
FOOT = (
    "Biochemical: Davis et al. 2011 (Kd); cellular: Vasta et al. 2018 (NanoBRET, HEK293); "
    "phenotype: PRISM 20Q2 secondary (CC BY 4.0)",
)


def selectivity_counts(dv: pd.DataFrame) -> pd.DataFrame:
    occ = occupancy().assign(target=lambda x: x["kinase"].map(to_davis))
    rows = []
    for drug, dose in DOSES_M.items():
        o = occ[occ["compound_id"] == drug]
        d = dv[dv["compound_id"] == drug].set_index("target")
        shared = [k for k in o["target"] if k in d.index]
        o = o[o["target"].isin(shared)]
        kd = d.loc[shared]
        cell = int((o["mean"] > 0.5).sum())
        bio = int(((kd["kd_nm"] <= dose * 1e9 * (1 + 1e-9)) & (kd["censored"] == "")).sum())
        rows.append(
            {
                "compound": drug,
                "dose": tc.report.fmt_conc(dose),
                "kinases": len(shared),
                "engaged in cells (>50%)": cell,
                "predicted from Kd": bio,
                "paper (cells / Kd)": f"{PAPER[drug][0]} / {PAPER[drug][1]}",
            }
        )
    return pd.DataFrame(rows)


def phenotype_rung(drug: str, line: str, label: str) -> dict | None:
    drs = load_prism_secondary(compounds=[drug])
    dr = drs.get(drs.keys()["compound_id"].iloc[0], line)
    h = fit_4pl(dr.dose, dr.response)
    if not h.converged or not np.isfinite(h.r2) or h.r2 < 0.8:
        return None
    dmin, dmax = dr.dose_range
    row = {"rung": "phenotype", "assay_type": "IC50", "cell_system": label}
    if h.fit.top < 0.5:  # below 50% already at the lowest dose
        return {**row, "px": -np.log10(dmin), "censored": "<"}
    ic = h.fit.inverse(0.5)
    if not np.isfinite(ic) or ic > dmax:
        return {**row, "px": -np.log10(dmax), "censored": ">"}
    b = bootstrap_4pl(dr.dose, dr.response, n_boot=200)
    lo, hi = b.interval(lambda f: f.inverse(0.5))
    return {
        **row,
        "px": -np.log10(ic),
        "px_lo": -np.log10(hi),
        "px_hi": -np.log10(lo),
        "censored": "",
    }


def run() -> tuple[Path, bool]:
    rep = Report("cs6_engagement_ladder", "CS6: Engagement Ladder on real kinase inhibitors")
    rep.text(
        "Pre-specified in the script header. Data: Davis 2011 Kd (biochemical), "
        "NanoBRET single-dose occupancy in HEK293 (cellular; apparent IC50 = C·(1−f)/f, "
        "1:1 binding), PRISM viability (phenotype)."
    )
    dv = load_davis2011(["Crizotinib", "Dasatinib"]).data
    dv["compound_id"] = dv["compound_id"].str.lower()
    nb = load_nanobret().data
    # 1. selectivity
    counts = selectivity_counts(dv)
    rep.text("## 1. Selectivity: biochemical vs live-cell engagement", md_table(counts))
    for drug, panel in (("crizotinib", dv), ("crizotinib", nb)):
        tab = PotencyTable(panel[panel["compound_id"] == drug].copy())
        kind = "biochemical Kd (Davis 2011)" if panel is dv else "cellular apparent IC50 (NanoBRET)"
        res = tc.selectivity(
            tab,
            drug,
            ["ALK", "MET", "ROS1"] if panel is dv else ["ALK", "MET"],
            panel_size=518,
            panel_name="kinases",
            title=f"Selectivity Spectrum: crizotinib, {kind}",
            footer=FOOT,
        )
        rep.figure(f"selectivity_crizotinib_{'davis' if panel is dv else 'nanobret'}", res)
    # 2. ladders
    rows, pairs, skipped = [], [], []
    for drug, tgt, line, label in PAIRS:
        b = dv[(dv["compound_id"] == drug) & (dv["target"] == tgt)].iloc[0]
        rows.append(
            {
                "compound_id": drug,
                "target": tgt,
                "rung": "biochemical",
                "px": b["px"],
                "assay_type": "Kd",
                "censored": b["censored"],
                "cell_system": "cell-free (KINOMEscan)",
            }
        )
        c = nb[(nb["compound_id"] == drug) & (nb["target"] == tgt)].iloc[0]
        rows.append(
            {
                "compound_id": drug,
                "target": tgt,
                "rung": "cellular",
                "px": c["px"],
                "px_lo": c["px_lo"],
                "px_hi": c["px_hi"],
                "assay_type": "IC50_cell",
                "censored": c["censored"],
                "cell_system": "HEK293",
            }
        )
        ph = phenotype_rung(drug, line, label)
        if ph is None:
            skipped.append(f"{drug}/{label}")
        else:
            rows.append({"compound_id": drug, "target": tgt, **ph})
        pairs.append((drug, tgt))
    lad = tc.engagement_ladder(
        PotencyTable(pd.DataFrame(rows)),
        pairs,
        km_atp=KM_ATP,
        title="Engagement Ladder: dasatinib and crizotinib",
        footer=(*FOOT, KM_SOURCES),
    )
    rep.text("## 2. Engagement ladders")
    if skipped:
        rep.text(f"Phenotype rung skipped (fit R² < 0.8): {', '.join(skipped)}.")
    rep.figure("engagement_ladder", lad)
    results = lad.meta["results"]
    flag_tab = pd.DataFrame(
        [
            {
                "compound / target": f"{r.compound} / {r.target}",
                "shift biochemical→cellular (log)": f"{r.shift_bc_bound}{r.shift_bc:.1f}"
                if r.shift_bc is not None
                else "—",
                "off-target flag": "yes" if "off_target" in r.flags else "no",
                "expected ATP shift (log)": f"{r.atp_shift[0]:.1f}–{r.atp_shift[1]:.1f}"
                if r.atp_shift
                else "—",
                "ATP-shift flag": "not evaluated (no ATP Km)"
                if "shift" in r.not_evaluated
                else ("yes" if "shift" in r.flags else "no"),
            }
            for r in results
        ]
    )
    rep.text(md_table(flag_tab))
    rep.text(
        "**Reading the dasatinib flag.** K562 is driven by BCR-ABL, so its phenotype is "
        "expected to be on-target. The flag compares it with NanoBRET engagement of "
        "full-length ABL1 in HEK293 cells: a different protein context and cell system, "
        "which the ladder shows as a dashed segment and labels on each rung. The flag is "
        "doing its job (phenotype far more potent than the measured engagement), and the "
        "likely explanation here is the engagement system, not an off-target. That is what "
        "the flag text means by 'likely': it asks for a check, not a conclusion."
    )
    # checks
    crz = counts.set_index("compound").loc["crizotinib"]
    rep.text("## Ground-truth checks")
    ok1 = rep.check(
        "Crizotinib shows improved cellular selectivity vs biochemical",
        bool(crz["engaged in cells (>50%)"] < crz["predicted from Kd"]),
        f"{crz['engaged in cells (>50%)']} of {crz['kinases']} kinases engaged > 50% in cells "
        f"at 1 µM vs {crz['predicted from Kd']} predicted from Kd (paper: 13 vs 69).",
    )
    flagged = [f"{r.compound}/{r.target}" for r in results if "off_target" in r.flags]
    ok2 = rep.check(
        "At least one compound gets the off-target flag",
        bool(flagged),
        f"flagged: {', '.join(flagged) or 'none'} (phenotype > 10× more potent than cellular "
        "engagement of the intended target).",
    )
    shifted = [f"{r.compound}/{r.target}" for r in results if "shift" in r.flags]
    rep.text(
        "**ATP-shift (purple) flag, reported, not pass/fail:** "
        + (f"fires for {', '.join(shifted)}." if shifted else "fires for none.")
        + " Caveats: the three Km values come from different labs, constructs (ABL core "
        "with SH3/SH2; ALK and MET kinase domains) and assay formats; the cellular rung is "
        "an apparent IC50 from single-dose occupancy; and the biochemical rung is a binding "
        "Kd for an isolated kinase domain. A shift beyond the ATP wedge is therefore a prompt "
        "to check permeability, efflux or the construct difference, not a conclusion."
    )
    ok = ok1 and ok2
    rep.text("", f"**Overall (pre-specified checks 1–2): {'PASS' if ok else 'FAIL'}**")
    return rep.write(), ok


if __name__ == "__main__":
    path, ok = run()
    print(path, "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
