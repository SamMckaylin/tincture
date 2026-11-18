"""CS5: Combination Surface on GDSC full 7×7 matrices (HANDOVER §6.2; D5, Bashi et al. 2024).

Pre-specified (fixed before the first run, from the paper's text):
1. **Beneficial:** selumetinib + venetoclax in AML. The paper names six AML lines with strong
   combination activity: EoL-1-cell, ML-2, OCI-AML5, NOMO-1, KG-1, HL-60. Pass: a significant
   excess region (Bliss, the plot's default) appears in ≥ 4 of the 6.
2. **Not beneficial:** olaparib + AZD5153, one of six combinations with "no combination
   benefit (HSA < 0.1, combination Emax < 0.5) in the 755 cell lines tested". Evaluated in
   every AML line screened with it. Pass: no significant region in ≥ 80% of them.
3. Secondary (reported, not pass/fail): the same counts for benefit over the best single
   agent (HSA, the paper's own metric), and agreement of our mean Bliss excess with GDSC's
   ``Bliss_matrix`` across lines.
4. Diagnostic added after the first run (report only, never pass/fail): the number of cells
   whose *point estimate* of HSA excess is ≥ 0.1 with observed viability ≤ 0.5, i.e. the
   paper's activity rule applied per cell without any significance test.

Uncertainty: each matrix is a single replicate, so significance uses a parametric bootstrap
with the well-level noise SD set to the root mean square of GDSC's two single-agent fit
RMSEs for that cell line (stated per figure). Region = BH q < 0.05 and excess ≥ 10%.
Figures: ML-2 for both combinations (named in the paper; screened with both).

GDSC data: figures and report stay local (internal research use).

    python case_studies/cs5_gdsc_matrix.py
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
from tincture.io.gdsc_combo import bashi_results, load_bashi_matrices
from tincture.stats.synergy import synergy

BENEFICIAL = ("Selumetinib", "Venetoclax")
NAMED = ["EoL-1-cell", "ML-2", "OCI-AML5", "NOMO-1", "KG-1", "HL-60"]
NOT_BENEFICIAL = ("Olaparib", "AZD5153")
FIGURE_LINE = "ML-2"
PE = "HSA ≥ 0.1 & viab ≤ 0.5, point est. (cells)"
FOOT = ("GDSC matrix screens, Bashi et al., Cancer Discov 2024 (internal research use)",)


def evaluate(drug_1: str, drug_2: str, lines: list[str]) -> tuple[pd.DataFrame, dict]:
    mats = load_bashi_matrices(drug_1, drug_2, lines)
    rows = []
    for line in lines:
        cm = mats[line]
        sd = float(cm.meta["noise_sd"])
        b = synergy(cm, "bliss", noise_sd=sd)
        h = synergy(cm, "hsa", noise_sd=sd)
        g = cm.meta["gdsc_summary"]
        rows.append(
            {
                "cell line": line,
                "noise SD": round(sd, 3),
                "Bliss region (cells)": b.n_significant,
                "masked": int(b.masked.sum()),
                "mean Bliss excess": round(b.score, 3),
                "GDSC Bliss_matrix": round(float(g["Bliss_matrix"].mean()), 3),
                "HSA region (cells)": h.n_significant,
                PE: int(((h.excess >= 0.1) & (h.observed <= 0.5) & ~h.masked).sum()),
                "GDSC HSA_matrix": round(float(g["HSA_matrix"].mean()), 3),
                "GDSC combo_MaxE": round(float(g["combo_MaxE"].mean()), 2),
            }
        )
    return pd.DataFrame(rows), mats


def run() -> tuple[Path, bool]:
    rep = Report("cs5_gdsc_matrix", "CS5: Combination Surface on GDSC 7×7 matrices")
    rep.text(
        "Pre-specified in the script header. Data: Bashi et al. 2024 per-well matrix results "
        "(figshare project 163378). Default model: Bliss independence; benefit over the best "
        "single agent (HSA) reported alongside."
    )
    res = bashi_results()
    aml = res[res["CANCER_TYPE"] == "Acute Myeloid Leukemia"]
    nb = aml[(aml["lib1_name"].isin(NOT_BENEFICIAL)) & (aml["lib2_name"].isin(NOT_BENEFICIAL))]
    nb_lines = sorted(nb["CELL_LINE_NAME"].unique())

    t1, m1 = evaluate(*BENEFICIAL, NAMED)
    rep.text("## 1. Selumetinib + venetoclax (reported beneficial in AML)", md_table(t1))
    t2, m2 = evaluate(*NOT_BENEFICIAL, nb_lines)
    rep.text("## 2. Olaparib + AZD5153 (reported no benefit)", md_table(t2))

    for (d1, d2), mats in ((BENEFICIAL, m1), (NOT_BENEFICIAL, m2)):
        cm = mats[FIGURE_LINE]
        r = tc.combination_surface(
            cm,
            noise_sd=float(cm.meta["noise_sd"]),
            title=f"Combination Surface: {cm.drug_a} + {cm.drug_b}, {FIGURE_LINE} (AML)",
            footer=FOOT,
        )
        rep.figure(f"surface_{d1.lower()}_{d2.lower()}_{FIGURE_LINE}", r)
        z = tc.combination_surface(
            cm,
            model="zip",
            noise_sd=float(cm.meta["noise_sd"]),
            title=f"Combination Surface (ZIP): {cm.drug_a} + {cm.drug_b}, {FIGURE_LINE}",
            footer=FOOT,
        )
        rep.figure(f"surface_zip_{d1.lower()}_{d2.lower()}_{FIGURE_LINE}", z)

    both = pd.concat([t1, t2])
    rho = float(both[["mean Bliss excess", "GDSC Bliss_matrix"]].corr(method="spearman").iloc[0, 1])
    rep.text("## Ground-truth checks")
    k1 = int((t1["Bliss region (cells)"] > 0).sum())
    ok1 = rep.check(
        "Beneficial combination shows a significant excess region",
        k1 >= 4,
        f"selumetinib + venetoclax: Bliss region in {k1} of {len(t1)} named AML lines "
        f"(HSA benefit region in {int((t1['HSA region (cells)'] > 0).sum())} of {len(t1)}).",
    )
    k2 = int((t2["Bliss region (cells)"] == 0).sum())
    ok2 = rep.check(
        "Non-beneficial combination shows no region",
        k2 >= 0.8 * len(t2),
        f"olaparib + AZD5153: no Bliss region in {k2} of {len(t2)} AML lines "
        f"(no HSA region in {int((t2['HSA region (cells)'] == 0).sum())} of {len(t2)}).",
    )
    rep.text(
        f"Secondary: Spearman correlation between our mean Bliss excess and GDSC's "
        f"Bliss_matrix across the {len(both)} matrices: ρ = {rho:.2f}."
    )
    rep.text(
        "Diagnostic (added after the first run, not pass/fail): with the paper's point-"
        "estimate rule per cell (HSA excess ≥ 0.1 and viability ≤ 0.5, no significance), "
        f"{int((t1[PE] > 0).sum())} of {len(t1)} beneficial-combination lines and "
        f"{int((t2[PE] > 0).sum())} of {len(t2)} non-beneficial lines have at least "
        "one such cell. Single 7×7 replicates with fit RMSE 0.1–0.2 cannot make most of these "
        "excesses significant."
    )
    ok = ok1 and ok2
    rep.text("", f"**Overall (pre-specified checks 1–2): {'PASS' if ok else 'FAIL'}**")
    _ = np
    return rep.write(), ok


if __name__ == "__main__":
    path, ok = run()
    print(path, "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
