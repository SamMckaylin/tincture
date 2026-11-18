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

Uncertainty (version 1, the first pre-registration): each matrix is a single replicate, so
significance uses a parametric bootstrap with the well-level noise SD set to the root mean
square of GDSC's two single-agent fit RMSEs for that cell line. Region = BH q < 0.05 and
excess ≥ 10%. Result of version 1: check 1 FAILED (2 of 6), check 2 passed (9 of 10).

**Version 2 (pre-registered 2026-11-18, after version 1 failed, chosen by the user as
option b; written before any version-2 run).** Everything as above except the noise SD,
which is estimated from each matrix itself (``stats.synergy.matrix_noise_sd``): second
differences y₋₁ − 2y₀ + y₊₁ along rows and columns of the 7×7 grid, over triplets whose three
viabilities all lie in (0.1, 0.9), robust scale 1.4826·median|Δ²|/√6; with fewer than 10
such triplets no significance is claimed. Same pass thresholds as version 1 (≥ 4 of 6;
≥ 80% clean).
5. **Calibration (gate for version 2):** on every cell line × drug pair screened as ≥ 3
   replicate matrices in the same screens (gdsc-008, gdsc-010), compare the single-matrix
   estimate (median over its replicates) with the empirical replicate SD (pooled across
   combination wells whose mean viability lies in (0.1, 0.9); ≥ 5 such wells). Version 2 is
   valid only if the median ratio estimate / replicate SD is ≥ 0.8 (not anti-conservative);
   otherwise its checks are reported as invalid, whatever they show.

Figures (version 2): ML-2 for both combinations (named in the paper; screened with both).

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
from tincture.io.gdsc_combo import bashi_results, load_bashi_matrices, load_bashi_replicated
from tincture.stats.synergy import NOISE_RANGE, matrix_noise_sd, synergy

BENEFICIAL = ("Selumetinib", "Venetoclax")
NAMED = ["EoL-1-cell", "ML-2", "OCI-AML5", "NOMO-1", "KG-1", "HL-60"]
NOT_BENEFICIAL = ("Olaparib", "AZD5153")
FIGURE_LINE = "ML-2"
PE = "HSA ≥ 0.1 & viab ≤ 0.5, point est. (cells)"
FOOT = ("GDSC matrix screens, Bashi et al., Cancer Discov 2024 (internal research use)",)


def evaluate(mats: dict, lines: list[str], version: int) -> pd.DataFrame:
    rows = []
    for line in lines:
        cm = mats[line]
        est = matrix_noise_sd(cm) if version == 2 else float(cm.meta["noise_sd"])
        sd = est if est is not None else 0.0
        b = synergy(cm, "bliss", noise_sd=sd or None)
        h = synergy(cm, "hsa", noise_sd=sd or None)
        g = cm.meta["gdsc_summary"]
        rows.append(
            {
                "cell line": line,
                "noise SD": round(sd, 3) if sd else "n/a",
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
    return pd.DataFrame(rows)


def calibration() -> tuple[pd.DataFrame, float]:
    rows = []
    for cm in load_bashi_replicated(("gdsc-008", "gdsc-010"), min_replicates=3):
        c = cm.combo
        w = c.groupby(["conc_a", "conc_b"])["viability"].agg(["mean", "var", "size"])
        w = w[(w["size"] >= 3) & (w["mean"] > NOISE_RANGE[0]) & (w["mean"] < NOISE_RANGE[1])]
        if len(w) < 5:
            continue
        est = matrix_noise_sd(cm)
        if est is None:
            continue
        rep_sd = float(np.sqrt(w["var"].mean()))
        rows.append(
            {"cell line": cm.cell_line, "pair": f"{cm.drug_a} + {cm.drug_b}",
             "replicates": len(cm.replicates), "estimate": est, "replicate SD": rep_sd,
             "ratio": est / rep_sd}
        )  # fmt: skip
    t = pd.DataFrame(rows)
    return t, float(t["ratio"].median())


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

    m1 = load_bashi_matrices(*BENEFICIAL, NAMED)
    m2 = load_bashi_matrices(*NOT_BENEFICIAL, nb_lines)

    rep.text("## Version 1 (first pre-registration: noise = GDSC single-agent fit RMSE)")
    v1a, v1b = evaluate(m1, NAMED, 1), evaluate(m2, nb_lines, 1)
    rep.text("### Selumetinib + venetoclax", md_table(v1a))
    rep.text("### Olaparib + AZD5153", md_table(v1b))
    k1a = int((v1a["Bliss region (cells)"] > 0).sum())
    k1b = int((v1b["Bliss region (cells)"] == 0).sum())
    rep.text(
        f"Version 1: beneficial region in {k1a} of {len(v1a)} (needs ≥ 4): "
        f"{'pass' if k1a >= 4 else 'FAIL'}; non-beneficial clean in {k1b} of {len(v1b)}: "
        f"{'pass' if k1b >= 0.8 * len(v1b) else 'FAIL'}."
    )

    rep.text("## Version 2 calibration (gate)")
    cal, ratio = calibration()
    rep.text(
        f"{len(cal)} replicated matrices (≥ 3 replicates, gdsc-008 and gdsc-010). Ratio of "
        f"the single-matrix estimate to the empirical replicate SD: median {ratio:.2f} "
        f"(IQR {cal['ratio'].quantile(0.25):.2f}–{cal['ratio'].quantile(0.75):.2f}); "
        f"estimate median {cal['estimate'].median():.3f}, replicate SD median "
        f"{cal['replicate SD'].median():.3f}."
    )
    valid = ratio >= 0.8

    rep.text("## Version 2 (noise estimated from each matrix)")
    t1, t2 = evaluate(m1, NAMED, 2), evaluate(m2, nb_lines, 2)
    rep.text("### 1. Selumetinib + venetoclax (reported beneficial in AML)", md_table(t1))
    rep.text("### 2. Olaparib + AZD5153 (reported no benefit)", md_table(t2))

    for (d1, d2), mats in ((BENEFICIAL, m1), (NOT_BENEFICIAL, m2)):
        cm = mats[FIGURE_LINE]
        for model in ("bliss", "zip"):
            r = tc.combination_surface(
                cm,
                model=model,
                noise_sd="matrix",
                title=f"Combination Surface ({'Bliss' if model == 'bliss' else 'ZIP'}): "
                f"{cm.drug_a} + {cm.drug_b}, {FIGURE_LINE} (AML)",
                footer=FOOT,
            )
            tag = "" if model == "bliss" else "zip_"
            rep.figure(f"surface_{tag}{d1.lower()}_{d2.lower()}_{FIGURE_LINE}", r)

    both = pd.concat([t1, t2])
    rho = float(both[["mean Bliss excess", "GDSC Bliss_matrix"]].corr(method="spearman").iloc[0, 1])
    rep.text("## Ground-truth checks")
    rep.check(
        "Version 2 calibration gate (median ratio ≥ 0.8)",
        valid,
        f"median ratio {ratio:.2f} over {len(cal)} replicated matrices"
        + ("" if valid else "; version 2 checks are invalid"),
    )
    k1 = int((t1["Bliss region (cells)"] > 0).sum())
    ok1 = rep.check(
        "Version 2, beneficial combination shows a significant excess region",
        valid and k1 >= 4,
        f"selumetinib + venetoclax: Bliss region in {k1} of {len(t1)} named AML lines "
        f"(HSA benefit region in {int((t1['HSA region (cells)'] > 0).sum())} of {len(t1)}).",
    )
    k2 = int((t2["Bliss region (cells)"] == 0).sum())
    ok2 = rep.check(
        "Version 2, non-beneficial combination shows no region",
        valid and k2 >= 0.8 * len(t2),
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
    rep.text(
        "",
        f"**Overall (version 2, checks 1–2 with the calibration gate): "
        f"{'PASS' if ok else 'FAIL'}**. Version 1 result is kept above.",
    )
    _ = np
    return rep.write(), ok


if __name__ == "__main__":
    path, ok = run()
    print(path, "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
