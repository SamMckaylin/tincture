"""CS4: Combination Lift on GDSC anchored screens (HANDOVER §6.2; D4 Jaaks et al. 2022, D6
Vis et al. 2024).

Pre-specified (fixed before the first run, from Jaaks et al. 2022, where TOP1 inhibition was
screened as camptothecin in colon and the CHEK1/2 inhibitor was AZD7762):
* Data: original screen, colon, anchor AZD7762 (0.0625 and 0.25 µM) × library camptothecin;
  per cell line, mean over replicate plates of the viabilities at the highest library
  concentration (anchor alone, library Emax, combination Emax).
* Benefit: combination effect ≥ 20 percentage points above the best single agent.
* Groups: microsatellite status from Cell Model Passports (MSS vs MSI).
1. Pass if ≥ 20% of colon lines show benefit at one or both anchor concentrations.
2. Pass if the benefit rate is higher in MSS than MSI lines at both anchor concentrations
   (direction pre-specified; Fisher and Mann–Whitney p reported, not thresholded, because
   the groups are small: ~30 vs 15 lines).
Reported alongside (not pass/fail): Jaaks's own Bliss-based synergy calls (``Synergy``) per
group, and a second dataset (D6, Vis et al. 2024): trametinib + pictilisib, which the paper
highlights in lung (ECB rate 47%), lung vs other tissues.

GDSC data: figures and report stay local (internal research use).

    python case_studies/cs4_gdsc_anchored.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd
from _report import Report, md_table

import tincture as tc
from tincture.io.gdsc_combo import (
    anchored_lift_table,
    cell_model_annotations,
    load_jaaks2022,
    load_vis2024,
)

ANCHOR, LIBRARY = "AZD7762", "Camptothecin"
FOOT = ("GDSC anchored screens, Jaaks et al., Nature 2022 (internal research use)",)
BLISS_NOTE = (
    "Synergy in the source paper is a separate, Bliss-based call (ΔEmax ≥ 0.2 or 8-fold "
    "ΔIC50 over Bliss); the Lift shows benefit over the best single agent."
)


def run() -> tuple[Path, bool]:
    rep = Report("cs4_gdsc_anchored", "CS4: Combination Lift on GDSC anchored screens")
    rep.text("Pre-specified in the script header.")
    j = load_jaaks2022("original")
    col = j[j["Tissue"] == "Colon"]
    lt = anchored_lift_table(col, ANCHOR, LIBRARY)
    cmp_ = cell_model_annotations().set_index("model_id")["msi_status"]
    lt["MSI status"] = lt["sidm"].map(cmp_)
    n_missing = int(lt["MSI status"].isna().sum())
    lt = lt.dropna(subset=["MSI status"])
    res = tc.combination_lift(
        lt,
        ANCHOR,
        LIBRARY,
        group="MSI status",
        facet="anchor",
        model_note=BLISS_NOTE,
        title=f"Combination Lift: {ANCHOR} (anchor) + {LIBRARY}, colon (Jaaks 2022)",
        footer=FOOT,
    )
    rep.text(f"## 1. {ANCHOR} + {LIBRARY} in colon")
    if n_missing:
        rep.text(f"{n_missing} rows without MSI annotation dropped.")
    rep.figure("lift_azd7762_camptothecin_colon", res)
    st = res.meta["stats"]
    rows = []
    for facet, s in st.items():
        rows.append(
            {
                "anchor": facet,
                "lines": s["n"],
                "benefit (all)": f"{s['k']}/{s['n']} ({s['k'] / s['n']:.0%})",
                "MSS": f"{s['k_MSS']}/{s['n_MSS']} ({s['rate_MSS']:.0%})",
                "MSI": f"{s['k_MSI']}/{s['n_MSI']} ({s['rate_MSI']:.0%})",
                "Fisher p": f"{s['fisher_p']:.3g}",
                "Mann–Whitney p (ΔEmax)": f"{s['mannwhitney_p']:.3g}",
            }
        )
    rep.text(md_table(pd.DataFrame(rows)))
    # Jaaks's own synergy calls, per group (anchor orientation as above)
    m = col[(col["ANCHOR_NAME"] == ANCHOR) & (col["LIBRARY_NAME"] == LIBRARY)].copy()
    m["MSI status"] = m["SIDM"].map(cmp_)
    syn = (
        m.groupby(["ANCHOR_CONC", "MSI status", "CELL_LINE_NAME"])["Synergy"]
        .mean()
        .ge(0.5)
        .groupby(level=[0, 1])
        .agg(["sum", "size"])
        .reset_index()
    )
    syn["rate"] = (syn["sum"] / syn["size"]).map("{:.0%}".format)
    rep.text("Source-paper synergy calls (≥ half of replicates synergistic):", md_table(syn))

    # second dataset: Vis 2024
    v = load_vis2024()
    vt = anchored_lift_table(v, "Trametinib", "Pictilisib")
    tissue = v.drop_duplicates("CELL_LINE_NAME").set_index("CELL_LINE_NAME")["tissue"]
    vt["tissue"] = vt["cell_line"].map(tissue).fillna("unknown")
    vt["lung"] = vt["tissue"].str.contains("Lung", case=False).map({True: "lung", False: "other"})
    vres = tc.combination_lift(
        vt,
        "Trametinib",
        "Pictilisib",
        group="lung",
        facet="anchor",
        model_note=BLISS_NOTE,
        title="Combination Lift: trametinib (anchor) + pictilisib, pan-cancer (Vis 2024)",
        footer=("GDSC pan-cancer anchored screen, Vis et al., Cell Rep Med 2024",),
    )
    rep.text("## 2. Second dataset: trametinib + pictilisib (Vis 2024), lung vs other")
    rep.figure("lift_trametinib_pictilisib_pancancer", vres)

    rep.text("## Ground-truth checks")
    ok1 = rep.check(
        "Benefit in colon",
        any(s["k"] / s["n"] >= 0.2 for s in st.values()),
        "; ".join(f"{f}: {s['k']}/{s['n']}" for f, s in st.items()),
    )
    ok2 = rep.check(
        "Benefit enriched in MSS vs MSI",
        all(s["rate_MSS"] > s["rate_MSI"] for s in st.values()),
        "; ".join(
            f"{f}: MSS {s['rate_MSS']:.0%} vs MSI {s['rate_MSI']:.0%} (Fisher p "
            f"{s['fisher_p']:.2g})"
            for f, s in st.items()
        ),
    )
    ok = ok1 and ok2
    rep.text("", f"**Overall (pre-specified checks 1–2): {'PASS' if ok else 'FAIL'}**")
    return rep.write(), ok


if __name__ == "__main__":
    path, ok = run()
    print(path, "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
