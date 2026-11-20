"""M2 acceptance: the Retrievability gate works on sci-Plex3 (D1).

Checks (written to case_studies/out/m2_retrievability_sciplex/report.md):
1. Calibration — vehicle wells held out as pseudo-conditions (no phenotype) are called
   at about the nominal rate (p < 0.05 in ≈ 5%).
2. Dose — median mAP rises with dose in every cell line.
3. Biology — HDAC inhibitors (a strong, known transcriptional class in sci-Plex) are
   retrievable at ≥ 1 µM more often than other compounds.

    python case_studies/m2_retrievability_sciplex.py
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
from tincture.io.sciplex import load_sciplex3
from tincture.stats.map_retrieval import retrievability
from tincture.stats.signatures import replicate_profiles

FOCUS = "S1090|1e-06"  # abexinostat 1 µM


def calibration(prof, n_random: int = 60, seed: int = 0) -> np.ndarray:
    o = prof.obs
    X = np.asarray(prof.X)
    veh = np.flatnonzero(o["dose"].to_numpy() == 0)
    rng = np.random.default_rng(seed)
    ps = []
    for _ in range(n_random):
        a, b = rng.choice(veh, 2, replace=False)
        g = np.array(["ctrl"] * len(o), dtype=object)
        g[[a, b]] = "pseudo"
        ctrl = np.zeros(len(o), bool)
        ctrl[veh] = True
        ctrl[[a, b]] = False
        keep = ctrl | (g == "pseudo")
        ps.append(float(retrievability(X[keep], g[keep], ctrl[keep]).table["p"].iloc[0]))
    return np.array(ps)


def run() -> tuple[Path, bool]:
    rep = Report("m2_retrievability_sciplex", "M2: Retrievability gate on sci-Plex3 (24 h)")
    rep.text(
        "Profiles: per-well logCPM minus same-line vehicle mean, scaled by vehicle SD; 2,000 "
        "highest-variance genes chosen without using condition labels. "
        "Conditions = compound × dose "
        "(2 wells each); negatives = the line's 32 vehicle wells; permutation null (exact)."
    )
    rows, cal_rows = [], []
    dose_ok: dict[str, bool] = {}
    hdac_ok: dict[str, bool] = {}
    cal_ok: dict[str, bool] = {}
    for line in ["A549", "MCF7", "K562"]:
        prof = replicate_profiles(load_sciplex3(cell_lines=[line]))
        res = tc.retrievability(
            prof,
            focus=[FOCUS],
            title=f"Retrievability: sci-Plex3 {line}, 24 h",
            footer=("sci-Plex3, Srivatsan et al. Science 2020 (GSE139944)",),
        )
        rep.text(f"## {line}")
        rep.figure(f"retrievability_{line}", res)
        t = res.meta["table"].copy()
        t["compound_id"] = t["group"].str.split("|").str[0]
        t["dose"] = t["group"].str.split("|").str[1].astype(float)
        tgt = prof.obs.drop_duplicates("compound_id").set_index("compound_id")["target"].astype(str)
        t["hdac"] = t["compound_id"].map(tgt).str.contains("HDAC")
        by_dose = t.groupby("dose")["map"].median()
        hi = t[t["dose"] >= 1e-6]
        rows.append(
            {
                "cell line": line,
                "retrievable": f"{int(t['retrievable'].sum())}/{len(t)}",
                **{f"median mAP {tc.report.fmt_conc(d)}": round(v, 2) for d, v in by_dose.items()},
                "HDACi ≥1 µM": f"{hi[hi['hdac']]['retrievable'].mean():.0%}",
                "others ≥1 µM": f"{hi[~hi['hdac']]['retrievable'].mean():.0%}",
            }
        )
        p = calibration(prof)
        cal_rows.append(
            {
                "cell line": line,
                "pseudo-conditions": len(p),
                "p < 0.05": f"{(p < 0.05).mean():.0%}",
                "median p": round(float(np.median(p)), 2),
            }
        )
        dose_ok[line] = bool(
            np.all(np.diff(by_dose.to_numpy()) > -0.02) and by_dose.iloc[-1] > by_dose.iloc[0]
        )
        hdac_ok[line] = bool(
            hi[hi["hdac"]]["retrievable"].mean() > hi[~hi["hdac"]]["retrievable"].mean()
        )
        cal_ok[line] = bool((p < 0.05).mean() <= 0.15)
    summary = pd.DataFrame(rows)
    cal = pd.DataFrame(cal_rows)
    rep.text(
        "## Summary",
        md_table(summary),
        "## Calibration: vehicle wells as pseudo-conditions",
        md_table(cal),
        "## Checks",
    )

    def fmt(d: dict[str, bool]) -> str:
        return ", ".join(f"{k}: {'yes' if v else 'NO'}" for k, v in d.items())

    ok = rep.check(
        "Calibrated null (≤ 15% of vehicle pseudo-conditions at p < 0.05)",
        all(cal_ok.values()),
        fmt(cal_ok),
    )
    ok &= rep.check("Median mAP rises with dose", all(dose_ok.values()), fmt(dose_ok))
    ok &= rep.check(
        "HDAC inhibitors retrievable more often at ≥ 1 µM", all(hdac_ok.values()), fmt(hdac_ok)
    )
    path = rep.write()
    return path, bool(ok)


if __name__ == "__main__":
    p, ok = run()
    print(p, "ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
