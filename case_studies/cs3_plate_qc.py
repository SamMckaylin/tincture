"""CS3 — Plate QC on GDSC2 raw monotherapy data (release 8.5).

Ground-truth check (HANDOVER §6.2): the Z′ distribution across plates is shown, the
edge-vs-interior ratio is reported, and failing plates are flagged.

    python case_studies/cs3_plate_qc.py [--max-plates N]

Writes case_studies/out/cs3_plate_qc/report.md (git-ignored; GDSC data are for
non-commercial research use and are not redistributed).
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from _report import Report, md_table

import tincture as tc
from tincture.io.gdsc import load_gdsc2_plates
from tincture.stats.zprime import Selector

NEG = Selector(tags=("NC-1",))
POS = Selector(well_types=("blank",))
EDGE = Selector(tags=("UN-USED", "NC-1", "NC-0"))


def run(max_plates: int | None = None) -> tuple[Path, bool]:
    t0 = time.time()
    pt = load_gdsc2_plates(max_plates=max_plates)
    load_s = time.time() - t0
    res = tc.plate_atlas(
        pt,
        neg=NEG,
        pos=POS,
        edge=EDGE,
        title="GDSC2 plate QC",
        footer=("GDSC2 raw data release 8.5 (27 Oct 2023), Wellcome Sanger Institute",),
    )
    q = res.meta["qc"]
    rep = Report("cs3_plate_qc", "CS3: Plate QC on GDSC2 raw monotherapy data")
    rep.text(
        f"Data: GDSC2 raw data release 8.5; {len(q):,} plates ({pt.n_wells}-well), "
        f"{len(pt.data):,} reported wells; loaded in {load_s:.0f} s.",
        "Z′ uses NC-1 (negative) vs B blanks (full-kill reference), as in the GDSC curve-fitting "
        "pipeline; standalone DMSO-tagged wells are not used as controls "
        "(see docs/data_status.md).",
    )
    rep.figure("plate_atlas", res)
    z = q["z_prime"].dropna()
    qs = z.quantile([0.05, 0.25, 0.5, 0.75, 0.95])
    rep.text(
        "## Z′ distribution",
        "| quantile | 5% | 25% | 50% | 75% | 95% |",
        "|---|---|---|---|---|---|",
        "| Z′ | " + " | ".join(f"{v:.2f}" for v in qs.to_numpy()) + " |",
        f"Plates with Z′ < 0.5: **{int(q['fails'].sum()):,} of {len(z):,} "
        f"({q['fails'].mean():.1%})**; "
        f"Z′ not computable: {int(q['z_prime'].isna().sum())}.",
    )
    worst = q.dropna(subset=["z_prime"]).nsmallest(15, "z_prime")[
        ["plate", "cell_line", "z_prime", "edge_ratio"]
    ]
    rep.text("Worst plates:", md_table(worst))
    q.to_csv(rep.dir / "plate_qc.csv", index=False)
    rep.text("Per-plate table: [plate_qc.csv](plate_qc.csv)", "## Ground-truth checks")
    ok = True
    ok &= rep.check(
        "Z′ distribution shown",
        len(z) > 0 and np.isfinite(qs).all(),
        f"median {qs.loc[0.5]:.2f}, IQR {qs.loc[0.25]:.2f}–{qs.loc[0.75]:.2f}",
    )
    er = q["edge_ratio"].dropna()
    ok &= rep.check(
        "Edge-vs-interior ratio reported",
        len(er) > 0,
        f"median {er.median():.2f} (IQR {er.quantile(0.25):.2f}–{er.quantile(0.75):.2f}) over "
        f"{len(er):,} plates; {res.meta['qc_object'].n_unseeded_outer:,} plates leave the outer "
        "ring unseeded",
    )
    ok &= rep.check(
        "Failing plates flagged",
        "fails" in q and int(q["fails"].sum()) >= 0,
        f"{int(q['fails'].sum()):,} plates flagged (terracotta in the figure; listed above)",
    )
    path = rep.write()
    return path, bool(ok)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-plates", type=int, default=None)
    a = ap.parse_args()
    p, ok = run(a.max_plates)
    print(p, "ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
