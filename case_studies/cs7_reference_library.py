"""CS7: a LINCS reference library retrieves known mechanisms (HANDOVER §6.2).

Pre-specified (fixed before the first run):
* Library: LINCS L1000 phase II (GSE70138) Level 5, landmark genes; one consensus signature
  per compound = mean of its 10 µM, 24 h signatures in the 7 core cell lines (≥ 3 lines);
  MoA from the Drug Repurposing Hub fields of the PRISM files (exact label match on the
  Broad ID core).
* Queries: every library compound whose MoA label has ≥ 3 compounds, searched against the
  library **without itself**.
* **Primary (ground truth):** share of queries whose top neighbour has the same MoA label,
  vs chance = the exact expected share under random ranking. Significance: 1,000 MoA-label
  permutations. Pass = ≥ 2× chance and permutation p < 0.01.
* Reported, not scored: the same for highly active compounds (mean TAS ≥ 0.3); the share of
  queries with a *confident* neighbour (q < 0.05) and how often it shares the MoA, under the
  handover's all-pairs null and the unrelated-pairs null.

    python case_studies/cs7_reference_library.py
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
from tincture.io.lincs import SOURCE, load_lincs
from tincture.model.reference import library_from_compounds
from tincture.model.signature import Signature
from tincture.stats.neighbours import neighbours

CORE = ["MCF7", "A375", "PC3", "HT29", "HA1E", "YAPC", "HELA"]
DOSE = 1e-5
PRISM = Path(__file__).resolve().parents[1] / "data/raw/prism"
TREAT = "prism-repurposing-20q2-secondary-screen-replicate-collapsed-treatment-info.csv"
FEATURED = ["bortezomib", "vorinostat", "trametinib"]
FOOT = (SOURCE, "MoA: Drug Repurposing Hub via PRISM 20Q2 (CC BY 4.0)")


def build_library():
    a = load_lincs(cell_lines=CORE, build=True)
    a = a[np.isclose(a.obs["dose"].astype(float).to_numpy(), DOSE)].copy()
    o = a.obs
    t = pd.read_csv(PRISM / TREAT)
    t["core"] = t["broad_id"].astype(str).str[:13]
    moa = t.dropna(subset=["moa"]).drop_duplicates("core").set_index("core")["moa"]
    X = np.asarray(a.X, float)
    sigs, ann, reps = {}, {}, {}
    for core, idx in o.groupby("pert_core").indices.items():
        if o["cell_id"].iloc[idx].nunique() < 3 or core not in moa.index:
            continue
        sigs[core], reps[core] = X[idx].mean(0), X[idx]
        r = o.iloc[idx[0]]
        ann[core] = {
            "moa": str(moa[core]),
            "name": str(r["pert_iname"]),
            "smiles": str(r["canonical_smiles"]),
            "tas": f"{o['tas'].iloc[idx].astype(float).mean():.3f}",
        }
    lib = library_from_compounds(
        sigs,
        a.var_names.to_numpy(),
        ann,
        version=f"GSE70138, 10 µM 24 h, consensus of {len(CORE)} core lines",
        source="LINCS L1000",
    )
    return lib, sigs, reps, a.var_names.to_numpy()


def evaluate(lib, sigs):
    moas = {c: lib.annotations[c]["moa"] for c in lib.names}
    counts = pd.Series(moas).value_counts()
    queries = [c for c in lib.names if counts[moas[c]] >= 3]
    rows = []
    for c in queries:
        loo = lib.without(c)
        r_all = neighbours(sigs[c], loo, null_kind="all", n_boot=0)
        r_unr = neighbours(sigs[c], loo, null_kind="unrelated", n_boot=0)
        top = r_all.names[0]  # ranking does not depend on the null
        rows.append(
            {
                "query": c,
                "moa": moas[c],
                "tas": float(lib.annotations[c]["tas"]),
                "top": top,
                "top_moa": moas[top],
                "hit": moas[top] == moas[c],
                "chance": (counts[moas[c]] - 1) / (len(lib.names) - 1),
                "conf_all": r_all.q[0] < 0.05,
                "conf_unr": r_unr.q[0] < 0.05,
            }
        )
    return pd.DataFrame(rows)


def permutation_p(lib, df, sigs, n_perm=1000, seed=0):
    """Shuffle MoA labels over the library; top neighbours stay fixed (they do not depend on
    labels), so only the label match is recomputed."""
    names = lib.names
    lab = np.array([lib.annotations[c]["moa"] for c in names])
    pos = {c: i for i, c in enumerate(names)}
    qi = np.array([pos[c] for c in df["query"]])
    ti = np.array([pos[c] for c in df["top"]])
    obs = df["hit"].mean()
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_perm):
        p = rng.permutation(lab)
        ge += np.mean(p[qi] == p[ti]) >= obs
    return (ge + 1) / (n_perm + 1)


def run() -> tuple[Path, bool]:
    rep = Report("cs7_reference_library", "CS7: LINCS reference library retrieves known MoA")
    lib, sigs, reps, feats = build_library()
    rep.text(
        f"Library: {len(lib.names)} LINCS compounds with a Repurposing Hub MoA (consensus of "
        f"10 µM, 24 h signatures over ≥ 3 of {len(CORE)} core lines; 978 landmark genes)."
    )
    names = {lib.annotations[c]["name"].lower(): c for c in lib.names}
    for nm in FEATURED:
        c = names.get(nm)
        if c is None:
            rep.text(f"_{nm} is not in the library._")
            continue
        q = Signature(
            features=feats,
            effect=sigs[c],
            compound_id=c,
            meta={"compound_name": nm, "smiles": lib.annotations[c]["smiles"]},
        )
        for kind in ("all", "unrelated"):
            res = tc.neighbourhood(
                q,
                lib,
                query_replicates=reps[c],
                null_kind=kind,
                title=f"Neighbourhood Card: {nm} ({kind}-pairs null)",
                footer=FOOT,
            )
            rep.text(f"## {nm}: null = {kind} pairs")
            rep.figure(f"card_{nm}_{kind}", res)
    df = evaluate(lib, sigs)
    rate, chance = df["hit"].mean(), df["chance"].mean()
    p = permutation_p(lib, df, sigs)
    hi = df[df["tas"] >= 0.3]
    summ = pd.DataFrame(
        [
            {
                "set": "all queries",
                "n": len(df),
                "top-1 same MoA": f"{rate:.0%}",
                "chance": f"{chance:.1%}",
                "× chance": f"{rate / chance:.1f}",
            },
            {
                "set": "TAS ≥ 0.3",
                "n": len(hi),
                "top-1 same MoA": f"{hi['hit'].mean():.0%}",
                "chance": f"{hi['chance'].mean():.1%}",
                "× chance": f"{hi['hit'].mean() / hi['chance'].mean():.1f}",
            },
        ]
    )
    conf = pd.DataFrame(
        [
            {
                "null": k,
                "queries with a confident neighbour": f"{df[c].mean():.0%}",
                "of which top shares MoA": f"{df[df[c]]['hit'].mean():.0%}" if df[c].any() else "—",
            }
            for k, c in (("all pairs (handover)", "conf_all"), ("unrelated pairs", "conf_unr"))
        ]
    )
    by_moa = (
        df.groupby("moa")
        .agg(n=("hit", "size"), hits=("hit", "sum"))
        .assign(rate=lambda x: x["hits"] / x["n"])
        .sort_values(["rate", "n"], ascending=False)
        .head(12)
        .reset_index()
    )
    by_moa["rate"] = by_moa["rate"].map("{:.0%}".format)
    rep.text(
        "## Held-out retrieval",
        md_table(summ),
        "### Confidence calls",
        md_table(conf),
        "### Best-retrieved MoA classes",
        md_table(by_moa),
        "## Ground-truth check",
    )
    ok = rep.check(
        "Held-out known-MoA compounds retrieve their own MoA well above chance",
        bool(rate >= 2 * chance and p < 0.01),
        f"top-1 neighbour shares the MoA for {int(df['hit'].sum())}/{len(df)} queries "
        f"({rate:.0%}) vs {chance:.1%} expected by chance ({rate / chance:.1f}×); permutation "
        f"p = {p:.3g} (1,000 label shuffles; criterion ≥ 2× and p < 0.01).",
    )
    rep.text("", f"**Overall: {'PASS' if ok else 'FAIL'}**")
    return rep.write(), ok


if __name__ == "__main__":
    path, ok = run()
    print(path, "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
