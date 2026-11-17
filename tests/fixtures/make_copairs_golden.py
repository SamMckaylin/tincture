"""Regenerate tests/fixtures/copairs_golden.json with copairs (NOT a tincture dependency).

Run in a throwaway environment, never in the project venv:

    uv venv /tmp/v_cp && uv pip install --python /tmp/v_cp/bin/python copairs==0.6.0
    /tmp/v_cp/bin/python -I tests/fixtures/make_copairs_golden.py tests/fixtures/copairs_golden.json

The input profiles are SYNTHETIC and stored inside the JSON, so the test does not depend on
random-number streams. Set-up = copairs' standard "phenotypic activity" recipe: positives are
replicates of the same compound, negatives are controls only (each control gets a unique
reference index).
"""

import json
import sys
from importlib.metadata import version

import numpy as np
import pandas as pd
from copairs.map import average_precision, mean_average_precision

rng = np.random.default_rng(7)
n_feat = 30
rows, labels = [], []
for _ in range(12):  # 12 controls
    rows.append(rng.normal(0, 1, n_feat))
    labels.append("DMSO")
effects = {"strong": 3.0, "medium": 1.2, "weak": 0.4, "none": 0.0, "noisy": 1.0}
for name, size in effects.items():
    direction = rng.normal(0, 1, n_feat)
    n_rep = 4 if name != "noisy" else 3
    for _ in range(n_rep):
        noise = rng.normal(0, 1 if name != "noisy" else 2.5, n_feat)
        rows.append(size * direction + noise)
        labels.append(f"synthetic_{name}")
X = np.vstack(rows)
meta = pd.DataFrame({"compound": labels})
meta["ref"] = np.where(meta["compound"] == "DMSO", np.arange(len(meta)), -1)

ap = average_precision(
    meta,
    X,
    pos_sameby=["compound", "ref"],
    pos_diffby=[],
    neg_sameby=[],
    neg_diffby=["compound", "ref"],
    progress_bar=False,
)
ap = ap[ap["compound"] != "DMSO"].dropna(subset=["average_precision"])
mp = mean_average_precision(
    ap,
    sameby=["compound"],
    null_size=20000,
    threshold=0.05,
    seed=0,
    progress_bar=False,
    null="copied",
)
out = {
    "copairs_version": version("copairs"),
    "X": X.round(10).tolist(),
    "labels": labels,
    "ap": ap[["compound", "average_precision", "n_pos_pairs", "n_total_pairs"]]
    .assign(profile=ap.index)
    .to_dict(orient="list"),
    "map": mp[["compound", "mean_average_precision", "p_value", "corrected_p_value"]].to_dict(
        orient="list"
    ),
    "null_size": 20000,
}
with open(sys.argv[1], "w") as fh:
    json.dump(out, fh)
print(mp)
