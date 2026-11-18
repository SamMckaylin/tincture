"""Loader for Davis et al. 2011 kinase binding constants (Nat Biotechnol 29:1046; D8).

Facts this loader relies on were verified against the file (docs/data_status.md, D8):

* Supplementary Table 4, one sheet: columns ``Accession Number``, ``Entrez Gene Symbol``,
  ``Kinase`` (assay construct), then **72 inhibitors**; **442 assay rows**; values are
  **Kd in nM** (0.016 to 9,900).
* **Blank cells are censored, not missing:** the combination was tested, but binding was
  weak (Kd > 10 µM) or not detected in the 10 µM primary screen. They become
  ``censored = ">"`` at 10 µM (pX 5).
* Constructs: mutants (e.g. ``ABL1(T315I)``, ``EGFR(L858R)``), phospho-states
  (``ABL1-phosphorylated``), single domains (``JAK1(JH1domain-catalytic)``) and three
  non-human kinases (P. falciparum, M. tuberculosis). 379 unique human genes.

``load_davis2011(per_gene=True)`` keeps wild-type human constructs and takes the most potent
Kd per gene (the construct used is recorded); ``per_gene=False`` returns every assay row.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from tincture.datasets.cache import fetch
from tincture.datasets.registry import DAVIS2011_TABLE4
from tincture.model.potency import PotencyTable

SOURCE = "Davis et al., Nat Biotechnol 2011, Supplementary Table 4 (Kd, nM)"
CENSOR_NM = 10_000.0
_MUT = re.compile(r"\(([^)]*)\)")
_MUT_TOKEN = re.compile(r"^[A-Z]\d+[A-Z]$|ITD|del|ins|^[A-Z]\d+-[A-Z]\d+")


def classify(construct: str) -> str:
    """``"mutant"``, ``"nonhuman"`` or ``"wildtype"`` for a Davis assay construct name."""
    m = _MUT.search(construct)
    inner = m.group(1) if m else ""
    if "falciparum" in inner or "tuberculosis" in inner:
        return "nonhuman"
    tokens = [t.strip() for t in re.split(r"[,\s]+", inner) if t.strip()]
    if any(_MUT_TOKEN.search(t) for t in tokens):
        return "mutant"
    return "wildtype"


def read_table4(root: Path | None = None) -> pd.DataFrame:
    path = fetch(DAVIS2011_TABLE4, root)
    t = pd.read_excel(path, header=0)
    t = t.rename(
        columns={
            "Entrez Gene Symbol": "gene",
            "Kinase": "construct",
            "Accession Number": "accession",
        }
    )
    t["construct_kind"] = t["construct"].astype(str).map(classify)
    return t


def load_davis2011(
    compounds: list[str] | None = None, per_gene: bool = True, root: Path | None = None
) -> PotencyTable:
    """Davis 2011 Kd values as a :class:`PotencyTable` (assay type ``Kd``)."""
    t = read_table4(root)
    drugs = [c for c in t.columns[3:] if c != "construct_kind"]
    if compounds is not None:
        low = {c.lower(): c for c in drugs}
        missing = [c for c in compounds if c.lower() not in low]
        if missing:
            from tincture.errors import TinctureDataError

            raise TinctureDataError(f"Not in Davis 2011: {missing}.", fix=f"Use one of {drugs}.")
        drugs = [low[c.lower()] for c in compounds]
    long = t.melt(
        id_vars=["accession", "gene", "construct", "construct_kind"],
        value_vars=drugs,
        var_name="compound_id",
        value_name="kd_nm",
    )
    long["censored"] = np.where(long["kd_nm"].isna(), ">", "")
    long["kd_nm"] = long["kd_nm"].fillna(CENSOR_NM).astype(float)
    long["px"] = 9.0 - np.log10(long["kd_nm"])
    long["assay_type"] = "Kd"
    long["context"] = "biochemical binding (KINOMEscan, no ATP)"
    if per_gene:
        long = long[long["construct_kind"] == "wildtype"]
        # most potent measured value per gene; censored only if every construct is censored
        long = long.sort_values(
            ["compound_id", "gene", "censored", "px"], ascending=[True, True, True, False]
        )
        long = long.drop_duplicates(["compound_id", "gene"])
    out = long.rename(columns={"gene": "target"})[
        [
            "compound_id",
            "target",
            "px",
            "assay_type",
            "censored",
            "construct",
            "construct_kind",
            "kd_nm",
            "context",
        ]
    ].reset_index(drop=True)
    return PotencyTable(out, source=SOURCE)
