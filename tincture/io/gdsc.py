"""Loader for GDSC2 raw monotherapy screening data (release 8.5, 27 Oct 2023, Sanger).

Verified facts (docs/data_status.md):

* One row per well; 16,354 plates (``SCAN_ID``), all **1536-well**; ``POSITION`` is
  1-based **row-major on 32 × 48** (checked empirically from control layout).
* ``CONC`` is µM (library range 1e-5 to 2,000 µM); converted to molar.
* ``TAG`` families → tincture ``well_type``:

  ==================  ===========  =====================================================
  TAG                 well_type    meaning
  ==================  ===========  =====================================================
  ``L#-D#-S``         compound     library drug, dose step ``D#``
  ``R#-D#-S``         compound     reference drug, dose step
  ``T#-S``            compound     single-concentration treatment
  ``PC#-D#-S``/PC-#   pos_ctrl     positive-control drug (strong kill)
  ``NC-0``/``NC-1``   neg_ctrl     negative control (cells, no drug)
  ``DMSO``            see below    solvent annotation (drug id 4000, no concentration)
  ``B``               blank        media only, no cells (≈2% of NC-1 signal)
  ``UN-USED``         empty        unused for analysis (wells **do** contain cells)
  ==================  ===========  =====================================================

Unknown tags raise rather than being guessed.

**Duplicate annotations.** 12,946 of 16,354 plates list some positions twice: the same
physical well (identical intensity) annotated once as ``DMSO`` and once as its drug
(``L…``/``PC…``). Such wells are collapsed to the drug annotation (``solvent_paired=True``).
Standalone ``DMSO`` wells read ≈0.89 of NC-1, between negative controls (1.00) and drug
wells, so their role is ambiguous: they are typed ``empty`` (tag kept) and are never used
as controls. Z′ uses NC-1 (negative) vs B (blank) only. A well that still carries two
drug annotations is kept once with joined ids (``"a+b"``) and NaN dose.
"""

from __future__ import annotations

import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from tincture.datasets.cache import fetch
from tincture.datasets.registry import GDSC2_RAW
from tincture.errors import TinctureDataError
from tincture.model.measurements import round_dose
from tincture.model.plate import PlateTable, position_to_rowcol


def _drug_id_str(v: object) -> object:
    if isinstance(v, str):
        return v
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return np.nan
    return str(int(v))  # type: ignore[call-overload]


SOURCE = "GDSC2 raw data, release 8.5 (27 Oct 2023), Wellcome Sanger Institute"
N_WELLS = 1536

_TAG_RULES: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"^L\d+-D\d+-S$"), "compound"),
    (re.compile(r"^R\d+-D\d+-S$"), "compound"),
    (re.compile(r"^T\d+-S$"), "compound"),
    (re.compile(r"^PC\d+-D\d+-S$"), "pos_ctrl"),
    (re.compile(r"^PC-\d+$"), "pos_ctrl"),
    (re.compile(r"^NC-\d+$"), "neg_ctrl"),
    (re.compile(r"^DMSO$"), "empty"),
    (re.compile(r"^B$"), "blank"),
    (re.compile(r"^UN-USED$"), "empty"),
]

_USECOLS = [
    "BARCODE",
    "SCAN_ID",
    "SCAN_DATE",
    "CELL_LINE_NAME",
    "SANGER_MODEL_ID",
    "DRUGSET_ID",
    "SEEDING_DENSITY",
    "POSITION",
    "TAG",
    "DRUG_ID",
    "CONC",
    "INTENSITY",
]


def tag_to_well_type(tags: pd.Series) -> pd.Series:
    """Map GDSC ``TAG`` values to tincture well types; raise on unknown tags."""
    uniq = pd.Series(tags.unique())
    mapping: dict[str, str] = {}
    for t in uniq.astype(str):
        for pat, wt in _TAG_RULES:
            if pat.match(t):
                mapping[t] = wt
                break
    unknown = sorted(set(uniq.astype(str)) - set(mapping))
    if unknown:
        raise TinctureDataError(
            f"Unknown GDSC TAG value(s): {unknown[:10]}.",
            fix="This may be a different GDSC release; add the tag to tincture.io.gdsc._TAG_RULES "
            "after checking its meaning, and record it in docs/data_status.md.",
        )
    return tags.astype(str).map(mapping)


def collapse_duplicate_wells(raw: pd.DataFrame) -> pd.DataFrame:
    """Collapse wells listed twice (DMSO solvent row + drug row) into one row."""
    key = ["SCAN_ID", "POSITION"]
    n = raw.groupby(key)["TAG"].transform("size")
    is_dmso = raw["TAG"].astype(str) == "DMSO"
    paired_keys = raw.loc[is_dmso & (n > 1), key]
    out = raw[~(is_dmso & (n > 1))].copy()
    out["solvent_paired"] = pd.MultiIndex.from_frame(out[key]).isin(
        pd.MultiIndex.from_frame(paired_keys)
    )
    n2 = out.groupby(key)["TAG"].transform("size")
    if (n2 > 1).any():
        multi = out[n2 > 1]
        same = multi.groupby(key)["INTENSITY"].nunique()
        if (same > 1).any():
            raise TinctureDataError(
                f"{int((same > 1).sum())} GDSC wells are listed twice with different intensities.",
                fix="Inspect these wells; the release may be inconsistent.",
            )
        agg = multi.groupby(key, as_index=False).agg(
            {
                **{c: "first" for c in multi.columns if c not in key},
                "TAG": lambda s: "+".join(map(str, s)),
                "DRUG_ID": lambda s: "+".join(str(int(v)) for v in s.dropna()),
            }
        )
        agg["CONC"] = np.nan
        out = pd.concat([out[n2 == 1], agg], ignore_index=True)
    return out.reset_index(drop=True)


def read_raw(
    path: Path | None = None,
    plates: list[int] | None = None,
    cell_lines: list[str] | None = None,
    max_plates: int | None = None,
    chunksize: int = 2_000_000,
) -> pd.DataFrame:
    """Read (a filtered subset of) the raw CSV inside the release zip."""
    path = path or fetch(GDSC2_RAW)
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if n.endswith(".csv")]
        if len(names) != 1:
            raise TinctureDataError(f"Expected one CSV in {path.name}, found {names}.")
        parts: list[pd.DataFrame] = []
        seen: list[int] = []
        with z.open(names[0]) as fh:
            for ch in pd.read_csv(fh, usecols=_USECOLS, chunksize=chunksize, low_memory=False):
                if plates is not None:
                    ch = ch[ch["SCAN_ID"].isin(plates)]
                if cell_lines is not None:
                    ch = ch[ch["CELL_LINE_NAME"].isin(cell_lines)]
                if max_plates is not None:
                    for s in pd.unique(ch["SCAN_ID"]):
                        if s not in seen:
                            seen.append(int(s))
                    allowed = set(seen[:max_plates])
                    ch = ch[ch["SCAN_ID"].isin(allowed)]
                if len(ch):
                    parts.append(ch)
                if max_plates is not None and len(seen) > max_plates:
                    break
    if not parts:
        raise TinctureDataError(
            "No GDSC rows match the filters.", fix="Check plate ids / cell-line names."
        )
    return pd.concat(parts, ignore_index=True)


def load_gdsc2_plates(
    plates: list[int] | None = None,
    cell_lines: list[str] | None = None,
    max_plates: int | None = None,
    path: Path | None = None,
) -> PlateTable:
    """GDSC2 raw wells as a :class:`PlateTable` (value = raw luminescence intensity).

    Plates are ordered by scan date (``plate_order``), which is what control-drift
    analysis needs. ``tag`` keeps the original GDSC label (e.g. to choose NC-1 vs B for
    Z′). Doses are molar.
    """
    raw = collapse_duplicate_wells(
        read_raw(path, plates=plates, cell_lines=cell_lines, max_plates=max_plates)
    )
    counts = raw.groupby("SCAN_ID")["POSITION"].agg(["max", "size"])
    if (counts["max"] > N_WELLS).any():
        raise TinctureDataError("Found well positions > 1536; not a GDSC2 1536-well release.")
    row, col = position_to_rowcol(raw["POSITION"].to_numpy(), N_WELLS, order="row-major")
    first_tag = raw["TAG"].astype(str).str.split("+").str[0]
    well_type = tag_to_well_type(first_tag)
    dose = round_dose(raw["CONC"].to_numpy(dtype=float) / 1e6)
    dose = np.where(well_type.isin(["neg_ctrl", "blank", "empty"]), 0.0, dose)
    order = (
        raw.drop_duplicates("SCAN_ID")
        .assign(_d=pd.to_datetime(raw["SCAN_DATE"]))
        .sort_values(["_d", "SCAN_ID"])["SCAN_ID"]
        .reset_index(drop=True)
    )
    plate_order = pd.Series(np.arange(len(order)), index=order.to_numpy())
    data = pd.DataFrame(
        {
            "plate": raw["SCAN_ID"].astype(str),
            "row": row,
            "col": col,
            "well_type": well_type.to_numpy(),
            "value": raw["INTENSITY"].astype(float),
            "tag": raw["TAG"].astype(str),
            "compound_id": raw["DRUG_ID"].map(_drug_id_str),
            "dose": dose,
            "cell_line": raw["CELL_LINE_NAME"].astype(str),
            "model_id": raw["SANGER_MODEL_ID"].astype(str),
            "drugset_id": raw["DRUGSET_ID"],
            "seeding_density": raw["SEEDING_DENSITY"],
            "scan_date": raw["SCAN_DATE"],
            "plate_order": raw["SCAN_ID"].map(plate_order).to_numpy(),
            "solvent_paired": raw["solvent_paired"].to_numpy(),
        }
    )
    return PlateTable(
        data, n_wells=N_WELLS, source=SOURCE, value_label="luminescence (CellTiter-Glo)"
    )
