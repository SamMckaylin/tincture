# Milestone 1: data model and loaders for D1–D3

## What was built

**Data model (`tincture.model`)**
- `validate_long`: validates the canonical long table (§4.1).
  - Required columns, numeric molar doses, a guard for µM mistakes (doses > 1 M), `well_type` vocabulary, ±inf, and an optional vehicle requirement.
  - Every error carries a fix suggestion.
  - Data without a dose column must state one explicitly via `add_constant_dose`.
  - `round_dose` removes unit-conversion float noise.
- `PerturbationSet` wraps AnnData; obs must have `compound_id`, `dose` (molar) and `replicate`.
  - `.subset(**equals, features=…)`, `.pseudobulk(by=…)` (sparse-aware, tracks `n_cells`), `.from_long()` / `.to_long()`.
  - Properties: `compounds`, `is_vehicle`, `doses()`, `var_names`.
  - `.signatures()` arrives in M2 with the stats layer.
- `Signature`: effect, SE and stored bootstrap draws.
  - `.ci()` uses bootstrap percentiles, falls back to normal from SE, and raises "no uncertainty" if it has neither.
  - `.align()` and `.magnitude`.
- `FourPL` (drc LL.4 convention, `.inverse()`), `DoseResponse` (raw points, optional fit, `quality_flags()`) and `DoseResponseSet` (points and fits tables; `.get()` returns a `DoseResponse`).
  - The quality flags cover no fit, non-convergence, R² < 0.8, asymptotes outside [−0.5, 1.5], and EC50 outside the tested range.
- `PlateTable`: 96/384/1536 geometry, row-major/column-major position mapping, edge mask, per-plate matrices (absent wells stay NaN), and validation of duplicates, bounds and well types.
- `well_type` gains **`blank`** (media only), which the handover's list lacked. GDSC uses blanks as its full-kill reference.

**Datasets (`tincture.datasets`)**
- `fetch()` downloads to a temp file, verifies MD5 or SHA-256 and size, then moves the file into place. It never uses a file that fails the checksum. On failure it says exactly where to put the file by hand.
- Cache root: `$TINCTURE_DATA`, else `data/raw` in a checkout, else `~/.cache/tincture`.
- `registry.py` pins 13 files: 5 sci-Plex (SHA-256 from first download), 7 PRISM (figshare MD5s) and 1 GDSC (SHA-256).

**Loaders (`tincture.io`)**
- `sciplex.build_pseudobulk()` streams the 3.1 GB triplet matrix once (about 5 min) into a per-well pseudobulk of human genes. `load_sciplex3()` returns a `PerturbationSet` filtered by cell line, time and compounds.
- `prism.load_prism_secondary()` returns a `DoseResponseSet` with raw replicate viability points (2^LFC) and the published drc fits. It prefers MTS010 for each compound and drops `_FAILED_STR` lines.
- `gdsc.load_gdsc2_plates()` returns a `PlateTable` (1536-well, row-major), ordered by scan date. It collapses DMSO-plus-drug double annotations and raises on unknown tags.

**Tests:** 158 unit and guard tests (synthetic fixtures only) plus 4 opt-in `-m realdata` round-trips. These pass on Python 3.12 and 3.10; lint, format and mypy are clean.

## Acceptance criteria

| Criterion | Status |
|---|---|
| Each loader round-trips into `PerturbationSet` / `DoseResponse` / `PlateTable` | ✅ `tests/realdata/test_loaders_real.py`: sci-Plex 4,992 wells × 58,347 genes (plus a `to_long`/`from_long` round-trip); PRISM dabrafenib → `DoseResponse` with fit; GDSC 20 plates → 1536-well `PlateTable` |
| Unmatched or surprising contents documented | ✅ `docs/data_status.md` (summary below) |

## What was verified, and what differed from the handover

Full detail is in `docs/data_status.md`.

**sci-Plex**
- 188 compounds × 3 lines confirmed.
- 799,317 cells, not ~650k; 762,795 are hash-assigned and already QC'd by the authors.
- 4 doses (10 nM–10 µM). A549 also has 72 h data, but for 47 compounds only.
- The gene table mixes human and mouse genes.
- The matrix is tab-separated. I first misread it, and the loader caught it.

**PRISM**
- Doses are µM and differ slightly per compound.
- There are 4 screens; MTS010 has 736 lines, not 489.
- **The README swaps upper/lower asymptotes, and its `convergence` column does not exist.**
- Many fits are pathological (EC50 up to 1e304, R² down to −5.8e5). Only about 12% of dabrafenib curves pass the quality flags; most failures are flat, non-responding lines.
- 407 of 1,421 compound × line pairs have points but no published fit.
- The SMILES field is repeated per batch.
- Control barcodes duplicate `row_name`s.

**GDSC**
- **cancerrxgene.org is retired** (HTTP 410). The data now live in the Sanger Cell Model Passports, whose host is not on the allowlist, so **the licence could not be re-read**. Release 8.5 files on `cog.sanger.ac.uk` still serve.
- All plates are 1536-well, row-major (verified empirically from control layout).
- Strong edge effect: outer wells sit at about 0.5–0.6 of the interior.
- Not every well is listed (708–1,392 rows per plate).
- 12,946 plates double-annotate wells as DMSO plus drug. Standalone DMSO wells are ambiguous and are never used as controls.

## What failed or changed along the way

- Network access needed a new **tincture** cloud environment with 8 dataset hosts (you approved it). The NCBI website shows a CAPTCHA to scripts, so metadata came from the GEO FTP SOFT file. The environment's setup script is untested until a fresh session starts on it.
- **GitHub push is blocked (403)**: Claude has no GitHub access to `SamMckaylin/tincture`. M0 and M1 are committed locally only.

## Open questions

1. **GitHub access:** connect GitHub or install the Claude app on the repo so commits can be pushed and CI can run.
2. **GDSC licence:** may I add `cellmodelpassports.sanger.ac.uk` to the environment to re-read the terms? Or can you confirm them? Until then, nothing derived from GDSC per-plate data is published.
3. **Compound IDs:** sci-Plex uses Selleck catalog numbers and PRISM uses Broad `BRD-K` ids. Cross-dataset mapping (M5/M6) needs InChIKeys, which means the `[chem]` extra (RDKit wheel about 35 MB, under the 50 MB line). OK to install it when needed?
4. Still open from §9 and M0: palette fix for novel/good (`docs/design/theme_colorblind.md`), `copairs` vs re-implementing mAP (needed in M2), stress axis, synergy default, package name.
