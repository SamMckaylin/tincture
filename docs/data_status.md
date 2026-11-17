# Data status

Ground rules: no data is committed (`data/` is git-ignored), no data is invented, and every "verify" fact in HANDOVER §6.1 is checked against its source before a loader is written. Anything that differs from the handover is recorded here.

## Environment / access log

| Date | Finding |
|---|---|
| 2026-11-16 | The Default environment blocked all dataset hosts. A personal **tincture** environment was created (approved by the user) allowing `ftp.ncbi.nlm.nih.gov`, `www.ncbi.nlm.nih.gov`, `figshare.com`, `api.figshare.com`, `ndownloader.figshare.com`, `depmap.org`, `www.cancerrxgene.org` and `cog.sanger.ac.uk`. |
| 2026-11-16 | `www.ncbi.nlm.nih.gov` serves a reCAPTCHA page to scripts, so GEO HTML pages can't be scraped. `ftp.ncbi.nlm.nih.gov` serves **files** fine but intermittently refuses **directory listings** (server-side 403). Metadata was therefore taken from `GSE139944_family.soft.gz` on the FTP host. |
| 2026-11-16 | figshare downloads redirect to a Cloudflare R2 storage host, and this works through the proxy. |
| 2026-11-16 | **cancerrxgene.org is retired**: every page returns HTTP 410 with a notice that GDSC data now live in the Sanger **Cell Model Passports** (cellmodelpassports.sanger.ac.uk/downloads) and DepMap Miner. Those hosts are *not* allowed, so the current GDSC licence page could not be read. The release files on `cog.sanger.ac.uk` are still served. |
| 2026-11-17 | `cellmodelpassports.sanger.ac.uk` added to the environment (user-approved) to read the GDSC licence. |
| 2026-11-17 | `static.omnipathdb.org`, `omnipathdb.org` added (gene sets for pathway scores). |
| 2026-11-16 | Nature, Cell and ACS author-guideline sites unreachable; journal presets are unverified (see `docs/design/journal_presets.md`). |

## D1: sci-Plex (Srivatsan et al., Science 2020), GEO GSE139944

**Status: verified, raw files cached.**

| Fact | Handover | Found | Note |
|---|---|---|---|
| Accession | GSE139944 | ✅ GSE139944; screen is sample **GSM4150378** "sciPlex3_A549_MCF7_K562_screen" | sciPlex1/2/4 are other experiments (barnyard; 4-drug A549 pilot; HDACi + acetyl-CoA) |
| Compounds | 188 | ✅ 188 + vehicle (189 `product_name` values) | |
| Cell lines | A549, MCF7, K562 | ✅ | |
| Cells | ~650k | **799,317 rows in pData**; 36,522 have no hash assignment (NaN cell_type, dose etc.), leaving 762,795 assigned | ~650k is the paper's post-QC number; our QC must be explicit and reported |
| Doses | "multiple" | **4 doses: 10, 100, 1,000, 10,000** (+ vehicle 0) for every compound | units are **nM** per the paper (10 nM–10 µM); file has no unit column; loader converts to molar |
| Time | — | **24 h** for all three lines; **A549 also 72 h** (82,110 cells, **47 of 188 compounds**) | 72 h is relevant for M7 (time) but covers a quarter of the library |
| Replicates | — | `rep1`, `rep2` (culture-plate replicates) | |
| Annotations | — | `product_name`, `target`, `pathway`, `pathway_level_1`, `pathway_level_2`; SMILES for all 188 compounds in `hashTable_metadata` (CXSMILES with a `|…|` suffix that must be stripped) | 27 compounds have `target` containing "HDAC" (14 named *-inostat* etc.) |
| Matrix layout | verify | `UMI.count.matrix.gz`: headerless, **tab-separated** triplets `gene_index cell_index count`, 1-based, row order of `gene.annotations` and `cell.annotations` | |
| Genes | — | `gene.annotations` has 110,983 rows: **58,347 human (ENSG) + 52,636 mouse (ENSMUSG)**, a combined reference | loader keeps human genes only |
| Processed | — | per-well pseudobulk (4,992 wells × 58,347 human genes; 4,608 wells at 24 h; median 129 cells/well) built in about 5 min by `tincture.io.sciplex.build_pseudobulk` → `data/processed/sciplex3_pseudobulk_wells.h5ad` | the treated condition has 2 wells (rep1/rep2); vehicle has 32 wells per line at 24 h |
| Other files | — | `.RDS` (R) objects of 0.55–2.4 GB | **not used**: would need R |

Files cached in `data/raw/sciplex/` (SHA-256 from first download on 2026-11-16):

| File | Bytes | SHA-256 |
|---|---|---|
| GSM4150378_sciPlex3_A549_MCF7_K562_screen_UMI.count.matrix.gz | 3,096,224,606 | 7d632716aa6ed0fc1780996003d0abc440ff78340609bceaaa4c8ade9345d00a |
| GSM4150378_sciPlex3_A549_MCF7_K562_screen_cell.annotations.txt.gz | 2,770,264 | 5db11591747a49f07ceb90a070971d60816591911f1022171adba7d8d032f379 |
| GSM4150378_sciPlex3_A549_MCF7_K562_screen_gene.annotations.txt.gz | 806,168 | fbe43028cfb75dc5ebf383cbbbb100da6b51e7be8e2e74f6354f7c59b9fda7c2 |
| GSM4150378_sciPlex3_pData.txt.gz | 84,673,357 | 860ba5c21846e2cfe97c39c805895982ef9db4b108cdf8b841a4a91be65024c2 |
| GSM4150378_sciPlex3_A549_MCF7_K562_hashTable_metadata.txt.gz | 64,715 | 98a65191d343d149c5742cf6e7b734d403e853d165da800fa2daea41b978a2ff |
| GSM4150378_sciPlex3_hashSampleSheet.txt.gz | 6,441 | ac66d0dfdcc3c170af9beafe2843d0250a10f753f7bf99968aafe7a2e7279084 |

## D2: PRISM Repurposing 20Q2 secondary (Corsello et al., Nat Cancer 2020), figshare 10.6084/m9.figshare.20564034.v1

**Status: verified, files cached and MD5-checked.** Licence: **CC BY 4.0**.

| Fact | Handover | Found | Note |
|---|---|---|---|
| Compounds × lines | 1,448 × 489 | README says the same; data: 1,547 `broad_id` batches / 1,488 names across 4 screens | `broad_id` is drug-*batch*; map to compounds by name/InChIKey |
| Screens | prefer MTS010 | **HTS002** (1,396 batches, 482 lines), **MTS010** (192, 736 lines), **MTS006** (73), **MTS005** (2) | README: MTS010 = 192 oncology compounds (152 redos) × 770 lines; figshare page says 147 redos. MTS010 has a much larger panel than 489 |
| Dilution | 8-point, 4-fold from 10 µM | ✅ 8 doses per compound in HTS002; doses in **µM** and **vary slightly per compound** (actual plated concentration; max 12.66 µM) | do not assume a fixed grid |
| Replicates | — | mostly 3 per compound × dose (some 4–6; vehicle more) | |
| Curve params | ec50, ic50, auc, r2 | ✅ present: `upper_limit`, `lower_limit`, `slope`, `r2`, `auc`, `ec50`, `ic50` (753,778 curves) | **no `convergence` column**, although the README lists one |
| Asymptotes | — | `upper_limit` median 1.0 = **no-drug** asymptote; `lower_limit` = high-dose asymptote (drc LL.4 convention) | **the README's descriptions of upper/lower are swapped**; trust the data |
| Fit quality | — | many pathological fits: EC50 up to 1e304, R² down to −5.8e5, `lower_limit` down to −4e7; `ic50` NaN for 46% (curve never crosses 50%) | loader flags fits (R² < 0.8, asymptotes outside [−0.5, 1.5], EC50 outside tested range); Window guard must use the flags |
| Raw points | — | replicate-level `logfold_change` matrix (rows = `PR500_ACH-…` cell-line rows, cols = `column_name` in treatment info); viability = 2^LFC | |

Cached files (MD5 from the figshare API, verified): readme, cell-line-info, pooling-info, dose-response-curve-parameters (290 MB), replicate-treatment-info, replicate-collapsed-treatment-info, logfold-change (308 MB). See `tincture/datasets/registry.py`.

## D3: GDSC monotherapy raw data

**Status: verified format, file cached. Licence: could not be re-read (see access log).**

- **Release chosen: GDSC2, release 8.5 (27 Oct 2023)**, `https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_public_raw_data_27Oct23.csv.zip` (122,676,850 bytes; SHA-256 `7bc484eb0e859d58bfca30f72121cd353b26d35b28d18549464d496f3fc1bfc6`). Unzips to a 2.11 GB CSV.
- GDSC2 (not GDSC1) was chosen because it is the newer screen (CellTiter-Glo, 7-point) and is the one used in the GDSC combination papers. GDSC1 release 8.5 (66 MB zip) is also available.
- Columns: `RESEARCH_PROJECT, BARCODE, SCAN_ID, DATE_CREATED, SCAN_DATE, CELL_ID, MASTER_CELL_ID, COSMIC_ID, CELL_LINE_NAME, SANGER_MODEL_ID, SEEDING_DENSITY, DRUGSET_ID, ASSAY, DURATION, POSITION, TAG, DRUG_ID, CONC, INTENSITY`.
- 15,069,449 well rows; **16,354 plates, all 1536-well** (`POSITION` 1–1536); 299 drugs; 969 cell lines; `ASSAY` = Glo for all; `DURATION` = 4 for all (unit not stated in file).
- `TAG` values: `UN-USED`, `DMSO`, `NC-1`, `NC-0`, `B`, `PC1-D1..5-S`, `PC2-D1..5-S`, `PC-1`, `L1/L2-D1..D7-S` (library drug, 7 concentrations), and others.
- **Control semantics, checked empirically** on the first 600 plates (intensity relative to the per-plate NC-1 median): `B` = 0.02 (**blank**, no cells); `NC-0` ≈ `NC-1` = 1.00 (**negative controls**); `PC-1` = 0.14 and `PC-2` = 0.19 (**positive controls**, strong kill); `UN-USED` = 0.94 (**wells contain cells but no drug**, so not "empty"); library tags are `L<n>-D<k>-S` (drug slot n, dose k). Plate Atlas Z′ defaults to NC-1 (neg) vs B (full-kill reference) and states the pair used.
- **Well orientation, checked empirically:** `POSITION` is **row-major on a 32 × 48 grid** (A1…A48, B1…). Under this mapping the blanks form a contiguous block (column 1, rows 3–30) and all four edges are depleted (outer row/column medians 0.51–0.59 of interior, a strong **edge effect**). Column-major gives a scattered, implausible layout.
- **Not every well is listed:** 708–1,392 rows per plate (median 922). The PlateTable therefore has NaN cells for absent wells, and the Plate Atlas must show them as "not reported", not as zero.
- **Duplicate annotations:** 12,946 of 16,354 plates list some positions **twice**, with the same intensity: once as `DMSO` (drug id 4000, no concentration) and once as the drug (`L…`/`PC…`). So `DMSO` is a **solvent annotation**, not a vehicle-control label, in those wells. The loader collapses them to the drug annotation (`solvent_paired=True`). About 2.1 M **standalone** `DMSO` rows read median 0.89 of NC-1, between negative controls (1.00) and library wells (0.95). Their role is ambiguous, so they are typed `empty` and **never used as controls**; Z′ uses NC-1 vs B only. `T#-S` wells also appear paired (5,397 rows in a 3 M-row sample).
- **Loader cost:** filtering by plate or cell line scans the whole 2.1 GB CSV (about 15 s); the first N plates load in about 5 s.
- **Licence (checked 2026-11-17 via Cell Model Passports, where GDSC now lives).** The site states: "All data in the Cell Model Passports is publicly available and individual users are free to use the API for **non-commercial use**" and "**Commercial use** of the Cell Model Passports API and/or data **is not permitted without prior consent**" (contact depmap@sanger.ac.uk). The full "Terms & Conditions" page is on depmap.sanger.ac.uk (not allowed), so **redistribution terms remain unconfirmed**. Practice: raw and derived per-plate GDSC data stay local and are never committed or published; case-study figures are for internal research use. The Passports downloads page lists the same release (`GDSC2_public_raw_data_27Oct23`) under `cmp.cog.sanger.ac.uk`, plus `screened_compounds_rel_8.5_inchi.csv` (InChI per GDSC compound, useful for cross-dataset ID mapping in M5/M6).

## Gene-set libraries (pathway scores, M2)

Fetched through decoupler 2.3.1 from OmniPath (`static.omnipathdb.org`, `omnipathdb.org`, allowed 2026-11-17). The content can change between fetches, so the first fetch is cached at `data/raw/genesets/` with a JSON sidecar, and later loads verify it. The version string is printed on figures.

| Library | Sets | Rows | Fetched | SHA-256 |
|---|---|---|---|---|
| MSigDB Hallmark (human) | 50 | 7,318 | 2026-11-17 | `41042764841a61b2208d79cbf98dcc18460a32f79712851dfa2952d25a9dccb3` |
| PROGENy top-500 (human, signed weights) | 14 | 6,463 | 2026-11-17 | `bb3461bdd10b4245ba19ad7db5333a11c67abce26b5453044f30e268e47a8d50` |

Gene identifiers are symbols. sci-Plex features are Ensembl ids and are mapped through `gene_symbol` (first occurrence wins for duplicated symbols). Our ULM matches `decoupler.mt.ulm` to within 6e-13 on real signatures.

## D4–D12

Not started (M5–M7).

## Unmatched compound names

None yet (cross-dataset mapping starts in M5/M6).
