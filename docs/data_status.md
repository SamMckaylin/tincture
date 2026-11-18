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
| 2026-11-17 | **GDSC-derived image found in the repo.** `docs/img/restyle_before_after.png` (committed during M2) showed a GDSC plate map and per-plate QC, against the practice below. Cropped to its sci-Plex half in M4, and (user-approved) the branch history was rewritten so no commit contains the GDSC version. GitHub may keep the old unreferenced commits reachable by their hash for a while; removing those needs a request to GitHub Support. Gallery and README use sci-Plex or synthetic figures only. |
| 2026-11-17 | `static.omnipathdb.org`, `omnipathdb.org` added (gene sets for pathway scores). |
| 2026-11-18 | M5: `pmc.ncbi.nlm.nih.gov`, `www.nature.com`, `static-content.springer.com` added (user-approved environment change) to fetch the D7 and D8 supplementary tables. |
| 2026-11-18 | **D7 download blocked.** PMC (`pmc.ncbi.nlm.nih.gov`) answers file requests for PMC5814754 `mmc2.xlsx` with a proof-of-work browser challenge ("Preparing to download …"), not the file. Not worked around. PMC's old OA service (`www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi`) returns 404 and the PMC FTP no longer lists article packages. |
| 2026-11-18 | **D8 download blocked.** `www.nature.com/articles/nbt.1990` redirects scripts to the `idp.nature.com` login handshake; fetching the supplementary files directly was denied by this session's permission policy. Not worked around. |
| 2026-11-18 | `eutils.ncbi.nlm.nih.gov`, `pubmed.ncbi.nlm.nih.gov`, `cran.r-project.org` and `www.cell.com` are not allowed (proxy 403). |
| 2026-11-18 | M7: `gdsc-combinations.depmap.sanger.ac.uk` and `www.ebi.ac.uk` (Europe PMC) added (user-approved environment change). The GDSC Combinations portal is a Dash app whose download links are built by its JavaScript; querying its internal endpoints was refused by this session's permission policy, so only its public `Matrix_documentation.pdf` was read. All M7 data came from figshare and Europe PMC instead. |
| 2026-11-18 | Paper full texts for M7 (Jaaks 2022 PMC8891012, Bashi 2024 PMC11061612, Vis 2024 PMC11384948, McFarland 2020 PMC7453022) read through the Europe PMC REST API (`fullTextXML`), which serves open-access XML to scripts. |
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

## D7: NanoBRET live-cell kinase engagement (Vasta et al., Cell Chem Biol 2018)

**Status: verified (file supplied by the user on 2026-11-18, because the PMC download was blocked).** Placed at `data/raw/nanobret/vasta2018_mmc2.xlsx` and pinned by SHA-256 `b6348a34…77ac3` (26,112 bytes) in the registry as a manual file.
- Paper: Vasta et al., *Cell Chem Biol* 25(2):206–214.e11 (2018), PMC5814754.
- **The file holds occupancy, not IC50s.** Table S1 (one sheet, three header rows) has 180 kinase rows (LCK appears twice, once "Not Determined"; **178 unique kinases**, matching the paper).
  - Columns: kinase, protein ID, NanoLuc orientation, probe, probe concentration (header shows "mM"; a mangled µ is likely), Z′, probe EC50.
  - Then **% occupancy, two replicates each, for crizotinib and dasatinib**; values are integers, "No Binding" or "Not Determined".
- Conditions, verified from the article text: live **HEK293** cells expressing full-length kinase–NanoLuc fusions; single doses of **1 µM crizotinib** and **100 nM dasatinib**.
- The handover expects "NanoBRET IC50" for the cellular rung. The loader derives an **apparent IC50 = C·(1−f)/f** (1:1 binding, Hill 1) from occupancy f at dose C, labels it as model-based, and censors "No Binding".
- **Reproduces the paper:**
  - crizotinib: **13** kinases > 50% occupied in cells (paper: 13);
  - dasatinib: 23 (paper: 22; one kinase sits at the 50% boundary in the replicate means).

## D8: Davis et al. 2011 kinase Kd panel (Nat Biotechnol 29:1046)

**Status: verified (file supplied by the user on 2026-11-18).** Placed at `data/raw/davis/davis2011_table4.xls` and pinned by SHA-256 `e47821b2…498d7` (304,640 bytes). The file's metadata title is "Assay Matrix" and its author is P. Zarrinkar, a co-author.
- One sheet: `Accession Number`, `Entrez Gene Symbol`, `Kinase` (assay construct), then **72 inhibitors**. **442 assay rows** (the handover says ~440). **Kd in nM**, 0.016–9,900.
- Blank cells are **censored** (Kd > 10 µM, or not detected at 10 µM) and are loaded as censored, never as missing.
- Constructs: 13 ABL1 variants, 12 EGFR, 10 PIK3CA, 8 KIT, 7 FLT3, etc. These are mutants (flagged, excluded by default), phospho-states and single domains (wild type; the most potent per gene is kept), plus 3 non-human kinases (excluded). 379 unique genes.
- **Reproduces the paper:** with Kd ≤ dose as "> 50% occupied" (no ATP in the binding assay), crizotinib at 1 µM gives **69** of the NanoBRET kinases (paper: 69) and dasatinib at 100 nM gives 44 (paper: 44).
- PMC3230241, which a search returned for this paper, is a different study (Anastassiadis et al. 2011).

**Kinase name mapping D7 ↔ D8** (`tincture/io/kinase_names.py`, curated, each pair checked in both files): JNK1 = MAPK8, JNK3 = MAPK10, MLTK = ZAK, SNF1LK2 = SIK2, SIK3 = KIAA0999 (QSK). **Unmatched: PKN3** (not in the Davis panel).

## D9: LINCS L1000 (GEO GSE92742, GSE70138)

**Status: verified (listing, sizes); see below for the chosen subset.**
- Both series are on the GEO FTP host. GSE92742 (phase I) Level 5 `COMPZ.MODZ_n473647x12328.gctx.gz` is **20 GB** compressed (about 40 GB on disk unpacked), too large for this environment. **GSE70138** (phase II) Level 5 `COMPZ_n118050x12328_2017-03-06.gctx.gz` is **5.0 GB**, with `sig_info`, `pert_info`, `gene_info`, `cell_info` and `SHA512SUMS` alongside.
- Neither series ships MoA annotations. MoA comes from the Drug Repurposing Hub fields in the **PRISM** files already cached (CC BY 4.0), matched on the Broad ID core (`BRD-K12345678`).
- **Verified from the GSE70138 metadata (2026-11-18; SHA-512 checked against `GSE70138_SHA512SUMS.txt.gz`):**
  - `sig_info`: 118,050 signatures (columns `sig_id, pert_id, pert_iname, pert_type, cell_id, pert_idose, pert_itime, distil_id`):
    - 107,404 are `trt_cp` (compounds); the rest are vehicle, CRISPR, vector and untreated controls;
    - 1,796 compounds across **30 cell lines**; seven core lines carry about 10–12.7k signatures each (MCF7, A375, PC3, HT29, HA1E, YAPC, HELA), and A549 only 1,954;
    - doses are a 6-point series 0.04–10 µM (plus 20 µM for some); times are mostly 24 h (101,360), with some 3 h and 6 h.
  - `pert_info`: 2,170 perturbagens with `canonical_smiles` and `inchi_key`; the InChIKey is the cross-dataset ID.
  - `gene_info`: 12,328 genes, **978 landmark** (`pr_is_lm == 1`); only the landmarks are measured (the rest are inferred) and only they are used.
  - `sig_metrics` includes `tas` (transcriptional activity score) and `distil_cc_q75` (replicate correlation), used to gate weak signatures.
- **MoA coverage via PRISM:** 616 of 1,796 compounds match a PRISM `moa` by Broad ID core; 303 MoA labels, **63 with ≥ 3 compounds** (largest: EGFR inhibitor 18, HDAC inhibitor 14, topoisomerase inhibitor 12, MEK inhibitor 11, PI3K inhibitor 11).
- Level 5 matrix: downloading `GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz` (5.0 GB). Only the 978 landmark rows will be kept, in a local cache of about 0.5 GB.

## Kinase ATP Km values (for the Engagement Ladder ATP wedge)

**Status: three values verified against open-access full text (2026-11-18); no single consistent source available.** Vendor tables (Thermo SelectScreen, Reaction Biology, Carna, Eurofins) were unreachable (proxy 403 on every vendor host). Panel papers either assume a fixed Km (Anastassiadis 2011: "a Km,ATP of 10 µM for all kinases") or only state that ATP was run "at Km app" without numbers (a JCI 2019 SelectScreen supplement). Full texts read from the PMC open-access copy on `pmc-oa-opendata.s3.amazonaws.com`.

| Kinase | Km,ATP | Construct / conditions | Source |
|---|---|---|---|
| ABL1 | **9.78 ± 0.14 µM** | ABL "core" (N-cap–SH3–SH2–kinase, Sf9, His-tag); ADP Quest; 10 mM MgCl2; peptide at its Km | Grover et al., PLoS One 2015, PMC4519180: "For wild type ABL, we obtained Km values of 9.78 ± 0.14 μM … for ATP" |
| ALK | **0.134 ± 0.007 mM** (unphosphorylated), **0.159 ± 0.012 mM** (autophosphorylated) | kinase domain 1090–1416, His6, Sf9; 32P; 10 mM MgCl2, no Mn2+; peptide 1 mM; 25 °C | Bresler et al., Sci Transl Med 2011, PMC3319004, Table 1 ("Determined at 1mM peptide"). The authors note these are "higher than reported from Mn2+-containing assays". |
| MET | **60 ± 12 µM** (phosphorylated), **200 ± 38 µM** (unphosphorylated) | catalytic domain 956–1391, His6, Sf21; HTRF; 10 mM MgCl2; peptide 1 µM; 25 °C | Behshad et al., Curr Chem Genomics 2010, PMC2885599, Table 1 |

Caveats: three labs, three constructs and assay formats, so the values are not directly comparable; ALK and MET depend on phosphorylation state (MET 3.3-fold). Values for other kinases were not found or conflicted (KIT: 9 µM vs 370 µM). Not yet applied to CS6.

## D4: GDSC anchored combinations (Jaaks, Coker et al., Nature 2022)

**Status: verified, files cached and MD5-checked (figshare).**

- Sources (all named "Raw data supporting Jaaks and Coker et al., 2022"): fitted data, original screen (figshare 16843597, `Original screen_All tissues_fitted data.csv`, 198,893,607 bytes); fitted data, validation screen (16843600, 11,774,397 bytes); glossary (16895371, `Glossary.xlsx`). Raw plate data (19141916, 19141922) and day-1 plates exist but are not needed.
- **Licence conflict:** the figshare licence field says **CC BY 4.0**, but every item's description says "non-exclusive, non-transferable right to use data files for internal proprietary research and educational purposes … Excluded … resale … or provision of commercial services." Treated as the restrictive terms: outputs stay local, nothing GDSC-derived is committed.
- Handover facts confirmed: **296,707 measurements, 125 cell lines (breast 156,065 rows, colon 74,525, pancreas 66,117), 2,025 combinations**; **2 × 7 anchored design** (anchor at two concentrations; library 7-point, 1,000-fold, "two 2-fold dilution steps from the highest used concentration, all other dilution steps are 4-fold", Extended Data Fig. 1e).
- Columns (46, glossary checked): per anchor-concentration × library × cell line: `ANCHOR_VIABILITY`, `LIBRARY_EMAX`, `LIBRARY_XMID(_uM)`, `SYNERGY_OBS_EMAX`, `SYNERGY_EXP_EMAX` (**Bliss**), `SYNERGY_DELTA_EMAX` (= expected − observed), `SYNERGY_DELTA_XMID` (log2), `Synergy` (bool), plus fit RMSEs. Viabilities are 0–1. `ANCHOR_CONC` is text in the original screen because a few anchors are two-drug composites written as `0.125|0.003` and `0.5|0.03`; the loader parses it as a number (composites become NaN and are not used).
- **Synergy definition (paper, Methods):** synergistic if the combination IC50 < 2× the highest library concentration and ΔIC50 ≥ 3 (log2, 8-fold) or ΔEmax ≥ 0.2 over Bliss; replicates summarised as synergistic if ≥ half are.
- **CS4 ground truth (paper):** in colon, TOP1 inhibition was screened as **camptothecin** (irinotecan in pancreas), and AZD7762 (CHEK1/2) + camptothecin is synergistic, more in MSS than MSI lines (n = 31 vs 15) and in KRAS–TP53 double-mutant lines. In the data: AZD7762 (anchor 0.0625 and 0.25 µM) × camptothecin, and camptothecin (anchor 0.0025, 0.025 µM) × AZD7762, in 43–45 colon lines.
- **MSI status:** from Cell Model Passports `model_list_latest.csv.gz` (downloaded 2026-11-18, pinned by SHA-256): **30 MSS, 15 MSI** of the 45 colon lines; the paper reports 31 MSS (annotations may have been revised since).

## D5: GDSC full-matrix combinations (Bashi et al., Cancer Discov 2024)

**Status: verified, files cached and MD5-checked (figshare).**

- The paper's Data Availability names **figshare project 163378**; the GDSC Combinations portal is the browsing front-end. Licence on figshare: **CC BY 4.0** (no extra restriction in the descriptions). GDSC-derived outputs are still kept local, as for D3/D4.
- Files: `matrix_well_stats_GDSC_007_010_plus_Spiper1.zip` (figshare 25266697; 241,088,863 bytes; unzips to five per-screen CSVs, 2.1 GB, plus `Glossary_combo_well_stats.pdf`) = per-well library/anchor responses, combination response, HSA and Bliss expected values and excess; per-screen summaries `gdsc-007…010_matrix_results.csv`, `sandpiper-01_matrix_results.csv`.
- Per-well file semantics (from `Glossary_combo_well_stats.pdf`, inside the zip): one row per combination well; `lib1_conc`/`lib2_conc` in µM with dose tags D1 (highest) … D7; `inhibition` = 1 − viability of the combination well (**raw, not fitted**); `lib*_y` / `lib*_yhat` = measured / fitted **monotherapy** inhibition at that concentration; `HSA` = max(lib1_yhat, lib2_yhat); `Bliss_additivity` = lib1_yhat + lib2_yhat − lib1_yhat·lib2_yhat (checked numerically on the first rows); replicate matrices are separate `BARCODE` × `cmatrix`.
- Handover facts confirmed in the paper: **109 combinations, 755 cell lines, 7 × 7 matrix**, 1,000-fold range per drug. The portal's methods note (`Matrix_documentation.pdf`) defines **HSA excess = E1,2 − max(E1, E2)** and **Bliss excess = E1,2 − (E1 + E2 − E1E2)** on inhibition, states that there is **no fitted combination surface** (so no Loewe), and reports synergy over the full matrix and over the best 3 × 3 window.

## D6: GDSC pan-cancer anchored (Vis et al., Cell Rep Med 2024)

**Status: verified, file cached (SHA-256 computed on first download).**

- Data Availability: "supplied as supplemental data with this publication". Fetched as the open-access supplement zip from Europe PMC (53,745,336 bytes); Data S1 (`mmc7.zip`) contains `pancan_combi_fitted_data_jan21.tsv.gz` (same column scheme as D4) and `pancan_col_description.tsv`. Table S1A (`mmc2.xlsx`) is the Cell Model Passports annotation of the lines.
- Checked in the data (137,738 rows, 68 columns): **51 anchor–library pairs** (48 unordered drug pairs, as some are screened in both orientations), matching the paper's 51 combinations; **762 distinct `CELL_LINE_NAME`s** vs 757 in the paper (likely naming variants or lines dropped in QC; not resolved). Anchored design, same metric columns as D4 (`SYNERGY_DELTA_EMAX` etc.).

## D10: MIX-seq (McFarland, Paolella et al., Nat Commun 2020)

**Status: verified, files cached and MD5-checked (figshare 10298696 v3, CC BY 4.0).**

- **Time course exists for trametinib only:** experiment 5, `trametinib_tc_expt5.zip` (291,873,754 bytes; 10x matrix 1.2 GB unzipped), 24 cell lines, trametinib 0.1 µM vs DMSO, hash-tagged at **3, 6, 12, 24 and 48 h** (plus Untreated 48 h). Supplementary Table 2 (experiment metadata) lists 3, 6, 12, 24, 48 h; **Table 3 (drug metadata) lists 3, 6, 9, 12, 24, 48 h, but there is no 9 h hash tag in the data**, so the time points are those of Table 2.
- QC-passing singlets (`cell_quality == "normal"`) per hash: 915–2,256 cells per condition.
- Pseudobulk (CS8): UMI counts summed per cell line × hash tag over QC-passing singlets, streamed from the zip (19 s); cell counts per condition reproduce the QC table above (e.g. DMSO 3 h: 1,098). Cached under `data/processed/mixseq/`. Each time point has a matched DMSO group for 23–24 lines.
- Other drugs have 6 and 24 h only (bortezomib, idasanutlin, BRD3379: two time points) or 24 h only.

## D11–D12

Not started (backlog).

## Unmatched names

- Kinases: PKN3 (NanoBRET only).

## Unmatched compound names

None yet (cross-dataset mapping starts in M5/M6).
