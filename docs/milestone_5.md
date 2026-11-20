# Milestone 5: Window, Selectivity Spectrum, Neighbourhood Card; loaders for D7, D8, D9

## What was built

**Statistics**
- `stats.hill`: 4-parameter log-logistic fits.
  - Fitting is bounded, soft-L1 least squares on log dose, with six starts for the point estimate and warm starts for the bootstrap (8× faster).
  - Two effective-concentration definitions: `ec_relative` (e.g. EC90 = 90% of the curve's own span) and `ec_drop` (e.g. EC10 of a toxicity readout = 10% loss of viability; `inf` when never reached).
  - **Median-of-units curves:** a curve that summarises many cell lines is fitted to per-dose medians, and the bootstrap resamples whole lines. Raw wells from different lines are not replicates; pooling them gave R² 0.03.
- `stats.neighbours`: cosine neighbours of a query in an annotated library.
  - Empirical null from library pairs, giving a percentile and p per neighbour, then BH q across the library.
  - Bootstrap intervals over the query's replicate signatures.
  - Two null options, see *Decision needed*: `"all"` pairs (handover, default) and `"unrelated"` pairs.

**Data**
- **D9 LINCS L1000 phase II (GSE70138): verified and loaded.**
  - The 5 GB Level 5 file is checked against GEO's SHA-512 sums, reduced to the 978 measured landmark genes (a 0.5 GB cache, built in 17 s once unpacked), and the unpacked file is deleted.
  - MoA comes from the Drug Repurposing Hub fields in the PRISM files (CC BY), matched by Broad ID. Facts are in `docs/data_status.md`.
- **D7 NanoBRET (Vasta 2018) and D8 Davis 2011: blocked, not worked around.**
  - PMC answers the D7 spreadsheet request with a proof-of-work bot challenge.
  - nature.com redirects scripts to a login handshake, and fetching the D8 file directly was denied by this session's permission policy.
  - Paper facts were verified from the article pages. Two corrections to the handover: Davis has **442** assays (the handover says ~440), and blank cells are *censored* (> 10 µM), not missing. A wrong lead was also recorded: PMC3230241 is a different paper, not Davis.
  - Exact manual steps are in `docs/data_status.md`. **These loaders are not written**, because their table layouts cannot be verified without the files (ground rule).
- `model.potency.PotencyTable`: compound × target pX with assay type, context and censoring.

**Plots**
- `tc.window` (§5.7):
  - shared log-concentration axis with raw points;
  - EC90(on) and EC10(tox) with bootstrap intervals;
  - green window band; Cmax as a range band, drawn only when supplied (PRISM has none, and the caption says so);
  - median-of-lines toxicity curves shown as medians with interquartile ranges.
  - Guards (each tested): a rejected on-target fit (no convergence, R² < 0.8 or effect < 50%) means no band, with every reason stated; a flat toxicity curve is censored ("> top dose", open-ended band), not failed; a toxicity curve that declines but fits badly means no band; the most sensitive toxicity curve limits the window; assay context is printed on each curve.
- `tc.neighbourhood` (§5.9):
  - RDKit structure next to the top-k neighbours;
  - similarity bars with 95% intervals, coloured by MoA only when confident (q < 0.05);
  - null percentile lines; the query is always left out of its own library.
  - **Guard (tested):** "no confident neighbour" replaces the ranked list, and the caption gives the best neighbour's percentile and why it fails.
- `tc.selectivity` (§5.8):
  - potency strip, strongest first; intended target in terracotta;
  - off-targets above the free-Cmax band in peach, or, without exposure data, within a stated 10-fold rule;
  - censored "not detected" values hollow and counted; panel coverage ("n of ~518 kinases") in the footer.
  - **Guard (tested):** mixed assay types are refused unless `allow_mixed=True`, which gives distinct markers and a caption warning.
  - Tested on synthetic panels only (D7/D8 blocked).
- Spec and renderer: `ImageMark` for structures, per-reference library annotations, SHA-512 checksums.
- Gallery extended with Window (CS2), Neighbourhood Card (CS7: bortezomib with the unrelated-pairs null, labelled as such) and a synthetic Selectivity Spectrum.

**Tests:** 300 passing, plus 8 opt-in real-data tests. New baselines for Window, Neighbourhood and Selectivity, all colour-blind-checked.

## Acceptance

| Criterion | Status |
|---|---|
| CS2 report | ✅ `case_studies/out/cs2_prism_windows/report.md`. **Pre-specified:** dabrafenib; on-target = each line's curve; toxicity = median of 676 non-skin lines. **12/44 skin lines vs 9/676 non-skin lines have a window** (one-sided Fisher p = 1.1e-10). A375: 275-fold window (EC90 8.4 nM vs toxicity EC10 2.3 µM). Vemurafenib agrees (16% vs 0%). Thyroid ranks third; BRAF V600E is common there. |
| CS7 report | ✅ `case_studies/out/cs7_reference_library/report.md`. **Pre-specified:** 603-compound LINCS consensus library; 325 held-out queries in MoA classes with ≥ 3 members. **The top-1 neighbour shares the exact MoA label for 96/325 (30%) vs 1.0% by chance (30×, permutation p = 0.001).** 41% for highly active compounds; HDAC 12/13, MEK 9/11, proteasome 5/5. |
| Neighbourhood "no confident neighbour" path tested | ✅ `tests/test_neighbourhood.py` |
| Loaders for D7, D8, D9 | D9 ✅. **D7 and D8 ❌ blocked** (downloads); see above. |

## What failed or changed along the way

1. **Pooled toxicity curves.** Fitting raw wells from 676 lines gave R² 0.03, because between-line spread dominates. This led to the median-of-units mode with a line bootstrap.
2. **Rejection reasons were incomplete.** The fit checks stopped at the first failure, so "R² 0.79" could hide "effect < 50%". Every failing check is now listed.
3. **The atlas preview hid most lines.** It showed only the 25 with a usable fit. Now every lineage counts window / too little effect (a real negative) / poor fit (not judged).
   - This exposed something worth knowing: **the R² ≥ 0.8 guard rejects most PRISM curves** (25/44 skin lines are "poor fit"), because flat or noisy curves can't reach R² 0.8. The CS2 counts are therefore conservative.
4. **The handover's all-pairs neighbour null is conservative.** The pairs it should detect, same-MoA pairs, crowd its upper tail: in LINCS, 24 of the 39 pairs with cosine ≥ 0.9 share an MoA. Vorinostat's top neighbour, belinostat (cosine 0.92, above 99.9% of pairs), gets q = 0.061. See *Decision needed*.
5. **A misleading caption.** The "no confident neighbour" caption said the best similarity was "within what random pairs reach" even when it was above the 99.9th percentile. It now gives the percentile and says the failure comes from the multiple-comparison correction.
6. **An unverified duration label.** I almost printed "PRISM 5 d" on the curves; the PRISM readme does not state a duration, so the label is "PRISM secondary viability".
7. **Disk:** the LINCS Level 5 download is 5.4 GB and is kept, git-ignored, so the cache can be rebuilt without re-downloading. Free disk after M5 is about 17 GB.

## Decision (resolved 2026-11-18): keep the handover's all-pairs null as the default

**Neighbourhood null:** keep the handover's "all library pairs" (current default) or switch to "unrelated pairs" (pairs with different MoA labels)? CS7 measured the trade-off:

| null | queries with a confident neighbour | of which the top shares the exact MoA |
|---|---|---|
| all pairs (handover) | 6% | 89% |
| unrelated pairs | 14% | 60% |

"Unrelated" makes more calls but they are less precise. Under the default, even textbook cases get no confident neighbour: vorinostat's top neighbour belinostat (cosine 0.92) has q = 0.061, and bortezomib's top neighbour oprozomib (cosine 0.92) has q = 0.056. Both are in the top 0.1% of library pairs, but the correction for 603 comparisons needs p < 8e-5, i.e. fewer than about 14 of 181k null pairs above them, and same-class pairs fill that tail. Exact-label matching understates its accuracy, since related mechanisms often carry different labels (e.g. PI3K vs mTOR inhibitor), but not by enough to close the gap. I'd keep the handover default and offer `null_kind="unrelated"` as an option, which is the current state.

## Open questions

1. ~~**D7 / D8 files**~~ (supplied by the user on 2026-11-18; loaders written and CS6 run in M6).
   - D7: `mmc2.xlsx` from https://pmc.ncbi.nlm.nih.gov/articles/PMC5814754/, placed at `data/raw/nanobret/vasta2018_mmc2.xlsx`.
   - D8: Supplementary Table 4 from https://www.nature.com/articles/nbt.1990, placed at `data/raw/davis/davis2011_table4.xls`.
2. **The Window R² guard:** keep R² ≥ 0.8 as the handover specifies, or add a second test that recognises flat curves as "no effect" rather than "poor fit" (an F-test of the 4PL against a constant)? This only changes how rejected lines are labelled, not whether a window is drawn.
3. **Still open from §9:** package name (a PyPI check), and Bliss vs HSA as the default synergy model (M7).
