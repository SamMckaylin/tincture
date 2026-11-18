# Changelog

## Unreleased

### M7: Combination Surface, Combination Lift, Onset Timeline (2026-11-18)
- `model.combo.ComboMatrix` (design inferred: full / anchored / partial) and `stats.synergy`: Bliss independence (default), ZIP (Yadav et al. 2015) and HSA; ceiling mask (expected inhibition > 90%); per-replicate t-tests (≥ 2 replicates) or a parametric bootstrap with a stated noise SD (one replicate), BH across cells, significant region = q < 0.05 and excess ≥ 10%.
- `tc.combination_surface`: observed contours inside the grid only, significant region outlined, ceiling cells hatched, free-Cmax box, single-agent margins, excess heatmap; refuses anchored designs (points to the Lift). `tc.combination_lift`: best single agent vs combination, benefit band, groups with Fisher / Mann–Whitney, failed fits hollow, facets per anchor concentration.
- `stats.onset` + `tc.onset_timeline`: interval-censored onset, peak, resolution; Gantt rows ordered by onset; horizon style for ≥ 6 time points only.
- Loaders: `io.gdsc_combo` (Jaaks 2022, Bashi 2024 per-well 7×7, Vis 2024, Cell Model Passports annotations) and `io.mixseq` (streaming pseudobulk of the trametinib time course).
- CS4 PASS (AZD7762 + camptothecin benefit in colon, MSS 59% vs MSI 20% at the low anchor); CS5 (selumetinib + venetoclax vs olaparib + AZD5153) see docs/milestone_7.md; CS8 PASS (MAPK output down by 3 h; p53 from 12–24 h; Hallmark apoptosis never active).
- Kinase ATP Km values (ABL1, ALK, MET) found and verified in open-access papers; recorded, not yet applied.

### M6: Engagement Ladder (2026-11-18)
- `io.davis` (Davis 2011 Kd, 442 assays × 72 inhibitors; censored blanks; construct classes), `io.nanobret` (Vasta 2018 single-dose occupancy → apparent IC50, model-based), `io.kinase_names` (curated synonyms); both files pinned as manual downloads. Reproduces the paper: crizotinib 13 engaged in cells vs 69 bound biochemically.
- `tc.engagement_ladder`: rungs on a shared log axis, intervals, small multiples per target; off-target, ATP-shift (Cheng–Prusoff wedge, needs ATP Km) and exposure flags; dashed across cell systems; missing rungs skipped; `<` and `>` censoring used only when conclusive.
- CS6: crizotinib selectivity reproduced; dasatinib/ABL1 off-target flag fires (K562 vs HEK293 caveat stated); ATP-shift flag not evaluable without ATP Km.
- Selectivity Spectrum on real data (crizotinib, Davis and NanoBRET); stacked labels with leader lines. `PolygonMark`.

### M5: Window, Neighbourhood Card, Selectivity Spectrum (2026-11-18)
- `stats.hill`: bounded 4PL fits (soft-L1), EC90/EC10 definitions, warm-started bootstrap, median-of-units curves with unit bootstrap.
- `tc.window`: fit-quality guards (convergence, R² ≥ 0.8, effect ≥ 50%), censored toxicity, Cmax as a range band, assay context; classic overlaid curves. CS2 (PRISM, dabrafenib): skin lines have windows (12/44 vs 9/676, p = 1e-10).
- `io.lincs`: LINCS L1000 phase II Level 5 (GSE70138, SHA-512 pinned) reduced to 978 landmark genes; MoA from the Repurposing Hub via PRISM.
- `stats.neighbours` + `tc.neighbourhood`: cosine neighbours with a library-pair null (`"all"` per handover or `"unrelated"`), BH q, replicate bootstrap CIs, RDKit structure, "no confident neighbour" guard. CS7: held-out top-1 neighbour shares the MoA 30% of the time vs 1.0% by chance.
- `model.potency.PotencyTable` + `tc.selectivity`: potency strip, exposure band or stated 10-fold rule, censored values, mixed-assay refusal, panel coverage.
- `ImageMark`, `ReferenceLibrary.annotations`, `library_from_compounds`, SHA-512 checksums.
- D7 (NanoBRET) and D8 (Davis 2011) downloads blocked by publisher access controls; recorded with manual steps.

### M4: v0.1 polish (2026-11-17)
- One-command CS1 (`python case_studies/cs1_sciplex_hdac.py`) verified from a fresh clone with the core install; the sci-Plex cache builds on first use.
- Hallmark gene sets are fetched directly from OmniPath (no decoupler needed); identical to decoupler's table.
- `docs/gallery.md` (classic vs tincture for every plot), README before/after, generated `docs/api.md` with a staleness test.
- User-study kit: `case_studies/user_study.py` (pairs, answer key from results, counterbalanced assignments, scorer) and `docs/user_study.md`.
- Classic fallbacks are fairer: Compass PCA colours reference classes; PCA by dose marks and labels compounds.
- CS1 uses a pharmacological HDAC group (pan / class-I), fixed before re-running; writes `results.json`.
- Docs image with GDSC panels cropped to its sci-Plex half (licence practice).

### M3: Compass, Dose Paths, Fingerprint Strip, Effect Ladder (2026-11-17)
- `model.reference`: `ReferenceLibrary` (class-mean references, leave-one-out, coherence) and `library_from_signatures`.
- `stats.projection`: NNLS on unit references with noise-corrected magnitude and R², χ² signal test, bootstrap ellipses, seriation.
- `tc.compass`: coherence-chosen spokes, leave-one-out, rim clipping, gate, novelty only with signal, greedy labels; PCA fallback.
- `stats.paths` and `tc.dose_paths`: stress axis (Hallmark apoptosis + p53 + UPR) or cell loss; bend test on scaled axes (turn > 30° beyond noise, Bonferroni); class-mean reference path; small multiples above 8; `reference=`; PCA-by-dose fallback.
- `tc.fingerprint`: ULM scores, BH over the matrix, curated Hallmark order and groups, references first, stated row ordering, Jaccard dedup, focus tick; volcano fallback.
- `tc.effect_ladder`: top-N forest grouped by gene set, CI colouring, transcriptome-wide count strip, N and shrinkage printed, `pathway=` drill-down; volcano fallback.
- Spec/renderer: heatmap column gaps, rotated/top tick labels, unclipped rects, colormapped points with a colourbar.
- `Signature.align` vectorised (about 5× faster on sci-Plex).
- Case study CS1: sci-Plex3 A549 HDAC inhibitors, gate to genes.

### M2: statistics, Retrievability, Plate Atlas (2026-11-17)
- `stats.ebayes`: limma-style variance moderation (`fit_f_dist`, `squeeze_var`, with an optional trend).
- `stats.signatures`: logCPM signatures vs same-stratum vehicle, moderated t and normal-prior shrinkage; `PerturbationSet.signatures()`; label-free `replicate_profiles`.
- `stats.pathways`: ULM pathway scores (matches decoupler to 1e-12) and Jaccard de-duplication; gene sets (Hallmark, PROGENy) cached with checksums.
- `stats.map_retrieval`: AP, mAP, a label-permutation null (exact when feasible) and the copairs "copied" null (validated against a copairs 0.6.0 golden file), BH q-values, CAT curve.
- `tc.retrievability`: gate plot with guards; flags designs that cannot reach significance; gate mask stored for other plots.
- `stats.zprime` and `tc.plate_atlas`: Z′, layout-aware edge effect, drift within cell line; density view for 1536-well plates.
- Palette: colour-blind-safe muted purple/green; GDSC licence terms recorded.
- Case studies: CS3 (GDSC2 plate QC) and the sci-Plex Retrievability validation.

### M1: data model and loaders (2026-11-16)
- `tincture.model`: `validate_long`, `PerturbationSet` (AnnData), `Signature`, `FourPL`, `DoseResponse`, `DoseResponseSet`, `PlateTable`; `well_type` gains `blank`.
- `tincture.datasets`: checksum-verified download cache and a registry of 13 pinned files.
- Loaders: sci-Plex3 (streamed per-well pseudobulk), PRISM 20Q2 secondary (points plus published fits), GDSC2 raw release 8.5 (1536-well `PlateTable`).
- `docs/data_status.md`: every "verify" fact for D1–D3 checked against the source; differences recorded.
- 158 unit tests plus 4 opt-in real-data round-trips.

### M0: scaffold (2026-11-16)
- Packaging (hatchling) with core, `[omics]`, `[chem]`, `[combo]` and `[dev]` extras; ruff, mypy and pytest-mpl; GitHub Actions CI.
- `PlotSpec` dataclasses and the matplotlib renderer `render_mpl`.
- Theme: palette, diverging and sequential map rules, typography, journal presets, `check_colorblind`.
- Captions (templated), methods paragraphs and the classic-fallback registry.
- `tc.export` / `PlotResult.save` with 6 pt and colour-blind checks and a markdown sidecar.
- `tc.demo` end-to-end demo plot; 104 tests, including 2 image baselines.
- Licence changed to Apache-2.0.
