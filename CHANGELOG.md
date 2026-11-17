# Changelog

## Unreleased

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
