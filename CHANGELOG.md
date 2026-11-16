# Changelog

## Unreleased

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
