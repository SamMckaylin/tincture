# Changelog

## Unreleased

### M0: scaffold (2026-11-16)
- Packaging (hatchling) with core, `[omics]`, `[chem]`, `[combo]` and `[dev]` extras; ruff, mypy and pytest-mpl; GitHub Actions CI.
- `PlotSpec` dataclasses and the matplotlib renderer `render_mpl`.
- Theme: palette, diverging and sequential map rules, typography, journal presets, `check_colorblind`.
- Captions (templated), methods paragraphs and the classic-fallback registry.
- `tc.export` / `PlotResult.save` with 6 pt and colour-blind checks and a markdown sidecar.
- `tc.demo` end-to-end demo plot; 104 tests, including 2 image baselines.
- Licence changed to Apache-2.0.
