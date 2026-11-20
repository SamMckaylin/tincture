# Milestone 0: scaffold, theme, PlotSpec, renderer, export, captions

## What was built

- **Packaging:** `pyproject.toml` (hatchling). Core deps: numpy, pandas, scipy, matplotlib, anndata, statsmodels. Extras: `[omics]` (pydeseq2, decoupler, scanpy), `[chem]` (rdkit, mmpdb), `[combo]` (empty) and `[dev]`. Tooling: ruff (lint and format), mypy (non-strict, `check_untyped_defs`), pytest and pytest-mpl. `data/` and `case_studies/out/` are git-ignored.
- **CI** (`.github/workflows/ci.yml`) has three jobs:
  - lint: ruff check, ruff format --check, mypy;
  - test: Python 3.10 and 3.12, pixel comparison on 3.12 with matplotlib pinned to the baseline version;
  - omics-extra: installs `.[omics,dev]`, imports, tests.
- **`tincture.spec`** has plain dataclasses only:
  - Marks: Point, Line, Interval, Band, Rect, Heatmap, Contour, Ellipse, RefLine, Text, Notice.
  - Supporting classes: `Scale` (linear, log, symlog, categorical, polar angle/radius), `Panel`, `PlotSpec` (grid layout, ratios, caption, methods, footer, warnings, meta), legend and colourbar specs.
  - `PlotSpec.validate()` catches inconsistent specs.
  - Colours are semantic names ("focus", "reference", "novel" …).
- **`tincture.render.mpl.render_mpl`** is the only code that creates matplotlib artists. It:
  - renders at the exact preset width;
  - starts from matplotlib defaults plus the theme, so the user's own style can't leak in;
  - pins tick styling so ticks created late stay on theme;
  - tags categorical artists for the colour-blind checker;
  - reserves footer space for reference-set, version and n lines.
- **`tincture.theme`:**
  - palette and semantic colours from §3; 5-step diverging map and `BoundaryNorm`;
  - sequential maps restricted to viridis, cividis, magma, inferno and plasma (jet, rainbow, turbo etc. raise `TinctureStyleError`);
  - typography; journal presets;
  - `check_colorblind(fig)`: Machado 2009 simulation and CIEDE2000.
- **`tincture.report`:**
  - `Caption`: sentences are dropped when a statistic is missing; `.weak()` notes always print; the §5 templates are stored in `TEMPLATES`.
  - `Methods`: reusable blocks with citations plus software versions.
  - Classic-fallback registry.
  - Formatters `fmt_p` and `fmt_conc`.
- **`tincture.export` / `PlotResult.save`:** writes the figure (pdf, png, svg, tiff or eps), a markdown sidecar (caption, warnings, methods, footer, preset) and the classic fallback. It refuses text under 6 pt and colour-blind collisions.
- **`tc.demo`:** a dummy plot that exercises the whole pipeline end to end. It reads the canonical long table and builds no data itself.
- **Errors:** `TinctureDataError(message, fix=...)` always carries a fix suggestion. There are also `TinctureSpecError`, `TinctureStyleError` and `TinctureColorblindError`.
- **Licence:** Apache-2.0, replacing the GPL-3.0 placeholder as you asked me to choose. Apache-2.0 is permissive, which suits pharma and scverse users, and has an explicit patent grant. The text was copied from an installed package whose file matches the canonical apache.org MD5 (`3b83ef96…`).

## Acceptance criteria

| Criterion | Status |
|---|---|
| `pip install -e .[omics]` works | ✅ locally (pydeseq2 0.6.0, decoupler 2.3.1, scanpy 1.12.6); CI job added |
| A dummy plot exports at `nature_double` size | ✅ `tests/test_export.py::test_demo_exports_at_nature_double`: PNG width = 183 mm at 600 dpi, PDF MediaBox = 518.74 pt |
| CI green | ✅ locally (lint, format, mypy, 104 tests on Python 3.12 and 3.10). GitHub run: see the commit's checks |

## What was verified

- CIEDE2000 matches the Sharma et al. (2005) reference values to 1e-4.
- Output is byte-identical under a global `classic` style and the default style (via export).
- Late-materialised tick labels stay at or above 6 pt even when global rcParams say 3 pt.
- The suite passes on Python 3.10 (matplotlib 3.10.9; pixel comparison disabled, see below) and 3.12 (matplotlib 3.11.2, pixels compared).

## What failed or surprised

1. **Palette collisions** (details in `docs/design/theme_colorblind.md`). In the palette as specified, reference blue vs novel purple has ΔE₀₀ 2.1 under deuteranopia, and focus terracotta vs good green has ΔE₀₀ 2.1 under protanopia. The palette is kept as specified and the checker refuses figures that combine those pairs. A replacement (`#55309A`, `#1E5C3E`) is proposed. **Needs your decision.**
2. **Journal presets are unverified.** The publishers' sites are not reachable from this environment, so the presets are from memory (`docs/design/journal_presets.md`).
3. **Pixel baselines depend on the matplotlib version.** matplotlib 3.11 dropped Python 3.10, and 3.10.x lays out text differently (RMS 23–34). Baselines record their matplotlib version; comparison runs only when major.minor matches, and figures are still built and colour-blind checked otherwise.
4. The `[omics]` extra is about 600 MB installed (pyarrow, llvmlite/numba, scikit-learn). You approved it.

## Open questions

- Palette replacement for novel and good (point 1 above).
- §9 items still open: package name (PyPI check), default synergy model (Bliss vs HSA), stress axis for Dose Paths, `copairs` vs re-implementing mAP. None of these block M1.
- Decided this milestone: licence (Apache-2.0); standalone package that reads AnnData rather than a pertpy extension; sci-Plex limited to the sciPlex3 subset.
