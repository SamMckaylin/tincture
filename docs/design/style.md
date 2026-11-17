# Visual style (aligned with the design whiteboard)

The concept whiteboard (https://claude.ai/artifact/2agybZq3fFVDpULaPx1qwS) sets the look. After M2 the theme and renderer were brought in line with it. These rules apply to every plot.

| Rule | Implementation |
|---|---|
| **Typeface:** Inter, regular-weight titles | `FONT_STACK = Inter → Arial → Helvetica → Liberation Sans → DejaVu Sans`; `axes.titleweight = normal`. For journals that require Arial or Helvetica: `tc.theme.set_font("Arial")`. Tests pin DejaVu Sans so baselines don't depend on installed fonts. |
| **Direct labels, no legend boxes** | `LineMark.end_label` (text at the line end in its colour), `BandMark.inner_label`, `TextMark` labels placed by the plot builder; point labels take the point's colour. All direct labels get a thin paper-coloured halo. |
| **Bands:** light fill, dashed outline in a darker shade | `BandMark(edgecolor=..., edge_style="dashed")`, as in the whiteboard's window and novelty halos. |
| **Cells:** rounded tiles with gutters | `HeatmapMark(style="auto")` draws rounded tiles (corner 0.2, gutter 0.14 cell) up to 5,000 cells, flat cells above that. `circles=` draws chosen cells as circles (controls on plate maps). |
| **Not significant ≠ not measured** | `mask` → neutral grey; NaN → empty outline. The two can never look alike. |
| **Points:** crisp | Filled points get a darker edge of their own colour; dense runs (> 150 points) skip it so they don't merge into a dark band. |
| **Quiet axes** | Thin ink spines; short mid-grey (`#878787`) ticks; soft tick labels; no ticks on tile maps. |
| **Colour keys:** small | `ColorbarSpec(compact=True, ticks=..., tick_labels=...)`, e.g. the plate map's "÷2 · same · ×2" key, like the whiteboard's "down ▪▪▪▪▪ up". |
| **Semantic colours** | Terracotta = focus; blue = reference / signal down; peach and terracotta = up; mid-grey = vehicle; green = good / window (Z′ bars); purple = novel. Purple and green are the colour-blind-safe variants (`docs/design/theme_colorblind.md`). |

## Plate Atlas vs the whiteboard sketch

The sketch shows grey compound wells, light-blue edge wells (evaporation), terracotta and peach anomalies, controls as circles, and green Z′ bars with a terracotta failing plate and a dashed threshold. The implementation follows it:
- tiles coloured by log₂(signal ÷ plate's negative-control median) on the 5-step diverging scale (neutral = ±15%);
- controls as circles; wells absent from the source as empty outlines;
- Z′ bars in window green, failing plates in terracotta;
- above 1,000 plates, a binned median line and 5–95% band, with every failing plate still drawn as a point.

## Deliberate deviations

- **Retrievability** has no whiteboard sketch. It uses the same language: direct labels and a dashed no-phenotype band.
- **Novel and good** colours differ from the whiteboard's `#7C5CBF` and `#2E8B57` for colour-blind safety (user-approved).
