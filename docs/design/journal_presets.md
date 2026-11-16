# Journal export presets

| Preset | Width | Max height | Notes |
|---|---|---|---|
| `nature_single` | 89 mm | 247 mm | |
| `nature_double` | 183 mm | 247 mm | M0 acceptance target |
| `cell_single` (`cell`) | 85 mm | 225 mm | |
| `cell_onehalf` | 114 mm | 225 mm | |
| `cell_double` | 174 mm | 225 mm | |
| `acs_single` (`acs`) | 82.55 mm (3.25 in) | 241.3 mm (9.5 in) | |
| `acs_double` | 177.8 mm (7 in) | 241.3 mm | |
| `screen` | 160 mm | 400 mm | notebooks/slides; 8 pt base, 150 dpi |

All journal presets use a 7 pt base font, a **6 pt minimum** (enforced at export), 600 dpi raster output, and TrueType-embedded text (`pdf.fonttype = 42`) so text stays editable. Fonts are tried in order: Arial, Helvetica, Liberation Sans (metric-compatible with Arial), then DejaVu Sans.

**Verification status: not verified.** These numbers are from memory of the publishers' author guidelines. The attempt to re-check them online during M0 failed because nature.com, cell.com and pubs.acs.org are not reachable from the build environment. Re-check before a real submission; maximum heights are the least certain values.

Output is independent of the user's matplotlib style. Rendering starts from matplotlib defaults and then applies the theme, and `tc.export` saves inside the same context. A test confirms that a global `classic` style produces byte-identical exports.
