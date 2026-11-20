# Theme: colour-blind check and a palette finding

## How the check works

`tc.theme.check_colorblind(fig)`:

1. Collects the figure's **categorical** colours. Artists drawn by the tincture renderer are tagged (`gid="tc:cat:<hex>"`) when their mark has `categorical=True`. For arbitrary matplotlib figures, it falls back to line, patch and unmapped collection colours. Greys, black, white, fully transparent and colormapped artists are ignored.
2. Simulates **deuteranopia** and **protanopia** with the Machado, Oliveira & Fernandes (2009) matrices at severity 1.0, applied in linear RGB.
3. Computes **CIEDE2000** between every pair, under normal vision and both simulations. The implementation is tested against the Sharma, Wu & Dalal (2005) reference pairs.
4. Fails if any pair is below **ΔE₀₀ = 10**, the default and configurable threshold.

`tc.export` runs this check on every export (main and classic figures). Pass `check=False` to override it deliberately.

## Finding: the brief's palette had collisions, now fixed (M2)

As specified, two pairs were near-identical for colour-blind readers (ΔE₀₀ ≈ 2): **reference blue `#2F7CC4` vs novel purple `#7C5CBF`** under deuteranopia, and **focus terracotta `#D97757` vs good green `#2E8B57`** under protanopia. The light variants collided too (light blue vs light purple ΔE₀₀ 1.6).

**Decision (approved by the user: "colour-blind friendly, but attractive, no bright high-contrast pinks"):** keep the brief's blue, terracotta, peach, light blue and neutral, which carry the diverging scale and already pass. Re-choose only the purple and green families, in muted tones:

| Role | Brief | Now | Character |
|---|---|---|---|
| novel | `#7C5CBF` | **`#4B3F70`** | dusty indigo |
| novel (light) | `#CDB6EC` | **`#9993AC`** | lavender grey |
| good | `#2E8B57` | **`#3B5F44`** | forest |
| window / good (light) | `#B7DFA4` | **`#B6D7C1`** | soft sage |

How they were chosen: random search in CIE LCh, constrained to violet/indigo hues (280–305°, no magenta or pink), sage/forest greens (120–175°), and chroma no higher than the brief's own blue. Each candidate was scored by its worst-case CIEDE2000 under normal vision, deuteranopia and protanopia.

Result:
- dark set {blue, terracotta, indigo, forest}: min ΔE₀₀ = **27.4 / 18.0 / 18.4** (normal / deutan / protan);
- light set: min **20.2 / 13.6 / 14.9**;
- all 8 together: ≥ 10 in every vision type.

A livelier purple (`#5A4A8E`) was rejected because it fell to ΔE₀₀ 8.8 against light blue under deuteranopia.

Visual comparison (old vs new, with a deuteranopia simulation): `docs/img/palette_compare.png`. Tests: `test_tincture_palette_is_colourblind_safe`, `test_palette_is_muted` (chroma ≤ 50), and the original colours are kept as `BRIEF_PURPLE`/`BRIEF_GREEN` for regression tests.
