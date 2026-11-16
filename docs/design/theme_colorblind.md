# Theme: colour-blind check and a palette finding

## How the check works

`tc.theme.check_colorblind(fig)`:

1. Collects the figure's **categorical** colours. Artists drawn by the tincture renderer are tagged (`gid="tc:cat:<hex>"`) when their mark has `categorical=True`. For arbitrary matplotlib figures, it falls back to line, patch and unmapped collection colours. Greys, black, white, fully transparent and colormapped artists are ignored.
2. Simulates **deuteranopia** and **protanopia** with the Machado, Oliveira & Fernandes (2009) matrices at severity 1.0, applied in linear RGB.
3. Computes **CIEDE2000** between every pair, under normal vision and both simulations. The implementation is tested against the Sharma, Wu & Dalal (2005) reference pairs.
4. Fails if any pair is below **ΔE₀₀ = 10**, the default and configurable threshold.

`tc.export` runs this check on every export (main and classic figures). Pass `check=False` to override it deliberately.

## Finding: two pairs in the brief's palette collide (open decision)

| Pair | Meaning | Worst ΔE₀₀ | Under |
|---|---|---|---|
| blue `#2F7CC4` vs purple `#7C5CBF` | reference vs novel | **2.1** | deuteranopia (7.2 protanopia) |
| terracotta `#D97757` vs green `#2E8B57` | focus vs good/window | **2.1** | protanopia |
| light blue `#9CCDE8` vs light purple `#CDB6EC` | reference vs novel (light) | 1.6 | deuteranopia |
| peach `#F6BD78` vs light green `#B7DFA4` | focus-light vs window | 7.4 | protanopia |

Reference vs novel is the central contrast of the Compass and Neighbourhood Card, and focus vs window matters in the Window plot. As specified, these pairs are indistinguishable to roughly 1 in 12 men.

**The palette is implemented exactly as specified.** The checker refuses any figure that puts a colliding pair together, so plots in M3 onward must avoid those combinations.

Proposed replacement (not applied):

| Role | Current | Proposed | Min ΔE₀₀ vs blue/terracotta/other, all three visions |
|---|---|---|---|
| novel | `#7C5CBF` | `#55309A` (deeper indigo-purple) | ≥ 15.8 |
| good | `#2E8B57` | `#1E5C3E` (deeper green) | ≥ 15.8 |

With both swaps, {blue, terracotta, purple, green} stay at ΔE₀₀ ≥ 15 under normal vision, deuteranopia and protanopia. The light variants need the same treatment once the dark ones are settled.
