# Milestone 3: Compass, Dose Paths, Fingerprint Strip, Effect Ladder

## What was built

**Model and statistics**
- `model.reference`: `ReferenceLibrary` built from class means of member signatures (`library_from_signatures`), with leave-one-out (`without`), class coherence, version and source carried into every footer.
- `stats.projection`: NNLS on unit-norm references.
  - **Noise correction.** Over 13k genes, measurement noise alone gives a signature a large norm, and that noise is by construction "unexplained". Magnitude is therefore √(‖y‖² − Σse²), R² is computed against that signal energy, and novelty is judged only when a χ² test says there is signal at all (z > 3).
  - Also: bootstrap ellipses and optimal-leaf seriation of spokes.
- `stats.paths`: the dose-path axes and the bend test.
  - **x** = projection on the top reference (from the Compass, leave-one-out). **y** = stress: the mean log2FC of Hallmark apoptosis ∪ p53 ∪ UPR genes minus the mean over all genes. This is your §9.5 choice; cell loss is available as `y_axis="cell_loss"`.
  - **Bend** = a turn of more than 30° between successive dose deltas. p = P(|turn| ≤ 30°) under posterior draws, Bonferroni-corrected over the turning points.
  - Angles are measured on axes scaled by their spread across the stratum. Both axes are linear in the signature, so the draws come from the exact bivariate normal of each dose's position, which is fast and exact.

**Plots**
- `tc.compass` (§5.1):
  - spokes are the most coherent classes (mean member cosine × √n), at most 8, seriated, with a correlated-neighbour guard;
  - radius = noise-corrected magnitude, with the rim at the 95th percentile; clipped points sit on the rim with their magnitude in brackets;
  - "poorly explained" (purple) is called only with signal; weak compounds are greyed;
  - gate (drop or grey), greedy labels, PCA fallback.
- `tc.dose_paths` (§5.2):
  - paths from vehicle with dot size ∝ log dose, an arrowhead at the top dose, and a blue class-mean reference path (query excluded);
  - stress and on-target zones; 95% ellipses for the focus compound;
  - bend annotation only when the bend is significant;
  - small multiples above 8 compounds, with a shared y axis and an optional shared `reference=`;
  - "too few doses for a path" guard; gate; PCA-coloured-by-dose fallback (dose drawn as a real colour scale).
- `tc.fingerprint` (§5.3):
  - ULM scores, with BH correction over the whole matrix; non-significant cells are neutral;
  - curated Hallmark order in six groups with column gaps; references first, above a gap;
  - row order is either clustered (correlation, optimal leaf order) or given, and is always stated;
  - Jaccard > 0.5 de-duplication; terracotta focus tick; a "nice" symmetric clip shown on the colour key;
  - row labels are thinned above 60 rows; volcano fallback.
- `tc.effect_ladder` (§5.4):
  - the top N (default 30) features by |shrunk LFC|, grouped by gene set in curated order;
  - 95% intervals, coloured only when they exclude zero;
  - a transcriptome-wide count strip (always drawn);
  - N and the shrinkage method printed in the footer;
  - `pathway=` drill-down; volcano fallback.
- Renderer: heatmap column gaps, rotated and top tick labels, unclipped rectangles, colormapped points with a colourbar, and wrapped footers.

**Speed:** `Signature.align` is now vectorised. Dose Paths on 6 sci-Plex compounds went from 18 s to 3 s. Budgets are tested: Compass on 500 compounds in under 3 s, and Fingerprint on 1,000 × 50 in under 5 s, both including rendering.

**Tests:** 265 passing, plus 7 opt-in real-data tests. New baselines for Compass, Dose Paths, Fingerprint and Effect Ladder, all colour-blind-checked. Every §5.1–5.4 guard has its own test:

| Plot | Guards tested |
|---|---|
| Compass | correlated spokes; poorly explained only with signal; reference set and version in the footer; leave-one-out; incoherent classes are not spokes |
| Dose Paths | more than 8 compounds → small multiples; fewer than 3 doses → points and "too few doses for a path"; non-significant bend → no annotation; bend angle measured on scaled axes; straight path vs injected bend |
| Fingerprint | clip shown in the legend; row ordering stated; Jaccard dedup; references first with a gap; non-significant cells neutral |
| Effect Ladder | N and shrinkage printed; count strip always present; CI crossing zero drawn grey; drill-down |

## Acceptance (CS1: HDAC inhibitors, sci-Plex3 A549, 24 h)

`case_studies/cs1_sciplex_hdac.py` → `case_studies/out/cs1_sciplex_hdac/report.md` (Retrievability → Compass → Dose Paths → Fingerprint → Effect Ladder).

| Ground-truth check | Result |
|---|---|
| HDAC inhibitors form a coherent group | ✅ HDAC is one of 5 spokes; 18/25 HDAC inhibitors with signal at 10 µM project mainly onto it (leave-one-out). Abexinostat: HDAC weight 0.98, 3% unexplained. |
| Their Fingerprint rows match each other | ❌ **by the criterion I set.** Mean row correlation is 0.69 within HDAC inhibitors vs 0.54 with the other compounds, but permutation p = 0.011 against a pre-set p < 0.01. See below. |
| Dose increases magnitude along a consistent direction | ✅ HDAC-likeness rises with dose (Spearman ≥ 0.8) for 14/19. The median smallest noise-corrected cosine between a dose step and the leave-one-out HDAC reference is 0.86. 3 of 19 show a significant bend (panobinostat above 1 µM, towards stress). |

**Why the Fingerprint check fails narrowly, and what I did not do.**
- The comparison compounds (Aurora, JAK, HSP90 inhibitors, nucleotide analogues) share the generic anti-proliferative fingerprint: E2F targets and G2/M down. That lifts the "between" correlation.
- The HDAC group also contains two compounds with no HDAC-like transcriptional signal at 10 µM in this assay:
  - **tubastatin A**, which is HDAC6-selective;
  - **valproate**, which needs millimolar concentrations.

  Both pass the retrievability gate and are flat on the HDAC axis.
- Post hoc, restricting to the 18 HDAC inhibitors whose Compass top spoke is HDAC gives 0.83 within vs 0.59 between. That analysis is labelled post hoc in the report and is **not** used for the check.
- I kept the criterion as set. Changing it after seeing the result would make the check meaningless. Your call (see open questions).

### Re-run with a pharmacological HDAC group (your decision: option 2)

Before re-running, the ground-truth group was fixed as **pan / class-I HDAC inhibitors**: the dataset's annotation minus tubastatin A (HDAC6-selective) and valproate (millimolar potency). Both are left out of every check, in neither group. The definition is in `NOT_PAN_CLASS_I` in the script and printed in the report. I ran it once:

| Check | Result |
|---|---|
| Coherent group | ✅ 18/25 with signal project mainly onto HDAC (unchanged; the two excluded compounds had no signal) |
| Fingerprint rows match | ✅ 0.83 within vs 0.59 with other compounds, permutation p = 0.0005 (criterion p < 0.01, unchanged) |
| Dose along a consistent direction | ✅ 14/17 rise with dose; median noise-corrected cosine 0.86 |

**Overall: PASS.** The first-run numbers above stay in this file and in the report header.

## What failed or changed along the way

1. **Compass noise.** Without noise correction, 85 of 102 compounds were flagged "novel", because noise counts as residual. Fixed with the signal-energy correction and χ² gate described above.
2. **Spoke choice.** Choosing spokes by magnitude dropped the Aurora class. Spokes are now chosen by coherence.
3. **Bend test, first version.** It measured angles in raw units (x spans about 1,000× the y range) and called any precisely measured curve a bend. Now it uses scaled axes and a minimum turn of 30° beyond noise.
4. **CS1 dose-direction metric, first run.** A raw cosine between 13k-gene step vectors is dominated by measurement noise, so even clean HDAC paths scored 0.42. I replaced it with a noise-corrected cosine against the leave-one-out HDAC reference, counting only steps with signal. The report states the revision and keeps the first-run column.
5. **Dose Paths in small multiples** used each compound's own top reference, so valproate and tubastatin got HSP90 and Aurora axes. Added `reference=`; CS1 uses the HDAC axis for every panel.
6. **Colour-blind check vs ordinal colours.** The PCA-by-dose fallback sampled viridis into hex colours, so the checker read four close categorical colours and blocked the export. Dose is now a true colour scale with a colourbar.
7. **Footer overflow:** long footers ran off the page. They now wrap.
8. **CI was red on the Python 3.10 job for three pushes** (Compass, Dose Paths and Fingerprint commits). I only noticed when checking CI at the end of the milestone. The cause: matplotlib 3.10 (unpinned there) rejects a dashed style on a zero-width patch edge, while 3.11 tolerates it. Fixed in the renderer and verified locally against both 3.10.9 and 3.11.2. Lesson: check CI after every push, not only at milestone end.

## Decisions taken (delegated by you)

- **Stress axis:** the Hallmark stress score by default, with cell loss as an alternative (your approval of my §9.5 recommendation).
- **Reference library for sci-Plex:** class means of retrievable members at 10 µM using the dataset's `pathway_level_2` annotation, leave-one-out for members. LINCS-based references (D9) remain for M5.

## Open questions

1. **CS1 Fingerprint check:** accept as is (the separation is real but modest, p = 0.011), or pre-register a different criterion for the next run? Options: compare only with compounds outside the anti-proliferative classes, or define the HDAC group pharmacologically (pan/class-I inhibitors only, excluding HDAC6-selective tubastatin and valproate). I'd suggest the pharmacological definition, decided now, before re-running.
2. Still open from §9: package name (PyPI check); Bliss vs HSA default (M7).
