# Tincture: build handover for Claude Code

> **What this is.** A complete brief for building `tincture`, a Python package for question-first, publication-quality visualisation of drug-discovery data, and for validating every plot on real public data.
> **Who it's for.** Claude Code, working in a fresh repo. Put this file in the repo root (or rename it `CLAUDE.md`) and start with Milestone 0.
> **Design reference.** The concept whiteboard (example sketches of every plot) lives at https://claude.ai/artifact/2agybZq3fFVDpULaPx1qwS. You may not be able to open it; everything you need is in this file.

---

## 0. Ground rules for the build

1. **Work milestone by milestone** (section 8). Finish a milestone's acceptance criteria and tests before starting the next. Commit at the end of each milestone, with a short changelog entry.
2. **Never invent data.** Every case study uses real downloaded data. If a download fails (blocked network, moved URL, licence click-through), stop that case study, write what failed in `docs/data_status.md`, and ask the human to place the file in `data/raw/<dataset>/`. Do not substitute simulated data silently. Simulated data is fine **only** in unit tests and must live in `tests/fixtures/`, clearly named `synthetic_*`.
3. **Never commit raw data.** `data/` is git-ignored. Write download scripts with checksums and a local cache. Several sources (GDSC especially) allow internal research use but forbid redistribution.
4. **Verify dataset facts before coding loaders.** Section 6 marks each fact as *verified* (checked against the source during planning) or *verify* (from memory). For anything marked *verify*, open the source, confirm the file names, columns and contents, and record what you found in `docs/data_status.md`.
5. **Ask before adding heavy dependencies** (anything that pulls in R, Java, CUDA, or more than ~50 MB). Core dependencies are listed in section 3.
6. **Every plot ships with three things:** the figure, an auto-generated reading-guide caption (templated from computed statistics, no LLM), and a methods paragraph. Plus a `classic=True` fallback where one exists (see each spec).
7. **Red-team guards are requirements, not suggestions.** Each plot spec lists failure modes. Each guard needs a test.

---

## 1. Goal and principles

**Problem.** Drug-discovery papers recycle the same generic plots (PCA, UMAP, volcano, GSEA bars). They show the structure of the data, not answers to drug-discovery questions.

**Thesis.** Organise the package around the questions drug discovery asks, with the right statistics underneath:

| Question | Plot(s) |
|---|---|
| What does this compound look like? | Compass, Neighbourhood Card |
| Does the mechanism change with dose? | Dose Paths |
| What did it do biologically? | Fingerprint Strip, Effect Ladder |
| Is there a safe, effective dose? | Window |
| How specific is it? | Selectivity Spectrum |
| Does potency hold up from enzyme to cell to phenotype? | **Engagement Ladder** |
| Do two drugs work better together, at reachable doses? | Combination Surface, Combination Lift |
| What happens first, and what follows? | Onset Timeline |
| Is the effect real? | Retrievability (replicate concordance), Plate Atlas |
| How did chemistry change biology? | SAR Arrows |
| Where did compounds drop out? | Triage Flow |

**Principles (apply to every plot):**
1. Question-first API: `tc.selectivity(...)`, not `tc.scatter(...)`.
2. Dose and time are first-class axes.
3. Always show a reference (vehicle, known-MoA compound, or comparator).
4. Uncertainty is on by default.
5. Every axis has biological or pharmacological meaning.
6. Opinionated, colour-blind-safe theme; journal-sized export presets.
7. Degrade gracefully: missing dose, single replicate, NaNs give a clear message or a reduced plot, never a crash or a misleading figure.
8. **Say when a result is weak.** If a fit fails, a neighbour is not significant, or a bend is within noise, the plot shows that explicitly instead of drawing a confident-looking mark.

---

## 2. Scope

**v0.1 (this handover, Milestones 0–4):** data model, theme, export, captions; Compass, Dose Paths, Fingerprint Strip, Effect Ladder, Retrievability, Plate Atlas; case studies on sci-Plex, PRISM and GDSC monotherapy raw data.

**v0.2 (Milestones 5–7):** Window, Selectivity Spectrum, Neighbourhood Card, Engagement Ladder, Combination Surface and Lift, Onset Timeline; case studies on GDSC combinations, NanoBRET/Davis, LINCS, MIX-seq.

**v0.3 (backlog, do not start):** SAR Arrows, Triage Flow, Cell-line Atlas, Figure Audit, interactive backend, R bindings.

**Out of scope:** docking, ADMET prediction, structural biology, any model training beyond curve fitting and simple linear models.

---

## 3. Architecture and stack

- **Python ≥ 3.10.** Packaging with `pyproject.toml` (hatchling), `ruff` for lint and format, `mypy` (non-strict) on `tincture/`, `pytest` + `pytest-mpl` for tests.
- **Core dependencies:** `numpy`, `pandas`, `scipy`, `matplotlib`, `anndata`, `statsmodels`.
- **Optional extras:** `[omics]` = `pydeseq2`, `decoupler`, `scanpy`; `[chem]` = `rdkit`, `mmpdb`; `[combo]` = none beyond core.
- **Rendering:** matplotlib only for v0.1/v0.2, but every plot builds an intermediate **PlotSpec** (a plain dataclass of marks, scales, annotations) and a separate `render_mpl(spec)` draws it. This keeps an interactive backend possible later without rewriting the plots. Do not leak matplotlib objects into the stats layer.

```
tincture/
  __init__.py          # public API: tc.compass, tc.dose_paths, ...
  io/                  # loaders: anndata, deseq/limma tables, gdsc, prism, sciplex, lincs, nanobret, plate csv
  model/               # PerturbationSet, Signature, ReferenceLibrary, DoseResponse, ComboMatrix,
                       # PotencyTable, TimeCourse, PlateTable, StageLog
  stats/               # signatures, shrinkage, hill (4PL), projection (NNLS), similarity + nulls,
                       # synergy (Bliss/HSA/Loewe), map_retrieval, onset, zprime
  spec/                # PlotSpec dataclasses
  plots/               # one module per plot: build_spec(...) -> PlotSpec
  render/mpl.py        # PlotSpec -> matplotlib Figure
  theme/               # palettes, typography, journal presets (nature_single, nature_double, cell, acs)
  report/              # captions.py (templated), methods.py, classic.py (fallback plots)
  datasets/            # download + cache helpers with checksums (no data in repo)
tests/
  fixtures/            # tiny synthetic_* and tiny real subsets (only if licence allows; else synthetic)
  baseline/            # pytest-mpl image baselines
case_studies/          # one script per case study; writes figures + a short markdown report
docs/                  # data_status.md, design notes, gallery
```

### Theme

- **Diverging (down / up):** blue `#2F7CC4` / `#9CCDE8`, neutral `#D1D1D1`, peach `#F6BD78` / terracotta `#D97757`.
- **Focus colour:** terracotta `#D97757` is reserved for "the compound or thing the figure is about".
- **Reference / known-MoA:** blue. **Novel / unexplained:** purple `#7C5CBF` / `#CDB6EC`. **Good / window:** green `#2E8B57` / `#B7DFA4`.
- No rainbow or jet colormaps. Sequential maps from `matplotlib`'s perceptually uniform set only.
- Provide a `tc.theme.check_colorblind(fig)` utility that simulates deuteranopia/protanopia and fails if any two categorical colours used in the figure fall below a ΔE threshold.

---

## 4. Data model

One canonical long table covers most plots. Loaders convert every source into these objects.

### 4.1 Measurement table (long format)

| column | type | notes |
|---|---|---|
| `compound_id` | str | required; stable ID (prefer InChIKey when available) |
| `compound_name` | str | optional |
| `smiles` | str | optional (needed for SAR, Neighbourhood structures) |
| `dose` | float | molar; `0` for vehicle |
| `time` | float | hours; optional |
| `replicate` | str/int | required for uncertainty |
| `plate`, `well`, `well_type` | str | `well_type` ∈ {compound, neg_ctrl, pos_ctrl, empty} |
| `cell_line` | str | optional |
| `feature` | str | gene, pathway, morphology feature, or "viability" |
| `value` | float | raw or normalised |
| `se` | float | optional |

### 4.2 Core objects

- `PerturbationSet` wraps `AnnData` (obs = samples with the metadata columns above, var = features). Methods: `.signatures(method=...)`, `.pseudobulk(by=[...])`, `.subset(...)`.
- `Signature`: effect vector + per-feature SE + metadata (compound, dose, time, cell line, n replicates). Bootstrap draws are stored, not just point estimates, so downstream CIs are honest.
- `ReferenceLibrary`: named signatures with MoA/target annotations, a version string and source citation. Printed on every figure that uses it.
- `DoseResponse`: per compound × assay; raw points + 4-parameter log-logistic fit (bottom, top, log EC50, Hill slope), fit quality (R², convergence flag), bootstrap CI on EC50.
- `PotencyTable`: compound × target (or rung) × value (pIC50/pKd/pEC50) × assay type × assay context (cell line, ATP concentration, readout).
- `ComboMatrix`: drug A conc × drug B conc × viability, with single-agent rows/cols, replicates, design flag `full` (e.g. 7×7) or `anchored` (2×7).
- `TimeCourse`: signatures indexed by time (and optionally dose).
- `PlateTable`: well-level values with plate geometry (96/384/1536).
- `StageLog`: compound × stage reached × drop reason.

---

## 5. Plot specifications

For each plot: **question → inputs → statistics → encoding → red-team guards (with tests) → caption template → classic fallback.**

### 5.1 Compass (replaces PCA)
- **Question:** what known mechanism does this compound resemble, and how strongly?
- **Inputs:** `Signature`s for query compounds; `ReferenceLibrary` (same platform and ideally cell line).
- **Statistics:** fit each query as non-negative least squares on reference signatures → weights `w_k`, reconstruction R². **Residual = 1 − R²** (the novelty signal). Position = Σ w_k·u_k where u_k is the unit vector at reference k's angle; radius scaled by signature magnitude (norm of shrunk effects). CI = bootstrap over replicates → ellipse.
- **Encoding:** polar plot; references on spokes, ordered by seriation of the reference–reference similarity matrix so correlated references sit next to each other; concentric gridlines; query dots; CI ellipses; dashed halo whose radius ∝ residual.
- **Guards:**
  - References are correlated → draw a thin outer ring showing ref–ref similarity between neighbours; warn if any pair has cosine > 0.8 ("these spokes are not independent").
  - Reference set printed in the figure footer with version.
  - Residual always shown; a compound with residual > 0.5 is labelled "poorly explained".
  - Polar area distortion → radius is linear in magnitude, not area; documented.
- **Caption template:** "{compound} most resembles {top_ref} (weight {w:.2f}); {pct_unexplained:.0%} of its signature is unexplained by the {n_refs} references."
- **Classic fallback:** PCA of the same signatures.

### 5.2 Dose Paths (replaces PCA coloured by dose)
- **Question:** does a compound's mechanism change with dose?
- **Inputs:** signatures at ≥ 3 doses per compound + vehicle; optional viability readout.
- **Statistics:** project each dose's signature into a 2-D space. Default axes: x = projection on the compound's top reference (from Compass), y = stress score (decoupler pathway score of a curated stress/apoptosis set, or viability loss). Alternative: user-supplied axes. **Bend test:** angle between successive dose deltas vs a bootstrap null from replicates; flag a bend only if p < 0.05.
- **Encoding:** one path per compound from vehicle to top dose; dot size ∝ log dose; arrowhead at top dose; shaded zones for "on-target" and "stress"; optional CI tubes from bootstrap.
- **Guards:**
  - More than 8 compounds → automatic small multiples.
  - Fewer than 3 doses → draw points, no path, and say "too few doses for a path".
  - Non-significant bend → no annotation.
- **Caption:** "{compound} tracks {ref} up to {dose_ok}; above {dose_bend} it moves toward {stress_axis} (bend p = {p:.2g})."
- **Classic fallback:** PCA coloured by dose.

### 5.3 Fingerprint Strip (replaces volcano + GSEA bar charts)
- **Question:** what did each compound do, at pathway level, compared across many compounds?
- **Inputs:** signatures; gene-set library (default MSigDB Hallmark via `decoupler`; PROGENy as alternative).
- **Statistics:** pathway activity per signature (`decoupler` ULM or MLM); per-score p-value; deduplicate gene sets with Jaccard > 0.5 (keep the larger).
- **Encoding:** rows = compounds (references first, separated by a gap), columns = pathways in a **fixed curated order** grouped (Cell cycle / Stress / Metabolism / Signalling / Immune / Other); diverging 5-step colour; non-significant cells drawn neutral; focus compound marked with a terracotta tick on the left.
- **Guards:** colour clipping limits shown in the legend; row ordering method stated; option to drill down (`tc.effect_ladder(compound, pathway=...)`).
- **Caption:** "{n} compounds × {m} pathways. {focus} is most similar to {nearest_row}; strongest effect: {pathway} ({direction})."
- **Classic fallback:** volcano for one selected compound.

### 5.4 Effect Ladder (replaces the volcano)
- **Question:** which features moved, and how confident are we?
- **Inputs:** per-feature effect + SE (from `pydeseq2` with LFC shrinkage, or any limma/DESeq2 table).
- **Statistics:** rank by |shrunk LFC|, top N (default 30), grouped by gene set.
- **Encoding:** forest plot; dashed zero line; CI bars; confident effects coloured (blue/terracotta), CIs crossing zero grey; small "count strip" at top showing total up/down significant features across the whole transcriptome (fixes the volcano's lost overview).
- **Guards:** N and shrinkage method printed; count strip mandatory.
- **Classic fallback:** volcano.

### 5.5 Retrievability (rethought "Replicate Concordance")
*Why rethought:* a replicate scatter with r = 0.9 is usually driven by a handful of large effects and says nothing about whether **each compound's** signature is reproducible. The question drug discovery actually needs answered is: "which compounds have a real, reproducible phenotype?" This plot is the **gate** before every other plot.
- **Inputs:** replicate-level profiles (transcriptomic or morphological).
- **Statistics:** for each compound, mean average precision (mAP) of its replicates retrieving each other among all profiles (cosine similarity), with a permutation null → p and q values. This is the approach used by the JUMP Cell Painting consortium's `copairs` framework (*verify* and use `copairs` directly if licence and dependencies allow; otherwise re-implement and test against it). Secondary panel: CAT curve (correspondence at the top: fraction of top-k features shared between two replicates vs k) against a null band.
- **Encoding:** strip of compounds sorted by mAP, with the null band and a q = 0.05 line; compounds below the line greyed.
- **Integration:** `PerturbationSet.retrievable()` returns a mask; other plots accept `gate=True` to grey out or drop non-retrievable compounds and say so in the caption.
- **Classic fallback:** replicate correlation scatter.

### 5.6 Plate Atlas (replaces plain plate heatmaps)
- **Question:** can we trust this plate?
- **Inputs:** `PlateTable`.
- **Statistics:** Z′ = 1 − 3(σ_pos + σ_neg)/|μ_pos − μ_neg| per plate; edge effect = median(edge) / median(interior); control drift across plate order.
- **Encoding:** plate grid (controls as circles, compounds as squares), side bars of Z′ per plate with a 0.5 line, failing plates in terracotta. Auto-switch to a density view for 1536-well plates.
- **Positioning:** QC and supplementary figure, not a main-figure plot.

### 5.7 Window (v0.2)
- **Question:** is there a dose that works without the toxicity?
- **Inputs:** `DoseResponse` for an on-target assay and ≥ 1 toxicity or counter-screen assay; optional free Cmax **range**.
- **Statistics:** 4PL fits; EC90(on) and EC10(tox) with bootstrap CIs; window exists if EC90(on) < EC10(tox).
- **Encoding:** shared log-concentration axis; both curves with raw points; green band = window; Cmax as a range band, not a line.
- **Guards:** if either fit fails (no convergence, R² < 0.8, or top < 50% effect), no band is drawn and the caption says why. Assay context (cell line, readout, duration) is printed on each curve.

### 5.8 Selectivity Spectrum (v0.2)
- **Inputs:** `PotencyTable` across a target panel, one assay type only.
- **Encoding:** potency strip sorted strongest first; free-Cmax band; intended target terracotta; off-targets above the band in peach; panel coverage ("n of ~518 kinases tested") in the footer.
- **Guards:** refuses mixed assay types (Kd vs IC50) unless `allow_mixed=True`, in which case markers differ by type.

### 5.9 Neighbourhood Card (v0.2)
- **Inputs:** query signature; annotated `ReferenceLibrary`; optional SMILES.
- **Statistics:** similarity (cosine on shrunk signatures) + null distribution from all library pairs → percentile and q value.
- **Encoding:** query structure (RDKit) + top-k neighbours with similarity bars, CIs, MoA colour.
- **Guards:** if no neighbour passes q < 0.05, the card says "no confident neighbour" in place of a ranked list.

### 5.10 Engagement Ladder (v0.2) — priority plot, the human particularly likes this one
- **Question:** does potency hold up as we move from purified protein to live cells to phenotype, and where does it break?
- **Inputs:** `PotencyTable` with rungs, in order:
  1. **Biochemical** (Kd or enzymatic IC50)
  2. **Cellular target engagement** (e.g. NanoBRET IC50 in live cells)
  3. **Pathway** (e.g. phospho-substrate IC50), optional
  4. **Phenotype** (viability IC50/GI50, or other phenotypic EC50)
  5. **Exposure** (free Cmax), optional
- **Statistics:** potency shift between consecutive rungs (Δlog). Optional expected ATP shift for ATP-competitive kinase inhibitors using Cheng–Prusoff (needs kinase ATP Km and assumed cellular ATP of ~1–5 mM; show as a grey wedge, not a line).
- **Encoding:** x = log10 molar concentration (shared), y = rungs top to bottom; one polyline per compound (slope-chart style); CI whiskers per rung. Small multiples per target.
- **Automatic flags (the point of the plot):**
  - **Phenotype more potent than target engagement** (by > 10×) → "phenotype likely driven by something other than {target}" (off-target or polypharmacology). Terracotta flag.
  - **Large biochemical → cellular shift beyond the ATP wedge** → "permeability, efflux or binding issue". Purple flag.
  - **Exposure below the phenotype rung** → "unlikely to reach active concentration in vivo".
- **Guards:** rungs measured in different cell systems are joined by **dashed** segments, and the cell system is labelled on each rung; missing rungs are skipped, never interpolated.
- **Caption:** "{compound}: {shift_bc:.1f} log shift biochemical→cellular ({expected} expected from ATP); phenotype is {rel} target engagement."

### 5.11 Combination Surface and Combination Lift (v0.2) — rethought
*Why rethought:* isobolograms and synergy heatmaps conflate **synergy** with **benefit**. A combination can be synergistic only at concentrations where both drugs barely work, or only above clinically reachable exposure. The plot must answer: *is there a dose pair, reachable in patients, where the combination clearly beats the best single agent?*

Two views, because public data come in two designs:

**(a) Combination Surface** — one cell line, full matrix (e.g. 7×7).
- **Statistics:** expected response under Bliss (default) and HSA; excess = observed − expected with replicate CIs. Mask cells where expected inhibition > 90% (ceiling: excess is meaningless there). Loewe optional, flagged as model-dependent.
- **Encoding:** log–log concentration axes; filled contours of **observed** viability interpolated only **within** the measured grid (no extrapolation), measured grid points as dots; outlined region where excess over Bliss is significant; a rectangle marking the clinically reachable box (Cmax range of each drug) when supplied; single-agent curves in the margins.
- **Guards:** model name printed; masked ceiling region hatched; if the matrix is anchored (not full), refuse and point to Lift.

**(b) Combination Lift** — many cell lines, works with anchored (2×7) designs.
- **Statistics:** per cell line, best single-agent Emax vs combination Emax (and ΔIC50 shift where available).
- **Encoding:** x = best single agent effect, y = combination effect, diagonal = no benefit; points = cell lines coloured by tissue or a biomarker; benefit threshold band (e.g. ΔEmax ≥ 20%, as used in the GDSC combination paper).
- **Guards:** threshold and model printed; cell lines with failed fits shown hollow.

### 5.12 Onset Timeline (v0.2) — rethought "Time Horizon"
*Why rethought:* horizon charts take training to read, and most drug datasets have 2–4 time points, so a dense-time-series chart is the wrong default. The useful question is *what happens first, and what follows?* Early, direct target effects vs later secondary stress or death programmes.
- **Inputs:** `TimeCourse` with ≥ 3 time points (2 allowed with a warning).
- **Statistics:** per pathway: onset = first time the score is significant and |score| > threshold; peak; resolution (return within threshold). With sparse time points, onset is **interval-censored**: between the last non-significant and first significant time point.
- **Encoding:** Gantt-style rows (pathways ordered by onset) on a time axis; onset drawn as an interval bracket, peak as a dot, sustained effect as a bar; colour by direction.
- **Option:** `style="horizon"` for dense time series (≥ 6 points, e.g. live-cell imaging).
- **Guards:** never draw onset as a single precise point when the data cannot support it.

### 5.13 Backlog specs (v0.3, do not build yet)
- **SAR Arrows:** matched molecular pairs (`mmpdb`) as arrows in (potency, toxicity) space; arrows only when |Δ| > 2× assay SD; filter to cliffs plus user-selected pairs.
- **Triage Flow:** custom ribbon renderer (not `matplotlib.sankey`); counts on every ribbon; "multiple reasons" bucket.
- **Cell-line Atlas:** response across a cell panel split by genotype or expression biomarker (natural fit for PRISM/GDSC).
- **Figure Audit:** score an existing figure spec for missing reference, dose, or uncertainty.

---

## 6. Real-data validation plan

Each case study is a script in `case_studies/` that downloads (or reads cached) data, builds plots, and writes `case_studies/out/<name>/report.md` containing the figures, captions, and a **ground-truth check**: a known biological result the plot must make visible. If the plot does not show it, that is a failed test of the plot, not a reason to tune the data.

### 6.1 Datasets

| ID | Dataset | Use for | Key facts | Status |
|---|---|---|---|---|
| D1 | **sci-Plex** (Srivatsan et al., *Science* 2020), GEO **GSE139944** | Compass, Dose Paths, Fingerprint, Effect Ladder, Retrievability | 188 compounds, 3 cancer cell lines (A549, MCF7, K562), multiple doses, ~650k single-cell transcriptomes | verified (dataset, scale). Dose levels and file layout: verify |
| D2 | **PRISM Repurposing secondary screen** (Corsello et al., *Nat Cancer* 2020), figshare 10.6084/m9.figshare.20564034 / depmap.org/repurposing | Window (viability side), Engagement Ladder (phenotype rung), Cell-line Atlas | 1,448 compounds × 489 cell lines, 8-point 4-fold dilution from 10 µM; prefer `screen_id` MTS010 where present; curve params include ec50, ic50, auc, r2 | verified |
| D3 | **GDSC monotherapy raw data**, cancerrxgene.org/downloads | Plate Atlas, Window | Raw CSV with plate barcode, well positions, controls; fit with the GDSC non-linear mixed-effects approach or our 4PL | verified (format doc). Pick a release and record it |
| D4 | **GDSC combinations, anchored** (Jaaks, Coker et al., *Nature* 2022) | Combination Lift | 2,025 combinations, 125 breast/colon/pancreas lines; **2 × 7 anchored design** (anchor at 2 concentrations, library 7-point); raw data on figshare; licence = internal research/education, no resale | verified |
| D5 | **GDSC MatrixExplorer, full matrix** (Bashi et al., *Cancer Discov* 2024; portal gdsc-combinations.depmap.sanger.ac.uk) | **Combination Surface** | 109 combinations, ~755 pan-cancer lines, **7 × 7 concentration matrix** | verified (design, scale). Download path and licence: verify |
| D6 | GDSC pan-cancer anchored (Vis et al., *Cell Rep Med* 2024) | Combination Lift (second dataset) | 51 combinations, 757 cell lines, anchored | verified (scale). Download: verify |
| D7 | **NanoBRET live-cell kinase engagement** (Vasta et al., *Cell Chem Biol* 2018) | Engagement Ladder (cellular rung), Selectivity Spectrum (cellular) | clinically relevant kinase inhibitors profiled against 178 full-length kinases in live cells; supplementary tables | verified (study). Table format: verify |
| D8 | **Davis et al. 2011** kinase Kd panel (*Nat Biotechnol*) | Engagement Ladder (biochemical rung), Selectivity Spectrum | 72 inhibitors × ~440 kinases, Kd | verify |
| D9 | **LINCS L1000 / CMap**, GEO GSE92742 and GSE70138 (or clue.io) | ReferenceLibrary for Compass and Neighbourhood | large reference signatures with MoA annotations; use Drug Repurposing Hub annotations | verify |
| D10 | **MIX-seq** (McFarland et al., *Nat Commun* 2020) | Onset Timeline | single-cell responses across pooled cell lines, with cell hashing for **time points** and doses | verified (study design). Which drugs have time courses and where the data live: verify |
| D11 | **JUMP Cell Painting** (cpg0016 in the Cell Painting Gallery) | Retrievability on morphology, Compass on morphology | very large compound set; replicate structure suited to mAP | verify; use a small subset only |
| D12 | ChEMBL | SAR (backlog), Window counter-screens | | verify, backlog |

### 6.2 Case studies and ground-truth checks

| Case | Data | Plots | Ground-truth check (must be visible) |
|---|---|---|---|
| CS1 sci-Plex HDAC | D1 | Retrievability → Compass → Dose Paths → Fingerprint → Effect Ladder | HDAC inhibitors form a coherent group (the sci-Plex paper highlights them); their Fingerprint rows match each other; dose increases magnitude along a consistent direction |
| CS2 PRISM windows | D2 | Window (viability as tox side) + Cell-line Atlas preview | Pick a targeted drug with a known sensitive lineage; sensitive lines show a window that insensitive lines lack |
| CS3 Plate QC | D3 | Plate Atlas | Z′ distribution across plates; edge-vs-interior ratio reported; any failing plates flagged |
| CS4 GDSC anchored combos | D4 | Combination Lift | Irinotecan + CHEK1 inhibitor in colon lines shows benefit, enriched in microsatellite-stable or KRAS–TP53 double-mutant lines (reported in Jaaks et al. 2022) |
| CS5 GDSC full-matrix combos | D5 | Combination Surface | For a combination reported as beneficial in Bashi et al. 2024, the significant excess region is inside the measured grid; for a non-beneficial one, no region appears |
| CS6 Engagement Ladder | D7 + D8 + D2 (or GDSC IC50) | Engagement Ladder | Crizotinib shows improved cellular selectivity vs biochemical (reported by Vasta et al. 2018); at least one compound where phenotype is more potent than primary-target engagement gets the off-target flag |
| CS7 Reference library | D9 | Compass, Neighbourhood | Held-out known-MoA compounds retrieve their own MoA as top neighbour at a rate well above chance; report the rate |
| CS8 Time | D10 | Onset Timeline | For a targeted drug (e.g. a MEK inhibitor if available), direct pathway output (e.g. MAPK targets) has an earlier onset than stress/apoptosis programmes |

Compound ID mapping across datasets (CS6 especially) is a known pain point: map via InChIKey where possible, then curated name synonyms; write unmatched names to `docs/data_status.md`.

---

## 7. Testing strategy

1. **Statistics unit tests** (synthetic fixtures, known answers):
   - 4PL fit recovers known parameters within tolerance; reports non-convergence correctly.
   - Bliss and HSA excess match hand-computed values on a 3×3 example; ceiling mask applied.
   - NNLS projection: query built as 0.7·refA + 0.3·refB recovers weights; residual ≈ 0 for pure mixtures and ≈ 1 for an orthogonal signature.
   - mAP retrieval matches `copairs` (or a hand-computed example) on a toy set.
   - Z′ matches the formula on synthetic plates; edge effect detected when injected.
   - Onset interval-censoring: with time points [2, 6, 24] and significance first at 6, onset = (2, 6].
   - Bend test: straight synthetic path → p > 0.05; injected bend → p < 0.05.
2. **Guard tests:** every red-team guard in section 5 has a test (e.g. Window with a failed fit draws no band and the caption contains "fit failed").
3. **Degenerate inputs:** no dose column, single replicate, all-NaN feature, one compound, 10,000 compounds. Each either works or raises a clear `TinctureDataError` with a fix suggestion.
4. **Visual regression:** `pytest-mpl` baselines for every plot on synthetic fixtures (deterministic seeds). Regenerate baselines only intentionally.
5. **Accessibility:** colour-blind simulation check on every baseline; minimum font size 6 pt at journal export size.
6. **Performance budget:** Fingerprint Strip with 1,000 compounds × 50 pathways renders in < 5 s; Compass with 500 compounds in < 3 s.
7. **Case studies** (section 6) run as a separate, opt-in `pytest -m realdata` suite; they are skipped (not failed) when data are absent.

---

## 8. Milestones

| # | Deliverable | Acceptance criteria |
|---|---|---|
| M0 | Repo scaffold, CI (lint, type check, tests), theme, PlotSpec + mpl renderer, export presets, caption/methods framework | `pip install -e .[omics]` works; a dummy plot exports at `nature_double` size; CI green |
| M1 | Data model + loaders for D1, D2, D3; `docs/data_status.md` with verified facts | Each loader round-trips into `PerturbationSet`/`DoseResponse`/`PlateTable`; unmatched or surprising contents documented |
| M2 | Signatures, shrinkage, pathway scores, Retrievability, Plate Atlas | Stats tests pass; CS3 report generated; Retrievability gate works on D1 |
| M3 | Compass, Dose Paths, Fingerprint Strip, Effect Ladder (+ classic fallbacks) | Guard tests pass; CS1 report shows the HDAC ground truth |
| M4 | v0.1 polish: gallery page, README with before/after figures, API docs | A new user can reproduce CS1 from the README in one command |
| M5 | Window, Selectivity Spectrum, Neighbourhood Card; loaders for D7, D8, D9 | CS2 and CS7 reports; Neighbourhood "no confident neighbour" path tested |
| M6 | **Engagement Ladder** | CS6 report with both flags demonstrated on real compounds |
| M7 | Combination Surface + Lift (D4, D5, D6), Onset Timeline (D10) | CS4, CS5, CS8 reports; ceiling mask and anchored-design refusal tested |

After each milestone, write a short `docs/milestone_<n>.md`: what was built, what was verified, what failed, open questions.

---

## 9. Decisions still open (ask the human, don't guess)

1. **Package name:** "Tincture" is a working name; check PyPI availability before publishing.
2. **Relationship to pertpy/scverse:** standalone package that reads AnnData (current plan) vs a pertpy extension.
3. **Licence** for the package (MIT vs Apache-2.0), and whether any case-study outputs can be published given dataset licences (GDSC restricts commercial reuse).
4. **Default synergy model:** Bliss (current plan) vs HSA as the headline.
5. **Stress axis for Dose Paths:** curated gene set vs viability readout when both exist.
6. Whether to depend on `copairs` or re-implement mAP.

---

## 10. Red-team summary (keep in mind throughout)

- New plots will be resisted by reviewers: always export the classic fallback and a methods paragraph.
- Reference-anchored views are only as good as the reference set: always print it, always show the residual.
- Pretty can mislead: every plot must show uncertainty and must say when a result is weak.
- The value is the statistics + data model, not the styling. If time is short, cut plots, not tests or guards.
- Success metric beyond "it renders": a small user study (8–10 scientists) measuring time-to-correct-answer and error rate, classic vs Tincture, on the case-study figures. Prepare the figure pairs in M4 so the study can run.
