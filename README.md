# tincture

Question-first, publication-quality plots for drug-discovery data, built on proper statistics.

Drug-discovery papers recycle the same generic plots (PCA, UMAP, volcano, GSEA bars). Those show the *structure* of the data, not answers to drug-discovery questions. tincture organises plotting around the questions:

| Question | Plot | Status |
|---|---|---|
| What does this compound look like? | Compass (**available**), Neighbourhood Card (M5) | **available** / planned |
| Does the mechanism change with dose? | Dose Paths | **available** |
| What did it do biologically? | Fingerprint Strip, Effect Ladder | **available** |
| Is the effect real? | Retrievability, Plate Atlas | **available** |
| Is there a safe, effective dose? | Window | planned (M5) |
| How specific is it? | Selectivity Spectrum | planned (M5) |
| Does potency hold up from enzyme to cell to phenotype? | Engagement Ladder | planned (M6) |
| Do two drugs work better together, at reachable doses? | Combination Surface / Lift | planned (M7) |
| What happens first, and what follows? | Onset Timeline | planned (M7) |

Every plot ships with:
- the figure;
- a **reading-guide caption** templated from the computed statistics (no LLM);
- a **methods paragraph** with citations and software versions;
- where a conventional plot exists, a **classic fallback** (e.g. PCA, volcano) for reviewers.

Uncertainty is on by default. When a result is weak (a failed fit, no confident neighbour, a bend within noise), the figure says so.

> **Status: Milestone 3 of 7.** Done so far: the rendering core, the data model, verified loaders for sci-Plex, PRISM and GDSC, the statistics layer (signatures with empirical-Bayes moderation and shrinkage, pathway scores, NNLS projection, dose-path bend test) and six question-first plots: **Retrievability**, **Plate Atlas**, **Compass**, **Dose Paths**, **Fingerprint Strip** and **Effect Ladder**. See [`HANDOVER.md`](HANDOVER.md) §8 and the reports in [`docs/`](docs/).

## Install

```bash
pip install -e ".[dev]"
```

Optional extras: `[omics]` (pydeseq2, decoupler, scanpy; about 600 MB installed), `[chem]` (rdkit, mmpdb). Requires Python ≥ 3.10.

## Quick look

```python
import tincture as tc
from tincture.io.sciplex import load_sciplex3
from tincture.stats.signatures import replicate_profiles

# Is each compound's phenotype real? (the gate before every other plot)
profiles = replicate_profiles(load_sciplex3(cell_lines=["A549"]))
res = tc.retrievability(profiles, focus=["S1090|1e-06"])  # abexinostat 1 µM
print(res.caption)  # templated reading guide, says when a result is weak
tc.export(res, "figures/retrievability", preset="nature_double")
# writes .pdf/.png, a .md sidecar (caption, methods, warnings), and the classic replicate scatter
```

```python
from tincture.io.gdsc import load_gdsc2_plates
from tincture.stats.zprime import Selector

# Can we trust these plates? Z′, edge effect and control drift
plates = load_gdsc2_plates(max_plates=300)
qc = tc.plate_atlas(
    plates,
    neg=Selector(tags=("NC-1",)),
    pos=Selector(well_types=("blank",)),
    edge=Selector(tags=("UN-USED", "NC-1", "NC-0")),
)
```

```python
from tincture.datasets.genesets import load_geneset
from tincture.stats.paths import stress_gene_set

sigs = load_sciplex3(cell_lines=["A549"]).signatures()[0]  # shrunk logFC vs vehicle
hallmark = load_geneset("hallmark").table
# lib = library_from_signatures(...)  # reference classes, e.g. sci-Plex pathway annotation

tc.compass({c: sigs.get(c, 1e-5) for c in compounds}, lib, focus="S1090")  # which mechanism?
tc.dose_paths(sigs, ["S1090", "S1030"], lib, stress_gene_set(hallmark), focus="S1090")
tc.fingerprint(sigs, hallmark, focus="S1090", references=hdac_ids)  # pathway-level strip
tc.effect_ladder(sigs, "S1090", 1e-5, gene_sets=hallmark)  # which genes, how sure
```

The full worked example is `case_studies/cs1_sciplex_hdac.py` (HDAC inhibitors in sci-Plex3 A549; report in `case_studies/out/`, summary in [`docs/milestone_3.md`](docs/milestone_3.md)).

| | |
|---|---|
| ![Compass](docs/img/cs1_compass.png) | ![Dose Paths](docs/img/cs1_dose_paths.png) |
| ![Fingerprint Strip](docs/img/cs1_fingerprint.png) | ![Effect Ladder](docs/img/cs1_effect_ladder.png) |

*Real data: sci-Plex3 (Srivatsan et al. 2020, GSE139944), A549, 24 h. Abexinostat in terracotta.*

Every plot returns a `PlotResult` with `.caption`, `.methods`, `.warnings`, `.classic` and `.save()`.

![demo plot on synthetic data](docs/img/demo_synthetic.png)

*The M0 demo plot on **synthetic** test data. It exercises the pipeline only and is not a tincture analysis plot.*

## Design in one paragraph

Plot builders compute statistics and return a `PlotSpec`, a plain dataclass of panels, marks, scales, caption, methods and footer. `render_mpl(spec, preset=...)` is the only code that touches matplotlib, which keeps an interactive backend possible later. Colours are semantic: terracotta always means the compound the figure is about, blue is reference or known-MoA, purple is novel or unexplained, and green is good or window. Sequential colours are perceptually uniform maps only; jet and rainbow raise an error.

## Export rules

`tc.export` / `PlotResult.save`:
- renders at the exact journal width (`nature_single` 89 mm, `nature_double` 183 mm, `cell_*`, `acs_*`; see [`docs/design/journal_presets.md`](docs/design/journal_presets.md));
- embeds TrueType fonts;
- **refuses** to export a figure with text under 6 pt, or with categorical colours that become indistinguishable under simulated deuteranopia or protanopia (CIEDE2000 < 10).

Output does not depend on your matplotlib style.

## Data

Case studies run on real public datasets. **No data is stored in this repository.** Loaders download into a git-ignored cache (`$TINCTURE_DATA`, or `data/raw/` in a checkout) and verify every file against a pinned checksum, because several sources (GDSC in particular) forbid redistribution.

| Dataset | Loader | Returns |
|---|---|---|
| sci-Plex3 (GSE139944) | `tincture.io.sciplex.load_sciplex3()` | `PerturbationSet` of per-well pseudobulk counts (run `build_pseudobulk()` once; about 5 min) |
| PRISM Repurposing 20Q2 secondary | `tincture.io.prism.load_prism_secondary(compounds=[...])` | `DoseResponseSet` (raw viability + published fits with quality flags) |
| GDSC2 raw, release 8.5 | `tincture.io.gdsc.load_gdsc2_plates(max_plates=...)` | `PlateTable` (1536-well) |

```python
from tincture.io.prism import load_prism_secondary

prism = load_prism_secondary(compounds=["dabrafenib"])
curve = prism.get("BRD-K09951645", "ACH-000824")  # a DoseResponse
curve.quality_flags()  # why a fit should not be trusted, if at all
```

What was verified about each source, and where it differs from the papers or READMEs, is in [`docs/data_status.md`](docs/data_status.md). Examples: PRISM's README swaps the asymptotes, and GDSC wells are double-annotated.

## Development

```bash
pytest
```

This runs the unit, guard and visual-regression tests.

```bash
pytest -m realdata
```

This runs the case studies on downloaded data; they are skipped when data are absent.

```bash
pytest -k baseline --mpl-generate-path=tests/baseline
```

This regenerates image baselines; do it only intentionally.

```bash
ruff check . && ruff format --check . && mypy
```

Image baselines are tied to the matplotlib version in `tests/baseline/MATPLOTLIB_VERSION`. With any other version, pixel comparison is skipped.

## Licence

Apache-2.0. Dataset licences are separate and are noted per dataset in `docs/data_status.md`.
