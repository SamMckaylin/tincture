# Milestone 6: Engagement Ladder

## What was built

**Data (D7, D8): the files you supplied, verified and loaded.** Both are pinned by SHA-256 as manual downloads (`RemoteFile(manual=True)`): if a file is missing, the error says where to download it in a browser and where to put it.
- `io.davis`: Davis 2011 Supplementary Table 4.
  - **442 assays × 72 inhibitors**, Kd in nM (the handover said ~440).
  - **Blank cells are censored** (> 10 µM), not missing.
  - Constructs are classified as mutant / wild type / non-human; by default the loader keeps the most potent wild-type construct per gene.
- `io.nanobret`: Vasta 2018 Table S1. **The file holds single-dose occupancy, not IC50s**: duplicates for 178 kinases at 1 µM crizotinib and 100 nM dasatinib, in HEK293 cells (doses verified from the article text).
  - The loader derives an **apparent IC50 = C·(1−f)/f** (1:1 binding, Hill 1), labelled as model-based.
  - "No Binding" is censored, with the detection limit taken as 5% occupancy.
- `io.kinase_names`: curated NanoBRET ↔ Davis synonyms (JNK1 = MAPK8, MLTK = ZAK, SIK3 = KIAA0999, …). Only PKN3 is unmatched (recorded).
- **Both loaders reproduce the paper:** crizotinib engages **13** kinases above 50% in cells against **69** predicted from Kd, exactly the paper's numbers. Dasatinib gives 23 vs 44 (paper: 22 vs 44).

**`tc.engagement_ladder`** (§5.10)
- Layout: rungs biochemical → cellular → pathway → phenotype → exposure on a shared log-concentration axis; one polyline per compound; interval whiskers; small multiples per target, stacked so each keeps a wide concentration axis.
- Statistics: the shift between consecutive rungs (Δlog). The **Cheng–Prusoff ATP wedge**, log10(1 + [ATP]/Km,ATP) for ATP 1–5 mM, is drawn grey only when a kinase ATP Km is supplied.
- Flags (each tested):
  - **phenotype > 10× more potent than cellular engagement** → "likely driven by something other than {target}" (terracotta);
  - **biochemical→cellular shift beyond the wedge** → "permeability, efflux or binding issue" (purple); not evaluated without an ATP Km, and the caption says so;
  - **exposure below the phenotype rung**.
- Guards (each tested):
  - **dashed** segments between different cell systems, and every point carries its cell system;
  - **missing rungs skipped**, never interpolated;
  - censored values are hollow and used only when the bound is conclusive. `PotencyTable` gained `<` censoring (e.g. already below 50% viability at the lowest dose).
- Caption follows the handover template; classic fallback is grouped bars.
- **Selectivity Spectrum on real data** (crizotinib, biochemical and cellular). Labels are now stacked in a column with leader lines; at real panel sizes they piled up.
- Spec: `PolygonMark`.

**Tests:** 319 passing, plus 9 opt-in real-data tests, including the paper reproduction (13 vs 69). New baselines for Engagement Ladder and Selectivity, colour-blind-checked.

## Acceptance (CS6)

`case_studies/out/cs6_engagement_ladder/report.md`; everything below was pre-specified in the script header.

| Check | Result |
|---|---|
| Crizotinib shows improved cellular selectivity vs biochemical (Vasta 2018) | ✅ 13 of 177 shared kinases engaged > 50% in cells at 1 µM vs 69 predicted from Kd; the paper reports 13 vs 69 |
| At least one compound gets the off-target flag | ✅ **dasatinib/ABL1**: K562 phenotype IC50 < 0.61 nM (below PRISM's lowest dose) vs apparent ABL1 engagement 29 nM in HEK293, so ≥ 48× more potent |
| ATP-shift (purple) flag on real data | ❌ **not evaluable.** None of the sources gives kinase ATP Km values, so the wedge cannot be drawn. The flag works and is tested on synthetic ladders. |

Other ladders, unflagged: crizotinib/ALK (2.0-log shift into cells; KARPAS-299 phenotype 4× more potent than engagement) and crizotinib/MET (1.4-log shift; MKN45 phenotype in line with engagement).

**How to read the dasatinib flag.** K562 is driven by BCR-ABL, so its phenotype should be on-target. The flag compares it with NanoBRET engagement of *full-length ABL1 in HEK293*: a different protein context and cell system. The ladder marks this with a dashed segment and labels on each rung, and the report says so. The flag behaves as specified, but here the likely explanation is the engagement system, not an off-target. I report it as a demonstrated flag with that caveat, not as a biological finding.

## What failed or changed along the way

1. **NanoBRET data is occupancy, not IC50.** The handover expects a NanoBRET IC50 rung, but the public table has single-dose occupancy only. The cellular rung is therefore model-based (1:1 binding, Hill 1), with replicate min–max as its interval, and labelled so everywhere.
2. **Phenotype below the tested range.** Dasatinib in K562 is already at 38% viability at 0.6 nM. This needed "more potent than" censoring, and the flags only use a censored value when its bound is conclusive.
3. **Floating-point slip in CS6.** `1000 nM × 1e-9 ≤ 1e-6` is false in floating point, so the first run reported 67 instead of 69. Fixed by comparing in nM.
4. **Layout.** Three side-by-side ladders squeezed the concentration axis and the flag text overflowed; targets are now stacked. The ⚑ glyph is not in Inter, so flags say "Flag:".
5. **Gallery.** The Selectivity Spectrum entry now uses real data (crizotinib in live cells) instead of the synthetic panel, and the Engagement Ladder is added. Both are figures derived from published supplementary tables, with citations; the raw tables stay local.

## Decision (from you)

- **Neighbourhood null:** the handover's all-pairs null stays the default (recorded in `docs/milestone_5.md`).

## Open questions

1. **Purple flag on real data:** it needs ATP Km values for the kinases (e.g. ABL1, ALK, MET). Do you have a source you trust, or do you want me to look for one? It would need a new allowed host and a check of each value. Otherwise CS6 stands with this flag demonstrated on synthetic data only.
2. **Exposure rung:** free Cmax for crizotinib and dasatinib would enable the exposure flag on real data. Same question: is there a source you want used?
3. **Still open from §9:** package name (a PyPI check), and **Bliss vs HSA as the default synergy model, needed for M7**.
