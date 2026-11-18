# Milestone 7: Combination Surface, Combination Lift, Onset Timeline

## Your decisions applied

- **Synergy models: Bliss (default) and ZIP**, as you asked. HSA is also computed: the handover's tests require it, and the Surface caption reports benefit over the best single agent next to synergy whatever the headline model. Loewe is not offered (the GDSC matrix documentation states there is no fitted surface to support it).
- **ATP Km:** searched; see the "Decisions for you" section below.

## Data (all verified, details in `docs/data_status.md`)

| ID | Source | What was checked |
|---|---|---|
| D4 | Jaaks 2022, figshare (fitted data, original + validation, glossary) | 296,707 measurements, 125 lines, 2,025 combinations, 2 × 7 anchored; synergy definition read from the paper. **Licence conflict:** figshare says CC BY 4.0, but the descriptions restrict use to internal research, so it is treated as restricted. |
| D5 | Bashi 2024, figshare project 163378 (per-well 7 × 7 results, 241 MB zip) | 7 × 7, 1,000-fold range; per-well semantics from the glossary; GDSC's Bliss formula checked numerically |
| D6 | Vis 2024, open-access supplement via Europe PMC | 51 anchor–library pairs (matches); **762 cell-line names vs 757 in the paper** (not resolved) |
| D10 | MIX-seq, figshare (CC BY 4.0) | Time course exists for **trametinib only** (experiment 5: 3, 6, 12, 24 and 48 h, 24 lines). Supplementary Table 3 also lists 9 h, but no 9 h hash exists in the data. |
| — | Cell Model Passports model list | MSI status for the colon lines: 30 MSS / 15 MSI (the paper reports 31 MSS) |

Network: you approved `gdsc-combinations.depmap.sanger.ac.uk` and `www.ebi.ac.uk`. The GDSC Combinations portal builds its download links in JavaScript, and the permission policy refused queries to its internal endpoints. I stopped there and found the same data on figshare, as named in the paper's Data Availability statement.

## What was built

- **`ComboMatrix`**: the design (full / anchored / partial) is inferred from the concentration grid.
- **`stats.synergy`** (Bliss, ZIP, HSA):
  - Bliss and HSA are hand-checked on a 3 × 3 example (§7). ZIP follows Yadav et al. 2015: row and column curves are refitted with the partner drug's effect as the baseline, and δ is the averaged fit minus the Bliss expectation from fitted single agents.
  - **Ceiling mask:** cells with expected inhibition > 90% are hatched and never tested.
  - **Uncertainty is honest about replicates.** With ≥ 2 replicate matrices, the excess is computed within each replicate and tested with a t-test. With one replicate and a stated noise SD, a parametric bootstrap is used. With neither, there is no significance claim and the caption says so.
  - **The first version was overconfident.** It bootstrapped over 3 replicates, which gave false positives on a null surface. I replaced it before writing the tests, and the tests now include a null surface.
  - The significant region requires BH q < 0.05 **and** excess ≥ 10%, so a trivial but "significant" excess is not outlined.
- **`tc.combination_surface`:**
  - left panel: observed viability as contours, interpolated only inside the measured grid, with dots at the measured points; the significant region outlined in terracotta; ceiling cells hatched; a dashed free-Cmax box when Cmax values are supplied; single-agent curves in the margins;
  - right panel: the excess per cell, with values printed.
  - Anchored designs are refused with a pointer to the Lift.
- **`tc.combination_lift`:**
  - the best single agent against the combination, with the no-benefit diagonal and a benefit band (default ΔEmax ≥ 20% over the best single agent);
  - groups compared with Fisher and Mann–Whitney tests; one facet per anchor concentration; failed fits hollow and not counted.
- **`stats.onset` + `tc.onset_timeline`:**
  - onset is interval-censored (the handover test case [2, 6, 24] with first significance at 6 gives (2, 6]); peak and resolution are reported;
  - Gantt rows are ordered by onset; never-active pathways are listed in the footer;
  - `style="horizon"` is allowed only with ≥ 6 time points.
- **Loaders:** `io.gdsc_combo` (Jaaks, Bashi, Vis, Cell Model Passports) and `io.mixseq` (a streaming pseudobulk of the 1.2 GB matrix in 19 s, cached).
- **Tests:** 351 passing (32 new for M7), plus 13 opt-in real-data tests (4 new). New baselines for Surface, Lift and Onset, all colour-blind-checked.
- **API reference:** it had silently left out every module added since M5. `make_api.py` now lists them all (1,097 lines).

## Acceptance

### CS4: Combination Lift on Jaaks 2022 (AZD7762 + camptothecin, colon). PASS

| anchor | benefit, all lines | MSS | MSI | Fisher p |
|---|---|---|---|---|
| AZD7762 0.0625 µM | 20/44 (45%) | 17/29 (59%) | 3/15 (20%) | 0.025 |
| AZD7762 0.25 µM | 24/43 (56%) | 18/28 (64%) | 6/15 (40%) | 0.20 |

The paper's own Bliss-based synergy calls show the same pattern (MSS 52–54% vs MSI 13–27%).

Second dataset (D6, report only): trametinib + pictilisib shows benefit in 17% / 40% of 750 lines across the two anchors, **with no lung enrichment** (Fisher p = 0.17 and 0.85). The paper's lung result is about ECB, a different metric that also requires efficacy above 50%, so this does not contradict it, but it is not a replication either.

### CS5: Combination Surface on Bashi 2024. FAIL on check 1

| Pre-specified check | Result |
|---|---|
| Selumetinib + venetoclax (beneficial): Bliss region in ≥ 4 of the 6 AML lines the paper names | ❌ **2 of 6** (EoL-1-cell, OCI-AML5); HSA region also 2 of 6 |
| Olaparib + AZD5153 (no benefit): no region in ≥ 80% of AML lines | ✅ 9 of 10 |
| Agreement with GDSC's own Bliss score (secondary) | Spearman ρ = 0.72 over 16 matrices |

**Why it failed, as far as I can tell:**
- Every matrix is a **single replicate**.
- The noise SD I fixed in advance (GDSC's own single-agent fit RMSE) is 0.13–0.23 for the four lines that failed.
- The paper calls a line active from point estimates (HSA > 0.1 and Emax > 0.5) without any significance test.
- So the plot is doing what it promises, refusing to outline an excess it cannot distinguish from noise. But the pre-specified check failed, and I have not changed the rule after the fact.
- A report-only diagnostic (added after the first run, labelled as such) applies the paper's point-estimate rule per cell (HSA excess ≥ 0.1 at viability ≤ 0.5, no test). It flags **6 of 6** beneficial lines but also **6 of 10** non-beneficial lines. Point estimates are sensitive but not specific; the significance rule is specific (9 of 10 negatives clean) but under-powered on single replicates.

This needs a decision from you (below).

### CS8: Onset Timeline on MIX-seq trametinib. PASS

- **MAPK output (PROGENy) goes down by 3 h**, together with EGFR, NFkB and TNFa.
- Cell-cycle programmes (E2F, G2/M, MYC) follow in 6–12 h.
- PROGENy p53 rises in 12–24 h.
- Hallmark apoptosis, p53 pathway and unfolded protein response never become active within 48 h.

The pre-specified check (MAPK onset earlier than every stress/apoptosis programme) holds.

## What failed or changed along the way

1. **First uncertainty model was overconfident**: see above. It was caught by the null test before any real-data run.
2. **CS5 failed its first check**: see above. The rule was not retuned.
3. **The CS5 run is slow** (about 30 minutes). The parametric bootstrap refits the single-agent curves 200 times per matrix, and ZIP also refits 14 row and column curves each time. It works, but making it faster (vectorised fits, or fewer draws for ZIP) is worth doing before anyone runs hundreds of matrices.
4. **Jaaks `ANCHOR_CONC` is text**, because a few anchors are two-drug composites ("0.125|0.003"). The loader parses it as a number and the composites are excluded.
5. **Missing modules in the API reference since M5**: fixed (see above).
6. **Licence:** GDSC-derived figures (CS3, CS4, CS5) stay in `case_studies/out/` (git-ignored). The gallery shows the Combination Surface and Lift on **synthetic** data, labelled. The Onset Timeline gallery image is MIX-seq (CC BY 4.0).

## Decisions for you

1. **CS5 single-replicate significance.** Options:
   - **(a)** Keep the fixed rule, so CS5 stands as FAIL. This is the honest record that single 7 × 7 matrices with fit RMSE ~0.2 cannot support a per-cell significance claim.
   - **(b)** Pre-register a different noise estimate for single replicates, for example the SD of the 49 combination wells' deviations from a smooth surface, or plate-level control CVs if I can find them in the data. Then rerun as a new, labelled check.
   - **(c)** For single replicates, outline cells on point estimates (as the paper does), drawn dashed and captioned "not tested". The diagnostic argues against this: it would also outline 6 of 10 non-beneficial lines.

   I'd lean towards (b), with (a) kept on record, but it's your call.
2. **ATP Km for the CS6 purple flag.** I found and verified values for ABL1 (9.8 µM), ALK (0.13–0.16 mM) and MET (60 µM phosphorylated, 200 µM unphosphorylated), but from **three different papers, constructs and assay formats**; vendor tables were unreachable. Applying them, with the wedge for 1–5 mM ATP:
   - dasatinib/ABL1: shift 3.0 vs wedge 2.0–2.7, so the flag would fire;
   - crizotinib/ALK: shift 2.0 vs wedge 0.8–1.6, so it would fire;
   - crizotinib/MET: shift 1.4 vs wedge 0.8–1.9, so it would not.

   Should I wire these in (labelled per source), or leave the purple flag demonstrated on synthetic data only?
3. **Still open from §9:** package name (a PyPI check), and the remaining §9 items (pertpy relationship, Dose Paths stress axis, copairs dependency) if you want to close them before a release.

## Follow-up after review (2026-11-18)

You chose **(b) for CS5** and asked me to **wire in the ATP Km values**.

### CS6: ATP Km applied

- Values: ABL1 9.8 µM (Grover 2015); ALK 134–159 µM (Bresler 2011, unphosphorylated and phosphorylated); MET 60–200 µM (Behshad 2010, phosphorylated and unphosphorylated).
- Where two enzyme forms are reported, `atp_shift` now takes a (low, high) Km range. The wedge then spans every combination of ATP and Km, so the flag fires only beyond the widest expectation (tested).
- **Purple flag:**
  - dasatinib/ABL1: shift 3.0 vs 2.0–2.7 expected, fires;
  - crizotinib/ALK: 2.0 vs 0.9–1.6, fires;
  - crizotinib/MET: 1.4 vs 0.8–1.9, does not fire.
- It is reported, not pass/fail, with the caveats: three labs and three constructs; the cellular rung is an apparent IC50 from occupancy; the biochemical rung is a domain-binding Kd. The footer names every Km source. CS6's pre-specified checks still PASS.

### CS5 version 2 (option b): invalid, because its calibration gate failed

- **Estimator (fixed before the run):** within-matrix noise from robust second differences along rows and columns of unsaturated wells.
- **Gate (also fixed before the run):** on replicated matrices, the estimate must be at least 0.8× the real replicate SD.
- **Result:** on **384 replicated matrices** (gdsc-008 and gdsc-010, ≥ 3 replicates), the median ratio is **0.71** (IQR 0.56–0.87). The estimate median is 0.064 viability, against a replicate SD median of 0.088.
- **Interpretation:** a single matrix only sees within-plate noise, but replicate matrices also differ plate to plate. Version 2 would therefore be anti-conservative, and its checks are invalid as pre-registered.
- **For the record only (invalid, not used):** under version 2 the Bliss region appears in 3 of 6 beneficial lines (HSA benefit in 5 of 6), and no region in 10 of 10 non-beneficial lines.
- **CS5 therefore stands as FAIL** under both versions. Version 1 failed its sensitivity check; version 2 failed its own validity gate.

### What I'd suggest next (your decision)

- **(d) Calibrated noise:** scale the within-matrix estimate by the factor measured on the replicated matrices (1 / 0.71 ≈ 1.4, ideally estimated per screen), pre-registered as version 3 with a fresh gate. One caveat: these 384 matrices were used to *derive* the factor, so the gate would need held-out matrices (e.g. calibrate on gdsc-008, check on gdsc-010).
- **Or accept the result:** the outcome is informative as it stands. Single-replicate 7×7 matrices can support a point-estimate *screen* but not a per-cell significance claim. The Surface already says so in its caption when no noise estimate is given.
