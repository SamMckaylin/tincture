"""Build docs/gallery.md: every plot next to its classic counterpart.

    python case_studies/cs1_sciplex_hdac.py   # once (real sci-Plex3 figures)
    python docs/make_gallery.py

Real-data figures come from CS1 (sci-Plex3, GEO GSE139944). Plate Atlas is shown on the
**synthetic** test plates: GDSC-derived plate figures stay local (licence, see
docs/data_status.md). Images are downscaled copies in docs/img/gallery/.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
CS1 = ROOT / "case_studies" / "out" / "cs1_sciplex_hdac"
IMG = ROOT / "docs" / "img" / "gallery"
OUT = ROOT / "docs" / "gallery.md"
WIDTH = 1000

#: (stem, plot, question, classic name, what the question-first version adds)
PLOTS = [
    (
        "1_retrievability",
        "Retrievability",
        "Is the effect real?",
        "replicate scatter",
        [
            "every condition is tested against vehicle, with an FDR-controlled call per "
            "condition, instead of one correlation driven by a few large genes",
            "the null band shows what 'no phenotype' looks like; non-retrievable "
            "conditions are greyed and later plots can drop them (the gate)",
        ],
    ),
    (
        "2_compass",
        "Compass",
        "What known mechanism does this compound resemble, and how strongly?",
        "PCA (reference classes coloured)",
        [
            "axes are mechanisms (reference classes), not principal components",
            "distance = noise-corrected signal strength; the residual says how much is "
            "unexplained, and compounds the references cannot explain are flagged",
            "the reference set and its version are printed; members are fitted leave-one-out",
        ],
    ),
    (
        "3_dose_paths",
        "Dose Paths",
        "Does the mechanism change with dose?",
        "PCA coloured by dose",
        [
            "each compound is a path from vehicle through its doses, on axes with meaning "
            "(likeness to the reference; stress)",
            "bends are tested (turn > 30° beyond noise) and only significant ones are "
            "annotated; the reference class's own path is drawn for comparison",
        ],
    ),
    (
        "4_fingerprint",
        "Fingerprint Strip",
        "What did each compound do, at pathway level?",
        "volcano (one compound)",
        [
            "many compounds × curated pathway groups at once, references first",
            "non-significant cells are neutral (BH over the whole matrix); colour "
            "limits and row ordering are stated",
        ],
    ),
    (
        "5_effect_ladder",
        "Effect Ladder",
        "Which genes moved, and how sure are we?",
        "volcano",
        [
            "the top genes with 95% intervals, grouped by gene set; intervals crossing "
            "zero are grey",
            "the count strip keeps the volcano's overview: how many genes went up and "
            "down across the transcriptome",
        ],
    ),
]


OUT_CS = ROOT / "case_studies" / "out"

#: Further real-data plots from other case studies: (folder, stem, plot, question, classic,
#: data note, what changes). Skipped with a note when the case study has not been run.
EXTRA = [
    (
        "cs6_engagement_ladder",
        "engagement_ladder",
        "Engagement Ladder",
        "Does potency hold up from enzyme to cell to phenotype, and where does it break?",
        "grouped bar chart of potencies",
        "Biochemical Kd (Davis et al. 2011), live-cell engagement in HEK293 (Vasta et al. "
        "2018; apparent IC50 from single-dose occupancy) and PRISM viability (CS6). The "
        "dasatinib flag is discussed in the CS6 report: BCR-ABL in K562 vs full-length ABL1 "
        "in HEK293 is the likely reason.",
        [
            "every rung on one concentration axis, so the shift from protein to cell to "
            "phenotype is read directly",
            "automatic flags: phenotype far more potent than engagement (off-target?), shifts "
            "beyond the expected ATP effect, exposure below the active concentration",
            "dashed segments and labels wherever rungs come from different cell systems",
        ],
    ),
    (
        "cs6_engagement_ladder",
        "selectivity_crizotinib_nanobret",
        "Selectivity Spectrum",
        "How specific is this compound?",
        "bar chart of most potent targets",
        "Crizotinib engagement of 178 kinases in live HEK293 cells (Vasta et al. 2018, "
        "apparent IC50 from occupancy at 1 µM); the biochemical spectrum (Davis 2011) is in "
        "the CS6 report: 69 kinases bound biochemically, 13 engaged in cells.",
        [
            "the whole panel, strongest first, with censored (not detected) targets shown "
            "and counted",
            "off-targets judged against exposure (free Cmax band) when supplied, otherwise by "
            "a stated 10-fold rule; panel coverage in the footer",
        ],
    ),
    (
        "cs2_prism_windows",
        "window_ACH-000219",
        "Window",
        "Is there a dose that works without the toxicity?",
        "overlaid dose-response curves",
        "PRISM 20Q2 secondary (Corsello et al. 2020, CC BY 4.0): dabrafenib in A375 (BRAF "
        "V600E melanoma) against the median of 676 non-skin lines as the toxicity proxy "
        "(CS2).",
        [
            "EC90 (on-target) and EC10 (toxicity) with bootstrap intervals instead of IC50s",
            "the window is drawn only when both fits pass quality checks, and is open-ended "
            "when toxicity never appears in the tested range",
            "exposure (free Cmax) is a band, never a line, and only when supplied",
        ],
    ),
    (
        "cs7_reference_library",
        "card_bortezomib_unrelated",
        "Neighbourhood Card",
        "Which annotated compounds does this one resemble?",
        "ranked similarity bars",
        "LINCS L1000 phase II (GEO GSE70138) consensus signatures, 603 compounds with "
        "Repurposing Hub MoA; bortezomib left out of its own library (CS7). Shown with the "
        "**unrelated-pairs null**; under the handover's all-pairs null (the current default) "
        "the same card reports 'no confident neighbour' (top neighbour oprozomib, q = 0.056). "
        "Which null to use is an open decision (docs/milestone_5.md).",
        [
            "each neighbour gets a q value against the library's own null, and only "
            "confident neighbours are coloured; with none, the card says so",
            "similarity intervals from the query's own replicate signatures",
            "the query's structure sits next to its neighbours",
        ],
    ),
]


def _caption(md: Path) -> str:
    if not md.exists():
        return ""
    m = re.search(r"## Caption\n\n(.+?)\n\n", md.read_text(encoding="utf-8"), re.S)
    return m.group(1).strip() if m else ""


def _shrink(src: Path, dst: Path) -> None:
    im = Image.open(src)
    if im.width > WIDTH:
        im = im.resize((WIDTH, round(im.height * WIDTH / im.width)), Image.Resampling.LANCZOS)
    dst.parent.mkdir(parents=True, exist_ok=True)
    im.convert("RGB").save(dst, optimize=True)


def _plate_atlas_synthetic(tmp: Path) -> tuple[Path, Path, str]:
    sys.path.insert(0, str(ROOT))
    import tincture as tc
    from tests.fixtures.synthetic_plates import synthetic_plates

    res = tc.plate_atlas(
        synthetic_plates(n_plates=40, bad_plates=(7, 23), edge_factor=0.75),
        title="Plate Atlas (synthetic test plates)",
        footer=("Synthetic plates from tests/fixtures; not real data.",),
    )
    written = tc.export(res, tmp / "plate_atlas", formats=("png",))
    return written["png"], written["classic_png"], res.caption


def main() -> int:
    if not (CS1 / "2_compass.png").exists():
        print("Run python case_studies/cs1_sciplex_hdac.py first.", file=sys.stderr)
        return 1
    lines = [
        "# Gallery",
        "",
        "Each plot answers one drug-discovery question. The left column is the plot people "
        "usually draw (tincture exports it too, as the *classic fallback*); the right column "
        "is the question-first version. Real data: sci-Plex3 A549, 24 h "
        "(Srivatsan et al. 2020, GEO GSE139944); abexinostat, an HDAC inhibitor, is the "
        "focus compound (terracotta). Reproduce with `python case_studies/cs1_sciplex_hdac.py` "
        "then `python docs/make_gallery.py`.",
        "",
    ]
    for stem, plot, question, classic, adds in PLOTS:
        _shrink(CS1 / f"{stem}.png", IMG / f"{stem}.png")
        _shrink(CS1 / f"{stem}_classic.png", IMG / f"{stem}_classic.png")
        lines += [
            f"## {plot}: {question}",
            "",
            f"| Classic: {classic} | tincture: {plot} |",
            "|---|---|",
            f"| ![classic](img/gallery/{stem}_classic.png) | ![{plot}](img/gallery/{stem}.png) |",
            "",
            f"**Caption (generated).** {_caption(CS1 / f'{stem}.md')}",
            "",
            "**What changes:**",
            *[f"- {a}" for a in adds],
            "",
        ]
    for folder, stem, plot, question, classic, note, adds in EXTRA:
        src = OUT_CS / folder
        if not (src / f"{stem}.png").exists():
            lines += [f"## {plot}: {question}", "", f"_Run the {folder} case study first._", ""]
            continue
        _shrink(src / f"{stem}.png", IMG / f"{stem}.png")
        _shrink(src / f"{stem}_classic.png", IMG / f"{stem}_classic.png")
        lines += [
            f"## {plot}: {question}",
            "",
            note,
            "",
            f"| Classic: {classic} | tincture: {plot} |",
            "|---|---|",
            f"| ![classic](img/gallery/{stem}_classic.png) | ![{plot}](img/gallery/{stem}.png) |",
            "",
            f"**Caption (generated).** {_caption(src / f'{stem}.md')}",
            "",
            "**What changes:**",
            *[f"- {a}" for a in adds],
            "",
        ]
    tmp = ROOT / "case_studies" / "out" / "gallery_tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    pa, pa_classic, pa_cap = _plate_atlas_synthetic(tmp)
    _shrink(pa, IMG / "plate_atlas_synthetic.png")
    _shrink(pa_classic, IMG / "plate_atlas_synthetic_classic.png")
    lines += [
        "## Plate Atlas: can we trust this plate?",
        "",
        "Shown on **synthetic** test plates: GDSC-derived plate figures are kept local "
        "(licence; see `docs/data_status.md`). Run `case_studies/cs3_plate_qc.py` to see it "
        "on all 16,354 GDSC2 plates.",
        "",
        "| Classic: plate heatmap | tincture: Plate Atlas |",
        "|---|---|",
        "| ![classic](img/gallery/plate_atlas_synthetic_classic.png) "
        "| ![Plate Atlas](img/gallery/plate_atlas_synthetic.png) |",
        "",
        f"**Caption (generated).** {pa_cap}",
        "",
        "**What changes:**",
        "- Z′ per plate with the 0.5 line, edge effect and control drift next to the plate "
        "map, so a bad plate is visible at a glance",
        "- controls and absent wells are marked; the QC definitions are printed",
        "",
    ]
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {OUT} and {len(list(IMG.glob('*.png')))} images")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
