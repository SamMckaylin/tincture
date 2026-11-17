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
