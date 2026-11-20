"""User-study kit: classic vs tincture figure pairs, questions, answer key, scoring.

HANDOVER §10: a small study (8–10 scientists) measuring time-to-correct-answer and error
rate, classic vs tincture, on the case-study figures. Protocol: docs/user_study.md.

    python case_studies/cs1_sciplex_hdac.py         # once: figures + results.json
    python case_studies/user_study.py build         # -> case_studies/out/user_study/
    python case_studies/user_study.py score responses.csv

Answers come from the same computed results as the figures (CS1 results.json), never typed
in by hand. Each item says whether the classic figure can answer it *in principle*; the
primary analysis uses only items answerable from both versions.
"""

from __future__ import annotations

import csv
import json
import shutil
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
CS1 = ROOT / "out" / "cs1_sciplex_hdac"
GALLERY = ROOT.parent / "docs" / "img" / "gallery"
OUT = ROOT / "out" / "user_study"
N_PARTICIPANTS = 10
SEED = 20261117


@dataclass
class Item:
    id: str
    figure: str  # stem in the CS1 output folder (or gallery for the synthetic plate)
    question: str
    answer: str
    accept: list[str]  # lower-case strings any of which marks an answer correct
    both_answerable: bool  # can the classic figure answer it in principle?
    note: str = ""


def items(f: dict) -> list[Item]:
    top = f["fingerprint_top_pathways"][0]
    up, down = f["ladder_counts"]["up"], f["ladder_counts"]["down"]
    return [
        Item(
            "Q01",
            "1_retrievability",
            "Does abexinostat at 10 µM have a reproducible effect (yes / no)?",
            "yes" if f["retrievable_by_dose"]["10 µM"] else "no",
            ["yes"] if f["retrievable_by_dose"]["10 µM"] else ["no"],
            True,
        ),
        Item(
            "Q02",
            "2_compass",
            "Which known mechanism class does abexinostat most resemble?",
            str(f["compass_top_reference"]),
            [str(f["compass_top_reference"]).lower(), "histone deacetylase", "hdac inhibitor"],
            True,
            "Classic PCA colours the reference classes.",
        ),
        Item(
            "Q03",
            "2_compass",
            "Is abexinostat's signature mostly explained by the known mechanisms, or does "
            "a large part remain unexplained?",
            "mostly explained" if f["compass_unexplained"] < 0.5 else "largely unexplained",
            ["explained", "mostly explained"]
            if f["compass_unexplained"] < 0.5
            else ["unexplained", "largely unexplained"],
            False,
            "PCA has no residual; included as a secondary item.",
        ),
        Item(
            "Q04",
            "3_dose_paths",
            "Does abexinostat's response change direction as the dose increases (yes / no)?",
            "yes" if f["abexinostat_bends"] else "no",
            ["yes"] if f["abexinostat_bends"] else ["no"],
            True,
            "Classic PCA marks each compound and labels its top dose.",
        ),
        Item(
            "Q05",
            "3_dose_paths",
            "Does panobinostat change direction with dose? If so, above which dose?",
            f"yes, above {f['panobinostat_bend_above']}" if f["panobinostat_bend_above"] else "no",
            [str(f["panobinostat_bend_above"]).lower().replace("µ", "u"), "1 um", "1um"]
            if f["panobinostat_bend_above"]
            else ["no"],
            True,
        ),
        Item(
            "Q06",
            "4_fingerprint",
            "Which biological process is most strongly changed by abexinostat, and in which "
            "direction?",
            f"{top['pathway']} ({top['direction']})",
            [top["pathway"].lower(), "cell cycle", "proliferation"],
            False,
            "Classic is a volcano: needs gene-level knowledge. Score process and direction.",
        ),
        Item(
            "Q07",
            "4_fingerprint",
            "Which other compound has the most similar pathway profile to abexinostat?",
            str(f.get("fingerprint_nearest_name", f["fingerprint_nearest"])),
            [str(f.get("fingerprint_nearest_name", f["fingerprint_nearest"])).lower()],
            False,
            "Classic shows one compound only.",
        ),
        Item(
            "Q08",
            "5_effect_ladder",
            "Which single gene decreases the most?",
            f["ladder_top_gene"],
            [f["ladder_top_gene"].lower()],
            True,
        ),
        Item(
            "Q09",
            "5_effect_ladder",
            "Are more genes significantly decreased or increased?",
            "decreased" if down > up else "increased",
            ["decreased", "down"] if down > up else ["increased", "up"],
            True,
            f"{down:,} down vs {up:,} up (q < 0.05).",
        ),
    ]


def is_correct(item: Item, answer: str) -> bool:
    a = answer.strip().lower().replace("µ", "u")
    return any(k.replace("µ", "u") in a for k in item.accept)


def build() -> Path:
    facts_p = CS1 / "results.json"
    if not facts_p.exists():
        raise SystemExit("Run python case_studies/cs1_sciplex_hdac.py first.")
    facts = json.loads(facts_p.read_text(encoding="utf-8"))
    its = items(facts)
    if OUT.exists():
        shutil.rmtree(OUT)
    (OUT / "figures").mkdir(parents=True)
    for stem in sorted({i.figure for i in its}):
        for version, suffix in (("tincture", ""), ("classic", "_classic")):
            src = CS1 / f"{stem}{suffix}.png"
            shutil.copy(src, OUT / "figures" / f"{stem}_{version}.png")
    # answer key (for the experimenter only)
    (OUT / "answer_key.json").write_text(
        json.dumps([asdict(i) for i in its], indent=2, ensure_ascii=False), encoding="utf-8"
    )
    # counterbalanced assignment: each participant sees every item once; half the items in
    # each version, alternating between participants; item order shuffled per participant
    rng = np.random.default_rng(SEED)
    rows = []
    for p in range(1, N_PARTICIPANTS + 1):
        order = rng.permutation(len(its))
        for pos, k in enumerate(order):
            it = its[k]
            version = "tincture" if (k + p) % 2 == 0 else "classic"
            rows.append(
                {
                    "participant": f"P{p:02d}",
                    "position": pos + 1,
                    "item": it.id,
                    "version": version,
                    "figure": f"figures/{it.figure}_{version}.png",
                    "question": it.question,
                }
            )
    with open(OUT / "assignments.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(OUT / "responses_template.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["participant", "item", "version", "answer", "seconds"])
        for r in rows:
            w.writerow([r["participant"], r["item"], r["version"], "", ""])
    print(f"wrote {OUT}: {len(its)} items, {N_PARTICIPANTS} participants, {len(rows)} trials")
    return OUT


def score(responses: Path) -> None:
    its = {i["id"]: Item(**i) for i in json.loads((OUT / "answer_key.json").read_text())}
    data = []
    with open(responses, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if not r["answer"].strip():
                continue
            it = its[r["item"]]
            data.append(
                (
                    r["participant"],
                    it.id,
                    r["version"],
                    is_correct(it, r["answer"]),
                    float(r["seconds"]),
                    it.both_answerable,
                )
            )
    if not data:
        raise SystemExit("No answered trials.")
    print("Primary analysis: items answerable from both versions")
    for v in ("classic", "tincture"):
        d = [x for x in data if x[2] == v and x[5]]
        if d:
            err = 1 - np.mean([x[3] for x in d])
            t_ok = [x[4] for x in d if x[3]]
            med = np.median(t_ok) if t_ok else float("nan")
            print(
                f"  {v:8s} n={len(d):3d}  error rate {err:.0%}  median time-to-correct {med:.0f} s"
            )
    print("Secondary: all items")
    for v in ("classic", "tincture"):
        d = [x for x in data if x[2] == v]
        if d:
            print(f"  {v:8s} n={len(d):3d}  error rate {1 - np.mean([x[3] for x in d]):.0%}")


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "build":
        build()
    elif len(sys.argv) >= 3 and sys.argv[1] == "score":
        score(Path(sys.argv[2]))
    else:
        print(__doc__)
