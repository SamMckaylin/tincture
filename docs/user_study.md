# User study: classic vs tincture figures (protocol)

HANDOVER §10: "Success metric beyond 'it renders': a small user study (8–10 scientists) measuring time-to-correct-answer and error rate, classic vs tincture, on the case-study figures." M4 prepares everything the study needs. Running it needs people, so that part is yours.

## Materials

```bash
python case_studies/cs1_sciplex_hdac.py   # figures + results.json (answer source)
python case_studies/user_study.py build   # -> case_studies/out/user_study/
```

`case_studies/out/user_study/` contains:

| File | For | Content |
|---|---|---|
| `figures/*_classic.png`, `*_tincture.png` | participants | 5 figure pairs from CS1 (sci-Plex3 A549, 24 h) |
| `assignments.csv` | experimenter | per participant: item order, version shown, figure, question |
| `responses_template.csv` | experimenter | one row per trial: fill in `answer` and `seconds` |
| `answer_key.json` | experimenter only | correct answers, accepted wordings, answerability flag |

The classic figures are the best conventional versions, so the comparison is fair:
- the PCA used for Compass colours the reference classes and labels the focus compound;
- the PCA by dose marks each compound and labels its top dose;
- the volcano labels the most significant genes;
- the replicate scatter shows the focus condition.

## Items

Answers are computed from the same results as the figures (`results.json`); none are typed in by hand.

| Item | Figure | Question | Classic can answer in principle? |
|---|---|---|---|
| Q01 | Retrievability | Does abexinostat at 10 µM have a reproducible effect? | yes |
| Q02 | Compass | Which known mechanism class does abexinostat most resemble? | yes |
| Q03 | Compass | Is its signature mostly explained by known mechanisms? | no (PCA has no residual) |
| Q04 | Dose Paths | Does abexinostat's response change direction with dose? | yes |
| Q05 | Dose Paths | Does panobinostat change direction, and above which dose? | yes |
| Q06 | Fingerprint | Which process is most strongly changed, and in which direction? | no (volcano needs gene knowledge) |
| Q07 | Fingerprint | Which compound has the most similar pathway profile? | no (volcano shows one compound) |
| Q08 | Effect Ladder | Which single gene decreases the most? | yes |
| Q09 | Effect Ladder | Are more genes significantly decreased or increased? | yes |

## Design

- **Participants:** 8–10 scientists who read omics or screening figures in their work; none involved in building tincture.
- **Within-subject, counterbalanced:**
  - each participant answers every item once: about half the items on the classic figure, half on the tincture figure;
  - which version an item gets alternates between participants, so each item is seen equally often in both versions (5 + 5 with 10 participants);
  - item order is shuffled per participant (fixed seed, `assignments.csv`).
- **Procedure:**
  - show one figure with its question, and nothing else: no caption, no methods;
  - start the clock when the figure appears and stop it when the participant commits to an answer;
  - no time limit, but record "don't know" as an incorrect answer;
  - a short practice item on a figure not in the study (e.g. the synthetic Plate Atlas in `docs/img/gallery/`) comes first.
- **Captions off:** tincture figures ship with generated captions that state many answers outright, so they are withheld to compare the figures themselves. A follow-up arm can test figures with captions.

## Outcomes and analysis (fixed before data collection)

- **Primary:** error rate and median time-to-correct-answer, classic vs tincture, on the items the classic figure can answer in principle (Q01, Q02, Q04, Q05, Q08, Q09).
- **Secondary:** the same on all items. Items the classic cannot answer measure what the new plots add, not how easy they are to read, so report them separately.
- **Statistics:**
  - paired by participant;
  - error rate: exact sign / McNemar-type comparison per item pair;
  - time: Wilcoxon signed-rank on per-participant median times;
  - with n = 8–10, report effect sizes and intervals, not just p-values.
- **Scoring:** `python case_studies/user_study.py score responses.csv` gives a first pass by matching accepted wordings. An experimenter should review every automatically scored answer, because keyword matching can accept an answer that mentions the right word in the wrong sense.

## Known limitations

- The study uses one dataset (sci-Plex3) and one focus compound. Q01–Q09 are not independent: several are about abexinostat.
- Participants will recognise the tincture style after a few items, so blinding to version isn't possible. The randomised order spreads the learning effect.
- Plate Atlas isn't in the study: the real GDSC figures stay local (licence) and the synthetic example isn't real data.
