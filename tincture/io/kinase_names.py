"""Curated kinase name synonyms between panels (HANDOVER §6.2: map by curated synonyms
when no shared identifier exists; unmatched names are listed in docs/data_status.md).

Keys are NanoBRET (Vasta 2018) names; values are the Davis 2011 ``Entrez Gene Symbol``.
Each pair was checked against both files: e.g. MLTK = ZAK (= MAP3K20), SIK3 = QSK (Davis
gene symbol KIAA0999), SNF1LK2 = SIK2.
"""

from __future__ import annotations

NANOBRET_TO_DAVIS = {
    "JNK1": "MAPK8",
    "JNK3": "MAPK10",
    "MLTK": "ZAK",
    "SNF1LK2": "SIK2",
    "SIK3": "KIAA0999",
}

#: NanoBRET kinases with no counterpart in the Davis 2011 panel.
NANOBRET_NOT_IN_DAVIS = ("PKN3",)


def to_davis(name: str) -> str:
    return NANOBRET_TO_DAVIS.get(name, name)
