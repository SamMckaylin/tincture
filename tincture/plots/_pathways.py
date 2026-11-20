"""Curated pathway order and labels shared by the Fingerprint Strip and the Effect Ladder.

HANDOVER §5.3: columns follow a **fixed curated order**, grouped Cell cycle / Stress /
Metabolism / Signalling / Immune / Other. The grouping below covers all 50 MSigDB Hallmark
sets; any other set goes to "Other" in alphabetical order.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd

GROUPS: dict[str, list[str]] = {
    "Cell cycle": [
        "E2F_TARGETS",
        "G2M_CHECKPOINT",
        "MITOTIC_SPINDLE",
        "MYC_TARGETS_V1",
        "MYC_TARGETS_V2",
        "DNA_REPAIR",
    ],
    "Stress": [
        "P53_PATHWAY",
        "APOPTOSIS",
        "UNFOLDED_PROTEIN_RESPONSE",
        "UV_RESPONSE_UP",
        "UV_RESPONSE_DN",
        "HYPOXIA",
        "REACTIVE_OXYGEN_SPECIES_PATHWAY",
    ],
    "Metabolism": [
        "OXIDATIVE_PHOSPHORYLATION",
        "GLYCOLYSIS",
        "FATTY_ACID_METABOLISM",
        "CHOLESTEROL_HOMEOSTASIS",
        "ADIPOGENESIS",
        "PEROXISOME",
        "HEME_METABOLISM",
        "XENOBIOTIC_METABOLISM",
        "BILE_ACID_METABOLISM",
    ],
    "Signalling": [
        "MTORC1_SIGNALING",
        "PI3K_AKT_MTOR_SIGNALING",
        "KRAS_SIGNALING_UP",
        "KRAS_SIGNALING_DN",
        "TNFA_SIGNALING_VIA_NFKB",
        "IL6_JAK_STAT3_SIGNALING",
        "IL2_STAT5_SIGNALING",
        "TGF_BETA_SIGNALING",
        "WNT_BETA_CATENIN_SIGNALING",
        "NOTCH_SIGNALING",
        "HEDGEHOG_SIGNALING",
        "ESTROGEN_RESPONSE_EARLY",
        "ESTROGEN_RESPONSE_LATE",
        "ANDROGEN_RESPONSE",
    ],
    "Immune": [
        "INTERFERON_ALPHA_RESPONSE",
        "INTERFERON_GAMMA_RESPONSE",
        "INFLAMMATORY_RESPONSE",
        "ALLOGRAFT_REJECTION",
        "COMPLEMENT",
        "COAGULATION",
    ],
    "Other": [
        "EPITHELIAL_MESENCHYMAL_TRANSITION",
        "ANGIOGENESIS",
        "APICAL_JUNCTION",
        "APICAL_SURFACE",
        "PROTEIN_SECRETION",
        "MYOGENESIS",
        "PANCREAS_BETA_CELLS",
        "SPERMATOGENESIS",
    ],
}

_LABELS = {
    "E2F_TARGETS": "E2F targets",
    "G2M_CHECKPOINT": "G2/M checkpoint",
    "MYC_TARGETS_V1": "MYC targets v1",
    "MYC_TARGETS_V2": "MYC targets v2",
    "DNA_REPAIR": "DNA repair",
    "P53_PATHWAY": "p53 pathway",
    "UNFOLDED_PROTEIN_RESPONSE": "Unfolded protein resp.",
    "UV_RESPONSE_UP": "UV response up",
    "UV_RESPONSE_DN": "UV response down",
    "REACTIVE_OXYGEN_SPECIES_PATHWAY": "ROS pathway",
    "OXIDATIVE_PHOSPHORYLATION": "Oxidative phosphorylation",
    "MTORC1_SIGNALING": "mTORC1",
    "PI3K_AKT_MTOR_SIGNALING": "PI3K/AKT/mTOR",
    "KRAS_SIGNALING_UP": "KRAS up",
    "KRAS_SIGNALING_DN": "KRAS down",
    "TNFA_SIGNALING_VIA_NFKB": "TNFα via NF-κB",
    "IL6_JAK_STAT3_SIGNALING": "IL-6/JAK/STAT3",
    "IL2_STAT5_SIGNALING": "IL-2/STAT5",
    "TGF_BETA_SIGNALING": "TGF-β",
    "WNT_BETA_CATENIN_SIGNALING": "WNT/β-catenin",
    "NOTCH_SIGNALING": "Notch",
    "HEDGEHOG_SIGNALING": "Hedgehog",
    "INTERFERON_ALPHA_RESPONSE": "Interferon α",
    "INTERFERON_GAMMA_RESPONSE": "Interferon γ",
    "EPITHELIAL_MESENCHYMAL_TRANSITION": "EMT",
    "PANCREAS_BETA_CELLS": "Pancreas β cells",
}

_RANK = {s: (g, i) for g, (_, sets) in enumerate(GROUPS.items()) for i, s in enumerate(sets)}
GROUP_NAMES = list(GROUPS)


def bare(name: str) -> str:
    return str(name).removeprefix("HALLMARK_")


def pathway_label(name: str) -> str:
    b = bare(name)
    if b in _LABELS:
        return _LABELS[b]
    return b.replace("_", " ").capitalize()


def group_of(name: str) -> str:
    b = bare(name)
    return GROUP_NAMES[_RANK[b][0]] if b in _RANK else "Other"


def curated_order(names: Iterable[str]) -> list[str]:
    """Sort pathway names into the fixed curated order (unknown sets last, alphabetical)."""
    other = len(GROUPS) - 1

    def key(n: str) -> tuple[int, int, str]:
        b = bare(n)
        g, i = _RANK.get(b, (other, 10_000))
        return g, i, b

    return sorted(names, key=key)


def group_breaks(ordered: list[str]) -> list[int]:
    """Indices after which a new group starts (for heatmap column gaps)."""
    gs = [group_of(n) for n in ordered]
    return [i for i in range(len(gs) - 1) if gs[i] != gs[i + 1]]


def feature_sets(net: pd.DataFrame, features: np.ndarray) -> dict[str, str]:
    """Assign each feature to one gene set: the first that contains it in curated order."""
    order = curated_order(net["source"].astype(str).unique())
    rank = {s: k for k, s in enumerate(order)}
    sub = net[net["target"].astype(str).isin(set(np.asarray(features).astype(str)))].copy()
    sub["rank"] = sub["source"].astype(str).map(rank)
    first = sub.sort_values("rank").drop_duplicates("target")
    return dict(zip(first["target"].astype(str), first["source"].map(bare), strict=True))
