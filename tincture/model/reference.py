"""ReferenceLibrary: named reference signatures with mechanism annotations.

A reference is the mean signature of its member compounds (one signature per member,
normally at its top dose). Member signatures are kept so that a query compound that is
itself a member can be projected against references rebuilt *without* it
(:meth:`ReferenceLibrary.without`): a compound must never be explained by a reference that
contains it. The library's name, version and citation are printed on every figure that
uses it (HANDOVER §4.2).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

from tincture.errors import TinctureDataError


@dataclass
class ReferenceLibrary:
    names: list[str]
    features: np.ndarray
    member_ids: dict[str, list[str]]
    member_signatures: dict[str, np.ndarray]  # name -> (n_members, n_features)
    version: str
    citation: str = ""
    source: str = ""
    meta: dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.features = np.asarray(self.features).astype(str)
        if not self.names:
            raise TinctureDataError(
                "Reference library is empty.",
                fix="Provide at least two reference classes with ≥ 2 usable members each.",
            )
        for n in self.names:
            m = np.asarray(self.member_signatures[n], dtype=float)
            if m.ndim != 2 or m.shape[1] != self.features.size:
                raise TinctureDataError(f"Reference {n!r}: member matrix has the wrong shape.")
            if m.shape[0] != len(self.member_ids[n]):
                raise TinctureDataError(f"Reference {n!r}: ids and signatures differ in number.")
            self.member_signatures[n] = m

    @property
    def matrix(self) -> np.ndarray:
        """References × features (member means)."""
        return np.vstack([self.member_signatures[n].mean(axis=0) for n in self.names])

    def __len__(self) -> int:
        return len(self.names)

    def label(self) -> str:
        return f"{self.source or 'reference library'} {self.version}".strip()

    def similarity(self) -> np.ndarray:
        """Cosine similarity between references."""
        M = self.matrix
        U = M / np.linalg.norm(M, axis=1, keepdims=True)
        return U @ U.T

    def reference_of(self, compound_id: str) -> list[str]:
        return [n for n in self.names if compound_id in self.member_ids[n]]

    def without(self, compound_id: str) -> ReferenceLibrary:
        """The library with ``compound_id`` removed from every reference (leave-one-out).

        A reference left with no members is dropped (and recorded in ``meta``).
        """
        names, ids, sigs, dropped = [], {}, {}, []
        for n in self.names:
            keep = [i for i, c in enumerate(self.member_ids[n]) if c != compound_id]
            if not keep:
                dropped.append(n)
                continue
            names.append(n)
            ids[n] = [self.member_ids[n][i] for i in keep]
            sigs[n] = self.member_signatures[n][keep]
        out = ReferenceLibrary(
            names=names,
            features=self.features,
            member_ids=ids,
            member_signatures=sigs,
            version=self.version,
            citation=self.citation,
            source=self.source,
            meta={**self.meta, "left_out": compound_id, "dropped": dropped},
        )
        return out

    def reorder(self, order: Sequence[int]) -> ReferenceLibrary:
        names = [self.names[i] for i in order]
        return ReferenceLibrary(
            names=names,
            features=self.features,
            member_ids={n: self.member_ids[n] for n in names},
            member_signatures={n: self.member_signatures[n] for n in names},
            version=self.version,
            citation=self.citation,
            source=self.source,
            meta=dict(self.meta),
        )


def library_from_signatures(
    signatures: Mapping[str, np.ndarray],
    features: np.ndarray,
    classes: Mapping[str, Sequence[str]],
    version: str,
    source: str = "",
    citation: str = "",
    min_members: int = 2,
) -> ReferenceLibrary:
    """Build a library from per-compound signatures and class membership.

    ``signatures`` maps compound id → effect vector (aligned to ``features``); classes with
    fewer than ``min_members`` available members are skipped (recorded in ``meta``).
    """
    names, ids, sigs, skipped = [], {}, {}, {}
    for name, members in classes.items():
        have = [c for c in members if c in signatures]
        if len(have) < min_members:
            skipped[name] = len(have)
            continue
        names.append(name)
        ids[name] = list(have)
        sigs[name] = np.vstack([np.asarray(signatures[c], dtype=float) for c in have])
    if len(names) < 2:
        raise TinctureDataError(
            f"Only {len(names)} reference class(es) have ≥ {min_members} usable members.",
            fix="Lower min_members, add compounds, or relax the retrievability gate.",
        )
    return ReferenceLibrary(
        names=names,
        features=np.asarray(features),
        member_ids=ids,
        member_signatures=sigs,
        version=version,
        citation=citation,
        source=source,
        meta={"skipped": skipped},
    )
