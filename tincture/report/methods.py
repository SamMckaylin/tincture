"""Methods paragraphs.

Plots assemble a methods paragraph from reusable *blocks* (one per statistical step),
each with an optional citation. The paragraph ends with software versions so a figure
can be reproduced. Blocks for specific statistics are registered by the stats modules
as they are built; the generic ones live here.
"""

from __future__ import annotations

import importlib.metadata as md
from dataclasses import dataclass, field
from typing import Any

from tincture.errors import TinctureSpecError
from tincture.report.captions import template_fields


@dataclass(frozen=True)
class MethodBlock:
    key: str
    text: str
    citations: tuple[str, ...] = ()


BLOCKS: dict[str, MethodBlock] = {}


def register_block(key: str, text: str, citations: tuple[str, ...] = ()) -> MethodBlock:
    """Register (or replace) a reusable methods block."""
    b = MethodBlock(key, text, citations)
    BLOCKS[key] = b
    return b


register_block(
    "bootstrap_ci",
    "Uncertainty is shown as {level:.0%} percentile bootstrap confidence intervals "
    "({n_boot} resamples over {unit}).",
    ("Efron B, Tibshirani RJ. An Introduction to the Bootstrap. Chapman & Hall; 1993.",),
)
register_block(
    "bh_fdr",
    "P values were adjusted for multiple testing with the Benjamini–Hochberg procedure "
    "(q < {alpha}).",
    (
        "Benjamini Y, Hochberg Y. Controlling the false discovery rate. "
        "J R Stat Soc B. 1995;57:289–300.",
    ),
)
register_block(
    "colour",
    "Colours follow the tincture palette: terracotta marks the focus compound, blue "
    "references, purple unexplained signal; non-significant values are grey.",
)


def software_versions(
    packages: tuple[str, ...] = ("tincture", "numpy", "scipy", "matplotlib"),
) -> str:
    parts = []
    for pkg in packages:
        try:
            parts.append(f"{pkg} {md.version(pkg)}")
        except md.PackageNotFoundError:
            if pkg == "tincture":
                from tincture._version import __version__

                parts.append(f"tincture {__version__}")
    return ", ".join(parts)


@dataclass
class Methods:
    """Builder for a methods paragraph."""

    sentences: list[str] = field(default_factory=list)
    citations: list[str] = field(default_factory=list)
    _keys: set[str] = field(default_factory=set)

    def use(self, key: str, /, **fields: Any) -> Methods:
        """Add a registered block (once per key), filled with ``fields``."""
        if key not in BLOCKS:
            raise TinctureSpecError(f"Unknown methods block {key!r}.")
        if key in self._keys:
            return self
        block = BLOCKS[key]
        absent = [n for n in template_fields(block.text) if n not in fields]
        if absent:
            raise TinctureSpecError(f"Methods block {key!r} needs fields {absent}.")
        self.sentences.append(block.text.format(**fields))
        for c in block.citations:
            if c not in self.citations:
                self.citations.append(c)
        self._keys.add(key)
        return self

    def text(self, sentence: str) -> Methods:
        """Add a free-text sentence specific to this plot."""
        self.sentences.append(sentence)
        return self

    def render(self, versions: bool = True) -> str:
        body = " ".join(s.strip() for s in self.sentences if s.strip())
        if versions:
            body = f"{body} Software: {software_versions()}.".strip()
        if self.citations:
            refs = " ".join(f"[{i}] {c}" for i, c in enumerate(self.citations, 1))
            body = f"{body}\nReferences: {refs}"
        return body
