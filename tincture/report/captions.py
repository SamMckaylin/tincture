"""Templated reading-guide captions (no LLM).

A caption is built from sentences. Each sentence is a ``str.format`` template filled
from computed statistics. If a statistic is unavailable (``None`` or NaN) the
sentence is *dropped* rather than printed with a placeholder; a template that
references a field that was never supplied is a bug and raises.

Weak-result notes (principle 8) are added with :meth:`Caption.weak` and always appear,
prefixed with "Note:".
"""

from __future__ import annotations

import math
import string
from dataclasses import dataclass, field
from typing import Any

from tincture.errors import TinctureSpecError

#: Caption templates from the brief (HANDOVER section 5), kept in one place.
TEMPLATES: dict[str, str] = {
    "compass": (
        "{compound} most resembles {top_ref} (weight {w:.2f}); {pct_unexplained:.0%} of its "
        "signature is unexplained by the {n_refs} references."
    ),
    "dose_paths": (
        "{compound} tracks {ref} up to {dose_ok}; above {dose_bend} it moves toward "
        "{stress_axis} (bend p = {p:.2g})."
    ),
    "fingerprint": (
        "{n} compounds × {m} pathways. {focus} is most similar to {nearest_row}; strongest "
        "effect: {pathway} ({direction})."
    ),
    "engagement_ladder": (
        "{compound}: {shift_bc:.1f} log shift biochemical→cellular ({expected} expected from "
        "ATP); phenotype is {rel} target engagement."
    ),
}


def _missing(v: Any) -> bool:
    """True for None and for any real-number NaN (Python or numpy scalar)."""
    if v is None:
        return True
    if isinstance(v, (str, bytes, bool)):
        return False
    try:
        return math.isnan(float(v))
    except (TypeError, ValueError):
        return False


def template_fields(template: str) -> list[str]:
    """Top-level field names referenced by a format template."""
    names = []
    for _, fname, _, _ in string.Formatter().parse(template):
        if fname:
            names.append(fname.split(".")[0].split("[")[0])
    return names


def fmt_p(p: float | None, digits: int = 2) -> str:
    """Format a p/q value: ``p < 0.001`` style below 1e-3, else 2 significant digits."""
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return "n/a"
    if p < 1e-3:
        return "< 0.001"
    return f"= {p:.{digits}g}"


def fmt_conc(molar: float) -> str:
    """Format a molar concentration with an SI prefix (e.g. 3.2e-7 → '320 nM')."""
    if molar is None or not math.isfinite(molar) or molar <= 0:
        return "n/a"
    # Round to 3 significant figures first so 999.9 nM becomes "1 µM", not "1e+03 nM".
    molar = float(f"{molar:.3g}")
    for unit, scale in (("M", 1.0), ("mM", 1e-3), ("µM", 1e-6), ("nM", 1e-9), ("pM", 1e-12)):
        if molar >= scale * (1 - 1e-12):
            return f"{molar / scale:.3g} {unit}"
    return f"{molar / 1e-15:.3g} fM"


@dataclass
class Caption:
    """Builder for a reading-guide caption."""

    sentences: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    dropped: list[str] = field(default_factory=list)

    def add(self, template: str, /, **fields: Any) -> Caption:
        """Add a sentence; dropped (and recorded) if any referenced statistic is missing."""
        names = template_fields(template)
        absent = [n for n in names if n not in fields]
        if absent:
            raise TinctureSpecError(f"Caption template needs fields {absent}: {template!r}")
        if any(_missing(fields[n]) for n in names):
            self.dropped.append(template)
            return self
        self.sentences.append(template.format(**fields))
        return self

    def add_if(self, condition: bool, template: str, /, **fields: Any) -> Caption:
        return self.add(template, **fields) if condition else self

    def weak(self, text: str) -> Caption:
        """Record that a result is weak or the plot was reduced (always shown)."""
        if text not in self.notes:
            self.notes.append(text)
        return self

    def render(self) -> str:
        parts = list(self.sentences)
        parts += [f"Note: {n}" if not n.startswith("Note:") else n for n in self.notes]
        return " ".join(p.strip() for p in parts if p.strip())

    def __str__(self) -> str:
        return self.render()
