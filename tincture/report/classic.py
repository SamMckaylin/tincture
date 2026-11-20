"""Registry of classic fallback plots.

Reviewers resist new plot types, so every tincture plot that has a conventional
counterpart (PCA, volcano, replicate scatter ...) can export it alongside. A plot
module registers its fallback builder here; :class:`tincture.result.PlotResult`
exposes it as ``result.classic``.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

from tincture.spec import PlotSpec

ClassicBuilder = Callable[..., PlotSpec]

_REGISTRY: dict[str, tuple[str, ClassicBuilder]] = {}


def register_classic(plot: str, name: str) -> Callable[[ClassicBuilder], ClassicBuilder]:
    """Decorator: register ``builder`` as the classic fallback for ``plot``.

    ``name`` is the conventional plot name (e.g. "PCA", "volcano").
    """

    def deco(builder: ClassicBuilder) -> ClassicBuilder:
        _REGISTRY[plot] = (name, builder)
        return builder

    return deco


def classic_name(plot: str) -> str | None:
    entry = _REGISTRY.get(plot)
    return entry[0] if entry else None


def build_classic(plot: str, **inputs: Any) -> PlotSpec | None:
    """Build the fallback spec for ``plot`` from the same inputs, or ``None`` if none exists."""
    entry = _REGISTRY.get(plot)
    if entry is None:
        return None
    spec = entry[1](**inputs)
    spec.meta.setdefault("classic_of", plot)
    return spec


def registered() -> dict[str, str]:
    """Mapping of plot name → classic fallback name."""
    return {k: v[0] for k, v in _REGISTRY.items()}
