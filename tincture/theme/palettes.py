"""Colour palette and colormap rules (HANDOVER section 3, "Theme").

Colours are referenced *semantically* in PlotSpecs (``"focus"``, ``"reference"``,
``"novel"`` ...) and resolved here, so the meaning of a colour is fixed across the
whole package: terracotta always means "the thing this figure is about".
"""

from __future__ import annotations

import re
from typing import Final

from matplotlib.colors import BoundaryNorm, Colormap, ListedColormap, to_rgba

from tincture.errors import TinctureStyleError

# --- raw palette -------------------------------------------------------------------------

BLUE: Final = "#2F7CC4"
BLUE_LIGHT: Final = "#9CCDE8"
NEUTRAL: Final = "#D1D1D1"
PEACH: Final = "#F6BD78"
TERRACOTTA: Final = "#D97757"
PURPLE: Final = "#7C5CBF"
PURPLE_LIGHT: Final = "#CDB6EC"
GREEN: Final = "#2E8B57"
GREEN_LIGHT: Final = "#B7DFA4"

INK: Final = "#222222"  # text, axes
INK_SOFT: Final = "#6B6B6B"  # secondary text, gridlines labels
GRID: Final = "#E6E6E6"
PAPER: Final = "#FFFFFF"

#: Five-step diverging scale, strongest down -> strongest up.
DIVERGING_5: Final[tuple[str, ...]] = (BLUE, BLUE_LIGHT, NEUTRAL, PEACH, TERRACOTTA)

#: Semantic names usable anywhere a colour is expected in a PlotSpec.
SEMANTIC: Final[dict[str, str]] = {
    "focus": TERRACOTTA,
    "focus_light": PEACH,
    "reference": BLUE,
    "reference_light": BLUE_LIGHT,
    "known": BLUE,
    "novel": PURPLE,
    "novel_light": PURPLE_LIGHT,
    "good": GREEN,
    "window": GREEN_LIGHT,
    "good_light": GREEN_LIGHT,
    "neutral": NEUTRAL,
    "nonsig": NEUTRAL,
    "down": BLUE,
    "down_light": BLUE_LIGHT,
    "up": TERRACOTTA,
    "up_light": PEACH,
    "warn": TERRACOTTA,
    "ink": INK,
    "ink_soft": INK_SOFT,
    "grid": GRID,
    "paper": PAPER,
    "none": "none",
}

# Colours that carry categorical meaning (used by the colour-blind checker on specs).
CATEGORICAL_SEMANTIC: Final[frozenset[str]] = frozenset(
    {"focus", "reference", "known", "novel", "good", "down", "up", "warn"}
)

_HEX = re.compile(r"^#(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{6}|[0-9a-fA-F]{8})$")


def resolve_color(color: str) -> str:
    """Resolve a semantic colour name or hex string to a hex string.

    Raises
    ------
    TinctureStyleError
        If ``color`` is neither a semantic name nor a valid hex colour / matplotlib name.
    """
    if color in SEMANTIC:
        return SEMANTIC[color]
    if _HEX.match(color):
        return color
    try:
        to_rgba(color)
    except ValueError as exc:
        raise TinctureStyleError(
            f"Unknown colour {color!r}. Use a semantic name ({', '.join(sorted(SEMANTIC))}) "
            "or a hex string."
        ) from exc
    return color


# --- colormaps ---------------------------------------------------------------------------

#: Perceptually uniform sequential maps from matplotlib (the only ones allowed).
ALLOWED_SEQUENTIAL: Final[frozenset[str]] = frozenset(
    {"viridis", "cividis", "magma", "inferno", "plasma"}
)
#: Explicitly banned (non-uniform or rainbow) maps; listed so the error is explicit.
BANNED_CMAPS: Final[frozenset[str]] = frozenset(
    {
        "jet",
        "rainbow",
        "hsv",
        "gist_rainbow",
        "nipy_spectral",
        "gist_ncar",
        "spectral",
        "Spectral",
        "turbo",
        "gist_stern",
        "flag",
        "prism",
        "cubehelix_rainbow",
    }
)


def diverging_cmap(steps: int = 5) -> ListedColormap:
    """The tincture diverging map (blue = down, terracotta = up), discrete by default."""
    if steps != 5:
        from matplotlib.colors import LinearSegmentedColormap

        lin = LinearSegmentedColormap.from_list("tc_diverging_cont", list(DIVERGING_5))
        return ListedColormap([lin(i / (steps - 1)) for i in range(steps)], name="tc_diverging")
    return ListedColormap(list(DIVERGING_5), name="tc_diverging5")


def diverging_norm(limit: float, steps: int = 5) -> BoundaryNorm:
    """Symmetric boundary norm for the discrete diverging map with clipping at ``±limit``.

    The two outer bins are open-ended (values beyond ``limit`` are clipped into them);
    the legend must state the clip limit (Fingerprint Strip guard).
    """
    if limit <= 0:
        raise TinctureStyleError("diverging_norm limit must be positive.")
    import numpy as np

    edges = np.linspace(-limit, limit, steps + 1)
    return BoundaryNorm(edges, ncolors=steps, clip=True)


def resolve_cmap(name: str) -> Colormap | str:
    """Return a colormap allowed by the theme.

    ``"diverging"`` gives the tincture 5-step diverging map. Sequential maps must be
    perceptually uniform; rainbow/jet-style maps raise :class:`TinctureStyleError`.
    """
    if name in ("diverging", "tc_diverging", "tc_diverging5"):
        return diverging_cmap()
    base = name[:-2] if name.endswith("_r") else name
    if base in BANNED_CMAPS:
        raise TinctureStyleError(
            f"Colormap {name!r} is not perceptually uniform and is not allowed. "
            f"Use one of {sorted(ALLOWED_SEQUENTIAL)} or 'diverging'."
        )
    if base not in ALLOWED_SEQUENTIAL:
        raise TinctureStyleError(
            f"Colormap {name!r} is not in the allowed set. Use one of "
            f"{sorted(ALLOWED_SEQUENTIAL)} or 'diverging'."
        )
    return name
