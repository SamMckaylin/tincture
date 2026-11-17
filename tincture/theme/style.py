"""Typography and matplotlib rcParams for the tincture theme."""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from typing import Any, cast

import matplotlib as mpl
import matplotlib.style

from tincture.theme.palettes import GRID, INK, INK_SOFT, MID_GREY, PAPER
from tincture.theme.presets import JournalPreset, get_preset

#: Font preference order. Inter (SIL Open Font Licence) is the tincture look, matching the
#: design whiteboard; journals that insist on Arial/Helvetica can switch with
#: ``tc.theme.set_font("Arial")``. Liberation Sans is metric-compatible with Arial; DejaVu
#: Sans ships with matplotlib so is always available.
FONT_STACK: list[str] = ["Inter", "Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]

#: Used by tests so baselines do not depend on which system fonts are installed.
DETERMINISTIC_FONT = "DejaVu Sans"

_font_override: str | None = None


def set_font(family: str | None) -> None:
    """Force a single font family for all subsequent renders (``None`` restores the stack)."""
    global _font_override
    _font_override = family


def rc_params(preset: str | JournalPreset = "nature_double") -> dict[str, Any]:
    """rcParams implementing the tincture look at a given journal preset."""
    p = get_preset(preset)
    base = p.base_font_pt
    small = max(p.min_font_pt, base - 1)
    family = [_font_override] if _font_override else FONT_STACK
    return {
        "font.family": "sans-serif",
        "font.sans-serif": family,
        "font.size": base,
        "axes.titlesize": base + 1,
        "axes.titleweight": "normal",
        "axes.titlepad": 6.0,
        "axes.titlelocation": "left",
        "axes.labelsize": base,
        "axes.labelcolor": INK,
        "axes.edgecolor": INK,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.facecolor": PAPER,
        "axes.grid": False,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "xtick.labelsize": small,
        "ytick.labelsize": small,
        "xtick.color": MID_GREY,
        "ytick.color": MID_GREY,
        "xtick.labelcolor": INK_SOFT,
        "ytick.labelcolor": INK_SOFT,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.0,
        "ytick.major.size": 2.0,
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "legend.fontsize": small,
        "legend.title_fontsize": small,
        "legend.frameon": False,
        "lines.linewidth": 1.4,
        "lines.solid_capstyle": "round",
        "lines.solid_joinstyle": "round",
        "lines.markersize": 4.0,
        "patch.linewidth": 0.6,
        "figure.facecolor": PAPER,
        "figure.dpi": 100,
        "figure.constrained_layout.use": True,
        "savefig.dpi": p.dpi,
        "savefig.facecolor": PAPER,
        "savefig.bbox": "standard",
        "pdf.fonttype": 42,  # TrueType: text stays editable in Illustrator, as journals ask
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "text.color": INK,
        "image.cmap": "viridis",
        "polaraxes.grid": True,
        "hatch.color": INK_SOFT,
        "hatch.linewidth": 0.5,
    }


@contextmanager
def theme_context(preset: str | JournalPreset = "nature_double") -> Iterator[None]:
    """Context manager applying the tincture theme (does not leak into global state).

    Starts from matplotlib's defaults, so the output does not depend on the user's own
    style (e.g. ``classic`` or a seaborn style) — only on the preset.
    """
    with mpl.style.context("default"), mpl.rc_context(cast(Any, rc_params(preset))):
        yield


def use(preset: str | JournalPreset = "nature_double") -> None:
    """Apply the theme globally (for users who want it in their own plots)."""
    mpl.rcParams.update(cast(Any, rc_params(preset)))
