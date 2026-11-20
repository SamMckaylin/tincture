"""Theme: palette, typography, journal presets, colour-blind check."""

from tincture.theme.colorblind import ColorblindReport, check_colorblind, check_colours
from tincture.theme.palettes import (
    DIVERGING_5,
    SEMANTIC,
    diverging_cmap,
    resolve_cmap,
    resolve_color,
)
from tincture.theme.presets import PRESETS, JournalPreset, get_preset
from tincture.theme.style import rc_params, set_font, theme_context, use

__all__ = [
    "DIVERGING_5",
    "PRESETS",
    "SEMANTIC",
    "ColorblindReport",
    "JournalPreset",
    "check_colorblind",
    "check_colours",
    "diverging_cmap",
    "get_preset",
    "rc_params",
    "resolve_cmap",
    "resolve_color",
    "set_font",
    "theme_context",
    "use",
]
