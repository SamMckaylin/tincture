"""Journal export presets.

Widths/heights are taken from the publishers' author guidelines as recalled during
the build; they could **not** be re-checked online from the build environment (the
publishers' sites were not reachable). See ``docs/design/journal_presets.md`` and
re-check before submission. The minimum font size is tincture's own floor (6 pt),
which is at or above every publisher's minimum.
"""

from __future__ import annotations

from dataclasses import dataclass

from tincture.errors import TinctureStyleError

MM_PER_INCH = 25.4


@dataclass(frozen=True)
class JournalPreset:
    """Physical figure size and typography for one journal format."""

    name: str
    width_mm: float
    max_height_mm: float
    min_font_pt: float = 6.0
    base_font_pt: float = 7.0
    dpi: int = 600
    note: str = ""

    @property
    def width_in(self) -> float:
        return self.width_mm / MM_PER_INCH

    @property
    def max_height_in(self) -> float:
        return self.max_height_mm / MM_PER_INCH

    def figsize(self, height_mm: float | None = None, aspect: float = 0.62) -> tuple[float, float]:
        """Figure size in inches.

        Parameters
        ----------
        height_mm
            Explicit height. Must not exceed ``max_height_mm``.
        aspect
            Height / width, used when ``height_mm`` is not given. Clamped to the
            maximum height.
        """
        if height_mm is None:
            height_mm = min(self.width_mm * aspect, self.max_height_mm)
        if height_mm <= 0:
            raise TinctureStyleError("Figure height must be positive.")
        if height_mm > self.max_height_mm + 1e-9:
            raise TinctureStyleError(
                f"Height {height_mm:.0f} mm exceeds the {self.name} maximum of "
                f"{self.max_height_mm:.0f} mm. Split the figure or use a larger preset."
            )
        return (self.width_in, height_mm / MM_PER_INCH)


_UNVERIFIED = "recalled from author guidelines; not re-verified online"

PRESETS: dict[str, JournalPreset] = {
    p.name: p
    for p in (
        JournalPreset("nature_single", 89.0, 247.0, note=_UNVERIFIED),
        JournalPreset("nature_double", 183.0, 247.0, note=_UNVERIFIED),
        JournalPreset("cell_single", 85.0, 225.0, note=_UNVERIFIED),
        JournalPreset("cell_onehalf", 114.0, 225.0, note=_UNVERIFIED),
        JournalPreset("cell_double", 174.0, 225.0, note=_UNVERIFIED),
        JournalPreset("acs_single", 82.55, 241.3, note=_UNVERIFIED + " (3.25 in)"),
        JournalPreset("acs_double", 177.8, 241.3, note=_UNVERIFIED + " (7 in)"),
        # Screen preset for notebooks / slides; not a journal format.
        JournalPreset("screen", 160.0, 400.0, base_font_pt=8.0, dpi=150, note="on-screen use"),
    )
}

#: Short aliases named in the brief.
ALIASES: dict[str, str] = {"cell": "cell_single", "acs": "acs_single", "nature": "nature_single"}


def get_preset(preset: str | JournalPreset) -> JournalPreset:
    """Look up a preset by name (or pass one through)."""
    if isinstance(preset, JournalPreset):
        return preset
    key = ALIASES.get(preset, preset)
    if key not in PRESETS:
        raise TinctureStyleError(
            f"Unknown preset {preset!r}. Available: {sorted(PRESETS) + sorted(ALIASES)}."
        )
    return PRESETS[key]
