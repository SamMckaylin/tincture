"""PlotResult: what every public plot function returns."""

from __future__ import annotations

import io
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from tincture.spec import PlotSpec
from tincture.theme.presets import JournalPreset

if TYPE_CHECKING:
    from matplotlib.figure import Figure


@dataclass
class PlotResult:
    """A built plot: spec, caption, methods, warnings and an optional classic fallback.

    The figure is rendered lazily (and cached per preset/height) with
    :func:`tincture.render.mpl.render_mpl`.
    """

    spec: PlotSpec
    classic: PlotSpec | None = None
    _cache: dict[tuple[str, float | None], Figure] = field(default_factory=dict, repr=False)

    @property
    def caption(self) -> str:
        return self.spec.caption

    @property
    def methods(self) -> str:
        return self.spec.methods

    @property
    def warnings(self) -> list[str]:
        return self.spec.warnings

    @property
    def meta(self) -> dict:
        return self.spec.meta

    def figure(
        self, preset: str | JournalPreset = "nature_double", height_mm: float | None = None
    ) -> Figure:
        from tincture.render.mpl import render_mpl

        key = (preset if isinstance(preset, str) else preset.name, height_mm)
        if key not in self._cache:
            self._cache[key] = render_mpl(self.spec, preset, height_mm)
        return self._cache[key]

    def classic_figure(
        self, preset: str | JournalPreset = "nature_double", height_mm: float | None = None
    ) -> Figure | None:
        if self.classic is None:
            return None
        from tincture.render.mpl import render_mpl

        return render_mpl(self.classic, preset, height_mm)

    def save(
        self,
        path: str | Path,
        preset: str | JournalPreset = "nature_double",
        height_mm: float | None = None,
        formats: tuple[str, ...] = ("pdf", "png"),
        classic: bool = True,
        check: bool = True,
    ) -> dict[str, Path]:
        """Export figure, caption/methods sidecar and (optionally) the classic fallback."""
        from tincture.export import export

        return export(
            self,
            path,
            preset=preset,
            height_mm=height_mm,
            formats=formats,
            classic=classic,
            check=check,
        )

    def _repr_png_(self) -> bytes:  # Jupyter display
        from tincture.theme.style import theme_context

        buf = io.BytesIO()
        with theme_context("screen"):
            self.figure("screen").savefig(buf, format="png", dpi=150)
        return buf.getvalue()

    def __repr__(self) -> str:
        return (
            f"PlotResult(plot={self.spec.plot!r}, panels={len(self.spec.panels)}, "
            f"warnings={len(self.spec.warnings)}, classic={self.classic is not None})"
        )
