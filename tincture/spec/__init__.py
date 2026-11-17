"""PlotSpec: a backend-neutral description of a figure.

Plot builders (``tincture.plots.*.build_spec``) compute statistics and return a
:class:`PlotSpec` made only of plain data (numbers, strings, numpy arrays). The
renderer (``tincture.render.mpl.render_mpl``) is the only code that touches
matplotlib. This keeps an interactive backend possible later and keeps the stats
layer free of plotting objects.

Colours are semantic names from :mod:`tincture.theme.palettes` (``"focus"``,
``"reference"``, ``"novel"`` ...) or hex strings.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

from tincture.errors import TinctureSpecError

ArrayLike = Any  # sequence of numbers or numpy array; normalised with np.asarray

LineStyle = Literal["solid", "dashed", "dotted", "dashdot"]
Coords = Literal["data", "axes"]

# ----------------------------------------------------------------------------------------
# Scales
# ----------------------------------------------------------------------------------------


@dataclass
class Scale:
    """One axis.

    ``kind``:
      * ``linear`` / ``log`` / ``symlog`` — continuous.
      * ``categorical`` — positions ``0..n-1`` labelled by ``categories``.
      * ``angle`` / ``radius`` — the two axes of a polar panel (angle in radians).
    """

    kind: Literal["linear", "log", "symlog", "categorical", "angle", "radius"] = "linear"
    label: str = ""
    domain: tuple[float, float] | None = None
    ticks: list[float] | None = None
    tick_labels: list[str] | None = None
    categories: list[str] | None = None
    invert: bool = False
    visible: bool = True
    grid: bool = False
    tick_rotation: float = 0.0  # degrees; rotated labels are anchored at their tick
    tick_side: Literal["default", "top"] = "default"  # x only: ticks and label on top

    def __post_init__(self) -> None:
        if self.kind == "categorical" and not self.categories:
            raise TinctureSpecError("A categorical scale needs `categories`.")
        if self.domain is not None and self.kind == "log" and min(self.domain) <= 0:
            raise TinctureSpecError("A log scale domain must be positive.")
        if (
            self.ticks is not None
            and self.tick_labels is not None
            and len(self.ticks) != len(self.tick_labels)
        ):
            raise TinctureSpecError("ticks and tick_labels differ in length.")


# ----------------------------------------------------------------------------------------
# Marks
# ----------------------------------------------------------------------------------------


@dataclass
class Mark:
    """Base class for all marks.

    ``categorical=True`` marks the colour as carrying categorical meaning, so it is
    included in the colour-blind check. ``legend`` adds a legend entry.
    """

    z: float = 2.0
    legend: str | None = None
    categorical: bool = False

    def validate(self) -> None:  # pragma: no cover - overridden
        return None


def _same_len(name: str, **arrays: ArrayLike) -> int:
    lengths = {k: np.size(v) for k, v in arrays.items() if v is not None and np.ndim(v) > 0}
    if len(set(lengths.values())) > 1:
        raise TinctureSpecError(f"{name}: arrays have different lengths {lengths}.")
    return next(iter(lengths.values()), 0)


@dataclass
class PointMark(Mark):
    """Scatter points. ``size`` is marker area in pt²; scalars broadcast."""

    x: ArrayLike = field(default_factory=list)
    y: ArrayLike = field(default_factory=list)
    size: ArrayLike = 12.0
    color: str | list[str] = "ink"
    edgecolor: str | None = None
    marker: str = "o"
    hollow: ArrayLike = False
    alpha: float = 1.0
    labels: list[str] | None = None  # optional per-point text labels (drawn in the point colour)

    def validate(self) -> None:
        n = _same_len("PointMark", x=self.x, y=self.y, size=self.size, hollow=self.hollow)
        if isinstance(self.color, list) and len(self.color) != n:
            raise TinctureSpecError("PointMark: colour list length differs from points.")
        if self.labels is not None and len(self.labels) != n:
            raise TinctureSpecError("PointMark: labels length differs from points.")


@dataclass
class LineMark(Mark):
    """Polyline through (x, y). ``arrow_end`` draws an arrowhead on the last segment."""

    x: ArrayLike = field(default_factory=list)
    y: ArrayLike = field(default_factory=list)
    color: str = "ink"
    width: float = 1.0
    style: LineStyle = "solid"
    alpha: float = 1.0
    arrow_end: bool = False
    step: bool = False
    end_label: str | None = None  # direct label at the line's last point, in its colour

    def validate(self) -> None:
        _same_len("LineMark", x=self.x, y=self.y)


@dataclass
class IntervalMark(Mark):
    """Error bars / CI whiskers: at position ``pos`` an interval ``lo``..``hi``.

    ``orientation="h"`` draws horizontal intervals (pos is y), ``"v"`` vertical.
    ``color`` may be a list (one per interval), e.g. grey when a CI crosses zero.
    """

    pos: ArrayLike = field(default_factory=list)
    lo: ArrayLike = field(default_factory=list)
    hi: ArrayLike = field(default_factory=list)
    orientation: Literal["h", "v"] = "v"
    color: str | list[str] = "ink"
    width: float = 0.8
    cap: float = 0.0

    def validate(self) -> None:
        n = _same_len("IntervalMark", pos=self.pos, lo=self.lo, hi=self.hi)
        if isinstance(self.color, list) and len(self.color) != n:
            raise TinctureSpecError("IntervalMark: colour list length differs from intervals.")
        lo, hi = np.asarray(self.lo, float), np.asarray(self.hi, float)
        ok = ~(np.isnan(lo) | np.isnan(hi))
        if np.any(lo[ok] > hi[ok]):
            raise TinctureSpecError("IntervalMark: some lo > hi.")


@dataclass
class BandMark(Mark):
    """Shaded band between ``lo`` and ``hi`` along ``x`` (``orientation="v"``) or along
    ``y`` (``"h"``). With ``span=True`` and scalar lo/hi, shades a full-width strip
    across the axis (e.g. a Cmax range or a null band)."""

    x: ArrayLike = None
    lo: ArrayLike = 0.0
    hi: ArrayLike = 0.0
    orientation: Literal["v", "h"] = "v"
    span: bool = False
    color: str = "neutral"
    alpha: float = 0.35
    hatch: str | None = None
    edgecolor: str | None = None
    edge_style: LineStyle = "dashed"  # outline style when ``edgecolor`` is given
    inner_label: str | None = None  # direct label inside the band (top-left)
    z: float = 0.5

    def validate(self) -> None:
        if self.span:
            if np.ndim(self.lo) != 0 or np.ndim(self.hi) != 0:
                raise TinctureSpecError("BandMark(span=True) needs scalar lo and hi.")
        else:
            if self.x is None:
                raise TinctureSpecError("BandMark needs x unless span=True.")
            _same_len("BandMark", x=self.x, lo=self.lo, hi=self.hi)


@dataclass
class RectMark(Mark):
    """Axis-aligned rectangles (x0, y0, width, height), e.g. a reachable-dose box."""

    x0: ArrayLike = field(default_factory=list)
    y0: ArrayLike = field(default_factory=list)
    width: ArrayLike = field(default_factory=list)
    height: ArrayLike = field(default_factory=list)
    color: str | list[str] = "neutral"
    edgecolor: str | None = None
    linestyle: LineStyle = "solid"
    linewidth: float = 0.6
    hatch: str | None = None
    alpha: float = 1.0
    clip: bool = True  # False: may draw outside the axes (e.g. a row marker left of a heatmap)

    def validate(self) -> None:
        _same_len("RectMark", x0=self.x0, y0=self.y0, width=self.width, height=self.height)


@dataclass
class HeatmapMark(Mark):
    """A matrix drawn cell by cell at categorical/integer positions.

    ``values`` has shape (n_rows, n_cols); row 0 is drawn at the top. ``mask`` cells
    (True) are drawn in the neutral colour (e.g. non-significant). ``cmap`` is
    ``"diverging"`` or an allowed sequential map; ``clip`` sets symmetric diverging
    limits that must be reported in the legend/colourbar.
    """

    values: ArrayLike = field(default_factory=lambda: np.zeros((0, 0)))
    mask: ArrayLike = None
    cmap: str = "diverging"
    clip: float | None = None
    vmin: float | None = None
    vmax: float | None = None
    row_gap_after: list[int] = field(default_factory=list)
    col_gap_after: list[int] = field(default_factory=list)
    cell_edge: str | None = "paper"
    style: Literal["auto", "rounded", "flat"] = "auto"  # rounded cells with gutters (≤ 5,000)
    circles: ArrayLike = None  # bool (n_rows, n_cols): draw these cells as circles
    z: float = 1.0

    def validate(self) -> None:
        v = np.asarray(self.values)
        if v.ndim != 2:
            raise TinctureSpecError("HeatmapMark.values must be 2-D.")
        if self.mask is not None and np.shape(self.mask) != v.shape:
            raise TinctureSpecError("HeatmapMark.mask shape differs from values.")
        if self.cmap == "diverging" and self.clip is None:
            raise TinctureSpecError(
                "HeatmapMark with a diverging cmap needs `clip` (the legend must show it)."
            )


@dataclass
class ContourMark(Mark):
    """Contours of z on a grid. ``filled=True`` gives filled contours.

    The renderer never extrapolates: only the convex hull of the given grid is drawn.
    """

    x: ArrayLike = field(default_factory=list)
    y: ArrayLike = field(default_factory=list)
    z_values: ArrayLike = field(default_factory=lambda: np.zeros((0, 0)))
    levels: ArrayLike = 8
    cmap: str = "viridis"
    filled: bool = True
    line_color: str | None = None
    z: float = 1.0

    def validate(self) -> None:
        zz = np.asarray(self.z_values)
        if zz.shape != (np.size(self.y), np.size(self.x)):
            raise TinctureSpecError("ContourMark: z_values must have shape (len(y), len(x)).")


@dataclass
class EllipseMark(Mark):
    """Ellipses (e.g. bootstrap CI regions). Angles in degrees."""

    cx: ArrayLike = field(default_factory=list)
    cy: ArrayLike = field(default_factory=list)
    width: ArrayLike = field(default_factory=list)
    height: ArrayLike = field(default_factory=list)
    angle: ArrayLike = 0.0
    color: str | list[str] = "ink"
    fill: bool = False
    alpha: float = 0.25
    linestyle: LineStyle = "solid"
    linewidth: float = 0.6

    def validate(self) -> None:
        _same_len("EllipseMark", cx=self.cx, cy=self.cy, width=self.width, height=self.height)


@dataclass
class RefLineMark(Mark):
    """A full-length reference line at ``value`` (``axis="x"`` → vertical line)."""

    axis: Literal["x", "y"] = "y"
    value: float = 0.0
    color: str = "ink_soft"
    style: LineStyle = "dashed"
    width: float = 0.6
    label: str | None = None  # small text placed at the end of the line
    z: float = 1.5


@dataclass
class TextMark(Mark):
    """Free text at (x, y) in data or axes coordinates."""

    x: float = 0.0
    y: float = 0.0
    text: str = ""
    coords: Coords = "data"
    ha: Literal["left", "center", "right"] = "left"
    va: Literal["top", "center", "bottom", "baseline"] = "center"
    color: str = "ink"
    size: float | None = None  # pt; None → theme default
    weight: Literal["normal", "bold"] = "normal"
    rotation: float = 0.0
    z: float = 5.0


@dataclass
class NoticeMark(Mark):
    """An in-panel notice that a result is weak or a plot was reduced.

    Principle 8: say when a result is weak. The notice is drawn inside the panel in a
    subtle box so that it travels with the figure (captions can get separated).
    """

    text: str = ""
    loc: Literal["upper left", "upper right", "lower left", "lower right", "center"] = "upper right"
    level: Literal["info", "weak", "warning"] = "weak"
    z: float = 6.0


MarkT = (
    PointMark
    | LineMark
    | IntervalMark
    | BandMark
    | RectMark
    | HeatmapMark
    | ContourMark
    | EllipseMark
    | RefLineMark
    | TextMark
    | NoticeMark
)

# ----------------------------------------------------------------------------------------
# Legend / colourbar / panel / figure
# ----------------------------------------------------------------------------------------


@dataclass
class ColorbarSpec:
    """Colourbar for a heatmap/contour mark in the same panel."""

    label: str = ""
    clip_note: str | None = None  # e.g. "clipped at ±3"; required for clipped maps
    orientation: Literal["vertical", "horizontal"] = "vertical"
    ticks: list[float] | None = None
    tick_labels: list[str] | None = None
    compact: bool = False  # small key (whiteboard "down ▪▪▪▪▪ up" style)


@dataclass
class LegendSpec:
    loc: str = "best"
    title: str | None = None
    ncol: int = 1
    outside: bool = False


@dataclass
class Panel:
    """One axes."""

    marks: list[MarkT] = field(default_factory=list)
    x: Scale = field(default_factory=Scale)
    y: Scale = field(default_factory=Scale)
    projection: Literal["cartesian", "polar"] = "cartesian"
    title: str = ""
    legend: LegendSpec | None = None
    colorbar: ColorbarSpec | None = None
    aspect: Literal["auto", "equal"] = "auto"
    row: int = 0
    col: int = 0
    rowspan: int = 1
    colspan: int = 1
    frame: bool = True
    key: str = ""  # identifier used in tests / meta

    def add(self, *marks: MarkT) -> Panel:
        self.marks.extend(marks)
        return self


@dataclass
class PlotSpec:
    """A complete figure description.

    Attributes
    ----------
    plot
        Name of the plot that built this spec (e.g. ``"compass"``).
    panels
        Panels placed on an ``nrows × ncols`` grid.
    caption
        Auto-generated reading-guide caption (templated from computed statistics).
    methods
        Methods paragraph.
    footer
        Short lines printed under the figure (reference set + version, n, model ...).
    warnings
        Weak-result / degraded-plot messages; also shown in the caption sidecar.
    meta
        Computed statistics and flags, for tests and for downstream use.
    aspect
        Preferred height / width ratio when no explicit height is requested.
    """

    plot: str
    panels: list[Panel] = field(default_factory=list)
    nrows: int = 1
    ncols: int = 1
    width_ratios: list[float] | None = None
    height_ratios: list[float] | None = None
    title: str = ""
    caption: str = ""
    methods: str = ""
    footer: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    meta: dict[str, Any] = field(default_factory=dict)
    aspect: float = 0.62

    def __post_init__(self) -> None:
        # Convenience: several panels all left at (0, 0) are laid out row-major.
        if len(self.panels) > 1 and all(p.row == 0 and p.col == 0 for p in self.panels):
            for i, p in enumerate(self.panels):
                p.row, p.col = divmod(i, self.ncols)

    def validate(self) -> PlotSpec:
        """Check internal consistency; returns self so it can be chained."""
        if not self.panels:
            raise TinctureSpecError(f"PlotSpec {self.plot!r} has no panels.")
        if self.width_ratios is not None and len(self.width_ratios) != self.ncols:
            raise TinctureSpecError("width_ratios length must equal ncols.")
        if self.height_ratios is not None and len(self.height_ratios) != self.nrows:
            raise TinctureSpecError("height_ratios length must equal nrows.")
        occupied: set[tuple[int, int]] = set()
        for p in self.panels:
            if p.row + p.rowspan > self.nrows or p.col + p.colspan > self.ncols:
                raise TinctureSpecError(f"Panel {p.key or p.title!r} lies outside the grid.")
            cells = {
                (r, c)
                for r in range(p.row, p.row + p.rowspan)
                for c in range(p.col, p.col + p.colspan)
            }
            if cells & occupied:
                raise TinctureSpecError(f"Panel {p.key or p.title!r} overlaps another panel.")
            occupied |= cells
            if p.projection == "polar" and (p.x.kind != "angle" or p.y.kind != "radius"):
                raise TinctureSpecError("A polar panel needs x=Scale('angle'), y=Scale('radius').")
            for m in p.marks:
                m.validate()
        return self

    def panel(self, key: str) -> Panel:
        for p in self.panels:
            if p.key == key:
                return p
        raise KeyError(key)

    def marks(self, kind: type | None = None) -> list[MarkT]:
        """All marks across panels, optionally filtered by type (for tests)."""
        out: list[MarkT] = [m for p in self.panels for m in p.marks]
        return [m for m in out if kind is None or isinstance(m, kind)]


__all__ = [
    "BandMark",
    "ColorbarSpec",
    "ContourMark",
    "EllipseMark",
    "HeatmapMark",
    "IntervalMark",
    "LegendSpec",
    "LineMark",
    "Mark",
    "NoticeMark",
    "Panel",
    "PlotSpec",
    "PointMark",
    "RectMark",
    "RefLineMark",
    "Scale",
    "TextMark",
]
