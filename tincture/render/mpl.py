"""Render a :class:`~tincture.spec.PlotSpec` with matplotlib.

This is the only module in tincture that creates matplotlib artists for plots.
Rendering happens inside :func:`tincture.theme.style.theme_context`, so output does
not depend on (or modify) the user's global rcParams.
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import Ellipse, Patch, Rectangle
from matplotlib.text import Text

from tincture.errors import TinctureSpecError
from tincture.spec import (
    BandMark,
    ContourMark,
    EllipseMark,
    HeatmapMark,
    IntervalMark,
    LineMark,
    NoticeMark,
    Panel,
    PlotSpec,
    PointMark,
    RectMark,
    RefLineMark,
    Scale,
    TextMark,
)
from tincture.theme.palettes import INK_SOFT, NEUTRAL, PAPER, resolve_cmap, resolve_color
from tincture.theme.presets import JournalPreset, get_preset
from tincture.theme.style import theme_context

_LS = {"solid": "-", "dashed": "--", "dotted": ":", "dashdot": "-."}

#: Gap (in row units) inserted after rows listed in ``HeatmapMark.row_gap_after``.
ROW_GAP = 0.6


def _tag(artist: Any, color: str, categorical: bool) -> None:
    if categorical and color not in ("none",):
        artist.set_gid(f"tc:cat:{color}")


def _colors(c: str | list[str], n: int) -> list[str]:
    if isinstance(c, list):
        return [resolve_color(x) for x in c]
    return [resolve_color(c)] * n


def _bcast(v: Any, n: int) -> np.ndarray:
    a = np.asarray(v)
    return np.broadcast_to(a, (n,)) if a.ndim == 0 else a


def heatmap_row_positions(n_rows: int, gaps_after: list[int]) -> np.ndarray:
    """Centre y-position of each heatmap row, with gaps inserted after given rows."""
    pos = np.arange(n_rows, dtype=float)
    for g in sorted(gaps_after):
        pos[g + 1 :] += ROW_GAP
    return pos


# ----------------------------------------------------------------------------------------
# Mark drawers
# ----------------------------------------------------------------------------------------


def _draw_points(ax: Axes, m: PointMark, handles: list) -> None:
    x, y = np.asarray(m.x, float), np.asarray(m.y, float)
    n = x.size
    cols = _colors(m.color, n)
    sizes = _bcast(m.size, n).astype(float)
    hollow = _bcast(m.hollow, n).astype(bool)
    edge = resolve_color(m.edgecolor) if m.edgecolor else None
    # group by (colour, hollow) so each artist has one categorical colour
    keys = sorted({(c, bool(h)) for c, h in zip(cols, hollow, strict=True)})
    for c, h in keys:
        idx = np.array([cc == c and bool(hh) == h for cc, hh in zip(cols, hollow, strict=True)])
        sc = ax.scatter(
            x[idx],
            y[idx],
            s=sizes[idx],
            marker=m.marker,
            facecolors="none" if h else c,
            edgecolors=c if h else (edge or c),
            linewidths=0.8 if h else 0.4,
            alpha=m.alpha,
            zorder=m.z,
        )
        _tag(sc, c, m.categorical)
    if m.labels is not None:
        for xi, yi, lab in zip(x, y, m.labels, strict=True):
            ax.annotate(
                lab,
                (xi, yi),
                xytext=(3, 2),
                textcoords="offset points",
                fontsize=plt.rcParams["xtick.labelsize"],
                zorder=m.z + 0.1,
            )
    if m.legend:
        c = cols[0] if cols else resolve_color("ink")
        handles.append(
            Line2D(
                [],
                [],
                ls="",
                marker=m.marker,
                markerfacecolor="none" if hollow.any() and hollow.all() else c,
                markeredgecolor=c,
                label=m.legend,
            )
        )


def _draw_line(ax: Axes, m: LineMark, handles: list) -> None:
    x, y = np.asarray(m.x, float), np.asarray(m.y, float)
    c = resolve_color(m.color)
    (ln,) = ax.plot(
        x,
        y,
        color=c,
        lw=m.width,
        ls=_LS[m.style],
        alpha=m.alpha,
        zorder=m.z,
        drawstyle="steps-post" if m.step else "default",
    )
    _tag(ln, c, m.categorical)
    if m.arrow_end and x.size >= 2:
        ax.annotate(
            "",
            xy=(x[-1], y[-1]),
            xytext=(x[-2], y[-2]),
            arrowprops={"arrowstyle": "-|>", "color": c, "lw": m.width, "shrinkA": 0},
            zorder=m.z,
        )
    if m.legend:
        handles.append(Line2D([], [], color=c, lw=m.width, ls=_LS[m.style], label=m.legend))


def _draw_interval(ax: Axes, m: IntervalMark, handles: list) -> None:
    pos = np.asarray(m.pos, float)
    lo, hi = np.asarray(m.lo, float), np.asarray(m.hi, float)
    cols = _colors(m.color, pos.size)
    for c in sorted(set(cols)):
        idx = np.array([cc == c for cc in cols])
        if m.orientation == "v":
            art = ax.vlines(pos[idx], lo[idx], hi[idx], colors=c, lw=m.width, zorder=m.z)
        else:
            art = ax.hlines(pos[idx], lo[idx], hi[idx], colors=c, lw=m.width, zorder=m.z)
        _tag(art, c, m.categorical)
        if m.cap > 0:
            mk = "_" if m.orientation == "v" else "|"
            for ends in (lo[idx], hi[idx]):
                xy = (pos[idx], ends) if m.orientation == "v" else (ends, pos[idx])
                ax.plot(*xy, ls="", marker=mk, ms=m.cap, color=c, mew=m.width, zorder=m.z)
    if m.legend:
        handles.append(Line2D([], [], color=cols[0], lw=m.width, label=m.legend))


def _draw_band(ax: Axes, m: BandMark, handles: list) -> None:
    c = resolve_color(m.color)
    ec = resolve_color(m.edgecolor) if m.edgecolor else "none"
    kw: dict[str, Any] = {
        "color": c,
        "alpha": m.alpha,
        "zorder": m.z,
        "lw": 0.6 if m.edgecolor else 0,
        "hatch": m.hatch,
    }
    art: Any
    if m.span:
        lo, hi = float(m.lo), float(m.hi)
        art = ax.axvspan(lo, hi, **kw) if m.orientation == "v" else ax.axhspan(lo, hi, **kw)
        art.set_edgecolor(ec)
    elif m.orientation == "v":
        art = ax.fill_between(np.asarray(m.x, float), m.lo, m.hi, **kw)
    else:
        art = ax.fill_betweenx(np.asarray(m.x, float), m.lo, m.hi, **kw)
    _tag(art, c, m.categorical)
    if m.legend:
        handles.append(Patch(facecolor=c, alpha=m.alpha, hatch=m.hatch, label=m.legend))


def _draw_rects(ax: Axes, m: RectMark, handles: list) -> None:
    x0 = np.atleast_1d(np.asarray(m.x0, float))
    n = x0.size
    y0 = _bcast(np.asarray(m.y0, float), n)
    w = _bcast(np.asarray(m.width, float), n)
    h = _bcast(np.asarray(m.height, float), n)
    cols = _colors(m.color, n)
    for i in range(n):
        ec = resolve_color(m.edgecolor) if m.edgecolor else "none"
        r = Rectangle(
            (x0[i], y0[i]),
            w[i],
            h[i],
            facecolor=cols[i],
            edgecolor=ec,
            lw=m.linewidth,
            ls=_LS[m.linestyle],
            hatch=m.hatch,
            alpha=m.alpha,
            zorder=m.z,
        )
        ax.add_patch(r)
        _tag(r, cols[i] if cols[i] != "none" else ec, m.categorical)
    if m.legend:
        handles.append(
            Patch(
                facecolor=cols[0] if n else "none",
                edgecolor=resolve_color(m.edgecolor) if m.edgecolor else "none",
                hatch=m.hatch,
                label=m.legend,
            )
        )


def _draw_heatmap(ax: Axes, m: HeatmapMark) -> Any:
    v = np.asarray(m.values, float)
    n_rows, n_cols = v.shape
    mask = np.zeros_like(v, bool) if m.mask is None else np.asarray(m.mask, bool)
    mask = mask | np.isnan(v)
    if m.cmap == "diverging":
        from tincture.theme.palettes import diverging_cmap, diverging_norm

        cmap = diverging_cmap().copy()
        norm = diverging_norm(float(m.clip))  # type: ignore[arg-type]
        vmin = vmax = None
    else:
        resolved = resolve_cmap(m.cmap)
        cmap = plt.get_cmap(resolved).copy() if isinstance(resolved, str) else resolved.copy()
        norm = None
        vmin, vmax = m.vmin, m.vmax
    cmap = cmap.with_extremes(bad=NEUTRAL)
    data = np.ma.masked_array(v, mask=mask)
    ypos = heatmap_row_positions(n_rows, m.row_gap_after)
    edge = resolve_color(m.cell_edge) if m.cell_edge else "face"
    mappable = None
    # draw contiguous row blocks separately so gaps stay empty
    starts = [0] + [g + 1 for g in sorted(m.row_gap_after)]
    ends = [g + 1 for g in sorted(m.row_gap_after)] + [n_rows]
    xe = np.arange(n_cols + 1) - 0.5
    for s, e in zip(starts, ends, strict=True):
        if s >= e:
            continue
        ye = np.concatenate([ypos[s:e] - 0.5, [ypos[e - 1] + 0.5]])
        mappable = ax.pcolormesh(
            xe,
            ye,
            data[s:e],
            cmap=cmap,
            norm=norm,
            vmin=vmin,
            vmax=vmax,
            edgecolors=edge,
            linewidth=0.3 if m.cell_edge else 0,
            zorder=m.z,
        )
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(ypos[-1] + 0.5 if n_rows else 0.5, -0.5)
    return mappable


def _draw_contour(ax: Axes, m: ContourMark) -> Any:
    resolved = resolve_cmap(m.cmap)
    x, y, z = np.asarray(m.x, float), np.asarray(m.y, float), np.asarray(m.z_values, float)
    if m.filled:
        cs = ax.contourf(x, y, z, levels=m.levels, cmap=resolved, zorder=m.z)
    else:
        lc = resolve_color(m.line_color) if m.line_color else None
        cs = ax.contour(
            x, y, z, levels=m.levels, cmap=None if lc else resolved, colors=lc, zorder=m.z
        )
    return cs


def _draw_ellipses(ax: Axes, m: EllipseMark, handles: list) -> None:
    cx = np.atleast_1d(np.asarray(m.cx, float))
    n = cx.size
    cy = _bcast(np.asarray(m.cy, float), n)
    w = _bcast(np.asarray(m.width, float), n)
    h = _bcast(np.asarray(m.height, float), n)
    ang = _bcast(np.asarray(m.angle, float), n)
    cols = _colors(m.color, n)
    for i in range(n):
        fc = to_rgba(cols[i], m.alpha) if m.fill else "none"
        e = Ellipse(
            (cx[i], cy[i]),
            w[i],
            h[i],
            angle=ang[i],
            facecolor=fc,
            edgecolor=cols[i],
            lw=m.linewidth,
            ls=_LS[m.linestyle],
            zorder=m.z,
        )
        ax.add_patch(e)
        _tag(e, cols[i], m.categorical)
    if m.legend and n:
        handles.append(Patch(facecolor="none", edgecolor=cols[0], label=m.legend))


def _draw_refline(ax: Axes, m: RefLineMark, handles: list) -> None:
    c = resolve_color(m.color)
    kw: dict[str, Any] = {"color": c, "ls": _LS[m.style], "lw": m.width, "zorder": m.z}
    ln = ax.axvline(m.value, **kw) if m.axis == "x" else ax.axhline(m.value, **kw)
    _tag(ln, c, m.categorical)
    if m.label:
        fs = plt.rcParams["xtick.labelsize"]
        if m.axis == "x":
            ax.annotate(
                m.label,
                xy=(m.value, 1.0),
                xycoords=ax.get_xaxis_transform(),
                xytext=(2, -2),
                textcoords="offset points",
                ha="left",
                va="top",
                fontsize=fs,
                color=c,
            )
        else:
            ax.annotate(
                m.label,
                xy=(1.0, m.value),
                xycoords=ax.get_yaxis_transform(),
                xytext=(-2, 1.5),
                textcoords="offset points",
                ha="right",
                va="bottom",
                fontsize=fs,
                color=c,
            )
    if m.legend:
        handles.append(Line2D([], [], color=c, ls=_LS[m.style], lw=m.width, label=m.legend))


def _draw_text(ax: Axes, m: TextMark) -> None:
    c = resolve_color(m.color)
    t = ax.text(
        m.x,
        m.y,
        m.text,
        transform=ax.transAxes if m.coords == "axes" else ax.transData,
        ha=m.ha,
        va=m.va,
        color=c,
        fontsize=m.size if m.size is not None else plt.rcParams["font.size"],
        fontweight=m.weight,
        rotation=m.rotation,
        zorder=m.z,
    )
    _tag(t, c, m.categorical)


def _draw_notice(ax: Axes, m: NoticeMark) -> None:
    ink = {"info": INK_SOFT, "weak": INK_SOFT, "warning": resolve_color("warn")}[m.level]
    at = AnchoredText(
        m.text,
        loc=m.loc,
        frameon=True,
        prop={"fontsize": plt.rcParams["xtick.labelsize"], "color": ink, "style": "italic"},
        pad=0.3,
        borderpad=0.4,
    )
    at.patch.set_facecolor(PAPER)
    at.patch.set_edgecolor(NEUTRAL)
    at.patch.set_linewidth(0.5)
    at.patch.set_alpha(0.9)
    at.set_zorder(m.z)
    at.set_gid("tc:notice")
    ax.add_artist(at)


# ----------------------------------------------------------------------------------------
# Scales and panels
# ----------------------------------------------------------------------------------------


def _apply_scale(ax: Axes, s: Scale, which: str, cat_pos: np.ndarray | None) -> None:
    set_scale = ax.set_xscale if which == "x" else ax.set_yscale
    set_lim = ax.set_xlim if which == "x" else ax.set_ylim
    set_label = ax.set_xlabel if which == "x" else ax.set_ylabel
    set_ticks = ax.set_xticks if which == "x" else ax.set_yticks
    axis = ax.xaxis if which == "x" else ax.yaxis

    if s.kind in ("log", "symlog"):
        set_scale(s.kind)
    if s.kind == "categorical":
        assert s.categories is not None
        pos = cat_pos if cat_pos is not None else np.arange(len(s.categories), dtype=float)
        if len(pos) != len(s.categories):
            raise TinctureSpecError(
                f"Categorical {which}-scale has {len(s.categories)} categories but the "
                f"marks define {len(pos)} positions."
            )
        set_ticks(pos, s.categories)
        axis.set_tick_params(length=0)
    elif s.ticks is not None:
        set_ticks(s.ticks, s.tick_labels if s.tick_labels is not None else None)
    if s.domain is not None:
        set_lim(*s.domain)
    if s.invert:
        (ax.invert_xaxis if which == "x" else ax.invert_yaxis)()
    if s.label:
        set_label(s.label)
    if s.grid:
        axis.grid(True)
    if not s.visible:
        axis.set_visible(False)
        side = ("bottom", "top") if which == "x" else ("left", "right")
        for sp in side:
            ax.spines[sp].set_visible(False)


def _draw_panel(fig: Figure, ax: Axes, panel: Panel) -> None:
    handles: list = []
    mappable = None
    heatmap_rows: np.ndarray | None = None
    for m in sorted(panel.marks, key=lambda mk: mk.z):
        if isinstance(m, PointMark):
            _draw_points(ax, m, handles)
        elif isinstance(m, LineMark):
            _draw_line(ax, m, handles)
        elif isinstance(m, IntervalMark):
            _draw_interval(ax, m, handles)
        elif isinstance(m, BandMark):
            _draw_band(ax, m, handles)
        elif isinstance(m, RectMark):
            _draw_rects(ax, m, handles)
        elif isinstance(m, HeatmapMark):
            mappable = _draw_heatmap(ax, m)
            heatmap_rows = heatmap_row_positions(np.shape(m.values)[0], m.row_gap_after)
        elif isinstance(m, ContourMark):
            mappable = _draw_contour(ax, m)
        elif isinstance(m, EllipseMark):
            _draw_ellipses(ax, m, handles)
        elif isinstance(m, RefLineMark):
            _draw_refline(ax, m, handles)
        elif isinstance(m, TextMark):
            _draw_text(ax, m)
        elif isinstance(m, NoticeMark):
            _draw_notice(ax, m)
        else:  # pragma: no cover
            raise TinctureSpecError(f"Unknown mark type {type(m).__name__}.")

    if panel.projection == "cartesian":
        _apply_scale(ax, panel.x, "x", None)
        _apply_scale(ax, panel.y, "y", heatmap_rows)
        if heatmap_rows is not None and panel.y.kind == "categorical":
            ax.set_ylim(heatmap_rows[-1] + 0.5, -0.5)
        if panel.aspect == "equal":
            ax.set_aspect("equal")
        if not panel.frame:
            for sp in ax.spines.values():
                sp.set_visible(False)
    else:  # polar
        ax.set_theta_zero_location("N")  # type: ignore[attr-defined]
        ax.set_theta_direction(-1)  # type: ignore[attr-defined]
        if panel.x.ticks is not None:
            ax.set_xticks(panel.x.ticks, panel.x.tick_labels)
        if panel.y.domain is not None:
            ax.set_ylim(*panel.y.domain)
        if panel.y.ticks is not None:
            ax.set_yticks(panel.y.ticks, panel.y.tick_labels)
        ax.spines["polar"].set_color(NEUTRAL)

    if panel.title:
        ax.set_title(panel.title)
    if panel.colorbar is not None:
        if mappable is None:
            raise TinctureSpecError("Panel has a colorbar but no heatmap/contour mark.")
        label = panel.colorbar.label
        if panel.colorbar.clip_note:
            label = f"{label}\n({panel.colorbar.clip_note})" if label else panel.colorbar.clip_note
        cb = fig.colorbar(
            mappable, ax=ax, orientation=panel.colorbar.orientation, shrink=0.8, aspect=25
        )
        cb.set_label(label)
        cb.outline.set_linewidth(0.4)
    if panel.legend is not None and handles:
        lkw: dict[str, Any] = {"title": panel.legend.title, "ncol": panel.legend.ncol}
        if panel.legend.outside:
            lkw.update(loc="upper left", bbox_to_anchor=(1.01, 1.0))
        else:
            lkw["loc"] = panel.legend.loc
        ax.legend(handles=handles, **lkw)


def render_mpl(
    spec: PlotSpec,
    preset: str | JournalPreset = "nature_double",
    height_mm: float | None = None,
) -> Figure:
    """Draw a PlotSpec as a matplotlib Figure sized for ``preset``.

    Parameters
    ----------
    spec
        The plot specification.
    preset
        Journal preset (sets width, fonts, dpi). See :mod:`tincture.theme.presets`.
    height_mm
        Explicit figure height; default is ``width × spec.aspect`` (clamped).
    """
    spec.validate()
    p = get_preset(preset)
    with theme_context(p):
        figsize = p.figsize(height_mm, spec.aspect)
        fig = plt.figure(figsize=figsize, layout="constrained")
        gs = fig.add_gridspec(
            spec.nrows,
            spec.ncols,
            width_ratios=spec.width_ratios,
            height_ratios=spec.height_ratios,
        )
        for panel in spec.panels:
            sub = gs[panel.row : panel.row + panel.rowspan, panel.col : panel.col + panel.colspan]
            ax = fig.add_subplot(sub, projection="polar" if panel.projection == "polar" else None)
            ax.set_gid(f"tc:panel:{panel.key}")
            _draw_panel(fig, ax, panel)
        if spec.title:
            fig.suptitle(spec.title, x=0.01, ha="left", fontweight="bold")
        if spec.footer:
            _add_footer(fig, spec.footer, p)
        _freeze_ticks(fig)
        # Keep the figure out of pyplot's global state (no surprise display/leaks).
        plt.close(fig)
    return fig


def _freeze_ticks(fig: Figure) -> None:
    """Pin tick styling to the theme.

    Matplotlib creates tick labels lazily (often at draw/save time, after we have left
    the theme context), reading global rcParams. Storing the values in each axis's
    tick kwargs keeps late ticks on-theme.
    """
    rc: Any = plt.rcParams
    family = rc["font.sans-serif"][0] if rc["font.family"] == ["sans-serif"] else None
    for ax in fig.axes:
        for which in ("x", "y"):
            kw: dict[str, Any] = {
                "labelsize": rc[f"{which}tick.labelsize"],
                "labelcolor": rc[f"{which}tick.color"],
                "color": rc[f"{which}tick.color"],
                "width": rc[f"{which}tick.major.width"],
            }
            if family is not None:
                kw["labelfontfamily"] = _first_available(rc["font.sans-serif"])
            ax.tick_params(axis=which, which="major", **kw)
            ax.tick_params(axis=which, which="minor", labelsize=kw["labelsize"])


def _first_available(families: list[str]) -> str:
    from matplotlib import font_manager

    names = {f.name for f in font_manager.fontManager.ttflist}
    for fam in families:
        if fam in names:
            return fam
    return "DejaVu Sans"


def _add_footer(fig: Figure, lines: list[str], p: JournalPreset) -> None:
    """Reserve space at the bottom of the figure and print footer lines there."""
    fs = max(p.min_font_pt, p.base_font_pt - 1)
    line_h_in = fs * 1.35 / 72.0
    h_in = fig.get_figheight()
    frac = min(0.4, (len(lines) * line_h_in + 0.04) / h_in)
    engine = fig.get_layout_engine()
    if engine is not None:
        engine.set(rect=(0.0, frac, 1.0, 1.0 - frac))  # type: ignore[call-arg]
    fig.text(
        0.005,
        frac - 0.02 / h_in,
        "\n".join(lines),
        ha="left",
        va="top",
        fontsize=fs,
        color=INK_SOFT,
        linespacing=1.25,
        gid="tc:footer",
    )


def iter_texts(fig: Figure) -> list[Text]:
    """All visible, non-empty text artists in a figure (tick labels included)."""
    out: list[Text] = []
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text().strip():
            continue
        out.append(t)
    return out


def min_font_size(fig: Figure) -> float:
    """Smallest font size (pt) of any visible text in the figure."""
    fig.canvas.draw()  # ensure ticks are materialised
    sizes = [t.get_fontsize() for t in iter_texts(fig)]
    return float(min(sizes)) if sizes else float("inf")
