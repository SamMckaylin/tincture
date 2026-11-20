from __future__ import annotations

import io

import matplotlib as mpl
import numpy as np
import pytest

from tincture.errors import TinctureSpecError
from tincture.render.mpl import heatmap_row_positions, min_font_size, render_mpl
from tincture.spec import (
    BandMark,
    ColorbarSpec,
    ContourMark,
    EllipseMark,
    HeatmapMark,
    IntervalMark,
    LegendSpec,
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
from tincture.theme import check_colorblind, theme_context


def kitchen_sink() -> PlotSpec:
    """Every mark type at least once (synthetic numbers, unit test only)."""
    xs = np.array([1e-8, 1e-7, 1e-6, 1e-5])
    a = Panel(
        key="a",
        x=Scale("log", label="Concentration (M)"),
        y=Scale(label="Response"),
        legend=LegendSpec(),
    )
    a.add(
        BandMark(lo=2e-7, hi=2e-6, span=True, color="window", legend="window"),
        BandMark(x=xs, lo=[0.8, 0.7, 0.4, 0.2], hi=[1.0, 0.9, 0.6, 0.4], color="reference_light"),
        LineMark(x=xs, y=[0.9, 0.8, 0.5, 0.3], color="reference", arrow_end=True, categorical=True),
        IntervalMark(pos=xs, lo=[0.85, 0.7, 0.4, 0.2], hi=[0.95, 0.9, 0.6, 0.4], color="reference"),
        PointMark(
            x=xs,
            y=[0.9, 0.8, 0.5, 0.3],
            color="focus",
            hollow=[0, 0, 1, 1],
            categorical=True,
            legend="focus",
            labels=["a", "b", "c", "d"],
        ),
        RefLineMark(axis="y", value=1.0, label="vehicle"),
        RectMark(
            x0=[1e-8],
            y0=[0.0],
            width=[9e-8],
            height=[0.1],
            color="none",
            edgecolor="ink",
            linestyle="dashed",
        ),
        TextMark(x=0.4, y=0.05, text="note", coords="axes"),
        NoticeMark(text="too few doses for a path", loc="lower right"),
    )
    b = Panel(
        key="b",
        x=Scale("categorical", categories=["p1", "p2", "p3"]),
        y=Scale("categorical", categories=["r1", "r2", "c1", "c2"]),
        colorbar=ColorbarSpec(label="score", clip_note="clipped at ±2"),
    )
    b.add(
        HeatmapMark(
            values=np.array([[1, -1, 0], [2, 3, -3], [0.1, 0.2, 0.3], [np.nan, 1, 1]]),
            mask=np.array([[0, 0, 1], [0, 0, 0], [0, 0, 0], [0, 0, 0]], bool),
            clip=2,
            row_gap_after=[1],
        )
    )
    g = np.linspace(0, 1, 6)
    c = Panel(key="c", x=Scale(label="A"), y=Scale(label="B"), colorbar=ColorbarSpec(label="v"))
    c.add(
        ContourMark(x=g, y=g, z_values=np.add.outer(g, g), levels=5),
        EllipseMark(
            cx=[0.5],
            cy=[0.5],
            width=[0.3],
            height=[0.1],
            angle=[30],
            color="focus",
            fill=True,
            categorical=True,
        ),
    )
    d = Panel(
        key="d",
        projection="polar",
        x=Scale(
            "angle",
            ticks=list(np.linspace(0, 2 * np.pi, 4, endpoint=False)),
            tick_labels=["A", "B", "C", "D"],
        ),
        y=Scale("radius", domain=(0, 1)),
    )
    d.add(PointMark(x=[0.3, 2.0], y=[0.5, 0.8], color="focus", categorical=True))
    return PlotSpec(
        plot="kitchen_sink",
        panels=[a, b, c, d],
        nrows=2,
        ncols=2,
        title="Kitchen sink",
        footer=["Reference set: synthetic v0", "n = 4"],
        aspect=0.8,
    )


@pytest.mark.parametrize("preset", ["nature_single", "nature_double", "cell_double", "acs_single"])
def test_figure_is_exact_preset_width(preset):
    from tincture.theme import get_preset

    fig = render_mpl(kitchen_sink(), preset=preset)
    assert fig.get_figwidth() * 25.4 == pytest.approx(get_preset(preset).width_mm)


def test_height_mm_respected():
    fig = render_mpl(kitchen_sink(), "nature_double", height_mm=120)
    assert fig.get_figheight() * 25.4 == pytest.approx(120)


def test_min_font_6pt_after_save_outside_theme():
    """Ticks materialised at save time (outside the theme context) stay ≥ 6 pt."""
    fig = render_mpl(kitchen_sink(), "nature_double")
    with mpl.rc_context({"xtick.labelsize": 3, "ytick.labelsize": 3, "font.size": 3}):
        fig.savefig(io.BytesIO(), format="png", dpi=50)
        assert min_font_size(fig) >= 6


def test_does_not_leak_rcparams():
    before = dict(mpl.rcParams)
    render_mpl(kitchen_sink())
    assert mpl.rcParams["font.size"] == before["font.size"]
    assert mpl.rcParams["axes.spines.top"] == before["axes.spines.top"]


def test_footer_notice_and_clip_note_present():
    fig = render_mpl(kitchen_sink())
    texts = [t.get_text() for t in fig.findobj(mpl.text.Text)]
    assert any("Reference set: synthetic v0" in t for t in texts)
    assert any("too few doses for a path" in t for t in texts)
    assert any("clipped at ±2" in t for t in texts)


def test_kitchen_sink_passes_colorblind():
    check_colorblind(render_mpl(kitchen_sink()))


def test_heatmap_row_gaps():
    pos = heatmap_row_positions(4, [1])
    assert pos[2] - pos[1] > 1 and pos[1] - pos[0] == 1


def test_categorical_mismatch_raises():
    p = Panel(x=Scale("categorical", categories=["a", "b"]), y=Scale())
    p.add(PointMark(x=[0, 1], y=[0, 1]))
    p.y = Scale("categorical", categories=["only-one"])
    p.add(HeatmapMark(values=np.zeros((3, 2)), clip=1))
    with pytest.raises(TinctureSpecError, match="categories"):
        render_mpl(PlotSpec(plot="x", panels=[p]))


def test_colorbar_without_mappable_raises():
    p = Panel(colorbar=ColorbarSpec(label="x"))
    p.add(PointMark(x=[0], y=[0]))
    with pytest.raises(TinctureSpecError, match="colorbar"):
        render_mpl(PlotSpec(plot="x", panels=[p]))


def test_output_independent_of_user_style():
    """A user's global style (here 'classic') must not change exported output."""
    import matplotlib.style

    def png(style):
        with matplotlib.style.context(style):
            buf = io.BytesIO()
            fig = render_mpl(kitchen_sink())
            # tc.export / PlotResult.save save inside the theme context like this.
            with theme_context():
                fig.savefig(buf, format="png", dpi=40)
            return buf.getvalue()

    assert png("classic") == png("default")


def test_rounded_sequential_heatmap_is_not_categorical():
    """Regression: tiles of a sequential map must not be checked as categorical colours."""
    v = np.linspace(0, 1, 64).reshape(8, 8)
    p = Panel(colorbar=ColorbarSpec(label="v")).add(
        HeatmapMark(values=v, cmap="viridis", style="rounded")
    )
    fig = render_mpl(PlotSpec(plot="x", panels=[p]))
    check_colorblind(fig)  # must not raise on near-identical viridis tiles
