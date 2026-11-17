from __future__ import annotations

import matplotlib.pyplot as plt
import pytest
from matplotlib.image import imread

import tincture as tc
from tincture.errors import TinctureColorblindError, TinctureStyleError
from tincture.spec import Panel, PlotSpec, PointMark, TextMark
from tincture.theme import palettes as P


def test_demo_exports_at_nature_double(tmp_path, dose_table):
    """M0 acceptance: a dummy plot exports at nature_double size."""
    res = tc.demo(dose_table)
    out = tc.export(res, tmp_path / "fig1", preset="nature_double", formats=("pdf", "png", "svg"))
    assert {"pdf", "png", "svg", "md", "classic_pdf", "classic_png", "classic_svg"} <= set(out)
    img = imread(out["png"])
    dpi = tc.theme.get_preset("nature_double").dpi
    assert img.shape[1] == pytest.approx(183 / 25.4 * dpi, abs=1)
    pdf = out["pdf"].read_bytes()
    assert b"/MediaBox [ 0 0 518.74" in pdf  # 183 mm in points
    md = out["md"].read_text()
    assert "## Caption" in md and "## Methods" in md and "nature_double" in md
    assert res.caption in md


def test_extension_in_path_is_replaced(tmp_path, dose_table):
    out = tc.export(tc.demo(dose_table), tmp_path / "f.png", formats=("png",), classic=False)
    assert out["png"].name == "f.png" and "classic_png" not in out


def test_too_small_text_refused(tmp_path):
    p = Panel().add(PointMark(x=[0], y=[0]), TextMark(x=0, y=0, text="tiny", size=4))
    with pytest.raises(TinctureStyleError, match="below the 6 pt"):
        tc.export(PlotSpec(plot="x", panels=[p]), tmp_path / "x")


def test_colour_collision_refused(tmp_path):
    # The brief's original purple collapses onto the reference blue under deuteranopia.
    p = Panel().add(
        PointMark(x=[0], y=[0], color="reference", categorical=True),
        PointMark(x=[1], y=[1], color=P.BRIEF_PURPLE, categorical=True),
    )
    with pytest.raises(TinctureColorblindError):
        tc.export(PlotSpec(plot="x", panels=[p]), tmp_path / "x")
    # check=False lets the user override deliberately
    assert "png" in tc.export(PlotSpec(plot="x", panels=[p]), tmp_path / "y", check=False)


def test_reference_vs_novel_now_exports(tmp_path):
    p = Panel().add(
        PointMark(x=[0], y=[0], color="reference", categorical=True),
        PointMark(x=[1], y=[1], color="novel", categorical=True),
        PointMark(x=[2], y=[2], color="focus", categorical=True),
        PointMark(x=[3], y=[3], color="good", categorical=True),
    )
    assert "png" in tc.export(PlotSpec(plot="x", panels=[p]), tmp_path / "z", formats=("png",))


def test_plain_matplotlib_figure_resized(tmp_path):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot([0, 1], [0, 1], color=P.BLUE)
    ax.tick_params(labelsize=7)
    out = tc.export(fig, tmp_path / "plain", preset="nature_single", formats=("png",))
    assert fig.get_figwidth() * 25.4 == pytest.approx(89)
    assert fig.get_figheight() * 25.4 == pytest.approx(44.5)
    assert out["png"].exists()
    plt.close(fig)


def test_unsupported_format(tmp_path, dose_table):
    with pytest.raises(TinctureStyleError, match="Unsupported"):
        tc.export(tc.demo(dose_table), tmp_path / "x", formats=("jpg",))
