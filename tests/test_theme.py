from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pytest

from tincture.errors import TinctureColorblindError, TinctureStyleError
from tincture.theme import (
    PRESETS,
    check_colorblind,
    check_colours,
    get_preset,
    resolve_cmap,
    resolve_color,
)
from tincture.theme import palettes as P
from tincture.theme.colorblind import delta_e2000, simulate, srgb_to_lab


class TestPalette:
    def test_semantic_names_resolve(self):
        assert resolve_color("focus") == P.TERRACOTTA
        assert resolve_color("reference") == P.BLUE
        assert resolve_color("novel") == P.PURPLE
        assert resolve_color("#123456") == "#123456"

    def test_unknown_colour_raises(self):
        with pytest.raises(TinctureStyleError, match="Unknown colour"):
            resolve_color("not-a-colour")

    @pytest.mark.parametrize("name", ["jet", "rainbow", "turbo", "hsv", "jet_r", "nipy_spectral"])
    def test_rainbow_maps_refused(self, name):
        with pytest.raises(TinctureStyleError, match="not perceptually uniform"):
            resolve_cmap(name)

    def test_unlisted_map_refused(self):
        with pytest.raises(TinctureStyleError, match="not in the allowed set"):
            resolve_cmap("coolwarm")

    @pytest.mark.parametrize(
        "name", ["viridis", "cividis", "magma", "inferno", "plasma", "viridis_r"]
    )
    def test_uniform_maps_allowed(self, name):
        assert resolve_cmap(name) == name

    def test_diverging_map_is_the_palette(self):
        cm = resolve_cmap("diverging")
        assert cm.N == 5  # type: ignore[union-attr]

    def test_diverging_norm_clips(self):
        n = P.diverging_norm(2.0)
        assert n(-10) == 0 and n(10) == 4 and n(0) == 2
        with pytest.raises(TinctureStyleError):
            P.diverging_norm(0)


class TestPresets:
    @pytest.mark.parametrize(
        "name,width",
        [("nature_single", 89), ("nature_double", 183), ("cell", 85), ("acs", 82.55)],
    )
    def test_widths(self, name, width):
        assert get_preset(name).width_mm == pytest.approx(width)

    def test_all_presets_have_6pt_floor(self):
        assert all(p.min_font_pt >= 6 for p in PRESETS.values())

    def test_height_over_max_raises(self):
        with pytest.raises(TinctureStyleError, match="exceeds"):
            get_preset("nature_double").figsize(height_mm=300)

    def test_aspect_clamped_to_max(self):
        _, h = get_preset("nature_single").figsize(aspect=10)
        assert h * 25.4 == pytest.approx(247)

    def test_unknown_preset(self):
        with pytest.raises(TinctureStyleError, match="Unknown preset"):
            get_preset("science")


class TestColourScience:
    # Sharma, Wu & Dalal (2005) CIEDE2000 test data, pairs 1 and 17.
    @pytest.mark.parametrize(
        "lab1,lab2,expected",
        [
            ((50.0, 2.6772, -79.7751), (50.0, 0.0, -82.7485), 2.0425),
            ((50.0, 2.5, 0.0), (73.0, 25.0, -18.0), 27.1492),
            ((50.0, 2.5, 0.0), (50.0, 0.0, -2.5), 4.3065),
        ],
    )
    def test_ciede2000_reference_values(self, lab1, lab2, expected):
        assert delta_e2000(np.array(lab1), np.array(lab2)) == pytest.approx(expected, abs=1e-4)

    def test_lab_of_white_and_black(self):
        assert srgb_to_lab(np.array([1.0, 1.0, 1.0])) == pytest.approx([100, 0, 0], abs=1e-3)
        assert srgb_to_lab(np.array([0.0, 0.0, 0.0])) == pytest.approx([0, 0, 0], abs=1e-6)

    def test_simulation_preserves_greys(self):
        grey = np.array([0.5, 0.5, 0.5])
        for kind in ("deuteranopia", "protanopia"):
            assert simulate(grey, kind) == pytest.approx(grey, abs=0.01)

    def test_red_green_collapse_under_deuteranopia(self):
        r = check_colours(["#D62728", "#2CA02C"], threshold=20)
        assert r.min_delta["normal"] > 20
        assert not r.passed


class TestColorblindCheck:
    def test_core_pair_passes(self):
        assert check_colours([P.BLUE, P.TERRACOTTA]).passed

    def test_brief_palette_collisions_are_detected(self):
        # Documented finding (docs/milestone_0.md): reference blue vs novel purple
        # collapse under deuteranopia, focus terracotta vs good green under protanopia.
        r = check_colours([P.BLUE, P.PURPLE, P.TERRACOTTA, P.GREEN])
        kinds = {k for *_, k, _ in [(a, b, kind, d) for a, b, kind, d in r.failures]}
        assert {"deuteranopia", "protanopia"} <= kinds

    def test_check_on_figure_with_tagged_artists(self):
        fig, ax = plt.subplots()
        (a,) = ax.plot([0, 1], [0, 1], color=P.BLUE)
        (b,) = ax.plot([0, 1], [1, 0], color=P.PURPLE)
        a.set_gid(f"tc:cat:{P.BLUE}")
        b.set_gid(f"tc:cat:{P.PURPLE}")
        with pytest.raises(TinctureColorblindError, match="deuteranopia"):
            check_colorblind(fig)
        rep = check_colorblind(fig, raise_on_fail=False)
        assert not rep.passed
        plt.close(fig)

    def test_untagged_figure_ignores_greys_and_colormaps(self):
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1], color="#888888")
        ax.plot([0, 1], [0, 1], color="black")
        ax.imshow(np.random.default_rng(0).random((4, 4)), cmap="viridis")
        ax.plot([0, 1], [1, 0], color=P.BLUE)
        ax.plot([0, 1], [1, 1], color=P.TERRACOTTA)
        rep = check_colorblind(fig)
        assert rep.passed and len(rep.colours) == 2
        plt.close(fig)
