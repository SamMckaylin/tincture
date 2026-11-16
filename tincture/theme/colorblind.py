"""Colour-vision-deficiency check for figures.

``check_colorblind(fig)`` collects the *categorical* colours used in a figure,
simulates deuteranopia and protanopia (Machado, Oliveira & Fernandes 2009, severity
1.0, applied in linear RGB) and reports any pair of colours whose CIEDE2000 distance
falls below a threshold under any of normal vision, deuteranopia or protanopia.

What counts as categorical
--------------------------
* Artists created by the tincture renderer are tagged with a ``gid`` of the form
  ``"tc:cat:<colour>"``; when any tagged artist exists, only tagged colours are checked.
* Otherwise (any matplotlib figure) the colours of lines, patches and unmapped
  collections are used. Achromatic colours (greys, black, white), fully transparent
  colours and colormapped artists (images, mapped scatter) are ignored: sequential
  maps are checked by construction (only perceptually uniform maps are allowed).
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
from matplotlib.collections import Collection
from matplotlib.colors import to_hex, to_rgba
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from tincture.errors import TinctureColorblindError

#: Default minimum CIEDE2000 distance between two categorical colours.
DEFAULT_DELTA_E = 10.0

# Machado et al. 2009, severity 1.0 (linear RGB).
_MACHADO = {
    "protanopia": np.array(
        [
            [0.152286, 1.052583, -0.204868],
            [0.114503, 0.786281, 0.099216],
            [-0.003882, -0.048116, 1.051998],
        ]
    ),
    "deuteranopia": np.array(
        [
            [0.367322, 0.860646, -0.227968],
            [0.280085, 0.672501, 0.047413],
            [-0.011820, 0.042940, 0.968881],
        ]
    ),
}


def _srgb_to_linear(c: np.ndarray) -> np.ndarray:
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _linear_to_srgb(c: np.ndarray) -> np.ndarray:
    c = np.clip(c, 0.0, 1.0)
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1 / 2.4) - 0.055)


def simulate(rgb: np.ndarray, kind: str) -> np.ndarray:
    """Simulate a colour-vision deficiency on sRGB colours in [0, 1], shape (..., 3)."""
    if kind == "normal":
        return np.asarray(rgb, dtype=float)
    m = _MACHADO[kind]
    lin = _srgb_to_linear(np.asarray(rgb, dtype=float))
    return _linear_to_srgb(lin @ m.T)


def srgb_to_lab(rgb: np.ndarray) -> np.ndarray:
    """sRGB in [0, 1] -> CIE L*a*b* (D65)."""
    lin = _srgb_to_linear(np.asarray(rgb, dtype=float))
    m = np.array(
        [
            [0.4124564, 0.3575761, 0.1804375],
            [0.2126729, 0.7151522, 0.0721750],
            [0.0193339, 0.1191920, 0.9503041],
        ]
    )
    xyz = lin @ m.T
    white = np.array([0.95047, 1.0, 1.08883])
    t = xyz / white
    eps, kappa = 216 / 24389, 24389 / 27
    f = np.where(t > eps, np.cbrt(t), (kappa * t + 16) / 116)
    lab_l = 116 * f[..., 1] - 16
    lab_a = 500 * (f[..., 0] - f[..., 1])
    lab_b = 200 * (f[..., 1] - f[..., 2])
    return np.stack([lab_l, lab_a, lab_b], axis=-1)


def delta_e2000(lab1: np.ndarray, lab2: np.ndarray) -> float:
    """CIEDE2000 colour difference (Sharma, Wu & Dalal 2005)."""
    l1, a1, b1 = (float(v) for v in lab1)
    l2, a2, b2 = (float(v) for v in lab2)
    c1, c2 = np.hypot(a1, b1), np.hypot(a2, b2)
    cbar = (c1 + c2) / 2
    g = 0.5 * (1 - np.sqrt(cbar**7 / (cbar**7 + 25.0**7)))
    a1p, a2p = (1 + g) * a1, (1 + g) * a2
    c1p, c2p = np.hypot(a1p, b1), np.hypot(a2p, b2)
    h1p = np.degrees(np.arctan2(b1, a1p)) % 360
    h2p = np.degrees(np.arctan2(b2, a2p)) % 360
    dlp = l2 - l1
    dcp = c2p - c1p
    if c1p * c2p == 0:
        dhp = 0.0
    elif abs(h2p - h1p) <= 180:
        dhp = h2p - h1p
    elif h2p - h1p > 180:
        dhp = h2p - h1p - 360
    else:
        dhp = h2p - h1p + 360
    dhp_big = 2 * np.sqrt(c1p * c2p) * np.sin(np.radians(dhp / 2))
    lbarp = (l1 + l2) / 2
    cbarp = (c1p + c2p) / 2
    if c1p * c2p == 0:
        hbarp = h1p + h2p
    elif abs(h1p - h2p) <= 180:
        hbarp = (h1p + h2p) / 2
    elif h1p + h2p < 360:
        hbarp = (h1p + h2p + 360) / 2
    else:
        hbarp = (h1p + h2p - 360) / 2
    t = (
        1
        - 0.17 * np.cos(np.radians(hbarp - 30))
        + 0.24 * np.cos(np.radians(2 * hbarp))
        + 0.32 * np.cos(np.radians(3 * hbarp + 6))
        - 0.20 * np.cos(np.radians(4 * hbarp - 63))
    )
    dtheta = 30 * np.exp(-(((hbarp - 275) / 25) ** 2))
    rc = 2 * np.sqrt(cbarp**7 / (cbarp**7 + 25.0**7))
    sl = 1 + 0.015 * (lbarp - 50) ** 2 / np.sqrt(20 + (lbarp - 50) ** 2)
    sc = 1 + 0.045 * cbarp
    sh = 1 + 0.015 * cbarp * t
    rt = -np.sin(np.radians(2 * dtheta)) * rc
    return float(
        np.sqrt(
            (dlp / sl) ** 2
            + (dcp / sc) ** 2
            + (dhp_big / sh) ** 2
            + rt * (dcp / sc) * (dhp_big / sh)
        )
    )


def _is_achromatic(rgba: tuple[float, float, float, float], tol: float = 0.04) -> bool:
    r, g, b, _ = rgba
    return max(r, g, b) - min(r, g, b) < tol


def _collect_colours(fig: Figure) -> list[str]:
    tagged: set[str] = set()
    untagged: set[str] = set()

    def add(dest: set[str], c: object) -> None:
        try:
            arr = np.atleast_2d(np.asarray(to_rgba_array_safe(c)))
        except (ValueError, TypeError):
            return
        for row in arr:
            rgba = tuple(float(x) for x in row)
            if len(rgba) != 4 or rgba[3] <= 0.0:
                continue
            if _is_achromatic(rgba):
                continue
            dest.add(to_hex(rgba[:3]))

    for ax in fig.get_axes():
        for art in ax.get_children():
            gid = art.get_gid() or ""
            if gid.startswith("tc:cat:"):
                tagged.add(to_hex(to_rgba(gid.split(":", 2)[2])[:3]))
                continue
            if isinstance(art, Line2D):
                add(untagged, art.get_color())
                add(untagged, art.get_markerfacecolor())
            elif isinstance(art, Patch):
                if art is ax.patch:
                    continue
                add(untagged, art.get_facecolor())
                add(untagged, art.get_edgecolor())
            elif isinstance(art, Collection):
                if art.get_array() is not None:  # colormapped
                    continue
                add(untagged, art.get_facecolor())
                add(untagged, art.get_edgecolor())
    return sorted(tagged) if tagged else sorted(untagged)


def to_rgba_array_safe(c: object) -> np.ndarray:
    from matplotlib.colors import to_rgba_array

    if isinstance(c, str) and c == "none":
        return np.zeros((0, 4))
    return to_rgba_array(c)


@dataclass
class ColorblindReport:
    """Result of :func:`check_colorblind`."""

    colours: list[str]
    threshold: float
    failures: list[tuple[str, str, str, float]] = field(default_factory=list)
    min_delta: dict[str, float] = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return not self.failures

    def __str__(self) -> str:
        if self.passed:
            return f"Colour-blind check passed for {len(self.colours)} categorical colours."
        lines = [f"Colour-blind check FAILED (ΔE2000 < {self.threshold}):"]
        lines += [f"  {a} vs {b} under {kind}: ΔE = {d:.1f}" for a, b, kind, d in self.failures]
        return "\n".join(lines)


def check_colours(
    colours: list[str],
    threshold: float = DEFAULT_DELTA_E,
    kinds: tuple[str, ...] = ("normal", "deuteranopia", "protanopia"),
) -> ColorblindReport:
    """Check a list of colours pairwise under each simulated vision type."""
    uniq = sorted({to_hex(to_rgba(c)[:3]) for c in colours})
    report = ColorblindReport(colours=uniq, threshold=threshold)
    if len(uniq) < 2:
        return report
    rgb = np.array([to_rgba(c)[:3] for c in uniq])
    for kind in kinds:
        lab = srgb_to_lab(simulate(rgb, kind))
        worst = np.inf
        for i, j in itertools.combinations(range(len(uniq)), 2):
            d = delta_e2000(lab[i], lab[j])
            worst = min(worst, d)
            if d < threshold:
                report.failures.append((uniq[i], uniq[j], kind, d))
        report.min_delta[kind] = float(worst)
    return report


def check_colorblind(
    fig: Figure,
    threshold: float = DEFAULT_DELTA_E,
    raise_on_fail: bool = True,
) -> ColorblindReport:
    """Check that a figure's categorical colours stay distinguishable for colour-blind readers.

    Parameters
    ----------
    fig
        A matplotlib figure (tincture-rendered or not).
    threshold
        Minimum CIEDE2000 difference between any two categorical colours.
    raise_on_fail
        Raise :class:`TinctureColorblindError` if any pair collides (default), otherwise
        return the report.
    """
    report = check_colours(_collect_colours(fig), threshold=threshold)
    if raise_on_fail and not report.passed:
        raise TinctureColorblindError(str(report))
    return report
