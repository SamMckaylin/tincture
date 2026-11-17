"""tincture: question-first, publication-quality plots for drug-discovery data.

Public API (grows milestone by milestone)::

    import tincture as tc
    res = tc.demo(df)                 # every plot returns a PlotResult
    res.caption, res.methods          # templated reading guide + methods paragraph
    tc.export(res, "fig1", preset="nature_double")
"""

from tincture import theme
from tincture._version import __version__
from tincture.errors import TinctureDataError, TinctureError, TinctureWarning
from tincture.export import export
from tincture.plots.compass import plot as compass
from tincture.plots.demo import demo
from tincture.plots.dose_paths import plot as dose_paths
from tincture.plots.effect_ladder import plot as effect_ladder
from tincture.plots.fingerprint import plot as fingerprint
from tincture.plots.plate_atlas import plot as plate_atlas
from tincture.plots.retrievability import plot as retrievability
from tincture.render.mpl import render_mpl
from tincture.result import PlotResult

__all__ = [
    "PlotResult",
    "TinctureDataError",
    "TinctureError",
    "TinctureWarning",
    "__version__",
    "compass",
    "demo",
    "dose_paths",
    "effect_ladder",
    "export",
    "fingerprint",
    "plate_atlas",
    "render_mpl",
    "retrievability",
    "theme",
]
