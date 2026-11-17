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
from tincture.plots.demo import demo
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
    "demo",
    "export",
    "plate_atlas",
    "render_mpl",
    "retrievability",
    "theme",
]
