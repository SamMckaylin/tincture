"""Errors and warnings raised by tincture.

Every user-facing failure carries a *fix suggestion*: principle 7 of the brief is
"degrade gracefully", so when we cannot draw something we say why and what to do.
"""

from __future__ import annotations


class TinctureError(Exception):
    """Base class for tincture errors."""


class TinctureDataError(TinctureError, ValueError):
    """Input data cannot support the requested plot or statistic.

    Parameters
    ----------
    message
        What is wrong with the data.
    fix
        A concrete suggestion for how to fix it (always shown to the user).
    """

    def __init__(self, message: str, fix: str | None = None) -> None:
        self.message = message
        self.fix = fix
        text = message if fix is None else f"{message}\n  How to fix: {fix}"
        super().__init__(text)


class TinctureSpecError(TinctureError, ValueError):
    """A PlotSpec is internally inconsistent (a bug in a plot builder or user spec)."""


class TinctureStyleError(TinctureError, ValueError):
    """A styling request violates the theme rules (e.g. a rainbow colormap)."""


class TinctureColorblindError(TinctureError, AssertionError):
    """Raised by ``check_colorblind(..., raise_on_fail=True)`` when colours collide."""


class TinctureWarning(UserWarning):
    """A result is weak, reduced or degraded; the plot still renders but says so."""
