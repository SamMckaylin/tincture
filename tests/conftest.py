from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from pathlib import Path

import pytest

from tests.fixtures.synthetic_dose_table import synthetic_dose_table
from tincture.theme import set_font
from tincture.theme.style import DETERMINISTIC_FONT

# Baselines must not depend on which system fonts are installed.
set_font(DETERMINISTIC_FONT)


@pytest.fixture
def dose_table():
    return synthetic_dose_table()


@pytest.fixture
def dose_table_single_rep():
    return synthetic_dose_table(n_rep=1)


BASELINE_MPL = (Path(__file__).parent / "baseline" / "MATPLOTLIB_VERSION").read_text().strip()


def _major_minor(v: str) -> tuple[str, ...]:
    return tuple(v.split(".")[:2])


@pytest.hookimpl(tryfirst=True)
def pytest_configure(config):
    """Pixel baselines are only meaningful for the matplotlib version that made them.

    With another major.minor version the image tests still run (figure is built and
    colour-blind checked) but pixels are not compared.
    """
    if config.getoption("--mpl", default=False) and _major_minor(
        matplotlib.__version__
    ) != _major_minor(BASELINE_MPL):
        config.option.mpl = False
        config._tc_mpl_note = (
            f"matplotlib {matplotlib.__version__} != baseline {BASELINE_MPL}: "
            "image comparison disabled for this run."
        )


def pytest_report_header(config):
    return getattr(config, "_tc_mpl_note", f"tincture baselines: matplotlib {BASELINE_MPL}")
