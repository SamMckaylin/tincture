"""Opt-in real-data checks for M5 (pytest -m realdata)."""

from __future__ import annotations

import numpy as np
import pytest

from tincture.datasets.cache import data_home

pytestmark = pytest.mark.realdata


def test_lincs_landmark_cache():
    cache = data_home().parent / "processed" / "lincs_gse70138_landmark.h5ad"
    if not cache.exists():
        pytest.skip("LINCS landmark cache not built")
    from tincture.io.lincs import load_lincs

    a = load_lincs(cell_lines=["MCF7"])
    assert a.n_vars == 978 and a.n_obs > 5000
    assert set(a.obs["pert_type"]) == {"trt_cp"} and np.allclose(a.obs["time"].astype(float), 24)
    assert np.isfinite(np.asarray(a.X)).all()
    vor = a.obs["pert_iname"] == "vorinostat"
    assert vor.any() and a.obs.loc[vor, "dose"].astype(float).max() == pytest.approx(1e-5)
