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


def test_davis_and_nanobret_reproduce_vasta_2018():
    from tincture.datasets.cache import is_cached
    from tincture.datasets.registry import DAVIS2011_TABLE4, NANOBRET_VASTA2018

    if not (is_cached(DAVIS2011_TABLE4) and is_cached(NANOBRET_VASTA2018)):
        pytest.skip("kinase panel files not placed")
    from tincture.io.davis import load_davis2011, read_table4
    from tincture.io.nanobret import occupancy

    t = read_table4()
    assert len(t) == 442 and len(t.columns) - 4 == 72  # 3 id columns + kind
    occ = occupancy()
    assert occ["kinase"].nunique() == 178
    dv = load_davis2011(["Crizotinib"]).data.set_index("target")
    from tincture.io.kinase_names import to_davis

    o = occ[occ["compound_id"] == "crizotinib"].assign(t=lambda x: x["kinase"].map(to_davis))
    q = [k for k in o["t"] if k in dv.index]
    assert (o["mean"] > 0.5).sum() == 13  # paper: 13 engaged in cells at 1 µM
    kd = dv.loc[q]
    assert ((kd["kd_nm"] <= 1000) & (kd["censored"] == "")).sum() == 69  # paper: 69
