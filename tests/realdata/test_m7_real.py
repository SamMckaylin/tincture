"""Opt-in real-data checks for M7 (pytest -m realdata)."""

from __future__ import annotations

import numpy as np
import pytest

from tincture.datasets.cache import is_cached
from tincture.datasets.registry import BASHI2024, JAAKS2022, MIXSEQ, VIS2024_SUPPLEMENT

pytestmark = pytest.mark.realdata


def test_jaaks_facts():
    if not is_cached(JAAKS2022["original_fitted"]):
        pytest.skip("Jaaks 2022 files not cached")
    from tincture.io.gdsc_combo import anchored_lift_table, load_jaaks2022

    j = load_jaaks2022()
    assert len(j) == 296_707 and j["CELL_LINE_NAME"].nunique() == 125
    assert j["COMBI_ID"].nunique() == 2025
    lt = anchored_lift_table(j[j["Tissue"] == "Colon"], "AZD7762", "Camptothecin")
    assert set(lt["anchor_conc_um"].round(4)) == {0.0625, 0.25}


def test_vis_facts():
    if not is_cached(VIS2024_SUPPLEMENT):
        pytest.skip("Vis 2024 supplement not cached")
    from tincture.io.gdsc_combo import load_vis2024

    v = load_vis2024()
    assert v.groupby(["ANCHOR_ID", "LIBRARY_ID"]).ngroups == 51


def test_bashi_matrix_is_full_7x7():
    if not is_cached(BASHI2024["well_stats"]):
        pytest.skip("Bashi 2024 files not cached")
    from tincture.io.gdsc_combo import load_bashi_matrices

    cm = load_bashi_matrices("Selumetinib", "Venetoclax", ["ML-2"])["ML-2"]
    assert cm.design == "full" and len(cm.levels_a) == len(cm.levels_b) == 7
    assert np.isclose(cm.levels_b.max() / cm.levels_b.min(), 1000, rtol=0.01)


def test_mixseq_pseudobulk_matches_qc_counts():
    if not is_cached(MIXSEQ["trametinib_tc"]):
        pytest.skip("MIX-seq time course not cached")
    from tincture.io.mixseq import pseudobulk_trametinib_tc

    _, g = pseudobulk_trametinib_tc()
    n = g.groupby(["condition", "time_h"])["n_cells"].sum()
    assert n[("DMSO", 3.0)] == 1098 and n[("trametinib", 48.0)] == 1439
    assert sorted(g["time_h"].unique()) == [3.0, 6.0, 12.0, 24.0, 48.0]
