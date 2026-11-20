"""M2 on real data (opt-in: pytest -m realdata). Skipped when data are absent."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from tincture.datasets.cache import data_home, is_cached
from tincture.datasets.registry import GDSC2_RAW

pytestmark = pytest.mark.realdata
CASES = Path(__file__).resolve().parents[2] / "case_studies"


def test_cs3_plate_qc_report_subset():
    if not is_cached(GDSC2_RAW):
        pytest.skip("GDSC2 raw data not cached")
    sys.path.insert(0, str(CASES))
    import cs3_plate_qc

    path, ok = cs3_plate_qc.run(max_plates=300)
    assert ok and path.exists()
    text = path.read_text()
    assert "Z′ distribution" in text and "Edge-vs-interior ratio reported" in text


def test_retrievability_gate_on_sciplex_a549():
    cache = data_home().parent / "processed" / "sciplex3_pseudobulk_wells.h5ad"
    if not cache.exists():
        pytest.skip("sci-Plex3 pseudobulk not built")
    import tincture as tc
    from tincture.io.sciplex import load_sciplex3
    from tincture.stats.signatures import replicate_profiles

    prof = replicate_profiles(load_sciplex3(cell_lines=["A549"]))
    res = tc.retrievability(prof, focus=["S1090|1e-06"])
    t = res.meta["table"].set_index("group")
    assert t.loc["S1090|1e-06", "retrievable"]  # abexinostat 1 µM
    assert 0.2 < res.meta["n_retrievable"] / res.meta["n_conditions"] < 0.6
    assert prof.adata.uns["tincture_retrievable"]["retrievable_groups"]


def test_hdaci_pathways_on_sciplex():
    cache = data_home().parent / "processed" / "sciplex3_pseudobulk_wells.h5ad"
    if not cache.exists():
        pytest.skip("sci-Plex3 pseudobulk not built")
    import pandas as pd

    from tincture.datasets.genesets import load_geneset
    from tincture.io.sciplex import load_sciplex3
    from tincture.stats.pathways import symbols_for, ulm

    s = load_sciplex3(cell_lines=["A549"]).signatures()[0]
    Y = s.matrix()
    Y.columns = symbols_for(s.features, s.feature_symbols)
    Y = Y.loc[:, ~pd.Index(Y.columns).duplicated()]
    sc = ulm(Y, load_geneset("hallmark").table).score.loc["S1090@1e-06"]
    # HDAC inhibition arrests the cell cycle: E2F targets and G2/M strongly down
    assert sc["E2F_TARGETS"] < -5 and sc["G2M_CHECKPOINT"] < -5
