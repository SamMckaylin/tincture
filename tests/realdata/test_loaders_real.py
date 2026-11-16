"""Loader round-trips on the real downloaded data (opt-in: ``pytest -m realdata``).

Skipped, not failed, when the data are absent. Expected counts are the verified facts
recorded in docs/data_status.md.
"""

from __future__ import annotations

import numpy as np
import pytest

from tincture.datasets.cache import data_home, is_cached
from tincture.datasets.registry import GDSC2_RAW, PRISM_SECONDARY, SCIPLEX3

pytestmark = pytest.mark.realdata


def _need(*files):
    missing = [f.name for f in files if not is_cached(f)]
    if missing:
        pytest.skip(f"real data not cached: {missing}")


# --------------------------------------------------------------------------- D1 sci-Plex
def test_sciplex_metadata_matches_verified_facts():
    _need(SCIPLEX3["pdata"], SCIPLEX3["cell_annotations"], SCIPLEX3["hash_metadata"])
    from tincture.io.sciplex import _well_obs, read_cell_metadata

    cells = read_cell_metadata(download=False)
    assert len(cells) == 799_317 and int((~cells["assigned"]).sum()) == 36_522
    a = cells[cells["assigned"]]
    assert a.loc[a["dose"] > 0, "compound_id"].nunique() == 188
    assert sorted(a["dose"].unique()) == [0.0, 1e-8, 1e-7, 1e-6, 1e-5]
    assert set(a["cell_line"]) == {"A549", "MCF7", "K562"}
    wells = _well_obs(cells)
    assert len(wells) == 4992
    assert wells.loc[wells["dose"] > 0, "smiles"].notna().all()


def test_sciplex_pseudobulk_round_trip():
    cache = data_home().parent / "processed" / "sciplex3_pseudobulk_wells.h5ad"
    if not cache.exists():
        pytest.skip("sci-Plex3 pseudobulk not built (tincture.io.sciplex.build_pseudobulk)")
    from tincture.io.sciplex import load_sciplex3
    from tincture.model import PerturbationSet

    ps = load_sciplex3(time=None)
    assert isinstance(ps, PerturbationSet) and ps.shape == (4992, 58_347)
    assert ps.var_names.str.startswith("ENSG").all()  # human genes only
    # summed human counts can never exceed the per-well UMI total (which includes mouse)
    counts = np.asarray(ps.X.sum(axis=1)).ravel()
    assert (counts <= ps.obs["total_umi"].to_numpy() + 1e-6).all()
    assert counts.sum() > 0.9 * ps.obs["total_umi"].sum()
    ps24 = load_sciplex3(cell_lines=["A549"], compounds=["Abexinostat (PCI-24781)"])
    assert set(ps24.obs["time"]) == {24.0}
    assert ps24.compounds == ["S1090"] and ps24.is_vehicle.sum() == 32
    assert ps24.obs.loc[ps24.is_vehicle, "smiles"].isna().all()  # no "nan" strings
    small = ps24.subset(features=list(ps24.var_names[:20]))
    back = type(small).from_long(small.to_long())
    assert back.shape == small.shape


# --------------------------------------------------------------------------- D2 PRISM
def test_prism_round_trip():
    _need(*PRISM_SECONDARY.values())
    from tincture.io.prism import load_prism_secondary
    from tincture.model import DoseResponse

    s = load_prism_secondary(compounds=["dabrafenib"], download=False)
    assert set(s.points["screen_id"]) == {"MTS010"}  # MTS010 preferred when available
    assert s.points["dose"].max() <= 1.3e-5 and s.points["dose"].min() > 0
    assert s.points.groupby("context_id")["dose"].nunique().max() == 8
    k = s.fits.iloc[0]
    dr = s.get(k["compound_id"], k["context_id"])
    assert isinstance(dr, DoseResponse) and dr.fit is not None and dr.fit_source
    assert 0 < np.median(s.fits["top"]) < 1.5  # top = zero-dose asymptote (README swapped)
    assert not s.points["context_id"].str.endswith("FAILED_STR").any()


# --------------------------------------------------------------------------- D3 GDSC
def test_gdsc_round_trip():
    _need(GDSC2_RAW)
    from tincture.io.gdsc import load_gdsc2_plates
    from tincture.model import PlateTable

    pt = load_gdsc2_plates(max_plates=20)
    assert isinstance(pt, PlateTable) and pt.shape == (32, 48) and len(pt.plates) == 20
    d = pt.data
    blanks = d[d["well_type"] == "blank"]
    assert (blanks["col"] == 0).all()  # verified row-major layout: blanks in column 1
    rel = d["value"] / d["plate"].map(d[d["tag"] == "NC-1"].groupby("plate")["value"].median())
    assert rel[d["well_type"] == "blank"].median() < 0.1
    assert rel[d["well_type"] == "pos_ctrl"].median() < 0.4
    assert d["dose"].max() < 1e-2  # molar
