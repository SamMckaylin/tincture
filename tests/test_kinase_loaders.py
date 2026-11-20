"""Davis / NanoBRET loader logic that needs no data files."""

from __future__ import annotations

import numpy as np
import pytest

from tincture.datasets.cache import DatasetUnavailable, RemoteFile, fetch
from tincture.io.davis import classify
from tincture.io.kinase_names import to_davis
from tincture.io.nanobret import apparent_ic50


@pytest.mark.parametrize(
    "construct,kind",
    [
        ("ABL1(T315I)-phosphorylated", "mutant"),
        ("EGFR(L747-E749del, A750P)", "mutant"),
        ("FLT3(ITD)", "mutant"),
        ("ABL1-nonphosphorylated", "wildtype"),
        ("JAK1(JH1domain-catalytic)", "wildtype"),
        ("RSK1(Kin.Dom.2-C-terminal)", "wildtype"),
        ("PFCDPK1(P.falciparum)", "nonhuman"),
        ("AAK1", "wildtype"),
    ],
)
def test_construct_classification(construct, kind):
    assert classify(construct) == kind


def test_apparent_ic50_from_occupancy():
    # 1:1 binding: 50% occupancy at C → IC50 = C; 75% → C/3; clipped at 5–95%
    assert apparent_ic50(np.array([0.5]), 1e-6)[0] == pytest.approx(1e-6)
    assert apparent_ic50(np.array([0.75]), 1e-6)[0] == pytest.approx(1e-6 / 3)
    assert apparent_ic50(np.array([0.0]), 1e-6)[0] == pytest.approx(1e-6 * 0.95 / 0.05)


def test_synonyms():
    assert to_davis("MLTK") == "ZAK" and to_davis("SIK3") == "KIAA0999" and to_davis("ALK") == "ALK"


def test_manual_file_gives_instructions(tmp_path):
    f = RemoteFile("x", "y.xls", "https://example.invalid/page", sha256="0" * 64, manual=True)
    with pytest.raises(DatasetUnavailable, match="in a browser"):
        fetch(f, root=tmp_path)
