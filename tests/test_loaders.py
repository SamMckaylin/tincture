"""Loader logic on synthetic inputs (no downloads): cache/checksums, GDSC tags, PRISM screens."""

from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd
import pytest

from tincture.datasets.cache import DatasetUnavailable, RemoteFile, fetch, is_cached, verify
from tincture.errors import TinctureDataError
from tincture.io.gdsc import collapse_duplicate_wells, tag_to_well_type
from tincture.io.prism import choose_screens, compound_core
from tincture.io.sciplex import clean_smiles


# --------------------------------------------------------------------------- cache
@pytest.fixture
def source_file(tmp_path):
    src = tmp_path / "src" / "synthetic_payload.txt"
    src.parent.mkdir()
    src.write_bytes(b"synthetic payload\n")
    return src


def _rf(src, **kw):
    data = src.read_bytes()
    base = {"md5": hashlib.md5(data).hexdigest(), "size": len(data)}
    base.update(kw)
    return RemoteFile("ds", "payload.txt", src.as_uri(), **base)


def test_fetch_downloads_and_verifies(tmp_path, source_file):
    root = tmp_path / "cache"
    f = _rf(source_file)
    p = fetch(f, root=root)
    assert p.read_bytes() == source_file.read_bytes() and verify(f, p)
    assert is_cached(f, root, check=True)
    assert fetch(f, root=root) == p  # second call uses cache


def test_checksum_mismatch_never_used(tmp_path, source_file):
    root = tmp_path / "cache"
    f = _rf(source_file, md5="0" * 32)
    with pytest.raises(DatasetUnavailable, match="checksum does not match"):
        fetch(f, root=root)
    assert not (root / "ds" / "payload.txt").exists()  # nothing half-written left behind


def test_corrupt_cached_file_refused(tmp_path, source_file):
    root = tmp_path / "cache"
    f = _rf(source_file)
    (root / "ds").mkdir(parents=True)
    (root / "ds" / "payload.txt").write_bytes(b"truncated")
    with pytest.raises(DatasetUnavailable, match="exists but its checksum"):
        fetch(f, root=root)


def test_missing_without_download_says_where(tmp_path, source_file):
    f = _rf(source_file)
    with pytest.raises(DatasetUnavailable, match="place it at"):
        fetch(f, root=tmp_path / "cache", download=False)


def test_failed_download_says_what_to_do(tmp_path):
    f = RemoteFile("ds", "x.txt", (tmp_path / "missing.txt").as_uri(), md5="0" * 32)
    with pytest.raises(DatasetUnavailable, match="Download it manually"):
        fetch(f, root=tmp_path / "cache")


def test_checksum_required():
    with pytest.raises(ValueError, match="checksum"):
        RemoteFile("ds", "x", "https://example.org/x")


# --------------------------------------------------------------------------- GDSC
def test_gdsc_tag_mapping():
    tags = pd.Series(
        ["L12-D3-S", "R1-D7-S", "T2-S", "PC1-D5-S", "PC-1", "NC-1", "NC-0", "B", "UN-USED", "DMSO"]
    )
    assert tag_to_well_type(tags).tolist() == [
        "compound",
        "compound",
        "compound",
        "pos_ctrl",
        "pos_ctrl",
        "neg_ctrl",
        "neg_ctrl",
        "blank",
        "empty",
        "empty",
    ]
    with pytest.raises(TinctureDataError, match="Unknown GDSC TAG"):
        tag_to_well_type(pd.Series(["Q9"]))


def _raw(rows):
    return pd.DataFrame(
        rows, columns=["SCAN_ID", "POSITION", "TAG", "DRUG_ID", "CONC", "INTENSITY"]
    )


def test_gdsc_collapse_solvent_pairs():
    raw = _raw(
        [
            (1, 5, "DMSO", 4000, np.nan, 100),
            (1, 5, "L3-D2-S", 1034, 0.5, 100),
            (1, 6, "DMSO", 4000, np.nan, 90),
            (1, 7, "NC-1", np.nan, np.nan, 120),
        ]
    )
    out = collapse_duplicate_wells(raw)
    assert len(out) == 3
    w5 = out[out["POSITION"] == 5].iloc[0]
    assert w5["TAG"] == "L3-D2-S" and bool(w5["solvent_paired"])
    assert out.loc[out["POSITION"] == 6, "TAG"].item() == "DMSO"  # standalone kept as-is


def test_gdsc_conflicting_duplicates_raise():
    raw = _raw([(1, 5, "L1-D1-S", 1, 0.1, 100), (1, 5, "L2-D1-S", 2, 0.1, 50)])
    with pytest.raises(TinctureDataError, match="different intensities"):
        collapse_duplicate_wells(raw)


def test_gdsc_two_drug_well_joined():
    raw = _raw([(1, 5, "T1-S", 11, 1.0, 100), (1, 5, "T2-S", 12, 1.0, 100)])
    out = collapse_duplicate_wells(raw)
    assert len(out) == 1 and out["TAG"].item() == "T1-S+T2-S" and out["DRUG_ID"].item() == "11+12"
    assert np.isnan(out["CONC"].item())


# --------------------------------------------------------------------------- PRISM
def test_prism_screen_preference():
    t = pd.DataFrame(
        {
            "compound_id": ["A", "A", "B", "C"],
            "screen_id": ["HTS002", "MTS010", "HTS002", "MTS006"],
            "perturbation_type": "experimental_treatment",
        }
    )
    s = choose_screens(t, "auto")
    assert s.to_dict() == {"A": "MTS010", "B": "HTS002", "C": "MTS006"}
    assert set(choose_screens(t, "HTS002").index) == {"A", "B"}
    with pytest.raises(TinctureDataError, match="Unknown screen"):
        choose_screens(t, "XYZ")


def test_prism_compound_core():
    assert compound_core(pd.Series(["BRD-K09951645-001-03-6"])).item() == "BRD-K09951645"


# --------------------------------------------------------------------------- sci-Plex
@pytest.mark.parametrize(
    "s,out",
    [("CCO |r,c:1|", "CCO"), ("CCO", "CCO"), ("NA", None), (np.nan, None), ("", None)],
)
def test_clean_smiles(s, out):
    assert clean_smiles(s) == out
