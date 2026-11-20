"""LINCS metadata parsing (no data needed)."""

from __future__ import annotations

import math

import pytest

from tincture.datasets.cache import RemoteFile
from tincture.io.lincs import parse_dose, parse_time


@pytest.mark.parametrize(
    "raw,molar",
    [("10.0 um", 1e-5), ("0.04 um", 4e-8), ("500 nM", 5e-7), ("1 µM", 1e-6)],
)
def test_parse_dose(raw, molar):
    assert parse_dose(raw) == pytest.approx(molar)


def test_parse_dose_controls_and_time():
    assert math.isnan(parse_dose("-666")) and math.isnan(parse_dose("10 furlongs"))
    assert parse_time("24 h") == 24.0 and math.isnan(parse_time("-666"))


def test_remote_file_accepts_sha512_only(tmp_path):
    from tincture.datasets.cache import verify

    p = tmp_path / "x.txt"
    p.write_text("abc")
    import hashlib

    f = RemoteFile(
        "t", "x.txt", "https://example.invalid/x", sha512=hashlib.sha512(b"abc").hexdigest()
    )
    assert verify(f, p)
    bad = RemoteFile("t", "x.txt", "https://example.invalid/x", sha512="0" * 128)
    assert not verify(bad, p)
