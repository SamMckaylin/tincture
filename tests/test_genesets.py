"""Gene-set fetching without decoupler (core install)."""

from __future__ import annotations

import gzip
import io
import json

import pandas as pd

from tincture.datasets import genesets


class _Resp(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _fake_table() -> bytes:
    df = pd.DataFrame(
        {
            "geneset": ["HALLMARK_SYNTH_A", "HALLMARK_SYNTH_A", "HALLMARK_SYNTH_B"],
            "genesymbol": ["G1", "COMPLEX:G2_G3", "G1"],
            "collection": ["hallmark"] * 3,
        }
    )
    return gzip.compress(df.to_csv(sep="\t", index=False).encode())


def test_fetch_hallmark_matches_decoupler_processing(monkeypatch):
    monkeypatch.setattr("urllib.request.urlopen", lambda *a, **k: _Resp(_fake_table()))
    t = genesets.fetch_hallmark()
    assert set(map(tuple, t.values)) == {
        ("SYNTH_A", "G1"),
        ("SYNTH_A", "G2"),
        ("SYNTH_A", "G3"),
        ("SYNTH_B", "G1"),
    }


def test_load_geneset_caches_with_checksum_without_decoupler(monkeypatch, tmp_path):
    monkeypatch.setattr("urllib.request.urlopen", lambda *a, **k: _Resp(_fake_table()))
    lib = genesets.load_geneset("hallmark", root=tmp_path)
    assert lib.n_sets == 2 and "sha256" in lib.version
    meta = json.loads((tmp_path / "genesets" / "hallmark_human.json").read_text())
    assert meta["fetched_with"].startswith("urllib")
    # second load uses the cache (no network)
    monkeypatch.setattr("urllib.request.urlopen", lambda *a, **k: (_ for _ in ()).throw(OSError))
    assert genesets.load_geneset("hallmark", root=tmp_path).n_sets == 2
