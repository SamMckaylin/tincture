"""Signatures, EB moderation, shrinkage, pathway scores (synthetic, known answers)."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from tests.fixtures.synthetic_counts import synthetic_count_set
from tincture.datasets.cache import DatasetUnavailable
from tincture.datasets.genesets import load_geneset
from tincture.errors import TinctureDataError
from tincture.stats.ebayes import fit_f_dist, squeeze_var, trigamma, trigamma_inverse
from tincture.stats.pathways import dedup_gene_sets, symbols_for, ulm
from tincture.stats.signatures import logcpm, replicate_profiles, shrink_normal


# --------------------------------------------------------------------------- EB
@pytest.mark.parametrize("x", [0.01, 0.5, 3.0, 50.0])
def test_trigamma_inverse(x):
    assert float(trigamma(trigamma_inverse(x))) == pytest.approx(x, rel=1e-6)


def test_fit_f_dist_recovers_prior():
    rng = np.random.default_rng(1)
    d0, s0, df, n = 8.0, 0.5, 4, 20000
    sigma2 = s0 * d0 / rng.chisquare(d0, n)
    s2 = sigma2 * rng.chisquare(df, n) / df
    fit = fit_f_dist(s2, df)
    assert fit.d0 == pytest.approx(d0, rel=0.2)
    assert float(np.median(fit.s0_sq)) == pytest.approx(s0, rel=0.1)


def test_squeeze_shrinks_toward_prior():
    rng = np.random.default_rng(2)
    s2 = 0.5 * rng.chisquare(3, 5000) / 3
    sq = squeeze_var(s2, 3)
    lo, hi = np.minimum(s2, sq.var_prior), np.maximum(s2, sq.var_prior)
    assert np.all((sq.var_post >= lo - 1e-12) & (sq.var_post <= hi + 1e-12))
    assert sq.var_post.var() < s2.var()


def test_fit_f_dist_needs_data():
    with pytest.raises(ValueError):
        fit_f_dist(np.array([1.0, np.nan]), 3)


# --------------------------------------------------------------------------- signatures
def test_known_fold_change_recovered():
    ps = synthetic_count_set()
    sets = ps.signatures(shrink=False, min_cpm=1)
    assert len(sets) == 1
    s = sets[0]
    a = s.get("synthetic_cpdA", 1e-6)
    up = np.isin(a.features, [f"ENSG{i:05d}" for i in range(20)])
    assert np.median(a.effect[up]) == pytest.approx(1.0, abs=0.2)
    assert abs(np.median(a.effect[~up])) < 0.1
    half = s.get("synthetic_cpdA", 1e-7)
    assert np.median(half.effect[up]) == pytest.approx(0.5, abs=0.2)
    b = s.get("synthetic_cpdB", 1e-6)
    p = s.pvalues[("synthetic_cpdB", 1e-6)]
    assert (p < 0.01).mean() < 0.03  # null compound: ~no false positives
    assert b.n_replicates == 3 and b.meta["n_vehicle"] == 8


def test_shrinkage_shrinks_noise_more_than_signal():
    s = synthetic_count_set().signatures(min_cpm=1)[0]
    a = s.get("synthetic_cpdA", 1e-6)
    raw, post = a.meta["raw_effect"], a.effect
    assert np.all(np.abs(post) <= np.abs(raw) + 1e-12)
    up = np.isin(a.features, [f"ENSG{i:05d}" for i in range(20)])
    keep_signal = np.median(np.abs(post[up]) / np.abs(raw[up]))
    keep_noise = np.median(np.abs(post[~up]) / np.abs(raw[~up]))
    assert keep_signal > keep_noise
    assert a.ci()[0].shape == a.effect.shape  # uncertainty available


def test_strata_never_mix_cell_lines():
    ps = synthetic_count_set(cell_lines=("L1", "L2"))
    sets = ps.signatures(min_cpm=1)
    assert [x.stratum["cell_line"] for x in sets] == ["L1", "L2"]
    assert all(x.n_vehicle == 8 for x in sets)


def test_no_vehicle_raises():
    ps = synthetic_count_set()
    with pytest.raises(TinctureDataError, match="No vehicle"):
        ps.subset(obs=(ps.obs["dose"] > 0).to_numpy()).signatures()


def test_no_replicates_raises():
    ps = synthetic_count_set(n_rep=1, n_vehicle=1)
    with pytest.raises(TinctureDataError, match="No replicates"):
        ps.signatures(min_cpm=1)


def test_expression_filter_too_strict():
    with pytest.raises(TinctureDataError, match="pass the expression filter"):
        synthetic_count_set().signatures(min_cpm=1e9)


def test_logcpm_zero_library():
    with pytest.raises(TinctureDataError, match="zero total"):
        logcpm(np.zeros((2, 3)))


def test_shrink_normal_null_and_signal():
    rng = np.random.default_rng(0)
    se = np.full(2000, 1.0)
    null = rng.normal(0, 1, 2000)
    post, _, _ = shrink_normal(null, se)
    assert np.abs(post).mean() < 0.1 * np.abs(null).mean()  # pure noise shrinks to ~0
    big = rng.normal(0, 5, 2000) + rng.normal(0, 1, 2000)
    post2, _, tau2b = shrink_normal(big, se)
    assert tau2b == pytest.approx(25, rel=0.15)
    assert np.abs(post2).mean() > 0.9 * np.abs(big).mean()


def test_replicate_profiles_label_free_and_shape():
    ps = synthetic_count_set()
    prof = replicate_profiles(ps, min_cpm=1, n_features=100)
    assert prof.shape == (ps.shape[0], 100)
    veh = prof.X[prof.is_vehicle]
    assert np.abs(veh.mean(axis=0)).max() < 1e-9  # centred on vehicle
    with pytest.raises(TinctureDataError, match="vehicle samples"):
        one = ps.subset(obs=((ps.obs["dose"] > 0) | (ps.obs["replicate"] == "r0")).to_numpy())
        replicate_profiles(one, min_cpm=1)


# --------------------------------------------------------------------------- pathways
def _net():
    return pd.DataFrame(
        {
            "source": ["UP"] * 20 + ["NULL"] * 20 + ["TINY"] * 2,
            "target": [f"G{i}" for i in range(20)]
            + [f"G{i}" for i in range(100, 120)]
            + ["G1", "G2"],
        }
    )


def test_ulm_detects_planted_set_and_matches_ols():
    rng = np.random.default_rng(3)
    feats = [f"G{i}" for i in range(400)]
    y = rng.normal(0, 1, 400)
    y[:20] += 2
    Y = pd.DataFrame([y], columns=feats, index=["sig"])
    sc = ulm(Y, _net())
    assert "TINY" not in sc.score.columns  # min_size guard
    assert sc.score.loc["sig", "UP"] > 6 and abs(sc.score.loc["sig", "NULL"]) < 3
    import statsmodels.api as sm

    x = np.isin(feats, [f"G{i}" for i in range(20)]).astype(float)
    ols = sm.OLS(y, sm.add_constant(x)).fit()
    assert sc.score.loc["sig", "UP"] == pytest.approx(ols.tvalues[1], rel=1e-9)
    assert sc.pvalue.loc["sig", "UP"] == pytest.approx(ols.pvalues[1], rel=1e-6)


def test_ulm_no_overlap_raises():
    Y = pd.DataFrame([[1.0, 2.0]], columns=["X1", "X2"])
    with pytest.raises(TinctureDataError, match="symbols vs Ensembl"):
        ulm(Y, _net())


def test_ulm_matches_decoupler():
    dc = pytest.importorskip("decoupler")
    rng = np.random.default_rng(4)
    Y = pd.DataFrame(
        rng.normal(size=(3, 400)), columns=[f"G{i}" for i in range(400)], index=["a", "b", "c"]
    )
    net = _net()[lambda d: d["source"] != "TINY"]
    res = dc.mt.ulm(data=Y, net=net, tmin=5, verbose=False)
    theirs = res[0] if isinstance(res, tuple) else res
    mine = ulm(Y, net).score
    assert np.allclose(mine[theirs.columns].to_numpy(), theirs.loc[Y.index].to_numpy(), atol=1e-9)


def test_dedup_keeps_larger():
    net = pd.DataFrame(
        {
            "source": ["BIG"] * 10 + ["SMALL"] * 7 + ["OTHER"] * 5,
            "target": [f"G{i}" for i in range(10)]
            + [f"G{i}" for i in range(7)]
            + [f"H{i}" for i in range(5)],
        }
    )
    out, dropped = dedup_gene_sets(net)
    assert set(out["source"]) == {"BIG", "OTHER"}
    assert (
        dropped[0][0] == "SMALL" and dropped[0][1] == "BIG" and dropped[0][2] == pytest.approx(0.7)
    )


def test_symbols_for_handles_duplicates():
    out = symbols_for(np.array(["E1", "E2", "E3"]), np.array(["A", "A", "B"]))
    assert out.tolist() == ["A", "E2", "B"]


# --------------------------------------------------------------------------- gene-set cache
def test_geneset_cache_and_checksum(tmp_path):
    d = tmp_path / "genesets"
    d.mkdir()
    tsv = d / "hallmark_human.tsv"
    pd.DataFrame({"source": ["S"] * 5, "target": list("ABCDE")}).to_csv(tsv, sep="\t", index=False)
    import hashlib

    meta = {
        "sha256": hashlib.sha256(tsv.read_bytes()).hexdigest(),
        "fetched": "2026-01-01",
        "source": "synthetic",
    }
    (d / "hallmark_human.json").write_text(json.dumps(meta))
    lib = load_geneset("hallmark", root=tmp_path, download=False)
    assert lib.n_sets == 1 and "fetched 2026-01-01" in lib.version
    tsv.write_text("source\ttarget\nS\tZ\n")
    with pytest.raises(DatasetUnavailable, match="checksum"):
        load_geneset("hallmark", root=tmp_path, download=False)
    with pytest.raises(DatasetUnavailable, match="not cached"):
        load_geneset("progeny", root=tmp_path, download=False)
    with pytest.raises(TinctureDataError, match="Unknown gene-set"):
        load_geneset("kegg", root=tmp_path)
