"""Neighbourhood Card: statistics, null kinds, leave-one-out, the no-confident-neighbour
guard (HANDOVER §5.9), structure drawing, baseline."""

from __future__ import annotations

import numpy as np
import pytest

import tincture as tc
from tincture.errors import TinctureDataError
from tincture.model.reference import library_from_compounds
from tincture.model.signature import Signature
from tincture.spec import ImageMark, NoticeMark, RectMark, TextMark
from tincture.stats.neighbours import library_null, neighbours
from tincture.theme import check_colorblind

P = 200
FEATS = np.array([f"g{i}" for i in range(P)])
MOAS = ["synthetic MoA A", "synthetic MoA B", "synthetic MoA C", "synthetic MoA D"]


def _library(seed: int = 0, per_class: int = 6, singletons: int = 40):
    rng = np.random.default_rng(seed)
    centres = rng.normal(0, 1, (len(MOAS), P))
    sigs, ann, reps = {}, {}, {}
    for k, moa in enumerate(MOAS):
        for j in range(per_class):
            cid = f"synthetic_{k}_{j}"
            r = centres[k] + rng.normal(0, 0.6, (4, P))  # 4 "cell lines"
            reps[cid], sigs[cid] = r, r.mean(0)
            ann[cid] = {"moa": moa, "name": cid.replace("synthetic_", "cpd ")}
    for j in range(singletons):  # unrelated compounds, each its own MoA
        cid = f"synthetic_s{j}"
        r = rng.normal(0, 1, (4, P))
        reps[cid], sigs[cid] = r, r.mean(0)
        ann[cid] = {"moa": f"synthetic singleton {j}", "name": f"single {j}"}
    lib = library_from_compounds(sigs, FEATS, ann, version="synthetic v1", source="synthetic")
    return lib, sigs, reps


def _query(cid, sigs, smiles=None):
    return Signature(features=FEATS, effect=sigs[cid], compound_id=cid,
                     meta={"compound_name": cid, "smiles": smiles})  # fmt: skip


def test_null_kinds():
    lib, _, _ = _library()
    all_pairs = library_null(lib, "all")
    unrel = library_null(lib, "unrelated")
    n = len(lib.names)
    assert all_pairs.size == n * (n - 1) // 2
    assert unrel.size == all_pairs.size - len(MOAS) * 15  # 6 members → 15 same-MoA pairs each
    assert np.quantile(unrel, 0.999) < np.quantile(all_pairs, 0.999)  # related pairs removed
    with pytest.raises(ValueError):
        library_null(lib, "nope")


def test_member_query_finds_its_class_leave_one_out():
    lib, sigs, reps = _library()
    q = _query("synthetic_1_0", sigs)
    # handover default: null = all library pairs. Classmates rank first, but the class's own
    # pairs sit in the null's upper tail, so they are not called confident (conservative).
    strict = tc.neighbourhood(q, lib, query_replicates=reps["synthetic_1_0"])
    r = strict.meta["result"]
    assert "synthetic_1_0" not in r.names  # left out of its own library
    assert all(n.startswith("synthetic_1_") for n in r.names[:5])
    assert strict.meta["confident"] == []
    # null of unrelated pairs (different MoA): the same classmates are confident
    res = tc.neighbourhood(q, lib, query_replicates=reps["synthetic_1_0"], null_kind="unrelated")
    conf = res.meta["confident"]
    assert len(conf) == 5 and all(c.startswith("synthetic_1_") for c in conf)
    r2 = res.meta["result"]
    assert r2.ci is not None and np.all(r2.ci[:, 0] <= r2.similarity + 1e-9)
    assert "most resembles cpd 1_" in res.caption and "synthetic MoA B" in res.caption
    assert [m for m in res.spec.marks(RectMark) if m.categorical]  # coloured bars
    assert "different MoA" in res.spec.footer[0]


def test_no_confident_neighbour_replaces_the_list():
    lib, _, _ = _library()
    rng = np.random.default_rng(9)
    q = Signature(features=FEATS, effect=rng.normal(0, 1, P), compound_id="synthetic_query")
    res = tc.neighbourhood(q, lib, null_kind="unrelated")
    assert res.meta["confident"] == []
    assert any(isinstance(m, NoticeMark) and m.text == "no confident neighbour"
               for m in res.spec.marks())  # fmt: skip
    assert not [m for m in res.spec.marks(RectMark)]  # no ranked list
    assert "No neighbour passes q < 0.05" in res.caption
    assert "no interval" in res.caption  # no replicates, no SEs


def test_p_values_and_bh():
    lib, sigs, _ = _library()
    q = sigs["synthetic_0_0"]
    r = neighbours(q, lib.without("synthetic_0_0"), null_kind="unrelated")
    assert np.all(np.diff(r.similarity) <= 0)  # sorted
    assert np.all(r.q >= r.p - 1e-12) and np.all((r.p > 0) & (r.p <= 1))
    assert r.percentile[0] > 0.99
    with pytest.raises(TinctureDataError, match="all zeros"):
        neighbours(np.zeros(P), lib)


def test_structure_drawn_from_smiles_or_placeholder():
    lib, sigs, _ = _library()
    with_smiles = tc.neighbourhood(_query("synthetic_2_0", sigs, smiles="CC(=O)Nc1ccccc1"), lib)
    imgs = list(with_smiles.spec.marks(ImageMark))
    if imgs:  # RDKit installed
        assert np.asarray(imgs[0].image).ndim == 3
    none = tc.neighbourhood(_query("synthetic_2_0", sigs), lib)
    assert any(isinstance(m, TextMark) and m.text == "structure not available"
               for m in none.spec.marks())  # fmt: skip
    assert none.classic is not None and none.classic.plot == "neighbourhood_classic"


@pytest.mark.mpl_image_compare(tolerance=5, style="default", savefig_kwargs={"dpi": 150})
def test_neighbourhood_baseline():
    lib, sigs, reps = _library()
    fig = tc.neighbourhood(_query("synthetic_1_0", sigs), lib, null_kind="unrelated",
                           query_replicates=reps["synthetic_1_0"]).figure()  # fmt: skip
    check_colorblind(fig)
    return fig
