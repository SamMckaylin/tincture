"""Data model: long-table validation, PerturbationSet, Signature, dose-response, plates."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from tests.fixtures.synthetic_long import synthetic_long_table
from tincture.errors import TinctureDataError
from tincture.model import (
    DoseResponse,
    DoseResponseSet,
    FourPL,
    PerturbationSet,
    PlateTable,
    Signature,
    add_constant_dose,
    position_to_rowcol,
    row_label,
    validate_long,
)
from tincture.model.measurements import round_dose

# --------------------------------------------------------------------------- long table


class TestValidateLong:
    def test_ok_and_dtypes(self):
        d = validate_long(synthetic_long_table())
        assert d["dose"].dtype == float
        assert pd.api.types.is_string_dtype(d["replicate"])

    def test_missing_column_has_fix(self):
        with pytest.raises(TinctureDataError, match="How to fix"):
            validate_long(synthetic_long_table().drop(columns=["dose"]))

    def test_not_a_dataframe(self):
        with pytest.raises(TinctureDataError, match="DataFrame"):
            validate_long([1, 2, 3])  # type: ignore[arg-type]

    def test_non_numeric_dose(self):
        d = synthetic_long_table()
        d["dose"] = d["dose"].astype(object)
        d.loc[0, "dose"] = "1 uM"
        with pytest.raises(TinctureDataError, match="not numeric"):
            validate_long(d)

    def test_negative_dose(self):
        d = synthetic_long_table()
        d.loc[0, "dose"] = -1e-6
        with pytest.raises(TinctureDataError, match="Negative"):
            validate_long(d)

    def test_micromolar_mistake_caught(self):
        d = synthetic_long_table()
        d["dose"] = d["dose"] * 1e6
        with pytest.raises(TinctureDataError, match="looks like µM"):
            validate_long(d)

    def test_unknown_well_type(self):
        d = synthetic_long_table().assign(well_type="sample")
        with pytest.raises(TinctureDataError, match="well_type"):
            validate_long(d)

    def test_require_vehicle(self):
        d = synthetic_long_table()
        with pytest.raises(TinctureDataError, match="vehicle"):
            validate_long(d[d["dose"] > 0], require_vehicle=True)

    def test_inf_rejected(self):
        d = synthetic_long_table()
        d.loc[0, "value"] = np.inf
        with pytest.raises(TinctureDataError, match="Infinite"):
            validate_long(d)

    def test_nan_values_kept(self):
        d = synthetic_long_table()
        d.loc[d["feature"] == "g0", "value"] = np.nan  # all-NaN feature
        assert validate_long(d)["value"].isna().sum() == (d["feature"] == "g0").sum()

    def test_no_dose_column_must_be_stated(self):
        d = synthetic_long_table().drop(columns=["dose"])
        with pytest.raises(TinctureDataError, match="no dose"):
            validate_long(add_constant_dose(d))
        assert len(validate_long(add_constant_dose(d, 0.0)))

    def test_round_dose(self):
        assert round_dose(np.array([0.00977 / 1e6, 0.0]))[0] == 9.77e-09


# --------------------------------------------------------------------------- PerturbationSet


class TestPerturbationSet:
    def test_long_round_trip(self):
        long = synthetic_long_table()
        ps = PerturbationSet.from_long(long, source="synthetic")
        back = ps.to_long()
        key = ["compound_id", "dose", "replicate", "feature"]
        a = long.sort_values(key).reset_index(drop=True)
        b = back.sort_values(key).reset_index(drop=True)
        assert np.allclose(a["value"], b["value"])
        assert ps.shape == (long.groupby(["compound_id", "dose", "replicate"]).ngroups, 5)

    def test_repr_and_props(self):
        ps = PerturbationSet.from_long(synthetic_long_table())
        assert ps.compounds == ["synthetic_cpdA", "synthetic_cpdB"]
        assert ps.is_vehicle.sum() == 3
        assert list(ps.doses("synthetic_cpdA")) == [1e-7, 1e-6, 1e-5]
        assert "PerturbationSet(" in repr(ps)

    def test_requires_obs_columns(self):
        a = ad.AnnData(np.zeros((2, 2)), obs=pd.DataFrame({"compound_id": ["a", "b"]}))
        with pytest.raises(TinctureDataError, match="missing"):
            PerturbationSet(a)

    def test_dose_must_be_molar(self):
        obs = pd.DataFrame({"compound_id": ["a"], "dose": [10.0], "replicate": ["1"]})
        with pytest.raises(TinctureDataError, match="molar"):
            PerturbationSet(ad.AnnData(np.zeros((1, 2)), obs=obs))

    def test_duplicate_sample_feature(self):
        d = pd.concat([synthetic_long_table()] * 2)
        with pytest.raises(TinctureDataError, match="share the same sample"):
            PerturbationSet.from_long(d)

    def test_subset(self):
        ps = PerturbationSet.from_long(synthetic_long_table())
        sub = ps.subset(compound_id="synthetic_cpdA", features=["g1", "g2"])
        assert sub.shape == (9, 2)
        assert ps.subset(compound_id=["synthetic_cpdA", "vehicle"]).shape[0] == 12
        with pytest.raises(TinctureDataError, match="selects no samples"):
            ps.subset(compound_id="nope")
        with pytest.raises(TinctureDataError, match="unknown feature"):
            ps.subset(features=["zzz"])
        with pytest.raises(TinctureDataError, match="Unknown obs column"):
            ps.subset(colour="red")

    def test_pseudobulk_sums_and_counts(self):
        ps = PerturbationSet.from_long(synthetic_long_table())
        pb = ps.pseudobulk(["compound_id", "dose"])
        X = np.asarray(ps.X)
        sel = (ps.obs["compound_id"] == "synthetic_cpdA") & (ps.obs["dose"] == 1e-6)
        row = pb.obs.index[(pb.obs["compound_id"] == "synthetic_cpdA") & (pb.obs["dose"] == 1e-6)]
        assert np.allclose(np.asarray(pb.adata[row].X).ravel(), X[sel.to_numpy()].sum(0))
        assert pb.obs.loc[row, "n_cells"].item() == 3

    def test_pseudobulk_sparse(self):
        ps = PerturbationSet.from_long(synthetic_long_table())
        import scipy.sparse as sp

        ps.adata.X = sp.csr_matrix(np.nan_to_num(np.asarray(ps.X)))
        assert ps.pseudobulk(["compound_id", "dose", "replicate"]).shape == ps.shape

    def test_single_replicate_is_fine(self):
        d = synthetic_long_table()
        ps = PerturbationSet.from_long(d[d["replicate"] == "r1"])
        assert ps.shape[0] == 7


# --------------------------------------------------------------------------- Signature


class TestSignature:
    def test_shapes_validated(self):
        with pytest.raises(TinctureDataError, match="shape"):
            Signature(features=["a", "b"], effect=[1.0])
        with pytest.raises(TinctureDataError, match="unique"):
            Signature(features=["a", "a"], effect=[1.0, 2.0])
        with pytest.raises(TinctureDataError, match="non-negative"):
            Signature(features=["a"], effect=[1.0], se=[-1.0])
        with pytest.raises(TinctureDataError, match="boot"):
            Signature(features=["a", "b"], effect=[1.0, 2.0], boot=np.zeros((10, 3)))

    def test_ci_from_boot_and_se(self):
        rng = np.random.default_rng(0)
        boot = rng.normal([0, 5], 1, size=(2000, 2))
        s = Signature(features=["a", "b"], effect=[0.0, 5.0], boot=boot)
        lo, hi = s.ci(0.95)
        assert lo[1] == pytest.approx(5 - 1.96, abs=0.15) and hi[0] == pytest.approx(1.96, abs=0.15)
        s2 = Signature(features=["a"], effect=[1.0], se=[0.5])
        assert s2.ci(0.95)[0][0] == pytest.approx(1 - 1.96 * 0.5, abs=1e-3)

    def test_no_uncertainty_says_so(self):
        with pytest.raises(TinctureDataError, match="no uncertainty"):
            Signature(features=["a"], effect=[1.0]).ci()

    def test_align_and_magnitude(self):
        s = Signature(features=["a", "b"], effect=[3.0, 4.0], se=[1.0, 1.0])
        assert s.magnitude == 5.0
        al = s.align(["b", "c"])
        assert al.effect[0] == 4.0 and np.isnan(al.effect[1]) and np.isnan(al.se[1])


# --------------------------------------------------------------------------- dose-response


class TestDoseResponse:
    curve = FourPL(bottom=0.1, top=1.0, log10_ec50=-6.0, hill=1.0)

    def test_fourpl_midpoint_and_inverse(self):
        assert self.curve(1e-6) == pytest.approx(0.55)
        assert self.curve.inverse(0.55) == pytest.approx(1e-6)
        assert np.isnan(self.curve.inverse(0.05))  # never reached

    def _dr(self, **kw):
        doses = np.repeat(np.logspace(-8, -5, 4), 2)
        return DoseResponse("x", doses, self.curve(doses), **kw)

    def test_flags(self):
        assert self._dr().quality_flags() == ["no fit"]
        assert self._dr(fit=self.curve, r2=0.95).fit_ok
        assert any("R²" in f for f in self._dr(fit=self.curve, r2=0.3).quality_flags())
        bad = FourPL(bottom=-40.0, top=1.0, log10_ec50=-6, hill=1)
        assert any("bottom" in f for f in self._dr(fit=bad, r2=0.9).quality_flags())
        far = FourPL(bottom=0.1, top=1.0, log10_ec50=2.0, hill=1)
        assert "EC50 outside tested range" in self._dr(fit=far, r2=0.9).quality_flags()
        assert (
            "did not converge" in self._dr(fit=self.curve, r2=0.9, converged=False).quality_flags()
        )

    def test_dose_units_checked(self):
        with pytest.raises(TinctureDataError, match="molar"):
            DoseResponse("x", [10.0], [1.0])
        with pytest.raises(TinctureDataError, match="equal length"):
            DoseResponse("x", [1e-6, 1e-5], [1.0])

    def test_set_get(self):
        pts = pd.DataFrame(
            {
                "compound_id": "x",
                "context_id": "L1",
                "dose": [1e-7, 1e-6],
                "replicate": "0",
                "response": [0.9, 0.5],
            }
        )
        fits = pd.DataFrame(
            {
                "compound_id": ["x"],
                "context_id": ["L1"],
                "bottom": [0.1],
                "top": [1.0],
                "log10_ec50": [-6.0],
                "hill": [1.0],
                "r2": [0.9],
                "auc": [0.7],
            }
        )
        s = DoseResponseSet(pts, fits, contexts=pd.DataFrame({"tissue": ["lung"]}, index=["L1"]))
        dr = s.get("x", "L1")
        assert dr.fit == self.curve and dr.context["auc"] == 0.7 and dr.context["tissue"] == "lung"
        assert len(s) == 1
        with pytest.raises(TinctureDataError, match="No points"):
            s.get("x", "L2")
        with pytest.raises(TinctureDataError, match="missing columns"):
            DoseResponseSet(pts.drop(columns=["replicate"]))


# --------------------------------------------------------------------------- plates


class TestPlates:
    @pytest.mark.parametrize("pos,rc", [(1, (0, 0)), (48, (0, 47)), (49, (1, 0)), (1536, (31, 47))])
    def test_row_major_1536(self, pos, rc):
        r, c = position_to_rowcol(np.array([pos]), 1536)
        assert (int(r[0]), int(c[0])) == rc

    def test_col_major_and_errors(self):
        r, c = position_to_rowcol(np.array([2]), 96, order="col-major")
        assert (int(r[0]), int(c[0])) == (1, 0)
        with pytest.raises(TinctureDataError, match="outside"):
            position_to_rowcol(np.array([97]), 96)
        with pytest.raises(TinctureDataError, match="Unknown plate format"):
            position_to_rowcol(np.array([1]), 100)

    def test_row_labels(self):
        assert [row_label(i) for i in (0, 25, 26, 31)] == ["A", "Z", "AA", "AF"]

    def _plate(self, n=96):
        nr, nc = {96: (8, 12), 1536: (32, 48)}[n]
        rr, cc = np.meshgrid(np.arange(nr), np.arange(nc), indexing="ij")
        d = pd.DataFrame(
            {
                "plate": "p1",
                "row": rr.ravel(),
                "col": cc.ravel(),
                "well_type": "compound",
                "value": 1.0,
            }
        )
        return PlateTable(d, n)

    def test_edges_and_matrix(self):
        pt = self._plate(1536)
        assert pt.is_edge().sum() == 2 * 48 + 2 * 30
        pt2 = PlateTable(pt.data.iloc[10:].reset_index(drop=True), 1536)
        m = pt2.plate_matrix("p1")
        assert np.isnan(m[0, :10]).all() and m.shape == (32, 48)  # absent wells stay NaN
        assert pt.well_name().iloc[0] == "A1"

    def test_validation(self):
        pt = self._plate()
        with pytest.raises(TinctureDataError, match="duplicated"):
            PlateTable(pd.concat([pt.data, pt.data.iloc[:1]]), 96)
        with pytest.raises(TinctureDataError, match="outside"):
            PlateTable(pt.data.assign(row=pt.data["row"] + 1), 96)
        with pytest.raises(TinctureDataError, match="well_type"):
            PlateTable(pt.data.assign(well_type="x"), 96)
        with pytest.raises(TinctureDataError, match="missing"):
            PlateTable(pt.data.drop(columns=["value"]), 96)
        with pytest.raises(TinctureDataError, match="plate"):
            pt.plate_matrix("nope")
