"""
Unit tests for the pure-numpy parts of tacs.postprocess.vmt.

These tests use synthetic nodal data with hand-derived expectations and do
not solve a TACS problem or read a real f5 file. Sign convention under test
(axis tangent t, shear direction s, bending axis b = s x t):

    V = F . s,  M = Mom . b,  T = Mom . t,  so that dM/ds = V.
"""

import os
import re
import struct
import tempfile
import unittest
import warnings

import numpy as np

from tacs.postprocess import vmt
from tacs.postprocess.vmt import (
    F5Data,
    VMTResult,
    buildSilhouette,
    computeEnvelope,
    computeVMT,
    plotVMTEnvelope,
    plotVMTReport,
    scaleVMTResult,
    selectComponents,
    writeVMTCsv,
    writeVMTEnvelopeCsv,
)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover
    matplotlib = None

P = 1000.0
STRAIGHT_AXIS = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]


def makeF5Data(X, ltypes, ptr, conn, comps=None, componentNames=("comp",)):
    """Build an F5Data object from raw arrays (coordinates only)."""
    X = np.asarray(X, dtype=np.float64)
    comps = np.zeros(len(ltypes), dtype=np.intp) if comps is None else np.asarray(comps)
    return F5Data(
        fileName="synthetic",
        varNames=("X", "Y", "Z"),
        data=X,
        comps=comps,
        ltypes=np.asarray(ltypes, dtype=np.intp),
        ptr=np.asarray(ptr, dtype=np.intp),
        conn=np.asarray(conn, dtype=np.intp),
        componentNames=tuple(componentNames),
    )


class TestComputeVMT(unittest.TestCase):
    N_PROCS = 1

    def setUp(self):
        warnings.simplefilter("ignore", UserWarning)

    def test_tip_force_straight_axis(self):
        result = computeVMT(
            [[10.0, 0.0, 0.0]], [[0.0, 0.0, P]], STRAIGHT_AXIS, numStations=11
        )
        s = np.linspace(0.0, 10.0, 11)
        np.testing.assert_allclose(result.s, s)
        np.testing.assert_allclose(result.V, np.where(s < 10.0, -P, 0.0))
        np.testing.assert_allclose(result.M, P * (10.0 - s))
        np.testing.assert_allclose(result.T, 0.0)
        np.testing.assert_allclose(result.force[:-1], [[0.0, 0.0, P]] * 10)
        np.testing.assert_allclose(result.moment[:-1, 1], -P * (10.0 - s[:-1]))
        self.assertEqual(result.axisLength, 10.0)

    def test_moment_derivative_equals_shear(self):
        rng = np.random.default_rng(0)
        X = np.column_stack(
            [np.arange(1, 10, dtype=float), rng.normal(size=9), rng.normal(size=9)]
        )
        F = rng.normal(size=(9, 3))
        M = rng.normal(size=(9, 3))
        # Stations strictly between two neighbouring nodes, so no node is crossed
        stations = np.linspace(4.1, 4.9, 5)
        result = computeVMT(X, F, STRAIGHT_AXIS, M=M, stations=stations)
        dMds = np.diff(result.M) / np.diff(result.s)
        np.testing.assert_allclose(dMds, result.V[:-1], rtol=1e-10, atol=1e-10)

    def test_offset_force_produces_torque(self):
        c = 2.0
        result = computeVMT(
            [[10.0, c, 0.0]], [[0.0, 0.0, P]], STRAIGHT_AXIS, numStations=3
        )
        np.testing.assert_allclose(result.V, [-P, -P, 0.0])
        np.testing.assert_allclose(result.M, [10.0 * P, 5.0 * P, 0.0])
        np.testing.assert_allclose(result.T, [c * P, c * P, 0.0])

    def test_pure_moment(self):
        result = computeVMT(
            [[10.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0]],
            STRAIGHT_AXIS,
            M=[[5.0, 7.0, 11.0]],
            numStations=2,
        )
        np.testing.assert_allclose(result.V, 0.0)
        np.testing.assert_allclose(result.T, [5.0, 0.0])
        np.testing.assert_allclose(result.M, [-7.0, 0.0])
        np.testing.assert_allclose(result.bendingAxis[0], [0.0, -1.0, 0.0])

    def test_kinked_axis(self):
        axis = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0]]
        result = computeVMT([[10.0, 10.0, 0.0]], [[0.0, 0.0, P]], axis, numStations=5)
        np.testing.assert_allclose(result.s, [0.0, 5.0, 10.0, 15.0, 20.0])
        np.testing.assert_allclose(result.V, [-P, -P, -P, -P, 0.0])
        # First segment (t = +x): bending (10 - s) P, torque 10 P from the offset
        np.testing.assert_allclose(result.M[:2], [10.0 * P, 5.0 * P])
        np.testing.assert_allclose(result.T[:2], [10.0 * P, 10.0 * P])
        # Station at the kink uses the outboard segment (t = +y)
        np.testing.assert_allclose(result.tangent[2], [0.0, 1.0, 0.0])
        np.testing.assert_allclose(result.M[2:], [10.0 * P, 5.0 * P, 0.0])
        np.testing.assert_allclose(result.T[2:], 0.0, atol=1e-12)
        np.testing.assert_allclose(result.bendingAxis[2], [1.0, 0.0, 0.0])

    def test_shear_dir_scaling_and_orthogonalisation(self):
        X, F = [[10.0, 0.0, 0.0]], [[0.0, 0.0, P]]
        ref = computeVMT(X, F, STRAIGHT_AXIS, numStations=4)
        scaled = computeVMT(
            X, F, STRAIGHT_AXIS, shearDir=(0.0, 0.0, -5.0), numStations=4
        )
        np.testing.assert_allclose(scaled.V, ref.V)
        np.testing.assert_allclose(scaled.M, ref.M)
        skew = computeVMT(X, F, STRAIGHT_AXIS, shearDir=(0.3, 0.0, -1.0), numStations=4)
        np.testing.assert_allclose(np.linalg.norm(skew.shearAxis, axis=1), 1.0)
        np.testing.assert_allclose(
            np.einsum("ij,ij->i", skew.shearAxis, skew.tangent), 0.0, atol=1e-14
        )
        np.testing.assert_allclose(skew.V, ref.V)

    def test_shear_dir_parallel_to_axis_raises(self):
        with self.assertRaises(ValueError) as ctx:
            computeVMT(
                [[10.0, 0.0, 0.0]],
                [[0.0, 0.0, P]],
                STRAIGHT_AXIS,
                shearDir=(1.0, 0.0, 0.0),
            )
        self.assertIn("segment 0", str(ctx.exception))

    def test_projection_clamps_and_kink_tie(self):
        axis = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0]]
        X = [[-3.0, 1.0, 0.0], [10.0, 12.0, 0.5], [12.0, -2.0, 0.0], [5.0, 3.0, 0.0]]
        result = computeVMT(X, np.zeros((4, 3)), axis, numStations=2)
        self.assertEqual(result.nodeS[0], 0.0)
        self.assertEqual(result.nodeS[1], 20.0)
        # Wedge outside the 90 degree kink: equidistant from both segments
        self.assertEqual(result.nodeS[2], 10.0)
        self.assertAlmostEqual(result.nodeS[3], 5.0)

    def test_station_tolerance(self):
        L = 10.0
        tol = 1e-6 * L
        X = [[5.0, 0.0, 0.0], [5.0 + 2.0 * tol, 0.0, 0.0], [5.0 + 0.5 * tol, 0.0, 0.0]]
        F = [[0.0, 0.0, 1.0], [0.0, 0.0, 10.0], [0.0, 0.0, 100.0]]
        result = computeVMT(X, F, STRAIGHT_AXIS, stations=[5.0])
        np.testing.assert_allclose(result.V, [-10.0])
        # float32 round trip of a node ring at the root is still excluded
        X32 = np.array([[0.0, 0.3, 0.7], [10.0, 0.0, 0.0]], dtype=np.float32).astype(
            np.float64
        )
        result = computeVMT(
            X32, [[0.0, 0.0, -P], [0.0, 0.0, P]], STRAIGHT_AXIS, numStations=3
        )
        np.testing.assert_allclose(result.V, [-P, -P, 0.0])

    def test_root_reaction_rule(self):
        # Tip load plus the equilibrating root reaction force and moment
        X = [[10.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        F = [[0.0, 0.0, P], [0.0, 0.0, -P]]
        M = [[0.0, 0.0, 0.0], [0.0, 10.0 * P, 0.0]]
        result = computeVMT(X, F, STRAIGHT_AXIS, M=M, numStations=3)
        np.testing.assert_allclose(result.V, [-P, -P, 0.0])
        np.testing.assert_allclose(result.M, [10.0 * P, 5.0 * P, 0.0])
        extended = computeVMT(
            X, F, [[-2.0, 0.0, 0.0], [10.0, 0.0, 0.0]], M=M, numStations=3
        )
        np.testing.assert_allclose(extended.V, [0.0, -P, 0.0])
        np.testing.assert_allclose(extended.M[0], 0.0, atol=1e-12)

    def test_midpoint_stations(self):
        result = computeVMT(
            [[10.0, 0.0, 0.0]],
            [[0.0, 0.0, P]],
            STRAIGHT_AXIS,
            numStations=4,
            stationMode="midpoint",
        )
        np.testing.assert_allclose(result.s, [1.25, 3.75, 6.25, 8.75])
        np.testing.assert_allclose(result.V, -P)

    def test_node_mask(self):
        X = [[10.0, 0.0, 0.0], [5.0, 0.0, 0.0]]
        F = [[0.0, 0.0, P], [0.0, 0.0, 7.0 * P]]
        result = computeVMT(X, F, STRAIGHT_AXIS, numStations=2, nodeMask=[True, False])
        np.testing.assert_allclose(result.V, [-P, 0.0])
        self.assertEqual(result.nodeS.shape, (1,))

    def test_scale_result(self):
        result = computeVMT(
            [[10.0, 2.0, 0.0]], [[0.0, 0.0, P]], STRAIGHT_AXIS, numStations=11
        )
        self.assertIs(scaleVMTResult(result), result)
        scaled = scaleVMTResult(result, forceScale=2.0, lengthScale=3.0)
        np.testing.assert_allclose(scaled.s, 3.0 * result.s)
        np.testing.assert_allclose(scaled.nodeS, 3.0 * result.nodeS)
        np.testing.assert_allclose(scaled.V, 2.0 * result.V)
        np.testing.assert_allclose(scaled.M, 6.0 * result.M)
        np.testing.assert_allclose(scaled.T, 6.0 * result.T)
        np.testing.assert_allclose(scaled.force, 2.0 * result.force)
        np.testing.assert_allclose(scaled.moment, 6.0 * result.moment)
        np.testing.assert_allclose(scaled.stationPoints, 3.0 * result.stationPoints)
        np.testing.assert_allclose(scaled.axisPoints, 3.0 * result.axisPoints)
        self.assertAlmostEqual(scaled.axisLength, 30.0)
        # The unit frame is not a length and must be untouched
        np.testing.assert_allclose(scaled.tangent, result.tangent)
        np.testing.assert_allclose(scaled.shearAxis, result.shearAxis)
        # dM/ds = V still holds after scaling
        np.testing.assert_allclose(
            np.diff(scaled.M) / np.diff(scaled.s), scaled.V[:-1], atol=1e-9
        )
        for bad in (0.0, -1.0, np.inf, np.nan, "abc"):
            with self.assertRaises(ValueError):
                scaleVMTResult(result, forceScale=bad)
            with self.assertRaises(ValueError):
                scaleVMTResult(result, lengthScale=bad)

    def test_invalid_input(self):
        with self.assertRaises(ValueError):
            computeVMT([[1.0, 0.0, 0.0]], [[0.0, 0.0, 1.0]], [[0.0, 0.0, 0.0]])
        with self.assertRaises(ValueError):
            computeVMT(
                [[1.0, 0.0, 0.0]], [[0.0, 0.0, 1.0]], [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]
            )
        with self.assertRaises(ValueError):
            computeVMT(
                [[1.0, 0.0, 0.0]], [[0.0, 0.0, 1.0]], STRAIGHT_AXIS, numStations=0
            )
        with self.assertRaises(ValueError):
            computeVMT(
                [[1.0, 0.0, 0.0]], [[0.0, 0.0, 1.0]], STRAIGHT_AXIS, stations=[11.0]
            )
        with self.assertRaises(ValueError):
            computeVMT([[1.0, 0.0, 0.0]], [[0.0, 0.0, np.nan]], STRAIGHT_AXIS)
        with self.assertRaises(ValueError):
            computeVMT(
                [[1.0, 0.0, 0.0]],
                [[0.0, 0.0, 1.0]],
                STRAIGHT_AXIS,
                shearDir=(0.0, 0.0, 0.0),
            )
        with self.assertRaises(ValueError):
            computeVMT(
                [[1.0, 0.0, 0.0]], [[0.0, 0.0, 1.0]], STRAIGHT_AXIS, stationMode="nope"
            )


class TestSelectionAndSilhouette(unittest.TestCase):
    N_PROCS = 1

    def test_quad4_polygon_order(self):
        # TACS tensor-product ordering: (0,0), (1,0), (0,1), (1,1)
        X = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
        data = makeF5Data(X, [vmt.QUAD_ELEMENT], [0, 4], [0, 1, 2, 3])
        result = computeVMT(
            X, np.zeros((4, 3)), [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], numStations=2
        )
        sil = buildSilhouette(data, result)
        self.assertEqual(sil.quads.shape, (1, 4, 2))
        xy = sil.quads[0]
        area = 0.5 * np.sum(
            xy[:, 0] * np.roll(xy[:, 1], -1) - np.roll(xy[:, 0], -1) * xy[:, 1]
        )
        self.assertAlmostEqual(abs(area), 1.0)

    def test_quad9_corners_and_other_layouts(self):
        # 3x3 tensor grid on the unit square plus a line element and an RBE2
        grid = (
            np.array([[i, j, 0.0] for j in range(3) for i in range(3)], dtype=float)
            / 2.0
        )
        X = np.vstack([grid, [[2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [9.0, 9.0, 9.0]]])
        ltypes = [
            vmt.QUAD_QUADRATIC_ELEMENT,
            vmt.LINE_ELEMENT,
            vmt.RBE2_ELEMENT,
            vmt.POINT_ELEMENT,
        ]
        conn = list(range(9)) + [9, 10] + [9, 10, 11] + [8]
        ptr = [0, 9, 11, 14, 15]
        data = makeF5Data(
            X, ltypes, ptr, conn, comps=[0, 1, 1, 1], componentNames=("skin", "other")
        )
        result = computeVMT(
            X, np.zeros((12, 3)), [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], numStations=2
        )
        sil = buildSilhouette(data, result)
        self.assertEqual(sil.quads.shape, (1, 4, 2))
        np.testing.assert_allclose(
            sil.quads[0], [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
        )
        self.assertEqual(sil.lines.shape, (1, 2, 2))
        self.assertEqual(sil.tris.shape[0], 0)
        # Node selection: RBE2 [1 indep | 1 dep | 1 multiplier] drops node 11
        _, nodeMask = selectComponents(data)
        self.assertFalse(nodeMask[11])
        self.assertTrue(nodeMask[[9, 10, 8]].all())
        elemMask, nodeMask = selectComponents(data, "SK*")
        np.testing.assert_array_equal(elemMask, [True, False, False, False])
        self.assertEqual(nodeMask.sum(), 9)
        elemMask, _ = selectComponents(data, [1])
        np.testing.assert_array_equal(elemMask, [False, True, True, True])
        with self.assertRaises(ValueError):
            selectComponents(data, ["nope"])
        with self.assertRaises(ValueError):
            selectComponents(data, [7])

    def test_csv_round_trip(self):
        result = computeVMT(
            [[10.0, 1.0, 0.0]], [[0.0, 0.0, P]], STRAIGHT_AXIS, numStations=5
        )
        with tempfile.TemporaryDirectory() as tmp:
            fname = os.path.join(tmp, "vmt.csv")
            writeVMTCsv(result, fname, metadata={"case": "unit"})
            with open(fname) as fp:
                header = (
                    [line for line in fp if not line.startswith("#")][0]
                    .strip()
                    .split(",")
                )
            table = np.loadtxt(fname, delimiter=",", comments="#", skiprows=1 + 5)
        self.assertEqual(header[:4], ["s", "V", "M", "T"])
        np.testing.assert_allclose(table[:, 0], result.s)
        np.testing.assert_allclose(table[:, 1], result.V)
        np.testing.assert_allclose(table[:, 2], result.M)
        np.testing.assert_allclose(table[:, 3], result.T)
        self.assertIsInstance(result, VMTResult)


class TestEnvelope(unittest.TestCase):
    N_PROCS = 1

    def setUp(self):
        warnings.simplefilter("ignore", UserWarning)
        X = [[10.0, 0.0, 0.0], [10.0, 2.0, 0.0]]
        self.up = computeVMT(
            X, [[0.0, 0.0, P], [0.0, 0.0, 0.0]], STRAIGHT_AXIS, numStations=6
        )
        self.down = computeVMT(
            X, [[0.0, 0.0, -2.0 * P], [0.0, 0.0, 0.0]], STRAIGHT_AXIS, numStations=6
        )
        self.twist = computeVMT(
            X, [[0.0, 0.0, 0.0], [0.0, 0.0, P]], STRAIGHT_AXIS, numStations=6
        )
        self.cases = {"up": self.up, "down": self.down, "twist": self.twist}

    def test_envelope_values_and_governing_case(self):
        # With shearDir = -z the downward (-2P) case has the largest positive
        # shear and the most negative bending moment; "up" and "twist" carry
        # the same P load, but only "twist" applies it off the axis
        env = computeEnvelope(self.cases)
        self.assertEqual(env.caseNames, ("up", "down", "twist"))
        np.testing.assert_allclose(env.s, self.up.s)
        np.testing.assert_allclose(env.Vmax, self.down.V)
        np.testing.assert_allclose(env.Vmin, self.up.V)
        np.testing.assert_allclose(env.Mmin, self.down.M)
        np.testing.assert_allclose(env.Mmax, self.up.M)
        np.testing.assert_allclose(env.Tmax, self.twist.T)
        np.testing.assert_allclose(env.Tmin, 0.0)
        self.assertEqual(list(env.VmaxCase[:-1]), ["down"] * 5)
        self.assertEqual(list(env.MminCase[:-1]), ["down"] * 5)
        self.assertEqual(list(env.TmaxCase[:-1]), ["twist"] * 5)
        # Ties resolve to the first case in order
        self.assertEqual(list(env.VminCase[:-1]), ["up"] * 5)

    def test_envelope_from_sequence_and_mismatch(self):
        env = computeEnvelope([self.up, self.down])
        self.assertEqual(env.caseNames, ("case 1", "case 2"))
        other = computeVMT(
            [[10.0, 0.0, 0.0]], [[0.0, 0.0, P]], STRAIGHT_AXIS, numStations=5
        )
        with self.assertRaises(ValueError):
            computeEnvelope({"a": self.up, "b": other})
        with self.assertRaises(ValueError):
            computeEnvelope({})

    def test_envelope_csv(self):
        env = computeEnvelope(self.cases)
        with tempfile.TemporaryDirectory() as tmp:
            fname = os.path.join(tmp, "env.csv")
            writeVMTEnvelopeCsv(env, fname, metadata={"note": "unit"})
            with open(fname) as fp:
                numMeta = sum(1 for line in fp if line.startswith("#"))
            table = np.genfromtxt(
                fname,
                delimiter=",",
                names=True,
                dtype=None,
                skip_header=numMeta,
                encoding=None,
            )
        np.testing.assert_allclose(table["Vmin"], env.Vmin)
        np.testing.assert_allclose(table["Tmax"], env.Tmax)
        self.assertEqual(list(table["VmaxCase"][:-1]), ["down"] * 5)

    @unittest.skipIf(matplotlib is None, "matplotlib not installed")
    @unittest.skipIf(matplotlib is None, "matplotlib not installed")
    def test_plot_units_and_scale(self):
        fig, axes = plotVMTEnvelope(
            self.cases,
            forceUnit="lbf",
            lengthUnit="in",
            forceScale=0.5,
            lengthScale=2.0,
        )
        try:
            axPlan, axV, axM, axT = axes
            self.assertTrue(axV.get_ylabel().endswith("[lbf]"))
            self.assertTrue(axM.get_ylabel().endswith("[lbf·in]"))
            self.assertTrue(axT.get_ylabel().endswith("[lbf·in]"))
            self.assertTrue(axT.get_xlabel().endswith("[in]"))
            self.assertTrue(axPlan.get_xlabel().endswith("[in]"))
            self.assertTrue(axPlan.get_ylabel().endswith("[in]"))
            self.assertEqual(axV.get_xlim(), (0.0, 20.0))
            env = computeEnvelope(self.cases)
            # Case lines are drawn after the two envelope curves
            lineV = axV.get_lines()[-3]
            np.testing.assert_allclose(lineV.get_xdata(), 2.0 * self.up.s)
            np.testing.assert_allclose(lineV.get_ydata(), 0.5 * self.up.V)
            lineM = axM.get_lines()[-3]
            np.testing.assert_allclose(lineM.get_ydata(), self.up.M)
            from matplotlib.collections import PolyCollection

            band = [c for c in axV.collections if isinstance(c, PolyCollection)][0]
            verts = band.get_paths()[0].vertices
            self.assertAlmostEqual(verts[:, 1].min(), 0.5 * env.Vmin.min())
            self.assertAlmostEqual(verts[:, 1].max(), 0.5 * env.Vmax.max())
        finally:
            plt.close(fig)
        # Only a force unit: moments get no unit rather than a half label
        fig, axes = plotVMTEnvelope(self.cases, forceUnit="N")
        try:
            self.assertTrue(axes[1].get_ylabel().endswith("[N]"))
            self.assertFalse(axes[2].get_ylabel().endswith("]"))
            self.assertFalse(axes[3].get_xlabel().endswith("]"))
        finally:
            plt.close(fig)
        with self.assertRaises(ValueError):
            plotVMTEnvelope(self.cases, lengthScale=0.0)

    def test_envelope_plot_and_report(self):
        fig, axes = plotVMTEnvelope(self.cases)
        try:
            self.assertEqual(len(axes), 4)
            # Envelope band plus one line per case in the shear panel legend
            labels = [t.get_text() for t in axes[1].get_legend().get_texts()]
            self.assertEqual(labels, ["envelope (min/max)", "up", "down", "twist"])
        finally:
            plt.close(fig)
        with tempfile.TemporaryDirectory() as tmp:
            fname = os.path.join(tmp, "report.pdf")
            pages = plotVMTReport(self.cases, fname, title="synthetic")
            self.assertEqual(pages, 4)
            with open(fname, "rb") as fp:
                content = fp.read()
            self.assertEqual(len(re.findall(rb"/Type\s*/Page\b", content)), 4)
            single = plotVMTReport({"up": self.up}, os.path.join(tmp, "one.pdf"))
            self.assertEqual(single, 1)


class TestF5HeaderScan(unittest.TestCase):
    N_PROCS = 1

    @staticmethod
    def _writeZone(fp, name, varNames, dtype, array):
        zname = name.encode() + b"\0"
        vname = varNames.encode() + b"\0"
        fp.write(
            struct.pack(
                "5i", dtype, array.shape[0], array.shape[1], len(zname), len(vname)
            )
        )
        fp.write(zname)
        fp.write(vname)
        fp.write(array.tobytes())

    def _writeMinimalF5(self, fname, truncate=False, dropZone=None):
        with open(fname, "wb") as fp:
            fp.write(struct.pack("i", 1))
            fp.write(struct.pack("i", 5) + b"comp\0")
            ints = np.array([[0]], dtype=np.int32)
            for zone in ("components", "ltypes", "ptr", "connectivity"):
                if zone != dropZone:
                    self._writeZone(fp, zone, "", 0, ints)
            self._writeZone(
                fp,
                "continuous data t=0",
                "X,Y,Z",
                2,
                np.zeros((3, 3), dtype=np.float32),
            )
        if truncate:
            size = os.path.getsize(fname)
            with open(fname, "r+b") as fp:
                fp.truncate(size - 8)

    def test_scan_and_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            good = os.path.join(tmp, "good.f5")
            self._writeMinimalF5(good)
            names, zones = vmt._scanF5Header(good)
            self.assertEqual(names, ["comp"])
            self.assertEqual([z.name for z in zones][-1], "continuous data t=0")
            self.assertEqual(zones[-1].varNames, "X,Y,Z")

            junk = os.path.join(tmp, "junk.f5")
            with open(junk, "wb") as fp:
                fp.write(b"\x07" * 100)
            with self.assertRaises(ValueError):
                vmt._scanF5Header(junk)

            truncated = os.path.join(tmp, "trunc.f5")
            self._writeMinimalF5(truncated, truncate=True)
            with self.assertRaises(ValueError):
                vmt._scanF5Header(truncated)

            noconn = os.path.join(tmp, "noconn.f5")
            self._writeMinimalF5(noconn, dropZone="connectivity")
            with self.assertRaises(ValueError) as ctx:
                vmt.loadF5(noconn)
            self.assertIn("writeConnectivity", str(ctx.exception))

            with self.assertRaises(FileNotFoundError):
                vmt.loadF5(os.path.join(tmp, "missing.f5"))


if __name__ == "__main__":
    unittest.main()
