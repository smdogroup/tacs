"""
Round-trip tests for tacs.postprocess.vmt on f5 files written by pyTACS.

The beam model gives exact cantilever values, the clamped plate exercises
shell elements, pressure loads and reactions, and the coarse wingbox
exercises component filtering, RBE3/CONM2 handling, plotting and the CLI.
"""

import os
import re
import tempfile
import unittest
import warnings

import numpy as np
from mpi4py import MPI

from tacs import constitutive, elements, pytacs
from tacs.postprocess import (
    buildSilhouette,
    computeEnvelope,
    computeVMTCases,
    computeVMTFromF5,
    listComponents,
    loadF5,
    nodalLoads,
    plotVMT,
    selectComponents,
)
from tacs.postprocess.vmt import main

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover
    matplotlib = None

base_dir = os.path.dirname(os.path.abspath(__file__))
beam_bdf = os.path.join(base_dir, "./input_files/beam_model.bdf")
plate_bdf = os.path.join(base_dir, "./input_files/plate.bdf")
wingbox_bdf = os.path.join(base_dir, "./input_files/coarse_mdo_tutorial_wingbox.bdf")

QUIET = {"printLevel": 0}


def solveAndWrite(problems, outputDir):
    """Solve every problem tightly and write its f5 file."""
    for problem in problems:
        problem.setOption("printLevel", 0)
        problem.setOption("L2Convergence", 1e-20)
        problem.setOption("L2ConvergenceRel", 1e-20)
        problem.solve()
        problem.writeSolution(outputDir=outputDir)


def f5Path(outputDir, name):
    return os.path.join(outputDir, f"{name}_000.f5")


def makePlateProblems(comm):
    """Replicate the plate problems of test_shell_plate_quad.py."""
    fea_assembler = pytacs.pyTACS(plate_bdf, comm, options=QUIET)

    def elem_call_back(
        dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
    ):
        prop = constitutive.MaterialProperties(rho=2500.0, E=70e9, nu=0.3, ys=464.0e6)
        con = constitutive.IsoShellConstitutive(prop, t=0.005, tNum=dv_num)
        return elements.Quad4Shell(None, con), [100.0]

    fea_assembler.initialize(elem_call_back)
    point = fea_assembler.createStaticProblem(name="point_load")
    point.addLoadToNodes(
        81, np.array([0.0, 0.0, 1e4, 0.0, 0.0, 0.0]), nastranOrdering=True
    )
    pressure = fea_assembler.createStaticProblem(name="pressure")
    pressure.addPressureToComponents(
        fea_assembler.selectCompIDs(include="PLATE"), 100e3
    )
    return fea_assembler, [point, pressure]


class TestBeamModelVMT(unittest.TestCase):
    N_PROCS = 1

    def setUp(self):
        # testflo runs each test method in its own process and never calls
        # setUpClass, so the (cheap) solve is repeated per test
        self.tmp = tempfile.TemporaryDirectory()
        fea_assembler = pytacs.pyTACS(beam_bdf, MPI.COMM_WORLD, options=QUIET)
        fea_assembler.initialize()
        problems = fea_assembler.createTACSProbsFromBDF()
        solveAndWrite(problems.values(), self.tmp.name)

        self.axis = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        self.tol = dict(rtol=1e-6, atol=1e-7)

    def tearDown(self):
        self.tmp.cleanup()

    def test_z_shear(self):
        result = computeVMTFromF5(
            f5Path(self.tmp.name, "z-shear"), self.axis, numStations=6
        )
        np.testing.assert_allclose(result.V, [-1.0] * 5 + [0.0], **self.tol)
        np.testing.assert_allclose(result.M, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0], **self.tol)
        np.testing.assert_allclose(result.T, 0.0, **self.tol)
        self.assertTrue(result.includeReactions)

    def test_x_torsion(self):
        result = computeVMTFromF5(
            f5Path(self.tmp.name, "x-torsion"), self.axis, numStations=6
        )
        np.testing.assert_allclose(result.T, [1.0] * 5 + [0.0], **self.tol)
        np.testing.assert_allclose(result.V, 0.0, **self.tol)
        np.testing.assert_allclose(result.M, 0.0, **self.tol)

    def test_y_shear_and_shear_dir_override(self):
        fname = f5Path(self.tmp.name, "y-shear")
        result = computeVMTFromF5(fname, self.axis, numStations=6)
        np.testing.assert_allclose(result.V, 0.0, **self.tol)
        np.testing.assert_allclose(result.M, 0.0, **self.tol)
        np.testing.assert_allclose(result.T, 0.0, **self.tol)
        np.testing.assert_allclose(result.force[0], [0.0, 1.0, 0.0], **self.tol)
        result = computeVMTFromF5(
            fname, self.axis, shearDir=(0.0, -1.0, 0.0), numStations=6
        )
        np.testing.assert_allclose(result.V, [-1.0] * 5 + [0.0], **self.tol)
        np.testing.assert_allclose(result.M, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0], **self.tol)

    def test_x_axial(self):
        result = computeVMTFromF5(
            f5Path(self.tmp.name, "x-axial"), self.axis, numStations=6
        )
        np.testing.assert_allclose(result.V, 0.0, **self.tol)
        np.testing.assert_allclose(result.M, 0.0, **self.tol)
        np.testing.assert_allclose(result.T, 0.0, **self.tol)
        np.testing.assert_allclose(result.force[0], [1.0, 0.0, 0.0], **self.tol)

    def test_reactions_and_equilibrium(self):
        fname = f5Path(self.tmp.name, "z-shear")
        withR = computeVMTFromF5(fname, self.axis, numStations=6, includeReactions=True)
        without = computeVMTFromF5(
            fname, self.axis, numStations=6, includeReactions=False
        )
        np.testing.assert_allclose(withR.V, without.V, **self.tol)
        np.testing.assert_allclose(withR.M, without.M, **self.tol)
        extended = computeVMTFromF5(
            fname, [[-0.5, 0.0, 0.0], [1.0, 0.0, 0.0]], numStations=4
        )
        np.testing.assert_allclose(
            [extended.V[0], extended.M[0], extended.T[0]], 0.0, atol=1e-9
        )
        np.testing.assert_allclose(extended.V[1:], [-1.0, -1.0, 0.0], **self.tol)

    def test_line_elements_in_silhouette(self):
        if matplotlib is None:
            self.skipTest("matplotlib not installed")
        data = loadF5(f5Path(self.tmp.name, "z-shear"))
        result = computeVMTFromF5(data, self.axis, numStations=6)
        fig, axes = plotVMT(result, data)
        self.assertEqual(len(axes), 4)
        plt.close(fig)


class TestPlateVMT(unittest.TestCase):
    N_PROCS = 1

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        _, problems = makePlateProblems(MPI.COMM_WORLD)
        solveAndWrite(problems, self.tmp.name)

        self.axis = [[0.0, 0.5, 0.0], [1.0, 0.5, 0.0]]
        self.tol = dict(rtol=1e-5, atol=1e-1)

    def tearDown(self):
        self.tmp.cleanup()

    def test_point_load(self):
        # The plate is clamped on all four edges, so reactions act outboard of
        # every station; sum the applied loads only to get a hand-checkable result
        result = computeVMTFromF5(
            f5Path(self.tmp.name, "point_load"),
            self.axis,
            numStations=11,
            includeReactions=False,
        )
        # Node 81 sits at x = 0.5, exactly on the sixth station, so it counts
        # as inboard there
        np.testing.assert_allclose(result.V, [-1e4] * 5 + [0.0] * 6, **self.tol)
        np.testing.assert_allclose(
            result.M, [5000.0, 4000.0, 3000.0, 2000.0, 1000.0] + [0.0] * 6, **self.tol
        )
        np.testing.assert_allclose(result.T, 0.0, atol=1e-1)

    def test_pressure(self):
        result = computeVMTFromF5(
            f5Path(self.tmp.name, "pressure"),
            self.axis,
            numStations=11,
            includeReactions=False,
        )
        # 100 kPa over 0.1 m x 0.1 m elements = 1000 N per element, 250 N per
        # corner. Every perimeter node is clamped, so its load row is zeroed
        # and only the 9 x 9 interior nodes carry 1000 N each: the root shear
        # is 81 kN, not the 100 kN total pressure load.
        np.testing.assert_allclose(result.V[0], -81000.0, **self.tol)
        np.testing.assert_allclose(result.V[5], -36000.0, **self.tol)
        np.testing.assert_allclose(result.M[5], 9000.0, **self.tol)
        np.testing.assert_allclose(result.V[-1], 0.0, atol=1e-1)
        np.testing.assert_allclose(result.T, 0.0, atol=1e-1)

    def test_reactions_close_equilibrium(self):
        axis = [[-0.25, 0.5, 0.0], [1.25, 0.5, 0.0]]
        for name, rootShear in (("point_load", -1e4), ("pressure", -81000.0)):
            fname = f5Path(self.tmp.name, name)
            withR = computeVMTFromF5(fname, axis, numStations=4, includeReactions=True)
            np.testing.assert_allclose(
                [withR.V[0], withR.M[0], withR.T[0]], 0.0, atol=1e-1
            )
            without = computeVMTFromF5(
                fname, axis, numStations=4, includeReactions=False
            )
            np.testing.assert_allclose(without.V[0], rootShear, **self.tol)

    def test_component_filter(self):
        data = loadF5(f5Path(self.tmp.name, "point_load"))
        names = [name for name, _ in listComponents(data)]
        self.assertTrue(names)
        full = computeVMTFromF5(data, self.axis, numStations=11)
        filtered = computeVMTFromF5(data, self.axis, numStations=11, components=names)
        np.testing.assert_allclose(filtered.V, full.V)
        with self.assertRaises(ValueError):
            computeVMTFromF5(data, self.axis, components=["NOPE"])


class TestPlateVMTParallel(unittest.TestCase):
    N_PROCS = 2

    def test_partition_invariance(self):
        comm = MPI.COMM_WORLD
        tmp = tempfile.TemporaryDirectory() if comm.rank == 0 else None
        tmpName = comm.bcast(tmp.name if tmp is not None else None, root=0)
        _, problems = makePlateProblems(comm)
        solveAndWrite(problems[:1], tmpName)
        comm.Barrier()
        error = None
        if comm.rank == 0:
            try:
                result = computeVMTFromF5(
                    f5Path(tmpName, "point_load"),
                    [[0.0, 0.5, 0.0], [1.0, 0.5, 0.0]],
                    numStations=11,
                    includeReactions=False,
                )
                np.testing.assert_allclose(
                    result.V, [-1e4] * 5 + [0.0] * 6, rtol=1e-5, atol=1e-1
                )
                np.testing.assert_allclose(
                    result.M,
                    [5000.0, 4000.0, 3000.0, 2000.0, 1000.0] + [0.0] * 6,
                    rtol=1e-5,
                    atol=1e-1,
                )
            except Exception as exc:  # noqa: BLE001 - re-raised on every rank below
                error = f"{type(exc).__name__}: {exc}"
            finally:
                tmp.cleanup()
        error = comm.bcast(error, root=0)
        self.assertIsNone(error, error)


class TestWingboxVMT(unittest.TestCase):
    N_PROCS = 1

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fea_assembler = pytacs.pyTACS(wingbox_bdf, MPI.COMM_WORLD, options=QUIET)
            fea_assembler.initialize()
        pullup = fea_assembler.createStaticProblem("grav")
        pullup.addInertialLoad([0.0, -9.81 * 2.5, 0.0])
        pushover = fea_assembler.createStaticProblem("pushover")
        pushover.addInertialLoad([0.0, 9.81, 0.0])
        solveAndWrite([pullup, pushover], self.tmp.name)
        self.f5 = f5Path(self.tmp.name, "grav")
        self.f5Pushover = f5Path(self.tmp.name, "pushover")
        self.data = loadF5(self.f5)

        # Span runs along z, the load acts in -y
        self.axis = [[4.0, 0.0, 0.0], [4.0, 0.0, 13.8]]
        self.shearDir = (0.0, -1.0, 0.0)

    def tearDown(self):
        self.tmp.cleanup()

    def test_root_shear_matches_total_load(self):
        result = computeVMTFromF5(
            self.data, self.axis, shearDir=self.shearDir, numStations=20
        )
        self.assertTrue(np.all(np.isfinite(result.V)))
        self.assertTrue(np.all(np.isfinite(result.M)))
        self.assertTrue(np.all(np.isfinite(result.T)))
        F, _, _ = nodalLoads(self.data, includeReactions=False)
        np.testing.assert_allclose(result.V[0], -F[:, 1].sum(), rtol=1e-6)
        np.testing.assert_allclose(result.V[-1], 0.0)
        extended = computeVMTFromF5(
            self.data,
            [[4.0, 0.0, -1.0], [4.0, 0.0, 13.8]],
            shearDir=self.shearDir,
            numStations=5,
        )
        np.testing.assert_allclose(extended.V[0], 0.0, atol=1e-6 * abs(result.V[0]))

    def test_component_globs_and_rbe(self):
        names = listComponents(self.data)
        self.assertGreaterEqual(len(names), 90)
        full = computeVMTFromF5(
            self.data, self.axis, shearDir=self.shearDir, numStations=20
        )
        # Every structural node of this coarse box lies on a skin, so restrict
        # to the upper skin only: the lower-skin nodes are then excluded and
        # the outboard load must strictly decrease
        upper = computeVMTFromF5(
            self.data,
            self.axis,
            shearDir=self.shearDir,
            numStations=20,
            components=["WING_U_SKIN*"],
        )
        self.assertLess(abs(upper.V[0]), 0.9 * abs(full.V[0]))
        self.assertGreater(abs(upper.V[0]), 0.1 * abs(full.V[0]))
        # Whether the RBE3 / CONM2 pair of this model reaches the f5 file
        # depends on the pyTACS version; when it does, its trailing
        # Lagrange-multiplier node must be excluded from the load sum
        _, nodeMask = selectComponents(self.data)
        rbe = np.flatnonzero(self.data.ltypes == 25)
        if rbe.size:
            multiplier = self.data.conn[self.data.ptr[rbe[0] + 1] - 1]
            self.assertFalse(nodeMask[multiplier])
            self.assertEqual(nodeMask.sum(), self.data.numNodes - rbe.size)
        else:
            self.assertEqual(nodeMask.sum(), self.data.numNodes)

    @unittest.skipIf(matplotlib is None, "matplotlib not installed")
    def test_plot(self):
        result = computeVMTFromF5(
            self.data, self.axis, shearDir=self.shearDir, numStations=20
        )
        png = os.path.join(self.tmp.name, "vmt.png")
        fig, axes = plotVMT(result, self.data, shearDir=self.shearDir, fileName=png)
        try:
            self.assertEqual(len(axes), 4)
            self.assertTrue(os.path.getsize(png) > 0)
            polys = [
                c
                for c in axes[0].collections
                if c.__class__.__name__ == "PolyCollection"
            ]
            self.assertEqual(sum(len(c.get_paths()) for c in polys), 91)
        finally:
            plt.close(fig)

    @unittest.skipIf(matplotlib is None, "matplotlib not installed")
    def test_plot_units_and_scale(self):
        result = computeVMTFromF5(
            self.data, self.axis, shearDir=self.shearDir, numStations=20
        )
        silhouette = buildSilhouette(self.data, result, self.shearDir)
        fig, axes = plotVMT(
            result,
            silhouette=silhouette,
            shearDir=self.shearDir,
            forceUnit="N",
            lengthUnit="m",
            forceScale=1e-3,
            lengthScale=1e3,
        )
        try:
            axPlan, axV, axM, axT = axes
            self.assertTrue(axV.get_ylabel().endswith("[N]"))
            self.assertTrue(axM.get_ylabel().endswith("[N·m]"))
            self.assertTrue(axT.get_xlabel().endswith("[m]"))
            self.assertAlmostEqual(axV.get_xlim()[1], 1e3 * result.axisLength)
            np.testing.assert_allclose(axV.get_lines()[-1].get_ydata(), 1e-3 * result.V)
            np.testing.assert_allclose(axM.get_lines()[-1].get_ydata(), result.M)
            np.testing.assert_allclose(axM.get_lines()[-1].get_xdata(), 1e3 * result.s)
            # The planform silhouette is drawn in the scaled length unit
            polys = [
                c
                for c in axPlan.collections
                if c.__class__.__name__ == "PolyCollection"
            ]
            drawn = np.concatenate([p.vertices for c in polys for p in c.get_paths()])
            expected = silhouette.quads.reshape(-1, 2)
            self.assertAlmostEqual(
                np.ptp(drawn[:, 0]), 1e3 * np.ptp(expected[:, 0]), places=3
            )
        finally:
            plt.close(fig)

    def test_cases_and_envelope(self):
        results = computeVMTCases(
            [self.f5, self.f5Pushover],
            self.axis,
            shearDir=self.shearDir,
            numStations=20,
        )
        self.assertEqual(list(results), ["grav_000", "pushover_000"])
        single = computeVMTFromF5(
            self.data, self.axis, shearDir=self.shearDir, numStations=20
        )
        np.testing.assert_allclose(results["grav_000"].V, single.V)
        # A -1 g case is exactly -0.4 times the 2.5 g case (linear analysis)
        np.testing.assert_allclose(
            results["pushover_000"].V, -0.4 * single.V, rtol=1e-5, atol=1e-2
        )
        env = computeEnvelope(results)
        np.testing.assert_allclose(env.Vmax, single.V)
        np.testing.assert_allclose(env.Vmin, -0.4 * single.V, rtol=1e-5, atol=1e-2)
        self.assertEqual(list(env.VmaxCase[:-1]), ["grav_000"] * 19)

    def test_cli_multiple_cases(self):
        csv = os.path.join(self.tmp.name, "multi.csv")
        argv = [
            self.f5,
            self.f5Pushover,
            "--axis",
            "4",
            "0",
            "0",
            "4",
            "0",
            "13.8",
            "--shear-dir",
            "0",
            "-1",
            "0",
            "--csv",
            csv,
        ]
        if matplotlib is not None:
            pdf = os.path.join(self.tmp.name, "report.pdf")
            argv += ["--output", pdf]
        self.assertEqual(main(argv), 0)
        for suffix in ("grav_000", "pushover_000", "envelope"):
            self.assertTrue(
                os.path.getsize(os.path.join(self.tmp.name, f"multi_{suffix}.csv")) > 0
            )
        table = np.loadtxt(
            os.path.join(self.tmp.name, "multi_grav_000.csv"),
            delimiter=",",
            comments="#",
            skiprows=8,
        )
        single = computeVMTFromF5(
            self.data, self.axis, shearDir=self.shearDir, numStations=20
        )
        np.testing.assert_allclose(table[:, 1], single.V)
        if matplotlib is not None:
            with open(pdf, "rb") as fp:
                content = fp.read()
            self.assertEqual(len(re.findall(rb"/Type\s*/Page\b", content)), 3)
            png = os.path.join(self.tmp.name, "multi.png")
            self.assertEqual(main(argv[:-2] + ["--output", png]), 0)
            for suffix in ("grav_000", "pushover_000", "envelope"):
                self.assertTrue(
                    os.path.getsize(os.path.join(self.tmp.name, f"multi_{suffix}.png"))
                    > 0
                )

    def test_cli(self):
        csv = os.path.join(self.tmp.name, "vmt.csv")
        argv = [
            self.f5,
            "--axis",
            "4",
            "0",
            "0",
            "4",
            "0",
            "13.8",
            "--shear-dir",
            "0",
            "-1",
            "0",
            "--csv",
            csv,
        ]
        noOutput = argv[:-2]
        if matplotlib is not None:
            png = os.path.join(self.tmp.name, "cli.png")
            argv += ["--output", png]
        self.assertEqual(main(argv), 0)
        table = np.loadtxt(csv, delimiter=",", comments="#", skiprows=8)
        result = computeVMTFromF5(
            self.data, self.axis, shearDir=self.shearDir, numStations=20
        )
        np.testing.assert_allclose(table[:, 1], result.V)
        if matplotlib is not None:
            self.assertTrue(os.path.getsize(png) > 0)
            unitsPng = os.path.join(self.tmp.name, "units.png")
            unitsArgv = noOutput + [
                "--output",
                unitsPng,
                "--force-unit",
                "N",
                "--length-unit",
                "m",
                "--force-scale",
                "1e-3",
                "--length-scale",
                "1000",
            ]
            self.assertEqual(main(unitsArgv), 0)
            self.assertTrue(os.path.getsize(unitsPng) > 0)
        # Scale factors must be positive; argparse rejects them before any work
        for flag in ("--force-scale", "--length-scale"):
            for bad in ("0", "-2", "nan", "abc"):
                with self.assertRaises(SystemExit):
                    main(noOutput + [flag, bad])
        self.assertEqual(
            main(
                [
                    os.path.join(self.tmp.name, "missing.f5"),
                    "--axis",
                    "0",
                    "0",
                    "0",
                    "1",
                    "0",
                    "0",
                ]
            ),
            2,
        )
        junk = os.path.join(self.tmp.name, "junk.f5")
        with open(junk, "wb") as fp:
            fp.write(b"\x07" * 100)
        self.assertEqual(main([junk, "--axis", "0", "0", "0", "1", "0", "0"]), 2)
        self.assertEqual(main([self.f5, "--axis", "0", "0", "0"]), 2)


if __name__ == "__main__":
    unittest.main()
