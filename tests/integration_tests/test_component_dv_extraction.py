"""Tests for pyTACS.getComponentDesignVars and pyTACS.getComponentDesignVarNums."""

import os
import unittest

import numpy as np
from mpi4py import MPI

from tacs import TACS, constitutive, elements, functions, pyTACS
from tacs.utilities import Error

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PARTITIONED_PLATE_BDF = os.path.join(BASE_DIR, "input_files", "partitioned_plate.bdf")
QUAD_TRI_BDF = os.path.join(BASE_DIR, "input_files", "quad_tri_plate.bdf")

# Constructor thickness for each component of partitioned_plate.bdf; the components in
# ACTIVE_COMPONENTS are given active thickness DVs, the others are inactive
PLATE_THICKNESSES = {
    "PLATE.00": 0.010,
    "PLATE.01": 0.012,
    "PLATE.02": 0.014,
    "PLATE.03": 0.016,
}
ACTIVE_COMPONENTS = ["PLATE.00", "PLATE.02"]


def makeShellCon(t, tNum):
    prop = constitutive.MaterialProperties(rho=2780.0, E=73.1e9, nu=0.33, ys=324e6)
    return constitutive.IsoShellConstitutive(prop, t=t, tNum=tNum)


def setupPartitionedPlate(comm, thicknesses):
    """Create an initialized pyTACS instance for partitioned_plate.bdf with a mix of
    active and inactive thickness DVs."""
    fea = pyTACS(PARTITIONED_PLATE_BDF, comm=comm)

    def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
        if compDescript in ACTIVE_COMPONENTS:
            tNum = dvNum
        else:
            tNum = -1
        con = makeShellCon(thicknesses[compDescript], tNum)
        return elements.Quad4Shell(None, con)

    fea.initialize(elemCallBack)
    return fea


class ComponentDVExtractionTest(unittest.TestCase):
    N_PROCS = 2

    def setUp(self):
        self.comm = MPI.COMM_WORLD

    def test_component_dv_nums(self):
        fea = setupPartitionedPlate(self.comm, PLATE_THICKNESSES)
        compDVNums = fea.getComponentDesignVarNums()
        expected = {
            "PLATE.00": {"t": 0},
            "PLATE.01": {"t": -1},
            "PLATE.02": {"t": 1},
            "PLATE.03": {"t": -1},
        }
        self.assertEqual(compDVNums, expected)

    def test_initial_values(self):
        fea = setupPartitionedPlate(self.comm, PLATE_THICKNESSES)
        compDVs = fea.getComponentDesignVars()
        self.assertEqual(sorted(compDVs.keys()), sorted(PLATE_THICKNESSES.keys()))
        for descript, t in PLATE_THICKNESSES.items():
            self.assertEqual(list(compDVs[descript].keys()), ["t"])
            self.assertNotIsInstance(compDVs[descript]["t"], np.ndarray)
            np.testing.assert_allclose(compDVs[descript]["t"], t, rtol=1e-12)

    def test_duplicate_descripts_error(self):
        fea = setupPartitionedPlate(self.comm, PLATE_THICKNESSES)
        # Simulate a model whose BDF labels two property groups identically
        fea.compDescripts = list(fea.compDescripts)
        fea.compDescripts[1] = fea.compDescripts[0]
        with self.assertRaises(Error):
            fea.getComponentDesignVars()
        with self.assertRaises(Error):
            fea.getComponentDesignVarNums()

    def test_multiple_constitutive_error(self):
        fea = pyTACS(QUAD_TRI_BDF, comm=self.comm)

        def elemCallBack(
            dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs
        ):
            elems = []
            for descript in elemDescripts:
                # Deliberately create a NEW constitutive object per element type
                con = makeShellCon(0.01, -1)
                if descript == "CQUAD4":
                    elems.append(elements.Quad4Shell(None, con))
                elif descript == "CTRIA3":
                    elems.append(elements.Tri3Shell(None, con))
            return elems

        fea.initialize(elemCallBack)
        with self.assertRaises(Error):
            fea.getComponentDesignVars()

    def test_shared_constitutive_ok(self):
        fea = pyTACS(QUAD_TRI_BDF, comm=self.comm)

        def elemCallBack(
            dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs
        ):
            # ONE constitutive object shared by both element types is the intended usage
            con = makeShellCon(0.01, dvNum)
            elems = []
            for descript in elemDescripts:
                if descript == "CQUAD4":
                    elems.append(elements.Quad4Shell(None, con))
                elif descript == "CTRIA3":
                    elems.append(elements.Tri3Shell(None, con))
            return elems

        fea.initialize(elemCallBack)
        compDVs = fea.getComponentDesignVars()
        self.assertEqual(len(compDVs), 1)
        (groupValues,) = compDVs.values()
        np.testing.assert_allclose(groupValues["t"], 0.01, rtol=1e-12)

    def test_unimplemented_class_gives_empty_entry(self):
        class NoGroupShellConstitutive(constitutive.IsoShellConstitutive):
            """Simulates a constitutive class that predates the DV group API."""

            def getNumDesignVarGroups(self):
                return 0

            def getDesignVarGroups(self):
                return {}

            def getDesignVarGroupDVNums(self):
                return {}

        fea = pyTACS(PARTITIONED_PLATE_BDF, comm=self.comm)

        def elemCallBack(
            dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs
        ):
            prop = constitutive.MaterialProperties(
                rho=2780.0, E=73.1e9, nu=0.33, ys=324e6
            )
            con = NoGroupShellConstitutive(prop, t=0.01, tNum=dvNum)
            return elements.Quad4Shell(None, con)

        fea.initialize(elemCallBack)
        # The components have DVs but report no groups, so each entry is empty (and a
        # warning naming the class is printed on the root proc)
        compDVs = fea.getComponentDesignVars()
        for groupValues in compDVs.values():
            self.assertEqual(groupValues, {})

    def test_partially_implemented_class_does_not_raise(self):
        class PartialGroupShellConstitutive(constitutive.IsoShellConstitutive):
            """Simulates a class that reports a group but whose group DV nums don't cover
            all of its active design variables."""

            def getDesignVarGroups(self):
                return {"t": 0.01}

            def getDesignVarGroupDVNums(self):
                return {"t": -1}

        fea = pyTACS(PARTITIONED_PLATE_BDF, comm=self.comm)

        def elemCallBack(
            dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs
        ):
            prop = constitutive.MaterialProperties(
                rho=2780.0, E=73.1e9, nu=0.33, ys=324e6
            )
            con = PartialGroupShellConstitutive(prop, t=0.01, tNum=dvNum)
            return elements.Quad4Shell(None, con)

        fea.initialize(elemCallBack)
        # The class claims one group but its group DV nums cover none of its active DVs, so
        # a warning naming the class is printed on the root proc; the call must still succeed
        # and return the stale (constructor) group values.
        compDVs = fea.getComponentDesignVars()
        for groupValues in compDVs.values():
            np.testing.assert_allclose(groupValues["t"], 0.01, rtol=1e-12)


class ComponentDVRoundTripTest(unittest.TestCase):
    N_PROCS = 2

    def setUp(self):
        self.comm = MPI.COMM_WORLD
        # Function-value comparison tolerance between the original and rebuilt models
        self.rtol = 1e-11 if TACS.dtype is complex else 1e-9

    @staticmethod
    def setupProblem(fea):
        problem = fea.createStaticProblem("gravity")
        problem.addInertialLoad(np.array([0.0, 0.0, -9.81]))
        problem.addFunction("mass", functions.StructuralMass)
        problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=100.0)
        return problem

    def test_round_trip(self):
        fea1 = setupPartitionedPlate(self.comm, PLATE_THICKNESSES)
        problem1 = self.setupProblem(fea1)

        # Perturb the active DVs so the assembler state differs from the constructor values
        xLocal = problem1.getDesignVars()
        problem1.setDesignVars(1.5 * xLocal)
        problem1.solve()
        funcs1 = {}
        problem1.evalFunctions(funcs1)

        compDVs = fea1.getComponentDesignVars()

        # The dictionary must be identical on every rank
        allCompDVs = self.comm.allgather(compDVs)
        for otherCompDVs in allCompDVs:
            self.assertEqual(sorted(otherCompDVs.keys()), sorted(compDVs.keys()))
            for descript in compDVs:
                self.assertEqual(
                    sorted(otherCompDVs[descript].keys()),
                    sorted(compDVs[descript].keys()),
                )
                for name in compDVs[descript]:
                    np.testing.assert_array_equal(
                        otherCompDVs[descript][name], compDVs[descript][name]
                    )

        # Active components must show the perturbed values, inactive components the
        # constructor values
        for descript, t0 in PLATE_THICKNESSES.items():
            if descript in ACTIVE_COMPONENTS:
                expected = 1.5 * t0
            else:
                expected = t0
            np.testing.assert_allclose(compDVs[descript]["t"], expected, rtol=1e-12)

        # Rebuild a second model from the extracted dictionary, exactly as a user would in
        # a new TACS execution, and check it produces the same function values
        fea2 = pyTACS(PARTITIONED_PLATE_BDF, comm=self.comm)

        def elemCallBack(
            dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs
        ):
            con = makeShellCon(compDVs[compDescript]["t"], dvNum)
            return elements.Quad4Shell(None, con)

        fea2.initialize(elemCallBack)
        problem2 = self.setupProblem(fea2)
        problem2.solve()
        funcs2 = {}
        problem2.evalFunctions(funcs2)

        self.assertEqual(sorted(funcs1.keys()), sorted(funcs2.keys()))
        for key in funcs1:
            np.testing.assert_allclose(funcs1[key], funcs2[key], rtol=self.rtol)


class ComponentDVArrayGroupOverlayTest(unittest.TestCase):
    """Exercises the array-valued-group overlay branch of getComponentDesignVars.

    All other tests in this module use IsoShellConstitutive, whose only DV group is the
    scalar "t", so the isinstance(nums, np.ndarray) branch that overlays array groups
    (e.g. ply_fractions) is otherwise untested end-to-end under MPI.
    """

    N_PROCS = 2

    def setUp(self):
        self.comm = MPI.COMM_WORLD

    @staticmethod
    def makeSmearedCompositeCon(dvNum):
        prop = constitutive.MaterialProperties(
            rho=1550.0,
            specific_heat=921.096,
            E1=54e9,
            E2=18e9,
            nu12=0.25,
            G12=9e9,
            G13=9e9,
            G23=9e9,
            Xt=2410.0e6,
            Xc=1040.0e6,
            Yt=73.0e6,
            Yc=173.0e6,
            S12=71.0e6,
            alpha=24.0e-6,
            kappa=230.0,
        )
        ply = constitutive.OrthotropicPly(0.1, prop)
        plyAngles = np.deg2rad([0.0, 45.0, 90.0]).astype(TACS.dtype)
        plyFractions = np.array([0.5, 0.3, 0.2], dtype=TACS.dtype)
        thicknessDVNum = dvNum
        plyFractionDVNums = np.array([dvNum + 1, dvNum + 2, dvNum + 3], dtype=np.intc)
        return constitutive.SmearedCompositeShellConstitutive(
            [ply] * 3,
            0.01,
            plyAngles,
            plyFractions,
            thickness_dv_num=thicknessDVNum,
            ply_fraction_dv_nums=plyFractionDVNums,
        )

    def test_array_group_overlay(self):
        fea = pyTACS(PARTITIONED_PLATE_BDF, comm=self.comm)

        def elemCallBack(
            dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs
        ):
            con = self.makeSmearedCompositeCon(dvNum)
            return elements.Quad4Shell(None, con)

        fea.initialize(elemCallBack)
        problem = fea.createStaticProblem("perturbed")

        # Perturb every active DV (thickness and ply fractions alike) away from the
        # constructor values
        xLocal = problem.getDesignVars()
        perturbed = 1.2 * xLocal
        problem.setDesignVars(perturbed)

        compDVs = fea.getComponentDesignVars()

        # The dictionary must be identical on every rank
        allCompDVs = self.comm.allgather(compDVs)
        for otherCompDVs in allCompDVs:
            self.assertEqual(sorted(otherCompDVs.keys()), sorted(compDVs.keys()))
            for descript in compDVs:
                self.assertEqual(
                    sorted(otherCompDVs[descript].keys()),
                    sorted(compDVs[descript].keys()),
                )
                for name in compDVs[descript]:
                    np.testing.assert_array_equal(
                        otherCompDVs[descript][name], compDVs[descript][name]
                    )

        expectedThickness = 1.2 * 0.01
        expectedPlyFractions = 1.2 * np.array([0.5, 0.3, 0.2])
        for descript in PLATE_THICKNESSES:
            groupValues = compDVs[descript]
            self.assertNotIsInstance(groupValues["thickness"], np.ndarray)
            np.testing.assert_allclose(
                groupValues["thickness"], expectedThickness, rtol=1e-12
            )
            self.assertIsInstance(groupValues["ply_fractions"], np.ndarray)
            np.testing.assert_allclose(
                groupValues["ply_fractions"], expectedPlyFractions, rtol=1e-12
            )


if __name__ == "__main__":
    unittest.main()
