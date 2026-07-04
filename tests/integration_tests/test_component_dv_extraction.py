"""Tests for pyTACS.getComponentDesignVars and pyTACS.getComponentDesignVarNums."""

import os
import unittest

import numpy as np
from mpi4py import MPI

from tacs import constitutive, elements, pyTACS
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


if __name__ == "__main__":
    unittest.main()
