"""
Unit tests for how the MACH StructProblem exposes promoted global design
variables to pyoptsparse: the design variable names, the data types of the
design variable dictionaries, and the ``wrt`` lists of the constraints.

The model is the four component partitioned plate, with a gravity load and a
pressure load controlled by global design variables. No problems are solved.
"""

import os
import unittest

import numpy as np

from tacs import pyTACS, constitutive, elements, functions

try:
    from tacs.mach import StructProblem
except ImportError:
    StructProblem = None

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/partitioned_plate.bdf")

# Gravity vector and pressure magnitude controlled by global DVs
g = np.array([0.0, 0.0, -9.81])
P = 1e5


@unittest.skipIf(StructProblem is None, "MACH dependencies are not installed")
class StructProblemPyOptTest(unittest.TestCase):
    # Every StructProblem getter broadcasts its result, so the same assertions
    # run on every rank
    N_PROCS = 2

    def setUp(self):
        self.structProb = self.setup_struct_problem()

    def setup_struct_problem(self, name="load", **kwargs):
        """
        Create a StructProblem promoting every global DV. Any keyword arguments
        are passed on to the StructProblem.
        """
        FEAAssembler = pyTACS(bdf_file)

        # Global DVs: an array, an array of length 1 and a scalar
        gDVNums = FEAAssembler.addGlobalDV("g", g, lower=-20.0, upper=20.0)
        pDVNums = FEAAssembler.addGlobalDV("pressure", [P], lower=0.0, upper=2 * P)
        FEAAssembler.addGlobalDV("unused", 1.0, lower=0.0, upper=2.0)

        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
        ):
            prop = constitutive.MaterialProperties(
                rho=2500.0, E=70e9, nu=0.3, ys=464.0e6
            )
            con = constitutive.IsoShellConstitutive(
                prop, t=0.01, tNum=dv_num, tlb=0.001, tub=0.1
            )
            return elements.Quad4Shell(None, con)

        FEAAssembler.initialize(elem_call_back)

        staticProb = FEAAssembler.createStaticProblem(name)
        staticProb.addFunction("mass", functions.StructuralMass)
        staticProb.addFunction("compliance", functions.Compliance)
        staticProb.addInertialLoad(g, inertiaVecDVNums=gDVNums)
        compIDs = FEAAssembler.selectCompIDs()
        staticProb.addPressureToComponents(compIDs, P, pressureDVNums=pDVNums[0])

        structProb = StructProblem(
            staticProb,
            FEAAssembler,
            promotedGlobalDVs=["g", "pressure", "unused"],
            **kwargs,
        )

        # A linear constraint on the thicknesses, which depends on none of the
        # promoted DVs
        constr = FEAAssembler.createDVConstraint("thickness")
        constr.addConstraint("sum", lower=0.0, upper=1.0)
        structProb.addConstraint(constr)

        self.numComps = FEAAssembler.getNumComponents()
        return structProb

    def test_var_names(self):
        self.assertEqual(self.structProb.getVarName(), "struct")
        self.assertEqual(
            self.structProb.getPromotedDVNames(),
            ["g_load", "pressure_load", "unused_load"],
        )
        self.assertEqual(
            self.structProb.getVarNames(),
            ["g_load", "pressure_load", "unused_load", "struct"],
        )

    def test_design_var_dicts_are_float_arrays(self):
        """
        Every entry has to be a float64 array with the shape of the DV, whatever
        the numpy version. In particular an array DV of length 1 must stay an
        array, and not become a scalar.
        """
        shapes = {
            "g_load": (3,),
            "pressure_load": (1,),
            "unused_load": (),
            "struct": (self.numComps,),
        }
        lb, ub = self.structProb.getDesignVarRange()
        dvDicts = {
            "value": self.structProb.getOrigDesignVars(),
            "lower": lb,
            "upper": ub,
            "scale": self.structProb.getDesignVarScales(),
        }
        for dictName, dvDict in dvDicts.items():
            with self.subTest(dict=dictName):
                self.assertEqual(set(dvDict.keys()), set(shapes.keys()))
                for dvName, shape in shapes.items():
                    with self.subTest(dv=dvName):
                        self.assertIsInstance(dvDict[dvName], np.ndarray)
                        self.assertEqual(dvDict[dvName].dtype, np.float64)
                        self.assertEqual(dvDict[dvName].shape, shape)

        np.testing.assert_array_equal(dvDicts["value"]["g_load"], g)
        np.testing.assert_array_equal(dvDicts["value"]["pressure_load"], [P])
        np.testing.assert_array_equal(dvDicts["lower"]["g_load"], [-20.0] * 3)
        np.testing.assert_array_equal(dvDicts["upper"]["pressure_load"], [2 * P])


if __name__ == "__main__":
    unittest.main()
