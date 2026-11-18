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
from tacs._testing import fails_at_version

try:
    from tacs.mach import StructProblem
except ImportError:
    StructProblem = None

try:
    from pyoptsparse import Optimization
except ImportError:
    Optimization = None

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
        """
        Given a StructProblem that promotes every global DV,
        when getStructDVName, getPromotedDVNames and getVarNames are called,
        then they return the struct DV name, the promoted DV names suffixed with the problem name,
        and both together with the struct DV last.
        """
        self.assertEqual(self.structProb.getStructDVName(), "struct")
        self.assertEqual(
            self.structProb.getPromotedDVNames(),
            ["g_load", "pressure_load", "unused_load"],
        )
        self.assertEqual(
            self.structProb.getVarNames(),
            ["g_load", "pressure_load", "unused_load", "struct"],
        )

    @fails_at_version("3.15")
    def test_get_var_name_is_deprecated(self):
        """Assert that the getVarName alias raises a DeprecationWarning and returns the same name as getStructDVName."""
        with self.assertWarns(DeprecationWarning):
            varName = self.structProb.getVarName()
        self.assertEqual(varName, self.structProb.getStructDVName())

    def test_design_var_dicts_are_float_arrays(self):
        """
        Given a StructProblem with an array DV, an array DV of length 1 and a scalar DV promoted,
        when its value, lower bound, upper bound and scale dictionaries are retrieved,
        then every entry is a float64 array with the shape of its DV, whatever the numpy version,
        so the length-1 array DV stays an array, and the entries hold the values given to addGlobalDV.
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

    def test_excluded_dvs_are_not_in_var_names(self):
        """
        Given a StructProblem with a promoted DV and the struct DVs excluded,
        when getVarNames and getPromotedDVNames are called,
        then the excluded DVs are missing from the variable names but still listed as promoted DVs.
        """
        structProb = self.setup_struct_problem(excludedDVs=["g", "struct"])
        self.assertEqual(structProb.getVarNames(), ["pressure_load", "unused_load"])
        # The promoted DVs themselves are unchanged
        self.assertEqual(
            structProb.getPromotedDVNames(),
            ["g_load", "pressure_load", "unused_load"],
        )

    def test_unknown_excluded_dv_raises(self):
        """Assert that creating a StructProblem raises a ValueError when excludedDVs names a DV by its pyOpt key or names a DV that does not exist."""
        # Excluded DVs are named like the promoted DVs, not by their pyOpt key
        for name in ["g_load", "notADV"]:
            with self.subTest(name=name):
                with self.assertRaises(ValueError):
                    self.setup_struct_problem(excludedDVs=[name])

    @unittest.skipIf(Optimization is None, "pyoptsparse is not installed")
    def test_empty_promoted_blocks_are_not_in_constraint_wrt(self):
        """
        Given a StructProblem with a linear constraint on the thicknesses only,
        when its variables and constraints are added to a pyoptsparse Optimization,
        then every promoted DV is added as a variable but the constraint's wrt list holds only the struct DVs.
        """
        optProb = Optimization("opt", lambda x: ({}, False), comm=self.structProb.comm)
        self.structProb.addVariablesPyOpt(optProb)
        self.assertEqual(
            list(optProb.variables.keys()),
            ["g_load", "pressure_load", "unused_load", "struct"],
        )

        self.structProb.addConstraintsPyOpt(optProb)
        # The constraint only depends on the thicknesses, the promoted DVs
        # would only add empty Jacobian blocks
        self.assertEqual(optProb.constraints["thickness_sum"].wrt, ["struct"])

    @unittest.skipIf(Optimization is None, "pyoptsparse is not installed")
    def test_excluded_dvs_are_not_added_to_pyopt(self):
        """
        Given a StructProblem with a promoted DV and the struct DVs excluded,
        when its variables and constraints are added to each of two pyoptsparse Optimizations,
        then neither optimization has the excluded DVs as variables or in the constraint's wrt list.
        """
        structProb = self.setup_struct_problem(excludedDVs=["g", "struct"])
        # The exclusions belong to the problem, so they apply to every
        # optimization the problem is added to
        for _ in range(2):
            optProb = Optimization("opt", lambda x: ({}, False), comm=structProb.comm)
            structProb.addVariablesPyOpt(optProb)
            self.assertEqual(
                list(optProb.variables.keys()), ["pressure_load", "unused_load"]
            )

            structProb.addConstraintsPyOpt(optProb)
            self.assertEqual(optProb.constraints["thickness_sum"].wrt, [])

    def test_excluded_dvs_are_not_in_sensitivities(self):
        """
        Given StructProblems with no DVs excluded and with a promoted DV and the struct DVs excluded,
        when the function and linear constraint sensitivities are evaluated,
        then the sensitivity dictionaries hold exactly the DVs that are not excluded.
        """
        allKeys = {"g_load", "pressure_load", "unused_load", "struct"}
        for excluded in ([], ["g", "struct"]):
            with self.subTest(excludedDVs=excluded):
                structProb = self.setup_struct_problem(excludedDVs=excluded)
                expected = allKeys - {
                    "struct" if n == "struct" else f"{n}_load" for n in excluded
                }
                structProb.solve()
                funcsSens = {}
                structProb.evalFunctionsSens(funcsSens)
                for funcKey in ["load_mass", "load_compliance"]:
                    self.assertEqual(set(funcsSens[funcKey].keys()), expected)
                conSens = {}
                structProb.evalConstraintsSens(conSens, linear=True)
                self.assertEqual(set(conSens["thickness_sum"].keys()), expected)

    @unittest.skipIf(Optimization is None, "pyoptsparse is not installed")
    def test_pyopt_accepts_objective_gradient_with_excluded_dvs(self):
        """
        Given a StructProblem with a promoted DV excluded and its variables added to a finalized pyoptsparse Optimization,
        when its compliance sensitivities are passed to processObjectiveGradient,
        then pyoptsparse accepts them and returns a nonzero gradient of the optimization's size.

        pyoptsparse raises an error for a gradient entry of a variable group that was never added to the optimization.
        """
        structProb = self.setup_struct_problem(excludedDVs=["g"])
        optProb = Optimization("opt", lambda x: ({}, False), comm=structProb.comm)
        structProb.addVariablesPyOpt(optProb)
        optProb.addObj("load_compliance")
        optProb.finalize()

        structProb.solve()
        funcsSens = {}
        structProb.evalFunctionsSens(funcsSens)
        gobj = optProb.processObjectiveGradient(funcsSens)
        self.assertEqual(np.size(gobj), optProb.ndvs)
        self.assertTrue(np.any(gobj != 0.0))

    def test_excluded_dvs_can_still_be_set(self):
        """Assert that StructProblem.setDesignVars still updates a promoted DV that is excluded from the optimization."""
        structProb = self.setup_struct_problem(excludedDVs=["g"])
        structProb.setDesignVars({"g_load": np.array([0.0, 0.0, -5.0])})
        x = structProb.convertDesignVecToDict(structProb.staticProblem.getDesignVars())
        np.testing.assert_array_equal(np.real(x["g_load"]), [0.0, 0.0, -5.0])

    def test_is_zero_jac_block(self):
        """
        Assert that _isZeroJacBlock treats a sparse block as zero when it stores no nonzeros, whether or not the
        constraint is linear, and a dense block as zero only when it is all zeros and the constraint is linear.
        """
        import scipy.sparse as sparse

        isZero = self.structProb._isZeroJacBlock
        # Sparse blocks fix the sparsity pattern, so an empty one stays empty
        self.assertTrue(isZero(sparse.csr_matrix((2, 3)), False))
        self.assertFalse(isZero(sparse.csr_matrix(np.eye(3)), False))
        # A stored zero is not a nonzero
        self.assertTrue(isZero(sparse.csr_matrix(([0.0], ([0], [0])), (2, 3)), False))
        # A dense zero block of a nonlinear constraint may just be zero here
        self.assertFalse(isZero(np.zeros((2, 3)), False))
        self.assertTrue(isZero(np.zeros((2, 3)), True))
        self.assertFalse(isZero(np.ones((2, 3)), True))


if __name__ == "__main__":
    unittest.main()
