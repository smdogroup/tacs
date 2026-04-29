from copy import deepcopy
import unittest

import numpy as np
from mpi4py import MPI

from tacs import elements
from tacs import TACS


"""
This is a base class for running MACH StructProblem unit test cases.
This base class will test function evaluations, total
sensitivities, and f5 file writing for the user-specified problems
implemented by the child test case. When the user creates a new test
based on this class only one method, setup_struct_problems, is required
to be defined in the child class. See the virtual method
implementation for this method below for more details.

NOTE: The child class must NOT implement its own setUp method
for the unittest class. This is handled in the base class.
"""


class MACHStructProblemTestCase:
    class MACHStructProblemTest(unittest.TestCase):
        dtype = TACS.dtype

        N_PROCS = 1  # this is how many MPI processes to use for this TestCase.

        FUNC_REFS = {}

        def setUp(self):
            # Default fd/cs step size and tolerances
            # Can be overridden in child class
            if self.dtype == complex:
                self.rtol = 1e-11
                self.atol = 1e-6
                self.dh = 1e-50
            else:
                self.rtol = 1e-2
                self.atol = 1e-4
                self.dh = 1e-5

            # Set the MPI communicator
            if not hasattr(self, "comm"):
                self.comm = MPI.COMM_WORLD

            self.absolute_compare = False

            # Setup struct problems to be tested
            self.struct_probs = self.setup_struct_problems(self.comm)

            # Get the design variable values (assuming all problems share the same arrays)
            # Use the first problem to get the arrays
            first_prob = self.struct_probs[0]
            self.dv0 = first_prob.getOrigDesignVars()
            self.dv_geo = first_prob.getDVGeo()
            self.geo_dvs0 = self.dv_geo.getValues() if self.dv_geo is not None else {}

            # Create temporary vectors for doing fd/cs
            self.dv1 = deepcopy(self.dv0)
            self.dv_pert = {
                dv_name: np.zeros_like(self.dv0[dv_name], dtype=self.dtype)
                for dv_name in self.dv0.keys()
            }

            self.geo_dvs1 = deepcopy(self.geo_dvs0)
            self.geo_dvs_pert = {
                dv_name: np.zeros_like(self.geo_dvs0[dv_name], dtype=self.dtype)
                for dv_name in self.geo_dvs0.keys()
            }

            # Populate fd/cs perturbation vectors based on user-defined method
            self.setup_struct_vecs(self.struct_probs, self.dv_pert, self.geo_dvs_pert)

            # Seed random number generator in tacs for consistent test results
            elements.SeedRandomGenerator(0)

        def setup_struct_problems(self, comm):
            """
            Setup MACH StructProblem objects we will be testing.
            Must be defined in child class that inherits from this class.
            This method should return a list of StructProblem objects to be tested.
            """
            raise NotImplementedError(
                "Child class must implement a 'setup_struct_problems' method"
            )

        def setup_struct_vecs(self, struct_probs, dv_pert_vec, geo_dvs_pert_vec):
            """
            Setup user-defined vectors for analysis and fd/cs sensitivity verification
            """
            # Create temporary dv vec for doing fd/cs
            for dv_name in dv_pert_vec.keys():
                dv_pert_vec[dv_name] = np.ones_like(
                    dv_pert_vec[dv_name], dtype=self.dtype
                )

            # Define perturbation array that moves all geometric design variables
            for dv_name in geo_dvs_pert_vec.keys():
                geo_dvs_pert_vec[dv_name][:] = 1.0

            return

        def test_solve(self):
            """
            Test linear solve and function evaluations
            """
            # solve
            funcs = self.run_solve()

            # Test functions values against historical values
            for prob in self.struct_probs:
                with self.subTest(problem=prob.name):
                    func_list = prob.evalFuncs
                    for func_name in func_list:
                        with self.subTest(function=func_name):
                            func_key = prob.name + "_" + func_name
                            if func_key in self.FUNC_REFS:
                                # Convert to abs value if requested
                                if self.absolute_compare:
                                    funcs[func_key] = abs(funcs[func_key])
                                    self.FUNC_REFS[func_key] = abs(
                                        self.FUNC_REFS[func_key]
                                    )
                                # Test values
                                np.testing.assert_allclose(
                                    funcs[func_key],
                                    self.FUNC_REFS[func_key],
                                    rtol=self.rtol,
                                    atol=self.atol,
                                )

        def test_total_struct_dv_sensitivities(self):
            """
            Test total struct DV sensitivity through adjoint against fd/cs
            """
            # Skip if no struct DVs across all problems
            total_struct_dvs = sum(len(prob.structDVList) for prob in self.struct_probs)
            if total_struct_dvs == 0:
                return

            # Initial solve
            funcs = self.run_solve()

            # Compute sensitivities using adjoint
            func_sens = self.run_sensitivities()

            # Perturb only struct DV keys; leave mass DV keys at original values
            dv1 = deepcopy(self.dv0)
            for prob in self.struct_probs:
                struct_key = prob.varName
                if struct_key in self.dv_pert:
                    dv1[struct_key] = self.perturb_struct_vec(
                        self.dv0[struct_key], self.dv_pert[struct_key]
                    )

            # Run perturbed solution
            funcs_pert = self.run_solve(dv=dv1)

            # Tests cs/fd against sensitivity from adjoint
            for prob in self.struct_probs:
                with self.subTest(problem=prob.name):
                    struct_key = prob.varName
                    func_list = prob.evalFuncs
                    for func_name in func_list:
                        with self.subTest(function=func_name):
                            func_key = prob.name + "_" + func_name
                            if func_key in self.FUNC_REFS:
                                if struct_key not in func_sens[func_key]:
                                    continue
                                # Project exact sens onto struct DV perturbation
                                dfddv_proj = func_sens[func_key][struct_key].dot(
                                    self.dv_pert[struct_key]
                                )
                                # Compute approximate sens
                                fdv_sens_approx = self.compute_fdcs_approx(
                                    funcs_pert[func_key], funcs[func_key]
                                )
                                # Convert to abs value if requested
                                if self.absolute_compare:
                                    dfddv_proj = np.abs(dfddv_proj)
                                    fdv_sens_approx = np.abs(fdv_sens_approx)
                                # Test values
                                np.testing.assert_allclose(
                                    dfddv_proj,
                                    fdv_sens_approx,
                                    rtol=self.rtol,
                                    atol=self.atol,
                                )

        def test_total_promoted_dv_sensitivities(self):
            """
            Test total promoted global DV sensitivity through adjoint against fd/cs.
            Covers all DVs returned by getPromotedDVNames(), which includes mass DVs,
            aux DVs (e.g. load-factor variables), and any other global DVs the user
            chose to promote to problem level.
            """
            # Skip if no promoted DVs across all problems
            total_promoted_dvs = sum(
                len(prob.getPromotedDVNames()) for prob in self.struct_probs
            )
            if total_promoted_dvs == 0:
                return

            # Initial solve
            funcs = self.run_solve()

            # Compute sensitivities using adjoint
            func_sens = self.run_sensitivities()

            # Perturb only promoted DV keys; leave struct DV keys at original values
            dv1 = deepcopy(self.dv0)
            for prob in self.struct_probs:
                for promoted_key in prob.getPromotedDVNames():
                    if promoted_key in self.dv_pert:
                        dv1[promoted_key] = self.perturb_struct_vec(
                            self.dv0[promoted_key], self.dv_pert[promoted_key]
                        )

            # Run perturbed solution
            funcs_pert = self.run_solve(dv=dv1)

            # Tests cs/fd against sensitivity from adjoint
            for prob in self.struct_probs:
                with self.subTest(problem=prob.name):
                    promoted_keys = prob.getPromotedDVNames()
                    func_list = prob.evalFuncs
                    for func_name in func_list:
                        with self.subTest(function=func_name):
                            func_key = prob.name + "_" + func_name
                            if func_key in self.FUNC_REFS:
                                # Sum projected sens over all promoted DVs for this problem
                                dfddv_proj = sum(
                                    func_sens[func_key][promoted_key]
                                    * self.dv_pert[promoted_key]
                                    for promoted_key in promoted_keys
                                    if promoted_key in func_sens[func_key]
                                )
                                # Compute approximate sens
                                fdv_sens_approx = self.compute_fdcs_approx(
                                    funcs_pert[func_key], funcs[func_key]
                                )
                                # Convert to abs value if requested
                                if self.absolute_compare:
                                    dfddv_proj = np.abs(dfddv_proj)
                                    fdv_sens_approx = np.abs(fdv_sens_approx)
                                # Test values
                                np.testing.assert_allclose(
                                    dfddv_proj,
                                    fdv_sens_approx,
                                    rtol=self.rtol,
                                    atol=self.atol,
                                )

        def test_total_dvgeo_sensitivities(self):
            """
            Test total geometric design variable sensitivity through adjoint against fd/cs
            """
            # Skip this check if no dvs were added to model
            total_geo_dvs = self.dv_geo.getNDV()
            if total_geo_dvs == 0:
                return

            # Initial solve
            funcs = self.run_solve()

            # Compute the total derivative w.r.t. geometric design variables using adjoint
            func_sens = self.run_sensitivities()

            # Compute the total derivative w.r.t. geometric design variables using fd/cs
            for dv_name in self.geo_dvs_pert:
                self.geo_dvs1[dv_name] = self.perturb_struct_vec(
                    self.geo_dvs0[dv_name], self.geo_dvs_pert[dv_name]
                )

            # Run perturbed solution
            funcs_pert = self.run_solve(geo_dvs=self.geo_dvs1)

            # Tests cs/fd against sensitivity from adjoint
            for prob in self.struct_probs:
                with self.subTest(problem=prob.name):
                    func_list = prob.evalFuncs
                    for func_name in func_list:
                        with self.subTest(function=func_name):
                            func_key = prob.name + "_" + func_name
                            if func_key in self.FUNC_REFS:
                                dfddvgeo_proj = np.atleast_1d(
                                    np.zeros_like(funcs[func_key])
                                )
                                # Get geometric design variable name for this problem
                                for dv_name in prob.getDVGeo().getVarNames():
                                    if dv_name in func_sens[func_key]:
                                        # project exact sens
                                        dfddvgeo_proj += func_sens[func_key][
                                            dv_name
                                        ].dot(self.geo_dvs_pert[dv_name])
                                # Compute approximate sens
                                fdvgeo_sens_approx = self.compute_fdcs_approx(
                                    funcs_pert[func_key], funcs[func_key]
                                )
                                # Convert to abs value if requested
                                if self.absolute_compare:
                                    dfddvgeo_proj = np.abs(dfddvgeo_proj)
                                    fdvgeo_sens_approx = np.abs(fdvgeo_sens_approx)
                                # Test values
                                np.testing.assert_allclose(
                                    dfddvgeo_proj,
                                    fdvgeo_sens_approx,
                                    rtol=self.rtol,
                                    atol=self.atol,
                                )

        def run_solve(self, dv=None, geo_dvs=None):
            """
            Run a linear solve at specified design point and return functions of interest
            """
            funcs = {}
            for prob in self.struct_probs:
                # Set the design variables if provided
                if dv is not None:
                    prob.setDesignVars(dv)
                else:
                    prob.setDesignVars(self.dv0)

                # Set geometric design variables if provided
                if self.dv_geo is not None:
                    if geo_dvs is not None:
                        self.dv_geo.setDesignVars(geo_dvs)
                    else:
                        self.dv_geo.setDesignVars(self.geo_dvs0)

                # Solve problem
                prob.solve()
                # Evaluate functions
                prob.evalFunctions(funcs)
                # Evaluate constraints
                prob.evalConstraints(funcs)

            return funcs

        def run_sensitivities(self):
            """
            Run a sensitivity solve at specified design point and return sens of functions of interest
            """
            funcs_sens = {}
            funcs = {}
            for prob in self.struct_probs:
                # Set the design variables
                prob.setDesignVars(self.dv0)

                # Solve problem
                prob.solve()
                # Evaluate functions
                prob.evalFunctions(funcs)
                # Evaluate function sensitivities
                prob.evalFunctionsSens(funcs_sens)
                # Evaluate constraint sensitivities
                prob.evalConstraintsSens(funcs_sens)

            return funcs_sens

        def perturb_struct_vec(self, vec_in, vec_pert):
            """
            Perform fd/cs perturbation on struct vector as follows
            vec_out = vec_in + scale * vec_pert

            where:
                scale = dh * 1j, in complex mode
                scale = dh, in real mode
            """
            if self.dtype == complex:
                vec_out = vec_in + self.dh * 1j * vec_pert
            else:
                vec_out = vec_in + self.dh * vec_pert
            return vec_out

        def compute_fdcs_approx(self, vec1, vec0):
            """
            Perform fd/cs calculation to approximate sensitivities

            difference performed as follows:
                sens = imag(vec1) / dh, in complex mode
                sens = (vec1 - vec0) / dh, in real mode
            """
            if self.dtype == complex:
                sens_approx = np.imag(vec1) / self.dh
            else:
                sens_approx = (vec1 - vec0) / self.dh
            return sens_approx
