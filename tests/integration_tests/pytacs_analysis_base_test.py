import numpy as np
from tacs import TACS
import unittest
from mpi4py import MPI

'''
This is a base class for running pytacs unit test cases.
This base class will test function evaluations and total 
sensitivities for the user-specified problems implimented by
the child test case. When the user creates a new test based 
on this class four methods are required to be defined in the child class. 

    1. setup_pytacs
    2. setup_tacs_problems
    3. setup_tacs_vecs
    4. setup_funcs
    
See the virtual method implementations for each method 
below for more details.

NOTE: The child class must NOT implement its own setUp method 
for the unittest class. This is handled in the base class.
'''

class PyTACSTestCase:

    class PyTACSTest(unittest.TestCase):
        def setUp(self):
            self.dtype = TACS.dtype

            # Default fd/cs step size and tolerances
            # Can be overridden in child class
            if self.dtype == complex:
                self.rtol = 1e-11
                self.atol = 1e-8
                self.dh = 1e-50
            else:
                self.rtol = 1e-2
                self.atol = 1e-4
                self.dh = 1e-5

            # Set the MPI communicator
            if not hasattr(self, 'comm'):
                self.comm = MPI.COMM_WORLD

            # Setup user-specified assembler for this test
            self.fea_assembler = self.setup_pytacs(self.comm, self.dtype)

            # Get the design variable values
            self.dv0 = self.fea_assembler.getOrigDesignVars()
            # Initial nodal location vector
            self.xpts0 = self.fea_assembler.getOrigCoordinates()
            # Create temporary dv vec for doing fd/cs
            self.dv1 = np.zeros_like(self.dv0, dtype=self.dtype)
            self.dv_pert = np.zeros_like(self.dv0, dtype=self.dtype)
            self.xpts1 = np.zeros_like(self.xpts0, dtype=self.dtype)
            self.xpts_pert = np.zeros_like(self.xpts0, dtype=self.dtype)

            # Setup tacs problems to be tested
            self.tacs_probs = self.setup_tacs_problems(self.fea_assembler)
            # Populate fd/cs perturbation vectors based on user-defined method
            self.setup_tacs_vecs(self.fea_assembler, self.dv_pert, self.xpts_pert)
            # Create the function list
            self.func_list, self.func_ref = self.setup_funcs(self.fea_assembler, self.tacs_probs)

        def setup_pytacs(self, comm, dtype):
            """
            Setup pytacs object for problems we will be testing.
            Must be defined in child class that inherits from this class.
            """
            raise NotImplementedError("Child class must implement a 'setup_pytacs' method")
            return

        def setup_tacs_problems(self, fea_assembler):
            """
            Setup tacs problems objects that describe different problem types we will be testing.
            Must be defined in child class that inherits from this class.
            """
            raise NotImplementedError("Child class must implement a 'setup_tacs_problems' method")
            return

        def setup_tacs_vecs(self, fea_assembler, dv_pert_vec, xpts_pert_vec):
            """
            Setup user-defined vectors for analysis and fd/cs sensitivity verification.
            Must be defined in child class that inherits from this class.
            """
            raise NotImplementedError("Child class must implement a 'setup_tacs_vecs' method")
            return

        def setup_funcs(self, fea_assembler, problems):
            """
            Create a list of functions to be tested and their reference values for the problem.
            Must be defined in child class that inherits from this class.
            """
            raise NotImplementedError("Child class must implement a 'setup_funcs' method")
            return

        def test_solve(self):
            """
            Test linear solve and function evaluations
            """
            # solve
            funcs = self.run_solve()

            # Test functions values against historical values
            for prob in self.tacs_probs:
                with self.subTest(problem=prob.name):
                    for func_name in self.func_list:
                        with self.subTest(function=func_name):
                            func_key = prob.name + '_' + func_name
                            np.testing.assert_allclose(funcs[func_key], self.func_ref[func_key],
                                                       rtol=self.rtol, atol=self.atol)

        def test_total_dv_sensitivities(self):
            """
            Test total dv sensitivity through adjoint against fd/cs
            """
            # Initial solve
            funcs = self.run_solve()

            # Compute the total derivative w.r.t. material design variables using adjoint
            func_sens = self.run_sensitivities()

            # Compute the total derivative w.r.t. material design variables using fd/cs
            self.dv1 = self.perturb_tacs_vec(self.dv0, self.dv_pert)
            # Run perturbed solution
            funcs_pert = self.run_solve(dv=self.dv1)

            # Tests cs/fd against sensitivity from adjoint
            for prob in self.tacs_probs:
                with self.subTest(problem=prob.name):
                    for func_name in self.func_list:
                        with self.subTest(function=func_name):
                            func_key = prob.name + '_' + func_name
                            # project exact sens
                            dfddv_proj = func_sens[func_key]['struct'].dot(self.dv_pert)
                            # Get contribution across all procs
                            dfddv_proj = self.comm.allreduce(dfddv_proj, op=MPI.SUM)
                            # Compute approximate sens
                            fdv_sens_approx = self.compute_fdcs_approx(funcs_pert[func_key], funcs[func_key])
                            np.testing.assert_allclose(dfddv_proj, fdv_sens_approx,
                                                       rtol=self.rtol, atol=self.atol)

        def test_total_xpt_sensitivities(self):
            """
            Test total xpt sensitivity through adjoint against fd/cs
            """
            # Initial solve
            funcs = self.run_solve()

            # Compute the total derivative w.r.t. material design variables using adjoint
            func_sens = self.run_sensitivities()

            # Compute the total derivative w.r.t. nodal xpt locations using fd/cs
            self.xpts1 = self.perturb_tacs_vec(self.xpts0, self.xpts_pert)
            # Run perturbed solution
            funcs_pert = self.run_solve(xpts=self.xpts1)

            # Tests cs/fd against sensitivity from adjoint
            for prob in self.tacs_probs:
                with self.subTest(problem=prob.name):
                    for func_name in self.func_list:
                        with self.subTest(function=func_name):
                            func_key = prob.name + '_' + func_name
                            # project exact sens
                            dfdx_proj = func_sens[func_key]['Xpts'].dot(self.xpts_pert)
                            # Get contribution across all procs
                            dfdx_proj = self.comm.allreduce(dfdx_proj, op=MPI.SUM)
                            # Compute approximate sens
                            f_xpt_sens_approx = self.compute_fdcs_approx(funcs_pert[func_key], funcs[func_key])
                            np.testing.assert_allclose(dfdx_proj, f_xpt_sens_approx,
                                                       rtol=self.rtol, atol=self.atol)

        def run_solve(self, dv=None, xpts=None):
            """
            Run a linear solve at specified design point and return functions of interest
            """
            if dv is None:
                dv = self.dv0

            if xpts is None:
                xpts = self.xpts0

            # Solve each problem and evaluate functions
            funcs = {}
            for prob in self.tacs_probs:
                # Set the design variables
                prob.setDesignVars(dv)
                # Set node locations
                prob.setCoordinates(xpts)
                # Solve problem
                prob.solve()
                # Evaluate functions
                prob.evalFunctions(funcs)

            return funcs

        def run_sensitivities(self, dv=None, xpts=None):
            """
            Run a sensitivity solve at specified design point and return sens of  functions of interest
            """
            if dv is None:
                dv = self.dv0

            if xpts is None:
                xpts = self.xpts0

            # Solve each problem and evaluate function sens
            funcs_sens = {}
            funcs = {}
            for prob in self.tacs_probs:
                # Set the design variables
                prob.setDesignVars(dv)
                # Set node locations
                prob.setCoordinates(xpts)
                # Solve problem
                prob.solve()
                # Evaluate functions
                prob.evalFunctions(funcs)
                prob.evalFunctionsSens(funcs_sens)

            return funcs_sens

        def perturb_tacs_vec(self, vec_in, vec_pert):
            """
            Perform fd/cs perturbation on tacs vector as follows
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
