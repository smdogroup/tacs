import os
import tempfile
import unittest

import numpy as np
from mpi4py import MPI

from tacs import TACS
from tacs import problems

"""
This is a base class for running pytacs unit test cases.
This base class will test function evaluations, total
sensitivities, and f5 file writing for the user-specified problems
implemented by the child test case. When the user creates a new test
based on this class only one method, setup_tacs_problems, is required
to be defined in the child class. See the virtual method
implementation for this method below for more details.

NOTE: The child class must NOT implement its own setUp method
for the unittest class. This is handled in the base class.
"""


class PyTACSTestCase:
    class PyTACSTest(unittest.TestCase):
        dtype = TACS.dtype

        N_PROCS = 1  # this is how many MPI processes to use for this TestCase.

        FUNC_REFS = {}

        def setUp(self):
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
            if not hasattr(self, "comm"):
                self.comm = MPI.COMM_WORLD

            # Setup tacs problems to be tested
            self.tacs_probs, self.fea_assembler = self.setup_tacs_problems(self.comm)

            # Get the design variable values
            self.dv0 = self.fea_assembler.getOrigDesignVars()
            # Initial nodal location vector
            self.xpts0 = self.fea_assembler.getOrigNodes()
            # Create temporary dv vec for doing fd/cs
            self.dv1 = np.zeros_like(self.dv0, dtype=self.dtype)
            self.dv_pert = np.zeros_like(self.dv0, dtype=self.dtype)
            self.xpts1 = np.zeros_like(self.xpts0, dtype=self.dtype)
            self.xpts_pert = np.zeros_like(self.xpts0, dtype=self.dtype)

            # Populate fd/cs perturbation vectors based on user-defined method
            self.setup_tacs_vecs(self.fea_assembler, self.dv_pert, self.xpts_pert)

        def setup_tacs_problems(self, comm):
            """
            Setup pytacs object and problems we will be testing.
            Must be defined in child class that inherits from this class.
            This method should return a list of problems to be tested,
            and the corresponding pytacs assembler used to create these problems.
            """
            raise NotImplementedError(
                "Child class must implement a 'setup_tacs_problems' method"
            )
            return

        def setup_tacs_vecs(self, fea_assembler, dv_pert_vec, xpts_pert_vec):
            """
            Setup user-defined vectors for analysis and fd/cs sensitivity verification
            """
            # Create temporary dv vec for doing fd/cs
            dv_pert_vec[:] = 1.0

            # Define perturbation array that moves all nodes on shell
            xpts = fea_assembler.getOrigNodes()
            xpts_pert_vec[:] = xpts

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
                    func_list = prob.getFunctionKeys()
                    for func_name in func_list:
                        with self.subTest(function=func_name):
                            func_key = prob.name + "_" + func_name
                            if func_key in self.FUNC_REFS:
                                np.testing.assert_allclose(
                                    funcs[func_key],
                                    self.FUNC_REFS[func_key],
                                    rtol=self.rtol,
                                    atol=self.atol,
                                )

        def test_total_dv_sensitivities(self):
            """
            Test total dv sensitivity through adjoint against fd/cs
            """
            # Skip this check if no dvs were added to model
            num_dvs = self.fea_assembler.getTotalNumDesignVars()
            if num_dvs == 0:
                return

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
                    func_list = prob.getFunctionKeys()
                    for func_name in func_list:
                        with self.subTest(function=func_name):
                            func_key = prob.name + "_" + func_name
                            if func_key in self.FUNC_REFS:
                                # project exact sens
                                dfddv_proj = func_sens[func_key]["struct"].dot(
                                    self.dv_pert
                                )
                                # Get contribution across all procs
                                dfddv_proj = self.comm.allreduce(dfddv_proj, op=MPI.SUM)
                                # Compute approximate sens
                                fdv_sens_approx = self.compute_fdcs_approx(
                                    funcs_pert[func_key], funcs[func_key]
                                )
                                np.testing.assert_allclose(
                                    dfddv_proj,
                                    fdv_sens_approx,
                                    rtol=self.rtol,
                                    atol=self.atol,
                                )

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
                    func_list = prob.getFunctionKeys()
                    for func_name in func_list:
                        with self.subTest(function=func_name):
                            func_key = prob.name + "_" + func_name
                            if func_key in self.FUNC_REFS:
                                # project exact sens
                                dfdx_proj = func_sens[func_key]["Xpts"].dot(
                                    self.xpts_pert
                                )
                                # Get contribution across all procs
                                dfdx_proj = self.comm.allreduce(dfdx_proj, op=MPI.SUM)
                                # Compute approximate sens
                                f_xpt_sens_approx = self.compute_fdcs_approx(
                                    funcs_pert[func_key], funcs[func_key]
                                )
                                np.testing.assert_allclose(
                                    dfdx_proj,
                                    f_xpt_sens_approx,
                                    rtol=self.rtol,
                                    atol=self.atol,
                                )

        def test_write_solution(self):
            """
            Test f5 solution writing procedure
            """
            # Create temporary directory to write f5 file to (only on root)
            tmp_dir = None
            tmp_dir_name = None
            if self.comm.rank == 0:
                tmp_dir = tempfile.TemporaryDirectory()
                tmp_dir_name = tmp_dir.name
            # Broadcast temp directory name to other procs
            tmp_dir_name = self.comm.bcast(tmp_dir_name, root=0)

            # Loop through each problem
            for prob in self.tacs_probs:
                # Solve problem
                prob.solve()
                # Write solution
                prob.writeSolution(outputDir=tmp_dir_name)
                if isinstance(prob, problems.StaticProblem):
                    prob.writeSolutionHistory(outputDir=tmp_dir_name)

            if self.comm.rank == 0:
                # Loop through each problem and make sure solution file exists and history file exists if static nonlinear problem
                for prob in self.tacs_probs:
                    with self.subTest(problem=prob.name):
                        base_name = os.path.join(tmp_dir_name, f"{prob.name}_000")
                        if isinstance(prob, problems.StaticProblem):
                            f5_file = f"{base_name}.f5"
                            self.assertTrue(
                                os.path.exists(f5_file), msg=f"{f5_file} not found"
                            )
                            history_file = f"{base_name}.pkl"
                            if prob.isNonlinear:
                                self.assertTrue(
                                    os.path.exists(history_file),
                                    msg=f"{history_file} not found",
                                )
                            else:
                                self.assertFalse(
                                    os.path.exists(history_file),
                                    msg=f"{history_file} found but should not have been",
                                )
                        else:
                            if isinstance(prob, problems.TransientProblem):
                                num_steps = prob.getNumTimeSteps() + 1
                            else:  # ModalProblem or BucklingProblem
                                num_steps = prob.getNumEigs()
                            for i in range(num_steps):
                                f5_file = f"{base_name}_%3.3d.f5" % (i)
                                self.assertTrue(
                                    os.path.exists(f5_file), msg=f"{f5_file} exists"
                                )

                # delete all files in temp dir
                tmp_dir.cleanup()

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
                prob.setNodes(xpts)
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
                prob.setNodes(xpts)
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
