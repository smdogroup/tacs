import numpy as np
from tacs import TACS
import unittest
from mpi4py import MPI

"""
This is a base class for static problem unit test cases.
This base class will test function evaluations and total 
and partial sensitivities for the user-specified problem 
that inherits from it.
When the user creates a new test based on this class three 
methods are required to be defined in the child class. 

    1. setup_assembler
    2. setup_tacs_vecs
    3. setup_funcs
    
See the virtual method implementations for each method 
below for more details.

NOTE: The child class must NOT implement its own setUp method 
for the unittest class. This is handled in the base class.
"""


class StaticTestCase:
    class StaticTest(unittest.TestCase):
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
            if not hasattr(self, "comm"):
                self.comm = MPI.COMM_WORLD

            # Setup user-specified assembler for this test
            self.assembler = self.setup_assembler(self.comm, self.dtype)

            # Get the design variable values
            self.dv0 = self.assembler.createDesignVec()
            self.assembler.getDesignVars(self.dv0)

            # Create tacs vectors and Matrix for linear/adjoint solve
            self.res0 = self.assembler.createVec()
            self.ans0 = self.assembler.createVec()
            self.mat = self.assembler.createSchurMat()

            # Jacobian matrix factors
            self.alpha = 1.0
            self.beta = 0.0
            self.gamma = 0.0

            # Initial nodal location vector
            self.xpts0 = self.assembler.createNodeVec()
            self.assembler.getNodes(self.xpts0)

            # Create temporary dv vec for doing fd/cs
            self.dv1 = self.assembler.createDesignVec()
            self.ans1 = self.assembler.createVec()
            self.xpts1 = self.assembler.createNodeVec()

            # Setup force and perturbation vectors used for fd/cs projections
            self.f = self.assembler.createVec()
            self.dv_pert = self.assembler.createDesignVec()
            self.ans_pert = self.assembler.createVec()
            self.xpts_pert = self.assembler.createNodeVec()
            # Populate force and perturbation vectors based on user-defined method
            self.setup_tacs_vecs(
                self.assembler, self.f, self.dv_pert, self.ans_pert, self.xpts_pert
            )

            # Zero out any bc nodes in the state variable vec (if the user didn't already do this)
            self.assembler.applyBCs(self.ans_pert)

            # Create the preconditioner for the corresponding matrix
            self.pc = TACS.Pc(self.mat)
            # Create GMRES solver object
            subspace = 100
            restarts = 2
            atol = 1e-30
            rtol = 1e-12
            self.gmres = TACS.KSM(self.mat, self.pc, subspace, restarts)
            self.gmres.setTolerances(rtol, atol)

            # Create the function list
            self.func_list, self.func_ref = self.setup_funcs(self.assembler)

            self.dfdu_list = []
            self.adjoint_list = []
            self.dfddv_list = []
            self.dfdx_list = []
            for i in range(len(self.func_list)):
                self.dfdu_list.append(self.assembler.createVec())
                self.adjoint_list.append(self.assembler.createVec())
                self.dfddv_list.append(self.assembler.createDesignVec())
                self.dfdx_list.append(self.assembler.createNodeVec())

        def setup_assembler(self, comm, dtype):
            """
            Setup mesh and tacs assembler for problem we will be testing.
            Must be defined in child class that inherits from this class.
            """
            raise NotImplementedError(
                "Child class must implement a 'setup_assembler' method"
            )
            return

        def setup_tacs_vecs(
            self, assembler, force_vec, dv_pert_vec, ans_pert_vec, xpts_pert_vec
        ):
            """
            Setup user-defined vectors for analysis and fd/cs sensitivity verification.
            Must be defined in child class that inherits from this class.
            """
            raise NotImplementedError(
                "Child class must implement a 'setup_tacs_vecs' method"
            )
            return

        def setup_funcs(self, assembler):
            """
            Create a list of functions to be tested and their reference values for the problem.
            Must be defined in child class that inherits from this class.
            """
            raise NotImplementedError(
                "Child class must implement a 'setup_funcs' method"
            )
            return

        def test_solve(self):
            """
            Test linear solve and function evaluations
            """
            # Make sure vecs are initialized to zero
            self.zero_tacs_vecs()

            # solve
            func_vals = self.run_solve()

            # Test functions values against historical values
            np.testing.assert_allclose(
                func_vals, self.func_ref, rtol=self.rtol, atol=self.atol
            )

        def test_partial_dv_sensitivities(self):
            """
            Test partial dv sensitivity against fd/cs
            """
            # Make sure vecs are initialized to zero
            self.zero_tacs_vecs()

            # Initial solve
            func_vals = self.run_solve()

            # Compute the partial derivative w.r.t. material design variables
            self.assembler.addDVSens(self.func_list, self.dfddv_list, 1.0)
            # Accumulate sensitivity across all procs
            self.set_tacs_vec_values(self.dfddv_list)

            # Compute the total derivative w.r.t. material design variables using fd/cs
            self.perturb_tacs_vec(self.dv1, self.dv0, self.dv_pert)
            # Set the perturbed design variables
            self.assembler.setDesignVars(self.dv1)
            # Compute functions w/o resolving problem
            func_vals_pert = self.assembler.evalFunctions(self.func_list)
            # Compute approximate sens
            fdv_sens_approx = self.compute_fdcs_approx(func_vals_pert, func_vals)

            # Tests cs/fd against sensitivity from partial
            for i in range(len(self.func_list)):
                with self.subTest(function=self.func_list[i]):
                    dfddv_proj_i = self.dfddv_list[i].dot(self.dv_pert)
                    np.testing.assert_allclose(
                        dfddv_proj_i, fdv_sens_approx[i], rtol=self.rtol, atol=self.atol
                    )

        def test_partial_xpt_sensitivities(self):
            """
            Test partial xpt sensitivity against fd/cs
            """
            # Make sure vecs are initialized to zero
            self.zero_tacs_vecs()

            # Initial solve
            func_vals = self.run_solve()

            # Compute the total derivative w.r.t. nodal xpt locations
            self.assembler.addXptSens(self.func_list, self.dfdx_list, 1.0)
            # Accumulate sensitivity across all procs
            self.set_tacs_vec_values(self.dfdx_list)

            # Compute the total derivative w.r.t. nodal xpt locations using fd/cs
            self.perturb_tacs_vec(self.xpts1, self.xpts0, self.xpts_pert)
            # Set the perturbed node locations
            self.assembler.setNodes(self.xpts1)
            # Compute functions w/o resolving problem
            func_vals_pert = self.assembler.evalFunctions(self.func_list)
            # Compute approximate sens
            f_xpt_sens_approx = self.compute_fdcs_approx(func_vals_pert, func_vals)

            # Tests cs/fd against sensitivity from partial
            for i in range(len(self.func_list)):
                with self.subTest(function=self.func_list[i]):
                    dfdx_proj_i = self.dfdx_list[i].dot(self.xpts_pert)
                    np.testing.assert_allclose(
                        dfdx_proj_i,
                        f_xpt_sens_approx[i],
                        rtol=self.rtol,
                        atol=self.atol,
                    )

        def test_partial_sv_sensitivities(self):
            """
            Test partial sv sensitivity against fd/cs
            """
            # Make sure vecs are initialized to zero
            self.zero_tacs_vecs()

            # Initial solve
            func_vals = self.run_solve()

            # Compute the partial derivative w.r.t. state variables
            self.assembler.addSVSens(
                self.func_list, self.dfdu_list, self.alpha, self.beta, self.gamma
            )

            # Compute the total derivative w.r.t. material design variables using fd/cs
            self.perturb_tacs_vec(self.ans1, self.ans0, self.ans_pert)
            # Set the perturbed state variables
            self.assembler.setVariables(self.ans1)
            # Compute functions w/o resolving problem
            func_vals_pert = self.assembler.evalFunctions(self.func_list)
            # Compute approximate sens
            f_u_sens_approx = self.compute_fdcs_approx(func_vals_pert, func_vals)

            # Tests cs/fd against sensitivity from partial
            for i in range(len(self.func_list)):
                with self.subTest(function=self.func_list[i]):
                    dfdu_proj_i = self.dfdu_list[i].dot(self.ans_pert)
                    np.testing.assert_allclose(
                        dfdu_proj_i, f_u_sens_approx[i], rtol=self.rtol, atol=self.atol
                    )

        def test_total_dv_sensitivities(self):
            """
            Test total dv sensitivity through adjoint against fd/cs
            """
            # Make sure vecs are initialized to zero
            self.zero_tacs_vecs()

            # Initial solve
            func_vals = self.run_solve()

            # Compute the total derivative w.r.t. material design variables using adjoint
            self.run_adjoints()
            self.assembler.addDVSens(self.func_list, self.dfddv_list, 1.0)
            self.assembler.addAdjointResProducts(
                self.adjoint_list, self.dfddv_list, -1.0
            )
            # Accumulate sensitivity across all procs
            self.set_tacs_vec_values(self.dfddv_list)

            # Compute the total derivative w.r.t. material design variables using fd/cs
            self.perturb_tacs_vec(self.dv1, self.dv0, self.dv_pert)
            # Run perturbed solution
            func_vals_pert = self.run_solve(dv=self.dv1)
            # Compute approximate sens
            fdv_sens_approx = self.compute_fdcs_approx(func_vals_pert, func_vals)

            # Tests cs/fd against sensitivity from adjoint
            for i in range(len(self.func_list)):
                with self.subTest(function=self.func_list[i]):
                    dfddv_proj_i = self.dfddv_list[i].dot(self.dv_pert)
                    np.testing.assert_allclose(
                        dfddv_proj_i, fdv_sens_approx[i], rtol=self.rtol, atol=self.atol
                    )

        def test_total_xpt_sensitivities(self):
            """
            Test total xpt sensitivity through adjoint against fd/cs
            """
            # Make sure vecs are initialized to zero
            self.zero_tacs_vecs()

            # Initial solve
            func_vals = self.run_solve()

            # Compute the total derivative w.r.t. nodal xpt locations using adjoint
            self.run_adjoints()
            self.assembler.addXptSens(self.func_list, self.dfdx_list, 1.0)
            self.assembler.addAdjointResXptSensProducts(
                self.adjoint_list, self.dfdx_list, -1.0
            )
            # Accumulate sensitivity across all procs
            self.set_tacs_vec_values(self.dfdx_list)

            # Compute the total derivative w.r.t. nodal xpt locations using fd/cs
            self.perturb_tacs_vec(self.xpts1, self.xpts0, self.xpts_pert)
            # Run perturbed solution
            func_vals_pert = self.run_solve(xpts=self.xpts1)
            # Compute approximate sens
            f_xpt_sens_approx = self.compute_fdcs_approx(func_vals_pert, func_vals)

            # Tests cs/fd against sensitivity from adjoint
            for i in range(len(self.func_list)):
                with self.subTest(function=self.func_list[i]):
                    dfdx_proj_i = self.dfdx_list[i].dot(self.xpts_pert)
                    np.testing.assert_allclose(
                        dfdx_proj_i,
                        f_xpt_sens_approx[i],
                        rtol=self.rtol,
                        atol=self.atol,
                    )

        def run_solve(self, dv=None, xpts=None):
            """
            Run a linear solve at specified design point and return functions of interest
            """
            if dv is None:
                dv = self.dv0

            if xpts is None:
                xpts = self.xpts0

            # Set the design variables
            self.assembler.setDesignVars(dv)

            # Set node locations
            self.assembler.setNodes(xpts)

            # Assemble the stiffness matrix
            self.assembler.zeroVariables()
            self.assembler.assembleJacobian(
                self.alpha, self.beta, self.gamma, self.res0, self.mat
            )
            self.pc.factor()

            # zero out bc terms in force
            self.assembler.applyBCs(self.f)
            # add force vector to residual (R = Ku - f)
            self.res0.axpy(-1.0, self.f)

            # Solve the linear system
            self.gmres.solve(self.res0, self.ans0)
            self.ans0.scale(-1.0)

            # Update state variables with solution
            self.assembler.setVariables(self.ans0)

            func_vals = self.assembler.evalFunctions(self.func_list)

            return np.array(func_vals)

        def run_adjoints(self):
            """
            Run adjoint solves for each function of interest
            """
            # Set the design variables
            self.assembler.setDesignVars(self.dv0)

            # Set node locations
            self.assembler.setNodes(self.xpts0)

            # Assemble the transpose stiffness matrix
            self.assembler.assembleJacobian(
                self.alpha, self.beta, self.gamma, None, self.mat, TACS.TRANSPOSE
            )
            self.pc.factor()

            # Solve for the adjoint variables
            self.assembler.addSVSens(
                self.func_list, self.dfdu_list, self.alpha, self.beta, self.gamma
            )
            for i in range(len(self.func_list)):
                self.gmres.solve(self.dfdu_list[i], self.adjoint_list[i])

        def zero_tacs_vecs(self):
            """
            Reset all vectors associated with solution and adjoint
            """
            # Zero solution vector
            self.ans0.zeroEntries()

            # Zero residual
            self.res0.zeroEntries()

            # Set state vars to zero
            self.assembler.zeroVariables()

            # Zero dv sens for each function
            for dfddv in self.dfddv_list:
                dfddv.zeroEntries()

            # Zero xpt sens for each function
            for dfdx in self.dfdx_list:
                dfdx.zeroEntries()

            # Zero sv sens for each function
            for dfdu in self.dfdu_list:
                dfdu.zeroEntries()

        def set_tacs_vec_values(self, tacs_vecs):
            """
            Begin setting the values: Collective on the TACS communicator
            """
            for vec in tacs_vecs:
                vec.beginSetValues()
                vec.endSetValues()

        def perturb_tacs_vec(self, vec_out, vec_in, vec_pert):
            """
            Perform fd/cs perturbation on tacs vector as follows
            vec_out = vec_in + scale * vec_pert

            where:
                scale = dh * 1j, in complex mode
                scale = dh, in real mode
            """
            vec_out.copyValues(vec_in)
            if self.dtype == complex:
                vec_out.axpy(self.dh * 1j, vec_pert)
            else:
                vec_out.axpy(self.dh, vec_pert)

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
