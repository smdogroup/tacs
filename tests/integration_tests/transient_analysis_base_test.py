import numpy as np
from tacs import TACS
import unittest
from mpi4py import MPI

"""
This is a base class for transient problem unit test cases.
This base class will test function evaluations and total
sensitivities for the user-specified problem that inherits from it.
When the user creates a new test based on this class four 
methods are required to be defined in the child class. 

    1. setup_assembler
    2. setup_integrator
    3. setup_tacs_vecs
    4. setup_funcs
    
See the virtual method implementations for each method 
below for more details.

NOTE: The child class must NOT implement its own setUp method 
for the unittest class. This is handled in the base class.
"""


class TransientTestCase:
    class TransientTest(unittest.TestCase):
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

            # Initial nodal location vector
            self.xpts0 = self.assembler.createNodeVec()
            self.assembler.getNodes(self.xpts0)

            # Create temporary dv vec for doing fd/cs
            self.dv1 = self.assembler.createDesignVec()
            self.xpts1 = self.assembler.createNodeVec()

            # Setup user-defined integrator for this test
            self.integrator = self.setup_integrator(self.assembler)
            self.num_steps = self.integrator.getNumTimeSteps()

            # Setup force and perturbation vectors used for fd/cs projections
            self.f_list = [
                self.assembler.createVec() for i in range(self.num_steps + 1)
            ]
            self.dv_pert = self.assembler.createDesignVec()
            self.xpts_pert = self.assembler.createNodeVec()
            # Populate force and perturbation vectors based on user-defined method
            self.setup_tacs_vecs(
                self.assembler, self.f_list, self.dv_pert, self.xpts_pert
            )

            # Create the function list
            self.func_list, self.func_ref = self.setup_funcs(self.assembler)
            # Set functions for integrator
            self.integrator.setFunctions(self.func_list)

        def setup_assembler(self, comm, dtype):
            """
            Setup mesh and tacs assembler for problem we will be testing.
            Must be defined in child class that inherits from this class.
            """
            raise NotImplementedError(
                "Child class must implement a 'setup_assembler' method"
            )
            return

        def setup_integrator(self, assembler):
            """
            Setup tacs integrator responsible for solving transient problem we will be testing.
            Must be defined in child class that inherits from this class.
            """
            raise NotImplementedError(
                "Child class must implement a 'setup_integrator' method"
            )
            return

        def setup_tacs_vecs(self, assembler, force_vec, dv_pert_vec, xpts_pert_vec):
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
            Test transient solve and function evaluations
            """
            # Make sure vecs are initialized to zero
            self.zero_tacs_vecs()

            # solve
            func_vals = self.run_solve()

            # Test functions values against historical values
            np.testing.assert_allclose(
                func_vals, self.func_ref, rtol=self.rtol, atol=self.atol
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

            # Compute the total derivative w.r.t. material design variables using fd/cs
            self.perturb_tacs_vec(self.dv1, self.dv0, self.dv_pert)
            # Run perturbed solution
            func_vals_pert = self.run_solve(dv=self.dv1)
            # Compute approximate sens
            fdv_sens_approx = self.compute_fdcs_approx(func_vals_pert, func_vals)

            # Tests cs/fd against sensitivity from adjoint
            for i in range(len(self.func_list)):
                with self.subTest(function=self.func_list[i]):
                    dfddv = self.integrator.getGradient(i)
                    dfddv_proj_i = dfddv.dot(self.dv_pert)
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

            # Compute the total derivative w.r.t. nodal xpt locations using fd/cs
            self.perturb_tacs_vec(self.xpts1, self.xpts0, self.xpts_pert)
            # Run perturbed solution
            func_vals_pert = self.run_solve(xpts=self.xpts1)
            # Compute approximate sens
            f_xpt_sens_approx = self.compute_fdcs_approx(func_vals_pert, func_vals)

            # Tests cs/fd against sensitivity from adjoint
            for i in range(len(self.func_list)):
                with self.subTest(function=self.func_list[i]):
                    dfdx = self.integrator.getXptGradient(i)
                    dfdx_proj_i = dfdx.dot(self.xpts_pert)
                    np.testing.assert_allclose(
                        dfdx_proj_i,
                        f_xpt_sens_approx[i],
                        rtol=self.rtol,
                        atol=self.atol,
                    )

        def run_solve(self, dv=None, xpts=None):
            """
            Run a transient solve at specified design point and return functions of interest
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

            # Loop over every time step and solve transient problem
            for i in range(self.num_steps + 1):
                self.integrator.iterate(i, forces=self.f_list[i])

            func_vals = self.integrator.evalFunctions(self.func_list)

            return np.array(func_vals)

        def run_adjoints(self):
            """
            Run adjoint solves for each function of interest
            """
            # Set the design variables
            self.assembler.setDesignVars(self.dv0)

            # Set node locations
            self.assembler.setNodes(self.xpts0)

            # Solve adjoint backwards in time
            self.integrator.integrateAdjoint()

        def zero_tacs_vecs(self):
            """
            Reset all vectors associated with solution and adjoint
            """

            # Set state vars to zero
            self.assembler.zeroVariables()

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
