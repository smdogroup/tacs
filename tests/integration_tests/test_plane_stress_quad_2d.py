import numpy as np
from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions
import unittest

'''
Create a uniform plate under uniform plane stress 
and test KSFailure, StructuralMass, and Compliance functions and sensitivities
'''

FUNC_REFS = np.array([1.46051701859883, 25700.0, 1.75e+07])


class ProblemTest(unittest.TestCase):
    def setUp(self):
        self.dtype = TACS.dtype

        if self.dtype == complex:
            self.rtol = 1e-11
            self.atol = 1e-8
            self.dh = 1e-50
        else:
            self.rtol = 1e-2
            self.atol = 1e-4
            self.dh = 1e-5

        # Length of plate in x/y direction
        Lx = 10.0
        Ly = 10.0

        # Number of elements in x/y direction
        nx = 10
        ny = 10

        # running loads (N/m)
        Nx = 175e5
        Ny = 175e5
        Nxy = -175e5

        # Set the MPI communicator
        comm = MPI.COMM_WORLD

        # Create the stiffness object
        props = constitutive.MaterialProperties(rho=2570.0, E=70e9, nu=0.3, ys=350e6)
        stiff = constitutive.PlaneStressConstitutive(props, t=0.1, tNum=0)

        # Set up the basis function
        model = elements.LinearElasticity2D(stiff)
        basis = elements.LinearQuadBasis()
        elem = elements.Element2D(model, basis)

        # Allocate the TACSCreator object
        vars_per_node = model.getVarsPerNode()
        creator = TACS.Creator(comm, vars_per_node)

        if comm.rank == 0:
            num_elems = nx * ny
            num_nodes = (nx + 1) * (ny + 1)

            # discretize plate
            x = np.linspace(0, Lx, nx + 1, self.dtype)
            y = np.linspace(0, Ly, ny + 1, self.dtype)
            xyz = np.zeros([ny + 1, nx + 1, 3], self.dtype)
            xyz[:, :, 0], xyz[:, :, 1] = np.meshgrid(x, y)

            node_ids = np.arange(num_nodes).reshape(ny + 1, nx + 1)

            # Set connectivity for each element
            conn = []
            for j in range(nx):
                for i in range(ny):
                    conn.append([node_ids[i, j],
                                 node_ids[i + 1, j],
                                 node_ids[i, j + 1],
                                 node_ids[i + 1, j + 1]])

            conn = np.array(conn, dtype=np.intc).flatten()
            ptr = np.arange(0, 4 * num_elems + 1, 4, dtype=np.intc)
            comp_ids = np.zeros(num_elems, dtype=np.intc)
            self.num_components = 1

            creator.setGlobalConnectivity(num_nodes, ptr, conn, comp_ids)

            # Set up the boundary conditions (fixed at bottom left corner, roller on top right)
            bcnodes = np.array([node_ids[0, 0], node_ids[-1, 0]], dtype=np.intc)
            bcvars = np.array([0, 1, 0], dtype=np.intc)

            # Set the boundary condition pointers
            bcptr = np.array([0, 2, 3], dtype=np.intc)
            creator.setBoundaryConditions(bcnodes, bcptr, bcvars)

            # Set the node locations
            creator.setNodes(xyz.flatten())

        # Set the elements for each (only one) component
        element_list = [elem]
        creator.setElements(element_list)

        # Create the tacs assembler object
        self.assembler = creator.createTACS()

        # Get the design variable values
        self.dv0 = self.assembler.createDesignVec()
        self.assembler.getDesignVars(self.dv0)

        # Create tacs vectors and Matrix for linear/adjoint solve
        self.f = self.assembler.createVec()
        self.res0 = self.assembler.createVec()
        self.ans0 = self.assembler.createVec()
        self.mat = self.assembler.createSchurMat()

        # Jacobian matrix factors
        self.alpha = 1.0
        self.beta = 0.0
        self.gamma = 0.0

        # The nodes have been distributed across processors now
        # Let's find which nodes this processor owns
        self.xpts0 = self.assembler.createNodeVec()
        self.assembler.getNodes(self.xpts0)
        xpts0_array = self.xpts0.getArray()
        # Split node vector into numpy arrays for easier parsing of vectors
        local_num_nodes = len(xpts0_array) // 3
        local_xyz = xpts0_array.reshape(local_num_nodes, 3)
        local_x, local_y, local_z = local_xyz[:, 0], local_xyz[:, 1], local_xyz[:, 2]

        # Create temporary dv vec for doing fd/cs
        self.dv1 = self.assembler.createDesignVec()
        self.dv_pert = self.assembler.createDesignVec()
        dv_pert_array = self.dv_pert.getArray()
        dv_pert_array[:] = 1.0

        # Create temporary state variable vec for doing fd/cs
        self.ans1 = self.assembler.createVec()
        self.ans_pert = self.assembler.createVec()
        ans_pert_array = self.ans_pert.getArray()
        # Define perturbation array that uniformly moves all nodes on right edge of plate to the right
        ans_pert_array = ans_pert_array.reshape(local_num_nodes, vars_per_node)
        ans_pert_array[local_x == Lx, 0] = 1.0

        # Create temporary nodal vec for doing fd/cs
        self.xpts1 = self.assembler.createNodeVec()
        # Define perturbation array that uniformly moves all nodes on right edge of plate to the right
        self.xpts_pert = self.assembler.createNodeVec()
        xpts_pert_array = self.xpts_pert.getArray()
        xpts_pert_array = xpts_pert_array.reshape(local_num_nodes, 3)
        xpts_pert_array[local_x == Lx, 0] = 1.0

        # Create the preconditioner for the corresponding matrix
        self.pc = TACS.Pc(self.mat)
        # Create GMRES solver object
        subspace = 100
        restarts = 2
        self.gmres = TACS.KSM(self.mat, self.pc, subspace, restarts)

        f_array = self.f.getArray().reshape(num_nodes, vars_per_node)

        # Apply distributed forces on edges of plate
        # Apply Nxx
        f_array[local_x == Lx, 0] += (Nx * Ly) / ny
        f_array[local_x == 0.0, 0] += -(Nx * Ly) / ny

        # Apply Nyy
        f_array[local_y == Ly, 1] += (Ny * Lx) / nx
        f_array[local_y == 0.0, 1] += -(Ny * Lx) / nx

        # Apply Nxy
        f_array[local_y == Ly, 0] += (Nxy * Lx) / nx
        f_array[local_x == Lx, 1] += (Nxy * Ly) / ny
        f_array[local_y == 0.0, 0] += -(Nxy * Lx) / nx
        f_array[local_x == 0.0, 1] += -(Nxy * Ly) / ny

        # drop force at corners by half to avoid stress concentration
        f_array[np.logical_and(local_x == 0.0, local_y == 0.0), :] *= 0.5
        f_array[np.logical_and(local_x == 0.0, local_y == Ly), :] *= 0.5
        f_array[np.logical_and(local_x == Lx, local_y == Ly), :] *= 0.5
        f_array[np.logical_and(local_x == Lx, local_y == 0.0), :] *= 0.5

        # Create the function list
        # KS function weight
        ksweight = 10.0
        self.func_list = [functions.KSFailure(self.assembler, ksweight),
                          functions.StructuralMass(self.assembler),
                          functions.Compliance(self.assembler)]

        self.dfdu_list = []
        self.adjoint_list = []
        self.dfddv_list = []
        self.dfdx_list = []
        for i in range(len(self.func_list)):
            self.dfdu_list.append(self.assembler.createVec())
            self.adjoint_list.append(self.assembler.createVec())
            self.dfddv_list.append(self.assembler.createDesignVec())
            self.dfdx_list.append(self.assembler.createNodeVec())

    def test_solve(self):
        '''
        Test linear solve and function evaluations
        '''
        # Make sure vecs are initialized to zero
        self.zero_tacs_vecs()

        # solve
        func_vals = self.run_solve()

        # Test functions values against historical values
        np.testing.assert_allclose(func_vals, FUNC_REFS, rtol=self.rtol, atol=self.atol)

    def test_partial_dv_sensitivities(self):
        '''
        Test partial dv sensitivity against fd/cs
        '''
        # Make sure vecs are initialized to zero
        self.zero_tacs_vecs()

        # Initial solve
        func_vals = self.run_solve()

        # Compute the partial derivative w.r.t. material design variables
        self.assembler.addDVSens(self.func_list, self.dfddv_list, 1.0)

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
                np.testing.assert_allclose(dfddv_proj_i, fdv_sens_approx[i],
                                           rtol=self.rtol, atol=self.atol)

    def test_partial_xpt_sensitivities(self):
        '''
        Test partial xpt sensitivity against fd/cs
        '''
        # Make sure vecs are initialized to zero
        self.zero_tacs_vecs()

        # Initial solve
        func_vals = self.run_solve()

        # Compute the total derivative w.r.t. nodal xpt locations
        self.assembler.addXptSens(self.func_list, self.dfdx_list, 1.0)

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
                np.testing.assert_allclose(dfdx_proj_i, f_xpt_sens_approx[i], rtol=self.rtol, atol=self.atol)

    def test_partial_sv_sensitivities(self):
        '''
        Test partial sv sensitivity against fd/cs
        '''
        # Make sure vecs are initialized to zero
        self.zero_tacs_vecs()

        # Initial solve
        func_vals = self.run_solve()

        # Compute the partial derivative w.r.t. state variables
        self.assembler.addSVSens(self.func_list, self.dfdu_list, self.alpha, self.beta, self.gamma)

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
                np.testing.assert_allclose(dfdu_proj_i, f_u_sens_approx[i],
                                           rtol=self.rtol, atol=self.atol)

    def test_total_dv_sensitivities(self):
        '''
        Test total dv sensitivity through adjoint against fd/cs
        '''
        # Make sure vecs are initialized to zero
        self.zero_tacs_vecs()

        # Initial solve
        func_vals = self.run_solve()

        # Compute the total derivative w.r.t. material design variables using adjoint
        self.run_adjoints()
        self.assembler.addDVSens(self.func_list, self.dfddv_list, 1.0)
        self.assembler.addAdjointResProducts(self.adjoint_list, self.dfddv_list, -1.0)

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
                np.testing.assert_allclose(dfddv_proj_i, fdv_sens_approx[i],
                                           rtol=self.rtol, atol=self.atol)

    def test_total_xpt_sensitivities(self):
        '''
        Test total xpt sensitivity through adjoint against fd/cs
        '''
        # Make sure vecs are initialized to zero
        self.zero_tacs_vecs()

        # Initial solve
        func_vals = self.run_solve()

        # Compute the total derivative w.r.t. nodal xpt locations using adjoint
        self.run_adjoints()
        self.assembler.addXptSens(self.func_list, self.dfdx_list, 1.0)
        self.assembler.addAdjointResXptSensProducts(self.adjoint_list, self.dfdx_list, -1.0)

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
                np.testing.assert_allclose(dfdx_proj_i, f_xpt_sens_approx[i], rtol=self.rtol, atol=self.atol)

    def run_solve(self, dv=None, xpts=None):
        '''
        Run a linear solve at specified design point and return functions of interest
        '''
        if dv == None:
            dv = self.dv0

        if xpts == None:
            xpts = self.xpts0

        # Set the design variables
        self.assembler.setDesignVars(dv)

        # Set node locations
        self.assembler.setNodes(xpts)

        # Assemble the stiffness matrix
        self.assembler.zeroVariables()
        self.assembler.assembleJacobian(self.alpha, self.beta, self.gamma, self.res0, self.mat)
        self.pc.factor()

        # add force vector to residual (R = Ku - f)
        self.res0.axpy(-1.0, self.f)

        # zero out bc terms in res and solve
        self.assembler.applyBCs(self.res0)
        # Solve the linear system
        self.gmres.solve(self.res0, self.ans0)
        self.ans0.scale(-1.0)

        # Update state variables with solution
        self.assembler.setVariables(self.ans0)

        func_vals = self.assembler.evalFunctions(self.func_list)

        return np.array(func_vals)

    def run_adjoints(self):
        '''
        Run adjoint solves for each function of interest
        '''
        # Set the design variables
        self.assembler.setDesignVars(self.dv0)

        # Set node locations
        self.assembler.setNodes(self.xpts0)

        # Assemble the transpose stiffness matrix
        self.assembler.assembleJacobian(self.alpha, self.beta, self.gamma, None, self.mat, TACS.TRANSPOSE)
        self.pc.factor()

        # Solve for the adjoint variables
        self.assembler.addSVSens(self.func_list, self.dfdu_list, self.alpha, self.beta, self.gamma)
        for i in range(len(self.func_list)):
            self.gmres.solve(self.dfdu_list[i], self.adjoint_list[i])

    def zero_tacs_vecs(self):
        '''
        Reset all vectors associated with solution and adjoint
        '''
        # Zero solution vector
        self.ans0.zeroEntries()

        # Zero residual
        self.res0.zeroEntries()

        # Set state vars to zero
        self.assembler.setVariables(self.ans0)

        # Zero dv sens for each function
        for dfddv in self.dfddv_list:
            dfddv.zeroEntries()

        # Zero xpt sens for each function
        for dfdx in self.dfdx_list:
            dfdx.zeroEntries()

        # Zero sv sens for each function
        for dfdu in self.dfdu_list:
            dfdu.zeroEntries()

    def perturb_tacs_vec(self, vec_out, vec_in, vec_pert):
        '''
        Perform fd/cs perturbation on tacs vector as follows
        vec_out = vec_in + scale * vec_pert

        where:
            scale = dh * 1j, in complex mode
            scale = dh, in real mode
        '''
        vec_out.copyValues(vec_in)
        if self.dtype == complex:
            vec_out.axpy(self.dh * 1j, vec_pert)
        else:
            vec_out.axpy(self.dh, vec_pert)

    def compute_fdcs_approx(self, vec1, vec0):
        '''
        Perform fd/cs calculation to approximate sensitivities

        difference performed as follows:
            sens = imag(vec1) / dh, in complex mode
            sens = (vec1 - vec0) / dh, in real mode
        '''
        if self.dtype == complex:
            sens_approx = np.imag(vec1) / self.dh
        else:
            sens_approx = (vec1 - vec0) / self.dh
        return sens_approx
