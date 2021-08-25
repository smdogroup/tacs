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
            self.dh = 1e-4

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
        self.dv0_array = self.dv0.getArray()

        # Create temporary dv vec for doing fd/cs
        self.dv_pert = self.assembler.createDesignVec()
        self.dv_pert_array = self.dv_pert.getArray()

        # Create tacs vectors and Matrix
        self.f = self.assembler.createVec()
        self.res = self.assembler.createVec()
        self.ans = self.assembler.createVec()
        self.mat = self.assembler.createSchurMat()

        # Jacobian matrix factors
        self.alpha = 1.0
        self.beta = 0.0
        self.gamma = 0.0

        # The nodes have been distributed across processors now
        # Let's find which nodes this processor owns
        self.xpts0 = self.assembler.createNodeVec()
        self.assembler.getNodes(self.xpts0)
        self.xpts0_array = self.xpts0.getArray()
        # Split node vector into numpy arrays for easier parsing of vectors
        local_num_nodes = len(self.xpts0_array) // 3
        local_xyz = self.xpts0_array.reshape(local_num_nodes, 3)
        local_x, local_y, local_z = local_xyz[:, 0], local_xyz[:, 1], local_xyz[:, 2]

        # Create temporary nodal vec for doing fd/cs
        self.xpts_pert = self.assembler.createNodeVec()
        # Define perturbation array that uniformly moves all nodes on right edge of plate to the right
        self.pert = self.assembler.createNodeVec()
        pert_array = self.pert.getArray()
        pert_array = pert_array.reshape(local_num_nodes, 3)
        pert_array[local_x == Lx, 0] = 1.0

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
        func_vals = self.run_solve()
        np.testing.assert_allclose(func_vals, FUNC_REFS, rtol=self.rtol, atol=self.atol)

    def test_dv_sensitivities(self):
        '''
        Test total dv sensitivity through adjoint against fd/cs
        '''
        # Initial solve
        func_vals = self.run_solve()

        # Compute the total derivative w.r.t. material design variables using adjoint
        self.run_adjoints()
        self.assembler.addDVSens(self.func_list, self.dfddv_list, 1.0)
        self.assembler.addAdjointResProducts(self.adjoint_list, self.dfddv_list, -1.0)

        # Compute the total derivative w.r.t. material design variables using fd/cs
        if self.dtype == complex:
            self.dv_pert_array[:] = self.dv0_array[:] + self.dh * 1j
        else:
            self.dv_pert_array[:] = self.dv0_array[:] + self.dh
        func_vals_pert = self.run_solve(dv=self.dv_pert)
        if self.dtype == complex:
            fdvSens_approx = np.imag(func_vals_pert) / self.dh
        else:
            fdvSens_approx = (func_vals_pert - func_vals) / self.dh

        # Tests cs/fd against sensitivity from adjoint
        for i in range(len(self.func_list)):
            with self.subTest(function=self.func_list[i]):
                np.testing.assert_allclose(self.dfddv_list[i].getArray(), fdvSens_approx[i],
                                           rtol=self.rtol, atol=self.atol)

    def test_xpt_sensitivities(self):
        '''
        Test total xpt sensitivity through adjoint against fd/cs
        '''
        # Initial solve
        func_vals = self.run_solve()

        # Compute the total derivative w.r.t. nodal xpt locations using adjoint
        self.run_adjoints()
        self.assembler.addXptSens(self.func_list, self.dfdx_list, 1.0)
        self.assembler.addAdjointResXptSensProducts(self.adjoint_list, self.dfdx_list, -1.0)

        # Compute the total derivative w.r.t. nodal xpt locations using fd/cs
        self.xpts_pert.copyValues(self.xpts0)
        if self.dtype == complex:
            self.xpts_pert.axpy(self.dh * 1j, self.pert)
        else:
            self.xpts_pert.axpy(self.dh, self.pert)
        func_vals_pert = self.run_solve(xpts=self.xpts_pert)
        if self.dtype == complex:
            f_xpt_sens_approx = np.imag(func_vals_pert) / self.dh
        else:
            f_xpt_sens_approx = (func_vals_pert - func_vals) / self.dh

        # Tests cs/fd against sensitivity from adjoint
        for i in range(len(self.func_list)):
            with self.subTest(function=self.func_list[i]):
                dfdx_proj_i = self.dfdx_list[i].dot(self.pert)
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
        self.assembler.assembleJacobian(self.alpha, self.beta, self.gamma, self.res, self.mat)
        self.pc.factor()

        # add force vector to residual (R = Ku - f)
        self.res.axpy(-1.0, self.f)

        # zero out bc terms in res and solve
        self.assembler.applyBCs(self.res)
        # Solve the linear system
        self.gmres.solve(self.res, self.ans)
        self.ans.scale(-1.0)

        # Update state variables with solution
        self.assembler.setVariables(self.ans)

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
