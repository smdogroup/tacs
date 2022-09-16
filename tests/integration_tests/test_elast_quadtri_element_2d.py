import numpy as np
from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions
from static_analysis_base_test import StaticTestCase

"""
Create a uniform plate of quadratic triangles under uniform edge loading
and test KSFailure, StructuralMass, and Compliance functions and sensitivities
"""

FUNC_REFS = np.array([2.007845956565272, 25700.0, 17567887.317833334])

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

# KS function weight
ksweight = 10.0


class ProblemTest(StaticTestCase.StaticTest):

    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    def setup_assembler(self, comm, dtype):
        """
        Setup mesh and tacs assembler for problem we will be testing.
        """

        # Overwrite default tolerances from base class
        if dtype == complex:
            self.rtol = 1e-11
            self.atol = 1e-8
            self.dh = 1e-50
        else:
            self.rtol = 1e-2
            self.atol = 1e-4
            self.dh = 1e-6

        # Create the stiffness object
        props = constitutive.MaterialProperties(rho=2570.0, E=70e9, nu=0.3, ys=350e6)
        stiff = constitutive.PlaneStressConstitutive(props, t=0.1, tNum=0)

        # Set up the basis function
        model = elements.LinearElasticity2D(stiff)
        basis = elements.QuadraticTriangleBasis()
        elem = elements.Element2D(model, basis)

        # Allocate the TACSCreator object
        vars_per_node = model.getVarsPerNode()
        creator = TACS.Creator(comm, vars_per_node)

        if comm.rank == 0:

            # Set the nodes
            num_nodes = (2 * nx + 1) * (2 * ny + 1)
            num_elems = 2 * nx * ny

            # discretize plate
            x = np.linspace(0, Lx, 2 * nx + 1, dtype)
            y = np.linspace(0, Ly, 2 * ny + 1, dtype)
            xyz = np.zeros([2 * nx + 1, 2 * ny + 1, 3], dtype)
            xyz[:, :, 0], xyz[:, :, 1] = np.meshgrid(x, y, indexing="ij")
            node_ids = np.arange(num_nodes).reshape((2 * nx + 1, 2 * ny + 1))

            conn = []
            for j in range(ny):
                for i in range(nx):
                    # Append the first set of nodes
                    conn.append(
                        [
                            node_ids[2 * i, 2 * j],
                            node_ids[2 * i + 2, 2 * j],
                            node_ids[2 * i + 2, 2 * j + 2],
                            node_ids[2 * i + 1, 2 * j],
                            node_ids[2 * i + 2, 2 * j + 1],
                            node_ids[2 * i + 1, 2 * j + 1],
                        ]
                    )

                    # Append the second set of nodes
                    conn.append(
                        [
                            node_ids[2 * i, 2 * j + 2],
                            node_ids[2 * i, 2 * j],
                            node_ids[2 * i + 2, 2 * j + 2],
                            node_ids[2 * i, 2 * j + 1],
                            node_ids[2 * i + 1, 2 * j + 1],
                            node_ids[2 * i + 1, 2 * j + 2],
                        ]
                    )
            # Set the node pointers
            conn = np.array(conn, dtype=np.intc).flatten()
            ptr = np.arange(0, 6 * num_elems + 1, 6, dtype=np.intc)
            elem_ids = np.zeros(num_elems, dtype=np.intc)
            creator.setGlobalConnectivity(num_nodes, ptr, conn, elem_ids)

            # Set up the boundary conditions (fixed at bottom left corner, roller on top right)
            bcnodes = np.array([node_ids[0, 0], node_ids[0, -1]], dtype=np.intc)
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
        assembler = creator.createTACS()

        return assembler

    def setup_tacs_vecs(
        self, assembler, force_vec, dv_pert_vec, ans_pert_vec, xpts_pert_vec
    ):
        """
        Setup user-defined vectors for analysis and fd/cs sensitivity verification
        """
        local_num_nodes = assembler.getNumOwnedNodes()
        vars_per_node = assembler.getVarsPerNode()

        # The nodes have been distributed across processors now
        # Let's find which nodes this processor owns
        xpts0 = assembler.createNodeVec()
        assembler.getNodes(xpts0)
        xpts0_array = xpts0.getArray()
        # Split node vector into numpy arrays for easier parsing of vectors
        local_xyz = xpts0_array.reshape(local_num_nodes, 3)
        local_x, local_y, local_z = local_xyz[:, 0], local_xyz[:, 1], local_xyz[:, 2]

        # Create force vector
        f_array = force_vec.getArray().reshape(local_num_nodes, vars_per_node)

        # Apply distributed forces on edges of plate
        # Apply Nxx
        f_array[local_x == Lx, 0] += (Nx * Ly) / (2 * ny)
        f_array[local_x == 0.0, 0] += -(Nx * Ly) / (2 * ny)

        # Apply Nyy
        f_array[local_y == Ly, 1] += (Ny * Lx) / (2 * nx)
        f_array[local_y == 0.0, 1] += -(Ny * Lx) / (2 * nx)

        # Apply Nxy
        f_array[local_y == Ly, 0] += (Nxy * Lx) / (2 * nx)
        f_array[local_x == Lx, 1] += (Nxy * Ly) / (2 * ny)
        f_array[local_y == 0.0, 0] += -(Nxy * Lx) / (2 * nx)
        f_array[local_x == 0.0, 1] += -(Nxy * Ly) / (2 * ny)

        # drop force at corners by half to avoid stress concentration
        f_array[np.logical_and(local_x == 0.0, local_y == 0.0), :] *= 0.5
        f_array[np.logical_and(local_x == 0.0, local_y == Ly), :] *= 0.5
        f_array[np.logical_and(local_x == Lx, local_y == Ly), :] *= 0.5
        f_array[np.logical_and(local_x == Lx, local_y == 0.0), :] *= 0.5

        # Create temporary dv vec for doing fd/cs
        dv_pert_array = dv_pert_vec.getArray()
        dv_pert_array[:] = 1.0

        # Create temporary state variable vec for doing fd/cs
        ans_pert_array = ans_pert_vec.getArray()
        # Define perturbation array that uniformly moves all nodes on right edge of plate to the right
        ans_pert_array = ans_pert_array.reshape(local_num_nodes, vars_per_node)
        ans_pert_array[local_x == Lx, 0] = 1.0

        # Define perturbation array that uniformly moves all nodes on right edge of plate to the right
        xpts_pert_array = xpts_pert_vec.getArray()
        xpts_pert_array = xpts_pert_array.reshape(local_num_nodes, 3)
        # Define perturbation array that uniformly moves all nodes on right edge of plate to the right
        xpts_pert_array[local_x == Lx, 0] = 1.0

        return

    def setup_funcs(self, assembler):
        """
        Create a list of functions to be tested and their reference values for the problem
        """
        func_list = [
            functions.KSFailure(assembler, ksWeight=ksweight, safetyFactor=1.5),
            functions.StructuralMass(assembler),
            functions.Compliance(assembler),
        ]
        return func_list, FUNC_REFS
