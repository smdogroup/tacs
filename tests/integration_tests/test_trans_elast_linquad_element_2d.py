import numpy as np
from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions
from transient_analysis_base_test import TransientTestCase

"""
Create a uniform plate under uniform plane stress with sinusoidal loading in time
and test KSFailure, StructuralMass, and Compliance functions and sensitivities
"""

FUNC_REFS = np.array([1.333569238830022, 25700.0, 8740219.285287695, 4.061351738402616])

# Length of plate in x/y direction
Lx = 10.0
Ly = 10.0

# Number of elements in x/y direction
nx = 5
ny = 5

# running loads (N/m)
Nx = 175e5
Ny = 175e5
Nxy = -175e5

# Frequency of applied loads
fhz = 10.0

# Transient problem parameters
tinit = 0.0
tfinal = 1.0
num_steps = 50

# KS function weight
ksweight = 10.0


class ProblemTest(TransientTestCase.TransientTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    def setup_assembler(self, comm, dtype):
        """
        Setup mesh and tacs assembler for problem we will be testing.
        """

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
            x = np.linspace(0, Lx, nx + 1, dtype)
            y = np.linspace(0, Ly, ny + 1, dtype)
            xyz = np.zeros([nx + 1, ny + 1, 3], dtype)
            xyz[:, :, 0], xyz[:, :, 1] = np.meshgrid(x, y, indexing="ij")

            node_ids = np.arange(num_nodes).reshape(nx + 1, ny + 1)

            # Set connectivity for each element
            conn = []
            for i in range(nx):
                for j in range(ny):
                    conn.append(
                        [
                            node_ids[i, j],
                            node_ids[i + 1, j],
                            node_ids[i, j + 1],
                            node_ids[i + 1, j + 1],
                        ]
                    )

            conn = np.array(conn, dtype=np.intc).flatten()
            ptr = np.arange(0, 4 * num_elems + 1, 4, dtype=np.intc)
            comp_ids = np.zeros(num_elems, dtype=np.intc)

            creator.setGlobalConnectivity(num_nodes, ptr, conn, comp_ids)

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

    def setup_integrator(self, assembler):
        """
        Setup tacs integrator responsible for solving transient problem we will be testing.
        """
        # Create the BDF integrator solver
        order = 2

        # Set the file output format
        integrator = TACS.BDFIntegrator(
            assembler, tinit, tfinal, float(num_steps), order
        )

        return integrator

    def setup_tacs_vecs(self, assembler, force_history, dv_pert_vec, xpts_pert_vec):
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

        # Loop through the force vector for every time step and set the time-dependent load
        time_history = np.linspace(tinit, tfinal, num_steps + 1)
        for t, force_vec in zip(time_history, force_history):
            # Create force vector
            f_array = force_vec.getArray().reshape(local_num_nodes, vars_per_node)

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

            # Scale entire force magnitude by time-dependent scaling factor
            alpha = np.cos(2 * np.pi * fhz * t)
            force_vec.scale(alpha)

        # Create temporary dv vec for doing fd/cs
        dv_pert_array = dv_pert_vec.getArray()
        dv_pert_array[:] = 1.0

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
            functions.KSFailure(assembler, ksWeight=ksweight),
            functions.StructuralMass(assembler),
            functions.Compliance(assembler),
            functions.KSDisplacement(
                assembler, ksWeight=ksweight, direction=[100.0, 100.0]
            ),
        ]
        return func_list, FUNC_REFS
