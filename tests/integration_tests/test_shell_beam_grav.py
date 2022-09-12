import numpy as np
from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions
from static_analysis_base_test import StaticTestCase

"""
Create a cantilevered beam of linear triangular shells under a gravity load
and test KSFailure, StructuralMass, and Compliance functions and sensitivities
"""

FUNC_REFS = np.array([3.26470740e-01, 2.57000000e03, 5.26220416e03, 1.15263745e-01])

# Length of plate in x/y direction
Lx = 10.0
Ly = 1.0

# Number of elements in x/y direction
nx = 10
ny = 10

# gravity vector (m/s^2)
g = np.array([0.0, 0.0, -9.81], dtype=TACS.dtype)

# KS function weight
ksweight = 10.0


class ProblemTest(StaticTestCase.StaticTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    def setup_assembler(self, comm, dtype):
        """
        Setup mesh and tacs assembler for problem we will be testing.
        """

        # Overwrite default check values
        if dtype == complex:
            self.rtol = 5e-8
            self.dh = 1e-50
        else:
            self.rtol = 1e-1
            self.dh = 1e-7

        # Only check for relative tolerance
        self.atol = 1e99

        # Create the stiffness object
        props = constitutive.MaterialProperties(rho=2570.0, E=70e9, nu=0.3, ys=350e6)
        stiff = constitutive.IsoShellConstitutive(props, t=0.1, tNum=0)

        # Set up the element transform
        ref_axis = np.array([0.0, 1.0, 1.0], dtype=dtype)
        transform = elements.ShellRefAxisTransform(ref_axis)
        elem = elements.Tri3Shell(transform, stiff)

        # Allocate the TACSCreator object
        vars_per_node = elem.getVarsPerNode()
        creator = TACS.Creator(comm, vars_per_node)

        if comm.rank == 0:
            num_elems = 2 * nx * ny
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
                        [node_ids[i, j], node_ids[i + 1, j], node_ids[i + 1, j + 1]]
                    )
                    conn.append(
                        [node_ids[i + 1, j + 1], node_ids[i, j + 1], node_ids[i, j]]
                    )

            conn = np.array(conn, dtype=np.intc).flatten()
            ptr = np.arange(0, 3 * num_elems + 1, 3, dtype=np.intc)
            comp_ids = np.zeros(num_elems, dtype=np.intc)

            creator.setGlobalConnectivity(num_nodes, ptr, conn, comp_ids)

            # Set up the boundary conditions (fixed at left hand edge)
            bcnodes = np.array(node_ids[0, :], dtype=np.intc)
            creator.setBoundaryConditions(bcnodes)

            # Set the node locations
            creator.setNodes(xyz.flatten())

        # Set the elements for each (only one) component
        element_list = [elem]
        creator.setElements(element_list)

        # Create the tacs assembler object
        assembler = creator.createTACS()

        # Get number of elements on this processor
        local_num_elems = assembler.getNumElements()

        # Create object to hold pressures
        aux_elems = TACS.AuxElements()

        # Add gravity load to all elements
        grav = elem.createElementInertialForce(g)
        for elem_id in range(local_num_elems):
            aux_elems.addElement(elem_id, grav)

        # Set tractions in assembler
        assembler.setAuxElements(aux_elems)

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

        # Don't need to modify force vector, since we already set pressure

        # Create temporary dv vec for doing fd/cs
        dv_pert_array = dv_pert_vec.getArray()
        dv_pert_array[:] = 1.0

        # Create temporary state variable vec for doing fd/cs
        ans_pert_array = ans_pert_vec.getArray()
        # Define perturbation array that uniformly moves all nodes on right edge of plate to the right
        ans_pert_array = ans_pert_array.reshape(local_num_nodes, vars_per_node)
        ans_pert_array[local_x == Lx, 2] = 1.0

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
                assembler, ksWeight=ksweight, direction=[0.0, 0.0, 1.0]
            ),
        ]
        return func_list, FUNC_REFS
