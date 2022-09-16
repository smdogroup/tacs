import numpy as np
from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions
from static_analysis_base_test import StaticTestCase

"""
Create a uniform plate under distributed inertial loads due to gravity
and test KSFailure, StructuralMass, and Compliance functions and sensitivities
"""

FUNC_REFS = np.array(
    [0.9516087233995378, 25700.0, 2919804.598684517, 2.0769491952311663]
)

# Length of plate in x/y direction
Lx = 10.0
Ly = 10.0

# Number of elements in x/y direction
nx = 10
ny = 10

# 1000 G gravity vector (m/s^2)
g = np.array([0.0, -9810.0], dtype=TACS.dtype)

# KS function weight
ksweight = 10.0


class ProblemTest(StaticTestCase.StaticTest):
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
            x, y, z = xyz[:, :, 0], xyz[:, :, 1], xyz[:, :, 2]

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

            # Set up the boundary conditions ( clamp at y == 0  edge)
            bcnodes = np.array(node_ids[y == 0.0], dtype=np.intc)
            creator.setBoundaryConditions(bcnodes)

            # Set the node locations
            creator.setNodes(xyz.flatten())

        # Set the elements for each (only one) component
        element_list = [elem]
        creator.setElements(element_list)

        # Create the tacs assembler object
        assembler = creator.createTACS()

        # Use traction elements to set running loads
        aux_elems = TACS.AuxElements()

        # Add gravity load to all elements
        local_elements = assembler.getElements()
        for elem_id, elem in enumerate(local_elements):
            grav = elem.createElementInertialForce(g)
            aux_elems.addElement(elem_id, grav)

        # Set forces in assembler
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

        # No need to populate force vector, since we're using inertial forces

        # Create temporary dv vec for doing fd/cs
        dv_pert_array = dv_pert_vec.getArray()
        dv_pert_array[:] = 1.0

        # Create temporary state variable vec for doing fd/cs
        ans_pert_array = ans_pert_vec.getArray()
        # Define perturbation array that uniformly moves all nodes on right edge of plate to the right
        ans_pert_array = ans_pert_array.reshape(local_num_nodes, vars_per_node)
        ans_pert_array[local_y == 0.0, 0] = 1.0

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
                assembler, ksWeight=ksweight, direction=[-100.0, -100.0]
            ),
        ]
        return func_list, FUNC_REFS
