import numpy as np
from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions
from static_analysis_base_test import StaticTestCase

r"""
Create a two separate cantilevered plates connected by an RBE3 element.
Apply a load at the RBE3 center node and test KSFailure, StructuralMass, 
and Compliance functions and sensitivities
-----------        ----------- 
|          |\    /|          |
|          | \  / |          |
| Plate 1  |__\/__| Plate 2  |
|          |  /\  |          |
|          | /  \ |          |
|          |/    \|          |
------------       -----------
"""

FUNC_REFS = np.array(
    [1.2205891205367805, 51400.0, 3368332.5161940744, 2.748016507122232]
)

# Length of plate in x/y direction
Lx = 10.0
Ly = 10.0

# Number of elements in x/y direction for each plate
nx = 4
ny = 4

# applied force at center node
applied_force = np.array([1e8, 0.0, 1.0e6, 0.0, 0.0, 1e8])

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
            self.rtol = 1e-8
            self.atol = 1e-8
            self.dh = 1e-50
        else:
            self.rtol = 1e-1
            self.atol = 1e-4
            self.dh = 1e-7

        # Create the stiffness object
        props = constitutive.MaterialProperties(rho=2570.0, E=70e9, nu=0.3, ys=350e6)
        stiff = constitutive.IsoShellConstitutive(props, t=0.1, tNum=0)

        # Set up the element transform function
        transform = elements.ShellNaturalTransform()
        shell = elements.Quad4Shell(transform, stiff)
        num_rbe_nodes = 0

        # Allocate the TACSCreator object
        vars_per_node = shell.getVarsPerNode()
        creator = TACS.Creator(comm, vars_per_node)

        if comm.rank == 0:
            num_elems = nx * ny
            num_nodes = (nx + 1) * (ny + 1)

            # discretize (left) plate
            x = np.linspace(0, Lx, nx + 1, dtype)
            y = np.linspace(0, Ly, ny + 1, dtype)
            left_xyz = np.zeros([nx + 1, ny + 1, 3], dtype)
            left_xyz[:, :, 0], left_xyz[:, :, 1] = np.meshgrid(x, y, indexing="ij")

            left_node_ids = np.arange(num_nodes, dtype=np.intc).reshape(nx + 1, ny + 1)

            # Define right plate by copying left plate and shifting 2 m
            right_xyz = left_xyz.copy()
            right_xyz[:, :, 0] += 2.0 * Lx
            right_node_ids = left_node_ids + num_nodes

            # Double the node/element count
            num_nodes *= 2
            num_elems *= 2

            # Set connectivity for each plate element
            conn = []
            for i in range(nx):
                for j in range(ny):
                    conn.extend(
                        [
                            left_node_ids[i, j],
                            left_node_ids[i + 1, j],
                            left_node_ids[i, j + 1],
                            left_node_ids[i + 1, j + 1],
                        ]
                    )
                    conn.extend(
                        [
                            right_node_ids[i, j],
                            right_node_ids[i + 1, j],
                            right_node_ids[i, j + 1],
                            right_node_ids[i + 1, j + 1],
                        ]
                    )

            # Append connectivity for rbe element
            center_node_id = num_nodes
            center_node_xyz = np.array([1.5 * Lx, 0.5 * Ly, 0.0], dtype=dtype)
            dummy_node_id = num_nodes + 1
            dummy_node_xyz = np.zeros(3, dtype=dtype)
            num_nodes += 2
            # Add center node as dep rbe node
            rbe_conn = [center_node_id]
            # Add nodes on right edge of left plate as indep RBE nodes
            rbe_conn.extend(left_node_ids[-1, :])
            # Add nodes on left edge of right plate as indep RBE nodes
            rbe_conn.extend(right_node_ids[0, :])
            rbe_conn.append(dummy_node_id)

            # Add rbe to global connectivity
            num_rbe_nodes = len(rbe_conn)
            conn.extend(rbe_conn)
            num_elems += 1

            # Set element info for plates
            conn = np.array(conn, dtype=np.intc)
            ptr = np.arange(0, 4 * num_elems + 1, 4, dtype=np.intc)
            comp_ids = np.zeros(num_elems, dtype=np.intc)

            # Correct last entries for RBE
            ptr[-1] = ptr[-2] + num_rbe_nodes
            comp_ids[-1] = 1

            creator.setGlobalConnectivity(num_nodes, ptr, conn, comp_ids)

            # Set up the boundary conditions (fixed at left hand edge)
            bcnodes = np.append(left_node_ids[0, :], right_node_ids[-1, :])
            creator.setBoundaryConditions(bcnodes)

            # Set the node locations
            xyz = np.append(left_xyz.flatten(), right_xyz.flatten())
            xyz = np.append(xyz.flatten(), center_node_xyz)
            xyz = np.append(xyz.flatten(), dummy_node_xyz)
            creator.setNodes(xyz.flatten())

        # Set up rbe object
        num_rbe_nodes = comm.bcast(num_rbe_nodes, root=0)
        # Which dependent dofs are connected
        dep_dofs = np.array([1, 1, 1, 1, 1, 1], np.intc)
        # Which indep dof are connected
        indep_dofs = np.array([1, 1, 1, 1, 1, 1], np.intc)
        # Indep node weights, we'll just use uniform
        indep_weights = np.array([1.0])
        rbe = elements.RBE3(num_rbe_nodes, dep_dofs, indep_weights, indep_dofs)
        # Set the elements for each (only two) component
        element_list = [shell, rbe]
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

        # Apply distributed forces at tip of beam
        # Apply Qxx
        f_array[
            np.logical_and(local_x == 1.5 * Lx, local_y == 0.5 * Ly), :
        ] = applied_force

        # Create temporary dv vec for doing fd/cs
        dv_pert_array = dv_pert_vec.getArray()
        dv_pert_array[:] = 1.0

        # Create temporary state variable vec for doing fd/cs
        ans_pert_array = ans_pert_vec.getArray()
        # Define perturbation array that uniformly moves all nodes on right edge of left plate to the upward
        ans_pert_array = ans_pert_array.reshape(local_num_nodes, vars_per_node)
        ans_pert_array[local_x == Lx, 1] = 1.0

        # Define perturbation array that uniformly moves all nodes on right edge of plate to the right
        xpts_pert_array = xpts_pert_vec.getArray()
        xpts_pert_array = xpts_pert_array.reshape(local_num_nodes, 3)
        # Define perturbation array that uniformly moves all nodes on right edge of left plate to the right
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
                assembler, ksWeight=ksweight, direction=[1.0, 1.0, 1.0]
            ),
        ]
        return func_list, FUNC_REFS
