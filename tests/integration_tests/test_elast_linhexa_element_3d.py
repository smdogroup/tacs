import numpy as np
from tacs import TACS, elements, constitutive, functions
from static_analysis_base_test import StaticTestCase

"""
Create a uniform cube under distributed point loads
and test KSFailure, StructuralMass, CenterOfMass, MomentOfInertia, and Compliance functions and sensitivities
"""

FUNC_REFS = np.array(
    [
        1.1074017881777185,
        2570000.0,
        178995399.11928475,
        7.192551651904685,
        5.00000000e00,
        5.00000000e00,
        5.00000000e00,
        171333333.33333328,
        -64250000.00000006,
        -64250000.00000006,
        171333333.33333328,
        -64250000.00000006,
        171333333.33333328,
    ]
)

# Length of plate in x/y direction
Lx = 10.0
Ly = 10.0
Lz = 10.0

# Number of elements in x/y direction
nx = 3
ny = 3
nz = 3

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
            self.atol = 1e-6
            self.dh = 1e-50
        else:
            self.rtol = 1e-2
            self.atol = 1e-4
            self.dh = 1e-9

        # Create the stiffness object
        props = constitutive.MaterialProperties(rho=2570.0, E=70e9, nu=0.3, ys=350e6)
        stiff = constitutive.SolidConstitutive(props, t=1.0, tNum=0)

        # Set up the basis function
        model = elements.LinearElasticity3D(stiff)
        basis = elements.LinearHexaBasis()
        elem = elements.Element3D(model, basis)

        # Allocate the TACSCreator object
        vars_per_node = model.getVarsPerNode()
        creator = TACS.Creator(comm, vars_per_node)

        if comm.rank == 0:
            num_elems = nx * ny * nz
            num_nodes = (nx + 1) * (ny + 1) * (nz + 1)

            x = np.linspace(0, Lx, nx + 1, dtype)
            y = np.linspace(0, Ly, ny + 1, dtype)
            z = np.linspace(0, Lz, nz + 1, dtype)
            xyz = np.zeros([nx + 1, ny + 1, nz + 1, 3], dtype)
            xyz[:, :, :, 0], xyz[:, :, :, 1], xyz[:, :, :, 2] = np.meshgrid(
                x, y, z, indexing="ij"
            )
            x, y, z = xyz[:, :, :, 0], xyz[:, :, :, 1], xyz[:, :, :, 2]

            node_ids = np.arange(num_nodes).reshape(nx + 1, ny + 1, nz + 1)

            conn = []
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        conn.append(
                            [
                                node_ids[i, j, k],
                                node_ids[i + 1, j, k],
                                node_ids[i, j + 1, k],
                                node_ids[i + 1, j + 1, k],
                                node_ids[i, j, k + 1],
                                node_ids[i + 1, j, k + 1],
                                node_ids[i, j + 1, k + 1],
                                node_ids[i + 1, j + 1, k + 1],
                            ]
                        )

            conn = np.array(conn, dtype=np.intc).flatten()
            ptr = np.arange(0, 8 * num_elems + 1, 8, dtype=np.intc)
            comp_ids = np.zeros(num_elems, dtype=np.intc)

            creator.setGlobalConnectivity(num_nodes, ptr, conn, comp_ids)

            # Set up the boundary conditions ( clamp at x == 0  face)
            bcnodes = np.array(node_ids[x == 0.0], dtype=np.intc)
            creator.setBoundaryConditions(bcnodes)

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

        # Apply uniform distributed forces on all nodes
        f_array[:] = 1e8

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
            functions.KSFailure(assembler, ksWeight=ksweight),
            functions.StructuralMass(assembler),
            functions.Compliance(assembler),
            functions.KSDisplacement(
                assembler, ksWeight=ksweight, direction=[100.0, 100.0, 100.0]
            ),
            functions.CenterOfMass(assembler, direction=[1.0, 0.0, 0.0]),
            functions.CenterOfMass(assembler, direction=[0.0, 1.0, 0.0]),
            functions.CenterOfMass(assembler, direction=[0.0, 0.0, 1.0]),
            functions.MomentOfInertia(
                assembler,
                direction1=[1.0, 0.0, 0.0],
                direction2=[1.0, 0.0, 0.0],
                aboutCM=False,
            ),
            functions.MomentOfInertia(
                assembler,
                direction1=[0.0, 1.0, 0.0],
                direction2=[1.0, 0.0, 0.0],
                aboutCM=False,
            ),
            functions.MomentOfInertia(
                assembler,
                direction1=[0.0, 0.0, 1.0],
                direction2=[1.0, 0.0, 0.0],
                aboutCM=False,
            ),
            functions.MomentOfInertia(
                assembler,
                direction1=[0.0, 1.0, 0.0],
                direction2=[0.0, 1.0, 0.0],
                aboutCM=False,
            ),
            functions.MomentOfInertia(
                assembler,
                direction1=[0.0, 1.0, 0.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=False,
            ),
            functions.MomentOfInertia(
                assembler,
                direction1=[0.0, 0.0, 1.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=False,
            ),
        ]
        return func_list, FUNC_REFS
