import numpy as np
from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions
from static_analysis_base_test import StaticTestCase

"""
Create a cantilevered composite beam of linear quad shells with an 
unbalanced tip shear load on the right corner
and test KSFailure, StructuralMass, and Compliance functions and sensitivities
"""

FUNC_REFS = np.array(
    [5812.5, 63907185.558059536, 12.21799417804536, 174.71401901274177]
)

# Length of plate in x/y direction
Lx = 10.0
Ly = 1.0

# Number of elements in x/y direction
nx = 10
ny = 10

# Shear load (N)
Vx = 5e6

# KS function weight
ksweight = 10.0

# Layup parameters
nplies = 3
ply_thickness = 0.125
ply_angles = np.array([0.0, -45.0, 90.0]) * np.pi / 180.0
# Axis that defines 0 deg direction for layup
ref_axis = np.array([0.0, 1.0, 0.0])


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
            self.dh = 1e-5

        # Create the stiffness object
        # Create the orthotropic layup
        rho = 1550.0
        specific_heat = 921.096
        E1 = 54e9
        E2 = 18e9
        nu12 = 0.25
        G12 = 9e9
        G13 = 9e9
        Xt = 2410.0e6
        Xc = 1040.0e6
        Yt = 73.0e6
        Yc = 173.0e6
        S12 = 71.0e6
        cte = 24.0e-6
        kappa = 230.0
        ortho_prop = constitutive.MaterialProperties(
            rho=rho,
            specific_heat=specific_heat,
            E1=E1,
            E2=E2,
            nu12=nu12,
            G12=G12,
            G13=G13,
            G23=G13,
            Xt=Xt,
            Xc=Xc,
            Yt=Yt,
            Yc=Yc,
            S12=S12,
            cte=cte,
            kappa=kappa,
        )

        ortho_ply = constitutive.OrthotropicPly(ply_thickness, ortho_prop)
        ortho_layup = [ortho_ply] * nplies

        ply_thicknesses = np.array([ply_thickness] * nplies, dtype=self.dtype)
        stiff = constitutive.CompositeShellConstitutive(
            ortho_layup, ply_thicknesses.astype(dtype), ply_angles.astype(dtype)
        )

        # Set up the element transform function
        transform = elements.ShellRefAxisTransform(ref_axis)
        elem = elements.Quad4Shell(transform, stiff)

        # Allocate the TACSCreator object
        vars_per_node = elem.getVarsPerNode()
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
        f_array[np.logical_and(local_x == Lx, local_y == 0.0), 2] += Vx

        # There are no dvs associated with the composite class

        # Create temporary state variable vec for doing fd/cs
        ans_pert_array = ans_pert_vec.getArray()
        # Define perturbation array that uniformly moves all nodes on right edge of plate to upward
        ans_pert_array = ans_pert_array.reshape(local_num_nodes, vars_per_node)
        ans_pert_array[local_x == Lx, 2] = 1.0

        # Define perturbation array that uniformly moves all nodes on right edge of plate to the right
        xpts_pert_array = xpts_pert_vec.getArray()
        xpts_pert_array = xpts_pert_array.reshape(local_num_nodes, 3)
        xpts_pert_array[local_x == Lx, 0] = 1.0

        return

    def setup_funcs(self, assembler):
        """
        Create a list of functions to be tested and their reference values for the problem
        """
        func_list = [
            functions.StructuralMass(assembler),
            functions.Compliance(assembler),
            functions.KSDisplacement(
                assembler, ksWeight=ksweight, direction=[0.0, 0.0, 1.0]
            ),
            functions.KSFailure(assembler, ksWeight=ksweight, safetyFactor=1.5),
        ]
        return func_list, FUNC_REFS
