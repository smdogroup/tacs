import numpy as np
from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions
from static_analysis_base_test import StaticTestCase

"""
  The following example demonstrates the use of TACS on a temperature
  loaded plate.

  This code uses the TACSAssembler interface directly. Other creation
  code (TACSCreator/TACSMeshLoader) can also be used to generate
  TACSAssembler instances. Once a TACSAssembler instance has been
  created and initialized, it should be able to be used interchangeably
  withe

  Note: This code does not intelligently partition the mesh. You could
  use TACSCreator to perform the partitioning for you to achieve
  better results, but this code is designed to demonstrate the
  TACSAssembler interface itself.

  This test is based on the "tutorial" script under the examples directory.
"""

FUNC_REFS = np.array(
    [
        0.981072186947658,
        2700.0,
        30.508623116027273,
        56.92510895125915,
        1.377153264551044,
    ]
)

# Length of plate in x/y direction
Lx = 1.0
Ly = 1.0

# Number of elements in x/y direction
nx = 10
ny = 10

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
            self.dh = 1e-8

        # Get the MPI communicator size
        rank = comm.rank
        size = comm.size

        # We know in advance that the number of unknowns per node is
        # going to be equal to 3 (You can find this value by checking
        # with element->getVarsPerNode() which returns the number
        # of unknowns per node)
        vars_per_node = 3

        nodes_per_proc = ((nx + 1) * (ny + 1)) // size
        elems_per_proc = (nx * ny) // size

        num_owned_nodes = nodes_per_proc
        num_elements = elems_per_proc

        # On the last rank, adjust the ownership so we get the
        # total that we need
        if rank == size - 1:
            num_owned_nodes = (nx + 1) * (ny + 1) - nodes_per_proc * (size - 1)
            num_elements = nx * ny - elems_per_proc * (size - 1)

        # There are no dependent nodes in this problem
        num_dependent_nodes = 0
        assembler = TACS.Assembler.create(
            comm, vars_per_node, num_owned_nodes, num_elements, num_dependent_nodes
        )

        # Set the global element index for the first and last element
        # in the partition
        first_elem = rank * elems_per_proc
        first_node = rank * nodes_per_proc

        last_elem = (rank + 1) * elems_per_proc
        last_node = (rank + 1) * nodes_per_proc
        if rank == size - 1:
            last_elem = nx * ny
            last_node = (nx + 1) * (ny + 1)

        """
        The element connectivity defines the mapping between the element
        and its corresponding nodes. The node numbers are global. Since
        the number of nodes per element may vary, we also provide a
        pointer into the element connectivity array denoting the begining
        location of each element node list. This data is passed in to
        TACSAssembler directly.

        In this case we know that we only ever have 4 nodes per element.
        """

        # The elements are ordered as (i + j*nx)
        ptr = np.zeros(num_elements + 1, dtype=np.intc)
        conn = np.zeros(4 * num_elements, dtype=np.intc)

        ptr[0] = 0
        for k, elem in zip(range(num_elements), range(first_elem, last_elem)):
            # Back out the i, j coordinates from the corresponding
            # element number
            i = elem % nx
            j = elem // nx

            # Set the node connectivity
            conn[4 * k] = i + j * (nx + 1)
            conn[4 * k + 1] = i + 1 + j * (nx + 1)
            conn[4 * k + 2] = i + (j + 1) * (nx + 1)
            conn[4 * k + 3] = i + 1 + (j + 1) * (nx + 1)
            ptr[k + 1] = 4 * (k + 1)

        # Set the connectivity
        assembler.setElementConnectivity(ptr, conn)

        # Create the isotropic material class
        rho = 2700.0
        specific_heat = 921.096
        E = 70e3
        nu = 0.3
        ys = 270.0
        cte = 24.0e-6
        kappa = 230.0
        props = constitutive.MaterialProperties(
            rho=rho,
            specific_heat=specific_heat,
            E=E,
            nu=nu,
            ys=ys,
            alpha=cte,
            kappa=kappa,
        )

        """
        Create the element. This element class consists of a constitutive
        object, which stores information about the material, a model class
        which computes the variational form of the governing equations based
        on input from the basis class, which contains the basis functions and
        quadrature scheme for the element.
        """
        linear_basis = elements.LinearQuadBasis()

        element_list = []

        for elem in range(first_elem, last_elem):
            # Set the thickness design variable = the element number
            t = 1.0
            t_num = elem

            # Create the plane stress constitutive object
            stiff = constitutive.PlaneStressConstitutive(props, t=t, tNum=t_num)

            # Create the element
            model = elements.LinearThermoelasticity2D(stiff)
            elem = elements.Element2D(model, linear_basis)
            element_list.append(elem)

        # Set the elements into the mesh
        assembler.setElements(element_list)

        # Set the boundary conditions - this will only record the
        # boundary conditions on its own nodes
        for i in range(nx + 1):
            bc_nodes = np.array(
                [i, i + (nx + 1) * ny, i * (nx + 1), (i + 1) * (nx + 1) - 1],
                dtype=np.intc,
            )
            bc_values = np.array([0.0, 0.0, 6.0 * i], dtype=dtype)
            bc_vars = np.array([0, 1, 2], dtype=np.intc)
            assembler.addBCs(bc_nodes, bc_vars, bc_values)

        # Perform initialization - cannot add any more elements / vars etc
        assembler.initialize()

        # Create the node vector
        X = assembler.createNodeVec()

        # Get the local node locations
        Xpts = X.getArray()
        for k, node in zip(
            range(0, 3 * num_owned_nodes, 3), range(first_node, last_node)
        ):
            i = node % (nx + 1)
            j = node // (nx + 1)
            Xpts[k] = i * Lx / nx
            Xpts[k + 1] = j * Ly / ny

        # Reorder the vector if required
        assembler.reorderVec(X)

        # Set the node locations
        assembler.setNodes(X)

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

        # Set all components of the vector to 1.0 and apply boundary
        # conditions
        f_array = force_vec.getArray()
        f_array[:] = 6.0

        # Create dv perturbation dv for doing fd/cs
        dv_pert_array = dv_pert_vec.getArray()
        dv_pert_array[:] = 1.0

        # Define uniform random state variable perturbation array
        # Set the variable arrays
        np.random.seed(30)  # Seed random numbers for deterministic/repeatable tests
        rand_data = np.random.rand(vars_per_node * local_num_nodes).astype(self.dtype)
        ans_pert_vec.getArray()[:] = rand_data

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
            functions.AverageTemperature(assembler),
            functions.KSTemperature(assembler, ksWeight=ksweight),
            functions.KSDisplacement(
                assembler, ksWeight=ksweight, direction=[1e3, 1e3]
            ),
        ]
        func_list[0].setKSFailureType("continuous")
        return func_list, FUNC_REFS
