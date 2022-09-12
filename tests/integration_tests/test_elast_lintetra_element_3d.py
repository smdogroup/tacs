import numpy as np
from tacs import TACS, elements, constitutive, functions
from static_analysis_base_test import StaticTestCase
import os

"""
Load in a bdf file with tetrahedral elements, apply a load,
and test KSFailure, StructuralMass, and Compliance functions and sensitivities.

This test is based on the "tetrahedral" script under the examples directory.
"""

FUNC_REFS = np.array([1.2622791020763084, 172800.0, 16.257419866831018])

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/5x5x5_cube.bdf")

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

        # Create the mesh loader object on MPI_COMM_WORLD.The
        # TACSAssembler object will be created on the same comm MPI_Comm
        mesh = TACS.MeshLoader(comm)

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
        # Create the stiffness object
        stiff = constitutive.SolidConstitutive(props, t=1.0, tNum=0)

        # Create model (need class)
        model = elements.LinearElasticity3D(stiff)
        vars_per_node = model.getVarsPerNode()
        # Set up the basis function
        linear_basis = elements.LinearTetrahedralBasis()
        quad_basis = elements.QuadraticTetrahedralBasis()

        # Create the element type (need 3D element class)
        linear_element = elements.Element3D(model, linear_basis)
        quad_element = elements.Element3D(model, quad_basis)

        # Read in bdf file
        fail = mesh.scanBDFFile(bdf_file)

        if fail is True:
            raise IOError("Failed to read in the BDF file")

        # Add the elements to the mesh loader class
        for i in range(mesh.getNumComponents()):
            elem_descript = mesh.getElementDescript(i)

            if elem_descript in ["CTETRA", "CTETRA4"]:
                elem = linear_element
            elif elem_descript == "CTETRA10":
                elem = quad_element
            else:
                elem = None

            if elem is not None:
                mesh.setElement(i, elem)

        # Now, create the TACSAssembler object
        assembler = mesh.createTACS(vars_per_node)

        return assembler

    def setup_tacs_vecs(
        self, assembler, force_vec, dv_pert_vec, ans_pert_vec, xpts_pert_vec
    ):
        """
        Setup user-defined vectors for analysis and fd/cs sensitivity verification
        """
        local_num_nodes = assembler.getNumOwnedNodes()
        vars_per_node = assembler.getVarsPerNode()

        # Create force vector
        f_array = force_vec.getArray()

        # Set uniform force on all nodes
        f_array[:] = -10.0

        # Create temporary dv vec for doing fd/cs
        dv_pert_vec.getArray()[:] = 1.0

        # Create temporary state variable vec for doing fd/cs
        np.random.seed(30)  # Seed random numbers for deterministic/repeatable tests
        rand_data = np.random.rand(vars_per_node * local_num_nodes).astype(self.dtype)
        ans_pert_vec.getArray()[:] = rand_data

        # Define perturbation array that uniformly moves all nodes on right edge of plate to the right
        np.random.seed(30)  # Seed random numbers for deterministic/repeatable tests
        rand_data = np.random.rand(3 * local_num_nodes).astype(self.dtype)
        xpts_pert_vec.getArray()[:] = rand_data

        return

    def setup_funcs(self, assembler):
        """
        Create a list of functions to be tested and their reference values for the problem
        """
        func_list = [
            functions.KSFailure(assembler, ksWeight=ksweight),
            functions.StructuralMass(assembler),
            functions.Compliance(assembler),
        ]
        return func_list, FUNC_REFS
