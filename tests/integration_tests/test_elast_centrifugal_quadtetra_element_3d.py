import numpy as np
import os
from tacs import pytacs, TACS, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

"""
10m x 1m x 1m Square beam constructed from tetrahedral elements. The beam is cantilevered at
one end and rotated about its base at a constant angular velocity. This leads to a pure axial loading on the model.
A second case is run where the beam is hung under a uniform gravity load.

test StructuralMass and Compliance functions and sensitivities
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/solid_beam.bdf")

FUNC_REFS = {
    "Centrifugal_compliance": 21447856.343444135,
    "Centrifugal_ks_disp": 5.04476307317472,
    "Centrifugal_mass": 27000.00000000002,
}

omega = 2 * np.pi * np.array([0.0, -10.0, 0.0])
rotCenter = np.array([0.5, 0.5, 0.0])
ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    def setup_pytacs(self, comm, dtype):
        """
        Setup mesh and pytacs object for problem we will be testing.
        """

        # Overwrite default check values
        if dtype == complex:
            self.rtol = 1e-8
            self.atol = 1e-8
            self.dh = 1e-50
        else:
            self.rtol = 5e-6
            self.atol = 1e-2
            self.dh = 1e-6

        # Instantiate FEA Assembler
        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, special_dvs, **kwargs
        ):
            # Material properties
            rho = 2700.0  # density kg/m^3
            E = 70e9  # Youngs modulus Pa
            nu = 0.3  # poissons ratio
            ys = 270e6  # yield strength Pa

            # Setup property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            # Set one thickness dv for every component
            con = constitutive.SolidConstitutive(prop, t=1.0, tNum=dv_num)

            # For each element type in this component,
            # pass back the appropriate tacs element object
            elem_list = []
            model = elements.LinearElasticity3D(con)
            basis = elements.QuadraticTetrahedralBasis()
            elem = elements.Element3D(model, basis)
            return elem

        # Set up constitutive objects and elements
        fea_assembler.initialize(elem_call_back)

        return fea_assembler

    def setup_tacs_vecs(self, fea_assembler, dv_pert_vec, xpts_pert_vec):
        """
        Setup user-defined vectors for analysis and fd/cs sensitivity verification
        """
        # Create temporary dv vec for doing fd/cs
        dv_pert_vec[:] = 1.0

        # Define perturbation array that moves all nodes on plate
        xpts = fea_assembler.getOrigNodes()
        xpts_pert_vec[:] = xpts

        return

    def setup_funcs(self, fea_assembler, problems):
        """
        Create a list of functions to be tested and their reference values for the problem
        """
        # Add Functions
        for problem in problems:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("compliance", functions.Compliance)
            problem.addFunction(
                "ks_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[0.0, 0.0, 100.0],
            )
        func_list = ["mass", "compliance", "ks_disp"]
        return func_list, FUNC_REFS

    def setup_tacs_problems(self, fea_assembler):
        """
        Setup pytacs object for problems we will be testing.
        """
        # Add centrifugal load case
        static_prob = fea_assembler.createStaticProblem("Centrifugal")
        static_prob.addCentrifugalLoad(omega, rotCenter)
        return [static_prob]
