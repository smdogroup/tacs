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

omega = 2 * np.pi * np.array([0.0, -10.0, 0.0])
rotCenter = np.array([0.5, 0.5, 0.0])
ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "Centrifugal_compliance": 21447856.343444135,
        "Centrifugal_ks_disp": 5.04476307317472,
        "Centrifugal_mass": 27000.00000000002,
    }

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        """

        # Overwrite default check values
        if self.dtype == complex:
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

        # Add centrifugal load case
        static_prob = fea_assembler.createStaticProblem("Centrifugal")
        static_prob.addCentrifugalLoad(omega, rotCenter)

        # Add Functions
        static_prob.addFunction("mass", functions.StructuralMass)
        static_prob.addFunction("compliance", functions.Compliance)
        static_prob.addFunction(
            "ks_disp",
            functions.KSDisplacement,
            ksWeight=ksweight,
            direction=[0.0, 0.0, 100.0],
        )

        return [static_prob], fea_assembler
