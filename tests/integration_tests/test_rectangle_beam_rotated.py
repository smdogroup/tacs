import numpy as np
import os
from tacs import pytacs, elements, constitutive, functions
from pytacs_analysis_base_test import PyTACSTestCase

"""
This is the same test cases as `test_rectangle_beam_tractions.py`, but the beam has been rotated 
about the y-axis by 45 degrees, so that it lies in a slant in the xz plane. This test ensures that the beam solution 
is invariant under trivial transformation: 

6 noded beam model 1 meter long in x' direction.
The cross-section is a solid rectangle with the following properties:
    w = 0.1
    t = 0.05
We apply two load cases: a distributed gravity and distributed traction case.
We test KSDisplacement, KSFailure, StructuralMass, and Compliance functions and sensitivities.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/beam_model_skewed.bdf")

from test_rectangle_beam_tractions import ProblemTest as PT

ksweight = 10.0

# Define rotated coordinate frame axes
x_prime = np.sqrt(0.5) * np.array([1.0, 0.0, 1.0])
y_prime = np.array([0.0, 1.0, 0.0])
z_prime = np.sqrt(0.5) * np.array([-1.0, 0.0, 1.0])


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    # Set reference functions to match unrotated case
    FUNC_REFS = PT.FUNC_REFS

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        """

        # Overwrite default check values
        if self.dtype == complex:
            self.rtol = 1e-6
            self.atol = 1e-6
            self.dh = 1e-50
        else:
            self.rtol = 2e-1
            self.atol = 1e-3
            self.dh = 1e-6

        # Material properties
        rho = 2700.0  # density kg/m^3
        E = 70.0e3  # Young's modulus (Pa)
        nu = 0.3  # Poisson's ratio
        ys = 2.7e3  # yield stress

        # Beam dimensions
        t = 0.05  # m
        w = 0.1  # m

        # Callback function used to setup TACS element objects and DVs
        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
        ):
            # Setup (isotropic) property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            con = constitutive.IsoRectangleBeamConstitutive(
                prop, t=t, tNum=dv_num, w=w, wNum=dv_num + 1
            )
            refAxis = y_prime
            transform = elements.BeamRefAxisTransform(refAxis)
            # pass back the appropriate tacs element object
            elem = elements.Beam2(transform, con)
            return elem

        # Instantiate FEA Assembler
        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        # Set up constitutive objects and elements
        fea_assembler.initialize(elem_call_back)

        grav_prob = fea_assembler.createStaticProblem("gravity")
        grav_vec = -10.0 * x_prime + 3.0 * y_prime + 5.0 * z_prime
        grav_prob.addInertialLoad(grav_vec)

        trac_prob = fea_assembler.createStaticProblem("traction")
        trac_vec = 1.0 * x_prime - 2.0 * y_prime + 3.0 * z_prime
        trac_prob.addTractionToComponents([0], trac_vec)

        tacs_probs = [grav_prob, trac_prob]

        # Set convergence to be tight for test
        for problem in tacs_probs:
            problem.setOption("L2Convergence", 1e-20)
            problem.setOption("L2ConvergenceRel", 1e-20)

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("compliance", functions.Compliance)
            problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=ksweight)
            # Get displacements in rotated coordinate frame
            problem.addFunction(
                "x_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=10.0 * x_prime,
            )
            problem.addFunction(
                "y_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=10.0 * y_prime,
            )
            problem.addFunction(
                "z_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=10.0 * z_prime,
            )

        return tacs_probs, fea_assembler
