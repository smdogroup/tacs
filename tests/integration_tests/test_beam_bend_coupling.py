import os
from tacs import pytacs, functions
from pytacs_analysis_base_test import PyTACSTestCase

"""
6 noded beam model 1 meter long in x direction.
The cross-sectional properties of the beam are as follows:
    A = 0.1
    Iz = 0.2
    Iy = 0.3
    J = 0.4
    Iyz = 0.1
Because Iyz =/= 0.0, we expect some coupling to show up in y and z bending. 
We apply apply various tip loads test KSDisplacement, StructuralMass, MomentOfInertia, 
and Compliance functions and sensitivities.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/beam_model.bdf")

ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "z-shear_compliance": 1.3200259999099195,
        "z-shear_mass": 0.27,
        "z-shear_x_disp": 0.0,
        "z-shear_y_disp": -0.29391826170251795,
        "z-shear_z_disp": 12.141597029011256,
        "z-shear_I_xx": 0.0,
        "z-shear_I_xy": 0.0,
        "z-shear_I_xz": 0.0,
        "z-shear_I_yy": 0.8325,
        "z-shear_I_yz": 0.27,
        "z-shear_I_zz": 0.5625,
        "y-shear_compliance": 1.9800259997664378,
        "y-shear_mass": 0.27,
        "y-shear_x_disp": 0.0,
        "y-shear_y_disp": 18.327400292118664,
        "y-shear_z_disp": -0.2939182616984911,
        "y-shear_I_xx": 0.0,
        "y-shear_I_xy": 0.0,
        "y-shear_I_xz": 0.0,
        "y-shear_I_yy": 0.8325,
        "y-shear_I_yz": 0.27,
        "y-shear_I_zz": 0.5625,
        "x-axial_compliance": 10.000000000000027,
        "x-axial_mass": 0.27,
        "x-axial_x_disp": 95.543244182597,
        "x-axial_y_disp": 0.0,
        "x-axial_z_disp": 0.0,
        "x-axial_I_xx": 0.0,
        "x-axial_I_xy": 0.0,
        "x-axial_I_xz": 0.0,
        "x-axial_I_yy": 0.8325,
        "x-axial_I_yz": 0.27,
        "x-axial_I_zz": 0.5625,
        "x-torsion_compliance": 6.499999999999995,
        "x-torsion_mass": 0.27,
        "x-torsion_x_disp": 0.0,
        "x-torsion_y_disp": 0.0,
        "x-torsion_z_disp": 0.0,
        "x-torsion_I_xx": 0.0,
        "x-torsion_I_xy": 0.0,
        "x-torsion_I_xz": 0.0,
        "x-torsion_I_yy": 0.8325,
        "x-torsion_I_yz": 0.27,
        "x-torsion_I_zz": 0.5625,
    }

    def setup_tacs_problems(self, comm):
        """
        Setup mesh and pytacs object for problem we will be testing.
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

        # Instantiate FEA Assembler
        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        # Set up constitutive objects and elements
        fea_assembler.initialize()

        # Read in forces from BDF and create tacs struct problems
        tacs_probs = fea_assembler.createTACSProbsFromBDF()
        # Convert from dict to list
        tacs_probs = tacs_probs.values()
        # Set convergence to be tight for test
        for problem in tacs_probs:
            problem.setOption("L2Convergence", 1e-20)
            problem.setOption("L2ConvergenceRel", 1e-20)

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("compliance", functions.Compliance)
            problem.addFunction(
                "x_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[10.0, 0.0, 0.0],
            )
            problem.addFunction(
                "y_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[0.0, 10.0, 0.0],
            )
            problem.addFunction(
                "z_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[0.0, 0.0, 10.0],
            )
            problem.addFunction(
                "I_xx",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[1.0, 0.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_xy",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[0.0, 1.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_xz",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_yy",
                functions.MomentOfInertia,
                direction1=[0.0, 1.0, 0.0],
                direction2=[0.0, 1.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_yz",
                functions.MomentOfInertia,
                direction1=[0.0, 1.0, 0.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_zz",
                functions.MomentOfInertia,
                direction1=[0.0, 0.0, 1.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=True,
            )

        return tacs_probs, fea_assembler
