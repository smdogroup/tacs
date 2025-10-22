import os

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, functions

"""
This model features two adjacent linear hex elements, clamped on one face.
Six load conditions are considered: a pressure applied to each of the six faces using PLOAD4 cards.
This test verifies that pyTACS/TACS applies pressure load on the consistent face based on the PLOAD4 information.
 
tests KSDisplacement, StructuralMass, and Compliance functions and sensitivities.

We also add a enclosed volume constraint and test it's values/sensitivties.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/two_hexs.bdf")

ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "Clamp_+X_-x_disp": -0.2716543418951961,
        "Clamp_+X_-y_disp": 0.1944965504786221,
        "Clamp_+X_-z_disp": 0.1122082988741163,
        "Clamp_+X_compliance": 25788.7451098405,
        "Clamp_+X_mass": 5.400000000000001,
        "Clamp_+X_x_disp": 1.0274069094117648,
        "Clamp_+X_y_disp": 0.19449655047862183,
        "Clamp_+X_z_disp": 0.11220829887411632,
        "Clamp_+Y_-x_disp": 0.372199546911399,
        "Clamp_+Y_-y_disp": -0.2108822590307131,
        "Clamp_+Y_-z_disp": 0.07859767285648897,
        "Clamp_+Y_compliance": 11409.61367995888,
        "Clamp_+Y_mass": 5.400000000000001,
        "Clamp_+Y_x_disp": 0.10005467583046843,
        "Clamp_+Y_y_disp": 1.4357366919830612,
        "Clamp_+Y_z_disp": 0.07859767285648912,
        "Clamp_-X_-x_disp": 0.06931471805599453,
        "Clamp_-X_-y_disp": 0.06931471805599453,
        "Clamp_-X_-z_disp": 0.06931471805599453,
        "Clamp_-X_compliance": 0.0,
        "Clamp_-X_mass": 5.400000000000001,
        "Clamp_-X_x_disp": 0.06931471805599453,
        "Clamp_-X_y_disp": 0.06931471805599453,
        "Clamp_-X_z_disp": 0.06931471805599453,
        "Clamp_-Y_-x_disp": 0.37219954691139895,
        "Clamp_-Y_-y_disp": 1.435736691983061,
        "Clamp_-Y_-z_disp": 0.0785976728564891,
        "Clamp_-Y_compliance": 11409.61367995888,
        "Clamp_-Y_mass": 5.400000000000001,
        "Clamp_-Y_x_disp": 0.10005467583046848,
        "Clamp_-Y_y_disp": -0.21088225903071312,
        "Clamp_-Y_z_disp": 0.07859767285648898,
        "Clamp_+Z_-x_disp": 0.6483960432371099,
        "Clamp_+Z_-y_disp": 0.1708538959664814,
        "Clamp_+Z_-z_disp": -0.6532709558080134,
        "Clamp_+Z_compliance": 36425.924866297384,
        "Clamp_+Z_mass": 5.400000000000001,
        "Clamp_+Z_x_disp": 0.423860759605031,
        "Clamp_+Z_y_disp": 0.1708538959664828,
        "Clamp_+Z_z_disp": 2.6861025015727664,
        "Clamp_-Z_-x_disp": 0.6483960432371095,
        "Clamp_-Z_-y_disp": 0.17085389596648218,
        "Clamp_-Z_-z_disp": 2.6861025015727664,
        "Clamp_-Z_compliance": 36425.924866297384,
        "Clamp_-Z_mass": 5.400000000000001,
        "Clamp_-Z_x_disp": 0.4238607596050309,
        "Clamp_-Z_y_disp": 0.17085389596648243,
        "Clamp_-Z_z_disp": -0.6532709558080131,
        "volume_con": 2.0,
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
        tacs_probs = list(tacs_probs.values())
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
                "-x_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[-10.0, 0.0, 0.0],
            )
            problem.addFunction(
                "-y_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[0.0, -10.0, 0.0],
            )
            problem.addFunction(
                "-z_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[0.0, 0.0, -10.0],
            )

        constr = fea_assembler.createVolumeConstraint("volume")
        constr.addConstraint("con")
        tacs_probs.append(constr)

        return tacs_probs, fea_assembler
