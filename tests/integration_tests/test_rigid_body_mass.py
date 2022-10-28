import os
from tacs import pytacs, functions
from pytacs_analysis_base_test import PyTACSTestCase

"""
3 point masses are rigidly connected to one another using a RBE2 element.

We create a static analysis to compute various mass balance information of 
the system (StructuralMass, CenterOfMass, MomentOfInertia) and their sensitivities.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/rigid_point_mass.bdf")


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 1  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "rigid_body_Ixx": 760.2857142857142,
        "rigid_body_Ixy": 10.5,
        "rigid_body_Ixz": -210.5,
        "rigid_body_Iyy": 478.14285714285717,
        "rigid_body_Iyz": 72.85714285714286,
        "rigid_body_Izz": 543.1428571428571,
        "rigid_body_cgx": 3.0,
        "rigid_body_cgy": -1.1428571428571428,
        "rigid_body_cgz": 3.857142857142857,
        "rigid_body_mass": 17.5,
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
            self.rtol = 2e-1
            self.atol = 1e-3
            self.dh = 1e-6

        # Instantiate FEA Assembler
        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        # Set up constitutive objects and elements
        fea_assembler.initialize()

        # Create case 1 static problem
        problem = fea_assembler.createStaticProblem("rigid_body")

        problem.addFunction("mass", functions.StructuralMass)
        problem.addFunction("cgx", functions.CenterOfMass, direction=[1.0, 0.0, 0.0])
        problem.addFunction("cgy", functions.CenterOfMass, direction=[0.0, 1.0, 0.0])
        problem.addFunction("cgz", functions.CenterOfMass, direction=[0.0, 0.0, 1.0])
        problem.addFunction(
            "Ixx",
            functions.MomentOfInertia,
            direction1=[1.0, 0.0, 0.0],
            direction2=[1.0, 0.0, 0.0],
            aboutCM=True,
        )
        problem.addFunction(
            "Ixy",
            functions.MomentOfInertia,
            direction1=[0.0, 1.0, 0.0],
            direction2=[1.0, 0.0, 0.0],
            aboutCM=True,
        )
        problem.addFunction(
            "Ixz",
            functions.MomentOfInertia,
            direction1=[0.0, 0.0, 1.0],
            direction2=[1.0, 0.0, 0.0],
            aboutCM=True,
        )
        problem.addFunction(
            "Iyy",
            functions.MomentOfInertia,
            direction1=[0.0, 1.0, 0.0],
            direction2=[0.0, 1.0, 0.0],
            aboutCM=True,
        )
        problem.addFunction(
            "Iyz",
            functions.MomentOfInertia,
            direction1=[0.0, 0.0, 1.0],
            direction2=[0.0, 1.0, 0.0],
            aboutCM=True,
        )
        problem.addFunction(
            "Izz",
            functions.MomentOfInertia,
            direction1=[0.0, 0.0, 1.0],
            direction2=[0.0, 0.0, 1.0],
            aboutCM=True,
        )
        return [problem], fea_assembler
