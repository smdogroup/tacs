import numpy as np
import os
from tacs import pytacs, TACS, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

'''
A point mass element free in space. 
Case 1: We apply a unit force in every direction on the mass for 10 seconds. 
Case 2: We apply a gravity load in the z direction for 10 seconds.

The inertial values for the mass are given below:
mass = 20.0 kg
Ixx = 23.5 kg*m^2
Iyy = 32.6 kg*m^2 
Izz = 12.8 kg*m^2 
'''

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/point_mass.bdf")

FUNC_REFS = {'constant_force_mass': 200.0,
             'constant_force_x_disp': 2.343553337720806,
             'constant_force_y_disp': 2.343553337720806,
             'constant_force_z_disp': 2.343553337720806,
             'gravity_mass': 200.0,
             'gravity_x_disp': 0.23025850929940442,
             'gravity_y_disp': 0.23025850929940442,
             'gravity_z_disp': 490.27400177252474}

# Force to apply
f = np.ones(6)

ksweight = 10.0

class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 1  # this is how many MPI processes to use for this TestCase.
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
            self.rtol = 2e-1
            self.atol = 1e-3
            self.dh = 1e-6

        # Instantiate FEA Assembler
        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        # Set up constitutive objects and elements
        fea_assembler.initialize()

        return fea_assembler

    def setup_tacs_vecs(self, fea_assembler, dv_pert_vec, xpts_pert_vec):
        """
        Setup user-defined vectors for analysis and fd/cs sensitivity verification
        """
        # This problem has no dvs

        # Define perturbation array that moves all nodes on shell
        xpts = fea_assembler.getOrigNodes()
        xpts_pert_vec[:] = xpts

        return

    def setup_funcs(self, fea_assembler, problems):
        """
        Create a list of functions to be tested and their reference values for the problem
        """
        # Add Functions
        for problem in problems:
            problem.addFunction('mass', functions.StructuralMass)
            problem.addFunction('x_disp', functions.KSDisplacement,
                                ksWeight=ksweight, direction=[1.0, 0.0, 0.0])
            problem.addFunction('y_disp', functions.KSDisplacement,
                                ksWeight=ksweight, direction=[0.0, 1.0, 0.0])
            problem.addFunction('z_disp', functions.KSDisplacement,
                                ksWeight=ksweight, direction=[0.0, 0.0, 1.0])
        func_list = ['mass', 'x_disp', 'y_disp', 'z_disp']
        return func_list, FUNC_REFS

    def setup_tacs_problems(self, fea_assembler):
        """
        Setup pytacs object for problems we will be testing.
        """
        # List to hold all problems
        allProblems = []

        # Create case 1 transient problem
        problem = fea_assembler.createTransientProblem('constant_force', 0.0, 10.0, 100)
        timeSteps = problem.getTimeSteps()
        for step_i, time in enumerate(timeSteps):
            problem.addLoadToNodes(step_i, 0, f, nastranOrdering=False)
        allProblems.append(problem)

        # Create case 2 transient problem
        problem = fea_assembler.createTransientProblem('gravity', 0.0, 10.0, 100)
        g = np.array([0.0, 0.0, 9.81], dtype=TACS.dtype)
        for step_i, time in enumerate(timeSteps):
            problem.addInertialLoad(step_i, g)
        allProblems.append(problem)

        return allProblems