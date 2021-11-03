import numpy as np
import os
from tacs import pytacs, TACS, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

'''
 10m x 1m x 1m Square beam constructed from tetrahedral elements. The beam is cantilevered at
one end and loaded using a pressure on the opposite face. This leads to a pure axial loading on the model.

test StructuralMass and Compliance functions and sensitivities
'''

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/solid_beam.bdf")

FUNC_REFS = {'AxialPressure_compliance': 14226.813534942703,
             'AxialPressure_ks_disp': 0.3088436763967383,
             'AxialPressure_mass': 27000.00000000002}

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
            self.rtol = 1e-6
            self.atol = 1e-2
            self.dh = 1e-6

        # Instantiate FEA Solver
        struct_options = {}

        fea_solver = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        # Set up constitutive objects and elements
        fea_solver.createTACSAssembler()

        return fea_solver

    def setup_tacs_vecs(self, fea_solver, dv_pert_vec, xpts_pert_vec):
        """
        Setup user-defined vectors for analysis and fd/cs sensitivity verification
        """
        # Create temporary dv vec for doing fd/cs
        dv_pert_vec[:] = 1.0

        # Define perturbation array that moves all nodes on plate
        xpts = fea_solver.getCoordinates()
        xpts_pert_vec[:] = xpts

        return

    def setup_funcs(self, fea_solver):
        """
        Create a list of functions to be tested and their reference values for the problem
        """
        # Add Functions
        fea_solver.addFunction('mass', functions.StructuralMass)
        fea_solver.addFunction('compliance', functions.Compliance)
        fea_solver.addFunction('ks_disp', functions.KSDisplacement,
                               ksWeight=ksweight, direction=[-100.0, -100.0, -100.0])
        func_list = ['mass', 'compliance', 'ks_disp']
        return func_list, FUNC_REFS

    def setup_tacs_problems(self, fea_solver):
        """
        Setup pytacs object for problems we will be testing.
        """
        # Read in forces from BDF and create tacs struct problems
        tacs_probs = fea_solver.createTACSProbsFromBDF()
        # Convert from dict to list
        tacs_probs = tacs_probs.values()
        return tacs_probs