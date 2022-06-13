import os
from tacs import pytacs, functions
from pytacs_analysis_base_test import PyTACSTestCase

'''
6 noded beam model 1 meter long in x direction with a point mass attached at the end.
The beam-mass assembly is rotated about the origin which generates a centrifugal load.
The cross-sectional properties of the beam are as follows:
    A = 0.1
    Iz = 0.2
    Iy = 0.3
    J = 0.4
    Iyz = 0.1
    M_tip = 20.0
    omega = 1.0 rev/s
Because Iyz =/= 0.0, we expect some coupling to show up in y and z bending. 
We apply apply various tip loads test KSDisplacement, StructuralMass, MomentOfInertia, 
and Compliance functions and sensitivities.
'''

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/beam_tip_mass.bdf")

FUNC_REFS = {'centrifugal_compliance': 972.8974314464419,
             'centrifugal_mass': 22.7,
             'centrifugal_x_disp': 11.787207585847018,
             'centrifugal_y_disp': 0.06931471805599453,
             'centrifugal_z_disp': 0.06931471805599453,
             'centrifugal_I_xx': 23.5, 'centrifugal_I_xy': 0.0,                'centrifugal_I_xz': 0.0,
                                       'centrifugal_I_yy': 41.519713656387665, 'centrifugal_I_yz': 2.7,
                                                                               'centrifugal_I_zz': 19.019713656387665}

ksweight = 10.0

class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.
    def setup_pytacs(self, comm, dtype):
        """
        Setup mesh and pytacs object for problem we will be testing.
        """

        # Overwrite default check values
        if dtype == complex:
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

        return fea_assembler

    def setup_tacs_vecs(self, fea_assembler, dv_pert_vec, xpts_pert_vec):
        """
        Setup user-defined vectors for analysis and fd/cs sensitivity verification
        """
        # Create temporary dv vec for doing fd/cs
        dv_pert_vec[:] = 1.0

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
            problem.addFunction('compliance', functions.Compliance)
            problem.addFunction('x_disp', functions.KSDisplacement,
                                ksWeight=ksweight, direction=[10.0, 0.0, 0.0])
            problem.addFunction('y_disp', functions.KSDisplacement,
                                ksWeight=ksweight, direction=[0.0, 10.0, 0.0])
            problem.addFunction('z_disp', functions.KSDisplacement,
                                ksWeight=ksweight, direction=[0.0, 0.0, 10.0])
            problem.addFunction('I_xx', functions.MomentOfInertia, direction1=[1.0, 0.0, 0.0],
                                direction2=[1.0, 0.0, 0.0], aboutCM=True)
            problem.addFunction('I_xy', functions.MomentOfInertia, direction1=[1.0, 0.0, 0.0],
                                direction2=[0.0, 1.0, 0.0], aboutCM=True)
            problem.addFunction('I_xz', functions.MomentOfInertia, direction1=[1.0, 0.0, 0.0],
                                direction2=[0.0, 0.0, 1.0], aboutCM=True)
            problem.addFunction('I_yy', functions.MomentOfInertia, direction1=[0.0, 1.0, 0.0],
                                direction2=[0.0, 1.0, 0.0], aboutCM=True)
            problem.addFunction('I_yz', functions.MomentOfInertia, direction1=[0.0, 1.0, 0.0],
                                direction2=[0.0, 0.0, 1.0], aboutCM=True)
            problem.addFunction('I_zz', functions.MomentOfInertia, direction1=[0.0, 0.0, 1.0],
                                direction2=[0.0, 0.0, 1.0], aboutCM=True)
        func_list = ['mass', 'compliance', 'x_disp', 'y_disp', 'z_disp', 'I_xx', 'I_xy', 'I_xz', 'I_yy', 'I_yz', 'I_zz']
        return func_list, FUNC_REFS

    def setup_tacs_problems(self, fea_assembler):
        """
        Setup pytacs object for problems we will be testing.
        """
        # Read in forces from BDF and create tacs struct problems
        tacs_probs = fea_assembler.createTACSProbsFromBDF()
        # Convert from dict to list
        tacs_probs = tacs_probs.values()
        # Set convergence to be tight for test
        for problem in tacs_probs:
            problem.setOption('L2Convergence', 1e-20)
            problem.setOption('L2ConvergenceRel', 1e-20)

        return tacs_probs