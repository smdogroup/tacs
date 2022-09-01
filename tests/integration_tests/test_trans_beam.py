import os
import numpy as np
from tacs import pytacs, constitutive, elements, functions
from pytacs_analysis_base_test import PyTACSTestCase

'''
6 noded beam model 1 meter long in x direction with a two transient tip load cases:
    1. linear ramp
    2. sinusoidal
The transient loads are read in from the BDF using the createTACSProbsFromBDF method. 
We apply apply various tip loads test KSDisplacement and Compliance functions and sensitivities.
'''

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/transient_beam.bdf")

FUNC_REFS = {'ramp_compliance': 10.507315462331803,
             'ramp_x_disp': 0.06931471805599457,
             'ramp_y_disp': 12.186052999423167,
             'ramp_z_disp': 0.06931471805599457,

             'sinusoid_compliance': 22.093569888841497,
             'sinusoid_x_disp': 0.06931471805599457,
             'sinusoid_y_disp': 25.654265895487846,
             'sinusoid_z_disp': 0.06931471805599457}

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
            self.rtol = 2e-3
            self.atol = 1e-3
            self.dh = 1e-6

        # Instantiate FEA Assembler
        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        # Material properties
        rho = 27.0  # density kg/m^3
        E = 70.0e2  # Young's modulus (Pa)
        nu = 0.3  # Poisson's ratio
        ys = 2.7e-2  # yield stress

        # Beam dimensions
        t = 0.2     # m
        w = 0.5  # m

        # Callback function used to setup TACS element objects and DVs
        def elem_call_back(dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs):
            # Setup (isotropic) property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            con = constitutive.IsoRectangleBeamConstitutive(prop, t=t, tNum=dv_num, w=w, wNum=dv_num+1)
            refAxis = np.array([0.0, 1.0, 0.0])
            transform = elements.BeamRefAxisTransform(refAxis)
            # pass back the appropriate tacs element object
            elem = elements.Beam2(transform, con)
            return elem

        # Instantiate FEA Assembler
        fea_assembler = pytacs.pyTACS(bdf_file, comm)

        # Set up constitutive objects and elements
        fea_assembler.initialize(elem_call_back)

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
            problem.addFunction('compliance', functions.Compliance)
            problem.addFunction('x_disp', functions.KSDisplacement,
                                ksWeight=ksweight, direction=[10.0, 0.0, 0.0])
            problem.addFunction('y_disp', functions.KSDisplacement,
                                ksWeight=ksweight, direction=[0.0, 10.0, 0.0])
            problem.addFunction('z_disp', functions.KSDisplacement,
                                ksWeight=ksweight, direction=[0.0, 0.0, 10.0])
        func_list = ['compliance', 'x_disp', 'y_disp', 'z_disp']
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