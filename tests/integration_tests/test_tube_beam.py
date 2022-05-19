import numpy as np
import os
from tacs import pytacs, elements, constitutive, functions
from pytacs_analysis_base_test import PyTACSTestCase

'''
6 noded beam model 1 meter long in x direction.
The cross-section is a hollow tube with the following properties:
    d = 0.1
    t = 0.01
We apply apply various tip loads test KSDisplacement, KSFailure, 
StructuralMass, and Compliance functions and sensitivities.
'''

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/beam_model.bdf")

FUNC_REFS = {'z-shear_compliance': 2.1117936117524927,
             'z-shear_ks_vmfailure': 6.631852693984316,
             'z-shear_mass': 4.4532075864635265,
             'z-shear_x_disp': 0.0, 'z-shear_y_disp': 0.0, 'z-shear_z_disp': 19.57107469646702,

             'y-shear_compliance': 2.1117936117524927,
             'y-shear_ks_vmfailure': 6.631852693984316,
             'y-shear_mass': 4.4532075864635265,
             'y-shear_x_disp': 0.0, 'y-shear_y_disp': 19.57107469646702, 'y-shear_z_disp': 0.0,

             'x-axial_compliance': 0.008661493501599744,
             'x-axial_ks_vmfailure': 0.17273633763874147,
             'x-axial_mass': 4.4532075864635265,
             'x-axial_x_disp': 0.04641402830875186, 'x-axial_y_disp': 0.0, 'x-axial_z_disp': 0.0,

             'x-torsion_compliance': 8.151993883858584,
             'x-torsion_ks_vmfailure': 8.447664369986043,
             'x-torsion_mass': 4.4532075864635265,
             'x-torsion_x_disp': 0.0, 'x-torsion_y_disp': 0.0, 'x-torsion_z_disp': 0.0}

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

        # Material properties
        rho = 2700.0  # density kg/m^3
        E = 70.0e3  # Young's modulus (Pa)
        nu = 0.3  # Poisson's ratio
        ys = 2.7e3  # yield stress

        # Shell thickness
        t = 0.01  # m
        d = 0.1  # m

        # Callback function used to setup TACS element objects and DVs
        def elem_call_back(dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs):
            # Setup (isotropic) property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            con = constitutive.IsoTubeBeamConstitutive(prop, t=t, tNum=dv_num, d=d, dNum=dv_num+1)
            refAxis = np.array([0.0, 1.0, 0.0])
            transform = elements.BeamRefAxisTransform(refAxis)
            # pass back the appropriate tacs element object
            elem = elements.Beam2(transform, con)
            return elem

        # Instantiate FEA Assembler
        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

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
            problem.addFunction('mass', functions.StructuralMass)
            problem.addFunction('compliance', functions.Compliance)
            problem.addFunction('ks_vmfailure', functions.KSFailure,
                                ksWeight=ksweight)
            problem.addFunction('x_disp', functions.KSDisplacement,
                                ksWeight=ksweight, direction=[10.0, 0.0, 0.0])
            problem.addFunction('y_disp', functions.KSDisplacement,
                                ksWeight=ksweight, direction=[0.0, 10.0, 0.0])
            problem.addFunction('z_disp', functions.KSDisplacement,
                                ksWeight=ksweight, direction=[0.0, 0.0, 10.0])
        func_list = ['mass', 'compliance', 'x_disp', 'y_disp', 'z_disp']
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