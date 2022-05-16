import numpy as np
import os
from tacs import pytacs, elements, constitutive, functions
from pytacs_analysis_base_test import PyTACSTestCase

'''
6 noded beam model 1 meter long in x direction.
The cross-section is a solid rectangle with the following properties:
    w = 0.1
    t = 0.05
We apply two load cases: a distributed gravity and distributed traction case.
We test KSDisplacement, KSFailure, StructuralMass, and Compliance functions and sensitivities.
'''

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/beam_model.bdf")

FUNC_REFS = {'gravity_compliance': 1941.637970569529,
             'gravity_ks_vmfailure': 215.18186143357514,
             'gravity_mass': 13.500000000000004,
             'gravity_x_disp': -0.3750404568908798,
             'gravity_y_disp': 656.1378755065329,
             'gravity_z_disp': 275.45079293281793,

             'traction_compliance': 4.325878285719202,
             'traction_ks_vmfailure': 9.71446798485217,
             'traction_mass': 13.500000000000004,
             'traction_x_disp': 0.009518432861763253,
             'traction_y_disp': -0.7122094288145717,
             'traction_z_disp': 12.02223266609341}

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

        # Beam dimensions
        t = 0.05  # m
        w = 0.1  # m

        # Callback function used to setup TACS element objects and DVs
        def elem_call_back(dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs):
            # Setup (isotropic) property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            con = constitutive.IsoRectangleBeamConstitutive(prop, t=t, tNum=dv_num, w=w, wNum=dv_num+1)
            refAxis = np.array([0.0, 1.0, 0.0])
            transform = elements.BeamRefAxisTransform(refAxis)
            # pass back the appropriate tacs element object
            elem = elements.LinearBeam(transform, con)
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
        grav_prob = fea_assembler.createStaticProblem('gravity')
        grav_prob.addInertialLoad([-10.0, 3.0, 5.0])

        trac_prob = fea_assembler.createStaticProblem('traction')
        trac_prob.addTractionToComponents([0], [1.0, -2.0, 3.0])

        tacs_probs = [grav_prob, trac_prob]

        # Set convergence to be tight for test
        for problem in tacs_probs:
            problem.setOption('L2Convergence', 1e-20)
            problem.setOption('L2ConvergenceRel', 1e-20)

        return tacs_probs