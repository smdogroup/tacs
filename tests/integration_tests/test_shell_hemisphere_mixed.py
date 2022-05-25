import numpy as np
import os
from tacs import pytacs, TACS, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

'''
Hemispherical shell constructed from mized quad/tri shell elements. 
The shell is subjected to an inward pressure and is supported at the rim.
The loads are applied in two euivilent load cases through the bdf:
    1. Using a PLOAD2 card
    2. Using a PLOAD4 card

tests StructuralMass, KSFailure, KSDisplacement and Compliance functions and sensitivities
'''

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/hemisphere.bdf")

FUNC_REFS = {'PLOAD4_cg_x': 0.0009653731820888509, 'PLOAD4_cg_y': -9.14227063766091e-05, 'PLOAD4_cg_z': 0.49758219135768283,
             'PLOAD4_I_xx': 721.1880296374977, 'PLOAD4_I_xy': 0.02506752378214014, 'PLOAD4_I_xz': -0.27765576618522025,
                                               'PLOAD4_I_yy': 718.8517782438879,   'PLOAD4_I_yz': -0.033253217390545584,
                                                                                   'PLOAD4_I_zz': 1152.5803780307508,
             'PLOAD4_compliance': 279300158.48951936,
             'PLOAD4_ks_disp': 9.927842420503762,
             'PLOAD4_ks_vmfailure': 29.307629374994303,
             'PLOAD4_mass': 1737.357316694243,

             'PLOAD2_cg_x': 0.0009653731820888509, 'PLOAD2_cg_y': -9.14227063766091e-05, 'PLOAD2_cg_z': 0.49758219135768283,
             'PLOAD2_I_xx': 721.1880296374977, 'PLOAD2_I_xy': 0.02506752378214014, 'PLOAD2_I_xz': -0.27765576618522025,
                                               'PLOAD2_I_yy': 718.8517782438879,   'PLOAD2_I_yz': -0.033253217390545584,
                                                                                   'PLOAD2_I_zz': 1152.5803780307508,
             'PLOAD2_compliance': 279300158.48951936,
             'PLOAD2_ks_disp': 9.927842420503762,
             'PLOAD2_ks_vmfailure': 29.307629374994303,
             'PLOAD2_mass': 1737.357316694243}

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
            self.rtol = 2e-1
            self.atol = 1e-3
            self.dh = 1e-6

        # Instantiate FEA Assembler
        struct_options = {}

        def elem_call_back(dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs):
            # Material properties
            rho = 2780.0  # density kg/m^3
            E = 70e9  # Young's modulus (Pa)
            nu = 0.3  # Poisson's ratio
            ys = 270e6  # yield stress

            # Shell thickness
            t = 0.1  # m

            # Setup (isotropic) property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            # Set one thickness dv for every component
            con = constitutive.IsoShellConstitutive(prop, t=t, tNum=dv_num)

            # For each element type in this component,
            # pass back the appropriate tacs element object
            elem_list = []
            transform = None
            for elem_descript in elem_descripts:
                if elem_descript in ['CQUAD4', 'CQUADR']:
                    elem = elements.Quad4Shell(transform, con)
                elif elem_descript in ['CTRIA3', 'CTRIAR']:
                    elem = elements.Tri3Shell(transform, con)
                else:
                    print("Uh oh, '%s' not recognized" % (elem_descript))
                elem_list.append(elem)

            # Add scale for thickness dv
            scale = [100.0]
            return elem_list, scale

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
            problem.addFunction('ks_disp', functions.KSDisplacement,
                                ksWeight=ksweight, direction=[-100.0, -100.0, -100.0])
            problem.addFunction('ks_vmfailure', functions.KSFailure,
                                ksWeight=ksweight)
            problem.addFunction('cg_x', functions.CenterOfMass, direction=[1.0, 0.0, 0.0])
            problem.addFunction('cg_y', functions.CenterOfMass, direction=[0.0, 1.0, 0.0])
            problem.addFunction('cg_z', functions.CenterOfMass, direction=[0.0, 0.0, 1.0])
            problem.addFunction('I_xx', functions.MomentOfInertia, direction1=[1.0, 0.0, 0.0],
                                direction2=[1.0, 0.0, 0.0], aboutCG=True)
            problem.addFunction('I_xy', functions.MomentOfInertia, direction1=[1.0, 0.0, 0.0],
                                direction2=[0.0, 1.0, 0.0], aboutCG=True)
            problem.addFunction('I_xz', functions.MomentOfInertia, direction1=[1.0, 0.0, 0.0],
                                direction2=[0.0, 0.0, 1.0], aboutCG=True)
            problem.addFunction('I_yy', functions.MomentOfInertia, direction1=[0.0, 1.0, 0.0],
                                direction2=[0.0, 1.0, 0.0], aboutCG=True)
            problem.addFunction('I_yz', functions.MomentOfInertia, direction1=[0.0, 0.0, 1.0],
                                direction2=[0.0, 1.0, 0.0], aboutCG=True)
            problem.addFunction('I_zz', functions.MomentOfInertia, direction1=[0.0, 0.0, 1.0],
                                direction2=[0.0, 0.0, 1.0], aboutCG=True)
        func_list = ['mass', 'compliance', 'ks_disp', 'ks_vmfailure', 'cg_x', 'cg_y', 'cg_z',
                     'I_xx', 'I_xy', 'I_xz', 'I_yy', 'I_yz', 'I_zz']
        return func_list, FUNC_REFS

    def setup_tacs_problems(self, fea_assembler):
        """
        Setup pytacs object for problems we will be testing.
        """
        # Read in forces from BDF and create tacs struct problems
        tacs_probs = fea_assembler.createTACSProbsFromBDF()
        # Convert from dict to list
        tacs_probs = tacs_probs.values()
        return tacs_probs