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

FUNC_REFS = {'PLOAD4_compliance': 279300158.48951936, 'PLOAD4_ks_disp': 9.927842420503762,
             'PLOAD4_ks_vmfailure': 29.307629374994303, 'PLOAD4_mass': 1737.357316694243,

             'PLOAD2_compliance': 279300158.48951936, 'PLOAD2_ks_disp': 9.927842420503762,
             'PLOAD2_ks_vmfailure': 29.307629374994303, 'PLOAD2_mass': 1737.357316694243}


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

        # Instantiate FEA Solver
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

        fea_solver = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        # Set up constitutive objects and elements
        fea_solver.createTACSAssembler(elem_call_back)

        return fea_solver

    def setup_tacs_vecs(self, fea_solver, dv_pert_vec, xpts_pert_vec):
        """
        Setup user-defined vectors for analysis and fd/cs sensitivity verification
        """
        # Create temporary dv vec for doing fd/cs
        dv_pert_vec[:] = 1.0

        # Define perturbation array that moves all nodes on shell
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
        fea_solver.addFunction('ks_vmfailure', functions.KSFailure,
                               ksWeight=ksweight)
        func_list = ['mass', 'compliance', 'ks_disp', 'ks_vmfailure']
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