import numpy as np
import os
from tacs import pytacs, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

'''
The nominal case is a 1m x 1m flat plate under three load cases: 
a 10 kN point force at center, a 100kPa pressure applied to the surface, and a 100G gravity load. The
perimeter of the plate is fixed in all 6 degrees of freedom. The plate comprises
100 CQUAD4 elements and test KSFailure, StructuralMass, and Compliance functions and sensitivities
'''

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/plate.bdf")

FUNC_REFS = {'point_load_compliance': 683.8571611640772, 'point_load_ks_vmfailure': 0.5757488025913641, 'point_load_mass': 12.5,
             'pressure_compliance': 4679.345460326432, 'pressure_ks_vmfailure': 1.2938623156872926, 'pressure_mass': 12.5,
             'gravity_compliance': 70.36280588344383, 'gravity_ks_vmfailure': 0.11707320009742483, 'gravity_mass': 12.5}

# KS function weight
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
            self.atol = 1e-4
            self.dh = 1e-6

        # Instantiate FEA Assembler
        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        def elem_call_back(dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs):
            # Material properties
            rho = 2500.0  # density kg/m^3
            E = 70e9  # Young's modulus (Pa)
            nu = 0.3  # Poisson's ratio
            ys = 464.0e6  # yield stress

            # Plate geometry
            tplate = 0.005  # 5 mm

            # Set up property model
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            # Set up constitutive model
            con = constitutive.IsoShellConstitutive(prop, t=tplate, tNum=dv_num)
            transform = None
            # Set up element
            elem = elements.Quad4Shell(transform, con)
            scale = [100.0]
            return elem, scale

        # Set up constitutive objects and elements
        fea_assembler.initialize(elem_call_back)

        return fea_assembler

    def setup_tacs_vecs(self, fea_assembler, dv_pert_vec, xpts_pert_vec):
        """
        Setup user-defined vectors for analysis and fd/cs sensitivity verification
        """
        # Create temporary dv vec for doing fd/cs
        dv_pert_vec[:] = 1.0

        # Define perturbation array that moves all nodes on plate
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
            problem.addFunction('ks_vmfailure', functions.KSFailure,
                                   ksWeight=ksweight)
            problem.addFunction('compliance', functions.Compliance)
        func_list = ['mass', 'ks_vmfailure', 'compliance']
        return func_list, FUNC_REFS

    def setup_tacs_problems(self, fea_assembler):
        """
        Setup pytacs object for problems we will be testing.
        """
        tacs_probs = []

        # Add point force to node 81 (center of plate)
        sp = fea_assembler.createStaticProblem(name='point_load')
        F = np.array([0.0, 0.0, 1e4, 0.0, 0.0, 0.0])
        sp.addLoadToNodes(81, F, nastranOrdering=True)
        tacs_probs.append(sp)

        # Add pressure to entire plate
        sp = fea_assembler.createStaticProblem(name='pressure')
        P = 100e3  # Pa
        compIDs = fea_assembler.selectCompIDs(include='PLATE')
        sp.addPressureToComponents(compIDs, P)
        tacs_probs.append(sp)

        # Add pressure to entire plate
        sp = fea_assembler.createStaticProblem(name='gravity')
        g = np.array([0.0, 0.0, -981.0], dtype=self.dtype)
        sp.addInertialLoad(g)
        tacs_probs.append(sp)

        return tacs_probs
