import numpy as np
import os
from tacs import pytacs, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

""""
The nominal case is a 1m x 1m flat plate under three load cases: 
a 1 MN point distributed force, a 10MPa pressure, and a 1 MPa traction. The
perimeter of the plate is fixed in all 6 degrees of freedom. The plate comprises
100 CQUAD4 elements and test KSFailure, KSDisplacement, StructuralMass, and Compliance functions and sensitivities

The plate is broken up into 4 components in the bdf file with the names 'PLATE.00', 'PLATE.01', etc.
The plate domains ordered as below:
           |
    3      |       2
           |
-----------|------------
           |
   0       |      1
           |
We use pyTACS selectCompIDs method to apply loads and eval functions specifically to each quadrant
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/partitioned_plate.bdf")

FUNC_REFS = {'load_compliance': 155663.05822283815, 'load_ks_disp': 38.52680169648463,
             'load_ks_vmfailure': 20.170707090964743, 'load_mass': 0.78125,

             'pressure_compliance': 73362.09834795444, 'pressure_ks_disp': 21.6030915592335,
             'pressure_ks_vmfailure': 16.464960170504657, 'pressure_mass': 0.78125,

             'traction_compliance': 735.6034203895599, 'traction_ks_disp': 1.821814969038552,
             'traction_ks_vmfailure': 0.7396840173568021, 'traction_mass': 0.78125}


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

        # Instantiate FEA Solver
        struct_options = {}

        fea_solver = pytacs.pyTACS(bdf_file, comm, options=struct_options)

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
        fea_solver.createTACSAssembler(elem_call_back)

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
        # Evaluate mass of bottom left quadrant of plate
        compIDs = fea_solver.selectCompIDs(include='PLATE.00')
        fea_solver.addFunction('mass', functions.StructuralMass, compIDs=compIDs)
        # Evaluate failure of bottom right quadrant of plate
        compIDs = fea_solver.selectCompIDs(include='PLATE.01')
        fea_solver.addFunction('ks_vmfailure', functions.KSFailure,
                               ksWeight=ksweight, compIDs=compIDs)
        # Evaluate displacement of upper left quadrant of plate
        compIDs = fea_solver.selectCompIDs(include='PLATE.03')
        fea_solver.addFunction('ks_disp', functions.KSDisplacement,
                               ksWeight=ksweight, direction= [100.0, 100.0, 100.0], compIDs=compIDs)
        # Evaluate compliance of entire plate
        fea_solver.addFunction('compliance', functions.Compliance)
        func_list = ['mass', 'ks_vmfailure', 'ks_disp', 'compliance']
        return func_list, FUNC_REFS

    def setup_tacs_problems(self, fea_solver):
        """
        Setup pytacs object for problems we will be testing.
        """
        tacs_probs = []

        # Distribute point force over all nodes in bottom left quadrant of plate
        sp = problems.StaticProblem(name='load')
        F = np.array([0.0, 0.0, 1e6, 0.0, 0.0, 0.0])
        compIDs = fea_solver.selectCompIDs(include='PLATE.00')
        fea_solver.addLoadToComponents(sp, compIDs, F)
        tacs_probs.append(sp)

        # Add pressure to bottom right quadrant of plate
        sp = problems.StaticProblem(name='pressure')
        P = 100e5  # Pa
        compIDs = fea_solver.selectCompIDs(include='PLATE.01')
        fea_solver.addPressureToComponents(sp, compIDs, P)
        tacs_probs.append(sp)

        # Add traction to upper right quadrant of plate
        sp = problems.StaticProblem(name='traction')
        trac = [1e6, 1e6, 1e6]  # N/m^2
        compIDs = fea_solver.selectCompIDs(include='PLATE.02')
        fea_solver.addTractionToComponents(sp, compIDs, trac)
        tacs_probs.append(sp)

        return tacs_probs
