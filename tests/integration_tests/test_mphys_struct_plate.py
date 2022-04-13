import numpy as np
import os
from tacs import elements, constitutive, functions
import tacs.mphys
from openmdao_analysis_base_test import OpenMDAOTestCase
from mphys.multipoint import Multipoint
from mphys.scenario_structural import ScenarioStructural
import openmdao.api as om

'''
The nominal case is a heat conduction problem of a
1m radius plate with a Dirichilet boundary condition applied at the edges,
such that:
    T(theta) = T0 + dT * sin(2*theta)
    T0 = 70 C
    dT = 30 C
The problem is then to solve for the temperature within the boundaries of the plate.
The problem basically boils down to Laplaces problem:
    grad**2 T = 0
test KSTemperature, StructuralMass, and AverageTemperature functions and sensitivities
'''

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/debug_plate.bdf")

FUNC_REFS = {'analysis.mass': 55.6,
             'analysis.ks_vmfailure': 5.778130269059719}

wrt = ['mesh.fea_mesh.x_struct0', 'dv_struct']

# KS function weight
ksweight = 10.0


class ProblemTest(OpenMDAOTestCase.OpenMDAOTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    def setup_problem(self, dtype):
        """
        Setup mesh and pytacs object for problem we will be testing.
        """

        # Overwrite default tolerances
        if dtype == complex:
            self.rtol = 1e-8
            self.dh = 1e-50
        else:
            self.rtol = 1e-2
            self.dh = 1e-6

        # Callback function used to setup TACS element objects and DVs
        def element_callback(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
            rho = 2780.0  # density, kg/m^3
            E = 73.1e9  # elastic modulus, Pa
            nu = 0.33  # poisson's ratio
            ys = 324.0e6  # yield stress, Pa
            thickness = 0.012
            min_thickness = 0.002
            max_thickness = 0.05

            # Setup (isotropic) property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            # Set one thickness dv for every component
            con = constitutive.IsoShellConstitutive(prop, t=thickness, tNum=dvNum, tlb=min_thickness, tub=max_thickness)

            # For each element type in this component,
            # pass back the appropriate tacs element object
            transform = None
            elem = elements.Quad4Shell(transform, con)

            return elem

        def problem_setup(scenario_name, fea_assembler, problem):
            """
            Helper function to add fixed forces and eval functions
            to structural problems used in tacs builder
            """
            # Set convergence to be tight for test
            problem.setOption('L2Convergence', 1e-20)
            problem.setOption('L2ConvergenceRel', 1e-20)

            # Add TACS Functions
            problem.addFunction('mass', functions.StructuralMass)
            problem.addFunction('ks_vmfailure', functions.KSFailure, ksWeight=ksweight)

            # Add gravity load
            g = np.array([0.0, 0.0, -9.81])*100  # m/s^2
            problem.addInertialLoad(g)

        class Top(Multipoint):

            def setup(self):
                tacs_options = {'element_callback': element_callback,
                                'problem_setup': problem_setup,
                                'mesh_file': bdf_file}

                tacs_builder = tacs.mphys.TacsBuilder(tacs_options, check_partials=True, coupled=True, write_solution=False)
                tacs_builder.initialize(self.comm)
                ndv_struct = tacs_builder.get_ndv()

                dvs = self.add_subsystem('dvs', om.IndepVarComp(), promotes=['*'])
                dvs.add_output('dv_struct', np.array(ndv_struct * [0.01]))

                f_size = tacs_builder.get_ndof() * tacs_builder.get_number_of_nodes()
                forces = self.add_subsystem('forces', om.IndepVarComp(), promotes=['*'])
                forces.add_output('f_struct', val=np.ones(f_size), distributed=True)

                self.add_subsystem('mesh', tacs_builder.get_mesh_coordinate_subsystem())
                self.mphys_add_scenario('analysis', ScenarioStructural(struct_builder=tacs_builder))
                self.connect('mesh.x_struct0', 'analysis.x_struct0')
                self.connect('dv_struct', 'analysis.dv_struct')
                self.connect('f_struct', 'analysis.f_struct')

        prob = om.Problem()
        prob.model = Top()

        return prob

    def setup_funcs(self):
        """
        Create a list of functions to be tested and their reference values for the problem
        """
        return FUNC_REFS, wrt