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

FUNC_REFS = {'analysis.avg_temp': 0.29063076564606977,
             'analysis.ks_temp': 0.6819632852575701,
             'analysis.mass': 25000.0}

wrt = ['mesh.fea_mesh.x_struct0', 'dv_struct', 'q_conduct']

# Radius of plate
R = 1.0
# Area of plate
area = np.pi * R ** 2

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
            self.rtol = 1e-6
            self.dh = 1e-50
        else:
            self.rtol = 1e-2
            self.dh = 1e-6

        # Callback function used to setup TACS element objects and DVs
        def element_callback(dv_num, comp_id, comp_descript, elem_descripts, special_dvs, **kwargs):
            # Material properties
            rho = 2500.0  # density kg/m^3
            kappa = 230.0e3  # Thermal conductivity W/(m⋅K)

            # Plate geometry
            tplate = 0.005  # 5 mm

            # Setup property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, kappa=kappa)
            # Set one thickness dv for every component
            con = constitutive.PlaneStressConstitutive(prop, t=tplate, tNum=dv_num)

            # For each element type in this component,
            # pass back the appropriate tacs element object
            elem_list = []
            model = elements.HeatConduction2D(con)
            for elem_descript in elem_descripts:
                if elem_descript in ['CQUAD4', 'CQUADR']:
                    basis = elements.LinearQuadBasis()
                elif elem_descript in ['CTRIA3', 'CTRIAR']:
                    basis = elements.LinearTriangleBasis()
                else:
                    print("Uh oh, '%s' not recognized" % (elem_descript))
                elem = elements.Element2D(model, basis)
                elem_list.append(elem)

            # Add scale for thickness dv
            scale = [100.0]
            return elem_list, scale

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
            problem.addFunction('ks_temp', functions.KSTemperature,
                                ksWeight=ksweight)
            problem.addFunction('avg_temp', functions.AverageTemperature, volume=area)

        class Top(Multipoint):

            def setup(self):
                tacs_options = {'element_callback': element_callback,
                                'problem_setup': problem_setup,
                                'mesh_file': bdf_file}

                tacs_builder = tacs.mphys.TacsBuilder(tacs_options, check_partials=True, coupled=True,
                                                      conduction=True, write_solution=False)
                tacs_builder.initialize(self.comm)
                ndv_struct = tacs_builder.get_ndv()

                dvs = self.add_subsystem('dvs', om.IndepVarComp(), promotes=['*'])
                dvs.add_output('dv_struct', np.array(ndv_struct * [5.0]))

                q_size = tacs_builder.get_ndof() * tacs_builder.get_number_of_nodes()
                heat = self.add_subsystem('heat', om.IndepVarComp(), promotes=['*'])
                heat.add_output('q_conduct', val=q_size * [1e5], distributed=True)

                self.add_subsystem('mesh', tacs_builder.get_mesh_coordinate_subsystem())
                self.mphys_add_scenario('analysis', ScenarioStructural(struct_builder=tacs_builder))
                self.connect('mesh.x_struct0', 'analysis.x_struct0')
                self.connect('dv_struct', 'analysis.dv_struct')
                self.connect('q_conduct', 'analysis.q_conduct')

        prob = om.Problem()
        prob.model = Top()
        prob.setup()

        return prob

    def setup_funcs(self):
        """
        Create a list of functions to be tested and their reference values for the problem
        """
        return FUNC_REFS, wrt