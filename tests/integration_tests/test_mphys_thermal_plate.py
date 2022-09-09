import numpy as np
import os
from tacs import elements, constitutive, functions
import tacs.mphys
from openmdao_analysis_base_test import OpenMDAOTestCase
from mphys.multipoint import Multipoint
from mphys.scenario_structural import ScenarioStructural
import openmdao.api as om

"""
This is a simple 1m by 2m plate made up of four quad thermal elements.
The plate is thermally loaded under a unit heat flow, 
"q_conduct", applied on on every node. The mass, AverageTemperature, and KSTemperature 
of the plate are evaluated as outputs and have their partial and total sensitivities checked.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/debug_plate.bdf")

# Historical reference values for function outputs
FUNC_REFS = {
    "analysis.avg_temp": 0.4565217391304351,
    "analysis.ks_temp": 0.6819632852575701,
    "analysis.mass": 25000.0,
}

# Inputs to check total sensitivities wrt
wrt = ["mesh.fea_mesh.x_struct0", "dv_struct", "q_conduct"]

# Area of plate
area = 2.0

# KS function weight
ksweight = 10.0


class ProblemTest(OpenMDAOTestCase.OpenMDAOTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    def setup_problem(self, dtype):
        """
        Setup openmdao problem object we will be testing.
        """

        # Overwrite default tolerances
        if dtype == complex:
            self.rtol = 1e-6
            self.dh = 1e-50
        else:
            self.rtol = 1e-2
            self.dh = 1e-6

        # Callback function used to setup TACS element objects and DVs
        def element_callback(
            dv_num, comp_id, comp_descript, elem_descripts, special_dvs, **kwargs
        ):
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
            basis = elements.LinearQuadBasis()
            elem = elements.Element2D(model, basis)

            # Add scale for thickness dv
            scale = [100.0]
            return elem, scale

        def problem_setup(scenario_name, fea_assembler, problem):
            """
            Helper function to add fixed forces and eval functions
            to thermal problems used in tacs builder
            """
            # Set convergence to be tight for test
            problem.setOption("L2Convergence", 1e-20)
            problem.setOption("L2ConvergenceRel", 1e-20)

            # Add TACS Functions
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("ks_temp", functions.KSTemperature, ksWeight=ksweight)
            problem.addFunction("avg_temp", functions.AverageTemperature, volume=area)

        class Top(Multipoint):
            def setup(self):
                tacs_options = {
                    "element_callback": element_callback,
                    "problem_setup": problem_setup,
                    "mesh_file": bdf_file,
                }

                tacs_builder = tacs.mphys.TacsBuilder(
                    tacs_options,
                    check_partials=True,
                    coupled=True,
                    conduction=True,
                    write_solution=False,
                )
                tacs_builder.initialize(self.comm)
                ndv_struct = tacs_builder.get_ndv()

                dvs = self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])
                dvs.add_output("dv_struct", np.array(ndv_struct * [5.0]))

                q_size = tacs_builder.get_ndof() * tacs_builder.get_number_of_nodes()
                heat = self.add_subsystem("heat", om.IndepVarComp(), promotes=["*"])
                heat.add_output("q_conduct", val=q_size * [1e5], distributed=True)

                self.add_subsystem("mesh", tacs_builder.get_mesh_coordinate_subsystem())
                self.mphys_add_scenario(
                    "analysis", ScenarioStructural(struct_builder=tacs_builder)
                )
                self.connect("mesh.x_struct0", "analysis.x_struct0")
                self.connect("dv_struct", "analysis.dv_struct")
                self.connect("q_conduct", "analysis.q_conduct")

        prob = om.Problem()
        prob.model = Top()

        return prob

    def setup_funcs(self):
        """
        Create a dict of functions to be tested and a list of inputs
        to test their sensitivities with respect to.
        """
        return FUNC_REFS, wrt
