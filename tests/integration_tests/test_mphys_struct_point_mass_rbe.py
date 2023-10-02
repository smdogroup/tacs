import os

import numpy as np
import openmdao.api as om
from mphys.multipoint import Multipoint
from mphys.scenario_structural import ScenarioStructural

import tacs.mphys
from openmdao_analysis_base_test import OpenMDAOTestCase
from tacs import elements, constitutive, functions

"""
This is a simple 1m by 2m plate made up of four quad shell elements.
A point mass is added to the tip of the plate and is attached using a RBE2 element.
The plate is structurally loaded under a 100G gravity load and a unit force, 
"f_struct", is applied on on every node. The mass and KSFailure of the plate 
are evaluated as outputs and have their partial and total sensitivities checked.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/debug_plate_w_pt_mass.bdf")

# Historical reference values for function outputs
FUNC_REFS = {
    "analysis.mass": 68.1,
    "analysis.Ixx": 186.63379667,
    "analysis.ks_vmfailure": 10.92484626,
}

# Inputs to check total sensitivities wrt
wrt = ["mesh.fea_mesh.x_struct0", "dv_struct", "f_struct"]

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
            self.rtol = 1e-7
            self.dh = 1e-50
        else:
            self.rtol = 1e-2
            self.dh = 1e-8

        def assembler_setup(fea_assembler):
            # Add concentrated mass as dv
            fea_assembler.assignMassDV("tip_mass", 5)

        # Callback function used to setup TACS element objects and DVs
        def element_callback(
            dv_num, comp_id, comp_descript, elem_descripts, special_dvs, **kwargs
        ):
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
            con = constitutive.IsoShellConstitutive(
                prop, t=thickness, tNum=dv_num, tlb=min_thickness, tub=max_thickness
            )

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
            problem.setOption("L2Convergence", 1e-99)
            problem.setOption("L2ConvergenceRel", 1e-20)

            # Add TACS Functions
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=ksweight)
            problem.addFunction(
                "Ixx",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[1.0, 0.0, 0.0],
            )

            # Add gravity load
            g = np.array([0.0, 0.0, -9.81]) * 100  # m/s^2
            problem.addInertialLoad(g)

        class Top(Multipoint):
            def setup(self):
                tacs_builder = tacs.mphys.TacsBuilder(
                    mesh_file=bdf_file,
                    assembler_setup=assembler_setup,
                    element_callback=element_callback,
                    problem_setup=problem_setup,
                    check_partials=True,
                    coupled=True,
                    write_solution=False,
                    separate_mass_dvs=True,
                )
                tacs_builder.initialize(self.comm)

                dvs = self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])
                dvs.add_output("dv_struct", 0.01)
                dvs.add_output("dv_mass_tip_mass", 12.5)

                f_size = tacs_builder.get_ndof() * tacs_builder.get_number_of_nodes()
                forces = self.add_subsystem("forces", om.IndepVarComp(), promotes=["*"])
                forces.add_output("f_struct", val=np.ones(f_size), distributed=True)

                self.add_subsystem("mesh", tacs_builder.get_mesh_coordinate_subsystem())
                self.mphys_add_scenario(
                    "analysis", ScenarioStructural(struct_builder=tacs_builder)
                )
                self.connect("mesh.x_struct0", "analysis.x_struct0")
                self.connect("dv_struct", "analysis.dv_struct")
                self.connect("dv_mass_tip_mass", "analysis.dv_mass_tip_mass")
                self.connect("f_struct", "analysis.f_struct")

        prob = om.Problem()
        prob.model = Top()

        return prob

    def setup_funcs(self):
        """
        Create a dict of functions to be tested and a list of inputs
        to test their sensitivities with respect to.
        """
        return FUNC_REFS, wrt
