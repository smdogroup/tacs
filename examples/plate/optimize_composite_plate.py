"""
Mass minimization of uCRM wingbox subject to a constant vertical force
"""
import os

import openmdao.api as om
import numpy as np
from mphys import Multipoint
from mphys.scenario_structural import ScenarioStructural

from tacs import elements, constitutive, functions
from tacs.mphys import TacsBuilder

# BDF file containing mesh
bdf_file = os.path.join(os.path.dirname(__file__), "partitioned_plate.bdf")

# Material properties
rho = 1550.0
E1 = 54e9
E2 = 18e9
nu12 = 0.25
G12 = 9e9
G13 = 9e9
Xt = 2410.0e6
Xc = 1040.0e6
Yt = 73.0e6
Yc = 173.0e6
S12 = 71.0e6

# Shell thickness
ply_thickness = 1.25e-3  # m
plate_thickness = 0.05  # m
tMin = 0.002  # m
tMax = 0.05  # m

# Ply angles/initial ply fractions
ply_angles = np.deg2rad([0.0, 45.0, -45.0, 90.0])
ply_fractions = np.array([0.25, 0.25, 0.25, 0.25])

# Pressure load to apply to plate
P = 100e3


# Callback function used to setup TACS element objects and DVs
def element_callback(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
    # Create ply object
    ortho_prop = constitutive.MaterialProperties(
        rho=rho,
        E1=E1,
        E2=E2,
        nu12=nu12,
        G12=G12,
        G13=G13,
        G23=G13,
        Xt=Xt,
        Xc=Xc,
        Yt=Yt,
        Yc=Yc,
        S12=S12,
    )
    ortho_ply = constitutive.OrthotropicPly(ply_thickness, ortho_prop)
    # Create the layup list (one for each angle)
    ortho_layup = [ortho_ply, ortho_ply, ortho_ply, ortho_ply]
    # Assign each ply fraction a unique DV
    ply_fraction_dv_nums = np.array(
        [dvNum, dvNum + 1, dvNum + 2, dvNum + 3], dtype=np.intc
    )
    # Create smeared stiffness object based on ply angles/fractions
    con = constitutive.SmearedCompositeShellConstitutive(
        ortho_layup,
        plate_thickness,
        ply_angles,
        ply_fractions,
        ply_fraction_dv_nums=ply_fraction_dv_nums,
    )

    # Define reference axis to define local 0 deg direction
    refAxis = np.array([1.0, 0.0, 0.0])
    transform = elements.ShellRefAxisTransform(refAxis)

    # Pass back the appropriate tacs element object
    elem = elements.Quad4Shell(transform, con)

    return elem


def problem_setup(scenario_name, fea_assembler, problem):
    """
    Helper function to add fixed forces and eval functions
    to structural problems used in tacs builder
    """

    # Add TACS Functions
    problem.addFunction("compliance", functions.Compliance)

    # Add forces to static problem
    allComponents = fea_assembler.selectCompIDs()
    problem.addPressureToComponents(allComponents, P)


def constraint_setup(scenario_name, fea_assembler, constraint_list):
    """
    Helper function to setup tacs constraint classes
    """
    constr = fea_assembler.createDVConstraint("ply_fractions")
    allComponents = fea_assembler.selectCompIDs()
    constr.addConstraint(
        "sum", allComponents, dvIndices=[0, 1, 2, 3], dvWeights=[1.0, 1.0, 1.0, 1.0]
    )
    constraint_list.append(constr)


class PlateModel(Multipoint):
    def setup(self):
        struct_builder = TacsBuilder(
            mesh_file=bdf_file,
            element_callback=element_callback,
            problem_setup=problem_setup,
            constraint_setup=constraint_setup,
            coupled=False,
            check_partials=True,
        )
        struct_builder.initialize(self.comm)
        dv_array = struct_builder.get_initial_dvs()

        dvs = self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])
        dvs.add_output("dv_struct", dv_array)

        self.add_subsystem("mesh", struct_builder.get_mesh_coordinate_subsystem())
        self.mphys_add_scenario(
            "pressure_load", ScenarioStructural(struct_builder=struct_builder)
        )
        self.mphys_connect_scenario_coordinate_source("mesh", "pressure_load", "struct")

        self.connect("dv_struct", "pressure_load.dv_struct")


################################################################################
# OpenMDAO setup
################################################################################

prob = om.Problem()
prob.model = PlateModel()
model = prob.model

# Declare design variables, objective, and constraint
model.add_design_var("dv_struct", lower=0.0, upper=1.0)
model.add_objective("pressure_load.compliance")
model.add_constraint("pressure_load.ply_fractions.sum", equals=1.0, linear=True)

# Configure optimizer
prob.driver = om.ScipyOptimizeDriver(debug_print=["objs", "nl_cons"], maxiter=100)
prob.driver.options["optimizer"] = "SLSQP"

# Setup OpenMDAO problem
prob.setup()

# Output N2 representation of OpenMDAO model
om.n2(prob, show_browser=False, outfile="tacs_struct.html")

# Run optimization
prob.run_driver()
