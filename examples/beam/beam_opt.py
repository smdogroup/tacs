"""
This example demonstrates TACS structural optimization capabilities.
The beam model that we will be using for this problem is a rectangular beam,
cantilevered, with a shear load applied at the tip. The beam is discretized using
100 beam elements along it's span.

The optimization problem is as follows:
Minimize the mass of the beam with respect to the thickness of each element,
subject to a max stress constraint dictated by the materials yield stress.

By realizing that the stress at any spanwise cross-section in the beam
can be found independently using:
    sigma(x,y) = y*M(x)/I
An analytical solution for this problem can be shown to be:
    t(x) = sqrt(6*V*(L-x)/(w*sigma_y))

The optimization is setup using TACS' MPHYS module, which acts as a wrapper
for OpenMDAO.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

import openmdao.api as om

from mphys import Multipoint
from mphys.scenario_structural import ScenarioStructural
from tacs.mphys import TacsBuilder
from tacs import elements, constitutive, functions

bdf_file = os.path.join(os.path.dirname(__file__), "beam_opt.bdf")

# Beam thickness (initial)
t = 0.05  # m
# Beam width
w = 0.05  # m
# Length of beam
L = 1.0

# Material properties
rho = 2500.0  # density kg/m^3
E = 70.0e9  # Young's modulus (Pa)
nu = 0.0  # Poisson's ratio
ys = 350e6  # yield stress

# Shear force applied at tip
V = 1e3

# Callback function used to setup TACS element objects and DVs
def element_callback(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
    # Setup (isotropic) property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    # Set one thickness dv for every property group
    con = constitutive.IsoRectangleBeamConstitutive(prop, t=t, w=w, tNum=dvNum)

    refAxis = np.array([0.0, 1.0, 0.0])
    transform = elements.BeamRefAxisTransform(refAxis)

    # Pass back the appropriate tacs element object
    elem = elements.Beam2(transform, con)
    return elem


def problem_setup(scenario_name, fea_assembler, problem):
    """
    Helper function to add fixed forces and eval functions
    to structural problems used in tacs builder
    """

    # Add TACS Functions
    problem.addFunction("mass", functions.StructuralMass)
    problem.addFunction(
        "ks_vmfailure", functions.KSFailure, safetyFactor=1.0, ksWeight=100.0
    )

    # Add forces to static problem
    problem.addLoadToNodes(101, [0.0, V, 0.0, 0.0, 0.0, 0.0], nastranOrdering=True)


class Top(Multipoint):
    def setup(self):
        tacs_options = {
            "element_callback": element_callback,
            "problem_setup": problem_setup,
            "mesh_file": bdf_file,
        }

        # Initialize MPHYS builder for TACS
        struct_builder = TacsBuilder(tacs_options, coupled=False, write_solution=False)
        struct_builder.initialize(self.comm)
        dv_array = struct_builder.get_initial_dvs()

        # Create independent variable component to control design variables for problem
        dvs = self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])
        dvs.add_output("dv_struct", dv_array)

        # Add mesh component
        self.add_subsystem("mesh", struct_builder.get_mesh_coordinate_subsystem())
        self.mphys_add_scenario(
            "tip_shear", ScenarioStructural(struct_builder=struct_builder)
        )
        self.mphys_connect_scenario_coordinate_source("mesh", "tip_shear", "struct")

        # Connect dv component to input of structural scenario
        self.connect("dv_struct", "tip_shear.dv_struct")


################################################################################
# OpenMDAO setup
################################################################################
# Instantiate OpenMDAO problem
prob = om.Problem()
prob.model = Top()
model = prob.model

# Declare design variables, objective, and constraint
model.add_design_var("dv_struct", lower=0.001, upper=0.1, scaler=100.0)
model.add_objective("tip_shear.mass", scaler=1.0)
model.add_constraint("tip_shear.ks_vmfailure", lower=0.0, upper=1.0, scaler=1.0)

# Configure optimizer
prob.driver = om.ScipyOptimizeDriver(debug_print=["objs", "nl_cons"], maxiter=1000)
prob.driver.options["optimizer"] = "SLSQP"

# Setup OpenMDAO problem
prob.setup()

# Run optimization
prob.run_driver()

# Get optimized solution variables
x = prob.get_val("mesh.x_struct0", get_remote=True)[:-3:3]
t_opt = prob["dv_struct"]
m_opt = prob["tip_shear.mass"]

# Get analytical solution
t_exact = np.sqrt(6 * (L - x) * V / w / ys)

if __name__ == "__main__" and prob.comm.size == 1:
    # Output N2 representation of OpenMDAO model
    om.n2(prob, show_browser=False, outfile="beam_opt_n2.html")

    # Plot results for solution
    plt.plot(x, t_opt, "o", x, t_exact)
    plt.legend(["optimized", "analytical"])
    plt.ylabel("t(x)")
    plt.xlabel("x")
    plt.title("Optimal beam thickness profile")
    plt.show()
