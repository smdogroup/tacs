"""
Sean Engelstad, Febuary 2023
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example
"""

from tacs import caps2tacs
from mpi4py import MPI
import openmdao.api as om

# --------------------------------------------------------------#
# Setup CAPS Problem
# --------------------------------------------------------------#
comm = MPI.COMM_WORLD
caps_struct = caps2tacs.CapsStruct.build(csm_file="simple_naca_wing.csm", comm=comm)
tacs_model = caps_struct.tacs_model
caps_struct.egads_aim.set_mesh(  # need a refined-enough mesh for the derivative test to pass
    edge_pt_min=15,
    edge_pt_max=20,
    global_mesh_size=0.1,
    max_surf_offset=0.01,
    max_dihedral_angle=5,
).register_to(
    tacs_model
)

aluminum = caps2tacs.Isotropic.aluminum().register_to(tacs_model)

# setup the thickness design variables + automatic shell properties
nribs = int(tacs_model.get_config_parameter("nribs"))
nspars = int(tacs_model.get_config_parameter("nspars"))
for irib in range(1, nribs + 1):
    caps2tacs.ShellProperty(
        caps_group=f"rib{irib}", material=aluminum, membrane_thickness=0.05
    ).register_to(tacs_model)
for ispar in range(1, nspars + 1):
    caps2tacs.ShellProperty(
        caps_group=f"spar{ispar}", material=aluminum, membrane_thickness=0.05
    ).register_to(tacs_model)
caps2tacs.ShellProperty(
    caps_group="OML", material=aluminum, membrane_thickness=0.03
).register_to(tacs_model)

# register one shape variable rib_a1
caps2tacs.ShapeVariable("rib_a1").register_to(tacs_model)
caps2tacs.ShapeVariable("rib_a2").register_to(tacs_model)

# add constraints and loads
caps2tacs.PinConstraint("root").register_to(tacs_model)
caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=100).register_to(tacs_model)

# add analysis functions to the model
caps2tacs.AnalysisFunction.mass().register_to(tacs_model)

# run the pre analysis to build tacs input files
tacs_model.setup(include_aim=True)

# ---------------------------------------------------------------------------------#
# Setup OpenMDAO Problem
# ---------------------------------------------------------------------------------#

# setup the OpenMDAO problem object
om_problem = om.Problem()

# Create the OpenMDAO component
tacs_system = caps2tacs.TacsComponent(tacs_model=tacs_model)
om_problem.model.add_subsystem("tacsSystem", tacs_system)

# setup the optimization
om_problem.driver = om.ScipyOptimizeDriver()
om_problem.driver.options["optimizer"] = "SLSQP"
om_problem.driver.options["tol"] = 1.0e-9
om_problem.driver.options["disp"] = True

# add design variables to the model
om_problem.model.add_design_var("tacsSystem.rib_a1", lower=0.6, upper=1.4)
om_problem.model.add_design_var("tacsSystem.rib_a2", lower=-0.3, upper=0.3)

# add objectives to the model
om_problem.model.add_objective("tacsSystem.mass")
# om_problem.model.add_objective('tacsSystem.ks_vmfailure')

# Start the optimization
print("\n==> Starting Optimization...")
om_problem.setup()
om_problem.run_driver()

tacs_system._design_hdl.write("--> Optimized values:\n")
rib_a1_value = om_problem.get_val("tacsSystem.rib_a1")
rib_a2_value = om_problem.get_val("tacsSystem.rib_a2")
tacs_system._design_hdl.write(f"\trib_a1 = {rib_a1_value}\n")
tacs_system._design_hdl.write(f"\trib_a2 = {rib_a2_value}\n")
# close the design hdl file and cleanup
tacs_system._design_hdl.close()
om_problem.cleanup()
