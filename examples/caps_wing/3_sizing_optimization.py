"""
Sean Engelstad, Febuary 2023
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example
"""

# also requires pySNOPT
from tacs import caps2tacs
import openmdao.api as om
from mpi4py import MPI

# --------------------------------------------------------------#
# Setup CAPS Problem
# --------------------------------------------------------------#
comm = MPI.COMM_WORLD
# can also switch to large_naca_wing.csm file here if you want and it will automatically update the DVs
tacs_model = caps2tacs.TacsModel.build(csm_file="large_naca_wing.csm", comm=comm)
tacs_model.mesh_aim.set_mesh(  # need a refined-enough mesh for the derivative test to pass
    edge_pt_min=15,
    edge_pt_max=20,
    global_mesh_size=0.01,
    max_surf_offset=0.01,
    max_dihedral_angle=5,
).register_to(
    tacs_model
)

aluminum = caps2tacs.Isotropic.aluminum().register_to(tacs_model)

# setup the thickness design variables + automatic shell properties
nribs = int(tacs_model.get_config_parameter("nribs"))
nspars = int(tacs_model.get_config_parameter("nspars"))
nOML = nribs - 1
for irib in range(1, nribs + 1):
    caps2tacs.ThicknessVariable(
        caps_group=f"rib{irib}", value=0.03, material=aluminum
    ).register_to(tacs_model)
for ispar in range(1, nspars + 1):
    caps2tacs.ThicknessVariable(
        caps_group=f"spar{ispar}", value=0.03, material=aluminum
    ).register_to(tacs_model)
for iOML in range(1, nOML + 1):
    caps2tacs.ThicknessVariable(
        caps_group=f"OML{iOML}", value=0.1, material=aluminum
    ).register_to(tacs_model)

# add constraints and loads
caps2tacs.PinConstraint("root").register_to(tacs_model)
caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=10).register_to(tacs_model)

# add in analysis functions into the tacs aim
caps2tacs.AnalysisFunction.mass().register_to(tacs_model)
caps2tacs.AnalysisFunction.ksfailure(safetyFactor=1.5, ksWeight=50.0).register_to(
    tacs_model
)

# run the tacs model setup
tacs_model.setup(include_aim=True)

# --------------------------------------------------------------------------#
# Setup OpenMDAO Problem and Perform Sizing Optimization on the Wing
# --------------------------------------------------------------------------#

# setup the OpenMDAO Problem object
prob = om.Problem()

# Create the OpenMDAO component using the built-in Tacs Static Component
tacs_system = caps2tacs.TacsStaticComponent(
    tacs_model=tacs_model, write_f5=True, track_history=True
)
prob.model.add_subsystem("tacsSystem", tacs_system)

# setup the optimizer settings # COBYLA for auto-FDing
prob.driver = om.ScipyOptimizeDriver(optimizer="SLSQP", tol=1.0e-9, disp=True)

# add design variables to the model
for thick_var in tacs_model.thickness_variables:
    prob.model.add_design_var(
        f"tacsSystem.{thick_var.name}", lower=0.001, upper=0.5, scaler=100.0
    )

# add objectives & constraints to the model
prob.model.add_objective("tacsSystem.mass", scaler=1.0e-2)
prob.model.add_constraint("tacsSystem.ksfailure", upper=0.267)

# Start the optimization
print("\n==> Starting optimization...")
prob.setup()

debug = False
if debug:
    print("Checking partials...", flush=True)
    prob.check_partials(compact_print=True)

else:
    prob.run_driver()

    # report the final optimal design
    design_hdl = tacs_system._design_hdl
    for thick_var in tacs_model.thickness_variables:
        opt_value = prob.get_val(f"tacsSystem.{thick_var.name}")
        design_hdl.write(f"\t{thick_var.name} = {opt_value}")

    prob.cleanup()
    # close the design hdl file
    design_hdl.close()
