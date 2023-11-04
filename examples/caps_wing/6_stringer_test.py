"""
Sean Engelstad, Febuary 2023
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example
"""

from tacs import caps2tacs
from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt

# run a steady elastic structural analysis in TACS using the tacsAIM wrapper caps2tacs submodule
# -------------------------------------------------------------------------------------------------
# 1: build the tacs aim, egads aim wrapper classes

comm = MPI.COMM_WORLD
tacs_model = caps2tacs.TacsModel.build(csm_file="simple_naca_wing.csm", comm=comm)
tacs_model.mesh_aim.set_mesh(
    edge_pt_min=15,
    edge_pt_max=20,
    global_mesh_size=0.25,
    max_surf_offset=0.01,
    max_dihedral_angle=15,
).register_to(tacs_model)

aluminum = caps2tacs.Isotropic.aluminum()
aluminum_w_stringer = caps2tacs.Orthotropic.smeared_stringer(
    aluminum, area_ratio=0.5
).register_to(tacs_model)

# setup the thickness design variables + automatic shell properties
nribs = int(tacs_model.get_config_parameter("nribs"))
nspars = int(tacs_model.get_config_parameter("nspars"))
for irib in range(1, nribs + 1):
    caps2tacs.CompositeProperty.one_ply(
        caps_group=f"rib{irib}",
        material=aluminum_w_stringer,
        ply_angle=0,
        thickness=0.05,  # * 0.5
    ).register_to(tacs_model)
for ispar in range(1, nspars + 1):
    caps2tacs.CompositeProperty.one_ply(
        caps_group=f"spar{ispar}",
        material=aluminum_w_stringer,
        ply_angle=0,
        thickness=0.05,  # * 0.5
    ).register_to(tacs_model)
caps2tacs.CompositeProperty.one_ply(
    caps_group="OML", material=aluminum_w_stringer, ply_angle=0, thickness=0.03  # * 0.5
).register_to(tacs_model)

# add constraints and loads
caps2tacs.PinConstraint("root").register_to(tacs_model)
caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=100).register_to(tacs_model)

# add analysis functions to the model
caps2tacs.AnalysisFunction.ksfailure(ksWeight=50.0, safetyFactor=1.5).register_to(
    tacs_model
)
caps2tacs.AnalysisFunction.mass().register_to(tacs_model)

caps2tacs.ShapeVariable("rib_a1", value=1.0).register_to(tacs_model)

# run the pre analysis to build tacs input files
# alternative is to call tacs_aim.setup_aim().pre_analysis() with tacs_aim = tacs_model.tacs_aim
tacs_model.setup(include_aim=True)

# ----------------------------------------------------------------------------------
# 2. Run the TACS steady elastic structural analysis, forward + adjoint

# choose method 1 or 2 to demonstrate the analysis : 1 - fewer lines, 2 - more lines
method = 1

# show both ways of performing the structural analysis in more or less detail
tacs_model.pre_analysis()
tacs_model.run_analysis()

# view failure profile
# assembler =
# sx = np.array([0 for _ in range(8)])
# sx[0] = 1.0
# sy = np.array([0 for _ in range(8)])
# sy[1] = 1.0

# it builds the wrong Composite constitutive model namely CompositeShellConstitutive
# so that we can't easily do derivatives with the default tacs callback
# instead in funtofem we can make a custom callback there.
# thickness variables for SmearedCompositeShellConstitutive don't work
tacs_model.post_analysis()

print("Tacs model static analysis outputs...\n")
for func in tacs_model.analysis_functions:
    print(f"func {func.name} = {func.value}")

    for var in tacs_model.variables:
        derivative = func.get_derivative(var)
        print(f"\td{func.name}/d{var.name} = {derivative}")
    print("\n")
