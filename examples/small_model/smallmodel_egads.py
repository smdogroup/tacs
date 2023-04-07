from tacs import caps2tacs
from mpi4py import MPI

import os
import numpy as np
import sys


# run a steady elastic structural analysis in TACS using the tacsAIM wrapper caps2tacs submodule
# -------------------------------------------------------------------------------------------------
# 1: build the tacs aim, egads aim wrapper classes

comm = MPI.COMM_WORLD
tacs_model = caps2tacs.TacsModel.build(csm_file="small.csm", comm=comm, mesh="egads")
tacs_model.mesh_aim.set_mesh(
    edge_pt_min=20,
    edge_pt_max=30,
    global_mesh_size=0.25,
    max_surf_offset=0.01,
    max_dihedral_angle=15,
).register_to(tacs_model)

aluminumT6 = caps2tacs.Material(
    name="aluminumT6",
    material_type="Isotropic",
    young_modulus=1e7,
    poisson_ratio=0.33,
    density=0.0975,
    tension_allow=45000,
    compression_allow=45000,
    yield_allow=40000,
    thermExpCoeff=13.1e-6,
).register_to(tacs_model)


# setup the thickness design variables + automatic shell properties
caps2tacs.ThicknessVariable(
    caps_group="frame", value=0.005, material=aluminumT6
).register_to(tacs_model)

caps2tacs.ThicknessVariable(
    caps_group="pylon1", value=0.005, material=aluminumT6
).register_to(tacs_model)

caps2tacs.ThicknessVariable(
    caps_group="pylon2", value=0.005, material=aluminumT6
).register_to(tacs_model)

caps2tacs.ThicknessVariable(
    caps_group="pylon3", value=0.005, material=aluminumT6
).register_to(tacs_model)

caps2tacs.ThicknessVariable(
    caps_group="pylon4", value=0.005, material=aluminumT6
).register_to(tacs_model)

caps2tacs.ThicknessVariable(
    caps_group="endcap", value=0.005, material=aluminumT6
).register_to(tacs_model)

caps2tacs.ThicknessVariable(
    caps_group="endcap1", value=0.005, material=aluminumT6
).register_to(tacs_model)

caps2tacs.ThicknessVariable(
    caps_group="fixed", value=0.005, material=aluminumT6
).register_to(tacs_model)

# add constraints and loads
# for icon in range(4):
# caps2tacs.PinConstraint(f"fixed{icon}").register_to(tacs_model)

caps2tacs.PinConstraint("fixed").register_to(tacs_model)
caps2tacs.Pressure("pylon1", force=14.71).register_to(tacs_model)
caps2tacs.Pressure("pylon2", force=14.71).register_to(tacs_model)
caps2tacs.Pressure("pylon3", force=14.71).register_to(tacs_model)
caps2tacs.Pressure("pylon4", force=14.71).register_to(tacs_model)
# caps2tacs.GridForce("endcap1", direction=[0, 0, 1.0], magnitude=14.17).register_to(tacs_model)

# add analysis functions to the model
caps2tacs.AnalysisFunction.ksfailure(ksWeight=50.0, safetyFactor=1.5).register_to(
    tacs_model
)
caps2tacs.AnalysisFunction.mass().register_to(tacs_model)

# run the pre analysis to build tacs input files
# alternative is to call tacs_aim.setup_aim().pre_analysis() with tacs_aim = tacs_model.tacs_aim
tacs_model.setup(include_aim=True)


# tacs_aim = tacs_model.tacs_aim
# print(f"ESP/CAPS loads = {tacs_aim.aim.input.Load}", flush=True)
# sys.exit()

# ----------------------------------------------------------------------------------
# 2. Run the TACS steady elastic structural analysis, forward + adjoint

# choose method 1 or 2 to demonstrate the analysis : 1 - fewer lines, 2 - more lines

# show both ways of performing the structural analysis in more or less detail

tacs_model.pre_analysis()
tacs_model.run_analysis()
tacs_model.post_analysis()

print("Tacs model static analysis outputs...\n")

for func in tacs_model.analysis_functions:
    print(f"func {func.name} = {func.value}")

    for var in tacs_model.variables:
        derivative = func.get_derivative(var)
        print(f"\td{func.name}/d{var.name} = {derivative}")
        print("\n")
