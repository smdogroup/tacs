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

carbon_fiber = caps2tacs.Orthotropic.carbon_fiber().register_to(tacs_model)

# setup the thickness design variables + automatic shell properties
nribs = int(tacs_model.get_config_parameter("nribs"))
nspars = int(tacs_model.get_config_parameter("nspars"))
# makes unidirectional laminate composite properties automatically
for irib in range(1, nribs + 1):
    caps_group = f"rib{irib}"
    thick = 0.05
    caps2tacs.CompositeProperty(
        caps_group=caps_group,
        ply_materials=[carbon_fiber] * 4,
        ply_thicknesses=[thick / 4] * 4,
        ply_angles=[0, -45, 45, 90],
    ).register_to(tacs_model)
    # caps2tacs.ThicknessVariable(
    #     caps_group=f"rib{irib}",
    #     value=thick,
    # ).register_to(tacs_model)
for ispar in range(1, nspars + 1):
    caps_group = f"spar{ispar}"
    thick = 0.04
    caps2tacs.CompositeProperty(
        caps_group=caps_group,
        ply_materials=[carbon_fiber] * 4,
        ply_thicknesses=[thick / 4] * 4,
        ply_angles=[0, -45, 45, 90],
    ).register_to(tacs_model)
    # caps2tacs.ThicknessVariable(
    #     caps_group=caps_group,
    #     value=thick,
    # ).register_to(tacs_model)

caps_group = "OML"
thick = 0.03
caps2tacs.CompositeProperty(
    caps_group=caps_group,
    ply_materials=[carbon_fiber] * 4,
    ply_thicknesses=[thick / 4] * 4,
    ply_angles=[0, -45, 45, 90],
).register_to(tacs_model)
# caps2tacs.ThicknessVariable(
#     caps_group=caps_group,
#     value=thick,
# ).register_to(tacs_model)

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
if method == 1:  # less detail version
    tacs_model.pre_analysis()
    tacs_model.run_analysis()

    # assembler =
    # sx = np.array([0 for _ in range(8)])
    # sx[0] = 1.0
    # sy = np.array([0 for _ in range(8)])
    # sy[1] = 1.0

    tacs_model.post_analysis()

    print("Tacs model static analysis outputs...\n")
    for func in tacs_model.analysis_functions:
        print(f"func {func.name} = {func.value}")

        for var in tacs_model.variables:
            derivative = func.get_derivative(var)
            print(f"\td{func.name}/d{var.name} = {derivative}")
        print("\n")

elif method == 2:  # longer way that directly uses pyTACS & static problem routines
    SPs = tacs_model.createTACSProbs(addFunctions=True)

    # solve each structural analysis problem (in this case 1)
    tacs_funcs = {}
    tacs_sens = {}
    function_names = tacs_model.function_names
    for caseID in SPs:
        SPs[caseID].solve()
        SPs[caseID].evalFunctions(tacs_funcs, evalFuncs=function_names)
        SPs[caseID].evalFunctionsSens(tacs_sens, evalFuncs=function_names)
        SPs[caseID].writeSolution(
            baseName="tacs_output", outputDir=tacs_model.analysis_dir
        )

    # print the output analysis functions and sensitivities
    print("\nTacs Analysis Outputs...")
    for func_name in function_names:
        # find the associated tacs key (tacs key = (loadset) + (func_name))
        for tacs_key in tacs_funcs:
            if func_name in tacs_key:
                break

        func_value = tacs_funcs[tacs_key].real
        print(f"\tfunctional {func_name} = {func_value}")

        struct_sens = tacs_sens[tacs_key]["struct"]
        print(f"\t{func_name} structDV gradient = {struct_sens}")

        xpts_sens = tacs_sens[tacs_key]["Xpts"]
        print(f"\t{func_name} coordinate derivatives = {xpts_sens}")

        print("\n")
