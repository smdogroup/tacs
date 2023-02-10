"""
Sean Engelstad, January 2023
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example
"""

from tacs import caps2tacs, functions
from mpi4py import MPI

# run a steady elastic structural analysis in TACS using the tacsAIM wrapper caps2tacs submodule
# -------------------------------------------------------------------------------------------------
# 1: build the tacs aim, egads aim wrapper classes

comm = MPI.COMM_WORLD
caps_struct = caps2tacs.CapsStruct.build(csm_file="simple_naca_wing.csm", comm=comm)
tacs_aim = caps_struct.tacs_aim
caps_struct.egads_aim.set_mesh(
    edge_pt_min=15,
    edge_pt_max=20,
    global_mesh_size=0.25,
    max_surf_offset=0.01,
    max_dihedral_angle=15,
).register_to(tacs_aim)

aluminum = caps2tacs.Isotropic.aluminum().register_to(tacs_aim)

# setup the thickness design variables + automatic shell properties
nribs = int(tacs_aim.get_config_parameter("nribs"))
nspars = int(tacs_aim.get_config_parameter("nspars"))
for irib in range(1, nribs + 1):
    caps2tacs.ThicknessVariable(
        caps_group=f"rib{irib}", value=0.05, material=aluminum
    ).register_to(tacs_aim)
for ispar in range(1, nspars + 1):
    caps2tacs.ThicknessVariable(
        caps_group=f"spar{ispar}", value=0.05, material=aluminum
    ).register_to(tacs_aim)
caps2tacs.ThicknessVariable(
    caps_group="OML", value=0.03, material=aluminum
).register_to(tacs_aim)

# add constraints and loads
caps2tacs.PinConstraint("root").register_to(tacs_aim)
caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=100).register_to(tacs_aim)

# run the pre analysis to build tacs input files
tacs_aim.setup_aim().pre_analysis()

# ----------------------------------------------------------------------------------
# 2. Run the TACS steady elastic structural analysis, forward + adjoint

SPs = tacs_aim.createTACSProbs()

# add mass and aggregated stress constraint functions
for caseID in SPs:
    SPs[caseID].addFunction("mass", functions.StructuralMass)
    SPs[caseID].addFunction(
        "ks_vmfailure", functions.KSFailure, safetyFactor=1.5, ksWeight=1000.0
    )

# solve each structural analysis problem (in this case 1)
tacs_funcs = {}
tacs_sens = {}
base_name = "tacs_output"
function_names = ["mass", "ks_vmfailure"]
for caseID in SPs:
    SPs[caseID].solve()
    SPs[caseID].evalFunctions(tacs_funcs, evalFuncs=function_names)
    SPs[caseID].evalFunctionsSens(tacs_sens, evalFuncs=function_names)
    SPs[caseID].writeSolution(baseName=base_name, outputDir=tacs_aim.analysis_dir)

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
