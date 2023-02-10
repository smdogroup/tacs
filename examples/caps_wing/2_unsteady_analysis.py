"""
Sean Engelstad, January 2023
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example
"""

from tacs import functions, caps2tacs
import numpy as np
from mpi4py import MPI

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
# 2. Run the TACS unsteady elastic structural analysis, forward + adjoint

# create a transient problem in pyTACS
fea_solver = tacs_aim.fea_solver
fea_solver.initialize()
TP = fea_solver.createTransientProblem(
    "sinusoidalWing", tInit=0.0, tFinal=10.0, numSteps=100
)
timeSteps = TP.getTimeSteps()

# get the subcase of the problem and its load ID
for subcase in fea_solver.bdfInfo.subcases.values():
    if "LOAD" in subcase:
        loadsID = subcase["LOAD"][0]

# loop over each load
for itime, time in enumerate(timeSteps):
    # compute the load scale at this time step
    load_scale = np.sin(0.5 * time)

    # add the loads from a static problem
    TP.addLoadFromBDF(itime, loadsID, scale=load_scale)

# Add functions to the transient problem
function_names = ["mass", "ks_vmfailure"]
TP.addFunction("mass", functions.StructuralMass)
TP.addFunction("ks_vmfailure", functions.KSFailure, safetyFactor=1.5, ksWeight=1000.0)

# run the tacs unsteady analysis, forward and adjoint
print("Starting unsteady tacs wing analysis, takes about 5 minutes")
tacs_funcs = {}
tacs_sens = {}
TP.solve()
TP.evalFunctions(tacs_funcs, evalFuncs=function_names)
TP.evalFunctionsSens(tacs_sens, evalFuncs=function_names)
TP.writeSolution(baseName="tacs_output", outputDir=tacs_aim.analysisDir)

# print the function outputs
print(f"\n\nOutputs of TACS unsteady structural analysis...")
for tacs_key in tacs_funcs:
    # find associated function name for tacs_key=loadset+func_name
    for func_name in function_names:
        if func_name in tacs_key:
            break
    print(f"\t{tacs_key} function = {tacs_funcs[tacs_key]}")

    c_struct_sens = tacs_sens[tacs_key]["struct"]
    xpts_sens = tacs_sens[tacs_key]["Xpts"]
    print(f"\tTACS unsteady thickness derivatives = {c_struct_sens}")
    print(f"\tTACS unsteady coordinate derivatives = {xpts_sens}")

# tell the user how to view the structural analysis results
base_name = "tacs_output"
print("\n\nPlease convert the f5 files to vtk by running the f5tovtk shell script...")
print("bash f5tovtk.sh")
print("Then view the results in paraview as follows")
print(f"\t paraview {tacs_aim.analysisDir}/{base_name}_000_..vtk")
print("cd ../../../")
