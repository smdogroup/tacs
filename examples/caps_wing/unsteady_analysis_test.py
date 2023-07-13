"""
Sean Engelstad, Febuary 2023
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example
"""

from tacs import functions, caps2tacs
import numpy as np
from mpi4py import MPI

# -------------------------------------------------------------------------------------------------
# 1: build the tacs aim, egads aim wrapper classes
comm = MPI.COMM_WORLD
tacs_model = caps2tacs.TacsModel.build(csm_file="simple_naca_wing.csm", comm=comm)
tacs_model.egads_aim.set_mesh(
    edge_pt_min=15,
    edge_pt_max=20,
    global_mesh_size=0.25,
    max_surf_offset=0.01,
    max_dihedral_angle=15,
).register_to(tacs_model)

aluminum = caps2tacs.Isotropic.aluminum().register_to(tacs_model)

# setup the thickness design variables + automatic shell properties
nribs = int(tacs_model.get_config_parameter("nribs"))
nspars = int(tacs_model.get_config_parameter("nspars"))
for irib in range(1, nribs + 1):
    caps2tacs.ThicknessVariable(
        caps_group=f"rib{irib}", value=0.05, material=aluminum
    ).register_to(tacs_model)
for ispar in range(1, nspars + 1):
    caps2tacs.ThicknessVariable(
        caps_group=f"spar{ispar}", value=0.05, material=aluminum
    ).register_to(tacs_model)
caps2tacs.ThicknessVariable(
    caps_group="OML", value=0.03, material=aluminum
).register_to(tacs_model)

# add constraints and loads
caps2tacs.PinConstraint("root").register_to(tacs_model)
caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=100).register_to(tacs_model)

# run the pre analysis to build tacs input files
tacs_model.setup(include_aim=True)

# ----------------------------------------------------------------------------------
# 2. Run the TACS unsteady elastic structural analysis, forward + adjoint

# create a transient problem in pyTACS
fea_solver = tacs_model.fea_solver
fea_solver.initialize()
TP = fea_solver.createTransientProblem(
    "sinusoidalWing", tInit=0.0, tFinal=1.0, numSteps=10
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
TP.addFunction("ks_vmfailure", functions.KSFailure, safetyFactor=1.5, ksWeight=50.0)

# run the tacs unsteady analysis, forward and adjoint
print("Starting unsteady tacs wing analysis, takes about 5 minutes")
tacs_funcs = {}
tacs_sens = {}
TP.solve()
TP.evalFunctions(tacs_funcs, evalFuncs=function_names)
TP.evalFunctionsSens(tacs_sens, evalFuncs=function_names)
TP.writeSolution(baseName="tacs_output", outputDir=tacs_model.analysis_dir)

struct_derivs = tacs_sens["ks_vmfailure"]["struct"]
adjoint_TD = struct_derivs[0]

# complex step derivative test
assembler = fea_solver.assembler
xvec = assembler.createDesignVec()
assembler.getDesignVars(xvec)
xarray = xvec.getArray()

# This assumes that the TACS variables are not distributed and are set
# only on the tacs_comm root processor.
h = 1e-30
if comm.rank == 0:
    xarray[0] += 1j * h

assembler.setDesignVars(xvec)


tacs_funcs = {}
TP.evalFunctions(tacs_funcs)
ks_value = tacs_funcs["ks_vmfailure"]
complex_step_TD = ks_value.imag / h

relative_error = (adjoint_TD - complex_step_TD) / complex_step_TD

print(f"approximate derivative = {adjoint_TD}")
print(f"complex step derivative = {complex_step_TD}")
print(f"relative error = {relative_error}")
