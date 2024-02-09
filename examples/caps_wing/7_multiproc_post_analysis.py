"""
Sean Engelstad, November 2023
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example

To call mpi version use (e.g. with 4 procs):
    mpirun -n 4 python 7_multiproc_post_analysis.py
* For local machine I use conda-mpirun alias instead of mpirun to point to my mpich from conda-forge mpi4py

Timing data comparison for different # of procs to show speedup 
    of optimization iterations with shape derivative on multi-procs:
    but use different # of active procs in the TacsModel inputs
The # of parallel tacsAIMs is controlled through the active_procs input to the TacsModel i.e. if len(active_procs) == 3
    then there are 3 active tacsAIMs running in parallel. The post_analysis() scales by the # of shape variables due to finite
    differencing, so by running parallel tacsAIMs and distributing the shape variable FD responsibility, we speed up the post_analysis.
--------------------------------------------------------------------
post_analysis time (4 shape variables, 4 procs)
1 tacsAIM  - 2.611 mins
2 tacsAIMs - 1.252 mins
3 tacsAIMs - 1.157 mins
4 tacsAIMs - 0.591 mins

NOTE : with 2 tacsAIMs - both tacsAIMs have 2 shape vars while 
    with 3 tacsAIMs - 2 AIMs have 1 shape var and 1 has 2 shape vars
    so the runtime for 2 or 3 instances is the same. It's not until 4 instances
    that we reach full 4x speedup.
"""

import time
from tacs import caps2tacs
from mpi4py import MPI

start_time1 = time.time()

# --------------------------------------------------------------#
# Setup CAPS Problem
# --------------------------------------------------------------#
comm = MPI.COMM_WORLD
tacs_model = caps2tacs.TacsModel.build(
    csm_file="large_naca_wing.csm", comm=comm, active_procs=[0, 1]
)
tacs_model.mesh_aim.set_mesh(  # need a refined-enough mesh for the derivative test to pass
    edge_pt_min=15,
    edge_pt_max=20,
    global_mesh_size=0.01,
    max_surf_offset=0.01,
    max_dihedral_angle=5,
).register_to(
    tacs_model
)
tacs_aim = tacs_model.tacs_aim

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

# register one shape variable rib_a1
caps2tacs.ShapeVariable("rib_a1", value=1.0).register_to(tacs_model)
caps2tacs.ShapeVariable("rib_a2", value=0.0).register_to(tacs_model)
caps2tacs.ShapeVariable("spar_a1", value=1.0).register_to(tacs_model)
caps2tacs.ShapeVariable("spar_a2", value=0.0).register_to(tacs_model)

# add constraints and loads
caps2tacs.PinConstraint("root").register_to(tacs_model)
caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=10.0).register_to(
    tacs_model
)

# add analysis functions to the model
caps2tacs.AnalysisFunction.mass().register_to(tacs_model)
caps2tacs.AnalysisFunction.ksfailure(safetyFactor=1.5, ksWeight=50.0).register_to(
    tacs_model
)

# run the pre analysis to build tacs input files
print(f"design 1 - setup tacs model on rank {comm.rank}", flush=True)
tacs_model.setup(include_aim=True)

# run analysis with one design
print(f"design 1  on rank {comm.rank}", flush=True)
tacs_model.pre_analysis()
tacs_model.run_analysis()
comm.Barrier()
start_time = time.time()
tacs_model.post_analysis()
comm.Barrier()
post_analysis_time = time.time() - start_time

# print out the derivatives => make sure shape derivatives are nonzero
comm.Barrier()
print("Tacs model static analysis outputs...\n")
for func in tacs_model.analysis_functions:
    print(f"func {func.name} = {func.value}")

    for var in tacs_model.variables:
        derivative = func.get_derivative(var)
        print(f"\td{func.name}/d{var.name} = {derivative}")
    print("\n", flush=True)

# print out timing results (compare serial to parallel instances of tacsAIM)
mins_elapsed = (time.time() - start_time1) / 60.0
naims = len(tacs_model.active_procs)
if comm.rank == 0:
    print(f"Elapsed time with {naims} procs is {mins_elapsed} mins using 4 shape vars.")
    print(f"Also each post analysis took {post_analysis_time/60.0:.4f} min", flush=True)
