"""
Sean Engelstad, November 2023
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example

To call mpi version use (e.g. with 4 procs):
    mpirun -n 4 python 7_multiproc_post_analysis.py
* For local machine I use conda-mpirun alias instead of mpirun to point to my mpich from conda-forge mpi4py

Timing data comparison for different # of procs to show speedup 
    of optimization iterations with shape derivative on multi-procs:
    NOTE : I still run with 4 procs (so TACS analysis takes same amount of time)
    but use different # of active procs in the TacsModel input.
--------------------------------------------------------------------
1 proc - 
2 proc - 
3 proc - 
4 proc -
"""

import time
from tacs import caps2tacs
from mpi4py import MPI

start_time = time.time()

# --------------------------------------------------------------#
# Setup CAPS Problem
# --------------------------------------------------------------#
comm = MPI.COMM_WORLD
tacs_model = caps2tacs.TacsModel.build(csm_file="large_naca_wing.csm", comm=comm, active_procs=[0])
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
tacs_model.post_analysis()

# update the shape variables in the design
print(f"design 2 on rank {comm.rank}", flush=True)
tacs_model.update_design(input_dict={"rib_a1" : 0.8, "rib_a2" : 0.2, "spar_a1" : 0.9, "spar_a2" : -0.2})

# run analysis with the second design
tacs_model.pre_analysis()
tacs_model.run_analysis()
tacs_model.post_analysis()

# update the shape variables in the design
print(f"design 3 on rank {comm.rank}", flush=True)
tacs_model.update_design(input_dict={"rib_a1" : 1.2, "rib_a2" : -0.2, "spar_a1" : 1.0, "spar_a2" : 0.2})

# run analysis with the third design
tacs_model.pre_analysis()
tacs_model.run_analysis()
tacs_model.post_analysis()

# print out the derivatives => make sure shape derivatives are nonzero
comm.Barrier()
print("Tacs model static analysis outputs...\n")
for func in tacs_model.analysis_functions:
    print(f"func {func.name} = {func.value}")

    for var in tacs_model.variables:
        derivative = func.get_derivative(var)
        print(f"\td{func.name}/d{var.name} = {derivative}")
    print("\n")

# print out timing results (compare serial to parallel instances of tacsAIM)
mins_elapsed = (time.time() - start_time) / 60.0
naims = len(tacs_model.active_procs)
if comm.rank == 0:
    print(f"Elapsed time with {naims} procs is {mins_elapsed} mins using 4 shape vars.")