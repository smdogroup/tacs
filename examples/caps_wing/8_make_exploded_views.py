"""
Sean Engelstad, Febuary 2023
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example

December 1, 2023
Run this with conda-mpirun -n 3 python 8_make_exploded_view.py
    Then process the f5 files which are moved to the exploded folder.
"""

from tacs import caps2tacs
from mpi4py import MPI
import os, shutil

# run a steady elastic structural analysis in TACS using the tacsAIM wrapper caps2tacs submodule
# -------------------------------------------------------------------------------------------------
# 1: build the tacs aim, egads aim wrapper classes

comm = MPI.COMM_WORLD

print(f"proc on rank {comm.rank}")

tacs_model = caps2tacs.TacsModel.build(
    csm_file="large_naca_wing.csm",
    comm=comm,
    problem_name="capsExploded",
    active_procs=[_ for _ in range(3)],
)
tacs_model.mesh_aim.set_mesh(
    edge_pt_min=15,
    edge_pt_max=20,
    global_mesh_size=0.25,
    max_surf_offset=0.01,
    max_dihedral_angle=15,
).register_to(tacs_model)

# set the exploded views to 1,2,3 for each of the 3 separate procs
tacs_aim = tacs_model.tacs_aim
exploded_view = comm.rank + 1
tacs_aim.set_config_parameter("exploded", exploded_view)

aluminum = caps2tacs.Isotropic.aluminum().register_to(tacs_model)

# setup the thickness design variables + automatic shell properties
# choose pattern based dv1 values
nribs = int(tacs_model.get_config_parameter("nribs"))
nspars = int(tacs_model.get_config_parameter("nspars"))
nOML = nribs - 1

if exploded_view == 2:  # only make these internal struct variables on proc 1 tacs model
    for irib in range(1, nribs + 1):
        caps2tacs.ThicknessVariable(
            caps_group=f"rib{irib}", value=0.1 - 0.01 * irib, material=aluminum
        ).register_to(tacs_model)
    for ispar in range(1, nspars + 1):
        caps2tacs.ThicknessVariable(
            caps_group=f"spar{ispar}", value=0.1 - 0.01 * irib, material=aluminum
        ).register_to(tacs_model)

if exploded_view in [
    1,
    3,
]:  # only make these OML struct variables on procs 0,2 tacs model
    for iOML in range(1, nOML + 1):
        caps2tacs.ThicknessVariable(
            caps_group=f"OML{iOML}", value=0.1 - 0.01 * iOML, material=aluminum
        ).register_to(tacs_model)

# add constraints and loads
caps2tacs.PinConstraint("root").register_to(tacs_model)
caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=100).register_to(tacs_model)

# run the pre analysis to build tacs input files
# alternative is to call tacs_aim.setup_aim().pre_analysis() with tacs_aim = tacs_model.tacs_aim
tacs_model.setup()

comm.Barrier()

SPs = tacs_model.createTACSProbs(addFunctions=True)

# make a directory for the exploded views
exploded_dir = os.path.join(os.getcwd(), "exploded")
if not os.path.exists(exploded_dir) and comm.rank == 0:
    os.mkdir(exploded_dir)

# create a mesh of each of the exploded views
names = ["upperOML", "int-struct", "lowerOML"]
for rank in range(3):
    # create tacs problems in the tacs directory of this proc
    SPs = tacs_model.createTACSProbs(root=rank)
    # then write solution using all procs
    for caseID in SPs:
        SPs[caseID].writeSolution(
            baseName=names[rank], outputDir=tacs_model.analysis_dir(proc=rank), number=0
        )

    # move the file to the exploded view directory
    filename = names[rank] + "_000.f5"
    src = os.path.join(tacs_model.analysis_dir(proc=rank), filename)
    dest = os.path.join(exploded_dir, filename)

    if comm.rank == 0:
        shutil.copy(src, dest)

    comm.Barrier()

print("Done creating exploded views for the large_naca_wing.csm..")
print(
    "\tGoto the exploded/ folder and process the f5 files to create an exploded view visualization!"
)
