"""
Sean Engelstad, November 2023
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example

The purpose of this file is to show how to animate shape variables in the structure.
Once you run this script it will generate .f5 files in the capsStruct/Scratch/tacs work directory.
You'll need to use f5tovtk to convert to *.vtk files and open up each group of shape_var_..vtk files
using the .. shortcut which opens up a group of files at once. For each of these groups of vtk files,
click file -> save animation and save all the png files in a subfolder. Then, either use the GifWriter in the
caps2tacs module or copy the GifWriter script into that directory and use it to make a gif. Repeat for each shape variable.

NOTE : this is meant to be done on one processor, but could be parallelized better in the future.
"""

from tacs import caps2tacs
from mpi4py import MPI

# --------------------------------------------------------------#
# Setup CAPS Problem
# --------------------------------------------------------------#
comm = MPI.COMM_WORLD
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

# register one shape variable rib_a1
rib_a1 = caps2tacs.ShapeVariable("rib_a1", value=1.0).register_to(tacs_model)
caps2tacs.ShapeVariable("rib_a2", value=0.0).register_to(tacs_model)
spar_a1 = caps2tacs.ShapeVariable("spar_a1", value=1.0).register_to(tacs_model)
caps2tacs.ShapeVariable("spar_a2", value=0.0).register_to(tacs_model)

# add constraints and loads
caps2tacs.PinConstraint("root").register_to(tacs_model)
caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=10.0).register_to(
    tacs_model
)

# add analysis functions to the model
caps2tacs.AnalysisFunction.mass().register_to(tacs_model)

# run the pre analysis to build tacs input files
tacs_model.setup(include_aim=True)

shape_var_dict = {
    rib_a1: [_ * 0.1 for _ in range(6, 14)],
    spar_a1: [_ * 0.1 for _ in range(6, 14)],
}
tacs_model.animate_shape_vars(shape_var_dict)
