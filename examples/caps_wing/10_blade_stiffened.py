"""
Sean Engelstad, March 2024
GT SMDO Lab, Dr. Graeme Kennedy

Based off of Alasdair Christison Gray's Mach wing blade stiffened example in TACS : examples/mach_tutorial_wing
    Can only do analysis in caps2tacs module in TACS need FUNtoFEM + casp2tacs to perform full optimization
"""

from tacs import caps2tacs
import openmdao.api as om
from mpi4py import MPI
import numpy as np

from _blade_callback import blade_elemCallBack

tacs_dtype = TACS.dtype
comm = MPI.COMM_WORLD

# --------------------------------------------------------------#
# Setup CAPS Problem and FUNtoFEM model
# --------------------------------------------------------------#

# define the Tacs model
tacs_model = caps2tacs.TacsModel.build(
    csm_file="large_naca_wing.csm", comm=comm, callback=blade_elemCallBack
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

# setup the thickness design variables + automatic shell properties
# using Composite functions, this part has to go after all funtofem variables are defined...
nribs = int(tacs_model.get_config_parameter("nribs"))
nspars = int(tacs_model.get_config_parameter("nspars"))
nOML = nribs - 1

null_material = caps2tacs.Orthotropic.null().register_to(tacs_model)

for irib in range(1, nribs + 1):
    name = f"rib{irib}"
    caps2tacs.CompositeProperty.null(name, null_material).register_to(tacs_model)
    # caps2tacs.ThicknessVariable(
    #     caps_group=f"rib{irib}", value=init_thickness
    # ).register_to(tacs_model)

for ispar in range(1, nspars + 1):
    name = f"spar{ispar}"
    caps2tacs.CompositeProperty.null(name, null_material).register_to(tacs_model)
    # caps2tacs.ThicknessVariable(
    #     caps_group=f"rib{irib}", value=init_thickness
    # ).register_to(tacs_model)

for iOML in range(1, nOML + 1):
    name = f"OML{iOML}"
    caps2tacs.CompositeProperty.null(name, null_material).register_to(tacs_model)
    # caps2tacs.ThicknessVariable(
    #     caps_group=f"rib{irib}", value=init_thickness
    # ).register_to(tacs_model)

# add constraints and loads
caps2tacs.PinConstraint("root").register_to(tacs_model)
caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=10).register_to(tacs_model)

# add in analysis functions into the tacs aim
caps2tacs.AnalysisFunction.mass().register_to(tacs_model)
caps2tacs.AnalysisFunction.ksfailure(safetyFactor=1.5, ksWeight=50.0).register_to(
    tacs_model
)

# run the tacs model setup
tacs_model.setup(include_aim=True)

# build the solver manager, no tacs interface since built for each new shape
# in the tacs driver
tacs_model.run_analysis()
