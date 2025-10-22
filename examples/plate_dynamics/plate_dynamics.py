#!/usr/bin/python
"""
This python script is used to simulate the dynamics of a flexible
plate subject to specified initial and boundary conditions.

The initial conditions such as the velocity, angular velocity and
accerations can be specified in the file as GibbsVectors. The boundary
conditions come from the BDF file and thus can not be spefied here
directly.

The structural/material properties must be set in the file. The code
is setup to use higher order shell elements "CQUAD".

The time integration can be controlled using the type of integrator to
use, start time, end time and number of steps to take per second.

The user can control the output to tecplot compatible f5 file by
specifying the quantities, frequency of output and filename/path.
"""
# ---------------------------------------------------------------------!
# Import necessary modules
# ---------------------------------------------------------------------!

import sys
import numpy as np

from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions

# TACS Communicator
comm = MPI.COMM_WORLD

# ---------------------------------------------------------------------!
# Parse command line arguments, if any
# ---------------------------------------------------------------------!

for arg in sys.argv:
    print(arg)

bdfFileName = "plate.bdf"  # Specify the name of the file to load which
# contains the mesh

# ---------------------------------------------------------------------!
# Configure F5 output
# ---------------------------------------------------------------------!

flag = (
    TACS.ToFH5.NODES
    | TACS.ToFH5.DISPLACEMENTS
    | TACS.ToFH5.STRAINS
    | TACS.ToFH5.STRESSES
    | TACS.ToFH5.EXTRAS
)
f5_format = "output/plate_nbg_%04d.f5"
write_freq = 1  # Set to 0 if no output is sought

# ---------------------------------------------------------------------!
# Time integration
# ---------------------------------------------------------------------!

tinit = 0.00
tfinal = 0.10
num_steps_per_sec = 100

# ---------------------------------------------------------------------!
# Properties for dynamics (Initial Conditions)
# ---------------------------------------------------------------------!

gravity = elements.GibbsVector(0.0, 0.0, -9.81)
v0 = elements.GibbsVector(0.0, 0.0, 0.0)
w0 = elements.GibbsVector(0.0, 0.0, 0.0)

# ---------------------------------------------------------------------!
# Properties for the structure
# ---------------------------------------------------------------------!

rho = 2500.0  # density, kg/m^3
E = 70.0e9  # elastic modulus, Pa
nu = 0.3  # poisson's ratio
kcorr = 5.0 / 6.0  # shear correction factor
ys = 350.0e6  # yield stress, Pa
min_thickness = 0.01  # minimum thickness of elements in m
max_thickness = 0.10  # maximum thickness of elemetns in m
thickness = 0.05  # currrent thickness of elements in m

# ---------------------------------------------------------------------!
# Load input BDF, set properties and create TACS
# ---------------------------------------------------------------------!

mesh = TACS.MeshLoader(comm)
mesh.setConvertToCoordinate(0)
mesh.scanBDFFile(bdfFileName)

num_components = mesh.getNumComponents()
for i in range(num_components):
    descriptor = mesh.getElementDescript(i)
    stiff = constitutive.isoFSDT(
        rho, E, nu, kcorr, ys, thickness, i, min_thickness, max_thickness
    )
    element = None
    if descriptor in ["CQUAD", "CQUAD9"]:
        element = elements.MITC(stiff, gravity, v0, w0)
    mesh.setElement(i, element)

tacs = mesh.createTACS(8)

######################################################################
# Time integration                                                   #
######################################################################

# Configure F5 output if tecplot output is required
f5_format = "output/plate_%04d.f5"
flag = (
    TACS.ToFH5.NODES
    | TACS.ToFH5.DISPLACEMENTS
    | TACS.ToFH5.STRAINS
    | TACS.ToFH5.STRESSES
    | TACS.ToFH5.EXTRAS
)
f5 = TACS.ToFH5(tacs, TACS.PY_SHELL, flag)

# Create the BDF integrator solver
tfinal = 0.5
num_steps_per_second = 250.0
order = 2

# Set the file output format
solver = TACS.BDFIntegrator(tacs, 0.0, tfinal, num_steps_per_second, order)

solver.setRelTol(1e-8)
# solver.setPrintLevel(2)
# solver.setOrderingType(TACS.PY_NATURAL_ORDER)
# solver.setUseLapack(1)
solver.setMaxNewtonIters(20)
solver.configureF5Output(f5, 1, f5_format)
solver.integrate()
