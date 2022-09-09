#!/usr/bin/python
"""
This python script is used to test the adjoint implementation on plate
example.

The initial conditions such as the velocity, angular velocity and
accerations can be specified in the file as GibbsVectors. The boundary
conditions come from the BDF file and thus can not be spefied here
directly.

The structural/material properties must be set in the file. The code
is setup to use higher order shell elements "CQUAD".

The time integration can be controlled using the type of integrator to
use, start time, end time and number of steps to take per second.

A list of functions of interest is sent as input to the adjoint
solver. One can create new functions and append it to the list.

"""
# ---------------------------------------------------------------------!
# Import necessary modules
# ---------------------------------------------------------------------!

import sys
import numpy as np

# Import argparse
import argparse

from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions

# TACS Communicator
comm = MPI.COMM_WORLD

# ---------------------------------------------------------------------!
# Parse command line arguments, if any
# ---------------------------------------------------------------------!
# Parse the arguments
parser = argparse.ArgumentParser()
parser.add_argument("--dh", type=float, default=1.0e-8, help="FD/CSD step size")
parser.add_argument(
    "--convert_mesh", type=int, default=0, help="Converts nodes to coordinate ordering"
)
args = parser.parse_args()

# ---------------------------------------------------------------------!
# Set variable values based on command line arguments
# ---------------------------------------------------------------------!

dh = args.dh
convert_mesh = args.convert_mesh
bdfFileName = "plate.bdf"  # Specify the name of the file to load which
# contains the mesh

# ---------------------------------------------------------------------!
# Set the design variables
# ---------------------------------------------------------------------!

x = np.array([0.01], dtype=np.complex)

# ---------------------------------------------------------------------!
# Configure F5 output
# ---------------------------------------------------------------------!

# flag       = (TACS.ToFH5.NODES | TACS.ToFH5.DISPLACEMENTS|
#              TACS.ToFH5.STRAINS| TACS.ToFH5.STRESSES|
#              TACS.ToFH5.EXTRAS)
# f5_format  = "output/plate_nbg_%04d.f5"
# write_freq = 1  # Set to 0 if no output is sought

# ---------------------------------------------------------------------!
# Time integration
# ---------------------------------------------------------------------!

tinit = 0.00
tfinal = 0.1
num_steps_per_sec = 1000

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
nu = 0.30  # poisson's ratio
kcorr = 5.0 / 6.0  # shear correction factor
ys = 350.0e6  # yield stress, Pa
min_thickness = 0.01  # minimum thickness of elements in m
max_thickness = 0.10  # maximum thickness of elemetns in m
thickness = 0.05  # currrent thickness of elements in m

# ---------------------------------------------------------------------!
# Load input BDF, set properties and create TACS
# ---------------------------------------------------------------------!

mesh = TACS.MeshLoader(comm)
mesh.setConvertToCoordinate(convert_mesh)

mesh.scanBDFFile(bdfFileName)

num_components = mesh.getNumComponents()
for i in range(num_components):
    descriptor = mesh.getElementDescript(i)
    stiff = constitutive.isoFSDT(
        rho, E, nu, kcorr, ys, thickness, i, min_thickness, max_thickness
    )
    element = None
    if descriptor in ["CQUAD"]:
        element = elements.MITC(stiff, gravity, v0, w0)
        mesh.setElement(i, element)

tacs = mesh.createTACS(8)

# ---------------------------------------------------------------------!
# Create the function list for adjoint solve
# ---------------------------------------------------------------------!

funcs = []
funcs.append(functions.StructuralMass(tacs))
funcs.append(functions.Compliance(tacs))
funcs.append(functions.InducedFailure(tacs, 20.0))
funcs.append(functions.KSFailure(tacs, 100.0))


# ---------------------------------------------------------------------#
# Setup space for function values and their gradients
# ---------------------------------------------------------------------#

num_funcs = len(funcs)
num_design_vars = len(x)

fvals = np.zeros(num_funcs)
dfdx = np.zeros(num_funcs * num_design_vars)

fvals_fd = np.zeros(num_funcs)
dfdx_fd = np.zeros(num_funcs * num_design_vars)

# ---------------------------------------------------------------------#
# Function to print errors and optionally store data for more post
# processin
# ---------------------------------------------------------------------#

data = []


def print_details(method, function, index, stepsize, fvals, adj_dfdx, fd_dfdx):
    # data.append([method, function, stepsize, adj_dfdx, fd_dfdx])
    record = dict(
        method=method,
        function=function,
        index=index,
        dh=stepsize,
        adjoint=adj_dfdx,
        complex_step=fd_dfdx,
        error=fd_dfdx - adj_dfdx,
    )
    data.append(record)
    print(
        "%10s %20s %4d %25.16e %25.16e %25.16e %25.16e %25.16e"
        % (
            method,
            function,
            index,
            stepsize,
            fvals,
            adj_dfdx,
            fd_dfdx,
            fd_dfdx - adj_dfdx,
        )
    )


# ---------------------------------------------------------------------#
# BDF Integrator
# ---------------------------------------------------------------------#

for bdf_order in [1, 2, 3]:
    bdf = TACS.BDFIntegrator(tacs, tinit, tfinal, num_steps_per_sec, bdf_order)
    bdf.setPrintLevel(0)
    bdf.setJacAssemblyFreq(1)
    bdf.setFunction(funcs)
    bdf.getFuncGrad(num_design_vars, x, fvals, dfdx)
    bdf.getFDFuncGrad(num_design_vars, x, fvals_fd, dfdx_fd, dh)

    fnum = 0
    for func in funcs:
        print_details(
            "BDF" + str(bdf_order),
            func.__class__.__name__,
            fnum,
            dh,
            fvals[fnum],
            np.real(dfdx[fnum]),
            np.real(dfdx_fd[fnum]),
        )
        fnum += 1


# ---------------------------------------------------------------------#
# DIRK Integrator
# ---------------------------------------------------------------------#

for order in [2, 3, 4]:
    dirk = TACS.DIRKIntegrator(tacs, tinit, tfinal, num_steps_per_sec, order)
    dirk.setPrintLevel(0)
    dirk.setJacAssemblyFreq(1)
    dirk.setFunction(funcs)
    dirk.getFuncGrad(num_design_vars, x, fvals, dfdx)
    dirk.getFDFuncGrad(num_design_vars, x, fvals_fd, dfdx_fd, dh)
    fnum = 0
    for func in funcs:
        print_details(
            "DIRK" + str(order),
            func.__class__.__name__,
            fnum,
            dh,
            fvals[fnum],
            np.real(dfdx[fnum]),
            np.real(dfdx_fd[fnum]),
        )
        fnum += 1

# ---------------------------------------------------------------------#
# ABM Integrator
# ---------------------------------------------------------------------#

for abm_order in [1, 2, 3, 4, 5, 6]:
    abm = TACS.ABMIntegrator(tacs, tinit, tfinal, num_steps_per_sec, abm_order)
    abm.setPrintLevel(0)
    abm.setJacAssemblyFreq(1)
    abm.setFunction(funcs)
    abm.getFuncGrad(num_design_vars, x, fvals, dfdx)
    abm.getFDFuncGrad(num_design_vars, x, fvals_fd, dfdx_fd, dh)
    fnum = 0
    for func in funcs:
        print_details(
            "ABM" + str(abm_order),
            func.__class__.__name__,
            fnum,
            dh,
            fvals[fnum],
            np.real(dfdx[fnum]),
            np.real(dfdx_fd[fnum]),
        )
        fnum += 1

# ---------------------------------------------------------------------#
# NBG Integrator
# ---------------------------------------------------------------------#

nbg = TACS.NBGIntegrator(tacs, tinit, tfinal, num_steps_per_sec, 2)
nbg.setPrintLevel(0)
nbg.setJacAssemblyFreq(1)
nbg.setFunction(funcs)
nbg.getFuncGrad(num_design_vars, x, fvals, dfdx)
nbg.getFDFuncGrad(num_design_vars, x, fvals_fd, dfdx_fd, dh)
fnum = 0
for func in funcs:
    print_details(
        "NBG2",
        func.__class__.__name__,
        fnum,
        dh,
        fvals[fnum],
        np.real(dfdx[fnum]),
        np.real(dfdx_fd[fnum]),
    )
    fnum += 1

# print data
