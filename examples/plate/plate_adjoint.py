#!/usr/bin/python
'''
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

'''
#---------------------------------------------------------------------!
# Import necessary modules
#---------------------------------------------------------------------!

import sys
import numpy as np

from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions

# TACS Communicator
comm = MPI.COMM_WORLD

#---------------------------------------------------------------------!
# Parse command line arguments, if any
#---------------------------------------------------------------------!

for arg in sys.argv:
    print arg

bdfFileName = "plate.bdf" # Specify the name of the file to load which
                          # contains the mesh
                          
#---------------------------------------------------------------------!
# Set the design variables
#---------------------------------------------------------------------!

x = np.array([0.03])

#---------------------------------------------------------------------!
# Configure F5 output
#---------------------------------------------------------------------!

#flag       = (TACS.ToFH5.NODES | TACS.ToFH5.DISPLACEMENTS|
#              TACS.ToFH5.STRAINS| TACS.ToFH5.STRESSES|
#              TACS.ToFH5.EXTRAS)
#f5_format  = "output/plate_nbg_%04d.f5"
#write_freq = 1  # Set to 0 if no output is sought

#---------------------------------------------------------------------!
# Time integration
#---------------------------------------------------------------------!

tinit             = 0.00
tfinal            = 0.01
num_steps_per_sec = 1000

#---------------------------------------------------------------------!
# Properties for dynamics (Initial Conditions)
#---------------------------------------------------------------------!

gravity = elements.GibbsVector(np.array([0.0, 0.0, -9.81]))
v0      = elements.GibbsVector(np.array([0.0, 0.0, 0.0]))
w0      = elements.GibbsVector(np.array([0., 0., 0.]))

#---------------------------------------------------------------------!
# Properties for the structure
#---------------------------------------------------------------------!

rho           = 2500.0  # density, kg/m^3
E             = 70.0e9  # elastic modulus, Pa
nu            = 0.3     # poisson's ratio
kcorr         = 5.0/6.0 # shear correction factor
ys            = 350.0e6 # yield stress, Pa
min_thickness = 0.01    # minimum thickness of elements in m
max_thickness = 0.10    # maximum thickness of elemetns in m  
thickness     = 0.05    # currrent thickness of elements in m

#---------------------------------------------------------------------!
# Load input BDF, set properties and create TACS
#---------------------------------------------------------------------!

mesh            = TACS.MeshLoader(comm)
mesh.scanBDFFile(bdfFileName)

num_components  = mesh.getNumComponents()
for i in xrange(num_components):
    descriptor  = mesh.getElementDescript(i)
    stiff       = constitutive.isoFSDT(rho, E, nu, kcorr, ys, thickness,
                                      i, min_thickness, max_thickness)
    element     = None
    if descriptor in ["CQUAD"]:
        element = elements.MITC(stiff, gravity, v0, w0)        
        mesh.setElement(i, element)

tacs            = mesh.createTACS(8)

#---------------------------------------------------------------------!
# Create the function list for adjoint solve
#---------------------------------------------------------------------!

funcs = []
funcs.append(functions.ksfailure(tacs, 100.0))
funcs.append(functions.compliance(tacs))

#---------------------------------------------------------------------#
# Setup space for function values and their gradients
#---------------------------------------------------------------------#

num_funcs       = len(funcs)
num_design_vars = len(x)

fvals           = np.zeros(num_funcs)
dfdx            = np.zeros(num_funcs*num_design_vars)

fvals_fd        = np.zeros(num_funcs)
dfdx_fd         = np.zeros(num_funcs*num_design_vars)

#---------------------------------------------------------------------#
# Test all BDF Integrators for these functions
#---------------------------------------------------------------------#
# BDF Integrator
for bdf_order in [1,2,3]:
    bdf = TACS.BDFIntegrator(tacs, tinit, tfinal, num_steps_per_sec, bdf_order)
    bdf.setPrintLevel(0)
    bdf.setJacAssemblyFreq(1)
    bdf.setFunction(funcs)

    print("")
    print("Evaluating adjoint gradient for BDF method of order ", bdf_order)
    bdf.getFuncGrad(num_design_vars, x, fvals, dfdx)
    print(fvals, dfdx)

    print("Evaluating finite difference gradient BDF method of order ", bdf_order)
    bdf.getFDFuncGrad(num_design_vars, x, fvals_fd, dfdx_fd, 1.0e-6)
    print(fvals_fd, dfdx_fd)

    print("Error for BDF order ", bdf_order)
    print(fvals-fvals_fd, dfdx-dfdx_fd)

# DIRK Integrator
for num_stages in [1,2,3]:
    dirk = TACS.DIRKIntegrator(tacs, tinit, tfinal, num_steps_per_sec, num_stages)
    dirk.setPrintLevel(0)
    dirk.setJacAssemblyFreq(1)
    dirk.setFunction(funcs)

    print("")
    print("Evaluating adjoint gradient for DIRK method of order ", num_stages+1)
    dirk.getFuncGrad(num_design_vars, x, fvals, dfdx)
    print(fvals, dfdx)

    print("Evaluating finite difference gradient DIRK method of order ", num_stages+1)
    dirk.getFDFuncGrad(num_design_vars, x, fvals_fd, dfdx_fd, 1.0e-6)
    print(fvals_fd, dfdx_fd)

    print("Error for DIRK order ", num_stages+1)
    print(fvals-fvals_fd, dfdx-dfdx_fd)

# ABM Integrator
for abm_order in [1,2,3,4,5,6]:
    abm = TACS.ABMIntegrator(tacs, tinit, tfinal, num_steps_per_sec, abm_order)
    abm.setPrintLevel(0)
    abm.setJacAssemblyFreq(1)
    abm.setFunction(funcs)

    print("")
    print("Evaluating adjoint gradient for ABM method of order ", abm_order)
    abm.getFuncGrad(num_design_vars, x, fvals, dfdx)
    print(fvals, dfdx)

    print("Evaluating finite difference gradient ABM method of order ", abm_order)
    abm.getFDFuncGrad(num_design_vars, x, fvals_fd, dfdx_fd, 1.0e-6)
    print(fvals_fd, dfdx_fd)

    print("Error for ABM order ", abm_order)
    print(fvals-fvals_fd, dfdx-dfdx_fd)

