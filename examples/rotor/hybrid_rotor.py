#!/usr/bin/python
'''
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
'''
#---------------------------------------------------------------------!
# Import necessary modules
#---------------------------------------------------------------------!

import sys
import numpy as np

# Import argparse
import argparse

from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions

#---------------------------------------------------------------------!
# Print command line arguments, if any
#---------------------------------------------------------------------!

for arg in sys.argv:
    print arg
    
# TACS Communicator
comm = MPI.COMM_WORLD

#---------------------------------------------------------------------!
# Parse command line arguments, if any
#---------------------------------------------------------------------!

parser = argparse.ArgumentParser()

parser.add_argument('--dh', type=float, default=1.0e-12,
                    help='FD/CSD step size')

parser.add_argument('--convert_mesh', type=int, default=1,
                    help='Converts nodes to coordinate ordering for loading a gmsh BDF file')

parser.add_argument('--write_solution', type=int, default=10,
                    help='Write solution frequency')

parser.add_argument('--print_level', type=int, default=1,
                    help='Amount of printing')

parser.add_argument('--input_file', default="rotor.bdf",
                    help='Write solution frequency')

parser.add_argument('--num_steps_per_sec', type=int, default=1000,
                    help='Number of steps per second')

parser.add_argument('--tfinal', type=float, default=1.0,
                    help='Write solution frequency')

parser.add_argument('--order', type=int, default=2,
                    help='Order of integration')

args = parser.parse_args()

#---------------------------------------------------------------------!
# Set variable values based on command line arguments
#---------------------------------------------------------------------!

dh                = args.dh
convert_mesh      = args.convert_mesh
write_freq        = args.write_solution
print_level       = args.print_level
bdfFileName       = args.input_file
num_steps_per_sec = args.num_steps_per_sec
tfinal            = args.tfinal
order             = args.order

print "dh.................", dh
print "convert_mesh.......", convert_mesh
print "write_freq.........", write_freq
print "print_level........", print_level
print "bdfFileName........", bdfFileName
print "num_steps_per_sec..", num_steps_per_sec
print "tfinal.............", tfinal
print "order..............", order

#---------------------------------------------------------------------!
# Configure F5 output
#---------------------------------------------------------------------!

flag       = (TACS.ToFH5.NODES | TACS.ToFH5.DISPLACEMENTS|
              TACS.ToFH5.STRAINS| TACS.ToFH5.STRESSES|
              TACS.ToFH5.EXTRAS)

f5_format  = "output/hybrid_rotor_%04d.f5"

#---------------------------------------------------------------------!
# Properties for dynamics (Initial Conditions)
#---------------------------------------------------------------------!

gravity = elements.GibbsVector(0.0, 0.0, -9.81)
v0      = elements.GibbsVector(0.0, 0.0, 0.0)
w0      = elements.GibbsVector(0.0, 0.0, 0.0)

#---------------------------------------------------------------------!
# Properties for the structure
#---------------------------------------------------------------------!

rho           = 2500.0  # density, kg/m^3
E             = 70.0e9  # elastic modulus, Pa
nu            = 0.3     # poisson's ratio
kcorr         = 5.0/6.0 # shear correction factor
ys            = 350.0e6 # yield stress, Pa
min_thickness = 0.001    # minimum thickness of elements in m
max_thickness = 0.10    # maximum thickness of elemetns in m  
thickness     = 0.005    # currrent thickness of elements in m

#---------------------------------------------------------------------!
# Load input BDF, set properties and create TACS
#---------------------------------------------------------------------!

vars_per_node = 8

mesh = TACS.MeshLoader(comm)

mesh.setConvertToCoordinate(convert_mesh)
mesh.scanBDFFile(bdfFileName)

num_components = mesh.getNumComponents()
for i in xrange(num_components):
    print "Creating elements for component", i
    descriptor = mesh.getElementDescript(i)
    stiff      = constitutive.isoFSDT(rho, E, nu, kcorr, ys, thickness,
                                      i, min_thickness, max_thickness)
    element = None
    if descriptor in ["CQUAD", "CQUAD9"]:
        element = elements.MITC(stiff, gravity, v0, w0)
    mesh.setElement(i, element)

print "Mesh loaded into memory"
ptr, conn, xpts = mesh.getConnectivity()
num_nodes_per_elem = 9
num_elems = len(conn)/num_nodes_per_elem

print "ptr",  ptr, len(ptr)
print "conn", conn, len(conn), max(conn)
print "xpts", xpts, len(xpts)

# Print a sorted connectivity array for each element
for elem in range(0,num_elems):
    print elem, conn[num_nodes_per_elem*(elem):num_nodes_per_elem*(elem+1)]

#creator = TACS.Creator(comm, vars_per_node)
#creator.setBoundaryConditions(nodes, bcvars, ptr)
#creator.setGlobalConnectivity(num_nodes, elem_node_ptr, elem_node, elem_id_nums)
#creator->setNodes(xpts)
       
print "Creating TACS using mesh loader"
tacs = mesh.createTACS(vars_per_node)

#nelems = tacs.getNumElements()
#print "Number of elements", nelems

# Create an FH5 object for tecplot output
f5 = TACS.ToFH5(tacs, TACS.PY_SHELL, flag)

print "Time Integration"
tinit = 0.0
logFileName = bdfFileName+".log"
solver = TACS.BDFIntegrator(tacs, tinit, tfinal, num_steps_per_sec, order)
solver.setPrintLevel(print_level, logFileName)
solver.configureOutput(f5, write_freq, f5_format)
solver.integrate()
