"""
This problem will show how to use some of pytacs more advanced load setting procedures.
The nominal case is a 1m x 1m flat plate. The perimeter of the plate is fixed in
all 6 degrees of freedom. The plate comprises 900 CQUAD4 elements.
We consider three structural problems:
    1. A static case where a 10 kN point force is applied at the plate center
    2. A transient problem with a pressure load that varies in time and space
    given by:
        P(x,y,t) = Pmax * sin(2*pi*x/L) * sin(2*pi*y/L) * sin(2*pi*fhz*t)
        where:
            Pmax = 100 kPa
            fhz = 1.0 hz
            L = 1.0 m
    3. A modal problem where the dynamic eigenmodes are solved for the plate structure
"""
# ==============================================================================
# Standard Python modules
# ==============================================================================
from __future__ import print_function
import os

# ==============================================================================
# External Python modules
# ==============================================================================
from pprint import pprint

import numpy as np
from mpi4py import MPI

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs import functions, constitutive, elements, pyTACS

# Center load magnitude
Q = 1e4
# Pressure magnitude
Pmax = 100e3
# Pressure frequency
fhz = 1.0
# Length of square plate
L = 1.0

comm = MPI.COMM_WORLD

# Instantiate FEAAssembler
structOptions = {
    "printtiming": True,
}

bdfFile = os.path.join(os.path.dirname(__file__), "plate.bdf")
FEAAssembler = pyTACS(bdfFile, comm=comm, options=structOptions)


def elemCallBack(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
    # Material properties
    rho = 2500.0  # density kg/m^3
    E = 70e9  # Young's modulus (Pa)
    nu = 0.3  # Poisson's ratio
    ys = 464.0e6  # yield stress

    # Plate geometry
    tplate = 0.005  # 1 mm
    tMin = 0.0001  # 0.1 mm
    tMax = 0.05  # 5 cm

    # Set up property model
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    # Set up constitutive model
    con = constitutive.IsoShellConstitutive(
        prop, t=tplate, tNum=dvNum, tlb=tMin, tub=tMax
    )
    transform = None
    # Set up element
    elem = elements.Quad4Shell(transform, con)
    scale = [100.0]
    return elem, scale


# Set up constitutive objects and elements
FEAAssembler.initialize(elemCallBack)

# List to hold all problems
allProblems = []

# Setup static problem
staticProb = FEAAssembler.createStaticProblem(name="point_force")
# Add functions
staticProb.addFunction("mass", functions.StructuralMass)
staticProb.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=100.0)
# Add point force to node 481 (center of plate)
F = np.array([0.0, 0.0, Q, 0.0, 0.0, 0.0])
staticProb.addLoadToNodes(481, F, nastranOrdering=True)
# Add static problem to list
allProblems.append(staticProb)

# Setup transient problem
# turn on print for solver in options
transientOptions = {"printlevel": 1}
transientProb = FEAAssembler.createTransientProblem(
    name="pressure", tInit=0.0, tFinal=10.0, numSteps=50, options=transientOptions
)
# Add functions
transientProb.addFunction("mass", functions.StructuralMass)
transientProb.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=100.0)
# Add presure load over plate
# pynastran bdf object for parsing mesh info
bdfInfo = FEAAssembler.getBDFInfo()
# cross-reference bdf object to use some of pynastrans advanced features
bdfInfo.cross_reference()
# Loop through each element and get the spatial variation of the pressure
Pxy = []
eIDs = []
for eID in bdfInfo.elements:
    [x, y, z] = bdfInfo.elements[eID].Centroid()
    Pmag = Pmax * np.sin(2 * np.pi * x / L) * np.sin(2 * np.pi * y / L)
    Pxy.append(Pmag)
    eIDs.append(eID)
Pxy = np.array(Pxy)
# Loop through each time step and add time varying load
timeSteps = transientProb.getTimeSteps()
for step_i, time in enumerate(timeSteps):
    # Multiply by time factor
    Pxyt = Pxy * np.sin(2 * np.pi * fhz * time)
    # Add pressure to problem for timestep
    transientProb.addPressureToElements(step_i, eIDs, Pxyt, nastranOrdering=True)
# Add transient problem to list
allProblems.append(transientProb)

# Setup modal problem (Note: eigenvalues are in (rad/s)**2
modalProb = FEAAssembler.createModalProblem("dynamic_modes", sigma=1e4, numEigs=10)
# Add modal problem to list
allProblems.append(modalProb)

# Solve all problems and evaluate functions
funcs = {}
funcsSens = {}
for problem in allProblems:
    problem.solve()
    problem.evalFunctions(funcs)
    problem.evalFunctionsSens(funcsSens)
    problem.writeSolution(outputDir=os.path.dirname(__file__))

if comm.rank == 0:
    pprint(funcs)
