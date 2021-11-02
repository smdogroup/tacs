"""
The nominal case is a heat conduction problem of a
1m radius plate with a Dirichilet boundary condition applied at the edges,
such that:
    T(theta) = T0 + dT * sin(2*theta)
    T0 = 70 C
    dT = 30 C
The problem is then to solve for the temperature within the boundaries of the plate.
The problem basically boils down to Laplaces problem:
    grad**2 T = 0
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
from tacs import functions, constitutive, elements, TACS, pyTACS, problems
# Radius of plate
R = 1.0
# Area of plate
area = np.pi * R ** 2

comm = MPI.COMM_WORLD

# Instantiate FEASolver
structOptions = {
    'printtimings':True,
    # Specify what type of elements we want in the f5
    'writeSolution':True,
    'outputElement': TACS.PLANE_STRESS_ELEMENT,
}

bdfFile = os.path.join(os.path.dirname(__file__), 'circ-plate-dirichlet-bcs.bdf')
FEASolver = pyTACS(bdfFile, comm, options=structOptions)

def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
    # Material properties
    rho = 2500.0        # density kg/m^3
    kappa = 230.0       # Thermal conductivity W/(m⋅K)

    # Plate geometry
    tplate = 0.005    # 5 mm

    # Setup property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, kappa=kappa)
    # Set one thickness dv for every component
    con = constitutive.PlaneStressConstitutive(prop, t=tplate, tNum=dvNum)

    # For each element type in this component,
    # pass back the appropriate tacs element object
    elemList = []
    model = elements.HeatConduction2D(con)
    for elemDescript in elemDescripts:
        if elemDescript in ['CQUAD4', 'CQUADR']:
            basis = elements.LinearQuadBasis()
        elif elemDescript in ['CTRIA3', 'CTRIAR']:
            basis = elements.LinearTriangleBasis()
        else:
            print("Uh oh, '%s' not recognized" % (elemDescript))
        elem = elements.Element2D(model, basis)
        elemList.append(elem)

    # Add scale for thickness dv
    scale = [100.0]
    return elemList, scale

# Set up constitutive objects and elements
FEASolver.createTACSAssembler(elemCallBack)

# Add Functions
FEASolver.addFunction('mass', functions.StructuralMass)
FEASolver.addFunction('ks_temp', functions.KSTemperature,
                      ksWeight=100.0)
FEASolver.addFunction('avg_temp', functions.AverageTemperature, volume=area)

# Structural problem
evalFuncs = ['mass', 'ks_temp', 'avg_temp']
sp = problems.StaticProblem(name='plate', evalFuncs=evalFuncs)

# Solve state
FEASolver(sp)

# Evaluate functions
funcs = {}
FEASolver.evalFunctions(sp, funcs)
if comm.rank == 0:
    pprint(funcs)

funcsSens = {}
FEASolver.evalFunctionsSens(sp, funcsSens)
if comm.rank == 0:
    pprint(funcsSens)

FEASolver.writeSolution(outputDir=os.path.dirname(__file__))