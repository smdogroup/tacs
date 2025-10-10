"""
This example demonstrates TACS structural optimization capabilities.
The beam model that we will be using for this problem is a rectangular beam,
cantilevered, with a shear load applied at the tip. The beam is discretized using
1001 shell elements along it's span and depth.

The optimization problem is as follows:
Minimize the mass of the beam with respect to the depth of the cross-section along the span,
subject to a max stress constraint dictated by the materials yield stress.

In order to change the shape of the FEM we use a free-form deformation (FFD) volume
parmeterization scheme provided by the pyGeo library.

An aproximate analytical solution can be derived from beam theory,
by realizing that the stress at any spanwise cross-section in the beam
can be found independently using:
    sigma(x,y) = y*M(x)/I
An analytical solution for this problem can be shown to be:
    d(x) = sqrt(6*V*(L-x)/(t*sigma_y))

The optimization is setup using TACS' MPHYS module, which acts as a wrapper
for OpenMDAO.
"""

import numpy as np
import os

from pygeo import DVGeometry
from pyoptsparse import Optimization, SNOPT

from tacs.mach import StructProblem
from tacs import pyTACS
from tacs import elements, constitutive, functions

bdf_file = os.path.join(os.path.dirname(__file__), "Slender_Beam.bdf")
ffd_file = os.path.join(os.path.dirname(__file__), "ffd_8_linear.fmt")

# Beam thickness
t = 0.01  # m
# Length of beam
L = 1.0

# Material properties
rho = 2780.0  # kg /m^3
E = 70.0e9
nu = 0.0
ys = 420.0e6

# Shear force applied at tip
V = 2.5e4


# Callback function used to setup TACS element objects and DVs
def element_callback(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
    # Setup (isotropic) property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    con = constitutive.IsoShellConstitutive(prop, t=t, tNum=-1)
    # TACS shells are sometimes a little overly-rigid in shear
    # We can reduce this effect by decreasing the drilling regularization
    con.setDrillingRegularization(0.1)
    refAxis = np.array([1.0, 0.0, 0.0])
    transform = elements.ShellRefAxisTransform(refAxis)
    elem = elements.Quad4Shell(transform, con)
    return elem


FEAAssembler = pyTACS(bdf_file)
FEAAssembler.initialize(element_callback)

DVGeo = DVGeometry(fileName=ffd_file)
# Create reference axis
nRefAxPts = DVGeo.addRefAxis(name="centerline", alignIndex="i", yFraction=0.5)


# Set up global design variables
def depth(val, geo):
    for i in range(nRefAxPts):
        geo.scale_y["centerline"].coef[i] = val[i]


DVGeo.addGlobalDV(
    dvName="depth",
    value=np.ones(nRefAxPts),
    func=depth,
    lower=1e-3,
    upper=10.0,
    scale=20.0,
)

staticProb = FEAAssembler.createStaticProblem("tip_shear")
# Add TACS Functions
staticProb.addFunction("mass", functions.StructuralMass)
staticProb.addFunction(
    "ks_vmfailure", functions.KSFailure, safetyFactor=1.0, ksWeight=100.0
)
# Add forces to static problem
staticProb.addLoadToNodes(1112, [0.0, V, 0.0, 0.0, 0.0, 0.0], nastranOrdering=True)

structProb = StructProblem(staticProb, FEAAssembler, DVGeo=DVGeo)


def structObj(x):
    """Evaluate the objective and constraints"""
    funcs = {}
    structProb.setDesignVars(x)
    DVGeo.setDesignVars(x)
    structProb.solve()
    structProb.evalFunctions(funcs)
    # structProb.evalConstraints(funcs)
    structProb.writeSolution()
    if structProb.comm.rank == 0:
        print(x)
        print(funcs)

    return funcs, False


def structSens(x, funcs):
    """Evaluate the objective and constraint sensitivities"""
    funcsSens = {}
    structProb.evalFunctionsSens(funcsSens)
    # structProb.evalConstraintsSens(funcsSens)
    for func in funcsSens:
        funcsSens[func].pop("struct")
    return funcsSens, False


# Now we create the structural optimization problem:
optProb = Optimization("Mass min", structObj)
optProb.addObj("tip_shear_mass")
structProb.addVariablesPyOpt(optProb)
# structProb.addConstraintsPyOpt(optProb)
DVGeo.addVariablesPyOpt(optProb)
optProb.addCon("tip_shear_ks_vmfailure", upper=1.0)

optProb.printSparsity()

optOptions = {
    "Major feasibility tolerance": 1e-4,
    "Major optimality tolerance": 1e-4,
    "Major iterations limit": 200,
    "Minor iterations limit": 150000,
    "Iterations limit": 1000000,
    "Major step limit": 0.1,
    "Function precision": 1.0e-8,
    "Problem Type": "Minimize",
    "New superbasics limit": 500,
    "Penalty parameter": 1e3,
}

opt = SNOPT(options=optOptions)

# Finally run the actual optimization
sol = opt(optProb, sens=structSens, storeSens=False)
