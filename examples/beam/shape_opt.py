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

The optimization is setup using TACS' MACH module, which acts as a wrapper
for MDOLab's MACH library.
"""

# [docs:imports-start]
import numpy as np
import os

from pygeo import DVGeometry
from pyoptsparse import Optimization, OPT

from tacs.mach import StructProblem
from tacs import pyTACS
from tacs import elements, constitutive, functions

# [docs:imports-end]

# [docs:parameters-start]
bdf_file = os.path.join(os.path.dirname(__file__), "Slender_Beam.bdf")
ffd_file = os.path.join(os.path.dirname(__file__), "ffd_8_linear.fmt")

# Beam thickness
t = 0.01  # m
# Length of beam
L = 1.0  # m

# Material properties
rho = 2780.0  # kg /m^3
E = 70.0e9  # Pa
nu = 0.0
ys = 420.0e6

# Shear force applied at tip
V = 2.5e4  # N
# [docs:parameters-end]


# Callback function used to setup TACS element objects and DVs
# [docs:element-callback-start]
def element_callback(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
    # Setup (isotropic) property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    con = constitutive.IsoShellConstitutive(prop, t=t, tNum=-1)
    refAxis = np.array([1.0, 0.0, 0.0])
    transform = elements.ShellRefAxisTransform(refAxis)
    elem = elements.Quad4Shell(transform, con)
    return elem


# [docs:element-callback-end]


# [docs:pytacs-init-start]
FEAAssembler = pyTACS(bdf_file)
FEAAssembler.initialize(element_callback)
# [docs:pytacs-init-end]

# [docs:dvgeo-setup-start]
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
# [docs:dvgeo-setup-end]

# [docs:static-problem-start]
staticProb = FEAAssembler.createStaticProblem("tip_shear")
# Add TACS Functions
staticProb.addFunction("mass", functions.StructuralMass)
staticProb.addFunction("ks_vmfailure", functions.KSFailure, safetyFactor=1.0, ksWeight=100.0)
# Add forces to static problem
staticProb.addLoadToNodes(1112, [0.0, V, 0.0, 0.0, 0.0, 0.0], nastranOrdering=True)
# [docs:static-problem-end]

# [docs:struct-problem-start]
structProb = StructProblem(staticProb, FEAAssembler, DVGeo=DVGeo)
# [docs:struct-problem-end]


# [docs:struct-obj-start]
def structObj(x):
    """Evaluate the objective and constraints"""
    funcs = {}
    structProb.setDesignVars(x)
    DVGeo.setDesignVars(x)
    structProb.solve()
    structProb.evalFunctions(funcs)
    structProb.writeSolution()
    if structProb.comm.rank == 0:
        print(x)
        print(funcs)

    return funcs, False


# [docs:struct-obj-end]


# [docs:struct-sens-start]
def structSens(x, funcs):
    """Evaluate the objective and constraint sensitivities"""
    funcsSens = {}
    structProb.evalFunctionsSens(funcsSens)
    for func in funcsSens:
        funcsSens[func].pop("struct")
    return funcsSens, False


# [docs:struct-sens-end]


# [docs:opt-setup-start]
# Now we create the structural optimization problem:
optProb = Optimization("Mass min", structObj)
optProb.addObj("tip_shear_mass")
structProb.addVariablesPyOpt(optProb)
DVGeo.addVariablesPyOpt(optProb)
optProb.addCon("tip_shear_ks_vmfailure", upper=1.0)

optProb.printSparsity()

opt = OPT(
    "SLSQP",
    options={
        "MAXIT": 100,
        "IPRINT": 1,
        "IFILE": os.path.join(os.path.dirname(__file__), "SLSQP.out"),
    },
)
# [docs:opt-setup-end]

# [docs:run-opt-start]
# Finally run the actual optimization
sol = opt(optProb, sens=structSens, storeSens=False)
# [docs:run-opt-end]
