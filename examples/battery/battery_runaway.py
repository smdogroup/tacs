from __future__ import print_function
import os
from pprint import pprint

import numpy as np
from mpi4py import MPI

from tacs import functions, constitutive, elements, TACS, pyTACS

"""
This example demonstrates transient heating of a battery pack with one cell undergoing thermal runaway.
The domain is an aluminum battery pack with 9 cylindrical cells embedded in a grid-pattern. The cell in
the corner undergoes thermal runaway, releasing a large amount of heat for 2 seconds. The heat conduction
of the battery pack is then computed for 5 seconds in total. The maximum temperature within a cell is computed
over all time steps. This is done for 3 cells: the cell undergoing thermal runaway, the nearest cell adject
to it, and the cell on the diagonal near it. Computing the maximum temperature of these cells could be used to
prevent other cells in the pack from going into thermal runaway, leading to a cascading failure.

This example demonstrates a number of useful pyTACS features, including:
- Transient heat conduction physics
- Two materials modeled in the same mesh
- Time-specified loading: the right-hand-side is given as a user-defined function of time
- Evaluating multiple functions in different regions
- Easy domain-selection enabling the previous three items
"""

comm = MPI.COMM_WORLD

# Name of the bdf file to get the mesh
bdfFile = os.path.join(os.path.dirname(__file__), "battery_pack.bdf")

# Instantiate the pyTACS object
FEAAssembler = pyTACS(bdfFile, comm)

# Specify the plate thickness
tplate = 0.065

# Define material properties for two materials used in this problem
# Properties of the battery cells
battery_rho = 1460.0  # density kg/m^3
battery_kappa = 1.3  # Thermal conductivity W/(m⋅K)
battery_cp = 880.0  # Specific heat J/(kg⋅K)

# Properties of the battery pack (aluminum)
alum_rho = 2700.0  # density kg/m^3
alum_kappa = 204.0  # Thermal conductivity W/(m⋅K)
alum_cp = 883.0  # Specific heat J/(kg⋅K)

# The callback function to define the element properties
def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):

    # Setup property and constitutive objects
    if (
        compDescript == "Block"
    ):  # If the bdf file labels this component as "Block", then it is aluminum
        prop = constitutive.MaterialProperties(
            rho=alum_rho, kappa=alum_kappa, specific_heat=alum_cp
        )
    else:  # otherwise it is a battery
        prop = constitutive.MaterialProperties(
            rho=battery_rho, kappa=battery_kappa, specific_heat=battery_cp
        )

    # Set one thickness value for every component
    con = constitutive.PlaneStressConstitutive(prop, t=tplate, tNum=-1)

    # For each element type in this component,
    # pass back the appropriate tacs element object
    elemList = []
    model = elements.HeatConduction2D(con)
    for elemDescript in elemDescripts:
        if elemDescript in ["CQUAD4", "CQUADR"]:
            basis = elements.LinearQuadBasis()
        elif elemDescript in ["CTRIA3", "CTRIAR"]:
            basis = elements.LinearTriangleBasis()
        else:
            print("Element '%s' not recognized" % (elemDescript))
        elem = elements.Element2D(model, basis)
        elemList.append(elem)

    return elemList


# Set up constitutive objects and elements
FEAAssembler.initialize(elemCallBack)

# Create a transient problem that will represent time-varying heat conduction
transientProblem = FEAAssembler.createTransientProblem(
    "Transient", tInit=0.0, tFinal=5.0, numSteps=50, options={"printLevel": 1}
)

# Get the time steps and define the loads
timeSteps = transientProblem.getTimeSteps()
for i, t in enumerate(timeSteps):
    if t <= 2.0:  # only apply the load for the first 2 seconds (step function)
        # select the component of the battery undergoing thermal runaway
        compIDs = FEAAssembler.selectCompIDs(include=["Battery.00"])

        # Define the heat-flux: apply 6000 Watts spread out over the face of the cell undergoing thermal runaway
        transientProblem.addLoadToComponents(i, compIDs, [6000.0])

# Define the functions of interest as maximum temperature withing 3 different batteries
compIDs_00 = FEAAssembler.selectCompIDs(
    ["Battery.00"]
)  # battery undergoing thermal runaway
compIDs_01 = FEAAssembler.selectCompIDs(["Battery.01"])  # adjecent battery
compIDs_04 = FEAAssembler.selectCompIDs(["Battery.04"])  # diagonal battery

transientProblem.addFunction(
    "ks_temp_corner", functions.KSTemperature, ksWeight=100.0, compIDs=compIDs_00
)
transientProblem.addFunction(
    "ks_temp_adjacent", functions.KSTemperature, ksWeight=100.0, compIDs=compIDs_01
)
transientProblem.addFunction(
    "ks_temp_diagonal", functions.KSTemperature, ksWeight=100.0, compIDs=compIDs_04
)

# Solve state for each problem, evaluate functions and sensitivities
funcs = {}  # Store the function values
funcsSens = {}  # Store the function gradients
transientProblem.solve()
transientProblem.evalFunctions(funcs)
transientProblem.writeSolution()
transientProblem.evalFunctionsSens(funcsSens)

if comm.rank == 0:
    pprint(funcs)
    pprint(funcsSens)
