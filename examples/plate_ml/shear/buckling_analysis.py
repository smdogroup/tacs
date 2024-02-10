"""
This problem performs a buckling analysis on a 4m x 3m flat plate.
The perimeter of the plate is pinned and loaded in compression (20kN/m) on its horizontal edges.
We use TACS buckling eigenvalue solver through the pyTACS BucklingProblem interface.
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
from tacs import pyTACS, constitutive, elements

comm = MPI.COMM_WORLD

# Instantiate FEAAssembler
bdfFile = os.path.join(os.path.dirname(__file__), "plate.bdf")
FEAAssembler = pyTACS(bdfFile, comm=comm)


def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
    # Set constitutive properties
    rho = 4540.0  # density, kg/m^3
    E = 70e9  # elastic modulus, Pa 118e9
    nu = 0.33  # poisson's ratio
    ys = 1050e6  # yield stress, Pa
    kappa = 6.89
    specific_heat = 463.0

    # Plate geometry
    tplate = 0.007  # 5 mm

    # Setup property and constitutive objects
    mat = constitutive.MaterialProperties(
        rho=rho, specific_heat=specific_heat, kappa=kappa, E=E, nu=nu, ys=ys
    )

    # Set one thickness dv for every component
    con = constitutive.IsoShellConstitutive(mat, t=tplate)

    # For each element type in this component,
    # pass back the appropriate tacs element object
    elemList = []
    for descript in elemDescripts:
        transform = None
        if descript in ["CQUAD4", "CQUADR"]:
            elem = elements.Quad4Shell(transform, con)
        elif descript in ["CQUAD9", "CQUAD"]:
            elem = elements.Quad9Shell(transform, con)
        else:
            raise AssertionError("Non CQUAD4 Elements in this plate?")

        elemList.append(elem)

    # Add scale for thickness dv
    scale = [100.0]
    return elemList, scale


# Set up constitutive objects and elements
FEAAssembler.initialize(elemCallBack)

# Setup buckling problem
bucklingProb = FEAAssembler.createBucklingProblem(name="buckle", sigma=10.0, numEigs=5)
bucklingProb.setOption("printLevel", 2)
# Add Loads
# bucklingProb.addLoadFromBDF(loadID=1)

# exit()

# solve and evaluate functions/sensitivities
funcs = {}
funcsSens = {}
bucklingProb.solve()
bucklingProb.evalFunctions(funcs)
bucklingProb.evalFunctionsSens(funcsSens)
bucklingProb.writeSolution(outputDir=os.path.dirname(__file__))


if comm.rank == 0:
    pprint(funcs)
    pprint(funcsSens)
