"""
This problem features a circular disk, clamped at its edges and rotated about its center.
Due to centrifugal forces, plane stress will build up in the disk as a function of radial position.
Since this problem is steady in the bod-frame of the disk, we can treat it as a static problem.
This example shows how to setup TACS models featuring mixed quad/tri elements.

material = Aluminum
R = 1 m
t = 1 mm
omega = 100 rev/s
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

# Rotational velocity vector (rad/s)
omega = 2 * np.pi * np.array([0.0, 0.0, 100.0])
# Rotation center (center of disk
rotCenter = np.zeros(3)

comm = MPI.COMM_WORLD

# Instantiate FEAAssembler
structOptions = {
    "printtiming": True,
}

bdfFile = os.path.join(os.path.dirname(__file__), "disk.bdf")
FEAAssembler = pyTACS(bdfFile, comm=comm, options=structOptions)


def elemCallBack(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
    # Material properties
    rho = 2700.0  # density kg/m^3
    E = 70e9  # Young's modulus (Pa)
    nu = 0.3  # Poisson's ratio
    ys = 270.0e6  # yield stress

    # Plate geometry
    tplate = 0.001  # 1 mm
    tMin = 0.0001  # 0.1 mm
    tMax = 0.05  # 5 cm

    # Set up property model
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    # Set up constitutive model
    con = constitutive.IsoShellConstitutive(
        prop, t=tplate, tNum=dvNum, tlb=tMin, tub=tMax
    )
    transform = None
    # Set up quad/tri elements
    elemList = []
    for elemDescript in elemDescripts:
        if elemDescript == "CQUAD4":
            elem = elements.Quad4Shell(transform, con)
        elif elemDescript == "CTRIA3":
            elem = elements.Tri3Shell(transform, con)
        else:
            raise AttributeError(
                f"Element descript of type {elemDescript} not recognized."
            )
        elemList.append(elem)
    return elemList


# Set up constitutive objects and elements
FEAAssembler.initialize(elemCallBack)

# Setup static problem
staticProb = FEAAssembler.createStaticProblem(name="centrifugal")
# Add functions
staticProb.addFunction("mass", functions.StructuralMass)
staticProb.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=100.0)
# Centrifugal inertial loads about center of disk
staticProb.addCentrifugalLoad(omega, rotCenter)

# solve and evaluate functions/sensitivities
funcs = {}
funcsSens = {}
staticProb.solve()
staticProb.evalFunctions(funcs)
staticProb.evalFunctionsSens(funcsSens)
staticProb.writeSolution(outputDir=os.path.dirname(__file__))

if comm.rank == 0:
    pprint(funcs)
