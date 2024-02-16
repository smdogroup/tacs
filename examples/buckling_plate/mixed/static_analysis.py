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


def run_static_analysis(
    thickness=0.01, E=70e9, nu=0.33, write_soln=False
):
    comm = MPI.COMM_WORLD

    # Instantiate FEAAssembler
    bdfFile = os.path.join(os.path.dirname(__file__), "plate.bdf")
    FEAAssembler = pyTACS(bdfFile, comm=comm)

    def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
        # Set constitutive properties
        rho = 4540.0  # density, kg/m^3
        # E = 70e9  # elastic modulus, Pa 118e9
        # nu = 0.33  # poisson's ratio
        ys = 1050e6  # yield stress, Pa
        kappa = 6.89
        specific_heat = 463.0

        # Plate geometry
        tplate = thickness
        # tplate = 0.007  # 5 mm

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

    # debug the static problem first
    SP = FEAAssembler.createStaticProblem(name="static")
    SP.solve()
    if write_soln:
        SP.writeSolution(outputDir=os.path.dirname(__file__))

    # test the average stresses routine
    avgStresses = FEAAssembler.assembler.getAverageStresses()
    print(f"avg Stresses = {avgStresses}")
    return avgStresses


if __name__ == "__main__":
    run_static_analysis(
        thickness=0.01, E=70e9, nu=0.33, write_soln=True
    )
