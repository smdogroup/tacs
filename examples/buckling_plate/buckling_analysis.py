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


def run_buckling_analysis(
    thickness=0.007,
    E=70e9,
    nu=0.33,
    sigma=30.0,
    num_eig=5,
    write_soln=False,
    derivatives=False,
):
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
        tplate = thickness  # 5 mm

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
    bucklingProb = FEAAssembler.createBucklingProblem(
        name="buckle", sigma=sigma, numEigs=num_eig
    )
    bucklingProb.setOption("printLevel", 2)

    # exit()

    # solve and evaluate functions/sensitivities
    funcs = {}
    funcsSens = {}
    bucklingProb.solve()
    bucklingProb.evalFunctions(funcs)
    if derivatives:
        bucklingProb.evalFunctionsSens(funcsSens)
    if write_soln:
        bucklingProb.writeSolution(outputDir=os.path.dirname(__file__))

    if comm.rank == 0:
        pprint(funcs)
        # pprint(funcsSens)

    # return the eigenvalues here
    funcs_list = [funcs[key] for key in funcs]
    return funcs_list


if __name__ == "__main__":
    run_buckling_analysis(
        thickness=0.07,  # 0.01
        E=70e9,
        nu=0.33,
        sigma=30.0,  # 1.0
        num_eig=12,
        write_soln=True,
    )
