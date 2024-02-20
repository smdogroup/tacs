"""
Run buckling analysis of a general orthotropic plate
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
from tacs import pyTACS, constitutive, elements, utilities

comm = MPI.COMM_WORLD
dtype = utilities.BaseUI.dtype


def run_buckling_analysis(
    thickness,
    E11,
    nu12,
    E22=None,
    G12=None,
    _G23=None,
    _G13=None,
    sigma=30.0,
    num_eig=5,
    write_soln=False,
    derivatives=False,
):
    comm = MPI.COMM_WORLD

    # Instantiate FEAAssembler
    bdfFile = os.path.join(os.path.dirname(__file__), "plate.bdf")
    FEAAssembler = pyTACS(bdfFile, comm=comm)

    def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
        # Set constitutive properties
        # rho = 4540.0  # density, kg/m^3
        # # E = 70e9  # elastic modulus, Pa 118e9
        # # nu = 0.33  # poisson's ratio
        # ys = 1050e6  # yield stress, Pa
        # kappa = 6.89
        # specific_heat = 463.0

        # if E22 not provided, isotropic
        isotropic = E22 is None

        # Setup property and constitutive objects
        if isotropic:
            mat = constitutive.MaterialProperties(E=E11, nu=nu12)

            # Set one thickness dv for every component
            con = constitutive.IsoShellConstitutive(mat, t=thickness)

        else:  # orthotropic
            # assume G23, G13 = G12
            G23 = G12 if _G23 is None else _G23
            G13 = G12 if _G13 is None else _G13
            ortho_prop = constitutive.MaterialProperties(
                E1=E11, E2=E22, nu12=nu12, G12=G12, G23=G23, G13=G13
            )

            ortho_ply = constitutive.OrthotropicPly(thickness, ortho_prop)

            # one play composite constitutive model
            con = constitutive.CompositeShellConstitutive(
                [ortho_ply],
                np.array([thickness], dtype=dtype),
                np.array([0], dtype=dtype),
                tOffset=0.0,
            )

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
        cpath = os.path.dirname(__file__)
        buckling_folder = os.path.join(cpath, "buckling")
        if not os.path.exists(buckling_folder):
            os.mkdir(buckling_folder)
        bucklingProb.writeSolution(outputDir=buckling_folder)

    if comm.rank == 0:
        pprint(funcs)
        # pprint(funcsSens)

    # return the eigenvalues here
    return np.array([funcs[key] for key in funcs])


if __name__ == "__main__":
    # isotropic case
    run_buckling_analysis(
        thickness=0.07,  # 0.01
        E11=70e9,
        nu12=0.33,
        E22=None,
        G12=None,
        sigma=30.0,  # 1.0
        num_eig=12,
        write_soln=True,
    )
