"""
The nominal case is a 1m x 1m flat plate under a 1 kN pressure load. The
perimeter of the plate is fixed in all 6 degrees of freedom.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import os
import pathlib

# ==============================================================================
# External Python modules
# ==============================================================================
from mpi4py import MPI
import numpy as np

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs import pyTACS, constitutive, elements, functions


def setup(compLPModel="LamParamFull", useVMFailure=False, useModifiedTsaiWu=False):
    """Setup FEASolver and static problem.

    Parameters
    ----------
    useComp : bool, optional
        Whether to use composite material, by default False
    useVMFailure : bool, optional
        Whether to use von Mises failure criterion, by default False
    useModifiedTsaiWu : bool, optional
        Whether to use modified Tsai-Wu failure criterion, by default False

    Returns
    -------
    FEASolver : FEASolver
        The TACS solver object
    sp : StaticProblem
        The static problem object
    """

    comm = MPI.COMM_WORLD

    # Instantiate FEASolver
    structOptions = {
        "printTiming": True,
    }
    staticProblemOptions = {
        "writeSolution": True,
    }

    # fName = pathlib.Path(fName).stem
    bdfFile = os.path.join(os.path.dirname(__file__), "plate.bdf")
    FEAAssembler = pyTACS(bdfFile, comm=comm, options=structOptions)

    def elemCallBack(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
        # Material properties
        rho = 1550.0  # Density kg/m^3
        E1 = 117.9e9  # Longitudinal modulus (Pa)
        E2 = 9.7e9  # Transverse modulus (Pa)
        nu12 = 0.04  # Major Poisson's ratio
        G12 = 4.8e9  # In-plane shear modulus (Pa)
        G13 = 4.8e9  # Out-of-plane shear modulus (Pa)
        G23 = G13  # Transverse shear modulus (Pa)
        Xt = 1647.8e6  # Longitudinal tensile strength (Pa)
        Xc = 1034.2e6  # Longitudinal compressive strength (Pa)
        Yt = 64.1e6  # Transverse tensile strength (Pa)
        Yc = 227.5e6  # Transverse compressive strength (Pa)
        S12 = 71.0e6  # In-plane shear strength (Pa)

        # Plate geometry
        tplate = 0.005  # 1 mm
        tMin = 0.0001  # 0.1 mm
        tMax = 0.05  # 5 cm
        tply = 0.125e-3  # 0.125 mm

        # Set up constitutive model
        ortho_prop = constitutive.MaterialProperties(
            rho=rho,
            E1=E1,
            E2=E2,
            nu12=nu12,
            G12=G12,
            G13=G13,
            G23=G23,
            Xt=Xt,
            Xc=Xc,
            Yt=Yt,
            Yc=Yc,
            S12=S12,
        )

        ortho_ply = constitutive.OrthotropicPly(
            tply, ortho_prop, max_strain_criterion=useVMFailure
        )

        if useVMFailure:
            if comm.rank == 0:
                print("Setting to use von Mises failure criterion")
            ortho_ply.setUseMaxStrainCriterion()
        else:
            if comm.rank == 0:
                print(
                    f"Setting to use Tsai-Wu failure criterion (modified: {useModifiedTsaiWu})"
                )
            ortho_ply.setUseTsaiWuCriterion()
            ortho_ply.setUseModifiedTsaiWu(useModifiedTsaiWu)

        if compLPModel == "LamParamSmeared":
            pfNums = np.arange(0, 3, dtype=np.intc) + dvNum + 1
            lpNums = np.arange(0, 2, dtype=np.intc) + pfNums[-1] + 1
            con = constitutive.LamParamSmearedShellConstitutive(
                ortho_ply,
                t=tplate,
                t_num=dvNum,
                min_t=tMin,
                max_t=tMax,
                f0=0.25,  # ply thickness fractions
                f45=0.5,
                f90=0.25,
                f0_num=pfNums[0],  # dv numbers for ply fractions
                f45_num=pfNums[1],
                f90_num=pfNums[2],
                min_f0=0.1,  # minimum ply thickness fractions
                min_f45=0.1,
                min_f90=0.1,
                W1=0.6,
                W3=0.6,
                W1_num=lpNums[0],
                W3_num=lpNums[1],
                ksWeight=100,
                epsilon=0.0,
            )
            scale = [100.0] + [1.0] * len(pfNums) + [1.0] * len(lpNums)

        elif compLPModel == "LamParamFull":
            lpNums = np.arange(0, 6, dtype=np.intc) + dvNum + 1  # 6
            con = constitutive.LamParamFullShellConstitutive(
                ortho_ply, tplate, dvNum, tMin, tMax, lpNums, 100.0
            )
            lp = 0.6 * np.ones(6)
            con.setLaminationParameters(lp)
            con.setNumFailAngles(12)
            scale = [100.0] + [1.0] * len(lp)

        refAxis = np.array([1.0, 1.0, 0.0])
        transform = elements.ShellRefAxisTransform(refAxis)

        # Create the element
        elem = elements.Quad4Shell(transform, con)

        return elem, scale

    # Set up constitutive objects and elements
    FEAAssembler.initialize(elemCallBack)

    # Setup static problem
    staticProb = FEAAssembler.createStaticProblem(
        name="plate", options=staticProblemOptions
    )

    # Add functions
    KSWeight = 100.0
    safetyFactor = 1.5
    staticProb.addFunction("mass", functions.StructuralMass)
    staticProb.addFunction(
        "KSFailure", functions.KSFailure, ksWeight=KSWeight, safetyFactor=safetyFactor
    )

    F = np.array([0.0, 1e2, 1e3, 0.0, 0.0, 0.0])
    nodeID = 481  # Center of plate
    staticProb.addLoadToNodes(nodeID, F, nastranOrdering=True)

    # Create lamination parameter constraints
    constraint = FEAAssembler.createLamParamFullConstraint("lam_param_con_full")
    allCompIDs = FEAAssembler.selectCompIDs()
    constraint.addConstraint("ALL", compIDs=allCompIDs, dvIndices=np.arange(1, 7), dvWeights=np.ones(6))


    return FEAAssembler, staticProb, constraint


if __name__ == "__main__":
    import argparse

    comm = MPI.COMM_WORLD

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--writeSol", action="store_true", help="Write solution to file"
    )
    parser.add_argument(
        "--compLPModel",
        type=str,
        choices=[
            "LamParamFull",
            "LamParamSmeared",
        ],
        default="LamParamFull",
        help="What type of composite lamination parameter model to use",
    )
    parser.add_argument(
        "--useVMFailure",
        action="store_true",
        help="Set von Mises failure (Tsai-Wu is default)",
    )
    parser.add_argument(
        "--useModifiedTsaiWu",
        action="store_true",
        help="Use the modified TW failure criterion",
    )

    args = parser.parse_args()

    FEAAssembler, staticProb, constraint = setup(
        compLPModel=args.compLPModel,
        useVMFailure=args.useVMFailure,
        useModifiedTsaiWu=args.useModifiedTsaiWu,
    )
    # Solve state
    staticProb.solve()

    # Evaluate functions
    funcs = {}
    staticProb.evalFunctions(funcs)
    constraint.evalConstraints(funcs)

    print(funcs)

    if comm.rank == 0:
        for funcType in ["mass", "failure"]:
            print("\n%s FUNCTIONS:" % funcType.upper())
            for func in funcs:
                if funcType in func.lower():
                    print("%s = %f" % (func, funcs[func]))

    if args.writeSol:
        # Create the output folder if it doesn't exist
        outputDir = pathlib.Path("output")
        outputDir.mkdir(exist_ok=True)

        # Construct the output filename

        outname = f"plate_composite_{args.compLPModel}"
        if args.useVMFailure:
            outname += "_VMF"
        else:
            if args.useModifiedTsaiWu:
                outname += "_ModTW"
            else:
                outname += "_TW"

        outname = os.path.join(outputDir, outname)
        staticProb.writeSolution(baseName=outname)
