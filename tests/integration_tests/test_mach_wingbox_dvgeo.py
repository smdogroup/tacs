"""
This tests TACS MACH structural optimization capabilities.
The wingbox model that we will be using for this problem is based on the MDOLab tutorial wingbox example,
cantilevered, with a pressure load applied on the skin and a shear load applied at the tip. The wingbox is
discretized using one shell element per panel. This tests the KSFailure function and the panel length constraint.

This tests the MACH StructProblem object's DVGeo and design variable sensitivities.
"""

import os
import numpy as np
from mpi4py import MPI
import unittest

from tacs import pyTACS
from tacs import elements, constitutive, functions
from tacs.mach import StructProblem
from mach_struct_problem_base_test import MACHStructProblemTestCase

try:
    from pygeo import DVGeometry
except ImportError:
    DVGeometry = None

# Get file paths
base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/coarse_mdo_tutorial_wingbox.bdf")
ffd_file = os.path.join(base_dir, "./input_files/mdo_tutorial_ffd.fmt")


@unittest.skipIf(DVGeometry is None, "pygeo is not installed")
class TestMACHWingboxExample(MACHStructProblemTestCase.MACHStructProblemTest):
    """
    Test case for MACH StructProblem using wingbox shape optimization.
    """

    N_PROCS = 2

    # Reference values for regression testing
    FUNC_REFS = {
        "2.5gload_SKIN_ksFailure": 2.28686464844228,
        "2.5gload_SPAR_RIB_ksFailure": 0.8179439383682638,
        "panel_length_con_ALL": np.array(
            [
                0.03177571,
                -0.00062845,
                0.06500047,
                0.03166657,
                0.03166604,
                0.03166361,
                0.03166229,
                0.03166113,
                0.03166204,
                0.03166379,
                0.03166625,
                0.03166962,
                0.03167413,
                0.03168508,
                0.03168686,
                0.03169697,
                0.0317101,
                0.03172713,
                0.03177573,
                -0.00062845,
                0.06500047,
                0.03166656,
                0.03166605,
                0.0316636,
                0.03166229,
                0.03166113,
                0.03166204,
                0.03166379,
                0.03166625,
                0.03166963,
                0.03167413,
                0.03168508,
                0.03168685,
                0.03169697,
                0.0317101,
                0.03172714,
                0.2341235,
                0.1751916,
                0.0676358,
                0.05324,
                0.0393297,
                0.0252944,
                0.0111235,
                -0.0032017,
                -0.017696,
                -0.0323818,
                -0.0472824,
                -0.0624296,
                -0.0778564,
                -0.0936098,
                -0.1097423,
                -0.1263259,
                -0.1434478,
                -0.1612237,
                -0.0534736,
                -0.0133258,
                -0.0133258,
                -0.0303375,
                -0.0466459,
                -0.0629617,
                -0.0792843,
                -0.095616,
                -0.1119564,
                -0.128308,
                -0.1446702,
                -0.1610457,
                -0.1774342,
                -0.1938382,
                -0.2102597,
                -0.2266995,
                -0.2431603,
                -0.2596498,
                0.2141235,
                0.1551916,
                0.0476358,
                0.03324,
                0.0193297,
                0.0052944,
                -0.0088765,
                -0.0232017,
                -0.037696,
                -0.0523818,
                -0.0672824,
                -0.0824296,
                -0.0978564,
                -0.1136098,
                -0.1297423,
                -0.1463259,
                -0.1634478,
                -0.1812237,
                -0.1998083,
            ]
        ),
    }

    def setup_struct_problems(self, comm):
        """
        Setup MACH StructProblem objects for testing.
        """

        FEAAssembler = pyTACS(bdf_file, comm=comm)

        # ==============================================================================
        # Design variable values, bounds, and scaling factors
        # ==============================================================================
        # Panel length
        panelLengthMax = np.inf
        panelLengthMin = 0.0
        panelLengthScale = 1.0

        # Stiffener pitch
        stiffenerPitch = 0.15  # m
        stiffenerPitchMax = 0.5  # m
        stiffenerPitchMin = 0.05  # m
        stiffenerPitchScale = 1.0

        # Panel thickness
        panelThickness = 0.02  # m
        panelThicknessMax = 0.0065  # m
        panelThicknessMin = 0.002  # m
        panelThicknessScale = 100.0

        # Stiffener height
        stiffenerHeight = 0.05  # m
        stiffenerHeightMax = 0.1  # m
        stiffenerHeightMin = 0.002  # m
        stiffenerHeightScale = 10.0

        # Stiffener thickness
        stiffenerThickness = 0.006  # m
        stiffenerThicknessMax = 0.1  # m
        stiffenerThicknessMin = 0.002  # m
        stiffenerThicknessScale = 100.0

        FEAAssembler.addGlobalDV(
            "USkinStiffPitch",
            stiffenerPitch,
            stiffenerPitchMin,
            stiffenerPitchMax,
            stiffenerPitchScale,
        )
        FEAAssembler.addGlobalDV(
            "LSkinStiffPitch",
            stiffenerPitch,
            stiffenerPitchMin,
            stiffenerPitchMax,
            stiffenerPitchScale,
        )
        FEAAssembler.addGlobalDV(
            "LESparStiffPitch",
            stiffenerPitch,
            stiffenerPitchMin,
            stiffenerPitchMax,
            stiffenerPitchScale,
        )
        FEAAssembler.addGlobalDV(
            "TESparStiffPitch",
            stiffenerPitch,
            stiffenerPitchMin,
            stiffenerPitchMax,
            stiffenerPitchScale,
        )
        FEAAssembler.addGlobalDV(
            "RibStiffPitch",
            stiffenerPitch,
            stiffenerPitchMin,
            stiffenerPitchMax,
            stiffenerPitchScale,
        )

        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
        ):
            currDVNum = dv_num
            elemList = []
            DVScales = []

            # Setup (isotropic) property and constitutive objects
            rho = 2700.0
            E = 70e9
            nu = 0.3
            ys = 300e6
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            ply = constitutive.OrthotropicPly(1.0, prop)

            # Use a 0-deg biased layup for the skin and a +-45-deg biased layup spars and ribs.
            # Align the stiffeners in the skins with the trailing edge spar, and the stiffeners
            # in the spars and ribs vertically.
            # The panel length values I set here are approximate, to get the real values, you'd
            # need to run an optimization with panel length design variables and constraints.
            if "SKIN" in comp_descript:
                refAxis = np.array([0.0, 0.0, 1.0])
                plyAngles = np.array([0.0], dtype=self.dtype)
                panelPlyFractions = np.array([1.0], dtype=self.dtype)
                panelLength = 0.735
            else:
                plyAngles = np.array([0.0], dtype=self.dtype)
                panelPlyFractions = np.array([1.0], dtype=self.dtype)
                refAxis = np.array([0.0, 1.0, 0.0])
                if "RIB" in comp_descript:
                    panelLength = 0.38
                elif "SPAR" in comp_descript:
                    panelLength = 0.36

            # Always use the 0-deg biased layup for the stiffeners
            stiffenerPlyFractions = panelPlyFractions
            numPlies = len(plyAngles)

            # --- Setup DV numbering and scaling ---

            # The ordering of the DVs used by the BladeStiffenedShell model is:
            # - panel length
            # - stiffener pitch
            # - panel thickness
            # - panel ply fractions (not used in this case)
            # - stiffener height
            # - stiffener thickness
            # - stiffener ply fractions (not used in this case)
            panelLengthNum = currDVNum
            DVScales.append(panelLengthScale)
            currDVNum += 1

            if "U_SKIN" in comp_descript:
                stiffenerPitch = global_dvs["USkinStiffPitch"]["value"]
                stiffenerPitchNum = global_dvs["USkinStiffPitch"]["num"]
                stiffenerPitchMin = global_dvs["USkinStiffPitch"]["lowerBound"]
                stiffenerPitchMax = global_dvs["USkinStiffPitch"]["upperBound"]
            elif "L_SKIN" in comp_descript:
                stiffenerPitch = global_dvs["LSkinStiffPitch"]["value"]
                stiffenerPitchNum = global_dvs["LSkinStiffPitch"]["num"]
                stiffenerPitchMin = global_dvs["LSkinStiffPitch"]["lowerBound"]
                stiffenerPitchMax = global_dvs["LSkinStiffPitch"]["upperBound"]
            elif "LE_SPAR" in comp_descript:
                stiffenerPitch = global_dvs["LESparStiffPitch"]["value"]
                stiffenerPitchNum = global_dvs["LESparStiffPitch"]["num"]
                stiffenerPitchMin = global_dvs["LESparStiffPitch"]["lowerBound"]
                stiffenerPitchMax = global_dvs["LESparStiffPitch"]["upperBound"]
            elif "TE_SPAR" in comp_descript:
                stiffenerPitch = global_dvs["TESparStiffPitch"]["value"]
                stiffenerPitchNum = global_dvs["TESparStiffPitch"]["num"]
                stiffenerPitchMin = global_dvs["TESparStiffPitch"]["lowerBound"]
                stiffenerPitchMax = global_dvs["TESparStiffPitch"]["upperBound"]
            elif "RIB" in comp_descript:
                stiffenerPitch = global_dvs["RibStiffPitch"]["value"]
                stiffenerPitchNum = global_dvs["RibStiffPitch"]["num"]
                stiffenerPitchMin = global_dvs["RibStiffPitch"]["lowerBound"]
                stiffenerPitchMax = global_dvs["RibStiffPitch"]["upperBound"]
            else:
                raise ValueError(f"Invalid component descriptor: {comp_descript}")

            panelThicknessNum = currDVNum
            DVScales.append(panelThicknessScale)
            currDVNum += 1

            stiffenerHeightNum = currDVNum
            DVScales.append(stiffenerHeightScale)
            currDVNum += 1

            stiffenerThicknessNum = currDVNum
            DVScales.append(stiffenerThicknessScale)
            currDVNum += 1

            con = constitutive.BladeStiffenedShellConstitutive(
                panelPly=ply,
                stiffenerPly=ply,
                panelLength=panelLength,
                stiffenerPitch=stiffenerPitch,
                panelThick=panelThickness,
                panelPlyAngles=plyAngles,
                panelPlyFracs=panelPlyFractions,
                stiffenerHeight=stiffenerHeight,
                stiffenerThick=stiffenerThickness,
                stiffenerPlyAngles=plyAngles,
                stiffenerPlyFracs=stiffenerPlyFractions,
                panelLengthNum=panelLengthNum,
                stiffenerPitchNum=stiffenerPitchNum,
                panelThickNum=panelThicknessNum,
                stiffenerHeightNum=stiffenerHeightNum,
                stiffenerThickNum=stiffenerThicknessNum,
            )
            con.setStiffenerPitchBounds(stiffenerPitchMin, stiffenerPitchMax)
            con.setPanelThicknessBounds(panelThicknessMin, panelThicknessMax)
            con.setStiffenerHeightBounds(stiffenerHeightMin, stiffenerHeightMax)
            con.setStiffenerThicknessBounds(
                stiffenerThicknessMin, stiffenerThicknessMax
            )

            # --- Create reference axis transform to define the stiffener direction ---
            transform = elements.ShellRefAxisTransform(refAxis)

            for elemDescript in elem_descripts:
                if elemDescript == "CQUAD4":
                    elem = elements.Quad4Shell(transform, con)
                    elemList.append(elem)
                elif elemDescript == "CTRIA3":
                    elem = elements.Tri3Shell(transform, con)
                    elemList.append(elem)

            return elemList, DVScales

        FEAAssembler.initialize(elem_call_back)

        staticProblem = FEAAssembler.createStaticProblem("2.5gload")

        F = np.array([0.0, 3e5, 0.0, 0.0, 0.0, 0.0])
        compIDs = FEAAssembler.selectCompIDs(["WING_RIB_018"])
        staticProblem.addLoadToComponents(compIDs, F, averageLoad=True)

        compIDs = FEAAssembler.selectCompIDs(["WING_L_SKIN"])
        staticProblem.addPressureToComponents(compIDs, -10e3)

        failureGroups = {"SPAR_RIB": ["SPAR", "RIB", "CAP"], "SKIN": ["SKIN"]}

        for failName in failureGroups:
            compIDs = FEAAssembler.selectCompIDs(
                include=failureGroups[failName], exclude="CAP"
            )
            staticProblem.addFunction(
                f"{failName}_ksFailure",
                functions.KSFailure,
                compIDs=compIDs,
                safetyFactor=1.5,
                ksWeight=10.0,
            )
        allCompIDs = FEAAssembler.selectCompIDs()
        panelLengthCon = FEAAssembler.createPanelLengthConstraint("panel_length_con")
        panelLengthCon.addConstraint("ALL", compIDs=allCompIDs)
        # ==============================================================================
        # Create DVGeometry
        # ==============================================================================
        DVGeo = DVGeometry(fileName=ffd_file, isComplex=self.dtype == complex)

        # Create reference axis
        nRefAxPts = DVGeo.addRefAxis("wing", xFraction=0.25, alignIndex="k")
        nTwist = nRefAxPts - 1

        # Set up global design variables
        def twist(val, geo):
            for i in range(1, nRefAxPts):
                geo.rot_z["wing"].coef[i] = val[i - 1]

        DVGeo.addGlobalDV(
            dvName="twist", value=[0] * nTwist, func=twist, lower=-10, upper=10, scale=1
        )

        # ==============================================================================
        # Create MACH StructProblem
        # ==============================================================================
        structProb = StructProblem(staticProblem, FEAAssembler, DVGeo=DVGeo)
        structProb.addConstraint(panelLengthCon)

        return [structProb]
