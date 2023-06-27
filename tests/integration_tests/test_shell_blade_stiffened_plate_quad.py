import os

import numpy as np

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive, functions

"""
This case is identical to the shell_plate_quad.py case, except that it uses the BladeStiffenedShell constitutive model
and involves a panel length computation.

The nominal case is a 1m x 1m flat plate under three load cases:
a 10 kN point force at center, a 100kPa pressure applied to the surface, and a 100G gravity load. The
perimeter of the plate is fixed in all 6 degrees of freedom. The plate comprises
100 CQUAD4 elements and test KSFailure, StructuralMass, CenterOfMass, MomentOfInertia,
and Compliance functions and sensitivities. Finally, a modal analysis is performed
and the eigenvalues tested.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/plate.bdf")

# KS function weight
ksweight = 10.0

# Stiffener angle w.r.t. the x-axis
STIFFENER_ANGLE = np.deg2rad(30.0)
TRUE_PANEL_LENGTH = 1.0 / np.cos(STIFFENER_ANGLE)


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "point_load_Ixx": 1.4589036458332965,
        "point_load_Ixy": 3.907985046680551e-14,
        "point_load_Ixz": 0.0,
        "point_load_Iyy": 1.4589036458333053,
        "point_load_Iyz": 0.0,
        "point_load_Izz": 2.916666666666586,
        "point_load_cgx": 0.500000000000004,
        "point_load_cgy": 0.500000000000004,
        "point_load_cgz": 0.0,
        "point_load_compliance": 71.87161795633577,
        "point_load_ks_vmfailure": 0.3304517634314118,
        "point_load_mass": 17.5,
        "pressure_Ixx": 1.4589036458332965,
        "pressure_Ixy": 3.907985046680551e-14,
        "pressure_Ixz": 0.0,
        "pressure_Iyy": 1.4589036458333053,
        "pressure_Iyz": 0.0,
        "pressure_Izz": 2.916666666666586,
        "pressure_cgx": 0.500000000000004,
        "pressure_cgy": 0.500000000000004,
        "pressure_cgz": 0.0,
        "pressure_compliance": 377.11604579572935,
        "pressure_ks_vmfailure": 0.9085085771277308,
        "pressure_mass": 17.5,
        "gravity_Ixx": 1.4589036458332965,
        "gravity_Ixy": 3.907985046680551e-14,
        "gravity_Ixz": 0.0,
        "gravity_Iyy": 1.4589036458333053,
        "gravity_Iyz": 0.0,
        "gravity_Izz": 2.916666666666586,
        "gravity_cgx": 0.500000000000004,
        "gravity_cgy": 0.500000000000004,
        "gravity_cgz": 0.0,
        "gravity_compliance": 11.114479357783475,
        "gravity_ks_vmfailure": 0.0908244150089233,
        "gravity_mass": 17.5,
        "modal_eigsm.0": 729196.8586597344,
        "modal_eigsm.1": 1593583.2300563648,
        "modal_eigsm.2": 3815858.4718147167,
        "modal_eigsm.3": 5110483.515481135,
        "panel_length_con_ALL": [0.0],
    }

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        """

        # Overwrite default check values
        if self.dtype == complex:
            self.rtol = 1e-8
            self.atol = 1e-8
            self.dh = 1e-200
        else:
            self.rtol = 1e-3
            self.atol = 1e-4
            self.dh = 1e-8

        # Instantiate FEA Assembler
        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
        ):
            # Material properties
            rho = 2500.0  # density kg/m^3
            E = 70e9  # Young's modulus (Pa)
            nu = 0.3  # Poisson's ratio
            ys = 464.0e6  # yield stress

            # Plate geometry
            tPlate = 0.005  # 5 mm
            tStiff = tPlate  # Stiffener thickness 5 mm
            hStiff = 5 * tStiff  # Stiffener height
            pStiff = 5 * hStiff  # Stiffener pitch
            lPanel = TRUE_PANEL_LENGTH  # Panel length

            # Set up property model
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            ply = constitutive.OrthotropicPly(1.0, prop)

            # Set up constitutive model
            plyAngles = np.zeros(1)
            plyFractions = stiffenerPlyFractions = np.ones(1)

            con = constitutive.BladeStiffenedShellConstitutive(
                panelPly=ply,
                stiffenerPly=ply,
                kcorr=5.0 / 6.0,
                panelLength=lPanel,
                panelLengthNum=dv_num,
                stiffenerPitch=pStiff,
                stiffenerPitchNum=dv_num + 1,
                panelThick=tPlate,
                panelThickNum=dv_num + 2,
                numPanelPlies=len(plyAngles),
                panelPlyAngles=plyAngles,
                panelPlyFracs=plyFractions,
                panelPlyFracNums=np.array([dv_num + 3], dtype=np.intc),
                stiffenerHeight=hStiff,
                stiffenerHeightNum=dv_num + 4,
                stiffenerThick=tStiff,
                stiffenerThickNum=dv_num + 5,
                numStiffenerPlies=1,
                stiffenerPlyAngles=plyAngles,
                stiffenerPlyFracs=plyFractions,
                stiffenerPlyFracNums=np.array([dv_num + 6], dtype=np.intc),
            )

            # Set up reference axis along stiffener direction
            refAxis = np.array([np.cos(STIFFENER_ANGLE), np.sin(STIFFENER_ANGLE), 0.0])
            transform = elements.ShellRefAxisTransform(refAxis)

            # Set up element
            elem = elements.Quad4Shell(transform, con)
            return elem

        # Set up constitutive objects and elements
        fea_assembler.initialize(elem_call_back)

        tacs_probs = []

        # Add point force to node 81 (center of plate)
        sp = fea_assembler.createStaticProblem(name="point_load")
        F = np.array([0.0, 0.0, 1e4, 0.0, 0.0, 0.0])
        sp.addLoadToNodes(81, F, nastranOrdering=True)
        tacs_probs.append(sp)

        # Add pressure to entire plate
        sp = fea_assembler.createStaticProblem(name="pressure")
        P = 100e3  # Pa
        compIDs = fea_assembler.selectCompIDs(include="PLATE")
        sp.addPressureToComponents(compIDs, P)
        tacs_probs.append(sp)

        # Add pressure to entire plate
        sp = fea_assembler.createStaticProblem(name="gravity")
        g = np.array([0.0, 0.0, -981.0], dtype=self.dtype)
        sp.addInertialLoad(g)
        tacs_probs.append(sp)

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=ksweight)
            problem.addFunction("compliance", functions.Compliance)
            problem.addFunction(
                "cgx", functions.CenterOfMass, direction=[1.0, 0.0, 0.0]
            )
            problem.addFunction(
                "cgy", functions.CenterOfMass, direction=[0.0, 1.0, 0.0]
            )
            problem.addFunction(
                "cgz", functions.CenterOfMass, direction=[0.0, 0.0, 1.0]
            )
            problem.addFunction(
                "Ixx",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[1.0, 0.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "Ixy",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[0.0, 1.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "Ixz",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=True,
            )
            problem.addFunction(
                "Iyy",
                functions.MomentOfInertia,
                direction1=[0.0, 1.0, 0.0],
                direction2=[0.0, 1.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "Iyz",
                functions.MomentOfInertia,
                direction1=[0.0, 0.0, 1.0],
                direction2=[0.0, 1.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "Izz",
                functions.MomentOfInertia,
                direction1=[0.0, 0.0, 1.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=True,
            )

        # Add modal problem
        mp = fea_assembler.createModalProblem("modal", 7e4, 10)
        mp.setOption("printLevel", 2)
        tacs_probs.append(mp)

        # Create panel length constraint
        allCompIDs = fea_assembler.selectCompIDs()
        constraint = fea_assembler.createPanelLengthConstraint("panel_length_con")
        constraint.addConstraint("ALL", compIDs=allCompIDs)
        tacs_probs.append(constraint)

        return tacs_probs, fea_assembler
