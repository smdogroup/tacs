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
        "point_load_Ixx": 1.4589304315475893,
        "point_load_Ixy": 3.552713678800501e-14,
        "point_load_Ixz": 4.163336342344337e-17,
        "point_load_Iyy": 1.4589304315475804,
        "point_load_Iyz": -1.3877787807814457e-17,
        "point_load_Izz": 2.9166666666665915,
        "point_load_cgx": 0.5000000000000038,
        "point_load_cgy": 0.5000000000000041,
        "point_load_cgz": -0.0035714285714285696,
        "point_load_compliance": 72.89406192461088,
        "point_load_ks_failure": 1.5271816916601026,
        "point_load_mass": 17.49999999999987,
        "pressure_Ixx": 1.4589304315475893,
        "pressure_Ixy": 3.552713678800501e-14,
        "pressure_Ixz": 4.163336342344337e-17,
        "pressure_Iyy": 1.4589304315475804,
        "pressure_Iyz": -1.3877787807814457e-17,
        "pressure_Izz": 2.9166666666665915,
        "pressure_cgx": 0.5000000000000038,
        "pressure_cgy": 0.5000000000000041,
        "pressure_cgz": -0.0035714285714285696,
        "pressure_compliance": 386.08932419333286,
        "pressure_ks_failure": 2.0650026254455125,
        "pressure_mass": 17.49999999999987,
        "gravity_Ixx": 1.4589304315475893,
        "gravity_Ixy": 3.552713678800501e-14,
        "gravity_Ixz": 4.163336342344337e-17,
        "gravity_Iyy": 1.4589304315475804,
        "gravity_Iyz": -1.3877787807814457e-17,
        "gravity_Izz": 2.9166666666665915,
        "gravity_cgx": 0.5000000000000038,
        "gravity_cgy": 0.5000000000000041,
        "gravity_cgz": -0.0035714285714285696,
        "gravity_compliance": 11.378942561176144,
        "gravity_ks_failure": 0.35047628416141674,
        "gravity_mass": 17.49999999999987,
        "no_load_Ixx": 1.4589304315475893,
        "no_load_Ixy": 3.552713678800501e-14,
        "no_load_Ixz": 4.163336342344337e-17,
        "no_load_Iyy": 1.4589304315475804,
        "no_load_Iyz": -1.3877787807814457e-17,
        "no_load_Izz": 2.9166666666665915,
        "no_load_cgx": 0.5000000000000038,
        "no_load_cgy": 0.5000000000000041,
        "no_load_cgz": -0.0035714285714285696,
        "no_load_compliance": 0j,
        "no_load_ks_failure": 0.027465307216701942,
        "no_load_mass": 17.49999999999987,
        "modal_eigsm.0": 721146.5245587686,
        "modal_eigsm.1": 1566963.722841484,
        "modal_eigsm.2": 3742182.9051853963,
        "modal_eigsm.3": 5039286.163250271,
        "panel_length_con_ALL": np.array([0.0]),
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
                panelLength=self.dtype(lPanel),
                stiffenerPitch=self.dtype(pStiff),
                panelThick=self.dtype(tPlate),
                panelPlyAngles=plyAngles.astype(self.dtype),
                panelPlyFracs=plyFractions.astype(self.dtype),
                stiffenerHeight=self.dtype(hStiff),
                stiffenerThick=self.dtype(tStiff),
                stiffenerPlyAngles=plyAngles.astype(self.dtype),
                stiffenerPlyFracs=plyFractions.astype(self.dtype),
                panelLengthNum=dv_num,
                stiffenerPitchNum=dv_num + 1,
                panelThickNum=dv_num + 2,
                panelPlyFracNums=np.array([dv_num + 3], dtype=np.intc),
                stiffenerHeightNum=dv_num + 4,
                stiffenerThickNum=dv_num + 5,
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

        # No load
        sp = fea_assembler.createStaticProblem(name="no_load")
        tacs_probs.append(sp)

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("ks_failure", functions.KSFailure, ksWeight=ksweight)
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

    def test_panel_length_constraint(self):
        """
        Test the panel length constraint
        """
        constraint = self.tacs_probs[-1]
        self.assertFalse(
            constraint.isLinear,
            "Panel length constraint should not be linear wrt nodes",
        )
