import os

import numpy as np

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive, functions

"""

Test the lamination parameters constitutive model consisting of all lamination parameters.

The nominal case is a 1m x 1m flat plate under two load cases:
a point force at center with 1kN in z and 0.1kN in y and 10kPa pressure applied to the surface.
The perimeter of the plate is fixed in all 6 degrees of freedom. The plate comprises
100 CQUAD4 elements and test Displacement, KSFailure, StructuralMass, CenterOfMass, MomentOfInertia,
and Compliance functions and sensitivities.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/plate.bdf")

# KS function weight for failure and displacement
ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "point_load_Ixx": 0.6458494791666745,
        "point_load_Ixy": -9.769962616701378e-15,
        "point_load_Ixz": 0.0,
        "point_load_Iyy": 0.645849479166678,
        "point_load_Iyz": 0.0,
        "point_load_Izz": 1.2916666666666847,
        "point_load_cgx": 0.49999999999999767,
        "point_load_cgy": 0.49999999999999756,
        "point_load_cgz": 0.0,
        "point_load_compliance": 16.45700369546179,
        "point_load_ks_failure": 0.21387130020636574,
        "point_load_mass": 7.750000000000039,
        "point_load_x_disp": -1.4294835340049015e-07,
        "point_load_y_disp": 7.015914611511476e-07,
        "point_load_z_disp": 0.046852109681568097,
        "pressure_Ixx": 0.6458494791666745,
        "pressure_Ixy": -9.769962616701378e-15,
        "pressure_Ixz": 0.0,
        "pressure_Iyy": 0.645849479166678,
        "pressure_Iyz": 0.0,
        "pressure_Izz": 1.2916666666666847,
        "pressure_cgx": 0.49999999999999767,
        "pressure_cgy": 0.49999999999999756,
        "pressure_cgz": 0.0,
        "pressure_compliance": 112.03681803257064,
        "pressure_ks_failure": 0.5981207921407652,
        "pressure_mass": 7.750000000000039,
        "pressure_x_disp": -7.438494264988578e-16,
        "pressure_y_disp": -7.438494264988578e-16,
        "pressure_z_disp": 0.17668717737166387,
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
        fea_assembler = pytacs.pyTACS(bdf_file, comm, options={})

        # Material properties
        rho = 1550.0
        E1 = 54e9
        E2 = 18e9
        nu12 = 0.25
        G12 = 9e9
        G13 = 9e9
        G23 = G13
        Xt = 2410.0e6
        Xc = 1040.0e6
        Yt = 73.0e6
        Yc = 173.0e6
        S12 = 71.0e6

        # Plate geometry
        tplate = 0.005
        tMin = 0.0001
        tMax = 0.05
        ply_thickness = 0.125e-3  # 0.125 mm

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

        ortho_ply = constitutive.OrthotropicPly(ply_thickness, ortho_prop, max_strain_criterion=False)

        def elemCallBack(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
            # Lamination parameter design variable numbers
            lpNums = np.arange(0, 6, dtype=np.intc) + dvNum + 1
            con = constitutive.LamParamAllShellConstitutive(ortho_ply, tplate, dvNum, tMin, tMax, lpNums, 100.0)

            # Set initial lamination parameters and number of failure angles
            lp = 0.6 * np.ones(6)
            con.setLaminationParameters(lp)
            con.setNumFailAngles(12)

            refAxis = np.array([1.0, 1.0, 0.0])
            transform = elements.ShellRefAxisTransform(refAxis)

            # Set up element
            elem = elements.Quad4Shell(transform, con)

            return elem

        # Set up constitutive objects and elements
        fea_assembler.initialize(elemCallBack)

        tacs_probs = []

        # Add point force to node 81 (center of plate)
        sp = fea_assembler.createStaticProblem(name="point_load")
        F = np.array([0.0, 1e2, 1e3, 0.0, 0.0, 0.0])
        sp.addLoadToNodes(81, F, nastranOrdering=True)
        tacs_probs.append(sp)

        # Add pressure to entire plate
        sp = fea_assembler.createStaticProblem(name="pressure")
        P = 10e3  # Pa
        compIDs = fea_assembler.selectCompIDs(include="PLATE")
        sp.addPressureToComponents(compIDs, P)
        tacs_probs.append(sp)

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("ks_failure", functions.KSFailure, ksWeight=ksweight, safetyFactor=1.25)
            problem.addFunction("compliance", functions.Compliance)
            problem.addFunction(
                "x_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[10.0, 0.0, 0.0],
            )
            problem.addFunction(
                "y_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[0.0, 10.0, 0.0],
            )
            problem.addFunction(
                "z_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[0.0, 0.0, 10.0],
            )
            problem.addFunction("cgx", functions.CenterOfMass, direction=[1.0, 0.0, 0.0])
            problem.addFunction("cgy", functions.CenterOfMass, direction=[0.0, 1.0, 0.0])
            problem.addFunction("cgz", functions.CenterOfMass, direction=[0.0, 0.0, 1.0])
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

        return tacs_probs, fea_assembler
