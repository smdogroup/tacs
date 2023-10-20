import os

import numpy as np

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive, functions

"""
Tests a smeared laminate shell model with the following layup: [0, 45, 30].
Two load cases are tested: an in-plane tension and out-of-plane shear.
This test is identical to test_shell_comp_unbalanced.py except since the laminate 
is smeared all stacking sequence dependence is neglected.
tests KSDisplacement, KSFailure, StructuralMass, CenterOfMass, MomentOfInertia, and Compliance functions
and sensitivities.
We also test a ply fraction summation constraint using the DVConstraint class.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/comp_plate.bdf")

ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "Tension_I_xx": 0.3875000544921874,
        "Tension_I_xy": 5.551115123125783e-17,
        "Tension_I_xz": 0.0,
        "Tension_I_yy": 0.02421880449218751,
        "Tension_I_yz": 0.0,
        "Tension_I_zz": 0.4117187499999999,
        "Tension_cgx": 0.25000000000000006,
        "Tension_cgy": 1.0000000000000002,
        "Tension_cgz": 0.0,
        "Tension_compliance": 8433.69262638969,
        "Tension_ks_TsaiWufailure": 4.35234046725213,
        "Tension_mass": 1.1624999999999999,
        "Tension_x_disp": 0.04604283873354243,
        "Tension_y_disp": -0.014684823847867035,
        "Tension_z_disp": 0.0,
        "VertShear_I_xx": 0.3875000544921874,
        "VertShear_I_xy": 5.551115123125783e-17,
        "VertShear_I_xz": 0.0,
        "VertShear_I_yy": 0.02421880449218751,
        "VertShear_I_yz": 0.0,
        "VertShear_I_zz": 0.4117187499999999,
        "VertShear_cgx": 0.25000000000000006,
        "VertShear_cgy": 1.0000000000000002,
        "VertShear_cgz": 0.0,
        "VertShear_compliance": 0.00010817284888043632,
        "VertShear_ks_TsaiWufailure": 0.22626387488408603,
        "VertShear_mass": 1.1624999999999999,
        "VertShear_x_disp": 0.0,
        "VertShear_y_disp": 0.0,
        "VertShear_z_disp": 0.0054598659927803965,
    }

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        """

        # Overwrite default check values
        if self.dtype == complex:
            self.rtol = 1e-6
            self.atol = 1e-6
            self.dh = 1e-50
        else:
            self.rtol = 1e-3
            self.atol = 1e-3
            self.dh = 1e-7

        fea_assembler = pytacs.pyTACS(bdf_file, comm)

        # Material properties
        rho = 1550.0
        specific_heat = 921.096
        E1 = 54e9
        E2 = 18e9
        nu12 = 0.25
        G12 = 9e9
        G13 = 9e9
        Xt = 2410.0e6
        Xc = 1040.0e6
        Yt = 73.0e6
        Yc = 173.0e6
        S12 = 71.0e6
        cte = 24.0e-6
        kappa = 230.0
        ply_thickness = 1.25e-3

        ortho_prop = constitutive.MaterialProperties(
            rho=rho,
            specific_heat=specific_heat,
            E1=E1,
            E2=E2,
            nu12=nu12,
            G12=G12,
            G13=G13,
            G23=G13,
            Xt=Xt,
            Xc=Xc,
            Yt=Yt,
            Yc=Yc,
            S12=S12,
            alpha=cte,
            kappa=kappa,
        )
        ortho_ply = constitutive.OrthotropicPly(ply_thickness, ortho_prop)

        # Shell thickness
        t = 7.5e-4  # m

        # Define reference axis for local shell stresses
        refAxis = np.array([1.0, 0.0, 0.0])
        transform = elements.ShellRefAxisTransform(refAxis)

        # Callback function used to setup TACS element objects and DVs
        def elemCallBack(
            dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs
        ):
            # Create the orthotropic layup
            ortho_layup = [ortho_ply] * 3
            ply_angles = np.deg2rad([0.0, 45.0, 30.0]).astype(self.dtype)
            ply_fractions = np.array([1.0 / 3.0] * 3, dtype=self.dtype)
            thickness_dv_num = dvNum
            ply_fraction_dv_nums = np.array(
                [dvNum + 1, dvNum + 2, dvNum + 3], dtype=np.intc
            )

            con = constitutive.SmearedCompositeShellConstitutive(
                ortho_layup,
                t,
                ply_angles,
                ply_fractions,
                thickness_dv_num=thickness_dv_num,
                ply_fraction_dv_nums=ply_fraction_dv_nums,
            )

            # For each element type in this component,
            # pass back the appropriate tacs element object
            elemList = []
            for elemDescript in elemDescripts:
                if elemDescript in ["CQUAD4", "CQUADR"]:
                    elem = elements.Quad4Shell(transform, con)
                elif elemDescript in ["CTRIA3", "CTRIAR"]:
                    elem = elements.Tri3Shell(transform, con)
                else:
                    print("Uh oh, '%s' not recognized" % (elemDescript))
                elemList.append(elem)

            return elemList

        # Set up constitutive objects and elements
        fea_assembler.initialize(elemCallBack)

        # Read in forces from BDF and create tacs struct problems
        tacs_probs = fea_assembler.createTACSProbsFromBDF()
        # Convert from dict to list
        tacs_probs = tacs_probs.values()
        # Set convergence to be tight for test
        for problem in tacs_probs:
            problem.setOption("L2Convergence", 1e-20)
            problem.setOption("L2ConvergenceRel", 1e-20)

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("compliance", functions.Compliance)
            problem.addFunction(
                "ks_TsaiWufailure",
                functions.KSFailure,
                ksWeight=ksweight,
                ftype="discrete",
            )
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
                "I_xx",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[1.0, 0.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_xy",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[0.0, 1.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_xz",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_yy",
                functions.MomentOfInertia,
                direction1=[0.0, 1.0, 0.0],
                direction2=[0.0, 1.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_yz",
                functions.MomentOfInertia,
                direction1=[0.0, 1.0, 0.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_zz",
                functions.MomentOfInertia,
                direction1=[0.0, 0.0, 1.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=True,
            )

        tacs_probs = list(tacs_probs)

        return tacs_probs, fea_assembler
