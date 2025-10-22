import os

import numpy as np

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive, functions

"""
This is the same geometry as `test_shell_plate_quad.py`, but the plate is offset 
in the z direction, so that the shell plane no longer aligns with the nodes of the model.
Tests a smeared laminate shell model with the following layup: [0, 45, 30].
Agravity load in the x direction is applied.
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
        "gravity_I_xx": 0.38750005449218716,
        "gravity_I_xy": 5.551115123125783e-17,
        "gravity_I_xz": -2.710505431213761e-20,
        "gravity_I_yy": 0.02421880449218755,
        "gravity_I_yz": -2.168404344971009e-19,
        "gravity_I_zz": 0.4117187499999999,
        "gravity_cgx": 0.25000000000000006,
        "gravity_cgy": 1.0000000000000002,
        "gravity_cgz": -0.00037500000000000006,
        "gravity_compliance": 3.4205235434811567e-06,
        "gravity_ks_TsaiWufailure": 0.2259320352899525,
        "gravity_mass": 1.1624999999999999,
        "gravity_x_disp": 1.625976981601372e-06,
        "gravity_y_disp": -1.2736463024905141e-07,
        "gravity_z_disp": -0.000920864663486706,
        "centrifugal_I_xx": 0.38750005449218716,
        "centrifugal_I_xy": 5.551115123125783e-17,
        "centrifugal_I_xz": -2.710505431213761e-20,
        "centrifugal_I_yy": 0.02421880449218755,
        "centrifugal_I_yz": -2.168404344971009e-19,
        "centrifugal_I_zz": 0.4117187499999999,
        "centrifugal_cgx": 0.25000000000000006,
        "centrifugal_cgy": 1.0000000000000002,
        "centrifugal_cgz": -0.00037500000000000006,
        "centrifugal_compliance": 0.23430251098271168,
        "centrifugal_ks_TsaiWufailure": 0.24632761512178053,
        "centrifugal_mass": 1.1624999999999999,
        "centrifugal_x_disp": 0.0001570586015183079,
        "centrifugal_y_disp": 8.748943337216315e-05,
        "centrifugal_z_disp": 0.0966983775511979,
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
                t_offset=0.5,
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
        grav_prob = fea_assembler.createStaticProblem("gravity")
        grav_prob.addInertialLoad(1000 * [9.81, 0.0, 0.0])

        rot_prob = fea_assembler.createStaticProblem("centrifugal")
        rot_prob.addCentrifugalLoad([10.0, 3.0, 5.0], [0.25, 1.0, -0.000375])

        probs = [grav_prob, rot_prob]

        # Set convergence to be tight for test
        for problem in probs:
            problem.setOption("L2Convergence", 1e-15)
            problem.setOption("L2ConvergenceRel", 1e-15)

        # Add Functions
        for problem in probs:
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

        return probs, fea_assembler
