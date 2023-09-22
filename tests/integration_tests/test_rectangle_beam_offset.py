import os

import numpy as np

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive, functions

"""
This is the same test cases as `test_rectangle_beam_tractions.py`, but the beam element is offset 
in the y and z directions, so that the beam axis no longer aligns with the nodes of the model. 
This test ensures that the beam solution is invariant under trivial transformation: 

6 noded beam model 1 meter long in x' direction.
The cross-section is a solid rectangle with the following properties:
    w = 0.1
    t = 0.05
We apply a distributed gravity load and a centrifugal load.
We test KSDisplacement, KSFailure, StructuralMass, CenterOfMass, MomentOfInertia, and Compliance 
functions and sensitivities.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/beam_model.bdf")

ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "gravity_I_xx": 0.0,
        "gravity_I_xy": 0.0,
        "gravity_I_xz": 0.0,
        "gravity_I_yy": 1.13625,
        "gravity_I_yz": 0.0,
        "gravity_I_zz": 1.1278125,
        "gravity_cgx": 0.5,
        "gravity_cgy": -0.025,
        "gravity_cgz": 0.05,
        "gravity_compliance": 2000.6857974986858,
        "gravity_ks_vmfailure": 219.37317664566595,
        "gravity_mass": 13.5,
        "gravity_x_disp": -0.7780125287429587,
        "gravity_y_disp": 656.1378755069089,
        "gravity_z_disp": 275.45079293282816,
        "centrifugal_I_xx": 0.0,
        "centrifugal_I_xy": 0.0,
        "centrifugal_I_xz": 0.0,
        "centrifugal_I_yy": 1.13625,
        "centrifugal_I_yz": 0.0,
        "centrifugal_I_zz": 1.1278125,
        "centrifugal_cgx": 0.5,
        "centrifugal_cgy": -0.025,
        "centrifugal_cgz": 0.05,
        "centrifugal_compliance": 9718.610964493391,
        "centrifugal_ks_vmfailure": 402.6071884894102,
        "centrifugal_mass": 13.5,
        "centrifugal_x_disp": 10.454782137288209,
        "centrifugal_y_disp": -18.840495836112435,
        "centrifugal_z_disp": -8.065586771046936,
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
            self.rtol = 2e-1
            self.atol = 1e-3
            self.dh = 1e-6

        # Material properties
        rho = 2700.0  # density kg/m^3
        E = 70.0e3  # Young's modulus (Pa)
        nu = 0.3  # Poisson's ratio
        ys = 2.7e3  # yield stress

        # Beam dimensions
        t = 0.05  # m
        w = 0.1  # m

        # Callback function used to setup TACS element objects and DVs
        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
        ):
            # Setup (isotropic) property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            con = constitutive.IsoRectangleBeamConstitutive(
                prop,
                t=t,
                tNum=dv_num,
                w=w,
                wNum=dv_num + 1,
                tOffset=0.5,
                wOffset=-0.5,
            )
            refAxis = np.array([0.0, 1.0, 0.0])
            transform = elements.BeamRefAxisTransform(refAxis)
            # pass back the appropriate tacs element object
            elem = elements.Beam2(transform, con)
            return elem

        # Instantiate FEA Assembler
        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        # Set up constitutive objects and elements
        fea_assembler.initialize(elem_call_back)

        grav_prob = fea_assembler.createStaticProblem("gravity")
        grav_prob.addInertialLoad([-10.0, 3.0, 5.0])

        rot_prob = fea_assembler.createStaticProblem("centrifugal")
        rot_prob.addCentrifugalLoad([10.0, 3.0, 5.0], [0.5, -0.025, 0.05])

        probs = [grav_prob, rot_prob]

        # Set convergence to be tight for test
        for problem in probs:
            problem.setOption("L2Convergence", 1e-15)
            problem.setOption("L2ConvergenceRel", 1e-15)

        # Add Functions
        for problem in probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("compliance", functions.Compliance)
            problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=ksweight)
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
