import os

import numpy as np

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive, functions

"""
6 noded beam model 1 meter long in x direction.
The cross-section is a solid rectangle with the following properties:
    w = 0.1
    t = 0.05
We apply two load cases: a distributed gravity and distributed traction case.
We apply apply various tip loads test KSDisplacement, StructuralMass, MomentOfInertia, 
and Compliance functions and sensitivities.
We also apply a constraint on the difference between the width and thickness dvs of the cross-section.
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
        "gravity_compliance": 2000.6857974986824,
        "gravity_ks_vmfailure": 219.37317664566586,
        "gravity_mass": 13.5,
        "gravity_x_disp": -0.7780125287445699,
        "gravity_y_disp": 656.1378755073663,
        "gravity_z_disp": 275.45079293279235,
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

        # Set convergence to be tight for test
        grav_prob.setOption("L2Convergence", 1e-20)
        grav_prob.setOption("L2ConvergenceRel", 1e-20)

        # Add Functions
        grav_prob.addFunction("mass", functions.StructuralMass)
        grav_prob.addFunction("compliance", functions.Compliance)
        grav_prob.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=ksweight)
        grav_prob.addFunction(
            "x_disp",
            functions.KSDisplacement,
            ksWeight=ksweight,
            direction=[10.0, 0.0, 0.0],
        )
        grav_prob.addFunction(
            "y_disp",
            functions.KSDisplacement,
            ksWeight=ksweight,
            direction=[0.0, 10.0, 0.0],
        )
        grav_prob.addFunction(
            "z_disp",
            functions.KSDisplacement,
            ksWeight=ksweight,
            direction=[0.0, 0.0, 10.0],
        )
        grav_prob.addFunction("cgx", functions.CenterOfMass, direction=[1.0, 0.0, 0.0])
        grav_prob.addFunction("cgy", functions.CenterOfMass, direction=[0.0, 1.0, 0.0])
        grav_prob.addFunction("cgz", functions.CenterOfMass, direction=[0.0, 0.0, 1.0])
        grav_prob.addFunction(
            "I_xx",
            functions.MomentOfInertia,
            direction1=[1.0, 0.0, 0.0],
            direction2=[1.0, 0.0, 0.0],
            aboutCM=True,
        )
        grav_prob.addFunction(
            "I_xy",
            functions.MomentOfInertia,
            direction1=[1.0, 0.0, 0.0],
            direction2=[0.0, 1.0, 0.0],
            aboutCM=True,
        )
        grav_prob.addFunction(
            "I_xz",
            functions.MomentOfInertia,
            direction1=[1.0, 0.0, 0.0],
            direction2=[0.0, 0.0, 1.0],
            aboutCM=True,
        )
        grav_prob.addFunction(
            "I_yy",
            functions.MomentOfInertia,
            direction1=[0.0, 1.0, 0.0],
            direction2=[0.0, 1.0, 0.0],
            aboutCM=True,
        )
        grav_prob.addFunction(
            "I_yz",
            functions.MomentOfInertia,
            direction1=[0.0, 1.0, 0.0],
            direction2=[0.0, 0.0, 1.0],
            aboutCM=True,
        )
        grav_prob.addFunction(
            "I_zz",
            functions.MomentOfInertia,
            direction1=[0.0, 0.0, 1.0],
            direction2=[0.0, 0.0, 1.0],
            aboutCM=True,
        )

        return [grav_prob], fea_assembler
