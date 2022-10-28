import numpy as np
import os
from tacs import pytacs, TACS, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

"""
Hemispherical shell constructed from mixed quad/tri shell elements. 
The shell is subjected to an inward pressure and is supported at the rim.
The loads are applied in two equivilent load cases through the bdf:
    1. Using a PLOAD2 card
    2. Using a PLOAD4 card
    
A third load case, not specified in the bdf, is also added where the sturcture 
is spun around its center at a constant angular velocity causing a centrifugal load.

tests StructuralMass, MomentOfInertia, KSFailure, KSDisplacement and Compliance functions and sensitivities
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/hemisphere.bdf")

omega = 2 * np.pi * np.array([0.0, 0.0, -10.0])
rotCenter = np.zeros(3)

ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "PLOAD4_cg_x": 0.0009653731820888509,
        "PLOAD4_cg_y": -9.14227063766091e-05,
        "PLOAD4_cg_z": 0.49758219135768283,
        "PLOAD4_I_xx": 721.1210796251873,
        "PLOAD4_I_xy": 0.02281438889615433,
        "PLOAD4_I_xz": -0.13557311929923765,
        "PLOAD4_I_yy": 718.9187282561999,
        "PLOAD4_I_yz": 0.037711775302945186,
        "PLOAD4_I_zz": 1152.580386827468,
        "PLOAD4_compliance": 279300158.48951936,
        "PLOAD4_ks_disp": 9.927842420503762,
        "PLOAD4_ks_vmfailure": 29.307629374994303,
        "PLOAD4_mass": 1737.357316694243,
        "PLOAD2_cg_x": 0.0009653731820888509,
        "PLOAD2_cg_y": -9.14227063766091e-05,
        "PLOAD2_cg_z": 0.49758219135768283,
        "PLOAD2_I_xx": 721.1210796251873,
        "PLOAD2_I_xy": 0.02281438889615433,
        "PLOAD2_I_xz": -0.13557311929923765,
        "PLOAD2_I_yy": 718.9187282561999,
        "PLOAD2_I_yz": 0.037711775302945186,
        "PLOAD2_I_zz": 1152.580386827468,
        "PLOAD2_compliance": 279300158.48951936,
        "PLOAD2_ks_disp": 9.927842420503762,
        "PLOAD2_ks_vmfailure": 29.307629374994303,
        "PLOAD2_mass": 1737.357316694243,
        "Centrifugal_cg_x": 0.0009653731820888509,
        "Centrifugal_cg_y": -9.14227063766091e-05,
        "Centrifugal_cg_z": 0.49758219135768283,
        "Centrifugal_I_xx": 721.1210796251873,
        "Centrifugal_I_xy": 0.022814388896140014,
        "Centrifugal_I_xz": -0.13557311929923765,
        "Centrifugal_I_yy": 718.9187282561999,
        "Centrifugal_I_yz": 0.037711775302945186,
        "Centrifugal_I_zz": 1152.580386827468,
        "Centrifugal_compliance": 303.4866002859211,
        "Centrifugal_ks_disp": 0.18580458183836876,
        "Centrifugal_ks_vmfailure": 0.21309095365567654,
        "Centrifugal_mass": 1737.357316694243,
    }

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        """

        # Overwrite default check values
        if self.dtype == complex:
            self.rtol = 1e-8
            self.atol = 1e-8
            self.dh = 1e-50
        else:
            self.rtol = 2e-1
            self.atol = 1e-3
            self.dh = 1e-6

        # Instantiate FEA Assembler
        struct_options = {}

        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
        ):
            # Material properties
            rho = 2780.0  # density kg/m^3
            E = 70e9  # Young's modulus (Pa)
            nu = 0.3  # Poisson's ratio
            ys = 270e6  # yield stress

            # Shell thickness
            t = 0.1  # m

            # Setup (isotropic) property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            # Set one thickness dv for every component
            con = constitutive.IsoShellConstitutive(prop, t=t, tNum=dv_num)

            # For each element type in this component,
            # pass back the appropriate tacs element object
            elem_list = []
            transform = None
            for elem_descript in elem_descripts:
                if elem_descript in ["CQUAD4", "CQUADR"]:
                    elem = elements.Quad4Shell(transform, con)
                elif elem_descript in ["CTRIA3", "CTRIAR"]:
                    elem = elements.Tri3Shell(transform, con)
                else:
                    print("Uh oh, '%s' not recognized" % (elem_descript))
                elem_list.append(elem)

            # Add scale for thickness dv
            scale = [100.0]
            return elem_list, scale

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        # Set up constitutive objects and elements
        fea_assembler.initialize(elem_call_back)

        # Read in forces from BDF and create tacs struct problems
        tacs_probs = fea_assembler.createTACSProbsFromBDF()
        # Convert from dict to list
        tacs_probs = list(tacs_probs.values())
        static_prob = fea_assembler.createStaticProblem("Centrifugal")
        static_prob.addCentrifugalLoad(omega, rotCenter)
        tacs_probs.append(static_prob)

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("compliance", functions.Compliance)
            problem.addFunction(
                "ks_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[-100.0, -100.0, -100.0],
            )
            problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=ksweight)
            problem.addFunction(
                "cg_x", functions.CenterOfMass, direction=[1.0, 0.0, 0.0]
            )
            problem.addFunction(
                "cg_y", functions.CenterOfMass, direction=[0.0, 1.0, 0.0]
            )
            problem.addFunction(
                "cg_z", functions.CenterOfMass, direction=[0.0, 0.0, 1.0]
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
                direction1=[0.0, 0.0, 1.0],
                direction2=[0.0, 1.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_zz",
                functions.MomentOfInertia,
                direction1=[0.0, 0.0, 1.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=True,
            )

        return tacs_probs, fea_assembler
