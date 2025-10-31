import os
import unittest

import numpy as np

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive, functions, TACS

"""
|---------o---------o---------      | g    ---------o---------o---------|
  CompID 0| CompID 1| CompID 2      \/      CompID 2| CompID 1| CompID 0

2 3 noded beam models, each 1 meter long in x direction (shown above).
The beams are discretized into 3 separate cross sections that grow linearly in thickness
The cross-section is a solid rectangle with the following properties:
    w = 0.1
    t0 = 0.05
    deltat = 0.005
A distributed gravity load case is applied.
We apply apply various tip loads test KSDisplacement, StructuralMass, MomentOfInertia,
and Compliance functions and sensitivities.
We also apply an adjacency constraint on the difference between the thickness dvs of each cross-section.
"""

TACS_IS_COMPLEX = TACS.dtype == complex

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/two_beam.bdf")

ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "gravity_I_xx": 0.0,
        "gravity_I_xy": 0.0,
        "gravity_I_xz": 0.0,
        "gravity_I_yy": 30.99975033000031,
        "gravity_I_yz": 0.0,
        "gravity_I_zz": 30.982610954932227,
        "gravity_compliance": 5073.790903996116,
        "gravity_ks_vmfailure": 218.95861025129074,
        "gravity_mass": 29.700000000000003,
        "gravity_x_disp": -0.3633512018122108,
        "gravity_y_disp": 692.9103268396991,
        "gravity_z_disp": 300.4311098863023,
        "adjcon_beam": np.array([-0.005, -0.005]),
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
        t0 = 0.05  # m
        delta_t = 0.005
        w = 0.1  # m

        # Callback function used to setup TACS element objects and DVs
        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
        ):
            # Setup (isotropic) property and constitutive objects
            t = t0 + delta_t * comp_id
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            con = constitutive.IsoRectangleBeamConstitutive(
                prop, t=t, tNum=dv_num, w=w, wNum=dv_num + 1
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

        tacs_probs = [grav_prob]

        # Set convergence to be tight for test
        for problem in tacs_probs:
            problem.setOption("L2Convergence", 1e-20)
            problem.setOption("L2ConvergenceRel", 1e-20)

        # Add Functions
        for problem in tacs_probs:
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

        # Add constraint on difference between thickness dvsof each cross-section
        constr = fea_assembler.createAdjacencyConstraint("adjcon")
        allComps = fea_assembler.selectCompIDs()
        constr.addConstraint(
            "beam",
            dvIndex=1,
            compIDs=allComps,
        )
        tacs_probs.append(constr)

        return tacs_probs, fea_assembler

    # We have to skip these tests in complex mode because the beam
    # element uses complex step to approximate the Jacobian and this
    # leads to issues with complex stepping the sensitivities.
    @unittest.skipIf(TACS_IS_COMPLEX, "Skipping test_total_dv_sensitivities")
    def test_total_dv_sensitivities(self):
        super().test_total_dv_sensitivities()

    @unittest.skipIf(TACS_IS_COMPLEX, "Skipping test_total_xpt_sensitivities")
    def test_total_xpt_sensitivities(self):
        super().test_total_xpt_sensitivities()

    def test_adjacency(self):
        adjCon = self.tacs_probs[-1]
        expectedAdjacency = [(0, 1), (1, 2)]
        self.assertEqual(sorted(expectedAdjacency), sorted(adjCon.adjacentComps))


if __name__ == "__main__":
    unittest.main()
