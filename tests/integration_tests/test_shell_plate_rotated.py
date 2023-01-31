import os

import numpy as np

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive, functions

"""
This is the same test cases as `test_shell_plate_quad.py`, but the plate is been rotated
about the y-axis by 45 degrees, so that it lies in a slant in the xz plane. This test ensures that the plate solution
is invariant under trivial transformation:
a 10 kN point force at center, a 100kPa pressure applied to the surface, and a 100G gravity load. The
perimeter of the plate is fixed in all 6 degrees of freedom. The plate comprises
100 CQUAD4 elements and test KSFailure, StructuralMass, CenterOfMass, MomentOfInertia,
and Compliance functions and sensitivities
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/slanted_plate.bdf")

from test_shell_plate_quad import ProblemTest as PT, ksweight

# Define rotated coordinate frame axes
x_prime = np.sqrt(0.5) * np.array([1.0, 0.0, 1.0])
y_prime = np.array([0.0, 1.0, 0.0])
z_prime = np.sqrt(0.5) * np.array([-1.0, 0.0, 1.0])


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = PT.FUNC_REFS

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
            self.atol = 1e-4
            self.dh = 1e-6

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
            tplate = 0.005  # 5 mm

            # Set up property model
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            # Set up constitutive model
            con = constitutive.IsoShellConstitutive(prop, t=tplate, tNum=dv_num)
            transform = None
            # Set up element
            elem = elements.Quad4Shell(transform, con)
            scale = [100.0]
            return elem, scale

        # Set up constitutive objects and elements
        fea_assembler.initialize(elem_call_back)

        tacs_probs = []

        # Add point force to node 81 (center of plate)
        sp = fea_assembler.createStaticProblem(name="point_load")
        F = np.zeros(6)
        F[:3] = 1e4 * z_prime
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
        g = -981.0 * z_prime
        sp.addInertialLoad(g)
        tacs_probs.append(sp)

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=ksweight)
            problem.addFunction("compliance", functions.Compliance)
            # Calculate cg and MOI in rotated coordinate frame
            problem.addFunction("cgx", functions.CenterOfMass, direction=x_prime)
            problem.addFunction("cgy", functions.CenterOfMass, direction=y_prime)
            problem.addFunction("cgz", functions.CenterOfMass, direction=z_prime)
            problem.addFunction(
                "Ixx",
                functions.MomentOfInertia,
                direction1=x_prime,
                direction2=x_prime,
                aboutCM=True,
            )
            problem.addFunction(
                "Ixy",
                functions.MomentOfInertia,
                direction1=x_prime,
                direction2=y_prime,
                aboutCM=True,
            )
            problem.addFunction(
                "Ixz",
                functions.MomentOfInertia,
                direction1=x_prime,
                direction2=z_prime,
                aboutCM=True,
            )
            problem.addFunction(
                "Iyy",
                functions.MomentOfInertia,
                direction1=y_prime,
                direction2=y_prime,
                aboutCM=True,
            )
            problem.addFunction(
                "Iyz",
                functions.MomentOfInertia,
                direction1=y_prime,
                direction2=z_prime,
                aboutCM=True,
            )
            problem.addFunction(
                "Izz",
                functions.MomentOfInertia,
                direction1=z_prime,
                direction2=z_prime,
                aboutCM=True,
            )

        return tacs_probs, fea_assembler

    def test_residual_scaling(self):
        """Test that the load scaling is working correctly for the point and pressure loads."""
        res = self.fea_assembler.createVec(asBVec=True)
        scaledRes = self.fea_assembler.createVec(asBVec=True)
        for problem in self.tacs_probs:
            with self.subTest(problem=problem.name):
                np.random.seed(1)
                # Check that the residual is zero if the states and load scale are both zero
                problem.loadScale = 0.0
                problem.zeroVariables()
                problem.getResidual(res)

                self.assertEqual(np.real(res.norm()), 0.0)

                # Check that the loadScale does linearly scale the external loads
                fullRes = problem.assembler.createVec()
                problem.setLoadScale(1.0)
                problem.getResidual(fullRes)

                loadScale = np.random.rand()
                problem.setLoadScale(loadScale)
                problem.getResidual(scaledRes)

                # scaledRes -= loadScale*fullRes should = 0
                scaledRes.axpy(-loadScale, fullRes)
                np.testing.assert_almost_equal(np.real(scaledRes.norm()), 0.0, 1e-12)
