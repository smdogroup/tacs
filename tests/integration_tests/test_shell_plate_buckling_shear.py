import os

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive, TACS
import unittest

complex_mode = TACS.dtype == complex

""""
The nominal case is a 1m x 0.7m flat plate under a buckling analysis. The
perimeter of the plate is clamped and loaded in shear on its horizontal edges. 
This tests the eigenvalues and eigenvalue sensitivities.
This analysis was verified against Abaqus and closed-form solution with a finer mesh.
    However, this test is scaled down to a 5x5 mesh to make sure the analysis result is identical.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/plate_shear_buckle.bdf")


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "buckling_eigsb.0": -170.74198051427817,
        "buckling_eigsb.1": 170.7419805143329,
        "buckling_eigsb.2": 199.91780267055614,
        "buckling_eigsb.3": -199.91780267064914,
        "buckling_eigsb.4": -377.1783008654174,
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
            self.atol = 1e-4
            self.dh = 1e-5

        # turn on absolute value comparison since +- shear mode eigenvalues can switch order
        self.absolute_compare = True

        # Instantiate FEA Assembler
        fea_assembler = pytacs.pyTACS(bdf_file, comm)

        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
        ):
            # Material properties
            rho = 2500.0  # density kg/m^3
            E = 70e9  # Young's modulus (Pa)
            nu = 0.33  # Poisson's ratio
            ys = 464.0e6  # yield stress

            # Plate geometry
            tplate = 0.07  # 20 mm

            # Set up property model
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            # Set up constitutive model
            con = constitutive.IsoShellConstitutive(prop, t=tplate, tNum=dv_num)
            transform = None
            # Set up element
            elem = elements.Quad4Shell(transform, con)
            elem.setComplexStepGmatrix(True)
            scale = [100.0]
            return elem, scale

        # Set up constitutive objects and elements
        fea_assembler.initialize(elem_call_back)

        buckle_prob = fea_assembler.createBucklingProblem("buckling", 10.0, 10)
        buckle_prob.setOption("L2Convergence", 1e-20)
        buckle_prob.setOption("L2ConvergenceRel", 1e-20)
        # no loads just displacement control

        return [buckle_prob], fea_assembler

    @unittest.skipIf(
        not complex_mode,
        "test with Gmatrix only in complex mode until analytic one implemented",
    )
    def test_total_dv_sensitivities(self):
        """
        Test total dv sensitivity through adjoint against fd/cs
        """
        PyTACSTestCase.PyTACSTest.test_total_dv_sensitivities(self)

    @unittest.skipIf(
        not complex_mode,
        "test with Gmatrix only in complex mode until analytic one implemented",
    )
    def test_total_xpt_sensitivities(self):
        """
        Test total xpt sensitivity through adjoint against fd/cs
        """
        PyTACSTestCase.PyTACSTest.test_total_xpt_sensitivities(self)


if __name__ == "__main__":
    import unittest

    unittest.main()
