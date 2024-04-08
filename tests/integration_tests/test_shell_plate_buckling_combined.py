import os

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive, TACS
import unittest

""""
The nominal case is a 1m x 0.7m flat plate under a buckling analysis. The
perimeter of the plate is pinned and loaded in combined axial and shear on its horizontal edges. 
This tests the eigenvalues and eigenvalue sensitivities
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/plate_combined_buckle.bdf")


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "buckle_eigsb.0": 52.61641147328781,
        "buckle_eigsb.1": 58.487138500669104,
        "buckle_eigsb.2": 117.21577758483404,
        "buckle_eigsb.3": 129.37327029030266,
        "buckle_eigsb.4": 143.2423322677979,
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
            scale = [100.0]
            return elem, scale

        # Set up constitutive objects and elements
        fea_assembler.initialize(elem_call_back)

        buckle_prob = fea_assembler.createBucklingProblem("buckling", 10.0, 10)
        buckle_prob.setOption("L2Convergence", 1e-20)
        buckle_prob.setOption("L2ConvergenceRel", 1e-20)
        # no loads just displacement control

        return [buckle_prob], fea_assembler


if __name__ == "__main__":
    import unittest

    unittest.main()
