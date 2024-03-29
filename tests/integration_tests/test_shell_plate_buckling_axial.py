import os

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive

""""
The nominal case is a 4m x 3m flat plate under a buckling analysis. The
perimeter of the plate is pinned and loaded in compression on its horizontal edges. 
This tests the eigenvalues and eigenvalue sensitivities
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/ss_plate.bdf")


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "buckling_eigsb.0": 20.18364503170502,
        "buckling_eigsb.1": 48.93915410646456,
        "buckling_eigsb.2": 82.32683532479277,
        "buckling_eigsb.3": 86.45944604770253,
        "buckling_eigsb.4": 126.77071579688763,
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
            E = 205e9  # Young's modulus (Pa)
            nu = 0.3  # Poisson's ratio
            ys = 464.0e6  # yield stress

            # Plate geometry
            tplate = 0.020  # 20 mm

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
        buckle_prob.addLoadFromBDF(loadID=1)

        return [buckle_prob], fea_assembler
