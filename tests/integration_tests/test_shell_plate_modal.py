import numpy as np
import os
from tacs import pytacs, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

""""
The nominal case is a 1m x 1m flat plate under a modal analysis. The
perimeter of the plate is fixed in all 6 degrees of freedom. The plate comprises
100 CQUAD4 elements and test the eigenvalues and eigenvalue sensitivities

The plate is broken up into 4 components in the bdf file with the names 'PLATE.00', 'PLATE.01', etc.
The plate domains ordered as below:
           |
    3      |       2
           |
-----------|------------
           |
   0       |      1
           |
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/partitioned_plate.bdf")


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "modal_eigsm.0": 1396464.9023496218,
        "modal_eigsm.1": 6329530.425709786,
        "modal_eigsm.2": 6329530.425710541,
        "modal_eigsm.3": 13800023.12391336,
        "modal_eigsm.4": 23935675.376059294,
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

        modal_prob = fea_assembler.createModalProblem("modal", 1.0e6, 10)

        return [modal_prob], fea_assembler
