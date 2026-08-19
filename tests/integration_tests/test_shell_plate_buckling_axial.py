import os

import numpy as np
from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive

""""
The nominal case is a 4m x 3m flat plate under a buckling analysis. The
perimeter of the plate is pinned and loaded in compression on its horizontal edges.
This tests the eigenvalues and eigenvalue sensitivities

A second buckling problem adds an in-plane traction and a pressure whose
magnitudes are both driven by global design variables. This covers the
auxiliary-element load design variables on the buckling problem, which the
assembler only updates if the auxiliary elements are attached before the design
variables are set.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/ss_plate.bdf")

# In-plane traction applied to the DV-driven buckling problem, aligned with the
# compression direction so that it perturbs the prebuckling stress state (N/m^2)
trac = [0.0, -1000.0, 0.0]
# Pressure applied to the DV-driven buckling problem (Pa)
P = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "buckling_eigsb.0": 20.183645031683817,
        "buckling_eigsb.1": 48.9391541064531,
        "buckling_eigsb.2": 82.32683532475838,
        "buckling_eigsb.3": 86.45944604744642,
        "buckling_eigsb.4": 126.77071579672322,
        "buckling_load_dv_eigsb.0": 18.668016675879374,
        "buckling_load_dv_eigsb.1": 45.27733031028185,
        "buckling_load_dv_eigsb.2": 74.6746903353058,
        "buckling_load_dv_eigsb.3": 80.39165954694316,
        "buckling_load_dv_eigsb.4": 116.36148871460612,
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

        # Register global DVs for the auxiliary load magnitudes before
        # initializing the assembler
        traction_dv_nums = fea_assembler.addGlobalDV("traction", trac)
        pressure_dv_num = fea_assembler.addGlobalDV("pressure", P)

        # Set up constitutive objects and elements
        fea_assembler.initialize(elem_call_back)

        # The load DVs were registered without bounds, so initialize must give
        # them the unbounded defaults rather than the zeros of an empty design
        # vector, and every initial value must sit inside its own bounds
        lb, ub = fea_assembler.getDesignVarRange()
        lb = np.concatenate(comm.allgather(np.real(lb)))
        ub = np.concatenate(comm.allgather(np.real(ub)))
        load_dv_nums = np.append(traction_dv_nums, pressure_dv_num)
        np.testing.assert_array_equal(lb[load_dv_nums], -1e20)
        np.testing.assert_array_equal(ub[load_dv_nums], 1e20)
        x0 = np.concatenate(comm.allgather(np.real(fea_assembler.getOrigDesignVars())))
        self.assertTrue(np.all(lb <= x0))
        self.assertTrue(np.all(x0 <= ub))

        buckle_prob = fea_assembler.createBucklingProblem("buckling", 10.0, 10)
        buckle_prob.addLoadFromBDF(loadID=1)
        buckle_prob.setOption("L2Convergence", 1e-20)
        buckle_prob.setOption("L2ConvergenceRel", 1e-20)

        # Same problem, with auxiliary loads whose magnitudes are design
        # variables
        all_comp_ids = fea_assembler.selectCompIDs()
        dv_prob = fea_assembler.createBucklingProblem("buckling_load_dv", 10.0, 10)
        dv_prob.addLoadFromBDF(loadID=1)
        dv_prob.addTractionToComponents(
            all_comp_ids, trac, tractionDVNums=traction_dv_nums
        )
        dv_prob.addPressureToComponents(all_comp_ids, P, pressureDVNums=pressure_dv_num)
        dv_prob.setOption("L2Convergence", 1e-20)
        dv_prob.setOption("L2ConvergenceRel", 1e-20)

        return [buckle_prob, dv_prob], fea_assembler


if __name__ == "__main__":
    import unittest

    unittest.main()
