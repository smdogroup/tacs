"""
Tests for tractions that have fewer components than the element expects.

A thermoelastic solid has 4 variables per node, so its traction has 4
components: 3 mechanical tractions and a heat flux. A traction with only 3
components (which is what a PLOAD4 with a direction vector gives) is taken to
apply no heat flux, and the user is warned once.
"""

import contextlib
import io
import os
import unittest

import numpy as np
from mpi4py import MPI

from tacs import pytacs, constitutive, elements

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/two_hexs.bdf")

WARNING_TEXT = "missing components were set to zero"


class ShortTractionTest(unittest.TestCase):
    N_PROCS = 2

    def setUp(self):
        self.comm = MPI.COMM_WORLD
        self.fea_assembler = pytacs.pyTACS(bdf_file, self.comm)
        self.trac_dv_nums = self.fea_assembler.addGlobalDV("trac", [1.0, 2.0, 3.0])

        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
        ):
            props = constitutive.MaterialProperties(
                rho=2700.0,
                specific_heat=921.0,
                E=70e9,
                nu=0.3,
                ys=270e6,
                alpha=24e-6,
                kappa=230.0,
            )
            con = constitutive.SolidConstitutive(props, t=1.0, tNum=dv_num)
            model = elements.LinearThermoelasticity3D(con)
            return elements.Element3D(model, elements.LinearHexaBasis())

        self.fea_assembler.initialize(elem_call_back)
        self.problem = self.fea_assembler.createStaticProblem("traction")
        # The model has two elements
        self.elem_ids = np.arange(2)
        self.comp_ids = self.fea_assembler.selectCompIDs()

    def capture_warnings(self, func, *args, **kwargs):
        """
        Call a function and return what it printed on the root rank, which is
        where the pyTACS warnings are printed, on every rank. The box the
        warnings are wrapped in is removed.
        """
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            func(*args, **kwargs)
        output = self.comm.bcast(stdout.getvalue(), root=0)
        return " ".join(output.replace("|", " ").split())

    def test_element_expects_four_components(self):
        for elem_id in self.elem_ids:
            elem = self.fea_assembler.meshLoader.getElementObjectForElemID(elem_id)
            self.assertEqual(elem.getNumTractionComponents(), 4)

    def test_short_traction_warns_once(self):
        output = self.capture_warnings(
            self.problem.addTractionToElements,
            self.elem_ids,
            [1.0, 2.0, 3.0],
            tractionDVNums=self.trac_dv_nums,
        )
        self.assertEqual(output.count(WARNING_TEXT), 1)
        self.assertIn("Tractions with 3 components", output)
        self.assertIn(f"added to {len(self.elem_ids)} elements", output)
        self.assertIn("expect 4", output)

    def test_short_traction_to_components_warns_once(self):
        output = self.capture_warnings(
            self.problem.addTractionToComponents, self.comp_ids, [1.0, 2.0, 3.0]
        )
        self.assertEqual(output.count(WARNING_TEXT), 1)

    def test_full_traction_does_not_warn(self):
        output = self.capture_warnings(
            self.problem.addTractionToElements, self.elem_ids, [1.0, 2.0, 3.0, 4.0]
        )
        self.assertNotIn(WARNING_TEXT, output)

    def test_short_traction_dvs(self):
        """
        The design variables of the short traction are exactly the ones given,
        the padded heat flux is not a design variable.
        """
        self.capture_warnings(
            self.problem.addTractionToElements,
            self.elem_ids,
            [1.0, 2.0, 3.0],
            tractionDVNums=self.trac_dv_nums,
        )
        x = self.problem.getDesignVars()
        if self.comm.rank == 0:
            x[self.trac_dv_nums] = [4.0, 5.0, 6.0]
        self.problem.setDesignVars(x)
        x = np.concatenate(self.comm.allgather(np.real(self.problem.getDesignVars())))
        np.testing.assert_array_equal(x[self.trac_dv_nums], [4.0, 5.0, 6.0])


if __name__ == "__main__":
    unittest.main()
