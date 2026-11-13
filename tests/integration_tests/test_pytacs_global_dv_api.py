import os
import unittest

import numpy as np

from mpi4py import MPI

from tacs import pytacs, constitutive, elements
from tacs.utilities import Error

try:
    import openmdao.api as om
    from tacs.mphys.dv import TacsDVComp
except ImportError:
    om = None

"""
Unit tests for the pyTACS addGlobalDV/assignMassDV API, covering both the
scalar and array-valued design variable paths (dv numbering, bounds/scale
broadcasting, error handling, and the mphys mass dv split), and for how the
global design variables are distributed across MPI ranks.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/point_mass.bdf")
plate_bdf_file = os.path.join(base_dir, "./input_files/partitioned_plate.bdf")


class GlobalDVAPITest(unittest.TestCase):
    N_PROCS = 1

    def setUp(self):
        self.fea_assembler = pytacs.pyTACS(bdf_file)

    def test_scalar_dv(self):
        dvNum = self.fea_assembler.addGlobalDV("alpha", 1.0, lower=0.0, upper=2.0)
        self.assertIsInstance(dvNum, int)
        self.assertEqual(dvNum, 0)
        dv_dict = self.fea_assembler.getGlobalDVs()["alpha"]
        self.assertIsInstance(dv_dict["num"], int)
        self.assertEqual(self.fea_assembler.getGlobalDVNums(), [0])
        self.assertEqual(self.fea_assembler.getTotalNumGlobalDVs(), 1)
        self.assertEqual(len(self.fea_assembler.scaleList), 1)

    def test_array_dv(self):
        dvNums = self.fea_assembler.addGlobalDV("g", [0.0, 0.0, -9.81])
        self.assertIsInstance(dvNums, np.ndarray)
        np.testing.assert_array_equal(dvNums, [0, 1, 2])
        self.assertEqual(self.fea_assembler.getGlobalDVNums(), [0, 1, 2])
        self.assertEqual(self.fea_assembler.getTotalNumGlobalDVs(), 3)
        self.assertEqual(len(self.fea_assembler.scaleList), 3)
        # Next dv should be numbered after the array
        nextNum = self.fea_assembler.addGlobalDV("alpha", 1.0)
        self.assertEqual(nextNum, 3)

    def test_array_dv_returned_nums_are_a_copy(self):
        dvNums = self.fea_assembler.addGlobalDV("g", [0.0, 0.0, -9.81])
        dvNums[0] = 99
        dv_dict = self.fea_assembler.getGlobalDVs()["g"]
        np.testing.assert_array_equal(dv_dict["num"], [0, 1, 2])

    def test_array_dv_scalar_bounds_broadcast(self):
        self.fea_assembler.addGlobalDV(
            "g", [0.0, 0.0, -9.81], lower=-20.0, upper=20.0, scale=0.1
        )
        dv_dict = self.fea_assembler.getGlobalDVs()["g"]
        np.testing.assert_array_equal(dv_dict["lowerBound"], [-20.0] * 3)
        np.testing.assert_array_equal(dv_dict["upperBound"], [20.0] * 3)
        self.assertEqual(self.fea_assembler.scaleList, [0.1] * 3)

    def test_array_dv_array_bounds(self):
        self.fea_assembler.addGlobalDV(
            "g",
            [0.0, 0.0, -9.81],
            lower=[-1.0, -2.0, -20.0],
            upper=[1.0, 2.0, 0.0],
            scale=[1.0, 1.0, 0.1],
        )
        dv_dict = self.fea_assembler.getGlobalDVs()["g"]
        np.testing.assert_array_equal(dv_dict["lowerBound"], [-1.0, -2.0, -20.0])
        np.testing.assert_array_equal(dv_dict["upperBound"], [1.0, 2.0, 0.0])
        self.assertEqual(self.fea_assembler.scaleList, [1.0, 1.0, 0.1])

    def test_array_dv_unbounded(self):
        self.fea_assembler.addGlobalDV("g", [0.0, 0.0, -9.81])
        dv_dict = self.fea_assembler.getGlobalDVs()["g"]
        self.assertIsNone(dv_dict["lowerBound"])
        self.assertIsNone(dv_dict["upperBound"])

    def test_empty_array_value_raises(self):
        with self.assertRaises(Error):
            self.fea_assembler.addGlobalDV("g", [])

    def test_multidim_array_value_raises(self):
        with self.assertRaises(Error):
            self.fea_assembler.addGlobalDV("g", np.zeros((2, 3)))

    def test_bound_length_mismatch_raises(self):
        with self.assertRaises(Error):
            self.fea_assembler.addGlobalDV("g", [0.0, 0.0, -9.81], lower=[-1.0, -2.0])

    def test_scale_length_mismatch_raises(self):
        with self.assertRaises(Error):
            self.fea_assembler.addGlobalDV("g", [0.0, 0.0, -9.81], scale=[1.0, 2.0])

    def test_array_bound_with_scalar_value_raises(self):
        with self.assertRaises(Error):
            self.fea_assembler.addGlobalDV("alpha", 1.0, lower=[0.0])

    def test_assign_mass_dv_array_entry(self):
        dvNums = self.fea_assembler.addGlobalDV("point_dvs", [20.0, 30.0])
        self.fea_assembler.assignMassDV("point_dvs", 1, "m", index=1)
        dv_dict = self.fea_assembler.getGlobalDVs()["point_dvs"]
        np.testing.assert_array_equal(dv_dict["isMassDV"], [False, True])
        # The mass element should be wired to the dv num/value of entry 1
        self.assertEqual(self.fea_assembler.massDVs[1]["mNum"], dvNums[1])
        self.assertEqual(self.fea_assembler.massDVs[1]["m"], 30.0)

    def test_assign_mass_dv_array_entry_bounds(self):
        self.fea_assembler.addGlobalDV(
            "point_dvs", [20.0, 30.0], lower=[1.0, 2.0], upper=[100.0, 200.0]
        )
        self.fea_assembler.assignMassDV("point_dvs", 1, "m", index=0)
        self.assertEqual(self.fea_assembler.massDVs[1]["mlb"], 1.0)
        self.assertEqual(self.fea_assembler.massDVs[1]["mub"], 100.0)

    def test_assign_mass_dv_array_requires_index(self):
        self.fea_assembler.addGlobalDV("point_dvs", [20.0, 30.0])
        with self.assertRaises(Error):
            self.fea_assembler.assignMassDV("point_dvs", 1, "m")

    def test_assign_mass_dv_array_index_out_of_range(self):
        self.fea_assembler.addGlobalDV("point_dvs", [20.0, 30.0])
        with self.assertRaises(Error):
            self.fea_assembler.assignMassDV("point_dvs", 1, "m", index=2)

    def test_assign_mass_dv_scalar_rejects_index(self):
        self.fea_assembler.addGlobalDV("point_mass", 20.0)
        with self.assertRaises(Error):
            self.fea_assembler.assignMassDV("point_mass", 1, "m", index=0)

    def test_assign_mass_dv_missing_key_rejects_index(self):
        with self.assertRaises(Error):
            self.fea_assembler.assignMassDV("point_mass", 1, "m", index=0)

    def test_assign_mass_dv_scalar_still_works(self):
        self.fea_assembler.addGlobalDV("point_mass", 20.0)
        self.fea_assembler.assignMassDV("point_mass", 1, "m")
        dv_dict = self.fea_assembler.getGlobalDVs()["point_mass"]
        self.assertTrue(dv_dict["isMassDV"])


class GlobalDVBoundsAfterInitializeTest(unittest.TestCase):
    """
    The bounds and values registered with addGlobalDV are only pushed into the
    design vectors during initialize, so these checks have to run on an
    initialized assembler. Zero-valued entries are the interesting case: they
    used to be dropped on the way into the vector, leaving the bound at whatever
    the elements reported.
    """

    N_PROCS = 1

    def test_scalar_zero_bounds_are_applied(self):
        fea_assembler = pytacs.pyTACS(bdf_file)
        # A mass dv gets its bounds from the constitutive object (+/-1e20 by
        # default), so a dropped zero bound here is visible
        dvNum = fea_assembler.addGlobalDV("point_mass", 20.0, lower=0.0, upper=0.0)
        fea_assembler.assignMassDV("point_mass", 1, "m")
        fea_assembler.initialize()

        xlb, xub = fea_assembler.getDesignVarRange()
        self.assertEqual(np.real(xlb[dvNum]), 0.0)
        self.assertEqual(np.real(xub[dvNum]), 0.0)

    def test_array_zero_bounds_and_values_are_applied(self):
        fea_assembler = pytacs.pyTACS(bdf_file)
        dvNums = fea_assembler.addGlobalDV(
            "g", [0.0, 0.0, -9.81], lower=[0.0, -1.0, -20.0], upper=[0.0, 1.0, 0.0]
        )
        fea_assembler.initialize()

        xlb, xub = fea_assembler.getDesignVarRange()
        np.testing.assert_array_equal(np.real(xlb[dvNums]), [0.0, -1.0, -20.0])
        np.testing.assert_array_equal(np.real(xub[dvNums]), [0.0, 1.0, 0.0])
        # A zero initial value must survive the same round trip
        x0 = fea_assembler.getOrigDesignVars()
        np.testing.assert_array_equal(np.real(x0[dvNums]), [0.0, 0.0, -9.81])

    def test_unheld_dv_without_bounds_gets_unbounded_defaults(self):
        """
        A global dv consumed only by aux element loads is held by no regular
        element at initialize, so getDesignVarRange leaves its entries at the
        zeros of createDesignVec. The unbounded defaults have to be filled in
        explicitly, otherwise the initial value violates its own bounds.
        """
        fea_assembler = pytacs.pyTACS(bdf_file)
        dvNum = fea_assembler.addGlobalDV("pressure", 100e5)
        fea_assembler.initialize()

        xlb, xub = fea_assembler.getDesignVarRange()
        self.assertEqual(np.real(xlb[dvNum]), -1e20)
        self.assertEqual(np.real(xub[dvNum]), 1e20)
        x0 = np.real(fea_assembler.getOrigDesignVars())
        self.assertLessEqual(np.real(xlb[dvNum]), x0[dvNum])
        self.assertLessEqual(x0[dvNum], np.real(xub[dvNum]))

    def test_unheld_dv_with_one_sided_bound_defaults_the_other_side(self):
        fea_assembler = pytacs.pyTACS(bdf_file)
        dvNum = fea_assembler.addGlobalDV("pressure", 100e5, lower=0.0)
        fea_assembler.initialize()

        xlb, xub = fea_assembler.getDesignVarRange()
        self.assertEqual(np.real(xlb[dvNum]), 0.0)
        self.assertEqual(np.real(xub[dvNum]), 1e20)

    def test_partially_held_array_dv_without_bounds(self):
        # Entry 0 is held by the point mass element, which reports its own
        # range (mass constitutive lower bound is 0, not -1e20), while entry 1
        # is held by nothing and needs the unbounded defaults
        fea_assembler = pytacs.pyTACS(bdf_file)
        dvNums = fea_assembler.addGlobalDV("point_dvs", [20.0, 30.0])
        fea_assembler.assignMassDV("point_dvs", 1, "m", index=0)
        fea_assembler.initialize()

        xlb, xub = fea_assembler.getDesignVarRange()
        np.testing.assert_array_equal(np.real(xlb[dvNums]), [0.0, -1e20])
        np.testing.assert_array_equal(np.real(xub[dvNums]), [1e20, 1e20])

    def test_zero_overwrites_a_nonzero_entry(self):
        """
        The two tests above only pin the observable behaviour; they cannot tell
        an insert from a non-zero insert, because an unclaimed global dv entry
        already defaults to zero. This checks the insert directly: a zero value
        has to overwrite whatever the elements put in the vector, which is what
        breaks if the values are inserted with INSERT_NONZERO_VALUES.
        """
        fea_assembler = pytacs.pyTACS(bdf_file)
        dvNum = fea_assembler.addGlobalDV("point_mass", 20.0)
        fea_assembler.assignMassDV("point_mass", 1, "m")
        fea_assembler.initialize()

        vec = fea_assembler.assembler.createDesignVec()
        vec.getArray()[:] = -1e20
        fea_assembler._setGlobalDVValues(vec, [dvNum], [0.0])
        self.assertEqual(np.real(vec.getArray()[dvNum]), 0.0)


class GlobalDVDistributionTest(unittest.TestCase):
    """
    The classes above run on one rank, where rank 0 owns every design variable
    and nothing is ever communicated. These tests run the same model on two
    ranks, so that the global dv values and bounds have to be inserted into a
    distributed design vector, and the load elements on the rank that does not
    own a global dv have to receive its value through the external dv
    distribution.

    Every assertion is made on quantities that are reduced over all ranks, so
    that a failure is raised on every rank.
    """

    N_PROCS = 2

    P = 1e5

    def setUp(self):
        self.comm = MPI.COMM_WORLD
        self.fea_assembler = pytacs.pyTACS(plate_bdf_file, self.comm)
        # The zero entries are the ones an INSERT_NONZERO_VALUES insert drops
        self.g_dv_nums = self.fea_assembler.addGlobalDV(
            "g", [0.0, 0.0, -9.81], lower=[0.0, -1.0, -20.0], upper=[0.0, 1.0, 0.0]
        )
        self.p_dv_num = self.fea_assembler.addGlobalDV("pressure", self.P)

        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
        ):
            prop = constitutive.MaterialProperties(
                rho=2500.0, E=70e9, nu=0.3, ys=464.0e6
            )
            con = constitutive.IsoShellConstitutive(
                prop, t=0.01, tNum=dv_num, tlb=0.001, tub=0.1
            )
            return elements.Quad4Shell(None, con)

        self.fea_assembler.initialize(elem_call_back)
        self.comp_ids = self.fea_assembler.selectCompIDs()

    def gather(self, local_array):
        """Assemble the global array from the entries owned by each rank"""
        return np.concatenate(self.comm.allgather(np.real(local_array)))

    def state_norm(self, local_states):
        return np.sqrt(
            self.comm.allreduce(np.sum(np.real(local_states) ** 2), op=MPI.SUM)
        )

    def double_pressure_dv(self, problem):
        """Double the pressure dv, which is owned by rank 0 like every global dv"""
        x = problem.getDesignVars()
        if self.comm.rank == 0:
            x[self.p_dv_num] *= 2.0
        problem.setDesignVars(x)

    def test_both_ranks_own_elements(self):
        # Otherwise nothing below depends on the dv distribution
        num_owned = self.comm.allgather(self.fea_assembler.assembler.getNumElements())
        self.assertTrue(all(num > 0 for num in num_owned))

    def test_values_and_bounds_are_inserted(self):
        xlb, xub = self.fea_assembler.getDesignVarRange()
        xlb = self.gather(xlb)
        xub = self.gather(xub)
        x0 = self.gather(self.fea_assembler.getOrigDesignVars())

        num_dvs = self.fea_assembler.getTotalNumDesignVars()
        self.assertEqual(len(x0), num_dvs)

        np.testing.assert_array_equal(xlb[self.g_dv_nums], [0.0, -1.0, -20.0])
        np.testing.assert_array_equal(xub[self.g_dv_nums], [0.0, 1.0, 0.0])
        np.testing.assert_array_equal(x0[self.g_dv_nums], [0.0, 0.0, -9.81])
        # The pressure dv has no bounds and is held by no regular element
        self.assertEqual(xlb[self.p_dv_num], -1e20)
        self.assertEqual(xub[self.p_dv_num], 1e20)
        self.assertEqual(x0[self.p_dv_num], self.P)
        # Inserting the global dvs must not disturb the element dvs, which are
        # owned by both ranks
        num_global_dvs = self.fea_assembler.getTotalNumGlobalDVs()
        np.testing.assert_array_equal(xlb[num_global_dvs:], 0.001)
        np.testing.assert_array_equal(xub[num_global_dvs:], 0.1)
        np.testing.assert_array_equal(x0[num_global_dvs:], 0.01)

    def test_static_load_dv_reaches_every_rank(self):
        problem = self.fea_assembler.createStaticProblem("pressure")
        problem.addPressureToComponents(
            self.comp_ids, self.P, pressureDVNums=self.p_dv_num
        )
        problem.solve()
        norm0 = self.state_norm(problem.getVariables())
        self.assertGreater(norm0, 0.0)

        # The response is linear in the pressure, and every element carries the
        # load, so the norm only doubles if every rank received the new value
        self.double_pressure_dv(problem)
        problem.solve()
        norm1 = self.state_norm(problem.getVariables())
        np.testing.assert_allclose(norm1, 2.0 * norm0, rtol=1e-8)

    def test_transient_load_dv_reaches_every_rank(self):
        # A transient problem owns one set of load elements per time step, and
        # updates them without going through the assembler
        num_steps = 3
        problem = self.fea_assembler.createTransientProblem(
            "pressure", tInit=0.0, tFinal=1e-3, numSteps=num_steps
        )
        for step in range(num_steps + 1):
            problem.addPressureToComponents(
                step, self.comp_ids, self.P, pressureDVNums=self.p_dv_num
            )
        problem.solve()
        norm0 = self.state_norm(problem.getVariables(num_steps)[1])
        self.assertGreater(norm0, 0.0)

        self.double_pressure_dv(problem)
        problem.solve()
        norm1 = self.state_norm(problem.getVariables(num_steps)[1])
        np.testing.assert_allclose(norm1, 2.0 * norm0, rtol=1e-8)


@unittest.skipIf(om is None, "openmdao is not installed")
class TacsDVCompArrayMassTest(unittest.TestCase):
    N_PROCS = 1

    def test_array_mass_dv_split(self):
        fea_assembler = pytacs.pyTACS(bdf_file)
        fea_assembler.addGlobalDV("point_dvs", [20.0, 9.81])
        fea_assembler.assignMassDV("point_dvs", 1, "m", index=0)
        fea_assembler.initialize()

        vals = np.real(fea_assembler.getOrigDesignVars()).astype(float)
        prob = om.Problem()
        prob.model.add_subsystem(
            "dv_comp",
            TacsDVComp(
                fea_assembler=fea_assembler,
                initial_dv_vals=vals,
                separate_mass_dvs=True,
            ),
        )
        prob.setup()
        # Only the mass-assigned entry should be split out of the struct vector
        np.testing.assert_array_equal(prob.get_val("dv_comp.dv_mass_point_dvs"), [20.0])
        np.testing.assert_array_equal(prob.get_val("dv_comp.dv_struct"), [9.81])
        # Full dv vector should be reassembled in the original tacs ordering
        prob.set_val("dv_comp.dv_mass_point_dvs", 25.0)
        prob.run_model()
        np.testing.assert_array_equal(prob.get_val("dv_comp.tacs_dvs"), [25.0, 9.81])


if __name__ == "__main__":
    unittest.main()
