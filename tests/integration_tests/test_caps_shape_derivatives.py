"""
Sean Engelstad, October 2022
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example
"""

import unittest, os, numpy as np, importlib
from tacs import functions, caps2tacs, TACS
from mpi4py import MPI

caps_loader = importlib.util.find_spec("pyCAPS")
complex_mode = TACS.dtype == complex

# only run the test if pyCAPS can be imported
@unittest.skipIf(caps_loader is None or complex_mode, "skipping ESP/CAPS test without pyCAPS module or in real mode")
class TestCaps2Tacs(unittest.TestCase):
    def _build_tacs_aim(self):
        comm = MPI.COMM_WORLD
        csm_path = os.path.join("input_files", "simple_naca_wing.csm")
        caps_struct = caps2tacs.CapsStruct.build(csm_file=csm_path, comm=comm)
        tacs_aim = caps_struct.tacs_aim
        caps_struct.egads_aim.set_mesh(  # need a refined-enough mesh for the derivative test to pass
            edge_pt_min=15,
            edge_pt_max=20,
            global_mesh_size=0.1,
            max_surf_offset=0.01,
            max_dihedral_angle=5,
        ).register_to(
            tacs_aim
        )

        aluminum = caps2tacs.Isotropic.aluminum().register_to(tacs_aim)

        # setup the thickness design variables + automatic shell properties
        nribs = int(tacs_aim.get_config_parameter("nribs"))
        nspars = int(tacs_aim.get_config_parameter("nspars"))
        for irib in range(1, nribs + 1):
            caps2tacs.ShellProperty(
                caps_group=f"rib{irib}", material=aluminum, membrane_thickness=0.05
            ).register_to(tacs_aim)
        for ispar in range(1, nspars + 1):
            caps2tacs.ShellProperty(
                caps_group=f"spar{ispar}", material=aluminum, membrane_thickness=0.05
            ).register_to(tacs_aim)
        caps2tacs.ShellProperty(
            caps_group="OML", material=aluminum, membrane_thickness=0.03
        ).register_to(tacs_aim)

        # register one shape variable rib_a1
        caps2tacs.ShapeVariable("rib_a1").register_to(tacs_aim)

        # add constraints and loads
        caps2tacs.PinConstraint("root").register_to(tacs_aim)
        caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=100).register_to(
            tacs_aim
        )

        # run the pre analysis to build tacs input files
        self.tacs_aim_wrapper = tacs_aim.setup_aim()
        self.tacs_aim = tacs_aim.aim  # get actual aim from underneath the wrapper

    def _run_analysis(self):
        """
        run a complete forward and adjoint analysis
        """
        self.tacs_aim.preAnalysis()
        SPs = self.tacs_aim_wrapper.createTACSProbs()
        for caseID in SPs:
            SPs[caseID].addFunction("mass", functions.StructuralMass)
            SPs[caseID].addFunction(
                "ks_vmfailure", functions.KSFailure, safetyFactor=1.5, ksWeight=50.0
            )

        # solve the forward and adjoint analysis for each struct problem
        self._func_names = ["mass", "ks_vmfailure"]
        for caseID in SPs:
            SPs[caseID].solve()
            SPs[caseID].writeSensFile(
                evalFuncs=self._func_names,
                tacsAim=self.tacs_aim_wrapper,
            )

        # compute the shape derivatives in ESP/CAPS which reads the sens file
        self.tacs_aim.postAnalysis()

        # functions and gradients are stored in the tacs AIM dynout method

    def test_mass_shape_derivatives(self):
        """
        test the shape derivatives from ESP/CAPS into TACS forward & adjoint analysis
        """

        test_functions = ["mass", "ks_vmfailure"]
        dLdf = {func_key: np.random.rand() for func_key in test_functions}

        # total derivative using adjoint & coordinate derivatives
        self._build_tacs_aim()
        self._run_analysis()
        adjoint_TD = 0.0
        initial = {}
        for func in test_functions:
            adjoint_TD += dLdf[func] * self.tacs_aim.dynout[func].deriv("rib_a1")
            initial[func] = self.tacs_aim.dynout[func].value

        # total derivative with finite difference
        h = 1.0e-5
        self.tacs_aim.geometry.despmtr["rib_a1"].value += h
        self._run_analysis()
        finite_diff_TD = 0.0
        for func in test_functions:
            deriv = (self.tacs_aim.dynout[func].value - initial[func]) / h
            finite_diff_TD += dLdf[func] * deriv

        # relative error of total derivatives
        rel_error = (adjoint_TD - finite_diff_TD) / finite_diff_TD
        print("\nFD shape derivative test with d(mass,ksfailure)/drib_a1...")
        print(f"\tAdjoint TD = {adjoint_TD}")
        print(f"\tFinite Diff TD = {finite_diff_TD}")
        print(f"\trelative error = {rel_error}")

        self.assertTrue(abs(rel_error) < 1.0e-4)


if __name__ == "__main__":
    unittest.main()
