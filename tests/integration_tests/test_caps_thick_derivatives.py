"""
Sean Engelstad, October 2022
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example
"""

import unittest, os, numpy as np, importlib
from tacs import functions, caps2tacs
from mpi4py import MPI

caps_loader = importlib.util.find_spec("pyCAPS")

# only run the test if pyCAPS can be imported
@unittest.skipIf(caps_loader is None, "skipping ESP/CAPS test without pyCAPS module")
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
            caps2tacs.ThicknessVariable(
                caps_group=f"rib{irib}", value=0.05, material=aluminum
            ).register_to(tacs_aim)
        for ispar in range(1, nspars + 1):
            caps2tacs.ThicknessVariable(
                caps_group=f"spar{ispar}", value=0.05, material=aluminum
            ).register_to(tacs_aim)
        caps2tacs.ThicknessVariable(
            caps_group="OML", value=0.03, material=aluminum
        ).register_to(tacs_aim)

        # add constraints and loads
        caps2tacs.PinConstraint("root").register_to(tacs_aim)
        caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=100).register_to(
            tacs_aim
        )
        # run the pre analysis to build tacs input files
        self.tacs_aim = tacs_aim.setup_aim().pre_analysis()
        self.fea_solver = self.tacs_aim.fea_solver.initialize()
        self.SPs = self.fea_solver.createTACSProbsFromBDF()
        for caseID in self.SPs:
            self.SPs[caseID].addFunction("mass", functions.StructuralMass)

    def _run_analysis(self):
        """
        run a complete forward and adjoint analysis
        """

        # solve the forward and adjoint analysis for each struct problem
        self._func_names = ["mass"]
        tacs_funcs = {}
        tacs_sens = {}
        for caseID in self.SPs:
            self.SPs[caseID].solve()
            self.SPs[caseID].evalFunctions(tacs_funcs, evalFuncs=self._func_names)
            self.SPs[caseID].evalFunctionsSens(tacs_sens, evalFuncs=self._func_names)
            self.SPs[caseID].writeSolution(
                baseName="tacs_output", outputDir=self.tacs_aim.analysis_dir
            )

        # functions and gradients are stored in the tacs AIM dynout method
        self._functions = {}
        self._gradients = {}
        tacs_key = list(tacs_funcs.keys())[0]
        self._functions["mass"] = tacs_funcs[tacs_key]
        self._gradients["mass"] = {}
        struct_sens = tacs_sens[tacs_key]["struct"]
        for ithick, thick_var in enumerate(self.tacs_aim.thickness_variables):
            self._gradients["mass"][thick_var.name] = struct_sens[ithick]

    def test_mass_thickness_derivatives(self):
        """
        test the shape derivatives from ESP/CAPS into TACS forward & adjoint analysis
        """

        # total derivative using adjoint & coordinate derivatives
        self._build_tacs_aim()
        self._run_analysis()

        # random dvar/ds contravariant tensor for the complex step perturbation
        dvar_ds = {
            thick_var.name: np.random.rand()
            for thick_var in self.tacs_aim.thickness_variables
        }
        adjoint_TD = 0.0
        for thick_var in self.tacs_aim.thickness_variables:
            adjoint_TD += (
                self._gradients["mass"][thick_var.name] * dvar_ds[thick_var.name]
            )

        # total derivative with finite difference
        h = 1.0e-5
        for caseID in self.SPs:
            xarray = self.SPs[
                caseID
            ].x.getArray()  # gets the xarray of last Struct Problem
            for ithick, thick_var in enumerate(self.tacs_aim.thickness_variables):
                xarray[ithick] += 1j * dvar_ds[thick_var.name] * h

        self._run_analysis()
        complex_step_TD = self._functions["mass"].imag / h

        # relative error of total derivatives
        rel_error = (adjoint_TD - complex_step_TD) / complex_step_TD
        print(
            "\Complex Step thickness derivative test with d(mass,ksfailure)/drib_a1..."
        )
        print(f"\tAdjoint TD = {adjoint_TD}")
        print(f"\tComplex Step TD = {complex_step_TD}")
        print(f"\trelative error = {rel_error}")

        tol = 1.0e-7  # tolerance not as strict since we only have FD not complex step available
        acceptable_error = abs(rel_error) < tol
        print(f"\ttest passed = {acceptable_error}")
        self.assertTrue(acceptable_error)


if __name__ == "__main__":
    unittest.main()
