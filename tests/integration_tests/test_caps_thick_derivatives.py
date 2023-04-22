"""
Sean Engelstad, Febuary 2022
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example
"""

import unittest, os, numpy as np, importlib
from tacs import caps2tacs, TACS
from mpi4py import MPI

caps_loader = importlib.util.find_spec("pyCAPS")

base_dir = os.path.dirname(os.path.abspath(__file__))
csm_path = os.path.join(base_dir, "input_files", "simple_naca_wing.csm")


# only run the test if pyCAPS can be imported
# runs on github workflow in real mode or offline in any mode
@unittest.skipIf(
    caps_loader is None,
    "skipping ESP/CAPS test without pyCAPS module",
)
class TestCaps2TacsSizing(unittest.TestCase):
    def test_thickness_derivatives(self):
        """
        test the thickness derivatives in TACS setup from ESP/CAPS
        """

        # build the tacs model with constraints, loads, properties, analysis functions, mesh, etc.
        comm = MPI.COMM_WORLD
        tacs_model = caps2tacs.TacsModel.build(
            csm_file=csm_path, comm=comm, problem_name="capsStruct2"
        )
        tacs_model.mesh_aim.set_mesh(  # need a refined-enough mesh for the derivative test to pass
            edge_pt_min=15,
            edge_pt_max=20,
            global_mesh_size=0.1,
            max_surf_offset=0.01,
            max_dihedral_angle=5,
        ).register_to(
            tacs_model
        )
        aluminum = caps2tacs.Isotropic.aluminum().register_to(tacs_model)

        # setup the thickness design variables + automatic shell properties
        nribs = int(tacs_model.get_config_parameter("nribs"))
        nspars = int(tacs_model.get_config_parameter("nspars"))
        for irib in range(1, nribs + 1):
            caps2tacs.ThicknessVariable(
                caps_group=f"rib{irib}", value=0.05, material=aluminum
            ).register_to(tacs_model)
        for ispar in range(1, nspars + 1):
            caps2tacs.ThicknessVariable(
                caps_group=f"spar{ispar}", value=0.05, material=aluminum
            ).register_to(tacs_model)
        caps2tacs.ThicknessVariable(
            caps_group="OML", value=0.03, material=aluminum
        ).register_to(tacs_model)

        # add constraints and loads
        caps2tacs.PinConstraint("root").register_to(tacs_model)
        caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=100).register_to(
            tacs_model
        )

        # add analysis functions
        caps2tacs.AnalysisFunction.mass().register_to(tacs_model)
        # caps2tacs.AnalysisFunction.ksfailure().register_to(tacs_model)

        # setup the tacs model
        tacs_model.setup()

        # perform the derivative test with these random covariant and contravariant tensors
        dLdf = {func_key: np.random.rand() for func_key in tacs_model.function_names}
        dxds = {var.name: np.random.rand() for var in tacs_model.variables}

        # total derivative using adjoint & coordinate derivatives
        tacs_model.update_design()  # update design with nothing to write in shape vars for next design
        tacs_model.pre_analysis()
        tacs_model.run_analysis()
        tacs_model.post_analysis()  # calls tacsAim.postAnalysis and reads tacsAim.dynout under the hood
        adjoint_TD = 0.0
        initial = {}
        for func in tacs_model.analysis_functions:
            initial[func.name] = func.value
            for var in tacs_model.variables:
                derivative = func.get_derivative(var)
                adjoint_TD += dLdf[func.name] * derivative * dxds[var.name]

        # total derivative with finite difference
        h = 1.0e-5
        for (
            var
        ) in tacs_model.variables:  # perturb the variables to affect update design
            var.value += dxds[var.name] * h
        tacs_model.update_design()
        tacs_model.pre_analysis()
        tacs_model.run_analysis()
        tacs_model.post_analysis()  # calls tacsAim.postAnalysis and reads tacsAim.dynout under the hood

        finite_diff_TD = 0.0
        for func in tacs_model.analysis_functions:
            deriv = (func.value - initial[func.name]) / h
            finite_diff_TD += dLdf[func.name] * deriv

        # relative error of total derivatives
        rel_error = (adjoint_TD - finite_diff_TD) / finite_diff_TD
        print(
            "\nFD shape derivative test with d(mass,ksfailure)/d(length,width,stiffHeight)..."
        )
        print(f"\tAdjoint TD = {adjoint_TD}")
        print(f"\tFinite Diff TD = {finite_diff_TD}")
        print(f"\trelative error = {rel_error}")

        self.assertTrue(abs(rel_error) < 1.0e-4)


if __name__ == "__main__":
    unittest.main()
