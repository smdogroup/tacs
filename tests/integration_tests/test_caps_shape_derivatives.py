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

# use random seed to avoid random outlier failures in shape derivatives
np.random.seed(1234567)


# only run the test if pyCAPS can be imported
# runs on github workflow in real mode or offline in any mode
@unittest.skipIf(
    caps_loader is None,
    "skipping ESP/CAPS test without pyCAPS module",
)
class TestCaps2TacsShape(unittest.TestCase):
    def test_wing_shape_derivatives(self):
        """
        test the shape derivatives from ESP/CAPS into TACS forward & adjoint analysis
        """

        # build the tacs model with constraints, loads, properties, analysis functions, mesh, etc.
        comm = MPI.COMM_WORLD
        tacs_model = caps2tacs.TacsModel.build(
            csm_file=csm_path, comm=comm, problem_name="capsStruct1"
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
            caps2tacs.ShellProperty(
                caps_group=f"rib{irib}", material=aluminum, membrane_thickness=0.05
            ).register_to(tacs_model)
        for ispar in range(1, nspars + 1):
            caps2tacs.ShellProperty(
                caps_group=f"spar{ispar}", material=aluminum, membrane_thickness=0.05
            ).register_to(tacs_model)
        caps2tacs.ShellProperty(
            caps_group="OML", material=aluminum, membrane_thickness=0.03
        ).register_to(tacs_model)

        # include shape variables
        caps2tacs.ShapeVariable("rib_a1").register_to(tacs_model)
        # caps2tacs.ShapeVariable("rib_a2").register_to(tacs_model)

        # add constraints and loads
        caps2tacs.PinConstraint("root").register_to(tacs_model)
        caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=100).register_to(
            tacs_model
        )

        # add analysis functions
        caps2tacs.AnalysisFunction.mass().register_to(tacs_model)
        caps2tacs.AnalysisFunction.ksfailure(ksWeight=5.0).register_to(tacs_model)

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
        h = 1.0e-4
        for (
            shape_var
        ) in (
            tacs_model.shape_variables
        ):  # perturb the variables to affect update design
            shape_var.value += dxds[var.name] * h
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
        print("\nFD shape derivative test with d(mass, ksfailure)/drib_a1")
        print(f"\tAdjoint TD = {adjoint_TD}")
        print(f"\tFinite Diff TD = {finite_diff_TD}")
        print(f"\trelative error = {rel_error}")

        self.assertTrue(abs(rel_error) < 1.0e-4)


if __name__ == "__main__":
    unittest.main()
