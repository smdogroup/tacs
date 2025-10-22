"""
==============================================================================

==============================================================================
@File    :   test_mphys_struct_hemisphere_nonlinear.py
@Date    :   2023/10/29
@Author  :   Alasdair Christison Gray
@Description : This integration test verifies that the TACS MPhys wrapper works correctly for nonlinear structural problems.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import os

# ==============================================================================
# External Python modules
# ==============================================================================
import openmdao.api as om
from mphys.core import Multipoint
from mphys.scenarios import ScenarioStructural
from mphys.core import MPhysVariables

import tacs.mphys
from openmdao_analysis_base_test import OpenMDAOTestCase

# ==============================================================================
# Extension modules
# ==============================================================================
# We import the setup functions and reference output values from the existing nonlinear
# hemisphere test so that we are comparing the output from the TACS MPhys wrapper to
# the output from standalone pyTACS.
from test_shell_hemisphere_nonlinear import (
    bdf_file,
    setupHemisphereProblem,
    elemCallBack,
    hemisphereProbRefFuncs,
)


# We need to rename the reference functions to match the names used in the TACS MPhys wrapper
FUNC_REFS = {
    name.replace("_", "."): value for name, value in hemisphereProbRefFuncs.items()
}

wrt = ["dv_struct", f"mesh.fea_mesh.{MPhysVariables.Structures.Mesh.COORDINATES}"]


class ProblemTest(OpenMDAOTestCase.OpenMDAOTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    def setup_problem(self, dtype):
        """
        Setup openmdao problem object we will be testing.
        """

        # Overwrite default tolerances
        if dtype == complex:
            self.rtol = 1e-7
            self.dh = 1e-200
        else:
            self.rtol = 1e-4
            self.dh = 1e-8

        def problem_setup(scenario_name, fea_assembler, problem):
            """
            Helper function to add fixed forces and eval functions
            to structural problems used in tacs builder
            """
            # Set convergence to be tight for test
            problem.setOption("L2Convergence", 1e-20)
            problem.setOption("L2ConvergenceRel", 1e-20)
            problem.setOption("resetBeforeSolve", True)
            return setupHemisphereProblem(fea_assembler, problem)

        class Top(Multipoint):
            def setup(self):
                tacs_builder = tacs.mphys.TacsBuilder(
                    mesh_file=bdf_file,
                    element_callback=elemCallBack,
                    problem_setup=problem_setup,
                    check_partials=True,
                    write_solution=False,
                )
                tacs_builder.initialize(self.comm)

                dvs = self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])
                structDVs = tacs_builder.get_initial_dvs()
                dvs.add_output("dv_struct", structDVs)
                lb, ub = tacs_builder.get_dv_bounds()
                structDVScaling = tacs_builder.get_dv_scalers()
                self.add_design_var(
                    "dv_struct", lower=lb, upper=ub, scaler=structDVScaling
                )

                self.add_subsystem("mesh", tacs_builder.get_mesh_coordinate_subsystem())
                self.mphys_add_scenario(
                    "RadialForces", ScenarioStructural(struct_builder=tacs_builder)
                )
                self.connect(
                    f"mesh.{MPhysVariables.Structures.Mesh.COORDINATES}",
                    f"RadialForces.{MPhysVariables.Structures.COORDINATES}",
                )
                self.connect("dv_struct", "RadialForces.dv_struct")

        prob = om.Problem()
        prob.model = Top()

        return prob

    def setup_funcs(self):
        """
        Create a dict of functions to be tested and a list of inputs
        to test their sensitivities with respect to.
        """
        return FUNC_REFS, wrt


if __name__ == "__main__":
    import unittest

    unittest.main()
