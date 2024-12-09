import os

import numpy as np
import openmdao.api as om
from mphys.core import Multipoint, MPhysVariables
from mphys.scenarios import ScenarioStructural

import tacs.mphys
from openmdao_analysis_base_test import OpenMDAOTestCase
from tacs import elements, constitutive, functions

"""
This is a simple 1m by 2m plate made up of four quad shell elements.
The plate is structurally loaded under a 100G gravity load and a unit force,
"f_struct", is applied on on every node. The mass and KSFailure of the plate
are evaluated as outputs and have their partial and total sensitivities checked.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/debug_plate.bdf")

# Historical reference values for function outputs
FUNC_REFS = {
    "analysis.mass": 55.6,
    "analysis.Ixx": 74.13379667,
    "analysis.ks_vmfailure": 5.778130269059719,
    "analysis.adjacency.PANEL": [0.0, 0.0, 0.0, 0.0],
}

# Inputs to check total sensitivities wrt
wrt = [
    f"mesh.fea_mesh.{MPhysVariables.Structures.Mesh.COORDINATES}",
    "dv_struct",
    MPhysVariables.Structures.Loads.AERODYNAMIC,
]

# KS function weight
ksweight = 10.0

# Adjacency constraint bounds
adj_lb = -2.5e-3
adj_ub = 2.5e-3


class ProblemTest(OpenMDAOTestCase.OpenMDAOTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    def setup_problem(self, dtype):
        """
        Setup openmdao problem object we will be testing.
        """

        # Overwrite default tolerances
        if dtype == complex:
            self.rtol = 1e-7
            self.dh = 1e-50
        else:
            self.rtol = 1e-2
            self.dh = 1e-8

        # Callback function used to setup TACS element objects and DVs
        def element_callback(
            dv_num, comp_id, comp_descript, elem_descripts, special_dvs, **kwargs
        ):
            rho = 2780.0  # density, kg/m^3
            E = 73.1e9  # elastic modulus, Pa
            nu = 0.33  # poisson's ratio
            ys = 324.0e6  # yield stress, Pa
            thickness = 0.01
            min_thickness = 0.002
            max_thickness = 0.05

            # Setup (isotropic) property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            # Set one thickness dv for every component
            con = constitutive.IsoShellConstitutive(
                prop, t=thickness, tNum=dv_num, tlb=min_thickness, tub=max_thickness
            )

            # For each element type in this component,
            # pass back the appropriate tacs element object
            transform = None
            elem = elements.Quad4Shell(transform, con)

            return elem

        def problem_setup(scenario_name, fea_assembler, problem):
            """
            Helper function to add fixed forces and eval functions
            to structural problems used in tacs builder
            """
            # Set convergence to be tight for test
            problem.setOption("L2Convergence", 1e-20)
            problem.setOption("L2ConvergenceRel", 1e-20)

            # Add TACS Functions
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=ksweight)
            problem.addFunction(
                "Ixx",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[1.0, 0.0, 0.0],
            )

            # Add gravity load
            g = np.array([0.0, 0.0, -9.81]) * 100  # m/s^2
            problem.addInertialLoad(g)

        def constraint_setup(scenario_name, fea_assembler, constraints):
            """
            Helper function to setup tacs constraint classes
            """
            # Setup adjacency constraints for panel thicknesses
            constr = fea_assembler.createAdjacencyConstraint("adjacency")
            constr.addConstraint("PANEL", lower=-2.5e-3, upper=2.5e-3)
            constr_list = [constr]
            return constr_list

        tacs_builder = tacs.mphys.TacsBuilder(
            mesh_file=bdf_file,
            element_callback=element_callback,
            problem_setup=problem_setup,
            constraint_setup=constraint_setup,
            check_partials=True,
            coupling_loads=MPhysVariables.Structures.Loads.AERODYNAMIC,
            write_solution=False,
        )
        self.tacs_builder = tacs_builder

        class Top(Multipoint):
            def setup(self):
                tacs_builder.initialize(self.comm)

                dvs = self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])
                structDVs = tacs_builder.get_initial_dvs()
                dvs.add_output("dv_struct", structDVs)
                lb, ub = tacs_builder.get_dv_bounds()
                structDVScaling = tacs_builder.get_dv_scalers()
                self.add_design_var(
                    "dv_struct", lower=lb, upper=ub, scaler=structDVScaling
                )

                f_size = tacs_builder.get_ndof() * tacs_builder.get_number_of_nodes()
                forces = self.add_subsystem("forces", om.IndepVarComp(), promotes=["*"])
                forces.add_output(
                    MPhysVariables.Structures.Loads.AERODYNAMIC,
                    val=np.ones(f_size),
                    distributed=True,
                )

                self.add_subsystem("mesh", tacs_builder.get_mesh_coordinate_subsystem())
                self.mphys_add_scenario(
                    "analysis", ScenarioStructural(struct_builder=tacs_builder)
                )
                self.connect(
                    f"mesh.{MPhysVariables.Structures.Mesh.COORDINATES}",
                    f"analysis.{MPhysVariables.Structures.COORDINATES}",
                )
                self.connect("dv_struct", "analysis.dv_struct")
                self.connect(
                    MPhysVariables.Structures.Loads.AERODYNAMIC,
                    f"analysis.{MPhysVariables.Structures.Loads.AERODYNAMIC}",
                )

            def configure(self):
                tacs.mphys.utils.add_tacs_constraints(self.analysis)
                self.add_constraint("analysis.ks_vmfailure", upper=1.0)
                self.add_objective("analysis.mass")

        prob = om.Problem()
        prob.model = Top()

        return prob

    def setup_funcs(self):
        """
        Create a dict of functions to be tested and a list of inputs
        to test their sensitivities with respect to.
        """
        return FUNC_REFS, wrt

    # def test_add_tacs_constraints(self):
    #     # prob = self.setup_problem(dtype=float)
    #     # prob.setup()
    #     prob = self.prob
    #     constraints = prob.model.get_constraints()
    #     self.assertIn("analysis.adjacency.PANEL", constraints)
    #     adjacency = constraints["analysis.adjacency.PANEL"]
    #     self.assertTrue(adjacency["linear"])
    #     lower_bound = adjacency["lower"]
    #     upper_bound = adjacency["upper"]
    #     np.testing.assert_equal(lower_bound, adj_lb)
    #     np.testing.assert_equal(upper_bound, adj_ub)

    def test_get_tagged_indices(self):
        """
        Test the get_tagged_indices method
        """

        # We want to test that:
        # - For each comp_id, get_tagged_indices returns the same indices as `getLocalNodeIDsForComps`
        # - For a random set of the NASTRAN node IDs, get_tagged_indices returns the corresponding local indices
        # - For a mix of comp_ids and NASTRAN node IDs, get_tagged_indices returns the correct local indices
        FEAAssembler = self.tacs_builder.get_fea_assembler()
        compIDs = FEAAssembler.selectCompIDs()[0]
        meshloader = FEAAssembler.meshLoader

        for compID in compIDs:
            trueNodeIDs = FEAAssembler.getLocalNodeIDsForComps([compID])
            compName = FEAAssembler.getCompNames(compID)
            taggedIndIDs = self.tacs_builder.get_tagged_indices(compName)
            self.assertEqual(sorted(trueNodeIDs), sorted(taggedIndIDs))

            nastranIDs = meshloader.getGlobalNodeIDsForComps(
                [compID], nastranOrdering=True
            )
            taggedIndIDs = self.tacs_builder.get_tagged_indices(nastranIDs)
            self.assertEqual(sorted(trueNodeIDs), sorted(taggedIndIDs))

        # now test a mix of comp_ids and NASTRAN node IDs, we'll use the name of the first compID and the NASTRAN node
        # IDs of the last compID
        tags = FEAAssembler.getCompNames(
            compIDs[0]
        ) + meshloader.getGlobalNodeIDsForComps([compIDs[-1]], nastranOrdering=True)
        trueNodeIDs = FEAAssembler.getLocalNodeIDsForComps([compIDs[0], compIDs[-1]])
        taggedIndIDs = self.tacs_builder.get_tagged_indices(tags)
        self.assertEqual(sorted(trueNodeIDs), sorted(taggedIndIDs))


if __name__ == "__main__":
    import unittest

    unittest.main()
