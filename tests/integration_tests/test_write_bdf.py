import os
import tempfile
import unittest
import numpy as np
from pyNastran.bdf.bdf import BDF
from mpi4py import MPI

from tacs import pytacs

"""
This file tests pyTACS's `writeBDF` method.
We first instantiate pyTACS and structural problems from a provided BDF, as usual.
Using that instance of pyTACS we then export a BDF file using the `writeBDF` method, once with no scaling
arguments, testing the defaults, and a second time with scaling factors.
We then load both files and compare that they are consistent when accounting for the scaling factors.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))


class ProblemTest(unittest.TestCase):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    # Define test cases as class attribute
    TEST_CASES = [
        {
            "name": "cylinder",
            "bdf_file": os.path.join(base_dir, "./input_files/cylinder.bdf"),
            "xyz_scale": 1000.0,  # m to mm
            "time_scale": 100.0,  # s to cs
            "mass_scale": 1000.0,  # kg to g
        },
        {
            "name": "plate",
            "bdf_file": os.path.join(base_dir, "./input_files/plate.bdf"),
            "xyz_scale": 1000.0,  # m to mm
            "time_scale": 100.0,  # s to cs
            "mass_scale": 1.0,  # no mass scaling
        },
        {
            "name": "composite plate",
            "bdf_file": os.path.join(base_dir, "./input_files/comp_plate.bdf"),
            "xyz_scale": 1000.0,  # m to mm
            "time_scale": 100.0,  # s to cs
            "mass_scale": 1.0,  # no mass scaling
        },
        {
            "name": "beam",
            "bdf_file": os.path.join(base_dir, "./input_files/beam_model.bdf"),
            "xyz_scale": 100.0,  # Different scaling
            "time_scale": 100.0,  # s to cs
            "mass_scale": 1.0,  # no mass scaling
        },
    ]

    def setUp(self):
        self.comm = MPI.COMM_WORLD

    def setup_tacs(self, comm, orig_bdf_file):
        """
        Setup pytacs object for a given input BDF file
        """
        # Create a temporary file on the root proc to hold our intermediate bdf
        if comm.rank == 0:
            temp_file = tempfile.NamedTemporaryFile(suffix=".bdf", delete=False)
            new_bdf_file = temp_file.name
            scaled_bdf_file = new_bdf_file.replace(".bdf", "_scaled.bdf")
            temp_file.close()
        else:
            new_bdf_file = None
            scaled_bdf_file = None

        # Broadcast the temp file name to other procs
        new_bdf_file = comm.bcast(new_bdf_file, root=0)
        scaled_bdf_file = comm.bcast(scaled_bdf_file, root=0)

        # Instantiate FEA Assembler
        fea_assembler = pytacs.pyTACS(orig_bdf_file, comm)
        # Set up constitutive objects and elements
        fea_assembler.initialize()
        # Read in forces from BDF and create tacs struct problems
        tacs_probs = fea_assembler.createTACSProbsFromBDF()
        # Convert from dict to list
        tacs_probs = list(tacs_probs.values())
        # We have to solve each problem before writing bdf file
        for problem in tacs_probs:
            problem.solve()

        # Store the assemblers and file paths for scaling tests
        self.fea_assembler = fea_assembler
        self.tacs_probs = tacs_probs
        self.new_bdf_file = new_bdf_file
        self.scaled_bdf_file = scaled_bdf_file

    def test_bdf_scaling(self):
        """Test scaling for several BDF files"""
        for test_case in self.TEST_CASES:
            print(f"\nRunning subTest for BDF: {test_case['name']}")
            with self.subTest(bdf=test_case["name"]):
                self._run_scaling_test(
                    test_case["bdf_file"],
                    test_case.get("xyz_scale", 1.0),
                    test_case.get("mass_scale", 1.0),
                    test_case.get("time_scale", 1.0),
                )

    def _run_scaling_test(
        self,
        orig_bdf_file,
        xyz_scale,
        mass_scale,
        time_scale,
    ):
        """Run scaling test for a single BDF file"""
        # Read in the original BDF, setup tacs and solver problems
        self.setup_tacs(self.comm, orig_bdf_file)

        # Write out a BDF file
        self.fea_assembler.writeBDF(self.new_bdf_file, self.tacs_probs)

        # Write a scaled BDF
        self.fea_assembler.writeBDF(
            self.scaled_bdf_file,
            self.tacs_probs,
            xyz_scale=xyz_scale,
            mass_scale=mass_scale,
            time_scale=time_scale,
        )

        # Only run the remaining comparison on root proc
        if self.comm.rank != 0:
            return

        # Read back the original and scaled BDF with pynastran
        new_bdf_info = BDF(debug=False)
        new_bdf_info.read_bdf(self.new_bdf_file)

        scaled_bdf_info = BDF(debug=False)
        scaled_bdf_info.read_bdf(self.scaled_bdf_file)

        # Run all comparisons
        self._compare_nodes(new_bdf_info, scaled_bdf_info, xyz_scale)
        self._compare_loads(
            new_bdf_info,
            scaled_bdf_info,
            xyz_scale,
            mass_scale=mass_scale,
            time_scale=time_scale,
        )
        self._compare_properties(new_bdf_info, scaled_bdf_info, xyz_scale)
        self._compare_materials(new_bdf_info, scaled_bdf_info, xyz_scale, mass_scale)

    def _compare_nodes(self, new_bdf, scaled_bdf, xyz_scale):
        """Compare nodal coordinates"""
        for node_id in new_bdf.nodes:
            new_node = new_bdf.nodes[node_id]
            scaled_node = scaled_bdf.nodes[node_id]

            new_xyz = np.array(new_node.xyz)
            scaled_xyz = np.array(scaled_node.xyz)

            # Verify scaling was applied correctly
            np.testing.assert_allclose(
                scaled_xyz,
                new_xyz * xyz_scale,
                rtol=1e-10,
                err_msg=f"Node {node_id} coordinates not scaled correctly",
            )

    def _compare_loads(self, new_bdf, scaled_bdf, xyz_scale, mass_scale, time_scale):
        """Compare forces and moments"""
        for load_id in new_bdf.loads:
            new_loads = new_bdf.loads[load_id]
            scaled_loads = scaled_bdf.loads[load_id]

            for new_load, scaled_load in zip(new_loads, scaled_loads):
                if new_load.type == "FORCE":
                    force_scale = mass_scale * xyz_scale / time_scale**2
                    # Forces should scale as mass * length / time^2
                    expected = new_load.mag * force_scale
                    np.testing.assert_allclose(
                        scaled_load.mag,
                        expected,
                        rtol=1e-10,
                        err_msg=f"Force load {new_load} magnitude not scaled correctly",
                    )
                elif new_load.type == "MOMENT":
                    # Moments should scale as force * length
                    moment_scale = (mass_scale * xyz_scale / time_scale**2) * xyz_scale
                    expected = new_load.mag * moment_scale
                    np.testing.assert_allclose(
                        scaled_load.mag,
                        expected,
                        rtol=1e-10,
                        err_msg=f"Moment load {new_load} magnitude not scaled correctly",
                    )

    def _compare_properties(self, new_bdf, scaled_bdf, xyz_scale):
        """Compare shell thicknesses and composite ply thicknesses"""
        for prop_id in new_bdf.properties:
            new_prop = new_bdf.properties[prop_id]
            scaled_prop = scaled_bdf.properties[prop_id]

            if new_prop.type == "PSHELL":
                if new_prop.t is not None:
                    np.testing.assert_allclose(
                        scaled_prop.t,
                        new_prop.t * xyz_scale,
                        rtol=1e-10,
                        err_msg=f"PSHELL {prop_id} thickness not scaled correctly",
                    )
            elif new_prop.type in ["PCOMP", "PCOMPG"]:
                # Check ply thicknesses
                for i, (new_t, scaled_t) in enumerate(
                    zip(new_prop.thicknesses, scaled_prop.thicknesses)
                ):
                    np.testing.assert_allclose(
                        scaled_t,
                        new_t * xyz_scale,
                        rtol=1e-10,
                        err_msg=f"PCOMP {prop_id} ply {i} thickness not scaled correctly",
                    )

    def _compare_materials(self, new_bdf, scaled_bdf, xyz_scale, mass_scale):
        """
        Compare material densities.

        Density has units of mass/volume = mass/length^3
        So when scaling: rho_scaled = rho * (mass_scale / xyz_scale^3)
        """
        for mat_id in new_bdf.materials:
            new_mat = new_bdf.materials[mat_id]
            scaled_mat = scaled_bdf.materials[mat_id]

            if hasattr(new_mat, "rho") and new_mat.rho is not None:
                # Density scales as: mass / length^3
                expected_rho = new_mat.rho * mass_scale / (xyz_scale**3)
                np.testing.assert_allclose(
                    scaled_mat.rho,
                    expected_rho,
                    rtol=1e-10,
                    err_msg=f"Material {mat_id} density not scaled correctly. "
                    f"Expected {expected_rho}, got {scaled_mat.rho}",
                )

    def tearDown(self):
        """Clean up temporary files"""
        if self.comm.rank == 0:
            if hasattr(self, "new_bdf_file") and os.path.exists(self.new_bdf_file):
                os.remove(self.new_bdf_file)

            if hasattr(self, "scaled_bdf_file") and os.path.exists(
                self.scaled_bdf_file
            ):
                os.remove(self.scaled_bdf_file)


if __name__ == "__main__":
    unittest.main()
