import os

import tempfile

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs

"""
A unit cube with each face constructed from a quad shell element is tested. 
Two of the faces have flipped normals (pointing inwards).
This test verifies the robustness of the VolumeConstraint class for shell elements.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/unit_cube.bdf")


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "volume_con": 1.0,
    }

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        """

        # Instantiate FEA Assembler
        fea_assembler = pytacs.pyTACS(bdf_file, comm)

        # Set up constitutive objects and elements
        fea_assembler.initialize()

        # Create tacs volume constraint
        constr = fea_assembler.createVolumeConstraint("volume")
        allIDs = fea_assembler.selectCompIDs()
        constr.addConstraint("con", compIDs=allIDs)

        return [constr], fea_assembler

    def test_is_linear(self):
        """
        Test that the volume constraint is linear wrt nodes/dvs.
        """
        constr = self.tacs_probs[0]
        self.assertFalse(
            constr.isLinear, "Volume constraint should not be linear wrt dvs"
        )

    def test_write_visualization(self):
        """
        Test tecplot visualization writer.
        """
        # Create temporary directory to write f5 file to (only on root)
        tmp_dir = None
        tmp_dir_name = None
        if self.comm.rank == 0:
            tmp_dir = tempfile.TemporaryDirectory()
            tmp_dir_name = tmp_dir.name
        # Broadcast temp directory name to other procs
        tmp_dir_name = self.comm.bcast(tmp_dir_name, root=0)
        constr = self.tacs_probs[0]
        f = {}
        constr.evalConstraints(f)
        constr.writeVisualization(outputDir=tmp_dir_name)
        if self.comm.rank == 0:
            # delete all files in temp dir
            tmp_dir.cleanup()
