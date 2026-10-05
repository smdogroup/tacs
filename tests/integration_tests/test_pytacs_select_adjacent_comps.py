"""
Tests for the ``adjacentTo`` argument of pyTACS.selectCompIDs.

Two meshes are checked:

- partitioned_plate.bdf: a flat shell plate split into 4 quadrant components (PLATE.00 - PLATE.03),
  ordered as below, which are adjacent when they share an element edge.

               |
        3      |       2
               |
    -----------|------------
               |
       0       |      1
               |

- multi_beam_frame.bdf: a frame of 4 beam components ("Property group 1" - "Property group 4"),
  which are adjacent when they share an end node. Component i is "Property group i+1".

The tests run on multiple procs to check that every proc gets the same answer.
"""

import os
import unittest

import numpy as np
from mpi4py import MPI

from tacs import pytacs

base_dir = os.path.dirname(os.path.abspath(__file__))
PLATE_BDF_FILE = os.path.join(base_dir, "./input_files/partitioned_plate.bdf")
FRAME_BDF_FILE = os.path.join(base_dir, "./input_files/multi_beam_frame.bdf")

PLATE_ADJACENT_COMP_IDS = {0: [1, 3], 1: [0, 2], 2: [1, 3], 3: [0, 2]}
FRAME_ADJACENT_COMP_IDS = {0: [1, 3], 1: [0, 2, 3], 2: [1], 3: [0, 1]}


class SelectAdjacentCompIDsTest(unittest.TestCase):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    def setUp(self):
        self.comm = MPI.COMM_WORLD

    def selectCompIDs(self, bdfFile, **kwargs) -> list[int]:
        """Create a pyTACS object for the bdf and return the sorted, single-group selection."""
        feaAssembler = pytacs.pyTACS(bdfFile, self.comm)
        compIDs = feaAssembler.selectCompIDs(**kwargs)
        self.assertEqual(len(compIDs), 1)
        return sorted(compIDs[0])

    def test_adjacent_to_single_shell_comp(self):
        for compID, expected in PLATE_ADJACENT_COMP_IDS.items():
            with self.subTest(compID=compID):
                # Select by component ID
                self.assertEqual(
                    self.selectCompIDs(PLATE_BDF_FILE, adjacentTo=compID), expected
                )
                # Select by (case-insensitive) component name
                self.assertEqual(
                    self.selectCompIDs(
                        PLATE_BDF_FILE, adjacentTo=f"plate.{compID:02d}"
                    ),
                    expected,
                )

    def test_adjacent_to_single_beam_comp(self):
        for compID, expected in FRAME_ADJACENT_COMP_IDS.items():
            with self.subTest(compID=compID):
                self.assertEqual(
                    self.selectCompIDs(
                        FRAME_BDF_FILE, adjacentTo=f"Property group {compID + 1}"
                    ),
                    expected,
                )

    def test_adjacent_to_multiple_comps(self):
        # Neighbours of either comp 2 or comp 3 -> union of {1} and {0, 1}
        self.assertEqual(
            self.selectCompIDs(FRAME_BDF_FILE, adjacentTo=["Property group 3", 3]),
            [0, 1],
        )

    def test_adjacent_to_combined_with_include(self):
        # Only PLATE.01 is both named PLATE.01 and adjacent to PLATE.00
        self.assertEqual(
            self.selectCompIDs(
                PLATE_BDF_FILE, include="PLATE.01", adjacentTo="PLATE.00"
            ),
            [1],
        )
        # PLATE.02 is not adjacent to PLATE.00, so nothing is selected
        self.assertEqual(
            self.selectCompIDs(
                PLATE_BDF_FILE, include="PLATE.02", adjacentTo="PLATE.00"
            ),
            [],
        )

    def test_adjacent_to_combined_with_exclude(self):
        self.assertEqual(
            self.selectCompIDs(
                PLATE_BDF_FILE, adjacentTo="PLATE.00", exclude="PLATE.03"
            ),
            [1],
        )

    def test_adjacent_to_does_not_include_reference_comps(self):
        # 'PLATE' matches every component, and every component has a neighbour, so this selects
        # everything. The reference components are only kept because they neighbour each other.
        self.assertEqual(
            self.selectCompIDs(PLATE_BDF_FILE, adjacentTo="PLATE"), [0, 1, 2, 3]
        )
        # Comp 2 only neighbours comp 1, so comp 2 itself is not selected
        self.assertEqual(
            self.selectCompIDs(FRAME_BDF_FILE, adjacentTo="Property group 3"), [1]
        )

    def test_adjacent_to_with_groups(self):
        feaAssembler = pytacs.pyTACS(FRAME_BDF_FILE, self.comm)
        compIDs = feaAssembler.selectCompIDs(adjacentTo="Property group 2", nGroup=-1)
        self.assertEqual(sorted(compIDs), [[0], [2], [3]])

    def test_adjacent_to_accepts_selected_comp_ids(self):
        # The output of one selectCompIDs call can be fed back into adjacentTo / include
        feaAssembler = pytacs.pyTACS(PLATE_BDF_FILE, self.comm)
        allCompIDs = feaAssembler.selectCompIDs()[0]
        self.assertTrue(all(isinstance(compID, int) for compID in allCompIDs))
        refCompIDs = feaAssembler.selectCompIDs(include="PLATE.00")[0]
        self.assertEqual(
            sorted(feaAssembler.selectCompIDs(adjacentTo=refCompIDs)[0]), [1, 3]
        )
        self.assertEqual(
            sorted(feaAssembler.selectCompIDs(include=allCompIDs)[0]), [0, 1, 2, 3]
        )
        # numpy integers are also accepted as component IDs
        self.assertEqual(
            sorted(
                feaAssembler.selectCompIDs(
                    include=np.array([0, 1, 2]), adjacentTo=np.int64(0)
                )[0]
            ),
            [1],
        )

    def test_adjacent_to_none_has_no_effect(self):
        self.assertEqual(
            self.selectCompIDs(PLATE_BDF_FILE, include="PLATE", adjacentTo=None),
            [0, 1, 2, 3],
        )


if __name__ == "__main__":
    unittest.main()
