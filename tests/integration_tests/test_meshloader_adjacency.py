"""
Tests for the component adjacency queries on pyMeshLoader.

Two meshes are checked:

- partitioned_plate.bdf: a flat shell plate split into 4 quadrant components, ordered as below,
  which are adjacent when they share an element edge.

               |
        3      |       2
               |
    -----------|------------
               |
       0       |      1
               |

- multi_beam_frame.bdf: a frame of 4 beam components, which are adjacent when they share an end node.

A small mesh is also built in memory with one component that is not connected to any other.

The tests run on multiple procs to check that every proc gets the same answer.
"""

import os
import unittest

import pyNastran.bdf.bdf as pn
from mpi4py import MPI

from tacs.pymeshloader import pyMeshLoader
from tacs.utilities import Error as TACSError

base_dir = os.path.dirname(os.path.abspath(__file__))
PLATE_BDF_FILE = os.path.join(base_dir, "./input_files/partitioned_plate.bdf")
FRAME_BDF_FILE = os.path.join(base_dir, "./input_files/multi_beam_frame.bdf")

PLATE_ADJACENT_COMP_IDS = {0: [1, 3], 1: [0, 2], 2: [1, 3], 3: [0, 2]}
FRAME_ADJACENT_COMP_IDS = {0: [1, 3], 1: [0, 2, 3], 2: [1], 3: [0, 1]}


def createIsolatedCompBDF() -> pn.BDF:
    """
    Create a shell mesh with 3 components, one CQUAD4 each.

    Components 0 and 1 share an edge, component 2 is not connected to either of them.
    """
    bdf = pn.BDF()
    bdf.add_mat1(1, 70e9, None, 0.3)
    for propID in [1, 2, 3]:
        bdf.add_pshell(propID, mid1=1, t=0.01, mid2=1, mid3=1)

    # Nodes 1-6 form a 2x1 strip of quads, nodes 7-10 form a separate quad
    nodeXYZ = {
        1: [0.0, 0.0, 0.0],
        2: [1.0, 0.0, 0.0],
        3: [2.0, 0.0, 0.0],
        4: [0.0, 1.0, 0.0],
        5: [1.0, 1.0, 0.0],
        6: [2.0, 1.0, 0.0],
        7: [5.0, 0.0, 0.0],
        8: [6.0, 0.0, 0.0],
        9: [6.0, 1.0, 0.0],
        10: [5.0, 1.0, 0.0],
    }
    for nodeID, xyz in nodeXYZ.items():
        bdf.add_grid(nodeID, xyz)

    bdf.add_cquad4(1, 1, [1, 2, 5, 4])
    bdf.add_cquad4(2, 2, [2, 3, 6, 5])
    bdf.add_cquad4(3, 3, [7, 8, 9, 10])
    return bdf


class MeshLoaderAdjacencyTest(unittest.TestCase):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    def setUp(self):
        self.comm = MPI.COMM_WORLD

    def loadMesh(self, bdf: str | pn.BDF) -> pyMeshLoader:
        """Create a pyMeshLoader and scan the given bdf file or object with it."""
        meshLoader = pyMeshLoader(self.comm)
        meshLoader.scanBdfFile(bdf)
        return meshLoader

    def assertAdjacentCompIDs(
        self, meshLoader: pyMeshLoader, expectedAdjacentCompIDs: dict[int, list[int]]
    ):
        """Check the adjacent component IDs of every component in the mesh."""
        self.assertEqual(meshLoader.getNumComponents(), len(expectedAdjacentCompIDs))
        for compID, expected in expectedAdjacentCompIDs.items():
            with self.subTest(compID=compID):
                self.assertEqual(meshLoader.getAdjacentCompIDs(compID), expected)

    def test_shell_adjacency(self):
        meshLoader = self.loadMesh(PLATE_BDF_FILE)
        self.assertAdjacentCompIDs(meshLoader, PLATE_ADJACENT_COMP_IDS)

    def test_beam_adjacency(self):
        meshLoader = self.loadMesh(FRAME_BDF_FILE)
        self.assertAdjacentCompIDs(meshLoader, FRAME_ADJACENT_COMP_IDS)

    def test_isolated_component(self):
        meshLoader = self.loadMesh(createIsolatedCompBDF())
        self.assertAdjacentCompIDs(meshLoader, {0: [1], 1: [0], 2: []})

    def test_rescan_resets_adjacency(self):
        meshLoader = self.loadMesh(PLATE_BDF_FILE)
        # Build the adjacency for the plate before scanning the frame
        meshLoader.getAdjacentCompIDs(0)
        meshLoader.scanBdfFile(FRAME_BDF_FILE)
        self.assertAdjacentCompIDs(meshLoader, FRAME_ADJACENT_COMP_IDS)

    def test_returned_list_is_a_copy(self):
        meshLoader = self.loadMesh(PLATE_BDF_FILE)
        adjacentCompIDs = meshLoader.getAdjacentCompIDs(0)
        adjacentCompIDs.append(2)
        self.assertEqual(meshLoader.getAdjacentCompIDs(0), PLATE_ADJACENT_COMP_IDS[0])

    def test_invalid_component_id(self):
        meshLoader = self.loadMesh(PLATE_BDF_FILE)
        numComps = meshLoader.getNumComponents()
        for badCompID in [-1, numComps]:
            with self.subTest(compID=badCompID):
                with self.assertRaises(TACSError):
                    meshLoader.getAdjacentCompIDs(badCompID)


if __name__ == "__main__":
    unittest.main()
