import unittest

import numpy as np
from mpi4py import MPI

from tacs import TACS, constitutive, elements

"""
Regression test for TACSBeamElement::getAverageStresses.

getAverageStresses is a forward-only value computation exposed through
Assembler.getAverageStresses (not a sensitivity), so the TacsTestElement*
FD/CS harnesses do not apply. This test builds a single straight 2-node
beam element directly via the low-level TACSAssembler interface (a single
element does not need mesh partitioning), prescribes a known uniform
axial-strain state directly via the state vector (bypassing solving, since
only the forward stress computation is under test), and checks the
assembler-level average stress against the closed-form axial stress
resultant s[0] = E*A*eps.
"""


class BeamAverageStressesTest(unittest.TestCase):
    def setUp(self):
        if TACS.dtype is complex:
            self.rtol = 1e-10
            self.atol = 1e-10
        else:
            self.rtol = 1e-10
            self.atol = 1e-8

        self.comm = MPI.COMM_WORLD

        # Material properties
        self.rho = 2700.0
        self.E = 70e3
        self.nu = 0.3
        self.ys = 270.0
        self.props = constitutive.MaterialProperties(
            rho=self.rho, E=self.E, nu=self.nu, ys=self.ys
        )

        # Tube cross-section (inner diameter "d", wall thickness "t" --
        # matches TACSIsoTubeBeamConstitutive.cpp's own inner/wall naming)
        self.d = 1.0
        self.t = 0.1
        d0 = self.d + 2.0 * self.t
        d1 = self.d
        self.A = np.pi * (d0**2 - d1**2) / 4.0

        self.con = constitutive.IsoTubeBeamConstitutive(self.props, d=self.d, t=self.t)

        ref_axis = np.array([0.0, 1.0, 1.0], dtype=TACS.dtype)
        self.transform = elements.BeamRefAxisTransform(ref_axis)

        # Length of the beam (arbitrary, along global x)
        self.L = 2.0
        # Prescribed uniform axial strain
        self.eps = 0.01

    def _build_assembler(self, element):
        """
        Build a minimal single-element, 2-node assembler directly via the
        TACSAssembler interface (no Creator/BDF -- a single element does
        not need mesh partitioning).
        """
        vars_per_node = element.getVarsPerNode()
        num_nodes = 2
        num_elements = 1

        assembler = TACS.Assembler.create(
            self.comm, vars_per_node, num_nodes, num_elements
        )

        ptr = np.array([0, num_nodes], dtype=np.intc)
        conn = np.array([0, 1], dtype=np.intc)
        assembler.setElementConnectivity(ptr, conn)
        assembler.setElements([element])
        assembler.initialize()

        # Straight beam along global x, from the origin to (L, 0, 0)
        Xvec = assembler.createNodeVec()
        Xarr = Xvec.getArray()
        Xarr[:] = 0.0
        Xarr[0:3] = [0.0, 0.0, 0.0]
        Xarr[3:6] = [self.L, 0.0, 0.0]
        assembler.setNodes(Xvec)

        return assembler, vars_per_node

    def _set_uniform_axial_state(self, assembler, vars_per_node):
        """
        Prescribe u_x(node) = eps * X(node), all other DOFs zero -- a
        linear axial displacement ramp along a straight beam gives a
        uniform axial strain e[0] = eps exactly (independent of the
        transverse reference-axis choice, which only resolves the beam's
        transverse frame, not its tangent direction).
        """
        varsVec = assembler.createVec()
        varsArr = varsVec.getArray()
        varsArr[:] = 0.0
        varsArr[0] = self.eps * 0.0
        varsArr[vars_per_node] = self.eps * self.L
        assembler.setVariables(varsVec)
        return varsVec

    def test_beam2_average_stresses_uniform_axial(self):
        element = elements.Beam2(self.transform, self.con)
        assembler, vars_per_node = self._build_assembler(element)
        self._set_uniform_axial_state(assembler, vars_per_node)

        avg = assembler.getAverageStresses(0)

        expected_axial = self.E * self.A * self.eps
        self.assertAlmostEqual(
            complex(avg[0]).real,
            expected_axial,
            delta=self.atol + self.rtol * abs(expected_axial),
        )
        # No bending/torsion/shear should appear for a pure axial state.
        for i in range(1, 6):
            self.assertAlmostEqual(complex(avg[i]).real, 0.0, delta=self.atol)
        # The beam's stress vector has only 6 components; avgStresses[6:9]
        # (shell-only drill/extra slots) are left untouched. Since
        # Assembler.getAverageStresses allocates a fresh zero buffer, this
        # assertion cannot distinguish "left untouched" from "explicitly
        # zeroed", but the implementation never writes to indices 6..8.
        for i in range(6, 9):
            self.assertEqual(complex(avg[i]).real, 0.0)

    def test_beam3_average_stresses_uniform_axial(self):
        element = elements.Beam3(self.transform, self.con)
        vars_per_node = element.getVarsPerNode()
        num_nodes = 3
        num_elements = 1

        assembler = TACS.Assembler.create(
            self.comm, vars_per_node, num_nodes, num_elements
        )
        ptr = np.array([0, num_nodes], dtype=np.intc)
        conn = np.array([0, 1, 2], dtype=np.intc)
        assembler.setElementConnectivity(ptr, conn)
        assembler.setElements([element])
        assembler.initialize()

        Xvec = assembler.createNodeVec()
        Xarr = Xvec.getArray()
        Xarr[:] = 0.0
        # Beam3's mid-node ordering follows the node ordering convention
        # used elsewhere in this test file's sibling element tests
        # (getNodePoint's parametric location -1, 1, 0 for a 3-node line
        # element) -- node 2 is the mid-node at x = L/2.
        Xarr[0:3] = [0.0, 0.0, 0.0]
        Xarr[3:6] = [self.L, 0.0, 0.0]
        Xarr[6:9] = [0.5 * self.L, 0.0, 0.0]
        assembler.setNodes(Xvec)

        varsVec = assembler.createVec()
        varsArr = varsVec.getArray()
        varsArr[:] = 0.0
        varsArr[0 * vars_per_node] = self.eps * 0.0
        varsArr[1 * vars_per_node] = self.eps * self.L
        varsArr[2 * vars_per_node] = self.eps * 0.5 * self.L
        assembler.setVariables(varsVec)

        avg = assembler.getAverageStresses(0)

        expected_axial = self.E * self.A * self.eps
        self.assertAlmostEqual(
            complex(avg[0]).real,
            expected_axial,
            delta=self.atol + self.rtol * abs(expected_axial),
        )
        for i in range(1, 6):
            self.assertAlmostEqual(complex(avg[i]).real, 0.0, delta=self.atol)
        for i in range(6, 9):
            self.assertEqual(complex(avg[i]).real, 0.0)


if __name__ == "__main__":
    unittest.main()
