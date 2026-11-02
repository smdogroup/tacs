"""
Check that design variable values reach the auxiliary load elements.

A single shell element with a thickness DV (DV 0) carries an auxiliary pressure
element whose pressure is DV 1. The tests set a new pressure in the design
vector and check the value the pressure element then holds, both for the default
design vector (one value per DV) and for a design vector with two values per DV.
"""

import unittest

import numpy as np
from mpi4py import MPI

from tacs import TACS, constitutive, elements

DTYPE = TACS.dtype
THICKNESS_DV_NUM = 0
PRESSURE_DV_NUM = 1
NUM_SHELL_NODES = 4
SHELL_XPTS = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0]
INITIAL_THICKNESS = 0.01
INITIAL_PRESSURE = 1.0e3
NEW_PRESSURE = 2.5e5
DVS_PER_NODE = 2


def createAssembler(designVarsPerNode=1):
    """Return a one-shell assembler and the aux pressure element attached to it.

    Parameters
    ----------
    designVarsPerNode : int
        Number of values per DV in the design vector, set through Assembler.setDesignNodeMap.

    Returns
    -------
    tuple
        The initialized assembler, the AuxElements object and the pressure element.
    """
    prop = constitutive.MaterialProperties(rho=2700.0, E=70e9, nu=0.3, ys=270e6)
    con = constitutive.IsoShellConstitutive(
        prop, t=INITIAL_THICKNESS, tNum=THICKNESS_DV_NUM
    )
    shell = elements.Quad4Shell(elements.ShellNaturalTransform(), con)

    assembler = TACS.Assembler.create(
        MPI.COMM_SELF, shell.getVarsPerNode(), NUM_SHELL_NODES, 1
    )
    assembler.setElementConnectivity(
        np.array([0, NUM_SHELL_NODES], dtype=np.intc),
        np.arange(NUM_SHELL_NODES, dtype=np.intc),
    )
    assembler.setElements([shell])

    pressureElem = shell.createElementPressure(0, INITIAL_PRESSURE, PRESSURE_DV_NUM)
    auxElems = TACS.AuxElements()
    auxElems.addElement(0, pressureElem)
    # Attach before initialize() so the pressure DV is counted in the design vector
    assembler.setAuxElements(auxElems)
    if designVarsPerNode != 1:
        assembler.setDesignNodeMap(designVarsPerNode, None)
    assembler.initialize()

    nodes = assembler.createNodeVec()
    nodes.getArray()[:] = np.array(SHELL_XPTS, dtype=DTYPE)
    assembler.setNodes(nodes)
    return assembler, auxElems, pressureElem


class AuxElementSetDesignVarsTest(unittest.TestCase):
    """Design variable updates of auxiliary load elements."""

    def test_assembler_set_design_vars_updates_aux_element(self):
        """
        Given an assembler with an aux pressure element whose pressure is a DV,
        when Assembler.setDesignVars is called with a new pressure in the design vector,
        then the pressure element holds the new pressure.
        """
        assembler, _, pressureElem = createAssembler()
        x = assembler.createDesignVec()
        assembler.getDesignVars(x)
        x.getArray()[PRESSURE_DV_NUM] = NEW_PRESSURE

        assembler.setDesignVars(x)

        np.testing.assert_array_equal(pressureElem.getDesignVars(0), [NEW_PRESSURE])

    def test_block_size_two_design_vector(self):
        """
        Given an assembler with two design values per DV and an aux pressure element whose pressure is DV 1,
        when AuxElements.setDesignVars or Assembler.setDesignVars is called with two different values in block 1,
        then the pressure element holds the first value of block 1.
        """
        # The heap overflow that this test guards against only fails under
        # AddressSanitizer; a real build can pass even with the overflow
        updateFunctions = {
            "AuxElements.setDesignVars": lambda assembler, auxElems, x: (
                auxElems.setDesignVars(x)
            ),
            "Assembler.setDesignVars": lambda assembler, auxElems, x: (
                assembler.setDesignVars(x)
            ),
        }
        for name, updateDesignVars in updateFunctions.items():
            with self.subTest(update=name):
                assembler, auxElems, pressureElem = createAssembler(DVS_PER_NODE)
                x = assembler.createDesignVec()
                xArray = x.getArray()
                self.assertEqual(len(xArray), DVS_PER_NODE * (PRESSURE_DV_NUM + 1))
                xArray[:] = INITIAL_THICKNESS
                pressureBlockStart = DVS_PER_NODE * PRESSURE_DV_NUM
                xArray[pressureBlockStart] = NEW_PRESSURE
                # A different second value, so that reading the wrong entry of the block is detected
                xArray[pressureBlockStart + 1] = -NEW_PRESSURE

                updateDesignVars(assembler, auxElems, x)

                np.testing.assert_array_equal(
                    pressureElem.getDesignVars(0), [NEW_PRESSURE]
                )


if __name__ == "__main__":
    unittest.main()
