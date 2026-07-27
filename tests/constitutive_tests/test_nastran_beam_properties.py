"""
Fast, numbers-only unit tests for the pure Nastran beam-card translation helpers in
``tacs.nastran.properties``.

These tests exercise the helper functions directly with plain Python/numpy inputs, so they
require no BDF file and no MPI. The end-to-end Nastran-to-TACS beam translation (building a
BDF, running pyTACS, and comparing against Nastran reference results) is covered separately by
``tests/integration_tests/test_nastran_beam_compat.py``.
"""

import unittest

import numpy as np

from tacs.nastran.properties import (
    averageStationProps,
    cowperHollowCircleShearFactor,
    isoTubeBeamDims,
    shearCentreOffset,
)


class CowperHollowCircleShearFactorTest(unittest.TestCase):
    def test_solidCircleLimit(self):
        """Assert that `cowperHollowCircleShearFactor(0.0, nu)` equals the solid-circle
        correction factor 6(1+nu)/(7+6nu).
        """
        nu = 0.33
        expected = 6.0 * (1.0 + nu) / (7.0 + 6.0 * nu)
        self.assertAlmostEqual(cowperHollowCircleShearFactor(0.0, nu), expected)

    def test_thinWallLimit(self):
        """Assert that `cowperHollowCircleShearFactor(1.0, nu)` equals the thin-wall
        correction factor 2(1+nu)/(4+3nu).
        """
        nu = 0.33
        expected = 2.0 * (1.0 + nu) / (4.0 + 3.0 * nu)
        self.assertAlmostEqual(cowperHollowCircleShearFactor(1.0, nu), expected)

    def test_midValueBetweenLimits(self):
        """
        Given a diameter ratio strictly between the solid-circle (m=0) and thin-wall (m=1)
        cases,
        when cowperHollowCircleShearFactor is evaluated at m=0.5,
        then the result lies strictly between the two limiting values.
        """
        nu = 0.33
        solidLimit = cowperHollowCircleShearFactor(0.0, nu)
        thinWallLimit = cowperHollowCircleShearFactor(1.0, nu)
        midValue = cowperHollowCircleShearFactor(0.5, nu)
        self.assertGreater(midValue, min(solidLimit, thinWallLimit))
        self.assertLess(midValue, max(solidLimit, thinWallLimit))


class IsoTubeBeamDimsTest(unittest.TestCase):
    def test_rod(self):
        """Assert that a ROD's single radius dim maps to zero inner diameter and a wall
        thickness equal to the radius.
        """
        radius = 0.05
        innerDiameter, wallThickness = isoTubeBeamDims("ROD", [radius])
        self.assertAlmostEqual(innerDiameter, 0.0)
        self.assertAlmostEqual(wallThickness, radius)

    def test_tube(self):
        """Assert that TUBE dims [r_outer, r_inner] map to innerDiameter=2*r_inner and
        wallThickness=r_outer-r_inner.
        """
        outerRadius = 0.05
        innerRadius = 0.03
        innerDiameter, wallThickness = isoTubeBeamDims(
            "TUBE", [outerRadius, innerRadius]
        )
        self.assertAlmostEqual(innerDiameter, 2.0 * innerRadius)
        self.assertAlmostEqual(wallThickness, outerRadius - innerRadius)

    def test_tube2(self):
        """Assert that TUBE2 dims [r_outer, wall] map to wallThickness=wall and
        innerDiameter=2*(r_outer-wall).
        """
        outerRadius = 0.05
        wall = 0.02
        innerDiameter, wallThickness = isoTubeBeamDims("TUBE2", [outerRadius, wall])
        self.assertAlmostEqual(wallThickness, wall)
        self.assertAlmostEqual(innerDiameter, 2.0 * (outerRadius - wall))

    def test_nonCircularSectionRaisesValueError(self):
        """Assert that a non-circular section type (e.g. BAR) raises a ValueError."""
        with self.assertRaises(ValueError):
            isoTubeBeamDims("BAR", [0.05, 0.02])

    def test_tubeAndTube2AgreeOnSamePhysicalTube(self):
        """
        Given the same physical hollow tube described once via TUBE dims and once via TUBE2
        dims,
        when isoTubeBeamDims is called with each section type,
        then both calls return the same (innerDiameter, wallThickness) pair.
        """
        outerRadius = 0.05
        innerRadius = 0.03
        wall = outerRadius - innerRadius

        tubeDims = isoTubeBeamDims("TUBE", [outerRadius, innerRadius])
        tube2Dims = isoTubeBeamDims("TUBE2", [outerRadius, wall])

        self.assertAlmostEqual(tubeDims[0], tube2Dims[0])
        self.assertAlmostEqual(tubeDims[1], tube2Dims[1])


class AverageStationPropsTest(unittest.TestCase):
    def test_singleStationReturnsSoleValue(self):
        """Assert that with a single station (xxb=[0.0]) each property collapses to its
        sole array element.
        """
        result = averageStationProps({"A": [3.0]}, [0.0])
        self.assertEqual(result, {"A": 3.0})

    def test_twoStationUniformAveragesToConstant(self):
        """Assert that a two-station property that is constant across stations averages to
        that same constant value.
        """
        result = averageStationProps({"A": [2.0, 2.0]}, [0.0, 1.0])
        self.assertAlmostEqual(result["A"], 2.0)

    def test_twoStationLinearTaperAveragesToArithmeticMean(self):
        """
        Given a property that varies linearly between two stations,
        when averageStationProps integrates it via trapezoid over xxb in [0, 1],
        then the result equals the simple arithmetic mean (a+b)/2, because xxb spans a unit
        interval so the trapezoid integral coincides with the length-average. Multiple keys
        are checked at once to confirm every entry in the dict is averaged this way.
        """
        result = averageStationProps({"A": [1.0, 3.0], "B": [4.0, 6.0]}, [0.0, 1.0])
        self.assertAlmostEqual(result["A"], (1.0 + 3.0) / 2.0)
        self.assertAlmostEqual(result["B"], (4.0 + 6.0) / 2.0)


class ShearCentreOffsetTest(unittest.TestCase):
    def _buildCbar(self, wa=None, wb=None):
        """Build a minimal in-memory pyNastran BDF with 2 grids, a MAT1, a PBAR, and a
        single CBAR, optionally with WA/WB offsets, and return the cross-referenced CBAR
        element card together with the BDF object.
        """
        from pyNastran.bdf.bdf import BDF

        bdfInfo = BDF()
        bdfInfo.add_grid(1, [0.0, 0.0, 0.0])
        bdfInfo.add_grid(2, [1.0, 0.0, 0.0])
        bdfInfo.add_mat1(1, 7.0e10, None, 0.33)
        bdfInfo.add_pbar(1, 1, A=1.0, i1=1.0, i2=1.0, i12=0.0, j=1.0)
        bdfInfo.add_cbar(1, 1, [1, 2], x=[0.0, 1.0, 0.0], g0=None, wa=wa, wb=wb)
        bdfInfo.cross_reference()
        elem0 = bdfInfo.elements[1]
        return elem0, bdfInfo

    def test_zeroOffsetGivesNoShearCentreOffset(self):
        """Assert that a CBAR with no WA/WB offset yields zero projected offsets and
        hasShearCenterOffset == False.
        """
        elem0, bdfInfo = self._buildCbar()
        (
            shearCenterYOffset,
            shearCenterZOffset,
            hasShearCenterOffset,
            _yElem,
            _zElem,
            _offsetVector,
        ) = shearCentreOffset(elem0, bdfInfo)
        self.assertAlmostEqual(shearCenterYOffset, 0.0)
        self.assertAlmostEqual(shearCenterZOffset, 0.0)
        self.assertFalse(hasShearCenterOffset)

    def test_knownOffsetProjectsOntoSectionAxes(self):
        """
        Given a CBAR with identical WA/WB offset vectors,
        when shearCentreOffset projects the averaged offset onto the element's local y/z
        section axes,
        then it returns those known projected components and hasShearCenterOffset == True.
        """
        offset = [0.0, 0.05, 0.02]
        elem0, bdfInfo = self._buildCbar(wa=offset, wb=offset)
        (
            shearCenterYOffset,
            shearCenterZOffset,
            hasShearCenterOffset,
            yElem,
            zElem,
            offsetVector,
        ) = shearCentreOffset(elem0, bdfInfo)
        self.assertAlmostEqual(shearCenterYOffset, np.dot(yElem, offset))
        self.assertAlmostEqual(shearCenterZOffset, np.dot(zElem, offset))
        np.testing.assert_allclose(offsetVector, offset)
        self.assertTrue(hasShearCenterOffset)


if __name__ == "__main__":
    unittest.main()
