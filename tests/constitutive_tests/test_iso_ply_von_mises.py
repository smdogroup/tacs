import unittest

import numpy as np

from tacs import TACS, constitutive

"""
Verify that the Tsai-Wu failure criterion reduces to the von Mises criterion
for isotropic materials.

TACSOrthotropicPly documents that its default TSAI_WU_MODIFIED criterion
"becomes equivalent to the von-Mises criterion when the material properties are
isotropic". These tests pin that claim down by comparing the failure envelope of
a single-ply smeared composite shell against an isotropic shell built from the
same material. Both constitutive models aggregate failure over the same two
through-thickness points with the same KS weight, so their envelopes coincide
only if the underlying pointwise failure criteria agree.
"""

# Indices of the in-plane stress resultants in the shell stress vector
NXX_INDEX = 0
NYY_INDEX = 1
NXY_INDEX = 2

# Number of stress components in a shell constitutive model
NUM_SHELL_STRESSES = 9


class IsotropicPlyVonMisesTest(unittest.TestCase):
    def setUp(self):
        self.rho = 2700.0
        self.specific_heat = 921.096
        self.E = 70e9
        self.nu = 0.3
        self.ys = 270e6
        self.thickness = 0.01
        self.num_envelope_points = 25

        # Failure envelopes are in stress-resultant space (force per unit
        # length), so the natural scale is ys * thickness.
        self.rtol = 1e-6
        self.atol = 1e-6 * self.ys * self.thickness

        self.props = constitutive.MaterialProperties(
            rho=self.rho,
            specific_heat=self.specific_heat,
            E=self.E,
            nu=self.nu,
            ys=self.ys,
        )
        ply = constitutive.OrthotropicPly(self.thickness, self.props)

        self.iso_con = constitutive.IsoShellConstitutive(self.props, t=self.thickness)
        self.ply_con = constitutive.SmearedCompositeShellConstitutive(
            [ply],
            self.thickness,
            np.array([0.0], dtype=TACS.dtype),
            np.array([1.0], dtype=TACS.dtype),
        )

    def assert_envelopes_match(self, first_index, second_index):
        """
        Compare the failure envelopes of the two constitutive models in the
        plane spanned by two stress-resultant directions.

        Parameters
        ----------
        first_index : int
            Index of the stress resultant swept along the envelope's x axis.
        second_index : int
            Index of the stress resultant swept along the envelope's y axis.
        """
        sx = np.zeros(NUM_SHELL_STRESSES, dtype=TACS.dtype)
        sy = np.zeros(NUM_SHELL_STRESSES, dtype=TACS.dtype)
        sx[first_index] = 1.0
        sy[second_index] = 1.0

        iso_x, iso_y = self.iso_con.getFailureEnvelope(
            sx, sy, npts=self.num_envelope_points
        )
        ply_x, ply_y = self.ply_con.getFailureEnvelope(
            sx, sy, npts=self.num_envelope_points
        )

        np.testing.assert_allclose(
            np.real(ply_x), np.real(iso_x), rtol=self.rtol, atol=self.atol
        )
        np.testing.assert_allclose(
            np.real(ply_y), np.real(iso_y), rtol=self.rtol, atol=self.atol
        )

    def test_shear_strength_is_von_mises(self):
        """An isotropic material yields in pure shear at ys / sqrt(3)."""
        S12 = self.props.getMaterialProperties()["S12"]
        expected = self.ys / np.sqrt(3.0)
        self.assertAlmostEqual(np.real(S12) / expected, 1.0, places=12)

    def test_shear_envelope_matches_von_mises(self):
        """
        Sweep the (Nxx, Nxy) plane. Nyy = 0 there, so the F12 interaction term
        drops out and this isolates the F66 shear term.
        """
        self.assert_envelopes_match(NXX_INDEX, NXY_INDEX)

    def test_biaxial_envelope_matches_von_mises(self):
        """
        Sweep the (Nxx, Nyy) plane. Nxy = 0 there, so the F66 shear term drops
        out and this isolates the F12 interaction term.
        """
        self.assert_envelopes_match(NXX_INDEX, NYY_INDEX)


if __name__ == "__main__":
    unittest.main()
