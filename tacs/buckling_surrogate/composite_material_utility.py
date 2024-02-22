__all__ = ["CompositeMaterialUtility"]

import numpy as np


class CompositeMaterialUtility:
    """
    utility class for computing composite material properties at an orientation
    Ref: Fiber Reinforced Composites TB by P.K. Mallick
    """

    @classmethod
    def from_fiber_matrix(cls, Ef, Em, nuf, num, vf):
        Gf = Ef / 2.0 / (1 + nuf)
        Gm = Em / 2.0 / (1 + num)
        vm = 1.0 - vf

        E11 = Ef * vf + Em * vm
        E22 = Ef * Em / (Ef * vm + Em * vf)
        nu12 = nuf * vf + num * vm
        G12 = Gf * Gm / (Gf * vm + Gm * nuf)
        return cls(E11=E11, E22=E22, nu12=nu12, G12=G12)

    def __init__(self, E11, E22, nu12, G12):
        self.E11 = E11
        self.E22 = E22
        self.nu12 = nu12
        self.G12 = G12

    def rotate_ply(self, angle=0.0):
        """
        compute the properties at a general angle
        assume the input angle is in degrees
        """
        # copy previous properties at 0 deg into temp variables
        E11 = self.E11 * 1.0
        E22 = self.E22 * 1.0
        nu12 = self.nu12 * 1.0
        G12 = self.G12 * 1.0

        # perform rotation to new axis system
        angle_rad = np.deg2rad(angle)
        C = np.cos(angle_rad)
        S = np.sin(angle_rad)
        C2 = np.cos(2 * angle_rad)
        S2 = np.sin(2 * angle_rad)
        self.E11 = (
            C**4 / E11 + S**4 / E22 + 0.25 * (1.0 / G12 - 2 * nu12 / E11) * S2**2
        ) ** (-1)
        self.E22 = (
            S**4 / E11 + C**4 / E22 + 0.25 * (1.0 / G12 - 2 * nu12 / E11) * S2**2
        ) ** (-1)
        _temp1 = 1.0 / E11 + 2.0 * nu12 / E11 + 1.0 / E22
        self.G12 = (_temp1 - (_temp1 - 1.0 / G12) * C2**2) ** (-1)
        self.nu12 = self.E11 * (nu12 / E11 - 0.25 * (_temp1 - 1.0 / G12) * S2**2)
        return self

    def __str__(self):
        return (
            f"E11 = {self.E11}, E22 = {self.E22}, nu12 = {self.nu12}, G12 = {self.G12}"
        )
