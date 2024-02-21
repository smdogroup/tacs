"""
Sean Engelstad, Feb 2024
GT SMDO Lab
"""
from tacs import buckling_surrogate
import numpy as np
import unittest

# 3 main verification cases for the buckling analysis
# 1, 2, 3
case = 3
run_static = False


class TestPlateCases(unittest.TestCase):
    # use the same flat plate geometry and material through all three test cases
    flat_plate = buckling_surrogate.FlatPlateAnalysis(
        bdf_file="plate.bdf",
        a=1.0,
        b=0.7,
        h=0.07,
        E11=70e9,
        nu12=0.33,
        E22=None,  # set to None if isotropic
        G12=None,  # set to None if isotropic
    )

    def test_uniaxial_compression(self):
        # case 1 - simply supported, uniaxial compression buckling

        self.flat_plate.generate_bdf(
            nx=30,
            ny=20,
            exx=0.001,
            eyy=0.0,
            exy=0.0,
            clamped=False,
        )

        # flat_plate.run_static_analysis(write_soln=True)

        tacs_eigvals = self.flat_plate.run_buckling_analysis(
            sigma=30.0, num_eig=12, write_soln=True
        )

        # eigenvalues from Abaqus for comparison
        abaqus_eigvals = np.array([36.083, 38.000, 51.634, 72.896, 96.711, 113.94])
        rel_error = (tacs_eigvals[:6] - abaqus_eigvals) / abaqus_eigvals

        print(
            f"\n\nVerification of case 1 - uniaxial compression,\n\tsimply supported plate buckling modes\n"
        )
        for i in range(6):
            print(f"mode {i+1} eigenvalues:")
            print(
                f"\ttacs = {tacs_eigvals[i]:.4f}, abaqus = {abaqus_eigvals[i]:.4f}, rel err = {rel_error[i]:.4f}"
            )

        assert np.max(np.abs(rel_error)) < 0.1

    def test_pure_shear_clamped(self):
        # case 2 - pure shear, clamped plate
        self.flat_plate.generate_bdf(
            nx=30,
            ny=20,
            exx=0.0,
            eyy=0.0,
            exy=0.001,
            clamped=True,
        )

        # flat_plate.run_static_analysis(write_soln=True)

        tacs_eigvals = self.flat_plate.run_buckling_analysis(
            sigma=30.0, num_eig=12, write_soln=True
        )

        # every other shear mode has negative eigenvalue (since reverse shear load will still cause failure by sym)
        pos_tacs_eigvals = [eigval for eigval in tacs_eigvals if eigval > 0.0]

        # eigenvalues from Abaqus for comparison
        abaqus_eigvals = np.array([111.79, 115.45, 169.71, 181.02, 236.06, 242.07])
        rel_error = (pos_tacs_eigvals[:6] - abaqus_eigvals) / abaqus_eigvals

        print(
            f"\n\nVerification of case 2 - pure shear,\n\tclamped plate buckling modes\n"
        )
        for i in range(6):
            print(f"mode {i+1} eigenvalues:")
            print(
                f"\ttacs = {pos_tacs_eigvals[i]:.4f}, abaqus = {abaqus_eigvals[i]:.4f}, rel err = {rel_error[i]:.4f}"
            )

        assert np.max(np.abs(rel_error)) < 0.1

    def test_combined_axial_shear_clamped(self):
        # case 3 - mixed shear and uniaxial compression case, clamped plate
        self.flat_plate.generate_bdf(
            nx=30,
            ny=20,
            exx=0.001,
            eyy=0.0,
            exy=0.001,
            clamped=True,
        )

        # flat_plate.run_static_analysis(write_soln=True)

        tacs_eigvals = self.flat_plate.run_buckling_analysis(
            sigma=30.0, num_eig=12, write_soln=True
        )

        # eigenvalues from Abaqus for comparison
        abaqus_eigvals = np.array([42.237, 42.290, 59.616, 68.999, 88.243, 91.976])
        rel_error = (tacs_eigvals[:6] - abaqus_eigvals) / abaqus_eigvals

        print(
            f"\n\nVerification of case 3 - mixed compression + shear,\n\tclamped plate buckling modes\n"
        )
        for i in range(6):
            print(f"mode {i+1} eigenvalues:")
            print(
                f"\ttacs = {tacs_eigvals[i]:.4f}, abaqus = {abaqus_eigvals[i]:.4f}, rel err = {rel_error[i]:.4f}"
            )

        assert np.max(np.abs(rel_error)) < 0.1


if __name__ == "__main__":
    test = "all"

    # run all cases with this
    if test == "all":
        unittest.main()

    elif test == 1:
        # run individual cases here
        tester = TestPlateCases()
        tester.test_uniaxial_compression()

    elif test == 2:
        tester = TestPlateCases()
        tester.test_pure_shear_clamped()

    elif test == 3:
        tester = TestPlateCases()
        tester.test_combined_axial_shear_clamped()
