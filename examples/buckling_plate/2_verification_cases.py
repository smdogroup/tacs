"""
Sean Engelstad, Feb 2024
GT SMDO Lab
"""
from tacs import buckling_surrogate
import numpy as np
import unittest
from mpi4py import MPI

comm = MPI.COMM_WORLD

# 3 main verification cases for the buckling analysis
# 1, 2, 3
case = 3
run_static = False


class TestPlateCases(unittest.TestCase):
    # use the same flat plate geometry and material through all three test cases
    flat_plate = buckling_surrogate.FlatPlateAnalysis(
        comm=comm,
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

    def test_affine_simply_supported(self):
        flat_plate = buckling_surrogate.FlatPlateAnalysis(
            comm=comm,
            bdf_file="plate.bdf",
            a=1.0,
            b=1.0,
            h=0.005,  # very slender => so near thin plate limit
            E11=70e9,
            nu12=0.33,
            E22=None,  # set to None if isotropic
            G12=None,  # set to None if isotropic
        )

        flat_plate.generate_bdf(
            nx=20,
            ny=20,
            exx=flat_plate.affine_exx,
            eyy=0.0,
            exy=0.0,
            clamped=False,
        )

        # avg_stresses = flat_plate.run_static_analysis(write_soln=True)

        tacs_eigvals, _ = flat_plate.run_buckling_analysis(
            sigma=10.0, num_eig=12, write_soln=True
        )

        # expect to get ~4.5
        # since k0-2D* = (m a_0/b_0)^2 + (b_0/a_0/m)^2
        # in affine space and D*=1 and k0-2D* = 2.5 in Brunelle paper (buckling-analysis section)
        # "Generic Buckling Curves For Specially Orthotropic Rectangular Plates"
        # but only works in thin plate limit (very thin)

        kx_0 = tacs_eigvals[0]  # exx_affine was scaled so that lambda = kx_0
        # the non-dimensional buckling coefficient

        kx_0_exact = 2.5 + 2.0 * flat_plate.Dstar
        rel_err = (kx_0 - kx_0_exact) / kx_0_exact

        print(f"Case 4 - Compare affine curve simply supported uniaxial compression..")
        print(
            f"\tDstar = {flat_plate.Dstar}, b/h = {flat_plate.slenderness} (near thin plate limit)"
        )
        print(f"\tkx0 pred = {kx_0}")
        print(f"\tkx0 exact = {kx_0_exact}")
        print(f"\tkx0 rel error = {rel_err}")

        assert abs(rel_err) < 0.05  # less than 5% error


if __name__ == "__main__":
    test = 4

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

    elif test == 4:
        tester = TestPlateCases()
        tester.test_affine_simply_supported()
