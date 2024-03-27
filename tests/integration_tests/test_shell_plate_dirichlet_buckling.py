"""
Sean Engelstad, Feb 2024
GT SMDO Lab
"""
import numpy as np
import unittest
from mpi4py import MPI

comm = MPI.COMM_WORLD

import numpy as np
from tacs import pyTACS, constitutive, elements, utilities
import os
from pprint import pprint

dtype = utilities.BaseUI.dtype

class FlatPlateAnalysis:
    def __init__(
        self,
        comm,
        bdf_file,
        a,
        b,
        h,
        E11,
        nu12,
        E22=None,
        G12=None,
        _G23=None,
        _G13=None,
        material_name=None,
        ply_angle=None,
        plate_name=None,  # use the plate name to differentiate plate folder names
    ):
        self.comm = comm

        # geometry properties
        self.a = a  # Lx
        self.b = b  # Ly
        self.h = h  # plate thickness

        # material properties
        self.material_name = material_name
        self.ply_angle = ply_angle
        self.E11 = E11
        self.nu12 = nu12
        self._E22 = E22
        self._G12 = G12
        self._G23 = _G23
        self._G13 = _G13

        self._plate_name = plate_name

        self._bdf_file = bdf_file

    @property
    def static_folder_name(self) -> str:
        if self._plate_name:
            return "static-" + self._plate_name
        else:
            return "static"

    @property
    def buckling_folder_name(self) -> str:
        if self._plate_name:
            return "buckling-" + self._plate_name
        else:
            return "buckling"

    @property
    def bdf_file(self) -> str:
        return self._bdf_file

    @bdf_file.setter
    def bdf_file(self, new_file: str):
        self._bdf_file = new_file

    @property
    def num_elements(self) -> int:
        return self._nx * self._ny

    @property
    def nu21(self) -> float:
        """reversed 12 Poisson's ratio"""
        return self.nu12 * self.E22 / self.E11

    @property
    def E22(self) -> float:
        if self._E22 is None:
            return self.E11
        else:
            return self._E22

    @property
    def G12(self) -> float:
        if self._G12 is None:
            return self.E11 / 2.0 / (1 + self.nu12)
        else:
            return self._G12

    @property
    def D11(self) -> float:
        """1-axial bending stiffness"""
        return self.E11 * self.h**3 / 12.0 / (1 - self.nu12 * self.nu21)

    @property
    def D22(self) -> float:
        """2-axial bending stiffness"""
        return self.E22 * self.h**3 / 12.0 / (1 - self.nu12 * self.nu21)

    @property
    def D12(self) -> float:
        return self.nu12 * self.D22

    @property
    def D66(self) -> float:
        return self.G12 * self.h**3 / 12.0

    @property
    def affine_exx(self):
        """
        get the exx so that lambda = kx_0 the affine buckling coefficient for pure axial load
        out of the buckling analysis!
        """
        exx_T = (
            np.pi**2 * np.sqrt(self.D11 * self.D22) / self.b**2 / self.h / self.E11
        )
        return exx_T

    @property
    def affine_eyy(self):
        """TODO : write this eqn out"""
        return None

    @property
    def affine_exy(self):
        """
        get the exy so that lambda = kx_0y_0 the affine buckling coefficient for pure shear load
        out of the buckling analysis!
        """
        option = 2
        # option 1 - based on self derivation (but didn't match data well)
        if option == 1:
            exy_T = np.pi**2 * (self.D11 * self.D22)**0.5 / self.a / self.b / self.h / self.G12
        # option 2 - based on NASA non-dimensional buckling parameter derivation (much better)
        elif option == 2:
            exy_T = (
                np.pi**2
                * (self.D11 * self.D22**3) ** 0.25
                / self.b**2
                / self.h
                / self.G12
            )
        return exy_T

    @property
    def aspect_ratio(self):
        """a/b ratio of the plate usually between 5 and 200"""
        return self.a / self.b

    @property
    def affine_aspect_ratio(self):
        """a0/b0 aspect ratio in the affine space"""
        return (self.D22 / self.D11) ** 0.25 * self.aspect_ratio

    @property
    def Dstar(self):
        """
        return Dstar the generalized rigidity from the affine transformatin of orthotropic CPT (Classical Plate Theory)
        TODO : may need to add isotropic case to this with Dstar = 1.0
        """
        return (self.D12 + 2 * self.D66) / np.sqrt(self.D11 * self.D22)

    @property
    def slenderness(self):
        """
        slenderness ratio b/h. As the slenderness ratio decreases (plate is thicker) the CPT (Classical Plate Theory)
        which doesn't include shear strain energy is less exact. And FSDT (First Order Shear Deformation Theory) or Reissler-Mindlin plate
        theory which includes shear is required.
        """
        return self.b / self.h

    def generate_bdf(self, nx=30, ny=30, exx=0.0, eyy=0.0, exy=0.0, clamped=True):
        """
        # Generate a plate mesh with CQUAD4 elements
        create pure axial, pure shear, or combined loading displacement control
        of a flat plate
        """

        nodes = np.arange(1, (nx + 1) * (ny + 1) + 1, dtype=np.int32).reshape(
            nx + 1, ny + 1
        )

        self._nx = nx
        self._ny = ny

        # num nodes in each direction
        self._M = self._nx + 1
        self._N = self._ny + 1

        x = self.a * np.linspace(0.0, 1.0, nx + 1)
        y = self.b * np.linspace(0.0, 1.0, ny + 1)

        # copy coordinates for use later in the modal assurance criterion
        self._x = x
        self._y = y

        if self.comm.rank == 0:
            fp = open(self.bdf_file, "w")
            fp.write("$ Input file for a square axial/shear-disp BC plate\n")
            fp.write("SOL 103\nCEND\nBEGIN BULK\n")

            # Write the grid points to a file
            for j in range(ny + 1):
                for i in range(nx + 1):
                    # Write the nodal data
                    spc = " "
                    coord_disp = 0
                    coord_id = 0
                    seid = 0

                    fp.write(
                        "%-8s%16d%16d%16.9e%16.9e*       \n"
                        % ("GRID*", nodes[i, j], coord_id, x[i], y[j])
                    )
                    fp.write(
                        "*       %16.9e%16d%16s%16d        \n"
                        % (0.0, coord_disp, spc, seid)
                    )

            # Output 2nd order elements
            elem = 1
            part_id = 1
            for j in range(ny):
                for i in range(nx):
                    # Write the connectivity data
                    # CQUAD4 elem id n1 n2 n3 n4
                    fp.write(
                        "%-8s%8d%8d%8d%8d%8d%8d\n"
                        % (
                            "CQUAD4",
                            elem,
                            part_id,
                            nodes[i, j],
                            nodes[i + 1, j],
                            nodes[i + 1, j + 1],
                            nodes[i, j + 1],
                        )
                    )
                    elem += 1

            # Set up the plate BCs so that it has u = uhat, for shear disp control
            # u = eps * y, v = eps * x, w = 0
            for j in range(ny + 1):
                for i in range(nx + 1):
                    u = exy * y[j]
                    v = exy * x[i]

                    if i == nx or exy != 0:
                        u -= exx * x[i]
                    elif j == ny:
                        v -= eyy * y[j]
                    elif i == 0 or j == 0:
                        pass

                    # check on boundary
                    if i == 0 or j == 0 or i == nx or j == ny:
                        if clamped or (i == 0 and j == 0):
                            fp.write(
                                "%-8s%8d%8d%8s%8.6f\n"
                                % ("SPC", 1, nodes[i, j], "3456", 0.0)
                            )  # w = theta_x = theta_y
                        else:
                            fp.write(
                                "%-8s%8d%8d%8s%8.6f\n"
                                % ("SPC", 1, nodes[i, j], "36", 0.0)
                            )  # w = theta_x = theta_y
                        if exy != 0 or i == 0 or i == nx:
                            fp.write(
                                "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "1", u)
                            )  # u = eps_xy * y
                        if exy != 0.0 or j == 0:
                            fp.write(
                                "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "2", v)
                            )  # v = eps_xy * x

            # # plot the mesh to make sure it makes sense
            # X, Y = np.meshgrid(x, y)
            # W = X * 0.0
            # for j in range(ny + 1):
            #     for i in range(nx + 1):
            #         if i == 0 or i == nx or j == 0 or j == ny:
            #             W[i, j] = 1.0

            # plt.scatter(X,Y)
            # plt.contour(X,Y,W, corner_mask=True, antialiased=True)
            # plt.show()

            fp.write("ENDDATA")
            fp.close()

        self.comm.Barrier()

    def _elemCallBack(self):
        def elemCallBack(
            dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs
        ):
            # Set constitutive properties
            # rho = 4540.0  # density, kg/m^3
            # # E = 70e9  # elastic modulus, Pa 118e9
            # # nu = 0.33  # poisson's ratio
            # ys = 1050e6  # yield stress, Pa
            # kappa = 6.89
            # specific_heat = 463.0

            # if E22 not provided, isotropic
            isotropic = self._E22 is None

            # Setup property and constitutive objects
            if isotropic:
                mat = constitutive.MaterialProperties(E=self.E11, nu=self.nu12)

                # Set one thickness dv for every component
                con = constitutive.IsoShellConstitutive(mat, t=self.h)

            else:  # orthotropic
                # assume G23, G13 = G12
                G23 = self.G12 if self._G23 is None else self._G23
                G13 = self.G12 if self._G13 is None else self._G13
                ortho_prop = constitutive.MaterialProperties(
                    E1=self.E11,
                    E2=self.E22,
                    nu12=self.nu12,
                    G12=self.G12,
                    G23=G23,
                    G13=G13,
                )

                ortho_ply = constitutive.OrthotropicPly(self.h, ortho_prop)

                # one play composite constitutive model
                con = constitutive.CompositeShellConstitutive(
                    [ortho_ply],
                    np.array([self.h], dtype=dtype),
                    np.array([0], dtype=dtype),
                    tOffset=0.0,
                )
            # For each element type in this component,
            # pass back the appropriate tacs element object
            elemList = []
            for descript in elemDescripts:
                transform = None
                if descript in ["CQUAD4", "CQUADR"]:
                    elem = elements.Quad4Shell(transform, con)
                elif descript in ["CQUAD9", "CQUAD"]:
                    elem = elements.Quad9Shell(transform, con)
                else:
                    raise AssertionError("Non CQUAD4 Elements in this plate?")

                elemList.append(elem)

            # Add scale for thickness dv
            scale = [100.0]
            return elemList, scale
        return elemCallBack

    def run_static_analysis(self, base_path=None, write_soln=False):
        """
        run a linear static analysis on the flat plate with either isotropic or composite materials
        return the average stresses in the plate => to compute in-plane loads Nx, Ny, Nxy
        """

        # Instantiate FEAAssembler
        FEAAssembler = pyTACS(self.bdf_file, comm=self.comm)

        # Set up constitutive objects and elements
        FEAAssembler.initialize(self._elemCallBack())

        # set complex step Gmatrix into all elements through assembler
        FEAAssembler.assembler.setComplexStepGmatrix(True)

        # debug the static problem first
        SP = FEAAssembler.createStaticProblem(name="static")
        SP.solve()
        if write_soln:
            if base_path is None:
                base_path = os.getcwd()
            static_folder = os.path.join(base_path, self.static_folder_name)
            if not os.path.exists(static_folder):
                os.mkdir(static_folder)
            SP.writeSolution(outputDir=static_folder)

        # test the average stresses routine
        avgStresses = FEAAssembler.assembler.getAverageStresses()
        return avgStresses

    def run_buckling_analysis(
        self,
        sigma=30.0,
        num_eig=5,
        write_soln=False,
        derivatives=False,
        base_path=None,
    ):
        """
        run a linear buckling analysis on the flat plate with either isotropic or composite materials
        return the sorted eigenvalues of the plate => would like to include M
        """

        # Instantiate FEAAssembler
        FEAAssembler = pyTACS(self.bdf_file, comm=self.comm)

        # Set up constitutive objects and elements
        FEAAssembler.initialize(self._elemCallBack())

        # set complex step Gmatrix into all elements through assembler
        FEAAssembler.assembler.setComplexStepGmatrix(True)

        # Setup buckling problem
        bucklingProb = FEAAssembler.createBucklingProblem(
            name="buckle", sigma=sigma, numEigs=num_eig
        )
        bucklingProb.setOption("printLevel", 2)

        # exit()

        # solve and evaluate functions/sensitivities
        funcs = {}
        funcsSens = {}
        bucklingProb.solve()
        bucklingProb.evalFunctions(funcs)
        if derivatives:
            bucklingProb.evalFunctionsSens(funcsSens)
        if write_soln:
            if base_path is None:
                base_path = os.getcwd()
            buckling_folder = os.path.join(base_path, self.buckling_folder_name)
            if not os.path.exists(buckling_folder):
                os.mkdir(buckling_folder)
            bucklingProb.writeSolution(outputDir=buckling_folder)

        # save the eigenvectors for MAC and return errors from function
        self._eigenvectors = []
        self._eigenvalues = []
        self._num_modes = num_eig
        errors = []
        for imode in range(num_eig):
            eigval, eigvec = bucklingProb.getVariables(imode)
            self._eigenvectors += [eigvec]
            self._eigenvalues += [eigval]
            error = bucklingProb.getModalError(imode)
            errors += [error]

        if self.comm.rank == 0:
            pprint(funcs)
            # pprint(funcsSens)

        self._solved_buckling = True
        self._alphas = {}

        # return the eigenvalues here
        return np.array([funcs[key] for key in funcs]), np.array(errors)


class TestPlateCases(unittest.TestCase):
    # use the same flat plate geometry and material through all three test cases
    flat_plate = FlatPlateAnalysis(
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
            nx=10,
            ny=10,
            exx=0.001,
            eyy=0.0,
            exy=0.0,
            clamped=False,
        )

        # flat_plate.run_static_analysis(write_soln=True)

        tacs_eigvals,_ = self.flat_plate.run_buckling_analysis(
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

        # relaxed the max relative error requirement due to very coarse mesh
        assert np.max(np.abs(rel_error)) < 0.2

    def test_pure_shear_clamped(self):
        # case 2 - pure shear, clamped plate
        self.flat_plate.generate_bdf(
            nx=10,
            ny=10,
            exx=0.0,
            eyy=0.0,
            exy=0.001,
            clamped=True,
        )

        # flat_plate.run_static_analysis(write_soln=True)

        tacs_eigvals,_ = self.flat_plate.run_buckling_analysis(
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

        # relaxed the max relative error requirement due to very coarse mesh
        assert np.max(np.abs(rel_error)) < 0.2

    def test_combined_axial_shear_clamped(self):
        # case 3 - mixed shear and uniaxial compression case, clamped plate
        self.flat_plate.generate_bdf(
            nx=10,
            ny=10,
            exx=0.001,
            eyy=0.0,
            exy=0.001,
            clamped=True,
        )

        # flat_plate.run_static_analysis(write_soln=True)

        tacs_eigvals,_ = self.flat_plate.run_buckling_analysis(
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

        # relaxed the max relative error requirement due to very coarse mesh
        assert np.max(np.abs(rel_error)) < 0.2

    def test_affine_simply_supported(self):
        flat_plate = FlatPlateAnalysis(
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
            nx=10,
            ny=10,
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

        # relaxed the max relative error requirement due to very coarse mesh
        assert np.max(np.abs(rel_err)) < 0.2


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

    elif test == 4:
        tester = TestPlateCases()
        tester.test_affine_simply_supported()
