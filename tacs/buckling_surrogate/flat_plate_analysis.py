__all__ = ["FlatPlateAnalysis", "exp_kernel1"]

import numpy as np
from tacs import pyTACS, constitutive, elements, utilities
import os
from pprint import pprint
from .composite_material_utility import CompositeMaterialUtility
from typing_extensions import Self

dtype = utilities.BaseUI.dtype


def exp_kernel1(xp, xq, sigma_f, L):
    # xp, xq are Nx1, Nx1 vectors
    return sigma_f**2 * np.exp(-0.5 * (xp - xq).T @ (xp - xq) / L**2)


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

        # temp
        self._nx = None
        self._ny = None
        self._M = None
        self._N = None
        self._num_modes = None
        self._eigenvectors = None
        self._eigenvalues = None
        self._alphas = None
        self._solved_buckling = False
        self._saved_alphas = False

    MAC_THRESHOLD = 0.1  # 0.6

    @classmethod
    def mac_permutation(
        cls, nominal_plate: Self, new_plate: Self, num_modes: int
    ) -> dict:
        """
        compute the permutation of modes in the new plate that correspond to the modes in the nominal plate
        using 2D Discrete Fourier transform in a model assurance criterion
        """
        eigenvalues = [None for _ in range(num_modes)]
        permutation = {}
        nominal_interp_modes = nominal_plate.interpolate_eigenvectors(
            X_test=new_plate.nondim_X
        )
        new_modes = new_plate.eigenvectors

        _debug = False
        if _debug:
            for imode, interp_mode in enumerate(nominal_interp_modes):
                interp_mat = new_plate._vec_to_plate_matrix(interp_mode)
                import matplotlib.pyplot as plt

                plt.imshow(interp_mat.astype(np.double))
                plt.show()

        for imode, nominal_mode in enumerate(nominal_interp_modes):
            if imode >= num_modes:  # if larger than number of nominal modes to compare
                break
            nominal_mode_unit = nominal_mode / np.linalg.norm(nominal_mode)

            similarity_list = []
            for new_mode in new_modes:
                new_mode_unit = new_mode / np.linalg.norm(new_mode)
                # compute cosine similarity with the unit vectors
                similarity_list += [
                    abs(np.dot(nominal_mode_unit, new_mode_unit).astype(np.double))
                ]

            # compute the maximum similarity index
            if _debug:
                print(f"similarity list imode {imode} = {similarity_list}")
            jmode_star = np.argmax(np.array(similarity_list))
            permutation[imode] = jmode_star
            if similarity_list[jmode_star] > cls.MAC_THRESHOLD:  # similarity threshold
                eigenvalues[imode] = new_plate.eigenvalues[jmode_star]

        # check the permutation map is one-to-one

        # print to the terminal about the modal criterion
        print(f"0-based mac criterion permutation map")
        print(
            f"\tbetween nominal plate {nominal_plate._plate_name} and new plate {new_plate._plate_name}"
        )
        print(f"\tthe permutation map is the following::\n")
        for imode in range(num_modes):
            print(f"\t nominal {imode} : new {permutation[imode]}")

        return eigenvalues, permutation

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

    @classmethod
    def get_materials(cls):
        return [
            cls.solvay5320,
            cls.solvayMTM45,
            cls.torayBT250E,
            cls.hexcelIM7,
            cls.victrexAE,
        ]
    
    @classmethod
    def get_material_from_str(cls, mat_name:str):
        method_names = [_.__qualname__ for _ in cls.get_materials()]
        materials = cls.get_materials()
        _method = None
        for i,method_name in enumerate(method_names):
            if mat_name in method_name:
                _method = materials[i]
        assert _method is not None
        return _method

    # MATERIALS CLASS METHODS
    # -----------------------------------------------------------

    # NIAR composite materials

    @classmethod
    def solvay5320(cls, comm, bdf_file, a, b, h, ply_angle=0.0):
        """
        NIAR dataset - Solvay 5320-1 material (thermoset)
        Fiber: T650 unitape, Resin: Cycom 5320-1
        Room Temperature Dry (RTD) mean properties shown below
        units in Pa, ND
        """
        comp_utility = CompositeMaterialUtility(
            E11=138.461e9, E22=9.177e9, nu12=0.326, G12=4.957e9
        )
        comp_utility.rotate_ply(ply_angle)

        return cls(
            comm=comm,
            bdf_file=bdf_file,
            a=a,
            b=b,
            h=h,
            material_name="solvay5320",
            ply_angle=ply_angle,
            E11=comp_utility.E11,
            E22=comp_utility.E22,
            nu12=comp_utility.nu12,
            G12=comp_utility.G12,
        )

    @classmethod
    def solvayMTM45(cls, comm, bdf_file, a, b, h, ply_angle=0.0):
        """
        NIAR dataset - Solvay MTM45 material (thermoset)
        Style: 12K AS4 Unidirectional
        Room Temperature Dry (RTD) mean properties shown below
        units in Pa, ND
        """
        comp_utility = CompositeMaterialUtility(
            E11=129.5e9, E22=7.936e9, nu12=0.313, G12=4.764e9
        )
        comp_utility.rotate_ply(ply_angle)

        return cls(
            comm=comm,
            bdf_file=bdf_file,
            a=a,
            b=b,
            h=h,
            material_name="solvayMTM45",
            ply_angle=ply_angle,
            E11=comp_utility.E11,
            E22=comp_utility.E22,
            nu12=comp_utility.nu12,
            G12=comp_utility.G12,
        )

    @classmethod
    def torayBT250E(cls, comm, bdf_file, a, b, h, ply_angle=0.0):
        """
        NIAR dataset - Toray (formerly Tencate) BT250E-6 S2 Unitape Gr 284 material (thermoset)
        Room Temperature Dry (RTD) mean properties shown below
        units in Pa, ND
        """
        comp_utility = CompositeMaterialUtility(
            E11=44.74e9, E22=11.36e9, nu12=0.278, G12=3.77e9
        )
        comp_utility.rotate_ply(ply_angle)

        return cls(
            comm=comm,
            bdf_file=bdf_file,
            a=a,
            b=b,
            h=h,
            material_name="torayBT250E",
            ply_angle=ply_angle,
            E11=comp_utility.E11,
            E22=comp_utility.E22,
            nu12=comp_utility.nu12,
            G12=comp_utility.G12,
        )

    @classmethod
    def victrexAE(cls, comm, bdf_file, a, b, h, ply_angle=0.0):
        """
        NIAR dataset - Victrex AE 250 LMPAEK (thermoplastic)
        Room Temperature Dry (RTD) mean properties shown below
        units in Pa, ND
        """
        comp_utility = CompositeMaterialUtility(
            E11=131.69e9, E22=9.694e9, nu12=0.3192, G12=4.524e9
        )
        comp_utility.rotate_ply(ply_angle)

        return cls(
            comm=comm,
            bdf_file=bdf_file,
            a=a,
            b=b,
            h=h,
            material_name="victrexAE",
            ply_angle=ply_angle,
            E11=comp_utility.E11,
            E22=comp_utility.E22,
            nu12=comp_utility.nu12,
            G12=comp_utility.G12,
        )

    @classmethod
    def hexcelIM7(cls, comm, bdf_file, a, b, h, ply_angle=0.0):
        """
        NIAR dataset - Hexcel 8552 IM7 Unidirectional Prepreg (thermoset)
        Room Temperature Dry (RTD) mean properties shown below
        units in Pa, ND
        """
        comp_utility = CompositeMaterialUtility(
            E11=158.51e9, nu12=0.316, E22=8.96e9, G12=4.688e9
        )
        comp_utility.rotate_ply(ply_angle)

        return cls(
            comm=comm,
            bdf_file=bdf_file,
            a=a,
            b=b,
            h=h,
            material_name="hexcelIM7",
            ply_angle=ply_angle,
            E11=comp_utility.E11,
            E22=comp_utility.E22,
            nu12=comp_utility.nu12,
            G12=comp_utility.G12,
        )

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

    def _vec_to_plate_matrix(self, vec):
        """
        build a matrix for transforming column vectors MN x 1 to M x N matrices
        for the mesh for Modal assurance criterion
        """
        return np.array(
            [[vec[self._M * j + i] for j in range(self._N)] for i in range(self._M)],
            dtype=dtype,
        )

    def get_eigenvector(self, imode):
        # convert eigenvectors to w coordinates only, 6 dof per shell
        # print(f"ndof in eigenvector = {self._eigenvectors[imode].shape[0]}")
        return self._eigenvectors[imode][2::6]

    @property
    def eigenvectors(self):
        return [self.get_eigenvector(imode) for imode in range(self.num_modes)]

    @property
    def eigenvalues(self):
        return self._eigenvalues

    @property
    def num_nodes(self) -> int:
        return self._N * self._M

    @property
    def nondim_X(self):
        """non-dimensional X matrix for Gaussian Process model"""
        return np.concatenate(
            [np.expand_dims(self._xi, axis=-1), np.expand_dims(self._eta, axis=-1)],
            axis=1,
        )

    def interpolate_eigenvectors(self, X_test, compute_covar=False):
        """
        interpolate the eigenvector from this object the nominal plate to a new mesh in non-dim coordinates
        """
        X_train = self.nondim_X
        num_train = self.num_nodes
        num_test = X_test.shape[0]
        # default hyperparameters
        sigma_n = 1e-4
        sigma_f = 1.0
        L = 0.4
        _kernel = lambda xp, xq: exp_kernel1(xp, xq, sigma_f=sigma_f, L=L)
        K_train = sigma_n**2 * np.eye(num_train) + np.array(
            [
                [_kernel(X_train[i, :], X_train[j, :]) for i in range(num_train)]
                for j in range(num_train)
            ]
        )
        K_test = np.array(
            [
                [_kernel(X_train[i, :], X_test[j, :]) for i in range(num_train)]
                for j in range(num_test)
            ]
        )

        if not compute_covar:
            _interpolated_eigenvectors = []
            for imode in range(self.num_modes):
                phi = self.get_eigenvector(imode)
                if self._saved_alphas:  # skip linear solve in this case
                    alpha = self._alphas[imode]
                else:
                    alpha = np.linalg.solve(K_train, phi)
                    self._alphas[imode] = alpha
                phi_star = K_test @ alpha
                _interpolated_eigenvectors += [phi_star]
            self._saved_alphas = True
            return _interpolated_eigenvectors
        else:
            raise AssertionError(
                "Haven't written part of extrapolate eigenvector to get the conditional covariance yet."
            )

    # decided no longer to do the discrete fourier transform approach for modal assurance criterion
    # -----------------------------------------------------------------------------------------------
    #
    # def _dft_matrix(self, k, l):
    #     """compute the Fourier matrix np.exp(-2*pi*sqrt(-1) * (km/M + ln/N)) in inline for loops (faster)"""
    #     return np.array([[np.exp(-2 * np.pi*1j * (k*m/self._M + l*n/self._N)) for n in range(self._N)] for m in range(self._M)])

    # def twod_DFT(self, eig_matrix):
    #     """
    #     perform the 2D Discrete Fourier Transform on the eigenvector
    #     in a meshgrid / matrix form (since this is a 2D mesh)
    #     report the amplitudes of the first 25 Fourier entries
    #     """

    #     dft_amplitudes = []
    #     for k in range(5):
    #         for l in range(5):
    #             # compute Phi_kl using DFT
    #             _debug = True
    #             if _debug: # see the eigenvectors
    #                 import matplotlib.pyplot as plt
    #                 plt.title(f"k = {k}, l = {l}")
    #                 plt.imshow(np.real(self._dft_matrix(k,l)))
    #                 plt.show()
    #             Phi_kl = 1.0 / self._M / self._N * np.sum(eig_matrix * self._dft_matrix(k,l))
    #             dft_amplitudes += [np.abs(Phi_kl)]
    #     return np.array(dft_amplitudes)
    # def get_fourier_eigenvector(self, imode):
    #     """
    #     get the 2D DFT fourier amplitude vector for each mode
    #     used for cosine similarity among different flat plates
    #     """
    #     eigvec = self.get_eigenvector(imode)
    #     eigvec_matrix = self._vec_to_plate_matrix(eigvec)
    #     _debug = False
    #     if _debug: # see the eigenvectors
    #         import matplotlib.pyplot as plt
    #         plt.imshow(eigvec_matrix.astype(np.double))
    #         plt.show()
    #     return self.twod_DFT(eigvec_matrix)

    @property
    def num_modes(self) -> int:
        """number of eigenvalues or modes that were recorded"""
        return self._num_modes

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
        # nodes count across the x-axis first then loop over y-axis from bot-left corner
        self._xi = [
            self._x[i % self._M] / self.a for i in range(self.num_nodes)
        ]  # non-dim coordinates
        self._eta = [self._y[int(i / self._M)] / self.b for i in range(self.num_nodes)]

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

    def run_static_analysis(self, base_path=None, write_soln=False):
        """
        run a linear static analysis on the flat plate with either isotropic or composite materials
        return the average stresses in the plate => to compute in-plane loads Nx, Ny, Nxy
        """

        # Instantiate FEAAssembler
        FEAAssembler = pyTACS(self.bdf_file, comm=self.comm)

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

        # Set up constitutive objects and elements
        FEAAssembler.initialize(elemCallBack)

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
        print(f"avg Stresses = {avgStresses}")
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

        # Set up constitutive objects and elements
        FEAAssembler.initialize(elemCallBack)

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
