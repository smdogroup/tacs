__all__ = ["FlatPlateAnalysis"]

import numpy as np
import matplotlib.pyplot as plt
from tacs import pyTACS, constitutive, elements, utilities
import os
from pprint import pprint

dtype = utilities.BaseUI.dtype


class FlatPlateAnalysis:
    def __init__(
        self, comm, bdf_file, a, b, h, E11, nu12, E22=None, G12=None, _G23=None, _G13=None
    ):
        self.comm = comm

        # geometry properties
        self.a = a  # Lx
        self.b = b  # Ly
        self.h = h  # plate thickness

        # material properties
        self.E11 = E11
        self.nu12 = nu12
        self._E22 = E22
        self._G12 = G12
        self._G23 = _G23
        self._G13 = _G13

        self._bdf_file = bdf_file

    @property
    def bdf_file(self) -> str:
        return self._bdf_file

    @bdf_file.setter
    def bdf_file(self, new_file: str):
        self._bdf_file = new_file

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
        get the exx so that lambda = kx_0 the affine buckling coefficient
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
        """TODO : write this eqn out"""
        return None
    
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

        x = self.a * np.linspace(0.0, 1.0, nx + 1)
        y = self.b * np.linspace(0.0, 1.0, ny + 1)

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
                                "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "36", 0.0)
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
            static_folder = os.path.join(base_path, "static")
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
            buckling_folder = os.path.join(base_path, "buckling")
            if not os.path.exists(buckling_folder):
                os.mkdir(buckling_folder)
            bucklingProb.writeSolution(outputDir=buckling_folder)

        errors = []
        for imode in range(num_eig):
            eigval,eigvec = bucklingProb.getVariables(imode)
            error = bucklingProb.getModalError(imode)
            errors += [error]

        if self.comm.rank == 0:
            pprint(funcs)
            # pprint(funcsSens)

        # return the eigenvalues here
        return np.array([funcs[key] for key in funcs]), np.array(errors)
