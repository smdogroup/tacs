from tacs import TACS, constitutive, elements
import numpy as np
import unittest


class ElementTest(unittest.TestCase):
    def setUp(self):
        max_nodes = 64
        max_vars_per_nodes = 8
        max_vars = max_nodes * max_vars_per_nodes

        # fd/cs step size
        if TACS.dtype is complex:
            self.dh = 1e-50
            self.rtol = 1e-10
        else:
            self.dh = 1e-5
            self.rtol = 1e-2
        self.dtype = TACS.dtype

        # Basically, only check relative tolerance
        self.atol = 1e99
        self.print_level = 0

        # Set element index
        self.elem_index = 0

        # Set the simulation time
        self.time = 0.0

        # Set the variable arrays
        np.random.seed(30)  # Seed random numbers for deterministic/repeatable tests
        self.xpts = np.random.rand(3 * max_nodes).astype(self.dtype)
        np.random.seed(30)  # Seed random numbers for deterministic/repeatable tests
        self.vars = np.random.rand(max_vars).astype(self.dtype)
        self.dvars = self.vars.copy()
        self.ddvars = self.vars.copy()

        # Create the isotropic material
        rho = 2700.0
        specific_heat = 921.096
        E = 70e3
        nu = 0.3
        ys = 270.0
        cte = 24.0e-6
        kappa = 230.0
        self.props = constitutive.MaterialProperties(
            rho=rho,
            specific_heat=specific_heat,
            E=E,
            nu=nu,
            ys=ys,
            cte=cte,
            kappa=kappa,
        )

        ref_axis = np.array([0.0, 1.0, 1.0], dtype=self.dtype)
        self.transforms = [
            elements.ShellNaturalTransform(),
            elements.ShellRefAxisTransform(ref_axis),
        ]

        # TACS shell elements of various orders and types
        self.elements = [
            elements.Tri3Shell,
            elements.Quad4Shell,
            elements.Quad9Shell,
            elements.Quad16Shell,
            elements.Tri3ThermalShell,
            elements.Quad4ThermalShell,
            elements.Quad9ThermalShell,
            elements.Quad16ThermalShell,
            elements.Quad4NonlinearShell,
            elements.Quad9NonlinearShell,
            elements.Quad16NonlinearShell,
            elements.Tri3NonlinearShell,
            elements.Quad4NonlinearThermalShell,
            elements.Quad9NonlinearThermalShell,
            elements.Quad16NonlinearThermalShell,
            elements.Tri3NonlinearThermalShell,
        ]

        # The thermal elements will not pass the residual test since they are not derived
        # from Lagrange's equations due to the presence of the thermal coupling equations.
        # They also fail the mat_sv_sens test, so we skip them there as well
        self.thermal_elements = [
            elements.Tri3ThermalShell,
            elements.Quad4ThermalShell,
            elements.Quad9ThermalShell,
            elements.Quad16ThermalShell,
            elements.Quad4NonlinearThermalShell,
            elements.Quad9NonlinearThermalShell,
            elements.Quad16NonlinearThermalShell,
            elements.Tri3NonlinearThermalShell,
        ]

        # Create stiffness (need class)
        self.con = constitutive.IsoShellConstitutive(self.props, t=1.0, tNum=0)

        # Set matrix types
        self.matrix_types = [
            TACS.STIFFNESS_MATRIX,
            TACS.MASS_MATRIX,
            TACS.GEOMETRIC_STIFFNESS_MATRIX,
        ]

        # Seed random number generator in tacs for consistent test results
        elements.SeedRandomGenerator(0)

    def test_element_residual(self):
        # Here we have to overwrite the step size rtol,
        # because TestElementResidual only supports FD testing right now
        dh = 1e-5
        rtol = 1e-2
        # Loop through every combination of transform type and shell element class and test residual
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for element_handle in self.elements:
                    if not (element_handle in self.thermal_elements):
                        with self.subTest(element=element_handle):
                            element = element_handle(transform, self.con)
                            fail = elements.TestElementResidual(
                                element,
                                self.elem_index,
                                self.time,
                                self.xpts,
                                self.vars,
                                self.dvars,
                                self.ddvars,
                                dh,
                                self.print_level,
                                self.atol,
                                rtol,
                            )
                            self.assertFalse(fail)

    def test_element_jacobian(self):
        # Loop through every combination of transform type and shell element class and test Jacobian
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for element_handle in self.elements:
                    with self.subTest(element=element_handle):
                        element = element_handle(transform, self.con)
                        fail = elements.TestElementJacobian(
                            element,
                            self.elem_index,
                            self.time,
                            self.xpts,
                            self.vars,
                            self.dvars,
                            self.ddvars,
                            -1,
                            self.dh,
                            self.print_level,
                            self.atol,
                            self.rtol,
                        )
                        self.assertFalse(fail)

    def test_adj_res_product(self):
        # Loop through every combination of transform type and shell element class and test adjoint residual-dvsens product
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for element_handle in self.elements:
                    with self.subTest(element=element_handle):
                        element = element_handle(transform, self.con)
                        dvs = element.getDesignVars(self.elem_index)
                        fail = elements.TestAdjResProduct(
                            element,
                            self.elem_index,
                            self.time,
                            self.xpts,
                            self.vars,
                            self.dvars,
                            self.ddvars,
                            dvs,
                            self.dh,
                            self.print_level,
                            self.atol,
                            self.rtol,
                        )
                        self.assertFalse(fail)

    def test_adj_res_xpt_product(self):
        # Loop through every combination of transform type and shell element class and test adjoint residual-xptsens product
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for element_handle in self.elements:
                    with self.subTest(element=element_handle):
                        element = element_handle(transform, self.con)
                        fail = elements.TestAdjResXptProduct(
                            element,
                            self.elem_index,
                            self.time,
                            self.xpts,
                            self.vars,
                            self.dvars,
                            self.ddvars,
                            self.dh,
                            self.print_level,
                            self.atol,
                            self.rtol,
                        )
                        self.assertFalse(fail)

    def test_element_mat_dv_sens(self):
        # Loop through every combination of transform type and shell element class and element matrix inner product sens
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for element_handle in self.elements:
                    with self.subTest(element=element_handle):
                        element = element_handle(transform, self.con)
                        dvs = element.getDesignVars(self.elem_index)
                        for matrix_type in self.matrix_types:
                            with self.subTest(matrix_type=matrix_type):
                                fail = elements.TestElementMatDVSens(
                                    element,
                                    matrix_type,
                                    self.elem_index,
                                    self.time,
                                    self.xpts,
                                    self.vars,
                                    dvs,
                                    self.dh,
                                    self.print_level,
                                    self.atol,
                                    self.rtol,
                                )
                                self.assertFalse(fail)

    def test_element_mat_xpt_sens(self):
        # Loop through every combination of transform type and shell element class and element matrix inner product sens
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for element_handle in self.elements:
                    with self.subTest(element=element_handle):
                        if element_handle not in self.thermal_elements:
                            element = element_handle(transform, self.con)
                            for matrix_type in self.matrix_types:
                                with self.subTest(matrix_type=matrix_type):
                                    fail = elements.TestElementMatXptSens(
                                        element,
                                        matrix_type,
                                        self.elem_index,
                                        self.time,
                                        self.xpts,
                                        self.vars,
                                        self.dh,
                                        self.print_level,
                                        self.atol,
                                        self.rtol,
                                    )
                                    self.assertFalse(fail)

    def test_element_mat_sv_sens(self):
        # Loop through every combination of model and basis class and test element matrix inner product sens
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for element_handle in self.elements:
                    with self.subTest(element=element_handle):
                        if element_handle not in self.thermal_elements:
                            element = element_handle(transform, self.con)
                            for matrix_type in self.matrix_types:
                                with self.subTest(matrix_type=matrix_type):
                                    fail = elements.TestElementMatSVSens(
                                        element,
                                        matrix_type,
                                        self.elem_index,
                                        self.time,
                                        self.xpts,
                                        self.vars,
                                        self.dh,
                                        self.print_level,
                                        self.atol,
                                        self.rtol,
                                    )
                                    self.assertFalse(fail)
