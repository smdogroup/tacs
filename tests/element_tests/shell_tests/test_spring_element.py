from tacs import TACS, constitutive, elements
import numpy as np
import unittest


class ElementTest(unittest.TestCase):
    def setUp(self):
        num_nodes = 2
        vars_per_node = 6
        num_vars = num_nodes * vars_per_node

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
        self.xpts = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0], dtype=self.dtype)
        np.random.seed(30)  # Seed random numbers for deterministic/repeatable tests
        self.vars = np.random.rand(num_vars).astype(self.dtype)
        self.dvars = self.vars.copy()
        self.ddvars = self.vars.copy()

        ref_axis1 = np.array([0.0, 1.0, 0.0], dtype=self.dtype)
        ref_axis2 = np.array([1.0, 0.0, 0.0], dtype=self.dtype)
        self.transforms = [
            elements.SpringIdentityTransform(),
            elements.SpringRefAxisTransform(ref_axis1),
            elements.SpringRefFrameTransform(ref_axis1, ref_axis2),
        ]

        # Create stiffness (need class)
        K = np.ones(21, dtype=self.dtype)
        k = np.arange(1, 7, dtype=self.dtype)
        self.con_objects = [
            constitutive.GeneralSpringConstitutive(K=K),
            constitutive.DOFSpringConstitutive(k=k),
        ]

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
        # Loop through every combination of transform type and spring constitutive class and test residual
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for con in self.con_objects:
                    with self.subTest(con=con):
                        element = elements.SpringElement(transform, con)
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
        # Loop through every combination of transform type and spring constitutive class and test Jacobian
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for con in self.con_objects:
                    with self.subTest(con=con):
                        element = elements.SpringElement(transform, con)
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
        # Loop through every combination of transform type and spring constitutive class and test adjoint residual-dvsens product
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for con in self.con_objects:
                    with self.subTest(con=con):
                        element = elements.SpringElement(transform, con)
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
        # Loop through every combination of transform type and spring constitutive class and test adjoint residual-xptsens product
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for con in self.con_objects:
                    with self.subTest(con=con):
                        element = elements.SpringElement(transform, con)
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
        # Loop through every combination of transform type and spring constitutive class and element matrix inner product sens
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for con in self.con_objects:
                    with self.subTest(con=con):
                        element = elements.SpringElement(transform, con)
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
        # Loop through every combination of transform type and spring constitutive class and element matrix inner product sens
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for con in self.con_objects:
                    with self.subTest(con=con):
                        element = elements.SpringElement(transform, con)
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
        # Loop through every combination of transform type and spring constitutive class
        # and test element matrix inner product sens
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for con in self.con_objects:
                    with self.subTest(con=con):
                        for matrix_type in self.matrix_types:
                            with self.subTest(matrix_type=matrix_type):
                                element = elements.SpringElement(transform, con)
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
