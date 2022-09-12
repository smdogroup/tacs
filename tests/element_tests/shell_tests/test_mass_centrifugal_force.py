from tacs import TACS, elements, constitutive
import numpy as np
import unittest


class ElementTest(unittest.TestCase):
    def setUp(self):
        num_nodes = 1
        vars_per_nodes = 6
        num_vars = num_nodes * vars_per_nodes

        # fd/cs step size
        if TACS.dtype is complex:
            self.dh = 1e-50
            self.rtol = 1e-10
        else:
            self.dh = 1e-6
            self.rtol = 1e-2
        self.dtype = TACS.dtype

        # Basically, only check relative tolerance
        self.atol = 1e99
        self.print_level = 0

        self.xpts = np.zeros(3, dtype=self.dtype)

        # Set element index
        self.elem_index = 0
        # Set the simulation time
        self.time = 0.0

        # Set the variable arrays
        np.random.seed(30)  # Seed random numbers for deterministic/repeatable tests
        self.vars = np.random.rand(num_vars).astype(self.dtype)
        self.dvars = self.vars.copy()
        self.ddvars = self.vars.copy()

        # General 6 dof mass matrix
        M = np.arange(21, dtype=self.dtype)
        # General 6 dof mass matrix for point mass
        m = 2.0
        I11 = I22 = I33 = 1.0
        I12 = I13 = I23 = 0.5

        # Rot velocity vector
        self.omega = np.array([1.0, 2.0, 3.0], dtype=self.dtype)
        self.rotCenter = np.array([4.0, -5.0, 6.0], dtype=self.dtype)

        # Create constitutive classes
        self.con_objects = [
            constitutive.GeneralMassConstitutive(M=M),
            constitutive.PointMassConstitutive(
                m=m, I11=I11, I22=I22, I33=I33, I12=I12, I13=I13, I23=I23
            ),
        ]

        # Set matrix types
        self.matrix_types = [
            TACS.STIFFNESS_MATRIX,
            TACS.MASS_MATRIX,
            TACS.GEOMETRIC_STIFFNESS_MATRIX,
        ]

        # Seed random number generator in tacs for consistent test results
        elements.SeedRandomGenerator(0)

    def test_element_jacobian(self):
        # Loop through each combination of consstitutive object and test Jacobian
        for con in self.con_objects:
            with self.subTest(con=con.getObjectName()):
                element = elements.MassElement(con)
                force = element.createElementCentrifugalForce(
                    self.omega, self.rotCenter
                )
                fail = elements.TestElementJacobian(
                    force,
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
        # Loop through each combination of consstitutive object and test adjoint residual-dvsens product
        for con in self.con_objects:
            with self.subTest(con=con.getObjectName()):
                element = elements.MassElement(con)
                force = element.createElementCentrifugalForce(
                    self.omega, self.rotCenter
                )
                dvs = force.getDesignVars(self.elem_index)
                fail = elements.TestAdjResProduct(
                    force,
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
        # Loop through each combination of constitutive object and test adjoint residual-xptsens product
        for con in self.con_objects:
            with self.subTest(con=con.getObjectName()):
                element = elements.MassElement(con)
                force = element.createElementCentrifugalForce(
                    self.omega, self.rotCenter
                )
                fail = elements.TestAdjResXptProduct(
                    force,
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
        # Loop through each combination of consstitutive object and test element matrix inner product sens
        for con in self.con_objects:
            with self.subTest(con=con.getObjectName()):
                element = elements.MassElement(con)
                force = element.createElementCentrifugalForce(
                    self.omega, self.rotCenter
                )
                dvs = force.getDesignVars(self.elem_index)
                for matrix_type in self.matrix_types:
                    with self.subTest(matrix_type=matrix_type):
                        fail = elements.TestElementMatDVSens(
                            force,
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

    def test_element_mat_sv_sens(self):
        # Loop through each combination of constitutive object and test element matrix inner product sens
        for con in self.con_objects:
            with self.subTest(con=con.getObjectName()):
                element = elements.MassElement(con)
                force = element.createElementCentrifugalForce(
                    self.omega, self.rotCenter
                )
                fail = elements.TestElementMatSVSens(
                    force,
                    TACS.GEOMETRIC_STIFFNESS_MATRIX,
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
