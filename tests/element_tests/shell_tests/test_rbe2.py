from tacs import TACS, elements
import numpy as np
import unittest


class ElementTest(unittest.TestCase):
    def setUp(self):
        num_indep_nodes = 1  # Always 1 for rbe2
        num_dep_nodes = 10
        num_dummy_nodes = num_dep_nodes
        self.num_nodes = num_dep_nodes + num_indep_nodes + num_dummy_nodes

        vars_per_nodes = 6
        num_vars = self.num_nodes * vars_per_nodes

        # Set artificial stiffness term to zero,
        # since we know it adds error to the sensitivities
        self.C1 = 1e6
        self.C2 = 0.0

        # fd/cs step size
        if TACS.dtype is complex:
            self.dh = 1e-50
            self.rtol = 1e-10
        else:
            self.dh = 1e-6
            self.rtol = 1e-2
        self.dtype = TACS.dtype

        dep_xpts = np.zeros([num_dep_nodes, 3], self.dtype)
        dep_xpts[:, 0] = np.linspace(-1.0, 1.0, num_dep_nodes)
        dep_xpts[:5, 1] = 1.0
        dep_xpts[5:, 1] = -1.0
        dep_xpts[:5, 2] = 1.0
        dep_xpts[5:, 2] = -1.0

        indep_xpts = np.zeros([num_indep_nodes, 3], self.dtype)
        dummy_xpts = np.zeros([num_dummy_nodes, 3], self.dtype)

        xpts = np.append(indep_xpts, dep_xpts)
        self.xpts = np.append(xpts, dummy_xpts)

        # Basically, only check relative tolerance
        self.atol = 1e99
        self.print_level = 0

        # Set element index
        self.elem_index = 0
        # Set the simulation time
        self.time = 0.0

        # Set the variable arrays
        np.random.seed(30)  # Seed random numbers for deterministic/repeatable tests
        self.vars = np.random.rand(num_vars).astype(self.dtype)
        self.dvars = self.vars.copy()
        self.ddvars = self.vars.copy()

        # Specify dofs for dependent nodes
        self.dep_dofs_constrained = [
            np.array([1, 1, 1, 1, 1, 1], np.intc),
            np.array([1, 1, 1, 0, 0, 0], np.intc),
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
        # Loop through each combination of dof constraints and test Jacobian
        for dep_dofs in self.dep_dofs_constrained:
            with self.subTest(dep_dofs=dep_dofs):
                element = elements.RBE2(self.num_nodes, dep_dofs, self.C1, self.C2)
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
        # Loop through each combination of dof constraints and test adjoint residual-dvsens product
        for dep_dofs in self.dep_dofs_constrained:
            with self.subTest(dep_dofs=dep_dofs):
                element = elements.RBE2(self.num_nodes, dep_dofs, self.C1, self.C2)
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
        # Loop through each combination of dof constraints and test adjoint residual-xptsens product
        for dep_dofs in self.dep_dofs_constrained:
            with self.subTest(dep_dofs=dep_dofs):
                element = elements.RBE2(self.num_nodes, dep_dofs, self.C1, self.C2)
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
        # Loop through each combination of dof constraints and element matrix inner product sens
        for dep_dofs in self.dep_dofs_constrained:
            with self.subTest(dep_dofs=dep_dofs):
                element = elements.RBE2(self.num_nodes, dep_dofs, self.C1, self.C2)
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
        # Loop through each combination of dof constraints and element matrix inner product sens
        for dep_dofs in self.dep_dofs_constrained:
            with self.subTest(dep_dofs=dep_dofs):
                element = elements.RBE2(self.num_nodes, dep_dofs, self.C1, self.C2)
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
        # Loop through each combination of dof constraints and test element matrix inner product sens
        for dep_dofs in self.dep_dofs_constrained:
            with self.subTest(dep_dofs=dep_dofs):
                element = elements.RBE2(self.num_nodes, dep_dofs, self.C1, self.C2)
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
