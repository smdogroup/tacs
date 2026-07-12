import unittest

import numpy as np

from tacs import TACS, constitutive, elements


class ElementTest(unittest.TestCase):
    def setUp(self):
        max_nodes = 64
        max_vars_per_nodes = 8
        max_vars = max_nodes * max_vars_per_nodes

        # fd/cs step size
        if TACS.dtype is complex:
            self.dh = 1e-200
            self.rtol = 1e-10
        else:
            self.dh = 1e-5
            self.rtol = 1e-2
        self.dtype = TACS.dtype

        self.atol = np.clip(1e-5 * self.rtol, 1e-14, 1e-8)
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
            alpha=cte,
            kappa=kappa,
        )

        ref_axis = np.array([0.0, 1.0, 1.0], dtype=self.dtype)
        self.transforms = [elements.BeamRefAxisTransform(ref_axis)]

        # TACS beam elements of various orders and types
        self.elements = [elements.Beam2, elements.Beam3]

        t = 0.1
        d = 1.0
        dNum = 0
        tNum = 1

        # Create stiffness (need class)
        self.con = constitutive.IsoTubeBeamConstitutive(
            self.props, t=t, tNum=tNum, d=d, dNum=dNum
        )

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
        # Loop through every combination of transform type and beam element class and test residual
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for element_handle in self.elements:
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
        # Loop through every combination of transform type and beam element class and test Jacobian
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
        # Loop through every combination of transform type and beam element class and test adjoint residual-dvsens product
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
        # Loop through every combination of transform type and beam element class and test adjoint residual-xptsens product
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

    def test_quantity_xpt_sens(self):
        # Loop through every combination of transform type and beam element
        # class and test the point-quantity Xpts-sensitivity (SPEC.md sec
        # 3.1, Phase 3 Task 3.2). TACS_ELEMENT_MOMENT_OF_INERTIA is the only
        # quantity type this feature brings in scope for
        # addPointQuantityXptSens; the other quantity types
        # (TACS_ELEMENT_DENSITY_MOMENT, TACS_FAILURE_INDEX,
        # TACS_STRAIN_ENERGY_DENSITY) were already analytic before this
        # feature and are not re-verified here.
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for element_handle in self.elements:
                    with self.subTest(element=element_handle):
                        element = element_handle(transform, self.con)
                        fail = elements.TestElementQuantityXptSens(
                            element,
                            self.elem_index,
                            elements.ELEMENT_MOMENT_OF_INERTIA,
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
        # Loop through every combination of transform type and beam element class and element matrix inner product sens
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

    def test_element_mat_sv_sens_nonlinear_director(self):
        # Task 7.3 (SPEC-phase-7.md sec 4.1/6.3): new test-surface work
        # bringing a nonlinear-director beam typedef (ModRot, i.e.
        # TACSQuadraticRotation) into this file's matrix-sens coverage --
        # a prerequisite for G1/G2's own oracle to exist at all, per
        # SPEC-phase-7.md sec 6.3. A separate, narrow element list (not
        # self.elements) since test_element_jacobian/other tests in this
        # file are NOT yet analytic for this director class (Phase 2 scope)
        # and adding ModRot to the shared self.elements list would break
        # them. No TACSBeam*Quaternion typedef exists in this codebase
        # (confirmed by grep) -- TACSQuaternionRotation's addDirectorHessian*
        # hooks are covered at the director-hook level instead
        # (src/elements/shell/tests/test_director_hessian_second_order.cpp).
        #
        # CURRENTLY TRIVIALLY-PASSING FD-vs-FD for TACS_STIFFNESS_MATRIX/
        # TACS_MASS_MATRIX specifically (both still forward to the base
        # FD/CS fallback for this director class -- see
        # TACSBeamElement::getMatSVSensInnerProduct's own documented-failure
        # comment, Task 7.3: a genuine attempt found the third-derivative
        # recipe SPEC-phase-7.md sec 4.1 describes is missing at least one
        # term, discovered via this exact test). GEOMETRIC_STIFFNESS_MATRIX
        # is genuinely analytic here (Task 5.4, TACSLinearizedRotation-only,
        # unaffected by ModRot -- this subtest exercises ModRot's own
        # base-FD-fallback path for that matType too, still a legitimate,
        # if not novel, regression check). This test's own value: it is the
        # regression net that will make TACS_STIFFNESS_MATRIX/TACS_MASS_MATRIX
        # genuine the moment a future session resolves the documented gap.
        modrot_elements = [elements.Beam2ModRot, elements.Beam3ModRot]
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for element_handle in modrot_elements:
                    with self.subTest(element=element_handle):
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

    def test_element_mat_sv_sens(self):
        # Loop through every combination of model and basis class and test element matrix inner product sens
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for element_handle in self.elements:
                    with self.subTest(element=element_handle):
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
