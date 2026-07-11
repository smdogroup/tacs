import unittest

import numpy as np

from tacs import TACS, constitutive, elements


def generateTestFailMessage(element, transform, matType=None):
    message = f"Failed test with element {element.__class__.__name__}, transform {transform.__class__.__name__}"
    if matType is not None:
        if matType == TACS.STIFFNESS_MATRIX:
            matName = "stiffness"
        if matType == TACS.MASS_MATRIX:
            matName = "mass"
        if matType == TACS.GEOMETRIC_STIFFNESS_MATRIX:
            matName = "geometric stiffness"
        message += f", {matName} matrix"
    return message


class ElementTest(unittest.TestCase):
    def setUp(self):
        max_nodes = 64
        max_vars_per_nodes = 8
        max_vars = max_nodes * max_vars_per_nodes

        # fd/cs step size
        if TACS.dtype is complex:
            self.dh = 1e-30
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
            elements.Quad4NonlinearShellModRot,
            elements.Quad9NonlinearShellModRot,
            elements.Quad16NonlinearShellModRot,
            elements.Tri3NonlinearShellModRot,
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

        # Set quantity types
        self.quantity_types = [
            elements.FAILURE_INDEX,
            elements.STRAIN_ENERGY_DENSITY,
            elements.ELEMENT_DENSITY,
            elements.ELEMENT_DISPLACEMENT,
            elements.ELEMENT_DENSITY_MOMENT,
            elements.ELEMENT_MOMENT_OF_INERTIA,
            elements.ELEMENT_ENCLOSED_VOLUME,
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
                    if element_handle not in self.thermal_elements:
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
                            self.assertFalse(
                                fail,
                                msg=generateTestFailMessage(
                                    element=element, transform=transform
                                ),
                            )

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
                        self.assertFalse(
                            fail,
                            msg=generateTestFailMessage(
                                element=element, transform=transform
                            ),
                        )

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
                        self.assertFalse(
                            fail,
                            msg=generateTestFailMessage(
                                element=element, transform=transform
                            ),
                        )

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
                        self.assertFalse(
                            fail,
                            msg=generateTestFailMessage(
                                element=element, transform=transform
                            ),
                        )

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
                                self.assertFalse(
                                    fail,
                                    msg=generateTestFailMessage(
                                        element=element, transform=transform
                                    ),
                                )

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
                                    self.assertFalse(
                                        fail,
                                        msg=generateTestFailMessage(
                                            element=element,
                                            transform=transform,
                                            matType=matrix_type,
                                        ),
                                    )

    def test_point_quantity_xpt_sens(self):
        # Here we use a tighter dh than self.dh (real mode only; complex mode
        # is unaffected since it always uses a fixed tiny complex-step size).
        # TACS_ELEMENT_ENCLOSED_VOLUME's quantity is (X.n0)/3, which can land
        # near zero for some element/random-seed/Xpts-DOF combinations (e.g.
        # Quad9Shell + ShellRefAxisTransform); at self.dh=1e-5 the harness's
        # plain forward-difference reference has enough truncation error at
        # that near-zero component to exceed self.rtol, even though the
        # analytic value there agrees with the base class's own (much
        # smaller-step) internal FD to ~6 significant figures. A tighter dh
        # removes the truncation artifact without touching atol/rtol.
        dh = self.dh if TACS.dtype is complex else 1e-6
        # Loop through every combination of transform type and shell element class and
        # test the point quantity output node coordinate sensitivities
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for element_handle in self.elements:
                    with self.subTest(element=element_handle):
                        element = element_handle(transform, self.con)
                        for quantity_type in self.quantity_types:
                            with self.subTest(quantity_type=quantity_type):
                                fail = elements.TestElementQuantityXptSens(
                                    element,
                                    self.elem_index,
                                    quantity_type,
                                    self.time,
                                    self.xpts,
                                    self.vars,
                                    self.dvars,
                                    self.ddvars,
                                    dh,
                                    self.print_level,
                                    self.atol,
                                    self.rtol,
                                )
                                self.assertFalse(
                                    fail,
                                    msg=generateTestFailMessage(
                                        element=element, transform=transform
                                    ),
                                )

    def test_point_quantity_xpt_sens_degenerate_ref_axis(self):
        # Regression test for a NaN in ShellRefAxisTransform::addTransformSens
        # (TACSShellElementTransform.h) that test_point_quantity_xpt_sens above
        # never exercised: that test's fixed ref_axis=[0, 1, 1] essentially
        # never lands bit-exactly parallel to a random element's normal, so it
        # never hit the t1 == [0, 0, 0] degenerate branch (ref axis exactly
        # parallel to the shell normal -- the case
        # TACSShellRefAxisTransform.computeTransform's own "ill-conditioned"
        # warning anticipates). In that branch, computing t1n = t1/|t1|
        # produced 0.0*inf = NaN, and that NaN then survived being multiplied
        # by an exactly-zero incoming seed (NaN*0.0 = NaN, not 0.0 in IEEE
        # 754), corrupting the Xpts sensitivity of every quantity type whose
        # addPointQuantityXptSens branch calls addTransformSens with a
        # zeroed dT (TACS_ELEMENT_DENSITY/DISPLACEMENT/DENSITY_MOMENT/
        # ENCLOSED_VOLUME) even though the correct contribution is exactly
        # zero. FAILURE_INDEX/STRAIN_ENERGY_DENSITY are excluded here since
        # their forward evalPointQuantity itself calls computeTransform,
        # which has its own (separate, pre-existing, unguarded) singularity
        # at this exact ref-axis/normal alignment -- out of scope for this
        # regression.
        #
        # Flattening every node's z-coordinate to exactly 0.0 forces every
        # element (any order/type) to have its normal land bit-exactly on
        # (0, 0, 1) or (0, 0, -1): the interpolated Xxi is a linear
        # combination of node coordinates, so an all-zero z-input produces an
        # exactly-zero z-component in Xxi, and hence in the normal from
        # crossProduct. Pairing that with ref_axis=[0, 0, 1] reproduces the
        # exact degenerate alignment.
        flat_xpts = self.xpts.copy()
        flat_xpts[2::3] = 0.0

        degenerate_transform = elements.ShellRefAxisTransform(
            np.array([0.0, 0.0, 1.0], dtype=self.dtype)
        )
        geometry_only_quantity_types = [
            elements.ELEMENT_DENSITY,
            elements.ELEMENT_DISPLACEMENT,
            elements.ELEMENT_DENSITY_MOMENT,
            elements.ELEMENT_ENCLOSED_VOLUME,
        ]
        # At this exact flat/degenerate configuration, the true Xpts
        # sensitivity of detXd (and hence of these quantities) w.r.t. a
        # node's out-of-plane (z) coordinate is analytically exactly zero to
        # first order (bending a flat sheet changes its area only at second
        # order). A plain forward difference at self.dh=1e-6 (or even 1e-8)
        # leaks enough second-order truncation error into that near-zero
        # component to exceed self.atol on its own (same class of artifact
        # already noted in test_point_quantity_xpt_sens above for
        # ENCLOSED_VOLUME, just sharper here since the zero is exact rather
        # than merely small); a tighter dh removes it without touching
        # atol/rtol. dh=1e-9 was verified empirically to clear self.atol for
        # every element/quantity-type combination below without yet being
        # small enough for floating-point roundoff to dominate.
        dh = self.dh if TACS.dtype is complex else 1e-9
        for element_handle in self.elements:
            with self.subTest(element=element_handle):
                element = element_handle(degenerate_transform, self.con)
                for quantity_type in geometry_only_quantity_types:
                    with self.subTest(quantity_type=quantity_type):
                        fail = elements.TestElementQuantityXptSens(
                            element,
                            self.elem_index,
                            quantity_type,
                            self.time,
                            flat_xpts,
                            self.vars,
                            self.dvars,
                            self.ddvars,
                            dh,
                            self.print_level,
                            self.atol,
                            self.rtol,
                        )
                        self.assertFalse(
                            fail,
                            msg=generateTestFailMessage(
                                element=element, transform=degenerate_transform
                            ),
                        )

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
                                    self.assertFalse(
                                        fail,
                                        msg=generateTestFailMessage(
                                            element=element,
                                            transform=transform,
                                            matType=matrix_type,
                                        ),
                                    )


if __name__ == "__main__":
    unittest.main()
