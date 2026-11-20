"""
Check the adjoint-residual product of every auxiliary load element that has
load design variables.

Each load element is attached to a single-element assembler whose state
variables are zero, so the regular element contributes nothing to the residual
and the residual comes from the load alone. For every design variable, the
entry of Assembler.addAdjointResProducts is compared against an independent
derivative of psi^T R obtained from Assembler.assembleRes:

- complex builds use a complex step, so the comparison is exact to round-off,
- real builds use a central difference with a step relative to the DV value,
  which is accurate to round-off for these residuals (they are at most
  quadratic in the load DVs and polynomial in the constitutive DVs).

Load magnitudes are deliberately large (10 MPa pressures, 1e6 m/s^2 accelerations,
1e4 rad/s spin rates) because a fixed absolute finite-difference step loses
accuracy as the DV value grows.
"""

import unittest

import numpy as np
from mpi4py import MPI

from tacs import TACS, constitutive, elements

DTYPE = TACS.dtype
IS_COMPLEX = DTYPE is complex
COMPLEX_STEP = 1e-30
RELATIVE_FD_STEP = 1e-4
RTOL = 1e-12 if IS_COMPLEX else 1e-7

# 2D quad, 3D hexa, shell quad and beam node coordinates (x, y, z per node)
QUAD_XPTS = [0.0, 0.0, 0.0, 1.1, 0.1, 0.0, -0.1, 0.9, 0.0, 1.0, 1.2, 0.0]
HEXA_XPTS = [
    0.0, 0.0, 0.0,
    1.1, 0.1, -0.1,
    0.1, 0.9, 0.05,
    1.0, 1.1, 0.0,
    -0.05, 0.1, 1.0,
    0.9, 0.0, 1.1,
    0.0, 1.05, 0.95,
    1.1, 0.95, 1.05,
]  # fmt: skip
SHELL_XPTS = [0.0, 0.0, 0.0, 1.1, 0.1, 0.05, -0.1, 0.9, -0.05, 1.0, 1.2, 0.1]
BEAM_XPTS = [0.0, 0.0, 0.0, 1.0, 0.2, 0.1]


def createMaterial():
    """Return an isotropic aluminium material."""
    return constitutive.MaterialProperties(rho=2700.0, E=70e9, nu=0.3, ys=270e6)


def createElement2D():
    """Return a linear quad plane-stress element with a thickness DV (DV 0)."""
    con = constitutive.PlaneStressConstitutive(createMaterial(), t=0.01, tNum=0)
    model = elements.LinearElasticity2D(con)
    return elements.Element2D(model, elements.LinearQuadBasis())


def createElement3D():
    """Return a linear hexa solid element with a 'thickness' DV (DV 0)."""
    con = constitutive.SolidConstitutive(createMaterial(), t=0.7, tNum=0)
    model = elements.LinearElasticity3D(con)
    return elements.Element3D(model, elements.LinearHexaBasis())


def createShell():
    """Return a 4-node shell element with a thickness DV (DV 0)."""
    con = constitutive.IsoShellConstitutive(createMaterial(), t=0.01, tNum=0)
    return elements.Quad4Shell(elements.ShellNaturalTransform(), con)


def createBeam():
    """Return a 2-node beam element with thickness (DV 0) and width (DV 1) DVs."""
    con = constitutive.IsoRectangleBeamConstitutive(
        createMaterial(), t=0.01, w=0.05, tNum=0, wNum=1
    )
    transform = elements.BeamRefAxisTransform(np.array([0.0, 1.0, 1.0], dtype=DTYPE))
    return elements.Beam2(transform, con)


def asScalars(values):
    """Convert a list of values to an array of the TACS scalar type."""
    return np.array(values, dtype=DTYPE)


def asDVNums(values):
    """Convert a list of DV numbers to the integer type TACS expects."""
    return np.array(values, dtype=np.intc)


OMEGA = [1.0e4, -2.0e4, 1.5e4]
ROT_CENTER = [0.3, -0.2, 0.5]

# (name, element factory, node coordinates, function creating the load element)
LOAD_CASES = [
    (
        "Traction2D",
        createElement2D,
        QUAD_XPTS,
        lambda e: e.createElementTraction(0, asScalars([3e6, -2e6]), asDVNums([1, 2])),
    ),
    (
        "Pressure2D",
        createElement2D,
        QUAD_XPTS,
        lambda e: e.createElementPressure(0, 1e7, 1),
    ),
    (
        "InertialForce2D",
        createElement2D,
        QUAD_XPTS,
        lambda e: e.createElementInertialForce(
            asScalars([1e6, -2e6]), asDVNums([1, 2])
        ),
    ),
    (
        "Traction3D",
        createElement3D,
        HEXA_XPTS,
        lambda e: e.createElementTraction(
            0, asScalars([3e6, -2e6, 1e6]), asDVNums([1, 2, 3])
        ),
    ),
    (
        "Pressure3D",
        createElement3D,
        HEXA_XPTS,
        lambda e: e.createElementPressure(0, 1e7, 1),
    ),
    (
        "InertialForce3D",
        createElement3D,
        HEXA_XPTS,
        lambda e: e.createElementInertialForce(
            asScalars([1e6, -2e6, 5e5]), asDVNums([1, 2, 3])
        ),
    ),
    (
        "CentrifugalForce3D",
        createElement3D,
        HEXA_XPTS,
        lambda e: e.createElementCentrifugalForce(
            asScalars(OMEGA),
            asScalars(ROT_CENTER),
            omegaDVNums=asDVNums([1, 2, 3]),
            rotCenterDVNums=asDVNums([4, 5, 6]),
        ),
    ),
    (
        "ShellTraction",
        createShell,
        SHELL_XPTS,
        lambda e: e.createElementTraction(
            0, asScalars([3e6, -2e6, 1e6]), asDVNums([1, 2, 3])
        ),
    ),
    (
        "ShellPressure",
        createShell,
        SHELL_XPTS,
        lambda e: e.createElementPressure(0, 1e7, 1),
    ),
    (
        "ShellInertialForce",
        createShell,
        SHELL_XPTS,
        lambda e: e.createElementInertialForce(
            asScalars([1e6, -2e6, 5e5]), asDVNums([1, 2, 3])
        ),
    ),
    (
        "ShellCentrifugalForce",
        createShell,
        SHELL_XPTS,
        lambda e: e.createElementCentrifugalForce(
            asScalars(OMEGA),
            asScalars(ROT_CENTER),
            omegaDVNums=asDVNums([1, 2, 3]),
            rotCenterDVNums=asDVNums([4, 5, 6]),
        ),
    ),
    (
        "BeamTraction",
        createBeam,
        BEAM_XPTS,
        lambda e: e.createElementTraction(
            0, asScalars([3e5, -2e5, 1e5]), asDVNums([2, 3, 4])
        ),
    ),
    (
        "BeamInertialForce",
        createBeam,
        BEAM_XPTS,
        lambda e: e.createElementInertialForce(
            asScalars([1e6, -2e6, 5e5]), asDVNums([2, 3, 4])
        ),
    ),
    (
        "BeamCentrifugalForce",
        createBeam,
        BEAM_XPTS,
        lambda e: e.createElementCentrifugalForce(
            asScalars(OMEGA),
            asScalars(ROT_CENTER),
            omegaDVNums=asDVNums([2, 3, 4]),
            rotCenterDVNums=asDVNums([5, 6, 7]),
        ),
    ),
]


def createAssembler(elem, loadElem, xpts):
    """Return an initialized single-element assembler with loadElem attached as an aux element."""
    numNodes = elem.getNumNodes()
    assembler = TACS.Assembler.create(MPI.COMM_SELF, elem.getVarsPerNode(), numNodes, 1)
    assembler.setElementConnectivity(
        asDVNums([0, numNodes]), np.arange(numNodes, dtype=np.intc)
    )
    assembler.setElements([elem])
    auxElems = TACS.AuxElements()
    auxElems.addElement(0, loadElem)
    # Attach before initialize() so the load DVs are counted in the design vector
    assembler.setAuxElements(auxElems)
    assembler.initialize()
    nodes = assembler.createNodeVec()
    nodes.getArray()[:] = asScalars(xpts)
    assembler.setNodes(nodes)
    return assembler


def computeResidualProduct(assembler, adjoint, x):
    """Return adjoint^T R(x) for the design variable values x."""
    xVec = assembler.createDesignVec()
    xVec.getArray()[:] = x
    assembler.setDesignVars(xVec)
    res = assembler.createVec()
    assembler.assembleRes(res, applyBCs=False)
    return np.sum(adjoint * res.getArray())


def computeReferenceProduct(assembler, adjoint, x):
    """Return d(adjoint^T R)/dx for every DV via complex step or a relative central difference."""
    ref = np.zeros(len(x))
    for k in range(len(x)):
        if IS_COMPLEX:
            xPert = x.copy()
            xPert[k] += 1j * COMPLEX_STEP
            ref[k] = (
                np.imag(computeResidualProduct(assembler, adjoint, xPert))
                / COMPLEX_STEP
            )
        else:
            step = RELATIVE_FD_STEP * abs(x[k])
            xPlus = x.copy()
            xPlus[k] += step
            xMinus = x.copy()
            xMinus[k] -= step
            ref[k] = (
                computeResidualProduct(assembler, adjoint, xPlus)
                - computeResidualProduct(assembler, adjoint, xMinus)
            ) / (2.0 * step)
    # Restore the unperturbed design variables
    computeResidualProduct(assembler, adjoint, x)
    return ref


class AuxElementDVAdjointTest(unittest.TestCase):
    """Adjoint-residual products of auxiliary load elements w.r.t. their DVs."""

    def test_adjoint_residual_product_matches_reference(self):
        """
        Given a single-element assembler with one auxiliary load element whose load components are DVs,
        when Assembler.addAdjointResProducts evaluates the adjoint-residual product for a random adjoint,
        then every DV entry matches the complex-step (complex builds) or relative central-difference (real builds)
        derivative of adjoint^T R.
        """
        rng = np.random.default_rng(0)
        for name, createBaseElement, xpts, createLoadElement in LOAD_CASES:
            with self.subTest(load=name):
                elem = createBaseElement()
                assembler = createAssembler(elem, createLoadElement(elem), xpts)
                x0Vec = assembler.createDesignVec()
                assembler.getDesignVars(x0Vec)
                x0 = x0Vec.getArray().copy()

                adjointVec = assembler.createVec()
                adjointVec.getArray()[:] = rng.uniform(
                    -1.0, 1.0, adjointVec.getArray().shape
                )
                adjoint = adjointVec.getArray().copy()

                dfdx = assembler.createDesignVec()
                assembler.addAdjointResProducts([adjointVec], [dfdx])
                analytic = np.real(dfdx.getArray())

                ref = computeReferenceProduct(assembler, adjoint, x0)
                # Entries differ by orders of magnitude between DVs, so compare each one relatively;
                # DVs the residual does not depend on are exactly zero in both
                np.testing.assert_allclose(analytic, ref, rtol=RTOL, atol=0.0)


if __name__ == "__main__":
    unittest.main()
