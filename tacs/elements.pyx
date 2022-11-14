#  This file is part of TACS: The Toolkit for the Analysis of Composite
#  Structures, a parallel finite-element code for structural and
#  multidisciplinary design optimization.
#
#  Copyright (C) 2014 Georgia Tech Research Corporation
#
#  TACS is licensed under the Apache License, Version 2.0 (the
#  "License"); you may not use this software except in compliance with
#  the License.  You may obtain a copy of the License at
#
#  http://www.apache.org/licenses/LICENSE-2.0

#distutils: language=c++

# For the use of MPI
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy
cimport numpy as np
import numpy as np

# Ensure that numpy is initialized
np.import_array()

# Import the definition required for const strings
from libc.string cimport const_char
from libc.stdlib cimport malloc, free

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

# Include the definitions
include "TacsDefs.pxi"

# Import the definitions
from TACS cimport *
from constitutive cimport *
from elements cimport *

# Include the mpi4py header
cdef extern from "mpi-compat.h":
    pass

# TACSElement quantity types
ELEMENT_DENSITY = TACS_ELEMENT_DENSITY
STRAIN_ENERGY_DENSITY = TACS_STRAIN_ENERGY_DENSITY
FAILURE_INDEX = TACS_FAILURE_INDEX
HEAT_FLUX = TACS_HEAT_FLUX
TEMPERATURE = TACS_TEMPERATURE
TOTAL_STRAIN_ENERGY_DENSITY = TACS_TOTAL_STRAIN_ENERGY_DENSITY
ELEMENT_DISPLACEMENT = TACS_ELEMENT_DISPLACEMENT
ELEMENT_STRAIN = TACS_ELEMENT_STRAIN
ELEMENT_STRESS = TACS_ELEMENT_STRESS

# Flags for the thermomechanical model
STEADY_STATE_MECHANICAL = TACS_STEADY_STATE_MECHANICAL
STEADY_STATE_THERMAL = TACS_STEADY_STATE_THERMAL

def TestElementBasis(ElementBasis basis, double dh=1e-6,
                     int test_print_level=2, double atol=1e-30,
                     double rtol=1e-5):
    return TacsTestElementBasis(basis.ptr, dh, test_print_level, atol, rtol)

def TestElementBasisFunctions(ElementBasis basis, double dh=1e-6,
                                  int test_print_level=2, double atol=1e-5,
                                  double rtol=1e-5):
    return TacsTestElementBasisFunctions(basis.ptr, dh, test_print_level, atol, rtol)

def TestElementBasisFaceNormals(ElementBasis basis, double dh=1e-6,
                                    int test_print_level=2, double atol=1e-5,
                                    double rtol=1e-5):
    return TacsTestElementBasisFaceNormals(basis.ptr, dh, test_print_level, atol, rtol)

def TestElementBasisJacobianTransform(ElementBasis basis, double dh=1e-6,
                                      int test_print_level=2, double atol=1e-5,
                                      double rtol=1e-5):
    return TacsTestElementBasisJacobianTransform(basis.ptr, dh, test_print_level, atol, rtol)

def TestElementModelJacobian(ElementModel model, int elem_index, double time, double dh=1e-6,
                             int test_print_level=2, double atol=1e-5,
                             double rtol=1e-5):
    return TacsTestElementModelJacobian(model.ptr, elem_index, time, dh, test_print_level, atol, rtol)

def TestElementModelAdjXptSensProduct(ElementModel model, int elem_index, double time, double dh=1e-6,
                                      int test_print_level=2, double atol=1e-5,
                                      double rtol=1e-5):
    return TacsTestElementModelAdjXptSensProduct(model.ptr, elem_index, time, dh, test_print_level, atol, rtol)

def TestElementResidual(Element element, int elem_index,
                        double time,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] xpts,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] vars,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] dvars,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] ddvars,
                        double dh=1e-6,
                        int test_print_level=2,
                        double atol=1e-5, double rtol=1e-5):
    num_nodes = element.getNumNodes()
    num_vars = element.getNumVariables()

    # Make sure input arrays are large enough for element to avoid segfault
    assert len(xpts) >= 3 * num_nodes
    assert len(vars) >= num_vars
    assert len(dvars) >= num_vars
    assert len(ddvars) >= num_vars

    return TacsTestElementResidual(element.ptr, elem_index, time,
                                   <TacsScalar*>xpts.data,
                                   <TacsScalar*>vars.data,
                                   <TacsScalar*>dvars.data,
                                   <TacsScalar*>ddvars.data,
                                   dh,
                                   test_print_level, atol, rtol)

def TestElementJacobian(Element element, int elem_index,
                        double time,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] xpts,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] vars,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] dvars,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] ddvars,
                        int col,
                        double dh=1e-6,
                        int test_print_level=2,
                        double atol=1e-5, double rtol=1e-5):
    num_nodes = element.getNumNodes()
    num_vars = element.getNumVariables()

    # Make sure input arrays are large enough for element to avoid segfault
    assert len(xpts) >= 3 * num_nodes
    assert len(vars) >= num_vars
    assert len(dvars) >= num_vars
    assert len(ddvars) >= num_vars

    return TacsTestElementJacobian(element.ptr, elem_index, time,
                                   <TacsScalar*>xpts.data,
                                   <TacsScalar*>vars.data,
                                   <TacsScalar*>dvars.data,
                                   <TacsScalar*>ddvars.data,
                                   col, dh,
                                   test_print_level, atol, rtol)

def TestAdjResProduct(Element element, int elem_index,
                        double time,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] xpts,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] vars,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] dvars,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] ddvars,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] design_vars,
                        double dh=1e-6,
                        int test_print_level=2,
                        double atol=1e-5, double rtol=1e-5):
    num_nodes = element.getNumNodes()
    num_vars = element.getNumVariables()
    num_dvs = len(element.getDesignVars(elem_index))

    # Make sure input arrays are large enough for element to avoid segfault
    assert len(xpts) >= 3 * num_nodes
    assert len(vars) >= num_vars
    assert len(dvars) >= num_vars
    assert len(ddvars) >= num_vars
    assert len(design_vars) >= num_dvs

    return TacsTestAdjResProduct(element.ptr, elem_index, time,
                                   <TacsScalar*>xpts.data,
                                   <TacsScalar*>vars.data,
                                   <TacsScalar*>dvars.data,
                                   <TacsScalar*>ddvars.data,
                                   num_dvs, <TacsScalar*>design_vars.data, dh,
                                   test_print_level, atol, rtol)

def TestAdjResXptProduct(Element element, int elem_index,
                        double time,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] xpts,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] vars,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] dvars,
                        np.ndarray[TacsScalar, ndim=1, mode='c'] ddvars,
                        double dh=1e-6,
                        int test_print_level=2,
                        double atol=1e-5, double rtol=1e-5):
    num_nodes = element.getNumNodes()
    num_vars = element.getNumVariables()

    # Make sure input arrays are large enough for element to avoid segfault
    assert len(xpts) >= 3 * num_nodes
    assert len(vars) >= num_vars
    assert len(dvars) >= num_vars
    assert len(ddvars) >= num_vars

    return TacsTestAdjResXptProduct(element.ptr, elem_index, time,
                                   <TacsScalar*>xpts.data,
                                   <TacsScalar*>vars.data,
                                   <TacsScalar*>dvars.data,
                                   <TacsScalar*>ddvars.data, dh,
                                   test_print_level, atol, rtol)

def TestElementMatDVSens(Element element, ElementMatrixType mat_type,
                         int elem_index, double time,
                         np.ndarray[TacsScalar, ndim=1, mode='c'] xpts,
                         np.ndarray[TacsScalar, ndim=1, mode='c'] vars,
                         np.ndarray[TacsScalar, ndim=1, mode='c'] design_vars,
                         double dh=1e-6,
                         int test_print_level=2,
                         double atol=1e-5, double rtol=1e-5):
    num_nodes = element.getNumNodes()
    num_vars = element.getNumVariables()
    num_dvs = len(element.getDesignVars(elem_index))

    # Make sure input arrays are large enough for element to avoid segfault
    assert len(xpts) >= 3 * num_nodes
    assert len(vars) >= num_vars
    assert len(design_vars) >= num_dvs

    return TacsTestElementMatDVSens(element.ptr, mat_type, elem_index, time,
                                   <TacsScalar*>xpts.data,
                                   <TacsScalar*>vars.data,
                                   num_dvs, <TacsScalar*>design_vars.data, dh,
                                   test_print_level, atol, rtol)

def TestElementMatXptSens(Element element, ElementMatrixType mat_type,
                          int elem_index, double time,
                          np.ndarray[TacsScalar, ndim=1, mode='c'] xpts,
                          np.ndarray[TacsScalar, ndim=1, mode='c'] vars,
                          double dh=1e-6,
                          int test_print_level=2,
                          double atol=1e-5, double rtol=1e-5):
    num_nodes = element.getNumNodes()
    num_vars = element.getNumVariables()

    # Make sure input arrays are large enough for element to avoid segfault
    assert len(xpts) >= 3 * num_nodes
    assert len(vars) >= num_vars

    return TacsTestElementMatXptSens(element.ptr, mat_type, elem_index, time,
                                    <TacsScalar*>xpts.data,
                                    <TacsScalar*>vars.data, dh,
                                    test_print_level, atol, rtol)

def TestElementMatSVSens(Element element, ElementMatrixType mat_type,
                         int elem_index, double time,
                         np.ndarray[TacsScalar, ndim=1, mode='c'] xpts,
                         np.ndarray[TacsScalar, ndim=1, mode='c'] vars,
                         double dh=1e-6,
                         int test_print_level=2,
                         double atol=1e-5, double rtol=1e-5):
    num_nodes = element.getNumNodes()
    num_vars = element.getNumVariables()

    # Make sure input arrays are large enough for element to avoid segfault
    assert len(xpts) >= 3 * num_nodes
    assert len(vars) >= num_vars

    return TacsTestElementMatSVSens(element.ptr, mat_type, elem_index, time,
                                   <TacsScalar*>xpts.data,
                                   <TacsScalar*>vars.data, dh,
                                   test_print_level, atol, rtol)

# Function for setting random seeds for all element tests
def SeedRandomGenerator(int seed=0):
    TacsSeedRandomGenerator(seed)
    return

cdef class LinearTetrahedralBasis(ElementBasis):
    """
    Basis class for a linear tetrahedral 3D element.
    """
    def __cinit__(self):
        self.ptr = new TACSLinearTetrahedralBasis()
        self.ptr.incref()

cdef class QuadraticTetrahedralBasis(ElementBasis):
    """
    Basis class for a quadratic tetrahedral 3D element.
    """
    def __cinit__(self):
        self.ptr = new TACSQuadraticTetrahedralBasis()
        self.ptr.incref()

cdef class CubicTetrahedralBasis(ElementBasis):
    """
    Basis class for a cubic tetrahedral 3D element.
    """
    def __cinit__(self):
        self.ptr = new TACSCubicTetrahedralBasis()
        self.ptr.incref()


cdef class LinearHexaBasis(ElementBasis):
    """
    Basis class for a linear hexahedral 3D element.
    """
    def __cinit__(self):
        self.ptr = new TACSLinearHexaBasis()
        self.ptr.incref()

cdef class QuadraticHexaBasis(ElementBasis):
    """
    Basis class for a quadratic hexahedral 3D element.
    """
    def __cinit__(self):
        self.ptr = new TACSQuadraticHexaBasis()
        self.ptr.incref()

cdef class CubicHexaBasis(ElementBasis):
    """
    Basis class for a cubic hexahedral 3D element.
    """
    def __cinit__(self):
        self.ptr = new TACSCubicHexaBasis()
        self.ptr.incref()


cdef class LinearQuadBasis(ElementBasis):
    """
    Basis class for a linear quad 2D element.
    """
    def __cinit__(self):
        self.ptr = new TACSLinearQuadBasis()
        self.ptr.incref()

cdef class QuadraticQuadBasis(ElementBasis):
    """
    Basis class for a quadratic quad 2D element.
    """
    def __cinit__(self):
        self.ptr = new TACSQuadraticQuadBasis()
        self.ptr.incref()

cdef class CubicQuadBasis(ElementBasis):
    """
    Basis class for a cubic quad 2D element.
    """
    def __cinit__(self):
        self.ptr = new TACSCubicQuadBasis()
        self.ptr.incref()

cdef class QuarticQuadBasis(ElementBasis):
    """
    Basis class for a quartic quad 2D element.
    """
    def __cinit__(self):
        self.ptr = new TACSQuarticQuadBasis()
        self.ptr.incref()

cdef class QuinticQuadBasis(ElementBasis):
    """
    Basis class for a quintic quad 2D element.
    """
    def __cinit__(self):
        self.ptr = new TACSQuinticQuadBasis()
        self.ptr.incref()


cdef class LinearTriangleBasis(ElementBasis):
    """
    Basis class for a linear triangular 2D element
    """
    def __cinit__(self):
        self.ptr = new TACSLinearTriangleBasis()
        self.ptr.incref()

cdef class QuadraticTriangleBasis(ElementBasis):
    """
    Basis class for a quadratic triangular 2D element
    """
    def __cinit__(self):
        self.ptr = new TACSQuadraticTriangleBasis()
        self.ptr.incref()

cdef class CubicTriangleBasis(ElementBasis):
    """
    Basis class for a cubic triangular 2D element
    """
    def __cinit__(self):
        self.ptr = new TACSCubicTriangleBasis()
        self.ptr.incref()

cdef class HeatConduction2D(ElementModel):
    """
    Model class for 2D heat conduction element.

    .. note::
        varsPerNode: 1

    Args:
        con (PlaneStressConstitutive): Material constitutive properties.
    """
    def __cinit__(self, PlaneStressConstitutive con):
        self.ptr = new TACSHeatConduction2D(con.cptr)
        self.ptr.incref()

cdef class PCMHeatConduction2D(ElementModel):
    """
    Model class for 2D phase change material heat conduction element.

    .. note::
        varsPerNode: 1

    Args:
        con (PhaseChangeMaterialConstitutive): Material constitutive properties.
    """
    def __cinit__(self, PhaseChangeMaterialConstitutive con):
        self.ptr = new TACSPCMHeatConduction2D(con.cptr)
        self.ptr.incref()

cdef class LinearElasticity2D(ElementModel):
    """
    Model class for 2D linear elasticity element.

    .. note::
        varsPerNode: 2

    Args:
        con (PlaneStressConstitutive): Material constitutive properties.
    """
    def __cinit__(self, PlaneStressConstitutive con):
        self.ptr = new TACSLinearElasticity2D(con.cptr, TACS_LINEAR_STRAIN)
        self.ptr.incref()

cdef class LinearThermoelasticity2D(ElementModel):
    """
    Model class for 2D linear thermoelasticity element.

    .. note::
        varsPerNode: 3

    Args:
        con (PlaneStressConstitutive): Material constitutive properties.
        steady_flag (bool, optional): Steady state flag. Defaults to False.
    """
    def __cinit__(self, PlaneStressConstitutive con, int steady_flag=0):
        self.ptr = new TACSLinearThermoelasticity2D(con.cptr, TACS_LINEAR_STRAIN,
                                                    steady_flag)
        self.ptr.incref()


cdef class HeatConduction3D(ElementModel):
    """
    Model class for 3D heat conduction element.

    .. note::
        varsPerNode: 1

    Args:
        con (SolidConstitutive): Material constitutive properties.
    """
    def __cinit__(self, SolidConstitutive con):
        self.ptr = new TACSHeatConduction3D(con.cptr)
        self.ptr.incref()

cdef class LinearElasticity3D(ElementModel):
    """
    Model class for 3D linear elasticity element.

    .. note::
        varsPerNode: 3

    Args:
        con (SolidConstitutive): Material constitutive properties.
    """
    cdef TACSLinearElasticity3D *leptr
    def __cinit__(self, SolidConstitutive con):
        self.leptr = new TACSLinearElasticity3D(con.cptr, TACS_LINEAR_STRAIN)
        self.ptr = self.leptr
        self.ptr.incref()

    def getConstitutive(self):
        if self.leptr:
            scon = SolidConstitutive()
            scon.ptr = self.leptr.getConstitutive()
            scon.ptr.incref()
            return scon
        return None

cdef class LinearThermoelasticity3D(ElementModel):
    """
    Model class for 3D linear thermoelasticity element.

    .. note::
        varsPerNode: 4

    Args:
        con (SolidConstitutive): Material constitutive properties.
        steady_flag (bool, optional): Steady state flag. Defaults to False.
    """
    def __cinit__(self, SolidConstitutive con, int steady_flag=0):
        self.ptr = new TACSLinearThermoelasticity3D(con.cptr, TACS_LINEAR_STRAIN,
                                                    steady_flag)
        self.ptr.incref()

cdef class PlateModel(ElementModel):
    def __cinit__(self, ShellConstitutive con):
        self.ptr = new TACSPlateModel(con.cptr)
        self.ptr.incref()

cdef class ThermoelasticPlateModel(ElementModel):
    def __cinit__(self, ShellConstitutive con):
        self.ptr = new TACSThermoelasticPlateModel(con.cptr)
        self.ptr.incref()

cdef class Element2D(Element):
    """
    General element class appropriate for 2D analysis.

    .. note::
        varsPerNode: depends on ``ElementModel``

        outputElement: ``TACS.PLANE_STRESS_ELEMENT``

    Args:
        model (ElementModel): Physics model for element.
        basis (ElementBasis): Basis for element.
    """
    def __cinit__(self, ElementModel model, ElementBasis basis):
        self.ptr = new TACSElement2D(model.ptr, basis.ptr)
        self.ptr.incref()

cdef class Element3D(Element):
    """
    General element class appropriate for 3D analysis.

    .. note::
        varsPerNode: depends on ``ElementModel``

        outputElement: ``TACS.SOLID_ELEMENT``

    Args:
        model (ElementModel): Physics model for element.
        basis (ElementBasis): Basis for element.
    """
    def __cinit__(self, ElementModel model, ElementBasis basis):
        self.ptr = new TACSElement3D(model.ptr, basis.ptr)
        self.ptr.incref()

cdef class Traction2D(Element):
    def __cinit__(self, int varsPerNode, int faceIndex,
                  ElementBasis basis, list traction, coordComp=True):
        cdef TacsScalar trac[16]
        cdef int tractionCoordinateComponent = 0
        if coordComp:
            tractionCoordinateComponent = 1
        for i in range(min(16, len(traction))):
            trac[i] = traction[i]
        self.ptr = new TACSTraction2D(varsPerNode, faceIndex, basis.ptr,
                                      trac, tractionCoordinateComponent)
        self.ptr.incref()

cdef class Traction3D(Element):
    def __cinit__(self, int varsPerNode, int faceIndex,
                  ElementBasis basis, list traction, coordComp=True):
        cdef TacsScalar trac[24]
        cdef int tractionCoordinateComponent = 0
        if coordComp:
            tractionCoordinateComponent = 1
        for i in range(min(24, len(traction))):
            trac[i] = traction[i]
        self.ptr = new TACSTraction3D(varsPerNode, faceIndex, basis.ptr,
                                      trac, tractionCoordinateComponent)
        self.ptr.incref()

cdef class ConvectiveTraction2D(Element):
    def __cinit__(self, int varsPerNode, int faceIndex,
                  int fieldIndex, TacsScalar alpha,
                  TacsScalar refValue,
                  ElementBasis basis):
        self.ptr = new TACSConvectiveTraction2D(varsPerNode,
                                                faceIndex,
                                                fieldIndex,
                                                alpha,
                                                refValue,
                                                basis.ptr)
        self.ptr.incref()

cdef class ConvectiveTraction3D(Element):
    def __cinit__(self, int varsPerNode, int faceIndex,
                  int fieldIndex, TacsScalar alpha,
                  TacsScalar refValue,
                  ElementBasis basis):
        self.ptr = new TACSConvectiveTraction3D(varsPerNode,
                                                faceIndex,
                                                fieldIndex,
                                                alpha,
                                                refValue,
                                                basis.ptr)
        self.ptr.incref()

cdef class ShellTransform:
    cdef TACSShellTransform *ptr
    def __cinit__(self):
        self.ptr = NULL
        return

    def __dealloc__(self):
        if self.ptr != NULL:
            self.ptr.decref()

cdef class ShellNaturalTransform(ShellTransform):
    """
    Class for computing the transformation from the global coordinates
    to the shell/material local coordinates (used for defining stresses).
    This class uses shell "natural" coordinate system (i.e. local :math:`x` is aligned with element's first edge.)
    This is appropriate for isotropic shells, who's stiffness properties don't depend on orientation.

    .. seealso:: :ref:`theory/shell_element:Natural transform`
    """
    def __cinit__(self):
        self.ptr = new TACSShellNaturalTransform()
        self.ptr.incref()

cdef class ShellRefAxisTransform(ShellTransform):
    """
    Class for computing the transformation from the global coordinates
    to the shell local coordinates (used for defining stresses).
    This class uses a projection of a user-supplied reference axis to define the shell/material coordinate system
    (i.e. local :math:`x` is aligned with the reference axis.)
    This is appropriate for orthotropic and anisotropic shells, who's stiffness properties depend on orientation.

    .. warning:: The reference direction cannot be normal to the shell surface.

    .. seealso:: :ref:`theory/shell_element:Reference axis projection transform`

    Args:
        axis (array-like): Reference axis.
    """
    def __cinit__(self, axis):
        cdef TacsScalar a[3]
        a[0] = axis[0]
        a[1] = axis[1]
        a[2] = axis[2]
        self.ptr = new TACSShellRefAxisTransform(a)
        self.ptr.incref()

cdef class Quad4Shell(Element):
    """
    A 4-node quad shell element for general linear elastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSQuad4Shell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Quad4NonlinearShell(Element):
    """
    A 4-node quad shell element for general geometric nonlinear elastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSQuad4NonlinearShell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Quad9Shell(Element):
    """
    A 9-node quad shell element for general linear elastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSQuad9Shell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Quad9NonlinearShell(Element):
    """
    A 9-node quad shell element for general geometric nonlinear elastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSQuad9NonlinearShell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Quad16Shell(Element):
    """
    A 16-node quad shell element for general linear elastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSQuad16Shell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Quad16NonlinearShell(Element):
    """
    A 16-node quad shell element for general geometric nonlinear elastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSQuad16NonlinearShell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Tri3Shell(Element):
    """
    A 3-node triangular shell element for general linear elastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSTri3Shell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Tri3NonlinearShell(Element):
    """
    A 3-node triangular shell element for geometric elastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSTri3NonlinearShell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Quad4ThermalShell(Element):
    """
    A 4-node quad shell element for general linear thermoelastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 7

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSQuad4ThermalShell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Quad4NonlinearThermalShell(Element):
    """
    A 4-node quad shell element for general linear thermoelastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 7

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSQuad4NonlinearThermalShell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Quad4ShellQuaternion(Element):
    """
    A 4-node quad shell element for general linear thermoelastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 8

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSQuad4ShellQuaternion(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Quad4ShellModRot(Element):
    """
    A 4-node quad shell element for general linear thermoelastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: ??

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSQuad4ShellModRot(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Quad9ThermalShell(Element):
    """
    A 9-node quad shell element for general linear thermoelastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 7

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSQuad9ThermalShell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Quad9NonlinearThermalShell(Element):
    """
    A 9-node quad shell element for general geometric nonlinear thermoelastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 7

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSQuad9NonlinearThermalShell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Quad16ThermalShell(Element):
    """
    A 16-node quad shell element for general linear thermoelastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 7

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSQuad16ThermalShell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Quad16NonlinearThermalShell(Element):
    """
    A 16-node quad shell element for general geometric nonlinear thermoelastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 7

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSQuad16NonlinearThermalShell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Tri3ThermalShell(Element):
    """
    A 3-node triangular shell element for general linear thermoelastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 7

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSTri3ThermalShell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Tri3NonlinearThermalShell(Element):
    """
    A 3-node triangular shell element for general geometric nonlinear thermoelastic analysis.

    This element employs a mixed interpolation of tensorial (strain)
    components (MITC) method to avoid shear locking problems.

    .. seealso:: :ref:`theory/shell_element:Mixed Interpolation of Tensorial Components`

    .. note::
        varsPerNode: 7

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (ShellTransform or None): Shell transform object.
          ``None`` is equivalent to :class:`~ShellNaturalTransform`.
        con (ShellConstitutive): Shell constitutive object.
    """
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        if transform is None:
            transform = ShellNaturalTransform()
        self.ptr = new TACSTri3NonlinearThermalShell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class BeamTransform:
    cdef TACSBeamTransform *ptr
    def __cinit__(self):
        self.ptr = NULL
        return

    def __dealloc__(self):
        if self.ptr != NULL:
            self.ptr.decref()

cdef class BeamRefAxisTransform(BeamTransform):
    """
    Class for computing the transformation from the global coordinates
    to the beam local coordinates (used for defining stresses).
    This class uses a projection of a user-supplied reference axis to define the beam coordinate system
    (i.e. local :math:`y` is aligned with the reference axis.)

    .. warning:: The reference direction cannot be parallel to the beam axis.

    Args:
        axis (array-like): Reference axis.
    """
    def __cinit__(self, axis):
        cdef TacsScalar a[3]
        a[0] = axis[0]
        a[1] = axis[1]
        a[2] = axis[2]
        self.ptr = new TACSBeamRefAxisTransform(a)
        self.ptr.incref()

cdef class Beam2(Element):
    """
    A 2-node Timoshenko beam element for general linear elastic analysis.

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (BeamTransform): Beam transform object.
        con (BeamConstitutive): Beam constitutive object.
    """
    def __cinit__(self, BeamTransform transform, BeamConstitutive con):
        self.ptr = new TACSBeam2(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Beam3(Element):
    """
    A 3-node Timoshenko beam element for general linear elastic analysis.

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (BeamTransform): Beam transform object.
        con (BeamConstitutive): Beam constitutive object.
    """
    def __cinit__(self, BeamTransform transform, BeamConstitutive con):
        self.ptr = new TACSBeam3(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Beam2ModRot(Element):
    """
    A 2-node Timoshenko beam element for general nonlinear elastic analysis
    with moderate rotations.

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (BeamTransform): Beam transform object.
        con (BeamConstitutive): Beam constitutive object.
    """
    def __cinit__(self, BeamTransform transform, BeamConstitutive con):
        self.ptr = new TACSBeam2ModRot(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class Beam3ModRot(Element):
    """
    A 3-node Timoshenko beam element for general nonlinear elastic analysis
    with moderate rotations.

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.BEAM_OR_SHELL_ELEMENT``

    Args:
        transform (BeamTransform): Beam transform object.
        con (BeamConstitutive): Beam constitutive object.
    """
    def __cinit__(self, BeamTransform transform, BeamConstitutive con):
        self.ptr = new TACSBeam3ModRot(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class SpringTransform:
    cdef TACSSpringTransform *ptr
    def __cinit__(self):
        self.ptr = NULL
        return

    def __dealloc__(self):
        if self.ptr != NULL:
            self.ptr.decref()

cdef class SpringIdentityTransform(SpringTransform):
    """
    Class for computing the transformation from the global coordinates
    to the spring/material local coordinates (used for defining spring stifness directions).
    This class uses an identity transform (i.e. the spring stiffnesses are aligned with the global axes.)
    """
    def __cinit__(self):
        self.ptr = new TACSSpringIdentityTransform()
        self.ptr.incref()

cdef class SpringRefAxisTransform(SpringTransform):
    """
    Class for computing the transformation from the global coordinates
    to the spring/material local coordinates (used for defining spring stiffness directions).
    This class uses a projection of a user-supplied reference axis to define the spring stiffness coordinate system.
    The local `x` always aligns with the element. The local `z` is given by the cross product of `x` and `axis`.
    The local `y` is given by the cross product of local `z` and `x`.

    .. warning:: Not appropriate for elements with coincident nodes. Use :class:`~SpringIdentityTransform` or
      :class:`~SpringRefFrameTransform`.

    Args:
        axis (array-like): Reference axis.
    """
    def __cinit__(self, axis):
        cdef TacsScalar a[3]
        a[0] = axis[0]
        a[1] = axis[1]
        a[2] = axis[2]
        self.ptr = new TACSSpringRefAxisTransform(a)
        self.ptr.incref()

cdef class SpringRefFrameTransform(SpringTransform):
    """
    Class for computing the transformation from the global coordinates
    to the spring/material local coordinates (used for defining spring stifness directions).
    This class uses an user-defined arbritrary reference frame for the transform (i.e. the spring stiffnesses are aligned with the global axes.)
    The reference is defined using two reference vectors. `axis1` defines the local `x` axes.
    The local `z` is given by the cross product of `x` and `axis2`.
    The local `y` is given by the cross product of local `z` and `x`.
    """
    def __cinit__(self, axis1, axis2):
        cdef TacsScalar a1[3], a2[3]
        a1[0] = axis1[0]
        a1[1] = axis1[1]
        a1[2] = axis1[2]
        a2[0] = axis2[0]
        a2[1] = axis2[1]
        a2[2] = axis2[2]
        self.ptr = new TACSSpringRefFrameTransform(a1, a2)
        self.ptr.incref()

cdef class SpringElement(Element):
    """
    A 6 DOF spring element.

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.SPRING_ELEMENT``

    Args:
        transform (SpringTransform or None): Spring transform object.
          ``None`` is equivalent to :class:`~SpringIdentityTransform`.
        con (GeneralSpringConstitutive): Spring constitutive object.
    """
    def __cinit__(self, SpringTransform transform, GeneralSpringConstitutive con):
        if transform is None:
            transform = SpringIdentityTransform()
        self.ptr = new TACSSpringElement(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class GibbsVector:
    cdef TACSGibbsVector *ptr
    def __cinit__(self, x, y, z):
        self.ptr = new TACSGibbsVector(x, y, z)
        self.ptr.incref()

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

cdef class RefFrame:
    cdef TACSRefFrame *ptr
    def __cinit__(self, GibbsVector r0, GibbsVector r1, GibbsVector r2):
        self.ptr = new TACSRefFrame(r0.ptr, r1.ptr, r2.ptr)
        self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

cdef class RigidBodyViz:
    cdef TACSRigidBodyViz *ptr
    def __cinit__(self,
                  int npts=0, int nelems=0,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] xpts=None,
                  np.ndarray[int, ndim=1, mode='c'] conn=None,
                  GibbsVector vref=None,
                  TacsScalar Lx=1.0, TacsScalar Ly=1.0, TacsScalar Lz=1.0):
        cdef TACSGibbsVector *vptr = NULL
        if vref is not None:
            vptr = vref.ptr
        if xpts is not None and conn is not None:
            self.ptr = new TACSRigidBodyViz(npts, nelems,
                                            <TacsScalar*>xpts.data,
                                            <int*>conn.data, vref.ptr)
        else:
            self.ptr = new TACSRigidBodyViz(Lx, Ly, Lz)

        self.ptr.incref()

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

cdef class RBE2(Element):
    """
    A rigid body connected to an arbitrary number of grid points.
    The independent degrees-of-freedom are the six components of motion at
    a single grid point. The dependent degrees-of-freedom at the other grid
    points.

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.RIGID_ELEMENT``

    Args:
        num_nodes (int): Total number of nodes associated with the element.
        constrained_dofs (numpy.ndarray[int]): Flags to determine which
          dependent node dof's are attached to the eleemnt.
    """
    cdef TACSRBE2 *cptr
    def __cinit__(self, int num_nodes,
                  np.ndarray[int, ndim=1, mode='c'] constrained_dofs,
                  double C1=1e3, double C2=1e-3):
        num_dep = (num_nodes - 1) / 2

        assert len(constrained_dofs) == 6 or len(constrained_dofs) == 6 * num_dep

        # Pad the dependent array, if necessary
        if len(constrained_dofs) == 6:
            constrained_dofs = np.tile(constrained_dofs, num_dep)

        self.cptr = new TACSRBE2(num_nodes, <int*>constrained_dofs.data, C1, C2)
        # Increase the reference count to the underlying object
        self.ptr = self.cptr
        self.ptr.incref()
        return

cdef class RBE3(Element):
    """
    The RBE3 element is a powerful tool for distributing applied
    loads and mass in a model. Unlike the RBAR and RBE2 elements,
    the RBE3 doesn’t add additional stiffness to your structure.
    Forces and moments applied to reference points are distributed
    to a set of independent degrees of freedom based on the RBE3
    geometry and local weight factors.

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.RIGID_ELEMENT``

    Args:
        num_nodes (int): Total number of nodes associated with the element.
        dep_constrained_dofs (numpy.ndarray[int]): Flags to determine which
          dependent node dof's are attached to the eleemnt.
        weights (numpy.ndarray[float]): RBE weighting factor for each independent node.
        indep_constrained_dofs (numpy.ndarray[int]): Flags to determine which
          independent node dof's are attached to the eleemnt.
    """
    cdef TACSRBE3 *cptr
    def __cinit__(self, int num_nodes,
                  np.ndarray[int, ndim=1, mode='c'] dep_constrained_dofs,
                  np.ndarray[double, ndim=1, mode='c'] weights,
                  np.ndarray[int, ndim=1, mode='c'] indep_constrained_dofs,
                  double C1=1e3, double C2=1e-3):
        num_indep = num_nodes - 2

        assert len(dep_constrained_dofs) == 6
        assert len(weights) == 1 or len(weights) == num_indep
        assert len(indep_constrained_dofs) == 6 or len(indep_constrained_dofs) == 6 * num_indep

        # Pad the independent arrays, if necessary
        if len(weights) == 1:
            weights = np.tile(weights, num_indep)
        if len(indep_constrained_dofs) == 6:
            indep_constrained_dofs = np.tile(indep_constrained_dofs, num_indep)

        self.cptr = new TACSRBE3(num_nodes, <int*>dep_constrained_dofs.data,
                                 <double*>weights.data, <int*>indep_constrained_dofs.data, C1, C2)
        # Increase the reference count to the underlying object
        self.ptr = self.cptr
        self.ptr.incref()
        return

cdef class MassElement(Element):
    """
    A 6 DOF point mass element.

    .. note::
        varsPerNode: 6

        outputElement: ``TACS.MASS_ELEMENT``

    Args:
        con (GeneralMassConstitutive): Point mass constitutive object.
    """
    cdef TACSMassElement *cptr
    def __cinit__(self, GeneralMassConstitutive con):

        self.cptr = new TACSMassElement(con.cptr)
        # Increase the reference count to the underlying object
        self.ptr = self.cptr
        self.ptr.incref()
        return

# cdef class RigidBody(Element):
#     cdef TACSRigidBody *cptr
#     def __cinit__(self, RefFrame frame, TacsScalar mass,
#                   np.ndarray[TacsScalar, ndim=1, mode='c'] cRef,
#                   np.ndarray[TacsScalar, ndim=1, mode='c'] JRef,
#                   GibbsVector r0,
#                   GibbsVector v0, GibbsVector omega0, GibbsVector g,
#                   int mdv=-1,
#                   np.ndarray[int, ndim=1, mode='c'] cdvs=None,
#                   np.ndarray[int, ndim=1, mode='c'] Jdvs=None):
#         cdef int *_cdvs = NULL
#         cdef int *_Jdvs = NULL

#         # Assign the the variable numbers if they are supplied by the
#         # user
#         if cdvs is not None:
#             _cdvs = <int*>cdvs.data
#         if Jdvs is not None:
#             _Jdvs = <int*>Jdvs.data

#         # Allocate the rigid body object and set the design variables
#         self.cptr = new TACSRigidBody(frame.ptr, mass,
#                                       <TacsScalar*>cRef.data,
#                                       <TacsScalar*>JRef.data, r0.ptr,
#                                       v0.ptr, omega0.ptr, g.ptr)
#         self.cptr.setDesignVarNums(mdv, _cdvs, _Jdvs)

#         # Increase the reference count to the underlying object
#         self.ptr = self.cptr
#         self.ptr.incref()

#     def setVisualization(self, RigidBodyViz viz):
#         self.cptr.setVisualization(viz.ptr)

# cdef class FixedConstraint(Element):
#     def __cinit__(self,
#                   GibbsVector point,
#                   RigidBody bodyA):
#         self.ptr = new TACSFixedConstraint(bodyA.cptr,
#                                            point.ptr)
#         self.ptr.incref()

# cdef class SphericalConstraint(Element):
#     def __cinit__(self,
#                   GibbsVector point,
#                   RigidBody bodyA, RigidBody bodyB=None):
#         if bodyB is None:
#             self.ptr = new TACSSphericalConstraint(bodyA.cptr,
#                                                    point.ptr)
#         else:
#             self.ptr = new TACSSphericalConstraint(bodyA.cptr, bodyB.cptr,
#                                                    point.ptr)
#         self.ptr.incref()

# cdef class RevoluteConstraint(Element):
#     def __cinit__(self, GibbsVector point, GibbsVector eA,
#                   int fixed_ref_point=0, int inertial_rev_axis=0,
#                   RigidBody bodyA=None, RigidBody bodyB=None):
#         if bodyA is not None and bodyB is not None:
#             self.ptr = new TACSRevoluteConstraint(bodyA.rbptr, bodyB.rbptr,
#                                                   point.ptr, eA.ptr)
#             self.ptr.incref()
#         elif bodyA is not None and bodyB is None:
#             self.ptr = new TACSRevoluteConstraint(bodyA.rbptr,
#                                                   point.ptr, eA.ptr)
#             self.ptr.incref()
#         else:
#             self.ptr = new TACSRevoluteConstraint(fixed_ref_point,
#                                                   point.ptr, eA.ptr,
#                                                   inertial_rev_axis)
#             self.ptr.incref()
#         return

# cdef class CylindricalConstraint(Element):
#     def __cinit__(self, GibbsVector point, GibbsVector eA,
#                   RigidBody bodyA, RigidBody bodyB=None):
#         if bodyB is None:
#             self.ptr = new TACSCylindricalConstraint(bodyA.cptr,
#                                                      point.ptr, eA.ptr)
#         else:
#             self.ptr = new TACSCylindricalConstraint(bodyA.cptr, bodyB.cptr,
#                                                      point.ptr, eA.ptr)
#         self.ptr.incref()

# cdef class RigidLink(Element):
#     def __cinit__(self, RigidBody bodyA):
#         self.ptr = new TACSRigidLink(bodyA.cptr)
#         self.ptr.incref()
#         return

# cdef class RevoluteDriver(Element):
#     def __cinit__(self, GibbsVector rev, TacsScalar omega):
#         self.ptr = new TACSRevoluteDriver(rev.ptr, omega)
#         self.ptr.incref()
#         return

# cdef class MotionDriver(Element):
#     def __cinit__(self, GibbsVector dir, TacsScalar omega,
#                   arrest_rot=False):
#         if arrest_rot is False:
#             self.ptr = new TACSMotionDriver(dir.ptr, omega, 0)
#         else:
#             self.ptr = new TACSMotionDriver(dir.ptr, omega, 1)
#         self.ptr.incref()

# cdef class AverageConstraint(Element):
#     def __cinit__(self, RigidBody body, GibbsVector point,
#                   RefFrame frame, int moment_flag=0):
#         self.ptr = new TACSAverageConstraint(body.cptr, point.ptr,
#                                              frame.ptr, moment_flag)
#         self.ptr.incref()

# cdef class MITCBeam(Element):
#     def __cinit__(self, BeamConstitutive stiff,
#                   GibbsVector gravity=None,
#                   GibbsVector vInit=None,
#                   GibbsVector omegaInit=None):
#         cdef TACSBeamConstitutive *con = _dynamicBeamConstitutive(stiff.ptr)
#         if omegaInit is not None:
#             self.ptr = new MITC3(con, gravity.ptr,
#                                  vInit.ptr, omegaInit.ptr)
#         elif vInit is not None:
#             self.ptr = new MITC3(con, gravity.ptr, vInit.ptr, NULL)
#         elif gravity is not None:
#             self.ptr = new MITC3(con, gravity.ptr, NULL, NULL)
#         else:
#             self.ptr = new MITC3(con, NULL, NULL, NULL)
#         self.ptr.incref()

# cdef inplace_array_1d(int nptype, int dim1, void *data_ptr):
#     '''Return a numpy version of the array'''
#     # Set the shape of the array
#     cdef int size = 1
#     cdef np.npy_intp shape[1]
#     cdef np.ndarray ndarray

#     # Set the first entry of the shape array
#     shape[0] = <np.npy_intp>dim1

#     # Create the array itself - Note that this function will not
#     # delete the data once the ndarray goes out of scope
#     ndarray = np.PyArray_SimpleNewFromData(size, shape,
#                                            nptype, data_ptr)

#     return ndarray

# cdef int getMultiplierIndex(void *self_ptr):
#     return (<object>self_ptr).getMultiplierIndex()

# cdef void getInitConditions(void *self_ptr, int elem_index, int num_nodes,
#                             const TacsScalar *Xpts, int num_vars,
#                             TacsScalar *vars,
#                             TacsScalar *dvars,
#                             TacsScalar *ddvars):
#     _Xpts = inplace_array_1d(TACS_NPY_SCALAR, 3*num_nodes, <void*>Xpts)
#     _vars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>vars)
#     _dvars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>dvars)
#     _ddvars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>ddvars)
#     (<object>self_ptr).getInitConditions(elem_index, _Xpts, _vars, _dvars, _ddvars)
#     return

# cdef void addResidual(void *self_ptr, int elem_index, double time,
#                       int num_nodes, const TacsScalar *Xpts,
#                       int num_vars, const TacsScalar *vars,
#                       const TacsScalar *dvars, const TacsScalar *ddvars,
#                       TacsScalar *res):
#     _Xpts = inplace_array_1d(TACS_NPY_SCALAR, 3*num_nodes, <void*>Xpts)
#     _vars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>vars)
#     _dvars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>dvars)
#     _ddvars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>ddvars)
#     _res = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>res)
#     (<object>self_ptr).addResidual(elem_index, time, _Xpts,
#                                    _vars, _dvars, _ddvars, _res)
#     return

# cdef void addJacobian(void *self_ptr, int elem_index, double time,
#                       TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
#                       int num_nodes, const TacsScalar *Xpts,
#                       int num_vars, const TacsScalar *vars,
#                       const TacsScalar *dvars, const TacsScalar *ddvars,
#                       TacsScalar *res, TacsScalar *mat):
#     _Xpts = inplace_array_1d(TACS_NPY_SCALAR, 3*num_nodes, <void*>Xpts)
#     _vars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>vars)
#     _dvars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>dvars)
#     _ddvars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>ddvars)
#     _res = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>res)
#     _mat = inplace_array_1d(TACS_NPY_SCALAR, num_vars*num_vars, <void*>mat)
#     (<object>self_ptr).addJacobian(elem_index, time, alpha, beta, gamma, _Xpts,
#                                    _vars, _dvars, _ddvars, _res, _mat)
#     return

# cdef class pyElement(Element):
#     def __cinit__(self, int vars_per_node, int num_nodes, *args, **kwargs):
#         cdef TACSElementWrapper *pointer
#         pointer = new TACSElementWrapper(<PyObject*>self, vars_per_node,
#                                          num_nodes)
#         pointer.incref()

#         # Set the function pointers
#         pointer.getmultiplierindex = getMultiplierIndex
#         pointer.getinitconditions = getInitConditions
#         pointer.addresidual = addResidual
#         pointer.addjacobian = addJacobian

#         self.ptr = pointer
