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

# distutils: language=c++

from __future__ import print_function, division
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

def TestElementBasis(ElementBasis basis, double dh=1e-6,
                     int test_print_level=2, double atol=1e-30,
                     double rtol=1e-5):
    return TacsTestElementBasis(basis.ptr, dh, test_print_level, atol, rtol)

cdef class LinearTetrahedralBasis(ElementBasis):
    def __cinit__(self):
        self.ptr = new TACSLinearTetrahedralBasis()
        self.ptr.incref()

cdef class QuadraticTetrahedralBasis(ElementBasis):
    def __cinit__(self):
        self.ptr = new TACSQuadraticTetrahedralBasis()
        self.ptr.incref()

cdef class CubicTetrahedralBasis(ElementBasis):
    def __cinit__(self):
        self.ptr = new TACSCubicTetrahedralBasis()
        self.ptr.incref()


cdef class LinearHexaBasis(ElementBasis):
    def __cinit__(self):
        self.ptr = new TACSLinearHexaBasis()
        self.ptr.incref()

cdef class QuadraticHexaBasis(ElementBasis):
    def __cinit__(self):
        self.ptr = new TACSQuadraticHexaBasis()
        self.ptr.incref()

cdef class CubicHexaBasis(ElementBasis):
    def __cinit__(self):
        self.ptr = new TACSCubicHexaBasis()
        self.ptr.incref()


cdef class LinearQuadBasis(ElementBasis):
    def __cinit__(self):
        self.ptr = new TACSLinearQuadBasis()
        self.ptr.incref()

cdef class QuadraticQuadBasis(ElementBasis):
    def __cinit__(self):
        self.ptr = new TACSQuadraticQuadBasis()
        self.ptr.incref()

cdef class CubicQuadBasis(ElementBasis):
    def __cinit__(self):
        self.ptr = new TACSCubicQuadBasis()
        self.ptr.incref()

cdef class QuarticQuadBasis(ElementBasis):
    def __cinit__(self):
        self.ptr = new TACSQuarticQuadBasis()
        self.ptr.incref()

cdef class QuinticQuadBasis(ElementBasis):
    def __cinit__(self):
        self.ptr = new TACSQuinticQuadBasis()
        self.ptr.incref()


cdef class LinearTriangleBasis(ElementBasis):
    def __cinit__(self):
        self.ptr = new TACSLinearTriangleBasis()
        self.ptr.incref()

cdef class QuadraticTriangleBasis(ElementBasis):
    def __cinit__(self):
        self.ptr = new TACSQuadraticTriangleBasis()
        self.ptr.incref()

cdef class CubicTriangleBasis(ElementBasis):
    def __cinit__(self):
        self.ptr = new TACSCubicTriangleBasis()
        self.ptr.incref()


cdef class HeatConduction2D(ElementModel):
    def __cinit__(self, PlaneStressConstitutive con):
        self.ptr = new TACSHeatConduction2D(con.cptr)
        self.ptr.incref()

cdef class LinearElasticity2D(ElementModel):
    def __cinit__(self, PlaneStressConstitutive con):
        self.ptr = new TACSLinearElasticity2D(con.cptr, TACS_LINEAR_STRAIN)
        self.ptr.incref()

cdef class LinearThermoelasticity2D(ElementModel):
    def __cinit__(self, PlaneStressConstitutive con):
        self.ptr = new TACSLinearThermoelasticity2D(con.cptr, TACS_LINEAR_STRAIN)
        self.ptr.incref()


cdef class HeatConduction3D(ElementModel):
    def __cinit__(self, SolidConstitutive con):
        self.ptr = new TACSHeatConduction3D(con.cptr)
        self.ptr.incref()

cdef class LinearElasticity3D(ElementModel):
    def __cinit__(self, SolidConstitutive con):
        self.ptr = new TACSLinearElasticity3D(con.cptr, TACS_LINEAR_STRAIN)
        self.ptr.incref()

cdef class LinearThermoelasticity3D(ElementModel):
    def __cinit__(self, SolidConstitutive con):
        self.ptr = new TACSLinearThermoelasticity3D(con.cptr, TACS_LINEAR_STRAIN)
        self.ptr.incref()


cdef class Element2D(Element):
    def __cinit__(self, ElementModel model, ElementBasis basis):
        self.ptr = new TACSElement2D(model.ptr, basis.ptr)
        self.ptr.incref()

cdef class Element3D(Element):
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


# cdef class GibbsVector:
#     cdef TACSGibbsVector *ptr
#     def __cinit__(self, x, y, z):
#         self.ptr = new TACSGibbsVector(x, y, z)
#         self.ptr.incref()

#     def __dealloc__(self):
#         self.ptr.decref()
#         return

# cdef class RefFrame:
#     cdef TACSRefFrame *ptr
#     def __cinit__(self, GibbsVector r0, GibbsVector r1, GibbsVector r2):
#         self.ptr = new TACSRefFrame(r0.ptr, r1.ptr, r2.ptr)
#         self.ptr.incref()
#         return

#     def __dealloc__(self):
#         self.ptr.decref()
#         return

# cdef class RigidBody(Element):
#     cdef TACSRigidBody *rbptr
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
#         self.rbptr = new TACSRigidBody(frame.ptr, mass,
#                                        <TacsScalar*>cRef.data,
#                                        <TacsScalar*>JRef.data, r0.ptr,
#                                        v0.ptr, omega0.ptr, g.ptr)
#         self.rbptr.setDesignVarNums(mdv, _cdvs, _Jdvs)

#         # Increase the reference count to the underlying object
#         self.ptr = self.rbptr
#         self.ptr.incref()
#         return

#     def setVisualization(self, RigidBodyViz viz):
#         self.rbptr.setVisualization(viz.ptr)

#     def __dealloc__(self):
#         self.ptr.decref()
#         return

#     def getNumNodes(self):
#         return self.ptr.numNodes()

#     def setComponentNum(self, int comp_num):
#         self.ptr.setComponentNum(comp_num)
#         return

# cdef class FixedConstraint(Element):
#     def __cinit__(self,
#                   GibbsVector point,
#                   RigidBody bodyA):
#         self.ptr = new TACSFixedConstraint(bodyA.rbptr,
#                                            point.ptr)
#         self.ptr.incref()
#         return

#     def __dealloc__(self):
#         self.ptr.decref()
#         return

#     def getNumNodes(self):
#         return self.ptr.numNodes()

#     def setComponentNum(self, int comp_num):
#         self.ptr.setComponentNum(comp_num)
#         return

# cdef class SphericalConstraint(Element):
#     def __cinit__(self,
#                   GibbsVector point,
#                   RigidBody bodyA, RigidBody bodyB=None):
#         if bodyB is None:
#             self.ptr = new TACSSphericalConstraint(bodyA.rbptr,
#                                                    point.ptr)
#         else:
#             self.ptr = new TACSSphericalConstraint(bodyA.rbptr, bodyB.rbptr,
#                                                    point.ptr)
#         self.ptr.incref()
#         return
#     def __dealloc__(self):
#         self.ptr.decref()
#         return
#     def numNodes(self):
#         return self.ptr.numNodes()
#     def setComponentNum(self, int comp_num):
#         self.ptr.setComponentNum(comp_num)
#         return

# cdef class RevoluteConstraint(Element):
#     def __cinit__(self, GibbsVector point, GibbsVector eA,
#                   RigidBody bodyA, RigidBody bodyB=None):
#         if bodyB is None:
#             self.ptr = new TACSRevoluteConstraint(bodyA.rbptr,
#                                                   point.ptr, eA.ptr)
#         else:
#             self.ptr = new TACSRevoluteConstraint(bodyA.rbptr, bodyB.rbptr,
#                                                   point.ptr, eA.ptr)
#         self.ptr.incref()
#         return
#     def __dealloc__(self):
#         self.ptr.decref()
#         return
#     def numNodes(self):
#         return self.ptr.numNodes()
#     def setComponentNum(self, int comp_num):
#         self.ptr.setComponentNum(comp_num)
#         return

# cdef class CylindricalConstraint(Element):
#     def __cinit__(self, GibbsVector point, GibbsVector eA,
#                   RigidBody bodyA, RigidBody bodyB=None):
#         if bodyB is None:
#             self.ptr = new TACSCylindricalConstraint(bodyA.rbptr,
#                                                      point.ptr, eA.ptr)
#         else:
#             self.ptr = new TACSCylindricalConstraint(bodyA.rbptr, bodyB.rbptr,
#                                                      point.ptr, eA.ptr)
#         self.ptr.incref()
#         return
#     def __dealloc__(self):
#         self.ptr.decref()
#         return
#     def numNodes(self):
#         return self.ptr.numNodes()
#     def setComponentNum(self, int comp_num):
#         self.ptr.setComponentNum(comp_num)
#         return

# cdef class RigidLink(Element):
#     def __cinit__(self, RigidBody bodyA):
#         self.ptr = new TACSRigidLink(bodyA.rbptr)
#         self.ptr.incref()
#         return
#     def __dealloc__(self):
#         self.ptr.decref()
#         return
#     def numNodes(self):
#         return self.ptr.numNodes()
#     def setComponentNum(self, int comp_num):
#         self.ptr.setComponentNum(comp_num)
#         return

# cdef class RevoluteDriver(Element):
#     def __cinit__(self, GibbsVector rev, TacsScalar omega):
#         self.ptr = new TACSRevoluteDriver(rev.ptr, omega)
#         self.ptr.incref()
#         return
#     def __dealloc__(self):
#         self.ptr.decref()
#         return
#     def numNodes(self):
#         return self.ptr.numNodes()
#     def setComponentNum(self, int comp_num):
#         self.ptr.setComponentNum(comp_num)
#         return

# cdef class MotionDriver(Element):
#     def __cinit__(self, GibbsVector dir, TacsScalar omega,
#                   arrest_rot=False):
#         if arrest_rot is False:
#             self.ptr = new TACSMotionDriver(dir.ptr, omega, 0)
#         else:
#             self.ptr = new TACSMotionDriver(dir.ptr, omega, 1)
#         self.ptr.incref()
#         return
#     def __dealloc__(self):
#         self.ptr.decref()
#         return
#     def numNodes(self):
#         return self.ptr.numNodes()
#     def setComponentNum(self, int comp_num):
#         self.ptr.setComponentNum(comp_num)
#         return

# cdef class AverageConstraint(Element):
#     def __cinit__(self, RigidBody body, GibbsVector point,
#                   RefFrame frame, int moment_flag=0):
#         self.ptr = new TACSAverageConstraint(body.rbptr, point.ptr,
#                                              frame.ptr, moment_flag)
#         self.ptr.incref()
#         return
#     def __dealloc__(self):
#         self.ptr.decref()
#         return

cdef inplace_array_1d(int nptype, int dim1, void *data_ptr):
    '''Return a numpy version of the array'''
    # Set the shape of the array
    cdef int size = 1
    cdef np.npy_intp shape[1]
    cdef np.ndarray ndarray

    # Set the first entry of the shape array
    shape[0] = <np.npy_intp>dim1

    # Create the array itself - Note that this function will not
    # delete the data once the ndarray goes out of scope
    ndarray = np.PyArray_SimpleNewFromData(size, shape,
                                           nptype, data_ptr)

    return ndarray

cdef void getInitConditions(void *self_ptr, int elem_index, int num_nodes,
                            const TacsScalar *Xpts, int num_vars,
                            TacsScalar *vars,
                            TacsScalar *dvars,
                            TacsScalar *ddvars):
    _Xpts = inplace_array_1d(TACS_NPY_SCALAR, 3*num_nodes, <void*>Xpts)
    _vars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>vars)
    _dvars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>dvars)
    _ddvars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>ddvars)
    (<object>self_ptr).getInitConditions(elem_index, _Xpts, _vars, _dvars, _ddvars)
    return

cdef void addResidual(void *self_ptr, int elem_index, double time,
                      int num_nodes, const TacsScalar *Xpts,
                      int num_vars, const TacsScalar *vars,
                      const TacsScalar *dvars, const TacsScalar *ddvars,
                      TacsScalar *res):
    _Xpts = inplace_array_1d(TACS_NPY_SCALAR, 3*num_nodes, <void*>Xpts)
    _vars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>vars)
    _dvars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>dvars)
    _ddvars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>ddvars)
    _res = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>res)
    (<object>self_ptr).addResidual(elem_index, time, _Xpts,
                                   _vars, _dvars, _ddvars, _res)
    return

cdef void addJacobian(void *self_ptr, int elem_index, double time,
                      TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                      int num_nodes, const TacsScalar *Xpts,
                      int num_vars, const TacsScalar *vars,
                      const TacsScalar *dvars, const TacsScalar *ddvars,
                      TacsScalar *res, TacsScalar *mat):
    _Xpts = inplace_array_1d(TACS_NPY_SCALAR, 3*num_nodes, <void*>Xpts)
    _vars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>vars)
    _dvars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>dvars)
    _ddvars = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>ddvars)
    _res = inplace_array_1d(TACS_NPY_SCALAR, num_vars, <void*>res)
    _mat = inplace_array_1d(TACS_NPY_SCALAR, num_vars*num_vars, <void*>mat)
    (<object>self_ptr).addJacobian(elem_index, time, alpha, beta, gamma, _Xpts,
                                   _vars, _dvars, _ddvars, _res, _mat)
    return

cdef class pyElement(Element):
    def __cinit__(self, int vars_per_node, int num_nodes, *args, **kwargs):
        cdef TACSElementWrapper *pointer
        pointer = new TACSElementWrapper(<PyObject*>self, vars_per_node,
                                         num_nodes)
        pointer.incref()

        # Set the function pointers
        pointer.getinitconditions = getInitConditions
        pointer.addresidual = addResidual
        pointer.addjacobian = addJacobian

        self.ptr = pointer
