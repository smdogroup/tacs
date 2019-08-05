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

# For the use of MPI
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy
import numpy as np
cimport numpy as np

# Ensure that numpy is initialized
np.import_array()

# Import the definition required for const strings from libc.string cimport const_char

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

# Import the definitions
from TACS cimport *
from constitutive cimport *
from elements cimport *

# Include the definitions
include "TacsDefs.pxi"

# Include the mpi4py header
cdef extern from "mpi-compat.h":
    pass

def setElementFDStepSize(double dh):
    TACSSetElementFDStepSize(dh)

# Wrap the element behavior types
PY_LINEAR = LINEAR
PY_NONLINEAR = NONLINEAR
PY_LARGE_ROTATION = LARGE_ROTATION

# A generic wrapper class for the TACSElement object
cdef class Element:
    '''Base element class'''
    def __cinit__(self):
        self.ptr = NULL
        return

    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return

    def numNodes(self):
        return self.ptr.numNodes()

cdef class GibbsVector:
    cdef TACSGibbsVector *ptr
    def __cinit__(self, x, y, z):
        self.ptr = new TACSGibbsVector(x, y, z)
        self.ptr.incref()

    def __dealloc__(self):
        self.ptr.decref()
        return

cdef class RefFrame:
    cdef TACSRefFrame *ptr
    def __cinit__(self, GibbsVector r0, GibbsVector r1, GibbsVector r2):
        self.ptr = new TACSRefFrame(r0.ptr, r1.ptr, r2.ptr)
        self.ptr.incref()
        return

    def __dealloc__(self):
        self.ptr.decref()
        return

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
        return

    def __dealloc__(self):
        self.ptr.decref()
        return

cdef class RigidBody(Element):
    cdef TACSRigidBody *rbptr
    def __cinit__(self, RefFrame frame, TacsScalar mass,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] cRef,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] JRef,
                  GibbsVector r0,
                  GibbsVector v0, GibbsVector omega0, GibbsVector g,
                  int mdv=-1,
                  np.ndarray[int, ndim=1, mode='c'] cdvs=None,
                  np.ndarray[int, ndim=1, mode='c'] Jdvs=None):
        cdef int *_cdvs = NULL
        cdef int *_Jdvs = NULL

        # Assign the the variable numbers if they are supplied by the
        # user
        if cdvs is not None:
            _cdvs = <int*>cdvs.data
        if Jdvs is not None:
            _Jdvs = <int*>Jdvs.data

        # Allocate the rigid body object and set the design variables
        self.rbptr = new TACSRigidBody(frame.ptr, mass,
                                       <TacsScalar*>cRef.data,
                                       <TacsScalar*>JRef.data, r0.ptr,
                                       v0.ptr, omega0.ptr, g.ptr)
        self.rbptr.setDesignVarNums(mdv, _cdvs, _Jdvs)

        # Increase the reference count to the underlying object
        self.ptr = self.rbptr
        self.ptr.incref()
        return

    def setVisualization(self, RigidBodyViz viz):
        self.rbptr.setVisualization(viz.ptr)
    def __dealloc__(self):
        self.ptr.decref()
        return

    def numNodes(self):
        return self.ptr.numNodes()

    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return

cdef class FixedConstraint(Element):
    def __cinit__(self,
                  GibbsVector point,
                  RigidBody bodyA):
        self.ptr = new TACSFixedConstraint(bodyA.rbptr,
                                           point.ptr)
        self.ptr.incref()
        return
    def __dealloc__(self):
        self.ptr.decref()
        return
    def numNodes(self):
        return self.ptr.numNodes()
    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return

cdef class SphericalConstraint(Element):
    def __cinit__(self,
                  GibbsVector point,
                  RigidBody bodyA, RigidBody bodyB=None):
        if bodyB is None:
            self.ptr = new TACSSphericalConstraint(bodyA.rbptr,
                                                   point.ptr)
        else:
            self.ptr = new TACSSphericalConstraint(bodyA.rbptr, bodyB.rbptr,
                                                   point.ptr)
        self.ptr.incref()
        return
    def __dealloc__(self):
        self.ptr.decref()
        return
    def numNodes(self):
        return self.ptr.numNodes()
    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return

cdef class RevoluteConstraint(Element):
    def __cinit__(self, GibbsVector point, GibbsVector eA,
                  RigidBody bodyA, RigidBody bodyB=None):
        if bodyB is None:
            self.ptr = new TACSRevoluteConstraint(bodyA.rbptr,
                                                  point.ptr, eA.ptr)
        else:
            self.ptr = new TACSRevoluteConstraint(bodyA.rbptr, bodyB.rbptr,
                                                  point.ptr, eA.ptr)
        self.ptr.incref()
        return
    def __dealloc__(self):
        self.ptr.decref()
        return
    def numNodes(self):
        return self.ptr.numNodes()
    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return

cdef class CylindricalConstraint(Element):
    def __cinit__(self, GibbsVector point, GibbsVector eA,
                  RigidBody bodyA, RigidBody bodyB=None):
        if bodyB is None:
            self.ptr = new TACSCylindricalConstraint(bodyA.rbptr,
                                                     point.ptr, eA.ptr)
        else:
            self.ptr = new TACSCylindricalConstraint(bodyA.rbptr, bodyB.rbptr,
                                                     point.ptr, eA.ptr)
        self.ptr.incref()
        return
    def __dealloc__(self):
        self.ptr.decref()
        return
    def numNodes(self):
        return self.ptr.numNodes()
    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return

cdef class PrismaticConstraint(Element):
    def __cinit__(self, GibbsVector point, GibbsVector eA,
                  RigidBody bodyA, RigidBody bodyB=None):
        if bodyB is None:
            self.ptr = new TACSPrismaticConstraint(bodyA.rbptr,
                                                   point.ptr, eA.ptr)
        else:
            self.ptr = new TACSPrismaticConstraint(bodyA.rbptr, bodyB.rbptr,
                                                   point.ptr, eA.ptr)
        self.ptr.incref()
        return
    def __dealloc__(self):
        self.ptr.decref()
        return
    def numNodes(self):
        return self.ptr.numNodes()
    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return

cdef class SlidingPivotConstraint(Element):
    def __cinit__(self, GibbsVector point, GibbsVector eA,
                  RigidBody bodyA, RigidBody bodyB=None):
        if bodyB is None:
            self.ptr = new TACSSlidingPivotConstraint(bodyA.rbptr,
                                                      point.ptr, eA.ptr)
        else:
            self.ptr = new TACSSlidingPivotConstraint(bodyA.rbptr, bodyB.rbptr,
                                                      point.ptr, eA.ptr)
        self.ptr.incref()
        return
    def __dealloc__(self):
        self.ptr.decref()
        return
    def numNodes(self):
        return self.ptr.numNodes()
    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return

cdef class RigidLink(Element):
    def __cinit__(self, RigidBody bodyA):
        self.ptr = new TACSRigidLink(bodyA.rbptr)
        self.ptr.incref()
        return
    def __dealloc__(self):
        self.ptr.decref()
        return
    def numNodes(self):
        return self.ptr.numNodes()
    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return

cdef class RevoluteDriver(Element):
    def __cinit__(self, GibbsVector rev, TacsScalar omega):
        self.ptr = new TACSRevoluteDriver(rev.ptr, omega)
        self.ptr.incref()
        return
    def __dealloc__(self):
        self.ptr.decref()
        return
    def numNodes(self):
        return self.ptr.numNodes()
    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return

cdef class MotionDriver(Element):
    def __cinit__(self, GibbsVector dir, TacsScalar omega,
                  arrest_rot=False):
        if arrest_rot is False:
            self.ptr = new TACSMotionDriver(dir.ptr, omega, 0)
        else:
            self.ptr = new TACSMotionDriver(dir.ptr, omega, 1)
        self.ptr.incref()
        return
    def __dealloc__(self):
        self.ptr.decref()
        return
    def numNodes(self):
        return self.ptr.numNodes()
    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return

cdef class LinearizedMotionDriver(Element):
    def __cinit__(self, GibbsVector dir, TacsScalar omega,
                  arrest_rot=False):
        if arrest_rot is False:
            self.ptr = new TACSLinearizedMotionDriver(dir.ptr, omega, 0)
        else:
            self.ptr = new TACSLinearizedMotionDriver(dir.ptr, omega, 1)
        self.ptr.incref()
        return
    def __dealloc__(self):
        self.ptr.decref()
        return
    def numNodes(self):
        return self.ptr.numNodes()
    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return

cdef class AverageConstraint(Element):
    def __cinit__(self, RigidBody body, GibbsVector point,
                  RefFrame frame, int moment_flag=0):
        self.ptr = new TACSAverageConstraint(body.rbptr, point.ptr,
                                             frame.ptr, moment_flag)
        self.ptr.incref()
        return
    def __dealloc__(self):
        self.ptr.decref()
        return

cdef class PlaneQuad(Element):
    def __cinit__(self, int order, PlaneStress stiff,
                  ElementBehaviorType elem_type=LINEAR,
                  int component_num=0):
        '''
        Wrap the PlaneStressQuad element class for order 2,3,4
        '''
        cdef PlaneStressStiffness *con = _dynamicPlaneStress(stiff.ptr)
        if order == 2:
            self.ptr = new PlaneStressQuad2(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new PlaneStressQuad3(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new PlaneStressQuad4(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 5:
            self.ptr = new PlaneStressQuad5(con, elem_type, component_num)
            self.ptr.incref()
        return

    def __dealloc__(self):
        self.ptr.decref()
        return

cdef class PSQuadTraction(Element):
    def __cinit__(self, int surf,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] tx,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] ty):
        if len(tx) != len(ty):
            errmsg = 'Traction lengths must be equal'
            raise ValueError(errmsg)
        if len(tx) < 2 or len(tx) > 5:
            errmsg = 'Traction lengths must be between 2 and 4'
        cdef int order = len(tx)
        self.ptr = NULL
        if order == 2:
            self.ptr = new PSQuadTraction2(surf, <TacsScalar*>tx.data,
                                           <TacsScalar*>ty.data)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new PSQuadTraction3(surf, <TacsScalar*>tx.data,
                                           <TacsScalar*>ty.data)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new PSQuadTraction4(surf, <TacsScalar*>tx.data,
                                           <TacsScalar*>ty.data)
            self.ptr.incref()
        elif order == 5:
            self.ptr = new PSQuadTraction5(surf, <TacsScalar*>tx.data,
                                           <TacsScalar*>ty.data)
            self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

    def numNodes(self):
        return self.ptr.numNodes()

cdef class Traction3D(Element):
    def __cinit__(self, int order, int surf,
                  TacsScalar tx, TacsScalar ty, TacsScalar tz,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] box=None):
        if order < 2 or order > 5:
            errmsg = 'Traction3D order must be between 2 and 5'
            raise ValueError(errmsg)
        if surf < 0 or surf >= 6:
            errmsg = 'Traction3D surf must be between 0 and 5'
            raise ValueError(errmsg)
        self.ptr = NULL
        if order == 2:
            if box is None:
                self.ptr = new TACS3DTraction2(surf, tx, ty, tz)
            else:
                self.ptr = new TACS3DBoundingTraction2(surf, tx, ty, tz,
                                                       <TacsScalar*>box.data)

            self.ptr.incref()
        elif order == 3:
            self.ptr = new TACS3DTraction3(surf, tx, ty, tz)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new TACS3DTraction4(surf, tx, ty, tz)
            self.ptr.incref()
        elif order == 5:
            self.ptr = new TACS3DTraction5(surf, tx, ty, tz)
            self.ptr.incref()
        return
    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return
    def numNodes(self):
        return self.ptr.numNodes()

cdef class PressureTraction3D(Element):
    def __cinit__(self, int order, int surf,
                  TacsScalar pressure):
        if order < 2 or order > 5:
            errmsg = 'PressureTraction3D order must be between 2 and 5'
            raise ValueError(errmsg)
        if surf < 0 or surf >= 6:
            errmsg = 'PressureTraction3D surf must be between 0 and 5'
            raise ValueError(errmsg)
        self.ptr = NULL
        if order == 2:
            self.ptr = new TACS3DPressureTraction2(surf, pressure)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new TACS3DPressureTraction3(surf, pressure)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new TACS3DPressureTraction4(surf, pressure)
            self.ptr.incref()
        elif order == 5:
            self.ptr = new TACS3DPressureTraction5(surf, pressure)
            self.ptr.incref()
        return
    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return
    def numNodes(self):
        return self.ptr.numNodes()

cdef class PlaneTri6(Element):
    def __cinit__(self, PlaneStress stiff,
                  ElementBehaviorType elem_type=LINEAR,
                  int component_num=0):
        '''
        Wrap the PlaneStressTri6 element class
        '''
        cdef PlaneStressStiffness *con = _dynamicPlaneStress(stiff.ptr)
        self.ptr = new PlaneStressTri6(con, elem_type, component_num)
        self.ptr.incref()
        return

    def __dealloc__(self):
        self.ptr.decref()
        return

    def numNodes(self):
        return self.ptr.numNodes()

cdef void _shell_evalf(void *_self, const TacsScalar *X, TacsScalar *T):
    cdef list Xp
    Xp = [X[0], X[1], X[2]]
    t = (<object>_self)(Xp)
    T[0] = t[0]
    T[1] = t[1]
    T[2] = t[2]
    return

cdef class ShellTraction(Element):
    def __cinit__(self, int order,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] tx=None,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] ty=None,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] tz=None,
                  evalf=None):
        self.ptr = NULL
        if order < 2 or order > 5:
            errmsg = 'ShellTraction order must be between 2 and 4'
            raise ValueError(errmsg)
        if evalf is not None:
            if order == 2:
                self.ptr = new TACSShellTraction2(<void*>evalf, _shell_evalf)
            elif order == 3:
                self.ptr = new TACSShellTraction3(<void*>evalf, _shell_evalf)
            elif order == 4:
                self.ptr = new TACSShellTraction4(<void*>evalf, _shell_evalf)
            elif order == 5:
                self.ptr = new TACSShellTraction5(<void*>evalf, _shell_evalf)
        else:
            if order == 2:
                self.ptr = new TACSShellTraction2(<TacsScalar*>tx.data,
                                                  <TacsScalar*>ty.data,
                                                  <TacsScalar*>tz.data)
            elif order == 3:
                self.ptr = new TACSShellTraction3(<TacsScalar*>tx.data,
                                                  <TacsScalar*>ty.data,
                                                  <TacsScalar*>tz.data)
            elif order == 4:
                self.ptr = new TACSShellTraction4(<TacsScalar*>tx.data,
                                                  <TacsScalar*>ty.data,
                                                  <TacsScalar*>tz.data)
            elif order == 5:
                self.ptr = new TACSShellTraction5(<TacsScalar*>tx.data,
                                                  <TacsScalar*>ty.data,
                                                  <TacsScalar*>tz.data)
        if self.ptr is not NULL:
            self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

    def numNodes(self):
        return self.ptr.numNodes()

cdef class MITCShell(Element):
    def __cinit__(self, int order, FSDT stiff,
                  ElementBehaviorType elem_type=LINEAR,
                  int component_num=0, int tying_order=-1):
        '''
        Wrap the MITCShell element class for order 2, 3, 4, and 5
        '''
        cdef FSDTStiffness *con = _dynamicFSDT(stiff.ptr)
        self.ptr = NULL
        if order == 2:
            self.ptr = new MITCShell2(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 3:
            if tying_order == 2:
                self.ptr = new MITCShell32(con, elem_type, component_num)
            else:
                self.ptr = new MITCShell3(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 4:
            if tying_order == 3:
                self.ptr = new MITCShell43(con, elem_type, component_num)
            else:
                self.ptr = new MITCShell4(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 5:
            if tying_order == 4:
                self.ptr = new MITCShell54(con, elem_type, component_num)
            else:
                self.ptr = new MITCShell5(con, elem_type, component_num)
            self.ptr.incref()

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

    def numNodes(self):
        return self.ptr.numNodes()

cdef class Solid(Element):
    def __cinit__(self, int order, SolidStiff stiff,
                  ElementBehaviorType elem_type=LINEAR,
                  int component_num=0):
        '''
        Wrap the Solid element class for order 2,3,4
        '''
        cdef SolidStiffness *con = _dynamicSolid(stiff.ptr)
        if order == 2:
            self.ptr = new Solid2(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new Solid3(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new Solid4(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 5:
            self.ptr = new Solid5(con, elem_type, component_num)
            self.ptr.incref()

    def __dealloc__(self):
        self.ptr.decref()
        return

    def numNodes(self):
        return self.ptr.numNodes()

cdef class MITC(Element):
    def __cinit__(self, FSDT stiff, GibbsVector gravity=None,
                  GibbsVector vInit=None, GibbsVector omegaInit=None):
        cdef FSDTStiffness *con = _dynamicFSDT(stiff.ptr)
        if omegaInit is not None:
            self.ptr = new MITC9(con, gravity.ptr,
                                 vInit.ptr, omegaInit.ptr)
        elif vInit is not None:
            self.ptr = new MITC9(con, gravity.ptr, vInit.ptr, NULL)
        elif gravity is not None:
            self.ptr = new MITC9(con, gravity.ptr, NULL, NULL)
        else:
            self.ptr = new MITC9(con, NULL, NULL, NULL)
        self.ptr.incref()
        return

    def __dealloc__(self):
        self.ptr.decref()
        return

    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return

    def numNodes(self):
        return self.ptr.numNodes()

cdef class MITCBeam(Element):
    def __cinit__(self, Timoshenko stiff, GibbsVector gravity=None,
                  GibbsVector vInit=None, GibbsVector omegaInit=None):
        cdef TimoshenkoStiffness *con = _dynamicTimoshenko(stiff.ptr)
        if omegaInit is not None:
            self.ptr = new MITC3(con, gravity.ptr,
                                 vInit.ptr, omegaInit.ptr)
        elif vInit is not None:
            self.ptr = new MITC3(con, gravity.ptr, vInit.ptr, NULL)
        elif gravity is not None:
            self.ptr = new MITC3(con, gravity.ptr, NULL, NULL)
        else:
            self.ptr = new MITC3(con, NULL, NULL, NULL)
        self.ptr.incref()
        return

    def __dealloc__(self):
        self.ptr.decref()
        return

    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return

    def numNodes(self):
        return self.ptr.numNodes()

cdef TacsScalar _poisson_evalf(void *_self, const TacsScalar *X):
    cdef list Xp
    Xp = [X[0], X[1], X[2]]
    return (<object>_self)(Xp)

cdef class PoissonQuad(Element):
    def __cinit__(self, int order,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] fx=None,
                  evalf=None):
        self.ptr = NULL
        if order < 2 or order > 5:
            errmsg = 'PoissonQuad order must be between 2 and 4'
            raise ValueError(errmsg)
        if evalf is not None:
            if order == 2:
                self.ptr = new PoissonQuad2(<void*>evalf, _poisson_evalf)
            elif order == 3:
                self.ptr = new PoissonQuad3(<void*>evalf, _poisson_evalf)
            elif order == 4:
                self.ptr = new PoissonQuad4(<void*>evalf, _poisson_evalf)
            elif order == 5:
                self.ptr = new PoissonQuad5(<void*>evalf, _poisson_evalf)
        else:
            if order == 2:
                self.ptr = new PoissonQuad2(<TacsScalar*>fx.data)
            elif order == 3:
                self.ptr = new PoissonQuad3(<TacsScalar*>fx.data)
            elif order == 4:
                self.ptr = new PoissonQuad4(<TacsScalar*>fx.data)
            elif order == 5:
                self.ptr = new PoissonQuad5(<TacsScalar*>fx.data)
        if self.ptr is not NULL:
            self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr is not NULL:
            self.ptr.decref()
        return

cdef class PSThermoelasticQuad(Element):
    def __cinit__(self, int order, CoupledPlaneStress stiff,
                  ElementBehaviorType elem_type=LINEAR,
                  int component_num=0):
        '''
        Wrap the PSThermoQuad element class for order 2,3,4
        '''
        cdef CoupledThermoPlaneStressStiffness *con = _dynamicPSThermo(stiff.ptr)
        if order == 2:
            self.ptr = new PSThermoQuad2(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new PSThermoQuad3(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new PSThermoQuad4(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 5:
            self.ptr = new PSThermoQuad5(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 6:
            self.ptr = new PSThermoQuad6(con, elem_type, component_num)
            self.ptr.incref()
        else:
            print("Order %d not implemented"%(order))
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

cdef class PSThermoQuadTraction(Element):
    def __cinit__(self, int surf,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] tx,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] ty):
        if len(tx) != len(ty):
            errmsg = 'Traction lengths must be equal'
            raise ValueError(errmsg)
        if len(tx) < 2 or len(tx) > 7:
            errmsg = 'Traction lengths must be between 2 and 7'
        cdef int order = len(tx)
        self.ptr = NULL
        if order == 2:
            self.ptr = new PSThermoQuadTraction2(surf, <TacsScalar*>tx.data,
                                                 <TacsScalar*>ty.data)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new PSThermoQuadTraction3(surf, <TacsScalar*>tx.data,
                                                 <TacsScalar*>ty.data)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new PSThermoQuadTraction4(surf, <TacsScalar*>tx.data,
                                                 <TacsScalar*>ty.data)
            self.ptr.incref()
        elif order == 5:
            self.ptr = new PSThermoQuadTraction5(surf, <TacsScalar*>tx.data,
                                                 <TacsScalar*>ty.data)
            self.ptr.incref()
        elif order == 6:
            self.ptr = new PSThermoQuadTraction6(surf, <TacsScalar*>tx.data,
                                                 <TacsScalar*>ty.data)
            self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

    def numNodes(self):
        return self.ptr.numNodes()

cdef class PSThermoQuadHF(Element):
    def __cinit__(self, int surf,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] tx,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] ty):
        if len(tx) != len(ty):
            errmsg = 'Traction lengths must be equal'
            raise ValueError(errmsg)
        if len(tx) < 2 or len(tx) > 5:
            errmsg = 'Traction lengths must be between 2 and 4'
        cdef int order = len(tx)
        self.ptr = NULL
        if order == 2:
            self.ptr = new PSThermoQuadHF2(surf, <TacsScalar*>tx.data,
                                           <TacsScalar*>ty.data)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new PSThermoQuadHF3(surf, <TacsScalar*>tx.data,
                                           <TacsScalar*>ty.data)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new PSThermoQuadHF4(surf, <TacsScalar*>tx.data,
                                           <TacsScalar*>ty.data)
            self.ptr.incref()
        elif order == 5:
            self.ptr = new PSThermoQuadHF5(surf, <TacsScalar*>tx.data,
                                           <TacsScalar*>ty.data)
            self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

    def numNodes(self):
        return self.ptr.numNodes()


cdef class SolidThermo(Element):
    def __cinit__(self, int order, CoupledSolid stiff,
                  ElementBehaviorType elem_type=LINEAR,
                  int component_num=0):
        '''
        Wrap the SolidThermo element class for order 2,3,4
        '''
        cdef CoupledThermoSolidStiffness *con = _dynamicSolidThermo(stiff.ptr)
        if order == 2:
            self.ptr = new SolidThermo2(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new SolidThermo3(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new SolidThermo4(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 5:
            self.ptr = new SolidThermo5(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 6:
            self.ptr = new SolidThermo6(con, elem_type, component_num)
            self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

cdef class TACS3DThermoTraction(Element):
    def __cinit__(self, int order, int surf,
                  TacsScalar tx, TacsScalar ty, TacsScalar tz):
        if order < 2 or order > 7:
            errmsg = 'ThermoTraction3D order must be between 2 and 4'
            raise ValueError(errmsg)
        if surf < 0 or surf >= 6:
            errmsg = 'ThermoTraction3D surf must be between 0 and 5'
            raise ValueError(errmsg)
        self.ptr = NULL
        if order == 2:
            self.ptr = new TACS3DThermoTraction2(surf, tx, ty, tz)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new TACS3DThermoTraction3(surf, tx, ty, tz)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new TACS3DThermoTraction4(surf, tx, ty, tz)
            self.ptr.incref()
        elif order == 5:
            self.ptr = new TACS3DThermoTraction5(surf, tx, ty, tz)
            self.ptr.incref()
        elif order == 6:
            self.ptr = new TACS3DThermoTraction6(surf, tx, ty, tz)
            self.ptr.incref()

        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

    def numNodes(self):
        return self.ptr.numNodes()

cdef class TACS3DThermoPressureTraction(Element):
    def __cinit__(self, int order, int surf,
                  TacsScalar pressure):
        if order < 2 or order > 5:
            errmsg = 'ThermoPressureTraction3D order must be between 2 and 5'
            raise ValueError(errmsg)
        if surf < 0 or surf >= 6:
            errmsg = 'ThermoPressureTraction3D surf must be between 0 and 5'
            raise ValueError(errmsg)
        self.ptr = NULL
        if order == 2:
            self.ptr = new TACS3DThermoPressureTraction2(surf, pressure)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new TACS3DThermoPressureTraction3(surf, pressure)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new TACS3DThermoPressureTraction4(surf, pressure)
            self.ptr.incref()
        elif order == 5:
            self.ptr = new TACS3DThermoPressureTraction5(surf, pressure)
            self.ptr.incref()

        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

    def numNodes(self):
        return self.ptr.numNodes()

cdef class TACS3DThermoHF(Element):
    def __cinit__(self, int order, int surf,
                  TacsScalar tx, TacsScalar ty, TacsScalar tz):
        if order < 2 or order > 5:
            errmsg = 'HeatFlux3D order must be between 2 and 4'
            raise ValueError(errmsg)
        if surf < 0 or surf >= 6:
            errmsg = 'HeatFlux3D surf must be between 0 and 5'
            raise ValueError(errmsg)
        self.ptr = NULL
        if order == 2:
            self.ptr = new TACS3DThermoHF2(surf, tx, ty, tz)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new TACS3DThermoHF3(surf, tx, ty, tz)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new TACS3DThermoHF4(surf, tx, ty, tz)
            self.ptr.incref()
        elif order == 5:
            self.ptr = new TACS3DThermoHF5(surf, tx, ty, tz)
            self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

    def numNodes(self):
        return self.ptr.numNodes()

cdef class TACS3DThermoNormalHF(Element):
    def __cinit__(self, int order, int surf,
                  TacsScalar qn):
        if order < 2 or order > 4:
            errmsg = 'HeatFlux3D order must be between 2 and 4'
            raise ValueError(errmsg)
        if surf < 0 or surf >= 6:
            errmsg = 'HeatFlux3D surf must be between 0 and 5'
            raise ValueError(errmsg)
        self.ptr = NULL
        if order == 2:
            self.ptr = new TACS3DThermoNormalHF2(surf, qn)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new TACS3DThermoNormalHF3(surf, qn)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new TACS3DThermoNormalHF4(surf, qn)
            self.ptr.incref()
        elif order == 5:
            self.ptr = new TACS3DThermoNormalHF5(surf, qn)
            self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

    def numNodes(self):
        return self.ptr.numNodes()

# This wraps a C++ array with a numpy array for later useage
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

cdef inplace_array_2d(int nptype, int dim1, int dim2, void *data_ptr):
    '''Return a numpy version of the array'''
    # Set the shape of the array
    cdef int size = 2
    cdef np.npy_intp shape[2]
    cdef np.ndarray ndarray

    # Set the first entry of the shape array
    shape[0] = <np.npy_intp>dim1
    shape[1] = <np.npy_intp>dim2

    # Create the array itself - Note that this function will not
    # delete the data once the ndarray goes out of scope
    ndarray = np.PyArray_SimpleNewFromData(size, shape,
                                           nptype, data_ptr)

    return ndarray

cdef void getinitconditions(void * _self, int nvars, int num_nodes,
                            TacsScalar * vars,
                            TacsScalar * dvars,
                            TacsScalar * ddvars,
                            const TacsScalar * Xpts):
    '''Get the initial conditions'''
    _vars = inplace_array_1d(TACS_NPY_SCALAR, nvars, <void*>vars)
    _dvars = inplace_array_1d(TACS_NPY_SCALAR, nvars, <void*>dvars)
    _ddvars = inplace_array_1d(TACS_NPY_SCALAR, nvars, <void*>ddvars)
    _Xpts = inplace_array_1d(TACS_NPY_SCALAR, 3*num_nodes, <void*>Xpts)
    (<object>_self).getInitConditions(_vars, _dvars, _ddvars, _Xpts)
    return

cdef void addresidual(void * _self, int nvars, int num_nodes,
                      double time, TacsScalar * res,
                      const TacsScalar * Xpts,
                      const TacsScalar * vars,
                      const TacsScalar * dvars,
                      const TacsScalar * ddvars):
    '''Add the residual'''
    _res = inplace_array_1d(TACS_NPY_SCALAR, nvars, <void*>res)
    _Xpts = inplace_array_1d(TACS_NPY_SCALAR, 3*num_nodes, <void*>Xpts)
    _vars = inplace_array_1d(TACS_NPY_SCALAR, nvars, <void*>vars)
    _dvars = inplace_array_1d(TACS_NPY_SCALAR, nvars, <void*>dvars)
    _ddvars = inplace_array_1d(TACS_NPY_SCALAR, nvars, <void*>ddvars)
    (<object>_self).addResidual(time, _res, _Xpts, _vars, _dvars, _ddvars)
    return

cdef void addjacobian(void * _self, int nvars, int num_nodes,
                      double time, TacsScalar J[],
                      double alpha, double beta, double gamma,
                      const TacsScalar Xpts[],
                      const TacsScalar vars[],
                      const TacsScalar dvars[],
                      const TacsScalar ddvars[]):
    '''Add the Jacobian'''
    _J = inplace_array_2d(TACS_NPY_SCALAR, nvars, nvars, <void*>J)
    _Xpts = inplace_array_1d(TACS_NPY_SCALAR, 3*num_nodes, <void*>Xpts)
    _vars = inplace_array_1d(TACS_NPY_SCALAR, nvars, <void*>vars)
    _dvars = inplace_array_1d(TACS_NPY_SCALAR, nvars, <void*>dvars)
    _ddvars = inplace_array_1d(TACS_NPY_SCALAR, nvars, <void*>ddvars)
    (<object>_self).addJacobian(time, _J, alpha, beta, gamma,
                                _Xpts, _vars, _dvars, _ddvars)
    return

cdef void addadjresproduct(void * _self, int nvars, int num_nodes, double time,
                           TacsScalar scale, TacsScalar dvSens[], int dvLen,
                           const TacsScalar psi[], const TacsScalar Xpts[],
                           const TacsScalar vars[], const TacsScalar dvars[],
                           const TacsScalar ddvars[]):
    '''Add the derivative of the adjoint-residual product'''
    _dvSens = inplace_array_1d(TACS_NPY_SCALAR, dvLen, <void*>dvSens)
    _psi = inplace_array_1d(TACS_NPY_SCALAR, nvars, <void*>psi)
    _Xpts = inplace_array_1d(TACS_NPY_SCALAR, 3*num_nodes, <void*>Xpts)
    _vars = inplace_array_1d(TACS_NPY_SCALAR, nvars, <void*>vars)
    _dvars = inplace_array_1d(TACS_NPY_SCALAR, nvars, <void*>dvars)
    _ddvars = inplace_array_1d(TACS_NPY_SCALAR, nvars, <void*>ddvars)

    (<object>_self).addAdjResProduct(time, scale, _dvSens, _psi, _Xpts,
                                     _vars, _dvars, _ddvars)

    return

cdef class pyElement(Element):
    def __cinit__(self, int num_nodes, int num_displacements, *args, **kwargs):
        cdef TACSElementWrapper *pointer
        pointer = new TACSElementWrapper(<PyObject*>self,
                                         num_nodes, num_displacements)
        pointer.incref()

        # Set the function pointers
        pointer.getinitconditions = getinitconditions
        pointer.addresidual = addresidual
        pointer.addjacobian = addjacobian
        pointer.addadjresproduct = addadjresproduct

        self.ptr = pointer

    def __dealloc__(self):
        self.ptr.decref()
        return

    def addAdjResProduct(self, time, scale, dvSens, psi,
                         Xpts, vars, dvars, ddvars):
        '''Raise error if this function is called and hasn\'t been overrided'''
        raise NotImplementedError()
