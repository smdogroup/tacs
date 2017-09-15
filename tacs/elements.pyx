# For the use of MPI
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy 
import numpy as np
cimport numpy as np

# Ensure that numpy is initialized
np.import_array()

# Import the definition required for const strings
from libc.string cimport const_char

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
    def __cinit__(self, GibbsVector orig, GibbsVector rev, TacsScalar omega):
        self.ptr = new TACSRevoluteDriver(orig.ptr, rev.ptr, omega)
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
        if len(tx) < 2 or len(tx) > 4:
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
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

cdef class Traction3D(Element):
    def __cinit__(self, int order, int surf, 
                  TacsScalar tx, TacsScalar ty, TacsScalar tz):
        if order < 2 or order > 4:
            errmsg = 'Traction3D order must be between 2 and 4'
            raise ValueError(errmsg)
        if surf < 0 or surf >= 6:
            errmsg = 'Traction3D surf must be between 0 and 5'
            raise ValueError(errmsg)
        self.ptr = NULL
        if order == 2:
            self.ptr = new TACS3DTraction2(surf, tx, ty, tz)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new TACS3DTraction3(surf, tx, ty, tz)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new TACS3DTraction4(surf, tx, ty, tz)
            self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

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

cdef class ShellTraction(Element):
    def __cinit__(self, int order, 
                  TacsScalar tx, TacsScalar ty, TacsScalar tz):
        if order < 2 or order > 4:
            errmsg = 'ShellTraction order must be between 2 and 4'
            raise ValueError(errmsg)
        self.ptr = NULL
        if order == 2:
            self.ptr = new TACSShellTraction2(tx, ty, tz)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new TACSShellTraction3(tx, ty, tz)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new TACSShellTraction4(tx, ty, tz)
            self.ptr.incref()
        return

    def __dealloc__(self):
        self.ptr.decref()
        return

cdef class MITCShell(Element):
    def __cinit__(self, int order, FSDT stiff, 
                  ElementBehaviorType elem_type=LINEAR,
                  int component_num=0):
        '''
        Wrap the MITCShell element class for order 2,3,4
        '''
        cdef FSDTStiffness *con = _dynamicFSDT(stiff.ptr)
        if order == 2:
            self.ptr = new MITCShell2(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new MITCShell3(con, elem_type, component_num)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new MITCShell4(con, elem_type, component_num)
            self.ptr.incref()
            
    def __dealloc__(self):
        self.ptr.decref()
        return

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
            
    def __dealloc__(self):
        self.ptr.decref()
        return
        
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
