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

cdef extern from "mpi-compat.h":
   pass

# A generic wrapper class for the TACSElement object
cdef class Element:
   def __cinit__(self):
       self.ptr = NULL
       return

cdef class TestElem:
   cdef TestElement *ptr
   def __cinit__(self, Element elem):
      self.ptr = new TestElement(elem.ptr, NULL)
      self.ptr.incref()

   def __dealloc__(self):
      self.ptr.decref()
      return

   def setFailTolerances(self, double rtol, double atol):
      self.ptr.setFailTolerances(rtol, atol)
      return

   def setPrintLevel(self, int lev):
      self.ptr.setPrintLevel(lev)
      return

   def setStepSize(self, double dh):
      self.ptr.setStepSize(dh)
      return

   def testResidual(self):
      return self.ptr.testResidual()

   def testJacobian(self, int col=-1):
      return self.ptr.testJacobian(col)

   def testAdjResProduct(self, np.ndarray[TacsScalar, ndim=1, mode='c'] x):
      return self.ptr.testAdjResProduct(<TacsScalar*>x.data, len(x))

   def testStrainSVSens(self, np.ndarray[double, ndim=1, mode='c'] pt):
      return self.ptr.testStrainSVSens(<double*>pt.data)

   def testJacobianXptSens(self, np.ndarray[double, ndim=1, mode='c'] pt):
      return self.ptr.testJacobianXptSens(<double*>pt.data)

cdef class GibbsVector:
   cdef TACSGibbsVector *ptr
   def __cinit__(self, np.ndarray[TacsScalar, ndim=1, mode='c'] x):
      assert(len(x) == 3)
      self.ptr = new TACSGibbsVector(<TacsScalar*>x.data)
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
      assert(len(tx) == len(ty))
      cdef int order = len(tx)
      if order == 2:
         self.ptr = new PSQuadTraction2(surf, <TacsScalar*>tx.data,
                                        <TacsScalar*>ty.data)
         self.ptr.incref()
      elif order == 3:
         self.ptr = new PSQuadTraction3(surf, <TacsScalar*>tx.data,
                                        <TacsScalar*>ty.data)
         self.ptr.incref()
      elif order == 4:
         self.ptr = new PSQuadTraction2(surf, <TacsScalar*>tx.data,
                                        <TacsScalar*>ty.data)
         self.ptr.incref()
      return

   def __dealloc__(self):
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

cdef class MITCShell(Element):
   def __cinit__(self, int order, FSDT stiff, ElementBehaviorType elem_type=LINEAR,
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
   def __cinit__(self, int order, solid stiff, ElementBehaviorType elem_type=LINEAR,
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
