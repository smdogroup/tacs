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

cdef extern from "mpi-compat.h":
   pass

# A generic wrapper class for the TACSConstitutive object
cdef class Constitutive:
   def __cinit__(self):
      self.ptr = NULL
      return
   def __dealloc__(self):
      if self.ptr:
         self.ptr.decref()
      return

cdef class FSDT(Constitutive):
   def __cinit__(self):
      self.ptr = NULL
      return

cdef class isoFSDT(Constitutive):
   def __cinit__(self, rho, E, nu, kcorr, ys, t, tNum, minT, maxT):
      '''
      Wraps the isoFSDTStiffness class that is used with shell elements
      '''
      self.ptr = new isoFSDTStiffness(rho, E, nu, kcorr, ys, t, tNum, minT, maxT)
      self.ptr.incref()
      return

cdef class PlaneStress(Constitutive):
   def __cinit__(self, rho, E, nu):
      '''
      Wraps the PlaneStressStiffness class that is used with 2D elements
      '''
      self.ptr = new PlaneStressStiffness(rho, E, nu)
      self.ptr.incref()
      return

   ## def getPlaneStress(self):
   ##    return _dynamicPlaneStress(self.ptr)
   
cdef class Solid(Constitutive):
   def __cinit__(self, rho, E, nu):
      '''
      Wraps the SolidStiffness class that is used with 3D elements
      '''
      self.ptr = new SolidStiffness(rho, E, nu)
      self.ptr.incref()
      return
