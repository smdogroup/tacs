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

# This class wraps a C++ array with a numpy array for later useage
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
   def __cinit__(self, *args, **kwargs):
      self.ptr = NULL
      return

cdef class isoFSDT(FSDT):
   def __cinit__(self, rho, E, nu, kcorr, ys, t, tNum, minT, maxT):
      '''
      Wraps the isoFSDTStiffness class that is used with shell elements
      '''
      self.ptr = new isoFSDTStiffness(rho, E, nu, kcorr, ys, t, tNum, minT, maxT)
      self.ptr.incref()
      return

cdef class PlaneStress(Constitutive):
   def __cinit__(self, rho, E, nu, *args, **kwargs):
      '''
      Wraps the PlaneStressStiffness class that is used with 2D elements
      '''
      self.ptr = new PlaneStressStiffness(rho, E, nu)
      self.ptr.incref()
      return
   
cdef class solid(Constitutive):
   def __cinit__(self, rho, E, nu, *args, **kwargs):
      '''
      Wraps the SolidStiffness class that is used with 3D elements
      '''
      self.ptr = new SolidStiffness(rho, E, nu)
      self.ptr.incref()
      return

# Wrapper for the PlaneStressStiffness object
cdef void ps_calculatestress(void *_self, const double *pt,
                             const TacsScalar *strain,
                             TacsScalar *stress):
   # Extract the strain
   e = np.zeros(3)
   e[0] = strain[0]
   e[1] = strain[1]
   e[2] = strain[2]

   # Extract the parametric point
   p = np.zeros(2)
   p[0] = pt[0]
   p[1] = pt[1]
   
   s = (<object>_self).calculateStress(p, e)

   # Conver the stress values
   stress[0] = s[0]
   stress[1] = s[1]
   stress[2] = s[2]
   return

cdef void ps_addstressdvsens(void *_self, const double *pt,
                             const TacsScalar *strain,
                             TacsScalar alpha, const TacsScalar *psi,
                             TacsScalar *fdvSens, int dvLen):
   # Extract the strain
   e = np.zeros(3)
   e[0] = strain[0]
   e[1] = strain[1]
   e[2] = strain[2]

   # Extract the sensitivity input vector
   ps = np.zeros(3)
   ps[0] = psi[0]
   ps[1] = psi[1]
   ps[2] = psi[2]

   # Extract the parametric point
   p = np.zeros(2)
   p[0] = pt[0]
   p[1] = pt[1]

   # Wrap the python array
   fdvs = inplace_array_1d(np.NPY_DOUBLE, dvLen, <void*>fdvSens)

   s = (<object>_self).addStressDVSens(p, e, alpha, ps, fdvs)

   return

cdef void ps_getpointwisemass(void *_self, const double *pt,
                              TacsScalar *mass):
   # Extract the parametric point
   p = np.zeros(2)
   p[0] = pt[0]
   p[1] = pt[1]

   mass[0] = (<object>_self).getPointwiseMass(p)
   return

cdef void ps_addpointwisemassdvsens(void *_self, const double *pt,
                                    const TacsScalar *alpha,
                                    TacsScalar *fdvSens, int dvLen):
   # Extract the parametric point
   p = np.zeros(2)
   p[0] = pt[0]
   p[1] = pt[1]

   # Set the scalar
   cdef TacsScalar ascalar = alpha[0]

   # Wrap the python array
   fdvs = inplace_array_1d(np.NPY_DOUBLE, dvLen, <void*>fdvSens)

   (<object>_self).addPointwiseMassDVSens(p, ascalar, fdvs)
   
   return

cdef TacsScalar ps_fail(void *_self, const double *pt,
                        const TacsScalar *strain):
   cdef TacsScalar fval = 0.0

   # Extract the strain
   e = np.zeros(3)
   e[0] = strain[0]
   e[1] = strain[1]
   e[2] = strain[2]

   # Extract the parametric point
   p = np.zeros(2)
   p[0] = pt[0]
   p[1] = pt[1]
   
   fval = (<object>_self).failure(p, e)

   return fval

cdef void ps_failstrainsens(void *_self, const double *pt,
                            const TacsScalar *strain,
                            TacsScalar *sens):
   # Extract the strain
   e = np.zeros(3)
   e[0] = strain[0]
   e[1] = strain[1]
   e[2] = strain[2]

   # Extract the parametric point
   p = np.zeros(2)
   p[0] = pt[0]
   p[1] = pt[1]
   
   fsens = (<object>_self).failureStrainSens(p, e)

   sens[0] = fsens[0]
   sens[1] = fsens[1]
   sens[2] = fsens[2]
   return

cdef void ps_addfaildvsens(void *_self, const double *pt,
                           const TacsScalar *strain, 
                           TacsScalar alpha,
                           TacsScalar *fdvSens, int dvLen):
   # Extract the strain
   e = np.zeros(3)
   e[0] = strain[0]
   e[1] = strain[1]
   e[2] = strain[2]

   # Extract the parametric point
   p = np.zeros(2)
   p[0] = pt[0]
   p[1] = pt[1]

   # Wrap the python array
   fdvs = inplace_array_1d(np.NPY_DOUBLE, dvLen, <void*>fdvSens)
   
   (<object>_self).addFailureDVSens(p, e, alpha, fdvs)

   return

cdef class pyPlaneStress(PlaneStress):
   def __cinit__(self, rho, E, nu, *args, **kwargs):
      cdef PSStiffnessWrapper *pointer
      pointer = new PSStiffnessWrapper()
      pointer.incref()
      
      # Set the function pointers
      pointer.self_ptr = <void*>self
      pointer.calculatestress = ps_calculatestress
      pointer.addstressdvsens = ps_addstressdvsens
      pointer.getpointwisemass = ps_getpointwisemass
      pointer.addpointwisemassdvsens = ps_addpointwisemassdvsens
      pointer.fail = ps_fail
      pointer.failstrainsens = ps_failstrainsens
      pointer.addfaildvsens = ps_addfaildvsens

      self.ptr = pointer
      return

   def calculateStress(self, pt, e):
      return np.zeros(3)

   def addStressDVSens(self, pt, e, alpha, psi, fdvs):
      return

   def getPointwiseMass(self, pt):
      return 0.0

   def addPointwiseMassDVSens(p, ascale, fdvs):
      return

   def failure(self, pt, e):
      return 0.0

   def failureStrainSens(self, pt, e):
      return np.zeros(3)

   def addFailureDVSens(self, pt, e, alpha, fdvs):
      return
