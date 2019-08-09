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

# Import the definition required for const strings
from libc.string cimport const_char

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

# Import the definitions
from TACS cimport *
from constitutive cimport *

# Include the definitions
include "TacsDefs.pxi"

# Include the mpi4py header
cdef extern from "mpi-compat.h":
    pass

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

cdef class Timoshenko(Constitutive):
    ## def __cinit__(self,
    ##               np.ndarray[TacsScalar, ndim=1, mode='c'] axis,
    ##               EA, EIy, EIz, EIyz,
    ##               GJ, kGy, kGz, kGyz,
    ##               rho, Iyy, Izz, Iyz,
    ##               xm2, xm3,
    ##               xc2, xc3,
    ##               xk2, xk3,
    ##               muS):        
    ##     self.ptr = new TimoshenkoStiffness(<TacsScalar*>axis.data,
    ##                                        EA, EIy, EIz, EIyz,
    ##                                        GJ, kGy, kGz, kGyz,
    ##                                        rho, Iyy, Izz, Iyz,
    ##                                        xm2, xm3,
    ##                                        xc2, xc3,
    ##                                        xk2, xk3,
    ##                                        muS)
    ##     return
                  
    def __cinit__(self, rhoA, rhoIy, rhoIz, rhoIyz,
                      EA, GJ, EIy, EIz, kGAy, kGAz,
                      np.ndarray[TacsScalar, ndim=1, mode='c'] axis):
        self.ptr = new TimoshenkoStiffness(rhoA, rhoIy, rhoIz, rhoIyz,
                                           EA, GJ, EIy, EIz, kGAy, kGAz,
                                           <TacsScalar*>axis.data)
        return

cdef class FSDT(Constitutive):
    def __cinit__(self, *args, **kwargs):
        self.ptr = NULL
        return

    def setRefAxis(self, np.ndarray[TacsScalar, ndim=1, mode='c'] axis):
        cdef FSDTStiffness *stiff = NULL
        stiff = _dynamicFSDT(self.ptr)
        if stiff:
            stiff.setRefAxis(<TacsScalar*>axis.data)
        return

    def printStiffness(self):
        stiff = _dynamicFSDT(self.ptr)
        if stiff:
            stiff.printStiffness()
        return

    def setDrillingRegularization(self, double krel):
        cdef FSDTStiffness *stiff = NULL
        stiff = _dynamicFSDT(self.ptr)
        if stiff:
            stiff.setDrillingRegularization(krel)
        return

cdef class isoFSDT(FSDT):
    def __cinit__(self, rho, E, nu, kcorr, ys, t, tNum, minT, maxT):
        '''
        Wraps the isoFSDTStiffness class that is used with shell elements
        '''
        self.ptr = new isoFSDTStiffness(rho, E, nu, kcorr, ys, 
                                        t, tNum, minT, maxT)
        self.ptr.incref()
        return

cdef class PlaneStress(Constitutive):
    def __cinit__(self, *args, **kwargs):
        self.ptr = NULL
        return

cdef class CoupledPlaneStress(Constitutive):
    def __cinit__(self, *args, **kwargs):
        self.ptr = NULL
        return

cdef class CoupledSolid(Constitutive):
    def __cinit__(self, *args, **kwargs):
        self.ptr = NULL
        return
    
cdef class SimplePlaneStress(PlaneStress):
    def __cinit__(self, *args, **kwargs):
        '''
        Wraps the PlaneStressStiffness class that is used with 2D elements
        '''
        if len(args) == 3:
            rho = args[0]
            E = args[1]
            nu = args[2]
            self.ptr = new PlaneStressStiffness(rho, E, nu)
        else:
            self.ptr = new PlaneStressStiffness()
            
        self.ptr.incref()
        return

cdef class SolidStiff(Constitutive):
    def __cinit__(self, *args, **kwargs):
        '''
        Wraps the SolidStiffness class that is used with 3D elements
        '''
        self.ptr = NULL
        return
    
cdef class isoSolidStiff(SolidStiff):
    def __cinit__(self, rho, E, nu, ys=1.0, eNum=-1, *args, **kwargs):
        '''
        Wraps the SolidStiffness class that is used with 3D elements
        '''
        self.ptr = new SolidStiffness(rho, E, nu, ys, eNum)
        self.ptr.incref()
        return

cdef void getdesignvars(void *_self, TacsScalar *x, int dvLen):
    '''Get the design variable values'''
    xvals = inplace_array_1d(np.NPY_DOUBLE, dvLen, <void*>x)
    (<object>_self).getDesignVars(xvals)
    return

cdef void setdesignvars(void *_self, const TacsScalar *x, int dvLen):
    '''Set the design variable values'''
    cdef np.ndarray xvals = np.zeros(dvLen)
    for i in range(dvLen):
        xvals[i] = x[i]
    (<object>_self).setDesignVars(xvals)
    return

cdef void getdesignvarrange(void *_self, TacsScalar *lb,
                            TacsScalar *ub, int dvLen):
    '''Get the design variable range'''
    xlb = inplace_array_1d(np.NPY_DOUBLE, dvLen, <void*>lb)
    xub = inplace_array_1d(np.NPY_DOUBLE, dvLen, <void*>ub)
    (<object>_self).getDesignVarRange(xlb, xub)
    return

# Callbacks for the python-level implementation of the
# PlaneStressStiffness object.
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

# Python-level interface for the plane stress constitutive object
cdef class pyPlaneStress(PlaneStress):
    def __cinit__(self, *args, **kwargs):
        cdef PSStiffnessWrapper *pointer
        pointer = new PSStiffnessWrapper(<PyObject*>self)
        pointer.incref()
        
        # Set the function pointers
        pointer.setdesignvars = setdesignvars
        pointer.getdesignvars = getdesignvars
        pointer.getdesignvarrange = getdesignvarrange
        pointer.calculatestress = ps_calculatestress
        pointer.addstressdvsens = ps_addstressdvsens
        pointer.getpointwisemass = ps_getpointwisemass
        pointer.addpointwisemassdvsens = ps_addpointwisemassdvsens
        pointer.fail = ps_fail
        pointer.failstrainsens = ps_failstrainsens
        pointer.addfaildvsens = ps_addfaildvsens

        self.ptr = pointer
        return

    def __dealloc__(self):
        self.ptr.decref()
        return

    def setDesignVars(self, x):
        return

    def getDesignVars(self, x):
        return

    def getDesignVarRange(self, lb, ub):
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

# Callbacks for the python-level implementation of the
# PlaneStressStiffness object.
cdef TacsScalar fsdt_getstiffness(void *_self, const double *pt,
                                  TacsScalar *a, TacsScalar *b,
                                  TacsScalar *d, TacsScalar *as):
    # Extract the parametric point
    p = np.zeros(2)
    p[0] = pt[0]
    p[1] = pt[1]

    # Get the numpy arrays for the stiffness object
    A, B, D, As, rot = (<object>_self).getStiffness(p)

    # Copy over the stiffness values
    for i in range(6):
        a[i] = A[i]
        b[i] = B[i]
        d[i] = D[i]

    # Set the out-of-plane shear stiffness
    as[0] = As[0]
    as[1] = As[1]
    as[2] = As[2]
    
    return rot

cdef void fsdt_addstiffnessdvsens(void *_self, const double *pt,
                                  const TacsScalar *strain,
                                  const TacsScalar *psi,
                                  TacsScalar rotPsi,
                                  TacsScalar *fdvSens, int dvLen):
    # Extract the strain
    cdef np.ndarray e = np.zeros(8)
    for i in range(8):
        e[i] = strain[i]
    cdef np.ndarray ps = np.zeros(8)
    for i in range(8):
        ps[i] = psi[i]
        
    # Extract the parametric point
    cdef np.ndarray p = np.zeros(2)
    p[0] = pt[0]
    p[1] = pt[1]

    # Wrap the python array
    fdvs = inplace_array_1d(TACS_NPY_SCALAR, dvLen, <void*>fdvSens)

    # Add the derivative to the array
    (<object>_self).addStiffnessDVSens(p, e, ps, rotPsi, fdvs)

    return

cdef void fsdt_getpointwisemass(void *_self, const double *pt,
                                TacsScalar *mass):
    # Extract the parametric point
    cdef np.ndarray p = np.zeros(2)
    p[0] = pt[0]
    p[1] = pt[1]

    rho = (<object>_self).getPointwiseMass(p)
    mass[0] = rho[0]
    mass[1] = rho[1]

    return

cdef void fsdt_addpointwisemassdvsens(void *_self, const double *pt,
                                      const TacsScalar *alpha,
                                      TacsScalar *fdvSens, int dvLen):
    # Extract the parametric point
    cdef np.ndarray p = np.zeros(2)
    p[0] = pt[0]
    p[1] = pt[1]

    # Set the scalar
    cdef np.ndarray avals = np.zeros(2)
    avals[0] = alpha[0]
    avals[1] = alpha[1]

    # Wrap the python array
    fdvs = inplace_array_1d(TACS_NPY_SCALAR, dvLen, <void*>fdvSens)

    (<object>_self).addPointwiseMassDVSens(p, avals, fdvs)
    
    return

cdef TacsScalar fsdt_fail(void *_self, const double *pt,
                          const TacsScalar *strain):
    cdef TacsScalar fval = 0.0

    # Extract the strain
    cdef np.ndarray e = np.zeros(8)
    for i in range(8):
        e[i] = strain[i]
        
    # Extract the parametric point
    cdef np.ndarray p = np.zeros(2)
    p[0] = pt[0]
    p[1] = pt[1]
    
    fval = (<object>_self).failure(p, e)

    return fval

cdef void fsdt_failstrainsens(void *_self, const double *pt,
                              const TacsScalar *strain,
                              TacsScalar *sens):
    # Extract the strain
    cdef np.ndarray e = np.zeros(8)
    for i in range(8):
        e[i] = strain[i]

    # Extract the parametric point
    cdef np.ndarray p = np.zeros(2)
    p[0] = pt[0]
    p[1] = pt[1]
    
    fsens = (<object>_self).failureStrainSens(p, e)

    for i in range(8):
        sens[i] = fsens[i]

    return

cdef void fsdt_addfaildvsens(void *_self, const double *pt,
                             const TacsScalar *strain, 
                             TacsScalar alpha,
                             TacsScalar *fdvSens, int dvLen):
    # Extract the strain
    cdef np.ndarray e = np.zeros(8)
    for i in range(8):
        e[i] = strain[i]
    
    # Extract the parametric point
    cdef np.ndarray p = np.zeros(2)
    p[0] = pt[0]
    p[1] = pt[1]

    # Wrap the python array
    fdvs = inplace_array_1d(TACS_NPY_SCALAR, dvLen, <void*>fdvSens)
    
    (<object>_self).addFailureDVSens(p, e, alpha, fdvs)

    return

# Python-level interface for the plane stress constitutive object
cdef class pyFSDT(FSDT):
    def __cinit__(self, *args, **kwargs):
        cdef FSDTStiffnessWrapper *pointer
        pointer = new FSDTStiffnessWrapper(<PyObject*>self)
        pointer.incref()
        
        # Set the function pointers
        pointer.setdesignvars = setdesignvars
        pointer.getdesignvars = getdesignvars
        pointer.getdesignvarrange = getdesignvarrange
        pointer.getstiffness = fsdt_getstiffness
        pointer.addstiffnessdvsens = fsdt_addstiffnessdvsens
        pointer.getpointwisemass = fsdt_getpointwisemass
        pointer.addpointwisemassdvsens = fsdt_addpointwisemassdvsens
        pointer.fail = fsdt_fail
        pointer.failstrainsens = fsdt_failstrainsens
        pointer.addfaildvsens = fsdt_addfaildvsens

        self.ptr = pointer
        return

    def __dealloc__(self):
        self.ptr.decref()
        return

    def setDesignVars(self, x):
        return

    def getDesignVars(self, x):
        return

    def getDesignVarRange(self, lb, ub):
        return

    def getStiffness(self, pt):
        '''Return the A, B, D, As matrices and the rotational stiffness'''
        return np.zeros(6), np.zeros(6), np.zeros(6), np.zeros(3), 0.0

    def addStiffnessDVSens(self, pt, e, psi, rot, fdvs):
        '''Add the derivative of the inner product to the fdvs array'''
        return

    def getPointwiseMass(self, pt):
        return np.zeros(2)

    def addPointwiseMassDVSens(p, alpha, fdvs):
        return

    def failure(self, pt, e):
        return 0.0

    def failureStrainSens(self, pt, e):
        return np.zeros(8)

    def addFailureDVSens(self, pt, e, alpha, fdvs):
        return
