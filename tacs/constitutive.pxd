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

# For MPI capabilities
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy
import numpy as np
cimport numpy as np
from libc.string cimport const_char

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

# Import from TACS for definitions
from TACS cimport *

cdef extern from "FSDTStiffness.h":
    cdef cppclass FSDTStiffness(TACSConstitutive):
        FSDTStiffness()
        void setRefAxis(TacsScalar*)
        void printStiffness()
        void setDrillingRegularization(double)

cdef extern from "isoFSDTStiffness.h":
    cdef cppclass isoFSDTStiffness(FSDTStiffness):
        isoFSDTStiffness(TacsScalar rho, TacsScalar E, TacsScalar nu, 
                         TacsScalar kcorr, TacsScalar yieldStress, 
                         TacsScalar thickness, int tNum,
                         TacsScalar minThickness, 
                         TacsScalar maxThickness)

cdef extern from "TimoshenkoStiffness.h":
    cdef cppclass TimoshenkoStiffness(TACSConstitutive):
        
        ## TimoshenkoStiffness(const TacsScalar* axis,
        ##                     TacsScalar EA, 
        ##                     TacsScalar EI22, TacsScalar EI33, TacsScalar EI23,
        ##                     TacsScalar GJ,
        ##                     TacsScalar kG22, TacsScalar kG33, TacsScalar kG23,
        ##                     TacsScalar m00,
        ##                     TacsScalar m11, TacsScalar m22, TacsScalar m33,
        ##                     TacsScalar xm2, TacsScalar xm3,
        ##                     TacsScalar xc2, TacsScalar xc3,
        ##                     TacsScalar xk2, TacsScalar xk3,
        ##                     TacsScalar muS);
                       
        TimoshenkoStiffness(TacsScalar, TacsScalar, TacsScalar, TacsScalar,
                            TacsScalar, TacsScalar, TacsScalar, TacsScalar,
                            TacsScalar, TacsScalar, const TacsScalar*)

cdef extern from "PlaneStressStiffness.h":
    cdef cppclass PlaneStressStiffness(TACSConstitutive):
        PlaneStressStiffness()
        PlaneStressStiffness(TacsScalar rho, TacsScalar E, TacsScalar nu)

cdef extern from "CoupledThermoPlaneStressStiffness.h":
    cdef cppclass CoupledThermoPlaneStressStiffness(PlaneStressStiffness):
        CoupledThermoPlaneStressStiffness()
        CoupledThermoPlaneStressStiffness( TacsScalar, TacsScalar, TacsScalar,
                                           TacsScalar, TacsScalar, TacsScalar )
        
cdef extern from "SolidStiffness.h":
    cdef cppclass SolidStiffness(TACSConstitutive):
        SolidStiffness()
        SolidStiffness(TacsScalar rho, TacsScalar E, TacsScalar nu,
                       TacsScalar ys, int eNum)
cdef extern from "CoupledThermoSolidStiffness.h":
    cdef cppclass CoupledThermoSolidStiffness(SolidStiffness):
        CoupledThermoSolidStiffness()
        CoupledThermoSolidStiffness( TacsScalar, TacsScalar, TacsScalar,
                                     TacsScalar, TacsScalar, TacsScalar )

cdef extern from "TACSConstitutiveWrapper.h":
    cdef cppclass PSStiffnessWrapper(PlaneStressStiffness):
        PSStiffnessWrapper(PyObject *obj)

        # Member functions
        void (*setdesignvars)(void*, const TacsScalar*, int)
        void (*getdesignvars)(void*, TacsScalar*, int)
        void (*getdesignvarrange)(void*, TacsScalar*, TacsScalar*, int)
        void (*calculatestress)(void*, const double*, 
                                const TacsScalar*, TacsScalar*)
        void (*addstressdvsens)(void*, const double*, const TacsScalar*,
                                TacsScalar, const TacsScalar*,
                                TacsScalar*, int)
        void (*getpointwisemass)(void*, const double*, TacsScalar*)
        void (*addpointwisemassdvsens)(void*, const double*,
                                       const TacsScalar*, TacsScalar*, int )
        TacsScalar (*fail)(void*, const double*, const TacsScalar*)
        void (*failstrainsens)(void*, const double*,
                               const TacsScalar*, TacsScalar*)
        void (*addfaildvsens)(void*, const double*, const TacsScalar*, 
                              TacsScalar, TacsScalar*, int)

    cdef cppclass FSDTStiffnessWrapper(FSDTStiffness):
        FSDTStiffnessWrapper(PyObject *obj)

        # Member functions
        void (*setdesignvars)(void*, const TacsScalar*, int)
        void (*getdesignvars)(void*, TacsScalar*, int)
        void (*getdesignvarrange)(void*, TacsScalar*, TacsScalar*, int)
        TacsScalar (*getstiffness)(void*, const double*, 
                                   TacsScalar*, TacsScalar*,
                                   TacsScalar*, TacsScalar*)
        void (*addstiffnessdvsens)(void*, const double*, const TacsScalar*,
                                   const TacsScalar*, TacsScalar,
                                   TacsScalar*, int)
        void (*getpointwisemass)(void*, const double*, TacsScalar*)
        void (*addpointwisemassdvsens)(void*, const double*,
                                       const TacsScalar*, TacsScalar*, int )
        TacsScalar (*fail)(void*, const double*, const TacsScalar*)
        void (*failstrainsens)(void*, const double*,
                               const TacsScalar*, TacsScalar*)
        void (*addfaildvsens)(void*, const double*, const TacsScalar*, 
                              TacsScalar, TacsScalar*, int)

cdef class Timoshenko(Constitutive):
    pass

cdef class FSDT(Constitutive):
    pass

cdef class PlaneStress(Constitutive):
    pass
   
cdef class SolidStiff(Constitutive):
    pass

cdef class CoupledPlaneStress(Constitutive):
    pass

cdef class CoupledSolid(Constitutive):
    pass

# Special functions required for converting pointers
cdef extern from "":
    PlaneStressStiffness* _dynamicPlaneStress"dynamic_cast<PlaneStressStiffness*>"(TACSConstitutive*)
    FSDTStiffness* _dynamicFSDT"dynamic_cast<FSDTStiffness*>"(TACSConstitutive*)
    SolidStiffness* _dynamicSolid"dynamic_cast<SolidStiffness*>"(TACSConstitutive*)
    TimoshenkoStiffness* _dynamicTimoshenko"dynamic_cast<TimoshenkoStiffness*>"(TACSConstitutive*)
    CoupledThermoPlaneStressStiffness* _dynamicPSThermo"dynamic_cast<CoupledThermoPlaneStressStiffness*>"(TACSConstitutive*)
    CoupledThermoSolidStiffness* _dynamicSolidThermo"dynamic_cast<CoupledThermoSolidStiffness*>"(TACSConstitutive*)
