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
import numpy as np
cimport numpy as np

# Ensure that numpy is initialized
np.import_array()

# Import the definition required for const strings
from libc.string cimport const_char
from libc.stdlib cimport malloc, free

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

cdef class MaterialProperties:
    def __cinit__(self, **kwargs):
        # Set the density and specific heat value
        cdef TacsScalar rho = 2700.0
        cdef TacsScalar specific_heat = 921.0

        # Stiffness properties
        cdef TacsScalar E = 70e9
        cdef TacsScalar E1 = 70e9
        cdef TacsScalar E2 = 70e9
        cdef TacsScalar E3 = 70e9

        cdef TacsScalar nu = 0.3
        cdef TacsScalar nu12 = 0.3
        cdef TacsScalar nu13 = 0.3
        cdef TacsScalar nu23 = 0.3

        cdef TacsScalar G12 = 0.5*70e9/(1.0 + 0.3)
        cdef TacsScalar G13 = 0.5*70e9/(1.0 + 0.3)
        cdef TacsScalar G23 = 0.5*70e9/(1.0 + 0.3)

        # Strength properties
        cdef TacsScalar ys = 270e6
        cdef TacsScalar T1 = 0.0
        cdef TacsScalar C1 = 0.0
        cdef TacsScalar T2 = 0.0
        cdef TacsScalar C2 = 0.0
        cdef TacsScalar T3 = 0.0
        cdef TacsScalar C3 = 0.0
        cdef TacsScalar S12 = 0.0
        cdef TacsScalar S13 = 0.0
        cdef TacsScalar S23 = 0.0

        # Set the thermal coefficient of expansion
        cdef TacsScalar alpha = 24e-6
        cdef TacsScalar alpha1 = 24e-6
        cdef TacsScalar alpha2 = 24e-6
        cdef TacsScalar alpha3 = 24e-6

        # Set the thermal conductivity
        cdef TacsScalar kappa = 230.0
        cdef TacsScalar kappa1 = 230.0
        cdef TacsScalar kappa2 = 230.0
        cdef TacsScalar kappa3 = 230.0

        if 'rho' in kwargs:
            rho = kwargs['rho']
        if 'specific_heat' in kwargs:
            specific_heat = kwargs['specific_heat']

        # Set the isotropic properties if they are defined
        if 'E' in kwargs:
            E = kwargs['E']
            E1 = E
            E2 = E
            E3 = E
        if 'nu' in kwargs:
            nu = kwargs['nu']
            nu12 = nu
            nu23 = nu
            nu13 = nu
        if 'ys' in kwargs:
            ys = kwargs['ys']
        if 'alpha' in kwargs:
            alpha = kwargs['alpha']
            alpha1 = alpha
            alpha2 = alpha
            alpha3 = alpha
        if 'kappa' in kwargs:
            kappa = kwargs['kappa']
            kappa1 = kappa
            kappa2 = kappa
            kappa3 = kappa

        # Stiffness properties
        if 'E1' in kwargs:
            E1 = kwargs['E1']
        if 'E2' in kwargs:
            E2 = kwargs['E2']
        if 'E3' in kwargs:
            E3 = kwargs['E3']
        if 'nu12' in kwargs:
            nu12 = kwargs['nu12']
        if 'nu13' in kwargs:
            nu13 = kwargs['nu13']
        if 'nu23' in kwargs:
            nu23 = kwargs['nu23']
        if 'G12' in kwargs:
            G12 = kwargs['G12']
        if 'G13' in kwargs:
            G13 = kwargs['G13']
        if 'G23' in kwargs:
            G23 = kwargs['G23']

        # Strength properties
        if 'T1' in kwargs:
            T1 = kwargs['T1']
        if 'Xt' in kwargs:
            T1 = kwargs['Xt']
        if 'C1' in kwargs:
            C1 = kwargs['C1']
        if 'Xc' in kwargs:
            C1 = kwargs['Xc']

        if 'T2' in kwargs:
            T2 = kwargs['T2']
        if 'Yt' in kwargs:
            T2 = kwargs['Yt']
        if 'C2' in kwargs:
            C2 = kwargs['C2']
        if 'Yc' in kwargs:
            C2 = kwargs['Yc']

        if 'T3' in kwargs:
            T3 = kwargs['T3']
        if 'C3' in kwargs:
            C3 = kwargs['C3']

        if 'S12' in kwargs:
            S12 = kwargs['S12']
        if 'S' in kwargs:
            S12 = kwargs['S']
        if 'S13' in kwargs:
            S13 = kwargs['S13']
        if 'S23' in kwargs:
            S23 = kwargs['S23']

        # Set the thermal coefficient of expansion
        if 'alpha1' in kwargs:
            alpha1 = kwargs['alpha1']
        if 'alpha2' in kwargs:
            alpha2 = kwargs['alpha2']
        if 'alpha3' in kwargs:
            alpha3 = kwargs['alpha3']

        # Set the thermal conductivity
        if 'kappa1' in kwargs:
            kappa1 = kwargs['kappa1']
        if 'kappa2' in kwargs:
            kappa2 = kwargs['kappa2']
        if 'kappa3' in kwargs:
            kappa3 = kwargs['kappa3']

        if 'E1' in kwargs:
            self.ptr = new TACSMaterialProperties(rho, specific_heat, E1, E2, E3,
                                                  nu12, nu13, nu23, G12, G13, G23,
                                                  T1, C1, T2, C2, T3, C3, S12, S13, S23,
                                                  alpha1, alpha2, alpha3,
                                                  kappa1, kappa2, kappa3)
        else:
            self.ptr = new TACSMaterialProperties(rho, specific_heat,
                                                  E, nu, ys, alpha, kappa)
        self.ptr.incref()

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getMaterialProperties(self):
        """
        Return a dictionary of the material properties
        """

        cdef TacsScalar E1, E2, E3, nu12, nu13, nu23, G12, G13, G23
        cdef TacsScalar T1, C1, T2, C2, T3, C3, S12, S13, S23
        cdef TacsScalar alpha1, alpha2, alpha3
        cdef TacsScalar kappa1, kappa2, kappa3

        self.ptr.getOrthotropicProperties(&E1, &E2, &E3, &nu12, &nu13, &nu23, &G12, &G13, &G23)
        self.ptr.getStrengthProperties(&T1, &C1, &T2, &C2, &T3, &C3,
                                       &S12, &S13, &S23)
        self.ptr.getCoefThermalExpansion(&alpha1, &alpha2, &alpha3)
        self.ptr.getThermalConductivity(&kappa1, &kappa2, &kappa3)

        mat = {}
        mat['E1'] = E1
        mat['E2'] = E2
        mat['E3'] = E3
        mat['nu12'] = nu12
        mat['nu13'] = nu13
        mat['nu23'] = nu23
        mat['G12'] = G12
        mat['G13'] = G13
        mat['G23'] = G23

        mat['T1'] = T1
        mat['C1'] = C1
        mat['T2'] = T2
        mat['C2'] = C2
        mat['T3'] = T3
        mat['C3'] = C3
        mat['S12'] = S12
        mat['S13'] = S13
        mat['S23'] = S23

        mat['alpha1'] = alpha1
        mat['alpha2'] = alpha2
        mat['alpha3'] = alpha3

        mat['kappa1'] = kappa1
        mat['kappa2'] = kappa2
        mat['kappa3'] = kappa3

        return mat

    def setDensity(self, TacsScalar rho):
        self.ptr.setDensity(rho)

    def setSpecificHeat(self, TacsScalar specific_heat):
        self.ptr.setSpecificHeat(specific_heat)

cdef class OrthotropicPly:
    cdef TACSOrthotropicPly *ptr
    def __cinit__(self, TacsScalar ply_thickness, MaterialProperties props,
                  max_strain_criterion=False):
        self.ptr = new TACSOrthotropicPly(ply_thickness, props.ptr)
        self.ptr.incref()
        if max_strain_criterion:
            self.ptr.setUseMaxStrainCriterion()
        else:
            self.ptr.setUseTsaiWuCriterion()

    def __dealloc__(self):
        self.ptr.decref()

cdef class PlaneStressConstitutive(Constitutive):
    def __cinit__(self, *args, **kwargs):
        cdef TACSMaterialProperties *props = NULL
        cdef TacsScalar t = 1.0
        cdef int tNum = -1
        cdef TacsScalar tlb = 0.0
        cdef TacsScalar tub = 10.0

        if len(args) >= 1:
            props = (<MaterialProperties>args[0]).ptr
        if 't' in kwargs:
            t = kwargs['t']
        if 'tNum' in kwargs:
            tNum = kwargs['tNum']
        if 'tlb' in kwargs:
            tlb = kwargs['tlb']
        if 'tub' in kwargs:
            tub = kwargs['tub']

        if props is not NULL:
            self.cptr = new TACSPlaneStressConstitutive(props, t, tNum,
                                                        tlb, tub)
            self.ptr = self.cptr
            self.ptr.incref()
        else:
            self.ptr = NULL
            self.cptr = NULL

cdef class SolidConstitutive(Constitutive):
    def __cinit__(self, *args, **kwargs):
        cdef TACSMaterialProperties *props = NULL
        cdef TacsScalar t = 1.0
        cdef int tNum = -1
        cdef TacsScalar tlb = 0.0
        cdef TacsScalar tub = 10.0

        if len(args) >= 1:
            props = (<MaterialProperties>args[0]).ptr
        if 't' in kwargs:
            t = kwargs['t']
        if 'tNum' in kwargs:
            tNum = kwargs['tNum']
        if 'tlb' in kwargs:
            tlb = kwargs['tlb']
        if 'tub' in kwargs:
            tub = kwargs['tub']

        if props is not NULL:
            self.cptr = new TACSSolidConstitutive(props, t, tNum,
                                                  tlb, tub)
            self.ptr = self.cptr
            self.ptr.incref()
        else:
            self.ptr = NULL
            self.cptr = NULL

    def getMaterialProperties(self):
        if self.cptr:
            return _init_MaterialProperties(self.cptr.getMaterialProperties())
        return None

cdef class IsoShellConstitutive(ShellConstitutive):
    def __cinit__(self, *args, **kwargs):
        cdef TACSMaterialProperties *props = NULL
        cdef TacsScalar t = 1.0
        cdef int tNum = -1
        cdef TacsScalar tlb = 0.0
        cdef TacsScalar tub = 10.0

        if len(args) >= 1:
            props = (<MaterialProperties>args[0]).ptr
        if 't' in kwargs:
            t = kwargs['t']
        if 'tNum' in kwargs:
            tNum = kwargs['tNum']
        if 'tlb' in kwargs:
            tlb = kwargs['tlb']
        if 'tub' in kwargs:
            tub = kwargs['tub']

        if props is not NULL:
            self.cptr = new TACSIsoShellConstitutive(props, t, tNum,
                                                     tlb, tub)
            self.ptr = self.cptr
            self.ptr.incref()
        else:
            self.ptr = NULL
            self.cptr = NULL

cdef class CompositeShellConstitutive(ShellConstitutive):
    def __cinit__(self, ply_list,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] ply_thicknesses,
                  np.ndarray[TacsScalar, ndim=1, mode='c'] ply_angles,
                  TacsScalar kcorr=5.0/6.0):

        num_plies = len(ply_list)
        if len(ply_thicknesses) != num_plies:
            raise ValueError('Ply thickness array must match length of ply list')
        if len(ply_angles) != num_plies:
            raise ValueError('Ply angle array must match length of ply list')

        # Allocate the array of TACSOrthotropicPly pointers
        cdef TACSOrthotropicPly **plys
        plys = <TACSOrthotropicPly**>malloc(num_plies*sizeof(TACSOrthotropicPly*))
        if plys is NULL:
            raise MemoryError()

        for i in range(num_plies):
            plys[i] = (<OrthotropicPly>ply_list[i]).ptr

        self.cptr = new TACSCompositeShellConstitutive(num_plies, plys,
                                                       <TacsScalar*>ply_thicknesses.data,
                                                       <TacsScalar*>ply_angles.data, kcorr)
        self.ptr = self.cptr
        self.ptr.incref()

        # Free the allocated array
        free(plys)

cdef class LamParamShellConstitutive(ShellConstitutive):
    def __cinit__(self, OrthotropicPly ply, **kwargs):
        cdef TacsScalar t = 1.0
        cdef int t_num = -1
        cdef TacsScalar min_t = 1.0
        cdef TacsScalar max_t = 1.0
        cdef TacsScalar f0 = 0.25
        cdef TacsScalar f45 = 0.5
        cdef TacsScalar f90 = 0.25
        cdef int f0_num = -1
        cdef int f45_num = -1
        cdef int f90_num = -1
        cdef TacsScalar min_f0 = 0.125
        cdef TacsScalar min_f45 = 0.125
        cdef TacsScalar min_f90 = 0.125
        cdef TacsScalar W1 = 0.0
        cdef TacsScalar W3 = 0.0
        cdef int W1_num = -1
        cdef int W3_num = -3
        cdef TacsScalar ksWeight = 30.0
        cdef TacsScalar epsilon = 0.0

        if 't' in kwargs:
            t = kwargs['t']
        if 't_num' in kwargs:
            t_num = kwargs['t_num']
        if 'min_t' in kwargs:
            min_t = kwargs['min_t']
        if 'max_t' in kwargs:
            max_t = kwargs['max_t']

        if 'f0' in kwargs:
            f0 = kwargs['f0']
        if 'f45' in kwargs:
            f45 = kwargs['f45']
        if 'f90' in kwargs:
            f90 = kwargs['f90']

        if 'f0_num' in kwargs:
            f0_num = kwargs['f0_num']
        if 'f45_num' in kwargs:
            f45_num = kwargs['f45_num']
        if 'f90_num' in kwargs:
            f90_num = kwargs['f90_num']

        if 'min_f0' in kwargs:
            min_f0 = kwargs['min_f0']
        if 'min_f45' in kwargs:
            min_f45 = kwargs['min_f45']
        if 'min_f90' in kwargs:
            min_f90 = kwargs['min_f90']

        if 'W1' in kwargs:
            W1 = kwargs['W1']
        if 'W3' in kwargs:
            W3 = kwargs['W3']
        if 'W1_num' in kwargs:
            W1_num = kwargs['W1_num']
        if 'W3_num' in kwargs:
            W3_num = kwargs['W3_num']
        if 'ksWeight' in kwargs:
            ksWeight = kwargs['ksWeight']
        if 'epsilon' in kwargs:
            epsilon = kwargs['epsilon']

        self.cptr = new TACSLamParamShellConstitutive(ply.ptr, t, t_num, min_t, max_t,
                                                      f0, f45, f90, f0_num, f45_num, f90_num,
                                                      min_f0, min_f45, min_f90,
                                                      W1, W3, W1_num, W3_num, ksWeight, epsilon)
        self.ptr = self.cptr
        self.ptr.incref()

cdef class TimoshenkoConstitutive(Constitutive):
    def __cinit__(self, rhoA, rhoIy, rhoIz, rhoIyz,
                      EA, GJ, EIy, EIz, kGAy, kGAz,
                      np.ndarray[TacsScalar, ndim=1, mode='c'] axis):
        self.cptr = new TACSTimoshenkoConstitutive(rhoA, rhoIy, rhoIz, rhoIyz,
                                                  EA, GJ, EIy, EIz, kGAy, kGAz,
                                                  <TacsScalar*>axis.data)
        self.ptr = self.cptr
        self.ptr.incref()

def TestConstitutive(Constitutive con, int elemIndex=0, double dh=1e-6,
                     int test_print_level=2, double atol=1e-30,
                     double rtol=1e-5):
    return TacsTestConstitutive(con.ptr, elemIndex, dh,
                                test_print_level, atol, rtol)

def TestConstitutiveDensity(Constitutive con, int elem_index,
                            np.ndarray[double, ndim=1, mode='c'] pt,
                            np.ndarray[TacsScalar, ndim=1, mode='c'] x,
                            np.ndarray[TacsScalar, ndim=1, mode='c'] dvs,
                            double dh=1e-6, int test_print_level=2, double atol=1e-30,
                            double rtol=1e-5):
    ndvs = len(dvs)
    assert len(pt) == 3
    assert len(x) == 3
    return TacsTestConstitutiveDensity(con.ptr, elem_index, <double*>pt.data, <TacsScalar*>x.data, ndvs,
                                       <TacsScalar*>dvs.data, dh, test_print_level, atol, rtol)

def TestConstitutiveSpecificHeat(Constitutive con, int elem_index,
                                 np.ndarray[double, ndim=1, mode='c'] pt,
                                 np.ndarray[TacsScalar, ndim=1, mode='c'] x,
                                 np.ndarray[TacsScalar, ndim=1, mode='c'] dvs,
                                 double dh=1e-6, int test_print_level=2, double atol=1e-30,
                                 double rtol=1e-5):
    ndvs = len(dvs)
    assert len(pt) == 3
    assert len(x) == 3
    return TacsTestConstitutiveSpecificHeat(con.ptr, elem_index, <double*>pt.data, <TacsScalar*>x.data, ndvs,
                                            <TacsScalar*>dvs.data, dh, test_print_level, atol, rtol)

def TestConstitutiveHeatFlux(Constitutive con, int elem_index,
                             np.ndarray[double, ndim=1, mode='c'] pt,
                             np.ndarray[TacsScalar, ndim=1, mode='c'] x,
                             np.ndarray[TacsScalar, ndim=1, mode='c'] dvs,
                             double dh=1e-6, int test_print_level=2, double atol=1e-30,
                             double rtol=1e-5):
    ndvs = len(dvs)
    assert len(pt) == 3
    assert len(x) == 3
    return TacsTestConstitutiveHeatFlux(con.ptr, elem_index, <double*>pt.data, <TacsScalar*>x.data, ndvs,
                                        <TacsScalar*>dvs.data, dh, test_print_level, atol, rtol)

def TestConstitutiveStress(Constitutive con, int elem_index,
                           np.ndarray[double, ndim=1, mode='c'] pt,
                           np.ndarray[TacsScalar, ndim=1, mode='c'] x,
                           np.ndarray[TacsScalar, ndim=1, mode='c'] dvs,
                           double dh=1e-6, int test_print_level=2, double atol=1e-30,
                           double rtol=1e-5):
    ndvs = len(dvs)
    assert len(pt) == 3
    assert len(x) == 3
    return TacsTestConstitutiveStress(con.ptr, elem_index, <double*>pt.data, <TacsScalar*>x.data, ndvs,
                                      <TacsScalar*>dvs.data, dh, test_print_level, atol, rtol)

def TestConstitutiveThermalStrain(Constitutive con, int elem_index,
                                  np.ndarray[double, ndim=1, mode='c'] pt,
                                  np.ndarray[TacsScalar, ndim=1, mode='c'] x,
                                  np.ndarray[TacsScalar, ndim=1, mode='c'] dvs,
                                  double dh=1e-6, int test_print_level=2, double atol=1e-30,
                                  double rtol=1e-5):
    ndvs = len(dvs)
    assert len(pt) == 3
    assert len(x) == 3
    return TacsTestConstitutiveThermalStrain(con.ptr, elem_index, <double*>pt.data, <TacsScalar*>x.data, ndvs,
                                             <TacsScalar*>dvs.data, dh, test_print_level, atol, rtol)

def TestConstitutiveFailure(Constitutive con, int elem_index,
                            np.ndarray[double, ndim=1, mode='c'] pt,
                            np.ndarray[TacsScalar, ndim=1, mode='c'] x,
                            np.ndarray[TacsScalar, ndim=1, mode='c'] dvs,
                            double dh=1e-6, int test_print_level=2, double atol=1e-30,
                            double rtol=1e-5):
    ndvs = len(dvs)
    assert len(pt) == 3
    assert len(x) == 3
    return TacsTestConstitutiveFailure(con.ptr, elem_index, <double*>pt.data, <TacsScalar*>x.data, ndvs,
                                             <TacsScalar*>dvs.data, dh, test_print_level, atol, rtol)

def TestConstitutiveFailureStrainSens(Constitutive con, int elem_index,
                                      np.ndarray[double, ndim=1, mode='c'] pt,
                                      np.ndarray[TacsScalar, ndim=1, mode='c'] x,
                                      double dh=1e-6, int test_print_level=2, double atol=1e-30,
                                      double rtol=1e-5):
    assert len(pt) == 3
    assert len(x) == 3
    return TacsTestConstitutiveFailureStrainSens(con.ptr, elem_index, <double*>pt.data, <TacsScalar*>x.data,
                                                 dh, test_print_level, atol, rtol)