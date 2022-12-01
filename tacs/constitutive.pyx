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
    """
    This class stores the mechanical and thermal material properties for
    isotropic and anisotropic materials.

    The goal of this class is to store a set of material properties
    that can be queried by constitutive classes for beams, shells,
    plane stress and solid elements. The minimum set of properties
    consists of isotroic mechanical properties, with zero thermal
    expansion.

    There are several different ways of specifying material
    properties. The following describes several of the possible
    ways and the appropriate situations:

        'rho' + 'specific_heat' + 'kappa':

            Specifies the density, specific heat, and thermal conductivity for an isotropic material.
            This is appropriate for 2D or 3D heat conduction analysis.

        'rho' + 'specific_heat' + 'kappa1' + 'kappa2':

            Specifies the density, specific heat, and thermal conductivity for an orthotropic material.
            This is appropriate for 2D heat conduction analysis.

        'rho' + 'specific_heat' + 'kappa1' + 'kappa2' + 'kappa3':

            Specifies the density, specific heat, and thermal conductivity for an orthotropic material.
            This is appropriate for 3D heat conduction analysis.

        'rho' + 'E' + 'nu' + 'ys':

            Specifies the density, Youngs' modulus, Poisson's ratio, and yield strength for an isotropic material.
            This is appropriate for 2D or 3D elastic structural analysis.

        'rho' + 'E1' + 'E2' + 'nu12' + 'G12' + 'T1' + 'T2' + 'C1' + 'C2' + 'S12':

            Specifies the density, Youngs' moduli, Poisson's ratios, and strength values for an orthotropic material.
            This is appropriate for 2D elastic structural analysis.

        'rho' + 'E1' + 'E2' + 'E3' + 'nu12' + 'nu13' + 'nu23' + 'G12' + 'G13' + 'G23'
        + 'T1' + 'T2' + 'T3' + 'C1' + 'C2' + 'C3' + 'S12' + 'S13' + 'S23':

            Specifies the density, Youngs' moduli, Poisson's ratios, and strength values for an orthotropic material.
            This is appropriate for 3D elastic structural analysis.

        'rho'  + 'specific_heat' + 'kappa' + 'alpha' + 'E' + 'nu' + 'ys':

            Specifies the density, specific heat, thermal conductivity, thermal expansion, Youngs' modulus,
            Poisson's ratio, and yield strength for an isotropic material.
            This is appropriate for 2D or 3D thermoelastic structural analysis.

        'rho'  + 'specific_heat' + 'kappa1' + 'kappa2' + 'alpha1' + 'alpha2'
        + 'E1' + 'E2' + 'nu12' + 'G12' + 'T1' + 'T2' + 'C1' + 'C2' + 'S12':

            Specifies the density, specific heat, thermal conductivity, thermal expansion, Youngs' moduli,
            Poisson's ratios, and strength values for an orthotropic material.
            This is appropriate for 2D thermoelastic structural analysis.

        'rho'  + 'specific_heat' + 'kappa1' + 'kappa2' + 'kappa3' + 'alpha1' + 'alpha2' + 'alpha3'
        + 'E1' + 'E2' + 'E3' + 'nu12' + 'nu13' + 'nu23' + 'G12' + 'G13' + 'G23'
        + 'T1' + 'T2' + 'T3' + 'C1' + 'C2' + 'C3' + 'S12' + 'S13' + 'S23':

            Specifies the density, specific heat, thermal conductivity, thermal expansion, Youngs' moduli,
            Poisson's ratios, and strength values for an orthotropic material.
            This is appropriate for 3D thermoelastic structural analysis.

    All parameters are optional and specified using a keyword arg format.

    Args:
        rho (float or complex, optional): The material density (keyword argument).

        specific_heat (float or complex, optional): The material specific heat (keyword argument).

        kappa (float or complex, optional): The material isotropic thermal conductivity (keyword argument).
        kappa1 (float or complex, optional): The material orthotropic thermal conductivity
            in the '1' direction (keyword argument).
        kappa2 (float or complex, optional): The material orthotropic thermal conductivity
           in the '2' direction (keyword argument).
        kappa3 (float or complex, optional): The material orthotropic thermal conductivity
           in the '3' direction (keyword argument).

        alpha (float or complex, optional): The material isotropic thermal expansion coeeficient (keyword argument).
        alpha1 (float or complex): The material orthotropic thermal expansion coefficient in the '1' direction
            (keyword argument).
        alpha2 (float or complex, optional): The material orthotropic thermal expansion coefficient in the '2' direction
            (keyword argument).
        alpha3 (float or complex, optional): The material orthotropic thermal expansion coefficient in the '3' direction
            (keyword argument).

        E (float or complex, optional): The material isotropic Youngs' modulus (keyword argument).
        E1 (float or complex, optional): The material orthotropic Youngs' modulus
            in the '1' direction (keyword argument).
        E2 (float or complex, optional): The material orthotropic Youngs' modulus
            in the '2' direction (keyword argument).
        E3 (float or complex, optional): The material orthotropic Youngs' modulus
            in the '3' direction (keyword argument).

        nu (float or complex, optional): The material isotropic Poisson's ratio (keyword argument).
        nu12 (float or complex, optional): The material orthotropic Poisson's ratio
            in the '1-2' plane (keyword argument).
        nu13 (float or complex, optional): The material orthotropic Poisson's ratio
            in the '1-3' plane (keyword argument).
        nu23 (float or complex, optional): The material orthotropic Poisson's ratio
            in the '2-3' plane (keyword argument).

        G (float or complex, optional): The material isotropic shear modulus (keyword argument).
        G12 (float or complex, optional): The material orthotropic shear modulus in the '1-2' plane (keyword argument).
        G13 (float or complex, optional): The material orthotropic shear modulus in the '1-3' plane (keyword argument).
        G23 (float or complex, optional): The material orthotropic shear modulus in the '2-3' plane (keyword argument).

        ys (float or complex, optional): The material isotropic yield strength (keyword argument).
        T1 (float or complex, optional): The material orthotropic tension strength
            in the '1' direction (keyword argument).
        T2 (float or complex, optional): The material orthotropic tension strength
            in the '2' direction (keyword argument).
        T3 (float or complex, optional): The material orthotropic tension strength
            in the '3' direction (keyword argument).
        C1 (float or complex, optional): The material orthotropic compression strength
            in the '1' direction (keyword argument).
        C2 (float or complex, optional): The material orthotropic compression strength
            in the '2' direction (keyword argument).
        C3 (float or complex, optional): The material orthotropic compression strength
            in the '3' direction (keyword argument).
        S12 (float or complex, optional): The material orthotropic shear strength in the '1-2' plane (keyword argument).
        S13 (float or complex, optional): The material orthotropic shear strength in the '1-3' plane (keyword argument).
        S23 (float or complex, optional): The material orthotropic shear strength in the '2-3' plane (keyword argument).
    """
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

        Returns:
            mat (dict): Dictionary holding material property information
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
        """
        Set the density property values

        Args:
            rho (float or complex): The material density
        """
        self.ptr.setDensity(rho)

    def setSpecificHeat(self, TacsScalar specific_heat):
        """
        Set the density property values

        Args:
            rho (float or complex): The material specific heat
        """
        self.ptr.setSpecificHeat(specific_heat)

cdef class OrthotropicPly:
    """
      The following class holds the material stiffness and strength
      properties for an orthotropic ply. This class is used by several
      constitutive classes within TACS.

      The interaction coefficient for the Tsai-Wu failure criterion is set
      to zero by default. If a value of C, the failure stress under
      combined in-plane loading, is supplied, the interaction coefficient
      is determined. Be careful - the value can easily fall outside
      acceptible bounds - these are tested during initialization.

      Args:
          ply_thickness (float or complex): The ply thickness.
          props (MaterialProperties): The ply material property.
          max_strain_criterion (bool): Flag to determine if max strain strength criterion is to be used.
            Defaults to False (i.e. use max strength).
    """
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
    """
    This is the base class for the plane stress constitutive objects.
    All objects performing plane stress analysis should utilize this class.

    Args:
        props (MaterialProperties): The material property.
        t (float or complex, optional): The material thickness (keyword argument). Defaults to 1.0.
        tNum (int, optional): Design variable number to assign to thickness (keyword argument). Defaults to -1
            (i.e. no design variable).
        tlb (float or complex, optional): Thickness varaible lower bound (keyword argument). Defaults to 0.0.
        tub (float or complex, optional): Thickness varaible upper bound (keyword argument). Defaults to 10.0.
    """
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

cdef class PhaseChangeMaterialConstitutive(Constitutive):
    """
    This is the base class for the phase change material constitutive objects.

    Args:
        solid_props (MaterialProperties): The material property of the solid phase.
        liquid_props (MaterialProperties): The material property of the liquid phase.
        lh (float or complex): The specific latent heat of the material.
        mt (float or complex): The melting temperature of the material.
        t (float or complex, optional): The material thickness (keyword argument). Defaults to 1.0.
        tNum (int, optional): Design variable number to assign to thickness (keyword argument). Defaults to -1
            (i.e. no design variable).
        tlb (float or complex, optional): Thickness varaible lower bound (keyword argument). Defaults to 0.0.
        tub (float or complex, optional): Thickness varaible upper bound (keyword argument). Defaults to 10.0.
    """
    def __cinit__(self, *args, **kwargs):
        cdef TACSMaterialProperties *solid_props = NULL
        cdef TACSMaterialProperties *liquid_props = NULL
        cdef TacsScalar lh = 0.0
        cdef TacsScalar Tm = 0.0
        cdef TacsScalar dT = 10.0
        cdef TacsScalar t = 1.0
        cdef int tNum = -1
        cdef TacsScalar tlb = 0.0
        cdef TacsScalar tub = 10.0

        if len(args) >= 2:
            solid_props = (<MaterialProperties>args[0]).ptr
            liquid_props = (<MaterialProperties>args[1]).ptr
        if 'lh' in kwargs:
            lh = kwargs['lh']
        if 'Tm' in kwargs:
            Tm = kwargs['Tm']
        if 'dT' in kwargs:
            dT = kwargs['dT']
        if 't' in kwargs:
            t = kwargs['t']
        if 'tNum' in kwargs:
            tNum = kwargs['tNum']
        if 'tlb' in kwargs:
            tlb = kwargs['tlb']
        if 'tub' in kwargs:
            tub = kwargs['tub']

        if solid_props is not NULL:
            self.cptr = new TACSPhaseChangeMaterialConstitutive(solid_props,
                                                                liquid_props,
                                                                lh, Tm, dT, t,
                                                                tNum, tlb, tub)
            self.ptr = self.cptr
            self.ptr.incref()
        else:
            self.ptr = NULL
            self.cptr = NULL

cdef class SolidConstitutive(Constitutive):
    """
    This is the base class for the solid constitutive objects.
    All objects performing solid elastic analysis should utilize this class.

    Args:
        props (MaterialProperties): The material property.
        t (float or complex, optional): The material "thickness" (keyword argument).
            Weighting factor used for topology optimization. 0.0 corresponds to void,
            1.0 corresponds material fully present, values in between are intermediate. Defaults to 1.0.
        tNum (int, optional): Design variable number to assign to "thickness" (keyword argument). Defaults to -1
            (i.e. no design variable).
        tlb (float or complex, optional): "Thickness" varaible lower bound (keyword argument). Defaults to 0.0.
        tub (float or complex, optional): "Thickness" varaible upper bound (keyword argument). Defaults to 10.0.
    """
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

cdef class ShellConstitutive(Constitutive):
    """
    This is the base class for the shell constitutive objects.
    All objects performing shell elastic analysis should utilize this class.
    """

    def setDrillingRegularization(self, double kpenalty=10.0):
        """
        Update regularization parameter used to stiffen shell in drilling rotation dof.

        Args:
            kpenalty (float): Drilling regularization parameter. Defaults to 10.0.
        """
        if self.cptr:
            self.cptr.setDrillingRegularization(kpenalty)

cdef class IsoShellConstitutive(ShellConstitutive):
    """
    This constitutive class defines the stiffness properties for a
    isotropic first-order shear deformation theory shell type element.

    Args:
        props (MaterialProperties): The material property.
        t (float or complex, optional): The material thickness (keyword argument). Defaults to 1.0.
        tNum (int, optional): Design variable number to assign to thickness (keyword argument). Defaults to -1
            (i.e. no design variable).
        tlb (float or complex, optional): Thickness variable lower bound (keyword argument). Defaults to 0.0.
        tub (float or complex, optional): Thickness variable upper bound (keyword argument). Defaults to 10.0.
    """
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
    """
    This constitutive class defines the stiffness properties for a
    composite laminate first-order shear deformation theory (FSDT) shell type element.

    Args:
       ply_list (list[OrthotropicPly]): List of ply properties in layup.
       ply_thicknesses (numpy.ndarray[float or complex]): Array of ply thicknesses in layup.
       ply_angles (numpy.ndarray[float or complex]): Array of ply angles (in radians) in layup.
       kcorr (float or complex, optional): FSDT shear correction factor. Defaults to 5.0/6.0.
    """
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

cdef class BasicBeamConstitutive(BeamConstitutive):
    """
    Timoshenko theory based constitutive object for an general beam.

     .. note::
        TACS uses a negative sign convention in the product of inertia definition, for example:
            Iyz = -int[y * z * dA]

        The moments of area are always positive, as usual:
            Iyy = int[(z^2 * dA]

    Args:
        props (MaterialProperties): The material property.
        A (float or complex): Beam cross-sectional area (keyword argument). Defaults to 0.0.
        J (float or complex): Beam polar moment of area about x-axis (keyword argument). Defaults to 0.0.
        Iy (float or complex): Beam area moment of area about y-axis (keyword argument). Defaults to 0.0.
        Iz (float or complex): Beam area moment of area about z-axis (keyword argument). Defaults to 0.0.
        Iyz (float or complex): Beam product of area in yz-plane (keyword argument). Defaults to 0.0.
        ky (float or complex): Shear correction factor in y-direction (keyword argument). Defaults to 5/6.
        kz (float or complex): Shear correction factor in z-direction (keyword argument). Defaults to 5/6.
    """
    def __cinit__(self, *args, **kwargs):
        cdef TACSMaterialProperties *props = NULL
        cdef TacsScalar A = 0.0
        cdef TacsScalar J = 0.0
        cdef TacsScalar Iy = 0.0
        cdef TacsScalar Iz = 0.0
        cdef TacsScalar Iyz = 0.0
        cdef TacsScalar ky = 5.0/6.0
        cdef TacsScalar kz = 5.0/6.0

        if len(args) >= 1:
            props = (<MaterialProperties>args[0]).ptr
        if 'A' in kwargs:
            A = kwargs['A']
        if 'J' in kwargs:
            J = kwargs['J']
        if 'Iy' in kwargs:
            Iy = kwargs['Iy']
        if 'Iz' in kwargs:
            Iz = kwargs['Iz']
        if 'Iyz' in kwargs:
            Iyz = kwargs['Iyz']
        if 'ky' in kwargs:
            ky = kwargs['ky']
        if 'kz' in kwargs:
            kz = kwargs['kz']

        self.cptr = new TACSBasicBeamConstitutive(props, A, J, Iy, Iz, Iyz, ky, kz)
        self.ptr = self.cptr
        self.ptr.incref()

cdef class IsoTubeBeamConstitutive(BeamConstitutive):
    """
    Timoshenko theory based constitutive object for a hollow circular beam.

    Args:
        props (MaterialProperties): The material property.
        d (float or complex, optional): Tube inner diameter (keyword argument). Defaults to 1.0.
        dNum (int, optional): Design variable number to assign to tube diameter (keyword argument). Defaults to -1
            (i.e. no design variable).
        dlb (float or complex, optional): Lower bound on tube diameter (keyword argument). Defaults to 0.0.
        dub (float or complex, optional): Upper bound on tube diameter (keyword argument). Defaults to 10.0.
        t (float or complex, optional): Tube wall thickness (keyword argument). Defaults to 0.1.
        tNum (int, optional): Design variable number to assign to wall thickness (keyword argument). Defaults to -1
            (i.e. no design variable).
        tlb (float or complex, optional): Lower bound on wall thickness (keyword argument). Defaults to 0.0.
        tub (float or complex, optional): Upper bound on wall thickness (keyword argument). Defaults to 10.0.
    """
    def __cinit__(self, *args, **kwargs):
        cdef TACSMaterialProperties *props = NULL
        cdef TacsScalar d = 1.0
        cdef int dNum = -1
        cdef TacsScalar dlb = 0.0
        cdef TacsScalar dub = 10.0
        cdef TacsScalar t = 0.1
        cdef int tNum = -1
        cdef TacsScalar tlb = 0.0
        cdef TacsScalar tub = 10.0

        if len(args) >= 1:
            props = (<MaterialProperties>args[0]).ptr
        if 'd' in kwargs:
            d = kwargs['d']
        if 'dNum' in kwargs:
            dNum = kwargs['dNum']
        if 'dlb' in kwargs:
            dlb = kwargs['dlb']
        if 'dub' in kwargs:
            dub = kwargs['dub']
        if 't' in kwargs:
            t = kwargs['t']
        if 'tNum' in kwargs:
            tNum = kwargs['tNum']
        if 'tlb' in kwargs:
            tlb = kwargs['tlb']
        if 'tub' in kwargs:
            tub = kwargs['tub']

        if props is not NULL:
            self.cptr = new TACSIsoTubeBeamConstitutive(props, d, t, dNum, tNum,
                                                        dlb, dub, tlb, tub)
            self.ptr = self.cptr
            self.ptr.incref()
        else:
            self.ptr = NULL
            self.cptr = NULL

cdef class IsoRectangleBeamConstitutive(BeamConstitutive):
    """
    Timoshenko theory based constitutive object for a solid rectangular beam.

    The thickness dimension is assumed to measured along the beam's local reference axis,
    the width is perpindicular to the reference axis.

    Args:
        props (MaterialProperties): The material property.
        width (float or complex, optional): Cross-section width (keyword argument). Defaults to 1.0.
        wNum (int, optional): Design variable number to assign to width (keyword argument). Defaults to -1
            (i.e. no design variable).
        wlb (float or complex, optional): Lower bound on width (keyword argument). Defaults to 0.0.
        wub (float or complex, optional): Upper bound on width diameter (keyword argument). Defaults to 10.0.
        t (float or complex, optional): Cross-section thickness (keyword argument). Defaults to 0.1.
        tNum (int, optional): Design variable number to assign to thickness (keyword argument). Defaults to -1
            (i.e. no design variable).
        tlb (float or complex, optional): Lower bound on thickness (keyword argument). Defaults to 0.0.
        tub (float or complex, optional): Upper bound on thickness (keyword argument). Defaults to 10.0.
    """
    def __cinit__(self, *args, **kwargs):
        cdef TACSMaterialProperties *props = NULL
        cdef TacsScalar w = 1.0
        cdef int wNum = -1
        cdef TacsScalar wlb = 0.0
        cdef TacsScalar wub = 10.0
        cdef TacsScalar t = 0.1
        cdef int tNum = -1
        cdef TacsScalar tlb = 0.0
        cdef TacsScalar tub = 10.0

        if len(args) >= 1:
            props = (<MaterialProperties>args[0]).ptr
        if 'w' in kwargs:
            w = kwargs['w']
        if 'wNum' in kwargs:
            wNum = kwargs['wNum']
        if 'wlb' in kwargs:
            wlb = kwargs['wlb']
        if 'wub' in kwargs:
            wub = kwargs['wub']
        if 't' in kwargs:
            t = kwargs['t']
        if 'tNum' in kwargs:
            tNum = kwargs['tNum']
        if 'tlb' in kwargs:
            tlb = kwargs['tlb']
        if 'tub' in kwargs:
            tub = kwargs['tub']

        if props is not NULL:
            self.cptr = new TACSIsoRectangleBeamConstitutive(props, w, t, wNum, tNum,
                                                             wlb, wub, tlb, tub)
            self.ptr = self.cptr
            self.ptr.incref()
        else:
            self.ptr = NULL
            self.cptr = NULL

cdef class GeneralMassConstitutive(Constitutive):
    """
    This is the base class for the fully general point mass constitutive objects.
    Assumes 6 dofs (3 translations + 3 rotations).
    The mass properties of this object are specified using a symmetric 6 x 6 mass matrix, as shown below.

        | M[ 0] M[ 1] M[ 2] M[ 3] M[ 4] M[ 5]
        | M[ 1] M[ 6] M[ 7] M[ 8] M[ 9] M[10]
        | M[ 2] M[ 7] M[11] M[12] M[13] M[14]
        | M[ 3] M[ 8] M[12] M[15] M[16] M[17]
        | M[ 4] M[ 9] M[13] M[16] M[18] M[19]
        | M[ 5] M[10] M[14] M[17] M[19] M[20]

    Args:
        M (array-like[float or complex]): Flattened array containing one side of symmetric mass matrix entries.
    """
    def __cinit__(self, **kwargs):
        cdef TacsScalar M[21]
        if 'M' in kwargs:
            _M = kwargs['M']
            assert isinstance(_M, list) or isinstance(_M, np.ndarray)
            assert len(_M) == 21
            for i in range(21):
                M[i] = _M[i]
        else:
            for i in range(21):
                M[i] = 0.0
        self.cptr = new TACSGeneralMassConstitutive(M)

        self.ptr = self.cptr
        self.ptr.incref()

cdef class PointMassConstitutive(GeneralMassConstitutive):
    """
    This is the base class for the traditional point mass constitutive objects with no translation-rotation coupling.
    Assumes 6 dofs.

     .. note::
        TACS uses a negative sign convention in the product of inertia definition, for example:
            I12 = -int[x1 * x2 * dm]

        The moments of inertia are always positive, as usual:
            I11 = int[(x2^2 + x3^2) * dm]

    Args:
        m (float or complex, optional): Mass value (keyword argument). Defaults to 0.0.
        I11 (float or complex, optional): Moment of inertia in '1' direction (keyword argument). Defaults to 0.0.
        I22 (float or complex, optional): Moment of inertia in '2' direction (keyword argument). Defaults to 0.0.
        I33 (float or complex, optional): Moment of inertia in '3' direction (keyword argument). Defaults to 0.0.
        I12 (float or complex, optional): Moment of inertia in '1-2' plane (keyword argument). Defaults to 0.0.
        I13 (float or complex, optional): Moment of inertia in '1-3' plane (keyword argument). Defaults to 0.0.
        I23 (float or complex, optional): Moment of inertia in '2-3' plane (keyword argument). Defaults to 0.0.
    """
    def __cinit__(self, **kwargs):
        cdef TacsScalar m = 0.0
        cdef TacsScalar I11 = 0.0
        cdef TacsScalar I22 = 0.0
        cdef TacsScalar I33 = 0.0
        cdef TacsScalar I12 = 0.0
        cdef TacsScalar I13 = 0.0
        cdef TacsScalar I23 = 0.0

        if 'm' in kwargs:
            m = kwargs['m']
        if 'I11' in kwargs:
            I11 = kwargs['I11']
        if 'I22' in kwargs:
            I22 = kwargs['I22']
        if 'I33' in kwargs:
            I33 = kwargs['I33']
        if 'I12' in kwargs:
            I12 = kwargs['I12']
        if 'I13' in kwargs:
            I13 = kwargs['I13']
        if 'I23' in kwargs:
            I23 = kwargs['I23']

        self.cptr = new TACSPointMassConstitutive(m, I11, I22, I33, I12, I13, I23)

        self.ptr = self.cptr
        self.ptr.incref()

cdef class GeneralSpringConstitutive(Constitutive):
    """
    This is the base class for the fully general spring constitutive objects.
    Assumes 6 dofs (3 translations + 3 rotations).
    The spring properties of this object are specified using a symmetric 6 x 6 stiffness matrix, as shown below.

        | K[ 0] K[ 1] K[ 2] K[ 3] K[ 4] K[ 5]
        | K[ 1] K[ 6] K[ 7] K[ 8] K[ 9] K[10]
        | K[ 2] K[ 7] K[11] K[12] K[13] K[14]
        | K[ 3] K[ 8] K[12] K[15] K[16] K[17]
        | K[ 4] K[ 9] K[13] K[16] K[18] K[19]
        | K[ 5] K[10] K[14] K[17] K[19] K[20]

    Args:
        K (array-like[float or complex]): Flattened array containing one side of symmetric stiffness matrix entries.
    """
    def __cinit__(self, **kwargs):
        cdef TacsScalar K[21]
        if 'K' in kwargs:
            _K = kwargs['K']
            assert isinstance(_K, list) or isinstance(_K, np.ndarray)
            assert len(_K) == 21
            for i in range(21):
                K[i] = _K[i]
        else:
            for i in range(21):
                K[i] = 0.0
        self.cptr = new TACSGeneralSpringConstitutive(K)

        self.ptr = self.cptr
        self.ptr.incref()

cdef class DOFSpringConstitutive(GeneralSpringConstitutive):
    """
    This is the base class for the traditional spring constitutive objects with no dof coupling.
    Assumes 6 dofs.

    Args:
        k (array-like[float or complex]): Stiffness values for all 6 dofs.
    """
    def __cinit__(self, **kwargs):
        cdef TacsScalar k[6]
        if 'k' in kwargs:
            _k = kwargs['k']
            assert isinstance(_k, list) or isinstance(_k, np.ndarray)
            assert len(_k) == 6
            for i in range(6):
                k[i] = _k[i]
        else:
            for i in range(6):
                k[i] = 0.0
        self.cptr = new TACSDOFSpringConstitutive(k)

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