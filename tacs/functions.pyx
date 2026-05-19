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

import warnings
from enum import IntEnum

# Import the definition required for const strings
from libc.string cimport const_char
from libc.stdlib cimport malloc, free

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

# Import the definitions
from tacs.TACS cimport *

# Include the definitions
include "TacsDefs.pxi"

# Include the mpi4py header
cdef extern from "mpi-compat.h":
    pass

cdef class StructuralMass(Function):
    """
    Evaluates the structural mass of the elements.

    Args:
        assembler (Assembler): TACS Assembler object that will evaluating this function.
    """
    def __cinit__(self, Assembler assembler):
        """
        Wrap the function StructuralMass
        """
        self.ptr = new TACSStructuralMass(assembler.ptr)
        self.ptr.incref()
        return

cdef class CenterOfMass(Function):
    """
    Evaluates the center of mass of the elements along specified axis.

    .. warning:: This function is only appropriate for static analyses and
        may give inconsistent results for transient problems.

    Args:
        assembler (Assembler): TACS Assembler object that will evaluating this function.
        direction (array-like[double], optional):
          3d vector specifying which direction to project cg position onto (keyword argument).
          Defaults to [0.0, 0.0, 0.0].
    """
    def __cinit__(self, Assembler assembler, **kwargs):
        """
        Wrap the function CenterOfMass
        """
        cdef double d[3]
        d[0] = d[1] = d[2] = 0.0

        if 'direction' in kwargs:
            dir = kwargs['direction']
            # Check if dir is a list or numpy array
            if isinstance(dir, list) or isinstance(dir, np.ndarray):
                dim = min(3, len(dir))
                for i in range(dim):
                    d[i] = dir[i]

        self.ptr = new TACSCenterOfMass(assembler.ptr, d)
        self.ptr.incref()
        return

cdef class MomentOfInertia(Function):
    """
    Evaluates the moment of inertia of the elements about origin or center of mass projected onto two input vectors:

        I_out = vec1^T * I_tensor * vec2

    Where I_tensor is the moment of inertia tensor in the global axis given below:

        |   Ixx  Ixy  Ixz
        |   Iyx  Iyy  Iyz
        |   Izx  Izy  Izz

    .. note::
        TACS uses a negative sign convention in the product of inertia definition, for example:
            Ixy = -int[x * y * dm]

        The moments of inertia are always positive, as usual:
            Ixx = int[(y^2 + z^2) * dm]

    .. warning:: This function is only appropriate for static analyses and
        may give inconsistent results for transient problems.

    Args:
        assembler (Assembler): TACS Assembler object that will evaluating this function.
        direction1 (array-like[double], optional):
          3d vector specifying first direction to project moment of inertia tensor onto (keyword argument).
          Defaults to [0.0, 0.0, 0.0].
        direction2 (array-like[double], optional):
          3d vector specifying second direction to project moment of inertia tensor onto (keyword argument).
          Defaults to [0.0, 0.0, 0.0].
        aboutCM (bool): Flag specifying whether moment of inertia should be taken
          about origin (False) or center of mass (True) (keyword argument). Defaults to False.
    """
    def __cinit__(self, Assembler assembler, **kwargs):
        """
        Wrap the function MomentOfInertia
        """
        cdef int cmFlag;
        cdef double d1[3], d2[3]
        d1[0] = d1[1] = d1[2] = 0.0
        d2[0] = d2[1] = d2[2] = 0.0

        if 'direction1' in kwargs:
            dir = kwargs['direction1']
            # Check if dir is a list or numpy array
            if isinstance(dir, list) or isinstance(dir, np.ndarray):
                dim = min(3, len(dir))
                for i in range(dim):
                    d1[i] = dir[i]

        if 'direction2' in kwargs:
            dir = kwargs['direction2']
            # Check if dir is a list or numpy array
            if isinstance(dir, list) or isinstance(dir, np.ndarray):
                dim = min(3, len(dir))
                for i in range(dim):
                    d2[i] = dir[i]

        cmFlag = kwargs.get('aboutCM', False);

        self.ptr = new TACSMomentOfInertia(assembler.ptr, d1, d2, cmFlag)
        self.ptr.incref()
        return

cdef class EnclosedVolume(Function):
    """
    Evaluates the volume enclosed by the elements.

    Args:
        assembler (Assembler): TACS Assembler object that will evaluating this function.
    """
    def __cinit__(self, Assembler assembler):
        """
        Wrap the function EnclosedVolume
        """
        self.ptr = new TACSEnclosedVolume(assembler.ptr)
        self.ptr.incref()
        return

cdef class Compliance(Function):
    """
    Evaluate the compliance of the structure.

    This evaluates the compliance within the structure based on an
    integration of the strain energy in each element, not by the product
    of the load vector with the displacement.

    Args:
        assembler (Assembler): TACS Assembler object that will evaluating this function.
    """
    cdef TACSCompliance *cptr
    def __cinit__(self, Assembler assembler):
        """
        Wrap the function Compliance
        """
        self.cptr = new TACSCompliance(assembler.ptr)
        self.ptr = self.cptr
        self.ptr.incref()
        return

    def setComplianceType(self, int compliance_type):
        """
        Set the type of compliance value to use
        """
        self.cptr.setComplianceType(compliance_type)
        return

cdef class AverageTemperature(Function):
    """
    Evaluate the spatial average of the temperature of the model.

    Args:
        assembler (Assembler): TACS Assembler object that will evaluating this function.
        volume (float, optional): Normalization factor used in averaging process. Defaults to 1.0.
    """
    def __cinit__(self, Assembler assembler, **kwargs):
        """
        Wrap the function AverageTemperature
        """
        cdef double volume = 1.0

        if 'volume' in kwargs:
            volume = kwargs['volume']

        self.ptr = new TACSAverageTemperature(assembler.ptr, volume)
        self.ptr.incref()
        return

class KSAggregationType(IntEnum):
    """
    Aggregation type for KS functions. Mirrors ``KSAggregationType`` from ``TACSFunction.h``.

    ``KS_DISCRETE_AVERAGE`` is only valid for :class:`KSFailure`.
    """
    KS_DISCRETE = _KSAGG_DISCRETE
    KS_CONTINUOUS = _KSAGG_CONTINUOUS
    PNORM_DISCRETE = _KSAGG_PNORM_DISCRETE
    PNORM_CONTINUOUS = _KSAGG_PNORM_CONTINUOUS
    KS_DISCRETE_AVERAGE = _KSAGG_DISCRETE_AVERAGE


cdef class KSTemperature(Function):
    """
    The following class implements the methods necessary
    to calculate the Kreisselmeier–Steinhauser (KS) aggregation
    of temperature over the domain of some finite element model.
    The KS aggregation gives a smooth and differentiable approximation to the
    maximum value.

    Args:
        assembler (Assembler): TACS Assembler object that will evaluating this function.
        ksWeight (float, optional): The ks weight used in the calculation (keyword argument). Defaults to 80.0.
        ksAggregationType (functions.KSAggregationType, optional): The type of KS aggregation to be used.
            Defaults to ``functions.KSAggregationType.KS_CONTINUOUS``. ``DISCRETE_AVERAGE`` is not supported.
        ftype (str, optional): Deprecated. Use ``ksAggregationType=functions.KSAggregationType.<VALUE>`` instead.
    """

    cdef TACSKSTemperature *kstptr
    def __cinit__(self, Assembler assembler, **kwargs):
        cdef double ksWeight = 80.0
        cdef double alpha = 1.0

        if 'ksWeight' in kwargs:
            ksWeight = kwargs['ksWeight']

        if 'alpha' in kwargs:
            alpha = kwargs['alpha']

        self.kstptr = new TACSKSTemperature(assembler.ptr, ksWeight, alpha)
        self.ptr = self.kstptr
        self.ptr.incref()

        if 'ftype' in kwargs:
            # Deprecated in v3.12, remove in v3.14
            warnings.warn(
                "The 'ftype' string kwarg is deprecated as of v3.12 and will be removed "
                "in v3.14. Use 'ksAggregationType=functions.KSAggregationType.<VALUE>' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            if 'ksAggregationType' in kwargs:
                raise ValueError(
                    "Cannot specify both 'ksAggregationType' and deprecated 'ftype' kwarg."
                )

            _str_map = {
                'discrete': KSAggregationType.KS_DISCRETE,
                'continuous': KSAggregationType.KS_CONTINUOUS,
                'pnorm-discrete': KSAggregationType.PNORM_DISCRETE,
                'pnorm-continuous': KSAggregationType.PNORM_CONTINUOUS,
            }
            ftype = kwargs['ftype']
            try:
                ksAggregationType = _str_map[ftype.lower()]
            except KeyError:
                raise ValueError(
                    f"Unknown ftype {ftype!r}. Valid options: {sorted(_str_map)}"
                ) from None
        else:
            ksAggregationType = kwargs.get('ksAggregationType', KSAggregationType.KS_CONTINUOUS)
        self.setKSAggregationType(ksAggregationType)

    def setKSAggregationType(self, ksAggregationType):
        """
        Set the type of KS aggregation.

        Args:
            ksAggregationType (functions.KSAggregationType):
                The aggregation type.

        Raises:
            ValueError: If ``ksAggregationType`` is ``KSAggregationType.KS_DISCRETE_AVERAGE``.
        """

        ksAggregationType = KSAggregationType(ksAggregationType)
        if ksAggregationType == KSAggregationType.KS_DISCRETE_AVERAGE:
            raise ValueError(
                "KSAggregationType.KS_DISCRETE_AVERAGE is only valid for KSFailure, not KSTemperature."
            )
        self.kstptr.setKSAggregationType(<_CKSAggregationType><int>ksAggregationType)

    def setParameter(self, double ksparam):
        self.kstptr.setParameter(ksparam)


cdef class KSFailure(Function):
    """
    The following class implements the methods necessary to calculate
    the Kreisselmeier–Steinhauser (KS) aggregation of either a stress
    or strain failure criteria over the domain of some finite element model.
    The KS aggregation gives a smooth and differentiable approximation to the
    maximum value.

    The failure load is calculated using the strain-based failure
    criteria from the base :class:`~TACS.Constitutive` class which requires linear and
    constant components of the strain to determine the failure load.

    For most elements, unity is generally considered to be the threshold value for failure.
    Meaning if this function returns a value > 1.0, at least one of the elements has exceeded
    its strength criteria. While values < 1.0 implies all elements are within their strength criteria.

    Args:
        assembler (Assembler): TACS Assembler object that will evaluating this function.
        ksWeight (float, optional): The ks weight used in the calculation (keyword argument). Defaults to 80.0.
        safetyFactor (float, optional):
            The safety factor to apply to loads before computing the failure (keyword argument). Defaults to 1.0.
        ksAggregationType (functions.KSAggregationType, optional): The type of KS aggregation to be used.
            Defaults to ``functions.KSAggregationType.KS_CONTINUOUS``.
        ftype (str, optional): Deprecated. Use ``ksAggregationType=functions.KSAggregationType.<VALUE>`` instead.
    """

    cdef TACSKSFailure *ksptr
    def __cinit__(self, Assembler assembler, **kwargs):
        cdef double ksWeight = 80.0
        cdef double alpha = 1.0
        cdef double safetyFactor = 1.0

        if 'ksWeight' in kwargs:
            ksWeight = kwargs['ksWeight']

        if 'alpha' in kwargs:
            alpha = kwargs['alpha']

        if 'safetyFactor' in kwargs:
            safetyFactor = kwargs['safetyFactor']

        self.ksptr = new TACSKSFailure(assembler.ptr, ksWeight, alpha, safetyFactor)
        self.ptr = self.ksptr
        self.ptr.incref()

        if 'ftype' in kwargs:
            # Deprecated in v3.12, remove in v3.14
            warnings.warn(
                "The 'ftype' string kwarg is deprecated as of v3.12 and will be removed "
                "in v3.14. Use 'ksAggregationType=functions.KSAggregationType.<VALUE>' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            if 'ksAggregationType' in kwargs:
                raise ValueError(
                    "Cannot specify both 'ksAggregationType' and deprecated 'ftype' kwarg."
                )

            _str_map = {
                'discrete': KSAggregationType.KS_DISCRETE,
                'continuous': KSAggregationType.KS_CONTINUOUS,
                'pnorm-discrete': KSAggregationType.PNORM_DISCRETE,
                'pnorm-continuous': KSAggregationType.PNORM_CONTINUOUS,
                'discrete-average': KSAggregationType.KS_DISCRETE_AVERAGE,
            }
            ftype = kwargs['ftype']
            try:
                ksAggregationType = _str_map[ftype.lower()]
            except KeyError:
                raise ValueError(
                    f"Unknown ftype {ftype!r}. Valid options: {sorted(_str_map)}"
                ) from None
        else:
            ksAggregationType = kwargs.get('ksAggregationType', KSAggregationType.KS_CONTINUOUS)
        self.setKSAggregationType(ksAggregationType)

    def setKSAggregationType(self, ksAggregationType):
        """
        Set the type of KS aggregation.

        Args:
            ksAggregationType (functions.KSAggregationType):
                The aggregation type.
        """
        ksAggregationType = KSAggregationType(ksAggregationType)
        self.ksptr.setKSAggregationType(<_CKSAggregationType><int>ksAggregationType)

    def setParameter(self, double ksparam):
        self.ksptr.setParameter(ksparam)


cdef class KSDisplacement(Function):
    """
    The following class implements the methods to calculate the
    Kreisselmeier–Steinhauser (KS) aggregation of the displacement in
    a particular direction over the domain of some finite element model.
    The KS aggregation gives a smooth and differentiable approximation to the
    maximum value.

    Args:
        assembler (Assembler): TACS Assembler object that will evaluating this function.
        ksWeight (float, optional): The ks weight used in the calculation (keyword argument). Defaults to 80.0.
        direction (array-like[double], optional):
          3d vector specifying which direction to project displacements in for KS aggregation (keyword argument).
          Defaults to [0.0, 0.0, 0.0].
        ksAggregationType (functions.KSAggregationType, optional): The type of KS aggregation to be used.
          Defaults to ``functions.KSAggregationType.KS_CONTINUOUS``. ``DISCRETE_AVERAGE`` is not supported.
        ftype (str, optional): Deprecated. Use ``ksAggregationType=functions.KSAggregationType.<VALUE>`` instead.
    """

    cdef TACSKSDisplacement *ksptr
    def __cinit__(self, Assembler assembler, **kwargs):
        cdef double ksWeight = 80.0
        cdef double alpha = 1.0
        cdef double d[3]
        d[0] = d[1] = d[2] = 0.0

        if 'ksWeight' in kwargs:
            ksWeight = kwargs['ksWeight']

        if 'alpha' in kwargs:
            alpha = kwargs['alpha']

        if 'direction' in kwargs:
            dir = kwargs['direction']
            # Check if dir is a list or numpy array
            if isinstance(dir, list) or isinstance(dir, np.ndarray):
                dim = min(3, len(dir))
                for i in range(dim):
                    d[i] = dir[i]

        self.ksptr = new TACSKSDisplacement(assembler.ptr, ksWeight, d, alpha)
        self.ptr = self.ksptr
        self.ptr.incref()

        if 'ftype' in kwargs:
            # Deprecated in v3.12, remove in v3.14
            warnings.warn(
                "The 'ftype' string kwarg is deprecated as of v3.12 and will be removed "
                "in v3.14. Use 'ksAggregationType=functions.KSAggregationType.<VALUE>' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            if 'ksAggregationType' in kwargs:
                raise ValueError(
                    "Cannot specify both 'ksAggregationType' and deprecated 'ftype' kwarg."
                )

            _str_map = {
                'discrete': KSAggregationType.KS_DISCRETE,
                'continuous': KSAggregationType.KS_CONTINUOUS,
                'pnorm-discrete': KSAggregationType.PNORM_DISCRETE,
                'pnorm-continuous': KSAggregationType.PNORM_CONTINUOUS,
            }
            ftype = kwargs['ftype']
            try:
                ksAggregationType = _str_map[ftype.lower()]
            except KeyError:
                raise ValueError(
                    f"Unknown ftype {ftype!r}. Valid options: {sorted(_str_map)}"
                ) from None
        else:
            ksAggregationType = kwargs.get('ksAggregationType', KSAggregationType.KS_CONTINUOUS)
        self.setKSAggregationType(ksAggregationType)

    def setKSAggregationType(self, ksAggregationType):
        """
        Set the type of KS aggregation.

        Args:
            ksAggregationType (functions.KSAggregationType):
                The aggregation type.

        Raises:
            ValueError: If ``ksAggregationType`` is ``KSAggregationType.KS_DISCRETE_AVERAGE``.
        """

        ksAggregationType = KSAggregationType(ksAggregationType)
        if ksAggregationType == KSAggregationType.KS_DISCRETE_AVERAGE:
            raise ValueError(
                "KSAggregationType.KS_DISCRETE_AVERAGE is only valid for KSFailure, not KSDisplacement."
            )
        self.ksptr.setKSAggregationType(<_CKSAggregationType><int>ksAggregationType)

    def setParameter(self, double ksparam):
        self.ksptr.setParameter(ksparam)

# cdef class InducedFailure(Function):
#     cdef TACSInducedFailure *iptr
#     def __cinit__(self, Assembler assembler, double P):
#         """
#         Wrap the function InducedFailure
#         """
#         self.iptr = new TACSInducedFailure(assembler.ptr, P)
#         self.ptr = self.iptr
#         self.ptr.incref()
#         return

#     def setInducedType(self, ftype='exponential'):
#         if ftype == 'exponential':
#             self.iptr.setInducedType(INDUCED_EXPONENTIAL)
#         elif ftype == 'power':
#             self.iptr.setInducedType(INDUCED_POWER)
#         elif ftype == 'exponential-squared':
#             self.iptr.setInducedType(INDUCED_EXPONENTIAL_SQUARED)
#         elif ftype == 'power-squared':
#             self.iptr.setInducedType(INDUCED_POWER_SQUARED)
#         elif ftype == 'discrete-exponential':
#             self.iptr.setInducedType(INDUCED_DISCRETE_EXPONENTIAL)
#         elif ftype == 'discrete-power':
#             self.iptr.setInducedType(INDUCED_DISCRETE_POWER)
#         elif ftype == 'discrete-exponential-squared':
#             self.iptr.setInducedType(INDUCED_DISCRETE_EXPONENTIAL_SQUARED)
#         elif ftype == 'discrete-power-squared':
#             self.iptr.setInducedType(INDUCED_DISCRETE_POWER_SQUARED)

#     def setParameter(self, double param):
#         self.iptr.setParameter(param)

cdef class HeatFlux(Function):
    cdef TACSHeatFlux *hptr
    def __cinit__(self, Assembler assembler, list elem_index,
                  list surfaces):
        cdef int num_elems = len(elem_index)
        cdef int *elem_ind = NULL
        cdef int *surf = NULL

        elem_ind = <int*>malloc(num_elems*sizeof(int));
        surf = <int*>malloc(num_elems*sizeof(int));

        for i in range(num_elems):
            elem_ind[i] = <int>elem_index[i]
            surf[i] = <int>surfaces[i]
        self.hptr = new TACSHeatFlux(assembler.ptr, elem_ind, surf,
                                     num_elems)
        self.ptr = self.hptr
        self.ptr.incref()

        free(elem_ind)
        free(surf)
        return

# cdef class DisplacementIntegral(Function):
#     cdef TACSDisplacementIntegral *dptr
#     def __cinit__(self, Assembler assembler, dirs):
#         """
#         Wrap the function KSFailure
#         """
#         cdef TacsScalar _dirs[3]
#         _dirs[0] = dirs[0]
#         _dirs[1] = dirs[1]
#         _dirs[2] = dirs[2]
#         self.dptr = new TACSDisplacementIntegral(assembler.ptr, _dirs)
#         self.ptr = self.dptr
#         self.ptr.incref()
#         return
