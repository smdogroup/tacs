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
from libc.string cimport const_char
from libcpp cimport bool

# Import numpy
cimport numpy as np
import numpy as np

# Import TACS c++ headers
from tacs.cpp_headers.TACS cimport *

cdef inline char* convert_to_chars(s):
    if isinstance(s, unicode):
        s = (<unicode>s).encode('utf8')
    return s

cdef inline convert_bytes_to_str(bytes s):
    return s.decode('utf8')

cdef class BcMap:
    cdef TACSBcMap *ptr

cdef inline _init_BcMap(TACSBcMap *ptr):
    bcmap = BcMap()
    bcmap.ptr = ptr
    bcmap.ptr.incref()
    return bcmap

cdef class NodeMap:
    cdef TACSNodeMap *ptr

cdef inline _init_NodeMap(TACSNodeMap *ptr):
    vmap = NodeMap()
    vmap.ptr = ptr
    vmap.ptr.incref()
    return vmap

cdef class VecIndices:
    cdef TACSBVecIndices *ptr

cdef inline _init_VecIndices(TACSBVecIndices *ptr):
    indices = VecIndices()
    indices.ptr = ptr
    indices.ptr.incref()
    return indices

cdef class VecInterp:
    cdef TACSBVecInterp *ptr

cdef inline _init_VecInterp(TACSBVecInterp *ptr):
    interp = VecInterp()
    interp.ptr = ptr
    interp.ptr.incref()
    return interp

cdef class Vec:
    cdef TACSVec *ptr
    cdef TACSBVec* getBVecPtr(self)

cdef inline _init_Vec(TACSVec *ptr):
    vec = Vec()
    vec.ptr = ptr
    vec.ptr.incref()
    return vec

cdef class Mat:
    cdef TACSMat *ptr

cdef inline _init_Mat(TACSMat *ptr):
    mat = Mat()
    mat.ptr = ptr
    mat.ptr.incref()
    return mat

cdef class Pc:
    cdef TACSPc *ptr

cdef inline _init_Pc(TACSPc *ptr):
    pc = Pc()
    pc.ptr = ptr
    pc.ptr.incref()
    return pc

cdef class Mg(Pc):
    cdef TACSMg *mg

cdef inline _init_Mg(TACSMg *ptr):
    mg = Mg()
    mg.mg = ptr
    mg.ptr = ptr
    mg.mg.incref()
    return mg

cdef class KSM:
    cdef TACSKsm *ptr

cdef class ElementBasis:
    cdef TACSElementBasis *ptr

cdef inline _init_ElementBasis(TACSElementBasis *ptr):
    basis = ElementBasis()
    basis.ptr = ptr
    basis.ptr.incref()
    return basis

cdef class ElementModel:
    cdef TACSElementModel *ptr
    cdef object con

cdef inline _init_ElementModel(TACSElementModel *ptr):
    model = ElementModel()
    model.ptr = ptr
    model.ptr.incref()
    return model

cdef class Element:
    cdef TACSElement *ptr
    cdef object con
    cdef object transform

cdef inline _init_Element(TACSElement *ptr):
    elem = Element()
    elem.ptr = ptr
    elem.ptr.incref()
    return elem

cdef class Function:
    cdef TACSFunction *ptr

cdef class Constitutive:
    cdef TACSConstitutive *ptr
    cdef object props
    cdef int nastranID

cdef inline _init_Constitutive(TACSConstitutive *ptr):
    cons = Constitutive()
    cons.ptr = ptr
    cons.ptr.incref()
    return cons

cdef class Assembler:
    cdef TACSAssembler *ptr

cdef inline _init_Assembler(TACSAssembler *ptr):
    tacs = Assembler()
    tacs.ptr = ptr
    tacs.ptr.incref()
    return tacs

cdef class JacobiDavidsonOperator:
    cdef TACSJacobiDavidsonOperator *ptr
