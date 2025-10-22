/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_OBJECT_H
#define TACS_OBJECT_H

/*!
  The following class implements a basic reference counting scheme
  for memory allocation/deallocation.

  The following assumptions are made about references:
  - References returned from functions are always borrowed
  (except for constructors which own the new reference -- this is
  the way I'd like to do it, but SWIG doesn't seem to support it easily)
  - References passed into functions have no effect -
  references are never stolen.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "TacsComplexStep.h"
#include "mpi.h"

extern MPI_Op TACS_MPI_MIN;
extern MPI_Op TACS_MPI_MAX;

/*
  Use the cplx type for TacsComplex
*/
typedef std::complex<double> TacsComplex;
typedef double TacsReal;

/*
  Define the basic scalar type TacsScalar
*/
#ifdef TACS_USE_COMPLEX
#define TACS_MPI_TYPE MPI_DOUBLE_COMPLEX
typedef TacsComplex TacsScalar;
#else
#define TACS_MPI_TYPE MPI_DOUBLE
typedef double TacsScalar;
#endif

/*
  Define the macro to add flop counts. This does not work for threaded
  implementations. Don't use it in threaded code!

  Note that this may only log FLOPs from certain parts of the code.
  I've tried to log all operations within the linear algebra portions of
  the code, but the analysis (e.g. residual/matrix computation) is much
  more difficult.
*/
extern double tacs_local_flop_count;

// Zero the number of counted flops
void TacsZeroNumFlops();

// Retrieve the total number of counted flops
double TacsGetNumFlops();

// Macro to record the number of flops in an operation
#ifdef TACS_LOG_FLOPS
#define TacsAddFlops(flop) (tacs_local_flop_count += (flop));
#else
#define TacsAddFlops(flop)
#endif

/*
  Set up the define statements for beginning/ending a namespace
  environment
*/
#define TACS_BEGIN_NAMESPACE(a) namespace a {
#define TACS_END_NAMESPACE }

/**
  Initialize some static variables related to MPI
*/
void TacsInitialize();

/**
  Check if TacsInitialize has been called
*/
int TacsIsInitialized();

/**
  Clean up from the Tacs initialization
*/
void TacsFinalize();

/**
  TACSObject: The base class for all TACS objects to enable reference
  counting. In most cases this is sufficient to free any allocated
  memory.
*/
class TACSObject {
 public:
  TACSObject();
  virtual ~TACSObject();

  /**
    Increase the reference count function
  */
  void incref();

  /**
    Decrease the reference count
  */
  void decref();

  /**
    Return the reference count
  */
  int refcount();

  /**
    Return the name of the object
  */
  virtual const char *getObjectName();

 private:
  int ref_count;
  static const char *tacsDefault;
};

/**
  Information about the number of threads to use in a computation.
  This should only be allocated by the TACSAssembler object. The
  number of threads is volitile in the sense that it can change
  between subsequent calls.
*/
class TACSThreadInfo : public TACSObject {
 public:
  static const int TACS_MAX_NUM_THREADS = 16;

  TACSThreadInfo(int _num_threads);
  ~TACSThreadInfo() {}

  void setNumThreads(int _num_threads);
  int getNumThreads();

 private:
  int num_threads;
};

#endif
