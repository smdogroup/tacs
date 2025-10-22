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

#include "TACSObject.h"

/*
  Implementation of the reference counting TACSObject as well as
  initialization of MPI data for complex arithmetic
*/

double tacs_local_flop_count = 0.0;

double TacsGetNumFlops() {
#ifdef TACS_LOG_FLOPS
  return tacs_local_flop_count;
#else
  fprintf(stderr, "Warning: TACS not compiled with -DTACS_LOG_FLOPS\n");
  return 0.0;
#endif
}

void TacsZeroNumFlops() { tacs_local_flop_count = 0.0; }

/*
  These definite the min/max operations for complex values
*/
#ifdef TACS_USE_COMPLEX
void TacsMPIComplexMax(void *_in, void *_out, int *count, MPI_Datatype *data) {
  if (*data == MPI_DOUBLE_COMPLEX) {
    TacsScalar *in = (TacsScalar *)_in;
    TacsScalar *out = (TacsScalar *)_out;

    // Compare the real parts of the array
    for (int i = 0; i < *count; i++) {
      if (TacsRealPart(in[i]) >= TacsRealPart(out[i])) {
        out[i] = in[i];
      }
    }
  }
}

void TacsMPIComplexMin(void *_in, void *_out, int *count, MPI_Datatype *data) {
  if (*data == MPI_DOUBLE_COMPLEX) {
    TacsScalar *in = (TacsScalar *)_in;
    TacsScalar *out = (TacsScalar *)_out;

    // Compare the real parts of the array
    for (int i = 0; i < *count; i++) {
      if (TacsRealPart(in[i]) < TacsRealPart(out[i])) {
        out[i] = in[i];
      }
    }
  }
}
#endif

// Static flag to test if TacsInitialize has been called
static int TacsInitialized = 0;

MPI_Op TACS_MPI_MIN = MPI_MAX;
MPI_Op TACS_MPI_MAX = MPI_MIN;

void TacsInitialize() {
  if (!TacsInitialized) {
#ifdef TACS_USE_COMPLEX
    // Try to add the MPI reduction operator for MPI_DOUBLE_COMPLEX
    int commute = 1;
    MPI_Op_create(TacsMPIComplexMax, commute, &TACS_MPI_MAX);
    MPI_Op_create(TacsMPIComplexMin, commute, &TACS_MPI_MIN);
#else
    TACS_MPI_MAX = MPI_MAX;
    TACS_MPI_MIN = MPI_MIN;
#endif
  }
  TacsInitialized++;
}

int TacsIsInitialized() { return TacsInitialized; }

void TacsFinalize() {
  TacsInitialized--;

  if (TacsInitialized == 0) {
#ifdef TACS_USE_COMPLEX
    MPI_Op_free(&TACS_MPI_MAX);
    MPI_Op_free(&TACS_MPI_MIN);
#endif
  }
}

TACSObject::TACSObject() { ref_count = 0; }

TACSObject::~TACSObject() {}

/*
  Increase the reference count functions
*/
void TACSObject::incref() { ref_count++; }

/*
  Decrease the reference count
*/
void TACSObject::decref() {
  ref_count--;

  if (ref_count == 0) {
#ifdef TACS_DEBUG
    fprintf(stderr, "Deleting object %s\n", this->getObjectName());
#endif
    delete this;
  } else if (ref_count < 0) {
    fprintf(stderr, "Encountered a negative reference count for %s\n",
            this->getObjectName());
  }
}

//! Return the reference count
int TACSObject::refcount() { return ref_count; }

//! Return the name of the object
const char *TACSObject::getObjectName() { return tacsDefault; }

const char *TACSObject::tacsDefault = "TACSObject";

/*
  Implementation of the TACSThreadInfo object
*/
TACSThreadInfo::TACSThreadInfo(int _num_threads) {
  if (_num_threads > 1) {
    num_threads = _num_threads;
  } else {
    num_threads = 1;
  }
}

void TACSThreadInfo::setNumThreads(int _num_threads) {
  if (_num_threads > 1 && _num_threads <= TACS_MAX_NUM_THREADS) {
    num_threads = _num_threads;
  } else if (_num_threads <= 1) {
    num_threads = 1;
  } else if (_num_threads >= TACS_MAX_NUM_THREADS) {
    num_threads = TACS_MAX_NUM_THREADS;
  }
}

int TACSThreadInfo::getNumThreads() { return num_threads; }
