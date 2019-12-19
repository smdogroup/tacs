/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_UTILITIES_H
#define TACS_UTILITIES_H

#include "TACSObject.h"

int TacsIntegerComparator( const void *a, const void *b );

/**
  Sort an array and return the number of unique design variables within
  that array start is the index of the first non-negative value and the
  return value is the number of non-negative design variables
  &dvNums[start] is the start of the new array of unique design vars...
*/
int TacsUniqueSort( int len, int *array );

/**
  Locate a value within a sorted array
*/
inline int* TacsSearchArray( int value, int len, const int *array ){
  return (int*)bsearch(&value, array, len, sizeof(int), TacsIntegerComparator);
}

/**
  Merge two sorted arrays into a single sorted array, in place.

  Two part algorithm:
  1. Find the number of duplicates between a and b
  2. Run through the list backwards placing elements into a[]
  when appropriate.

  Memory requirements: note that len(a) >= na + nb
*/
int TacsMergeSortedArrays( int na, int *a, int nb, const int *b );

/**
  Find the interval such that the given index satisfies:

  intv[k] <= index < intv[k+1]

  The intervals must be non-decreasing. Note that len is equal
  to the length of the intv array which is one more than the
  total number of intervals.
*/
int TacsFindInterval( int index, int len, const int intv[] );

/**
  Match the intervals in a list of sorted variables.

  Given the intervals, the range of the intervals and
  the list of sorted variables, local a point such that

  vars[ext_ptr[n]] <= ownerRange[n]
*/
void TacsMatchIntervals( int size, const int range[],
                         int nvars, const int vars[], int ext_ptr[] );

/*!
  Extend an array and copy values from the old to the new array
*/
void TacsExtendArray( int **_array, int oldlen, int newlen );
void TacsExtendArray( TacsScalar **_array, int oldlen, int newlen );

/*!
  Sort and uniquify a CSR data structure.

  nvars:        the number of variables
  rowp:         pointer into the rows of of the matrix
  cols:         the column indices
  remove_diag:  remove any diagonal entries encountered
*/
void TacsSortAndUniquifyCSR( int nvars, int *rowp, int *cols,
                             int remove_diag=0 );

/*
  Reorder locally based on RCM-reordering
*/
int TacsComputeRCMOrder( const int nvars, const int *rowp, const int *cols,
                         int *rcm_order, int root, int n_rcm_iters );

/*
  Compute the multicolor order via a greedy algorithm
*/
int TacsComputeSerialMultiColor( const int nvars, const int *rowp,
                                 const int *cols, int *colors,
                                 int *new_vars );

#endif // TACS_UTILITIES_H
