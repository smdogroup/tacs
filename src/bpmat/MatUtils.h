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

#ifndef TACS_MAT_UTILS_H
#define TACS_MAT_UTILS_H

/*
  The following are some useful functions that may be used
  to manipulate graph data structures.
*/

#include "TACSObject.h"

TACS_BEGIN_NAMESPACE(matutils)

/*!
  Extend an array and copy values from the old to the new array
*/
void ExtendArray( int **_array, int oldlen, int newlen );  
void ExtendArray( TacsScalar **_array, int oldlen, int newlen );

/*!
  Sort and uniquify a CSR data structure.

  nvars:        the number of variables
  rowp:         pointer into the rows of of the matrix
  cols:         the column indices
  remove_diag:  remove any diagonal entries encountered
*/
void SortAndUniquifyCSR( int nvars, int *rowp, int *cols, 
                         int remove_diag=0 );

/*
  Reorder locally based on RCM-reordering
*/
int ComputeRCMOrder( const int nvars, const int *rowp, const int *cols,
                     int *rcm_order, int root, int n_rcm_iters );
int ComputeRCMLevSetOrder( const int nvars, const int *rowp, 
                           const int *cols, int *rcm_vars, int *levset, 
                           int root );

/*
  Compute the multicolor order via a greedy algorithm
*/
int ComputeSerialMultiColor( const int nvars, const int *rowp,
                             const int *cols, int *colors, 
                             int *new_vars );

TACS_END_NAMESPACE

#endif
