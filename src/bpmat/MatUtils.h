#ifndef TACS_MAT_UTILS_H
#define TACS_MAT_UTILS_H

/*
  The following are some useful functions that may be used
  to manipulate graph data structures.

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

#include "TACSObject.h"

TACS_BEGIN_NAMESPACE(matutils)
/*!
  Extend an array and copy values from the old to the new array
*/
void ExtendArray( int ** _array, int oldlen, int newlen );  
void ExtendArray( TacsScalar ** _array, int oldlen, int newlen );

/*!
  Sort and uniquify a CSR data structure.
  nvars/rowp/cols == Input/output CSR data structure
  nodiag == flag, remove the diagonal entry from the array.
*/
void SortAndUniquifyCSR( int nvars, int * rowp, int * cols, int nodiag );

/*
  Reorder based on RCM-reordering
*/
int ComputeRCMOrder( const int nvars, const int * rowp, const int * cols,
                     int * rcm_order, int root, int n_rcm_iters );
int ComputeRCMLevSetOrder( const int nvars, const int * rowp, 
                           const int * cols, int * rcm_vars, int * levset, 
                           int root );

void ElementListToVarList( const int * elems, const int * elemptr, int nelems,
                           int lower, int upper,
                           int * _varelems, int * _varptr );
  
void ElementListToExtVarList( const int * elems, const int * elemptr, 
                              int nelems, const int * extvars, int nextvars,
                              int lower, int upper,
                              int * _varelems, int * _varptr );

TACS_END_NAMESPACE

#endif
