#include "FElibrary.h"
#include "MatUtils.h"

/*
  Use MatUtils

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

TACS_BEGIN_NAMESPACE(matutils)

/*!
  Extend the length of an integer array to a new length of array
*/
void ExtendArray( int ** _array, int oldlen, int newlen ){
  int * oldarray = *_array;
  int * newarray = new int[ newlen ];
  memcpy(newarray, oldarray, oldlen*sizeof(int));
  delete [] *_array;
  *_array = newarray;
}

/*
  Extend the length of a TacsScalar array to a new length
*/
void ExtendArray( TacsScalar ** _array, int oldlen, int newlen ){
  TacsScalar * oldarray = *_array;
  TacsScalar * newarray = new TacsScalar[ newlen ];
  memcpy(newarray, oldarray, oldlen*sizeof(TacsScalar));
  delete [] *_array;
  *_array = newarray;
}

/*!  
  Given an unsorted CSR data structure, with duplicate entries,
  sort/uniquify each array, and recompute the rowp values on the
  fly.

  For each row, sort the array and remove duplicates.  Copy the
  values if required and skip the diagonal.  Note that the copy may
  be overlapping so memcpy cannot be used.
*/
void SortAndUniquifyCSR( int nvars, int * rowp, 
			 int * cols, int nodiag ){
  // Uniquify each column of the array
  int old_start = 0;
  int new_start = 0;
  for ( int i = 0; i < nvars; i++ ){
    // sort cols[start:rowp[i]] 
    int rsize = FElibrary::uniqueSort(&cols[old_start], rowp[i+1]-old_start);
      
    if (nodiag){
      int end = old_start + rsize;
      for ( int j = old_start, k = new_start; j < end; j++, k++ ){
        if (cols[j] == i){
          rsize--;
          k--;	  
        }
        else if (j != k){
          cols[k] = cols[j];
        }
      }
    }
    else if (old_start != new_start){
      int end = old_start + rsize;
      for ( int j = old_start, k = new_start; j < end; j++, k++ ){
        cols[k] = cols[j];
      }
    }
      
    old_start = rowp[i+1];
    rowp[i] = new_start;
    new_start += rsize;    
  }

  rowp[nvars] = new_start;
}


/*!
  Perform Reverse Cuthill-McKee ordering.

  Input:
  ------
  nvars, rowp, cols == The CSR data structure containing the 
  graph representation of the matrix

  root: The root node to perform the reordering from
  Returns:
  --------
  rcm_vars == The new variable ordering 

  The number of variables ordered in this pass of the RCM reordering
*/
int ComputeRCMOrder( const int nvars, const int * rowp, const int * cols,
                     int * rcm_vars, int root, int n_rcm_iters ){    
  if (n_rcm_iters < 1){
    n_rcm_iters = 1;
  }

  int * levset = new int[ nvars ];
   
  int rvars = 0;
  for ( int k = 0; k < n_rcm_iters; k++ ){
    rvars = ComputeRCMLevSetOrder(nvars, rowp, cols, 
                                  rcm_vars, levset, root);
    if (nvars != rvars){
      return rvars;
    }

    root = rcm_vars[0];
  }

  delete [] levset;

  return rvars;
}

/*
  Compute the RCM level sets andreordering of the graph given by the
  symmetric CSR data structure rowp/cols.

  rowp/cols represents the non-zero structure of the matrix to be
  re-ordered

  Here levset is a unique, 0 to nvars array containing the level
  sets
*/
int ComputeRCMLevSetOrder( const int nvars, const int * rowp, 
                           const int * cols, int * rcm_vars, int * levset, 
                           int root ){
  int start = 0; // The start of the current level
  int end   = 0; // The end of the current level

  // Set all the new variable numbers to -1
  for ( int k = 0; k < nvars; k++ ){
    rcm_vars[k] = -1;
  }

  int var_num = 0;
  while ( var_num < nvars ){
    // If the current level set is empty, find any un-ordered variables
    if ( end-start == 0 ){
      if ( rcm_vars[root] >= 0 ){
        // Find an appropriate root
        int i = 0;
        for ( ; i < nvars; i++ ){
          if ( rcm_vars[i] < 0 ){
            root = i;
            break;
          }
        }
        if ( i >= nvars ){
          return var_num;
        }
      }
    
      levset[end] = root;
      rcm_vars[root] = var_num;
      var_num++;
      end++;
    }

    while ( start < end ){
      int next = end;
      // Iterate over the nodes added to the previous level set
      for ( int current = start; current < end; current++ ){
        int node = levset[current];

        // Add all the nodes in the next level set
        for ( int j = rowp[node]; j < rowp[node+1]; j++ ){
          int next_node = cols[j];
	    
          if ( rcm_vars[next_node] < 0 ){
            rcm_vars[next_node] = var_num;
            levset[next] = next_node;
            var_num++;
            next++;	
          }      
        }
      }    

      start = end;
      end   = next;
    }
  }

  // Go through and reverse the ordering of all the variables
  for ( int i = 0; i < var_num; i++ ){
    int node = levset[i];
    rcm_vars[node] = var_num-1 - i;
  }

  return var_num;
}

/*!  
  The following function takes a list of elements and the
  corresponding variables and creates a list of variables that point
  to each variable.
*/
void ElementListToVarList( const int * elems, const int * elemptr, int nelems,
                           int lower, int upper,
                           int ** _varelems, int ** _varptr ){
  int N = upper - lower;
  int * varptr = new int[ N+1 ];

  for ( int i = 0; i < N+1; i++ ){
    varptr[i] = 0;
  }
  
  for ( int i = 0; i < nelems; i++ ){
    for ( int j = elemptr[i]; j < elemptr[i+1]; j++ ){
      if ( elems[j] >= lower && elems[j] < upper ){
        int var = elems[j] - lower;
        varptr[var+1]++;
      }
    }
  }

  for ( int i = 0; i < N; i++ ){
    varptr[i+1] += varptr[i];
  }

  int * varelems = new int[ varptr[N] ];

  for ( int i = 0; i < nelems; i++ ){
    for ( int j = elemptr[i]; j < elemptr[i+1]; j++ ){
      if ( elems[j] >= lower && elems[j] < upper ){
        int var   = elems[j] - lower;
        varelems[ varptr[var] ] = i;
        varptr[var]++;
      }
    }
  }

  // Put the pointers back to their proper position
  for ( int i = 0; i < N; i++ ){
    varptr[ N-i ] = varptr[ N-i-1 ];
  }
  varptr[0] = 0;

  *_varptr = varptr;
  *_varelems = varelems;
}

/*!  
  Given a list of elements -> variables, an ownership range
  defined by lower <= var < upper and an external variable to local
  map defined by extvars[i], compute the variable to element map
*/
void ElementListToExtVarList( const int * elems, const int * elemptr, int nelems,
                              const int * extvars, int nextvars,
                              int lower, int upper,
                              int ** _varelems, int ** _varptr ){
  // Construct a variable -> element scheme for the off-processor variables  
  int * varptr = new int[ nextvars+1 ];
  for ( int i = 0; i < nextvars+1; i++ ){
    varptr[i] = 0;
  }

  for ( int i = 0; i < nelems; i++ ){
    for ( int j = elemptr[i]; j < elemptr[i+1]; j++ ){
      if ( elems[j] < lower || elems[j] >= upper ){
        // search ext_vars for this variable
        int * item = (int*)bsearch(&elems[j], extvars, nextvars, 
                                   sizeof (int), FElibrary::comparator);

        if ( item != NULL ){
          int var = item - extvars;
          varptr[var+1]++;
        }  
      }
    }
  }

  for ( int i = 0; i < nextvars; i++ ){
    varptr[i+1] += varptr[i];
  }

  int * varelems = new int[ varptr[nextvars] ];
  for ( int i = 0; i < nelems; i++ ){
    for ( int j = elemptr[i]; j < elemptr[i+1]; j++ ){
      if ( elems[j] < lower || elems[j] >= upper ){
        // search the ext_vars for this variable
        int * item = (int*) bsearch(&elems[j], extvars, nextvars, 
                                    sizeof (int), FElibrary::comparator);

        if ( item != NULL ){
          int var = item - extvars; // index into ext_vars corresponding to elems[j]
          varelems[ varptr[var] ] = i;
          varptr[var]++;
        }
      }
    }
  }

  // Put the pointers back to their proper position
  for ( int i = 0; i < nextvars; i++ ){
    varptr[ nextvars-i ] = varptr[ nextvars-i-1 ];
  }
  varptr[0] = 0;

  *_varptr = varptr;
  *_varelems = varelems;
}

TACS_END_NAMESPACE
