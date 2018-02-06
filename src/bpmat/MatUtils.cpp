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

#include "FElibrary.h"
#include "MatUtils.h"

TACS_BEGIN_NAMESPACE(matutils)

/*!
  Extend the length of an integer array to a new length of array
*/
void ExtendArray( int **_array, int oldlen, int newlen ){
  int *oldarray = *_array;
  int *newarray = new int[ newlen ];
  memcpy(newarray, oldarray, oldlen*sizeof(int));
  delete [] *_array;
  *_array = newarray;
}

/*
  Extend the length of a TacsScalar array to a new length
*/
void ExtendArray( TacsScalar **_array, int oldlen, int newlen ){
  TacsScalar *oldarray = *_array;
  TacsScalar *newarray = new TacsScalar[ newlen ];
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
void SortAndUniquifyCSR( int nvars, int *rowp, 
                         int *cols, int remove_diag ){
  // Uniquify each column of the array
  int old_start = 0;
  int new_start = 0;
  for ( int i = 0; i < nvars; i++ ){
    // sort cols[start:rowp[i]] 
    int rsize = FElibrary::uniqueSort(&cols[old_start], rowp[i+1]-old_start);
      
    if (remove_diag){
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
int ComputeRCMOrder( const int nvars, const int *rowp, const int *cols,
                     int *rcm_vars, int root, int n_rcm_iters ){    
  if (n_rcm_iters < 1){
    n_rcm_iters = 1;
  }

  int *levset = new int[ nvars ];
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
int ComputeRCMLevSetOrder( const int nvars, const int *rowp, 
                           const int *cols, int *rcm_vars, int *levset, 
                           int root ){
  int start = 0; // The start of the current level
  int end   = 0; // The end of the current level

  // Set all the new variable numbers to -1
  for ( int k = 0; k < nvars; k++ ){
    rcm_vars[k] = -1;
  }

  int var_num = 0;
  while (var_num < nvars){
    // If the current level set is empty, find any un-ordered variables
    if (end-start == 0){
      if (rcm_vars[root] >= 0){
        // Find an appropriate root
        int i = 0;
        for ( ; i < nvars; i++ ){
          if (rcm_vars[i] < 0){
            root = i;
            break;
          }
        }
        if (i >= nvars){
          return var_num;
        }
      }
    
      levset[end] = root;
      rcm_vars[root] = var_num;
      var_num++;
      end++;
    }

    while (start < end){
      int next = end;
      // Iterate over the nodes added to the previous level set
      for ( int current = start; current < end; current++ ){
        int node = levset[current];

        // Add all the nodes in the next level set
        for ( int j = rowp[node]; j < rowp[node+1]; j++ ){
          int next_node = cols[j];
            
          if (rcm_vars[next_node] < 0){
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

/*
  Multicolor code for a single process using a greedy algorithm
*/
int ComputeSerialMultiColor( const int nvars, const int *rowp,
                             const int *cols, int *colors, 
                             int *new_vars ){
  // Allocate a temporary array to store the
  int *tmp = new int[ nvars+1 ];
  for ( int i = 0; i < nvars; i++ ){
    tmp[i] = -1;
    colors[i] = -1;
  }
  
  int num_colors = 0;
  for ( int i = 0; i < nvars; i++ ){
    // Find the minimum color that is not referred to by any adjacent
    // node. 
    for ( int jp = rowp[i]; jp < rowp[i+1]; jp++ ){
      int j = cols[jp];
      if (colors[j] >= 0){
        tmp[colors[j]] = i;
      }
    }

    // Set the color for this variable if it already exists
    int flag = 1;
    for ( int k = 0; k < num_colors; k++ ){
      if (tmp[k] != i){
        colors[i] = k;
        flag = 0;
        break;
      }
    }

    // Create a new color
    if (flag){
      colors[i] = num_colors;
      num_colors++;
    }
  }

  // Now that all the nodes have been colored, order them
  memset(tmp, 0, (num_colors+1)*sizeof(int));

  // Count up the number of nodes for each color
  for ( int i = 0; i < nvars; i++ ){
    tmp[colors[i]+1]++;
  }

  // Set tmp as an offset for each color
  for ( int i = 1; i < num_colors+1; i++ ){
    tmp[i] += tmp[i-1];
  }
  
  // Create the new color variables
  for ( int i = 0; i < nvars; i++ ){
    new_vars[i] = tmp[colors[i]];
    tmp[colors[i]]++;
  }

  delete [] tmp;

  return num_colors;
}

TACS_END_NAMESPACE
