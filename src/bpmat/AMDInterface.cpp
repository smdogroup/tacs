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

#include "AMDInterface.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
  Compare two integers. This is used when sorting an array.
*/
static int compare_ints(const void *a, const void *b) {
  return (*(int *)a - *(int *)b);
}

/*
  Sort an array of length len, then remove duplicate entries and
  entries with values -1.
*/
static int amd_remove_duplicates(int *array, int len) {
  qsort(array, len, sizeof(int), compare_ints);

  // Remove any negative numbers
  int i = 0;  // location to take entries from
  int j = 0;  // location to place entries

  while (i < len && array[i] < 0) i++;

  for (; i < len; i++, j++) {
    while ((i < len - 1) && (array[i] == array[i + 1])) i++;

    if (i != j) {
      array[j] = array[i];
    }
  }

  return j;  // The new length of the array
}

/*
  Check the formatting of the AMD data structure to ensure that things
  are still ordered correctly.

  If there is a problem with one of the rows of the data structure,
  return row+1, otherwise, return 0;
*/
static int amd_check_format(int nvars, int *rowp, int *cols, int *elen,
                            int *alen) {
  int flag = 0;
  for (int i = 0; i < nvars; i++) {
    for (int j = rowp[i]; j < rowp[i] + elen[i] - 1; j++) {
      if (cols[j] < 0 || cols[j + 1] < 0 || cols[j + 1] <= cols[j]) {
        flag = i + 1;
        break;
      }
    }
    if (flag) {
      break;
    }

    for (int j = rowp[i] + elen[i]; j < rowp[i] + alen[i] - 1; j++) {
      if (cols[j] < 0 || cols[j + 1] < 0 || cols[j + 1] <= cols[j]) {
        flag = i + 1;
        break;
      }
    }
    if (flag) {
      break;
    }
  }

  return flag;
}

/*
  Find 'required_space' extra locations for the row entry by moving
  the entries in the array to different positions.

  Input:
  nvars, rowp, cols, alen: The current entries in the AMD data structure

  Procedure:

  First scan through the rows starting from r, and proceeding to 0,
  counting up the space that would be freed by compressing them.
  Next, copy those rows, and modifying rowp/cols such that the bound
  alen is tight - note that elen/alen values are not modified.  If
  that isn't enough space, count up the space from r, to nvars-1.
  Compress the required number of rows to free up the required amount
  of space.
*/
static void amd_add_space(int nvars, int *rowp, int *cols, const int *alen,
                          const int r, int required_space) {
  // First, try and collect the space required from the rows preceeding r
  int new_space = 0;

  if (r > 0) {
    int start = r - 1;
    // Count up the space freed by compressing rows start through r
    for (; new_space < required_space && start >= 0; start--) {
      new_space += (rowp[start + 1] - rowp[start]) - alen[start];
    }
    start++;

    int j = rowp[start] + alen[start];  // Where new entries will be placed
    for (int i = start + 1; i <= r; i++) {
      int k = rowp[i];
      int new_rowp = j;
      int end = k + alen[i];
      for (; k < end; k++, j++) {
        cols[j] = cols[k];
      }
      rowp[i] = new_rowp;
    }
  }

  // If not enough space has been found, use the remaining columns
  if (new_space < required_space) {
    int start = r + 1;
    for (; new_space < required_space && start < nvars; start++) {
      new_space += (rowp[start + 1] - rowp[start]) - alen[start];
    }
    start--;

    // Cannot exceed the size of the matrix - print an error here
    if (start >= nvars) {
      start = nvars - 1;
      fprintf(stderr, "Error, not enough memory found\n");
    }

    // Start from the end of the entries
    int j = rowp[start + 1] - 1;  // Where new entries will be placed
    for (int i = start; i > r; i--) {
      int end = rowp[i];
      int k = end + alen[i] - 1;
      for (; k >= end; k--, j--) {
        cols[j] = cols[k];
      }
      rowp[i] = rowp[i + 1] - alen[i];  // Tight bound on the entries for i
    }
  }
}

/*
  Compare two variables to see if they are indistinguishable

  This checks to see if the variables i and j form the same adjacency
  structure. If so, then they will be used to form a supervariable.
  This function first checks if the nodes are the same, and then
  checks that the { adj[j], i } is equal to { adj[i], j }. This makes
  things a bit more interesting.

  Input:
  i, j: the nodes to be compared

  elen, alen: The length of the elements, and total number variables
  and nodes respectively

  rowp, cols: The quotient graph data structure.
*/
static int amd_compare_variables(const int i, const int j, const int *elen,
                                 const int *alen, const int *rowp,
                                 const int *cols) {
  // First, check if they are the same length
  if (i == j) {
    return 0;  // The same node, this should be avoided
  } else if ((alen[i] != alen[j]) ||
             (elen[i] != elen[j])) {  // The lengths must be the same
    return 0;
  } else {
    // Compare cols[rowp[i] ... rowp[i] + alen[i]-1] and
    //         cols[rowp[j] ... rowp[j] + alen[j]-1]
    int size = alen[i];
    int ip = rowp[i], iend = ip + alen[i];
    int jp = rowp[j], jend = jp + alen[j];

    for (int k = 0; k < size; k++, jp++, ip++) {
      if (cols[ip] == j) {
        ip++;
      }
      if (cols[jp] == i) {
        jp++;
      }
      if (cols[ip] != cols[jp]) {
        break;
      }
    }

    return (ip == iend && jp == jend);
  }
}

/*
  Remove the references to a variable from the data structure.

  This is used when removing an indistinguishable variable.
*/
void amd_remove_variable(int var, int *elen, int *alen, int *rowp, int *cols,
                         int nvars) {
  int i = rowp[var];

  // First, visit all the elements pointed to by var
  for (; i < rowp[var] + elen[var]; i++) {
    int e = cols[i];  // The element number

    int j = rowp[e];
    int jend = (rowp[e] + alen[e]) - 1;
    for (; j < jend; j++) {
      if (cols[j] == var) {
        break;
      }
    }

    // cols[j] == var: This should always be true
    for (; j < jend; j++) {
      cols[j] = cols[j + 1];
    }
    alen[e]--;
  }

  // Remove the variable from the reference
  for (; i < rowp[var] + alen[var]; i++) {
    int v = cols[i];

    int j = rowp[v] + elen[v];
    int jend = (rowp[v] + alen[v]) - 1;
    for (; j < jend; j++) {
      if (cols[j] == var) {
        break;
      }
    }

    // cols[j] == var: This should always be true
    for (; j < jend; j++) {
      cols[j] = cols[j + 1];
    }
    alen[v]--;
  }

  elen[var] = 0;
  alen[var] = 0;
}

/*
  The following function performs an approximate minimum degree
  reordering of the input matrix.  This code follows the description
  of (Amestoy, Davis and Duff, An Approximate Minimum Degree
  Reordering Algorithm, SIAM J. Mat. Anl. 1996). The advantage of
  the method implemented below, is that the ordering of the
  interface unknowns, used in the global Schur complement system,
  can be delayed until all other variables are ordered. This allows
  a the full Schur method to be used efficiently within TACS.

  This uses a quotient graph representation of the elimination graph.
  The quotient graph representation uses a known amount of memory -
  less than or equal to the storage of A itself. The implementation
  below uses the arrays provided and as a result, the contents of
  these arrays are destroyed during the call.

  Input data description:
  -----------------------

  nvars: The dimensions of the square input matrix

  rowp: An integer that points to the beginnning entry for each row in
  cols

  cols: An array of the entries in A. During the computation this is
  modified. The output of the matrix is meaningless. All entries are
  destroyed.

  Output:
  -------

  perm: The permutation array for the new ordering.

  Data allocated within the function:
  -----------------------------------

  Lp: A temporary array of size nvars used to collect the variables
  associated pivot

  degrees: An integer array storing the current degree estimates

  elem_degrees: An integer array used during the aprpoximation of
  the variable degrees - represents |Le \ Lp|

  state: An integer array denoting the current status of a variable.
  The states are: interface variable (cannot be selected as a
  pivot), variable (not yet eliminated), element (an eliminated
  variable) and inactive element (no longer required by the algorithm)

  elen: The number of elements in the row. The elements are stored
  first within the row in ascending order.

  alen: The number of elements + variables in the row. The variables
  are stored immediately after the elements.

  slen: The length of the supernode!

  hash: The hash table values used to compare nodes to see if they
  are indistinguishable

  As a result:
  cols[rowp[i] ... rowp[i] + elen[i]] == The elements in the row

  cols[rowp[i] + elen[i] ... rowp[i] + alen[i]] == The
  variables in each row

  On initialization elen[:] = 0, alen[i] = rowp[i+1] - rowp[i]
*/
void amd_order_interface(int nvars, int *rowp, int *cols, int *perm,
                         int *interface_nodes, int ninterface_nodes,
                         int ndep_vars, const int *dep_vars,
                         const int *indep_ptr, const int *indep_vars,
                         int use_exact_degree) {
  int *alen = new int[nvars];  // Number of entries in a row
  int *elen = new int[nvars];  // Number of elements in a row
  int *slen = new int[nvars];  // The length of each supernode

  int *degree = new int[nvars];       // Degree (exact or estimate) of the node
  int *elem_degree = new int[nvars];  // The degree of the element
  int *Lp = new int[nvars];           // A set of variables
  int *state = new int[nvars];        // Are we a variable, or element?

  // Allocate more space for the dependent vars - if any
  int *dep_count = NULL;
  int *vars_to_dep_ptr = NULL;
  int *vars_to_dep = NULL;

  // Build the information required for the dependent-nodes.
  // All independent nodes associated with a dependent node must
  // be ordered before the dependent.
  if (ndep_vars > 0) {
    // Count up number of times that the indepdnent variables are
    // referenced by each dependent variable
    dep_count = new int[nvars];
    memset(dep_count, 0, nvars * sizeof(int));
    for (int i = 0; i < ndep_vars; i++) {
      dep_count[dep_vars[i]] += indep_ptr[i + 1] - indep_ptr[i];
    }

    // Set up the data for the independent to dependent counter
    vars_to_dep_ptr = new int[nvars + 1];
    memset(vars_to_dep_ptr, 0, (nvars + 1) * sizeof(int));
    for (int i = 0; i < indep_ptr[ndep_vars]; i++) {
      vars_to_dep_ptr[indep_vars[i] + 1]++;
    }

    // Set up the pointer so that it offsets to the dep_pointer
    for (int i = 1; i <= nvars; i++) {
      vars_to_dep_ptr[i] += vars_to_dep_ptr[i - 1];
    }

    // Allocate the variable array
    vars_to_dep = new int[vars_to_dep_ptr[nvars]];
    for (int i = 0; i < ndep_vars; i++) {
      for (int j = indep_ptr[i]; j < indep_ptr[i + 1]; j++) {
        int index = indep_vars[j];

        // Set the dependent variable
        vars_to_dep[vars_to_dep_ptr[index]] = dep_vars[i];
        vars_to_dep_ptr[index]++;
      }
    }

    // Reset the vars_to_dep_ptr array
    int offset = 0;
    for (int i = 0; i <= nvars; i++) {
      int tmp = vars_to_dep_ptr[i];
      vars_to_dep_ptr[i] = offset;
      offset = tmp;
    }
  }

  for (int i = 0; i < nvars; i++) {
    perm[i] = i;
    degree[i] = rowp[i + 1] - rowp[i];
    alen[i] = rowp[i + 1] - rowp[i];
    elen[i] = 0;
    state[i] = 1;  // Set everything to variables
    slen[i] = 1;
  }

  // Sort the interface nodes
  if (ninterface_nodes !=
      amd_remove_duplicates(interface_nodes, ninterface_nodes)) {
    fprintf(stderr, "amd_order_interace: Interface nodes must be unique\n");
  }

  // Set the interface node states
  for (int i = 0; i < ninterface_nodes; i++) {
    int node = interface_nodes[i];
    state[node] = 2;
  }

  if (amd_check_format(nvars, rowp, cols, elen, alen) != 0) {
    fprintf(stderr,
            "amd_order_interface: Detected problem with \
initial data structure\n");
    return;
  }

  // Perform the elimination
  int nsvars = 0;  // Keep track of the number of supervariables
  for (int i = 0; i < nvars; nsvars++) {
    // Select the pivots
    int piv = -1;
    int min_degree = nvars + 1;

    if (dep_count) {
      // Find the minimum degree with dep_count[j] == 0
      if (i < nvars - ninterface_nodes) {
        // Find the minimum degree variable
        for (int j = 0; j < nvars; j++) {
          if (state[j] == 1 && dep_count[j] == 0 && degree[j] < min_degree) {
            piv = j;
            min_degree = degree[j];
          }
        }
      } else {
        // Find the minimum interface variable
        for (int j = 0; j < nvars; j++) {
          if (state[j] == 2 && dep_count[j] == 0 && degree[j] < min_degree) {
            piv = j;
            min_degree = degree[j];
          }
        }
      }

      if (piv < 0) {
        fprintf(stderr,
                "amd_order_interface: Negative pivot encountered at %d\n", i);
        fflush(stderr);
      }

      perm[nsvars] = piv;
      state[piv] = 0;  // Eliminated variable

      // Update the counts
      for (int j = vars_to_dep_ptr[piv]; j < vars_to_dep_ptr[piv + 1]; j++) {
        dep_count[vars_to_dep[j]]--;
      }
    } else {
      if (i < nvars - ninterface_nodes) {
        // Find the minimum degree variable
        for (int j = 0; j < nvars; j++) {
          if (state[j] == 1 && degree[j] < min_degree) {
            piv = j;
            min_degree = degree[j];
          }
        }
      } else {
        // Find the minimum interface variable
        for (int j = 0; j < nvars; j++) {
          if (state[j] == 2 && degree[j] < min_degree) {
            piv = j;
            min_degree = degree[j];
          }
        }
      }

      perm[nsvars] = piv;
      state[piv] = 0;  // Eliminated variable
    }

    // Determine the set of variables Lp U { piv } - the non-zeros
    // in the column of L
    int lenlp = 0;
    Lp[lenlp] = piv;
    lenlp++;

    // Add the contributions from the row piv
    for (int j = rowp[piv] + elen[piv]; j < rowp[piv] + alen[piv]; j++) {
      if (piv != cols[j]) {
        Lp[lenlp] = cols[j];
        lenlp++;
      }
    }

    // Add the non-zero pattern to Lp
    // First add the contributions from the elements
    for (int j = rowp[piv]; j < rowp[piv] + elen[piv]; j++) {
      int e = cols[j];
      for (int k = rowp[e] + elen[e]; k < rowp[e] + alen[e]; k++) {
        if (lenlp >= nvars) {
          // Remove duplicates - this has to free up enough space for Lp
          lenlp = amd_remove_duplicates(Lp, lenlp);
        }
        if (cols[k] != piv) {
          Lp[lenlp] = cols[k];
          lenlp++;
        }
      }
    }

    // Sort and remove any duplicates from the list Lp
    lenlp = amd_remove_duplicates(Lp, lenlp);

    // Update the non-zero pattern in cols
    for (int j = 0; j < lenlp; j++) {
      int lj = Lp[j];

      if (lj == piv) {
        continue;
      }

      // Absorb elements that are no longer required
      int nre = 0;  // Number of removed elements
      int end = rowp[piv] + elen[piv];
      for (int k = rowp[piv], ck = rowp[lj];
           ((k < end) && (ck < rowp[lj] + elen[lj])); ck++) {
        while ((k < end) && (cols[k] < cols[ck])) {
          k++;
        }
        if (k < end && cols[k] == cols[ck]) {
          cols[ck] = -1;  // Remove the element from the list
          nre++;
        }
      }

      // Remove Lp[j], Lp[k] if it exists in cols
      int nrv = 0;  // Number of removed variables
      // Remove Lp U { piv } from the columns
      for (int k = 0, ck = rowp[lj] + elen[lj];
           (k < lenlp) && (ck < rowp[lj] + alen[lj]); ck++) {
        while ((k < lenlp) && (Lp[k] < cols[ck])) {
          k++;
        }
        if (k < lenlp && cols[ck] == Lp[k]) {
          cols[ck] = -1;  // Remove it from the list
          nrv++;
        }
      }

      // Remove negative entries from the element list
      if (nre > 0 || nrv > 0) {
        int end = rowp[lj] + alen[lj];
        int k = rowp[lj];
        int nk = k;
        for (; k < end; k++, nk++) {
          while (k < end && cols[k] == -1) {
            k++;
          }
          if (k < end) {
            cols[nk] = cols[k];
          }
        }

        // Adjust the number of variables and elements within the list
        elen[lj] = elen[lj] - nre;
        alen[lj] = alen[lj] - nre - nrv;
      }
    }

    // Now, add piv to the elements in rows Lp \ {piv}
    for (int j = 0; j < lenlp; j++) {
      int lj = Lp[j];

      if (lj == piv) {
        continue;
      }

      if (alen[lj] == rowp[lj + 1] - rowp[lj]) {
        amd_add_space(nvars, rowp, cols, alen, lj, 1);
      }

      // Now, find where piv should go in [rowp[lj] ... rowp[lj] + elen[lj]-1]
      int k = rowp[lj];
      int end = k + elen[lj];
      while (k < end && cols[k] < piv) {
        k++;
      }

      if (cols[k] != piv) {
        int p = rowp[lj] + alen[lj] - 1;
        while (p >= k) {
          cols[p + 1] = cols[p];
          p--;
        }
        cols[k] = piv;

        alen[lj] += 1;
        elen[lj] += 1;
      }
    }

    // Remove the rows associated with the elements of piv
    for (int j = rowp[piv]; j < rowp[piv] + elen[piv]; j++) {
      int e = cols[j];
      // Remove row e
      alen[e] = 0;
      elen[e] = 0;
      state[e] = -1;  // This element has been entirely removed
    }

    // Copy Lp to the row piv
    int rsize = rowp[piv + 1] - rowp[piv];

    // Test to see if new space is requried
    if (lenlp - 1 > rsize) {
      amd_add_space(nvars, rowp, cols, alen, piv,
                    lenlp - 1 - (rowp[piv + 1] - rowp[piv]));
    }

    elen[piv] = 0;
    alen[piv] = lenlp - 1;
    for (int j = rowp[piv], k = 0; k < lenlp; k++, j++) {
      if (Lp[k] == piv) {
        k++;
        if (k == lenlp) {
          break;
        }
      }
      cols[j] = Lp[k];
    }

    if (use_exact_degree) {
      // Update the degrees for each variable in the row piv
      for (int j = rowp[piv]; j < rowp[piv] + alen[piv]; j++) {
        // This row should be entirely variables by definition
        if (state[cols[j]] <= 0) {
          fprintf(stderr,
                  "Error, the pivot row %d should contain only \
variables, not element %d\n",
                  piv, cols[j]);
        }

        lenlp = 0;
        int var = cols[j];

        // Add the contributions from A
        for (int k = rowp[var] + elen[var]; k < rowp[var] + alen[var]; k++) {
          Lp[lenlp] = cols[k];
          lenlp++;
        }

        // Add the sets from the elements in this row
        for (int k = rowp[var]; k < rowp[var] + elen[var]; k++) {
          int e = cols[k];  // Take the elements belonging to var

          // Now, add in the variables corresponding to the element e
          for (int p = rowp[e] + elen[e]; p < rowp[e] + alen[e]; p++) {
            if (cols[p] != var) {
              if (lenlp >= nvars) {
                lenlp = amd_remove_duplicates(Lp, lenlp);
              }
              Lp[lenlp] = cols[p];
              lenlp++;
            }
          }
        }

        lenlp = amd_remove_duplicates(Lp, lenlp);

        int deg = 0;
        for (int k = 0; k < lenlp; k++) {
          deg += slen[Lp[k]];
        }
        degree[var] = deg;
      }
    } else {  // Approximate degree
      // The worst cases are:
      // The trailing submatrix is dense: degree = N-i
      // All nodes in k, result in new fill in: d^{i+1} = d^{i} + lenlp
      // The approximate degree is:
      //          |A_{i} \ i| + |L_p \ i| + sum_{e} | L_{e} \ L_p|

      // Determine the degrees of the un-eliminated elements
      for (int j = 0; j < nvars; j++) {
        elem_degree[j] = -1;
      }

      // For each supervariable in Lp
      for (int j = 0; j < lenlp; j++) {
        int lj = Lp[j];
        int end = rowp[lj] + elen[lj];
        // For all elements pointed to by row Lp[j] = lj
        for (int k = rowp[lj]; k < end; k++) {
          // Find all the elements
          int e = cols[k];
          if (elem_degree[e] < 0) {
            // Calculate the element degree
            elem_degree[e] = 0;
            for (int p = rowp[e]; p < rowp[e] + alen[e]; p++) {
              elem_degree[e] += slen[cols[p]];
            }
          }
          elem_degree[e] -= slen[lj];
        }
      }

      // Compute | Lp |
      int deg_Lp = 0;
      for (int j = 0; j < lenlp; j++) {
        if (Lp[j] != piv) {
          deg_Lp += slen[Lp[j]];
        }
      }

      for (int j = 0; j < lenlp; j++) {
        int lj = Lp[j];

        if (lj == piv) {
          continue;
        }

        // Add the contributions from A
        int deg_estimate = deg_Lp;
        for (int k = rowp[lj] + elen[lj]; k < rowp[lj] + alen[lj]; k++) {
          deg_estimate += slen[cols[k]];
        }

        // If lj is in A
        int start = rowp[lj] + elen[lj];
        int len = alen[lj] - elen[lj];
        if (len > 0) {
          if (bsearch(&lj, &cols[start], len, sizeof(int), compare_ints)) {
            deg_estimate -= slen[lj];
          }
        }

        // If lj is in Lp U piv
        if (bsearch(&lj, Lp, lenlp, sizeof(int), compare_ints)) {
          deg_estimate -= slen[lj];
        }

        for (int k = rowp[lj]; k < rowp[lj] + elen[lj]; k++) {
          int e = cols[k];
          if (e == piv) {
            continue;
          }

          if (elem_degree[e] < 0) {
            elem_degree[e] = 0;
            for (int p = rowp[e]; p < rowp[e] + alen[e]; p++) {
              elem_degree[e] += slen[cols[p]];
            }
            deg_estimate += elem_degree[e];
          } else {
            deg_estimate += elem_degree[e];
          }
        }

        // Now, compute the degree estimate
        deg_estimate = (deg_estimate < nvars - i ? deg_estimate : nvars - i);
        degree[lj] = (deg_estimate < degree[lj] + deg_Lp ? deg_estimate
                                                         : degree[lj] + deg_Lp);
      }
    }

    // Supervariable detection and construction
    if (i < nvars - ninterface_nodes) {
      // If not all interior nodes have been found, only create
      // supervariables from interior nodes
      for (int j = 0; j < lenlp; j++) {
        int lj = Lp[j];

        // Check if this is a dependent variable that has not
        // been ordered yet
        int dep_check = 0;
        if (dep_count) {
          dep_check = (dep_count[lj] > 0);
        }
        if (dep_check || lj == piv || slen[lj] < 0 ||
            bsearch(&lj, interface_nodes, ninterface_nodes, sizeof(int),
                    compare_ints)) {
          continue;
        }

        for (int k = j + 1; k < lenlp; k++) {
          int lk = Lp[k];

          // Check if this is a dependent variable
          dep_check = 0;
          if (dep_count) {
            dep_check = (dep_count[lk] > 0);
          }
          if (dep_check || lk == piv || slen[lk] < 0 ||
              bsearch(&lk, interface_nodes, ninterface_nodes, sizeof(int),
                      compare_ints)) {
            continue;
          }

          // Quick check to see if the nodes are the same
          if (amd_compare_variables(lj, lk, elen, alen, rowp, cols)) {
            // Update the dependent node counts
            if (dep_count && slen[lk] == 1) {
              for (int p = vars_to_dep_ptr[lk]; p < vars_to_dep_ptr[lk + 1];
                   p++) {
                dep_count[vars_to_dep[p]]--;
              }
            }

            // Merge lk into lj
            slen[lj] += slen[lk];
            degree[lj] -= slen[lk];
            slen[lk] = -(lj + 1);
            state[lk] = -1;

            // Remove lj from the quotient graph
            amd_remove_variable(lk, elen, alen, rowp, cols, nvars);
          }
        }
      }
    } else {
      // All interface nodes have been detected.. go nuts!
      for (int j = 0; j < lenlp; j++) {
        int lj = Lp[j];

        // Check if this is a dependent variable
        int dep_check = 0;
        if (dep_count) {
          dep_check = (dep_count[lj] > 0);
        }
        if (dep_check || lj == piv || slen[lj] < 0) {
          continue;
        }

        for (int k = j + 1; k < lenlp; k++) {
          int lk = Lp[k];
          // Check if this is a dependent variable
          int dep_check = 0;
          if (dep_count) {
            dep_check = (dep_count[lk] > 0);
          }
          if (dep_check || lk == piv || slen[lk] < 0) {
            continue;
          }

          // Quick check to see if the nodes are the same
          if (amd_compare_variables(lj, lk, elen, alen, rowp, cols)) {
            // Update the dependent node counts
            if (dep_count && slen[lk] == 1) {
              for (int p = vars_to_dep_ptr[lk]; p < vars_to_dep_ptr[lk + 1];
                   p++) {
                dep_count[vars_to_dep[p]]--;
              }
            }

            // Merge lk into lj
            slen[lj] += slen[lk];
            degree[lj] -= slen[lk];
            slen[lk] = -(lj + 1);
            state[lk] = -1;

            // Remove lj from the quotient graph
            amd_remove_variable(lk, elen, alen, rowp, cols, nvars);
          }
        }
      }
    }

    // Reduce the number of variables to eliminate
    i = i + slen[piv];
  }

  // We have the following situation...
  // perm[0:nsvars] = piv, contains the supernodal pivots
  // slen[i] > 0 means i is a principal supervariable
  // slen[i] < 0 means variable i was collapsed into supervariable (slen[i]-1)

  // First arrange the principal supervariables in the permutation array
  int end = nvars;
  for (int j = nsvars - 1; j >= 0; j--) {
    int sv = perm[j];     // Get the supervariable selected last
    end -= slen[sv];      // Figure out where it should start
    perm[end] = sv;       // Place the principal variable at the beginning
    state[sv] = end + 1;  // Increment the pointer into the array
  }

  // Fill in the non-principal variables
  for (int k = 0; k < nvars; k++) {
    if (slen[k] < 0) {  // non-principal variable
      // Back-track until we find a princiapl supervariable
      int j = k;
      while (slen[j] < 0) {  //
        j = -(slen[j] + 1);  // j was eliminated by -(slen[j]+1) and so on...
      }
      // j should now be a principal supervariable
      perm[state[j]] = k;
      state[j]++;
    }
  }

  delete[] alen;
  delete[] elen;
  delete[] degree;
  delete[] elem_degree;
  delete[] Lp;
  delete[] state;
  delete[] slen;

  // If the dependent variables are allocated, delete them now
  if (ndep_vars > 0) {
    delete[] dep_count;
    delete[] vars_to_dep_ptr;
    delete[] vars_to_dep;
  }
}
