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

#include "TacsUtilities.h"

#include <stdint.h>
#include <stdlib.h>

int TacsIntegerComparator(const void *a, const void *b) {
  return (*(int *)a - *(int *)b);
}

int TacsSort(int len, int *array) {
  qsort(array, len, sizeof(int), TacsIntegerComparator);
  return len;
}

/*
  Data and comparison function for the argsort
*/
static const TacsScalar *tacs_arg_sort_list = NULL;

static int TacsCompareArgSort(const void *a, const void *b) {
  if (TacsRealPart(tacs_arg_sort_list[*(int *)a]) <
      TacsRealPart(tacs_arg_sort_list[*(int *)b])) {
    return -1;
  } else if (TacsRealPart(tacs_arg_sort_list[*(int *)a]) >
             TacsRealPart(tacs_arg_sort_list[*(int *)b])) {
    return 1;
  }
  return 0;
}

/*
  Sort the values by argument
*/
int TacsArgSort(int len, const TacsScalar *values, int *array) {
  for (int i = 0; i < len; i++) {
    array[i] = i;
  }
  tacs_arg_sort_list = values;
  qsort(array, len, sizeof(int), TacsCompareArgSort);
  tacs_arg_sort_list = NULL;
  return len;
}

/*
  Sort an array and remove duplicate entries from the array. Negative
  values are removed from the array.

  This function is useful when trying to determine the number of
  unique design variable from a set of unsorted variable numbers.

  The algorithm proceeds by first sorting the array, then scaning from
  the beginning to the end of the array, copying values back only once
  duplicate entries are discovered.

  input:
  len:     the length of the array
  array:   the array of values to be sorted

  returns:
  the size of the unique list <= len
*/
int TacsUniqueSort(int len, int *array) {
  // Sort the array
  qsort(array, len, sizeof(int), TacsIntegerComparator);

  int i = 0;  // The location from which to take the entires
  int j = 0;  // The location to place the entries

  // Remove the negative entries
  while (i < len && array[i] < 0) i++;

  for (; i < len; i++, j++) {
    while ((i < len - 1) && (array[i] == array[i + 1])) {
      i++;
    }

    if (i != j) {
      array[j] = array[i];
    }
  }

  return j;
}

/*!
  Merge two sorted arrays into a single sorted array, in place.
  This relies on both arrays having unique elements independently (ie
  a cannot contain duplicates and b cannot contain duplicates, but a
  can contain some of the same elements as b).

  Two part algorithm:
  1. Find the number of duplicates between a and b
  2. Run through the list backwards placing elements into a[]
  when appropriate.
*/
int TacsMergeSortedArrays(int na, int *a, int nb, const int *b) {
  int ndup = 0;

  int j = 0, i = 0;
  for (; i < na; i++) {
    while ((j < nb) && b[j] < a[i]) {
      j++;
    }
    if (j >= nb) {
      break;
    }
    if (a[i] == b[j]) {
      ndup++;
    }
  }

  int len = na + nb - ndup;  // End of the array
  int end = len - 1;

  j = nb - 1;
  i = na - 1;
  while (i >= 0 && j >= 0) {
    if (a[i] > b[j]) {
      a[end] = a[i];
      end--, i--;
    } else if (b[j] > a[i]) {
      a[end] = b[j];
      end--, j--;
    } else {  // b[j] == a[i]
      a[end] = a[i];
      end--, j--, i--;
    }
  }

  // Only need to copy over remaining elements from b - if any
  while (j >= 0) {
    a[j] = b[j];
    j--;
  }

  return len;
}

/*!
  Find the interval such that the given index satisfies:

  intv[k] <= index < intv[k+1]

  The intervals must be non-decreasing. Note that len is equal to the
  length of the intv array which is one more than the total number of
  intervals.
*/
int TacsFindInterval(int index, int len, const int intv[]) {
  // Check that the basic conditions are satisfied
  if (len <= 0) {
    return 0;
  } else if (index < intv[0]) {
    return -1;
  } else if (index > intv[len - 1]) {
    return len - 1;  // The number of intervals is equal to len-1
  }

  int low = 0;
  int high = len - 1;
  int mid = low + (int)((high - low) / 2);

  // By construction the following condition always applies:
  // intv[low] <= index < intv[high]

  while (low != high) {
    if (index < intv[low + 1]) {  // Check the low interval
      return low;
    } else if (intv[high - 1] <= index) {  // Check the high interval
      return high - 1;
    }

    if (intv[mid] <= index) {
      low = mid;
    } else {
      high = mid;
    }

    mid = low + (int)((high - low) / 2);
  }

  return low;
}

/*!
  Match the intervals in a list of sorted variables.
*/
void TacsMatchIntervals(int mpiSize, const int ownerRange[], int nvars,
                        const int vars[], int ext_ptr[]) {
  // ext_ptr[n] should be the greatest integer such that
  // vars[ext_ptr[n]] <= ownerRange[n]
  // all variables on [exp_ptr[m],ext_ptr[m+1]) belong to
  // processor m

  // Nothing to do
  if (nvars == 0) {
    for (int n = 0; n < mpiSize + 1; n++) {
      ext_ptr[n] = 0;
    }

    return;
  }

  for (int n = 0; n < mpiSize + 1; n++) {
    // First check the lower bound
    if (ownerRange[n] <= vars[0]) {
      ext_ptr[n] = 0;
    } else if (ownerRange[n] > vars[nvars - 1]) {  // No more variables
      ext_ptr[n] = nvars;
    } else {  // Determine the interval using a binary search
      int low = 0;
      int high = nvars - 1;
      int mid = low + (int)((high - low) / 2);

      // maintain that the variable is in the interval (vars[low],vars[high])
      // note that if high-low=1, then mid = high
      while (high != mid) {
        if (vars[mid] == ownerRange[n]) {
          break;
        }

        if (ownerRange[n] < vars[mid]) {
          high = mid;
        } else {
          low = mid;
        }

        mid = high - (int)((high - low) / 2);
      }

      ext_ptr[n] = mid;
    }
  }
}

/*!
  Extend the length of an integer array to a new length of array
*/
void TacsExtendArray(int **_array, int oldlen, int newlen) {
  int *oldarray = *_array;
  int *newarray = new int[newlen];
  memcpy(newarray, oldarray, oldlen * sizeof(int));
  delete[] * _array;
  *_array = newarray;
}

/*
  Extend the length of a TacsScalar array to a new length
*/
void TacsExtendArray(TacsScalar **_array, int oldlen, int newlen) {
  TacsScalar *oldarray = *_array;
  TacsScalar *newarray = new TacsScalar[newlen];
  memcpy(newarray, oldarray, oldlen * sizeof(TacsScalar));
  delete[] * _array;
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
void TacsSortAndUniquifyCSR(int nvars, int *rowp, int *cols, int remove_diag) {
  // Uniquify each column of the array
  int old_start = 0;
  int new_start = 0;
  for (int i = 0; i < nvars; i++) {
    // sort cols[start:rowp[i]]
    int rsize = TacsUniqueSort(rowp[i + 1] - old_start, &cols[old_start]);

    if (remove_diag) {
      int end = old_start + rsize;
      for (int j = old_start, k = new_start; j < end; j++, k++) {
        if (cols[j] == i) {
          rsize--;
          k--;
        } else if (j != k) {
          cols[k] = cols[j];
        }
      }
    } else if (old_start != new_start) {
      int end = old_start + rsize;
      for (int j = old_start, k = new_start; j < end; j++, k++) {
        cols[k] = cols[j];
      }
    }

    old_start = rowp[i + 1];
    rowp[i] = new_start;
    new_start += rsize;
  }

  rowp[nvars] = new_start;
}

/*
  Compute the RCM level sets andreordering of the graph given by the
  symmetric CSR data structure rowp/cols.

  rowp/cols represents the non-zero structure of the matrix to be
  re-ordered

  Here levset is a unique, 0 to nvars array containing the level
  sets
*/
static int TacsComputeRCMLevSetOrder(const int nvars, const int *rowp,
                                     const int *cols, int *rcm_vars,
                                     int *levset, int root) {
  int start = 0;  // The start of the current level
  int end = 0;    // The end of the current level

  // Set all the new variable numbers to -1
  for (int k = 0; k < nvars; k++) {
    rcm_vars[k] = -1;
  }

  int var_num = 0;
  while (var_num < nvars) {
    // If the current level set is empty, find any un-ordered variables
    if (end - start == 0) {
      if (rcm_vars[root] >= 0) {
        // Find an appropriate root
        int i = 0;
        for (; i < nvars; i++) {
          if (rcm_vars[i] < 0) {
            root = i;
            break;
          }
        }
        if (i >= nvars) {
          return var_num;
        }
      }

      levset[end] = root;
      rcm_vars[root] = var_num;
      var_num++;
      end++;
    }

    while (start < end) {
      int next = end;
      // Iterate over the nodes added to the previous level set
      for (int current = start; current < end; current++) {
        int node = levset[current];

        // Add all the nodes in the next level set
        for (int j = rowp[node]; j < rowp[node + 1]; j++) {
          int next_node = cols[j];

          if (rcm_vars[next_node] < 0) {
            rcm_vars[next_node] = var_num;
            levset[next] = next_node;
            var_num++;
            next++;
          }
        }
      }

      start = end;
      end = next;
    }
  }

  // Go through and reverse the ordering of all the variables
  for (int i = 0; i < var_num; i++) {
    int node = levset[i];
    rcm_vars[node] = var_num - 1 - i;
  }

  return var_num;
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
int TacsComputeRCMOrder(const int nvars, const int *rowp, const int *cols,
                        int *rcm_vars, int root, int n_rcm_iters) {
  if (n_rcm_iters < 1) {
    n_rcm_iters = 1;
  }

  int *levset = new int[nvars];
  int rvars = 0;
  for (int k = 0; k < n_rcm_iters; k++) {
    rvars =
        TacsComputeRCMLevSetOrder(nvars, rowp, cols, rcm_vars, levset, root);
    if (nvars != rvars) {
      return rvars;
    }
    root = rcm_vars[0];
  }

  delete[] levset;
  return rvars;
}

/*
  Multicolor code for a single process using a greedy algorithm
*/
int TacsComputeSerialMultiColor(const int nvars, const int *rowp,
                                const int *cols, int *colors, int *new_vars) {
  // Allocate a temporary array to store the
  int *tmp = new int[nvars + 1];
  for (int i = 0; i < nvars; i++) {
    tmp[i] = -1;
    colors[i] = -1;
  }

  int num_colors = 0;
  for (int i = 0; i < nvars; i++) {
    // Find the minimum color that is not referred to by any adjacent
    // node.
    for (int jp = rowp[i]; jp < rowp[i + 1]; jp++) {
      int j = cols[jp];
      if (colors[j] >= 0) {
        tmp[colors[j]] = i;
      }
    }

    // Set the color for this variable if it already exists
    int flag = 1;
    for (int k = 0; k < num_colors; k++) {
      if (tmp[k] != i) {
        colors[i] = k;
        flag = 0;
        break;
      }
    }

    // Create a new color
    if (flag) {
      colors[i] = num_colors;
      num_colors++;
    }
  }

  // Now that all the nodes have been colored, order them
  memset(tmp, 0, (num_colors + 1) * sizeof(int));

  // Count up the number of nodes for each color
  for (int i = 0; i < nvars; i++) {
    tmp[colors[i] + 1]++;
  }

  // Set tmp as an offset for each color
  for (int i = 1; i < num_colors + 1; i++) {
    tmp[i] += tmp[i - 1];
  }

  // Create the new color variables
  for (int i = 0; i < nvars; i++) {
    new_vars[i] = tmp[colors[i]];
    tmp[colors[i]]++;
  }

  delete[] tmp;

  return num_colors;
}

/*
  Hash function definitions
*/
#define TACS_MAT_ROT_INT(x, k) (((x) << (k)) | ((x) >> (32 - (k))))
#define TACS_MAT_MIX_INT(a, b, c)                                              \
  ((a -= c, a ^= TACS_MAT_ROT_INT(c, 4), c += b, b -= a,                       \
    b ^= TACS_MAT_ROT_INT(a, 6), a += c, c -= b, c ^= TACS_MAT_ROT_INT(b, 8),  \
    b += a, a -= c, a ^= TACS_MAT_ROT_INT(c, 16), c += b, b -= a,              \
    b ^= TACS_MAT_ROT_INT(a, 19), a += c, c -= b, c ^= TACS_MAT_ROT_INT(b, 4), \
    b += a))
#define TACS_MAT_FINAL_HASH(a, b, c)                                           \
  ((c ^= b, c -= TACS_MAT_ROT_INT(b, 14), a ^= c,                              \
    a -= TACS_MAT_ROT_INT(c, 11), b ^= a, b -= TACS_MAT_ROT_INT(a, 25),        \
    c ^= b, c -= TACS_MAT_ROT_INT(b, 16), a ^= c, a -= TACS_MAT_ROT_INT(c, 4), \
    b ^= a, b -= TACS_MAT_ROT_INT(a, 14), c ^= b,                              \
    c -= TACS_MAT_ROT_INT(b, 24)))

/*
  Create a hash value for pairs of usigned integers
*/
inline uint32_t TACSMatIntegerPairHash(uint32_t u, uint32_t v) {
  uint32_t w = 0;
  TACS_MAT_MIX_INT(u, v, w);
  return TACS_MAT_FINAL_HASH(u, v, w);
}

inline uint32_t TACSIndexHashFunction(uint32_t x) {
  x ^= x >> 16;
  x *= UINT32_C(0x7feb352d);
  x ^= x >> 15;
  x *= UINT32_C(0x846ca68b);
  x ^= x >> 16;
  return x;
}

TACSIndexHash::TACSIndexHash(int approx_size, int _increment_size) {
  // Keep track of the total number of entries
  num_entries = 0;

  // Set the table size based on the non-zero estimate. This
  // is important and will
  table_size = approx_size;
  if (table_size < 1024) {
    table_size = 1024;
  }

  // Set the number of non-zeros
  if (_increment_size < 1) {
    increment_size = approx_size / 2;
  } else {
    increment_size = _increment_size;
  }
  if (increment_size < 1024) {
    increment_size = 1024;
  }

  // Allocate the table and set the entries
  table = new HashEntry *[table_size];
  memset(table, 0, table_size * sizeof(HashEntry *));

  // Allocate the memory nodes
  mem = new MemNode();
  mem->current = 0;
  mem->size = approx_size;
  mem->array = new HashEntry[approx_size];
  mem->next = NULL;
  mem_root = mem;
}

TACSIndexHash::~TACSIndexHash() {
  while (mem_root) {
    delete[] mem_root->array;
    MemNode *temp = mem_root;
    mem_root = mem_root->next;
    delete temp;
  }
  delete[] table;
}

/*
  Add an entry into the matrix.
*/
int TACSIndexHash::addEntry(int i) {
  if (i >= 0) {
    if (num_entries >= 4 * table_size) {
      reset(2 * (table_size + 1) - 1);
    }

    // Set the node attributes
    HashEntry node;
    node.i = i;
    node.next = NULL;

    // Compute the index
    uint32_t index = TACSIndexHashFunction(i);
    index = index % table_size;

    if (table[index]) {
      // Loop until finding the location where to insert the entry
      HashEntry *ptr = table[index];
      while (ptr) {
        // If the entry already exists, don't add it
        if (ptr->i == i) {
          return 0;
        }

        // If the next list element does not exist, find/allocate space
        // for the new entry
        if (!ptr->next) {
          if (mem->current >= mem->size) {
            mem->next = new MemNode();
            mem = mem->next;

            // Finish allocating the memory required
            mem->current = 0;
            mem->size = increment_size;
            mem->array = new HashEntry[increment_size];
            mem->next = NULL;
          }

          // Add the new entry and update the allocation
          mem->array[mem->current] = node;
          ptr->next = &(mem->array[mem->current]);
          mem->current++;
          num_entries++;
          return 1;
        } else {
          ptr = ptr->next;
        }
      }
    } else {
      if (mem->current >= mem->size) {
        mem->next = new MemNode();
        mem = mem->next;

        // Finish allocating the memory required
        mem->current = 0;
        mem->size = increment_size;
        mem->array = new HashEntry[increment_size];
        mem->next = NULL;
      }

      // Add the new entry and update the allocation
      mem->array[mem->current] = node;
      table[index] = &(mem->array[mem->current]);
      mem->current++;
      num_entries++;
      return 1;
    }
  }

  return 0;
}

/*
  Convert the hash table to an array
*/
void TACSIndexHash::toArray(int *_size, int **_array) {
  if (_size) {
    *_size = num_entries;
  }
  if (_array) {
    int *array = new int[num_entries];
    int i = 0;

    // Reset the locations for all the
    MemNode *node = mem_root;
    while (node) {
      for (int k = 0; k < node->current; k++, i++) {
        array[i] = node->array[k].i;
      }
      node = node->next;
    }
    *_array = array;
  }
}

void TACSIndexHash::reset(int new_table_size) {
  // Allocate the entries for the new table
  HashEntry **new_table = new HashEntry *[new_table_size];
  memset(new_table, 0, new_table_size * sizeof(HashEntry *));

  // Reset the locations for all the
  MemNode *node = mem_root;
  while (node) {
    for (int k = 0; k < node->current; k++) {
      // Reset the new entries
      HashEntry *entry = &node->array[k];
      entry->next = NULL;

      uint32_t index = TACSIndexHashFunction(entry->i);
      index = index % new_table_size;

      if (new_table[index]) {
        HashEntry *ptr = new_table[index];
        while (ptr->next) {
          ptr = ptr->next;
        }
        ptr->next = entry;
      } else {
        new_table[index] = entry;
      }
    }

    node = node->next;
  }

  // Reset the table
  delete[] table;
  table = new_table;
  table_size = new_table_size;
}

TACSMatrixHash::TACSMatrixHash(int approx_num_nonzero, int _increment_size) {
  // Keep track of the total number of entries
  num_entries = 0;

  // Set the table size based on the non-zero estimate. This
  // is important and will
  int temp = approx_num_nonzero;
  int exponent = 0;
  while ((temp >>= 1)) {
    exponent++;
  }

  table_size = (1 << exponent) - 1;
  if (table_size < 1023) {
    table_size = 1023;
  }

  // Set the number of non-zeros
  if (_increment_size < 1) {
    increment_size = approx_num_nonzero / 2;
  } else {
    increment_size = _increment_size;
  }
  if (increment_size < 1024) {
    increment_size = 1024;
  }

  // Allocate the table and set the entries
  table = new HashEntry *[table_size];
  memset(table, 0, table_size * sizeof(HashEntry *));

  // Allocate the memory nodes
  mem = new MemNode();
  mem->current = 0;
  mem->size = approx_num_nonzero;
  mem->array = new HashEntry[approx_num_nonzero];
  mem->next = NULL;
  mem_root = mem;
}

TACSMatrixHash::~TACSMatrixHash() {
  while (mem_root) {
    delete[] mem_root->array;
    MemNode *temp = mem_root;
    mem_root = mem_root->next;
    delete temp;
  }
  delete[] table;
}

/*
  Add an entry into the matrix.
*/
int TACSMatrixHash::addEntry(int i, int j) {
  if (i >= 0 && j >= 0) {
    if (num_entries >= 2 * table_size) {
      reset(2 * (table_size + 1) - 1);
    }

    // Set the node attributes
    HashEntry node;
    node.i = i;
    node.j = j;
    node.next = NULL;

    // Compute the index
    uint32_t index = TACSMatIntegerPairHash(i, j);
    index = index % table_size;

    if (table[index]) {
      // Loop until finding the location where to insert the entry
      HashEntry *ptr = table[index];
      while (ptr) {
        // If the entry already exists, don't add it
        if (ptr->i == i && ptr->j == j) {
          return 0;
        }

        // If the next list element does not exist, find/allocate space
        // for the new entry
        if (!ptr->next) {
          if (mem->current >= mem->size) {
            mem->next = new MemNode();
            mem = mem->next;

            // Finish allocating the memory required
            mem->current = 0;
            mem->size = increment_size;
            mem->array = new HashEntry[increment_size];
            mem->next = NULL;
          }

          // Add the new entry and update the allocation
          mem->array[mem->current] = node;
          ptr->next = &(mem->array[mem->current]);
          mem->current++;
          num_entries++;
          return 1;
        } else {
          ptr = ptr->next;
        }
      }
    } else {
      if (mem->current >= mem->size) {
        mem->next = new MemNode();
        mem = mem->next;

        // Finish allocating the memory required
        mem->current = 0;
        mem->size = increment_size;
        mem->array = new HashEntry[increment_size];
        mem->next = NULL;
      }

      // Add the new entry and update the allocation
      mem->array[mem->current] = node;
      table[index] = &(mem->array[mem->current]);
      mem->current++;
      num_entries++;
      return 1;
    }
  }

  return 0;
}

/*
  Convert the hash table to a non-zero CSR pattern
*/
void TACSMatrixHash::tocsr(int *_nrows, int **_rows, int **_rowp, int **_cols) {
  // Add in all the rows
  TACSIndexHash *row_hash = new TACSIndexHash(1 << 13);
  row_hash->incref();

  // Find a unique list of rows
  MemNode *node = mem_root;
  while (node) {
    for (int k = 0; k < node->current; k++) {
      row_hash->addEntry(node->array[k].i);
      row_hash->addEntry(node->array[k].j);
    }
    node = node->next;
  }

  // Get the array
  int nrows, *rows;
  row_hash->toArray(&nrows, &rows);
  row_hash->decref();

  TacsUniqueSort(nrows, rows);

  // Allocate space for the row pointer array
  int *rowp = new int[nrows + 1];
  memset(rowp, 0, (nrows + 1) * sizeof(int));

  node = mem_root;
  while (node) {
    for (int k = 0; k < node->current; k++) {
      int row = node->array[k].i;
      int *item = TacsSearchArray(row, nrows, rows);
      if (item) {
        int index = item - rows;
        node->array[k].i = index;
        rowp[index + 1]++;
      }

      int col = node->array[k].j;
      item = TacsSearchArray(col, nrows, rows);
      if (item) {
        int index = item - rows;
        node->array[k].j = index;
      }
    }
    node = node->next;
  }

  // Sum up the contributions so that rowp indexes into cols
  for (int k = 0; k < nrows; k++) {
    rowp[k + 1] += rowp[k];
  }

  // Allocate space to store the column indices
  int *cols = new int[num_entries];

  node = mem_root;
  while (node) {
    for (int k = 0; k < node->current; k++) {
      int index = node->array[k].i;
      cols[rowp[index]] = node->array[k].j;
      rowp[index]++;
    }
    node = node->next;
  }

  // Reset the columns array
  for (int k = nrows; k > 0; k--) {
    rowp[k] = rowp[k - 1];
  }
  rowp[0] = 0;

  // Restore the index values
  node = mem_root;
  while (node) {
    for (int k = 0; k < node->current; k++) {
      node->array[k].i = rows[node->array[k].i];
      node->array[k].j = rows[node->array[k].j];
    }
    node = node->next;
  }

  // Set the output
  if (_nrows) {
    *_nrows = nrows;
  }
  if (_rows) {
    *_rows = rows;
  } else {
    delete[] rows;
  }
  if (_rowp) {
    *_rowp = rowp;
  } else {
    delete[] rowp;
  }
  if (_cols) {
    *_cols = cols;
  } else {
    delete[] cols;
  }
}

void TACSMatrixHash::reset(int new_table_size) {
  // Allocate the entries for the new table
  HashEntry **new_table = new HashEntry *[new_table_size];
  memset(new_table, 0, new_table_size * sizeof(HashEntry *));

  // Reset the locations for all the
  MemNode *node = mem_root;
  while (node) {
    for (int k = 0; k < node->current; k++) {
      // Reset the new entries
      HashEntry *entry = &node->array[k];
      entry->next = NULL;

      uint32_t index = TACSMatIntegerPairHash(entry->i, entry->j);
      index = index % new_table_size;

      if (new_table[index]) {
        HashEntry *ptr = new_table[index];
        while (ptr->next) {
          ptr = ptr->next;
        }
        ptr->next = entry;
      } else {
        new_table[index] = entry;
      }
    }

    node = node->next;
  }

  // Reset the table
  delete[] table;
  table = new_table;
  table_size = new_table_size;
}
