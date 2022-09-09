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

int TacsIntegerComparator(const void *a, const void *b);

/**
  Sort an array of integers
*/
int TacsSort(int len, int *array);

/**
  Sort a list by argument

  The array indices are sorted such that values[array[i]] are
  in increasing order. Note that this is not thread-safe.

  @param len The length of the array
  @param values The array of values
  @param array The array of indices sorted so that values[array[i]] is ascending
*/
int TacsArgSort(int len, const TacsScalar *values, int *array);

/**
  Sort an array and return the number of unique design variables within
  that array start is the index of the first non-negative value and the
  return value is the number of non-negative design variables
  &dvNums[start] is the start of the new array of unique design vars...
*/
int TacsUniqueSort(int len, int *array);

/**
  Locate a value within a sorted array
*/
inline int *TacsSearchArray(int value, int len, const int *array) {
  return (int *)bsearch(&value, array, len, sizeof(int), TacsIntegerComparator);
}

/**
  Merge two sorted arrays into a single sorted array, in place.

  Two part algorithm:
  1. Find the number of duplicates between a and b
  2. Run through the list backwards placing elements into a[]
  when appropriate.

  Memory requirements: note that len(a) >= na + nb
*/
int TacsMergeSortedArrays(int na, int *a, int nb, const int *b);

/**
  Find the interval such that the given index satisfies:

  intv[k] <= index < intv[k+1]

  The intervals must be non-decreasing. Note that len is equal
  to the length of the intv array which is one more than the
  total number of intervals.
*/
int TacsFindInterval(int index, int len, const int intv[]);

/**
  Match the intervals in a list of sorted variables.

  Given the intervals, the range of the intervals and
  the list of sorted variables, local a point such that

  vars[ptr[n]] <= ownerRange[n]
*/
void TacsMatchIntervals(int size, const int range[], int nvars,
                        const int vars[], int ptr[]);

/*!
  Extend an array and copy values from the old to the new array
*/
void TacsExtendArray(int **_array, int oldlen, int newlen);
void TacsExtendArray(TacsScalar **_array, int oldlen, int newlen);

/**
  Sort and uniquify a CSR data structure.

  @param nvars The number of variables
  @param rowp Pointer into the rows of of the matrix
  @param cols The column indices
  @param remove_diag Flag indicating whether to remove diagonal entries
*/
void TacsSortAndUniquifyCSR(int nvars, int *rowp, int *cols,
                            int remove_diag = 0);

/*
  Reorder locally based on RCM-reordering
*/
int TacsComputeRCMOrder(const int nvars, const int *rowp, const int *cols,
                        int *rcm_order, int root, int n_rcm_iters);

/*
  Compute the multicolor order via a greedy algorithm
*/
int TacsComputeSerialMultiColor(const int nvars, const int *rowp,
                                const int *cols, int *colors, int *new_vars);

/*
  Hash table implementation for storing a set of integers
*/
class TACSIndexHash : public TACSObject {
 public:
  TACSIndexHash(int approx_size, int _increment_size = -1);
  ~TACSIndexHash();

  /*
    Add an entry to the hash table
  */
  int addEntry(int i);

  /*
    Convert the hash table to an array
  */
  void toArray(int *_size, int **_array);

 private:
  void reset(int new_table_size);

  // The hash entry corresponding to an (i,j) entry in the matrix
  class HashEntry {
   public:
    int i;
    HashEntry *next;
  };

  // An entry in the linked list of hash nodes
  class MemNode {
   public:
    int size;
    int current;
    HashEntry *array;
    MemNode *next;
  };

  // Number of total entries in the table
  int num_entries;

  // The number of slots avaiable in the hash table data structure.
  // This is based on the estimated number of non-zero entries in the matrix.
  int table_size;
  HashEntry **table;

  // The size of the memory increments allocated
  int increment_size;

  // Memory root and head of the linked list of memory
  MemNode *mem_root, *mem;
};

/*
  Hash table implementation for storing pairs of indices.

  The matrix pairs can be extracted to obtain a compressed sparse row matrix.
*/
class TACSMatrixHash : public TACSObject {
 public:
  TACSMatrixHash(int approx_size, int _increment_size = -1);
  ~TACSMatrixHash();

  /*
    Add an entry to the hash table
  */
  int addEntry(int i, int j);

  /*
    Convert the hash table to an array
  */
  void tocsr(int *_nrows, int **_rows, int **_rowp, int **_cols);

 private:
  void reset(int new_table_size);

  // The hash entry corresponding to an (i,j) entry in the matrix
  class HashEntry {
   public:
    int i, j;
    HashEntry *next;
  };

  // An entry in the linked list of hash nodes
  class MemNode {
   public:
    int size;
    int current;
    HashEntry *array;
    MemNode *next;
  };

  // Number of total entries in the table
  int num_entries;

  // The number of slots avaiable in the hash table data structure.
  // This is based on the estimated number of non-zero entries in the matrix.
  int table_size;
  HashEntry **table;

  // The size of the memory increments allocated
  int increment_size;

  // Memory root and head of the linked list of memory
  MemNode *mem_root, *mem;
};

#endif  // TACS_UTILITIES_H
