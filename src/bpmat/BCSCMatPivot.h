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

#ifndef TACS_BCSC_MAT_PIVOT_H
#define TACS_BCSC_MAT_PIVOT_H

/*
  Matrix classes based on compressed-sparse-column formats with
  partial pivoting
*/

#include "TACSObject.h"

/*
  A block-CSC matrix for general purpose storage. Note that this
  object cannot be used directly for matrix factorization. Instead,
  the object can be used in conjunction with the BCSCMatPivot object
  to obtain the factorization. This interface is different than the
  interface provided in TACSMat. That interface is designed to work
  efficiently with row-oriented data structures. This may cause issues
  in the future. For instance, zeroing rows will be an issue in this
  new code since that information cannot easily be derived from the
  existing data structure in the CSC matrix format.

  This matrix uses a blocked-column scheme in which the row indices
  for each blocked column are stored together. This reduces the memory
  requirements over an element-based scheme, but is less efficient
  than a full block-matrix approach (as is currently used by the BCSR
  format.)

  The input for the CSR-based data structure assumes that the input is
  a blocked-matrix where the row/columns correspond to a block (bsize
  x bsize) arrays. Note, however, that the addValues takes in the
  actual rows/ indices.
*/
class BCSCMat : public TACSObject {
 public:
  // Take in the non-zero structure of the CSC matrix
  BCSCMat(MPI_Comm _comm, int bsize, int _nrows, int _ncols, int **_colp,
          int **_rows, TacsScalar **_A = NULL);

  // Take in a CSR-type format and use it to generate the CSC data
  BCSCMat(MPI_Comm _comm, int bsize, int num_block_rows, int num_block_cols,
          const int *block_rowp, const int *block_cols);

  // Take in CSR-type format with a block-size array
  BCSCMat(MPI_Comm _comm, const int bsize[], int num_block_rows,
          int num_block_cols, const int *block_rowp, const int *block_cols);

  // Take in directly a block CSCtype format
  BCSCMat(MPI_Comm _comm, const int bsize[], int num_block_rows,
          int num_block_cols, const int *block_colp, const int *block_rows,
          int rows_sorted);

  // Create the matrix from a binary file
  BCSCMat(MPI_Comm _comm, const char *filename);

  // Copy over the values from the matrix
  BCSCMat(BCSCMat *mat);

  // Deallocate the matrix
  ~BCSCMat();

  // Add/set values in the matrix
  // ----------------------------
  void zeroEntries();
  void addMatBlockValues(int num_rows, const int brows[], int num_cols,
                         const int bcols[], const TacsScalar *mat, int ldmat,
                         int is_transpose = 0);

  // Retrieve the MPI communicator
  // -----------------------------
  MPI_Comm getMPIComm();

  // Compute multiple matrix-vector products with dense vectors
  // ----------------------------------------------------------
  void mult(TacsScalar *X, TacsScalar *Y, int vec_bsize);
  void multAdd(TacsScalar *X, TacsScalar *Z, TacsScalar *Y, int vec_bsize);

  // Retrieve the arrays from the underlying data structure
  // ------------------------------------------------------
  int getMaxBlockSize();
  void getArrays(int *_nrows, int *_ncols, int *_nblock_cols, const int **_bptr,
                 const int **_aptr, const int **_colp, const int **_rows,
                 TacsScalar **_A);

  // Archive the matrix in a file
  // ----------------------------
  void writeToFile(const char *filename);

 private:
  // Create the internal CSC data structures from the given (row and
  // column) blocked CSR data
  void initCSCfromCSR(int num_block_rows, int num_block_cols,
                      const int *block_rowp, const int *block_cols);

  MPI_Comm comm;
  int rows_sorted;  // Are the row indices sorted in each column?
  int nrows;        // The number of rows
  int ncols;        // The number of columns
  int nblock_cols;  // The number of block columns

  int max_block_size;  // Maximum block size in the matrix
  int *bptr;           // Index into the input vector (column) for each block
  int *rows;           // The non-zero rows in the data structure
  int *colp;           // Pointer to the beginning of each column

  int *aptr;      // Pointer to the initial index of each new column
  TacsScalar *A;  // Pointer to the matrix of data

  int temp_array_size;     // The temporary array size
  TacsScalar *temp_array;  // The temporary array (not always allocated)
};

/*
  This is a sparse matrix class that implements a relative simple
  partial pivoting scheme. The class uses a column-oriented storage
  format where each block column is stored in row-major order. Most
  operations are implemented using BLAS, with the most time consuming
  portions of the matrix factorization - the column updates - are
  implemented using BLAS level-3 GEMM calls.  While it is difficult to
  make a direct comparison, this code is generally 2 - 10 times slower
  than the static pivoting scheme that is used by the BCSRMat code in
  TACS. This is due to a number of factors, in part the extra overhead
  associated with the pivoting operations that generally require
  additional memory addressing and memory movement above and beyond
  the static pivoting scheme. However, this must be measured against
  the numerical stability properties associated with the partial
  pivoting approach. In general, this will likely be better for
  problems with non-symmetric matrices that are not
  diagonally-dominant.

  This data storage format is more complex than the CSR or block-CSR
  type data storage formats. This is due to the complications of
  pivoting and storing the L/U factors in a contiguous storage format
  that makes the best use of BLAS routines.

  There are two pieces of information that are required for every
  matrix element: the row and the column index. When special structure
  exists within the matrix, the rows and columns of certain portions
  of the factored matrix can be inferred enabling a reduction in the
  amount of memory required to define the non-zero pattern of the
  matrix.

  In this storage format, each row is indexed by its original row
  index to avoid having to re-number everything on the fly which would
  be cumbersome and slow. Instead, the original row index is retained
  and the position of the row, in either the L or U portions of the
  factored matrix, can be determined from the location it is stored in
  the L or U data structures. The permutation array stores the order
  of the pivots employed in the factorization.  Application of the
  factorization involves first permuting the input vector to the
  appropriate ordering, then applying the factorization.

  There are two levels of blocking that are employed in this
  factorization code that are both designed to utilize the most BLAS
  level-3 operations as possible. The first level of blocking is by
  column. Due to the structure of matrices arrising from the
  discretization of finite-element problems, the columns associated
  with a node will all have the same initial non-zero structure. If
  only partial row pivoting is employed, this identical column
  structure will be maintained throughout the factorization process.
  Therefore, it is advantageous to store the rows of the LU
  factorization in a manner that exploits this structure. Furthmore,
  this reduces the amount of memory required to index the non-zero
  space within the matrices.
*/
class BCSCMatPivot : public TACSObject {
 public:
  BCSCMatPivot(BCSCMat *C);
  ~BCSCMatPivot();

  // Factor the matrix and store the result
  double factor(double _fill = -1.0);

  // Apply the factorization to a right-hand-side vector
  void applyFactor(TacsScalar *X, int vec_bsize);

 private:
  // Apply the lower portion of the factorization
  void applyLower(TacsScalar *B, int vec_bsize, TacsScalar *temp,
                  int temp_size);

  // Apply the upper portion of the factorization
  void applyUpper(TacsScalar *B, int vec_bsize, TacsScalar *temp,
                  int temp_size);

  // Apply the column update, spa = L[node]^{-1}*spa, to the columns
  void applyNodeUpdate(int node, int node_dim, TacsScalar *spa, int spa_width,
                       TacsScalar *temp_block, int temp_block_size,
                       TacsScalar *temp_cols, int temp_cols_size);

  // Apply the column update, spa = U[node]^{-1}*spa, to the columns
  void applyNodeUpperUpdate(int node, int node_dim, TacsScalar *spa,
                            int spa_width, TacsScalar *temp_block,
                            int temp_block_size, TacsScalar *temp_cols,
                            int temp_cols_size);

  // Factor the panel matrix associated with a node
  void factorNode(int node, TacsScalar *column, int node_size, int *rows,
                  int num_rows, int diag_index);

  // Compute the non-zero pattern of a sparse right-hand-side
  // by performing a depth-first search of the graph G(L^{T})
  // to get a topological ordering of the verticies reachable from
  // the non-zero entries in rhs_rows
  void computeTopologicalOrder(int iteration, const int *rhs_rows,
                               const int num_rhs_rows, int *topo_order,
                               int *_num_node_labels, int *node_stack,
                               int *node_labels);

  // Given the topological ordering of the non-zero pattern, find the
  // positions in L[:iteration, iteration:] that update the variables
  // on the right-hand-side
  void computeColNzPattern(int iteration, const int *rhs_rows,
                           const int num_rhs_rows, const int *update_nodes,
                           const int num_nodes, int *var_labels, int *rhs_vars,
                           int *_num_vars);

  // The matrix we're going to factorize
  BCSCMat *mat;

  // The expected fill-in in the matrix
  double fill;

  // Buffer size info
  int temp_array_size;
  TacsScalar *temp_array;

  // Temporary index array information
  int *temp_iarray;

  // The block column that servese as the sparse accumulator
  TacsScalar *temp_column;

  // Data to store the LU factorization
  // ----------------------------------
  int nrows;           // The number of rows in the matrix
  int ncols;           // The number of columns in the matrix
  int nblock_cols;     // The number of block columns
  int max_block_size;  // The maximum block size
  const int *bptr;     // Pointer to the block-pointer in mat

  // Store the permutation array
  int *perm;  // row i in A -> row perm[i] in P*A

  // Store information about the node -> var, and var -> node data
  int *var_to_node;       // Super node stored for each variable (max. ncols)
  int *node_to_vars_ptr;  // Pointer to the first element of snode_to_vars
  // for each super node (ncols+1)
  int *node_to_vars;  // Super node to variable array (ncols)

  // The storage for the L/U factorization
  int max_lu_size;  // The predicted LU size based on the fill
  TacsScalar *LU;   // The storage for the L/U factors

  int max_lu_rows_size;  // Predicted final size of the lu_rows array
  int *lu_rows;          // Pointer to the rows for each node

  // Fixed-sized arrays required to point into the LU data
  int *lu_aptr;      // Pointer into the numerical columns stored in A
  int *lu_diag_ptr;  // Pointer to the diagonal of the matrix
  int *lu_col_ptr;   // Pointer to the indices of each column
};

#endif  // TACS_BCSC_MAT_PIVOT_H
