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

#ifndef TACS_BLOCK_CYCLIC_MAT_H
#define TACS_BLOCK_CYCLIC_MAT_H

#include "TACSObject.h"

/*!
  Parallel partially dense matrix format.

  This code implements a parallel partially dense matrix in block
  format. The blocks are stored in column-major fortran order. The
  matrix multiplications, and factorizations are performed in
  parallel, and in place using LAPACK/BLAS for all block
  operations. At present each block matrix is considered fully dense,
  but a sparse format could be implemented in future to store fewer
  zero off-diagonal terms.

  The parallelism is based on a 2D block cyclic format. This format is
  more difficult to implement but results in better parallelism for
  factorization. The 2D block-cyclic format repeats over blocks of
  rows and columns in the matrix. For instance, an 8-block cycle is
  arranged as follows

   __ __ __ __
  |1 |2 |3 |4 |
  |__|__|__|__|
  |5 |6 |7 |8 |
  |__|__|__|__|

  The information required in matrix operations are:
  - the processor assigned to entry (i,j)
  - the processors in the block format column that own row i
  - the processors in the block format row that own column j
  - the processor block column/row

  The factorizations are performed without pivoting for
  stability. This simplifies the algorithmic requirements, but may
  result in numerical instability. The effect of round-off errors can
  be assessed after the factorization is complete. In practice there
  should be few difficulties because of the nice properties that
  result from the use of positive definite matrices.

  The main purpose of this matrix class is to be used for an interface
  problem (the Schur-complement problem) in a precondition or direct
  solve. Therefore, there is no need to save the original values in
  the matrix. As a result, the factorization is done in place, using
  the existing non-zero pattern.

  The main requirements for the code are:

  1. Initialize the non-zero pattern from the distributed contributions
  2. Transfer the values from the distributed contributions to the
  appropriate places.
  3. Compute the factorization fill-in.
  4. Compute the factorization in parallel.
  5. Perform back-solves in parallel.
*/
class TACSBlockCyclicMat : public TACSObject {
 public:
  // Create a sparse matrix
  TACSBlockCyclicMat(MPI_Comm _comm, int csr_m, int csr_n, int csr_bsize,
                     const int *csr_vars, int nvars, const int *csr_rowp,
                     const int *csr_cols, int csr_blocks_per_block,
                     int reorder_blocks, int max_grid_size = -1);

  // Create a dense matrix
  TACSBlockCyclicMat(MPI_Comm _comm, int _nrows, int _ncols);
  ~TACSBlockCyclicMat();

  // Functions for various parts of the matrix
  // -----------------------------------------
  void getSize(int *nr, int *nc);
  void getProcessGridSize(int *_nprows, int *_npcols);
  void setMonitorFactorFlag(int flag);
  int getLocalVecSize() { return xbptr[nrows]; }

  // Get block pointers to the columns
  // ---------------------------------
  void getBlockPointers(int *_nrows, int *_ncols, const int **_bptr,
                        const int **_xbptr, const int **_perm,
                        const int **_iperm, const int **_orig_bptr);

  // Get right-hand-side index of the given global variable number
  // -------------------------------------------------------------
  int getVecIndex(int index) {
    int ioff = 0, ib = 0;
    if (orig_bptr) {
      ib = get_block_num(index, orig_bptr);
      ioff = index - orig_bptr[ib];
      ib = iperm[ib];
    } else {
      ib = get_block_num(index, bptr);
      ioff = index - bptr[ib];
    }

    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);
    if (mpi_rank == get_block_owner(ib, ib)) {
      return ioff + xbptr[ib];
    }
    return -1;
  }

  // Functions for setting values into the matrix
  // --------------------------------------------
  void zeroEntries();
  void addAllValues(int csr_bsize, int nvars, const int *vars,
                    const int *csr_rowp, const int *csr_cols, TacsScalar *vals);
  void addAlltoallValues(int csr_bsize, int nvars, const int *vars,
                         const int *csr_rowp, const int *csr_cols,
                         TacsScalar *vals);
  void setRand();

  // Matrix operations - note that factorization is in-place
  // -------------------------------------------------------
  void mult(TacsScalar *x, TacsScalar *y);
  void applyFactor(TacsScalar *x);
  void factor();

  // Given the i/j location within the matrix, determine the owner
  // -------------------------------------------------------------
  int get_block_owner(int i, int j) const {
    i = i % nprows;
    j = j % npcols;
    return proc_grid[j + i * npcols];
  }

 private:
  void init_proc_grid(int size);
  void init_nz_arrays();
  void init_row_counts();
  void merge_nz_pattern(int root, int *rowp, int *cols, int reorder_blocks);
  void compute_symbolic_factor(int **_rowp, int **_cols, int max_size);
  void init_ptr_arrays(int *rowp, int *cols);
  int get_block_num(int var, const int *ptr);
  int add_values(int rank, int i, int j, int csr_bsize, int csr_i, int csr_j,
                 TacsScalar *b);

  // Helper functions for applying the lower-triangular back-solve
  void lower_column_update(int col, TacsScalar *x, TacsScalar *xsum,
                           TacsScalar *xlocal, int *row_sum_count,
                           int *row_sum_recvd);
  void add_lower_row_sum(int row, TacsScalar *x, TacsScalar *xsum,
                         TacsScalar *xlocal, int *row_sum_count,
                         int *row_sum_recvd);

  // Helper functions for applying the upper-triangular back-solve
  void upper_column_update(int col, TacsScalar *x, TacsScalar *xsum,
                           TacsScalar *xlocal, int *row_sum_count,
                           int *row_sum_recv);
  void add_upper_row_sum(int row, TacsScalar *x, TacsScalar *xsum,
                         TacsScalar *xlocal, int *row_sum_count,
                         int *row_sum_recv);

  // Get the process row, of the provided matrix row
  int get_proc_row(int row) const { return row % nprows; }

  // Get the process column of the provided matrix column
  int get_proc_column(int col) const { return col % npcols; }

  // Get the process column and row of the given rank process
  // Return 1 upon success
  int get_proc_row_column(int rank, int *proc_row, int *proc_col) const {
    for (int i = 0; i < nprows; i++) {
      for (int j = 0; j < npcols; j++) {
        if (proc_grid[j + i * npcols] == rank) {
          *proc_row = i;
          *proc_col = j;
          return 1;
        }
      }
    }

    *proc_row = -1;
    *proc_col = -1;
    return 0;
  }

  static const int BACKSOLVE_COLUMN_SIZE = 16;
  static const int BACKSOLVE_BUFF_SIZE = 64;

  // Retrieve the block matrix at the specified entry
  TacsScalar *get_block(int rank, int i, int j);

  // The communicator for this matrix
  MPI_Comm comm;

  // This data controls how the data is assigned to the processors
  int npcols,
      nprows;  // How many processors are assigned for each grid location
  // npcols == number of columns in the process grid
  // nprows == number of rows in the process grid
  int *proc_grid;  // The processors assigned for each part of the grid

  // The global non-zero block-pattern of the matrix
  // This is global information that is duplicated on all procs
  int nrows, ncols;

  // The block sizes for the matrix
  int *bptr;      // len(bptr) = max(nrows, ncols)+1
  int max_bsize;  // max_bsize = max(bsize)

  // Pointer into the column space
  int *cbptr;

  // Pointer into the row space
  int *rbptr;

  // Pointer into the input/output vector
  int *xbptr;  // Pointer to the beginning/end of the vector

  // Additional variables if reordering has been performed
  int *orig_bptr;
  int *perm, *iperm;  // The permutation arrays

  // The non-zero off-diagonal contributions
  int *Lcolp, *Lrows;
  int *Urowp, *Ucols;

  // The locally stored components of the matrix
  TacsScalar *Dvals, *Lvals, *Uvals;
  int *dval_offset, *lval_offset, *uval_offset;
  int dval_size, uval_size, lval_size;

  // Store information about the size of the buffers required for the
  // factorization.

  // The maximum size of the recieve buffers during the factorization
  // max_ubuff_size <= nrows/npcol
  // max_lbuff_size <= nrows/nprow
  int max_ubuff_size, max_lbuff_size;

  // Monitor the time spent in the factorization process
  int monitor_factor;

  // Store information about the back-solve
  int lower_block_count, upper_block_count;
  int *lower_row_sum_count, *lower_row_sum_recv;
  int *upper_row_sum_count, *upper_row_sum_recv;
};

#endif  // TACS_BLOCK_CYCLIC_MAT_H
