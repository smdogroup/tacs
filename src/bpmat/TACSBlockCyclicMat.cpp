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

#include "TACSBlockCyclicMat.h"

#include <stdlib.h>

#include "TacsUtilities.h"
#include "tacslapack.h"

#ifdef TACS_HAS_AMD_LIBRARY
#include "amd.h"
#else
#include "AMDInterface.h"
#endif  // TACS_HAS_AMD_LIBRARY

/*
  Assemble the non-zero pattern of the matrix and pass it to the
  required processes.

  Purpose:
  --------

  This constructor accepts the distributed non-zero pattern of a
  matrix that will be assembled at a future point in time.  The
  numerical values of the entries are not yet known.  The provided
  non-zero pattern is in the form of a block-CSR format where the
  block size is provided (csr_bsize). The non-zero pattern of the
  matrix is determined by adding the contribution from all
  processors. The non-zero contribution is first collected to the root
  processor, then distributed to all remaining processors.

  input:
  comm:                  MPI communicator for the matrix
  csr_m, csr_n:          number of block-CSR rows and columns
  csr_bsize:             input block-CSR block size
  csr_vars:              global block-CSR variable numbers
  csr_nvars:             number of CSR variables
  csr_rowp:              CSR row pointer
  csr_cols:              global non-zero column indices
  csr_blocks_per_block:  number of CSR blocks per block

  Procedure:
  ----------
  1. Determine the block layout structure, populate bptr.
  2. Determine the CSR structure from the block CSR structure
  provided on each processor.
  3. Pass the non-zero patterns to the root processor.
  4. Determine the fill-ins required for the factorization process on
  the root process.
  5. Pass the block structure back to all processes.
  6. Clean up time.
*/
TACSBlockCyclicMat::TACSBlockCyclicMat(MPI_Comm _comm, int csr_m, int csr_n,
                                       int csr_bsize, const int *csr_vars,
                                       int csr_nvars, const int *csr_rowp,
                                       const int *csr_cols,
                                       int csr_blocks_per_block,
                                       int reorder_blocks, int max_grid_size) {
  comm = _comm;
  monitor_factor = 0;
  perm = iperm = orig_bptr = NULL;

  int rank = 0, size = 0;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  // Determine the process grid
  if (max_grid_size <= 0 || max_grid_size > size) {
    max_grid_size = size;
  }
  init_proc_grid(max_grid_size);

  // The total size of the matrix
  int m = csr_m * csr_bsize;
  int n = csr_n * csr_bsize;

  // Set up the bptr array
  int bsize = csr_blocks_per_block * csr_bsize;
  nrows = m / bsize;
  ncols = n / bsize;
  if (m % bsize > 0) {
    nrows++;
  }
  if (n % bsize > 0) {
    ncols++;
  }

  int len_bptr = (nrows > ncols ? nrows : ncols) + 1;
  bptr = new int[len_bptr];
  max_bsize = bsize;

  if (m >= n) {
    bptr[0] = 0;
    for (int i = 0; i < nrows; i++) {
      bptr[i + 1] = bptr[i] + bsize;
    }
    bptr[nrows] = m;
    if (nrows > 0 && bptr[nrows] - bptr[nrows - 1] > max_bsize) {
      max_bsize = bptr[nrows] - bptr[nrows - 1];
    }
  } else {
    bptr[0] = 0;
    for (int i = 0; i < ncols; i++) {
      bptr[i + 1] = bptr[i] + bsize;
    }
    bptr[ncols] = n;
    if (ncols > 0 && bptr[ncols] - bptr[ncols - 1] > max_bsize) {
      max_bsize = bptr[ncols] - bptr[ncols - 1];
    }
  }

  // Determine the block-CSR format for the block-cyclic matrix.

  // The following approach allocates more memory than is actually
  // required. However, it is fixed to at most len(csr_rowp) + nvars+1
  // new allocations. This is probably not the end of the world.
  // A tighter bound could be found.
  int *rowp = new int[nrows + 1];
  memset(rowp, 0, (nrows + 1) * sizeof(int));

  for (int ip = 0; ip < csr_nvars; ip++) {
    int i = csr_vars[ip];
    int ib = get_block_num(csr_bsize * i, bptr);

    for (int jp = csr_rowp[ip]; jp < csr_rowp[ip + 1]; jp++) {
      rowp[ib + 1]++;
    }
  }

  for (int k = 0; k < nrows; k++) {
    rowp[k + 1] += rowp[k];
  }

  int col_size = rowp[nrows];
  int *cols = new int[col_size];

  for (int ip = 0; ip < csr_nvars; ip++) {
    int i = csr_vars[ip];
    int ib = get_block_num(csr_bsize * i, bptr);

    for (int jp = csr_rowp[ip]; jp < csr_rowp[ip + 1]; jp++) {
      int j = csr_cols[jp];
      int jb = get_block_num(csr_bsize * j, bptr);
      cols[rowp[ib]] = jb;
      rowp[ib]++;
    }
  }

  // Sort and uniqify the result
  for (int k = nrows; k > 0; k--) {
    rowp[k] = rowp[k - 1];
  }
  rowp[0] = 0;

  int nodiag = 0;
  TacsSortAndUniquifyCSR(nrows, rowp, cols, nodiag);

  int root = size - 1;

  // This initializes Urowp/Ucols and Lcolp/Lrows on all procs
  merge_nz_pattern(root, rowp, cols, reorder_blocks);

  delete[] rowp;
  delete[] cols;

  // Initialize the component arrays
  init_nz_arrays();

  // Initialze data for the back-solves
  init_row_counts();
}

/*
  Create a dense matrix.
*/
TACSBlockCyclicMat::TACSBlockCyclicMat(MPI_Comm _comm, int _nrows, int _ncols) {
  comm = _comm;
  monitor_factor = 0;
  perm = iperm = orig_bptr = NULL;

  int rank = 0, size = 0;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  init_proc_grid(size);

  // Set up the data structure - for now, use a fully-dense matrix
  nrows = _nrows;
  ncols = _ncols;

  int len_bsize = (nrows > ncols ? nrows : ncols);
  bptr = new int[len_bsize + 1];
  bptr[0] = 0;
  max_bsize = 0;

  for (int i = 0; i < len_bsize; i++) {
    if (i % 3 == 0) {
      bptr[i + 1] = bptr[i] + 24;
    } else if (i % 3 == 1) {
      bptr[i + 1] = bptr[i] + 36;
    } else if (i % 3 == 2) {
      bptr[i + 1] = bptr[i] + 30;
    }
    if (bptr[i + 1] - bptr[i] > max_bsize) {
      max_bsize = bptr[i + 1] - bptr[i];
    }
  }

  // Set the offset to the solution vector itself
  xbptr = new int[len_bsize + 1];
  xbptr[0] = 0;
  for (int i = 0; i < len_bsize; i++) {
    if (rank == get_block_owner(i, i)) {
      xbptr[i + 1] = xbptr[i] + (bptr[i + 1] - bptr[i]);
    } else {
      xbptr[i + 1] = xbptr[i];
    }
  }

  // Get the process row/column
  int proc_row, proc_col;
  get_proc_row_column(rank, &proc_row, &proc_col);

  // Allocate cbptr and rbptr
  cbptr = new int[ncols + 1];
  cbptr[0] = 0;
  for (int j = 0; j < ncols; j++) {
    if (rank == get_block_owner(proc_row, j)) {
      cbptr[j + 1] = cbptr[j] + bptr[j + 1] - bptr[j];
    } else {
      cbptr[j + 1] = cbptr[j];
    }
  }
  rbptr = new int[nrows + 1];
  rbptr[0] = 0;
  for (int i = 0; i < nrows; i++) {
    if (rank == get_block_owner(i, proc_col)) {
      rbptr[i + 1] = rbptr[i] + bptr[i + 1] - bptr[i];
    } else {
      rbptr[i + 1] = rbptr[i];
    }
  }

  // Set Urowp and Ucols and Lrowp and Lcols
  Urowp = new int[nrows + 1];
  int unnz = ((nrows + 1) * nrows) / 2;
  Ucols = new int[unnz];

  Urowp[0] = 0;
  for (int i = 0, n = 0; i < nrows; i++) {
    for (int j = i + 1; j < ncols; j++, n++) {
      Ucols[n] = j;
    }
    Urowp[i + 1] = n;
  }

  // Set Lcolp and Lrows
  Lcolp = new int[ncols + 1];
  int lnnz = ((ncols + 1) * ncols) / 2;
  Lrows = new int[lnnz];

  Lcolp[0] = 0;
  for (int j = 0, n = 0; j < ncols; j++) {
    for (int i = j + 1; i < nrows; i++, n++) {
      Lrows[n] = i;
    }
    Lcolp[j + 1] = n;
  }

  // Initialize the component arrays
  init_nz_arrays();

  // Initialze data for the back-solves
  init_row_counts();
}

TACSBlockCyclicMat::~TACSBlockCyclicMat() {
  // Delete the process grid information
  delete[] proc_grid;

  // Delete the pointer to the block starting locations
  delete[] bptr;
  delete[] xbptr;
  delete[] rbptr;
  delete[] cbptr;

  // Delete the permutation information if it exists
  if (orig_bptr) {
    delete[] perm;
    delete[] iperm;
    delete[] orig_bptr;
  }

  // Delete the diagonal components
  delete[] Dvals;
  delete[] dval_offset;

  // Delete the upper triangular components
  delete[] Urowp;
  delete[] Ucols;
  delete[] uval_offset;
  delete[] Uvals;

  // Delete the lower triangular components
  delete[] Lcolp;
  delete[] Lrows;
  delete[] lval_offset;
  delete[] Lvals;

  // Delete arrays for the back-solves
  if (lower_row_sum_count) {
    delete[] lower_row_sum_count;
    delete[] lower_row_sum_recv;
    delete[] upper_row_sum_count;
    delete[] upper_row_sum_recv;
  }
}

/*
  Merge the non-zero pattern of the matrices to the root process.

  At this point rowp/cols point to the local contribution to the
  non-zero pattern. Send these contributions from all processors to
  the root process, at which point the full set of fill ins can be
  computed. Finally, the arrays Urowp/Ucols and Lcolp/Lrows can be
  allocated and set.
*/
void TACSBlockCyclicMat::merge_nz_pattern(int root, int *rowp, int *cols,
                                          int reorder_blocks) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // Count up the maximum size of the CSR array
  int col_size = rowp[nrows];
  int root_col_size = 0;
  MPI_Reduce(&col_size, &root_col_size, 1, MPI_INT, MPI_SUM, root, comm);

  int *root_rowp = NULL;
  int *root_cols = NULL;
  int *root_all_rowp = NULL;

  int *recv_row = NULL;
  int *recv_ptr = NULL;
  int *recv_count = NULL;

  if (rank == root) {
    root_rowp = new int[nrows + 1];
    root_rowp[0] = 0;

    root_cols = new int[root_col_size];
    root_all_rowp = new int[(nrows + 1) * size];

    recv_row = new int[size * ncols];
    recv_ptr = new int[size + 1];
    recv_count = new int[size];
  }

  MPI_Gather(rowp, nrows + 1, MPI_INT, root_all_rowp, (nrows + 1), MPI_INT,
             root, comm);

  for (int i = 0; i < nrows; i++) {
    if (rank == root) {
      recv_ptr[0] = 0;
      for (int k = 0; k < size; k++) {
        recv_count[k] = (root_all_rowp[i + 1 + k * (nrows + 1)] -
                         root_all_rowp[i + k * (nrows + 1)]);
        recv_ptr[k + 1] = recv_ptr[k] + recv_count[k];
      }
    }

    int row_size = rowp[i + 1] - rowp[i];
    MPI_Gatherv(&cols[rowp[i]], row_size, MPI_INT, recv_row, recv_count,
                recv_ptr, MPI_INT, root, comm);

    if (rank == root) {
      // Sort and uniquify this row
      int row_len = recv_ptr[size];
      row_len = TacsUniqueSort(row_len, recv_row);

      // Now add this row to the CSR data structure
      root_rowp[i + 1] = root_rowp[i] + row_len;
      for (int j = 0, k = root_rowp[i]; j < row_len; j++, k++) {
        root_cols[k] = recv_row[j];
      }
    }
  }

  if (reorder_blocks && nrows == ncols) {
    // Allocate the permutation array
    perm = new int[nrows];
    iperm = new int[ncols];
    orig_bptr = new int[nrows + 1];
  } else {
    reorder_blocks = 0;
  }

  if (rank == root) {
    delete[] root_all_rowp;
    delete[] recv_row;
    delete[] recv_ptr;
    delete[] recv_count;

    int m = bptr[nrows], n = bptr[ncols];
    if (reorder_blocks && n == m) {
      if (m && n) {
#ifdef TACS_HAS_AMD_LIBRARY
        // Use AMD to compute the reordering of the variables.
        double control[AMD_CONTROL], info[AMD_INFO];
        amd_defaults(control);  // Use the default values
        amd_order(nrows, root_rowp, root_cols, perm, control, info);
#else
        // Allocate temporary space
        int *tmp_rowp = new int[nrows + 1];
        int *tmp_cols = new int[root_rowp[nrows]];
        memcpy(tmp_rowp, root_rowp, (nrows + 1) * sizeof(int));
        memcpy(tmp_cols, root_cols, root_rowp[nrows] * sizeof(int));

        // This call destroys data in tmp_rowp/tmp_cols
        int use_exact_degree = 0;
        amd_order_interface(nrows, tmp_rowp, tmp_cols, perm, NULL, 0, 0, NULL,
                            NULL, NULL, use_exact_degree);

        // Free the temporary space
        delete[] tmp_rowp;
        delete[] tmp_cols;
#endif  // TACS_HAS_AMD_LIBRARY
        // Reorder the non-zero pattern to correspond to the new ordering
        // perm:  new variable i -> old variable perm[i]
        // iperm: old variable i -> new variable iperm[i]
        for (int i = 0; i < nrows; i++) {
          iperm[perm[i]] = i;
        }

        int *temp_rowp = new int[nrows + 1];
        int *temp_cols = new int[root_col_size];

        temp_rowp[0] = 0;
        for (int i = 0, p = 0; i < nrows; i++) {
          int row = perm[i];
          int size = root_rowp[row + 1] - root_rowp[row];
          for (int jp = root_rowp[row]; jp < root_rowp[row + 1]; jp++, p++) {
            temp_cols[p] = iperm[root_cols[jp]];
          }
          // Sort the row
          temp_rowp[i + 1] = p;
          size = TacsUniqueSort(size, &temp_cols[temp_rowp[i]]);
          if (size != temp_rowp[i + 1] - temp_rowp[i]) {
            printf(
                "[%d] TACSBlockCyclicMat: problem with the "
                "permutation array\n",
                rank);
          }
        }

        delete[] root_rowp;
        delete[] root_cols;

        root_rowp = temp_rowp;
        root_cols = temp_cols;
      }
    }

    // perform the symbolic factorization
    int init_nnz = root_rowp[nrows];
    compute_symbolic_factor(&root_rowp, &root_cols, root_col_size);
    int final_nnz = root_rowp[nrows];

    if (m && n) {
      printf(
          "[%d] TACSBlockCyclicMat: (%d,%d) Initial density: %4.3f "
          "factor fill in: %4.3f\n",
          root, m, n, (1.0 * init_nnz) / (nrows * ncols),
          (1.0 * (final_nnz - init_nnz)) / init_nnz);
    }

    // create the Urowp/Ucols and Lcolp/Lrows arrays
    init_ptr_arrays(root_rowp, root_cols);

    delete[] root_rowp;
    delete[] root_cols;
  }

  // Reorder the blocks
  if (reorder_blocks) {
    MPI_Bcast(perm, nrows, MPI_INT, root, comm);

    if (rank != root) {
      for (int i = 0; i < nrows; i++) {
        iperm[perm[i]] = i;
      }
    }

    // Set the bptr array with the new permutation
    for (int i = 0; i < nrows; i++) {
      int ip = perm[i];
      orig_bptr[i + 1] = bptr[ip + 1] - bptr[ip];
    }
    orig_bptr[0] = 0;
    for (int i = 0; i < nrows; i++) {
      orig_bptr[i + 1] += orig_bptr[i];
    }

    // Flip the two
    int *temp = bptr;
    bptr = orig_bptr;
    orig_bptr = temp;
  }

  // Set the offset to the vector
  xbptr = new int[nrows + 1];
  xbptr[0] = 0;
  for (int i = 0; i < nrows; i++) {
    if (rank == get_block_owner(i, i)) {
      xbptr[i + 1] = xbptr[i] + (bptr[i + 1] - bptr[i]);
    } else {
      xbptr[i + 1] = xbptr[i];
    }
  }

  // Get the process row/column
  int proc_row, proc_col;
  get_proc_row_column(rank, &proc_row, &proc_col);

  // Allocate cbptr and rbptr
  cbptr = new int[ncols + 1];
  cbptr[0] = 0;
  for (int j = 0; j < ncols; j++) {
    if (rank == get_block_owner(proc_row, j)) {
      cbptr[j + 1] = cbptr[j] + bptr[j + 1] - bptr[j];
    } else {
      cbptr[j + 1] = cbptr[j];
    }
  }
  rbptr = new int[nrows + 1];
  rbptr[0] = 0;
  for (int i = 0; i < nrows; i++) {
    if (rank == get_block_owner(i, proc_col)) {
      rbptr[i + 1] = rbptr[i] + bptr[i + 1] - bptr[i];
    } else {
      rbptr[i + 1] = rbptr[i];
    }
  }

  // Broadcast the nonzero pattern
  if (rank != root) {
    Urowp = new int[nrows + 1];
    Lcolp = new int[ncols + 1];
  }

  MPI_Bcast(Urowp, nrows + 1, MPI_INT, root, comm);
  MPI_Bcast(Lcolp, ncols + 1, MPI_INT, root, comm);

  // Allocate space for the remaining arrays
  if (rank != root) {
    Ucols = new int[Urowp[nrows]];
    Lrows = new int[Lcolp[ncols]];
  }

  MPI_Bcast(Ucols, Urowp[nrows], MPI_INT, root, comm);
  MPI_Bcast(Lrows, Lcolp[ncols], MPI_INT, root, comm);
}

/*
  Perform the symbolic factorization of the non-zero pattern in
  rowp/cols. This computation is performed in place.

  This may take extra time because of the number of copies that are
  required. The array must be shifted to account for the new entries
  generated at each step in the factorization process.  This shift is
  delayed until all the new entries for the new row are processed.
*/
void TACSBlockCyclicMat::compute_symbolic_factor(int **_rowp, int **_cols,
                                                 int max_size) {
  int *rowp = *_rowp;
  int *cols = *_cols;

  int *diag = new int[nrows];
  int *rcols = new int[ncols];

  for (int i = 0; i < nrows; i++) {
    int nr = 0;  // Number of entries in the current row

    int diag_flag = 0;
    for (int j = rowp[i]; j < rowp[i + 1]; j++) {
      if (cols[j] == i) {
        diag_flag = 1;
      }
      rcols[nr] = cols[j];
      nr++;
    }

    if (!diag_flag) {
      nr = TacsMergeSortedArrays(nr, rcols, 1, &i);
    }

    // Perform the symbolic factorization - generating new entries
    for (int j = 0; rcols[j] < i; j++) {
      int p = j + 1;                   // The index into rcols
      int k_end = rowp[rcols[j] + 1];  // the end of row number cols[j]

      // Start with the first entry after the diagonal in row, cols[j]
      // k is the index into cols for row cols[j]
      for (int k = diag[rcols[j]] + 1; k < k_end; k++) {
        // Increment p to an entry where we may have cols[k] == rcols[p]
        while ((p < nr) && (rcols[p] < cols[k])) {
          p++;
        }

        if (p == nr) {
          // Insert the new element at the end of the row
          rcols[p] = cols[k];
          nr++;
        } else if (rcols[p] != cols[k]) {
          // Insert the new entry into the list, but keep the list sorted
          for (int n = nr; n > p; n--) {
            rcols[n] = rcols[n - 1];
          }

          rcols[p] = cols[k];
          nr++;
        }
      }
    }

    // Look for the diagonal entry
    for (int j = 0; j < nr; j++) {
      if (rcols[j] == i) {
        diag[i] = rowp[i] + j;
        break;
      }
    }

    // Check if the size will be exceeded by adding the new elements
    int current_size = rowp[nrows];
    int new_entries = nr - (rowp[i + 1] - rowp[i]);
    if (current_size + new_entries > max_size) {
      int mat_ext = 2 * current_size;
      if (nr > mat_ext) {
        mat_ext = nr;
      }
      max_size = max_size + mat_ext;
      TacsExtendArray(&cols, current_size, max_size);
    }

    // Move everthing forward so that there's enough space
    int k = current_size + new_entries - 1;  // Place new entries here
    for (int j = nrows - 1; j >= i; j--) {
      int start = rowp[j];
      int end = rowp[j + 1] - 1;
      rowp[j + 1] = k + 1;  // The position of the new end of the array

      for (; end >= start; end--, k--) {
        cols[k] = cols[end];
      }
    }

    // Now, put the new entries into the cols/levs arrays
    for (int j = rowp[i], k = 0; k < nr; j++, k++) {
      cols[j] = rcols[k];
    }
  }

  *_rowp = rowp;
  *_cols = cols;

  delete[] diag;
  delete[] rcols;
}

/*
  Given the non-zero pattern in rowp/cols create Urowp/Ucols and
  Lcolp/Lrows arrays. Note that rowp/cols must be sorted row-wise.

  Note that the row indices in Urowp/Ucols and the column indices
  in Lcolp/Lrows will be sorted if rowp/cols are sorted - which
  they must be!
*/
void TACSBlockCyclicMat::init_ptr_arrays(int *rowp, int *cols) {
  Urowp = new int[nrows + 1];
  Lcolp = new int[ncols + 1];
  memset(Urowp, 0, (nrows + 1) * sizeof(int));
  memset(Lcolp, 0, (ncols + 1) * sizeof(int));

  for (int i = 0; i < nrows; i++) {
    for (int jp = rowp[i]; jp < rowp[i + 1]; jp++) {
      int j = cols[jp];
      if (j > i) {  // Upper part
        Urowp[i + 1]++;
      } else if (j < i) {  // Lower part
        Lcolp[j + 1]++;
      }
    }
  }

  // Count up the contributions to the arrays
  for (int i = 0; i < nrows; i++) {
    Urowp[i + 1] += Urowp[i];
  }

  for (int i = 0; i < ncols; i++) {
    Lcolp[i + 1] += Lcolp[i];
  }

  // Allocate the arrays
  Ucols = new int[Urowp[nrows]];
  Lrows = new int[Lcolp[ncols]];

  for (int i = 0; i < nrows; i++) {
    for (int jp = rowp[i]; jp < rowp[i + 1]; jp++) {
      int j = cols[jp];
      if (j > i) {  // Upper part
        Ucols[Urowp[i]] = j;
        Urowp[i]++;
      } else if (j < i) {  // Lower part
        Lrows[Lcolp[j]] = i;
        Lcolp[j]++;
      }
    }
  }

  // Fix the pointer arrays
  for (int i = nrows; i > 0; i--) {
    Urowp[i] = Urowp[i - 1];
  }
  Urowp[0] = 0;

  for (int i = ncols; i > 0; i--) {
    Lcolp[i] = Lcolp[i - 1];
  }
  Lcolp[0] = 0;
}

/*
  Set up the process grid. This has to assign a unique processor to
  each process block. Therefore, in general there may be processes
  left without blocks. For np = 1, 2, 3, 4, 6, 8, 9, 12, 15, 16,...
  it should work fine. More flexibility with the block assignment
  should be considered for future versions of the code.
*/
void TACSBlockCyclicMat::init_proc_grid(int size) {
  // Special cases
  if (size == 1) {
    nprows = 1;
    npcols = 1;
  } else {
    // Favor more columns than rows
    int n = 1;
    while (n * (n + 1) > 0) {
      if (n * (n + 1) <= size && size < (n + 1) * (n + 2)) {
        break;
      }
      n++;
    }

    nprows = n;
    npcols = size / n;
  }

  // Allocate the process grid and assign the processors
  proc_grid = new int[nprows * npcols];

  for (int i = 0, p = 0; i < nprows; i++) {
    for (int j = 0; j < npcols; j++, p++) {
      proc_grid[j + i * npcols] = p;
    }
  }
}

/*
  Get the size of the matrix
*/
void TACSBlockCyclicMat::getSize(int *nr, int *nc) {
  if (nr) {
    *nr = bptr[nrows];
  }
  if (nc) {
    *nc = bptr[ncols];
  }
}

/*
  Retrieve the size of the process grid.
*/
void TACSBlockCyclicMat::getProcessGridSize(int *_nprows, int *_npcols) {
  if (_nprows) {
    *_nprows = nprows;
  }
  if (_npcols) {
    *_npcols = npcols;
  }
}

/*
  Retrieve the block pointers
*/
void TACSBlockCyclicMat::getBlockPointers(int *_nrows, int *_ncols,
                                          const int **_bptr, const int **_xbptr,
                                          const int **_perm, const int **_iperm,
                                          const int **_orig_bptr) {
  if (_nrows) {
    *_nrows = nrows;
  }
  if (_ncols) {
    *_ncols = ncols;
  }
  if (_bptr) {
    *_bptr = bptr;
  }
  if (_xbptr) {
    *_xbptr = xbptr;
  }
  if (_perm) {
    *_perm = perm;
  }
  if (_iperm) {
    *_iperm = iperm;
  }
  if (_orig_bptr) {
    *_orig_bptr = orig_bptr;
  }
}

/*
  Set the flag that prints out the factorization time
*/
void TACSBlockCyclicMat::setMonitorFactorFlag(int flag) {
  monitor_factor = flag;
}

/*
  This function performs several initialization tasks, including
  determining the number of matrix elements that are stored locally,
  and computing pointers into the arrays that store these values.
  This function can only be used after the process grid is set up, the
  bptr array is initialized AND the non-zero pattern is set in
  Lcolp/Lrows and Urowp/Ucols.

  This function initializes:
  - dval_offset/uval_offset/lval_offset: Integer pointers into the
  arrays Dvals, Uvals and Lvals
  - Dvals, Uvals and Lvals: The arrays storing the matrix components
  - max_ubuff_size, max_lbuff_size: Maximum size required for
  buffers during factorization
*/
void TACSBlockCyclicMat::init_nz_arrays() {
  int rank;
  MPI_Comm_rank(comm, &rank);

  dval_offset = new int[nrows + 1];
  uval_offset = new int[Urowp[nrows] + 1];
  lval_offset = new int[Lcolp[ncols] + 1];

  // Calculate the diagonal offsets
  dval_offset[0] = 0;
  for (int i = 0; i < nrows; i++) {
    int bi = bptr[i + 1] - bptr[i];
    if (rank == get_block_owner(i, i)) {
      dval_offset[i + 1] = dval_offset[i] + bi * bi;
    } else {
      dval_offset[i + 1] = dval_offset[i];
    }
  }

  // Allocate the diagonal values
  dval_size = dval_offset[nrows];
  Dvals = new TacsScalar[dval_size];
  memset(Dvals, 0, dval_size * sizeof(TacsScalar));

  // Calculate the off-diagonal offsets
  max_ubuff_size = 0;
  uval_offset[0] = 0;
  for (int i = 0, n = 0; i < nrows; i++) {
    int bi = bptr[i + 1] - bptr[i];
    int ubuff_size = 0;
    for (int jp = Urowp[i]; jp < Urowp[i + 1]; jp++, n++) {
      int j = Ucols[jp];
      int bj = bptr[j + 1] - bptr[j];
      if (rank == get_block_owner(i, j)) {
        uval_offset[n + 1] = uval_offset[n] + bi * bj;
        ubuff_size += bi * bj;
      } else {
        uval_offset[n + 1] = uval_offset[n];
      }
    }
    if (ubuff_size > max_ubuff_size) {
      max_ubuff_size = ubuff_size;
    }
  }

  // Allocate the upper triangular components
  uval_size = uval_offset[Urowp[nrows]];
  Uvals = new TacsScalar[uval_size];
  memset(Uvals, 0, uval_size * sizeof(TacsScalar));

  // Calculate the lower triangular size
  lval_offset[0] = 0;
  max_lbuff_size = 0;
  for (int j = 0, n = 0; j < ncols; j++) {
    int bj = bptr[j + 1] - bptr[j];
    int lbuff_size = 0;
    for (int ip = Lcolp[j]; ip < Lcolp[j + 1]; ip++, n++) {
      int i = Lrows[ip];
      int bi = bptr[i + 1] - bptr[i];
      if (rank == get_block_owner(i, j)) {
        lval_offset[n + 1] = lval_offset[n] + bi * bj;
        lbuff_size += bi * bj;
      } else {
        lval_offset[n + 1] = lval_offset[n];
      }
    }
    if (lbuff_size > max_lbuff_size) {
      max_lbuff_size = lbuff_size;
    }
  }

  // Allocate space for the lower triangular components
  lval_size = lval_offset[Lcolp[ncols]];
  Lvals = new TacsScalar[lval_size];
  memset(Lvals, 0, lval_size * sizeof(TacsScalar));

  int max_buff[2], max_buff_all[2];
  max_buff[0] = max_lbuff_size;
  max_buff[1] = max_ubuff_size;

  MPI_Allreduce(max_buff, max_buff_all, 2, MPI_INT, MPI_MAX, comm);

  max_lbuff_size = max_buff_all[0];
  max_ubuff_size = max_buff_all[1];
}

/*
  Zero all the matrix entries.
*/
void TACSBlockCyclicMat::zeroEntries() {
  memset(Dvals, 0, dval_size * sizeof(TacsScalar));
  memset(Lvals, 0, lval_size * sizeof(TacsScalar));
  memset(Uvals, 0, uval_size * sizeof(TacsScalar));
}

/*
  Add values into the matrix. This function is collective on all
  processes in the block-cyclic comm.

  First, add the on-process components. Next, add the off-process
  components from each process. This approach tries to minimize the
  amount of extra memory required for the off-process contributions.
  This code uses MPI_Gatherv for each process rather than a single
  call to MPI_Alltoallv which may be faster, but requires more memory.
*/
void TACSBlockCyclicMat::addAllValues(int csr_bsize, int nvars,
                                      const int *csr_vars, const int *csr_rowp,
                                      const int *csr_cols, TacsScalar *vals) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int b2 = csr_bsize * csr_bsize;

  // Count up the number of recvs for each processor
  int *recv_count = new int[size];
  int *recv_ptr = new int[size + 1];

  int csr_size = csr_rowp[nvars];
  int *jblock = new int[csr_size];

  int *send_counts = new int[size];
  memset(send_counts, 0, size * sizeof(int));

  // Go through and add any contributions from the local process
  for (int ip = 0; ip < nvars; ip++) {
    int i = csr_bsize * csr_vars[ip];
    int ioff = 0, ib = 0;
    if (orig_bptr) {
      ib = get_block_num(i, orig_bptr);
      ioff = i - orig_bptr[ib];
      ib = iperm[ib];
    } else {
      ib = get_block_num(i, bptr);
      ioff = i - bptr[ib];
    }

    // Determine the local block for vars
    for (int jp = csr_rowp[ip]; jp < csr_rowp[ip + 1]; jp++) {
      int j = csr_bsize * csr_cols[jp];
      int joff = 0;

      // Get the block numbers
      if (orig_bptr) {
        int jb = get_block_num(j, orig_bptr);
        joff = j - orig_bptr[jb];
        jblock[jp] = iperm[jb];
      } else {
        jblock[jp] = get_block_num(j, bptr);
        joff = j - bptr[jblock[jp]];
      }

      int owner = get_block_owner(ib, jblock[jp]);
      if (owner == rank) {
        add_values(rank, ib, jblock[jp], csr_bsize, ioff, joff, &vals[b2 * jp]);
      } else {
        send_counts[owner]++;
      }
    }
  }

  int max_send_size = send_counts[0];
  for (int k = 1; k < size; k++) {
    if (send_counts[k] > max_send_size) {
      max_send_size = send_counts[k];
    }
  }

  delete[] send_counts;

  int *send_ivals = new int[max_send_size];
  int *send_jvals = new int[max_send_size];
  TacsScalar *send_vals = new TacsScalar[b2 * max_send_size];

  // Iterate over each process
  for (int k = 0; k < size; k++) {
    int send_count = 0;

    if (rank != k) {
      // Figure out what to add where
      for (int ip = 0; ip < nvars; ip++) {
        int i = csr_bsize * csr_vars[ip];
        int ib;
        if (orig_bptr) {
          ib = get_block_num(i, orig_bptr);
          ib = iperm[ib];
        } else {
          ib = get_block_num(i, bptr);
        }

        // Determine the local block for vars
        for (int jp = csr_rowp[ip]; jp < csr_rowp[ip + 1]; jp++) {
          int j = csr_bsize * csr_cols[jp];

          if (get_block_owner(ib, jblock[jp]) == k) {
            send_ivals[send_count] = i;
            send_jvals[send_count] = j;
            memcpy(&send_vals[b2 * send_count], &vals[b2 * jp],
                   b2 * sizeof(TacsScalar));
            send_count++;
          }
        }
      }
    }

    MPI_Gather(&send_count, 1, MPI_INT, recv_count, 1, MPI_INT, k, comm);

    // Count up the reciving information
    int nrecv = 0;
    int *recv_ivals = NULL;
    int *recv_jvals = NULL;
    TacsScalar *recv_vals = NULL;

    if (rank == k) {
      recv_ptr[0] = 0;
      for (int n = 0; n < size; n++) {
        recv_ptr[n + 1] = recv_ptr[n] + recv_count[n];
      }

      nrecv = recv_ptr[size];

      recv_ivals = new int[nrecv];
      recv_jvals = new int[nrecv];
      recv_vals = new TacsScalar[b2 * nrecv];
    }

    MPI_Gatherv(send_ivals, send_count, MPI_INT, recv_ivals, recv_count,
                recv_ptr, MPI_INT, k, comm);
    MPI_Gatherv(send_jvals, send_count, MPI_INT, recv_jvals, recv_count,
                recv_ptr, MPI_INT, k, comm);

    // Reset the array sizes to send the blocks
    send_count *= b2;
    for (int n = 0; n < size; n++) {
      recv_count[n] *= b2;
      recv_ptr[n] *= b2;
    }

    // Send the values
    MPI_Gatherv(send_vals, send_count, TACS_MPI_TYPE, recv_vals, recv_count,
                recv_ptr, TACS_MPI_TYPE, k, comm);

    if (rank == k) {
      // Allocate the sending/receiving arrays
      for (int n = 0; n < nrecv; n++) {
        int i = recv_ivals[n];
        int j = recv_jvals[n];

        int ib, jb;  // The blocks
        int ioff, joff;
        if (orig_bptr) {
          ib = get_block_num(i, orig_bptr);
          jb = get_block_num(j, orig_bptr);
          ioff = i - orig_bptr[ib];
          joff = j - orig_bptr[jb];
          ib = iperm[ib];
          jb = iperm[jb];
        } else {
          ib = get_block_num(i, bptr);
          jb = get_block_num(j, bptr);
          ioff = i - bptr[ib];
          joff = j - bptr[jb];
        }

        // Add this to the correct block ...
        add_values(rank, ib, jb, csr_bsize, ioff, joff, &recv_vals[b2 * n]);
      }

      delete[] recv_ivals;
      delete[] recv_jvals;
      delete[] recv_vals;
    }
  }

  delete[] jblock;

  delete[] recv_count;
  delete[] recv_ptr;

  delete[] send_ivals;
  delete[] send_jvals;
  delete[] send_vals;
}

/*
  Add values into the matrix. This function is collective on all matrix
  processes.

  This function requires more memory than the function above, but may
  be faster - especially when more processors are used.

  This function first adds the contributions to the on-processor
  blocks, while computing the number of off-processor blocks for each
  proc. Next, buffers are allocated that can store the entire
  off-process contribution to all non-local processes and a receive
  buffer is allocated to store all incoming contributions. Next,
  MPI_Alltoallv calls are made to pass the information to all required
  processes. The receive buffers are then added to the local components.
  All allocated memory is freed.
*/
void TACSBlockCyclicMat::addAlltoallValues(int csr_bsize, int nvars,
                                           const int *csr_vars,
                                           const int *csr_rowp,
                                           const int *csr_cols,
                                           TacsScalar *vals) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int b2 = csr_bsize * csr_bsize;

  // Count up the number of recvs for each processor
  int *recv_counts = new int[size];
  int *recv_ptr = new int[size + 1];

  int *send_counts = new int[size];
  int *send_ptr = new int[size + 1];
  memset(send_counts, 0, size * sizeof(int));

  // Go through and add any contributions from the local process
  for (int ip = 0; ip < nvars; ip++) {
    int i = csr_bsize * csr_vars[ip];
    int ioff = 0, ib = 0;
    if (orig_bptr) {
      ib = get_block_num(i, orig_bptr);
      ioff = i - orig_bptr[ib];
      ib = iperm[ib];
    } else {
      ib = get_block_num(i, bptr);
      ioff = i - bptr[ib];
    }

    // Determine the local block for vars
    for (int jp = csr_rowp[ip]; jp < csr_rowp[ip + 1]; jp++) {
      int j = csr_bsize * csr_cols[jp];
      int joff = 0, jb = 0;

      // Get the block numbers
      if (orig_bptr) {
        jb = get_block_num(j, orig_bptr);
        joff = j - orig_bptr[jb];
        jb = iperm[jb];
      } else {
        jb = get_block_num(j, bptr);
        joff = j - bptr[jb];
      }

      int owner = get_block_owner(ib, jb);
      if (owner == rank) {
        add_values(rank, ib, jb, csr_bsize, ioff, joff, &vals[b2 * jp]);
      } else {
        send_counts[owner]++;
      }
    }
  }

  // Send the counts to all processors
  MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);

  // Sum up the send and recv pointers
  send_ptr[0] = 0;
  recv_ptr[0] = 0;
  for (int k = 0; k < size; k++) {
    send_ptr[k + 1] = send_ptr[k] + send_counts[k];
    recv_ptr[k + 1] = recv_ptr[k] + recv_counts[k];
  }

  // Create a buffer that is large enough to send everything at once
  // Zero out the send counter
  memset(send_counts, 0, size * sizeof(int));
  int *send_index = new int[2 * send_ptr[size]];
  TacsScalar *send_vals = new TacsScalar[b2 * send_ptr[size]];

  // Go back through and copy over values to pass to all the processes
  for (int ip = 0; ip < nvars; ip++) {
    int i = csr_bsize * csr_vars[ip];
    int ib = 0;
    if (orig_bptr) {
      ib = get_block_num(i, orig_bptr);
      ib = iperm[ib];
    } else {
      ib = get_block_num(i, bptr);
    }

    // Determine the local block for vars
    for (int jp = csr_rowp[ip]; jp < csr_rowp[ip + 1]; jp++) {
      int j = csr_bsize * csr_cols[jp];
      int jb = 0;

      // Get the block numbers
      if (orig_bptr) {
        jb = get_block_num(j, orig_bptr);
        jb = iperm[jb];
      } else {
        jb = get_block_num(j, bptr);
      }

      int owner = get_block_owner(ib, jb);
      if (owner != rank) {
        int sc = send_ptr[owner] + send_counts[owner];
        send_counts[owner]++;

        // Set the index values
        send_index[2 * sc] = i;
        send_index[2 * sc + 1] = j;

        // Copy the values to the send buffer
        memcpy(&send_vals[b2 * sc], &vals[b2 * jp], b2 * sizeof(TacsScalar));
      }
    }
  }

  // Send the indices and the values
  int *send_counts2 = new int[size];
  int *send_ptr2 = new int[size];
  int *recv_counts2 = new int[size];
  int *recv_ptr2 = new int[size];

  // Adjust the pointers/counts to receive the (i,j) indices
  for (int k = 0; k < size; k++) {
    send_counts2[k] = 2 * send_counts[k];
    send_ptr2[k] = 2 * send_ptr[k];
    recv_counts2[k] = 2 * recv_counts[k];
    recv_ptr2[k] = 2 * recv_ptr[k];
  }

  int *recv_index = new int[2 * recv_ptr[size]];
  TacsScalar *recv_vals = new TacsScalar[b2 * recv_ptr[size]];

  MPI_Alltoallv(send_index, send_counts2, send_ptr2, MPI_INT, recv_index,
                recv_counts2, recv_ptr2, MPI_INT, comm);

  // Adjust the pointers to receive the values
  for (int k = 0; k < size; k++) {
    send_counts2[k] = b2 * send_counts[k];
    send_ptr2[k] = b2 * send_ptr[k];
    recv_counts2[k] = b2 * recv_counts[k];
    recv_ptr2[k] = b2 * recv_ptr[k];
  }

  MPI_Alltoallv(send_vals, send_counts2, send_ptr2, TACS_MPI_TYPE, recv_vals,
                recv_counts2, recv_ptr2, TACS_MPI_TYPE, comm);

  // Free all the pointers that are not required anymore
  delete[] send_counts2;
  delete[] send_ptr2;
  delete[] recv_counts2;
  delete[] recv_ptr2;
  delete[] send_counts;
  delete[] send_ptr;

  // Free the send buffer indices and values
  delete[] send_index;
  delete[] send_vals;

  int nrecv = recv_ptr[size];

  // Add the values from the receiving arrays
  for (int n = 0; n < nrecv; n++) {
    int i = recv_index[2 * n];
    int j = recv_index[2 * n + 1];

    int ib, jb;  // The block indices
    int ioff, joff;
    if (orig_bptr) {
      ib = get_block_num(i, orig_bptr);
      jb = get_block_num(j, orig_bptr);
      ioff = i - orig_bptr[ib];
      joff = j - orig_bptr[jb];
      ib = iperm[ib];
      jb = iperm[jb];
    } else {
      ib = get_block_num(i, bptr);
      jb = get_block_num(j, bptr);
      ioff = i - bptr[ib];
      joff = j - bptr[jb];
    }

    // Add this to the correct block ...
    add_values(rank, ib, jb, csr_bsize, ioff, joff, &recv_vals[b2 * n]);
  }

  delete[] recv_index;
  delete[] recv_vals;
  delete[] recv_counts;
  delete[] recv_ptr;
}

/*
  Determine the block number i such that var is within the interval:

  ptr[i] <= var < ptr[i+1]
*/
int TACSBlockCyclicMat::get_block_num(int var, const int *ptr) {
  int high = (ncols > nrows ? ncols : nrows);
  int low = 0;

  if ((var < ptr[low]) || (var >= ptr[high])) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    fprintf(stderr, "[%d] TACSBlockCyclicMat::get_block_num(%d) out of range\n",
            rank, var);
    return -1;
  }

  if (ptr[low] == var) {
    return low;
  }

  while (high - low > 1) {
    int mid = low + (int)((high - low) / 2);

    if (ptr[mid] == var) {
      return mid;  // quick return
    } else if (var > ptr[mid]) {
      low = mid;
    } else {
      high = mid;
    }
  }

  return low;
}

/*
  Given the block a[] from the input CSR format, add a[] to the given
  block (i, j) within the cyclic block format.  This function checks
  to see whether the process owns the current block and whether it
  exists within the non-zero format.

  The input block matrix a[] is in row-major order (C-order), and the
  storage format is in column-major (fortran-order).
*/
int TACSBlockCyclicMat::add_values(int rank, int i, int j, int csr_bsize,
                                   int ioff, int joff, TacsScalar *a) {
  TacsScalar *A = get_block(rank, i, j);

  if (A) {
    int bi = bptr[i + 1] - bptr[i];
    int bj = bptr[j + 1] - bptr[j];

    if ((ioff >= 0 && ioff + csr_bsize <= bi) &&
        (joff >= 0 && joff + csr_bsize <= bj)) {
      for (int m = 0; m < csr_bsize; m++) {
        for (int n = 0; n < csr_bsize; n++) {
          A[(ioff + m) + bi * (joff + n)] += a[csr_bsize * m + n];
        }
      }

      return 1;
    }
  } else {
    fprintf(stderr,
            "[%d] TACSBlockCyclicMat: Error, (%d, %d) not in nz-pattern\n",
            rank, i, j);
  }

  return 0;
}

/*
  Assign randomly generated entries to the matrix.

  Note that this uses the rand() function from stdlib. The matrix
  entries lie within the unit interval.
*/
void TACSBlockCyclicMat::setRand() {
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Fill the matrix with randomly generated entries
  for (int i = 0; i < nrows; i++) {
    int bi = bptr[i + 1] - bptr[i];
    if (rank == get_block_owner(i, i)) {
      int np = dval_offset[i];
      TacsScalar *a = &Dvals[np];
      for (int k = 0; k < bi * bi; k++) {
        a[k] = (1.0 * rand()) / RAND_MAX;
      }
    } else {
      for (int k = 0; k < bi * bi; k++) {
        rand();
      }
    }
  }

  // Compute the upper triangular part
  for (int i = 0; i < nrows; i++) {
    int bi = bptr[i + 1] - bptr[i];
    for (int jp = Urowp[i]; jp < Urowp[i + 1]; jp++) {
      int j = Ucols[jp];
      int bj = bptr[j + 1] - bptr[j];

      if (rank == get_block_owner(i, j)) {
        int np = uval_offset[jp];
        TacsScalar *a = &Uvals[np];
        for (int k = 0; k < bi * bj; k++) {
          a[k] = (1.0 * rand()) / RAND_MAX;
        }
      } else {
        for (int k = 0; k < bi * bj; k++) {
          rand();
        }
      }
    }
  }

  // Compute the lower triangular part
  for (int j = 0; j < ncols; j++) {
    int bj = bptr[j + 1] - bptr[j];
    for (int ip = Lcolp[j]; ip < Lcolp[j + 1]; ip++) {
      int i = Lrows[ip];
      int bi = bptr[i + 1] - bptr[i];

      if (rank == get_block_owner(i, j)) {
        int np = lval_offset[ip];
        TacsScalar *a = &Lvals[np];

        for (int k = 0; k < bi * bj; k++) {
          a[k] = (1.0 * rand()) / RAND_MAX;
        }
      } else {
        for (int k = 0; k < bi * bj; k++) {
          rand();
        }
      }
    }
  }
}

/*
  Multiply y = A*x

  Note that this code assumes that x/y are on all processors.  This
  function will not work after the matrix is factored since the
  factorization over-writes the original matrix entries.
*/
void TACSBlockCyclicMat::mult(TacsScalar *x, TacsScalar *y) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Get the location of rank on the process grid
  int proc_row = -1, proc_col = -1;
  get_proc_row_column(rank, &proc_row, &proc_col);

  memset(y, 0, xbptr[nrows] * sizeof(TacsScalar));

  // Allocate temporary space (if needed)
  TacsScalar *tx = new TacsScalar[cbptr[ncols]];
  TacsScalar *ty = new TacsScalar[rbptr[nrows]];
  memset(ty, 0, rbptr[nrows] * sizeof(TacsScalar));

  if (proc_row >= 0) {
    // Allocate an array of requests
    MPI_Request *sreq = new MPI_Request[nprows];

    // Compute the on-process parts for the diagonal contributions
    for (int j = 0; j < ncols; j++) {
      // Compute the block size
      int bj = bptr[j + 1] - bptr[j];

      if (rank == get_block_owner(j, j)) {
        // Send the x-values to the other processors in this process
        // column
        int nj = xbptr[j];

        // Send the vector components to the other processors in this
        // process column
        for (int p = 0, k = 0; p < nprows; p++) {
          int dest = proc_grid[proc_col + p * npcols];
          if (dest != rank) {
            MPI_Isend(&x[nj], bj, TACS_MPI_TYPE, dest, p, comm, &sreq[k]);
            k++;
          } else {
            memcpy(&tx[cbptr[j]], &x[nj], bj * sizeof(TacsScalar));
          }
        }

        if (j < nrows) {
          // Multiply by the diagonal part
          int np = dval_offset[j];
          TacsScalar alpha = 1.0, beta = 1.0;
          int one = 1;
          BLASgemv("N", &bj, &bj, &alpha, &Dvals[np], &bj, &x[nj], &one, &beta,
                   &y[nj], &one);
          TacsAddFlops(2 * bj * bj);

          // Compute the lower triangular part
          for (int ip = Lcolp[j]; ip < Lcolp[j + 1]; ip++) {
            int i = Lrows[ip];

            if (rank == get_block_owner(i, j)) {
              // Compute the size of the block
              int bi = bptr[i + 1] - bptr[i];

              // Compute where this should go in the row-space
              int ni = rbptr[i];

              int np = lval_offset[ip];
              TacsScalar alpha = 1.0, beta = 1.0;
              int one = 1;
              BLASgemv("N", &bi, &bj, &alpha, &Lvals[np], &bi, &x[nj], &one,
                       &beta, &ty[ni], &one);
              TacsAddFlops(2 * bi * bj);
            }
          }
        }

        // Wait for all the sends to complete
        MPI_Waitall(nprows - 1, sreq, MPI_STATUSES_IGNORE);
      } else if (proc_col == get_proc_column(j)) {
        // Recv the column values
        int source = get_block_owner(j, j);

        // Get the column owner
        int nj = cbptr[j];

        // Recv and store the column value
        MPI_Recv(&tx[nj], bj, TACS_MPI_TYPE, source, proc_row, comm,
                 MPI_STATUS_IGNORE);

        // Compute the lower triangular part
        for (int ip = Lcolp[j]; ip < Lcolp[j + 1]; ip++) {
          int i = Lrows[ip];

          if (rank == get_block_owner(i, j)) {
            // Compute the size of the block
            int bi = bptr[i + 1] - bptr[i];

            // Compute where this should go in the row-space
            int ni = rbptr[i];

            // Compute the product
            int np = lval_offset[ip];
            TacsScalar alpha = 1.0, beta = 1.0;
            int one = 1;
            BLASgemv("N", &bi, &bj, &alpha, &Lvals[np], &bi, &tx[nj], &one,
                     &beta, &ty[ni], &one);
            TacsAddFlops(2 * bi * bj);
          }
        }
      }
    }

    // Free the send request array
    delete[] sreq;

    // Loop over the rows in the matrix
    for (int i = 0; i < nrows; i++) {
      int bi = bptr[i + 1] - bptr[i];
      int ni = rbptr[i];

      // Compute the product from this row
      if (proc_row == get_proc_row(i)) {
        for (int jp = Urowp[i]; jp < Urowp[i + 1]; jp++) {
          int j = Ucols[jp];
          int bj = bptr[j + 1] - bptr[j];

          if (rank == get_block_owner(i, j)) {
            // Set the pointer to the correct components of x
            TacsScalar *xp = NULL;
            if (i == j) {
              int nj = xbptr[j];
              xp = &x[nj];
            } else {
              int nj = cbptr[j];
              xp = &tx[nj];
            }

            // Compute the product
            int np = uval_offset[jp];
            TacsScalar alpha = 1.0, beta = 1.0;
            int one = 1;
            BLASgemv("N", &bi, &bj, &alpha, &Uvals[np], &bi, xp, &one, &beta,
                     &ty[ni], &one);
            TacsAddFlops(2 * bi * bj);
          }
        }

        // Send the product from this processor from the current
        // iteration
        int owner = get_block_owner(i, i);
        if (owner != rank) {
          int tag = proc_col;
          MPI_Send(&ty[ni], bi, TACS_MPI_TYPE, owner, tag, comm);
        } else if (owner == rank) {
          // Add up the contributions from this processor
          int nj = xbptr[i];
          for (int j = 0; j < bi; j++) {
            y[nj + j] += ty[ni + j];
          }

          for (int k = 0; k < npcols; k++) {
            int source = get_block_owner(proc_row, k);
            if (source != owner) {
              int tag = k;
              MPI_Recv(&ty[ni], bi, TACS_MPI_TYPE, source, tag, comm,
                       MPI_STATUS_IGNORE);

              for (int j = 0; j < bi; j++) {
                y[nj + j] += ty[ni + j];
              }
            }
          }
        }
      }
    }

    // Free the temporary data
    delete[] tx;
    delete[] ty;
  }

  MPI_Barrier(comm);
}

/*
  Based on the non-zero structure of the LU-factorization, determine
  in advance the sizes of the data structure required for the
  back-solves.

  This function allocates the following data:

  lower_row_sum_recv: the number of row sum recvs for L^{-1} for row i
  upper_row_sum_recv: the number of row sum recvs for U^{-1} for row i

  lower_row_sum_count: the number of local updates for L^{-1} for row i
  upper_row_sum_count: the number of local updates for U^{-1} for row i
*/
void TACSBlockCyclicMat::init_row_counts() {
  int rank;
  MPI_Comm_rank(comm, &rank);

  lower_block_count = 0;
  upper_block_count = 0;

  lower_row_sum_count = NULL;
  upper_row_sum_count = NULL;
  lower_row_sum_recv = NULL;
  upper_row_sum_recv = NULL;

  // Get the location of rank on the process grid
  int proc_row = -1, proc_col = -1;

  // Check if this process is on the process grid - if not, then
  // it does not participate in the computation
  if (get_proc_row_column(rank, &proc_row, &proc_col)) {
    // Compute the expected number of recvs for this processor
    int *row_sum_recv = new int[nrows * npcols];
    memset(row_sum_recv, 0, nrows * npcols * sizeof(int));

    // Compute the number of local updates for row i in L^{-1}
    lower_row_sum_count = new int[nrows];
    memset(lower_row_sum_count, 0, nrows * sizeof(int));
    for (int j = 0; j < ncols; j++) {
      for (int ip = Lcolp[j]; ip < Lcolp[j + 1]; ip++) {
        int i = Lrows[ip];
        if (i < nrows) {
          int owner = get_block_owner(i, j);
          if (rank == owner) {
            lower_row_sum_count[i]++;
          }

          // If this is the diagonal owner but not the
          // owner of L[i,j], then this will contribute a
          // row sum to the proc. 'owner'
          int diag_owner = get_block_owner(i, i);
          if (rank == diag_owner && rank != owner) {
            row_sum_recv[npcols * i + (j % npcols)] = 1;
          }
        }
      }
    }

    // Count up the number of recvs from each row
    lower_row_sum_recv = new int[nrows];
    memset(lower_row_sum_recv, 0, nrows * sizeof(int));
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < npcols; j++) {
        lower_row_sum_recv[i] += row_sum_recv[i * npcols + j];
      }
      if (rank == get_block_owner(proc_row, i) &&
          rank != get_block_owner(i, i)) {
        lower_block_count++;
      }
      if (rank == get_block_owner(i, i)) {
        lower_block_count += lower_row_sum_recv[i];
      }
    }

    // Reset the row sum recv count to zero
    memset(row_sum_recv, 0, nrows * npcols * sizeof(int));

    // Compute the number of for updates for row i in U^{-1}
    upper_row_sum_count = new int[nrows];
    memset(upper_row_sum_count, 0, nrows * sizeof(int));
    for (int i = 0; i < nrows; i++) {
      for (int jp = Urowp[i]; jp < Urowp[i + 1]; jp++) {
        int j = Ucols[jp];
        int owner = get_block_owner(i, j);
        if (rank == owner) {
          upper_row_sum_count[i]++;
        }

        // If this is the diagonal owner but not the
        // owner of L[i,j], then this will contribute a
        // row sum to the proc. 'owner'
        int diag_owner = get_block_owner(i, i);
        if (rank == diag_owner && rank != owner) {
          row_sum_recv[npcols * i + (j % npcols)] = 1;
        }
      }
    }

    // Count up the number of recvs from each row
    upper_row_sum_recv = new int[nrows];
    memset(upper_row_sum_recv, 0, nrows * sizeof(int));
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < npcols; j++) {
        upper_row_sum_recv[i] += row_sum_recv[i * npcols + j];
      }
      if (rank == get_block_owner(proc_row, i) &&
          rank != get_block_owner(i, i)) {
        upper_block_count++;
      }
      upper_block_count += upper_row_sum_recv[i];
    }

    delete[] row_sum_recv;
  }
}

/*
  Apply the LU factorization to x

  x <-- U^{-1} L^{-1} x

  The back-solve is performed in place. No check is performed to
  ensure that you've factored the matrix first, so beware!

  Procedure:

  This is an event-driven implementation of the back-solve that should
  work better in parallel than a row-wise or a coordinated
  variant. Note that the back-solve is more challenging to implement
  in parallel because the communication to computation ratio is
  higher. The event-driven implementation works in two stages.  First,
  the L^{-1} application, then the U^{-1} application. Both work in a
  similar manner.

  Algorithm for the application of L^{-1}:

  Initialize:

  for i in range(0, nrows):
  .   if lower_row_sum_count[i] == 0 and lower_row_sum_recv[i] == 0:
  .       compute x[i] = L[i,i]^{-1}x[i]
  .       send x[i] to processors in column i

  Until completed do:

  1. Receive x[i] and add the product with column xsum[:] += L[:,i]*x[i]
  for j in Lrows[Lcolp[i]:Lcolp[i+1]]:
  .   xsum[j] += L[j,i]*x[i]
  .   row_sum_count[j] += 1.
  .   if row_sum_count[j] == lower_row_sum_count[j]
  .       send xsum[j] to the diagonal owner of L[j,j]

  2. Receive xsum for row i and add to x[i]
  x[i] += xsum
  row_sum_recv[i] += 1
  If (row_sum_recv[i] == lower_row_sum_recv[i] and
  .   row_sum_count[i] == lower_row_sum_count[i]):
  .     compute x[i] = L[i,i]^{-1}(x[i] - xsum[i]),
  .     send x[i] to the j-th columns
*/
void TACSBlockCyclicMat::applyFactor(TacsScalar *x) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Space required for the factorization
  TacsScalar *tx = NULL;
  TacsScalar *xsum = NULL;
  TacsScalar *xlocal = NULL;

  // Row sum and recv counters
  int *row_sum_count = NULL;
  int *row_sum_recv = NULL;

  // Get the location of rank on the process grid
  int proc_row = -1, proc_col = -1;
  get_proc_row_column(rank, &proc_row, &proc_col);

  // Check if this process is on the process grid - if not, then
  // it does not participate in the computation
  if (proc_row >= 0) {
    // The temporary array used to store in-coming values
    xlocal = new TacsScalar[max_bsize];

    // Temporary vector of x-values on this processor
    tx = new TacsScalar[cbptr[nrows]];
    memset(tx, 0, cbptr[nrows] * sizeof(TacsScalar));

    // The locally stored sum of the products of L[i:0:i]*x[0:i]
    xsum = new TacsScalar[rbptr[nrows]];
    memset(xsum, 0, rbptr[nrows] * sizeof(TacsScalar));

    // Allocate space for the counters
    row_sum_count = new int[nrows];
    row_sum_recv = new int[nrows];
    memset(row_sum_count, 0, nrows * sizeof(int));
    memset(row_sum_recv, 0, nrows * sizeof(int));

    // Initiate the back-solve by initiating all factorizations with
    // a row-sum of zero.
    for (int row = 0; row < nrows; row++) {
      if (lower_row_sum_count[row] == 0 && lower_row_sum_recv[row] == 0 &&
          rank == get_block_owner(row, row)) {
        // Finish this entry
        add_lower_row_sum(row, x, tx, xsum, row_sum_count, row_sum_recv);
      }
    }

    for (int k = 0; k < lower_block_count; k++) {
      // Probe the message to see what size it is
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);

      // Even tags indicate a new value of x[col]
      // Odd tags indicate a newly completed row sum
      if (status.MPI_TAG % 2 == 0) {
        // The message contains the new value of x[col]
        int col = status.MPI_TAG / 2;

        // Find the local block size and location
        int ni = cbptr[col];
        int bi = bptr[col + 1] - bptr[col];

        // Receieve x[ni]
        MPI_Recv(&tx[ni], bi, TACS_MPI_TYPE, status.MPI_SOURCE, status.MPI_TAG,
                 comm, &status);

        // Compute the update to the column
        lower_column_update(col, x, tx, xsum, row_sum_count, row_sum_recv);
      } else {
        // The message contains the completed row sum: xsum[row]
        int row = (status.MPI_TAG - 1) / 2;

        // Find the local block size and location
        int ni = rbptr[row];
        int bi = bptr[row + 1] - bptr[row];

        // Receieve xsum[i] into the buffer xlocal
        MPI_Recv(xlocal, bi, TACS_MPI_TYPE, status.MPI_SOURCE, status.MPI_TAG,
                 comm, &status);

        // Add the contribution from the non-local block
        TacsScalar alpha = 1.0;
        int one = 1;
        BLASaxpy(&bi, &alpha, xlocal, &one, &xsum[ni], &one);
        TacsAddFlops(2 * bi);

        // Increment the number of received rows and the block count
        row_sum_recv[row]++;

        // Finalize the entry if all sums have been computed
        if (row_sum_recv[row] == lower_row_sum_recv[row] &&
            row_sum_count[row] == lower_row_sum_count[row]) {
          add_lower_row_sum(row, x, tx, xsum, row_sum_count, row_sum_recv);
        }
      }
    }
  }

  // Barrier required to prevent interference between
  // communication since we're using an event-driven
  // scheme

  MPI_Barrier(comm);

  if (proc_row >= 0) {
    // Set the partial row sums to zero again
    memset(xsum, 0, rbptr[nrows] * sizeof(TacsScalar));

    // Now, apply the upper factorization
    memset(row_sum_count, 0, nrows * sizeof(int));
    memset(row_sum_recv, 0, nrows * sizeof(int));

    for (int row = nrows - 1; row >= 0; row--) {
      if (upper_row_sum_count[row] == 0 && upper_row_sum_recv[row] == 0 &&
          rank == get_block_owner(row, row)) {
        // Entry can be completed
        add_upper_row_sum(row, x, tx, xsum, row_sum_count, row_sum_recv);
      }
    }

    for (int k = 0; k < upper_block_count; k++) {
      // Receieve the next task to perform
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);

      // Even tags indicate a new value of x[col]
      // Odd tags indicate a newly completed row sum
      if (status.MPI_TAG % 2 == 0) {
        // The message contains the new value of x[col]
        int col = status.MPI_TAG / 2;

        // Find the local block size and location
        int ni = cbptr[col];
        int bi = bptr[col + 1] - bptr[col];

        // Receieve x[ni]
        MPI_Recv(&tx[ni], bi, TACS_MPI_TYPE, status.MPI_SOURCE, status.MPI_TAG,
                 comm, &status);

        // Add the update to the column
        upper_column_update(col, x, tx, xsum, row_sum_count, row_sum_recv);
      } else {
        // The message contains the completed row sum: xsum[row]
        int row = (status.MPI_TAG - 1) / 2;

        // Find the local block size and location
        int ni = rbptr[row];
        int bi = bptr[row + 1] - bptr[row];

        // Receieve xsum[i] into the buffer xlocal
        MPI_Recv(xlocal, bi, TACS_MPI_TYPE, status.MPI_SOURCE, status.MPI_TAG,
                 comm, &status);

        // Add the contributions from the non-local blocks
        TacsScalar alpha = 1.0;
        int one = 1;
        BLASaxpy(&bi, &alpha, xlocal, &one, &xsum[ni], &one);
        TacsAddFlops(2 * bi);

        // Increment the number of receieved rows
        row_sum_recv[row]++;

        // Finalize the entry if all sums have been computed
        if (row_sum_recv[row] == upper_row_sum_recv[row] &&
            row_sum_count[row] == upper_row_sum_count[row]) {
          add_upper_row_sum(row, x, tx, xsum, row_sum_count, row_sum_recv);
        }
      }
    }

    delete[] xlocal;
    delete[] tx;
    delete[] xsum;
    delete[] row_sum_count;
    delete[] row_sum_recv;
  }

  // Barrier required to prevent interference between U^{-1}
  // application and the re-distribution of the vector components
  MPI_Barrier(comm);
}

/*
  Apply the update to the columns of the lower-diagonal matrix.

  Note that this is a doubly-recursive function with
  add_lower_row_sum.

  Given the factored matrix S = L*U, compute

  for i = col+1:n
  .  xsum[i] += L[i,col]*x[col]
  .  row_sum_count[i] += 1
  .  if (row_sum_count[i] == lower_row_sum_count[i] and
  .      row_sum_recv[i] == lower_row_sum_recv[i])
  .      send xsum[i] to proc_owner(i,i)

  Note that if the destination is the local processor,
  add_lower_row_sum is called immediately.

  input:
  col:           the column to update
  x:             the solution vector on all procs
  xsum:          the column sums on this processor
  row_sum_count: the number of updates required for row[i] from L[i,0:i]
  row_sum_recv:  the number of received row sums
*/
void TACSBlockCyclicMat::lower_column_update(int col, TacsScalar *x,
                                             TacsScalar *tx, TacsScalar *xsum,
                                             int *row_sum_count,
                                             int *row_sum_recv) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Keep track of the send requests
  MPI_Request req_list[BACKSOLVE_BUFF_SIZE];
  int nreq = 0;

  // Compute the location/number of columns in the block
  int dj = xbptr[col];                 // Variable location in x
  int nj = cbptr[col];                 // Variable location in tx
  int bj = bptr[col + 1] - bptr[col];  // block size

  // Loop over the column and compute xsum[j] += L[j,col]*x[col]
  for (int jp = Lcolp[col]; jp < Lcolp[col + 1]; jp++) {
    int row = Lrows[jp];

    // If this is the block owner, multiply by L[row, col]*x[col]
    if (rank == get_block_owner(row, col)) {
      int ni = rbptr[row];
      int bi = bptr[row + 1] - bptr[row];
      TacsScalar *L = &Lvals[lval_offset[jp]];

      // Set the pointer for where to take xp from...
      TacsScalar *xp = &tx[nj];
      if (rank == get_block_owner(col, col)) {
        xp = &x[dj];
      }

      TacsScalar alpha = 1.0, beta = 1.0;
      int one = 1;
      // xsum[i] = xsum[i] + L[i,j]*x[j] for
      BLASgemv("N", &bi, &bj, &alpha, L, &bi, xp, &one, &beta, &xsum[ni], &one);
      TacsAddFlops(2 * bi * bj);

      // Update row_sum_count[row]
      row_sum_count[row]++;

      // The local summation for this row is complete
      if (row_sum_count[row] == lower_row_sum_count[row]) {
        // Send the result to the diagonal owner
        int dest = get_block_owner(row, row);

        if (rank != dest) {
          // If the number of requests exceeds the buffer size,
          // wait for them to complete
          if (nreq >= BACKSOLVE_BUFF_SIZE) {
            MPI_Waitall(nreq, req_list, MPI_STATUSES_IGNORE);
            nreq = 0;
          }

          // Send the row sum to the other processor
          MPI_Isend(&xsum[ni], bi, TACS_MPI_TYPE, dest, 2 * row + 1, comm,
                    &req_list[nreq]);
          nreq++;
        } else if (rank == dest &&
                   row_sum_recv[row] == lower_row_sum_recv[row]) {
          // If all remaining blocks have been received, finish this entry
          add_lower_row_sum(row, x, tx, xsum, row_sum_count, row_sum_recv);
        }
      }
    }
  }

  // Wait for the pending communication requests
  MPI_Waitall(nreq, req_list, MPI_STATUSES_IGNORE);
}

/*
  Add the row sum from the lower portion of the factorization.

  This is a doubly-recursive function with lower_column_update.

  Add the row sum stored in xlocal to xsum. If the total number of
  updates to the row == the number of required updates

  input:
  col:           the integer for the column
  x:             the solution vector
  xsum:          the sum of all rows
  xlocal:        a temporary array
  row_sum_count: the number of updates required for the row[i]
  row_sum_recv:  the number of received row sums
*/
void TACSBlockCyclicMat::add_lower_row_sum(int col, TacsScalar *x,
                                           TacsScalar *tx, TacsScalar *xsum,
                                           int *row_sum_count,
                                           int *row_sum_recv) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Keep track of the number of the send requests
  MPI_Request req_list[BACKSOLVE_COLUMN_SIZE];
  int nreq = 0;

  // Find the local block size and location
  int di = xbptr[col];  // Variable location in x
  int ni = rbptr[col];  // Variable location in xsum
  int bi = bptr[col + 1] - bptr[col];

  // Add the contributions from the non-local blocks
  // x[i] <- x[i] - xsum[i]
  TacsScalar alpha = -1.0;
  int one = 1;
  BLASaxpy(&bi, &alpha, &xsum[ni], &one, &x[di], &one);
  TacsAddFlops(2 * bi);

  // Send the result to the processors in this column
  for (int ii = 0; ii < nprows; ii++) {
    int dest = get_block_owner(ii, col);

    // Send the result to the column processors
    if (dest != rank) {
      // If the number of requests exceeds the buffer size,
      // wait for them to complete
      if (nreq >= BACKSOLVE_COLUMN_SIZE) {
        MPI_Waitall(nreq, req_list, MPI_STATUSES_IGNORE);
        nreq = 0;
      }

      // Send x[i] to all the processes in the column
      MPI_Isend(&x[di], bi, TACS_MPI_TYPE, dest, 2 * col, comm,
                &req_list[nreq]);
      nreq++;
    }
  }

  // Now, initiate the block-diagonal solve
  lower_column_update(col, x, tx, xsum, row_sum_count, row_sum_recv);

  // Wait for the pending communication requests
  MPI_Waitall(nreq, req_list, MPI_STATUSES_IGNORE);
}

/*
  Apply the update to the columns of the upper-diagonal matrix.

  Note that this is a doubly-recursive function with
  add_upper_row_sum.

  Given the factored matrix S = L*U, compute

  for i = col+1:n
  .  xsum[i] += L[i,col]*x[col]
  .  row_sum_count[i] -= 1
  .  if (row_sum_count[i] == 0)
  .      send xsum[i] to proc_owner(i,i)

  Note that if the destination is the local processor,
  add_upper_row_sum is called immediately.

  input:
  col:           the column to update
  x:             the solution vector on all procs
  tx:            temporary local x-values
  xsum:          the column sums on this processor
  row_sum_count: the number of updates required for row[i] from L[i,0:i]
*/
void TACSBlockCyclicMat::upper_column_update(int col, TacsScalar *x,
                                             TacsScalar *tx, TacsScalar *xsum,
                                             int *row_sum_count,
                                             int *row_sum_recv) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Keep track of the number of buffers
  MPI_Request req_list[BACKSOLVE_BUFF_SIZE];
  int nreq = 0;

  // Compute the location/number of columns in the block
  int dj = xbptr[col];                 // Variable location in x
  int nj = cbptr[col];                 // Variable location in tx
  int bj = bptr[col + 1] - bptr[col];  // block size

  for (int row = col - 1; row >= 0; row--) {
    if (rank == get_block_owner(row, col)) {
      // Find the location where Ucols[jp] == col
      int jp = Urowp[row];
      int size = Urowp[row + 1] - jp;
      int *item = TacsSearchArray(col, size, &Ucols[jp]);

      // Check if the column is in the array
      if (item) {
        jp += item - &Ucols[jp];

        // Find the row-sum information
        int ni = rbptr[row];
        int bi = bptr[row + 1] - bptr[row];
        TacsScalar *U = &Uvals[uval_offset[jp]];

        // Determine whether this is locally owned or not
        TacsScalar *xp = &tx[nj];
        if (rank == get_block_owner(col, col)) {
          xp = &x[dj];
        }

        TacsScalar alpha = 1.0, beta = 1.0;
        int one = 1;
        // xsum[i] <-- xsum[i] + U[i,j]*x[j] for
        BLASgemv("N", &bi, &bj, &alpha, U, &bi, xp, &one, &beta, &xsum[ni],
                 &one);
        TacsAddFlops(2 * bi * bj);

        // Update row_sum_count[row]
        row_sum_count[row]++;

        if (row_sum_count[row] == upper_row_sum_count[row]) {
          // Send the result to the diagonal owner
          int dest = get_block_owner(row, row);

          if (rank != dest) {
            // If the number of requests exceeds the buffer size,
            // wait for them to complete
            if (nreq >= BACKSOLVE_BUFF_SIZE) {
              MPI_Waitall(nreq, req_list, MPI_STATUSES_IGNORE);
              nreq = 0;
            }

            // Send the row sum to processors in the row
            MPI_Isend(&xsum[ni], bi, TACS_MPI_TYPE, dest, 2 * row + 1, comm,
                      &req_list[nreq]);
            nreq++;
          } else if (rank == dest &&
                     row_sum_recv[row] == upper_row_sum_recv[row]) {
            // Update this row
            add_upper_row_sum(row, x, tx, xsum, row_sum_count, row_sum_recv);
          }
        }
      }
    }
  }

  // Wait the pending sends
  MPI_Waitall(nreq, req_list, MPI_STATUSES_IGNORE);
}

/*
  Add the row sum from the upper portion of the factorization.

  This is a doubly-recursive function with upper_column_update.

  Add the row sum stored in xlocal to xsum. If the total number of
  updates to the row == the number of required updates, perform

  input:
  row:           the integer for the row sum
  x:             the solution vector
  tx:            a temporary array of x-local values
  xsum:          the sum of all rows
  row_sum_count: the number of updates required for the row[i]
*/
void TACSBlockCyclicMat::add_upper_row_sum(int row, TacsScalar *x,
                                           TacsScalar *tx, TacsScalar *xsum,
                                           int *row_sum_count,
                                           int *row_sum_recv) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Keep track of the number of the send requests
  MPI_Request req_list[BACKSOLVE_COLUMN_SIZE];
  int nreq = 0;

  // Find the local block size and location
  int di = xbptr[row];
  int ni = rbptr[row];
  int bi = bptr[row + 1] - bptr[row];

  // Add the contributions from the non-local blocks
  // xsum[i] <- xsum[i] - x[i]
  TacsScalar alpha = -1.0;
  int one = 1;
  BLASaxpy(&bi, &alpha, &x[di], &one, &xsum[ni], &one);
  TacsAddFlops(2 * bi);

  // Assign x <- -D*(xsum - x)
  int np = dval_offset[row];
  alpha = -1.0;
  TacsScalar beta = 0.0;
  BLASgemv("N", &bi, &bi, &alpha, &Dvals[np], &bi, &xsum[ni], &one, &beta,
           &x[di], &one);
  TacsAddFlops(2 * bi * bi);

  // Send the result to the processors in this column
  for (int ii = 0; ii < nprows; ii++) {
    int dest = get_block_owner(ii, row);

    if (dest != rank) {
      // If the number of requests exceeds the buffer size,
      // wait for them to complete
      if (nreq >= BACKSOLVE_COLUMN_SIZE) {
        MPI_Waitall(nreq, req_list, MPI_STATUSES_IGNORE);
        nreq = 0;
      }

      // Send x[i] to all processors in this column
      MPI_Isend(&x[di], bi, TACS_MPI_TYPE, dest, 2 * row, comm,
                &req_list[nreq]);
      nreq++;
    }
  }

  // Perform the column update
  upper_column_update(row, x, tx, xsum, row_sum_count, row_sum_recv);

  // Wait the pending communication requests
  MPI_Waitall(nreq, req_list, MPI_STATUSES_IGNORE);
}

/*
  Retrieve the pointer to the (i,j) block. If the block is not owned
  by this process or if the block is not in the non-zero pattern,
  return NULL.

  This function first checks whether the block is in the L/D/U blocks.
  Then, it uses either the column or row indices to locate the row or
  column it is in. These searches rely on the fact that the
  rows/column indices are sorted.
*/
TacsScalar *TACSBlockCyclicMat::get_block(int rank, int i, int j) {
  TacsScalar *A = NULL;

  if (rank == get_block_owner(i, j)) {
    if (i > j) {  // L
      int lp = Lcolp[j];
      int size = Lcolp[j + 1] - Lcolp[j];

      // Look for row i
      int *item = TacsSearchArray(i, size, &Lrows[lp]);
      if (item) {
        lp = lp + (item - &Lrows[lp]);
        A = &Lvals[lval_offset[lp]];
      }
    } else if (i == j) {  // D
      A = &Dvals[dval_offset[i]];
    } else {  // i < j : U
      int up = Urowp[i];
      int size = Urowp[i + 1] - Urowp[i];

      // Look for column j
      int *item = TacsSearchArray(j, size, &Ucols[up]);
      if (item) {
        up = up + (item - &Ucols[up]);
        A = &Uvals[uval_offset[up]];
      }
    }
  }

  return A;
}

/*
  Factor the matrix in-place, in parallel.

  This code is intended for use only if the matrix is nearly dense and
  the block sizes are chosen sufficiently large. The optimal block
  size will depend on the computer architecture. Some test runs may be
  required to tune the code.

  This code computes an L/U matrix factorization of the almost dense
  matrix.

  The L and U factors are defined such that,

  [ A[i,i]      |  A[i+1,i:n]     ]
  [ A[i+1:n,i]  |  A[i+1:n,i+1:n] ]
  =
  [ 1           |   0             ][ U[i,i]  |  U[i,i+1:n]     ]
  [ L[i+1:n,i]  |  L[i+1:n,i+1:n] ][ 0       |  U[i+1:n,i+1:n] ]

  As a result, the following relationships hold:

  U[i,i] = A[i,i]
  U[i,i+1:n] = A[i,i+1:n]
  L[i+1:n,i] = A[i+1:n,i]*U[i,i]^{-1}

  With the rank-bi update:

  A[i+1:n,i+1:n] <-- A[i+1:n,i+1:n] - L[+1i:n,i]*U[i,i+1:n]

  The computational steps involved in the factorization are:

  For i = 1,...,n

  1. Compute the diagonal block factorization of A[i,i]
  2. Compute L[i+1:n,i] = A[i+1:n,i]*U[i,i]^{-1}
  3. Compute the update
  A[i+1:n,i+1:n] <-- A[i+1:n,i+1:n] - L[i+1:n,i]*U[i,i+1:n]
*/
void TACSBlockCyclicMat::factor() {
  int rank;
  MPI_Comm_rank(comm, &rank);

  int proc_row, proc_col;  // Get the location of rank on the process grid
  if (!get_proc_row_column(rank, &proc_row, &proc_col)) {
    // This process is not on the process grid - does not participate
    return;
  }

  int *temp_piv = new int[max_bsize];
  TacsScalar *temp_diag = new TacsScalar[max_bsize * max_bsize];
  TacsScalar *temp_block = new TacsScalar[max_bsize * max_bsize];
  int lwork = 128 * max_bsize;
  TacsScalar *work = new TacsScalar[lwork];

  // Buffers to handle the recieves information
  TacsScalar *Ubuff = new TacsScalar[max_ubuff_size];
  TacsScalar *Lbuff = new TacsScalar[max_lbuff_size];

  // Send information for rows owning U
  MPI_Request *U_send_request = new MPI_Request[nprows - 1];
  MPI_Status *U_send_status = new MPI_Status[nprows - 1];

  // Send information for columns owning L
  MPI_Request *L_send_request = new MPI_Request[npcols - 1];
  MPI_Status *L_send_status = new MPI_Status[npcols - 1];

  double t_update = 0.0;
  double t_recv_wait = 0.0;
  double t_send_wait = 0.0;
  int n_gemm = 0;

  for (int i = 0; i < nrows; i++) {
    int bi = bptr[i + 1] - bptr[i];

    // The diagonal factor of A and its pivot
    TacsScalar *d_diag = NULL;
    int diag_owner = get_block_owner(i, i);

    // Get the owner for the diagonal block
    if (rank == diag_owner) {
      // Determine the address of the diagonal block
      int nd = dval_offset[i];
      d_diag = &Dvals[nd];

      // Compute the inverse of the diagonal block
      int info;
      LAPACKgetrf(&bi, &bi, d_diag, &bi, temp_piv, &info);
      LAPACKgetri(&bi, d_diag, &bi, temp_piv, work, &lwork, &info);
      // Add flops from the inversion
      TacsAddFlops(1.333333 * bi * bi * bi);

      // Send the factor to the column processes
      for (int p = 0; p < nprows; p++) {
        int dest = proc_grid[proc_col + p * npcols];
        if (rank != dest) {
          MPI_Send(d_diag, bi * bi, TACS_MPI_TYPE, dest, p, comm);
        }
      }
    }

    // Receive U[i,i]^{-1}
    if (rank != diag_owner && proc_col == get_proc_column(i)) {
      MPI_Status status;
      MPI_Recv(temp_diag, bi * bi, TACS_MPI_TYPE, diag_owner, proc_row, comm,
               &status);
      d_diag = temp_diag;
    }

    // Before computing the entries for L[i:n,i], post all the recieve
    // information for U and L.

    // Determine the size of the incoming/outgoing U
    int ubuff_size = 0;
    for (int jp = Urowp[i]; jp < Urowp[i + 1]; jp++) {
      int j = Ucols[jp];
      int bj = bptr[j + 1] - bptr[j];

      if (get_proc_column(j) == proc_col) {
        ubuff_size += bi * bj;
      }
    }

    // Set the U values to the row processes that need it
    MPI_Request U_recv_request;
    int source_proc_row = get_proc_row(i);
    if (source_proc_row == proc_row) {
      // The sending processes
      int offset = uval_offset[Urowp[i]];
      for (int p = 0, k = 0; p < nprows; p++) {
        int dest = proc_grid[proc_col + p * npcols];
        if (rank != dest) {
          int tag = 2 * i;
          MPI_Isend(&Uvals[offset], ubuff_size, TACS_MPI_TYPE, dest, tag, comm,
                    &U_send_request[k]);
          k++;
        }
      }
    } else {
      // The receiving processes
      int source = proc_grid[proc_col + source_proc_row * npcols];
      int tag = 2 * i;
      MPI_Irecv(Ubuff, ubuff_size, TACS_MPI_TYPE, source, tag, comm,
                &U_recv_request);
    }

    // Determine the size of the incoming/outgoing L
    int lbuff_size = 0;
    for (int jp = Lcolp[i]; jp < Lcolp[i + 1]; jp++) {
      int j = Lrows[jp];
      int bj = bptr[j + 1] - bptr[j];

      if (get_proc_row(j) == proc_row) {
        lbuff_size += bi * bj;
      }
    }

    // Compute L[i+1:n,i] = A[i+1:n,i]*U[i,i]^{-1} and send the results to
    // the destination immediately
    if (proc_col == get_proc_column(i)) {
      for (int jp = Lcolp[i]; jp < Lcolp[i + 1]; jp++) {
        int j = Lrows[jp];
        int bj = bptr[j + 1] - bptr[j];

        if (rank == get_block_owner(j, i)) {
          int np = lval_offset[jp];

          // Compute L[i+1:n,i] = A[i+1:n,i]*U[i,i]^{-1}
          // L in bj x bi
          // A in bj x bi
          // U^{-1} in bi x bi
          TacsScalar alpha = 1.0, beta = 0.0;
          BLASgemm("N", "N", &bj, &bi, &bi, &alpha, &Lvals[np], &bj, d_diag,
                   &bi, &beta, temp_block, &bj);
          n_gemm++;
          TacsAddFlops(2 * bi * bi * bj);

          memcpy(&Lvals[np], temp_block, bi * bj * sizeof(TacsScalar));
        }
      }
    }

    // Set the L values to the row processes that need it
    MPI_Request L_recv_request;
    int source_proc_column = get_proc_column(i);
    if (source_proc_column == proc_col) {
      // Send the proc_rows requiring everything
      for (int p = 0, k = 0; p < npcols; p++) {
        int dest = proc_grid[p + proc_row * npcols];
        if (rank != dest) {
          int tag = 2 * i + 1;
          int offset = lval_offset[Lcolp[i]];
          MPI_Isend(&Lvals[offset], lbuff_size, TACS_MPI_TYPE, dest, tag, comm,
                    &L_send_request[k]);
          k++;
        }
      }
    } else {
      // The receiving processes
      int source = proc_grid[source_proc_column + proc_row * npcols];
      int tag = 2 * i + 1;
      MPI_Irecv(Lbuff, lbuff_size, TACS_MPI_TYPE, source, tag, comm,
                &L_recv_request);
    }

    // There are four cases:
    // 1. The processor owns both the row and column elements required.
    // This is true of only the root processor.
    // 2. The processor owns the column but not the row and must wait for
    // any of the U_recv requests.
    // 3. The processor owns the row but not the column and must wait for
    // any of the L_recv requests.
    // 4. The processor owns neither the column or the row and must wait
    // for combinations of requests to begin.

    if (monitor_factor) {
      t_send_wait -= MPI_Wtime();
    }

    // Wait for the remaining sends to complete
    if (source_proc_row == proc_row) {
      MPI_Waitall(nprows - 1, U_send_request, U_send_status);
    }
    if (source_proc_column == proc_col) {
      MPI_Waitall(npcols - 1, L_send_request, L_send_status);
    }

    if (monitor_factor) {
      t_send_wait += MPI_Wtime();
      t_recv_wait -= MPI_Wtime();
    }

    // Wait for the receive to completex
    if (source_proc_column != proc_col) {
      MPI_Status L_recv_status;
      MPI_Wait(&L_recv_request, &L_recv_status);
    }
    if (source_proc_row != proc_row) {
      MPI_Status U_recv_status;
      MPI_Wait(&U_recv_request, &U_recv_status);
    }

    // Compute the bi-rank update to the remainder of the matrix
    // A[i+1:n,i+1:n] = A[i+1:n,i+1:n] - L[i:n,i]*U[i,i:n]

    if (monitor_factor) {
      t_recv_wait += MPI_Wtime();
      t_update -= MPI_Wtime();
    }

    // Initialize the L-pointer
    TacsScalar *L = Lbuff;
    if (source_proc_column == proc_col) {
      L = &Lvals[lval_offset[Lcolp[i]]];
    }

    for (int iip = Lcolp[i]; iip < Lcolp[i + 1]; iip++) {
      // Skip rows not locally owned
      int ii = Lrows[iip];
      int bii = bptr[ii + 1] - bptr[ii];
      if (get_proc_row(ii) != proc_row) {
        continue;
      }

      // Initialize the U-pointer
      TacsScalar *U = Ubuff;
      if (source_proc_row == proc_row) {
        U = &Uvals[uval_offset[Urowp[i]]];
      }

      for (int jjp = Urowp[i]; jjp < Urowp[i + 1]; jjp++) {
        // Skip columns not locally owned
        int jj = Ucols[jjp];
        int bjj = bptr[jj + 1] - bptr[jj];
        if (get_proc_column(jj) != proc_col) {
          continue;
        }

        TacsScalar *A = get_block(rank, ii, jj);

        if (A) {
          // Perform the GEMM update
          // A[i+1:n,i+1:n] = A[i+1:n,i+1:n] - L[i:n,i]*U[i,i:n]
          // A in bii x bjj
          // L in bii x bi
          // U in bi x bjj

          TacsScalar alpha = -1.0, beta = 1.0;
          BLASgemm("N", "N", &bii, &bjj, &bi, &alpha, L, &bii, U, &bi, &beta, A,
                   &bii);
          n_gemm++;
          TacsAddFlops(2 * bii * bjj * bi);
        }

        U += bi * bjj;
      }

      L += bi * bii;
    }

    if (monitor_factor) {
      t_update += MPI_Wtime();
    }
  }

  if (monitor_factor) {
    printf("[%d] Number of GEMM updates: %d\n", rank, n_gemm);
    printf("[%d] Update time:      %15.8f\n", rank, t_update);
    printf("[%d] Recv wait time:   %15.8f\n", rank, t_recv_wait);
    printf("[%d] Send wait time:   %15.8f\n", rank, t_send_wait);
  }

  // Remove memory for the block factorization
  delete[] temp_piv;
  delete[] temp_diag;
  delete[] temp_block;
  delete[] work;

  // Release memory for the data transfer
  delete[] Ubuff;
  delete[] Lbuff;

  delete[] U_send_request;
  delete[] L_send_request;
  delete[] U_send_status;
  delete[] L_send_status;
}
