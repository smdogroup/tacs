#ifndef TACS_PD_MAT_PIVOT_H
#define TACS_PD_MAT_PIVOT_H

#include "TACSObject.h"

/*!
  Parallel sparse block-cyclic matrix format with pivoting.

  Copyright (c) 2013 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.

  This code implements a parallel sparse block-cylic matrix
  format. The blocks are stored in column-major fortran order. The
  matrix multiplications, and factorizations are performed in
  parallel, and in place using LAPACK/BLAS for all block
  operations. Each block is full with no sub-block non-zero pattern
  taken into account.

  The parallelism is based on a 2D block cyclic format. This format is
  more difficult to implement but results in better parallelism for
  factorization. The 2D block-cyclic format repeats over blocks of
  rows and columns in the matrix. 

  The information required in matrix operations are:
  - the processor assigned to entry (i,j)
  - the processors in the block format column that own row i
  - the processors in the block format row that own column j
  - the processor block column/row
  
  The factorizations are performed with partial pivoting in the rows
  for stability. The main purpose of this matrix class is to be used
  for an interface problem (the Schur-complement problem) in a
  precondition or direct solve. Therefore, there is no need to save
  the original values in the matrix. As a result, the factorization is
  done in place, using the existing non-zero pattern.

  The main requirements for the code are:

  1. Initialize the non-zero pattern from the distributed contributions
  2. Transfer the values from the distributed contributions to the
  appropriate places.
  3. Compute the factorization fill-in.
  4. Compute the factorization in parallel.
  5. Perform back-solves in parallel.
*/
class PDMatPivot : public TACSObject {
 public:
  // Create a sparse matrix with the specified block pattern
  PDMatPivot( MPI_Comm _comm, int _num_block_rows, int _num_block_cols,
              const int *_block_ptr, const int * in_rowp, const int * in_cols,
              double fill );
  ~PDMatPivot();

  // Functions for various parts of the matrix
  // -----------------------------------------
  void getSize( int * nr, int * nc );
  void getProcessGridSize( int * _nprows, int * _npcols );
  void setMonitorFactorFlag( int flag );

  // Functions for setting values into the matrix
  // --------------------------------------------
  void zeroEntries();
  /*
  void addAllValues( int csr_bsize, int nvars, const int * vars,
                     const int * csr_rowp, const int * csr_cols, 
                     TacsScalar * vals );
  void addAlltoallValues( int csr_bsize, int nvars, const int * vars,
			  const int * csr_rowp, const int * csr_cols, 
			  TacsScalar * vals );
  */
  void setRand();

  // Matrix operations - note that factorization is in-place
  // -------------------------------------------------------
  void mult( TacsScalar * x, TacsScalar * y );
  void applyFactor( TacsScalar * x );
  // void applyFactor( int local_size, TacsScalar * x ); // This is faster
  void factor();

 private:
  void init_proc_grid();
  void init_mpi_type();

  // void init_nz_arrays();
  // void merge_nz_pattern( int root, int * rowp, int * cols,
  // int reorder_blocks );
  // void init_ptr_arrays( int * rowp, int * cols );

  void init_symbolic_factor( int _num_block_rows, 
                             int _num_block_cols,
                             const int *_block_ptr,
                             const int * in_rowp, 
                             const int * in_cols,
                             double fill,
                             int extra_nnz_per_row );

  int get_block_num( int var, const int * ptr );
  
  int add_values( int rank, int i, int j, 
		  int csr_bsize, int csr_i, int csr_j, 
		  TacsScalar * b );

  void compress_row( int row );
  void insert_row_nonzeros( int row, const int * new_cols,
			    int num_nnz );
  void swap_numerical_values( int row,
			      const int * pivot_sequence,
			      int num_pivots, 
			      TacsScalar * send_buff, 
			      TacsScalar * recv_buff );
  void add_non_zero_pattern( int row, 
			     const int * pivot, 
			     int num_pivots,
			     const int * col_rows, int num_col_rows,
			     int * marker, int * temp );

  void factor_panel( int col, const int * col_rows, int num_cols, 
                     TacsScalar * temp, const int * current_col_offset );

  // Given the i/j location within the matrix, determine the owner
  int get_block_owner( int i, int j ) const {
    i = i % nprows;
    j = j % npcols;
    return proc_grid[j + i*npcols];
  }

  // Get the process row, of the provided matrix row
  int get_proc_row( int row ) const {
    return row % nprows;
  }

  // Get the process column of the provided matrix column
  int get_proc_column( int col ) const {
    return col % npcols;
  }

  // Get the process column and row of the given rank process
  // Return 1 upon success 
  int get_proc_row_column( int rank, int * _proc_row, int * _proc_col ) const {
    for ( int i = 0; i < nprows; i++ ){
      for ( int j = 0; j < npcols; j++ ){
        if (proc_grid[j + i*npcols] == rank){
          *_proc_row = i;
          *_proc_col = j;
          return 1;
        }        
      }
    }

    *_proc_row = -1;
    *_proc_col = -1;
    return 0;
  }

  static const int BACKSOLVE_COLUMN_SIZE = 16;
  static const int BACKSOLVE_BUFF_SIZE = 64;

  // Retrieve the block matrix at the specified entry
  TacsScalar * get_block( int rank, int i, int j );

  // The communicator for this matrix
  MPI_Comm comm; 

  // This data controls how the data is assigned to the processors
  int npcols, nprows; // How many processors are assigned for each grid location
  int proc_col, proc_row; // The process column and row
  // npcols == number of columns in the process grid
  // nprows == number of rows in the process grid
  int *proc_grid;     // The processors assigned for each part of the grid

  // The block sizes for the matrix
  int *block_ptr; // len(bptr) = max(nrows, ncols)+1
  int max_block_size; // max(block_ptr[i+1] - block_ptr[i])
  int num_new_blocks; // Allocate new data size: num_new_blocks*max_block_size^2

  // The pivot sequence: block row/block index pairs
  int *pivots; // Note that this is a 2*block_ptr[num_block_rows]

  // Additional variables if reordering has been performed
  // int *orig_bptr;
  // int *perm, *iperm; // The permutation arrays

  // The storage information for the non-zero blocks
  int num_block_rows, num_block_cols;
  int max_block_array_size; // The maximum size of the arrays rows/cols

  // The compressed-row storage information
  int *rowp, *cols, *diag_offset, *row_len; // Non-zero data
  TacsScalar ** block_row_ptr; // Pointer to block data
  
  // The data storage infrastructure
  // This stores the data in a linked list. New entries are appended
  // to the list while the pointers to the blocks do not have to be
  // updated to match.
  class PDMatData {
  public:
    PDMatData( int _data_size ){ 
      data_size = _data_size;
      data = new TacsScalar[ data_size ];
      memset(data, 0, data_size*sizeof(TacsScalar));
      data_ptr = data;
      next = NULL;
    }
    ~PDMatData(){
      if (data){ delete [] data; }
    }
    TacsScalar * data; // The pointer to the beginning of the data
    TacsScalar * data_ptr; // The pointer to the next un-allocated mem location
    int data_size; // The total size of the allocated space
    PDMatData * next;
  } *data_root, *data_top;

  // Monitor the time spent in the factorization process
  int monitor_factor;

  // The data type that is used to communicate the location of the
  // maximum entry in each column
  static MPI_Datatype mat_block_type;
};

#endif
