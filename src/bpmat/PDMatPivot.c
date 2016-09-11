#include <stdlib.h>
#include "PDMatPivot.h"
#include "tacslapack.h"
#include "FElibrary.h"
#include "MatUtils.h"
#include "amd.h"

/*
  Copyright (c) 2013 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

/*
  This defines the structure type used in comparison operations for
  the maximum entry in the given column. This is used to reduce the
  data to the root processor or all processors during the
  factorization. 
  
  The data consists of the following:
  entry:       the value of the entry in the column
  block_index: the block index (row) that contains the entry
  index:       the index of the entry within the block
  rank:        the rank of the owner block
*/
typedef struct {
  TacsScalar entry;
  int block_index;
  int index;
  int rank;
} MatBlockIndex;

MPI_Datatype PDMatPivot::mat_block_type = NULL;

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

  Input:
  ------

  comm: the MPI communicator that defines the matrix

  csr_m, csr_n: The number of block-CSR rows and columns
  
  csr_bsize: The input block-CSR block size

  csr_vars: The global block-CSR variable numbers

  csr_nvars, csr_rowp, csr_cols: The CSR non-zero pattern
  
  csr_blocks_per_block: The number of CSR blocks per block of the
  non-zero pattern in the PDMatPivot class.

  Procedure:
  ----------
  
  1. Determine the block layout structure for PDMatPivot: populate bptr.
  2. Determine the PDMatPivot CSR structure from the block CSR structure
  provided on each processor.
  3. Pass the non-zero patterns to the root processor.
  4. Determine the fill-ins required for the factorization process on
  the root process.
  5. Pass the block structure back to all processes.
  6. Clean up time.
*/
/*
PDMatPivot::PDMatPivot( MPI_Comm _comm, int csr_m, int csr_n, 
                        int csr_bsize, const int * csr_vars, 
                        int csr_nvars, const int * csr_rowp, 
                        const int * csr_cols,
                        int csr_blocks_per_block, int reorder_blocks ){
  comm = _comm;
  monitor_factor = 0;
  perm = iperm = orig_block_ptr = NULL;

  int rank = 0, size = 0;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  
  // Determine the process grid
  init_proc_grid(size);

  // The total size of the matrix
  int m = csr_m*csr_bsize;
  int n = csr_n*csr_bsize;

  // Set up the block_ptr array
  int bsize = csr_blocks_per_block*csr_bsize;
  nrows = m/bsize; 
  ncols = n/bsize;
  if (m % bsize > 0){
    nrows++;
  }
  if (n % bsize > 0){
    ncols++;
  }

  int len_block_ptr = (nrows > ncols ? nrows : ncols) + 1;
  block_ptr = new int[len_block_ptr];
  max_bsize = bsize;

  if ( m >= n ){
    block_ptr[0] = 0;
    for ( int i = 0; i < nrows; i++ ){
      block_ptr[i+1] = block_ptr[i] + bsize;
    }
    block_ptr[nrows] = m;
    if (nrows > 0 && block_ptr[nrows] - block_ptr[nrows-1] > max_bsize){
      max_bsize = block_ptr[nrows] - block_ptr[nrows-1];
    }
  }
  else {
    block_ptr[0] = 0;
    for ( int i = 0; i < ncols; i++ ){
      block_ptr[i+1] = block_ptr[i] + bsize;
    }
    block_ptr[ncols] = n;
    if (ncols > 0 && block_ptr[ncols] - block_ptr[ncols-1] > max_bsize){
      max_bsize = block_ptr[ncols] - block_ptr[ncols-1];
    }
  }

  // Determine the block-CSR format for PDMatPivot.

  // The following approach allocates more memory than is actually
  // required. However, it is fixed to at most len(csr_rowp) + nvars+1
  // new allocations. This is probably not the end of the world.
  // A tighter bound could be found.
  int * rowp = new int[ nrows+1 ];
  memset(rowp, 0, (nrows+1)*sizeof(int));

  for ( int ip = 0; ip < csr_nvars; ip++ ){
    int i = csr_vars[ip];
    int ib = get_block_num(csr_bsize*i, block_ptr);

    for ( int jp = csr_rowp[ip]; jp < csr_rowp[ip+1]; jp++ ){
      rowp[ib+1]++;
    }
  }
  
  for ( int k = 0; k < nrows; k++ ){
    rowp[k+1] += rowp[k];
  }

  int col_size = rowp[nrows];
  int * cols = new int[col_size];
  
  for ( int ip = 0; ip < csr_nvars; ip++ ){
    int i = csr_vars[ip];
    int ib = get_block_num(csr_bsize*i, block_ptr);

    for ( int jp = csr_rowp[ip]; jp < csr_rowp[ip+1]; jp++ ){
      int j = csr_vars[csr_cols[jp]];
      int jb = get_block_num(csr_bsize*j, block_ptr);
      cols[rowp[ib]] = jb;
      rowp[ib]++;
    }
  }

  // Sort and uniqify the result
  for ( int k = nrows; k > 0; k-- ){
    rowp[k] = rowp[k-1];
  }
  rowp[0] = 0;

  int nodiag = 0;
  matutils::SortAndUniquifyCSR(nrows, rowp, cols, nodiag);
  
  int root = size-1;

  // This initializes Urowp/Ucols and Lcolp/Lrows on all procs
  merge_nz_pattern(root, rowp, cols, reorder_blocks);

  delete [] rowp;
  delete [] cols;

  // Initialize the component arrays
  init_nz_arrays();

  // Initialze data for the back-solves
  init_row_counts();
}
*/

/*
  Create a dense matrix.
*/
PDMatPivot::PDMatPivot( MPI_Comm _comm, 
                        int _num_block_rows, int _num_block_cols,
                        const int *_block_ptr, 
                        const int * in_rowp, const int * in_cols,
                        double fill ){
  comm = _comm;
  monitor_factor = 0; // Don't monitor the factorization by default
  
  // Set the number of non-zero blocks allocated each time
  // the allocated memory runs out
  num_new_blocks = 100; 

  // Set some of the data arrays to NULL
  block_ptr = NULL;
  rowp = cols = diag_offset = row_len = NULL;
  block_row_ptr = NULL;

  // Set the pivot rows/indices to NULL before factorization
  pivots = NULL;

  // Set the linked list of storage to NULL
  data_root = data_top = NULL;
  
  // Initialize the process grid
  init_proc_grid();
 
  // Perform the symbolic factorization
  // Extra non-zero entries for each row - note that this does not increase
  // the size of the memory in storage. Initially the memory allocated is
  // exactly that predicted by static pivoting.
  int extra_nnz_per_row = 0; 
  init_symbolic_factor(_num_block_rows, _num_block_cols,
                       _block_ptr, in_rowp, in_cols, fill,
                       extra_nnz_per_row);
  
  // Initialize the derived MPI data type 
  init_mpi_type();
}

/*
  Deallocate all the memory that this object allocated 
*/
PDMatPivot::~PDMatPivot(){
  // Delete the process grid information
  delete [] proc_grid;

  // Delete the pointer to the block starting locations
  delete [] block_ptr;

  // Delete the pivot sequence
  if (pivots){ delete [] pivots; }

  // Delete the row-oriented storage scheme
  if (rowp){ delete [] rowp; }
  if (cols){ delete [] cols; }
  if (row_len){ delete [] row_len; }
  if (diag_offset){ delete [] diag_offset; }
  if (block_row_ptr){ delete [] block_row_ptr; }
  
  // Delete the data allocated during the factorization
  while (data_root){
    PDMatData * temp = data_root;
    data_root = data_root->next;
    delete temp;
  }
}

/*
  Set up the MPI data type required to pass the pivot information
  between processors. This includes the maximum pivot entry, the block
  index, index within the block and the rank of the owner.  Setting up
  this type avoids having to make multiple MPI calls with different
  data types.  
*/
void PDMatPivot::init_mpi_type(){
  if (!mat_block_type){
    // Set up the MPI data type to pass information back and forth between
    MatBlockIndex max_index;
    int block[] = {1, 1, 1, 1};
    MPI_Aint disp[4];
    MPI_Datatype types[] = {TACS_MPI_TYPE, MPI_INT, MPI_INT, MPI_INT};
    
    MPI_Get_address(&max_index, disp);
    MPI_Get_address(&max_index.block_index, disp+1);
    MPI_Get_address(&max_index.index, disp+2);
    MPI_Get_address(&max_index.rank, disp+3);
    for ( int i = 3; i >= 0; i-- ){
      disp[i] -= disp[0];
    }
    
    MPI_Type_struct(4, block, disp, types, &mat_block_type);
    MPI_Type_commit(&mat_block_type);
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
/*
void PDMatPivot::merge_nz_pattern( int root, int * rowp, int * cols, 
                                   int reorder_blocks ){
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // Count up the maximum size of the CSR array
  int col_size = rowp[nrows];
  int root_col_size = 0;
  MPI_Reduce(&col_size, &root_col_size, 1, MPI_INT, MPI_SUM, root, comm);

  int * root_rowp = NULL;
  int * root_cols = NULL; 
  int * root_all_rowp = NULL;

  int * recv_row = NULL;
  int * recv_ptr = NULL;
  int * recv_count = NULL;

  if (rank == root){
    root_rowp = new int[nrows+1];
    root_rowp[0] = 0;

    root_cols = new int[root_col_size];
    root_all_rowp = new int[(nrows+1)*size];

    recv_row = new int[size*ncols];
    recv_ptr = new int[size+1];
    recv_count = new int[size];
  }

  MPI_Gather(rowp, nrows+1, MPI_INT, 
	     root_all_rowp, (nrows+1), MPI_INT, root, comm);

  for ( int i = 0; i < nrows; i++ ){
    if (rank == root){
      recv_ptr[0] = 0;
      for ( int k = 0; k < size; k++ ){	
	recv_count[k] = (root_all_rowp[i+1 + k*(nrows+1)] - 
			 root_all_rowp[i + k*(nrows+1)]);
	recv_ptr[k+1] = recv_ptr[k] + recv_count[k];
      }
    }
    
    int row_size = rowp[i+1] - rowp[i];
    MPI_Gatherv(&cols[rowp[i]], row_size, MPI_INT,
		recv_row, recv_count, recv_ptr, MPI_INT, root, comm);

    if (rank == root){
      // Sort and uniquify this row
      int row_len = recv_ptr[size];
      row_len = FElibrary::uniqueSort(recv_row, row_len);

      // Now add this row to the CSR data structure
      root_rowp[i+1] = root_rowp[i] + row_len;
      for ( int j = 0, k = root_rowp[i]; j < row_len; j++, k++ ){
	root_cols[k] = recv_row[j];
      }
    }
  }

  if (reorder_blocks && nrows == ncols){
    // Allocate the permutation array
    perm = new int[ nrows ];
    iperm = new int[ ncols ];
    orig_block_ptr = new int[ nrows+1 ];
  }
  else{
    reorder_blocks = 0;
    if (rank == 0){
      printf("[%d] PDMatPivot: Cannot re-order blocks with nrows != ncols\n",
             rank);
    }
  }

  if (rank == root){
    delete [] root_all_rowp;
    delete [] recv_row;
    delete [] recv_ptr;
    delete [] recv_count;

    int m = block_ptr[nrows], n = block_ptr[ncols];
    if (reorder_blocks && n == m){
      if (m && n){
        // Use AMD to compute the reordering of the variables.
        double control[AMD_CONTROL], info[AMD_INFO];
        amd_defaults(control); // Use the default values
        amd_order(nrows, root_rowp, root_cols, perm, 
                  control, info);
        
        // Reorder the non-zero pattern to correspond to the new ordering
        // perm:  new variable i -> old variable perm[i]
        // iperm: old variable i -> new variable iperm[i]
        for ( int i = 0; i < nrows; i++ ){
          iperm[perm[i]] = i;
        }

        int * temp_rowp = new int[nrows+1];
        int * temp_cols = new int[root_col_size];
        
        temp_rowp[0] = 0;
        for ( int i = 0, p = 0; i < nrows; i++ ){
          int row = perm[i];
          int size = root_rowp[row+1] - root_rowp[row];
          for ( int jp = root_rowp[row]; jp < root_rowp[row+1]; jp++, p++ ){
            temp_cols[p] = iperm[root_cols[jp]];
          }
          // Sort the row
          temp_rowp[i+1] = p;
          size = FElibrary::uniqueSort(&temp_cols[temp_rowp[i]], size);
          if (size != temp_rowp[i+1] - temp_rowp[i]){
            printf("[%d] PDMatPivot: problem with the permutation array\n",
                   rank);
          }
        }
        
        delete [] root_rowp;
        delete [] root_cols;
        
        root_rowp = temp_rowp;
        root_cols = temp_cols;
      }
    }

    // perform the symbolic factorization
    int init_nnz = root_rowp[nrows];
    compute_symbolic_factor(&root_rowp, &root_cols, root_col_size);
    int final_nnz = root_rowp[nrows];    
    
    printf("[%d] PDMatPivot factorization: M: %d N: %d\n",
           root, m, n);
    if (m && n){
      printf("[%d] initial density: %4.3f factor fill in: %4.3f\n",
             root, (1.0*init_nnz)/(nrows*ncols), 
             (1.0*(final_nnz - init_nnz))/init_nnz);
    }

    // create the Urowp/Ucols and Lcolp/Lrows arrays
    init_ptr_arrays(root_rowp, root_cols);

    delete [] root_rowp;
    delete [] root_cols;
  }

  if (reorder_blocks){
    MPI_Bcast(perm, nrows, MPI_INT, root, comm);

    if (rank != root){
      for ( int i = 0; i < nrows; i++ ){
        iperm[perm[i]] = i;
      }
    }

    // Set the block_ptr array with the new permutation
    for ( int i = 0; i < nrows; i++ ){
      int ip = perm[i];
      orig_block_ptr[i+1] = block_ptr[ip+1] - block_ptr[ip];
    }
    orig_block_ptr[0] = 0;
    for ( int i = 0; i < nrows; i++ ){
      orig_block_ptr[i+1] += orig_block_ptr[i];
    }

    // Flip the two
    int * temp = block_ptr;
    block_ptr = orig_block_ptr;
    orig_block_ptr = temp;    
  }
  
  // Broadcast the nonzero pattern
  if (rank != root){
    Urowp = new int[nrows+1];
    Lcolp = new int[ncols+1];
  }

  MPI_Bcast(Urowp, nrows+1, MPI_INT, root, comm);
  MPI_Bcast(Lcolp, ncols+1, MPI_INT, root, comm);

  // Allocate space for the remaining arrays
  if (rank != root){
    Ucols = new int[Urowp[nrows]];
    Lrows = new int[Lcolp[ncols]];
  }

  MPI_Bcast(Ucols, Urowp[nrows], MPI_INT, root, comm);
  MPI_Bcast(Lrows, Lcolp[ncols], MPI_INT, root, comm);
}
*/

/*
  Perform the symbolic factorization of the non-zero pattern.

  This takes in the initial compressed-sparse row format of the
  matrix.  It then computes the fill-in due to the factorization and
  the static factored fill pattern. This function then allocates the
  memory that is required for the storage of the non-zero
  components. This allocates exactly the required amount of space, so
  any additional fill during the numerical factorization will have to
  be allocated in a new buffer. Finally, the initial non-zero pattern
  of the matrix is restored in the row_len/cols array.

  input:
  fill:           the expected fill-in during the factorization
  num_block_rows: the number of block rows
  num_block_cols: the number of block columns
  in_rowp:        the row pointer for the initial non-zero pattern
  in_cols:        the column indices for the initial non-zero pattern  
*/
void PDMatPivot::init_symbolic_factor( int _num_block_rows, 
                                       int _num_block_cols,
                                       const int *_block_ptr,
                                       const int * in_rowp, 
                                       const int * in_cols,
                                       double fill,
                                       int extra_nnz_per_row ){
  // Copy over the number of rows/columns
  num_block_rows = _num_block_rows;
  num_block_cols = _num_block_cols;

  // Initialize the block pointer
  block_ptr = new int[ num_block_rows+1 ];
  memcpy(block_ptr, _block_ptr, (num_block_rows+1)*sizeof(int));

  // Compute the maximum block pointer size
  max_block_size = 0;
  for ( int i = 0; i < num_block_rows; i++ ){
    if (block_ptr[i+1] - block_ptr[i] > max_block_size){
      max_block_size = block_ptr[i+1] - block_ptr[i];
    }
  }

  if (proc_row < 0 || proc_col < 0){
    rowp = row_len = diag_offset = cols = NULL;
    block_row_ptr = NULL;
    pivots = NULL;
    data_root = data_top = NULL;
    return;
  }

  // Create a temporary array for the column markers
  int * marker = new int[ num_block_cols ];
  memset(marker, 0, num_block_cols*sizeof(int));

  // Create a temporary array of the column indices
  int * tcols = new int[ num_block_cols ];

  // Compute the approximate amount of fill-in expected during
  // the factorization
  max_block_array_size = 
    (fill*in_rowp[num_block_rows] + 
     extra_nnz_per_row);

  // Allocate the space for the row-oriented storage
  rowp = new int[ num_block_rows+1 ];
  row_len = new int[ num_block_rows ];
  diag_offset = new int[ num_block_rows ];
  cols = new int[ max_block_array_size ];

  // Initialize the row pointer array
  rowp[0] = 0;

  for ( int i = 0; i < num_block_rows; i++ ){
    // Check if the columns
    int num_cols = 0;
    for ( int j = in_rowp[i]; j < in_rowp[i+1]; j++ ){
      int col = in_cols[j]; 
      marker[col] = i+1;
      tcols[num_cols] = col;
      num_cols++;
    }
    
    // Merge the diagonal into the array
    if (marker[i] != i+1){
      marker[i] = i+1;
      tcols[num_cols] = i;
      num_cols++;
    }

    // Perform the symbolic factorization - generating new entries
    for ( int j = 0; j < num_cols; j++ ){
      // Search through the current column if it is less than
      // the current row i
      int col = tcols[j];
      if (col < i){
	// The end of the col row
	int kp_end = rowp[col] + row_len[col];

	// Start with the first entry after the diagonal in row, cols[j] 
	// k is the index into cols for row cols[j]
	int kp = rowp[col] + diag_offset[col] + 1;
	for ( ; kp < kp_end; kp++ ){
	  int k = cols[kp];
	  // Add the marker 
	  if (marker[k] != i+1){
	    marker[k] = i+1;
	    tcols[num_cols] = k;
	    num_cols++;
	  }
	}
      }
    }
      
    // Sort the row
    qsort(tcols, num_cols, sizeof(int), FElibrary::comparator);

    // Look for the diagonal entry
    int * item = (int*)bsearch(&i, tcols, num_cols, sizeof(int), 
			       FElibrary::comparator);
    diag_offset[i] = item - tcols;

    // Check if the size will be exceeded by adding the new elements
    if (num_cols + extra_nnz_per_row + rowp[i] > max_block_array_size){
      // Estimate the number of non-zeros that will be requried for
      // the remainder of the matrix - and multiply by 2
      max_block_array_size = max_block_array_size + 
        2*(num_block_rows - i)*(rowp[i]/i);

      int * temp = cols;
      cols = new int[ max_block_array_size ];
      memcpy(cols, temp, rowp[i]*sizeof(int));
      delete [] temp;
    }

    row_len[i] = num_cols;
    rowp[i+1] = rowp[i] + num_cols + extra_nnz_per_row;
    memcpy(&cols[rowp[i]], tcols, num_cols*sizeof(int));
    for ( int ip = rowp[i] + num_cols; ip < rowp[i+1]; ip++ ){
      cols[ip] = -1;
    }
  }

  // Free the memory 
  delete [] marker;
  delete [] tcols;

  // Restore the initial non-zero pattern into the array
  // But leave the space for the additional non-zeros that will be
  // allocated later. As a result, the array rowp is not touched.
  for ( int i = 0; i < num_block_rows; i++ ){
    int num_cols = 0;
    diag_offset[i] = -1;
    for ( int ip = rowp[i], jp = in_rowp[i]; 
	  jp < in_rowp[i+1]; ip++, jp++, num_cols++ ){
      cols[ip] = in_cols[jp];
      if (cols[ip] == i){
	diag_offset[i] = ip - rowp[i];
      }
    }
    row_len[i] = num_cols;
  }

  // Allocate the pointers to the blocks of data storage
  block_row_ptr = new TacsScalar*[ max_block_array_size ];
  memset(block_row_ptr, 0, max_block_array_size*sizeof(TacsScalar*));

  // Count up the size of the space required to store all the non-zeros 
  // in the initial non-zero pattern of the matrix
  int data_size = 0;
  for ( int i = proc_row; i < num_block_rows; i += nprows ){
    int bi = block_ptr[i+1] - block_ptr[i];

    int jp_end = rowp[i] + row_len[i];
    for ( int jp = rowp[i]; jp < jp_end; jp++ ){
      int j = cols[jp];
      if (proc_col == get_proc_column(j)){
        int bj = block_ptr[j+1] - block_ptr[j];
        data_size += bi*bj;
      }
    }
  }

  // Make sure that the data is all destroyed
  while (data_root){
    PDMatData * temp = data_root;
    data_root = data_root->next;
    delete temp;
  }

  // Allocate the root of the linked list
  data_root = new PDMatData(data_size);

  // Set the pointer into the block of the data storage
  TacsScalar * data_ptr = data_root->data;

  for ( int i = proc_row; i < num_block_rows; i += nprows ){
    int bi = block_ptr[i+1] - block_ptr[i];

    int jp_end = rowp[i] + row_len[i];
    for ( int jp = rowp[i]; jp < jp_end; jp++ ){
      int j = cols[jp];
      if (proc_col == get_proc_column(j)){
        int bj = block_ptr[j+1] - block_ptr[j];

        // Set the pointer into the data storage
        block_row_ptr[jp] = data_ptr;
        data_ptr += bi*bj;
      }
    }
  }

  data_root->data_ptr = data_ptr;
  data_top = data_root;
}

/*
  Set up the process grid. 

  This function assigns a unique processor to each process
  block. Therefore, in general there may be processes left without
  blocks. For np = 1, 2, 3, 4, 6, 8, 9, 12, 15, 16,...  it should work
  fine. More flexibility with the block assignment should be
  considered for future versions of the code.  

  input:
  size: the MPI size of the comm
*/
void PDMatPivot::init_proc_grid(){
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  // Check for the nearest perfect square  
  int n = 1;
  while (n*n > 0){
    if (n*n <= size && size < (n+1)*(n+1)){
      break;
    }
    n++;
  }

  // Set the size of the process grid based on the
  // search for the perfect square
  nprows = n;
  npcols = size/n;

  proc_grid = new int[ nprows*npcols ];  
  for ( int i = 0, p = 0; i < nprows; i++ ){
    for ( int j = 0; j < npcols; j++, p++ ){
      proc_grid[j + i*npcols] = p;
    }
  }

  // Set the process row/column for this processor. Note
  // that negative values indicate that this is not in
  // the process grid and should not participate in the
  // computation
  proc_row = -1;
  proc_col = -1;
  get_proc_row_column(rank, &proc_row, &proc_col);
}

/*
  Get the size of the matrix
*/
void PDMatPivot::getSize( int * nr, int * nc ){
  *nr = block_ptr[num_block_rows];
  *nc = block_ptr[num_block_cols];
}

/*
  Retrieve the size of the process grid.
*/
void PDMatPivot::getProcessGridSize( int * _nprows, int * _npcols ){
  *_nprows = nprows;
  *_npcols = npcols;
}

/*
  Set the flag that prints out the factorization time
*/
void PDMatPivot::setMonitorFactorFlag( int flag ){
  monitor_factor = flag;
}

/*
  Zero all the matrix entries.
*/
void PDMatPivot::zeroEntries(){
  // Zero all the un-allocated data in each of the 
  // allocated blocks
  PDMatData * current = data_root;
  
  while (current){
    int size = current->data_ptr - current->data;
    memset(current->data, 0, size*sizeof(TacsScalar));
    current = current->next;
  }
}

/*
  Add values into the matrix. This function is collective on all
  processes in the PDMatPivot comm.

  First, add the on-process components. Next, add the off-process
  components from each process. This approach tries to minimize the
  amount of extra memory required for the off-process contributions.
  This code uses MPI_Gatherv for each process rather than a single
  call to MPI_Alltoallv which may be faster, but requires more memory.
*/
/*
void PDMatPivot::addAllValues( int csr_bsize, int nvars, const int * vars,
                          const int * csr_rowp, const int * csr_cols, 
                          TacsScalar * vals ){
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  
  int b2 = csr_bsize*csr_bsize;

  // Count up the number of recvs for each processor 
  int * recv_count = new int[ size ];
  int * recv_ptr = new int[ size+1 ];

  int csr_size = csr_rowp[nvars];
  int * jblock = new int[csr_size];

  int * send_counts = new int[ size ];
  memset(send_counts, 0, size*sizeof(int));
  
  // Go through and add any contributions from the local process
  for ( int ip = 0; ip < nvars; ip++ ){
    int i = csr_bsize*vars[ip];
    int ioff = 0, ib = 0;
    if (orig_block_ptr){
      ib = get_block_num(i, orig_block_ptr);
      ioff = i - orig_block_ptr[ib];
      ib = iperm[ib];
    }
    else {
      ib = get_block_num(i, block_ptr);
      ioff = i - block_ptr[ib];
    }
        
    // Determine the local block for vars
    for ( int jp = csr_rowp[ip]; jp < csr_rowp[ip+1]; jp++ ){
      int j = csr_bsize*vars[csr_cols[jp]];
      int joff = 0;
          
      // Get the block numbers
      if (orig_block_ptr){
        int jb = get_block_num(j, orig_block_ptr);
        joff = j - orig_block_ptr[jb];
        jblock[jp] = iperm[jb];
      }
      else {
        jblock[jp] = get_block_num(j, block_ptr);
        joff = j - block_ptr[jblock[jp]];
      }
      
      int owner = get_block_owner(ib, jblock[jp]);
      if (owner == rank){
        add_values(rank, ib, jblock[jp],
		   csr_bsize, ioff, joff, &vals[b2*jp]);       
      }
      else {
        send_counts[owner]++;
      }
    }
  }

  int max_send_size = send_counts[0];
  for ( int k = 1; k < size; k++ ){
    if (send_counts[k] > max_send_size){
      max_send_size = send_counts[k];
    }
  }

  delete [] send_counts;

  int * send_ivals = new int[max_send_size];
  int * send_jvals = new int[max_send_size];
  TacsScalar * send_vals = new TacsScalar[b2*max_send_size];

  // Iterate over each process
  for ( int k = 0; k < size; k++ ){
    int send_count = 0;

    if (rank != k){
      // Figure out what to add where
      for ( int ip = 0; ip < nvars; ip++ ){
        int i = csr_bsize*vars[ip];
        int ib;
        if (orig_block_ptr){
          ib = get_block_num(i, orig_block_ptr);
          ib = iperm[ib];
        }
        else {
          ib = get_block_num(i, block_ptr);
        }
        
        // Determine the local block for vars
        for ( int jp = csr_rowp[ip]; jp < csr_rowp[ip+1]; jp++ ){
          int j = csr_bsize*vars[csr_cols[jp]];
         
          if (get_block_owner(ib, jblock[jp]) == k){
            send_ivals[send_count] = i;
            send_jvals[send_count] = j;
            memcpy(&send_vals[b2*send_count], 
                   &vals[b2*jp], b2*sizeof(TacsScalar));
            send_count++;
          }
        }
      }
    }

    MPI_Gather(&send_count, 1, MPI_INT,
               recv_count, 1, MPI_INT, k, comm);
    
    // Count up the reciving information
    int nrecv = 0;
    int * recv_ivals = NULL;
    int * recv_jvals = NULL;
    TacsScalar * recv_vals = NULL;

    if (rank == k){
      recv_ptr[0] = 0;
      for ( int n = 0; n < size; n++ ){
        recv_ptr[n+1] = recv_ptr[n] + recv_count[n];
      }

      nrecv = recv_ptr[size];
    
      recv_ivals = new int[nrecv];
      recv_jvals = new int[nrecv];
      recv_vals = new TacsScalar[b2*nrecv];
    }

    MPI_Gatherv(send_ivals, send_count, MPI_INT,
                recv_ivals, recv_count, recv_ptr, MPI_INT, k, comm);
    MPI_Gatherv(send_jvals, send_count, MPI_INT,
                recv_jvals, recv_count, recv_ptr, MPI_INT, k, comm);
    
    // Reset the array sizes to send the blocks
    send_count *= b2;
    for ( int n = 0; n < size; n++ ){
      recv_count[n] *= b2;
      recv_ptr[n] *= b2;
    }

    // Send the values
    MPI_Gatherv(send_vals, send_count, TACS_MPI_TYPE,
                recv_vals, recv_count, recv_ptr, TACS_MPI_TYPE, k, comm);

    if (rank == k){
      // Allocate the sending/receiving arrays
      for ( int n = 0; n < nrecv; n++ ){
        int i = recv_ivals[n];
        int j = recv_jvals[n];
        
        int ib, jb; // The blocks
        int ioff, joff;
        if (orig_block_ptr){
          ib = get_block_num(i, orig_block_ptr);
          jb = get_block_num(j, orig_block_ptr);
          ioff = i - orig_block_ptr[ib];
          joff = j - orig_block_ptr[jb];
          ib = iperm[ib];
          jb = iperm[jb];
        }
        else {
          ib = get_block_num(i, block_ptr);
          jb = get_block_num(j, block_ptr);
          ioff = i - block_ptr[ib];
          joff = j - block_ptr[jb];
        }
        
        // Add this to the correct block ...
        add_values(rank, ib, jb,
                   csr_bsize, ioff, joff, &recv_vals[b2*n]);
      }

      delete [] recv_ivals;
      delete [] recv_jvals;
      delete [] recv_vals;
    }
  }

  delete [] jblock;

  delete [] recv_count;
  delete [] recv_ptr;

  delete [] send_ivals;
  delete [] send_jvals;
  delete [] send_vals;
}
*/

/*
  Add values into the matrix. This function is collective on all PDMatPivot
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
/*  
void PDMatPivot::addAlltoallValues( int csr_bsize, int nvars, const int * vars,
			       const int * csr_rowp, const int * csr_cols, 
			       TacsScalar * vals ){
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  
  int b2 = csr_bsize*csr_bsize;

  // Count up the number of recvs for each processor 
  int * recv_counts = new int[ size ];
  int * recv_ptr = new int[ size+1 ];

  int * send_counts = new int[ size ];
  int * send_ptr = new int[ size+1 ];
  memset(send_counts, 0, size*sizeof(int));
 
  // Go through and add any contributions from the local process
  for ( int ip = 0; ip < nvars; ip++ ){
    int i = csr_bsize*vars[ip];
    int ioff = 0, ib = 0;
    if (orig_block_ptr){
      ib = get_block_num(i, orig_block_ptr);
      ioff = i - orig_block_ptr[ib];
      ib = iperm[ib];
    }
    else {
      ib = get_block_num(i, block_ptr);
      ioff = i - block_ptr[ib];
    }
        
    // Determine the local block for vars
    for ( int jp = csr_rowp[ip]; jp < csr_rowp[ip+1]; jp++ ){
      int j = csr_bsize*vars[csr_cols[jp]];
      int joff = 0, jb = 0;
          
      // Get the block numbers
      if (orig_block_ptr){
        jb = get_block_num(j, orig_block_ptr);
        joff = j - orig_block_ptr[jb];
        jb = iperm[jb];
      }
      else {
        jb = get_block_num(j, block_ptr);
        joff = j - block_ptr[jb];
      }
      
      int owner = get_block_owner(ib, jb);
      if (owner == rank){
        add_values(rank, ib, jb,
		   csr_bsize, ioff, joff, &vals[b2*jp]);
      }
      else {
        send_counts[owner]++;
      }
    }
  }

  // Send the counts to all processors
  MPI_Alltoall(send_counts, 1, MPI_INT,
	       recv_counts, 1, MPI_INT, comm);
  
  // Sum up the send and recv pointers
  send_ptr[0] = 0;
  recv_ptr[0] = 0;
  for ( int k = 0; k < size; k++ ){
    send_ptr[k+1] = send_ptr[k] + send_counts[k];
    recv_ptr[k+1] = recv_ptr[k] + recv_counts[k];
  }

  // Create a buffer that is large enough to send everything at once
  // Zero out the send counter
  memset(send_counts, 0, size*sizeof(int));
  int * send_index = new int[2*send_ptr[size]];
  TacsScalar * send_vals = new TacsScalar[b2*send_ptr[size]];

  // Go back through and copy over values to pass to all the processes
  for ( int ip = 0; ip < nvars; ip++ ){
    int i = csr_bsize*vars[ip];
    int ib = 0;
    if (orig_block_ptr){
      ib = get_block_num(i, orig_block_ptr);
      ib = iperm[ib];
    }
    else {
      ib = get_block_num(i, block_ptr);
    }
        
    // Determine the local block for vars
    for ( int jp = csr_rowp[ip]; jp < csr_rowp[ip+1]; jp++ ){
      int j = csr_bsize*vars[csr_cols[jp]];
      int jb = 0;
          
      // Get the block numbers
      if (orig_block_ptr){
        jb = get_block_num(j, orig_block_ptr);
        jb = iperm[jb];
      }
      else {
        jb = get_block_num(j, block_ptr);
      }
      
      int owner = get_block_owner(ib, jb);
      if (owner != rank){
	int sc = send_ptr[owner] + send_counts[owner];
	send_counts[owner]++;

	// Set the index values
	send_index[2*sc]   = i;
	send_index[2*sc+1] = j;

	// Copy the values to the send buffer
	memcpy(&send_vals[b2*sc],
	       &vals[b2*jp], b2*sizeof(TacsScalar));
      }
    }
  }

  // Send the indices and the values
  int * send_counts2 = new int[ size ];
  int * send_ptr2 = new int[ size ];
  int * recv_counts2 = new int[ size ];
  int * recv_ptr2 = new int[ size ];

  // Adjust the pointers/counts to receive the (i,j) indices 
  for ( int k = 0; k < size; k++ ){
    send_counts2[k] = 2*send_counts[k];
    send_ptr2[k] = 2*send_ptr[k];
    recv_counts2[k] = 2*recv_counts[k];
    recv_ptr2[k] = 2*recv_ptr[k];
  }

  int * recv_index = new int[ 2*recv_ptr[size] ];
  TacsScalar * recv_vals = new TacsScalar[ b2*recv_ptr[size] ];

  MPI_Alltoallv(send_index, send_counts2, send_ptr2, MPI_INT,
		recv_index, recv_counts2, recv_ptr2, MPI_INT, comm);

  // Adjust the pointers to receive the values
  for ( int k = 0; k < size; k++ ){
    send_counts2[k] = b2*send_counts[k];
    send_ptr2[k] = b2*send_ptr[k];
    recv_counts2[k] = b2*recv_counts[k];
    recv_ptr2[k] = b2*recv_ptr[k];
  }

  MPI_Alltoallv(send_vals, send_counts2, send_ptr2, TACS_MPI_TYPE,
		recv_vals, recv_counts2, recv_ptr2, TACS_MPI_TYPE, comm);

  // Free all the pointers that are not required anymore
  delete [] send_counts2;
  delete [] send_ptr2;
  delete [] recv_counts2;
  delete [] recv_ptr2;
  delete [] send_counts;
  delete [] send_ptr;

  // Free the send buffer indices and values
  delete [] send_index;
  delete [] send_vals;

  int nrecv = recv_ptr[size];

  // Add the values from the receiving arrays
  for ( int n = 0; n < nrecv; n++ ){
    int i = recv_index[2*n];
    int j = recv_index[2*n+1];
        
    int ib, jb; // The block indices
    int ioff, joff;
    if (orig_block_ptr){
      ib = get_block_num(i, orig_block_ptr);
      jb = get_block_num(j, orig_block_ptr);
      ioff = i - orig_block_ptr[ib];
      joff = j - orig_block_ptr[jb];
      ib = iperm[ib];
      jb = iperm[jb];
    }
    else {
      ib = get_block_num(i, block_ptr);
      jb = get_block_num(j, block_ptr);
      ioff = i - block_ptr[ib];
      joff = j - block_ptr[jb];
    }
    
    // Add this to the correct block ...
    add_values(rank, ib, jb,
	       csr_bsize, ioff, joff, &recv_vals[b2*n]);
  }
  
  delete [] recv_index;
  delete [] recv_vals;
  delete [] recv_counts;
  delete [] recv_ptr;
}
*/
/*
  Determine the block number i such that var is within the interval:

  ptr[i] <= var < ptr[i+1]

  This uses a bisection search method since the block_ptr array must
  be sorted in ascending order (otherwise there would be negative
  block sizes!). 

  input:
  var:   the variable to search for
  ptr:   the block->variable range array

  returns:
  the interval satisfying the criterion above
*/
int PDMatPivot::get_block_num( int var, const int * ptr ){
  int high = num_block_cols;
  if (num_block_rows > high){
    high = num_block_rows; 
  }

  int low = 0;

  if ((var < ptr[low]) || (var >= ptr[high])){
    int rank;
    MPI_Comm_rank(comm, &rank);
    fprintf(stderr, "[%d] PDMatPivot::get_block_num(%d) out of range\n",
	    rank, var);
    return -1;
  }

  if (ptr[low] == var){
    return low;
  }
  
  while (high - low > 1){
    int mid = low + (int)((high - low)/2);

    if (ptr[mid] == var){
      return mid; // quick return
    }
    else if (var > ptr[mid]){
      low = mid;
    }
    else {
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
int PDMatPivot::add_values( int rank, int i, int j, 
                            int csr_bsize, int ioff, int joff, 
                            TacsScalar * a ){
  TacsScalar * A = get_block(rank, i, j);

  if (A){
    int bi = block_ptr[i+1] - block_ptr[i];
    int bj = block_ptr[j+1] - block_ptr[j];

    if ((ioff >= 0 && ioff + csr_bsize <= bi) &&
	(joff >= 0 && joff + csr_bsize <= bj)){
      for ( int m = 0; m < csr_bsize; m++ ){
	for ( int n = 0; n < csr_bsize; n++ ){
	  A[(ioff + m) + bi*(joff + n)] += a[csr_bsize*m + n];
	}
      }
      
      return 1;
    }
  }
  else {
    fprintf(stderr, "[%d] PDMatPivot: Error, (%d, %d) not in nz-pattern\n",
            rank, i, j);
  }

  return 0;
}
  
/*
  Assign randomly generated entries to the matrix.

  Note that this uses the rand() function from stdlib. The matrix
  entries lie within the interval [-0.5, 0.5].
*/
void PDMatPivot::setRand(){
  // Set random values everywhere
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // Keep track of the number of entries set 
  int nvars = 0;
  if (rank > 0){
    MPI_Recv(&nvars, 1, MPI_INT, rank-1, 0, 
	     comm, MPI_STATUS_IGNORE);
  }

  // Call rand size times
  for ( int i = 0; i < nvars; i++ ){ rand(); }

  PDMatData * current = data_root;
  while (current){
    int local_nvars = current->data_ptr - current->data;
    for ( int i = 0; i < local_nvars; i++ ){
      current->data[i] = -0.5 + (1.0*rand())/RAND_MAX;
    }
    current = current->next;
    nvars += local_nvars;
  }

  if (rank < size-1){
    MPI_Send(&nvars, 1, MPI_INT, rank+1, 0, comm);
  }
}

/*
  Multiply y = A*x

  Note that this code assumes that x/y are on all processors.  This
  function will be less than useless after the matrix is factored
  since the factorization over-writes the original matrix entries.  
*/
void PDMatPivot::mult( TacsScalar * x, TacsScalar * y ){
  int size = block_ptr[num_block_rows];
  TacsScalar * ty = new TacsScalar[size];
  memset(ty, 0, size*sizeof(TacsScalar));

  if (proc_row >= 0 && proc_col >= 0){
    // Compute the on-process parts for the diagonal contributions
    for ( int i = proc_row; i < num_block_rows; i += nprows ){
      int bi = block_ptr[i+1] - block_ptr[i];
      
      int jp_end = rowp[i] + row_len[i];
      for ( int jp = rowp[i]; jp < jp_end; jp++ ){
	int j = cols[jp];
	if (proc_col == get_proc_column(j)){
	  int bj = block_ptr[j+1] - block_ptr[j];
	  TacsScalar * A = block_row_ptr[jp];
	  
	  TacsScalar alpha = 1.0, beta = 1.0;
	  int one = 1;
	  BLASgemv("N", &bi, &bj, &alpha, A, &bi,
		   &x[block_ptr[j]], &one, 
		   &beta, &ty[block_ptr[i]], &one);
	}
      }
    }
  }

  MPI_Allreduce(ty, y, block_ptr[num_block_rows], 
		TACS_MPI_TYPE, MPI_SUM, comm);

  delete [] ty;
}

/*
  Apply the LU factorization to x

  x = U^{-1} L^{-1} x

  The factorization is performed in place. No check is performed to
  ensure that you've factored the matrix first, so beware!

  The pivots must be applied during the application of the lower
  factor, otherwise the results will be wrong. This requires a
  fine-grain parallelism that doesn't achieve ideal performance.
  Small gains could be achieved by performing the permutations in a
  batch. The degree to which this would improve the computational time
  is not clear, especially in the context of the entire factorization
  process.
*/
void PDMatPivot::applyFactor( TacsScalar * x ){
  int rank;
  MPI_Comm_rank(comm, &rank);

  TacsScalar * x_local = new TacsScalar[max_block_size];

  // Apply the lower part of the factorization
  for ( int i = 0; i < num_block_rows; i++ ){
    int ni = block_ptr[i];
    int bi = block_ptr[i+1] - ni;
    
    // Get the owner of the diagonal block
    int diag_owner = get_block_owner(i, i);

    // Perform the permutation of the entries of x
    for ( int ii = block_ptr[i]; ii < block_ptr[i+1]; ii++ ){
      // Swap entry ii with block 
      int swap_row = pivots[2*ii];
      int swap_index = pivots[2*ii+1];
      int jj = block_ptr[swap_row] + swap_index;
    
      // Get the owner of the swap block
      int swap_owner = get_block_owner(swap_row, i);

      // Swap the entries
      if (diag_owner == swap_owner){
	TacsScalar t = x[ii];
	x[ii] = x[jj];
	x[jj] = t;
      }
      else {
	if (rank == diag_owner){
	  TacsScalar t;
	  MPI_Request send_req;
	  MPI_Isend(&x[ii], 1, TACS_MPI_TYPE, swap_owner, i, comm, &send_req);
	  MPI_Recv(&t, 1, TACS_MPI_TYPE, swap_owner, i, comm, MPI_STATUS_IGNORE);
	  MPI_Wait(&send_req, MPI_STATUS_IGNORE);
	  x[ii] = t;
	}
	else if (rank == swap_owner){
	  TacsScalar t;
	  MPI_Request send_req;
	  MPI_Isend(&x[jj], 1, TACS_MPI_TYPE, diag_owner, i, comm, &send_req);
	  MPI_Recv(&t, 1, TACS_MPI_TYPE, swap_owner, i, comm, MPI_STATUS_IGNORE);
	  MPI_Wait(&send_req, MPI_STATUS_IGNORE);
	  x[jj] = t;
	}
      }
    }

    // Solve L*x = b
    if (rank == diag_owner){
      TacsScalar * D = block_row_ptr[rowp[i] + diag_offset[i]];
      int incx = 1;
      BLAStrsv("L", "N", "U", &bi, D, &bi, &x[ni], &incx);
    
      // Send the result to the other processors in the
      // column group
      for ( int k = 0; k < nprows; k++ ){
	int dest = get_block_owner(k, i);
	if (dest != rank){
	  MPI_Send(&x[ni], bi, TACS_MPI_TYPE, dest, i, comm);
	}
      }
    }

    // Cycle over the columns
    if (proc_col == get_proc_column(i)){
      if (rank != diag_owner){
	MPI_Recv(&x[ni], bi, TACS_MPI_TYPE, diag_owner, i, comm, MPI_STATUS_IGNORE);
      }

      int j = i+1;
      for ( ; j < num_block_rows; j++ ){
        int nj = block_ptr[j];
        int bj = block_ptr[j+1] - nj;

        // Find the column such that cols[jp] == i
        int jp = rowp[j];
        int jp_end = rowp[j] + diag_offset[j];
        while (jp < jp_end && cols[jp] < i){ 
          jp++; 
        }

        // x_sum_{i} = - L_{ij}*sum_{j}
        if (jp < jp_end && cols[jp] == i){
          TacsScalar * A = block_row_ptr[jp];
          
          TacsScalar alpha = -1.0, beta = 1.0;
          int one = 1;
          BLASgemv("N", &bj, &bi, &alpha, A, &bj,
                   &x[ni], &one, &beta, &x[nj], &one);
        }
      }
    }  

    // MPI_Reduce(x_local, &x[ni], bi, TACS_MPI_TYPE, 
    // MPI_SUM, diag_owner, comm);
    
    // MPI_Bcast(&x[ni], bi, TACS_MPI_TYPE, diag_owner, comm);
  }

  // Apply the upper part of the factorization
  // x_{i} <-- U_{i,i}^{-1} ( x_{i} - sum_{j}( U_{ij} * x_{j} ) )
  for ( int i = num_block_rows-1; i >= 0; i-- ){
    int ni = block_ptr[i];
    int bi = block_ptr[i+1] - ni;

    int diag_owner = get_block_owner(i, i);
    if (rank == diag_owner){
      memcpy(x_local, &x[ni], bi*sizeof(TacsScalar));
    }
    else {
      memset(x_local, 0, bi*sizeof(TacsScalar));
    }

    // U_{i,i} * x_sum_{i} = sum_{j}( - U_{ij} x_{j} )
    int jp_end = rowp[i] + row_len[i];
    int jp = rowp[i] + diag_offset[i]+1;
    for ( ; jp < jp_end; jp++ ){
      int j = cols[jp];
      if (rank == get_block_owner(i, j)){
	int nj = block_ptr[j];
	int bj = block_ptr[j+1] - nj;
	TacsScalar * A = block_row_ptr[jp];
	
	TacsScalar alpha = -1.0, beta = 1.0;
        int one = 1;
        BLASgemv("N", &bi, &bj, &alpha, A, &bi,
                 &x[nj], &one, &beta, x_local, &one);        
      }
    }

    // Reduce the sum to the diagonal owner
    MPI_Reduce(x_local, &x[ni], bi, TACS_MPI_TYPE, MPI_SUM, diag_owner, comm);
    
    // Solve U*x[ni] = x[ni] on the diagonal processor
    if (rank == diag_owner){
      TacsScalar * D = block_row_ptr[rowp[i] + diag_offset[i]];
      int incx = 1;
      BLAStrsv("U", "N", "N", &bi, D, &bi, &x[ni], &incx);
    }

    // Broadcast the result to the remaining processors
    MPI_Bcast(&x[ni], bi, TACS_MPI_TYPE, diag_owner, comm);
  }

  delete [] x_local;
}

/*
  Retrieve the pointer to the (i,j) block. If the block is not owned
  by this process or if the block is not in the non-zero pattern,
  return NULL.

  This function first checks that the processor owns the block, next
  it does a bisection search on the column j to find row i. If the
  entry is found, it returns the block pointer, otherwise, it returns
  NULL.
*/
TacsScalar * PDMatPivot::get_block( int rank, int i, int j ){
  TacsScalar * A = NULL;

  if (rank == get_block_owner(i, j)){
    // Look for row i in column j   
    int size = row_len[i];
    int * item = (int*)bsearch(&j, &cols[rowp[i]], size, 
                               sizeof(int), FElibrary::comparator);

    // If the entry is found, look up the storage location
    if (item){
      int ip = item - cols;
      A = block_row_ptr[ip];
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

  The L and U factors are defined such that,

  [ A[i,i]      |  A[i,i+1:n]     ]
  [ A[i+1:n,i]  |  A[i+1:n,i+1:n] ]
  = 
  [ L[i,i]      |   0             ][ U[i,i]  |  U[i,i+1:n]     ] 
  [ L[i+1:n,i]  |  L[i+1:n,i+1:n] ][ 0       |  U[i+1:n,i+1:n] ]

  As a result, the following relationships hold:

  L[i,i]*U[i,i] = A[i,i] 
  U[i,i+1:n]    = L[i,i]^{-1}*A[i,i+1:n]
  L[i+1:n,i]    = A[i+1:n,i]*U[i,i]^{-1}

  With the update:

  A[i+1:n,i+1:n] <-- A[i+1:n,i+1:n] - L[+1i:n,i]*U[i,i+1:n]

  The computational steps involved in the factorization are:
  
  For i = 1,...,n
  
  1. Compute the panel factorization by pivoting on rows in block i

  2. Apply the pivots to rows not in i, in turn generating new
  non-zero fill-ins.

  3. Compute the updates to the row and columns:
  L[i+1:n,i] = A[i+1:n,i]*U[i,i]^{-1}  and 
  U[i,i+1:n] = L[i,i]^{-1}*A[i,i+1:n]

  3. Compute the full update to the remainder of the matrix
  A[i+1:n,i+1:n] <-- A[i+1:n,i+1:n] - L[i+1:n,i]*U[i,i+1:n]
*/
void PDMatPivot::factor(){  
  int rank;
  MPI_Comm_rank(comm, &rank);

  if (proc_row < 0 || proc_col < 0){
    // This process is not on the process grid - does not participate
    return;
  }

  // Compute the absolute maximum size of the row and column buffers
  int max_row_buff_size = 0;
  for ( int i = proc_col; i < num_block_cols; i += npcols ){
    max_row_buff_size += block_ptr[i+1] - block_ptr[i];
  }
  max_row_buff_size *= max_block_size;

  int max_col_buff_size = 0;
  for ( int i = proc_row; i < num_block_rows; i += nprows ){
    max_col_buff_size += block_ptr[i+1] - block_ptr[i];
  }
  max_col_buff_size *= max_block_size;

  // Set the pivots array
  if (pivots){ delete [] pivots; }
  pivots = new int[ 2*block_ptr[num_block_rows] ];

  // Buffers to handle the recieves information
  TacsScalar * row_buff = new TacsScalar[max_row_buff_size];
  TacsScalar * col_buff = new TacsScalar[max_col_buff_size];
  TacsScalar * temp_diag = new TacsScalar[ max_block_size*max_block_size ];

  // Send information for rows owning U
  MPI_Request * row_send_request = new MPI_Request[nprows-1];
  MPI_Status * row_send_status = new MPI_Status[nprows-1];

  // Send information for columns owning L
  MPI_Request * col_send_request = new MPI_Request[npcols-1];
  MPI_Status * col_send_status = new MPI_Status[npcols-1];

  // The marker array denoting the non-zero entries in a row
  int * marker_array = new int[ num_block_cols ];
  memset(marker_array, 0, num_block_cols*sizeof(int));

  // The temporary storage array
  int * temp_array = new int[ max_block_size + num_block_cols ];

  // Keep track of the current column in the matrix
  int * current_col_offset = new int[ num_block_rows ];
  memset(current_col_offset, 0, num_block_rows*sizeof(int));
  int * col_rows = new int[ num_block_rows ];

  for ( int i = 0; i < num_block_rows; i++ ){
    // Get the block pointer
    int bi = block_ptr[i+1] - block_ptr[i];

    // The process row and column that actually own the
    // row and column in the 2D block cyclic structure
    int source_proc_row = get_proc_row(i);
    int source_proc_column = get_proc_column(i);

    // Determine which entries in this column have a non-zero
    // block associated with them. Each row with a non-zero block
    // is stored in col_rows[j] for 0 <= j < num_rows, and the offset
    // into the array is stored in current_col_offset[col_rows[j]]
    int num_rows = 0;
  
    for ( int j = i; j < num_block_rows; j++ ){
      int jp = rowp[j] + current_col_offset[j];
      while (cols[jp] < i){ jp++; }
      if (cols[jp] == i){
	col_rows[num_rows] = j;
	num_rows++;
      }
      current_col_offset[j] = jp - rowp[j];
    }

    // Get the owner of the diagonal block
    int diag_owner = get_block_owner(i, i);       
    TacsScalar * D = NULL;
    if (rank == diag_owner){
      D = block_row_ptr[rowp[i] + diag_offset[i]];
    }

    // Step 1: Panel factorization for column i
    if (proc_col == source_proc_column){
      factor_panel(i, col_rows, num_rows, row_buff,
		   current_col_offset);

      // Send the pivot sequence to the processors not in
      // the proc column
      int size = 2*(block_ptr[i+1] - block_ptr[i]);
      for ( int p = 0; p < npcols; p++ ){
	int dest = proc_grid[p + proc_row*npcols];
	if (rank != dest){
	  MPI_Send(&pivots[2*block_ptr[i]], size, MPI_INT,
		   dest, i, comm);
	}
      }
    }
    else {
      // Receive the pivot sequence from the proc column
      int source = proc_grid[source_proc_column + proc_row*npcols];
      int size = 2*(block_ptr[i+1] - block_ptr[i]);
      MPI_Recv(&pivots[2*block_ptr[i]], size, MPI_INT,
	       source, i, comm, MPI_STATUS_IGNORE); 
    }

    // Step 2: Update the non-zero pattern and swap the pivot sequence
    // Add the non-zero pattern from the swaps plus the GEMM update
    add_non_zero_pattern(i, &pivots[2*block_ptr[i]], bi, 
			 col_rows, num_rows, marker_array, temp_array);

    // Swap the rows in the remainder of the matrix
    int buff_offset = max_row_buff_size/max_block_size;
    swap_numerical_values(i, &pivots[2*block_ptr[i]], bi,
			  row_buff, &row_buff[buff_offset]);
   
    // Step 3: Compute the updates to the row/column i
    if (rank == diag_owner){
      // Send the factor to the column processes
      for ( int p = 0; p < npcols; p++ ){
	int dest = proc_grid[p + proc_row*npcols];
	if (rank != dest){
	  MPI_Send(D, bi*bi, TACS_MPI_TYPE, dest, p, comm);
        }
      }
    }
    
    // Receive the factorization of the diagonal: L[i,i]*U[i,i]
    if (proc_row == get_proc_row(i)){
      if (rank != diag_owner){
        MPI_Status status;
        MPI_Recv(temp_diag, bi*bi, TACS_MPI_TYPE,
                 diag_owner, proc_col, comm, &status);
        D = temp_diag;
      }
    }

    // Copy the required values from the column
    int col_buff_size = 0;
    if (proc_col == source_proc_column){
      // Loop through column i and for each non-zero block
      for ( int jp = 1; jp < num_rows; jp++ ){
	int j = col_rows[jp];
                
        if (proc_row == get_proc_row(j)){
          int bj = block_ptr[j+1] - block_ptr[j];
          int jpp = rowp[j] + current_col_offset[j];
          TacsScalar * A = block_row_ptr[jpp];

          memcpy(&col_buff[col_buff_size], A, bi*bj*sizeof(TacsScalar));
          col_buff_size += bi*bj;
	}
      }
    }
    else {
      // Compute the size of the incoming column
      for ( int jp = 1; jp < num_rows; jp++ ){
	int j = col_rows[jp];
        if (proc_row == get_proc_row(j)){
          col_buff_size += bi*(block_ptr[j+1] - block_ptr[j]);
        }
      }
    }

    // Apply L^{-1} to row i
    int row_buff_size = 0;
    if (proc_row == source_proc_row){
      // Loop through row i and for each non-zero block
      // compute L[i,i]^{-1}*A[i, j]
      int jp_end = rowp[i] + row_len[i];
      int jp = rowp[i] + diag_offset[i]+1;
      for ( ; jp < jp_end; jp++ ){
        int j = cols[jp];

        if (proc_col == get_proc_column(j)){
          int bj = block_ptr[j+1] - block_ptr[j];
          TacsScalar * A = block_row_ptr[jp];

          // Solve L[i,i]*A[i, j] = A[i, j]
          TacsScalar alpha = 1.0;
          BLAStrsm("L", "L", "N", "U", &bi, &bj, &alpha,
                   D, &bi, A, &bi);

          memcpy(&row_buff[row_buff_size], A, bi*bj*sizeof(TacsScalar));
          row_buff_size += bi*bj;
        }
      }
    }
    else {
      // Compute the size of the incoming row
      int jp_end = rowp[i] + row_len[i];
      int jp = rowp[i] + diag_offset[i]+1;
      for ( ; jp < jp_end; jp++ ){
	int j = cols[jp];

        if (proc_col == get_proc_column(j)){
          row_buff_size += bi*(block_ptr[j+1] - block_ptr[j]);
        }
      }
    }

    // Set up for step 4: Send and receive the rows to the column
    // processes and the column to the row processes.
    MPI_Request row_recv_request, col_recv_request;

    // Send the column buffer to all the rows
    if (proc_col == source_proc_column){
      for ( int p = 0, k = 0; p < npcols; p++ ){
        int dest = proc_grid[p + proc_row*npcols];
        if (rank != dest){
          int tag = 2*i+1;
          MPI_Isend(col_buff, col_buff_size, TACS_MPI_TYPE,
                    dest, tag, comm, &col_send_request[k]);
	  k++;
        }
      }
    }
    else {
      // The receiving processes
      int source = proc_grid[source_proc_column + proc_row*npcols];
      int tag = 2*i+1;
      MPI_Irecv(col_buff, col_buff_size, TACS_MPI_TYPE, 
                source, tag, comm, &col_recv_request);
    }

    // Send the row buffer to all the columns
    if (proc_row == source_proc_row){
      // The sending processes
      for ( int p = 0, k = 0; p < nprows; p++ ){
        int dest = proc_grid[proc_col + p*npcols];
        if (rank != dest){
          int tag = 2*i;
          MPI_Isend(row_buff, row_buff_size, TACS_MPI_TYPE,
                    dest, tag, comm, &row_send_request[k]);
	  k++;
        }
      }
    }
    else {
      // The receiving processes
      int source = proc_grid[proc_col + source_proc_row*npcols];
      int tag = 2*i;
      MPI_Irecv(row_buff, row_buff_size, TACS_MPI_TYPE, 
                source, tag, comm, &row_recv_request);
    }

    // Wait for the sends to complete
    if (source_proc_row == proc_row){
      MPI_Waitall(nprows-1, row_send_request, row_send_status);
    }
    if (source_proc_column == proc_col){
      MPI_Waitall(npcols-1, col_send_request, col_send_status);
    }
    
    // Wait for any receives to complete
    if (source_proc_column != proc_col){
      MPI_Status col_recv_status;
      MPI_Wait(&col_recv_request, &col_recv_status);
    }
    if (source_proc_row != proc_row){
      MPI_Status row_recv_status;
      MPI_Wait(&row_recv_request, &row_recv_status);
    }

    // Step 4: Apply the matrix GEMM update
    TacsScalar * L = col_buff;

    // Skip the diagonal entry by starting from iip = 1
    for ( int iip = 1; iip < num_rows; iip++ ){
      // Skip rows not locally owned
      int ii = col_rows[iip];
      if (proc_row == get_proc_row(ii)){

        // Compute the size of the ii-block
        int bii = block_ptr[ii+1] - block_ptr[ii];

        // Initialize the U pointer
        TacsScalar * U = row_buff;

	// Keep track of the destination row
	int kp = rowp[ii] + current_col_offset[ii];
	int kp_end = rowp[ii] + row_len[ii]; 

	// Set the pointer into the row
	int jjp = rowp[i] + diag_offset[i]+1;
        int jjp_end = rowp[i] + row_len[i];
        for ( ; jjp < jjp_end; jjp++ ){
          // Skip columns not locally owned
          int jj = cols[jjp];

	  // Increment the location of the destination until
	  // cols[kp] == jj
	  while (cols[kp] < jj && kp < kp_end){
	    kp++;
	  }

	  // If there are no more possible destination entries for this
	  // row, then we are done
	  if (kp == kp_end){
	    break;
	  }

	  // If the column is owned by this processor, we must
	  // increment the U-block, regardless of whether there is
	  // a destination block
          if (proc_col == get_proc_column(jj)){ 
            // Compute the size of the block
            int bjj = block_ptr[jj+1] - block_ptr[jj];

	    // If the destination matches, and this processor is the
	    // owner, then perform the GEMM update
	    if (jj == cols[kp] && 
		rank == get_block_owner(ii, jj)){
	      TacsScalar * A = block_row_ptr[kp];

              TacsScalar alpha = -1.0, beta = 1.0;
              BLASgemm("N", "N", &bii, &bjj, &bi,
                       &alpha, L, &bii, 
                       U, &bi, &beta, A, &bii);
            }
            
            // Increment the buffer pointer
            U += bi*bjj;
          }
        }

        // Increment the buffer pointer
        L += bi*bii;
      }
    }
  }

  // Release memory for the data transfer
  delete [] marker_array;
  delete [] temp_array;
  delete [] current_col_offset;
  delete [] col_rows;
  delete [] row_buff;
  delete [] col_buff;
  delete [] temp_diag;
  delete [] row_send_request;
  delete [] col_send_request;
  delete [] row_send_status;
  delete [] col_send_status;
}

/*
  Factor the given panel with partial pivoting.

  This function computes the factorization of a panel: a portion of
  the matrix that has a much smaller width than depth. This function
  effectively computes the following:

  [ D ] = [ L_{D} ]*U_{D} 
  [ B ]   [ L_{B} ]

  where L_{D}*U_{D} is the diagonal factorization, B is the matrix
  below the diagonal, and L_{B} is the lower diagonal after the
  factorization. Note that this formula does not reflect the partial
  pivoting that this function performs.

  This function generates the pivot sequence for the partial pivoting.
  The panel factorization is simplified by the fact that all the rows
  that participate have the same number of column entries (they span
  the column), therefore no new non-zero entries are generated during
  this stage. Instead, new non-zero entries are generated when
  completing the pivot sequence on the remainder of the rows outside
  the column block.
  
  input:
  col:      the fully summed column that can be factorized nex
  col_rows: the non-zero rows in the column 
  num_rows: the number of rows in the column
  temp:     a temporary buffer twice the panel width

  side-effects:
  Note that this function inserts the pivot sequence locally into
  the member pivots array
*/
void PDMatPivot::factor_panel( int col, 
                               const int * col_rows, int num_rows,
                               TacsScalar * temp,
			       const int * current_col_offset ){
  // Check if the processor is in the column
  int rank;
  MPI_Comm_rank(comm, &rank);

  MatBlockIndex max_index;
  int diag_owner = get_block_owner(col, col);

  // Set the diagonal amtrix if this is the diagonal pointer
  TacsScalar * D = NULL;
  if (rank == diag_owner){
    D = block_row_ptr[rowp[col] + diag_offset[col]];
  }

  // Compute the width of the panel
  int bi = block_ptr[col+1] - block_ptr[col];
  
  for ( int i = 0; i < bi; i++ ){
    if (i > 0){
      if (rank == diag_owner){
        // Solve L[:i, :i]*y = D[i]
        int incx = 1;
        BLAStrsv("L", "N", "U", &i, D, &bi, &D[bi*i], &incx);

        // Copy the result into temp
        memcpy(temp, &D[bi*i], i*sizeof(TacsScalar));

        // Send the result to the remaining processors in 
        // this column
        for ( int p = 0; p < nprows; p++ ){          
          int proc = proc_grid[proc_col + p*npcols];
          if (proc != rank){
            MPI_Send(temp, i, TACS_MPI_TYPE, proc, i, comm);            
          }
        }

        // Compute the update to the next column i from this
        // new portion of the factorization
        TacsScalar alpha = -1.0, beta = 1.0;
        int size = bi - i;
        BLASgemv("N", &size, &i, &alpha, &D[i], &bi, 
                 &D[bi*i], &incx, &beta, &D[(bi+1)*i], &incx);
      }
      else {
        // Receive the result into the array temp
        MPI_Recv(temp, i, TACS_MPI_TYPE, diag_owner, i, 
                 comm, MPI_STATUS_IGNORE);
      }

      // Start in the row after the diagonal pivot
      for ( int jp = 1; jp < num_rows; jp++ ){
        // Compute the row to perform the update on
        int j = col_rows[jp];

        // Perform the update on this block
        if (proc_row == get_proc_row(j)){
          int jpp = rowp[j] + current_col_offset[j];
          TacsScalar * A = block_row_ptr[jpp];
          int bj = block_ptr[j+1] - block_ptr[j];
          
          // Compute the update to the next column on this
          // block
          TacsScalar alpha = -1.0, beta = 1.0;
          int incx = 1;
          BLASgemv("N", &bj, &i, &alpha, A, &bj, 
                   temp, &incx, &beta, &A[bj*i], &incx);
        }
      }
    }

    // Compute the maximum index on each processor
    TacsScalar * A_max = NULL;
    max_index.entry = 0.0;
    max_index.block_index = -1;
    max_index.index = -1;
    max_index.rank = rank;

    // First, deal with the diagonal entry
    if (rank == diag_owner){
      TacsScalar * d = &D[(bi+1)*i];

      A_max = D;
      max_index.entry = d[0];
      max_index.block_index = col;
      max_index.index = i;
      d++;

      for ( int k = i+1; k < bi; k++ ){
        if (fabs(RealPart(d[0])) > fabs(RealPart(max_index.entry))){
          max_index.entry = d[0];
          max_index.index = k;
        }
        d++;
      }
    }
    
    // Find the maximum entry that is stored locally
    // on this processor within this column
    for ( int jp = 1; jp < num_rows; jp++ ){
      // Compute the row to perform the update on
      int j = col_rows[jp];

      // Perform the update on this block
      if (proc_row == get_proc_row(j)){
        int jpp = rowp[j] + current_col_offset[j];
        TacsScalar * A = block_row_ptr[jpp];
        int bj = block_ptr[j+1] - block_ptr[j];
        
        TacsScalar * a = &A[bj*i];
        for ( int k = 0; k < bj; k++ ){
          if (fabs(RealPart(a[0])) > fabs(RealPart(max_index.entry))){
            A_max = A;
            max_index.entry = a[0];
            max_index.block_index = j;
            max_index.index = k;
          }
          a++;
        }
      }
    }

    // Now that we have the local maximum entry, perform a reduction
    // across all the processors in this column by sending the result
    // to the root 
    if (rank == diag_owner){
      MatBlockIndex temp;
      for ( int p = 0; p < nprows; p++ ){
        int proc = proc_grid[proc_col + p*npcols];
        if (rank != proc){
          MPI_Recv(&temp, 1, mat_block_type, proc, i, 
                   comm, MPI_STATUS_IGNORE);

          // Check if which is the larger entry
          if (fabs(RealPart(temp.entry)) > 
              fabs(RealPart(max_index.entry))){
            max_index.entry = temp.entry;
            max_index.block_index = temp.block_index;
            max_index.index = temp.index;
            max_index.rank = temp.rank;
          }
        }
      }

      for ( int p = 0; p < nprows; p++ ){
        int proc = proc_grid[proc_col + p*npcols];
        if (rank != proc){
          MPI_Send(&max_index, 1, mat_block_type, proc, i, comm);
        }
      }
    }
    else {
      // Send the result to the root processor, wait for it to receive
      // all the results, compute the maximum and return the result back
      MPI_Send(&max_index, 1, mat_block_type, diag_owner, i, comm);
      MPI_Recv(&max_index, 1, mat_block_type, diag_owner, i, 
               comm, MPI_STATUS_IGNORE);
    }

    // Set the pivot sequence into the global array
    pivots[2*(block_ptr[col]+i)] = max_index.block_index;
    pivots[2*(block_ptr[col]+i)+1] = max_index.index;

    // Send the result to the diagonal owner, 
    // First, check to see if the maximum entry is on the same processor
    // If it is, this is fortuitous and we don't need to perform any 
    // communication
    int max_rank = max_index.rank;
    if (diag_owner == max_rank){
      if (rank == diag_owner &&
          !(max_index.block_index == col &&
            max_index.index == i)){
        // Swap the values directly since they are on the same processor
        int index = max_index.index;
        TacsScalar * a = &A_max[index];
	int j = max_index.block_index;
	int bj = block_ptr[j+1] - block_ptr[j];
        TacsScalar * d = &D[i];

        for ( int k = 0; k < bi; k++ ){
          TacsScalar t = a[0];
          a[0] = d[0];
          d[0] = t;
          d += bi;
          a += bj;
        }
      }
    }
    else {
      // Unfortunately, the max entry is not on the same processor, therefore,
      // we need to communicate the max row in the panel to the diagonal owner
      // and the diagonal owner row to the former max entry row
      if (diag_owner == rank){
        // Post the recieve request
        MPI_Request send_request, recv_request;
        MPI_Irecv(&temp[bi], bi, TACS_MPI_TYPE, max_rank, i, 
                  comm, &recv_request);

        // Copy the current from the diagonal entry
        TacsScalar * t = &temp[0];
        TacsScalar * d = &D[i];
        for ( int k = 0; k < bi; k++ ){
          t[0] = d[0];
          d += bi;  t++;
        }

        // Send the buffer and wait for the receive to complete
        MPI_Isend(&temp[0], bi, TACS_MPI_TYPE, max_rank, i, 
                  comm, &send_request);
        MPI_Wait(&recv_request, MPI_STATUS_IGNORE);

        // Copy the recieved array back into the diagonal entry
        t = &temp[bi];
        d = &D[i];
        for ( int k = 0; k < bi; k++ ){
          d[0] = t[0];
          d += bi;  t++;
        }

        // Wait for the send to finish
        MPI_Wait(&send_request, MPI_STATUS_IGNORE);
      }
      if (rank == max_rank){
	// Compute the block size
	int j = max_index.block_index;
	int bj = block_ptr[j+1] - block_ptr[j];

        // Post the recieve request
        MPI_Request send_request, recv_request;
        MPI_Irecv(&temp[bi], bi, TACS_MPI_TYPE, diag_owner, i, 
                  comm, &recv_request);

        // Copy the max row into the send buffer
        int index = max_index.index;
        TacsScalar * a = &A_max[index];
        TacsScalar * t = &temp[0];
        for ( int k = 0; k < bi; k++ ){
          t[0] = a[0];
          a += bj;  t++;
        }

        // Send the buffer and wait for the receive to complete
        MPI_Isend(&temp[0], bi, TACS_MPI_TYPE, diag_owner, i, 
                  comm, &send_request);
        MPI_Wait(&recv_request, MPI_STATUS_IGNORE);

        // Copy the recieve into the original buffer
	a = &A_max[index];
	t = &temp[bi];
        for ( int k = 0; k < bi; k++ ){
          a[0] = t[0];
          a += bj;  t++;
        }        

        // Wait for the send to finish
        MPI_Wait(&send_request, MPI_STATUS_IGNORE);
      }
    }

    // Now, divide the entries in column i below the diagonal 
    // by the maximum entry that is now on the diagonal.
    // First, deal with the diagonal entry:
    if (rank == diag_owner){
      TacsScalar * d = &D[i*bi + i+1];
      for ( int k = i+1; k < bi; k++ ){
        d[0] = d[0]/max_index.entry;
        d++;
      }
    }

    // Divide the remainder of the column by the index
    // Again, skip the first entry (the diagonal)
    for ( int jp = 1; jp < num_rows; jp++ ){
      // Compute the row to perform the update on
      int j = col_rows[jp];
      
      // Perform the update on this block
      if (proc_row == get_proc_row(j)){
        int jpp = rowp[j] + current_col_offset[j];
        TacsScalar * A = block_row_ptr[jpp];
        int bj = block_ptr[j+1] - block_ptr[j];
        
        TacsScalar * a = &A[bj*i];
        for ( int k = 0; k < bj; k++ ){
          a[0] = a[0]/max_index.entry;
          a++;
        }
      }
    }
  }
}

/*
  Compress the sparse matrix storage format to maximize the amount
  of avaiable space for the next row.

  This code moves the current space to the end of the array to
  facilitate adding entries during the merging of two rows. If enough
  free space is not made avaiable during this compression, then
  additional space might need to be allocated to the sparsity pattern
  storage scheme (but not the sparse matrix data itself.)

  This code modifies the values of rowp[row] and moves
  all the column indices so that rowp[row] = rowp[row-1] +
  row_len[row]. This equation is not necessarily satisfied by the
  non-compressed rows. This function also updates the block_row_ptr
  array to maintain consistency with the original data. By calling
  this function for each row, the code will have the maximum amount of
  space available at the during the computation.  Note that row_len
  remains unchanged after this operation.

  input:
  row:  the row to compress
*/
void PDMatPivot::compress_row( int row ){
  if (row < num_block_rows-1){
    int old_rowp = rowp[row+1];
    rowp[row+1] = rowp[row] + row_len[row];
    
    int * r = &cols[old_rowp];
    int * n = &cols[rowp[row+1]];
    TacsScalar ** rd = &block_row_ptr[old_rowp];
    TacsScalar ** nd = &block_row_ptr[rowp[row+1]];

    for ( int i = 0; i < row_len[row+1]; i++ ){
      n[0] = r[0];    n++;  r++;
      nd[0] = rd[0];  nd++; rd++;
    }
    
    // Now, assign null values to the remaining components
    for ( int ip = rowp[row+1] + row_len[row+1]; ip < rowp[row+2]; ip++ ){
      cols[ip] = -1;
      block_row_ptr[ip] = NULL;
    }
  }
}

/*
  Insert new entries into the non-zero pattern.

  This assumes that the new entries are not duplicates of the entries
  that exist in this array already and that the entries in the
  new_cols array are sorted in ascending order. These requirements
  make adding the additional rows easier and less cumbersome, and
  likely more efficient - although I make no claim that this is the
  best possible way to do this! 

  For each new non-zero added to the data structure, additional space
  is allocated to store the corresponding non-zero matrix
  entries. These additional entries are automatically zeroed upon
  their allocation. In addition, the diag pointer is updated and the
  block_row_ptr is updated as well.

  input:
  row:      the row to add the new entries to
  new_cols: the new columns to add to the non-zero pattern
  num_nnz:  the number of non-zero columns to add
*/
void PDMatPivot::insert_row_nonzeros( int row, const int * new_cols,
                                      int num_nnz ){
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Set the new length of the array
  int bi = block_ptr[row+1] - block_ptr[row];
  int col_ptr = rowp[row] + row_len[row]-1;
  int ptr = rowp[row] + row_len[row] + num_nnz-1;

  // Check if updating the row length will exceed the total length
  // of the array - if it will, then we need to either compress the
  // existing data structure or allocate new storage space
  int new_space_req = num_nnz - (rowp[row+1] - rowp[row] - row_len[row]);
  
  // The search for space proceeds in two stages:
  // First, we try and partially compress the next few rows to 
  // free up enough space to store data required for the new
  // non-zeros. If that does not create enough room, we re-allocate
  // the arrays and copy over the results into the data structure.
  if (new_space_req > 0){
    // We need to leave at least this much extra space per row
    // so that we have enough space for the future
    int min_space_per_row = 5; 

    // A running total of the space we could create by compression
    int compress_space = 0;

    // First, count up whether it is feasible to obtain the required space
    int final_row = row+1; 
    for ( ; final_row < num_block_rows; final_row++ ){
      int extra_space = rowp[final_row+1] - rowp[final_row] - row_len[final_row];
      if (extra_space > min_space_per_row){
	compress_space += extra_space - min_space_per_row;
      }
      if (compress_space >= new_space_req){
	break;
      }
    }
    
    // We have enough space to compress and avoid reallocation for now
    if (compress_space >= new_space_req){
      // Compress the data structure 
      int new_ptr = rowp[final_row];
      int extra_space = (rowp[final_row+1] - 
			 rowp[final_row] - row_len[final_row]);
      
      // Adjust the new pointer position
      if (extra_space > min_space_per_row){
	new_ptr += extra_space - min_space_per_row;
      }

      // First check if additional values should be copied 
      for ( int i = final_row; i > row; i-- ){
	// Set the new 
	for ( int k = row_len[i]-1; k >= 0; k-- ){
	  cols[new_ptr + k] = cols[rowp[i] + k];
	  block_row_ptr[new_ptr + k] = block_row_ptr[rowp[i] + k];
	}

	// Calculate the extra space in the next row
	extra_space = rowp[i] - rowp[i-1] - row_len[i-1];
	int temp = new_ptr;
	new_ptr -= (rowp[i] - rowp[i-1]);

	// Set the new location of the last pointer
	rowp[i] = temp;

	// Adjust the new pointer position to account for the extra space
	if (extra_space > min_space_per_row){
	  new_ptr += extra_space - min_space_per_row;
	}
      }

      // Set the space request to zero - since we've now added
      // enough new space
      new_space_req = 0;
    }

    // The compression didn't work, we need to allocate more memory and
    // expand the existing arrays to match
    if (new_space_req > 0){
      int new_len = rowp[num_block_rows] + new_space_req*(num_block_rows - row);

      // Allocate temporary arrays for the storage
      int * new_cols = new int[ new_len ];
      TacsScalar ** new_br_ptr = new TacsScalar*[ new_len ];

      // Copy over the existing data
      memcpy(new_cols, cols, rowp[row]*sizeof(int));
      memcpy(new_br_ptr, block_row_ptr, rowp[row]*sizeof(TacsScalar*));

      // In each row, copy over the values in [rowp[i]:rowp[i]+row_len[i]]
      int ptr = rowp[row];
      for ( int i = row; i < num_block_rows; i++ ){
	for ( int j = 0; j < row_len[i]; j++ ){
	  new_cols[ptr + j] = cols[rowp[i] + j];
	  new_br_ptr[ptr + j] = block_row_ptr[rowp[i] + j];
	}
	int temp = ptr;
	ptr += rowp[i+1] - rowp[i] + new_space_req;
	rowp[i] = temp;
      }
      rowp[num_block_rows] = ptr;

      // Delete the old data and set the pointers to the new arrays
      delete [] cols;
      delete [] block_row_ptr;

      cols = new_cols;
      block_row_ptr = new_br_ptr;
    }
  }

  // Update the length of the current row
  row_len[row] += num_nnz;

  // Set the num_nnz to point to the last entry
  num_nnz--;
  
  while (num_nnz >= 0 && col_ptr >= rowp[row]){
    // Copy the maximum entry from the old array to the
    // new array
    if (cols[col_ptr] > new_cols[num_nnz]){
      cols[ptr] = cols[col_ptr];
      block_row_ptr[ptr] = block_row_ptr[col_ptr];

      // Check if this is the diagonal entry and update it
      // if it is required
      if (cols[ptr] == row){
        diag_offset[row] = ptr - rowp[row];
      }

      ptr--;  col_ptr--;
    }
    else {
      int col = new_cols[num_nnz];
      cols[ptr] = col;

      if (rank == get_block_owner(row, col)){
        // Check if allocating a block of bi*bj will exceed
        // the available free space in the current block
        int bj = block_ptr[col+1] - block_ptr[col];
        int current_size = data_top->data_ptr - data_top->data;
        if (current_size + bi*bj > data_top->data_size){
          // Allocate new space for the non-zero terms
          int new_size = max_block_size*max_block_size*num_new_blocks;
          data_top->next = new PDMatData(new_size);
          data_top = data_top->next;
        }

        block_row_ptr[ptr] = data_top->data_ptr;
        data_top->data_ptr += bi*bj;
      }
      else {
        block_row_ptr[ptr] = NULL;
      }

      ptr--;  num_nnz--;
    }
  }

  // Finish off the row - this is only required when there are 
  // additional terms that have to be added to the non-zero pattern
  // at the beginning of the cols[rowp[row]] array, otherwise, the
  // old pattern can remain in place
  while (num_nnz >= 0){
    int col = new_cols[num_nnz];
    cols[ptr] = col;

    if (rank == get_block_owner(row, col)){
      // Check if allocating a block of bi*bj will exceed
      // the available free space in the current block
      int bj = block_ptr[col+1] - block_ptr[col];
      int current_size = data_top->data_ptr - data_top->data;
      if (current_size + bi*bj > data_top->data_size){
        // Allocate new space for the non-zero terms
        int new_size = max_block_size*max_block_size*num_new_blocks;
        data_top->next = new PDMatData(new_size);
        data_top = data_top->next;
      }
      
      block_row_ptr[ptr] = data_top->data_ptr;
      data_top->data_ptr += bi*bj;
    }
    else {
      block_row_ptr[ptr] = NULL;
    }

    ptr--;  num_nnz--;
  }
}

/*
  This function adds the non-zero pattern to the matrix requried
  for the swapping of rows and the trailing-matrix update.

  This function works in two stages: first, the effect of the swap on
  the symbolic structure is computed, based on the non-zero entries in
  the rows participating in the swap. Second, the non-zero pattern of
  the pivot rows are determined and the effect on the trailing matrix is
  computed.

  The non-zero pattern in the pivot row is the union of the initial
  non-zero pattern and the non-zero patterns of all rows that swap
  with this row. The non-zero pattern of the trailing matrix update
  can be thought of as the result of an outer product of the current
  row and column and must take into consideration the fill-ins due to
  pivoting.

  Note that this approach makes it very difficult to swap just the
  original numerical values.  Instead, the entire non-zero pattern is
  swapped for two participating rows - which is much easier.  This
  results in an additional overhead due to the zero values that are
  passed between processors.

  Note that a conservative estimate for the size of the temporary
  integer array is: max_block_size + num_block_cols

  input:
  row:           the row that we're working on
  pivot:         the block row/index pairs for each pivot
  num_pivots:    the number of rows we're swapping (local block size)
  col_rows:      the row indices in the current column
  num_col_rows:  the number of non-zero block rows in the current column
  marker:        an array used to mark entries to determine the nz patter
  temp:          a temporary integer array 
*/
void PDMatPivot::add_non_zero_pattern( int row, 
				       const int * pivot,
				       int num_pivots,
				       const int * col_rows, 
				       int num_col_rows,
				       int * marker, int * temp ){
  // "Uniquify" the rows so we don't have to do things twice (or more)
  int * sorted_rows = &temp[0];
  for ( int i = 0; i < num_pivots; i++ ){
    sorted_rows[i] = pivot[2*i];
  }
  int num_unique_rows = FElibrary::uniqueSort(sorted_rows, num_pivots);

  // Add the entries from the current row
  int ip_end = rowp[row] + row_len[row];
  for ( int ip = rowp[row]; ip < ip_end; ip++ ){
    int col = cols[ip];
    marker[col] = row+1;
  }

  // Keep track of the new columns
  int num_new_cols = 0;
  int * new_cols = &temp[num_pivots];

  // Keep track of how many columns will be added
  for ( int i = 0; i < num_unique_rows; i++ ){
    int new_row = sorted_rows[i];
    int jp_end = rowp[new_row] + row_len[new_row];
    int jp = rowp[new_row];
    while (jp < jp_end && cols[jp] <= row){ 
      jp++;
    }

    for ( ; jp < jp_end; jp++ ){
      int col = cols[jp];
      if (marker[col] != row+1){
        marker[col] = row+1;
        new_cols[num_new_cols] = col;
        num_new_cols++;
      }
    }
  }

  // Add the new rows to the non-zero pattern
  qsort(new_cols, num_new_cols, sizeof(int), FElibrary::comparator);
  if (num_new_cols > 0){
    insert_row_nonzeros(row, new_cols, num_new_cols);
  }

  // Now, loop through and find which columns must be added to each
  // new row. Note that we start the loop after the first entry
  // to skip the initial row.
  for ( int jp = 1; jp < num_col_rows; jp++ ){
    int j = col_rows[jp];

    // Keep track of the number of new columns to be added to row j
    num_new_cols = 0;
    
    int ip = rowp[row] + diag_offset[row]+1;
    int ip_end = rowp[row] + row_len[row];
    
    int kp = rowp[j];
    int kp_end = rowp[j] + row_len[j];
    for ( ; kp < kp_end && ip < ip_end; kp++ ){
      // We need to add the column if it is in cols[ip] but
      // not in cols[kp]
      while (ip < ip_end && 
	     cols[ip] < cols[kp]){
	new_cols[num_new_cols] = cols[ip];
	num_new_cols++;
	ip++;
      }
      
      // If the entries are equal, increment ip without
      // adding cols[ip] to the new entries
      if (ip < ip_end && cols[ip] == cols[kp]){
        ip++;
      }
    }

    while (ip < ip_end){
      new_cols[num_new_cols] = cols[ip];
      num_new_cols++;
      ip++;
    }

    if (num_new_cols > 0){
      insert_row_nonzeros(j, new_cols, num_new_cols);
    }
  }

  // Compress the row that we've just computed
  compress_row(row);
}

/*
  Swap the numerical values in the row in sequence. 
  
  This function performs the sequence of row swaps for for the rows
  that are in the trailing matrix (in python these are the entries in
  A[row:, row+1:].)  Note that the column is excluded since these
  swaps have already been performed when computing the pivot sequence
  to begin with.

  This code could be improved by permuting the swaps before hand,
  instead of computing each in sequence which would lead to better
  parallel performance.

  input:
  row:         the current row
  pivot:       the array of row/index pairs for the pivots
  num_pivots:  the number of pivots to perform
  send_buff:   temporary buffer to store the sends
  recv_buff:   temporary buffer to store the receives
*/
void PDMatPivot::swap_numerical_values( int row,
					const int * pivot,
					int num_pivots, 
					TacsScalar * send_buff, 
					TacsScalar * recv_buff ){
  // Determine the size of the buffers that are requried
  int buff_size = 0;
  int jp = rowp[row] + diag_offset[row] + 1;
  int jp_end = rowp[row] + row_len[row];
  for ( ; jp < jp_end; jp++ ){
    int j = cols[jp];
    if (proc_col == get_proc_column(j)){
      buff_size += block_ptr[j+1] - block_ptr[j];
    }
  }

  // Get the rank of this row
  int row_rank = get_block_owner(row, proc_col);
  
  // Compute the size of the block row
  int bi = block_ptr[row+1] - block_ptr[row];

  int rank;
  MPI_Comm_rank(comm, &rank);

  for ( int i = 0; i < num_pivots; i++ ){
    // Read in the pivot row and pivot index
    int piv_row = pivot[2*i];
    int index = pivot[2*i+1];

    // Compute the leading dimension of the pivot row
    int bpiv = block_ptr[piv_row+1] - block_ptr[piv_row];

    // Check where the rows need to be exchanged
    if (proc_row == get_proc_row(row) && 
	proc_row == get_proc_row(piv_row)){
      // Everything is on the same group of processors - no need for MPI
      if (i != index || piv_row != row){
	// Find the starting location for the row
	int kp_start = rowp[piv_row];
	int kp_end = kp_start + row_len[piv_row];
	while (kp_start < kp_end &&
	       cols[kp_start] <= row){
	  kp_start++;
	}
	// Now cols[kp_start] > row
	
	// Swap the rows
	int jp = rowp[row] + diag_offset[row] + 1;
	int jp_end = rowp[row] + row_len[row];
	for ( int kp = kp_start; 
	      (jp < jp_end && kp < kp_end); jp++, kp++ ){
	  // The column structure must be the same
	  // between the two rows for cols[] > row
	  int j = cols[jp];
	  if (proc_col == get_proc_column(j)){
	    int bj = block_ptr[j+1] - block_ptr[j];
	    TacsScalar * A = block_row_ptr[jp];
	    TacsScalar * a = &A[i];

	    TacsScalar * B = block_row_ptr[kp];
	    TacsScalar * b = &B[index];
	    
	    // Copy to the send buffer
	    for ( int k = 0; k < bj; k++ ){
	      TacsScalar t = a[0];
	      a[0] = b[0];
	      b[0] = t;
	      a += bi;
	      b += bpiv;
	    }
	  }
	}
      }
    }
    else {
      int piv_rank = get_block_owner(piv_row, proc_col);
      if (proc_row == get_proc_row(row)){
	// Post the recieve request
	MPI_Request send_request, recv_request;
	MPI_Irecv(recv_buff, buff_size, TACS_MPI_TYPE, piv_rank, i,
		  comm, &recv_request);
	
	// Copy all the values from the row to the send buffer
	TacsScalar * s = send_buff;
	int jp = rowp[row] + diag_offset[row]+1;
	int jp_end = rowp[row] + row_len[row];
	for ( ; jp < jp_end; jp++ ){
	  int j = cols[jp];
	  
	  if (proc_col == get_proc_column(j)){
	    int bj = block_ptr[j+1] - block_ptr[j];
	    TacsScalar * A = block_row_ptr[jp];
	    TacsScalar * a = &A[i];
	    
	    // Copy to the send buffer
	    for ( int k = 0; k < bj; k++ ){
	      s[0] = a[0];
	      a += bi;
	      s++;
	    }
	  }
	}
	
        // Send the buffer and wait for the receive to complete
        MPI_Isend(send_buff, buff_size, TACS_MPI_TYPE, piv_rank, i, 
                  comm, &send_request);
        MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
	
	// Receieve the result and place into row i locally
	TacsScalar * r = recv_buff;
	
	// Set the jp pointer to the entry after the diagonal
	jp = rowp[row] + diag_offset[row] + 1;
	for ( ; jp < jp_end; jp++ ){
	  int j = cols[jp];

	  if (proc_col == get_proc_column(j)){
	    int bj = block_ptr[j+1] - block_ptr[j];
	    TacsScalar * A = block_row_ptr[jp];
	    TacsScalar * a = &A[i];
	    
	    // Copy to the send buffer
	    for ( int k = 0; k < bj; k++ ){
	      a[0] = r[0];
	      a += bi;
	      r++;
	    }
	  }
	}

	MPI_Wait(&send_request, MPI_STATUS_IGNORE);
      }
      if (proc_row == get_proc_row(piv_row)){
	// Post the recieve request
	MPI_Request send_request, recv_request;
	MPI_Irecv(recv_buff, buff_size, TACS_MPI_TYPE, row_rank, i,
		  comm, &recv_request);

	// Find the offset into the row we're taking the pivot from
	int jp_start = rowp[piv_row];
	int jp_end = jp_start + row_len[piv_row];

	// Find the location such that jp_start > row
	while (jp_start < jp_end && cols[jp_start] <= row){
	  jp_start++;
	}

	// Compute the end location of the current row
	int kp_end = rowp[row] + row_len[row];

	// Copy all the values from the row to the send buffer
	TacsScalar * s = send_buff;
	int kp = rowp[row] + diag_offset[row]+1;
	for ( int jp = jp_start; kp < kp_end; kp++ ){
	  int j = cols[kp];
	  if (proc_col == get_proc_column(j)){
	    // Increment jp until we have j = cols[jp]
	    while (jp < jp_end && cols[jp] < j){
	      jp++;
	    }

	    if (jp < jp_end && cols[jp] == j){
	      int bj = block_ptr[j+1] - block_ptr[j];
	      TacsScalar * A = block_row_ptr[jp];
	      TacsScalar * a = &A[index];
	  
	      // Copy to the send buffer
	      for ( int k = 0; k < bj; k++ ){
		s[0] = a[0];
		a += bpiv;
		s++;
	      }
	    }
	  }
	}

        // Send the buffer and wait for the receive to complete
        MPI_Isend(send_buff, buff_size, TACS_MPI_TYPE, row_rank, i, 
                  comm, &send_request);
        MPI_Wait(&recv_request, MPI_STATUS_IGNORE);

	// Receieve the result and place into row i locally
	TacsScalar * r = recv_buff;
	kp = rowp[row] + diag_offset[row]+1;
	for ( int jp = jp_start; jp < jp_end; jp++ ){
	  int j = cols[kp];
	  if (proc_col == get_proc_column(j)){
	    // Increment jp until we have j = cols[jp]
	    while (jp < jp_end && cols[jp] < j){
	      jp++;
	    }

	    if (jp < jp_end && cols[jp] == j){
	      int bj = block_ptr[j+1] - block_ptr[j];
	      TacsScalar * A = block_row_ptr[jp];
	      TacsScalar * a = &A[index];
	    
	      // Copy to the send buffer
	      for ( int k = 0; k < bj; k++ ){
		a[0] = r[0];
		a += bpiv;
		r++;
	      }
	    }
	  }
	}

	MPI_Wait(&send_request, MPI_STATUS_IGNORE);
      }
    }
  }
}
