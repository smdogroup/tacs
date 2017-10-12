#include "BCSRMatImpl.h"

/*
  The following file contains the specific implementation for block size = 3.

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

/*!
  Compute the matrix-vector product: y = A * x
*/
void BCSRMatVecMult3( BCSRMatData * data,
                      TacsScalar * x, TacsScalar * y ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;
  const TacsScalar * a = data->A;

  for ( int i = 0; i < nrows; i++ ){
    y[0] = 0.0;
    y[1] = 0.0;
    y[2] = 0.0;

    int end = rowp[i+1];
    for ( int k = rowp[i]; k < end; k++ ){
      int j = 3*cols[k];
        
      y[0] += a[0]*x[j] + a[1]*x[j+1] + a[2]*x[j+2];
      y[1] += a[3]*x[j] + a[4]*x[j+1] + a[5]*x[j+2];
      y[2] += a[6]*x[j] + a[7]*x[j+1] + a[8]*x[j+2];
      a += 9;
    }
    y += 3;
  }
}

/*!
  Compute the matrix vector product plus addition: z = A * x + y
*/
void BCSRMatVecMultAdd3( BCSRMatData * data,
                         TacsScalar * x, TacsScalar * y, TacsScalar * z ){

  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;
  const TacsScalar * a = data->A;

  for ( int i = 0; i < nrows; i++ ){
    y[0] = z[0];
    y[1] = z[1];
    y[2] = z[2];

    int end = rowp[i+1];
    for ( int k = rowp[i]; k < end; k++ ){
      int j = 3*cols[k];
        
      y[0] += a[0]*x[j] + a[1]*x[j+1] + a[2]*x[j+2];
      y[1] += a[3]*x[j] + a[4]*x[j+1] + a[5]*x[j+2];
      y[2] += a[6]*x[j] + a[7]*x[j+1] + a[8]*x[j+2];
      a += 9;
    }
    y += 3;
    z += 3;
  }
}

/*!
  Apply the lower factorization y = L^{-1} x
*/
void BCSRMatApplyLower3( BCSRMatData * data,
                         TacsScalar * x, TacsScalar * y ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;
  const int * diag = data->diag;
  const TacsScalar * A = data->A;

  TacsScalar * z = y;

  for ( int i = 0; i < nrows; i++ ){
    z[0] = x[0];
    z[1] = x[1];
    z[2] = x[2];

    int end = diag[i];
    int k   = rowp[i];
    const TacsScalar * a = &A[9*k];
    for ( ; k < end; k++ ){
      int j = 3*cols[k];

      z[0] -= a[0]*y[j] + a[1]*y[j+1] + a[2]*y[j+2];
      z[1] -= a[3]*y[j] + a[4]*y[j+1] + a[5]*y[j+2];
      z[2] -= a[6]*y[j] + a[7]*y[j+1] + a[8]*y[j+2];
      a += 9;
    }
    z += 3;
    x += 3;
  }
}

/*!
  Apply the upper factorization y = U^{-1} x
*/
void BCSRMatApplyUpper3( BCSRMatData * data,
                         TacsScalar * x, TacsScalar * y ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;
  const int * diag = data->diag;
  const TacsScalar * A = data->A;

  TacsScalar y0, y1, y2;

  x = &x[3*(nrows-1)];
  for ( int i = nrows-1; i >= 0; i-- ){
    y0 = x[0];
    y1 = x[1];
    y2 = x[2];

    int end = rowp[i+1];
    int k   = diag[i]+1;
    const TacsScalar * a = &A[9*k];
    for ( ; k < end; k++ ){
      int j = 3*cols[k];

      y0 -= a[0]*y[j] + a[1]*y[j+1] + a[2]*y[j+2];
      y1 -= a[3]*y[j] + a[4]*y[j+1] + a[5]*y[j+2];
      y2 -= a[6]*y[j] + a[7]*y[j+1] + a[8]*y[j+2];
      a += 9;
    }

    int bi = 3*i;
    a = &A[9*diag[i]];
    y[bi  ] = a[0]*y0 + a[1]*y1 + a[2]*y2;
    y[bi+1] = a[3]*y0 + a[4]*y1 + a[5]*y2;
    y[bi+2] = a[6]*y0 + a[7]*y1 + a[8]*y2;
    
    x -= 3;
  }
}

/*!
  Apply a portion of the lower factorization x = L^{-1} x
*/

void BCSRMatApplyPartialLower3( BCSRMatData * data, TacsScalar * x, 
                                int var_offset ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;
  const int * diag = data->diag;
  const TacsScalar * A = data->A;

  TacsScalar * xx = &x[3];
  int off = 3*var_offset;

  for ( int i = var_offset+1; i < nrows; i++ ){
    int end = diag[i];
    int k   = rowp[i];
    while ( cols[k] < var_offset ) k++;

    const TacsScalar * a = &A[9*k];
    for ( ; k < end; k++ ){
      int j = 3*cols[k] - off;

      xx[0] -= a[0]*x[j] + a[1]*x[j+1] + a[2]*x[j+2];
      xx[1] -= a[3]*x[j] + a[4]*x[j+1] + a[5]*x[j+2];
      xx[2] -= a[6]*x[j] + a[7]*x[j+1] + a[8]*x[j+2];
      a += 9;
    }
    xx += 3;
  }
}


/*!
  Apply a portion of he upper factorization x = U^{-1} x
*/

void BCSRMatApplyPartialUpper3( BCSRMatData * data, TacsScalar * x, 
                                int var_offset ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;
  const int * diag = data->diag;
  const TacsScalar * A = data->A;

  TacsScalar y0, y1, y2;
  TacsScalar * xx = &x[3*(nrows-var_offset-1)];
  int off = 3*var_offset;

  for ( int i = nrows-1; i >= var_offset; i-- ){
    y0 = xx[0];
    y1 = xx[1];
    y2 = xx[2];

    int end = rowp[i+1];
    int k   = diag[i]+1;
    const TacsScalar * a = &A[9*k];
    for ( ; k < end; k++ ){
      int j = 3*cols[k] - off;

      y0 -= a[0]*x[j] + a[1]*x[j+1] + a[2]*x[j+2];
      y1 -= a[3]*x[j] + a[4]*x[j+1] + a[5]*x[j+2];
      y2 -= a[6]*x[j] + a[7]*x[j+1] + a[8]*x[j+2];
      a += 9;
    }

    a = &A[9*diag[i]];
    xx[0] = a[0]*y0 + a[1]*y1 + a[2]*y2;
    xx[1] = a[3]*y0 + a[4]*y1 + a[5]*y2;
    xx[2] = a[6]*y0 + a[7]*y1 + a[8]*y2;
    xx -= 3;
  }
}

/*!
  Function for the approximate Schur preconditioner
*/

void BCSRMatApplyFactorSchur3( BCSRMatData * data, TacsScalar * x, 
                               int var_offset ){
  const int * rowp = data->rowp;
  const int * cols = data->cols;
  const int * diag = data->diag;
  const TacsScalar * A = data->A;

  TacsScalar y0, y1, y2;
  TacsScalar * xx = &x[3*(var_offset-1)];

  for ( int i = var_offset-1; i >= 0; i-- ){
    y0 = xx[0];
    y1 = xx[1];
    y2 = xx[2];

    int end = rowp[i+1];
    int k   = diag[i]+1;
    const TacsScalar * a = &A[9*k];
    for ( ; k < end; k++ ){
      int j = 3*cols[k];

      y0 -= a[0]*x[j] + a[1]*x[j+1] + a[2]*x[j+2];
      y1 -= a[3]*x[j] + a[4]*x[j+1] + a[5]*x[j+2];
      y2 -= a[6]*x[j] + a[7]*x[j+1] + a[8]*x[j+2];
      a += 9;
    }

    a = &A[9*diag[i]];
    xx[0] = a[0]*y0 + a[1]*y1 + a[2]*y2;
    xx[1] = a[3]*y0 + a[4]*y1 + a[5]*y2;
    xx[2] = a[6]*y0 + a[7]*y1 + a[8]*y2;
    xx -= 3;
  }
}

/*!
  Perform a matrix-matrix multiplication
*/

void BCSRMatMatMultAdd3( double alpha, BCSRMatData * Adata, 
                         BCSRMatData * Bdata, BCSRMatData * Cdata ){

  // Retrieve the data required from the matrix
  const int nrows_a = Adata->nrows;
  const int * arowp = Adata->rowp;
  const int * acols = Adata->cols;
  const TacsScalar * A = Adata->A;

  const int * browp = Bdata->rowp;
  const int * bcols = Bdata->cols;
  const TacsScalar * B = Bdata->A;

  // The matrix being written to
  const int * crowp = Cdata->rowp;
  const int * ccols = Cdata->cols;
  TacsScalar * C = Cdata->A;

  if ( alpha == 1.0 ){
    // C_{ik} = A_{ij} B_{jk}
    for ( int i = 0; i < nrows_a; i++ ){
      for ( int jp = arowp[i]; jp < arowp[i+1]; jp++ ){
        int j = acols[jp];
        const TacsScalar * a = &A[9*jp];
        
        int kp     = browp[j];
        int kp_end = browp[j+1];
        const TacsScalar * b = &B[9*kp];
        
        int cp     = crowp[i];
        int cp_end = crowp[i+1];
        TacsScalar * c = &C[9*cp];
        
        for ( ; kp < kp_end; kp++ ){
          while ( ( cp < cp_end ) && ( ccols[cp] < bcols[kp] ) ){ cp++; c+= 9; }
          if ( cp >= cp_end ){ break; }

          if ( bcols[kp] == ccols[cp] ){
            c[0] += a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
            c[3] += a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
            c[6] += a[6]*b[0] + a[7]*b[3] + a[8]*b[6];

            c[1] += a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
            c[4] += a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
            c[7] += a[6]*b[1] + a[7]*b[4] + a[8]*b[7];

            c[2] += a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
            c[5] += a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
            c[8] += a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
          }
          b += 9;
        }
      }
    }
  }
  else if ( alpha == -1.0 ){
    // C_{ik} = A_{ij} B_{jk}
    for ( int i = 0; i < nrows_a; i++ ){
      for ( int jp = arowp[i]; jp < arowp[i+1]; jp++ ){
        int j = acols[jp];
        const TacsScalar * a = &A[9*jp];
        
        int kp     = browp[j];
        int kp_end = browp[j+1];
        const TacsScalar * b = &B[9*kp];
        
        int cp     = crowp[i];
        int cp_end = crowp[i+1];
        TacsScalar * c = &C[9*cp];
        
        for ( ; kp < kp_end; kp++ ){
          while ( ( cp < cp_end ) && ( ccols[cp] < bcols[kp] ) ){ cp++; c += 9; }
          if ( cp >= cp_end ){ break; }

          if ( bcols[kp] == ccols[cp] ){
            c[0] -= a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
            c[3] -= a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
            c[6] -= a[6]*b[0] + a[7]*b[3] + a[8]*b[6];

            c[1] -= a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
            c[4] -= a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
            c[7] -= a[6]*b[1] + a[7]*b[4] + a[8]*b[7];

            c[2] -= a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
            c[5] -= a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
            c[8] -= a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
          }

          b += 9;
        }
      }
    }
  }
  else {

    // C_{ik} = A_{ij} B_{jk}
    for ( int i = 0; i < nrows_a; i++ ){
      for ( int jp = arowp[i]; jp < arowp[i+1]; jp++ ){
        int j = acols[jp];
        const TacsScalar * a = &A[9*jp];
        
        int kp     = browp[j];
        int kp_end = browp[j+1];
        const TacsScalar * b = &B[9*kp];
        
        int cp     = crowp[i];
        int cp_end = crowp[i+1];
        TacsScalar * c = &C[9*cp];
        
        for ( ; kp < kp_end; kp++ ){
          while ( ( cp < cp_end ) && ( ccols[cp] < bcols[kp] ) ){ cp++; c += 9; }
          if ( cp >= cp_end ){ break; }

          if ( bcols[kp] == ccols[cp] ){
            c[0] += alpha*( a[0]*b[0] + a[1]*b[3] + a[2]*b[6] );
            c[3] += alpha*( a[3]*b[0] + a[4]*b[3] + a[5]*b[6] );
            c[6] += alpha*( a[6]*b[0] + a[7]*b[3] + a[8]*b[6] );

            c[1] += alpha*( a[0]*b[1] + a[1]*b[4] + a[2]*b[7] );
            c[4] += alpha*( a[3]*b[1] + a[4]*b[4] + a[5]*b[7] );
            c[7] += alpha*( a[6]*b[1] + a[7]*b[4] + a[8]*b[7] );

            c[2] += alpha*( a[0]*b[2] + a[1]*b[5] + a[2]*b[8] );
            c[5] += alpha*( a[3]*b[2] + a[4]*b[5] + a[5]*b[8] );
            c[8] += alpha*( a[6]*b[2] + a[7]*b[5] + a[8]*b[8] );
          }

          b += 9;
        }
      }
    }
  }
}

/*!
  Apply a given number of steps of SOR to the system A*x = b.
*/
void BCSRMatApplySOR3( BCSRMatData * data, TacsScalar * Adiag,
                       TacsScalar omega, int iters, TacsScalar * b, 
                       TacsScalar * x ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;

  TacsScalar t1, t2, t3;

  for ( int iter = 0; iter < iters; iter++ ){
    for ( int i = 0; i < nrows; i++ ){
      // Copy the right-hand-side to the temporary vector
      // for this row
      t1 = b[3*i];
      t2 = b[3*i+1];
      t3 = b[3*i+2];
        
      // Set the pointer to the beginning of the current
      // row
      TacsScalar * a = &data->A[9*rowp[i]];

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = rowp[i+1];
      for ( int k = rowp[i]; k < end; k++ ){
        int j = cols[k];
        TacsScalar *y = &x[3*j];

        if (i != j){
          t1 -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2];
          t2 -= a[3]*y[0] + a[4]*y[1] + a[5]*y[2];
          t3 -= a[6]*y[0] + a[7]*y[1] + a[8]*y[2];
        }
        
        // Increment the block pointer by bsize^2
        a += 9;
      }

      // Set a pointer to the inverse of the diagonal
      TacsScalar * d = &Adiag[9*i];

      // Compute the first term in the update:
      // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
      x[3*i]   = (1.0 - omega)*x[3*i]   + omega*(d[0]*t1 + d[1]*t2 + d[2]*t3);
      x[3*i+1] = (1.0 - omega)*x[3*i+1] + omega*(d[3]*t1 + d[4]*t2 + d[5]*t3);
      x[3*i+2] = (1.0 - omega)*x[3*i+2] + omega*(d[6]*t1 + d[7]*t2 + d[8]*t3);
    }
  }  
}

/*!
  Apply a given number of steps of symmetric SOR to the system A*x = b.
*/
void BCSRMatApplySSOR3( BCSRMatData * data, TacsScalar * Adiag,
                        TacsScalar omega, int iters, TacsScalar * b, 
                        TacsScalar * x ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;

  TacsScalar t1, t2, t3;

  for ( int iter = 0; iter < iters; iter++ ){
    // Go through the matrix with the forward ordering
    for ( int i = 0; i < nrows; i++ ){
      // Copy the right-hand-side to the temporary vector
      // for this row
      t1 = b[3*i];
      t2 = b[3*i+1];
      t3 = b[3*i+2];
        
      // Set the pointer to the beginning of the current
      // row
      TacsScalar * a = &data->A[9*rowp[i]];

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = rowp[i+1];
      for ( int k = rowp[i]; k < end; k++ ){
        int j = cols[k];
        TacsScalar * y = &x[3*j];

        if (i != j){
          t1 -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2];
          t2 -= a[3]*y[0] + a[4]*y[1] + a[5]*y[2];
          t3 -= a[6]*y[0] + a[7]*y[1] + a[8]*y[2];
        }
        
        // Increment the block pointer by bsize^2
        a += 9;
      }

      // Set a pointer to the inverse of the diagonal
      TacsScalar * d = &Adiag[9*i];

      // Compute the first term in the update:
      // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
      x[3*i]   = (1.0 - omega)*x[3*i]   + omega*(d[0]*t1 + d[1]*t2 + d[2]*t3);
      x[3*i+1] = (1.0 - omega)*x[3*i+1] + omega*(d[3]*t1 + d[4]*t2 + d[5]*t3);
      x[3*i+2] = (1.0 - omega)*x[3*i+2] + omega*(d[6]*t1 + d[7]*t2 + d[8]*t3);
    }

    // Go through the matrix with the reverse orderin
    for ( int i = nrows-1; i >= 0; i-- ){
      // Copy the right-hand-side to the temporary vector
      // for this row
      t1 = b[3*i];
      t2 = b[3*i+1];
      t3 = b[3*i+2];
        
      // Set the pointer to the beginning of the current
      // row
      TacsScalar * a = &data->A[9*rowp[i]];

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = rowp[i+1];
      for ( int k = rowp[i]; k < end; k++ ){
        int j = cols[k];
        TacsScalar * y = &x[3*j];

        if (i != j){
          t1 -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2];
          t2 -= a[3]*y[0] + a[4]*y[1] + a[5]*y[2];
          t3 -= a[6]*y[0] + a[7]*y[1] + a[8]*y[2];
        }
        
        // Increment the block pointer by bsize^2
        a += 9;
      }

      // Set a pointer to the inverse of the diagonal
      TacsScalar * d = &Adiag[9*i];

      // Compute the first term in the update:
      // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
      x[3*i]   = (1.0 - omega)*x[3*i]   + omega*(d[0]*t1 + d[1]*t2 + d[2]*t3);
      x[3*i+1] = (1.0 - omega)*x[3*i+1] + omega*(d[3]*t1 + d[4]*t2 + d[5]*t3);
      x[3*i+2] = (1.0 - omega)*x[3*i+2] + omega*(d[6]*t1 + d[7]*t2 + d[8]*t3);
    }
  }  
}

/*!
  Apply a given number of steps of SOR to the system A*x = b.
*/
void BCSRMatApplySOR3( BCSRMatData *data, TacsScalar *Adiag,
                       const int *pairs, int npairs,
                       TacsScalar omega, int iters, TacsScalar *b, 
                       TacsScalar *x ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;

  for ( int iter = 0; iter < iters; iter++ ){
    const TacsScalar *D = Adiag;
    
    for ( int i = 0, p = 0; i < nrows; i++ ){  
      if (p < npairs && i == pairs[p]){        
        // Copy the right-hand-side to the temporary vector
        // for this row
        TacsScalar t1, t2, t3, t4, t5, t6;
        t1 = b[3*i];
        t2 = b[3*i+1];
        t3 = b[3*i+2];
        t4 = b[3*i+3];
        t5 = b[3*i+4];
        t6 = b[3*i+5];
        
        // Set the pointer to the beginning of the current
        // row
        const TacsScalar *a = &data->A[9*rowp[i]];
        
        // Scan through the row and compute the result:
        // tx <- b_i - A_{ij}*x_{j} for j != i
        int end = rowp[i+1];
        for ( int k = rowp[i]; k < end; k++ ){
          int j = cols[k];
          TacsScalar *y = &x[3*j];

          if (j != i && j != i+1){
            t1 -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2];
            t2 -= a[3]*y[0] + a[4]*y[1] + a[5]*y[2];
            t3 -= a[6]*y[0] + a[7]*y[1] + a[8]*y[2];
          }
          
          // Increment the block pointer by bsize^2
          a += 9;
        }

        end = rowp[i+2];
        for ( int k = rowp[i+1]; k < end; k++ ){
          int j = cols[k];
          TacsScalar *y = &x[3*j];

          if (j != i && j != i+1){
            t4 -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2];
            t5 -= a[3]*y[0] + a[4]*y[1] + a[5]*y[2];
            t6 -= a[6]*y[0] + a[7]*y[1] + a[8]*y[2];
          }
          
          // Increment the block pointer by bsize^2
          a += 9;
        }

        // Compute the first term in the update:
        // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
        x[3*i]   = (1.0 - omega)*x[3*i]   + omega*(D[0]*t1 + D[1]*t2 + D[2]*t3 + D[3]*t4 + D[4]*t5 + D[5]*t6);
        x[3*i+1] = (1.0 - omega)*x[3*i+1] + omega*(D[6]*t1 + D[7]*t2 + D[8]*t3 + D[9]*t4 + D[10]*t5 + D[11]*t6);
        x[3*i+2] = (1.0 - omega)*x[3*i+2] + omega*(D[12]*t1 + D[13]*t2 + D[14]*t3 + D[15]*t4 + D[16]*t5 + D[17]*t6);
        x[3*i+3] = (1.0 - omega)*x[3*i+3] + omega*(D[18]*t1 + D[19]*t2 + D[20]*t3 + D[21]*t4 + D[22]*t5 + D[23]*t6);
        x[3*i+4] = (1.0 - omega)*x[3*i+4] + omega*(D[24]*t1 + D[25]*t2 + D[26]*t3 + D[27]*t4 + D[28]*t5 + D[29]*t6);
        x[3*i+5] = (1.0 - omega)*x[3*i+5] + omega*(D[30]*t1 + D[31]*t2 + D[32]*t3 + D[33]*t4 + D[34]*t5 + D[35]*t6);

        // Increment the pointer
        p++;
        i++;;
        D += 36;
      }
      else {
        // Copy the right-hand-side to the temporary vector
        // for this row
        TacsScalar t1, t2, t3;
        t1 = b[3*i];
        t2 = b[3*i+1];
        t3 = b[3*i+2];
        
        // Set the pointer to the beginning of the current
        // row
        const TacsScalar *a = &data->A[9*rowp[i]];
        
        // Scan through the row and compute the result:
        // tx <- b_i - A_{ij}*x_{j} for j != i
        int end = rowp[i+1];
        for ( int k = rowp[i]; k < end; k++ ){
          int j = cols[k];
          TacsScalar *y = &x[3*j];

          if (i != j){
            t1 -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2];
            t2 -= a[3]*y[0] + a[4]*y[1] + a[5]*y[2];
            t3 -= a[6]*y[0] + a[7]*y[1] + a[8]*y[2];
          }
          
          // Increment the block pointer by bsize^2
          a += 9;
        }
        
        // Compute the first term in the update:
        // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
        x[3*i]   = (1.0 - omega)*x[3*i]   + omega*(D[0]*t1 + D[1]*t2 + D[2]*t3);
        x[3*i+1] = (1.0 - omega)*x[3*i+1] + omega*(D[3]*t1 + D[4]*t2 + D[5]*t3);
        x[3*i+2] = (1.0 - omega)*x[3*i+2] + omega*(D[6]*t1 + D[7]*t2 + D[8]*t3);
        D += 9;
      }
    }
  }
}

/*!
  Apply a given number of steps of symmetric SOR to the system A*x = b.
*/
void BCSRMatApplySSOR3( BCSRMatData *data, TacsScalar *Adiag,
                        const int *pairs, int npairs,
                        TacsScalar omega, int iters, TacsScalar * b, 
                        TacsScalar * x ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;

  for ( int iter = 0; iter < iters; iter++ ){
    // Pointer to the last pair
    const TacsScalar *D = Adiag;
    
    for ( int i = 0, p = 0; i < nrows; i++ ){  
      if (p < npairs && i == pairs[p]){
        // Copy the right-hand-side to the temporary vector
        // for this row
        TacsScalar t1, t2, t3, t4, t5, t6;
        t1 = b[3*i];
        t2 = b[3*i+1];
        t3 = b[3*i+2];
        t4 = b[3*i+3];
        t5 = b[3*i+4];
        t6 = b[3*i+5];
        
        // Set the pointer to the beginning of the current
        // row
        const TacsScalar *a = &data->A[9*rowp[i]];
        
        // Scan through the row and compute the result:
        // tx <- b_i - A_{ij}*x_{j} for j != i
        int end = rowp[i+1];
        for ( int k = rowp[i]; k < end; k++ ){
          int j = cols[k];
          TacsScalar *y = &x[3*j];

          if (j != i && j != i+1){
            t1 -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2];
            t2 -= a[3]*y[0] + a[4]*y[1] + a[5]*y[2];
            t3 -= a[6]*y[0] + a[7]*y[1] + a[8]*y[2];
          }
          
          // Increment the block pointer by bsize^2
          a += 9;
        }

        end = rowp[i+2];
        for ( int k = rowp[i+1]; k < end; k++ ){
          int j = cols[k];
          TacsScalar *y = &x[3*j];

          if (j != i && j != i+1){
            t4 -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2];
            t5 -= a[3]*y[0] + a[4]*y[1] + a[5]*y[2];
            t6 -= a[6]*y[0] + a[7]*y[1] + a[8]*y[2];
          }
          
          // Increment the block pointer by bsize^2
          a += 9;
        }

        // Compute the first term in the update:
        // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
        x[3*i]   = (1.0 - omega)*x[3*i]   + omega*(D[0]*t1 + D[1]*t2 + D[2]*t3 + D[3]*t4 + D[4]*t5 + D[5]*t6);
        x[3*i+1] = (1.0 - omega)*x[3*i+1] + omega*(D[6]*t1 + D[7]*t2 + D[8]*t3 + D[9]*t4 + D[10]*t5 + D[11]*t6);
        x[3*i+2] = (1.0 - omega)*x[3*i+2] + omega*(D[12]*t1 + D[13]*t2 + D[14]*t3 + D[15]*t4 + D[16]*t5 + D[17]*t6);
        x[3*i+3] = (1.0 - omega)*x[3*i+3] + omega*(D[18]*t1 + D[19]*t2 + D[20]*t3 + D[21]*t4 + D[22]*t5 + D[23]*t6);
        x[3*i+4] = (1.0 - omega)*x[3*i+4] + omega*(D[24]*t1 + D[25]*t2 + D[26]*t3 + D[27]*t4 + D[28]*t5 + D[29]*t6);
        x[3*i+5] = (1.0 - omega)*x[3*i+5] + omega*(D[30]*t1 + D[31]*t2 + D[32]*t3 + D[33]*t4 + D[34]*t5 + D[35]*t6);

        // Increment the rows/pair number and the pointer to the
        // diagonal block
        p++;
        i++;
        D += 36;
      }
      else {
        // Copy the right-hand-side to the temporary vector
        // for this row
        TacsScalar t1, t2, t3;
        t1 = b[3*i];
        t2 = b[3*i+1];
        t3 = b[3*i+2];
        
        // Set the pointer to the beginning of the current
        // row
        const TacsScalar *a = &data->A[9*rowp[i]];
        
        // Scan through the row and compute the result:
        // tx <- b_i - A_{ij}*x_{j} for j != i
        int end = rowp[i+1];
        for ( int k = rowp[i]; k < end; k++ ){
          int j = cols[k];
          TacsScalar *y = &x[3*j];

          if (i != j){
            t1 -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2];
            t2 -= a[3]*y[0] + a[4]*y[1] + a[5]*y[2];
            t3 -= a[6]*y[0] + a[7]*y[1] + a[8]*y[2];
          }
          
          // Increment the block pointer by bsize^2
          a += 9;
        }
        
        // Compute the first term in the update:
        // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
        x[3*i]   = (1.0 - omega)*x[3*i]   + omega*(D[0]*t1 + D[1]*t2 + D[2]*t3);
        x[3*i+1] = (1.0 - omega)*x[3*i+1] + omega*(D[3]*t1 + D[4]*t2 + D[5]*t3);
        x[3*i+2] = (1.0 - omega)*x[3*i+2] + omega*(D[6]*t1 + D[7]*t2 + D[8]*t3);
        D += 9;
      }
    }

    for ( int i = nrows-1, p = npairs-1; i >= 0; i-- ){
      if (p >= 0 && i-1 == pairs[p]){
        // Copy the right-hand-side to the temporary vector
        // for this row
        TacsScalar t1, t2, t3, t4, t5, t6;
        t1 = b[3*i-3];
        t2 = b[3*i-2];
        t3 = b[3*i-1];
        t4 = b[3*i];
        t5 = b[3*i+1];
        t6 = b[3*i+2];
        
        // Set the pointer to the beginning of the current
        // row
        const TacsScalar *a = &data->A[9*rowp[i-1]];
        
        // Scan through the row and compute the result:
        // tx <- b_i - A_{ij}*x_{j} for j != i
        int end = rowp[i];
        for ( int k = rowp[i-1]; k < end; k++ ){
          int j = cols[k];
          TacsScalar *y = &x[3*j];

          if (j != i-1 && j != i){
            t1 -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2];
            t2 -= a[3]*y[0] + a[4]*y[1] + a[5]*y[2];
            t3 -= a[6]*y[0] + a[7]*y[1] + a[8]*y[2];
          }
          
          // Increment the block pointer by bsize^2
          a += 9;
        }

        end = rowp[i+1];
        for ( int k = rowp[i]; k < end; k++ ){
          int j = cols[k];
          TacsScalar *y = &x[3*j];

          if (j != i-1 && j != i){
            t4 -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2];
            t5 -= a[3]*y[0] + a[4]*y[1] + a[5]*y[2];
            t6 -= a[6]*y[0] + a[7]*y[1] + a[8]*y[2];
          }
          
          // Increment the block pointer by bsize^2
          a += 9;
        }

        // Compute the first term in the update:
        // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
        D -= 36;
        x[3*i-3] = (1.0 - omega)*x[3*i-3] + omega*(D[0]*t1 + D[1]*t2 + D[2]*t3 + D[3]*t4 + D[4]*t5 + D[5]*t6);
        x[3*i-2] = (1.0 - omega)*x[3*i-2] + omega*(D[6]*t1 + D[7]*t2 + D[8]*t3 + D[9]*t4 + D[10]*t5 + D[11]*t6);
        x[3*i-1] = (1.0 - omega)*x[3*i-1] + omega*(D[12]*t1 + D[13]*t2 + D[14]*t3 + D[15]*t4 + D[16]*t5 + D[17]*t6);
        x[3*i]   = (1.0 - omega)*x[3*i]   + omega*(D[18]*t1 + D[19]*t2 + D[20]*t3 + D[21]*t4 + D[22]*t5 + D[23]*t6);
        x[3*i+1] = (1.0 - omega)*x[3*i+1] + omega*(D[24]*t1 + D[25]*t2 + D[26]*t3 + D[27]*t4 + D[28]*t5 + D[29]*t6);
        x[3*i+2] = (1.0 - omega)*x[3*i+2] + omega*(D[30]*t1 + D[31]*t2 + D[32]*t3 + D[33]*t4 + D[34]*t5 + D[35]*t6);

        // Increment the rows/pair number and the pointer to the
        // diagonal block
        p--;
        i--;
      }
      else {
        // Copy the right-hand-side to the temporary vector
        // for this row
        TacsScalar t1, t2, t3;
        t1 = b[3*i];
        t2 = b[3*i+1];
        t3 = b[3*i+2];
        
        // Set the pointer to the beginning of the current
        // row
        const TacsScalar *a = &data->A[9*rowp[i]];
        
        // Scan through the row and compute the result:
        // tx <- b_i - A_{ij}*x_{j} for j != i
        int end = rowp[i+1];
        for ( int k = rowp[i]; k < end; k++ ){
          int j = cols[k];
          TacsScalar *y = &x[3*j];

          if (i != j){
            t1 -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2];
            t2 -= a[3]*y[0] + a[4]*y[1] + a[5]*y[2];
            t3 -= a[6]*y[0] + a[7]*y[1] + a[8]*y[2];
          }
          
          // Increment the block pointer by bsize^2
          a += 9;
        }
        
        // Compute the first term in the update:
        // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
        D -= 9;
        x[3*i]   = (1.0 - omega)*x[3*i]   + omega*(D[0]*t1 + D[1]*t2 + D[2]*t3);
        x[3*i+1] = (1.0 - omega)*x[3*i+1] + omega*(D[3]*t1 + D[4]*t2 + D[5]*t3);
        x[3*i+2] = (1.0 - omega)*x[3*i+2] + omega*(D[6]*t1 + D[7]*t2 + D[8]*t3);
      }
    }
  }
}
