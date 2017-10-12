#include "BCSRMatImpl.h"

/*
  The following file contains the specific implementation for block size = 4.

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

/*!
  Compute the matrix-vector product: y = A * x
*/
void BCSRMatVecMult4( BCSRMatData * data,
                      TacsScalar * x, TacsScalar * y ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;
  const TacsScalar * a = data->A;

  for ( int i = 0; i < nrows; i++ ){
    y[0] = 0.0;
    y[1] = 0.0;
    y[2] = 0.0;
    y[3] = 0.0;

    int end = rowp[i+1];
    for ( int k = rowp[i]; k < end; k++ ){
      int j = 4*cols[k];
        
      y[0] += a[0 ]*x[j] + a[1 ]*x[j+1] + a[2 ]*x[j+2] + a[3 ]*x[j+3];
      y[1] += a[4 ]*x[j] + a[5 ]*x[j+1] + a[6 ]*x[j+2] + a[7 ]*x[j+3];
      y[2] += a[8 ]*x[j] + a[9 ]*x[j+1] + a[10]*x[j+2] + a[11]*x[j+3];
      y[3] += a[12]*x[j] + a[13]*x[j+1] + a[14]*x[j+2] + a[15]*x[j+3];
      a += 16;
    }
    y += 4;
  }
}

/*!
  Compute the matrix vector product plus addition: z = A * x + y
*/
void BCSRMatVecMultAdd4( BCSRMatData * data,
                         TacsScalar * x, TacsScalar * y, 
                         TacsScalar * z ){

  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;
  const TacsScalar * a = data->A;

  for ( int i = 0; i < nrows; i++ ){
    y[0] = z[0];
    y[1] = z[1];
    y[2] = z[2];
    y[3] = z[3];

    int end = rowp[i+1];
    for ( int k = rowp[i]; k < end; k++ ){
      int j = 4*cols[k];
      y[0] += a[0 ]*x[j] + a[1 ]*x[j+1] + a[2 ]*x[j+2] + a[3 ]*x[j+3];
      y[1] += a[4 ]*x[j] + a[5 ]*x[j+1] + a[6 ]*x[j+2] + a[7 ]*x[j+3];
      y[2] += a[8 ]*x[j] + a[9 ]*x[j+1] + a[10]*x[j+2] + a[11]*x[j+3];
      y[3] += a[12]*x[j] + a[13]*x[j+1] + a[14]*x[j+2] + a[15]*x[j+3];

      a += 16;
    }
    y += 4;
    z += 4;
  }
}

/*!
  Apply the lower factorization y = L^{-1} x
*/
void BCSRMatApplyLower4( BCSRMatData * data,
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
    z[3] = x[3];

    int end = diag[i];
    int k   = rowp[i];
    const TacsScalar * a = &A[16*k];
    for ( ; k < end; k++ ){
      int j = 4*cols[k];

      z[0] -= a[0 ]*y[j] + a[1 ]*y[j+1] + a[2 ]*y[j+2] + a[3 ]*y[j+3];
      z[1] -= a[4 ]*y[j] + a[5 ]*y[j+1] + a[6 ]*y[j+2] + a[7 ]*y[j+3];
      z[2] -= a[8 ]*y[j] + a[9 ]*y[j+1] + a[10]*y[j+2] + a[11]*y[j+3];
      z[3] -= a[12]*y[j] + a[13]*y[j+1] + a[14]*y[j+2] + a[15]*y[j+3];
      a += 16;
    }
    z += 4;
    x += 4;
  }
}

/*!
  Apply the upper factorization y = U^{-1} x
*/
void BCSRMatApplyUpper4( BCSRMatData * data,
                         TacsScalar * x, TacsScalar * y ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;
  const int * diag = data->diag;
  const TacsScalar * A = data->A;

  TacsScalar y0, y1, y2, y3;

  x = &x[4*(nrows-1)];
  for ( int i = nrows-1; i >= 0; i-- ){
    y0 = x[0];
    y1 = x[1];
    y2 = x[2];
    y3 = x[3];
    
    int end = rowp[i+1];
    int k   = diag[i]+1;
    const TacsScalar * a = &A[16*k];
    for ( ; k < end; k++ ){
      int j = 4*cols[k];
      y0 -= a[0 ]*y[j] + a[1 ]*y[j+1] + a[2 ]*y[j+2] + a[3 ]*y[j+3];
      y1 -= a[4 ]*y[j] + a[5 ]*y[j+1] + a[6 ]*y[j+2] + a[7 ]*y[j+3];
      y2 -= a[8 ]*y[j] + a[9 ]*y[j+1] + a[10]*y[j+2] + a[11]*y[j+3];
      y3 -= a[12]*y[j] + a[13]*y[j+1] + a[14]*y[j+2] + a[15]*y[j+3];
      a += 16;
    }

    int bi = 4*i;
    a = &A[16*diag[i]];
    
    y[bi  ] = a[0 ]*y0 + a[1 ]*y1 + a[2 ]*y2 + a[3 ]*y3;
    y[bi+1] = a[4 ]*y0 + a[5 ]*y1 + a[6 ]*y2 + a[7 ]*y3;
    y[bi+2] = a[8 ]*y0 + a[9 ]*y1 + a[10]*y2 + a[11]*y3;
    y[bi+3] = a[12]*y0 + a[13]*y1 + a[14]*y2 + a[15]*y3;

    x -= 4;
  }
}

/*!
  Apply a portion of the lower factorization x = L^{-1} x
*/

void BCSRMatApplyPartialLower4( BCSRMatData * data, 
                                TacsScalar * x, 
                                int var_offset ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;
  const int * diag = data->diag;
  const TacsScalar * A = data->A;

  TacsScalar * xx = &x[4];
  int off = 4*var_offset;

  for ( int i = var_offset+1; i < nrows; i++ ){
    int end = diag[i];
    int k   = rowp[i];
    while ( cols[k] < var_offset ) k++;

    const TacsScalar * a = &A[16];
    for ( ; k < end; k++ ){
      int j = 4*cols[k] - off;

      xx[0] -= a[0 ]*x[j] + a[1 ]*x[j+1] + a[2 ]*x[j+2] 
        + a[3 ]*x[j+3];
      xx[1] -= a[4 ]*x[j] + a[5 ]*x[j+1] + a[6 ]*x[j+2] 
        + a[7 ]*x[j+3];
      xx[2] -= a[8 ]*x[j] + a[9 ]*x[j+1] + a[10]*x[j+2] 
        + a[11]*x[j+3];
      xx[3] -= a[12]*x[j] + a[13]*x[j+1] + a[14]*x[j+2] 
        + a[15]*x[j+3];

      a += 16;
    }
    xx += 4;
  }
}


/*!
  Apply a portion of he upper factorization x = U^{-1} x
*/

void BCSRMatApplyPartialUpper4( BCSRMatData * data, 
                                TacsScalar * x, 
                                int var_offset ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;
  const int * diag = data->diag;
  const TacsScalar * A = data->A;

  TacsScalar y0, y1, y2, y3;
  TacsScalar * xx = &x[4*(nrows-var_offset-1)];
  int off = 4*var_offset;

  for ( int i = nrows-1; i >= var_offset; i-- ){
    y0 = xx[0];
    y1 = xx[1];
    y2 = xx[2];
    y3 = xx[3];

    int end = rowp[i+1];
    int k   = diag[i]+1;
    const TacsScalar * a = &A[16*k];
    for ( ; k < end; k++ ){
      int j = 4*cols[k] - off;
      y0 -= a[0 ]*x[j] + a[1 ]*x[j+1] + a[2 ]*x[j+2] + a[3 ]*x[j+3];
      y1 -= a[4 ]*x[j] + a[5 ]*x[j+1] + a[6 ]*x[j+2] + a[7 ]*x[j+3];
      y2 -= a[8 ]*x[j] + a[9 ]*x[j+1] + a[10]*x[j+2] + a[11]*x[j+3];
      y3 -= a[12]*x[j] + a[13]*x[j+1] + a[14]*x[j+2] + a[15]*x[j+3];
      a += 16;
    }

    a = &A[16*diag[i]];
    xx[0] = a[0 ]*y0 + a[1 ]*y1 + a[2 ]*y2 + a[3 ]*y3;
    xx[1] = a[4 ]*y0 + a[5 ]*y1 + a[6 ]*y2 + a[7 ]*y3;
    xx[2] = a[8 ]*y0 + a[9 ]*y1 + a[10]*y2 + a[11]*y3;
    xx[3] = a[12]*y0 + a[13]*y1 + a[14]*y2 + a[15]*y3;
   
    xx -= 4;
  }
}

/*!
  Function for the approximate Schur preconditioner
*/

void BCSRMatApplyFactorSchur4( BCSRMatData * data, TacsScalar * x, 
                               int var_offset ){
  const int * rowp = data->rowp;
  const int * cols = data->cols;
  const int * diag = data->diag;
  const TacsScalar * A = data->A;

  TacsScalar y0, y1, y2, y3;
  TacsScalar * xx = &x[4*(var_offset-1)];

  for ( int i = var_offset-1; i >= 0; i-- ){
    y0 = xx[0];
    y1 = xx[1];
    y2 = xx[2];
    y3 = xx[3];

    int end = rowp[i+1];
    int k   = diag[i]+1;
    const TacsScalar * a = &A[16*k];
    for ( ; k < end; k++ ){
      int j = 4*cols[k];
      y0 -= a[0 ]*x[j] + a[1 ]*x[j+1] + a[2 ]*x[j+2] + a[3 ]*x[j+3];
      
      y0 -= a[4 ]*x[j] + a[5 ]*x[j+1] + a[6 ]*x[j+2] + a[7 ]*x[j+3];
      y0 -= a[8 ]*x[j] + a[9 ]*x[j+1] + a[10]*x[j+2] + a[11]*x[j+3];
      y0 -= a[12]*x[j] + a[13]*x[j+1] + a[14]*x[j+2] + a[15]*x[j+3];

      a += 16;
    }

    a = &A[16*diag[i]];

    xx[0] = a[0 ]*y0 + a[1 ]*y1 + a[2 ]*y2 + a[3 ]*y3;
    xx[1]= a[4 ]*y0 + a[5 ]*y1 + a[6 ]*y2 + a[7 ]*y3;
    xx[2] = a[8 ]*y0 + a[9 ]*y1 + a[10]*y2 + a[11]*y3;
    xx[3] = a[12]*y0 + a[13]*y1 + a[14]*y2 + a[15]*y3;
  
    xx -= 4;
  }
}

/*!
  Perform a matrix-matrix multiplication
*/

void BCSRMatMatMultAdd4( double alpha, BCSRMatData * Adata, 
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
        const TacsScalar * a = &A[16*jp];
        
        int kp     = browp[j];
        int kp_end = browp[j+1];
        const TacsScalar * b = &B[16*kp];
        
        int cp     = crowp[i];
        int cp_end = crowp[i+1];
        TacsScalar * c = &C[16*cp];
        
        for ( ; kp < kp_end; kp++ ){
          while ( ( cp < cp_end ) && ( ccols[cp] < bcols[kp] ) ){ cp++; c+= 16; }
          if ( cp >= cp_end ){ break; }

          if ( bcols[kp] == ccols[cp] ){
            TacsScalar b0, b1, b2, b3;
            b0 = b[0 ]; b1 = b[4 ]; b2 = b[8 ]; b3 = b[12];
            c[0 ] += a[0 ]*b0 + a[1 ]*b1 + a[2 ]*b2 + a[3 ]*b3;
            c[4 ] += a[4 ]*b0 + a[5 ]*b1 + a[6 ]*b2 + a[7 ]*b3;
            c[8 ] += a[8 ]*b0 + a[9 ]*b1 + a[10]*b2 + a[11]*b3;
            c[12] += a[12]*b0 + a[13]*b1 + a[14]*b2 + a[15]*b3;
            
            b0 = b[1 ]; b1 = b[5 ]; b2 = b[9 ]; b3 = b[13];
            c[1 ] += a[0 ]*b0 + a[1 ]*b1 + a[2 ]*b2 + a[3 ]*b3;
            c[5 ] += a[4 ]*b0 + a[5 ]*b1 + a[6 ]*b2 + a[7 ]*b3;
            c[9 ] += a[8 ]*b0 + a[9 ]*b1 + a[10]*b2 + a[11]*b3;
            c[13] += a[12]*b0 + a[13]*b1 + a[14]*b2 + a[15]*b3;

            b0 = b[2 ]; b1 = b[6 ]; b2 = b[10]; b3 = b[14];
            c[2 ] += a[0 ]*b0 + a[1 ]*b1 + a[2 ]*b2 + a[3 ]*b3;
            c[6 ] += a[4 ]*b0 + a[5 ]*b1 + a[6 ]*b2 + a[7 ]*b3;
            c[10] += a[8 ]*b0 + a[9 ]*b1 + a[10]*b2 + a[11]*b3;
            c[14] += a[12]*b0 + a[13]*b1 + a[14]*b2 + a[15]*b3;

            b0 = b[3 ]; b1 = b[7 ]; b2 = b[11]; b3 = b[15];
            c[3 ] += a[0 ]*b0 + a[1 ]*b1 + a[2 ]*b2 + a[3 ]*b3;
            c[7 ] += a[4 ]*b0 + a[5 ]*b1 + a[6 ]*b2 + a[7 ]*b3;
            c[11] += a[8 ]*b0 + a[9 ]*b1 + a[10]*b2 + a[11]*b3;
            c[15] += a[12]*b0 + a[13]*b1 + a[14]*b2 + a[15]*b3;

          }
          b += 16;
        }
      }
    }
  }
  else if ( alpha == -1.0 ){
    // C_{ik} = A_{ij} B_{jk}
    for ( int i = 0; i < nrows_a; i++ ){
      for ( int jp = arowp[i]; jp < arowp[i+1]; jp++ ){
        int j = acols[jp];
        const TacsScalar * a = &A[16*jp];
        
        int kp     = browp[j];
        int kp_end = browp[j+1];
        const TacsScalar * b = &B[16*kp];
        
        int cp     = crowp[i];
        int cp_end = crowp[i+1];
        TacsScalar * c = &C[16*cp];
        
        for ( ; kp < kp_end; kp++ ){
          while ( ( cp < cp_end ) && ( ccols[cp] < bcols[kp] ) ){ cp++; c += 16; }
          if ( cp >= cp_end ){ break; }

          if ( bcols[kp] == ccols[cp] ){
            TacsScalar b0, b1, b2, b3;
            b0 = b[0 ]; b1 = b[4 ]; b2 = b[8 ]; b3 = b[12];
            c[0 ] -= a[0 ]*b0 + a[1 ]*b1 + a[2 ]*b2 + a[3 ]*b3;
            c[4 ] -= a[4 ]*b0 + a[5 ]*b1 + a[6 ]*b2 + a[7 ]*b3;
            c[8 ] -= a[8 ]*b0 + a[9 ]*b1 + a[10]*b2 + a[11]*b3;
            c[12] -= a[12]*b0 + a[13]*b1 + a[14]*b2 + a[15]*b3;
            
            b0 = b[1 ]; b1 = b[5 ]; b2 = b[9 ]; b3 = b[13];
            c[1 ] -= a[0 ]*b0 + a[1 ]*b1 + a[2 ]*b2 + a[3 ]*b3;
            c[5 ] -= a[4 ]*b0 + a[5 ]*b1 + a[6 ]*b2 + a[7 ]*b3;
            c[9 ] -= a[8 ]*b0 + a[9 ]*b1 + a[10]*b2 + a[11]*b3;
            c[13] -= a[12]*b0 + a[13]*b1 + a[14]*b2 + a[15]*b3;

            b0 = b[2 ]; b1 = b[6 ]; b2 = b[10]; b3 = b[14];
            c[2 ] -= a[0 ]*b0 + a[1 ]*b1 + a[2 ]*b2 + a[3 ]*b3;
            c[6 ] -= a[4 ]*b0 + a[5 ]*b1 + a[6 ]*b2 + a[7 ]*b3;
            c[10] -= a[8 ]*b0 + a[9 ]*b1 + a[10]*b2 + a[11]*b3;
            c[14] -= a[12]*b0 + a[13]*b1 + a[14]*b2 + a[15]*b3;

            b0 = b[3 ]; b1 = b[7 ]; b2 = b[11]; b3 = b[15];
            c[3 ] -= a[0 ]*b0 + a[1 ]*b1 + a[2 ]*b2 + a[3 ]*b3;
            c[7 ] -= a[4 ]*b0 + a[5 ]*b1 + a[6 ]*b2 + a[7 ]*b3;
            c[11] -= a[8 ]*b0 + a[9 ]*b1 + a[10]*b2 + a[11]*b3;
            c[15] -= a[12]*b0 + a[13]*b1 + a[14]*b2 + a[15]*b3;
          }

          b += 16;
        }
      }
    }
  }
  else {

    // C_{ik} = A_{ij} B_{jk}
    for ( int i = 0; i < nrows_a; i++ ){
      for ( int jp = arowp[i]; jp < arowp[i+1]; jp++ ){
        int j = acols[jp];
        const TacsScalar * a = &A[16*jp];
        
        int kp     = browp[j];
        int kp_end = browp[j+1];
        const TacsScalar * b = &B[16*kp];
        
        int cp     = crowp[i];
        int cp_end = crowp[i+1];
        TacsScalar * c = &C[16*cp];
        
        for ( ; kp < kp_end; kp++ ){
          while ( ( cp < cp_end ) && ( ccols[cp] < bcols[kp] ) ){ cp++; c += 16; }
          if ( cp >= cp_end ){ break; }

          if ( bcols[kp] == ccols[cp] ){
            TacsScalar b0, b1, b2, b3;
            b0 = b[0 ]; b1 = b[4 ]; b2 = b[8 ]; b3 = b[12];
            c[0 ] += alpha*(a[0 ]*b0 + a[1 ]*b1 + a[2 ]*b2 + a[3 ]*b3);
            c[4 ] += alpha*(a[4 ]*b0 + a[5 ]*b1 + a[6 ]*b2 + a[7 ]*b3);
            c[8 ] += alpha*(a[8 ]*b0 + a[9 ]*b1 + a[10]*b2 + a[11]*b3);
            c[12] += alpha*(a[12]*b0 + a[13]*b1 + a[14]*b2 + a[15]*b3);
            
            c[1 ] += alpha*(a[0 ]*b0 + a[1 ]*b1 + a[2 ]*b2 + a[3 ]*b3);
            c[5 ] += alpha*(a[4 ]*b0 + a[5 ]*b1 + a[6 ]*b2 + a[7 ]*b3);
            c[9 ] += alpha*(a[8 ]*b0 + a[9 ]*b1 + a[10]*b2 + a[11]*b3);
            c[13] += alpha*(a[12]*b0 + a[13]*b1 + a[14]*b2 + a[15]*b3);

            c[2 ] += alpha*(a[0 ]*b0 + a[1 ]*b1 + a[2 ]*b2 + a[3 ]*b3);
            c[6 ] += alpha*(a[4 ]*b0 + a[5 ]*b1 + a[6 ]*b2 + a[7 ]*b3);
            c[10] += alpha*(a[8 ]*b0 + a[9 ]*b1 + a[10]*b2 + a[11]*b3);
            c[14] += alpha*(a[12]*b0 + a[13]*b1 + a[14]*b2 + a[15]*b3);

            c[3 ] += alpha*(a[0 ]*b0 + a[1 ]*b1 + a[2 ]*b2 + a[3 ]*b3);
            c[7 ] += alpha*(a[4 ]*b0 + a[5 ]*b1 + a[6 ]*b2 + a[7 ]*b3);
            c[11] += alpha*(a[8 ]*b0 + a[9 ]*b1 + a[10]*b2 + a[11]*b3);
            c[15] += alpha*(a[12]*b0 + a[13]*b1 + a[14]*b2 + a[15]*b3);

          }

          b += 16;
        }
      }
    }
  }
}

/*!
  Apply a given number of steps of SOR to the system A*x = b.
*/
void BCSRMatApplySOR4( BCSRMatData * data, TacsScalar * Adiag,
                       TacsScalar omega, int iters, TacsScalar * b, 
                       TacsScalar * x ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;

  TacsScalar t1, t2, t3, t4;

  for ( int iter = 0; iter < iters; iter++ ){
    for ( int i = 0; i < nrows; i++ ){
      // Copy the right-hand-side to the temporary vector
      // for this row
      t1 = b[4*i];
      t2 = b[4*i+1];
      t3 = b[4*i+2];
      t4 = b[4*i+3];
        
      // Set the pointer to the beginning of the current
      // row
      TacsScalar * a = &data->A[16*rowp[i]];

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = rowp[i+1];
      for ( int k = rowp[i]; k < end; k++ ){
        int j = cols[k];
        TacsScalar *y = &x[4*j];

        if (i != j){
          t1 -= a[0 ]*y[0] + a[1 ]*y[1] + a[2 ]*y[2] + a[3 ]*y[3];
          t2 -= a[4 ]*y[0] + a[5 ]*y[1] + a[6 ]*y[2] + a[7 ]*y[3];
          t3 -= a[8 ]*y[0] + a[9 ]*y[1] + a[10]*y[2] + a[11]*y[3];
          t4 -= a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3];
        }
        
        // Increment the block pointer by bsize^2
        a += 16;
      }

      // Set a pointer to the inverse of the diagonal
      TacsScalar * d = &Adiag[16*i];

      // Compute the first term in the update:
      // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
      
      x[4*i]   = (1.0 - omega)*x[4*i]   + omega*(d[0 ]*t1 + d[1 ]*t2 + d[2 ]*t3 + d[3 ]*t4);
      x[4*i+1] = (1.0 - omega)*x[4*i+1] + omega*(d[4 ]*t1 + d[5 ]*t2 + d[6 ]*t3 + d[7 ]*t4);
      x[4*i+2] = (1.0 - omega)*x[4*i+2] + omega*(d[8 ]*t1 + d[9 ]*t2 + d[10]*t3 + d[11]*t4);
      x[4*i+3] = (1.0 - omega)*x[4*i+3] + omega*(d[12]*t1 + d[13]*t2 + d[14]*t3 + d[15]*t4);
   
    }
  }  
}

/*!
  Apply a given number of steps of symmetric SOR to the system A*x = b.
*/
void BCSRMatApplySSOR4( BCSRMatData * data, TacsScalar * Adiag,
                        TacsScalar omega, int iters, TacsScalar * b,
                        TacsScalar * x ){
  const int nrows = data->nrows;
  const int * rowp = data->rowp;
  const int * cols = data->cols;

  TacsScalar t1, t2, t3, t4;

  for ( int iter = 0; iter < iters; iter++ ){
    // Go through the matrix with the forward ordering
    for ( int i = 0; i < nrows; i++ ){
      // Copy the right-hand-side to the temporary vector
      // for this row
      t1 = b[4*i];
      t2 = b[4*i+1];
      t3 = b[4*i+2];
      t4 = b[4*i+3];
     
      // Set the pointer to the beginning of the current
      // row
      TacsScalar * a = &data->A[16*rowp[i]];

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = rowp[i+1];
      for ( int k = rowp[i]; k < end; k++ ){
        int j = cols[k];
        TacsScalar * y = &x[4*j];

        if (i != j){
          t1 -= a[0 ]*y[0] + a[1 ]*y[1] + a[2 ]*y[2] + a[3 ]*y[3];
          t2 -= a[4 ]*y[0] + a[5 ]*y[1] + a[6 ]*y[2] + a[7 ]*y[3];
          t3 -= a[8 ]*y[0] + a[9 ]*y[1] + a[10]*y[2] + a[11]*y[3];
          t4 -= a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3];

        }
        
        // Increment the block pointer by bsize^2
        a += 16;
      }

      // Set a pointer to the inverse of the diagonal
      TacsScalar * d = &Adiag[16*i];

      // Compute the first term in the update:
      // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
      x[4*i]   = (1.0 - omega)*x[4*i]   + omega*(d[0 ]*t1 + d[1 ]*t2 + d[2 ]*t3 + d[3 ]*t4);
      x[4*i+1] = (1.0 - omega)*x[4*i+1] + omega*(d[4 ]*t1 + d[5 ]*t2 + d[6 ]*t3 + d[7 ]*t4);
      x[4*i+2] = (1.0 - omega)*x[4*i+2] + omega*(d[8 ]*t1 + d[9 ]*t2 + d[10]*t3 + d[11]*t4);
      x[4*i+3] = (1.0 - omega)*x[4*i+3] + omega*(d[12]*t1 + d[13]*t2 + d[14]*t3 + d[15]*t4);

    }

    // Go through the matrix with the reverse orderin
    for ( int i = nrows-1; i >= 0; i-- ){
      // Copy the right-hand-side to the temporary vector
      // for this row
      t1 = b[4*i];
      t2 = b[4*i+1];
      t3 = b[4*i+2];
      t4 = b[4*i+3];
      // Set the pointer to the beginning of the current
      // row
      TacsScalar * a = &data->A[16*rowp[i]];

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = rowp[i+1];
      for ( int k = rowp[i]; k < end; k++ ){
        int j = cols[k];
        TacsScalar * y = &x[4*j];

        if (i != j){
          t1 -= a[0 ]*y[0] + a[1 ]*y[1] + a[2 ]*y[2] + a[3 ]*y[3];
          t2 -= a[4 ]*y[0] + a[5 ]*y[1] + a[6 ]*y[2] + a[7 ]*y[3];
          t3 -= a[8 ]*y[0] + a[9 ]*y[1] + a[10]*y[2] + a[11]*y[3];
          t4 -= a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3];
        }
        
        // Increment the block pointer by bsize^2
        a += 16;
      }

      // Set a pointer to the inverse of the diagonal
      TacsScalar * d = &Adiag[16*i];

      // Compute the first term in the update:
      // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
      x[4*i]   = (1.0 - omega)*x[4*i]   + omega*(d[0 ]*t1 + d[1 ]*t2 + d[2 ]*t3 + d[3 ]*t4);
      x[4*i+1] = (1.0 - omega)*x[4*i+1] + omega*(d[4 ]*t1 + d[5 ]*t2 + d[6 ]*t3 + d[7 ]*t4);
      x[4*i+2] = (1.0 - omega)*x[4*i+2] + omega*(d[8 ]*t1 + d[9 ]*t2 + d[10]*t3 + d[11]*t4);
      x[4*i+3] = (1.0 - omega)*x[4*i+3] + omega*(d[12]*t1 + d[13]*t2 + d[14]*t3 + d[15]*t4);
     
    }
  }
}

/*!
/*   Apply a given number of steps of SOR to the system A*x = b. */
/* *\/ */
/* void BCSRMatApplySOR4( BCSRMatData *data, TacsScalar *Adiag, */
/*                        const int *pairs, int npairs, */
/*                        TacsScalar omega, int iters, TacsScalar *b,  */
/*                        TacsScalar *x ){ */
/*   const int nrows = data->nrows; */
/*   const int * rowp = data->rowp; */
/*   const int * cols = data->cols; */

/*   for ( int iter = 0; iter < iters; iter++ ){ */
/*     const TacsScalar *D = Adiag; */
    
/*     for ( int i = 0, p = 0; i < nrows; i++ ){   */
/*       if (p < npairs && i == pairs[p]){         */
/*         // Copy the right-hand-side to the temporary vector */
/*         // for this row */
/*         TacsScalar t1, t2, t3, t4, t5, t6, t7, t8; */
/*         t1 = b[4*i]; */
/*         t2 = b[4*i+1]; */
/*         t3 = b[4*i+2]; */
/*         t4 = b[4*i+3]; */
/*         t5 = b[4*i+4]; */
/*         t6 = b[4*i+5]; */
/*         t7 = b[4*i+6]; */
/*         t8 = b[4*i+7]; */
        
/*         // Set the pointer to the beginning of the current */
/*         // row */
/*         const TacsScalar *a = &data->A[16*rowp[i]]; */
        
/*         // Scan through the row and compute the result: */
/*         // tx <- b_i - A_{ij}*x_{j} for j != i */
/*         int end = rowp[i+1]; */
/*         for ( int k = rowp[i]; k < end; k++ ){ */
/*           int j = cols[k]; */
/*           TacsScalar *y = &x[4*j]; */

/*           if (j != i && j != i+1){ */
/*             t1 -= a[0 ]*y[0] + a[1 ]*y[1] + a[2 ]*y[2] + a[3 ]*y[3]; */
/*             t2 -= a[4 ]*y[0] + a[5 ]*y[1] + a[6 ]*y[2] + a[7 ]*y[3]; */
/*             t3 -= a[8 ]*y[0] + a[9 ]*y[1] + a[10]*y[2] + a[11]*y[3]; */
/*             t4 -= a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3]; */
/*           } */
          
/*           // Increment the block pointer by bsize^2 */
/*           a += 16; */
/*         } */

/*         end = rowp[i+3]; */
/*         for ( int k = rowp[i+1]; k < end; k++ ){ */
/*           int j = cols[k]; */
/*           TacsScalar *y = &x[4*j]; */

/*           if (j != i && j != i+1){ */
/*              t5 -= a[0 ]*y[0] + a[1 ]*y[1] + a[2 ]*y[2] + a[3 ]*y[3]; */
/*              t6 -= a[4 ]*y[0] + a[5 ]*y[1] + a[6 ]*y[2] + a[7 ]*y[3]; */
/*              t7 -= a[8 ]*y[0] + a[9 ]*y[1] + a[10]*y[2] + a[11]*y[3]; */
/*              t8 -= a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3];           */
/*           } */
          
/*           // Increment the block pointer by bsize^2 */
/*           a += 16; */
/*         } */

/*         // Compute the first term in the update: */
/*         // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx */
/*         x[4*i  ]  = (1.0 - omega)*x[4*i  ] + omega*(D[0]*t1 + D[1]*t2 + D[2]*t3 + D[3]*t4 + D[4]*t5 + D[5]*t6 + D[6]*t7 + D[7]*t8); */
/*         x[4*i+1]  = (1.0 - omega)*x[4*i+1] + omega*(D[8]*t1 + D[9]*t2 + D[10]*t3 + D[11]*t4 + D[12]*t5 + D[13]*t6 + D[14]*t7 + D[15]*t8); */
/*         x[4*i+2]  = (1.0 - omega)*x[4*i+2] + omega*(D[16]*t1 + D[17]*t2 + D[18]*t3 + D[19]*t4 + D[20]*t5 + D[21]*t6 + D[22]*t7 + D[23]*t8); */
/*         x[4*i+3]  = (1.0 - omega)*x[4*i+3] + omega*(D[24]*t1 + D[25]*t2 + D[26]*t3 + D[27]*t4 + D[28]*t5 + D[29]*t6 + D[30]*t7 + D[31]*t8); */
/*         x[4*i+4]  = (1.0 - omega)*x[4*i+3] + omega*(D[32]*t1 + D[33]*t2 + D[34]*t3 + D[35]*t4 + D[36]*t5 + D[37]*t6 + D[38]*t7 + D[39]*t8); */
/*         x[4*i+5]  = (1.0 - omega)*x[4*i+3] + omega*(D[40]*t1 + D[41]*t2 + D[42]*t3 + D[43]*t4 + D[44]*t5 + D[45]*t6 + D[46]*t7 + D[47]*t8); */
/*         x[4*i+6]  = (1.0 - omega)*x[4*i+3] + omega*(D[48]*t1 + D[49]*t2 + D[50]*t3 + D[51]*t4 + D[52]*t5 + D[53]*t6 + D[54]*t7 + D[55]*t8); */
/*         x[4*i+7]  = (1.0 - omega)*x[4*i+3] + omega*(D[56]*t1 + D[57]*t2 + D[58]*t3 + D[59]*t4 + D[60]*t5 + D[61]*t6 + D[62]*t7 + D[63]*t8); */

/*         // Increment the pointer */
/*         p++; */
/*         i++;; */
/*         D += 64; */
/*       } */
/*       else { */
/*         // Copy the right-hand-side to the temporary vector */
/*         // for this row */
/*         TacsScalar t1, t2, t3, t4; */
/*         t1 = b[4*i]; */
/*         t2 = b[4*i+1]; */
/*         t3 = b[4*i+2]; */
/*         t4 = b[4*i+3]; */
        
/*         // Set the pointer to the beginning of the current */
/*         // row */
/*         const TacsScalar *a = &data->A[16*rowp[i]]; */
        
/*         // Scan through the row and compute the result: */
/*         // tx <- b_i - A_{ij}*x_{j} for j != i */
/*         int end = rowp[i+2]; */
/*         for ( int k = rowp[i]; k < end; k++ ){ */
/*           int j = cols[k]; */
/*           TacsScalar *y = &x[4*j]; */

/*           if (i != j){ */
/*             t1 -= a[0 ]*y[0] + a[1 ]*y[1] + a[2 ]*y[2] + a[3 ]*y[3]; */
/*             t2 -= a[4 ]*y[0] + a[5 ]*y[1] + a[6 ]*y[2] + a[7 ]*y[3]; */
/*             t3 -= a[8 ]*y[0] + a[9 ]*y[1] + a[10]*y[2] + a[11]*y[3]; */
/*             t4 -= a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3]; */
            
/*           } */
          
/*           // Increment the block pointer by bsize^2 */
/*           a += 16; */
/*         } */
        
/*         // Compute the first term in the update: */
/*         // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx */
/*         x[4*i]   = (1.0 - omega)*x[4*i]  + omega*(D[0 ]*t1 + D[1 ]*t2 + D[2 ]*t3 + D[3 ]*t4); */
/*         x[4*i+1] = (1.0 - omega)*x[4*i+1] + omega*(D[4 ]*t1 + D[5 ]*t2 + D[6 ]*t3 + D[7 ]*t4); */
/*         x[4*i+2] = (1.0 - omega)*x[4*i+2] + omega*(D[8 ]*t1 + D[9 ]*t2 + D[10]*t3 + D[11]*t4); */
/*         x[4*i+3] = (1.0 - omega)*x[4*i+3] + omega*(D[12]*t1 + D[13]*t2 + D[14]*t3 + D[15]*t4); */
/*         D += 16; */
/*       } */
/*     } */
/*   } */
/* } */

/* /\*! */
/*   Apply a given number of steps of symmetric SOR to the system A*x = b. */
/* *\/ */
/* void BCSRMatApplySSOR4( BCSRMatData *data, TacsScalar *Adiag, */
/*                         const int *pairs, int npairs, */
/*                         TacsScalar omega, int iters, TacsScalar * b,  */
/*                         TacsScalar * x ){ */
/*   const int nrows = data->nrows; */
/*   const int * rowp = data->rowp; */
/*   const int * cols = data->cols; */

/*   for ( int iter = 0; iter < iters; iter++ ){ */
/*     // Pointer to the last pair */
/*     const TacsScalar *D = Adiag; */
    
/*     for ( int i = 0, p = 0; i < nrows; i++ ){   */
/*       if (p < npairs && i == pairs[p]){ */
/*         // Copy the right-hand-side to the temporary vector */
/*         // for this row */
/*         TacsScalar t1, t2, t3, t4, t5, t6, t7, t8; */
/*         t1 = b[4*i]; */
/*         t2 = b[4*i+1]; */
/*         t3 = b[4*i+2]; */
/*         t4 = b[4*i+3]; */
/*         t5 = b[4*i+4]; */
/*         t6 = b[4*i+5]; */
/*         t7 = b[4*i+6]; */
/*         t8 = b[4*i+7]; */
        
/*         // Set the pointer to the beginning of the current */
/*         // row */
/*         const TacsScalar *a = &data->A[16*rowp[i]]; */
        
/*         // Scan through the row and compute the result: */
/*         // tx <- b_i - A_{ij}*x_{j} for j != i */
/*         int end = rowp[i+1]; */
/*         for ( int k = rowp[i]; k < end; k++ ){ */
/*           int j = cols[k]; */
/*           TacsScalar *y = &x[4*j]; */

/*           if (j != i && j != i+1){ */
/*             t1 -= a[0 ]*y[0] + a[1 ]*y[1] + a[2 ]*y[2] + a[3 ]*y[3]; */
/*             t2 -= a[4 ]*y[0] + a[5 ]*y[1] + a[6 ]*y[2] + a[7 ]*y[3]; */
/*             t3 -= a[8 ]*y[0] + a[9 ]*y[1] + a[10]*y[2] + a[11]*y[3]; */
/*             t4 -= a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3]; */
/*           } */
          
/*           // Increment the block pointer by bsize^2 */
/*           a += 16; */
/*         } */

/*         end = rowp[i+3]; */
/*         for ( int k = rowp[i+1]; k < end; k++ ){ */
/*           int j = cols[k]; */
/*           TacsScalar *y = &x[4*j]; */

/*           if (j != i && j != i+1){ */
/*             t5 -= a[0 ]*y[0] + a[1 ]*y[1] + a[2 ]*y[2] + a[3 ]*y[3]; */
/*             t6 -= a[4 ]*y[0] + a[5 ]*y[1] + a[6 ]*y[2] + a[7 ]*y[3]; */
/*             t7 -= a[8 ]*y[0] + a[9 ]*y[1] + a[10]*y[2] + a[11]*y[3]; */
/*             t8 -= a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3];    */
            
/*           } */
          
/*           // Increment the block pointer by bsize^2 */
/*           a += 16; */
/*         } */

/*         // Compute the first term in the update: */
/*         // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx */
/*         x[4*i  ]  = (1.0 - omega)*x[4*i  ] + omega*(D[0]*t1 + D[1]*t2 + D[2]*t3 + D[3]*t4 + D[4]*t5 + D[5]*t6 + D[6]*t7 + D[7]*t8); */
/*         x[4*i+1]  = (1.0 - omega)*x[4*i+1] + omega*(D[8]*t1 + D[9]*t2 + D[10]*t3 + D[11]*t4 + D[12]*t5 + D[13]*t6 + D[14]*t7 + D[15]*t8); */
/*         x[4*i+2]  = (1.0 - omega)*x[4*i+2] + omega*(D[16]*t1 + D[17]*t2 + D[18]*t3 + D[19]*t4 + D[20]*t5 + D[21]*t6 + D[22]*t7 + D[23]*t8); */
/*         x[4*i+3]  = (1.0 - omega)*x[4*i+3] + omega*(D[24]*t1 + D[25]*t2 + D[26]*t3 + D[27]*t4 + D[28]*t5 + D[29]*t6 + D[30]*t7 + D[31]*t8); */
/*         x[4*i+4]  = (1.0 - omega)*x[4*i+3] + omega*(D[32]*t1 + D[33]*t2 + D[34]*t3 + D[35]*t4 + D[36]*t5 + D[37]*t6 + D[38]*t7 + D[39]*t8); */
/*         x[4*i+5]  = (1.0 - omega)*x[4*i+3] + omega*(D[40]*t1 + D[41]*t2 + D[42]*t3 + D[43]*t4 + D[44]*t5 + D[45]*t6 + D[46]*t7 + D[47]*t8); */
/*         x[4*i+6]  = (1.0 - omega)*x[4*i+3] + omega*(D[48]*t1 + D[49]*t2 + D[50]*t3 + D[51]*t4 + D[52]*t5 + D[53]*t6 + D[54]*t7 + D[55]*t8); */
/*         x[4*i+7]  = (1.0 - omega)*x[4*i+3] + omega*(D[56]*t1 + D[57]*t2 + D[58]*t3 + D[59]*t4 + D[60]*t5 + D[61]*t6 + D[62]*t7 + D[63]*t8); */
       
/*         // Increment the rows/pair number and the pointer to the */
/*         // diagonal block */
/*         p++; */
/*         i++; */
/*         D += 64; */
/*       } */
/*       else { */
/*         // Copy the right-hand-side to the temporary vector */
/*         // for this row */
/*         TacsScalar t1, t2, t3, t4; */
/*         t1 = b[4*i]; */
/*         t2 = b[4*i+1]; */
/*         t3 = b[4*i+2]; */
/*         t4 = b[4*i+3]; */
/*         // Set the pointer to the beginning of the current */
/*         // row */
/*         const TacsScalar *a = &data->A[16*rowp[i]]; */
        
/*         // Scan through the row and compute the result: */
/*         // tx <- b_i - A_{ij}*x_{j} for j != i */
/*         int end = rowp[i+2]; */
/*         for ( int k = rowp[i]; k < end; k++ ){ */
/*           int j = cols[k]; */
/*           TacsScalar *y = &x[4*j]; */

/*           if (i != j){ */
/*             t1 -= a[0 ]*y[0] + a[1 ]*y[1] + a[2 ]*y[2] + a[3 ]*y[3]; */
/*             t2 -= a[4 ]*y[0] + a[5 ]*y[1] + a[6 ]*y[2] + a[7 ]*y[3]; */
/*             t3 -= a[8 ]*y[0] + a[9 ]*y[1] + a[10]*y[2] + a[11]*y[3]; */
/*             t4 -= a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3]; */
/*           } */
          
/*           // Increment the block pointer by bsize^2 */
/*           a += 16; */
/*         } */
        
/*         // Compute the first term in the update: */
/*         // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx */
/*         x[4*i]   = (1.0 - omega)*x[4*i]  + omega*(D[0 ]*t1 + D[1 ]*t2 + D[2 ]*t3 + D[3 ]*t4); */
/*         x[4*i+1] = (1.0 - omega)*x[4*i+1] + omega*(D[4 ]*t1 + D[5 ]*t2 + D[6 ]*t3 + D[7 ]*t4); */
/*         x[4*i+2] = (1.0 - omega)*x[4*i+2] + omega*(D[8 ]*t1 + D[9 ]*t2 + D[10]*t3 + D[11]*t4); */
/*         x[4*i+3] = (1.0 - omega)*x[4*i+3] + omega*(D[12]*t1 + D[13]*t2 + D[14]*t3 + D[15]*t4); */
/*         D += 16; */
/*       } */
/*     } */

/*     for ( int i = nrows-1, p = npairs-1; i >= 0; i-- ){ */
/*       if (p >= 0 && i-1 == pairs[p]){ */
/*         // Copy the right-hand-side to the temporary vector */
/*         // for this row */
/*         TacsScalar t1, t2, t3, t4, t5, t6, t7, t8; */
/*         t1 = b[4*i-4]; */
/*         t2 = b[4*i-3]; */
/*         t3 = b[4*i-2]; */
/*         t4 = b[4*i-1]; */
/*         t5 = b[4*i]; */
/*         t6 = b[4*i+1]; */
/*         t7 = b[4*i+2]; */
/*         t8 = b[4*i+3]; */
        
/*         // Set the pointer to the beginning of the current */
/*         // row */
/*         const TacsScalar *a = &data->A[16*rowp[i-1]]; */
        
/*         // Scan through the row and compute the result: */
/*         // tx <- b_i - A_{ij}*x_{j} for j != i */
/*         int end = rowp[i]; */
/*         // -------------- Change? */
/*         for ( int k = rowp[i-1]; k < end; k++ ){ */
/*           int j = cols[k]; */
/*           TacsScalar *y = &x[4*j]; */

/*           if (j != i-1 && j != i){ */
/*             t1 -= a[0 ]*y[0] + a[1 ]*y[1] + a[2 ]*y[2] + a[3 ]*y[3]; */
/*             t2 -= a[4 ]*y[0] + a[5 ]*y[1] + a[6 ]*y[2] + a[7 ]*y[3]; */
/*             t3 -= a[8 ]*y[0] + a[9 ]*y[1] + a[10]*y[2] + a[11]*y[3]; */
/*             t4 -= a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3]; */
            
/*           } */
          
/*           // Increment the block pointer by bsize^2 */
/*           a += 16; */
/*         } */

/*         end = rowp[i+2]; */
/*         for ( int k = rowp[i]; k < end; k++ ){ */
/*           int j = cols[k]; */
/*           TacsScalar *y = &x[4*j]; */

/*           if (j != i-1 && j != i){ */
/*             t5 -= a[0 ]*y[0] + a[1 ]*y[1] + a[2 ]*y[2] + a[3 ]*y[3]; */
/*             t6 -= a[4 ]*y[0] + a[5 ]*y[1] + a[6 ]*y[2] + a[7 ]*y[3]; */
/*             t7 -= a[8 ]*y[0] + a[9 ]*y[1] + a[10]*y[2] + a[11]*y[3]; */
/*             t8 -= a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3];    */
            
/*           } */
          
/*           // Increment the block pointer by bsize^2 */
/*           a += 16; */
/*         } */

/*         // Compute the first term in the update: */
/*         // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx */
/*         D -= 64; */
/*         x[4*i-4]  = (1.0 - omega)*x[4*i-4] + omega*(D[0]*t1 + D[1]*t2 + D[2]*t3 + D[3]*t4 + D[4]*t5 + D[5]*t6 + D[6]*t7 + D[7]*t8); */
/*         x[4*i-3]  = (1.0 - omega)*x[4*i-3] + omega*(D[8]*t1 + D[9]*t2 + D[10]*t3 + D[11]*t4 + D[12]*t5 + D[13]*t6 + D[14]*t7 + D[15]*t8); */
/*         x[4*i-2]  = (1.0 - omega)*x[4*i-2] + omega*(D[16]*t1 + D[17]*t2 + D[18]*t3 + D[19]*t4 + D[20]*t5 + D[21]*t6 + D[22]*t7 + D[23]*t8); */
/*         x[4*i-1]  = (1.0 - omega)*x[4*i-1] + omega*(D[24]*t1 + D[25]*t2 + D[26]*t3 + D[27]*t4 + D[28]*t5 + D[29]*t6 + D[30]*t7 + D[31]*t8); */
/*         x[4*i  ]  = (1.0 - omega)*x[4*i  ] + omega*(D[32]*t1 + D[33]*t2 + D[34]*t3 + D[35]*t4 + D[36]*t5 + D[37]*t6 + D[38]*t7 + D[39]*t8); */
/*         x[4*i+1]  = (1.0 - omega)*x[4*i+1] + omega*(D[40]*t1 + D[41]*t2 + D[42]*t3 + D[43]*t4 + D[44]*t5 + D[45]*t6 + D[46]*t7 + D[47]*t8); */
/*         x[4*i+2]  = (1.0 - omega)*x[4*i+2] + omega*(D[48]*t1 + D[49]*t2 + D[50]*t3 + D[51]*t4 + D[52]*t5 + D[53]*t6 + D[54]*t7 + D[55]*t8); */
/*         x[4*i+3]  = (1.0 - omega)*x[4*i+3] + omega*(D[56]*t1 + D[57]*t2 + D[58]*t3 + D[59]*t4 + D[60]*t5 + D[61]*t6 + D[62]*t7 + D[63]*t8); */
        
/*         // Increment the rows/pair number and the pointer to the */
/*         // diagonal block */
/*         p--; */
/*         i--; */
/*       } */
/*       else { */
/*         // Copy the right-hand-side to the temporary vector */
/*         // for this row */
/*         TacsScalar t1, t2, t3, t4; */
/*         t1 = b[4*i]; */
/*         t2 = b[4*i+1]; */
/*         t3 = b[4*i+2]; */
/*         t4 = b[4*i+3]; */
/*         // Set the pointer to the beginning of the current */
/*         // row */
/*         const TacsScalar *a = &data->A[16*rowp[i]]; */
        
/*         // Scan through the row and compute the result: */
/*         // tx <- b_i - A_{ij}*x_{j} for j != i */
/*         int end = rowp[i+2]; */
/*         for ( int k = rowp[i]; k < end; k++ ){ */
/*           int j = cols[k]; */
/*           TacsScalar *y = &x[4*j]; */

/*           if (i != j){ */
/*             t1 -= a[0 ]*y[0] + a[1 ]*y[1] + a[2 ]*y[2] + a[3 ]*y[3]; */
/*             t2 -= a[4 ]*y[0] + a[5 ]*y[1] + a[6 ]*y[2] + a[7 ]*y[3]; */
/*             t3 -= a[8 ]*y[0] + a[9 ]*y[1] + a[10]*y[2] + a[11]*y[3]; */
/*             t4 -= a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3]; */
/*           } */
          
/*           // Increment the block pointer by bsize^2 */
/*           a += 16; */
/*         } */
        
/*         // Compute the first term in the update: */
/*         // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx */
/*         D -= 16; */
/*         x[4*i]   = (1.0 - omega)*x[4*i]  + omega*(D[0 ]*t1 + D[1 ]*t2 + D[2 ]*t3 + D[3 ]*t4); */
/*         x[4*i+1] = (1.0 - omega)*x[4*i+1] + omega*(D[4 ]*t1 + D[5 ]*t2 + D[6 ]*t3 + D[7 ]*t4); */
/*         x[4*i+2] = (1.0 - omega)*x[4*i+2] + omega*(D[8 ]*t1 + D[9 ]*t2 + D[10]*t3 + D[11]*t4); */
/*         x[4*i+3] = (1.0 - omega)*x[4*i+3] + omega*(D[12]*t1 + D[13]*t2 + D[14]*t3 + D[15]*t4); */
        
/*       } */
/*     } */
/*   } */
/* } */
