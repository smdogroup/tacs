#ifndef TACS_PARALLEL_BCSR_MAT_H
#define TACS_PARALLEL_BCSR_MAT_H

/*
  Parallel matrix code

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

#include "BVec.h"
#include "BVecDist.h"
#include "BCSRMat.h"
#include "KSM.h"

/*!
  Parallel matrix class with Block Compressed-Sparse row format.

  Each processor contains the following information:

  1. The local matrix that acts on the local variables (A_loc)
  2. The external-coupling matrix (B_ext)
  3. A map that takes the global vector x and maps it to the local external
  vector components: x_ext

  Matrix-vector products are done row-wise:

  y_loc = A_loc * x_loc + B_ext * x_ext

  Matrix-vector products
  ----------------------
  1. Begin scatter operation from global to local-external unknowns
  2. Perform local matrix-vector product
  3. End scatter operation 
  4. Perform local, exteral matrix-vector product  
*/
class PMat : public TACSMat {
 public:
  PMat( TACSVarMap *rmap,
        BCSRMat *_Aloc, BCSRMat *_Bext,
        TACSBVecDistribute *_col_map );
  ~PMat();

  // Functions for setting values in the matrix
  // ------------------------------------------
  void zeroEntries();                // Zero the matrix values 
  void applyBCs( TACSBcMap *bcmap ); // Apply the boundary conditions

  // Functions required for solving linear systems
  // ---------------------------------------------
  void getSize( int *_nr, int *_nc );          // Get the local dimensions
  void mult( TACSVec *x, TACSVec *y );         // y <- A*x 
  TACSVec *createVec();                        // Create a vector
  void copyValues( TACSMat *mat );             // Copy matrix entries
  void scale( TacsScalar alpha );              // Scale the matrix
  void axpy( TacsScalar alpha, TACSMat *mat ); // Compute y <- y + alpha*x
  void axpby( TacsScalar alpha, 
              TacsScalar beta, TACSMat *mat ); // Compute y <- alpha*x + beta*y

  void addDiag( TacsScalar alpha );

  // Other miscelaneous functions
  // ----------------------------
  void getBCSRMat( BCSRMat ** A, BCSRMat ** B );  // Access the underlying mats
  void getRowMap( int *bs, int *_N, int *_Nc );
  void getColMap( int *bs, int *_M );
  TACSVarMap *getRowMap(){ return rmap; }
  void getExtColMap( TACSBVecDistribute **ext_map ); // Access the column map  
  void printNzPattern( const char *fileName ); // Print the non-zero pattern

  const char *TACSObjectName();
  MPI_Comm getMPIComm(){ return rmap->getMPIComm(); }

 protected:
  PMat();

  // Common initialization routine
  void init( TACSVarMap *_rmap, 
             BCSRMat *_Aloc, BCSRMat *_Bext,
             TACSBVecDistribute *_col_map );

  // Local entries for the matrix
  BCSRMat *Aloc, *Bext;
    
  // Map the local entries into the global data
  TACSVarMap *rmap;
  TACSBVecDistribute *ext_dist; 
  TACSBVecDistCtx *ctx;

  // Sizes/dimensions of the matrix
  int bsize; // The block size
  int N; // Number of local rows
  int Nc; // Number of equations that are coupled to other processors
  int Np; // The number of local-only equations Np + Nc = N

 private: 
  // External values - used for matrix-vector products
  TacsScalar *x_ext; 
  int ext_offset;

  static const char *matName;
};

/*
  Parallel successive over relaxation (PSOR)
  ------------------------------------------
  
  Set up involves factoring the diagonal matrix
*/
class PSOR : public TACSPc {
 public:
  PSOR( PMat *_mat, int _zero_guess, TacsScalar _omega, int _iters, 
        int _isSymmetric, int *pairs=NULL, int npairs=0 );
  ~PSOR();

  void factor();
  void applyFactor( TACSVec *xvec, TACSVec *yvec );
  void getMat( TACSMat **_mat );

 private:
  // Parallel matrix pointer
  PMat *mat;

  // Information about how to handle the smoother
  int iters;
  TacsScalar omega;
  int zero_guess, isSymmetric;

  // Pointers to the local/external matrix
  BCSRMat *Aloc, *Bext;
  TACSBVec *bvec;

  // Parallel data for the PSOR object
  TACSBVecDistribute *ext_dist;
  TACSBVecDistCtx *ctx;
  int ext_offset;
  TacsScalar *yext;
};

/*
  Additive Schwarz Method (ASM)
  -----------------------------

  Set up involves factoring the diagonal portion of the matrix
  
  ASM: Apply preconditioner
  Apply the local preconditioner to the local components of the residual. 
*/
class AdditiveSchwarz : public TACSPc {
 public:
  AdditiveSchwarz( PMat *mat, int levFill, double fill );
  ~AdditiveSchwarz();

  void setDiagShift( TacsScalar _alpha );
  void factor();
  void applyFactor( TACSVec *xvec, TACSVec *yvec );
  void applyFactor( TACSVec *yvec );
  void getMat( TACSMat **_mat );

 private:
  PMat *mat;
  BCSRMat *Aloc;
  TacsScalar alpha;
  BCSRMat *Apc;
};

/*
  Approximate Schur 
  -----------------

  Set up involves factoring the diagonal portion of the matrix

  AS: Apply preconditioner
  1. Restrict to the interface unknowns
  2. Solve a GMRES-accelerated, Jacobi-preconditioned problem for the interface
  unknowns
  3. Determine the solution at the internal interface unknowns
*/
class GlobalSchurMat : public TACSMat {
 public:
  GlobalSchurMat( PMat *mat, BCSRMat *Apc );
  ~GlobalSchurMat();

  // Functions used to solve the linear system
  // -----------------------------------------
  void getSize( int *_nr, int *_nc );
  void mult( TACSVec *x, TACSVec *y );
  TACSVec *createVec();

  // Multiply y <- Bext *x_ext 
  void multOffDiag( TACSBVec *x, TACSBVec *y );

 private:
  TACSVarMap *rmap; // The variable map for the interface variables
  BCSRMat *Apc; // The factored diagonal matrix 
  BCSRMat *Bext; // The off-diagonal part
  TACSBVecDistribute *ext_dist; 
  TACSBVecDistCtx *ctx;

  int nvars, varoffset;
  TacsScalar *x_ext;
};

/*!
  The approximate Schur preconditioner
*/
class ApproximateSchur : public TACSPc {
 public:
  ApproximateSchur( PMat *mat, int levFill, double fill, 
                    int inner_gmres_iters, double inner_rtol=1e-3, 
                    double inner_atol=1e-30 );
  ~ApproximateSchur();

  void setDiagShift( TacsScalar _alpha );
  void setMonitor( KSMPrint *ksm_print );
  void factor();
  void applyFactor( TACSVec *xvec, TACSVec *yvec );
  void getMat( TACSMat **_mat );
  void printNzPattern( const char *fileName );

 private:
  TACSBVec *rvec, *wvec;
  PMat *mat;
  BCSRMat *Aloc, *Apc;
  TacsScalar alpha;
 
  // Offsets into the array
  int start, end, var_offset;

  // Global Schur matrix and its associated KSM object
  GlobalSchurMat *gsmat;
  TACSKsm *inner_ksm;
};

#endif
