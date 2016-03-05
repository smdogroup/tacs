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
  PMat( VarMap * rmap,
	BCSRMat * _Aloc, BCSRMat * _Bext,
	BVecDistribute * _col_map,
	BCMap * _bcs = NULL );
  ~PMat();

  // Functions for setting values in the matrix
  // ------------------------------------------
  void zeroEntries();  // Zero the matrix values 
  void applyBCs();     // Apply the boundary conditions to the matrix

  // Functions required for solving linear systems
  // ---------------------------------------------
  void getSize( int * _nr, int * _nc ); // Get the local dimensions of the matrix
  void mult( TACSVec * x, TACSVec * y );   // y <- A*x 
  TACSVec * createVec(); // Create a vector to operate on

  void copyValues( TACSMat * mat ); // Copy the values from another PMat
  void scale( TacsScalar alpha ); // Scale the matrix by a value
  void axpy( TacsScalar alpha, TACSMat * mat ); // Compute y <- y + alpha*x
  void axpby( TacsScalar alpha, 
	      TacsScalar beta, TACSMat * mat ); // Compute y <- alpha*x + beta*y

  void addDiag( TacsScalar alpha );

  // Other miscelaneous functions
  // ----------------------------
  void getBCSRMat( BCSRMat ** A, BCSRMat ** B ); // Access the underlying mats
  void getRowMap( int * bs, int * _N, int * _Nc );
  void getColMap( int * bs, int * _M );
  VarMap * getRowMap(){ return rmap; }
  void getExtColMap( BVecDistribute ** ext_map );  // Access the column map  
  void printNzPattern( const char * fileName );    // Print the non-zero pattern

  const char * TACSObjectName();
  MPI_Comm getMPIComm(){ return rmap->getMPIComm(); }

 protected:
  PMat();

  void init( VarMap * _rmap, 
	     BCSRMat * _Aloc, BCSRMat * _Bext,
	     BVecDistribute * _col_map,
	     BCMap * _bcs = NULL );

  // Matrix data
  // -----------
  VarMap * rmap; // Map for the rows 
  BCSRMat * Aloc;
  BCSRMat * Bext;
  BVecDistribute * col_map; 

  // Variable data
  // -------------
  int bsize; // The block size
  int N;     // The number of local rows
  int Nc;    // The number of equations that are coupled to other processors
  int Np;    // The number of local-only equations Np + Nc = N

 private: 

  // BC information
  // --------------
  BCMap * bcs;
    
  TacsScalar * x_ext; // External values - used for matrix-vector products
  int ext_offset;

  static const char * matName;
};

/*
  Parallel successive over relaxation (PSOR)
  ------------------------------------------
  
  Set up involves factoring the diagonal matrix
*/
class PSOR : public TACSPc {
 public:
  PSOR( PMat * mat, int _zero_guess, 
	TacsScalar _omega, int _iters, int _isSymmetric );
  ~PSOR();

  void factor();
  void applyFactor( TACSVec * xvec, TACSVec * yvec );

 private:
  int zero_guess;
  TacsScalar omega;
  int iters;
  int isSymmetric;

  BCSRMat * Aloc;
  BCSRMat * Bext;

  BVec * bvec;

  BVecDistribute * col_map;
  int ext_offset;
  TacsScalar * yext;
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
  AdditiveSchwarz( PMat * mat, int levFill, double fill );
  ~AdditiveSchwarz();

  void setDiagShift( TacsScalar _alpha );
  void factor();
  void applyFactor( TACSVec * xvec, TACSVec * yvec );
  void applyFactor( TACSVec * yvec );

 private:
  BCSRMat * Aloc;
  TacsScalar alpha;
  BCSRMat * Apc;
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
  GlobalSchurMat( PMat * mat, BCSRMat * Apc );
  ~GlobalSchurMat();

  // These functions have no effect for this matrix
  // ----------------------------------------------
  void zeroEntries(){}
  void addValues( int nrow, const int * row, int ncol, const int * col,
		  int nv, int mv, const TacsScalar * values ){} 
  void applyBCs(){}
  void beginAssembly(){}
  void endAssembly(){}

  // Functions used to solve the linear system
  // -----------------------------------------
  void getSize( int * _nr, int * _nc );  // Get the local dimensions of the matrix
  void mult( TACSVec * x, TACSVec * y );  // y <- A * x 
  TACSVec * createVec(); // Create a vector to operate on

  // Multiply y <- Bext * x_ext 
  void multOffDiag( BVec * x, BVec * y );

 private:
  VarMap * rmap;  // The variable map for the interface variables
  BCSRMat * Apc;  // The factored diagonal matrix 
  BCSRMat * Bext; // The off-diagonal part
  BVecDistribute * col_map; 

  int varoffset;
  int nvars;
  TacsScalar * x_ext;
};

/*!
  The approximate Schur preconditioner
*/
class ApproximateSchur : public TACSPc {
 public:
  ApproximateSchur( PMat * mat, int levFill, double fill, 
		    int inner_gmres_iters, double inner_rtol = 1e-3, 
		    double inner_atol = 1e-30 );
  ~ApproximateSchur();

  void setDiagShift( TacsScalar _alpha );
  void setMonitor( KSMPrint * ksm_print );
  void factor();
  void applyFactor( TACSVec * xvec, TACSVec * yvec );
  void printNzPattern( const char * fileName );

 private:
  BVec *rvec, *wvec;
  PMat * mat;
  BCSRMat * Aloc;
  TacsScalar alpha;
  BCSRMat * Apc;
 
  int start;
  int end;
  int var_offset;

  GlobalSchurMat * gsmat;
  TACSKsm * inner_ksm;
};

#endif
