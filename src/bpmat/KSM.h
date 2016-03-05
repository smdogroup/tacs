#ifndef TACS_KSM_H
#define TACS_KSM_H

/*
  Solve linear systems using Krylov-subspace methods
  
  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

#include <math.h>
#include "TACSObject.h"

/*!
  The abstract vector class. 

  All TACS KSM objects use this interface for matrix-vector products
  and preconditioning operations. In practice, within TACS itself,
  most of the matricies require the BVec base class. However, when
  coupling to other disciplines it is frequently convenient to define
  a hybrid vector consisting of multiple components. In these cases,
  it is handy to have a base vector class. 
*/
class TACSVec : public TACSObject {
 public:
  virtual ~TACSVec(){}

  // Functions required for Krylov methods
  // -------------------------------------
  virtual TacsScalar norm() = 0; // Compute the Cartesian 2 norm
  virtual void scale( TacsScalar alpha ) = 0; // Scale the vector by a value
  virtual TacsScalar dot( TACSVec * x ) = 0;  // Compute x^{T} * y
  virtual void mdot( TACSVec ** x, TacsScalar * ans, int m ); // Multiple dot product
  virtual void axpy( TacsScalar alpha, TACSVec * x ) = 0; // Compute y <- y + alpha * x
  virtual void copyValues( TACSVec * x ) = 0; // Copy values from x to this
  virtual void axpby( TacsScalar alpha, TacsScalar beta, 
                      TACSVec * x ) = 0; // Compute y <- alpha * x + beta * y 
  virtual void zeroEntries() = 0; // Zero all the entries

  // Additional useful member functions
  // ----------------------------------
  virtual void setRand( double lower = -1.0, double upper = 1.0 ){}
  virtual void applyBCs(){}
};

/*!  
  The abstract matrix base class for all TACS matrices. All of
  these operations must be defined to be utilized by TACSAssembler.

  zeroEntries(): Zeros the values of all entries in the matrix.

  addValues(): Adds a small, dense matrix to the rows set in row[] and
  the columns set in col[].

  applyBCs(): Applies the Dirichlet boundary conditions to the matrix
  by settin the associated diagonal elements to 1.

  beginAssembly()/endAssembly(): Begin/end the parallel assembly of
  the matrix.

  getSize(): Return the number of rows and number of columns in the
  matrix.

  createVec(): Return a vector for output from the matrix - typically
  only used when the matrix is square.

  mult(x, y): Perform the matrix multiplication y = A*x

  copyValues(): Copy the values from the matrix mat to this matrix.

  scale(alpha): Scale all entries in the matrix by alpha.

  axpy(alpha, mat): Compute this = this + alpha*mat - only works on
  the same matrix types.
*/
class TACSMat : public TACSObject {
 public:
  virtual ~TACSMat(){}
  
  // Operations for assembly/setting values
  // --------------------------------------  
  virtual void zeroEntries(){}
  virtual void addValues( int nrow, const int * row, 
			  int ncol, const int * col,
			  int nv, int mv, const TacsScalar * values ){}
  virtual void addWeightValues( int nvars, const int *varp, const int *vars,
				const TacsScalar *weights,
				int nv, int mv, const TacsScalar *values ){}
  virtual void applyBCs(){}
  virtual void beginAssembly(){}
  virtual void endAssembly(){}

  // Create vectors/retrieve sizes
  // -----------------------------
  virtual void getSize( int * _nr, int * _nc ){}
  virtual TACSVec * createVec() = 0;

  // Operations required for solving problems
  // ----------------------------------------
  virtual void mult( TACSVec * x, TACSVec * y ) = 0;
  virtual void copyValues( TACSMat * mat ){}
  virtual void scale( TacsScalar alpha ){}
  virtual void axpy( TacsScalar alpha, TACSMat * mat ){}

  // Return the name of the object
  // -----------------------------
  const char * TACSObjectName();

 private:
  static const char * matName;
};

/*!
  The abstract pre-conditioner class

  applyFactor(): Compute y = M^{-1}x, where M^{-1} is the
  preconditioner.

  factor(): Factor the preconditioner based on values in the matrix
  associated with the preconditioner
*/
class TACSPc : public TACSObject {
 public:
  virtual ~TACSPc(){}
  // Apply the preconditioner to x, to produce y
  // -------------------------------------------
  virtual void applyFactor( TACSVec * x, TACSVec * y ) = 0;

  // Factor (or set up) the preconditioner 
  // -------------------------------------
  virtual void factor() = 0;

  // Retrieve the object name
  // ------------------------
  const char * TACSObjectName();

 private:
  static const char * pcName;
};

/*!
  The abstract residual history class

  Monitor the residual or other interesting quantity and print out the
  result to the screen or file.
*/
class KSMPrint : public TACSObject {
 public:
  virtual ~KSMPrint(){}

  virtual void printResidual( int iter, TacsScalar res ) = 0;
  virtual void print( const char * cstr ) = 0;
  const char * TACSObjectName();

 private:
  static const char * printName;
};

/*!
  The abstract Krylov-subspace method class
  
  Solve the linear system A*x = b with a Krylov subspace method,
  possibly using some preconditioner.

  createVec(): Create a vector for the matrix

  setOperations(): Set the matrix (required) and the preconditioner
  (optional)

  getOperators(): Get the matrix and the preconditioner

  solve(): Solve the linear system to the specified tolerance using
  the Krylov subspace method.

  setTolerances(rtol, atol): Set the relative and absolute stopping
  tolerances for the method

  setMonitor(): Set the monitor - possibly NULL - that will be used
 */
class TACSKsm : public TACSObject {
 public:
  virtual ~TACSKsm(){}

  virtual TACSVec * createVec() = 0;
  virtual void setOperators( TACSMat * _mat, TACSPc * _pc ) = 0;
  virtual void getOperators( TACSMat ** _mat, TACSPc ** _pc ) = 0;
  virtual void solve( TACSVec * b, TACSVec * x, int zero_guess = 1 ) = 0;
  virtual void setTolerances( double _rtol, double _atol ) = 0;
  virtual void setMonitor( KSMPrint * _monitor ) = 0;
  const char * TACSObjectName();

 private:
  static const char * ksmName;
};

/*
  The following classes define a consistent way to write
  residual histories to a file/to stdout
*/
class KSMPrintStdout : public KSMPrint {
 public:
  KSMPrintStdout( const char * _descript, int _rank, int _freq );
  ~KSMPrintStdout();
  
  void printResidual( int iter, TacsScalar res );
  void print( const char * cstr );

 private:
  char * descript;
  int rank;
  int freq;
};

class KSMPrintFile : public KSMPrint {
  KSMPrintFile( const char * filename, const char * _descript, 
		int _rank, int _freq );
  ~KSMPrintFile();
  
  void printResidual( int iter, TacsScalar res );
  void print( const char * cstr );

 private:
  FILE * fp;
  char * descript;
  int rank;
  int freq;
};

/*!
  Pre-conditioned conjugate gradient method
*/
class PCG : public TACSKsm {
 public:
  PCG( TACSMat * _mat, TACSPc * _pc, int reset, int _nouter );
  ~PCG();

  TACSVec * createVec(){ return mat->createVec(); }
  void solve( TACSVec * b, TACSVec * x, int zero_guess = 1 );
  void setOperators( TACSMat * _mat, TACSPc * _pc );
  void getOperators( TACSMat ** _mat, TACSPc ** _pc );
  void setTolerances( double _rtol, double _atol );
  void setMonitor( KSMPrint * _print );

 private:
  TACSMat * mat;
  TACSPc * pc;

  double rtol;
  double atol;
  
  int nouter;
  int reset;
  
  TACSVec * work;
  TACSVec * R;
  TACSVec * P;
  TACSVec * Z;

  KSMPrint * monitor;
};

/*!
  Right-preconditioned GMRES 

  This class handles both right-preconditioned GMRES and the flexible
  variant of GMRES. Either the classical or modified Gram Schmidt
  orthogonalization algorithm may be used. Modified Gram Schmidt has
  better numerical stability properties, but requires more
  communication overhead. Modified Gram Schmidt is the default.

  The input parameters are:
  -------------------------
  mat: the matrix handle
  pc: (optional - but highly recommended!) the preconditioner
  m: the size of the Krylov-subspace to use before restarting
  nrestart: the number of restarts to use nrestart = 0 results in

  isflexible: flag to indicate whether to use a flexible variant of
  GMRES or not

  Notes:
  ------
  - The preconditioner is used as is without any additional
  calls. Ensure that it is properly factored by calling pc->factor().

  - When calling solve, set zero_guess = 0 (False) to use the values
  in 'x' as the initial guess.

  - You can change the matrix and preconditioner operators by making
  calls to setOperators - this is useful if you want to re-use the
  vectors in the Krylov subspace, but have two different problems.
*/
class GMRES : public TACSKsm {
 public:
  enum OrthoType { CLASSICAL_GRAM_SCHMIDT, MODIFIED_GRAM_SCHMIDT };
  GMRES( TACSMat * _mat, TACSPc * _pc, int _m, int _nrestart, int _isFlexible );
  GMRES( TACSMat * _mat, int _m, int _nrestart );
  ~GMRES();

  TACSVec * createVec(){ return mat->createVec(); }  
  void solve( TACSVec * b, TACSVec * x, int zero_guess = 1 );
  void setOperators( TACSMat * _mat, TACSPc * _pc );
  void getOperators( TACSMat ** _mat, TACSPc ** _pc );
  void setTolerances( double _rtol, double _atol );
  void setMonitor( KSMPrint * _monitor );
  void setOrthoType( enum OrthoType otype );
  void setTimeMonitor();

  const char * TACSObjectName();

 private:
  // Initialize the class
  void init( TACSMat * _mat, TACSPc * _pc, int _m, 
	     int _nrestart, int _isFlexible );

  // Orthogonalize a vector against a set of vectors
  void (*orthogonalize)(TacsScalar *, TACSVec *, TACSVec **, int);

  TACSMat * mat; 
  TACSPc * pc; 
  int msub;
  int nrestart;
  int isFlexible;

  TACSVec ** W;   // The Arnoldi vectors that span the Krylov subspace
  TACSVec ** Z;   // An additional subspace of vectors -- for the flexible variant
  TACSVec * work; // A work vector

  int * Hptr;   // A vector to make accessing the elements of the matrix easier!
  TacsScalar * H; // The Hessenberg matrix

  double rtol;
  double atol;

  TacsScalar * Qsin;
  TacsScalar * Qcos;
  TacsScalar * res;

  int monitor_time;
  KSMPrint * monitor;

  static const char * gmresName;
};

/*!
  A simplified and flexible variant of GCROT - from Hicken and Zingg

  This is a simplified and flexible variant of the GCROT algorithm
  suggested by Hicken and Zingg. This Krylov subspace method - GCR
  with outer truncation - has better restart properties than F/GMRES.

  The input parameters are:
  mat: The matrix operator
  pc: The preconditioner (optional but recommended)
  outer: The number of outer subspace vectors to use
  max_outer >= outer: The maximum number of outer iterations
  msub: The size of the GMRES subspace to use
  isflexible: Parameter to indicate whether to use a flexible variant
*/
class GCROT : public TACSKsm {
 public:
  GCROT( TACSMat * _mat, TACSPc * _pc, int _outer, int _max_outer, 
	 int _msub, int _isFlexible );
  GCROT( TACSMat * _mat, int _outer, int _max_outer, int _msub );
  ~GCROT();

  TACSVec * createVec(){ return mat->createVec(); }
  void solve( TACSVec * b, TACSVec * x, int zero_guess = 1 );
  void setOperators( TACSMat * _mat, TACSPc * _pc );
  void getOperators( TACSMat ** _mat, TACSPc ** _pc );
  void setTolerances( double _rtol, double _atol );
  void setMonitor( KSMPrint * _monitor );

  const char * TACSObjectName();

 private:

  void init( TACSMat * _mat, TACSPc * _pc, int _outer, int _max_outer, 
	     int _msub, int _isFlexible );

  TACSMat * mat; 
  TACSPc * pc; 
  int msub;             // Size of the GMRES subspace
  int outer, max_outer; // Number of outer vectors
  int isFlexible;

  TACSVec ** W;     // The Arnoldi vectors that span the Krylov subspace
  TACSVec ** Z;     // An additional subspace of vectors - for the flexible variant
  TACSVec **U, **C; // The subspaces for the outer GCR iterations
  TACSVec *u_hat, *c_hat, *R;

  int * Hptr;   // A vector to make accessing the elements of the matrix easier!
  TacsScalar * H; // The Hessenberg matrix
  TacsScalar * B; // The matrix that stores C^{T} * A * W

  double rtol;
  double atol;

  TacsScalar * Qsin;
  TacsScalar * Qcos;
  TacsScalar * res;

  KSMPrint * monitor;

  static const char * gcrotName;
};

/*!
  Right-preconditioned GMRES with a constraint Cx = 0

  Determine the solution of a system of linear equations using GMRES
  subject to a set of constraints that can be set later.
*/
class ConGMRES : public TACSKsm {
 public:
  enum OrthoType { CLASSICAL_GRAM_SCHMIDT, MODIFIED_GRAM_SCHMIDT };
  ConGMRES( TACSMat * _mat, TACSPc * _pc, 
            int _m, int _nrestart, int _isFlexible );
  ConGMRES( TACSMat * _mat, int _m, int _nrestart );
  ~ConGMRES();

  TACSVec * createVec(){ return mat->createVec(); }
  void solve( TACSVec * b, TACSVec * x, int zero_guess = 1 ); 
  void setOperators( TACSMat * _mat, TACSPc * _pc );
  void getOperators( TACSMat ** _mat, TACSPc ** _pc );
  void setTolerances( double _rtol, double _atol );
  void setOrthoType( enum OrthoType otype );
  void setMonitor( KSMPrint * _monitor );
  void setTimeMonitor();

  void setConstraints( TACSVec ** _C, int nconstr );

  const char * TACSObjectName();

 private:
  // Initialize the class
  void init( TACSMat * _mat, TACSPc * _pc, int _m, 
	     int _nrestart, int _isFlexible );

  // Orthogonalize a vector against a set of vectors
  void (*orthogonalize)(TacsScalar *, TACSVec *, TACSVec **, int);

  TACSMat * mat; 
  TACSPc * pc; 
  int msub;
  int nrestart;
  int isFlexible;

  // Constraint information
  int nconstr;     // The number of constraints
  TACSVec ** C;       // The set of constraints (that are orthonormalized on input)
  TacsScalar * lambda;

  void con_mat_mult( TACSVec * x, TacsScalar * y, 
		     TACSVec * out, TacsScalar * out_lambda );
  TacsScalar con_norm( TACSVec * x, TacsScalar * y );
  TacsScalar con_dot( TACSVec * x, TacsScalar * y,
		      TACSVec * X, TacsScalar * Y );
  void con_axpy( TACSVec * x, TacsScalar * X, 
		 TacsScalar alpha,
		 TACSVec * y, TacsScalar * Y );
  void con_scale( TACSVec * x, TacsScalar * y, TacsScalar alpha );

  TACSVec ** W;   // The Arnoldi vectors that span the Krylov subspace
  TacsScalar * Wl; // The portion of the Arnoldi vector for lambda
  TACSVec ** Z;   // An additional subspace of vectors -- for the flexible variant
  TACSVec * work; // A work vector

  int * Hptr;   // A vector to make accessing the elements of the matrix easier!
  TacsScalar * H; // The Hessenberg matrix

  double rtol;
  double atol;

  TacsScalar * Qsin;
  TacsScalar * Qcos;
  TacsScalar * res;

  int monitor_time;
  KSMPrint * monitor;

  static const char * gmresName;
};

/*
  Create a Krylov-subspace class that is just a preconditioner.

  This is used in cases where a preconditioner may consist of an
  internal subspace method. In these cases, it is sometimes useful to
  have a KSM class that just applies the preconditioner itself. This
  class fills this role.
*/
class KsmPreconditioner : public TACSKsm {
 public:
  KsmPreconditioner( TACSMat * _mat, TACSPc * _pc );
  ~KsmPreconditioner();

  TACSVec * createVec();
  void setOperators( TACSMat * _mat, TACSPc * _pc );
  void getOperators( TACSMat ** _mat, TACSPc ** _pc );
  void solve( TACSVec * b, TACSVec * x, int zero_guess = 1 );
  void setTolerances( double _rtol, double _atol );
  void setMonitor( KSMPrint * _monitor );

 private:
  TACSMat * mat;
  TACSPc * pc;
};

#endif
