#include <stdio.h>
#include <math.h>
#include "KSM.h"

/*
  Implementation of various Krylov-subspace methods

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

/*
  The default implementation for performing multiple dot-products
  simultaneously. This operation can be implemented more efficiently
  in parallel by performing all local operations, then performing the
  communication step.

  input:
  x:   the vectors
  ans: an array of the dot product results
  m:   the number of vectors in x
*/
void TACSVec::mdot( TACSVec ** x, TacsScalar * ans, int m ){
  for ( int k = 0; k < m; k++ ){
    ans[k] = dot(x[k]);
  }
}

const char * TACSMat::TACSObjectName(){ return matName; }
const char * TACSMat::matName = "TACSMat";

const char * TACSPc::TACSObjectName(){ return pcName; }
const char * TACSPc::pcName = "TACSPc";

const char * KSMPrint::TACSObjectName(){ return printName; }
const char * KSMPrint::printName = "TACSPrint";

const char * TACSKsm::TACSObjectName(){ return ksmName; }
const char * TACSKsm::ksmName = "TACSKsm";

/*
  The KSMPrint object controls how the data it output to the screen
  (if at all) 

  input:
  descript: prepend this descriptor to the output
  rank:     print out the result only if rank == 0
  freq:     print out the result with this frequency
*/
KSMPrintStdout::KSMPrintStdout( const char * _descript, int _rank, int _freq ){
  rank = _rank;
  freq = _freq;
  if (freq < 1){ freq = 1; }

  // Copy the description string to a local array
  size_t n = strlen(_descript);
  descript = new char[ n+1 ];
  strcpy(descript, _descript);
}

KSMPrintStdout::~KSMPrintStdout(){
  delete [] descript;
}

/*
  Print the residual norm for the given iteration on the root processor

  input:
  iter:  the iteration count
  res:   the residual normal
*/
void KSMPrintStdout::printResidual( int iter, TacsScalar res ){
  if ((iter == 0 || iter % freq == 0) && rank == 0){
    printf("%s[%3d]: %15.8e\n", descript, iter, RealPart(res));
  }
}

/*
  Print out the information string only on the root processor

  input:
  cstr: the information string to print out
*/
void KSMPrintStdout::print( const char * cstr ){
  if (rank == 0){
    printf("%s", cstr);
  }
}

/*
  A KSM print object for output to a file directly

  input:
  filename: the name of the file
  descript: the description to append to any residual output
  rank:     output/create file if rank == 0
  freq:     the frequency to generate output
*/
KSMPrintFile::KSMPrintFile( const char * filename, const char * _descript, 
			    int _rank, int _freq ){
  rank = _rank;
  freq = _freq;
  if (freq <= 0){ freq = 1; }

  // Copy the description to a local array
  size_t n = strlen(_descript);
  descript = new char[ n+1 ];
  strcpy(descript, _descript);

  fp = NULL;
  if (rank == 0){
    fp = fopen(filename, "w");
  }
}

KSMPrintFile::~KSMPrintFile(){
  if (fp){ fclose(fp); }
  delete [] descript;
}

/*
  Print the residual out to the file
*/
void KSMPrintFile::printResidual( int iter, TacsScalar res ){
  if ((iter == 0 || iter % freq == 0) && rank == 0){
    if (fp){ fprintf(fp, "%s[%d]: %15.8e\n", descript, iter, RealPart(res)); }
  }
}

/*
  Print a string out to the file
*/
void KSMPrintFile::print( const char * cstr ){
  if (fp){ fprintf(fp, "%s", cstr); }
}

/*
  The preconditioned conjugate gradient method

  This object is used to perform the preconditioned conjugate gradient
  method for a linear system. The preconditioner cannot be flexible.

  input:
  mat:    the matrix operator
  pc:     the preconditioner operator
  reset:  reset the CG iterations every 'reset' iterations
  nouter: the number of resets to try before giving up
*/
PCG::PCG( TACSMat * _mat, TACSPc * _pc, int _reset, int _nouter ){
  monitor = NULL;

  mat = _mat;
  pc = _pc;
  
  mat->incref();
  pc->incref();

  reset = _reset;
  nouter = _nouter;

  // Set default absolute and relative tolerances
  rtol = 1e-8;
  atol = 1e-30;

  // Create the vectors required 
  work = mat->createVec();
  R = mat->createVec();
  Z = mat->createVec();
  P = mat->createVec();

  work->incref();
  R->incref();
  Z->incref();
  P->incref();  
}

PCG::~PCG(){
  mat->decref();
  pc->decref();

  work->decref();
  R->decref();
  Z->decref();
  P->decref();  
}

/*
  Set the operators for the preconditioned conjugate gradient method.
*/
void PCG::setOperators( TACSMat * _mat, TACSPc * _pc ){
  if (_mat){
    _mat->incref();
    if (mat){ mat->decref(); }
    mat = _mat;
  }
  if (_pc){
    _pc->incref();  
    if (pc){ pc->decref(); }
    pc = _pc;
  }
}

void PCG::getOperators( TACSMat ** _mat, TACSPc ** _pc ){
  *_mat = mat;
  *_pc = pc;
}

void PCG::setTolerances( double _rtol, double _atol ){
  rtol = _rtol;
  atol = _atol;
}

void PCG::setMonitor( KSMPrint * _monitor ){
  _monitor->incref();
  if ( monitor ){
    monitor->decref();
  }
  monitor = _monitor;
}


/*
  Solve the linear system with the preconditioned conjugate gradient
  method 

  input:
  b:          the right-hand-side
  x:          the solution vector
  zero_guess: flag to indicate whether to start with x = 0
*/
void PCG::solve( TACSVec * b, TACSVec * x, int zero_guess ){
  int solve_flag = 0;
  TacsScalar rhs_norm = 0.0;
  // R, Z and P are work-vectors
  // R == the residual

  for ( int count = 0; count < nouter; count++ ){
    // R is the residual
    if (zero_guess && count == 0){
      // If the initial guess is zero
      x->zeroEntries();  // Set x = 0
      R->copyValues(b);  // R = b
    }
    else {
      mat->mult(x, R);        // R = A*x
      R->axpby(1.0, -1.0, b); // R = b - A*x
    }    
    
    if (count == 0){
      rhs_norm = R->norm();
    }

    if (monitor && count == 0){ 
      monitor->printResidual(0, rhs_norm);
    }

    if ( rhs_norm > atol ){
      // Apply the preconditioner
      pc->applyFactor(R, Z);
      
      // Copy Z to P, ie. P = Z
      P->copyValues(Z);
      
      for ( int i = 0; i < reset; i++ ){
	mat->mult(P, work);                     // work = A*P
	TacsScalar temp = R->dot(Z);            // (R,Z)
	TacsScalar alpha = temp/(work->dot(P)); // alpha = (R,Z)/(A*P,P)
	x->axpy(alpha, P);                      // x = x + alpha * P
	R->axpy(-alpha, work);                  // R' = R - alpha * A*P
	pc->applyFactor(R, Z);                  // Z' = M^{-1} R
	TacsScalar beta = R->dot(Z)/temp;       // beta = (R',Z')/(R,Z)
        P->axpby(1.0, beta, Z);                 // P' = Z' + beta * P 

	TacsScalar norm = R->norm();
	
	if (monitor){
	  monitor->printResidual(i+1, norm);
	}
	
	if (norm < atol || norm < rhs_norm*rtol){
	  solve_flag = 1;
	  break;
	}
      }
    }

    if ( solve_flag ){
      break;
    }
  }
}

/*
  Classical Gram-Schmidt orthogonalization.
  
  Given an input vector q, and the set of orthogonal vectors w. Make q
  orthogonal to the set of vectors w.
  
  Return an array h = w^{T} q.
  q' = q - w h
  w^{T} q' = w^{T} q - w^{T} w w^{T} q 
*/
static void ClassicalGramSchmidt( TacsScalar * h, TACSVec * q, 
				  TACSVec ** w, int nvecs ){
  q->mdot(w, h, nvecs);
  for ( int j = 0; j < nvecs; j++ ){
    q->axpy(-h[j], w[j]);
  }
}

/*
  Modified Gram-Schmidt orthogonalization
*/
static void ModifiedGramSchmidt( TacsScalar * h, TACSVec * q, 
				 TACSVec ** w, int nvecs ){
  for ( int j = 0; j < nvecs; j++ ){
    h[j] = w[j]->dot(q);
    q->axpy(-h[j], w[j]);
  }
}

/*
  Create a GMRES object for solving a linear system.

  This automatically allocates the requried Krylov subspace on
  initialization.

  input:
  mat:        the matrix operator
  pc:         the preconditioner
  m:          the size of the Krylov subspace
  nrestart:   the number of restarts before we give up
  isFlexible: is the preconditioner actually flexible? If so use FGMRES
*/
GMRES::GMRES( TACSMat * _mat, TACSPc * _pc, int _m, int _nrestart, 
	      int _isFlexible ){
  init(_mat, _pc, _m, _nrestart, _isFlexible);
}

/*
  Initialize GMRES without a preconditioner. This is used for the
  ApproximateSchur preconditioner.

  input:
  mat:      the matrix operator
  m:        the size of the Krylov subspace
  nrestart: try this many times before giving up
*/
GMRES::GMRES( TACSMat * _mat, int _m, int _nrestart ){
  init(_mat, NULL, _m, _nrestart, 0);
}

/*
  Initialize the underlying GMRES data structure.

  This is called by both of the two constructors above. 
*/
void GMRES::init( TACSMat * _mat, TACSPc * _pc, 
		  int _m, int _nrestart, int _isFlexible ){

  orthogonalize = ModifiedGramSchmidt;
  monitor = NULL;
  monitor_time = 0;
  msub = _m;
  nrestart = (_nrestart >= 0 ? _nrestart : 0);
  isFlexible = _isFlexible;

  mat = _mat;
  pc = _pc;
  mat->incref();

  if (pc){
    pc->incref();
  }
  else {
    // If there's no Pc then it doesn't have to be a flexible variant
    isFlexible = 0; 
  }

  // Set default absolute and relative tolerances
  rtol = 1e-8;
  atol = 1e-30;

  // Allocate the subspace of vectors
  work = NULL;  // Required only for regular GMRES
  Z = NULL;     // Required only for flexible GMRES
  W = new TACSVec*[ msub+1 ];

  for ( int i = 0; i < msub+1; i++ ){
    W[i] = mat->createVec();
    W[i]->incref();
  }

  if ( isFlexible ){
    Z = new TACSVec*[ msub ];

    for ( int i = 0; i < msub; i++ ){
      Z[i] = mat->createVec();
      Z[i]->incref();
    }
  }
  else if ( pc ){
    // Allocate the work array
    work = mat->createVec();
    work->incref();
  }

  // Allocate space for the Hessenberg matrix
  // This is a (msub+1) x msub matrix with non-zeros 
  // on and above the first diagonal below the main diagonal
  Hptr = new int[msub + 1];
  Hptr[0] = 0;

  for ( int i = 0; i < msub; i++ ){
    Hptr[i+1] = Hptr[i] + i+2;
  }
  
  int size = Hptr[msub];
  H = new TacsScalar[size];        // The Hessenberg matrix
  res = new TacsScalar[ msub+1 ];  // The residual

  memset(H, 0, size*sizeof(TacsScalar));
  memset(res, 0, (msub+1)*sizeof(TacsScalar));

  // Allocate the terms that represent the unitary Q matrix
  // in the QR factorixation of H
  Qsin = new TacsScalar[ msub ];
  Qcos = new TacsScalar[ msub ];

  memset(Qsin, 0, msub*sizeof(TacsScalar));
  memset(Qcos, 0, msub*sizeof(TacsScalar));
}

/*
  Free the data/memory allocated by GMRES
*/
GMRES::~GMRES(){
  mat->decref();

  if (pc){   pc->decref(); }
  if (work){ work->decref(); }
  
  if (Z){
    for ( int i = 0; i < msub; i++ ){
      Z[i]->decref();
    }
    delete [] Z;
  }

  for ( int i = 0; i < msub+1; i++ ){
    W[i]->decref();
  }
  delete [] W;

  if (monitor){ monitor->decref(); }

  delete [] H;
  delete [] Hptr;
  delete [] res;
  delete [] Qsin;
  delete [] Qcos;
}

/*
  Set the matrix and/or preconditioner operators. If either (or both)
  are NULL, then the operator is not replaced.

  input:
  mat: the new matrix operator (possibly NULL)
  pc:  the new preconditioner operator (possibly NULL)
*/
void GMRES::setOperators( TACSMat * _mat, TACSPc * _pc ){
  if (_mat){
    _mat->incref();
    if (mat){ mat->decref(); }
    mat = _mat;
  }
  if (_pc){
    _pc->incref();  
    if (pc){ pc->decref(); }
    pc = _pc;
  }
}

/*
  Retrieve the preconditioner operators from GMRES

  output:
  mat: the matrix
  pc:  the preconditioner
*/
void GMRES::getOperators( TACSMat ** _mat, TACSPc ** _pc ){
  *_mat = mat;
  *_pc = pc;
}

/*
  Set the relative and absolute tolerances used for the stopping
  criterion.

  input:
  rtol: the relative tolerance ||r_k|| < rtol*||r_0||
  atol: the absolute tolerancne ||r_k|| < atol
*/
void GMRES::setTolerances( double _rtol, double _atol ){
  rtol = _rtol;
  atol = _atol;
}

/*
  Set the object to control how the convergence history is displayed
  (if at all)

  input:
  monitor: the KSMPrint monitor object
*/
void GMRES::setMonitor( KSMPrint * _monitor ){
  _monitor->incref();
  if (monitor){
    monitor->decref();
  }
  monitor = _monitor;
}

/*
  Set the type of orthogonalization to use. This will be either
  classical Gram-Schmidt or modified Gram-Schmidt (default).
  
  Unless you have a good reason, you should use modified Gram-Schmidt.
*/
void GMRES::setOrthoType( enum OrthoType otype ){
  if (otype == CLASSICAL_GRAM_SCHMIDT){
    orthogonalize = ClassicalGramSchmidt;
  }
  else {
    orthogonalize = ModifiedGramSchmidt;
  }
}

/*
  Set a flag to also monitor the time spent in various operations
  internally. 
*/
void GMRES::setTimeMonitor(){
  monitor_time = 1;
}

const char * GMRES::TACSObjectName(){
  return gmresName;
}

const char * GMRES::gmresName = "GMRES";

/*
  Try to solve the linear system using GMRES.

  The following code tries to solve the linear system using GMRES (or
  FGMRES if the preconditioner is flexible.) 

  input:
  b:          the right-hand-side
  x:          the solution vector (with possibly significant entries)
  zero_guess: flag to indicate whether to zero entries of x before solution
*/
void GMRES::solve( TACSVec * b, TACSVec * x, int zero_guess ){
  TacsScalar rhs_norm = 0.0;
  int solve_flag = 0; 

  double t_pc = 0.0, t_ortho = 0.0;
  double t_total = 0.0;

  if (monitor_time){ t_total = MPI_Wtime(); }

  for ( int count = 0 ; count < nrestart+1; count++ ){
    // Compute the residual
    if (zero_guess && count == 0){
      // If the initial guess is zero
      x->zeroEntries();    // Set x = 0
      W[0]->copyValues(b); // W[0] = b

      res[0] = W[0]->norm();
      W[0]->scale(1.0/res[0]); // W[0] = b/|| b ||
    }
    else {
      // If the initial guess is non-zero or restarting
      mat->mult(x, W[0]);
      W[0]->axpy(-1.0, b); // W[0] = A*x - b

      res[0] = W[0]->norm();
      W[0]->scale(-1.0/res[0]); // W[0] = (b - A*x)/|| b - A*x ||    
    }

    if (monitor){
      monitor->printResidual(0, fabs(res[0]));
    }

    if (count == 0){
      rhs_norm = res[0]; // The initial residual 
    }  

    int niters = 0; // Keep track of the size of the Hessenberg matrix  

    if (res[0] < atol){
      break;
    }

    for ( int i = 0; i < msub; i++ ){
      if (monitor_time){ t_pc -= MPI_Wtime(); }
      if (isFlexible){
        // Apply the preconditioner, Z[i] = M^{-1} W[i] 
	pc->applyFactor(W[i], Z[i]); 
	mat->mult(Z[i], W[i+1]);  // W[i+1] = A * Z[i] = A * M^{-1} * W[i]
      }
      else {
	if (pc){
          // Apply the preconditioner, work = M^{-1} W[i] 
	  pc->applyFactor(W[i], work); 
	  mat->mult(work, W[i+1]); // W[i+1] = A * work = A * M^{-1} * W[i]
	}
	else {
	  mat->mult(W[i], W[i+1]); // Compute W[i+1] = A * W[i] 
	}
      }
      if (monitor_time){ 
	double t0 = MPI_Wtime();
	t_pc += t0; t_ortho -= t0; 
      }

      // Build the orthogonal factorization MGS
      orthogonalize(&H[Hptr[i]], W[i+1], W, i+1);
      if (monitor_time){ t_ortho += MPI_Wtime(); }

      H[i+1 + Hptr[i]] = W[i+1]->norm(); // H[i+1,i] = || W[i+1] ||
      W[i+1]->scale(1.0/H[i+1 + Hptr[i]]); // W[i+1] = W[i+1]/|| W[i+1] ||
      
      // Apply the existing part of Q to the new components of 
      // the Hessenberg matrix
      TacsScalar h1, h2;
      for ( int k = 0; k < i; k++ ){
	h1 = H[k   + Hptr[i]];
	h2 = H[k+1 + Hptr[i]];
	H[k   + Hptr[i]] =  h1*Qcos[k] + h2*Qsin[k];
	H[k+1 + Hptr[i]] = -h1*Qsin[k] + h2*Qcos[k];
      }
      
      // Now, compute the rotation for the new column that was just added
      h1 = H[i   + Hptr[i]];
      h2 = H[i+1 + Hptr[i]];
      TacsScalar sq = sqrt(h1*h1 + h2*h2);
      
      Qcos[i] = h1/sq;
      Qsin[i] = h2/sq;
      H[i   + Hptr[i]] =  h1*Qcos[i] + h2*Qsin[i];
      H[i+1 + Hptr[i]] = -h1*Qsin[i] + h2*Qcos[i];
      
      // Update the residual
      h1 = res[i];
      // h2 = res[i+1]; = 0
      res[i]   =   h1*Qcos[i];
      res[i+1] = - h1*Qsin[i];

      if (monitor){
	monitor->printResidual(i+1, fabs(res[i+1]));
      }
      
      niters++;
      
      if (fabs(res[i+1]) < atol ||
          fabs(res[i+1]) < rtol*rhs_norm){
	// Set the solve flag 
	solve_flag = 1;
	
	break;
      }
    }

    // Now, compute the solution - the linear combination of the
    // Arnoldi vectors. H is upper triangular

    // Compute the weights
    for ( int i = niters-1; i >= 0; i-- ){
      for ( int j = i+1; j < niters; j++ ){ // 
	res[i] = res[i] - H[i + Hptr[j]]*res[j];
      }
      res[i] = res[i]/H[i + Hptr[i]];
    }
    
    // Compute the linear combination
    if (isFlexible){ // Flexible variant
      for ( int i = 0; i < niters; i++ ){
	x->axpy(res[i], Z[i]);
      }
    }
    else if (!pc){   // If there's no pc
      for ( int i = 0; i < niters; i++ ){
	x->axpy(res[i], W[i]);
      }
    }
    else {             // If the pc isn't flexible
      work->zeroEntries();
      for ( int i = 0; i < niters; i++ ){
	work->axpy(res[i], W[i]);
      }
    
      // Apply M^{-1} to the linear combination
      pc->applyFactor(work, W[0]);
      x->axpy(1.0, W[0]);
    }
    
    if (solve_flag){
      break;
    }
  }
  
  if (monitor_time && monitor){
    t_total = MPI_Wtime() - t_total;
    char str_mat[80], str_ort[80], str_tot[80];
    sprintf(str_mat, "pc-mat time %10.6f\n", t_pc);
    sprintf(str_ort, "ortho time  %10.6f\n", t_ortho);
    sprintf(str_tot, "total time  %10.6f\n", t_total);
    monitor->print(str_mat);
    monitor->print(str_ort);
    monitor->print(str_tot);
  }
}

/*
  Create the GCROT linear system solver

  This constructor creates an object that uses a simplified variant of
  GCROT described by Hicken and Zingg.

  input:
  mat:        the matrix operator
  pc:         the preconditioner
  outer:      the number of outer vectors
  max_outer:  the maximum number of outer iterations before we give up
  msub:       the size of the underlying GMRES (FGMRES) subspace
  isFlexible: flag to indicate required use of flexible GCROT
*/
GCROT::GCROT( TACSMat * _mat, TACSPc * _pc, int _outer, int _max_outer, 
	      int _msub, int _isFlexible ){
  init(_mat, _pc, _outer, _max_outer, _msub, _isFlexible);
}

/*
  Create the GCROT linear system solver without a preconditioner. Note
  that this variant is never flexible.
*/
GCROT::GCROT( TACSMat * _mat, int _outer, int _max_outer, int _msub ){
  init(_mat, NULL, _outer, _max_outer, _msub, 0);
}

/*
  Initialize the GCROT data.

  This allocates all the vectors required for the underlying subspace.
  After initialization, no new significant memory allocation is
  required.  
*/
void GCROT::init( TACSMat * _mat, TACSPc * _pc, int _outer, int _max_outer, 
		  int _msub, int _isFlexible ){
  monitor = NULL;
  msub = _msub;             // Size of the F/GMRES subspace
  outer = _outer;           // Number of outer iterations
  max_outer = _max_outer;   // Maximum number of outer iterations
  isFlexible = _isFlexible; // Flexible variant?

  mat = _mat;
  pc = _pc;
  mat->incref();

  if (pc){
    pc->incref();
  }
  else {
    // If there's no Pc then it doesn't have to be a flexible variant
    isFlexible = 0; 
  }

  // Set default absolute and relative tolerances
  rtol = 1e-8;
  atol = 1e-30;

  // Allocate the subspace of vectors
  Z = NULL;     // Required only for flexible GMRES
  W = new TACSVec*[ msub+1 ];

  for ( int i = 0; i < msub+1; i++ ){
    W[i] = mat->createVec();
    W[i]->incref();
  }

  if (isFlexible){
    Z = new TACSVec*[ msub ];

    for ( int i = 0; i < msub; i++ ){
      Z[i] = mat->createVec();
      Z[i]->incref();
    }
  }

  R = mat->createVec();  
  c_hat = mat->createVec();
  u_hat = mat->createVec();

  R->incref();
  c_hat->incref();
  u_hat->incref();

  U = new TACSVec*[ outer ];
  C = new TACSVec*[ outer ];

  for ( int i = 0; i < outer; i++ ){
    U[i] = mat->createVec();
    U[i]->incref();
    C[i] = mat->createVec();
    C[i]->incref();
  }

  // Allocate space for the Hessenberg matrix
  // This is a (msub+1) x msub matrix with non-zeros 
  // on and above the first diagonal below the main diagonal
  Hptr = new int[msub + 1];
  Hptr[0] = 0;

  for ( int i = 0; i < msub; i++ ){
    Hptr[i+1] = Hptr[i] + i+2;
  }
  
  int size = Hptr[msub];
  H   = new TacsScalar[size];        // The Hessenberg matrix
  res = new TacsScalar[ msub+1 ];  // The residual
  B   = new TacsScalar[ msub * outer ];

  memset(H, 0, size*sizeof(TacsScalar));
  memset(res, 0, (msub+1)*sizeof(TacsScalar));
  memset(B, 0, msub*outer*sizeof(TacsScalar));

  // Allocate the terms that represent the unitary Q matrix
  // in the QR factorixation of H
  Qsin = new TacsScalar[ msub ];
  Qcos = new TacsScalar[ msub ];

  memset(Qsin, 0, msub*sizeof(TacsScalar));
  memset(Qcos, 0, msub*sizeof(TacsScalar));
}

/*
  Delete the object and free all the data
*/
GCROT::~GCROT(){
  mat->decref();

  if (pc){ pc->decref(); }

  R->decref();
  u_hat->decref();
  c_hat->decref();

  for ( int i = 0; i < outer; i++ ){
    if (U[i]){ U[i]->decref(); }
    if (C[i]){ C[i]->decref(); }
  }

  delete [] U;
  delete [] C;

  if (Z){
    for ( int i = 0; i < msub; i++ ){
      Z[i]->decref();
    }
    delete [] Z;
  }

  for ( int i = 0; i < msub+1; i++ ){
    W[i]->decref();
  }
  delete [] W;

  if (monitor){ monitor->decref(); }

  delete [] H;
  delete [] Hptr;
  delete [] res;
  delete [] B;
  delete [] Qsin;
  delete [] Qcos;
}

/*
  Set the matrix/preconditioner operators used for GCROT
*/ 
void GCROT::setOperators( TACSMat * _mat, TACSPc * _pc ){
  if (_mat){
    _mat->incref();
    if (mat){ mat->decref(); }
    mat = _mat;
  }
  if (_pc){
    _pc->incref();  
    if (pc){ pc->decref(); }
    pc = _pc;
  }
}

/*
  Retrieve the matrix/preconditioner operators set in the GCROT object
*/
void GCROT::getOperators( TACSMat ** _mat, TACSPc ** _pc ){
  *_mat = mat;
  *_pc = pc;
}

/*
  Set the relative and absolute convergence tolerances for GCROT
*/
void GCROT::setTolerances( double _rtol, double _atol ){
  rtol = _rtol;
  atol = _atol;
}

/*
  Set the residual/solution monitor object
*/
void GCROT::setMonitor( KSMPrint * _monitor ){
  _monitor->incref();
  if (monitor){
    monitor->decref();
  }
  monitor = _monitor;
}

const char * GCROT::TACSObjectName(){
  return gcrotName;
}

const char * GCROT::gcrotName = "GCROT";

/*
  Solve the linear system defined by the matrix operator with the
  given (possibly flexible) preconditioner 

  input:
  b:          the input right-hand-side
  x:          the solution vector
  zero_guess: flag to treat x as an initial guess or zero
*/
void GCROT::solve( TACSVec * b, TACSVec * x, int zero_guess ){
  TacsScalar rhs_norm = 0.0;
  int solve_flag = 0; 
  int mat_iters = 0;

  // Compute the residual
  if (zero_guess){
    // If the initial guess is zero
    x->zeroEntries();    // Set x = 0
    R->copyValues(b);  // R = b
  }
  else {
    // If the initial guess is non-zero or restarting
    mat->mult(x, R);      // R = A*x
    mat_iters++;
    R->axpy(-1.0, b);     // R = A*x - b
    R->scale(-1.0);       // R = b - A*x0
  }

  rhs_norm = R->norm(); // The initial residual 

  if (rhs_norm < atol){
    return;
  }

  for ( int count = 0 ; count < max_outer; count++ ){
    // Size of the U/C subspaces
    int outer_size = ( count < outer ? count : outer ); // min( count, outer ); 
    int niters = 0; // The size of the Hessenberg matrix  

    W[0]->copyValues(R); // W[0] = R
    res[0] = W[0]->norm();
    W[0]->scale(1.0/res[0]); // W[0] = b/|| b ||

    if ( monitor ){
      monitor->printResidual(mat_iters, fabs(res[0]));
    }

    // The inner F/GMRES loop
    for ( int i = 0; i < msub; i++ ){
      if (isFlexible){
        // Apply the preconditioner, Z[i] = M^{-1} W[i] 
	pc->applyFactor(W[i], Z[i]);
	mat->mult( Z[i], W[i+1] ); // W[i+1] = A * Z[i] = A * M^{-1} * W[i]
      }
      else {
	if (pc){
	  // Use u_hat here as a temporary array
          // Apply the preconditioner, work = M^{-1} W[i] 
	  pc->applyFactor(W[i], u_hat); 
	  mat->mult(u_hat, W[i+1]); // W[i+1] = A * work = A * M^{-1} * W[i]
	}
	else {
	  mat->mult(W[i], W[i+1]); // Compute W[i+1] = A * W[i] 
	}
      }
      mat_iters++;

      // First, orthonormalize W[i+1] w.r.t. 
      // the basis C[j] j = 0 .. outer_size-1 B is (outer_size, msub)
      for ( int j = 0; j < outer_size; j++ ){
	B[j*msub + i] = W[i+1]->dot(C[j]); // B[j,i] = dot( W[i+1], C[j] )  
	W[i+1]->axpy(-B[j*msub + i], C[j]); // W[i+1] = W[i+1] - B[j,i] * C[j]
      }
      
      // Build the orthogonal factorization MGS
      for ( int j = i; j >= 0; j-- ){
	H[j + Hptr[i]] = W[i+1]->dot(W[j]); // H[j,i] = dot( W[i+1], W[i] )
	W[i+1]->axpy(-H[j + Hptr[i]], W[j]); // W[i+1] = W[i+1] - H[j,i] * W[j]
      }
      
      H[i+1 + Hptr[i]] = W[i+1]->norm(); // H[i+1,i] = || W[i+1] ||
      W[i+1]->scale(1.0/H[i+1 + Hptr[i]]); // W[i+1] = W[i+1]/|| W[i+1] ||
       
      // Apply the existing part of Q to the new components of 
      // the Hessenberg matrix
      TacsScalar h1, h2;
      for ( int k = 0; k < i; k++ ){
	h1 = H[k   + Hptr[i]];
	h2 = H[k+1 + Hptr[i]];
	H[k   + Hptr[i]] =  h1*Qcos[k] + h2*Qsin[k];
	H[k+1 + Hptr[i]] = -h1*Qsin[k] + h2*Qcos[k];
      }
      
      // Now, compute the rotation for the new column that was just added
      h1 = H[i   + Hptr[i]];
      h2 = H[i+1 + Hptr[i]];
      TacsScalar sq = sqrt(h1*h1 + h2*h2);
      
      Qcos[i] = h1/sq;
      Qsin[i] = h2/sq;
      H[i   + Hptr[i]] =  h1*Qcos[i] + h2*Qsin[i];
      H[i+1 + Hptr[i]] = -h1*Qsin[i] + h2*Qcos[i];
      
      // Update the residual
      h1 = res[i];
      // h2 = res[i+1]; = 0
      res[i]   =   h1*Qcos[i];
      res[i+1] = - h1*Qsin[i];

      if ( monitor ){
	monitor->printResidual(mat_iters, fabs(res[i+1]));
      }
      
      niters++;
      
      if (fabs(res[i+1]) < atol ||
          fabs(res[i+1]) < rtol*rhs_norm){
	// Set the solve flag 
	solve_flag = 1;
	
	break;
      }
    }

    // Now, compute the solution - the linear combination of the
    // Arnoldi vectors
    // H is now upper triangular

    // Compute the weights
    for ( int i = niters-1; i >= 0; i-- ){
      // Compute res[i] = res[i] - H[i,j] * res[j];
      for ( int j = i+1; j < niters; j++ ){
	res[i] = res[i] - H[i + Hptr[j]] * res[j];
      }
      res[i] = res[i]/H[i + Hptr[i]];
    }

    // u_hat = Z * y - U * B * y
    // where y = the weights computed from the inner GMRES loop above

    // Compute the linear combination
    if ( isFlexible ){ // Flexible variant
      u_hat->zeroEntries();
      for ( int i = 0; i < niters; i++ ){
	u_hat->axpy(res[i], Z[i]);
      }
    }
    else if ( !pc ){   // If there's no pc
      u_hat->zeroEntries();
      for ( int i = 0; i < niters; i++ ){
	u_hat->axpy(res[i], W[i]);
      }
    }
    else {             // If the pc isn't flexible
      // Use c_hat here as a temporary array
      c_hat->zeroEntries();
      for ( int i = 0; i < niters; i++ ){
	c_hat->axpy(res[i], W[i]);
      }
    
      // Apply u_hat = M^{-1} c_hat to the linear combination
      pc->applyFactor(c_hat, u_hat);
    }

    // Now complete u_hat by computing u_hat = u_hat - U * B * y
    // Here y = res
    // Note that B is (outer_size,msub)
    for ( int i = 0; i < outer_size; i++ ){
      TacsScalar bsum = TacsScalar(0.0);
      for ( int j = 0; j < niters; j++ ){
	bsum = bsum + B[i*msub + j] * res[j];
      }

      u_hat->axpy(-bsum, U[i]);
    }


    c_hat->zeroEntries();
    
    // Compute c_hat = W * H * y
    // The Hessenberg matrix has been modified by pre-multiplying by Q
    // now, undo these rotations and multiply by Q^{T}
    for ( int i = 0; i < niters; i++ ){
      TacsScalar h1, h2;
      for ( int k = i; k >= 0; k-- ){
	h1 = H[k   + Hptr[i]];
	h2 = H[k+1 + Hptr[i]];
	H[k   + Hptr[i]] =  h1*Qcos[k] - h2*Qsin[k];
	H[k+1 + Hptr[i]] =  h1*Qsin[k] + h2*Qcos[k];
      }
    }

    // The matrix H is (niters+1) by (niters) // H[i,j] = H[i + Hptr[j]]
    // y = the variables res is niters
    for ( int i = 0; i < niters+1; i++ ){
      TacsScalar hsum = TacsScalar(0.0);
      int j = i-1;
      if ( i == 0 ){ j = 0; }

      // hsum = hsum + H[i,j] * res[j]
      for ( ; j < niters; j++ ){
	hsum += H[i + Hptr[j]] * res[j];
      }

      c_hat->axpy(hsum, W[i]);
    }

    TacsScalar cnorm = c_hat->norm();
    c_hat->scale(1.0/cnorm);
    u_hat->scale(1.0/cnorm);

    // Update the residual and solution
    TacsScalar alpha = R->dot(c_hat);
    R->axpy(-alpha, c_hat);
    x->axpy(alpha, u_hat); 

    // Compute the residual and compare... 
    // mat->mult( x, R );      // R = A*x
    // R->axpy( -1.0, b );     // R = A*x - b
    // R->scale( -1.0 );

    if ( solve_flag ){
      break;
    }

    // Now, prepend c_hat and u_hat to C and U respectively.
    // All this is done with pointers - no copying actually occurs.
    // Otherwise, discard the oldest vector - 
    // but use it for the 'new' u_hat/c_hat
    TACSVec * u = u_hat;
    TACSVec * c = c_hat;
    
    // This is the size of the subspace on the next iteration
    outer_size = ( count+1 < outer ? count+1 : outer );

    for ( int i = 0; i < outer_size; i++ ){
      TACSVec * tu = U[i];
      TACSVec * tc = C[i];
      
      U[i] = u;
      C[i] = c;	
      
      u = tu;
      c = tc;
    }
      
    u_hat = u;
    c_hat = c;          
  }
}

/*
  Constrainted GMRES code.

  Perform GMRES but also satisfy a set of constraints at the same time.

  These constraints may be dense and are added as a set of vectors at
  some point after the creation of the object. This code could
  probably be incorperated into the regular GMRES code for more
  flexibility later, but for now will just be included separately.
*/
ConGMRES::ConGMRES( TACSMat * _mat, TACSPc * _pc, int _m, int _nrestart, 
		    int _isFlexible ){
  init(_mat, _pc, _m, _nrestart, _isFlexible);
}

ConGMRES::ConGMRES( TACSMat * _mat, int _m, int _nrestart ){
  init(_mat, NULL, _m, _nrestart, 0);
}

void ConGMRES::init( TACSMat * _mat, TACSPc * _pc, int _m, int _nrestart, 
		     int _isFlexible ){

  orthogonalize = ModifiedGramSchmidt;
  monitor = NULL;
  monitor_time = 0;
  msub = _m;
  nrestart = (_nrestart >= 0 ? _nrestart : 0);
  isFlexible = _isFlexible;

  mat = _mat;
  pc = _pc;
  mat->incref();

  if (pc){
    pc->incref();
  }
  else {
    // If there's no Pc then it doesn't have to be a flexible variant
    isFlexible = 0;
  }

  // Set default absolute and relative tolerances
  rtol = 1e-8;
  atol = 1e-30;

  // Allocate the subspace of vectors
  work = NULL;   // Required only for regular GMRES
  Z = NULL;      // Required only for flexible GMRES
  W = new TACSVec*[ msub+1 ]; // The orthonormal Krylov subspace
  C = NULL;      // The constraint vectors
  nconstr = 0;   // The number of constraints
  lambda = NULL; // The Lagrange multipliers
  Wl = NULL;     // The basis for lambda

  for ( int i = 0; i < msub+1; i++ ){
    W[i] = mat->createVec();
    W[i]->incref();
  }

  if ( isFlexible ){
    Z = new TACSVec*[ msub ];

    for ( int i = 0; i < msub; i++ ){
      Z[i] = mat->createVec();
      Z[i]->incref();
    }
  }
  else if ( pc ){
    // Allocate the work array
    work = mat->createVec();
    work->incref();
  }

  // Allocate space for the Hessenberg matrix
  // This is a (msub+1) x msub matrix with non-zeros 
  // on and above the first diagonal below the main diagonal
  Hptr = new int[msub + 1];
  Hptr[0] = 0;

  for ( int i = 0; i < msub; i++ ){
    Hptr[i+1] = Hptr[i] + i+2;
  }
  
  int size = Hptr[msub];
  H = new TacsScalar[size];        // The Hessenberg matrix
  res = new TacsScalar[ msub+1 ];  // The residual

  memset(H, 0, size*sizeof(TacsScalar));
  memset(res, 0, (msub+1)*sizeof(TacsScalar));

  // Allocate the terms that represent the unitary Q matrix
  // in the QR factorixation of H
  Qsin = new TacsScalar[ msub ];
  Qcos = new TacsScalar[ msub ];

  memset(Qsin, 0, msub*sizeof(TacsScalar));
  memset(Qcos, 0, msub*sizeof(TacsScalar));
}

ConGMRES::~ConGMRES(){
  mat->decref();

  if (pc){   pc->decref(); }
  if (work){ work->decref(); }
  
  if (Z){
    for ( int i = 0; i < msub; i++ ){
      Z[i]->decref();
    }
    delete [] Z;
  }

  for ( int i = 0; i < msub+1; i++ ){
    W[i]->decref();
  }
  delete [] W;

  if (C){
    for ( int i = 0; i < nconstr; i++ ){
      C[i]->decref();
    }
    delete [] C;
  }

  if (Wl){ delete [] Wl; }
  if (lambda){ delete [] lambda; }
  if (monitor){ monitor->decref(); }

  delete [] H;
  delete [] Hptr;
  delete [] res;
  delete [] Qsin;
  delete [] Qcos;
}

void ConGMRES::setOperators( TACSMat * _mat, TACSPc * _pc ){
  if (_mat){
    _mat->incref();
    if (mat){ mat->decref(); }
    mat = _mat;
  }
  if (_pc){
    _pc->incref();
    if (pc){ pc->decref(); }
    pc = _pc;
  }
}

void ConGMRES::getOperators( TACSMat ** _mat, TACSPc ** _pc ){
  *_mat = mat;
  *_pc = pc;
}

void ConGMRES::setTolerances( double _rtol, double _atol ){
  rtol = _rtol;
  atol = _atol;
}

void ConGMRES::setMonitor( KSMPrint * _monitor ){
  _monitor->incref();
  if ( monitor ){
    monitor->decref();
  }
  monitor = _monitor;
}

void ConGMRES::setOrthoType( enum OrthoType otype ){
  if ( otype == CLASSICAL_GRAM_SCHMIDT ){
    orthogonalize = ClassicalGramSchmidt;
  }
  else {
    orthogonalize = ModifiedGramSchmidt;
  }
}

void ConGMRES::setTimeMonitor(){
  monitor_time = 1;
}

const char * ConGMRES::TACSObjectName(){
  return gmresName;
}

const char * ConGMRES::gmresName = "ConGMRES";


/*!
  Add the constraints to the KSM object.

  These constraints are immediately othonormalized and should not be 
  accessed again unless this object is later destroyed or new constraints
  are set.
*/

void ConGMRES::setConstraints( TACSVec ** _C, int _nconstr ){
  // Read in the constraints and set them to be equal 
  if (C){
    for ( int i = 0; i < nconstr; i++ ){
      C[i]->decref();
    }

    delete [] C;
  }

  nconstr = _nconstr;
  C = new TACSVec*[ nconstr ];

  for ( int i = 0; i < nconstr; i++ ){
    C[i] = _C[i];
    C[i]->incref();
  }

  if (lambda){ delete [] lambda; }
  lambda = new TacsScalar[ nconstr ];
  if (Wl){ delete [] Wl; }
  Wl = new TacsScalar[ nconstr*(msub+1) ];
}

/*
  Perform vector operations but add an additional constraint to the
  linear system. This can be used to impose additional condditions.
*/

/*
  Perform the matrix-vector multiplication
  out        = A*x + C*y
  out_lambda = C*x
*/
void ConGMRES::con_mat_mult( TACSVec * x, TacsScalar * y, 
			     TACSVec * out, TacsScalar * out_lambda ){
  mat->mult(x, out);
  for ( int k = 0; k < nconstr; k++ ){
    out->axpy(y[k], C[k]);
    out_lambda[k] = C[k]->dot(x);
  }
}

// Compute the norm with the constraints
TacsScalar ConGMRES::con_norm( TACSVec * x, TacsScalar * y ){
  TacsScalar norm = x->norm();
  norm *= norm;
  for ( int k = 0; k < nconstr; k++ ){
    norm += y[k]*y[k];
  }

  return sqrt(norm);
}

// Compute the dot product
TacsScalar ConGMRES::con_dot( TACSVec * x, TacsScalar * y,
			      TACSVec * X, TacsScalar * Y ){
  TacsScalar dot_product = 0.0;
  for ( int k = 0; k < nconstr; k++ ){
    dot_product += y[k]*Y[k];
  }
  
  return x->dot(X) + dot_product;
}

  // Compute an axpy including constraints
void ConGMRES::con_axpy( TACSVec * x, TacsScalar * X, 
			 TacsScalar alpha,
			 TACSVec * y, TacsScalar * Y ){
  x->axpy(alpha, y);
  for ( int k = 0; k < nconstr; k++ ){
    X[k] += alpha*Y[k];
  }
}

// Compute the constrained scaling
void ConGMRES::con_scale( TACSVec * x, TacsScalar * y, TacsScalar alpha ){
  x->scale(alpha);
  for ( int k = 0; k < nconstr; k++ ){
    y[k] *= alpha;
  }
}

/*!
  Constrained GMRES iteration

  Compute (x,lambda), such that

  A x + C lambda = b
  C^{T} x = d 

  where C is a small set of vectors, lambda is a set of Lagrange
  multipliers and b is the right hand side.
  
  Approach:
  
  1. Compute the C = QR factorization of the set of vectors, where Q
  is orthonormal and R is upper triangular.  This is stored in place
  in the vectors C and the matrix R.
  
  2. Compute the right hand side r = (I - Q*Q^{T})b
  
  3. Compute the solution to the problem: (I - Q Q^{T})Ax = r

  (Optional)
  4. Compute the Lagrange multipliers lambda = R^{-1}(Q^{T}Ax - b)  
*/
void ConGMRES::solve( TACSVec * b, TACSVec * x, int zero_guess ){
  TacsScalar rhs_norm = 0.0;
  int solve_flag = 0; 

  double t_pc = 0.0, t_ortho = 0.0;
  double t_total = 0.0;

  if ( monitor_time ){ t_total = MPI_Wtime(); }

  for ( int count = 0 ; count < nrestart+1; count++ ){
    // Compute the residual
    if ( count == 0 && zero_guess ){
      // If the initial guess is zero
      x->zeroEntries();    // Set x = 0
      W[0]->copyValues(b); // W[0] = b

      memset(lambda, 0, nconstr*sizeof(TacsScalar));
      memset(&Wl[0], 0, nconstr*sizeof(TacsScalar));

      // W[0] = b/|| b ||
      res[0] = con_norm(W[0], &Wl[0]);
      con_scale(W[0], &Wl[0], 1.0/res[0]);
    }
    else {
      // If the initial guess is non-zero or restarting
      con_mat_mult(x, lambda, W[0], &Wl[0]);
      W[0]->axpy(-1.0, b); // W[0] = A*x + C*lambda - b

      // W[0] = b/|| b ||
      res[0] = con_norm(W[0], &Wl[0]);
      con_scale(W[0], &Wl[0], -1.0/res[0]);
    }

    if ( monitor ){
      monitor->printResidual(0, fabs(res[0]));
    }

    if ( count == 0 ){
      rhs_norm = res[0]; // The initial residual 
    }  

    int niters = 0; // Keep track of the size of the Hessenberg matrix  

    if ( res[0] < atol ){
      break;
    }

    for ( int i = 0; i < msub; i++ ){
      if ( monitor_time ){ t_pc -= MPI_Wtime(); }
      if ( isFlexible ){
	pc->applyFactor(W[i], Z[i]); // Apply the preconditioner
	con_mat_mult(Z[i], &Wl[i*nconstr], 
		     W[i+1], &Wl[(i+1)*nconstr]); // W[i+1] = A*Z[i]
      }
      else {
	if ( pc ){
	  // Apply the preconditioner, work = M^{-1} W[i] 
	  pc->applyFactor(W[i], work); 
	  con_mat_mult(work, &Wl[i*nconstr], 
		       W[i+1], &Wl[(i+1)*nconstr]); 
	}
	else {
	  con_mat_mult(W[i], &Wl[i*nconstr], 
		       W[i+1], &Wl[(i+1)*nconstr]); 
	}
      }
      if ( monitor_time ){ 
	double t0 = MPI_Wtime();
	t_pc += t0; t_ortho -= t0; 
      }

      // Build the orthogonal factorization using MGS
      for ( int j = 0; j < i+1; j++ ){
	H[j + Hptr[i]] = con_dot(W[i+1], &Wl[(i+1)*nconstr], 
				 W[j], &Wl[j*nconstr]);

	con_axpy(W[i+1], &Wl[(i+1)*nconstr], 
		 -H[j + Hptr[i]],
		 W[j], &Wl[j*nconstr]);
      }

      if ( monitor_time ){ t_ortho += MPI_Wtime(); }

      // compute W[i+1] = W[i+1]/|| W[i+1] ||
      H[i+1 + Hptr[i]] = con_norm(W[i+1], &Wl[(i+1)*nconstr]);
      con_scale(W[i+1], &Wl[(i+1)*nconstr], 1.0/H[i+1 + Hptr[i]]); 
      
      // Apply the existing part of Q to the new components of 
      // the Hessenberg matrix
      TacsScalar h1, h2;
      for ( int k = 0; k < i; k++ ){
	h1 = H[k   + Hptr[i]];
	h2 = H[k+1 + Hptr[i]];
	H[k   + Hptr[i]] =  h1*Qcos[k] + h2*Qsin[k];
	H[k+1 + Hptr[i]] = -h1*Qsin[k] + h2*Qcos[k];
      }
      
      // Now, compute the rotation for the new column that was just added
      h1 = H[i   + Hptr[i]];
      h2 = H[i+1 + Hptr[i]];
      TacsScalar sq = sqrt(h1*h1 + h2*h2);
      
      Qcos[i] = h1/sq;
      Qsin[i] = h2/sq;
      H[i   + Hptr[i]] =  h1*Qcos[i] + h2*Qsin[i];
      H[i+1 + Hptr[i]] = -h1*Qsin[i] + h2*Qcos[i];
      
      // Update the residual
      h1 = res[i];
      // h2 = res[i+1]; = 0
      res[i]   =   h1*Qcos[i];
      res[i+1] = - h1*Qsin[i];

      if ( monitor ){
	monitor->printResidual(i+1, fabs(res[i+1]));
      }
      
      niters++;
      
      if ( fabs(res[i+1]) < atol ||
	   fabs(res[i+1]) < rtol*rhs_norm ){
	// Set the solve flag 
	solve_flag = 1;
	
	break;
      }
    }

    // Now, compute the solution - the linear combination of the
    // Arnoldi vectors H is upper triangular

    // Compute the weights
    for ( int i = niters-1; i >= 0; i-- ){
      for ( int j = i+1; j < niters; j++ ){ 
	res[i] = res[i] - H[i + Hptr[j]] * res[j];
      }
      res[i] = res[i]/H[i + Hptr[i]];
    }
    
    // Compute the linear combination
    if ( isFlexible ){ // Flexible variant
      for ( int i = 0; i < niters; i++ ){
	con_axpy(x, lambda, 
		 res[i], Z[i], &Wl[i*nconstr]);
      }
    }
    else if ( !pc ){   // If there's no pc
      for ( int i = 0; i < niters; i++ ){
	con_axpy(x, lambda, 
		 res[i], W[i], &Wl[i*nconstr]);
      }
    }
    else {             // If the pc isn't flexible
      work->zeroEntries();
      for ( int i = 0; i < niters; i++ ){
	con_axpy(work, lambda, 
		 res[i], W[i], &Wl[i*nconstr]);
      }
    
      // Apply M^{-1} to the linear combination
      pc->applyFactor(work, W[0]);
      x->axpy(1.0, W[0]);
    }

    if ( solve_flag ){
      break;
    }
  }

  if ( monitor ){
    // If the initial guess is non-zero or restarting
    con_mat_mult(x, lambda, W[0], &Wl[0]);
    W[0]->axpy(-1.0, b); // W[0] = A*x + C*y - b
    TacsScalar norm = con_norm(W[0], &Wl[0]);
    char str_norm[80];
    sprintf(str_norm, "|A*x + C*y - b|: %15.8e\n", RealPart(norm));
    monitor->print(str_norm);
  }
  
  if ( monitor_time && monitor ){
    t_total = MPI_Wtime() - t_total;
    char str_mat[80], str_ort[80], str_tot[80];
    sprintf(str_mat, "pc-mat time %10.6f\n", t_pc);
    sprintf(str_ort, "ortho time  %10.6f\n", t_ortho);
    sprintf(str_tot, "total time  %10.6f\n", t_total);
    monitor->print(str_mat);
    monitor->print(str_ort);
    monitor->print(str_tot);
  }
}

/*
  Create the preconditioner class with the specified
  matrix/preconditioner pair 
*/
KsmPreconditioner::KsmPreconditioner( TACSMat * _mat, 
				      TACSPc * _pc ){
  mat = _mat;
  pc = _pc;
  mat->incref();
  pc->incref();
}

KsmPreconditioner::~KsmPreconditioner(){
  mat->decref();
  pc->decref();
}

/*
  Create a solution vector using the direct preconditioner
*/
TACSVec * KsmPreconditioner::createVec(){ 
  return mat->createVec(); 
}

/*
  Set the operators for the direct preconditioner
*/
void KsmPreconditioner::setOperators( TACSMat * _mat, TACSPc * _pc ){
  if (_mat){ _mat->incref(); }
  mat->decref();
  mat = _mat;
  
  if (_pc){ _pc->incref(); }
  pc->decref();
  pc = _pc;
}

void KsmPreconditioner::getOperators( TACSMat ** _mat, TACSPc ** _pc ){
  *_mat = mat;
  *_pc = pc;
}
 
/*
  'Solve' the problem - this is only really a solution if a direct
  solver is used. This is often the case within TACS.
*/
void KsmPreconditioner::solve( TACSVec * b, TACSVec * x, 
			       int zero_guess ){
  pc->applyFactor(b, x);
}

/*
  The solution tolerances have no effect.
*/
void KsmPreconditioner::setTolerances( double _rtol, double _atol ){}

/*
  Set a solution monitor - don't monitor anything since there are 
  no real iterations
*/
void KsmPreconditioner::setMonitor( KSMPrint * monitor ){
  // Delete the monitor if it's allocated as an argument
  monitor->incref();
  monitor->decref();
}
