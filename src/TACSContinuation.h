#ifndef TACS_CONTINUATION_H
#define TACS_CONTINUATION_H

/*!
  Continuation methods for bifurcation analysis/path following.

  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include "TACSAssembler.h"

/*!
  Matrix type required for continuation methods.
  
  This solves a system of equations that takes the following
  form:
  
  [ A | r ][x]
  .        [y] = [b]

  subject to a single constraint:

  [t^{T} | s][x]
  .          [y] = 0

  The advantage of this approach is that the constraint is satisified
  exactly, even if the system of equations is only solved approximately.
*/
class PathMat : public TACSMat {
 public:
  PathMat( TACSMat * _A, 
	   TACSVec * _r, TACSVec * _t, TacsScalar s );
  ~PathMat();

  // Return the vectors associated with the constraint
  // -------------------------------------------------
  void getVectors( TACSVec ** _r, TACSVec ** _t );

  // Reset the constraint values based on the value s
  // ------------------------------------------------
  void resetConstraint( TacsScalar s );

  // Multiply x <-- Qx, return the value of the n+1-th row
  // -----------------------------------------------------
  TacsScalar extract( TACSVec * x );

  // The required TACSMat implementation - performs no operations
  // ------------------------------------------------------------
  void zeroEntries(){}
  void addValues( int nrow, const int * row, int ncol, const int * col,
			  int nv, int mv, const TacsScalar * values ){}
  void applyBCs(){}
  void beginAssembly(){}
  void endAssembly(){}

  // Specific implementation for PathMath
  // ------------------------------------
  void getSize( int * _nr, int * _nc );
  void mult( TACSVec * x, TACSVec * y );
  TACSVec * createVec();

 private:
  TACSMat *A;
  TACSVec *r, *t, *xtmp;
  TacsScalar tn, wn;
  TacsScalar s, w;
};

/*!
  Callback class for the continuation solver.  This outputs a monitor
  file and/or solution file at every iteration. Passing in NULL file
  names/object values is legal.  
*/
class TACSContCallback : public TACSObject {
 public:
  TACSContCallback( MPI_Comm _comm,
		    const char * file_prefix,
		    const char * monitor_file,
		    int * _var_nums, int _nvars );
  ~TACSContCallback();

  void setUnitBCDisplacement( BVec * unit );

  void initialize();
  void iteration( int iter, TacsScalar lambda, TacsScalar dlambda_ds, 
		  int loadCase, TACSAssembler * tacs,
		  BVec * ans );
  void finalize();
 private:
  MPI_Comm comm;
  char * file_prefix;
  char * monitor_file;

  // What variables to output to the file
  int nvars; 
  int *var_nums;
  TacsScalar *var_values;

  FILE * mfp;

  // Information for computing the force exherted by an imposed
  // boundary displacement
  BVec * unit_bc_displacement;
  BVec * res_no_bcs;
};

/*!
  Arc-length continuation algorithm. This method is based on a
  constant displacement per iteration. The initial increment in
  displacement is determined based on a linearized estimate of the
  equilibrium equations. All subsequent iterations match the initial
  displacement.

  There are two exit conditions that can be set, or no exit condition
  at all (that defaults to a fixed number of iterations).

  The following exit conditions are implemented such that if any are
  satisfied, the iterations terminate.

  1. The path derivative d(lambda)/ds falls below a preset value
  (typically set to zero). d(lambda)/ds corresponds to structural
  collapse. Setting the d(lambda)/ds to a large, negative value
  removes this constraint.

  2. A failure constraint is violated. The failure constraint is given
  by a value of the TACSFunction passed to TACSContinuation algorithm.

  3. The maximum number of iterations is exceeded.
*/
class TACSContinuation : public TACSObject { 
 public:
  TACSContinuation( TACSAssembler * _tacs, int _max_continuation_ters );
  ~TACSContinuation();

  // Set the callback for monitoring purposes
  // ----------------------------------------
  void setCallback( TACSContCallback * _callback );

  // Set iteration limits/tolerances for the modified Newton iteration
  // -----------------------------------------------------------------
  void setCorrectionIters( int _max_correction_iters,
			   int _max_correction_restarts );
  void setCorrectionTol( double rtol, double dtol );

  // Set the termination conditions
  // ------------------------------
  void setTermFunction( TACSFunction * func, 
			TacsScalar term_value );
  void setTermLambdaRate( TacsScalar term_dlambda_ds );

  // Perform a continuation solve using a linearized arc-length constraint
  // ---------------------------------------------------------------------
  void solve_tangent( int loadCase,
		      TACSMat * mat, TACSPc * pc, 
		      TACSKsm * ksm,
                      TacsScalar lambda_init,
		      TacsScalar target_delta_lambda,
                      KSMPrint * ksm_print = NULL );

  // Retrieve information about the solve
  // ------------------------------------
  int getNumIterations();
  void getSolution( int iter, TacsScalar * lambda, 
                    TacsScalar * dlambda_ds );
  
 private:
  // Solution parameters
  int max_continuation_iters;
  int max_correction_iters;
  int max_correction_restarts;

  // relative and divergence tolerances for the correction
  double correction_rtol;
  double correction_dtol;

  // relative and absolute tolerances for the Krylov updates
  double krylov_rtol, krylov_atol;

  // relative and absolute tolerances for the tangent computation
  double tangent_rtol, tangent_atol;
  
  // Solution variables
  TACSAssembler * tacs; 

  // Callback used to output information at every iteration
  TACSContCallback * callback;

  // Information to store the iteration history
  int iteration_count; // The number of iterations actually used
  TacsScalar * lambda_history;     // The history of the parameter
  TacsScalar * dlambda_ds_history; // The history of dlambda/ds

  // Termination condition information
  TACSFunction * term_function;
  TacsScalar term_function_value;  
  TacsScalar dlambda_ds_term_value;
};

#endif
