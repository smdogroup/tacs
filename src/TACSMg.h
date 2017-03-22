#ifndef TACS_MG_H
#define TACS_MG_H

/*
  Multigrid code for solution/preconditioning problems.

  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include "TACSAssembler.h"
#include "BVec.h"
#include "BVecInterp.h"

/*
  This class implements a geometric multi-grid solution method as
  either a preconditioner or a solution method.

  This uses geometric multigrid where a number of related
  finite-element models are used to assemble a related series of
  linear systems. The interpolation/prolongation operators transfer
  the solution between different grid levels. 
  
  This automates the recursive application of multi-grid, but the
  requirement to formulate the nested series of TACS models remains
  with the user. This is automated somewhat with the TACSMesh class.

  This code can be used as either a stand-alone solver, or as a
  preconditioner for a Krylov-subspace method. As such, this could be
  included within the block-CSR matrix factorization code directory,
  but is included here instead due to the nested TACSAssembler
  formation requirements.

  Note that the lowest multi-grid level should be solved with a direct
  solver. By default, this is the parallel direct solver implemented
  in TACS. This solver has good parallel scalability, but poor
  scalability with the numbers of degrees of freedom in the
  model. This should be considered when deciding on the number of
  levels of multi-grid to use.
*/
class TACSMg : public TACSPc {
 public:
  TACSMg( MPI_Comm comm, int _nlevels, double _sor_omega, 
	  int _sor_iters, int _sor_symmetric );
  ~TACSMg();

  // Set the data for the multi-grid level
  // -------------------------------------
  void setLevel( int level, TACSAssembler *_tacs,
                 TACSBVecInterp *interp,
		 int _iters=1 );
    
  // Set the state/design variables of all lower finite-element models
  // -----------------------------------------------------------------
  void setVariables( TACSBVec * vec ); 
  void setDesignVars( const TacsScalar dvs[], int numDVs );

  // Assemble the given finite-element matrix at all levels
  // ------------------------------------------------------
  void assembleJacobian( double alpha, double beta, double gamma,
                         TACSBVec *res=NULL, 
                         MatrixOrientation matOr=NORMAL );
  void assembleMatType( ElementMatrixType matType=STIFFNESS_MATRIX, 
			MatrixOrientation matOr=NORMAL );

  // Methods required by the TACSPc class
  // ------------------------------------
  void applyFactor( TACSVec *x, TACSVec *y );
  void factor();

  // Solve the problem using the full multi-grid method
  // --------------------------------------------------
  void solve( TACSBVec *bvec, TACSBVec *xvec, int max_iters=200, 
	      double rtol=1e-8, double atol=1e-30 );

  // Retrieve the matrix from the specified level
  // --------------------------------------------
  TACSMat *getMat( int level );  

  // Set the solution monitor context
  // --------------------------------
  void setMonitor( KSMPrint *_monitor );

 private:
  // Recursive function to apply multi-grid at each level
  void applyMg( int level ); 

  // The MPI communicator for this object
  MPI_Comm comm;

  // Monitor the solution
  KSMPrint *monitor;

  // The SOR data 
  int sor_iters, sor_symmetric;
  double sor_omega;

  // The number of multi-grid levels
  int nlevels;

  // The TACSAssembler object for each level
  TACSAssembler **tacs;
  int *iters;

  // The solution, right-hand-side and residual on each level
  TACSBVec **x, **b, **r;

  // The interpolation operators
  TACSBVecInterp **interp;

  // The matrices/preconditioner objects required for multigrid
  TACSMat *root_mat; // The root matrix 
  TACSPc *root_pc; // The root direct solver
  TACSMat **mat; // The matrices associated with each level
  TACSPc **pc; // The smoothers for all but the lowest level
};

#endif
