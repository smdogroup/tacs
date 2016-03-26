#include "TACSMg.h"

/*
  Implementation of Multigrid within TACS

  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

/*
  Set up the TACS multi-grid object with a given number of multi-grid
  levels.

  Each level is initialized with the given successive over relaxation
  (SOR) factor (sor_omega), and the number of relaxations to perform
  at each level. The flag sor_symmetric indicates whether or not to
  use a symmetric SOR iteration.

  input:
  comm:          MPI communicator
  nlevels:       number of multi-grid levels
  sor_omega:     SOR factor
  sor_iters:     number of SOR iterations to perform at each level
  sor_symmetric: symmetric SOR flag
*/
TACSMg::TACSMg( MPI_Comm _comm, int _nlevels, 
		double _sor_omega, int _sor_iters, 
		int _sor_symmetric ){
  // Copy over the data 
  comm = _comm;
  nlevels = _nlevels;
  sor_omega = _sor_omega;
  sor_iters = _sor_iters;
  sor_symmetric = _sor_symmetric;

  if (nlevels < 2){
    fprintf(stderr, "Multigrid with fewer than 2 levels does \
not make any sense!\n");
  }

  // Create the list of pointers to the TACS objects at each level 
  tacs = new TACSAssembler*[ nlevels ];
  
  // Number of smoothing operations to be performed at each iteration
  iters = new int[ nlevels ];    

  // The solution, residual and right-hand-side at each iteration
  x = new BVec*[ nlevels ];
  b = new BVec*[ nlevels ];
  r = new BVec*[ nlevels ];

  // Initialie the data in the arrays
  for ( int i = 0; i < nlevels; i++ ){ 
    iters[i] = 1; // defaults to one - a V cycle
    tacs[i] = NULL;
    r[i] = NULL;
    x[i] = NULL;
    b[i] = NULL;
  }

  // Create the pointers to the matrices
  mat = new PMat*[ nlevels-1 ];
 
  // should re-assemble this matrix
  restrct = new BVecInterp*[ nlevels-1 ];
  interp = new BVecInterp*[ nlevels-1 ];
  psor = new PSOR*[ nlevels-1 ]; 

  for ( int i = 0; i < nlevels-1; i++ ){
    restrct[i] = NULL;
    interp[i] = NULL;
    mat[i] = NULL;
    psor[i] = NULL;
  }

  monitor = NULL;
  root_mat = NULL;
  root_pc = NULL;
}

/*
  Deallocate/dereference the data stored internally
*/
TACSMg::~TACSMg(){
  for ( int i = 0; i < nlevels; i++ ){
    if (tacs[i]){ tacs[i]->decref(); }
    if (r[i]){ r[i]->decref(); }
    if (b[i]){ b[i]->decref(); }
    if (x[i]){ x[i]->decref(); }
  }

  for ( int i = 0; i < nlevels-1; i++ ){
    if (mat[i]){ mat[i]->decref(); }
    if (restrct[i]){ restrct[i]->decref(); }
    if (interp[i]){ interp[i]->decref(); }
    if (psor[i]){ psor[i]->decref(); }
  }

  if (monitor){ monitor->decref(); }
  if (root_mat){ root_mat->decref(); }
  if (root_pc){ root_pc->decref(); }

  delete [] tacs;
  delete [] mat;
  delete [] iters;
  delete [] x;
  delete [] r;
  delete [] b;

  delete [] restrct;
  delete [] interp;
  delete [] psor;
}

/*
  Set the data for the given multi-grid level.

  This consists of setting the TACS finite-element model, and the
  restriction and interpolation operators for this level. (These are
  not requried for the lowest level)

  input:
  level:     the multigrid level
  tacs:      the TACSAssembler object
  restrict:  the restriction operator
  interp:    the interpolation operator
  iters:     the number of iterations to take at this level
*/
void TACSMg::setLevel( int level, TACSAssembler * _tacs,
		       BVecInterp * _restrct, BVecInterp * _interp, 
		       int _iters ){
  tacs[level] = _tacs;
  tacs[level]->incref();

  iters[level] = 1;
  if (_iters > 0){
    iters[level] = _iters;
  }

  // Only define the restriction/interpolation 
  // operators for level < nlevels-1
  if (level < nlevels-1){
    if (!_restrct && !_interp){
      fprintf(stderr, "TACSMg: Must define restriction and prolongation\
 operators for all but the coarsest problem\n");
    }

    restrct[level] = _restrct;    
    restrct[level]->incref();

    interp[level] = _interp;
    interp[level]->incref();
    
    mat[level] = tacs[level]->createMat();
    mat[level]->incref();
  }
  else {
    // Set up the root matrix
    root_mat = tacs[level]->createFEMat();
    root_mat->incref();

    // Set up the root preconditioner/solver
    int lev = 10000;
    double fill = 15.0;
    int reorder_schur = 1; 
    root_pc = new PcScMat(root_mat, lev, fill, reorder_schur);
    root_pc->incref();
  }
  
  if (level > 0){
    x[level] = tacs[level]->createVec();
    x[level]->incref();

    b[level] = tacs[level]->createVec();
    b[level]->incref();
  }

  r[level] = tacs[level]->createVec();
  r[level]->incref();
}

/*
  Set the state variables for this, and all subsequent levels.

  input:
  vec:      the input vector of state variables
*/
void TACSMg::setVariables( BVec * vec ){
  tacs[0]->setVariables(vec);

  for ( int i = 0; i < nlevels-1; i++ ){
    restrct[i]->mult(x[i], x[i+1]);
    x[i+1]->applyBCs();
    tacs[i+1]->setVariables(x[i+1]);
  }
}

/*
  Set the design variables for all multi-grid levels.

  This call ensures that all the TACSAssembler objects, and the
  objects they reference, share the same set of design variables.

  input:
  dvs:     the design variable values
  numDVs:  the number of design variables 
*/
void TACSMg::setDesignVars( const TacsScalar dvs[], int numDVs ){
  for ( int i = 0; i < nlevels; i++ ){
    if (tacs[i]){
      tacs[i]->setDesignVars(dvs, numDVs);
    }
  }
}

/*
  Factor the SOR objects for each level of multi-grid, and 
  the direct solver for the lowest level.
*/
void TACSMg::factor(){
  for ( int i = 0; i < nlevels-1; i++ ){
    if (psor[i]){ psor[i]->factor(); }
  }

  if (root_pc){ root_pc->factor(); }
}

/*
  Set up the multi-grid data by computing the matrices at each
  multi-grid level within the problem. 
*/
void TACSMg::assembleJacobian( BVec *res, 
                               double alpha, double beta, double gamma,
                               MatrixOrientation matOr ){
  // Assemble the matrices if they are locally owned, otherwise assume
  // that they have already been assembled
  if (tacs[0]){
    tacs[0]->assembleJacobian(res, mat[0], 
                              alpha, beta, gamma, matOr);
  }

  for ( int i = 0; i < nlevels-1; i++ ){
    if (tacs[i]){
      tacs[i]->assembleJacobian(NULL, mat[i], 
                                alpha, beta, gamma, matOr);
    }
  }

  // Assemble the coarsest problem
  if (tacs[nlevels-1]){
    tacs[nlevels-1]->assembleJacobian(NULL, root_mat, 
                                      alpha, beta, gamma, matOr);
  }

  // For all but the lowest level, set up the SOR object
  for ( int i = 0; i < nlevels-1; i++ ){
    if (!psor[i]){
      // Do not zero the initial guess for the PSOR object
      int zero_guess = 0; 
      psor[i] = new PSOR(mat[i], zero_guess, sor_omega, 
			 sor_iters, sor_symmetric);
      psor[i]->incref();
    }
  }
}

/*
  Set up the multi-grid data by computing the matrices at each
  multi-grid level within the problem. 
*/
void TACSMg::assembleMatType( ElementMatrixType matType, 
			      MatrixOrientation matOr ){
  // Assemble the matrices if they are locally owned, otherwise assume
  // that they have already been assembled
  for ( int i = 0; i < nlevels-1; i++ ){
    if (tacs[i]){
      tacs[i]->assembleMatType(matType, mat[i], matOr);
    }
  }

  // Assemble the coarsest problem 
  if (tacs[nlevels-1]){
    tacs[nlevels-1]->assembleMatType(matType, root_mat, matOr);
  }

  // For all but the lowest level, set up the SOR object
  for ( int i = 0; i < nlevels-1; i++ ){
    if (!psor[i]){
      // Do not zero the initial guess for the PSOR object
      int zero_guess = 0; 
      psor[i] = new PSOR(mat[i], zero_guess, sor_omega, 
			 sor_iters, sor_symmetric);
      psor[i]->incref();
    }
  }
}

/*
  Retrieve the matrix at the specified multigrid level.

  Return NULL if no such matrix exists.
*/
TACSMat * TACSMg::getMat( int level ){
  if (level >= 0 && level < nlevels-1){
    return mat[level];
  }
  
  return NULL;
}

/*
  Set the monitor to use internally for printing out convergence data.
*/
void TACSMg::setMonitor( KSMPrint * _monitor ){
  if (_monitor){
    _monitor->incref();
  }
  if (monitor){
    monitor->decref();
  }
  monitor = _monitor;
}

/*
  Repeatedly apply the multi-grid method until the problem is solved
*/
void TACSMg::solve( BVec * bvec, BVec * xvec, int max_iters,
		    double rtol, double atol ){
  b[0] = bvec; // Set the RHS at the finest level
  x[0] = xvec; // Set the solution at the finest level
  
  // Compute the initial residual and multiply
  TacsScalar rhs_norm = 0.0;
  for ( int i = 0; i < max_iters; i++ ){
    TacsScalar norm = applyMg(0);    
    if (monitor){ monitor->printResidual(i, norm); }
    if (i == 0){ 
      rhs_norm = norm; 
    }

    if (norm < atol || norm < rtol*rhs_norm){
      break;
    }
  } 

  b[0] = NULL;
  x[0] = NULL;
}

/*  
  Apply the multi-grid preconditioner to try and solve the problem
  x = A^{-1} b

  Assume an initial guess of zero. 
*/
void TACSMg::applyFactor( TACSVec * bvec, TACSVec * xvec ){
  b[0] = dynamic_cast<BVec*>(bvec); // Set the RHS at the finest level
  x[0] = dynamic_cast<BVec*>(xvec); // Set the solution at the finest level

  if (b[0] && x[0]){
    x[0]->zeroEntries();
    applyMg(0);
  }
  else {
    fprintf(stderr, "TACSMg type error: Input/output must be BVec\n");
  }
  
  b[0] = NULL;
  x[0] = NULL;
}

/*
  Apply the multigrid algorithm from McCormick and Runge 1983 - but
  with added post-smoothing. 

  This function applies multigrid recursively by smoothing the
  residual, restricting to the next level, applying multigrid, then
  post-smoothing.
*/
TacsScalar TACSMg::applyMg( int level ){
  // Pre-smooth at the current level
  psor[level]->applyFactor(b[level], x[level]);  

  // Compute r[level] = b[level] - A*x[level]
  mat[level]->mult(x[level], r[level]);
  r[level]->axpby(1.0, -1.0, b[level]);

  TacsScalar norm = r[level]->norm();

  // Restrict the residual to the next lowest level 
  // to form the RHS at that level
  restrct[level]->mult(r[level], b[level+1]);
  b[level+1]->applyBCs();

  // If we've made it to the lowest level, apply the direct solver
  // otherwise, perform multigrid on the next-lowest level
  if (level+1 == nlevels-1){
    // Perform a direct solve on the smallest grid
    root_pc->applyFactor(b[nlevels-1], x[nlevels-1]); 
  }
  else {
    x[level+1]->zeroEntries();

    // Perform iters[level] cycle at the next lowest level
    for ( int k = 0; k < iters[level]; k++ ){
      applyMg(level+1);
    }
  }

  // Interpolate back from the next lowest level
  interp[level]->multAdd(x[level+1], x[level], x[level]);
  x[level]->applyBCs();

  // Post-Smooth the residual
  psor[level]->applyFactor(b[level], x[level]);  

  return norm;
}
