/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TACSMg.h"

/*
  Implementation of geometric multigrid using TACS
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
  x = new TACSBVec*[ nlevels ];
  b = new TACSBVec*[ nlevels ];
  r = new TACSBVec*[ nlevels ];

  // Initialie the data in the arrays
  for ( int i = 0; i < nlevels; i++ ){
    iters[i] = 1; // defaults to one - a V cycle
    tacs[i] = NULL;
    r[i] = NULL;
    x[i] = NULL;
    b[i] = NULL;
  }

  // Create the pointers to the matrices
  mat = new TACSMat*[ nlevels-1 ];

  // should re-assemble this matrix
  interp = new TACSBVecInterp*[ nlevels-1 ];
  pc = new TACSPc*[ nlevels-1 ];

  for ( int i = 0; i < nlevels-1; i++ ){
    interp[i] = NULL;
    mat[i] = NULL;
    pc[i] = NULL;
  }

  monitor = NULL;
  root_mat = NULL;
  root_pc = NULL;

  // Total time for smoothing and interpolation for each level
  cumulative_level_time = new double[ nlevels ];
  memset(cumulative_level_time, 0, nlevels*sizeof(double));
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
    if (interp[i]){ interp[i]->decref(); }
    if (pc[i]){ pc[i]->decref(); }
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
  delete [] interp;
  delete [] pc;
  delete [] cumulative_level_time;
}

/*
  Set the data for the given multi-grid level.

  This consists of setting the TACS finite-element model, and the
  restriction and interpolation operators for this level. (These are
  not requried for the lowest level)

  input:
  level:     the multigrid level
  tacs:      the TACSAssembler object
  interp:    the interpolation operator
  iters:     the number of iterations to take at this level
*/
void TACSMg::setLevel( int level, TACSAssembler *_tacs,
                       TACSBVecInterp *_interp,
                       int _iters, TACSMat *_mat,
                       TACSPc *_smoother ){
  tacs[level] = _tacs;
  tacs[level]->incref();

  if ((!_mat && _smoother) || (!_smoother && _mat)){
    fprintf(stderr,
            "TACSMg: You must define both the matrix and preconditioner\n");
    _mat = NULL;
    _smoother = NULL;
  }

  iters[level] = 1;
  if (_iters > 0){
    iters[level] = _iters;
  }

  // Only define the restriction/interpolation
  // operators for level < nlevels-1
  if (level < nlevels-1){
    if (!_interp){
      fprintf(stderr, "TACSMg: Must define prolongation\
 operators for all but the coarsest problem\n");
    }
    interp[level] = _interp;
    interp[level]->incref();

    if (_mat){
      _mat->incref();
      if (mat[level]){
        mat[level]->decref();
      }
      mat[level] = _mat;

      _smoother->incref();
      if (pc[level]){
        pc[level]->decref();
      }
      pc[level] = _smoother;
    }
    else {
      TACSPMat *pmat = tacs[level]->createMat();
      mat[level] = pmat;
      mat[level]->incref();

      // Do not zero the initial guess for the PSOR object
      int zero_guess = 0;
      pc[level] = new TACSGaussSeidel(pmat, zero_guess, sor_omega,
                                      sor_iters, sor_symmetric);
      pc[level]->incref();
    }
  }
  else {
    if (_mat){
      _mat->incref();
      if (root_mat){
        root_mat->decref();
      }
      root_mat = _mat;

      _smoother->incref();
      if (root_pc){
        root_pc->decref();
      }
      root_pc = _smoother;
    }
    else {
      // Set up the root matrix
      FEMat *femat = tacs[level]->createFEMat();
      root_mat = femat;
      root_mat->incref();

      // Set up the root preconditioner/solver
      int lev = 10000;
      double fill = 15.0;
      int reorder_schur = 1;
      PcScMat *_pc = new PcScMat(femat, lev, fill, reorder_schur);
      // _pc->setMonitorFactorFlag(1);
      // _pc->setMonitorBackSolveFlag(1);
      root_pc = _pc;
      root_pc->incref();
    }
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
void TACSMg::setVariables( TACSBVec *vec ){
  tacs[0]->setVariables(vec);

  for ( int i = 0; i < nlevels-1; i++ ){
    if (i == 0){
      interp[i]->multWeightTranspose(vec, x[i+1]);
    }
    else {
      interp[i]->multWeightTranspose(x[i], x[i+1]);
    }
    x[i+1]->applyBCs(tacs[i+1]->getBcMap());
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
    if (pc[i]){ pc[i]->factor(); }
  }
  if (root_pc){ root_pc->factor(); }
}

/*
  Set up the multi-grid data by computing the matrices at each
  multi-grid level within the problem.
*/
void TACSMg::assembleJacobian( double alpha, double beta, double gamma,
                               TACSBVec *res,
                               MatrixOrientation matOr ){
  // Assemble the matrices if they are locally owned, otherwise assume
  // that they have already been assembled
  if (tacs[0]){
    tacs[0]->assembleJacobian(alpha, beta, gamma,
                              res, mat[0], matOr);
  }

  for ( int i = 1; i < nlevels-1; i++ ){
    if (tacs[i]){
      tacs[i]->assembleJacobian(alpha, beta, gamma,
                                NULL, mat[i], matOr);
    }
  }

  // Assemble the coarsest problem
  if (tacs[nlevels-1]){
    tacs[nlevels-1]->assembleJacobian(alpha, beta, gamma,
                                      NULL, root_mat, matOr);
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
}

/*
  Assemble a linear combination of matrices with the specified linear
  combination
*/
void TACSMg::assembleMatCombo( ElementMatrixType matTypes[],
                               TacsScalar scale[], int nmats,
                               MatrixOrientation matOr ){
  for ( int i = 0; i < nlevels-1; i++ ){
    if (tacs[i]){
      tacs[i]->assembleMatCombo(matTypes, scale, nmats, mat[i], matOr);
    }
  }

  // Assemble the coarsest problem
  if (tacs[nlevels-1]){
    tacs[nlevels-1]->assembleMatCombo(matTypes, scale, nmats, root_mat, matOr);
  }
}

/*
  Get the matrix operator associated with the highest
*/
void TACSMg::getMat( TACSMat **_mat ){
  *_mat = mat[0];
}

/*
  Retrieve the matrix at the specified multigrid level.

  Return NULL if no such matrix exists.
*/
TACSMat *TACSMg::getMat( int level ){
  if (level >= 0 && level < nlevels-1){
    return mat[level];
  }

  return NULL;
}

/*
  Retrieve the TACSAssembler model
*/
TACSAssembler *TACSMg::getTACS( int level ){
  if (level >= 0 && level < nlevels-1){
    return tacs[level];
  }
  return NULL;
}

/*
  Retrieve the TACSBVecInterp object for the given level
*/
TACSBVecInterp *TACSMg::getInterpolation( int level ){
  if (level >= 0 && level < nlevels-1){
    return interp[level];
  }
  return NULL;
}

/*
  Set the monitor to use internally for printing out convergence data.
*/
void TACSMg::setMonitor( KSMPrint *_monitor ){
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
void TACSMg::solve( TACSBVec *bvec, TACSBVec *xvec, int max_iters,
                    double rtol, double atol ){
  b[0] = bvec; // Set the RHS at the finest level
  x[0] = xvec; // Set the solution at the finest level

  // Compute the initial residual and multiply
  TacsScalar rhs_norm = 0.0;
  for ( int i = 0; i < max_iters; i++ ){
    applyMg(0);
    TacsScalar norm = r[0]->norm();
    if (monitor){ monitor->printResidual(i, norm); }
    if (i == 0){
      rhs_norm = norm;
    }

    if (TacsRealPart(norm) < atol ||
        TacsRealPart(norm) < rtol*TacsRealPart(rhs_norm)){
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
void TACSMg::applyFactor( TACSVec *bvec, TACSVec *xvec ){
  // Set the RHS at the finest level
  b[0] = dynamic_cast<TACSBVec*>(bvec);

  // Set the solution at the finest level
  x[0] = dynamic_cast<TACSBVec*>(xvec);

  if (b[0] && x[0]){
    x[0]->zeroEntries();
    if (monitor){
      memset(cumulative_level_time, 0, nlevels*sizeof(double));
    }
    applyMg(0);
    if (monitor){
      for ( int k = 0; k < nlevels; k++ ){
        char descript[128];
        sprintf(descript, "TACSMg cumulative level %2d time %15.8e\n",
                k, cumulative_level_time[k]);
        monitor->print(descript);
      }
    }
  }
  else {
    fprintf(stderr, "TACSMg type error: Input/output must be TACSBVec\n");
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
void TACSMg::applyMg( int level ){
  // If we've made it to the lowest level, apply the direct solver
  // otherwise, perform multigrid on the next-lowest level
  if (level == nlevels-1){
    double t1 = 0.0;
    if (monitor){ t1 = MPI_Wtime(); }
    root_pc->applyFactor(b[level], x[level]);
    if (monitor){ cumulative_level_time[level] += MPI_Wtime() - t1; }
    return;
  }

  // Perform iters[level] cycle at the next lowest level
  double t1 = 0.0;
  if (monitor){ t1 = MPI_Wtime(); }
  for ( int k = 0; k < iters[level]; k++ ){
    // Pre-smooth at the current level
    pc[level]->applyFactor(b[level], x[level]);

    // Compute r[level] = b[level] - A*x[level]
    mat[level]->mult(x[level], r[level]);
    r[level]->axpby(1.0, -1.0, b[level]);

    // Restrict the residual to the next lowest level
    // to form the RHS at that level
    interp[level]->multTranspose(r[level], b[level+1]);
    b[level+1]->applyBCs(tacs[level+1]->getBcMap());
    x[level+1]->zeroEntries();

    applyMg(level+1);

    // Interpolate back from the next lowest level
    interp[level]->multAdd(x[level+1], x[level], x[level]);
    x[level]->applyBCs(tacs[level]->getBcMap());
  }

  // Post-Smooth the residual
  pc[level]->applyFactor(b[level], x[level]);

  if (monitor){
    cumulative_level_time[level] += MPI_Wtime() - t1;
  }
}
