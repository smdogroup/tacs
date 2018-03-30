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

#include "TACSBuckling.h"
#include "tacslapack.h"

/*
  Implementation of the buckling/frequency analysis
*/

/*
  Linear buckling analysis object.

  This object uses a shift and invert strategy in conjunction with a
  full-orthogonalization Lanczos method to compute the buckling
  eigenvalues and eigenvectors. The method requires objects for the
  stiffness and geometric stiffness matrices as well as an auxiliary
  matrix that is used to store a linear combination of the two. The
  solver must be associated with the auxiliary matrix.

  Note: all the matrices supplied must be of the same type and support
  copy/axpy/axpby etc. operations.

  input:
  tacs:         The TACS model corresponding to the analysis problem
  sigma:        The spectral shift
  gmat:         The geometric stiffness matrix
  kmat:         The stiffness matrix
  aux_mat:      The auxiliary matrix associated with the solver
  solver:       Whatever KSM object you create
  max_lanczos:  Maximum size of the projected subspace
  num_eigvals:  Number of converged eigenvalues required
  eig_tol:      Tolerance of the eigenvalues
*/
TACSLinearBuckling::TACSLinearBuckling( TACSAssembler *_tacs,
                                        TacsScalar _sigma,
                                        TACSMat *_gmat, TACSMat *_kmat,
                                        TACSMat *_aux_mat, TACSKsm *_solver,
                                        int _max_lanczos_vecs,
                                        int _num_eigvals,
                                        double _eig_tol ){
  // Copy pointer to the TACS assembler object
  tacs = _tacs;
  tacs->incref();

  // Store the matrices required
  aux_mat = _aux_mat;
  gmat = _gmat;
  kmat = _kmat;
  aux_mat->incref();
  gmat->incref();
  kmat->incref();

  // The Krylov subspace method (KSM) associated with the solver
  solver = _solver;
  solver->incref();

  // Get the operators
  TACSMat *mat;
  solver->getOperators(&mat, &pc);

  // Check that the matrix associated with the solver object is the auxiliary
  // matrix. If not, complain about it.
  if (mat != aux_mat){
    fprintf(stderr, "TACSBuckling: Error, solver must be associated with the \
auxiliary matrix\n");
  }

  // Check that the auxiliary matrix, geometric stiffness and stiffness
  // matrices are different objects
  if (aux_mat == kmat){
    fprintf(stderr, "TACSBuckling: Error, stiffness and auxiliary matrices \
must be different instances\n");
  }
  if (aux_mat == gmat){
    fprintf(stderr, "TACSBuckling: Error, geometric stiffness and auxiliary matrices \
must be different instances\n");
  }
  if (gmat == kmat){
    fprintf(stderr, "TACSBuckling: Error, geometric stiffness and stiffness matrices \
must be different instances\n");
  }

  // Check if the preconditioner is actually a multigrid object. If
  // so, then we have to allocate extra data to store things for each
  // multigrid level.
  mg = dynamic_cast<TACSMg*>(pc);

  // Store the spectral shift info
  sigma = _sigma;

  // The eigenvalues solver options
  max_lanczos_vecs = _max_lanczos_vecs;
  num_eigvals = _num_eigvals;
  eig_tol = _eig_tol;

  // Allocate the shift and invert inner product
  ep_op = new EPBucklingShiftInvert(sigma, solver, kmat);
  ep_op->incref();

  // Allocate the eigenvalue solver
  sep = new SEP(ep_op, max_lanczos_vecs, SEP::FULL, tacs->getBcMap());
  sep->incref();
  sep->setTolerances(eig_tol, SEP::SMALLEST_MAGNITUDE,
                     num_eigvals);

  // Allocate temporary local vectors
  res = tacs->createVec();
  update = tacs->createVec();
  eigvec = tacs->createVec();
  path = tacs->createVec();
  res->incref();
  update->incref();
  eigvec->incref();
  path->incref();
}

/*
  Destructor object for the buckling object
*/
TACSLinearBuckling::~TACSLinearBuckling(){
  // Dereference the matrix objects
  aux_mat->decref();
  gmat->decref();
  kmat->decref();

  // Dereference the solver/tacs
  tacs->decref();
  solver->decref();

  // Deallocate the solvers
  ep_op->decref();
  sep->decref();

  // Deallocate the vectors
  path->decref();
  res->decref();
  update->decref();
  eigvec->decref();
}

/*
  Retrieve the value of sigma
*/
TacsScalar TACSLinearBuckling::getSigma(){
  return sigma;
}

/*
  Set a different value for the shift parameter
*/
void TACSLinearBuckling::setSigma( TacsScalar _sigma ){
  sigma = _sigma;
  ep_op->setSigma(sigma);
}

/*
  Solve the linearized buckling problem about x = 0.

  This code determine the lowest magnitude eigenvalue such that,

  K x = - lambda * G x

  where K is the stiffness matrix and G is the geometric stiffness
  matrix.  The difficulty is that the geometric stiffness matrix is
  not positive definite necessarily. As a a result, a modified
  eigenvalue problem must be solved instead. In order to attain better
  convergence behaviour from the eigensolver, it is necessary to
  employ a shift and invert strategy to isolate a segment of the
  eigen-spectrum that is of primary interest in the analysis.  For
  buckling problems, the lowest buckling modes must be determined.
  The shift-invert strategy is employed such that the eigenvector
  remains unmodified, while the eigenvalue is shifted.

  (K + sigma G) x = - (lambda - sigma) G x = (lambda - sigma)/lambda K x

  The modified problem is then,

  (K + sigma G)^{-1} K x = lambda/(lambda - sigma) x
*/
void TACSLinearBuckling::solve( TACSVec *rhs, KSMPrint *ksm_print ){
  // Zero the variables
  tacs->zeroVariables();

  if (mg){
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    mg->assembleJacobian(alpha, beta, gamma, res);
    mg->factor();

    // Add the right-hand-side due to external forces
    if (rhs){
      //res->axpy(1.0, rhs);
      tacs->applyBCs(res);
    }

    // Solve for the load path
    solver->solve(res, path);
    path->scale(-1.0);
    tacs->setVariables(path);

    // Assemble the linear combination of the stiffness matrix
    // and geometric stiffness matrix
    ElementMatrixType matTypes[2] = {STIFFNESS_MATRIX,
                                     GEOMETRIC_STIFFNESS_MATRIX};
    TacsScalar scale[2] = {1.0, sigma};
    mg->assembleMatCombo(matTypes, scale, 2);

    // Assemble the geometric stiffness matrix and the stiffness matrix itself
    tacs->assembleMatType(STIFFNESS_MATRIX, kmat);
    tacs->assembleMatType(GEOMETRIC_STIFFNESS_MATRIX, gmat);
  }
  else{
    // Compute the stiffness matrix and copy the values to the
    // auxiliary matrix used to solve for the load path.
    tacs->assembleMatType(STIFFNESS_MATRIX, kmat);
    aux_mat->copyValues(kmat);

    pc->factor();
    tacs->assembleRes(res);

    // If need to add rhs
    if (rhs){
      res->axpy(-1.0, rhs);
      tacs->applyBCs(res);
    }

    // Solve for the load path and set the variables
    solver->solve(res, path);
    path->scale(-1.0);
    tacs->setVariables(path);

    // Assemble the stiffness and geometric stiffness matrix
    tacs->assembleMatType(GEOMETRIC_STIFFNESS_MATRIX, gmat);

    // Form the shifted operator and factor it
    aux_mat->axpy(sigma, gmat);
    aux_mat->applyBCs(tacs->getBcMap());
  }

  // Factor the preconditioner
  pc->factor();

  // Solve the symmetric eigenvalue problem
  sep->solve(ksm_print);
}

/*!
  Extract the eigenvalue from the analysis.
*/
TacsScalar TACSLinearBuckling::extractEigenvalue( int n,
                                                  TacsScalar *error ){
  return sep->extractEigenvalue(n, error);
}

/*!
  Extract the eigenvector and eigenvalue from the eigenvalue analysis
*/
TacsScalar TACSLinearBuckling::extractEigenvector( int n, TACSBVec *ans,
                                                   TacsScalar *error ){
  return sep->extractEigenvector(n, ans, error);
}

/*
  Return ||I - Q^{T}Q ||_{F}
*/
TacsScalar TACSLinearBuckling::checkOrthogonality(){
  return sep->checkOrthogonality();
}

/*
  Print the components of the matrix Q^{T}Q
*/
void TACSLinearBuckling::printOrthogonality(){
  sep->printOrthogonality();
}

/*!
  Check the actual residual for the given eigenvalue
*/
void TACSLinearBuckling::checkEigenvector( int n ){
  // Test the eignevalue
  TACSBVec *t1 = tacs->createVec();
  TACSBVec *t2 = tacs->createVec();
  t1->incref();
  t2->incref();

  // Extract the eigenvalue and eigenvectors
  TacsScalar eig, error;
  eig = extractEigenvector(n, eigvec, &error);

  // Multiply to get the t
  kmat->mult(eigvec, t1);
  gmat->mult(eigvec, t2);

  // Get the rank so that we only print this once
  int rank;
  MPI_Comm_rank(tacs->getMPIComm(), &rank);

  // Compute the norms of the products
  TacsScalar t1n = t1->norm();
  TacsScalar t2n = t2->norm();

  // Print out the norms of the products K*eigvec and G*eigvec
  if (rank == 0){
    printf("|K*e| = %15.5e  \n|G*e| = %15.5e \n",
           TacsRealPart(t1n), TacsRealPart(t2n));
  }

  // Add the two vectors together and print out the error
  t1->axpy(eig, t2);
  t1n = t1->norm();
  if (rank == 0){
    printf("||K*e + lambda*G*e|| = %15.5e \n", TacsRealPart(t1n));
  }

  // Decref the vectors to free the memory
  t1->decref();
  t2->decref();
}

/*
  The function computes the derivatives of the buckling eigenvalues.

  Compute the derivative of the eignevalues w.r.t. the design
  variables. This function must be called after the solve function has
  been called. The stiffness matrix and geometric stiffness matrix
  cannot be modified from the previous call to solve.

  The original eigenvalue problem is

  K*u + lambda*G*u = 0

  The derivative of the eigenvalue problem is given as follows:

  d(lambda)/dx = - u^{T}*(dK/dx + lambda*dG/dx)*u/(u^{T}*G*u)

  The difficulty is that the load path is determined by solving an
  auxiliary linear system:

  K*path = f

  Since the geometric stiffness matrix is a function of the path, we
  must compute the total derivative of the inner product of the
  geometric stiffness matrix as follows:

  d(u^{T}*G*u)/dx = [ p(u^{T}*G*u)/px - psi*d(K*path)/dx ]

  where the adjoint variables psi are found by solving the linear
  system:

  K*psi = d(u^{T}*G*u)/d(path)
*/
void TACSLinearBuckling::evalEigenDVSens( int n,
                                          TacsScalar fdvSens[],
                                          int numDVs ){
  // Zero the derivative
  memset(fdvSens, 0, numDVs*sizeof(TacsScalar));

  // Copy over the values of the stiffness matrix, factor
  // the stiffness matrix.
  aux_mat->copyValues(kmat);
  pc->factor();

  // Get the eigenvalue and eigenvector
  TacsScalar error;
  TacsScalar eig = extractEigenvector(n, eigvec, &error);

  // Evaluate the partial derivative for the stiffness matrix
  tacs->addMatDVSensInnerProduct(1.0, STIFFNESS_MATRIX,
                                 eigvec, eigvec, fdvSens, numDVs);
  int mpi_rank;
  MPI_Comm_rank(tacs->getMPIComm(),&mpi_rank);

  // Evaluate the derivative of the geometric stiffness matrix
  tacs->addMatDVSensInnerProduct(TacsRealPart(eig), GEOMETRIC_STIFFNESS_MATRIX,
                                 eigvec, eigvec, fdvSens, numDVs);

  // Evaluate derivative of the inner product with respect to
  // the path variables
  tacs->evalMatSVSensInnerProduct(GEOMETRIC_STIFFNESS_MATRIX,
                                  eigvec, eigvec, res);

  // Solve for the adjoint vector and evaluate the derivative of
  // the adjoint-residual inner product
  solver->solve(res, update);
  tacs->addAdjointResProducts(-TacsRealPart(eig), &update, 1,
                              fdvSens, numDVs);

  // Now compute the inner product: u^{T}*G*u
  gmat->mult(eigvec, res);
  TacsScalar scale = res->dot(eigvec);

  // Prepare to scale the final result
  scale = -1.0/scale;

  // Scale the gradient to complete the calculation
  for ( int i = 0; i < numDVs; i++ ){
    fdvSens[i] *= scale;
  }
}

/*!
  The following code computes the eigenvalues and eigenvectors
  for the natural frequency eigenproblem:

  K u = lambda M u

  The method uses a shift and invert strategy in conjunction with the
  Lanczos method with full orthogonalization.

  Input:
  tacs:        the TACS assembler object
  sigma:       the initial value of the shift-invert
  mmat:        the mass matrix object
  kmat:        the stiffness matrix object
  solver:      the Krylov subspace method associated with the kmat
  max_lanczos: the maximum number of Lanczos vectors to use
  num_eigvals: the number of eigenvalues to use
  eig_tol:     the eigenproblem tolerance
*/
TACSFrequencyAnalysis::TACSFrequencyAnalysis( TACSAssembler *_tacs,
                                              TacsScalar _sigma,
                                              TACSMat *_mmat,
                                              TACSMat *_kmat,
                                              TACSKsm *_solver,
                                              int max_lanczos,
                                              int num_eigvals,
                                              double eig_tol ){
  // Store the TACSAssembler pointer
  tacs = _tacs;
  tacs->incref();

  // Set the shift value
  sigma = _sigma;

  // Store the stiffness/mass matrices
  mmat = _mmat;
  kmat = _kmat;
  mmat->incref();
  kmat->incref();

  // Store the pointer to the KSM solver and ensure that the solver
  // is associated with the stiffness matrix object
  TACSMat *mat;
  solver = _solver;
  solver->incref();
  solver->getOperators(&mat, &pc);

  if (mat != kmat){
    fprintf(stderr, "Error, solver must be associated with the \
stiffness matrix\n");
  }

  // Check if the preconditioner is actually a multigrid object. If
  // so, then we have to allocate extra data to store things for each
  // multigrid level.
  mg = dynamic_cast<TACSMg*>(pc);

  // Allocate vectors that are required for the eigenproblem
  eigvec = tacs->createVec();
  res = tacs->createVec();
  eigvec->incref();
  res->incref();

  // Allocate the eigenproblem operator
  ep_op = new EPGeneralizedShiftInvert(sigma, solver, mmat);
  ep_op->incref();

  // Allocate the symmetric eigenproblem solver
  sep = new SEP(ep_op, max_lanczos, SEP::FULL, tacs->getBcMap());
  sep->incref();
  sep->setTolerances(eig_tol, SEP::SMALLEST_MAGNITUDE,
                     num_eigvals);
}

/*
  Deallocate all of the stored data
*/
TACSFrequencyAnalysis::~TACSFrequencyAnalysis(){
  tacs->decref();
  eigvec->decref();
  res->decref();
  solver->decref();
  ep_op->decref();
  sep->decref();
}

/*
  Retrieve the value of sigma
*/
TacsScalar TACSFrequencyAnalysis::getSigma(){
  return sigma;
}

/*
  Set the eigenvalue estimate.
*/
void TACSFrequencyAnalysis::setSigma( TacsScalar _sigma ){
  sigma = _sigma;
  ep_op->setSigma(_sigma);
}

/*
  Solve the eigenvalue problem
*/
void TACSFrequencyAnalysis::solve( KSMPrint *ksm_print ){
  // Zero the variables
  tacs->zeroVariables();

  if (mg){
    // Assemble the mass matrix
    ElementMatrixType matTypes[2] = {STIFFNESS_MATRIX, MASS_MATRIX};
    TacsScalar scale[2] = {1.0, -sigma};

    // Assemble the mass matrix
    tacs->assembleMatType(MASS_MATRIX, mmat);

    // Assemble the linear combination
    mg->assembleMatCombo(matTypes, scale, 2);
  }
  else {
    // Assemble the stiffness and mass matrices
    tacs->assembleMatType(STIFFNESS_MATRIX, kmat);
    tacs->assembleMatType(MASS_MATRIX, mmat);

    // Form the shifted operator and factor it
    kmat->axpy(-sigma, mmat);
    kmat->applyBCs(tacs->getBcMap());
  }

  // Factor the preconditioner
  pc->factor();

  // Solve the symmetric eigenvalue problem
  sep->solve(ksm_print);
}

/*!
  Extract the eigenvalue from the analysis
*/
TacsScalar TACSFrequencyAnalysis::extractEigenvalue( int n,
                                                     TacsScalar *error ){
  return sep->extractEigenvalue(n, error);
}

/*!
  Extract the eigenvector and eigenvalue from the eigenvalue analysis
*/
TacsScalar TACSFrequencyAnalysis::extractEigenvector( int n, TACSBVec *ans,
                                                      TacsScalar *error ){
  return sep->extractEigenvector(n, ans, error);
}

TacsScalar TACSFrequencyAnalysis::checkOrthogonality(){
  return sep->checkOrthogonality();
}

/*!
  Compute the derivative of the eignevalues w.r.t. the design variables

  The original eigenvalue problem is,

  K*u = lambda*M*u

  The derivative of the eigenvalue problem is given as follows,

  dK/dx*u + K*du/dx =
  d(lambda)/dx*M*u + lambda*dM/dx*u + lambda*M*du/dx

  Since M = M^{T} and K = K^{T}, pre-multiplying by u^{T} gives,

  u^{T}*dK/dx*u = d(lambda)/dx*(u^{T}*M*u) + lambda*u^{T}*dM/dx*u

  Rearranging gives,

  (u^{T}*M*u)*d(lambda)/dx = u^{T}*(dK/dx - lambda*dM/dx)*u
*/
void TACSFrequencyAnalysis::evalEigenDVSens( int n,
                                             TacsScalar fdvSens[],
                                             int numDVs ){
  // Zero the derivative
  memset(fdvSens, 0, numDVs*sizeof(TacsScalar));

  // Extract the eigenvalue and eigenvector
  TacsScalar error;
  TacsScalar eig = extractEigenvector(n, eigvec, &error);

  // Evaluate the partial derivative for the stiffness matrix
  tacs->addMatDVSensInnerProduct(1.0, STIFFNESS_MATRIX,
                                 eigvec, eigvec, fdvSens, numDVs);

  // Evaluate the derivative of the geometric stiffness matrix
  tacs->addMatDVSensInnerProduct(-TacsRealPart(eig), MASS_MATRIX,
                                 eigvec, eigvec, fdvSens, numDVs);

  // Finish computing the derivative
  mmat->mult(eigvec, res);
  TacsScalar scale = 1.0/res->dot(eigvec);

  // Finish computing the derivative
  for ( int i = 0; i < numDVs; i++ ){
    fdvSens[i] *= scale;
  }
}

/*!
  Check the actual residual for the given eigenvalue
*/
void TACSFrequencyAnalysis::checkEigenvector( int n ){
  // Assemble the stiffness/mass matrices
  tacs->zeroVariables();
  tacs->assembleMatType(STIFFNESS_MATRIX, kmat);
  tacs->assembleMatType(MASS_MATRIX, mmat);

  // Create temporary arrays required
  TACSBVec *t1 = tacs->createVec();
  TACSBVec *t2 = tacs->createVec();
  t1->incref();
  t2->incref();

  // Extract the eigenvalue and eigenvectors
  TacsScalar eig, error;
  eig = extractEigenvector(n, eigvec, &error);

  // Multiply to get the t
  kmat->mult(eigvec, t1);
  mmat->mult(eigvec, t2);

  // Print out the norms of the products K*eigvec and G*eigvec
  printf("|K*e| = %15.5e  \n|M*e| = %15.5e \n",
         TacsRealPart(t1->norm()), TacsRealPart(t2->norm()));

  // Add the two vectors together and print out the error
  t1->axpy(-eig, t2);
  printf("||K*e - lambda*M*e|| = %15.5e \n", TacsRealPart(t1->norm()));

  // Decref the vectors to free the memory
  t1->decref();
  t2->decref();
}


/*
  Constraint
*/
TACSFrequencyConstraint::TACSFrequencyConstraint( TACSAssembler *_tacs,
                                                  TacsScalar _min_lambda,
                                                  TACSMat *_mmat,
                                                  TACSMat *_kmat,
                                                  TACSPc *_pc ){
  tacs = _tacs;
  tacs->incref();

  // Set the matrices
  mmat = _mmat;
  kmat = _kmat;
  pc = _pc;

  min_lambda = _min_lambda;

  // The maximum size of the Jacobi-Davidson subspace
  max_jd_size = 25;

  // The maximum size of the deflation space
  max_eigen_vectors = 10;

  // The eigen tolerance
  eigtol = 1e-6;

  // Allocate space for the vectors
  V = new TACSBVec*[ max_jd_size+1 ];
  Q = new TACSBVec*[ max_eigen_vectors+1 ];

  // The maximum number of gmres iterations
  max_gmres_size = 30;

  // The residual tolerances for the GMRES iterations
  rtol = 1e-3;
  atol = 1e-30;

  // Create and populate the pointer into the H columns
  Hptr = new int[ max_gmres_size+1 ];
  Hptr[0] = 0;
  for ( int i = 0; i < max_gmres_size; i++ ){
    Hptr[i+1] = Hptr[i] + i+2;
  }

  // Allocate the Hessenberg matrix
  H = new TacsScalar[ Hptr[max_gmres_size] ];
  res = new TacsScalar[ max_gmres_size+1 ];

  // Allocate space for the unitary matrix Q
  Qsin = new TacsScalar[ max_gmres_size ];
  Qcos = new TacsScalar[ max_gmres_size ];

  // Allocate space for the vectors
  W = new TACSBVec*[ max_gmres_size+1 ];

  Z = NULL;
  if (is_flexible){
    Z = new TACSBVec*[ max_gmres_size ];
    for ( int i = 0; i < max_gmres_size; i++ ){

    }
  }
}

TACSFrequencyConstraint::~TACSFrequencyConstraint(){
  // Free the data to solve the GMRES problem
  delete [] H;
  delete [] Hptr;
  delete [] res;
  delete [] Qsin;
  delete [] Qcos;

  // Free the


  // Free the GMRES subspace vectors
  W[0]->decref();
  for ( int i = 0; i < max_gmres_size; i++ ){
    W[i+1]->decref();
  }
  delete [] W;

  if (Z){
    for ( int i = 0; i < max_gmres_size; i++ ){
      Z[i]->decref();
    }
    delete [] Z;
  }
}

void TACSFrequencyConstraint::mult( TACSBVec *x, TACSBVec *y ){
  kmat->mult(x, y); // mmat->mult(x, y);
}

void TACSFrequencyConstraint::applyPc( TACSBVec *x, TACSBVec *y ){
  pc->applyFactor(x, y);
}

void TACSFrequencyConstraint::solve(){
  // Keep track of the current subspace
  V[0]->setRand(-1.0, 1.0);

  // The maximum size of the Jacobi--Davidson subspace
  const int m = max_jd_size;
  memset(M, 0, m*m*sizeof(TacsScalar));

  // Allocate the space for the real work vector
  int lwork = 16*max_jd_size;
  double *rwork = new double[ lwork ];

  // Store the number of converged eigenvalues
  int nconverged = 0;

  for ( int k = 0; k < m; k++ ){
    // Orthogonalize V[k] against all other vectors in the current
    // solution subspace
    for ( int i = 0; i < k; i++ ){
      TacsScalar h = V[k]->dot(V[i]);
      V[k]->axpy(-h, V[i]);
    }

    // Normalize the vector so that it is orthonormal
    TacsScalar vnorm = V[k]->norm();
    V[k]->scale(1.0/vnorm);

    // Compute work = A*V[k]
    mult(V[k], work);

    // Complete the entries in the symmetric matrix M that is formed
    // by M = V^{T}*A*V
    for ( int i = 0; i <= k; i++ ){
      M[k*m + i] = V[i]->dot(work);
      M[i*m + k] = M[k*m + i];
    }

    // Compute the eigenvalues/eigenvectors of the M matrix. Copy over
    // the values from the M matrix into the eigvecs array.
    for ( int j = 0; j <= k; j++ ){
      for ( int i = 0; i <= k; i++ ){
        eigvecs[i + (k+1)*j] = TacsRealPart(M[i + m*j]);
      }
    }

    // The input arguments required for LAPACK
    const char *jobz = "V", *uplo = "U";

    // Compute the eigenvalues of a symmetric matrix
    int info;
    int n = k+1, ldm = k+1;
    LAPACKdsyev(jobz, uplo, &n, eigvecs, &ldm, eigvals,
                rwork, &lwork, &info);

    // Assemble the th predicted eigenvector
    Q[nconverged]->zeroEntries();
    for ( int j = 0; j <= k; j++ ){
      Q[nconverged]->axpy(eigvecs[j], V[j]);
    }
    TacsScalar qnorm = Q[nconverged]->norm();
    Q[nconverged]->scale(1.0/qnorm);

    // Compute the residual and store it in the work vector
    mult(Q[nconverged], work);
    work->axpy(-eigvals[0], Q[nconverged]);

    // Compute the norm of the eigenvalue to check if it has converged
    if (work->norm() < eigtol){
      // This eigenvalue has converged, store it in Q[nconverged] and
      // modify the remaining entries of V. Now we need to compute the
      // new starting vector for the next eigenvalue
      work->zeroEntries();
      for ( int j = 0; j <= k; j++ ){
        work->axpy(eigvecs[k+1 + j], V[j]);
      }

      // Orthogonalize the new starting vector against all other
      // converged eigenvalues
      for ( int j = 0; j < nconverged; j++ ){
        TacsScalar h = work->dot(Q[j]);
        work->axpy(-h, Q[j]);
      }
      nconverged++;

      // Check if we should quit the loop or continue
      if (nconverged >= max_eigen_vectors){
        break;
      }

      // Reset the iteration loop and continue
      k = 0;
      V[0]->copyValues(work);
      continue;
    }

    // Now solve the system (K - min_lambda*M)*t = -work
    // Keep track of the number of iterations in GMRES
    int niters = 0;

    // Copy the residual to the first work vector
    W[0]->copyValues(work);
    res[0] = W[0]->norm();
    W[0]->scale(1.0/res[0]); // W[0] = b/|| b ||

    // Keep track of the initial norm of the right-hand-side
    double beta = TacsRealPart(res[0]);

    // Using GMRES, solve for the update equation
    for ( int i = 0; i < max_gmres_size; i++ ){
      if (is_flexible){
        // Apply the preconditioner, Z[i] = M^{-1} W[i]
        applyPc(W[i], Z[i]);
        mult(Z[i], W[i+1]); // W[i+1] = A*Z[i] = A*M^{-1}*W[i]
      }
      else {
        // W[i+1] = A*work = A*M^{-1}*W[i]
        applyPc(W[i], work);
        mult(work, W[i+1]);
      }

      // First, orthonormalize W[i+1] with respect to the basis from
      // the deflation sub-space. This ensures that the solution
      // remains orthogonal to the basis vectors Q[j] which are
      // orthogonal.
      for ( int j = 0; j <= nconverged; j++ ){
        TacsScalar h = W[i+1]->dot(Q[j]); // h = dot( W[i+1], C[j] )
        W[i+1]->axpy(-h, Q[j]);           // W[i+1] = W[i+1] - h*C[j]
      }

      // Build the orthogonal basis using MGS
      for ( int j = i; j >= 0; j-- ){
        H[j + Hptr[i]] = W[i+1]->dot(W[j]); // H[j,i] = dot( W[i+1], W[i] )
        W[i+1]->axpy(-H[j + Hptr[i]], W[j]); // W[i+1] = W[i+1] - H[j,i]*W[j]
      }

      // Complete the basis
      H[i+1 + Hptr[i]] = W[i+1]->norm(); // H[i+1,i] = || W[i+1] ||
      W[i+1]->scale(1.0/H[i+1 + Hptr[i]]); // W[i+1] = W[i+1]/|| W[i+1] ||

      // Apply the existing part of Q to the new components of the
      // Hessenberg matrix
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

      // Evaluate the sin/cos of the rotation matrix
      Qcos[i] = h1/sq;
      Qsin[i] = h2/sq;
      H[i   + Hptr[i]] =  h1*Qcos[i] + h2*Qsin[i];
      H[i+1 + Hptr[i]] = -h1*Qsin[i] + h2*Qcos[i];

      // Update the residual
      h1 = res[i];
      res[i]   =   h1*Qcos[i];
      res[i+1] = - h1*Qsin[i];

      // if (monitor){
      //   monitor->printResidual(mat_iters, fabs(TacsRealPart(res[i+1])));
      // }

      niters++;

      if (fabs(TacsRealPart(res[i+1])) < atol ||
          fabs(TacsRealPart(res[i+1])) < rtol*beta){
        break;
      }
    }

    // Now, compute the solution. The linear combination of the
    // Arnoldi vectors. The matrix H is now upper triangular. First
    // compute the weights for each basis vector.
    for ( int i = niters-1; i >= 0; i-- ){
      for ( int j = i+1; j < niters; j++ ){
        res[i] = res[i] - H[i + Hptr[j]]*res[j];
      }
      res[i] = res[i]/H[i + Hptr[i]];
    }

    // Zero the next basis vector for the outer Jacobi--Davidson basis
    V[k+1]->zeroEntries();

    // Build the solution based on the minimum residual solution. Note
    // that a factor of -1.0 is applied here since we should solve (K
    // - min_lambda*M)*t = -r = -(A*x - lambda*x). Since the factor is
    // not included above, we include it here.
    double fact = -1.0;
    if (is_flexible){
      for ( int i = 0; i < niters; i++ ){
        V[k+1]->axpy(fact*res[i], Z[i]);
      }
    }
    else {
      // Apply M^{-1} to the linear combination
      work->zeroEntries();
      for ( int i = 0; i < niters; i++ ){
        work->axpy(fact*res[i], W[i]);
      }
      applyPc(work, V[k+1]);
    }
  }

  delete [] rwork;
}
