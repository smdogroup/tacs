#include "TACSBuckling.h"
#include "tacslapack.h"

/*
  Implementation of the buckling/frequency analysis

  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
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

  // Check that the matrix associated with the solver object
  // is the auxiliary matrix. If not, complain about it.
  if (mat != aux_mat){
    fprintf(stderr, "Error, solver must be associated with the \
auxiliary matrix\n");
  } 
    
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
  tacs->zeroVariables();
  tacs->assembleMatType(STIFFNESS_MATRIX, kmat);

  // Copy values from kmat to aux_mat
  // This requires a pointer to either ScMat or PMat 
  aux_mat->copyValues(kmat);
  
  // Determine the tangent to the solution path at the origin
  pc->factor();
  tacs->assembleRes(res);
  // If need to add rhs
  if (rhs){
    res->axpy(-1.0, rhs);
    tacs->applyBCs(res);
  }
  solver->solve(res, path);
  path->scale(-1.0);

  // Assemble the geometric stiffness matrix
  tacs->setVariables(path);
  tacs->assembleMatType(GEOMETRIC_STIFFNESS_MATRIX, gmat);

  // Set up the eigenvalue problem
  aux_mat->axpy(sigma, gmat);
  aux_mat->applyBCs(tacs->getBcMap());
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

  // All reduce across the processors
  MPI_Allreduce(MPI_IN_PLACE, fdvSens, numDVs, MPI_INT,
                MPI_SUM, tacs->getMPIComm());

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

  dK/dx * u + K * du/dx = 
  d lambda/dx M u + lambda * dM/dx * u + lambda M * du/dx

  Since M = M^{T} and K = K^{T}, pre-multiplying by u^{T} gives,

  u^{T} * dK/dx * u = d lambda/dx ( u^{T} * M * u ) + lambda * u^{T} * dM/dx * u

  Rearranging gives, 

  ( u^{T} * M * u ) [ d lambda/dx ] = u^{T} * ( dK/dx - lambda * dM/dx ) * u
*/
void TACSFrequencyAnalysis::evalEigenDVSens( int n,
                                             TacsScalar fdvSens[], 
                                             int numDVs ){
  // Allocate extra space for the derivative
  TacsScalar *temp = new TacsScalar[ numDVs ];

  // Extract the eigenvalue and eigenvector
  TacsScalar error;
  TacsScalar eig = extractEigenvector(n, eigvec, &error);
  
  // Evaluate the partial derivative for the stiffness matrix
  tacs->addMatDVSensInnerProduct(1.0, STIFFNESS_MATRIX,
                                 eigvec, eigvec, fdvSens, numDVs);

  // Evaluate the derivative of the geometric stiffness matrix
  tacs->addMatDVSensInnerProduct(-TacsRealPart(eig), MASS_MATRIX,
                                 eigvec, eigvec, temp, numDVs);
  
  // Finish computing the derivative
  mmat->mult(eigvec, res);
  TacsScalar scale = 1.0/res->dot(eigvec);

  // Finish computing the derivative
  for ( int i = 0; i < numDVs; i++ ){
    fdvSens[i] *= scale;
  }

  delete [] temp;
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
