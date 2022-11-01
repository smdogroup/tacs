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
  assembler:    The TACS model corresponding to the analysis problem
  sigma:        The spectral shift
  gmat:         The geometric stiffness matrix
  kmat:         The stiffness matrix
  aux_mat:      The auxiliary matrix associated with the solver
  solver:       Whatever KSM object you create
  max_lanczos:  Maximum size of the projected subspace
  num_eigvals:  Number of converged eigenvalues required
  eig_tol:      Tolerance of the eigenvalues
*/
TACSLinearBuckling::TACSLinearBuckling(TACSAssembler *_assembler,
                                       TacsScalar _sigma, TACSMat *_gmat,
                                       TACSMat *_kmat, TACSMat *_aux_mat,
                                       TACSKsm *_solver, int _max_lanczos_vecs,
                                       int _num_eigvals, double _eig_tol) {
  // Copy pointer to the TACS assembler object
  assembler = _assembler;
  assembler->incref();

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
  if (mat != aux_mat) {
    fprintf(stderr,
            "TACSBuckling: Error, solver must be associated with the "
            "auxiliary matrix\n");
  }

  // Check that the auxiliary matrix, geometric stiffness and stiffness
  // matrices are different objects
  if (aux_mat == kmat) {
    fprintf(stderr,
            "TACSBuckling: Error, stiffness and auxiliary matrices "
            "must be different instances\n");
  }
  if (aux_mat == gmat) {
    fprintf(stderr,
            "TACSBuckling: Error, geometric stiffness and auxiliary "
            "matrices must be different instances\n");
  }
  if (gmat == kmat) {
    fprintf(stderr,
            "TACSBuckling: Error, geometric stiffness and stiffness "
            "matrices must be different instances\n");
  }

  // Check if the preconditioner is actually a multigrid object. If
  // so, then we have to allocate extra data to store things for each
  // multigrid level.
  mg = dynamic_cast<TACSMg *>(pc);

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
  sep = new SEP(ep_op, max_lanczos_vecs, SEP::FULL, assembler->getBcMap());
  sep->incref();
  sep->setTolerances(eig_tol, SEP::SMALLEST_MAGNITUDE, num_eigvals);

  // Allocate temporary local vectors
  res = assembler->createVec();
  update = assembler->createVec();
  eigvec = assembler->createVec();
  path = assembler->createVec();
  res->incref();
  update->incref();
  eigvec->incref();
  path->incref();
}

/*
  Destructor object for the buckling object
*/
TACSLinearBuckling::~TACSLinearBuckling() {
  // Dereference the matrix objects
  aux_mat->decref();
  gmat->decref();
  kmat->decref();

  // Dereference the solver/tacs
  assembler->decref();
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
TacsScalar TACSLinearBuckling::getSigma() { return sigma; }

/*
  Set a different value for the shift parameter
*/
void TACSLinearBuckling::setSigma(TacsScalar _sigma) {
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
void TACSLinearBuckling::solve(TACSVec *rhs, KSMPrint *ksm_print) {
  // Zero the variables
  assembler->zeroVariables();

  if (mg) {
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    mg->assembleJacobian(alpha, beta, gamma, res);
    mg->factor();

    // Add the right-hand-side due to external forces
    if (rhs) {
      // res->axpy(1.0, rhs);
      assembler->applyBCs(res);
    }

    // Solve for the load path
    solver->solve(res, path);
    path->scale(-1.0);
    assembler->setVariables(path);

    // Assemble the linear combination of the stiffness matrix
    // and geometric stiffness matrix
    ElementMatrixType matTypes[2] = {TACS_STIFFNESS_MATRIX,
                                     TACS_GEOMETRIC_STIFFNESS_MATRIX};
    TacsScalar scale[2] = {1.0, sigma};
    mg->assembleMatCombo(matTypes, scale, 2);

    // Assemble the geometric stiffness matrix and the stiffness matrix itself
    assembler->assembleMatType(TACS_STIFFNESS_MATRIX, kmat);
    assembler->assembleMatType(TACS_GEOMETRIC_STIFFNESS_MATRIX, gmat);
  } else {
    // Compute the stiffness matrix and copy the values to the
    // auxiliary matrix used to solve for the load path.
    assembler->assembleMatType(TACS_STIFFNESS_MATRIX, kmat);
    aux_mat->copyValues(kmat);

    pc->factor();
    assembler->assembleRes(res);

    // If need to add rhs
    if (rhs) {
      res->axpy(-1.0, rhs);
      assembler->applyBCs(res);
    }

    // Solve for the load path and set the variables
    solver->solve(res, path);
    path->scale(-1.0);
    assembler->setVariables(path);

    // Assemble the stiffness and geometric stiffness matrix
    assembler->assembleMatType(TACS_GEOMETRIC_STIFFNESS_MATRIX, gmat);

    // Form the shifted operator and factor it
    aux_mat->axpy(sigma, gmat);
    aux_mat->applyBCs(assembler->getBcMap());
  }

  // Factor the preconditioner
  pc->factor();

  // Solve the symmetric eigenvalue problem
  sep->solve(ksm_print);
}

/*!
  Extract the eigenvalue from the analysis.
*/
TacsScalar TACSLinearBuckling::extractEigenvalue(int n, TacsScalar *error) {
  return sep->extractEigenvalue(n, error);
}

/*!
  Extract the eigenvector and eigenvalue from the eigenvalue analysis
*/
TacsScalar TACSLinearBuckling::extractEigenvector(int n, TACSBVec *ans,
                                                  TacsScalar *error) {
  return sep->extractEigenvector(n, ans, error);
}

/*
  Return ||I - Q^{T}Q ||_{F}
*/
TacsScalar TACSLinearBuckling::checkOrthogonality() {
  return sep->checkOrthogonality();
}

/*
  Print the components of the matrix Q^{T}Q
*/
void TACSLinearBuckling::printOrthogonality() { sep->printOrthogonality(); }

/*!
  Check the actual residual for the given eigenvalue
*/
void TACSLinearBuckling::checkEigenvector(int n) {
  // Test the eignevalue
  TACSBVec *t1 = assembler->createVec();
  TACSBVec *t2 = assembler->createVec();
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
  MPI_Comm_rank(assembler->getMPIComm(), &rank);

  // Compute the norms of the products
  TacsScalar t1n = t1->norm();
  TacsScalar t2n = t2->norm();

  // Print out the norms of the products K*eigvec and G*eigvec
  if (rank == 0) {
    printf("|K*e| = %15.5e  \n|G*e| = %15.5e \n", TacsRealPart(t1n),
           TacsRealPart(t2n));
  }

  // Add the two vectors together and print out the error
  t1->axpy(eig, t2);
  t1n = t1->norm();
  if (rank == 0) {
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
void TACSLinearBuckling::evalEigenDVSens(int n, TACSBVec *dfdx) {
  // Zero the derivative
  dfdx->zeroEntries();

  // Copy over the values of the stiffness matrix, factor
  // the stiffness matrix.
  aux_mat->copyValues(kmat);
  pc->factor();

  // Get the eigenvalue and eigenvector
  TacsScalar error;
  TacsScalar eig = extractEigenvector(n, eigvec, &error);

  // Evaluate the partial derivative for the stiffness matrix
  assembler->addMatDVSensInnerProduct(1.0, TACS_STIFFNESS_MATRIX, eigvec,
                                      eigvec, dfdx);

  // Evaluate the derivative of the geometric stiffness matrix
  assembler->addMatDVSensInnerProduct(
      TacsRealPart(eig), TACS_GEOMETRIC_STIFFNESS_MATRIX, eigvec, eigvec, dfdx);

  // Evaluate derivative of the inner product with respect to
  // the path variables
  assembler->evalMatSVSensInnerProduct(TACS_GEOMETRIC_STIFFNESS_MATRIX, eigvec,
                                       eigvec, res);

  // Solve for the adjoint vector and evaluate the derivative of
  // the adjoint-residual inner product
  solver->solve(res, update);
  assembler->addAdjointResProducts(-TacsRealPart(eig), 1, &update, &dfdx);

  // Now compute the inner product: u^{T}*G*u
  gmat->mult(eigvec, res);
  TacsScalar scale = res->dot(eigvec);

  dfdx->beginSetValues(TACS_ADD_VALUES);
  dfdx->endSetValues(TACS_ADD_VALUES);

  // Scale the final result
  dfdx->scale(-1.0 / scale);
}

/*!
  The following code computes the eigenvalues and eigenvectors
  for the natural frequency eigenproblem:

  K u = lambda M u

  The method uses a shift and invert strategy in conjunction with the
  Lanczos method with full orthogonalization.

  Input:
  assembler:   the TACS assembler object
  sigma:       the initial value of the shift-invert
  mmat:        the mass matrix object
  kmat:        the stiffness matrix object
  solver:      the Krylov subspace method associated with the kmat
  max_lanczos: the maximum number of Lanczos vectors to use
  num_eigvals: the number of eigenvalues to use
  eig_tol:     the eigenproblem tolerance
*/
TACSFrequencyAnalysis::TACSFrequencyAnalysis(TACSAssembler *_assembler,
                                             TacsScalar _sigma, TACSMat *_mmat,
                                             TACSMat *_kmat, TACSKsm *_solver,
                                             int max_lanczos, int num_eigvals,
                                             double eig_tol) {
  // Store the TACSAssembler pointer
  assembler = _assembler;
  assembler->incref();

  // Set the shift value
  sigma = _sigma;

  // Store the stiffness/mass matrices
  mmat = _mmat;
  kmat = _kmat;
  if (mmat) {
    mmat->incref();
  }
  kmat->incref();

  // Store the pointer to the KSM solver and ensure that the solver
  // is associated with the stiffness matrix object
  TACSMat *mat;
  solver = _solver;
  solver->incref();
  solver->getOperators(&mat, &pc);

  if (mat != kmat) {
    fprintf(stderr,
            "Error, solver must be associated with the "
            "stiffness matrix\n");
  }

  // Check if the preconditioner is actually a multigrid object. If
  // so, then we have to allocate extra data to store things for each
  // multigrid level.
  mg = dynamic_cast<TACSMg *>(pc);

  // Allocate vectors that are required for the eigenproblem
  eigvec = assembler->createVec();
  res = assembler->createVec();
  eigvec->incref();
  res->incref();

  // Allocate the eigenproblem operator
  if (mmat) {
    ep_op = new EPGeneralizedShiftInvert(sigma, solver, mmat);
    ep_op->incref();
    simple_ep_op = NULL;
  } else {
    simple_ep_op = new EPShiftInvert(sigma, solver);
    simple_ep_op->incref();
    ep_op = NULL;
  }

  // Allocate the symmetric eigenproblem solver
  if (mmat) {
    sep = new SEP(ep_op, max_lanczos, SEP::FULL, assembler->getBcMap());
  } else {
    sep = new SEP(simple_ep_op, max_lanczos, SEP::FULL, assembler->getBcMap());
  }
  sep->incref();
  sep->setTolerances(eig_tol, SEP::SMALLEST_MAGNITUDE, num_eigvals);

  // Set unallocated objects to NULL
  pcmat = NULL;
  jd_op = NULL;
  jd = NULL;
}

/*!
  The following code computes the eigenvalues and eigenvectors
  for the natural frequency eigenproblem:

  K u = lambda M u

  The method uses the Jacobi-Davidson method

  Input:
  assembler:   the TACS assembler object
  init_eig:    the initial eigenvalue estimate
  mmat:        the mass matrix object
  kmat:        the stiffness matrix object
  pcmat:       the preconditioner matrix object
  max_jd_size: the maximum number of vectors to use
  fgmres_size: the number of FGMRES iteration
  num_eigvals: the number of eigenvalues to use
  eig_tol:     the eigenproblem tolerance
*/
TACSFrequencyAnalysis::TACSFrequencyAnalysis(
    TACSAssembler *_assembler, TacsScalar _init_eig, TACSMat *_mmat,
    TACSMat *_kmat, TACSMat *_pcmat, TACSPc *_pc, int max_jd_size,
    int fgmres_size, int num_eigvals, double eigtol, double eig_rtol,
    double eig_atol, int num_recycle, JDRecycleType recycle_type) {
  // Store the TACSAssembler pointer
  assembler = _assembler;
  assembler->incref();

  // Set the initial eigenvalue estimate
  sigma = _init_eig;

  // Store the stiffness/mass/preconditioner matrices
  mmat = _mmat;
  kmat = _kmat;
  pcmat = _pcmat;
  pc = _pc;

  if (!pcmat) {
    fprintf(stderr,
            "TACSFrequency: Error, the preconditioner matrix associated "
            "with the Jacobi-Davidson method cannot be NULL\n");
  }
  mmat->incref();
  kmat->incref();
  pcmat->incref();

  // Check if the preconditioner is actually a multigrid object. If
  // so, then we have to allocate extra data to store things for each
  // multigrid level.
  mg = dynamic_cast<TACSMg *>(pc);

  // Allocate vectors that are required for the eigenproblem
  eigvec = assembler->createVec();
  res = assembler->createVec();
  eigvec->incref();
  res->incref();

  // Allocate the Jacobi-Davidson operator
  if (mg) {
    jd_op = new TACSJDFrequencyOperator(assembler, kmat, mmat, pcmat, mg);
  } else {
    jd_op = new TACSJDFrequencyOperator(assembler, kmat, mmat, pcmat, pc);
  }
  jd_op->incref();

  // Allocate the Jacobi-Davidson solver
  jd = new TACSJacobiDavidson(jd_op, num_eigvals, max_jd_size, fgmres_size);
  jd->incref();

  // Set unallocated objects to NULL
  ep_op = NULL;
  sep = NULL;
  solver = NULL;

  // Set the tolerance to the Jacobi-Davidson solver
  jd->setTolerances(eigtol, 1e-30, eig_rtol, eig_atol);

  // Set the number of eigenvectors to recycle
  jd->setRecycle(num_recycle, recycle_type);
}

/*
  Deallocate all of the stored data
*/
TACSFrequencyAnalysis::~TACSFrequencyAnalysis() {
  assembler->decref();
  eigvec->decref();
  res->decref();

  if (jd) {
    jd_op->decref();
    jd->decref();
  } else {
    solver->decref();
    sep->decref();
    if (ep_op) {
      ep_op->decref();
    }
    if (simple_ep_op) {
      simple_ep_op->decref();
    }
  }
}

/*
  Retrieve the value of sigma
*/
TacsScalar TACSFrequencyAnalysis::getSigma() { return sigma; }

/*
  Set the eigenvalue estimate.
*/
void TACSFrequencyAnalysis::setSigma(TacsScalar _sigma) {
  sigma = _sigma;
  if (ep_op) {
    ep_op->setSigma(_sigma);
  } else if (jd_op) {
    jd_op->setEigenvalueEstimate(TacsRealPart(_sigma));
  }
}

/*
  Solve the eigenvalue problem
*/
void TACSFrequencyAnalysis::solve(KSMPrint *ksm_print, int print_level) {
  // Zero the variables
  assembler->zeroVariables();
  if (jd) {
    if (mg) {
      // Assemble the mass matrix
      assembler->assembleMatType(TACS_MASS_MATRIX, mmat);
      assembler->assembleMatType(TACS_STIFFNESS_MATRIX, kmat);

      // Assemble the linear combination
      mg->assembleMatType(TACS_STIFFNESS_MATRIX);
      mg->factor();
    } else {
      assembler->assembleMatType(TACS_MASS_MATRIX, mmat);
      assembler->assembleMatType(TACS_STIFFNESS_MATRIX, kmat);
    }

    // Keep track of the computational time
    double t0 = 0.0;
    if (ksm_print && print_level > 0) {
      t0 = MPI_Wtime();
    }

    // Solve the problem using Jacobi-Davidson
    jd->solve(ksm_print, print_level);

    if (ksm_print && print_level > 0) {
      t0 = MPI_Wtime() - t0;

      char line[256];
      sprintf(line, "JD computational time: %15.6f\n", t0);
      ksm_print->print(line);
    }
  } else {
    if (mg) {
      // Assemble the mass matrix
      ElementMatrixType matTypes[2] = {TACS_STIFFNESS_MATRIX, TACS_MASS_MATRIX};
      TacsScalar scale[2] = {1.0, -sigma};

      // Assemble the mass matrix
      if (mmat) {
        assembler->assembleMatType(TACS_MASS_MATRIX, mmat);
      }

      // Assemble the linear combination
      mg->assembleMatCombo(matTypes, scale, 2);
    } else {
      // Assemble the stiffness and mass matrices
      assembler->assembleMatType(TACS_STIFFNESS_MATRIX, kmat);
      if (mmat) {
        assembler->assembleMatType(TACS_MASS_MATRIX, mmat);
      }

      // Form the shifted operator and factor it
      if (mmat) {
        kmat->axpy(-sigma, mmat);
      }
      kmat->applyBCs(assembler->getBcMap());
    }

    // Keep track of the computational time
    double t0 = 0.0;
    if (ksm_print && print_level > 0) {
      t0 = MPI_Wtime();
    }

    // Factor the preconditioner
    pc->factor();

    // Solve the problem using Jacobi-Davidson
    sep->solve(ksm_print);

    if (ksm_print && print_level > 0) {
      t0 = MPI_Wtime() - t0;

      char line[256];
      sprintf(line, "Lanczos computational time: %15.6f\n", t0);
      ksm_print->print(line);
    }
  }
}

/*!
  Extract the eigenvalue from the analysis
*/
TacsScalar TACSFrequencyAnalysis::extractEigenvalue(int n, TacsScalar *error) {
  if (sep) {
    return sep->extractEigenvalue(n, error);
  } else {
    // Error should be NULL unless needed
    return jd->extractEigenvalue(n, error);
  }
}

/*!
  Extract the eigenvector and eigenvalue from the eigenvalue analysis
*/
TacsScalar TACSFrequencyAnalysis::extractEigenvector(int n, TACSBVec *ans,
                                                     TacsScalar *error) {
  if (sep) {
    return sep->extractEigenvector(n, ans, error);
  } else {
    // Error should be NULL unless needed
    return jd->extractEigenvector(n, ans, error);
  }
}

/*
  Check the orthogonality of the Lanczos subspace
*/
TacsScalar TACSFrequencyAnalysis::checkOrthogonality() {
  if (sep) {
    return sep->checkOrthogonality();
  } else {
    fprintf(stderr,
            "TACSFrequency: No orthogonality check for Jacobi-Davidson\n");
    return 0.0;
  }
}

/*!
  Compute the derivative of the eigenvalues w.r.t. the design variables

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
void TACSFrequencyAnalysis::evalEigenDVSens(int n, TACSBVec *dfdx) {
  // Zero the derivative
  dfdx->zeroEntries();

  // Extract the eigenvalue and eigenvector
  TacsScalar error;
  TacsScalar eig = extractEigenvector(n, eigvec, &error);

  // Evaluate the partial derivative for the stiffness matrix
  assembler->addMatDVSensInnerProduct(1.0, TACS_STIFFNESS_MATRIX, eigvec,
                                      eigvec, dfdx);

  // Evaluate the derivative of the geometric stiffness matrix
  assembler->addMatDVSensInnerProduct(-TacsRealPart(eig), TACS_MASS_MATRIX,
                                      eigvec, eigvec, dfdx);

  // Finish computing the derivative
  if (mmat) {
    mmat->mult(eigvec, res);
  }
  TacsScalar scale = 1.0 / res->dot(eigvec);

  dfdx->beginSetValues(TACS_ADD_VALUES);
  dfdx->endSetValues(TACS_ADD_VALUES);
  dfdx->scale(scale);
}

/*!
  Compute the derivative of the eigenvalues w.r.t. the nodal coordinates

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
void TACSFrequencyAnalysis::evalEigenXptSens(int n, TACSBVec *dfdXpt) {
  // Zero the derivative
  dfdXpt->zeroEntries();

  // Extract the eigenvalue and eigenvector
  TacsScalar error;
  TacsScalar eig = extractEigenvector(n, eigvec, &error);

  // Evaluate the partial derivative for the stiffness matrix
  assembler->addMatXptSensInnerProduct(1.0, TACS_STIFFNESS_MATRIX, eigvec,
                                       eigvec, dfdXpt);

  // Evaluate the derivative of the geometric stiffness matrix
  assembler->addMatXptSensInnerProduct(-TacsRealPart(eig), TACS_MASS_MATRIX,
                                       eigvec, eigvec, dfdXpt);

  // Finish computing the derivative
  if (mmat) {
    mmat->mult(eigvec, res);
  }
  TacsScalar scale = 1.0 / res->dot(eigvec);

  dfdXpt->beginSetValues(TACS_ADD_VALUES);
  dfdXpt->endSetValues(TACS_ADD_VALUES);
  dfdXpt->scale(scale);
}

/*!
  Check the actual residual for the given eigenvalue
*/
void TACSFrequencyAnalysis::checkEigenvector(int n) {
  // Assemble the stiffness/mass matrices
  assembler->zeroVariables();
  assembler->assembleMatType(TACS_STIFFNESS_MATRIX, kmat);
  if (mmat) {
    assembler->assembleMatType(TACS_MASS_MATRIX, mmat);
  }

  // Create temporary arrays required
  TACSBVec *t1 = assembler->createVec();
  TACSBVec *t2 = assembler->createVec();
  t1->incref();
  t2->incref();

  // Extract the eigenvalue and eigenvectors
  TacsScalar eig, error;
  eig = extractEigenvector(n, eigvec, &error);

  // Multiply to get the t
  kmat->mult(eigvec, t1);
  if (mmat) {
    mmat->mult(eigvec, t2);
  }

  // Print out the norms of the products K*eigvec and G*eigvec
  printf("|K*e| = %15.5e  \n|M*e| = %15.5e \n", TacsRealPart(t1->norm()),
         TacsRealPart(t2->norm()));

  // Add the two vectors together and print out the error
  t1->axpy(-eig, t2);
  printf("||K*e - lambda*M*e|| = %15.5e \n", TacsRealPart(t1->norm()));

  // Decref the vectors to free the memory
  t1->decref();
  t2->decref();
}
