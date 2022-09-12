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

#include "TACSMatrixFreeMat.h"

/*
  Implementation of geometric multigrid using TACS
*/

/**
  Set up the TACS multi-grid object with a given number of multi-grid
  levels.

  Each level is initialized with the given successive over relaxation
  (SOR) factor (sor_omega), and the number of relaxations to perform
  at each level. The flag sor_symmetric indicates whether or not to
  use a symmetric SOR iteration.

  @param comm The MPI communicator for the multigrid object
  @param nlevels The number of multigrid levels
  @param sor_omega Default successive over-relaxation factor
  @param sor_iters Default number of SOR iterations at each level
  @param sor_symmetric Default flag for using symmetric SOR
*/
TACSMg::TACSMg(MPI_Comm _comm, int _nlevels, double _sor_omega, int _sor_iters,
               int _sor_symmetric) {
  // Copy over the data
  comm = _comm;
  nlevels = _nlevels;
  sor_omega = _sor_omega;
  sor_iters = _sor_iters;
  sor_symmetric = _sor_symmetric;

  if (nlevels < 2) {
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);
    fprintf(stderr,
            "[%d] TACSMg: Multigrid object cannot "
            "be used less than 2 levels\n",
            mpi_rank);
  }

  // Create the list of pointers to the TACS objects at each level
  assembler = new TACSAssembler *[nlevels];
  memset(assembler, 0, nlevels * sizeof(TACSAssembler *));

  // Number of smoothing operations to be performed at each iteration
  iters = new int[nlevels];
  use_galerkin = new int[nlevels];

  // The solution, residual and right-hand-side at each iteration
  x = new TACSBVec *[nlevels];
  b = new TACSBVec *[nlevels];
  r = new TACSBVec *[nlevels];

  // Initialie the data in the arrays
  for (int i = 0; i < nlevels; i++) {
    iters[i] = 1;         // defaults to one - a V cycle
    use_galerkin[i] = 0;  // Use Galerkin projection
    assembler[i] = NULL;
    r[i] = NULL;
    x[i] = NULL;
    b[i] = NULL;
  }

  // Create the pointers to the matrices
  mat = new TACSMat *[nlevels - 1];

  // should re-assemble this matrix
  interp = new TACSBVecInterp *[nlevels - 1];
  pc = new TACSPc *[nlevels - 1];

  for (int i = 0; i < nlevels - 1; i++) {
    interp[i] = NULL;
    mat[i] = NULL;
    pc[i] = NULL;
  }

  monitor = NULL;
  root_mat = NULL;
  root_pc = NULL;

  // Total time for smoothing and interpolation for each level
  cumulative_level_time = new double[nlevels];
  memset(cumulative_level_time, 0, nlevels * sizeof(double));
}

/**
  Deallocate the data stored internally
*/
TACSMg::~TACSMg() {
  for (int i = 0; i < nlevels; i++) {
    if (assembler[i]) {
      assembler[i]->decref();
    }
    if (r[i]) {
      r[i]->decref();
    }
    if (b[i]) {
      b[i]->decref();
    }
    if (x[i]) {
      x[i]->decref();
    }
  }

  for (int i = 0; i < nlevels - 1; i++) {
    if (mat[i]) {
      mat[i]->decref();
    }
    if (interp[i]) {
      interp[i]->decref();
    }
    if (pc[i]) {
      pc[i]->decref();
    }
  }

  if (monitor) {
    monitor->decref();
  }
  if (root_mat) {
    root_mat->decref();
  }
  if (root_pc) {
    root_pc->decref();
  }

  delete[] assembler;
  delete[] mat;
  delete[] iters;
  delete[] x;
  delete[] r;
  delete[] b;
  delete[] interp;
  delete[] pc;
  delete[] cumulative_level_time;
}

/**
  Set the data for the given multi-grid level.

  This consists of setting the TACS finite-element model, and the
  interpolation operator for this level that interpolates from the
  next coarsest level to this one. If this is the coarsest mesh, no
  interpolation operator is needed.

  A custom matrix may be provided. If it is, then the smoother must
  also be provided as well.

  If the use_galerkin flag is used at a given level, the matrix is formed
  from Galerkin projection: A[i] = P[i-1]^{T}*A[i-1]*P[i-1]. This only works
  when the matrix on the finer level is of type TACSParallelMat. If not, then
  the galerkin flag has no effect. Also note, that the use_galerkin flag
  has no effect on the finest level other than to force the use of a
  TACSParallelMat

  @param level The multigrid level
  @param assembler The assembler defined for this level
  @param interp The interpolation object from the coarser level to this one
  @param iters The number of iterations to take at this level
  @param use_galerkin Use Galerkin projection at this level
  @param mat The matrix at this level (may be NULL)
  @param smoother The smoother defined at this level (may be NULL)
*/
void TACSMg::setLevel(int level, TACSAssembler *_assembler,
                      TACSBVecInterp *_interp, int _iters, int _use_galerkin,
                      TACSMat *_mat, TACSPc *_smoother) {
  assembler[level] = _assembler;
  assembler[level]->incref();

  if ((!_mat && _smoother) || (!_smoother && _mat)) {
    fprintf(stderr,
            "TACSMg: You must define both the matrix and preconditioner\n");
    _mat = NULL;
    _smoother = NULL;
  }

  iters[level] = 1;
  if (_iters > 0) {
    iters[level] = _iters;
  }

  // Only define the restriction/interpolation
  // operators for level < nlevels-1
  if (level < nlevels - 1) {
    if (!_interp) {
      fprintf(stderr,
              "TACSMg: Must define prolongation "
              "operators for all but the coarsest problem\n");
    }
    interp[level] = _interp;
    interp[level]->incref();

    // Check to see if we can actually use the Galerkin approach.
    // If not, no matrix or smoother will be defined.
    if (_use_galerkin) {
      if (level == 0) {
        TACSParallelMat *pmat = assembler[level]->createMat();
        mat[level] = pmat;
        mat[level]->incref();

        int zero_guess = 0;
        pc[level] = new TACSGaussSeidel(pmat, zero_guess, sor_omega, sor_iters,
                                        sor_symmetric);
        pc[level]->incref();
      } else {
        TACSParallelMat *fine_mat =
            dynamic_cast<TACSParallelMat *>(mat[level - 1]);
        TACSParallelMat *coarse_mat = NULL;

        if (fine_mat && interp[level - 1]) {
          interp[level - 1]->computeGalerkinNonZeroPattern(fine_mat,
                                                           &coarse_mat);
          mat[level] = coarse_mat;
          mat[level]->incref();

          int zero_guess = 0;
          pc[level] = new TACSGaussSeidel(coarse_mat, zero_guess, sor_omega,
                                          sor_iters, sor_symmetric);
          pc[level]->incref();

          // This matrix is computed using Galerkin projection
          use_galerkin[level] = 1;
        } else {
          fprintf(stderr,
                  "TACSMg: Must define operators from finest to "
                  "to coarsest levels\n");
        }
      }
    }

    // Galerkin matrix was not created. Use either the provided matrix
    // and smoother or create defaults.
    if (!mat[level]) {
      if (_mat) {
        _mat->incref();
        mat[level] = _mat;

        _smoother->incref();
        pc[level] = _smoother;
      } else {
        TACSParallelMat *pmat = assembler[level]->createMat();
        mat[level] = pmat;
        mat[level]->incref();

        // Do not zero the initial guess for the PSOR object
        int zero_guess = 0;
        pc[level] = new TACSGaussSeidel(pmat, zero_guess, sor_omega, sor_iters,
                                        sor_symmetric);
        pc[level]->incref();
      }
    }
  } else {
    // On the coarsest level, use a different approach
    if (_use_galerkin) {
      TACSParallelMat *fine_mat =
          dynamic_cast<TACSParallelMat *>(mat[level - 1]);
      TACSParallelMat *coarse_mat = NULL;

      if (fine_mat && interp[level - 1]) {
        interp[level - 1]->computeGalerkinNonZeroPattern(fine_mat, &coarse_mat);
        root_mat = coarse_mat;
        root_mat->incref();

        root_pc = new TACSBlockCyclicPc(coarse_mat);
        root_pc->incref();

        // Use Galerkin projection to create the coarsest problem
        use_galerkin[level] = 1;
      } else {
        fprintf(stderr,
                "TACSMg: Must define operators from finest to "
                "to coarsest levels\n");
      }
    }

    if (!root_mat) {
      if (_mat) {
        _mat->incref();
        root_mat = _mat;

        _smoother->incref();
        root_pc = _smoother;
      } else {
        // Set up the root matrix
        TACSSchurMat *schur_mat = assembler[level]->createSchurMat();
        root_mat = schur_mat;
        root_mat->incref();

        // Set up the root preconditioner/solver
        int lev = 10000;
        double fill = 15.0;
        int reorder_schur = 1;
        TACSSchurPc *_pc = new TACSSchurPc(schur_mat, lev, fill, reorder_schur);
        root_pc = _pc;
        root_pc->incref();
      }
    }
  }

  if (level > 0) {
    x[level] = assembler[level]->createVec();
    x[level]->incref();

    b[level] = assembler[level]->createVec();
    b[level]->incref();
  }

  r[level] = assembler[level]->createVec();
  r[level]->incref();
}

/**
  Set the state variables for and levels

  @param vec The vector of state variables on the finest mesh
*/
void TACSMg::setVariables(TACSBVec *vec) {
  assembler[0]->setVariables(vec);

  for (int i = 0; i < nlevels - 1; i++) {
    if (i == 0) {
      interp[i]->multWeightTranspose(vec, x[i + 1]);
    } else {
      interp[i]->multWeightTranspose(x[i], x[i + 1]);
    }
    x[i + 1]->applyBCs(assembler[i + 1]->getBcMap());
    assembler[i + 1]->setVariables(x[i + 1]);
  }
}

/**
  Set the design variables for all multi-grid levels.

  This call ensures that all the TACSAssembler objects, and the
  objects they reference, share the same set of design variables.
  This will only work if the same design vector can be shared across
  all multigrid levels.

  @param x The design variable values
*/
void TACSMg::setDesignVars(TACSBVec *x) {
  for (int i = 0; i < nlevels; i++) {
    if (assembler[i]) {
      assembler[i]->setDesignVars(x);
    }
  }
}

/**
  Factor the smoother objects for each level and
  the direct solver for the lowest level.
*/
void TACSMg::factor() {
  for (int i = 0; i < nlevels - 1; i++) {
    if (pc[i]) {
      pc[i]->factor();
    }
  }
  if (root_pc) {
    root_pc->factor();
  }
}

/**
  Compute the Jacobian matrices on each level of the multigrid problem.

  @param alpha The coefficient for the state variable Jacobian
  @param beta The coefficient for the first time derivative Jacobian
  @param gamma The coefficient for the second time derivative Jacobian
  @param res The residual vector, can be NULL
  @param matOr The matrix orientation (either TACS_MAT_NORMAL or
  TACS_MAT_TRANSPOSE)
*/
void TACSMg::assembleJacobian(double alpha, double beta, double gamma,
                              TACSBVec *res, MatrixOrientation matOr) {
  // Assemble the matrices
  if (assembler[0]) {
    TACSMatrixFreeMat *mat_free = dynamic_cast<TACSMatrixFreeMat *>(mat[0]);
    if (mat_free) {
      mat_free->assembleMatrixFreeData(TACS_JACOBIAN_MATRIX, alpha, beta,
                                       gamma);
    } else {
      assembler[0]->assembleJacobian(alpha, beta, gamma, res, mat[0], matOr);
    }
  }

  // Compute the number
  for (int i = 1; i < nlevels - 1; i++) {
    if (use_galerkin[i]) {
      TACSParallelMat *fine_mat = dynamic_cast<TACSParallelMat *>(mat[i - 1]);
      TACSParallelMat *coarse_mat = dynamic_cast<TACSParallelMat *>(mat[i]);
      if (fine_mat && coarse_mat) {
        interp[i - 1]->computeGalerkin(fine_mat, coarse_mat);
        assembler[i]->applyBCs(coarse_mat);
      }
    } else {
      TACSMatrixFreeMat *mat_free = dynamic_cast<TACSMatrixFreeMat *>(mat[i]);
      if (mat_free) {
        mat_free->assembleMatrixFreeData(TACS_JACOBIAN_MATRIX, alpha, beta,
                                         gamma);
      } else if (assembler[i]) {
        assembler[i]->assembleJacobian(alpha, beta, gamma, NULL, mat[i], matOr);
      }
    }
  }

  // Assemble the coarsest problem
  if (use_galerkin[nlevels - 1]) {
    TACSParallelMat *fine_mat =
        dynamic_cast<TACSParallelMat *>(mat[nlevels - 2]);
    TACSParallelMat *coarse_mat = dynamic_cast<TACSParallelMat *>(root_mat);
    if (fine_mat && coarse_mat) {
      interp[nlevels - 2]->computeGalerkin(fine_mat, coarse_mat);
      assembler[nlevels - 1]->applyBCs(coarse_mat);
    }
  } else if (assembler[nlevels - 1]) {
    assembler[nlevels - 1]->assembleJacobian(alpha, beta, gamma, NULL, root_mat,
                                             matOr);
  }
}

/**
  Set up the multigrid data by computing the matrices at each
  multigrid level. This utilizes the Galkerin projection where appropriate.

  @param matType The type of matrix to assemble
  @param matOr The matrix orientation
*/
void TACSMg::assembleMatType(ElementMatrixType matType,
                             MatrixOrientation matOr) {
  if (assembler[0]) {
    TACSMatrixFreeMat *mat_free = dynamic_cast<TACSMatrixFreeMat *>(mat[0]);
    if (mat_free) {
      double alpha = 1.0, beta = 0.0, gamma = 0.0;
      mat_free->assembleMatrixFreeData(matType, alpha, beta, gamma);
    } else {
      assembler[0]->assembleMatType(matType, mat[0], matOr);
    }
  }

  // Compute the number
  for (int i = 1; i < nlevels - 1; i++) {
    if (use_galerkin[i]) {
      TACSParallelMat *fine_mat = dynamic_cast<TACSParallelMat *>(mat[i - 1]);
      TACSParallelMat *coarse_mat = dynamic_cast<TACSParallelMat *>(mat[i]);
      if (fine_mat && coarse_mat) {
        interp[i - 1]->computeGalerkin(fine_mat, coarse_mat);
        assembler[i]->applyBCs(coarse_mat);
      }
    } else {
      TACSMatrixFreeMat *mat_free = dynamic_cast<TACSMatrixFreeMat *>(mat[i]);
      if (mat_free) {
        double alpha = 1.0, beta = 0.0, gamma = 0.0;
        mat_free->assembleMatrixFreeData(matType, alpha, beta, gamma);
      } else if (assembler[i]) {
        assembler[i]->assembleMatType(matType, mat[i], matOr);
      }
    }
  }

  // Assemble the coarsest problem
  if (use_galerkin[nlevels - 1]) {
    TACSParallelMat *fine_mat =
        dynamic_cast<TACSParallelMat *>(mat[nlevels - 2]);
    TACSParallelMat *coarse_mat = dynamic_cast<TACSParallelMat *>(root_mat);
    if (fine_mat && coarse_mat) {
      interp[nlevels - 2]->computeGalerkin(fine_mat, coarse_mat);
      assembler[nlevels - 1]->applyBCs(coarse_mat);
    }
  } else if (assembler[nlevels - 1]) {
    assembler[nlevels - 1]->assembleMatType(matType, root_mat, matOr);
  }
}

/**
  Assemble a linear combination of matrices with the specified linear
  combination

  @param matTypes The array of the types of matrices
  @param scale The scalar factors to apply to each matrix
  @param nmats The number of matrices in the combination
  @param matOr The orientation of the matrices
*/
void TACSMg::assembleMatCombo(ElementMatrixType matTypes[], TacsScalar scale[],
                              int nmats, MatrixOrientation matOr) {
  if (assembler[0]) {
    assembler[0]->assembleMatCombo(matTypes, scale, nmats, mat[0], matOr);
  }

  // Compute the number
  for (int i = 1; i < nlevels - 1; i++) {
    if (use_galerkin[i]) {
      TACSParallelMat *fine_mat = dynamic_cast<TACSParallelMat *>(mat[i - 1]);
      TACSParallelMat *coarse_mat = dynamic_cast<TACSParallelMat *>(mat[i]);
      if (fine_mat && coarse_mat) {
        interp[i - 1]->computeGalerkin(fine_mat, coarse_mat);
        assembler[i]->applyBCs(coarse_mat);
      }
    } else if (assembler[i]) {
      assembler[i]->assembleMatCombo(matTypes, scale, nmats, mat[i], matOr);
    }
  }

  // Assemble the coarsest problem
  if (use_galerkin[nlevels - 1]) {
    TACSParallelMat *fine_mat =
        dynamic_cast<TACSParallelMat *>(mat[nlevels - 2]);
    TACSParallelMat *coarse_mat = dynamic_cast<TACSParallelMat *>(root_mat);
    if (fine_mat && coarse_mat) {
      interp[nlevels - 2]->computeGalerkin(fine_mat, coarse_mat);
      assembler[nlevels - 1]->applyBCs(coarse_mat);
    }
  } else if (assembler[nlevels - 1]) {
    assembler[nlevels - 1]->assembleMatCombo(matTypes, scale, nmats, root_mat,
                                             matOr);
  }
}

/**
  Assemble the multigrid preconditioner matrix with Galerkin flag, this can be
  used to assemble the matrix for mg when mat gets modified outside the class
  via getMat()

  @return fail flag, true if the Galerkin flag is false
*/
int TACSMg::assembleGalerkinMat() {
  int fail = 0;

  // Compute the number
  for (int i = 1; i < nlevels - 1; i++) {
    if (use_galerkin[i]) {
      TACSParallelMat *fine_mat = dynamic_cast<TACSParallelMat *>(mat[i - 1]);
      TACSParallelMat *coarse_mat = dynamic_cast<TACSParallelMat *>(mat[i]);
      if (fine_mat && coarse_mat) {
        interp[i - 1]->computeGalerkin(fine_mat, coarse_mat);
        assembler[i]->applyBCs(coarse_mat);
      }
    } else if (assembler[i]) {
      fail = 1;
    }
  }

  // Assemble the coarsest problem
  if (use_galerkin[nlevels - 1]) {
    TACSParallelMat *fine_mat =
        dynamic_cast<TACSParallelMat *>(mat[nlevels - 2]);
    TACSParallelMat *coarse_mat = dynamic_cast<TACSParallelMat *>(root_mat);
    if (fine_mat && coarse_mat) {
      interp[nlevels - 2]->computeGalerkin(fine_mat, coarse_mat);
      assembler[nlevels - 1]->applyBCs(coarse_mat);
    }
  } else if (assembler[nlevels - 1]) {
    fail = 1;
  }
  return fail;
}

/**
  Get the matrix operator associated with the finest level

  @param _mat A pointer to the matrix
*/
void TACSMg::getMat(TACSMat **_mat) { *_mat = mat[0]; }

/**
  Retrieve the matrix at the specified multigrid level.

  Return NULL if no such matrix exists.

  @param level The multigrid level
  @return The matrix at the specified level (NULL if invalid)
*/
TACSMat *TACSMg::getMat(int level) {
  if (level >= 0 && level < nlevels - 1) {
    return mat[level];
  }
  return NULL;
}

/**
  Retrieve the TACSAssembler at the specified level

  @param level The multgrid level
  @return The TACSAssembler object
*/
TACSAssembler *TACSMg::getAssembler(int level) {
  if (level >= 0 && level < nlevels - 1) {
    return assembler[level];
  }
  return NULL;
}

/**
  Retrieve the TACSBVecInterp object for the given level

  @param level The multgrid level
  @return The TACSBVecInterp object associated with the level
*/
TACSBVecInterp *TACSMg::getInterpolation(int level) {
  if (level >= 0 && level < nlevels - 1) {
    return interp[level];
  }
  return NULL;
}

/**
  Set the monitor to use internally for printing out convergence data.

  @param monitor The print monitor object
*/
void TACSMg::setMonitor(KSMPrint *_monitor) {
  if (_monitor) {
    _monitor->incref();
  }
  if (monitor) {
    monitor->decref();
  }
  monitor = _monitor;
}

/**
  Repeatedly apply the multi-grid method until the problem is solved
*/
void TACSMg::solve(TACSBVec *bvec, TACSBVec *xvec, int max_iters, double rtol,
                   double atol) {
  b[0] = bvec;  // Set the RHS at the finest level
  x[0] = xvec;  // Set the solution at the finest level

  // Compute the initial residual and multiply
  TacsScalar rhs_norm = 0.0;
  for (int i = 0; i < max_iters; i++) {
    applyMg(0);
    TacsScalar norm = r[0]->norm();
    if (monitor) {
      monitor->printResidual(i, norm);
    }
    if (i == 0) {
      rhs_norm = norm;
    }

    if (TacsRealPart(norm) < atol ||
        TacsRealPart(norm) < rtol * TacsRealPart(rhs_norm)) {
      break;
    }
  }

  b[0] = NULL;
  x[0] = NULL;
}

/**
  Apply the multi-grid preconditioner to try and solve the problem
  x = A^{-1} b

  Assume an initial guess of zero.

  @param bvec The input right-hand-side
  @param xvec The output vector
*/
void TACSMg::applyFactor(TACSVec *bvec, TACSVec *xvec) {
  // Set the RHS at the finest level
  b[0] = dynamic_cast<TACSBVec *>(bvec);

  // Set the solution at the finest level
  x[0] = dynamic_cast<TACSBVec *>(xvec);

  if (b[0] && x[0]) {
    x[0]->zeroEntries();
    if (monitor) {
      memset(cumulative_level_time, 0, nlevels * sizeof(double));
    }
    applyMg(0);
    if (monitor) {
      for (int k = 0; k < nlevels; k++) {
        char descript[128];
        sprintf(descript, "TACSMg cumulative level %2d time %15.8e\n", k,
                cumulative_level_time[k]);
        monitor->print(descript);
      }
    }
  } else {
    fprintf(stderr, "TACSMg type error: Input/output must be TACSBVec\n");
  }

  b[0] = NULL;
  x[0] = NULL;
}

/**
  Apply the multigrid algorithm from McCormick and Runge 1983 - but
  with added post-smoothing.

  This function applies multigrid recursively by smoothing the
  residual, restricting to the next level, applying multigrid, then
  post-smoothing.

  @param level Apply a cycle of multigrid at this level
*/
void TACSMg::applyMg(int level) {
  // If we've made it to the lowest level, apply the direct solver
  // otherwise, perform multigrid on the next-lowest level
  if (level == nlevels - 1) {
    double t1 = 0.0;
    if (monitor) {
      t1 = MPI_Wtime();
    }
    root_pc->applyFactor(b[level], x[level]);
    if (monitor) {
      cumulative_level_time[level] += MPI_Wtime() - t1;
    }
    return;
  }

  // Perform iters[level] cycle at the next lowest level
  double t1 = 0.0;
  if (monitor) {
    t1 = MPI_Wtime();
  }
  for (int k = 0; k < iters[level]; k++) {
    // Pre-smooth at the current level
    pc[level]->applyFactor(b[level], x[level]);

    // Compute r[level] = b[level] - A*x[level]
    mat[level]->mult(x[level], r[level]);
    r[level]->axpby(1.0, -1.0, b[level]);

    // Restrict the residual to the next lowest level
    // to form the RHS at that level
    interp[level]->multTranspose(r[level], b[level + 1]);
    b[level + 1]->applyBCs(assembler[level + 1]->getBcMap());
    x[level + 1]->zeroEntries();

    applyMg(level + 1);

    // Interpolate back from the next lowest level
    interp[level]->multAdd(x[level + 1], x[level], x[level]);
    x[level]->applyBCs(assembler[level]->getBcMap());
  }

  // Post-Smooth the residual
  pc[level]->applyFactor(b[level], x[level]);

  if (monitor) {
    cumulative_level_time[level] += MPI_Wtime() - t1;
  }
}
