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

#ifndef TACS_BUCKLING_H
#define TACS_BUCKLING_H

/*
  Implementation of buckling and frequency analysis and sensitivity
  analysis of eigenvalues.
*/

#include "GSEP.h"
#include "JacobiDavidson.h"
#include "TACSAssembler.h"
#include "TACSMg.h"

/*
  Linearized buckling analysis code.

  The efficient solution of generalized eigenvalue problems requires a
  shift and invert operation that involves among other things, copying
  and axpy operations on matrices.  These operations are not supported
  with the TACSMat interface because it would be difficult to do this
  for matrices of different alternate types. This intermediate object
  maintains consistency between matrix types involved in the operation
  without exposing the underlying matrix type.
*/
class TACSLinearBuckling : public TACSObject {
 public:
  TACSLinearBuckling(TACSAssembler *_assembler, TacsScalar _sigma,
                     TACSMat *_gmat, TACSMat *_kmat, TACSMat *_aux_mat,
                     TACSKsm *_solver, int _max_lanczos_vecs, int _num_eigvals,
                     double _eig_tol);
  ~TACSLinearBuckling();

  // Retrieve the instance of TACSAssembler
  // --------------------------------------
  TACSAssembler *getAssembler() { return assembler; }

  // Functions to set the shift value
  // --------------------------------
  TacsScalar getSigma();
  void setSigma(TacsScalar sigma);

  // Solve the eigenvalue problem
  // ----------------------------
  void solve(TACSVec *rhs = NULL, KSMPrint *ksm_print = NULL);
  void evalEigenDVSens(int n, TACSBVec *dfdx);

  // Extract the eigenvalue or check the solution
  // --------------------------------------------
  TacsScalar extractEigenvalue(int n, TacsScalar *error);
  TacsScalar extractEigenvector(int n, TACSBVec *ans, TacsScalar *error);
  void checkEigenvector(int n);
  TacsScalar checkOrthogonality();
  void printOrthogonality();

 private:
  // Data for the eigenvalue analysis
  TacsScalar sigma;

  EPBucklingShiftInvert *ep_op;
  SEP *sep;

  // The tacs object
  TACSAssembler *assembler;

  // Tolerances/required number of eigenvalues
  int max_lanczos_vecs;
  int num_eigvals;
  double eig_tol;

  // These are used by the eigenvalue solver and to solve the linear systems
  // for the path determination problem
  TACSPc *pc;
  TACSKsm *solver;
  TACSMat *aux_mat, *kmat, *gmat;

  // Vectors used in the analysis
  TACSBVec *path;  // The solution path
  TACSBVec *res, *update, *eigvec;

  // The multigrid object -- only defined if a multigrid
  // preconditioner is used
  TACSMg *mg;
};

/*!
  The following class performs frequency analysis and gradient
  evaluation of a TACS finite-element model.

  The code computes eigenvalues and eigenvectors of the generalized
  eigenproblem:

  K u = lambda M u

  The natural frequency of vibration are determined where lambda =
  omega^2.

  The code uses a Lanczos eigenproblem solver with full
  orthogonalization.  The full orthogonalization ensures that the
  Lanczos basis is linearly independent to the required precision. The
  derivatives of the eigenvalues are obtained using an efficient
  method for computing the derivative of the inner product of two
  vectors and the corresponding matrix.
*/
class TACSFrequencyAnalysis : public TACSObject {
 public:
  TACSFrequencyAnalysis(TACSAssembler *_assembler, TacsScalar _sigma,
                        TACSMat *_mmat, TACSMat *_kmat, TACSKsm *_solver,
                        int max_lanczos, int num_eigvals, double _eig_tol);

  TACSFrequencyAnalysis(TACSAssembler *_assembler, TacsScalar _sigma,
                        TACSMat *_mmat, TACSMat *_kmat, TACSMat *_pcmat,
                        TACSPc *_pc, int max_jd_size, int fgmres_size,
                        int num_eigvals, double eigtol = 1e-9,
                        double eig_rtol = 1e-9, double eig_atol = 1e-30,
                        int num_recycle = 0,
                        JDRecycleType recycle_type = JD_NUM_RECYCLE);

  ~TACSFrequencyAnalysis();

  // Retrieve the instance of TACSAssembler
  // --------------------------------------
  TACSAssembler *getAssembler() { return assembler; }

  // Solve the generalized eigenvalue problem
  // ----------------------------------------
  TacsScalar getSigma();
  void setSigma(TacsScalar _sigma);
  void solve(KSMPrint *ksm_print = NULL, int print_level = 0);
  void evalEigenDVSens(int n, TACSBVec *dfdx);
  void evalEigenXptSens(int n, TACSBVec *dfdX);

  // Extract and check the solution
  // ------------------------------
  TacsScalar extractEigenvalue(int n, TacsScalar *error);
  TacsScalar extractEigenvector(int n, TACSBVec *ans, TacsScalar *error);
  void checkEigenvector(int n);
  TacsScalar checkOrthogonality();

 private:
  // The TACS assembler object
  TACSAssembler *assembler;

  // The matrices used in the analysis
  TACSMat *mmat;    // The mass matrix
  TACSMat *kmat;    // The stiffness matrix
  TACSMat *pcmat;   // Matrix associated with the preconditioner (JD)
  TACSKsm *solver;  // Associated with kmat
  TACSPc *pc;       // The preconditioner

  // The multigrid object -- only defined if a multigrid
  // preconditioner is used
  TACSMg *mg;

  // The eigen solver
  TacsScalar sigma;
  EPGeneralizedShiftInvert *ep_op;
  EPShiftInvert *simple_ep_op;
  SEP *sep;

  // Objects associated with the Jacobi-Davidson method
  TACSJDFrequencyOperator *jd_op;
  TACSJacobiDavidson *jd;

  // Vectors required for eigen-sensitivity analysis
  TACSBVec *eigvec, *res;
};

#endif  // TACS_BUCKLING_H
