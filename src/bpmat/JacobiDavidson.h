/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_JACOBI_DAVIDSON_H
#define TACS_JACOBI_DAVIDSON_H

#include "KSM.h"
#include "TACSAssembler.h"

enum JDRecycleType { JD_SUM_TWO, JD_NUM_RECYCLE };

/*
  The following code implements a Jacobi-Davidson method for
  simple and generalized eigenvalue problems.
*/

/*
  This is the abstract base class that are required to define the
  operators needed in the Jacobi-Davidson method.
*/
class TACSJacobiDavidsonOperator : public TACSObject {
 public:
  // Get the MPI_Comm
  virtual MPI_Comm getMPIComm() = 0;

  // Create a vector
  virtual TACSVec *createVec() = 0;

  // Set the eigenvalue estimate (and reset the factorization)
  virtual void setEigenvalueEstimate(double estimate) = 0;

  // Apply the preconditioner
  virtual void applyFactor(TACSVec *x, TACSVec *y) = 0;

  // Apply boundary conditions associated with the matrix
  virtual void applyBCs(TACSVec *x) {}

  // Perform a dot-product of two vectors in the operator space e.g.
  // output <- x^{T}*B*y
  virtual TacsScalar dot(TACSVec *x, TACSVec *y) { return x->dot(y); }

  // Matrix-vector products with either the A or B operators.
  virtual void multA(TACSVec *x, TACSVec *y) = 0;
  virtual void multB(TACSVec *x, TACSVec *y) { y->copyValues(x); }
};

/*
  The Jacobi-Davidson operator for simple eigenvalue problems
*/
class TACSJDSimpleOperator : public TACSJacobiDavidsonOperator {
 public:
  TACSJDSimpleOperator(TACSAssembler *_assembler, TACSMat *_mat, TACSPc *_pc);
  ~TACSJDSimpleOperator();

  // Get the MPI_Comm
  MPI_Comm getMPIComm();

  // Create a vector
  TACSVec *createVec();

  // Set the eigenvalue estimate (and reset the factorization)
  // pc = K in this case
  void setEigenvalueEstimate(double estimate);

  // Apply the preconditioner
  void applyFactor(TACSVec *x, TACSVec *y);

  // Apply boundary conditions associated with the matrix
  void applyBCs(TACSVec *x);

  // Perform a dot-product of two vectors in the operator space e.g.
  // output <- x^{T}*B*y
  TacsScalar dot(TACSVec *x, TACSVec *y);

  // Matrix-vector products with either the A or B operators.
  void multA(TACSVec *x, TACSVec *y);
  void multB(TACSVec *x, TACSVec *y);

 private:
  TACSAssembler *assembler;
  TACSMat *mat;
  TACSPc *pc;
};

/*
  The Jacobi-Davidson operator for the natural frequency problem
*/
class TACSJDFrequencyOperator : public TACSJacobiDavidsonOperator {
 public:
  TACSJDFrequencyOperator(TACSAssembler *_assembler, TACSMat *_kmat,
                          TACSMat *_mmat, TACSMat *_pc_mat, TACSPc *_pc);
  ~TACSJDFrequencyOperator();

  // Get the MPI_Comm
  MPI_Comm getMPIComm();

  // Create a vector
  TACSVec *createVec();

  // Set the eigenvalue estimate (and reset the factorization)
  // pc = K in this case
  void setEigenvalueEstimate(double estimate);

  // Apply the preconditioner
  void applyFactor(TACSVec *x, TACSVec *y);

  // Apply boundary conditions associated with the matrix
  void applyBCs(TACSVec *x);

  // Perform a dot-product of two vectors in the operator space e.g.
  // output <- x^{T}*B*y
  TacsScalar dot(TACSVec *x, TACSVec *y);

  // Matrix-vector products with either the A or B operators.
  void multA(TACSVec *x, TACSVec *y);
  void multB(TACSVec *x, TACSVec *y);

 private:
  TACSAssembler *assembler;
  TACSMat *kmat, *mmat, *pc_mat;
  TACSPc *pc;
  TACSVec *work;
};

/*
  The Jacobi-Davidson solver
*/
class TACSJacobiDavidson : public TACSObject {
 public:
  TACSJacobiDavidson(TACSJacobiDavidsonOperator *oper,
                     int _max_eigen_vectors = 10, int _max_jd_size = 20,
                     int _max_gmres_size = 30);
  ~TACSJacobiDavidson();

  // Get the MPI_Comm
  MPI_Comm getMPIComm();

  // Extract the eigenvalues and eigenvectors
  int getNumConvergedEigenvalues();
  TacsScalar extractEigenvalue(int n, TacsScalar *error);
  TacsScalar extractEigenvector(int n, TACSVec *ans, TacsScalar *error);

  // Solve the eigenvalue problem
  void solve(KSMPrint *ksm_print = NULL, int print_level = 0);

  // Set tolerances to FGMRES
  void setTolerances(double _eig_rtol, double _eig_atol, double _rtol,
                     double _atol);

  // Set the paramter that decide whether the convergence check relies more on
  // eig_atol or more on eig_rtol
  void setThetaCutoff(double _theta_cutoff);

  // Set the number of vectors to recycle
  void setRecycle(int _recycle, JDRecycleType _recycle_type);

 private:
  // The operator class that defines the eigenproblem
  TACSJacobiDavidsonOperator *oper;

  // The recycle flag that indicates the number of converged eigenvectors to
  // recycle in the next solve
  JDRecycleType recycle_type;
  int num_recycle_vecs;

  // Generic work vector
  TACSVec *work;

  // The maximum size of the Jacobi--Davidson subspace
  int max_jd_size;

  // The desired number of eigenvectors
  int max_eigen_vectors;

  // The relative and absolute eigenvalue tolerances
  double eig_rtol, eig_atol;

  // The paramter that decide whether the convergence check relies more on
  // eig_atol or more on eig_rtol
  double theta_cutoff;

  // The matrix of variables
  TacsScalar *M;
  double *ritzvecs, *ritzvals;

  // The eigenvalues
  TacsScalar *eigvals;   // The eigenvalues
  TacsScalar *eigerror;  // The error associated with the eigenvalues
  int *eigindex;  // Argsorted list such that eigvals[argindex[i]] is sorted

  // The vectors for the eigenvalue space
  TACSVec **V;

  // The number of converged eigenvectors/eigenvalues
  int nconverged;

  // The vectors for the deflation space
  TACSVec **Q, **P;

  // The data for the FGMRES code
  int max_gmres_size;
  double rtol, atol;

  // Data for the Hessenberg matrix
  int *Hptr;
  TacsScalar *H;
  TacsScalar *Qcos, *Qsin;
  TacsScalar *res;

  // Subspace data for GMRES whether its flexible or not
  TACSVec **W, **Z;
};

#endif  // TACS_JACOBI_DAVIDSON_H
