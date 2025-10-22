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

#ifndef TACS_MG_H
#define TACS_MG_H

/*
  Multigrid code for solution/preconditioning problems.
*/

#include "TACSAssembler.h"
#include "TACSBVec.h"
#include "TACSBVecInterp.h"

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
  TACSMg(MPI_Comm comm, int _nlevels, double _sor_omega = 1.0,
         int _sor_iters = 1, int _sor_symmetric = 0);
  ~TACSMg();

  // Set the data for the multi-grid level
  // -------------------------------------
  void setLevel(int level, TACSAssembler *_assembler,
                TACSBVecInterp *interp = NULL, int _iters = 1,
                int _use_galerkin = 0, TACSMat *_mat = NULL,
                TACSPc *_smoother = NULL);

  // Set the state/design variables of all lower finite-element models
  // -----------------------------------------------------------------
  void setVariables(TACSBVec *vec);
  void setDesignVars(TACSBVec *x);

  // Assemble the given finite-element matrix at all levels
  // ------------------------------------------------------
  void assembleJacobian(double alpha, double beta, double gamma,
                        TACSBVec *res = NULL,
                        MatrixOrientation matOr = TACS_MAT_NORMAL);
  void assembleMatType(ElementMatrixType matType = TACS_STIFFNESS_MATRIX,
                       MatrixOrientation matOr = TACS_MAT_NORMAL);
  void assembleMatCombo(ElementMatrixType matTypes[], TacsScalar scale[],
                        int nmats, MatrixOrientation matOr = TACS_MAT_NORMAL);
  int assembleGalerkinMat();

  // Methods required by the TACSPc class
  // ------------------------------------
  void applyFactor(TACSVec *x, TACSVec *y);
  void factor();

  // Solve the problem using the full multi-grid method
  // --------------------------------------------------
  void solve(TACSBVec *bvec, TACSBVec *xvec, int max_iters = 200,
             double rtol = 1e-8, double atol = 1e-30);

  // Retrieve the matrix from the specified level
  // --------------------------------------------
  void getMat(TACSMat **_mat);
  TACSMat *getMat(int level);
  TACSAssembler *getAssembler(int level);
  TACSBVecInterp *getInterpolation(int level);

  // Set the solution monitor context
  // --------------------------------
  void setMonitor(KSMPrint *_monitor);

 private:
  // Recursive function to apply multi-grid at each level
  void applyMg(int level);

  // The MPI communicator for this object
  MPI_Comm comm;

  // Monitor the solution
  KSMPrint *monitor;

  // The SOR data
  int sor_iters, sor_symmetric;
  double sor_omega;

  // Flag to indicate whether to form the coarse grid operators
  // via Galerkin projection coarse = P^{T}*A*P
  int *use_galerkin;

  // The number of multi-grid levels
  int nlevels;

  // The TACSAssembler object for each level
  TACSAssembler **assembler;
  int *iters;

  // The solution, right-hand-side and residual on each level
  TACSBVec **x, **b, **r;

  // The interpolation operators
  TACSBVecInterp **interp;

  // Time spent on each level
  double *cumulative_level_time;

  // The matrices/preconditioner objects required for multigrid
  TACSMat *root_mat;  // The root matrix
  TACSPc *root_pc;    // The root direct solver
  TACSMat **mat;      // The matrices associated with each level
  TACSPc **pc;        // The smoothers for all but the lowest level
};

#endif  // TACS_MG_H
