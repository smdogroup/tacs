#include "TACSContinuation.h"

/*
  Implementation of buckling and frequency analysis and sensitivity
  analysis of eigenvalues.

  Copyright (c) 2015 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

/**
  The following is a special matrix object that may be used to
  compute the tangent and compute corrector steps for the continuation
  method.

  The objective is to solve an underdetermined system of equations
  with an additional constraint imposed to make the system non-singular.
  The equations are:

  [ A | r ][x]
  .        [y] = [b]    (1)

  With A in R^{nxn} b,r,x in R^{n} and y in R. The full system is nx(n+1).

  The following constraint is imposed on (x,y):

  t^{T} x + s y = 0     (2)

  This constraint is handled by forcing the Krylov iterate to satisfy
  (2) exactly, while approximating solving (1). We compute,

  range(Q) = [t]^{perp}
  .          [s]

  such that [t]^{T} Q = 0
  .         [s]

  where Q in R^{(n+1)xn}. Furthermore ||Qy||_2 = ||y||_{2}.

  An appropriate choice for Q is to use:

  Q[x] = [x - 2v v^{T}x] = [x - 2t^{T}x/(tn*tn) t ]
  .      [  - 2w v^{T}x]   [ -2wn t^{T}x/(tn*tn)  ]

  with v = t/tn, w = wn/tn

  wn = s + sign(s)*sqrt( ||t||^2 + s^2 )
  where tn = sqrt( ||t||^2 + wn^2 )

  The solution to the whole system of equations is,

  [x] = Q[x']
  [y]
*/
TACSContinuationPathMat::TACSContinuationPathMat(TACSMat *_A, TACSVec *_r,
                                                 TACSVec *_t, TacsScalar s) {
  A = _A;
  A->incref();
  r = _r;
  r->incref();
  t = _t;
  t->incref();
  xtmp = A->createVec();
  xtmp->incref();
  resetConstraint(s);
}

TACSContinuationPathMat::~TACSContinuationPathMat() {
  A->decref();
  r->decref();
  t->decref();
  xtmp->decref();
}

void TACSContinuationPathMat::getVectors(TACSVec **_r, TACSVec **_t) {
  if (_r) {
    *_r = r;
  }
  if (_t) {
    *_t = t;
  }
}

void TACSContinuationPathMat::resetConstraint(TacsScalar s) {
  TacsScalar tnorm = t->norm();
  TacsScalar t2 = tnorm * tnorm;

  wn = s + (TacsRealPart(s) >= 0.0 ? 1.0 : -1.0) * sqrt(t2 + s * s);
  tn = sqrt(t2 + wn * wn);
}

// Multiply x <-- Qx, return the value of the n+1-th row
TacsScalar TACSContinuationPathMat::extract(TACSVec *x) {
  TacsScalar tTx = t->dot(x);
  x->axpy(-(2.0 * tTx) / (tn * tn), t);
  return -(2.0 * wn * tTx) / (tn * tn);
}

void TACSContinuationPathMat::getSize(int *_nr, int *_nc) {
  A->getSize(_nr, _nc);
}

void TACSContinuationPathMat::mult(TACSVec *x, TACSVec *y) {
  xtmp->copyValues(x);

  TacsScalar tTx = t->dot(x);
  xtmp->axpy(-(2.0 * tTx) / (tn * tn), t);

  A->mult(xtmp, y);
  y->axpy(-(2.0 * wn * tTx) / (tn * tn), r);
}

TACSVec *TACSContinuationPathMat::createVec() { return A->createVec(); }

/**
  Callback class for monitoring the continuation algorithm
*/
TACSContinuationCallback::TACSContinuationCallback(MPI_Comm _comm,
                                                   const char *filename) {
  comm = _comm;

  int rank;
  MPI_Comm_rank(comm, &rank);

  if (filename && rank == 0) {
    fp = fopen(filename, "w");

    if (fp) {
      fprintf(fp, "Variables = iter, lambda, dlambda_ds\n");
      fprintf(fp, "Zone T=\"Continuation points\"\n");
    }
  }
}

TACSContinuationCallback::~TACSContinuationCallback() {
  if (fp) {
    fclose(fp);
  }
}

void TACSContinuationCallback::iteration(int iter, TACSBVec *vars,
                                         TacsScalar lambda,
                                         TacsScalar dlambda_ds,
                                         TACSAssembler *assembler) {
  if (fp) {
    fprintf(fp, "%2d %15.6e %15.6e\n", iter + 1, TacsRealPart(lambda),
            TacsRealPart(dlambda_ds));
    fflush(fp);
  }
}

/**
  Initialize the arc-length continuation class
*/
TACSContinuation::TACSContinuation(TACSAssembler *_assembler,
                                   int _max_continuation_iters,
                                   int _max_correction_iters,
                                   int _max_correction_restarts,
                                   double _corr_rtol, double _corr_dtol,
                                   double _krylov_rtol, double _krylov_atol,
                                   double _tangent_rtol, double _tangent_atol) {
  assembler = _assembler;
  assembler->incref();

  max_continuation_iters = _max_continuation_iters;
  max_correction_iters = _max_correction_iters;
  max_correction_restarts = _max_correction_restarts;

  correction_rtol = 1e-8;
  if (_corr_rtol < 0.1) {
    correction_rtol = _corr_rtol;
  }
  correction_dtol = 1e3;
  if (_corr_dtol > 10.0) {
    correction_dtol = _corr_dtol;
  }

  krylov_rtol = _krylov_rtol;
  krylov_atol = _krylov_atol;

  tangent_rtol = _tangent_rtol;
  tangent_atol = _tangent_atol;

  iteration_count = 0;
  lambda_history = new TacsScalar[max_continuation_iters];
  dlambda_ds_history = new TacsScalar[max_continuation_iters];

  memset(lambda_history, 0, max_continuation_iters * sizeof(TacsScalar));
  memset(dlambda_ds_history, 0, max_continuation_iters * sizeof(TacsScalar));

  // Set the default termination information
  term_function = NULL;
  term_function_value = 1.0;
  dlambda_ds_term_value = -1e20;
}

TACSContinuation::~TACSContinuation() {
  assembler->decref();

  delete[] lambda_history;
  delete[] dlambda_ds_history;

  if (term_function) {
    term_function->decref();
  }
}

/**
  Set the termination conditions
*/
void TACSContinuation::setTermFunction(TACSFunction *func,
                                       TacsScalar term_value) {
  func->incref();
  if (term_function) {
    term_function->decref();
  }
  term_function = func;
  term_function_value = term_value;
}

void TACSContinuation::setTermLambdaRate(TacsScalar term_dlambda_ds) {
  dlambda_ds_term_value = term_dlambda_ds;
}

/**
  Retrieve information about the solve
*/
int TACSContinuation::getNumIterations() { return iteration_count; }

void TACSContinuation::getSolution(int iter, TacsScalar *lambda,
                                   TacsScalar *dlambda_ds) {
  if (iter >= 0 && iter < iteration_count) {
    *lambda = lambda_history[iter];
    *dlambda_ds = dlambda_ds_history[iter];
  } else {
    *lambda = 0.0;
    *dlambda_ds = 0.0;
  }
}

/**
  Perform arc-length continuation to track the equilibrium path past
  limit points (folds).

  The equations are given by

  R(u, lambda) = r(u) - lambda * load

  Here r(u) are the residuals stored in the TACSAssembler object, load
  is a specified load vector containing the loads and boundary conditions
  from the problem. The load path parameter lambda controls the distance
  along the load path.

  The solution proceeds as follows:

  1. Compute the tangent vector at the current point and compute the
  first order update.

  The solution curve is given by (u,lambda), and represents solutions
  of the following system of equations:

  J * tangent - load * delta_lambda = 0 (1)

  where J = dr/du is the Jacobian matrix from the TACSAssembler object.

  The tangent to the curve at a point is given by tangent = J^{-1} load.

  An additional equation is added to (1) to parameterize the curve
  through the arc-length. This is performed by employing a
  linearization of the arc-length constraint:

  t^{T}(u_{k+1} - u_{k}) + (lambda_{k+1} - lambda_k) = delta_s/(dlambda_ds)

  2. Perform corrector iterations to bring the solution back to the
  equilibrium path. These corrector iterations employ a Newton method
  to bring the point back onto the equilibrium path.
*/
void TACSContinuation::solve_tangent(TACSMat *mat, TACSPc *pc, TACSKsm *ksm,
                                     TACSBVec *load, TacsScalar lambda_init,
                                     TacsScalar target_delta_lambda,
                                     KSMPrint *ksm_print,
                                     TACSContinuationCallback *callback) {
  TACSBVec *vars = assembler->createVec();
  TACSBVec *old_vars = assembler->createVec();
  TACSBVec *temp = assembler->createVec();
  TACSBVec *tangent = assembler->createVec();
  TACSBVec *update = assembler->createVec();
  TACSBVec *res = assembler->createVec();

  vars->incref();
  old_vars->incref();
  temp->incref();
  tangent->incref();
  update->incref();
  res->incref();

  TACSContinuationPathMat *path_mat =
      new TACSContinuationPathMat(mat, load, tangent, 0.0);
  path_mat->incref();

  TacsScalar lambda = 0.0;          // The load factor
  TacsScalar lambda_old = 0.0;      // The previous load factor
  TacsScalar target_delta_r = 0.0;  // Target change in r = (u - u_k)^T(u - u_k)

  double t0 = MPI_Wtime();

  if (lambda_init != 0.0) {
    lambda = lambda_init;

    if (ksm_print) {
      char line[256];
      ksm_print->print("Performing initial Newton iterations\n");
      sprintf(line, "%5s %9s %10s\n", "Iter", "t", "|R|");
      ksm_print->print(line);
    }

    // Set the tolerances for the tangent computation
    ksm->setTolerances(tangent_rtol, tangent_atol);

    // Compute and factor the stiffness matrix
    TacsScalar alpha = 1.0, beta = 0.0, gamma = 0.0;
    assembler->assembleJacobian(alpha, beta, gamma, res, mat);
    pc->factor();

    // Compute the tangent vector
    ksm->solve(load, tangent);

    // Update the variables based on the tangent computation
    vars->axpy(lambda, tangent);

    // Compute the initial norm based on the tangent approximation
    TacsScalar res_norm_init = 1.0;

    // Now perform a Newton iteration until convergence
    for (int k = 0; k < max_continuation_iters; k++) {
      // Set the variables
      assembler->setVariables(vars);

      // Assemble the residuals at the current point
      assembler->assembleRes(res);
      res->axpy(-lambda, load);

      TacsScalar res_norm = res->norm();
      if (ksm_print) {
        char line[256];
        sprintf(line, "%5d %9.4f %10.4e\n", k + 1, MPI_Wtime() - t0,
                TacsRealPart(res_norm));
        ksm_print->print(line);
      }
      if (k == 0) {
        res_norm_init = res_norm;
      }
      if (TacsRealPart(res_norm) <
          correction_rtol * TacsRealPart(res_norm_init)) {
        break;
      }

      // Solve for the update and update the state variables
      ksm->solve(res, update);
      vars->axpy(-1.0, update);
    }
  }

  if (ksm_print) {
    ksm_print->print("Beginning arc-length continuation method\n");
  }

  // The rate of change of the lambda w.r.t. the arc-length
  TacsScalar dlambda_ds = 0.0;

  for (iteration_count = 0; iteration_count < max_continuation_iters;
       iteration_count++) {
    // Copy the current values to a vector
    assembler->setVariables(vars);

    // Assemble the stiffness matrix at the current iteration
    TacsScalar alpha = 1.0, beta = 0.0, gamma = 0.0;
    assembler->assembleJacobian(alpha, beta, gamma, res, mat);

    // Compute the residual as r(u, lambda) = R(u) - lambda*load
    res->axpy(-lambda, load);

    // Factor the preconditioner
    pc->factor();

    // Compute the change in r = (u - u_k)^{T}(u - u_k)
    TacsScalar delta_s = 1.0;  // this will be over-written later

    // Set the tolerances for the tangent computation
    ksm->setTolerances(tangent_rtol, tangent_atol);

    if (iteration_count == 0) {
      // Compute the initial tangent vector to the solution path
      ksm->setOperators(mat, pc);
      ksm->solve(load, tangent);

      TacsScalar tnorm = tangent->norm();
      dlambda_ds = 1.0 / sqrt(1.0 + tnorm * tnorm);

      tangent->scale(dlambda_ds);
    } else {
      // Set the ksm to use the path_mat object
      ksm->setOperators(path_mat, pc);

      // compute res = -(Kmat*tangent - load*dlambda_ds)
      mat->mult(tangent, res);
      res->axpy(-dlambda_ds, load);
      res->scale(-1.0);

      // Set new values back into the matrix
      path_mat->resetConstraint(dlambda_ds);
      ksm->solve(res, temp);
      dlambda_ds += path_mat->extract(temp);

      // tangent = tangent + temp
      tangent->axpy(1.0, temp);
    }

    // Save values and update the iteration counter
    lambda_history[iteration_count] = lambda;

    // Set the derivative of the load variable w.r.t. path parameter
    dlambda_ds_history[iteration_count] = dlambda_ds;

    // Check for the termination conditions
    if (term_function) {
      TacsScalar value;
      assembler->evalFunctions(1, &term_function, &value);
      if (TacsRealPart(value) < TacsRealPart(term_function_value)) {
        iteration_count++;  // make the counter match the number of iterations
        if (ksm_print) {
          ksm_print->print(
              "TACSContinuation::solve_tangent: "
              "Terminating due to function value\n");
        }
        break;
      }
    }
    if (TacsRealPart(dlambda_ds) < TacsRealPart(dlambda_ds_term_value)) {
      iteration_count++;
      if (ksm_print) {
        ksm_print->print(
            "TACSContinuation::solve_tangent: "
            "Terminating for collapse condition\n");
      }
      break;
    }

    if (iteration_count == 0) {
      delta_s = target_delta_lambda / dlambda_ds;

      TacsScalar tnorm = tangent->norm();

      // Assign the target change in the r = (u - u_k)^{T}(u - u_k)
      // dr = sqrt( theta - dlambda_ds^2 ) = ||t||_{2}/sqrt( 1 + ||t||_{2}^2 )
      target_delta_r =
          delta_s * tnorm / sqrt(dlambda_ds * dlambda_ds + tnorm * tnorm);
    } else {
      // Find ds based on target_delta_r
      TacsScalar tnorm = tangent->norm();

      if (tnorm != 0.0 && dlambda_ds != 0.0) {
        delta_s = target_delta_r *
                  sqrt(dlambda_ds * dlambda_ds + tnorm * tnorm) / (tnorm);
      } else {
        ksm_print->print("Encountered error with step size selection\n");
      }
    }

    if (ksm_print) {
      char line[256];
      sprintf(line, "Outer iteration %3d: t: %9.4f dp_ds: %10.4e\n",
              iteration_count, MPI_Wtime() - t0, TacsRealPart(dlambda_ds));
      ksm_print->print(line);
      sprintf(line, "%5s %9s %10s %10s %10s\n", "Iter", "t", "|R|", "lambda",
              "|u|");
      ksm_print->print(line);
    }

    // Store the values of u, lambda
    old_vars->copyValues(vars);
    lambda_old = lambda;

    // Set the tolerances for the correction update iterations
    ksm->setTolerances(krylov_rtol, krylov_atol);

    // Try using the current solution again
    int fail_flag = 1;
    int nrestarts = 0;
    for (; fail_flag && (nrestarts < max_correction_restarts); nrestarts++) {
      // Perform an update based on the calculated value of ds
      // This ensures that the step lenght constraint is satisfied
      vars->axpy(delta_s, tangent);
      lambda = lambda + dlambda_ds * delta_s;

      // Set the ksm to use the path_mat object
      ksm->setOperators(path_mat, pc);

      // Now compute the next iteration
      TacsScalar init_res_norm = 0.0;
      for (int j = 0; j < max_correction_iters; j++) {
        // Compute the residual at the current value of (u, lambda)
        assembler->setVariables(vars);
        assembler->assembleRes(res);
        res->axpy(-lambda, load);

        TacsScalar res_norm = res->norm();
        if (ksm_print) {
          char line[256];
          sprintf(line, "%5d %9.4f %10.3e %10.3e %10.3e\n", j, MPI_Wtime() - t0,
                  TacsRealPart(res_norm), TacsRealPart(lambda),
                  TacsRealPart(vars->norm()));
          ksm_print->print(line);
        }

        // Set the initial norm or check the rtol/dtol
        if (j == 0) {
          init_res_norm = res_norm;
        } else if (TacsRealPart(res_norm) <
                   correction_rtol * TacsRealPart(init_res_norm)) {
          fail_flag = 0;
          break;
        } else if (TacsRealPart(res_norm) >
                   correction_dtol * TacsRealPart(init_res_norm)) {
          break;
        }

        // Set new values back into the matrix
        path_mat->resetConstraint(dlambda_ds);
        ksm->solve(res, temp);

        TacsScalar delta_lambda = path_mat->extract(temp);
        lambda = lambda - delta_lambda;
        vars->axpy(-1.0, temp);
      }

      // The corrector has failed. Try again with a smaller step size.
      if (fail_flag) {
        vars->copyValues(old_vars);
        lambda = lambda_old;
        delta_s = 0.5 * delta_s;

        if (ksm_print) {
          char line[256];
          sprintf(line,
                  "Failed to converge, retrying with step size = %10.3e\n",
                  TacsRealPart(delta_s));
          ksm_print->print(line);
        }
      }
    }

    if (nrestarts >= max_correction_restarts) {
      break;
    }

    if (callback) {
      callback->iteration(iteration_count, vars, lambda, dlambda_ds, assembler);
    }
  }

  // Deallocate temporary variables
  vars->decref();
  old_vars->decref();
  temp->decref();
  tangent->decref();
  update->decref();
  res->decref();
  path_mat->decref();
}
