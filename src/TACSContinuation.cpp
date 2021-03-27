#include "TACSContinuation.h"

/*!
  Implementation of buckling and frequency analysis and sensitivity
  analysis of eigenvalues.

  Copyright (c) 2015 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/


/*!
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
PathMat::PathMat( TACSMat * _A, 
                  TACSVec * _r, TACSVec * _t, TacsScalar s ){
  A = _A; A->incref();
  r = _r; r->incref();
  t = _t; t->incref();   
  xtmp = A->createVec(); xtmp->incref();
  resetConstraint(s);
}

PathMat::~PathMat(){
  A->decref();
  r->decref();
  t->decref();
  xtmp->decref();
}

void PathMat::getVectors( TACSVec ** _r, TACSVec ** _t ){
  if (_r){ *_r = r; }
  if (_t){ *_t = t; }
}

void PathMat::resetConstraint( TacsScalar s ){
  TacsScalar tnorm = t->norm();
  TacsScalar t2 = tnorm*tnorm;
  
  wn = s + (s >= 0.0 ? 1.0 : -1.0)*sqrt(t2 + s*s);
  tn = sqrt(t2 + wn*wn);
}

// Multiply x <-- Qx, return the value of the n+1-th row
TacsScalar PathMat::extract( TACSVec * x ){
  TacsScalar tTx = t->dot(x);
  x->axpy(-(2.0*tTx)/(tn*tn), t);
  return -(2.0*wn*tTx)/(tn*tn);
}

void PathMat::getSize( int * _nr, int * _nc ){ A->getSize( _nr, _nc ); }

void PathMat::mult( TACSVec * x, TACSVec * y ){
  xtmp->copyValues(x);
  
  TacsScalar tTx = t->dot(x);
  xtmp->axpy(-(2.0*tTx)/(tn*tn), t);
  
  A->mult(xtmp, y);
  y->axpy(-(2.0*wn*tTx)/(tn*tn), r);
}

TACSVec * PathMat::createVec(){ return A->createVec(); }

/*!
  Callback class for monitoring the continuation algorithm
*/
TACSContCallback::TACSContCallback( MPI_Comm _comm,
                                    const char * _file_prefix,
                                    const char * _monitor_file,
                                    int * _var_nums, int _nvars ){
  comm = _comm;

  if ( _file_prefix ){
    file_prefix = new char[strlen(_file_prefix)+1];
    strcpy(file_prefix, _file_prefix);
  }
  else {
    file_prefix = NULL;
  }

  if ( _monitor_file ){
    monitor_file = new char[strlen(_monitor_file)+1];
    strcpy(monitor_file, _monitor_file);
  }
  else {
    monitor_file = NULL;
  }
  mfp = NULL;

  nvars = _nvars;
  var_nums = new int[ nvars ];
  for ( int k = 0; k < nvars; k++ ){
    var_nums[k] = _var_nums[k];
  }

  var_values = new TacsScalar[ nvars ];
  memset(var_values, 0, sizeof(TacsScalar)*nvars);

  unit_bc_displacement = NULL;
  res_no_bcs = NULL;
}

TACSContCallback::~TACSContCallback(){
  if ( file_prefix ){ delete [] file_prefix; }
  if ( monitor_file ){ delete [] monitor_file; }
  delete [] var_nums;
  delete [] var_values;
  if ( mfp ){ fclose(mfp); }
  if (unit_bc_displacement){ unit_bc_displacement->decref(); }
  if (res_no_bcs){ res_no_bcs->decref(); }
}

void TACSContCallback::setUnitBCDisplacement( BVec * unit ){
  unit->incref();
  if (unit_bc_displacement){
    unit_bc_displacement->decref();
  }
  unit_bc_displacement = unit;
  if (res_no_bcs){
    res_no_bcs->decref();
    res_no_bcs = NULL;
  }
}

void TACSContCallback::initialize(){
  int rank;
  MPI_Comm_rank(comm, &rank);

  if (rank == 0){
    if (mfp){
      fclose(mfp);
    }

    if (monitor_file){
      mfp = fopen(monitor_file, "w");
    }

    if (mfp){    
      fprintf(mfp, "Variables = iter, lambda, dlambda_ds, Fbc");
      for ( int k = 0; k < nvars; k++ ){
        fprintf(mfp, ", d%d", k+1);
      }
      fprintf(mfp, "\n");
      fprintf(mfp, "Zone T=\"Continuation points\"\n");      
    }
  }
  else {
    mfp = NULL;
  }
}

void TACSContCallback::iteration( int iter, TacsScalar lambda,
                                  TacsScalar dlambda_ds,
                                  int loadCase, TACSAssembler * tacs,
                                  BVec * vars ){

  // Get the values on the root process
  TacsScalar * temp = new TacsScalar[nvars];
  memset(var_values, 0, nvars*sizeof(TacsScalar));
  memset(temp, 0, nvars*sizeof(TacsScalar));
 
  VarMap * rmap = tacs->getVarMap();
  int rank, size;
  const int * range;
  int bs = rmap->getBlockSize();

  rmap->getOwnerRange(&range, &rank, &size);
  TacsScalar * xvars;
  vars->getArray(&xvars);

  for ( int k = 0; k < nvars; k++ ){
    if (var_nums[k] >= bs*range[rank] && 
        var_nums[k] < bs*range[rank+1]){
      temp[k] = xvars[var_nums[k] - bs*range[rank]];
    }
  }

  MPI_Reduce(temp, var_values, nvars, TACS_MPI_TYPE, 
             MPI_SUM, 0, comm);
  
  TacsScalar Fbc = 0.0;
  if (unit_bc_displacement){
    if (!res_no_bcs){
      res_no_bcs = tacs->createVec();
      res_no_bcs->incref();
    }
    tacs->assembleResNoBCs(loadCase, res_no_bcs);
    Fbc = res_no_bcs->dot(unit_bc_displacement);
  }

  if (mfp){
    fprintf(mfp, "%2d %15.6e %15.6e %15.6e ", 
            iter+1, RealPart(lambda), RealPart(dlambda_ds),
	    RealPart(Fbc));
    for ( int k = 0; k < nvars; k++ ){
      fprintf(mfp, "%15.6e ", RealPart(var_values[k]));
    }
    fprintf(mfp, "\n");
    fflush(mfp);    
  } 

  if (file_prefix){
    char * filename = new char[strlen(file_prefix)+15];
    sprintf(filename, "%s%03d.bin", file_prefix, iter);
    vars->writeToFile(filename);
    delete [] filename;
  }
}

void TACSContCallback::finalize(){
  if (mfp){ fclose(mfp); }
  mfp = NULL;
}

/*!
  Initialize the arc-length continuation class
*/

TACSContinuation::TACSContinuation( TACSAssembler * _tacs, 
                                    int _max_continuation_iters ){
  tacs = _tacs;
  tacs->incref();

  max_continuation_iters  = _max_continuation_iters;
  max_correction_iters    = 25;
  max_correction_restarts = 4;

  tangent_rtol = 1e-8;
  tangent_atol = 1e-30;

  krylov_rtol = 1e-3;
  krylov_atol = 1e-30;

  correction_rtol = 1e-8;
  correction_dtol = 1e3;

  iteration_count = 0;
  lambda_history = new TacsScalar[ max_continuation_iters ];
  dlambda_ds_history = new TacsScalar[ max_continuation_iters ];

  memset(lambda_history, 0, max_continuation_iters*sizeof(TacsScalar));
  memset(dlambda_ds_history, 0, max_continuation_iters*sizeof(TacsScalar));

  // Set the default termination information
  term_function = NULL;
  term_function_value = 1.0;
  dlambda_ds_term_value = -1e20;

  callback = NULL;
}

TACSContinuation::~TACSContinuation(){
  tacs->decref();

  delete [] lambda_history;
  delete [] dlambda_ds_history;

  if ( term_function ){ term_function->decref(); }
  if ( callback ){ callback->decref(); }
}

/*!
  Set the callback object to be used
*/

void TACSContinuation::setCallback( TACSContCallback * _callback ){
  _callback->incref();
  if (callback){ callback->decref(); }
  callback = _callback;
}

/*!
  Set the iteration limits/tolerances
*/

void TACSContinuation::setCorrectionIters( int _max_correction_iters,
                                           int _max_correction_restarts ){
  max_correction_iters    = _max_correction_iters;  
  max_correction_restarts = _max_correction_restarts;  
}

void TACSContinuation::setCorrectionTol( double rtol, double dtol ){  
  if (rtol < 0.1){
    correction_rtol = rtol;
  }
  if (dtol > 10.0){
    correction_dtol = dtol;
  }
}

/*!
  Set the termination conditions
*/

void TACSContinuation::setTermFunction( TACSFunction * func, 
                                        TacsScalar term_value ){
  func->incref();
  if ( term_function ){ term_function->decref(); }
  term_function = func;
  term_function_value = term_value;
}

void TACSContinuation::setTermLambdaRate( TacsScalar term_dlambda_ds ){
  dlambda_ds_term_value = term_dlambda_ds;
}

/*!
  Retrieve information about the solve
*/
int TACSContinuation::getNumIterations(){
  return iteration_count;
}

void TACSContinuation::getSolution( int iter, TacsScalar * lambda, 
                                    TacsScalar * dlambda_ds ){
  if ( iter >= 0 && iter < iteration_count ){
    *lambda = lambda_history[iter];
    *dlambda_ds = dlambda_ds_history[iter];
  }
  else {
    *lambda = 0.0;
    *dlambda_ds = 0.0;
  }
}

/*!
  Perform arc-length continuation to track the equilibrium path past
  limit points (folds).

  1. Compute the tangent vector at the current point and compute the
  first order update.

  The solution curve is given by (u,lambda), and represents solutions
  of the following system of equations:

  R' * du + R,lambda * lambda = 0 (1)

  The tangent to the curve at a point is given by t = [R']^{-1} R,lambda,
  where R,lambda == the load vector.

  An additional equation is added to (1) to parameterize the curve
  through the arc-length. This is performed by employing a
  linearization of the arc-length constraint:
  
  t^{T}(u_{k+1} - u_{k}) + (lambda_{k+1} - lambda_k) = delta_s/(dlambda_ds)

  2. Perform corrector iterations to bring the solution back to the
  equilibrium path. These corrector iterations employ a Newton method
  to bring the point back onto the equilibrium path. 
*/
void TACSContinuation::solve_tangent( int loadCase,
                                      TACSMat * mat, TACSPc * pc, 
                                      TACSKsm * ksm,
                                      TacsScalar lambda_init,
                                      TacsScalar target_delta_lambda,
                                      KSMPrint * ksm_print ){

  BVec * vars = tacs->createVec();   
  BVec * old_vars = tacs->createVec();      
  BVec * temp = tacs->createVec();      
  BVec * tangent = tacs->createVec();
  BVec * update = tacs->createVec();
  BVec * res = tacs->createVec();
  BVec * load = tacs->createVec();

  vars->incref();
  old_vars->incref();
  temp->incref();
  tangent->incref();
  update->incref();
  res->incref();
  load->incref();

  if ( callback ){ callback->initialize(); }
  
  PathMat * path_mat = new PathMat(mat, load, tangent, 0.0);
  path_mat->incref();                               

  TacsScalar lambda = 0.0;         // The load factor
  TacsScalar lambda_old = 0.0;     // The previous load factor
  TacsScalar target_delta_r = 0.0; // Target change in r = (u - u_k)^T(u - u_k)

  double t0 = MPI_Wtime();

  if (lambda_init != 0.0){
    lambda = lambda_init;
 
    if (ksm_print){
      char line[256];
      ksm_print->print("Performing initial Newton iterations\n");
      sprintf(line, "%5s %9s %10s\n", "Iter", "t", "|R|");
      ksm_print->print(line);
    }

    // Set the tolerances for the tangent computation
    ksm->setTolerances(tangent_rtol, tangent_atol);

    // Compute and factor the stiffness matrix
    tacs->setLoadFactor(loadCase, 0.0);
    tacs->assembleMatType(loadCase, mat, 1.0, STIFFNESS_MATRIX, NORMAL);
    pc->factor();

    // Compute the tangent vector
    tacs->assembleResLoadFactor(loadCase, load);
    ksm->solve(load, tangent);

    // Update the variables based on the tangent computation
    vars->axpy(-lambda, tangent);
    tacs->setLoadFactor(loadCase, lambda);

    // Compute the initial norm based on the tangent approximation
    TacsScalar res_norm_init = 1.0;
    
    // Now perform a Newton iteration until convergence
    for ( int k = 0; k < max_continuation_iters; k++ ){
      tacs->setVariables(loadCase, vars);
      tacs->assembleRes(loadCase, res);

      TacsScalar res_norm = res->norm();
      if (ksm_print){
        char line[256];
        sprintf(line, "%5d %9.4f %10.4e\n",
		k+1, MPI_Wtime() - t0, RealPart(res_norm));
        ksm_print->print(line);
      }
      if (k == 0){
	res_norm_init = res_norm;
      }
      if (res_norm < correction_rtol*res_norm_init){
        break;
      }

      ksm->solve(res, update);
      vars->axpy(-1.0, update);
    }
  }

  if (ksm_print){
    ksm_print->print("Beginning arc-length continuation method\n");
  }

  // The rate of change of the lambda w.r.t. the arc-length
  TacsScalar dlambda_ds = 0.0;

  for ( iteration_count = 0; iteration_count < max_continuation_iters; 
        iteration_count++ ){

    // Copy the current values to a vector
    tacs->setLoadFactor(loadCase, lambda);
    tacs->setVariables(loadCase, vars);
    
    // Assemble the stiffness matrix at the current iteration
    tacs->assembleMatType(loadCase, mat, 1.0, STIFFNESS_MATRIX, NORMAL);
    pc->factor(); 
    tacs->assembleResLoadFactor(loadCase, load);

    // Compute the change in r = (u - u_k)^{T}(u - u_k)
    TacsScalar delta_s = 1.0; // this will be over-written later

    // Set the tolerances for the tangent computation
    ksm->setTolerances(tangent_rtol, tangent_atol);

    if (iteration_count == 0){
      // Compute the initial tangent vector to the solution path
      ksm->setOperators(mat, pc);
      ksm->solve(load, tangent);

      TacsScalar tnorm = tangent->norm();
      dlambda_ds = 1.0/sqrt(1.0 + tnorm*tnorm);    

      tangent->scale(-dlambda_ds);
    }
    else {
      // Set the ksm to use the path_mat object
      ksm->setOperators(path_mat, pc);

      // compute res = -(Kmat*tangent + load*dlambda_ds)
      mat->mult(tangent, res);
      res->axpy(dlambda_ds, load);
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
    if (term_function){
      if (tacs->evalFunction(loadCase, term_function) < term_function_value){
        iteration_count++; // make the counter match the number of iterations
        if (ksm_print){
          ksm_print->print("Terminating due to function value\n");
        }
        break;
      }
    }
    if (dlambda_ds < dlambda_ds_term_value){
      iteration_count++;
      if (ksm_print){
        ksm_print->print("Terminating for collapse condition\n");
      }
      break;
    }

    if (iteration_count == 0){
      delta_s = target_delta_lambda/dlambda_ds;
      
      TacsScalar tnorm = tangent->norm();

      // Assign the target change in the r = (u - u_k)^{T}(u - u_k) 
      // dr = sqrt( theta - dlambda_ds^2 ) = ||t||_{2}/sqrt( 1 + ||t||_{2}^2 )
      target_delta_r = delta_s*tnorm/sqrt(dlambda_ds*dlambda_ds + 
                                          tnorm*tnorm);
    }
    else {
      // Find ds based on target_delta_r
      TacsScalar tnorm = tangent->norm();
        
      if (tnorm != 0.0 && dlambda_ds != 0.0){
        delta_s = target_delta_r*sqrt(dlambda_ds*dlambda_ds +
                                      tnorm*tnorm)/(tnorm);
      }
      else {
        ksm_print->print("Encountered error with step size selection\n");
      }
    }

    if (ksm_print){
      char line[256];
      sprintf(line, "Outer iteration %3d: t: %9.4f dp_ds: %10.4e\n",
              iteration_count, MPI_Wtime() - t0, RealPart(dlambda_ds));
      ksm_print->print(line);
      sprintf(line, "%5s %9s %10s %10s %10s\n",
              "Iter", "t", "|R|", "lambda", "|u|");
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
    for ( ; fail_flag && (nrestarts < max_correction_restarts); 
          nrestarts++ ){
      // Perform an update based on the calculated value of ds
      // This ensures that the step lenght constraint is satisfied
      vars->axpy(delta_s, tangent);
      lambda = lambda + dlambda_ds*delta_s;
    
      // Set the ksm to use the path_mat object
      ksm->setOperators(path_mat, pc);

      // Now compute the next iteration
      TacsScalar init_res_norm = 0.0;
      for ( int j = 0; j < max_correction_iters; j++ ){
        // Compute the residual at the current value of (u, lambda)
        tacs->setLoadFactor(loadCase, lambda);
        tacs->setVariables(loadCase, vars);
	tacs->assembleRes(loadCase, res);
	
	TacsScalar res_norm = res->norm();
        if (ksm_print){
          char line[256];
          sprintf(line, "%5d %9.4f %10.3e %10.3e %10.3e\n",
                  j, MPI_Wtime() - t0, RealPart(res_norm), RealPart(lambda), 
                  RealPart(vars->norm()));
          ksm_print->print(line);
        }
        
        // Set the initial norm or check the rtol/dtol 
        if (j == 0){
          init_res_norm = res_norm;
        }
        else if (res_norm < correction_rtol*init_res_norm){
          fail_flag = 0;
          break;
        }
        else if (res_norm > correction_dtol*init_res_norm){
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
      if (fail_flag){
        vars->copyValues(old_vars);
        lambda = lambda_old;
        delta_s = 0.5*delta_s;

        if (ksm_print){
          char line[256];
          sprintf(line, 
                  "Failed to converge, retrying with step size = %10.3e\n", 
                  RealPart(delta_s));
          ksm_print->print(line);
        }
      }
    }

    if (nrestarts >= max_correction_restarts){
      break;
    }

    if (callback){ 
      callback->iteration(iteration_count, lambda, dlambda_ds, 
                          loadCase, tacs, vars);
    }
  }

  if (callback){ callback->finalize(); }

  // Clean up the temporary vectors
  vars->decref();
  old_vars->decref();
  temp->decref();
  tangent->decref();
  update->decref();
  res->decref();
  load->decref();

  path_mat->decref();
} 

/*!
  Perform arc-length continuation to track the equilibrium path past
  limit points (folds).

  1. Compute the tangent vector at the current point and compute the
  first order update.

  The solution curve is given by (u,lambda), and represents solutions
  of the following system of equations:

  R' * du + f * lambda = 0 (1)

  The tangent to the curve at a point is given by t = [R']^{-1}R,lambda.

  An additional equation is added to (1) to parameterize the curve
  through the arc-length. This is performed by enforcing the following
  additional cylindrical constraint, outlined by Crisfield,

  (u_{k+1} - u_{k})^2 = (s_{k+1} - s_{k})^2 

  2. Perform corrector iterations to bring the solution back to the
  equilibrium path. These corrector iterations employ a Newton method
  to bring the point back onto the equilibrium path. The cylindrical
  arc-length constraint is imposed simultaneously.  
*/
/*!
void TACSContinuation::solve_cylindrical( int loadCase,
                                          TACSMat * mat, TACSPc * pc, 
                                          TACSKsm * ksm,
                                          TacsScalar target_delta_lambda ){
*/

