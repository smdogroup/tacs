#include "TACSIntegrator.h"
#include <math.h>
#include "tacslapack.h"

/* 
   Static factory method that returns an instance of the concrete
   child class. Note that the base class is abstract (contains
   unimplemented pure virtual functions), thus it can not be
   instantiated.
   
   input:
   tacs              : the TACS assembler object
   tinit             : start time of simulation
   tfinal            : end time of simulation
   num_steps_per_sec : the number of steps to take per second
   type              : type of integrator (BDF1-3, DIRK2-4, ABM1-6, NBG)
*/
TACSIntegrator* TACSIntegrator::getInstance( TACSAssembler * _tacs,
                                             double _tinit, double _tfinal, 
                                             int _num_steps_per_sec, 
                                             enum IntegratorType type ){
  if ( type == DIRK2 ) {
    fprintf(stdout, ">> TACSIntegrator: Instantiating DIRK Order 2 integrator...\n");
    return new TACSDIRKIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 1);
  } else if ( type == DIRK3 ) {
    fprintf(stdout, ">> TACSIntegrator: Instantiating DIRK Order 3 integrator...\n");
    return new TACSDIRKIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 2);
  } else if ( type == DIRK4 ) {
    fprintf(stdout, ">> TACSIntegrator: Instantiating DIRK Order 4 integrator...\n");
    return new TACSDIRKIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 3);

  } else if ( type == BDF1 ) {
    fprintf(stdout, ">> TACSIntegrator: Instantiating BDF Order 1 integrator...\n");
    return new TACSBDFIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 1);
  } else if ( type == BDF2 ) {
    fprintf(stdout, ">> TACSIntegrator: Instantiating BDF Order 2 integrator...\n");
    return new TACSBDFIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 2);
  } else if ( type == BDF3 ) {
    fprintf(stdout, ">> TACSIntegrator: Instantiating BDF Order 3 integrator...\n");
    return new TACSBDFIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 3);

  } else if ( type == ABM1 ) {
    fprintf(stdout, ">> TACSIntegrator: Instantiating ABM Order 1 integrator...\n");
    return new TACSABMIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 1);
  } else if ( type == ABM2 ) {
    fprintf(stdout, ">> TACSIntegrator: Instantiating ABM Order 2 integrator...\n");
    return new TACSABMIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 2);
  } else if ( type == ABM3 ) {
    fprintf(stdout, ">> TACSIntegrator: Instantiating ABM Order 3 integrator...\n");
    return new TACSABMIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 3);
  } else if ( type == ABM4 ) {
    fprintf(stdout, ">> TACSIntegrator: Instantiating ABM Order 4 integrator...\n");
    return new TACSABMIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 4);
  } else if ( type == ABM5 ) {
    fprintf(stdout, ">> TACSIntegrator: Instantiating ABM Order 5 integrator...\n");
    return new TACSABMIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 5);
  } else if ( type == ABM6 ) {
    fprintf(stdout, ">> TACSIntegrator: Instantiating ABM Order 6 integrator...\n");
    return new TACSABMIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 6);

    // Newmark Beta Gamma Method
    } else if ( type == NBG ) {
    fprintf(stdout, ">> TACSIntegrator: Instantiating NBG integrator...\n");
    return new TACSNBGIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec);
  } else { 
    // Default NBG integrator
    fprintf(stdout, ">> TACSIntegrator: Instantiating NBG integrator...\n");
    return new TACSNBGIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec);
  }  
}

/*
  Base class constructor for integration schemes.

  input:
  tacs              : the tacs assembler object
  tInit             : the initial time
  tfinal            : the final time
  num_steps_per_sec : the number of steps to take for each second
*/
TACSIntegrator::TACSIntegrator( TACSAssembler * _tacs,
                                double _tinit, double _tfinal, 
                                int _num_steps_per_sec ){
  // Copy over the input parameters
  tacs              = _tacs; tacs->incref();
  tinit             = _tinit;
  tfinal            = _tfinal;
  num_steps_per_sec = _num_steps_per_sec;

  // MPI information
  MPI_Comm_rank(tacs->getMPIComm(), &mpiRank);
  MPI_Comm_size(tacs->getMPIComm(), &mpiSize);
  
  // Compute the step size
  h = 1.0/num_steps_per_sec;

  // compute the total number of time steps
  num_time_steps = int(num_steps_per_sec*(tfinal-tinit)) + 1;

  // Default print level and logging control
  print_level = 1;
  logfp       = stdout;
  logfilename = NULL;
  if (mpiRank != 0){logfp = NULL;}   // Only root writes log

  // State variables that store the entire time history
  q     = new TACSBVec*[num_time_steps];
  qdot  = new TACSBVec*[num_time_steps];
  qddot = new TACSBVec*[num_time_steps];

  // create state vectors for TACS during each timestep
  for ( int k = 0; k < num_time_steps; k++ ) {
    q[k]     = tacs->createVec(); q[k]->incref(); 
    qdot[k]  = tacs->createVec(); qdot[k]->incref(); 
    qddot[k] = tacs->createVec(); qddot[k]->incref(); 
  }

  // Get the number of state variables and store into the class variable
  q[0]->getSize(&num_state_vars);
  
  // Store time
  time = new double[num_time_steps];
  memset(time, 0, num_time_steps*sizeof(double));

  //------------------------------------------------------------------//
  //                     Newton's method                              //
  //------------------------------------------------------------------//

  // Frequency of Jacobian recomputation during nonlinear solve
  jac_comp_freq = 1;

  // Set the default LINEAR solver
  use_lapack = 0;
    
  // Default parameters for Newton solve
  max_newton_iters = 25;
  atol = 1.0e-12;
  rtol = 1.0e-10;

  // Create vector for storing the residual at each Newton iteration
  res = tacs->createVec();
  res->incref();
  
  // Create a vector for storing the Newton updates
  update = tacs->createVec();
  update->incref();
  
  // Create a matrix for storing the Jacobian
  D = tacs->createFEMat();
  D->incref();

  // Associate the maxtrix with FEMatrix
  mat = D;
  mat->incref();
    
  // Allocate the factorization
  int lev = 4500; double fill = 10.0; int reorder_schur = 1;
  pc = new PcScMat(D, lev, fill, reorder_schur);
  pc->incref();
    
  // The Krylov subspace method (KSM) associated with the solver
  int gmres_iters = 10, num_restarts = 0, is_flexible = 0;
  ksm = new GMRES(D, pc, gmres_iters, num_restarts, is_flexible);
  ksm->incref();
  
  ksm->setTolerances(rtol, atol);

  // Variables used in adjoint solve (use setFunction(...) to set these 
  num_funcs = 0;
  funcs     = NULL;

  // Tecplot solution export (use configureOutput(...) to set these
  f5_write_freq = 0;
  f5_file_fmt   = NULL;
  f5            = NULL;

  // Allocate space for kinetic and potential energies
  energies = new TacsScalar[ 2 ];
  energies[0] = 0.0;
  energies[1] = 0.0;

  init_energy = 0.0;
}

/*
  Destructor for base class resources
*/
TACSIntegrator::~TACSIntegrator(){
  // Close any open file pointers
  if (logfp != stdout && logfp){
    fclose(logfp);
  }
  
  // Dereference position, velocity and acceleration states
  for ( int k = 0; k < num_time_steps; k++ ) {
    q[k]->decref();
    qdot[k]->decref();
    qddot[k]->decref();
  }
  
  // Dereference Newton's method objects
  res->decref();
  update->decref();
  D->decref();
  mat->decref();
  pc->decref();
  ksm->decref();
  
  if (time)     { delete [] time; }
  if (q)        { delete [] q; }
  if (qdot)     { delete [] qdot; }
  if (qddot)    { delete [] qddot; }

  if (energies) { delete [] energies; }
  
  // Dereference TACS
  tacs->decref();
}

/*
  Drives the residual R(t,q,qdot,qddot) to zero using Newton's method

  Input: 
  The guessed (initial) state variable values q, qdot, qddot are
  supplied
  
  Output: q, qdot, qddot updated iteratively until the corresponding
  residual R <= tolerance
  
  alpha: multiplier for derivative of Residual wrt to q
  beta : multiplier for derivative of Residual wrt to qdot
  gamma: multiplier for derivative of Residual wrt to qddot
*/
void TACSIntegrator::newtonSolve( double alpha, double beta, double gamma,
                                  double t, TACSBVec *u, TACSBVec *udot, 
                                  TACSBVec *uddot ){
  // Initialize the norms
  init_norm = 0.0;
  norm = 0.0;

  if ( logfp && print_level >= 2){
    fprintf(logfp, "%12s %8s %12s %12s %12s %12s\n",
            "time", "NIter", "tcpu", "|R|", "|R|/|R0|", "delta");
  }
  
  double t0 = MPI_Wtime();

  // Iterate until max iters or R <= tol
  double delta = 0.0;
  for ( niter = 0; niter < max_newton_iters; niter++ ){
    // Set the supplied initial input states into TACS
    setTACSStates(t, u, udot, uddot);

    // Use a globalization based on a dual-type time step method
    /*
      double delta = alpha;
      if (n > 0){
      double frac = 10.0*RealPart(norm/(init_norm + rtol));
      if (frac < 1.0){
      delta = frac*alpha;
      }
      else {
      delta = alpha;
      }
      }
    */
    // Assemble the Jacobian matrix once in five newton iterations
    if (niter % jac_comp_freq == 0){
      tacs->assembleJacobian(alpha, beta + delta, gamma + delta,
                             res, mat, NORMAL);
    }
    else {
      tacs->assembleRes(res);
    }
   
    // Compute the L2-norm of the residual
    norm = res->norm();
    
    // Record the residual norm at the first Newton iteration
    if (niter == 0){
      init_norm = norm;
    }

    // Write a summary    
    if( logfp && print_level >= 2) {
      fprintf(logfp, "%12.5e %8d %12.5e %12.5e %12.5e %12.5e\n",
              t, niter, MPI_Wtime()-t0, 
              RealPart(norm),  (niter == 0) ? 1.0 : RealPart(norm/init_norm), delta);
    }

    // Check if the norm of the residuals is a NaN
    if (norm != norm){ 
      fprintf(stderr,
              "Newton iteration %d, failed with NaN residual norm\n", niter);
      break;
    }
    
    // Check if the Newton convergence tolerance is satisfied
    if (norm < rtol*init_norm || norm < atol){
      break;
    }
    
    if (!use_lapack) {    
      // LU Factor the matrix when needed
      if (niter % jac_comp_freq == 0){
	pc->factor();
      }  
      // Solve for update using KSM
      ksm->solve(res, update);
    }
    else {
      // Perform the linear solve using LAPACK (serial only)
      lapackLinearSolve(res, mat, update);
    }

    // Update the state variables using the solution
    uddot->axpy(-gamma, update);
    udot->axpy(-beta, update);
    u->axpy(-alpha, update);

    // Check whether the Newton iteration was successful
    if (niter == max_newton_iters && norm >= rtol*init_norm){
      fprintf(stderr,"Newton iteration failed to converge in %d iters\n", niter);
      break;
    }
  }
}

/*
  Function that writes time, q, qdot, qddot to file
*/
void TACSIntegrator::writeSolution( const char *filename ) {
  if (mpiRank != 0) { return; }
  // Temporary variables to access the states at each time
  TacsScalar *qvals, *qdotvals, *qddotvals;

  // Open a new file
  FILE *fp = fopen(filename, "w");
 
  for ( int k = 0; k < num_time_steps; k++ ){    
    // Copy over the state values from TACSBVec
    q[k]->getArray(&qvals);
    qdot[k]->getArray(&qdotvals);
    qddot[k]->getArray(&qddotvals);
  
    // Write the time and states to file
    fprintf(fp, "%e ", time[k]);
    for ( int j = 0; j < num_state_vars; j++ ){
      fprintf(fp, "%e %e %e ", RealPart(qvals[j]), 
              RealPart(qdotvals[j]), RealPart(qddotvals[j]));
    }
    fprintf(fp, "\n");
  }

  // Close the file
  fclose(fp);
}

/*
  Don't write if f5_write_freq is set to 0. If write output is sought,
  force writing initial(k=0) and final (k=num_time_steps-1) time steps
  as they are the most important outputs that an user would need
*/
int TACSIntegrator::getWriteFlag( int k ){
  int write_now;
  if (f5_write_freq > 0) {
    if ( k == 0 || k == num_time_steps-1 ){
      write_now = 1;
    } else {
      write_now = (k % f5_write_freq == 0);
    }
  }
  return write_now;
}

/*
  Writes the output as an f5 file. Note the states and time must be
  set appropriately before calling this function.
*/
void TACSIntegrator::writeStepToF5( int k ){
  if (mpiRank != 0) { return; }
 
  if(f5 && f5_file_fmt && getWriteFlag(k)) {
    // Create a buffer for filename 
    char buffer[128];
    // Format the buffer based on the time step
    getString(buffer, f5_file_fmt, k);      
    // Write the f5 file for this time step
    f5->writeToFile(buffer);
  }
}

/*
  Creates an f5 file for each time step and writes the data. This call
  is distributed in time. This may be faster than calling f5 output at each step
*/
void TACSIntegrator::writeSolutionToF5(){

  // Take write_flag and element type as inputs too.
  
  int is, ie, idec;
  idec = num_time_steps/mpiSize;
  is   = idec*mpiRank + 1;
  ie   = idec*(mpiRank + 1);
  if ( mpiRank == mpiSize - 1) { ie = num_time_steps; }

  //Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);

  // Create a viewer
  TACSToFH5 * f5 = new TACSToFH5(tacs, SHELL, write_flag);
  f5->incref();

  for ( int k = is; k < ie; k++ ){
    // Set the current states into TACS
    setTACSStates(time[k], q[k], qdot[k], qddot[k]);

    // Make a filename
    char fname[128];
    sprintf(fname, "output/solution_%04d.f5", k);

    // Write the solution
    f5->writeToFile(fname);

  }
  // Delete the viewer
  f5->decref();
}

/*
  Configure the F5 output 
 */
void TACSIntegrator::configureOutput(TACSToFH5 *_viewer, 
                                     int _write_freq, 
                                     char *_f5_file_fmt ) {
  f5             = _viewer;
  f5_write_freq  = _write_freq;
  f5_file_fmt    = _f5_file_fmt;
}

/*
  Integate forward in time using the initial conditions retrieved from
  TACS
*/
void TACSIntegrator::integrate( ) {
  // Get the initial condition
  current_time_step = 0;
  tacs->getInitConditions(q[0], qdot[0]);
  tacs->setVariables(q[0], qdot[0]);
  
  // Perform logging, tecplot export, etc.
  doEachTimeStep(0);

  for (int k = 1; k < num_time_steps; k++){
    // Set the class variable
    current_time_step = k;

    // Advance time
    time[k] = time[k-1] + h;
    
    // Approximate states and their derivatives using ABM formula
    approxStates();

    // Determine the coefficients for linearizing the Residual
    double alpha, beta, gamma;   
    getLinearizationCoeffs(&alpha, &beta, &gamma);

    // Solve the nonlinear system of stage equations starting with the approximated states
    newtonSolve(alpha, beta, gamma, time[k], q[k], qdot[k], qddot[k]);

    // Perform logging, tecplot export, etc.
    doEachTimeStep(k);
  }
}

/*
  Compute the gradient for the given functions of interest with
  respect to the design variable.
*/
void TACSIntegrator::getFuncGrad( int _num_dv,
                                  TacsScalar *_x,
                                  TacsScalar *_fvals,
                                  TacsScalar *_dfdx ) {
  // Copy the inputs
  num_design_vars = _num_dv;
  tacs->setDesignVars(_x, num_design_vars);

  // Check whether the function has been set properly
  if (num_funcs == 0 || funcs == NULL) {
    fprintf(stderr, "TACS Warning: Function is not set\n");
    return;
  }
  
  // Point the class variable to the memory allocated by the optimizer
  fvals = _fvals; dfdx = _dfdx;
  memset(fvals, 0, num_funcs*sizeof(TacsScalar));
  memset(dfdx, 0, num_funcs*num_design_vars*sizeof(TacsScalar));

  // Integrate forward in time to solve for the states (q, qdot, qddot)
  integrate();
  
  // March backwards in time and solve for the adjoint variables
  marchBackwards();

  // Reduce the gradients on all processors
  MPI_Allreduce(MPI_IN_PLACE, dfdx, num_funcs*num_design_vars, TACS_MPI_TYPE, 
                MPI_SUM, tacs->getMPIComm());
}

/*
  Function that returns the finite-difference/complex-step gradient
  (used for testing purposes)
*/
void TACSIntegrator::getFDFuncGrad( int _num_dv, TacsScalar *_x,
                                        TacsScalar *_fvals, TacsScalar *_dfdx, 
                                        double dh ) {
  // Copy the inputs
  num_design_vars = _num_dv;
  tacs->setDesignVars(_x, num_design_vars);

  // Check whether the function has been set properly
  if (num_funcs == 0 || funcs == NULL) {
    fprintf(stderr, "TACS Warning: Function is not set, skipping adjoint solve. \n");
    return;
  }
  
  // Point the class variable to the memory allocated by the optimizer
  fvals = _fvals; dfdx = _dfdx;
  memset(fvals, 0, num_funcs*sizeof(TacsScalar));
  memset(dfdx, 0, num_funcs*num_design_vars*sizeof(TacsScalar));

  TacsScalar *ftmp = new TacsScalar[num_funcs];
  memset(ftmp, 0, num_funcs*sizeof(TacsScalar));

  // Perform the forward integration if we're using a
  // finite-difference approximation

#ifndef TACS_USE_COMPLEX
  
  // Integrate forward in time
  integrate();

  // Evaluate the functions
  evalTimeAvgFunctions(funcs, num_funcs, fvals); 

#endif

  // Find a finite-difference (or complex-step) approximation of the
  // total derivative
  for ( int k = 0; k < num_design_vars; k++ ){
    TacsScalar xtmp = _x[k];
    
#ifdef TACS_USE_COMPLEX
    
    // Perturb the DV
    _x[k] = xtmp + TacsScalar(0.0, dh);

    // Set the perturbed vector into TACS
    tacs->setDesignVars(_x, num_design_vars);
    
    // Integrate with perturbed x
    integrate();
   
    // Evaluate the functions
    evalTimeAvgFunctions(funcs, num_funcs, fvals); 

    // Evaluate the CS derivative
    for ( int j = 0; j < num_funcs; j++ ){
      dfdx[k+j*num_design_vars] = ImagPart(fvals[j])/dh;
    }

#else
    
    // Perturn the DV
    _x[k] = xtmp + dh;

    // Set the perturbed vector into TACS
    tacs->setDesignVars(_x, num_design_vars);
    
    // Integrate with perturbed x
    integrate();

    // Evaluate the functions at the current time 
    evalTimeAvgFunctions(funcs, num_funcs, ftmp); 

    // Evaluate the FD derivative
    for ( int j = 0; j < num_funcs; j++ ){
      dfdx[k+j*num_design_vars] = (ftmp[j] - fvals[j])/dh;
    }

#endif // TACS_USE_COMPLEX

    // Restrore the perturbed vector
    _x[k] = xtmp;
  }

  // Make sure the DV's in TACS are the same as before
  tacs->setDesignVars(_x, num_design_vars);

  // Reduce the gradients on all processors
  MPI_Allreduce(MPI_IN_PLACE, dfdx, num_funcs*num_design_vars, TACS_MPI_TYPE, 
                MPI_SUM, tacs->getMPIComm());

  delete [] ftmp;
}

/*
  Evaluate time average of the function value using discretization
  from the integration scheme (default implementation)
*/
void TACSIntegrator::evalTimeAvgFunctions( TACSFunction **funcs, 
                                           int numFuncs, 
                                           TacsScalar *funcVals) {
  memset(funcVals, 0, numFuncs*sizeof(TacsScalar));
 
  TacsScalar *ftmp = new TacsScalar[numFuncs];  
  memset(ftmp, 0, numFuncs*sizeof(TacsScalar));
  
  // Loop over time steps
  for ( int k = 1; k < num_time_steps; k++ ) {
    // Set the states into TACS
    setTACSStates(time[k], q[k], qdot[k], qddot[k]);
    
    // Evaluate the functions

    tacs->evalFunctions(funcs, numFuncs, ftmp); 
    
    // Compute the mean function value (h is the weight)
    for ( int j = 0; j < numFuncs; j++ ) {
      funcVals[j] += h*ftmp[j];
    }        
  }
  
  delete [] ftmp;
}


/*
  Add up the contributions to the total derivative from this stage/stage
*/
void TACSIntegrator::addToTotalDerivative( double scale, TACSBVec **adjoint ) {
  TacsScalar *tmp = new TacsScalar[ num_design_vars*num_funcs ];
  // Add the partial derivative w.r.t. the design variables
  memset(tmp, 0, num_design_vars*num_funcs*sizeof(TacsScalar));
  tacs->addDVSens(1.0, funcs, num_funcs, tmp, num_design_vars);
  for ( int m = 0; m < num_design_vars*num_funcs; m++ ){
    dfdx[m] += scale*tmp[m];
  }

  // Add the Adjoint-Residual partial derivative contributions
  memset(tmp, 0, num_design_vars*num_funcs*sizeof(TacsScalar));
  tacs->addAdjointResProducts(1.0, adjoint, num_funcs, tmp, num_design_vars);
  for ( int m = 0; m < num_design_vars*num_funcs; m++ ){
    dfdx[m] += scale*tmp[m];
  }
  delete [] tmp;
}

/*
  Solves the linear system Ax=b using LAPACK. The execution should be
  in serial mode.
*/
void TACSIntegrator::lapackLinearSolve( TACSBVec *res, TACSMat *mat, TACSBVec *update ) {
  // Serial only usage for debugging
  // Get the right hand side as an array
  TacsScalar *R;
  res->getArray(&R);
      
  // The following code retrieves a dense column-major 
  // matrix from the FEMat matrix
  TacsScalar *J;
  FEMat *femat = dynamic_cast<FEMat*>(mat);
  if (femat){
    int bsize, nrows;
    BCSRMat *B;
    femat->getBCSRMat(&B, NULL, NULL, NULL);
    B->getArrays(&bsize, &nrows, NULL, NULL, NULL, NULL);
	
    J = new TacsScalar[ bsize*bsize*nrows*nrows ];

    // Get the matrix in LAPACK column major order
    B->getDenseColumnMajor(J);
  }
  else {
    fprintf(stderr,"Casting to FEMat failed\n");
    exit(-1);
  }

  // Perform the linear solve
  int size = num_state_vars;
  int *dpiv = new int[size];
  int info = 0, one = 1;

  // Perform LU factorization
  LAPACKgetrf(&size, &size, J, &size, dpiv, &info);
  if (info){
    fprintf(stderr,"LAPACK GETRF output error %d\n", info);
    if (info < 0) {
      fprintf(stderr,"LAPACK GETRF: %d-th argument had an illegal value\n", info);
    } else {
      fprintf(stderr,"LAPACK GETRF: The factorization has been completed, \
but the factor U(%d,%d) is exactly singular, and division by zero will occur \
if it is used to solve a system of equations\n", info, info);
    }
    exit(-1);
  } 
  
  // Apply factorization
  LAPACKgetrs("N", &size, &one, J, &size, dpiv, R, &size, &info);
  if (info){
    fprintf(stderr,"LAPACK GETRS output error %d\n", info);
    exit(-1);
  }

  delete [] dpiv;

  // Set the lapack solution into the distributed vector
  TacsScalar *resvals;
  update->getArray(&resvals);
  memcpy(resvals, R, num_state_vars*sizeof(TacsScalar));

  if (J) { delete [] J; }
}

/*
  Function that returns a buffer that is formatted with the supplied
  format specifiers and arguments. Note that this function takes
  variable number of arguments. The argument list and format needs to
  be consistent.
*/
void TACSIntegrator::getString(char *buffer, const char * format, ... ){
  va_list args;
  va_start(args, format);
  vsprintf(buffer, format, args); 
  va_end(args);
}

/*
  Implement all the tasks to perform during each time step
*/
void TACSIntegrator::doEachTimeStep( int current_step ) {
  // Write the tecplot output to disk if sought
  writeStepToF5(current_step);    

  if ( current_step == 0) {
    // Log information
    if ( logfp && print_level >= 1){
      fprintf(logfp, "Variables=\n%12s %8s %12s %12s %15s %15s %15s\n", 
              "time", "NItrs", "|R|", "|R|/|R0|",
              "KineticEnrgy", "PotentialEnrgy",
              "EInit-E");
      
      // Compute the initial energy
      tacs->evalEnergies(&energies[0], &energies[1]);
      init_energy = energies[0] + energies[1];

      // Log the details
      fprintf(logfp, "%12.5e %8d %12.5e %12.5e %15.7e %15.7e %15.7e\n",
              time[0], 0, 0.0, 1.0,
              RealPart(energies[0]), RealPart(energies[1]),  0.0);
    }
  } else {
    // Print out the time step summary
    if ( logfp && print_level >= 1){
      tacs->evalEnergies(&energies[0], &energies[1]);
      fprintf(logfp, "%12.5e %8d %12.5e %12.5e %15.7e %15.7e %15.7e\n",
	      time[current_step], niter, RealPart(norm), RealPart(norm/(rtol + init_norm)),
	      RealPart(energies[0]), RealPart(energies[1]), 
	      RealPart((init_energy - (energies[0] + energies[1]))));
    }
  }
}

/*
  Implement all the tasks to perform during each nonlinear solve
*/
void TACSIntegrator::doEachNonLinearIter( int iter_num ) {}

/*
  Use this to set the values of class variables. Might be handy to
  pass a python dictionary/map containing the parameters read from an
  input file/user supplied. Yet to be implemented.
*/
void TACSIntegrator::setParameters( ... ){}

/*
  Set the relative tolerance for GMRES solver
*/
void TACSIntegrator::setRelTol( double _rtol ){ 
  rtol = _rtol; 
}

/*
  Set the absolute tolerance for GMRES solver
*/
void TACSIntegrator::setAbsTol( double _atol ){ 
  atol = _atol; 
}

/*
  Set the maximum number of allowed Newton iterations for the solution
  of nonlinear system
*/
void TACSIntegrator::setMaxNewtonIters( int _max_newton_iters ){
  max_newton_iters = _max_newton_iters;
}

/*
  Control the amount of information printed to the console and the
  logging stream
*/
void TACSIntegrator::setPrintLevel( int _print_level, char *_logfilename ){ 
  print_level = _print_level;
  if ( mpiRank == 0 && _logfilename && !( strlen(_logfilename) == 0) ) {
    // Set the log file name
    logfilename = _logfilename;
    
    // Close any opened non stdout logstreams
    if (logfp != stdout && logfp){
      fclose(logfp);
    }
    
    // Open a new file for logstream
    logfp = fopen(logfilename, "w");
  }
}

/*
  Number of times the Jacobian is recomputed/assembled during each
  nonlinear solve
*/
void TACSIntegrator::setJacAssemblyFreq( int _jac_comp_freq ){ 
  jac_comp_freq = _jac_comp_freq; 
}

/*
  Set whether or not to use LAPACK for linear solve
*/
void TACSIntegrator::setUseLapack( int _use_lapack ) {
  use_lapack = _use_lapack;
}

/*
  Set the functions of interest that take part in the adjoint solve.
*/
void TACSIntegrator::setFunction( TACSFunction **_funcs, int _num_funcs ) {
  num_funcs = _num_funcs;
  funcs = _funcs;
}

/*
  Set whether the LU factorization was done on the Jacobian
*/
void TACSIntegrator::setIsFactorized( int flag ){
  factorized = flag;
}

/*
  Update TACS states with the supplied ones (q, qdot, qddot)
*/
void TACSIntegrator::setTACSStates( double time, TACSBVec *q, TACSBVec *qdot, TACSBVec * qddot ){
  tacs->setSimulationTime(time);
  tacs->setVariables(q, qdot, qddot);
}

/*
  Constructor for BDF Integration scheme

  Input:
  tinit: the initial time
  tfinal: the final time
  num_steps_per_sec: the number of steps to take for each second
  max_bdf_order: global order of accuracy
*/
TACSBDFIntegrator::TACSBDFIntegrator( TACSAssembler * _tacs, 
                                      double _tinit, double _tfinal, 
                                      int _num_steps_per_sec, 
                                      int _max_bdf_order ):
TACSIntegrator(_tacs, _tinit,  _tfinal,  _num_steps_per_sec){		
  // copy over the variables
  max_bdf_order = _max_bdf_order;

  // Truncate the maximum order to 3rd order
  max_bdf_order = (max_bdf_order <= 3 ? 
		   max_bdf_order : 3);

  // Number of first and second order BDF coefficients
  nbdf = 0;
  nbddf = 0;
    
  // As many RHS as the number of second derivative coeffs
  num_adjoint_rhs = (2*max_bdf_order+1)+1;  
}

/*
  Destructor for TACSBDFIntegrator
*/
TACSBDFIntegrator::~TACSBDFIntegrator(){
}

/*
  Approximate states (q, qdot, qddot) at the current time step using
  the BDF coefficients and previous time step values of the states q,
  qdot and qddot.
  
  Input:
  pointers to the global states q, qdot, qddot
*/
void TACSBDFIntegrator::approxStates(){
  int k = current_time_step;

  // Zero the current states (these may not be zero when integrate() is
  // called for the second time)
  q[k]->zeroEntries();
  qdot[k]->zeroEntries();
  qddot[k]->zeroEntries();
  
  // get the BDF coefficients
  get2ndBDFCoeff(k, bdf_coeff, &nbdf, bddf_coeff, &nbddf, max_bdf_order);
  
  // Extrapolate to next time step: q[k] = q[k-1] + h*qdot[k-1] +
  // h^2/2*qddot[k-1] (helps reduce the initial residual)
  q[k]->copyValues(q[k-1]);
  q[k]->axpy(h, qdot[k-1]);
  q[k]->axpy(h*h/2.0, qddot[k-1]);

  // approximate qdot using BDF formula
  for ( int i = 0; i < nbdf; i++ ){
    double scale = bdf_coeff[i]/h;
    qdot[k]->axpy(scale, q[k-i]);
  }

  // approximate qddot using BDF formula
  for ( int i = 0; i < nbddf; i++ ){
    double scale = bddf_coeff[i]/(h*h);
    qddot[k]->axpy(scale, q[k-i]);
  }

  // If required, add the contribution to the second derivative
  // from the initial values of the first derivatives
  if (k == nbdf-1){
    double scale = bdf_coeff[nbdf-1]/h;
    qddot[k]->axpy(scale, qdot[0]);
  }
}

/*
  This code computes the BDF coefficients for first and second
  derivative approximations.
  
  input:
  k:         the integration time step
  max_order: the maximum order to use
  
  output:
  bdf:    the first derivative approximation
  nbdf:   the number first derivative of coefficients
  bddf:   the second derivative approximation
  nbddf:  the number second derivative of coefficients
*/
void TACSBDFIntegrator::get2ndBDFCoeff( const int k, 
					double bdf[], int *nbdf,
					double bddf[], int *nbddf,
					const int max_order ){
  // Construct the second-order BDF scheme
  memset(bddf, 0, (2*max_order+1)*sizeof(double));
  memset(bdf, 0, (max_order+1)*sizeof(double));

  // For the first time step, set the first coefficient to 1.0, 
  // but set the number of coefficients to zero
  if (k == 0){
    bdf[0] = 1.0;
    *nbdf = 0;
    *nbddf = 0;
    return;
  }

  // Get the first-order coefficients - one greater than the maximum
  // order of coefficients
  int order = (k < max_order ? k : max_order);
  *nbdf = getBDFCoeff(bdf, order)+1;
  *nbddf = 2*(*nbdf)-1;
  if (*nbddf > k+1){
    *nbddf = k+1;
  }

  // Now, compute the second-order coefficients
  for ( int j = 0; j < *nbdf && (k - j > 0); j++ ){
    // order >= 1 always
    int order = (k-j < max_order ? k-j : max_order);
    double bdf2[5];
    int len = getBDFCoeff(bdf2, order)+1;

    for ( int i = 0; i < len; i++ ){
      bddf[j + i] += bdf[j]*bdf2[i];
    }
  }
}

/*
  Setup the BDF coefficients
*/
void TACSBDFIntegrator::setupCoeffs(){
  
}

/*
  Get the first-order BDF coefficients of order <= 3 

  input:
  order:  order of the backwards-difference coefficients

  output: 
  bdf:    the backwards difference coefficients
*/ 
int TACSBDFIntegrator::getBDFCoeff( double bdf[], int order ){
  if (order <= 1){
    bdf[0] = 1.0;
    bdf[1] = -1.0;
    return 1;
  }
  else if (order == 2){
    bdf[0] =  1.5;
    bdf[1] = -2.0;
    bdf[2] =  0.5;
    return 2;
  }
  else if (order >= 3){
    bdf[0] =  35.0/24.0;
    bdf[1] = -15.0/8.0;
    bdf[2] =  3.0/8.0;
    bdf[3] =  1.0/24.0;
    return 3;
  }
  return 0;
}

/*
  Return the coefficients for linearizing the Residual using BDF method
*/
void TACSBDFIntegrator::getLinearizationCoeffs( double *alpha, double *beta, double *gamma ) {
  *gamma = bddf_coeff[0]/(h*h);
  *beta  = bdf_coeff[0]/h;
  *alpha = 1.0;
}

/*
  March backward in time and solve for the adjoint variables and find
  the total derivatives
*/
void TACSBDFIntegrator::marchBackwards( ) {
  // Adjoint variables for each function of interest
  TACSBVec **psi = new TACSBVec*[ num_funcs ];
  for ( int n = 0; n < num_funcs; n++ ){
    psi[n] = tacs->createVec();
    psi[n]->incref();
  }
  
  // Right hand sides for adjoint linear-system
  TACSBVec **rhs = new TACSBVec*[ num_funcs*num_adjoint_rhs ];
  for ( int n = 0; n < num_funcs*num_adjoint_rhs; n++ ){
    rhs[n] = tacs->createVec();
    rhs[n]->incref();
  }
 
  // March backwards in time (initial condition not evaluated)
  for ( int k = num_time_steps-1; k >=1 ; k--){
    current_time_step = k;
    
    // Get the BDF coefficients at this time step
    get2ndBDFCoeff(k, bdf_coeff, &nbdf, bddf_coeff, &nbddf, max_bdf_order);

    // Determine the linearization coefficients for Jacobian Assembly
    double gamma = bddf_coeff[0]/(h*h);
    double beta  = bdf_coeff[0]/h;
    double alpha = 1.0;

    // Set the stages
    setTACSStates(time[k], q[k], qdot[k], qddot[k]);
    
    // Find the adjoint index
    int adj_index = k % num_adjoint_rhs;

    // Setup the adjoint RHS
    TacsScalar ftmp;
    for ( int n = 0; n < num_funcs; n++ ){
      // Add the contribution to function value from this stage
      ftmp = 0.0;
      tacs->evalFunctions(&funcs[n], 1, &ftmp);
      fvals[n] += h*ftmp;

      // Add up the contribution from function state derivative to RHS
      psi[n]->zeroEntries();
      tacs->addSVSens(1.0, 0.0, 0.0, &funcs[n], 1, &psi[n]);

      // Add the contributions to the current adjoint RHS
      rhs[adj_index*num_funcs+n]->axpy(1.0, psi[n]);
      rhs[adj_index*num_funcs+n]->scale(-1.0);
    }
    
    // Setup the Jacobian
    tacs->assembleJacobian(alpha, beta, gamma, NULL, mat, TRANSPOSE);

    // LU factorization of the Jacobian
    pc->factor();
    
    // Apply the factorization for all right hand sides and solve for
    // the adjoint variables
    for ( int n = 0; n < num_funcs; n++ ){
      ksm->solve(rhs[adj_index*num_funcs + n], psi[n]);
      rhs[adj_index*num_funcs+n]->zeroEntries();
    }

    // Add total derivative contributions from this step to all
    // functions
    addToTotalDerivative(h, psi);

    // Drop the contributions from this step to other right hand sides
    for ( int ii = 1; (ii < nbdf || ii < nbddf); ii++ ){
      int rhs_index = (k - ii) % num_adjoint_rhs;
      double beta = 0.0, gamma = 0.0;
      if (ii < nbdf){
        beta = bdf_coeff[ii]/h;
      }
      if (ii < nbddf){
        gamma = bddf_coeff[ii]/(h*h);
      }
      for ( int n = 0; n < num_funcs; n++ ){      
	tacs->addJacobianVecProduct(1.0, 0.0, beta, gamma,
				    psi[n], rhs[rhs_index*num_funcs+n], TRANSPOSE);
      }
    }
  }
  // Freeup objects

  // Adjoint variables for each function of interest
  for ( int n = 0; n < num_funcs; n++ ){
    psi[n]->decref();
  }
  delete [] psi;

  // Right hand sides for adjoint linear-system
  for ( int n = 0; n < num_funcs*num_adjoint_rhs; n++ ){
    rhs[n]->decref();
  }
  delete [] rhs;
}

/*
  Constructor for TACSDIRKIntegrator

  Input:

  num_stages:        the number of Runge-Kutta stages
  tinit:             the initial time
  tfinal:            the final time
  num_steps_per_sec: the number of steps to take for each second
  num_stages:        the number of stages in DIRK
*/
TACSDIRKIntegrator::TACSDIRKIntegrator( TACSAssembler * _tacs, 
                                        double _tinit, double _tfinal, 
                                        int _num_steps_per_sec,
                                        int _num_stages ): 
TACSIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec){   
  // copy over the variables
  num_stages = _num_stages;
  
  // allocate space for stage state variables
  qS = new TACSBVec*[num_stages*num_time_steps];
  qdotS = new TACSBVec*[num_stages*num_time_steps];
  qddotS = new TACSBVec*[num_stages*num_time_steps];

  // store stage time
  tS = new double[num_stages*num_time_steps];
  memset(tS, 0, num_stages*num_time_steps*sizeof(double));
  
  // create state vectors for TACS during each timestep
  for ( int k = 0; k < num_stages*num_time_steps; k++ ){
    qS[k] = tacs->createVec(); 
    qS[k]->incref(); 
    
    qdotS[k] = tacs->createVec(); 
    qdotS[k]->incref(); 

    qddotS[k] = tacs->createVec(); 
    qddotS[k]->incref(); 
  }
  
  // Allocate space for Butcher tableau
  A = new double[num_stages*(num_stages+1)/2];
  B = new double[num_stages];
  C = new double[num_stages];
  
  // Set the Butcher tableau entries to zero
  memset(A, 0, num_stages*(num_stages+1)/2*sizeof(double));
  memset(B, 0, num_stages*sizeof(double));
  memset(C, 0, num_stages*sizeof(double));

  // Add entries into the Butcher tableau
  setupCoeffs();

  // As many RHS as the number of stages
  num_adjoint_rhs = num_stages;
}

/*
  Destructor for TACSDIRKIntegrator
*/
TACSDIRKIntegrator::~TACSDIRKIntegrator(){
  // Cleanup Butcher's Tableau
  delete [] A;
  delete [] B;
  delete [] C;

  // Cleanup stage stages
  for ( int i = 0; i < num_stages*num_time_steps; i++ ){
    qS[i]->decref();
    qdotS[i]->decref();
    qddotS[i]->decref();
  }

  // cleanup stage time
  delete [] tS;

  delete [] qS;
  delete [] qdotS;
  delete [] qddotS;
}

/*
  Function that puts the entries into Butcher tableau
*/
void TACSDIRKIntegrator::setupCoeffs(){
  if (num_stages == 1){
    // Implicit mid-point rule (A-stable)
    A[0] = 0.5;
    B[0] = 1.0;
    C[0] = 0.5;
    order = 2;
  } 
  else if (num_stages == 2){
    // Crouzeix formula (A-stable)
    double tmp = 1.0/(2.0*sqrt(3.0));

    A[0] = 0.5 + tmp;
    A[1] = -1.0/sqrt(3.0);
    A[2] = A[0];

    B[0] = 0.5;
    B[1] = 0.5;

    C[0] = 0.5 + tmp;
    C[1] = 0.5 - tmp;

    order = 3;     
  } 
  else if (num_stages == 3){
    // Crouzeix formula (A-stable)
    double alpha = 2.0*cos(M_PI/18.0)/sqrt(3.0);
    
    A[0] = (1.0 + alpha)*0.5;
    A[1] = -0.5*alpha;
    A[2] = A[0];
    A[3] = 1.0 + alpha;
    A[4] = -(1.0 + 2.0*alpha);
    A[5] = A[0];    

    B[0] = 1.0/(6.0*alpha*alpha);
    B[1] = 1.0 - 1.0/(3.0*alpha*alpha);
    B[2] = B[0];

    C[0] = 0.5*(1.0+alpha);
    C[1] = 0.5;
    C[2] = 0.5*(1.0-alpha);

    order = 4;
  }
  else {
    fprintf(stderr, "ERROR: Invalid number of stages %d\n", num_stages);
    num_stages = 1;
    setupCoeffs();
  }

  // check for the consistency of butcher tableau entries
  checkButcherTableau();
}

/*
  Function that checks the consistency of Butcher tableau values
*/
void TACSDIRKIntegrator::checkButcherTableau(){
  double tmp;

  // Check #1: sum(A(i,:)) = C(i)  
  int idx = -1;
  for ( int i = 0; i < num_stages; i++ ){
    tmp = 0.0;
    for ( int j = 0; j <= i; j++ ){
      idx++;
      tmp += A[idx];
    }

    // Check the difference
    if (fabs(C[i] - tmp) >= 1.0e-6) {
      fprintf(stderr, "WARNING: Sum A[%d,:] != C[%d] i.e. %f != %f \n", 
              i, i, C[i], tmp);
    }  
  }
  
  // Check #2: sum(B) = 1.0
  tmp = 0.0;
  for ( int i = 0; i < num_stages; i++ ){
    tmp += B[i];
  }
  if (fabs(1.0 - tmp) >= 1.0e-6) {
    fprintf(stderr, "WARNING: Sum B != 1.0 \n");
  }
}

/*
  Return the coefficients for linearizing the Residual using NBG method
*/
void TACSDIRKIntegrator::getLinearizationCoeffs( double *alpha, double *beta, double *gamma ) {
  int k   = current_time_step; // Retrieve the current time step
  int i   = current_stage;     // Retrieve the current stage
  int idx = getRowIdx(i);      // Starting entry of Butcher Tableau for this stage
  
  // Compute the coefficients
  *gamma = 1.0;
  *beta  = h*A[idx+i];           // h Aii
  *alpha = h*A[idx+i]*h*A[idx+i];  // h Aii h Aii
}

/*
  Function that computes the stage values at each time step. This
  function uses the global states and time at the previous time-step
  and sets the stage states tS, qS, qdotS and qdotS.
*/
void TACSDIRKIntegrator::approxStates(){
  int k = current_time_step;
  int i = current_stage;

  // Pointer to current stage
  int toffset = k*num_stages;

  // Zero the states (these may not be zero when integrate() is
  // called for the second time)
  qS[toffset+i]->zeroEntries();
  qdotS[toffset+i]->zeroEntries();
  qddotS[toffset+i]->zeroEntries();

  // Initial guess for qddotS
  if (i == 0){
    qddotS[toffset+i]->copyValues(qddot[k-1]);
  }
  else {
    qddotS[toffset+i]->copyValues(qddotS[toffset+i-1]);
  }

  // Compute qdotS
  qdotS[toffset+i]->copyValues(qdot[k-1]);
  int idx = getRowIdx(i);
  for ( int j = 0; j <= i; j++ ){
    qdotS[toffset+i]->axpy(h*A[idx+j], qddotS[toffset+j]);
  }

  // Compute qS
  qS[toffset+i]->copyValues(q[k-1]);
  idx = getRowIdx(i);
  for ( int j = 0; j <= i; j++ ){
    qS[toffset+i]->axpy(h*A[idx+j], qdotS[toffset+j]);
  }    
}

/*
  Start index of the Butcher Tableau A for the supplied stage
*/
int TACSDIRKIntegrator::getRowIdx( int stageNum ){
  return stageNum*(stageNum+1)/2;
}

/*
  Returns the coefficients for inter-stage contributions from stage
  adjoint variables. Note that there is no inter stage contributions
  for one-stage DIRK.
*/
void TACSDIRKIntegrator::getCoeffsInterStage( int current_stage, int target_stage, 
					      double *alpha, double *beta, double *gamma ){
  int i = current_stage;
  int j = target_stage;
  
  *gamma = 0.0;
  
  double h2 = h*h;

  if ( i == 2 &&  j == 0) { // contribution from third stage to first rhs

    *beta = B[2]*h*A[3];

    *alpha  = h2*B[2]*(A[0]*A[3]);
    *alpha += h2*B[2]*(A[1]*A[4]);
    *alpha += h2*B[2]*(A[3]*A[5]);

  } else if ( i == 1 &&  j == 0) { // contribution from second stage to first rhs

    *beta = B[1]*h*A[1];

    *alpha  = h2*B[1]*(A[0]*A[1]);
    *alpha += h2*B[1]*(A[1]*A[2]);
      
  } else if ( i == 2 &&  j == 1) {  // contribution from third stage to second rhs

    *beta = B[2]*h*A[4];

    *alpha  = h2*B[2]*(A[2]*A[4]);
    *alpha += h2*B[2]*(A[4]*A[5]);
        
  } else {
    printf("TACSIntegrator::Incorrect access to inter-stage coefficients...\n");
    exit(-1);
  }
}

/*
  Function that advances the global states q, qdot, qddot and time to
  next time step
  
  Input:
  The pointers to the global state variables q, qdot, qddot
  
  Output:
  Updated global state variables q, qdot, qddot at the time step
*/
void TACSDIRKIntegrator::computeTimeStepStates( int current_step, TACSBVec **q, TACSBVec **qdot, TACSBVec **qddot ){
  int k = current_step;
  int toffset  = num_stages*k;

  // Zero the states (these may not be zero when integrete() is
  // called for the second time)
  q[k]->zeroEntries();
  qdot[k]->zeroEntries();
  qddot[k]->zeroEntries();

  // advance the position state
  q[k]->copyValues(q[k-1]);
  for ( int j = 0; j < num_stages; j++ ){
    q[k]->axpy(h*B[j], qdotS[toffset+j]);
  }
  
  // advance the velocity state
  qdot[k]->copyValues(qdot[k-1]);
  for ( int j = 0; j < num_stages; j++ ){
    qdot[k]->axpy(h*B[j], qddotS[toffset+j]);
  }
  
  // advance the acceleration state
  for ( int j = 0; j < num_stages; j++ ){
    qddot[k]->axpy(B[j], qddotS[toffset+j]);
  }
}

/*
  Integration logic of DIRK. Use this function to march in time. The
  solution over time is set into the class variables q, qdot and qddot
  and time.
*/
void TACSDIRKIntegrator::integrate( ) {
  current_time_step = 0;
  
  // Get the initial condition
  tacs->getInitConditions(q[0], qdot[0]);
  tacs->setVariables(q[0], qdot[0]);

  // Perform logging, tecplot export, etc.
  doEachTimeStep(0);   

  for ( int k = 1; k < num_time_steps; k++ ){
    current_time_step = k;
    
    // Find the offset to current state values
    int toffset = k*num_stages;           

    // Compute the stage states qS, qdotS, qddotS based on the DIRK formula
    for ( int i = 0; i < num_stages; i++){
      current_stage = i;
      
      // Compute the stage time
      tS[toffset+i] = time[k-1] + C[i]*h;

      // Approximate the stage states using the DIRK formula
      approxStates();

      // Determine the coefficients for Jacobian assembly
      double alpha, beta, gamma;
      getLinearizationCoeffs(&alpha, &beta, &gamma);
      
      // Solve the nonlinear system of stage equations starting with the approximated states
      newtonSolve(alpha, beta, gamma, tS[toffset+i], 
		  qS[toffset+i], qdotS[toffset+i], qddotS[toffset+i]);
    }
    
    // Advance the time
    time[k] = time[k-1] + h;

    // Compute the state varialbes at the current time step using the
    // intermediate stage states
    computeTimeStepStates(k, q, qdot, qddot);

    // Perform logging, tecplot export, etc.
    doEachTimeStep(k);
  }
}
/*
  March backward in time and solve for the adjoint variables
*/
void TACSDIRKIntegrator::marchBackwards( ) {
  // Inter-step adjoint variables for each function of interest
  TACSBVec **psi = new TACSBVec*[ num_funcs ];
  TACSBVec **phi = new TACSBVec*[ num_funcs ];

  TACSBVec **psiTmp = new TACSBVec*[ num_funcs ];
  TACSBVec **phiTmp = new TACSBVec*[ num_funcs ];

  for ( int n = 0; n < num_funcs; n++ ){
    psi[n] = tacs->createVec();
    psi[n]->incref();

    phi[n] = tacs->createVec();
    phi[n]->incref();

    psiTmp[n] = tacs->createVec();
    psiTmp[n]->incref();
    
    phiTmp[n] = tacs->createVec();
    phiTmp[n]->incref();
  }
  
  // Stage right hand sides and adjoint variables for each function of interest
  TACSBVec **lambda = new TACSBVec*[ num_funcs*num_adjoint_rhs ];
  TACSBVec **rhs = new TACSBVec*[ num_funcs*num_adjoint_rhs ];
  TACSBVec **dfdq = new TACSBVec*[ num_funcs*num_adjoint_rhs ];
  for ( int n = 0; n < num_funcs*num_adjoint_rhs; n++ ){
    lambda[n] = tacs->createVec();
    lambda[n]->incref();

    rhs[n] = tacs->createVec();
    rhs[n]->incref();

    dfdq[n] = tacs->createVec();
    dfdq[n]->incref();
  }

  // March backwards in time (initial condition not evaluated)
  for ( int k = num_time_steps-1; k >= 1; k-- ) {
    current_time_step = k;
    
    // March backwards in stage
    int toffset = k*num_stages;

    // Set the pointer to the stage values
    double *ts = &tS[toffset];
    TACSBVec **qs = &qS[toffset];
    TACSBVec **qdots = &qdotS[toffset];
    TACSBVec **qddots = &qddotS[toffset];
    
    for ( int i = num_stages-1; i >= 0; i-- ) {
      current_stage = i;

      //-------------------------------------------------------------//
      // Solve for the adjoint variables 
      //-------------------------------------------------------------//

      // Compute the time
      double t = time[k-1] + C[i]*h;

      // Get the starting index for the corresponding Butcher's tableau row
      int idx = getRowIdx(i);

      // Determine the coefficients for linearizing the Residual
      double gamma = B[i];
      double beta  = B[i]*h*A[idx+i]; 
      double alpha = beta*h*A[idx+i];

      // Set the stages
      setTACSStates(ts[i], qs[i], qdots[i], qddots[i]);

      //--------------------------------------------------------------//
      // Assemble the right hand side
      //--------------------------------------------------------------//
      
      // Add function contributions
      TacsScalar ftmp;
      for ( int n = 0; n < num_funcs; n++ ){
        // Add up the function contribution (dfdq)
        ftmp = 0.0;
	tacs->evalFunctions(&funcs[n], 1, &ftmp);
        fvals[n] += h*B[i]*ftmp;

        // Add up the contribution from its derivative to RHS (drdq.lam)
	dfdq[i*num_funcs+n]->zeroEntries();
        tacs->addSVSens(1.0, 0.0, 0.0, &funcs[n], 1, &dfdq[i*num_funcs+n]);
	rhs[i*num_funcs+n]->axpy(alpha, dfdq[i*num_funcs+n]);

        // Add up the contribution from PSI to this RHS
        rhs[i*num_funcs+n]->axpy(B[i], phi[n]);

        // Add up the contribution from PSI to this RHS
        TacsScalar scale = 0.0;
        for ( int jj = i; jj < num_stages; jj++ ){
          idx = getRowIdx(jj);
          scale += B[jj]*h*A[idx+i];
        }
        rhs[i*num_funcs+n]->axpy(scale, psi[n]);
      
        // Negate the RHS
	rhs[i*num_funcs+n]->scale(-1.0);
      }
     
      // Setup the Jacobian
      tacs->assembleJacobian(alpha, beta, gamma, NULL, mat, TRANSPOSE);

      // LU factorization of the Jacobian
      pc->factor();
    
      // Apply the factorization for all right hand sides and solve for
      // the adjoint variables
      for ( int n = 0; n < num_funcs; n++ ){
	ksm->solve(rhs[i*num_funcs+n], lambda[i*num_funcs+n]);
	rhs[i*num_funcs+n]->zeroEntries();
      }

      // Add total derivative contributions from this step to all
      // functions
      addToTotalDerivative(h*B[i], &lambda[i*num_funcs]);
      
      //--------------------------------------------------------------//
      // Put the contributions from this stage to the right hand sides
      // of upcoming stages
      //--------------------------------------------------------------//
      
      for ( int j = i-1; j >= 0; j-- ){
	// Determine the corresponding coefficients
	getCoeffsInterStage(i, j, &alpha, &beta, &gamma);
	for ( int n = 0; n < num_funcs; n++ ){
	  // Put the contributions to the rhs from the function state
	  // variable sensitivity
	  rhs[j*num_funcs+n]->axpy(alpha, dfdq[i*num_funcs+n]);
	
	  // Put the contributions from Adjoint-Residual state variable
	  // sensitivity
	  tacs->addJacobianVecProduct(1.0, alpha, beta, gamma, 
				      lambda[i*num_funcs+n], rhs[j*num_funcs+n], 
				      TRANSPOSE);
	}
      }

      for ( int n = 0; n < num_funcs; n++ ){

        //---------------------------------------------------------//
        // Add the contributions to PHI from this stage
        //---------------------------------------------------------//

        gamma = 0.0;
        beta  = h*B[i];
        alpha = beta*h*C[i];

        // Add function SV sens contributions
        phiTmp[n]->axpy(alpha, dfdq[i*num_funcs+n]);
	  
        // Add the contributions from Adjoint-Residual state variable
        // sensitivity
        tacs->addJacobianVecProduct(1.0, alpha, beta, gamma,
                                    lambda[i*num_funcs+n], phiTmp[n], TRANSPOSE);

        //---------------------------------------------------------//
        // Add the contributions to PSI from this stage
        //---------------------------------------------------------//

        gamma = 0.0; beta = 0.0; alpha = h*B[i];

        // Add function SV sens contributions
        psiTmp[n]->axpy(alpha, dfdq[i*num_funcs+n]);
      
        // Add the contributions from Adjoint-Residual state variable
        // sensitivity
        tacs->addJacobianVecProduct(1.0, alpha, beta, gamma,
                                    lambda[i*num_funcs+n], psiTmp[n], TRANSPOSE);
        
      }
    } // stage

    // Add up the remaining PHI contributions
    for ( int n = 0; n < num_funcs; n++ ){
      phi[n]->axpy(h, psi[n]);       // phi = phi + h * psi
      phi[n]->axpy(1.0, phiTmp[n]);  // phi = phi + phiTmp
      phiTmp[n]->zeroEntries();       // Reset the stage accumulations now
    }

    // Add up the remaining PSI contributions
    for ( int n = 0; n < num_funcs; n++ ){
      psi[n]->axpy(1.0, psiTmp[n]);  // psi = psi + psiTmp
      psiTmp[n]->zeroEntries();       // Reset the stage contributions now
    }
 
  } // time
  // Freeup objects
  for ( int n = 0; n < num_funcs; n++ ){
    psi[n]->decref();
    phi[n]->decref();
    psiTmp[n]->decref();
    phiTmp[n]->decref();
  }

  for ( int n = 0; n < num_funcs*num_stages; n++ ){
    lambda[n]->decref();
    rhs[n]->decref();
    dfdq[n]->decref();
  }

  delete [] psi;
  delete [] psiTmp;
  delete [] phi;
  delete [] phiTmp;
  delete [] lambda;
  delete [] rhs;
  delete [] dfdq;
}

/*
  Evaluate the time average of functions: F = \sum_{k=0}^N h_k \sum_{j=0}^s b_j f(q_{k,j},t)
*/
void TACSDIRKIntegrator::evalTimeAvgFunctions( TACSFunction **funcs, 
					       int numFuncs, 
					       TacsScalar *funcVals) {
  memset(funcVals, 0, numFuncs*sizeof(TacsScalar));

  TacsScalar *ftmp = new TacsScalar[numFuncs];  
  memset(ftmp, 0, numFuncs*sizeof(TacsScalar));
  
  // Loop over time steps
  for ( int k = 1; k < num_time_steps; k++ ) {
    int toffset = k*num_stages;

    // Loop over all stages
    for ( int j = 0; j < num_stages; j++ ){

      // Set the states into TACS
      setTACSStates(tS[toffset+j], qS[toffset+j], qdotS[toffset+j], qddotS[toffset+j]);
      
      // Evaluate the functions
      tacs->evalFunctions(funcs, numFuncs, ftmp); 
      
      // Compute the mean function value (h is the weight)
      for ( int i = 0; i < numFuncs; i++ ) {
	funcVals[i] += h*B[j]*ftmp[i];
      }        
    }
  }
  delete [] ftmp;
}

/*
  Constructor for ABM Integration scheme

  Input:
  tinit: the initial time
  tfinal: the final time
  num_steps_per_sec: the number of steps to take for each second
  max_abm_order: the maximum global order of accuracy
*/
TACSABMIntegrator::TACSABMIntegrator( TACSAssembler * _tacs, 
                                      double _tinit, double _tfinal, 
                                      int _num_steps_per_sec, 
                                      int _max_abm_order ):
TACSIntegrator(_tacs, _tinit,  _tfinal,  _num_steps_per_sec){		
  // copy over the variables
  max_abm_order = _max_abm_order;

  // Sanity check on the input order
  max_abm_order = ((max_abm_order <= 6 && max_abm_order >= 1) ?
		   max_abm_order : 2);
  
  // Setup ABM Coefficents
  A = new double[max_abm_order*(max_abm_order+1)/2];
  memset(A, 0, max_abm_order*(max_abm_order+1)/2*sizeof(double)); // Lower triangular matrix
  setupCoeffs();
  checkABMCoeffs();

  // As many RHS as the number of second derivative coeffs
  num_adjoint_rhs = max_abm_order + 1;
}

/*
  Destructor for TACSABMIntegrator
*/
TACSABMIntegrator::~TACSABMIntegrator(){
  delete [] A;
}

/*
  Setup the coefficients needed for the ABM integration based on the
  maximum requested order and stores in the array
*/
void TACSABMIntegrator::setupCoeffs( ) {
  for ( int order = 1; order <= max_abm_order; order++ ){
    // Set the coefficients
    if ( order == 1 ) {  
      A[0] = 1.0;
    } else if ( order == 2 ) {
      A[1] = 1.0/2.0;
      A[2] = 1.0/2.0;
    } else if ( order == 3 ) {
      A[3] = 5.0/12.0;
      A[4] = 8.0/12.0;
      A[5] = -1.0/12.0;
    } else if ( order == 4 ) {
      A[6] = 9.0/24.0;
      A[7] = 19.0/24.0;
      A[8] = -5.0/24.0;
      A[9] = 1.0/24.0;
    } else if ( order == 5 ) {
      A[10] =  251.0/720.0;
      A[11] =  646.0/720.0;
      A[12] = -264.0/720.0;
      A[13] =  106.0/720.0;
      A[14] = -19.0/720.0;
    } else if ( order == 6 ) {
      A[15] =  475.0/1440.0;
      A[16] = 1427.0/1440.0;
      A[17] = -798.0/1440.0;
      A[18] =  482.0/1440.0;
      A[19] = -173.0/1440.0;
      A[20] =   27.0/1440.0;
    } else {
      fprintf(stderr, "WARNING: Wrong ABM order %d out of %d\n", order, max_abm_order);
    }
  }
}

/*
  Approximate states (q, qdot, qddot) at the current time step using
  the ABM coefficients and previous time step values of the states q,
  qdot and qddot. The accelration states are the unknowns from which
  the velocity states and position states are obtained.
  
  input:
  pointers to the global states q, qdot, qddot

  output:
  the state vectors (q, qdot, qddot) are prepared for nonlinear solve
*/
void TACSABMIntegrator::approxStates(){
  int k    = current_time_step;
  int m    = getOrder(k);
  int ridx = getRowIdx(m-1);

  // Zero the current states (these may not be zero when integrate() is
  // called for the second time)
  q[k]->zeroEntries();
  qdot[k]->zeroEntries();
  qddot[k]->zeroEntries();
  
  // Initial guess for qddot -- copy over previous accelearation states
  qddot[k]->copyValues(qddot[k-1]);

  // Approximate qdot using qdot and qddot:
  // qdot[k] = qdot[k-1] + \sum _{i=0}^{m-1} h A_i qddot[k-i]
  qdot[k]->copyValues(qdot[k-1]);
  for ( int i = 0; i <= m-1; i++ ){
    qdot[k]->axpy(h*A[ridx+i], qddot[k-i]);
  }

  // Approximate q using q and qdot:
  // q[k] = q[k-1] + \sum _{i=0}^{m-1} h A_i qdot[k-i]
  q[k]->copyValues(q[k-1]);
  for ( int i = 0; i <= m-1; i++ ){
    q[k]->axpy(h*A[ridx+i], qdot[k-i]);
  }
}

/*
  Return the  coefficients for linearizing the Residual using ABM method
*/
void TACSABMIntegrator::getLinearizationCoeffs( double *alpha, double *beta, double *gamma ) {
  int k   = current_time_step; // Retrieve the current time step
  int m   = getOrder(k);       // Determine the order of approximation
  int idx = getRowIdx(m-1);    // Order starts with 1 but the table starts at 0
  
  // Compute the coefficients
  *gamma = 1.0;
  *beta  = h*A[idx]; 
  *alpha = h*A[idx]*h*A[idx];
}

/*
  March backwards in time to solve for adjoint variables and computing
  total derivatives
*/
void TACSABMIntegrator::marchBackwards( ){
  int num_adjoint_rhs = 2; // ABM currently defaulted to one step

  // Adjoint variables
  TACSBVec **psi    = new TACSBVec*[ num_funcs ];
  TACSBVec **phi    = new TACSBVec*[ num_funcs ];
  TACSBVec **lambda = new TACSBVec*[ num_funcs ];

  TACSBVec **dfdq   = new TACSBVec*[ num_funcs ];
  TACSBVec **rhs    = new TACSBVec*[ num_funcs*num_adjoint_rhs ];

  for ( int n = 0; n < num_funcs*num_adjoint_rhs; n++ ){
    if (n < num_funcs) {    
      psi[n] = tacs->createVec();
      psi[n]->incref();

      phi[n] = tacs->createVec();
      phi[n]->incref();
      
      lambda[n] = tacs->createVec();
      lambda[n]->incref();
      
      dfdq[n] = tacs->createVec();
      dfdq[n]->incref();
    }
    rhs[n] = tacs->createVec();
    rhs[n]->incref();
  }

  // March backwards in time (initial condition not evaluated)
  for ( int k = num_time_steps-1; k >=1 ; k-- ){
    current_time_step = k;

    // Determine the coefficients for Jacobian assembly
    double alpha, beta, gamma;
    getLinearizationCoeffs(&alpha, &beta, &gamma);

    // Set the stages
    setTACSStates(time[k], q[k], qdot[k], qddot[k]);
    
    // Find the adjoint index
    int adj_index = k % num_adjoint_rhs;

    //---------------------------------------------------------------//
    // Setup the adjoint RHS
    //---------------------------------------------------------------//

    TacsScalar ftmp;
    for ( int n = 0; n < num_funcs; n++ ){
      // Evaluate the function
      tacs->evalFunctions(&funcs[n], 1, &ftmp);
      fvals[n] += h*ftmp;

      // Add up the contribution from function state derivative to RHS
      dfdq[n]->zeroEntries();
      tacs->addSVSens(1.0, 0.0, 0.0, &funcs[n], 1, &dfdq[n]);

      // Add the contributions to the current adjoint RHS
      rhs[adj_index*num_funcs+n]->axpy(alpha, dfdq[n]);

      rhs[adj_index*num_funcs+n]->scale(-1.0);
    }

    // Setup the Jacobian
    tacs->assembleJacobian(alpha, beta, gamma, NULL, mat, TRANSPOSE);

    // LU factorization of the Jacobian
    pc->factor();
    
    // Apply the factorization for all right hand sides and solve for
    // the adjoint variables
    for ( int n = 0; n < num_funcs; n++ ){
      ksm->solve(rhs[adj_index*num_funcs+n], lambda[n]);
      rhs[adj_index*num_funcs+n]->zeroEntries();
    }

    // Add total derivative contributions from this step for all
    // functions
    addToTotalDerivative(h, lambda);

    //-------------------------------------------------------------//
    // Put the contribution from this step to the next adjoint RHS //
    //-------------------------------------------------------------//

    for ( int ii = 1; ii < num_adjoint_rhs ; ii++ ){
      int rhs_index = (k - ii) % num_adjoint_rhs;
      for ( int n = 0; n < num_funcs; n++ ){

        //-----------------------------------------------------------//        
        // Find the new PHI
        //-----------------------------------------------------------//

        gamma = 0.0; beta  = 0.0; alpha = h;
        
        phi[n]->axpy(alpha, dfdq[n]);
        tacs->addJacobianVecProduct(1.0, alpha, beta, gamma, 
                                    lambda[n], phi[n], 
                                    TRANSPOSE);

        //-----------------------------------------------------------//
        // Find the new PSI
        //-----------------------------------------------------------//

        psi[n]->axpy(h*A[0], phi[n]);

        gamma = 0.0; beta  = h;  alpha = h*h*A[0];

        psi[n]->axpy(alpha, dfdq[n]);
        tacs->addJacobianVecProduct(1.0, alpha, beta, gamma, 
                                    lambda[n], psi[n], 
                                    TRANSPOSE);



        // Add contribution from psi
        rhs[rhs_index*num_funcs+n]->axpy(A[0], psi[n]);

        /*
          gamma = 0.0;
          beta  = 1.0/h;
          alpha = (0.5 + GAMMA);

          // Add function contribution
          rhs[rhs_index*num_funcs+n]->axpy(alpha, dfdq[n]);

          // Add Residual-Adjoint Product
          tacs->addJacobianVecProduct(1.0, alpha, beta, gamma, lambda[n], rhs[rhs_index*num_funcs+n], TRANSPOSE);
        */

        // Add the cumulative contributions
        //  rhs[rhs_index*num_funcs+n]->axpy(beta/h, phi[n]);

      }
    }
  }
  // Freeup objects
  for ( int n = 0; n < num_funcs*num_adjoint_rhs; n++ ){
    if (n < num_funcs) {    
      psi[n]->decref();
      phi[n]->decref();     
      lambda[n]->decref();
      dfdq[n]->decref();
    }
    rhs[n]->decref();
  }
  delete [] psi;
  delete [] phi;
  delete [] lambda;
  delete [] rhs;
  delete [] dfdq;
}

/*
  Start index of the the row of ABM coefficient matrix.
*/
int TACSABMIntegrator::getRowIdx( int row ){
  return row*(row+1)/2;
}

/*
  Start index of the the row of ABM coefficient matrix.
*/
int TACSABMIntegrator::getOrder( int k ){
  int order = k;
  if (order > max_abm_order ) order = max_abm_order;
  return order;
}

/*
  Sanity check on the ABM coefficient values 
*/
void TACSABMIntegrator::checkABMCoeffs(){
  for ( int order = 1; order <= max_abm_order; order++ ){
    int start_idx = getRowIdx(order-1);
    int end_idx = getRowIdx(order);
    double sum = 0.0;
    for ( int i = start_idx; i < end_idx; i++ ){
      sum += A[i];
    }
    if (abs(sum-1.0) > 1.e-15) {
      fprintf(stderr, "WARNING: Wrong ABM coefficients for \
order %d max order %d sum %f\n", order, max_abm_order, sum);
    }
  }
}

//----------------------------------------------------------//
//               Newmark Beta Gamma Method (NBG)            //  
//----------------------------------------------------------//

/*
  Constructor for NBG Integration scheme

  input:
  tacs              : the TACS assembler object
  tinit             : the initial time
  tfinal            : the final time
  num_steps_per_sec : the number of steps to take for each second
*/
TACSNBGIntegrator::TACSNBGIntegrator( TACSAssembler * _tacs, 
                                      double _tinit,
                                      double _tfinal, 
                                      int _num_steps_per_sec ):
TACSIntegrator(_tacs, _tinit,  _tfinal,  _num_steps_per_sec){
  // Setup the NBG coefficients
  setupCoeffs();
}

/*
  Destructor for TACSNBGIntegrator
*/
TACSNBGIntegrator::~TACSNBGIntegrator(){}

/*
  Setup the coefficients for the the integration scheme
*/
void TACSNBGIntegrator::setupCoeffs(){
  BETA  = 0.25;
  GAMMA = 0.50;
  
  //Other popular sets of values in Newmark family of integrators
  /*
  // Fox & Goodwin  (third order & conditionally stable rho=2.45)
  BETA  = 1.0/12.0;
  GAMMA = 0.50;

  // Linear Acceleration (second order & conditionally stable rho=3.46)
  BETA  = 1.0/6.0;
  GAMMA = 0.50;

  // Central Difference (second order & conditionally stable rho=2)
  BETA  = 1.0/2.0;
  GAMMA = 0.50;

  // Purely Explicit (firstorder & conditionally stable rho=0)
  BETA  = 0.0;
  GAMMA = 0.0;
  */
}

/*
  Approximate states (q, qdot, qddot) at the current time step using
  the NBG formula and previous time step values of the states q, qdot
  and qddot. The acceleration states are the unknowns from which the
  velocity states and position states are obtained.
  
  Use the appropriate class variables to prepare the states (q[k],
  qdot[k], qddot[k]) at the k-th timestep for nonlinear solution. Make
  sure to setup the coefficients, zero the states.
*/
void TACSNBGIntegrator::approxStates(){
  int k = current_time_step;

  double scale;
  // Zero the current states (these may not be zero when integrate() is
  // called for the second time)
  q[k]->zeroEntries();
  qdot[k]->zeroEntries();
  qddot[k]->zeroEntries();
  
  // Initial guess for qddot -- copy over previous acceleration states
  qddot[k]->copyValues(qddot[k-1]);

  // Approximate qdot using qdot and qddot:
  // qdot[k] = qdot[k-1] + (1-GAMMA) h qddot[k-1]) + GAMMA h qddot[k]
  qdot[k]->copyValues(qdot[k-1]);

  scale = (1.0-GAMMA)*h;
  qdot[k]->axpy(h*0.5, qddot[k-1]);

  scale = GAMMA*h;
  qdot[k]->axpy(h*0.5, qddot[k]);

  // Approximate q using q, qdot and qddot:
  // q[k] = q[k-1] + h qdot[k-1] + h^2(1-2*BETA)/2 qddot[k-1] + h^2 BETA qddot[k])
  q[k]->copyValues(q[k-1]);
  q[k]->axpy(h, qdot[k-1]);

  scale = h*h*(1.0-2.0*BETA)/2.0;
  q[k]->axpy(scale, qddot[k-1]);

  scale = h*h*BETA;
  q[k]->axpy(scale, qddot[k]);
}

/*
  Return the coefficients for linearizing the Residual using NBG method
*/
void TACSNBGIntegrator::getLinearizationCoeffs( double *alpha, double *beta, double *gamma ) {
  *gamma = 1.0;
  *beta  = GAMMA*h; 
  *alpha = BETA*h*h;
}


/*
  March backwards in time to solve for adjoint variables and computing
  total derivatives
*/
void TACSNBGIntegrator::marchBackwards( ){
  int num_adjoint_rhs = 2; // NBG is a one step method (uses
                           // information from previous and current
                           // steps)
  TACSBVec **psi    = new TACSBVec*[ num_funcs ];
  TACSBVec **phi    = new TACSBVec*[ num_funcs ];
  TACSBVec **lambda = new TACSBVec*[ num_funcs ];
  TACSBVec **dfdq   = new TACSBVec*[ num_funcs ];
  TACSBVec **rhs    = new TACSBVec*[ num_funcs*num_adjoint_rhs ];

  for ( int n = 0; n < num_funcs*num_adjoint_rhs; n++ ){
    if (n < num_funcs) {    
      psi[n] = tacs->createVec();
      psi[n]->incref();

      phi[n] = tacs->createVec();
      phi[n]->incref();
      
      lambda[n] = tacs->createVec();
      lambda[n]->incref();
      
      dfdq[n] = tacs->createVec();
      dfdq[n]->incref();
    }
    rhs[n] = tacs->createVec();
    rhs[n]->incref();
  }

  // March backwards in time (initial condition not evaluated)
  for ( int k = num_time_steps-1; k >=1 ; k-- ){
    current_time_step = k;
    
    // Determine the coefficients for Jacobian assembly
    double gamma = 1.0/(h*h);
    double beta  = GAMMA/h; 
    double alpha = BETA;

    // Set the stages
    setTACSStates(time[k], q[k], qdot[k], qddot[k]);
    
    // Find the adjoint index
    int adj_index = k % num_adjoint_rhs;

    //---------------------------------------------------------------//
    // Setup the adjoint RHS
    //---------------------------------------------------------------//

    TacsScalar ftmp;
    for ( int n = 0; n < num_funcs; n++ ){
      // Evaluate the function
      ftmp = 0.0;
      tacs->evalFunctions(&funcs[n], 1, &ftmp);
      fvals[n] += h*ftmp;

      // Add up the contribution from function state derivative to RHS
      dfdq[n]->zeroEntries();
      tacs->addSVSens(1.0, 0.0, 0.0, &funcs[n], 1, &dfdq[n]);

      // Add the contributions to the current adjoint RHS
      rhs[adj_index*num_funcs+n]->axpy(alpha, dfdq[n]);
      rhs[adj_index*num_funcs+n]->scale(-1.0);
    }

    // Setup the Jacobian
    tacs->assembleJacobian(alpha, beta, gamma, NULL, mat, TRANSPOSE);

    // LU factorization of the Jacobian
    pc->factor();
    
    // Apply the factorization for all right hand sides and solve for
    // the adjoint variables
    for ( int n = 0; n < num_funcs; n++ ){
      ksm->solve(rhs[adj_index*num_funcs+n], lambda[n]);
      rhs[adj_index*num_funcs+n]->zeroEntries();
    }

    // Add total derivative contributions from this step for all
    // functions
    addToTotalDerivative(h, lambda);

    //-------------------------------------------------------------//
    // Put the contribution from this step to the next adjoint RHS //
    //-------------------------------------------------------------//

    for ( int ii = 1; ii < num_adjoint_rhs ; ii++ ){
      int rhs_index = (k - ii) % num_adjoint_rhs;
      for ( int n = 0; n < num_funcs; n++ ){
        gamma = 0.0;
        beta  = 1.0/h;
        alpha = (0.5 + GAMMA);

        // Add function contribution
        rhs[rhs_index*num_funcs+n]->axpy(alpha, dfdq[n]);

        // Add Residual-Adjoint Product
        tacs->addJacobianVecProduct(1.0, alpha, beta, gamma, 
                                    lambda[n], rhs[rhs_index*num_funcs+n], 
                                    TRANSPOSE);

        // Add the cumulative contributions
        rhs[rhs_index*num_funcs+n]->axpy(beta/h, psi[n]);
        rhs[rhs_index*num_funcs+n]->axpy(alpha/h, phi[n]);
        
        // Find the new PSI
        gamma = 0.0; beta  = h;  alpha = h*h;

        psi[n]->axpy(h, phi[n]);
        psi[n]->axpy(alpha, dfdq[n]);
        tacs->addJacobianVecProduct(1.0, alpha, beta, gamma, 
                                    lambda[n], psi[n], 
                                    TRANSPOSE);
        
        // Find the new PHI
        gamma = 0.0; beta  = 0.0; alpha = h;

        phi[n]->axpy(alpha, dfdq[n]);
        tacs->addJacobianVecProduct(1.0, alpha, beta, gamma, 
                                    lambda[n], phi[n], 
                                    TRANSPOSE);

      }
    }
  }

  // Freeup objects
  for ( int n = 0; n < num_funcs*num_adjoint_rhs; n++ ){
    if (n < num_funcs) {    
      psi[n]->decref();
      phi[n]->decref();     
      lambda[n]->decref();
      dfdq[n]->decref();
    }
    rhs[n]->decref();
  }
  delete [] psi;
  delete [] phi;
  delete [] lambda;
  delete [] rhs;
  delete [] dfdq;
}
