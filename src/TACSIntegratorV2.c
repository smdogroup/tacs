#include "TACSIntegratorV2.h"
#include <math.h>
#include "tacslapack.h"

/*
  Base class constructor for integration schemes.

  input:
  tacs              : the tacs assembler object
  tInit             : the initial time
  tfinal            : the final time
  num_steps_per_sec : the number of steps to take for each second
*/
TACSIntegrator::TACSIntegrator( TACSAssembler * _tacs,
                                double _tinit, 
                                double _tfinal, 
                                double _num_steps_per_sec ){
  // Copy over the input parameters
  tacs = _tacs; 
  tacs->incref();
  tinit = _tinit;
  tfinal = _tfinal;
  num_steps_per_sec = _num_steps_per_sec;

  // MPI information
  MPI_Comm_rank(tacs->getMPIComm(), &mpiRank);
  MPI_Comm_size(tacs->getMPIComm(), &mpiSize);
  
  // Compute the step size
  h = 1.0/num_steps_per_sec;

  // compute the total number of time steps
  num_time_steps = int(num_steps_per_sec*(tfinal-tinit)) + 1;

  // Default print level and logging control
  print_level = 2;
  logfp = NULL;
  if (mpiRank == 0){ 
    logfp = stdout; 
  }

  // State variables that store the entire time history
  q = new TACSBVec*[num_time_steps];
  qdot = new TACSBVec*[num_time_steps];
  qddot = new TACSBVec*[num_time_steps];

  // create state vectors for TACS during each timestep
  for ( int k = 0; k < num_time_steps; k++ ) {
    q[k] = tacs->createVec(); q[k]->incref();
    qdot[k] = tacs->createVec(); qdot[k]->incref();
    qddot[k] = tacs->createVec(); qddot[k]->incref();
  }

  // Get the number of state variables and store into the class variable
  q[0]->getSize(&num_state_vars);
  
  // Store physical time of simulation
  time = new double[num_time_steps];
  memset(time, 0, num_time_steps*sizeof(double));

  //------------------------------------------------------------------//
  //                     Newton's method                              //
  //------------------------------------------------------------------//

  // Frequency of Jacobian recomputation during nonlinear solve
  jac_comp_freq = 1;

  // Set the default LINEAR solver
  use_lapack = 0;
  eigensolve = 0;
  use_femat = 1;

  // AMD reordering parameters
  lev = 100000; 
  fill = 10.0;  
  reorder_schur = 1;

  // KSM parameters
  gmres_iters  = 10;
  num_restarts = 0;
  is_flexible  = 0;
    
  // Default parameters for Newton solve
  max_newton_iters = 25;
  init_newton_delta = 0.0;
  atol = 1.0e-6;
  rtol = 1.0e-9;

  // Create vector for storing the residual at each Newton iteration
  res = tacs->createVec();
  res->incref();
  
  // Create a vector for storing the Newton updates
  update = tacs->createVec();
  update->incref();

  // NULL the different KSM/solver objects
  D = NULL;
  mat = NULL;
  pc = NULL;
  ksm = NULL;

  // Variables used in adjoint solve (use setFunction(...) to set these 
  num_funcs = 0;
  funcs = NULL;

  // Tecplot solution export
  f5_write_freq = 0; 
  f5_newton_freq = 0; 

  // Set the rigid and shell visualization objects to NULL
  rigidf5 = NULL;
  shellf5 = NULL;
  beamf5 = NULL;
  
  // Set kinetic and potential energies
  energies[0] = 0.0;
  energies[1] = 0.0;
  init_energy = 0.0;

  // Termination flag for newton solve
  newton_term = 0;

  // Force vector
  forces = NULL;
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
  if (D   && !isAdaptiveInstance()) { D->decref(); }
  if (mat && !isAdaptiveInstance()) { mat->decref(); }
  if (pc  && !isAdaptiveInstance()) { pc->decref(); }
  if (ksm && !isAdaptiveInstance()) { ksm->decref(); }
  
  if (time)     { delete [] time; }
  if (q)        { delete [] q; }
  if (qdot)     { delete [] qdot; }
  if (qddot)    { delete [] qddot; }

  // Dereference TACS
  if (tacs){ tacs->decref(); }
  
  if (rigidf5){ rigidf5->decref();}
  if (shellf5){ shellf5->decref();}
  if (beamf5){ beamf5->decref();}

  // External forces
  if (forces){forces->decref();}
}

/*
  March one step and exit the integration. This is primarily for use 
  with FUNtoFEM. 
*/
void TACSIntegrator::iterate( int k, TACSBVec *forces=NULL ){
  if (k == 0){
    // Perform initialization tasks
    initialize();

    // Retrieve the initial conditions and set into TACS
    tacs->getInitConditions(q[0], qdot[0], qddot[0]);
    tacs->setVariables(q[0], qdot[0], qddot[0]);

    // Output the results at the initial condition if configured
    doEachTimeStep(current_time_step);
    return;
  }
 
  // Advance time
  time[k] = time[k-1] + h;
  
  // Extrapolate to next time step: q[k] = q[k-1] + h*qdot[k-1] +
  // h^2/2*qddot[k-1] (helps reduce the initial residual)
  q[k]->copyValues(q[k-1]);
  q[k]->axpy(h, qdot[k-1]);
  q[k]->axpy(0.5*h*h, qddot[k-1]);

  // Get the BDF coefficients
  get2ndBDFCoeff(k, bdf_coeff, &nbdf, bddf_coeff, &nbddf, max_bdf_order);
  
  // approximate qdot using BDF formula
  qdot[k]->zeroEntries();
  for ( int i = 0; i < nbdf; i++ ){
    double scale = bdf_coeff[i]/h;
    qdot[k]->axpy(scale, q[k-i]);
  }

  // approximate qddot using BDF formula
  qddot[k]->zeroEntries();
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
  
  // Determine the coefficients for linearizing the Residual
  double alpha, beta, gamma;   
  getLinearizationCoeffs(&alpha, &beta, &gamma);
  
  // Solve the nonlinear system of stage equations starting with the
  // approximated states
  newton_term = newtonSolve(alpha, beta, gamma,
                            time[k], q[k], qdot[k], qddot[k],
                            forces, NULL);
  if (newton_term < 0){
    exit(-1);
  }

  // Tecplot output and print related stuff as configured
  doEachTimeStep(k);  
}

/*
  Integate forward in time using the initial conditions retrieved from
  TACS
*/
void TACSIntegrator::integrate( ){
  // Keep track of the time taken for foward mode
  double t0 = MPI_Wtime();
  time_forward = 0.0;
  time_fwd_assembly = 0.0;
  time_fwd_factor = 0.0;
  time_fwd_apply_factor = 0.0;
  time_newton = 0.0;

  // March forward in time
  printOptionSummary(logfp);
  for (int k = 0; k < num_time_steps; k++){
    iterate(k);
  }
  printOptionSummary(logfp);

  // Keep track of the time taken for foward mode
  time_forward += MPI_Wtime() - t0;
}

/*
  Compute the gradient for the given functions of interest with
  respect to the design variable.
*/
void TACSIntegrator::getFuncGrad( int _num_dv,
                                  TacsScalar *_x,
                                  TacsScalar *_fvals,
                                  TacsScalar *_dfdx ){
  // Copy the inputs
  num_design_vars = _num_dv;
  tacs->setDesignVars(_x, num_design_vars);

  //tacs->getDesignVars(xvals, num_design_vars);
  //MPI_Allreduce(MPI_IN_PLACE, xvals, num_design_vars, TACS_MPI_TYPE,
  //TACS_MPI_MAX, comm);

  // Check whether the function has been set properly
  if (num_funcs == 0 || funcs == NULL) {
    fprintf(stderr, "TACS Warning: Function is not set\n");
    return;
  }
  
  // Point the class variable to the memory allocated by the optimizer
  fvals = _fvals; 
  dfdx = _dfdx;
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
                                    double dh ){
  // Copy the inputs
  num_design_vars = _num_dv;
  tacs->setDesignVars(_x, num_design_vars);

  //tacs->getDesignVars(xvals, num_design_vars);
  //MPI_Allreduce(MPI_IN_PLACE, xvals, num_design_vars, TACS_MPI_TYPE,
  //TACS_MPI_MAX, comm);

  // Check whether the function has been set properly
  if (num_funcs == 0 || funcs == NULL) {
    fprintf(stderr, "TACS Warning: Function is not set, skipping adjoint solve. \n");
    return;
  }
  
  // Point the class variable to the memory allocated by the optimizer
  fvals = _fvals; 
  dfdx = _dfdx;
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
  this->evalFunctions(funcs, fvals, num_funcs);

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
    this->evalFunctions(funcs, fvals, num_funcs);

    // Evaluate the CS derivative
    for ( int j = 0; j < num_funcs; j++ ){
      dfdx[k+j*num_design_vars] = TacsImagPart(fvals[j])/dh;
    }

#else
    
    // Perturn the DV
    _x[k] = xtmp + dh;

    // Set the perturbed vector into TACS
    tacs->setDesignVars(_x, num_design_vars);
    
    // Integrate with perturbed x
    integrate();

    // Evaluate the functions at the current time
    this->evalFunctions(funcs, ftmp, num_funcs); 

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

  delete [] ftmp;
}

/*
  Drives the residual R(t,q,qdot,qddot) to zero using Newton's method

  Input: 
  The guessed (initial) state variable values q, qdot, qddot are
  supplied
  
  Output: q, qdot, qddot updated iteratively until the corresponding
  residual ||R|| <= tolerance
  
  alpha: multiplier for derivative of Residual wrt to q
  beta : multiplier for derivative of Residual wrt to qdot
  gamma: multiplier for derivative of Residual wrt to qddot

  forces: contains the additional contributions to the RHS
*/
int TACSIntegrator::newtonSolve( double alpha, double beta, double gamma,
                                 double t, TACSBVec *u, TACSBVec *udot, 
                                 TACSBVec *uddot, TACSBVec *forces,
                                 TACSBcMap *addBcs ){
  if (!mat || !ksm){
    // Set the D matrix to NULL
    if (use_femat){
      // Create a matrix for storing the Jacobian
      D = tacs->createFEMat(TACSAssembler::TACS_AMD_ORDER);
      D->incref();
      
      // Allocate the factorization
      pc = new PcScMat(D, lev, fill, reorder_schur);
      pc->incref();
      
      // Associate the maxtrix with FEMatrix
      mat = D;
      mat->incref();
    } 
    else {
      SerialBCSCMat *A = tacs->createSerialBCSCMat();
      pc = new SerialBCSCPc(A);
      pc->incref();
      
      mat = A;
      mat->incref();
      if (mpiSize > 1) {
        fprintf(stderr, "TACSIntegrator error: Using SerialBCSCMat in parallel\n");
      }
    }
  
    // The Krylov subspace method (KSM) associated with the solver
    ksm = new GMRES(mat, pc, gmres_iters, num_restarts, is_flexible);
    ksm->incref();
  }

  // ksm->setMonitor(new KSMPrintStdout("GMRES", 0, 1));
  ksm->setTolerances(rtol, atol);

  // Initialize the update norms
  update_norm = 1.0e99; 

  // Initialize the residual norms
  init_res_norm = 0.0;
  res_norm = 0.0;
  int newton_exit_flag = 0;

  if (logfp && print_level >= 2){
    fprintf(logfp, "%12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
            "time step", "#iters", "|R|", "|R|/|R0|", "|dq|", 
            "alpha", "beta", "gamma","delta");
  }

  // Track the time taken for newton solve at each time step
  double tnewton = MPI_Wtime();

  // Iterate until max iters or R <= tol
  double delta = 0.0;
  for ( niter = 0; niter < max_newton_iters; niter++ ){    
    // Set the supplied initial input states into TACS
    setTACSStates(t, u, udot, uddot);

    // Write the tecplot output to disk if sought
    writeNewtonIterToF5(current_time_step, niter);

    // Assemble the Jacobian matrix once in Newton iterations
    double t0 = MPI_Wtime();
    if ((niter % jac_comp_freq) == 0){
      delta = init_newton_delta*gamma;
      if (niter > 0 && (TacsRealPart(res_norm) < TacsRealPart(init_res_norm))){
        delta *= TacsRealPart(res_norm/init_res_norm);
      }

      tacs->assembleJacobian(alpha, beta, gamma + delta,
                             res, mat, NORMAL);
      if (addBcs){
        res->applyBCs(addBcs);
        mat->applyBCs(addBcs);
      }
    }
    else {
      tacs->assembleRes(res);
      if (addBcs){
        res->applyBCs(addBcs);
      }
    }

    // Add the forces into the residual    
    if (forces){
      tacs->applyBCs(forces);
      res->axpy(-1.0, forces);
      tacs->applyBCs(res);
      if (addBcs){
        res->applyBCs(addBcs);
      }
    }
    
    time_fwd_assembly += MPI_Wtime() - t0;

    // Compute the L2-norm of the residual
    res_norm = res->norm();
    
    // Record the residual norm at the first Newton iteration
    if (niter == 0){
      init_res_norm = res_norm;
    }

    // Write a summary    
    if(logfp && print_level >= 2){
      if (niter == 0){
        fprintf(logfp, "%12d %12d %12.5e %12.5e %12s %12.5e %12.5e %12.5e %12.5e\n",
                current_time_step, niter, TacsRealPart(res_norm),  
                (niter == 0) ? 1.0 : TacsRealPart(res_norm/init_res_norm), 
                " ", alpha, beta, gamma, delta);
      }
      else {
        fprintf(logfp, "%12d %12d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
                current_time_step, niter, TacsRealPart(res_norm),  
                (niter == 0) ? 1.0 : TacsRealPart(res_norm/init_res_norm), 
                TacsRealPart(update_norm),  
                alpha, beta, gamma, delta);
      }
    }

    // Check if the norm of the residuals is a NaN
    if (res_norm != res_norm || update_norm != update_norm ){ 
      if (logfp) {
        fprintf(stderr,"[%d] Newton iteration %d, failed with NaN residual norm\n", 
          mpiRank, niter);
      }
      newton_exit_flag = -2;
      break;
    }
    
    // Check whether the update is sufficiently small
    if (TacsRealPart(update_norm) < atol){
      newton_exit_flag = 2;
      break;
    }

    // Check if the Newton convergence tolerance is satisfied
    if (TacsRealPart(res_norm) < atol){
      newton_exit_flag = 1;
      break;
    }

    // Check for relative reduction in residual magnitude
    if (TacsRealPart(res_norm) < rtol*TacsRealPart(rtol + init_res_norm)){
      newton_exit_flag = 3;
      break;
    }

    if (use_lapack){
      if (mpiSize > 1){
        fprintf(stderr, "TACSIntegrator:: Using LAPACK in parallel!\n");
      }
      // Perform the linear solve using LAPACK (serial only)
      lapackLinearSolve(res, mat, update);
    }
    else {
      // LU Factor the matrix when needed
      double t1 = MPI_Wtime();
      if ((niter % jac_comp_freq) == 0){
        pc->factor();
      }
      time_fwd_factor += MPI_Wtime() - t1;      

      // Solve for update using KSM
      double t2 = MPI_Wtime();
      ksm->solve(res, update);
      time_fwd_apply_factor += MPI_Wtime() - t2;
    }

    // Solve for the eigenvalues      
    if (eigensolve){
      lapackEigenSolve(mat);
    }

    // Find the norm of the displacement update
    update_norm = update->norm()*alpha;

    // Update the state variables using the solution
    uddot->axpy(-gamma, update);
    udot->axpy(-beta, update);
    u->axpy(-alpha, update);
  }

  // Failed nonlinear solution
  if (niter == max_newton_iters) {
    newton_exit_flag = -1;
  }

  // Record the time taken for nonlinear solution
  time_newton = MPI_Wtime() - tnewton;

  // Return the exit flag
  return newton_exit_flag;
}

/*
  Function that writes time, q, qdot, qddot to file
*/
void TACSIntegrator::writeSolution( const char *filename, int format ) {
  FILE *fp = fopen(filename, "w");
  TacsScalar *qvals, *qdotvals, *qddotvals;

  if (format == 1){

    // Plain format with t q[0], q[1],...,qdot[0], qdot[1],...,qddot[0], qddot[1],...
    for ( int k = 0; k < num_time_steps; k++ ){    
      // Copy over the state values from TACSBVec
      q[k]->getArray(&qvals);
      qdot[k]->getArray(&qdotvals);
      qddot[k]->getArray(&qddotvals);
  
      // Write the time and states to file
      fprintf(fp, "%12.5e ", time[k]);  
      for ( int j = 0; j < num_state_vars; j++ ){
        fprintf(fp, "%12.5e %12.5e %12.5e ", TacsRealPart(qvals[j]), 
                TacsRealPart(qdotvals[j]), TacsRealPart(qddotvals[j]));
      }
      fprintf(fp, "\n");
    }

  } else {

    if (format > 0) {

      /*
        A little more easily readable formatting
        . time
        . DOF number
        . DOF of element
        . Element number
        . q
        . qdot
        . qddot
      */
      for ( int k = 0; k < num_time_steps; k++ ){       
        // Copy over the state values from TACSBVec
        q[k]->getArray(&qvals);
        qdot[k]->getArray(&qdotvals);
        qddot[k]->getArray(&qddotvals);

        fprintf(fp, "time=%e \n", time[k]);
    
        // Write the time and states to file
        int elem_ctr = 0;     
        for ( int j = 0; j < num_state_vars; j++ ){
          fprintf(fp, "%12d %3d %5d %12.5e %12.5e %12.5e \n", 
                  j, // global DOF number
                  j % 8, // DOF number of each element
                  (j % 8 == 7) ? (elem_ctr++) : elem_ctr, 
                  TacsRealPart(qvals[j]), // q
                  TacsRealPart(qdotvals[j]), // qdots
                  TacsRealPart(qddotvals[j])); // qddots
        }
        fprintf(fp, "\n");
      }

    } else {

      // Write the DOFS on user specified element number in final ordering
      
      // Plain format with t q[0], q[1],...,qdot[0], qdot[1],...,qddot[0], qddot[1],...
      for ( int k = 0; k < num_time_steps; k++ ){    
        // Copy over the state values from TACSBVec
        q[k]->getArray(&qvals);
        qdot[k]->getArray(&qdotvals);
        qddot[k]->getArray(&qddotvals);
  
        // Write the time and states to file
        fprintf(fp, "%12.5e ", time[k]);  
        int elem_ctr = 0;
        for ( int j = 0; j < num_state_vars; j++ ){
          // Write if we have found the sought element
          if (elem_ctr == -format ) {
            fprintf(fp, " %12.5e %12.5e %12.5e ",
                    TacsRealPart(qvals[j]), 
                    TacsRealPart(qdotvals[j]), 
                    TacsRealPart(qddotvals[j]));
          }
          if (j % 8 == 7) {
            elem_ctr++;
          }
        }
        fprintf(fp, "\n");
      }

    }
    
  }
  // Close the output file safely
  fclose(fp);
}


/*
  Writes the output as an f5 file. Note the states and time must be
  set appropriately before calling this function.
*/
void TACSIntegrator::writeNewtonIterToF5( int k, int n ){
  if(rigidf5 && getWriteFlag(k, f5_newton_freq)){
    // Create a buffer for filename 
    char rbuffer[256];

    // Format the buffer based on the time step
    getString(rbuffer, "results/rigid_%06d_%03d.f5", k, n);
    
    // Write the f5 file for this time step
    rigidf5->writeToFile(rbuffer);

  }

  if(shellf5 && getWriteFlag(k, f5_newton_freq)){

    // Create a buffer for shell filename 
    char sbuffer[256];

    // Format the buffer based on the time step
    getString(sbuffer, "results/shell_%06d_%03d.f5", k, n);
    
    // Write the f5 file for this time step
    shellf5->writeToFile(sbuffer);
  }

  if(beamf5 && getWriteFlag(k, f5_newton_freq)){

    // Create a buffer for beam filename 
    char sbuffer[256];

    // Format the buffer based on the time step
    getString(sbuffer, "results/beam_%06d_%03d.f5", k, n);
    
    // Write the f5 file for this time step
    beamf5->writeToFile(sbuffer);
  }
}

/*
  Writes the output as an f5 file. Note the states and time must be
  set appropriately before calling this function.
*/
void TACSIntegrator::writeStepToF5( int k ){
  if(rigidf5 && getWriteFlag(k, f5_write_freq)){
    // Create a buffer for filename 
    char rbuffer[256];

    // Format the buffer based on the time step
    getString(rbuffer, "results/rigid_%06d_%04d.f5", k, iter);
    
    // Write the f5 file for this time step
    rigidf5->writeToFile(rbuffer);
  }

  if(shellf5 && getWriteFlag(k, f5_write_freq)){
    // Create a buffer for shell filename 
    char sbuffer[256];

    // Format the buffer based on the time step
    getString(sbuffer, "results/shell_%06d_%04d.f5", k, iter);
    
    // Write the f5 file for this time step
    shellf5->writeToFile(sbuffer);
  }

  if(beamf5 && getWriteFlag(k, f5_write_freq)){
    // Create a buffer for beam filename 
    char sbuffer[256];

    // Format the buffer based on the time step
    getString(sbuffer, "results/beam_%06d.f5", k);
    
    // Write the f5 file for this time step
    beamf5->writeToFile(sbuffer);
  }
}

/*
  Creates an f5 file for each time step and writes the data.
*/
void TACSIntegrator::writeSolutionToF5( int _f5_write_freq ){
  /* Determine parallelism
     int is, ie, idec;
     idec = num_time_steps/mpiSize;
     is   = idec*mpiRank;
     ie   = idec*(mpiRank + 1);
     if ( mpiRank == mpiSize - 1) { ie = num_time_steps-1; }
  */
  // Loop through all timesteps
  for ( int k = 0; k < num_time_steps; k++ ){
    // Set the current states into TACS
    setTACSStates(time[k], q[k], qdot[k], qddot[k]);
      
    // Write RIGID body if set
    if(rigidf5 && getWriteFlag(k, _f5_write_freq)){
      char fname[128];
      getString(fname, "results/rigid_%06d.f5", k);
      rigidf5->writeToFile(fname);
    }

    // Write SHELL body if set
    if(shellf5 && getWriteFlag(k, _f5_write_freq)){
      char fname2[128];
      getString(fname2, "results/shell_%06d.f5", k);
      shellf5->writeToFile(fname2);
    }

    // Write BEAM body if set
    if(beamf5 && getWriteFlag(k, _f5_write_freq)){
      char fname2[128];
      getString(fname2, "results/beam_%06d.f5", k);
      beamf5->writeToFile(fname2);
    }
  }
}

/*
  Prints the wall time taken during operations in TACSIntegrator
   
  input:
  level: controls the level of detail requested in timing
  t0   : reference time to normalize the times calculated within TACSIntegrator
*/
void TACSIntegrator::printWallTime( double t0, int level ){
  if(level >= 0) { 
    fprintf(logfp, "[%d] Total            : %8.2f %6.2f\n", 
      mpiRank, t0, t0/t0); 
    fprintf(logfp, "[%d] Integrator       : %8.2f %6.2f\n", 
      mpiRank, time_forward + time_reverse, (time_forward +time_reverse)/t0); 
  }

  if (level >= 1) { 
    fprintf(logfp, ".[%d] Forward         :  %8.2f %6.2f\n", 
      mpiRank, time_forward, time_forward/t0); 
  }

  if (level >= 2) {
    fprintf(logfp, "..[%d] Assembly       :   %8.2f %6.2f\n", 
      mpiRank, time_fwd_assembly, time_fwd_assembly/t0);
    fprintf(logfp, "..[%d] Factor         :   %8.2f %6.2f\n", 
      mpiRank, time_fwd_factor, time_fwd_factor/t0);
    fprintf(logfp, "..[%d] ApplyFac       :   %8.2f %6.2f\n", 
      mpiRank, time_fwd_apply_factor, time_fwd_apply_factor/t0);
  }

  if (level >= 1) { 
    fprintf(logfp, ".[%d] Reverse         :  %8.2f %6.2f\n", 
      mpiRank, time_reverse, time_reverse/t0); 
  }

  if (level >= 2) {
    fprintf(logfp, "..[%d] Assembly       :   %8.2f %6.2f\n", 
      mpiRank, time_rev_assembly, time_rev_assembly/t0);
    fprintf(logfp, "..[%d] Factor         :   %8.2f %6.2f\n", 
      mpiRank, time_rev_factor, time_rev_factor/t0);
    fprintf(logfp, "..[%d] ApplyFac       :   %8.2f %6.2f\n", 
      mpiRank, time_rev_apply_factor, time_rev_apply_factor/t0);
    fprintf(logfp, "..[%d] JacVecPdt      :   %8.2f %6.2f\n", 
      mpiRank, time_rev_jac_pdt, time_rev_jac_pdt/t0);
  }
}

/*
  Configure the F5 output 
*/
void TACSIntegrator::setOutputFrequency( int _write_freq, int _newton_freq ){
  f5_write_freq = _write_freq;
  f5_newton_freq = _newton_freq;
}

/*
  Set whether RIGID body components are a part of the output
*/
void TACSIntegrator::setRigidOutput( int flag ){
  // Create an TACSToFH5 object for writing output to files
  unsigned int rigid_write_flag = (TACSElement::OUTPUT_NODES |  
                                   TACSElement::OUTPUT_DISPLACEMENTS);
  // Create a new instance if it does not exist and the flag is
  // greater than zero
  if (!rigidf5 && flag>0){
    // Create an instance of f5 and store to the class variable
    rigidf5 = new TACSToFH5(tacs, TACS_RIGID, rigid_write_flag);
    rigidf5->incref();
  }
}

/*
  Set whether SHELL body components are a part of the output
*/
void TACSIntegrator::setShellOutput( int flag ){
  // Create an TACSToFH5 object for writing output to files
  unsigned int shell_write_flag = (TACSElement::OUTPUT_NODES |
                                   TACSElement::OUTPUT_DISPLACEMENTS |
                                   TACSElement::OUTPUT_STRAINS |
                                   TACSElement::OUTPUT_STRESSES |
                                   TACSElement::OUTPUT_EXTRAS);
  // Create a new instance if it does not exist and the flag is
  // greater than zero
  if (!shellf5 && flag>0){
    shellf5 = new TACSToFH5(tacs, TACS_SHELL, shell_write_flag);
    shellf5->incref();
  }
}

/*
  Set whether BEAM components are a part of the output
*/
void TACSIntegrator::setBeamOutput( int flag ){
  // Create an TACSToFH5 object for writing output to files
  unsigned int beam_write_flag = (TACSElement::OUTPUT_NODES |
                                  TACSElement::OUTPUT_DISPLACEMENTS |
                                  TACSElement::OUTPUT_STRAINS |
                                  TACSElement::OUTPUT_STRESSES |
                                  TACSElement::OUTPUT_EXTRAS);
  // Create a new instance if it does not exist and the flag is
  // greater than zero
  if (!beamf5 && flag>0){
    beamf5 = new TACSToFH5(tacs, TACS_TIMOSHENKO_BEAM, beam_write_flag);
    beamf5->incref();
  }
}

/*
  Don't write if f5_write_freq is set to 0. If write output is sought,
  force writing initial(k=0) and final (k=num_time_steps-1) time steps
  as they are the most important outputs that an user would need
*/
int TACSIntegrator::getWriteFlag( int k, int f5_write_freq ){
  int write_now = 0;
  if (f5_write_freq > 0) {
    if ( k == 0 || k == num_time_steps-1 ){
      write_now = 1;
    } 
    else {
      write_now = (k % f5_write_freq == 0);
    }
  }
  return write_now;
}

/*
  Add up the contributions to the total derivative from this stage/stage
*/
void TACSIntegrator::addToTotalDerivative( double scale, TACSBVec **adjoint ) {
  // Add df/dx
  tacs->addDVSens(scale, funcs, num_funcs, dfdx, num_design_vars);

  // Add lam.dR/dx
  tacs->addAdjointResProducts(scale, adjoint, num_funcs, dfdx, num_design_vars);
}

/*
  Performs an eigen solve of the Jacobian matrix
 */
void TACSIntegrator::lapackEigenSolve( TACSMat *mat ) {    
  /*
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

  // Set up LAPACK call
  int size = num_state_vars;
  TacsScalar *wr = new TacsScalar[size];
  TacsScalar *wi = new TacsScalar[size];
  int lwork = 8*size;
  TacsScalar *work = new TacsScalar[lwork];
  TacsScalar *vl = new TacsScalar[size*size];
  TacsScalar *vr = new TacsScalar[size*size];
  int info = 0; 
  
  // Call lapack to solve the eigenvalue problem
  LAPACKdgeev("N", "V", &size, J, &size,
              wr, wi,
              vl, &size, vr, &size, 
              work, &lwork,
              &info);

  // Print the eigenvalues
  for (int i = 0; i < size; i++){
    printf("%d %12.5e %12.5e\n", i, wr[i], wi[i]);
  }

  // Print the eigenvector
  for (int i = 0; i < size; i++){
    printf("\n%d Eigenvector for %12.5e\n", i, wr[i]);
    for (int j = 0; j < size; j++){
      if (wr[i] > 0.0){
        printf("[%12.5e, %12.5e] ", vr[i*size+j], vr[i*size+j+1]);
      } else if (wr[i] < 0.0){
        printf("[%12.5e, %12.5e] ", vr[i*size+j-1], vr[i*size+j]);
      } else {
        printf("%12.5e ", vr[i*size+j]);
      }
    }
  }

  if (info){
    fprintf(stderr,"LAPACK DGEEV output error %d\n", info);
    if (info < 0){
      fprintf(stderr,"LAPACK : %d-th argument had an illegal value\n", info);
    } else {
      fprintf(stderr,"LAPACK DGEEV: The QR algorithm failed to compute \
all the eigenvalues, and no eigenvectors have been computed");
    }
    exit(-1);
  }
  
  delete [] wr;
  delete [] wi;
  delete [] work;
  delete [] vl;
  delete [] vr;
  if (J) { delete [] J; }
  */
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

  // Set the lapack solution into the distributed vector
  TacsScalar *ans;
  update->getArray(&ans);
  memcpy(ans, R, num_state_vars*sizeof(TacsScalar));
      
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
  LAPACKgetrs("N", &size, &one, J, &size, dpiv, ans, &size, &info);
  if (info){
    fprintf(stderr,"LAPACK GETRS output error %d\n", info);
    exit(-1);
  }

  delete [] dpiv;

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

  // Evaluate the energies
  tacs->evalEnergies(&energies[0], &energies[1]);

  if (current_step == 0){
    // Log information
    if (logfp && print_level >= 1){
      fprintf(logfp, "%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n", 
              "status", "time", "tnewton", "#iters", "NwtnFlg", "|R|", "|R|/|R0|", "|dq|",
              "KE", "PE", "E0-E", "|F|");
      
      // Compute the initial energy
      init_energy = energies[0] + energies[1];

      // Log the details
      fprintf(logfp, "%6d/%-6d %12.5e %12.5e %12d %12d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
              0, num_time_steps,
              time[0], time_newton, 0, 0, 0.0, 0.0, 0.0,
              TacsRealPart(energies[0]), TacsRealPart(energies[1]),  
              0.0, 
              (forces) ? TacsRealPart(forces->norm()):0.0);
    }
  } 
  else {
    // Print out the time step summary
    if (logfp && print_level >= 1){

      // Need a title for total summary as details of Newton iteration
      // will overshadow this one line summary
      if (print_level == 2){
        fprintf(logfp, "%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
                "status",
                "time", "tnewton", "#iters", "NwtnFlg", "|R|", "|R|/|R0|", "|dq|",
                "KE", "PE", "E0-E", "|F|");
      }
      
      fprintf(logfp, "%6d/%-6d %12.5e %12.5e %12d %12d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
              current_step, num_time_steps,
              time[current_step], time_newton,
              niter, newton_term,
              TacsRealPart(res_norm), 
              TacsRealPart(res_norm/(rtol + init_res_norm)),
              TacsRealPart(update_norm),
              TacsRealPart(energies[0]), TacsRealPart(energies[1]), 
              TacsRealPart((init_energy - (energies[0] + energies[1]))), 
              (forces) ? TacsRealPart(forces->norm()):0.0);
    }
  }
}

/*
  Implement all the tasks to perform during each nonlinear solve
*/
void TACSIntegrator::doEachNonLinearIter( int iter_num ) {}

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
  Write out all of the options that have been set to a output stream.
*/
void TACSIntegrator::printOptionSummary( FILE *fp ) {
  if (fp && print_level>0){
    fprintf(fp, "===============================================\n");
    fprintf(fp, "TACSIntegrator: Parameter values\n");
    fprintf(fp, "===============================================\n");
    fprintf(fp, "%-30s %15s\n", "integrator type", getIntegratorType(mytype));
    fprintf(fp, "%-30s %15g\n", "tinit", tinit );
    fprintf(fp, "%-30s %15g\n", "tfinal", tfinal );
    fprintf(fp, "%-30s %15g\n", "num_steps_per_sec", num_steps_per_sec );
    fprintf(fp, "%-30s %15g\n", "step size", h);
    fprintf(fp, "%-30s %15d\n", "num_time_steps", num_time_steps );
    fprintf(fp, "%-30s %15d\n", "num_state_vars", num_state_vars);
    fprintf(fp, "%-30s %15d\n", "mpiSize", mpiSize);
    fprintf(fp, "%-30s %15d\n", "print_level", print_level);

    fprintf(fp, "===============================================\n");
    fprintf(fp, "Nonlinear Solver: Parameter values\n");
    fprintf(fp, "===============================================\n");    
    fprintf(fp, "%-30s %15d\n", "max_newton_iters", max_newton_iters);
    fprintf(fp, "%-30s %15g\n", "absolute_tolerance", atol);
    fprintf(fp, "%-30s %15g\n", "relative_tolerance", rtol);
    fprintf(fp, "%-30s %15d\n", "jac_comp_freq", jac_comp_freq);
    fprintf(fp, "%-30s %15d\n", "use_line_search", use_line_search);

    fprintf(fp, "===============================================\n");
    fprintf(fp, "Linear Solver: Parameter values\n");
    fprintf(fp, "===============================================\n");
    fprintf(fp, "%-30s %15d\n", "use_lapack", use_lapack);
    fprintf(fp, "%-30s %15d\n", "use_femat", use_femat);
    fprintf(fp, "%-30s %15d\n", "lev", lev);
    fprintf(fp, "%-30s %15d\n", "fill", fill);
    fprintf(fp, "%-30s %15d\n", "reorder_schur", reorder_schur);    
    fprintf(fp, "%-30s %15d\n", "gmres_iters", gmres_iters);
    fprintf(fp, "%-30s %15d\n", "num_gmres_restarts", num_restarts);
    fprintf(fp, "%-30s %15d\n", "is_flexible", is_flexible);
    fprintf(fp, "===============================================\n");
  }
}

/*
  Print the adjoint options before marching backwards in time
*/
void TACSIntegrator::printAdjointOptionSummary( FILE *fp ) {
  if (fp && print_level>0){
    fprintf(fp, "===============================================\n");
    fprintf(fp, "Adjoint Mode : Parameter values\n");
    fprintf(fp, "===============================================\n");
    fprintf(fp, "%-30s %15d\n", "num_funcs", num_funcs);
    fprintf(fp, "%-30s %15d\n", "num_design_vars", num_design_vars);
    fprintf(fp, "===============================================\n");
  }
}

/*
  Control the amount of information printed to the console and the
  logging stream
*/
void TACSIntegrator::setPrintLevel( int _print_level, 
                                    const char *logfilename ){ 
  print_level = _print_level;
  if (logfilename){
    // Close any opened non stdout logstreams
    if (logfp && (logfp != stdout)){
      fclose(logfp);
    }

    // Open a new file for logstream
    if (mpiRank == 0){
      logfp = fopen(logfilename, "w");
      // Use stdout as output stream if the filename is empty
      if (!logfp){
        logfp = stdout;
      }
    }
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
void TACSIntegrator::setUseLapack( int _use_lapack, int _eigensolve ) {
  use_lapack = _use_lapack;
  eigensolve = _eigensolve;
  if (use_lapack == 1 || eigensolve == 1){
    // LAPACK works with natural ordering only
    this->setOrderingType(TACSAssembler::NATURAL_ORDER);
    // LAPACK needs FEMat
    this->setUseFEMat(1);
  }
}

/*
  Set FEMat usage flag into TACS
*/
void TACSIntegrator::setUseFEMat( int _use_femat ){
  use_femat = _use_femat;
}

/*
  Use line search for the solution of nonlinear problem
*/
void TACSIntegrator::setUseLineSearch( int _use_line_search ) {
  use_line_search = _use_line_search;
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
void TACSIntegrator::setTACSStates( double time, TACSBVec *q, 
                                    TACSBVec *qdot, TACSBVec * qddot ){
  tacs->setSimulationTime(time);
  tacs->setVariables(q, qdot, qddot);
}

/*
  Set external forces into the class variable
*/
void TACSIntegrator::setLoads( TACSBVec * forces ){
  if (this->forces){this->forces->decref();}
  if (forces){
    this->forces = forces;
    this->forces->incref();
  }
}

/*
  Update TACS states with the supplied ones (q, qdot, qddot)
*/
double TACSIntegrator::getTACSStates( TACSBVec *q, TACSBVec *qdot, TACSBVec * qddot ){
  tacs->getVariables(q, qdot, qddot);
  return tacs->getSimulationTime();
}

/*
  Set the ordering type within TACSIntegrator
*/
void TACSIntegrator::setOrderingType( TACSAssembler::OrderingType _type ){
  ordering_type = _type;
}

/*
  Set the initial fraction for the delta parameter within the Newton
  globalization strategy.
*/
void TACSIntegrator::setInitNewtonDeltaFraction( double frac ){
  if (frac >= 0.0){
    init_newton_delta = frac;
  }
  else {
    init_newton_delta = 0.0;
  }
}

/*
  Set if the current instance is due to adaptive time stepping
*/
void TACSIntegrator::setAdaptiveInstance( int flag ){
  adaptive_instance = flag;
}

/*
  Know if the current instance is due to adpative stepping
*/
int TACSIntegrator::isAdaptiveInstance(){
  return adaptive_instance;
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
                                      double _num_steps_per_sec, 
                                      int _max_bdf_order):
TACSIntegrator(_tacs, _tinit,  _tfinal,  _num_steps_per_sec){
  // copy over the variables
  max_bdf_order = _max_bdf_order;

  // Truncate the maximum order to 3rd order
  max_bdf_order = (max_bdf_order <= 3 ? 
                   max_bdf_order : 3);

  // Set the type of integrator
  if ( max_bdf_order == 3) {
    mytype = BDF3;
  } else if ( max_bdf_order == 2) {
    mytype = BDF2;
  } else if ( max_bdf_order == 1) {
    mytype = BDF1;
  }

  // Number of first and second order BDF coefficients
  nbdf = 0;
  nbddf = 0;
    
  // As many RHS as the number of second derivative coeffs
  num_adjoint_rhs = (2*max_bdf_order+1)+1; 
}

void TACSBDFIntegrator::initialize(){
  // Keep track of the time taken for foward mode
  time_forward = 0.0;
  time_fwd_assembly = 0.0;
  time_fwd_factor = 0.0;
  time_fwd_apply_factor = 0.0;
  time_newton = 0.0;
  double t0 = MPI_Wtime();
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
void TACSBDFIntegrator::getLinearizationCoeffs( double *alpha, 
                                                double *beta, 
                                                double *gamma ) {
  *gamma = bddf_coeff[0]/(h*h);
  *beta  = bdf_coeff[0]/h;
  *alpha = 1.0;
}

/*
  March backward in time and solve for the adjoint variables and find
  the total derivatives
*/
void TACSBDFIntegrator::marchBackwards( ) {
  // Print adjoint mode summary before maching backwards
  printAdjointOptionSummary(logfp);
    
  current_time_step = num_time_steps;

  time_rev_assembly     = 0.0;
  time_rev_factor       = 0.0;
  time_rev_apply_factor = 0.0;
  time_reverse          = 0.0;

  double t0 = MPI_Wtime();

  // Evaluate the function values
  this->evalFunctions(funcs, fvals, num_funcs);

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
  for ( int k = num_time_steps-1; k >= 1 ; k--){
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

    double tassembly = MPI_Wtime();

    // Setup the adjoint RHS
    for ( int n = 0; n < num_funcs; n++ ){
      // Add up the contribution from function state derivative to RHS
      psi[n]->zeroEntries();
      tacs->addSVSens(1.0, 0.0, 0.0, &funcs[n], 1, &psi[n]);

      // Add the contributions to the current adjoint RHS
      rhs[adj_index*num_funcs+n]->axpy(1.0, psi[n]);
      rhs[adj_index*num_funcs+n]->scale(-1.0);
    }
 
    // Setup the Jacobian
    tacs->assembleJacobian(alpha, beta, gamma, NULL, mat, TRANSPOSE);

    time_rev_assembly += MPI_Wtime() - tassembly;
     
    // LU factorization of the Jacobian
    double tfactor = MPI_Wtime();
    pc->factor();
    time_rev_factor += MPI_Wtime() - tfactor;

    // Apply the factorization for all right hand sides and solve for
    // the adjoint variables
    double tapply = MPI_Wtime();
    for ( int n = 0; n < num_funcs; n++ ){
      ksm->solve(rhs[adj_index*num_funcs + n], psi[n]);
      rhs[adj_index*num_funcs+n]->zeroEntries();
    }
    time_rev_apply_factor += MPI_Wtime() - tapply;

    // Add total derivative contributions from this step to all
    // functions
    double jacpdt = MPI_Wtime();
    addToTotalDerivative(h, psi);    
    time_rev_jac_pdt += MPI_Wtime() - jacpdt;

    // Drop the contributions from this step to other right hand sides
    double tassembly2 = MPI_Wtime();
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
    time_rev_assembly += MPI_Wtime() - tassembly2;
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

  // Print adjoint mode summary before maching backwards
  printAdjointOptionSummary(logfp);

  // Keep track of the time taken for foward mode
  time_reverse += MPI_Wtime() - t0;
}
