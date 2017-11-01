#include "TACSIntegrator.h"
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
TACSIntegrator::TACSIntegrator( TACSAssembler *_tacs,
                                double tinit,
                                double tfinal, 
                                double num_steps_per_sec ){
  // Copy over the input parameters
  tacs = _tacs; 
  tacs->incref();

  // Set the default prefix = results
  sprintf(prefix, "results");

  // Allocate the total number of time steps
  num_time_steps = int(num_steps_per_sec*(tfinal-tinit));

  // Store physical time of simulation
  time = new double[ num_time_steps+1 ];
  memset(time, 0, num_time_steps*sizeof(double));
  for ( int k = 0; k < num_time_steps+1; k++ ){
    time[k] = tinit + k*(tfinal - tinit)/num_time_steps;
  }

  // MPI information
  MPI_Comm_rank(tacs->getMPIComm(), &mpiRank);
  MPI_Comm_size(tacs->getMPIComm(), &mpiSize);
  
  // Default print level and logging control
  print_level = 2;
  logfp = NULL;
  if (mpiRank == 0){ 
    logfp = stdout; 
  }

  // State variables that store the entire time history
  q = new TACSBVec*[ num_time_steps+1 ];
  qdot = new TACSBVec*[ num_time_steps+1 ];
  qddot = new TACSBVec*[ num_time_steps+1 ];

  // create state vectors for TACS during each timestep
  for ( int k = 0; k < num_time_steps+1; k++ ) {
    q[k] = tacs->createVec(); q[k]->incref();
    qdot[k] = tacs->createVec(); qdot[k]->incref();
    qddot[k] = tacs->createVec(); qddot[k]->incref();
  }

  // Objects to store information about the functions of interest
  funcs = NULL;
  num_funcs = 0;
  fvals = NULL;
  dfdx = NULL;
  num_design_vars = 0;

  // Linear algebra objects for Newton's method
  res = tacs->createVec();
  update = tacs->createVec();
  res->incref();
  update->incref();

  // NULL the different KSM/solver objects
  mat = NULL;
  pc = NULL;
  ksm = NULL;

  // Default parameters for Newton's method
  max_newton_iters = 25;
  atol = 1.0e-6;
  rtol = 1.0e-9;
  init_newton_delta = 0.0;
  jac_comp_freq = 1;

  // Set the default LINEAR solver
  use_lapack = 0;
  use_femat = 1;
  order_type = TACSAssembler::TACS_AMD_ORDER;

  // AMD reordering parameters
  lev = 100000; 
  fill = 10.0;  
  reorder_schur = 1;

  // KSM parameters
  gmres_iters  = 10;
  num_restarts = 0;
  is_flexible  = 0;

  // Tecplot solution export
  f5_write_freq = 0; 

  // Set the rigid and shell visualization objects to NULL
  rigidf5 = NULL;
  shellf5 = NULL;
  beamf5 = NULL;
  
  // Set kinetic and potential energies
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
  for ( int k = 0; k < num_time_steps+1; k++ ) {
    q[k]->decref();
    qdot[k]->decref();
    qddot[k]->decref();
  }
  
  // Dereference Newton's method objects
  res->decref();
  update->decref();
  if (mat){ mat->decref(); }
  if (pc){ pc->decref(); }
  if (ksm){ ksm->decref(); }
  
  if (time){ delete [] time; }
  if (q){ delete [] q; }
  if (qdot){ delete [] qdot; }
  if (qddot){ delete [] qddot; }

  // Dereference TACS
  if (tacs){ tacs->decref(); }
  
  if (rigidf5){ rigidf5->decref();}
  if (shellf5){ shellf5->decref();}
  if (beamf5){ beamf5->decref();}
}

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
void TACSIntegrator::setUseLapack( int _use_lapack ){
  use_lapack = _use_lapack;
  if (use_lapack){
    // LAPACK works with natural ordering only
    setUseFEMat(1, TACSAssembler::NATURAL_ORDER);
  }
}

/*
  Set FEMat usage flag into TACS
*/
void TACSIntegrator::setUseFEMat( int _use_femat, 
                                  TACSAssembler::OrderingType _order_type ){
  use_femat = _use_femat;
  order_type = _order_type;
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
  Set the functions of interest that take part in the adjoint solve.
*/
void TACSIntegrator::setFunctions( TACSFunction **_funcs, int _num_funcs,
                                   int _num_design_vars,
                                   int _start_step, int _end_step ){
  // Set the time window
  start_step = _start_step;
  if (start_step < 0){
    start_step = 0;
  }
  end_step = _end_step;
  if (end_step < start_step){
    end_step = num_time_steps;
  }

  // Set the number of design variables
  num_design_vars = _num_design_vars;

  // Increase the reference counts to all functions to be set
  for ( int i = 0; i < _num_funcs; i++ ){
    _funcs[i]->incref();
  }

  // Free functions that have already been set
  if (funcs){
    for ( int i = 0; i < num_funcs; i++ ){
      funcs[i]->decref();
    }
    delete [] funcs;
    delete [] dfdx;
  }

  // Set the number of functions
  num_funcs = _num_funcs;
  funcs = new TACSFunction*[ num_funcs ];
  for ( int i = 0; i < num_funcs; i++ ){
    funcs[i] = _funcs[i];
  }
  dfdx = new TacsScalar[ num_funcs*num_design_vars ];
  memset(dfdx, 0, num_funcs*num_design_vars*sizeof(TacsScalar));
}

/*
  Set the output directory prefix
*/
void TACSIntegrator::setOutputPrefix( const char *_prefix ){
  strncpy(prefix, _prefix, sizeof(prefix));
}

/*
  Function that writes time, q, qdot, qddot to file
*/
void TACSIntegrator::writeSolution( const char *filename, int format ){
  FILE *fp = fopen(filename, "w");
  TacsScalar *qvals, *qdotvals, *qddotvals;

  if (format == 1){
    // Plain format with t q[0], q[1],...,qdot[0], qdot[1],...,qddot[0], qddot[1],...
    for ( int k = 0; k < num_time_steps; k++ ){    
      // Copy over the state values from TACSBVec
      int num_state_vars = q[k]->getArray(&qvals);
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
  } 
  else {
    if (format > 0){
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
        int num_state_vars = q[k]->getArray(&qvals);
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
    } 
    else {
      // Write the DOFS on user specified element number in final ordering  
      // Plain format with t q[0], q[1],...,qdot[0], qdot[1],...,qddot[0], qddot[1],...
      for ( int k = 0; k < num_time_steps; k++ ){    
        // Copy over the state values from TACSBVec
        int num_state_vars = q[k]->getArray(&qvals);
        qdot[k]->getArray(&qdotvals);
        qddot[k]->getArray(&qddotvals);
  
        // Write the time and states to file
        fprintf(fp, "%12.5e ", time[k]);  
        int elem_ctr = 0;
        for ( int j = 0; j < num_state_vars; j++ ){
          // Write if we have found the sought element
          if (elem_ctr == -format){
            fprintf(fp, " %12.5e %12.5e %12.5e ",
                    TacsRealPart(qvals[j]), 
                    TacsRealPart(qdotvals[j]), 
                    TacsRealPart(qddotvals[j]));
          }
          if (j % 8 == 7){
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
  Creates an f5 file for each time step and writes the data.
*/
void TACSIntegrator::writeSolutionToF5( int step_num ){
  if (step_num >= 0 && step_num < num_time_steps){
    // Loop through all timesteps
    for ( int k = 0; k < num_time_steps; k++ ){
      // Set the current states into TACS
      tacs->setVariables(q[k], qdot[k], qddot[k]);
      tacs->setSimulationTime(time[k]);
      
      // Write RIGID body if set
      if (rigidf5){
        char fname[256];
        sprintf(fname, "%s/rigid_%06d.f5", prefix, step_num);
        rigidf5->writeToFile(fname);
      }

      // Write SHELL body if set
      if (shellf5){
        char fname[256];
        sprintf(fname, "%s/rigid_%06d.f5", prefix, step_num);
        shellf5->writeToFile(fname);
      }

      // Write BEAM body if set
      if (beamf5){
        char fname[256];
        sprintf(fname, "%s/rigid_%06d.f5", prefix, step_num);
        beamf5->writeToFile(fname);
      }
    }
  }
}

/*
  Configure the F5 output 
*/
void TACSIntegrator::setOutputFrequency( int _write_freq ){
  f5_write_freq = _write_freq;
}

/*
  Set whether RIGID body components are a part of the output
*/
void TACSIntegrator::setRigidOutput( TACSToFH5 *_rigidf5 ){
  if (_rigidf5){
    _rigidf5->incref();
  }
  if (rigidf5){
    rigidf5->decref();
  }
  rigidf5 = _rigidf5;
}

/*
  Set whether SHELL body components are a part of the output
*/
void TACSIntegrator::setShellOutput( TACSToFH5 *_shellf5 ){
  if (_shellf5){
    _shellf5->incref();
  }
  if (shellf5){
    shellf5->decref();
  }
  shellf5 = _shellf5;
}

/*
  Set whether SHELL body components are a part of the output
*/
void TACSIntegrator::setBeamOutput( TACSToFH5 *_beamf5 ){
  if (_beamf5){
    _beamf5->incref();
  }
  if (beamf5){
    beamf5->decref();
  }
  beamf5 = _beamf5;
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
  Write out all of the options that have been set to a output stream.
*/
void TACSIntegrator::printOptionSummary(){
  if (logfp && print_level>0){
    fprintf(logfp, "===============================================\n");
    fprintf(logfp, "TACSIntegrator: Parameter values\n");
    fprintf(logfp, "===============================================\n");
    fprintf(logfp, "%-30s %15g\n", "tinit", time[0] );
    fprintf(logfp, "%-30s %15g\n", "tfinal", time[num_time_steps-1]);
    fprintf(logfp, "%-30s %15d\n", "num_time_steps", num_time_steps );
    fprintf(logfp, "%-30s %15d\n", "mpiSize", mpiSize);
    fprintf(logfp, "%-30s %15d\n", "print_level", print_level);

    fprintf(logfp, "===============================================\n");
    fprintf(logfp, "Nonlinear Solver: Parameter values\n");
    fprintf(logfp, "===============================================\n");    
    fprintf(logfp, "%-30s %15d\n", "max_newton_iters", max_newton_iters);
    fprintf(logfp, "%-30s %15g\n", "absolute_tolerance", atol);
    fprintf(logfp, "%-30s %15g\n", "relative_tolerance", rtol);
    fprintf(logfp, "%-30s %15d\n", "jac_comp_freq", jac_comp_freq);

    fprintf(logfp, "===============================================\n");
    fprintf(logfp, "Linear Solver: Parameter values\n");
    fprintf(logfp, "===============================================\n");
    fprintf(logfp, "%-30s %15d\n", "use_lapack", use_lapack);
    fprintf(logfp, "%-30s %15d\n", "use_femat", use_femat);
    fprintf(logfp, "%-30s %15d\n", "lev", lev);
    fprintf(logfp, "%-30s %15g\n", "fill", fill);
    fprintf(logfp, "%-30s %15d\n", "reorder_schur", reorder_schur);    
    fprintf(logfp, "%-30s %15d\n", "gmres_iters", gmres_iters);
    fprintf(logfp, "%-30s %15d\n", "num_gmres_restarts", num_restarts);
    fprintf(logfp, "%-30s %15d\n", "is_flexible", is_flexible);
    fprintf(logfp, "===============================================\n");
  }
}

/*
  Print the adjoint options before marching backwards in time
*/
void TACSIntegrator::printAdjointOptionSummary(){
  if (logfp && print_level>0){
    fprintf(logfp, "===============================================\n");
    fprintf(logfp, "Adjoint Mode : Parameter values\n");
    fprintf(logfp, "===============================================\n");
    fprintf(logfp, "%-30s %15d\n", "num_funcs", num_funcs);
    fprintf(logfp, "%-30s %15d\n", "num_design_vars", num_design_vars);
    fprintf(logfp, "===============================================\n");
  }
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

  Returns: Integer: Termination of nonlinear solver 
  1: |R| < atol; 
  2: |dq| < atol
  3: |R|/|R0| < rtol
  -1: max_newton_iters 
  -2: Nan
*/
int TACSIntegrator::newtonSolve( double alpha, double beta, double gamma,
                                 double t, TACSBVec *u, TACSBVec *udot, 
                                 TACSBVec *uddot, TACSBVec *forces ){
  if (!mat || !ksm){
    // Set the D matrix to NULL
    if (use_femat){
      // Create a matrix for storing the Jacobian
      FEMat *D = tacs->createFEMat(order_type);
      
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
    fprintf(logfp, "%12s %12s %12s %12s %12s %12s %12s %12s\n",
            "#iters", "|R|", "|R|/|R0|", "|dq|", 
            "alpha", "beta", "gamma","delta");
  }

  // Track the time taken for newton solve at each time step
  double tnewton = MPI_Wtime();

  // Iterate until max iters or R <= tol
  double delta = 0.0;
  for ( niter = 0; niter < max_newton_iters; niter++ ){    
    // Set the supplied initial input states into TACS
    tacs->setSimulationTime(t);
    tacs->setVariables(u, udot, uddot);

    // Assemble the Jacobian matrix once in Newton iterations
    double t0 = MPI_Wtime();
    if ((niter % jac_comp_freq) == 0){
      delta = init_newton_delta*gamma;
      if (niter > 0 && (TacsRealPart(res_norm) < TacsRealPart(init_res_norm))){
        delta *= TacsRealPart(res_norm/init_res_norm);
      }

      tacs->assembleJacobian(alpha, beta, gamma + delta,
                             res, mat, NORMAL);
    }
    else {
      tacs->assembleRes(res);
    }

    // Add the forces into the residual    
    if (forces){
      tacs->applyBCs(forces);
      res->axpy(-1.0, forces);
      tacs->applyBCs(res);
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
        fprintf(logfp, "%12d %12.5e %12.5e %12s %12.5e %12.5e %12.5e %12.5e\n",
                niter, TacsRealPart(res_norm),  
                (niter == 0) ? 1.0 : TacsRealPart(res_norm/init_res_norm), 
                " ", alpha, beta, gamma, delta);
      }
      else {
        fprintf(logfp, "%12d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
                niter, TacsRealPart(res_norm),  
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
  int num_state_vars = update->getArray(&ans);
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
  Implement all the tasks to perform during each time step
*/
void TACSIntegrator::logTimeStep( int step_num ){
  if (step_num == 0){
    // Keep track of the time taken for foward mode
    time_forward = MPI_Wtime();
    time_fwd_assembly = 0.0;
    time_fwd_factor = 0.0;
    time_fwd_apply_factor = 0.0;
    time_newton = 0.0;
  }
  if (step_num == num_time_steps-1){
    time_forward = MPI_Wtime() - time_forward;
  }

  // Write the tecplot output to disk if sought
  if (f5_write_freq > 0 && step_num % f5_write_freq == 0){
    writeSolutionToF5(step_num);
  }

  // Evaluate the energies
  TacsScalar energies[2];
  tacs->evalEnergies(&energies[0], &energies[1]);

  if (step_num == 0){
    // Log information
    if (logfp && print_level >= 1){
      fprintf(logfp, "%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n", 
              "status", "time", "tnewton", "#iters",
              "|R|", "|R|/|R0|", "|dq|", "KE", "PE", "E0-E");
      
      // Compute the initial energy
      init_energy = energies[0] + energies[1];

      // Log the details
      fprintf(logfp, "%6d/%-6d %12.5e %12.5e %12d %12.5e \
%12.5e %12.5e %12.5e %12.5e %12.5e\n",
              0, num_time_steps,
              time[0], time_newton, 0, 0.0, 0.0, 0.0,
              TacsRealPart(energies[0]), TacsRealPart(energies[1]),  
              0.0);
    }
  } 
  else {
    // Print out the time step summary
    if (logfp && print_level >= 1){
      // Need a title for total summary as details of Newton iteration
      // will overshadow this one line summary
      if (print_level == 2){
        fprintf(logfp, "%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
                "status",
                "time", "tnewton", "#iters",
                "|R|", "|R|/|R0|", "|dq|", "KE", "PE", "E0-E");
      }
      
      fprintf(logfp, "%6d/%-6d %12.5e %12.5e %12d %12.5e \
%12.5e %12.5e %12.5e %12.5e %12.5e\n",
              step_num, num_time_steps,
              time[step_num], time_newton, niter,
              TacsRealPart(res_norm), 
              TacsRealPart(res_norm/(rtol + init_res_norm)),
              TacsRealPart(update_norm),
              TacsRealPart(energies[0]), TacsRealPart(energies[1]), 
              TacsRealPart((init_energy - (energies[0] + energies[1]))));
    }
  }
}

/*
  Get the state variables from TACS at the given step
*/
double TACSIntegrator::getStates( int step_num, 
                                  TACSBVec *_q, TACSBVec *_qdot, TACSBVec *_qddot){
  if (_q){
    _q->copyValues(q[step_num]);
  }

  if (_qdot){
    _qdot->copyValues(qdot[step_num]);
  }

  if (_qddot){
    _qddot->copyValues(qddot[step_num]);
  }
  return time[step_num];
}

/*
  Function that returns the finite-difference/complex-step gradient
  (used for testing purposes)
*/
void TACSIntegrator::checkGradients( double dh ){
 // Check whether the function has been set properly
  if (num_funcs == 0 || funcs == NULL) {
    fprintf(stderr, "TACS Warning: Function is not set, skipping adjoint solve. \n");
    return;
  }

  // Get the design variable values
  TacsScalar *x = new TacsScalar[ num_design_vars ];
  tacs->getDesignVars(x, num_design_vars);

  // Create a temporary vector of design variable values
  TacsScalar *xtmp = new TacsScalar[ num_design_vars ];
  memcpy(xtmp, x, num_design_vars*sizeof(TacsScalar)); 
 
  // Allocate an array of function values
  TacsScalar *ftmp1 = new TacsScalar[ num_funcs ];
  TacsScalar *ftmp2 = new TacsScalar[ num_funcs ];
  
  // Form a finite-difference or complex step approximation of the
  // total derivative
#ifdef TACS_USE_COMPLEX
  for ( int k = 0; k < num_design_vars; k++ ){
    xtmp[k] = x[k] + TacsScalar(0.0, dh);
  }
  // Set the design variables
  tacs->setDesignVars(xtmp, num_design_vars);

  // Integrate forward in time
  integrate();

  // Evaluate the functions
  evalFunctions(ftmp2);

  // Compute the total derivative
  for ( int k = 0; k < num_funcs; k++ ){
    ftmp2[k] = TacsImagPart(ftmp2[k])/dh;
  }

  // Evaluate the derivative
  integrateAdjoint();
#else // REAL code
  // Integrate forward in time
  integrate();

  // Evaluate the functions
  evalFunctions(ftmp1);

  // Evaluate the derivative
  integrateAdjoint();

  for ( int k = 0; k < num_design_vars; k++ ){
    xtmp[k] = x[k] + dh;
  }

  // Set the design variables
  tacs->setDesignVars(xtmp, num_design_vars);

  // Integrate forward in time
  integrate();

  // Evaluate the functions
  evalFunctions(ftmp2);

  // Compute the total derivative
  for ( int k = 0; k < num_funcs; k++ ){
    ftmp2[k] = (ftmp2[k] - ftmp1[k])/dh;
  }
#endif // TACS_USE_COMPLEX

  // Compute the total projected derivative
  for ( int i = 0; i < num_funcs; i++ ){
    ftmp1[i] = 0.0;
    for ( int j = 0; j < num_design_vars; j++ ){
      ftmp1[i] += dfdx[num_design_vars*i + j];
    }
  }

  printf("%3s %15s %15s %15s\n", "Fn", "Adjoint", "FD", "Relative Err");
  for ( int i = 0; i < num_funcs; i++ ){
    double relerr = fabs(TacsRealPart((ftmp1[i] - ftmp2[i])/ftmp2[i]));
    printf("%3d %15.8e %15.8e %15.8e\n",
           i, TacsRealPart(ftmp1[i]), TacsRealPart(ftmp2[i]), relerr);
  }

  delete [] x;
  delete [] ftmp1;
  delete [] ftmp2;
  delete [] xtmp;
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

  // Set the adjoint variables and right-hand-sides to NULL
  rhs = NULL;
  psi = NULL;
    
  // As many RHS as the number of second derivative coeffs
  num_adjoint_rhs = (2*max_bdf_order+1)+1; 
}

/*
  Free any data that is allocated by this class
*/
TACSBDFIntegrator::~TACSBDFIntegrator(){
  if (rhs){
    for ( int i = 0; i < num_funcs; i++ ){
      psi[i]->decref();
    } 
    for ( int i = 0; i < num_funcs*num_adjoint_rhs; i++ ){
      rhs[i]->decref();
    }
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
  March one step and exit the integration. This is primarily for use 
  with FUNtoFEM. 
*/
int TACSBDFIntegrator::iterate( int k, TACSBVec *forces ){
  if (k == 0){    
    // Retrieve the initial conditions and set into TACS
    tacs->getInitConditions(q[0], qdot[0], qddot[0]);
    tacs->setVariables(q[0], qdot[0], qddot[0]);

    // Output the results at the initial condition if configured
    printOptionSummary();
    logTimeStep(k);
    return 0;
  }
 
  // Advance time
  double h = time[k] - time[k-1];
  
  // Extrapolate to next time step: q[k] = q[k-1] + h*qdot[k-1] +
  // h^2/2*qddot[k-1] (helps reduce the initial residual)
  q[k]->copyValues(q[k-1]);
  q[k]->axpy(h, qdot[k-1]);
  q[k]->axpy(0.5*h*h, qddot[k-1]);

  // Get the BDF coefficients
  int nbdf, nbddf;
  double bdf_coeff[4];
  double bddf_coeff[9];
  get2ndBDFCoeff(k, bdf_coeff, &nbdf, 
                 bddf_coeff, &nbddf, max_bdf_order);
  
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
  double alpha = 1.0;
  double beta = bdf_coeff[0]/h;
  double gamma = bddf_coeff[0]/(h*h);
  
  // Solve the nonlinear system of stage equations starting with the
  // approximated states
  int newton_term = newtonSolve(alpha, beta, gamma,
                                time[k], q[k], qdot[k], qddot[k],
                                forces);

  // Tecplot output and print related stuff as configured
  logTimeStep(k);  

  // Return a non-zero flag when the Newton iteration fails
  int fail = 0;
  if (newton_term < 0){
    fail = 1;
  }
  return fail;
}

/*
  Evaluate the functions of interest
*/
void TACSBDFIntegrator::evalFunctions( TacsScalar *fvals ){
  // Check whether these are two-stage or single-stage functions
  int twoStage = 0;
  for ( int n = 0; n < num_funcs; n++ ){
    if (funcs[n]->getStageType() == TACSFunction::TWO_STAGE){
      twoStage = 1;
      break;
    }
  }

  // Initialize the function if had already not been initialized
  if (twoStage){
    // First stage
    for ( int n = 0; n < num_funcs; n++ ){
      funcs[n]->initEvaluation(TACSFunction::INITIALIZE);
    }
    
    for ( int k = 0; k < num_time_steps+1; k++ ){
      // Set the stages
      tacs->setSimulationTime(time[k]);
      tacs->setVariables(q[k], qdot[k], qddot[k]);

      double tcoeff = 0.0;
      if (k > 0){
        tcoeff += 0.5*(time[k] - time[k-1]);
      }
      if (k < num_time_steps){
        tcoeff += 0.5*(time[k+1] - time[k]);
      }
      tacs->integrateFunctions(tcoeff, TACSFunction::INITIALIZE,
                               funcs, num_funcs);
    }

    for ( int n = 0; n < num_funcs; n++ ){
      funcs[n]->finalEvaluation(TACSFunction::INITIALIZE);      
    }
  }

  // Second stage
  for ( int n = 0; n < num_funcs; n++ ){
    funcs[n]->initEvaluation(TACSFunction::INTEGRATE);
  }
    
  for ( int k = 0; k < num_time_steps+1; k++ ){
    tacs->setSimulationTime(time[k]);
    tacs->setVariables(q[k], qdot[k], qddot[k]);
      
    double tcoeff = 0.0;
    if (k > 0){
      tcoeff += 0.5*(time[k] - time[k-1]);
    }
    if (k < num_time_steps){
      tcoeff += 0.5*(time[k+1] - time[k]);
    }
    tacs->integrateFunctions(tcoeff, TACSFunction::INTEGRATE,
                             funcs, num_funcs);
  }
       
  for ( int n = 0; n < num_funcs; n++ ){
    funcs[n]->finalEvaluation(TACSFunction::INTEGRATE);    
  }

  // Retrieve the function values
  for ( int n = 0; n < num_funcs; n++ ){
    fvals[n] = funcs[n]->getFunctionValue();
  }
}

/*
  Initialize the right-hand-side contributions to the adjoint, compute
  the Jacobian and factorize it.
*/
void TACSBDFIntegrator::initAdjoint( int k ){
  // Adjoint variables for each function of interest
  if (!psi){
    psi = new TACSBVec*[ num_funcs ];
    for ( int i = 0; i < num_funcs; i++ ){
      psi[i] = tacs->createVec();
      psi[i]->incref();
    }
  }

  if (!rhs){
    // Right hand sides for adjoint linear-system
    rhs = new TACSBVec*[ num_funcs*num_adjoint_rhs ];
    for ( int i = 0; i < num_funcs*num_adjoint_rhs; i++ ){
      rhs[i] = tacs->createVec();
      rhs[i]->incref();
    }
  }

  if (k == num_time_steps){
    // Print adjoint mode summary before maching backwards
    printAdjointOptionSummary();
    
    double t0 = MPI_Wtime();
    time_rev_assembly = 0.0;
    time_rev_factor = 0.0;
    time_rev_apply_factor = 0.0;
    time_reverse = t0;
    
    // Zero the right-hand-sides and adjoint
    for ( int i = 0; i < num_funcs; i++ ){
      psi[i]->zeroEntries();
    }
    for ( int i = 0; i < num_funcs*num_adjoint_rhs; i++ ){
      rhs[i]->zeroEntries();
    }
  } 

  // Set the simulation time
  tacs->setSimulationTime(time[k]);
  tacs->setVariables(q[k], qdot[k], qddot[k]);

  if (k > 0){
    // Get the BDF coefficients at this time step
    int nbdf, nbddf;
    double bdf_coeff[4];
    double bddf_coeff[9];
    get2ndBDFCoeff(k, bdf_coeff, &nbdf, bddf_coeff, 
                   &nbddf, max_bdf_order);

    // Keep track of the multiplier on the function value (divided by
    // the time step)
    double tfact = 0.5;
    if (k < num_time_steps){
      tfact += 0.5;
    }

    // Compute the time interval
    double h = time[k] - time[k-1];

    // Determine the linearization coefficients for Jacobian Assembly
    double alpha = 1.0;
    double beta  = bdf_coeff[0]/h;
    double gamma = bddf_coeff[0]/(h*h);
    
    // Find the adjoint index
    int adj_index = k % num_adjoint_rhs;
    
    // Setup the adjoint RHS
    for ( int n = 0; n < num_funcs; n++ ){
      // Add up the contribution from function state derivative to RHS
      psi[n]->zeroEntries();
      tacs->addSVSens(tfact, 0.0, 0.0, &funcs[n], 1, &psi[n]);
      
      // Add the contributions to the current adjoint RHS
      rhs[adj_index*num_funcs + n]->axpy(1.0, psi[n]);
      rhs[adj_index*num_funcs + n]->scale(-1.0);
    }
 
    // Setup the Jacobian
    double tassembly = MPI_Wtime();
    tacs->assembleJacobian(alpha, beta, gamma, NULL, mat, TRANSPOSE);
    time_rev_assembly += MPI_Wtime() - tassembly;
    
    // LU factorization of the Jacobian
    double tfactor = MPI_Wtime();
    pc->factor();
    time_rev_factor += MPI_Wtime() - tfactor;
  }
}

/*
  Solve the adjoint equations, including any extra external
  contributions to the right-hand-side
*/
void TACSBDFIntegrator::iterateAdjoint( int k, TACSBVec **adj_rhs ){
  // Find the adjoint index
  int adj_index = k % num_adjoint_rhs;

  // Apply the factorization for all right hand sides and solve for
  // the adjoint variables
  double tapply = MPI_Wtime();
  for ( int n = 0; n < num_funcs; n++ ){
    if (adj_rhs){
      // Use the residual vector as a temp vector for the purposes
      // of this computation
      res->copyValues(adj_rhs[n]);
      res->axpy(1.0, rhs[adj_index*num_funcs + n]);
      ksm->solve(res, psi[n]);
    }
    else {
      ksm->solve(rhs[adj_index*num_funcs + n], psi[n]);
    }
  }
  time_rev_apply_factor += MPI_Wtime() - tapply;
}

/*
  Add the additional terms for the right-hand-side to the
  right-hand-side
*/
void TACSBDFIntegrator::postAdjoint( int k ){
  // Find the adjoint index
  int adj_index = k % num_adjoint_rhs;

  // Zero the right-hand-sides that we just solved
  for ( int n = 0; n < num_funcs; n++ ){
    rhs[adj_index*num_funcs + n]->zeroEntries();
  }

  double tcoeff = 0.0;
  if (k > 0){
    tcoeff += 0.5*(time[k] - time[k-1]);
  }
  if (k < num_time_steps){
    tcoeff += 0.5*(time[k+1] - time[k]);
  }
  tacs->addDVSens(tcoeff, funcs, num_funcs, dfdx, num_design_vars);

  if (k > 0){
    // Get the BDF coefficients at this time step
    int nbdf, nbddf;
    double bdf_coeff[4];
    double bddf_coeff[9];
    get2ndBDFCoeff(k, bdf_coeff, &nbdf, bddf_coeff, 
                   &nbddf, max_bdf_order);

    // Compute the time interval
    double h = time[k] - time[k-1];
    
    // Determine the linearization coefficients for Jacobian Assembly
    double alpha = 1.0;
    double beta  = bdf_coeff[0]/h;
    double gamma = bddf_coeff[0]/(h*h);
    
    // Add total derivative contributions from this step to all
    // functions
    double jacpdt = MPI_Wtime();
    tacs->addAdjointResProducts(h, psi, num_funcs, dfdx, num_design_vars);
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
                                    psi[n], rhs[rhs_index*num_funcs+n], 
                                    TRANSPOSE);
      }
    }
    time_rev_assembly += MPI_Wtime() - tassembly2;
  }

  if (k == 0){
    // Keep track of the time taken for foward mode
    time_reverse = MPI_Wtime() - time_reverse;
  }  
}

/*
  In this case, only the adjoint at the current iteration is stored.
*/
void TACSBDFIntegrator::getAdjoint( int step_num, int func_num, 
                                    TACSBVec **adjoint ){
  *adjoint = psi[func_num];
}
