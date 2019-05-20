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

#include "TACSIntegrator.h"
#include <math.h>
#include "tacslapack.h"
#include "TACSMg.h"

/*
  Base class constructor for integration schemes.

  input:
  tacs:               the tacs assembler object
  tInit:              the initial time
  tfinal:             the final time
  num_steps_per_sec:  the number of steps to take for each second
*/
TACSIntegrator::TACSIntegrator( TACSAssembler *_tacs,
                                double tinit,
                                double tfinal,
                                double num_steps ){
  // Copy over the input parameters
  tacs = _tacs;
  tacs->incref();

  // Set the default prefix = results
  sprintf(prefix, "results");

  // Allocate the total number of time steps
  num_time_steps = int(num_steps);

  // Store physical time of simulation
  time = new double[ num_time_steps+1 ];
  memset(time, 0, num_time_steps*sizeof(double));
  for ( int k = 0; k < num_time_steps+1; k++ ){
    time[k] = tinit + double(k)*(tfinal - tinit)/double(num_time_steps);
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
  start_plane = 0;
  end_plane = num_time_steps;
  num_funcs = 0;
  fvals = NULL;
  dfdx = NULL;
  num_design_vars = 0;
  dfdXpt = NULL;

  // Linear algebra objects for Newton's method
  res = tacs->createVec();
  update = tacs->createVec();
  res->incref();
  update->incref();

  // NULL the different KSM/solver objects
  mat = NULL;
  pc = NULL;
  ksm = NULL;
  linear_solver_initialized = 0;

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
  solidf5 = NULL;

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
  if (solidf5){ solidf5->decref();}
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
  Set the liner solver
 */
void TACSIntegrator::setKrylovSubspaceMethod( TACSKsm *_ksm ){
  if (_ksm){
    _ksm->incref();
  }
  if (ksm){
    ksm->decref();
  }
  ksm = _ksm;
}

/*
  Set the time interval for the simulation
*/
void TACSIntegrator::setTimeInterval( double tinit, double tfinal ){
  for ( int k = 0; k < num_time_steps+1; k++ ){
    time[k] = tinit + double(k)*(tfinal - tinit)/double(num_time_steps);
  }
}

/*
  Set the functions of interest that take part in the adjoint solve.
*/
void TACSIntegrator::setFunctions( TACSFunction **_funcs, int _num_funcs,
                                   int _num_design_vars,
                                   int _start_plane, int _end_plane ){
  // Set the time window
  start_plane = _start_plane;
  if (start_plane < 0){
    start_plane = 0;
  }
  end_plane = _end_plane;
  if (end_plane <= start_plane){
    end_plane = num_time_steps;
  }
  if (end_plane > num_time_steps){
    end_plane = num_time_steps;
  }

  // Set the number of design variables
  num_design_vars = _num_design_vars;

  // Increase the reference counts to all functions to be set
  for ( int i = 0; i < _num_funcs; i++ ){
    if (_funcs[i]){
      _funcs[i]->incref();
    }
  }

  // Free functions that have already been set
  if (funcs){
    for ( int i = 0; i < num_funcs; i++ ){
      if (funcs[i]){
        funcs[i]->decref();
      }
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

  dfdXpt = new TACSBVec*[ num_funcs ];
  for ( int i = 0; i < num_funcs; i++ ){
    dfdXpt[i] = tacs->createNodeVec();
    dfdXpt[i]->incref();
  }
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
void TACSIntegrator::writeRawSolution( const char *filename, int format ){
  FILE *fp = fopen(filename, "w");
  TacsScalar *qvals, *qdotvals, *qddotvals;

  if (format == 1){
    for ( int k = 0; k < num_time_steps + 1; k++ ){
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
      for ( int k = 0; k < num_time_steps + 1; k++ ){
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
      for ( int k = 0; k < num_time_steps + 1; k++ ){
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
  Creates f5 files for each time step.
*/
void TACSIntegrator::writeSolutionToF5(){
  for ( int k = 0; k < num_time_steps + 1; k++ ){
    writeStepToF5(k);
  }
}

/*
  Creates an f5 file for each time step and writes the data.
*/
void TACSIntegrator::writeStepToF5( int step_num ){
  // Set the current states into TACS
  tacs->setVariables(q[step_num], qdot[step_num], qddot[step_num]);
  tacs->setSimulationTime(time[step_num]);

  // Write RIGID body if set
  if (rigidf5){
    char fname[256];
    sprintf(fname, "%s/rigid_%06d.f5", prefix, step_num);
    rigidf5->writeToFile(fname);
  }

  // Write SHELL body if set
  if (shellf5){
    char fname[256];
    sprintf(fname, "%s/shell_%06d.f5", prefix, step_num);
    shellf5->writeToFile(fname);
  }

  // Write BEAM body if set
  if (beamf5){
    char fname[256];
    sprintf(fname, "%s/beam_%06d.f5", prefix, step_num);
    beamf5->writeToFile(fname);
  }

  // Write solid body if set
  if (solidf5){
    char fname[256];
    sprintf(fname, "%s/solid_%06d.f5", prefix, step_num);
    solidf5->writeToFile(fname);
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
  Set whether SOLID body components are a part of the output
*/
void TACSIntegrator::setSolidOutput( TACSToFH5 *_solidf5 ){
  if (_solidf5){
    _solidf5->incref();
  }
  if (solidf5){
    solidf5->decref();
  }
  solidf5 = _solidf5;
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
    fprintf(logfp, "%-30s %15g\n", "tfinal", time[num_time_steps]);
    fprintf(logfp, "%-30s %15g\n", "step_size", time[1]-time[0] );
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
  The number of time steps configured during instantiation
*/
int TACSIntegrator::getNumTimeSteps(){
  return num_time_steps;
}

/*
  Perform a nonlinear solve to obtain accelerations with q and qdot
  held constant as obtained from intial conditions.
*/
int TACSIntegrator::initAccelerationSolve(TACSBVec *forces){
 double force_norm = 0.0;
  if(forces){
    force_norm = TacsRealPart(forces->norm());
  }

  // initialize the linear solver
  initializeLinearSolver();

  // Create KSM
  TACSMg *mg = dynamic_cast<TACSMg*>(pc);
  if (mg){
    mg->incref();
  };

  if (logfp && print_level >= 2){
    fprintf(logfp, "%12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
            "#iters", "|R|", "|R|/|R0|", "|dq|",
            "alpha", "beta", "gamma","delta", "|F|");
  }

  // Assemble the residual and jacobian
  double alpha = 0.0;
  double beta  = 0.0;
  double gamma = 1.0;
    
  // Iterate accelerations until residual is zero
  for ( niter = 0; niter < max_newton_iters; niter++ ){

    // Assemble residual and jacobian
    if (mg){
      mg->assembleJacobian(alpha, beta, gamma, res, NORMAL);
    }
    else {
      tacs->assembleJacobian(alpha, beta, gamma, res, mat, NORMAL);
    }

    // Add the forces into the residual
    if (forces){
      tacs->applyBCs(forces);
      res->axpy(-1.0, forces);
      tacs->applyBCs(res);
    }   
    res_norm = res->norm();    

    // Record the residual norm at the first Newton iteration
    if (niter == 0){
      init_res_norm = res_norm;
    }

    // Write a summary
    if(logfp && print_level >= 2){
      if (niter == 0){
        fprintf(logfp,
                "%12d %12.5e %12.5e %12s %12.5e %12.5e %12.5e %12.5e %12.5e\n",
                niter,
                TacsRealPart(res_norm),
                (niter == 0) ? 1.0 : TacsRealPart(res_norm/init_res_norm),
                " ",
                alpha, beta, gamma, 0.0, force_norm);
      }
      else {
        fprintf(logfp,
                "%12d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
                niter,
                TacsRealPart(res_norm),
                (niter == 0) ? 1.0 : TacsRealPart(res_norm/init_res_norm),
                TacsRealPart(update_norm),
                alpha, beta, gamma, 0.0, force_norm);
      }
    }

    // Check if the Newton convergence tolerance is satisfied
    if (TacsRealPart(res_norm) < atol){
      break;
    }

    // LU Factor the matrix when needed
    if ((niter % jac_comp_freq) == 0){
      pc->factor();
    }

    // Solve for update using KSM
    ksm->solve(res, update);
     
    // Update the state variables using the solution
    qddot[0]->axpy(-1.0, update);

    // Set states into TACS (new acceleration)
    tacs->setVariables(q[0], qdot[0], qddot[0]);

    // Check whether the update is sufficiently small
    if (TacsRealPart(update_norm) < atol){
      break;
    }

  }

  if (mg){
    mg->decref();
  };
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
  // Compute the norm of the forces if supplied
  double force_norm = 0.0;
  if(forces){
    force_norm = TacsRealPart(forces->norm());
  }

  // Create KSM
  initializeLinearSolver();

  // Initialize the update norms
  update_norm = 1.0e99;

  // Initialize the residual norms
  init_res_norm = 0.0;
  res_norm = 0.0;
  int newton_exit_flag = 0;

  if (logfp && print_level >= 2){
    fprintf(logfp, "%12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
            "#iters", "|R|", "|R|/|R0|", "|dq|",
            "alpha", "beta", "gamma","delta", "|F|");
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
      if (niter > 0 &&
          (TacsRealPart(res_norm) < TacsRealPart(init_res_norm))){
        delta *= TacsRealPart(res_norm/init_res_norm);
      }

      // Try to down cast the preconditioner to a multigrid pc
      TACSMg *mg = dynamic_cast<TACSMg*>(pc);
      if (mg){
        mg->assembleJacobian(alpha, beta, gamma + delta,
                             res, NORMAL);
      }
      else {
        tacs->assembleJacobian(alpha, beta, gamma + delta,
                               res, mat, NORMAL);
      }
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
        fprintf(logfp,
                "%12d %12.5e %12.5e %12s %12.5e %12.5e %12.5e %12.5e %12.5e\n",
                niter,
                TacsRealPart(res_norm),
                (niter == 0) ? 1.0 : TacsRealPart(res_norm/init_res_norm),
                " ",
                alpha, beta, gamma, delta, force_norm);
      }
      else {
        fprintf(logfp,
                "%12d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
                niter,
                TacsRealPart(res_norm),
                (niter == 0) ? 1.0 : TacsRealPart(res_norm/init_res_norm),
                TacsRealPart(update_norm),
                alpha, beta, gamma, delta, force_norm);
      }
    }

    // Check if the norm of the residuals is a NaN
    if (res_norm != res_norm || update_norm != update_norm ){
      if (logfp) {
        fprintf(stderr,
                "[%d] Newton iteration %d, failed with NaN residual norm\n",
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
  Obtain a column-major array containing the matrix values.
*/
void TACSIntegrator::getRawMatrix( TACSMat *mat, TacsScalar *mat_vals) {
  FEMat *femat = dynamic_cast<FEMat*>(mat);
  if (femat){
    int bsize, nrows;
    BCSRMat *B;
    femat->getBCSRMat(&B, NULL, NULL, NULL);
    B->getArrays(&bsize, &nrows, NULL, NULL, NULL, NULL);
    B->getDenseColumnMajor(mat_vals);
  }
  else {
    fprintf(stderr,"Casting to FEMat failed\n");
  }
}

/*
  Performs an eigen solve of the Jacobian matrix
*/
int TACSIntegrator::lapackNaturalFrequencies( int use_gyroscopic,
                                              TACSBVec *q,
                                              TACSBVec *qdot,
                                              TACSBVec *qddot,
                                              TacsScalar *freq,
                                              TacsScalar *modes ){
  // TACSVec for mode
  TACSBVec *mode = tacs->createVec();
  mode->incref();
  TacsScalar *mode_vals;
  mode->getArray(&mode_vals);

  if (mpiSize > 1){
    fprintf(stderr, "TACSIntegrator: Natural frequencies \
      can only be determined in serial\n");
    return -1;
  }

  // Determine the size of state vector
  int num_state_vars;
  q->getSize(&num_state_vars);

  // Set the (steady-state) state variables into TACS
  tacs->setVariables(q, qdot, qddot);

  // Create K matrix
  FEMat *DK = tacs->createFEMat();  DK->incref();
  FEMat *DG = tacs->createFEMat();  DG->incref();
  FEMat *DM = tacs->createFEMat();  DM->incref();
  tacs->assembleJacobian(1.0, 0.0, 0.0, NULL, DK);
  tacs->assembleJacobian(0.0, 1.0, 0.0, NULL, DG);
  tacs->assembleJacobian(0.0, 0.0, 1.0, NULL, DM);

  // Get the dense column-major orientations of the matrix
  BCSRMat *Kbcsr, *Gbcsr, *Mbcsr;
  DK->getBCSRMat(&Kbcsr, NULL, NULL, NULL);
  DG->getBCSRMat(&Gbcsr, NULL, NULL, NULL);
  DM->getBCSRMat(&Mbcsr, NULL, NULL, NULL);

  // Retrieve the values from the matrices
  TacsScalar *Kvals = new TacsScalar[ num_state_vars*num_state_vars ];
  TacsScalar *Gvals = new TacsScalar[ num_state_vars*num_state_vars ];
  TacsScalar *Mvals = new TacsScalar[ num_state_vars*num_state_vars ];
  Kbcsr->getDenseColumnMajor(Kvals);
  Gbcsr->getDenseColumnMajor(Gvals);
  Mbcsr->getDenseColumnMajor(Mvals);

  int index = 0;
  if (use_gyroscopic){
    // Create an A/B matrix for the generalized eigenproblem
    int n = 2*num_state_vars;
    double *A = new double[ n*n ];
    double *B = new double[ n*n ];
    memset(A, 0, n*n*sizeof(double));
    memset(B, 0, n*n*sizeof(double));

    // Call LAPACK
    double *alphar = new double[ n ];
    double *alphai = new double[ n ];
    double *beta = new double[ n ];
    double *vl = NULL, *vr = NULL;
    int lwork = 20*n;
    double *work = new double[ lwork ];
    int info = 0;

    // A*x = lambda*B*x
    const int nvars = num_state_vars;
    for ( int j = 0; j < nvars; j++ ){
      for ( int i = 0; i < nvars; i++ ){
        // Set the A-matrix:
        // [ 0  -K ]
        // [ K   G ]
        A[i+nvars + n*j] = TacsRealPart(Kvals[i + nvars*j]);
        A[i + n*(j + nvars)] = -TacsRealPart(Kvals[i + nvars*j]);
        A[i+nvars + n*(j + nvars)] = TacsRealPart(Gvals[i + nvars*j]);

        // Set the B-matrix:
        // [ K  0 ]
        // [ 0  M ]
        B[i + n*j] = TacsRealPart(Kvals[i + nvars*j]);
        B[i+nvars + n*(j + nvars)] = TacsRealPart(Mvals[i + nvars*j]);
      }
    }

    // Call lapack to solve the eigenvalue problem
    if (modes){
      vr = new double[ n*n];
      LAPACKdggev("N", "V", &n, A, &n, B, &n,
                  alphar, alphai, beta,
                  vl, &n, vr, &n, work, &lwork, &info);
    }
    else {
      LAPACKdggev("N", "N", &n, A, &n, B, &n,
                  alphar, alphai, beta,
                  vl, &n, vr, &n, work, &lwork, &info);
    }

    // Print the eigenvalues
    for ( int i = 0; i < n; i++ ){
      if (fabs(beta[i]) > 1e-14 && alphai[i] > 0.0){
        freq[index] = alphai[i]/beta[i];
        //  Get the corresponding eigenvector. If the j-th eigenvalue
        // is real, then v(j) = VR(:,j), the j-th column of VR. If the
        // j-th and (j+1)-th eigenvalues form a complex conjugate
        // pair, then v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) =
        // VR(:,j)-i*VR(:,j+1).
        if (modes){
          for ( int k = 0; k < n; k++ ){
            modes[index*n+k] = vr[i*n+k];
            mode_vals[k] = vr[i*n+k];
          }
          // Write the mode to disk as f5
          tacs->setVariables(mode, mode, mode);
          if (beamf5){
            char fname[256];
            sprintf(fname, "mode_freq_%g.f5", TacsRealPart(freq[index]));
            beamf5->writeToFile(fname);
          }
        }
        index++;
        i++;
      }
    }

    delete [] A;
    delete [] B;
    delete [] alphar;
    delete [] alphai;
    delete [] beta;
    delete [] work;
  }
  else {
    int n = num_state_vars;
    double *A = new double[ n*n ];
    double *B = new double[ n*n ];
    memset(A, 0, n*n*sizeof(double));
    memset(B, 0, n*n*sizeof(double));

    // Call LAPACK
    double *alphar = new double[ n ];
    double *alphai = new double[ n ];
    double *beta = new double[ n ];
    double *vl = NULL, *vr = NULL;
    int lwork = 20*n;
    double *work = new double[ lwork ];
    int info = 0;

    for ( int j = 0; j < n; j++ ){
      for ( int i = 0; i < n; i++ ){
        A[i + n*j] = TacsRealPart(Kvals[i + n*j]);
        B[i + n*j] = TacsRealPart(Mvals[i + n*j]);
      }

    }

    // Call lapack to solve the eigenvalue problem
    if (modes){
      vr = new double[ n*n];
      LAPACKdggev("N", "V", &n, A, &n, B, &n,
                  alphar, alphai, beta,
                  vl, &n, vr, &n, work, &lwork, &info);
    }
    else {
      LAPACKdggev("N", "N", &n, A, &n, B, &n,
                  alphar, alphai, beta,
                  vl, &n, vr, &n, work, &lwork, &info);
    }

    // Print the eigenvalues K v = lam M v
    for ( int i = 0; i < n; i++ ){
      if (fabs(beta[i]) > 1e-14 && alphar[i] > 0.0){
        freq[index] = sqrt(alphar[i]/beta[i]);

        // Get the eigenvector corresponding to this eigenvalue. If
        // the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th
        // column of VR
        if (modes){
          for ( int k = 0; k < n; k++ ){
            modes[index*n+k] = vr[i*n+k];
            mode_vals[k] = vr[i*n+k];
          }
          // Write the mode to disk as f5
          tacs->setVariables(mode, mode, mode);
          if (beamf5){
            char fname[256];
            sprintf(fname, "mode_freq_%g.f5", TacsRealPart(freq[index]));
            beamf5->writeToFile(fname);
          }
        }
        index++;
      }
    }

    if (vr){ delete [] vr; };
    delete [] A;
    delete [] B;
    delete [] alphar;
    delete [] alphai;
    delete [] beta;
    delete [] work;
  }

  mode->decref();
  DK->decref();
  DG->decref();
  DM->decref();

  delete [] Kvals;
  delete [] Gvals;
  delete [] Mvals;

  return index;
}

/*
  Solves the linear system Ax=b using LAPACK. The execution should be
  in serial mode.
*/
void TACSIntegrator::lapackLinearSolve( TACSBVec *res,
                                        TACSMat *mat,
                                        TACSBVec *update ) {
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
      fprintf(stderr,"LAPACK GETRF: %d-th argument had an illegal value\n",
              info);
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
  if (step_num == num_time_steps){
    time_forward = MPI_Wtime() - time_forward;
  }

  // Write the tecplot output to disk if sought
  if (f5_write_freq > 0 && step_num % f5_write_freq == 0){
    writeStepToF5(step_num);
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
                                  TACSBVec **_q,
                                  TACSBVec **_qdot,
                                  TACSBVec **_qddot){
  if (_q){
    *_q = q[step_num];
  }
  if (_qdot){
    *_qdot = qdot[step_num];
  }
  if (_qddot){
    *_qddot = qddot[step_num];
  }
  return time[step_num];
}

/*
  Creates mat, ksm and pc objects
*/
void TACSIntegrator::initializeLinearSolver( ){

  // Return if already initialized
  if (linear_solver_initialized == 1){
    return;  
  };
  
  if (!ksm){
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
        fprintf(stderr,
                "TACSIntegrator error: Using SerialBCSCMat in parallel\n");
      }
    }

    // The Krylov subspace method (KSM) associated with the solver
    ksm = new GMRES(mat, pc, gmres_iters, num_restarts, is_flexible);
    ksm->incref();
  }
  else {
    ksm->getOperators(&mat, &pc);
  }

  // ksm->setMonitor(new KSMPrintStdout("GMRES", 0, 1));
  ksm->setTolerances(0.1*rtol, 1.0e-30);

  // Set the global variable to initialize linear solver
  linear_solver_initialized = 1;
}

/*
  Function that returns the finite-difference/complex-step gradient
  (used for testing purposes)
*/
void TACSIntegrator::checkGradients( double dh ){
 // Check whether the function has been set properly
  if (num_funcs == 0 || funcs == NULL) {
    fprintf(stderr,
            "TACS Warning: Function is not set, skipping adjoint solve. \n");
    return;
  }

  // Get the design variable values
  TacsScalar *x = new TacsScalar[ num_design_vars ];
  tacs->getDesignVars(x, num_design_vars);

  // Create a temporary vector of design variable values
  TacsScalar *xtmp = new TacsScalar[ num_design_vars ];
  memcpy(xtmp, x, num_design_vars*sizeof(TacsScalar));

  // Allocate an array of function values
  TacsScalar *fvs = new TacsScalar[ num_funcs ];
  TacsScalar *ftmp = new TacsScalar[ num_funcs ];
  TacsScalar *dfp = new TacsScalar[ num_funcs ];

  // Allocate a temporary vector
  TACSBVec *Xpts = tacs->createNodeVec();
  TACSBVec *Xdup = tacs->createNodeVec();
  TACSBVec *Xtmp = tacs->createNodeVec();
  Xpts->incref();
  Xdup->incref();
  Xtmp->incref();

  // Get the nodal vector
  tacs->getNodes(Xdup);
  tacs->getNodes(Xpts);

  // Create a random vector
  Xtmp->setRand(-1.0, 1.0);

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
  evalFunctions(ftmp);

  // Compute the total derivative
  for ( int k = 0; k < num_funcs; k++ ){
    ftmp[k] = TacsImagPart(ftmp[k])/dh;
  }

  // Evaluate the derivative
  integrateAdjoint();
#else // REAL code
  // Integrate forward in time
  integrate();

  // Evaluate the functions
  evalFunctions(fvs);

  // Evaluate the derivative
  integrateAdjoint();

  // Set the design variables
  for ( int k = 0; k < num_design_vars; k++ ){
    xtmp[k] = x[k] + dh;
  }
  tacs->setDesignVars(xtmp, num_design_vars);

  // Integrate forward in time
  integrate();

  // Evaluate the functions
  evalFunctions(ftmp);

  // Compute the total derivative
  for ( int k = 0; k < num_funcs; k++ ){
    ftmp[k] = (ftmp[k] - fvs[k])/dh;
  }
#endif // TACS_USE_COMPLEX
  tacs->setDesignVars(x, num_design_vars);

  // Compute the total projected derivative
  for ( int i = 0; i < num_funcs; i++ ){
    dfp[i] = 0.0;
    for ( int j = 0; j < num_design_vars; j++ ){
      dfp[i] += dfdx[num_design_vars*i + j];
    }
  }

  printf("%3s %15s %15s %15s\n", "Fn", "Adjoint", "FD", "Relative Err");
  for ( int i = 0; i < num_funcs; i++ ){
    double relerr = fabs(TacsRealPart((dfp[i] - ftmp[i])/ftmp[i]));
    printf("%3d %15.8e %15.8e %15.8e\n",
           i, TacsRealPart(dfp[i]), TacsRealPart(ftmp[i]), relerr);
  }

  // Now test the derivative w.r.t. node locations
#ifdef TACS_USE_COMPLEX
  // Perturb the nodes
  Xpts->axpy(TacsScalar(0.0, dh), Xtmp);
  tacs->setNodes(Xpts);

  // Integrate forward in time
  integrate();

  // Evaluate the functions
  evalFunctions(ftmp);

  // Compute the total derivative
  for ( int k = 0; k < num_funcs; k++ ){
    ftmp[k] = TacsImagPart(ftmp[k])/dh;
  }
#else // REAL code
  // Perturb the nodes
  Xpts->axpy(dh, Xtmp);
  tacs->setNodes(Xpts);

  // Integrate forward in time
  integrate();

  // Evaluate the functions
  evalFunctions(ftmp);

  // Compute the total derivative
  for ( int k = 0; k < num_funcs; k++ ){
    ftmp[k] = (ftmp[k] - fvs[k])/dh;
  }
#endif // TACS_USE_COMPLEX

  // Add the contribution from the nodal derivatives
  for ( int i = 0; i < num_funcs; i++ ){
    dfp[i] = dfdXpt[i]->dot(Xtmp);
  }

  printf("%3s %15s %15s %15s\n", "Fn", "Adjoint", "FD", "Relative Err");
  for ( int i = 0; i < num_funcs; i++ ){
    double relerr = fabs(TacsRealPart((dfp[i] - ftmp[i])/ftmp[i]));
    printf("%3d %15.8e %15.8e %15.8e\n",
           i, TacsRealPart(dfp[i]), TacsRealPart(ftmp[i]), relerr);
  }

  // Reset the nodes/design variable values
  tacs->setNodes(Xdup);

  // Free the allocated data
  Xpts->decref();
  Xdup->decref();
  Xtmp->decref();
  delete [] x;
  delete [] fvs;
  delete [] ftmp;
  delete [] dfp;
  delete [] xtmp;
}

/*
  Constructor for BDF Integration scheme

  Input:
  tinit: the initial time
  tfinal: the final time
  num_steps: the number of steps to take for each second
  max_bdf_order: global order of accuracy
*/
TACSBDFIntegrator::TACSBDFIntegrator( TACSAssembler * _tacs,
                                      double _tinit, double _tfinal,
                                      double _num_steps,
                                      int _max_bdf_order):
TACSIntegrator(_tacs, _tinit,  _tfinal,  _num_steps){
  if (mpiRank == 0){
    fprintf(logfp, "[%d] Creating TACSIntegrator of type %s order %d\n",
            mpiRank, "BDF", _max_bdf_order);
  }

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
  *nbdf = getBDFCoeff(k, bdf, order)+1;
  *nbddf = 2*(*nbdf)-1;
  if (*nbddf > k+1){
    *nbddf = k+1;
  }

  // Now, compute the second-order coefficients
  for ( int j = 0; j < *nbdf && (k - j > 0); j++ ){
    // order >= 1 always
    int order = (k-j < max_order ? k-j : max_order);
    double bdf2[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    int len = getBDFCoeff(k - j, bdf2, order)+1;

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
int TACSBDFIntegrator::getBDFCoeff( const int k, double bdf[], int order ){
  if (order <= 1){
    double h = time[k] - time[k-1];
    bdf[0] = 1.0/h;
    bdf[1] = -bdf[0];
    return 1;
  }
  else if (order == 2){
    double h = time[k] - time[k-1];
    double h1 = time[k-1] - time[k-2];
    bdf[0] = (2.0*h + h1)/(h*(h + h1));
    bdf[1] = -(h + h1)/(h*h1);
    bdf[2] = h/(h1*(h + h1));
    return 2;
  }
  else if (order >= 3){
    double h = time[k] - time[k-1];
    double h1 = time[k-1] - time[k-2];
    double h2 = time[k-2] - time[k-3];
    bdf[0] = (h1*h2 + 2*h*h2 + h1*h1 + 4*h*h1 + 3*h*h)/(
              (h*(h1 + h)*(h2 + h1 + h)));
    bdf[1] = -(h1*h2 + h*h2 + h1*h1 + 2*h*h1 + h*h)/(
               (h*h1*(h2 + h1)));
    bdf[2] = (h*h2 + h*h1 + h*h)/(h1*(h1 + h)*h2);
    bdf[3] = -(h*h1 + h*h)/(h2*(h2 + h1)*(h2 + h1 + h));
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
    // Output the results at the initial condition if configured
    printOptionSummary();

    // Retrieve the initial conditions and set into TACS
    tacs->getInitConditions(q[0], qdot[0], qddot[0]);
    tacs->setVariables(q[0], qdot[0], qddot[0]);

    // Solve for acceleration and set into TACS
    logTimeStep(k);
    initAccelerationSolve(forces);

    return 0;
  }

  // Extrapolate to next time step: q[k] = q[k-1] + h*qdot[k-1] +
  // h^2/2*qddot[k-1] (helps reduce the initial residual)
  double h = time[k] - time[k-1];
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
    double scale = bdf_coeff[i];
    qdot[k]->axpy(scale, q[k-i]);
  }

  // approximate qddot using BDF formula
  qddot[k]->zeroEntries();
  for ( int i = 0; i < nbddf; i++ ){
    double scale = bddf_coeff[i];
    qddot[k]->axpy(scale, q[k-i]);
  }

  // If required, add the contribution to the second derivative
  // from the initial values of the first derivatives
  if (k == nbdf-1){
    double scale = bdf_coeff[nbdf-1];
    qddot[k]->axpy(scale, qdot[0]);
  }

  // Determine the coefficients for linearizing the Residual
  double alpha = 1.0;
  double beta = bdf_coeff[0];
  double gamma = bddf_coeff[0];

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
  //TODO: integrate functions as the forward problem is integrated
  // Check whether these are two-stage or single-stage functions
  int twoStage = 0;
  for ( int n = 0; n < num_funcs; n++ ){
    if (funcs[n] && funcs[n]->getStageType() == TACSFunction::TWO_STAGE){
      twoStage = 1;
      break;
    }
  }

  // Initialize the function if had already not been initialized
  if (twoStage){
    // First stage
    for ( int n = 0; n < num_funcs; n++ ){
      if (funcs[n]){
        funcs[n]->initEvaluation(TACSFunction::INITIALIZE);
      }
    }

    for ( int k = start_plane; k <= end_plane; k++ ){
      // Set the stages
      tacs->setSimulationTime(time[k]);
      tacs->setVariables(q[k], qdot[k], qddot[k]);

      double tcoeff = 0.0;
      if (k > start_plane && k <= end_plane){
        tcoeff += 0.5*(time[k] - time[k-1]);
      }
      if (k >= start_plane && k < end_plane){
        tcoeff += 0.5*(time[k+1] - time[k]);
      }
      tacs->integrateFunctions(tcoeff, TACSFunction::INITIALIZE,
                               funcs, num_funcs);
    }

    for ( int n = 0; n < num_funcs; n++ ){
      if (funcs[n]){
        funcs[n]->finalEvaluation(TACSFunction::INITIALIZE);
      }
    }
  }

  // Second stage
  for ( int n = 0; n < num_funcs; n++ ){
    if (funcs[n]){
      funcs[n]->initEvaluation(TACSFunction::INTEGRATE);
    }
  }

  for ( int k = start_plane; k <= end_plane; k++ ){
    tacs->setSimulationTime(time[k]);
    tacs->setVariables(q[k], qdot[k], qddot[k]);

    double tcoeff = 0.0;
    if (k > start_plane && k <= end_plane){
      tcoeff += 0.5*(time[k] - time[k-1]);
    }
    if (k >= start_plane && k < end_plane){
      tcoeff += 0.5*(time[k+1] - time[k]);
    }
    tacs->integrateFunctions(tcoeff, TACSFunction::INTEGRATE,
                             funcs, num_funcs);
  }

  for ( int n = 0; n < num_funcs; n++ ){
    if (funcs[n]){
      funcs[n]->finalEvaluation(TACSFunction::INTEGRATE);
    }
  }

  // Retrieve the function values
  for ( int n = 0; n < num_funcs; n++ ){
    fvals[n] = 0.0;
    if (funcs[n]){
      fvals[n] = funcs[n]->getFunctionValue();
    }
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

  // Zero the right-hand-sides at the last time-step
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

    // Zero the derivative!
    memset(dfdx, 0, num_funcs*num_design_vars*sizeof(TacsScalar));

    // Zero the derivatives w.r.t. the node locations
    for ( int i = 0; i < num_funcs; i++ ){
      dfdXpt[i]->zeroEntries();
    }

    // Initialize linear solver
    initializeLinearSolver();
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

    // Compute the time coefficient for the integral of
    // the functional
    double tcoeff = 0.0;
    if (k > start_plane && k <= end_plane){
      tcoeff += 0.5*(time[k] - time[k-1]);
    }
    if (k >= start_plane && k < end_plane){
      tcoeff += 0.5*(time[k+1] - time[k]);
    }

    // Determine the linearization coefficients for Jacobian Assembly
    double alpha = 1.0;
    double beta = bdf_coeff[0];
    double gamma = bddf_coeff[0];

    // Find the adjoint index
    int adj_index = k % num_adjoint_rhs;

    // Setup the adjoint RHS
    if (k >= start_plane && k <= end_plane){
      for ( int n = 0; n < num_funcs; n++ ){
        if (funcs[n]){
          // Add up the contribution from function state derivative to RHS
          tacs->addSVSens(tcoeff, 0.0, 0.0, &funcs[n], 1,
                          &rhs[adj_index*num_funcs + n]);
        }
      }
    }

    // Scale the right-hand-side by -1
    for ( int n = 0; n < num_funcs; n++ ){
      rhs[adj_index*num_funcs + n]->scale(-1.0);
    }

    // Setup the Jacobian
    double tassembly = MPI_Wtime();

    // Try to downcast to a multigrid pc
    TACSMg *mg = dynamic_cast<TACSMg*>(pc);
    if (mg){
      mg->assembleJacobian(alpha, beta, gamma, NULL, TRANSPOSE);
    }
    else {
      tacs->assembleJacobian(alpha, beta, gamma, NULL, mat, TRANSPOSE);
    }
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
      tacs->applyBCs(res);
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
  if (k > start_plane && k <= end_plane){
    tcoeff += 0.5*(time[k] - time[k-1]);
  }
  if (k >= start_plane && k < end_plane){
    tcoeff += 0.5*(time[k+1] - time[k]);
  }
  if (k >= start_plane && k <= end_plane){
    tacs->addDVSens(tcoeff, funcs, num_funcs, dfdx, num_design_vars);
    tacs->addXptSens(tcoeff, funcs, num_funcs, dfdXpt);
  }

  if (k > 0){
    // Get the BDF coefficients at this time step
    int nbdf, nbddf;
    double bdf_coeff[4];
    double bddf_coeff[9];
    get2ndBDFCoeff(k, bdf_coeff, &nbdf, bddf_coeff,
                   &nbddf, max_bdf_order);

    // Add total derivative contributions from this step to all
    // functions
    double jacpdt = MPI_Wtime();
    tacs->addAdjointResProducts(1.0, psi, num_funcs,
                                dfdx, num_design_vars);
    tacs->addAdjointResXptSensProducts(1.0, psi, num_funcs, dfdXpt);
    time_rev_jac_pdt += MPI_Wtime() - jacpdt;

    // Drop the contributions from this step to other right hand sides
    double tassembly2 = MPI_Wtime();
    for ( int ii = 1; (ii < nbdf || ii < nbddf); ii++ ){
      int rhs_index = (k - ii) % num_adjoint_rhs;
      double beta = 0.0, gamma = 0.0;
      if (ii < nbdf){
        beta = bdf_coeff[ii];
      }
      if (ii < nbddf){
        gamma = bddf_coeff[ii];
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
    // Finally sum up all of the results across all processors
    for ( int n = 0; n < num_funcs; n++ ){
      dfdXpt[n]->beginSetValues(TACS_ADD_VALUES);
    }
    for ( int n = 0; n < num_funcs; n++ ){
      dfdXpt[n]->endSetValues(TACS_ADD_VALUES);
    }

    // All reduce the contributions across processors
    MPI_Allreduce(MPI_IN_PLACE, dfdx, num_funcs*num_design_vars,
                  TACS_MPI_TYPE, MPI_SUM, tacs->getMPIComm());

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

/*
  Constructor for TACSDIRKIntegrator

  Input:

  num_stages:        the number of Runge-Kutta stages
  tinit:             the initial time
  tfinal:            the final time
  num_steps:         the number of steps to take
  order:             order of the truncation error
*/
TACSDIRKIntegrator::TACSDIRKIntegrator( TACSAssembler * _tacs,
                                        double _tinit, double _tfinal,
                                        double _num_steps,
                                        int _num_stages ):
TACSIntegrator(_tacs, _tinit, _tfinal, _num_steps){
  if (mpiRank == 0){
    fprintf(logfp, "[%d] Creating TACSIntegrator of type %s stages %d\n",
            mpiRank, "DIRK", _num_stages);
  }
  // Set the number of stages
  num_stages = _num_stages;

  // allocate space for stage state variables
  qS = new TACSBVec*[ num_stages*num_time_steps ];
  qdotS = new TACSBVec*[ num_stages*num_time_steps ];
  qddotS = new TACSBVec*[ num_stages*num_time_steps ];

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
  a = new double[num_stages*(num_stages+1)/2];
  b = new double[num_stages];
  c = new double[num_stages];
  A = new double[num_stages*(num_stages+1)/2];
  B = new double[num_stages];

  // Set the Butcher tableau entries to zero
  memset(a, 0, num_stages*(num_stages+1)/2*sizeof(double));
  memset(b, 0, num_stages*sizeof(double));
  memset(c, 0, num_stages*sizeof(double));
  memset(A, 0, num_stages*(num_stages+1)/2*sizeof(double));
  memset(B, 0, num_stages*sizeof(double));

  // Add entries into the Butcher tableau
  setupDefaultCoeffs();

  // Set the second-order coefficients
  setupSecondCoeffs();

  // NULL-out the adjoint stuff
  rhs = NULL;
  lambda = NULL;
  omega = NULL;
  domega = NULL;
  phi = NULL;
  psi = NULL;
}

/*
  Destructor for TACSDIRKIntegrator
*/
TACSDIRKIntegrator::~TACSDIRKIntegrator(){
  // Cleanup Butcher's Tableau
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] A;
  delete [] B;

  // Cleanup stage states
  for ( int i = 0; i < num_stages*num_time_steps; i++ ){
    qS[i]->decref();
    qdotS[i]->decref();
    qddotS[i]->decref();
  }

  delete [] qS;
  delete [] qdotS;
  delete [] qddotS;

  // Free the data that was allocated
  if (lambda){
    rhs->decref();

    for ( int i = 0; i < num_funcs; i++ ){
      phi[i]->decref();
      psi[i]->decref();
    }

    for ( int i = 0; i < num_funcs*num_stages; i++ ){
      lambda[i]->decref();
      omega[i]->decref();
      domega[i]->decref();
    }

    delete [] lambda;
    delete [] phi;
    delete [] psi;
    delete [] omega;
    delete [] domega;
  }
}

/*
  Function that puts the entries into Butcher tableau
*/
void TACSDIRKIntegrator::setupDefaultCoeffs(){
  if (num_stages == 1){
    // Implicit mid-point rule (A-stable)
    a[0] = 0.5;
    b[0] = 1.0;
    c[0] = 0.5;
  }
  else if (num_stages == 2){
    // Crouzeix formula (A-stable)
    double tmp = 1.0/(2.0*sqrt(3.0));
    a[0] = 0.5 + tmp;
    a[1] = -1.0/sqrt(3.0);
    a[2] = a[0];

    b[0] = 0.5;
    b[1] = 0.5;

    c[0] = 0.5 + tmp;
    c[1] = 0.5 - tmp;
  }
  else if (num_stages == 3){
    // Crouzeix formula (A-stable)
    double alpha = 2.0*cos(M_PI/18.0)/sqrt(3.0);

    a[0] = (1.0 + alpha)*0.5;
    a[1] = -0.5*alpha;
    a[2] = a[0];
    a[3] = 1.0 + alpha;
    a[4] = -(1.0 + 2.0*alpha);
    a[5] = a[0];

    b[0] = 1.0/(6.0*alpha*alpha);
    b[1] = 1.0 - 1.0/(3.0*alpha*alpha);
    b[2] = b[0];

    c[0] = 0.5*(1.0+alpha);
    c[1] = 0.5;
    c[2] = 0.5*(1.0-alpha);
  }
  else {
    fprintf(stderr, "ERROR: Invalid number of stages %d\n", num_stages);
    num_stages = 1;
    setupDefaultCoeffs();
  }

  // check for the consistency of butcher tableau entries
  checkButcherTableau();
}

/*
  Get the stage coefficient from the tableau using the full index
*/
double TACSDIRKIntegrator::getACoeff( const int i, const int j ){
  if (i < num_stages && j < num_stages){
    if (i <= j){
      int index = getRowIndex(i);
      return a[index + j];
    }
  }
  return 0.0;
}

/*
  Set the second-order coefficients based on the Butcher tableau values
*/
void TACSDIRKIntegrator::setupSecondCoeffs(){
  // Set the values of the B coefficients
  for ( int i = 0; i < num_stages; i++ ){
    B[i] = 0.0;

    // Loop over the rows in the tableau
    for ( int j = 0; j < num_stages; j++ ){
      B[i] += b[j]*getACoeff(j, i);
    }
  }

  // Set the values of the A coefficients
  for ( int i = 0; i < num_stages; i++ ){
    int index = getRowIndex(i);
    for ( int j = 0; j <= i; j++ ){
      A[index + j] = 0.0;
      for ( int k = 0; k < num_stages; k++ ){
        A[index + j] += getACoeff(i, k)*getACoeff(k, j);
      }
    }
  }
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
      tmp += a[idx];
    }

    // Check the difference
    if (fabs(c[i] - tmp) >= 1.0e-6) {
      fprintf(stderr, "WARNING: Sum A[%d,:] != c[%d] i.e. %f != %f \n",
              i, i, c[i], tmp);
    }
  }

  // Check #2: sum(B) = 1.0
  tmp = 0.0;
  for ( int i = 0; i < num_stages; i++ ){
    tmp += b[i];
  }
  if (fabs(1.0 - tmp) >= 1.0e-6) {
    fprintf(stderr, "WARNING: Sum B != 1.0 \n");
  }
}

/*
  Return the coefficients for linearizing the Residual using NBG method
*/
void TACSDIRKIntegrator::getLinearizationCoeffs( const int stage,
                                                 const double h,
                                                 double *alpha,
                                                 double *beta,
                                                 double *gamma ){
  // Starting entry of Butcher Tableau for this stage
  int index = getRowIndex(stage);

  // Compute the coefficients
  *gamma = 1.0;
  *beta = h*a[index + stage];
  *alpha = h*h*A[index + stage];
}

/*
  Start index of the Butcher Tableau A for the supplied stage
*/
int TACSDIRKIntegrator::getRowIndex( int stageNum ){
  return stageNum*(stageNum+1)/2;
}

/*
  Integration logic of DIRK. Use this function to march in time. The
  solution over time is set into the class variables q, qdot and qddot
  and time.
*/
int TACSDIRKIntegrator::iterate( int k, TACSBVec *forces ){
  if (k == 0){
    // Output the results at the initial condition if configured
    printOptionSummary();

    // Retrieve the initial conditions and set into TACS
    tacs->getInitConditions(q[0], qdot[0], qddot[0]);
    tacs->setVariables(q[0], qdot[0], qddot[0]);

    // Solve for acceleration and set into TACS
    logTimeStep(k);
    initAccelerationSolve(forces);

    return 0;
  }

  // Compute the time step
  double h = time[k] - time[k-1];

  // Compute the stage states (qS, qdotS, qddotS) based on the DIRK
  // formula
  for ( int stage = 0; stage < num_stages; stage++ ){
    // Compute the stage time
    double tS = time[k-1] + c[stage]*h;

    // Compute the offset to this stage
    int offset = (k-1)*num_stages + stage;

    // Approximate the stage states using the DIRK formula
    if (stage == 0){
      qddotS[offset]->copyValues(qddot[k-1]);
    }
    else {

      qddotS[offset]->copyValues(qddotS[offset-1]);
    }

    //
    qddotS[offset]->zeroEntries();

    // Compute approximations for qS
    qS[offset]->copyValues(q[k-1]);
    qS[offset]->axpy(h*c[stage], qdot[k-1]);

    // Compute the approximations for qdotS
    qdotS[offset]->copyValues(qdot[k-1]);

    int index = getRowIndex(stage);
    for ( int j = 0; j <= stage; j++ ){
      int prev = (k-1)*num_stages + j;
      qS[offset]->axpy(h*h*A[index+j], qddotS[prev]);
      qdotS[offset]->axpy(h*a[index+j], qddotS[prev]);
    }

    // Determine the coefficients for Jacobian assembly
    double alpha, beta, gamma;
    getLinearizationCoeffs(stage, h, &alpha, &beta, &gamma);

    // Solve the nonlinear system of stage equations starting with
    // the approximated states
    int newton_term = newtonSolve(alpha, beta, gamma, tS,
                                  qS[offset], qdotS[offset], qddotS[offset],
                                  forces);

    // Check the flag set from Newton's method to see if we have a
    // problem at this point..
    int fail = 0;
    if (newton_term < 0){
      fail = 1;
      return fail;
    }
  }

  // Compute the state varialbes at the current time step using the
  // intermediate stage states
  q[k]->copyValues(q[k-1]);
  q[k]->axpy(h, qdot[k-1]);

  qdot[k]->copyValues(qdot[k-1]);
  qddot[k]->zeroEntries();

  // Compute the new values
  for ( int stage = 0; stage < num_stages; stage++ ){
    int offset = (k-1)*num_stages + stage;
    q[k]->axpy(h*h*B[stage], qddotS[offset]);
    qdot[k]->axpy(h*b[stage], qddotS[offset]);
    qddot[k]->axpy(b[stage], qddotS[offset]);
  }

  // Perform logging, tecplot export, etc.
  logTimeStep(k);

  return 0;
}

/*
  Evaluate the functions of interest
*/
void TACSDIRKIntegrator::evalFunctions( TacsScalar *fvals ){
  // Check whether these are two-stage or single-stage functions
  int twoStage = 0;
  for ( int n = 0; n < num_funcs; n++ ){
    if (funcs[n] && funcs[n]->getStageType() == TACSFunction::TWO_STAGE){
      twoStage = 1;
      break;
    }
  }

  // Initialize the function if had already not been initialized
  if (twoStage){
    // First stage
    for ( int n = 0; n < num_funcs; n++ ){
      if (funcs[n]){
        funcs[n]->initEvaluation(TACSFunction::INITIALIZE);
      }
    }

    for ( int k = start_plane; k < end_plane; k++ ){
      // Compute the time-step
      double h = time[k+1] - time[k];

      for ( int stage = 0; stage < num_stages; stage++ ){
        double tS = time[k] + c[stage]*h;
        tacs->setSimulationTime(tS);

        // Set the offset into the stage variable vectors
        int offset = k*num_stages + stage;
        tacs->setVariables(qS[offset], qdotS[offset], qddotS[offset]);

        // Set the time-integration coefficient
        double tcoeff = h*b[stage];

        // Integrate the functions
        tacs->integrateFunctions(tcoeff, TACSFunction::INITIALIZE,
                                 funcs, num_funcs);
      }
    }

    for ( int n = 0; n < num_funcs; n++ ){
      if (funcs[n]){
        funcs[n]->finalEvaluation(TACSFunction::INITIALIZE);
      }
    }
  }

  // Second stage
  for ( int n = 0; n < num_funcs; n++ ){
    if (funcs[n]){
      funcs[n]->initEvaluation(TACSFunction::INTEGRATE);
    }
  }

  for ( int k = start_plane; k < end_plane; k++ ){
    // Compute the time-step
    double h = time[k+1] - time[k];

    for ( int stage = 0; stage < num_stages; stage++ ){
      double tS = time[k] + c[stage]*h;
      tacs->setSimulationTime(tS);

      // Set the offset into the stage variable vectors
      int offset = k*num_stages + stage;
      tacs->setVariables(qS[offset], qdotS[offset], qddotS[offset]);

      // Set the time-integration coefficient
      double tcoeff = h*b[stage];

      // Integrate the functions
      tacs->integrateFunctions(tcoeff, TACSFunction::INTEGRATE,
                               funcs, num_funcs);
    }
  }

  for ( int n = 0; n < num_funcs; n++ ){
    if (funcs[n]){
      funcs[n]->finalEvaluation(TACSFunction::INTEGRATE);
    }
  }

  // Retrieve the function values
  for ( int n = 0; n < num_funcs; n++ ){
    fvals[n] = 0.0;
    if (funcs[n]){
      fvals[n] = funcs[n]->getFunctionValue();
    }
  }
}

/*
  Set-up right-hand-sides for the adjoint equations
*/
void TACSDIRKIntegrator::initAdjoint( int step_num ){
  // Adjoint variables for each function of interest
  if (!lambda){
    rhs = tacs->createVec();
    rhs->incref();

    // Allocate the adjoint stage vectors
    lambda = new TACSBVec*[ num_funcs*num_stages ];
    for ( int i = 0; i < num_funcs*num_stages; i++ ){
      lambda[i] = tacs->createVec();
      lambda[i]->incref();
    }

    // Allocate the temporary adjoint stage-vectors
    omega = new TACSBVec*[ num_funcs*num_stages ];
    domega = new TACSBVec*[ num_funcs*num_stages ];
    for ( int i = 0; i < num_funcs*num_stages; i++ ){
      omega[i] = tacs->createVec();
      domega[i] = tacs->createVec();
      omega[i]->incref();
      domega[i]->incref();
    }

    // Allocate the inter-stage adjoint vectors
    phi = new TACSBVec*[ num_funcs ];
    psi = new TACSBVec*[ num_funcs ];
    for ( int i = 0; i < num_funcs; i++ ){
      phi[i] = tacs->createVec();
      psi[i] = tacs->createVec();
      phi[i]->incref();
      psi[i]->incref();
    }
  }

  // Zero the entries at the final step
  if (step_num == num_time_steps){
    for ( int i = 0; i < num_funcs; i++ ){
      psi[i]->zeroEntries();
      phi[i]->zeroEntries();
    }

    // Zero the derivative!
    memset(dfdx, 0, num_funcs*num_design_vars*sizeof(TacsScalar));

    // Zero the derivatives w.r.t. the node locations
    for ( int i = 0; i < num_funcs; i++ ){
      dfdXpt[i]->zeroEntries();
    }

    // Initialize linear solver
    initializeLinearSolver();
  }

  for ( int i = 0; i < num_funcs*num_stages; i++ ){
    omega[i]->zeroEntries();
    domega[i]->zeroEntries();
  }
}

/*
  Iterate to find a solution of the adjoint equations
*/
void TACSDIRKIntegrator::iterateAdjoint( int k, TACSBVec **adj_rhs ){
  if (k == 0){
    // Retrieve the initial conditions and set into TACS
    tacs->getInitConditions(q[0], qdot[0], qddot[0]);
    tacs->setVariables(q[0], qdot[0], qddot[0]);

    // Output the results at the initial condition if configured
    printOptionSummary();
    return;
  }

  // Compute the time step
  double h = time[k] - time[k-1];

  // Iterate in reverse through the stage equations
  for ( int stage = num_stages-1; stage >= 0; stage-- ){
    // Compute the stage time
    double tS = time[k-1] + c[stage]*h;

    // Compute the offset to this stage
    int offset = (k-1)*num_stages + stage;

    // Set the time step and stage variables
    tacs->setSimulationTime(tS);
    tacs->setVariables(qS[offset], qdotS[offset], qddotS[offset]);

    // Determine the coefficients for Jacobian assembly
    double alpha, beta, gamma;
    getLinearizationCoeffs(stage, h, &alpha, &beta, &gamma);

    // Setup the Jacobian
    TACSMg *mg = dynamic_cast<TACSMg*>(pc);
    if (mg){
      mg->assembleJacobian(alpha, beta, gamma, NULL, TRANSPOSE);
    }
    else {
      tacs->assembleJacobian(alpha, beta, gamma, NULL, mat, TRANSPOSE);
    }

    // Factor the preconditioner
    pc->factor();

    // Compute the derivatives and store them
    if (k > start_plane && k <= end_plane){
      tacs->addSVSens(1.0, 0.0, 0.0, funcs, num_funcs, &omega[num_funcs*stage]);
      tacs->addSVSens(0.0, 1.0, 0.0, funcs, num_funcs, &domega[num_funcs*stage]);
    }

    // Compute all the contributions to the right-hand-side
    for ( int i = 0; i < num_funcs; i++ ){
      rhs->copyValues(phi[i]);
      rhs->axpy(h*B[stage]/b[stage], psi[i]);

      // Add the derivatives to the right-hand-side
      int index = getRowIndex(stage);
      rhs->axpy(-h*h*A[index + stage], omega[num_funcs*stage + i]);
      rhs->axpy(-h*a[index + stage], domega[num_funcs*stage + i]);

      // Add up the values to get the right-hand-side
      for ( int j = stage+1; j < num_stages; j++ ){
        // a[j,stage]*b[j]/b[stage]
        index = getRowIndex(j);

        // Compute the coefficient for the domega terms
        double alpha = a[index + stage]*b[j]/b[stage];
        rhs->axpy(-h*alpha, domega[num_funcs*j + i]);

        // Compute the coefficient for the omega terms
        double beta = A[index + stage]*b[j]/b[stage];
        rhs->axpy(-h*h*beta, omega[num_funcs*j + i]);
      }

      // Apply the factorization for all right hand sides and solve
      // for the adjoint variables
      ksm->solve(rhs, lambda[num_funcs*stage + i]);
    }

    // Add the products to omega and domega
    for ( int i = 0; i < num_funcs; i++ ){
      tacs->addJacobianVecProduct(1.0, 1.0, 0.0, 0.0,
                                  lambda[num_funcs*stage + i],
                                  omega[num_funcs*stage + i],
                                  TRANSPOSE);
    }

    // Add the products of omega and domega
    for ( int i = 0; i < num_funcs; i++ ){
      tacs->addJacobianVecProduct(1.0, 0.0, 1.0, 0.0,
                                  lambda[num_funcs*stage + i],
                                  domega[num_funcs*stage + i],
                                  TRANSPOSE);
    }
  }
}

/*
  Add the contributions to the total derivative from the time-step
*/
void TACSDIRKIntegrator::postAdjoint( int k ){
  if (k > 0){
    // Compute the time step
    double h = time[k] - time[k-1];

    // Add the contributions to the total derivatives
    for ( int stage = num_stages-1; stage >= 0; stage-- ){
      // Compute the stage time
      double tS = time[k-1] + c[stage]*h;

      // Compute the offset to this stage
      int offset = (k-1)*num_stages + stage;

      // Set the time step and stage variables
      tacs->setSimulationTime(tS);
      tacs->setVariables(qS[offset], qdotS[offset], qddotS[offset]);

      double tcoeff = h*b[stage];
      if (k > start_plane && k <= end_plane){
        tacs->addDVSens(tcoeff, funcs, num_funcs, dfdx, num_design_vars);
        tacs->addXptSens(tcoeff, funcs, num_funcs, dfdXpt);
      }

      // Add the derivative of the product of the adjoint to the
      // output vector
      tacs->addAdjointResProducts(tcoeff,
                                  &lambda[num_funcs*stage], num_funcs,
                                  dfdx, num_design_vars);
      tacs->addAdjointResXptSensProducts(tcoeff,
                                         &lambda[num_funcs*stage], num_funcs,
                                         dfdXpt);
    }

    // Sum up the contributions to the phi vectors
    for ( int i = 0; i < num_funcs; i++ ){
      // Integrate phi over the last time step
      phi[i]->axpy(h, psi[i]);
      for ( int stage = 0; stage < num_stages; stage++ ){
        phi[i]->axpy(-h*b[stage], domega[num_funcs*stage + i]);
        phi[i]->axpy(-h*h*b[stage]*c[stage], omega[num_funcs*stage + i]);
      }

      // Integrate psi over the last time step
      for ( int stage = 0; stage < num_stages; stage++ ){
        psi[i]->axpy(-h*b[stage], omega[num_funcs*stage + i]);
      }
    }
  }

  if (k == 0){
    // Finally sum up all of the results across all processors
    for ( int i = 0; i < num_funcs; i++ ){
      dfdXpt[i]->beginSetValues(TACS_ADD_VALUES);
    }
    for ( int i = 0; i < num_funcs; i++ ){
      dfdXpt[i]->endSetValues(TACS_ADD_VALUES);
    }

    // All reduce the contributions across processors
    MPI_Allreduce(MPI_IN_PLACE, dfdx, num_funcs*num_design_vars,
                  TACS_MPI_TYPE, MPI_SUM, tacs->getMPIComm());
  }
}

/*
  Get the adjoint value for the given function
*/
void TACSDIRKIntegrator::getAdjoint( int step_num, int func_num,
                                     TACSBVec **adjoint ){
  *adjoint = lambda[func_num];
}
