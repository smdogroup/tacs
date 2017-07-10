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
   type              : type of integrator (BDF1-3, DIRK2-4, ABM1-6, NBG[E,2-3])
*/
TACSIntegrator* TACSIntegrator::getInstance( TACSAssembler * _tacs,
                                             double _tinit, double _tfinal, 
                                             double _num_steps_per_sec, 
                                             enum IntegratorType type ){
  if ( type == DIRK2 ) {
    return new TACSDIRKIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 2);
  } else if ( type == DIRK3 ) {
    return new TACSDIRKIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 3);
  } else if ( type == DIRK4 ) {
    return new TACSDIRKIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 4);

  } else if ( type == BDF1 ) {
    return new TACSBDFIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 1);
  } else if ( type == BDF2 ) {
    return new TACSBDFIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 2);
  } else if ( type == BDF3 ) {
    return new TACSBDFIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 3);

  } else if ( type == ABM1 ) {
    return new TACSABMIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 1);
  } else if ( type == ABM2 ) {
    return new TACSABMIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 2);
  } else if ( type == ABM3 ) {
    return new TACSABMIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 3);
  } else if ( type == ABM4 ) {
    return new TACSABMIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 4);
  } else if ( type == ABM5 ) {
    return new TACSABMIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 5);
  } else if ( type == ABM6 ) {
    return new TACSABMIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 6);

    // Newmark Beta Gamma Method
  } else if ( type == NBGE ) {
    return new TACSNBGIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 1);
  } else if ( type == NBG2 ) {
    return new TACSNBGIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 2);
  } else if ( type == NBG3 ) {
    return new TACSNBGIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 3);
    
  } else { 
    // Default BDF2 integrator
    return new TACSBDFIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 2);
  }  
}

/*
  String for enum types of integrator
*/
inline const char* getIntegratorType( enum IntegratorType type ) {
  switch (type) {
  case DIRK2: return "DIRK2";
  case DIRK3: return "DIRK3";
  case DIRK4: return "DIRK4";

  case BDF1:  return "BDF1";
  case BDF2:  return "BDF2";
  case BDF3:  return "BDF3";

  case ABM1:  return "ABM1";
  case ABM2:  return "ABM2";
  case ABM3:  return "ABM3";
  case ABM4:  return "ABM4";
  case ABM5:  return "ABM5";
  case ABM6:  return "ABM6";

  case NBGE:  return "NBGE";
  case NBG2:  return "NBG2";
  case NBG3:  return "NBG3";

  default:    return "UNKNOWN";
  }
}

/*
  String for enum types of ordering type
*/
inline const char* getOrderingType( enum TACSAssembler::OrderingType otype ) {
  switch (otype) {
  case TACSAssembler::NATURAL_ORDER:  return "NATURAL_ORDER";
  case TACSAssembler::RCM_ORDER:      return "RCM_ORDER";
  case TACSAssembler::AMD_ORDER:      return "AMD_ORDER";
  case TACSAssembler::ND_ORDER:       return "ND_ORDER";
  case TACSAssembler::TACS_AMD_ORDER: return "TACS_AMD_ORDER";
  default:                            return "UNKNOWN";
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
                                double _num_steps_per_sec ){
  // Mark that this is the global instance of integrator
  adaptive_instance    = 0;
  num_adaptive_retry   = 0;
  adaptive_step_factor = 10;

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
  print_level = 2;
  logfp = NULL;
  if (mpiRank == 0){ 
    logfp = stdout; 
  }

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
  
  // Store physical time of simulation
  time = new double[num_time_steps];
  memset(time, 0, num_time_steps*sizeof(double));

  //------------------------------------------------------------------//
  //                     Newton's method                              //
  //------------------------------------------------------------------//

  // Frequency of Jacobian recomputation during nonlinear solve
  jac_comp_freq = 1;

  // Set the default LINEAR solver
  use_lapack      = 0;
  eigensolve      = 0;
  use_line_search = 0;
  use_femat       = 1;

  // AMD reordering parameters
  lev           = 100000; 
  fill          = 10.0;  
  reorder_schur = 1;

  // KSM parameters
  gmres_iters  = 10;
  num_restarts = 0;
  is_flexible  = 0;
    
  // Default parameters for Newton solve
  max_newton_iters = 25;
  init_newton_delta = 0.0;
  atol = 1.0e-13;
  rtol = 1.0e-13;

  // Create vector for storing the residual at each Newton iteration
  res = tacs->createVec();
  res->incref();
  
  // Create a vector for storing the Newton updates
  update = tacs->createVec();
  update->incref();

  // Use TACS_AMD_ORDER by default
  ordering_type = TACSAssembler::TACS_AMD_ORDER;

  // Variables used in adjoint solve (use setFunction(...) to set these 
  num_funcs = 0;
  funcs = NULL;

  // NULL the different KSM/solver objects
  D = NULL;
  mat = NULL;
  pc = NULL;
  ksm = NULL;

  // Tecplot solution export (use configureOutput(...) to set these
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
  
  if (rigidf5){ rigidf5->incref();}
  if (shellf5){ shellf5->incref();}
  if (beamf5){ beamf5->incref();}
}

/*
  Logic for adaptive time marching when the original time solution
  fails due to failure flag from newton solver
*/
int TACSIntegrator::doAdaptiveMarching( ) {
  int k = current_time_step;
  
  // Set the old state into TACS
  this->setTACSStates(time[k-1], q[k-1], qdot[k-1], qddot[k-1]);
      
  // Create new integrator to march from tstart to tend
  double tstart        = time[k-1];
  double tend          = time[k-1] + h;
  int    steps_per_sec = num_steps_per_sec*adaptive_step_factor;
  TACSIntegrator *child = TACSIntegrator::getInstance(tacs,
                                                      tstart, tend,
                                                      steps_per_sec,
                                                      mytype);
  child->incref();
  
  // Setting parameters as relevant to the child integrator
  child->setAdaptiveInstance(1);
  child->energies[0] = energies[0];
  child->energies[1] = energies[1];
  child->init_energy = init_energy;

  // Use matrices (these are not deleted by the child)
  child->D   = D;
  child->pc  = pc;
  child->mat = mat;
  child->ksm = ksm;
  
  // Integrate to the next global step using reduced time step h/reduction_factor
  child->integrate();
  
  // Store final adaptively marched state variables into TACS
  TACSBVec  *qtmp, *qdottmp, *qddottmp;  
  qtmp     = tacs->createVec(); qtmp->incref();
  qdottmp  = tacs->createVec(); qdottmp->incref();
  qddottmp = tacs->createVec(); qddottmp->incref();

  // Copy the values into the global TACS instance
  time[k] = child->getTACSStates(qtmp, qdottmp, qddottmp);  
  q[k]->copyValues(qtmp);
  qdot[k]->copyValues(qdottmp);
  qddot[k]->copyValues(qddottmp);
        
  // Set the new state into TACS
  this->setTACSStates(time[k], q[k], qdot[k], qddot[k]);

  qtmp->decref();
  qdottmp->decref();
  qddottmp->decref();
                        
  // Replace the global nonlinear termination flag
  res_norm    = child->res_norm;
  update_norm = child->update_norm;
  energies[0] = child->energies[0];
  energies[1] = child->energies[1];
  
  int newton_flag = child->getNewtonTermFlag();
              
  child->decref();

  return newton_flag;
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
        fprintf(stderr,"[%d] Newton iteration %d, failed with NaN residual norm\n", mpiRank, niter);
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
                  (j % 8 == 7) ? (elem_ctr++) : elem_ctr, // body number (each 7 belongs to a node)
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
  Prints the wall time taken during operations in TACSIntegrator
   
  input:
  level: controls the level of detail requested in timing
  t0   : reference time to normalize the times calculated within TACSIntegrator
*/
void TACSIntegrator::printWallTime( double t0, int level ){
  if(level >= 0) { 
    fprintf(logfp, "[%d] Total            : %8.2f %6.2f\n", mpiRank, t0, t0/t0); 
    fprintf(logfp, "[%d] Integrator       : %8.2f %6.2f\n", mpiRank, time_forward + time_reverse, (time_forward +time_reverse)/t0); 
  }

  if (level >= 1) { 
    fprintf(logfp, ".[%d] Forward         :  %8.2f %6.2f\n", mpiRank, time_forward, time_forward/t0); 
  }

  if (level >= 2) {
    fprintf(logfp, "..[%d] Assembly       :   %8.2f %6.2f\n", mpiRank, time_fwd_assembly, time_fwd_assembly/t0);
    fprintf(logfp, "..[%d] Factor         :   %8.2f %6.2f\n", mpiRank, time_fwd_factor, time_fwd_factor/t0);
    fprintf(logfp, "..[%d] ApplyFac       :   %8.2f %6.2f\n", mpiRank, time_fwd_apply_factor, time_fwd_apply_factor/t0);
  }

  if (level >= 1) { 
    fprintf(logfp, ".[%d] Reverse         :  %8.2f %6.2f\n", mpiRank, time_reverse, time_reverse/t0); 
  }

  if (level >= 2) {
    fprintf(logfp, "..[%d] Assembly       :   %8.2f %6.2f\n", mpiRank, time_rev_assembly, time_rev_assembly/t0);
    fprintf(logfp, "..[%d] Factor         :   %8.2f %6.2f\n", mpiRank, time_rev_factor, time_rev_factor/t0);
    fprintf(logfp, "..[%d] ApplyFac       :   %8.2f %6.2f\n", mpiRank, time_rev_apply_factor, time_rev_apply_factor/t0);
    fprintf(logfp, "..[%d] JacVecPdt      :   %8.2f %6.2f\n", mpiRank, time_rev_jac_pdt, time_rev_jac_pdt/t0);
  }
}

/*
  Returns the termination flag of newton iterations
*/
int TACSIntegrator::getNewtonTermFlag(){
  return newton_term;
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
    getString(rbuffer, "results/rigid_%06d.f5", k);
    
    // Write the f5 file for this time step
    rigidf5->writeToFile(rbuffer);
  }

  if(shellf5 && getWriteFlag(k, f5_write_freq)){
    // Create a buffer for shell filename 
    char sbuffer[256];

    // Format the buffer based on the time step
    getString(sbuffer, "results/shell_%06d.f5", k);
    
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
  March one step and exit the integration. This is primarily for use 
  with FUNtoFEM. 
*/
void TACSIntegrator::marchOneStep( int k, TACSBVec *forces ){
  // Retrieve the initial conditions and set into TACS
  if (k == 1){
    tacs->getInitConditions(q[0], qdot[0], qddot[0]);
    tacs->setVariables(q[0], qdot[0], qddot[0]);

    /*
    // Perform an initial solve to get initial conditions compatible
    // with the initial accelerations/velocities
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    newtonSolve(alpha, beta, gamma, 0.0, q[0], qdot[0], qddot[0], 
                NULL, tacs->getInitBcMap());
    */
  }
  
  // Set the class variable
  current_time_step = k;  

  // Advance time
  time[k] = time[k-1] + h;
  
  // Approximate states and their derivatives using state approximations
  approxStates();
  
  // Determine the coefficients for linearizing the Residual
  double alpha, beta, gamma;   
  getLinearizationCoeffs(&alpha, &beta, &gamma);
  
  // Solve the nonlinear system of stage equations starting with the approximated states
  newton_term = newtonSolve(alpha, beta, gamma, 
                            time[k], q[k], qdot[k], qddot[k], forces);  
  if (newton_term < 0){
    exit(-1);
  }

  doEachTimeStep(k);  
}

/*
  Integate forward in time using the initial conditions retrieved from
  TACS
*/
void TACSIntegrator::integrate(){
  // Keep track of the time taken for foward mode
  time_forward = 0.0;
  time_fwd_assembly = 0.0;
  time_fwd_factor = 0.0;
  time_fwd_apply_factor = 0.0;
  time_newton = 0.0;
  double t0 = MPI_Wtime();

  // Perform logging, tecplot export, etc
  current_time_step = 0;
  
  // Set the initial conditions
  tacs->getInitConditions(q[0], qdot[0], qddot[0]);
  tacs->setVariables(q[0], qdot[0], qddot[0]);

  // Output the results at the initial condition
  doEachTimeStep(current_time_step);  

  // March forward in time
  for (int k = 1; k < num_time_steps; k++){
    marchOneStep(k, NULL);
  }

  // Perform logging, tecplot export, etc.
  doEachTimeStep(current_time_step);  

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
                                  TacsScalar *_dfdx ) {
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
                                    double dh ) {
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

  // Loop over time steps
  for ( int k = 1; k < num_time_steps; k++ ) {
    // Set the states into TACS
    setTACSStates(time[k], q[k], qdot[k], qddot[k]);
    
    // Evaluate and add the function values from this step
    this->addFunctions(h, funcs, numFuncs, funcVals); 
  }
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
  Evalute the functions using the integration coefficient
*/
void TACSIntegrator::addFunctions( double tcoeff, TACSFunction **funcs,
                                   int numFuncs, TacsScalar *funcVals ){
  // Check whether these are two-stage or single-stage functions
  int twoStage = 0;
  for ( int n = 0; n < numFuncs; n++ ){
    if (funcs[n]->getStageType() == TACSFunction::TWO_STAGE){
      twoStage = 1;
      break;
    }
  }

  // Initialize the function if had already not been initialized
  if (twoStage){
    for ( int n = 0; n < numFuncs; n++ ){
      funcs[n]->initEvaluation(TACSFunction::INITIALIZE);
    }
    tacs->integrateFunctions(tcoeff, TACSFunction::INITIALIZE,
                             funcs, numFuncs);
    for ( int n = 0; n < numFuncs; n++ ){
      funcs[n]->finalEvaluation(TACSFunction::INITIALIZE);
    }
  }
  
  // Perform the integration required to evaluate the function
  for ( int n = 0; n < numFuncs; n++ ){
    funcs[n]->initEvaluation(TACSFunction::INTEGRATE);
  }
  tacs->integrateFunctions(tcoeff, TACSFunction::INTEGRATE,
                           funcs, numFuncs);
  for ( int n = 0; n < numFuncs; n++ ){
    funcs[n]->finalEvaluation(TACSFunction::INTEGRATE);
  }
  
  // Retrieve the function values
  for ( int n = 0; n < numFuncs; n++ ){
    funcVals[n] += funcs[n]->getFunctionValue();
  }
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
      fprintf(logfp, "%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n", 
              "time", "tnewton", "#iters", "NwtnFlg", "|R|", "|R|/|R0|", "|dq|",
              "KE", "PE", "E0-E");
      
      // Compute the initial energy
      init_energy = energies[0] + energies[1];

      // Log the details
      fprintf(logfp, "%12.5e %12.5e %12d %12d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
              time[0], time_newton, 0, 0, 0.0, 0.0, 0.0,
              TacsRealPart(energies[0]), TacsRealPart(energies[1]),  0.0);
    }
  } 
  else {
    // Print out the time step summary
    if (logfp && print_level >= 1){

      // Need a title for total summary as details of Newton iteration
      // will overshadow this one line summary
      if (print_level == 2){
        fprintf(logfp, "%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n", 
                "time", "tnewton", "#iters", "NwtnFlg", "|R|", "|R|/|R0|", "|dq|",
                "KE", "PE", "E0-E");
      }
      
      fprintf(logfp, "%12.5e %12.5e %12d %12d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	      time[current_step], time_newton,
              niter, newton_term,
              TacsRealPart(res_norm), 
              TacsRealPart(res_norm/(rtol + init_res_norm)),
              TacsRealPart(update_norm),
              TacsRealPart(energies[0]), TacsRealPart(energies[1]), 
              TacsRealPart((init_energy - (energies[0] + energies[1]))));
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
  Configure adaptive time step recovery
*/
void TACSIntegrator::configureAdaptiveMarch( int factor, int num_retry ) {
  if ( factor > 1) {
    adaptive_step_factor = factor;
  }
  if ( num_retry > 0 ) {
    num_adaptive_retry = num_retry;
  }
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
    fprintf(fp, "%-30s %15d\n", "num_adaptive_retry", num_adaptive_retry);
    fprintf(fp, "%-30s %15d\n", "adaptive_step_factor", adaptive_step_factor);
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
    fprintf(fp, "%-30s %15s\n", "ordering_type", getOrderingType(ordering_type));
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

  // get the BDF coefficients
  get2ndBDFCoeff(k, bdf_coeff, &nbdf, bddf_coeff, &nbddf, max_bdf_order);
  
  // Extrapolate to next time step: q[k] = q[k-1] + h*qdot[k-1] +
  // h^2/2*qddot[k-1] (helps reduce the initial residual)
  q[k]->copyValues(q[k-1]);
  q[k]->axpy(h, qdot[k-1]);
  q[k]->axpy(0.5*h*h, qddot[k-1]);

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
  // Print adjoint mode summary before maching backwards
  printAdjointOptionSummary(logfp);
    
  current_time_step = num_time_steps;

  time_rev_assembly     = 0.0;
  time_rev_factor       = 0.0;
  time_rev_apply_factor = 0.0;
  time_reverse          = 0.0;

  double t0 = MPI_Wtime();
  
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

    // Evalute the function
    this->addFunctions(h, funcs, num_funcs, fvals);

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

/*
  Constructor for TACSDIRKIntegrator

  Input:

  num_stages:        the number of Runge-Kutta stages
  tinit:             the initial time
  tfinal:            the final time
  num_steps_per_sec: the number of steps to take for each second
  order:             order of the truncation error
*/
TACSDIRKIntegrator::TACSDIRKIntegrator( TACSAssembler * _tacs, 
                                        double _tinit, double _tfinal, 
                                        double _num_steps_per_sec,
                                        int _order ): 
TACSIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec){   
  // copy over the variables
  num_stages = _order - 1;

  // Set the type of integrator
  if ( num_stages == 3){
    mytype = DIRK4;
  } else if ( num_stages == 2){
    mytype = DIRK3;
  } else if ( num_stages == 1){
    mytype = DIRK2;
  } else {
    fprintf(stderr, "ERROR: Invalid number of stages %d\n", num_stages);
  }
  
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
  int i   = current_stage;     // Retrieve the current stage
  int idx = getRowIdx(i);      // Starting entry of Butcher Tableau for this stage
  
  // Compute the coefficients
  *gamma = 1.0;
  *beta  = h*A[idx+i];             // h Aii
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

  // Initial guess for qddotS
  if (i == 0){
    qddotS[toffset+i]->zeroEntries(); //copyValues(qddot[k-1]);
  }
  else {
    qddotS[toffset+i]->zeroEntries(); //copyValues(qddotS[toffset+i-1]);
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

  if ( i == 2 &&  j == 0) { 
    // contribution from third stage to first rhs

    *beta = B[2]*h*A[3];

    *alpha  = h2*B[2]*(A[0]*A[3]);
    *alpha += h2*B[2]*(A[1]*A[4]);
    *alpha += h2*B[2]*(A[3]*A[5]);

  } else if ( i == 1 &&  j == 0) { 
    // contribution from second stage to first rhs

    *beta = B[1]*h*A[1];

    *alpha  = h2*B[1]*(A[0]*A[1]);
    *alpha += h2*B[1]*(A[1]*A[2]);
      
  } else if ( i == 2 &&  j == 1) {  
    // contribution from third stage to second rhs

    *beta = B[2]*h*A[4];

    *alpha  = h2*B[2]*(A[2]*A[4]);
    *alpha += h2*B[2]*(A[4]*A[5]);
        
  } else {
    fprintf(stderr, "TACSIntegrator::Incorrect access to inter-stage coefficients...\n");
    exit(-1);
  }
}

void TACSDIRKIntegrator::marchOneStep( int step_num, TACSBVec *forces){ 
  printf("DIRK marchOneStep unimplemented");
}

/*
  Function that advances the global states q, qdot, qddot and time to
  next time step
  
  Input:
  The pointers to the global state variables q, qdot, qddot
  
  Output:
  Updated global state variables q, qdot, qddot at the time step
*/
void TACSDIRKIntegrator::computeTimeStepStates( int current_step, 
                                                TACSBVec **q, TACSBVec **qdot, TACSBVec **qddot ){
  int k       = current_step;
  int toffset = num_stages*k;

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
  qddot[k]->zeroEntries();
  for ( int j = 0; j < num_stages; j++ ){
    qddot[k]->axpy(B[j], qddotS[toffset+j]);
  }
}

/*
  Integration logic of DIRK. Use this function to march in time. The
  solution over time is set into the class variables q, qdot and qddot
  and time.
*/
void TACSDIRKIntegrator::integrate(){
  // Print a summary of the object
  if (!isAdaptiveInstance()){
    printOptionSummary(logfp);
  }

  current_time_step = 0;

  // Keep track of the time taken for foward mode
  time_forward          = 0.0;
  time_fwd_assembly     = 0.0;
  time_fwd_factor       = 0.0;
  time_fwd_apply_factor = 0.0;

  double t0 = MPI_Wtime();
  
  // Get the initial condition
  tacs->getInitConditions(q[0], qdot[0], qddot[0]);
  tacs->setVariables(q[0], qdot[0], qddot[0]);

  // Perform logging, tecplot export, etc.
  doEachTimeStep(current_time_step);   

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
      
      // Solve the nonlinear system of stage equations starting with
      // the approximated states
      newton_term = newtonSolve(alpha, beta, gamma, tS[toffset+i], 
                                qS[toffset+i], qdotS[toffset+i],
                                qddotS[toffset+i], NULL);
      if (newton_term < 0){
        exit(-1);
      }
    }
    
    // Advance the time
    time[k] = time[k-1] + h;

    // Compute the state varialbes at the current time step using the
    // intermediate stage states
    computeTimeStepStates(k, q, qdot, qddot);

    // Perform logging, tecplot export, etc.
    doEachTimeStep(k);
  }

  // Print a summary of the object
  if (!isAdaptiveInstance()){
    printOptionSummary(logfp);
  }

  // Keep track of the time taken for foward mode
  time_reverse += MPI_Wtime() - t0;
}
/*
  March backward in time and solve for the adjoint variables
*/
void TACSDIRKIntegrator::marchBackwards( ) {
  // Print adjoint mode summary before maching backwards
  printAdjointOptionSummary(logfp);

  current_time_step = num_time_steps;

  time_rev_assembly     = 0.0;
  time_rev_factor       = 0.0;
  time_rev_apply_factor = 0.0;
  time_reverse          = 0.0;

  double t0 = MPI_Wtime();

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
  TACSBVec **rhs    = new TACSBVec*[ num_funcs*num_adjoint_rhs ];
  TACSBVec **dfdq   = new TACSBVec*[ num_funcs*num_adjoint_rhs ];
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

      double tassembly = MPI_Wtime();

      // Add the contribution to function value from this stage
      this->addFunctions(h*B[i], funcs, num_funcs, fvals);
      
      // Add function contributions
      for ( int n = 0; n < num_funcs; n++ ){
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

      time_rev_assembly += MPI_Wtime() - tassembly;

      // LU factorization of the Jacobian
      double tfactor = MPI_Wtime();
      pc->factor();
      time_rev_factor += MPI_Wtime() - tfactor;

      // Apply the factorization for all right hand sides and solve for
      // the adjoint variables
      double tapply = MPI_Wtime();
      for ( int n = 0; n < num_funcs; n++ ){
	ksm->solve(rhs[i*num_funcs+n], lambda[i*num_funcs+n]);
	rhs[i*num_funcs+n]->zeroEntries();
      }
      time_rev_apply_factor += MPI_Wtime() - tapply;

      // Add total derivative contributions from this step to all
      // functions
      double jacpdt = MPI_Wtime();        
      addToTotalDerivative(h*B[i], &lambda[i*num_funcs]);
      time_rev_jac_pdt += MPI_Wtime() - jacpdt;
      
      //--------------------------------------------------------------//
      // Put the contributions from this stage to the right hand sides
      // of upcoming stages
      //--------------------------------------------------------------//

      double tassembly2 = MPI_Wtime();

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
      time_rev_assembly += MPI_Wtime() - tassembly2;
    } // stage

    double tassembly3 = MPI_Wtime();
    
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
    time_rev_assembly += MPI_Wtime() - tassembly3;
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

  // Print adjoint mode summary before maching backwards
  printAdjointOptionSummary(logfp);

  // Keep track of the time taken for foward mode
  time_reverse += MPI_Wtime() - t0;
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
      this->addFunctions(h*B[j], funcs, numFuncs, funcVals);
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
                                      double _num_steps_per_sec, 
                                      int _max_abm_order ):
TACSIntegrator(_tacs, _tinit,  _tfinal,  _num_steps_per_sec){		
  // copy over the variables
  max_abm_order = _max_abm_order;

  // Set the type of integrator
  if ( max_abm_order == 6) {
    mytype = ABM6;
  } else if ( max_abm_order == 5) {
    mytype = ABM5;
  } else if ( max_abm_order == 4) {
    mytype = ABM4;
  } else if ( max_abm_order == 3) {
    mytype = ABM3;
  } else if ( max_abm_order == 2) {
    mytype = ABM2;
  } else if ( max_abm_order == 1) {
    mytype = ABM1;
  }

  // Sanity check on the input order
  max_abm_order = ((max_abm_order <= 6 && max_abm_order >= 1) ?
		   max_abm_order : 2);
  
  // Setup ABM Coefficents
  A = new double[max_abm_order*(max_abm_order+1)/2];
  memset(A, 0, max_abm_order*(max_abm_order+1)/2*sizeof(double)); // Lower triangular matrix
  setupCoeffs();
  checkABMCoeffs();

  // As many RHS as the number of second derivative coeffs
  num_adjoint_rhs = max_abm_order;
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
  qdot and qddot. The acceleration states are the unknowns from which
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
 
  // Initial guess for qddot -- copy over previous accelearation states
  qddot[k]->zeroEntries(); //copyValues(qddot[k-1]);

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
  Return the  coefficients for linearizing the Residual using ABM method.
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
  The coefficients are in matrix form A[m,i] in row major order:
  
  [0,0]                            IDX = 0
  [1,0] [1,1]                      IDX = 1   
  [2,0] [2,1] [2,2]                IDX = 3
  [3,0] [3,1] [3,2] [3,3]          IDX = 6
  [4,0] [4,1] [4,2] [4,3] [4,4]    IDX = 9
  
*/
int TACSABMIntegrator::getCoeffIndex( int time_step ){
  // Determine the order of approximation from time step
  int m   = getOrder(time_step);

  // Get the starting index for the order from the table of coeffs
  // Order starts with 1 but the table starts at 0
  int idx = getRowIdx(m-1);   

  return idx;   
}

/*
  March backwards in time to solve for adjoint variables and computing
  total derivatives
*/
void TACSABMIntegrator::marchBackwards(){
  // Print adjoint mode summary before maching backwards
  printAdjointOptionSummary(logfp);

  // Set the current step as the last step
  current_time_step = num_time_steps;

  time_rev_assembly     = 0.0;
  time_rev_factor       = 0.0;
  time_rev_apply_factor = 0.0;
  time_reverse          = 0.0;

  double t0 = MPI_Wtime();

  // Temporary vectors for adjoint 
  TACSBVec **psibin = new TACSBVec*[ num_funcs*num_adjoint_rhs ];  
  TACSBVec **phibin = new TACSBVec*[ num_funcs*num_adjoint_rhs ]; 
  TACSBVec **rhsbin = new TACSBVec*[ num_funcs*num_adjoint_rhs ]; 

  // Adjoint vectors
  TACSBVec **psi    = new TACSBVec*[ num_funcs ];
  TACSBVec **phi    = new TACSBVec*[ num_funcs ];
  TACSBVec **lambda = new TACSBVec*[ num_funcs ];

  // Allocate vectors
  for ( int n = 0; n < num_funcs*num_adjoint_rhs; n++ ){
    if (n < num_funcs){
      phi[n] = tacs->createVec();
      phi[n]->incref();
      
      psi[n] = tacs->createVec();
      psi[n]->incref();    
      
      lambda[n] = tacs->createVec();
      lambda[n]->incref();      
    }
    phibin[n] = tacs->createVec();
    phibin[n]->incref();

    psibin[n] = tacs->createVec();
    psibin[n]->incref();     

    rhsbin[n] = tacs->createVec();
    rhsbin[n]->incref();
  }

  // March backwards in time (initial condition not evaluated)
  for ( int k = num_time_steps-1; k > 0; k-- ){

    // Set the current step as k
    current_time_step = k;

    // Integration order for this time step
    int p   = getOrder(k);

    // Pointer index into the coefficient matrix
    int idx = getCoeffIndex(k);

    // Determine the coefficients for Jacobian assembly
    double alpha, beta, gamma;
    getLinearizationCoeffs(&alpha, &beta, &gamma);

    // Set the states into the tacs instance
    setTACSStates(time[k], q[k], qdot[k], qddot[k]);
    
    // Find the adjoint index pointer
    int adj_index = k % num_adjoint_rhs;

    //---------------------------------------------------------------//
    // 1. Solve for PHI (no linear solves here)
    //---------------------------------------------------------------//

    for ( int n = 0; n < num_funcs; n++ ){
      phi[n]->copyValues(phibin[adj_index*num_funcs+n]);
    }
 
    // Zero the current phibin at adjoint index for further use next time
    for ( int n = 0; n < num_funcs; n++ ){
      phibin[adj_index*num_funcs+n]->zeroEntries();
    }

    //---------------------------------------------------------------//
    // 2. Solve for PSI (no linear solves here)
    //---------------------------------------------------------------//

    for ( int n = 0; n < num_funcs; n++ ){
      psi[n]->copyValues(psibin[adj_index*num_funcs+n]);
    }
 
    // Zero the current psibin at adjoint index for further use next time
    for ( int n = 0; n < num_funcs; n++ ){
      psibin[adj_index*num_funcs+n]->zeroEntries();
    }

    //---------------------------------------------------------------//
    // 3. Setup the adjoint RHS (involves linear solve)
    //---------------------------------------------------------------//
   
    double tassembly = MPI_Wtime();

    // Evaluate all functions (should be done before all sens calls)
    this->addFunctions(h, funcs, num_funcs, fvals);

    // Add the contribution from dfdq to RHS of the corresponding
    // adjoint index
    tacs->addSVSens(alpha, beta, gamma, 
                    funcs, num_funcs, 
                    &rhsbin[adj_index*num_funcs]);
  
    // Add the contributions from PSI
    for ( int n = 0; n < num_funcs; n++ ){
      rhsbin[adj_index*num_funcs+n]->axpy(beta/h, psi[n]);
    }
    
    // Add the contributions for PHI 
    for ( int n = 0; n < num_funcs; n++ ){
      rhsbin[adj_index*num_funcs+n]->axpy(alpha/h, phi[n]);
    }
    
    // Negate the right hand side of the adjoint index
    for ( int n = 0; n < num_funcs; n++ ){
      rhsbin[adj_index*num_funcs+n]->scale(-1.0);
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
      ksm->solve(rhsbin[adj_index*num_funcs+n], lambda[n]);
    }
    time_rev_apply_factor += MPI_Wtime() - tapply;

    // Zero the current rhs at adjoint index for further use next time
    for ( int n = 0; n < num_funcs; n++ ){
      rhsbin[adj_index*num_funcs+n]->zeroEntries();
    }

    // Add total derivative contributions from this step using adjoint
    double jacpdt = MPI_Wtime();
    addToTotalDerivative(h, lambda);
    time_rev_jac_pdt += MPI_Wtime() - jacpdt;
    
    //----------------------------------------------------------------//
    // A. Put the contribution from this step to rhsbin               //
    //----------------------------------------------------------------//
    double tassembly2 = MPI_Wtime();
    for ( int ii = 1; ii < p ; ii++ ){

      // Find the adjoint index to which the current contributions are
      // to be added
      int rhs_index = (k - ii) % num_adjoint_rhs;

      double alphatmp, betatmp, gammatmp;
      gammatmp = 0.0;
      betatmp  = h*A[idx+ii];
      alphatmp = h*A[idx]*h*A[idx+ii];

      // Add the contributions from PSI into RHS
      for ( int n = 0; n < num_funcs; n++ ){
        rhsbin[rhs_index*num_funcs+n]->axpy(betatmp/h, psi[n]);
      }

      // Add the contributions for PHI into RHS
      for ( int n = 0; n < num_funcs; n++ ){
        rhsbin[rhs_index*num_funcs+n]->axpy(alphatmp/h, phi[n]);
      }

      // Add the function and residual adjoint contributions into the
      // other rhs
      addVectorTransProducts(&rhsbin[rhs_index*num_funcs], 
                             alphatmp, betatmp, gammatmp, 
                             num_funcs, funcs, 
                             lambda);          
    }

    //----------------------------------------------------------------//
    // B. Put the contributions from this step to PSI bins            //
    //----------------------------------------------------------------//

    if ( max_abm_order > 0) {

      int rhs_index = (k-1) % num_adjoint_rhs;

      double alphatmp, betatmp, gammatmp;
      gammatmp = 0.0;
      betatmp  = h;
      alphatmp = h*h*A[idx];

      // Add the contributions for PHI into psibin
      for ( int n = 0; n < num_funcs; n++ ){
        psibin[rhs_index*num_funcs+n]->axpy(alphatmp/h, phi[n]);
      }

      // Add the contributions for PSI into psibin
      for ( int n = 0; n < num_funcs; n++ ){
        psibin[rhs_index*num_funcs+n]->axpy(betatmp/h, psi[n]);
      }      

      // Add the function and residual adjoint contributions into the
      // other psibin
      addVectorTransProducts(&psibin[rhs_index*num_funcs], 
                             alphatmp, betatmp, gammatmp,
                             num_funcs, funcs,
                             lambda);
    }

    for ( int ii = 1; ii < p ; ii++ ){

      int rhs_index = (k - ii) % num_adjoint_rhs;

      double alphatmp, betatmp, gammatmp;
      gammatmp = 0.0;
      betatmp  = 0.0;
      alphatmp = h*h*A[idx+ii];

      // Add the contributions for PHI into psibin
      for ( int n = 0; n < num_funcs; n++ ){
        psibin[rhs_index*num_funcs+n]->axpy(alphatmp/h, phi[n]);
      }

      // Add the function and residual adjoint contributions into the
      // other psibin
      addVectorTransProducts(&psibin[rhs_index*num_funcs], 
                             alphatmp, betatmp, gammatmp,
                             num_funcs, funcs,
                             lambda);               
    }

    //----------------------------------------------------------------//
    // C. Put the contributions from this step to PHI bins            //
    //----------------------------------------------------------------//

    if ( max_abm_order > 0) {

      int rhs_index = (k-1) % num_adjoint_rhs;

      double alphatmp, betatmp, gammatmp;
      gammatmp = 0.0;
      betatmp  = 0.0;
      alphatmp = h;

      // Add the contributions for PHI into phibin
      for ( int n = 0; n < num_funcs; n++ ){
        phibin[rhs_index*num_funcs+n]->axpy(alphatmp/h, phi[n]);
      }  

      // Add the function and residual adjoint contributions into the
      // other psibin
      addVectorTransProducts(&phibin[rhs_index*num_funcs], 
                             alphatmp, betatmp, gammatmp,
                             num_funcs, funcs,
                             lambda);
    }
    time_rev_assembly += MPI_Wtime() - tassembly2;

  } // end time loop backwards

  // Freeup objects
  for ( int n = 0; n < num_funcs*num_adjoint_rhs; n++ ){
    if (n < num_funcs) {  
      phi[n]->decref();
      psi[n]->decref();  
      lambda[n]->decref();
    }
    phibin[n]->decref();
    psibin[n]->decref();
    rhsbin[n]->decref();
  }
  delete [] psi;
  delete [] phi;
  delete [] lambda;

  delete [] phibin;
  delete [] psibin;
  delete [] rhsbin;

  // Print adjoint mode summary before maching backwards
  printAdjointOptionSummary(logfp);

  // Keep track of the time taken for foward mode
  time_reverse += MPI_Wtime() - t0;
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
                                      double _num_steps_per_sec,
                                      int _order):
TACSIntegrator(_tacs, _tinit,  _tfinal,  _num_steps_per_sec){

  // Set the integration order
  order = _order;
  
  // Set the type of integrator
  if ( order == 1) {
    mytype = NBGE;
  } else if ( order == 2) {
    mytype = NBG2;
  } else if ( order == 3) {
    mytype = NBG3;
  }

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

  if (mytype == NBGE) {

    // Purely Explicit (firstorder & conditionally stable rho=0)
    BETA  = 0.0;
    GAMMA = 0.0;

  } else if (mytype == NBG2) {

    // Average constant acceleration
    //    BETA = 1.0/4.0;
    //    GAMMA = 0.50;

    BETA  = 1.0/2.0;
    GAMMA = 0.50;

    //    BETA  = 1.0/8.0;
    //    GAMMA = 0.50;

    // Central Difference (second order & conditionally stable rho=2)
    //    BETA  = 1.0/12.0;
    //    GAMMA = 0.50;

  } else if (mytype == NBG3) {

    // Fox & Goodwin  (third order & conditionally stable rho=2.45)
    BETA  = 1.0/12.0;
    GAMMA = 0.50;

  }

}

/*
  Approximate states (q, qdot, qddot) at the current time step using
  the NBG formula and previous time step values of the states q, qdot
  and qddot. The acceleration states are the unknowns from which the
  velocity states and position states are obtained.
  
  Use the appropriate class variables to prepare the states (q[k],
  qdot[k], qddot[k]) at the k-th timestep for nonlinear solution. Make
  sure to setup the coefficients before calling this function.
*/
void TACSNBGIntegrator::approxStates(){
  int k = current_time_step;

  // Initial guess for qddot -- copy over previous accelearation states
  qddot[k]->zeroEntries(); //copyValues(qddot[k-1]);

  // Approximate qdot using qdot and qddot:
  // qdot[k] = qdot[k-1] + (1-GAMMA) h qddot[k-1]) + GAMMA h qddot[k]
  qdot[k]->copyValues(qdot[k-1]);
  qdot[k]->axpy((1.0-GAMMA)*h, qddot[k-1]);
  qdot[k]->axpy(GAMMA*h, qddot[k]);

  // Approximate q using q, qdot and qddot:
  // q[k] = q[k-1] + h qdot[k-1] + h^2(1-2*BETA)/2 qddot[k-1] + h^2 BETA qddot[k])
  q[k]->copyValues(q[k-1]);
  q[k]->axpy(h, qdot[k-1]);
  q[k]->axpy(h*h*(1.0-2.0*BETA)/2.0, qddot[k-1]);
  q[k]->axpy(h*h*BETA, qddot[k]);
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
void TACSNBGIntegrator::marchBackwards(){
  // Print adjoint mode summary before maching backwards
  printAdjointOptionSummary(logfp);

  current_time_step = num_time_steps;
  
  time_rev_assembly     = 0.0;
  time_rev_factor       = 0.0;
  time_rev_apply_factor = 0.0;
  time_reverse          = 0.0;

  int num_adjoint_rhs = 1; // NBG is a one step method (uses
                           // information from previous and current
                           // steps)
  double t0 = MPI_Wtime();
  
  // Temporary vectors for adjoint 
  TACSBVec **psibin = new TACSBVec*[ num_funcs*num_adjoint_rhs ];  
  TACSBVec **phibin = new TACSBVec*[ num_funcs*num_adjoint_rhs ]; 
  TACSBVec **rhsbin = new TACSBVec*[ num_funcs*num_adjoint_rhs ]; 

  // Adjoint vectors
  TACSBVec **psi    = new TACSBVec*[ num_funcs ];
  TACSBVec **phi    = new TACSBVec*[ num_funcs ];
  TACSBVec **lambda = new TACSBVec*[ num_funcs ];

  // Allocate vectors
  for ( int n = 0; n < num_funcs*num_adjoint_rhs; n++ ){
    if (n < num_funcs){
      phi[n] = tacs->createVec();
      phi[n]->incref();
      
      psi[n] = tacs->createVec();
      psi[n]->incref();    
      
      lambda[n] = tacs->createVec();
      lambda[n]->incref();      
    }

    phibin[n] = tacs->createVec();
    phibin[n]->incref();

    psibin[n] = tacs->createVec();
    psibin[n]->incref();     

    rhsbin[n] = tacs->createVec();
    rhsbin[n]->incref();
  }

  // March backwards in time (initial condition not evaluated)
  for ( int k = num_time_steps-1; k > 0; k-- ){

    // Set the current step as k
    current_time_step = k;

    // Determine the coefficients for Jacobian assembly
    double alpha, beta, gamma;
    getLinearizationCoeffs(&alpha, &beta, &gamma);

    // Set the states into the tacs instance
    setTACSStates(time[k], q[k], qdot[k], qddot[k]);
    
    // Find the adjoint index pointer
    int adj_index = k % num_adjoint_rhs;

    //---------------------------------------------------------------//
    // 1. Solve for PHI (no linear solves here)
    //---------------------------------------------------------------//

    for ( int n = 0; n < num_funcs; n++ ){
      phi[n]->copyValues(phibin[adj_index*num_funcs+n]);
    }
 
    // Zero the current phibin at adjoint index for further use next time
    for ( int n = 0; n < num_funcs; n++ ){
      phibin[adj_index*num_funcs+n]->zeroEntries();
    }

    //---------------------------------------------------------------//
    // 2. Solve for PSI (no linear solves here)
    //---------------------------------------------------------------//

    for ( int n = 0; n < num_funcs; n++ ){
      psi[n]->copyValues(psibin[adj_index*num_funcs+n]);
    }
 
    // Zero the current psibin at adjoint index for further use next time
    for ( int n = 0; n < num_funcs; n++ ){
      psibin[adj_index*num_funcs+n]->zeroEntries();
    }

    //---------------------------------------------------------------//
    // 3. Setup the adjoint RHS (involves linear solve)
    //---------------------------------------------------------------//
   
    double tassembly = MPI_Wtime();

    // Evaluate all functions (should be done before all sens calls)
    this->addFunctions(h, funcs, num_funcs, fvals);

    // Add the contribution from dfdq to RHS of the corresponding
    // adjoint index
    tacs->addSVSens(alpha, beta, gamma, 
                    funcs, num_funcs, 
                    &rhsbin[adj_index*num_funcs]);
    
    // Add the contributions from PSI
    for ( int n = 0; n < num_funcs; n++ ){
      rhsbin[adj_index*num_funcs+n]->axpy(beta/h, psi[n]);
    }
    
    // Add the contributions for PHI 
    for ( int n = 0; n < num_funcs; n++ ){
      rhsbin[adj_index*num_funcs+n]->axpy(alpha/h, phi[n]);
    }
    
    // Negate the right hand side of the adjoint index
    for ( int n = 0; n < num_funcs; n++ ){
      rhsbin[adj_index*num_funcs+n]->scale(-1.0);
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
      ksm->solve(rhsbin[adj_index*num_funcs+n], lambda[n]);
    }
    time_rev_apply_factor += MPI_Wtime() - tapply;

    // Zero the current rhs at adjoint index for further use next time
    for ( int n = 0; n < num_funcs; n++ ){
      rhsbin[adj_index*num_funcs+n]->zeroEntries();
    }

    // Add total derivative contributions from this step using adjoint
    double jacpdt = MPI_Wtime();
    addToTotalDerivative(h, lambda);
    time_rev_jac_pdt += MPI_Wtime() - jacpdt;
        
    //----------------------------------------------------------------//
    // A. Put the contribution from this step to rhsbin               //
    //----------------------------------------------------------------//

    double tassembly2 = MPI_Wtime();

    if (1) {

      // Find the adjoint index to which the current contributions are
      // to be added
      int rhs_index = k % num_adjoint_rhs;

      double alphatmp, betatmp, gammatmp;
      gammatmp = 0.0;
      betatmp  = (1.0-GAMMA)*h;
      alphatmp = 0.5*(1.0-2.0*BETA)*h*h;

      // Add the contributions from PSI into RHS
      for ( int n = 0; n < num_funcs; n++ ){
        rhsbin[rhs_index*num_funcs+n]->axpy(betatmp/h, psi[n]);
      }

      // Add the contributions for PHI into RHS
      for ( int n = 0; n < num_funcs; n++ ){
        rhsbin[rhs_index*num_funcs+n]->axpy(alphatmp/h, phi[n]);
      }

      // Add the function and residual adjoint contributions into the
      // other rhs
      addVectorTransProducts(&rhsbin[rhs_index*num_funcs], 
                             alphatmp, betatmp, gammatmp, 
                             num_funcs, funcs, 
                             lambda);          
    }

    //----------------------------------------------------------------//
    // B. Put the contributions from this step to PSI bins            //
    //----------------------------------------------------------------//

    if (1) {

      int rhs_index = k % num_adjoint_rhs;

      double alphatmp, betatmp, gammatmp;
      gammatmp = 0.0;
      betatmp  = h;
      alphatmp = h*h;

      // Add the contributions for PHI into psibin
      for ( int n = 0; n < num_funcs; n++ ){
        psibin[rhs_index*num_funcs+n]->axpy(alphatmp/h, phi[n]);
      }

      // Add the contributions for PSI into psibin
      for ( int n = 0; n < num_funcs; n++ ){
        psibin[rhs_index*num_funcs+n]->axpy(betatmp/h, psi[n]);
      }      

      // Add the function and residual adjoint contributions into the
      // other psibin
      addVectorTransProducts(&psibin[rhs_index*num_funcs], 
                             alphatmp, betatmp, gammatmp,
                             num_funcs, funcs,
                             lambda);
    }

    //----------------------------------------------------------------//
    // C. Put the contributions from this step to PHI bins            //
    //----------------------------------------------------------------//

    if (1) {

      int rhs_index = k % num_adjoint_rhs;

      double alphatmp, betatmp, gammatmp;
      gammatmp = 0.0;
      betatmp  = 0.0;
      alphatmp = h;

      // Add the contributions for PHI into phibin
      for ( int n = 0; n < num_funcs; n++ ){
        phibin[rhs_index*num_funcs+n]->axpy(alphatmp/h, phi[n]);
      }  

      // Add the function and residual adjoint contributions into the
      // other psibin
      addVectorTransProducts(&phibin[rhs_index*num_funcs], 
                             alphatmp, betatmp, gammatmp,
                             num_funcs, funcs,
                             lambda);
    }

    time_rev_assembly += MPI_Wtime() - tassembly2;
  }

  // Freeup objects
  for ( int n = 0; n < num_funcs*num_adjoint_rhs; n++ ){
    if (n < num_funcs) {  
      phi[n]->decref();
      psi[n]->decref();  
      lambda[n]->decref();
    }
    phibin[n]->decref();
    psibin[n]->decref();
    rhsbin[n]->decref();
  }
  delete [] psi;
  delete [] phi;
  delete [] lambda;

  delete [] phibin;
  delete [] psibin;
  delete [] rhsbin;

  // Print adjoint mode summary before maching backwards
  printAdjointOptionSummary(logfp);

  // Keep track of the time taken for foward mode
  time_reverse += MPI_Wtime() - t0;
}

/*
  Add Residual-Vector product and function state variable sensitivity
  into the ans vector
*/
void TACSIntegrator::addVectorTransProducts( TACSBVec **ans, 
                                             double alpha, double beta, double gamma,
                                             int num_funcs, TACSFunction **funcs,
                                             TACSBVec **input ){

  // Add the Residual-Adj vec product for each adjoint vector
  for ( int n = 0; n < num_funcs; n++ ){     
    tacs->addJacobianVecProduct(1.0, 
                                alpha, beta, gamma, 
                                input[n], ans[n], 
                                TRANSPOSE);     
  }
  
  // Add the function SV sensitivities
  tacs->addSVSens(alpha, beta, gamma, 
                  funcs, num_funcs, 
                  ans);
}
