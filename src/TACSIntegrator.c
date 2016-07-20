#include "TACSIntegrator.h"
#include <math.h>
#include "tacslapack.h"

/* 
   Static factory method that returns an instance of the concrete
   child class. Note that the base class is abstract, thus it can not
   be instantiated.
   
   Input:
   tacs              : the TACS assembler object
   tinit             : start time of simulation
   tfinal            : end time of simulation
   num_steps_per_sec : the number of steps to take per second
   type              : type of integrator (BDF1,BDF2,BDF3, DIRK2, DIRK3, DIRK4)
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

  } else { // default
    fprintf(stdout, ">> TACSIntegrator: Instantiating the default DIRK Order 2 integrator...\n");
    return new TACSDIRKIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec, 1);
  }  
}

/*
  Base class constructor for integration schemes. This base class
  contains common methods and variables pertaining to the integration
  schemes used in TACS.

  input:
  tacs: tacs assembler object
  tInit: the initial time
  tfinal: the final time
  num_steps_per_sec: the number of steps to take for each second
*/
TACSIntegrator::TACSIntegrator( TACSAssembler * _tacs,
                                double _tinit, double _tfinal, 
                                int _num_steps_per_sec ){
  // copy over the input parameters
  tacs = _tacs;
  tacs->incref();

  tinit = _tinit;
  tfinal = _tfinal;
  num_steps_per_sec = _num_steps_per_sec;

  // compute the step size
  h = 1.0/num_steps_per_sec;

  // compute the total number of time steps
  num_time_steps = int(num_steps_per_sec*(tfinal-tinit)) + 1;

  // Default print level
  print_level = 1;

  //------------------------------------------------------------------//
  //                     Time history of states                       //
  //------------------------------------------------------------------//

  // state variables that store the entire time history
  q = new BVec*[num_time_steps];
  qdot = new BVec*[num_time_steps];
  qddot = new BVec*[num_time_steps];

  // create state vectors for TACS during each timestep
  for ( int k = 0; k < num_time_steps; k++ ) {
    q[k] = tacs->createVec(); 
    q[k]->incref(); 

    qdot[k] = tacs->createVec(); 
    qdot[k]->incref(); 

    qddot[k] = tacs->createVec(); 
    qddot[k]->incref(); 
  }

  // Get the number of state variables and store into the class variable
  q[0]->getSize(&num_state_vars);
  
  // store time
  time = new double[num_time_steps];
  memset(time, 0, num_time_steps*sizeof(double));

  //------------------------------------------------------------------//
  //                     Newton's method                              //
  //------------------------------------------------------------------//

  // Frequency of Jacobian recomputation during nonlinear solve
  jac_comp_freq = 3;

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

  //------------------------------------------------------------------//
  // Variables used in adjoint solve
  //------------------------------------------------------------------//
  
  num_funcs = 0;
  funcs = NULL;
}

/*
  Destructor for base class
*/
TACSIntegrator::~TACSIntegrator(){
  // Dereference TACS
  tacs->decref();

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
  
  delete [] time;
  delete [] q;
  delete [] qdot;
  delete [] qddot;
}

/*
  Update TACS states with the supplied ones (q, qdot, qddot)
*/
void TACSIntegrator::setTACSStates( BVec *q, BVec *qdot, BVec * qddot ){
  tacs->setVariables(q);
  tacs->setDotVariables(qdot);
  tacs->setDDotVariables(qddot);
}

/*
  Adds up the contribution to the function values from this time stage/time
*/
void TACSIntegrator::addFunctionContribution( double scale ) {  
  TacsScalar ftmp;
  for ( int j = 0; j < num_funcs; j++ ) {
    tacs->evalFunctions(&funcs[j], 1, &ftmp);
    fvals[j] += scale*ftmp;
  }
}

/*
  Add up the contributions to the total derivative from this stage/stage
*/
void TACSIntegrator::addTotalDerivative( double scale, BVec **adjoint ) {
  TacsScalar *tmp = new TacsScalar[ num_design_vars*num_funcs ];
  // Add the partial derivative w.r.t. the design variables
  tacs->evalDVSens(funcs, num_funcs, tmp, num_design_vars);
  for ( int m = 0; m < num_design_vars*num_funcs; m++ ){
    dfdx[m] += scale*tmp[m];
  }

  // Add the Adjoint-Residual partial derivative contributions
  tacs->evalAdjointResProducts(adjoint, num_funcs, tmp, num_design_vars);
  for ( int m = 0; m < num_design_vars*num_funcs; m++ ){
    dfdx[m] += scale*tmp[m];
  }
  delete [] tmp;
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
                                  double t, BVec *u, BVec *udot, 
                                  BVec *uddot ){
  // Initialize the norms
  TacsScalar init_norm = 0.0;
  TacsScalar norm = 0.0;

  if (print_level > 0){
    fprintf(stdout, "%12s %8s %12s %12s %12s\n",
            "time", "Newton", "tcpu", "|R|", "|R|/|R0|");
  }
  
  double t0 = MPI_Wtime();

  // Iterate until max iters or R <= tol
  double delta = 0.0;
  for ( int n = 0; n < max_newton_iters; n++ ){
    // Set the supplied initial input states into TACS
    setTACSStates(u, udot, uddot);

    // Assemble the Jacobian matrix once in five newton iterations
    if (n % jac_comp_freq == 0){
      tacs->assembleJacobian(res, mat, alpha, beta, 
                             gamma + delta, NORMAL);
    }
    else {
      tacs->assembleRes(res);
    }
    
    // Compute the L2-norm of the residual
    norm = res->norm();
    
    // Record the residual norm at the first Newton iteration
    if (n == 0){
      init_norm = norm;
    }

    // Write a summary
    if(print_level > 0) {
      if (n == 0){
        fprintf(stdout, "%12.5e %8d %12.5e %12.5e %12.5e \n",
                t, n, MPI_Wtime()-t0, 
		RealPart(norm), RealPart(norm/init_norm));
      }
      else {
        fprintf(stdout, "%12s %8d %12.5e %12.5e %12.5e\n",
                " ", n, MPI_Wtime()-t0, 
		RealPart(norm), RealPart(norm/init_norm));
      }
    }
           
    // Check if the norm of the residuals is a NaN
    if (norm != norm){ 
      fprintf(stderr,
              "Newton iteration %d, failed with NaN residual norm\n", n);
      break;
    }
    
    // Check if the Newton convergence tolerance is satisfied
    if (norm < rtol*init_norm || norm < atol){
      break;
    }
    
    if (!use_lapack) {    
      // LU Factor the matrix when needed
      if (n % jac_comp_freq == 0){
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
    if (n == max_newton_iters && norm >= rtol*init_norm){
      fprintf(stderr,"Newton iteration failed to converge in %d iters\n", n);
      break;
    }
  }
}

/*
  Function that writes time, q, qdot, qddot to file
*/
void TACSIntegrator::writeSolution( const char *filename ) {
  // Temporary variables to access the states at each time
  TacsScalar *qvals, *qdotvals, *qddotvals;

  // Open a new file
  FILE *fp = fopen(filename, "w");
 
  for ( int k = 0; k < num_time_steps; k++ ){    
    // Copy over the state values from BVec
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
  Creates an f5 file for each time step and writes the data
*/
void TACSIntegrator::writeSolutionToF5(){
  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);

  // Create a viewer
  TACSToFH5 * f5 = new TACSToFH5(tacs, SHELL, write_flag);
  f5->incref();

  for ( int k = 0; k < num_time_steps; k++ ){    
    // Set the current states into TACS
    setTACSStates( q[k], qdot[k], qddot[k]);

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
  Compute the gradient for the given functions of interest with
  respect to the design variable.
*/
void TACSIntegrator::getFuncGrad( int _num_dv, TacsScalar *_x,
				  TacsScalar *_fvals, TacsScalar *_dfdx ) {
  // Copy the inputs
  this->num_design_vars = _num_dv;
  tacs->setDesignVars(_x, num_design_vars);

  // Check whether the function has been set properly
  if (this->num_funcs == 0 || this->funcs == NULL) {
    fprintf(stdout, 
            "TACS Warning: Function is not set, skipping adjoint solve. \n");
    return;
  }
  
  // Point the class variable to the memory allocated by the optimizer
  fvals = _fvals; dfdx = _dfdx;
  memset(fvals, 0, num_funcs*sizeof(TacsScalar));
  memset(dfdx, 0, num_funcs*num_design_vars*sizeof(TacsScalar));

  // Integrate forward in time to solve for the states (q, qdot, qddot)
  this->integrate();
  
  // March backwards in time and solve for the adjoint variables
  this->marchBackwards();
}

/*
  Function that returns the finite-difference/complex-step gradient
  (used for testing purposes)
*/
void TACSIntegrator::getFDFuncGrad( int _num_dv, TacsScalar *_x,
				    TacsScalar *_fvals, TacsScalar *_dfdx, 
				    double dh ) {
  // Copy the inputs
  this->num_design_vars = _num_dv;
  tacs->setDesignVars(_x, num_design_vars);

  // Check whether the function has been set properly
  if (this->num_funcs == 0 || this->funcs == NULL) {
    fprintf(stdout, "TACS Warning: Function is not set, skipping adjoint solve. \n");
    return;
  }
  
  // Point the class variable to the memory allocated by the optimizer
  fvals = _fvals; dfdx = _dfdx;
  memset(fvals, 0, num_funcs*sizeof(TacsScalar));
  memset(dfdx, 0, num_funcs*num_design_vars*sizeof(TacsScalar));

  TacsScalar *ftmp = new TacsScalar[num_funcs];
  memset(ftmp, 0, num_funcs*sizeof(TacsScalar));

  // Perform the forward integration if we're using a finite-difference approximation

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

  delete [] ftmp;
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
  Control the amount of information printed to the console
*/
void TACSIntegrator::setPrintLevel( int _print_level ){ 
  print_level = _print_level; 
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
  this->num_funcs = _num_funcs;
  this->funcs = _funcs;
}

/*
  Configure the F5 output 
 */
void TACSIntegrator::configureOutput(TACSToFH5 *_viewer, 
                                     int _write_freq, 
                                     char *_f5_file_fmt ) {
  this->f5             = _viewer;
  this->f5_write_freq  = _write_freq;
  this->f5_file_fmt    = _f5_file_fmt;
}

/*
  Solves the linear system Ax=b using LAPACK. The execution should be
  in serial mode.
*/
void TACSIntegrator::lapackLinearSolve( BVec *res, TACSMat *mat, BVec *update ) {
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
  Perform sanity checks on the right hand side of adjoint linear
  system
*/
void TACSIntegrator::checkAdjointRHS( BVec *rhs ) {
  // The norm of the rhs vector shouldn't be zero or shouldn't contain
  // Nan, Inf's etc
  TacsScalar norm = rhs->norm();
  if (norm != norm){ 
    fprintf(stdout, "TACS Warning: Invalid entries detected in \
the adjoint RHS vector\n");
  } 
  else if (norm <= 1.0e-15) { 
    fprintf(stdout, "TACS Warning: Zero RHS in adjoint linear system. \
The adjoint variables will be zero. Check the state variables/SVSens \
implementation.\n");
  }
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
  Constructor for BDF Integration scheme

  Input:
  tinit: the initial time
  tfinal: the final time
  num_steps_per_sec: the number of steps to take for each second
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
void TACSBDFIntegrator::approxStates( int current_step ){
  int k = current_step;

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
  Integration logic of BDF. Use this function to march in time. The
  solution over time is set into the class variables q, qdot and qddot
  and time.
*/
void TACSBDFIntegrator::integrate( ){
  // Get the initial condition
  tacs->getInitConditions(q[0], qdot[0]);
    
  for ( int k = 1; k < num_time_steps; k++ ){
    // Advance time
    time[k] = time[k-1] + h;

    // Approximate states and their derivatives using BDF formula
    approxStates(k);
    
    // Determine the coefficients for Jacobian Assembly
    double gamma = bddf_coeff[0]/(h*h);
    double beta  = bdf_coeff[0]/h;
    double alpha = 1.0;

    // Solve the nonlinear system of equations. Note that the states
    // will be advanced at the end of Newton solve
    newtonSolve(alpha, beta, gamma, time[k], q[k], qdot[k], qddot[k]);
    
    // Write the tecplot output to disk if sought
    if(f5 && f5_write_freq > 0 && k % f5_write_freq == 0){
      // Create a buffer for filename 
      char buffer[128];
      // Format the buffer based on the time step
      getString(buffer, f5_file_fmt, k);      
      // Write the f5 file for this time step
      f5->writeToFile(buffer);
    }
  }
}

/*
  Forward solve test function
*/
TacsScalar TACSBDFIntegrator::forward( const TacsScalar *x,
                                       int num_design_vars,
                                       TACSFunction *func ){
  tacs->setDesignVars(x, num_design_vars);

  // Get the initial condition
  tacs->getInitConditions(q[0], qdot[0]);
  qddot[0]->zeroEntries();

  TacsScalar fval = 0.0;

  for ( int k = 1; k < num_time_steps; k++ ){
    time[k] = time[k-1] + h;

    // Get the BDF coefficients at this time step
    get2ndBDFCoeff(k, bdf_coeff, &nbdf, 
                   bddf_coeff, &nbddf, max_bdf_order);

    // Approximate the next values
    q[k]->copyValues(q[k-1]);
    q[k]->axpy(h, qdot[k-1]);
    q[k]->axpy(h*h/2.0, qddot[k-1]);

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

    // Determine the linearization coefficients for Jacobian Assembly
    double gamma = bddf_coeff[0]/(h*h);
    double beta = bdf_coeff[0]/h;
    double alpha = 1.0;

    newtonSolve(alpha, beta, gamma, time[k],
                q[k], qdot[k], qddot[k]);


    tacs->setVariables(q[k]);
    tacs->setDotVariables(qdot[k]);
    tacs->setDDotVariables(qddot[k]);

    TacsScalar ftmp;
    tacs->evalFunctions(&func, 1, &ftmp);
    fval += h*ftmp;
  }

  return fval;
}

/*
  Adjoint solve test function
 */
void TACSBDFIntegrator::reverse( TacsScalar *dfdx,
                                 int num_design_vars,
                                 TACSFunction *func ){
  int num_rhs = (2*max_bdf_order+1) + 1;
  BVec **rhs = new BVec*[ num_rhs ];
  for ( int i = 0; i < num_rhs; i++ ){
    rhs[i] = tacs->createVec();
    rhs[i]->incref();
  }

  BVec *adj = tacs->createVec();
  adj->incref();

  TacsScalar *tmp = new TacsScalar[ num_design_vars ];
  
  memset(dfdx, 0, num_design_vars*sizeof(TacsScalar));

  for ( int k = num_time_steps-1; k >= 1; k-- ){
    // Get the BDF coefficients at this time step
    get2ndBDFCoeff(k, bdf_coeff, &nbdf, bddf_coeff, &nbddf, max_bdf_order);

    // Determine the linearization coefficients for Jacobian Assembly
    double gamma = bddf_coeff[0]/(h*h);
    double beta = bdf_coeff[0]/h;
    double alpha = 1.0;

    tacs->setVariables(q[k]);
    tacs->setDotVariables(qdot[k]);
    tacs->setDDotVariables(qddot[k]);
    
    // Find the adjoint index
    int adj_index = k % num_rhs;

    TacsScalar ftmp;
    tacs->evalFunctions(&func, 1, &ftmp);
    tacs->evalSVSens(func, adj);
    rhs[adj_index]->axpy(-1.0, adj);

    // Add the terms from the partial derivative w.r.t. 
    // the design variables
    tacs->evalDVSens(&func, 1, tmp, num_design_vars);
    for ( int i = 0; i < num_design_vars; i++ ){
      dfdx[i] += h*tmp[i];
    }

    // Solve the adjoint
    tacs->assembleJacobian(NULL, mat, alpha, beta, gamma, TRANSPOSE);
    pc->factor();
    ksm->solve(rhs[adj_index], adj);

    // Add the result from the total derivative
    tacs->evalAdjointResProducts(&adj, 1, tmp, num_design_vars);
    for ( int i = 0; i < num_design_vars; i++ ){
      dfdx[i] += h*tmp[i];
    }

    // Zero the entries from the unused 
    rhs[adj_index]->zeroEntries();

    // Add the terms to the remaining adjoint variables
    for ( int ii = 1; (ii < nbdf || ii < nbddf); ii++ ){
      int rhs_index = (k - ii) % num_rhs;
      double beta = 0.0, gamma = 0.0;
      if (ii < nbdf){
        beta = -bdf_coeff[ii]/h;
      }
      if (ii < nbddf){
        gamma = -bddf_coeff[ii]/(h*h);
      }
      tacs->addJacobianVecProduct(1.0, 0.0, beta, gamma,
                                  adj, rhs[rhs_index], TRANSPOSE);
    }
  }
}

/*
  March backward in time and solve for the adjoint variables and find
  the total derivatives
*/
void TACSBDFIntegrator::marchBackwards( ) {
  // Adjoint variables for each function of interest
  BVec **psi = new BVec*[ num_funcs ];
  for ( int n = 0; n < num_funcs; n++ ){
    psi[n] = tacs->createVec();
    psi[n]->incref();
  }
  
  // Right hand sides for adjoint linear-system
  BVec **rhs = new BVec*[ num_funcs*num_adjoint_rhs ];
  for ( int n = 0; n < num_funcs*num_adjoint_rhs; n++ ){
    rhs[n] = tacs->createVec();
    rhs[n]->incref();
  }
 
  // March backwards in time (initial condition not evaluated)
  for ( int k = num_time_steps-1; k >=1 ; k-- ){
    // Get the BDF coefficients at this time step
    this->get2ndBDFCoeff(k, bdf_coeff, &nbdf, bddf_coeff, &nbddf, max_bdf_order);

    // Determine the linearization coefficients for Jacobian Assembly
    double gamma = bddf_coeff[0]/(h*h);
    double beta  = bdf_coeff[0]/h;
    double alpha = 1.0;

    // Set the stages
    this->setTACSStates(q[k], qdot[k], qddot[k]);
    
    // Find the adjoint index
    int adj_index = k % num_adjoint_rhs;

    // Setup the adjoint RHS
    TacsScalar ftmp;
    for ( int n = 0; n < num_funcs; n++ ){
      // Add the contribution to function value from this stage
      tacs->evalFunctions(&funcs[n], 1, &ftmp);
      fvals[n] += h*ftmp;

      // Add up the contribution from function state derivative to RHS
      tacs->evalSVSens(funcs[n], psi[n]);

      // Add the contributions to the current adjoint RHS
      rhs[adj_index*num_funcs+n]->axpy(1.0, psi[n]);
      rhs[adj_index*num_funcs+n]->scale(-1.0);

      // Sanity check on the RHS
      //checkAdjointRHS(rhs[adj_index*num_funcs+n]);
    }
    
    // Setup the Jacobian
    tacs->assembleJacobian(NULL, mat, alpha, beta, gamma, TRANSPOSE);

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
    this->addTotalDerivative(h, psi);

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
  max_newton_iters:  the max number of Newton iterations
*/
TACSDIRKIntegrator::TACSDIRKIntegrator( TACSAssembler * _tacs, 
                                        double _tinit, double _tfinal, 
                                        int _num_steps_per_sec,
                                        int _num_stages ): 
TACSIntegrator(_tacs, _tinit, _tfinal, _num_steps_per_sec){   
  // copy over the variables
  num_stages = _num_stages;
  
  // allocate space for stage state variables
  qS = new BVec*[num_stages*num_time_steps];
  qdotS = new BVec*[num_stages*num_time_steps];
  qddotS = new BVec*[num_stages*num_time_steps];

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
  setupButcherTableau();

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
void TACSDIRKIntegrator::setupButcherTableau(){
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
    setupButcherTableau();
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
  Integrate forward in time and solve for states
*/
TacsScalar TACSDIRKIntegrator::forward( const TacsScalar *x, 
                                        int num_design_vars,
                                        TACSFunction *func ){
  tacs->setDesignVars(x, num_design_vars);
  
  // Get the initial condition
  tacs->getInitConditions(q[0], qdot[0]);
  qddot[0]->zeroEntries();
  
  TacsScalar fval = 0.0;
  for ( int k = 1; k < num_time_steps; k++ ){
    // Set the pointer to the stage values
    BVec **qs = &qS[num_stages*k];
    BVec **qdots = &qdotS[num_stages*k];
    BVec **qddots = &qddotS[num_stages*k];

    q[k]->zeroEntries();
    qdot[k]->zeroEntries();
    qddot[k]->zeroEntries();

    for ( int i = 0; i < num_stages; i++ ){
      // Zero the states (these may not be zero when integrate() is
      // called for the second time)
      qs[i]->zeroEntries();
      qdots[i]->zeroEntries();
      qddots[i]->zeroEntries();
      
      // Compute the stage time
      double t = time[k-1] + C[i]*h;
      
      // Initial guess for qddotS
      if (i == 0){
        qddots[i]->copyValues(qddot[k-1]);
      }
      else {
        qddots[i]->copyValues(qddots[i-1]);
      }
      
      // Compute qdotS
      int idx = getIdx(i);
      qdots[i]->copyValues(qdot[k-1]);
      for ( int j = 0; j <= i; j++ ){
        qdots[i]->axpy(h*A[idx+j], qddots[j]);
      }
      
      // Compute qS
      qs[i]->copyValues(q[k-1]);
      for ( int j = 0; j <= i; j++ ){
        qs[i]->axpy(h*A[idx+j], qdots[j]);
      }
    
      // Determine the coefficients for linearizing the Residual
      double gamma = 1.0;
      double beta  = h*A[idx+i]; 
      double alpha = beta*beta;

      // Solve the nonlinear system of stage equations
      newtonSolve(alpha, beta, gamma, t,
                  qs[i], qdots[i], qddots[i]);

      // Add the contribution to the function value
      TacsScalar ftmp;
      tacs->evalFunctions(&func, 1, &ftmp);
      fval += h*B[i]*ftmp;
    }

    // advance the time
    time[k] = time[k-1] + h;
    
    // advance the position state
    q[k]->copyValues(q[k-1]);
    for ( int j = 0; j < num_stages; j++ ){
      q[k]->axpy(h*B[j], qdots[j]);
    }

    // advance the velocity state
    qdot[k]->copyValues(qdot[k-1]);
    for ( int j = 0; j < num_stages; j++ ){
      qdot[k]->axpy(h*B[j], qddots[j]);
    }
    
    // advance the acceleration state
    qddot[k]->zeroEntries();
    for ( int j = 0; j < num_stages; j++ ){
      qddot[k]->axpy(B[j], qddots[j]);
    }

  }

  return fval;
}
/*
  March backwards in time and find the derivatives
 */
void TACSDIRKIntegrator::reverse( TacsScalar *dfdx, 
                                  int num_design_vars,
                                  TACSFunction *func ){}

/*
  Function that computes the stage values at each time step. This
  function uses the global states and time at the previous time-step
  and sets the stage states tS, qS, qdotS and qdotS.
*/
void TACSDIRKIntegrator::approxStates( int current_step, int current_stage ){
  int k = current_step;
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
  int idx = getIdx(i);
  for ( int j = 0; j <= i; j++ ){
    qdotS[toffset+i]->axpy(h*A[idx+j], qddotS[toffset+j]);
  }

  // Compute qS
  qS[toffset+i]->copyValues(q[k-1]);
  idx = getIdx(i);
  for ( int j = 0; j <= i; j++ ){
    qS[toffset+i]->axpy(h*A[idx+j], qdotS[toffset+j]);
  }    
}

/*
  Start index of the Butcher Tableau A for the supplied stage
*/
int TACSDIRKIntegrator::getIdx( int stageNum ){
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
void TACSDIRKIntegrator::computeTimeStepStates( int current_step, BVec **q, BVec **qdot, BVec **qddot ){
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
void TACSDIRKIntegrator::integrate( ){
  // Get the initial condition
  tacs->getInitConditions(q[0], qdot[0]);

  for ( int k = 1; k < num_time_steps; k++ ){    
    // Find the offset to current state values
    int toffset = k*num_stages;           

    // Compute the stage states qS, qdotS, qddotS based on the DIRK formula
    for ( int i = 0; i < num_stages; i++ ){
      // Compute the stage time
      tS[toffset+i] = time[k-1] + C[i]*h;

      // Approximate the stage states using the DIRK formula
      approxStates(k, i);

      // Determine the coefficients for linearizing the Residual
      int idx = getIdx(i);

      double gamma = 1.0;
      double beta  = h*A[idx+i]; 
      double alpha = beta*h*A[idx+i];

      // Solve the nonlinear system of stage equations starting with the approximated states
      newtonSolve(alpha, beta, gamma, tS[toffset+i], 
		  qS[toffset+i], qdotS[toffset+i], qddotS[toffset+i]);
    }
    
    // Advance the time
    time[k] = time[k-1] + h;

    // Compute the state varialbes at the current time step using the
    // intermediate stage states
    computeTimeStepStates(k, q, qdot, qddot);

    // Write the tecplot output to disk if sought
    if(f5 && f5_write_freq > 0 && k % f5_write_freq == 0){
      // Create a buffer for filename 
      char buffer[128];
      // Format the buffer based on the time step
      getString(buffer, f5_file_fmt, k);      
      // Write the f5 file for this time step
      f5->writeToFile(buffer);
    }
  }
}
/*
  March backward in time and solve for the adjoint variables
*/
void TACSDIRKIntegrator::marchBackwards( ) {
  // Inter-step adjoint variables for each function of interest
  BVec **psi = new BVec*[ num_funcs ];
  BVec **phi = new BVec*[ num_funcs ];

  BVec **psiTmp = new BVec*[ num_funcs ];
  BVec **phiTmp = new BVec*[ num_funcs ];

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
  BVec **lambda = new BVec*[ num_funcs*num_adjoint_rhs ];
  BVec **rhs = new BVec*[ num_funcs*num_adjoint_rhs ];
  BVec **dfdq = new BVec*[ num_funcs*num_adjoint_rhs ];
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
    // March backwards in stage
    int toffset = k*num_stages;

    // Set the pointer to the stage values
    BVec **qs = &qS[toffset];
    BVec **qdots = &qdotS[toffset];
    BVec **qddots = &qddotS[toffset];
    
    for ( int i = num_stages-1; i >= 0; i-- ) {

      //-------------------------------------------------------------//
      // Solve for the adjoint variables 
      //-------------------------------------------------------------//

      // Compute the time
      double t = time[k-1] + C[i]*h;

      // Get the starting index for the corresponding Butcher's tableau row
      int idx = getIdx(i);

      // Determine the coefficients for linearizing the Residual
      double gamma = B[i];
      double beta  = B[i]*h*A[idx+i]; 
      double alpha = beta*h*A[idx+i];

      // Set the stages
      this->setTACSStates(qs[i], qdots[i], qddots[i]);

      //--------------------------------------------------------------//
      // Assemble the right hand side
      //--------------------------------------------------------------//
      
      // Add function contributions
      TacsScalar ftmp;
      for ( int n = 0; n < num_funcs; n++ ){
        // Add up the function contribution (dfdq)
	tacs->evalFunctions(&funcs[n], 1, &ftmp);
        fvals[n] += h*B[i]*ftmp;

        // Add up the contribution from its derivative to RHS (drdq.lam)
	tacs->evalSVSens(funcs[n], dfdq[i*num_funcs+n]);
	rhs[i*num_funcs+n]->axpy(alpha, dfdq[i*num_funcs+n]);

        // Add up the contribution from PSI to this RHS
        rhs[i*num_funcs+n]->axpy(B[i], phi[n]);

        // Add up the contribution from PSI to this RHS
        TacsScalar scale = 0.0;
        for ( int jj = i; jj < num_stages; jj++ ){
          idx = getIdx(jj);
          scale += B[jj]*h*A[idx+i];
        }
        rhs[i*num_funcs+n]->axpy(scale, psi[n]);
      
        // Negate the RHS
	rhs[i*num_funcs+n]->scale(-1.0);

        // Sanity check on the RHS
        //checkAdjointRHS(rhs[i*num_funcs+n]);
      }
     
      // Setup the Jacobian
      tacs->assembleJacobian(NULL, mat, alpha, beta, gamma, TRANSPOSE);

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
      this->addTotalDerivative(h*B[i], &lambda[i*num_funcs]);
      
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
  Evaluate the time average of functions: F = \sum_{k=0}^N h f(q,t)
*/
void TACSBDFIntegrator::evalTimeAvgFunctions( TACSFunction **funcs, 
					      int numFuncs, 
					      TacsScalar *funcVals) {
  memset(funcVals, 0, numFuncs*sizeof(TacsScalar));
  
  TacsScalar *ftmp = new TacsScalar[numFuncs];  
  
  // Loop over time steps
  for ( int k = 1; k < num_time_steps; k++ ) {
    // Set the states into TACS
    setTACSStates(q[k], qdot[k], qddot[k]);

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
  Evaluate the time average of functions: F = \sum_{k=0}^N h_k \sum_{j=0}^s b_j f(q_{k,j},t)
*/
void TACSDIRKIntegrator::evalTimeAvgFunctions( TACSFunction **funcs, 
					       int numFuncs, 
					       TacsScalar *funcVals) {
  memset(funcVals, 0, numFuncs*sizeof(TacsScalar));

  TacsScalar *ftmp = new TacsScalar[numFuncs];  

  // Loop over time steps
  for ( int k = 1; k < num_time_steps; k++ ) {
    int toffset = k*num_stages;

    // Loop over all stages
    for ( int j = 0; j < num_stages; j++ ){

      // Set the states into TACS
      setTACSStates(qS[toffset+j], qdotS[toffset+j], qddotS[toffset+j]);
      
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
