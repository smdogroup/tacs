#include <math.h>
#include "TACSIntegrator.h"

/*
  Base class constructor for integration schemes. This base class
  contains common methods and variables pertaining to the integration
  schemes used in TACS.

  Input:
  tInit: the initial time
  tfinal: the final time
  num_steps_per_sec: the number of steps to take for each second
*/
TacsIntegrator::TacsIntegrator( TACSAssembler * _tacs,
                                double _tinit, double _tfinal, 
                                int _num_steps_per_sec ){
  // copy over the input parameters
  tacs = _tacs;
  tacs->incref();

  tinit = _tinit;
  tfinal = _tfinal;
  num_steps_per_sec = _num_steps_per_sec;

  // compute the step size
  h = 1.0/double(num_steps_per_sec);

  // compute the total number of time steps
  num_time_steps = int(double(num_steps_per_sec)*(tfinal-tinit)) + 1;

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
  for (int k = 0; k < num_time_steps; k++) {
    q[k] = tacs->createVec(); 
    q[k]->incref(); 

    qdot[k] = tacs->createVec(); 
    qdot[k]->incref(); 

    qddot[k] = tacs->createVec(); 
    qddot[k]->incref(); 
  }

  // store time
  time = new double[num_time_steps];
  memset(time, 0, num_time_steps*sizeof(double));

  //------------------------------------------------------------------//
  //                     Newton's method                              //
  //------------------------------------------------------------------//

  // Default parameters for Newton solve
  max_newton_iters = 25;
  atol = 1.0e-30;
  rtol = 1.0e-12;

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
}

/*
  Destructor
*/
TacsIntegrator::~TacsIntegrator(){
  // Dereference TACS
  tacs->decref();

  // Dereference position, velocity and acceleration states
  for (int k = 0; k < num_time_steps; k++) {
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

  Drives the residual R(t,q,qdot,qddot) to zero using Newton's method
  with KSM solver

  Input: 
  The guessed (initial) state variable values q, qdot, qddot are
  supplied
  
  Output: q, qdot, qddot updated iteratively until the corresponding
  residual R = 0
  
  alpha: multiplier for derivative of Residual wrt to q
  beta : multiplier for derivative of Residual wrt to qdot
  gamma: multiplier for derivative of Residual wrt to qddot

*/
void TacsIntegrator::newtonSolve( double alpha, double beta, double gamma,
                                  double t, BVec *q, BVec *qdot, 
                                  BVec *qddot ){
  // Initialize the norms
  TacsScalar init_norm = 0.0;
  TacsScalar norm = 0.0;

  // Iterate until max iters or R <= tol
  for ( int n = 0; n < max_newton_iters; n++ ){
    // Set the supplied initial input states into TACS
    tacs->setVariables(q);
    tacs->setDotVariables(qdot);
    tacs->setDDotVariables(qddot);

    // Assemble the Jacobian matrix once in five newton iterations
    if (1){ // n % 5 == 0){
      tacs->assembleJacobian(res, mat, alpha, beta, gamma, NORMAL);
    } 
    else {
      tacs->assembleRes(res);
    }    

    // Compute the L2-norm of the residual
    norm = res->norm();
    
    // Write a summary
    if(print_level > 0) {
      printf("# Newton iters=%d, |R|=%e \n", n, RealPart(norm));
    }
    // Record the residual norm at the first Newton iteration
    if (n == 0){
      init_norm = norm;
    }
           
    // Check if the norm of the residuals is a NaN
    if (norm != norm){ 
      /* fprintf(stderr,  */
      /*         "Newton iteration %d, failed with NaN residual norm\n", n); */
      break;
    }
    
    // Check if the Newton convergence tolerance is satisfied
    if (norm < rtol*init_norm || norm < atol){
      break;
    }

    // Factor the preconditioner
    if (1){ //n % 5 == 0){
      pc->factor();
    }
   
    // Solve for update using KSM
    ksm->solve(res, update);

    // Update the state variables using the solution
    qddot->axpy(-gamma, update);
    qdot->axpy(-beta, update);
    q->axpy(-alpha, update);

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
void TacsIntegrator::writeSolution( const char *filename ){
  // Temporary variables to access the states at each time
  TacsScalar *qvals, *qdotvals, *qddotvals;

  // Figure out the number of variables
  int num_vars;
  q[0]->getSize(&num_vars);

  // Open a new file
  FILE *fp = fopen(filename, "w");
 
  for ( int k = 0; k < num_time_steps; k++ ){    
    // Copy over the state values from BVec
    q[k]->getArray(&qvals);
    qdot[k]->getArray(&qdotvals);
    qddot[k]->getArray(&qddotvals);
  
    // Write the time
    fprintf(fp, "%e ", time[k]);

    // Write the states (q, qdot, qddot) to file
    for ( int j = 0; j < num_vars; j++ ){
      fprintf(fp, "%e %e %e ", RealPart(qvals[j]), 
              RealPart(qdotvals[j]), RealPart(qddotvals[j]));
    }
    fprintf(fp, "\n");
  }

  // Close the file
  fclose(fp);
}

/*
  Constructor for BDF Integration scheme

  Input:
  tinit: the initial time
  tfinal: the final time
  num_steps_per_sec: the number of steps to take for each second
*/

TacsBDFIntegrator:: TacsBDFIntegrator( TACSAssembler * _tacs, 
                                       double _tinit, double _tfinal, 
                                       int _num_steps_per_sec, 
                                       int _max_bdf_order ):
TacsIntegrator(_tacs, _tinit,  _tfinal,  _num_steps_per_sec){		
  // copy over the variables
  max_bdf_order = _max_bdf_order;

  // Truncate the maximum order to 3rd order
  max_bdf_order = (max_bdf_order <= 3 ? 
		   max_bdf_order : 3);
}

/*
  Destructor for TACSBDFIntegrator
*/
TacsBDFIntegrator::~TacsBDFIntegrator(){}

/*
  Integration logic of BDF. Use this function to march in time. The
  solution over time is set into the class variables q, qdot and qddot
  and time.
*/
void TacsBDFIntegrator::integrate(){
  current_time_step = 0;

  // Initial condition
  tacs->getInitConditions(q[0], qdot[0]);

  for ( int k = 1; k < num_time_steps; k++ ){
    current_time_step++;

    // Approximate states and their derivatives using BDF formula
    approxStates(q, qdot, qddot);
    
    // Determine the coefficients for Jacobian Assembly
    int gamma = bddf_coeff[0]/h/h;
    int beta = bdf_coeff[0]/h;
    int alpha = 1.0;

    // Solve the nonlinear system of equations
    newtonSolve(alpha, beta, gamma, time[k], q[k], qdot[k], qddot[k]);
    
    // Advance time (states are already advanced at the end of Newton solve)
    time[k] = time[k-1] + h;
  }
}

/*
  Approximate states (q, qdot, qddot) at the current time step using
  the BDF coefficients and previous time step values of the states q,
  qdot and qddot.
  
  Input:
  pointers to the global states q, qdot, qddot
*/
void TacsBDFIntegrator::approxStates( BVec **q, BVec **qdot, BVec **qddot ){
  int k = current_time_step;
  
  // get the BDF coefficients
  int nbdf, nbddf;
  get2ndBDFCoeff(k, bdf_coeff, &nbdf, bddf_coeff, &nbddf, max_bdf_order);

  // copy over the q values from previous time step into the current q
  q[k]->copyValues(q[k-1]);

  // approximate qdot using BDF formula
  for ( int i = 0; i < nbdf; i++ ){
    double scale = bdf_coeff[i]/h;
    qdot[k]->axpy(scale, q[k-i]);
  }

  // approximate qddot using BDF formula
  for ( int i = 0; i < nbddf; i++ ){
    double scale = bddf_coeff[i]/h/h;
    qddot[k]->axpy(scale, q[k-i]);
  }

  // If required, add the contribution to the second derivative
  // from the initial values of the first derivatives
  if (k == nbdf-1){
    double scale = bdf_coeff[nbdf-1]/h;
    qddot[k]->axpy(scale, qdot[k-1]);
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
void TacsBDFIntegrator::get2ndBDFCoeff( const int k, 
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
int TacsBDFIntegrator::getBDFCoeff( double bdf[], int order ){
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
  Constructor for TacsDIRKIntegrator

  Input:

  num_stages:        the number of Runge-Kutta stages
  tinit:             the initial time
  tfinal:            the final time
  num_steps_per_sec: the number of steps to take for each second
  max_newton_iters:  the max number of Newton iterations
*/
TacsDIRKIntegrator::TacsDIRKIntegrator( TACSAssembler * _tacs, 
                                        double _tinit, double _tfinal, 
                                        int _num_steps_per_sec,
                                        int _num_stages ): 
TacsIntegrator(_tacs, _tinit,  _tfinal,  _num_steps_per_sec){   
  // copy over the variables
  num_stages = _num_stages;
  
  // allocate space for stage state variables
  qS = new BVec*[num_stages];
  qdotS = new BVec*[num_stages];
  qddotS = new BVec*[num_stages];

  // store stage time
  tS = new double[num_stages];
  memset(tS, 0, num_stages*sizeof(double));
  
  // create state vectors for TACS during each timestep
  for ( int k = 0; k < num_stages; k++ ){
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
}

/*
  Destructor for TacsDIRKIntegrator
*/
TacsDIRKIntegrator::~TacsDIRKIntegrator(){
  // Cleanup Butcher's Tableau
  delete [] A;
  delete [] B;
  delete [] C;

  // Cleanup stage stages
  for ( int i = 0; i < num_stages; i++ ){
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
void TacsDIRKIntegrator::setupButcherTableau(){
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
void TacsDIRKIntegrator::checkButcherTableau(){
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
  Integration logic of DIRK. Use this function to march in time. The
  solution over time is set into the class variables q, qdot and qddot
  and time.
*/
void TacsDIRKIntegrator::integrate(){
  current_time_step = 0;
  
  // Initial condition
  tacs->getInitConditions(q[0], qdot[0]);

  for ( int k = 1; k < num_time_steps; k++ ){
    current_time_step++;
       
    // Compute the stage states qS, qdotS, qddotS based on the DIRK formula
    computeStageValues();
    
    // Advance the global state to the next time
    timeMarch(time, q, qdot, qddot);
  }
}

/*
  Function that computes the stage values at each time step. This
  function uses the global states and time at the previous time-step
  and sets the stage states tS, qS, qdotS and qdotS.
*/
void TacsDIRKIntegrator::computeStageValues(){
  int k = current_time_step;

  // Set the stage values tS, qS, qdotS, qddotS to zero before
  // evaluting them at the current time step
  resetStageValues();

  for ( int i = 0; i < num_stages; i++ ){
    // Compute the stage time
    tS[i] = time[k-1] + C[i]*h;

    // Initial guess for qddotS
    if (i == 0){
      qddotS[i]->copyValues(qddot[k-1]);
    }
    else {
      qddotS[i]->copyValues(qddotS[i-1]);
    }

    // Compute qdotS
    int idx1 = getIdx(i);
    for ( int j = 0; j <= i; j++ ){
      qdotS[i]->axpy(h*A[idx1], qddotS[j]);
      idx1++;
    }
    qdotS[i]->axpy(1.0, qdot[k-1]);

    // Compute qS
    int idx2 = getIdx(i);
    for ( int j = 0; j <= i; j++ ){
      qS[i]->axpy(h*A[idx2], qdotS[j]);
      idx2++;
    }
    qS[i]->axpy(1.0, q[k-1]);
    
    // Determine the coefficients for linearizing the Residual
    double gamma = 1.0;
    double beta  = h*A[0]; 
    double alpha = h*A[0]*h*A[0];

    // Solve the nonlinear system of stage equations
    newtonSolve(alpha, beta, gamma, tS[i], qS[i], qdotS[i], qddotS[i]);
  }
}

/*
  Start index of the Butcher Tableau A for the supplied stage
*/
int TacsDIRKIntegrator::getIdx( int stageNum ){
  return stageNum*(stageNum+1)/2;
}

/*
  Function that advances the global states q, qdot, qddot and time to
  next time step
  
  Input:
  The pointers to the global state variables q, qdot, qddot
  
  Output:
  Updated global state variables q, qdot, qddot at the time step
*/
void TacsDIRKIntegrator::timeMarch( double *time, 
				    BVec **q, BVec **qdot, BVec **qddot ){
  int k = current_time_step;
  
  // advance the time
  time[k] = time[k-1] + h;
  
  // advance the position state
  for ( int j = 0; j < num_stages; j++ ){
    q[k]->axpy(h*B[j], qdotS[j]);
  }
  q[k]->axpy(1.0, q[k-1]);

  // advance the velocity state
  for ( int j = 0; j < num_stages; j++ ){
    qdot[k]->axpy(h*B[j], qddotS[j]);
  }
  qdot[k]->axpy(1.0, qdot[k-1]);

  // advance the acceleration state
  for ( int j = 0; j < num_stages; j++ ){
    qddot[k]->axpy(B[j], qddotS[j]);
  }
}

/*
  Function that sets the arrays used during previous time marching
  step to zero, to make way for the next time step
*/
void TacsDIRKIntegrator::resetStageValues(){
  memset(tS, 0, num_stages*sizeof(double));

  for (int i=0; i < num_stages; i++){
    qS[i]->zeroEntries();
    qdotS[i]->zeroEntries();
    qddotS[i]->zeroEntries();
  }
}
