#include "TACSIntegrator.h"
#include "tacslapack.h"

/*
  Abstract class for integration schemes to extend

  Input:

  numVars : the number of varibles/degrees of freedom in the system
  tInit: the initial time
  tFinal: the final time
  numStepsPerSec: the number of steps to take for each second
  maxNewtonIters: the max number of Newton iterations
  
*/
TacsIntegrator::TacsIntegrator(int _numVars, double _tInit, 
			       double _tFinal, int _numStepsPerSec, 
			       int _maxNewtonIters, double _atol, double _rtol) {
  
  // copy over the input parameters
  numVars = _numVars;
  tInit = _tInit;
  tFinal = _tFinal;
  numStepsPerSec = _numStepsPerSec;
  maxNewtonIters = _maxNewtonIters;

  // compute the step size
  h = 1.0/double(numStepsPerSec);

  // compute the total number of time steps
  numSteps = int(double(numStepsPerSec)*(tFinal-tInit)) + 1;

  // create an instance of a nonlinear solver
  nonlinearSolver =  new TacsNewtonSolver(_numVars, _maxNewtonIters, _atol, _rtol);
  nonlinearSolver->incref();

}

/*
  Base class destructor for the integration schemes
*/
TacsIntegrator::~TacsIntegrator(){

  nonlinearSolver->decref();

}

/*
  Constructor for BDF Integration scheme

  Input:

  numVars : the number of varibles/degrees of freedom in the system
  tInit: the initial time
  tFinal: the final time
  numStepsPerSec: the number of steps to take for each second
  maxNewtonIters: the max number of Newton iterations
  
*/
TacsBDFIntegrator:: TacsBDFIntegrator(int _numVars, double _tInit, 
				      double _tFinal, int _numStepsPerSec, 
				      int _maxNewtonIters, double _atol, double _rtol,
				      int _maxBDFOrder) 
: TacsIntegrator(_numVars, _tInit,  _tFinal,  _numStepsPerSec, _maxNewtonIters, _atol, _rtol) {
  
  // copy over the variables
  maxBDFOrder = _maxBDFOrder;

  // Truncate the maximum order to 3rd order
  maxBDFOrder = (maxBDFOrder <= 3 ? 
		 maxBDFOrder : 3);

  // allocate space for the residuals
  res = new TacsScalar[numVars];
  memset(res, 0, numVars*sizeof(TacsScalar));
    
  // allocate space for the jacobian
  D = new TacsScalar[numVars*numVars];
  memset(D, 0, numVars*numVars*sizeof(TacsScalar));
  
}

/*
  Destructor for TACSBDFIntegrator
*/
TacsBDFIntegrator::~TacsBDFIntegrator(){

  delete [] res;
  delete [] D;

}

/*
  Set the coefficients for Jacobian linearization and evaluation
*/
void TacsBDFIntegrator::setCoeffs(double _alpha, double _beta, double _gamma){
  alpha = _alpha;
  beta  = _beta;
  gamma = _gamma;
}

/*
  Integration logic of BDF
  
  Input:
  
  time: pointer to global time
  q: pointer to the position states (initial value should be set)
  qdot: pointer to the velocity states (initial value should be set)
  qddot: pointer to the acceleration states

*/
void TacsBDFIntegrator::integrate(TacsScalar *time, 
				  TacsScalar *q, 
				  TacsScalar *qdot, 
				  TacsScalar *qddot){
  
  currentTimeStep = 0;
  
  for (int k = 1; k < numSteps; k++) {
    
    currentTimeStep ++;

    // approximate states and their derivatives using BDF formula
    approxStates(q, qdot, qddot);
    
    // get the coeffcients for residual linearization during the nonlinear solve
    setCoeffs(bddf_coeff[0]/h/h, bdf_coeff[0]/h, 1.0);
    
    // solve the nonlinear system of equations for q
    nonlinearSolver->solve(alpha, beta, gamma,
			   time[k], &q[k*numVars], &qdot[k*numVars], &qddot[k*numVars]);
    
    time[k] = time[k-1] + h;
    
  }
  
}

/*
  Approximate states at the current time step using the BDF
  coefficients and previous time step values of the states

  Input:   
  pointers to the global states q, qdot, qddot
  
*/
void TacsBDFIntegrator::approxStates(TacsScalar *q, 
				     TacsScalar *qdot, 
				     TacsScalar *qddot){
  
  int k = currentTimeStep;
  int idx =  k*numVars;
  
  // get the BDF coefficients
  int nbdf, nbddf;
  get2ndBDFCoeff(k, bdf_coeff, &nbdf, bddf_coeff, &nbddf, 
		 maxBDFOrder);
  
  // Copy the values of q from the previous time step
  memcpy(&q[idx], &q[numVars*(k-1)], numVars*sizeof(TacsScalar));
  
  // approximate qdot using BDF formula
  for ( int i = 0; i < nbdf; i++ ){
    const TacsScalar *qi = &q[numVars*(k-i)];
    double scale = bdf_coeff[i]/h;
    for ( int j = 0; j < numVars; j++ ){
      qdot[idx+j] += scale*qi[j];
    }
  }

  // approximate qddot using BDF formula
  for ( int i = 0; i < nbddf; i++ ){
    const TacsScalar *qi = &q[numVars*(k-i)];
    double scale = bddf_coeff[i]/h/h;
    for ( int j = 0; j < numVars; j++ ){
      qddot[idx+j] += scale*qi[j];
    }
  }

  // If required, add the contribution to the second derivative
  // from the initial values of the first derivatives
  if (k == nbdf-1){
    double scale = bdf_coeff[nbdf-1]/h;
    for ( int j = 0; j < numVars; j++ ){
      qddot[idx+j] += scale*qdot[(k-1)*numVars+j];
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
int TacsBDFIntegrator::getBDFCoeff( double bdf[], 
				    int order ){
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

  numStages: the number of Runge-Kutta stages
  numVars : the number of varibles/degrees of freedom in the system
  tInit: the initial time
  tFinal: the final time
  numStepsPerSec: the number of steps to take for each second
  maxNewtonIters: the max number of Newton iterations

*/
TacsDIRKIntegrator::TacsDIRKIntegrator(int _numVars,  
				       double _tInit, double _tFinal, 
				       int _numStepsPerSec,
				       int _maxNewtonIters, double _atol, double _rtol,
				       int _numStages) 
: TacsIntegrator(_numVars, _tInit,  _tFinal,  _numStepsPerSec, _maxNewtonIters, _atol, _rtol) {
  
  // copy over the variables
  numStages = _numStages;
  
  // allocate space for stage variables
  tS = new TacsScalar[numStages];
  qS = new TacsScalar[numStages*numVars];
  qdotS = new TacsScalar[numStages*numVars];
  qddotS = new TacsScalar[numStages*numVars];

  // initialize states to zero
  memset(tS, 0, numStages*sizeof(TacsScalar));
  memset(qS, 0, numStages*numVars*sizeof(TacsScalar));
  memset(qdotS, 0, numStages*numVars*sizeof(TacsScalar));
  memset(qddotS, 0, numStages*numVars*sizeof(TacsScalar));
  
  // allocate space for Butcher tableau
  A = new double[numStages*(numStages+1)/2];
  B = new double[numStages];
  C = new double[numStages];

  // set the Butcher tableau entries to zero
  memset(A, 0, numStages*(numStages+1)/2*sizeof(double));
  memset(B, 0, numStages*sizeof(double));
  memset(C, 0, numStages*sizeof(double));

  // now add entries into the Butcher tableau
  setupButcherTableau();

  // determine the integration coeffs
  alpha = 1.0;
  beta  = h*A[0]; 
  gamma = h*A[0]*h*A[0];

}

/*
  Destructor for TACSSIRKIntegrator
*/
TacsDIRKIntegrator::~TacsDIRKIntegrator(){
  
  delete [] A;
  delete [] B;
  delete [] C;

  delete [] tS;
  delete [] qS;
  delete [] qdotS;
  delete [] qddotS;

}

/*
  Function that adds the entries into Butcher tableau
*/
void TacsDIRKIntegrator::setupButcherTableau(){
  
  if (numStages == 1) {

    // Implicit mid-point rule (A-stable)

    A[0] = 0.5;
    B[0] = 1.0;
    C[0] = 0.5;

    order = 2;
    
  } else if (numStages == 2) {

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
    
  } else if (numStages == 3) {

    // Crouzeix formula (A-stable)
    double PI  = 22.0/7.0;
    double alpha = 2.0*cos(PI/18.0)/sqrt(3.0);
    
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

  } else {

    fprintf(stderr, "ERROR: Invalid number of stages %d\n", numStages);
    exit(-1);
    
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
  for (int i = 0; i < numStages; i++) {
    
    tmp = 0.0;
    for (int j = 0; j <= i; j++) {
      idx ++;
      tmp += A[idx];
    }

    // check the difference
    if (fabs(C[i] - tmp) >= 1.0e-6) {
      fprintf(stderr, "WARNING: Sum A[%d,:] != C[%d] i.e. %f != %f \n", 
	      i, i, C[i], tmp);
    }
   
  }
  
  // Check #2: sum(B) = 1.0
  tmp = 0.0;
  for (int i = 0; i < numStages; i++) {
    tmp += B[i];
  }
  if (fabs(1.0 - tmp) >= 1.0e-6) {
    fprintf(stderr, "WARNING: Sum B != 1.0 \n");
  }
  
}

/*
  Integration logic of DIRK
  
  Input:
  
  time: pointer to global time
  q: pointer to the position states (initial value should be set)
  qdot: pointer to the velocity states (initial value should be set)
  qddot: pointer to the acceleration states

*/
void TacsDIRKIntegrator::integrate(TacsScalar *time, 
				   TacsScalar *q, 
				   TacsScalar *qdot, 
				   TacsScalar *qddot){
  
  currentTimeStep = 0;

  for (int k = 1; k < numSteps; k++) {

    currentTimeStep ++;
    
    // set the stage values tS, qS, qdotS, qddotS to zero
    resetStageValues();
    
    // pass in the global states at previous time step and compute the
    // stage states
    computeStageValues(time[k-1], 
		       &q[(k-1)*numVars], 
		       &qdot[(k-1)*numVars], 
		       &qddot[(k-1)*numVars]);
    
    // advance the global state to the next time
    timeMarch(k, time, q, qdot, qddot);

  }
  
}

/*
  Function that computes the stage values at each time step. 

  Input: 
  The global states and time at the previous time-step

  Output:
  The stage states tS, qS, qdotS and qdotS
*/
void TacsDIRKIntegrator::computeStageValues(TacsScalar tk, 
					    TacsScalar *qk, 
					    TacsScalar *qdotk, 
					    TacsScalar *qddotk){

  // global counter to the current stage
  int currentStage = 0;
  
  // stage offset for q values
  int stageOffCtr= 0;

  for (int i = 0; i < numStages; i++) {

    // compute the stage time
    tS[i] = tk + C[i]*h;

    // initial guess for qddotS
    for (int n = 0; n < numVars; n++) {
      qddotS[stageOffCtr+n] = 1.0;
    }

    // compute qdotS
    for (int n = 0; n < numVars; n++) {  
      int ctr2 = 0;
      TacsScalar tmp = 0.0;
      int idx1 = getIdx(i);
      for (int j = 0; j <= i; j++) {
	tmp += A[idx1]*qddotS[ctr2 + n];
	ctr2 +=  numVars;
	idx1++;
      }
      qdotS[stageOffCtr+n] = qdotk[n] + h*tmp;
    }
    
    // compute qS
    for (int n = 0; n < numVars; n++) {  
      int ctr2 = 0;
      TacsScalar tmp = 0.0;
      int idx2 = getIdx(i);
      for (int j = 0; j <= i; j++) {
	tmp += A[idx2]*qdotS[ctr2 + n];
	ctr2 +=  numVars;
	idx2++;
      }
      qS[stageOffCtr+n] = qk[n] + h*tmp;
    }

    // solve the nonlinear system of stage equations for qddotS
    nonlinearSolver->solve(alpha, beta, gamma,
			   tS[i], &qS[stageOffCtr], &qdotS[stageOffCtr], &qddotS[stageOffCtr]);
    
    // recompute the offsets for next stage
    stageOffCtr +=  numVars;

    // increament the global stage counter
    currentStage ++;
   
  }
  
}

/*
  Start index of the Butcher Tableau A for the supplied stage
*/
int TacsDIRKIntegrator::getIdx(int stageNum){
  return stageNum*(stageNum+1)/2;
}

/*
  Function that advances the states and time to next time step
  
  Input:
  The pointers to the global state variables q, qdot, qddot
  
  Output:
  Updated global state variables q, qdot, qddot at the time step
  
*/
void TacsDIRKIntegrator::timeMarch(int k, TacsScalar *time, TacsScalar *q, 
				   TacsScalar *qdot, TacsScalar *qddot){
  
  // advance the time
  time[k] = time[k-1] + h;
  
  // advance the state
  for (int i = 0; i < numVars; i++) {
    TacsScalar tmp = 0.0;
    int ctr = 0;
    for (int j = 0; j < numStages; j++) {
      tmp += B[j]*qdotS[ctr+i];
      ctr += numVars;
    }
    q[k*numVars+i] = q[(k-1)*numVars+i] + h*tmp;
  }

  // advance the velocity state
  for (int i = 0; i < numVars; i++) {
    TacsScalar tmp = 0.0;
    int ctr = 0;
    for (int j = 0; j < numStages; j++) {
      tmp += B[j]*qddotS[ctr+i];
      ctr += numVars;
    }
    qdot[k*numVars+i] = qdot[(k-1)*numVars+i] + h*tmp;
  }

  // advance the acceleration state
  for (int i = 0; i < numVars; i++) {
    TacsScalar tmp = 0.0;
    int ctr = 0;
    for (int j = 0; j < numStages; j++) {
      tmp += B[j]*qddotS[ctr+i];
      ctr += numVars;
    }
    qddot[k*numVars+i] = qddot[(k-1)*numVars+i] + tmp;
  }

  
}

/*
  Function that sets the arrays used during previous time marching
  step to zero, to make way for the next time step
*/
void TacsDIRKIntegrator::resetStageValues(){

  memset(tS, 0, numStages*sizeof(TacsScalar));

  memset(qS, 0, numStages*numVars*sizeof(TacsScalar));
  memset(qdotS, 0, numStages*numVars*sizeof(TacsScalar));
  memset(qddotS, 0, numStages*numVars*sizeof(TacsScalar));

}







