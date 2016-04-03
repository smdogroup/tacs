#include "TACSIntegrator.h"
#include "tacslapack.h"

/*
  Abstract class for integration schemes to extend
*/
TacsIntegrator::TacsIntegrator() {
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
TacsDIRKIntegrator::TacsDIRKIntegrator(int _numStages, int _numVars,  
				       double _tInit, double _tFinal, 
				       int _numStepsPerSec,
				       int _maxNewtonIters) 
: TacsIntegrator() {

  // copy over the input parameters
  numStages = _numStages;
  numVars = _numVars;
  tInit = _tInit;
  tFinal = _tFinal;
  numStepsPerSec = _numStepsPerSec;

  // create an instance of a nonlinear solver
  nonlinearSolver =  new TacsNewtonSolver(_numVars, _maxNewtonIters, 1e-30, 1.0e-8);
  nonlinearSolver->incref();
  
  // compute the step size
  h = 1.0/double(numStepsPerSec);

  // compute the total number of time steps
  numSteps = int(double(numStepsPerSec)*(tFinal-tInit)) + 1;

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

  nonlinearSolver->decref();

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
  Function that wraps the integration logic.
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







