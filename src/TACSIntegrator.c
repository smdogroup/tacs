#include "TACSIntegrator.h"
#include "tacslapack.h"

/*
  Constructor for TacsDIRKIntegrator

  Input:

  numStages: the number of Runge-Kutta stages
  numVars : the number of varibles/degrees of freedom in the system
  tInit: the initial time
  tFinal: the final time
  numStepsPerSec: the number of steps to take for each second

*/
TacsDIRKIntegrator::TacsDIRKIntegrator(int _numStages, int _numVars,  
				       double _tInit, double _tFinal, 
				       int _numStepsPerSec,
				       int _max_newton_iters){

  // copy over the input parameters
  numStages = _numStages;
  numVars = _numVars;
  tInit = _tInit;
  tFinal = _tFinal;
  numStepsPerSec = _numStepsPerSec;
  max_newton_iters = _max_newton_iters;

  // sanity check the input params

  // compute the step size
  h = 1.0/numStepsPerSec;

  // compute the total number of time steps
  numSteps = int(numStepsPerSec*(tFinal-tInit)) + 1;

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
  Function that wraps the integration logic.
*/
void TacsDIRKIntegrator::integrate(TacsScalar *time, 
				   TacsScalar *q, 
				   TacsScalar *qdot, 
				   TacsScalar *qddot){
  
  for (int k = 1; k <= numSteps; k++) {
    
    // pass in the global states at previous time step
    computeStageValues(time[k-1], &q[(k-1)*numVars], &qdot[(k-1)*numVars], 
			 &qddot[(k-1)*numVars]);
    
    // advance the state to the current step
    timeMarch(k, time, q, qdot, qddot);

    // set the stage values to zero
    resetStageValues();
    
  }
  
}

/*
  Function that computes the values at each time step. 

  Input: 
  The global states and time for at the prvious time-step

  Output:
  The state states tS, qS, qdotS and qdotS
*/
void TacsDIRKIntegrator::computeStageValues(TacsScalar tk, 
					      TacsScalar *qk, 
					      TacsScalar *qdotk, 
					      TacsScalar *qddotk){

  currentStage = 0;

  // stage offset for q values
  stageOffCtr= 0;

  for (int i = 0; i < numStages; i++) {

    // compute the stage time
    tS[i] = tk + C[i]*h;

    // initial guess for qddotS
    for (int n = 0; n < numVars; n++) {
      qddotS[stageOffCtr+n] = 1.0;
    }

    // prepare qdotS
    int stageOffCtr2 = 0;
    for (int n = 0; n < numVars; n++) {  
      TacsScalar tmp = 0.0;
      for (int j = 0; j <= i; j++) {
	int idx1 = getIdx(j);
	tmp += A[idx1]*qddotS[stageOffCtr2 + n];
	stageOffCtr2 +=  numVars;
	idx1++;
      }
      qdotS[stageOffCtr+n] = qdotk[n] + h*tmp;
    }
    
    // prepare qS
    stageOffCtr2 = 0;
    for (int n = 0; n < numVars; n++) {  
      TacsScalar tmp = 0.0;
      for (int j = 0; j <= i; j++) {
	int idx2 = getIdx(j);
	tmp += A[idx2]*qdotS[stageOffCtr2 + n];
	stageOffCtr2 +=  numVars;
	idx2++;
      }
      qS[stageOffCtr+n] = qk[n] + h*tmp;
    }

    // solve the nonlinear system of stage equations for qddotS
    nonlinearSolve(tS[i], &qS[stageOffCtr], &qdotS[stageOffCtr], &qddotS[stageOffCtr]);

    // recompute the offsets
    stageOffCtr +=  numVars;

    // increament the global stage counter
    currentStage ++;
    
  }
  
}

/*
  Start index for the supplied stage  
*/
int TacsDIRKIntegrator::getIdx(int stageNum){
  
  if (stageNum == 0) {
    return 0;
  } else if (stageNum == 1) {
    return 1;
  } else if (stageNum == 2) {
    return 3;
  }  else if (stageNum == 3) {
    return 6;
  }

}

/*
  Solve the nonlinear stage equations using Newton's method.

  Input: 
  
  The stage-state variable values, where qddotS is a guessed
  solution -- qdotS, qS are computed from qddotS 
  
  Output: 
  
  qS, qdotS, qddotS solved iteratively until the corresponding
  residual R = 0
*/
void TacsDIRKIntegrator::nonlinearSolve(TacsScalar t, TacsScalar *q, 
				      TacsScalar *qdot, TacsScalar *qddot){

  // size of the linear system
  int size = numVars;
  
  // make space for residual and the jacobian
  TacsScalar *res = new TacsScalar[numVars];
  TacsScalar *D = new TacsScalar[numVars*numVars];

  // used by lapack
  int *dpiv = new int[numVars];

  // some tolerances 
  double atol = 1e-30; 
  double rtol = 1.0e-8;

  // initialize the norms
  TacsScalar init_norm = 0.0;
  TacsScalar norm = 0.0;

  // iterate until max iters or R <= tol
  int n;
  for (n = 0; n < max_newton_iters; n++) {
    
    // make sure the residual and jacobian are zeroed
    memset(res, 0, numVars*sizeof(TacsScalar));
    memset(D, 0, numVars*numVars*sizeof(TacsScalar));

    // get the residual
    computeResidual(res, t, q, qdot, qddot);
    
    // get the jacobian
    computeJacobian(D, 1.0, h*A[0], h*A[0]*h*A[0],
		     t, q, qdot, qddot);
    
    // Compute the l2 norm of the residual
    norm = 0.0;
    for ( int i = 0; i < numVars; i++ ){
      norm += res[i]*res[i];
    }
    norm = sqrt(norm);
    
    // Check if the norm of the residuals is a NaN - if so then
    // quit the integration
    if (norm != norm){ 
      fprintf(stderr, "Newton iteration %d, failed with NaN residual norm\n", n);
      break;
    }

    // Record the residual norm at the first Newton iteration
    if (n == 0){
      init_norm = norm;
    }

    // Check if the Newton convergence tolerance is satisfied
    if (norm < rtol*init_norm || norm < atol){
      break;
    }

    // Check whether the Newton iteration was successful
    if (n == max_newton_iters && norm >= rtol*init_norm){
      fprintf(stderr,
	      "Newton iteration failed to converge in %d iters \n", n);
      break;
    }

    // call lapack
    int info = 0;
    LAPACKgetrf(&size, &size, D, &size, dpiv, &info);
    if (info){
      fprintf(stderr,"LAPACK GETRF output error %d\n", info);
    }

    int one = 1;
    LAPACKgetrs("T", &size, &one, D, &size, dpiv, res, &size, &info);
    if (info){
      fprintf(stderr,"LAPACK GETRS output error %d\n", info);
    }

    // update the solution at the current iteration
    updateState(res, q, qdot, qddot);
        
  }

  // write a summary
  printf("Stage=%d, Stage Time=%f, # Newton iters=%d, |R|=%e \n", 
	 currentStage, t, n, norm);

  delete [] res;
  delete [] D;
  delete [] dpiv;

}

/*
  Update the state variable values after each Newton iteration

  Input: 
  The corresponding stage-states (q, qdot, qddot)

  Output:
  Updated stage-states (q, qdot, qddot)

*/

void TacsDIRKIntegrator::updateState(TacsScalar * res, 
				      TacsScalar *q, 
				      TacsScalar *qdot, 
				      TacsScalar *qddot) {

  // update qddot
  for (int i = 0; i < numVars; i++) {
    qddot[i] = qddot[i] - res[i];
  }

  // update qdot
  for (int i = 0; i < numVars; i++) {
    qdot[i] = qdot[i] - h*A[0]*res[i];
  }

  // update q
  for (int i = 0; i < numVars; i++) {
    q[i] = q[i] - h*A[0]*h*A[0]*res[i];
  }

}

/*
  Function that advances the states and time to next time step
  
  Input:
  The pointers to the global state variables q, qdot, qddot
  
  Output:
  Updated global state variables q, qdot, qddot at the time step
  
*/
void TacsDIRKIntegrator::timeMarch(int k, TacsScalar *time, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot){
  
  // advance the time
  time[k] = time[k-1] + h;
  
  // advance the state
  int ctr = 0;
  for (int i = 0; i < numVars; i++) {
    TacsScalar tmp = 0.0;
    for (int j = 0; j < numStages; j++) {
      tmp += B[j]*qdotS[ctr+j];
      ctr += numVars;
    }
    q[k*numVars+i] = q[(k-1)*numVars+i] + h*tmp;
  }

  // advance the velocity state
  ctr = 0;
  for (int i = 0; i < numVars; i++) {
    TacsScalar tmp = 0.0;
    for (int j = 0; j < numStages; j++) {
      tmp += B[j]*qddotS[ctr+j];
      ctr += numVars;
    }
    qdot[k*numVars+i] = qdot[(k-1)*numVars+i] + h*tmp;
  }
  
  // qddot?

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

/*
  The Residual implementation
  
  Input:
  The pointers at the current stage i.e. qS, qdotS and qddotS

  Output:
  res is filled with computed residual values

*/
void TacsDIRKIntegrator::computeResidual(TacsScalar *res, 
					  TacsScalar t, 
					  TacsScalar *q,
					  TacsScalar *qdot, 
					  TacsScalar *qddot){

  res[0] = qddot[0] + 0.02*qdot[0]*qdot[1] + 5.0*q[0];
  res[1] = qddot[1] - 0.05*qdot[1]*qdot[0] + q[1]*q[0];

}


/*
  The Jacobian implementation

  Input:
  The pointers at the current stage i.e. qS, qdotS and qddotS

  Output:
  J is filled with computed jacobian values

*/
void TacsDIRKIntegrator::computeJacobian(TacsScalar *J,
					  double alpha, double beta, double gamma, 
					  TacsScalar t, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot){

  // derivative wrt qddot
  J[0] = alpha*1.0;
  J[1] = alpha*0.0;
  J[2] = alpha*0.0;
  J[3] = alpha*1.0;

  // derivative wrt qdot
  J[0] += beta*0.02*qdot[1];
  J[1] += beta*0.02*qdot[0];
  J[2] += beta*-0.05*qdot[1];
  J[3] += beta*-0.05*qdot[0];

  // derivative wrt q
  J[0] += gamma*5.0;
  J[1] += gamma*0.0;
  J[2] += gamma*q[1];
  J[3] += gamma*q[0];

}








