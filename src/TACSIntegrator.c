#include "TACSIntegrator.h"
#include "tacslapack.h"

/*
  Constructor for TacsDIRKIntegrator
  
  input:
  tInit: initial time
  h: time step size
  num_stages: number of stages in implicit Runge Kutta

*/
TacsDIRKIntegrator::TacsDIRKIntegrator(int _numStages, int _numVars,  double _tInit, double _tFinal, int _numStepsPerSec){

  // copy over the input parameters
  numStages = _numStages;
  numVars = _numVars;
  tInit = _tInit;
  tFinal = _tFinal;
  numStepsPerSec = _numStepsPerSec;

  // sanity check the input params

  // compute the step size
  h = 1.0/numStepsPerSec;

  // compute the total number of time steps
  numSteps = int(numStepsPerSec*(tFinal-tInit));

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
  
  // allocate space for the stage residual and jacobian. Since we
  // solve for each stage successively, so we don't have to make space
  // for all the stages
  RS = new TacsScalar[numVars];
  JS = new TacsScalar[numVars*numVars];

  // initialize stage residual and jacobian
  memset(RS, 0, numVars*sizeof(TacsScalar));
  memset(JS, 0, numVars*numVars*sizeof(TacsScalar));

  // allocate space for Butcher tableau
  A = new double[numStages*(numStages+1)/2];
  B = new double[numStages];
  C = new double[numStages];

  // set the Butcher tableau entries to zero
  memset(A, 0, numStages*(numStages+1)/2*sizeof(double));
  memset(B, 0, numStages*sizeof(double));
  memset(C, 0, numStages*sizeof(double));

  // now add entries into the Butcher tableau
  setup_butcher_tableau();
  
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

  delete [] RS;
  delete [] JS;

}

/*
  Function that adds the entries into Butcher tableau
*/
void TacsDIRKIntegrator::setup_butcher_tableau(){
  
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
  check_butcher_tableau();

}

/*
  Function that checks the consistency of Butcher tableau values
*/
void TacsDIRKIntegrator::check_butcher_tableau(){

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
      fprintf(stderr, "WARNING: Sum A[%d,:] != C[%d] i.e. %f != %f \n", i, i, C[i], tmp);
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
void TacsDIRKIntegrator::integrate(TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot){
  
  for (int k = 1; k <= numSteps; k++) {
    
    // find the stage derivatives at the current step
    int idx = (k-1)*numSteps;

    // pass in the global states at previous time step
    compute_stage_values(double(k)*h, &q[idx], &qdot[idx], &qddot[idx]);
    
    // advance the state to the current step
    time_march(k, q, qdot, qddot);

    // set the stage values to zero
    reset_stage_values();
    
  }
  
}

/*
  Function that computes the values at each time step. The values that
  are computed are qS, tS and qdotS. These are used in subsequently in
  the call to time_march function.
*/
void TacsDIRKIntegrator::compute_stage_values(TacsScalar tk, TacsScalar *qk, TacsScalar *qdotk, TacsScalar *qddotk){

  currentStage = 0;

  // stage offset for q values
  int soffset = 0;
  
  for (int i = 0; i < numStages; i++) {

    // printf("time=%f, stage = %d, soffset=%d\n", tk, i, soffset);

    // compute the stage time
    tS[i] = time + C[i]*h;

    // initial guess for qddotS
    for (int n = 0; n < numVars; n++) {
      qddotS[soffset+n] = 1.0;
    }

    // prepare qdotS
    for (int n = 0; n < numVars; n++) {  
      TacsScalar tmp = 0.0;
      int idx1 = getIdx(i+1);
      for (int j = 0; j <= i; j++) {
	tmp += A[idx1]*qddotS[n];
	idx1++;
      }
      qdotS[soffset+n] = qdotk[n] + h* tmp;
    }
    
    // prepare qS
    for (int n = 0; n < numVars; n++) {  
      TacsScalar tmp = 0.0;
      int idx2 = getIdx(i+1);
      for (int j = 0; j <= i; j++) {
	tmp += A[idx2]*qdotS[n];
	idx2++;
      }
      qS[soffset+n] = qk[n] + h* tmp;
    }
    
    // solve the nonlinear system of stage equations for qddotS
    newton_solve(tS[i], &qdotS[soffset], &qdotS[soffset], &qddotS[soffset]);

    // recompute the offsets
    soffset +=  numVars;

    // increament the global stage counter
    currentStage ++;
    
  }
  
}

/*
  Start index for the supplied stage  
*/
int TacsDIRKIntegrator::getIdx(int stageNum){
  
  if (stageNum == 1) {
    return 0;
  } else if (stageNum == 2) {
    return 1;
  } else if (stageNum == 3) {
    return 3;
  }  else if (stageNum == 4) {
    return 6;
  }

}

/*
  Solve the stage equations using Newton's method
*/
void TacsDIRKIntegrator::newton_solve(TacsScalar tS, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot){

  int size = numVars;
  
  TacsScalar *res = new TacsScalar[numVars];
  TacsScalar *D = new TacsScalar[numVars*numVars];

  int *dpiv = new int[numVars];

  double atol = 1e-30; 
  double rtol = 1.0e-8;

  TacsScalar init_norm = 0.0;
  TacsScalar norm = 0.0;

  for (int n = 0; n < max_newton_iters; n++) {
    
    // get the residual
    compute_residual(res, tS, qS, qdotS, qddotS);
    
    // get the jacobian
    compute_jacobian(D, tS, qS, qdotS, qddotS);

    // Compute the l2 norm of the residual
    double norm = 0.0;
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
    if (n >= max_newton_iters && norm >= rtol*init_norm){
      fprintf(stderr,
	      "Newton iteration for time step %d failed to converge\n",
	      n);
      break;
    }

    // call lapack
    int info = 0;
    LAPACKgetrf(&size, &size, D, &size, dpiv, &info);
    if (info){
      fprintf(stderr,"LAPACK GETRF output error %d\n", info);
    }

    int one = 1;
    LAPACKgetrs("N", &size, &one, D, &size, dpiv, res, &size, &info);
    if (info){
      fprintf(stderr,"LAPACK GETRS output error %d\n", info);
    }

    // update the solution at the current iteration
    state_update(&res[0]);

  }

  delete [] res;
  delete [] D;
  delete [] dpiv;

}

/*
  Update the state variable values after each Newton iteration
*/

void TacsDIRKIntegrator::state_update(TacsScalar * res) {

  // update qddot
  for (int i = 0; i < numVars; i++) {
    qddotS[i] = qddotS[i] - res[i];
  }

  // update qdot
  for (int i = 0; i < numVars; i++) {
    qdotS[i] = qdotS[i] - h*A[0]*res[i];
  }

  // update q
  for (int i = 0; i < numVars; i++) {
    qS[i] = qS[i] - h*A[0]*h*A[0]*res[i];
  }

}

/*
  Function that advances the states and time to next time step
*/
void TacsDIRKIntegrator::time_march(int k, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot){
  
  // advance the time
  time = time + h;
  
  // compute the summation
  TacsScalar tmp = 0.0;
  for (int j = 0; j < numStages; j++) {
    tmp += B[j]*qdotS[j];
  }
  
  // advance the state
  q[k] = q[k-1] + h*tmp;

  // compute the summation
  tmp = 0.0;
  for (int j = 0; j < numStages; j++) {
    tmp += B[j]*qddotS[j];
  }

  // advance the velocity states
  qdot[k] = qdot[k-1] + h*tmp;

}

/*
  Function that sets the arrays used during previous time marching
  step to zero, to make way for the next time step
*/
void TacsDIRKIntegrator::reset_stage_values(){
  
  memset(tS, 0, numStages*sizeof(TacsScalar));
  memset(qS, 0, numStages*sizeof(TacsScalar));
  memset(qdotS, 0, numStages*sizeof(TacsScalar));
  memset(qddotS, 0, numStages*sizeof(TacsScalar));

  memset(RS, 0, numStages*sizeof(TacsScalar));
  memset(JS, 0, numStages*numStages*sizeof(TacsScalar));

}

/*
  Residual implementation
*/
void TacsDIRKIntegrator::compute_residual(TacsScalar *R, 
					  TacsScalar t, 
					  TacsScalar *q,
					  TacsScalar *qdot, 
					  TacsScalar *qddot){

}


/*
  Jacobian implementation
*/
void TacsDIRKIntegrator::compute_jacobian(TacsScalar *J, 
					  TacsScalar t, 
					  TacsScalar *q, 
					  TacsScalar *qdot, 
					  TacsScalar *qddot){
}

					  
					  


					  










