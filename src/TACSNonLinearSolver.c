#include "TACSNonLinearSolver.h"
#include "tacslapack.h"

/*
  Default constructor for NonLinearSolver
 */
TacsNonLinearSolver::TacsNonLinearSolver(){
  
}

/*
  Constructor for Newton solver
 */
TacsNewtonSolver::TacsNewtonSolver(int _numVars, 
				   int _maxNewtonIters, 
				   double _atol, double _rtol) 
: TacsNonLinearSolver() {
  
  // copy over constructor arguments
  maxNewtonIters = _maxNewtonIters;
  numVars = _numVars;
  atol = _atol;
  rtol = _rtol;

}

/*
  Destructor for Newton solver
*/
TacsNewtonSolver::~TacsNewtonSolver(){

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
void TacsNewtonSolver::solve(double alpha, double beta, double gamma,
			     TacsScalar t, TacsScalar *q, TacsScalar *qdot, 
			     TacsScalar *qddot){

  int size = numVars;

  // make space for residual and the jacobian
  TacsScalar *res = new TacsScalar[size];
  TacsScalar *D = new TacsScalar[size*size];

  // used by lapack
  int *dpiv = new int[size];
  int info = 0;
  int one = 1;
 
  // initialize the norms
  TacsScalar init_norm = 0.0;
  TacsScalar norm = 0.0;

  // iterate until max iters or R <= tol
  int n;
  for (n = 0; n < maxNewtonIters; n++) {
    
    // make sure the residual and jacobian are zeroed
    memset(res, 0, size*sizeof(TacsScalar));
    memset(D, 0, size*size*sizeof(TacsScalar));

    // get the residual
    computeResidual(res, t, q, qdot, qddot);
    
    // get the jacobian
    computeJacobian(D, alpha, beta, gamma, t, q, qdot, qddot);

    // Compute the l2 norm of the residual
    norm = 0.0;
    for ( int i = 0; i < size; i++ ){
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
    if (n == maxNewtonIters && norm >= rtol*init_norm){
      fprintf(stderr,
	      "Newton iteration failed to converge in %d iters \n", n);
      break;
    }

    LAPACKgesv(&size, &one, D, &size, dpiv, res, &size, &info);
    if (info){
      fprintf(stderr,"LAPACK DGESV output error %d\n", info);
      break;
    }

    /*      
    // call lapack
    LAPACKgetrf(&size, &size, D, &size, dpiv, &info);
    if (info){
    fprintf(stderr,"LAPACK GETRF output error %d\n", info);
    break;
    }

 
    LAPACKgetrs("N", &size, &one, D, &size, dpiv, res, &size, &info);
    if (info){
    fprintf(stderr,"LAPACK GETRS output error %d\n", info);
    break;
    }
    */

    // update the state variables using the solution at the current
    // iteration
    updateState(res, alpha, beta, gamma, q, qdot, qddot);

  }

  // write a summary
  //  printf("Step=%d, Stage=%d, Stage Time=%f, # Newton iters=%d, |R|=%e \n", 
  //  currentTimeStep, currentStage, t, n, norm);

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

void TacsNewtonSolver::updateState(TacsScalar * res, 
				   double alpha, 
				   double beta, 
				   double gamma,
				   TacsScalar *q, 
				   TacsScalar *qdot, 
				   TacsScalar *qddot) {

  // update qddot
  for (int i = 0; i < numVars; i++) {
    qddot[i] = qddot[i] - alpha*res[i];
  }

  // update qdot
  for (int i = 0; i < numVars; i++) {
    qdot[i] = qdot[i] - beta*res[i];
  }

  // update q
  for (int i = 0; i < numVars; i++) {
    q[i] = q[i] - gamma*res[i];
  }

}


/*
  The Residual implementation
  
  Input:
  The pointers at the current stage i.e. qS, qdotS and qddotS

  Output:
  res is filled with computed residual values

*/
void TacsNewtonSolver::computeResidual(TacsScalar *res, 
					  TacsScalar t, 
					  TacsScalar *q,
					  TacsScalar *qdot, 
					  TacsScalar *qddot){

  res[0] = qddot[0] + 0.02*qdot[0]*qdot[1] + 5.0*q[0];
  res[1] = qddot[1] - 0.05*qdot[0]*qdot[1] + q[0]*q[1];
  
  // res[0] = qddot[0] + 0.02*qdot[0] + 5.0*q[0];

}


/*
  The Jacobian implementation

  Input:
  The pointers at the current stage i.e. qS, qdotS and qddotS

  Output:
  J is filled with computed jacobian values

*/
void TacsNewtonSolver::computeJacobian(TacsScalar *J,
					  double alpha, double beta, double gamma, 
					  TacsScalar t, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot){

  //  res[0] = qddot[0] + 0.02*qdot[0] + 5.0*q[0];

  /*
    J[0] = alpha*1.0;
    J[0] += beta*0.02;
    J[0] += gamma*5.0;
  */
  
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

