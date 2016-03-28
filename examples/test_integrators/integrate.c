#include "TACSIntegrator.h"

/*
  Function that tests the integration logic in TACS
*/
int main( int argc, char *argv[] ){
  
  MPI_Init(&argc, &argv);

  int numStages = 3;
  int numVars = 2;
  
  double tInit = 0.0;
  double tFinal = 1.0;

  int numStepsPerSec = 10;
  int nSteps = numStepsPerSec*(tFinal-tInit) + 1;

  int max_newton_iters = 25;

  // create space for states
  TacsScalar *q = new TacsScalar[numVars];
  TacsScalar *qdot = new TacsScalar[numVars];
  TacsScalar *qddot = new TacsScalar[numVars];
  TacsScalar *time = new TacsScalar[nSteps];
  
  // initialize the states with zeroes
  memset(q, 0, numVars*sizeof(TacsScalar));
  memset(qdot, 0, numVars*sizeof(TacsScalar));
  memset(qddot, 0, numVars*sizeof(TacsScalar));

  // set initial condition
  time[0] =0.0;

  // initial q values for all variables
  q[0] = 1.0;
  q[1] = 2.0;

  // initial qdot values for all variables
  qdot[0] = 0.0;
  qdot[1] = 0.0;
  
  // create an integrator object
  TacsDIRKIntegrator *dirk = new TacsDIRKIntegrator(numStages, numVars, 
						    tInit, tFinal, 
						    numStepsPerSec,
						    max_newton_iters);
  dirk->incref();

  dirk->integrate(time, q, qdot, qddot);

  // write results
  FILE *fp = fopen("solution.dat", "w");
  
  int cnt = 0;
  for (int k = 0; k < nSteps; k++) {
    fprintf(fp, "%f ", time[k]);
    for (int j = 0; j < numVars; j++) {
      fprintf(fp, "%f ", q[cnt]);
      cnt++;
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  dirk->decref();

  delete [] q;
  delete [] qdot;
  delete [] qddot;

  MPI_Finalize();

  return (0);
}




