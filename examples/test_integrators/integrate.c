#include "TACSIntegrator.h"

/*
  Function that tests the integration logic in TACS
*/
int main( int argc, char *argv[] ){
  
  MPI_Init(&argc, &argv);

  int numStages = 3;
  int numVars = 5;
  
  double tInit = 0.0;
  double tFinal = 0.2;
  int numStepsPerSec = 10;
  
  // create space for states
  TacsScalar *q = new TacsScalar[numVars];
  TacsScalar *qdot = new TacsScalar[numVars];
  TacsScalar *qddot = new TacsScalar[numVars];

  // initialize the states with zeroes
  memset(q, 0, numVars*sizeof(TacsScalar));
  memset(qdot, 0, numVars*sizeof(TacsScalar));
  memset(qddot, 0, numVars*sizeof(TacsScalar));

  // set initial condition
  q[0] = 1.0;
  qdot[0] = 0.0;
  
  // create an integrator object
  TacsDIRKIntegrator *dirk = new TacsDIRKIntegrator(numStages, numVars, tInit, tFinal, numStepsPerSec);
  
  dirk->integrate(q, qdot, qddot);

  // print results
  for (int k = 0; k < numVars; k++) {
    printf("%f \n", q[k]);
  }
  
  delete dirk;

  delete [] q;
  delete [] qdot;
  delete [] qddot;

  MPI_Finalize();

  return (0);
}




