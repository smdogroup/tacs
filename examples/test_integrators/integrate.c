#include "TACSIntegrator.h"

/*
  Function that tests the integration logic in TACS
*/
int main( int argc, char *argv[] ){
  
  MPI_Init(&argc, &argv);

  int numStages = 1;
  int numVars = 2;
  
  double tInit = 0.0;
  double tFinal = 25.0;

  int numStepsPerSec = 10;
  int nSteps = int(double(numStepsPerSec)*(tFinal-tInit)) + 1;

  int max_newton_iters = 25;
  double atol = 1.0e-30;
  double rtol = 1.0e-8;
  
  // create space for states
  TacsScalar *q = new TacsScalar[numVars*nSteps];
  TacsScalar *qdot = new TacsScalar[numVars*nSteps];
  TacsScalar *qddot = new TacsScalar[numVars*nSteps];
  TacsScalar *time = new TacsScalar[nSteps];
  
  // initialize the states with zeroes
  memset(q, 0, numVars*nSteps*sizeof(TacsScalar));
  memset(qdot, 0, numVars*nSteps*sizeof(TacsScalar));
  memset(qddot, 0, numVars*nSteps*sizeof(TacsScalar));
  memset(time, 0, nSteps*sizeof(TacsScalar));

  // set initial condition
  time[0] = 0.0;

  // initial q values for all variables
  q[0] = 1.0;
  q[1] = 2.0;

  // initial qdot values for all variables
  qdot[0] = 0.0;
  qdot[1] = 0.0;

  printf("Integrating using DIRK \n");  

  // create an integrator object
  TacsIntegrator *dirk = new TacsDIRKIntegrator(numVars, 
						tInit, tFinal, 
						numStepsPerSec,
						max_newton_iters, 
						atol, rtol,
						numStages);

  dirk->incref();
  dirk->integrate(time, q, qdot, qddot);

  // write results
  FILE *fp = fopen("solution_dirk.dat", "w");
 
  int cnt = 0;
  for (int k = 0; k < nSteps; k++) {
    fprintf(fp, "%e ", time[k]);
    for (int j = 0; j < numVars; j++) {
      fprintf(fp, "%e %e %e ", q[cnt], qdot[cnt], qddot[cnt]);
      cnt++;
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  printf("DIRK complete \n");  
  
  // initialize the states with zeroes
  memset(q, 0, numVars*nSteps*sizeof(TacsScalar));
  memset(qdot, 0, numVars*nSteps*sizeof(TacsScalar));
  memset(qddot, 0, numVars*nSteps*sizeof(TacsScalar));
  memset(time, 0, nSteps*sizeof(TacsScalar));

  // set initial condition
  time[0] = 0.0;

  // initial q values for all variables
  q[0] = 1.0;
  q[1] = 2.0;

  // initial qdot values for all variables
  qdot[0] = 0.0;
  qdot[1] = 0.0;

  printf("Integrating using BDF\n");  

  int maxBDFOrder = 2;

  // create an integrator object
  TacsIntegrator *bdf = new TacsBDFIntegrator(numVars, 
					      tInit, tFinal, 
					      numStepsPerSec,
					      max_newton_iters, 
					      atol, rtol,
					      maxBDFOrder);
  
  bdf->incref();
  bdf->integrate(time, q, qdot, qddot);
  
  // write results
  fp = fopen("solution_bdf.dat", "w");
 
  cnt = 0;
  for (int k = 0; k < nSteps; k++) {
    fprintf(fp, "%e ", time[k]);
    for (int j = 0; j < numVars; j++) {
      fprintf(fp, "%e %e %e ", q[cnt], qdot[cnt], qddot[cnt]);
      cnt++;
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  printf("BDF complete \n");

  delete [] time;
  delete [] q;
  delete [] qdot;
  delete [] qddot;
 
  dirk->decref();
  bdf->decref();

  MPI_Finalize();

  return (0);
}




