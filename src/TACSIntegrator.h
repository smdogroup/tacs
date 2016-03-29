#ifndef TACS_INTEGRATOR_H
#define TACS_INTEGRATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "TACSObject.h"

/*
  A DIRK integration scheme for TACS
*/
class TacsDIRKIntegrator : public TACSObject {

 public:

  // constructor for DIRK object
  TacsDIRKIntegrator(int numStages, int numVars, double tInit, 
		     double tFinal, int numStepsPerSec, 
		     int max_newton_iters);
  
  // destructor for DIRK object
  ~TacsDIRKIntegrator();

  // integrate call
  void integrate(TacsScalar *time, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);
  
 private:

  // private variables
  double h;
  double tInit, tFinal;

  int numStepsPerSec;
  int numSteps;

  int numVars;
  int numStages, order;

  int max_newton_iters;
  
  // variables for Butcher tableau
  double *A, *B, *C;

  // stage values (computed at each time stage for each time step)
  double  *tS, *qS, *qdotS, *qddotS;
  
  // global index of stage
  int stageOffCtr;
  
  // variables to keep track of current stage and time index
  int currentStage;
  int currentTimeStep;

  // private functions

  void setupButcherTableau();
  void checkButcherTableau();

  // returns the starting index of Butcher tableau
  int getIdx(int stageNum);

  void computeStageValues(TacsScalar tk, TacsScalar *qk, TacsScalar *qdotk, TacsScalar *qddotk);
  void timeMarch(int k, TacsScalar *time, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);

  void resetStageValues();

  void computeResidual(TacsScalar *Rki, TacsScalar tS, TacsScalar *qS,
			TacsScalar *qdotS, TacsScalar *qddotS);  
  void computeJacobian(TacsScalar *J,
			double alpha, double beta, double gamma,
			TacsScalar t, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);

  void nonlinearSolve(TacsScalar t, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);
  void updateState(TacsScalar * sol, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);
 
};

#endif
