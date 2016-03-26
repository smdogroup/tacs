#ifndef TACS_INTEGRATOR_H
#define TACS_INTEGRATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "TACSObject.h"

/*
  A DIRK integration scheme for TACS
*/
class TacsDIRKIntegrator {

 public:

  // constructor for DIRK object
  TacsDIRKIntegrator(int numStages, int numVars, double tInit, double tFinal, int numStepsPerSec);
  
  // dectructor for DIRK object
  ~TacsDIRKIntegrator();

  // integrate call
  void integrate(TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);
  
 private:

  // private variables

  double h;
  double tInit, tFinal;
  double time;

  int numStepsPerSec;
  int numSteps;

  int numVars;
  
  int currentStage;

  int numStages, order;
  
  // variables for Butcher tableau
  double *A, *B, *C;

  // stage values (computed at each time step)
  double  *tS, *qS, *qdotS, *qddotS;
  
  // stage residual and jacobianx
  double *RS, *JS;
  
  int max_newton_iters = 25;
  
  // private functions
  void setup_butcher_tableau();
  void check_butcher_tableau();

  // returns the starting index of Butcher tableau
  int getIdx(int stageNum);

  void compute_stage_values(TacsScalar tk, TacsScalar *qk, TacsScalar *qdotk, TacsScalar *qddotk);
  void time_march(int k, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);
  void reset_stage_values();

  // external function and residual assembly as needed by the scheme
  void compute_residual(TacsScalar *Rki, TacsScalar tS, TacsScalar *qS,
			TacsScalar *qdotS, TacsScalar *qddotS);
  void compute_jacobian(TacsScalar *Jki, TacsScalar tS, TacsScalar *qS, 
			TacsScalar *qdotS, TacsScalar *qddotS);

  // external interface for solving the linear system
  void newton_solve(TacsScalar t, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);
  void state_update(TacsScalar * res);

  // function writes the time history to a file
  void writeSolution();

};

#endif







