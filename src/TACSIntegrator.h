#ifndef TACS_INTEGRATOR_H
#define TACS_INTEGRATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "TACSObject.h"
#include "TACSNonLinearSolver.h"

/*
  Abstract class for integration schemes to extend
*/

class TacsIntegrator : public TACSObject {

 public:
 
  TacsIntegrator();
  virtual ~TacsIntegrator(){}
  virtual void integrate(TacsScalar *time, TacsScalar *q, 
			 TacsScalar *qdot, TacsScalar *qddot){}
  
 protected:
  
  int maxNewtonIters;
  
  double h;
  double tInit, tFinal;

  // integration coefficients
  double alpha, beta, gamma;

  int numStepsPerSec;
  int numSteps;
  
  int numVars;

  int currentTimeStep;

  TacsNonLinearSolver *nonlinearSolver;

};

/*
  A DIRK integration scheme for TACS
*/
class TacsDIRKIntegrator : public TacsIntegrator {

 public:

  // constructor for DIRK object
  TacsDIRKIntegrator(int numStages, int numVars, double tInit, 
		     double tFinal, int numStepsPerSec, 
		     int maxNewtonIters);
  
  // destructor for DIRK object
  ~TacsDIRKIntegrator();

  // integrate call
  void integrate(TacsScalar *time, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);
  
 private:

  //------------------------------------------------------------------//
  // private variables
  //------------------------------------------------------------------//

  // the number of stage in RK scheme
  int numStages;
  
  // the order of accuracy of the scheme
  int order;
  
  // variables for Butcher tableau
  double *A, *B, *C;
  
  // stage values (computed at each time stage for each time step)
  double  *tS, *qS, *qdotS, *qddotS;

  //------------------------------------------------------------------//
  // private functions
  //------------------------------------------------------------------//

  void setupButcherTableau();
  void checkButcherTableau();

  // returns the starting index of Butcher tableau
  int getIdx(int stageNum);

  void computeStageValues(TacsScalar tk, TacsScalar *qk, TacsScalar *qdotk, TacsScalar *qddotk);
  void timeMarch(int k, TacsScalar *time, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);
  void resetStageValues();

};

/*
  BDF integration scheme for TACS
*/

class TacsBDFIntegrator : public TacsIntegrator {

 public:
  
  // constructor for BDF object
  TacsBDFIntegrator(int numVars, double tInit, 
		    double tFinal, int numStepsPerSec, 
		    int maxNewtonIters, int maxBDFOrder);
  
  // destructor for BDF object
  ~TacsBDFIntegrator();
  
  // integrate call
  void integrate(TacsScalar *time, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);
  
 private:
  
  int maxBDFOrder;

};

#endif

