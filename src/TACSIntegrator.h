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
  
  TacsIntegrator(int numVars, double tInit, 
		 double tFinal, int numStepsPerSec, 
		 int maxNewtonIters, double _atol, double _rtol);
  
  ~TacsIntegrator();
  
  virtual void integrate(TacsScalar *time, TacsScalar *q, 
			 TacsScalar *qdot, TacsScalar *qddot){}
  
 protected:
 
  int numVars;
 
  double h;
  double tInit, tFinal;

  int numStepsPerSec;
  int numSteps; 
  int currentTimeStep;

  double alpha, beta, gamma;
  
  // nonlinear solver
  int maxNewtonIters;
  TacsNonLinearSolver *nonlinearSolver;

};

/*
  A DIRK integration scheme for TACS
*/
class TacsDIRKIntegrator : public TacsIntegrator {

 public:

  // constructor for DIRK object
  TacsDIRKIntegrator(int numVars, double tInit, 
		     double tFinal, int numStepsPerSec, 
		     int maxNewtonIters, double atol, double rtol,
		     int numStages);
  
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
  double *A=NULL, *B=NULL, *C=NULL;
  
  // stage values (computed at each time stage for each time step)
  double  *tS=NULL, *qS=NULL, *qdotS=NULL, *qddotS=NULL;

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
		    int maxNewtonIters, double atol, double rtol,
		    int maxBDFOrder);
  
  // destructor for BDF object
  ~TacsBDFIntegrator();
  
  // integrate call
  void integrate(TacsScalar *time, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);
  
 private:
  
  TacsScalar *res = NULL;
  TacsScalar *D = NULL;
    
  // The BDF coefficients for time integration
  double bdf_coeff[4], bddf_coeff[9];
  
  int maxBDFOrder;

  // Retireve the BDF coefficients
  int getBDFCoeff(double bdf[], int order );
  void get2ndBDFCoeff(const int k, double bdf[], int *nbdf,
		      double bddf[], int *nbddf,
		      const int max_order);
  
  // approximate derivatives using BDF stencil
  void approxStates(TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);
  void setCoeffs(double alpha, double beta, double gamma);

};

#endif

