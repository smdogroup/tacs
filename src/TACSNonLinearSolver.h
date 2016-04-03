#ifndef TACS_NONLINEAR_SOLVER_H
#define TACS_NONLINEAR_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "TACSObject.h"

/*
  Interface for nonlinear solvers in TACS
*/

class TacsNonLinearSolver : public TACSObject {

 public:

  TacsNonLinearSolver();
  virtual ~TacsNonLinearSolver(){}
  
  virtual void solve(double alpha, double beta, double gamma,
			      TacsScalar t, TacsScalar *q, TacsScalar *qdot, 
			      TacsScalar *qddot){}

  virtual void updateState(TacsScalar * res, 
			   double alpha, double beta,  double gamma,
			   TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot){}
};

/*
  A Newton-Raphson implementation of the non-linear solver
*/

class TacsNewtonSolver : public TacsNonLinearSolver {

 public:
  
  TacsNewtonSolver(int _numVars, int maxNewtonIters, double atol, double rtol);
  
  ~TacsNewtonSolver();
  
  void solve(double alpha, double beta, double gamma,
	     TacsScalar t, TacsScalar *q, TacsScalar *qdot, 
	     TacsScalar *qddot);
  
  void updateState(TacsScalar * res, 
		   double alpha, double beta,  double gamma,
		   TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);

  void computeResidual(TacsScalar *res, 
		       TacsScalar t, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);

  void computeJacobian(TacsScalar *J,
		       double alpha, double beta, double gamma, 
		       TacsScalar t, TacsScalar *q, TacsScalar *qdot, TacsScalar *qddot);
  
 private:

  int maxNewtonIters;
  int numVars;

  // some tolerances 
  double atol;
  double rtol;

};

#endif








