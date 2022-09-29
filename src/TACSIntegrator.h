/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_INTEGRATOR_H
#define TACS_INTEGRATOR_H

#include "KSM.h"
#include "TACSAssembler.h"
#include "TACSObject.h"
#include "TACSToFH5.h"

/*
  Abstract base class for integration schemes.

  This base class contains common functions and variables pertaining
  to the integration and adjoint schemes used in TACS.

  Certain functions are pure virtual and therefore the integrator can
  be instantiated only after those functions are provided a full
  implementation in an extending child class.
*/

class TACSIntegrator : public TACSObject {
 public:
  TACSIntegrator(TACSAssembler *assembler, double tinit, double tfinal,
                 double num_steps);
  ~TACSIntegrator();

  // Define members of the class
  // ---------------------------
  void setRelTol(double _rtol);
  void setAbsTol(double _atol);
  void setMaxNewtonIters(int _max_newton_iters);
  void setPrintLevel(int _print_level, const char *logfilename = NULL);
  void setJacAssemblyFreq(int _jac_comp_freq);
  void setUseLapack(int _use_lapack);
  void setUseSchurMat(int _use_schur_mat, TACSAssembler::OrderingType _type);
  void setInitNewtonDeltaFraction(double frac);
  void setKrylovSubspaceMethod(TACSKsm *_ksm);

  // Set (or reset) the time interval
  // --------------------------------
  void setTimeInterval(double tinit, double tfinal);

  // Set the functions to integrate
  //--------------------------------
  void setFunctions(int num_funcs, TACSFunction **funcs, int start_plane = -1,
                    int end_plane = -1);

  // Solve for time-step t with
  virtual int iterate(int step_num, TACSBVec *forces) = 0;

  // Integrate the equations of motion forward in time
  virtual int integrate();

  // Evaluate the functions of interest
  virtual void evalFunctions(TacsScalar *fvals) = 0;

  // Set-up right-hand-sides for the adjoint equations
  virtual void initAdjoint(int step_num) {}

  // Iterate to find a solution of the adjoint equations
  virtual void iterateAdjoint(int step_num, TACSBVec **adj_rhs) {}

  // Add the contributions to the total derivative from the time-step
  virtual void postAdjoint(int step_num) {}

  // Integrate the adjoint and add the total derivative from all
  // time-steps
  virtual void integrateAdjoint();

  // Get the adjoint vector for the given function
  virtual void getAdjoint(int step_num, int func_num, TACSBVec **adjoint) = 0;

  // Copy out the function
  virtual void getGradient(int func_num, TACSBVec **_dfdx) {
    *_dfdx = NULL;
    if (func_num >= 0 && func_num < num_funcs) {
      *_dfdx = dfdx[func_num];
    }
  }

  // Get the vector of the shape derivatives
  virtual void getXptGradient(int func_num, TACSBVec **_dfdXpt) {
    *_dfdXpt = NULL;
    if (func_num >= 0 && func_num < num_funcs) {
      *_dfdXpt = dfdXpt[func_num];
    }
  }

  // Retrieve the internal states
  double getStates(int step_num, TACSBVec **q, TACSBVec **qdot,
                   TACSBVec **qddot);

  // Initialize linear solver
  //--------------------------
  void initializeLinearSolver();

  // Check the gradient using finite-difference
  void checkGradients(double dh);

  // Functions to export the solution in raw and tecplot binary forms
  //-----------------------------------------------------------------
  void setOutputPrefix(const char *prefix);
  void setOutputFrequency(int _write_step);
  void setFH5(TACSToFH5 *_f5);
  void writeRawSolution(const char *filename, int format = 2);
  void writeSolutionToF5();
  void writeStepToF5(int step_num);
  void printWallTime(double t0, int level = 1);
  void printOptionSummary();
  void printAdjointOptionSummary();

  // Returns the number of time steps configured during instantiation
  //-----------------------------------------------------------------
  int getNumTimeSteps();
  int lapackNaturalFrequencies(int use_gyroscopic, TACSBVec *q, TACSBVec *qdot,
                               TACSBVec *qddot, TacsScalar *eigvals,
                               TacsScalar *modes = NULL);
  void getRawMatrix(TACSMat *mat, TacsScalar *mat_vals);

 protected:
  // Functions for solutions to linear and nonlinear problems
  int newtonSolve(double alpha, double beta, double gamma, double t,
                  TACSBVec *q, TACSBVec *qdot, TACSBVec *qddot,
                  TACSBVec *forces = NULL);
  void lapackLinearSolve(TACSBVec *res, TACSMat *mat, TACSBVec *update);

  // Variables that keep track of time
  double time_fwd_assembly;
  double time_fwd_factor;
  double time_fwd_apply_factor;
  double time_forward;
  double time_newton;

  double time_rev_assembly;
  double time_rev_factor;
  double time_rev_apply_factor;
  double time_rev_jac_pdt;
  double time_reverse;

  // Log the time step information
  void logTimeStep(int time_step);

  // TACSAssembler information
  TACSAssembler *assembler;  // Instance of TACSAssembler

  // The step information
  int num_time_steps;  // Total number of time steps
  double *time;        // Stores the time values
  TACSBVec **q;        // state variables across all time steps
  TACSBVec **qdot;     // first time derivative of ''
  TACSBVec **qddot;    // second time derivative of ''

  // Objects that store information about the functions of interest
  int start_plane, end_plane;  // Time-window for the functions of interest
  int num_funcs;               // The number of objective functions
  TACSFunction **funcs;        // List of functions
  TacsScalar *fvals;           // Function values
  TACSBVec **dfdx;             // Derivative values
  TACSBVec **dfdXpt;           // Derivatives w.r.t. node locations

  // Linear algebra objects and parameters associated with the Newton solver
  TACSBVec *res, *update;         // Residual and Newton update
  TACSMat *mat;                   // Jacobian matrix
  TACSPc *pc;                     // Preconditioner
  TACSKsm *ksm;                   // KSM solver
  int linear_solver_initialized;  // flag to indicate whether the linear solver
                                  // is initialized
  int mpiRank;                    // rank of the processor
  int mpiSize;                    // number of processors
  FILE *logfp;                    // Pointer to the output filename

  // Newton solver parameters
  int max_newton_iters;      // The max number of nonlinear iterations
  double atol;               // Absolute tolerance
  double rtol;               // Relative tolerance
  double init_newton_delta;  // Initial value of delta for globalization
  int jac_comp_freq;         // Frequency of Jacobian factorization
  int use_schur_mat;         // use the Schur matrix type for parallel execution
  TACSAssembler::OrderingType order_type;
  int use_lapack;  // Flag to switch to LAPACK for linear solve

  int lev;
  double fill;
  int reorder_schur;
  int gmres_iters;
  int num_restarts;
  int is_flexible;

 private:
  char prefix[256];  // Output prefix

  // Information for visualization/logging purposes
  int print_level;    // 0 = off;
                      // 1 = summary per time step;
                      // 2 = summary per Newton solve iteration
  TACSToFH5 *f5;      // F5 output visualization
  int f5_write_freq;  // Frequency for output during time marching

  int niter;                 // Newton iteration number
  TacsScalar res_norm;       // residual norm
  TacsScalar init_res_norm;  // Initial norm of the residual
  TacsScalar update_norm;    // Norm of the update

  TacsScalar init_energy;  // The energy during time = 0
};

/*
  Backwards Difference Formulas integration scheme
*/
class TACSBDFIntegrator : public TACSIntegrator {
 public:
  TACSBDFIntegrator(TACSAssembler *_tacs, double _tinit, double _tfinal,
                    double _num_steps, int max_bdf_order);
  ~TACSBDFIntegrator();

  // Iterate through the forward solution
  int iterate(int k, TACSBVec *forces);

  // Set-up right-hand-sides for the adjoint equations
  void initAdjoint(int step_num);

  // Iterate to find a solution of the adjoint equations
  void iterateAdjoint(int step_num, TACSBVec **adj_rhs);

  // Add the contributions to the total derivative from the time-step
  void postAdjoint(int step_num);

  // Get the adjoint value for the given function
  void getAdjoint(int step_num, int func_num, TACSBVec **adjoint);

  // Evaluate the functions of interest
  void evalFunctions(TacsScalar *fvals);

 private:
  void get2ndBDFCoeff(const int k, double bdf[], int *nbdf, double bddf[],
                      int *nbddf, const int max_order);
  int getBDFCoeff(const int k, double bdf[], int order);

  int max_bdf_order;  // Maximum order of the BDF integration scheme

  // Adjoint information
  int num_adjoint_rhs;  // the number of right hand sides allocated
  TACSBVec **rhs;       // storage vector for the right hand sides
  TACSBVec **psi;       // adjoint variable accumulating qdot dependance
};

/*
  DIRK integration scheme for TACS
*/
class TACSDIRKIntegrator : public TACSIntegrator {
 public:
  TACSDIRKIntegrator(TACSAssembler *_tacs, double _tinit, double _tfinal,
                     double _num_steps, int _num_stages);
  ~TACSDIRKIntegrator();

  // Iterate through the forward solution
  int iterate(int k, TACSBVec *forces);

  // Iterate through the forward solution - per stage
  int iterateStage(int k, int s, TACSBVec *forces);

  // Retrieve the internal states - per stage
  double getStageStates(int step, int stage, TACSBVec **qS, TACSBVec **qdotS,
                        TACSBVec **qddotS);

  // Set-up right-hand-sides for the adjoint equations
  void initAdjoint(int step_num);

  // Iterate to find a solution of the adjoint equations
  void iterateAdjoint(int step_num, TACSBVec **adj_rhs);

  // Add the contributions to the total derivative from the time-step
  void postAdjoint(int step_num);

  // Get the adjoint value for the given function
  void getAdjoint(int step_num, int func_num, TACSBVec **adjoint);

  // Evaluate the functions of interest
  void evalFunctions(TacsScalar *fvals);

 private:
  // Set the default coefficients
  void setupDefaultCoeffs();

  // Set the second-order coefficients based on the first-order values
  void setupSecondCoeffs();

  // Get the stage coefficient from the Tableau
  double getACoeff(const int i, const int j);

  // Check the Butcher tableau for consistency
  void checkButcherTableau();

  // Get the row index for the stage
  int getRowIndex(int stageNum);

  // Get the linearization coefficients for the given stage/step
  void getLinearizationCoeffs(const int stage, const double h, double *alpha,
                              double *beta, double *gamma);

  // The number of stages for this method
  int num_stages;

  // States at each stage
  TACSBVec **qS, **qdotS, **qddotS;

  // The Butcher coefficients for the integration scheme
  double *a, *b, *c;

  // The second-order coefficients for the integration scheme
  double *A, *B;

  // Right-hand-side vector
  TACSBVec *rhs;

  // Store the adjoint-stage vectors
  TACSBVec **lambda;

  // Save the stage right-hand-side integrator
  TACSBVec **omega, **domega;

  // Save the right-hand-side
  TACSBVec **psi, **phi;
};

/*
  ESDIRK integration scheme for TACS. *No adjoint implementation yet*
*/
class TACSESDIRKIntegrator : public TACSIntegrator {
 public:
  // Constructor
  TACSESDIRKIntegrator(TACSAssembler *_tacs, double _tinit, double _tfinal,
                       double _num_steps, int _num_stages);

  // Destructor
  ~TACSESDIRKIntegrator();

  // Iterate through the forward solution
  int iterate(int k, TACSBVec *forces);

  // Iterate through the forward solution - per stage
  int iterateStage(int k, int s, TACSBVec *forces);

  // Retrieve the internal states - per stage
  double getStageStates(int step_num, int stage_num, TACSBVec **qS,
                        TACSBVec **qdotS, TACSBVec **qddotS);

  // Evaluate the functions of interest
  void evalFunctions(TacsScalar *fvals);

  // Get the adjoint value for the given function - adjoint not implemented yet
  void getAdjoint(int step_num, int func_num, TACSBVec **adjoint);

 private:
  // set the first-order descirption integration coefficients
  void setupDefaultCoeffs();

  // set the second-order description integration coefficients
  void setupSecondCoeffs();

  // return the stage coefficients from the Butcher Tableau
  double getACoeff(const int i, const int j);

  // check the Butcher Tableau for consistency
  void checkButcherTableau();

  // get the row index of the a/A coefficient matrix for the stage
  int getRowIndex(int stage_num);

  // get the linearization coefficients for the given step/stage
  void getLinearizationCoeffs(const int stage, const double h, double *alpha,
                              double *beta, double *gamma);

  // the number of stages for this method
  int num_stages;

  // the states at every stage
  TACSBVec **qS, **qdotS, **qddotS;

  // the (first-order) Butcher Tableau coefficients for the integration scheme
  double *a, *b, *c;

  // the second order coefficients for the integration scheme
  double *A, *B;
};

/*
  Adams-Bashforth-Moulton integration scheme for TACS
*/
/*
class TACSABMIntegrator : public TACSIntegrator {
 public:
  TACSABMIntegrator( TACSAssembler *tacs,
                     double tinit,
                     double tfinal,
                     double num_steps,
                     int max_abm_order );
  ~TACSABMIntegrator();

  // Iterate through the forward solution
  int iterate( int k, TACSBVec *forces );

  // Set-up right-hand-sides for the adjoint equations
  void initAdjoint( int step_num );

  // Iterate to find a solution of the adjoint equations
  void iterateAdjoint( int step_num, TACSBVec **adj_rhs );

  // Add the contributions to the total derivative from the time-step
  void postAdjoint( int step_num );

  // Get the adjoint value for the given function
  void getAdjoint( int step_num, int func_num,
                   TACSBVec **adjoint );

  // Evaluate the functions of interest
  void evalFunctions( TacsScalar *fvals );

 private:
  void setupCoeffs();
  void checkABMCoeffs();
  int getCoeffIndex( int time_step );
  int getRowIdx( int row_num );
  int getOrder( int k );

  int max_abm_order; // Order of integration
  double *A;         // Array of ABM coefficients
};
*/

/*
  Newmark Beta Gamma Method (NBG) integration scheme for TACS
*/
/*
class TACSNBGIntegrator : public TACSIntegrator {
 public:
  TACSNBGIntegrator( TACSAssembler *tacs,
                     double tinit,
                     double tfinal,
                     double num_steps,
                     int order );
  ~TACSNBGIntegrator();

  // Iterate through the forward solution
  int iterate( int k, TACSBVec *forces );

  // Set-up right-hand-sides for the adjoint equations
  void initAdjoint( int step_num );

  // Iterate to find a solution of the adjoint equations
  void iterateAdjoint( int step_num, TACSBVec **adj_rhs );

  // Add the contributions to the total derivative from the time-step
  void postAdjoint( int step_num );

  // Get the adjoint value for the given function
  void getAdjoint( int step_num, int func_num,
                   TACSBVec **adjoint );

  // Evaluate the functions of interest
  void evalFunctions( TacsScalar *fvals );

 private:
  void setupCoeffs();
  int order;          // order of integration
  double BETA, GAMMA; // Newmark coefficients
};
*/
#endif  // TACS_INTEGRATOR_H
