#ifndef TACS_INTEGRATOR_H
#define TACS_INTEGRATOR_H

#include <stdarg.h>
#include "TACSObject.h"
#include "TACSAssembler.h"
#include "KSM.h"
#include "TACSToFH5.h"

// Type of integrator to use. The following are the supported methods.
//--------------------------------------------------------------------
enum IntegratorType { BDF1, BDF2, BDF3,                   // Backward-difference methods
                      ABM1, ABM2, ABM3, ABM4, ABM5, ABM6, // Adams-Bashforth-method
                      DIRK2, DIRK3, DIRK4,                // Diagonally-Implicit-Runge-Kutta methods
                      NBG};                               // Newmark Beta Gamma Method

/*
  Base class for integration schemes. This base class contains common
  methods and variables pertaining to the integration schemes used in
  TACS.
  
  Copyright (c) 2010-2016 Graeme Kennedy. All rights reserved. 
*/

class TACSIntegrator : public TACSObject {
 public:
  // Factory method for intantiating integrators of supplied type
  //-------------------------------------------------------------
  static TACSIntegrator* getInstance( TACSAssembler *_tacs, 
                                      double _tinit, double _tfinal, 
                                      int _num_steps_per_sec, 
                                      enum IntegratorType type );

  // Useful setters for class variables. These setters must be called
  // before calling the 'integrate()' function. These won't have any
  // effect after the integrate call has been made as default values
  // would have been used during time integration.
  // -----------------------------------------------------------------
  void setRelTol( double _rtol );
  void setAbsTol( double _atol );
  void setMaxNewtonIters( int _max_newton_iters );
  void setPrintLevel( int _print_level );
  void setJacAssemblyFreq( int _jac_comp_freq );
  void setUseLapack( int _use_lapack );

  // Set the objective/constraint functions of interest and increment.
  // -----------------------------------------------------------------
  void setFunction( TACSFunction **_func, int _num_funcs );

  // Call this function after integrating to write the solution to
  // file in ASCII text format (might be slower for bigger problems)
  // ------------------------------------------------------------------
  void writeSolution( const char *filename );

  // Write a single step to F5 file created based on the step number
  //----------------------------------------------------------------
  void writeStepToF5( int step=0 );
  
  // Call this function after integrating to write the solution to f5
  // file. Since it is written in binary form, it is faster.
  // ------------------------------------------------------------------
  void writeSolutionToF5();

  // Configure the F5 file output
  //-----------------------------
  void configureOutput(TACSToFH5 *_viewer, int _write_freq, char *_f5_file_fmt);
  
  // Get the finite-difference/complex-step gradient for testing purposes
  //---------------------------------------------------------------------
  void getFDFuncGrad( int num_dv, TacsScalar *x,
		      TacsScalar *fvals, TacsScalar *dfdx, double dh );
  
  // Function for returning the adjoint derivative for the functions
  //----------------------------------------------------------------
  void getFuncGrad( int num_dv, TacsScalar *x,
		    TacsScalar *fvals, TacsScalar *dfdx );

  // Pure virtual function that the derived classes must override/implement
  //-----------------------------------------------------------------------
  virtual void integrate() = 0;

 protected:
  // Parent class constructor for integration schemes
  //-------------------------------------------------
  TACSIntegrator( TACSAssembler *_tacs, 
                  double _tinit, double _tfinal, 
                  int _num_steps_per_sec ); 

  // Destructor
  //-----------
  ~TACSIntegrator();

  // Pure virtual function for marching backwards in stage and time
  //---------------------------------------------------------------
  virtual void marchBackwards() = 0;

  // Evaluate time average of the function value using discretization
  // from the integration scheme
  //---------------------------------------------------------------------
  virtual void evalTimeAvgFunctions( TACSFunction **funcs, int numFuncs, TacsScalar *funcVals) = 0;

  // Update TACS states with the supplied ones (q, qdot, qddot)
  void setTACSStates( double time, BVec *q, BVec *qdot, BVec *qddot );

  // Add up the function contribution from each time step
  //-----------------------------------------------------
  void addFunctionContribution( double scale );

  // Add up the total derivative after solving for adjoint variables
  // ----------------------------------------------------------------
  void addTotalDerivative( double scale, BVec **adjoint );
  
  //Method to solve the non-linear system using Newton's method
  //------------------------------------------------------------
  void newtonSolve( double alpha, double beta, double gamma,
                    double t, BVec *q, BVec *qdot, 
                    BVec *qddot );

  // Calls LAPACK for the the solution of the linear system Ax =
  // b. Serial execution mode only.
  //-------------------------------------------------------------
  void lapackLinearSolve( BVec *res, TACSMat *mat, BVec *update );
  
  // Sanity checks on the RHS of the adjoint linear system
  //------------------------------------------------------
  static void checkAdjointRHS( BVec *rhs );
    
  // Get a formatted filename based on the current time step for
  // tecplot output
  //-----------------------------------------------------------
  static void getString(char *buffer, const char * format, ... );

  //-------------------------------------------------------------------------------------//
  //                                 Variables
  //-------------------------------------------------------------------------------------//

  // Instance of TACS
  TACSAssembler *tacs; 
  
  // Store the history of states over time
  BVec **q, **qdot, **qddot;
  double *time;

  // Number of state variables
  int num_state_vars;
  
  // Number of design variables
  int num_design_vars;
  
  // Class variables used to manage/monitor time marching in
  // integration schemes
  int num_time_steps, num_steps_per_sec; 
  double h, tinit, tfinal;

  // Print and output options
  int print_level;

  // Variables controlling the nonlinear solution
  int max_newton_iters;
  double atol, rtol;
  
  // Frequency of Jacobian recomputation during nonlinear solve
  int jac_comp_freq;

  // Flag to switch to LAPACK for linear solve
  int use_lapack;

  // Matrices and vectors for the nonlinear solution
  BVec *res, *update;  // Residual and Newton update
  FEMat *D;            // Matrix associated with Preconditioner
  TACSMat *mat;        // Jacobian matrix
  TACSPc *pc;          // Preconditioner
  TACSKsm *ksm;        // KSM solver

  // The objective and contraint functions
  int num_funcs;
  TACSFunction **funcs;

  // Pointers to the memory for storing function and derivative values
  TacsScalar *fvals;
  TacsScalar *dfdx;

  // Number of right hand sides we store during adjoint solve
  int num_adjoint_rhs; 

  // Tecplot output related variables
  TACSToFH5 *f5;
  char* f5_file_fmt;
  int f5_write_freq;
};

/*
  DIRK integration scheme for TACS
*/
class TACSDIRKIntegrator : public TACSIntegrator {
 public:
  // Constructor for DIRK object
  //----------------------------
  TACSDIRKIntegrator( TACSAssembler * _tacs, 
                      double _tinit, double _tfinal, int _num_steps_per_sec, 
                      int num_stages );
  
  // Desctructor
  //------------
  ~TACSDIRKIntegrator();

  // Integrate forward in time
  //--------------------------
  TacsScalar forward( const TacsScalar *x, int num_design_vars,
                      TACSFunction *func );

  // March backwards and find total derivative
  //------------------------------------------
  void reverse( TacsScalar *dfdx, int num_design_vars,
                TACSFunction *func );
  
  // Function to call to integrate in time
  //--------------------------------------
  void integrate();

 protected:
  // Evaluate time average of the function value using discretization
  // from the integration scheme
  //---------------------------------------------------------------------
  void evalTimeAvgFunctions( TACSFunction **funcs, int numFuncs, TacsScalar *funcVals);

 private:  
  // the order of accuracy of the scheme
  int num_stages, order;
  
  // variables for Butcher tableau
  double *A, *B, *C;
  
  // stage values (computed at each time stage for each time step)
  double *tS;
  BVec **qS, **qdotS, **qddotS;

  // Functions related to Butcher Tableau
  //-------------------------------------
  void setupButcherTableau();
  void checkButcherTableau();

  // Returns the starting index of corresponding Butcher tableau row
  //----------------------------------------------------------------
  int getRowIdx( int stageNum );

  // Approximate derivatives using DIRK formulae
  //--------------------------------------------
  void approxStates( int current_step, int current_stage );

  // Advance the time and states to next step
  //-----------------------------------------
  void computeTimeStepStates(int k,  BVec **q, BVec **qdot, BVec **qddot);

  // Function for marching backwards in stage and time
  //---------------------------------------------------
  void marchBackwards();
  
  // Get the coefficients for adding inter-stage contributions during adjoint solve
  //-------------------------------------------------------------------------------
  void getCoeffsInterStage( int i, int j, double *alpha, double *beta, double *gamma );
};

/*
  BDF integration scheme for TACS which extends TACSIntegrator
*/

class TACSBDFIntegrator : public TACSIntegrator {
 public:
  // Constructor for BDF object
  //---------------------------
  TACSBDFIntegrator(TACSAssembler * _tacs, 
		    double _tinit, double _tfinal, int _num_steps_per_sec, 
		    int max_bdf_order);
  
  // Destructor for BDF object
  //--------------------------
  ~TACSBDFIntegrator();
  
  // Integrate forward in time and find the solution
  //-------------------------------------------------
  TacsScalar forward( const TacsScalar *x, int num_design_vars,
                      TACSFunction *func );

  // March backwards in time and find the derivatives
  //-------------------------------------------------
  void reverse( TacsScalar *dfdx, int num_design_vars,
                TACSFunction *func );
  
  // Function that integrates forward in time
  //-----------------------------------------
  void integrate();
  
 protected:
  // Evaluate time average of the function value using discretization
  // from the integration scheme
  //---------------------------------------------------------------------
  void evalTimeAvgFunctions( TACSFunction **funcs, int numFuncs, TacsScalar *funcVals);
  
 private:  
  // Maximum order of the BDF integration scheme
  int max_bdf_order;

  // Number of first and second order BDF coefficients
  int nbdf, nbddf;

  // Class variable to store BDF coefficients
  double bdf_coeff[4], bddf_coeff[9];

  // Retrieve the first order BDF coefficients
  //------------------------------------------
  int getBDFCoeff(double bdf[], int order );

  // Retrieve the second order BDF coefficients
  //------------------------------------------
  void get2ndBDFCoeff(const int k, double bdf[], int *nbdf,
		      double bddf[], int *nbddf,
		      const int max_order);
  
  // approximate derivatives using BDF stencil
  //------------------------------------------
  void approxStates( int current_step );

  // Function for marching backwards in stage and time
  //---------------------------------------------------
  void marchBackwards();
};

/*
  Adams-Bashforth-Moulton integration scheme for TACS
*/
class TACSABMIntegrator : public TACSIntegrator {
 public:
  // Constructor for ABM object
  //---------------------------
  TACSABMIntegrator(TACSAssembler * _tacs, 
		    double _tinit, double _tfinal, int _num_steps_per_sec, 
		    int max_abm_order);
  
  // Destructor for ABM object
  //--------------------------
  ~TACSABMIntegrator();
  
  // Function that integrates forward in time
  //-----------------------------------------
  void integrate();
  
 private:   
  // Start index of the the row of ABM coefficient matrix.
  // -----------------------------------------------------
  int getRowIdx( int stageNum );

  // Order of integration to use for the current time step k
  //--------------------------------------------------------
  int getOrder( int k );

  // Create a table of coeffs for ABM
  void setupABMCoeffs( int max_order, double *A);
  void checkABMCoeffs();

  // Approximate derivatives using ABM stencil
  //------------------------------------------
  void approxStates( int current_step, int current_order );
  
  // Evaluate time average of the function value using discretization
  // from the integration scheme
  //---------------------------------------------------------------------
  void evalTimeAvgFunctions( TACSFunction **funcs, int numFuncs, TacsScalar *funcVals);

  // Function for marching backwards in stage and time
  //---------------------------------------------------
  void marchBackwards();
  
  // Maximum order of the ABM integration scheme
  int max_abm_order;

  // ABM Coefficients
  double *A;
};

/*
  Newmark Beta Gamma Method (NBG) integration scheme for TACS
*/
class TACSNBGIntegrator : public TACSIntegrator {
 public:
  // Constructor for NBG object
  //---------------------------
  TACSNBGIntegrator(TACSAssembler * _tacs, 
		    double _tinit, double _tfinal, 
                    int _num_steps_per_sec);
  
  // Destructor for NBG object
  //--------------------------
  ~TACSNBGIntegrator();
  
  // Function that integrates forward in time
  //-----------------------------------------
  void integrate();
  
 private:
  // Approximate derivatives using NBG formula
  //------------------------------------------
  void approxStates( int current_step );
  
  // Evaluate time average of the function value using discretization
  // from the integration scheme
  //---------------------------------------------------------------------
  void evalTimeAvgFunctions( TACSFunction **funcs, int numFuncs, TacsScalar *funcVals);

  // Function for marching backwards in stage and time
  //---------------------------------------------------
  void marchBackwards();
  
  // Newmark Coefficients
  static const double BETA  = 0.25;
  static const double GAMMA = 0.50;
};

#endif
