#ifndef TACS_INTEGRATOR_H
#define TACS_INTEGRATOR_H

#include "TACSObject.h"
#include "TACSAssembler.h"
#include "KSM.h"
#include "TACSToFH5.h"

// Type of integrator to use. The following are the supported methods.
//--------------------------------------------------------------------
enum IntegratorType { 
  BDF1, BDF2, BDF3,                   // Backward-difference methods
  ABM1, ABM2, ABM3, ABM4, ABM5, ABM6, // Adams-Bashforth-method
  DIRK2, DIRK3, DIRK4,                // Diagonally-Implicit-Runge-Kutta
  NBGE, NBG2, NBG3 };                 // Newmark Beta Gamma Method

/*
  Abstract base class for integration schemes.
  
  This base class contains common functions and variables pertaining
  to the integration and adjoint schemes used in TACS.

  Certain functions are pure virtual and therefore the integrator can
  be instantiated only after those functions are provided a full
  implementation in an extending child class.

  Copyright (c) 2010-2016 Graeme Kennedy. All rights reserved.
*/


class TACSIntegrator : public TACSObject {
 public:
  TACSIntegrator( TACSAssembler *tacs,
                  double tinit,
                  double tfinal, 
                  double num_steps_per_sec ); 
  ~TACSIntegrator();

  // Set the functions to integrate
  //--------------------------------
  void setFunctions( TACSFunction **funcs, int num_funcs,
                     int start_step=-1, int end_step=-1 );

  // Initialize the TACSIntegrator object
  //-------------------------------------
  void initialize();
  
  // Solve for time-step t with 
  virtual void iterate( int step_num, TACSBVec *forces=NULL ) = 0;

  // Integrate the equations of motion forward in time
  virtual void integrate(){
    for ( int i = 0; i < num_steps; i++ ){
      iterate(i);
    }
  }

  // Set-up right-hand-sides for the adjoint equations
  virtual void initAdjoint( int step_num ) = 0;

  // Iterate to find a solution of the adjoint equations
  virtual void iterateAdjoint( int step_num, TACSBVec *psi=NULL ) = 0;

  // Add the contributions to the total derivative from the time-step 
  virutal void addAdjointDerivative( int step_num ) = 0;

  // Integrate the adjoint and add the total derivative from all
  // time-steps
  virtual void integrateAdjoint(){
    for ( int i = num_steps-1; i >= 0; i-- ){
      initAdjoint(i);
      iterateAdjoint(i);
      addAdjointDerivative(i);
    }
  }

  // Get the adjoint vector
  //-----------------------
  void getAdjoint( int step_num, TACSBVec **adjoint );
  virtual void evalFunctions( TacsScalar *fvals );
  virtual void getGradient( TacsScalar *dfdx );

  // Retrieve the internal states
  //------------------------------
  void getStates( int step_num, 
                  TACSBVec **q, TACSBVec **qdot, TACSBVec **qddot );

  // Define members of the class
  // ---------------------------
  void setRelTol( double _rtol );
  void setAbsTol( double _atol );
  void setMaxNewtonIters( int _max_newton_iters );
  void setPrintLevel( int _print_level, const char *logfilename=NULL );
  void setJacAssemblyFreq( int _jac_comp_freq );
  void setUseLapack( int _use_lapack, int _eigensolve=0 );
  void setUseLineSearch( int _use_line_search );
  void setUseFEMat( int _use_femat );
  void setInitNewtonDeltaFraction( double frac );

  // Functions to export the solution in raw and tecplot binary forms
  //-----------------------------------------------------------------
  void writeSolution( const char *filename, int format=2 );
  void writeSolutionToF5( int _write_step );
  void writeStepToF5( int step );
  void writeNewtonIterToF5( int step, int newton );
  void setOutputFrequency( int _write_step, int _write_newton = 0);
  void setRigidOutput( int flag );    
  void setShellOutput( int flag );
  void setBeamOutput( int flag );
  void printWallTime( double t0, int level=1 );
  void printOptionSummary( FILE *fp );
  void printAdjointOptionSummary( FILE *fp );

 protected:

  // Functions for solutions to linear and nonlinear problems
  // --------------------------------------------------------
  int newtonSolve( double alpha, double beta, double gamma,
                   double t, TACSBVec *q, TACSBVec *qdot, TACSBVec *qddot, 
                   TACSBVec *forces=NULL, TACSBcMap *addBcs=NULL );
  void lapackLinearSolve( TACSBVec *res, TACSMat *mat, TACSBVec *update );
  void lapackEigenSolve( TACSMat *mat );

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
  
  // Functions to aid adjoint solve
  // -------------------------------
  void addVectorTransProducts( TACSBVec **ans, 
                               double alpha, double beta, double gamma,
                               int num_funcs, TACSFunction **funcs,
                               TACSBVec **input );

  // // Virtual functions for forward mode
  // // -----------------------------------
  // virtual void setupCoeffs() = 0;
  // virtual void approxStates() = 0;
  // virtual void getLinearizationCoeffs( double *alpha, double *beta, 
  //                                      double *gamma ) = 0;
  
  // // Virtual functions for reverse mode
  // // ----------------------------------
  // virtual void marchBackwards() = 0;
  // virtual void evalFunctions( TACSFunction **funcs, TacsScalar *funcVals,
  //                             int numFuncs );  

  // These functions are to be used during adjoint solve and total
  // derivative computations. Once the state variables are set, the
  // current state can not only beq used to find the adjoint variables
  // at the current step but also used to compute as much as possible
  // with the current snapshot of tacs(t, q, qdot, qddot). These
  // possible operations are packed into each one of the following
  // functions.
  //----------------------------------------------------------------
  /*
  void addFuncSVSensToAdjointRHS( double alpha, double beta, double gamma, int index ) ;
  void addSVSensToAdjointRHS( double scale, double alpha, double beta, double gamma,
                              int adj_index, TACSFunction **funcs,
                              int num_funcs, TacsScalar *fvals );
  void adjointSolve( int adj_index, TACSBVec **adjoint );
  */
  // void addToTotalDerivative( double scale, TACSBVec **adjoint );
    
  // Functions to pack common logics (file io, output) in one place
  // and use them in appropriate places
  // ------------------------------------------------------------------
  void doEachTimeStep( int current_step );
  void doEachNonLinearIter( int iter_num );
  int  doAdaptiveMarching();
                    
  //------------------------------------------------------------------//
  //                  Protected variables
  //------------------------------------------------------------------//

  TACSAssembler *tacs;        // Instance of TACS
  enum IntegratorType mytype; // Will be set the by child class
  FILE *logfp;                // Pointer to the output filename

  TACSBVec *forces;           // forces applied to the RHS
  TACSBVec **q;               // state variables across all time steps
  TACSBVec **qdot;            // first time derivative of ''
  TACSBVec **qddot;           // second time derivative of ''
  
  double *time;               // Stores the time values
  int num_time_steps;         // Total number of time steps to take to go from initial time to final time
  double num_steps_per_sec;   // Number of time steps to take to advance one second
  int current_time_step;      // Tracks the current time step
  double h;                   // Time step size
  double tinit;               // Initial simulation time
  double tfinal;              // Final simulation time
  
  TACSFunction **funcs;       // List of functions
  int num_funcs;              // The number of objective functions
  TacsScalar *fvals;          // Function values
  TacsScalar *dfdx;           // Derivative values
  int num_state_vars;         // Number of state variables
  int num_design_vars;        // Number of design variables

  int newton_term;            // Termination of nonlinear solver 
                              // 1: |R| < atol; 2: |dq| < atol
                              // 3: |R|/|R0| < rtol
                              // -1: max_newton_iters // -2: Nan
  TACSMat *mat;               // Jacobian matrix
  TACSPc *pc;                 // Preconditioner
  TACSBVec *res, *update;     // Residual and Newton update
  TACSKsm *ksm;               // KSM solver

  int num_adjoint_rhs;        // the number of right hand sides allocated in adjoint loop  
  TACSBVec **rhs;             // storage vector for the right hand sides
  TACSBVec **psi;             // adjoint variable accumulating qdot dependance
  TACSBVec **phi;             // adjoint variable accumulating q dependance
  TACSBVec **lambda;          // adjoint variable qddot
  TACSBVec **dfdq;            // storage vector for statevariable sensitivities

 private:

  static void getString( char *buffer, const char * format, ... );
  int getWriteFlag( int step, int freq );

  int print_level;          // 0 = off;
                            // 1 = summary per time step;
                            // 2 = summary per Newton solve iteration
  TACSToFH5 *rigidf5;       // F5 file for rigid body visualization
  TACSToFH5 *shellf5;       // F5 file for shell visualization
  TACSToFH5 *beamf5;        // F5 file for beam visualization
  int f5_write_freq;        // How frequent to write the output during time marching
  int f5_newton_freq;       // How frequent to write the output during nonlinear solve
  
  int max_newton_iters;     // The max number of nonlinear iterations
  double atol;              // Absolute tolerance
  double rtol;              // Relative tolerance
  int jac_comp_freq;        // Frequency of Jacobian factorization
  int use_line_search;      // Flag to make use of line search for nonlinear root finding

  int use_lapack;           // Flag to switch to LAPACK for linear solve
  int use_femat;            // use femet for parallel execution
  int eigensolve;           // Solve for eigenvalues
  int lev;
  int fill;
  int reorder_schur;
  int gmres_iters;
  int num_restarts;
  int is_flexible;

  FEMat *D;                 // Matrix associated with Preconditioner
  int factorized;           // Set whether the matrix is factorized
  int niter;                // Newton iteration number
  TacsScalar res_norm;      // residual norm
  TacsScalar init_res_norm; // Initial norm of the residual
  TacsScalar update_norm;   // Norm of the update                            
  double init_newton_delta; // Initial value of delta used in the globalization strategy
 
  TacsScalar energies[2];   // Keep track of energies
  TacsScalar init_energy;   // The energy during time = 0
  int mpiRank;              // rank of the processor
  int mpiSize;              // number of processors
};

/*
  Backwards Difference Formulas integration scheme
*/
class TACSBDFIntegrator : public TACSIntegrator {
 public:
  TACSBDFIntegrator( TACSAssembler * _tacs, 
		     double _tinit,
                     double _tfinal,
                     double _num_steps_per_sec, 
		     int max_bdf_order );

 protected:
  void iterate( int k, TACSBVec *forces=NULL );

 private:  
  void getLinearizationCoeffs( double *alpha, double *beta, double *gamma );
  void get2ndBDFCoeff( const int k, double bdf[], int *nbdf,
		       double bddf[], int *nbddf,
		       const int max_order );
  int getBDFCoeff( double bdf[], int order );

  int max_bdf_order;    // Maximum order of the BDF integration scheme
  int nbdf, nbddf;      // Number of first and second order BDF coefficients
  double bdf_coeff[4];  // Store first order BDF coefficients
  double bddf_coeff[9]; // Store second order BDF coefficients
};

#endif
