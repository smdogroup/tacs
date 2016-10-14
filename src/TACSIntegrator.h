#ifndef TACS_INTEGRATOR_H
#define TACS_INTEGRATOR_H

#include <stdarg.h>
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
  NBG };                              // Newmark Beta Gamma Method

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
                  int num_steps_per_sec ); 
  ~TACSIntegrator();

  // Factory method for intantiating integrators of supplied type
  //-------------------------------------------------------------
  static TACSIntegrator *getInstance( TACSAssembler *_tacs, 
                                      double _tinit,
                                      double _tfinal, 
                                      int _num_steps_per_sec, 
                                      enum IntegratorType type );
  
  // Function to solve for the states in time
  //-----------------------------------------
  virtual void integrate();

  // Function for returning the derivatives for the functions
  //----------------------------------------------------------
  void getFuncGrad( int num_dv, TacsScalar *x, 
                    TacsScalar *fvals, TacsScalar *dfdx );
  void getFDFuncGrad( int num_dv, TacsScalar *x, 
                      TacsScalar *fvals, TacsScalar *dfdx, 
                      double dh=1.0e-8);
 
  // Useful setters for class variables. These setters must be called
  // before calling the 'integrate()' function. These won't have any
  // effect after the integrate call has been made as default values
  // would have been used during time integration.
  // -----------------------------------------------------------------
  void setParameters( ... );
  void setRelTol( double _rtol );
  void setAbsTol( double _atol );
  void setMaxNewtonIters( int _max_newton_iters );
  void setPrintLevel( int _print_level, const char *logfilename=NULL );
  void setJacAssemblyFreq( int _jac_comp_freq );
  void setUseLapack( int _use_lapack );
  void setFunction( TACSFunction **_func, int _num_funcs );
  void setIsFactorized( int flag );
  void setTACSStates( double time, TACSBVec *q, 
                      TACSBVec *qdot, TACSBVec *qddot );
  void setOrderingType( TACSAssembler::OrderingType _type );

  // Functions to export the solution in raw and tecplot binary forms
  //-----------------------------------------------------------------
  void writeSolution( const char *filename );
  void writeSolutionToF5( TACSToFH5 *f5, int _write_freq, 
                          const char *_f5_file_fmt );
  void writeStepToF5( int step = 0 );          
  void configureOutput( TACSToFH5 *_viewer, int _write_freq, 
                        const char *_f5_file_fmt );
  void printWallTime( double t0, int level=1 );

  // Variables that keep track of time
  double time_fwd_assembly;
  double time_fwd_factor;
  double time_fwd_apply_factor;
  double time_forward;

  double time_rev_assembly;
  double time_rev_factor;
  double time_rev_apply_factor;
  double time_rev_jac_pdt;
  double time_reverse;
 protected:

  //------------------------------------------------------------------//
  //                  Protected functions
  //------------------------------------------------------------------//

  // Functions for solutions to linear and nonlinear problems
  // --------------------------------------------------------
  void newtonSolve( double alpha, double beta, double gamma,
                    double t, TACSBVec *q, TACSBVec *qdot, 
                    TACSBVec *qddot );
  void lapackLinearSolve( TACSBVec *res, TACSMat *mat, TACSBVec *update );
    
  // Virtual functions for forward mode
  // -----------------------------------
  virtual void setupCoeffs() = 0;
  virtual void approxStates() = 0;
  virtual void getLinearizationCoeffs( double *alpha, double *beta, 
                                       double *gamma ) = 0;
  
  // Virtual functions for reverse mode
  // ----------------------------------
  virtual void marchBackwards() = 0;
  virtual void evalTimeAvgFunctions( TACSFunction **funcs, 
                                     int numFuncs, TacsScalar *funcVals);
  
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
  void addToTotalDerivative( double scale, TACSBVec **adjoint );
  void addFunctions( int time_step, double tcoeff, TACSFunction **funcs, 
                     int numFuncs, TacsScalar *funcVals);
    
  // Functions to pack common logics (file io, output) in one place
  // and use them in appropriate places
  // ------------------------------------------------------------------
  void doEachTimeStep( int current_step );
  void doEachNonLinearIter( int iter_num );
                    
  //------------------------------------------------------------------//
  //                  Protected variables
  //------------------------------------------------------------------//
  
  TACSAssembler  *tacs;               // Instance of TACS
  TACSBVec      **q, **qdot, **qddot; // Store the history of states over time
  
  double         *time;               // Stores the time values
  int             num_time_steps;     // Total number of time steps to take to go from initial time to final time
  int             num_steps_per_sec;  // Number of time steps to take to advance one second
  int             current_time_step;  // Tracks the current time step
  double          h, tinit, tfinal;   // Time step size, initial and final times
  
  TACSFunction  **funcs;              // List of objective functions
  int             num_funcs;          // The number of objective functions
  TacsScalar     *fvals;              // Function values (returned at the end of adjoint solve)
  TacsScalar     *dfdx;               // Derivative values (returned at the end of adjoint solve)

  int             num_state_vars;     // Number of state variables
  int             num_design_vars;    // Number of design variables
  
  TACSMat        *mat;                // Jacobian matrix
  TACSPc         *pc;                 // Preconditioner
  TACSBVec       *res, *update;       // Residual and Newton update
  TACSKsm        *ksm;                // KSM solver
  int             num_adjoint_rhs;    // the number of right hand sides allocated in adjoint loop
  
  TACSBVec      **rhs;                // storage vector for the right hand sides
  TACSBVec      **psi;                // adjoint variable accumulating qdot dependance
  TACSBVec      **phi;                // adjoint variable accumulating q dependance
  TACSBVec      **lambda;             // adjoint variable qddot
  TACSBVec      **dfdq;               // storage vector for statevariable sensitivities
 private:

  //-----------------------------------------------------------------//
  //                   Private functions
  //-----------------------------------------------------------------//
  
  static void getString( char *buffer, const char * format, ... );
  int getWriteFlag( int step, int f5_write_freq );

  //-----------------------------------------------------------------//
  //                   Private variables
  //-----------------------------------------------------------------//

  int print_level;      // 0 = off;
                        // 1 = summary per time step;
                        // 2 = summary per Newton solve iteration
  FILE *logfp;          // Pointer to the output filename
  TACSToFH5 *f5;        // Pointer to FH5 output object
  char *f5_file_fmt;    // Formatting of the output filename
  int f5_write_freq;    // How frequent to write the output
  
  int max_newton_iters; // The max number of nonlinear iterations
  double atol, rtol;
  
  int jac_comp_freq;    // Frequency of Jacobian factorization
  int use_lapack;       // Flag to switch to LAPACK for linear solve
  FEMat *D;             // Matrix associated with Preconditioner
  int factorized;       // Set whether the matrix is factorized
  int niter;            // Newton iteration number
  TacsScalar norm, init_norm, update_norm;  // Norms and initial norm
  
  TacsScalar energies[2]; // Keep track of energies
  TacsScalar init_energy; // The energy during time = 0
  int mpiRank, mpiSize;   // MPI information

  TACSAssembler::OrderingType ordering_type; // The type of ordering to use
};

/*
  DIRK integration scheme for TACS
*/
class TACSDIRKIntegrator : public TACSIntegrator {
 public:
  TACSDIRKIntegrator( TACSAssembler * _tacs, 
                      double _tinit,
                      double _tfinal,
                      int _num_steps_per_sec, 
                      int num_stages );
  ~TACSDIRKIntegrator();

  // Overriding the default integration logic
  void integrate();
  
 protected:
  void setupCoeffs();
  void approxStates();
  void getLinearizationCoeffs( double *alpha, double *beta, double *gamma );
  
  void marchBackwards();

  // Overridden procedure
  void evalTimeAvgFunctions( TACSFunction **funcs, int numFuncs, TacsScalar *funcVals);

 private:
  void checkButcherTableau();
  int getRowIdx( int stageNum );
  void computeTimeStepStates( int k,  TACSBVec **q, TACSBVec **qdot, TACSBVec **qddot );
  void getCoeffsInterStage( int i, int j, double *alpha, double *beta, double *gamma );
  
  double    *tS;                    // Time at each stage
  TACSBVec **qS, **qdotS, **qddotS; // States at each stage
  int        num_stages, order;     // The order of accuracy of the scheme
  double    *A, *B, *C;             // Variables for Butcher tableau
  int        current_stage;         // The current stage during iteration
};

/*
  Backwards Difference Formulas integration scheme
*/
class TACSBDFIntegrator : public TACSIntegrator {
 public:
  TACSBDFIntegrator( TACSAssembler * _tacs, 
		     double _tinit,
                     double _tfinal,
                     int _num_steps_per_sec, 
		     int max_bdf_order );
  ~TACSBDFIntegrator();
 protected:
  void setupCoeffs();
  void approxStates();
  void getLinearizationCoeffs( double *alpha, double *beta, double *gamma );
  void marchBackwards();
 private:  
  int getBDFCoeff( double bdf[], int order );
  void get2ndBDFCoeff( const int k, double bdf[], int *nbdf,
		       double bddf[], int *nbddf,
		       const int max_order );

  int    max_bdf_order; // Maximum order of the BDF integration scheme
  int    nbdf, nbddf;   // Number of first and second order BDF coefficients
  double bdf_coeff[4];  // Store first order BDF coefficients
  double bddf_coeff[9]; // Store second order BDF coefficients
};

/*
  Adams-Bashforth-Moulton integration scheme for TACS
*/
class TACSABMIntegrator : public TACSIntegrator {
 public:
  TACSABMIntegrator( TACSAssembler *tacs, 
		     double         tinit,
                     double         tfinal,
                     int            num_steps_per_sec, 
		     int            max_abm_order );
  ~TACSABMIntegrator();
 protected:
  void setupCoeffs();
  void approxStates();
  void getLinearizationCoeffs( double *alpha, double *beta, double *gamma );
  void marchBackwards();
 private:
  int getRowIdx( int row_num );
  int getOrder( int k );
  void checkABMCoeffs();

  int max_abm_order; // Order of integration
  double *A;         // Array of ABM coefficients
};

/*
  Newmark Beta Gamma Method (NBG) integration scheme for TACS
*/
class TACSNBGIntegrator : public TACSIntegrator {
 public:
  TACSNBGIntegrator( TACSAssembler *tacs, 
		     double         tinit,
                     double         tfinal, 
                     int            num_steps_per_sec );
  ~TACSNBGIntegrator();
  
 protected:
  void setupCoeffs();
  void approxStates();
  void getLinearizationCoeffs( double *alpha, double *beta, double *gamma );
  void marchBackwards();
 private:
  double BETA, GAMMA; // Newmark coefficients
};

#endif
