#ifndef TACS_INTEGRATOR_H
#define TACS_INTEGRATOR_H

#include "TACSObject.h"
#include "TACSAssembler.h"
#include "KSM.h"
#include "TACSToFH5.h"

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

  // Define members of the class
  // ---------------------------
  void setRelTol( double _rtol );
  void setAbsTol( double _atol );
  void setMaxNewtonIters( int _max_newton_iters );
  void setPrintLevel( int _print_level, const char *logfilename=NULL );
  void setJacAssemblyFreq( int _jac_comp_freq );
  void setUseLapack( int _use_lapack );
  void setUseFEMat( int _use_femat, TACSAssembler::OrderingType _type );
  void setInitNewtonDeltaFraction( double frac );

  // Set the functions to integrate
  //--------------------------------
  void setFunctions( TACSFunction **funcs, int num_funcs,
                     int num_design_vars,
                     int start_step=-1, int end_step=-1 );
  
  // Solve for time-step t with 
  virtual int iterate( int step_num, TACSBVec *forces ) = 0;

  // Integrate the equations of motion forward in time
  virtual void integrate(){
    for ( int i = 0; i < num_time_steps+1; i++ ){
      iterate(i, NULL);
    }
  }

  // Evaluate the functions of interest
  virtual void evalFunctions( TacsScalar *fvals ) = 0;

  // Set-up right-hand-sides for the adjoint equations
  virtual void initAdjoint( int step_num ){}

  // Iterate to find a solution of the adjoint equations
  virtual void iterateAdjoint( int step_num, TACSBVec **adj_rhs ){}

  // Add the contributions to the total derivative from the time-step 
  virtual void postAdjoint( int step_num ){}

  // Integrate the adjoint and add the total derivative from all
  // time-steps
  virtual void integrateAdjoint(){
    for ( int i = num_time_steps; i >= 0; i-- ){
      initAdjoint(i);
      iterateAdjoint(i, NULL);
      postAdjoint(i);
    }
  }

  // Get the adjoint vector for the given function
  virtual void getAdjoint( int step_num, int func_num, 
                           TACSBVec **adjoint ) = 0;

  // Copy out the function
  virtual void getGradient( TacsScalar *_dfdx ){
    memcpy(_dfdx, dfdx, num_funcs*num_design_vars*sizeof(TacsScalar));
  }

  // Retrieve the internal states
  double getStates( int step_num, 
                    TACSBVec *q, TACSBVec *qdot, TACSBVec *qddot );

  // Check the gradient using finite-difference
  void checkGradients( double dh );

  // Functions to export the solution in raw and tecplot binary forms
  //-----------------------------------------------------------------
  void setOutputPrefix( const char *prefix );
  void setOutputFrequency( int _write_step );
  void setRigidOutput( TACSToFH5 *_rigidf5 );
  void setShellOutput( TACSToFH5 *_shellf5 );
  void setBeamOutput( TACSToFH5 *_beamf5 );
  void writeSolution( const char *filename, int format=2 );
  void writeSolutionToF5( int step_num );
  void printWallTime( double t0, int level=1 );
  void printOptionSummary();
  void printAdjointOptionSummary();

 protected:
  // Functions for solutions to linear and nonlinear problems
  int newtonSolve( double alpha, double beta, double gamma,
                   double t, TACSBVec *q, TACSBVec *qdot, TACSBVec *qddot, 
                   TACSBVec *forces=NULL );
  void lapackLinearSolve( TACSBVec *res, TACSMat *mat, TACSBVec *update );

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
  void logTimeStep( int time_step );

  // TACSAssembler information
  TACSAssembler *tacs;        // Instance of TACS

  // The step information
  int num_time_steps;         // Total number of time steps
  double *time;               // Stores the time values
  TACSBVec **q;               // state variables across all time steps
  TACSBVec **qdot;            // first time derivative of ''
  TACSBVec **qddot;           // second time derivative of ''

  // Objects that store information about the functions of interest
  int start_step, end_step;   // Time-window for the functions of interest
  TACSFunction **funcs;       // List of functions
  int num_funcs;              // The number of objective functions
  TacsScalar *fvals;          // Function values
  TacsScalar *dfdx;           // Derivative values
  int num_design_vars;        // Number of design variables

  // Linear algebra objects and parameters associated with the Newton solver  
  TACSBVec *res, *update;     // Residual and Newton update
  TACSMat *mat;               // Jacobian matrix
  TACSPc *pc;                 // Preconditioner
  TACSKsm *ksm;               // KSM solver

 private:
  char prefix[256];           // Output prefix
  FILE *logfp;                // Pointer to the output filename

  // Newton solver parameters
  int max_newton_iters;     // The max number of nonlinear iterations
  double atol;              // Absolute tolerance
  double rtol;              // Relative tolerance
  double init_newton_delta; // Initial value of delta for globalization
  int jac_comp_freq;        // Frequency of Jacobian factorization
  int use_femat;            // use femet for parallel execution
  TACSAssembler::OrderingType order_type;
  int use_lapack;           // Flag to switch to LAPACK for linear solve

  int lev;
  double fill;
  int reorder_schur;
  int gmres_iters;
  int num_restarts;
  int is_flexible;

  // Information for visualization/logging purposes
  int print_level;          // 0 = off;
                            // 1 = summary per time step;
                            // 2 = summary per Newton solve iteration
  TACSToFH5 *rigidf5;       // F5 file for rigid body visualization
  TACSToFH5 *shellf5;       // F5 file for shell visualization
  TACSToFH5 *beamf5;        // F5 file for beam visualization
  int f5_write_freq;        // Frequency for output during time marching
  
  int niter;                // Newton iteration number
  TacsScalar res_norm;      // residual norm
  TacsScalar init_res_norm; // Initial norm of the residual
  TacsScalar update_norm;   // Norm of the update                            
 
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
  ~TACSBDFIntegrator();

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
  void get2ndBDFCoeff( const int k, double bdf[], int *nbdf,
		       double bddf[], int *nbddf,
		       const int max_order );
  int getBDFCoeff( double bdf[], int order );

  int max_bdf_order;    // Maximum order of the BDF integration scheme

  // Adjoint information
  int num_adjoint_rhs;  // the number of right hand sides allocated
  TACSBVec **rhs;  // storage vector for the right hand sides
  TACSBVec **psi;  // adjoint variable accumulating qdot dependance
};

/*
  DIRK integration scheme for TACS
*/
/*
class TACSDIRKIntegrator : public TACSIntegrator {
 public:
  TACSDIRKIntegrator( TACSAssembler * _tacs, 
                      double _tinit,
                      double _tfinal,
                      double _num_steps_per_sec, 
                      int _order );
  ~TACSDIRKIntegrator();

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
  void checkButcherTableau();
  int getRowIdx( int stageNum );

  double    *tS;                    // Time at each stage
  TACSBVec **qS, **qdotS, **qddotS; // States at each stage
  int        num_stages, order;     // The order of accuracy of the scheme
  double    *A, *B, *C;             // Variables for Butcher tableau
  int        current_stage;         // The current stage during iteration
};
*/

/*
  Adams-Bashforth-Moulton integration scheme for TACS
*/
/*
class TACSABMIntegrator : public TACSIntegrator {
 public:
  TACSABMIntegrator( TACSAssembler *tacs, 
		     double tinit,
                     double tfinal,
                     double num_steps_per_sec, 
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
                     double num_steps_per_sec,
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
#endif // TACS_INTEGRATOR_H
