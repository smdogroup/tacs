#ifndef TACS_OSCILLATOR_ELEMENT_H
#define TACS_OSCILLATOR_ELEMENT_H

#include "TACSElement.h"

/*
  A test implementation of TACSElement class with a two-DOF aero
  elastic oscillator system. The degrees of freedom are pitching and
  plunging states.

  Reference:
  Zhao, L. and Yang, Z., ``Chaotic motions of an airfoil with
  non-linear stiffness in incompressible flow,'' Journal of Sound and
  Vibration, Vol. 138, No. 2, 1990, pp. 245–254.
*/
class TACSOscillatorElement : public TACSElement {
 public:
  TACSOscillatorElement(){
    Q = 8.0;
  }
  ~TACSOscillatorElement(){}
  
  // Retrieve information about the name and quantity of variables
  // -------------------------------------------------------------
  const char * elementName() { return  "TACSOscillatorODE"; } 
  const char * displacementName( int i ) {}
  const char * stressName( int i ) {}
  const char * strainName( int i ) {}
  const char * extraName( int i ) {}
  
  // Get the number of displacements, stresses, nodes, etc.
  // ------------------------------------------------------
  int numDisplacements() { return 2; } // should equal num_state_vars
  int numStresses() { return 0; }
  int numNodes() { return 1; }
  int numVariables() { return 2; } // should equal num_state_vars
  int numExtras(){ return 0; }
  enum ElementType getElementType(){return SHELL;}

  // Return the name of the element
  // ------------------------------
  const char * TACSObjectName(){ return this->elementName(); }

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies( double time,
			TacsScalar *_Te, 
			TacsScalar *_Pe,
			const TacsScalar Xpts[],
			const TacsScalar vars[],
			const TacsScalar dvars[] ) {}

  // Compute the Residual of the governing equations
  // -----------------------------------------------
  void addResidual( double time, TacsScalar res[],
		    const TacsScalar Xpts[],
		    const TacsScalar q[],
		    const TacsScalar qdot[],
		    const TacsScalar qddot[] ) {
    // First equation
    res[0] += qddot[0] + 0.25*qdot[1] + 0.1*qdot[0] 
      + 0.2*q[0] + 0.1*Q*q[1];

    // Second equation
    res[1] += 0.25*qddot[0] + 0.5*qddot[1] 
      + 0.1*qdot[1] + 0.5*q[1] + 20.0*q[1]*q[1]*q[1] - 0.1*Q*q[1];
  }

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian( double time, TacsScalar J[],
		    double alpha, double beta, double gamma,
		    const TacsScalar Xpts[],
		    const TacsScalar q[],
		    const TacsScalar qdot[],
		    const TacsScalar qddot[] ) {

    //-------------------------------------------------------------//
    // Add dR/dUDDOT
    //-------------------------------------------------------------//

    // derivative of first equation

    J[0] += gamma*1.0;
    J[1] += gamma*0.0;

    // derivative of second equation

    J[2] += gamma*0.25;
    J[3] += gamma*0.50;


    //-------------------------------------------------------------//
    // Add dR/dUDOT
    //-------------------------------------------------------------//

    // derivative of first equation

    J[0] += beta*0.10;
    J[1] += beta*0.25;

    // derivative of second equation

    J[2] += beta*0.0;
    J[3] += beta*0.1;


    //-------------------------------------------------------------//
    // Add dR/du
    //-------------------------------------------------------------//

    // derivative of first equation

    J[0] += alpha*0.2;
    J[1] += alpha*0.1*Q;

    // derivative of second equation

    J[2] += alpha*0.0;
    J[3] += alpha*(0.5 + 60.0*q[1]*q[1] - 0.1*Q);
  }

  // Retrieve the initial values of the state variables
  // --------------------------------------------------
  void getInitCondition( TacsScalar vars[],
			 TacsScalar dvars[],
			 const TacsScalar X[] ){
    vars[0] = 2.0;
    vars[1] = 1.0;

    dvars[0] = 0.2;
    dvars[1] = -0.1;
  }  

  // Add the product of the adjoint with the derivative of the design variables
  // --------------------------------------------------------------------------
  void addAdjResProduct( double scale,
			 TacsScalar dvSens[], int dvLen,
			 const TacsScalar psi[],
			 const TacsScalar Xpts[],
			 const TacsScalar vars[],
			 const TacsScalar dvars[],
			 const TacsScalar ddvars[] ){}

  // Add the product of the adjoint with the derivative of the design variables
  // --------------------------------------------------------------------------  
  void getAdjResXptProduct( TacsScalar XptSens[],
			    const TacsScalar psi[],
			    const TacsScalar Xpts[],
			    const TacsScalar vars[],
			    const TacsScalar dvars[],
			    const TacsScalar ddvars[] ){}
  
  // Retrieve a specific time-independent matrix from the element
  // ------------------------------------------------------------
  void getMatType( ElementMatrixType matType, 
		   TacsScalar mat[], 
		   const TacsScalar Xpts[],
		   const TacsScalar vars[] ) {}

  // Compute the derivative of the inner product w.r.t. design variables
  // -------------------------------------------------------------------
  void addMatDVSensInnerProduct( ElementMatrixType matType, 
				 double scale,
				 TacsScalar dvSens[], int dvLen,
				 const TacsScalar psi[], 
				 const TacsScalar phi[],
				 const TacsScalar Xpts[],
				 const TacsScalar vars[] ){}

  // Compute the derivative of the inner product w.r.t. vars[]
  // ---------------------------------------------------------
  void getMatSVSensInnerProduct( ElementMatrixType matType, 
				 TacsScalar res[],
				 const TacsScalar psi[], 
				 const TacsScalar phi[],
				 const TacsScalar Xpts[],
				 const TacsScalar vars[] ){}

  // Member functions for evaluating global functions of interest
  // ------------------------------------------------------------
  TACSConstitutive * getConstitutive(){ return NULL; }

  // Get the number of Gauss quadrature points
  // -----------------------------------------
  int getNumGaussPts(){ return 0; }
  
  // Get the quadrature points and weights
  // -------------------------------------
  double getGaussWtsPts( const int num, double * pt ) { return 0.0; }

  // Get the shape functions from the element
  // ----------------------------------------
  void getShapeFunctions( const double pt[], double N[] ) {}
  
  // Return the determinant of the Jacobian at this point
  // ----------------------------------------------------
  TacsScalar getDetJacobian( const double * pt, 
			     const TacsScalar Xpts[] ){ return 0.0; }
  
  // Return the determinant of the Jacobian and its sensitivity at this point
  // ------------------------------------------------------------------------
  TacsScalar getDetJacobianXptSens( TacsScalar * hXptSens, 
				    const double * pt, 
				    const TacsScalar Xpts[] ){ return 0.0; }

  // This function returns the strain evaluated at pt
  // ------------------------------------------------
  void getStrain( TacsScalar strain[], 
		  const double pt[], 
		  const TacsScalar Xpts[],
		  const TacsScalar vars[] ) {}

  // This function returns the sensitivity of the strain w.r.t. Xpts
  // ---------------------------------------------------------------
  void addStrainXptSens( TacsScalar strainXptSens[],
			 const double pt[], 
			 const TacsScalar scale,
			 const TacsScalar strainSens[], 
			 const TacsScalar Xpts[],
			 const TacsScalar vars[] ){}
  
  // This function adds the sensitivity of the strain to the state variables
  // -----------------------------------------------------------------------
  void addStrainSVSens( TacsScalar strainSVSens[], 
			const double pt[], 
			const TacsScalar scale,
			const TacsScalar strainSens[], 
			const TacsScalar Xpts[],
			const TacsScalar vars[] ){}

  // Functions for retrieving data from the element
  // ----------------------------------------------    
  void setComponentNum( int comp_num ){}
  int getComponentNum(){}
  void addOutputCount( int * nelems, int * nnodes, int * ncsr ) {}
  void getOutputData( unsigned int out_type,
		      double * data, int ld_data,
		      const TacsScalar Xpts[],
		      const TacsScalar vars[] ) {}
  void getOutputConnectivity( int * con, int start_node ) {}
 private:
  // Design variables
  TacsScalar Q; // reduced dynamic pressure
};

#endif // TACS_OSCILLATOR_ELEMENT_H
