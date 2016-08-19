#ifndef TACS_SMD_ELEMENT_H
#define TACS_SMD_ELEMENT_H

#include "TACSElement.h"

/*
  A test implementation of TACSElement class with a Spring mass
  damper system ODE.
*/
class TACSSMDElement : public TACSElement {
 public:
  TACSSMDElement(){}
  ~TACSSMDElement(){}
  
  // Design variables
  static const TacsScalar M = 1.00;
  static const TacsScalar C = 0.02;
  static const TacsScalar K = 5.00;

  // Retrieve information about the name and quantity of variables
  // -------------------------------------------------------------
  const char * elementName() { return  "TACSSpringMassDamperODE"; } 
  const char * displacementName( int i ) {}
  const char * stressName( int i ) {}
  const char * strainName( int i ) {}
  const char * extraName( int i ) {}
  
  // Get the number of displacements, stresses, nodes, etc.
  // ------------------------------------------------------
  int numDisplacements() { return 1; } // should equal num_state_vars
  int numStresses() { return 0; }
  int numNodes() { return 1; }
  int numVariables() { return 1; } // should equal num_state_vars
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

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void getResidual( double time, 
		    TacsScalar res[],
		    const TacsScalar Xpts[],
		    const TacsScalar q[],
		    const TacsScalar qdot[],
		    const TacsScalar qddot[] ) {
    exit(-1);
  }

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void getJacobian( double time, 
		    TacsScalar J[],
		    double alpha, double beta, double gamma,
		    const TacsScalar Xpts[],
		    const TacsScalar q[],
		    const TacsScalar qdot[],
		    const TacsScalar qddot[] ) {
    exit(-1);
  }

  // Compute the Residual of the governing equations
  // -----------------------------------------------
  void addResidual( double time, TacsScalar res[],
		    const TacsScalar Xpts[],
		    const TacsScalar q[],
		    const TacsScalar qdot[],
		    const TacsScalar qddot[] ) {
    res[0] += M*qddot[0] + C*qdot[0] + K*q[0];
  }

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian( double time, TacsScalar J[],
		    double alpha, double beta, double gamma,
		    const TacsScalar Xpts[],
		    const TacsScalar q[],
		    const TacsScalar qdot[],
		    const TacsScalar qddot[] ) {
    J[0] += M*gamma + C*beta + K*alpha;
  }

  // Retrieve the initial values of the state variables
  // --------------------------------------------------
  void getInitCondition( TacsScalar vars[],
			 TacsScalar dvars[],
			 const TacsScalar X[] ){
    vars[0]  = 0.5;
    dvars[0] = 0.1;
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
};

#endif // TACS_SMD_ELEMENT_H
