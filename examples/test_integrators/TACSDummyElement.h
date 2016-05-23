#ifndef TACS_DUMMY_ELEMENT_H
#define TACS_DUMMY_ELEMENT_H

#include "TACSElement.h"

/*
  A dummy implementation of TACSElement class used to test integration
  schemes (BDF and DIRK)
*/
class TACSDummyElement : public TACSElement {
 public:
  TACSDummyElement(){}
  ~TACSDummyElement(){}
  
  // Retrieve information about the name and quantity of variables
  // -------------------------------------------------------------
  const char * elementName() { return  "TACSDummy"; } 
  const char * displacementName( int i ) {}
  const char * stressName( int i ) {}
  const char * strainName( int i ) {}
  const char * extraName( int i ) {}
  
  // Get the number of displacements, stresses, nodes, etc.
  // ------------------------------------------------------
  int numDisplacements() { return 2; }
  int numStresses() { return 0; }
  int numNodes() { return 1; }
  int numVariables() { return 2; }
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

    res[0] = qddot[0] + 0.02*qdot[0]*qdot[1] + 5.0*q[0];
    res[1] = qddot[1] - 0.05*qdot[0]*qdot[1] + q[0]*q[1];

  }

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian( double time, TacsScalar J[],
		    double alpha, double beta, double gamma,
		    const TacsScalar Xpts[],
		    const TacsScalar q[],
		    const TacsScalar qdot[],
		    const TacsScalar qddot[] ) {

    // derivative wrt qddot
    J[0] = gamma*1.0;
    J[1] = gamma*0.0;
    J[2] = gamma*0.0;
    J[3] = gamma*1.0;

    // derivative wrt qdot
    J[0] += beta*0.02*qdot[1];
    J[1] += beta*0.02*qdot[0];
    J[2] += beta*-0.05*qdot[1];
    J[3] += beta*-0.05*qdot[0];

    // derivative wrt q
    J[0] += alpha*5.0;
    J[1] += alpha*0.0;
    J[2] += alpha*q[1];
    J[3] += alpha*q[0];
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
  TACSConstitutive * getConstitutive(){}

  // Get the number of Gauss quadrature points
  // -----------------------------------------
  int getNumGaussPts(){}
  
  // Get the quadrature points and weights
  // -------------------------------------
  double getGaussWtsPts( const int num, double * pt ) {}

  // Get the shape functions from the element
  // ----------------------------------------
  void getShapeFunctions( const double pt[], double N[] ) {}
  
  // Return the determinant of the Jacobian at this point
  // ----------------------------------------------------
  TacsScalar getJacobian( const double * pt, 
			  const TacsScalar Xpts[] ){}
  
  // Return the determinant of the Jacobian and its sensitivity at this point
  // ------------------------------------------------------------------------
  TacsScalar getJacobianXptSens( TacsScalar * hXptSens, 
				 const double * pt, 
				 const TacsScalar Xpts[] ){}

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

#endif // TACS_DUMMY_ELEMENT_H
