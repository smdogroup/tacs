#ifndef TACS_ELEMENT_H
#define TACS_ELEMENT_H

/*
  Basic Element and ElementTraction definitions

  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.

  The purpose of this file is to provide an interface for creating and storing 
  different instances of the finite elements that will be used by TACS. 
  This is what should be extended when including more elements and not the 
  underlying TACS implementation itself. 
*/

#include "TACSObject.h"
#include "TACSConstitutive.h"

enum ElementType { POINT_ELEMENT,
		   EULER_BEAM,
		   TIMOSHENKO_BEAM,
		   PLANE_STRESS,
		   SHELL, 
		   SOLID,
		   Q3D_ELEMENT };

// The different element matrix types
enum ElementMatrixType { STIFFNESS_MATRIX, 
			 MASS_MATRIX, 
			 GEOMETRIC_STIFFNESS_MATRIX, 
			 STIFFNESS_PRODUCT_DERIVATIVE };

//! Provide the matrix in either the normal or the transpose 
enum MatrixOrientation { NORMAL, TRANSPOSE };


/*
  This is the base Element class that all other elements inherit.
  The functions in this class are broken into a few different groups:

  Functions inherited from TACSOptObject:
  ---------------------------------------

  int ownsDesignVar( const int dvNum ) const;
  int getNumDesignVars() const;
  int getDesignVarNums( int * dvNums, int *dvIndex, int dvLen ) const;
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs ) const;
  void getDesignVarRange( TacsScalar lowerBound[], 
			  TacsScalar upperBound[], int numDVs ) const;
  int getDesignVarIndex( int dvNum );

  These functions are used to set and retrieve variable information 
  and set up internal data for sensitivity calculations. Further details
  can be found in TACSObject.h

  Information about sizes and names of variables:
  -----------------------------------------------

  const char * elementName() const;
  const char * displacementName( int i ) const;
  const char * stressName( int i ) const;
  const char * strainName( int i ) const;
  const char * extraName( int i ) const;
  int numDisplacements() const;
  int numStresses() const;
  int numNodes() const;
  int numVariables() const;
  int numExtras() const;
  ElementType getElementType();

  The functions returning char* provide a name for the
  displacement/strain etc.  components that is useful for writing
  output files with named fields. numDisplacements() etc. provides the
  number of components in each field. getElementType() provides an
  enumerated type that idicates the underlying type of element (BEAM,
  SHELL etc.). Note that different element classes can return the same
  ElementType. 

  Functions for analysis:
  -----------------------

  getRes(): Returns the residual of the strain energy components for
  the element. For a linear element this is Ke * ue (where Ke and ue
  are the element stiffness matrix and ue are the element
  variables). Note that surface traction forces and nodal loads are
  applied elsewhere.

  getMat(): Returns the element stiffness matrix and residual (the
  same as above) for given values of the local variables. For a linear
  element this returns Ke, while for a nonlinear element this returns
  the tangent stiffness matrix.

  getMatType(): Return an element matrix of a given type as specified
  by the ElementMatrixType enumeration. This can be used to evaluate
  mass and geometric stiffness matrices. Note that not all element
  classes implement all the matrix types.

  Functions for sensitivity analysis:
  -----------------------------------
  
  There are two types of sensitivities that are considered: one for
  material design variables and the other for geometric design
  variables.  One important assumption is that the geometry variables
  only change the nodal coordinates. Thus, all geometric sensivity
  information can be expanded through a derivative of the element
  residuals w.r.t. the nodal coordinates.

  There are two types of ways sensitivity information is returned:
  For the residual calculations, the derivative of the residual w.r.t.
  the design variable or all nodal coordinates are returned. 
  For the matrix sensitivity routines for the nodal sensitivities,
  the directional derivative is returned: ie, the derivative times the
  components of the sensitivity. (This avoids trying to return a 3rd
  order tensor.)

  getResDVSens(): Return the derivative of the residual w.r.t. the
  currently set design variable. These are usually just the material
  design variables.

  getResXptSens(): Return the (transpose of the) derivative of the
  residual with respect to the nodal coordinates. (Note that the
  transpose is returned because it is easier to compute it in this
  manner. Thus res[i] i = 0 ... nvars-1 represents the derivative
  w.r.t. the first nodal coordinate Xpt[0])

  getMatTypeDVSens(): Return the derivative of a given matrix type
  w.r.t. the set design variable.

  getMatTypeXptSens(): Return the derivative of the matrix w.r.t. the
  nodal coordinates times the sensitivity of the nodal coordinates.
  This is a projected derivative. Note that these computations are
  often much more computationally expensive and can take a
  considerable amount of time.

  Post-processing functions:
  --------------------------

  There are several post-processing calculations that are useful for
  evaluating structural functions of interest. These may be used to
  obtain the underlying consitutive class, the Gauss points and
  weights, the determinant of the Jacobian and the the strain and
  derivative of the strain. These are used by the TACSFunction classes
  for evaluating the compliance and the K.S. function.

  TACSConstitutive * getConstitutive(): return the constitutive relationship

  getJacobian(): Get the determinant of the Jacobian evaluated at the 
  specified parametric location

  getJacobianXptSens(): Return the derivative of the determinant of 
  the Jacobian w.r.t. all nodal coordinates.

  getPtwiseStrain(): Retrieve the strain evaluated at a parametric location 

  getPtwiseStrainXptSens(): The derivative of the strain w.r.t. all nodal
  locations

  addPtwiseStrainSVSens(): The derivative of the strain w.r.t. all nodal
  displacements  

  Functions for aerostructural coupling
  -------------------------------------

  The following functions are used to transmit a set of forces to the
  structure through a series of rigid links. The absolution location
  of the aerodynamic points are provided (Xa). These functions then
  evaluate the rigid displacements (RD) and forces (RF) by evaluating
  the structural surface at the provided parametric location. 
  
  getClosestPt(): Get the parametric locatio of the closest point on
  the element to the physical location provided. This is required for
  pre-processing the closest parametric locations to all the
  aerodynamic points in advance.

  The displacement extrapolation functions:
  getRD(): Extrapolate the displacement to the aerodynamic surface

  getRDXptSens(): Take the sensitivity of the extranpolation to the
  aerodynamic surface w.r.t. the nodal locations - this provides a
  projected derivative.

  getRDTranspose(): The transpose operation - this is required for
  computing the aerostructural adjoint

  The consistent force calculation:
  getRF(): Get the force for the element given the aerodynamic force
  and the parametric location 
  
  getRFXptSens(): Get the derivative of the force w.r.t. the nodal 
  locations - a projected derivative.

  getRFTranspose(): Get the transpose of the force transfer

  Functions for generating output
  -------------------------------
  
  Visualization is a critical component of any finite-element analysis.
  The following functions are used to retrieve the element-level data
  for post-analysis visualization. 

  addOutputCount(): Add the number of nodes and number of csr entries
  that are requried for this element

  getOutputData(): Get the output data of a given type from the element.

  getOutputConnectivity(): Get the connectivity of the data
*/
class TACSElement : public TACSOptObject {
 public: 
  // These constants are used to denote which output to obtain
  // ---------------------------------------------------------
  static const unsigned int OUTPUT_NODES = 1;
  static const unsigned int OUTPUT_DISPLACEMENTS = 2;
  static const unsigned int OUTPUT_STRAINS = 4;
  static const unsigned int OUTPUT_STRESSES = 8;
  static const unsigned int OUTPUT_EXTRAS = 16;
  static const unsigned int OUTPUT_COORDINATES = 32;
  static const unsigned int OUTPUT_DESIGN_VARIABLES = 64;

  TACSElement( int _componentNum = 0 ){
    componentNum = _componentNum;
  }
  virtual ~TACSElement(){}

  // Retrieve information about the name and quantity of variables
  // -------------------------------------------------------------
  virtual const char * elementName() const = 0;
  virtual const char * displacementName( int i ) const = 0;
  virtual const char * stressName( int i ) const = 0;
  virtual const char * strainName( int i ) const = 0;
  virtual const char * extraName( int i ) const = 0;

  // Get the number of displacements, stresses, nodes, etc.
  // ------------------------------------------------------
  virtual int numDisplacements() const = 0;
  virtual int numStresses() const = 0;
  virtual int numNodes() const = 0;
  virtual int numVariables() const = 0;
  virtual int numExtras() const = 0;
  virtual enum ElementType getElementType() = 0;

  // Return the name of the element
  // ------------------------------
  virtual const char * TACSObjectName(){ return this->elementName(); }

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  virtual void computeEnergies( TacsScalar *_Te, 
				TacsScalar *_Pe,
				const TacsScalar Xpts[],
				const TacsScalar vars[],
				const TacsScalar dvars[] ) = 0;

  // Compute the residual of the governing equations
  // -----------------------------------------------
  virtual void getResidual( TacsScalar res[],
			    const TacsScalar Xpts[],
			    const TacsScalar vars[],
			    const TacsScalar dvars[],
			    const TacsScalar ddvars[] ) = 0;

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  virtual void getJacobian( TacsScalar J[],
			    double alpha, double beta, double gamma,
			    const TacsScalar Xpts[],
			    const TacsScalar vars[],
			    const TacsScalar dvars[],
			    const TacsScalar ddvars[] ) = 0;

  // Add the product of the adjoint with the derivative of the design variables
  // --------------------------------------------------------------------------
  virtual void addAdjResProduct( double scale,
				 TacsScalar dvSens[], int dvLen,
				 const TacsScalar psi[],
				 const TacsScalar Xpts[],
				 const TacsScalar vars[],
				 const TacsScalar dvars[],
				 const TacsScalar ddvars[] ) = 0;

  // Add the product of the adjoint with the derivative of the design variables
  // --------------------------------------------------------------------------  
  virtual void getAdjResXptProduct( TacsScalar XptSens[],
				    const TacsScalar psi[],
				    const TacsScalar Xpts[],
				    const TacsScalar vars[],
				    const TacsScalar dvars[],
				    const TacsScalar ddvars[] ) = 0;

  // Retrieve a specific time-independent matrix from the element
  // ------------------------------------------------------------
  virtual void getMatType( ElementMatrixType matType, 
			   TacsScalar mat[], 
			   const TacsScalar Xpts[],
			   const TacsScalar vars[] ){}

  // Compute the derivative of the inner product w.r.t. design variables
  // -------------------------------------------------------------------
  virtual void addMatDVSensInnerProduct( ElementMatrixType matType, 
					 double scale,
					 TacsScalar dvSens[], int dvLen,
					 const TacsScalar psi[], 
					 const TacsScalar phi[],
					 const TacsScalar Xpts[],
					 const TacsScalar vars[] ){}

  // Compute the derivative of the inner product w.r.t. vars[]
  // ---------------------------------------------------------
  virtual void getMatSVSensInnerProduct( ElementMatrixType matType, 
					 TacsScalar res[],
					 const TacsScalar psi[], 
					 const TacsScalar phi[],
					 const TacsScalar Xpts[],
					 const TacsScalar vars[] ){}

  // Member functions for evaluating global functions of interest
  // ------------------------------------------------------------
  virtual TACSConstitutive * getConstitutive() = 0;

  // Get the number of Gauss quadrature points
  // -----------------------------------------
  virtual int getNumGaussPts() = 0;

  // Get the quadrature points and weights
  // -------------------------------------
  virtual double getGaussWtsPts( const int num, double * pt ) = 0;

  // Get the shape functions from the element
  // ----------------------------------------
  virtual void getShapeFunctions( const double pt[], double N[] ){}
  
  // Return the determinant of the Jacobian at this point
  // ----------------------------------------------------
  virtual TacsScalar getJacobian( const double * pt, 
				  const TacsScalar Xpts[] ) = 0;

  // Return the determinant of the Jacobian and its sensitivity at this point
  // ------------------------------------------------------------------------
  virtual TacsScalar getJacobianXptSens( TacsScalar * hXptSens, 
                                         const double * pt, 
					 const TacsScalar Xpts[] );

  // This function returns the strain evaluated at pt
  // ------------------------------------------------
  virtual void getStrain( TacsScalar strain[], 
			  const double pt[], 
			  const TacsScalar Xpts[],
			  const TacsScalar vars[] ) = 0;

  // This function returns the sensitivity of the strain w.r.t. Xpts
  // ---------------------------------------------------------------
  virtual void addStrainXptSens( TacsScalar strainXptSens[],
				 const double pt[], 
				 const TacsScalar scale,
				 const TacsScalar strainSens[], 
				 const TacsScalar Xpts[],
				 const TacsScalar vars[] );
  
  // This function adds the sensitivity of the strain to the state variables
  // -----------------------------------------------------------------------
  virtual void addStrainSVSens( TacsScalar strainSVSens[], 
				const double pt[], 
				const TacsScalar scale,
				const TacsScalar strainSens[], 
				const TacsScalar Xpts[],
				const TacsScalar vars[] );

  // Functions for retrieving data from the element
  // ----------------------------------------------    
  void setComponentNum( int comp_num ){ componentNum = comp_num; }
  int getComponentNum(){ return componentNum; }
  virtual void addOutputCount( int * nelems, int * nnodes, int * ncsr ) = 0;
  virtual void getOutputData( unsigned int out_type,
			      double * data, int ld_data,
			      const TacsScalar Xpts[],
			      const TacsScalar vars[] );
  virtual void getOutputConnectivity( int * con, int start_node ) = 0;

 private: 
  int componentNum;
};

/*
  This class is used to test the element implementation

  Each function tests a single function within the Element class for
  internal self-consistency. These functions test three types of 
  quantities:

  1. Consistency between the residual calculation and the calculation
  of the Jacobian of the residuals. 

  2. Consistency between the element calculations and derivatives
  with respect to the nodal coordinates (XptSens functions.)

  3. Consistency between element calculations and the derivatives
  with respect to material design variables (DVSens functions.)

  The error in each component is measured against the global relative
  tolerance fail_rtol and the global absolute tolerance fail_atol. If
  the component in any component of the derivatives is greater than
  either of these tolerances, the function returns a non-zero fail
  flag (indicating failure.) Otherwise the function returns 0 to
  indicate that no failure has occured. It is important to note that
  the absolute tolerance check is sensitive to the order of magnitude
  of the quantities. Setting the print level to greater than 0 prints
  all components.

  Note that if no variables are supplied to the class, it generates a
  random set of variables on the interval [-1, 1]. (Such large
  displacements often produce large residuals.)
*/
class TestElement : public TACSObject {
 public:
  TestElement( TACSElement * _element, 
	       const TacsScalar _Xpts[] );
  ~TestElement();

  // Set parameters within the test object
  // -------------------------------------
  void setFailTolerances( double _fail_rtol, double _fail_atol ){
    fail_rtol = _fail_rtol;
    fail_atol = _fail_atol;
  }

  void setPrintLevel( int _flag ){
    print_level = _flag;
  }

  void setStepSize( TacsScalar _dh ){
    dh = _dh;
  }

  // Tests for consistency amongst the analysis functions
  // ----------------------------------------------------
  int testStiffnessMat( int col = -1 );
  int testMatDVSens( ElementMatrixType type );
  int testStrainSVSens( const double pt[] );
  
  // Tests for the sensitivities w.r.t. nodal coordinates
  // ----------------------------------------------------
  int testJacobianXptSens( const double pt[] );
  int testStrainXptSens( const double pt[] );
  int testResXptSens( TacsScalar alpha = 1.0 );

  // Design variable sensitivity tests
  // ---------------------------------
  int testResDVSens();
  // int testForceTransferDVSens();
  // int testStiffnessMatDVSens();
  // int testMassMatDVSens();

 private:
  TacsScalar dh; // Step size

  // print_level: 0 - print nothing, 1 - print summary, 2 - print everything
  int print_level; 

  // A test fails if the relative or absolute tolerances are greater than
  // these values
  double fail_rtol, fail_atol;

  TACSElement * element;
  TacsScalar *vars, *Xpts;
};

/*
  Test the implementation of the constitutive class
*/

class TestConstitutive : public TACSObject {
 public:
  TestConstitutive( TACSConstitutive * _con );
  ~TestConstitutive();

  // Set parameters within the test object
  // -------------------------------------
  void setFailTolerances( double _fail_rtol, double _fail_atol ){
    fail_rtol = _fail_rtol;
    fail_atol = _fail_atol;
  }

  void setPrintLevel( int _flag ){
    print_level = _flag;
  }

  void setStepSize( TacsScalar _dh ){
    dh = _dh;
  }
  
  // Test the failure load implementation
  // ------------------------------------
  int testFailStrainSens( const double pt[] );
  int testFailDVSens( const double pt[] );

  // Test the buckling implementation
  // --------------------------------
  int testBucklingStrainSens();
  int testBucklingDVSens();

  // Test the mass implementation
  // ----------------------------
  int testMassDVSens( const double pt[] );
  
 private:
  void compute_strain( TacsScalar strain[], 
                       const double pt[],
                       const TacsScalar stress[] );

  // The constitutive relationship to test
  TACSConstitutive * con;

  TacsScalar dh; // Step size

  // print_level: 0 - print nothing, 1 - print summary, 2 - print everything
  int print_level; 

  // A test fails if the relative or absolute tolerances are greater than
  // these values
  double fail_rtol, fail_atol;  
};

#endif
