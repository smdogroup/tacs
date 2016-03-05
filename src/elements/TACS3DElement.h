#ifndef TACS_3D_ELEMENT_H
#define TACS_3D_ELEMENT_H

/*
  The following file contains the general definition of a
  three-dimensional element that can be used in TACS.

  Copyright (c) 2014 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include "TACSElement.h"
#include "SolidStiffness.h"

/*
  The following class defines a generic three-dimensional element
  without defining the strain expressions, shape functions or
  quadrature scheme. This class can be used to implement multiple 3D
  elements that could be either linear/nonlinear later on.  This does
  not significantly impact the computational performance since the
  cost of the element computations is consumed in the inner product of
  the B-matrix with the constitutive matrix.
*/
template<int NUM_NODES>
class TACS3DElement : public TACSElement {
 public:
  // Define some constants for this element type
  static const int NUM_DISPS = 3;
  static const int NUM_STRESSES = 6;
  static const int NUM_EXTRAS = 4;
  static const int NUM_VARIABLES = 3*NUM_NODES;

  TACS3DElement( SolidStiffness * _stiff,
		 int _linear_strain, int _component );
  ~TACS3DElement();

  // Get the Gauss point and weight information
  // ------------------------------------------
  virtual int getNumQuadraturePoints() = 0;
  virtual double getQuadraturePoint( int n, double pt[] ) = 0;

  // Retrieve the shape functions
  // ----------------------------
  virtual void getShapeFunctions( const double pt[], double N[],
				  double Na[], double Nb[], double Nc[] ) = 0;

  // Compute the position vector and its derivative
  // ----------------------------------------------
  void solidJacobian( TacsScalar X[], TacsScalar Xa[],
		      const double N[], const double Na[], const double Nb[], 
		      const double Nc[], const TacsScalar Xpts[] );
  
  // Compute the displacement gradient
  // --------------------------------
  void getDeformGradient( TacsScalar Ud[], const TacsScalar J[],
			  const double Na[], const double Nb[], const double Nc[],
			  const TacsScalar vars[] );
  void getDeformGradientSens( TacsScalar Ud[], TacsScalar UdSens[],
			      const TacsScalar J[], const TacsScalar JSens[],
			      const double Na[], const double Nb[], const double Nc[],
			      const TacsScalar vars[] );

  // Compute the strain in the element
  // ---------------------------------
  void getStrain( TacsScalar strain[], const TacsScalar J[],
		  const double Na[], const double Nb[], const double Nc[],
		  const TacsScalar vars[] );
 
  // Compute the derivative of the strain with respect to vars
  // ---------------------------------------------------------
  void getBmat( TacsScalar B[], const TacsScalar J[],
		const double Na[], const double Nb[], const double Nc[],
		const TacsScalar vars[] ); 

  // Add the second derivatives of the strain times the stress to the
  // upper portion of the matrix
  // ----------------------------------------------------------------
  void addGeoStiffness( TacsScalar kmat[], TacsScalar h, const TacsScalar stress[],
			const TacsScalar J[], const double Na[], 
			const double Nb[], const double Nc[] );
  
  // Compute the derivative of the strain with respect to the nodal coordinates
  // --------------------------------------------------------------------------
  void getStrainXptSens( TacsScalar sens[], const TacsScalar J[], const TacsScalar Xa[], 
			 const double Na[],  const double Nb[], const double Nc[],
			 const TacsScalar vars[] );

  // The design variable query functions
  // -----------------------------------
  int ownsDesignVar( const int dvNum ) const;
  int getNumDesignVars() const;
  int getDesignVarNums( int * dvNums, int *dvIndex, int dvLen ) const;
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs ) const;
  void getDesignVarRange( TacsScalar lowerBound[], 
			  TacsScalar upperBound[], int numDVs ) const;
  
  // Get the variable information
  // ----------------------------
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
  
  // Functions for computing the residual/stiffness matrix
  // -----------------------------------------------------
  void getRes( TacsScalar * res, const TacsScalar vars[], 
	       const TacsScalar Xpts[] );
  void getMat( TacsScalar * mat, TacsScalar * res, 
	       const TacsScalar vars[], 
	       const TacsScalar Xpts[], MatrixOrientation matOr );
  void getMatType( ElementMatrixTypes matType, TacsScalar scaleFactor, 
		   TacsScalar * mat, const TacsScalar vars[], const TacsScalar Xpts[], 
		   MatrixOrientation matOr );

  // Functions required to determine the derivatives w.r.t. the design variables
  // ---------------------------------------------------------------------------
  void getResXptSens( TacsScalar * res, const TacsScalar vars[], 
		      const TacsScalar Xpts[] );
  void getResDVSens( int dvNum, TacsScalar * res, const TacsScalar vars[], 
		     const TacsScalar Xpts[] );  
  void addAdjResDVSensProduct( TacsScalar alpha, TacsScalar dvSens[], 
			       int dvLen, const TacsScalar psi[],
			       const TacsScalar vars[], const TacsScalar Xpts[] );
  void addMatDVSensInnerProduct( ElementMatrixTypes matType, TacsScalar alpha,
				 TacsScalar dvSens[], int dvLen,
				 const TacsScalar psi[], const TacsScalar phi[],
				 const TacsScalar vars[], const TacsScalar Xpts[] );
  void getMatTypeDVSens( int dvNum, 
                         ElementMatrixTypes matType, TacsScalar scaleFactor, 
			 TacsScalar * mat, const TacsScalar vars[], 
			 const TacsScalar Xpts[], MatrixOrientation matOr );

  // Functions for evaluating global functionals of interest
  // -------------------------------------------------------
  TACSConstitutive * getConstitutive(){ return stiff; }

  // Evaluate the determinant of the Jacobian and its derivative
  // -----------------------------------------------------------
  TacsScalar getJacobian( const double * pt, const TacsScalar Xpts[] );
  TacsScalar getJacobianXptSens( TacsScalar * sh, const double * pt, 
				 const TacsScalar Xpts[] );

  // Compute the point-wise strain and its derivative
  // ------------------------------------------------
  void getPtwiseStrain( TacsScalar strain[], const double pt[],
			const TacsScalar vars[], const TacsScalar Xpts[] );
  void getPtwiseStrainXptSens( TacsScalar strain[], TacsScalar sens[],
			       const double pt[], const TacsScalar vars[], 
			       const TacsScalar Xpts[] );
  void addPtwiseStrainSVSens( TacsScalar sens[], const double pt[], 
			      const TacsScalar scaleFactor, 
			      const TacsScalar strainSens[],
			      const TacsScalar vars[], const TacsScalar Xpts[] );
  
  // Set how the aerodynamic forces are computed
  // -------------------------------------------
  static int USE_RIGID_MOMENT;
  static void setRigidMomentFlag( int _flag );

  // Functions used mostly for aerostructural coupling
  // -------------------------------------------------
  void getDisplacement( TacsScalar U[], const double pt[],
			const TacsScalar vars[] );
  void getRD( TacsScalar U[], 
	      const TacsScalar Xa[], const double pt[], 
	      const TacsScalar vars[], const TacsScalar Xpts[] );
  void getRDXptSens( TacsScalar U[], 
		     const TacsScalar Xa[], const TacsScalar XaSens[], 
		     const double pt[], const TacsScalar vars[], 
		     const TacsScalar Xpts[], const TacsScalar XptSens[] );
  void getRDTranspose( TacsScalar elemAdj[], const TacsScalar Uaero[], 
		       const TacsScalar Xa[], const double pt[], 
		       const TacsScalar Xpts[] );
  void getRF( TacsScalar res[], const TacsScalar Fnode[], 
	      const TacsScalar Xa[], const double pt[], 
	      const TacsScalar Xpts[] );
  void getRFXptSens( TacsScalar res[], const TacsScalar F[], 
		     const TacsScalar Xa[], const TacsScalar XaSens[], 
		     const double pt[], const TacsScalar Xpts[], 
		     const TacsScalar XptSens[] );
  void getRFTranspose( TacsScalar Fnode[], const TacsScalar Xa[], 
		       const double pt[], const TacsScalar res[], 
		       const TacsScalar Xpts[] );

 protected:
  int linear_strain;
  SolidStiffness * stiff;

 private:
  static const char * dispNames[NUM_DISPS];
  static const char * stressNames[NUM_STRESSES];
  static const char * strainNames[NUM_STRESSES];
  static const char * extraNames[NUM_EXTRAS];
};

/*
  The constructor for the 3D element matrix
*/
template <int NUM_NODES>
TACS3DElement<NUM_NODES>::TACS3DElement( SolidStiffness * _stiff,
					 int _linear_strain, 
					 int component ) :
TACSElement(component){
  linear_strain = _linear_strain;
  stiff = _stiff;
  stiff->incref();
}
  
template <int NUM_NODES>
TACS3DElement<NUM_NODES>::~TACS3DElement(){
  stiff->decref();
}

/*
  Provide the names for the different components of the displacements
  and stresses
*/ 
template <int NUM_NODES>
const char * TACS3DElement<NUM_NODES>::dispNames[] = { "u", "v", "w" };

template <int NUM_NODES>
const char * TACS3DElement<NUM_NODES>::stressNames[] = { "sxx", "syy", "szz", 
							 "syz", "sxz", "sxy" };
 
template <int NUM_NODES>
const char * TACS3DElement<NUM_NODES>::strainNames[] = { "exx", "eyy", "ezz", 
							 "eyz", "exz", "exy" };

template <int NUM_NODES>
const char * TACS3DElement<NUM_NODES>::extraNames[] = { "lambda", "buckling",
							"dv1", "dv2" };

/*
  Get the names of the displacements/stress etc.
*/
template <int NUM_NODES>
const char * TACS3DElement<NUM_NODES>::displacementName( int i ) const {
  if (i >= 0 && i < NUM_DISPS){
    return dispNames[i];
  }
  return NULL;
}
  
template <int NUM_NODES>
const char * TACS3DElement<NUM_NODES>::stressName( int i ) const { 
  if (i >= 0 && i < NUM_STRESSES){
    return stressNames[i];
  }
  return NULL;
}

template <int NUM_NODES>
const char * TACS3DElement<NUM_NODES>::strainName( int i ) const {
  if (i >= 0 && i < NUM_STRESSES){
    return strainNames[i];
  }
  return NULL;
}

template <int NUM_NODES>
const char * TACS3DElement<NUM_NODES>::extraName( int i ) const {
  if (i >= 0 && i < NUM_EXTRAS){
    return extraNames[i];
  }
  return NULL;
}

/*
  Retrieve information about the number of displacements/stress etc.
*/
template <int NUM_NODES>
int TACS3DElement<NUM_NODES>::numDisplacements() const { 
  return NUM_DISPS; 
}

template <int NUM_NODES>
int TACS3DElement<NUM_NODES>::numStresses() const { 
  return NUM_STRESSES; 
}

template <int NUM_NODES>
int TACS3DElement<NUM_NODES>::numNodes() const { 
  return NUM_NODES; 
}

template <int NUM_NODES>
int TACS3DElement<NUM_NODES>::numVariables() const { 
  return NUM_VARIABLES; 
}

template <int NUM_NODES>
int TACS3DElement<NUM_NODES>::numExtras() const { 
  return NUM_EXTRAS; 
}

template <int NUM_NODES>
ElementType TACS3DElement<NUM_NODES>::getElementType(){ 
  return SOLID; 
}

/*
  Is the element modified by this element?
*/
template <int NUM_NODES>
int TACS3DElement<NUM_NODES>::ownsDesignVar( const int dvNum ) const { 
  return stiff->ownsDesignVar(dvNum);
}

/*
  Return the design variable numbers that modify this element
*/
template <int NUM_NODES>
int TACS3DElement<NUM_NODES>::getDesignVarNums( int * dvNums, 
						int * dvIndex,
						int dvLen ) const {
  return stiff->getDesignVarNums(dvNums, dvIndex, dvLen);
}

/*
  Retrun the number of design variables modifying this element
*/
template <int NUM_NODES>
int TACS3DElement<NUM_NODES>::getNumDesignVars() const {
  return stiff->getNumDesignVars();
}

/*
  Set the design variable values
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::setDesignVars( const TacsScalar dvs[],
					      int numDVs ){ 
  stiff->setDesignVars(dvs, numDVs);
}

/*
  Retrive the design variable numbers
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getDesignVars( TacsScalar dvs[], 
					      int numDVs ) const {
  stiff->getDesignVars(dvs, numDVs);
}

/*
  Set the design variable lower/upper bounds in the provided arrays
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getDesignVarRange( TacsScalar lowerBound[], 
						  TacsScalar upperBound[], 
						  int numDVs ) const {
  stiff->getDesignVarRange(lowerBound, upperBound, numDVs);
}

/*
  Compute the product of the shape functions with the element shape
  functions.
  
  output:
  X:   the x,y,z location within the structure
  Xa:  the derivative of x,y,z with respect to the parametric location

  input:
  N:           the shape functions
  Na, Nb, Nc:  the derivative of the shape functions
  Xpts:        the nodal locations
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::solidJacobian( TacsScalar X[], 
					      TacsScalar Xa[],
					      const double N[], 
					      const double Na[], 
					      const double Nb[], 
					      const double Nc[],
					      const TacsScalar Xpts[] ){
  X[0] = X[1] = X[2] = 0.0;
  Xa[0] = Xa[1] = Xa[2] = 0.0;
  Xa[3] = Xa[4] = Xa[5] = 0.0;
  Xa[6] = Xa[7] = Xa[8] = 0.0;

  for ( int i = 0; i < NUM_NODES; i++ ){
    X[0] += Xpts[0]*N[0];
    X[1] += Xpts[1]*N[0];
    X[2] += Xpts[2]*N[0];

    Xa[0] += Xpts[0]*Na[0];
    Xa[1] += Xpts[0]*Nb[0];
    Xa[2] += Xpts[0]*Nc[0];
    
    Xa[3] += Xpts[1]*Na[0];
    Xa[4] += Xpts[1]*Nb[0];
    Xa[5] += Xpts[1]*Nc[0];
    
    Xa[6] += Xpts[2]*Na[0];
    Xa[7] += Xpts[2]*Nb[0];
    Xa[8] += Xpts[2]*Nc[0];
    
    N++;
    Na++; Nb++; Nc++;
    Xpts += 3;
  }
}

/*
  Compute the displacement gradient using the provided basis functions.

  The displacement gradient is computed as follows:
  Ud = Ua*J =
  [ Ua[0]  Ua[1]  Ua[2] ][ J[0]  J[1]  J[2] ]
  [ Ua[3]  Ua[4]  Ua[5] ][ J[3]  J[4]  J[5] ]
  [ Ua[6]  Ua[7]  Ua[8] ][ J[6]  J[7]  J[8] ]

  output:
  Ud:    the displacement gradient

  input:
  J:          the transformation J = [X,a]^{-1}
  Na, Nb, Nc: the derivatives of the shape functions
  vars:       the element variables
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getDeformGradient( TacsScalar Ud[],
						  const TacsScalar J[],
						  const double Na[],
						  const double Nb[], 
						  const double Nc[],
						  const TacsScalar vars[] ){
  TacsScalar Ua[9];
  Ua[0] = Ua[1] = Ua[2] = 
    Ua[3] = Ua[4] = Ua[5] = 
    Ua[6] = Ua[7] = Ua[8] = 0.0;

  // Compute the derivative of the u,v,w displacements with
  // respect to the parametric locations within the element
  for ( int i = 0; i < NUM_NODES; i++ ){
    Ua[0] += vars[0]*Na[0];
    Ua[1] += vars[0]*Nb[0];
    Ua[2] += vars[0]*Nc[0];
    
    Ua[3] += vars[1]*Na[0];
    Ua[4] += vars[1]*Nb[0];
    Ua[5] += vars[1]*Nc[0];
    
    Ua[6] += vars[2]*Na[0];
    Ua[7] += vars[2]*Nb[0];
    Ua[8] += vars[2]*Nc[0];
    
    Na++; Nb++; Nc++;
    vars += 3;
  }

  // Compute the displacement gradient: Ud = Ua*J
  Ud[0] = Ua[0]*J[0] + Ua[1]*J[3] + Ua[2]*J[6];
  Ud[3] = Ua[3]*J[0] + Ua[4]*J[3] + Ua[5]*J[6];
  Ud[6] = Ua[6]*J[0] + Ua[7]*J[3] + Ua[8]*J[6];

  Ud[1] = Ua[0]*J[1] + Ua[1]*J[4] + Ua[2]*J[7];
  Ud[4] = Ua[3]*J[1] + Ua[4]*J[4] + Ua[5]*J[7];
  Ud[7] = Ua[6]*J[1] + Ua[7]*J[4] + Ua[8]*J[7];

  Ud[2] = Ua[0]*J[2] + Ua[1]*J[5] + Ua[2]*J[8];
  Ud[5] = Ua[3]*J[2] + Ua[4]*J[5] + Ua[5]*J[8];
  Ud[8] = Ua[6]*J[2] + Ua[7]*J[5] + Ua[8]*J[8];
}

/*
  Compute the displacement gradient and the derivative of the
  displacement gradient using the provided basis functions.

  The displacement gradient is computed as follows:
  Ud = Ua*J =
  [ Ua[0]  Ua[1]  Ua[2] ][ J[0]  J[1]  J[2] ]
  [ Ua[3]  Ua[4]  Ua[5] ][ J[3]  J[4]  J[5] ]
  [ Ua[6]  Ua[7]  Ua[8] ][ J[6]  J[7]  J[8] ]

  output:
  Ud:       the displacement gradient
  UdSens:   the derivative of the displacement gradient

  input:
  J:          the transformation J = [X,a]^{-1}
  JSens:      the derivative of the inverse Jacobian
  Na, Nb, Nc: the derivatives of the shape functions
  vars:       the element variables
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getDeformGradientSens( TacsScalar Ud[], 
						      TacsScalar UdSens[],
						      const TacsScalar J[],
						      const TacsScalar JSens[],
						      const double Na[],
						      const double Nb[], 
						      const double Nc[],
						      const TacsScalar vars[] ){
  TacsScalar Ua[9];
  Ua[0] = Ua[1] = Ua[2] = 
    Ua[3] = Ua[4] = Ua[5] = 
    Ua[6] = Ua[7] = Ua[8] = 0.0;

  // Compute the derivative of the u,v,w displacements with
  // respect to the parametric locations within the element
  for ( int i = 0; i < NUM_NODES; i++ ){
    Ua[0] += vars[0]*Na[0];
    Ua[1] += vars[0]*Nb[0];
    Ua[2] += vars[0]*Nc[0];
    
    Ua[3] += vars[1]*Na[0];
    Ua[4] += vars[1]*Nb[0];
    Ua[5] += vars[1]*Nc[0];
    
    Ua[6] += vars[2]*Na[0];
    Ua[7] += vars[2]*Nb[0];
    Ua[8] += vars[2]*Nc[0];
    
    Na++; Nb++; Nc++;
    vars += 3;
  }

  // Compute the displacement gradient: Ud = Ua*J
  Ud[0] = Ua[0]*J[0] + Ua[1]*J[3] + Ua[2]*J[6];
  Ud[3] = Ua[3]*J[0] + Ua[4]*J[3] + Ua[5]*J[6];
  Ud[6] = Ua[6]*J[0] + Ua[7]*J[3] + Ua[8]*J[6];

  Ud[1] = Ua[0]*J[1] + Ua[1]*J[4] + Ua[2]*J[7];
  Ud[4] = Ua[3]*J[1] + Ua[4]*J[4] + Ua[5]*J[7];
  Ud[7] = Ua[6]*J[1] + Ua[7]*J[4] + Ua[8]*J[7];

  Ud[2] = Ua[0]*J[2] + Ua[1]*J[5] + Ua[2]*J[8];
  Ud[5] = Ua[3]*J[2] + Ua[4]*J[5] + Ua[5]*J[8];
  Ud[8] = Ua[6]*J[2] + Ua[7]*J[5] + Ua[8]*J[8];

  // Compute the derivative of the displacement gradient
  UdSens[0] = Ua[0]*JSens[0] + Ua[1]*JSens[3] + Ua[2]*JSens[6];
  UdSens[3] = Ua[3]*JSens[0] + Ua[4]*JSens[3] + Ua[5]*JSens[6];
  UdSens[6] = Ua[6]*JSens[0] + Ua[7]*JSens[3] + Ua[8]*JSens[6];

  UdSens[1] = Ua[0]*JSens[1] + Ua[1]*JSens[4] + Ua[2]*JSens[7];
  UdSens[4] = Ua[3]*JSens[1] + Ua[4]*JSens[4] + Ua[5]*JSens[7];
  UdSens[7] = Ua[6]*JSens[1] + Ua[7]*JSens[4] + Ua[8]*JSens[7];

  UdSens[2] = Ua[0]*JSens[2] + Ua[1]*JSens[5] + Ua[2]*JSens[8];
  UdSens[5] = Ua[3]*JSens[2] + Ua[4]*JSens[5] + Ua[5]*JSens[8];
  UdSens[8] = Ua[6]*JSens[2] + Ua[7]*JSens[5] + Ua[8]*JSens[8];
}

/*
  Compute the strain using the specified transformation, shape
  functions and variables

  Note: The strain is computed as follows:

  strain = 0.5*(Ud + Ud^{T} + Ud^{T}*Ud)
  =
  [ Ud[0]  0.5*(Ud[1] + Ud[3])  0.5*(Ud[2] + Ud[6]) ]
  [             Ud[4]           0.5*(Ud[5] + Ud[7]) ] + 
  [                                  Ud[8]          ]

  [ Ud[0]  Ud[3]  Ud[6] ][ Ud[0]  Ud[1]  Ud[2] ]
  [ Ud[1]  Ud[4]  Ud[7] ][ Ud[3]  Ud[4]  Ud[5] ]
  [ Ud[2]  Ud[5]  Ud[8] ][ Ud[6]  Ud[7]  Ud[8] ]

  output:
  strain:   the strain
  
  input:
  J:          the Jacobian of the transformation
  Na, Nb, Nc: the derivatives of the basis functions
  vars:       the variables
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getStrain( TacsScalar strain[], 
					  const TacsScalar J[],
					  const double Na[], 
					  const double Nb[], 
					  const double Nc[], 
					  const TacsScalar vars[] ){
  // Compute the displacement gradient
  TacsScalar Ud[9];
  getDeformGradient(Ud, J, Na, Nb, Nc, vars);

  // Compute the strain using either linear or nonlinear expression
  if (linear_strain){
    strain[0] = Ud[0];
    strain[1] = Ud[4];
    strain[2] = Ud[8];

    strain[3] = Ud[5] + Ud[7];
    strain[4] = Ud[2] + Ud[6];
    strain[5] = Ud[1] + Ud[3];
  }
  else {
    strain[0] = Ud[0] + 0.5*(Ud[0]*Ud[0] + Ud[3]*Ud[3] + Ud[6]*Ud[6]);
    strain[1] = Ud[4] + 0.5*(Ud[1]*Ud[1] + Ud[4]*Ud[4] + Ud[7]*Ud[7]);
    strain[2] = Ud[8] + 0.5*(Ud[2]*Ud[2] + Ud[5]*Ud[5] + Ud[8]*Ud[8]);
    
    strain[3] = Ud[5] + Ud[7] + (Ud[1]*Ud[2] + Ud[4]*Ud[5] + Ud[7]*Ud[8]);
    strain[4] = Ud[2] + Ud[6] + (Ud[0]*Ud[2] + Ud[3]*Ud[5] + Ud[6]*Ud[8]);
    strain[5] = Ud[1] + Ud[3] + (Ud[0]*Ud[1] + Ud[3]*Ud[4] + Ud[6]*Ud[7]);
  }
}
  
/*
  Compute the derivative of the strain with respect to vars
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getBmat( TacsScalar B[], 
					const TacsScalar J[],
					const double Na[], 
					const double Nb[], 
					const double Nc[],
					const TacsScalar vars[] ){
  if (linear_strain){
    // If this is a linear element, then things are relative easy to
    // deal with - we just compute B alternatively by row
    for ( int i = 0; i < NUM_NODES; i++ ){
      TacsScalar Dx = Na[0]*J[0] + Nb[0]*J[3] + Nc[0]*J[6];
      TacsScalar Dy = Na[0]*J[1] + Nb[0]*J[4] + Nc[0]*J[7];
      TacsScalar Dz = Na[0]*J[2] + Nb[0]*J[5] + Nc[0]*J[8];

      B[0] = Dx;
      B[1] = 0.0;
      B[2] = 0.0;
      B[3] = 0.0;
      B[4] = Dz;
      B[5] = Dy;
      B += 6;

      B[0] = 0.0;
      B[1] = Dy;
      B[2] = 0.0;
      B[3] = Dz;
      B[4] = 0.0;
      B[5] = Dx;
      B += 6;

      B[0] = 0.0;
      B[1] = 0.0;
      B[2] = Dz;
      B[3] = Dy;
      B[4] = Dx;
      B[5] = 0.0;
      B += 6;

      Na++; Nb++; Nc++;
    }
  }
  else {
    // Compute the displacement gradient: Ud = Ua*J
    TacsScalar Ud[9];
    getDeformGradient(Ud, J, Na, Nb, Nc, vars);

    // Compute the derivative of the strain with respect to
    // the nodal displacements
    for ( int i = 0; i < NUM_NODES; i++ ){
      TacsScalar Dx = Na[0]*J[0] + Nb[0]*J[3] + Nc[0]*J[6];
      TacsScalar Dy = Na[0]*J[1] + Nb[0]*J[4] + Nc[0]*J[7];
      TacsScalar Dz = Na[0]*J[2] + Nb[0]*J[5] + Nc[0]*J[8];

      B[0] = Dx + Ud[0]*Dx;
      B[1] = Ud[1]*Dy;
      B[2] = Ud[2]*Dz;
      B[3] = Dy*Ud[2] + Ud[1]*Dz;
      B[4] = Dz + (Dx*Ud[2] + Ud[0]*Dz);
      B[5] = Dy + (Dx*Ud[1] + Ud[0]*Dy);
      B += 6;
      
      B[0] = Ud[3]*Dx;
      B[1] = Dy + Ud[4]*Dy;
      B[2] = Ud[5]*Dz;
      B[3] = Dz + (Dy*Ud[5] + Ud[4]*Dz);
      B[4] = Dx*Ud[5] + Ud[3]*Dz;
      B[5] = Dx + (Dx*Ud[4] + Ud[3]*Dy);
      B += 6;

      B[0] = Ud[6]*Dx;
      B[1] = Ud[7]*Dy;
      B[2] = Dz + Ud[8]*Dz;
      B[3] = Dy + (Dy*Ud[8] + Ud[7]*Dz);
      B[4] = Dx + (Dx*Ud[8] + Ud[6]*Dz);
      B[5] = Dx*Ud[7] + Ud[6]*Dy;
      B += 6;

      Na++; Nb++; Nc++;
    }
  }
}

/*
  Add the second derivatives of the strain times the stress to the
  upper portion of the matrix 
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::addGeoStiffness( TacsScalar mat[],
						TacsScalar h,
						const TacsScalar stress[],
						const TacsScalar J[], 
						const double Na[], 
						const double Nb[], 
						const double Nc[] ){
  if (!linear_strain){
    for ( int j = 0; j < NUM_NODES; j++ ){
      TacsScalar Dxj = Na[j]*J[0] + Nb[j]*J[3] + Nc[j]*J[6];
      TacsScalar Dyj = Na[j]*J[1] + Nb[j]*J[4] + Nc[j]*J[7];
      TacsScalar Dzj = Na[j]*J[2] + Nb[j]*J[5] + Nc[j]*J[8];
      
      for ( int i = 0; i < NUM_NODES; i++ ){
	TacsScalar Dxi = Na[i]*J[0] + Nb[i]*J[3] + Nc[i]*J[6];
	TacsScalar Dyi = Na[i]*J[1] + Nb[i]*J[4] + Nc[i]*J[7];
	TacsScalar Dzi = Na[i]*J[2] + Nb[i]*J[5] + Nc[i]*J[8];

	// Add the contributions to the stiffness matrix
	mat[3*i + 3*j*NUM_VARIABLES] += h*stress[0]*Dxi*Dxj;
	mat[3*i+1 + 3*j*NUM_VARIABLES] += h*stress[5]*(Dxi*Dyj + Dxj*Dyi);
	mat[3*i+2 + 3*j*NUM_VARIABLES] += h*stress[4]*(Dxi*Dzj + Dxj*Dzi);

	mat[3*i + (3*j+1)*NUM_VARIABLES] += h*stress[4]*(Dxi*Dzj + Dxj*Dzi);
	mat[3*i+1 + (3*j+1)*NUM_VARIABLES] += h*stress[1]*Dyi*Dyj;
	mat[3*i+2 + (3*j+1)*NUM_VARIABLES] += h*stress[3]*(Dyi*Dzj + Dyj*Dzi);

	mat[3*i + (3*j+2)*NUM_VARIABLES] += h*stress[4]*(Dxi*Dzj + Dxj*Dzi);
	mat[3*i+1 + (3*j+2)*NUM_VARIABLES] += h*stress[3]*(Dyi*Dzj + Dyj*Dzi);
	mat[3*i+2 + (3*j+2)*NUM_VARIABLES] += h*stress[2]*Dzi*Dzj;
      }
    }
  }
}
  
/*
  Compute the derivative of the strain with respect to the nodal
  coordinates

  Note that the derivative of the Jacobian matrix can be computed as
  follows:

  dJ/dx = - J*d(Xa)/dx*J

  where d(Xa)/dx is the derivative of the parametric derivative with
  respect to the nodal position. This derivative takes the following
  form:

  dJ/dXa = 
  [ J[0]  J[1]  J[2] ][ Na[0]  Nb[0]  Nc[0] ][ J[0]  J[1]  J[2] ]
  [ J[3]  J[4]  J[5] ][     0      0      0 ][ J[3]  J[4]  J[5] ]
  [ J[6]  J[7]  J[8] ][     0      0      0 ][ J[6]  J[7]  J[8] ]

  (where the middle term is the derivative d(Xa)/dx). These
  derivatives take a special form as follows:

  For the x-derivative:
  dJ/dXa = 
  [ J[0]  J[1]  J[2] ][ d1  d2  d3 ]
  [ J[3]  J[4]  J[5] ][  0   0   0 ]
  [ J[6]  J[7]  J[8] ][  0   0   0 ]

  = 
  [ d1*J[0]  d2*J[0]  d3*J[0] ]
  [ d1*J[3]  d2*J[3]  d3*J[3] ]
  [ d1*J[6]  d2*J[6]  d3*J[6] ]

  where:
  d1 = Na[0]*J[0] + Nb[0]*J[3] + Nc[0]*J[6]
  d2 = Na[0]*J[1] + Nb[0]*J[4] + Nc[0]*J[7]
  d3 = Na[0]*J[2] + Nb[0]*J[5] + Nc[0]*J[8]
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getStrainXptSens( TacsScalar sens[],
						 const TacsScalar J[], 
						 const TacsScalar Xa[],
						 const double Na[], 
						 const double Nb[], 
						 const double Nc[],
						 const TacsScalar vars[] ){
  // Compute the derivative of the u,v,w displacements with
  // respect to the parametric locations within the element
  TacsScalar Ua[9];
  Ua[0] = Ua[1] = Ua[2] = 
    Ua[3] = Ua[4] = Ua[5] = 
    Ua[6] = Ua[7] = Ua[8] = 0.0;

  const double *na = Na, *nb = Nb, *nc = Nc;
  for ( int i = 0; i < NUM_NODES; i++ ){
    Ua[0] += vars[0]*na[0];
    Ua[1] += vars[0]*nb[0];
    Ua[2] += vars[0]*nc[0];
    
    Ua[3] += vars[1]*na[0];
    Ua[4] += vars[1]*nb[0];
    Ua[5] += vars[1]*nc[0];
    
    Ua[6] += vars[2]*na[0];
    Ua[7] += vars[2]*nb[0];
    Ua[8] += vars[2]*nc[0];
    
    na++; nb++; nc++;
    vars += 3;
  }

  if (linear_strain){
    for ( int i = 0; i < NUM_NODES; i++ ){
      // JSens = -J*d(Xa)/dx*J
      TacsScalar d1 = -(Na[0]*J[0] + Nb[0]*J[3] + Nc[0]*J[6]);
      TacsScalar d2 = -(Na[0]*J[1] + Nb[0]*J[4] + Nc[0]*J[7]);
      TacsScalar d3 = -(Na[0]*J[2] + Nb[0]*J[5] + Nc[0]*J[8]);

      // The derivative of the inverse transformation matrix
      // in each coordinate direction
      TacsScalar Ja[9], Jb[9], Jc[9];
      Ja[0] = d1*J[0];   Ja[1] = d2*J[0];   Ja[2] = d3*J[0];
      Ja[3] = d1*J[3];   Ja[4] = d2*J[3];   Ja[5] = d3*J[3];
      Ja[6] = d1*J[6];   Ja[7] = d2*J[6];   Ja[8] = d3*J[6];

      Jb[0] = d1*J[1];   Jb[1] = d2*J[1];   Jb[2] = d3*J[1];
      Jb[3] = d1*J[4];   Jb[4] = d2*J[4];   Jb[5] = d3*J[4];
      Jb[6] = d1*J[7];   Jb[7] = d2*J[7];   Jb[8] = d3*J[7];

      Jc[0] = d1*J[2];   Jc[1] = d2*J[2];   Jc[2] = d3*J[2];
      Jc[3] = d1*J[5];   Jc[4] = d2*J[5];   Jc[5] = d3*J[5];
      Jc[6] = d1*J[8];   Jc[7] = d2*J[8];   Jc[8] = d3*J[8];

      // Compute the derivative of the displacement gradient
      // with respect to each of the three coordinate directions
      // displacement gradient: Ud = Ua*J
      TacsScalar Uda[9];
      Uda[0] = Ua[0]*Ja[0] + Ua[1]*Ja[3] + Ua[2]*Ja[6];
      Uda[3] = Ua[3]*Ja[0] + Ua[4]*Ja[3] + Ua[5]*Ja[6];
      Uda[6] = Ua[6]*Ja[0] + Ua[7]*Ja[3] + Ua[8]*Ja[6];

      Uda[1] = Ua[0]*Ja[1] + Ua[1]*Ja[4] + Ua[2]*Ja[7];
      Uda[4] = Ua[3]*Ja[1] + Ua[4]*Ja[4] + Ua[5]*Ja[7];
      Uda[7] = Ua[6]*Ja[1] + Ua[7]*Ja[4] + Ua[8]*Ja[7];

      Uda[2] = Ua[0]*Ja[2] + Ua[1]*Ja[5] + Ua[2]*Ja[8];
      Uda[5] = Ua[3]*Ja[2] + Ua[4]*Ja[5] + Ua[5]*Ja[8];
      Uda[8] = Ua[6]*Ja[2] + Ua[7]*Ja[5] + Ua[8]*Ja[8];

      TacsScalar Udb[9];
      Udb[0] = Ua[0]*Jb[0] + Ua[1]*Jb[3] + Ua[2]*Jb[6];
      Udb[3] = Ua[3]*Jb[0] + Ua[4]*Jb[3] + Ua[5]*Jb[6];
      Udb[6] = Ua[6]*Jb[0] + Ua[7]*Jb[3] + Ua[8]*Jb[6];

      Udb[1] = Ua[0]*Jb[1] + Ua[1]*Jb[4] + Ua[2]*Jb[7];
      Udb[4] = Ua[3]*Jb[1] + Ua[4]*Jb[4] + Ua[5]*Jb[7];
      Udb[7] = Ua[6]*Jb[1] + Ua[7]*Jb[4] + Ua[8]*Jb[7];

      Udb[2] = Ua[0]*Jb[2] + Ua[1]*Jb[5] + Ua[2]*Jb[8];
      Udb[5] = Ua[3]*Jb[2] + Ua[4]*Jb[5] + Ua[5]*Jb[8];
      Udb[8] = Ua[6]*Jb[2] + Ua[7]*Jb[5] + Ua[8]*Jb[8];

      TacsScalar Udc[9];
      Udc[0] = Ua[0]*Jc[0] + Ua[1]*Jc[3] + Ua[2]*Jc[6];
      Udc[3] = Ua[3]*Jc[0] + Ua[4]*Jc[3] + Ua[5]*Jc[6];
      Udc[6] = Ua[6]*Jc[0] + Ua[7]*Jc[3] + Ua[8]*Jc[6];

      Udc[1] = Ua[0]*Jc[1] + Ua[1]*Jc[4] + Ua[2]*Jc[7];
      Udc[4] = Ua[3]*Jc[1] + Ua[4]*Jc[4] + Ua[5]*Jc[7];
      Udc[7] = Ua[6]*Jc[1] + Ua[7]*Jc[4] + Ua[8]*Jc[7];

      Udc[2] = Ua[0]*Jc[2] + Ua[1]*Jc[5] + Ua[2]*Jc[8];
      Udc[5] = Ua[3]*Jc[2] + Ua[4]*Jc[5] + Ua[5]*Jc[8];
      Udc[8] = Ua[6]*Jc[2] + Ua[7]*Jc[5] + Ua[8]*Jc[8];

      sens[0] = Uda[0];
      sens[1] = Uda[4];
      sens[2] = Uda[8];
      sens[3] = (Uda[5] + Uda[7]);
      sens[4] = (Uda[2] + Uda[6]);
      sens[5] = (Uda[1] + Uda[3]);
      sens += 6;

      sens[0] = Udb[0];
      sens[1] = Udb[4];
      sens[2] = Udb[8];
      sens[3] = (Udb[5] + Udb[7]);
      sens[4] = (Udb[2] + Udb[6]);
      sens[5] = (Udb[1] + Udb[3]);
      sens += 6;

      sens[0] = Udc[0];
      sens[1] = Udc[4];
      sens[2] = Udc[8];
      sens[3] = (Udc[5] + Udc[7]);
      sens[4] = (Udc[2] + Udc[6]);
      sens[5] = (Udc[1] + Udc[3]);
      sens += 6;

      Na++; Nb++; Nc++;
    }
  }
  else {
    // Compute the displacement gradient: Ud = Ua*J
    TacsScalar Ud[9];
    Ud[0] = Ua[0]*J[0] + Ua[1]*J[3] + Ua[2]*J[6];
    Ud[3] = Ua[3]*J[0] + Ua[4]*J[3] + Ua[5]*J[6];
    Ud[6] = Ua[6]*J[0] + Ua[7]*J[3] + Ua[8]*J[6];

    Ud[1] = Ua[0]*J[1] + Ua[1]*J[4] + Ua[2]*J[7];
    Ud[4] = Ua[3]*J[1] + Ua[4]*J[4] + Ua[5]*J[7];
    Ud[7] = Ua[6]*J[1] + Ua[7]*J[4] + Ua[8]*J[7];

    Ud[2] = Ua[0]*J[2] + Ua[1]*J[5] + Ua[2]*J[8];
    Ud[5] = Ua[3]*J[2] + Ua[4]*J[5] + Ua[5]*J[8];
    Ud[8] = Ua[6]*J[2] + Ua[7]*J[5] + Ua[8]*J[8];

    for ( int i = 0; i < NUM_NODES; i++ ){
      // JSens = -J*d(Xa)/dx*J
      TacsScalar d1 = -(Na[0]*J[0] + Nb[0]*J[3] + Nc[0]*J[6]);
      TacsScalar d2 = -(Na[0]*J[1] + Nb[0]*J[4] + Nc[0]*J[7]);
      TacsScalar d3 = -(Na[0]*J[2] + Nb[0]*J[5] + Nc[0]*J[8]);

      // The derivative of the inverse transformation matrix
      // in each coordinate direction
      TacsScalar Ja[9], Jb[9], Jc[9];
      Ja[0] = d1*J[0];   Ja[1] = d2*J[0];   Ja[2] = d3*J[0];
      Ja[3] = d1*J[3];   Ja[4] = d2*J[3];   Ja[5] = d3*J[3];
      Ja[6] = d1*J[6];   Ja[7] = d2*J[6];   Ja[8] = d3*J[6];

      Jb[0] = d1*J[1];   Jb[1] = d2*J[1];   Jb[2] = d3*J[1];
      Jb[3] = d1*J[4];   Jb[4] = d2*J[4];   Jb[5] = d3*J[4];
      Jb[6] = d1*J[7];   Jb[7] = d2*J[7];   Jb[8] = d3*J[7];

      Jc[0] = d1*J[2];   Jc[1] = d2*J[2];   Jc[2] = d3*J[2];
      Jc[3] = d1*J[5];   Jc[4] = d2*J[5];   Jc[5] = d3*J[5];
      Jc[6] = d1*J[8];   Jc[7] = d2*J[8];   Jc[8] = d3*J[8];

      // Compute the derivative of the displacement gradient
      // with respect to each of the three coordinate directions
      // displacement gradient: Ud = Ua*J
      TacsScalar Uda[9];
      Uda[0] = Ua[0]*Ja[0] + Ua[1]*Ja[3] + Ua[2]*Ja[6];
      Uda[3] = Ua[3]*Ja[0] + Ua[4]*Ja[3] + Ua[5]*Ja[6];
      Uda[6] = Ua[6]*Ja[0] + Ua[7]*Ja[3] + Ua[8]*Ja[6];

      Uda[1] = Ua[0]*Ja[1] + Ua[1]*Ja[4] + Ua[2]*Ja[7];
      Uda[4] = Ua[3]*Ja[1] + Ua[4]*Ja[4] + Ua[5]*Ja[7];
      Uda[7] = Ua[6]*Ja[1] + Ua[7]*Ja[4] + Ua[8]*Ja[7];

      Uda[2] = Ua[0]*Ja[2] + Ua[1]*Ja[5] + Ua[2]*Ja[8];
      Uda[5] = Ua[3]*Ja[2] + Ua[4]*Ja[5] + Ua[5]*Ja[8];
      Uda[8] = Ua[6]*Ja[2] + Ua[7]*Ja[5] + Ua[8]*Ja[8];

      TacsScalar Udb[9];
      Udb[0] = Ua[0]*Jb[0] + Ua[1]*Jb[3] + Ua[2]*Jb[6];
      Udb[3] = Ua[3]*Jb[0] + Ua[4]*Jb[3] + Ua[5]*Jb[6];
      Udb[6] = Ua[6]*Jb[0] + Ua[7]*Jb[3] + Ua[8]*Jb[6];

      Udb[1] = Ua[0]*Jb[1] + Ua[1]*Jb[4] + Ua[2]*Jb[7];
      Udb[4] = Ua[3]*Jb[1] + Ua[4]*Jb[4] + Ua[5]*Jb[7];
      Udb[7] = Ua[6]*Jb[1] + Ua[7]*Jb[4] + Ua[8]*Jb[7];

      Udb[2] = Ua[0]*Jb[2] + Ua[1]*Jb[5] + Ua[2]*Jb[8];
      Udb[5] = Ua[3]*Jb[2] + Ua[4]*Jb[5] + Ua[5]*Jb[8];
      Udb[8] = Ua[6]*Jb[2] + Ua[7]*Jb[5] + Ua[8]*Jb[8];

      TacsScalar Udc[9];
      Udc[0] = Ua[0]*Jc[0] + Ua[1]*Jc[3] + Ua[2]*Jc[6];
      Udc[3] = Ua[3]*Jc[0] + Ua[4]*Jc[3] + Ua[5]*Jc[6];
      Udc[6] = Ua[6]*Jc[0] + Ua[7]*Jc[3] + Ua[8]*Jc[6];

      Udc[1] = Ua[0]*Jc[1] + Ua[1]*Jc[4] + Ua[2]*Jc[7];
      Udc[4] = Ua[3]*Jc[1] + Ua[4]*Jc[4] + Ua[5]*Jc[7];
      Udc[7] = Ua[6]*Jc[1] + Ua[7]*Jc[4] + Ua[8]*Jc[7];

      Udc[2] = Ua[0]*Jc[2] + Ua[1]*Jc[5] + Ua[2]*Jc[8];
      Udc[5] = Ua[3]*Jc[2] + Ua[4]*Jc[5] + Ua[5]*Jc[8];
      Udc[8] = Ua[6]*Jc[2] + Ua[7]*Jc[5] + Ua[8]*Jc[8];

      sens[0] = Uda[0] + (Ud[0]*Uda[0] + Ud[3]*Uda[3] + Ud[6]*Uda[6]);
      sens[1] = Uda[4] + (Ud[1]*Uda[1] + Ud[4]*Uda[4] + Ud[7]*Uda[7]);
      sens[2] = Uda[8] + (Ud[2]*Uda[2] + Ud[5]*Uda[5] + Ud[8]*Uda[8]);
      sens[3] = Uda[5] + Uda[7] + (Uda[1]*Ud[2] + Uda[4]*Ud[5] + Uda[7]*Ud[8] +
				   Ud[1]*Uda[2] + Ud[4]*Uda[5] + Ud[7]*Uda[8]);
      sens[4] = Uda[2] + Uda[6] + (Uda[0]*Ud[2] + Uda[3]*Ud[5] + Uda[6]*Ud[8] +
				   Ud[0]*Uda[2] + Ud[3]*Uda[5] + Ud[6]*Uda[8]);
      sens[5] = Uda[1] + Uda[3] + (Uda[0]*Ud[1] + Uda[3]*Ud[4] + Uda[6]*Ud[7] +
				   Ud[0]*Uda[1] + Ud[3]*Uda[4] + Ud[6]*Uda[7]);
      sens += 6;

      sens[0] = Udb[0] + (Ud[0]*Udb[0] + Ud[3]*Udb[3] + Ud[6]*Udb[6]);
      sens[1] = Udb[4] + (Ud[1]*Udb[1] + Ud[4]*Udb[4] + Ud[7]*Udb[7]);
      sens[2] = Udb[8] + (Ud[2]*Udb[2] + Ud[5]*Udb[5] + Ud[8]*Udb[8]);
      sens[3] = Udb[5] + Udb[7] + (Udb[1]*Ud[2] + Udb[4]*Ud[5] + Udb[7]*Ud[8] +
				   Ud[1]*Udb[2] + Ud[4]*Udb[5] + Ud[7]*Udb[8]);
      sens[4] = Udb[2] + Udb[6] + (Udb[0]*Ud[2] + Udb[3]*Ud[5] + Udb[6]*Ud[8] +
				   Ud[0]*Udb[2] + Ud[3]*Udb[5] + Ud[6]*Udb[8]);
      sens[5] = Udb[1] + Udb[3] + (Udb[0]*Ud[1] + Udb[3]*Ud[4] + Udb[6]*Ud[7] +
				   Ud[0]*Udb[1] + Ud[3]*Udb[4] + Ud[6]*Udb[7]);
      sens += 6;

      sens[0] = Udc[0] + (Ud[0]*Udc[0] + Ud[3]*Udc[3] + Ud[6]*Udc[6]);
      sens[1] = Udc[4] + (Ud[1]*Udc[1] + Ud[4]*Udc[4] + Ud[7]*Udc[7]);
      sens[2] = Udc[8] + (Ud[2]*Udc[2] + Ud[5]*Udc[5] + Ud[8]*Udc[8]);
      sens[3] = Udc[5] + Udc[7] + (Udc[1]*Ud[2] + Udc[4]*Ud[5] + Udc[7]*Ud[8] +
				   Ud[1]*Udc[2] + Ud[4]*Udc[5] + Ud[7]*Udc[8]);
      sens[4] = Udc[2] + Udc[6] + (Udc[0]*Ud[2] + Udc[3]*Ud[5] + Udc[6]*Ud[8] +
				   Ud[0]*Udc[2] + Ud[3]*Udc[5] + Ud[6]*Udc[8]);
      sens[5] = Udc[1] + Udc[3] + (Udc[0]*Ud[1] + Udc[3]*Ud[4] + Udc[6]*Ud[7] +
				   Ud[0]*Udc[1] + Ud[3]*Udc[4] + Ud[6]*Udc[7]);
      sens += 6;

      Na++; Nb++; Nc++;
    }
  }
}

/*
  Compute the residual contribution from this element not including
  any externally applied loads.

  output:
  res:     the element residual
  
  input:
  vars:    the element variables
  Xpts:    the element nodal locations in R^{3}
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getRes( TacsScalar * res, 
				       const TacsScalar vars[], 
				       const TacsScalar Xpts[] ){
  memset(res, 0, NUM_VARIABLES*sizeof(TacsScalar));

  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  
  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  // Get the number of quadrature points
  int numGauss = getNumQuadraturePoints();

  for ( int n = 0; n < numGauss; n++ ){
    // Retrieve the quadrature points and weight
    double pt[3];
    double weight = getQuadraturePoint(n, pt);

    // Compute the element shape functions
    getShapeFunctions(pt, N, Na, Nb, Nc);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[9];
    solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[9];
    TacsScalar h = FElibrary::jacobian3d(Xa, J);
    h = h*weight;

    // Compute the strain
    TacsScalar strain[NUM_STRESSES];
    getStrain(strain, J, Na, Nb, Nc, vars);
 
    // Compute the corresponding stress
    TacsScalar stress[NUM_STRESSES];
    stiff->calculateStress(pt, strain, stress);
       
    // Get the derivative of the strain with respect to the nodal
    // displacements
    getBmat(B, J, Na, Nb, Nc, vars);

    TacsScalar *b = B;
    for ( int i = 0; i < NUM_VARIABLES; i++ ){
      res[i] += h*(b[0]*stress[0] + b[1]*stress[1] + b[2]*stress[2] +
		   b[3]*stress[3] + b[4]*stress[4] + b[5]*stress[5]);
      b += 6;
    }
  }
}

/*
  Get the element tangent stiffness matrix - the exact Jacobian of the
  residual expressions.

  output:
  mat:     the element tangent stiffness matrix
  res:     the element residual
  
  input:
  vars:    the element variables
  Xpts:    the element nodal locations in R^{3}
  matOr:   the matrix orientation (NORMAL or TRANSPOSE)
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getMat( TacsScalar *mat, 
				       TacsScalar *res, 
				       const TacsScalar vars[], 
				       const TacsScalar Xpts[],
				       MatrixOrientation matOr ){
  memset(mat, 0, NUM_VARIABLES*NUM_VARIABLES*sizeof(TacsScalar));
  memset(res, 0, NUM_VARIABLES*sizeof(TacsScalar));

  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  
  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  // Get the number of quadrature points
  int numGauss = getNumQuadraturePoints();

  for ( int n = 0; n < numGauss; n++ ){
    // Retrieve the quadrature points and weight
    double pt[3];
    double weight = getQuadraturePoint(n, pt);

    // Compute the element shape functions
    getShapeFunctions(pt, N, Na, Nb, Nc);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[9];
    solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[9];
    TacsScalar h = FElibrary::jacobian3d(Xa, J);
    h = h*weight;

    // Compute the strain
    TacsScalar strain[NUM_STRESSES];
    getStrain(strain, J, Na, Nb, Nc, vars);
 
    // Compute the corresponding stress
    TacsScalar stress[NUM_STRESSES];
    stiff->calculateStress(pt, strain, stress);
       
    // Add the stress times the second derivative of the strain
    addGeoStiffness(mat, h, stress, J, Na, Nb, Nc);

    // Get the derivative of the strain with respect to the nodal
    // displacements
    getBmat(B, J, Na, Nb, Nc, vars);

    TacsScalar *b = B;
    for ( int i = 0; i < NUM_VARIABLES; i++ ){
      res[i] += h*(b[0]*stress[0] + b[1]*stress[1] + b[2]*stress[2] +
		   b[3]*stress[3] + b[4]*stress[4] + b[5]*stress[5]);
      b += 6;
    }

    // Fill-in the upper portion of the matrix
    TacsScalar *bj = B;
    for ( int j = 0; j < NUM_VARIABLES; j++ ){
      // Compute the stress at the given point
      TacsScalar bs[NUM_STRESSES];
      stiff->calculateStress(pt, bj, bs);

      TacsScalar *bi = B;
      for ( int i = 0; i <= j; i++ ){
	mat[i + j*NUM_VARIABLES] += h*(bi[0]*bs[0] + bi[1]*bs[1] + bi[2]*bs[2] + 
				       bi[3]*bs[3] + bi[4]*bs[4] + bi[5]*bs[5]);
	bi += 6;
      }

      bj += 6;
    }
  }

  // Apply symmetry to the matrix
  for ( int j = 0; j < NUM_VARIABLES; j++ ){
    for ( int i = 0; i < j; i++ ){
      mat[j + i*NUM_VARIABLES] = mat[i + j*NUM_VARIABLES];
    }
  }
}

/*
  Evaluate the derivative of the element residuals with respect to the
  nodal coordinates e.g res = dR/dXpts

  output:
  res:  the derivative of the residuals w.r.t. the element nodes
  
  input:
  vars:    the element variables
  Xpts:    the element nodal locations
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getResXptSens( TacsScalar * res, 
					      const TacsScalar vars[],
					      const TacsScalar Xpts[] ){
  memset(res, 0, 3*NUM_NODES*NUM_VARIABLES*sizeof(TacsScalar));
}

/*
  Compute the derivative of the residuals with respect to one of the
  material design variables.

  output:
  res:    the derivative of the residual w.r.t. the material design var
  
  input:
  dvNum:   the design variable number 
  vars:    the variables
  Xpts:    the element nodal locations
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getResDVSens( int dvNum, 
					     TacsScalar * res, 
					     const TacsScalar vars[], 
					     const TacsScalar Xpts[] ){
  memset(res, 0, NUM_VARIABLES*sizeof(TacsScalar));

  if (dvNum < 0 || !stiff->ownsDesignVar(dvNum)){
    return;
  }

  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  
  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  // Get the number of quadrature points
  int numGauss = getNumQuadraturePoints();

  for ( int n = 0; n < numGauss; n++ ){
    // Retrieve the quadrature points and weight
    double pt[3];
    double weight = getQuadraturePoint(n, pt);

    // Compute the element shape functions
    getShapeFunctions(pt, N, Na, Nb, Nc);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[9];
    solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[9];
    TacsScalar h = FElibrary::jacobian3d(Xa, J);
    h = h*weight;

    // Compute the strain
    TacsScalar strain[NUM_STRESSES];
    getStrain(strain, J, Na, Nb, Nc, vars);
 
    // Compute the corresponding stress
    TacsScalar stress[NUM_STRESSES];
    stiff->calculateStressDVSens(dvNum, pt, strain, stress);
       
    // Get the derivative of the strain with respect to the nodal
    // displacements
    getBmat(B, J, Na, Nb, Nc, vars);

    TacsScalar *b = B;
    for ( int i = 0; i < NUM_VARIABLES; i++ ){
      res[i] += h*(b[0]*stress[0] + b[1]*stress[1] + b[2]*stress[2] +
		   b[3]*stress[3] + b[4]*stress[4] + b[5]*stress[5]);
      b += 6;
    }
  }
}

/*
  Add the product of the adjoint vector times the derivative of the
  residuals multiplied by a scalar to the given derivative vector.

  output:
  mat:     the element tangent stiffness matrix
  res:     the element residual
  
  input:
  psi:     the element adjoint variables
  vars:    the element variables
  Xpts:    the element nodal locations in R^{3}
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::addAdjResDVSensProduct( TacsScalar alpha,
						       TacsScalar dvSens[], 
						       int dvLen,
						       const TacsScalar psi[],
						       const TacsScalar vars[],
						       const TacsScalar Xpts[] ){
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  
  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  // Get the number of quadrature points
  int numGauss = getNumQuadraturePoints();

  for ( int n = 0; n < numGauss; n++ ){
    // Retrieve the quadrature points and weights
    double pt[3];
    double weight = getQuadraturePoint(n, pt);

    // Compute the element shape functions
    getShapeFunctions(pt, N, Na, Nb, Nc);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[9];
    solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[9];
    TacsScalar h = FElibrary::jacobian3d(Xa, J);
    h = h*weight;

    // Compute the strain
    TacsScalar strain[NUM_STRESSES];
    getStrain(strain, J, Na, Nb, Nc, vars);
        
    // Get the derivative of the strain with respect to the nodal
    // displacements
    getBmat(B, J, Na, Nb, Nc, vars);

    // Compute the product of psi^{T}*B^{T}
    TacsScalar bpsi[NUM_STRESSES];
    memset(bpsi, 0, NUM_STRESSES*sizeof(TacsScalar));

    TacsScalar *b = B;
    const TacsScalar *ps = psi;
    for ( int i = 0; i < NUM_VARIABLES; i++ ){
      bpsi[0] += ps[0]*b[0];
      bpsi[1] += ps[0]*b[1];
      bpsi[2] += ps[0]*b[2];
      bpsi[3] += ps[0]*b[3];
      bpsi[4] += ps[0]*b[4];
      bpsi[5] += ps[0]*b[5];
      b += 6;
      ps++;
    }
    
    // Add the term: alpha*psi^{T}*B^{T}*dC/dx*strain to the vector
    // dvSens - Note that this is much more efficient than computing
    // the terms component by component
    stiff->addStressDVSens(pt, strain, alpha*h, bpsi, dvSens, dvLen);
  }
}

/*
  Add the derivative of the inner product of the stiffness or mass
  matrix with respect to the design variables to a design variable
  vector. This is much more efficient than computing the derivative of
  the stiffness/mass matrix, then computing the product for each 
  design variable.

  input:
  matType:     the matrix type (e.g. MASS_MATRIX)
  alpha:       the scaling factor
  dvLen:       the length of the design variable vector
  psi:         the left inner-product vector
  phi:         the right inner-product vector
  vars:        the state variable values
  Xpts:        the nodal locations

  output:
  dvSens:      vector of the design sensitivity
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::addMatDVSensInnerProduct( ElementMatrixTypes matType,
							 TacsScalar alpha,
							 TacsScalar dvSens[], int dvLen,
							 const TacsScalar psi[], 
							 const TacsScalar phi[],
							 const TacsScalar vars[],
							 const TacsScalar Xpts[] ){
  if (matType == STIFFNESS_MATRIX){
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
    
    // The derivative of the stress with respect to the strain
    TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

    // Get the number of quadrature points
    int numGauss = getNumQuadraturePoints();

    for ( int n = 0; n < numGauss; n++ ){
      // Retrieve the quadrature points and weights
      double pt[3];
      double weight = getQuadraturePoint(n, pt);
      
      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb, Nc);
      
      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);
      
      // Compute the determinant of Xa and the transformation
      TacsScalar J[9];
      TacsScalar h = FElibrary::jacobian3d(Xa, J);
      h = h*weight;
      
      // Compute the strain
      TacsScalar strain[NUM_STRESSES];
      getStrain(strain, J, Na, Nb, Nc, vars);
      
      // Get the derivative of the strain with respect to the nodal
      // displacements
      getBmat(B, J, Na, Nb, Nc, vars);
      
      // Compute the product of psi^{T}*B^{T}
      TacsScalar bpsi[NUM_STRESSES], bphi[NUM_STRESSES];
      memset(bpsi, 0, NUM_STRESSES*sizeof(TacsScalar));
      memset(bphi, 0, NUM_STRESSES*sizeof(TacsScalar));
      
      TacsScalar *b = B;
      const TacsScalar *ps = psi, *ph = phi;
      for ( int i = 0; i < NUM_VARIABLES; i++ ){
	bpsi[0] += ps[0]*b[0];  bpsi[1] += ps[0]*b[1];
	bpsi[2] += ps[0]*b[2];  bpsi[3] += ps[0]*b[3];
	bpsi[4] += ps[0]*b[4];  bpsi[5] += ps[0]*b[5];

	bphi[0] += ph[0]*b[0];  bphi[1] += ph[0]*b[1];
	bphi[2] += ph[0]*b[2];  bphi[3] += ph[0]*b[3];
	bphi[4] += ph[0]*b[4];  bphi[5] += ph[0]*b[5];

	b += 6;
	ps++;  ph++;
      }

      // Add the result to the design variable vector
      stiff->addStressDVSens(pt, bphi, alpha*h, bpsi, dvSens, dvLen);
    }
  }
  else if (matType == MASS_MATRIX){
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  
    // Get the number of quadrature points
    int numGauss = getNumQuadraturePoints();
    
    for ( int n = 0; n < numGauss; n++ ){
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getQuadraturePoint(n, pt);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);
      
      // Compute the determinant of Xa and the transformation
      TacsScalar J[9];
      TacsScalar h = FElibrary::jacobian3d(Xa, J);
      h = h*weight;
      
      // Compute the nodal accelerations at the quadrature point
      TacsScalar upsi[3], uphi[3];
      upsi[0] = upsi[1] = upsi[2] = 0.0;
      uphi[0] = uphi[1] = uphi[2] = 0.0;

      double *ns = N;
      const TacsScalar *ps = psi, *ph = phi;
      for ( int i = 0; i < NUM_NODES; i++ ){
	upsi[0] += ns[0]*ps[0];
	upsi[1] += ns[0]*ps[1];  
	upsi[2] += ns[0]*ps[2];
	
	uphi[0] += ns[0]*ph[0];  
	uphi[1] += ns[0]*ph[1];  
	uphi[2] += ns[0]*ph[2];
	
	ps += 3; ph += 3; ns++;
      }

      // Add the result to the design variable vector
      TacsScalar rho_alpha = alpha*h*(upsi[0]*uphi[0] + upsi[1]*uphi[1] + 
				      upsi[2]*uphi[2]);
      stiff->addPointwiseMassDVSens(pt, &rho_alpha, dvSens, dvLen);
    }
  }
}

/*
  Get the element matrix of the specified type (e.g. mass matrix)
  from the element. 

  output:
  mat:         the element matrix of the specified type 

  input:
  matType:     the matrix type (e.g. MASS_MATRIX)
  scaleFactor: scale factor such that mat = scaleFactor*M
  vars:        the element variables
  Xpts:        the nodal coordinates in R^{3}
  matOr:       the matrix orientation either NORMAL or TRANSPOSE
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getMatType( ElementMatrixTypes matType, 
					   TacsScalar scaleFactor, 
					   TacsScalar * mat, 
					   const TacsScalar vars[], 
					   const TacsScalar Xpts[], 
					   MatrixOrientation matOr ){
  memset(mat, 0, NUM_VARIABLES*NUM_VARIABLES*sizeof(TacsScalar));

  // The mass matrix
  if (matType == MASS_MATRIX){
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  
    // Get the number of quadrature points
    int numGauss = getNumQuadraturePoints();
    
    for ( int n = 0; n < numGauss; n++ ){
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getQuadraturePoint(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb, Nc);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);
      
      // Compute the determinant of Xa and the transformation
      TacsScalar J[9];
      TacsScalar h = FElibrary::jacobian3d(Xa, J);
      h = h*weight;
      
      // Get the pointwise mass
      TacsScalar ptmass[3];
      stiff->pointwiseMass(pt, ptmass);

      // Fill-in the upper-portion of the matrix
      for ( int j = 0; j < NUM_NODES; j++ ){
	for ( int i = 0; i <= j; i++ ){
	  double d = h*ptmass[0]*N[i]*N[j];

	  mat[3*i + 3*j*NUM_VARIABLES] += d;
	  mat[3*i+1 + (3*j+1)*NUM_VARIABLES] += d;
	  mat[3*i+2 + (3*j+2)*NUM_VARIABLES] += d;
	}
      }
    }

    // Apply symmetry to the matrix
    for ( int j = 0; j < NUM_VARIABLES; j++ ){
      for ( int i = 0; i < j; i++ ){
	mat[j + i*NUM_VARIABLES] = mat[i + j*NUM_VARIABLES];
      }
    }
  }
}

/*
  Compute the derivative of the element matrix of the specified type

  output:
  mat:         the element matrix of the specified type 

  input:
  dvNum:       the derivative number to take
  matType:     the matrix type (e.g. MASS_MATRIX)
  scaleFactor: scale factor such that mat = scaleFactor*M
  vars:        the element variables
  Xpts:        the nodal coordinates in R^{3}
  matOr:       the matrix orientation either NORMAL or TRANSPOSE
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getMatTypeDVSens( int dvNum, 
						 ElementMatrixTypes matType, 
						 TacsScalar scaleFactor, 
						 TacsScalar * mat, 
						 const TacsScalar vars[],
						 const TacsScalar Xpts[], 
						 MatrixOrientation matOr ){
  memset(mat, 0, NUM_VARIABLES*NUM_VARIABLES*sizeof(TacsScalar));

  // The mass matrix
  if (dvNum >= 0 && stiff->ownsDesignVar(dvNum)){
    if (matType == MASS_MATRIX){
      
      // The shape functions associated with the element
      double N[NUM_NODES];
      double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  
      // Get the number of quadrature points
      int numGauss = getNumQuadraturePoints();
    
      for ( int n = 0; n < numGauss; n++ ){
	// Retrieve the quadrature points and weight
	double pt[3];
	double weight = getQuadraturePoint(n, pt);

	// Compute the element shape functions
	getShapeFunctions(pt, N, Na, Nb, Nc);

	// Compute the derivative of X with respect to the
	// coordinate directions
	TacsScalar X[3], Xa[9];
	solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);
      
	// Compute the determinant of Xa and the transformation
	TacsScalar J[9];
	TacsScalar h = FElibrary::jacobian3d(Xa, J);
	h = h*weight;
      
	// Get the pointwise mass
	TacsScalar ptmass[3];
	stiff->pointwiseMassDVSens(dvNum, pt, ptmass);

	// Fill-in the upper-portion of the matrix
	for ( int j = 0; j < NUM_NODES; j++ ){
	  for ( int i = 0; i <= j; i++ ){
	    double d = h*ptmass[0]*N[i]*N[j];

	    mat[3*i + 3*j*NUM_VARIABLES] += d;
	    mat[3*i+1 + (3*j+1)*NUM_VARIABLES] += d;
	    mat[3*i+2 + (3*j+2)*NUM_VARIABLES] += d;
	  }
	}
      }

      // Apply symmetry to the matrix
      for ( int j = 0; j < NUM_VARIABLES; j++ ){
	for ( int i = 0; i < j; i++ ){
	  mat[j + i*NUM_VARIABLES] = mat[i + j*NUM_VARIABLES];
	}
      }
    }
  }
}

/*
  Evaluate the determinant of the Jacobian for numerical integration

  returns: the determinant of the Jacobian
 
  input:
  pt:    the parametric point within the element
  Xpts:  the element nodes
*/
template <int NUM_NODES>
TacsScalar TACS3DElement<NUM_NODES>::getJacobian( const double pt[], 
						  const TacsScalar Xpts[] ){
  // Compute the element shape functions
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of the shape functions w.r.t. the 
  // parametric locations
  TacsScalar X[3], Xa[9];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

  return FElibrary::jacobian3d(Xa);
}

/*
  Evaluate the derivative of the determinant of the Jacobian with
  respect to the element nodal locations
  
  output:
  hXptSens: the derivative of the determinant w.r.t. the nodal locations

  returns:  the determinant of the Jacobian
 
  input:
  pt:    the parametric point within the element
  Xpts:  the element nodes
*/
template <int NUM_NODES>
TacsScalar TACS3DElement<NUM_NODES>::getJacobianXptSens( TacsScalar * hXptSens, 
							 const double * pt,
							 const TacsScalar Xpts[] ){
  // Compute the element shape functions
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of the shape functions w.r.t. the 
  // parametric locations
  TacsScalar X[3], Xa[9];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

  // Evaluate the determinant of the Jacobian
  TacsScalar J[9];
  TacsScalar h = FElibrary::jacobian3d(Xa, J);
  
  for ( int i = 0; i < NUM_NODES; i++ ){
    for ( int k = 0; k < 3; k++ ){
      TacsScalar XaSens[9];
      XaSens[0] = XaSens[1] = XaSens[2] = 0.0;
      XaSens[3] = XaSens[4] = XaSens[5] = 0.0;
      XaSens[6] = XaSens[7] = XaSens[8] = 0.0;      
      XaSens[k] = Na[i];
      XaSens[k+3] = Nb[i];
      XaSens[k+6] = Nc[i];

      FElibrary::jacobian3dSens(Xa, XaSens, &hXptSens[0]);
      hXptSens++;
    }
  }

  return h;
}

/*
  Evaluate the strain at the specified point using the provided set of
  variables

  output:
  strain:   the strain evaluate at the specific parametric point
  
  input:
  vars:     the element variable values
  Xpts:     the element nodal locations
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getPtwiseStrain( TacsScalar strain[], 
						const double pt[],
						const TacsScalar vars[], 
						const TacsScalar Xpts[] ){
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

  // Compute the element shape functions
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of X with respect to the coordinate directions
  TacsScalar X[3], Xa[9];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

  // Compute the determinant of Xa and the transformation
  TacsScalar J[9];
  FElibrary::jacobian3d(Xa, J);

  // Compute the strain
  getStrain(strain, J, Na, Nb, Nc, vars);
}

/*
  Compute the derivative of the point-wise strain multiplied by a
  specified vector with respect to the element variables and add the
  result to an array. This can be used to evaluate the derivative of a
  function of interest with respect to the element variables.

  output:
  sens:        the output array - same length as the number of elem variables
 
  input:
  pt:          parametric point used to evaluate the derivative [-1, 1]^{3}
  scaleFactor: scale ther result by this scalar
  strainSens:  the sensitivity of each straint component 
  vars:        the element variables
  Xpts:        the element nodal locations
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::addPtwiseStrainSVSens( TacsScalar sens[], 
						      const double pt[], 
						      const TacsScalar scaleFactor, 
						      const TacsScalar strainSens[],
						      const TacsScalar vars[], 
						      const TacsScalar Xpts[] ){
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  // Compute the element shape functions
  getShapeFunctions(pt, N, Na, Nb, Nc);
  
  // Compute the derivative of X with respect to the coordinate
  // directions
  TacsScalar X[3], Xa[9];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);
  
  // Compute the determinant of Xa and the transformation
  TacsScalar J[9];
  FElibrary::jacobian3d(Xa, J);
        
  // Get the derivative of the strain with respect to the nodal
  // displacements
  getBmat(B, J, Na, Nb, Nc, vars);

  TacsScalar *b = B;
  for ( int i = 0; i < NUM_VARIABLES; i++ ){
    sens[i] += scaleFactor*(strainSens[0]*b[0] + strainSens[1]*b[1] + 
			    strainSens[2]*b[2] + strainSens[3]*b[3] +
			    strainSens[4]*b[4] + strainSens[5]*b[5]);
    b += 6;
  }
}

/*
  Compute the strain and the derivative of the strain with respect to
  the nodal locations.

  output:
  strain:        the strain evaluate at the pamametric point pt
  strainXptSens: the derivative of the straint w.r.t. the nodal locations

  input:
  pt:        the parametric point within the element
  vars:      the element variables
  Xpts:      the nodal locations
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getPtwiseStrainXptSens( TacsScalar strain[], 
						       TacsScalar strainXptSens[], 
						       const double * pt,
						       const TacsScalar vars[], 
						       const TacsScalar Xpts[] ){

  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  // Compute the element shape functions
  getShapeFunctions(pt, N, Na, Nb, Nc);
  
  // Compute the derivative of X with respect to the coordinate
  // directions
  TacsScalar X[3], Xa[9];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);
  
  // Compute the determinant of Xa and the transformation
  TacsScalar J[9];
  FElibrary::jacobian3d(Xa, J);

  // Compute the strain
  getStrain(strain, J, Na, Nb, Nc, vars);

  // Compute the derivative of the strain w.r.t. nocal coordinates
  getStrainXptSens(strainXptSens, J, Xa, 
		   Na, Nb, Nc, vars);
}

/*
  Set the rigid-moment displacement extrapolation flag. Only
  extrapolate the displacements and rotations if this flag is
  activated.
*/
template <int NUM_NODES>
int TACS3DElement<NUM_NODES>::USE_RIGID_MOMENT = 1;

/*
  Set the rigid moment flag for all TACS3DElements. This flag can be
  turned either on-or off (since it is a static member it is the same
  for all TACS3DElements - (on this processor!))
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::setRigidMomentFlag( int _flag ){
  USE_RIGID_MOMENT = _flag;
}

/*
  Get the extrapolation of the displacement of a point in R^{3} where
  the extrapolation may include a rigid link term that accounts for
  the rotational displacement of the body. Note that the displacement
  may not lie directly on the surface of the structure, therefore
  displacement extrapolation is required.

  output:
  Uaero:  the displacement of the aerodynamic point
  
  input:
  Xa:     the location of the aerodynamic point in the undeformed config.
  pt:     the parametric location of the attachment point
  vars:   the element variable values
  Xpts:   the element nodal locations
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getRD( TacsScalar Uaero[], 
				      const TacsScalar Xaero[], 
				      const double pt[], 
				      const TacsScalar vars[], 
				      const TacsScalar Xpts[] ){
  // Compute the element shape functions
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of the position w.r.t. the parametric
  // coordinates 
  TacsScalar X[3], Xa[3];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);
  
  // Compute the inverse of the transformation Xa
  TacsScalar J[9];
  FElibrary::jacobian3d(Xa, J);

  // Compute the displacement gradient
  TacsScalar Ud[9];
  getDeformGradient(Ud, J, Na, Nb, Nc, vars);

  // Compute the displacement
  TacsScalar U[3] = {0.0, 0.0, 0.0};
  double *n = N;
  const TacsScalar *u = vars; 
  for ( int i = 0; i < NUM_NODES; i++ ){
    U[0] += n[0]*u[0];
    U[1] += n[0]*u[1];
    U[2] += n[0]*u[2];
    u += 3; n++;
  }

  // Compute the rigid
  TacsScalar R[3];
  R[0] = Xaero[0] - X[0];
  R[1] = Xaero[1] - X[1];
  R[2] = Xaero[2] - X[2];

  // Compute the rotation: The skew symmetric part of the
  // displacement gradient
  TacsScalar rot[3];
  rot[0] = 0.5*(Ud[5] - Ud[7]);
  rot[1] = 0.5*(Ud[6] - Ud[2]);
  rot[2] = 0.5*(Ud[1] - Ud[3]);

  // Compute the displacements at the aerodynamic point:
  // Uaero = U + cross(rot, R)
  Uaero[0] = U[0] + rot[1]*R[2] - rot[2]*R[1];
  Uaero[1] = U[1] + rot[2]*R[0] - rot[0]*R[2];
  Uaero[2] = U[2] + rot[0]*R[1] - rot[1]*R[0];
}

/*
  Compute the derivative of the extrapolated displacement with respect
  to the element nodal locations. This is the derivative of the code
  getRD displacements with respect to the input Xpts.  This
  computation is performed against an input vector XptSens which
  provides the sensitivities of nodal locations w.r.t. the specified
  derivative

  output:
  Uaero:   the displacement of the aerodynamic point
  
  input:
  Xa:      the location of the aerodynamic point in the undeformed config.
  XaSens:  the derivative of Xa with respect to the nodal locations
  pt:      the parametric location of the attachment point
  vars:    the element variable values
  Xpts:    the element nodal locations
  XptSens: the sensitivity of the nodal locations
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getRDXptSens( TacsScalar UaeroSens[], 
					     const TacsScalar Xaero[], 
					     const TacsScalar XaeroSens[], 
					     const double pt[], 
					     const TacsScalar vars[], 
					     const TacsScalar Xpts[], 
					     const TacsScalar XptSens[] ){ 
  // Compute the element shape functions
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of the position w.r.t. the parametric
  // coordinates 
  TacsScalar X[3], Xa[3];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

  // Compute the derivatives of the provided positions at the parametric
  // coordinate specified
  TacsScalar XSens[3], XaSens[3];
  solidJacobian(XSens, XaSens, N, Na, Nb, Nc, Xpts);

  // Compute the inverse of the transformation Xa and the derivative
  // of the inverse transformation
  TacsScalar sh;
  TacsScalar J[9], JSens[9];
  FElibrary::jacobian3dSens(Xa, XaSens, J, JSens, &sh);
  
  // Get the displacement gradient and its derivative
  TacsScalar Ud[9], UdSens[9];
  getDeformGradientSens(Ud, UdSens, 
			J, JSens, Na, Nb, Nc, vars);

  // Compute the rigid link vector and its sensitivity
  TacsScalar R[3], RSens[3];
  R[0] = Xaero[0] - X[0];
  R[1] = Xaero[1] - X[1];
  R[2] = Xaero[2] - X[2];
  RSens[0] = Xaero[0] - XSens[0];
  RSens[1] = Xaero[1] - XSens[1];
  RSens[2] = Xaero[2] - XSens[2];

  TacsScalar rot[3], rotSens[3];
  rot[0] = 0.5*(Ud[5] - Ud[7]);
  rot[1] = 0.5*(Ud[6] - Ud[2]);
  rot[2] = 0.5*(Ud[1] - Ud[3]);
  rotSens[0] = 0.5*(UdSens[5] - UdSens[7]);
  rotSens[1] = 0.5*(UdSens[6] - UdSens[2]);
  rotSens[2] = 0.5*(UdSens[1] - UdSens[3]);

  UaeroSens[0] = (rotSens[1]*R[2] - rotSens[2]*R[1] + 
		  rot[1]*RSens[2] - rot[2]*RSens[1]);
  UaeroSens[1] = (rotSens[2]*R[0] - rotSens[0]*R[2] + 
		  rot[2]*RSens[0] - rot[0]*RSens[2]);
  UaeroSens[2] = (rotSens[0]*R[1] - rotSens[1]*R[0] + 
		  rot[0]*RSens[1] - rot[1]*RSens[0]);
}

/*
  Compute the product of the transpose of the derivative of the rigid
  extrapolation with an input adjoint vector

  output:
  elemAdj:  the transpose product of the derivative w.r.t. the adjoint input
  
  input:
  Uaero:    the extrapolated displacement adjoint sensitivity
  Xa:       the location of the aerodynamic point in the undeformed config.
  pt:       the parametric location of the attachment point
  Xpts:     the element nodal locations
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getRDTranspose( TacsScalar elemAdj[], 
					       const TacsScalar Uaero[], 
					       const TacsScalar Xaero[], 
					       const double pt[], 
					       const TacsScalar Xpts[] ){
  // Compute the element shape functions
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of the position w.r.t. the parametric
  // coordinates 
  TacsScalar X[3], Xa[3];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);
  
  // Compute the inverse of the transformation Xa
  TacsScalar J[9];
  FElibrary::jacobian3d(Xa, J);

  // Compute the rigid link between the aerodynamic point and the
  // structural surface
  TacsScalar R[3];
  R[0] = Xaero[0] - X[0];
  R[1] = Xaero[1] - X[1];
  R[2] = Xaero[2] - X[2];

  double *n = N, *na = Na, *nb = Nb, *nc = Nc;
  for ( int i = 0; i < NUM_NODES; i++ ){
    // Compute the derivative of the displacement gradient
    // with respect to the nodal variables
    TacsScalar Dn[9];
    Dn[0] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
    Dn[3] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
    Dn[6] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
    
    Dn[1] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
    Dn[4] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
    Dn[7] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
    
    Dn[2] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
    Dn[5] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
    Dn[8] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];

    // Compute the derivative of rot w.r.t. the u-variables
    TacsScalar rot[3];
    rot[0] = 0.0;
    rot[1] = -0.5*Dn[2];
    rot[2] =  0.5*Dn[1];
    elemAdj[0] = (Uaero[0]*n[0] +
		  Uaero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
		  Uaero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
		  Uaero[2]*(rot[0]*R[1] - rot[1]*R[0]));

    rot[0] =  0.5*Dn[5];
    rot[1] =  0.0;
    rot[2] = -0.5*Dn[3];
    elemAdj[1] = (Uaero[1]*n[0] +
		  Uaero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
		  Uaero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
		  Uaero[2]*(rot[0]*R[1] - rot[1]*R[0]));

    rot[0] = -0.5*Dn[7];
    rot[1] =  0.5*Dn[6];
    rot[2] = 0.0;
    elemAdj[2] = (Uaero[2]*n[0] +
		  Uaero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
		  Uaero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
		  Uaero[2]*(rot[0]*R[1] - rot[1]*R[0]));
    
    elemAdj += 3;
    n++; na++; nb++; nc++;
  }
}

/*
  Add the consistent force applied to the element due to a point load
  applied from a rigid link. Note that the work due to the aerodynamic
  force through the rigid link can be written as:

  Work =
  Faero[0]*(U[0] + rot[1]*R[2] - rot[2]*R[1] ) +
  Faero[1]*(U[1] + rot[2]*R[0] - rot[0]*R[2] ) +
  Faero[2]*(U[2] + rot[0]*R[1] - rot[1]*R[0] )

  output:
  res:   the consistent force vector

  input:
  Faero:   the applied force
  Xaero:   the point of application
  pt:      the parametric location
  Xpts:    the element nodal locations
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getRF( TacsScalar res[],
				      const TacsScalar Faero[],
				      const TacsScalar Xaero[],
				      const double pt[], 
				      const TacsScalar Xpts[] ){
  // Compute the element shape functions
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of the position w.r.t. the parametric
  // coordinates 
  TacsScalar X[3], Xa[3];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);
  
  // Compute the inverse of the Jacobian
  TacsScalar J[9];
  FElibrary::jacobian3d(Xa, J);

  // Compute the rigid link
  TacsScalar R[3];
  R[0] = Xaero[0] - X[0];
  R[1] = Xaero[1] - X[1];
  R[2] = Xaero[2] - X[2];

  // Compute the contribution from each node
  double *n = N, *na = Na, *nb = Nb, *nc = Nc;
  for ( int i = 0; i < NUM_NODES; i++ ){
    // Compute the derivative of the displacement gradient
    // with respect to the nodal variables
    if (USE_RIGID_MOMENT){
      TacsScalar Dn[9];
      Dn[0] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      Dn[3] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      Dn[6] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      
      Dn[1] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      Dn[4] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      Dn[7] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      
      Dn[2] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
      Dn[5] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
      Dn[8] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
      
      TacsScalar rot[3];
      rot[0] = 0.0;
      rot[1] = -0.5*Dn[2];
      rot[2] =  0.5*Dn[1];
      res[0] = - (Faero[0]*n[0] +
		  Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
		  Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
		  Faero[2]*(rot[0]*R[1] - rot[1]*R[0]));
      
      rot[0] =  0.5*Dn[5];
      rot[1] =  0.0;
      rot[2] = -0.5*Dn[3];
      res[1] = - (Faero[1]*n[0] + 
		  Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
		  Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
		  Faero[2]*(rot[0]*R[1] - rot[1]*R[0]));
      
      rot[0] = -0.5*Dn[7];
      rot[1] =  0.5*Dn[6];
      rot[2] = 0.0;
      res[2] = - (Faero[2]*n[0] +
		  Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
		  Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
		  Faero[2]*(rot[0]*R[1] - rot[1]*R[0]));
    }
    else {
      res[0] = -Faero[0]*n[0];
      res[1] = -Faero[1]*n[1];
      res[2] = -Faero[2]*n[2];
    }

    // Increment the res pointers
    res += 3;    
    n++; na++; nb++; nc++;
  }
}

/*
  Determine the consistent force applied to the element due to a 
  point load applied from a rigid link. 

  F  == the force vector
  Rd == the relative position vector Rd = Xa - Xelem 
  pt == the parametric location of the point
*/

template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getRFXptSens( TacsScalar res[], 
					     const TacsScalar Faero[], 
					     const TacsScalar Xaero[], 
					     const TacsScalar XaeroSens[],
					     const double pt[], 
					     const TacsScalar Xpts[], 
					     const TacsScalar XptSens[] ){
  if (USE_RIGID_MOMENT){
    // Compute the element shape functions
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
    getShapeFunctions(pt, N, Na, Nb, Nc);
    
    // Compute the derivative of the position w.r.t. the parametric
    // coordinates 
    TacsScalar X[3], Xa[3];
    solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);
    
    // Compute the sensitivity of the points
    TacsScalar XSens[3], XaSens[9];
    solidJacobian(XSens, XaSens, N, Na, Nb, Nc, XptSens);

    // Compute the derivativve of the Jacobian w.r.t. the
    // perturbation XptSens
    TacsScalar sh;
    TacsScalar J[9], JSens[9];
    FElibrary::jacobian3dSens(Xa, XaSens, J, JSens, &sh);
    
    TacsScalar R[3], RSens[3];
    R[0] = Xaero[0] - X[0];
    R[1] = Xaero[1] - X[1];
    R[2] = Xaero[2] - X[2];
    RSens[0] = XaeroSens[0] - XSens[0];
    RSens[1] = XaeroSens[1] - XSens[1];
    RSens[2] = XaeroSens[2] - XSens[2];
    
    // Compute the contribution from each node
    double *n = N, *na = Na, *nb = Nb, *nc = Nc;
    for ( int i = 0; i < NUM_NODES; i++ ){
      // Compute the derivative of the displacement gradient
      // with respect to the nodal variables
      TacsScalar Dn[9];
      Dn[0] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      Dn[3] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      Dn[6] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      
      Dn[1] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      Dn[4] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      Dn[7] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      
      Dn[2] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
      Dn[5] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
      Dn[8] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];

      TacsScalar DnSens[9];
      DnSens[0] = na[0]*JSens[0] + nb[0]*JSens[3] + nc[0]*JSens[6];
      DnSens[3] = na[0]*JSens[0] + nb[0]*JSens[3] + nc[0]*JSens[6];
      DnSens[6] = na[0]*JSens[0] + nb[0]*JSens[3] + nc[0]*JSens[6];
      
      DnSens[1] = na[0]*JSens[1] + nb[0]*JSens[4] + nc[0]*JSens[7];
      DnSens[4] = na[0]*JSens[1] + nb[0]*JSens[4] + nc[0]*JSens[7];
      DnSens[7] = na[0]*JSens[1] + nb[0]*JSens[4] + nc[0]*JSens[7];
      
      DnSens[2] = na[0]*JSens[2] + nb[0]*JSens[5] + nc[0]*JSens[8];
      DnSens[5] = na[0]*JSens[2] + nb[0]*JSens[5] + nc[0]*JSens[8];
      DnSens[8] = na[0]*JSens[2] + nb[0]*JSens[5] + nc[0]*JSens[8];
      
      TacsScalar rot[3], rotSens[3];
      rot[0] = 0.0;
      rot[1] = -0.5*Dn[2];
      rot[2] =  0.5*Dn[1];
      rotSens[0] = 0.0;
      rotSens[1] = -0.5*DnSens[2];
      rotSens[2] =  0.5*DnSens[1];
      res[0] = - (Faero[0]*(rotSens[1]*R[2] - rotSens[2]*R[1] +
			    rot[1]*RSens[2] - rot[2]*RSens[1]) +
		  Faero[1]*(rotSens[2]*R[0] - rotSens[0]*R[2] +
			    rot[2]*RSens[0] - rot[0]*RSens[2]) +
		  Faero[2]*(rotSens[0]*R[1] - rotSens[1]*R[0] +
			    rot[0]*RSens[1] - rot[1]*RSens[0]));
      
      rot[0] =  0.5*Dn[5];
      rot[1] =  0.0;
      rot[2] = -0.5*Dn[3];
      rotSens[0] =  0.5*DnSens[5];
      rotSens[1] =  0.0;
      rotSens[2] = -0.5*DnSens[3];
      res[1] = - (Faero[0]*(rotSens[1]*R[2] - rotSens[2]*R[1] +
			    rot[1]*RSens[2] - rot[2]*RSens[1]) +
		  Faero[1]*(rotSens[2]*R[0] - rotSens[0]*R[2] +
			    rot[2]*RSens[0] - rot[0]*RSens[2]) +
		  Faero[2]*(rotSens[0]*R[1] - rotSens[1]*R[0] +
			    rot[0]*RSens[1] - rot[1]*RSens[0]));
      
      rot[0] = -0.5*Dn[7];
      rot[1] =  0.5*Dn[6];
      rot[2] = 0.0;
      rotSens[0] = -0.5*DnSens[7];
      rotSens[1] =  0.5*DnSens[6];
      rotSens[2] = 0.0;
      res[2] = - (Faero[0]*(rotSens[1]*R[2] - rotSens[2]*R[1] +
			    rot[1]*RSens[2] - rot[2]*RSens[1]) +
		  Faero[1]*(rotSens[2]*R[0] - rotSens[0]*R[2] +
			    rot[2]*RSens[0] - rot[0]*RSens[2]) +
		  Faero[2]*(rotSens[0]*R[1] - rotSens[1]*R[0] +
			    rot[0]*RSens[1] - rot[1]*RSens[0]));
    }
  }
  else {
    memset(res, 0, NUM_VARIABLES*sizeof(TacsScalar));
  }
}

/*
  Determine the consistent force applied to the element due to a 
  point load applied from a rigid link. 

  F  == the force vector
  Xa == the aerodynamic surface node
  pt == the parametric location of the point
  res == the element residual
*/

template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getRFTranspose( TacsScalar Faero[], 
					       const TacsScalar Xaero[],
					       const double pt[], 
					       const TacsScalar res[], 
					       const TacsScalar Xpts[] ){
  // Compute the element shape functions
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of the position w.r.t. the parametric
  // coordinates 
  TacsScalar X[3], Xa[3];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

  // Compute the inverse of the transformation Xa
  TacsScalar J[9];
  FElibrary::jacobian3d(Xa, J);

  // Compute the displacement gradient
  TacsScalar Ud[9];
  getDeformGradient(Ud, J, Na, Nb, Nc, res);

  // Compute the displacement
  TacsScalar U[3] = {0.0, 0.0, 0.0};
  double *n = N;
  for ( int i = 0; i < NUM_NODES; i++ ){
    U[0] += n[0]*res[0];
    U[1] += n[0]*res[1];
    U[2] += n[0]*res[2];
    res += 3; n++;
  }

  // Compute the vector from the structural surface to the aerodynamic
  // point
  TacsScalar R[3];
  R[0] = Xaero[0] - X[0];
  R[1] = Xaero[1] - X[1];
  R[2] = Xaero[2] - X[2];

  if ( USE_RIGID_MOMENT){
    TacsScalar rot[3];
    rot[0] = 0.5*(Ud[5] - Ud[7]);
    rot[1] = 0.5*(Ud[6] - Ud[2]);
    rot[2] = 0.5*(Ud[1] - Ud[3]);
    
    Faero[0] = -(U[0] + rot[1]*R[2] - rot[2]*R[1]);
    Faero[1] = -(U[1] + rot[2]*R[0] - rot[0]*R[2]);
    Faero[2] = -(U[2] + rot[0]*R[1] - rot[1]*R[0]);
  }
  else {
    Faero[0] = -U[0];
    Faero[1] = -U[1];
    Faero[2] = -U[2];
  }
}

#endif
