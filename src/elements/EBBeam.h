#ifndef TACS_EULER_BERNOULLI_BEAM_H
#define TACS_EULER_BERNOULLI_BEAM_H

/*
  The Euler-Bernoulli beam element.

  Copyright (c) 2011 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.

  This element takes as input a beam stiffness object. This object
  includes the reference direction used to define the strong or weak
  axis of the beam. 
*/

#include "TACSElement.h"
#include "EBStiffness.h"
#include "KSM.h"

class EBBeam : public TACSElement {
 public:

  EBBeam( EBStiffness * _beamStiff);
  ~EBBeam();

  // Functions to access the design variables
  // ----------------------------------------
  int ownsDesignVar( const int dvNum ) const;
  int getNumDesignVars() const;
  int getDesignVarNums( int * dvNums, int * dvIndex, int dvLen ) const;
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs ) const;
  void getDesignVarRange( TacsScalar lowerBound[], 
			  TacsScalar upperBound[], int numDVs ) const;

  // Get the element properties and names
  // ------------------------------------
  const char * elementName() const;
  const char * displacementName( int i ) const;
  const char * stressName( int i ) const;
  const char * strainName( int i ) const;
  const char * extraName( int i ) const;
  int numDisplacements();
  int numStresses() const;
  int numNodes();
  int numVariables() const;
  int numExtras() const;
  enum ElementType getElementType();

  // Functions for analysis
  // ----------------------
  void getRes( TacsScalar * res, const TacsScalar vars[], 
	       const TacsScalar Xpts[] );
  void getMat( TacsScalar * mat, TacsScalar * res, 
	       const TacsScalar vars[], const TacsScalar Xpts[], 
	       MatrixOrientation matOr );
  void getMatType( ElementMatrixType matType, TacsScalar scaleFactor, 
		   TacsScalar * mat, const TacsScalar vars[], 
		   const TacsScalar Xpts[], MatrixOrientation matOr );

  void addResidual( double time, TacsScalar res[],
                            const TacsScalar Xpts[],
                            const TacsScalar vars[],
                            const TacsScalar dvars[],
                            const TacsScalar ddvars[] );

  // Functions required to determine the derivatives w.r.t. the design variables
  // ---------------------------------------------------------------------------
  void getResXptSens( TacsScalar * res, const TacsScalar vars[], 
		      const TacsScalar Xpts[] );  
  void getResDVSens( int dvNum, TacsScalar * res, 
		     const TacsScalar vars[], const TacsScalar Xpts[] );  

  // Functions used mostly for aerostructural coupling
  // -------------------------------------------------
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

  // Functions for evaluating global functionals of interest
  // -------------------------------------------------------
  TACSConstitutive * getConstitutive(){ return beamStiff; }
  int getGaussPtScheme();
  int getNumGaussPts( const int scheme = -1 );
  TacsScalar getGaussWtsPts( const int scheme, const int num, double * pt );

  TacsScalar getClosestPt( double pt[], int * fflag, 
			   const TacsScalar Xp[], const TacsScalar Xpts[] );
  void getPoint( TacsScalar Xp[], const double pt[], const TacsScalar Xpts[] );
  void getShapeFunctions( const double pt[], double N[] );
  
  TacsScalar getJacobian( const double * pt, const TacsScalar Xpts[] );
  TacsScalar getJacobianXptSens( TacsScalar hXptSens[], const double * pt, 
				 const TacsScalar Xpts[] );

  void getPtwiseStrain( TacsScalar strain[], const double * pt,
			const TacsScalar vars[],
			const TacsScalar Xpts[] );

  void getPtwiseStrainXptSens( TacsScalar strain[], TacsScalar strainSens[],
			       const double * pt, const TacsScalar vars[], 
			       const TacsScalar Xpts[] );
   
  void addPtwiseStrainSVSens( TacsScalar elementSens[], const double * pt, 
			      const TacsScalar scaleFactor, 
			      const TacsScalar strainSens[],
			      const TacsScalar vars[],
			      const TacsScalar Xpts[] );

  // Functions for post-processing
  // -----------------------------
  void addOutputCount( int * nelems, int * nnodes, int * ncsr );
  void getOutputData( unsigned int out_type, double * data, int ld_data,
		      const TacsScalar vars[], const TacsScalar Xpts[] );
  void getOutputConnectivity( int * con, int node );

  // Functions related to the ref axis
  enum EBBeamReferenceDirection getRefDirTyp(){return ref_dir_type;}
  const TacsScalar * getRefAxis(){ 
    return ref_dir;
  }

 private:

  // Compute the transformation between reference frames
  // ---------------------------------------------------
  TacsScalar computeTransform( TacsScalar t[],
                               const TacsScalar Xpts[] );

  TacsScalar computeTransformXptSens( TacsScalar t[],
                                      TacsScalar tSens[],
                                      TacsScalar * LSens,
                                      int component,
                                      const TacsScalar Xpts[] );

  // Transform the variables, residual and stiffness matrix
  // ------------------------------------------------------  
  void transformVarsGlobalToLocal( TacsScalar vars[],
                                   const TacsScalar t[] );
  void transformVarsGlobalToLocal( TacsScalar vars[],
                                   const TacsScalar elemVars[],
                                   const TacsScalar t[] );

  void transformResLocalToGlobal( TacsScalar res[],
                                  const TacsScalar t[] );
  void transformResLocalToGlobalSens( TacsScalar res[],
                                      const TacsScalar resSens[],
                                      const TacsScalar t[],
                                      const TacsScalar tSens[] );

  void transformStiffnessMat( TacsScalar mat[],
                              const TacsScalar t[] );

  inline TacsScalar strain_product( const TacsScalar a[], 
				    const TacsScalar b[] ){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
  }

  // Function pointers to the strain expressions
  void evalStrain( TacsScalar h, const double Na[], 
                   const double Nahp[], const double Naahp[],
                   const TacsScalar vars[], TacsScalar strain[] );
  void evalBvec( TacsScalar h, const double Na[], 
                 const double Nahp[], const double Naahp[],
                 TacsScalar strain[] );

  void evalStrainSens( TacsScalar h, TacsScalar hSens,
                       const double Na[], const double Nahp[], 
                       const double Naahp[],
                       const TacsScalar vars[], const TacsScalar varsSens[],
                       TacsScalar strainSens[] );
  void evalBvecSens( TacsScalar h, TacsScalar hSens, const double Na[], 
                     const double Nahp[], const double Naahp[],
                     TacsScalar strain[] );

  static const int NUM_DISPS = 6;
  static const int NUM_STRESSES = 4;
  static const int NUM_NODES = 2;
  static const int NUM_VARIABLES = 6*NUM_NODES;
  static const int NUM_EXTRAS = 2;

  int numGauss;
  const double * gaussWts;
  const double * gaussPts;

  static const char * elemName;
  static const char * dispNames[6];
  static const char * stressNames[4];
  static const char * strainNames[4];
  static const char * extraNames[2];

  // Information about the reference axis
  enum EBBeamReferenceDirection ref_dir_type;
  TacsScalar ref_dir[3];

  int stiffFlag; // Flag if the beamStiff dvs are set
  EBStiffness * beamStiff;
};

#endif

