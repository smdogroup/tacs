#ifndef TACS_DYNAMICS_ELEMENT_H
#define TACS_DYNAMICS_ELEMENT_H

/*
  Rigid-body dynamics routines for TACS

  Copyright (c) 2015-2016 Graeme Kennedy. All rights reserved. 
*/

#include "TACSObject.h"

/*
  Abstract interface for TACSDynElement. 

  This is an interface that enables flexible multibody dynamic
  analysis of finite-element models. The interface is designed for
  elements which are inherantly geometric/inertially nonlinear and may
  be coupled to other bodies via kinematic constraints.  

  The interface is also designed to enable the computation of
  derivatives of functions of interest. These derivatives may be used
  for gradient-based design optimization. 
*/
class TACSDynElement : public TACSObject {
 public:
  // Get the initial variable values
  // -------------------------------
  virtual void getInitVariables( TacsScalar q[] ) = 0;

  // Compute the kinetic and potential energies associated with the element
  // ----------------------------------------------------------------------
  virtual void computeEnergies( TacsScalar *_Te, 
				TacsScalar *_Pe,
				const TacsScalar X[],
				const TacsScalar vars[],
				const TacsScalar dvars[] ) = 0;

  // Compute the residual of the governing equations
  // -----------------------------------------------
  virtual void getResidual( TacsScalar res[],
			    const TacsScalar X[],
			    const TacsScalar vars[],
			    const TacsScalar dvars[],
			    const TacsScalar ddvars[] ) = 0;

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  virtual void getJacobian( TacsScalar J[],
			    double alpha, double beta, double gamma,
			    const TacsScalar X[],
			    const TacsScalar vars[],
			    const TacsScalar dvars[],
			    const TacsScalar ddvars[] ) = 0;
};

/*
  A Mixed-Interpolation of Tensorial Components element for dynamic
  analysis.

  The following class implements a geometrically nonlinear
  finite-element shell for large displacement/rotation problems.  The
  element permits arbitrary rotation/displacement rigid body motion.
  The rotational parametrization is based on the quaternions with an
  added constraint at each node.

  The shell formulation utilizes through-thickness strain/kinematic
  assumptions that are based on first-order shear deformation
  theory. The theory takes into account the nonlinear rotational
  kinematics that are required to obtain strain-free rotation of the
  elements.

  The drilling degree of freedom is handled through the use of a penalty
  term that penalizes the discrepancy between the in-plane rotations
  predicted from nonlinear shell theory and those predicted by stress-
*/
class MITC9 : public TACSDynElement {
 public:
  static const int ORDER = 3;
  static const int NUM_NODES = ORDER*ORDER;

  MITC9( TacsScalar density, TacsScalar E, 
	 TacsScalar nu, TacsScalar t, TacsScalar _g[] );
  ~MITC9();

  // Code to test the residual and Jacobian implementations
  // ------------------------------------------------------
  void testResidual( double dh, 
		     const TacsScalar X[],
		     const TacsScalar vars[],
		     const TacsScalar dvars[],
		     const TacsScalar ddvars[] );
  void testJacobian( double dh, 
		     double alpha, double beta, double gamma,
		     const TacsScalar X[],
		     const TacsScalar vars[],
		     const TacsScalar dvars[],
		     const TacsScalar ddvars[] );  

  // Retrieve the initial values of the state variables
  // --------------------------------------------------
  void getInitVariables( TacsScalar q[] );

  // Compute the kinetic and potential energies
  // ------------------------------------------
  void computeEnergies( TacsScalar *_Te, TacsScalar *_Pe,
			const TacsScalar X[],
			const TacsScalar vars[],
			const TacsScalar dvars[] );

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void getResidual( TacsScalar res[],
		    const TacsScalar X[],
		    const TacsScalar vars[],
		    const TacsScalar dvars[],
		    const TacsScalar ddvars[] );

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void getJacobian( TacsScalar J[],
		    double alpha, double beta, double gamma,
		    const TacsScalar X[],
		    const TacsScalar vars[],
		    const TacsScalar dvars[],
		    const TacsScalar ddvars[] );

  // Test the strain implementation
  // ------------------------------
  void testStrain( const TacsScalar X[] );

  // Get the strain and the parametric location from the element
  // -----------------------------------------------------------
  TacsScalar getStrain( TacsScalar Xp[], TacsScalar e[],
			const double u, const double v,
			const TacsScalar X[],
			const TacsScalar vars[] );

 private:
  // Helper functions required for analysis
  void computeAngularVelocity( TacsScalar omega[],
			       const TacsScalar vars[],
			       const TacsScalar dvars[] );

  // Compute the angular acceleration at the nodes
  void computeAngularAccel( TacsScalar domega[],
			    const TacsScalar vars[],
			    const TacsScalar dvars[],
			    const TacsScalar ddvars[] );

  // Compute the reference frames at each node of the element
  void computeFrames( TacsScalar Xr[], const TacsScalar X[] );

  // Compute the directors for each node
  void computeDirectors( TacsScalar d[],
			 const TacsScalar vars[],
			 const TacsScalar Xr[] );

  // Compute the derivative of the directors w.r.t. the variables
  void computeDirectorDeriv( TacsScalar dirdq[],
			     const TacsScalar vars[],
			     const TacsScalar Xr[] );

  // Evaluate the strain
  void evalStrain( TacsScalar e[],
		   const TacsScalar Ur[], 
		   const TacsScalar dr[],
		   const TacsScalar Xdinv[],
		   const TacsScalar zXdinv[],
		   const TacsScalar T[] );

  // Evaluate the derivative of the strain w.r.t. the element variables 
  void evalBmat( TacsScalar e[], TacsScalar B[],
		 const double N[], const double Na[], const double Nb[],
		 const TacsScalar Ur[], const TacsScalar dr[],
		 const TacsScalar Xdinv[], const TacsScalar zXdinv[],
		 const TacsScalar T[], const TacsScalar dirdq[] );

  // Add the interpolated strain to the strain vector
  void addTyingStrain( TacsScalar e[],
		       const double N13[], const double N23[],
		       const TacsScalar g13[], const TacsScalar g23[],
		       const TacsScalar Xdinv[], const TacsScalar T[] );

  // Add the contribution from the tying strain to the b-matrix
  void addTyingBmat( TacsScalar B[],
		     const double N13[], const double N23[],
		     const TacsScalar b13[], const TacsScalar b23[],
		     const TacsScalar Xdinv[], const TacsScalar T[] );

  // Compute the shear strain at the tying points
  void computeTyingStrain( TacsScalar g13[], TacsScalar g23[],
			   const TacsScalar X[], const TacsScalar Xr[],
			   const TacsScalar vars[], 
			   const TacsScalar dir[] );

  // Compute the derivative of the strain at the tying points
  void computeTyingBmat( TacsScalar g13[], TacsScalar g23[],
			 TacsScalar b13[], TacsScalar b23[],
			 const TacsScalar X[], const TacsScalar Xr[],
			 const TacsScalar vars[], const TacsScalar dir[],
			 const TacsScalar dirdq[] );

  // Add the terms from the geometric stiffness matrix
  void addGmat( TacsScalar J[], 
		const TacsScalar scale, const TacsScalar s[],
		const double N[], const double Na[], const double Nb[],
		const TacsScalar Ur[], const TacsScalar dr[],
		const TacsScalar Xdinv[], const TacsScalar zXdinv[],
		const TacsScalar T[], const TacsScalar Xr[],
		const TacsScalar dirdq[] );

  // Add the term from the drilling rotation
  TacsScalar computeRotPenalty( const double N[],
				const TacsScalar Xa[],
				const TacsScalar Xb[],
				const TacsScalar Ua[],
				const TacsScalar Ub[],
				const TacsScalar vars[] );

  // Compute the derivative of the rotation term 
  TacsScalar computeBRotPenalty( TacsScalar brot[], const double N[], 
				 const double Na[], const double Nb[],
				 const TacsScalar Xa[], 
				 const TacsScalar Xb[],
				 const TacsScalar Ua[], 
				 const TacsScalar Ub[],
				 const TacsScalar vars[] );

  // Add the geometric stiffness term from the rotation
  void addGRotMat( TacsScalar J[], const TacsScalar scale, 
		   const double N[], const double Na[], const double Nb[],
		   const TacsScalar Xa[], const TacsScalar Xb[],
		   const TacsScalar Ua[], const TacsScalar Ub[],
		   const TacsScalar vars[] );

  // Add the geometric stiffness matrix from the tying strain
  void addTyingGmat( TacsScalar J[],
		     const TacsScalar w13[], const TacsScalar w23[],
		     const TacsScalar X[], const TacsScalar Xr[],
		     const TacsScalar vars[], const TacsScalar dir[],
		     const TacsScalar dirdq[] );

  // Add to the weights required to compute the geometric stiffness
  void addTyingGmatWeights( TacsScalar w13[], TacsScalar w23[],
			    const TacsScalar scalar,
			    const TacsScalar s[],
			    const double N13[], const double N23[],
			    const TacsScalar Xdinv[],
			    const TacsScalar T[] );

  // Compute the stress based on the strain values
  inline void computeStress( TacsScalar s[], const TacsScalar e[] ){
    s[0] = A[0]*e[0] + A[1]*e[1] + A[2]*e[2];
    s[1] = A[1]*e[0] + A[3]*e[1] + A[4]*e[2];
    s[2] = A[2]*e[0] + A[4]*e[1] + A[5]*e[2];
    
    s[3] = D[0]*e[3] + D[1]*e[4] + D[2]*e[5];
    s[4] = D[1]*e[3] + D[3]*e[4] + D[4]*e[5];
    s[5] = D[2]*e[3] + D[4]*e[4] + D[5]*e[5];

    s[6] = As[0]*e[6] + As[1]*e[7];
    s[7] = As[1]*e[6] + As[2]*e[7];
  }

  // Compute the product of the stress and the strain
  inline TacsScalar strainProduct( const TacsScalar s[], 
				   const TacsScalar e[] ){
    return (e[0]*s[0] + e[1]*s[1] + e[2]*s[2] +
	    e[3]*s[3] + e[4]*s[4] + e[5]*s[5] + 
	    e[6]*s[6] + e[7]*s[7]);
  }

  // Set the pointers to quadrature points/weights
  const double *gaussPts, *gaussWts;

  // The areal density/areal moment of inertia of the element
  TacsScalar rho, rhoI;

  // The stiffness matrices for the element
  TacsScalar A[6], D[6], As[3];
  TacsScalar kpenalty;

  // The gravity vector
  TacsScalar g[3];
};

/*
  Apply constraints between one or two nodes.

  This function applies constraints between one or two nodes and
  should be re-implemented as a constraint.
*/
class NodeConstraint : public TACSDynElement {
 public:
  static const double SCALE = 1000.0;

  NodeConstraint( int _num_nodes ){
    num_nodes = _num_nodes;
  }

  void getInitVariables( TacsScalar q[] ){
    memset(q, 0, 8*num_nodes*sizeof(TacsScalar));
    for ( int i = 1; i < num_nodes; i++ ){
      q[8*i+3] = 1.0;
    } 
  }

  // Compute the kinetic and potential energies
  // ------------------------------------------
  void computeEnergies( TacsScalar *_Te, TacsScalar *_Pe,
			const TacsScalar X[],
			const TacsScalar vars[],
			const TacsScalar dvars[] ){
    *_Te = 0.0;
    *_Pe = 0.0;
  }

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void getResidual( TacsScalar res[],
		    const TacsScalar X[],
		    const TacsScalar vars[],
		    const TacsScalar dvars[],
		    const TacsScalar ddvars[] ){
    memset(res, 0, 8*num_nodes*sizeof(TacsScalar));

    const TacsScalar *lamb = &vars[0];
    const TacsScalar *u1 = &vars[8]; 
    
    if (num_nodes == 2){
      res[0] = SCALE*u1[0];
      res[1] = SCALE*u1[1];
      res[2] = SCALE*u1[2];
      res[3] = SCALE*vars[3];
      res[4] = SCALE*vars[4];
      res[5] = SCALE*vars[5];
      res[6] = SCALE*vars[6];
      res[7] = SCALE*vars[7];
    
      res += 8;
      res[0] = SCALE*lamb[0];
      res[1] = SCALE*lamb[1];
      res[2] = SCALE*lamb[2];
      res[3] = res[4] = res[5] = res[6] = res[7] = 0.0;
    }
    else if (num_nodes == 3){
      const TacsScalar *u2 = &vars[16]; 
      res[0] = SCALE*(u1[0] - u2[0]);
      res[1] = SCALE*(u1[1] - u2[1]);
      res[2] = SCALE*(u1[2] - u2[2]);
      res[3] = SCALE*vars[3];
      res[4] = SCALE*vars[4];
      res[5] = SCALE*vars[5];
      res[6] = SCALE*vars[6];
      res[7] = SCALE*vars[7];
    
      res += 8;
      res[0] = SCALE*lamb[0];
      res[1] = SCALE*lamb[1];
      res[2] = SCALE*lamb[2];
      res[3] = res[4] = res[5] = res[6] = res[7] = 0.0;

      res += 8;
      res[0] = -SCALE*lamb[0];
      res[1] = -SCALE*lamb[1];
      res[2] = -SCALE*lamb[2];
      res[3] = res[4] = res[5] = res[6] = res[7] = 0.0;
    }
  }

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void getJacobian( TacsScalar J[],
		    double alpha, double beta, double gamma,
		    const TacsScalar X[],
		    const TacsScalar vars[],
		    const TacsScalar dvars[],
		    const TacsScalar ddvars[] ){
    memset(J, 0, 64*num_nodes*num_nodes*sizeof(TacsScalar));
    
    // Set the entries in the matrix
    for ( int k = 3; k < 8; k++ ){
      J[(8*num_nodes+1)*k] = SCALE*alpha;
    }

    for ( int k = 0; k < 3; k++ ){
      J[(8*num_nodes)*(8 + k) + k] = SCALE*alpha;
      J[(8*num_nodes)*k + 8 + k] = SCALE*alpha;
    }

    if (num_nodes == 3){
      for ( int k = 0; k < 3; k++ ){
	J[(8*num_nodes)*(16 + k) + k] = -SCALE*alpha;
	J[(8*num_nodes)*k + 16 + k] = -SCALE*alpha;
      }
    }
  }

  // The number of nodes (either 2 or 3)
  int num_nodes;
};

#endif // TACS_DYNAMICS_ELEMENT_H
