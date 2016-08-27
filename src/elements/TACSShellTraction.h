#ifndef TACS_SHELL_TRACTION_H
#define TACS_SHELL_TRACTION_H

/*
  Shell element traction class

  Copyright (c) 2010-2016 Graeme Kennedy. All rights reserved.
*/

#include "TACSElement.h"
#include "ShellUtils.h"

/*
  TACSShellTraction class

  This class defines a general shell traction. The traction may apply
  a force either normal to or in-plane of the shell. 
*/
template <int order>
class TACSShellTraction : public TACSElement {
 public:
  // Set the number of nodes
  static const int NUM_NODES = order*order;

  // Constructor for the shell element
  // ---------------------------------
  TACSShellTraction( TacsScalar _tx[], 
                     TacsScalar _ty[],
                     TacsScalar _tz[] ){
    memcpy(tx, _tx, NUM_NODES*sizeof(TacsScalar));
    memcpy(ty, _ty, NUM_NODES*sizeof(TacsScalar));
    memcpy(tz, _tz, NUM_NODES*sizeof(TacsScalar));
  }
  TACSShellTraction( TacsScalar _tx,
                     TacsScalar _ty,
                     TacsScalar _tz ){
    for ( int i = 0; i < NUM_NODES; i++ ){
      tx[i] = _tx;  ty[i] = _ty;  tz[i] = _tz;
    }
  }
                     

  // Return the number of displacements, stresses and nodes
  // ------------------------------------------------------
  int numDisplacements(){ return 6; }
  int numNodes(){ return NUM_NODES; }
  int numStresses(){ return 0; }
  
  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies( double time, 
                        TacsScalar *Te, TacsScalar *Pe,
                        const TacsScalar Xpts[],
                        const TacsScalar vars[],
                        const TacsScalar dvars[] ){
    *Te = 0.0, *Pe = 0.0;
  }

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){
    // Get the quadrature points and weights
    const double *gaussPts, *gaussWts;
    FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);

    // Add the residual due to the shell traction
    for ( int m = 0; m < order; m++ ){
      for ( int n = 0; n < order; n++ ){
        // Set the quadrature point
        double pt[2];
	pt[0] = gaussPts[n];
	pt[1] = gaussPts[m];	  
	
        // Compute X, Xd, N, Na and Nb
	double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
        TacsScalar X[3], Xd[9];
        shellutils::shell_jacobian(order, X, Xd, N, Na, Nb, pt, Xpts);

        // Determine the normal direction and normalize it
	Tensor::crossProduct3D(&Xd[6], &Xd[0], &Xd[3]);
	Tensor::normalize3D(&Xd[6]);
        	  
        // Compute the determinant of the Jacobian
	TacsScalar h = FElibrary::jacobian3d(Xd);
	h *= gaussWts[n]*gaussWts[m];

	// Evaluate the traction force evaluated at the
        // quadrature point within the element
        TacsScalar Tx = 0.0, Ty = 0.0, Tz = 0.0;
	for ( int i = 0; i < NUM_NODES; i++ ){
	  Tx += tx[i]*N[i];
	  Ty += ty[i]*N[i];
	  Tz += tz[i]*N[i];
	}
	
        // Add the contribution to the residual - the minus sign
        // is due to the fact that this is a work term
        TacsScalar *r = res;
	for ( int i = 0; i < NUM_NODES; i++ ){
	  r[0] -= h*Tx*N[i];
	  r[1] -= h*Ty*N[i];
	  r[2] -= h*Tz*N[i];
          r += 6;
	}
      }
    }
  }

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian( double time, TacsScalar J[],
                    double alpha, double beta, double gamma,
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){}

 private:
  TacsScalar tx[NUM_NODES], ty[NUM_NODES], tz[NUM_NODES];
};

/*
  TACSShellPressure class

  This class defines a surface pressure. The pressure is always normal
  to the surface.
*/
template <int order>
class TACSShellPressure : public TACSElement {
 public:
  // Set the number of nodes
  static const int NUM_NODES = order*order;

  // Constructor for the shell element
  // ---------------------------------
  TACSShellPressure( TacsScalar _p[] ){
    memcpy(p, _p, NUM_NODES*sizeof(TacsScalar));
  }
  TACSShellPressure( TacsScalar _p ){
    for ( int i = 0; i < NUM_NODES; i++ ){
      p[i] = _p;
    }
  }

  // Return the number of displacements, stresses and nodes
  // ------------------------------------------------------
  int numDisplacements(){ return 6; }
  int numNodes(){ return NUM_NODES; }
  int numStresses(){ return 0; }
  
  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies( double time, 
                        TacsScalar *Te, TacsScalar *Pe,
                        const TacsScalar Xpts[],
                        const TacsScalar vars[],
                        const TacsScalar dvars[] ){
    *Te = 0.0, *Pe = 0.0;
  }

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){
    // Get the quadrature points and weights
    const double *gaussPts, *gaussWts;
    FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);

    // Add the residual due to the shell traction
    for ( int m = 0; m < order; m++ ){
      for ( int n = 0; n < order; n++ ){
        // Set the quadrature point
        double pt[2];
	pt[0] = gaussPts[n];
	pt[1] = gaussPts[m];	  
	
        // Compute X, Xd, N, Na and Nb
	double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
        TacsScalar X[3], Xd[9];
        shellutils::shell_jacobian(order, X, Xd, N, Na, Nb, pt, Xpts);

        // Determine the normal direction and normalize it
	Tensor::crossProduct3D(&Xd[6], &Xd[0], &Xd[3]);
	Tensor::normalize3D(&Xd[6]);
        	  
        // Compute the determinant of the Jacobian
	TacsScalar h = FElibrary::jacobian3d(Xd);
	h *= gaussWts[n]*gaussWts[m];

	// Evaluate the traction force evaluated at the
        // quadrature point within the element
        TacsScalar P = 0.0;
	for ( int i = 0; i < NUM_NODES; i++ ){
	  P += p[i]*N[i];
	}
	
        // Add the contribution to the residual - the minus sign
        // is due to the fact that this is a work term
        TacsScalar *r = res;
	for ( int i = 0; i < NUM_NODES; i++ ){
	  r[0] -= h*P*Xd[6]*N[i];
	  r[1] -= h*P*Xd[7]*N[i];
	  r[2] -= h*P*Xd[8]*N[i];
          r += 6;
	}
      }
    }
  }

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian( double time, TacsScalar J[],
                    double alpha, double beta, double gamma,
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){}

 private:
  TacsScalar p[NUM_NODES];
};

#endif // TACS_SHELL_TRACTION
