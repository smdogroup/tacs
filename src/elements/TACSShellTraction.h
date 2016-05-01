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
    memcpy(tx, _tx, order*order*sizeof(TacsScalar));
    memcpy(ty, _ty, order*order*sizeof(TacsScalar));
    memcpy(tz, _tz, order*order*sizeof(TacsScalar));
  }

  // Return the number of displacements, stresses and nodes
  // ------------------------------------------------------
  int numDisplacements(){ return 6; }
  int numNodes(){ return order*order; }
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
    FElibrary::int getGaussPtsWts(order, &gaussPts, &gaussWts);

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
        	  
        // Compute the determinant of the Jacobian
	TacsScalar h = FElibrary::jacobian3d(Xd);
	h *= gaussWts[n]*gaussWts[m];

	// Evaluate the traction force evaluated at the
        // quadrature point within the element
        TacsScalar Tx = 0.0, Ty = 0.0, Tz = 0.0;
	for ( int i = 0; i < order*order; i++ ){
	  Tx += tx[i]*N[i];
	  Ty += ty[i]*N[i];
	  Tz += tz[i]*N[i];
	}
	
        // Add the contribution to the residual - the minus sign
        // is due to the fact that this is a work term
        TacsScalar *r = res;
	for ( int i = 0; i < order*order; i++ ){
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

#endif // TACS_SHELL_TRACTION
