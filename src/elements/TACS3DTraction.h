#ifndef TACS_3D_TRACTION_H
#define TACS_3D_TRACTION_H

#include "TACSElement.h"
#include "TensorToolbox.h"
#include "FElibrary.h"

/*
  The surface traction class for 3D elements
*/
template <int order>
class TACS3DTraction : public TACSElement {
 public:
  static const int NUM_NODES = order*order*order;

  static const int U_NEGATIVE = 0;
  static const int U_POSITIVE = 1;
  static const int V_NEGATIVE = 2;
  static const int V_POSITIVE = 3;
  static const int W_NEGATIVE = 4;
  static const int W_POSITIVE = 5;
  
  // Constructor for the shell element
  TACS3DTraction( int _surface, 
                  TacsScalar _tx[], 
                  TacsScalar _ty[],
                  TacsScalar _tz[] ){
    surface = surface;
    memcpy(tx, _tx, order*order*sizeof(TacsScalar));
    memcpy(ty, _ty, order*order*sizeof(TacsScalar));
    memcpy(tz, _tz, order*order*sizeof(TacsScalar));
  }
  TACS3DTraction( int _surface, 
                  TacsScalar _tx,
                  TacsScalar _ty,
                  TacsScalar _tz ){
    surface = _surface;
    for ( int i = 0; i < order*order; i++ ){
      tx[i] = _tx;  ty[i] = _ty;  tz[i] = _tz;
    }
  }
                     
  // Return the number of displacements, stresses and nodes
  // ------------------------------------------------------
  int numDisplacements(){ return 3; }
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
        double pt[3];
        pt[0] = gaussPts[n];
        pt[1] = gaussPts[m];

        // Compute X, Xd, N, Na and Nb
	double N[order*order];
        double Na[order*order], Nb[order*order];
        FElibrary::biLagrangeSF(N, Na, Nb, pt, order);

        TacsScalar Xa[3], Xb[3];
        Xa[0] = Xa[1] = Xa[2] = 0.0;
        Xb[0] = Xb[1] = Xb[2] = 0.0;

        if (surface < 2){
          const int i = (order-1)*(surface % 2);
          for ( int k = 0; k < order; k++ ){
            for ( int j = 0; j < order; j++ ){
              const int node = i + j*order + k*order*order;
              
              Xa[0] += Na[j + k*order]*Xpts[3*node];
              Xa[1] += Na[j + k*order]*Xpts[3*node+1];
              Xa[2] += Na[j + k*order]*Xpts[3*node+2];

              Xb[0] += Nb[j + k*order]*Xpts[3*node];
              Xb[1] += Nb[j + k*order]*Xpts[3*node+1];
              Xb[2] += Nb[j + k*order]*Xpts[3*node+2];
            }
          }
        }
        else if (surface < 4){
          const int j = (order-1)*(surface % 2);
          for ( int k = 0; k < order; k++ ){
            for ( int i = 0; i < order; i++ ){
              const int node = i + j*order + k*order*order;
              
              Xa[0] += Na[i + k*order]*Xpts[3*node];
              Xa[1] += Na[i + k*order]*Xpts[3*node+1];
              Xa[2] += Na[i + k*order]*Xpts[3*node+2];

              Xb[0] += Nb[i + k*order]*Xpts[3*node];
              Xb[1] += Nb[i + k*order]*Xpts[3*node+1];
              Xb[2] += Nb[i + k*order]*Xpts[3*node+2];
            }
          }
        }
        else {
          const int k = (order-1)*(surface % 2);
          for ( int j = 0; j < order; j++ ){
            for ( int i = 0; i < order; i++ ){
              const int node = i + j*order + k*order*order;
              
              Xa[0] += Na[i + j*order]*Xpts[3*node];
              Xa[1] += Na[i + j*order]*Xpts[3*node+1];
              Xa[2] += Na[i + j*order]*Xpts[3*node+2];

              Xb[0] += Nb[i + j*order]*Xpts[3*node];
              Xb[1] += Nb[i + j*order]*Xpts[3*node+1];
              Xb[2] += Nb[i + j*order]*Xpts[3*node+2];
            }
          }
        }

        // Determine the normal direction
        TacsScalar normal[3];
	Tensor::crossProduct3D(normal, Xa, Xb);
        TacsScalar h = sqrt(normal[0]*normal[0] +
                            normal[1]*normal[1] +
                            normal[2]*normal[2]);
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
        if (surface < 2){
          const int i = (order-1)*(surface % 2);
          for ( int k = 0; k < order; k++ ){
            for ( int j = 0; j < order; j++ ){
              const int node = i + j*order + k*order*order;
              res[3*node] -= h*Tx*N[j + k*order];
              res[3*node+1] -= h*Ty*N[j + k*order];
              res[3*node+2] -= h*Tz*N[j + k*order];
            }
          }
        }
        else if (surface < 4){
          const int j = (order-1)*(surface % 2);
          for ( int k = 0; k < order; k++ ){
            for ( int i = 0; i < order; i++ ){
              const int node = i + j*order + k*order*order;
              res[3*node] -= h*Tx*N[i + k*order];
              res[3*node+1] -= h*Ty*N[i + k*order];
              res[3*node+2] -= h*Tz*N[i + k*order];

            }
          }
        }
        else {
          const int k = (order-1)*(surface % 2);
          for ( int j = 0; j < order; j++ ){
            for ( int i = 0; i < order; i++ ){
              const int node = i + j*order + k*order*order;
              res[3*node] -= h*Tx*N[i + j*order];
              res[3*node+1] -= h*Ty*N[i + j*order];
              res[3*node+2] -= h*Tz*N[i + j*order];
            }
          }
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
  int surface;
  TacsScalar tx[order*order];
  TacsScalar ty[order*order];
  TacsScalar tz[order*order];
};


#endif // TACS_3D_TRACTION_H
