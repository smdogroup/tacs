#ifndef TACS_3D_COUPLED_THERMO_TRACTION_H
#define TACS_3D_COUPLED_THERMO_TRACTION_H

#include "TACSElement.h"
#include "TensorToolbox.h"
#include "FElibrary.h"

/*
  The surface traction class for 3D elements
*/
template <int order>
class TACS3DThermoTraction : public TACSElement {
 public:
  static const int NUM_NODES = order*order*order;

  static const int U_NEGATIVE = 0;
  static const int U_POSITIVE = 1;
  static const int V_NEGATIVE = 2;
  static const int V_POSITIVE = 3;
  static const int W_NEGATIVE = 4;
  static const int W_POSITIVE = 5;
  
  // Constructor for the shell element
  TACS3DThermoTraction( int _surface, 
                  TacsScalar _tx[], 
                  TacsScalar _ty[],
                  TacsScalar _tz[] ){
    surface = surface;
    memcpy(tx, _tx, order*order*sizeof(TacsScalar));
    memcpy(ty, _ty, order*order*sizeof(TacsScalar));
    memcpy(tz, _tz, order*order*sizeof(TacsScalar));

    // Set the knot locations
    if (order == 2){
      knots[0] = -1.0;
      knots[1] = 1.0;
    }
    else if (order == 3){
      knots[0] = -1.0;
      knots[1] = 0.0;
      knots[2] = 1.0;
    }
    else {
      // Set a co-sine spacing for the knot locations
      for ( int k = 0; k < order; k++ ){
        knots[k] = -cos(M_PI*k/(order-1));
      }
    }  
  
  }
  TACS3DThermoTraction( int _surface, 
                  TacsScalar _tx,
                  TacsScalar _ty,
                  TacsScalar _tz ){
    surface = _surface;
    for ( int i = 0; i < order*order; i++ ){
      tx[i] = _tx;  ty[i] = _ty;  tz[i] = _tz;
    }
    // Set the knot locations
    if (order == 2){
      knots[0] = -1.0;
      knots[1] = 1.0;
    }
    else if (order == 3){
      knots[0] = -1.0;
      knots[1] = 0.0;
      knots[2] = 1.0;
    }
    else {
      // Set a co-sine spacing for the knot locations
      for ( int k = 0; k < order; k++ ){
        knots[k] = -cos(M_PI*k/(order-1));
      }
    }  
  }
                     
  // Return the number of displacements, stresses and nodes
  // ------------------------------------------------------
  int numDisplacements(){ return 4; }
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

  /*
    Compute the shape functions and their derivatives w.r.t. the
    parametric element location 
  */
  void getShapeFunctions( const double pt[], double N[],
                          double Na[], double Nb[] ){
    double na[order], nb[order];
    double dna[order], dnb[order];
    FElibrary::lagrangeSFKnots(na, dna, pt[0], knots, order);
    FElibrary::lagrangeSFKnots(nb, dnb, pt[1], knots, order);
    for ( int j = 0; j < order; j++ ){
      for ( int i = 0; i < order; i++ ){
        N[0] = na[i]*nb[j];
        Na[0] = dna[i]*nb[j];
        Nb[0] = na[i]*dnb[j];
        N++;
        Na++;  Nb++;          
      }
    }
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
        getShapeFunctions(pt, N, Na, Nb);

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
        if (surface < 2){
          const int i = (order-1)*(surface % 2);
          for ( int k = 0; k < order; k++ ){
            for ( int j = 0; j < order; j++ ){
              const int node = i + j*order + k*order*order;
              res[4*node] -= h*Tx*N[j + k*order];
              res[4*node+1] -= h*Ty*N[j + k*order];
              res[4*node+2] -= h*Tz*N[j + k*order];
            }
          }
        }
        else if (surface < 4){
          const int j = (order-1)*(surface % 2);
          for ( int k = 0; k < order; k++ ){
            for ( int i = 0; i < order; i++ ){
              const int node = i + j*order + k*order*order;
              res[4*node] -= h*Tx*N[i + k*order];
              res[4*node+1] -= h*Ty*N[i + k*order];
              res[4*node+2] -= h*Tz*N[i + k*order];

            }
          }
        }
        else {
          const int k = (order-1)*(surface % 2);
          for ( int j = 0; j < order; j++ ){
            for ( int i = 0; i < order; i++ ){
              const int node = i + j*order + k*order*order;
              res[4*node] -= h*Tx*N[i + j*order];
              res[4*node+1] -= h*Ty*N[i + j*order];
              res[4*node+2] -= h*Tz*N[i + j*order];
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
  // The knot locations for the basis functions
  double knots[order];
};

/*
  The surface heat flux class for 3D elements
*/
template <int order>
class TACS3DHeatFluxTraction : public TACSElement {
 public:
  static const int NUM_NODES = order*order*order;

  static const int U_NEGATIVE = 0;
  static const int U_POSITIVE = 1;
  static const int V_NEGATIVE = 2;
  static const int V_POSITIVE = 3;
  static const int W_NEGATIVE = 4;
  static const int W_POSITIVE = 5;
  
  // Constructor for the shell element
  TACS3DHeatFluxTraction( int _surface, 
                          TacsScalar _tx[], 
                          TacsScalar _ty[],
                          TacsScalar _tz[] ){
    surface = surface;
    memcpy(tx, _tx, order*order*sizeof(TacsScalar));
    memcpy(ty, _ty, order*order*sizeof(TacsScalar));
    memcpy(tz, _tz, order*order*sizeof(TacsScalar));
    // Set the knot locations
    if (order == 2){
      knots[0] = -1.0;
      knots[1] = 1.0;
    }
    else if (order == 3){
      knots[0] = -1.0;
      knots[1] = 0.0;
      knots[2] = 1.0;
    }
    else {
      // Set a co-sine spacing for the knot locations
      for ( int k = 0; k < order; k++ ){
        knots[k] = -cos(M_PI*k/(order-1));
      }
    }
  }
  TACS3DHeatFluxTraction( int _surface, 
                          TacsScalar _tx,
                          TacsScalar _ty,
                          TacsScalar _tz ){
    surface = _surface;
    for ( int i = 0; i < order*order; i++ ){
      tx[i] = _tx;  ty[i] = _ty;  tz[i] = _tz;
    }
    // Set the knot locations
    if (order == 2){
      knots[0] = -1.0;
      knots[1] = 1.0;
    }
    else if (order == 3){
      knots[0] = -1.0;
      knots[1] = 0.0;
      knots[2] = 1.0;
    }
    else {
      // Set a co-sine spacing for the knot locations
      for ( int k = 0; k < order; k++ ){
        knots[k] = -cos(M_PI*k/(order-1));
      }
    }
  }
                     
  // Return the number of displacements, stresses and nodes
  // ------------------------------------------------------
  int numDisplacements(){ return 4; }
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
  /*
    Compute the shape functions and their derivatives w.r.t. the
    parametric element location 
  */
  void getShapeFunctions( const double pt[], double N[],
                          double Na[], double Nb[] ){
    double na[order], nb[order];
    double dna[order], dnb[order];
    FElibrary::lagrangeSFKnots(na, dna, pt[0], knots, order);
    FElibrary::lagrangeSFKnots(nb, dnb, pt[1], knots, order);
    for ( int j = 0; j < order; j++ ){
      for ( int i = 0; i < order; i++ ){
        N[0] = na[i]*nb[j];
        Na[0] = dna[i]*nb[j];
        Nb[0] = na[i]*dnb[j];
        N++;
        Na++;  Nb++;          
      }
    }
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
        getShapeFunctions(pt, N, Na, Nb);

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

        // Evaluate the heat flux force evaluated at the
        // quadrature point within the element
        TacsScalar Qx = 0.0, Qy = 0.0, Qz = 0.0;
        for ( int i = 0; i < order*order; i++ ){
          Qx += tx[i]*N[i];
          Qy += ty[i]*N[i];
          Qz += tz[i]*N[i];
        }

        // Compute the heat flux outward normal to the edge
        TacsScalar n_dir[3];
        Tensor::normalize3D(n_dir,normal);
        TacsScalar qn = Qx*n_dir[0]+Qy*n_dir[1]+Qz*n_dir[2];

        // Add the contribution to the residual - the minus sign
        // is due to the fact that this is a work term
        if (surface < 2){
          const int i = (order-1)*(surface % 2);
          for ( int k = 0; k < order; k++ ){
            for ( int j = 0; j < order; j++ ){
              const int node = i + j*order + k*order*order;
              res[4*node+3] += h*qn*N[j + k*order];
            }
          }
        }
        else if (surface < 4){
          const int j = (order-1)*(surface % 2);
          for ( int k = 0; k < order; k++ ){
            for ( int i = 0; i < order; i++ ){
              const int node = i + j*order + k*order*order;
              res[4*node+3] += h*qn*N[j + k*order];
            }
          }
        }
        else {
          const int k = (order-1)*(surface % 2);
          for ( int j = 0; j < order; j++ ){
            for ( int i = 0; i < order; i++ ){
              const int node = i + j*order + k*order*order;
              res[4*node+3] += h*qn*N[j + k*order];
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
  // The knot locations for the basis functions
  double knots[order];
};

/*
  The surface heat source/sink class for 3D elements; treat it like a
  body force
*/
template <int order>
class TACS3DHeatSourceSink : public TACSElement {
 public:
  static const int NUM_NODES = order*order*order;

  static const int U_NEGATIVE = 0;
  static const int U_POSITIVE = 1;
  static const int V_NEGATIVE = 2;
  static const int V_POSITIVE = 3;
  static const int W_NEGATIVE = 4;
  static const int W_POSITIVE = 5;

  TACS3DHeatSourceSink( TacsScalar _Q){
    for ( int i = 0; i < NUM_NODES; i++){
      Q[i] = _Q;
    }
    // Set the knot locations
    if (order == 2){
      knots[0] = -1.0;
      knots[1] = 1.0;
    }
    else if (order == 3){
      knots[0] = -1.0;
      knots[1] = 0.0;
      knots[2] = 1.0;
    }
    else {
      // Set a co-sine spacing for the knot locations
      for ( int k = 0; k < order; k++ ){
        knots[k] = -cos(M_PI*k/(order-1));
      }
    }
    
    numGauss = FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);    
  }
  // Get the number of displacements/nodes
  int numDisplacements(){ return 4; } // u,v,w,dT
  void getShapeFunctions( const double pt[], double N[],
                          double Na[], double Nb[], double Nc[]){
    double na[order], nb[order], nc[order];
    double dna[order], dnb[order], dnc[order];
    FElibrary::lagrangeSFKnots(na, dna, pt[0], knots, order);
    FElibrary::lagrangeSFKnots(nb, dnb, pt[1], knots, order);
    FElibrary::lagrangeSFKnots(nc, dnc, pt[2], knots, order);
    for ( int k = 0; k < order; k++ ){
      for ( int j = 0; j < order; j++ ){
        for ( int i = 0; i < order; i++ ){
          N[0] = na[i]*nb[j]*nc[k];
          Na[0] = dna[i]*nb[j]*nc[k];
          Nb[0] = na[i]*dnb[j]*nc[k];
          Nc[0] = na[i]*nb[j]*dnc[k];
          N++;
          Na++;  Nb++;  Nc++;          
        }
      }
    }
  }
  int getNumGaussPts(){ return numGauss*numGauss*numGauss; }
  double getGaussWtsPts( int npoint, double pt[] ){
    // Compute the n/m/p indices of the Gauss quadrature scheme
    int p = (int)((npoint)/(numGauss*numGauss));
    int m = (int)((npoint - numGauss*numGauss*p)/numGauss);
    int n = npoint - numGauss*m - numGauss*numGauss*p;
  
    pt[0] = gaussPts[n];
    pt[1] = gaussPts[m];    
    pt[2] = gaussPts[p];
  
    return gaussWts[n]*gaussWts[m]*gaussWts[p];    
  }
  // Add the residual from the heat source or sink through a volume integral
  // ----------------------------------------------------------------------
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
    
    // Get the number of quadrature points
    int numGauss = getNumGaussPts();
    for ( int n = 0; n < numGauss; n++ ){
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getGaussWtsPts(n, pt);
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
      
      // Add the contribution to the residual - the minus sign
      // is due to the fact that this is a work term
      for ( int i = 0; i < NUM_NODES; i++){
        res[4*i+3] -= h*Q[i]*N[i];
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
  // Information on the heat source/sink
  TacsScalar Q[NUM_NODES];
  // The Gauss quadrature scheme
  int numGauss;
  const double *gaussWts, *gaussPts;
  // The knot locations for the basis functions
  double knots[order];
};

#endif // TACS_3D_COUPLED_THERMO_TRACTION_H
