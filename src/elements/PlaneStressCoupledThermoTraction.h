#ifndef PLANE_STRESS_COUPLED_TRACTION_H
#define PLANE_STRESS_COUPLED_TRACTION_H

#include "FElibrary.h"
/*
  The following file contains definitions for plane stress
  surface tractions for coupled thermoelastic elements
  
  Copyright (c) 2010-2016 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

/*
  The following class defines a traction for use with a quadrilateral
  element. 

  The surface argument is an integer that denotes the element surface
  over which the traction should be applied. The surfaces are indexed
  as follows in the xi/eta parameter space:

  ---- 3 ----
  |         |
  0         1
  |         |
  ---- 2 ----
  
  Note that the integration within the element occurs over each face.  
*/
template <int order>
class PSQuadThermoTraction : public TACSElement {
 public:
  PSQuadThermoTraction( int surf, 
                        TacsScalar _tx, TacsScalar _ty ){
    for ( int k = 0; k < order; k++ ){
      tx[k] = _tx;
      ty[k] = _ty;
    }
    initBaseDir(surf);
  }
  PSQuadThermoTraction( int surf, 
                        TacsScalar _tx[], TacsScalar _ty[] ){
    for ( int k = 0; k < order; k++ ){
      tx[k] = _tx[k];
      ty[k] = _ty[k];
    }
    initBaseDir(surf);
  }

  // Get the number of displacements/nodes
  int numDisplacements(){ return 3; }
  int numNodes(){ return order*order; }
  
  // Add the residual from the forces
  // --------------------------------
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){
    // Retrieve the quadrature scheme of the appropriate order
    const double *gaussPts, *gaussWts;
    int numGauss = FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);

    // Integrate over the specified element surface
    for ( int n = 0; n < numGauss; n++ ){
      double pt[2];
      pt[0] = gaussPts[n]*dir[0] + base[0]; 
      pt[1] = gaussPts[n]*dir[1] + base[1]; 

      // Evaluate the Lagrange basis in each direction
      double na[order], nb[order], dna[order], dnb[order];
      FElibrary::lagrangeSF(na, dna, pt[0], order);
      FElibrary::lagrangeSF(nb, dnb, pt[1], order);
    
      // Calcualte the Jacobian at the current point	
      const TacsScalar *x = Xpts;
      TacsScalar Xd[4] = {0.0, 0.0, 0.0, 0.0};
      for ( int j = 0; j < order; j++ ){
        for ( int i = 0; i < order; i++ ){
          Xd[0] += x[0]*dna[i]*nb[j];
          Xd[1] += x[0]*na[i]*dnb[j];
          
          Xd[2] += x[1]*dna[i]*nb[j];
          Xd[3] += x[1]*na[i]*dnb[j];
          x += 3;
        }
      }

      // Compute the derivative along each direction
      TacsScalar dx = Xd[0]*dir[0] + Xd[2]*dir[1];
      TacsScalar dy = Xd[1]*dir[0] + Xd[3]*dir[1];
      TacsScalar hsurf = gaussWts[n]*sqrt(dx*dx + dy*dy); 

      // Calculate the traction at the current point
      TacsScalar Tx = 0.0, Ty = 0.0;
      double N[order];
      FElibrary::lagrangeSF(N, gaussPts[n], order);
      for ( int i = 0; i < order; i++ ){
        Tx += N[i]*tx[i];
        Ty += N[i]*ty[i];
      }

      // Add the result to the element
      TacsScalar *r = res;
      for ( int j = 0; j < order; j++ ){
        for ( int i = 0; i < order; i++ ){
          r[0] -= hsurf*na[i]*nb[j]*Tx;
          r[1] -= hsurf*na[i]*nb[j]*Ty;
          r += 3;
        }
      }
    }
  }

 private:
  // Set the base point and direction
  // --------------------------------
  void initBaseDir( int surf ){
    // Determine the base point and integration direction
    if (surf == 0 || surf == 1){
      dir[0] = 0.0;
      dir[1] = 1.0;
    }
    else {
      dir[0] = 1.0;
      dir[1] = 0.0;
    }
    
    // Set the base point: The mid-point of the edge
    base[0] = base[1] = 0.0;
    if (surf == 0 || surf == 1){
      base[0] = -1.0 + 2.0*(surf % 2);
    }
    else {
      base[1] = -1.0 + 2.0*(surf % 2);
    }
  }

  // The traction information
  TacsScalar tx[order], ty[order];

  // The parametric base point and direction along which the
  // integration will occur.
  double base[2], dir[2];
};
/*
  Class for heat flux due to conduction, treat it like a traction
*/
template <int order>
class PSQuadHeatFluxTraction : public TACSElement {
 public:
  PSQuadHeatFluxTraction( int surf, 
                          TacsScalar _tx, 
                          TacsScalar _ty ){
    for ( int k = 0; k < order; k++ ){
      tx[k] = _tx;
      ty[k] = _ty;
    }
    initBaseDir(surf);
  }
  PSQuadHeatFluxTraction( int surf, 
                          TacsScalar _tx[], 
                          TacsScalar _ty[] ){
    for ( int k = 0; k < order; k++ ){
      tx[k] = _tx[k];
      ty[k] = _ty[k];
    }
    initBaseDir(surf);
  }
  // Get the number of displacements/nodes
  int numDisplacements(){ return 3; } // u,v,dT
  int numNodes(){ return order*order; }
  // Add the residual from the heat flux
  // ------------------------------------------------------
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){
    // Retrieve the quadrature scheme of the appropriate order
    const double *gaussPts, *gaussWts;
    int numGauss = FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);

    // Integrate over the specified element surface
    for ( int n = 0; n < numGauss; n++ ){
      double pt[2];
      pt[0] = gaussPts[n]*dir[0] + base[0]; 
      pt[1] = gaussPts[n]*dir[1] + base[1]; 

      // Evaluate the Lagrange basis in each direction
      double na[order], nb[order], dna[order], dnb[order];
      FElibrary::lagrangeSF(na, dna, pt[0], order);
      FElibrary::lagrangeSF(nb, dnb, pt[1], order);
    
      // Calcualte the Jacobian at the current point	
      const TacsScalar *x = Xpts;
      TacsScalar Xd[4] = {0.0, 0.0, 0.0, 0.0};
      for ( int j = 0; j < order; j++ ){
        for ( int i = 0; i < order; i++ ){
          Xd[0] += x[0]*dna[i]*nb[j];
          Xd[1] += x[0]*na[i]*dnb[j];
          
          Xd[2] += x[1]*dna[i]*nb[j];
          Xd[3] += x[1]*na[i]*dnb[j];
          x += 3;
        }
      }
      // ---------------------------------------------------
      // Compute the derivative along each direction
      TacsScalar dx = Xd[0]*dir[0] + Xd[2]*dir[1];
      TacsScalar dy = Xd[1]*dir[0] + Xd[3]*dir[1];
      TacsScalar hsurf = gaussWts[n]*sqrt(dx*dx + dy*dy);
      // --------------------------------------------------
      // Calculate the heat flux 
      TacsScalar qx = 0.0, qy = 0.0;
      double N[order];
      FElibrary::lagrangeSF(N, gaussPts[n], order);
      for ( int i = 0; i < order; i++ ){
        qx += N[i]*tx[i];
        qy += N[i]*ty[i];
      } // end for int i = 0; i < order
      // ---------------------------------------------------
      // Compute the heat flux outward normal to the edge
      TacsScalar nx = Xd[0]*norm_dir[0] + Xd[2]*norm_dir[1];
      TacsScalar ny = Xd[1]*norm_dir[0] + Xd[3]*norm_dir[1];
      // Normalize the vector
      nx /= sqrt(nx*nx + ny*ny);
      ny /= sqrt(nx*nx + ny*ny);
      //TacsScalar hsurf = gaussWts[n]*sqrt(nx*nx + ny*ny);
      // Compute flux that is outwards normal to the edge
      TacsScalar qn = qx*nx+qy*ny;
      
      // ---------------------------------------------(+ or -)
      // Add the result to the element
      TacsScalar *r = res;
      for ( int j = 0; j < order; j++ ){
        for ( int i = 0; i < order; i++ ){
          //printf("r: %e\n", r[2]);
          r[2] += hsurf*na[i]*nb[j]*qn;
          //printf("r: %e\n", hsurf*na[i]*nb[j]*qn);
          r += 3;
        }
      }
    } // end for int n = 0; n < numGauss
  }
  
 private:
  // Set the base point and direction
  // --------------------------------
  void initBaseDir( int surf ){
    // Determine the base point and integration direction
    if (surf == 0 || surf == 1){
      dir[0] = 0.0;
      dir[1] = 1.0;
      
      norm_dir[1] = 0.0;
      if (surf == 0){
        norm_dir[0] = -1.0;
      }
      else{
        norm_dir[1] = 1.0;
      }
    }
    else {
      dir[0] = 1.0;
      dir[1] = 0.0;
      
      norm_dir[0] = 0.0;
      if (surf == 2){
        norm_dir[1] = -1.0;
      }
      else{
        norm_dir[1] = 1.0;
      }
    }
    
    // Set the base point: The mid-point of the edge
    base[0] = base[1] = 0.0;
    if (surf == 0 || surf == 1){
      base[0] = -1.0 + 2.0*(surf % 2);
    }
    else {
      base[1] = -1.0 + 2.0*(surf % 2);
    }
  }

  // The heat flux information
  TacsScalar tx[order], ty[order];

  // The parametric base point and direction along which the
  // integration will occur as well as the normal direction of the
  // edge 
  double base[2], dir[2], norm_dir[2];
};
/*
  Class for heat source/sink in the plane quad element. Treat it like
  a body force
*/
template <int order>
class PSQuadHeatSourceSink : public TACSElement {
 public:
  PSQuadHeatSourceSink( TacsScalar _Q){
    for ( int i = 0; i < order*order; i++){
      Q[i] = _Q;
    }    
    numGauss = FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);
    
  }
  // Get the number of displacements/nodes
  int numDisplacements(){ return 3; } // u,v,dT
  int numNodes(){ return order*order; }
  void getShapeFunctions( const double pt[], double N[],
                          double Na[], double Nb[]){
    FElibrary::biLagrangeSF(N, Na, Nb, pt, order);
  }
  int getNumGaussPts(){ return numGauss*numGauss; }
  double getGaussWtsPts( int npoint, double pt[] ){
    // Compute the n/m/p indices of the Gauss quadrature scheme
    int m = (int)((npoint)/(numGauss));
    int n = npoint - numGauss*m;
    
    pt[0] = gaussPts[n];
    pt[1] = gaussPts[m];
    
    return gaussWts[n]*gaussWts[m];
  }
  // Add the residual from the heat source or sink through an area integral
  // ----------------------------------------------------------------------
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){
    // The shape functions associated with the element
    int NUM_NODES = order*order;
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES];
    
    // Get the number of quadrature points
    int numGauss = getNumGaussPts();
    
    for ( int n = 0; n < numGauss; n++ ){
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getGaussWtsPts(n, pt);
      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      planeJacobian(X, Xa, N, Na, Nb, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[4];
      TacsScalar h = FElibrary::jacobian2d(Xa, J);
      h = h*weight;
      // Add the contribution to the residual - the minus sign
      // is due to the fact that this is a work term
      for ( int i = 0; i < NUM_NODES; i++){
        res[3*i+2] -= h*Q[i]*N[i];
      }

    } // end for int n = 0; n < numGauss
  }
 private:
  // Information on the heat source/sink
  TacsScalar Q[order*order];
  // The Gauss quadrature scheme
  int numGauss;
  const double *gaussWts, *gaussPts;
};

#endif // PLANE_STRESS_COUPLED_TRACTION_H
