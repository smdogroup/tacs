/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TACSElementBasis.h"
#include "TACSElementAlgebra.h"

/*
  Get the field values at the specified quadrature point
*/
void TACSElementBasis::getFieldValues( const double pt[],
                                       const TacsScalar Xpts[],
                                       const int vars_per_node,
                                       const TacsScalar vars[],
                                       TacsScalar X[],
                                       TacsScalar U[] ){
  const int num_nodes = getNumNodes();

  // Compute the basis functions
  double N[MAX_BASIS_SIZE];
  computeBasis(pt, N);

  computeFieldValues(num_nodes, N, Xpts, vars_per_node, X, U);
}

void TACSElementBasis::computeFieldValues( const int num_nodes,
                                           const double N[],
                                           const TacsScalar Xpts[],
                                           const int vars_per_node,
                                           const TacsScalar vars[],
                                           TacsScalar X[],
                                           TacsScalar U[] ){
  // Zero the values of the coordinates and state variable values
  X[0] = X[1] = X[2] = 0.0;
  for ( int j = 0; j < vars_per_node; j++ ){
    U[j] = 0.0;
  }

  // Loop over each quadrature point for each basis function
  const double *n = N;
  for ( int i = 0; i < num_nodes; i++ ){
    // Add the contribution to the coordinates
    X[0] += n[0]*Xpts[0];
    X[1] += n[0]*Xpts[1];
    X[2] += n[0]*Xpts[2];
    n++;
    Xpts += 3;
  }

  // Re-set the shape function pointer
  n = N;
  for ( int i = 0; i < num_nodes; i++ ){
    // Add the contribution to the state variable values
    TacsScalar *u = U;
    for ( int j = 0; j < vars_per_node; j++ ){
      u[0] += n[0]*vars[0];
      vars++;
      u++;
    }

    // Increment the pointer to the shape functions
    n++
  }
}

/*
  Get the field values at the specified quadrature point
*/
void TACSElementBasis::getFieldValues( int n,
                                       const TacsScalar Xpts[],
                                       const int vars_per_node,
                                       const TacsScalar vars[],
                                       TacsScalar X[],
                                       TacsScalar U[] ){
  double pt[3];
  getQuadraturePoint(n, pt);
  getFieldValues(pt, Xpts, vars_per_node, X, U);
}

/*
  Get the gradient of the field at the quadrature point
*/
TacsScalar TACSElementBasis::getFieldGradient( const double pt[],
                                               const TacsScalar Xpts[],
                                               const int vars_per_node,
                                               const TacsScalar vars[],
                                               const TacsScalar dvars[],
                                               const TacsScalar ddvars[],
                                               TacsScalar X[],
                                               TacsScalar Xd[],
                                               TacsScalar J[],
                                               TacsScalar U[],
                                               TacsScalar Udot[],
                                               TacsScalar Uddot[],
                                               TacsScalar Ux[] ){
  const int num_nodes = getNumNodes();
  const int num_params = getNumParameters();
  double N[MAX_BASIS_SIZE], Nxi[3*MAX_BASIS_SIZE];
  computeBasisGradient(pt, N, Nxi);

}

TacsScalar TACSElementBasis::computeFieldGradient( const int num_params,
                                                   const int num_nodes,
                                                   const double N[],
                                                   const double Nxi[],
                                                   const TacsScalar Xpts[],
                                                   const int vars_per_node,
                                                   const TacsScalar vars[],
                                                   const TacsScalar dvars[],
                                                   const TacsScalar ddvars[],
                                                   TacsScalar X[],
                                                   TacsScalar Xd[],
                                                   TacsScalar J[],
                                                   TacsScalar U[],
                                                   TacsScalar Udot[],
                                                   TacsScalar Uddot[],
                                                   TacsScalar Ux[] ){
  if (num_params == 3){
    // Zero the values of the coordinate and its derivative
    X[0] = X[1] = X[2] = 0.0;
    Xd[0] = Xd[1] = Xd[2] = 0.0;
    Xd[3] = Xd[4] = Xd[5] = 0.0;
    Xd[6] = Xd[7] = Xd[8] = 0.0;

    // Zero the values of the displacements and their derivatives
    for ( int j = 0; j < vars_per_node; j++ ){
      U[j] = 0.0;
      Ux[3*j] = Ux[3*j+1] = Ux[3*j+2] = 0.0;
    }

    // Loop over each quadrature point for each basis function
    const double *n = N, *nxi = Nxi;
    for ( int i = 0; i < num_nodes; i++ ){
      X[0] += n[0]*Xpts[0];
      X[1] += n[0]*Xpts[1];
      X[2] += n[0]*Xpts[2];

      Xd[0] += nxi[0]*Xpts[0];
      Xd[1] += nxi[1]*Xpts[0];
      Xd[2] += nxi[2]*Xpts[0];

      Xd[3] += nxi[0]*Xpts[1];
      Xd[4] += nxi[1]*Xpts[1];
      Xd[5] += nxi[2]*Xpts[1];

      Xd[6] += nxi[0]*Xpts[2];
      Xd[7] += nxi[1]*Xpts[2];
      Xd[8] += nxi[2]*Xpts[2];
      Xpts += 3;
      n++;
      nxi += 3;
    }

    n = N, nxi = Nxi;
    for ( int i = 0; i < num_nodes; i++ ){
      // Add contributions to the derivatives of the displacements
      TacsScalar *u = U, *ux = Ux;
      for ( int j = 0; j < vars_per_node; j++ ){
        u[0] += n[0]*vars[0];

        ux[0] += nxi[0]*vars[0];
        ux[1] += nxi[1]*vars[0];
        ux[2] += nxi[2]*vars[0];

        vars++;
        u++;
        nx += 3;
      }

      n++;
      nxi += 3;
    }
  }
  else if (param_size == 2){
    // Zero the values of the coordinate and its derivative
    X[0] = X[1] = X[2] = 0.0;
    Xd[0] = Xd[1] = Xd[2] = Xd[3] = 0.0;

    // Zero the values of the displacements and their derivatives
    for ( int j = 0; j < vars_per_node; j++ ){
      U[j] = 0.0;
      Ux[2*j] = Ux[2*j+1] = 0.0;
    }

    // Loop over each quadrature point for each basis function
    const double *n = N, *nxi = Nxi;
    for ( int i = 0; i < num_nodes; i++ ){
      X[0] += n[0]*Xpts[0];
      X[1] += n[0]*Xpts[1];
      X[2] += n[0]*Xpts[2];

      X[0] += nxi[0]*Xpts[0];
      X[1] += nxi[1]*Xpts[0];

      X[2] += nxi[0]*Xpts[1];
      X[3] += nxi[1]*Xpts[1];
      Xpts += 3;

      // Add contributions to the derivatives of the displacements
      TacsScalar *u = U, *ux = Ux;
      for ( int j = 0; j < vars_per_node; j++ ){
        u[0] += n[0]*vars[0];

        ux[0] += nxi[0]*vars[0];
        ux[1] += nxi[1]*vars[0];

        vars++;
        u++;
        nx += 2;
      }

      n++;
      nxi += 2;
    }
  }
  else if (param_size == 1){
    // Zero the values of the coordinate and its derivative
    X[0] = X[1] = X[2] = 0.0;
    Xd[0] = Xd[1] = Xd[2] = 0.0;

    // Zero the values of the displacements and their derivatives
    for ( int j = 0; j < vars_per_node; j++ ){
      U[j] = 0.0;
      Ux[2*j] = 0.0;
    }

    // Loop over each quadrature point for each basis function
    const double *n = N, *nxi = Nxi;
    for ( int i = 0; i < num_nodes; i++ ){
      X[0] += n[0]*Xpts[0];
      X[1] += n[0]*Xpts[1];
      X[2] += n[0]*Xpts[2];

      Xd[0] += nxi[0]*Xpts[0];
      Xd[1] += nxi[1]*Xpts[0];

      Xd[2] += nxi[0]*Xpts[1];
      Xd[3] += nxi[1]*Xpts[1];

      Xd[4] += nxi[0]*Xpts[1];
      Xd[5] += nxi[1]*Xpts[1];
      Xpts += 3;

      // Add contributions to the derivatives of the displacements
      TacsScalar *u = U, *ux = Ux;
      for ( int j = 0; j < vars_per_node; j++ ){
        u[0] += n[0]*vars[0];

        ux[0] += nxi[0]*vars[0];
        ux[1] += nxi[1]*vars[0];

        vars++;
        u++;
        nxi;
      }

      n++;
      nxi += 1;
    }

  }
}

/*
  Add the weak form of the governing equations to the residual
*/
void TACSElementBasis::addWeakFormResidual( int n,
                                            const double pt[],
                                            TacsScalar weight,
                                            const TacsScalar J[],
                                            const int vars_per_node,
                                            const TacsScalar DUt[],
                                            const TacsScalar DUx[],
                                            TacsScalar res[] ){
  const int num_nodes = getNumNodes();
  double N[MAX_BASIS_SIZE], Nxi[3*MAX_BASIS_SIZE];
  TacsScalar Nx[3*MAX_BASIS_SIZE];

  // Compute the basis and its derivative
  computeBasisGradient(pt, N, Nxi);


}

void

  // Transform Nxi to Nx
  for ( int i = 0; i < num_nodes; i++ ){
    matMult(J, &Nxi[3*i], &Nx[3*i]);
  }

  // Add contributions from DUt
  const double *n = N;
  TacsScalar *r = res;
  for ( int i = 0; i < num_nodes; i++ ){
    for ( int ii = 0; ii < vars_per_node; ii++ ){
      r[0] += weight*n[0]*(DUt[3*ii] + DUt[3*ii+1] + DUt[3*ii+2]);
      r++;
    }
    n++;
  }

  // Add contributions from DUx
  n = N;
  const double *nx = Nx;
  r = res;

  if (getParameterSize() == 3){
    for ( int i = 0; i < num_nodes; i++ ){
      for ( int ii = 0; ii < vars_per_node; ii++ ){
        r[0] += weight*(n[0]*DUx[4*ii] + nx[0]*DUx[4*ii+1] +
                        nx[1]*DUx[4*ii+2] + nx[2]*DUx[4*ii+3]);
        r++;
      }
      n++;
      nx += 3;
    }
  }
  else if (getParameterSize() == 2){
    for ( int i = 0; i < num_nodes; i++ ){
      for ( int ii = 0; ii < vars_per_node; ii++ ){
        r[0] += weight*(n[0]*DUx[3*ii] + nx[0]*DUx[3*ii+1] + nx[1]*DUx[3*ii+2]);
        r++;
      }
      n++;
      nx += 2;
    }
  }
  else if (getParameterSize() == 1){
    for ( int i = 0; i < num_nodes; i++ ){
      for ( int ii = 0; ii < vars_per_node; ii++ ){
        r[0] += weight*(n[0]*DUx[2*ii] + nx[0]*DUx[2*ii+1]);
        r++;
      }
      n++;
      nx++;
    }
  }
}

/*
  Add the weak form of the governing equations to the residual
*/
void TACSElementBasis::addWeakFormJacobian( int n,
                                            const double pt[],
                                            TacsScalar weight,
                                            const TacsScalar J[],
                                            const int vars_per_node,
                                            const TacsScalar DUt[],
                                            const TacsScalar DUx[],
                                            double alpha, double beta, double gamma,
                                            const int DDUt_nnz,
                                            const int *DDUt_pairs,
                                            const TacsScalar *DDUt,
                                            const int DDUx_nnz,
                                            const int *DDUx_paris,
                                            const TacsScalar *DDUx,
                                            TacsScalar *res,
                                            TacsScalar *mat ){
  const int num_nodes = getNumNodes();
  double N[MAX_BASIS_SIZE], Nxi[3*MAX_BASIS_SIZE];
  TacsScalar Nx[3*MAX_BASIS_SIZE];

  // Compute the basis and its derivative
  computeBasisGradient(pt, N, Nxi);

  // Transform Nxi to Nx
  for ( int i = 0; i < num_nodes; i++ ){
    matMult(J, &Nxi[3*i], &Nx[3*i]);
  }

  // Transform Nxi to Nx
  for ( int i = 0; i < num_nodes; i++ ){
    matMult(J, &Nxi[3*i], &Nx[3*i]);
  }

  // Add contributions from DUt
  const double *n = N;
  TacsScalar *r = res;
  for ( int i = 0; i < num_nodes; i++ ){
    for ( int ii = 0; ii < vars_per_node; ii++ ){
      r[0] += weight*n[0]*(DUt[3*ii] + DUt[3*ii+1] + DUt[3*ii+2]);
      r++;
    }
    n++;
  }

    // Decode whether to use the dense or sparse representation
  if (DDUt_nnz < 0){
    for ( int ix = 0; ix < 3*vars_per_node; ix++ ){
      for ( int jx = 0; jx < 3*vars_per_node; jx++ ){
        if (DDUt[ix*3*vars_per_node + jx] != 0.0){
          TacsScalar scale = weight*DDUt[ix*3*vars_per_node + jx];
          if (jx % 3 == 0){ scale *= alpha; }
          else if (jx % 3 == 1){ scale *= beta; }
          else { scale *= gamma; }

          for ( int i = 0; i < num_nodes; i++ ){
            for ( int j = 0; j < vars_per_node; j++ ){
              mat[ vars_per_node*i
            }
          }
        }
      }
    }
  }
  else {
    // Use a sparse representation
    

    
  }  

  // Add contributions from DUx
  n = N;
  const double *nx = Nx;
  r = res;

  if (getParameterSize() == 3){
    for ( int i = 0; i < num_nodes; i++ ){
      for ( int ii = 0; ii < vars_per_node; ii++ ){
        r[0] += weight*(n[0]*DUx[4*ii] + nx[0]*DUx[4*ii+1] +
                        nx[1]*DUx[4*ii+2] + nx[2]*DUx[4*ii+3]);
        r++;
      }
      n++;
      nx += 3;
    }
  }
  else if (getParameterSize() == 2){
    for ( int i = 0; i < num_nodes; i++ ){
      for ( int ii = 0; ii < vars_per_node; ii++ ){
        r[0] += weight*(n[0]*DUx[3*ii] + nx[0]*DUx[3*ii+1] + nx[1]*DUx[3*ii+2]);
        r++;
      }
      n++;
      nx += 2;
    }
  }
  else if (getParameterSize() == 1){
    for ( int i = 0; i < num_nodes; i++ ){
      for ( int ii = 0; ii < vars_per_node; ii++ ){
        r[0] += weight*(n[0]*DUx[2*ii] + nx[0]*DUx[2*ii+1]);
        r++;
      }
      n++;
      nx++;
    }
  }
}

