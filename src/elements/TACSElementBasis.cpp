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

  computeFieldValues(num_nodes, N, Xpts, vars_per_node, vars, X, U);
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
    n++;
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
  getFieldValues(pt, Xpts, vars_per_node, vars, X, U);
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
  const int num_params = getNumParameters();
  const int num_nodes = getNumNodes();
  double N[MAX_BASIS_SIZE], Nxi[3*MAX_BASIS_SIZE];
  computeBasisGradient(pt, N, Nxi);

  return computeFieldGradient(num_params, num_nodes, N, Nxi, Xpts,
                              vars_per_node, vars, dvars, ddvars,
                              X, Xd, J, U, Udot, Uddot, Ux);
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
      Udot[j] = 0.0;
      Uddot[j] = 0.0;
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
      TacsScalar *u = U, *udot = Udot, *uddot = Uddot, *ux = Ux;
      for ( int j = 0; j < vars_per_node; j++ ){
        u[0] += n[0]*vars[0];
        udot[0] += n[0]*dvars[0];
        uddot[0] += n[0]*ddvars[0];

        ux[0] += nxi[0]*vars[0];
        ux[1] += nxi[1]*vars[0];
        ux[2] += nxi[2]*vars[0];

        vars++;
        dvars++;
        ddvars++;
        u++;
        udot++;
        uddot++;
        ux += 3;
      }

      n++;
      nxi += 3;
    }

    // Compute the Jacobian transformation
    TacsScalar detJ = inv3x3(Xd, J);

    // U,x = U,xi * J
    TacsScalar *ux = Ux;
    for ( int j = 0; j < vars_per_node; j++ ){
      TacsScalar ud[3];
      ud[0] = ux[0];
      ud[1] = ux[1];
      ud[2] = ux[2];

      mat3x3MultTrans(J, ud, ux);

      ux += 3;
    }

    return detJ;
  }
  else if (num_params == 2){
    // Zero the values of the coordinate and its derivative
    X[0] = X[1] = X[2] = 0.0;
    Xd[0] = Xd[1] = Xd[2] = Xd[3] = 0.0;

    // Zero the values of the displacements and their derivatives
    for ( int j = 0; j < vars_per_node; j++ ){
      U[j] = 0.0;
      Udot[j] = 0.0;
      Uddot[j] = 0.0;
      Ux[2*j] = Ux[2*j+1] = 0.0;
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
      Xpts += 3;

      // Add contributions to the derivatives of the displacements
      TacsScalar *u = U, *udot = Udot, *uddot = Uddot, *ux = Ux;
      for ( int j = 0; j < vars_per_node; j++ ){
        u[0] += n[0]*vars[0];
        udot[0] += n[0]*dvars[0];
        uddot[0] += n[0]*ddvars[0];

        ux[0] += nxi[0]*vars[0];
        ux[1] += nxi[1]*vars[0];

        vars++;
        dvars++;
        ddvars++;
        u++;
        udot++;
        uddot++;
        ux += 2;
      }

      n++;
      nxi += 2;
    }

    // Compute the Jacobian transformation
    TacsScalar detJ = inv2x2(Xd, J);

    // U,x = U,xi * J
    TacsScalar *ux = Ux;
    for ( int j = 0; j < vars_per_node; j++ ){
      TacsScalar ud[2];
      ud[0] = ux[0];
      ud[1] = ux[1];

      mat2x2MultTrans(J, ud, ux);

      ux += 2;
    }

    return detJ;
  }
  else if (num_params == 1){
    // Zero the values of the coordinate and its derivative
    X[0] = X[1] = X[2] = 0.0;
    Xd[0] = 0.0;

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
      Xpts += 3;

      // Add contributions to the derivatives of the displacements
      TacsScalar *u = U, *udot = Udot, *uddot = Uddot, *ux = Ux;
      for ( int j = 0; j < vars_per_node; j++ ){
        u[0] += n[0]*vars[0];
        udot[0] += n[0]*dvars[0];
        uddot[0] += n[0]*ddvars[0];
        ux[0] += nxi[0]*vars[0];

        vars++;
        dvars++;
        ddvars++;
        u++;
        udot++;
        uddot++;
        ux++;
      }

      n++;
      nxi++;
    }

    return Xd[0];
  }

  return 0.0;
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
  const int num_params = getNumParameters();
  const int num_nodes = getNumNodes();
  double N[MAX_BASIS_SIZE], Nxi[3*MAX_BASIS_SIZE];
  computeBasisGradient(pt, N, Nxi);

  addWeakFormResidual(num_params, num_nodes, N, Nxi, weight, J,
                      vars_per_node, DUt, DUx, res);
}

void TACSElementBasis::addWeakFormResidual( const int num_params,
                                            const int num_nodes,
                                            const double N[],
                                            const double Nxi[],
                                            TacsScalar weight,
                                            const TacsScalar J[],
                                            const int vars_per_node,
                                            const TacsScalar DUt[],
                                            const TacsScalar DUx[],
                                            TacsScalar res[] ){
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
  const double *nxi = Nxi;
  r = res;

  if (num_params == 3){
    for ( int i = 0; i < num_nodes; i++ ){
      TacsScalar nx[3];
      mat3x3MultTrans(J, nxi, nx);
      for ( int ii = 0; ii < vars_per_node; ii++ ){
        r[0] += weight*(n[0]*DUx[4*ii] + nx[0]*DUx[4*ii+1] +
                        nx[1]*DUx[4*ii+2] + nx[2]*DUx[4*ii+3]);
        r++;
      }
      n++;
      nxi += 3;
    }
  }
  else if (num_params == 2){
    for ( int i = 0; i < num_nodes; i++ ){
      TacsScalar nx[2];
      mat2x2MultTrans(J, nxi, nx);
      for ( int ii = 0; ii < vars_per_node; ii++ ){
        r[0] += weight*(n[0]*DUx[3*ii] + nx[0]*DUx[3*ii+1] + nx[1]*DUx[3*ii+2]);
        r++;
      }
      n++;
      nxi += 2;
    }
  }
  else if (num_params == 1){
    for ( int i = 0; i < num_nodes; i++ ){
      TacsScalar nx = J[0]*nxi[0];
      for ( int ii = 0; ii < vars_per_node; ii++ ){
        r[0] += weight*(n[0]*DUx[2*ii] + nx*DUx[2*ii+1]);
        r++;
      }
      n++;
      nxi++;
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
                                            double alpha,
                                            double beta,
                                            double gamma,
                                            const int DDUt_nnz,
                                            const int *DDUt_pairs,
                                            const TacsScalar *DDUt,
                                            const int DDUx_nnz,
                                            const int *DDUx_pairs,
                                            const TacsScalar *DDUx,
                                            TacsScalar *res,
                                            TacsScalar *mat ){
  const int num_params = getNumParameters();
  const int num_nodes = getNumNodes();
  double N[MAX_BASIS_SIZE], Nxi[3*MAX_BASIS_SIZE];
  computeBasisGradient(pt, N, Nxi);

  // Convert the array Nxi (derivative of the shape functions w.r.t. the
  // parametric coordinates) to Nx (the derivative w.r.t. the nodal
  // coordinates). In other words, apply the coordinate transform
  // to the derivatives of the shape functions to avoid repeatedly
  // applying the transformation later on. This does not incurr extra
  // overhead during the residual computation but does here..
#ifdef TACS_USE_COMPLEX
  TacsScalar Nx[3*MAX_BASIS_SIZE];

  if (num_params == 3){
    for ( int i = 0; i < num_nodes; i++ ){
      mat3x3MultTrans(J, &Nxi[3*i], &Nx[3*i]);
    }
  }
  else if (num_params == 2){
    for ( int i = 0; i < num_nodes; i++ ){
      mat2x2MultTrans(J, &Nxi[2*i], &Nx[2*i]);
    }
  }
  else {
    for ( int i = 0; i < num_nodes; i++ ){
      Nx[i] = J[0]*Nxi[0];
    }
  }

  addWeakFormJacobian(num_params, num_nodes, N, Nx, weight, J,
                      vars_per_node, DUt, DUx, alpha, beta, gamma,
                      DDUt_nnz, DDUt_pairs, DDUt,
                      DDUx_nnz, DDUx_pairs, DDUx, res, mat);
#else // Real code
  if (num_params == 3){
    for ( int i = 0; i < num_nodes; i++ ){
      TacsScalar nx[3];
      mat3x3MultTrans(J, &Nxi[3*i], nx);
      Nxi[3*i] = nx[0];
      Nxi[3*i+1] = nx[1];
      Nxi[3*i+2] = nx[2];
    }
  }
  else if (num_params == 2){
    for ( int i = 0; i < num_nodes; i++ ){
      TacsScalar nx[2];
      mat2x2MultTrans(J, &Nxi[2*i], nx);
      Nxi[2*i] = nx[0];
      Nxi[2*i+1] = nx[1];
    }
  }
  else {
    for ( int i = 0; i < num_nodes; i++ ){
      Nxi[i] *= J[0];
    }
  }

  addWeakFormJacobian(num_params, num_nodes, N, Nxi, weight, J,
                      vars_per_node, DUt, DUx, alpha, beta, gamma,
                      DDUt_nnz, DDUt_pairs, DDUt,
                      DDUx_nnz, DDUx_pairs, DDUx, res, mat);
#endif // TACS_USE_COMPLEX
}

void TACSElementBasis::addWeakFormJacobian( const int num_params,
                                            const int num_nodes,
                                            const double N[],
                                            const TacsScalar Nx[],
                                            TacsScalar weight,
                                            const TacsScalar J[],
                                            const int vars_per_node,
                                            const TacsScalar DUt[],
                                            const TacsScalar DUx[],
                                            double alpha,
                                            double beta,
                                            double gamma,
                                            const int DDUt_nnz,
                                            const int *DDUt_pairs,
                                            const TacsScalar *DDUt,
                                            const int DDUx_nnz,
                                            const int *DDUx_pairs,
                                            const TacsScalar *DDUx,
                                            TacsScalar *res,
                                            TacsScalar *mat ){
  const int num_vars = num_nodes*vars_per_node;
  const int DDUt_size = 3*vars_per_node;
  const int DDUx_size = (num_params+1)*vars_per_node;

  // Add contributions from DUt
  const double *n = N;
  if (res){
    TacsScalar *r = res;
    for ( int i = 0; i < num_nodes; i++ ){
      for ( int ii = 0; ii < vars_per_node; ii++ ){
        r[0] += weight*n[0]*(DUt[3*ii] + DUt[3*ii+1] + DUt[3*ii+2]);
        r++;
      }
      n++;
    }
  }

  // Decode whether to use the dense or sparse representation
  if (DDUt_nnz < 0){
    for ( int ix = 0; ix < DDUt_size; ix++ ){
      for ( int jx = 0; jx < DDUt_size; jx++ ){
        const int ddut_index = DDUt_size*ix + jx;
        if (DDUt[ddut_index] != 0.0){
          TacsScalar scale = weight*DDUt[ddut_index];

          if (jx % 3 == 0){ scale *= alpha; }
          else if (jx % 3 == 1){ scale *= beta; }
          else { scale *= gamma; }

          const int init_index = num_vars*(ix/3) + (jx/3);
          double *M = &mat[init_index];
          for ( int i = 0; i < num_nodes; i++ ){
            for ( int j = 0; j < num_nodes; j++ ){
              M[0] += scale*N[i]*N[j];
              M += vars_per_node;
            }

            M += num_vars*(vars_per_node-1);
          }
        }
      }
    }
  }
  else {
    // Use a sparse representation
    for ( int ii = 0; ii < DDUt_nnz; ii++ ){
      int ix = DDUt_pairs[2*ii];
      int jx = DDUt_pairs[2*ii+1];

      if (DDUt[ii] != 0.0){
        TacsScalar scale = weight*DDUt[ii];

        if (jx % 3 == 0){ scale *= alpha; }
        else if (jx % 3 == 1){ scale *= beta; }
        else { scale *= gamma; }

        const int init_index = num_vars*(ix/3) + (jx/3);
        double *M = &mat[init_index];
        for ( int i = 0; i < num_nodes; i++ ){
          for ( int j = 0; j < num_nodes; j++ ){
            M[0] += scale*N[i]*N[j];
            M += vars_per_node;
          }

          M += num_vars*(vars_per_node-1);
        }
      }
    }
  }

  // Add contributions from DUx
  n = N;
  const double *nx = Nx;

  if (res){
    TacsScalar *r = res;
    if (num_params == 3){
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
    else if (num_params == 2){
      for ( int i = 0; i < num_nodes; i++ ){
        for ( int ii = 0; ii < vars_per_node; ii++ ){
          r[0] += weight*(n[0]*DUx[3*ii] + nx[0]*DUx[3*ii+1] +
                          nx[1]*DUx[3*ii+2]);
          r++;
        }
        n++;
        nx += 2;
      }
    }
    else if (num_params == 1){
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

  // Decode whether to use the dense or sparse representation
  if (DDUx_nnz < 0){
    for ( int ix = 0; ix < DDUx_size; ix++ ){
      for ( int jx = 0; jx < DDUx_size; jx++ ){
        const int ddux_index = DDUx_size*ix + jx;
        if (DDUx[ddux_index] != 0.0){
          TacsScalar scale = alpha*weight*DDUx[ddux_index];

          // Compute the offset point in the matrix
          const int init_index = (num_vars*(ix/(num_params+1)) +
                                  (jx/(num_params+1)));
          double *M = &mat[init_index];

          // Evaluate which shape function to use
          int incr = 1;
          const TacsScalar *ni = NULL;
          if (ix % (num_params+1) == 0){
            ni = N;
          }
          else {
            incr = num_params;
            ni = &Nx[(ix % (num_params+1))-1];
          }

          for ( int i = 0; i < num_nodes; i++ ){
            int jncr = 1;
            const TacsScalar *nj = NULL;
            if (jx % (num_params+1) == 0){
              nj = N;
            }
            else {
              jncr = num_params;
              nj = &Nx[(jx % (num_params+1))-1];
            }

            for ( int j = 0; j < num_nodes; j++ ){
              M[0] += scale*ni[0]*nj[0];
              M += vars_per_node;
              nj += jncr;
            }
            ni += incr;

            M += num_vars*(vars_per_node-1);
          }
        }
      }
    }
  }
  else {
    // Use a sparse representation
    for ( int ii = 0; ii < DDUx_nnz; ii++ ){
      int ix = DDUx_pairs[2*ii];
      int jx = DDUx_pairs[2*ii+1];

      if (DDUx[ii] != 0.0){
        TacsScalar scale = alpha*weight*DDUx[ii];

        // Compute the offset point in the matrix
        const int init_index = (num_vars*(ix/(num_params+1)) +
                                (jx/(num_params+1)));
        double *M = &mat[init_index];

        // Evaluate which shape function to use
        int incr = 1;
        const TacsScalar *ni = NULL;
        if (ix % (num_params+1) == 0){
          ni = N;
        }
        else {
          incr = num_params;
          ni = &Nx[(ix % (num_params+1))-1];
        }

        for ( int i = 0; i < num_nodes; i++ ){
          int jncr = 1;
          const TacsScalar *nj = NULL;
          if (jx % (num_params+1) == 0){
            nj = N;
          }
          else {
            jncr = num_params;
            nj = &Nx[(jx % (num_params+1))-1];
          }

          for ( int j = 0; j < num_nodes; j++ ){
            M[0] += scale*ni[0]*nj[0];
            M += vars_per_node;
            nj += jncr;
          }
          ni += incr;

          M += num_vars*(vars_per_node-1);
        }
      }
    }
  }
}
