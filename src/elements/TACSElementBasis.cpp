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
  Get the layout type
*/
ElementLayout TACSElementBasis::getLayoutType(){
  return TACS_LAYOUT_NONE;
}

/*
  Get the visualization parametric point
*/
void TACSElementBasis::getVisPoint( int n, double pt[] ){
  if (getNumParameters() == 3){
    pt[0] = pt[1] = pt[2] = 0.0;
  }
  else if (getNumParameters() == 2){
    pt[0] = pt[1] = 0.0;
  }
  else {
    pt[0] = 0.0;
  }
}

TacsScalar TACSElementBasis::getFaceNormal( int face, int n,
                                            const TacsScalar Xpts[],
                                            TacsScalar X[],
                                            TacsScalar Xd[],
                                            TacsScalar normal[] ){
  const int num_params = getNumParameters();

  double pt[3];
  double tangents[2*3];
  getFaceQuadraturePoint(face, n, pt, tangents);

  // Compute the interpolation on the face
  interpFaceFields(face, n, pt, 3, Xpts, 1, X);
  interpFaceFieldsGrad(face, n, pt, 3, Xpts, Xd);

  if (num_params == 3){
    // Compute the tangent directions
    TacsScalar t1[3], t2[3];
    mat3x3Mult(Xd, &tangents[0], t1);
    mat3x3Mult(Xd, &tangents[3], t2);

    // Find the cross product
    crossProduct(1.0, t1, t2, normal);

    TacsScalar A = sqrt(vec3Dot(normal, normal));
    vec3Scale(1.0/A, normal);

    return A;
  }
  else if (num_params == 2){
    // Compute the tangent direction
    TacsScalar t[2];
    mat2x2Mult(Xd, tangents, t);

    // Compute the norm of the vector
    TacsScalar A = sqrt(vec2Dot(t, t));

    normal[0] = t[1]/A;
    normal[1] = -t[0]/A;

    return A;
  }
  return 0.0;
}

void TACSElementBasis::addFaceNormalXptSens( int face, int n,
                                             const TacsScalar A,
                                             const TacsScalar Xd[],
                                             const TacsScalar normal[],
                                             const TacsScalar dfdA,
                                             const TacsScalar dfdX[],
                                             const TacsScalar dfdXd[],
                                             const TacsScalar dfdn[],
                                             TacsScalar dfdXpts[] ){
  const int num_params = getNumParameters();

  double pt[3];
  double tangents[2*3];
  getFaceQuadraturePoint(face, n, pt, tangents);

  if (num_params == 3){
    TacsScalar Ainv = 1.0/A;

    // Compute dfda = 1.0/A*dfdn(I - (n*n^{T}))
    TacsScalar dfda[3];
    dfda[0] = Ainv*dfdn[0];
    dfda[1] = Ainv*dfdn[1];
    dfda[2] = Ainv*dfdn[2];
    vec3Axpy(dfdA - Ainv*vec3Dot(normal, dfdn), normal, dfda);

    // Compute the tangent directions
    TacsScalar t1[3], t2[3];
    mat3x3Mult(Xd, &tangents[0], t1);
    mat3x3Mult(Xd, &tangents[3], t2);

    TacsScalar dfdt1[3], dfdt2[3];
    crossProduct(1.0, t2, dfda, dfdt1);
    crossProduct(1.0, dfda, t1, dfdt2);

    TacsScalar t[9];
    if (dfdXd){
      memcpy(t, dfdXd, 9*sizeof(TacsScalar));
      vec3x3OuterAdd(1.0, dfdt1, &tangents[0], t);
      vec3x3OuterAdd(1.0, dfdt2, &tangents[3], t);
    }
    else {
      vec3x3Outer(dfdt1, &tangents[0], t);
      vec3x3OuterAdd(1.0, dfdt2, &tangents[3], t);
    }

    // Loop over the basis functions
    addInterpFaceFieldsGradTranspose(face, n, pt, 3, t, dfdXpts);

    if (dfdX){
      addInterpFaceFieldsTranspose(face, n, pt, 1, dfdX, 3, dfdXpts);
    }
  }
  else if (num_params == 2){
    TacsScalar Ainv = 1.0/A;

    TacsScalar dfda[2];
    dfda[0] = Ainv*dfdn[0];
    dfda[1] = Ainv*dfdn[1];
    vec2Axpy(dfdA - Ainv*vec2Dot(normal, dfdn), normal, dfda);

    TacsScalar dfdt[2];
    dfdt[0] = -dfda[1];
    dfdt[1] = dfda[0];

    TacsScalar t[6] = {0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0};
    if (dfdXd){
      memcpy(t, dfdXd, 4*sizeof(TacsScalar));
      vec2x2OuterAdd(1.0, dfdt, &tangents[0], t);
    }
    else {
      vec2x2Outer(dfdt, &tangents[0], t);
    }

    // Loop over the basis functions
    addInterpFaceFieldsGradTranspose(face, n, pt, 3, t, dfdXpts);

    if (dfdX){
      addInterpFaceFieldsTranspose(face, n, pt, 1, dfdX, 3, dfdXpts);
    }
  }
}

TacsScalar TACSElementBasis::getJacobianTransform( int n, const double pt[],
                                                   const TacsScalar Xpts[],
                                                   TacsScalar Xd[],
                                                   TacsScalar J[] ){
  const int num_params = getNumParameters();

  interpFieldsGrad(n, pt, 3, Xpts, Xd);

  if (num_params == 3){
    // Compute the Jacobian transformation
    TacsScalar det = inv3x3(Xd, J);

    return det;
  }
  else if (num_params == 2){
    // Compute the Jacobian transformation
    TacsScalar det = inv2x2(Xd, J);

    return det;
  }
  else if (num_params == 1){
    J[0] = 1.0/Xd[0];

    return Xd[0];
  }

  return 0.0;
}

void TACSElementBasis::addJacobianTransformXptSens( int n, const double pt[],
                                                    const TacsScalar Xd[],
                                                    const TacsScalar J[],
                                                    TacsScalar dfddetXd,
                                                    const TacsScalar dfdXd[],
                                                    const TacsScalar dfdJ[],
                                                    TacsScalar dfdXpts[] ){
  const int num_params = getNumParameters();

  if (num_params == 3){
    // Compute t = d(detXd)/d(Xd)
    TacsScalar t[9];
    det3x3Sens(Xd, t);

    // Multiply the derivative by df/d(detXd)
    TacsScalar *T = t;
    for ( int i = 0; i < 9; i++, T++ ){
      T[0] *= dfddetXd;
    }

    // If dfdXd is supplied, add it to the derivative
    if (dfdXd){
      for ( int i = 0; i < 9; i++ ){
        t[i] += dfdXd[i];
      }
    }

    // if df/d(J) is supplied, compute df/d(J)*d(J)/d(Xd) and add it
    // to the array t
    if (dfdJ){
      TacsScalar t2[9];
      inv3x3Sens(J, dfdJ, t2);
      for ( int i = 0; i < 9; i++ ){
        t[i] += t2[i];
      }
    }

    // Loop over each quadrature point for each basis function
    addInterpFieldsGradTranspose(n, pt, 3, t, dfdXpts);
  }
  else if (num_params == 2){
    // Compute t = d(detXd)/d(Xd)
    TacsScalar t[6] = {0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0};
    det2x2Sens(Xd, t);

    // Multiply the derivative by df/d(detXd)
    TacsScalar *T = t;
    for ( int i = 0; i < 4; i++, T++ ){
      T[0] *= dfddetXd;
    }

    // If dfdXd is supplied, add it to the derivative
    if (dfdXd){
      for ( int i = 0; i < 4; i++ ){
        t[i] += dfdXd[i];
      }
    }

    // if df/d(J) is supplied, compute df/d(J)*d(J)/d(Xd) and add it
    // to the array t
    if (dfdJ){
      TacsScalar t2[4];
      inv2x2Sens(J, dfdJ, t2);
      for ( int i = 0; i < 4; i++ ){
        t[i] += t2[i];
      }
    }

    // Loop over each quadrature point for each basis function
    addInterpFieldsGradTranspose(n, pt, 3, t, dfdXpts);
  }
}

/*
  Get the field values at the specified quadrature point
*/
void TACSElementBasis::getFieldValues( int n, const double pt[],
                                       const TacsScalar Xpts[],
                                       const int vars_per_node,
                                       const TacsScalar vars[],
                                       TacsScalar X[],
                                       TacsScalar U[] ){
  // Compute the interpolation point
  interpFields(n, pt, 3, Xpts, 1, X);
  interpFields(n, pt, vars_per_node, vars, 1, U);
}

/*
  Get the gradient of the field at the quadrature point
*/
TacsScalar TACSElementBasis::getFieldGradient( int n, const double pt[],
                                               const TacsScalar Xpts[],
                                               const int vars_per_node,
                                               const TacsScalar vars[],
                                               const TacsScalar dvars[],
                                               const TacsScalar ddvars[],
                                               TacsScalar X[],
                                               TacsScalar Xd[],
                                               TacsScalar J[],
                                               TacsScalar Ut[],
                                               TacsScalar Ud[],
                                               TacsScalar Ux[] ){
  const int num_params = getNumParameters();

  // Compute the interpolation point
  interpFields(n, pt, 3, Xpts, 1, X);

  // Compute derivative of the interpolation
  interpFieldsGrad(n, pt, 3, Xpts, Xd);

  if (vars && dvars && ddvars){
    interpFieldsGrad(n, pt, vars_per_node, vars, Ud);
    interpFields(n, pt, vars_per_node, vars, dvars, ddvars, Ut);
  }
  else {
    interpFieldsGrad(n, pt, vars_per_node, vars, Ud);
    interpFields(n, pt, vars_per_node, vars, 3, &Ut[0]);

    // Zero the time derivatives
    for ( int j = 0; j < vars_per_node; j++ ){
      Ut[3*j+1] = Ut[3*j+2] = 0.0;
    }
  }

  if (num_params == 3){
    // Compute the Jacobian transformation
    TacsScalar detXd = inv3x3(Xd, J);

    // U,x = U,xi * J
    TacsScalar *ux = Ux;
    for ( int j = 0; j < vars_per_node; j++ ){
      mat3x3MultTrans(J, &Ud[3*j], ux);
      ux += 3;
    }

    return detXd;
  }
  else if (num_params == 2){
    // Compute the Jacobian transformation
    TacsScalar detXd = inv2x2(Xd, J);

    // U,x = U,xi * J
    TacsScalar *ux = Ux;
    for ( int j = 0; j < vars_per_node; j++ ){
      mat2x2MultTrans(J, &Ud[2*j], ux);
      ux += 2;
    }

    return detXd;
  }
  else if (num_params == 1){
    J[0] = 1.0/Xd[0];
    for ( int j = 0; j < vars_per_node; j++ ){
      Ux[j] = J[0]*Ud[j];
    }

    return Xd[0];
  }

  return 0.0;
}

TacsScalar TACSElementBasis::getFieldGradient( int n, const double pt[],
                                               const TacsScalar Xpts[],
                                               const int vars_per_node,
                                               const TacsScalar vars[],
                                               const TacsScalar dvars[],
                                               const TacsScalar ddvars[],
                                               const TacsScalar psi[],
                                               TacsScalar X[],
                                               TacsScalar Xd[],
                                               TacsScalar J[],
                                               TacsScalar Ut[],
                                               TacsScalar Ud[],
                                               TacsScalar Ux[],
                                               TacsScalar Psi[],
                                               TacsScalar Psid[],
                                               TacsScalar Psix[] ){
  const int num_params = getNumParameters();

  // Compute the interpolation point
  interpFields(n, pt, 3, Xpts, 1, X);

  // Compute derivative of the interpolation
  interpFieldsGrad(n, pt, 3, Xpts, Xd);

  if (vars && dvars && ddvars && psi){
    interpFieldsGrad(n, pt, vars_per_node, vars, Ud);
    interpFieldsGrad(n, pt, vars_per_node, psi, Psid);
    interpFields(n, pt, vars_per_node, psi, 1, Psi);
    interpFields(n, pt, vars_per_node, vars, dvars, ddvars, Ut);
  }
  else {
    interpFieldsGrad(n, pt, vars_per_node, vars, Ud);
    interpFieldsGrad(n, pt, vars_per_node, psi, Psid);
    interpFields(n, pt, vars_per_node, psi, 1, Psi);
    interpFields(n, pt, vars_per_node, vars, 3, &Ut[0]);

    // Zero the time derivatives
    for ( int j = 0; j < vars_per_node; j++ ){
      Ut[3*j+1] = Ut[3*j+2] = 0.0;
    }
  }

  if (num_params == 3){
    // Compute the Jacobian transformation
    TacsScalar detXd = inv3x3(Xd, J);

    // U,x = U,xi * J
    TacsScalar *ux = Ux;
    for ( int j = 0; j < vars_per_node; j++ ){
      mat3x3MultTrans(J, &Ud[3*j], ux);
      ux += 3;
    }

    TacsScalar *px = Psix;
    for ( int j = 0; j < vars_per_node; j++ ){
      mat3x3MultTrans(J, &Psid[3*j], px);
      px += 3;
    }

    return detXd;
  }
  else if (num_params == 2){
    // Compute the Jacobian transformation
    TacsScalar detXd = inv2x2(Xd, J);

    // U,x = U,xi * J
    TacsScalar *ux = Ux;
    for ( int j = 0; j < vars_per_node; j++ ){
      mat2x2MultTrans(J, &Ud[2*j], ux);
      ux += 2;
    }

    TacsScalar *px = Psix;
    for ( int j = 0; j < vars_per_node; j++ ){
      mat2x2MultTrans(J, &Psid[2*j], px);
      px += 2;
    }

    return detXd;
  }
  else if (num_params == 1){
    J[0] = 1.0/Xd[0];
    for ( int j = 0; j < vars_per_node; j++ ){
      Ux[j] = J[0]*Ud[j];
    }

    return Xd[0];
  }

  return 0.0;
}

void TACSElementBasis::addFieldGradientSVSens( int n, const double pt[],
                                               const TacsScalar Xpts[],
                                               const int vars_per_node,
                                               const TacsScalar Xd[],
                                               const TacsScalar J[],
                                               const TacsScalar Ud[],
                                               const TacsScalar dfdUt[],
                                               TacsScalar dfdUx[],
                                               TacsScalar dfdu[] ){
  const int num_params = getNumParameters();

  if (num_params == 3){
    for ( int i = 0; i < vars_per_node; i++ ){
      TacsScalar jx[3];
      mat3x3Mult(J, &dfdUx[3*i], jx);
      dfdUx[3*i] = jx[0];
      dfdUx[3*i+1] = jx[1];
      dfdUx[3*i+2] = jx[2];
    }
  }
  else if (num_params == 2){
    for ( int i = 0; i < vars_per_node; i++ ){
      TacsScalar jx[2];
      mat2x2Mult(J, &dfdUx[2*i], jx);
      dfdUx[2*i] = jx[0];
      dfdUx[2*i+1] = jx[1];
    }
  }
  else {
    for ( int i = 0; i < vars_per_node; i++ ){
      dfdUx[i] *= J[0];
    }
  }

  addInterpFieldsTranspose(n, pt, 3, &dfdUt[0], vars_per_node, dfdu);
  addInterpFieldsTranspose(n, pt, 3, &dfdUt[1], vars_per_node, dfdu);
  addInterpFieldsTranspose(n, pt, 3, &dfdUt[2], vars_per_node, dfdu);
  addInterpFieldsGradTranspose(n, pt, vars_per_node, dfdUx, dfdu);
}

void TACSElementBasis::addFieldGradientXptSens( int n, const double pt[],
                                                const TacsScalar Xpts[],
                                                const int vars_per_node,
                                                const TacsScalar Xd[],
                                                const TacsScalar J[],
                                                const TacsScalar Ud[],
                                                const TacsScalar dfddetXd,
                                                const TacsScalar dfdX[],
                                                const TacsScalar dfdXd[],
                                                const TacsScalar dfdJ[],
                                                const TacsScalar dfdUx[],
                                                TacsScalar dfdXpts[] ){
  const int num_params = getNumParameters();

  if (num_params == 3){
    // Compute t = d(detXd)/d(Xd)
    TacsScalar t[9];
    det3x3Sens(Xd, t);

    // Multiply the derivative by df/d(detXd)
    for ( int i = 0; i < 9; i++ ){
      t[i] *= dfddetXd;
    }

    // If dfdXd is supplied, add it to the derivative
    if (dfdXd){
      for ( int i = 0; i < 9; i++ ){
        t[i] += dfdXd[i];
      }
    }

    // if df/d(J) is supplied, compute df/d(J)*d(J)/d(Xd) and add it
    // to the array t
    if (dfdJ && dfdUx){
      // dfdJ = dfdUx*dUx/dJ, where dUx/dJ = Ud
      TacsScalar df[9];
      memcpy(df, dfdJ, 9*sizeof(TacsScalar));
      for ( int k = 0; k < vars_per_node; k++ ){
        df[0] += Ud[3*k]*dfdUx[3*k];
        df[1] += Ud[3*k]*dfdUx[3*k+1];
        df[2] += Ud[3*k]*dfdUx[3*k+2];

        df[3] += Ud[3*k+1]*dfdUx[3*k];
        df[4] += Ud[3*k+1]*dfdUx[3*k+1];
        df[5] += Ud[3*k+1]*dfdUx[3*k+2];

        df[6] += Ud[3*k+2]*dfdUx[3*k];
        df[7] += Ud[3*k+2]*dfdUx[3*k+1];
        df[8] += Ud[3*k+2]*dfdUx[3*k+2];
      }

      TacsScalar t2[9];
      inv3x3Sens(J, df, t2);
      for ( int i = 0; i < 9; i++ ){
        t[i] += t2[i];
      }
    }
    else if (dfdJ){
      TacsScalar t2[9];
      inv3x3Sens(J, dfdJ, t2);
      for ( int i = 0; i < 9; i++ ){
        t[i] += t2[i];
      }
    }
    else if (dfdUx){
      // dfdJ = dfdUx*dUx/dJ, where dUx/dJ = Ud
      TacsScalar df[9];
      memset(df, 0, 9*sizeof(TacsScalar));
      for ( int k = 0; k < vars_per_node; k++ ){
        df[0] += Ud[3*k]*dfdUx[3*k];
        df[1] += Ud[3*k]*dfdUx[3*k+1];
        df[2] += Ud[3*k]*dfdUx[3*k+2];

        df[3] += Ud[3*k+1]*dfdUx[3*k];
        df[4] += Ud[3*k+1]*dfdUx[3*k+1];
        df[5] += Ud[3*k+1]*dfdUx[3*k+2];

        df[6] += Ud[3*k+2]*dfdUx[3*k];
        df[7] += Ud[3*k+2]*dfdUx[3*k+1];
        df[8] += Ud[3*k+2]*dfdUx[3*k+2];
      }

      TacsScalar t2[9];
      inv3x3Sens(J, df, t2);
      for ( int i = 0; i < 9; i++ ){
        t[i] += t2[i];
      }
    }

    // Loop over the basis functions
    addInterpFieldsGradTranspose(n, pt, 3, t, dfdXpts);

    if (dfdX){
      addInterpFieldsTranspose(n, pt, 1, dfdX, 3, dfdXpts);
    }
  }
  else if (num_params == 2){
    // Compute t = d(detXd)/d(Xd)
    TacsScalar t[6] = {0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0};
    det2x2Sens(Xd, t);

    // Multiply the derivative by df/d(detXd)
    TacsScalar *T = t;
    for ( int i = 0; i < 4; i++, T++ ){
      T[0] *= dfddetXd;
    }

    // If dfdXd is supplied, add it to the derivative
    if (dfdXd){
      for ( int i = 0; i < 4; i++ ){
        t[i] += dfdXd[i];
      }
    }

    // If df/d(J) is supplied, compute df/d(J)*d(J)/d(Xd) and add it
    // to the array t
    if (dfdJ && dfdUx){
      // dfdJ = dfdUx*dUx/dJ, where dUx/dJ = Ud
      TacsScalar df[4];
      memcpy(df, dfdJ, 4*sizeof(TacsScalar));
      for ( int k = 0; k < vars_per_node; k++ ){
        df[0] += Ud[2*k]*dfdUx[2*k];
        df[1] += Ud[2*k]*dfdUx[2*k+1];

        df[2] += Ud[2*k+1]*dfdUx[2*k];
        df[3] += Ud[2*k+1]*dfdUx[2*k+1];
      }

      TacsScalar t2[4];
      inv2x2Sens(J, df, t2);
      for ( int i = 0; i < 4; i++ ){
        t[i] += t2[i];
      }
    }
    else if (dfdJ){
      TacsScalar t2[4];
      inv2x2Sens(J, dfdJ, t2);
      for ( int i = 0; i < 4; i++ ){
        t[i] += t2[i];
      }
    }
    else if (dfdUx){
      TacsScalar df[4];
      memset(df, 0, 4*sizeof(TacsScalar));
      for ( int k = 0; k < vars_per_node; k++ ){
        df[0] += Ud[2*k]*dfdUx[2*k];
        df[1] += Ud[2*k]*dfdUx[2*k+1];

        df[2] += Ud[2*k+1]*dfdUx[2*k];
        df[3] += Ud[2*k+1]*dfdUx[2*k+1];
      }

      TacsScalar t2[4];
      inv2x2Sens(J, df, t2);
      for ( int i = 0; i < 4; i++ ){
        t[i] += t2[i];
      }
    }

    // Loop over each basis function
    addInterpFieldsGradTranspose(n, pt, 3, t, dfdXpts);

    if (dfdX){
      addInterpFieldsTranspose(n, pt, 1, dfdX, 3, dfdXpts);
    }
  }
}

void TACSElementBasis::addFieldGradientXptSens( int n, const double pt[],
                                                const TacsScalar Xpts[],
                                                const int vars_per_node,
                                                const TacsScalar Xd[],
                                                const TacsScalar J[],
                                                const TacsScalar Ud[],
                                                const TacsScalar Psid[],
                                                const TacsScalar dfddetXd,
                                                const TacsScalar dfdX[],
                                                const TacsScalar dfdXd[],
                                                const TacsScalar dfdJ[],
                                                const TacsScalar dfdUx[],
                                                const TacsScalar dfdPsix[],
                                                TacsScalar dfdXpts[] ){
  const int num_params = getNumParameters();

  if (num_params == 3){
    // Compute t = d(detXd)/d(Xd)
    TacsScalar t[9];
    det3x3Sens(Xd, t);

    // Multiply the derivative by df/d(detXd)
    for ( int i = 0; i < 9; i++ ){
      t[i] *= dfddetXd;
    }

    // If dfdXd is supplied, add it to the derivative
    if (dfdXd){
      for ( int i = 0; i < 9; i++ ){
        t[i] += dfdXd[i];
      }
    }

    // if df/d(J) is supplied, compute df/d(J)*d(J)/d(Xd) and add it
    // to the array t
    if (dfdJ && dfdUx && dfdPsix){
      // dfdJ = dfdUx*dUx/dJ, where dUx/dJ = Ud
      TacsScalar df[9];
      memcpy(df, dfdJ, 9*sizeof(TacsScalar));
      for ( int k = 0; k < vars_per_node; k++ ){
        df[0] += Ud[3*k]*dfdUx[3*k];
        df[1] += Ud[3*k]*dfdUx[3*k+1];
        df[2] += Ud[3*k]*dfdUx[3*k+2];

        df[3] += Ud[3*k+1]*dfdUx[3*k];
        df[4] += Ud[3*k+1]*dfdUx[3*k+1];
        df[5] += Ud[3*k+1]*dfdUx[3*k+2];

        df[6] += Ud[3*k+2]*dfdUx[3*k];
        df[7] += Ud[3*k+2]*dfdUx[3*k+1];
        df[8] += Ud[3*k+2]*dfdUx[3*k+2];
      }

      for ( int k = 0; k < vars_per_node; k++ ){
        df[0] += Psid[3*k]*dfdPsix[3*k];
        df[1] += Psid[3*k]*dfdPsix[3*k+1];
        df[2] += Psid[3*k]*dfdPsix[3*k+2];

        df[3] += Psid[3*k+1]*dfdPsix[3*k];
        df[4] += Psid[3*k+1]*dfdPsix[3*k+1];
        df[5] += Psid[3*k+1]*dfdPsix[3*k+2];

        df[6] += Psid[3*k+2]*dfdPsix[3*k];
        df[7] += Psid[3*k+2]*dfdPsix[3*k+1];
        df[8] += Psid[3*k+2]*dfdPsix[3*k+2];
      }

      TacsScalar t2[9];
      inv3x3Sens(J, df, t2);
      for ( int i = 0; i < 9; i++ ){
        t[i] += t2[i];
      }
    }
    else if (dfdUx && dfdPsix){
      // dfdJ = dfdUx*dUx/dJ, where dUx/dJ = Ud
      TacsScalar df[9];
      memset(df, 0, 9*sizeof(TacsScalar));
      for ( int k = 0; k < vars_per_node; k++ ){
        df[0] += Ud[3*k]*dfdUx[3*k];
        df[1] += Ud[3*k]*dfdUx[3*k+1];
        df[2] += Ud[3*k]*dfdUx[3*k+2];

        df[3] += Ud[3*k+1]*dfdUx[3*k];
        df[4] += Ud[3*k+1]*dfdUx[3*k+1];
        df[5] += Ud[3*k+1]*dfdUx[3*k+2];

        df[6] += Ud[3*k+2]*dfdUx[3*k];
        df[7] += Ud[3*k+2]*dfdUx[3*k+1];
        df[8] += Ud[3*k+2]*dfdUx[3*k+2];
      }

      for ( int k = 0; k < vars_per_node; k++ ){
        df[0] += Psid[3*k]*dfdPsix[3*k];
        df[1] += Psid[3*k]*dfdPsix[3*k+1];
        df[2] += Psid[3*k]*dfdPsix[3*k+2];

        df[3] += Psid[3*k+1]*dfdPsix[3*k];
        df[4] += Psid[3*k+1]*dfdPsix[3*k+1];
        df[5] += Psid[3*k+1]*dfdPsix[3*k+2];

        df[6] += Psid[3*k+2]*dfdPsix[3*k];
        df[7] += Psid[3*k+2]*dfdPsix[3*k+1];
        df[8] += Psid[3*k+2]*dfdPsix[3*k+2];
      }

      TacsScalar t2[9];
      inv3x3Sens(J, df, t2);
      for ( int i = 0; i < 9; i++ ){
        t[i] += t2[i];
      }
    }
    else if (dfdJ){
      TacsScalar t2[9];
      inv3x3Sens(J, dfdJ, t2);
      for ( int i = 0; i < 9; i++ ){
        t[i] += t2[i];
      }
    }

    // Loop over each basis function
    addInterpFieldsGradTranspose(n, pt, 3, t, dfdXpts);

    if (dfdX){
      addInterpFieldsTranspose(n, pt, 1, dfdX, 3, dfdXpts);
    }
  }
  else if (num_params == 2){
    // Compute t = d(detXd)/d(Xd)
    TacsScalar t[6] = {0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0};
    det2x2Sens(Xd, t);

    // Multiply the derivative by df/d(detXd)
    TacsScalar *T = t;
    for ( int i = 0; i < 4; i++, T++ ){
      T[0] *= dfddetXd;
    }

    // If dfdXd is supplied, add it to the derivative
    if (dfdXd){
      for ( int i = 0; i < 4; i++ ){
        t[i] += dfdXd[i];
      }
    }

    // If df/d(J) is supplied, compute df/d(J)*d(J)/d(Xd) and add it
    // to the array t
    if (dfdJ && dfdUx && dfdPsix){
      // dfdJ = dfdUx*dUx/dJ, where dUx/dJ = Ud
      TacsScalar df[4];
      memcpy(df, dfdJ, 4*sizeof(TacsScalar));
      for ( int k = 0; k < vars_per_node; k++ ){
        df[0] += Ud[2*k]*dfdUx[2*k];
        df[1] += Ud[2*k]*dfdUx[2*k+1];

        df[2] += Ud[2*k+1]*dfdUx[2*k];
        df[3] += Ud[2*k+1]*dfdUx[2*k+1];
      }

      for ( int k = 0; k < vars_per_node; k++ ){
        df[0] += Psid[2*k]*dfdPsix[2*k];
        df[1] += Psid[2*k]*dfdPsix[2*k+1];

        df[2] += Psid[2*k+1]*dfdPsix[2*k];
        df[3] += Psid[2*k+1]*dfdPsix[2*k+1];
      }

      TacsScalar t2[4];
      inv2x2Sens(J, df, t2);
      for ( int i = 0; i < 4; i++ ){
        t[i] += t2[i];
      }
    }
    else if (dfdUx && dfdPsix){
      TacsScalar df[4];
      memset(df, 0, 4*sizeof(TacsScalar));
      for ( int k = 0; k < vars_per_node; k++ ){
        df[0] += Ud[2*k]*dfdUx[2*k];
        df[1] += Ud[2*k]*dfdUx[2*k+1];

        df[2] += Ud[2*k+1]*dfdUx[2*k];
        df[3] += Ud[2*k+1]*dfdUx[2*k+1];
      }

      for ( int k = 0; k < vars_per_node; k++ ){
        df[0] += Psid[2*k]*dfdPsix[2*k];
        df[1] += Psid[2*k]*dfdPsix[2*k+1];

        df[2] += Psid[2*k+1]*dfdPsix[2*k];
        df[3] += Psid[2*k+1]*dfdPsix[2*k+1];
      }

      TacsScalar t2[4];
      inv2x2Sens(J, df, t2);
      for ( int i = 0; i < 4; i++ ){
        t[i] += t2[i];
      }
    }
    else if (dfdJ){
      TacsScalar t2[4];
      inv2x2Sens(J, dfdJ, t2);
      for ( int i = 0; i < 4; i++ ){
        t[i] += t2[i];
      }
    }

    // Loop over each basis function
    addInterpFieldsGradTranspose(n, pt, 3, t, dfdXpts);

    if (dfdX){
      addInterpFieldsTranspose(n, pt, 1, dfdX, 3, dfdXpts);
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
                                            TacsScalar DUt[],
                                            TacsScalar DUx[],
                                            TacsScalar res[] ){
  const int num_params = getNumParameters();

  // Add contributions from DUt
  for ( int i = 0; i < 3*vars_per_node; i++ ){
    DUt[i] *= weight;
  }
  addInterpFieldsTranspose(n, pt, 3, &DUt[0], vars_per_node, res);
  addInterpFieldsTranspose(n, pt, 3, &DUt[1], vars_per_node, res);
  addInterpFieldsTranspose(n, pt, 3, &DUt[2], vars_per_node, res);

  if (num_params == 3){
    for ( int i = 0; i < vars_per_node; i++ ){
      TacsScalar dx[3];
      mat3x3Mult(J, &DUx[3*i], dx);
      DUx[3*i] = weight*dx[0];
      DUx[3*i+1] = weight*dx[1];
      DUx[3*i+2] = weight*dx[2];
    }
  }
  else if (num_params == 2){
    for ( int i = 0; i < vars_per_node; i++ ){
      TacsScalar dx[2];
      mat2x2Mult(J, &DUx[2*i], dx);
      DUx[2*i] = weight*dx[0];
      DUx[2*i+1] = weight*dx[1];
    }
  }
  else {
    for ( int i = 0; i < vars_per_node; i++ ){
      DUx[i] *= J[0]*weight;
    }
  }

  addInterpFieldsGradTranspose(n, pt, vars_per_node, DUx, res);
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
                                            TacsScalar alpha,
                                            TacsScalar beta,
                                            TacsScalar gamma,
                                            const int Jac_nnz,
                                            const int *Jac_pairs,
                                            const TacsScalar *Jac,
                                            TacsScalar *res,
                                            TacsScalar *mat ){
  const int num_params = getNumParameters();
  const int num_nodes = getNumNodes();
  double N[256], Nxi[3*256];
  computeBasisGradient(pt, N, Nxi);

  // Convert the array Nxi (derivative of the shape functions w.r.t. the
  // parametric coordinates) to Nx (the derivative w.r.t. the nodal
  // coordinates). In other words, apply the coordinate transform
  // to the derivatives of the shape functions to avoid repeatedly
  // applying the transformation later on. This does not incurr extra
  // overhead during the residual computation but does here..
#ifdef TACS_USE_COMPLEX
  TacsScalar Nc[256];
  TacsScalar Nx[3*256];

  for ( int i = 0; i < num_nodes; i++ ){
    Nc[i] = N[i];
  }

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

  addWeakFormJacobian(num_params, num_nodes, Nc, Nx, weight, J,
                      vars_per_node, DUt, DUx, alpha, beta, gamma,
                      Jac_nnz, Jac_pairs, Jac, res, mat);
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
                      Jac_nnz, Jac_pairs, Jac, res, mat);
#endif // TACS_USE_COMPLEX
}

void TACSElementBasis::addWeakFormJacobian( const int num_params,
                                            const int num_nodes,
                                            const TacsScalar N[],
                                            const TacsScalar Nx[],
                                            TacsScalar weight,
                                            const TacsScalar J[],
                                            const int vars_per_node,
                                            const TacsScalar DUt[],
                                            const TacsScalar DUx[],
                                            TacsScalar alpha,
                                            TacsScalar beta,
                                            TacsScalar gamma,
                                            const int Jac_nnz,
                                            const int *Jac_pairs,
                                            const TacsScalar *Jac,
                                            TacsScalar *res,
                                            TacsScalar *mat ){
  const int num_vars = num_nodes*vars_per_node;
  const int Jac_size = (3 + num_params)*vars_per_node;

  // Add contributions from DUt
  const TacsScalar *n = N;
  const TacsScalar *nx = Nx;
  if (res){
    // Add the residual contribution from the time derivative components
    TacsScalar *r = res;
    for ( int i = 0; i < num_nodes; i++ ){
      for ( int ii = 0; ii < vars_per_node; ii++ ){
        r[0] += weight*n[0]*(DUt[3*ii] +
                             DUt[3*ii+1] +
                             DUt[3*ii+2]);
        r++;
      }
      n++;
    }

    // Add the residual contribution from the spatial derivative components
    r = res;
    if (num_params == 3){
      for ( int i = 0; i < num_nodes; i++ ){
        for ( int ii = 0; ii < vars_per_node; ii++ ){
          r[0] += weight*(nx[0]*DUx[3*ii] +
                          nx[1]*DUx[3*ii+1] +
                          nx[2]*DUx[3*ii+2]);
          r++;
        }
        nx += 3;
      }
    }
    else if (num_params == 2){
      for ( int i = 0; i < num_nodes; i++ ){
        for ( int ii = 0; ii < vars_per_node; ii++ ){
          r[0] += weight*(nx[0]*DUx[2*ii] +
                          nx[1]*DUx[2*ii+1]);
          r++;
        }
        nx += 2;
      }
    }
    else if (num_params == 1){
      for ( int i = 0; i < num_nodes; i++ ){
        for ( int ii = 0; ii < vars_per_node; ii++ ){
          r[0] += weight*(nx[0]*DUx[ii]);
          r++;
        }
        nx++;
      }
    }
  }

  // Decode whether to use the dense or sparse representation
  if (Jac_nnz < 0){
    for ( int ix = 0; ix < Jac_nnz; ix++ ){
      for ( int jx = 0; jx < Jac_nnz; jx++ ){
        const int Jac_index = Jac_size*ix + jx;
        if (Jac[Jac_index] != 0.0){
          TacsScalar init_scale = weight*Jac[Jac_index];

          // Compute the offset point in the matrix
          const int init_index = (num_vars*(ix/(num_params+3)) +
                                  (jx/(num_params+3)));
          TacsScalar *M = &mat[init_index];

          // Evaluate which shape function to use
          int incr = 1;
          const TacsScalar *ni = NULL;
          if (ix % (num_params + 3) == 0){
            ni = N;
          }
          else if (ix % (num_params + 3) == 1){
            ni = N;
          }
          else if (ix % (num_params + 3) == 2){
            ni = N;
          }
          else {
            incr = num_params;
            ni = &Nx[(ix % (num_params + 3)) - 3];
          }

          for ( int i = 0; i < num_nodes; i++ ){
            TacsScalar scale = 0.0;
            int jncr = 1;
            const TacsScalar *nj = NULL;
            if (jx % (num_params + 3) == 0){
              scale = alpha*init_scale;
              nj = N;
            }
            else if (jx % (num_params + 3) == 1){
              scale = beta*init_scale;
              nj = N;
            }
            else if (jx % (num_params + 3) == 2){
              scale = gamma*init_scale;
              nj = N;
            }
            else {
              scale = alpha*init_scale;
              jncr = num_params;
              nj = &Nx[(jx % (num_params + 3)) - 3];
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
    for ( int ii = 0; ii < Jac_nnz; ii++ ){
      int ix = Jac_pairs[2*ii];
      int jx = Jac_pairs[2*ii+1];

      if (Jac[ii] != 0.0){
        TacsScalar init_scale = weight*Jac[ii];

        // Compute the offset point in the matrix
        const int init_index = (num_vars*(ix/(num_params+3)) +
                                (jx/(num_params+3)));
        TacsScalar *M = &mat[init_index];

        // Evaluate which shape function to use
        int incr = 1;
        const TacsScalar *ni = NULL;
        if (ix % (num_params + 3) == 0){
          ni = N;
        }
        else if (ix % (num_params + 3) == 1){
          ni = N;
        }
        else if (ix % (num_params + 3) == 2){
          ni = N;
        }
        else {
          incr = num_params;
          ni = &Nx[(ix % (num_params + 3)) - 3];
        }

        for ( int i = 0; i < num_nodes; i++ ){
          TacsScalar scale = 0.0;
          int jncr = 1;
          const TacsScalar *nj = NULL;
          if (jx % (num_params + 3) == 0){
            scale = alpha*init_scale;
            nj = N;
          }
          else if (jx % (num_params + 3) == 1){
            scale = beta*init_scale;
            nj = N;
          }
          else if (jx % (num_params + 3) == 2){
            scale = gamma*init_scale;
            nj = N;
          }
          else {
            scale = alpha*init_scale;
            jncr = num_params;
            nj = &Nx[(jx % (num_params + 3)) - 3];
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
