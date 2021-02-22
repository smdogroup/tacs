#ifndef TACS_SHELL_ELEMENT_H
#define TACS_SHELL_ELEMENT_H

#include "TACSShellElementModel.h"
#include "TACSShellElementBasis.h"
#include "TACSGaussQuadrature.h"
#include "TACSElementAlgebra.h"
#include "TACSShellConstitutive.h"
#include "TACSElement.h"
#include "TACSElementTypes.h"

/*
  Compute the transformation from the local coordinates
  to
*/
class TACSShellTransform : public TACSObject {
 public:
  /*
    Given the local shell element reference frame Xf, compute the
    transformation from the global coordinates to the shell-aligned local axis.
  */
  virtual void computeTransform( const TacsScalar Xxi[], TacsScalar T[] ) = 0;
};

class TACSShellNaturalTransform : public TACSShellTransform {
 public:
  TACSShellNaturalTransform(){}

  void computeTransform( const TacsScalar Xxi[], TacsScalar T[] ){
    // Compute the transformation
    TacsScalar t1[3], t2[3];
    t1[0] = Xxi[0];
    t1[1] = Xxi[2];
    t1[2] = Xxi[4];

    t2[0] = Xxi[1];
    t2[1] = Xxi[3];
    t2[2] = Xxi[5];

    // Compute the normal direction
    TacsScalar n[3];
    crossProduct(t1, t2, n);

    // Normalize the normal direction
    TacsScalar invNorm = 1.0/sqrt(vec3Dot(n, n));
    vec3Scale(invNorm, n);

    // Normalize the 1-direction of the element
    TacsScalar inv = 1.0/sqrt(vec3Dot(t1, t1));
    vec3Scale(inv, t1);

    // Take the cross product to determine the 2-direction
    crossProduct(n, t1, t2);

    // Set the components of the transformation
    T[0] = t1[0];
    T[3] = t1[1];
    T[6] = t1[2];

    T[1] = t2[0];
    T[4] = t2[1];
    T[7] = t2[2];

    T[2] = n[0];
    T[5] = n[1];
    T[8] = n[2];
  }
};


class TACSShellRefAxisTransform : public TACSShellTransform {
 public:
  TACSShellRefAxisTransform( const TacsScalar _axis[] ){
    axis[0] = _axis[0];
    axis[1] = _axis[1];
    axis[2] = _axis[2];

    TacsScalar norm = sqrt(vec3Dot(axis, axis));
    TacsScalar invNorm = 0.0;
    if (norm != 0.0){
      invNorm = 1.0/norm;
    }
    vec3Scale(invNorm, axis);
  }

  void computeTransform( const TacsScalar Xxi[], TacsScalar T[] ){
    // Compute the transformation
    TacsScalar t1[3], t2[3];
    t1[0] = Xxi[0];
    t1[1] = Xxi[2];
    t1[2] = Xxi[4];

    t2[0] = Xxi[1];
    t2[1] = Xxi[3];
    t2[2] = Xxi[5];

    // Compute the normal direction
    TacsScalar n[3];
    crossProduct(t1, t2, n);

    // Normalize the normal direction
    TacsScalar invNorm = 1.0/sqrt(vec3Dot(n, n));
    vec3Scale(invNorm, n);

    // Compute the dot product with
    TacsScalar an = vec3Dot(axis, n);

    // Take the component of the reference axis perpendicular
    // to the surface
    t1[0] = axis[0] - an*n[0];
    t1[1] = axis[1] - an*n[1];
    t1[2] = axis[2] - an*n[2];

    // Normalize the new direction
    TacsScalar inv = 1.0/sqrt(vec3Dot(t1, t1));
    vec3Scale(inv, t1);

    // Take the cross product to determine the 2-direction
    crossProduct(n, t1, t2);

    // Set the components of the transformation
    T[0] = t1[0];
    T[3] = t1[1];
    T[6] = t1[2];

    T[1] = t2[0];
    T[4] = t2[1];
    T[7] = t2[2];

    T[2] = n[0];
    T[5] = n[1];
    T[8] = n[2];
  }

 private:
  TacsScalar axis[3];
};

/*
  The director class.

  Given a reference vector, t, from the element geometry, the director computes
  the exact or approximate rate of change of the displacement t.
*/
class TACSLinearizedRotation {
 public:

  static const int NUM_PARAMETERS = 3;

  /**
    Compute the director at a point.

    d = Q(q)*t = (C(q)^{T} - I)*t

    @param q The input rotation parametrization
    @param t The reference direction
    @param d The director values
  */
  static void computeDirector( const TacsScalar q[],
                               const TacsScalar t[],
                               TacsScalar d[] ){
    crossProduct(q, t, d);
  }

  /*
    Compute the director and rates at a point.

    d = Q(q)*t = (C(q)^{T} - I)*t
    ddot = d/dt(Q(q))*t
    dddot = d^2/dt^2(Q(q))*t

    @param q The input rotation parametrization
    @param t The reference direction
    @param d The director values
  */
  static void computeDirectorRates( const TacsScalar q[],
                                    const TacsScalar qdot[],
                                    const TacsScalar t[],
                                    TacsScalar d[],
                                    TacsScalar ddot[] ){
    crossProduct(q, t, d);
    crossProduct(qdot, t, ddot);
  }

  /*
    Compute the director and rates at a point.

    d = Q(q)*t = (C(q)^{T} - I)*t
    ddot = d/dt(Q(q))*t
    dddot = d^2/dt^2(Q(q))*t

    @param q The input rotation parametrization
    @param t The reference direction
    @param d The director values
  */
  static void computeDirectorRates( const TacsScalar q[],
                                    const TacsScalar qdot[],
                                    const TacsScalar qddot[],
                                    const TacsScalar t[],
                                    TacsScalar d[],
                                    TacsScalar ddot[],
                                    TacsScalar dddot[] ){
    crossProduct(q, t, d);
    crossProduct(qdot, t, ddot);
    crossProduct(qddot, t, dddot);
  }

  /*
    Given the derivatives of the kinetic energy expression with respect to time,
    add the contributions to the derivative of the

    Given the partial derivatives of the Lagrangian with respect to the
    director and the time derivative of the vector, compute

    ddtdTddot = d/dt(dT/d(ddot))
    dTddot = dT/d(dot)

    Compute:

    res += scale*(d/dt(dT/d(ddot))*d(ddot)/d(qdot) + dT/d(ddot)*d/dt(d(ddot)/d(qdot)))

  */
  static void addDirectorResidual( const TacsScalar q[],
                                    const TacsScalar qdot[],
                                    const TacsScalar qddot[],
                                    const TacsScalar t[],
                                    const TacsScalar dLdd[],
                                    TacsScalar res[] ){
    crossProductAdd(1.0, t, dLdd, res);
  }
};



/*
  The director class.

  Given a reference vector, t, from the element geometry, the director computes
  the exact or approximate rate of change of the displacement t.
*/
class TACSQuaternionDirector : public TACSObject {
 public:

  static const int NUM_PARAMETERS = 5;

  /**
    Compute the director at a point.

    d = Q(q)*t = (C(q)^{T} - I)*t

    @param q The input rotation parametrization
    @param t The reference direction
    @param d The director values
  */
  static void computeDirector( const TacsScalar q[],
                               const TacsScalar t[],
                               TacsScalar d[] ){
    // Compute Q = C^{T} - I
    TacsScalar Q[9];
    Q[0] =-2.0*(q[2]*q[2] + q[3]*q[3]);
    Q[1] = 2.0*(q[2]*q[1] - q[3]*q[0]);
    Q[2] = 2.0*(q[3]*q[1] + q[2]*q[0]);

    Q[3] = 2.0*(q[1]*q[2] + q[3]*q[0]);
    Q[4] =-2.0*(q[1]*q[1] + q[3]*q[3]);
    Q[5] = 2.0*(q[3]*q[2] - q[1]*q[0]);

    Q[6] = 2.0*(q[1]*q[3] - q[2]*q[0]);
    Q[7] = 2.0*(q[2]*q[3] + q[1]*q[0]);
    Q[8] =-2.0*(q[1]*q[1] + q[2]*q[2]);

    // Compute d = Q*t
    d[0] = Q[0]*t[0] + Q[1]*t[1] + Q[2]*t[2];
    d[1] = Q[3]*t[0] + Q[4]*t[1] + Q[5]*t[2];
    d[2] = Q[6]*t[0] + Q[7]*t[1] + Q[8]*t[2];
  }

  /*
    Compute the director and rates at a point.

    d = Q(q)*t = (C(q)^{T} - I)*t
    ddot = d/dt(Q(q))*t
    dddot = d^2/dt^2(Q(q))*t

    @param q The input rotation parametrization
    @param t The reference direction
    @param d The director values
  */
  // void computeDirectorRates( const TacsScalar q[],
  //                            const TacsScalar qdot[],
  //                            const TacsScalar qddot[],
  //                            const TacsScalar t[],
  //                            TacsScalar d[],
  //                            TacsScalar ddot[],
  //                            TacsScalar dddot[] ){
  //   oid TACSShellQuaternion::getAngularAcceleration( const int num_nodes,
  //                                                 const int vars_per_node,
  //                                                 const TacsScalar fn[],
  //                                                 const TacsScalar vars[],
  //                                                 const TacsScalar dvars[],
  //                                                 const TacsScalar ddvars[],
  //                                                 TacsScalar omega[],
  //                                                 TacsScalar domega[] )

  // for ( int i = 0; i < num_nodes; i++ ){
  //   TacsScalar eta = vars[0];
  //   const TacsScalar *eps = &vars[1];
  //   TacsScalar deta = dvars[0];
  //   const TacsScalar *deps = &dvars[1];
  //   TacsScalar ddeta = ddvars[0];
  //   const TacsScalar *ddeps = &ddvars[1];

  //   // omega = -2*eps^{x}*deps + 2*eta*deps - eps*deta
  //   TacsScalar omeg[3];
  //   crossProduct(-2.0, eps, deps, omega);
  //   vecAxpy(2.0*eta, deps, omega);
  //   vecAxpy(-2.0*deta, eps, omega);

  //   // domega = S(q)*ddot{q}
  //   TacsScalar domeg[3];
  //   crossProduct(-2.0, eps, ddeps, domeg);
  //   vecAxpy(2.0*eta, ddeps, domeg);
  //   vecAxpy(-2.0*ddeta, eps, domeg);

  //   TacsScalar tmp = 0.0;
  //   tmp = vecDot(omeg, fn);
  //   omega[0] = omeg[0] - tmp*fn[0];
  //   omega[1] = omeg[1] - tmp*fn[1];
  //   omega[2] = omeg[2] - tmp*fn[2];

  //   tmp = vecDot(domeg, fn);
  //   domega[0] = domeg[0] - tmp*fn[0];
  //   domega[1] = domeg[1] - tmp*fn[1];
  //   domega[2] = domeg[2] - tmp*fn[2];

  //   fn += 3;
  //   vars += vars_per_node;
  //   dvars += vars_per_node;
  //   ddvars += vars_per_node;
  //   omega += 3;
  //   domega += 3;
  // }

  // }

  /*
    Given the derivatives of the kinetic energy expression with respect to time,
    add the contributions to the derivative of the

    Given the partial derivatives of the Lagrangian with respect to the
    director and the time derivative of the vector, compute

    ddtdTddot = d/dt(dT/d(ddot))
    dTddot = dT/d(dot)

    Compute:

    res += scale*(d/dt(dT/d(ddot))*d(ddot)/d(qdot) + dT/d(ddot)*d/dt(d(ddot)/d(qdot)))

  */
//  void addResidual( TacsScalar scale,
//                    const TacsScalar ddtdTddot[],
//                    const TacsScalar dTddot[],
//                    const TacsScalar dLdd[],
//                    const TacsScalar t[],
//                    TacsScalar res[] ){
//     TacsScalar ddot


//     res[0] +=
//   }
};



template <class quadrature, class basis, class director, class model>
class TACSShellElement : public TACSElement {
 public:
  TACSShellElement( TACSShellTransform *_transform,
                    TACSShellConstitutive *_con ){
    transform = _transform;
    transform->incref();

    con = _con;
    con->incref();
  }

  int getVarsPerNode(){
    return 3 + director::NUM_PARAMETERS;
  }
  int getNumNodes(){
    return basis::NUM_NODES;
  }

  ElementLayout getLayoutType(){
    return basis::getLayoutType();
  }


  void computeEnergies( int elemIndex,
                        double time,
                        const TacsScalar Xpts[],
                        const TacsScalar vars[],
                        const TacsScalar dvars[],
                        TacsScalar *Te,
                        TacsScalar *Pe );

  void addResidual( int elemIndex,
                    double time,
                    const TacsScalar *Xpts,
                    const TacsScalar *vars,
                    const TacsScalar *dvars,
                    const TacsScalar *ddvars,
                    TacsScalar *res );

/*
  void addJacobian( int elemIndex, double time,
                    TacsScalar alpha,
                    TacsScalar beta,
                    TacsScalar gamma,
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[],
                    TacsScalar res[],
                    TacsScalar mat[] );
*/

  void getOutputData( int elemIndex,
                      ElementType etype,
                      int write_flag,
                      const TacsScalar Xpts[],
                      const TacsScalar vars[],
                      const TacsScalar dvars[],
                      const TacsScalar ddvars[],
                      int ld_data,
                      TacsScalar *data );

 private:
  TACSShellTransform *transform;
  TACSShellConstitutive *con;

  void getNodeNormals( const TacsScalar Xpts[],
                       TacsScalar fn[],
                       TacsScalar fnorm[]=NULL );

  static void assembleFrame( const TacsScalar Xxi[],
                             const TacsScalar n[],
                             TacsScalar Xd[] ){
    Xd[0] = Xxi[0];
    Xd[1] = Xxi[1];
    Xd[2] = n[0];

    Xd[3] = Xxi[2];
    Xd[4] = Xxi[3];
    Xd[5] = n[1];

    Xd[6] = Xxi[4];
    Xd[7] = Xxi[5];
    Xd[8] = n[2];
  }

  static void assembleFrame( const TacsScalar nxi[],
                             TacsScalar Xdz[] ){
    Xdz[0] = nxi[0];
    Xdz[1] = nxi[1];
    Xdz[2] = 0.0;

    Xdz[3] = nxi[2];
    Xdz[4] = nxi[3];
    Xdz[5] = 0.0;

    Xdz[6] = nxi[4];
    Xdz[7] = nxi[5];
    Xdz[8] = 0.0;
  }

  static void extractFrame( const TacsScalar Xd[],
                            TacsScalar Xxi[],
                            TacsScalar n[] ){
    Xxi[0] = Xd[0];
    Xxi[1] = Xd[1];
    n[0] = Xd[2];

    Xxi[2] = Xd[3];
    Xxi[3] = Xd[4];
    n[1] = Xd[5];

    Xxi[4] = Xd[6];
    Xxi[5] = Xd[7];
    n[2] = Xd[8];
  }

  static void extractFrame( const TacsScalar Xd[],
                            TacsScalar Xxi[] ){
    Xxi[0] = Xd[0];
    Xxi[1] = Xd[1];

    Xxi[2] = Xd[3];
    Xxi[3] = Xd[4];

    Xxi[4] = Xd[6];
    Xxi[5] = Xd[7];
  }
};

template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
  getNodeNormals( const TacsScalar Xpts[],
                  TacsScalar fn[],
                  TacsScalar fnorm[] ){
  for ( int i = 0; i < basis::NUM_NODES; i++ ){
    double pt[2];
    basis::getNodePoint(i, pt);

    // Compute the derivative X,xi at each node
    TacsScalar Xxi[6];
    basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);

    TacsScalar a[3], b[3];
    a[0] = Xxi[0];
    a[1] = Xxi[2];
    a[2] = Xxi[4];

    b[0] = Xxi[1];
    b[1] = Xxi[3];
    b[2] = Xxi[5];

    // Compute the normal direction at the point
    crossProduct(a, b, &fn[3*i]);

    // Compute the 2-norm of the vector in the normal direction
    TacsScalar norm = sqrt(vec3Dot(&fn[3*i], &fn[3*i]));

    // Save the 2-norm value if the fnorm argument is not NULL
    if (fnorm){
      fnorm[i] = norm;
    }

    // Scale the normal direction
    if (norm != 0.0){
      vec3Scale(1.0/norm, &fn[3*i]);
    }
  }
}

/*
  Add the residual to the provided vector
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
  computeEnergies( int elemIndex,
                   double time,
                   const TacsScalar *Xpts,
                   const TacsScalar *vars,
                   const TacsScalar *dvars,
                   TacsScalar *_Te, TacsScalar *_Ue ){
  // Zero the kinetic and potential energies
  TacsScalar Te = 0.0;
  TacsScalar Ue = 0.0;

  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();
  const int vars_per_node = 3 + director::NUM_PARAMETERS;

  // Compute the node normal directions
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals(Xpts, fn);

  // Compute the frame normal and directors at each node
  TacsScalar d[3*basis::NUM_NODES];
  TacsScalar ddot[3*basis::NUM_NODES];
  for ( int i = 0, offset = 3; i < basis::NUM_NODES; i++, offset += vars_per_node ){
    director::computeDirectorRates(&vars[offset], &dvars[offset],
                                   &fn[3*i], &d[3*i], &ddot[3*i]);
  }

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<basis>(Xpts, fn, vars_per_node, vars, d, ety);

  // Loop over each quadrature point and add the residual contribution
  for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // Compute X,xi = [dX/dxi1 ; dX/dxi2] and n,xi = [dn/dxi1; dn/dxi2]
    TacsScalar X[3], Xxi[6], n[3], nxi[6];
    basis::interpFields(pt, 3, Xpts, 3, X);
    basis::interpFields(pt, 3, fn, 3, n);
    basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);
    basis::interpFieldsGrad(pt, 3, fn, 3, nxi);

    // Compute the transformation at the quadrature point
    TacsScalar T[9];
    transform->computeTransform(Xxi, T);

    // Assemble the terms Xd = [Xxi; n] and Xdz
    TacsScalar Xd[9], Xdz[9];
    assembleFrame(Xxi, n, Xd);
    assembleFrame(nxi, Xdz);

    // Compute the inverse of the 3x3 Jacobian transformation
    TacsScalar Xdinv[9];
    TacsScalar detXd = inv3x3(Xd, Xdinv);
    detXd *= weight;

    // Compute negXdinvXdz = -Xdinv*Xdz
    TacsScalar negXdinvXdz[9];
    mat3x3MatMult(Xdinv, Xdz, negXdinvXdz);
    for ( int i = 0; i < 9; i++ ){
      negXdinvXdz[i] *= -1.0;
    }

    // Compute XdinvT = Xdinv*T
    TacsScalar XdinvT[9];
    mat3x3MatMult(Xdinv, T, XdinvT);

    // Compute the director field and the gradient of the director
    // field at the specified point
    TacsScalar d0[3], d0xi[6];
    basis::interpFields(pt, 3, d, 3, d0);
    basis::interpFieldsGrad(pt, 3, d, 3, d0xi);

    // Compute the gradient of the displacement solution at the quadrature points
    TacsScalar u0xi[6];
    basis::interpFieldsGrad(pt, vars_per_node, vars, 3, u0xi);

    // Input: u0xi, d0, d0xi, T, negXdinvXdz, XdinvT
    // Output: u0x, u1x

    // Set u0x = [u0,xi ; d]
    TacsScalar u0x[9];
    assembleFrame(u0xi, d0, u0x); // Use u0x to store [u0,xi; d0]

    // u1x = T^{T}*(u0d*(-Xdinv*Xdz) + u1d)*Xdinv*T
    TacsScalar u1x[9], tmp[9];
    assembleFrame(d0xi, u1x); // Use u1x to store [d0,xi; 0]
    mat3x3MatMultAdd(u0x, negXdinvXdz, u1x);
    mat3x3TransMatMult(T, u1x, tmp);
    mat3x3MatMult(tmp, XdinvT, u1x);

    // Compute the transformation u0x = T^{T}*ueta*Xdinv*T
    // u0x = T^{T}*u0d*Xdinv*T
    mat3x3MatMult(u0x, XdinvT, tmp);
    mat3x3TransMatMult(T, tmp, u0x);

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    model::template interpTyingStrain<basis>(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, e);

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, e, s);

    Ue += 0.5*detXd*(s[0]*e[0] + s[1]*e[1] + s[2]*e[2] +
                     s[3]*e[3] + s[4]*e[4] + s[5]*e[5] +
                     s[6]*e[6] + s[7]*e[7] + s[8]*e[8]);
  }

  for ( int i = 0; i < basis::NUM_NODES; i++ ){
    Ue += 0.5*vars[5 + 6*i]*vars[5 + 6*i];
  }

  *_Te = Te;
  *_Ue = Ue;
}

/*
  Add the residual to the provided vector
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
  addResidual( int elemIndex,
               double time,
               const TacsScalar *Xpts,
               const TacsScalar *vars,
               const TacsScalar *dvars,
               const TacsScalar *ddvars,
               TacsScalar *res ){
  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();
  const int vars_per_node = 3 + director::NUM_PARAMETERS;

  // Compute the node normal directions
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals(Xpts, fn);

  // Compute the frame normal and directors at each node
  TacsScalar d[3*basis::NUM_NODES];
  TacsScalar ddot[3*basis::NUM_NODES];
  TacsScalar dddot[3*basis::NUM_NODES];
  for ( int i = 0, offset = 3; i < basis::NUM_NODES; i++, offset += vars_per_node ){
    director::computeDirectorRates(&vars[offset], &dvars[offset], &ddvars[offset],
                                   &fn[3*i], &d[3*i], &ddot[3*i], &dddot[3*i]);
  }

  // Derivative of the director field
  TacsScalar dd[3*basis::NUM_NODES];
  memset(dd, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<basis>(Xpts, fn, vars_per_node, vars, d, ety);

  // Zero the contributions to the
  TacsScalar dety[basis::NUM_TYING_POINTS];
  memset(dety, 0, basis::NUM_TYING_POINTS*sizeof(TacsScalar));

  // Loop over each quadrature point and add the residual contribution
  for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // Compute X,xi = [dX/dxi1 ; dX/dxi2] and n,xi = [dn/dxi1; dn/dxi2]
    TacsScalar X[3], Xxi[6], n[3], nxi[6];
    basis::interpFields(pt, 3, Xpts, 3, X);
    basis::interpFields(pt, 3, fn, 3, n);
    basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);
    basis::interpFieldsGrad(pt, 3, fn, 3, nxi);

    // Compute the transformation at the quadrature point
    TacsScalar T[9];
    transform->computeTransform(Xxi, T);

    // Assemble the terms Xd = [Xxi; n] and Xdz
    TacsScalar Xd[9], Xdz[9];
    assembleFrame(Xxi, n, Xd);
    assembleFrame(nxi, Xdz);

    // Compute the inverse of the 3x3 Jacobian transformation
    TacsScalar Xdinv[9];
    TacsScalar detXd = inv3x3(Xd, Xdinv);
    detXd *= weight;

    // Compute negXdinvXdz = -Xdinv*Xdz
    TacsScalar negXdinvXdz[9];
    mat3x3MatMult(Xdinv, Xdz, negXdinvXdz);
    for ( int i = 0; i < 9; i++ ){
      negXdinvXdz[i] *= -1.0;
    }

    // Compute XdinvT = Xdinv*T
    TacsScalar XdinvT[9];
    mat3x3MatMult(Xdinv, T, XdinvT);

    // Compute the director field and the gradient of the director
    // field at the specified point
    TacsScalar d0[3], d0xi[6];
    basis::interpFields(pt, 3, d, 3, d0);
    basis::interpFieldsGrad(pt, 3, d, 3, d0xi);

    // Compute the gradient of the displacement solution at the quadrature points
    TacsScalar u0xi[6];
    basis::interpFieldsGrad(pt, vars_per_node, vars, 3, u0xi);

    // Input: u0xi, d0, d0xi, T, negXdinvXdz, XdinvT
    // Output: u0x, u1x

    // Set u0x = [u0,xi ; d]
    TacsScalar u0x[9];
    assembleFrame(u0xi, d0, u0x); // Use u0x to store [u0,xi; d0]

    // u1x = T^{T}*(u0d*(-Xdinv*Xdz) + u1d)*Xdinv*T
    TacsScalar u1x[9], tmp[9];
    assembleFrame(d0xi, u1x); // Use u1x to store [d0,xi; 0]
    mat3x3MatMultAdd(u0x, negXdinvXdz, u1x);
    mat3x3TransMatMult(T, u1x, tmp);
    mat3x3MatMult(tmp, XdinvT, u1x);

    // Compute the transformation u0x = T^{T}*ueta*Xdinv*T
    // u0x = T^{T}*u0d*Xdinv*T
    mat3x3MatMult(u0x, XdinvT, tmp);
    mat3x3TransMatMult(T, tmp, u0x);

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    model::template interpTyingStrain<basis>(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, e);

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, e, s);

    // Compute the derivative of the product of the stress and strain
    // with respect to u0x, u1x and e0ty
    TacsScalar du0x[9], du1x[9], de0ty[6];
    model::evalStrainSens(detXd, s, u0x, u1x, du0x, du1x, de0ty);

    // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
    TacsScalar dgty[6];
    mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);

    // Evaluate the tying strain
    model::template addInterpTyingStrainTranspose<basis>(pt, dgty, dety);

    // Compute du0d = T*(du0x*XdinvT^{T} + du1x*XdinvT^{T}*negXdinvXdz^{T})
    TacsScalar du0d[9];
    mat3x3MatTransMult(du1x, XdinvT, du0d);
    mat3x3MatTransMult(du0d, negXdinvXdz, tmp);
    mat3x3MatTransMultAdd(du0x, XdinvT, tmp);
    mat3x3MatMult(T, tmp, du0d);

    // Compute du1d = T*du1x*XdinvT^{T}
    TacsScalar du1d[9];
    mat3x3MatTransMult(du1x, XdinvT, tmp);
    mat3x3MatMult(T, tmp, du1d);

    // du0d = [du0xi; dd0]
    TacsScalar du0xi[6], dd0[3];
    extractFrame(du0d, du0xi, dd0);

    TacsScalar dd0xi[6];
    extractFrame(du1d, dd0xi);

    // Compute the director field and the gradient of the director
    // field at the specified point
    basis::addInterpFieldsTranspose(pt, 3, dd0, 3, dd);
    basis::addInterpFieldsGradTranspose(pt, 3, dd0xi, 3, dd);

    // Compute the gradient of the displacement solution at the quadrature points
    basis::addInterpFieldsGradTranspose(pt, 3, du0xi, vars_per_node, res);
  }

  // Set the total number of tying points needed for this element
  model::template addComputeTyingStrainTranspose<basis>(Xpts, fn, vars_per_node, vars, d, dety, res, dd);

  // Add the contributions to the director field
  for ( int i = 0, offset = 3; i < basis::NUM_NODES; i++, offset += vars_per_node ){
    director::addDirectorResidual(&vars[offset], &dvars[offset], &ddvars[offset],
                                   &fn[3*i], &dd[3*i], &res[offset]);
  }

  for ( int i = 0; i < basis::NUM_NODES; i++ ){
    res[5 + 6*i] += vars[5 + 6*i];
  }
}

/*
  Add the residual to the provided vector
*/
/*
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
  addJacobian( int elemIndex, double time,
               TacsScalar alpha,
               TacsScalar beta,
               TacsScalar gamma,
               const TacsScalar Xpts[],
               const TacsScalar vars[],
               const TacsScalar dvars[],
               const TacsScalar ddvars[],
               TacsScalar res[],
               TacsScalar mat[] ){

}
*/

/*
  Get the element data for the basis
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
  getOutputData( int elemIndex,
                 ElementType etype,
                 int write_flag,
                 const TacsScalar Xpts[],
                 const TacsScalar vars[],
                 const TacsScalar dvars[],
                 const TacsScalar ddvars[],
                 int ld_data,
                 TacsScalar *data ){
  // Get the number of nodes associated with the visualization
  int num_vis_nodes = TacsGetNumVisNodes(basis::getLayoutType());

  // Compute the number of quadrature points
  const int vars_per_node = 3 + director::NUM_PARAMETERS;

  // Compute the node normal directions
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals(Xpts, fn);

  // Compute the frame normal and directors at each node
  TacsScalar d[3*basis::NUM_NODES];
  TacsScalar ddot[3*basis::NUM_NODES];
  for ( int i = 0, offset = 3; i < basis::NUM_NODES; i++, offset += vars_per_node ){
    director::computeDirectorRates(&vars[offset], &dvars[offset],
                                   &fn[3*i], &d[3*i], &ddot[3*i]);
  }

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<basis>(Xpts, fn, vars_per_node, vars, d, ety);

  // Loop over each quadrature point and add the residual contribution
  for ( int index = 0; index < num_vis_nodes; index++ ){
    // Get the quadrature weight
    double pt[3];
    basis::getNodePoint(index, pt);

    // Compute X,xi = [dX/dxi1 ; dX/dxi2] and n,xi = [dn/dxi1; dn/dxi2]
    TacsScalar X[3], Xxi[6], n[3], nxi[6];
    basis::interpFields(pt, 3, Xpts, 3, X);
    basis::interpFields(pt, 3, fn, 3, n);
    basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);
    basis::interpFieldsGrad(pt, 3, fn, 3, nxi);

    // Compute the transformation at the quadrature point
    TacsScalar T[9];
    transform->computeTransform(Xxi, T);

    // Assemble the terms Xd = [Xxi; n] and Xdz
    TacsScalar Xd[9], Xdz[9];
    assembleFrame(Xxi, n, Xd);
    assembleFrame(nxi, Xdz);

    // Compute the inverse of the 3x3 Jacobian transformation
    TacsScalar Xdinv[9];
    inv3x3(Xd, Xdinv);

    // Compute negXdinvXdz = -Xdinv*Xdz
    TacsScalar negXdinvXdz[9];
    mat3x3MatMult(Xdinv, Xdz, negXdinvXdz);
    for ( int i = 0; i < 9; i++ ){
      negXdinvXdz[i] *= -1.0;
    }

    // Compute XdinvT = Xdinv*T
    TacsScalar XdinvT[9];
    mat3x3MatTransMult(Xdinv, T, XdinvT);

    // Compute the director field and the gradient of the director
    // field at the specified point
    TacsScalar d0[3], d0xi[6];
    basis::interpFields(pt, 3, d, 3, d0);
    basis::interpFieldsGrad(pt, 3, d, 3, d0xi);

    // Compute the gradient of the displacement solution at the quadrature points
    TacsScalar u0xi[6];
    basis::interpFieldsGrad(pt, vars_per_node, vars, 3, u0xi);

    // Input: u0xi, d0, d0xi, T, negXdinvXdz, XdinvT
    // Output: u0x, u1x

    // Set u0x = [u0,xi ; d]
    TacsScalar u0x[9];
    assembleFrame(u0xi, d0, u0x); // Use u0x to store [u0,xi; d0]

    // u1x = T^{T}*(u0d*(-Xdinv*Xdz) + u1d)*Xdinv*T
    TacsScalar u1x[9], tmp[9];
    assembleFrame(d0xi, u1x); // Use u1x to store [d0,xi; 0]
    mat3x3MatMultAdd(u0x, negXdinvXdz, u1x);
    mat3x3TransMatMult(T, u1x, tmp);
    mat3x3MatMult(tmp, XdinvT, u1x);

    // Compute the transformation u0x = T^{T}*ueta*Xdinv*T
    // u0x = T^{T}*u0d*Xdinv*T
    mat3x3MatMult(u0x, XdinvT, tmp);
    mat3x3TransMatMult(T, tmp, u0x);

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    model::template interpTyingStrain<basis>(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, e);

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, e, s);

    if (etype == TACS_BEAM_OR_SHELL_ELEMENT){
      if (write_flag & TACS_OUTPUT_NODES){
        data[0] = X[0];
        data[1] = X[1];
        data[2] = X[2];
        data += 3;
      }
      if (write_flag & TACS_OUTPUT_DISPLACEMENTS){
        int len = vars_per_node;
        if (len > 6){
          len = 6;
        }
        for ( int i = 0; i < len; i++ ){
          data[i] = vars[i + vars_per_node*index];
        }
        for ( int i = len; i < 6; i++ ){
          data[i] = 0.0;
        }
        data += 6;
      }
      if (write_flag & TACS_OUTPUT_STRAINS){
        for ( int i = 0; i < 9; i++ ){
          data[i] = e[i];
        }
        data += 9;
      }
      if (write_flag & TACS_OUTPUT_STRESSES){
        for ( int i = 0; i < 9; i++ ){
          data[i] = s[i];
        }
        data += 9;
      }
      if (write_flag & TACS_OUTPUT_EXTRAS){
        data[0] = con->evalFailure(elemIndex, pt, X, e);
        data[1] = con->evalDesignFieldValue(elemIndex, pt, X, 0);
        data[2] = con->evalDesignFieldValue(elemIndex, pt, X, 1);
        data[3] = con->evalDesignFieldValue(elemIndex, pt, X, 2);
        data += 4;
      }
    }
  }
}

#endif // TACS_SHELL_ELEMENT_H
