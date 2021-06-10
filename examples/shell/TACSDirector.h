#ifndef TACS_DIRECTOR_H
#define TACS_DIRECTOR_H

#include "TACSObject.h"

/*
  The director class.

  Given a reference vector, t, from the element geometry, the director computes
  the exact or approximate rate of change of the displacement t.
*/
class TACSLinearizedRotation {
 public:
  static const int NUM_PARAMETERS = 3;

  /**
    Compute the rotation matrices at each node

    @param vars The full variable vector
    @param C The rotation matrices at each point
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeRotationMat( const TacsScalar vars[],
                                  TacsScalar C[] ){
    const TacsScalar *q = &vars[offset];
    for ( int i = 0; i < num_nodes; i++ ){
      // C = I - q^{x}
      C[0] = C[4] = C[8] = 1.0;
      setMatSkew(-1.0, q, C);

      C += 9;
      q += vars_per_node;
    }
  }

  /*
    Compute the derivative of the rotation matrices at each node

  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeRotationMatDeriv( const TacsScalar vars[],
                                       const TacsScalar varsd[],
                                       TacsScalar C[],
                                       TacsScalar Cd[] ){
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qd = &varsd[offset];
    for ( int i = 0; i < num_nodes; i++ ){
      // C = I - q^{x}
      C[0] = C[4] = C[8] = 1.0;
      setMatSkew(-1.0, q, C);

      // Cd = - qd^{x}
      Cd[0] = Cd[4] = Cd[8] = 0.0;
      setMatSkew(-1.0, qd, Cd);

      C += 9;
      Cd += 9;

      q += vars_per_node;
      qd += vars_per_node;
    }

  }

  /*
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationMatResidual( const TacsScalar vars[],
                                      const TacsScalar dC[],
                                      TacsScalar res[] ){
    TacsScalar *r = &res[offset];

    for ( int i = 0; i < num_nodes; i++ ){
      r[0] += -(dC[7] - dC[5]);
      r[1] += -(dC[2] - dC[6]);
      r[2] += -(dC[3] - dC[1]);

      r += vars_per_node;
      dC += 9;
    }
  }

  /*
    Add the derivative w.r.t.
  */


  /**
    Compute the director and rates at all nodes.

    d = Q(q)*t = (C(q)^{T} - I)*t
    ddot = d/dt(Q(q))*t

    @param vars The full variable vector
    @param dvars The first time derivative of the variables
    @param ddvars The second derivatives of the variables
    @param t The reference directions
    @param d The director values
    @param ddot The first time derivative of the director
    @param dddot The second time derivative of the director
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeDirectorRates( const TacsScalar vars[],
                                    const TacsScalar dvars[],
                                    const TacsScalar t[],
                                    TacsScalar d[],
                                    TacsScalar ddot[] ){
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    for ( int i = 0; i < num_nodes; i++ ){
      crossProduct(q, t, d);
      crossProduct(qdot, t, ddot);

      t += 3;
      d += 3;
      ddot += 3;

      q += vars_per_node;
      qdot += vars_per_node;
    }
  }

  /**
    Compute the director and rates at all nodes.

    d = Q(q)*t = (C(q)^{T} - I)*t
    ddot = d/dt(Q(q))*t
    dddot = d^2/dt^2(Q(q))*t

    @param vars The full variable vector
    @param dvars The first time derivative of the variables
    @param ddvars The second derivatives of the variables
    @param t The reference directions
    @param d The director values
    @param ddot The first time derivative of the director
    @param dddot The second time derivative of the director
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeDirectorRates( const TacsScalar vars[],
                                    const TacsScalar dvars[],
                                    const TacsScalar ddvars[],
                                    const TacsScalar t[],
                                    TacsScalar d[],
                                    TacsScalar ddot[],
                                    TacsScalar dddot[] ){
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    const TacsScalar *qddot = &ddvars[offset];
    for ( int i = 0; i < num_nodes; i++ ){
      crossProduct(q, t, d);
      crossProduct(qdot, t, ddot);
      crossProduct(qddot, t, dddot);

      t += 3;
      d += 3;
      ddot += 3;
      dddot += 3;

      q += vars_per_node;
      qdot += vars_per_node;
      qddot += vars_per_node;
    }
  }

  /**
    Compute the director and rates at all nodes and the derivative.

    d = Q(q)*t = (C(q)^{T} - I)*t
    ddot = d/dt(Q(q))*t
    dddot = d^2/dt^2(Q(q))*t

    @param vars The full variable vector
    @param dvars The first time derivative of the variables
    @param ddvars The second derivatives of the variables
    @param varsd The full variable vector derivative
    @param t The reference directions
    @param C The rotation matrices at each point
    @param d The director values
    @param ddot The first time derivative of the director
    @param dddot The second time derivative of the director
    @param Cd The derivative of the rotation matrices at each point
    @param dd The derivator of the director values
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeDirectorRatesDeriv( const TacsScalar vars[],
                                         const TacsScalar dvars[],
                                         const TacsScalar ddvars[],
                                         const TacsScalar varsd[],
                                         const TacsScalar t[],
                                         TacsScalar d[],
                                         TacsScalar ddot[],
                                         TacsScalar dddot[],
                                         TacsScalar dd[] ){
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    const TacsScalar *qddot = &ddvars[offset];
    const TacsScalar *qd = &varsd[offset];
    for ( int i = 0; i < num_nodes; i++ ){
      crossProduct(q, t, d);
      crossProduct(qdot, t, ddot);
      crossProduct(qddot, t, dddot);

      // Cd = - qd^{x}
      crossProduct(qd, t, dd);

      t += 3;
      d += 3;
      dd += 3;
      ddot += 3;
      dddot += 3;

      q += vars_per_node;
      qdot += vars_per_node;
      qddot += vars_per_node;
      qd += vars_per_node;
    }
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
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorResidual( const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar ddvars[],
                                   const TacsScalar t[],
                                   const TacsScalar dd[],
                                   TacsScalar res[] ){
    TacsScalar *r = &res[offset];

    for ( int i = 0; i < num_nodes; i++ ){
      crossProductAdd(1.0, t, dd, r);

      r += vars_per_node;
      dd += 3;
      t += 3;
    }
  }

  /*
    Add terms from the Jacobian
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorJacobian( const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar ddvars[],
                                   const TacsScalar t[],
                                   const TacsScalar d2d[],
                                   const TacsScalar d2du[],
                                   TacsScalar mat[] ){
    // Add the derivative due to and d2d
    const int dsize = 3*num_nodes;
    const int nvars = vars_per_node*num_nodes;

    // d = crossProduct(q, t, d)
    const TacsScalar *ti = t;
    for ( int i = 0; i < num_nodes; i++, ti += 3 ){
      TacsScalar *jac1 = &mat[(offset + vars_per_node*i)*nvars + 3];
      TacsScalar *jac2 = &mat[(offset + vars_per_node*i + 1)*nvars + 3];
      TacsScalar *jac3 = &mat[(offset + vars_per_node*i + 2)*nvars + 3];

      const TacsScalar *tj = t;
      for ( int j = 0; j < num_nodes; j++, tj += 3 ){
        // Add the derivative
        TacsScalar d[9];
        d[0] = d2d[0];
        d[1] = d2d[1];
        d[2] = d2d[2];

        d[3] = d2d[dsize];
        d[4] = d2d[dsize+1];
        d[5] = d2d[dsize+2];

        d[6] = d2d[2*dsize];
        d[7] = d2d[2*dsize+1];
        d[8] = d2d[2*dsize+2];

        TacsScalar tmp[9];
        mat3x3SkewMatSkewTransform(ti, d, tj, tmp);

        jac1[0] -= tmp[0];
        jac1[1] -= tmp[1];
        jac1[2] -= tmp[2];

        jac2[0] -= tmp[3];
        jac2[1] -= tmp[4];
        jac2[2] -= tmp[5];

        jac3[0] -= tmp[6];
        jac3[1] -= tmp[7];
        jac3[2] -= tmp[8];

        jac1 += vars_per_node;
        jac2 += vars_per_node;
        jac3 += vars_per_node;
        d2d += 3;
      }

      d2d += 2*dsize;
    }

    for ( int i = 0; i < num_nodes; i++ ){
      for ( int j = 0; j < num_nodes; j++ ){
        // Add the derivative
        TacsScalar d[9];
        d[0] = d2du[0];
        d[1] = d2du[1];
        d[2] = d2du[2];

        d[3] = d2du[dsize];
        d[4] = d2du[dsize+1];
        d[5] = d2du[dsize+2];

        d[6] = d2du[2*dsize];
        d[7] = d2du[2*dsize+1];
        d[8] = d2du[2*dsize+2];

        TacsScalar tmp[9];
        mat3x3SkewMatTransform(&t[3*i], d, tmp);

        for ( int ii = 0; ii < 3; ii++ ){
          for ( int jj = 0; jj < 3; jj++ ){
            int index =
              (vars_per_node*i + ii + offset)*nvars +
              vars_per_node*j + jj;

            mat[index] += tmp[3*ii + jj];
          }
        }

        for ( int ii = 0; ii < 3; ii++ ){
          for ( int jj = 0; jj < 3; jj++ ){
            int index =
              (vars_per_node*j + jj)*nvars +
              vars_per_node*i + ii + offset;

            mat[index] += tmp[3*ii + jj];
          }
        }

        d2du += 3;
      }

      d2du += 2*dsize;
    }
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

#endif // TACS_DIRECTOR_H
