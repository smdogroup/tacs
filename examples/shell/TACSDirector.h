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
                                    TacsScalar C[],
                                    TacsScalar d[],
                                    TacsScalar ddot[] ){
    // C = I - q^{x}
    C[0] = C[4] = C[8] = 1.0;
    setMatSkew(-1.0, q, C);
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
                                    TacsScalar C[],
                                    TacsScalar d[],
                                    TacsScalar ddot[],
                                    TacsScalar dddot[] ){
    // C = I - q^{x}
    C[0] = C[4] = C[8] = 1.0;
    setMatSkew(-1.0, q, C);
    crossProduct(q, t, d);
    crossProduct(qdot, t, ddot);
    crossProduct(qddot, t, dddot);
  }

  /*
    Compute the directional derivative of the director
  */
  static void computeDirectorRatesDeriv( const TacsScalar q[],
                                         const TacsScalar qdot[],
                                         const TacsScalar qddot[],
                                         const TacsScalar qd[],
                                         const TacsScalar t[],
                                         TacsScalar C[],
                                         TacsScalar d[],
                                         TacsScalar ddot[],
                                         TacsScalar dddot[],
                                         TacsScalar Cd[],
                                         TacsScalar dd[] ){
    // C = I - q^{x}
    C[0] = C[4] = C[8] = 1.0;
    setMatSkew(-1.0, q, C);
    crossProduct(q, t, d);
    crossProduct(qdot, t, ddot);
    crossProduct(qddot, t, dddot);

    // Cd = - qd^{x}
    Cd[0] = Cd[4] = Cd[8] = 0.0;
    setMatSkew(-1.0, qd, Cd);
    crossProduct(qd, t, dd);
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
                                    const TacsScalar dC[],
                                    const TacsScalar dd[],
                                    TacsScalar res[] ){
    res[0] += -(dC[7] - dC[5]);
    res[1] += -(dC[2] - dC[6]);
    res[2] += -(dC[3] - dC[1]);
    crossProductAdd(1.0, t, dd, res);
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