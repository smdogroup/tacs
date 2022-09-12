#ifndef TACS_DIRECTOR_H
#define TACS_DIRECTOR_H

#include "TACSElementAlgebra.h"
#include "TACSElementVerification.h"

/*
  The director class.

  Given a reference vector, t, from the element geometry, the director
  computes the exact or approximate rate of change of the displacement
  through the thickness.
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
  static void computeRotationMat(const TacsScalar vars[], TacsScalar C[]) {
    const TacsScalar *q = &vars[offset];
    for (int i = 0; i < num_nodes; i++) {
      // C = I - q^{x}
      setMatSkew(-1.0, q, C);
      C[0] = C[4] = C[8] = 1.0;

      C += 9;
      q += vars_per_node;
    }
  }

  /**
    Compute the derivative of the rotation matrices at each node

    @param vars The full variable vector
    @param varsd The full variable vector
    @param C The rotation matrices at each point
    @param Cd The rotation matrices at each point
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeRotationMatDeriv(const TacsScalar vars[],
                                      const TacsScalar varsd[], TacsScalar C[],
                                      TacsScalar Cd[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qd = &varsd[offset];
    for (int i = 0; i < num_nodes; i++) {
      // C = I - q^{x}
      setMatSkew(-1.0, q, C);
      C[0] = C[4] = C[8] = 1.0;

      // Cd = - qd^{x}
      setMatSkew(-1.0, qd, Cd);

      C += 9;
      Cd += 9;

      q += vars_per_node;
      qd += vars_per_node;
    }
  }

  /**
    Add the contribution to the residual from the rotation matrix

    This code adds the contribution to the residual via the derivative

    d(tr(dC^{T}C(q)))/dq_{i} = d(tr(dC^{T}*(I - q^{x})))/dq_{i}

    @param vars The full variable vector
    @param dC The derivative w.r.t. the rotation matrix
    @param res The residual array
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationMatResidual(const TacsScalar vars[],
                                     const TacsScalar dC[], TacsScalar res[]) {
    TacsScalar *r = &res[offset];

    for (int i = 0; i < num_nodes; i++) {
      r[0] += -(dC[7] - dC[5]);
      r[1] += -(dC[2] - dC[6]);
      r[2] += -(dC[3] - dC[1]);

      r += vars_per_node;
      dC += 9;
    }
  }

  /*
    Add the Jacobian of the rotation matrix to the output

    @param alpha Scalar coefficient for the Jacobian matrix
    @param vars The variable values
    @param dC The derivative of the functional w.r.t. C
    @param d2C The second derivatives of the functional w.r.t. C
    @param res The residual
    @param mat The Jacobian matrix
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationMatJacobian(const TacsScalar alpha,
                                     const TacsScalar vars[],
                                     const TacsScalar dC[],
                                     const TacsScalar d2C[], TacsScalar res[],
                                     TacsScalar mat[]) {
    const int size = vars_per_node * num_nodes;
    const int csize = 9 * num_nodes;

    TacsScalar *r = NULL;
    if (res) {
      r = &res[offset];
    }
    TacsScalar *m = &mat[offset * size + offset];

    for (int i = 0; i < num_nodes; i++) {
      if (res) {
        r[0] += -(dC[7] - dC[5]);
        r[1] += -(dC[2] - dC[6]);
        r[2] += -(dC[3] - dC[1]);
        r += vars_per_node;
      }

      for (int j = 0; j < num_nodes; j++) {
        m[vars_per_node * j] += d2C[csize * (9 * i + 5) + 9 * j + 5] -
                                d2C[csize * (9 * i + 5) + 9 * j + 7] -
                                d2C[csize * (9 * i + 7) + 9 * j + 5] +
                                d2C[csize * (9 * i + 7) + 9 * j + 7];

        m[vars_per_node * j + 1] += -d2C[csize * (9 * i + 5) + 9 * j + 2] +
                                    d2C[csize * (9 * i + 5) + 9 * j + 6] +
                                    d2C[csize * (9 * i + 7) + 9 * j + 2] -
                                    d2C[csize * (9 * i + 7) + 9 * j + 6];

        m[vars_per_node * j + 2] += d2C[csize * (9 * i + 5) + 9 * j + 1] -
                                    d2C[csize * (9 * i + 5) + 9 * j + 3] -
                                    d2C[csize * (9 * i + 7) + 9 * j + 1] +
                                    d2C[csize * (9 * i + 7) + 9 * j + 3];

        m[vars_per_node * j + size] += -d2C[csize * (9 * i + 2) + 9 * j + 5] +
                                       d2C[csize * (9 * i + 2) + 9 * j + 7] +
                                       d2C[csize * (9 * i + 6) + 9 * j + 5] -
                                       d2C[csize * (9 * i + 6) + 9 * j + 7];

        m[vars_per_node * j + 1 + size] +=
            d2C[csize * (9 * i + 2) + 9 * j + 2] -
            d2C[csize * (9 * i + 2) + 9 * j + 6] -
            d2C[csize * (9 * i + 6) + 9 * j + 2] +
            d2C[csize * (9 * i + 6) + 9 * j + 6];

        m[vars_per_node * j + 2 + size] +=
            -d2C[csize * (9 * i + 2) + 9 * j + 1] +
            d2C[csize * (9 * i + 2) + 9 * j + 3] +
            d2C[csize * (9 * i + 6) + 9 * j + 1] -
            d2C[csize * (9 * i + 6) + 9 * j + 3];

        m[vars_per_node * j + 2 * size] +=
            d2C[csize * (9 * i + 1) + 9 * j + 5] -
            d2C[csize * (9 * i + 1) + 9 * j + 7] -
            d2C[csize * (9 * i + 3) + 9 * j + 5] +
            d2C[csize * (9 * i + 3) + 9 * j + 7];

        m[vars_per_node * j + 1 + 2 * size] +=
            -d2C[csize * (9 * i + 1) + 9 * j + 2] +
            d2C[csize * (9 * i + 1) + 9 * j + 6] +
            d2C[csize * (9 * i + 3) + 9 * j + 2] -
            d2C[csize * (9 * i + 3) + 9 * j + 6];

        m[vars_per_node * j + 2 + 2 * size] +=
            d2C[csize * (9 * i + 1) + 9 * j + 1] -
            d2C[csize * (9 * i + 1) + 9 * j + 3] -
            d2C[csize * (9 * i + 3) + 9 * j + 1] +
            d2C[csize * (9 * i + 3) + 9 * j + 3];
      }

      m += vars_per_node * size;
      dC += 9;
    }
  }

  /**
    The linearized rotation class is unconstrained
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationConstraint(const TacsScalar vars[], TacsScalar res[]) {
  }

  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationConstrJacobian(const TacsScalar alpha,
                                        const TacsScalar vars[],
                                        TacsScalar res[], TacsScalar mat[]) {}

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
  static void computeDirectorRates(const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar t[], TacsScalar d[],
                                   TacsScalar ddot[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    for (int i = 0; i < num_nodes; i++) {
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
  static void computeDirectorRates(const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar ddvars[],
                                   const TacsScalar t[], TacsScalar d[],
                                   TacsScalar ddot[], TacsScalar dddot[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    const TacsScalar *qddot = &ddvars[offset];
    for (int i = 0; i < num_nodes; i++) {
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
    @param psi The full variable vector derivative
    @param t The reference directions
    @param C The rotation matrices at each point
    @param d The director values
    @param ddot The first time derivative of the director
    @param dddot The second time derivative of the director
    @param dd The derivative of the director values
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeDirectorRatesDeriv(const TacsScalar vars[],
                                        const TacsScalar dvars[],
                                        const TacsScalar ddvars[],
                                        const TacsScalar psi[],
                                        const TacsScalar t[], TacsScalar d[],
                                        TacsScalar ddot[], TacsScalar dddot[],
                                        TacsScalar dpsi[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    const TacsScalar *qddot = &ddvars[offset];
    const TacsScalar *qpsi = &psi[offset];
    for (int i = 0; i < num_nodes; i++) {
      crossProduct(q, t, d);
      crossProduct(qdot, t, ddot);
      crossProduct(qddot, t, dddot);
      crossProduct(qpsi, t, dpsi);

      t += 3;
      d += 3;
      ddot += 3;
      dddot += 3;
      dpsi += 3;

      q += vars_per_node;
      qdot += vars_per_node;
      qddot += vars_per_node;
      qpsi += vars_per_node;
    }
  }

  /**
    Given the derivatives of the kinetic energy expression with
    respect to time, add the contributions to the derivative of the

    Given the partial derivatives of the Lagrangian with respect to the
    director and the time derivative of the vector, compute

    dd = d/dt(dT/d(dot{d})) - dL/dd

    In general, the residual contribution is:

    res +=
    dTdot*d(dot{d})/d(dot{q}) +
    dT/d(dot{d})*d/dt(dot{d})/d(dot{q}) +
    dd*d(d)/d(q)

    For the linearized rotation director these expressions are:

    d = q^{x} t
    dot{d} = - t^{x} dot{q}
    d(dot{d})/d(dot{q}) = - t^{x}
    d/dt(d(dot{d})/d(dot{q})) = 0

    @param vars The full variable vector
    @param dvars The first time derivative of the variables
    @param ddvars The second derivatives of the variables
    @param t The normal direction
    @param dd The contribution from the derivative of the director
    @param res The output residual
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorResidual(const TacsScalar vars[],
                                  const TacsScalar dvars[],
                                  const TacsScalar ddvars[],
                                  const TacsScalar t[], const TacsScalar dd[],
                                  TacsScalar res[]) {
    TacsScalar *r = &res[offset];

    for (int i = 0; i < num_nodes; i++) {
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
  static void addDirectorJacobian(
      TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
      const TacsScalar vars[], const TacsScalar dvars[],
      const TacsScalar ddvars[], const TacsScalar t[], const TacsScalar dd[],
      const TacsScalar d2Tdotd[], const TacsScalar d2Tdotu[],
      const TacsScalar d2d[], const TacsScalar d2du[], TacsScalar res[],
      TacsScalar mat[]) {
    // Add the derivative due to and d2d
    const int dsize = 3 * num_nodes;
    const int nvars = vars_per_node * num_nodes;

    // d = crossProduct(q, t, d)
    const TacsScalar *ti = t;
    for (int i = 0; i < num_nodes; i++, ti += 3) {
      TacsScalar *jac1 = &mat[(offset + vars_per_node * i) * nvars + offset];
      TacsScalar *jac2 =
          &mat[(offset + vars_per_node * i + 1) * nvars + offset];
      TacsScalar *jac3 =
          &mat[(offset + vars_per_node * i + 2) * nvars + offset];

      const TacsScalar *tj = t;
      for (int j = 0; j < num_nodes; j++, tj += 3) {
        // Add the derivative
        TacsScalar d[9];
        d[0] = d2d[0] + gamma * d2Tdotd[0];
        d[1] = d2d[1] + gamma * d2Tdotd[1];
        d[2] = d2d[2] + gamma * d2Tdotd[2];

        d[3] = d2d[dsize] + gamma * d2Tdotd[dsize];
        d[4] = d2d[dsize + 1] + gamma * d2Tdotd[dsize + 1];
        d[5] = d2d[dsize + 2] + gamma * d2Tdotd[dsize + 2];

        d[6] = d2d[2 * dsize] + gamma * d2Tdotd[2 * dsize];
        d[7] = d2d[2 * dsize + 1] + gamma * d2Tdotd[2 * dsize + 1];
        d[8] = d2d[2 * dsize + 2] + gamma * d2Tdotd[2 * dsize + 2];

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
        d2Tdotd += 3;
      }

      d2d += 2 * dsize;
      d2Tdotd += 2 * dsize;
    }

    for (int i = 0; i < num_nodes; i++) {
      for (int j = 0; j < num_nodes; j++) {
        // Add the derivative
        TacsScalar d[9];
        d[0] = d2du[0] + gamma * d2Tdotu[0];
        d[1] = d2du[1] + gamma * d2Tdotu[1];
        d[2] = d2du[2] + gamma * d2Tdotu[2];

        d[3] = d2du[dsize] + gamma * d2Tdotu[dsize];
        d[4] = d2du[dsize + 1] + gamma * d2Tdotu[dsize + 1];
        d[5] = d2du[dsize + 2] + gamma * d2Tdotu[dsize + 2];

        d[6] = d2du[2 * dsize] + gamma * d2Tdotu[2 * dsize];
        d[7] = d2du[2 * dsize + 1] + gamma * d2Tdotu[2 * dsize + 1];
        d[8] = d2du[2 * dsize + 2] + gamma * d2Tdotu[2 * dsize + 2];

        TacsScalar tmp[9];
        mat3x3SkewMatTransform(&t[3 * i], d, tmp);

        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            int index = (vars_per_node * i + ii + offset) * nvars +
                        vars_per_node * j + jj;

            mat[index] += tmp[3 * ii + jj];
          }
        }

        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            int index = (vars_per_node * j + jj) * nvars + vars_per_node * i +
                        ii + offset;

            mat[index] += tmp[3 * ii + jj];
          }
        }

        d2du += 3;
        d2Tdotu += 3;
      }

      d2du += 2 * dsize;
      d2Tdotu += 2 * dsize;
    }

    // Update residual
    TacsScalar *r = &res[offset];

    for (int i = 0; i < num_nodes; i++) {
      crossProductAdd(1.0, t, dd, r);

      r += vars_per_node;
      dd += 3;
      t += 3;
    }
  }

  /*
    Add the director contributions to the derivative of the normal

    Add the adjoint sensitivity of the reference normal (dt) based on
    the adjoint sensitivity of the director (dd).

    Given that the parametrization is d = (C^{T}(q) - I) * t, compute

    dt += d(dd^{T}d)/dt = dd^{T}*C^{T}(q) = C(q) * dd

    @param vars The full variable vector
    @param t The reference directions
    @param dd The adjoint sensitivities w.r.t. the director
    @param dt The adjoint sensitivity w.r.t. the reference directions
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorRefNormalSens(const TacsScalar vars[],
                                       const TacsScalar t[],
                                       const TacsScalar dd[], TacsScalar dt[]) {
    const TacsScalar *q = &vars[offset];

    for (int i = 0; i < num_nodes; i++) {
      crossProductAdd(1.0, dd, q, dt);

      t += 3;
      dd += 3;
      dt += 3;
      q += vars_per_node;
    }
  }

  /*
    Add the director contributions to the derivative of the normal

    Add the adjoint sensitivity of the reference normal (dt) based on
    the adjoint sensitivity of the director (dd) and the sensitivity
    of the derivative field (ddadj).

    Given that the parametrization is d = (C^{T}(q) - I) * t and the field
    dpsi = d(d)/dq^{T} * psi, compute

    dt += d(dd^{T}d)/dt = dd^{T}*C^{T}(q) = C(q) * dd

    dt += d(ddpsi^{T}*dpsi)/dt = [ d(C(q))/dq * psi ] * ddpsi

    @param vars The full variable vector
    @param psi The full variable vector derivative
    @param t The reference directions
    @param dd The adjoint sensitivities w.r.t. the director
    @param ddpsi The adjoint sensitivities w.r.t. the director derivative
    @param dt The adjoint sensitivity w.r.t. the reference directions
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorRefNormalSens(
      const TacsScalar vars[], const TacsScalar psi[], const TacsScalar t[],
      const TacsScalar dd[], const TacsScalar ddpsi[], TacsScalar dt[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qpsi = &psi[offset];

    for (int i = 0; i < num_nodes; i++) {
      crossProductAdd(1.0, dd, q, dt);
      crossProductAdd(1.0, ddpsi, qpsi, dt);

      t += 3;
      dd += 3;
      ddpsi += 3;
      dt += 3;
      q += vars_per_node;
      qpsi += vars_per_node;
    }
  }

  static TacsScalar evalDrillStrain(const TacsScalar u0x[],
                                    const TacsScalar Ct[]) {
    // Compute the rotational penalty
    return 0.5 * (Ct[3] + u0x[3] - Ct[1] - u0x[1]);
  }

  static void evalDrillStrainSens(TacsScalar scale, const TacsScalar u0x[],
                                  const TacsScalar Ct[], TacsScalar du0x[],
                                  TacsScalar dCt[]) {
    dCt[0] = 0.0;
    dCt[1] = -0.5 * scale;
    dCt[2] = 0.0;
    dCt[3] = 0.5 * scale;
    dCt[4] = 0.0;
    dCt[5] = 0.0;
    dCt[6] = 0.0;
    dCt[7] = 0.0;
    dCt[8] = 0.0;

    du0x[0] = 0.0;
    du0x[1] = -0.5 * scale;
    du0x[2] = 0.0;
    du0x[3] = 0.5 * scale;
    du0x[4] = 0.0;
    du0x[5] = 0.0;
    du0x[6] = 0.0;
    du0x[7] = 0.0;
    du0x[8] = 0.0;
  }

  static TacsScalar evalDrillStrainDeriv(const TacsScalar u0x[],
                                         const TacsScalar Ct[],
                                         const TacsScalar u0xd[],
                                         const TacsScalar Ctd[],
                                         TacsScalar *ed) {
    *ed = 0.5 * (Ctd[3] + u0xd[3] - Ctd[1] - u0xd[1]);

    // Compute the rotational penalty
    return 0.5 * (Ct[3] + u0x[3] - Ct[1] - u0x[1]);
  }

  static void evalDrillStrainHessian(TacsScalar d2et, const TacsScalar u0x[],
                                     const TacsScalar Ct[], TacsScalar d2u0x[],
                                     TacsScalar d2Ct[], TacsScalar d2Ctu0x[]) {
    memset(d2u0x, 0, 81 * sizeof(TacsScalar));
    memset(d2Ct, 0, 81 * sizeof(TacsScalar));
    memset(d2Ctu0x, 0, 81 * sizeof(TacsScalar));
  }
};

/*
  A quadratic approximation of the director

  The director matrix is given as

  C(q) = 1 - q^{x} + 0.5*q^{x}q^{x}
  .    = 1 - q^{x} + 0.5*(qq^{T} - q^{T}q*1)
  .    = 1*(1 - 0.5*q^{T}q) - q^{x} + 0.5*q*q^{T}
*/
class TACSQuadraticRotation {
 public:
  static const int NUM_PARAMETERS = 3;

  /**
    Compute the rotation matrices at each node

    @param vars The full variable vector
    @param C The rotation matrices at each point
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeRotationMat(const TacsScalar vars[], TacsScalar C[]) {
    const TacsScalar *q = &vars[offset];
    for (int i = 0; i < num_nodes; i++) {
      TacsScalar qTq = vec3Dot(q, q);
      setMatSkew(-1.0, q, C);
      C[0] = C[4] = C[8] = 1.0 - 0.5 * qTq;
      vec3x3OuterAdd(0.5, q, q, C);

      C += 9;
      q += vars_per_node;
    }
  }

  /**
    Compute the derivative of the rotation matrices at each node

    @param vars The full variable vector
    @param varsd The full variable vector
    @param C The rotation matrices at each point
    @param Cd The rotation matrices at each point
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeRotationMatDeriv(const TacsScalar vars[],
                                      const TacsScalar varsd[], TacsScalar C[],
                                      TacsScalar Cd[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qd = &varsd[offset];
    for (int i = 0; i < num_nodes; i++) {
      // Compute C
      TacsScalar qTq = vec3Dot(q, q);
      setMatSkew(-1.0, q, C);
      C[0] = C[4] = C[8] = 1.0 - 0.5 * qTq;
      vec3x3OuterAdd(0.5, q, q, C);

      // Compute Cd
      TacsScalar qTqd = vec3Dot(q, qd);
      setMatSkew(-1.0, qd, Cd);
      Cd[0] = Cd[4] = Cd[8] = -qTqd;
      vec3x3OuterAdd(0.5, qd, q, Cd);
      vec3x3OuterAdd(0.5, q, qd, Cd);

      C += 9;
      Cd += 9;

      q += vars_per_node;
      qd += vars_per_node;
    }
  }

  /**
    Add the residual rotation matrix to the output

    This code adds the contribution to the residual via the derivative

    d(tr(dC^{T}C(q)))/dq_{i}
    = d(tr(dC^{T}*(I - q^{x} + 0.5*(qq^{T} - q^{T}q*1))))/dq_{i}

    @param vars The full variable vector
    @param dC The derivative w.r.t. the rotation matrix
    @param res The residual array
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationMatResidual(const TacsScalar vars[],
                                     const TacsScalar dC[], TacsScalar res[]) {
    const TacsScalar *q = &vars[offset];
    TacsScalar *r = &res[offset];

    for (int i = 0; i < num_nodes; i++) {
      TacsScalar dCtr = (dC[0] + dC[4] + dC[8]);
      r[0] -= dC[7] - dC[5] + dCtr * q[0];
      r[1] -= dC[2] - dC[6] + dCtr * q[1];
      r[2] -= dC[3] - dC[1] + dCtr * q[2];

      TacsScalar e1[3], e2[3];
      mat3x3Mult(dC, q, e1);
      mat3x3MultTrans(dC, q, e2);

      r[0] += 0.5 * (e1[0] + e2[0]);
      r[1] += 0.5 * (e1[1] + e2[1]);
      r[2] += 0.5 * (e1[2] + e2[2]);

      r += vars_per_node;
      q += vars_per_node;
      dC += 9;
    }
  }

  /*
    Add the Jacobian of the rotation matrix to the output

    @param alpha Scalar coefficient for the Jacobian matrix
    @param vars The variable values
    @param dC The derivative of the functional w.r.t. C
    @param d2C The second derivatives of the functional w.r.t. C
    @param res The residual
    @param mat The Jacobian matrix
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationMatJacobian(const TacsScalar alpha,
                                     const TacsScalar vars[],
                                     const TacsScalar dC[],
                                     const TacsScalar d2C[], TacsScalar res[],
                                     TacsScalar mat[]) {
    const int size = vars_per_node * num_nodes;
    const int csize = 9 * num_nodes;

    const TacsScalar *q = &vars[offset];
    TacsScalar *r = NULL;
    if (res) {
      r = &res[offset];
    }
    TacsScalar *m = &mat[offset * size + offset];

    for (int i = 0; i < num_nodes; i++) {
      // Add the contribution to the residual
      if (res) {
        TacsScalar dCtr = (dC[0] + dC[4] + dC[8]);
        r[0] -= dC[7] - dC[5] + dCtr * q[0];
        r[1] -= dC[2] - dC[6] + dCtr * q[1];
        r[2] -= dC[3] - dC[1] + dCtr * q[2];

        TacsScalar e1[3], e2[3];
        mat3x3Mult(dC, q, e1);
        mat3x3MultTrans(dC, q, e2);

        r[0] += 0.5 * (e1[0] + e2[0]);
        r[1] += 0.5 * (e1[1] + e2[1]);
        r[2] += 0.5 * (e1[2] + e2[2]);
        r += vars_per_node;
      }

      const TacsScalar *qi = &vars[offset + i * vars_per_node];
      const TacsScalar *qj = &vars[offset];

      for (int j = 0; j < num_nodes; j++, qj += vars_per_node) {
        TacsScalar dfdC[27];
        for (int k = 0; k < 9; k++) {
          TacsScalar d2Ct[9];
          for (int jj = 0; jj < 9; jj++) {
            d2Ct[jj] = d2C[csize * (9 * i + k) + 9 * j + jj];
          }

          TacsScalar d2Ctr = (d2Ct[0] + d2Ct[4] + d2Ct[8]);
          dfdC[k] = -(d2Ct[7] - d2Ct[5] + d2Ctr * qj[0]);
          dfdC[9 + k] = -(d2Ct[2] - d2Ct[6] + d2Ctr * qj[1]);
          dfdC[18 + k] = -(d2Ct[3] - d2Ct[1] + d2Ctr * qj[2]);

          TacsScalar e1[3], e2[3];
          mat3x3Mult(d2Ct, qj, e1);
          mat3x3MultTrans(d2Ct, qj, e2);

          dfdC[k] += 0.5 * (e1[0] + e2[0]);
          dfdC[9 + k] += 0.5 * (e1[1] + e2[1]);
          dfdC[18 + k] += 0.5 * (e1[2] + e2[2]);
        }

        TacsScalar jac[9];
        for (int k = 0; k < 3; k++) {
          TacsScalar *d2Ct = &dfdC[9 * k];

          TacsScalar d2Ctr = (d2Ct[0] + d2Ct[4] + d2Ct[8]);
          jac[k] = -(d2Ct[7] - d2Ct[5] + d2Ctr * qi[0]);
          jac[3 + k] = -(d2Ct[2] - d2Ct[6] + d2Ctr * qi[1]);
          jac[6 + k] = -(d2Ct[3] - d2Ct[1] + d2Ctr * qi[2]);

          TacsScalar e1[3], e2[3];
          mat3x3Mult(d2Ct, qi, e1);
          mat3x3MultTrans(d2Ct, qi, e2);

          jac[k] += 0.5 * (e1[0] + e2[0]);
          jac[3 + k] += 0.5 * (e1[1] + e2[1]);
          jac[6 + k] += 0.5 * (e1[2] + e2[2]);
        }

        if (i == j) {
          jac[0] -= dC[4] + dC[8];
          jac[1] += 0.5 * (dC[1] + dC[3]);
          jac[2] += 0.5 * (dC[2] + dC[6]);

          jac[3] += 0.5 * (dC[1] + dC[3]);
          jac[4] -= dC[0] + dC[8];
          jac[5] += 0.5 * (dC[5] + dC[7]);

          jac[6] += 0.5 * (dC[2] + dC[6]);
          jac[7] += 0.5 * (dC[5] + dC[7]);
          jac[8] -= dC[0] + dC[4];
        }

        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            m[vars_per_node * j + ii * size + jj] += jac[3 * ii + jj];
          }
        }
      }

      m += vars_per_node * size;
      q += vars_per_node;
      dC += 9;
    }
  }

  /**
    The quadratic rotation matrix is unconstrained
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationConstraint(const TacsScalar vars[], TacsScalar res[]) {
  }

  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationConstrJacobian(const TacsScalar alpha,
                                        const TacsScalar vars[],
                                        TacsScalar res[], TacsScalar mat[]) {}

  /**
    Compute the director and rates at all nodes.

    d = Q(q)*t = (1 + 0.5*q^{x})*q^{x}t
    ddot = qdot^{x}*t + 0.5*q^{x}*qdot^{x}*t + 0.5*qdot^{x}*q^{x}*t

    @param vars The full variable vector
    @param dvars The first time derivative of the variables
    @param ddvars The second derivatives of the variables
    @param t The reference directions
    @param d The director values
    @param ddot The first time derivative of the director
    @param dddot The second time derivative of the director
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeDirectorRates(const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar t[], TacsScalar d[],
                                   TacsScalar ddot[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    for (int i = 0; i < num_nodes; i++) {
      TacsScalar qxt[3], qxtdot[3];

      // Compute d = q^{x}*t + 0.5*q^{x}*q^{x}*t
      crossProduct(q, t, qxt);
      d[0] = qxt[0];
      d[1] = qxt[1];
      d[2] = qxt[2];
      crossProductAdd(0.5, q, qxt, d);

      // Compute ddot =
      // qdot^{x}*t + 0.5*qdot^{x}*q^{x}*t + 0.5*q^{x}*qdot^{x}*t
      crossProduct(qdot, t, qxtdot);
      ddot[0] = qxtdot[0];
      ddot[1] = qxtdot[1];
      ddot[2] = qxtdot[2];
      crossProductAdd(0.5, qdot, qxt, ddot);
      crossProductAdd(0.5, q, qxtdot, ddot);

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
  static void computeDirectorRates(const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar ddvars[],
                                   const TacsScalar t[], TacsScalar d[],
                                   TacsScalar ddot[], TacsScalar dddot[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    const TacsScalar *qddot = &ddvars[offset];
    for (int i = 0; i < num_nodes; i++) {
      TacsScalar qxt[3], qxtdot[3], qxtddot[3];

      // Compute d = q^{x}*t + 0.5*q^{x}*q^{x}*t
      crossProduct(q, t, qxt);
      d[0] = qxt[0];
      d[1] = qxt[1];
      d[2] = qxt[2];
      crossProductAdd(0.5, q, qxt, d);

      // Compute ddot =
      // qdot^{x}*t + 0.5*qdot^{x}*q^{x}*t + 0.5*q^{x}*qdot^{x}*t
      crossProduct(qdot, t, qxtdot);
      ddot[0] = qxtdot[0];
      ddot[1] = qxtdot[1];
      ddot[2] = qxtdot[2];
      crossProductAdd(0.5, qdot, qxt, ddot);
      crossProductAdd(0.5, q, qxtdot, ddot);

      // Compute ddot =
      // qddot^{x}*t + 0.5*qddot^{x}*q^{x}*t +
      // qdot^{x}*qdot^{x}*t + 0.5*q^{x}*qddot^{x}*t
      crossProduct(qddot, t, qxtddot);
      dddot[0] = qxtddot[0];
      dddot[1] = qxtddot[1];
      dddot[2] = qxtddot[2];
      crossProductAdd(0.5, qddot, qxt, dddot);
      crossProductAdd(1.0, qdot, qxtdot, dddot);
      crossProductAdd(0.5, q, qxtddot, dddot);

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
    @param psi The full variable vector derivative
    @param t The reference directions
    @param C The rotation matrices at each point
    @param d The director values
    @param ddot The first time derivative of the director
    @param dddot The second time derivative of the director
    @param dd The derivative of the director values
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeDirectorRatesDeriv(const TacsScalar vars[],
                                        const TacsScalar dvars[],
                                        const TacsScalar ddvars[],
                                        const TacsScalar psi[],
                                        const TacsScalar t[], TacsScalar d[],
                                        TacsScalar ddot[], TacsScalar dddot[],
                                        TacsScalar dpsi[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    const TacsScalar *qddot = &ddvars[offset];
    const TacsScalar *qpsi = &psi[offset];
    for (int i = 0; i < num_nodes; i++) {
      TacsScalar qxt[3], qxtdot[3], qxtddot[3], qdxt[3];

      // Compute d = q^{x}*t + 0.5*q^{x}*q^{x}*t
      crossProduct(q, t, qxt);
      d[0] = qxt[0];
      d[1] = qxt[1];
      d[2] = qxt[2];
      crossProductAdd(0.5, q, qxt, d);

      // Compute ddot =
      // qdot^{x}*t + 0.5*qdot^{x}*q^{x}*t + 0.5*q^{x}*qdot^{x}*t
      crossProduct(qdot, t, qxtdot);
      ddot[0] = qxtdot[0];
      ddot[1] = qxtdot[1];
      ddot[2] = qxtdot[2];
      crossProductAdd(0.5, qdot, qxt, ddot);
      crossProductAdd(0.5, q, qxtdot, ddot);

      // Compute ddot =
      // qddot^{x}*t + 0.5*qddot^{x}*q^{x}*t +
      // qdot^{x}*qdot^{x}*t + 0.5*q^{x}*qddot^{x}*t
      crossProduct(qddot, t, qxtddot);
      dddot[0] = qxtddot[0];
      dddot[1] = qxtddot[1];
      dddot[2] = qxtddot[2];
      crossProductAdd(0.5, qddot, qxt, dddot);
      crossProductAdd(1.0, qdot, qxtdot, dddot);
      crossProductAdd(0.5, q, qxtddot, dddot);

      // Compute dpsi = (qpsi^{x} + 0.5*qpsi^{x}*q^{x} + 0.5*q^{x}*qpsi^{x})*t
      crossProduct(qpsi, t, qdxt);
      dpsi[0] = qdxt[0];
      dpsi[1] = qdxt[1];
      dpsi[2] = qdxt[2];
      crossProductAdd(0.5, qpsi, qxt, dpsi);
      crossProductAdd(0.5, q, qdxt, dpsi);

      t += 3;
      d += 3;
      ddot += 3;
      dddot += 3;
      dpsi += 3;

      q += vars_per_node;
      qdot += vars_per_node;
      qddot += vars_per_node;
      qpsi += vars_per_node;
    }
  }

  /**
    Given the derivatives of the kinetic energy expression with
    respect to time, add the contributions to the derivative of the

    Given the partial derivatives of the Lagrangian with respect to the
    director and the time derivative of the vector, compute

    dTdot = d/dt(dT/d(dot{d}))
    dd = -dL/dd

    @param vars The full variable vector
    @param dvars The first time derivative of the variables
    @param ddvars The second derivatives of the variables
    @param t The normal direction
    @param dd The contribution from the derivative of the director
    @param res The output residual
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorResidual(const TacsScalar vars[],
                                  const TacsScalar dvars[],
                                  const TacsScalar ddvars[],
                                  const TacsScalar t[], const TacsScalar dd[],
                                  TacsScalar res[]) {
    TacsScalar *r = &res[offset];
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];

    for (int i = 0; i < num_nodes; i++) {
      // D = (t^{x} + 0.5*(q^{x}t)^{x} - 0.5*t^{x}*q^{x})
      //   = (t^{x} - t^x*q^{x} + 0.5*q^{x}*t^{x})
      TacsScalar D[9];
      D[0] = 0.5 * (q[1] * t[1] + q[2] * t[2]);
      D[1] = -q[0] * t[1] + 0.5 * q[1] * t[0] - t[2];
      D[2] = -q[0] * t[2] + 0.5 * q[2] * t[0] + t[1];

      D[3] = 0.5 * q[0] * t[1] - q[1] * t[0] + t[2];
      D[4] = 0.5 * (q[0] * t[0] + q[2] * t[2]);
      D[5] = -q[1] * t[2] + 0.5 * q[2] * t[1] - t[0];

      D[6] = 0.5 * q[0] * t[2] - q[2] * t[0] - t[1];
      D[7] = 0.5 * q[1] * t[2] - q[2] * t[1] + t[0];
      D[8] = 0.5 * (q[0] * t[0] + q[1] * t[1]);

      // Add the contribution to the residual
      r[0] += D[0] * dd[0] + D[1] * dd[1] + D[2] * dd[2];
      r[1] += D[3] * dd[0] + D[4] * dd[1] + D[5] * dd[2];
      r[2] += D[6] * dd[0] + D[7] * dd[1] + D[8] * dd[2];

      r += vars_per_node;
      q += vars_per_node;
      qdot += vars_per_node;

      dd += 3;
      t += 3;
    }
  }

  /*
    Add terms from the Jacobian
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorJacobian(
      TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
      const TacsScalar vars[], const TacsScalar dvars[],
      const TacsScalar ddvars[], const TacsScalar t[], const TacsScalar dd[],
      const TacsScalar d2Tdotd[], const TacsScalar d2Tdotu[],
      const TacsScalar d2d[], const TacsScalar d2du[], TacsScalar res[],
      TacsScalar mat[]) {
    const int size = vars_per_node * num_nodes;
    const int dsize = 3 * num_nodes;

    // Pre-compute the Jacobian matrices
    const TacsScalar *ti = t;
    TacsScalar D[9 * num_nodes], Ddot[9 * num_nodes], Dddot[9 * num_nodes];
    for (int i = 0; i < num_nodes; i++, ti += 3) {
      TacsScalar *Di = &D[9 * i];
      TacsScalar *Didot = &Ddot[9 * i];
      TacsScalar *Diddot = &Dddot[9 * i];
      const TacsScalar *qi = &vars[offset + vars_per_node * i];
      const TacsScalar *qidot = &dvars[offset + vars_per_node * i];
      const TacsScalar *qiddot = &ddvars[offset + vars_per_node * i];

      Di[0] = 0.5 * (qi[1] * ti[1] + qi[2] * ti[2]);
      Di[1] = -qi[0] * ti[1] + 0.5 * qi[1] * ti[0] - ti[2];
      Di[2] = -qi[0] * ti[2] + 0.5 * qi[2] * ti[0] + ti[1];

      Di[3] = 0.5 * qi[0] * ti[1] - qi[1] * ti[0] + ti[2];
      Di[4] = 0.5 * (qi[0] * ti[0] + qi[2] * ti[2]);
      Di[5] = -qi[1] * ti[2] + 0.5 * qi[2] * ti[1] - ti[0];

      Di[6] = 0.5 * qi[0] * ti[2] - qi[2] * ti[0] - ti[1];
      Di[7] = 0.5 * qi[1] * ti[2] - qi[2] * ti[1] + ti[0];
      Di[8] = 0.5 * (qi[0] * ti[0] + qi[1] * ti[1]);

      Didot[0] = 0.5 * (qidot[1] * ti[1] + qidot[2] * ti[2]);
      Didot[1] = -qidot[0] * ti[1] + 0.5 * qidot[1] * ti[0];
      Didot[2] = -qidot[0] * ti[2] + 0.5 * qidot[2] * ti[0];

      Didot[3] = 0.5 * qidot[0] * ti[1] - qidot[1] * ti[0];
      Didot[4] = 0.5 * (qidot[0] * ti[0] + qidot[2] * ti[2]);
      Didot[5] = -qidot[1] * ti[2] + 0.5 * qidot[2] * ti[1];

      Didot[6] = 0.5 * qidot[0] * ti[2] - qidot[2] * ti[0];
      Didot[7] = 0.5 * qidot[1] * ti[2] - qidot[2] * ti[1];
      Didot[8] = 0.5 * (qidot[0] * ti[0] + qidot[1] * ti[1]);

      Diddot[0] = 0.5 * (qiddot[1] * ti[1] + qiddot[2] * ti[2]);
      Diddot[1] = -qiddot[0] * ti[1] + 0.5 * qiddot[1] * ti[0];
      Diddot[2] = -qiddot[0] * ti[2] + 0.5 * qiddot[2] * ti[0];

      Diddot[3] = 0.5 * qiddot[0] * ti[1] - qiddot[1] * ti[0];
      Diddot[4] = 0.5 * (qiddot[0] * ti[0] + qiddot[2] * ti[2]);
      Diddot[5] = -qiddot[1] * ti[2] + 0.5 * qiddot[2] * ti[1];

      Diddot[6] = 0.5 * qiddot[0] * ti[2] - qiddot[2] * ti[0];
      Diddot[7] = 0.5 * qiddot[1] * ti[2] - qiddot[2] * ti[1];
      Diddot[8] = 0.5 * (qiddot[0] * ti[0] + qiddot[1] * ti[1]);
    }

    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    TacsScalar *r = NULL;
    if (res) {
      r = &res[offset];
    }
    TacsScalar *m = &mat[offset * size + offset];

    // Re-set the pointer
    const TacsScalar *Di = D;
    for (int i = 0; i < num_nodes; i++, Di += 9) {
      if (res) {
        r[0] += D[0] * dd[0] + D[1] * dd[1] + D[2] * dd[2];
        r[1] += D[3] * dd[0] + D[4] * dd[1] + D[5] * dd[2];
        r[2] += D[6] * dd[0] + D[7] * dd[1] + D[8] * dd[2];

        r += vars_per_node;
      }

      const TacsScalar *Dj = D;
      const TacsScalar *Djdot = Ddot;
      const TacsScalar *Djddot = Dddot;
      for (int j = 0; j < num_nodes; j++, Dj += 9, Djdot += 9, Djddot += 9) {
        TacsScalar dfdq[9];
        for (int k = 0; k < 3; k++) {
          TacsScalar tmp[3];
          tmp[0] = d2d[dsize * (3 * i + k) + 3 * j] +
                   gamma * d2Tdotd[dsize * (3 * i + k) + 3 * j];
          tmp[1] = d2d[dsize * (3 * i + k) + 3 * j + 1] +
                   gamma * d2Tdotd[dsize * (3 * i + k) + 3 * j + 1];
          tmp[2] = d2d[dsize * (3 * i + k) + 3 * j + 2] +
                   gamma * d2Tdotd[dsize * (3 * i + k) + 3 * j + 2];

          dfdq[k] = Dj[0] * tmp[0] + Dj[1] * tmp[1] + Dj[2] * tmp[2];
          dfdq[3 + k] = Dj[3] * tmp[0] + Dj[4] * tmp[1] + Dj[5] * tmp[2];
          dfdq[6 + k] = Dj[6] * tmp[0] + Dj[7] * tmp[1] + Dj[8] * tmp[2];

          tmp[0] = beta * d2Tdotd[dsize * (3 * i + k) + 3 * j];
          tmp[1] = beta * d2Tdotd[dsize * (3 * i + k) + 3 * j + 1];
          tmp[2] = beta * d2Tdotd[dsize * (3 * i + k) + 3 * j + 2];

          dfdq[k] +=
              2.0 * (Djdot[0] * tmp[0] + Djdot[1] * tmp[1] + Djdot[2] * tmp[2]);
          dfdq[3 + k] +=
              2.0 * (Djdot[3] * tmp[0] + Djdot[4] * tmp[1] + Djdot[5] * tmp[2]);
          dfdq[6 + k] +=
              2.0 * (Djdot[6] * tmp[0] + Djdot[7] * tmp[1] + Djdot[8] * tmp[2]);

          tmp[0] = alpha * d2Tdotd[dsize * (3 * i + k) + 3 * j];
          tmp[1] = alpha * d2Tdotd[dsize * (3 * i + k) + 3 * j + 1];
          tmp[2] = alpha * d2Tdotd[dsize * (3 * i + k) + 3 * j + 2];

          dfdq[k] +=
              Djddot[0] * tmp[0] + Djddot[1] * tmp[1] + Djddot[2] * tmp[2];
          dfdq[3 + k] +=
              Djddot[3] * tmp[0] + Djddot[4] * tmp[1] + Djddot[5] * tmp[2];
          dfdq[6 + k] +=
              Djddot[6] * tmp[0] + Djddot[7] * tmp[1] + Djddot[8] * tmp[2];
        }

        TacsScalar jac[9];
        for (int k = 0; k < 3; k++) {
          jac[k] = dfdq[3 * k] * Di[0] + dfdq[3 * k + 1] * Di[1] +
                   dfdq[3 * k + 2] * Di[2];
          jac[3 + k] = dfdq[3 * k] * Di[3] + dfdq[3 * k + 1] * Di[4] +
                       dfdq[3 * k + 2] * Di[5];
          jac[6 + k] = dfdq[3 * k] * Di[6] + dfdq[3 * k + 1] * Di[7] +
                       dfdq[3 * k + 2] * Di[8];
        }

        if (i == j) {
          TacsScalar tmp[3];
          tmp[0] = alpha * dd[0];
          tmp[1] = alpha * dd[1];
          tmp[2] = alpha * dd[2];

          jac[0] -= tmp[1] * t[1] + tmp[2] * t[2];
          jac[1] += 0.5 * (tmp[0] * t[1] + tmp[1] * t[0]);
          jac[2] += 0.5 * (tmp[0] * t[2] + tmp[2] * t[0]);

          jac[3] += 0.5 * (tmp[0] * t[1] + tmp[1] * t[0]);
          jac[4] -= tmp[0] * t[0] + tmp[2] * t[2];
          jac[5] += 0.5 * (tmp[1] * t[2] + tmp[2] * t[1]);

          jac[6] += 0.5 * (tmp[0] * t[2] + tmp[2] * t[0]);
          jac[7] += 0.5 * (tmp[1] * t[2] + tmp[2] * t[1]);
          jac[8] -= tmp[0] * t[0] + tmp[1] * t[1];
        }

        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            m[vars_per_node * j + ii * size + jj] += jac[3 * ii + jj];
          }
        }
      }

      q += vars_per_node;
      qdot += vars_per_node;

      dd += 3;
      t += 3;
      m += vars_per_node * size;
    }

    Di = D;
    const TacsScalar *Didot = Ddot;
    const TacsScalar *Diddot = Dddot;
    for (int i = 0; i < num_nodes; i++, Di += 9, Didot += 9, Diddot += 9) {
      for (int j = 0; j < num_nodes; j++) {
        TacsScalar dfdq[9], dfdqT[9];
        for (int k = 0; k < 3; k++) {
          TacsScalar tmp[3];
          tmp[0] = d2du[dsize * (3 * i) + 3 * j + k] +
                   gamma * d2Tdotu[dsize * (3 * i) + 3 * j + k];
          tmp[1] = d2du[dsize * (3 * i + 1) + 3 * j + k] +
                   gamma * d2Tdotu[dsize * (3 * i + 1) + 3 * j + k];
          tmp[2] = d2du[dsize * (3 * i + 2) + 3 * j + k] +
                   gamma * d2Tdotu[dsize * (3 * i + 2) + 3 * j + k];

          dfdq[k] = Di[0] * tmp[0] + Di[1] * tmp[1] + Di[2] * tmp[2];
          dfdq[3 + k] = Di[3] * tmp[0] + Di[4] * tmp[1] + Di[5] * tmp[2];
          dfdq[6 + k] = Di[6] * tmp[0] + Di[7] * tmp[1] + Di[8] * tmp[2];

          tmp[0] = d2du[dsize * (3 * i) + 3 * j + k];
          tmp[1] = d2du[dsize * (3 * i + 1) + 3 * j + k];
          tmp[2] = d2du[dsize * (3 * i + 2) + 3 * j + k];

          dfdqT[k] = dfdq[k];
          dfdqT[3 + k] = dfdq[3 + k];
          dfdqT[6 + k] = dfdq[6 + k];

          tmp[0] = beta * d2Tdotu[dsize * (3 * i) + 3 * j + k];
          tmp[1] = beta * d2Tdotu[dsize * (3 * i + 1) + 3 * j + k];
          tmp[2] = beta * d2Tdotu[dsize * (3 * i + 2) + 3 * j + k];

          dfdqT[k] +=
              2.0 * (Didot[0] * tmp[0] + Didot[1] * tmp[1] + Didot[2] * tmp[2]);
          dfdqT[3 + k] +=
              2.0 * (Didot[3] * tmp[0] + Didot[4] * tmp[1] + Didot[5] * tmp[2]);
          dfdqT[6 + k] +=
              2.0 * (Didot[6] * tmp[0] + Didot[7] * tmp[1] + Didot[8] * tmp[2]);

          tmp[0] = alpha * d2Tdotu[dsize * (3 * i) + 3 * j + k];
          tmp[1] = alpha * d2Tdotu[dsize * (3 * i + 1) + 3 * j + k];
          tmp[2] = alpha * d2Tdotu[dsize * (3 * i + 2) + 3 * j + k];

          dfdqT[k] +=
              Diddot[0] * tmp[0] + Diddot[1] * tmp[1] + Diddot[2] * tmp[2];
          dfdqT[3 + k] +=
              Diddot[3] * tmp[0] + Diddot[4] * tmp[1] + Diddot[5] * tmp[2];
          dfdqT[6 + k] +=
              Diddot[6] * tmp[0] + Diddot[7] * tmp[1] + Diddot[8] * tmp[2];
        }

        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            int index = (vars_per_node * i + ii + offset) * size +
                        vars_per_node * j + jj;

            mat[index] += dfdq[3 * ii + jj];
          }
        }

        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            int index = (vars_per_node * j + jj) * size + vars_per_node * i +
                        ii + offset;

            mat[index] += dfdqT[3 * ii + jj];
          }
        }
      }
    }
  }

  /*
    Add the director contributions to the derivative of the normal

    Add the adjoint sensitivity of the reference normal (dt) based on
    the adjoint sensitivity of the director (dd).

    Given that the parametrization is d = (C^{T}(q) - I) * t, compute

    dt += d(dd^{T} * d)/dt = dd^{T} * C^{T}(q) = C(q) * dd
    .   = (0.5*q^{x} - 1)*q^{x} * dpsi

    @param vars The full variable vector
    @param t The reference directions
    @param dd The adjoint sensitivities w.r.t. the director
    @param dt The adjoint sensitivity w.r.t. the reference directions
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorRefNormalSens(const TacsScalar vars[],
                                       const TacsScalar t[],
                                       const TacsScalar dd[], TacsScalar dt[]) {
    const TacsScalar *q = &vars[offset];

    for (int i = 0; i < num_nodes; i++) {
      TacsScalar tmp[3];
      crossProduct(q, dd, tmp);
      dt[0] -= tmp[0];
      dt[1] -= tmp[1];
      dt[2] -= tmp[2];
      crossProductAdd(0.5, q, tmp, dt);

      t += 3;
      dd += 3;
      dt += 3;
      q += vars_per_node;
    }
  }

  /*
    Add the director contributions to the derivative of the normal

    Add the adjoint sensitivity of the reference normal (dt) based on
    the adjoint sensitivity of the director (dd) and the sensitivity
    of the derivative field (ddadj).

    Given that the parametrization is d = (C^{T}(q) - I) * t and the field
    dpsi = d(d)/dq^{T} * psi, compute

    dt += d(dd^{T}d)/dn = dd^{T} * C^{T}(q) = C(q) * dd
    .   = (0.5*q^{x} - 1)*q^{x} * dpsi

    dt += d(ddpsi^{T} * dpsi)/dt = [ d(C(q))/dq * psi ] * ddpsi
    .   = (qpsi^{x} + 0.5*qpsi^{x}*q^{x} + 0.5*q^{x}*qpsi^{x})^{T} * ddpsi
    .   = (0.5*qpsi^{x}*q^{x} + 0.5*q^{x}*qpsi^{x}) - qpsi^{x}) * ddpsi

    @param vars The full variable vector
    @param psi The full variable vector derivative
    @param t The reference directions
    @param dd The adjoint sensitivities w.r.t. the director
    @param ddpsi The adjoint sensitivities w.r.t. the director derivative
    @param dt The adjoint sensitivity w.r.t. the reference directions
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorRefNormalSens(
      const TacsScalar vars[], const TacsScalar psi[], const TacsScalar t[],
      const TacsScalar dd[], const TacsScalar ddpsi[], TacsScalar dt[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qpsi = &psi[offset];

    for (int i = 0; i < num_nodes; i++) {
      TacsScalar tmp[3];
      crossProduct(q, dd, tmp);
      dt[0] -= tmp[0];
      dt[1] -= tmp[1];
      dt[2] -= tmp[2];
      crossProductAdd(0.5, q, tmp, dt);

      // tmp = qpsi^{x} * ddpsi
      crossProduct(qpsi, ddpsi, tmp);
      dt[0] -= tmp[0];
      dt[1] -= tmp[1];
      dt[2] -= tmp[2];
      // Add dt += 0.5 * q^{x} * tmp = 0.5 * q^{x} * qpsi^{x} * ddpsi
      crossProductAdd(0.5, q, tmp, dt);

      // tmp = q^{x} * ddpsi
      crossProduct(q, ddpsi, tmp);
      crossProductAdd(0.5, qpsi, tmp, dt);

      t += 3;
      dd += 3;
      ddpsi += 3;
      dt += 3;
      q += vars_per_node;
      qpsi += vars_per_node;
    }
  }

  static TacsScalar evalDrillStrain(const TacsScalar u0x[],
                                    const TacsScalar Ct[]) {
    // Compute the rotational penalty
    return 0.5 * (Ct[3] + u0x[3] - Ct[1] - u0x[1]);
  }

  static void evalDrillStrainSens(TacsScalar scale, const TacsScalar u0x[],
                                  const TacsScalar Ct[], TacsScalar du0x[],
                                  TacsScalar dCt[]) {
    dCt[0] = 0.0;
    dCt[1] = -0.5 * scale;
    dCt[2] = 0.0;
    dCt[3] = 0.5 * scale;
    dCt[4] = 0.0;
    dCt[5] = 0.0;
    dCt[6] = 0.0;
    dCt[7] = 0.0;
    dCt[8] = 0.0;

    du0x[0] = 0.0;
    du0x[1] = -0.5 * scale;
    du0x[2] = 0.0;
    du0x[3] = 0.5 * scale;
    du0x[4] = 0.0;
    du0x[5] = 0.0;
    du0x[6] = 0.0;
    du0x[7] = 0.0;
    du0x[8] = 0.0;
  }

  static TacsScalar evalDrillStrainDeriv(const TacsScalar u0x[],
                                         const TacsScalar Ct[],
                                         const TacsScalar u0xd[],
                                         const TacsScalar Ctd[],
                                         TacsScalar *ed) {
    *ed = 0.5 * (Ctd[3] + u0xd[3] - Ctd[1] - u0xd[1]);

    // Compute the rotational penalty
    return 0.5 * (Ct[3] + u0x[3] - Ct[1] - u0x[1]);
  }

  static void evalDrillStrainHessian(TacsScalar d2et, const TacsScalar u0x[],
                                     const TacsScalar Ct[], TacsScalar d2u0x[],
                                     TacsScalar d2Ct[], TacsScalar d2Ctu0x[]) {
    memset(d2u0x, 0, 81 * sizeof(TacsScalar));
    memset(d2Ct, 0, 81 * sizeof(TacsScalar));
    memset(d2Ctu0x, 0, 81 * sizeof(TacsScalar));
  }
};

class TACSQuaternionRotation {
 public:
  static const int NUM_PARAMETERS = 5;

  /**
    Compute the rotation matrices at each node

    @param vars The full variable vector
    @param C The rotation matrices at each point
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeRotationMat(const TacsScalar vars[], TacsScalar C[]) {
    const TacsScalar *q = &vars[offset];
    for (int i = 0; i < num_nodes; i++) {
      C[0] = 1.0 - 2.0 * (q[2] * q[2] + q[3] * q[3]);
      C[1] = 2.0 * (q[1] * q[2] + q[3] * q[0]);
      C[2] = 2.0 * (q[1] * q[3] - q[2] * q[0]);

      C[3] = 2.0 * (q[2] * q[1] - q[3] * q[0]);
      C[4] = 1.0 - 2.0 * (q[1] * q[1] + q[3] * q[3]);
      C[5] = 2.0 * (q[2] * q[3] + q[1] * q[0]);

      C[6] = 2.0 * (q[3] * q[1] + q[2] * q[0]);
      C[7] = 2.0 * (q[3] * q[2] - q[1] * q[0]);
      C[8] = 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);

      C += 9;
      q += vars_per_node;
    }
  }

  /**
    Compute the derivative of the rotation matrices at each node

    @param vars The full variable vector
    @param varsd The full variable vector
    @param C The rotation matrices at each point
    @param Cd The rotation matrices at each point
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeRotationMatDeriv(const TacsScalar vars[],
                                      const TacsScalar varsd[], TacsScalar C[],
                                      TacsScalar Cd[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qd = &varsd[offset];

    for (int i = 0; i < num_nodes; i++) {
      C[0] = 1.0 - 2.0 * (q[2] * q[2] + q[3] * q[3]);
      C[1] = 2.0 * (q[1] * q[2] + q[3] * q[0]);
      C[2] = 2.0 * (q[1] * q[3] - q[2] * q[0]);

      C[3] = 2.0 * (q[2] * q[1] - q[3] * q[0]);
      C[4] = 1.0 - 2.0 * (q[1] * q[1] + q[3] * q[3]);
      C[5] = 2.0 * (q[2] * q[3] + q[1] * q[0]);

      C[6] = 2.0 * (q[3] * q[1] + q[2] * q[0]);
      C[7] = 2.0 * (q[3] * q[2] - q[1] * q[0]);
      C[8] = 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);

      Cd[0] = -4.0 * (q[2] * qd[2] + q[3] * qd[3]);
      Cd[1] = 2.0 * (q[1] * qd[2] + q[3] * qd[0] + qd[1] * q[2] + qd[3] * q[0]);
      Cd[2] = 2.0 * (q[1] * qd[3] - q[2] * qd[0] + qd[1] * q[3] - qd[2] * q[0]);

      Cd[3] = 2.0 * (q[2] * qd[1] - q[3] * qd[0] + qd[2] * q[1] - qd[3] * q[0]);
      Cd[4] = -4.0 * (q[1] * qd[1] + q[3] * qd[3]);
      Cd[5] = 2.0 * (q[2] * qd[3] + q[1] * qd[0] + qd[2] * q[3] + qd[1] * q[0]);

      Cd[6] = 2.0 * (q[3] * qd[1] + q[2] * qd[0] + qd[3] * q[1] + qd[2] * q[0]);
      Cd[7] = 2.0 * (q[3] * qd[2] - q[1] * qd[0] + qd[3] * q[2] - qd[1] * q[0]);
      Cd[8] = -4.0 * (q[1] * qd[1] + q[2] * qd[2]);

      C += 9;
      Cd += 9;

      q += vars_per_node;
      qd += vars_per_node;
    }
  }

  /**
    Add the residual rotation matrix to the output

    This code adds the contribution to the residual via the derivative

    d(tr(dC^{T}C(q)))/dq_{i}

    @param vars The full variable vector
    @param dC The derivative w.r.t. the rotation matrix
    @param res The residual array
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationMatResidual(const TacsScalar vars[],
                                     const TacsScalar dC[], TacsScalar res[]) {
    const TacsScalar *q = &vars[offset];
    TacsScalar *r = &res[offset];

    for (int i = 0; i < num_nodes; i++) {
      r[0] += 2.0 * (q[3] * (dC[1] - dC[3]) + q[2] * (dC[6] - dC[2]) +
                     q[1] * (dC[5] - dC[7]));
      r[1] += 2.0 * (q[0] * (dC[5] - dC[7]) - 2.0 * q[1] * (dC[4] + dC[8]) +
                     q[2] * (dC[1] + dC[3]) + q[3] * (dC[2] + dC[6]));
      r[2] += 2.0 * (q[0] * (dC[6] - dC[2]) + q[1] * (dC[1] + dC[3]) -
                     2.0 * q[2] * (dC[0] + dC[8]) + q[3] * (dC[7] + dC[5]));
      r[3] += 2.0 * (q[0] * (dC[1] - dC[3]) + q[1] * (dC[2] + dC[6]) +
                     q[2] * (dC[5] + dC[7]) - 2.0 * q[3] * (dC[0] + dC[4]));

      r += vars_per_node;
      q += vars_per_node;
      dC += 9;
    }
  }

  /*
    Add the Jacobian of the rotation matrix to the output

    @param alpha Scalar coefficient for the Jacobian matrix
    @param vars The variable values
    @param dC The derivative of the functional w.r.t. C
    @param d2C The second derivatives of the functional w.r.t. C
    @param res The residual
    @param mat The Jacobian matrix
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationMatJacobian(const TacsScalar alpha,
                                     const TacsScalar vars[],
                                     const TacsScalar dC[],
                                     const TacsScalar d2C[], TacsScalar res[],
                                     TacsScalar mat[]) {
    const int size = vars_per_node * num_nodes;
    const int csize = 9 * num_nodes;

    const TacsScalar *q = &vars[offset];
    TacsScalar *r = NULL;
    if (res) {
      r = &res[offset];
    }
    TacsScalar *m = &mat[offset * size + offset];

    for (int i = 0; i < num_nodes; i++) {
      if (res) {
        r[0] += 2.0 * (q[3] * (dC[1] - dC[3]) + q[2] * (dC[6] - dC[2]) +
                       q[1] * (dC[5] - dC[7]));
        r[1] += 2.0 * (q[0] * (dC[5] - dC[7]) - 2.0 * q[1] * (dC[4] + dC[8]) +
                       q[2] * (dC[1] + dC[3]) + q[3] * (dC[2] + dC[6]));
        r[2] += 2.0 * (q[0] * (dC[6] - dC[2]) + q[1] * (dC[1] + dC[3]) -
                       2.0 * q[2] * (dC[0] + dC[8]) + q[3] * (dC[7] + dC[5]));
        r[3] += 2.0 * (q[0] * (dC[1] - dC[3]) + q[1] * (dC[2] + dC[6]) +
                       q[2] * (dC[5] + dC[7]) - 2.0 * q[3] * (dC[0] + dC[4]));
        r += vars_per_node;
      }

      const TacsScalar *qi = &vars[offset + i * vars_per_node];
      const TacsScalar *qj = &vars[offset];

      for (int j = 0; j < num_nodes; j++, qj += vars_per_node) {
        TacsScalar dfdC[36];
        for (int k = 0; k < 9; k++) {
          TacsScalar d2Ct[9];
          for (int jj = 0; jj < 9; jj++) {
            d2Ct[jj] = d2C[csize * (9 * i + k) + 9 * j + jj];
          }

          dfdC[k] =
              2.0 * (qj[3] * (d2Ct[1] - d2Ct[3]) + qj[2] * (d2Ct[6] - d2Ct[2]) +
                     qj[1] * (d2Ct[5] - d2Ct[7]));
          dfdC[9 + k] =
              2.0 *
              (qj[0] * (d2Ct[5] - d2Ct[7]) - 2.0 * qj[1] * (d2Ct[4] + d2Ct[8]) +
               qj[2] * (d2Ct[1] + d2Ct[3]) + qj[3] * (d2Ct[2] + d2Ct[6]));
          dfdC[18 + k] =
              2.0 *
              (qj[0] * (d2Ct[6] - d2Ct[2]) + qj[1] * (d2Ct[1] + d2Ct[3]) -
               2.0 * qj[2] * (d2Ct[0] + d2Ct[8]) + qj[3] * (d2Ct[7] + d2Ct[5]));
          dfdC[27 + k] =
              2.0 *
              (qj[0] * (d2Ct[1] - d2Ct[3]) + qj[1] * (d2Ct[2] + d2Ct[6]) +
               qj[2] * (d2Ct[5] + d2Ct[7]) - 2.0 * qj[3] * (d2Ct[0] + d2Ct[4]));
        }

        TacsScalar jac[16];
        for (int k = 0; k < 4; k++) {
          TacsScalar *d2Ct = &dfdC[9 * k];

          jac[k] =
              2.0 * (qi[3] * (d2Ct[1] - d2Ct[3]) + qi[2] * (d2Ct[6] - d2Ct[2]) +
                     qi[1] * (d2Ct[5] - d2Ct[7]));
          jac[4 + k] =
              2.0 *
              (qi[0] * (d2Ct[5] - d2Ct[7]) - 2.0 * qi[1] * (d2Ct[4] + d2Ct[8]) +
               qi[2] * (d2Ct[1] + d2Ct[3]) + qi[3] * (d2Ct[2] + d2Ct[6]));
          jac[8 + k] =
              2.0 *
              (qi[0] * (d2Ct[6] - d2Ct[2]) + qi[1] * (d2Ct[1] + d2Ct[3]) -
               2.0 * qi[2] * (d2Ct[0] + d2Ct[8]) + qi[3] * (d2Ct[7] + d2Ct[5]));
          jac[12 + k] =
              2.0 *
              (qi[0] * (d2Ct[1] - d2Ct[3]) + qi[1] * (d2Ct[2] + d2Ct[6]) +
               qi[2] * (d2Ct[5] + d2Ct[7]) - 2.0 * qi[3] * (d2Ct[0] + d2Ct[4]));
        }

        if (i == j) {
          jac[1] += 2.0 * (dC[5] - dC[7]);
          jac[2] += 2.0 * (dC[6] - dC[2]);
          jac[3] += 2.0 * (dC[1] - dC[3]);

          jac[4] += 2.0 * (dC[5] - dC[7]);
          jac[5] -= 4.0 * (dC[4] + dC[8]);
          jac[6] += 2.0 * (dC[1] + dC[3]);
          jac[7] += 2.0 * (dC[2] + dC[6]);

          jac[8] += 2.0 * (dC[6] - dC[2]);
          jac[9] += 2.0 * (dC[1] + dC[3]);
          jac[10] -= 4.0 * (dC[0] + dC[8]);
          jac[11] += 2.0 * (dC[5] + dC[7]);

          jac[12] += 2.0 * (dC[1] - dC[3]);
          jac[13] += 2.0 * (dC[2] + dC[6]);
          jac[14] += 2.0 * (dC[5] + dC[7]);
          jac[15] -= 4.0 * (dC[0] + dC[4]);
        }

        for (int ii = 0; ii < 4; ii++) {
          for (int jj = 0; jj < 4; jj++) {
            m[vars_per_node * j + ii * size + jj] += jac[4 * ii + jj];
          }
        }
      }

      m += vars_per_node * size;
      q += vars_per_node;
      dC += 9;
    }
  }

  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationConstraint(const TacsScalar vars[], TacsScalar res[]) {
    const TacsScalar *q = &vars[offset];
    TacsScalar *r = &res[offset];
    for (int i = 0; i < num_nodes; i++) {
      TacsScalar lamb = q[4];

      // Add the result to the governing equations
      r[0] += 2.0 * q[0] * lamb;
      r[1] += 2.0 * q[1] * lamb;
      r[2] += 2.0 * q[2] * lamb;
      r[3] += 2.0 * q[3] * lamb;

      // Enforce the quaternion constraint
      r[4] += q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3] - 1.0;

      r += vars_per_node;
      q += vars_per_node;
    }
  }

  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationConstrJacobian(const TacsScalar alpha,
                                        const TacsScalar vars[],
                                        TacsScalar res[], TacsScalar mat[]) {
    const int size = vars_per_node * num_nodes;
    const TacsScalar *q = &vars[offset];
    TacsScalar *r = NULL;
    if (res) {
      r = &res[offset];
    }

    // Add the constribution to the constraints from the quaternions
    TacsScalar *m = &mat[offset * size + offset];
    for (int i = 0; i < num_nodes; i++) {
      TacsScalar lamb = q[4];

      if (res) {
        // Add the result to the governing equations
        r[0] += 2.0 * q[0] * lamb;
        r[1] += 2.0 * q[1] * lamb;
        r[2] += 2.0 * q[2] * lamb;
        r[3] += 2.0 * q[3] * lamb;

        // Enforce the quaternion constraint
        r[4] += q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3] - 1.0;

        r += vars_per_node;
      }

      // Add the constraint terms
      m[4] += 2.0 * alpha * q[0];
      m[4 + size] += 2.0 * alpha * q[1];
      m[4 + 2 * size] += 2.0 * alpha * q[2];
      m[4 + 3 * size] += 2.0 * alpha * q[3];

      // Enforce the quaternion constraint
      m[4 * size] += 2.0 * alpha * q[0];
      m[4 * size + 1] += 2.0 * alpha * q[1];
      m[4 * size + 2] += 2.0 * alpha * q[2];
      m[4 * size + 3] += 2.0 * alpha * q[3];

      // Add the terms to the diagonal
      m[0] += 2.0 * alpha * lamb;
      m[size + 1] += 2.0 * alpha * lamb;
      m[2 * (size + 1)] += 2.0 * alpha * lamb;
      m[3 * (size + 1)] += 2.0 * alpha * lamb;

      q += vars_per_node;

      // Increment to the next block diagonal entry
      m += vars_per_node * (size + 1);
    }
  }

  /**
    Compute the director and rates at all nodes.

    d = Q(q)*t = (1 + 0.5*q^{x})*q^{x}t
    ddot = qdot^{x}*t + 0.5*q^{x}*qdot^{x}*t + 0.5*qdot^{x}*q^{x}*t

    @param vars The full variable vector
    @param dvars The first time derivative of the variables
    @param ddvars The second derivatives of the variables
    @param t The reference directions
    @param d The director values
    @param ddot The first time derivative of the director
    @param dddot The second time derivative of the director
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeDirectorRates(const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar t[], TacsScalar d[],
                                   TacsScalar ddot[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];

    for (int i = 0; i < num_nodes; i++) {
      TacsScalar Q[9];
      Q[0] = -2.0 * (q[2] * q[2] + q[3] * q[3]);
      Q[1] = 2.0 * (q[2] * q[1] - q[3] * q[0]);
      Q[2] = 2.0 * (q[3] * q[1] + q[2] * q[0]);

      Q[3] = 2.0 * (q[1] * q[2] + q[3] * q[0]);
      Q[4] = -2.0 * (q[1] * q[1] + q[3] * q[3]);
      Q[5] = 2.0 * (q[3] * q[2] - q[1] * q[0]);

      Q[6] = 2.0 * (q[1] * q[3] - q[2] * q[0]);
      Q[7] = 2.0 * (q[2] * q[3] + q[1] * q[0]);
      Q[8] = -2.0 * (q[1] * q[1] + q[2] * q[2]);

      // Compute d = Q*t
      d[0] = Q[0] * t[0] + Q[1] * t[1] + Q[2] * t[2];
      d[1] = Q[3] * t[0] + Q[4] * t[1] + Q[5] * t[2];
      d[2] = Q[6] * t[0] + Q[7] * t[1] + Q[8] * t[2];

      TacsScalar Qdot[9];
      Qdot[0] = -4.0 * (q[2] * qdot[2] + q[3] * qdot[3]);
      Qdot[1] = 2.0 * (q[2] * qdot[1] - q[3] * qdot[0] + qdot[2] * q[1] -
                       qdot[3] * q[0]);
      Qdot[2] = 2.0 * (q[3] * qdot[1] + q[2] * qdot[0] + qdot[3] * q[1] +
                       qdot[2] * q[0]);

      Qdot[3] = 2.0 * (q[1] * qdot[2] + q[3] * qdot[0] + qdot[1] * q[2] +
                       qdot[3] * q[0]);
      Qdot[4] = -4.0 * (q[1] * qdot[1] + q[3] * qdot[3]);
      Qdot[5] = 2.0 * (q[3] * qdot[2] - q[1] * qdot[0] + qdot[3] * q[2] -
                       qdot[1] * q[0]);

      Qdot[6] = 2.0 * (q[1] * qdot[3] - q[2] * qdot[0] + qdot[1] * q[3] -
                       qdot[2] * q[0]);
      Qdot[7] = 2.0 * (q[2] * qdot[3] + q[1] * qdot[0] + qdot[2] * q[3] +
                       qdot[1] * q[0]);
      Qdot[8] = -4.0 * (q[1] * qdot[1] + q[2] * qdot[2]);

      // Compute d = Q*t
      ddot[0] = Qdot[0] * t[0] + Qdot[1] * t[1] + Qdot[2] * t[2];
      ddot[1] = Qdot[3] * t[0] + Qdot[4] * t[1] + Qdot[5] * t[2];
      ddot[2] = Qdot[6] * t[0] + Qdot[7] * t[1] + Qdot[8] * t[2];

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
  static void computeDirectorRates(const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar ddvars[],
                                   const TacsScalar t[], TacsScalar d[],
                                   TacsScalar ddot[], TacsScalar dddot[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    const TacsScalar *qddot = &ddvars[offset];

    for (int i = 0; i < num_nodes; i++) {
      TacsScalar Q[9];
      Q[0] = -2.0 * (q[2] * q[2] + q[3] * q[3]);
      Q[1] = 2.0 * (q[2] * q[1] - q[3] * q[0]);
      Q[2] = 2.0 * (q[3] * q[1] + q[2] * q[0]);

      Q[3] = 2.0 * (q[1] * q[2] + q[3] * q[0]);
      Q[4] = -2.0 * (q[1] * q[1] + q[3] * q[3]);
      Q[5] = 2.0 * (q[3] * q[2] - q[1] * q[0]);

      Q[6] = 2.0 * (q[1] * q[3] - q[2] * q[0]);
      Q[7] = 2.0 * (q[2] * q[3] + q[1] * q[0]);
      Q[8] = -2.0 * (q[1] * q[1] + q[2] * q[2]);

      // Compute d = Q*t
      d[0] = Q[0] * t[0] + Q[1] * t[1] + Q[2] * t[2];
      d[1] = Q[3] * t[0] + Q[4] * t[1] + Q[5] * t[2];
      d[2] = Q[6] * t[0] + Q[7] * t[1] + Q[8] * t[2];

      TacsScalar Qdot[9];
      Qdot[0] = -4.0 * (q[2] * qdot[2] + q[3] * qdot[3]);
      Qdot[1] = 2.0 * (q[2] * qdot[1] - q[3] * qdot[0] + qdot[2] * q[1] -
                       qdot[3] * q[0]);
      Qdot[2] = 2.0 * (q[3] * qdot[1] + q[2] * qdot[0] + qdot[3] * q[1] +
                       qdot[2] * q[0]);

      Qdot[3] = 2.0 * (q[1] * qdot[2] + q[3] * qdot[0] + qdot[1] * q[2] +
                       qdot[3] * q[0]);
      Qdot[4] = -4.0 * (q[1] * qdot[1] + q[3] * qdot[3]);
      Qdot[5] = 2.0 * (q[3] * qdot[2] - q[1] * qdot[0] + qdot[3] * q[2] -
                       qdot[1] * q[0]);

      Qdot[6] = 2.0 * (q[1] * qdot[3] - q[2] * qdot[0] + qdot[1] * q[3] -
                       qdot[2] * q[0]);
      Qdot[7] = 2.0 * (q[2] * qdot[3] + q[1] * qdot[0] + qdot[2] * q[3] +
                       qdot[1] * q[0]);
      Qdot[8] = -4.0 * (q[1] * qdot[1] + q[2] * qdot[2]);

      // Compute ddot = Qdot*t
      ddot[0] = Qdot[0] * t[0] + Qdot[1] * t[1] + Qdot[2] * t[2];
      ddot[1] = Qdot[3] * t[0] + Qdot[4] * t[1] + Qdot[5] * t[2];
      ddot[2] = Qdot[6] * t[0] + Qdot[7] * t[1] + Qdot[8] * t[2];

      TacsScalar Qddot[9];
      Qddot[0] = -4.0 * (q[2] * qddot[2] + q[3] * qddot[3] + qdot[2] * qdot[2] +
                         qdot[3] * qdot[3]);
      Qddot[1] =
          2.0 * (q[2] * qddot[1] - q[3] * qddot[0] + qddot[2] * q[1] -
                 qddot[3] * q[0] + qdot[2] * qdot[1] - qdot[3] * qdot[0] +
                 qdot[2] * qdot[1] - qdot[3] * qdot[0]);
      Qddot[2] =
          2.0 * (q[3] * qddot[1] + q[2] * qddot[0] + qddot[3] * q[1] +
                 qddot[2] * q[0] + qdot[3] * qdot[1] + qdot[2] * qdot[0] +
                 qdot[3] * qdot[1] + qdot[2] * qdot[0]);

      Qddot[3] =
          2.0 * (q[1] * qddot[2] + q[3] * qddot[0] + qddot[1] * q[2] +
                 qddot[3] * q[0] + qdot[1] * qdot[2] + qdot[3] * qdot[0] +
                 qdot[1] * qdot[2] + qdot[3] * qdot[0]);
      Qddot[4] = -4.0 * (q[1] * qddot[1] + q[3] * qddot[3] + qdot[1] * qdot[1] +
                         qdot[3] * qdot[3]);
      Qddot[5] =
          2.0 * (q[3] * qddot[2] - q[1] * qddot[0] + qddot[3] * q[2] -
                 qddot[1] * q[0] + qdot[3] * qdot[2] - qdot[1] * qdot[0] +
                 qdot[3] * qdot[2] - qdot[1] * qdot[0]);

      Qddot[6] =
          2.0 * (q[1] * qddot[3] - q[2] * qddot[0] + qddot[1] * q[3] -
                 qddot[2] * q[0] + qdot[1] * qdot[3] - qdot[2] * qdot[0] +
                 qdot[1] * qdot[3] - qdot[2] * qdot[0]);
      Qddot[7] =
          2.0 * (q[2] * qddot[3] + q[1] * qddot[0] + qddot[2] * q[3] +
                 qddot[1] * q[0] + qdot[2] * qdot[3] + qdot[1] * qdot[0] +
                 qdot[2] * qdot[3] + qdot[1] * qdot[0]);
      Qddot[8] = -4.0 * (q[1] * qddot[1] + q[2] * qddot[2] + qdot[1] * qdot[1] +
                         qdot[2] * qdot[2]);

      // Compute dddot = Qddot*t
      dddot[0] = Qddot[0] * t[0] + Qddot[1] * t[1] + Qddot[2] * t[2];
      dddot[1] = Qddot[3] * t[0] + Qddot[4] * t[1] + Qddot[5] * t[2];
      dddot[2] = Qddot[6] * t[0] + Qddot[7] * t[1] + Qddot[8] * t[2];

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
    @param dd The derivative of the director values
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeDirectorRatesDeriv(
      const TacsScalar vars[], const TacsScalar dvars[],
      const TacsScalar ddvars[], const TacsScalar varsd[], const TacsScalar t[],
      TacsScalar d[], TacsScalar ddot[], TacsScalar dddot[], TacsScalar dd[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    const TacsScalar *qddot = &ddvars[offset];
    const TacsScalar *qd = &varsd[offset];

    for (int i = 0; i < num_nodes; i++) {
      TacsScalar Q[9];
      Q[0] = -2.0 * (q[2] * q[2] + q[3] * q[3]);
      Q[1] = 2.0 * (q[2] * q[1] - q[3] * q[0]);
      Q[2] = 2.0 * (q[3] * q[1] + q[2] * q[0]);

      Q[3] = 2.0 * (q[1] * q[2] + q[3] * q[0]);
      Q[4] = -2.0 * (q[1] * q[1] + q[3] * q[3]);
      Q[5] = 2.0 * (q[3] * q[2] - q[1] * q[0]);

      Q[6] = 2.0 * (q[1] * q[3] - q[2] * q[0]);
      Q[7] = 2.0 * (q[2] * q[3] + q[1] * q[0]);
      Q[8] = -2.0 * (q[1] * q[1] + q[2] * q[2]);

      // Compute d = Q*t
      d[0] = Q[0] * t[0] + Q[1] * t[1] + Q[2] * t[2];
      d[1] = Q[3] * t[0] + Q[4] * t[1] + Q[5] * t[2];
      d[2] = Q[6] * t[0] + Q[7] * t[1] + Q[8] * t[2];

      TacsScalar Qdot[9];
      Qdot[0] = -4.0 * (q[2] * qdot[2] + q[3] * qdot[3]);
      Qdot[1] = 2.0 * (q[2] * qdot[1] - q[3] * qdot[0] + qdot[2] * q[1] -
                       qdot[3] * q[0]);
      Qdot[2] = 2.0 * (q[3] * qdot[1] + q[2] * qdot[0] + qdot[3] * q[1] +
                       qdot[2] * q[0]);

      Qdot[3] = 2.0 * (q[1] * qdot[2] + q[3] * qdot[0] + qdot[1] * q[2] +
                       qdot[3] * q[0]);
      Qdot[4] = -4.0 * (q[1] * qdot[1] + q[3] * qdot[3]);
      Qdot[5] = 2.0 * (q[3] * qdot[2] - q[1] * qdot[0] + qdot[3] * q[2] -
                       qdot[1] * q[0]);

      Qdot[6] = 2.0 * (q[1] * qdot[3] - q[2] * qdot[0] + qdot[1] * q[3] -
                       qdot[2] * q[0]);
      Qdot[7] = 2.0 * (q[2] * qdot[3] + q[1] * qdot[0] + qdot[2] * q[3] +
                       qdot[1] * q[0]);
      Qdot[8] = -4.0 * (q[1] * qdot[1] + q[2] * qdot[2]);

      // Compute ddot = Qdot*t
      ddot[0] = Qdot[0] * t[0] + Qdot[1] * t[1] + Qdot[2] * t[2];
      ddot[1] = Qdot[3] * t[0] + Qdot[4] * t[1] + Qdot[5] * t[2];
      ddot[2] = Qdot[6] * t[0] + Qdot[7] * t[1] + Qdot[8] * t[2];

      TacsScalar Qddot[9];
      Qddot[0] = -4.0 * (q[2] * qddot[2] + q[3] * qddot[3] + qdot[2] * qdot[2] +
                         qdot[3] * qdot[3]);
      Qddot[1] =
          2.0 * (q[2] * qddot[1] - q[3] * qddot[0] + qddot[2] * q[1] -
                 qddot[3] * q[0] + qdot[2] * qdot[1] - qdot[3] * qdot[0] +
                 qdot[2] * qdot[1] - qdot[3] * qdot[0]);
      Qddot[2] =
          2.0 * (q[3] * qddot[1] + q[2] * qddot[0] + qddot[3] * q[1] +
                 qddot[2] * q[0] + qdot[3] * qdot[1] + qdot[2] * qdot[0] +
                 qdot[3] * qdot[1] + qdot[2] * qdot[0]);

      Qddot[3] =
          2.0 * (q[1] * qddot[2] + q[3] * qddot[0] + qddot[1] * q[2] +
                 qddot[3] * q[0] + qdot[1] * qdot[2] + qdot[3] * qdot[0] +
                 qdot[1] * qdot[2] + qdot[3] * qdot[0]);
      Qddot[4] = -4.0 * (q[1] * qddot[1] + q[3] * qddot[3] + qdot[1] * qdot[1] +
                         qdot[3] * qdot[3]);
      Qddot[5] =
          2.0 * (q[3] * qddot[2] - q[1] * qddot[0] + qddot[3] * q[2] -
                 qddot[1] * q[0] + qdot[3] * qdot[2] - qdot[1] * qdot[0] +
                 qdot[3] * qdot[2] - qdot[1] * qdot[0]);

      Qddot[6] =
          2.0 * (q[1] * qddot[3] - q[2] * qddot[0] + qddot[1] * q[3] -
                 qddot[2] * q[0] + qdot[1] * qdot[3] - qdot[2] * qdot[0] +
                 qdot[1] * qdot[3] - qdot[2] * qdot[0]);
      Qddot[7] =
          2.0 * (q[2] * qddot[3] + q[1] * qddot[0] + qddot[2] * q[3] +
                 qddot[1] * q[0] + qdot[2] * qdot[3] + qdot[1] * qdot[0] +
                 qdot[2] * qdot[3] + qdot[1] * qdot[0]);
      Qddot[8] = -4.0 * (q[1] * qddot[1] + q[2] * qddot[2] + qdot[1] * qdot[1] +
                         qdot[2] * qdot[2]);

      // Compute dddot = Qddot*t
      dddot[0] = Qddot[0] * t[0] + Qddot[1] * t[1] + Qddot[2] * t[2];
      dddot[1] = Qddot[3] * t[0] + Qddot[4] * t[1] + Qddot[5] * t[2];
      dddot[2] = Qddot[6] * t[0] + Qddot[7] * t[1] + Qddot[8] * t[2];

      Q[0] = -4.0 * (q[2] * qd[2] + q[3] * qd[3]);
      Q[1] = 2.0 * (q[2] * qd[1] - q[3] * qd[0] + qd[2] * q[1] - qd[3] * q[0]);
      Q[2] = 2.0 * (q[3] * qd[1] + q[2] * qd[0] + qd[3] * q[1] + qd[2] * q[0]);

      Q[3] = 2.0 * (q[1] * qd[2] + q[3] * qd[0] + qd[1] * q[2] + qd[3] * q[0]);
      Q[4] = -4.0 * (q[1] * qd[1] + q[3] * qd[3]);
      Q[5] = 2.0 * (q[3] * qd[2] - q[1] * qd[0] + qd[3] * q[2] - qd[1] * q[0]);

      Q[6] = 2.0 * (q[1] * qd[3] - q[2] * qd[0] + qd[1] * q[3] - qd[2] * q[0]);
      Q[7] = 2.0 * (q[2] * qd[3] + q[1] * qd[0] + qd[2] * q[3] + qd[1] * q[0]);
      Q[8] = -4.0 * (q[1] * qd[1] + q[2] * qd[2]);

      // Compute d = Q*t
      dd[0] = Q[0] * t[0] + Q[1] * t[1] + Q[2] * t[2];
      dd[1] = Q[3] * t[0] + Q[4] * t[1] + Q[5] * t[2];
      dd[2] = Q[6] * t[0] + Q[7] * t[1] + Q[8] * t[2];

      t += 3;
      d += 3;
      ddot += 3;
      dddot += 3;
      dd += 3;

      q += vars_per_node;
      qdot += vars_per_node;
      qddot += vars_per_node;
      qd += vars_per_node;
    }
  }

  /**
    Given the derivatives of the kinetic energy expression with
    respect to time, add the contributions to the derivative of the

    Given the partial derivatives of the Lagrangian with respect to the
    director and the time derivative of the vector, compute

    dd =  d/dt(dT/d(dot{d})) - dL/dd

    In general, the residual contribution is:

    res +=
    dTdot*d(dot{d})/d(dot{q}) + dd*d(d)/d(q)

    @param vars The full variable vector
    @param dvars The first time derivative of the variables
    @param ddvars The second derivatives of the variables
    @param t The normal direction
    @param dd The contribution from the derivative of the director
    @param res The output residual
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorResidual(const TacsScalar vars[],
                                  const TacsScalar dvars[],
                                  const TacsScalar ddvars[],
                                  const TacsScalar t[], const TacsScalar dd[],
                                  TacsScalar res[]) {
    TacsScalar *r = &res[offset];
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];

    for (int i = 0; i < num_nodes; i++) {
      // D = d(Qdot*t)/d(qdot)
      TacsScalar D[12];
      D[0] = 2.0 * (q[2] * t[2] - q[3] * t[1]);
      D[1] = 2.0 * (q[2] * t[1] + q[3] * t[2]);
      D[2] = 2.0 * (-2.0 * q[2] * t[0] + q[1] * t[1] + q[0] * t[2]);
      D[3] = 2.0 * (-2.0 * q[3] * t[0] - q[0] * t[1] + q[1] * t[2]);

      D[4] = 2.0 * (q[3] * t[0] - q[1] * t[2]);
      D[5] = 2.0 * (q[2] * t[0] - 2.0 * q[1] * t[1] - q[0] * t[2]);
      D[6] = 2.0 * (q[1] * t[0] + q[3] * t[2]);
      D[7] = 2.0 * (q[0] * t[0] - 2.0 * q[3] * t[1] + q[2] * t[2]);

      D[8] = 2.0 * (q[1] * t[1] - q[2] * t[0]);
      D[9] = 2.0 * (q[3] * t[0] + q[0] * t[1] - 2.0 * q[1] * t[2]);
      D[10] = 2.0 * (q[3] * t[1] - q[0] * t[0] - 2.0 * q[2] * t[2]);
      D[11] = 2.0 * (q[1] * t[0] + q[2] * t[1]);

      r[0] += dd[0] * D[0] + dd[1] * D[4] + dd[2] * D[8];
      r[1] += dd[0] * D[1] + dd[1] * D[5] + dd[2] * D[9];
      r[2] += dd[0] * D[2] + dd[1] * D[6] + dd[2] * D[10];
      r[3] += dd[0] * D[3] + dd[1] * D[7] + dd[2] * D[11];

      r += vars_per_node;
      q += vars_per_node;
      qdot += vars_per_node;

      dd += 3;
      t += 3;
    }
  }

  /*
    Add the contributions to the Jacobian matrix from the director
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorJacobian(
      TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
      const TacsScalar vars[], const TacsScalar dvars[],
      const TacsScalar ddvars[], const TacsScalar t[], const TacsScalar dd[],
      const TacsScalar d2Tdotd[], const TacsScalar d2Tdotu[],
      const TacsScalar d2d[], const TacsScalar d2du[], TacsScalar res[],
      TacsScalar mat[]) {
    const int size = vars_per_node * num_nodes;
    const int dsize = 3 * num_nodes;

    // Pre-compute the Jacobian matrices
    const TacsScalar *ti = t;
    TacsScalar D[12 * num_nodes], Ddot[12 * num_nodes], Dddot[12 * num_nodes];
    for (int i = 0; i < num_nodes; i++, ti += 3) {
      // D = d(Qdot*t)/d(qdot)
      TacsScalar *Di = &D[12 * i];
      TacsScalar *Didot = &Ddot[12 * i];
      TacsScalar *Diddot = &Dddot[12 * i];
      const TacsScalar *qi = &vars[offset + vars_per_node * i];
      const TacsScalar *qidot = &dvars[offset + vars_per_node * i];
      const TacsScalar *qiddot = &ddvars[offset + vars_per_node * i];

      Di[0] = 2.0 * (qi[2] * ti[2] - qi[3] * ti[1]);
      Di[1] = 2.0 * (qi[2] * ti[1] + qi[3] * ti[2]);
      Di[2] = 2.0 * (-2.0 * qi[2] * ti[0] + qi[1] * ti[1] + qi[0] * ti[2]);
      Di[3] = 2.0 * (-2.0 * qi[3] * ti[0] - qi[0] * ti[1] + qi[1] * ti[2]);

      Di[4] = 2.0 * (qi[3] * ti[0] - qi[1] * ti[2]);
      Di[5] = 2.0 * (qi[2] * ti[0] - 2.0 * qi[1] * ti[1] - qi[0] * ti[2]);
      Di[6] = 2.0 * (qi[1] * ti[0] + qi[3] * ti[2]);
      Di[7] = 2.0 * (qi[0] * ti[0] - 2.0 * qi[3] * ti[1] + qi[2] * ti[2]);

      Di[8] = 2.0 * (qi[1] * ti[1] - qi[2] * ti[0]);
      Di[9] = 2.0 * (qi[3] * ti[0] + qi[0] * ti[1] - 2.0 * qi[1] * ti[2]);
      Di[10] = 2.0 * (qi[3] * ti[1] - qi[0] * ti[0] - 2.0 * qi[2] * ti[2]);
      Di[11] = 2.0 * (qi[1] * ti[0] + qi[2] * ti[1]);

      Didot[0] = 2.0 * (qidot[2] * ti[2] - qidot[3] * ti[1]);
      Didot[1] = 2.0 * (qidot[2] * ti[1] + qidot[3] * ti[2]);
      Didot[2] =
          2.0 * (-2.0 * qidot[2] * ti[0] + qidot[1] * ti[1] + qidot[0] * ti[2]);
      Didot[3] =
          2.0 * (-2.0 * qidot[3] * ti[0] - qidot[0] * ti[1] + qidot[1] * ti[2]);

      Didot[4] = 2.0 * (qidot[3] * ti[0] - qidot[1] * ti[2]);
      Didot[5] =
          2.0 * (qidot[2] * ti[0] - 2.0 * qidot[1] * ti[1] - qidot[0] * ti[2]);
      Didot[6] = 2.0 * (qidot[1] * ti[0] + qidot[3] * ti[2]);
      Didot[7] =
          2.0 * (qidot[0] * ti[0] - 2.0 * qidot[3] * ti[1] + qidot[2] * ti[2]);

      Didot[8] = 2.0 * (qidot[1] * ti[1] - qidot[2] * ti[0]);
      Didot[9] =
          2.0 * (qidot[3] * ti[0] + qidot[0] * ti[1] - 2.0 * qidot[1] * ti[2]);
      Didot[10] =
          2.0 * (qidot[3] * ti[1] - qidot[0] * ti[0] - 2.0 * qidot[2] * ti[2]);
      Didot[11] = 2.0 * (qidot[1] * ti[0] + qidot[2] * ti[1]);

      Diddot[0] = 2.0 * (qiddot[2] * ti[2] - qiddot[3] * ti[1]);
      Diddot[1] = 2.0 * (qiddot[2] * ti[1] + qiddot[3] * ti[2]);
      Diddot[2] = 2.0 * (-2.0 * qiddot[2] * ti[0] + qiddot[1] * ti[1] +
                         qiddot[0] * ti[2]);
      Diddot[3] = 2.0 * (-2.0 * qiddot[3] * ti[0] - qiddot[0] * ti[1] +
                         qiddot[1] * ti[2]);

      Diddot[4] = 2.0 * (qiddot[3] * ti[0] - qiddot[1] * ti[2]);
      Diddot[5] = 2.0 * (qiddot[2] * ti[0] - 2.0 * qiddot[1] * ti[1] -
                         qiddot[0] * ti[2]);
      Diddot[6] = 2.0 * (qiddot[1] * ti[0] + qiddot[3] * ti[2]);
      Diddot[7] = 2.0 * (qiddot[0] * ti[0] - 2.0 * qiddot[3] * ti[1] +
                         qiddot[2] * ti[2]);

      Diddot[8] = 2.0 * (qiddot[1] * ti[1] - qiddot[2] * ti[0]);
      Diddot[9] = 2.0 * (qiddot[3] * ti[0] + qiddot[0] * ti[1] -
                         2.0 * qiddot[1] * ti[2]);
      Diddot[10] = 2.0 * (qiddot[3] * ti[1] - qiddot[0] * ti[0] -
                          2.0 * qiddot[2] * ti[2]);
      Diddot[11] = 2.0 * (qiddot[1] * ti[0] + qiddot[2] * ti[1]);
    }

    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    TacsScalar *r = NULL;
    if (res) {
      r = &res[offset];
    }
    TacsScalar *m = &mat[offset * size + offset];

    // Re-set the pointer
    const TacsScalar *Di = D;
    for (int i = 0; i < num_nodes; i++, Di += 12) {
      if (res) {
        r[0] += dd[0] * Di[0] + dd[1] * Di[4] + dd[2] * Di[8];
        r[1] += dd[0] * Di[1] + dd[1] * Di[5] + dd[2] * Di[9];
        r[2] += dd[0] * Di[2] + dd[1] * Di[6] + dd[2] * Di[10];
        r[3] += dd[0] * Di[3] + dd[1] * Di[7] + dd[2] * Di[11];

        r += vars_per_node;
      }

      const TacsScalar *Dj = D;
      const TacsScalar *Djdot = Ddot;
      const TacsScalar *Djddot = Dddot;
      for (int j = 0; j < num_nodes; j++, Dj += 12, Djdot += 12, Djddot += 12) {
        TacsScalar dfdq[12];
        for (int k = 0; k < 3; k++) {
          TacsScalar tmp[3];
          tmp[0] = d2d[dsize * (3 * i + k) + 3 * j] +
                   gamma * d2Tdotd[dsize * (3 * i + k) + 3 * j];
          tmp[1] = d2d[dsize * (3 * i + k) + 3 * j + 1] +
                   gamma * d2Tdotd[dsize * (3 * i + k) + 3 * j + 1];
          tmp[2] = d2d[dsize * (3 * i + k) + 3 * j + 2] +
                   gamma * d2Tdotd[dsize * (3 * i + k) + 3 * j + 2];

          dfdq[k] = Dj[0] * tmp[0] + Dj[4] * tmp[1] + Dj[8] * tmp[2];
          dfdq[3 + k] = Dj[1] * tmp[0] + Dj[5] * tmp[1] + Dj[9] * tmp[2];
          dfdq[6 + k] = Dj[2] * tmp[0] + Dj[6] * tmp[1] + Dj[10] * tmp[2];
          dfdq[9 + k] = Dj[3] * tmp[0] + Dj[7] * tmp[1] + Dj[11] * tmp[2];

          tmp[0] = beta * d2Tdotd[dsize * (3 * i + k) + 3 * j];
          tmp[1] = beta * d2Tdotd[dsize * (3 * i + k) + 3 * j + 1];
          tmp[2] = beta * d2Tdotd[dsize * (3 * i + k) + 3 * j + 2];

          dfdq[k] +=
              2.0 * (Djdot[0] * tmp[0] + Djdot[4] * tmp[1] + Djdot[8] * tmp[2]);
          dfdq[3 + k] +=
              2.0 * (Djdot[1] * tmp[0] + Djdot[5] * tmp[1] + Djdot[9] * tmp[2]);
          dfdq[6 + k] += 2.0 * (Djdot[2] * tmp[0] + Djdot[6] * tmp[1] +
                                Djdot[10] * tmp[2]);
          dfdq[9 + k] += 2.0 * (Djdot[3] * tmp[0] + Djdot[7] * tmp[1] +
                                Djdot[11] * tmp[2]);

          tmp[0] = alpha * d2Tdotd[dsize * (3 * i + k) + 3 * j];
          tmp[1] = alpha * d2Tdotd[dsize * (3 * i + k) + 3 * j + 1];
          tmp[2] = alpha * d2Tdotd[dsize * (3 * i + k) + 3 * j + 2];

          dfdq[k] +=
              Djddot[0] * tmp[0] + Djddot[4] * tmp[1] + Djddot[8] * tmp[2];
          dfdq[3 + k] +=
              Djddot[1] * tmp[0] + Djddot[5] * tmp[1] + Djddot[9] * tmp[2];
          dfdq[6 + k] +=
              Djddot[2] * tmp[0] + Djddot[6] * tmp[1] + Djddot[10] * tmp[2];
          dfdq[9 + k] +=
              Djddot[3] * tmp[0] + Djddot[7] * tmp[1] + Djddot[11] * tmp[2];
        }

        TacsScalar jac[16];
        for (int k = 0; k < 4; k++) {
          jac[k] = dfdq[3 * k] * Di[0] + dfdq[3 * k + 1] * Di[4] +
                   dfdq[3 * k + 2] * Di[8];
          jac[4 + k] = dfdq[3 * k] * Di[1] + dfdq[3 * k + 1] * Di[5] +
                       dfdq[3 * k + 2] * Di[9];
          jac[8 + k] = dfdq[3 * k] * Di[2] + dfdq[3 * k + 1] * Di[6] +
                       dfdq[3 * k + 2] * Di[10];
          jac[12 + k] = dfdq[3 * k] * Di[3] + dfdq[3 * k + 1] * Di[7] +
                        dfdq[3 * k + 2] * Di[11];
        }

        if (i == j) {
          TacsScalar tmp[3];
          tmp[0] = alpha * dd[0];
          tmp[1] = alpha * dd[1];
          tmp[2] = alpha * dd[2];

          jac[1] += 2.0 * (tmp[2] * t[1] - tmp[1] * t[2]);
          jac[2] += 2.0 * (tmp[0] * t[2] - tmp[2] * t[0]);
          jac[3] += 2.0 * (tmp[1] * t[0] - tmp[0] * t[1]);

          jac[4] += 2.0 * (tmp[2] * t[1] - tmp[1] * t[2]);
          jac[5] -= 4.0 * (tmp[1] * t[1] + tmp[2] * t[2]);
          jac[6] += 2.0 * (tmp[0] * t[1] + tmp[1] * t[0]);
          jac[7] += 2.0 * (tmp[0] * t[2] + tmp[2] * t[0]);

          jac[8] += 2.0 * (tmp[0] * t[2] - tmp[2] * t[0]);
          jac[9] += 2.0 * (tmp[0] * t[1] + tmp[1] * t[0]);
          jac[10] -= 4.0 * (tmp[0] * t[0] + tmp[2] * t[2]);
          jac[11] += 2.0 * (tmp[1] * t[2] + tmp[2] * t[1]);

          jac[12] += 2.0 * (tmp[1] * t[0] - tmp[0] * t[1]);
          jac[13] += 2.0 * (tmp[0] * t[2] + tmp[2] * t[0]);
          jac[14] += 2.0 * (tmp[1] * t[2] + tmp[2] * t[1]);
          jac[15] -= 4.0 * (tmp[0] * t[0] + tmp[1] * t[1]);
        }

        for (int ii = 0; ii < 4; ii++) {
          for (int jj = 0; jj < 4; jj++) {
            m[vars_per_node * j + ii * size + jj] += jac[4 * ii + jj];
          }
        }
      }

      q += vars_per_node;
      qdot += vars_per_node;

      dd += 3;
      t += 3;
      m += vars_per_node * size;
    }

    Di = D;
    const TacsScalar *Didot = Ddot;
    const TacsScalar *Diddot = Dddot;
    for (int i = 0; i < num_nodes; i++, Di += 12, Didot += 12, Diddot += 12) {
      for (int j = 0; j < num_nodes; j++) {
        TacsScalar dfdq[12], dfdqT[12];
        for (int k = 0; k < 3; k++) {
          TacsScalar tmp[3];
          tmp[0] = d2du[dsize * (3 * i) + 3 * j + k] +
                   gamma * d2Tdotu[dsize * (3 * i) + 3 * j + k];
          tmp[1] = d2du[dsize * (3 * i + 1) + 3 * j + k] +
                   gamma * d2Tdotu[dsize * (3 * i + 1) + 3 * j + k];
          tmp[2] = d2du[dsize * (3 * i + 2) + 3 * j + k] +
                   gamma * d2Tdotu[dsize * (3 * i + 2) + 3 * j + k];

          dfdq[k] = Di[0] * tmp[0] + Di[4] * tmp[1] + Di[8] * tmp[2];
          dfdq[3 + k] = Di[1] * tmp[0] + Di[5] * tmp[1] + Di[9] * tmp[2];
          dfdq[6 + k] = Di[2] * tmp[0] + Di[6] * tmp[1] + Di[10] * tmp[2];
          dfdq[9 + k] = Di[3] * tmp[0] + Di[7] * tmp[1] + Di[11] * tmp[2];

          tmp[0] = d2du[dsize * (3 * i) + 3 * j + k];
          tmp[1] = d2du[dsize * (3 * i + 1) + 3 * j + k];
          tmp[2] = d2du[dsize * (3 * i + 2) + 3 * j + k];

          dfdqT[k] = dfdq[k];
          dfdqT[3 + k] = dfdq[3 + k];
          dfdqT[6 + k] = dfdq[6 + k];
          dfdqT[9 + k] = dfdq[9 + k];

          tmp[0] = beta * d2Tdotu[dsize * (3 * i) + 3 * j + k];
          tmp[1] = beta * d2Tdotu[dsize * (3 * i + 1) + 3 * j + k];
          tmp[2] = beta * d2Tdotu[dsize * (3 * i + 2) + 3 * j + k];

          dfdqT[k] +=
              2.0 * (Didot[0] * tmp[0] + Didot[4] * tmp[1] + Didot[8] * tmp[2]);
          dfdqT[3 + k] +=
              2.0 * (Didot[1] * tmp[0] + Didot[5] * tmp[1] + Didot[9] * tmp[2]);
          dfdqT[6 + k] += 2.0 * (Didot[2] * tmp[0] + Didot[6] * tmp[1] +
                                 Didot[10] * tmp[2]);
          dfdqT[9 + k] += 2.0 * (Didot[3] * tmp[0] + Didot[7] * tmp[1] +
                                 Didot[11] * tmp[2]);

          tmp[0] = alpha * d2Tdotu[dsize * (3 * i) + 3 * j + k];
          tmp[1] = alpha * d2Tdotu[dsize * (3 * i + 1) + 3 * j + k];
          tmp[2] = alpha * d2Tdotu[dsize * (3 * i + 2) + 3 * j + k];

          dfdqT[k] +=
              Diddot[0] * tmp[0] + Diddot[4] * tmp[1] + Diddot[8] * tmp[2];
          dfdqT[3 + k] +=
              Diddot[1] * tmp[0] + Diddot[5] * tmp[1] + Diddot[9] * tmp[2];
          dfdqT[6 + k] +=
              Diddot[2] * tmp[0] + Diddot[6] * tmp[1] + Diddot[10] * tmp[2];
          dfdqT[9 + k] +=
              Diddot[3] * tmp[0] + Diddot[7] * tmp[1] + Diddot[11] * tmp[2];
        }

        for (int ii = 0; ii < 4; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            int index = (vars_per_node * i + ii + offset) * size +
                        vars_per_node * j + jj;

            mat[index] += dfdq[3 * ii + jj];
          }
        }

        for (int ii = 0; ii < 4; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            int index = (vars_per_node * j + jj) * size + vars_per_node * i +
                        ii + offset;

            mat[index] += dfdqT[3 * ii + jj];
          }
        }
      }
    }
  }

  /*
    Add the director contributions to the derivative of the normal

    Add the adjoint sensitivity of the reference normal (dt) based on
    the adjoint sensitivity of the director (dd).

    Given that the parametrization is d = (C^{T}(q) - I) * t, compute

    dt += d(dd^{T} * d)/dt = dd^{T} * C^{T}(q) = C(q) * dd
    .   = (0.5*q^{x} - 1)*q^{x} * dpsi

    @param vars The full variable vector
    @param t The reference directions
    @param dd The adjoint sensitivities w.r.t. the director
    @param dt The adjoint sensitivity w.r.t. the reference directions
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorRefNormalSens(const TacsScalar vars[],
                                       const TacsScalar t[],
                                       const TacsScalar dd[], TacsScalar dt[]) {
    const TacsScalar *q = &vars[offset];

    for (int i = 0; i < num_nodes; i++) {
      TacsScalar C[9];
      C[0] = -2.0 * (q[2] * q[2] + q[3] * q[3]);
      C[1] = 2.0 * (q[1] * q[2] + q[3] * q[0]);
      C[2] = 2.0 * (q[1] * q[3] - q[2] * q[0]);

      C[3] = 2.0 * (q[2] * q[1] - q[3] * q[0]);
      C[4] = -2.0 * (q[1] * q[1] + q[3] * q[3]);
      C[5] = 2.0 * (q[2] * q[3] + q[1] * q[0]);

      C[6] = 2.0 * (q[3] * q[1] + q[2] * q[0]);
      C[7] = 2.0 * (q[3] * q[2] - q[1] * q[0]);
      C[8] = -2.0 * (q[1] * q[1] + q[2] * q[2]);

      dt[0] += C[0] * dd[0] + C[1] * dd[1] + C[2] * dd[2];
      dt[1] += C[3] * dd[0] + C[4] * dd[1] + C[5] * dd[2];
      dt[2] += C[6] * dd[0] + C[7] * dd[1] + C[8] * dd[2];

      t += 3;
      dd += 3;
      dt += 3;
      q += vars_per_node;
    }
  }

  /*
    Add the director contributions to the derivative of the normal

    Add the adjoint sensitivity of the reference normal (dt) based on
    the adjoint sensitivity of the director (dd) and the sensitivity
    of the derivative field (ddadj).

    Given that the parametrization is d = (C^{T}(q) - I) * t and the field
    dpsi = d(d)/dq^{T} * psi, compute

    dt += d(dd^{T}d)/dn = dd^{T} * (C^{T}(q) - I) = (C(q) - I) * dd

    dt += d(ddpsi^{T} * dpsi)/dt = [ d(C(q))/dq * psi ] * ddpsi

    @param vars The full variable vector
    @param psi The full variable vector derivative
    @param t The reference directions
    @param dd The adjoint sensitivities w.r.t. the director
    @param ddpsi The adjoint sensitivities w.r.t. the director derivative
    @param dt The adjoint sensitivity w.r.t. the reference directions
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorRefNormalSens(
      const TacsScalar vars[], const TacsScalar psi[], const TacsScalar t[],
      const TacsScalar dd[], const TacsScalar ddpsi[], TacsScalar dt[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qpsi = &psi[offset];

    for (int i = 0; i < num_nodes; i++) {
      TacsScalar C[9], Cd[9];
      C[0] = -2.0 * (q[2] * q[2] + q[3] * q[3]);
      C[1] = 2.0 * (q[1] * q[2] + q[3] * q[0]);
      C[2] = 2.0 * (q[1] * q[3] - q[2] * q[0]);

      C[3] = 2.0 * (q[2] * q[1] - q[3] * q[0]);
      C[4] = -2.0 * (q[1] * q[1] + q[3] * q[3]);
      C[5] = 2.0 * (q[2] * q[3] + q[1] * q[0]);

      C[6] = 2.0 * (q[3] * q[1] + q[2] * q[0]);
      C[7] = 2.0 * (q[3] * q[2] - q[1] * q[0]);
      C[8] = -2.0 * (q[1] * q[1] + q[2] * q[2]);

      dt[0] += C[0] * dd[0] + C[1] * dd[1] + C[2] * dd[2];
      dt[1] += C[3] * dd[0] + C[4] * dd[1] + C[5] * dd[2];
      dt[2] += C[6] * dd[0] + C[7] * dd[1] + C[8] * dd[2];

      Cd[0] = -4.0 * (q[2] * qpsi[2] + q[3] * qpsi[3]);
      Cd[1] = 2.0 * (q[1] * qpsi[2] + q[3] * qpsi[0] + qpsi[1] * q[2] +
                     qpsi[3] * q[0]);
      Cd[2] = 2.0 * (q[1] * qpsi[3] - q[2] * qpsi[0] + qpsi[1] * q[3] -
                     qpsi[2] * q[0]);

      Cd[3] = 2.0 * (q[2] * qpsi[1] - q[3] * qpsi[0] + qpsi[2] * q[1] -
                     qpsi[3] * q[0]);
      Cd[4] = -4.0 * (q[1] * qpsi[1] + q[3] * qpsi[3]);
      Cd[5] = 2.0 * (q[2] * qpsi[3] + q[1] * qpsi[0] + qpsi[2] * q[3] +
                     qpsi[1] * q[0]);

      Cd[6] = 2.0 * (q[3] * qpsi[1] + q[2] * qpsi[0] + qpsi[3] * q[1] +
                     qpsi[2] * q[0]);
      Cd[7] = 2.0 * (q[3] * qpsi[2] - q[1] * qpsi[0] + qpsi[3] * q[2] -
                     qpsi[1] * q[0]);
      Cd[8] = -4.0 * (q[1] * qpsi[1] + q[2] * qpsi[2]);

      dt[0] += Cd[0] * ddpsi[0] + Cd[1] * ddpsi[1] + Cd[2] * ddpsi[2];
      dt[1] += Cd[3] * ddpsi[0] + Cd[4] * ddpsi[1] + Cd[5] * ddpsi[2];
      dt[2] += Cd[6] * ddpsi[0] + Cd[7] * ddpsi[1] + Cd[8] * ddpsi[2];

      t += 3;
      dd += 3;
      ddpsi += 3;
      dt += 3;
      q += vars_per_node;
      qpsi += vars_per_node;
    }
  }

  static TacsScalar evalDrillStrain(const TacsScalar u0x[],
                                    const TacsScalar Ct[]) {
    // e2^{T}*Ct*(e1 + u_{,x}*e1) - e1^{T}*Ct*(e2 + u_{,x}*e2)
    return ((Ct[3] * (1.0 + u0x[0]) + Ct[4] * u0x[3] + Ct[5] * u0x[6]) -
            (Ct[0] * u0x[1] + Ct[1] * (1.0 + u0x[4]) + Ct[2] * u0x[7]));
  }

  static void evalDrillStrainSens(TacsScalar scale, const TacsScalar u0x[],
                                  const TacsScalar Ct[], TacsScalar du0x[],
                                  TacsScalar dCt[]) {
    // Derivative with respect to u0x
    du0x[0] = scale * Ct[3];
    du0x[1] = -scale * Ct[0];
    du0x[2] = 0.0;
    du0x[3] = scale * Ct[4];
    du0x[4] = -scale * Ct[1];
    du0x[5] = 0.0;
    du0x[6] = scale * Ct[5];
    du0x[7] = -scale * Ct[2];
    du0x[8] = 0.0;

    dCt[0] = -scale * u0x[1];
    dCt[1] = -scale * (1.0 + u0x[4]);
    dCt[2] = -scale * u0x[7];
    dCt[3] = scale * (1.0 + u0x[0]);
    dCt[4] = scale * u0x[3];
    dCt[5] = scale * u0x[6];
    dCt[6] = 0.0;
    dCt[7] = 0.0;
    dCt[8] = 0.0;
  }

  static TacsScalar evalDrillStrainDeriv(const TacsScalar u0x[],
                                         const TacsScalar Ct[],
                                         const TacsScalar u0xd[],
                                         const TacsScalar Ctd[],
                                         TacsScalar *ed) {
    // e2^{T}*Ct*(e1 + u_{,x}*e1) - e1^{T}*Ct*(e2 + u_{,x}*e2)
    *ed = ((Ctd[3] * (1.0 + u0x[0]) + Ctd[4] * u0x[3] + Ctd[5] * u0x[6]) +
           (Ct[3] * u0xd[0] + Ct[4] * u0xd[3] + Ct[5] * u0xd[6])) -
          ((Ctd[0] * u0x[1] + Ctd[1] * (1.0 + u0x[4]) + Ctd[2] * u0x[7]) +
           (Ct[0] * u0xd[1] + Ct[1] * u0xd[4] + Ct[2] * u0xd[7]));

    return ((Ct[3] * (1.0 + u0x[0]) + Ct[4] * u0x[3] + Ct[5] * u0x[6]) -
            (Ct[0] * u0x[1] + Ct[1] * (1.0 + u0x[4]) + Ct[2] * u0x[7]));
  }

  static void evalDrillStrainHessian(TacsScalar det, const TacsScalar u0x[],
                                     const TacsScalar Ct[], TacsScalar d2u0x[],
                                     TacsScalar d2Ct[], TacsScalar d2Ctu0x[]) {
    memset(d2u0x, 0, 81 * sizeof(TacsScalar));
    memset(d2Ct, 0, 81 * sizeof(TacsScalar));
    memset(d2Ctu0x, 0, 81 * sizeof(TacsScalar));

    // d2Ctu0x[0] = det*(0);
    // d2Ctu0x[1] = det*(-1);
    // d2Ctu0x[2] = det*(0);
    // d2Ctu0x[3] = det*(0);
    // d2Ctu0x[4] = det*(0);
    // d2Ctu0x[5] = det*(0);
    // d2Ctu0x[6] = det*(0);
    // d2Ctu0x[7] = det*(0);
    // d2Ctu0x[8] = det*(0);
    // d2Ctu0x[9] = det*(0);
    // d2Ctu0x[10] = det*(0);
    // d2Ctu0x[11] = det*(0);
    // d2Ctu0x[12] = det*(0);
    // d2Ctu0x[13] = det*(-1);
    // d2Ctu0x[14] = det*(0);
    // d2Ctu0x[15] = det*(0);
    // d2Ctu0x[16] = det*(0);
    // d2Ctu0x[17] = det*(0);
    // d2Ctu0x[18] = det*(0);
    // d2Ctu0x[19] = det*(0);
    // d2Ctu0x[20] = det*(0);
    // d2Ctu0x[21] = det*(0);
    // d2Ctu0x[22] = det*(0);
    // d2Ctu0x[23] = det*(0);
    // d2Ctu0x[24] = det*(0);
    // d2Ctu0x[25] = det*(-1);
    // d2Ctu0x[26] = det*(0);
    // d2Ctu0x[27] = det*(1);
    // d2Ctu0x[28] = det*(0);
    // d2Ctu0x[29] = det*(0);
    // d2Ctu0x[30] = det*(0);
    // d2Ctu0x[31] = det*(0);
    // d2Ctu0x[32] = det*(0);
    // d2Ctu0x[33] = det*(0);
    // d2Ctu0x[34] = det*(0);
    // d2Ctu0x[35] = det*(0);
    // d2Ctu0x[36] = det*(0);
    // d2Ctu0x[37] = det*(0);
    // d2Ctu0x[38] = det*(0);
    // d2Ctu0x[39] = det*(1);
    // d2Ctu0x[40] = det*(0);
    // d2Ctu0x[41] = det*(0);
    // d2Ctu0x[42] = det*(0);
    // d2Ctu0x[43] = det*(0);
    // d2Ctu0x[44] = det*(0);
    // d2Ctu0x[45] = det*(0);
    // d2Ctu0x[46] = det*(0);
    // d2Ctu0x[47] = det*(0);
    // d2Ctu0x[48] = det*(0);
    // d2Ctu0x[49] = det*(0);
    // d2Ctu0x[50] = det*(0);
    // d2Ctu0x[51] = det*(1);
    // d2Ctu0x[52] = det*(0);
    // d2Ctu0x[53] = det*(0);
    // d2Ctu0x[54] = det*(0);
    // d2Ctu0x[55] = det*(0);
    // d2Ctu0x[56] = det*(0);
    // d2Ctu0x[57] = det*(0);
    // d2Ctu0x[58] = det*(0);
    // d2Ctu0x[59] = det*(0);
    // d2Ctu0x[60] = det*(0);
    // d2Ctu0x[61] = det*(0);
    // d2Ctu0x[62] = det*(0);
    // d2Ctu0x[63] = det*(0);
    // d2Ctu0x[64] = det*(0);
    // d2Ctu0x[65] = det*(0);
    // d2Ctu0x[66] = det*(0);
    // d2Ctu0x[67] = det*(0);
    // d2Ctu0x[68] = det*(0);
    // d2Ctu0x[69] = det*(0);
    // d2Ctu0x[70] = det*(0);
    // d2Ctu0x[71] = det*(0);
    // d2Ctu0x[72] = det*(0);
    // d2Ctu0x[73] = det*(0);
    // d2Ctu0x[74] = det*(0);
    // d2Ctu0x[75] = det*(0);
    // d2Ctu0x[76] = det*(0);
    // d2Ctu0x[77] = det*(0);
    // d2Ctu0x[78] = det*(0);
    // d2Ctu0x[79] = det*(0);
    // d2Ctu0x[80] = det*(0);

    // d2Ct[0] = scale*(drill*u0x[1]*u0x[1]);
    // d2Ct[1] = scale*(drill*u0x[1]*(u0x[4] + 1.0));
    // d2Ct[2] = scale*(drill*u0x[1]*u0x[7]);
    // d2Ct[3] = scale*(-drill*u0x[1]*(u0x[0] + 1.0));
    // d2Ct[4] = scale*(-drill*u0x[1]*u0x[3]);
    // d2Ct[5] = scale*(-drill*u0x[1]*u0x[6]);
    // d2Ct[6] = 0.0;
    // d2Ct[7] = 0.0;
    // d2Ct[8] = 0.0;

    // d2Ct[9] = scale*(drill*u0x[1]*(u0x[4] + 1.0));
    // d2Ct[10] = scale*(0.5*drill*(u0x[4] + 1.0)*(2.0*u0x[4] + 2.0));
    // d2Ct[11] = scale*(drill*u0x[7]*(u0x[4] + 1.0));
    // d2Ct[12] = scale*(-0.5*drill*(u0x[0] + 1.0)*(2.0*u0x[4] + 2.0));
    // d2Ct[13] = scale*(-drill*u0x[3]*(u0x[4] + 1.0));
    // d2Ct[14] = scale*(-drill*u0x[6]*(u0x[4] + 1.0));
    // d2Ct[15] = 0.0;
    // d2Ct[16] = 0.0;
    // d2Ct[17] = 0.0;

    // d2Ct[18] = scale*(drill*u0x[1]*u0x[7]);
    // d2Ct[19] = scale*(drill*u0x[7]*(u0x[4] + 1.0));
    // d2Ct[20] = scale*(drill*u0x[7]*u0x[7]);
    // d2Ct[21] = scale*(-drill*u0x[7]*(u0x[0] + 1.0));
    // d2Ct[22] = scale*(-drill*u0x[3]*u0x[7]);
    // d2Ct[23] = scale*(-drill*u0x[6]*u0x[7]);
    // d2Ct[24] = 0.0;
    // d2Ct[25] = 0.0;
    // d2Ct[26] = 0.0;

    // d2Ct[27] = scale*(-drill*u0x[1]*(u0x[0] + 1.0));
    // d2Ct[28] = scale*(-0.5*drill*(2.0*u0x[0] + 2.0)*(u0x[4] + 1.0));
    // d2Ct[29] = scale*(-drill*u0x[7]*(u0x[0] + 1.0));
    // d2Ct[30] = scale*(0.5*drill*(u0x[0] + 1.0)*(2.0*u0x[0] + 2.0));
    // d2Ct[31] = scale*(drill*u0x[3]*(u0x[0] + 1.0));
    // d2Ct[32] = scale*(drill*u0x[6]*(u0x[0] + 1.0));
    // d2Ct[33] = 0.0;
    // d2Ct[34] = 0.0;
    // d2Ct[35] = 0.0;

    // d2Ct[36] = scale*(-drill*u0x[1]*u0x[3]);
    // d2Ct[37] = scale*(-drill*u0x[3]*(u0x[4] + 1.0));
    // d2Ct[38] = scale*(-drill*u0x[3]*u0x[7]);
    // d2Ct[39] = scale*(drill*u0x[3]*(u0x[0] + 1.0));
    // d2Ct[40] = scale*(drill*u0x[3]*u0x[3]);
    // d2Ct[41] = scale*(drill*u0x[3]*u0x[6]);
    // d2Ct[42] = 0.0;
    // d2Ct[43] = 0.0;
    // d2Ct[44] = 0.0;

    // d2Ct[45] = scale*(-drill*u0x[1]*u0x[6]);
    // d2Ct[46] = scale*(-drill*u0x[6]*(u0x[4] + 1.0));
    // d2Ct[47] = scale*(-drill*u0x[6]*u0x[7]);
    // d2Ct[48] = scale*(drill*u0x[6]*(u0x[0] + 1.0));
    // d2Ct[49] = scale*(drill*u0x[3]*u0x[6]);
    // d2Ct[50] = scale*(drill*u0x[6]*u0x[6]);
    // d2Ct[51] = 0.0;
    // d2Ct[52] = 0.0;
    // d2Ct[53] = 0.0;

    // d2Ct[54] = 0.0;
    // d2Ct[55] = 0.0;
    // d2Ct[56] = 0.0;
    // d2Ct[57] = 0.0;
    // d2Ct[58] = 0.0;
    // d2Ct[59] = 0.0;
    // d2Ct[60] = 0.0;
    // d2Ct[61] = 0.0;
    // d2Ct[62] = 0.0;

    // d2Ct[63] = 0.0;
    // d2Ct[64] = 0.0;
    // d2Ct[65] = 0.0;
    // d2Ct[66] = 0.0;
    // d2Ct[67] = 0.0;
    // d2Ct[68] = 0.0;
    // d2Ct[69] = 0.0;
    // d2Ct[70] = 0.0;
    // d2Ct[71] = 0.0;

    // d2Ct[72] = 0.0;
    // d2Ct[73] = 0.0;
    // d2Ct[74] = 0.0;
    // d2Ct[75] = 0.0;
    // d2Ct[76] = 0.0;
    // d2Ct[77] = 0.0;
    // d2Ct[78] = 0.0;
    // d2Ct[79] = 0.0;
    // d2Ct[80] = 0.0;

    // d2Ctu0x[0] = scale*(-Ct[3]*drill*u0x[1]);
    // d2Ctu0x[1] = scale*(drill*(2.0*Ct[0]*u0x[1] + Ct[1]*(u0x[4] + 1.0) +
    // Ct[2]*u0x[7] - Ct[3]*(u0x[0] + 1.0) - Ct[4]*u0x[3] - Ct[5]*u0x[6]));
    // d2Ctu0x[2] = 0.0;
    // d2Ctu0x[3] = scale*(-Ct[4]*drill*u0x[1]);
    // d2Ctu0x[4] = scale*(Ct[1]*drill*u0x[1]);
    // d2Ctu0x[5] = 0.0;
    // d2Ctu0x[6] = scale*(-Ct[5]*drill*u0x[1]);
    // d2Ctu0x[7] = scale*(Ct[2]*drill*u0x[1]);
    // d2Ctu0x[8] = 0.0;

    // d2Ctu0x[9] = scale*(-Ct[3]*drill*(u0x[4] + 1.0));
    // d2Ctu0x[10] = scale*(Ct[0]*drill*(u0x[4] + 1.0));
    // d2Ctu0x[11] = 0.0;
    // d2Ctu0x[12] = scale*(-Ct[4]*drill*(u0x[4] + 1.0));
    // d2Ctu0x[13] = scale*(drill*(Ct[0]*u0x[1] + Ct[1]*(u0x[4] + 1.0) +
    // 0.5*Ct[1]*(2.0*u0x[4] + 2.0) + Ct[2]*u0x[7] - Ct[3]*(u0x[0] + 1.0) -
    // Ct[4]*u0x[3] - Ct[5]*u0x[6])); d2Ctu0x[14] = 0.0; d2Ctu0x[15] =
    // scale*(-Ct[5]*drill*(u0x[4] + 1.0)); d2Ctu0x[16] =
    // scale*(Ct[2]*drill*(u0x[4] + 1.0)); d2Ctu0x[17] = 0.0;

    // d2Ctu0x[18] = scale*(-Ct[3]*drill*u0x[7]);
    // d2Ctu0x[19] = scale*(Ct[0]*drill*u0x[7]);
    // d2Ctu0x[20] = 0.0;
    // d2Ctu0x[21] = scale*(-Ct[4]*drill*u0x[7]);
    // d2Ctu0x[22] = scale*(Ct[1]*drill*u0x[7]);
    // d2Ctu0x[23] = 0.0;
    // d2Ctu0x[24] = scale*(-Ct[5]*drill*u0x[7]);
    // d2Ctu0x[25] = scale*(drill*(Ct[0]*u0x[1] + Ct[1]*(u0x[4] + 1.0)
    // + 2.0*Ct[2]*u0x[7] - Ct[3]*(u0x[0] + 1.0) - Ct[4]*u0x[3] -
    // Ct[5]*u0x[6])); d2Ctu0x[26] = 0.0;

    // d2Ctu0x[27] = scale*(drill*(-Ct[0]*u0x[1] - Ct[1]*(u0x[4] + 1.0) -
    // Ct[2]*u0x[7] + Ct[3]*(u0x[0] + 1.0) + 0.5*Ct[3]*(2.0*u0x[0] + 2.0) +
    // Ct[4]*u0x[3] + Ct[5]*u0x[6])); d2Ctu0x[28] = scale*(-Ct[0]*drill*(u0x[0]
    // + 1.0)); d2Ctu0x[29] = 0.0; d2Ctu0x[30] = scale*(Ct[4]*drill*(u0x[0]
    // + 1.0)); d2Ctu0x[31] = scale*(-Ct[1]*drill*(u0x[0] + 1.0)); d2Ctu0x[32] =
    // 0.0; d2Ctu0x[33] = scale*(Ct[5]*drill*(u0x[0] + 1.0)); d2Ctu0x[34] =
    // scale*(-Ct[2]*drill*(u0x[0] + 1.0)); d2Ctu0x[35] = 0.0;

    // d2Ctu0x[36] = scale*(Ct[3]*drill*u0x[3]);
    // d2Ctu0x[37] = scale*(-Ct[0]*drill*u0x[3]);
    // d2Ctu0x[38] = 0.0;
    // d2Ctu0x[39] = scale*(drill*(-Ct[0]*u0x[1] - Ct[1]*(u0x[4] + 1.0) -
    // Ct[2]*u0x[7] + Ct[3]*(u0x[0] + 1.0) + 2.0*Ct[4]*u0x[3] + Ct[5]*u0x[6]));
    // d2Ctu0x[40] = scale*(-Ct[1]*drill*u0x[3]);
    // d2Ctu0x[41] = 0.0;
    // d2Ctu0x[42] = scale*(Ct[5]*drill*u0x[3]);
    // d2Ctu0x[43] = scale*(-Ct[2]*drill*u0x[3]);
    // d2Ctu0x[44] = 0.0;

    // d2Ctu0x[45] = scale*(Ct[3]*drill*u0x[6]);
    // d2Ctu0x[46] = scale*(-Ct[0]*drill*u0x[6]);
    // d2Ctu0x[47] = 0.0;
    // d2Ctu0x[48] = scale*(Ct[4]*drill*u0x[6]);
    // d2Ctu0x[49] = scale*(-Ct[1]*drill*u0x[6]);
    // d2Ctu0x[50] = 0.0;
    // d2Ctu0x[51] = scale*(drill*(-Ct[0]*u0x[1] - Ct[1]*(u0x[4] + 1.0) -
    // Ct[2]*u0x[7] + Ct[3]*(u0x[0] + 1.0) + Ct[4]*u0x[3] + 2.0*Ct[5]*u0x[6]));
    // d2Ctu0x[52] = scale*(-Ct[2]*drill*u0x[6]);
    // d2Ctu0x[53] = 0.0;

    // d2Ctu0x[54] = 0.0;
    // d2Ctu0x[55] = 0.0;
    // d2Ctu0x[56] = 0.0;
    // d2Ctu0x[57] = 0.0;
    // d2Ctu0x[58] = 0.0;
    // d2Ctu0x[59] = 0.0;
    // d2Ctu0x[60] = 0.0;
    // d2Ctu0x[61] = 0.0;

    // d2Ctu0x[62] = 0.0;
    // d2Ctu0x[63] = 0.0;
    // d2Ctu0x[64] = 0.0;
    // d2Ctu0x[65] = 0.0;
    // d2Ctu0x[66] = 0.0;
    // d2Ctu0x[67] = 0.0;
    // d2Ctu0x[68] = 0.0;
    // d2Ctu0x[69] = 0.0;
    // d2Ctu0x[70] = 0.0;
    // d2Ctu0x[71] = 0.0;

    // d2Ctu0x[72] = 0.0;
    // d2Ctu0x[73] = 0.0;
    // d2Ctu0x[74] = 0.0;
    // d2Ctu0x[75] = 0.0;
    // d2Ctu0x[76] = 0.0;
    // d2Ctu0x[77] = 0.0;
    // d2Ctu0x[78] = 0.0;
    // d2Ctu0x[79] = 0.0;
    // d2Ctu0x[80] = 0.0;
  }
};

template <int vars_per_node, int offset, int num_nodes, class director>
int TacsTestDirector(double dh = 1e-7, int test_print_level = 2,
                     double test_fail_atol = 1e-5,
                     double test_fail_rtol = 1e-5) {
  const int size = vars_per_node * num_nodes;
  const int dsize = 3 * num_nodes;
  const int csize = 9 * num_nodes;

  // Generate random arrays for the state variables and their time derivatives
  TacsScalar vars[size], dvars[size], ddvars[size];
  TacsGenerateRandomArray(vars, size);
  TacsGenerateRandomArray(dvars, size);
  TacsGenerateRandomArray(ddvars, size);

  // Compute/normalize the normals
  TacsScalar t[dsize];
  TacsGenerateRandomArray(t, dsize);
  for (int i = 0; i < num_nodes; i++) {
    TacsScalar tnrm = sqrt(vec3Dot(&t[3 * i], &t[3 * i]));
    vec3Scale(1.0 / tnrm, &t[3 * i]);
  }

  // Random perturbation for the variables
  TacsScalar varsd[size];
  TacsGenerateRandomArray(varsd, size);

  // Create random arrays for testing the residual and Jacobian
  TacsScalar dC[csize], d2C[csize * csize];
  TacsGenerateRandomArray(dC, csize);
  TacsGenerateRandomArray(d2C, csize * csize);
  for (int i = 0; i < csize; i++) {
    for (int j = 0; j < i; j++) {
      d2C[j + i * csize] = d2C[i + j * csize];
    }
  }

  // Compute the rotation matrices
  TacsScalar C[csize];
  director::template computeRotationMat<vars_per_node, offset, num_nodes>(vars,
                                                                          C);

  // Compute the Jacobian and residual
  TacsScalar res[size], mat[size * size];
  memset(res, 0, size * sizeof(TacsScalar));
  memset(mat, 0, size * size * sizeof(TacsScalar));
  director::template addRotationMatJacobian<vars_per_node, offset, num_nodes>(
      1.0, vars, dC, d2C, res, mat);

  // Verify the implementation of the residual
  TacsScalar fd[size], C0 = 0.0;
  for (int i = 0; i < csize; i++) {
    C0 += dC[i] * C[i];
  }

  for (int k = 0; k < size; k++) {
    TacsScalar varst[size];
    memcpy(varst, vars, size * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    varst[k] = vars[k] + TacsScalar(0.0, dh);
#else
    varst[k] = vars[k] + dh;
#endif  // TACS_USE_COMPLEX

    TacsScalar Ct[csize];
    director::template computeRotationMat<vars_per_node, offset, num_nodes>(
        varst, Ct);

    TacsScalar C1 = 0.0;
    for (int i = 0; i < csize; i++) {
      C1 += dC[i] * Ct[i];
    }

#ifdef TACS_USE_COMPLEX
    fd[k] = TacsImagPart(C1) / dh;
#else
    fd[k] = (C1 - C0) / dh;
#endif  // TACS_USE_COMPLEX
  }

  // Variables to store the max error and indices
  int max_err_index, max_rel_index;
  double max_err, max_rel;

  // Keep track of the failure flag
  int fail = 0;

  // Compute the error
  max_err = TacsGetMaxError(res, fd, size, &max_err_index);
  max_rel = TacsGetMaxRelError(res, fd, size, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the director residual implementation\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "res", res, fd, size);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute the derivative of the rotation matrix
  TacsScalar Cd[csize];
  director::template computeRotationMatDeriv<vars_per_node, offset, num_nodes>(
      vars, varsd, C, Cd);

  TacsScalar q[size];
  for (int k = 0; k < size; k++) {
#ifdef TACS_USE_COMPLEX
    q[k] = vars[k] + varsd[k] * TacsScalar(0.0, dh);
#else
    q[k] = vars[k] + varsd[k] * dh;
#endif  // TACS_USE_COMPLEX
  }

  TacsScalar Ctemp[csize];
  director::template computeRotationMat<vars_per_node, offset, num_nodes>(
      q, Ctemp);

  TacsScalar fdC[csize];
  for (int k = 0; k < csize; k++) {
#ifdef TACS_USE_COMPLEX
    fdC[k] = TacsImagPart(Ctemp[k]) / dh;
#else
    fdC[k] = (Ctemp[k] - C[k]) / dh;
#endif  // TACS_USE_COMPLEX
  }

  // Compute the error
  max_err = TacsGetMaxError(Cd, fdC, csize, &max_err_index);
  max_rel = TacsGetMaxRelError(Cd, fdC, csize, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative of the rotation matrix\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "Cd", Cd, fdC, csize);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute the derivative of the rotation matrix residual
  TacsScalar fdmat[size * size];
  for (int k = 0; k < size; k++) {
    TacsScalar varst[size];
    memcpy(varst, vars, size * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    varst[k] = vars[k] + TacsScalar(0.0, dh);
#else
    varst[k] = vars[k] + dh;
#endif  // TACS_USE_COMPLEX

    TacsScalar Ct[csize];
    director::template computeRotationMat<vars_per_node, offset, num_nodes>(
        varst, Ct);

    // Add the contributions from the
    TacsScalar dCt[csize];
    for (int i = 0; i < csize; i++) {
      dCt[i] = dC[i];

      for (int j = 0; j < csize; j++) {
        dCt[i] += d2C[j + i * csize] * (Ct[j] - C[j]);
      }
    }

    TacsScalar rest[size];
    memset(rest, 0, size * sizeof(TacsScalar));
    director::template addRotationMatResidual<vars_per_node, offset, num_nodes>(
        varst, dCt, rest);

    for (int j = 0; j < size; j++) {
#ifdef TACS_USE_COMPLEX
      fdmat[k + size * j] = TacsImagPart(rest[j]) / dh;
#else
      fdmat[k + size * j] = (rest[j] - res[j]) / dh;
#endif  // TACS_USE_COMPLEX
    }
  }

  // Compute the error
  max_err = TacsGetMaxError(mat, fdmat, size * size, &max_err_index);
  max_rel = TacsGetMaxRelError(mat, fdmat, size * size, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr,
            "Testing the derivative of the rotation matrix w.r.t. vars\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "mat", mat, fdmat, size * size);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Check for consistency between the director and C
  director::template computeRotationMat<vars_per_node, offset, num_nodes>(vars,
                                                                          C);

  TacsScalar d[dsize], ddot[dsize];
  director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
      vars, dvars, t, d, ddot);

  TacsScalar dcal[dsize];
  for (int i = 0; i < num_nodes; i++) {
    const TacsScalar *C0 = &C[9 * i];
    TacsScalar Qt[9];
    Qt[0] = C0[0] - 1.0;
    Qt[1] = C0[1];
    Qt[2] = C0[2];
    Qt[3] = C0[3];
    Qt[4] = C0[4] - 1.0;
    Qt[5] = C0[5];
    Qt[6] = C0[6];
    Qt[7] = C0[7];
    Qt[8] = C0[8] - 1.0;

    mat3x3MultTrans(Qt, &t[3 * i], &dcal[3 * i]);
  }

  // Compute the error
  max_err = TacsGetMaxError(d, dcal, dsize, &max_err_index);
  max_rel = TacsGetMaxRelError(d, dcal, dsize, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the consistency of the director\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d", d, dcal, dsize);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Test the implementation of the director time derivative
  TacsScalar dddot[dsize];
  director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
      vars, dvars, ddvars, t, d, ddot, dddot);

  TacsScalar fddot[dsize], varst[size];
  for (int k = 0; k < size; k++) {
#ifdef TACS_USE_COMPLEX
    varst[k] = vars[k] + dvars[k] * TacsScalar(0.0, dh);
#else
    varst[k] = vars[k] + dh * dvars[k];
#endif  // TACS_USE_COMPLEX
  }

  TacsScalar dt[dsize], dtdot[dsize];
  director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
      varst, dvars, t, dt, dtdot);

  for (int k = 0; k < dsize; k++) {
#ifdef TACS_USE_COMPLEX
    fddot[k] = TacsImagPart(dt[k]) / dh;
#else
    fddot[k] = (dt[k] - d[k]) / dh;
#endif  // TACS_USE_COMPLEX
  }

  // Compute the error
  max_err = TacsGetMaxError(ddot, fddot, dsize, &max_err_index);
  max_rel = TacsGetMaxRelError(ddot, fddot, dsize, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the time derivative of the director\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "ddot", ddot, fddot, dsize);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  TacsScalar fdddot[dsize], dvarst[size];
  for (int k = 0; k < size; k++) {
#ifdef TACS_USE_COMPLEX
    varst[k] = vars[k] + dvars[k] * TacsScalar(0.0, dh) +
               0.5 * ddvars[k] * TacsScalar(0.0, dh * dh);
    dvarst[k] = dvars[k] + ddvars[k] * TacsScalar(0.0, dh);
#else
    varst[k] = vars[k] + dh * dvars[k] + 0.5 * dh * dh * ddvars[k];
    dvarst[k] = dvars[k] + dh * ddvars[k];
#endif  // TACS_USE_COMPLEX
  }

  director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
      varst, dvarst, t, dt, dtdot);

  for (int k = 0; k < dsize; k++) {
#ifdef TACS_USE_COMPLEX
    fdddot[k] = TacsImagPart(dtdot[k]) / dh;
#else
    fdddot[k] = (dtdot[k] - ddot[k]) / dh;
#endif  // TACS_USE_COMPLEX
  }

  // Compute the error
  max_err = TacsGetMaxError(dddot, fdddot, dsize, &max_err_index);
  max_rel = TacsGetMaxRelError(dddot, fdddot, dsize, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the second time derivative of the director\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "dddot", dddot, fdddot, dsize);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Set up the data for the director
  TacsScalar dTdot[dsize], dd[dsize];
  TacsScalar d2Tdotd[dsize * dsize], d2Tdotu[dsize * dsize];
  TacsScalar d2d[dsize * dsize], d2du[dsize * dsize];

  TacsScalar alpha = 0.0, beta = 0.0, gamma = 0.0;
  TacsGenerateRandomArray(&alpha, 1);
  TacsGenerateRandomArray(&beta, 1);
  TacsGenerateRandomArray(&gamma, 1);
  TacsGenerateRandomArray(dTdot, dsize);
  TacsGenerateRandomArray(dd, dsize);
  TacsGenerateRandomArray(d2Tdotd, dsize * dsize);
  TacsGenerateRandomArray(d2Tdotu, dsize * dsize);
  TacsGenerateRandomArray(d2d, dsize * dsize);
  TacsGenerateRandomArray(d2du, dsize * dsize);

  // Compute the director rates
  director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
      vars, dvars, ddvars, t, d, ddot, dddot);

  // Compute the director Jacobian matrix
  memset(res, 0, size * sizeof(TacsScalar));
  memset(mat, 0, size * size * sizeof(TacsScalar));

  for (int i = 0; i < dsize * dsize; i++) {
    d2d[i] *= alpha;
    d2du[i] *= alpha;
  }

  director::template addDirectorJacobian<vars_per_node, offset, num_nodes>(
      alpha, beta, gamma, vars, dvars, ddvars, t, dd, d2Tdotd, d2Tdotu, d2d,
      d2du, res, mat);

  if (alpha != 0.0) {
    for (int i = 0; i < dsize * dsize; i++) {
      d2d[i] *= 1.0 / alpha;
      d2du[i] *= 1.0 / alpha;
    }
  }

  for (int k = 0; k < size; k++) {
    TacsScalar varst[size], dvarst[size], ddvarst[size];
    memcpy(varst, vars, size * sizeof(TacsScalar));
    memcpy(dvarst, dvars, size * sizeof(TacsScalar));
    memcpy(ddvarst, ddvars, size * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    varst[k] = varst[k] + alpha * TacsScalar(0.0, dh);
    dvarst[k] = dvarst[k] + beta * TacsScalar(0.0, dh);
    ddvarst[k] = ddvarst[k] + gamma * TacsScalar(0.0, dh);
#else
    varst[k] = varst[k] + alpha * dh;
    dvarst[k] = dvarst[k] + beta * dh;
    ddvarst[k] = ddvarst[k] + gamma * dh;
#endif  // TACS_USE_COMPLEX

    // Compute the director rates at the perturbed value
    TacsScalar dt[dsize], ddott[dsize], dddott[dsize];
    director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
        varst, dvarst, ddvarst, t, dt, ddott, dddott);

    TacsScalar rest[size];
    memset(rest, 0, size * sizeof(TacsScalar));

    // Add the change in coefficient
    TacsScalar dTdott[dsize], ddt[dsize];
    for (int j = 0; j < dsize; j++) {
      dTdott[j] = dTdot[j];
      ddt[j] = dd[j];
      for (int i = 0; i < dsize; i++) {
        dTdott[j] += d2Tdotd[j * dsize + i] * (dddott[i] - dddot[i]);
        ddt[j] += d2d[j * dsize + i] * (dt[i] - d[i]);
      }
    }

    for (int j = 0; j < dsize; j++) {
      for (int ii = 0; ii < size; ii++) {
        if (ii % vars_per_node < 3) {
          int i = ii % vars_per_node + 3 * (ii / vars_per_node);
          dTdott[j] += d2Tdotu[j * dsize + i] * (ddvarst[ii] - ddvars[ii]);
          rest[ii] += d2Tdotu[j * dsize + i] * (dddott[j] - dddot[j]);

          ddt[j] += d2du[j * dsize + i] * (varst[ii] - vars[ii]);
          rest[ii] += d2du[j * dsize + i] * (dt[j] - d[j]);
        }
      }
    }

    director::template addDirectorResidual<vars_per_node, offset, num_nodes>(
        varst, dvarst, ddvarst, t, dTdott, ddt, rest);

#ifdef TACS_USE_COMPLEX
    for (int j = 0; j < size; j++) {
      fdmat[size * j + k] = TacsImagPart(rest[j]) / dh;
    }
#else
    for (int j = 0; j < size; j++) {
      fdmat[size * j + k] = (rest[j] - res[j]) / dh;
    }
#endif  // TACS_USE_COMPLEX
  }

  // Compute the error
  max_err = TacsGetMaxError(mat, fdmat, size * size, &max_err_index);
  max_rel = TacsGetMaxRelError(mat, fdmat, size * size, &max_err_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the Jacobian of the director\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "mat", mat, fdmat, size * size);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  return fail;
}

template <int dsize>
void TacsTestEvalDirectorEnergy(const TacsScalar Tlin[],
                                const TacsScalar Tquad[],
                                const TacsScalar Plin[],
                                const TacsScalar Pquad[], const TacsScalar d[],
                                const TacsScalar ddot[], TacsScalar *_T,
                                TacsScalar *_P) {
  TacsScalar T = 0.0;
  TacsScalar P = 0.0;

  for (int j = 0; j < dsize; j++) {
    T += Tlin[j] * (ddot[j] * ddot[j] + ddot[j] * d[j]);
    P += Plin[j] * d[j];

    for (int i = 0; i < dsize; i++) {
      T += Tquad[i + j * dsize] * ddot[i] * ddot[j];
      P += Pquad[i + j * dsize] * d[i] * d[j];
    }
  }

  *_T = T;
  *_P = P;
}

template <int dsize>
void TacsTestEvalDirectorEnergyDerivatives(
    const TacsScalar Tlin[], const TacsScalar Tquad[], const TacsScalar Plin[],
    const TacsScalar Pquad[], const TacsScalar d[], const TacsScalar ddot[],
    const TacsScalar dddot[], TacsScalar dTdot[], TacsScalar dd[]) {
  for (int j = 0; j < dsize; j++) {
    dTdot[j] = Tlin[j] * (2.0 * dddot[j] + ddot[j]);
    dd[j] = Plin[j] - Tlin[j] * ddot[j];
  }

  for (int j = 0; j < dsize; j++) {
    for (int i = 0; i < dsize; i++) {
      dTdot[j] += Tquad[i + j * dsize] * dddot[i];
      dTdot[i] += Tquad[i + j * dsize] * dddot[j];

      dd[j] += Pquad[i + j * dsize] * d[i];
      dd[i] += Pquad[i + j * dsize] * d[j];
    }
  }
}

template <int vars_per_node, int offset, int num_nodes, class director>
int TacsTestDirectorResidual(double dh = 1e-5, int test_print_level = 2,
                             double test_fail_atol = 1e-5,
                             double test_fail_rtol = 1e-5) {
  const int size = vars_per_node * num_nodes;
  const int dsize = 3 * num_nodes;

  // Generate random arrays for the state variables and their time derivatives
  TacsScalar vars[size], dvars[size], ddvars[size];
  TacsGenerateRandomArray(vars, size);
  TacsGenerateRandomArray(dvars, size);
  TacsGenerateRandomArray(ddvars, size);

  // Compute/normalize the normals
  TacsScalar t[dsize];
  TacsGenerateRandomArray(t, dsize);
  for (int i = 0; i < num_nodes; i++) {
    TacsScalar tnrm = sqrt(vec3Dot(&t[3 * i], &t[3 * i]));
    vec3Scale(1.0 / tnrm, &t[3 * i]);
  }

  // The kinetic energy is computed as
  TacsScalar Tlin[dsize], Plin[dsize];
  TacsGenerateRandomArray(Tlin, dsize);
  TacsGenerateRandomArray(Plin, dsize);

  TacsScalar Tquad[dsize * dsize], Pquad[dsize * dsize];
  TacsGenerateRandomArray(Tquad, dsize * dsize);
  TacsGenerateRandomArray(Pquad, dsize * dsize);

  // Compute the director rates
  TacsScalar d[dsize], ddot[dsize], dddot[dsize];
  director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
      vars, dvars, ddvars, t, d, ddot, dddot);

  // Compute the derivatives of the kinetic and potential energies
  TacsScalar dTdot[dsize], dd[dsize];
  TacsTestEvalDirectorEnergyDerivatives<dsize>(Tlin, Tquad, Plin, Pquad, d,
                                               ddot, dddot, dTdot, dd);

  // Compute the residual
  TacsScalar res[size];
  memset(res, 0, size * sizeof(TacsScalar));
  director::template addDirectorResidual<vars_per_node, offset, num_nodes>(
      vars, dvars, ddvars, t, dTdot, dd, res);

  // Compute the values of the variables at (t + dt)
  TacsScalar q[size], qdot[size];
  for (int i = 0; i < size; i++) {
    q[i] = vars[i] + dh * dvars[i] + 0.5 * dh * dh * ddvars[i];
    qdot[i] = dvars[i] + dh * ddvars[i];
  }

  // Evaluate the derivative w.r.t. dot{q}
  TacsScalar res1[size];
  for (int i = 0; i < size; i++) {
    // Evaluate the finite-difference for component i
    TacsScalar dqtmp = qdot[i];
#ifdef TACS_USE_COMPLEX
    TacsScalar T1, P1;
    qdot[i] = dqtmp + TacsScalar(0.0, dh);
    director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
        q, qdot, t, d, ddot);
    TacsTestEvalDirectorEnergy<dsize>(Tlin, Tquad, Plin, Pquad, d, ddot, &T1,
                                      &P1);
    res1[i] = TacsImagPart((T1 - P1)) / dh;
#else
    TacsScalar T1, P1, T2, P2;
    qdot[i] = dqtmp + dh;
    director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
        q, qdot, t, d, ddot);
    TacsTestEvalDirectorEnergy<dsize>(Tlin, Tquad, Plin, Pquad, d, ddot, &T1,
                                      &P1);

    qdot[i] = dqtmp - dh;
    director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
        q, qdot, t, d, ddot);
    TacsTestEvalDirectorEnergy<dsize>(Tlin, Tquad, Plin, Pquad, d, ddot, &T2,
                                      &P2);
    res1[i] = 0.5 * ((T1 - P1) - (T2 - P2)) / dh;
#endif  // TACS_USE_COMPLEX
    qdot[i] = dqtmp;
  }

  // Compute the values of the variables at (t - dt)
  for (int i = 0; i < size; i++) {
    q[i] = vars[i] - dh * dvars[i] + 0.5 * dh * dh * ddvars[i];
    qdot[i] = dvars[i] - dh * ddvars[i];
  }

  // Evaluate the derivative w.r.t. dot{q}
  TacsScalar res2[size];
  for (int i = 0; i < size; i++) {
    // Evaluate the finite-difference for component i
    TacsScalar dqtmp = qdot[i];
#ifdef TACS_USE_COMPLEX
    TacsScalar T1, P1;
    qdot[i] = dqtmp + TacsScalar(0.0, dh);
    director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
        q, qdot, t, d, ddot);
    TacsTestEvalDirectorEnergy<dsize>(Tlin, Tquad, Plin, Pquad, d, ddot, &T1,
                                      &P1);
    res2[i] = TacsImagPart((T1 - P1)) / dh;
#else
    TacsScalar T1, P1, T2, P2;
    qdot[i] = dqtmp + dh;
    director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
        q, qdot, t, d, ddot);
    TacsTestEvalDirectorEnergy<dsize>(Tlin, Tquad, Plin, Pquad, d, ddot, &T1,
                                      &P1);

    qdot[i] = dqtmp - dh;
    director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
        q, qdot, t, d, ddot);
    TacsTestEvalDirectorEnergy<dsize>(Tlin, Tquad, Plin, Pquad, d, ddot, &T2,
                                      &P2);
    res2[i] = 0.5 * ((T1 - P1) - (T2 - P2)) / dh;
#endif  // TACS_USE_COMPLEX
    qdot[i] = dqtmp;
  }

  // Evaluate the finite-difference for the first term in Largrange's
  // equations of motion
  TacsScalar fd[size];
  for (int i = 0; i < size; i++) {
    fd[i] = 0.5 * (res1[i] - res2[i]) / dh;
  }

  // Reset the values of q and dq at time t
  for (int i = 0; i < size; i++) {
    q[i] = vars[i];
    qdot[i] = dvars[i];
  }

  // Compute the contribution from dL/dq^{T}
  for (int i = 0; i < size; i++) {
    // Evaluate the finite-difference for component i
    TacsScalar qtmp = q[i];

#ifdef TACS_USE_COMPLEX
    TacsScalar T1, P1;
    q[i] = qtmp + TacsScalar(0.0, dh);
    director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
        q, qdot, t, d, ddot);
    TacsTestEvalDirectorEnergy<dsize>(Tlin, Tquad, Plin, Pquad, d, ddot, &T1,
                                      &P1);
    res1[i] = TacsImagPart((T1 - P1)) / dh;
#else
    TacsScalar T1, P1, T2, P2;
    q[i] = qtmp + dh;
    director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
        q, qdot, t, d, ddot);
    TacsTestEvalDirectorEnergy<dsize>(Tlin, Tquad, Plin, Pquad, d, ddot, &T1,
                                      &P1);

    q[i] = qtmp - dh;
    director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
        q, qdot, t, d, ddot);
    TacsTestEvalDirectorEnergy<dsize>(Tlin, Tquad, Plin, Pquad, d, ddot, &T2,
                                      &P2);

    // Compute and store the approximation
    res1[i] = 0.5 * ((T1 - P1) - (T2 - P2)) / dh;
#endif  // TACS_USE_COMPLEX
    q[i] = qtmp;
  }

  // Add the result to the finite-difference result
  for (int i = 0; i < size; i++) {
    fd[i] -= res1[i];
  }

  // Variables to store the max error and indices
  int max_err_index, max_rel_index;
  double max_err, max_rel;

  // Keep track of the failure flag
  int fail = 0;

  // Compute the error
  max_err = TacsGetMaxError(res, fd, size, &max_err_index);
  max_rel = TacsGetMaxRelError(res, fd, size, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the residual implementation\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "res", res, fd, size);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  return fail;
}

#endif  // TACS_DIRECTOR_H
