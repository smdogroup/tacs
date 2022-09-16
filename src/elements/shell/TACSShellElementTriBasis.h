#ifndef TACS_SHELL_ELEMENT_TRI_BASIS_H
#define TACS_SHELL_ELEMENT_TRI_BASIS_H

#include "TACSElementAlgebra.h"
#include "TACSElementTypes.h"
#include "TACSLagrangeInterpolation.h"
#include "TACSShellElementQuadrature.h"

class TACSShellTriLinearBasis {
 public:
  static const int NUM_NODES = 3;

  // Set the number of tying points for each of the 5 components
  // of the tying strain
  static const int NUM_TYING_POINTS = 4;

  static void getNodePoint(const int n, double pt[]) {
    if (n == 0) {
      pt[0] = 0.0;
      pt[1] = 0.0;
    } else if (n == 1) {
      pt[0] = 1.0;
      pt[1] = 0.0;
    } else if (n == 2) {
      pt[0] = 0.0;
      pt[1] = 1.0;
    }
  }
  static ElementLayout getLayoutType() { return TACS_TRI_ELEMENT; }

  static inline void computeBasis(const double pt[], double N[]) {
    N[0] = (1.0 - pt[0] - pt[1]);
    N[1] = pt[0];
    N[2] = pt[1];
  }

  static inline void computeBasisGradient(const double pt[], double Nxi[]) {
    Nxi[0] = -1.0;
    Nxi[1] = -1.0;
    Nxi[2] = 1.0;
    Nxi[3] = 0.0;
    Nxi[4] = 0.0;
    Nxi[5] = 1.0;
  }

  template <int vars_per_node, int m>
  static void interpFields(const double pt[], const TacsScalar values[],
                           TacsScalar field[]) {
    double N[NUM_NODES];
    computeBasis(pt, N);

    for (int k = 0; k < m; k++) {
      field[k] = 0.0;
    }

    for (int i = 0; i < NUM_NODES; i++) {
      for (int k = 0; k < m; k++) {
        field[k] += N[i] * values[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void addInterpFieldsTranspose(const double pt[],
                                       const TacsScalar field[],
                                       TacsScalar values[]) {
    double N[NUM_NODES];
    computeBasis(pt, N);

    for (int i = 0; i < NUM_NODES; i++) {
      for (int k = 0; k < m; k++) {
        values[k] += N[i] * field[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void interpFieldsGrad(const double pt[], const TacsScalar values[],
                               TacsScalar grad[]) {
    double Nxi[2 * NUM_NODES];
    computeBasisGradient(pt, Nxi);

    for (int k = 0; k < m; k++) {
      grad[2 * k] = 0.0;
      grad[2 * k + 1] = 0.0;
    }

    for (int i = 0; i < NUM_NODES; i++) {
      for (int k = 0; k < m; k++) {
        grad[2 * k] += Nxi[2 * i] * values[k];
        grad[2 * k + 1] += Nxi[2 * i + 1] * values[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void addInterpFieldsGradTranspose(const double pt[], TacsScalar grad[],
                                           TacsScalar values[]) {
    double Nxi[2 * NUM_NODES];
    computeBasisGradient(pt, Nxi);

    for (int i = 0; i < NUM_NODES; i++) {
      for (int k = 0; k < m; k++) {
        values[k] +=
            (Nxi[2 * i] * grad[2 * k] + Nxi[2 * i + 1] * grad[2 * k + 1]);
      }
      values += vars_per_node;
    }
  }

  /**
    Add the outer-product of the shape functions to the matrix

    jac[m*ix + jy] stores the derivative of the ix term with
    respect to the jx component.

    @param pt The parametric location of the quadrature point
    @param m The number of field components
    @param jac The 2m x 2m Jacobian matrix of coefficients
    @param vars_per_node The number of variables per node
    @param mat The element matrix
  */
  template <int nbrows, int nbcols, int njrows, int njcols>
  static void addInterpFieldsOuterProduct(const double pt[],
                                          const TacsScalar jac[],
                                          TacsScalar *mat) {
    double N[NUM_NODES];
    computeBasis(pt, N);

    const int ncols = NUM_NODES * nbcols;

    for (int jx = 0; jx < NUM_NODES; jx++) {
      const TacsScalar *jac1 = jac;
      for (int jm = 0; jm < njrows; jm++, jac1 += njcols) {
        for (int ix = 0; ix < NUM_NODES; ix++) {
          double Ni = N[jx] * N[ix];

          for (int im = 0; im < njcols; im++) {
            mat[im] += Ni * jac1[im];
          }

          mat += nbcols;
        }
      }

      mat += (nbrows - njrows) * ncols;
    }
  }

  /**
    Add the outer-product of the shape functions to the matrix

    jac[2*m*(2*ix + jx) + 2*iy + jy] stores the derivative of the
    2*ix + jx term with respect to the 2*ix + jx component.

    @param pt The parametric location of the quadrature point
    @param jac The 2*njrows x 2*njcols Jacobian matrix of coefficients
    @param mat The element matrix
  */
  template <int nbrows, int nbcols, int njrows, int njcols>
  static void addInterpGradOuterProduct(const double pt[],
                                        const TacsScalar jac[],
                                        TacsScalar *mat) {
    double Nxi[2 * NUM_NODES];
    computeBasisGradient(pt, Nxi);

    const int ncols = NUM_NODES * nbcols;

    for (int jx = 0; jx < NUM_NODES; jx++) {
      const TacsScalar *jac1 = jac;
      const TacsScalar *jac2 = &jac[2 * njcols];

      for (int jm = 0; jm < njrows;
           jm++, jac1 += 4 * njcols, jac2 += 4 * njcols) {
        for (int ix = 0; ix < NUM_NODES; ix++) {
          double Naa = Nxi[2 * jx] * Nxi[2 * ix];
          double Nab = Nxi[2 * jx] * Nxi[2 * ix + 1];
          double Nba = Nxi[2 * jx + 1] * Nxi[2 * ix];
          double Nbb = Nxi[2 * jx + 1] * Nxi[2 * ix + 1];

          for (int im = 0; im < njcols; im++) {
            mat[im] += (Naa * jac1[2 * im] + Nab * jac1[2 * im + 1] +
                        Nba * jac2[2 * im] + Nbb * jac2[2 * im + 1]);
          }

          mat += nbcols;
        }
      }

      mat += (nbrows - njrows) * ncols;
    }
  }

  /**
    Add the outer-product of the shape functions and their
    derivatives

    Here jac and jacT are the Jacobians of the coefficients of the
    Jacobian matrix of the displacement with respect to the gradient
    of the displacement.

    @param pt The parametric location of the quadrature point
    @param m The number of field components
    @param jac The m x 2m Jacobian matrix of coefficients
    @param jacT The m x 2m Jacobian matrix of coefficients
    @param vars_per_node The number of variables per node
    @param mat The element matrix
  */
  template <int nbrows, int nbcols, int njrows, int njcols>
  static void addInterpGradMixedOuterProduct(const double pt[],
                                             const TacsScalar jac[],
                                             const TacsScalar jacT[],
                                             TacsScalar *mat) {
    double N[NUM_NODES];
    computeBasis(pt, N);

    double Nxi[2 * NUM_NODES];
    computeBasisGradient(pt, Nxi);

    const int ncols = NUM_NODES * nbcols;

    if (jac && jacT) {
      for (int jx = 0; jx < NUM_NODES; jx++) {
        const TacsScalar *jac1 = jac;

        for (int jm = 0; jm < njrows; jm++, jac1 += 2 * njcols) {
          for (int ix = 0; ix < NUM_NODES; ix++) {
            double Na1 = N[jx] * Nxi[2 * ix];
            double Nb1 = N[jx] * Nxi[2 * ix + 1];
            double Na2 = N[ix] * Nxi[2 * jx];
            double Nb2 = N[ix] * Nxi[2 * jx + 1];

            const TacsScalar *jac2 = &jacT[2 * jm];
            for (int im = 0; im < njcols; im++, jac2 += 2 * njcols) {
              mat[im] += Na1 * jac1[2 * im] + Nb1 * jac1[2 * im + 1] +
                         Na2 * jac2[0] + Nb2 * jac2[1];
            }

            mat += nbcols;
          }

          mat += (nbrows - njrows) * ncols;
        }
      }
    } else if (jac) {
      for (int jx = 0; jx < NUM_NODES; jx++) {
        const TacsScalar *jac1 = jac;

        for (int jm = 0; jm < njrows; jm++, jac1 += 2 * njcols) {
          for (int ix = 0; ix < NUM_NODES; ix++) {
            double Na1 = N[jx] * Nxi[2 * ix];
            double Nb1 = N[jx] * Nxi[2 * ix + 1];

            for (int im = 0; im < njcols; im++) {
              mat[im] += Na1 * jac1[2 * im] + Nb1 * jac1[2 * im + 1];
            }

            mat += nbcols;
          }

          mat += (nbrows - njrows) * ncols;
        }
      }
    } else if (jacT) {
      for (int jx = 0; jx < NUM_NODES; jx++) {
        for (int jm = 0; jm < njrows; jm++) {
          for (int ix = 0; ix < NUM_NODES; ix++) {
            double Na2 = N[ix] * Nxi[2 * jx];
            double Nb2 = N[ix] * Nxi[2 * jx + 1];

            const TacsScalar *jac2 = &jacT[2 * jm];
            for (int im = 0; im < njcols; im++, jac2 += 2 * njcols) {
              mat[im] += Na2 * jac2[0] + Nb2 * jac2[1];
            }

            mat += nbcols;
          }

          mat += (nbrows - njrows) * ncols;
        }
      }
    }
  }

  /*
    Get the tying strain component
  */
  static inline TacsShellTyingStrainComponent getTyingField(int ty_index) {
    if (ty_index == 0 || ty_index == 1) {
      // g13
      return TACS_SHELL_G13_COMPONENT;
    } else {  // if (ty_index == 2 || ty_index == 3){
      // g23
      return TACS_SHELL_G23_COMPONENT;
    }
  }

  /*
    Get the parametric point within the element
  */
  static inline void getTyingPoint(int ty_index, double pt[]) {
    if (ty_index == 0) {  // g13
      pt[0] = 0.5;
      pt[1] = 0.0;
    } else if (ty_index == 1) {  // g13
      pt[0] = 0.5;
      pt[1] = 0.5;
    } else if (ty_index == 2) {  // g23
      pt[0] = 0.0;
      pt[1] = 0.5;
    } else if (ty_index == 3) {  // g23
      pt[0] = 0.5;
      pt[1] = 0.5;
    }
  }

  /**
    Evaluate the tensorial components of the strain tensor at the
    specific quadrature point

    gty = [g11  g12  g13]
    .     [sym  g22  g23]
    .     [sym  sym  g33]

    As a result: gty[0] = g11, gty[1] = g12, gty[2] = g13, gty[3] = g22
    and gty[4] = g23, with gty[5] = 0.0

    @param pt The quadrature point
    @param ety The strain computed at the tying points
    @param gty The interpolated tying strain
  */
  static inline void interpTyingStrain(const double pt[],
                                       const TacsScalar ety[],
                                       TacsScalar gty[]) {
    // Perform the interpolation
    // cc = g23[1] - g13[0] - g23[2] + g13[2];
    TacsScalar cc = ety[2] - ety[0] - ety[3] + ety[1];

    // Interpolate g13 == g13[0] + cc * pt[1]
    gty[2] = ety[0] + pt[1] * cc;

    // Interpolate g23 = g23[1] - cc * pt[0]
    gty[4] = ety[2] - pt[0] * cc;

    // Set all remaining entries to zero
    gty[0] = gty[1] = gty[3] = gty[5] = 0.0;
  }

  /**
    Add the derivative of the tying strain to the residual

    @param pt The quadrature point
    @param dgty The derivative of the interpolated strain
    @param dety The output derivative of the strain at the tying points
  */
  static inline void addInterpTyingStrainTranspose(const double pt[],
                                                   const TacsScalar dgty[],
                                                   TacsScalar dety[]) {
    dety[0] += dgty[2];
    dety[2] += dgty[4];

    TacsScalar dcc = pt[1] * dgty[2] - pt[0] * dgty[4];

    dety[2] += dcc;
    dety[0] -= dcc;
    dety[3] -= dcc;
    dety[1] += dcc;
  }

  /**
    Add the second derivative of the tying strain at the tying points

    @param pt The quadrature point
    @param d2gty The second derivative of the interpolated strain
    @param d2ety The second derivatives of the strain at the tying points
  */
  static inline void addInterpTyingStrainHessian(const double pt[],
                                                 const TacsScalar d2gty[],
                                                 TacsScalar d2ety[]) {
    TacsScalar e0 = (1.0 - pt[1]) * d2gty[2 * 6 + 2] + pt[0] * d2gty[2 * 6 + 4];
    TacsScalar f0 = (1.0 - pt[1]) * d2gty[4 * 6 + 2] + pt[0] * d2gty[4 * 6 + 4];

    TacsScalar e1 = pt[1] * d2gty[2 * 6 + 2] - pt[0] * d2gty[2 * 6 + 4];
    TacsScalar f1 = pt[1] * d2gty[4 * 6 + 2] - pt[0] * d2gty[4 * 6 + 4];

    TacsScalar e2 = pt[1] * d2gty[2 * 6 + 2] + (1.0 - pt[0]) * d2gty[2 * 6 + 4];
    TacsScalar f2 = pt[1] * d2gty[4 * 6 + 2] + (1.0 - pt[0]) * d2gty[4 * 6 + 4];

    d2ety[0] += (1.0 - pt[1]) * e0 + pt[0] * f0;
    d2ety[1] += (1.0 - pt[1]) * e1 + pt[0] * f1;
    d2ety[2] += (1.0 - pt[1]) * e2 + pt[0] * f2;
    d2ety[3] -= (1.0 - pt[1]) * e1 + pt[0] * f1;

    d2ety[4] += pt[1] * e0 - pt[0] * f0;
    d2ety[5] += pt[1] * e1 - pt[0] * f1;
    d2ety[6] += pt[1] * e2 - pt[0] * f2;
    d2ety[7] -= pt[1] * e1 - pt[0] * f1;

    d2ety[8] += pt[1] * e0 + (1.0 - pt[0]) * f0;
    d2ety[9] += pt[1] * e1 + (1.0 - pt[0]) * f1;
    d2ety[10] += pt[1] * e2 + (1.0 - pt[0]) * f2;
    d2ety[11] -= pt[1] * e1 + (1.0 - pt[0]) * f1;

    d2ety[12] -= pt[1] * e0 - pt[0] * f0;
    d2ety[13] -= pt[1] * e1 - pt[0] * f1;
    d2ety[14] -= pt[1] * e2 - pt[0] * f2;
    d2ety[15] += pt[1] * e1 - pt[0] * f1;
  }
};

class TACSShellTriQuadraticBasis {
 public:
  static const int NUM_NODES = 6;

  static void getNodePoint(const int n, double pt[]) {
    if (n == 0) {
      pt[0] = 0.0;
      pt[1] = 0.0;
    } else if (n == 1) {
      pt[0] = 1.0;
      pt[1] = 0.0;
    } else if (n == 2) {
      pt[0] = 0.0;
      pt[1] = 1.0;
    } else if (n == 3) {
      pt[0] = 0.5;
      pt[1] = 0.0;
    } else if (n == 4) {
      pt[0] = 0.5;
      pt[1] = 0.5;
    } else if (n == 5) {
      pt[0] = 0.0;
      pt[1] = 0.5;
    }
  }
  static ElementLayout getLayoutType() { return TACS_TRI_QUADRATIC_ELEMENT; }

  static inline void computeBasis(const double pt[], double N[]) {
    N[0] = (1.0 - pt[0] - pt[1]) * (1.0 - 2.0 * pt[0] - 2.0 * pt[1]);
    N[1] = pt[0] * (2.0 * pt[0] - 1.0);
    N[2] = pt[1] * (2.0 * pt[1] - 1.0);
    N[3] = 4.0 * pt[0] * (1.0 - pt[0] - pt[1]);
    N[4] = 4.0 * pt[0] * pt[1];
    N[5] = 4.0 * pt[1] * (1.0 - pt[0] - pt[1]);
  }

  static inline void computeBasisGradient(const double pt[], double Nxi[]) {
    Nxi[0] = 4.0 * pt[0] + 4.0 * pt[1] - 3.0;
    Nxi[1] = 4.0 * pt[0] + 4.0 * pt[1] - 3.0;
    Nxi[2] = 4.0 * pt[0] - 1.0;
    Nxi[3] = 0.0;
    Nxi[4] = 0.0;
    Nxi[5] = 4.0 * pt[1] - 1.0;
    Nxi[6] = 4.0 - 8.0 * pt[0] - 4.0 * pt[1];
    Nxi[7] = -4.0 * pt[0];
    Nxi[8] = 4.0 * pt[1];
    Nxi[9] = 4.0 * pt[0];
    Nxi[10] = -4.0 * pt[1];
    Nxi[11] = 4.0 - 4.0 * pt[0] - 8.0 * pt[1];
  }

  template <int vars_per_node, int m>
  static void interpFields(const double pt[], const TacsScalar values[],
                           TacsScalar field[]) {
    double N[NUM_NODES];
    computeBasis(pt, N);

    for (int k = 0; k < m; k++) {
      field[k] = 0.0;
    }

    for (int i = 0; i < NUM_NODES; i++) {
      for (int k = 0; k < m; k++) {
        field[k] += N[i] * values[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void addInterpFieldsTranspose(const double pt[],
                                       const TacsScalar field[],
                                       TacsScalar values[]) {
    double N[NUM_NODES];
    computeBasis(pt, N);

    for (int i = 0; i < NUM_NODES; i++) {
      for (int k = 0; k < m; k++) {
        values[k] += N[i] * field[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void interpFieldsGrad(const double pt[], const TacsScalar values[],
                               TacsScalar grad[]) {
    double Nxi[2 * NUM_NODES];
    computeBasisGradient(pt, Nxi);

    for (int k = 0; k < m; k++) {
      grad[2 * k] = 0.0;
      grad[2 * k + 1] = 0.0;
    }

    for (int i = 0; i < NUM_NODES; i++) {
      for (int k = 0; k < m; k++) {
        grad[2 * k] += Nxi[2 * i] * values[k];
        grad[2 * k + 1] += Nxi[2 * i + 1] * values[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void addInterpFieldsGradTranspose(const double pt[], TacsScalar grad[],
                                           TacsScalar values[]) {
    double Nxi[2 * NUM_NODES];
    computeBasisGradient(pt, Nxi);

    for (int i = 0; i < NUM_NODES; i++) {
      for (int k = 0; k < m; k++) {
        values[k] +=
            (Nxi[2 * i] * grad[2 * k] + Nxi[2 * i + 1] * grad[2 * k + 1]);
      }
      values += vars_per_node;
    }
  }

  /**
    Add the outer-product of the shape functions to the matrix

    jac[m*ix + jy] stores the derivative of the ix term with
    respect to the jx component.

    @param pt The parametric location of the quadrature point
    @param m The number of field components
    @param jac The 2m x 2m Jacobian matrix of coefficients
    @param vars_per_node The number of variables per node
    @param mat The element matrix
  */
  template <int nbrows, int nbcols, int njrows, int njcols>
  static void addInterpFieldsOuterProduct(const double pt[],
                                          const TacsScalar jac[],
                                          TacsScalar *mat) {
    double N[NUM_NODES];
    computeBasis(pt, N);

    const int ncols = NUM_NODES * nbcols;

    for (int jx = 0; jx < NUM_NODES; jx++) {
      const TacsScalar *jac1 = jac;
      for (int jm = 0; jm < njrows; jm++, jac1 += njcols) {
        for (int ix = 0; ix < NUM_NODES; ix++) {
          double Ni = N[jx] * N[ix];

          for (int im = 0; im < njcols; im++) {
            mat[im] += Ni * jac1[im];
          }

          mat += nbcols;
        }
      }

      mat += (nbrows - njrows) * ncols;
    }
  }

  /**
    Add the outer-product of the shape functions to the matrix

    jac[2*m*(2*ix + jx) + 2*iy + jy] stores the derivative of the
    2*ix + jx term with respect to the 2*ix + jx component.

    @param pt The parametric location of the quadrature point
    @param jac The 2*njrows x 2*njcols Jacobian matrix of coefficients
    @param mat The element matrix
  */
  template <int nbrows, int nbcols, int njrows, int njcols>
  static void addInterpGradOuterProduct(const double pt[],
                                        const TacsScalar jac[],
                                        TacsScalar *mat) {
    double Nxi[2 * NUM_NODES];
    computeBasisGradient(pt, Nxi);

    const int ncols = NUM_NODES * nbcols;

    for (int jx = 0; jx < NUM_NODES; jx++) {
      const TacsScalar *jac1 = jac;
      const TacsScalar *jac2 = &jac[2 * njcols];

      for (int jm = 0; jm < njrows;
           jm++, jac1 += 4 * njcols, jac2 += 4 * njcols) {
        for (int ix = 0; ix < NUM_NODES; ix++) {
          double Naa = Nxi[2 * jx] * Nxi[2 * ix];
          double Nab = Nxi[2 * jx] * Nxi[2 * ix + 1];
          double Nba = Nxi[2 * jx + 1] * Nxi[2 * ix];
          double Nbb = Nxi[2 * jx + 1] * Nxi[2 * ix + 1];

          for (int im = 0; im < njcols; im++) {
            mat[im] += (Naa * jac1[2 * im] + Nab * jac1[2 * im + 1] +
                        Nba * jac2[2 * im] + Nbb * jac2[2 * im + 1]);
          }

          mat += nbcols;
        }
      }

      mat += (nbrows - njrows) * ncols;
    }
  }

  /**
    Add the outer-product of the shape functions and their
    derivatives

    Here jac and jacT are the Jacobians of the coefficients of the
    Jacobian matrix of the displacement with respect to the gradient
    of the displacement.

    @param pt The parametric location of the quadrature point
    @param m The number of field components
    @param jac The m x 2m Jacobian matrix of coefficients
    @param jacT The m x 2m Jacobian matrix of coefficients
    @param vars_per_node The number of variables per node
    @param mat The element matrix
  */
  template <int nbrows, int nbcols, int njrows, int njcols>
  static void addInterpGradMixedOuterProduct(const double pt[],
                                             const TacsScalar jac[],
                                             const TacsScalar jacT[],
                                             TacsScalar *mat) {
    double N[NUM_NODES];
    computeBasis(pt, N);

    double Nxi[2 * NUM_NODES];
    computeBasisGradient(pt, Nxi);

    const int ncols = NUM_NODES * nbcols;

    if (jac && jacT) {
      for (int jx = 0; jx < NUM_NODES; jx++) {
        const TacsScalar *jac1 = jac;

        for (int jm = 0; jm < njrows; jm++, jac1 += 2 * njcols) {
          for (int ix = 0; ix < NUM_NODES; ix++) {
            double Na1 = N[jx] * Nxi[2 * ix];
            double Nb1 = N[jx] * Nxi[2 * ix + 1];
            double Na2 = N[ix] * Nxi[2 * jx];
            double Nb2 = N[ix] * Nxi[2 * jx + 1];

            const TacsScalar *jac2 = &jacT[2 * jm];
            for (int im = 0; im < njcols; im++, jac2 += 2 * njcols) {
              mat[im] += Na1 * jac1[2 * im] + Nb1 * jac1[2 * im + 1] +
                         Na2 * jac2[0] + Nb2 * jac2[1];
            }

            mat += nbcols;
          }

          mat += (nbrows - njrows) * ncols;
        }
      }
    } else if (jac) {
      for (int jx = 0; jx < NUM_NODES; jx++) {
        const TacsScalar *jac1 = jac;

        for (int jm = 0; jm < njrows; jm++, jac1 += 2 * njcols) {
          for (int ix = 0; ix < NUM_NODES; ix++) {
            double Na1 = N[jx] * Nxi[2 * ix];
            double Nb1 = N[jx] * Nxi[2 * ix + 1];

            for (int im = 0; im < njcols; im++) {
              mat[im] += Na1 * jac1[2 * im] + Nb1 * jac1[2 * im + 1];
            }

            mat += nbcols;
          }

          mat += (nbrows - njrows) * ncols;
        }
      }
    } else if (jacT) {
      for (int jx = 0; jx < NUM_NODES; jx++) {
        for (int jm = 0; jm < njrows; jm++) {
          for (int ix = 0; ix < NUM_NODES; ix++) {
            double Na2 = N[ix] * Nxi[2 * jx];
            double Nb2 = N[ix] * Nxi[2 * jx + 1];

            const TacsScalar *jac2 = &jacT[2 * jm];
            for (int im = 0; im < njcols; im++, jac2 += 2 * njcols) {
              mat[im] += Na2 * jac2[0] + Nb2 * jac2[1];
            }

            mat += nbcols;
          }

          mat += (nbrows - njrows) * ncols;
        }
      }
    }
  }

  // static int getNumTyingFields(){
  //   return 5;
  // }
  // static int getNumTyingPoints( const int field ){
  //   if (field == 0){ return NUM_G11_TYING_POINTS; }
  //   else if (field == 1){ return NUM_G22_TYING_POINTS; }
  //   else if (field == 2){ return NUM_G12_TYING_POINTS; }
  //   else if (field == 3){ return NUM_G13_TYING_POINTS; }
  //   else if (field == 4){ return NUM_G23_TYING_POINTS; }
  //   return 0;
  // }
  // static void getTyingPoint( const int field,
  //                            const int ty,
  //                            double pt[] ){
  //   const double s = 0.774596669241483;
  //   const double t = 0.577350269189626;
  //   const double s0 = 0.5 - 0.5*s;
  //   const double t0 = 0.5 - 0.5*t;
  //   const double t1 = 0.5 + 0.5*t;

  //   if (field == 0 || field == 4){ // g11 or g13
  //     if (ty == 0){
  //       pt[0] = t0;
  //       pt[1] = 0.0;
  //     }
  //     else if (ty == 1){
  //       pt[0] = t1;
  //       pt[1] = 0.0;
  //     }
  //     else if (ty == 2){
  //       pt[0] = t0;
  //       pt[1] = t;
  //     }
  //   }
  //   else if (field == 1 || field == 3){ // g22 or g23
  //     if (ty == 0){
  //       pt[0] = 0.0;
  //       pt[1] = t0;
  //     }
  //     else if (ty == 1){
  //       pt[0] = t;
  //       pt[1] = t0;
  //     }
  //     else if (ty == 2){
  //       pt[0] = 0.0;
  //       pt[1] = t1;
  //     }
  //   }
  //   else { // (field == 2) g12
  //     if (ty == 0){
  //       pt[0] = t0;
  //       pt[1] = t0;
  //     }
  //     else if (ty == 1){
  //       pt[0] = t1;
  //       pt[1] = t0;
  //     }
  //     else if (ty == 2){
  //       pt[0] = t0;
  //       pt[1] = t1;
  //     }
  //   }
  // }
  // static TacsScalar interpTying( const int field,
  //                                const double pt[],
  //                                const TacsScalar ety[] ){
  //   const double s = 0.774596669241483;
  //   const double t = 0.577350269189626;
  //   const double s0 = 0.5 - 0.5*s;
  //   const double t0 = 0.5 - 0.5*t;
  //   const double tinv = 1.0/t;

  //   double N[3];
  //   if (field == 0 || field == 4){
  //     N[0] = 1.0 - tinv*((pt[0] - t0) + pt[1]);
  //     N[1] = tinv*(pt[0] - t0);
  //     N[2] = tinv*pt[1];
  //   }
  //   else if (field == 1 || field == 3){
  //     N[0] = 1.0 - tinv*(pt[0] + (pt[1] - t0));
  //     N[1] = tinv*pt[0];
  //     N[2] = tinv*(pt[1] - t0);
  //   }
  //   else { // field == 2
  //     N[0] = 1.0 - tinv*((pt[0] - t0) + (pt[1] - t0));
  //     N[1] = tinv*(pt[0] - t0);
  //     N[2] = tinv*(pt[1] - t0);
  //   }

  //   TacsScalar value = 0.0;
  //   for ( int i = 0; i < 3; i++ ){
  //     value += N[i]*ety[i];
  //   }

  //   return value;
  // }

  // static void addInterpTyingTranspose( const int field,
  //                                      const double pt[],
  //                                      const TacsScalar value,
  //                                      TacsScalar ety[] ){
  //   const double s = 0.774596669241483;
  //   const double t = 0.577350269189626;
  //   const double s0 = 0.5 - 0.5*s;
  //   const double t0 = 0.5 - 0.5*t;
  //   const double tinv = 1.0/t;

  //   double N[3];
  //   if (field == 0 || field == 4){
  //     N[0] = 1.0 - tinv*((pt[0] - t0) + pt[1]);
  //     N[1] = tinv*(pt[0] - t0);
  //     N[2] = tinv*pt[1];
  //   }
  //   else if (field == 1 || field == 3){
  //     N[0] = 1.0 - tinv*(pt[0] + (pt[1] - t0));
  //     N[1] = tinv*pt[0];
  //     N[2] = tinv*(pt[1] - t0);
  //   }
  //   else { // field == 2
  //     N[0] = 1.0 - tinv*((pt[0] - t0) + (pt[1] - t0));
  //     N[1] = tinv*(pt[0] - t0);
  //     N[2] = tinv*(pt[1] - t0);
  //   }

  //   for ( int i = 0; i < 3; i++ ){
  //     ety[i] += N[i]*value;
  //   }
  // }

  // static void addInterpTyingOuterProduct( const int f1,
  //                                         const int f2,
  //                                         const double pt[],
  //                                         const TacsScalar value,
  //                                         TacsScalar d2ety[] ){
  //   const double s = 0.774596669241483;
  //   const double t = 0.577350269189626;
  //   const double s0 = 0.5 - 0.5*s;
  //   const double t0 = 0.5 - 0.5*t;
  //   const double tinv = 1.0/t;

  //   double N1[3];
  //   if (f1 == 0 || f1 == 4){
  //     N1[0] = 1.0 - tinv*((pt[0] - t0) + pt[1]);
  //     N1[1] = tinv*(pt[0] - t0);
  //     N1[2] = tinv*pt[1];
  //   }
  //   else if (f1 == 1 || f1 == 3){
  //     N1[0] = 1.0 - tinv*(pt[0] + (pt[1] - t0));
  //     N1[1] = tinv*pt[0];
  //     N1[2] = tinv*(pt[1] - t0);
  //   }
  //   else { // f1 == 2
  //     N1[0] = 1.0 - tinv*((pt[0] - t0) + (pt[1] - t0));
  //     N1[1] = tinv*(pt[0] - t0);
  //     N1[2] = tinv*(pt[1] - t0);
  //   }

  //   double N2[3];
  //   if (f2 == 0 || f2 == 4){
  //     N2[0] = 1.0 - tinv*((pt[0] - t0) + pt[1]);
  //     N2[1] = tinv*(pt[0] - t0);
  //     N2[2] = tinv*pt[1];
  //   }
  //   else if (f2 == 1 || f2 == 3){
  //     N2[0] = 1.0 - tinv*(pt[0] + (pt[1] - t0));
  //     N2[1] = tinv*pt[0];
  //     N2[2] = tinv*(pt[1] - t0);
  //   }
  //   else { // f2 == 2
  //     N2[0] = 1.0 - tinv*((pt[0] - t0) + (pt[1] - t0));
  //     N2[1] = tinv*(pt[0] - t0);
  //     N2[2] = tinv*(pt[1] - t0);
  //   }

  //   for ( int i = 0; i < 3; i++ ){
  //     for ( int j = 0; j < 3; j++, d2ety++ ){
  //       d2ety[0] += value*N1[i]*N2[j];
  //     }
  //   }
  // }
};

#endif  // TACS_SHELL_ELEMENT_TRI_BASIS_H
