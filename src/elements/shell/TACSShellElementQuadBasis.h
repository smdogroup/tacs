#ifndef TACS_SHELL_ELEMENT_QUAD_BASIS_H
#define TACS_SHELL_ELEMENT_QUAD_BASIS_H

#include "TACSElementAlgebra.h"
#include "TACSElementTypes.h"
#include "TACSLagrangeInterpolation.h"
#include "TACSShellElementQuadrature.h"

enum TacsShellTyingStrainComponent {
  TACS_SHELL_G11_COMPONENT = 0,
  TACS_SHELL_G22_COMPONENT = 1,
  TACS_SHELL_G12_COMPONENT = 2,
  TACS_SHELL_G23_COMPONENT = 3,
  TACS_SHELL_G13_COMPONENT = 4
};

template <int order>
inline void TacsLagrangeShapeFunction(const double u, const double knots[],
                                      double N[]) {
  // Loop over the shape functions
  for (int i = 0; i < order; i++) {
    N[i] = 1.0;
    for (int j = 0; j < order; j++) {
      if (i != j) {
        double d = 1.0 / (knots[i] - knots[j]);
        N[i] *= (u - knots[j]) * d;
      }
    }
  }
}

template <int order>
inline void TacsLagrangeShapeFuncDerivative(const double u,
                                            const double knots[], double N[],
                                            double Nd[]) {
  // Loop over the shape function knot locations
  for (int i = 0; i < order; i++) {
    N[i] = 1.0;
    Nd[i] = 0.0;

    // Loop over each point again, except for the current control
    // point, adding the contribution to the shape function
    for (int j = 0; j < order; j++) {
      if (i != j) {
        double d = 1.0 / (knots[i] - knots[j]);
        N[i] *= (u - knots[j]) * d;

        // Now add up the contribution to the derivative
        for (int k = 0; k < order; k++) {
          if (k != i && k != j) {
            d *= (u - knots[k]) / (knots[i] - knots[k]);
          }
        }

        // Add the derivative contribution
        Nd[i] += d;
      }
    }
  }
}

template <int order>
inline void TacsLagrangeLobattoShapeFunction(const double u, double *N) {
  if (order == 1) {
    N[0] = 1.0;
  } else if (order == 2) {
    N[0] = 0.5 * (1.0 - u);
    N[1] = 0.5 * (1.0 + u);
  } else if (order == 3) {
    N[0] = -0.5 * u * (1.0 - u);
    N[1] = (1.0 - u) * (1.0 + u);
    N[2] = 0.5 * (1.0 + u) * u;
  } else {
    const double *knots = TacsGaussLobattoPoints4;
    if (order == 5) {
      knots = TacsGaussLobattoPoints5;
    } else if (order == 6) {
      knots = TacsGaussLobattoPoints6;
    }

    TacsLagrangeShapeFunction<order>(u, knots, N);
  }
}

template <int order>
inline void TacsLagrangeLobattoShapeFuncDerivative(const double u, double *N,
                                                   double *Nd) {
  if (order == 1) {
    N[0] = 1.0;
  } else if (order == 2) {
    N[0] = 0.5 * (1.0 - u);
    N[1] = 0.5 * (1.0 + u);

    Nd[0] = -0.5;
    Nd[1] = 0.5;
  } else if (order == 3) {
    N[0] = -0.5 * u * (1.0 - u);
    N[1] = (1.0 - u) * (1.0 + u);
    N[2] = 0.5 * (1.0 + u) * u;

    Nd[0] = -0.5 + u;
    Nd[1] = -2.0 * u;
    Nd[2] = 0.5 + u;
  } else {
    const double *knots = TacsGaussLobattoPoints4;
    if (order == 5) {
      knots = TacsGaussLobattoPoints5;
    } else if (order == 6) {
      knots = TacsGaussLobattoPoints6;
    }

    TacsLagrangeShapeFuncDerivative<order>(u, knots, N, Nd);
  }
}

const double TacsShellLinearTyingPoints[2] = {-1.0, 1.0};

template <int order>
class TACSShellQuadBasis {
 public:
  static const int NUM_NODES = order * order;

  // Set the number of tying points for each of the 5 components
  // of the tying strain
  static const int NUM_G11_TYING_POINTS = order * (order - 1);
  static const int NUM_G22_TYING_POINTS = order * (order - 1);
  static const int NUM_G12_TYING_POINTS = (order - 1) * (order - 1);
  static const int NUM_G13_TYING_POINTS = order * (order - 1);
  static const int NUM_G23_TYING_POINTS = order * (order - 1);

  // Set the offsets used to interpolate the tying strain
  static const int G11_OFFSET = NUM_G11_TYING_POINTS;
  static const int G22_OFFSET = (NUM_G11_TYING_POINTS + NUM_G22_TYING_POINTS);
  static const int G12_OFFSET =
      (NUM_G11_TYING_POINTS + NUM_G22_TYING_POINTS + NUM_G12_TYING_POINTS);
  static const int G13_OFFSET = (NUM_G11_TYING_POINTS + NUM_G22_TYING_POINTS +
                                 NUM_G12_TYING_POINTS + NUM_G13_TYING_POINTS);

  // Set the total number of tying points
  static const int NUM_TYING_POINTS =
      NUM_G11_TYING_POINTS + NUM_G22_TYING_POINTS + NUM_G12_TYING_POINTS +
      NUM_G13_TYING_POINTS + NUM_G23_TYING_POINTS;

  /*
    Get the parametric points of each node in the element
  */
  static void getNodePoint(const int n, double pt[]) {
    pt[0] = -1.0 + (2.0 / (order - 1)) * (n % order);
    pt[1] = -1.0 + (2.0 / (order - 1)) * (n / order);
  }

  /*
    Get the layout type of the element
  */
  static ElementLayout getLayoutType() {
    if (order == 2) {
      return TACS_QUAD_ELEMENT;
    } else if (order == 3) {
      return TACS_QUAD_QUADRATIC_ELEMENT;
    } else if (order == 4) {
      return TACS_QUAD_CUBIC_ELEMENT;
    } else if (order == 5) {
      return TACS_QUAD_QUARTIC_ELEMENT;
    } else if (order == 6) {
      return TACS_QUAD_QUINTIC_ELEMENT;
    }

    return TACS_LAYOUT_NONE;
  }

  template <int vars_per_node, int m>
  static void interpFields(const double pt[], const TacsScalar values[],
                           TacsScalar field[]) {
    double na[order], nb[order];
    TacsLagrangeLobattoShapeFunction<order>(pt[0], na);
    TacsLagrangeLobattoShapeFunction<order>(pt[1], nb);

    for (int k = 0; k < m; k++) {
      field[k] = 0.0;
    }

    for (int j = 0; j < order; j++) {
      for (int i = 0; i < order; i++) {
        for (int k = 0; k < m; k++) {
          field[k] += na[i] * nb[j] * values[k];
        }
        values += vars_per_node;
      }
    }
  }

  template <int vars_per_node, int m>
  static void addInterpFieldsTranspose(const double pt[],
                                       const TacsScalar field[],
                                       TacsScalar values[]) {
    double na[order], nb[order];
    TacsLagrangeLobattoShapeFunction<order>(pt[0], na);
    TacsLagrangeLobattoShapeFunction<order>(pt[1], nb);

    for (int j = 0; j < order; j++) {
      for (int i = 0; i < order; i++) {
        for (int k = 0; k < m; k++) {
          values[k] += na[i] * nb[j] * field[k];
        }
        values += vars_per_node;
      }
    }
  }

  template <int vars_per_node, int m>
  static void interpFieldsGrad(const double pt[], const TacsScalar values[],
                               TacsScalar grad[]) {
    double na[order], dna[order];
    double nb[order], dnb[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[1], nb, dnb);

    for (int k = 0; k < m; k++) {
      grad[2 * k] = 0.0;
      grad[2 * k + 1] = 0.0;
    }

    for (int j = 0; j < order; j++) {
      for (int i = 0; i < order; i++) {
        for (int k = 0; k < m; k++) {
          grad[2 * k] += dna[i] * nb[j] * values[k];
          grad[2 * k + 1] += na[i] * dnb[j] * values[k];
        }
        values += vars_per_node;
      }
    }
  }

  template <int vars_per_node, int m>
  static void addInterpFieldsGradTranspose(const double pt[], TacsScalar grad[],
                                           TacsScalar values[]) {
    double na[order], dna[order];
    double nb[order], dnb[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[1], nb, dnb);

    for (int j = 0; j < order; j++) {
      for (int i = 0; i < order; i++) {
        for (int k = 0; k < m; k++) {
          values[k] +=
              (dna[i] * nb[j] * grad[2 * k] + na[i] * dnb[j] * grad[2 * k + 1]);
        }
        values += vars_per_node;
      }
    }
  }

  /**
    Add the outer-product of the shape functions to the matrix.

    The Jacobian matrix 'mat' is a block matrix with block size nbrows
    x nbcols. The input matrix 'jac' is of size njrows x njcols and
    jac[njcols*i + j] stores the derivative of the i-th term with
    respect to the j-th component.

    The number of column and row blocks in the Jacobian mat matrix is
    equal, however the matrix need not be square because nbrows may
    not be equal to nbcols. However, you must have that nbrows <=
    njrows and nbcols <= njcols.

    @param pt The parametric location of the quadrature point
    @param jac The njrows x njcols Jacobian matrix of coefficients
    @param mat The Jacobian matrix
  */
  template <int nbrows, int nbcols, int njrows, int njcols>
  static void addInterpFieldsOuterProduct(const double pt[],
                                          const TacsScalar jac[],
                                          TacsScalar *mat) {
    double na[order], nb[order];
    TacsLagrangeLobattoShapeFunction<order>(pt[0], na);
    TacsLagrangeLobattoShapeFunction<order>(pt[1], nb);

    const int ncols = NUM_NODES * nbcols;

    for (int jy = 0; jy < order; jy++) {
      for (int jx = 0; jx < order; jx++) {
        double Nj = na[jx] * nb[jy];

        const TacsScalar *jac1 = jac;
        for (int jm = 0; jm < njrows; jm++, jac1 += njcols) {
          for (int iy = 0; iy < order; iy++) {
            for (int ix = 0; ix < order; ix++) {
              double Ni = na[ix] * nb[iy];

              for (int im = 0; im < njcols; im++) {
                mat[im] += Ni * Nj * jac1[im];
              }

              mat += nbcols;
            }
          }
        }

        mat += (nbrows - njrows) * ncols;
      }
    }
  }

  /**
    Add the outer-product of the shape functions to the matrix

    The Jacobian matrix 'mat' is a block matrix with block size nbrows
    x nbcols. The input matrix 'jac' is of size 2*njrows x 2*njcols
    and jac[2*m*(2*ix + jx) + 2*iy + jy] stores the derivative of the
    2*ix + jx term with respect to the 2*ix + jx component.

    The number of column and row blocks in the Jacobian mat matrix is
    equal, however the matrix need not be square because nbrows may
    not be equal to nbcols. However, you must have that nbrows <=
    njrows and nbcols <= njcols.

    @param pt The parametric location of the quadrature point
    @param jac The 2*njrows x 2*njcols Jacobian matrix of coefficients
    @param mat The element matrix
  */
  template <int nbrows, int nbcols, int njrows, int njcols>
  static void addInterpGradOuterProduct(const double pt[],
                                        const TacsScalar jac[],
                                        TacsScalar *mat) {
    double na[order], dna[order];
    double nb[order], dnb[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[1], nb, dnb);

    const int ncols = NUM_NODES * nbcols;

    for (int jy = 0; jy < order; jy++) {
      for (int jx = 0; jx < order; jx++) {
        double Naj = dna[jx] * nb[jy];
        double Nbj = na[jx] * dnb[jy];

        const TacsScalar *jac1 = jac;
        const TacsScalar *jac2 = &jac[2 * njcols];
        for (int jm = 0; jm < njrows;
             jm++, jac1 += 4 * njcols, jac2 += 4 * njcols) {
          for (int iy = 0; iy < order; iy++) {
            for (int ix = 0; ix < order; ix++) {
              double Nai = dna[ix] * nb[iy];
              double Nbi = na[ix] * dnb[iy];
              double Naa = Naj * Nai;
              double Nab = Naj * Nbi;
              double Nba = Nbj * Nai;
              double Nbb = Nbj * Nbi;

              for (int im = 0; im < njcols; im++) {
                mat[im] += (Naa * jac1[2 * im] + Nab * jac1[2 * im + 1] +
                            Nba * jac2[2 * im] + Nbb * jac2[2 * im + 1]);
              }

              mat += nbcols;
            }
          }
        }

        mat += (nbrows - njrows) * ncols;
      }
    }
  }

  /**
    Add the outer product of the shape functions and their derivatives
    to a matrix with a rectangular layout.

    The Jacobian matrix 'mat' is a block matrix with block size nbrows
    x nbcols.

    jac is a njrows x 2*njcols Jacobian matrix
    jacT is a 2*njrows x njcols Jacobian matrix

    @param pt The parametric location of the quadrature point
    @param jac The 2m x 2m Jacobian matrix of coefficients
    @param mat The element matrix
  */
  template <int nbrows, int nbcols, int njrows, int njcols>
  static void addInterpGradMixedOuterProduct(const double pt[],
                                             const TacsScalar jac[],
                                             const TacsScalar jacT[],
                                             TacsScalar *mat) {
    double na[order], dna[order];
    double nb[order], dnb[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[1], nb, dnb);

    const int ncols = NUM_NODES * nbcols;

    if (jac && jacT) {
      for (int jy = 0; jy < order; jy++) {
        for (int jx = 0; jx < order; jx++) {
          double Nj = na[jx] * nb[jy];
          double Naj = dna[jx] * nb[jy];
          double Nbj = na[jx] * dnb[jy];

          const TacsScalar *jac1 = jac;
          for (int jm = 0; jm < njrows; jm++, jac1 += 2 * njcols) {
            for (int iy = 0; iy < order; iy++) {
              for (int ix = 0; ix < order; ix++) {
                double Ni = na[ix] * nb[iy];
                double Nai = dna[ix] * nb[iy];
                double Nbi = na[ix] * dnb[iy];

                double Na1 = Nj * Nai;
                double Nb1 = Nj * Nbi;
                double Na2 = Ni * Naj;
                double Nb2 = Ni * Nbj;

                const TacsScalar *jac2 = &jacT[2 * jm];
                for (int im = 0; im < njcols; im++, jac2 += 2 * njcols) {
                  mat[im] += Na1 * jac1[2 * im] + Nb1 * jac1[2 * im + 1] +
                             Na2 * jac2[0] + Nb2 * jac2[1];
                }

                mat += nbcols;
              }
            }
          }

          mat += (nbrows - njrows) * ncols;
        }
      }
    } else if (jac) {
      for (int jy = 0; jy < order; jy++) {
        for (int jx = 0; jx < order; jx++) {
          double Nj = na[jx] * nb[jy];

          const TacsScalar *jac1 = jac;
          for (int jm = 0; jm < njrows; jm++, jac1 += 2 * njcols) {
            for (int iy = 0; iy < order; iy++) {
              for (int ix = 0; ix < order; ix++) {
                double Nai = dna[ix] * nb[iy];
                double Nbi = na[ix] * dnb[iy];

                double Na1 = Nj * Nai;
                double Nb1 = Nj * Nbi;

                for (int im = 0; im < njcols; im++) {
                  mat[im] += Na1 * jac1[2 * im] + Nb1 * jac1[2 * im + 1];
                }

                mat += nbcols;
              }
            }
          }

          mat += (nbrows - njrows) * ncols;
        }
      }
    } else if (jacT) {
      for (int jy = 0; jy < order; jy++) {
        for (int jx = 0; jx < order; jx++) {
          double Naj = dna[jx] * nb[jy];
          double Nbj = na[jx] * dnb[jy];

          for (int jm = 0; jm < njrows; jm++) {
            for (int iy = 0; iy < order; iy++) {
              for (int ix = 0; ix < order; ix++) {
                double Ni = na[ix] * nb[iy];
                double Na2 = Ni * Naj;
                double Nb2 = Ni * Nbj;

                const TacsScalar *jac2 = &jacT[2 * jm];
                for (int im = 0; im < njcols; im++, jac2 += 2 * njcols) {
                  mat[im] += Na2 * jac2[0] + Nb2 * jac2[1];
                }

                mat += nbcols;
              }
            }
          }

          mat += (nbrows - njrows) * ncols;
        }
      }
    }
  }

  /**
    Get the tying strain field index associated with this tying point

    @param index The index of the tying point
    @return The tying strain field index
  */
  static inline TacsShellTyingStrainComponent getTyingField(int index) {
    if (index < G11_OFFSET) {
      return TACS_SHELL_G11_COMPONENT;
    } else if (index < G22_OFFSET) {
      return TACS_SHELL_G22_COMPONENT;
    } else if (index < G12_OFFSET) {
      return TACS_SHELL_G12_COMPONENT;
    } else if (index < G13_OFFSET) {
      return TACS_SHELL_G23_COMPONENT;
    } else {
      return TACS_SHELL_G13_COMPONENT;
    }
  }

  /*
    Get the knots associated with the tying points
  */
  static inline void getTyingKnots(const double **ty_knots_order,
                                   const double **ty_knots_reduced) {
    if (order == 2) {
      *ty_knots_order = TacsShellLinearTyingPoints;
      *ty_knots_reduced = TacsGaussQuadPts1;
    } else if (order == 3) {
      *ty_knots_order = TacsGaussQuadPts3;
      *ty_knots_reduced = TacsGaussQuadPts2;
    } else if (order == 4) {
      *ty_knots_order = TacsGaussQuadPts4;
      *ty_knots_reduced = TacsGaussQuadPts4;
    } else if (order == 5) {
      *ty_knots_order = TacsGaussQuadPts5;
      *ty_knots_reduced = TacsGaussQuadPts4;
    } else {  // order == 6
      *ty_knots_order = TacsGaussQuadPts6;
      *ty_knots_reduced = TacsGaussQuadPts5;
    }
  }

  /**
    Get the tying point parametric location associated with the tying point

    @param index The index of the tying point
    @param pt The parametric point associated with the tying point
  */
  static inline void getTyingPoint(int ty_index, double pt[]) {
    const double *ty_knots_order, *ty_knots_reduced;
    getTyingKnots(&ty_knots_order, &ty_knots_reduced);

    int field = 0, ty = 0;
    if (ty_index < G11_OFFSET) {
      field = TACS_SHELL_G11_COMPONENT;
      ty = ty_index;
    } else if (ty_index < G22_OFFSET) {
      field = TACS_SHELL_G22_COMPONENT;
      ty = ty_index - G11_OFFSET;
    } else if (ty_index < G12_OFFSET) {
      field = TACS_SHELL_G12_COMPONENT;
      ty = ty_index - G22_OFFSET;
    } else if (ty_index < G13_OFFSET) {
      field = TACS_SHELL_G23_COMPONENT;
      ty = ty_index - G12_OFFSET;
    } else {
      field = TACS_SHELL_G13_COMPONENT;
      ty = ty_index - G13_OFFSET;
    }

    if (field == TACS_SHELL_G11_COMPONENT ||
        field == TACS_SHELL_G13_COMPONENT) {  // g11 or g13
      pt[0] = ty_knots_reduced[ty % (order - 1)];
      pt[1] = ty_knots_order[ty / (order - 1)];
    } else if (field == TACS_SHELL_G22_COMPONENT ||
               field == TACS_SHELL_G23_COMPONENT) {  // g22 or g23
      pt[0] = ty_knots_order[ty % order];
      pt[1] = ty_knots_reduced[ty / order];
    } else {  // (field == TACS_SHELL_G12_COMPONENT) g12
      pt[0] = ty_knots_reduced[ty % (order - 1)];
      pt[1] = ty_knots_reduced[ty / (order - 1)];
    }
  }

  /*
    Evaluate the interpolation for all of the tying points
  */
  static void evalTyingInterp(const double pt[], double N[]) {
    const double *ty_knots_order, *ty_knots_reduced;
    getTyingKnots(&ty_knots_order, &ty_knots_reduced);

    // Evaluate the required shape functions
    double na[order], nb[order];
    TacsLagrangeShapeFunction<order>(pt[0], ty_knots_order, na);
    TacsLagrangeShapeFunction<order>(pt[1], ty_knots_order, nb);

    double nar[order - 1], nbr[order - 1];
    TacsLagrangeShapeFunction<order - 1>(pt[0], ty_knots_reduced, nar);
    TacsLagrangeShapeFunction<order - 1>(pt[1], ty_knots_reduced, nbr);

    // TACS_SHELL_G11_COMPONENT
    for (int j = 0; j < order; j++) {
      for (int i = 0; i < order - 1; i++, N++) {
        N[0] = nar[i] * nb[j];
      }
    }

    // TACS_SHELL_G22_COMPONENT
    for (int j = 0; j < order - 1; j++) {
      for (int i = 0; i < order; i++, N++) {
        N[0] = na[i] * nbr[j];
      }
    }

    // TACS_SHELL_G12_COMPONENT
    for (int j = 0; j < order - 1; j++) {
      for (int i = 0; i < order - 1; i++, N++) {
        N[0] = nar[i] * nbr[j];
      }
    }

    // TACS_SHELL_G23_COMPONENT
    for (int j = 0; j < order - 1; j++) {
      for (int i = 0; i < order; i++, N++) {
        N[0] = na[i] * nbr[j];
      }
    }

    // TACS_SHELL_G13_COMPONENT
    for (int j = 0; j < order; j++) {
      for (int i = 0; i < order - 1; i++, N++) {
        N[0] = nar[i] * nb[j];
      }
    }
  }

  /*
    Get the number of tying points associated with each field
  */
  static inline int getNumTyingPoints(const int field) {
    if (field == TACS_SHELL_G11_COMPONENT) {
      return NUM_G11_TYING_POINTS;
    } else if (field == TACS_SHELL_G22_COMPONENT) {
      return NUM_G22_TYING_POINTS;
    } else if (field == TACS_SHELL_G12_COMPONENT) {
      return NUM_G12_TYING_POINTS;
    } else if (field == TACS_SHELL_G23_COMPONENT) {
      return NUM_G23_TYING_POINTS;
    } else if (field == TACS_SHELL_G13_COMPONENT) {
      return NUM_G13_TYING_POINTS;
    }
    return 0;
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
    const int index[] = {0, 3, 1, 4, 2};
    const int num_tying_fields = 5;

    double N[NUM_TYING_POINTS];
    evalTyingInterp(pt, N);

    const double *N0 = N;
    for (int field = 0; field < num_tying_fields; field++) {
      const int npts = getNumTyingPoints(field);

      gty[index[field]] = 0.0;
      for (int k = 0; k < npts; k++, N0++, ety++) {
        gty[index[field]] += N0[0] * ety[0];
      }
    }

    // Set the last tying strain entry to zero
    gty[5] = 0.0;
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
    const int index[] = {0, 3, 1, 4, 2};
    const int num_tying_fields = 5;

    double N[NUM_TYING_POINTS];
    evalTyingInterp(pt, N);

    const double *N0 = N;
    for (int field = 0; field < num_tying_fields; field++) {
      const int npts = getNumTyingPoints(field);

      for (int k = 0; k < npts; k++, N0++, dety++) {
        dety[0] += dgty[index[field]] * N0[0];
      }
    }
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
    // Set the values into the strain tensor
    const int index[] = {0, 3, 1, 4, 2};
    const int num_strains = 6;
    const int num_tying_fields = 5;

    double N[NUM_TYING_POINTS];
    evalTyingInterp(pt, N);

    const double *N1 = N;
    for (int field1 = 0; field1 < num_tying_fields; field1++) {
      const int npts1 = getNumTyingPoints(field1);
      for (int k1 = 0; k1 < npts1; k1++, N1++) {
        const double *N2 = N;
        for (int field2 = 0; field2 < num_tying_fields; field2++) {
          const int npts2 = getNumTyingPoints(field2);
          const TacsScalar value =
              d2gty[num_strains * index[field1] + index[field2]];

          for (int k2 = 0; k2 < npts2; k2++, N2++, d2ety++) {
            d2ety[0] += N1[0] * N2[0] * value;
          }
        }
      }
    }
  }
};

#endif  // TACS_SHELL_ELEMENT_QUAD_BASIS_H
