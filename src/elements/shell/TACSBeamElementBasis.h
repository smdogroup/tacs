#ifndef TACS_BEAM_ELEMENT_BASIS_H
#define TACS_BEAM_ELEMENT_BASIS_H

#include "TACSElementAlgebra.h"
#include "TACSElementTypes.h"
#include "TACSGaussQuadrature.h"
#include "TACSShellElementQuadBasis.h"

enum TacsBeamTyingStrainComponent {
  TACS_BEAM_G12_COMPONENT = 0,
  TACS_BEAM_G13_COMPONENT = 1
};

template <int order>
class TACSBeamBasis {
 public:
  static const int NUM_NODES = order;
  static const int NUM_G12_TYING_POINTS = (order - 1);
  static const int NUM_G13_TYING_POINTS = (order - 1);
  static const int NUM_TYING_POINTS =
      (NUM_G12_TYING_POINTS + NUM_G13_TYING_POINTS);

  static const int G12_OFFSET = NUM_G12_TYING_POINTS;

  static void getNodePoint(const int n, double pt[]) { pt[0] = -1.0 + 2.0 * n; }
  static ElementLayout getLayoutType() {
    if (order == 2) {
      return TACS_LINE_ELEMENT;
    } else if (order == 3) {
      return TACS_LINE_QUADRATIC_ELEMENT;
    } else if (order == 4) {
      return TACS_LINE_CUBIC_ELEMENT;
    }

    return TACS_LAYOUT_NONE;
  }

  template <int vars_per_node, int m>
  static void interpFields(const double pt[], const TacsScalar values[],
                           TacsScalar field[]) {
    double na[order];
    TacsLagrangeLobattoShapeFunction<order>(pt[0], na);

    for (int k = 0; k < m; k++) {
      field[k] = 0.0;
    }

    for (int i = 0; i < order; i++) {
      for (int k = 0; k < m; k++) {
        field[k] += na[i] * values[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void addInterpFieldsTranspose(const double pt[],
                                       const TacsScalar field[],
                                       TacsScalar values[]) {
    double na[order];
    TacsLagrangeLobattoShapeFunction<order>(pt[0], na);

    for (int i = 0; i < order; i++) {
      for (int k = 0; k < m; k++) {
        values[k] += na[i] * field[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void interpFieldsGrad(const double pt[], const TacsScalar values[],
                               TacsScalar grad[]) {
    double na[order], dna[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);

    for (int k = 0; k < m; k++) {
      grad[k] = 0.0;
    }

    for (int i = 0; i < order; i++) {
      for (int k = 0; k < m; k++) {
        grad[k] += dna[i] * values[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void addInterpFieldsGradTranspose(const double pt[], TacsScalar grad[],
                                           TacsScalar values[]) {
    double na[order], dna[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);

    for (int i = 0; i < order; i++) {
      for (int k = 0; k < m; k++) {
        values[k] += dna[i] * grad[k];
      }
      values += vars_per_node;
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
    double na[order];
    TacsLagrangeLobattoShapeFunction<order>(pt[0], na);

    const int ncols = NUM_NODES * nbcols;

    for (int jx = 0; jx < order; jx++) {
      const TacsScalar *jac1 = jac;
      for (int jm = 0; jm < njrows; jm++, jac1 += njcols) {
        for (int ix = 0; ix < order; ix++) {
          for (int im = 0; im < njcols; im++) {
            mat[im] += na[ix] * na[jx] * jac1[im];
          }

          mat += nbcols;
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
    @param jac The 2m x 2m Jacobian matrix of coefficients
    @param mat The element matrix
  */
  template <int nbrows, int nbcols, int njrows, int njcols>
  static void addInterpGradOuterProduct(const double pt[],
                                        const TacsScalar jac[],
                                        TacsScalar *mat) {
    double na[order], dna[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);

    const int ncols = NUM_NODES * nbcols;

    for (int jx = 0; jx < order; jx++) {
      const TacsScalar *jac1 = jac;
      for (int jm = 0; jm < njrows; jm++, jac1 += njcols) {
        for (int ix = 0; ix < order; ix++) {
          for (int im = 0; im < njcols; im++) {
            mat[im] += dna[ix] * dna[jx] * jac1[im];
          }
          mat += nbcols;
        }

        mat += (nbrows - njrows) * ncols;
      }
    }
  }

  /*
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
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);

    const int ncols = NUM_NODES * nbcols;

    // if (jac && jacT){
    //   for ( int jx = 0; jx < order; jx++ ){
    //     const TacsScalar *jac1 = jac;
    //     for ( int jm = 0; jm < njrows; jm++, jac1 += njcols ){
    //       for ( int ix = 0; ix < order; ix++ ){
    //         const TacsScalar *jac2 = &jacT[2*jm];
    //         for ( int im = 0; im < njcols; im++, jac2 += njcols ){
    //           mat[im] += na[]
    //             Na1*jac1[2*im] + Nb1*jac1[2*im+1] +
    //             Na2*jac2[0] + Nb2*jac2[1];
    //         }

    //       mat += nbcols;
    //     }

    //     mat += (nbrows - njrows)*ncols;
    //   }
    // }
  }

  /**
    Get the tying strain field index associated with this tying point

    @param index The index of the tying point
    @return The tying strain field index
  */
  static inline TacsBeamTyingStrainComponent getTyingField(int index) {
    if (index < G12_OFFSET) {
      return TACS_BEAM_G12_COMPONENT;
    } else {
      return TACS_BEAM_G13_COMPONENT;
    }
  }

  /*
    Get the knots associated with the tying points
  */
  static inline void getTyingKnots(const double **ty_knots) {
    if (order == 2) {
      *ty_knots = TacsGaussQuadPts1;
    } else if (order == 3) {
      *ty_knots = TacsGaussQuadPts2;
    } else if (order == 4) {
      *ty_knots = TacsGaussQuadPts4;
    } else if (order == 5) {
      *ty_knots = TacsGaussQuadPts4;
    } else {  // order == 6
      *ty_knots = TacsGaussQuadPts5;
    }
  }

  /**
    Get the tying point parametric location associated with the tying point

    @param index The index of the tying point
    @param pt The parametric point associated with the tying point
  */
  static inline void getTyingPoint(int ty_index, double pt[]) {
    const double *ty_knots;
    getTyingKnots(&ty_knots);

    if (ty_index < G12_OFFSET) {
      pt[0] = ty_knots[ty_index];
    } else {
      pt[0] = ty_knots[ty_index - G12_OFFSET];
    }
  }

  /*
    Evaluate the interpolation for all of the tying points
  */
  static void evalTyingInterp(const double pt[], double N[]) {
    const double *ty_knots;
    getTyingKnots(&ty_knots);

    double nar[order - 1];
    TacsLagrangeShapeFunction<order - 1>(pt[0], ty_knots, nar);

    // TACS_SHELL_G12_COMPONENT
    for (int i = 0; i < order - 1; i++, N++) {
      N[0] = nar[i];
    }

    // TACS_SHELL_G13_COMPONENT
    for (int i = 0; i < order - 1; i++, N++) {
      N[0] = nar[i];
    }
  }

  /*
    Get the number of tying points associated with each field
  */
  static inline int getNumTyingPoints(const int field) {
    if (field == TACS_BEAM_G12_COMPONENT) {
      return NUM_G12_TYING_POINTS;
    } else if (field == TACS_BEAM_G13_COMPONENT) {
      return NUM_G13_TYING_POINTS;
    }
    return 0;
  }

  /**
    Evaluate the tensorial components of the strain tensor at the
    specific quadrature point

    gty = [g12  g13]

    @param pt The quadrature point
    @param ety The strain computed at the tying points
    @param gty The interpolated tying strain
  */
  static inline void interpTyingStrain(const double pt[],
                                       const TacsScalar ety[],
                                       TacsScalar gty[]) {
    const int num_tying_fields = 2;

    double N[NUM_TYING_POINTS];
    evalTyingInterp(pt, N);

    const double *N0 = N;
    for (int field = 0; field < num_tying_fields; field++) {
      const int npts = getNumTyingPoints(field);

      gty[field] = 0.0;
      for (int k = 0; k < npts; k++, N0++, ety++) {
        gty[field] += N0[0] * ety[0];
      }
    }
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
    const int num_tying_fields = 2;

    double N[NUM_TYING_POINTS];
    evalTyingInterp(pt, N);

    const double *N0 = N;
    for (int field = 0; field < num_tying_fields; field++) {
      const int npts = getNumTyingPoints(field);

      for (int k = 0; k < npts; k++, N0++, dety++) {
        dety[0] += dgty[field] * N0[0];
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
    const int num_tying_fields = 2;

    double N[NUM_TYING_POINTS];
    evalTyingInterp(pt, N);

    const double *N1 = N;
    for (int field1 = 0; field1 < num_tying_fields; field1++) {
      const int npts1 = getNumTyingPoints(field1);
      for (int k1 = 0; k1 < npts1; k1++, N1++) {
        const double *N2 = N;
        for (int field2 = 0; field2 < num_tying_fields; field2++) {
          const int npts2 = getNumTyingPoints(field2);
          const TacsScalar value = d2gty[num_tying_fields * field1 + field2];

          for (int k2 = 0; k2 < npts2; k2++, N2++, d2ety++) {
            d2ety[0] += N1[0] * N2[0] * value;
          }
        }
      }
    }
  }
};

#endif  // TACS_BEAM_ELEMENT_BASIS_H
