#ifndef TACS_BEAM_ELEMENT_BASIS_H
#define TACS_BEAM_ELEMENT_BASIS_H

#include "TACSGaussQuadrature.h"
#include "TACSElementAlgebra.h"
#include "TACSElementTypes.h"

/**
  Defines the quadrature over both the face and quadrature
*/
class TACSBeamLinearQuadrature {
 public:
  static int getNumParameters(){
    return 1;
  }
  static int getNumQuadraturePoints(){
    return 2;
  }
  static double getQuadratureWeight( int n ){
    return TacsGaussQuadWts2[n];
  }
  static double getQuadraturePoint( int n, double pt[] ){
    pt[0] = TacsGaussQuadPts2[n];

    return TacsGaussQuadWts2[n];
  }
  static int getNumFaces(){
    return 2;
  }
  static int getNumFaceQuadraturePoints( int face ){
    return 1;
  }
  static double getFaceQuadraturePoint( int face, int n,
                                        double pt[],
                                        double t[] ){
    pt[0] = -1.0 + 2.0*face;
    return 1.0;
  }
};

/**
  Defines the quadrature over both the face and quadrature
*/
class TACSBeamQuadraticQuadrature {
 public:
  static int getNumParameters(){
    return 1;
  }
  static int getNumQuadraturePoints(){
    return 3;
  }
  static double getQuadratureWeight( int n ){
    return TacsGaussQuadWts3[n];
  }
  static double getQuadraturePoint( int n, double pt[] ){
    pt[0] = TacsGaussQuadPts3[n];

    return TacsGaussQuadWts3[n];
  }
  static int getNumFaces(){
    return 2;
  }
  static int getNumFaceQuadraturePoints( int face ){
    return 1;
  }
  static double getFaceQuadraturePoint( int face, int n,
                                        double pt[],
                                        double t[] ){
    pt[0] = -1.0 + 2.0*face;
    return 1.0;
  }
};

template <int order>
class TACSBeamBasis {
 public:
  static const int NUM_NODES = order;
  static const int NUM_TYING_POINTS = 2*(order - 1);

  static void getNodePoint( const int n, double pt[] ){
    pt[0] = -1.0 + 2.0*n;
  }
  static ElementLayout getLayoutType(){
    return TACS_LINE_ELEMENT;
  }

  template <int vars_per_node, int m>
  static void interpFields( const double pt[],
                            const TacsScalar values[],
                            TacsScalar field[] ){
    double na[order];
    TacsLagrangeLobattoShapeFunction<order>(pt[0], na);

    for ( int k = 0; k < m; k++ ){
      field[k] = 0.0;
    }

    for ( int i = 0; i < 2; i++ ){
      for ( int k = 0; k < m; k++ ){
        field[k] += na[i]*values[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void addInterpFieldsTranspose( const double pt[],
                                        const TacsScalar field[],
                                        TacsScalar values[] ){
    double na[order];
    TacsLagrangeLobattoShapeFunction<order>(pt[0], na);

    for ( int i = 0; i < 2; i++ ){
      for ( int k = 0; k < m; k++ ){
        values[k] += na[i]*field[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void interpFieldsGrad( const double pt[],
                                const TacsScalar values[],
                                TacsScalar grad[] ){
    double na[order], dna[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);

    for ( int k = 0; k < m; k++ ){
      grad[k] = 0.0;
    }

    for ( int i = 0; i < 2; i++ ){
      for ( int k = 0; k < m; k++ ){
        grad[k] += dna[i]*values[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void addInterpFieldsGradTranspose( const double pt[],
                                            TacsScalar grad[],
                                            TacsScalar values[] ){
    double na[order], dna[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);

    for ( int i = 0; i < 2; i++ ){
      for ( int k = 0; k < m; k++ ){
        values[k] += dna[i]*grad[k];
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
  static void addInterpFieldsOuterProduct( const double pt[],
                                           const TacsScalar jac[],
                                           TacsScalar *mat ){
    double na[order];
    TacsLagrangeLobattoShapeFunction<order>(pt[0], na);

    const int ncols = NUM_NODES*nbcols;

    for ( int jx = 0; jx < order; jx++ ){
      const TacsScalar *jac1 = jac;
      for ( int jm = 0; jm < njrows; jm++, jac1 += njcols ){
        for ( int ix = 0; ix < order; ix++ ){
          for ( int im = 0; im < njcols; im++ ){
            mat[im] += na[ix]*na[jx]*jac1[im];
          }

          mat += nbcols;
        }

        mat += (nbrows - njrows)*ncols;
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
  static void addInterpGradOuterProduct( const double pt[],
                                         const TacsScalar jac[],
                                         TacsScalar *mat ){
    double na[order], dna[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);

    const int ncols = NUM_NODES*nbcols;

    for ( int jx = 0; jx < order; jx++ ){
      const TacsScalar *jac1 = jac;
      for ( int jm = 0; jm < njrows; jm++, jac1 += njcols ){
        for ( int ix = 0; ix < order; ix++ ){
          for ( int im = 0; im < njcols; im++ ){
            mat[im] += dna[ix]*dna[jx]*jac1[im];
          }
          mat += nbcols;
        }

        mat += (nbrows - njrows)*ncols;
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
  static void addInterpGradMixedOuterProduct( const double pt[],
                                              const TacsScalar jac[],
                                              const TacsScalar jacT[],
                                              TacsScalar *mat ){
    double na[order], dna[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);

    const int ncols = NUM_NODES*nbcols;

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

  static inline void getTyingKnots( const double **ty_knots ){
    if (order == 2){
      *ty_knots = TacsGaussQuadPts1;
    }
    else if (order == 3){
      *ty_knots = TacsGaussQuadPts2;
    }
    else if (order == 4){
      *ty_knots = TacsGaussQuadPts4;
    }
    else if (order == 5){
      *ty_knots = TacsGaussQuadPts4;
    }
    else { // order == 6
      *ty_knots = TacsGaussQuadPts5;
    }
  }
  static int getNumTyingFields(){
    return 2;
  }
  static int getNumTyingPoints( const int field ){
    if (field == 0 || field == 1){ return order-1; }
    return 0;
  }
  static void getTyingPoint( const int field,
                             const int ty,
                             double pt[] ){
    const double *ty_knots;
    getTyingKnots(&ty_knots);
    pt[0] = ty_knots[ty];
  }
  static TacsScalar interpTying( const int field,
                                 const double pt[],
                                 const TacsScalar ety[] ){
    const double *ty_knots;
    getTyingKnots(&ty_knots);

    TacsScalar value = 0.0;
    double na[order-1];
    TacsLagrangeShapeFunction<order-1>(pt[0], ty_knots, na);

    for ( int i = 0; i < order-1; i++, ety++ ){
      value += na[i]*ety[0];
    }

    return value;
  }

  static void addInterpTyingTranspose( const int field,
                                       const double pt[],
                                       const TacsScalar value,
                                       TacsScalar ety[] ){
    const double *ty_knots;
    getTyingKnots(&ty_knots);

    double na[order-1];
    TacsLagrangeShapeFunction<order-1>(pt[0], ty_knots, na);

    for ( int i = 0; i < order-1; i++, ety++ ){
      ety[0] += na[i]*value;
    }
  }

  static void addInterpTyingOuterProduct( const int f1,
                                          const int f2,
                                          const double pt[],
                                          const TacsScalar value,
                                          TacsScalar d2ety[] ){
    const double *ty_knots;
    getTyingKnots(&ty_knots);

    double na[order-1];
    TacsLagrangeShapeFunction<order-1>(pt[0], ty_knots, na);

    for ( int i = 0; i < order-1; i++ ){
      TacsScalar N = value*na[i];

      for ( int j = 0; j < order-1; j++ ){
        d2ety[0] += na[j]*N;
      }
    }
  }
};

#endif // TACS_BEAM_ELEMENT_BASIS_H