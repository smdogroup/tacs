#ifndef TACS_SHELL_ELEMENT_BASIS_H
#define TACS_SHELL_ELEMENT_BASIS_H

#include "TACSGaussQuadrature.h"
#include "TACSTriangleQuadrature.h"
#include "TACSLagrangeInterpolation.h"
#include "TACSElementAlgebra.h"
#include "TACSElementTypes.h"

/**
  Defines the quadrature over both the face and quadrature
*/
class TACSQuadLinearQuadrature {
 public:
  static const int NUM_QUADRATURE_POINTS = 4;

  static int getNumParameters(){
    return 2;
  }
  static int getNumQuadraturePoints(){
    return 4;
  }
  static double getQuadratureWeight( int n ){
    return TacsGaussQuadWts2[n % 2]*TacsGaussQuadWts2[n / 2];
  }
  static double getQuadraturePoint( int n, double pt[] ){
    pt[0] = TacsGaussQuadPts2[n % 2];
    pt[1] = TacsGaussQuadPts2[n / 2];

    return TacsGaussQuadWts2[n % 2]*TacsGaussQuadWts2[n/2];
  }
  static int getNumElementFaces(){
    return 4;
  }
  static int getNumFaceQuadraturePoints( int face ){
    return 2;
  }
  static double getFaceQuadraturePoint( int face, int n,
                                        double pt[],
                                        double t[] ){
    if (face/2 == 0){
      pt[0] = -1.0 + 2.0*(face % 2);
      pt[1] = TacsGaussQuadPts2[n];
    }
    else {
      pt[0] = TacsGaussQuadPts2[n];
      pt[1] = -1.0 + 2.0*(face % 2);
    }

    if (face == 0){
      // -X edge
      t[0] = 0.0;  t[1] = -1.0;
    }
    else if (face == 1){
      // +X edge
      t[0] = 0.0;  t[1] = 1.0;
    }
    else if (face == 2){
      // -Y edge
      t[0] = 1.0;  t[1] = 0.0;
    }
    else if (face == 3){
      // +Y edge
      t[0] = -1.0;  t[1] = 0.0;
    }

    return TacsGaussQuadWts2[n];
  }
};

class TACSQuadQuadraticQuadrature {
 public:
  static const int NUM_QUADRATURE_POINTS = 9;

  static int getNumParameters(){
    return 2;
  }
  static int getNumQuadraturePoints(){
    return 9;
  }
  static double getQuadratureWeight( int n ){
    return TacsGaussQuadWts3[n % 3]*TacsGaussQuadWts3[n / 3];
  }
  static double getQuadraturePoint( int n, double pt[] ){
    pt[0] = TacsGaussQuadPts3[n % 3];
    pt[1] = TacsGaussQuadPts3[n / 3];

    return TacsGaussQuadWts3[n % 3]*TacsGaussQuadWts3[n/3];
  }
  static int getNumElementFaces(){
    return 4;
  }
  static int getNumFaceQuadraturePoints( int face ){
    return 3;
  }
  static double getFaceQuadraturePoint( int face, int n,
                                        double pt[],
                                        double t[] ){
    if (face/2 == 0){
      pt[0] = -1.0 + 2.0*(face % 2);
      pt[1] = TacsGaussQuadPts3[n];
    }
    else {
      pt[0] = TacsGaussQuadPts3[n];
      pt[1] = -1.0 + 2.0*(face % 2);
    }

    if (face == 0){
      // -X edge
      t[0] = 0.0;  t[1] = -1.0;
    }
    else if (face == 1){
      // +X edge
      t[0] = 0.0;  t[1] = 1.0;
    }
    else if (face == 2){
      // -Y edge
      t[0] = 1.0;  t[1] = 0.0;
    }
    else if (face == 3){
      // +Y edge
      t[0] = -1.0;  t[1] = 0.0;
    }

    return TacsGaussQuadWts3[n];
  }
};


class TACSQuadCubicQuadrature {
 public:
  static const int NUM_QUADRATURE_POINTS = 16;

  static int getNumParameters(){
    return 2;
  }
  static int getNumQuadraturePoints(){
    return 16;
  }
  static double getQuadratureWeight( int n ){
    return TacsGaussQuadWts4[n % 4]*TacsGaussQuadWts4[n / 4];
  }
  static double getQuadraturePoint( int n, double pt[] ){
    pt[0] = TacsGaussQuadPts4[n % 4];
    pt[1] = TacsGaussQuadPts4[n / 4];

    return TacsGaussQuadWts4[n % 4]*TacsGaussQuadWts4[n / 4];
  }
  static int getNumElementFaces(){
    return 4;
  }
  static int getNumFaceQuadraturePoints( int face ){
    return 4;
  }
  static double getFaceQuadraturePoint( int face, int n,
                                        double pt[],
                                        double t[] ){
    if (face/2 == 0){
      pt[0] = -1.0 + 2.0*(face % 2);
      pt[1] = TacsGaussQuadPts4[n];
    }
    else {
      pt[0] = TacsGaussQuadPts4[n];
      pt[1] = -1.0 + 2.0*(face % 2);
    }

    if (face == 0){
      // -X edge
      t[0] = 0.0;  t[1] = -1.0;
    }
    else if (face == 1){
      // +X edge
      t[0] = 0.0;  t[1] = 1.0;
    }
    else if (face == 2){
      // -Y edge
      t[0] = 1.0;  t[1] = 0.0;
    }
    else if (face == 3){
      // +Y edge
      t[0] = -1.0;  t[1] = 0.0;
    }

    return TacsGaussQuadWts4[n];
  }
};

class TACSTriLinearQuadrature {
 public:
  static const int NUM_QUADRATURE_POINTS = 3;

  static int getNumParameters(){
    return 2;
  }
  static int getNumQuadraturePoints(){
    return 3;
  }
  static double getQuadratureWeight( int n ){
    return TacsTriangleWts3[n];
  }
  static double getQuadraturePoint( int n, double pt[] ){
    pt[0] = TacsTrianglePts3[2*n];
    pt[1] = TacsTrianglePts3[2*n+1];

    return TacsTriangleWts3[n];
  }
  static int getNumElementFaces(){
    return 3;
  }
  static int getNumFaceQuadraturePoints( int face ){
    return 2;
  }
  static double getFaceQuadraturePoint( int face, int n,
                                        double pt[],
                                        double t[] ){

    if (face == 0){
      t[0] = 1.0;
      t[1] = 0.0;
    }
    else if (face == 1){
      t[0] = -1.0;
      t[1] = 1.0;
    }
    else if (face == 2){
      t[0] = 0.0;
      t[1] = 1.0;
    }

    if (face == 0){
      pt[0] = 0.5*TacsGaussQuadPts2[n] + 0.5;
      pt[1] = 0.0;
      return 0.5*TacsGaussQuadWts2[n];
    }
    else if (face == 2){
      pt[0] = 0.0;
      pt[1] = 0.5*TacsGaussQuadPts2[n] + 0.5;
      return 0.5*TacsGaussQuadWts2[n];
    }
    else if (face == 1){
      pt[0] = 1.0 - (0.5*TacsGaussQuadPts2[n] + 0.5);
      pt[1] = 0.5*TacsGaussQuadPts2[n] + 0.5;
      return 0.5*sqrt(2.0)*TacsGaussQuadWts2[n];
    }

    return 0.0;
  }
};


class TACSTriQuadraticQuadrature {
 public:
  static const int NUM_QUADRATURE_POINTS = 6;

  static int getNumParameters(){
    return 2;
  }
  static int getNumQuadraturePoints(){
    return 4;
  }
  static double getQuadratureWeight( int n ){
    return TacsTriangleWts4[n];
  }
  static double getQuadraturePoint( int n, double pt[] ){
    pt[0] = TacsTrianglePts4[2*n];
    pt[1] = TacsTrianglePts4[2*n+1];

    return TacsTriangleWts4[n];
  }
  static int getNumElementFaces(){
    return 3;
  }
  static int getNumFaceQuadraturePoints( int face ){
    return 2;
  }
  static double getFaceQuadraturePoint( int face, int n,
                                        double pt[],
                                        double t[] ){
    if (face == 0){
      t[0] = 1.0;
      t[1] = 0.0;
    }
    else if (face == 1){
      t[0] = -1.0;
      t[1] = 1.0;
    }
    else if (face == 2){
      t[0] = 0.0;
      t[1] = 1.0;
    }

    if (face == 0){
      pt[0] = 0.5*TacsGaussQuadPts2[n] + 0.5;
      pt[1] = 0.0;
      return 0.5*TacsGaussQuadWts2[n];
    }
    else if (face == 2){
      pt[0] = 0.0;
      pt[1] = 0.5*TacsGaussQuadPts2[n] + 0.5;
      return 0.5*TacsGaussQuadWts2[n];
    }
    else if (face == 1){
      pt[0] = 1.0 - (0.5*TacsGaussQuadPts2[n] + 0.5);
      pt[1] = 0.5*TacsGaussQuadPts2[n] + 0.5;
      return 0.5*sqrt(2.0)*TacsGaussQuadWts2[n];
    }

    return 0.0;
  }
};


template <int order>
inline void TacsLagrangeShapeFunction( const double u,
                                       const double knots[],
                                       double N[] ){
  // Loop over the shape functions
  for ( int i = 0; i < order; i++ ){
    N[i] = 1.0;
    for ( int j = 0; j < order; j++ ){
      if (i != j){
        double d = 1.0/(knots[i] - knots[j]);
        N[i] *= (u - knots[j])*d;
      }
    }
  }
}

template <int order>
inline void TacsLagrangeShapeFuncDerivative( const double u,
                                             const double knots[],
                                             double N[],
                                             double Nd[] ){
  // Loop over the shape function knot locations
  for ( int i = 0; i < order; i++ ){
    N[i] = 1.0;
    Nd[i] = 0.0;

    // Loop over each point again, except for the current control
    // point, adding the contribution to the shape function
    for ( int j = 0; j < order; j++ ){
      if (i != j){
        double d = 1.0/(knots[i] - knots[j]);
        N[i] *= (u - knots[j])*d;

        // Now add up the contribution to the derivative
        for ( int k = 0; k < order; k++ ){
          if (k != i && k != j){
            d *= (u - knots[k])/(knots[i] - knots[k]);
          }
        }

        // Add the derivative contribution
        Nd[i] += d;
      }
    }
  }
}

template <int order>
inline void TacsLagrangeLobattoShapeFunction( const double u,
                                              double *N ){
  if (order == 1){
    N[0] = 1.0;
  }
  else if (order == 2){
    N[0] = 0.5*(1.0 - u);
    N[1] = 0.5*(1.0 + u);
  }
  else if (order == 3){
    N[0] = -0.5*u*(1.0 - u);
    N[1] = (1.0 - u)*(1.0 + u);
    N[2] = 0.5*(1.0 + u)*u;
  }
  else {
    const double *knots = TacsGaussLobattoPoints4;
    if (order == 5){
      knots = TacsGaussLobattoPoints5;
    }
    else if (order == 6){
      knots = TacsGaussLobattoPoints6;
    }

    TacsLagrangeShapeFunction<order>(u, knots, N);
  }
}

template <int order>
inline void TacsLagrangeLobattoShapeFuncDerivative( const double u,
                                                    double *N,
                                                    double *Nd ){
  if (order == 1){
    N[0] = 1.0;
  }
  else if (order == 2){
    N[0] = 0.5*(1.0 - u);
    N[1] = 0.5*(1.0 + u);

    Nd[0] = -0.5;
    Nd[1] = 0.5;
  }
  else if (order == 3){
    N[0] = -0.5*u*(1.0 - u);
    N[1] = (1.0 - u)*(1.0 + u);
    N[2] = 0.5*(1.0 + u)*u;

    Nd[0] = -0.5 + u;
    Nd[1] = -2.0*u;
    Nd[2] = 0.5 + u;
  }
  else {
    const double *knots = TacsGaussLobattoPoints4;
    if (order == 5){
      knots = TacsGaussLobattoPoints5;
    }
    else if (order == 6){
      knots = TacsGaussLobattoPoints6;
    }

    TacsLagrangeShapeFuncDerivative<order>(u, knots, N, Nd);
  }
}

const double TacsShellLinearTyingPoints[2] = {-1.0, 1.0};

template <int order>
class TACSShellQuadBasis {
 public:
  static const int NUM_NODES = order*order;

  // Set the number of tying points for each of the 5 components
  // of the tying strain
  static const int NUM_G11_TYING_POINTS = order*(order - 1);
  static const int NUM_G22_TYING_POINTS = order*(order - 1);
  static const int NUM_G12_TYING_POINTS = (order - 1)*(order - 1);
  static const int NUM_G13_TYING_POINTS = order*(order - 1);
  static const int NUM_G23_TYING_POINTS = order*(order - 1);

  static const int NUM_TYING_POINTS =
    NUM_G11_TYING_POINTS +
    NUM_G22_TYING_POINTS +
    NUM_G12_TYING_POINTS +
    NUM_G13_TYING_POINTS +
    NUM_G23_TYING_POINTS;

  static void getNodePoint( const int n, double pt[] ){
    pt[0] = -1.0 + (2.0/(order - 1))*(n % order);
    pt[1] = -1.0 + (2.0/(order - 1))*(n / order);
  }
  static ElementLayout getLayoutType(){
    if (order == 2){
      return TACS_QUAD_ELEMENT;
    }
    else if (order == 3){
      return TACS_QUAD_QUADRATIC_ELEMENT;
    }
    else if (order == 4){
      return TACS_QUAD_CUBIC_ELEMENT;
    }
    else if (order == 5){
      return TACS_QUAD_QUARTIC_ELEMENT;
    }
    else if (order == 6){
      return TACS_QUAD_QUINTIC_ELEMENT;
    }

    return TACS_LAYOUT_NONE;
  }

  template <int vars_per_node, int m>
  static void interpFields( const double pt[],
                            const TacsScalar values[],
                            TacsScalar field[] ){
    double na[order], nb[order];
    TacsLagrangeLobattoShapeFunction<order>(pt[0], na);
    TacsLagrangeLobattoShapeFunction<order>(pt[1], nb);

    for ( int k = 0; k < m; k++ ){
      field[k] = 0.0;
    }

    for ( int j = 0; j < order; j++ ){
      for ( int i = 0; i < order; i++ ){
        for ( int k = 0; k < m; k++ ){
          field[k] += na[i]*nb[j]*values[k];
        }
        values += vars_per_node;
      }
    }
  }

  template <int vars_per_node, int m>
  static void addInterpFieldsTranspose( const double pt[],
                                        const TacsScalar field[],
                                        TacsScalar values[] ){
    double na[order], nb[order];
    TacsLagrangeLobattoShapeFunction<order>(pt[0], na);
    TacsLagrangeLobattoShapeFunction<order>(pt[1], nb);

    for ( int j = 0; j < order; j++ ){
      for ( int i = 0; i < order; i++ ){
        for ( int k = 0; k < m; k++ ){
          values[k] += na[i]*nb[j]*field[k];
        }
        values += vars_per_node;
      }
    }
  }

  template <int vars_per_node, int m>
  static void interpFieldsGrad( const double pt[],
                                const TacsScalar values[],
                                TacsScalar grad[] ){
    double na[order], dna[order];
    double nb[order], dnb[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[1], nb, dnb);

    for ( int k = 0; k < m; k++ ){
      grad[2*k] = 0.0;
      grad[2*k+1] = 0.0;
    }

    for ( int j = 0; j < order; j++ ){
      for ( int i = 0; i < order; i++ ){
        for ( int k = 0; k < m; k++ ){
          grad[2*k]   += dna[i]*nb[j]*values[k];
          grad[2*k+1] += na[i]*dnb[j]*values[k];
        }
        values += vars_per_node;
      }
    }
  }

  template <int vars_per_node, int m>
  static void addInterpFieldsGradTranspose( const double pt[],
                                            TacsScalar grad[],
                                            TacsScalar values[] ){
    double na[order], dna[order];
    double nb[order], dnb[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[1], nb, dnb);

    for ( int j = 0; j < order; j++ ){
      for ( int i = 0; i < order; i++ ){
        for ( int k = 0; k < m; k++ ){
          values[k] += (dna[i]*nb[j]*grad[2*k] + na[i]*dnb[j]*grad[2*k+1]);
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
  static void addInterpFieldsOuterProduct( const double pt[],
                                           const TacsScalar jac[],
                                           TacsScalar *mat ){
    double na[order], nb[order];
    TacsLagrangeLobattoShapeFunction<order>(pt[0], na);
    TacsLagrangeLobattoShapeFunction<order>(pt[1], nb);

    const int ncols = NUM_NODES*nbcols;

    for ( int jy = 0; jy < order; jy++ ){
      for ( int jx = 0; jx < order; jx++ ){
        double Nj = na[jx]*nb[jy];

        const TacsScalar *jac1 = jac;
        for ( int jm = 0; jm < njrows; jm++, jac1 += njcols ){
          for ( int iy = 0; iy < order; iy++ ){
            for ( int ix = 0; ix < order; ix++ ){
              double Ni = na[ix]*nb[iy];

              for ( int im = 0; im < njcols; im++ ){
                mat[im] += Ni*Nj*jac1[im];
              }

              mat += nbcols;
            }
          }
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
    double nb[order], dnb[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[1], nb, dnb);

    const int ncols = NUM_NODES*nbcols;

    for ( int jy = 0; jy < order; jy++ ){
      for ( int jx = 0; jx < order; jx++ ){
        double Naj = dna[jx]*nb[jy];
        double Nbj = na[jx]*dnb[jy];

        const TacsScalar *jac1 = jac;
        const TacsScalar *jac2 = &jac[2*njcols];
        for ( int jm = 0; jm < njrows; jm++, jac1 += 4*njcols, jac2 += 4*njcols ){
          for ( int iy = 0; iy < order; iy++ ){
            for ( int ix = 0; ix < order; ix++ ){
              double Nai = dna[ix]*nb[iy];
              double Nbi = na[ix]*dnb[iy];
              double Naa = Naj*Nai;
              double Nab = Naj*Nbi;
              double Nba = Nbj*Nai;
              double Nbb = Nbj*Nbi;

              for ( int im = 0; im < njcols; im++ ){
                mat[im] += (Naa*jac1[2*im] + Nab*jac1[2*im + 1] +
                  Nba*jac2[2*im] + Nbb*jac2[2*im + 1]);
              }

              mat += nbcols;
            }
          }
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
    double nb[order], dnb[order];
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[0], na, dna);
    TacsLagrangeLobattoShapeFuncDerivative<order>(pt[1], nb, dnb);

    const int ncols = NUM_NODES*nbcols;

    if (jac && jacT){
      for ( int jy = 0; jy < order; jy++ ){
        for ( int jx = 0; jx < order; jx++ ){
          double Nj = na[jx]*nb[jy];
          double Naj = dna[jx]*nb[jy];
          double Nbj = na[jx]*dnb[jy];

          const TacsScalar *jac1 = jac;
          for ( int jm = 0; jm < njrows; jm++, jac1 += 2*njcols ){
            for ( int iy = 0; iy < order; iy++ ){
              for ( int ix = 0; ix < order; ix++ ){
                double Ni = na[ix]*nb[iy];
                double Nai = dna[ix]*nb[iy];
                double Nbi = na[ix]*dnb[iy];

                double Na1 = Nj*Nai;
                double Nb1 = Nj*Nbi;
                double Na2 = Ni*Naj;
                double Nb2 = Ni*Nbj;

                const TacsScalar *jac2 = &jacT[2*jm];
                for ( int im = 0; im < njcols; im++, jac2 += 2*njcols ){
                  mat[im] +=
                    Na1*jac1[2*im] + Nb1*jac1[2*im+1] +
                    Na2*jac2[0] + Nb2*jac2[1];
                }

                mat += nbcols;
              }
            }
          }

          mat += (nbrows - njrows)*ncols;
        }
      }
    }
    else if (jac){
      for ( int jy = 0; jy < order; jy++ ){
        for ( int jx = 0; jx < order; jx++ ){
          double Nj = na[jx]*nb[jy];

          const TacsScalar *jac1 = jac;
          for ( int jm = 0; jm < njrows; jm++, jac1 += 2*njcols ){

            for ( int iy = 0; iy < order; iy++ ){
              for ( int ix = 0; ix < order; ix++ ){
                double Nai = dna[ix]*nb[iy];
                double Nbi = na[ix]*dnb[iy];

                double Na1 = Nj*Nai;
                double Nb1 = Nj*Nbi;

                for ( int im = 0; im < njcols; im++ ){
                  mat[im] += Na1*jac1[2*im] + Nb1*jac1[2*im+1];
                }

                mat += nbcols;
              }
            }
          }

          mat += (nbrows - njrows)*ncols;
        }
      }
    }
    else if (jacT){
      for ( int jy = 0; jy < order; jy++ ){
        for ( int jx = 0; jx < order; jx++ ){
          double Naj = dna[jx]*nb[jy];
          double Nbj = na[jx]*dnb[jy];

          for ( int jm = 0; jm < njrows; jm++ ){
            for ( int iy = 0; iy < order; iy++ ){
              for ( int ix = 0; ix < order; ix++ ){
                double Ni = na[ix]*nb[iy];
                double Na2 = Ni*Naj;
                double Nb2 = Ni*Nbj;

                const TacsScalar *jac2 = &jacT[2*jm];
                for ( int im = 0; im < njcols; im++, jac2 += 2*njcols ){
                  mat[im] += Na2*jac2[0] + Nb2*jac2[1];
                }

                mat += nbcols;
              }
            }
          }

          mat += (nbrows - njrows)*ncols;
        }
      }
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
  static inline void interpTyingStrain( const double pt[],
                                        const TacsScalar ety[],
                                        TacsScalar gty[] ){
    // Set the values into the strain tensor
    const int index[] = {0, 3, 1, 4, 2};
    const int num_tying_fields = 5;
    for ( int field = 0; field < num_tying_fields; field++ ){
      gty[index[field]] = interpTying(field, pt, ety);
      ety += getNumTyingPoints(field);
    }
    gty[5] = 0.0;
  }

  /**
    Add the derivative of the tying strain to the residual

    @param pt The quadrature point
    @param dgty The derivative of the interpolated strain
    @param dety The output derivative of the strain at the tying points
  */
  static inline void addInterpTyingStrainTranspose( const double pt[],
                                                    const TacsScalar dgty[],
                                                    TacsScalar dety[] ){
    // Set the values into the strain tensor
    const int index[] = {0, 3, 1, 4, 2};
    const int num_tying_fields = 5;
    for ( int field = 0; field < num_tying_fields; field++ ){
      addInterpTyingTranspose(field, pt, dgty[index[field]], dety);
      dety += getNumTyingPoints(field);
    }
  }

  /**
    Add the second derivative of the tying strain at the tying points

    @param pt The quadrature point
    @param d2gty The second derivative of the interpolated strain
    @param d2ety The second derivatives of the strain at the tying points
  */
  static inline void addInterpTyingStrainHessian( const double pt[],
                                                  const TacsScalar d2gty[],
                                                  TacsScalar d2ety[] ){
    // Set the values into the strain tensor
    const int index[] = {0, 3, 1, 4, 2};
    const int num_strains = 6;
    const int num_tying_fields = 5;
    for ( int field1 = 0; field1 < num_tying_fields; field1++ ){
      for ( int field2 = 0; field2 < num_tying_fields; field2++ ){
        TacsScalar value = d2gty[num_strains*index[field1] + index[field2]];
        addInterpTyingOuterProduct(field1, field2, pt, value, d2ety);
        d2ety += getNumTyingPoints(field1)*getNumTyingPoints(field2);
      }
    }
  }

  static inline void getTyingKnots( const double **ty_knots_order,
                                    const double **ty_knots_reduced ){
    if (order == 2){
      *ty_knots_order = TacsShellLinearTyingPoints;
      *ty_knots_reduced = TacsGaussQuadPts1;
    }
    else if (order == 3){
      *ty_knots_order = TacsGaussQuadPts3;
      *ty_knots_reduced = TacsGaussQuadPts2;
    }
    else if (order == 4){
      *ty_knots_order = TacsGaussQuadPts4;
      *ty_knots_reduced = TacsGaussQuadPts4;
    }
    else if (order == 5){
      *ty_knots_order = TacsGaussQuadPts5;
      *ty_knots_reduced = TacsGaussQuadPts4;
    }
    else { // order == 6
      *ty_knots_order = TacsGaussQuadPts6;
      *ty_knots_reduced = TacsGaussQuadPts5;
    }
  }
  static int getNumTyingFields(){
    return 5;
  }
  static int getNumTyingPoints( const int field ){
    if (field == 0){ return NUM_G11_TYING_POINTS; }
    else if (field == 1){ return NUM_G22_TYING_POINTS; }
    else if (field == 2){ return NUM_G12_TYING_POINTS; }
    else if (field == 3){ return NUM_G13_TYING_POINTS; }
    else if (field == 4){ return NUM_G23_TYING_POINTS; }
    return 0;
  }
  static void getTyingPoint( const int field,
                             const int ty,
                             double pt[] ){
    const double *ty_knots_order, *ty_knots_reduced;
    getTyingKnots(&ty_knots_order, &ty_knots_reduced);

    if (field == 0 || field == 4){ // g11 or g13
      pt[0] = ty_knots_reduced[ty % (order - 1)];
      pt[1] = ty_knots_order[ty / (order - 1)];
    }
    else if (field == 1 || field == 3){ // g22 or g23
      pt[0] = ty_knots_order[ty % order];
      pt[1] = ty_knots_reduced[ty / order];
    }
    else { // (field == 2) g12
      pt[0] = ty_knots_reduced[ty % (order - 1)];
      pt[1] = ty_knots_reduced[ty / (order - 1)];
    }
  }
  static TacsScalar interpTying( const int field,
                                 const double pt[],
                                 const TacsScalar ety[] ){
    const double *ty_knots_order, *ty_knots_reduced;
    getTyingKnots(&ty_knots_order, &ty_knots_reduced);

    TacsScalar value = 0.0;
    if (field == 0 || field == 4){
      double na[order-1], nb[order];
      TacsLagrangeShapeFunction<order-1>(pt[0], ty_knots_reduced, na);
      TacsLagrangeShapeFunction<order>(pt[1], ty_knots_order, nb);

      for ( int j = 0; j < order; j++ ){
        for ( int i = 0; i < order-1; i++, ety++ ){
          value += na[i]*nb[j]*ety[0];
        }
      }
    }
    else if (field == 1 || field == 3){
      double na[order], nb[order-1];
      TacsLagrangeShapeFunction<order>(pt[0], ty_knots_order, na);
      TacsLagrangeShapeFunction<order-1>(pt[1], ty_knots_reduced, nb);

      for ( int j = 0; j < order-1; j++ ){
        for ( int i = 0; i < order; i++, ety++ ){
          value += na[i]*nb[j]*ety[0];
        }
      }
    }
    else { // field == 2
      double na[order-1], nb[order-1];
      TacsLagrangeShapeFunction<order-1>(pt[0], ty_knots_reduced, na);
      TacsLagrangeShapeFunction<order-1>(pt[1], ty_knots_reduced, nb);

      for ( int j = 0; j < order-1; j++ ){
        for ( int i = 0; i < order-1; i++, ety++ ){
          value += na[i]*nb[j]*ety[0];
        }
      }
    }

    return value;
  }

  static void addInterpTyingTranspose( const int field,
                                       const double pt[],
                                       const TacsScalar value,
                                       TacsScalar ety[] ){
    const double *ty_knots_order, *ty_knots_reduced;
    getTyingKnots(&ty_knots_order, &ty_knots_reduced);

    if (field == 0 || field == 4){
      double na[order-1], nb[order];
      TacsLagrangeShapeFunction<order-1>(pt[0], ty_knots_reduced, na);
      TacsLagrangeShapeFunction<order>(pt[1], ty_knots_order, nb);

      for ( int j = 0; j < order; j++ ){
        for ( int i = 0; i < order-1; i++, ety++ ){
          ety[0] += value*na[i]*nb[j];
        }
      }
    }
    else if (field == 1 || field == 3){
      double na[order], nb[order-1];
      TacsLagrangeShapeFunction<order>(pt[0], ty_knots_order, na);
      TacsLagrangeShapeFunction<order-1>(pt[1], ty_knots_reduced, nb);

      for ( int j = 0; j < order-1; j++ ){
        for ( int i = 0; i < order; i++, ety++ ){
          ety[0] += value*na[i]*nb[j];
        }
      }
    }
    else { // field == 2
      double na[order-1], nb[order-1];
      TacsLagrangeShapeFunction<order-1>(pt[0], ty_knots_reduced, na);
      TacsLagrangeShapeFunction<order-1>(pt[1], ty_knots_reduced, nb);

      for ( int j = 0; j < order-1; j++ ){
        for ( int i = 0; i < order-1; i++, ety++ ){
          ety[0] += value*na[i]*nb[j];
        }
      }
    }
  }

  static void addInterpTyingOuterProduct( const int f1,
                                          const int f2,
                                          const double pt[],
                                          const TacsScalar value,
                                          TacsScalar d2ety[] ){
    const double *ty_knots_order, *ty_knots_reduced;
    getTyingKnots(&ty_knots_order, &ty_knots_reduced);

    int ntu1, ntv1;
    double nu1[order], nv1[order];
    if (f1 == 0 || f1 == 4){
      ntu1 = order-1;
      ntv1 = order;
      TacsLagrangeShapeFunction<order-1>(pt[0], ty_knots_reduced, nu1);
      TacsLagrangeShapeFunction<order>(pt[1], ty_knots_order, nv1);
    }
    else if (f1 == 1 || f1 == 3){
      ntu1 = order;
      ntv1 = order-1;
      TacsLagrangeShapeFunction<order>(pt[0], ty_knots_order, nu1);
      TacsLagrangeShapeFunction<order-1>(pt[1], ty_knots_reduced, nv1);
    }
    else { // f1 == 2
      ntu1 = order-1;
      ntv1 = order-1;
      TacsLagrangeShapeFunction<order-1>(pt[0], ty_knots_reduced, nu1);
      TacsLagrangeShapeFunction<order-1>(pt[1], ty_knots_reduced, nv1);
    }

    int ntu2, ntv2;
    double nu2[order], nv2[order];
    if (f2 == 0 || f2 == 4){
      ntu2 = order-1;
      ntv2 = order;
      TacsLagrangeShapeFunction<order-1>(pt[0], ty_knots_reduced, nu2);
      TacsLagrangeShapeFunction<order>(pt[1], ty_knots_order, nv2);
    }
    else if (f2 == 1 || f2 == 3){
      ntu2 = order;
      ntv2 = order-1;
      TacsLagrangeShapeFunction<order>(pt[0], ty_knots_order, nu2);
      TacsLagrangeShapeFunction<order-1>(pt[1], ty_knots_reduced, nv2);
    }
    else { // f2 == 2
      ntu2 = order-1;
      ntv2 = order-1;
      TacsLagrangeShapeFunction<order-1>(pt[0], ty_knots_reduced, nu2);
      TacsLagrangeShapeFunction<order-1>(pt[1], ty_knots_reduced, nv2);
    }

    for ( int iv = 0; iv < ntv1; iv++ ){
      for ( int iu = 0; iu < ntu1; iu++ ){
        TacsScalar N = value*nu1[iu]*nv1[iv];

        for ( int jv = 0; jv < ntv2; jv++ ){
          for ( int ju = 0; ju < ntu2; ju++, d2ety++ ){
            d2ety[0] += nu2[ju]*nv2[jv]*N;
          }
        }
      }
    }
  }
};



class TACSShellTriQuadraticBasis {
 public:
  static const int NUM_NODES = 6;

  // Set the number of tying points for each of the 5 components
  // of the tying strain
  static const int NUM_G11_TYING_POINTS = 3;
  static const int NUM_G22_TYING_POINTS = 3;
  static const int NUM_G12_TYING_POINTS = 3;
  static const int NUM_G13_TYING_POINTS = 3;
  static const int NUM_G23_TYING_POINTS = 3;

  static const int NUM_TYING_POINTS =
    NUM_G11_TYING_POINTS +
    NUM_G22_TYING_POINTS +
    NUM_G12_TYING_POINTS +
    NUM_G13_TYING_POINTS +
    NUM_G23_TYING_POINTS;

  static void getNodePoint( const int n, double pt[] ){
    if (n == 0){
      pt[0] = 0.0;  pt[1] = 0.0;
    }
    else if (n == 1){
      pt[0] = 1.0;  pt[1] = 0.0;
    }
    else if (n == 2){
      pt[0] = 0.0;  pt[1] = 1.0;
    }
    else if (n == 3){
      pt[0] = 0.5;  pt[1] = 0.0;
    }
    else if (n == 4){
      pt[0] = 0.5;  pt[1] = 0.5;
    }
    else if (n == 5){
      pt[0] = 0.0;  pt[1] = 0.5;
    }
  }
  static ElementLayout getLayoutType(){
    return TACS_TRI_QUADRATIC_ELEMENT;
  }

  static inline void computeBasis( const double pt[], double N[] ){
    N[0] = (1.0 - pt[0] - pt[1])*(1.0 - 2.0*pt[0] - 2.0*pt[1]);
    N[1] = pt[0]*(2.0*pt[0] - 1.0);
    N[2] = pt[1]*(2.0*pt[1] - 1.0);
    N[3] = 4.0*pt[0]*(1.0 - pt[0] - pt[1]);
    N[4] = 4.0*pt[0]*pt[1];
    N[5] = 4.0*pt[1]*(1.0 - pt[0] - pt[1]);
  }

  static inline void computeBasisGradient( const double pt[], double Nxi[] ){
    Nxi[0] = 4.0*pt[0] + 4.0*pt[1] - 3.0;
    Nxi[1] = 4.0*pt[0] + 4.0*pt[1] - 3.0;
    Nxi[2] = 4.0*pt[0] - 1.0;
    Nxi[3] = 0.0;
    Nxi[4] = 0.0;
    Nxi[5] = 4.0*pt[1] - 1.0;
    Nxi[6] = 4.0 - 8.0*pt[0] - 4.0*pt[1];
    Nxi[7] = -4.0*pt[0];
    Nxi[8] = 4.0*pt[1];
    Nxi[9] = 4.0*pt[0];
    Nxi[10] = -4.0*pt[1];
    Nxi[11] = 4.0 - 4.0*pt[0] - 8.0*pt[1];
  }

  template <int vars_per_node, int m>
  static void interpFields( const double pt[],
                            const TacsScalar values[],
                            TacsScalar field[] ){
    double N[6];
    computeBasis(pt, N);

    for ( int k = 0; k < m; k++ ){
      field[k] = 0.0;
    }

    for ( int i = 0; i < NUM_NODES; i++ ){
      for ( int k = 0; k < m; k++ ){
        field[k] += N[i]*values[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void addInterpFieldsTranspose( const double pt[],
                                        const TacsScalar field[],
                                        TacsScalar values[] ){
    double N[6];
    computeBasis(pt, N);

    for ( int i = 0; i < NUM_NODES; i++ ){
      for ( int k = 0; k < m; k++ ){
        values[k] += N[i]*field[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void interpFieldsGrad( const double pt[],
                                const TacsScalar values[],
                                TacsScalar grad[] ){
    double Nxi[12];
    computeBasisGradient(pt, Nxi);

    for ( int k = 0; k < m; k++ ){
      grad[2*k] = 0.0;
      grad[2*k+1] = 0.0;
    }

    for ( int i = 0; i < NUM_NODES; i++ ){
      for ( int k = 0; k < m; k++ ){
        grad[2*k]   += Nxi[2*i]*values[k];
        grad[2*k+1] += Nxi[2*i+1]*values[k];
      }
      values += vars_per_node;
    }
  }

  template <int vars_per_node, int m>
  static void addInterpFieldsGradTranspose( const double pt[],
                                            TacsScalar grad[],
                                            TacsScalar values[] ){
    double Nxi[12];
    computeBasisGradient(pt, Nxi);

    for ( int i = 0; i < NUM_NODES; i++ ){
      for ( int k = 0; k < m; k++ ){
        values[k] += (Nxi[2*i]*grad[2*k] + Nxi[2*i+1]*grad[2*k+1]);
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
  template <int vars_per_node, int m>
  static void addInterpFieldsOuterProduct( const double pt[],
                                           const TacsScalar jac[],
                                           TacsScalar *mat ){
    double N[6];
    computeBasis(pt, N);

    const int nvars = 6*vars_per_node;

    for ( int jx = 0; jx < NUM_NODES; jx++ ){
      const TacsScalar *jac1 = jac;
      for ( int jm = 0; jm < m; jm++, jac1 += m ){
        for ( int ix = 0; ix < NUM_NODES; ix++ ){
          double Ni = N[jx]*N[ix];

          for ( int im = 0; im < m; im++ ){
            mat[im] += Ni*jac1[im];
          }

          mat += vars_per_node;
        }
      }

      mat += (vars_per_node - m)*nvars;
    }
  }

  /**
    Add the outer-product of the shape functions to the matrix

    jac[2*m*(2*ix + jx) + 2*iy + jy] stores the derivative of the
    2*ix + jx term with respect to the 2*ix + jx component.

    @param pt The parametric location of the quadrature point
    @param m The number of field components
    @param jac The 2m x 2m Jacobian matrix of coefficients
    @param vars_per_node The number of variables per node
    @param mat The element matrix
  */
  template <int vars_per_node, int m>
  static void addInterpGradOuterProduct( const double pt[],
                                         const TacsScalar jac[],
                                         TacsScalar *mat ){
    double Nxi[12];
    computeBasisGradient(pt, Nxi);

    const int nvars = 6*vars_per_node;

    for ( int jx = 0; jx < NUM_NODES; jx++ ){
      const TacsScalar *jac1 = jac;
      const TacsScalar *jac2 = &jac[2*m];

      for ( int jm = 0; jm < m; jm++, jac1 += 4*m, jac2 += 4*m ){
        for ( int ix = 0; ix < NUM_NODES; ix++ ){
          double Naa = Nxi[2*jx]*Nxi[2*ix];
          double Nab = Nxi[2*jx]*Nxi[2*ix+1];
          double Nba = Nxi[2*jx+1]*Nxi[2*ix];
          double Nbb = Nxi[2*jx+1]*Nxi[2*ix+1];

          for ( int im = 0; im < m; im++ ){
            mat[im] += (Naa*jac1[2*im] + Nab*jac1[2*im + 1] +
              Nba*jac2[2*im] + Nbb*jac2[2*im + 1]);
          }

          mat += vars_per_node;
        }
      }

      mat += (vars_per_node - m)*nvars;
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
  template <int vars_per_node, int m>
  static void addInterpGradMixedOuterProduct( const double pt[],
                                              const TacsScalar jac[],
                                              const TacsScalar jacT[],
                                              TacsScalar *mat ){
    double N[6];
    computeBasis(pt, N);

    double Nxi[12];
    computeBasisGradient(pt, Nxi);

    const int nvars = 6*vars_per_node;

    if (jac && jacT){
      for ( int jx = 0; jx < NUM_NODES; jx++ ){
        const TacsScalar *jac1 = jac;

        for ( int jm = 0; jm < m; jm++, jac1 += 2*m ){
          for ( int ix = 0; ix < NUM_NODES; ix++ ){
            double Na1 = N[jx]*Nxi[2*ix];
            double Nb1 = N[jx]*Nxi[2*ix+1];
            double Na2 = N[ix]*Nxi[2*jx];
            double Nb2 = N[ix]*Nxi[2*jx+1];

            const TacsScalar *jac2 = &jacT[2*jm];
            for ( int im = 0; im < m; im++, jac2 += 2*m ){
              mat[im] +=
                Na1*jac1[2*im] + Nb1*jac1[2*im+1] +
                Na2*jac2[0] + Nb2*jac2[1];
            }

            mat += vars_per_node;
          }

          mat += (vars_per_node - m)*nvars;
        }
      }
    }
    else if (jac){
      for ( int jx = 0; jx < NUM_NODES; jx++ ){
        const TacsScalar *jac1 = jac;

        for ( int jm = 0; jm < m; jm++, jac1 += 2*m ){
          for ( int ix = 0; ix < NUM_NODES; ix++ ){
            double Na1 = N[jx]*Nxi[2*ix];
            double Nb1 = N[jx]*Nxi[2*ix+1];

            for ( int im = 0; im < m; im++ ){
              mat[im] += Na1*jac1[2*im] + Nb1*jac1[2*im+1];
            }

            mat += vars_per_node;
          }

          mat += (vars_per_node - m)*nvars;
        }
      }
    }
    else if (jacT){
      for ( int jx = 0; jx < NUM_NODES; jx++ ){
        for ( int jm = 0; jm < m; jm++ ){
          for ( int ix = 0; ix < NUM_NODES; ix++ ){
            double Na2 = N[ix]*Nxi[2*jx];
            double Nb2 = N[ix]*Nxi[2*jx+1];

            const TacsScalar *jac2 = &jacT[2*jm];
            for ( int im = 0; im < m; im++, jac2 += 2*m ){
              mat[im] += Na2*jac2[0] + Nb2*jac2[1];
            }

            mat += vars_per_node;
          }

          mat += (vars_per_node - m)*nvars;
        }
      }
    }
  }

  static int getNumTyingFields(){
    return 5;
  }
  static int getNumTyingPoints( const int field ){
    if (field == 0){ return 3; }
    else if (field == 1){ return 3; }
    else if (field == 2){ return 3; }
    else if (field == 3){ return 3; }
    else if (field == 4){ return 3; }
    return 0;
  }
  static void getTyingPoint( const int field,
                             const int ty,
                             double pt[] ){
    const double s = 0.774596669241483;
    const double t = 0.577350269189626;
    const double s0 = 0.5 - 0.5*s;
    const double t0 = 0.5 - 0.5*t;
    const double t1 = 0.5 + 0.5*t;

    if (field == 0 || field == 4){ // g11 or g13
      if (ty == 0){
        pt[0] = t0;
        pt[1] = 0.0;
      }
      else if (ty == 1){
        pt[0] = t1;
        pt[1] = 0.0;
      }
      else if (ty == 2){
        pt[0] = t0;
        pt[1] = t;
      }
    }
    else if (field == 1 || field == 3){ // g22 or g23
      if (ty == 0){
        pt[0] = 0.0;
        pt[1] = t0;
      }
      else if (ty == 1){
        pt[0] = t;
        pt[1] = t0;
      }
      else if (ty == 2){
        pt[0] = 0.0;
        pt[1] = t1;
      }
    }
    else { // (field == 2) g12
      if (ty == 0){
        pt[0] = t0;
        pt[1] = t0;
      }
      else if (ty == 1){
        pt[0] = t1;
        pt[1] = t0;
      }
      else if (ty == 2){
        pt[0] = t0;
        pt[1] = t1;
      }
    }
  }
  static TacsScalar interpTying( const int field,
                                 const double pt[],
                                 const TacsScalar ety[] ){
    const double s = 0.774596669241483;
    const double t = 0.577350269189626;
    const double s0 = 0.5 - 0.5*s;
    const double t0 = 0.5 - 0.5*t;
    const double tinv = 1.0/t;

    double N[3];
    if (field == 0 || field == 4){
      N[0] = 1.0 - tinv*((pt[0] - t0) + pt[1]);
      N[1] = tinv*(pt[0] - t0);
      N[2] = tinv*pt[1];
    }
    else if (field == 1 || field == 3){
      N[0] = 1.0 - tinv*(pt[0] + (pt[1] - t0));
      N[1] = tinv*pt[0];
      N[2] = tinv*(pt[1] - t0);
    }
    else { // field == 2
      N[0] = 1.0 - tinv*((pt[0] - t0) + (pt[1] - t0));
      N[1] = tinv*(pt[0] - t0);
      N[2] = tinv*(pt[1] - t0);
    }

    TacsScalar value = 0.0;
    for ( int i = 0; i < 3; i++ ){
      value += N[i]*ety[i];
    }

    return value;
  }

  static void addInterpTyingTranspose( const int field,
                                       const double pt[],
                                       const TacsScalar value,
                                       TacsScalar ety[] ){
    const double s = 0.774596669241483;
    const double t = 0.577350269189626;
    const double s0 = 0.5 - 0.5*s;
    const double t0 = 0.5 - 0.5*t;
    const double tinv = 1.0/t;

    double N[3];
    if (field == 0 || field == 4){
      N[0] = 1.0 - tinv*((pt[0] - t0) + pt[1]);
      N[1] = tinv*(pt[0] - t0);
      N[2] = tinv*pt[1];
    }
    else if (field == 1 || field == 3){
      N[0] = 1.0 - tinv*(pt[0] + (pt[1] - t0));
      N[1] = tinv*pt[0];
      N[2] = tinv*(pt[1] - t0);
    }
    else { // field == 2
      N[0] = 1.0 - tinv*((pt[0] - t0) + (pt[1] - t0));
      N[1] = tinv*(pt[0] - t0);
      N[2] = tinv*(pt[1] - t0);
    }

    for ( int i = 0; i < 3; i++ ){
      ety[i] += N[i]*value;
    }
  }

  static void addInterpTyingOuterProduct( const int f1,
                                          const int f2,
                                          const double pt[],
                                          const TacsScalar value,
                                          TacsScalar d2ety[] ){
    const double s = 0.774596669241483;
    const double t = 0.577350269189626;
    const double s0 = 0.5 - 0.5*s;
    const double t0 = 0.5 - 0.5*t;
    const double tinv = 1.0/t;

    double N1[3];
    if (f1 == 0 || f1 == 4){
      N1[0] = 1.0 - tinv*((pt[0] - t0) + pt[1]);
      N1[1] = tinv*(pt[0] - t0);
      N1[2] = tinv*pt[1];
    }
    else if (f1 == 1 || f1 == 3){
      N1[0] = 1.0 - tinv*(pt[0] + (pt[1] - t0));
      N1[1] = tinv*pt[0];
      N1[2] = tinv*(pt[1] - t0);
    }
    else { // f1 == 2
      N1[0] = 1.0 - tinv*((pt[0] - t0) + (pt[1] - t0));
      N1[1] = tinv*(pt[0] - t0);
      N1[2] = tinv*(pt[1] - t0);
    }

    double N2[3];
    if (f2 == 0 || f2 == 4){
      N2[0] = 1.0 - tinv*((pt[0] - t0) + pt[1]);
      N2[1] = tinv*(pt[0] - t0);
      N2[2] = tinv*pt[1];
    }
    else if (f2 == 1 || f2 == 3){
      N2[0] = 1.0 - tinv*(pt[0] + (pt[1] - t0));
      N2[1] = tinv*pt[0];
      N2[2] = tinv*(pt[1] - t0);
    }
    else { // f2 == 2
      N2[0] = 1.0 - tinv*((pt[0] - t0) + (pt[1] - t0));
      N2[1] = tinv*(pt[0] - t0);
      N2[2] = tinv*(pt[1] - t0);
    }

    for ( int i = 0; i < 3; i++ ){
      for ( int j = 0; j < 3; j++, d2ety++ ){
        d2ety[0] += value*N1[i]*N2[j];
      }
    }
  }
};

#endif // TACS_SHELL_ELEMENT_BASIS_H
