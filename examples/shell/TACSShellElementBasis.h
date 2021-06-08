#ifndef TACS_SHELL_ELEMENT_BASIS_H
#define TACS_SHELL_ELEMENT_BASIS_H

#include "TACSGaussQuadrature.h"
#include "TACSTriangleQuadrature.h"
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

class TACSShellQuadLinearBasis {
 public:
  static const int NUM_NODES = 4;

  // Set the number of tying points for each of the 5 components
  // of the tying strain
  static const int NUM_G11_TYING_POINTS = 2;
  static const int NUM_G22_TYING_POINTS = 2;
  static const int NUM_G12_TYING_POINTS = 1;
  static const int NUM_G13_TYING_POINTS = 2;
  static const int NUM_G23_TYING_POINTS = 2;

  static const int NUM_TYING_POINTS =
    NUM_G11_TYING_POINTS +
    NUM_G22_TYING_POINTS +
    NUM_G12_TYING_POINTS +
    NUM_G13_TYING_POINTS +
    NUM_G23_TYING_POINTS;

  static void getNodePoint( const int n, double pt[] ){
    pt[0] = -1.0 + 2.0*(n % 2);
    pt[1] = -1.0 + 2.0*(n / 2);
  }
  static ElementLayout getLayoutType(){
    return TACS_QUAD_ELEMENT;
  }

  template <int vars_per_node, int m>
  static void interpFields( const double pt[],
                            const TacsScalar values[],
                            TacsScalar field[] ){
    double na[2];
    na[0] = 0.5*(1.0 - pt[0]);
    na[1] = 0.5*(1.0 + pt[0]);

    double nb[2];
    nb[0] = 0.5*(1.0 - pt[1]);
    nb[1] = 0.5*(1.0 + pt[1]);

    for ( int k = 0; k < m; k++ ){
      field[k] = 0.0;
    }

    for ( int j = 0; j < 2; j++ ){
      for ( int i = 0; i < 2; i++ ){
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
    double na[2];
    na[0] = 0.5*(1.0 - pt[0]);
    na[1] = 0.5*(1.0 + pt[0]);

    double nb[2];
    nb[0] = 0.5*(1.0 - pt[1]);
    nb[1] = 0.5*(1.0 + pt[1]);

    for ( int j = 0; j < 2; j++ ){
      for ( int i = 0; i < 2; i++ ){
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
    double na[2];
    na[0] = 0.5*(1.0 - pt[0]);
    na[1] = 0.5*(1.0 + pt[0]);

    double nb[2];
    nb[0] = 0.5*(1.0 - pt[1]);
    nb[1] = 0.5*(1.0 + pt[1]);

    double dna[2];
    dna[0] = -0.5;
    dna[1] = 0.5;

    double dnb[2];
    dnb[0] = -0.5;
    dnb[1] = 0.5;

    for ( int k = 0; k < m; k++ ){
      grad[2*k] = 0.0;
      grad[2*k+1] = 0.0;
    }

    for ( int j = 0; j < 2; j++ ){
      for ( int i = 0; i < 2; i++ ){
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
    double na[2];
    na[0] = 0.5*(1.0 - pt[0]);
    na[1] = 0.5*(1.0 + pt[0]);

    double nb[2];
    nb[0] = 0.5*(1.0 - pt[1]);
    nb[1] = 0.5*(1.0 + pt[1]);

    double dna[2];
    dna[0] = -0.5;
    dna[1] = 0.5;

    double dnb[2];
    dnb[0] = -0.5;
    dnb[1] = 0.5;

    for ( int j = 0; j < 2; j++ ){
      for ( int i = 0; i < 2; i++ ){
        for ( int k = 0; k < m; k++ ){
          values[k] += (dna[i]*nb[j]*grad[2*k] + na[i]*dnb[j]*grad[2*k+1]);
        }
        values += vars_per_node;
      }
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
    double na[2];
    na[0] = 0.5*(1.0 - pt[0]);
    na[1] = 0.5*(1.0 + pt[0]);

    double nb[2];
    nb[0] = 0.5*(1.0 - pt[1]);
    nb[1] = 0.5*(1.0 + pt[1]);

    const int nvars = 4*vars_per_node;

    for ( int jy = 0; jy < 2; jy++ ){
      for ( int jx = 0; jx < 2; jx++ ){
        double Nj = na[jx]*nb[jy];

        const TacsScalar *jac1 = jac;
        for ( int jm = 0; jm < m; jm++, jac1 += m ){
          for ( int iy = 0; iy < 2; iy++ ){
            for ( int ix = 0; ix < 2; ix++ ){
              double Ni = na[ix]*nb[iy];

              for ( int im = 0; im < m; im++ ){
                mat[im] += Ni*Nj*jac1[im];
              }

              mat += vars_per_node;
            }
          }
        }

        mat += (vars_per_node - m)*nvars;
      }
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
    double na[2];
    na[0] = 0.5*(1.0 - pt[0]);
    na[1] = 0.5*(1.0 + pt[0]);

    double nb[2];
    nb[0] = 0.5*(1.0 - pt[1]);
    nb[1] = 0.5*(1.0 + pt[1]);

    double dna[2];
    dna[0] = -0.5;
    dna[1] = 0.5;

    double dnb[2];
    dnb[0] = -0.5;
    dnb[1] = 0.5;

    const int nvars = 4*vars_per_node;

    for ( int jy = 0; jy < 2; jy++ ){
      for ( int jx = 0; jx < 2; jx++ ){
        double Naj = dna[jx]*nb[jy];
        double Nbj = na[jx]*dnb[jy];

        const TacsScalar *jac1 = jac;
        const TacsScalar *jac2 = &jac[2*m];
        for ( int jm = 0; jm < m; jm++, jac1 += 4*m, jac2 += 4*m ){
          for ( int iy = 0; iy < 2; iy++ ){
            for ( int ix = 0; ix < 2; ix++ ){
              double Nai = dna[ix]*nb[iy];
              double Nbi = na[ix]*dnb[iy];
              double Naa = Naj*Nai;
              double Nab = Naj*Nbi;
              double Nba = Nbj*Nai;
              double Nbb = Nbj*Nbi;

              for ( int im = 0; im < m; im++ ){
                mat[im] += (Naa*jac1[2*im] + Nab*jac1[2*im + 1] +
                  Nba*jac2[2*im] + Nbb*jac2[2*im + 1]);
              }

              mat += vars_per_node;
            }
          }
        }

        mat += (vars_per_node - m)*nvars;
      }
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
    double na[2];
    na[0] = 0.5*(1.0 - pt[0]);
    na[1] = 0.5*(1.0 + pt[0]);

    double nb[2];
    nb[0] = 0.5*(1.0 - pt[1]);
    nb[1] = 0.5*(1.0 + pt[1]);

    double dna[2];
    dna[0] = -0.5;
    dna[1] = 0.5;

    double dnb[2];
    dnb[0] = -0.5;
    dnb[1] = 0.5;

    const int nvars = 4*vars_per_node;

    if (jac && jacT){
      for ( int jy = 0; jy < 2; jy++ ){
        for ( int jx = 0; jx < 2; jx++ ){
          double Nj = na[jx]*nb[jy];
          double Naj = dna[jx]*nb[jy];
          double Nbj = na[jx]*dnb[jy];

          const TacsScalar *jac1 = jac;
          for ( int jm = 0; jm < m; jm++, jac1 += 2*m ){
            for ( int iy = 0; iy < 2; iy++ ){
              for ( int ix = 0; ix < 2; ix++ ){
                double Ni = na[ix]*nb[iy];
                double Nai = dna[ix]*nb[iy];
                double Nbi = na[ix]*dnb[iy];

                double Na1 = Nj*Nai;
                double Nb1 = Nj*Nbi;
                double Na2 = Ni*Naj;
                double Nb2 = Ni*Nbj;

                const TacsScalar *jac2 = &jacT[2*jm];
                for ( int im = 0; im < m; im++, jac2 += 2*m ){
                  mat[im] +=
                    Na1*jac1[2*im] + Nb1*jac1[2*im+1] +
                    Na2*jac2[0] + Nb2*jac2[1];
                }

                mat += vars_per_node;
              }
            }
          }

          mat += (vars_per_node - m)*nvars;
        }
      }
    }
    else if (jac){
      for ( int jy = 0; jy < 2; jy++ ){
        for ( int jx = 0; jx < 2; jx++ ){
          double Nj = na[jx]*nb[jy];

          const TacsScalar *jac1 = jac;
          for ( int jm = 0; jm < m; jm++, jac1 += 2*m ){
            for ( int iy = 0; iy < 2; iy++ ){
              for ( int ix = 0; ix < 2; ix++ ){
                double Nai = dna[ix]*nb[iy];
                double Nbi = na[ix]*dnb[iy];

                double Na1 = Nj*Nai;
                double Nb1 = Nj*Nbi;

                for ( int im = 0; im < m; im++ ){
                  mat[im] += Na1*jac1[2*im] + Nb1*jac1[2*im+1];
                }

                mat += vars_per_node;
              }
            }
          }

          mat += (vars_per_node - m)*nvars;
        }
      }
    }
    else if (jacT){
      for ( int jy = 0; jy < 2; jy++ ){
        for ( int jx = 0; jx < 2; jx++ ){
          double Naj = dna[jx]*nb[jy];
          double Nbj = na[jx]*dnb[jy];

          for ( int jm = 0; jm < m; jm++ ){
            for ( int iy = 0; iy < 2; iy++ ){
              for ( int ix = 0; ix < 2; ix++ ){
                double Ni = na[ix]*nb[iy];
                double Na2 = Ni*Naj;
                double Nb2 = Ni*Nbj;

                const TacsScalar *jac2 = &jacT[2*jm];
                for ( int im = 0; im < m; im++, jac2 += 2*m ){
                  mat[im] += Na2*jac2[0] + Nb2*jac2[1];
                }

                mat += vars_per_node;
              }
            }
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
    if (field == 0){ return 2; }
    else if (field == 1){ return 2; }
    else if (field == 2){ return 1; }
    else if (field == 3){ return 2; }
    else if (field == 4){ return 2; }
    return 0;
  }
  static void getTyingPoint( const int field,
                             const int ty,
                             double pt[] ){
    if (field == 0 || field == 4){ // g11 or g13
      if (ty == 0){
        pt[0] = 0.0;
        pt[1] = -1.0;
      }
      else {
        pt[0] = 0.0;
        pt[1] = 1.0;
      }
    }
    else if (field == 1 || field == 3){ // g22 or g23
      if (ty == 0){
        pt[0] = -1.0;
        pt[1] = 0.0;
      }
      else {
        pt[0] = 1.0;
        pt[1] = 0.0;
      }
    }
    else { // (field == 2) g12
      pt[0] = pt[1] = 0.0;
    }
  }
  static TacsScalar interpTying( const int field,
                                 const double pt[],
                                 const TacsScalar ety[] ){
    TacsScalar value = 0.0;
    if (field == 0 || field == 4){
      value = 0.5*(1.0 - pt[1])*ety[0] + 0.5*(1.0 + pt[1])*ety[1];
    }
    else if (field == 1 || field == 3){
      value = 0.5*(1.0 - pt[0])*ety[0] + 0.5*(1.0 + pt[0])*ety[1];
    }
    else { // field == 2
      value = ety[0];
    }

    return value;
  }

  static void addInterpTyingTranspose( const int field,
                                       const double pt[],
                                       const TacsScalar value,
                                       TacsScalar ety[] ){
    if (field == 0 || field == 4){
      ety[0] += value*0.5*(1.0 - pt[1]);
      ety[1] += value*0.5*(1.0 + pt[1]);
    }
    else if (field == 1 || field == 3){
      ety[0] += value*0.5*(1.0 - pt[0]);
      ety[1] += value*0.5*(1.0 + pt[0]);
    }
    else { // field == 2
      ety[0] += value;
    }
  }

  static void addInterpTyingOuterProduct( const int f1,
                                          const int f2,
                                          const double pt[],
                                          const TacsScalar value,
                                          TacsScalar d2ety[] ){
    int nf1 = 2, nf2 = 2;
    TacsScalar n1[2];
    double n2[2];
    if (f1 == 0 || f1 == 4){
      n1[0] = value*0.5*(1.0 - pt[1]);
      n1[1] = value*0.5*(1.0 + pt[1]);
    }
    else if (f1 == 1 || f1 == 3){
      n1[0] = value*0.5*(1.0 - pt[0]);
      n1[1] = value*0.5*(1.0 + pt[0]);
    }
    else { // f1 == 2
      nf1 = 1;
      n1[0] = value;
    }

    if (f2 == 0 || f2 == 4){
      n2[0] = 0.5*(1.0 - pt[1]);
      n2[1] = 0.5*(1.0 + pt[1]);
    }
    else if (f2 == 1 || f2 == 3){
      n2[0] = 0.5*(1.0 - pt[0]);
      n2[1] = 0.5*(1.0 + pt[0]);
    }
    else { // f2 == 2
      nf2 = 1;
      n2[0] = 1.0;
    }

    for ( int i = 0; i < nf1; i++ ){
      for ( int j = 0; j < nf2; j++, d2ety++ ){
        d2ety[0] += n1[i]*n2[j];
      }
    }
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

class TACSShellQuadQuadraticBasis {
 public:
  static const int NUM_NODES = 9;

  // Set the number of tying points for each of the 5 components
  // of the tying strain
  static const int NUM_G11_TYING_POINTS = 6;
  static const int NUM_G22_TYING_POINTS = 6;
  static const int NUM_G12_TYING_POINTS = 4;
  static const int NUM_G13_TYING_POINTS = 6;
  static const int NUM_G23_TYING_POINTS = 6;

  static const int NUM_TYING_POINTS =
    NUM_G11_TYING_POINTS +
    NUM_G22_TYING_POINTS +
    NUM_G12_TYING_POINTS +
    NUM_G13_TYING_POINTS +
    NUM_G23_TYING_POINTS;

  static void getNodePoint( const int n, double pt[] ){
    pt[0] = -1.0 + 1.0*(n % 3);
    pt[1] = -1.0 + 1.0*(n / 3);
  }
  static ElementLayout getLayoutType(){
    return TACS_QUAD_QUADRATIC_ELEMENT;
  }

  template <int vars_per_node, int m>
  static void interpFields( const double pt[],
                            const TacsScalar values[],
                            TacsScalar field[] ){
    double na[3];
    na[0] = -0.5*pt[0]*(1.0 - pt[0]);
    na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
    na[2] = 0.5*(1.0 + pt[0])*pt[0];

    double nb[3];
    nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
    nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
    nb[2] = 0.5*(1.0 + pt[1])*pt[1];

    for ( int k = 0; k < m; k++ ){
      field[k] = 0.0;
    }

    for ( int j = 0; j < 3; j++ ){
      for ( int i = 0; i < 3; i++ ){
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
    double na[3];
    na[0] = -0.5*pt[0]*(1.0 - pt[0]);
    na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
    na[2] = 0.5*(1.0 + pt[0])*pt[0];

    double nb[3];
    nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
    nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
    nb[2] = 0.5*(1.0 + pt[1])*pt[1];

    for ( int j = 0; j < 3; j++ ){
      for ( int i = 0; i < 3; i++ ){
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
    double na[3];
    na[0] = -0.5*pt[0]*(1.0 - pt[0]);
    na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
    na[2] = 0.5*(1.0 + pt[0])*pt[0];

    double nb[3];
    nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
    nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
    nb[2] = 0.5*(1.0 + pt[1])*pt[1];

    double dna[3];
    dna[0] = -0.5 + pt[0];
    dna[1] = -2.0*pt[0];
    dna[2] = 0.5 + pt[0];

    double dnb[3];
    dnb[0] = -0.5 + pt[1];
    dnb[1] = -2.0*pt[1];
    dnb[2] = 0.5 + pt[1];

    for ( int k = 0; k < m; k++ ){
      grad[2*k] = 0.0;
      grad[2*k+1] = 0.0;
    }

    for ( int j = 0; j < 3; j++ ){
      for ( int i = 0; i < 3; i++ ){
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
    double na[3];
    na[0] = -0.5*pt[0]*(1.0 - pt[0]);
    na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
    na[2] = 0.5*(1.0 + pt[0])*pt[0];

    double nb[3];
    nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
    nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
    nb[2] = 0.5*(1.0 + pt[1])*pt[1];

    double dna[3];
    dna[0] = -0.5 + pt[0];
    dna[1] = -2.0*pt[0];
    dna[2] = 0.5 + pt[0];

    double dnb[3];
    dnb[0] = -0.5 + pt[1];
    dnb[1] = -2.0*pt[1];
    dnb[2] = 0.5 + pt[1];

    for ( int j = 0; j < 3; j++ ){
      for ( int i = 0; i < 3; i++ ){
        for ( int k = 0; k < m; k++ ){
          values[k] += (dna[i]*nb[j]*grad[2*k] + na[i]*dnb[j]*grad[2*k+1]);
        }
        values += vars_per_node;
      }
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
    double na[3];
    na[0] = -0.5*pt[0]*(1.0 - pt[0]);
    na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
    na[2] = 0.5*(1.0 + pt[0])*pt[0];

    double nb[3];
    nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
    nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
    nb[2] = 0.5*(1.0 + pt[1])*pt[1];

    const int nvars = 9*vars_per_node;

    for ( int jy = 0; jy < 3; jy++ ){
      for ( int jx = 0; jx < 3; jx++ ){
        double Nj = na[jx]*nb[jy];

        const TacsScalar *jac1 = jac;
        for ( int jm = 0; jm < m; jm++, jac1 += m ){
          for ( int iy = 0; iy < 3; iy++ ){
            for ( int ix = 0; ix < 3; ix++ ){
              double Ni = na[ix]*nb[iy];

              for ( int im = 0; im < m; im++ ){
                mat[im] += Ni*Nj*jac1[im];
              }

              mat += vars_per_node;
            }
          }
        }

        mat += (vars_per_node - m)*nvars;
      }
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
    double na[3];
    na[0] = -0.5*pt[0]*(1.0 - pt[0]);
    na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
    na[2] = 0.5*(1.0 + pt[0])*pt[0];

    double nb[3];
    nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
    nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
    nb[2] = 0.5*(1.0 + pt[1])*pt[1];

    double dna[3];
    dna[0] = -0.5 + pt[0];
    dna[1] = -2.0*pt[0];
    dna[2] = 0.5 + pt[0];

    double dnb[3];
    dnb[0] = -0.5 + pt[1];
    dnb[1] = -2.0*pt[1];
    dnb[2] = 0.5 + pt[1];

    const int nvars = 9*vars_per_node;

    for ( int jy = 0; jy < 3; jy++ ){
      for ( int jx = 0; jx < 3; jx++ ){
        double Naj = dna[jx]*nb[jy];
        double Nbj = na[jx]*dnb[jy];

        const TacsScalar *jac1 = jac;
        const TacsScalar *jac2 = &jac[2*m];
        for ( int jm = 0; jm < m; jm++, jac1 += 4*m, jac2 += 4*m ){
          for ( int iy = 0; iy < 3; iy++ ){
            for ( int ix = 0; ix < 3; ix++ ){
              double Nai = dna[ix]*nb[iy];
              double Nbi = na[ix]*dnb[iy];
              double Naa = Naj*Nai;
              double Nab = Naj*Nbi;
              double Nba = Nbj*Nai;
              double Nbb = Nbj*Nbi;

              for ( int im = 0; im < m; im++ ){
                mat[im] += (Naa*jac1[2*im] + Nab*jac1[2*im + 1] +
                  Nba*jac2[2*im] + Nbb*jac2[2*im + 1]);
              }

              mat += vars_per_node;
            }
          }
        }

        mat += (vars_per_node - m)*nvars;
      }
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
    double na[3];
    na[0] = -0.5*pt[0]*(1.0 - pt[0]);
    na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
    na[2] = 0.5*(1.0 + pt[0])*pt[0];

    double nb[3];
    nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
    nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
    nb[2] = 0.5*(1.0 + pt[1])*pt[1];

    double dna[3];
    dna[0] = -0.5 + pt[0];
    dna[1] = -2.0*pt[0];
    dna[2] = 0.5 + pt[0];

    double dnb[3];
    dnb[0] = -0.5 + pt[1];
    dnb[1] = -2.0*pt[1];
    dnb[2] = 0.5 + pt[1];

    const int nvars = 9*vars_per_node;

    if (jac && jacT){
      for ( int jy = 0; jy < 3; jy++ ){
        for ( int jx = 0; jx < 3; jx++ ){
          double Nj = na[jx]*nb[jy];
          double Naj = dna[jx]*nb[jy];
          double Nbj = na[jx]*dnb[jy];

          const TacsScalar *jac1 = jac;
          for ( int jm = 0; jm < m; jm++, jac1 += 2*m ){
            for ( int iy = 0; iy < 3; iy++ ){
              for ( int ix = 0; ix < 3; ix++ ){
                double Ni = na[ix]*nb[iy];
                double Nai = dna[ix]*nb[iy];
                double Nbi = na[ix]*dnb[iy];

                double Na1 = Nj*Nai;
                double Nb1 = Nj*Nbi;
                double Na2 = Ni*Naj;
                double Nb2 = Ni*Nbj;

                const TacsScalar *jac2 = &jacT[2*jm];
                for ( int im = 0; im < m; im++, jac2 += 2*m ){
                  mat[im] +=
                    Na1*jac1[2*im] + Nb1*jac1[2*im+1] +
                    Na2*jac2[0] + Nb2*jac2[1];
                }

                mat += vars_per_node;
              }
            }
          }

          mat += (vars_per_node - m)*nvars;
        }
      }
    }
    else if (jac){
      for ( int jy = 0; jy < 3; jy++ ){
        for ( int jx = 0; jx < 3; jx++ ){
          double Nj = na[jx]*nb[jy];

          const TacsScalar *jac1 = jac;
          for ( int jm = 0; jm < m; jm++, jac1 += 2*m ){
            for ( int iy = 0; iy < 3; iy++ ){
              for ( int ix = 0; ix < 3; ix++ ){
                double Nai = dna[ix]*nb[iy];
                double Nbi = na[ix]*dnb[iy];

                double Na1 = Nj*Nai;
                double Nb1 = Nj*Nbi;

                for ( int im = 0; im < m; im++ ){
                  mat[im] += Na1*jac1[2*im] + Nb1*jac1[2*im+1];
                }

                mat += vars_per_node;
              }
            }
          }

          mat += (vars_per_node - m)*nvars;
        }
      }
    }
    else if (jacT){
      for ( int jy = 0; jy < 3; jy++ ){
        for ( int jx = 0; jx < 3; jx++ ){
          double Naj = dna[jx]*nb[jy];
          double Nbj = na[jx]*dnb[jy];

          for ( int jm = 0; jm < m; jm++ ){
            for ( int iy = 0; iy < 3; iy++ ){
              for ( int ix = 0; ix < 3; ix++ ){
                double Ni = na[ix]*nb[iy];
                double Na2 = Ni*Naj;
                double Nb2 = Ni*Nbj;

                const TacsScalar *jac2 = &jacT[2*jm];
                for ( int im = 0; im < m; im++, jac2 += 2*m ){
                  mat[im] += Na2*jac2[0] + Nb2*jac2[1];
                }

                mat += vars_per_node;
              }
            }
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
    if (field == 0){ return 6; }
    else if (field == 1){ return 6; }
    else if (field == 2){ return 4; }
    else if (field == 3){ return 6; }
    else if (field == 4){ return 6; }
    return 0;
  }
  static void getTyingPoint( const int field,
                             const int ty,
                             double pt[] ){
    const double s = 0.774596669241483;
    const double t = 0.577350269189626;

    if (field == 0 || field == 4){ // g11 or g13
      if (ty % 2 == 0){
        pt[0] = -t;
      }
      else {
        pt[0] = t;
      }
      if (ty / 2 == 0){
        pt[1] = -s;
      }
      else if (ty / 2 == 1){
        pt[1] = 0.0;
      }
      else {
        pt[1] = s;
      }
    }
    else if (field == 1 || field == 3){ // g22 or g23
      if (ty % 3 == 0){
        pt[0] = -s;
      }
      else if (ty % 3 == 1){
        pt[0] = 0.0;
      }
      else {
        pt[0] = s;
      }
      if (ty / 3 == 0){
        pt[1] = -t;
      }
      else {
        pt[1] = t;
      }
    }
    else { // (field == 2) g12
      if (ty % 2 == 0){
        pt[0] = -t;
      }
      else {
        pt[0] = t;
      }
      if (ty / 2 == 0){
        pt[1] = -t;
      }
      else {
        pt[1] = t;
      }
    }
  }
  static TacsScalar interpTying( const int field,
                                 const double pt[],
                                 const TacsScalar ety[] ){
    const double s = 0.774596669241483;
    const double t = 0.577350269189626;
    const double tinv = 1.0/t;
    const double sinv = 1.0/(s*s);

    TacsScalar value = 0.0;
    if (field == 0 || field == 4){
      double ntu[2];
      ntu[0] = 0.5*tinv*(t - pt[0]);
      ntu[1] = 0.5*tinv*(t + pt[0]);

      double nv[3];
      nv[0] = 0.5*sinv*pt[1]*(pt[1] - s);
      nv[1] = sinv*(s - pt[1])*(s + pt[1]);
      nv[2] = 0.5*sinv*pt[1]*(s + pt[1]);

      for ( int j = 0; j < 3; j++ ){
        for ( int i = 0; i < 2; i++ ){
          value += ntu[i]*nv[j]*ety[i + 2*j];
        }
      }
    }
    else if (field == 1 || field == 3){
      double nu[3];
      nu[0] = 0.5*sinv*pt[0]*(pt[0] - s);
      nu[1] = sinv*(s - pt[0])*(s + pt[0]);
      nu[2] = 0.5*sinv*pt[0]*(s + pt[0]);

      double ntv[2];
      ntv[0] = 0.5*tinv*(t - pt[1]);
      ntv[1] = 0.5*tinv*(t + pt[1]);

      for ( int j = 0; j < 2; j++ ){
        for ( int i = 0; i < 3; i++ ){
          value += nu[i]*ntv[j]*ety[i + 3*j];
        }
      }
    }
    else { // field == 2
      double ntu[2];
      ntu[0] = 0.5*tinv*(t - pt[0]);
      ntu[1] = 0.5*tinv*(t + pt[0]);

      double ntv[2];
      ntv[0] = 0.5*tinv*(t - pt[1]);
      ntv[1] = 0.5*tinv*(t + pt[1]);

      for ( int j = 0; j < 2; j++ ){
        for ( int i = 0; i < 2; i++ ){
          value += ntu[i]*ntv[j]*ety[i + 2*j];
        }
      }
    }

    return value;
  }

  static void addInterpTyingTranspose( const int field,
                                       const double pt[],
                                       const TacsScalar value,
                                       TacsScalar ety[] ){
    const double s = 0.774596669241483;
    const double t = 0.577350269189626;
    const double tinv = 1.0/t;
    const double sinv = 1.0/(s*s);

    if (field == 0 || field == 4){
      double ntu[2];
      ntu[0] = 0.5*tinv*(t - pt[0]);
      ntu[1] = 0.5*tinv*(t + pt[0]);

      double nv[3];
      nv[0] = 0.5*sinv*pt[1]*(pt[1] - s);
      nv[1] = sinv*(s - pt[1])*(s + pt[1]);
      nv[2] = 0.5*sinv*pt[1]*(s + pt[1]);

      for ( int j = 0; j < 3; j++ ){
        for ( int i = 0; i < 2; i++ ){
          ety[i + 2*j] += ntu[i]*nv[j]*value;
        }
      }
    }
    else if (field == 1 || field == 3){
      double nu[3];
      nu[0] = 0.5*sinv*pt[0]*(pt[0] - s);
      nu[1] = sinv*(s - pt[0])*(s + pt[0]);
      nu[2] = 0.5*sinv*pt[0]*(s + pt[0]);

      double ntv[2];
      ntv[0] = 0.5*tinv*(t - pt[1]);
      ntv[1] = 0.5*tinv*(t + pt[1]);

      for ( int j = 0; j < 2; j++ ){
        for ( int i = 0; i < 3; i++ ){
          ety[i + 3*j] += nu[i]*ntv[j]*value;
        }
      }
    }
    else { // field == 2
      double ntu[2];
      ntu[0] = 0.5*tinv*(t - pt[0]);
      ntu[1] = 0.5*tinv*(t + pt[0]);

      double ntv[2];
      ntv[0] = 0.5*tinv*(t - pt[1]);
      ntv[1] = 0.5*tinv*(t + pt[1]);

      for ( int j = 0; j < 2; j++ ){
        for ( int i = 0; i < 2; i++ ){
          ety[i + 2*j] += ntu[i]*ntv[j]*value;
        }
      }
    }
  }

  static void addInterpTyingOuterProduct( const int f1,
                                          const int f2,
                                          const double pt[],
                                          const TacsScalar value,
                                          TacsScalar d2ety[] ){
    const double s = 0.774596669241483;
    const double t = 0.577350269189626;
    const double tinv = 1.0/t;
    const double sinv = 1.0/(s*s);

    int ntu1, ntv1;
    double nu1[3], nv1[3];

    int ntu2, ntv2;
    double nu2[3], nv2[3];

    if (f1 == 0 || f1 == 4){
      ntu1 = 2;
      nu1[0] = 0.5*tinv*(t - pt[0]);
      nu1[1] = 0.5*tinv*(t + pt[0]);

      ntv1 = 3;
      nv1[0] = 0.5*sinv*pt[1]*(pt[1] - s);
      nv1[1] = sinv*(s - pt[1])*(s + pt[1]);
      nv1[2] = 0.5*sinv*pt[1]*(s + pt[1]);
    }
    else if (f1 == 1 || f1 == 3){
      ntu1 = 3;
      nu1[0] = 0.5*sinv*pt[0]*(pt[0] - s);
      nu1[1] = sinv*(s - pt[0])*(s + pt[0]);
      nu1[2] = 0.5*sinv*pt[0]*(s + pt[0]);

      ntv1 = 2;
      nv1[0] = 0.5*tinv*(t - pt[1]);
      nv1[1] = 0.5*tinv*(t + pt[1]);
    }
    else { // f1 == 2
      ntu1 = 2;
      nu1[0] = 0.5*tinv*(t - pt[0]);
      nu1[1] = 0.5*tinv*(t + pt[0]);

      ntv1 = 2;
      nv1[0] = 0.5*tinv*(t - pt[1]);
      nv1[1] = 0.5*tinv*(t + pt[1]);
    }

    if (f2 == 0 || f2 == 4){
      ntu2 = 2;
      nu2[0] = 0.5*tinv*(t - pt[0]);
      nu2[1] = 0.5*tinv*(t + pt[0]);

      ntv2 = 3;
      nv2[0] = 0.5*sinv*pt[1]*(pt[1] - s);
      nv2[1] = sinv*(s - pt[1])*(s + pt[1]);
      nv2[2] = 0.5*sinv*pt[1]*(s + pt[1]);
    }
    else if (f2 == 1 || f2 == 3){
      ntu2 = 3;
      nu2[0] = 0.5*sinv*pt[0]*(pt[0] - s);
      nu2[1] = sinv*(s - pt[0])*(s + pt[0]);
      nu2[2] = 0.5*sinv*pt[0]*(s + pt[0]);

      ntv2 = 2;
      nv2[0] = 0.5*tinv*(t - pt[1]);
      nv2[1] = 0.5*tinv*(t + pt[1]);
    }
    else { // f2 == 2
      ntu2 = 2;
      nu2[0] = 0.5*tinv*(t - pt[0]);
      nu2[1] = 0.5*tinv*(t + pt[0]);

      ntv2 = 2;
      nv2[0] = 0.5*tinv*(t - pt[1]);
      nv2[1] = 0.5*tinv*(t + pt[1]);
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