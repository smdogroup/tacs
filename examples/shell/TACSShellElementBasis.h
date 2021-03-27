#ifndef TACS_SHELL_ELEMENT_BASIS_H
#define TACS_SHELL_ELEMENT_BASIS_H

#include "TACSGaussQuadrature.h"
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

  static void interpFields( const double pt[],
                            const int vars_per_node,
                            const TacsScalar values[],
                            const int m,
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

  static void addInterpFieldsTranspose( const double pt[],
                                        const int m,
                                        const TacsScalar field[],
                                        const int vars_per_node,
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

  static void interpFieldsGrad( const double pt[],
                                const int vars_per_node,
                                const TacsScalar values[],
                                const int m,
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

  static void addInterpFieldsGradTranspose( const double pt[],
                                            const int m,
                                            TacsScalar grad[],
                                            const int vars_per_node,
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

    mat[row_incr*i + col_incr*j] += scale*N[i]*N[j]

    @param pt The parametric location of the quadrature point
    @param weight The weight factor added to the matrix
    @param row_incr The row increment applied to the matrix
    @param col_incr The column increment applied to the matrix
    @param mat The element matrix
  */
  virtual void addInterpGradGradOuterProduct( const double pt[],
                                              TacsScalar jac[],
                                              const int row_incr,
                                              const int col_incr,
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

    for ( int jy = 0; jy < 2; jy++ ){
      for ( int jx = 0; jx < 2; jx++, mat += row_incr ){
        double Naj = dna[jx]*nb[jy];
        double Nbj = na[jx]*dnb[jy];

        for ( int iy = 0; iy < 2; iy++ ){
          for ( int ix = 0; ix < 2; ix++, mat += col_incr ){
            double Nai = dna[ix]*nb[iy];
            double Nbi = na[ix]*dnb[iy];

            mat[0] +=
              (Naj*Nai*jac[0] + Naj*Nbi*jac[1] +
               Nbj*Nai*jac[2] + Nbj*Nbi*jac[3]);
          }
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

  static void interpFields( const double pt[],
                            const int vars_per_node,
                            const TacsScalar values[],
                            const int m,
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

  static void addInterpFieldsTranspose( const double pt[],
                                        const int m,
                                        const TacsScalar field[],
                                        const int vars_per_node,
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

  static void interpFieldsGrad( const double pt[],
                                const int vars_per_node,
                                const TacsScalar values[],
                                const int m,
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

  static void addInterpFieldsGradTranspose( const double pt[],
                                            const int m,
                                            TacsScalar grad[],
                                            const int vars_per_node,
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

    mat[row_incr*i + col_incr*j] += scale*N[i]*N[j]

    @param pt The parametric location of the quadrature point
    @param weight The weight factor added to the matrix
    @param row_incr The row increment applied to the matrix
    @param col_incr The column increment applied to the matrix
    @param mat The element matrix
  */
  virtual void addInterpGradGradOuterProduct( const double pt[],
                                              TacsScalar jac[],
                                              const int row_incr,
                                              const int col_incr,
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

    for ( int jy = 0; jy < 3; jy++ ){
      for ( int jx = 0; jx < 3; jx++, mat += row_incr ){
        double Naj = dna[jx]*nb[jy];
        double Nbj = na[jx]*dnb[jy];

        for ( int iy = 0; iy < 3; iy++ ){
          for ( int ix = 0; ix < 3; ix++, mat += col_incr ){
            double Nai = dna[ix]*nb[iy];
            double Nbi = na[ix]*dnb[iy];

            mat[0] +=
              (Naj*Nai*jac[0] + Naj*Nbi*jac[1] +
               Nbj*Nai*jac[2] + Nbj*Nbi*jac[3]);
          }
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
};

#endif // TACS_SHELL_ELEMENT_BASIS_H