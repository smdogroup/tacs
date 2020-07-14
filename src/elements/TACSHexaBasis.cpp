/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TACSHexaBasis.h"
#include "TACSGaussQuadrature.h"
#include "TACSLagrangeInterpolation.h"
#include "TACSBasisMacros.h"

static void getFaceTangents( int face, double t[] ){
  if (face == 0){
    // Z - Y plane make a -X normal direction
    t[0] = 0.0;  t[1] = 0.0;  t[2] = 1.0;
    t[3] = 0.0;  t[4] = 1.0;  t[5] = 0.0;
  }
  else if (face == 1){
    // Y - Z plane makes a X normal direction
    t[0] = 0.0;  t[1] = 1.0;  t[2] = 0.0;
    t[3] = 0.0;  t[4] = 0.0;  t[5] = 1.0;
  }
  else if (face == 2){
    // X - Z plane makes a -Y normal direction
    t[0] = 1.0;  t[1] = 0.0;  t[2] = 0.0;
    t[3] = 0.0;  t[4] = 0.0;  t[5] = 1.0;
  }
  else if (face == 3){
    // Z - X plane makes a Y normal direction
    t[0] = 0.0;  t[1] = 0.0;  t[2] = 1.0;
    t[3] = 1.0;  t[4] = 0.0;  t[5] = 0.0;
  }
  else if (face == 4){
    // Y - X plane makes a -Z normal direction
    t[0] = 0.0;  t[1] = 1.0;  t[2] = 0.0;
    t[3] = 1.0;  t[4] = 0.0;  t[5] = 0.0;
  }
  else if (face == 5){
    // X - Y plane makes a Z normal direction
    t[0] = 1.0;  t[1] = 0.0;  t[2] = 0.0;
    t[3] = 0.0;  t[4] = 0.0;  t[5] = 1.0;
  }
}

/*
  Linear Hexa basis class functions
*/
ElementLayout TACSLinearHexaBasis::getLayoutType(){
  return TACS_HEXA_ELEMENT;
}

void TACSLinearHexaBasis::getVisPoint( int n, double pt[] ){
  pt[0] = -1.0 + 2.0*(n % 2);
  pt[1] = -1.0 + 2.0*((n % 4) / 2);
  pt[2] = -1.0 + 2.0*(n / 4);
}

int TACSLinearHexaBasis::getNumNodes(){
  return 8;
}

int TACSLinearHexaBasis::getNumParameters(){
  return 3;
}

int TACSLinearHexaBasis::getNumQuadraturePoints(){
  return 8;
}

double TACSLinearHexaBasis::getQuadratureWeight( int n ){
  return (TacsGaussQuadWts2[n % 2]*
          TacsGaussQuadWts2[(n % 4)/2]*
          TacsGaussQuadWts2[n/4]);
}

double TACSLinearHexaBasis::getQuadraturePoint( int n,
                                                double pt[] ){
  pt[0] = TacsGaussQuadPts2[n % 2];
  pt[1] = TacsGaussQuadPts2[(n % 4)/2];
  pt[2] = TacsGaussQuadPts2[n/4];

  return (TacsGaussQuadWts2[n % 2]*
          TacsGaussQuadWts2[(n % 4)/2]*
          TacsGaussQuadWts2[n/4]);
}

int TACSLinearHexaBasis::getNumElementFaces(){
  return 6;
}

int TACSLinearHexaBasis::getNumFaceQuadraturePoints( int face ){
  return 4;
}

double TACSLinearHexaBasis::getFaceQuadraturePoint( int face,
                                                    int n,
                                                    double pt[],
                                                    double t[] ){
  if (face/2 == 0){
    pt[0] = -1.0 + 2.0*(face % 2);
    pt[1] = TacsGaussQuadPts2[n % 2];
    pt[2] = TacsGaussQuadPts2[n/2];
  }
  else if (face/2 == 1){
    pt[0] = TacsGaussQuadPts2[n % 2];
    pt[1] = -1.0 + 2.0*(face % 2);
    pt[2] = TacsGaussQuadPts2[n/2];
  }
  else {
    pt[0] = TacsGaussQuadPts2[n % 2];
    pt[1] = TacsGaussQuadPts2[n/2];
    pt[2] = -1.0 + 2.0*(face % 2);
  }

  getFaceTangents(face, t);

  return (TacsGaussQuadWts2[n % 2]*TacsGaussQuadWts2[n/2]);
}

void TACSLinearHexaBasis::computeBasis( const double pt[],
                                        double N[] ){
  N[0] = 0.125*(1.0 - pt[0])*(1.0 - pt[1])*(1.0 - pt[2]);
  N[1] = 0.125*(1.0 + pt[0])*(1.0 - pt[1])*(1.0 - pt[2]);
  N[2] = 0.125*(1.0 - pt[0])*(1.0 + pt[1])*(1.0 - pt[2]);
  N[3] = 0.125*(1.0 + pt[0])*(1.0 + pt[1])*(1.0 - pt[2]);
  N[4] = 0.125*(1.0 - pt[0])*(1.0 - pt[1])*(1.0 + pt[2]);
  N[5] = 0.125*(1.0 + pt[0])*(1.0 - pt[1])*(1.0 + pt[2]);
  N[6] = 0.125*(1.0 - pt[0])*(1.0 + pt[1])*(1.0 + pt[2]);
  N[7] = 0.125*(1.0 + pt[0])*(1.0 + pt[1])*(1.0 + pt[2]);
}

void TACSLinearHexaBasis::computeBasisGradient( const double pt[],
                                                double N[],
                                                double Nxi[] ){
  N[0] = 0.125*(1.0 - pt[0])*(1.0 - pt[1])*(1.0 - pt[2]);
  N[1] = 0.125*(1.0 + pt[0])*(1.0 - pt[1])*(1.0 - pt[2]);
  N[2] = 0.125*(1.0 - pt[0])*(1.0 + pt[1])*(1.0 - pt[2]);
  N[3] = 0.125*(1.0 + pt[0])*(1.0 + pt[1])*(1.0 - pt[2]);
  N[4] = 0.125*(1.0 - pt[0])*(1.0 - pt[1])*(1.0 + pt[2]);
  N[5] = 0.125*(1.0 + pt[0])*(1.0 - pt[1])*(1.0 + pt[2]);
  N[6] = 0.125*(1.0 - pt[0])*(1.0 + pt[1])*(1.0 + pt[2]);
  N[7] = 0.125*(1.0 + pt[0])*(1.0 + pt[1])*(1.0 + pt[2]);

  Nxi[0] = -0.125*(1.0 - pt[1])*(1.0 - pt[2]);
  Nxi[1] = -0.125*(1.0 - pt[0])*(1.0 - pt[2]);
  Nxi[2] = -0.125*(1.0 - pt[0])*(1.0 - pt[1]);

  Nxi[3] = 0.125*(1.0 - pt[1])*(1.0 - pt[2]);
  Nxi[4] = -0.125*(1.0 + pt[0])*(1.0 - pt[2]);
  Nxi[5] = -0.125*(1.0 + pt[0])*(1.0 - pt[1]);

  Nxi[6] = -0.125*(1.0 + pt[1])*(1.0 - pt[2]);
  Nxi[7] = 0.125*(1.0 - pt[0])*(1.0 - pt[2]);
  Nxi[8] = -0.125*(1.0 - pt[0])*(1.0 + pt[1]);

  Nxi[9] = 0.125*(1.0 + pt[1])*(1.0 - pt[2]);
  Nxi[10] = 0.125*(1.0 + pt[0])*(1.0 - pt[2]);
  Nxi[11] = -0.125*(1.0 + pt[0])*(1.0 + pt[1]);

  Nxi[12] = -0.125*(1.0 - pt[1])*(1.0 + pt[2]);
  Nxi[13] = -0.125*(1.0 - pt[0])*(1.0 + pt[2]);
  Nxi[14] = 0.125*(1.0 - pt[0])*(1.0 - pt[1]);

  Nxi[15] = 0.125*(1.0 - pt[1])*(1.0 + pt[2]);
  Nxi[16] = -0.125*(1.0 + pt[0])*(1.0 + pt[2]);
  Nxi[17] = 0.125*(1.0 + pt[0])*(1.0 - pt[1]);

  Nxi[18] = -0.125*(1.0 + pt[1])*(1.0 + pt[2]);
  Nxi[19] = 0.125*(1.0 - pt[0])*(1.0 + pt[2]);
  Nxi[20] = 0.125*(1.0 - pt[0])*(1.0 + pt[1]);

  Nxi[21] = 0.125*(1.0 + pt[1])*(1.0 + pt[2]);
  Nxi[22] = 0.125*(1.0 + pt[0])*(1.0 + pt[2]);
  Nxi[23] = 0.125*(1.0 + pt[0])*(1.0 + pt[1]);
}

/*
  Quadratic Hexa basis class functions
*/
ElementLayout TACSQuadraticHexaBasis::getLayoutType(){
  return TACS_HEXA_QUADRATIC_ELEMENT;
}

void TACSQuadraticHexaBasis::getVisPoint( int n, double pt[] ){
  pt[0] = -1.0 + 2.0*(n % 3);
  pt[1] = -1.0 + 2.0*((n % 9) / 3);
  pt[2] = -1.0 + 2.0*(n / 9);
}

int TACSQuadraticHexaBasis::getNumNodes(){
  return 27;
}

int TACSQuadraticHexaBasis::getNumParameters(){
  return 3;
}

int TACSQuadraticHexaBasis::getNumQuadraturePoints(){
  return 27;
}

double TACSQuadraticHexaBasis::getQuadratureWeight( int n ){
  return (TacsGaussQuadWts3[n % 3]*
          TacsGaussQuadWts3[(n % 9)/3]*
          TacsGaussQuadWts3[n/9]);
}

double TACSQuadraticHexaBasis::getQuadraturePoint( int n,
                                                   double pt[] ){
  pt[0] = TacsGaussQuadPts3[n % 3];
  pt[1] = TacsGaussQuadPts3[(n % 9)/3];
  pt[2] = TacsGaussQuadPts3[n/9];

  return (TacsGaussQuadWts3[n % 3]*
          TacsGaussQuadWts3[(n % 9)/3]*
          TacsGaussQuadWts3[n/9]);
}

int TACSQuadraticHexaBasis::getNumElementFaces(){
  return 6;
}

int TACSQuadraticHexaBasis::getNumFaceQuadraturePoints( int face ){
  return 9;
}

double TACSQuadraticHexaBasis::getFaceQuadraturePoint( int face,
                                                       int n,
                                                       double pt[],
                                                       double t[] ){
  if (face/2 == 0){
    pt[0] = -1.0 + 2.0*(face % 2);
    pt[1] = TacsGaussQuadPts3[n % 3];
    pt[2] = TacsGaussQuadPts3[n/3];
  }
  else if (face/2 == 1){
    pt[0] = TacsGaussQuadPts3[n % 3];
    pt[1] = -1.0 + 2.0*(face % 2);
    pt[2] = TacsGaussQuadPts3[n/3];
  }
  else {
    pt[0] = TacsGaussQuadPts3[n % 3];
    pt[1] = TacsGaussQuadPts3[n/3];
    pt[2] = -1.0 + 2.0*(face % 2);
  }

  getFaceTangents(face, t);

  return (TacsGaussQuadWts3[n % 3]*TacsGaussQuadWts3[n/3]);
}

void TACSQuadraticHexaBasis::computeBasis( const double pt[],
                                           double N[] ){
  double na[3];
  na[0] = -0.5*pt[0]*(1.0 - pt[0]);
  na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
  na[2] = 0.5*(1.0 + pt[0])*pt[0];

  double nb[3];
  nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
  nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
  nb[2] = 0.5*(1.0 + pt[1])*pt[1];

  double nc[3];
  nc[0] = -0.5*pt[2]*(1.0 - pt[2]);
  nc[1] = (1.0 - pt[2])*(1.0 + pt[2]);
  nc[2] = 0.5*(1.0 + pt[2])*pt[2];

  for ( int k = 0; k < 3; k++ ){
    for ( int j = 0; j < 3; j++ ){
      for ( int i = 0; i < 3; i++ ){
        N[0] = na[i]*nb[j]*nc[k];
        N++;
      }
    }
  }
}

void TACSQuadraticHexaBasis::computeBasisGradient( const double pt[],
                                                    double N[],
                                                    double Nxi[] ){
  double na[3];
  na[0] = -0.5*pt[0]*(1.0 - pt[0]);
  na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
  na[2] = 0.5*(1.0 + pt[0])*pt[0];

  double nb[3];
  nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
  nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
  nb[2] = 0.5*(1.0 + pt[1])*pt[1];

  double nc[3];
  nc[0] = -0.5*pt[2]*(1.0 - pt[2]);
  nc[1] = (1.0 - pt[2])*(1.0 + pt[2]);
  nc[2] = 0.5*(1.0 + pt[2])*pt[2];

  double dna[3];
  dna[0] = -0.5 + pt[0];
  dna[1] = -2.0*pt[0];
  dna[2] = 0.5 + pt[0];

  double dnb[3];
  dnb[0] = -0.5 + pt[1];
  dnb[1] = -2.0*pt[1];
  dnb[2] = 0.5 + pt[1];

  double dnc[3];
  dnc[0] = -0.5 + pt[2];
  dnc[1] = -2.0*pt[2];
  dnc[2] = 0.5 + pt[2];

  for ( int k = 0; k < 3; k++ ){
    for ( int j = 0; j < 3; j++ ){
      for ( int i = 0; i < 3; i++ ){
        N[0] = na[i]*nb[j]*nc[k];
        Nxi[0] = dna[i]*nb[j]*nc[k];
        Nxi[1] = na[i]*dnb[j]*nc[k];
        Nxi[2] = na[i]*nb[j]*dnc[k];
        N++;
        Nxi += 3;
      }
    }
  }
}

void TACSQuadraticHexaBasis::interpFields( const int n,
                                           const double pt[],
                                           const int m,
                                           const TacsScalar v[],
                                           const int incr,
                                           TacsScalar u[] ){
  double na[3];
  na[0] = -0.5*pt[0]*(1.0 - pt[0]);
  na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
  na[2] = 0.5*(1.0 + pt[0])*pt[0];

  double nb[3];
  nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
  nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
  nb[2] = 0.5*(1.0 + pt[1])*pt[1];

  double nc[3];
  nc[0] = -0.5*pt[2]*(1.0 - pt[2]);
  nc[1] = (1.0 - pt[2])*(1.0 + pt[2]);
  nc[2] = 0.5*(1.0 + pt[2])*pt[2];

  for ( int i = 0; i < m; i++, u += incr, v++ ){
    u[0] = TACS_BASIS_EVAL_TENSOR3D_ORDER3(na, nb, nc, m, v);
  }
}

void TACSQuadraticHexaBasis::addInterpFieldsTranspose( const int n,
                                                       const double pt[],
                                                       const int incr,
                                                       const TacsScalar u[],
                                                       const int m,
                                                       TacsScalar v[] ){
  double na[3];
  na[0] = -0.5*pt[0]*(1.0 - pt[0]);
  na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
  na[2] = 0.5*(1.0 + pt[0])*pt[0];

  double nb[3];
  nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
  nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
  nb[2] = 0.5*(1.0 + pt[1])*pt[1];

  double nc[3];
  nc[0] = -0.5*pt[2]*(1.0 - pt[2]);
  nc[1] = (1.0 - pt[2])*(1.0 + pt[2]);
  nc[2] = 0.5*(1.0 + pt[2])*pt[2];

  for ( int i = 0; i < m; i++, u += incr, v++ ){
    TacsScalar temp;
    TACS_BASIS_TRANSPOSE_TENSOR3D_ORDER3(na, nb, nc, u[0], temp, m, v);
  }
}

void TACSQuadraticHexaBasis::interpFieldsGrad( const int n,
                                               const double pt[],
                                               const int m,
                                               const TacsScalar v[],
                                               TacsScalar g[] ){
  double na[3];
  na[0] = -0.5*pt[0]*(1.0 - pt[0]);
  na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
  na[2] = 0.5*(1.0 + pt[0])*pt[0];

  double nb[3];
  nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
  nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
  nb[2] = 0.5*(1.0 + pt[1])*pt[1];

  double nc[3];
  nc[0] = -0.5*pt[2]*(1.0 - pt[2]);
  nc[1] = (1.0 - pt[2])*(1.0 + pt[2]);
  nc[2] = 0.5*(1.0 + pt[2])*pt[2];

  double dna[3];
  dna[0] = -0.5 + pt[0];
  dna[1] = -2.0*pt[0];
  dna[2] = 0.5 + pt[0];

  double dnb[3];
  dnb[0] = -0.5 + pt[1];
  dnb[1] = -2.0*pt[1];
  dnb[2] = 0.5 + pt[1];

  double dnc[3];
  dnc[0] = -0.5 + pt[2];
  dnc[1] = -2.0*pt[2];
  dnc[2] = 0.5 + pt[2];

  for ( int i = 0; i < m; i++, g += 3, v++ ){
    g[0] = TACS_BASIS_EVAL_TENSOR3D_ORDER3(dna, nb, nc, m, v);
    g[1] = TACS_BASIS_EVAL_TENSOR3D_ORDER3(na, dnb, nc, m, v);
    g[2] = TACS_BASIS_EVAL_TENSOR3D_ORDER3(na, nb, dnc, m, v);
  }
}

void TACSQuadraticHexaBasis::addInterpFieldsGradTranspose( int n,
                                                           const double pt[],
                                                           const int m,
                                                           const TacsScalar g[],
                                                           TacsScalar v[] ){
  double na[3];
  na[0] = -0.5*pt[0]*(1.0 - pt[0]);
  na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
  na[2] = 0.5*(1.0 + pt[0])*pt[0];

  double nb[3];
  nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
  nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
  nb[2] = 0.5*(1.0 + pt[1])*pt[1];

  double nc[3];
  nc[0] = -0.5*pt[2]*(1.0 - pt[2]);
  nc[1] = (1.0 - pt[2])*(1.0 + pt[2]);
  nc[2] = 0.5*(1.0 + pt[2])*pt[2];

  double dna[3];
  dna[0] = -0.5 + pt[0];
  dna[1] = -2.0*pt[0];
  dna[2] = 0.5 + pt[0];

  double dnb[3];
  dnb[0] = -0.5 + pt[1];
  dnb[1] = -2.0*pt[1];
  dnb[2] = 0.5 + pt[1];

  double dnc[3];
  dnc[0] = -0.5 + pt[2];
  dnc[1] = -2.0*pt[2];
  dnc[2] = 0.5 + pt[2];

  for ( int i = 0; i < m; i++, g += 3, v++ ){
    TacsScalar temp;
    TACS_BASIS_TRANSPOSE_TENSOR3D_ORDER3(dna, nb, nc, g[0], temp, m, v);
    TACS_BASIS_TRANSPOSE_TENSOR3D_ORDER3(na, dnb, nc, g[1], temp, m, v);
    TACS_BASIS_TRANSPOSE_TENSOR3D_ORDER3(na, nb, dnc, g[2], temp, m, v);
  }
}

/*
  Cubic Hexa basis class functions
*/
ElementLayout TACSCubicHexaBasis::getLayoutType(){
  return TACS_HEXA_CUBIC_ELEMENT;
}

void TACSCubicHexaBasis::getVisPoint( int n, double pt[] ){
  pt[0] = -1.0 + 2.0*(n % 4);
  pt[1] = -1.0 + 2.0*((n % 16) / 4);
  pt[2] = -1.0 + 2.0*(n / 16);
}

int TACSCubicHexaBasis::getNumNodes(){
  return 64;
}

int TACSCubicHexaBasis::getNumParameters(){
  return 3;
}

int TACSCubicHexaBasis::getNumQuadraturePoints(){
  return 64;
}

double TACSCubicHexaBasis::getQuadratureWeight( int n ){
  return (TacsGaussQuadWts4[n % 4]*
          TacsGaussQuadWts4[(n % 16)/4]*
          TacsGaussQuadWts4[n/16]);
}

double TACSCubicHexaBasis::getQuadraturePoint( int n,
                                               double pt[] ){
  pt[0] = TacsGaussQuadPts4[n % 4];
  pt[1] = TacsGaussQuadPts4[(n % 16)/4];
  pt[2] = TacsGaussQuadPts4[n/16];

  return (TacsGaussQuadWts4[n % 4]*
          TacsGaussQuadWts4[(n % 16)/4]*
          TacsGaussQuadWts4[n/16]);
}

int TACSCubicHexaBasis::getNumElementFaces(){
  return 6;
}

int TACSCubicHexaBasis::getNumFaceQuadraturePoints( int face ){
  return 16;
}

double TACSCubicHexaBasis::getFaceQuadraturePoint( int face,
                                                   int n,
                                                   double pt[],
                                                   double t[] ){
  if (face/2 == 0){
    pt[0] = -1.0 + 2.0*(face % 2);
    pt[1] = TacsGaussQuadPts4[n % 4];
    pt[2] = TacsGaussQuadPts4[n/4];
  }
  else if (face/2 == 1){
    pt[0] = TacsGaussQuadPts3[n % 4];
    pt[1] = -1.0 + 2.0*(face % 2);
    pt[2] = TacsGaussQuadPts3[n/4];
  }
  else {
    pt[0] = TacsGaussQuadPts3[n % 4];
    pt[1] = TacsGaussQuadPts3[n/4];
    pt[2] = -1.0 + 2.0*(face % 2);
  }

  getFaceTangents(face, t);

  return (TacsGaussQuadWts4[n % 4]*TacsGaussQuadWts4[n/4]);
}

void TACSCubicHexaBasis::computeBasis( const double pt[],
                                       double N[] ){
  double na[4];
  na[0] = -(2.0/3.0)*(0.5 + pt[0])*(0.5 - pt[0])*(1.0 - pt[0]);
  na[1] = (4.0/3.0)*(1.0 + pt[0])*(0.5 - pt[0])*(1.0 - pt[0]);
  na[2] = (4.0/3.0)*(1.0 + pt[0])*(0.5 + pt[0])*(1.0 - pt[0]);
  na[3] = -(2.0/3.0)*(1.0 + pt[0])*(0.5 + pt[0])*(0.5 - pt[0]);

  double nb[4];
  nb[0] = -(2.0/3.0)*(0.5 + pt[1])*(0.5 - pt[1])*(1.0 - pt[1]);
  nb[1] = (4.0/3.0)*(1.0 + pt[1])*(0.5 - pt[1])*(1.0 - pt[1]);
  nb[2] = (4.0/3.0)*(1.0 + pt[1])*(0.5 + pt[1])*(1.0 - pt[1]);
  nb[3] = -(2.0/3.0)*(1.0 + pt[1])*(0.5 + pt[1])*(0.5 - pt[1]);

  double nc[4];
  nc[0] = -(2.0/3.0)*(0.5 + pt[2])*(0.5 - pt[2])*(1.0 - pt[2]);
  nc[1] = (4.0/3.0)*(1.0 + pt[2])*(0.5 - pt[2])*(1.0 - pt[2]);
  nc[2] = (4.0/3.0)*(1.0 + pt[2])*(0.5 + pt[2])*(1.0 - pt[2]);
  nc[3] = -(2.0/3.0)*(1.0 + pt[2])*(0.5 + pt[2])*(0.5 - pt[2]);

  for ( int k = 0; k < 4; k++ ){
    for ( int j = 0; j < 4; j++ ){
      for ( int i = 0; i < 4; i++ ){
        N[0] = na[i]*nb[j]*nc[k];
        N++;
      }
    }
  }
}

void TACSCubicHexaBasis::computeBasisGradient( const double pt[],
                                               double N[],
                                               double Nxi[] ){
  double na[4];
  na[0] = -(2.0/3.0)*(0.5 + pt[0])*(0.5 - pt[0])*(1.0 - pt[0]);
  na[1] = (4.0/3.0)*(1.0 + pt[0])*(0.5 - pt[0])*(1.0 - pt[0]);
  na[2] = (4.0/3.0)*(1.0 + pt[0])*(0.5 + pt[0])*(1.0 - pt[0]);
  na[3] = -(2.0/3.0)*(1.0 + pt[0])*(0.5 + pt[0])*(0.5 - pt[0]);

  double nb[4];
  nb[0] = -(2.0/3.0)*(0.5 + pt[1])*(0.5 - pt[1])*(1.0 - pt[1]);
  nb[1] = (4.0/3.0)*(1.0 + pt[1])*(0.5 - pt[1])*(1.0 - pt[1]);
  nb[2] = (4.0/3.0)*(1.0 + pt[1])*(0.5 + pt[1])*(1.0 - pt[1]);
  nb[3] = -(2.0/3.0)*(1.0 + pt[1])*(0.5 + pt[1])*(0.5 - pt[1]);

  double nc[4];
  nc[0] = -(2.0/3.0)*(0.5 + pt[2])*(0.5 - pt[2])*(1.0 - pt[2]);
  nc[1] = (4.0/3.0)*(1.0 + pt[2])*(0.5 - pt[2])*(1.0 - pt[2]);
  nc[2] = (4.0/3.0)*(1.0 + pt[2])*(0.5 + pt[2])*(1.0 - pt[2]);
  nc[3] = -(2.0/3.0)*(1.0 + pt[2])*(0.5 + pt[2])*(0.5 - pt[2]);

  double dna[4];
  dna[0] = -2.0*pt[0]*pt[0] + (4.0/3.0)*pt[0] + 1.0/6.0;
  dna[1] = 4.0*pt[0]*pt[0] - (4.0/3.0)*pt[0] - 4.0/3.0;
  dna[2] = -4.0*pt[0]*pt[0] - (4.0/3.0)*pt[0] + 4.0/3.0;
  dna[3] = 2.0*pt[0]*pt[0] + (4.0/3.0)*pt[0] - 1.0/6.0;

  double dnb[4];
  dnb[0] = -2.0*pt[1]*pt[1] + (4.0/3.0)*pt[1] + 1.0/6.0;
  dnb[1] = 4.0*pt[1]*pt[1] - (4.0/3.0)*pt[1] - 4.0/3.0;
  dnb[2] = -4.0*pt[1]*pt[1] - (4.0/3.0)*pt[1] + 4.0/3.0;
  dnb[3] = 2.0*pt[1]*pt[1] + (4.0/3.0)*pt[1] - 1.0/6.0;

  double dnc[4];
  dnc[0] = -2.0*pt[2]*pt[2] + (4.0/3.0)*pt[2] + 1.0/6.0;
  dnc[1] = 4.0*pt[2]*pt[2] - (4.0/3.0)*pt[2] - 4.0/3.0;
  dnc[2] = -4.0*pt[2]*pt[2] - (4.0/3.0)*pt[2] + 4.0/3.0;
  dnc[3] = 2.0*pt[2]*pt[2] + (4.0/3.0)*pt[2] - 1.0/6.0;

  for ( int k = 0; k < 4; k++ ){
    for ( int j = 0; j < 4; j++ ){
      for ( int i = 0; i < 4; i++ ){
        N[0] = na[i]*nb[j]*nc[k];
        Nxi[0] = dna[i]*nb[j]*nc[k];
        Nxi[1] = na[i]*dnb[j]*nc[k];
        Nxi[2] = na[i]*nb[j]*dnc[k];
        N++;
        Nxi += 3;
      }
    }
  }
}

void TACSCubicHexaBasis::interpFields( const int n,
                                       const double pt[],
                                       const int m,
                                       const TacsScalar v[],
                                       const int incr,
                                       TacsScalar u[] ){
  double na[4];
  na[0] = -(2.0/3.0)*(0.5 + pt[0])*(0.5 - pt[0])*(1.0 - pt[0]);
  na[1] = (4.0/3.0)*(1.0 + pt[0])*(0.5 - pt[0])*(1.0 - pt[0]);
  na[2] = (4.0/3.0)*(1.0 + pt[0])*(0.5 + pt[0])*(1.0 - pt[0]);
  na[3] = -(2.0/3.0)*(1.0 + pt[0])*(0.5 + pt[0])*(0.5 - pt[0]);

  double nb[4];
  nb[0] = -(2.0/3.0)*(0.5 + pt[1])*(0.5 - pt[1])*(1.0 - pt[1]);
  nb[1] = (4.0/3.0)*(1.0 + pt[1])*(0.5 - pt[1])*(1.0 - pt[1]);
  nb[2] = (4.0/3.0)*(1.0 + pt[1])*(0.5 + pt[1])*(1.0 - pt[1]);
  nb[3] = -(2.0/3.0)*(1.0 + pt[1])*(0.5 + pt[1])*(0.5 - pt[1]);

  double nc[4];
  nc[0] = -(2.0/3.0)*(0.5 + pt[2])*(0.5 - pt[2])*(1.0 - pt[2]);
  nc[1] = (4.0/3.0)*(1.0 + pt[2])*(0.5 - pt[2])*(1.0 - pt[2]);
  nc[2] = (4.0/3.0)*(1.0 + pt[2])*(0.5 + pt[2])*(1.0 - pt[2]);
  nc[3] = -(2.0/3.0)*(1.0 + pt[2])*(0.5 + pt[2])*(0.5 - pt[2]);

  for ( int i = 0; i < m; i++, u += incr, v++ ){
    u[0] = TACS_BASIS_EVAL_TENSOR3D_ORDER4(na, nb, nc, m, v);
  }
}

void TACSCubicHexaBasis::addInterpFieldsTranspose( const int n,
                                                   const double pt[],
                                                   const int incr,
                                                   const TacsScalar u[],
                                                   const int m,
                                                   TacsScalar v[] ){
  double na[4];
  na[0] = -(2.0/3.0)*(0.5 + pt[0])*(0.5 - pt[0])*(1.0 - pt[0]);
  na[1] = (4.0/3.0)*(1.0 + pt[0])*(0.5 - pt[0])*(1.0 - pt[0]);
  na[2] = (4.0/3.0)*(1.0 + pt[0])*(0.5 + pt[0])*(1.0 - pt[0]);
  na[3] = -(2.0/3.0)*(1.0 + pt[0])*(0.5 + pt[0])*(0.5 - pt[0]);

  double nb[4];
  nb[0] = -(2.0/3.0)*(0.5 + pt[1])*(0.5 - pt[1])*(1.0 - pt[1]);
  nb[1] = (4.0/3.0)*(1.0 + pt[1])*(0.5 - pt[1])*(1.0 - pt[1]);
  nb[2] = (4.0/3.0)*(1.0 + pt[1])*(0.5 + pt[1])*(1.0 - pt[1]);
  nb[3] = -(2.0/3.0)*(1.0 + pt[1])*(0.5 + pt[1])*(0.5 - pt[1]);

  double nc[4];
  nc[0] = -(2.0/3.0)*(0.5 + pt[2])*(0.5 - pt[2])*(1.0 - pt[2]);
  nc[1] = (4.0/3.0)*(1.0 + pt[2])*(0.5 - pt[2])*(1.0 - pt[2]);
  nc[2] = (4.0/3.0)*(1.0 + pt[2])*(0.5 + pt[2])*(1.0 - pt[2]);
  nc[3] = -(2.0/3.0)*(1.0 + pt[2])*(0.5 + pt[2])*(0.5 - pt[2]);

  for ( int i = 0; i < m; i++, u += incr, v++ ){
    TacsScalar temp;
    TACS_BASIS_TRANSPOSE_TENSOR3D_ORDER4(na, nb, nc, u[0], temp, m, v);
  }
}

void TACSCubicHexaBasis::interpFieldsGrad( const int n,
                                           const double pt[],
                                           const int m,
                                           const TacsScalar v[],
                                           TacsScalar g[] ){
  double na[4];
  na[0] = -(2.0/3.0)*(0.5 + pt[0])*(0.5 - pt[0])*(1.0 - pt[0]);
  na[1] = (4.0/3.0)*(1.0 + pt[0])*(0.5 - pt[0])*(1.0 - pt[0]);
  na[2] = (4.0/3.0)*(1.0 + pt[0])*(0.5 + pt[0])*(1.0 - pt[0]);
  na[3] = -(2.0/3.0)*(1.0 + pt[0])*(0.5 + pt[0])*(0.5 - pt[0]);

  double nb[4];
  nb[0] = -(2.0/3.0)*(0.5 + pt[1])*(0.5 - pt[1])*(1.0 - pt[1]);
  nb[1] = (4.0/3.0)*(1.0 + pt[1])*(0.5 - pt[1])*(1.0 - pt[1]);
  nb[2] = (4.0/3.0)*(1.0 + pt[1])*(0.5 + pt[1])*(1.0 - pt[1]);
  nb[3] = -(2.0/3.0)*(1.0 + pt[1])*(0.5 + pt[1])*(0.5 - pt[1]);

  double nc[4];
  nc[0] = -(2.0/3.0)*(0.5 + pt[2])*(0.5 - pt[2])*(1.0 - pt[2]);
  nc[1] = (4.0/3.0)*(1.0 + pt[2])*(0.5 - pt[2])*(1.0 - pt[2]);
  nc[2] = (4.0/3.0)*(1.0 + pt[2])*(0.5 + pt[2])*(1.0 - pt[2]);
  nc[3] = -(2.0/3.0)*(1.0 + pt[2])*(0.5 + pt[2])*(0.5 - pt[2]);

  double dna[4];
  dna[0] = -2.0*pt[0]*pt[0] + (4.0/3.0)*pt[0] + 1.0/6.0;
  dna[1] = 4.0*pt[0]*pt[0] - (4.0/3.0)*pt[0] - 4.0/3.0;
  dna[2] = -4.0*pt[0]*pt[0] - (4.0/3.0)*pt[0] + 4.0/3.0;
  dna[3] = 2.0*pt[0]*pt[0] + (4.0/3.0)*pt[0] - 1.0/6.0;

  double dnb[4];
  dnb[0] = -2.0*pt[1]*pt[1] + (4.0/3.0)*pt[1] + 1.0/6.0;
  dnb[1] = 4.0*pt[1]*pt[1] - (4.0/3.0)*pt[1] - 4.0/3.0;
  dnb[2] = -4.0*pt[1]*pt[1] - (4.0/3.0)*pt[1] + 4.0/3.0;
  dnb[3] = 2.0*pt[1]*pt[1] + (4.0/3.0)*pt[1] - 1.0/6.0;

  double dnc[4];
  dnc[0] = -2.0*pt[2]*pt[2] + (4.0/3.0)*pt[2] + 1.0/6.0;
  dnc[1] = 4.0*pt[2]*pt[2] - (4.0/3.0)*pt[2] - 4.0/3.0;
  dnc[2] = -4.0*pt[2]*pt[2] - (4.0/3.0)*pt[2] + 4.0/3.0;
  dnc[3] = 2.0*pt[2]*pt[2] + (4.0/3.0)*pt[2] - 1.0/6.0;

  for ( int i = 0; i < m; i++, g += 3, v++ ){
    g[0] = TACS_BASIS_EVAL_TENSOR3D_ORDER4(dna, nb, nc, m, v);
    g[1] = TACS_BASIS_EVAL_TENSOR3D_ORDER4(na, dnb, nc, m, v);
    g[2] = TACS_BASIS_EVAL_TENSOR3D_ORDER4(na, nb, dnc, m, v);
  }
}

void TACSCubicHexaBasis::addInterpFieldsGradTranspose( int n,
                                                         const double pt[],
                                                         const int m,
                                                         const TacsScalar g[],
                                                         TacsScalar v[] ){
  double na[4];
  na[0] = -(2.0/3.0)*(0.5 + pt[0])*(0.5 - pt[0])*(1.0 - pt[0]);
  na[1] = (4.0/3.0)*(1.0 + pt[0])*(0.5 - pt[0])*(1.0 - pt[0]);
  na[2] = (4.0/3.0)*(1.0 + pt[0])*(0.5 + pt[0])*(1.0 - pt[0]);
  na[3] = -(2.0/3.0)*(1.0 + pt[0])*(0.5 + pt[0])*(0.5 - pt[0]);

  double nb[4];
  nb[0] = -(2.0/3.0)*(0.5 + pt[1])*(0.5 - pt[1])*(1.0 - pt[1]);
  nb[1] = (4.0/3.0)*(1.0 + pt[1])*(0.5 - pt[1])*(1.0 - pt[1]);
  nb[2] = (4.0/3.0)*(1.0 + pt[1])*(0.5 + pt[1])*(1.0 - pt[1]);
  nb[3] = -(2.0/3.0)*(1.0 + pt[1])*(0.5 + pt[1])*(0.5 - pt[1]);

  double nc[4];
  nc[0] = -(2.0/3.0)*(0.5 + pt[2])*(0.5 - pt[2])*(1.0 - pt[2]);
  nc[1] = (4.0/3.0)*(1.0 + pt[2])*(0.5 - pt[2])*(1.0 - pt[2]);
  nc[2] = (4.0/3.0)*(1.0 + pt[2])*(0.5 + pt[2])*(1.0 - pt[2]);
  nc[3] = -(2.0/3.0)*(1.0 + pt[2])*(0.5 + pt[2])*(0.5 - pt[2]);

  double dna[4];
  dna[0] = -2.0*pt[0]*pt[0] + (4.0/3.0)*pt[0] + 1.0/6.0;
  dna[1] = 4.0*pt[0]*pt[0] - (4.0/3.0)*pt[0] - 4.0/3.0;
  dna[2] = -4.0*pt[0]*pt[0] - (4.0/3.0)*pt[0] + 4.0/3.0;
  dna[3] = 2.0*pt[0]*pt[0] + (4.0/3.0)*pt[0] - 1.0/6.0;

  double dnb[4];
  dnb[0] = -2.0*pt[1]*pt[1] + (4.0/3.0)*pt[1] + 1.0/6.0;
  dnb[1] = 4.0*pt[1]*pt[1] - (4.0/3.0)*pt[1] - 4.0/3.0;
  dnb[2] = -4.0*pt[1]*pt[1] - (4.0/3.0)*pt[1] + 4.0/3.0;
  dnb[3] = 2.0*pt[1]*pt[1] + (4.0/3.0)*pt[1] - 1.0/6.0;

  double dnc[4];
  dnc[0] = -2.0*pt[2]*pt[2] + (4.0/3.0)*pt[2] + 1.0/6.0;
  dnc[1] = 4.0*pt[2]*pt[2] - (4.0/3.0)*pt[2] - 4.0/3.0;
  dnc[2] = -4.0*pt[2]*pt[2] - (4.0/3.0)*pt[2] + 4.0/3.0;
  dnc[3] = 2.0*pt[2]*pt[2] + (4.0/3.0)*pt[2] - 1.0/6.0;

  for ( int i = 0; i < m; i++, g += 3, v++ ){
    TacsScalar temp;
    TACS_BASIS_TRANSPOSE_TENSOR3D_ORDER4(dna, nb, nc, g[0], temp, m, v);
    TACS_BASIS_TRANSPOSE_TENSOR3D_ORDER4(na, dnb, nc, g[1], temp, m, v);
    TACS_BASIS_TRANSPOSE_TENSOR3D_ORDER4(na, nb, dnc, g[2], temp, m, v);
  }
}
/*
  Quartic hexahedral basis class functions
*/
TACSQuarticHexaBasis::TACSQuarticHexaBasis(){
  for ( int i = 0; i < 5; i++ ){
    TacsLagrangeShapeFuncDerivative(5, TacsGaussQuadPts5[i], cosine_pts,
                                    &Nf[5*i], &Nfxi[5*i]);
  }
}

const double TACSQuarticHexaBasis::cosine_pts[5] =
  {-1.0, -0.7071067811865475, 0.0, 0.7071067811865475, 1.0};

ElementLayout TACSQuarticHexaBasis::getLayoutType(){
  return TACS_HEXA_QUARTIC_ELEMENT;
}

void TACSQuarticHexaBasis::getVisPoint( int n, double pt[] ){
  pt[0] = cosine_pts[(n % 25) % 5];
  pt[1] = cosine_pts[(n / 25) % 5];
  pt[2] = cosine_pts[n / 25];
}

int TACSQuarticHexaBasis::getNumNodes(){
  return 125;
}

int TACSQuarticHexaBasis::getNumParameters(){
  return 3;
}

int TACSQuarticHexaBasis::getNumQuadraturePoints(){
  return 125;
}

double TACSQuarticHexaBasis::getQuadratureWeight( int n ){
  return
    (TacsGaussQuadWts5[n % 5]*
     TacsGaussQuadWts5[(n % 25) / 5]*
     TacsGaussQuadWts5[n / 25]);
}

double TACSQuarticHexaBasis::getQuadraturePoint( int n,
                                                 double pt[] ){
  pt[0] = TacsGaussQuadPts5[n % 5];
  pt[1] = TacsGaussQuadPts5[(n % 25) / 5];
  pt[2] = TacsGaussQuadPts5[n / 25];

  return
    (TacsGaussQuadWts5[n % 5]*
     TacsGaussQuadWts5[(n % 25) / 5]*
     TacsGaussQuadWts5[n / 25]);
}

int TACSQuarticHexaBasis::getNumElementFaces(){
  return 6;
}

int TACSQuarticHexaBasis::getNumFaceQuadraturePoints( int face ){
  return 25;
}

double TACSQuarticHexaBasis::getFaceQuadraturePoint( int face,
                                                     int n,
                                                     double pt[],
                                                     double t[] ){
  if (face/2 == 0){
    pt[0] = -1.0 + 2.0*(face % 2);
    pt[1] = TacsGaussQuadPts5[n % 5];
    pt[2] = TacsGaussQuadPts5[n / 5];
  }
  else if (face/2 == 1){
    pt[0] = TacsGaussQuadPts5[n % 5];
    pt[1] = -1.0 + 2.0*(face % 2);
    pt[2] = TacsGaussQuadPts5[n / 5];
  }
  else {
    pt[0] = TacsGaussQuadPts5[n % 5];
    pt[1] = TacsGaussQuadPts5[n / 5];
    pt[2] = -1.0 + 2.0*(face % 2);
  }

  getFaceTangents(face, t);

  return (TacsGaussQuadWts5[n % 5]*TacsGaussQuadWts5[n / 5]);
}

void TACSQuarticHexaBasis::computeBasis( const double pt[],
                                         double N[] ){
  double na[5], nb[5], nc[5];
  TacsLagrangeShapeFunctions(5, pt[0], cosine_pts, na);
  TacsLagrangeShapeFunctions(5, pt[1], cosine_pts, nb);
  TacsLagrangeShapeFunctions(5, pt[2], cosine_pts, nc);

  for ( int k = 0; k < 5; k++ ){
    for ( int j = 0; j < 5; j++ ){
      for ( int i = 0; i < 5; i++ ){
        N[i + 5*j + 25*k] = na[i]*nb[j]*nc[k];
      }
    }
  }
}

void TACSQuarticHexaBasis::computeBasisGradient( const double pt[],
                                                 double N[],
                                                 double Nxi[] ){
  double na[5], nb[5], nc[5];
  double dna[5], dnb[5], dnc[5];
  TacsLagrangeShapeFuncDerivative(5, pt[0], cosine_pts, na, dna);
  TacsLagrangeShapeFuncDerivative(5, pt[1], cosine_pts, nb, dnb);
  TacsLagrangeShapeFuncDerivative(5, pt[2], cosine_pts, nc, dnc);

  for ( int k = 0; k < 5; k++ ){
    for ( int j = 0; j < 5; j++ ){
      for ( int i = 0; i < 5; i++ ){
        N[i + 5*j + 25*k] = na[i]*nb[j]*nc[k];
        Nxi[3*(i + 5*j + 25*k)] = dna[i]*nb[j]*nc[k];
        Nxi[3*(i + 5*j + 25*k)+1] = na[i]*dnb[j]*nc[k];
        Nxi[3*(i + 5*j + 25*k)+2] = na[i]*nb[j]*dnc[k];
      }
    }
  }
}

void TACSQuarticHexaBasis::interpFields( const int n,
                                         const double pt[],
                                         const int m,
                                         const TacsScalar v[],
                                         const int incr,
                                         TacsScalar u[] ){
  if (n >= 0){
    const double *n1 = &Nf[5*(n % 5)];
    const double *n2 = &Nf[5*((n % 25)/5)];
    const double *n3 = &Nf[5*(n / 25)];

    for ( int k = 0; k < 5; k++ ){
      for ( int j = 0; j < 5; j++ ){
        TacsScalar n23 = n2[j]*n3[k];
        for ( int p = 0; p < m; p++ ){
          TacsScalar t1 = (n1[0]*v[0]   + n1[1]*v[m]   + n1[2]*v[2*m] +
                           n1[3]*v[3*m] + n1[4]*v[4*m]);
          u[p*incr] += n23*t1;
          v++;
        }
        v += 4*m;
      }
    }
  }
  else {
    double n1[5], n2[5], n3[5];
    TacsLagrangeShapeFunctions(5, pt[0], cosine_pts, n1);
    TacsLagrangeShapeFunctions(5, pt[1], cosine_pts, n2);
    TacsLagrangeShapeFunctions(5, pt[2], cosine_pts, n3);

    for ( int k = 0; k < 5; k++ ){
      for ( int j = 0; j < 5; j++ ){
        TacsScalar n23 = n2[j]*n3[k];
        for ( int p = 0; p < m; p++ ){
          TacsScalar t1 = (n1[0]*v[0]   + n1[1]*v[m]   + n1[2]*v[2*m] +
                           n1[3]*v[3*m] + n1[4]*v[4*m]);
          u[p*incr] += n23*t1;
          v++;
        }
        v += 4*m;
      }
    }
  }
}

void TACSQuarticHexaBasis::addInterpFieldsTranspose( const int n,
                                                     const double pt[],
                                                     const int incr,
                                                     const TacsScalar u[],
                                                     const int m,
                                                     TacsScalar v[] ){
  if (n >= 0){
    const double *n1 = &Nf[5*(n % 5)];
    const double *n2 = &Nf[5*((n % 25)/5)];
    const double *n3 = &Nf[5*(n / 25)];

    for ( int k = 0; k < 5; k++ ){
      for ( int j = 0; j < 5; j++ ){
        TacsScalar n23 = n2[j]*n3[k];
        for ( int p = 0; p < m; p++ ){
          TacsScalar a = n23*u[p*incr];
          v[0] +=   a*n1[0];
          v[m] +=   a*n1[1];
          v[2*m] += a*n1[2];
          v[3*m] += a*n1[3];
          v[4*m] += a*n1[4];
          v++;
        }
        v += 4*m;
      }
    }
  }
  else {
    double n1[5], n2[5], n3[5];
    TacsLagrangeShapeFunctions(5, pt[0], cosine_pts, n1);
    TacsLagrangeShapeFunctions(5, pt[1], cosine_pts, n2);
    TacsLagrangeShapeFunctions(5, pt[2], cosine_pts, n3);

    for ( int k = 0; k < 5; k++ ){
      for ( int j = 0; j < 5; j++ ){
        TacsScalar n23 = n2[j]*n3[k];
        for ( int p = 0; p < m; p++ ){
          TacsScalar a = n23*u[p*incr];
          v[0] +=   a*n1[0];
          v[m] +=   a*n1[1];
          v[2*m] += a*n1[2];
          v[3*m] += a*n1[3];
          v[4*m] += a*n1[4];
          v++;
        }
        v += 4*m;
      }
    }
  }
}

void TACSQuarticHexaBasis::interpFieldsGrad( const int n,
                                             const double pt[],
                                             const int m,
                                             const TacsScalar v[],
                                             TacsScalar grad[] ){
  if (n >= 0){
    const double *n1 = &Nf[5*(n % 5)];
    const double *n2 = &Nf[5*((n % 25)/5)];
    const double *n3 = &Nf[5*(n / 25)];
    const double *n1x = &Nfxi[5*(n % 5)];
    const double *n2x = &Nfxi[5*((n % 25)/5)];
    const double *n3x = &Nfxi[5*(n / 25)];

    memset(grad, 0, 3*m*sizeof(TacsScalar));

    for ( int k = 0; k < 5; k++ ){
      for ( int j = 0; j < 5; j++ ){
        TacsScalar n23  = n2[j]*n3[k];
        TacsScalar n2x3 = n2x[j]*n3[k];
        TacsScalar n23x = n2[j]*n3x[k];

        TacsScalar *g = grad;
        for ( int p = 0; p < m; p++ ){
          g[0] += n23*(n1x[0]*v[0]   + n1x[1]*v[m]   + n1x[2]*v[2*m] +
                       n1x[3]*v[3*m] + n1x[4]*v[4*m]);
          TacsScalar t1 = (n1[0]*v[0]   + n1[1]*v[m]   + n1[2]*v[2*m] +
                           n1[3]*v[3*m] + n1[4]*v[4*m]);
          g[1] += n2x3*t1;
          g[2] += n23x*t1;
          g += 3;
          v++;
        }
        v += 4*m;
      }
    }
  }
  else {
    double n1[5], n2[5], n3[5], n1x[5], n2x[5], n3x[5];
    TacsLagrangeShapeFuncDerivative(5, pt[0], cosine_pts, n1, n1x);
    TacsLagrangeShapeFuncDerivative(5, pt[1], cosine_pts, n2, n2x);
    TacsLagrangeShapeFuncDerivative(5, pt[2], cosine_pts, n3, n3x);

    for ( int k = 0; k < 5; k++ ){
      for ( int j = 0; j < 5; j++ ){
        TacsScalar n23  = n2[j]*n3[k];
        TacsScalar n2x3 = n2x[j]*n3[k];
        TacsScalar n23x = n2[j]*n3x[k];

        TacsScalar *g = grad;
        for ( int p = 0; p < m; p++ ){
          g[0] += n23*(n1x[0]*v[0]   + n1x[1]*v[m]   + n1x[2]*v[2*m] +
                       n1x[3]*v[3*m] + n1x[4]*v[4*m]);
          TacsScalar t1 = (n1[0]*v[0]   + n1[1]*v[m]   + n1[2]*v[2*m] +
                           n1[3]*v[3*m] + n1[4]*v[4*m]);
          g[1] += n2x3*t1;
          g[2] += n23x*t1;
          g += 3;
          v++;
        }
        v += 4*m;
      }
    }
  }
}

void TACSQuarticHexaBasis::addInterpFieldsGradTranspose( int n,
                                                         const double pt[],
                                                         const int m,
                                                         const TacsScalar grad[],
                                                         TacsScalar v[] ){
  if (n >= 0){
    const double *n1 = &Nf[5*(n % 5)];
    const double *n2 = &Nf[5*((n % 25)/5)];
    const double *n3 = &Nf[5*(n / 25)];
    const double *n1x = &Nfxi[5*(n % 5)];
    const double *n2x = &Nfxi[5*((n % 25)/5)];
    const double *n3x = &Nfxi[5*(n / 25)];

    for ( int k = 0; k < 5; k++ ){
      for ( int j = 0; j < 5; j++ ){
        TacsScalar n23  = n2[j]*n3[k];
        TacsScalar n2x3 = n2x[j]*n3[k];
        TacsScalar n23x = n2[j]*n3x[k];

        const TacsScalar *g = grad;
        for ( int p = 0; p < m; p++ ){
          TacsScalar a = n2x3*g[1] + n23x*g[2];
          TacsScalar b = n23*g[0];
          v[0] +=   a*n1[0] + b*n1x[0];
          v[m] +=   a*n1[1] + b*n1x[1];
          v[2*m] += a*n1[2] + b*n1x[2];
          v[3*m] += a*n1[3] + b*n1x[3];
          v[4*m] += a*n1[4] + b*n1x[4];
          g += 3;
          v++;
        }
        v += 4*m;
      }
    }
  }
  else {
    double n1[5], n2[5], n3[5], n1x[5], n2x[5], n3x[5];
    TacsLagrangeShapeFuncDerivative(5, pt[0], cosine_pts, n1, n1x);
    TacsLagrangeShapeFuncDerivative(5, pt[1], cosine_pts, n2, n2x);
    TacsLagrangeShapeFuncDerivative(5, pt[2], cosine_pts, n3, n3x);

    for ( int k = 0; k < 5; k++ ){
      for ( int j = 0; j < 5; j++ ){
        TacsScalar n23  = n2[j]*n3[k];
        TacsScalar n2x3 = n2x[j]*n3[k];
        TacsScalar n23x = n2[j]*n3x[k];

        const TacsScalar *g = grad;
        for ( int p = 0; p < m; p++ ){
          TacsScalar a = n2x3*g[1] + n23x*g[2];
          TacsScalar b = n23*g[0];
          v[0] +=   a*n1[0] + b*n1x[0];
          v[m] +=   a*n1[1] + b*n1x[1];
          v[2*m] += a*n1[2] + b*n1x[2];
          v[3*m] += a*n1[3] + b*n1x[3];
          v[4*m] += a*n1[4] + b*n1x[4];
          g += 3;
          v++;
        }
        v += 4*m;
      }
    }
  }
}

void TACSQuarticHexaBasis::interpAllFieldsGrad( const int m,
                                                const TacsScalar values[],
                                                TacsScalar out[] ){
  memset(out, 0, 4*m*125*sizeof(TacsScalar));

  for ( int n = 0; n < 125; n++ ){
    const double *n1 = &Nf[5*(n % 5)];
    const double *n2 = &Nf[5*((n % 25)/5)];
    const double *n3 = &Nf[5*(n / 25)];
    const double *n1x = &Nfxi[5*(n % 5)];
    const double *n2x = &Nfxi[5*((n % 25)/5)];
    const double *n3x = &Nfxi[5*(n / 25)];

    const TacsScalar *v = values;

    if (m == 1){
      for ( int k = 0; k < 5; k++ ){
        for ( int j = 0; j < 5; j++ ){
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          out[1] += n23*(n1x[0]*v[0] + n1x[1]*v[1] + n1x[2]*v[2] +
                         n1x[3]*v[3] + n1x[4]*v[4]);
          TacsScalar t1 = (n1[0]*v[0] + n1[1]*v[1] + n1[2]*v[2] +
                           n1[3]*v[3] + n1[4]*v[4]);
          out[0] += n23*t1;
          out[2] += n2x3*t1;
          out[3] += n23x*t1;
          v += 5;
        }
      }
    }
    else if (m == 3){
      for ( int k = 0; k < 5; k++ ){
        for ( int j = 0; j < 5; j++ ){
          TacsScalar t1;
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          out[3] += n23*(n1x[0]*v[0] + n1x[1]*v[3]  + n1x[2]*v[6] +
                         n1x[3]*v[9] + n1x[4]*v[12]);
          t1 = (n1[0]*v[0] + n1[1]*v[3]  + n1[2]*v[6] +
                n1[3]*v[9] + n1[4]*v[12]);
          out[0] += n23*t1;
          out[4] += n2x3*t1;
          out[5] += n23x*t1;

          out[6] += n23*(n1x[0]*v[1]  + n1x[1]*v[4]  + n1x[2]*v[7] +
                         n1x[3]*v[10] + n1x[4]*v[13]);
          t1 = (n1[0]*v[1]  + n1[1]*v[4]  + n1[2]*v[7] +
                n1[3]*v[10] + n1[4]*v[13]);
          out[1] += n23*t1;
          out[7] += n2x3*t1;
          out[8] += n23x*t1;

          out[9] += n23*(n1x[0]*v[2]  + n1x[1]*v[5]  + n1x[2]*v[8] +
                         n1x[3]*v[11] + n1x[4]*v[14]);
          t1 = (n1[0]*v[2]  + n1[1]*v[5]   + n1[2]*v[8] +
                n1[3]*v[11] + n1[4]*v[14]);
          out[2] += n23*t1;
          out[10] += n2x3*t1;
          out[11] += n23x*t1;

          v += 15;
        }
      }
    }
    else if (m == 4){
      for ( int k = 0; k < 5; k++ ){
        for ( int j = 0; j < 5; j++ ){
          TacsScalar t1;
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          out[4] += n23*(n1x[0]*v[0]  + n1x[1]*v[4]  + n1x[2]*v[8] +
                         n1x[3]*v[12] + n1x[4]*v[16]);
          t1 = (n1[0]*v[0] + n1[1]*v[4]  + n1[2]*v[8] +
                n1[3]*v[12] + n1[4]*v[16]);
          out[0] += n23*t1;
          out[5] += n2x3*t1;
          out[6] += n23x*t1;

          out[7] += n23*(n1x[0]*v[1]  + n1x[1]*v[5]  + n1x[2]*v[9] +
                         n1x[3]*v[13] + n1x[4]*v[17]);
          t1 = (n1[0]*v[1]  + n1[1]*v[5]  + n1[2]*v[9] +
                n1[3]*v[13] + n1[4]*v[17]);
          out[1] += n23*t1;
          out[8] += n2x3*t1;
          out[9] += n23x*t1;

          out[10] += n23*(n1x[0]*v[2]  + n1x[1]*v[6]  + n1x[2]*v[10] +
                          n1x[3]*v[14] + n1x[4]*v[18]);
          t1 = (n1[0]*v[2]  + n1[1]*v[6]  + n1[2]*v[10] +
                n1[3]*v[14] + n1[4]*v[18]);
          out[2] += n23*t1;
          out[11] += n2x3*t1;
          out[12] += n23x*t1;

          out[13] += n23*(n1x[0]*v[3]  + n1x[1]*v[7]  + n1x[2]*v[11] +
                         n1x[3]*v[15] + n1x[4]*v[19]);
          t1 = (n1[0]*v[3] + n1[1]*v[7]  + n1[2]*v[11] +
                n1[3]*v[15] + n1[4]*v[19]);
          out[3] += n23*t1;
          out[14] += n2x3*t1;
          out[15] += n23x*t1;

          v += 20;
        }
      }
    }
    else {
      for ( int k = 0; k < 5; k++ ){
        for ( int j = 0; j < 5; j++ ){
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          TacsScalar *g = &out[m];
          for ( int p = 0; p < m; p++ ){
            g[0] += n23*(n1x[0]*v[0]   + n1x[1]*v[m]   + n1x[2]*v[2*m] +
                         n1x[3]*v[3*m] + n1x[4]*v[4*m]);
            TacsScalar t1 = (n1[0]*v[0]   + n1[1]*v[m]   + n1[2]*v[2*m] +
                             n1[3]*v[3*m] + n1[4]*v[4*m]);
            out[p] += n23*t1;
            g[1] += n2x3*t1;
            g[2] += n23x*t1;
            g += 3;
            v++;
          }
          v += 4*m;
        }
      }
    }

    out += 4*m;
  }
}

void TACSQuarticHexaBasis::addInterpAllFieldsGradTranspose( const int m,
                                                            const TacsScalar in[],
                                                            TacsScalar values[] ){
  for ( int n = 0; n < 125; n++ ){
    const double *n1 = &Nf[5*(n % 5)];
    const double *n2 = &Nf[5*((n % 25)/5)];
    const double *n3 = &Nf[5*(n / 25)];
    const double *n1x = &Nfxi[5*(n % 5)];
    const double *n2x = &Nfxi[5*((n % 25)/5)];
    const double *n3x = &Nfxi[5*(n / 25)];

    TacsScalar *v = values;

    if (m == 1){
      for ( int k = 0; k < 5; k++ ){
        for ( int j = 0; j < 5; j++ ){
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          TacsScalar a = n23*in[0] + n2x3*in[2] + n23x*in[3];
          TacsScalar b = n23*in[1];
          v[0] += a*n1[0] + b*n1x[0];
          v[1] += a*n1[1] + b*n1x[1];
          v[2] += a*n1[2] + b*n1x[2];
          v[3] += a*n1[3] + b*n1x[3];
          v[4] += a*n1[4] + b*n1x[4];
          v += 5;
        }
      }
    }
    else if (m == 3){
      for ( int k = 0; k < 5; k++ ){
        for ( int j = 0; j < 5; j++ ){
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          TacsScalar a1 = n23*in[0] + n2x3*in[4] + n23x*in[5];
          TacsScalar b1 = n23*in[3];

          TacsScalar a2 = n23*in[1] + n2x3*in[7] + n23x*in[8];
          TacsScalar b2 = n23*in[6];

          TacsScalar a3 = n23*in[2] + n2x3*in[10] + n23x*in[11];
          TacsScalar b3 = n23*in[9];

          v[0] +=  a1*n1[0] + b1*n1x[0];
          v[1] +=  a2*n1[0] + b2*n1x[0];
          v[2] +=  a3*n1[0] + b3*n1x[0];

          v[3] +=  a1*n1[1] + b1*n1x[1];
          v[4] +=  a2*n1[1] + b2*n1x[1];
          v[5] +=  a3*n1[1] + b3*n1x[1];

          v[6] +=  a1*n1[2] + b1*n1x[2];
          v[7] +=  a2*n1[2] + b2*n1x[2];
          v[8] +=  a3*n1[2] + b3*n1x[2];

          v[9] +=  a1*n1[3] + b1*n1x[3];
          v[10] += a2*n1[3] + b2*n1x[3];
          v[11] += a3*n1[3] + b3*n1x[3];

          v[12] += a1*n1[4] + b1*n1x[4];
          v[13] += a2*n1[4] + b2*n1x[4];
          v[14] += a3*n1[4] + b3*n1x[4];

          v += 15;
        }
      }
    }
    else if (m == 4){
      for ( int k = 0; k < 5; k++ ){
        for ( int j = 0; j < 5; j++ ){
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          TacsScalar a1 = n23*in[0] + n2x3*in[5] + n23x*in[6];
          TacsScalar b1 = n23*in[4];

          TacsScalar a2 = n23*in[1] + n2x3*in[8] + n23x*in[9];
          TacsScalar b2 = n23*in[7];

          TacsScalar a3 = n23*in[2] + n2x3*in[11] + n23x*in[12];
          TacsScalar b3 = n23*in[10];

          TacsScalar a4 = n23*in[3] + n2x3*in[14] + n23x*in[15];
          TacsScalar b4 = n23*in[13];

          v[0] +=  a1*n1[0] + b1*n1x[0];
          v[1] +=  a2*n1[0] + b2*n1x[0];
          v[2] +=  a3*n1[0] + b3*n1x[0];
          v[3] +=  a4*n1[0] + b4*n1x[0];

          v[4] +=  a1*n1[1] + b1*n1x[1];
          v[5] +=  a2*n1[1] + b2*n1x[1];
          v[6] +=  a3*n1[1] + b3*n1x[1];
          v[7] +=  a4*n1[1] + b4*n1x[1];

          v[8] +=  a1*n1[2] + b1*n1x[2];
          v[9] +=  a2*n1[2] + b2*n1x[2];
          v[10] += a3*n1[2] + b3*n1x[2];
          v[11] += a4*n1[2] + b4*n1x[2];

          v[12] += a1*n1[3] + b1*n1x[3];
          v[13] += a2*n1[3] + b2*n1x[3];
          v[14] += a3*n1[3] + b3*n1x[3];
          v[15] += a4*n1[3] + b4*n1x[3];

          v[16] += a1*n1[4] + b1*n1x[4];
          v[17] += a2*n1[4] + b2*n1x[4];
          v[18] += a3*n1[4] + b3*n1x[4];
          v[19] += a4*n1[4] + b4*n1x[4];

          v += 20;
        }
      }
    }
    else {
      for ( int k = 0; k < 5; k++ ){
        for ( int j = 0; j < 5; j++ ){
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          const TacsScalar *g = &in[m];
          for ( int p = 0; p < m; p++ ){
            TacsScalar a = n23*in[p] + n2x3*g[1] + n23x*g[2];
            TacsScalar b = n23*g[0];
            v[0] +=   a*n1[0] + b*n1x[0];
            v[m] +=   a*n1[1] + b*n1x[1];
            v[2*m] += a*n1[2] + b*n1x[2];
            v[3*m] += a*n1[3] + b*n1x[3];
            v[4*m] += a*n1[4] + b*n1x[4];
            g += 3;
            v++;
          }
          v += 4*m;
        }
      }
    }

    in += 4*m;
  }
}

/*
  Quintic Quad basis class functions
*/
TACSQuinticHexaBasis::TACSQuinticHexaBasis(){
  for ( int i = 0; i < 6; i++ ){
    TacsLagrangeShapeFuncDerivative(6, TacsGaussQuadPts6[i], cosine_pts,
                                    &Nf[6*i], &Nfxi[6*i]);
  }
}

const double TACSQuinticHexaBasis::cosine_pts[6] =
  {-1.0, -0.8090169943749475, -0.30901699437494745,
   0.30901699437494745, 0.8090169943749475, 1.0};

ElementLayout TACSQuinticHexaBasis::getLayoutType(){
  return TACS_HEXA_QUINTIC_ELEMENT;
}

void TACSQuinticHexaBasis::getVisPoint( int n, double pt[] ){
  pt[0] = cosine_pts[n % 6];
  pt[1] = cosine_pts[n/6];
}

int TACSQuinticHexaBasis::getNumNodes(){
  return 216;
}

int TACSQuinticHexaBasis::getNumParameters(){
  return 3;
}

int TACSQuinticHexaBasis::getNumQuadraturePoints(){
  return 216;
}

double TACSQuinticHexaBasis::getQuadratureWeight( int n ){
  return
    (TacsGaussQuadWts6[n % 6]*
     TacsGaussQuadWts6[(n % 36) / 6]*
     TacsGaussQuadWts6[n / 36]);
}

double TACSQuinticHexaBasis::getQuadraturePoint( int n,
                                                 double pt[] ){
  pt[0] = TacsGaussQuadPts6[n % 6];
  pt[1] = TacsGaussQuadPts6[(n % 36) / 6];
  pt[2] = TacsGaussQuadPts6[n / 36];

  return
    (TacsGaussQuadWts6[n % 6]*
     TacsGaussQuadWts6[(n % 36) / 6]*
     TacsGaussQuadWts6[n / 36]);
}

int TACSQuinticHexaBasis::getNumElementFaces(){
  return 6;
}

int TACSQuinticHexaBasis::getNumFaceQuadraturePoints( int face ){
  return 36;
}

double TACSQuinticHexaBasis::getFaceQuadraturePoint( int face,
                                                     int n,
                                                     double pt[],
                                                     double t[] ){
  if (face/2 == 0){
    pt[0] = -1.0 + 2.0*(face % 2);
    pt[1] = TacsGaussQuadPts6[n % 6];
    pt[2] = TacsGaussQuadPts6[n / 6];
  }
  else if (face/2 == 1){
    pt[0] = TacsGaussQuadPts6[n % 6];
    pt[1] = -1.0 + 2.0*(face % 2);
    pt[2] = TacsGaussQuadPts6[n / 6];
  }
  else {
    pt[0] = TacsGaussQuadPts6[n % 6];
    pt[1] = TacsGaussQuadPts6[n / 6];
    pt[2] = -1.0 + 2.0*(face % 2);
  }

  getFaceTangents(face, t);

  return (TacsGaussQuadWts6[n % 6]*TacsGaussQuadWts6[n / 6]);
}

void TACSQuinticHexaBasis::computeBasis( const double pt[],
                                         double N[] ){
  double na[6], nb[6], nc[6];
  TacsLagrangeShapeFunctions(6, pt[0], cosine_pts, na);
  TacsLagrangeShapeFunctions(6, pt[1], cosine_pts, nb);
  TacsLagrangeShapeFunctions(6, pt[2], cosine_pts, nc);

  for ( int k = 0; k < 6; k++ ){
    for ( int j = 0; j < 6; j++ ){
      for ( int i = 0; i < 6; i++ ){
        N[i + 6*j + 36*k] = na[i]*nb[j]*nc[k];
      }
    }
  }
}

void TACSQuinticHexaBasis::computeBasisGradient( const double pt[],
                                                 double N[],
                                                 double Nxi[] ){
  double na[6], nb[6], nc[6];
  double dna[6], dnb[6], dnc[6];
  TacsLagrangeShapeFuncDerivative(6, pt[0], cosine_pts, na, dna);
  TacsLagrangeShapeFuncDerivative(6, pt[1], cosine_pts, nb, dnb);
  TacsLagrangeShapeFuncDerivative(6, pt[2], cosine_pts, nc, dnc);

  for ( int k = 0; k < 6; k++ ){
    for ( int j = 0; j < 6; j++ ){
      for ( int i = 0; i < 6; i++ ){
        N[i + 6*j + 36*k] = na[i]*nb[j]*nc[k];
        Nxi[3*(i + 6*j + 36*k)] = dna[i]*nb[j]*nc[k];
        Nxi[3*(i + 6*j + 36*k)+1] = na[i]*dnb[j]*nc[k];
        Nxi[3*(i + 6*j + 36*k)+2] = na[i]*nb[j]*dnc[k];
      }
    }
  }
}

void TACSQuinticHexaBasis::interpFields( const int n,
                                         const double pt[],
                                         const int m,
                                         const TacsScalar v[],
                                         const int incr,
                                         TacsScalar u[] ){
  if (n >= 0){
    const double *n1 = &Nf[6*(n % 6)];
    const double *n2 = &Nf[6*((n % 36)/6)];
    const double *n3 = &Nf[6*(n / 36)];

    for ( int k = 0; k < 6; k++ ){
      for ( int j = 0; j < 6; j++ ){
        TacsScalar n23 = n2[j]*n3[k];
        for ( int p = 0; p < m; p++ ){
          TacsScalar t1 = (n1[0]*v[0]   + n1[1]*v[m]   + n1[2]*v[2*m] +
                           n1[3]*v[3*m] + n1[4]*v[4*m] + n1[5]*v[5*m]);
          u[p*incr] += n23*t1;
          v++;
        }
        v += 5*m;
      }
    }
  }
  else {
    double n1[6], n2[6], n3[6];
    TacsLagrangeShapeFunctions(6, pt[0], cosine_pts, n1);
    TacsLagrangeShapeFunctions(6, pt[1], cosine_pts, n2);
    TacsLagrangeShapeFunctions(6, pt[2], cosine_pts, n3);

    for ( int k = 0; k < 6; k++ ){
      for ( int j = 0; j < 6; j++ ){
        TacsScalar n23 = n2[j]*n3[k];
        for ( int p = 0; p < m; p++ ){
          TacsScalar t1 = (n1[0]*v[0]   + n1[1]*v[m]   + n1[2]*v[2*m] +
                           n1[3]*v[3*m] + n1[4]*v[4*m] + n1[5]*v[5*m]);
          u[p*incr] += n23*t1;
          v++;
        }
        v += 5*m;
      }
    }
  }
}

void TACSQuinticHexaBasis::addInterpFieldsTranspose( const int n,
                                                     const double pt[],
                                                     const int incr,
                                                     const TacsScalar u[],
                                                     const int m,
                                                     TacsScalar v[] ){
  if (n >= 0){
    const double *n1 = &Nf[6*(n % 6)];
    const double *n2 = &Nf[6*((n % 36)/6)];
    const double *n3 = &Nf[6*(n / 36)];

    for ( int k = 0; k < 6; k++ ){
      for ( int j = 0; j < 6; j++ ){
        TacsScalar n23  = n2[j]*n3[k];
        for ( int p = 0; p < m; p++ ){
          TacsScalar a = n23*u[p*incr];
          v[0] +=   a*n1[0];
          v[m] +=   a*n1[1];
          v[2*m] += a*n1[2];
          v[3*m] += a*n1[3];
          v[4*m] += a*n1[4];
          v[5*m] += a*n1[5];
          v++;
        }
        v += 5*m;
      }
    }
  }
  else {
    double n1[6], n2[6], n3[6];
    TacsLagrangeShapeFunctions(6, pt[0], cosine_pts, n1);
    TacsLagrangeShapeFunctions(6, pt[1], cosine_pts, n2);
    TacsLagrangeShapeFunctions(6, pt[2], cosine_pts, n3);

    for ( int k = 0; k < 6; k++ ){
      for ( int j = 0; j < 6; j++ ){
        TacsScalar n23 = n2[j]*n3[k];
        for ( int p = 0; p < m; p++ ){
          TacsScalar a = n23*u[p*incr];
          v[0] +=   a*n1[0];
          v[m] +=   a*n1[1];
          v[2*m] += a*n1[2];
          v[3*m] += a*n1[3];
          v[4*m] += a*n1[4];
          v[5*m] += a*n1[5];
          v++;
        }
        v += 5*m;
      }
    }
  }
}

void TACSQuinticHexaBasis::interpFieldsGrad( const int n,
                                             const double pt[],
                                             const int m,
                                             const TacsScalar v[],
                                             TacsScalar grad[] ){
  if (n >= 0){
    const double *n1 = &Nf[6*(n % 6)];
    const double *n2 = &Nf[6*((n % 36)/6)];
    const double *n3 = &Nf[6*(n / 36)];
    const double *n1x = &Nfxi[6*(n % 6)];
    const double *n2x = &Nfxi[6*((n % 36)/6)];
    const double *n3x = &Nfxi[6*(n / 36)];

    memset(grad, 0, 3*m*sizeof(TacsScalar));

    for ( int k = 0; k < 6; k++ ){
      for ( int j = 0; j < 6; j++ ){
        TacsScalar n23  = n2[j]*n3[k];
        TacsScalar n2x3 = n2x[j]*n3[k];
        TacsScalar n23x = n2[j]*n3x[k];

        TacsScalar *g = grad;
        for ( int p = 0; p < m; p++ ){
          g[0] += n23*(n1x[0]*v[0]   + n1x[1]*v[m]   + n1x[2]*v[2*m] +
                       n1x[3]*v[3*m] + n1x[4]*v[4*m] + n1x[5]*v[5*m]);
          TacsScalar t1 = (n1[0]*v[0]   + n1[1]*v[m]   + n1[2]*v[2*m] +
                           n1[3]*v[3*m] + n1[4]*v[4*m] + n1[5]*v[5*m]);
          g[1] += n2x3*t1;
          g[2] += n23x*t1;
          g += 3;
          v++;
        }
        v += 5*m;
      }
    }
  }
  else {
    double n1[6], n2[6], n3[6], n1x[6], n2x[6], n3x[6];
    TacsLagrangeShapeFuncDerivative(6, pt[0], cosine_pts, n1, n1x);
    TacsLagrangeShapeFuncDerivative(6, pt[1], cosine_pts, n2, n2x);
    TacsLagrangeShapeFuncDerivative(6, pt[2], cosine_pts, n3, n3x);

    for ( int k = 0; k < 6; k++ ){
      for ( int j = 0; j < 6; j++ ){
        TacsScalar n23  = n2[j]*n3[k];
        TacsScalar n2x3 = n2x[j]*n3[k];
        TacsScalar n23x = n2[j]*n3x[k];

        TacsScalar *g = grad;
        for ( int p = 0; p < m; p++ ){
          g[0] += n23*(n1x[0]*v[0]   + n1x[1]*v[m]   + n1x[2]*v[2*m] +
                       n1x[3]*v[3*m] + n1x[4]*v[4*m] + n1x[5]*v[5*m]);
          TacsScalar t1 = (n1[0]*v[0]   + n1[1]*v[m]   + n1[2]*v[2*m] +
                           n1[3]*v[3*m] + n1[4]*v[4*m] + n1[5]*v[5*m]);
          g[1] += n2x3*t1;
          g[2] += n23x*t1;
          g += 3;
          v++;
        }
        v += 5*m;
      }
    }
  }
}

void TACSQuinticHexaBasis::addInterpFieldsGradTranspose( int n,
                                                         const double pt[],
                                                         const int m,
                                                         const TacsScalar grad[],
                                                         TacsScalar v[] ){
  if (n >= 0){
    const double *n1 = &Nf[6*(n % 6)];
    const double *n2 = &Nf[6*((n % 36)/6)];
    const double *n3 = &Nf[6*(n / 36)];
    const double *n1x = &Nfxi[6*(n % 6)];
    const double *n2x = &Nfxi[6*((n % 36)/6)];
    const double *n3x = &Nfxi[6*(n / 36)];

    for ( int k = 0; k < 6; k++ ){
      for ( int j = 0; j < 6; j++ ){
        TacsScalar n23  = n2[j]*n3[k];
        TacsScalar n2x3 = n2x[j]*n3[k];
        TacsScalar n23x = n2[j]*n3x[k];

        const TacsScalar *g = grad;
        for ( int p = 0; p < m; p++ ){
          TacsScalar a = n2x3*g[1] + n23x*g[2];
          TacsScalar b = n23*g[0];
          v[0] +=   a*n1[0] + b*n1x[0];
          v[m] +=   a*n1[1] + b*n1x[1];
          v[2*m] += a*n1[2] + b*n1x[2];
          v[3*m] += a*n1[3] + b*n1x[3];
          v[4*m] += a*n1[4] + b*n1x[4];
          v[5*m] += a*n1[5] + b*n1x[5];
          g += 3;
          v++;
        }
        v += 5*m;
      }
    }
  }
  else {
    double n1[6], n2[6], n3[6], n1x[6], n2x[6], n3x[6];
    TacsLagrangeShapeFuncDerivative(6, pt[0], cosine_pts, n1, n1x);
    TacsLagrangeShapeFuncDerivative(6, pt[1], cosine_pts, n2, n2x);
    TacsLagrangeShapeFuncDerivative(6, pt[2], cosine_pts, n3, n3x);

    for ( int k = 0; k < 6; k++ ){
      for ( int j = 0; j < 6; j++ ){
        TacsScalar n23  = n2[j]*n3[k];
        TacsScalar n2x3 = n2x[j]*n3[k];
        TacsScalar n23x = n2[j]*n3x[k];

        const TacsScalar *g = grad;
        for ( int p = 0; p < m; p++ ){
          TacsScalar a = n2x3*g[1] + n23x*g[2];
          TacsScalar b = n23*g[0];
          v[0] +=   a*n1[0] + b*n1x[0];
          v[m] +=   a*n1[1] + b*n1x[1];
          v[2*m] += a*n1[2] + b*n1x[2];
          v[3*m] += a*n1[3] + b*n1x[3];
          v[4*m] += a*n1[4] + b*n1x[4];
          v[5*m] += a*n1[5] + b*n1x[5];
          g += 3;
          v++;
        }
        v += 5*m;
      }
    }
  }
}

void TACSQuinticHexaBasis::interpAllFieldsGrad( const int m,
                                                const TacsScalar values[],
                                                TacsScalar out[] ){
  memset(out, 0, 4*m*216*sizeof(TacsScalar));

  for ( int n = 0; n < 216; n++ ){
    const double *n1 = &Nf[6*(n % 6)];
    const double *n2 = &Nf[6*((n % 36)/6)];
    const double *n3 = &Nf[6*(n / 36)];
    const double *n1x = &Nfxi[6*(n % 6)];
    const double *n2x = &Nfxi[6*((n % 36)/6)];
    const double *n3x = &Nfxi[6*(n / 36)];

    const TacsScalar *v = values;

    if (m == 1){
      for ( int k = 0; k < 6; k++ ){
        for ( int j = 0; j < 6; j++ ){
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          out[1] += n23*(n1x[0]*v[0] + n1x[1]*v[1] + n1x[2]*v[2] +
                         n1x[3]*v[3] + n1x[4]*v[4] + n1x[5]*v[5]);
          TacsScalar t1 = (n1[0]*v[0] + n1[1]*v[1] + n1[2]*v[2] +
                           n1[3]*v[3] + n1[4]*v[4] + n1[5]*v[5]);
          out[0] += n23*t1;
          out[2] += n2x3*t1;
          out[3] += n23x*t1;
          v += 6;
        }
      }
    }
    else if (m == 3){
      for ( int k = 0; k < 6; k++ ){
        for ( int j = 0; j < 6; j++ ){
          TacsScalar t1;
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          out[3] += n23*(n1x[0]*v[0] + n1x[1]*v[3]  + n1x[2]*v[6] +
                         n1x[3]*v[9] + n1x[4]*v[12] + n1x[5]*v[15]);
          t1 = (n1[0]*v[0] + n1[1]*v[3]  + n1[2]*v[6] +
                n1[3]*v[9] + n1[4]*v[12] + n1[5]*v[15]);
          out[0] += n23*t1;
          out[4] += n2x3*t1;
          out[5] += n23x*t1;

          out[6] += n23*(n1x[0]*v[1]  + n1x[1]*v[4]  + n1x[2]*v[7] +
                         n1x[3]*v[10] + n1x[4]*v[13] + n1x[5]*v[16]);
          t1 = (n1[0]*v[1]  + n1[1]*v[4]  + n1[2]*v[7] +
                n1[3]*v[10] + n1[4]*v[13] + n1[5]*v[16]);
          out[1] += n23*t1;
          out[7] += n2x3*t1;
          out[8] += n23x*t1;

          out[9] += n23*(n1x[0]*v[2]  + n1x[1]*v[5]  + n1x[2]*v[8] +
                         n1x[3]*v[11] + n1x[4]*v[14] + n1x[5]*v[17]);
          t1 = (n1[0]*v[2]  + n1[1]*v[5]   + n1[2]*v[8] +
                n1[3]*v[11] + n1[4]*v[14] + n1[5]*v[17]);
          out[2] += n23*t1;
          out[10] += n2x3*t1;
          out[11] += n23x*t1;

          v += 18;
        }
      }
    }
    else if (m == 4){
      for ( int k = 0; k < 6; k++ ){
        for ( int j = 0; j < 6; j++ ){
          TacsScalar t1;
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          out[4] += n23*(n1x[0]*v[0]  + n1x[1]*v[4]  + n1x[2]*v[8] +
                         n1x[3]*v[12] + n1x[4]*v[16] + n1x[5]*v[20]);
          t1 = (n1[0]*v[0] + n1[1]*v[4]  + n1[2]*v[8] +
                n1[3]*v[12] + n1[4]*v[16] + n1[5]*v[20]);
          out[0] += n23*t1;
          out[5] += n2x3*t1;
          out[6] += n23x*t1;

          out[7] += n23*(n1x[0]*v[1]  + n1x[1]*v[5]  + n1x[2]*v[9] +
                         n1x[3]*v[13] + n1x[4]*v[17] + n1x[5]*v[21]);
          t1 = (n1[0]*v[1]  + n1[1]*v[5]  + n1[2]*v[9] +
                n1[3]*v[13] + n1[4]*v[17] + n1[5]*v[21]);
          out[1] += n23*t1;
          out[8] += n2x3*t1;
          out[9] += n23x*t1;

          out[10] += n23*(n1x[0]*v[2]  + n1x[1]*v[6]  + n1x[2]*v[10] +
                          n1x[3]*v[14] + n1x[4]*v[18] + n1x[5]*v[22]);
          t1 = (n1[0]*v[2]  + n1[1]*v[6]  + n1[2]*v[10] +
                n1[3]*v[14] + n1[4]*v[18] + n1[5]*v[22]);
          out[2] += n23*t1;
          out[11] += n2x3*t1;
          out[12] += n23x*t1;

          out[13] += n23*(n1x[0]*v[3]  + n1x[1]*v[7]  + n1x[2]*v[11] +
                         n1x[3]*v[15] + n1x[4]*v[19] + n1x[5]*v[23]);
          t1 = (n1[0]*v[3] + n1[1]*v[7]  + n1[2]*v[11] +
                n1[3]*v[15] + n1[4]*v[19] + n1[5]*v[23]);
          out[3] += n23*t1;
          out[14] += n2x3*t1;
          out[15] += n23x*t1;

          v += 24;
        }
      }
    }
    else {
      for ( int k = 0; k < 6; k++ ){
        for ( int j = 0; j < 6; j++ ){
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          TacsScalar *g = &out[m];
          for ( int p = 0; p < m; p++ ){
            g[0] += n23*(n1x[0]*v[0]   + n1x[1]*v[m]   + n1x[2]*v[2*m] +
                         n1x[3]*v[3*m] + n1x[4]*v[4*m] + n1x[5]*v[5*m]);
            TacsScalar t1 = (n1[0]*v[0]   + n1[1]*v[m]   + n1[2]*v[2*m] +
                             n1[3]*v[3*m] + n1[4]*v[4*m] + n1[5]*v[5*m]);
            out[p] += n23*t1;
            g[1] += n2x3*t1;
            g[2] += n23x*t1;
            g += 3;
            v++;
          }
          v += 5*m;
        }
      }
    }

    out += 4*m;
  }
}

void TACSQuinticHexaBasis::addInterpAllFieldsGradTranspose( const int m,
                                                            const TacsScalar in[],
                                                            TacsScalar values[] ){
  for ( int n = 0; n < 216; n++ ){
    const double *n1 = &Nf[6*(n % 6)];
    const double *n2 = &Nf[6*((n % 36)/6)];
    const double *n3 = &Nf[6*(n / 36)];
    const double *n1x = &Nfxi[6*(n % 6)];
    const double *n2x = &Nfxi[6*((n % 36)/6)];
    const double *n3x = &Nfxi[6*(n / 36)];

    TacsScalar *v = values;

    if (m == 1){
      for ( int k = 0; k < 6; k++ ){
        for ( int j = 0; j < 6; j++ ){
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          TacsScalar a = n23*in[0] + n2x3*in[2] + n23x*in[3];
          TacsScalar b = n23*in[1];
          v[0] += a*n1[0] + b*n1x[0];
          v[1] += a*n1[1] + b*n1x[1];
          v[2] += a*n1[2] + b*n1x[2];
          v[3] += a*n1[3] + b*n1x[3];
          v[4] += a*n1[4] + b*n1x[4];
          v[5] += a*n1[5] + b*n1x[5];
          v += 6;
        }
      }
    }
    else if (m == 3){
      for ( int k = 0; k < 6; k++ ){
        for ( int j = 0; j < 6; j++ ){
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          TacsScalar a1 = n23*in[0] + n2x3*in[4] + n23x*in[5];
          TacsScalar b1 = n23*in[3];

          TacsScalar a2 = n23*in[1] + n2x3*in[7] + n23x*in[8];
          TacsScalar b2 = n23*in[6];

          TacsScalar a3 = n23*in[2] + n2x3*in[10] + n23x*in[11];
          TacsScalar b3 = n23*in[9];

          v[0] +=  a1*n1[0] + b1*n1x[0];
          v[1] +=  a2*n1[0] + b2*n1x[0];
          v[2] +=  a3*n1[0] + b3*n1x[0];

          v[3] +=  a1*n1[1] + b1*n1x[1];
          v[4] +=  a2*n1[1] + b2*n1x[1];
          v[5] +=  a3*n1[1] + b3*n1x[1];

          v[6] +=  a1*n1[2] + b1*n1x[2];
          v[7] +=  a2*n1[2] + b2*n1x[2];
          v[8] +=  a3*n1[2] + b3*n1x[2];

          v[9] +=  a1*n1[3] + b1*n1x[3];
          v[10] += a2*n1[3] + b2*n1x[3];
          v[11] += a3*n1[3] + b3*n1x[3];

          v[12] += a1*n1[4] + b1*n1x[4];
          v[13] += a2*n1[4] + b2*n1x[4];
          v[14] += a3*n1[4] + b3*n1x[4];

          v[15] += a1*n1[5] + b1*n1x[5];
          v[16] += a2*n1[5] + b2*n1x[5];
          v[17] += a3*n1[5] + b3*n1x[5];

          v += 18;
        }
      }
    }
    else if (m == 4){
      for ( int k = 0; k < 6; k++ ){
        for ( int j = 0; j < 6; j++ ){
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          TacsScalar a1 = n23*in[0] + n2x3*in[5] + n23x*in[6];
          TacsScalar b1 = n23*in[4];

          TacsScalar a2 = n23*in[1] + n2x3*in[8] + n23x*in[9];
          TacsScalar b2 = n23*in[7];

          TacsScalar a3 = n23*in[2] + n2x3*in[11] + n23x*in[12];
          TacsScalar b3 = n23*in[10];

          TacsScalar a4 = n23*in[3] + n2x3*in[14] + n23x*in[15];
          TacsScalar b4 = n23*in[13];

          v[0] +=  a1*n1[0] + b1*n1x[0];
          v[1] +=  a2*n1[0] + b2*n1x[0];
          v[2] +=  a3*n1[0] + b3*n1x[0];
          v[3] +=  a4*n1[0] + b4*n1x[0];

          v[4] +=  a1*n1[1] + b1*n1x[1];
          v[5] +=  a2*n1[1] + b2*n1x[1];
          v[6] +=  a3*n1[1] + b3*n1x[1];
          v[7] +=  a4*n1[1] + b4*n1x[1];

          v[8] +=  a1*n1[2] + b1*n1x[2];
          v[9] +=  a2*n1[2] + b2*n1x[2];
          v[10] += a3*n1[2] + b3*n1x[2];
          v[11] += a4*n1[2] + b4*n1x[2];

          v[12] += a1*n1[3] + b1*n1x[3];
          v[13] += a2*n1[3] + b2*n1x[3];
          v[14] += a3*n1[3] + b3*n1x[3];
          v[15] += a4*n1[3] + b4*n1x[3];

          v[16] += a1*n1[4] + b1*n1x[4];
          v[17] += a2*n1[4] + b2*n1x[4];
          v[18] += a3*n1[4] + b3*n1x[4];
          v[19] += a4*n1[4] + b4*n1x[4];

          v[20] += a1*n1[5] + b1*n1x[5];
          v[21] += a2*n1[5] + b2*n1x[5];
          v[22] += a3*n1[5] + b3*n1x[5];
          v[23] += a4*n1[5] + b4*n1x[5];

          v += 24;
        }
      }
    }
    else {
      for ( int k = 0; k < 6; k++ ){
        for ( int j = 0; j < 6; j++ ){
          TacsScalar n23  = n2[j]*n3[k];
          TacsScalar n2x3 = n2x[j]*n3[k];
          TacsScalar n23x = n2[j]*n3x[k];

          const TacsScalar *g = &in[m];
          for ( int p = 0; p < m; p++ ){
            TacsScalar a = n23*in[p] + n2x3*g[1] + n23x*g[2];
            TacsScalar b = n23*g[0];
            v[0] +=   a*n1[0] + b*n1x[0];
            v[m] +=   a*n1[1] + b*n1x[1];
            v[2*m] += a*n1[2] + b*n1x[2];
            v[3*m] += a*n1[3] + b*n1x[3];
            v[4*m] += a*n1[4] + b*n1x[4];
            v[5*m] += a*n1[5] + b*n1x[5];
            g += 3;
            v++;
          }
          v += 5*m;
        }
      }
    }

    in += 4*m;
  }
}
