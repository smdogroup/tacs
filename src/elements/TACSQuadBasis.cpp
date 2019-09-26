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

#include "TACSQuadBasis.h"
#include "TACSGaussQuadrature.h"

static void getEdgeTangent( int edge, double t[] ){
  if (edge == 0){
    // -X edge
    t[0] = 0.0;  t[1] = -1.0;
  }
  else if (edge == 1){
    // +X edge
    t[0] = 0.0;  t[1] = 1.0;
  }
  else if (edge == 2){
    // -Y edge
    t[0] = 1.0;  t[1] = 0.0;
  }
  else if (edge == 3){
    // +Y edge
    t[0] = -1.0;  t[1] = 0.0;
  }
}

/*
  Linear Quad basis class functions
*/
ElementLayout TACSLinearQuadBasis::getLayoutType(){
  return TACS_QUAD_ELEMENT;
}

void TACSLinearQuadBasis::getVisPoint( int n, double pt[] ){
  pt[0] = -1.0 + 2.0*(n % 2);
  pt[1] = -1.0 + 2.0*(n/2);
}

int TACSLinearQuadBasis::getNumNodes(){
  return 4;
}

int TACSLinearQuadBasis::getNumParameters(){
  return 2;
}

int TACSLinearQuadBasis::getNumQuadraturePoints(){
  return 4;
}

double TACSLinearQuadBasis::getQuadratureWeight( int n ){
  return TacsGaussQuadWts2[n % 2]*TacsGaussQuadWts2[n/2];
}

double TACSLinearQuadBasis::getQuadraturePoint( int n,
                                                double pt[] ){
  pt[0] = TacsGaussQuadPts2[n % 2];
  pt[1] = TacsGaussQuadPts2[n/2];

  return TacsGaussQuadWts2[n % 2]*TacsGaussQuadWts2[n/2];
}

int TACSLinearQuadBasis::getNumElementFaces(){
  return 4;
}

int TACSLinearQuadBasis::getNumFaceQuadraturePoints( int face ){
  return 2;
}

double TACSLinearQuadBasis::getFaceQuadraturePoint( int face,
                                                    int n,
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

  getEdgeTangent(face, t);

  return TacsGaussQuadWts2[n];
}

void TACSLinearQuadBasis::computeBasis( const double pt[],
                                        double N[] ){
  N[0] = 0.25*(1.0 - pt[0])*(1.0 - pt[1]);
  N[1] = 0.25*(1.0 + pt[0])*(1.0 - pt[1]);
  N[2] = 0.25*(1.0 - pt[0])*(1.0 + pt[1]);
  N[3] = 0.25*(1.0 + pt[0])*(1.0 + pt[1]);
}

void TACSLinearQuadBasis::computeBasisGradient( const double pt[],
                                                double N[],
                                                double Nxi[] ){
  N[0] = 0.25*(1.0 - pt[0])*(1.0 - pt[1]);
  N[1] = 0.25*(1.0 + pt[0])*(1.0 - pt[1]);
  N[2] = 0.25*(1.0 - pt[0])*(1.0 + pt[1]);
  N[3] = 0.25*(1.0 + pt[0])*(1.0 + pt[1]);

  Nxi[0] = -0.25*(1.0 - pt[1]);
  Nxi[1] = -0.25*(1.0 - pt[0]);

  Nxi[2] = 0.25*(1.0 - pt[1]);
  Nxi[3] = -0.25*(1.0 + pt[0]);

  Nxi[4] = -0.25*(1.0 + pt[1]);
  Nxi[5] = 0.25*(1.0 - pt[0]);

  Nxi[6] = 0.25*(1.0 + pt[1]);
  Nxi[7] = 0.25*(1.0 + pt[0]);
}

/*
  Quadratic Quad basis class functions
*/
ElementLayout TACSQuadraticQuadBasis::getLayoutType(){
  return TACS_QUAD_QUADRATIC_ELEMENT;
}

void TACSQuadraticQuadBasis::getVisPoint( int n, double pt[] ){
  pt[0] = -1.0 + 1.0*(n % 3);
  pt[1] = -1.0 + 1.0*(n/3);
}

int TACSQuadraticQuadBasis::getNumNodes(){
  return 9;
}

int TACSQuadraticQuadBasis::getNumParameters(){
  return 2;
}

int TACSQuadraticQuadBasis::getNumQuadraturePoints(){
  return 9;
}

double TACSQuadraticQuadBasis::getQuadratureWeight( int n ){
  return TacsGaussQuadWts3[n % 3]*TacsGaussQuadWts3[n/3];
}

double TACSQuadraticQuadBasis::getQuadraturePoint( int n,
                                                   double pt[] ){
  pt[0] = TacsGaussQuadPts3[n % 3];
  pt[1] = TacsGaussQuadPts3[n/3];

  return TacsGaussQuadWts3[n % 3]*TacsGaussQuadWts3[n/3];
}

int TACSQuadraticQuadBasis::getNumElementFaces(){
  return 4;
}

int TACSQuadraticQuadBasis::getNumFaceQuadraturePoints( int face ){
  return 3;
}

double TACSQuadraticQuadBasis::getFaceQuadraturePoint( int face,
                                                       int n,
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

  getEdgeTangent(face, t);

  return TacsGaussQuadWts3[n];
}

void TACSQuadraticQuadBasis::computeBasis( const double pt[],
                                           double N[] ){
  double na[3];
  na[0] = -0.5*pt[0]*(1.0 - pt[0]);
  na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
  na[2] = 0.5*(1.0 + pt[0])*pt[0];

  double nb[3];
  nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
  nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
  nb[2] = 0.5*(1.0 + pt[1])*pt[1];

  N[0] = na[0]*nb[0];
  N[1] = na[1]*nb[0];
  N[2] = na[2]*nb[0];
  N[3] = na[0]*nb[1];
  N[4] = na[1]*nb[1];
  N[5] = na[2]*nb[1];
  N[6] = na[0]*nb[2];
  N[7] = na[1]*nb[2];
  N[8] = na[2]*nb[2];
}

void TACSQuadraticQuadBasis::computeBasisGradient( const double pt[],
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

  N[0] = na[0]*nb[0];
  N[1] = na[1]*nb[0];
  N[2] = na[2]*nb[0];
  N[3] = na[0]*nb[1];
  N[4] = na[1]*nb[1];
  N[5] = na[2]*nb[1];
  N[6] = na[0]*nb[2];
  N[7] = na[1]*nb[2];
  N[8] = na[2]*nb[2];

  double dna[3];
  dna[0] = -0.5 + pt[0];
  dna[1] = -2.0*pt[0];
  dna[2] = 0.5 + pt[0];

  double dnb[3];
  dnb[0] = -0.5 + pt[1];
  dnb[1] = -2.0*pt[1];
  dnb[2] = 0.5 + pt[1];

  Nxi[0] = dna[0]*nb[0];
  Nxi[1] = na[0]*dnb[0];
  Nxi[2] = dna[1]*nb[0];
  Nxi[3] = na[1]*dnb[0];
  Nxi[4] = dna[2]*nb[0];
  Nxi[5] = na[2]*dnb[0];
  Nxi[6] = dna[0]*nb[1];
  Nxi[7] = na[0]*dnb[1];
  Nxi[8] = dna[1]*nb[1];
  Nxi[9] = na[1]*dnb[1];
  Nxi[10] = dna[2]*nb[1];
  Nxi[11] = na[2]*dnb[1];
  Nxi[12] = dna[0]*nb[2];
  Nxi[13] = na[0]*dnb[2];
  Nxi[14] = dna[1]*nb[2];
  Nxi[15] = na[1]*dnb[2];
  Nxi[16] = dna[2]*nb[2];
  Nxi[17] = na[2]*dnb[2];
}

/*
  Cubic Quad basis class functions
*/
ElementLayout TACSCubicQuadBasis::getLayoutType(){
  return TACS_QUAD_CUBIC_ELEMENT;
}

void TACSCubicQuadBasis::getVisPoint( int n, double pt[] ){
  static const double pts[4] = {-1.0, -0.5, 0.5, 1.0};
  pt[0] = pts[n % 4];
  pt[1] = pts[n/4];
}

int TACSCubicQuadBasis::getNumNodes(){
  return 9;
}

int TACSCubicQuadBasis::getNumParameters(){
  return 2;
}

int TACSCubicQuadBasis::getNumQuadraturePoints(){
  return 9;
}

double TACSCubicQuadBasis::getQuadratureWeight( int n ){
  return TacsGaussQuadWts4[n % 4]*TacsGaussQuadWts4[n/4];
}

double TACSCubicQuadBasis::getQuadraturePoint( int n,
                                               double pt[] ){
  pt[0] = TacsGaussQuadPts4[n % 4];
  pt[1] = TacsGaussQuadPts4[n/4];

  return TacsGaussQuadWts4[n % 4]*TacsGaussQuadWts4[n/4];
}

int TACSCubicQuadBasis::getNumElementFaces(){
  return 4;
}

int TACSCubicQuadBasis::getNumFaceQuadraturePoints( int face ){
  return 4;
}

double TACSCubicQuadBasis::getFaceQuadraturePoint( int face,
                                                   int n,
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

  getEdgeTangent(face, t);

  return TacsGaussQuadWts4[n];
}

void TACSCubicQuadBasis::computeBasis( const double pt[],
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

  N[0] = na[0]*nb[0];
  N[1] = na[1]*nb[0];
  N[2] = na[2]*nb[0];
  N[3] = na[3]*nb[0];
  N[4] = na[0]*nb[1];
  N[5] = na[1]*nb[1];
  N[6] = na[2]*nb[1];
  N[7] = na[3]*nb[1];
  N[8] = na[0]*nb[2];
  N[9] = na[1]*nb[2];
  N[10] = na[2]*nb[2];
  N[11] = na[3]*nb[2];
  N[12] = na[0]*nb[3];
  N[13] = na[1]*nb[3];
  N[14] = na[2]*nb[3];
  N[15] = na[3]*nb[3];
}

void TACSCubicQuadBasis::computeBasisGradient( const double pt[],
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

  N[0] = na[0]*nb[0];
  N[1] = na[1]*nb[0];
  N[2] = na[2]*nb[0];
  N[3] = na[3]*nb[0];
  N[4] = na[0]*nb[1];
  N[5] = na[1]*nb[1];
  N[6] = na[2]*nb[1];
  N[7] = na[3]*nb[1];
  N[8] = na[0]*nb[2];
  N[9] = na[1]*nb[2];
  N[10] = na[2]*nb[2];
  N[11] = na[3]*nb[2];
  N[12] = na[0]*nb[3];
  N[13] = na[1]*nb[3];
  N[14] = na[2]*nb[3];
  N[15] = na[3]*nb[3];

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

  N[0] = dna[0]*nb[0];
  N[1] = na[0]*dnb[0];

  N[2] = dna[1]*nb[0];
  N[3] = na[1]*dnb[0];

  N[4] = dna[2]*nb[0];
  N[5] = na[2]*dnb[0];

  N[6] = dna[3]*nb[0];
  N[7] = na[3]*dnb[0];

  N[8] = dna[0]*nb[1];
  N[9] = na[0]*dnb[1];

  N[10] = dna[1]*nb[1];
  N[11] = na[1]*dnb[1];

  N[12] = dna[2]*nb[1];
  N[13] = na[2]*dnb[1];

  N[14] = dna[3]*nb[1];
  N[15] = na[3]*dnb[1];

  N[16] = dna[0]*nb[2];
  N[17] = na[0]*dnb[2];

  N[18] = dna[1]*nb[2];
  N[19] = na[1]*dnb[2];

  N[20] = dna[2]*nb[2];
  N[21] = na[2]*dnb[2];

  N[22] = dna[3]*nb[2];
  N[23] = na[3]*dnb[2];

  N[24] = dna[0]*nb[3];
  N[25] = na[0]*dnb[3];

  N[26] = dna[1]*nb[3];
  N[27] = na[1]*dnb[3];

  N[28] = dna[2]*nb[3];
  N[29] = na[2]*dnb[3];

  N[30] = dna[3]*nb[3];
  N[31] = na[3]*dnb[3];
}
