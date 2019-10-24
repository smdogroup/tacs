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
    pt[1] = -1.0 + 2.0*(face % 2);
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
    pt[1] = -1.0 + 2.0*(face % 2);
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
  return 4;
}

int TACSCubicHexaBasis::getNumFaceQuadraturePoints( int face ){
  return 4;
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
    pt[1] = -1.0 + 2.0*(face % 2);
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
