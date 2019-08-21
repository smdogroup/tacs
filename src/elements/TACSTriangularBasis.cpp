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

#ifndef TACS_TRIANGULAR_BASIS_H
#define TACS_TRIANGULAR_BASIS_H

/*
  Linear Triangle basis class functions
*/
int TACSLinearTriangleBasis::getNumNodes(){
  return 3;
}

int TACSLinearTriangleBasis::getNumParameters(){
  return 2;
}

int TACSLinearTriangleBasis::getNumQuadraturePoints(){
  return 1;
}

double TACSLinearTriangleBasis::getQuadratureWeight( int n ){
  return 1.0;
}

double TACSLinearTriangleBasis::getQuadraturePoint( int n,
                                                    double pt[] ){
  return 1.0;
}

int TACSLinearTriangleBasis::getNumElementFaces(){
  return 3;
}

int TACSLinearTriangleBasis::getNumFaceQuadraturePoints( int face ){
  return 1;
}

double TACSLinearTriangleBasis::getFaceQuadraturePoint( int face,
                                                        int n,
                                                        double pt[] ){

  return 0.1;
}

void TACSLinearTriangleBasis::computeBasis( const double pt[],
                                            double N[] ){
  

}

void TACSLinearTriangleBasis::computeBasisGradient( const double pt[],
                                                    double N[],
                                                    double Nxi[] ){
  

}

/*
  Quadratic Triangle basis class functions
*/
int TACSQuadraticTriangleBasis::getNumNodes(){
  return 3;
}

int TACSQuadraticTriangleBasis::getNumParameters(){
  return 2;
}

int TACSQuadraticTriangleBasis::getNumQuadraturePoints(){
  return 1;
}

double TACSQuadraticTriangleBasis::getQuadratureWeight( int n ){
  return 1.0;
}

double TACSQuadraticTriangleBasis::getQuadraturePoint( int n,
                                                       double pt[] ){
  return 1.0;
}

int TACSQuadraticTriangleBasis::getNumElementFaces(){
  return 3;
}

int TACSQuadraticTriangleBasis::getNumFaceQuadraturePoints( int face ){
  return 1;
}

double TACSQuadraticTriangleBasis::getFaceQuadraturePoint( int face,
                                                        int n,
                                                        double pt[] ){

  return 0.1;
}

void TACSQuadraticTriangleBasis::computeBasis( const double pt[],
                                            double N[] ){
  

}

void TACSQuadraticTriangleBasis::computeBasisGradient( const double pt[],
                                                    double N[],
                                                    double Nxi[] ){
  

}

/*
  Cubic Triangle basis class functions
*/
int TACSCubicTriangleBasis::getNumNodes(){
  return 3;
}

int TACSCubicTriangleBasis::getNumParameters(){
  return 2;
}

int TACSCubicTriangleBasis::getNumQuadraturePoints(){
  return 1;
}

double TACSCubicTriangleBasis::getQuadratureWeight( int n ){
  return 1.0;
}

double TACSCubicTriangleBasis::getQuadraturePoint( int n,
                                                   double pt[] ){
  return 1.0;
}

int TACSCubicTriangleBasis::getNumElementFaces(){
  return 3;
}

int TACSCubicTriangleBasis::getNumFaceQuadraturePoints( int face ){
  return 1;
}

double TACSCubicTriangleBasis::getFaceQuadraturePoint( int face,
                                                       int n,
                                                       double pt[] ){

  return 0.1;
}

void TACSCubicTriangleBasis::computeBasis( const double pt[],
                                           double N[] ){
  

}

void TACSCubicTriangleBasis::computeBasisGradient( const double pt[],
                                                   double N[],
                                                   double Nxi[] ){
  

}
