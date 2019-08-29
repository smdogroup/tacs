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


/*
  Linear Tetrahedral basis class functions
*/
int TACSLinearTetrahedralBasis::getNumNodes(){
  return 3;
}

int TACSLinearTetrahedralBasis::getNumParameters(){
  return 2;
}

int TACSLinearTetrahedralBasis::getNumQuadraturePoints(){
  return 1;
}

double TACSLinearTetrahedralBasis::getQuadratureWeight( int n ){
  return 1.0;
}

double TACSLinearTetrahedralBasis::getQuadraturePoint( int n,
                                                       double pt[] ){
  return 1.0;
}

int TACSLinearTetrahedralBasis::getNumElementFaces(){
  return 3;
}

int TACSLinearTetrahedralBasis::getNumFaceQuadraturePoints( int face ){
  return 1;
}

double TACSLinearTetrahedralBasis::getFaceQuadraturePoint( int face,
                                                           int n,
                                                           double pt[] ){

  return 0.1;
}

void TACSLinearTetrahedralBasis::computeBasis( const double pt[],
                                               double N[] ){


}

void TACSLinearTetrahedralBasis::computeBasisGradient( const double pt[],
                                                       double N[],
                                                       double Nxi[] ){


}

/*
  Quadratic Tetrahedral basis class functions
*/
int TACSQuadraticTetrahedralBasis::getNumNodes(){
  return 3;
}

int TACSQuadraticTetrahedralBasis::getNumParameters(){
  return 2;
}

int TACSQuadraticTetrahedralBasis::getNumQuadraturePoints(){
  return 1;
}

double TACSQuadraticTetrahedralBasis::getQuadratureWeight( int n ){
  return 1.0;
}

double TACSQuadraticTetrahedralBasis::getQuadraturePoint( int n,
                                                          double pt[] ){
  return 1.0;
}

int TACSQuadraticTetrahedralBasis::getNumElementFaces(){
  return 3;
}

int TACSQuadraticTetrahedralBasis::getNumFaceQuadraturePoints( int face ){
  return 1;
}

double TACSQuadraticTetrahedralBasis::getFaceQuadraturePoint( int face,
                                                              int n,
                                                              double pt[] ){

  return 0.1;
}

void TACSQuadraticTetrahedralBasis::computeBasis( const double pt[],
                                                  double N[] ){


}

void TACSQuadraticTetrahedralBasis::computeBasisGradient( const double pt[],
                                                          double N[],
                                                          double Nxi[] ){


}

/*
  Cubic Tetrahedral basis class functions
*/
int TACSCubicTetrahedralBasis::getNumNodes(){
  return 3;
}

int TACSCubicTetrahedralBasis::getNumParameters(){
  return 2;
}

int TACSCubicTetrahedralBasis::getNumQuadraturePoints(){
  return 1;
}

double TACSCubicTetrahedralBasis::getQuadratureWeight( int n ){
  return 1.0;
}

double TACSCubicTetrahedralBasis::getQuadraturePoint( int n,
                                                      double pt[] ){
  return 1.0;
}

int TACSCubicTetrahedralBasis::getNumElementFaces(){
  return 3;
}

int TACSCubicTetrahedralBasis::getNumFaceQuadraturePoints( int face ){
  return 1;
}

double TACSCubicTetrahedralBasis::getFaceQuadraturePoint( int face,
                                                          int n,
                                                          double pt[] ){

  return 0.1;
}

void TACSCubicTetrahedralBasis::computeBasis( const double pt[],
                                              double N[] ){


}

void TACSCubicTetrahedralBasis::computeBasisGradient( const double pt[],
                                                      double N[],
                                                      double Nxi[] ){


}
