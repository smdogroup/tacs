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

#ifndef TACS_QUAD_BASIS_H
#define TACS_QUAD_BASIS_H

#include "TACSMixedInterpElementBasis.h"

class TACSLinearMITCQuadBasis : public TACSMixedInterpElementBasis {
 public:
  ElementLayout getLayoutType();
  void getVisPoint( int n, double pt[] );
  int getNumNodes();
  int getNumParameters();
  int getNumQuadraturePoints();
  double getQuadratureWeight( int n );
  double getQuadraturePoint( int n, double pt[] );
  int getNumElementFaces();
  int getNumFaceQuadraturePoints( int face );

  int getNumTyingPoints(){
    return 4;
  }
  int getNumTyingFieldValues(){
    return 2;
  }
  void getTyingPoint( int ty, double pt ){
    if (ty == 0){

    }
    else if (ty == 1){

    }
    else if (ty == 2){

    }
    else if (ty == 3){

    }
  }
  void getTyingFieldValue( int n, const double pt[],
                           const TacsScalar qtys[], TacsScalar Uty ){
    
  }

  double getFaceQuadraturePoint( int face, int n, double pt[], double t[] );
  void computeBasis( const double pt[], double N[] );
  void computeBasisGradient( const double pt[], double N[], double Nxi[] );
};
