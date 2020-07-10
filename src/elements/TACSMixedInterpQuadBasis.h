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

  int getNumTyingFields(){
    return 2;
  }
  int getNumTyingPoints( int field ){
    return 2;
  }
  void getTyingPoint( int field, int ty, double pt[] ){
    if (field == 0){
      if (ty == 0){
        pt[0] = -1.0;
        pt[1] = 0.0;
      }
      else if (ty == 1){
        pt[0] = 1.0;
        pt[1] = 0.0;
      }
    }
    else if (field == 1){
      if (ty == 0){
        pt[0] = 0.0;
        pt[1] = -1.0;
      }
      else if (ty == 1){
        pt[0] = 0.0;
        pt[1] = 1.0;
      }
    }
  }

  void getTyingFieldValue( int n, const double pt[],
                           const TacsScalar qtys[], TacsScalar Uty ){

  }

  double getFaceQuadraturePoint( int face, int n, double pt[], double t[] );
  void computeBasis( const double pt[], double N[] );
  void computeBasisGradient( const double pt[], double N[], double Nxi[] );
};
