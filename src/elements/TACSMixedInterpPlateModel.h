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

#include "TACSMixedInterpElementModel.h"

class TACSMixedInterpPlateModel : public TACSMixedInterpElementModel {
 public:
  TACSMixedInterpPlateModel();

  // Get the problem dimensions
  int getSpatialDim();
  int getVarsPerNode();
  int getDesignVarsPerNode();

  /**
    Retrieve the global design variable numbers associated with this element
  */
  int getDesignVarNums( int elemIndex, int dvLen, int dvNums[] );

  /**
    Set the element design variables from the design vector
  */
  int setDesignVars( int elemIndex, int dvLen, const TacsScalar dvs[] );

  /**
    Get the element design variables values
  */
  int getDesignVars( int elemIndex, int dvLen, TacsScalar dvs[] );

  /**
    Get the lower and upper bounds for the design variable values
  */
  int getDesignVarRange( int elemIndex, int dvLen,
                         TacsScalar lb[], TacsScalar ub[] );

  void evalPlateStrain( const TacsScalar Ut[],
                        const TacsScalar Ux[],
                        TacsScalar e[] ){
    //   0,  1,  2   3   4   5      6       7      8      9
    // u,x u,y v,x v,y w,x w,y rotx,x  rotx,y roty,x roty,y
    e[0] = Ux[0]; // exx = u,x
    e[1] = Ux[3]; // eyy = v,y
    e[2] = Ux[1] + Ux[2]; // exy = u,y + v,x

    e[3] =  Ux[8]; // kxx = roty,x
    e[4] = -Ux[7]; // kyy = -rotx,y
    e[5] = Ux[9] - Ux[6]; // kxy = roty,y - rotx,x

    e[8] = 0.0;
  }

  /*
    Evaluate the tying point quantity for the specified field value
    at the given tying point
  */
  TacsScalar evalTyingQuantity( int field, int typt, const double pt[],
                                const TacsScalar Xd[],
                                const TacsScalar U[], const TacsScalar Ud[] ){
    if (field == 0){
      return Ud[4]*Xd[2] + Ud[5]*Xd[3] - U[3]; // e23 = w,2 - rotx
    }
    else if (field == 1){
      return Ud[4]*Xd[0] + Ud[5]*Xd[1] + U[4]; // e13 = w,1 + roty
    }
  }

  void evalWeakIntegrand( int elemIndex,
                          const double time,
                          int n, const double pt[],
                          const TacsScalar X[],
                          const TacsScalar Xd[],
                          const TacsScalar Ut[],
                          const TacsScalar Ux[],
                          TacsScalar DUt[],
                          TacsScalar DUx[] ){

  }



};
