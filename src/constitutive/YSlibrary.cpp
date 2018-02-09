/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0 
*/

#include "YSlibrary.h"

/*
  Code for predicting the yield stress using various failure models
*/

void TestVonMises3D( double tol, double dh, 
		     TacsScalar stress[], 
		     TacsScalar ys ){
  
  printf("Testing the sensitivities of the 3D vonMises \
yield-stress calculations\n");
  
  TacsScalar sens[6];

  VonMisesFailure3DStressSens(sens, stress, ys);

  for ( int i = 0; i < 6; i++ ){
    stress[i] = stress[i] + dh;
    TacsScalar pf = VonMisesFailure3D(stress, ys);
    stress[i] = stress[i] - 2.0 * dh;
    TacsScalar pb = VonMisesFailure3D(stress, ys);
    stress[i] = stress[i] + dh;

    TacsScalar fd = 0.5*(pf - pb)/dh;
    printf("VonMises [%d] Value: %15.5f FD %15.5f Rel error: %10.3e \n", i, 
           TacsRealPart(sens[i]), TacsRealPart(fd), 
           TacsRealPart((sens[i] - fd)/fd));
  }
} 

void TestVonMisesPlaneStress( double tol, double dh, 
			      TacsScalar stress[], 
			      TacsScalar ys ){
  
  printf("Testing the sensitivities of the plane stress \
vonMises yield-stress calculations\n");
  
  TacsScalar sens[3];

  VonMisesFailurePlaneStressSens(sens, stress, ys);

  for ( int i = 0; i < 3; i++ ){
    stress[i] = stress[i] + dh;
    TacsScalar pf = VonMisesFailurePlaneStress(stress, ys);
    stress[i] = stress[i] - 2.0 * dh;
    TacsScalar pb = VonMisesFailurePlaneStress(stress, ys);
    stress[i] = stress[i] + dh;

    TacsScalar fd = 0.5*(pf - pb)/dh;

    printf("VonMises [%d] Value: %15.5f FD %15.5f Rel error: %10.3e \n", i, 
           TacsRealPart(sens[i]), TacsRealPart(fd), 
           TacsRealPart((sens[i] - fd)/fd));
  }
} 

/*!
  The von Mises failure criteria is the following:
  
  (sx - sy)^2 + (sx - sz)^2 + (sy - sz)^2 +
  6.0 * ( sxy^2 + sxz^2 + syz^2 ) = 2.0 * ys^2

  sx  = s[0]
  sy  = s[1]
  sz  = s[2]
  syz = s[3]
  sxz = s[4]
  sxy = s[5]
*/

TacsScalar VonMisesFailure3D( const TacsScalar s[], 
                              TacsScalar ys ){

  TacsScalar fail = sqrt(0.5*((s[0] - s[1])*(s[0] - s[1]) +
                              (s[0] - s[2])*(s[0] - s[2]) +
                              (s[1] - s[2])*(s[1] - s[2]) + 
                              6.0*(s[3]*s[3] + s[4]*s[4] + s[5]*s[5])))/ys;
  return fail;
}

TacsScalar VonMisesFailure3DStressSens( TacsScalar sens[], 
                                        const TacsScalar s[], 
                                        TacsScalar ys ){

  TacsScalar fail = sqrt(0.5*((s[0] - s[1])*(s[0] - s[1]) +
                              (s[0] - s[2])*(s[0] - s[2]) +
                              (s[1] - s[2])*(s[1] - s[2]) + 
                              6.0*(s[3]*s[3] + s[4]*s[4] + s[5]*s[5])));
  TacsScalar fact = 0.5/(ys*fail);

  sens[0] = fact*(2.0*s[0] - s[1] - s[2]);
  sens[1] = fact*(2.0*s[1] - s[0] - s[2]);
  sens[2] = fact*(2.0*s[2] - s[0] - s[1]);
  sens[3] = 6.0*fact*s[3];
  sens[4] = 6.0*fact*s[4];
  sens[5] = 6.0*fact*s[5];

  return fail;
}

/*
  Compute the von Mises failure criteria:

  (sx^2 + sy^2 - sx*sy + 3.0*sxy^2)/ys^2 < 1.0

  sx  = s[0]
  sy  = s[1]
  sxy = s[2]
*/

TacsScalar VonMisesFailurePlaneStress( const TacsScalar s[], TacsScalar ys ){
  TacsScalar fail = sqrt(s[0]*s[0] + s[1]*s[1] - s[0]*s[1] + 3.0*s[2]*s[2])/ys;
  return fail;
}

TacsScalar VonMisesFailurePlaneStressSens( TacsScalar sens[], 
					   const TacsScalar s[], 
					   TacsScalar ys ){
  TacsScalar fail = sqrt(s[0]*s[0] + s[1]*s[1] - s[0]*s[1] + 3.0*s[2]*s[2]);

  if (fail != 0.0){
    sens[0] = (s[0] - 0.5*s[1])/(fail*ys);
    sens[1] = (s[1] - 0.5*s[0])/(fail*ys);
    sens[2] = (3.0*s[2])/(fail*ys);
  }
  else {
    sens[0] = sens[1] = sens[2] = 0.0;
  }
  fail = fail/ys;

  return fail;
}

