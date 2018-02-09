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

#ifndef TACS_YS_LIBRARY_H
#define TACS_YS_LIBRARY_H

/*
  The following library implements several methods used in different
  failure criteria. All the code here is based on an assumption of
  linear material behaviour.
*/

#include <math.h>
#include <stdlib.h>
#include "TACSObject.h"

/*
  The von Mises failure criteria for three-dimensional stress problems
*/

TacsScalar VonMisesFailure3D( const TacsScalar stress[], 
                              TacsScalar ys );
TacsScalar VonMisesFailure3DStressSens( TacsScalar sens[], 
                                        const TacsScalar stress[], 
                                        TacsScalar ys );

/*
  The von Mises failure criteria for plane stress problems
*/

TacsScalar VonMisesFailurePlaneStress( const TacsScalar stress[], 
                                       TacsScalar ys );
TacsScalar VonMisesFailurePlaneStressSens( TacsScalar sens[], 
                                           const TacsScalar stress[], 
                                           TacsScalar ys );

void TestVonMises3D( double tol, double dh, 
		     TacsScalar stress[], 
		     TacsScalar ys );
void TestVonMisesPlaneStress( double tol, double dh, 
			      TacsScalar stress[], 
			      TacsScalar ys );

#endif

