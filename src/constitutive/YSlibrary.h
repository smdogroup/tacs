#ifndef TACS_YS_LIBRARY_H
#define TACS_YS_LIBRARY_H

/*
  The following library implements several methods used in different
  failure criteria. All the code here is based on an assumption of
  linear material behaviour.

  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
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

