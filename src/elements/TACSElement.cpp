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

#include <stdlib.h>
#include "TACSElement.h"
#include "tacslapack.h"

/*
  Finds the finite-difference based Jacobian of the element. This is
  the default Jacobian implementation for any TACSElement. The user
  can override this function and provide an analytic Jacobian
  implemention in descendant classes.
*/
void TACSElement::addJacobian( int elemIndex,
                               double time,
                               TacsScalar alpha,
                               TacsScalar beta,
                               TacsScalar gamma,
                               const TacsScalar Xpts[],
                               const TacsScalar vars[],
                               const TacsScalar dvars[],
                               const TacsScalar ddvars[],
                               TacsScalar res[],
                               TacsScalar J[] ){
  // Call the residual implementation
  addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);

  // Get the number of variables
  int nvars = getNumVariables();

  // The step length
#ifdef TACS_USE_COMPLEX
  const double dh = 1e-30;
#else
  const double dh = 1e-7;
#endif // TACS_USE_COMPLEX

  // Original and perturbed residual vectors
  TacsScalar *Rtmp1 = new TacsScalar[nvars];
  TacsScalar *Rtmp2 = new TacsScalar[nvars];

  // Perturbed state vectors
  // TacsScalar *pstate = new TacsScalar[nvars];
  TacsScalar *qTmp = new TacsScalar[nvars];
  TacsScalar *qdotTmp = new TacsScalar[nvars];
  TacsScalar *qddotTmp = new TacsScalar[nvars];

  // Copy the state variables into pstate
  memcpy(qTmp, vars, nvars*sizeof(TacsScalar));
  memcpy(qdotTmp, dvars, nvars*sizeof(TacsScalar));
  memcpy(qddotTmp, ddvars, nvars*sizeof(TacsScalar));

  // Perturb each state variable and find the residual
  for ( int i = 0; i < nvars; i++ ){
#ifdef TACS_USE_COMPLEX
    qTmp[i] = vars[i] + TacsScalar(0.0, dh);

    // Assemble the unperturbed residual
    memset(Rtmp1, 0, nvars*sizeof(TacsScalar));
    addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp1);

    // Find the approximated jacobian
    for ( int j = 0; j < nvars; j++ ){
      J[j*nvars+i] += alpha*TacsImagPart(Rtmp1[j])/dh;
    }
#else
    // Perturb the i-th variable
    qTmp[i] = vars[i] + dh;

    // Assemble the unperturbed residual
    memset(Rtmp1, 0, nvars*sizeof(TacsScalar));
    addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp1);

    // Perturb the i-th variable
    qTmp[i] = vars[i] - dh;

    // Assemble the unperturbed residual
    memset(Rtmp2, 0, nvars*sizeof(TacsScalar));
    addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp2);

    // Find the approximated jacobian
    for ( int j = 0; j < nvars; j++ ){
      J[j*nvars+i] += 0.5*alpha*(Rtmp1[j] - Rtmp2[j])/dh;
    }
#endif // TACS_USE_COMPLEX
    // Restore the i-th variable
    qTmp[i] = vars[i];
  }

  // Perturb each state variable and find the residual
  for ( int i = 0; i < nvars; i++ ){
#ifdef TACS_USE_COMPLEX
    qdotTmp[i] = dvars[i] + TacsScalar(0.0, dh);

    // Assemble the unperturbed residual
    memset(Rtmp1, 0, nvars*sizeof(TacsScalar));
    addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp1);

    // Find the approximated jacobian
    for ( int j = 0; j < nvars; j++ ){
      J[j*nvars+i] += beta*TacsImagPart(Rtmp1[j])/dh;
    }
#else
    // Perturb the i-th variable
    qdotTmp[i] = dvars[i] + dh;

    // Assemble the unperturbed residual
    memset(Rtmp1, 0, nvars*sizeof(TacsScalar));
    addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp1);

    // Perturb the i-th variable
    qdotTmp[i] = dvars[i] - dh;

    // Assemble the unperturbed residual
    memset(Rtmp2, 0, nvars*sizeof(TacsScalar));
    addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp2);

    // Find the approximated jacobian
    for ( int j = 0; j < nvars; j++ ){
      J[j*nvars+i] += 0.5*beta*(Rtmp1[j] - Rtmp2[j])/dh;
    }
#endif // TACS_USE_COMPLEX
    // Restore the i-th variable
    qdotTmp[i] = dvars[i];
  }

  // Perturb each state variable and find the residual
  for ( int i = 0; i < nvars; i++ ){
#ifdef TACS_USE_COMPLEX
    qddotTmp[i] = ddvars[i] + TacsScalar(0.0, dh);

    // Assemble the unperturbed residual
    memset(Rtmp1, 0, nvars*sizeof(TacsScalar));
    addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp1);

    // Find the approximated jacobian
    for ( int j = 0; j < nvars; j++ ){
      J[j*nvars+i] += gamma*TacsImagPart(Rtmp1[j])/dh;
    }
#else
    // Perturb the i-th variable
    qddotTmp[i] = ddvars[i] + dh;

    // Assemble the unperturbed residual
    memset(Rtmp1, 0, nvars*sizeof(TacsScalar));
    addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp1);

    // Perturb the i-th variable
    qddotTmp[i] = ddvars[i] - dh;

    // Assemble the unperturbed residual
    memset(Rtmp2, 0, nvars*sizeof(TacsScalar));
    addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp2);

    // Find the approximated jacobian
    for ( int j = 0; j < nvars; j++ ){
      J[j*nvars+i] += 0.5*gamma*(Rtmp1[j] - Rtmp2[j])/dh;
    }
#endif // TACS_USE_COMPLEX
    // Restore the i-th variable
    qddotTmp[i] = ddvars[i];
  }

  delete [] Rtmp1;
  delete [] Rtmp2;
  delete [] qTmp;
  delete [] qdotTmp;
  delete [] qddotTmp;
}
