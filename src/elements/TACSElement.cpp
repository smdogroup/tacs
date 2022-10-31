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

#include "TACSElement.h"

#include <stdlib.h>

#include "tacslapack.h"

/*
  Default finite difference order for real analysis
*/
int TACSElement::fdOrder = 2;

/*
  Allow users to set default finite difference order for real analysis
*/
void TACSElement::setFiniteDifferenceOrder(int order) { fdOrder = order; }

/*
  Finds the finite-difference based Jacobian of the element. This is
  the default Jacobian implementation for any TACSElement. The user
  can override this function and provide an analytic Jacobian
  implemention in descendant classes.
*/
void TACSElement::addJacobian(int elemIndex, double time, TacsScalar alpha,
                              TacsScalar beta, TacsScalar gamma,
                              const TacsScalar Xpts[], const TacsScalar vars[],
                              const TacsScalar dvars[],
                              const TacsScalar ddvars[], TacsScalar res[],
                              TacsScalar J[]) {
  // The step length
#ifdef TACS_USE_COMPLEX
  const double dh = 1e-30;
#else
  const double dh = 1e-7;
#endif  // TACS_USE_COMPLEX

  // Get the number of variables
  int nvars = getNumVariables();

  // Call the residual implementation
  if (res) {
    addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);
  }

  // Original and perturbed residual vectors
  TacsScalar *Rtmp1 = new TacsScalar[nvars];
  TacsScalar *Rtmp2 = new TacsScalar[nvars];

  // Perturbed state vectors
  TacsScalar *qTmp = new TacsScalar[nvars];
  TacsScalar *qdotTmp = new TacsScalar[nvars];
  TacsScalar *qddotTmp = new TacsScalar[nvars];

  // Copy the state variables into pstate
  memcpy(qTmp, vars, nvars * sizeof(TacsScalar));
  memcpy(qdotTmp, dvars, nvars * sizeof(TacsScalar));
  memcpy(qddotTmp, ddvars, nvars * sizeof(TacsScalar));

  // Perturb each state variable and find the residual
  if (TacsRealPart(alpha) != 0.0) {
    for (int i = 0; i < nvars; i++) {
#ifdef TACS_USE_COMPLEX
      qTmp[i] = vars[i] + TacsScalar(0.0, dh);

      // Assemble the unperturbed residual
      memset(Rtmp1, 0, nvars * sizeof(TacsScalar));
      addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp1);

      // Find the approximated jacobian
      for (int j = 0; j < nvars; j++) {
        J[j * nvars + i] += alpha * TacsImagPart(Rtmp1[j]) / dh;
      }
#else
      // Perturb the i-th variable
      qTmp[i] = vars[i] + dh;

      // Assemble the unperturbed residual
      memset(Rtmp1, 0, nvars * sizeof(TacsScalar));
      addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp1);

      // Perturb the i-th variable
      qTmp[i] = vars[i] - dh;

      // Assemble the unperturbed residual
      memset(Rtmp2, 0, nvars * sizeof(TacsScalar));
      addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp2);

      // Find the approximated jacobian
      for (int j = 0; j < nvars; j++) {
        J[j * nvars + i] += 0.5 * alpha * (Rtmp1[j] - Rtmp2[j]) / dh;
      }
#endif  // TACS_USE_COMPLEX
      // Restore the i-th variable
      qTmp[i] = vars[i];
    }
  }

  if (TacsRealPart(beta) != 0.0) {
    // Perturb each state variable and find the residual
    for (int i = 0; i < nvars; i++) {
#ifdef TACS_USE_COMPLEX
      qdotTmp[i] = dvars[i] + TacsScalar(0.0, dh);

      // Assemble the unperturbed residual
      memset(Rtmp1, 0, nvars * sizeof(TacsScalar));
      addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp1);

      // Find the approximated jacobian
      for (int j = 0; j < nvars; j++) {
        J[j * nvars + i] += beta * TacsImagPart(Rtmp1[j]) / dh;
      }
#else
      // Perturb the i-th variable
      qdotTmp[i] = dvars[i] + dh;

      // Assemble the unperturbed residual
      memset(Rtmp1, 0, nvars * sizeof(TacsScalar));
      addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp1);

      // Perturb the i-th variable
      qdotTmp[i] = dvars[i] - dh;

      // Assemble the unperturbed residual
      memset(Rtmp2, 0, nvars * sizeof(TacsScalar));
      addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp2);

      // Find the approximated jacobian
      for (int j = 0; j < nvars; j++) {
        J[j * nvars + i] += 0.5 * beta * (Rtmp1[j] - Rtmp2[j]) / dh;
      }
#endif  // TACS_USE_COMPLEX
      // Restore the i-th variable
      qdotTmp[i] = dvars[i];
    }
  }

  if (TacsRealPart(gamma) != 0.0) {
    // Perturb each state variable and find the residual
    for (int i = 0; i < nvars; i++) {
#ifdef TACS_USE_COMPLEX
      qddotTmp[i] = ddvars[i] + TacsScalar(0.0, dh);

      // Assemble the unperturbed residual
      memset(Rtmp1, 0, nvars * sizeof(TacsScalar));
      addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp1);

      // Find the approximated jacobian
      for (int j = 0; j < nvars; j++) {
        J[j * nvars + i] += gamma * TacsImagPart(Rtmp1[j]) / dh;
      }
#else
      // Perturb the i-th variable
      qddotTmp[i] = ddvars[i] + dh;

      // Assemble the unperturbed residual
      memset(Rtmp1, 0, nvars * sizeof(TacsScalar));
      addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp1);

      // Perturb the i-th variable
      qddotTmp[i] = ddvars[i] - dh;

      // Assemble the unperturbed residual
      memset(Rtmp2, 0, nvars * sizeof(TacsScalar));
      addResidual(elemIndex, time, Xpts, qTmp, qdotTmp, qddotTmp, Rtmp2);

      // Find the approximated jacobian
      for (int j = 0; j < nvars; j++) {
        J[j * nvars + i] += 0.5 * gamma * (Rtmp1[j] - Rtmp2[j]) / dh;
      }
#endif  // TACS_USE_COMPLEX
      // Restore the i-th variable
      qddotTmp[i] = ddvars[i];
    }
  }

  delete[] Rtmp1;
  delete[] Rtmp2;
  delete[] qTmp;
  delete[] qdotTmp;
  delete[] qddotTmp;
}

void TACSElement::addAdjResProduct(
    int elemIndex, double time, TacsScalar scale, const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int dvLen, TacsScalar dfdx[]) {
  // The step length
#ifdef TACS_USE_COMPLEX
  const double dh = 1e-30;
#else
  const double dh = 1e-7;
#endif  // TACS_USE_COMPLEX

  TacsScalar *x = new TacsScalar[dvLen];
  getDesignVars(elemIndex, dvLen, x);

  int nvars = getNumVariables();
  TacsScalar *res = new TacsScalar[nvars];
  TacsScalar *tmp = new TacsScalar[nvars];

  memset(res, 0, nvars * sizeof(TacsScalar));
  addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);

  for (int k = 0; k < dvLen; k++) {
    TacsScalar xt = x[k];

#ifdef TACS_USE_COMPLEX
    x[k] = xt + TacsScalar(0.0, dh);
#else
    x[k] = xt + dh;
#endif  // TACS_USE_COMPLEX
    setDesignVars(elemIndex, dvLen, x);

    memset(tmp, 0, nvars * sizeof(TacsScalar));
    addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, tmp);

    TacsScalar product = 0.0;
#ifdef TACS_USE_COMPLEX
    for (int i = 0; i < nvars; i++) {
      product += psi[i] * TacsImagPart(tmp[i]) / dh;
    }
#else
    if (fdOrder < 2) {
      // Use first-order forward differencing
      for (int i = 0; i < nvars; i++) {
        product += psi[i] * (tmp[i] - res[i]) / dh;
      }
    } else {
      // Use second-order central differencing
      x[k] = xt - dh;  //  backward step
      setDesignVars(elemIndex, dvLen, x);
      memset(res, 0, nvars * sizeof(TacsScalar));
      addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);
      // Central difference
      for (int i = 0; i < nvars; i++) {
        product += psi[i] * (tmp[i] - res[i]) / (2.0 * dh);
      }
    }
#endif  // TACS_USE_COMPLEX

    dfdx[k] += scale * product;
    x[k] = xt;
  }

  // Reset the design variable values
  setDesignVars(elemIndex, dvLen, x);

  delete[] x;
  delete[] res;
  delete[] tmp;
}

void TACSElement::addAdjResXptProduct(
    int elemIndex, double time, TacsScalar scale, const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar fXptSens[]) {
  // The step length
#ifdef TACS_USE_COMPLEX
  const double dh = 1e-30;
#else
  const double dh = 1e-7;
#endif  // TACS_USE_COMPLEX

  int nnodes = getNumNodes();
  TacsScalar *X = new TacsScalar[3 * nnodes];
  memcpy(X, Xpts, 3 * nnodes * sizeof(TacsScalar));

  int nvars = getNumVariables();
  TacsScalar *res = new TacsScalar[nvars];
  TacsScalar *tmp = new TacsScalar[nvars];

  memset(res, 0, nvars * sizeof(TacsScalar));
  addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);

  for (int k = 0; k < 3 * nnodes; k++) {
#ifdef TACS_USE_COMPLEX
    X[k] = Xpts[k] + TacsScalar(0.0, dh);
#else
    X[k] = Xpts[k] + dh;
#endif  // TACS_USE_COMPLEX
    memset(tmp, 0, nvars * sizeof(TacsScalar));
    addResidual(elemIndex, time, X, vars, dvars, ddvars, tmp);

    TacsScalar product = 0.0;
#ifdef TACS_USE_COMPLEX
    for (int i = 0; i < nvars; i++) {
      product += psi[i] * TacsImagPart(tmp[i]) / dh;
    }
#else
    if (fdOrder < 2) {
      // Use first-order forward differencing
      for (int i = 0; i < nvars; i++) {
        product += psi[i] * (tmp[i] - res[i]) / (dh);
      }
    } else {
      // Use second-order central differencing
      X[k] = Xpts[k] - dh;  //  backward step
      memset(res, 0, nvars * sizeof(TacsScalar));
      addResidual(elemIndex, time, X, vars, dvars, ddvars, res);
      // Central difference
      for (int i = 0; i < nvars; i++) {
        product += psi[i] * (tmp[i] - res[i]) / (2.0 * dh);
      }
    }
#endif  // TACS_USE_COMPLEX

    fXptSens[k] += scale * product;

    X[k] = Xpts[k];
  }

  delete[] X;
  delete[] res;
  delete[] tmp;
}

void TACSElement::addPointQuantityDVSens(
    int elemIndex, int quantityType, double time, TacsScalar scale, int n,
    double pt[], const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[],
    const TacsScalar dfdq[], int dvLen, TacsScalar dfdx[]) {
  // The step length
#ifdef TACS_USE_COMPLEX
  const double dh = 1e-30;
#else
  const double dh = 1e-7;
#endif  // TACS_USE_COMPLEX

  TacsScalar *x = new TacsScalar[dvLen];
  getDesignVars(elemIndex, dvLen, x);

  TacsScalar detXd = 0.0;
  int nvals = evalPointQuantity(elemIndex, quantityType, time, n, pt, Xpts,
                                vars, dvars, ddvars, &detXd, NULL);
  if (nvals >= 1) {
    TacsScalar *q0 = new TacsScalar[nvals];
    TacsScalar *q1 = new TacsScalar[nvals];
    TacsScalar *fd = new TacsScalar[nvals];
    evalPointQuantity(elemIndex, quantityType, time, n, pt, Xpts, vars, dvars,
                      ddvars, &detXd, q0);

    for (int k = 0; k < dvLen; k++) {
      TacsScalar xt = x[k];

#ifdef TACS_USE_COMPLEX
      x[k] = xt + TacsScalar(0.0, dh);
#else
      x[k] = xt + dh;
#endif  // TACS_USE_COMPLEX
      setDesignVars(elemIndex, dvLen, x);

      evalPointQuantity(elemIndex, quantityType, time, n, pt, Xpts, vars, dvars,
                        ddvars, &detXd, q1);

#ifdef TACS_USE_COMPLEX
      for (int i = 0; i < nvals; i++) {
        fd[i] = TacsImagPart(q1[i]) / dh;
      }
#else
      if (fdOrder < 2) {
        // Use first-order forward differencing
        for (int i = 0; i < nvals; i++) {
          fd[i] = (q1[i] - q0[i]) / dh;
        }
      } else {
        // Use second-order central differencing
        x[k] = xt - dh;  //  backward step
        setDesignVars(elemIndex, dvLen, x);
        evalPointQuantity(elemIndex, quantityType, time, n, pt, Xpts, vars,
                          dvars, ddvars, &detXd, q0);
        // Central difference
        for (int i = 0; i < nvals; i++) {
          fd[i] = (q1[i] - q0[i]) / (2.0 * dh);
        }
      }
#endif  // TACS_USE_COMPLEX
      for (int i = 0; i < nvals; i++) {
        dfdx[k] += scale * dfdq[i] * fd[i];
      }

      x[k] = xt;
    }

    // Reset the design variable values
    setDesignVars(elemIndex, dvLen, x);

    delete[] q0;
    delete[] q1;
    delete[] fd;
  }
  delete[] x;
}

void TACSElement::addPointQuantitySVSens(
    int elemIndex, int quantityType, double time, TacsScalar alpha,
    TacsScalar beta, TacsScalar gamma, int n, double pt[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], const TacsScalar dfdq[], TacsScalar dfdu[]) {
  // The step length
#ifdef TACS_USE_COMPLEX
  const double dh = 1e-30;
#else
  const double dh = 1e-7;
#endif  // TACS_USE_COMPLEX

  TacsScalar detXd = 0.0;
  int nvals = evalPointQuantity(elemIndex, quantityType, time, n, pt, Xpts,
                                vars, dvars, ddvars, &detXd, NULL);
  if (nvals >= 1) {
    int nvars = getNumVariables();
    TacsScalar *v = new TacsScalar[nvars];
    TacsScalar *dv = new TacsScalar[nvars];
    TacsScalar *ddv = new TacsScalar[nvars];
    memcpy(v, vars, nvars * sizeof(TacsScalar));
    memcpy(dv, dvars, nvars * sizeof(TacsScalar));
    memcpy(ddv, ddvars, nvars * sizeof(TacsScalar));

    TacsScalar *q0 = new TacsScalar[nvals];
    TacsScalar *q1 = new TacsScalar[nvals];
    TacsScalar *fd = new TacsScalar[nvals];
    evalPointQuantity(elemIndex, quantityType, time, n, pt, Xpts, vars, dvars,
                      ddvars, &detXd, q0);

    for (int k = 0; k < nvars; k++) {
#ifdef TACS_USE_COMPLEX
      v[k] = vars[k] + alpha * TacsScalar(0.0, dh);
      dv[k] = dvars[k] + beta * TacsScalar(0.0, dh);
      ddv[k] = ddvars[k] + gamma * TacsScalar(0.0, dh);
#else
      v[k] = vars[k] + alpha * dh;
      dv[k] = dvars[k] + beta * dh;
      ddv[k] = ddvars[k] + gamma * dh;
#endif  // TACS_USE_COMPLEX

      evalPointQuantity(elemIndex, quantityType, time, n, pt, Xpts, v, dv, ddv,
                        &detXd, q1);

#ifdef TACS_USE_COMPLEX
      for (int i = 0; i < nvals; i++) {
        fd[i] = TacsImagPart(q1[i]) / dh;
      }
#else
      if (fdOrder < 2) {
        // Use first-order forward differencing
        for (int i = 0; i < nvals; i++) {
          fd[i] = (q1[i] - q0[i]) / dh;
        }
      } else {
        // Use second-order central differencing
        v[k] = vars[k] - alpha * dh;      //  backward step
        dv[k] = dvars[k] - beta * dh;     //  backward step
        ddv[k] = ddvars[k] - gamma * dh;  //  backward step
        evalPointQuantity(elemIndex, quantityType, time, n, pt, Xpts, v, dv,
                          ddv, &detXd, q0);
        // Central difference
        for (int i = 0; i < nvals; i++) {
          fd[i] = (q1[i] - q0[i]) / (2.0 * dh);
        }
      }
#endif  // TACS_USE_COMPLEX

      for (int i = 0; i < nvals; i++) {
        dfdu[k] += dfdq[i] * fd[i];
      }

      v[k] = vars[k];
      dv[k] = dvars[k];
      ddv[k] = ddvars[k];
    }

    delete[] q0;
    delete[] q1;
    delete[] fd;
    delete[] v;
    delete[] dv;
    delete[] ddv;
  }
}

void TACSElement::addPointQuantityXptSens(
    int elemIndex, int quantityType, double time, TacsScalar scale, int n,
    double pt[], const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[],
    const TacsScalar dfddetXd, const TacsScalar dfdq[], TacsScalar dfdXpts[]) {
  // The step length
#ifdef TACS_USE_COMPLEX
  const double dh = 1e-30;
#else
  const double dh = 1e-7;
#endif  // TACS_USE_COMPLEX

  TacsScalar detXd = 0.0;
  int nvals = evalPointQuantity(elemIndex, quantityType, time, n, pt, Xpts,
                                vars, dvars, ddvars, &detXd, NULL);
  if (nvals >= 1) {
    int nnodes = getNumNodes();
    TacsScalar *X = new TacsScalar[3 * nnodes];
    memcpy(X, Xpts, 3 * nnodes * sizeof(TacsScalar));

    TacsScalar detXd0 = 0.0;
    TacsScalar *q0 = new TacsScalar[nvals];
    TacsScalar *q1 = new TacsScalar[nvals];
    TacsScalar *fd = new TacsScalar[nvals];
    evalPointQuantity(elemIndex, quantityType, time, n, pt, Xpts, vars, dvars,
                      ddvars, &detXd0, q0);

    for (int k = 0; k < 3 * nnodes; k++) {
#ifdef TACS_USE_COMPLEX
      X[k] = Xpts[k] + TacsScalar(0.0, dh);
#else
      X[k] = Xpts[k] + dh;
#endif  // TACS_USE_COMPLEX
      TacsScalar detXd1 = 0.0;
      evalPointQuantity(elemIndex, quantityType, time, n, pt, X, vars, dvars,
                        ddvars, &detXd1, q1);

      TacsScalar fddetXd = 0.0;
#ifdef TACS_USE_COMPLEX
      for (int i = 0; i < nvals; i++) {
        fd[i] = TacsImagPart(q1[i]) / dh;
      }
      fddetXd = TacsImagPart(detXd1) / dh;
#else
      if (fdOrder < 2) {
        // Use first-order forward differencing
        for (int i = 0; i < nvals; i++) {
          fd[i] = (q1[i] - q0[i]) / dh;
        }
        fddetXd = (detXd1 - detXd0) / dh;
      } else {
        // Use second-order central differencing
        X[k] = Xpts[k] - dh;  //  backward step
        evalPointQuantity(elemIndex, quantityType, time, n, pt, X, vars, dvars,
                          ddvars, &detXd0, q0);
        // Central difference
        for (int i = 0; i < nvals; i++) {
          fd[i] = (q1[i] - q0[i]) / (2.0 * dh);
        }
        fddetXd = (detXd1 - detXd0) / (2.0 * dh);
      }
#endif  // TACS_USE_COMPLEX
      for (int i = 0; i < nvals; i++) {
        dfdXpts[k] += scale * dfdq[i] * fd[i];
      }
      dfdXpts[k] += scale * dfddetXd * fddetXd;

      X[k] = Xpts[k];
    }

    delete[] X;
    delete[] q0;
    delete[] q1;
    delete[] fd;
  }
}
