/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2018 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "KSMatTemperature.h"

#include "CoupledThermoPlaneStressStiffness.h"
#include "CoupledThermoSolidStiffness.h"
#include "TACSAssembler.h"

/*
  The context for the TACSKSTemperature function
*/
class KSTemperatureCtx : public TACSFunctionCtx {
 public:
  KSTemperatureCtx(TACSFunction *func, int maxNodes) {
    // Allocate the working array
    N = new double[maxNodes];
    ksSum = 0.0;
    maxValue = new TacsScalar[5];
    for (int i = 0; i < 5; i++) {
      maxValue[i] = -1e20;
    }
  }
  ~KSTemperatureCtx() {
    delete[] N;
    delete[] maxValue;
  }

  // Data to be used for the function computation
  TacsScalar ksSum;
  TacsScalar *maxValue;
  double *N;
};

/*
  Initialize the TACSKSMatTemperature class properties
*/
TACSKSMatTemperature::TACSKSMatTemperature(
    TACSAssembler *_tacs, double _ksWeight,
    TACSKSTemperature::KSTemperatureType _ksType, int _nmats)
    : TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN, TACSFunction::TWO_STAGE,
                   0) {
  maxNumNodes = _tacs->getMaxElementNodes();
  ksWeight = _ksWeight;
  ksType = _ksType;
  nmats = _nmats;
  maxValue = new TacsScalar[nmats];
  for (int i = 0; i < nmats; i++) {
    maxValue[i] = -1e20;
  }
  ksSum = 0.0;
  invPnorm = 0.0;
  int maxNumStrains = _tacs->getMaxElementStrains();
  if (maxNumStrains == 6) {
    is_3d = 1;
    is_2d = 0;
  } else if (maxNumStrains == 3) {
    is_2d = 1;
    is_3d = 0;
  }
}

TACSKSMatTemperature::~TACSKSMatTemperature() {}

/*
  TACSKSMatTemperature function name
*/
const char *TACSKSMatTemperature::funcName = "TACSKSMatTemperature";

/*
  Return the function name
*/
const char *TACSKSMatTemperature::functionName() { return funcName; }

/*
  Set the displacement aggregate type
*/
void TACSKSMatTemperature::setKSDispType(
    TACSKSTemperature::KSTemperatureType _ksType) {
  ksType = _ksType;
}

/*
  Retrieve the function value
*/
TacsScalar TACSKSMatTemperature::getFunctionValue() {
  // Compute the final value of the KS function on all processors
  if (ksType == TACSKSTemperature::CONTINUOUS ||
      ksType == TACSKSTemperature::DISCRETE) {
    // TacsScalar ksValue = maxValue + log(ksSum)/ksWeight;
    TacsScalar ksValue = log(ksSum) / ksWeight;
    for (int i = 0; i < nmats; i++) {
      ksValue += maxValue[i];
    }
    int mpi_rank;
    MPI_Comm_rank(tacs->getMPIComm(), &mpi_rank);
    if (mpi_rank == 0) {
      printf("KS temperature value: %25.10e\n", TacsRealPart(ksValue));
      for (int i = 0; i < nmats; i++) {
        printf("Max temperature value[%d]: %25.10e\n", i,
               TacsRealPart(maxValue[i]));
      }
    }
    return ksValue;
  } else {
    return maxValue[0] *
           pow(ksSum,
               1.0 / ksWeight);  // return maxValue*pow(ksSum, 1.0/ksWeight);
  }
}

/*
  Allocate and return the function-specific context
*/
TACSFunctionCtx *TACSKSMatTemperature::createFunctionCtx() {
  return new KSTemperatureCtx(this, maxNumNodes);
}

/*
  Initialize the internal values stored within the KS function
*/
void TACSKSMatTemperature::initEvaluation(EvaluationType ftype) {
  if (ftype == TACSFunction::INITIALIZE) {
    for (int i = 0; i < nmats; i++) {
      maxValue[i] = -1e20;
    }
  } else if (ftype == TACSFunction::INTEGRATE) {
    ksSum = 0.0;
  }
}

/*
  Reduce the function values across all MPI processes
*/
void TACSKSMatTemperature::finalEvaluation(EvaluationType ftype) {
  if (ftype == TACSFunction::INITIALIZE) {
    // Distribute the values of the KS function computed on this domain
    for (int i = 0; i < nmats; i++) {
      TacsScalar temp = maxValue[i];
      MPI_Allreduce(&temp, &maxValue[i], 1, TACS_MPI_TYPE, TACS_MPI_MAX,
                    tacs->getMPIComm());
    }
  } else {
    // Find the sum of the ks contributions from all processes
    TacsScalar temp = ksSum;
    MPI_Allreduce(&temp, &ksSum, 1, TACS_MPI_TYPE, MPI_SUM, tacs->getMPIComm());

    // Compute the P-norm quantity if needed
    invPnorm = 0.0;
    if (ksType == TACSKSTemperature::PNORM_DISCRETE ||
        ksType == TACSKSTemperature::PNORM_CONTINUOUS) {
      if (ksSum != 0.0) {
        invPnorm = pow(ksSum, (1.0 - ksWeight) / ksWeight);
      }
    }
  }
}

/*
  Initialize the context for either integration or initialization
*/
void TACSKSMatTemperature::initThread(const double tcoef, EvaluationType ftype,
                                      TACSFunctionCtx *fctx) {
  KSTemperatureCtx *ctx = dynamic_cast<KSTemperatureCtx *>(fctx);
  if (ctx) {
    if (ftype == TACSFunction::INITIALIZE) {
      for (int i = 0; i < nmats; i++) {
        ctx->maxValue[i] = -1e20;
      }
      ctx->ksSum = 0.0;
    }
  }
}

/*
  Perform the element-wise evaluation of the TACSKSMatTemperature function.
*/
void TACSKSMatTemperature::elementWiseEval(
    EvaluationType ftype, TACSElement *element, int elemNum,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TACSFunctionCtx *fctx) {
  KSTemperatureCtx *ctx = dynamic_cast<KSTemperatureCtx *>(fctx);
  if (ctx) {
    // Get the number of quadrature points for this element
    const int numGauss = element->getNumGaussPts();
    const int numDisps = element->numDisplacements();
    const int numNodes = element->numNodes();
    TacsScalar *dT = new TacsScalar[nmats];
    memset(dT, 0.0, nmats * sizeof(TacsScalar));
    if (ftype == TACSFunction::INITIALIZE) {
      // With the first iteration, find the maximum over the domain
      for (int i = 0; i < numGauss; i++) {
        // Get the Gauss points one at a time
        double pt[3];
        element->getGaussWtsPts(i, pt);
        element->getShapeFunctions(pt, ctx->N);

        // Evaluate the dot-product with the displacements
        const double *N = ctx->N;
        const TacsScalar *d = vars;

        TacsScalar value = 0.0;
        TacsScalar ns = 0.0;
        for (int j = 0; j < numNodes; j++) {
          ns += N[0];
          value += N[0] * d[numDisps - 1];
          d += numDisps;
          N++;
        }

        // --------------------------------------------------------
        // Get the constitutive object for this element
        TacsScalar value1 = value;
        TACSConstitutive *constitutive = element->getConstitutive();
        if (is_3d) {
          CoupledThermoSolidStiffness *con =
              dynamic_cast<CoupledThermoSolidStiffness *>(constitutive);
          if (con) {
            for (int i = 0; i < nmats; i++) {
              con->maxtemp(pt, value1, &dT[i], i);
              if (TacsRealPart(dT[i]) > TacsRealPart(ctx->maxValue[i])) {
                ctx->maxValue[i] = dT[i];
              }
            }
          }
        } else {
          CoupledThermoPlaneStressStiffness *con =
              dynamic_cast<CoupledThermoPlaneStressStiffness *>(constitutive);
          if (con) {
            for (int i = 0; i < nmats; i++) {
              con->maxtemp(pt, value1, &dT[i], i);
              if (TacsRealPart(dT[i]) > TacsRealPart(ctx->maxValue[i])) {
                ctx->maxValue[i] = dT[i];
              }
            }
          }
        }
      }
    } else {
      // With the first iteration, find the maximum over the domain
      for (int i = 0; i < numGauss; i++) {
        // Get the Gauss points one at a time
        double pt[3];
        double weight = element->getGaussWtsPts(i, pt);
        element->getShapeFunctions(pt, ctx->N);

        // Evaluate the dot-product with the displacements
        const double *N = ctx->N;
        const TacsScalar *d = vars;

        TacsScalar value = 0.0;
        for (int j = 0; j < numNodes; j++) {
          value += N[0] * d[numDisps - 1];
          d += numDisps;
          N++;
        }
        // Add up the contribution from the quadrature
        TacsScalar h = element->getDetJacobian(pt, Xpts);

        // --------------------------------------------------------
        // Get the constitutive object for this element
        TacsScalar value1 = value;
        TACSConstitutive *constitutive = element->getConstitutive();
        if (is_3d) {
          CoupledThermoSolidStiffness *con =
              dynamic_cast<CoupledThermoSolidStiffness *>(constitutive);
          if (con) {
            for (int i = 0; i < nmats; i++) {
              con->maxtemp(pt, value1, &dT[i], i);
            }
          }
        } else {
          CoupledThermoPlaneStressStiffness *con =
              dynamic_cast<CoupledThermoPlaneStressStiffness *>(constitutive);
          if (con) {
            for (int i = 0; i < nmats; i++) {
              con->maxtemp(pt, value1, &dT[i], i);
            }
          }
        }

        // ---------------------------------------------------------
        if (ksType == TACSKSTemperature::CONTINUOUS) {
          for (int i = 0; i < nmats; i++) {
            ctx->ksSum += h * weight * exp(ksWeight * (dT[i] - maxValue[i]));
          }
        }
        // exit(0);
        //  else if (ksType == DISCRETE){
        //    for
        //      ctx->ksSum += exp(ksWeight*(value - maxValue[0]));
        //  }
        //  else if (ksType == PNORM_CONTINUOUS){
        //    ctx->ksSum +=
        //      h*weight*pow(fabs(TacsRealPart(value/maxValue[0])), ksWeight);
        //  }
        //  else if (ksType == PNORM_DISCRETE){
        //    ctx->ksSum += pow(fabs(TacsRealPart(value/maxValue[0])),
        //    ksWeight);
        //  }
      }
    }
  }
}

/*
  For each thread used to evaluate the function, call the
  post-evaluation code once.
*/
void TACSKSMatTemperature::finalThread(const double tcoef, EvaluationType ftype,
                                       TACSFunctionCtx *fctx) {
  KSTemperatureCtx *ctx = dynamic_cast<KSTemperatureCtx *>(fctx);
  if (ctx) {
    if (ftype == TACSFunction::INITIALIZE) {
      for (int i = 0; i < nmats; i++) {
        if (TacsRealPart(ctx->maxValue[i]) > TacsRealPart(maxValue[i])) {
          maxValue[i] = ctx->maxValue[i];
        }
      }
    } else {
      ksSum += tcoef * ctx->ksSum;
    }
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.
*/
void TACSKSMatTemperature::getElementSVSens(
    double alpha, double beta, double gamma, TacsScalar *elemSVSens,
    TACSElement *element, int elemNum, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TACSFunctionCtx *fctx) {
  KSTemperatureCtx *ctx = dynamic_cast<KSTemperatureCtx *>(fctx);

  // Zero the derivative of the function w.r.t. the element state
  // variables
  int numVars = element->numVariables();
  memset(elemSVSens, 0, numVars * sizeof(TacsScalar));
  TacsScalar *dT = new TacsScalar[nmats];
  memset(dT, 0.0, nmats * sizeof(TacsScalar));
  if (ctx) {
    // Get the number of quadrature points for this element
    const int numGauss = element->getNumGaussPts();
    const int numDisps = element->numDisplacements();
    const int numNodes = element->numNodes();

    // With the first iteration, find the maximum over the domain
    for (int i = 0; i < numGauss; i++) {
      // Get the Gauss points one at a time
      double pt[3];
      double weight = element->getGaussWtsPts(i, pt);
      element->getShapeFunctions(pt, ctx->N);

      // Evaluate the dot-product with the displacements
      const double *N = ctx->N;
      const TacsScalar *d = vars;

      TacsScalar value = 0.0;
      for (int j = 0; j < numNodes; j++) {
        value += N[0] * d[numDisps - 1];
        d += numDisps;
        N++;
      }
      // Get the constitutive object for this element
      TacsScalar value1 = value;
      TACSConstitutive *constitutive = element->getConstitutive();
      if (is_3d) {
        CoupledThermoSolidStiffness *con =
            dynamic_cast<CoupledThermoSolidStiffness *>(constitutive);
        if (con) {
          for (int i = 0; i < nmats; i++) {
            con->maxtemp(pt, value1, &dT[i], i);
          }
        }
      } else {
        CoupledThermoPlaneStressStiffness *con =
            dynamic_cast<CoupledThermoPlaneStressStiffness *>(constitutive);
        if (con) {
          for (int i = 0; i < nmats; i++) {
            con->maxtemp(pt, value1, &dT[i], i);
          }
        }
      }
      // Add up the contribution from the quadrature
      TacsScalar h = element->getDetJacobian(pt, Xpts);
      TacsScalar ptWeight = 0.0;

      if (ksType == TACSKSTemperature::CONTINUOUS) {
        for (int i = 0; i < nmats; i++) {
          ptWeight = alpha * h * weight *
                     exp(ksWeight * (dT[i] - maxValue[i])) / ksSum;
          // Reset the shape function pointer and run through the
          // element nodes again to set the derivative
          // Get the weights from design variables
          TacsScalar wx[] = {0.0};
          if (is_3d) {
            CoupledThermoSolidStiffness *con =
                dynamic_cast<CoupledThermoSolidStiffness *>(constitutive);
            if (con) {
              con->maxtempStrainSens(pt, value1, wx, i);
            }
          } else {
            CoupledThermoPlaneStressStiffness *con =
                dynamic_cast<CoupledThermoPlaneStressStiffness *>(constitutive);
            if (con) {
              con->maxtempStrainSens(pt, value1, wx, i);
            }
          }
          N = ctx->N;
          TacsScalar *s = elemSVSens;
          for (int j = 0; j < numNodes; j++) {
            s[numDisps - 1] += ptWeight * N[0] * wx[0];
            s += numDisps;
            N++;
          }
        }
      }
    }
  }
  // else if (ksType == DISCRETE){
  //   ptWeight = alpha*exp(ksWeight*(value - maxValue))/ksSum;
  // }
  // else if (ksType == PNORM_CONTINUOUS){
  //   ptWeight = value*pow(fabs(TacsRealPart(value/maxValue)), ksWeight-2.0);
  //   ptWeight *= alpha*h*weight*invPnorm;
  // }
  // else if (ksType == PNORM_DISCRETE){
  //   ptWeight = value*pow(fabs(TacsRealPart(value/maxValue)), ksWeight-2.0);
  //   ptWeight *= alpha*ksWeight*invPnorm;
  // }
  /*// Get the weights from design variables
  TacsScalar wx[] = {0.0};
  if (is_3d){
    CoupledThermoSolidStiffness *con =
      dynamic_cast<CoupledThermoSolidStiffness*>(constitutive);
    if (con){
      con->maxtempStrainSens(pt, value1, wx);
    }
  }
  else {
    CoupledThermoPlaneStressStiffness *con =
      dynamic_cast<CoupledThermoPlaneStressStiffness*>(constitutive);
    if (con){
      con->maxtempStrainSens(pt, value1, wx);
    }
  }
  // Reset the shape function pointer and run through the
  // element nodes again to set the derivative
  N = ctx->N;
  TacsScalar *s = elemSVSens;
  for ( int j = 0; j < numNodes; j++ ){
    //s[numDisps-1] += ptWeight*N[0];
    s[numDisps-1] += ptWeight*N[0]*wx[0];
    s += numDisps;
    N++;
  }
}
}*/
}

/*
  Determine the derivative of the function with respect to
  the element nodal locations
*/
void TACSKSMatTemperature::getElementXptSens(
    const double tcoef, TacsScalar fXptSens[], TACSElement *element,
    int elemNum, const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[],
    TACSFunctionCtx *fctx) {}

/*
  Determine the derivative of the function with respect to
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void TACSKSMatTemperature::addElementDVSens(
    const double tcoef, TacsScalar *fdvSens, int numDVs, TACSElement *element,
    int elemNum, const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[],
    TACSFunctionCtx *fctx) {
  KSTemperatureCtx *ctx = dynamic_cast<KSTemperatureCtx *>(fctx);
  // Get the constitutive object for this element
  TACSConstitutive *constitutive = element->getConstitutive();
  TacsScalar *dT = new TacsScalar[nmats];
  memset(dT, 0.0, nmats * sizeof(TacsScalar));
  if (ctx) {
    // Get the number of quadrature points for this element
    const int numGauss = element->getNumGaussPts();
    const int numDisps = element->numDisplacements();
    const int numNodes = element->numNodes();
    // With the first iteration, find the maximum over the domain
    for (int i = 0; i < numGauss; i++) {
      // Get the Gauss points one at a time
      double pt[3];
      double weight = element->getGaussWtsPts(i, pt);
      element->getShapeFunctions(pt, ctx->N);

      // Evaluate the dot-product with the displacements
      const double *N = ctx->N;
      const TacsScalar *d = vars;

      TacsScalar value = 0.0;
      for (int j = 0; j < numNodes; j++) {
        value += N[0] * d[numDisps - 1];
        d += numDisps;
        N++;
      }
      // --------------------------------------------------------
      // Get the constitutive object for this element
      TacsScalar value1 = value;
      TACSConstitutive *constitutive = element->getConstitutive();
      if (is_3d) {
        CoupledThermoSolidStiffness *con =
            dynamic_cast<CoupledThermoSolidStiffness *>(constitutive);
        if (con) {
          for (int i = 0; i < nmats; i++) {
            con->maxtemp(pt, value1, &dT[i], i);
          }
        }
      } else {
        CoupledThermoPlaneStressStiffness *con =
            dynamic_cast<CoupledThermoPlaneStressStiffness *>(constitutive);
        if (con) {
          for (int i = 0; i < nmats; i++) {
            con->maxtemp(pt, value1, &dT[i], i);
          }
        }
      }
      // Add up the contribution from the quadrature
      TacsScalar h = element->getDetJacobian(pt, Xpts);
      TacsScalar ptWeight = 0.0;

      if (ksType == TACSKSTemperature::CONTINUOUS) {
        for (int i = 0; i < nmats; i++) {
          ptWeight = h * weight * exp(ksWeight * (dT[i] - maxValue[i])) / ksSum;
          // Add contribution of the relaxation to the design sensitivity
          if (is_3d) {
            CoupledThermoSolidStiffness *con =
                dynamic_cast<CoupledThermoSolidStiffness *>(constitutive);
            if (con) {
              con->addMaxTempDVSens(pt, dT[i], tcoef * ptWeight, fdvSens,
                                    numDVs, i);
            }
          } else {
            CoupledThermoPlaneStressStiffness *con =
                dynamic_cast<CoupledThermoPlaneStressStiffness *>(constitutive);
            if (con) {
              con->addMaxTempDVSens(pt, dT[i], tcoef * ptWeight, fdvSens,
                                    numDVs, i);
            }
          }
        }
      }
      // else if (ksType == DISCRETE){
      //   ptWeight = exp(ksWeight*(value - maxValue))/ksSum;
      // }
      // else if (ksType == PNORM_CONTINUOUS){
      //   ptWeight = value*pow(fabs(TacsRealPart(value/maxValue)),
      //   ksWeight-2.0); ptWeight *= h*weight*invPnorm;
      // }
      // else if (ksType == PNORM_DISCRETE){
      //   ptWeight = value*pow(fabs(TacsRealPart(value/maxValue)),
      //   ksWeight-2.0); ptWeight *= ksWeight*invPnorm;
      // }
      /*// Add contribution of the relaxation to the design sensitivity
      if (is_3d){
        CoupledThermoSolidStiffness *con =
          dynamic_cast<CoupledThermoSolidStiffness*>(constitutive);
        if (con){
          con->addMaxTempDVSens(pt, value1, tcoef*ptWeight,
                                fdvSens, numDVs);
        }
      }
      else {
        CoupledThermoPlaneStressStiffness *con =
          dynamic_cast<CoupledThermoPlaneStressStiffness*>(constitutive);
        if (con){
          con->addMaxTempDVSens(pt, value1, tcoef*ptWeight,
                                fdvSens, numDVs);
        }
      }*/
    }
  }
}
