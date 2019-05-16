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

#include "ThermalKSFailure.h"
#include "TACSAssembler.h"
#include "FElibrary.h"
#include "ThermoElements.h"
#include "CoupledThermoSolidStiffness.h"
#include "CoupledThermoPlaneStressStiffness.h"

/*
  The context for the TACSKSFailure function
*/
class KSFunctionCtx : public TACSFunctionCtx {
 public:
  KSFunctionCtx( TACSFunction *ks,
                 int maxStrains,
                 int maxNodes ){
    maxFail = -1e20;
    ksFailSum = 0.0;

    // Allocate the working array
    work = new TacsScalar[2*maxStrains + 3*maxNodes];

    // Set the pointers into the work array
    strain = &work[0];
    failSens = &work[maxStrains];
    hXptSens = &work[2*maxStrains];
  }
  ~KSFunctionCtx(){
    delete [] work;
  }

  // Data to be used for the function computation
  TacsScalar maxFail;
  TacsScalar ksFailSum;
  TacsScalar *strain;
  TacsScalar *failSens;
  TacsScalar *hXptSens;
  TacsScalar *work;
};

/*
  Initialize the TACSThermalKSFailure class properties
*/
TACSThermalKSFailure::TACSThermalKSFailure( TACSAssembler *_tacs,
                                            double _ksWeight,
                                            TACSKSFailure::KSConstitutiveFunction func,
                                            double _alpha ):
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN,
             TACSFunction::TWO_STAGE, 0){
  ksWeight = _ksWeight;
  alpha = _alpha;
  conType = func;
  ksType = TACSKSFailure::CONTINUOUS;
  loadFactor = 1.0;

  // Initialize the maximum failure value and KS sum to default values
  // that will be overwritten later.
  maxFail = -1e20;
  ksFailSum = 0.0;

  // Get the max number of nodes/stresses
  maxNumNodes = _tacs->getMaxElementNodes();
  maxNumStrains = _tacs->getMaxElementStrains();
  if (maxNumStrains == 6){
    is_3d = 1;
    is_2d = 0;
  }
  else if (maxNumStrains == 3){
    is_2d = 1;
    is_3d = 0;
  }
}

TACSThermalKSFailure::~TACSThermalKSFailure(){}

/*
  TACSThermalKSFailure function name
*/
const char * TACSThermalKSFailure::funcName = "TACSThermalKSFailure";

/*
  Set the KS aggregation type
*/
void TACSThermalKSFailure::setKSFailureType( enum TACSKSFailure::KSFailureType type ){
  ksType = type;
}

/*
  Retrieve the KS aggregation weight
*/
double TACSThermalKSFailure::getParameter(){
  return ksWeight;
}

/*
  Set the KS aggregation parameter
*/
void TACSThermalKSFailure::setParameter( double _ksWeight ){
  ksWeight = _ksWeight;
}

/*
  Set the load factor to some value greater than or equal to 1.0
*/
void TACSThermalKSFailure::setLoadFactor( TacsScalar _loadFactor ){
  if (TacsRealPart(_loadFactor) >= 1.0){
    loadFactor = _loadFactor;
  }
}

/*
  Return the function name
*/
const char *TACSThermalKSFailure::functionName(){
  return funcName;
}

/*
  Retrieve the function value
*/
TacsScalar TACSThermalKSFailure::getFunctionValue(){
  // Compute the final value of the KS function on all processors
  TacsScalar ksFail = maxFail + log(ksFailSum/alpha)/ksWeight;

  return ksFail;
}

/*
  Allocate and return the function-specific context
*/
TACSFunctionCtx *TACSThermalKSFailure::createFunctionCtx(){
  return new KSFunctionCtx(this, maxNumStrains, maxNumNodes);
}

/*
  Initialize the internal values stored within the KS function
*/
void TACSThermalKSFailure::initEvaluation( EvaluationType ftype ){
  if (ftype == TACSFunction::INITIALIZE){
    maxFail = -1e20;
  }
  else if (ftype == TACSFunction::INTEGRATE){
    ksFailSum = 0.0;
  }
}

/*
  Reduce the function values across all MPI processes
*/
void TACSThermalKSFailure::finalEvaluation( EvaluationType ftype ){
  if (ftype == TACSFunction::INITIALIZE){
    // Distribute the values of the KS function computed on this domain
    TacsScalar temp = maxFail;
    MPI_Allreduce(&temp, &maxFail, 1, TACS_MPI_TYPE,
                  TACS_MPI_MAX, tacs->getMPIComm());
  }
  else {
    // Find the sum of the ks contributions from all processes
    TacsScalar temp = ksFailSum;
    MPI_Allreduce(&temp, &ksFailSum, 1, TACS_MPI_TYPE,
                  MPI_SUM, tacs->getMPIComm());
  }
}

/*
  Initialize the context for either integration or initialization
*/
void TACSThermalKSFailure::initThread( const double tcoef,
                                EvaluationType ftype,
                                TACSFunctionCtx *fctx ){
  KSFunctionCtx *ctx = dynamic_cast<KSFunctionCtx*>(fctx);
  if (ctx){
    if (ftype == TACSFunction::INITIALIZE){
      ctx->maxFail = -1e20;
      ctx->ksFailSum = 0.0;
    }
  }
}

/*
  Perform the element-wise evaluation of the TACSThermalKSFailure function.
*/
void TACSThermalKSFailure::elementWiseEval( EvaluationType ftype,
                                            TACSElement *element, int elemNum,
                                            const TacsScalar Xpts[],
                                            const TacsScalar vars[],
                                            const TacsScalar dvars[],
                                            const TacsScalar ddvars[],
                                            TACSFunctionCtx *fctx ){
  KSFunctionCtx *ctx = dynamic_cast<KSFunctionCtx*>(fctx);

  if (ctx){
    // Retrieve the number of stress components for this element
    int numStresses = element->numStresses();
    int numNodes = element->numNodes();
    
    // Get the number of quadrature points for this element
    int numGauss = element->getNumGaussPts();
    double N[numNodes];
   
    // Get the constitutive object for this element
    TACSConstitutive *constitutive = element->getConstitutive();
    if (constitutive){
      // Set the strain buffer
      TacsScalar *strain = ctx->strain;
      if (ftype == TACSFunction::INITIALIZE){
        // With the first iteration, find the maximum over the domain
        for ( int i = 0; i < numGauss; i++ ){
          // Get the Gauss points one at a time
          double pt[3];
          TacsScalar T[] = {0.0};
          element->getGaussWtsPts(i, pt);
          // Get the strain B*u and temperature dT
          if (is_3d){
            ThermoSolid *elem = dynamic_cast<ThermoSolid*>(element);
            if (elem){
              elem->getShapeFunctions(pt, N);
              elem->getStrain(strain, pt, Xpts, vars);
              elem->getTemperature(T, N, vars);
            }
          }
          else{
            ThermoQuad* elem = dynamic_cast<ThermoQuad*>(element);
            if (elem){
              elem->getShapeFunctions(pt, N);
              elem->getStrain(strain, pt, Xpts, vars);
              elem->getTemperature(T, N, vars);
            }
          }
          
          // Scale the strain by the load factor
          for ( int k = 0; k < numStresses; k++ ){
            strain[k] *= loadFactor;
          }

          T[0] *= loadFactor;
          // Determine the failure criteria
          TacsScalar fail;
          // Test constitutive type
          if (is_3d){
            CoupledThermoSolidStiffness *con =
              dynamic_cast<CoupledThermoSolidStiffness*>(constitutive);

            if (conType == TACSKSFailure::TACSKSFailure::FAILURE && con){
              con->failure(pt, T, strain, &fail);
            }
          }
          else {
            CoupledThermoPlaneStressStiffness *con =
              dynamic_cast<CoupledThermoPlaneStressStiffness*>(constitutive);

            if (conType == TACSKSFailure::FAILURE && con){
              con->failure(pt, T, strain, &fail);
            }
          }
          // Set the maximum failure load
          if (TacsRealPart(fail) > TacsRealPart(ctx->maxFail)){
            ctx->maxFail = fail;
          }
        }
      }
      else {
        for ( int i = 0; i < numGauss; i++ ){
          // Get the Gauss points one at a time
          double pt[3];
          double weight = element->getGaussWtsPts(i, pt);
          TacsScalar T[] = {0.0};

          // Get the strain B*u and temperature dT
          if (is_3d){
            ThermoSolid *elem = dynamic_cast<ThermoSolid*>(element);
            if (elem){
              elem->getShapeFunctions(pt, N);
              elem->getStrain(strain, pt, Xpts, vars);
              elem->getTemperature(T, N, vars);
            }
          }
          else{
            ThermoQuad* elem = dynamic_cast<ThermoQuad*>(element);
            if (elem){
              elem->getShapeFunctions(pt, N);
              elem->getStrain(strain, pt, Xpts, vars);
              elem->getTemperature(T, N, vars);
            }
          }
          // Scale the strain by the load factor
          for ( int k = 0; k < numStresses; k++ ){
            strain[k] *= loadFactor;
          }
          T[0] *= loadFactor;
          // Determine the failure criteria again
          // Test constitutive type
          TacsScalar fail;
          if (is_3d){
            CoupledThermoSolidStiffness *con =
              dynamic_cast<CoupledThermoSolidStiffness*>(constitutive);
            if (conType == TACSKSFailure::FAILURE && con){
              con->failure(pt, T, strain, &fail);
            }
          }
          else {
            CoupledThermoPlaneStressStiffness *con =
              dynamic_cast<CoupledThermoPlaneStressStiffness*>(constitutive);

            if (conType == TACSKSFailure::FAILURE && con){
              con->failure(pt, T, strain, &fail);
            }
          }
          // Add the failure load to the sum
          TacsScalar fexp = exp(ksWeight*(fail - maxFail));
          if (ksType == TACSKSFailure::DISCRETE){
            ctx->ksFailSum += fexp;
          }
          else {
            // Get the determinant of the Jacobian
            TacsScalar h = element->getDetJacobian(pt, Xpts);
            ctx->ksFailSum += h*weight*fexp;
          }
        }
      }
    }
  }
}

/*
  For each thread used to evaluate the function, call the
  post-evaluation code once.
*/
void TACSThermalKSFailure::finalThread( const double tcoef,
                                        EvaluationType ftype,
                                        TACSFunctionCtx *fctx ){
  KSFunctionCtx *ctx = dynamic_cast<KSFunctionCtx*>(fctx);

  if (ctx){
    if (ftype == TACSFunction::INITIALIZE){
      if (TacsRealPart(ctx->maxFail) > TacsRealPart(maxFail)){
        maxFail = ctx->maxFail;
      }
    }
    else {
      ksFailSum += tcoef*ctx->ksFailSum;
    }
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.
*/
void TACSThermalKSFailure::getElementSVSens( double alpha, double beta, double gamma,
                                             TacsScalar *elemSVSens,
                                             TACSElement *element, int elemNum,
                                             const TacsScalar Xpts[],
                                             const TacsScalar vars[],
                                             const TacsScalar dvars[],
                                             const TacsScalar ddvars[],
                                             TACSFunctionCtx *fctx ){
  KSFunctionCtx *ctx = dynamic_cast<KSFunctionCtx*>(fctx);

  // Zero the derivative of the function w.r.t. the element state
  // variables
  int numVars = element->numVariables();
  memset(elemSVSens, 0, numVars*sizeof(TacsScalar));
  if (ctx){
    // Get the number of stress components and total number of variables
    // for this element.
    int numStresses = element->numStresses();
    int numNodes = element->numNodes();
    // Get the quadrature scheme information
    int numGauss = element->getNumGaussPts();

    double N[numNodes];
    int nvars = 1;
    
    // Get the constitutive object
    TACSConstitutive *constitutive = element->getConstitutive();
    if (constitutive){
      // Set pointers into the buffer
      TacsScalar *strain = ctx->strain;
      TacsScalar *failSens = ctx->failSens;
      if (is_3d){
        CoupledThermoSolidStiffness *con =
          dynamic_cast<CoupledThermoSolidStiffness*>(constitutive);
        if (con){
          nvars = con->getVarsPerNode();
        }
      }
      else {
        CoupledThermoPlaneStressStiffness *con =
          dynamic_cast<CoupledThermoPlaneStressStiffness*>(constitutive);
        if (con){
          nvars = con->getVarsPerNode();
        }
      }
      for ( int i = 0; i < numGauss; i++ ){
        double pt[3];
        TacsScalar T[] = {0.0};
        double weight = element->getGaussWtsPts(i, pt);

        // Get the strain B*u and temperature dT
        if (is_3d){
          ThermoSolid *elem = dynamic_cast<ThermoSolid*>(element);
          if (elem){
            elem->getShapeFunctions(pt, N);
            elem->getStrain(strain, pt, Xpts, vars);
            elem->getTemperature(T, N, vars);
          }
        }
        else{
          ThermoQuad* elem = dynamic_cast<ThermoQuad*>(element);
          if (elem){
            elem->getShapeFunctions(pt, N);
            elem->getStrain(strain, pt, Xpts, vars);
            elem->getTemperature(T, N, vars);
          }
        }

        for ( int k = 0; k < numStresses; k++ ){
          strain[k] *= loadFactor;
        }
        T[0] *= loadFactor;

        // Compute the failure term
        TacsScalar fail;
        if (is_3d){
          CoupledThermoSolidStiffness *con =
            dynamic_cast<CoupledThermoSolidStiffness*>(constitutive);
          if (conType == TACSKSFailure::FAILURE && con){
            con->failure(pt, T, strain, &fail);
          }
        }
        else {
          CoupledThermoPlaneStressStiffness *con =
            dynamic_cast<CoupledThermoPlaneStressStiffness*>(constitutive);

          if (conType == TACSKSFailure::FAILURE && con){
            con->failure(pt, T, strain, &fail);
          }
        }
        if (conType == TACSKSFailure::FAILURE){
          TacsScalar ksPtWeight = 0.0;
          if (ksType == TACSKSFailure::DISCRETE){
            // d(log(ksFailSum))/du = 1/(ksFailSum)*d(fail)/du
            ksPtWeight = loadFactor*exp(ksWeight*(fail - maxFail))/ksFailSum;
          }
          else {
            // Get the determinant of the Jacobian
            TacsScalar h = element->getDetJacobian(pt, Xpts);

            ksPtWeight = h*weight*loadFactor*exp(ksWeight*(fail - maxFail))/ksFailSum;
          }

          if (nvars == 1){
            if (is_3d){
              CoupledThermoSolidStiffness *con =
                dynamic_cast<CoupledThermoSolidStiffness*>(constitutive);
              if (con){
                con->failureStrainSens(pt, T, strain, failSens, 0);
              }
              ThermoSolid *elem = dynamic_cast<ThermoSolid*>(element);
              if (elem){
                elem->addEffStrainSVSens(elemSVSens, pt, alpha*ksPtWeight,
                                         failSens, Xpts, vars, 0);
              }
            }
            else {
              CoupledThermoPlaneStressStiffness *con =
                dynamic_cast<CoupledThermoPlaneStressStiffness*>(constitutive);

              if (con){
                con->failureStrainSens(pt, T, strain, failSens, 0);
              }
              ThermoQuad* elem = dynamic_cast<ThermoQuad*>(element);
              if (elem){
                elem->addEffStrainSVSens(elemSVSens, pt, alpha*ksPtWeight,
                                           failSens, Xpts, vars, 0);
              }
            }
          }
          else {
            for ( int j = 1; j < nvars; j++ ){
              // Evaluate the derivative of the failure criteria for i-th
              // material
              if (is_3d){
                CoupledThermoSolidStiffness *con =
                  dynamic_cast<CoupledThermoSolidStiffness*>(constitutive);
                if (con){
                  con->failureStrainSens(pt, T, strain, failSens, j);
                }
                // Evaluate the derivative of the i-th strain w.r.t. to the state
                // variables
                ThermoSolid *elem = dynamic_cast<ThermoSolid*>(element);
                if (elem){
                  elem->addEffStrainSVSens(elemSVSens, pt, alpha*ksPtWeight,
                                           failSens, Xpts, vars, j);
                }
              }
              else {
                CoupledThermoPlaneStressStiffness *con =
                  dynamic_cast<CoupledThermoPlaneStressStiffness*>(constitutive);

                if (con){
                  con->failureStrainSens(pt, T, strain, failSens, j);
                }
                // Evaluate the derivative of the i-th strain w.r.t. to the state
                // variables
                ThermoQuad* elem = dynamic_cast<ThermoQuad*>(element);
                if (elem){
                  elem->addEffStrainSVSens(elemSVSens, pt, alpha*ksPtWeight,
                                             failSens, Xpts, vars, j);
                }
              }              
            }// end j < nvars
          }// else
        }
      }
    }
  }
}

/*
  Determine the derivative of the function with respect to
  the element nodal locations
*/
void TACSThermalKSFailure::getElementXptSens( const double tcoef,
                                              TacsScalar fXptSens[],
                                              TACSElement *element, int elemNum,
                                              const TacsScalar Xpts[],
                                              const TacsScalar vars[],
                                              const TacsScalar dvars[],
                                              const TacsScalar ddvars[],
                                              TACSFunctionCtx *fctx ){
  KSFunctionCtx *ctx = dynamic_cast<KSFunctionCtx*>(fctx);

  // Zero the sensitivity w.r.t. the nodes
  int numNodes = element->numNodes();
  memset(fXptSens, 0, 3*numNodes*sizeof(TacsScalar));

  if (ctx){
    // Get the number of stress components, the total number of
    // variables, and the total number of nodes
    int numStresses = element->numStresses();
    int numNodes = element->numNodes();
    // Get the quadrature scheme information
    int numGauss = element->getNumGaussPts();
    double N[numNodes];
    // Get the constitutive object for this element
    TACSConstitutive *constitutive = element->getConstitutive();
    if (constitutive){
      // Set pointers into the buffer
      TacsScalar *strain = ctx->strain;
      TacsScalar *failSens = ctx->failSens;
      TacsScalar *hXptSens = ctx->hXptSens;

      for ( int i = 0; i < numGauss; i++ ){
        // Get the gauss point
        double pt[3];
        TacsScalar T[] = {0.0};
        double weight = element->getGaussWtsPts(i, pt);

        // Get the strain
        if (is_3d){
          ThermoSolid *elem = dynamic_cast<ThermoSolid*>(element);
          if (elem){
            elem->getShapeFunctions(pt, N);
            elem->getStrain(strain, pt, Xpts, vars);
            elem->getTemperature(T, N, vars);
          }
        }
        else {
          ThermoQuad* elem = dynamic_cast<ThermoQuad*>(element);
          if (elem){
            elem->getShapeFunctions(pt, N);
            elem->getStrain(strain, pt, Xpts, vars);
            elem->getTemperature(T, N, vars);
          }
        }
        
        // Multiply by the load factor
        for ( int k = 0; k < numStresses; k++ ){
          strain[k] *= loadFactor;
        }
        T[0] *= loadFactor;
        // Determine the strain failure criteria
        TacsScalar fail;
        // Test constitutive type
        CoupledThermoSolidStiffness *con =
          dynamic_cast<CoupledThermoSolidStiffness*>(constitutive);
        if (conType == TACSKSFailure::FAILURE && con){
          // Determine the sensitivity of the failure criteria to
          // the design variables and stresses
          con->failure(pt, T, strain, &fail);
          con->failureStrainSens(pt, T, strain, failSens);
        }
        else {
          constitutive->buckling(strain, &fail);
          constitutive->bucklingStrainSens(strain, failSens);
        }

        // Compute the sensitivity contribution
        if (ksType == TACSKSFailure::DISCRETE){
          // d(log(ksFailSum))/dx = 1/(ksFailSum)*d(fail)/dx
          TacsScalar ksPtWeight =
            loadFactor*exp(ksWeight*(fail - maxFail))/ksFailSum;

          element->addStrainXptSens(fXptSens, pt, tcoef*ksPtWeight, failSens,
                                    Xpts, vars);
        }
        else {
          // Get the derivative of the determinant of the Jacobian
          // w.r.t. the nodes
          TacsScalar h = element->getDetJacobianXptSens(hXptSens, pt, Xpts);

          // Compute the derivative of the KS functional
          TacsScalar ksExp = exp(ksWeight*(fail - maxFail))/ksFailSum;
          TacsScalar ksHptWeight = tcoef*weight*ksExp/ksWeight;
          TacsScalar ksPtWeight = h*weight*loadFactor*ksExp;

          for ( int j = 0; j < 3*numNodes; j++ ){
            fXptSens[j] += ksHptWeight*hXptSens[j];
          }

          element->addStrainXptSens(fXptSens, pt, tcoef*ksPtWeight, failSens,
                                    Xpts, vars);
        }
      }
    }
  }
}

/*
  Determine the derivative of the function with respect to
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void TACSThermalKSFailure::addElementDVSens( const double tcoef,
                                             TacsScalar *fdvSens, int numDVs,
                                             TACSElement *element, int elemNum,
                                             const TacsScalar Xpts[],
                                             const TacsScalar vars[],
                                             const TacsScalar dvars[],
                                             const TacsScalar ddvars[],
                                             TACSFunctionCtx *fctx ){
  KSFunctionCtx *ctx = dynamic_cast<KSFunctionCtx*>(fctx);

  // Get the constitutive object for this element
  TACSConstitutive *constitutive = element->getConstitutive();
  if (ctx && constitutive){
    // Get the number of stress components, the total number of
    // variables, and the total number of nodes
    int numStresses = element->numStresses();
    int numNodes = element->numNodes();
    // Get the quadrature scheme information
    int numGauss = element->getNumGaussPts();
    double N[numNodes];
    // Set pointers into the buffer
    TacsScalar *strain = ctx->strain;

    for ( int i = 0; i < numGauss; i++ ){
      // Get the gauss point
      double pt[3];
      double weight = element->getGaussWtsPts(i, pt);
      TacsScalar T[] = {0.0};

      // Get the strain B*u and temperature dT
      if (is_3d){
        ThermoSolid *elem = dynamic_cast<ThermoSolid*>(element);
        if (elem){
          elem->getShapeFunctions(pt, N);
          elem->getStrain(strain, pt, Xpts, vars);
          elem->getTemperature(T, N, vars);
        }
      }
      else {
        ThermoQuad* elem = dynamic_cast<ThermoQuad*>(element);
        if (elem){
          elem->getShapeFunctions(pt, N);
          elem->getStrain(strain, pt, Xpts, vars);
          elem->getTemperature(T, N, vars);
        }
      }      

      for ( int k = 0; k < numStresses; k++ ){
        strain[k] *= loadFactor;
      }
      T[0] *= loadFactor;
      // Determine the strain failure criteria
      TacsScalar fail;
      // Test constitutive type
      if (is_3d){
        CoupledThermoSolidStiffness *con =
          dynamic_cast<CoupledThermoSolidStiffness*>(constitutive);
        if (conType == TACSKSFailure::FAILURE && con){
          con->failure(pt, T, strain, &fail);
        }
      }
      else {
        CoupledThermoPlaneStressStiffness *con =
          dynamic_cast<CoupledThermoPlaneStressStiffness*>(constitutive);

        if (conType == TACSKSFailure::FAILURE && con){
          con->failure(pt, T, strain, &fail);
        }
      }

      // Add contribution from the design variable sensitivity
      // of the failure calculation
      // Compute the sensitivity contribution
      TacsScalar ksPtWeight = 0.0;
      if (ksType == TACSKSFailure::DISCRETE){
        // d(log(ksFailSum))/dx = 1/(ksFailSum)*d(fail)/dx
        ksPtWeight = exp(ksWeight*(fail - maxFail))/ksFailSum;
      }
      else {
        // Get the determinant of the Jacobian
        TacsScalar h = element->getDetJacobian(pt, Xpts);
        ksPtWeight = h*weight*exp(ksWeight*(fail - maxFail))/ksFailSum;
      }

      // Add the derivative of the criteria w.r.t. design variables
      if (is_3d){
        CoupledThermoSolidStiffness *con =
          dynamic_cast<CoupledThermoSolidStiffness*>(constitutive);
        if (conType == TACSKSFailure::FAILURE && con){
          con->addFailureDVSens(pt, T, strain, tcoef*ksPtWeight,
                                fdvSens, numDVs);
        }
      }
      else {
        CoupledThermoPlaneStressStiffness *con =
          dynamic_cast<CoupledThermoPlaneStressStiffness*>(constitutive);

        if (conType == TACSKSFailure::FAILURE && con){
          con->addFailureDVSens(pt, T, strain, tcoef*ksPtWeight,
                                fdvSens, numDVs);
        }
      }
    }
  }
}
