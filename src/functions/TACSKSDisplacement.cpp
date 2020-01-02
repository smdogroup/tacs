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

#include "TACSKSDisplacement.h"
#include "TACSAssembler.h"

/*
  Initialize the TACSKSDisplacement class properties
*/
TACSKSDisplacement::TACSKSDisplacement( TACSAssembler *_assembler,
                                        double _ksWeight,
                                        double _alpha ):
  TACSFunction(_assembler, TACSFunction::ENTIRE_DOMAIN,
               TACSFunction::TWO_STAGE, 0){
  ksWeight = _ksWeight;
  alpha = _alpha;

  // Initialize the maximum failure value and KS sum to default values
  // that will be overwritten later.
  maxDisp = -1e20;
  ksDispSum = 0.0;
  invPnorm = 0.0;
}

TACSKSDisplacement::~TACSKSDisplacement(){}

/*
  TACSKSDisplacement function name
*/
const char * TACSKSDisplacement::funcName = "TACSKSDisplacement";

/*
  Retrieve the KS aggregation weight
*/
double TACSKSDisplacement::getParameter(){
  return ksWeight;
}

/*
  Set the KS aggregation parameter
*/
void TACSKSDisplacement::setParameter( double _ksWeight ){
  ksWeight = _ksWeight;
}

/*
  Return the function name
*/
const char *TACSKSDisplacement::getObjectName(){
  return funcName;
}

/*
  Retrieve the function value
*/
TacsScalar TACSKSDisplacement::getFunctionValue(){
  // Compute the final value of the KS function on all processors
  TacsScalar ksDisp = maxDisp + log(ksDispSum/alpha)/ksWeight;

  return ksDisp;
}

/*
  Initialize the internal values stored within the KS function
*/
void TACSKSDisplacement::initEvaluation( EvaluationType ftype ){
  if (ftype == TACSFunction::INITIALIZE){
    maxDisp = -1e20;
  }
  else if (ftype == TACSFunction::INTEGRATE){
    ksDispSum = 0.0;
  }
}

/*
  Reduce the function values across all MPI processes
*/
void TACSKSDisplacement::finalEvaluation( EvaluationType ftype ){
  if (ftype == TACSFunction::INITIALIZE){
    // Distribute the values of the KS function computed on this domain
    TacsScalar temp = maxDisp;
    MPI_Allreduce(&temp, &maxDisp, 1, TACS_MPI_TYPE,
                  TACS_MPI_MAX, assembler->getMPIComm());
  }
  else {
    // Find the sum of the ks contributions from all processes
    TacsScalar temp = ksDispSum;
    MPI_Allreduce(&temp, &ksDispSum, 1, TACS_MPI_TYPE,
                  MPI_SUM, assembler->getMPIComm());
  }
}

/*
  Perform the element-wise evaluation of the TACSKSDisplacement function.
*/
void TACSKSDisplacement::elementWiseEval( EvaluationType ftype,
                                          int elemIndex,
                                          TACSElement *element,
                                          double time,
                                          TacsScalar scale,
                                          const TacsScalar Xpts[],
                                          const TacsScalar vars[],
                                          const TacsScalar dvars[],
                                          const TacsScalar ddvars[] ){
  // Retrieve the number of stress components for this element
  TACSElementBasis *basis = element->getElementBasis();

  if (basis){
    for ( int i = 0; i < basis->getNumQuadraturePoints(); i++ ){
      double pt[3];
      double weight = basis->getQuadraturePoint(i, pt);

      // Evaluate the failure index, and check whether it is an
      // undefined quantity of interest on this element
      TacsScalar fail = 0.0;
      int count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DISPLACEMENT,
                                             time, i, pt,
                                             Xpts, vars, dvars, ddvars,
                                             &fail);

      // Check whether the quantity requested is defined or not
      if (count >= 1){
        if (ftype == TACSFunction::INITIALIZE){
          // Set the maximum failure load
          if (TacsRealPart(fail) > TacsRealPart(maxDisp)){
            maxDisp = fail;
          }
        }
        else {
          // Evaluate the determinant of the Jacobian
          TacsScalar Xd[9], J[9];
          TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);

          TacsScalar fexp = exp(ksWeight*(fail - maxDisp));
          ksDispSum += scale*weight*detJ*fexp;

        }
      }
    }
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.
*/
void TACSKSDisplacement::getElementSVSens( int elemIndex, TACSElement *element,
                                           double time,
                                           TacsScalar alpha,
                                           TacsScalar beta,
                                           TacsScalar gamma,
                                           const TacsScalar Xpts[],
                                           const TacsScalar vars[],
                                           const TacsScalar dvars[],
                                           const TacsScalar ddvars[],
                                           TacsScalar dfdu[] ){
  // Zero the derivative of the function w.r.t. the element state
  // variables
  int numVars = element->getNumVariables();
  memset(dfdu, 0, numVars*sizeof(TacsScalar));

  // Get the element basis class
  TACSElementBasis *basis = element->getElementBasis();

  if (basis){
    for ( int i = 0; i < basis->getNumQuadraturePoints(); i++ ){
      double pt[3];
      double weight = basis->getQuadraturePoint(i, pt);

      TacsScalar fail = 0.0;
      int count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DISPLACEMENT,
                                             time, i, pt,
                                             Xpts, vars, dvars, ddvars,
                                             &fail);

      if (count >= 1){
        // Evaluate the determinant of the Jacobian
        TacsScalar Xd[9], J[9];
        TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);

        // Compute the sensitivity contribution
        TacsScalar ksPtWeight = 0.0;
        ksPtWeight = exp(ksWeight*(fail - maxDisp))/ksDispSum;
        ksPtWeight *= weight*detJ;

        TacsScalar dfdq = ksPtWeight;
        element->addPointQuantitySVSens(elemIndex, TACS_ELEMENT_DISPLACEMENT, time,
                                        alpha, beta, gamma,
                                        i, pt, Xpts, vars, dvars, ddvars,
                                        &dfdq, dfdu);
      }
    }
  }
}

/*
  Determine the derivative of the function with respect to
  the element nodal locations
*/
void TACSKSDisplacement::getElementXptSens( int elemIndex,
                                            TACSElement *element,
                                            double time,
                                            TacsScalar scale,
                                            const TacsScalar Xpts[],
                                            const TacsScalar vars[],
                                            const TacsScalar dvars[],
                                            const TacsScalar ddvars[],
                                            TacsScalar dfdXpts[] ){
  // Zero the derivative of the function w.r.t. the element node
  // locations
  int numNodes = element->getNumNodes();
  memset(dfdXpts, 0, 3*numNodes*sizeof(TacsScalar));

  // Get the element basis class
  TACSElementBasis *basis = element->getElementBasis();

  if (basis){
    for ( int i = 0; i < basis->getNumQuadraturePoints(); i++ ){
      double pt[3];
      double weight = basis->getQuadraturePoint(i, pt);

      TacsScalar fail = 0.0;
      int count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DISPLACEMENT, 
                                             time, i, pt,
                                             Xpts, vars, dvars, ddvars,
                                             &fail);

      if (count >= 1){
        // Evaluate the determinant of the Jacobian
        TacsScalar Xd[9], J[9];
        TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);

        // Compute the sensitivity contribution
        TacsScalar dfdq = 0.0;
        TacsScalar dfddetJ = 0.0;
        TacsScalar expfact = exp(ksWeight*(fail - maxDisp))/ksDispSum;
        dfddetJ = weight*expfact/ksWeight;
        dfdq = weight*detJ*expfact;
        element->addPointQuantityXptSens(elemIndex, TACS_ELEMENT_DISPLACEMENT, time,
                                         scale, i, pt, Xpts, vars, dvars, ddvars,
                                         &dfdq, dfdXpts);
        if (dfddetJ != 0.0){
          basis->addJacobianTransformXptSens(pt, Xd, J, scale*dfddetJ,
                                             NULL, NULL, dfdXpts);
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
void TACSKSDisplacement::addElementDVSens( int elemIndex,
                                           TACSElement *element,
                                           double time,
                                           TacsScalar scale,
                                           const TacsScalar Xpts[],
                                           const TacsScalar vars[],
                                           const TacsScalar dvars[],
                                           const TacsScalar ddvars[],
                                           int dvLen, TacsScalar dfdx[] ){
  // Get the element basis class
  TACSElementBasis *basis = element->getElementBasis();

  if (basis){
    for ( int i = 0; i < basis->getNumQuadraturePoints(); i++ ){
      double pt[3];
      double weight = basis->getQuadraturePoint(i, pt);

      TacsScalar fail = 0.0;
      int count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DISPLACEMENT, 
                                             time, i, pt,
                                             Xpts, vars, dvars, ddvars,
                                             &fail);

      if (count >= 1){
        // Evaluate the determinant of the Jacobian
        TacsScalar Xd[9], J[9];
        TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);

        // Compute the sensitivity contribution
        TacsScalar dfdq = 0.0;
        TacsScalar expfact = exp(ksWeight*(fail - maxDisp))/ksDispSum;
        dfdq = weight*detJ*expfact;

        element->addPointQuantityDVSens(elemIndex, TACS_ELEMENT_DISPLACEMENT, 
                                        time, scale, i, pt,
                                        Xpts, vars, dvars, ddvars,
                                        &dfdq, dvLen, dfdx);
      }
    }
  }
}
