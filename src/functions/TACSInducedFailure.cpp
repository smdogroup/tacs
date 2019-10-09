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

#include "TACSInducedFailure.h"
#include "TACSAssembler.h"

/*
  Initialize the TACSInducedFailure class properties

  Evaluate the Induced only on the elements specified
*/
TACSInducedFailure::TACSInducedFailure( TACSAssembler *_assembler,
                                        double _P ):
  TACSFunction(_assembler, TACSFunction::ENTIRE_DOMAIN,
               TACSFunction::TWO_STAGE, 0){
  // Set the penalization information
  P = _P;
  normType = EXPONENTIAL;

  // Initialize values
  maxFail = 0.0;
  failNumer = 0.0;
  failDenom = 0.0;
}

/*
  Delete all the allocated data
*/
TACSInducedFailure::~TACSInducedFailure(){}

/*
  The name of the function class
*/
const char * TACSInducedFailure::funcName = "TACSInducedFailure";

/*
  Set the value of P
*/
void TACSInducedFailure::setParameter( double _P ){
  P = _P;
}

/*
  Retrieve the value of P
*/
double TACSInducedFailure::getParameter(){
  return P;
}

/*
  Set the type of p-norm to use
*/
void TACSInducedFailure::setInducedType( enum InducedNormType type ){
  normType = type;
}

/*
  Retrieve the function name
*/
const char *TACSInducedFailure::getObjectName(){ return funcName; }

/*
  Retrieve the function value
*/
TacsScalar TACSInducedFailure::getFunctionValue(){
  // Compute the final value of the function on all processors
  return maxFail*failNumer/failDenom;
}

/*
  Initialize the internal values stored within the KS function
*/
void TACSInducedFailure::initEvaluation( EvaluationType ftype ){
  if (ftype == TACSFunction::INITIALIZE){
    maxFail = -1e20;
  }
  else if (ftype == TACSFunction::INTEGRATE){
    failNumer = 0.0;
    failDenom = 0.0;
  }
}

/*
  Reduce the function values across all MPI processes
*/
void TACSInducedFailure::finalEvaluation( EvaluationType ftype ){
  if (ftype == TACSFunction::INITIALIZE){
    // Distribute the values of the KS function computed on this domain
    TacsScalar temp = maxFail;
    MPI_Allreduce(&temp, &maxFail, 1, TACS_MPI_TYPE,
                  TACS_MPI_MAX, assembler->getMPIComm());
  }
  else if (ftype == TACSFunction::INTEGRATE){
    // Find the sum of the ks contributions from all processes
    TacsScalar in[2], out[2];
    in[0] = failNumer;
    in[1] = failDenom;

    MPI_Allreduce(in, out, 2, TACS_MPI_TYPE, MPI_SUM, assembler->getMPIComm());

    failNumer = out[0];
    failDenom = out[1];
  }
}

/*
  Perform the element-wise evaluation of the TACSInducedFailure
  function.
*/
void TACSInducedFailure::elementWiseEval( EvaluationType ftype,
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
      int count = element->evalPointQuantity(elemIndex, TACS_FAILURE_INDEX,
                                             time, i, pt,
                                             Xpts, vars, dvars, ddvars,
                                             &fail);

      // Check whether the quantity requested is defined or not
      if (count >= 1){
        if (ftype == TACSFunction::INITIALIZE){
          // Set the maximum failure load
          if (TacsRealPart(fail) > TacsRealPart(maxFail)){
            maxFail = fail;
          }
        }
        else {
          // Evaluate the determinant of the Jacobian
          TacsScalar Xd[9], J[9];
          TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);

          if (normType == POWER){
            TacsScalar fp = pow(fabs(fail/maxFail), P);
            failNumer += scale*weight*detJ*(fail/maxFail)*fp;
            failDenom += weight*detJ*fp;
          }
          else if (normType == DISCRETE_POWER){
            TacsScalar fp = pow(fabs(fail/maxFail), P);
            failNumer += scale*(fail/maxFail)*fp;
            failDenom += fp;
          }
          else if (normType == POWER_SQUARED){
            TacsScalar fp = pow(fabs(fail/maxFail), P);
            failNumer += scale*weight*detJ*(fail*fail/maxFail)*fp;
            failDenom += weight*detJ*fp;
          }
          else if (normType == DISCRETE_POWER_SQUARED){
            TacsScalar fp = pow(fabs(fail/maxFail), P);
            failNumer += scale*(fail*fail/maxFail)*fp;
            failDenom += fp;
          }
          else if (normType == EXPONENTIAL){
            TacsScalar efp = exp(P*(fail - maxFail));
            failNumer += scale*weight*detJ*(fail/maxFail)*efp;
            failDenom += weight*detJ*efp;
          }
          else if (normType == DISCRETE_EXPONENTIAL){
            TacsScalar efp = exp(P*(fail - maxFail));
            failNumer += scale*(fail/maxFail)*efp;
            failDenom += efp;
          }
          else if (normType == EXPONENTIAL_SQUARED){
            TacsScalar efp = exp(P*(fail - maxFail));
            failNumer += scale*weight*detJ*(fail*fail/maxFail)*efp;
            failDenom += weight*detJ*efp;
          }
          else if (normType == DISCRETE_EXPONENTIAL_SQUARED){
            TacsScalar efp = exp(P*(fail - maxFail));
            failNumer += scale*(fail*fail/maxFail)*efp;
            failDenom += efp;
          }
        }
      }
    }
  }
}

/*
  Determine the derivative of the P-norm function w.r.t. the state
  variables over this element.
*/
void TACSInducedFailure::getElementSVSens( int elemIndex,
                                           TACSElement *element,
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
      int count = element->evalPointQuantity(elemIndex, TACS_FAILURE_INDEX,
                                             time, i, pt,
                                             Xpts, vars, dvars, ddvars,
                                             &fail);

      if (count >= 1){
        // Evaluate the determinant of the Jacobian
        TacsScalar Xd[9], J[9];
        TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);

        // Compute the derivative of the induced aggregation with
        // respect to the failure function
        TacsScalar dfdq = 0.0;
        if (normType == POWER){
          TacsScalar g = fail/maxFail;
          TacsScalar fp = pow(fabs(g), P);

          dfdq = weight*detJ*((1.0 + P)*g*failDenom -
                              P*failNumer)*fp/(g*failDenom*failDenom);
        }
        else if (normType == DISCRETE_POWER){
          TacsScalar g = fail/maxFail;
          TacsScalar fp = pow(fabs(g), P);

          dfdq = ((1.0 + P)*g*failDenom -
                  P*failNumer)*fp/(g*failDenom*failDenom);
        }
        else if (normType == POWER_SQUARED){
          TacsScalar g = fail/maxFail;
          TacsScalar fp = pow(fabs(g), P);

          dfdq = weight*detJ*((2.0 + P)*fail*g*failDenom -
                              P*failNumer)*fp/(g*failDenom*failDenom);
        }
        else if (normType == DISCRETE_POWER_SQUARED){
          TacsScalar g = fail/maxFail;
          TacsScalar fp = pow(fabs(g), P);

          dfdq = ((2.0 + P)*fail*g*failDenom -
                  P*failNumer)*fp/(g*failDenom*failDenom);
        }
        else if (normType == EXPONENTIAL){
          TacsScalar efp = exp(P*(fail - maxFail));

          dfdq = weight*detJ*(((1.0 + P*fail)*failDenom -
                               P*maxFail*failNumer)*efp)/(failDenom*failDenom);
        }
        else if (normType == DISCRETE_EXPONENTIAL){
          TacsScalar efp = exp(P*(fail - maxFail));

          dfdq = (((1.0 + P*fail)*failDenom -
                   P*maxFail*failNumer)*efp)/(failDenom*failDenom);
        }
        else if (normType == EXPONENTIAL_SQUARED){
          TacsScalar efp = exp(P*(fail - maxFail));

          dfdq = weight*detJ*(((2.0 + P*fail)*fail*failDenom -
                               P*maxFail*failNumer)*efp)/(failDenom*failDenom);
        }
        else if (normType == DISCRETE_EXPONENTIAL_SQUARED){
          TacsScalar efp = exp(P*(fail - maxFail));

          dfdq = (((2.0 + P*fail)*fail*failDenom -
                   P*maxFail*failNumer)*efp)/(failDenom*failDenom);
        }

        element->addPointQuantitySVSens(elemIndex, TACS_FAILURE_INDEX, time,
                                        alpha, beta, gamma,
                                        i, pt, Xpts, vars, dvars, ddvars,
                                        &dfdq, dfdu);
      }
    }
  }
}

/*
  Determine the derivative of the function with respect to the element
  nodal locations
*/
void TACSInducedFailure::getElementXptSens( int elemIndex,
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
      int count = element->evalPointQuantity(elemIndex, TACS_FAILURE_INDEX,
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
        if (normType == POWER){
          TacsScalar g = fail/maxFail;
          TacsScalar fp = pow(fabs(g), P);

          dfdq = weight*detJ*((1.0 + P)*g*failDenom -
                              P*failNumer)*fp/(g*failDenom*failDenom);
          dfddetJ = ((fail*failDenom -
                      maxFail*failNumer)*weight*fp)/(failDenom*failDenom);
        }
        else if (normType == DISCRETE_POWER){
          TacsScalar g = fail/maxFail;
          TacsScalar fp = pow(fabs(g), P);

          dfdq = ((1.0 + P)*g*failDenom -
                  P*failNumer)*fp/(g*failDenom*failDenom);
        }
        else if (normType == POWER_SQUARED){
          TacsScalar g = fail/maxFail;
          TacsScalar fp = pow(fabs(g), P);

          dfdq = weight*detJ*((2.0 + P)*fail*g*failDenom -
                              P*failNumer)*fp/(g*failDenom*failDenom);
          dfddetJ = ((fail*fail*failDenom -
                      maxFail*failNumer)*weight*fp)/(failDenom*failDenom);
        }
        else if (normType == DISCRETE_POWER_SQUARED){
          TacsScalar g = fail/maxFail;
          TacsScalar fp = pow(fabs(g), P);

          dfdq = ((2.0 + P)*fail*g*failDenom -
                  P*failNumer)*fp/(g*failDenom*failDenom);
        }
        else if (normType == EXPONENTIAL){
          TacsScalar efp = exp(P*(fail - maxFail));

          dfdq = weight*detJ*(((1.0 + P*fail)*failDenom -
                               P*maxFail*failNumer)*efp)/(failDenom*failDenom);
          dfddetJ = ((fail*failDenom -
                      maxFail*failNumer)*weight*efp)/(failDenom*failDenom);
        }
        else if (normType == DISCRETE_EXPONENTIAL){
          TacsScalar efp = exp(P*(fail - maxFail));

          dfdq = (((1.0 + P*fail)*failDenom -
                   P*maxFail*failNumer)*efp)/(failDenom*failDenom);
        }
        else if (normType == EXPONENTIAL_SQUARED){
          TacsScalar efp = exp(P*(fail - maxFail));

          dfdq = weight*detJ*(((2.0 + P*fail)*fail*failDenom -
                               P*maxFail*failNumer)*efp)/(failDenom*failDenom);
          dfddetJ = ((fail*fail*failDenom -
                      maxFail*failNumer)*weight*efp)/(failDenom*failDenom);
        }
        else if (normType == DISCRETE_EXPONENTIAL_SQUARED){
          TacsScalar efp = exp(P*(fail - maxFail));

          dfdq = (((2.0 + P*fail)*fail*failDenom -
                   P*maxFail*failNumer)*efp)/(failDenom*failDenom);
        }

        element->addPointQuantityXptSens(elemIndex, TACS_FAILURE_INDEX, time,
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
  Determine the derivative of the function with respect to the design
  variables defined by the element - usually just the
  constitutive/material design variables.
*/
void TACSInducedFailure::addElementDVSens( int elemIndex,
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
      int count = element->evalPointQuantity(elemIndex, TACS_FAILURE_INDEX,
                                             time, i, pt,
                                             Xpts, vars, dvars, ddvars,
                                             &fail);

      if (count >= 1){
        // Evaluate the determinant of the Jacobian
        TacsScalar Xd[9], J[9];
        TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);

        // Compute the derivative of the induced aggregation with
        // respect to the failure function
        TacsScalar dfdq = 0.0;
        if (normType == POWER){
          TacsScalar g = fail/maxFail;
          TacsScalar fp = pow(fabs(g), P);

          dfdq = weight*detJ*((1.0 + P)*g*failDenom -
                              P*failNumer)*fp/(g*failDenom*failDenom);
        }
        else if (normType == DISCRETE_POWER){
          TacsScalar g = fail/maxFail;
          TacsScalar fp = pow(fabs(g), P);

          dfdq = ((1.0 + P)*g*failDenom -
                  P*failNumer)*fp/(g*failDenom*failDenom);
        }
        else if (normType == POWER_SQUARED){
          TacsScalar g = fail/maxFail;
          TacsScalar fp = pow(fabs(g), P);

          dfdq = weight*detJ*((2.0 + P)*fail*g*failDenom -
                              P*failNumer)*fp/(g*failDenom*failDenom);
        }
        else if (normType == DISCRETE_POWER_SQUARED){
          TacsScalar g = fail/maxFail;
          TacsScalar fp = pow(fabs(g), P);

          dfdq = ((2.0 + P)*fail*g*failDenom -
                  P*failNumer)*fp/(g*failDenom*failDenom);
        }
        else if (normType == EXPONENTIAL){
          TacsScalar efp = exp(P*(fail - maxFail));

          dfdq = weight*detJ*(((1.0 + P*fail)*failDenom -
                               P*maxFail*failNumer)*efp)/(failDenom*failDenom);
        }
        else if (normType == DISCRETE_EXPONENTIAL){
          TacsScalar efp = exp(P*(fail - maxFail));

          dfdq = (((1.0 + P*fail)*failDenom -
                   P*maxFail*failNumer)*efp)/(failDenom*failDenom);
        }
        else if (normType == EXPONENTIAL_SQUARED){
          TacsScalar efp = exp(P*(fail - maxFail));

          dfdq = weight*detJ*(((2.0 + P*fail)*fail*failDenom -
                               P*maxFail*failNumer)*efp)/(failDenom*failDenom);
        }
        else if (normType == DISCRETE_EXPONENTIAL_SQUARED){
          TacsScalar efp = exp(P*(fail - maxFail));

          dfdq = (((2.0 + P*fail)*fail*failDenom -
                   P*maxFail*failNumer)*efp)/(failDenom*failDenom);
        }

        element->addPointQuantityDVSens(elemIndex, TACS_FAILURE_INDEX,
                                        time, scale, i, pt,
                                        Xpts, vars, dvars, ddvars,
                                        &dfdq, dvLen, dfdx);
      }
    }
  }
}
