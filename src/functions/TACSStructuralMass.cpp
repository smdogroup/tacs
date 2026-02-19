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

#include "TACSStructuralMass.h"

#include "TACSAssembler.h"

/*
  Allocate the structural mass function
*/
TACSStructuralMass::TACSStructuralMass(TACSAssembler *_assembler)
    : TACSFunction(_assembler) {
  totalMass = 0.0;
}

/*
  Destructor for the structural mass
*/
TACSStructuralMass::~TACSStructuralMass() {}

const char *TACSStructuralMass::funcName = "StructuralMass";

/*
  The structural mass function name
*/
const char *TACSStructuralMass::getObjectName() { return funcName; }

/*
  Get the function name
*/
TacsScalar TACSStructuralMass::getFunctionValue() { return totalMass; }

/*
  Initialize the mass to zero
*/
void TACSStructuralMass::initEvaluation(EvaluationType ftype) {
  totalMass = 0.0;
}

/*
  Sum the mass across all MPI processes
*/
void TACSStructuralMass::finalEvaluation(EvaluationType ftype) {
  TacsScalar temp = totalMass;
  MPI_Allreduce(&temp, &totalMass, 1, TACS_MPI_TYPE, MPI_SUM,
                assembler->getMPIComm());
}

/*
  Perform the element-wise evaluation of the TACSKSFailure function.
*/
void TACSStructuralMass::elementWiseEval(
    EvaluationType ftype, int elemIndex, TACSElement *element, double time,
    TacsScalar scale, const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[]) {
  for (int i = 0; i < element->getNumQuadraturePoints(); i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    // Evaluate the failure index, and check whether it is an
    // undefined quantity of interest on this element
    TacsScalar density = 0.0, detXd = 0.0;
    int count =
        element->evalPointQuantity(elemIndex, TACS_ELEMENT_DENSITY, time, i, pt,
                                   Xpts, vars, dvars, ddvars, &detXd, &density);

    if (count >= 1) {
      totalMass += scale * weight * detXd * density;
    }
  }
}

/*
  Determine the derivative of the mass w.r.t. the material
  design variables
*/
void TACSStructuralMass::addElementDVSens(
    int elemIndex, TACSElement *element, double time, TacsScalar scale,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int dvLen, TacsScalar dfdx[]) {
  for (int i = 0; i < element->getNumQuadraturePoints(); i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    TacsScalar density = 0.0, detXd = 0.0;
    int count =
        element->evalPointQuantity(elemIndex, TACS_ELEMENT_DENSITY, time, i, pt,
                                   Xpts, vars, dvars, ddvars, &detXd, &density);

    if (count >= 1) {
      TacsScalar dfdq = scale * weight * detXd;
      element->addPointQuantityDVSens(elemIndex, TACS_ELEMENT_DENSITY, time,
                                      1.0, i, pt, Xpts, vars, dvars, ddvars,
                                      &dfdq, dvLen, dfdx);
    }
  }
}

/*
  Determine the derivative of the mass w.r.t. the element nodal
  locations.
*/
void TACSStructuralMass::getElementXptSens(
    int elemIndex, TACSElement *element, double time, TacsScalar scale,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar dfdXpts[]) {
  // Zero the derivative of the function w.r.t. the node locations
  int numNodes = element->getNumNodes();
  memset(dfdXpts, 0, 3 * numNodes * sizeof(TacsScalar));

  for (int i = 0; i < element->getNumQuadraturePoints(); i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    TacsScalar density = 0.0, detXd = 0.0;
    int count =
        element->evalPointQuantity(elemIndex, TACS_ELEMENT_DENSITY, time, i, pt,
                                   Xpts, vars, dvars, ddvars, &detXd, &density);

    if (count >= 1) {
      TacsScalar dfdq = scale * weight * detXd;
      TacsScalar dfddetXd = scale * weight * density;
      element->addPointQuantityXptSens(elemIndex, TACS_ELEMENT_DENSITY, time,
                                       1.0, i, pt, Xpts, vars, dvars, ddvars,
                                       dfddetXd, &dfdq, dfdXpts);
    }
  }
}
