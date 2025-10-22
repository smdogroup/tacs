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

#include "TACSCenterOfMass.h"

#include "TACSAssembler.h"

/*
  Allocate the structural mass function
*/
TACSCenterOfMass::TACSCenterOfMass(TACSAssembler *_assembler,
                                   const double _dir[])
    : TACSFunction(_assembler) {
  totalMass = 0.0;
  massMoment = 0.0;

  dir[0] = _dir[0];
  dir[1] = _dir[1];
  dir[2] = _dir[2];
}

/*
  Destructor for the structural mass
*/
TACSCenterOfMass::~TACSCenterOfMass() {}

const char *TACSCenterOfMass::funcName = "CenterOfMass";

/*
  The structural mass function name
*/
const char *TACSCenterOfMass::getObjectName() { return funcName; }

/*
  Get the function name
*/
TacsScalar TACSCenterOfMass::getFunctionValue() {
  return massMoment / totalMass;
}

/*
  Initialize the mass to zero
*/
void TACSCenterOfMass::initEvaluation(EvaluationType ftype) {
  totalMass = 0.0;
  massMoment = 0.0;
}

/*
  Sum the mass across all MPI processes
*/
void TACSCenterOfMass::finalEvaluation(EvaluationType ftype) {
  TacsScalar temp = totalMass;
  MPI_Allreduce(&temp, &totalMass, 1, TACS_MPI_TYPE, MPI_SUM,
                assembler->getMPIComm());
  temp = massMoment;
  MPI_Allreduce(&temp, &massMoment, 1, TACS_MPI_TYPE, MPI_SUM,
                assembler->getMPIComm());
}

/*
  Perform the element-wise evaluation of the TACSKSFailure function.
*/
void TACSCenterOfMass::elementWiseEval(
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
    TacsScalar densityMoment[3];
    count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DENSITY_MOMENT,
                                       time, i, pt, Xpts, vars, dvars, ddvars,
                                       &detXd, densityMoment);

    for (int j = 0; j < count; j++) {
      massMoment += scale * weight * detXd * densityMoment[j] * dir[j];
    }
  }
}

/*
  Determine the derivative of the mass w.r.t. the material
  design variables
*/
void TACSCenterOfMass::addElementDVSens(
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
      TacsScalar dfdq =
          -massMoment / totalMass / totalMass * scale * weight * detXd;
      element->addPointQuantityDVSens(elemIndex, TACS_ELEMENT_DENSITY, time,
                                      1.0, i, pt, Xpts, vars, dvars, ddvars,
                                      &dfdq, dvLen, dfdx);
    }

    TacsScalar densityMoment[3];
    count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DENSITY_MOMENT,
                                       time, i, pt, Xpts, vars, dvars, ddvars,
                                       &detXd, densityMoment);

    if (count >= 1) {
      TacsScalar dfdq[3] = {0.0};
      for (int j = 0; j < count; j++) {
        dfdq[j] = scale * weight * detXd * dir[j] / totalMass;
      }

      element->addPointQuantityDVSens(elemIndex, TACS_ELEMENT_DENSITY_MOMENT,
                                      time, 1.0, i, pt, Xpts, vars, dvars,
                                      ddvars, dfdq, dvLen, dfdx);
    }
  }
}

/*
  Determine the derivative of the mass w.r.t. the element nodal
  locations.
*/
void TACSCenterOfMass::getElementXptSens(
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
      TacsScalar dfdq =
          -massMoment / totalMass / totalMass * scale * weight * detXd;
      TacsScalar dfddetXd =
          -massMoment / totalMass / totalMass * scale * weight * density;
      element->addPointQuantityXptSens(elemIndex, TACS_ELEMENT_DENSITY, time,
                                       1.0, i, pt, Xpts, vars, dvars, ddvars,
                                       dfddetXd, &dfdq, dfdXpts);
    }

    TacsScalar densityMoment[3];
    count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DENSITY_MOMENT,
                                       time, i, pt, Xpts, vars, dvars, ddvars,
                                       &detXd, densityMoment);

    if (count >= 1) {
      TacsScalar dfdq[3] = {0.0};
      TacsScalar dfddetXd = 0.0;
      for (int j = 0; j < count; j++) {
        dfdq[j] = scale * weight * detXd * dir[j] / totalMass;
        dfddetXd += scale * weight * densityMoment[j] * dir[j] / totalMass;
      }
      element->addPointQuantityXptSens(elemIndex, TACS_ELEMENT_DENSITY_MOMENT,
                                       time, 1.0, i, pt, Xpts, vars, dvars,
                                       ddvars, dfddetXd, dfdq, dfdXpts);
    }
  }
}
