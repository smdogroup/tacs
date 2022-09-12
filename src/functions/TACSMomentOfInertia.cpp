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

#include "TACSMomentOfInertia.h"

#include "TACSAssembler.h"
#include "TACSElementAlgebra.h"

/*
  Allocate the structural mass function
*/
TACSMomentOfInertia::TACSMomentOfInertia(TACSAssembler *_assembler,
                                         const double _dir1[],
                                         const double _dir2[], int _cmFlag)
    : TACSFunction(_assembler) {
  totalMass = 0.0;
  massMoment[0] = massMoment[1] = massMoment[2] = 0.0;
  I0 = 0.0;

  dir1[0] = TacsScalar(_dir1[0]);
  dir1[1] = TacsScalar(_dir1[1]);
  dir1[2] = TacsScalar(_dir1[2]);

  dir2[0] = TacsScalar(_dir2[0]);
  dir2[1] = TacsScalar(_dir2[1]);
  dir2[2] = TacsScalar(_dir2[2]);

  cmFlag = _cmFlag;
}

/*
  Destructor for the structural mass
*/
TACSMomentOfInertia::~TACSMomentOfInertia() {}

const char *TACSMomentOfInertia::funcName = "MomentOfInertia";

/*
  The structural mass function name
*/
const char *TACSMomentOfInertia::getObjectName() { return funcName; }

/*
  Get the function name
*/
TacsScalar TACSMomentOfInertia::getFunctionValue() {
  // Compute moment of inertia about origin
  TacsScalar value = I0;

  // Use reverse parallel axis theorem to move moment of inertia to cm
  if (cmFlag) {
    TacsScalar term2 = vec3Dot(massMoment, massMoment) * vec3Dot(dir1, dir2) -
                       vec3Dot(massMoment, dir1) * vec3Dot(massMoment, dir2);
    value -= term2 / totalMass;
  }

  return value;
}

/*
  Initialize the mass/moments to zero
*/
void TACSMomentOfInertia::initEvaluation(EvaluationType ftype) {
  totalMass = I0 = 0.0;
  massMoment[0] = massMoment[1] = massMoment[2] = 0.0;
}

/*
  Sum the mass/moments across all MPI processes
*/
void TACSMomentOfInertia::finalEvaluation(EvaluationType ftype) {
  TacsScalar temp = totalMass;
  MPI_Allreduce(&temp, &totalMass, 1, TACS_MPI_TYPE, MPI_SUM,
                assembler->getMPIComm());
  TacsScalar temp2[3] = {massMoment[0], massMoment[1], massMoment[2]};
  MPI_Allreduce(temp2, massMoment, 3, TACS_MPI_TYPE, MPI_SUM,
                assembler->getMPIComm());
  temp = I0;
  MPI_Allreduce(&temp, &I0, 1, TACS_MPI_TYPE, MPI_SUM, assembler->getMPIComm());
}

/*
  Perform the element-wise evaluation of the TACSMomentOfInertia function.
*/
void TACSMomentOfInertia::elementWiseEval(
    EvaluationType ftype, int elemIndex, TACSElement *element, double time,
    TacsScalar scale, const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[]) {
  for (int i = 0; i < element->getNumQuadraturePoints(); i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    // Evaluate parallel axis theorem portion of function
    TacsScalar density = 0.0, detXd = 0.0;
    int count;
    if (cmFlag) {
      count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DENSITY, time,
                                         i, pt, Xpts, vars, dvars, ddvars,
                                         &detXd, &density);

      if (count >= 1) {
        totalMass += scale * weight * detXd * density;
      }

      TacsScalar densityMoment[3];
      count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DENSITY_MOMENT,
                                         time, i, pt, Xpts, vars, dvars, ddvars,
                                         &detXd, densityMoment);

      for (int j = 0; j < count; j++) {
        massMoment[j] += scale * weight * detXd * densityMoment[j];
      }
    }

    // Evaluate moment of inertia about origin
    TacsScalar I0_elem[6];
    count = element->evalPointQuantity(
        elemIndex, TACS_ELEMENT_MOMENT_OF_INERTIA, time, i, pt, Xpts, vars,
        dvars, ddvars, &detXd, I0_elem);

    // Compute inner product of I0 (i.e. v1^T.I0.v2)
    TacsScalar ip[6];
    getInnerProductFactor(count, ip);
    for (int j = 0; j < count; j++) {
      I0 += scale * weight * detXd * I0_elem[j] * ip[j];
    }
  }
}

/*
  Determine the derivative of the mass w.r.t. the material
  design variables
*/
void TACSMomentOfInertia::addElementDVSens(
    int elemIndex, TACSElement *element, double time, TacsScalar scale,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int dvLen, TacsScalar dfdx[]) {
  for (int i = 0; i < element->getNumQuadraturePoints(); i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    TacsScalar density = 0.0, detXd = 0.0;
    int count;
    if (cmFlag) {
      count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DENSITY, time,
                                         i, pt, Xpts, vars, dvars, ddvars,
                                         &detXd, &density);

      if (count >= 1) {
        TacsScalar term2 =
            vec3Dot(massMoment, massMoment) * vec3Dot(dir1, dir2) -
            vec3Dot(massMoment, dir1) * vec3Dot(massMoment, dir2);
        TacsScalar dfdq =
            term2 / totalMass / totalMass * scale * weight * detXd;
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
          TacsScalar termj = 2.0 * vec3Dot(dir1, dir2) * massMoment[j] -
                             vec3Dot(massMoment, dir1) * dir2[j] -
                             vec3Dot(massMoment, dir2) * dir1[j];
          dfdq[j] = -scale * weight * detXd / totalMass * termj;
        }

        element->addPointQuantityDVSens(elemIndex, TACS_ELEMENT_DENSITY_MOMENT,
                                        time, 1.0, i, pt, Xpts, vars, dvars,
                                        ddvars, dfdq, dvLen, dfdx);
      }
    }

    TacsScalar I0_elem[6];
    count = element->evalPointQuantity(
        elemIndex, TACS_ELEMENT_MOMENT_OF_INERTIA, time, i, pt, Xpts, vars,
        dvars, ddvars, &detXd, I0_elem);

    if (count >= 1) {
      TacsScalar dfdq[6] = {0.0};
      TacsScalar ip[6];
      getInnerProductFactor(count, ip);
      for (int j = 0; j < count; j++) {
        dfdq[j] = scale * weight * detXd * ip[j];
      }

      element->addPointQuantityDVSens(elemIndex, TACS_ELEMENT_MOMENT_OF_INERTIA,
                                      time, 1.0, i, pt, Xpts, vars, dvars,
                                      ddvars, dfdq, dvLen, dfdx);
    }
  }
}

/*
  Determine the derivative of the mass w.r.t. the element nodal
  locations.
*/
void TACSMomentOfInertia::getElementXptSens(
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
    int count;
    if (cmFlag) {
      count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DENSITY, time,
                                         i, pt, Xpts, vars, dvars, ddvars,
                                         &detXd, &density);

      if (count >= 1) {
        TacsScalar term2 =
            vec3Dot(massMoment, massMoment) * vec3Dot(dir1, dir2) -
            vec3Dot(massMoment, dir1) * vec3Dot(massMoment, dir2);
        TacsScalar dfdq =
            term2 / totalMass / totalMass * scale * weight * detXd;
        TacsScalar dfddetXd =
            term2 / totalMass / totalMass * scale * weight * density;
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
          TacsScalar termj = 2.0 * vec3Dot(dir1, dir2) * massMoment[j] -
                             vec3Dot(massMoment, dir1) * dir2[j] -
                             vec3Dot(massMoment, dir2) * dir1[j];
          dfdq[j] = -scale * weight * detXd / totalMass * termj;
          dfddetXd += -scale * weight / totalMass * termj * densityMoment[j];
        }
        element->addPointQuantityXptSens(elemIndex, TACS_ELEMENT_DENSITY_MOMENT,
                                         time, 1.0, i, pt, Xpts, vars, dvars,
                                         ddvars, dfddetXd, dfdq, dfdXpts);
      }
    }

    TacsScalar I0_elem[6];
    count = element->evalPointQuantity(
        elemIndex, TACS_ELEMENT_MOMENT_OF_INERTIA, time, i, pt, Xpts, vars,
        dvars, ddvars, &detXd, I0_elem);

    if (count >= 1) {
      TacsScalar dfdq[6] = {0.0};
      TacsScalar dfddetXd = 0.0;
      TacsScalar ip[6];
      getInnerProductFactor(count, ip);
      for (int j = 0; j < count; j++) {
        dfdq[j] = scale * weight * detXd * ip[j];
        dfddetXd += scale * weight * ip[j] * I0_elem[j];
      }

      element->addPointQuantityXptSens(
          elemIndex, TACS_ELEMENT_MOMENT_OF_INERTIA, time, 1.0, i, pt, Xpts,
          vars, dvars, ddvars, dfddetXd, dfdq, dfdXpts);
    }
  }
}

/*
  Get 2D/3D inner product factors
*/
void TACSMomentOfInertia::getInnerProductFactor(int count, TacsScalar ip[]) {
  if (count == 3) {  // 2D inner product
    ip[0] = dir1[0] * dir2[0];
    ip[1] = dir1[0] * dir2[1] + dir1[1] * dir2[0];
    ip[2] = dir1[1] * dir2[1];
  } else {  // 3D inner product
    ip[0] = dir1[0] * dir2[0];
    ip[1] = dir1[0] * dir2[1] + dir1[1] * dir2[0];
    ip[2] = dir1[0] * dir2[2] + dir1[2] * dir2[0];
    ip[3] = dir1[1] * dir2[1];
    ip[4] = dir1[1] * dir2[2] + dir1[2] * dir2[1];
    ip[5] = dir1[2] * dir2[2];
  }
}
