/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2019 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TACSHeatFlux.h"

#include "TACSAssembler.h"
#include "TACSElementAlgebra.h"

TACSHeatFlux::TACSHeatFlux(TACSAssembler *_assembler, int *elem_index,
                           int *face_index, int num_elems)
    : TACSFunction(_assembler) {
  // Set the element indices
  for (int i = 0; i < num_elems; i++) {
    int value = 0;
    if (face_index[i] >= 0 && face_index[i] < MAX_SURFACE_INDEX) {
      value = 1 << face_index[i];
    }

    std::map<int, int>::iterator it = element_to_face_key.find(elem_index[i]);

    // Perform a bitwise or operation if the element has already
    // been added to the domain
    if (it != element_to_face_key.end()) {
      value = it->second | value;
    }

    // Set the value in the map
    element_to_face_key[elem_index[i]] = value;
  }

  // Set the domain of the function
  setDomain(num_elems, elem_index);
}

TACSHeatFlux::~TACSHeatFlux() {}

/*
  HeatFluxIntegral function name
*/
const char *TACSHeatFlux::funcName = "TACSHeatFlux";

/*
  Return the function name
*/
const char *TACSHeatFlux::getObjectName() { return funcName; }

/*
  Retrieve the function value
*/
TacsScalar TACSHeatFlux::getFunctionValue() { return value; }

/*
  Reduce the function values across all MPI processes
*/
void TACSHeatFlux::initEvaluation(EvaluationType ftype) { value = 0.0; }

/*
  Reduce the function values across all MPI processes
*/
void TACSHeatFlux::finalEvaluation(EvaluationType ftype) {
  TacsScalar temp = value;
  MPI_Allreduce(&temp, &value, 1, TACS_MPI_TYPE, MPI_SUM,
                assembler->getMPIComm());
}

/*
  Perform the element-wise evaluation of the TACSDisplacementIntegral function.
*/
void TACSHeatFlux::elementWiseEval(EvaluationType ftype, int elemIndex,
                                   TACSElement *element, double time,
                                   TacsScalar scale, const TacsScalar Xpts[],
                                   const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar ddvars[]) {
  // Retrieve the number of stress components for this element
  TACSElementBasis *basis = element->getElementBasis();

  if (basis) {
    // Get the surface index
    int face_key = element_to_face_key[elemIndex];

    for (int face = 0; (face_key && face < MAX_SURFACE_INDEX); face++) {
      // Check if this is a surface that we need to integrate
      if (1 << face & face_key) {
        // Clear the bit from the surface index we just integrated
        face_key &= ~(1 << face);

        for (int i = 0; i < basis->getNumFaceQuadraturePoints(face); i++) {
          double pt[3], tangents[6];
          double weight = basis->getFaceQuadraturePoint(face, i, pt, tangents);

          // Evaluate the heat flux at the quadrature point
          TacsScalar flux[3], detXd = 0.0;
          const int not_a_quadrature_pt = -1;
          int count = element->evalPointQuantity(
              elemIndex, TACS_HEAT_FLUX, time, not_a_quadrature_pt, pt, Xpts,
              vars, dvars, ddvars, &detXd, flux);

          // Compute the component of the flux normal to the surface
          if (count > 0) {
            TacsScalar X[3], Xd[9], normal[3];
            TacsScalar Area =
                basis->getFaceNormal(face, i, Xpts, X, Xd, normal);
            if (count == 2) {
              value += scale * weight * Area * vec2Dot(flux, normal);
            } else if (count == 3) {
              value += scale * weight * Area * vec3Dot(flux, normal);
            }
          }
        }
      }
    }
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.
*/
void TACSHeatFlux::getElementSVSens(
    int elemIndex, TACSElement *element, double time, TacsScalar alpha,
    TacsScalar beta, TacsScalar gamma, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar dfdu[]) {
  // Zero the derivative of the function w.r.t. the element state
  // variables
  int numVars = element->getNumVariables();
  memset(dfdu, 0, numVars * sizeof(TacsScalar));

  // Get the element basis class
  TACSElementBasis *basis = element->getElementBasis();

  if (basis) {
    // Get the surface index
    int face_key = element_to_face_key[elemIndex];

    for (int face = 0; (face_key && face < MAX_SURFACE_INDEX); face++) {
      // Check if this is a surface that we need to integrate
      if (1 << face & face_key) {
        // Clear the bit from the surface index we just integrated
        face_key &= ~(1 << face);

        for (int i = 0; i < basis->getNumFaceQuadraturePoints(face); i++) {
          double pt[3], tangents[6];
          double weight = basis->getFaceQuadraturePoint(face, i, pt, tangents);

          // Evaluate the heat flux at the quadrature point
          TacsScalar flux[3], detXd = 0.0;
          ;
          const int not_a_quadrature_pt = -1;
          int count = element->evalPointQuantity(
              elemIndex, TACS_HEAT_FLUX, time, not_a_quadrature_pt, pt, Xpts,
              vars, dvars, ddvars, &detXd, flux);

          if (count > 0) {
            TacsScalar X[3], Xd[9], normal[3];
            TacsScalar Area =
                basis->getFaceNormal(face, i, Xpts, X, Xd, normal);

            TacsScalar dfdq[3] = {0.0, 0.0, 0.0};
            if (count == 2) {
              dfdq[0] = weight * Area * normal[0];
              dfdq[1] = weight * Area * normal[1];
            } else if (count == 3) {
              dfdq[0] = weight * Area * normal[0];
              dfdq[1] = weight * Area * normal[1];
              dfdq[2] = weight * Area * normal[2];
            }

            element->addPointQuantitySVSens(
                elemIndex, TACS_HEAT_FLUX, time, alpha, beta, gamma,
                not_a_quadrature_pt, pt, Xpts, vars, dvars, ddvars, dfdq, dfdu);
          }
        }
      }
    }
  }
}

/*
  Determine the derivative of the function with respect to
  the element nodal locations
*/
void TACSHeatFlux::getElementXptSens(
    int elemIndex, TACSElement *element, double time, TacsScalar scale,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar dfdXpts[]) {
  // Zero the derivative of the function w.r.t. the element node
  // locations
  int numNodes = element->getNumNodes();
  memset(dfdXpts, 0, 3 * numNodes * sizeof(TacsScalar));

  // Get the element basis class
  TACSElementBasis *basis = element->getElementBasis();

  if (basis) {
    // Get the surface index
    int face_key = element_to_face_key[elemIndex];

    for (int face = 0; (face_key && face < MAX_SURFACE_INDEX); face++) {
      // Check if this is a surface that we need to integrate
      if (1 << face & face_key) {
        // Clear the bit from the surface index we just integrated
        face_key &= ~(1 << face);

        for (int i = 0; i < basis->getNumFaceQuadraturePoints(face); i++) {
          double pt[3], tangents[6];
          double weight = basis->getFaceQuadraturePoint(face, i, pt, tangents);

          // Evaluate the heat flux at the quadrature point
          TacsScalar flux[3], detXd = 0.0;
          const int not_a_quadrature_pt = -1;
          int count = element->evalPointQuantity(
              elemIndex, TACS_HEAT_FLUX, time, not_a_quadrature_pt, pt, Xpts,
              vars, dvars, ddvars, &detXd, flux);

          if (count > 0) {
            TacsScalar X[3], Xd[9], normal[3];
            TacsScalar Area =
                basis->getFaceNormal(face, i, Xpts, X, Xd, normal);

            TacsScalar dfdq[3] = {0.0, 0.0, 0.0};
            TacsScalar dfdn[3] = {0.0, 0.0, 0.0};
            TacsScalar dfdA = 0.0;
            if (count == 2) {
              dfdq[0] = weight * Area * normal[0];
              dfdq[1] = weight * Area * normal[1];

              dfdn[0] = scale * weight * Area * flux[0];
              dfdn[1] = scale * weight * Area * flux[1];

              dfdA = scale * weight * vec2Dot(flux, normal);
            } else if (count == 3) {
              dfdq[0] = weight * Area * normal[0];
              dfdq[1] = weight * Area * normal[1];
              dfdq[2] = weight * Area * normal[2];

              dfdn[0] = scale * weight * Area * flux[0];
              dfdn[1] = scale * weight * Area * flux[1];
              dfdn[2] = scale * weight * Area * flux[2];

              dfdA = scale * weight * vec3Dot(flux, normal);
            }
            TacsScalar dfddetXd = 0.0;
            element->addPointQuantityXptSens(
                elemIndex, TACS_HEAT_FLUX, time, scale, not_a_quadrature_pt, pt,
                Xpts, vars, dvars, ddvars, dfddetXd, dfdq, dfdXpts);

            basis->addFaceNormalXptSens(face, i, Area, Xd, normal, dfdA, NULL,
                                        NULL, dfdn, dfdXpts);
          }
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
void TACSHeatFlux::addElementDVSens(
    int elemIndex, TACSElement *element, double time, TacsScalar scale,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int dvLen, TacsScalar dfdx[]) {
  // Get the element basis class
  TACSElementBasis *basis = element->getElementBasis();

  if (basis) {
    // Get the surface index
    int face_key = element_to_face_key[elemIndex];

    for (int face = 0; (face_key && face < MAX_SURFACE_INDEX); face++) {
      // Check if this is a surface that we need to integrate
      if (1 << face & face_key) {
        // Clear the bit from the surface index we just integrated
        face_key &= ~(1 << face);

        for (int i = 0; i < basis->getNumFaceQuadraturePoints(face); i++) {
          double pt[3], tangents[6];
          double weight = basis->getFaceQuadraturePoint(face, i, pt, tangents);

          // Evaluate the heat flux at the quadrature point
          TacsScalar flux[3], detXd = 0.0;
          const int not_a_quadrature_pt = -1;
          int count = element->evalPointQuantity(
              elemIndex, TACS_HEAT_FLUX, time, not_a_quadrature_pt, pt, Xpts,
              vars, dvars, ddvars, &detXd, flux);

          if (count > 0) {
            TacsScalar X[3], Xd[9], normal[3];
            TacsScalar Area =
                basis->getFaceNormal(face, i, Xpts, X, Xd, normal);

            TacsScalar dfdq[3] = {0.0, 0.0, 0.0};
            if (count == 2) {
              dfdq[0] = weight * Area * normal[0];
              dfdq[1] = weight * Area * normal[1];
            } else if (count == 3) {
              dfdq[0] = weight * Area * normal[0];
              dfdq[1] = weight * Area * normal[1];
              dfdq[2] = weight * Area * normal[2];
            }

            element->addPointQuantityDVSens(
                elemIndex, TACS_HEAT_FLUX, time, scale, not_a_quadrature_pt, pt,
                Xpts, vars, dvars, ddvars, dfdq, dvLen, dfdx);
          }
        }
      }
    }
  }
}
