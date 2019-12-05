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
TACSStructuralMass::TACSStructuralMass( TACSAssembler *_assembler ):
TACSFunction(_assembler){
  totalMass = 0.0;
}

/*
  Destructor for the structural mass
*/
TACSStructuralMass::~TACSStructuralMass(){}

const char *TACSStructuralMass::funcName = "StructuralMass";

/*
  The structural mass function name
*/
const char* TACSStructuralMass::getObjectName(){
  return funcName;
}

/*
  Get the function name
*/
TacsScalar TACSStructuralMass::getFunctionValue(){
  return totalMass;
}

/*
  Initialize the mass to zero
*/
void TACSStructuralMass::initEvaluation( EvaluationType ftype ){
  totalMass = 0.0;
}

/*
  Sum the mass across all MPI processes
*/
void TACSStructuralMass::finalEvaluation( EvaluationType ftype ){
  TacsScalar temp = totalMass;
  MPI_Allreduce(&temp, &totalMass, 1, TACS_MPI_TYPE,
                MPI_SUM, assembler->getMPIComm());
}

/*
  Perform the element-wise evaluation of the TACSKSFailure function.
*/
void TACSStructuralMass::elementWiseEval( EvaluationType ftype,
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
      TacsScalar density = 0.0;
      int count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DENSITY,
                                             time, i, pt,
                                             Xpts, vars, dvars, ddvars,
                                             &density);

      if (count >= 1){
        // Evaluate the determinant of the Jacobian
        TacsScalar Xd[9], J[9];
        TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);
        totalMass += scale*weight*detJ*density;
      }
    }
  }
}

/*
  Determine the derivative of the mass w.r.t. the element nodal
  locations.
*/
void TACSStructuralMass::getElementXptSens( int elemIndex,
                                            TACSElement *element,
                                            double time,
                                            TacsScalar scale,
                                            const TacsScalar Xpts[],
                                            const TacsScalar vars[],
                                            const TacsScalar dvars[],
                                            const TacsScalar ddvars[],
                                            TacsScalar dfdXpts[] ){
  // Zero the derivative of the function w.r.t. the node locations
  int numNodes = element->getNumNodes();
  memset(dfdXpts, 0, 3*numNodes*sizeof(TacsScalar));

  // Get the element basis class
  TACSElementBasis *basis = element->getElementBasis();

  if (basis){
    for ( int i = 0; i < basis->getNumQuadraturePoints(); i++ ){
      double pt[3];
      double weight = basis->getQuadraturePoint(i, pt);

      TacsScalar density = 0.0;
      int count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DENSITY,
                                             time, i, pt,
                                             Xpts, vars, dvars, ddvars,
                                             &density);

      if (count >= 1){
        // Evaluate the determinant of the Jacobian
        TacsScalar Xd[9], J[9];
        basis->getJacobianTransform(pt, Xpts, Xd, J);

        // Compute the sensitivity contribution
        TacsScalar dfddetJ = density*weight;
        basis->addJacobianTransformXptSens(pt, Xd, J, scale*dfddetJ,
                                           NULL, NULL, dfdXpts);
      }
    }
  }
}

/*
  Determine the derivative of the mass w.r.t. the material
  design variables
*/
void TACSStructuralMass::addElementDVSens( int elemIndex,
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

      TacsScalar density = 0.0;
      int count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DENSITY,
                                             time, i, pt,
                                             Xpts, vars, dvars, ddvars,
                                             &density);
      if (count >= 1){
        // Evaluate the determinant of the Jacobian
        TacsScalar Xd[9], J[9];
        TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);
        TacsScalar dfdq = weight*detJ;

        element->addPointQuantityDVSens(elemIndex, TACS_ELEMENT_DENSITY,
                                        time, scale, i, pt,
                                        Xpts, vars, dvars, ddvars,
                                        &dfdq, dvLen, dfdx);
      }
    }
  }
}
