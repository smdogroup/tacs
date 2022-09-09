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

#include "PlaneStressTri6.h"

#include "FElibrary.h"

PlaneStressTri6::PlaneStressTri6(PlaneStressStiffness *_stiff,
                                 ElementBehaviorType type, int _componentNum)
    : TACS2DElement<6>(_stiff, type, _componentNum) {}

PlaneStressTri6::~PlaneStressTri6() {}

const char *PlaneStressTri6::elemName = "PlaneStressTri6";

/*
  Evaluates the shape function at pt = [xi, eta], and its derivatives
*/
void PlaneStressTri6::getShapeFunctions(const double pt[], double N[],
                                        double Na[], double Nb[]) {
  // shape function values
  N[0] = 2.0 * (1.0 - (pt[0] + pt[1])) * (0.5 - (pt[0] + pt[1]));
  N[1] = 2.0 * pt[0] * (pt[0] - 0.5);
  N[2] = 2.0 * pt[1] * (pt[1] - 0.5);
  N[3] = 4.0 * pt[0] * (1.0 - (pt[0] + pt[1]));
  N[4] = 4.0 * pt[0] * pt[1];
  N[5] = 4.0 * pt[1] * (1.0 - (pt[0] + pt[1]));

  // derivative of shape function with respect to xi
  Na[0] = 4.0 * (pt[0] + pt[1] - 0.75);
  Na[1] = 4.0 * pt[0] - 1.0;
  Na[2] = 0.0;
  Na[3] = 4.0 * (1.0 - 2.0 * pt[0] - pt[1]);
  Na[4] = 4.0 * pt[1];
  Na[5] = -4.0 * pt[1];

  // derivative of shape function with respect to eta
  Nb[0] = 4.0 * (pt[0] + pt[1] - 0.75);
  Nb[1] = 0.0;
  Nb[2] = 4.0 * pt[1] - 1.0;
  Nb[3] = -4.0 * pt[0];
  Nb[4] = 4.0 * pt[0];
  Nb[5] = 4.0 * (1.0 - pt[0] - 2.0 * pt[1]);
}

/*
  Evaluates the shape function at pt = [xi, eta]
*/
void PlaneStressTri6::getShapeFunctions(const double pt[], double N[]) {
  // shape function values
  N[0] = 2.0 * (1.0 - (pt[0] + pt[1])) * (0.5 - (pt[0] + pt[1]));
  N[1] = 2.0 * pt[0] * (pt[0] - 0.5);
  N[2] = 2.0 * pt[1] * (pt[1] - 0.5);
  N[3] = 4.0 * pt[0] * (1.0 - (pt[0] + pt[1]));
  N[4] = 4.0 * pt[0] * pt[1];
  N[5] = 4.0 * pt[1] * (1.0 - (pt[0] + pt[1]));
}

/*
  Get the number of Gauss points in the quadrature scheme
*/
int PlaneStressTri6::getNumGaussPts() { return 4; }

/*
  Get the quadrature points
*/
double PlaneStressTri6::getGaussWtsPts(const int num, double pt[]) {
  // Set coordinates of point specified by num input in pt[]
  // Return weight of point as TacsScalar output
  switch (num) {
    case 0:
      pt[0] = 1.0 / 3.0;
      pt[1] = 1.0 / 3.0;
      return -27.0 / 48.0;
    case 1:
      pt[0] = 1.0 / 5.0;
      pt[1] = 3.0 / 5.0;
      return 25.0 / 48.0;
    case 2:
      pt[0] = 1.0 / 5.0;
      pt[1] = 1.0 / 5.0;
      return 25.0 / 48.0;
    case 3:
      pt[0] = 3.0 / 5.0;
      pt[1] = 1.0 / 5.0;
      return 25.0 / 48.0;
    default:
      break;
  }

  return 0.0;
}

/*
  Get the number of elemens/nodes and CSR size of the contributed by
  this element.
*/
void PlaneStressTri6::addOutputCount(int *nelems, int *nnodes, int *ncsr) {
  *nelems += 3;
  *nnodes += 6;
  *ncsr += 12;
}

/*
  Get the output data from this element and place it in a real
  array for visualization later. The values generated for visualization
  are determined by a bit-wise selection variable 'out_type' which is
  can be used to simultaneously write out different data. Note that this
  is why the bitwise operation & is used below.

  The output may consist of the following:
  - the nodal locations
  - the displacements and rotations
  - the strains or strains within the element
  - extra variables that are used for optimization

  output:
  data:     the data to write to the file (eventually)

  input:
  out_type: the bit-wise variable used to specify what data to generate
  vars:     the element variables
  Xpts:     the element nodal locations
*/
void PlaneStressTri6::getOutputData(unsigned int out_type, double *data,
                                    int ld_data, const TacsScalar Xpts[],
                                    const TacsScalar vars[]) {
  // Set the nodal parametric coordinates
  double pt[][2] = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
                    {0.5, 0.0}, {0.5, 0.5}, {0.0, 0.5}};

  for (int n = 0; n < 6; n++) {
    int index = 0;
    if (out_type & TACSElement::OUTPUT_NODES) {
      for (int k = 0; k < 3; k++) {
        data[index + k] = TacsRealPart(Xpts[3 * n + k]);
      }
      index += 3;
    }
    if (out_type & TACSElement::OUTPUT_DISPLACEMENTS) {
      for (int k = 0; k < 2; k++) {
        data[index + k] = TacsRealPart(vars[2 * n + k]);
      }
      index += 2;
    }

    // Compute the shape functions
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES];
    getShapeFunctions(pt[n], N, Na, Nb);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[4];
    this->planeJacobian(X, Xa, N, Na, Nb, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[4];
    FElibrary::jacobian2d(Xa, J);

    // Compute the strain
    TacsScalar strain[3];
    this->evalStrain(strain, J, Na, Nb, vars);

    if (out_type & TACSElement::OUTPUT_STRAINS) {
      for (int k = 0; k < 3; k++) {
        data[index + k] = TacsRealPart(strain[k]);
      }
      index += 3;
    }
    if (out_type & TACSElement::OUTPUT_STRESSES) {
      // Calculate the strain at the current point
      TacsScalar stress[3];
      this->stiff->calculateStress(pt[n], strain, stress);

      for (int k = 0; k < 3; k++) {
        data[index + k] = TacsRealPart(stress[k]);
      }
      index += 3;
    }
    if (out_type & TACSElement::OUTPUT_EXTRAS) {
      // Compute the failure value
      TacsScalar lambda;
      this->stiff->failure(pt[n], strain, &lambda);
      data[index] = TacsRealPart(lambda);

      this->stiff->buckling(strain, &lambda);
      data[index + 1] = TacsRealPart(lambda);

      data[index + 2] = TacsRealPart(this->stiff->getDVOutputValue(0, pt[n]));
      data[index + 3] = TacsRealPart(this->stiff->getDVOutputValue(1, pt[n]));

      index += this->NUM_EXTRAS;
    }

    data += ld_data;
  }
}

/*
  Get the element connectivity for visualization purposes. Since each
  element consists of a series of sub-elements used for visualization,
  we also need the connectivity of these visualization elements.

  output:
  con:  the connectivity of the local visualization elements contributed
  by this finite-element

  input:
  node:  the node offset number - so that this connectivity is more or
  less global
*/
void PlaneStressTri6::getOutputConnectivity(int *con, int node) {
  con[0] = node;
  con[1] = node + 3;
  con[2] = node + 4;
  con[3] = node + 5;
  con += 4;

  con[0] = node + 3;
  con[1] = node + 1;
  con[2] = node + 4;
  con[3] = node + 4;
  con += 4;

  con[0] = node + 5;
  con[1] = node + 4;
  con[2] = node + 2;
  con[3] = node + 2;
}
