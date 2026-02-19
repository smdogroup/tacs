#ifndef TACS_PLANE_STRESS_COUPLED_THERMO_MULTI_QUAD_H
#define TACS_PLANE_STRESS_COUPLED_THERMO_MULTI_QUAD_H

/*
  Plane stress element with coupled thermoelastic analysis

  Copyright (c) 2010-2015 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.

  The following code uses templates to allow for arbitrary order elements.
*/
#include "FElibrary.h"
#include "TACS2DCoupledThermoElement.h"
#include "TACSElement.h"

template <int order>
class PlaneStressCoupledThermoQuad
    : public TACS2DCoupledThermoElement<order * order> {
 public:
  PlaneStressCoupledThermoQuad(CoupledThermoPlaneStressStiffness *_stiff,
                               ElementBehaviorType type = LINEAR,
                               int _componentNum = 0);
  ~PlaneStressCoupledThermoQuad();

  // Return the name of this element
  // -------------------------------
  const char *elementName() { return elemName; }

  // Retrieve the shape functions
  // ----------------------------
  void getShapeFunctions(const double pt[], double N[]);
  void getShapeFunctions(const double pt[], double N[], double Na[],
                         double Nb[]);

  // Retrieve the Gauss points/weights
  // ---------------------------------
  int getNumGaussPts();
  double getGaussWtsPts(const int num, double pt[]);

  // Functions for post-processing
  // -----------------------------
  void addOutputCount(int *nelems, int *nnodes, int *ncsr);
  void getOutputData(unsigned int out_type, double *data, int ld_data,
                     const TacsScalar Xpts[], const TacsScalar vars[]);
  void getOutputConnectivity(int *con, int node);

 private:
  // The knot locations for the basis functions
  double knots[order];

  static const int NUM_NODES = order * order;

  // The Gauss quadrature scheme
  int numGauss;
  const double *gaussWts, *gaussPts;

  // Store the name of the element
  static const char *elemName;
};

template <int order>
PlaneStressCoupledThermoQuad<order>::PlaneStressCoupledThermoQuad(
    CoupledThermoPlaneStressStiffness *_stiff, ElementBehaviorType type,
    int _componentNum)
    : TACS2DCoupledThermoElement<order * order>(_stiff, type, _componentNum) {
  numGauss = FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);
  // Set the knot locations
  if (order == 2) {
    knots[0] = -1.0;
    knots[1] = 1.0;
  } else if (order == 3) {
    knots[0] = -1.0;
    knots[1] = 0.0;
    knots[2] = 1.0;
  } else {
    // Set a co-sine spacing for the knot locations
    for (int k = 0; k < order; k++) {
      knots[k] = -cos(M_PI * k / (order - 1));
    }
  }
}

template <int order>
PlaneStressCoupledThermoQuad<order>::~PlaneStressCoupledThermoQuad() {}

template <int order>
const char *PlaneStressCoupledThermoQuad<order>::elemName =
    "PlaneStressCoupledThermoQuad";
/*
  Get the number of Gauss points in the Gauss quadrature scheme
*/
template <int order>
int PlaneStressCoupledThermoQuad<order>::getNumGaussPts() {
  return numGauss * numGauss;
}

/*
  Get the Gauss points
*/
template <int order>
double PlaneStressCoupledThermoQuad<order>::getGaussWtsPts(int npoint,
                                                           double pt[]) {
  // Compute the n/m/p indices of the Gauss quadrature scheme
  int m = (int)((npoint) / (numGauss));
  int n = npoint - numGauss * m;

  pt[0] = gaussPts[n];
  pt[1] = gaussPts[m];

  return gaussWts[n] * gaussWts[m];
}

/*
  Evaluate the shape functions and their derivatives
*/
template <int order>
void PlaneStressCoupledThermoQuad<order>::getShapeFunctions(const double pt[],
                                                            double N[]) {
  double na[order], nb[order];
  FElibrary::lagrangeSFKnots(na, pt[0], knots, order);
  FElibrary::lagrangeSFKnots(nb, pt[1], knots, order);
  for (int j = 0; j < order; j++) {
    for (int i = 0; i < order; i++) {
      N[i + j * order] = na[i] * nb[j];
    }
  }
}

/*
  Compute the shape functions and their derivatives w.r.t. the
  parametric element location
*/
template <int order>
void PlaneStressCoupledThermoQuad<order>::getShapeFunctions(const double pt[],
                                                            double N[],
                                                            double Na[],
                                                            double Nb[]) {
  double na[order], nb[order];
  double dna[order], dnb[order];
  FElibrary::lagrangeSFKnots(na, dna, pt[0], knots, order);
  FElibrary::lagrangeSFKnots(nb, dnb, pt[1], knots, order);
  for (int j = 0; j < order; j++) {
    for (int i = 0; i < order; i++) {
      N[i + j * order] = na[i] * nb[j];
      Na[i + j * order] = dna[i] * nb[j];
      Nb[i + j * order] = na[i] * dnb[j];
    }
  }
}

/*
  Get the number of elemens/nodes and CSR size of the contributed by
  this element.
*/
template <int order>
void PlaneStressCoupledThermoQuad<order>::addOutputCount(int *nelems,
                                                         int *nnodes,
                                                         int *ncsr) {
  *nelems += (order - 1) * (order - 1);
  *nnodes += order * order;
  *ncsr += 4 * (order - 1) * (order - 1);
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
template <int order>
void PlaneStressCoupledThermoQuad<order>::getOutputData(
    unsigned int out_type, double *data, int ld_data, const TacsScalar Xpts[],
    const TacsScalar vars[]) {
  for (int m = 0; m < order; m++) {
    for (int n = 0; n < order; n++) {
      int p = n + order * m;
      int index = 0;
      // Set the parametric point to extract the data
      double pt[2];
      pt[0] = knots[n];
      pt[1] = knots[m];
      // Compute the shape functions
      double N[NUM_NODES];
      double Na[NUM_NODES], Nb[NUM_NODES];
      getShapeFunctions(pt, N, Na, Nb);

      if (out_type & TACSElement::OUTPUT_NODES) {
        for (int k = 0; k < 3; k++) {
          TacsScalar X = 0.0;
          for (int i = 0; i < NUM_NODES; i++) {
            X += N[i] * Xpts[3 * p + k];
          }
          data[index + k] = TacsRealPart(X);
        }
        index += 3;
      }
      if (out_type & TACSElement::OUTPUT_DISPLACEMENTS) {
        for (int k = 0; k < 3; k++) {
          TacsScalar u = 0.0;
          for (int i = 0; i < NUM_NODES; i++) {
            u += N[i] * vars[3 * p + k];
          }
          data[index + k] = TacsRealPart(u);
        }
        index += 3;
      }

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
        // Calculate the stress D*B*u at the current point
        TacsScalar stress[3];
        this->stiff->calculateStress(pt, strain, stress);
        for (int k = 0; k < 3; k++) {
          data[index + k] = TacsRealPart(stress[k]);
        }
        index += 3;
      }
      if (out_type & TACSElement::OUTPUT_EXTRAS) {
        // Get the temperature
        TacsScalar T[] = {0.0};
        ThermoQuad *elem = dynamic_cast<ThermoQuad *>(this);
        if (elem) {
          elem->getTemperature(T, N, vars);
        }
        // Compute the failure value
        TacsScalar lambda = 0.0;
        CoupledThermoPlaneStressStiffness *con =
            dynamic_cast<CoupledThermoPlaneStressStiffness *>(this->stiff);
        if (con) {
          con->failure(pt, T, strain, &lambda);
        }
        data[index] = TacsRealPart(lambda);

        this->stiff->buckling(strain, &lambda);
        data[index + 1] = TacsRealPart(lambda);

        data[index + 2] = TacsRealPart(this->stiff->getDVOutputValue(0, pt));
        data[index + 3] = TacsRealPart(this->stiff->getDVOutputValue(1, pt));

        index += this->NUM_EXTRAS;
      }

      data += ld_data;
    }
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
template <int order>
void PlaneStressCoupledThermoQuad<order>::getOutputConnectivity(int *con,
                                                                int node) {
  int p = 0;
  for (int m = 0; m < order - 1; m++) {
    for (int n = 0; n < order - 1; n++) {
      con[4 * p] = node + n + m * order;
      con[4 * p + 1] = node + n + 1 + m * order;
      con[4 * p + 2] = node + n + 1 + (m + 1) * order;
      con[4 * p + 3] = node + n + (m + 1) * order;
      p++;
    }
  }
}

#endif  // TACS_PLANE_STRESS_THERMO_QUAD_H
