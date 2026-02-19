#ifndef TACS_COUPLED_THERMO_SOLID_H
#define TACS_COUPLED_THERMO_SOLID_H

/*
  Solid element implementation

  Copyright (c) 2010-2015 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.

  The following code uses templates to allow for arbitrary order elements.
*/

#include "FElibrary.h"
#include "TACS3DCoupledThermoElement.h"
#include "TACSElement.h"

template <int order>
class CoupledThermoSolid
    : public TACS3DCoupledThermoElement<order * order * order> {
 public:
  CoupledThermoSolid(CoupledThermoSolidStiffness *_stiff,
                     ElementBehaviorType type = LINEAR, int _componentNum = 0);
  ~CoupledThermoSolid();

  // Return the name of this element
  // -------------------------------
  const char *elementName() { return elemName; }

  // Retrieve the shape functions
  // ----------------------------
  void getShapeFunctions(const double pt[], double N[]);
  void getShapeFunctions(const double pt[], double N[], double Na[],
                         double Nb[], double Nc[]);

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
  // The number of nodes in the element
  static const int NUM_NODES = order * order * order;

  // The knot locations for the basis functions
  double knots[order];

  // The Gauss quadrature scheme
  int numGauss;
  const double *gaussWts, *gaussPts;

  // Store the name of the element
  static const char *elemName;
};

template <int order>
CoupledThermoSolid<order>::CoupledThermoSolid(
    CoupledThermoSolidStiffness *_stiff, ElementBehaviorType type,
    int _componentNum)
    : TACS3DCoupledThermoElement<order * order * order>(_stiff, type,
                                                        _componentNum) {
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
CoupledThermoSolid<order>::~CoupledThermoSolid() {}

template <int order>
const char *CoupledThermoSolid<order>::elemName = "CoupledThermoSolid";

/*
  Get the number of Gauss points in the scheme
*/
template <int order>
int CoupledThermoSolid<order>::getNumGaussPts() {
  return numGauss * numGauss * numGauss;
}

/*
  Get the Gauss points
*/
template <int order>
double CoupledThermoSolid<order>::getGaussWtsPts(int npoint, double pt[]) {
  // Compute the n/m/p indices of the Gauss quadrature scheme
  int p = (int)((npoint) / (numGauss * numGauss));
  int m = (int)((npoint - numGauss * numGauss * p) / numGauss);
  int n = npoint - numGauss * m - numGauss * numGauss * p;

  pt[0] = gaussPts[n];
  pt[1] = gaussPts[m];
  pt[2] = gaussPts[p];

  return gaussWts[n] * gaussWts[m] * gaussWts[p];
}

/*
  Evaluate the shape functions and their derivatives
*/
template <int order>
void CoupledThermoSolid<order>::getShapeFunctions(const double pt[],
                                                  double N[]) {
  double na[order], nb[order], nc[order];
  FElibrary::lagrangeSFKnots(na, pt[0], knots, order);
  FElibrary::lagrangeSFKnots(nb, pt[1], knots, order);
  FElibrary::lagrangeSFKnots(nc, pt[2], knots, order);
  for (int k = 0; k < order; k++) {
    for (int j = 0; j < order; j++) {
      for (int i = 0; i < order; i++) {
        N[0] = na[i] * nb[j] * nc[k];
        N++;
      }
    }
  }
}

/*
  Compute the shape functions and their derivatives w.r.t. the
  parametric element location
*/
template <int order>
void CoupledThermoSolid<order>::getShapeFunctions(const double pt[], double N[],
                                                  double Na[], double Nb[],
                                                  double Nc[]) {
  double na[order], nb[order], nc[order];
  double dna[order], dnb[order], dnc[order];
  FElibrary::lagrangeSFKnots(na, dna, pt[0], knots, order);
  FElibrary::lagrangeSFKnots(nb, dnb, pt[1], knots, order);
  FElibrary::lagrangeSFKnots(nc, dnc, pt[2], knots, order);
  for (int k = 0; k < order; k++) {
    for (int j = 0; j < order; j++) {
      for (int i = 0; i < order; i++) {
        N[0] = na[i] * nb[j] * nc[k];
        Na[0] = dna[i] * nb[j] * nc[k];
        Nb[0] = na[i] * dnb[j] * nc[k];
        Nc[0] = na[i] * nb[j] * dnc[k];
        N++;
        Na++;
        Nb++;
        Nc++;
      }
    }
  }
}

/*
  Get the number of elemens/nodes and CSR size of the contributed by
  this element.
*/
template <int order>
void CoupledThermoSolid<order>::addOutputCount(int *nelems, int *nnodes,
                                               int *ncsr) {
  *nelems += (order - 1) * (order - 1) * (order - 1);
  *nnodes += order * order * order;
  *ncsr += 8 * (order - 1) * (order - 1) * (order - 1);
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
void CoupledThermoSolid<order>::getOutputData(unsigned int out_type,
                                              double *data, int ld_data,
                                              const TacsScalar Xpts[],
                                              const TacsScalar vars[]) {
  for (int p = 0; p < order; p++) {
    for (int m = 0; m < order; m++) {
      for (int n = 0; n < order; n++) {
        int node = n + m * order + p * order * order;
        int index = 0;
        // Set the parametric point where to evaluate the
        // stresses/strains
        double pt[3];
        pt[0] = knots[n];
        pt[1] = knots[m];
        pt[2] = knots[p];

        // Compute the shape functions
        double N[NUM_NODES];
        double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
        getShapeFunctions(pt, N, Na, Nb, Nc);

        if (out_type & TACSElement::OUTPUT_NODES) {
          for (int k = 0; k < 3; k++) {
            TacsScalar X = 0.0;
            for (int i = 0; i < NUM_NODES; i++) {
              X += N[i] * Xpts[3 * node + k];
            }
            data[index + k] = TacsRealPart(X);
          }
          index += 3;
        }
        if (out_type & TACSElement::OUTPUT_DISPLACEMENTS) {
          for (int k = 0; k < 4; k++) {
            TacsScalar u = 0.0;
            for (int i = 0; i < NUM_NODES; i++) {
              u += N[i] * vars[4 * node + k];
            }
            data[index + k] = TacsRealPart(u);
          }
          index += 4;
        }

        // Compute the derivative of X with respect to the
        // coordinate directions
        TacsScalar X[3], Xa[9];
        this->solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

        // Compute the determinant of Xa and the transformation
        TacsScalar J[9];
        FElibrary::jacobian3d(Xa, J);

        // Compute the strain B*u
        TacsScalar strain[6];
        this->evalStrain(strain, J, Na, Nb, Nc, vars);

        if (out_type & TACSElement::OUTPUT_STRAINS) {
          for (int k = 0; k < 6; k++) {
            data[index + k] = TacsRealPart(strain[k]);
          }
          index += 6;
        }
        if (out_type & TACSElement::OUTPUT_STRESSES) {
          // Calculate the stress D*B*u at the current point
          TacsScalar stress[6];
          this->stiff->calculateStress(pt, strain, stress);

          for (int k = 0; k < 6; k++) {
            data[index + k] = TacsRealPart(stress[k]);
          }
          index += 6;
        }
        if (out_type & TACSElement::OUTPUT_EXTRAS) {
          // Get the temperature
          TacsScalar T[] = {0.0};
          ThermoSolid *elem = dynamic_cast<ThermoSolid *>(this);
          if (elem) {
            elem->getTemperature(T, N, vars);
          }
          // Compute the failure value
          TacsScalar lambda = 0.0;
          CoupledThermoSolidStiffness *con =
              dynamic_cast<CoupledThermoSolidStiffness *>(this->stiff);
          if (con) {
            con->failure(pt, T, strain, &lambda);
          }
          data[index] = TacsRealPart(lambda);
          lambda = 0.0;
          if (con) {
            con->maxtemp(pt, T[0], &lambda);
            // this->stiff->buckling(strain, &lambda);
          }
          data[index + 1] = TacsRealPart(lambda);

          data[index + 2] = TacsRealPart(this->stiff->getDVOutputValue(0, pt));
          data[index + 3] = TacsRealPart(this->stiff->getDVOutputValue(1, pt));

          index += this->NUM_EXTRAS;
        }

        data += ld_data;
      }
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
void CoupledThermoSolid<order>::getOutputConnectivity(int *con, int node) {
  int j = 0;
  for (int p = 0; p < order - 1; p++) {
    for (int m = 0; m < order - 1; m++) {
      for (int n = 0; n < order - 1; n++) {
        con[8 * j] = node + n + m * order + p * order * order;
        con[8 * j + 1] = node + n + 1 + m * order + p * order * order;
        con[8 * j + 2] = node + n + 1 + (m + 1) * order + p * order * order;
        con[8 * j + 3] = node + n + (m + 1) * order + p * order * order;
        con[8 * j + 4] = node + n + m * order + (p + 1) * order * order;
        con[8 * j + 5] = node + n + 1 + m * order + (p + 1) * order * order;
        con[8 * j + 6] =
            node + n + 1 + (m + 1) * order + (p + 1) * order * order;
        con[8 * j + 7] = node + n + (m + 1) * order + (p + 1) * order * order;
        j++;
      }
    }
  }
}

#endif  // TACS_SOLID_H
