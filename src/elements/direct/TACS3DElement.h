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

#ifndef TACS_3D_ELEMENT_H
#define TACS_3D_ELEMENT_H

/*
  The following file contains the general definition of a
  three-dimensional element that can be used in TACS.
*/

#include "FElibrary.h"
#include "SolidStiffness.h"
#include "TACSElement.h"

/*
  The following class defines a generic three-dimensional element
  without defining the strain expressions, shape functions or
  quadrature scheme. This class can be used to implement multiple 3D
  elements that could be either linear/nonlinear later on.  This does
  not significantly impact the computational performance since the
  cost of the element computations is consumed in the inner product of
  the B-matrix with the constitutive matrix.
*/
template <int NUM_NODES>
class TACS3DElement : public TACSElement {
 public:
  // Define some constants for this element type
  static const int NUM_DISPS = 3;
  static const int NUM_STRESSES = 6;
  static const int NUM_EXTRAS = 4;
  static const int NUM_VARIABLES = 3 * NUM_NODES;

  TACS3DElement(SolidStiffness *_stiff, ElementBehaviorType type,
                int _component);
  ~TACS3DElement();

  // Retrieve the shape functions
  // ----------------------------
  virtual void getShapeFunctions(const double pt[], double N[], double Na[],
                                 double Nb[], double Nc[]) = 0;

  // Compute the position vector and its derivative
  // ----------------------------------------------
  void solidJacobian(TacsScalar X[], TacsScalar Xa[], const double N[],
                     const double Na[], const double Nb[], const double Nc[],
                     const TacsScalar Xpts[]);

  // Compute the displacement gradient
  // ---------------------------------
  void getDisplacement(TacsScalar U[], const double N[],
                       const TacsScalar vars[]);
  void getDisplGradient(TacsScalar Ud[], const TacsScalar J[],
                        const double Na[], const double Nb[], const double Nc[],
                        const TacsScalar vars[]);
  void getDisplGradientSens(TacsScalar Ud[], TacsScalar UdSens[],
                            const TacsScalar J[], const TacsScalar JSens[],
                            const double Na[], const double Nb[],
                            const double Nc[], const TacsScalar vars[]);

  // Compute the strain in the element
  // ---------------------------------
  void evalStrain(TacsScalar strain[], const TacsScalar J[], const double Na[],
                  const double Nb[], const double Nc[],
                  const TacsScalar vars[]);

  // Compute the derivative of the strain with respect to vars
  // ---------------------------------------------------------
  void getBmat(TacsScalar B[], const TacsScalar J[], const double Na[],
               const double Nb[], const double Nc[], const TacsScalar vars[]);

  // Add the second derivatives of the strain times the stress to the
  // upper portion of the matrix
  // ----------------------------------------------------------------
  void addGeoStiffness(TacsScalar kmat[], TacsScalar h,
                       const TacsScalar stress[], const TacsScalar J[],
                       const double Na[], const double Nb[], const double Nc[]);

  // Compute the derivative of the strain with respect to the nodal coordinates
  // --------------------------------------------------------------------------
  void addStrainXptSens(TacsScalar sens[], TacsScalar scale,
                        const TacsScalar strainSens[], const TacsScalar J[],
                        const TacsScalar Xa[], const double Na[],
                        const double Nb[], const double Nc[],
                        const TacsScalar vars[]);

  // The design variable query functions
  // -----------------------------------
  void setDesignVars(const TacsScalar dvs[], int numDVs);
  void getDesignVars(TacsScalar dvs[], int numDVs);
  void getDesignVarRange(TacsScalar lowerBound[], TacsScalar upperBound[],
                         int numDVs);

  // Get the variable information
  // ----------------------------
  const char *displacementName(int i);
  const char *stressName(int i);
  const char *strainName(int i);
  const char *extraName(int i);
  int numDisplacements();
  int numStresses();
  int numNodes();
  int numVariables();
  int numExtras();
  ElementType getElementType();

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies(double time, TacsScalar *_Te, TacsScalar *_Pe,
                       const TacsScalar Xpts[], const TacsScalar vars[],
                       const TacsScalar dvars[]);

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual(double time, TacsScalar res[], const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[]);

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian(double time, TacsScalar J[], double alpha, double beta,
                   double gamma, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[]);

  // Add the product of the adjoint with the derivative of the design variables
  // --------------------------------------------------------------------------
  void addAdjResProduct(double time, double scale, TacsScalar dvSens[],
                        int dvLen, const TacsScalar psi[],
                        const TacsScalar Xpts[], const TacsScalar vars[],
                        const TacsScalar dvars[], const TacsScalar ddvars[]);

  // Add the product of the adjoint with the derivative of the design variables
  // --------------------------------------------------------------------------
  void addAdjResXptProduct(double time, double scale, TacsScalar fXptSens[],
                           const TacsScalar psi[], const TacsScalar Xpts[],
                           const TacsScalar vars[], const TacsScalar dvars[],
                           const TacsScalar ddvars[]);

  // Retrieve a specific time-independent matrix from the element
  // ------------------------------------------------------------
  void getMatType(ElementMatrixType matType, TacsScalar mat[],
                  const TacsScalar Xpts[], const TacsScalar vars[]);

  // Compute the derivative of the inner product w.r.t. design variables
  // -------------------------------------------------------------------
  void addMatDVSensInnerProduct(ElementMatrixType matType, double scale,
                                TacsScalar dvSens[], int dvLen,
                                const TacsScalar psi[], const TacsScalar phi[],
                                const TacsScalar Xpts[],
                                const TacsScalar vars[]);

  // Compute the derivative of the inner product w.r.t. vars[]
  // ---------------------------------------------------------
  void getMatSVSensInnerProduct(ElementMatrixType matType, TacsScalar res[],
                                const TacsScalar psi[], const TacsScalar phi[],
                                const TacsScalar Xpts[],
                                const TacsScalar vars[]);

  // Functions for evaluating global functionals of interest
  // -------------------------------------------------------
  TACSConstitutive *getConstitutive() { return stiff; }

  // Evaluate the determinant of the Jacobian and its derivative
  // -----------------------------------------------------------
  TacsScalar getDetJacobian(const double pt[], const TacsScalar Xpts[]);
  TacsScalar getDetJacobianXptSens(TacsScalar hXptSens[], const double pt[],
                                   const TacsScalar Xpts[]);

  // Compute the point-wise strain and its derivative
  // ------------------------------------------------
  void getStrain(TacsScalar strain[], const double pt[],
                 const TacsScalar Xpts[], const TacsScalar vars[]);
  void addStrainXptSens(TacsScalar strainXptSens[], const double pt[],
                        const TacsScalar scale, const TacsScalar strainSens[],
                        const TacsScalar Xpts[], const TacsScalar vars[]);
  void addStrainSVSens(TacsScalar strainSVSens[], const double pt[],
                       const TacsScalar scale, const TacsScalar strainSens[],
                       const TacsScalar Xpts[], const TacsScalar vars[]);

 protected:
  ElementBehaviorType strain_type;
  SolidStiffness *stiff;

 private:
  static const char *dispNames[NUM_DISPS];
  static const char *stressNames[NUM_STRESSES];
  static const char *strainNames[NUM_STRESSES];
  static const char *extraNames[NUM_EXTRAS];
};

/*
  The constructor for the 3D element matrix
*/
template <int NUM_NODES>
TACS3DElement<NUM_NODES>::TACS3DElement(SolidStiffness *_stiff,
                                        ElementBehaviorType type, int component)
    : TACSElement(component) {
  strain_type = type;
  stiff = _stiff;
  stiff->incref();
}

template <int NUM_NODES>
TACS3DElement<NUM_NODES>::~TACS3DElement() {
  stiff->decref();
}

/*
  Provide the names for the different components of the displacements
  and stresses
*/
template <int NUM_NODES>
const char *TACS3DElement<NUM_NODES>::dispNames[] = {"u", "v", "w"};

template <int NUM_NODES>
const char *TACS3DElement<NUM_NODES>::stressNames[] = {"sxx", "syy", "szz",
                                                       "syz", "sxz", "sxy"};

template <int NUM_NODES>
const char *TACS3DElement<NUM_NODES>::strainNames[] = {"exx", "eyy", "ezz",
                                                       "eyz", "exz", "exy"};

template <int NUM_NODES>
const char *TACS3DElement<NUM_NODES>::extraNames[] = {"lambda", "buckling",
                                                      "dv1", "dv2"};

/*
  Get the names of the displacements/stress etc.
*/
template <int NUM_NODES>
const char *TACS3DElement<NUM_NODES>::displacementName(int i) {
  if (i >= 0 && i < NUM_DISPS) {
    return dispNames[i];
  }
  return NULL;
}

template <int NUM_NODES>
const char *TACS3DElement<NUM_NODES>::stressName(int i) {
  if (i >= 0 && i < NUM_STRESSES) {
    return stressNames[i];
  }
  return NULL;
}

template <int NUM_NODES>
const char *TACS3DElement<NUM_NODES>::strainName(int i) {
  if (i >= 0 && i < NUM_STRESSES) {
    return strainNames[i];
  }
  return NULL;
}

template <int NUM_NODES>
const char *TACS3DElement<NUM_NODES>::extraName(int i) {
  if (i >= 0 && i < NUM_EXTRAS) {
    return extraNames[i];
  }
  return NULL;
}

/*
  Retrieve information about the number of displacements/stress etc.
*/
template <int NUM_NODES>
int TACS3DElement<NUM_NODES>::numDisplacements() {
  return NUM_DISPS;
}

template <int NUM_NODES>
int TACS3DElement<NUM_NODES>::numStresses() {
  return NUM_STRESSES;
}

template <int NUM_NODES>
int TACS3DElement<NUM_NODES>::numNodes() {
  return NUM_NODES;
}

template <int NUM_NODES>
int TACS3DElement<NUM_NODES>::numVariables() {
  return NUM_VARIABLES;
}

template <int NUM_NODES>
int TACS3DElement<NUM_NODES>::numExtras() {
  return NUM_EXTRAS;
}

template <int NUM_NODES>
ElementType TACS3DElement<NUM_NODES>::getElementType() {
  return TACS_SOLID;
}

/*
  Set the design variable values
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::setDesignVars(const TacsScalar dvs[],
                                             int numDVs) {
  stiff->setDesignVars(dvs, numDVs);
}

/*
  Retrive the design variable numbers
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getDesignVars(TacsScalar dvs[], int numDVs) {
  stiff->getDesignVars(dvs, numDVs);
}

/*
  Set the design variable lower/upper bounds in the provided arrays
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getDesignVarRange(TacsScalar lowerBound[],
                                                 TacsScalar upperBound[],
                                                 int numDVs) {
  stiff->getDesignVarRange(lowerBound, upperBound, numDVs);
}

/*
  Compute the product of the shape functions with the element shape
  functions.

  output:
  X:   the x,y,z location within the structure
  Xa:  the derivative of x,y,z with respect to the parametric location

  input:
  N:           the shape functions
  Na, Nb, Nc:  the derivative of the shape functions
  Xpts:        the nodal locations
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::solidJacobian(
    TacsScalar X[], TacsScalar Xa[], const double N[], const double Na[],
    const double Nb[], const double Nc[], const TacsScalar Xpts[]) {
  X[0] = X[1] = X[2] = 0.0;
  Xa[0] = Xa[1] = Xa[2] = 0.0;
  Xa[3] = Xa[4] = Xa[5] = 0.0;
  Xa[6] = Xa[7] = Xa[8] = 0.0;

  for (int i = 0; i < NUM_NODES; i++) {
    X[0] += Xpts[0] * N[0];
    X[1] += Xpts[1] * N[0];
    X[2] += Xpts[2] * N[0];

    Xa[0] += Xpts[0] * Na[0];
    Xa[1] += Xpts[0] * Nb[0];
    Xa[2] += Xpts[0] * Nc[0];

    Xa[3] += Xpts[1] * Na[0];
    Xa[4] += Xpts[1] * Nb[0];
    Xa[5] += Xpts[1] * Nc[0];

    Xa[6] += Xpts[2] * Na[0];
    Xa[7] += Xpts[2] * Nb[0];
    Xa[8] += Xpts[2] * Nc[0];

    N++;
    Na++;
    Nb++;
    Nc++;
    Xpts += 3;
  }
}

/*
  Compute the displacement given the provided values of the shape functions
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getDisplacement(TacsScalar U[], const double N[],
                                               const TacsScalar vars[]) {
  // Zero the displacement values
  U[0] = U[1] = U[2] = 0.0;

  for (int i = 0; i < NUM_NODES; i++) {
    U[0] += N[0] * vars[0];
    U[1] += N[0] * vars[1];
    U[2] += N[0] * vars[2];

    N++;
    vars += 3;
  }
}

/*
  Compute the displacement gradient using the provided basis functions.

  The displacement gradient is computed as follows:
  Ud = Ua*J =
  [ Ua[0]  Ua[1]  Ua[2] ][ J[0]  J[1]  J[2] ]
  [ Ua[3]  Ua[4]  Ua[5] ][ J[3]  J[4]  J[5] ]
  [ Ua[6]  Ua[7]  Ua[8] ][ J[6]  J[7]  J[8] ]

  output:
  Ud:    the displacement gradient

  input:
  J:          the transformation J = [X,a]^{-1}
  Na, Nb, Nc: the derivatives of the shape functions
  vars:       the element variables
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getDisplGradient(
    TacsScalar Ud[], const TacsScalar J[], const double Na[], const double Nb[],
    const double Nc[], const TacsScalar vars[]) {
  TacsScalar Ua[9];
  Ua[0] = Ua[1] = Ua[2] = Ua[3] = Ua[4] = Ua[5] = Ua[6] = Ua[7] = Ua[8] = 0.0;

  // Compute the derivative of the u,v,w displacements with
  // respect to the parametric locations within the element
  for (int i = 0; i < NUM_NODES; i++) {
    Ua[0] += vars[0] * Na[0];
    Ua[1] += vars[0] * Nb[0];
    Ua[2] += vars[0] * Nc[0];

    Ua[3] += vars[1] * Na[0];
    Ua[4] += vars[1] * Nb[0];
    Ua[5] += vars[1] * Nc[0];

    Ua[6] += vars[2] * Na[0];
    Ua[7] += vars[2] * Nb[0];
    Ua[8] += vars[2] * Nc[0];

    Na++;
    Nb++;
    Nc++;
    vars += 3;
  }

  // Compute the displacement gradient: Ud = Ua*J
  Ud[0] = Ua[0] * J[0] + Ua[1] * J[3] + Ua[2] * J[6];
  Ud[3] = Ua[3] * J[0] + Ua[4] * J[3] + Ua[5] * J[6];
  Ud[6] = Ua[6] * J[0] + Ua[7] * J[3] + Ua[8] * J[6];

  Ud[1] = Ua[0] * J[1] + Ua[1] * J[4] + Ua[2] * J[7];
  Ud[4] = Ua[3] * J[1] + Ua[4] * J[4] + Ua[5] * J[7];
  Ud[7] = Ua[6] * J[1] + Ua[7] * J[4] + Ua[8] * J[7];

  Ud[2] = Ua[0] * J[2] + Ua[1] * J[5] + Ua[2] * J[8];
  Ud[5] = Ua[3] * J[2] + Ua[4] * J[5] + Ua[5] * J[8];
  Ud[8] = Ua[6] * J[2] + Ua[7] * J[5] + Ua[8] * J[8];
}

/*
  Compute the displacement gradient and the derivative of the
  displacement gradient using the provided basis functions.

  The displacement gradient is computed as follows:
  Ud = Ua*J =
  [ Ua[0]  Ua[1]  Ua[2] ][ J[0]  J[1]  J[2] ]
  [ Ua[3]  Ua[4]  Ua[5] ][ J[3]  J[4]  J[5] ]
  [ Ua[6]  Ua[7]  Ua[8] ][ J[6]  J[7]  J[8] ]

  output:
  Ud:       the displacement gradient
  UdSens:   the derivative of the displacement gradient

  input:
  J:          the transformation J = [X,a]^{-1}
  JSens:      the derivative of the inverse Jacobian
  Na, Nb, Nc: the derivatives of the shape functions
  vars:       the element variables
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getDisplGradientSens(
    TacsScalar Ud[], TacsScalar UdSens[], const TacsScalar J[],
    const TacsScalar JSens[], const double Na[], const double Nb[],
    const double Nc[], const TacsScalar vars[]) {
  TacsScalar Ua[9];
  Ua[0] = Ua[1] = Ua[2] = Ua[3] = Ua[4] = Ua[5] = Ua[6] = Ua[7] = Ua[8] = 0.0;

  // Compute the derivative of the u,v,w displacements with
  // respect to the parametric locations within the element
  for (int i = 0; i < NUM_NODES; i++) {
    Ua[0] += vars[0] * Na[0];
    Ua[1] += vars[0] * Nb[0];
    Ua[2] += vars[0] * Nc[0];

    Ua[3] += vars[1] * Na[0];
    Ua[4] += vars[1] * Nb[0];
    Ua[5] += vars[1] * Nc[0];

    Ua[6] += vars[2] * Na[0];
    Ua[7] += vars[2] * Nb[0];
    Ua[8] += vars[2] * Nc[0];

    Na++;
    Nb++;
    Nc++;
    vars += 3;
  }

  // Compute the displacement gradient: Ud = Ua*J
  Ud[0] = Ua[0] * J[0] + Ua[1] * J[3] + Ua[2] * J[6];
  Ud[3] = Ua[3] * J[0] + Ua[4] * J[3] + Ua[5] * J[6];
  Ud[6] = Ua[6] * J[0] + Ua[7] * J[3] + Ua[8] * J[6];

  Ud[1] = Ua[0] * J[1] + Ua[1] * J[4] + Ua[2] * J[7];
  Ud[4] = Ua[3] * J[1] + Ua[4] * J[4] + Ua[5] * J[7];
  Ud[7] = Ua[6] * J[1] + Ua[7] * J[4] + Ua[8] * J[7];

  Ud[2] = Ua[0] * J[2] + Ua[1] * J[5] + Ua[2] * J[8];
  Ud[5] = Ua[3] * J[2] + Ua[4] * J[5] + Ua[5] * J[8];
  Ud[8] = Ua[6] * J[2] + Ua[7] * J[5] + Ua[8] * J[8];

  // Compute the derivative of the displacement gradient
  UdSens[0] = Ua[0] * JSens[0] + Ua[1] * JSens[3] + Ua[2] * JSens[6];
  UdSens[3] = Ua[3] * JSens[0] + Ua[4] * JSens[3] + Ua[5] * JSens[6];
  UdSens[6] = Ua[6] * JSens[0] + Ua[7] * JSens[3] + Ua[8] * JSens[6];

  UdSens[1] = Ua[0] * JSens[1] + Ua[1] * JSens[4] + Ua[2] * JSens[7];
  UdSens[4] = Ua[3] * JSens[1] + Ua[4] * JSens[4] + Ua[5] * JSens[7];
  UdSens[7] = Ua[6] * JSens[1] + Ua[7] * JSens[4] + Ua[8] * JSens[7];

  UdSens[2] = Ua[0] * JSens[2] + Ua[1] * JSens[5] + Ua[2] * JSens[8];
  UdSens[5] = Ua[3] * JSens[2] + Ua[4] * JSens[5] + Ua[5] * JSens[8];
  UdSens[8] = Ua[6] * JSens[2] + Ua[7] * JSens[5] + Ua[8] * JSens[8];
}

/*
  Compute the strain using the specified transformation, shape
  functions and variables

  Note: The strain is computed as follows:

  strain = 0.5*(Ud + Ud^{T} + Ud^{T}*Ud)
  =
  [ Ud[0]  0.5*(Ud[1] + Ud[3])  0.5*(Ud[2] + Ud[6]) ]
  [             Ud[4]           0.5*(Ud[5] + Ud[7]) ] +
  [                                  Ud[8]          ]

  [ Ud[0]  Ud[3]  Ud[6] ][ Ud[0]  Ud[1]  Ud[2] ]
  [ Ud[1]  Ud[4]  Ud[7] ][ Ud[3]  Ud[4]  Ud[5] ]
  [ Ud[2]  Ud[5]  Ud[8] ][ Ud[6]  Ud[7]  Ud[8] ]

  output:
  strain:   the strain

  input:
  J:          the Jacobian of the transformation
  Na, Nb, Nc: the derivatives of the basis functions
  vars:       the variables
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::evalStrain(TacsScalar strain[],
                                          const TacsScalar J[],
                                          const double Na[], const double Nb[],
                                          const double Nc[],
                                          const TacsScalar vars[]) {
  // Compute the displacement gradient
  TacsScalar Ud[9];
  getDisplGradient(Ud, J, Na, Nb, Nc, vars);

  // Compute the strain using either linear or nonlinear expression
  if (strain_type == LINEAR) {
    strain[0] = Ud[0];
    strain[1] = Ud[4];
    strain[2] = Ud[8];

    strain[3] = Ud[5] + Ud[7];
    strain[4] = Ud[2] + Ud[6];
    strain[5] = Ud[1] + Ud[3];
  } else {
    strain[0] = Ud[0] + 0.5 * (Ud[0] * Ud[0] + Ud[3] * Ud[3] + Ud[6] * Ud[6]);
    strain[1] = Ud[4] + 0.5 * (Ud[1] * Ud[1] + Ud[4] * Ud[4] + Ud[7] * Ud[7]);
    strain[2] = Ud[8] + 0.5 * (Ud[2] * Ud[2] + Ud[5] * Ud[5] + Ud[8] * Ud[8]);

    strain[3] = Ud[5] + Ud[7] + (Ud[1] * Ud[2] + Ud[4] * Ud[5] + Ud[7] * Ud[8]);
    strain[4] = Ud[2] + Ud[6] + (Ud[0] * Ud[2] + Ud[3] * Ud[5] + Ud[6] * Ud[8]);
    strain[5] = Ud[1] + Ud[3] + (Ud[0] * Ud[1] + Ud[3] * Ud[4] + Ud[6] * Ud[7]);
  }
}

/*
  Compute the derivative of the strain with respect to vars
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getBmat(TacsScalar B[], const TacsScalar J[],
                                       const double Na[], const double Nb[],
                                       const double Nc[],
                                       const TacsScalar vars[]) {
  if (strain_type == LINEAR) {
    // If this is a linear element, then things are relative easy to
    // deal with - we just compute B alternatively by row
    for (int i = 0; i < NUM_NODES; i++) {
      TacsScalar Dx = Na[0] * J[0] + Nb[0] * J[3] + Nc[0] * J[6];
      TacsScalar Dy = Na[0] * J[1] + Nb[0] * J[4] + Nc[0] * J[7];
      TacsScalar Dz = Na[0] * J[2] + Nb[0] * J[5] + Nc[0] * J[8];

      B[0] = Dx;
      B[1] = 0.0;
      B[2] = 0.0;
      B[3] = 0.0;
      B[4] = Dz;
      B[5] = Dy;
      B += 6;

      B[0] = 0.0;
      B[1] = Dy;
      B[2] = 0.0;
      B[3] = Dz;
      B[4] = 0.0;
      B[5] = Dx;
      B += 6;

      B[0] = 0.0;
      B[1] = 0.0;
      B[2] = Dz;
      B[3] = Dy;
      B[4] = Dx;
      B[5] = 0.0;
      B += 6;

      Na++;
      Nb++;
      Nc++;
    }
  } else {
    // Compute the displacement gradient: Ud = Ua*J
    TacsScalar Ud[9];
    getDisplGradient(Ud, J, Na, Nb, Nc, vars);

    // Compute the derivative of the strain with respect to
    // the nodal displacements
    for (int i = 0; i < NUM_NODES; i++) {
      TacsScalar Dx = Na[0] * J[0] + Nb[0] * J[3] + Nc[0] * J[6];
      TacsScalar Dy = Na[0] * J[1] + Nb[0] * J[4] + Nc[0] * J[7];
      TacsScalar Dz = Na[0] * J[2] + Nb[0] * J[5] + Nc[0] * J[8];

      B[0] = Dx + Ud[0] * Dx;
      B[1] = Ud[1] * Dy;
      B[2] = Ud[2] * Dz;
      B[3] = Dy * Ud[2] + Ud[1] * Dz;
      B[4] = Dz + (Dx * Ud[2] + Ud[0] * Dz);
      B[5] = Dy + (Dx * Ud[1] + Ud[0] * Dy);
      B += 6;

      B[0] = Ud[3] * Dx;
      B[1] = Dy + Ud[4] * Dy;
      B[2] = Ud[5] * Dz;
      B[3] = Dz + (Dy * Ud[5] + Ud[4] * Dz);
      B[4] = Dx * Ud[5] + Ud[3] * Dz;
      B[5] = Dx + (Dx * Ud[4] + Ud[3] * Dy);
      B += 6;

      B[0] = Ud[6] * Dx;
      B[1] = Ud[7] * Dy;
      B[2] = Dz + Ud[8] * Dz;
      B[3] = Dy + (Dy * Ud[8] + Ud[7] * Dz);
      B[4] = Dx + (Dx * Ud[8] + Ud[6] * Dz);
      B[5] = Dx * Ud[7] + Ud[6] * Dy;
      B += 6;

      Na++;
      Nb++;
      Nc++;
    }
  }
}

/*
  Add the second derivatives of the strain times the stress to the
  upper portion of the matrix
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::addGeoStiffness(TacsScalar mat[], TacsScalar h,
                                               const TacsScalar stress[],
                                               const TacsScalar J[],
                                               const double Na[],
                                               const double Nb[],
                                               const double Nc[]) {
  if (!(strain_type == LINEAR)) {
    for (int j = 0; j < NUM_NODES; j++) {
      TacsScalar Dxj = Na[j] * J[0] + Nb[j] * J[3] + Nc[j] * J[6];
      TacsScalar Dyj = Na[j] * J[1] + Nb[j] * J[4] + Nc[j] * J[7];
      TacsScalar Dzj = Na[j] * J[2] + Nb[j] * J[5] + Nc[j] * J[8];

      for (int i = 0; i < NUM_NODES; i++) {
        TacsScalar Dxi = Na[i] * J[0] + Nb[i] * J[3] + Nc[i] * J[6];
        TacsScalar Dyi = Na[i] * J[1] + Nb[i] * J[4] + Nc[i] * J[7];
        TacsScalar Dzi = Na[i] * J[2] + Nb[i] * J[5] + Nc[i] * J[8];

        // Add the contributions to the stiffness matrix
        TacsScalar scale =
            h * (stress[0] * Dxi * Dxj + stress[1] * Dyi * Dyj +
                 stress[2] * Dzi * Dzj + stress[3] * (Dyi * Dzj + Dyj * Dzi) +
                 stress[4] * (Dxi * Dzj + Dxj * Dzi) +
                 stress[5] * (Dxi * Dyj + Dxj * Dyi));

        mat[3 * i + 3 * j * NUM_VARIABLES] += scale;
        mat[3 * i + 1 + (3 * j + 1) * NUM_VARIABLES] += scale;
        mat[3 * i + 2 + (3 * j + 2) * NUM_VARIABLES] += scale;
      }
    }
  }
}

/*
  Compute the derivative of the strain with respect to the nodal
  coordinates

  Note that the derivative of the Jacobian matrix can be computed as
  follows:

  dJ/dx = - J*d(Xa)/dx*J

  where d(Xa)/dx is the derivative of the parametric derivative with
  respect to the nodal position. This derivative takes the following
  form:

  dJ/dXa =
  [ J[0]  J[1]  J[2] ][ Na[0]  Nb[0]  Nc[0] ][ J[0]  J[1]  J[2] ]
  [ J[3]  J[4]  J[5] ][     0      0      0 ][ J[3]  J[4]  J[5] ]
  [ J[6]  J[7]  J[8] ][     0      0      0 ][ J[6]  J[7]  J[8] ]

  (where the middle term is the derivative d(Xa)/dx). These
  derivatives take a special form as follows:

  For the x-derivative:
  dJ/dXa =
  [ J[0]  J[1]  J[2] ][ d1  d2  d3 ]
  [ J[3]  J[4]  J[5] ][  0   0   0 ]
  [ J[6]  J[7]  J[8] ][  0   0   0 ]

  =
  [ d1*J[0]  d2*J[0]  d3*J[0] ]
  [ d1*J[3]  d2*J[3]  d3*J[3] ]
  [ d1*J[6]  d2*J[6]  d3*J[6] ]

  where:
  d1 = Na[0]*J[0] + Nb[0]*J[3] + Nc[0]*J[6]
  d2 = Na[0]*J[1] + Nb[0]*J[4] + Nc[0]*J[7]
  d3 = Na[0]*J[2] + Nb[0]*J[5] + Nc[0]*J[8]
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::addStrainXptSens(
    TacsScalar sens[], TacsScalar scale, const TacsScalar strainSens[],
    const TacsScalar J[], const TacsScalar Xa[], const double Na[],
    const double Nb[], const double Nc[], const TacsScalar vars[]) {
  // Compute the derivative of the u,v,w displacements with
  // respect to the parametric locations within the element
  TacsScalar Ua[9];
  Ua[0] = Ua[1] = Ua[2] = Ua[3] = Ua[4] = Ua[5] = Ua[6] = Ua[7] = Ua[8] = 0.0;

  const double *na = Na, *nb = Nb, *nc = Nc;
  for (int i = 0; i < NUM_NODES; i++) {
    Ua[0] += vars[0] * na[0];
    Ua[1] += vars[0] * nb[0];
    Ua[2] += vars[0] * nc[0];

    Ua[3] += vars[1] * na[0];
    Ua[4] += vars[1] * nb[0];
    Ua[5] += vars[1] * nc[0];

    Ua[6] += vars[2] * na[0];
    Ua[7] += vars[2] * nb[0];
    Ua[8] += vars[2] * nc[0];

    na++;
    nb++;
    nc++;
    vars += 3;
  }

  // Compute the scaled strain sensitivity
  TacsScalar eSens[6];
  eSens[0] = scale * strainSens[0];
  eSens[1] = scale * strainSens[1];
  eSens[2] = scale * strainSens[2];
  eSens[3] = scale * strainSens[3];
  eSens[4] = scale * strainSens[4];
  eSens[5] = scale * strainSens[5];

  if (strain_type == LINEAR) {
    for (int i = 0; i < NUM_NODES; i++) {
      // JSens = -J*d(Xa)/dx*J
      TacsScalar d1 = -(Na[0] * J[0] + Nb[0] * J[3] + Nc[0] * J[6]);
      TacsScalar d2 = -(Na[0] * J[1] + Nb[0] * J[4] + Nc[0] * J[7]);
      TacsScalar d3 = -(Na[0] * J[2] + Nb[0] * J[5] + Nc[0] * J[8]);

      // The derivative of the inverse transformation matrix
      // in each coordinate direction
      TacsScalar Ja[9], Jb[9], Jc[9];
      Ja[0] = d1 * J[0];
      Ja[1] = d2 * J[0];
      Ja[2] = d3 * J[0];
      Ja[3] = d1 * J[3];
      Ja[4] = d2 * J[3];
      Ja[5] = d3 * J[3];
      Ja[6] = d1 * J[6];
      Ja[7] = d2 * J[6];
      Ja[8] = d3 * J[6];

      Jb[0] = d1 * J[1];
      Jb[1] = d2 * J[1];
      Jb[2] = d3 * J[1];
      Jb[3] = d1 * J[4];
      Jb[4] = d2 * J[4];
      Jb[5] = d3 * J[4];
      Jb[6] = d1 * J[7];
      Jb[7] = d2 * J[7];
      Jb[8] = d3 * J[7];

      Jc[0] = d1 * J[2];
      Jc[1] = d2 * J[2];
      Jc[2] = d3 * J[2];
      Jc[3] = d1 * J[5];
      Jc[4] = d2 * J[5];
      Jc[5] = d3 * J[5];
      Jc[6] = d1 * J[8];
      Jc[7] = d2 * J[8];
      Jc[8] = d3 * J[8];

      // Compute the derivative of the displacement gradient
      // with respect to each of the three coordinate directions
      // displacement gradient: Ud = Ua*J
      TacsScalar Uda[9];
      Uda[0] = Ua[0] * Ja[0] + Ua[1] * Ja[3] + Ua[2] * Ja[6];
      Uda[3] = Ua[3] * Ja[0] + Ua[4] * Ja[3] + Ua[5] * Ja[6];
      Uda[6] = Ua[6] * Ja[0] + Ua[7] * Ja[3] + Ua[8] * Ja[6];

      Uda[1] = Ua[0] * Ja[1] + Ua[1] * Ja[4] + Ua[2] * Ja[7];
      Uda[4] = Ua[3] * Ja[1] + Ua[4] * Ja[4] + Ua[5] * Ja[7];
      Uda[7] = Ua[6] * Ja[1] + Ua[7] * Ja[4] + Ua[8] * Ja[7];

      Uda[2] = Ua[0] * Ja[2] + Ua[1] * Ja[5] + Ua[2] * Ja[8];
      Uda[5] = Ua[3] * Ja[2] + Ua[4] * Ja[5] + Ua[5] * Ja[8];
      Uda[8] = Ua[6] * Ja[2] + Ua[7] * Ja[5] + Ua[8] * Ja[8];

      TacsScalar Udb[9];
      Udb[0] = Ua[0] * Jb[0] + Ua[1] * Jb[3] + Ua[2] * Jb[6];
      Udb[3] = Ua[3] * Jb[0] + Ua[4] * Jb[3] + Ua[5] * Jb[6];
      Udb[6] = Ua[6] * Jb[0] + Ua[7] * Jb[3] + Ua[8] * Jb[6];

      Udb[1] = Ua[0] * Jb[1] + Ua[1] * Jb[4] + Ua[2] * Jb[7];
      Udb[4] = Ua[3] * Jb[1] + Ua[4] * Jb[4] + Ua[5] * Jb[7];
      Udb[7] = Ua[6] * Jb[1] + Ua[7] * Jb[4] + Ua[8] * Jb[7];

      Udb[2] = Ua[0] * Jb[2] + Ua[1] * Jb[5] + Ua[2] * Jb[8];
      Udb[5] = Ua[3] * Jb[2] + Ua[4] * Jb[5] + Ua[5] * Jb[8];
      Udb[8] = Ua[6] * Jb[2] + Ua[7] * Jb[5] + Ua[8] * Jb[8];

      TacsScalar Udc[9];
      Udc[0] = Ua[0] * Jc[0] + Ua[1] * Jc[3] + Ua[2] * Jc[6];
      Udc[3] = Ua[3] * Jc[0] + Ua[4] * Jc[3] + Ua[5] * Jc[6];
      Udc[6] = Ua[6] * Jc[0] + Ua[7] * Jc[3] + Ua[8] * Jc[6];

      Udc[1] = Ua[0] * Jc[1] + Ua[1] * Jc[4] + Ua[2] * Jc[7];
      Udc[4] = Ua[3] * Jc[1] + Ua[4] * Jc[4] + Ua[5] * Jc[7];
      Udc[7] = Ua[6] * Jc[1] + Ua[7] * Jc[4] + Ua[8] * Jc[7];

      Udc[2] = Ua[0] * Jc[2] + Ua[1] * Jc[5] + Ua[2] * Jc[8];
      Udc[5] = Ua[3] * Jc[2] + Ua[4] * Jc[5] + Ua[5] * Jc[8];
      Udc[8] = Ua[6] * Jc[2] + Ua[7] * Jc[5] + Ua[8] * Jc[8];

      sens[0] += (Uda[0] * eSens[0] + eSens[1] * Uda[4] + Uda[8] * eSens[2] +
                  (Uda[5] + Uda[7]) * eSens[3] + (Uda[2] + Uda[6]) * eSens[4] +
                  (Uda[1] + Uda[3]) * eSens[5]);
      sens++;

      sens[0] += (Udb[0] * eSens[0] + Udb[4] * eSens[1] + Udb[8] * eSens[2] +
                  (Udb[5] + Udb[7]) * eSens[3] + (Udb[2] + Udb[6]) * eSens[4] +
                  (Udb[1] + Udb[3]) * eSens[5]);
      sens++;

      sens[0] += (Udc[0] * eSens[0] + Udc[4] * eSens[1] + Udc[8] * eSens[2] +
                  (Udc[5] + Udc[7]) * eSens[3] + (Udc[2] + Udc[6]) * eSens[4] +
                  (Udc[1] + Udc[3]) * eSens[5]);
      sens++;

      Na++;
      Nb++;
      Nc++;
    }
  } else {
    // Compute the displacement gradient: Ud = Ua*J
    TacsScalar Ud[9];
    Ud[0] = Ua[0] * J[0] + Ua[1] * J[3] + Ua[2] * J[6];
    Ud[3] = Ua[3] * J[0] + Ua[4] * J[3] + Ua[5] * J[6];
    Ud[6] = Ua[6] * J[0] + Ua[7] * J[3] + Ua[8] * J[6];

    Ud[1] = Ua[0] * J[1] + Ua[1] * J[4] + Ua[2] * J[7];
    Ud[4] = Ua[3] * J[1] + Ua[4] * J[4] + Ua[5] * J[7];
    Ud[7] = Ua[6] * J[1] + Ua[7] * J[4] + Ua[8] * J[7];

    Ud[2] = Ua[0] * J[2] + Ua[1] * J[5] + Ua[2] * J[8];
    Ud[5] = Ua[3] * J[2] + Ua[4] * J[5] + Ua[5] * J[8];
    Ud[8] = Ua[6] * J[2] + Ua[7] * J[5] + Ua[8] * J[8];

    for (int i = 0; i < NUM_NODES; i++) {
      // JSens = -J*d(Xa)/dx*J
      TacsScalar d1 = -(Na[0] * J[0] + Nb[0] * J[3] + Nc[0] * J[6]);
      TacsScalar d2 = -(Na[0] * J[1] + Nb[0] * J[4] + Nc[0] * J[7]);
      TacsScalar d3 = -(Na[0] * J[2] + Nb[0] * J[5] + Nc[0] * J[8]);

      // The derivative of the inverse transformation matrix
      // in each coordinate direction
      TacsScalar Ja[9], Jb[9], Jc[9];
      Ja[0] = d1 * J[0];
      Ja[1] = d2 * J[0];
      Ja[2] = d3 * J[0];
      Ja[3] = d1 * J[3];
      Ja[4] = d2 * J[3];
      Ja[5] = d3 * J[3];
      Ja[6] = d1 * J[6];
      Ja[7] = d2 * J[6];
      Ja[8] = d3 * J[6];

      Jb[0] = d1 * J[1];
      Jb[1] = d2 * J[1];
      Jb[2] = d3 * J[1];
      Jb[3] = d1 * J[4];
      Jb[4] = d2 * J[4];
      Jb[5] = d3 * J[4];
      Jb[6] = d1 * J[7];
      Jb[7] = d2 * J[7];
      Jb[8] = d3 * J[7];

      Jc[0] = d1 * J[2];
      Jc[1] = d2 * J[2];
      Jc[2] = d3 * J[2];
      Jc[3] = d1 * J[5];
      Jc[4] = d2 * J[5];
      Jc[5] = d3 * J[5];
      Jc[6] = d1 * J[8];
      Jc[7] = d2 * J[8];
      Jc[8] = d3 * J[8];

      // Compute the derivative of the displacement gradient
      // with respect to each of the three coordinate directions
      // displacement gradient: Ud = Ua*J
      TacsScalar Uda[9];
      Uda[0] = Ua[0] * Ja[0] + Ua[1] * Ja[3] + Ua[2] * Ja[6];
      Uda[3] = Ua[3] * Ja[0] + Ua[4] * Ja[3] + Ua[5] * Ja[6];
      Uda[6] = Ua[6] * Ja[0] + Ua[7] * Ja[3] + Ua[8] * Ja[6];

      Uda[1] = Ua[0] * Ja[1] + Ua[1] * Ja[4] + Ua[2] * Ja[7];
      Uda[4] = Ua[3] * Ja[1] + Ua[4] * Ja[4] + Ua[5] * Ja[7];
      Uda[7] = Ua[6] * Ja[1] + Ua[7] * Ja[4] + Ua[8] * Ja[7];

      Uda[2] = Ua[0] * Ja[2] + Ua[1] * Ja[5] + Ua[2] * Ja[8];
      Uda[5] = Ua[3] * Ja[2] + Ua[4] * Ja[5] + Ua[5] * Ja[8];
      Uda[8] = Ua[6] * Ja[2] + Ua[7] * Ja[5] + Ua[8] * Ja[8];

      TacsScalar Udb[9];
      Udb[0] = Ua[0] * Jb[0] + Ua[1] * Jb[3] + Ua[2] * Jb[6];
      Udb[3] = Ua[3] * Jb[0] + Ua[4] * Jb[3] + Ua[5] * Jb[6];
      Udb[6] = Ua[6] * Jb[0] + Ua[7] * Jb[3] + Ua[8] * Jb[6];

      Udb[1] = Ua[0] * Jb[1] + Ua[1] * Jb[4] + Ua[2] * Jb[7];
      Udb[4] = Ua[3] * Jb[1] + Ua[4] * Jb[4] + Ua[5] * Jb[7];
      Udb[7] = Ua[6] * Jb[1] + Ua[7] * Jb[4] + Ua[8] * Jb[7];

      Udb[2] = Ua[0] * Jb[2] + Ua[1] * Jb[5] + Ua[2] * Jb[8];
      Udb[5] = Ua[3] * Jb[2] + Ua[4] * Jb[5] + Ua[5] * Jb[8];
      Udb[8] = Ua[6] * Jb[2] + Ua[7] * Jb[5] + Ua[8] * Jb[8];

      TacsScalar Udc[9];
      Udc[0] = Ua[0] * Jc[0] + Ua[1] * Jc[3] + Ua[2] * Jc[6];
      Udc[3] = Ua[3] * Jc[0] + Ua[4] * Jc[3] + Ua[5] * Jc[6];
      Udc[6] = Ua[6] * Jc[0] + Ua[7] * Jc[3] + Ua[8] * Jc[6];

      Udc[1] = Ua[0] * Jc[1] + Ua[1] * Jc[4] + Ua[2] * Jc[7];
      Udc[4] = Ua[3] * Jc[1] + Ua[4] * Jc[4] + Ua[5] * Jc[7];
      Udc[7] = Ua[6] * Jc[1] + Ua[7] * Jc[4] + Ua[8] * Jc[7];

      Udc[2] = Ua[0] * Jc[2] + Ua[1] * Jc[5] + Ua[2] * Jc[8];
      Udc[5] = Ua[3] * Jc[2] + Ua[4] * Jc[5] + Ua[5] * Jc[8];
      Udc[8] = Ua[6] * Jc[2] + Ua[7] * Jc[5] + Ua[8] * Jc[8];

      sens[0] +=
          ((Uda[0] + (Ud[0] * Uda[0] + Ud[3] * Uda[3] + Ud[6] * Uda[6])) *
               eSens[0] +
           (Uda[4] + (Ud[1] * Uda[1] + Ud[4] * Uda[4] + Ud[7] * Uda[7])) *
               eSens[1] +
           (Uda[8] + (Ud[2] * Uda[2] + Ud[5] * Uda[5] + Ud[8] * Uda[8])) *
               eSens[2] +
           (Uda[5] + Uda[7] +
            (Uda[1] * Ud[2] + Uda[4] * Ud[5] + Uda[7] * Ud[8] + Ud[1] * Uda[2] +
             Ud[4] * Uda[5] + Ud[7] * Uda[8])) *
               eSens[3] +
           (Uda[2] + Uda[6] +
            (Uda[0] * Ud[2] + Uda[3] * Ud[5] + Uda[6] * Ud[8] + Ud[0] * Uda[2] +
             Ud[3] * Uda[5] + Ud[6] * Uda[8])) *
               eSens[4] +
           (Uda[1] + Uda[3] +
            (Uda[0] * Ud[1] + Uda[3] * Ud[4] + Uda[6] * Ud[7] + Ud[0] * Uda[1] +
             Ud[3] * Uda[4] + Ud[6] * Uda[7])) *
               eSens[5]);
      sens++;

      sens[0] +=
          ((Udb[0] + (Ud[0] * Udb[0] + Ud[3] * Udb[3] + Ud[6] * Udb[6])) *
               eSens[0] +
           (Udb[4] + (Ud[1] * Udb[1] + Ud[4] * Udb[4] + Ud[7] * Udb[7])) *
               eSens[1] +
           (Udb[8] + (Ud[2] * Udb[2] + Ud[5] * Udb[5] + Ud[8] * Udb[8])) *
               eSens[2] +
           (Udb[5] + Udb[7] +
            (Udb[1] * Ud[2] + Udb[4] * Ud[5] + Udb[7] * Ud[8] + Ud[1] * Udb[2] +
             Ud[4] * Udb[5] + Ud[7] * Udb[8])) *
               eSens[3] +
           (Udb[2] + Udb[6] +
            (Udb[0] * Ud[2] + Udb[3] * Ud[5] + Udb[6] * Ud[8] + Ud[0] * Udb[2] +
             Ud[3] * Udb[5] + Ud[6] * Udb[8])) *
               eSens[4] +
           (Udb[1] + Udb[3] +
            (Udb[0] * Ud[1] + Udb[3] * Ud[4] + Udb[6] * Ud[7] + Ud[0] * Udb[1] +
             Ud[3] * Udb[4] + Ud[6] * Udb[7])) *
               eSens[5]);
      sens++;

      sens[0] +=
          ((Udc[0] + (Ud[0] * Udc[0] + Ud[3] * Udc[3] + Ud[6] * Udc[6])) *
               eSens[0] +
           (Udc[4] + (Ud[1] * Udc[1] + Ud[4] * Udc[4] + Ud[7] * Udc[7])) *
               eSens[1] +
           (Udc[8] + (Ud[2] * Udc[2] + Ud[5] * Udc[5] + Ud[8] * Udc[8])) *
               eSens[2] +
           (Udc[5] + Udc[7] +
            (Udc[1] * Ud[2] + Udc[4] * Ud[5] + Udc[7] * Ud[8] + Ud[1] * Udc[2] +
             Ud[4] * Udc[5] + Ud[7] * Udc[8])) *
               eSens[3] +
           (Udc[2] + Udc[6] +
            (Udc[0] * Ud[2] + Udc[3] * Ud[5] + Udc[6] * Ud[8] + Ud[0] * Udc[2] +
             Ud[3] * Udc[5] + Ud[6] * Udc[8])) *
               eSens[4] +
           (Udc[1] + Udc[3] +
            (Udc[0] * Ud[1] + Udc[3] * Ud[4] + Udc[6] * Ud[7] + Ud[0] * Udc[1] +
             Ud[3] * Udc[4] + Ud[6] * Udc[7])) *
               eSens[5]);
      sens++;

      Na++;
      Nb++;
      Nc++;
    }
  }
}

/*
  Compute the kinetic and strain energies in the element

  output:
  Te:      the kinetic energy
  Pe:      the potential energy (including strain energy)

  input:
  Xpts:    the nodal locations
  vars:    the nodal variables
  dvars:   the time derivatives of the nodal variables

 */
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::computeEnergies(double time, TacsScalar *_Te,
                                               TacsScalar *_Pe,
                                               const TacsScalar Xpts[],
                                               const TacsScalar vars[],
                                               const TacsScalar dvars[]) {
  // Compute the kinetic and potential energy
  TacsScalar Te = 0.0, Pe = 0.0;

  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

  // Get the number of quadrature points
  int numGauss = getNumGaussPts();

  for (int n = 0; n < numGauss; n++) {
    // Retrieve the quadrature points and weight
    double pt[3];
    double weight = getGaussWtsPts(n, pt);

    // Compute the element shape functions
    getShapeFunctions(pt, N, Na, Nb, Nc);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[9];
    solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[9];
    TacsScalar h = FElibrary::jacobian3d(Xa, J);
    h = h * weight;

    // Compute the strain
    TacsScalar strain[NUM_STRESSES];
    evalStrain(strain, J, Na, Nb, Nc, vars);

    // Compute the corresponding stress
    TacsScalar stress[NUM_STRESSES];
    stiff->calculateStress(pt, strain, stress);

    // Compute the contribution from the potential energy
    Pe +=
        0.5 * h *
        (stress[0] * strain[0] + stress[1] * strain[1] + stress[2] * strain[2] +
         stress[3] * strain[3] + stress[4] * strain[4] + stress[5] * strain[5]);

    // Get value of the mass/area at this point
    TacsScalar mass;
    stiff->getPointwiseMass(pt, &mass);

    // Compute the contribution from the kinetic energy
    TacsScalar dUdt[3];
    getDisplacement(dUdt, N, dvars);
    Te += 0.5 * h * mass *
          (dUdt[0] * dUdt[0] + dUdt[1] * dUdt[1] + dUdt[2] * dUdt[2]);
  }

  // Set the output values
  *_Te = Te;
  *_Pe = Pe;
}

/*
  Compute the residual contribution from this element not including
  any externally applied loads.

  output:
  res:     the element residual

  input:
  vars:    the element variables
  Xpts:    the element nodal locations in R^{3}
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::addResidual(double time, TacsScalar res[],
                                           const TacsScalar Xpts[],
                                           const TacsScalar vars[],
                                           const TacsScalar dvars[],
                                           const TacsScalar ddvars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

  // Get the number of quadrature points
  int numGauss = getNumGaussPts();

  for (int n = 0; n < numGauss; n++) {
    // Retrieve the quadrature points and weight
    double pt[3];
    double weight = getGaussWtsPts(n, pt);

    // Compute the element shape functions
    getShapeFunctions(pt, N, Na, Nb, Nc);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[9];
    solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[9];
    TacsScalar h = FElibrary::jacobian3d(Xa, J);
    h = h * weight;

    // Compute the strain
    TacsScalar strain[NUM_STRESSES];
    evalStrain(strain, J, Na, Nb, Nc, vars);

    // Compute the corresponding stress
    TacsScalar stress[NUM_STRESSES];
    stiff->calculateStress(pt, strain, stress);

    // Get the derivative of the strain with respect to the nodal
    // displacements
    getBmat(B, J, Na, Nb, Nc, vars);

    TacsScalar *b = B;
    for (int i = 0; i < NUM_VARIABLES; i++) {
      res[i] += h * (b[0] * stress[0] + b[1] * stress[1] + b[2] * stress[2] +
                     b[3] * stress[3] + b[4] * stress[4] + b[5] * stress[5]);
      b += NUM_STRESSES;
    }

    // Get value of the mass/area at this point
    TacsScalar mass;
    stiff->getPointwiseMass(pt, &mass);

    // Add the contribution from the inertial terms
    TacsScalar d2Udt2[3];
    getDisplacement(d2Udt2, N, ddvars);
    for (int i = 0; i < NUM_NODES; i++) {
      res[3 * i] += h * mass * N[i] * d2Udt2[0];
      res[3 * i + 1] += h * mass * N[i] * d2Udt2[1];
      res[3 * i + 2] += h * mass * N[i] * d2Udt2[2];
    }
  }
}

/*
  Get the Jacobian of the governing equations- the exact Jacobian of the
  residual expressions.

  output:
  mat:     the element Jacobian

  input:
  alpha:   coefficient of the time-independent terms
  beta:    coefficient of the first time derivatives
  gamma:   coefficient of the second time derivatives
  Xpts:    the element nodal locations in R^{3}
  vars:    the element variables
  dvars:   time derivative of the element variables
  ddvars:  second time derivative of the element variables
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::addJacobian(
    double time, TacsScalar mat[], double alpha, double beta, double gamma,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

  // Get the number of quadrature points
  int numGauss = getNumGaussPts();

  for (int n = 0; n < numGauss; n++) {
    // Retrieve the quadrature points and weight
    double pt[3];
    double weight = getGaussWtsPts(n, pt);

    // Compute the element shape functions
    getShapeFunctions(pt, N, Na, Nb, Nc);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[9];
    solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[9];
    TacsScalar h = FElibrary::jacobian3d(Xa, J);
    h = h * weight;

    if (alpha != 0.0) {
      // Compute the strain
      TacsScalar strain[NUM_STRESSES];
      evalStrain(strain, J, Na, Nb, Nc, vars);

      // Compute the corresponding stress
      TacsScalar stress[NUM_STRESSES];
      stiff->calculateStress(pt, strain, stress);

      // Add the stress times the second derivative of the strain
      addGeoStiffness(mat, alpha * h, stress, J, Na, Nb, Nc);

      // Get the derivative of the strain with respect to the nodal
      // displacements
      getBmat(B, J, Na, Nb, Nc, vars);

      // Fill-in the upper portion of the matrix
      TacsScalar *bj = B;
      for (int j = 0; j < NUM_VARIABLES; j++) {
        // Compute the stress at the given point
        TacsScalar bs[NUM_STRESSES];
        stiff->calculateStress(pt, bj, bs);

        TacsScalar *bi = B;
        for (int i = 0; i <= j; i++) {
          mat[i + j * NUM_VARIABLES] +=
              alpha * h *
              (bi[0] * bs[0] + bi[1] * bs[1] + bi[2] * bs[2] + bi[3] * bs[3] +
               bi[4] * bs[4] + bi[5] * bs[5]);
          bi += NUM_STRESSES;
        }
        bj += NUM_STRESSES;
      }
    }

    if (gamma != 0.0) {
      // Get value of the mass/area at this point
      TacsScalar mass;
      stiff->getPointwiseMass(pt, &mass);

      // Add the contributions from the stiffness matrix
      TacsScalar scale = gamma * h * mass;
      for (int j = 0; j < NUM_NODES; j++) {
        for (int i = 0; i <= j; i++) {
          mat[3 * i + 3 * j * NUM_VARIABLES] += scale * N[i] * N[j];
          mat[3 * i + 1 + (3 * j + 1) * NUM_VARIABLES] += scale * N[i] * N[j];
          mat[3 * i + 2 + (3 * j + 2) * NUM_VARIABLES] += scale * N[i] * N[j];
        }
      }
    }
  }

  // Apply symmetry to the matrix
  for (int j = 0; j < NUM_VARIABLES; j++) {
    for (int i = 0; i < j; i++) {
      mat[j + i * NUM_VARIABLES] = mat[i + j * NUM_VARIABLES];
    }
  }
}

/*
  Add the product of the adjoint vector times the derivative of the
  residuals multiplied by a scalar to the given derivative vector.

  output:
  dvSens:  the derivative of the values w.r.t. the nodes

  input:
  dvLen:   the design variable array length
  scale:   scale the derivative by this scalar
  psi:     the element adjoint variables
  Xpts:    the element nodal locations in R^{3}
  vars:    the element variables
  dvars:   time derivative of the element variables
  ddvars:  second time derivative of the element variables
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::addAdjResProduct(
    double time, double scale, TacsScalar dvSens[], int dvLen,
    const TacsScalar psi[], const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

  // Get the number of quadrature points
  int numGauss = getNumGaussPts();

  for (int n = 0; n < numGauss; n++) {
    // Retrieve the quadrature points and weights
    double pt[3];
    double weight = getGaussWtsPts(n, pt);

    // Compute the element shape functions
    getShapeFunctions(pt, N, Na, Nb, Nc);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[9];
    solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[9];
    TacsScalar h = FElibrary::jacobian3d(Xa, J);
    h = h * weight;

    // Compute the strain
    TacsScalar strain[NUM_STRESSES];
    evalStrain(strain, J, Na, Nb, Nc, vars);

    // Get the derivative of the strain with respect to the nodal
    // displacements
    getBmat(B, J, Na, Nb, Nc, vars);

    // Compute the product of psi^{T}*B^{T}
    TacsScalar bpsi[NUM_STRESSES];
    memset(bpsi, 0, NUM_STRESSES * sizeof(TacsScalar));

    TacsScalar *b = B;
    const TacsScalar *ps = psi;
    for (int i = 0; i < NUM_VARIABLES; i++) {
      bpsi[0] += ps[0] * b[0];
      bpsi[1] += ps[0] * b[1];
      bpsi[2] += ps[0] * b[2];
      bpsi[3] += ps[0] * b[3];
      bpsi[4] += ps[0] * b[4];
      bpsi[5] += ps[0] * b[5];
      b += NUM_STRESSES;
      ps++;
    }

    // Add the term: alpha*psi^{T}*B^{T}*dC/dx*strain to the vector
    // dvSens - Note that this is much more efficient than computing
    // the terms component by component
    stiff->addStressDVSens(pt, strain, scale * h, bpsi, dvSens, dvLen);
  }
}

/*
  Evaluate the derivative of the element residuals with respect to the
  nodal coordinates such that fXptSens += scale*psi^{T}*dR/d(Xpts)

  To compute this term, first consider the inner product

  psi^{T}*R = sum (h*psi^{T}*B^{T}*D*strain)

  This derivative can be broken down into three contributions:

  psi^{T}*d(R)/d(Xpt) =
  sum( d(h)/d(Xpt)*(psi^{T}*B^{T}*D*strain) +
  .    psi^{T}*B^{T}*D*(d(strain)/d(Xpt)) +
  .    strain^{T}*D*d(psi^{T}*B)/d(Xpt))

  which involve the derivatives: d(h)/d(Xpt), d(v^{T}*strain)/d(Xpt)
  and d(w^{T}*B*v)/d(Xpt)

  output:
  fXptSens:  the derivative of the residuals w.r.t. the element nodes

  input:
  time:      the simulation time
  scale:     scale the result by this value
  psi:       the adjoint variables
  Xpts:      the element nodal locations
  vars:      the element variables
  dvars:     the first time derivative of the element variables
  ddvars:    the second time derivative of the element variables
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::addAdjResXptProduct(
    double time, double scale, TacsScalar fXptSens[], const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[]) {
  /*
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  // Get the number of quadrature points
  int numGauss = getNumGaussPts();

  for ( int n = 0; n < numGauss; n++ ){
    // Retrieve the quadrature points and weight
    double pt[3];
    double weight = getGaussWtsPts(n, pt);

    // Compute the element shape functions
    getShapeFunctions(pt, N, Na, Nb, Nc);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[9];
    solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[9];
    TacsScalar h = FElibrary::jacobian3d(Xa, J);
    h = h*weight;

    // Compute the derivative of the determinant w.r.t. nodes
    TacsScalar hXptSens[3*NUM_NODES];
    getDetJacobianXptSens(hXptSens, pt, Xpts);

    // Compute the strain
    TacsScalar strain[NUM_STRESSES];
    evalStrain(strain, J, Na, Nb, Nc, vars);

    // Compute the corresponding stress
    TacsScalar stress[NUM_STRESSES];
    stiff->calculateStress(pt, strain, stress);

    // Get the derivative of the strain with respect to the nodal
    // displacements
    getBmat(B, J, Na, Nb, Nc, vars);

    TacsScalar *b = B;
    for ( int i = 0; i < NUM_VARIABLES; i++ ){
      res[i] += h*(b[0]*stress[0] + b[1]*stress[1] + b[2]*stress[2] +
                   b[3]*stress[3] + b[4]*stress[4] + b[5]*stress[5]);
      b += NUM_STRESSES;
    }

    // Get value of the mass/area at this point
    TacsScalar mass;
    stiff->getPointwiseMass(pt, &mass);

    // Add the contribution from the inertial terms
    TacsScalar d2Udt2[3];
    getDisplacement(d2Udt2, N, ddvars);
    for ( int i = 0; i < NUM_NODES; i++ ){
      res[3*i] += h*mass*N[i]*d2Udt2[0];
      res[3*i+1] += h*mass*N[i]*d2Udt2[1];
      res[3*i+2] += h*mass*N[i]*d2Udt2[2];
    }
  }
  */
}

/*
  Add the derivative of the inner product of the stiffness or mass
  matrix with respect to the design variables to a design variable
  vector. This is much more efficient than computing the derivative of
  the stiffness/mass matrix, then computing the product for each
  design variable.

  input:
  matType:     the matrix type (e.g. MASS_MATRIX)
  scale:       the scaling factor
  dvLen:       the length of the design variable vector
  psi:         the left inner-product vector
  phi:         the right inner-product vector
  Xpts:        the nodal locations
  vars:        the state variable values

  output:
  dvSens:      vector of the design sensitivity
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::addMatDVSensInnerProduct(
    ElementMatrixType matType, double scale, TacsScalar dvSens[], int dvLen,
    const TacsScalar psi[], const TacsScalar phi[], const TacsScalar Xpts[],
    const TacsScalar vars[]) {
  if (matType == STIFFNESS_MATRIX) {
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

    // The derivative of the stress with respect to the strain
    TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

    // Get the number of quadrature points
    int numGauss = getNumGaussPts();

    for (int n = 0; n < numGauss; n++) {
      // Retrieve the quadrature points and weights
      double pt[3];
      double weight = getGaussWtsPts(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb, Nc);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[9];
      TacsScalar h = FElibrary::jacobian3d(Xa, J);
      h = h * weight;

      // Get the derivative of the strain with respect to the nodal
      // displacements
      getBmat(B, J, Na, Nb, Nc, vars);

      // Compute the product of psi^{T}*B^{T}
      TacsScalar bpsi[NUM_STRESSES], bphi[NUM_STRESSES];
      memset(bpsi, 0, NUM_STRESSES * sizeof(TacsScalar));
      memset(bphi, 0, NUM_STRESSES * sizeof(TacsScalar));

      TacsScalar *b = B;
      const TacsScalar *ps = psi, *ph = phi;
      for (int i = 0; i < 3 * NUM_NODES; i++) {
        bpsi[0] += ps[0] * b[0];
        bpsi[1] += ps[0] * b[1];
        bpsi[2] += ps[0] * b[2];
        bpsi[3] += ps[0] * b[3];
        bpsi[4] += ps[0] * b[4];
        bpsi[5] += ps[0] * b[5];

        bphi[0] += ph[0] * b[0];
        bphi[1] += ph[0] * b[1];
        bphi[2] += ph[0] * b[2];
        bphi[3] += ph[0] * b[3];
        bphi[4] += ph[0] * b[4];
        bphi[5] += ph[0] * b[5];

        b += 6;
        ps++;
        ph++;
      }

      // Add the result to the design variable vector
      stiff->addStressDVSens(pt, bphi, scale * h, bpsi, dvSens, dvLen);
    }
  } else if (matType == GEOMETRIC_STIFFNESS_MATRIX) {
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

    // Get the number of quadrature points
    int numGauss = getNumGaussPts();

    for (int n = 0; n < numGauss; n++) {
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getGaussWtsPts(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb, Nc);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[9];
      TacsScalar h = FElibrary::jacobian3d(Xa, J);
      h = h * weight;

      // Compute the strain derived from load path
      TacsScalar strain[NUM_STRESSES];
      evalStrain(strain, J, Na, Nb, Nc, vars);

      // Compute dN/dx*phi
      // Compute the product of psi^{T}*G^{T} and G*phi
      TacsScalar gpsi[9], gphi[9];
      memset(gpsi, 0, 9 * sizeof(TacsScalar));
      memset(gphi, 0, 9 * sizeof(TacsScalar));
      for (int j = 0; j < NUM_NODES; j++) {
        TacsScalar Dx = Na[j] * J[0] + Nb[j] * J[3] + Nc[j] * J[6];
        TacsScalar Dy = Na[j] * J[1] + Nb[j] * J[4] + Nc[j] * J[7];
        TacsScalar Dz = Na[j] * J[2] + Nb[j] * J[5] + Nc[j] * J[8];

        gpsi[0] += Dx * psi[3 * j];
        gphi[0] += Dx * phi[3 * j];
        gpsi[1] += Dx * psi[3 * j + 1];
        gphi[1] += Dx * phi[3 * j + 1];
        gpsi[2] += Dx * psi[3 * j + 2];
        gphi[2] += Dx * phi[3 * j + 2];

        gpsi[3] += Dy * psi[3 * j];
        gphi[3] += Dy * phi[3 * j];
        gpsi[4] += Dy * psi[3 * j + 1];
        gphi[4] += Dy * phi[3 * j + 1];
        gpsi[5] += Dy * psi[3 * j + 2];
        gphi[5] += Dy * phi[3 * j + 2];

        gpsi[6] += Dz * psi[3 * j];
        gphi[6] += Dz * phi[3 * j];
        gpsi[7] += Dz * psi[3 * j + 1];
        gphi[7] += Dz * phi[3 * j + 1];
        gpsi[8] += Dz * psi[3 * j + 2];
        gphi[8] += Dz * phi[3 * j + 2];
      }

      TacsScalar sumN[NUM_STRESSES];
      memset(sumN, 0, NUM_STRESSES * sizeof(TacsScalar));
      TacsScalar *gs = gpsi, *gh = gphi;
      for (int j = 0; j < 3; j++) {
        sumN[j] = gs[0] * gh[0] + gs[1] * gh[1] + gs[2] * gh[2];
        gs += 3;
        gh += 3;
      }
      for (int j = 0; j < 3; j++) {
        sumN[3] += gpsi[3 + j] * gphi[6 + j] + gpsi[6 + j] * gphi[3 + j];
      }
      for (int j = 0; j < 3; j++) {
        sumN[4] += gpsi[j] * gphi[6 + j] + gpsi[6 + j] * gphi[j];
      }
      for (int j = 0; j < 3; j++) {
        sumN[5] += gpsi[j] * gphi[3 + j] + gpsi[3 + j] * gphi[j];
      }

      // Add the result to the design variable vector
      stiff->addStressDVSens(pt, sumN, scale * h, strain, dvSens, dvLen);
    }
  } else if (matType == MASS_MATRIX) {
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

    // Get the number of quadrature points
    int numGauss = getNumGaussPts();

    for (int n = 0; n < numGauss; n++) {
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getGaussWtsPts(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb, Nc);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[9];
      TacsScalar h = FElibrary::jacobian3d(Xa, J);
      h = h * weight;

      // Compute the nodal accelerations at the quadrature point
      TacsScalar upsi[3], uphi[3];
      upsi[0] = upsi[1] = upsi[2] = 0.0;
      uphi[0] = uphi[1] = uphi[2] = 0.0;

      double *ns = N;
      const TacsScalar *ps = psi, *ph = phi;
      for (int i = 0; i < NUM_NODES; i++) {
        upsi[0] += ns[0] * ps[0];
        upsi[1] += ns[0] * ps[1];
        upsi[2] += ns[0] * ps[2];

        uphi[0] += ns[0] * ph[0];
        uphi[1] += ns[0] * ph[1];
        uphi[2] += ns[0] * ph[2];

        ps += 3;
        ph += 3;
        ns++;
      }

      // Add the result to the design variable vector
      TacsScalar rho_alpha =
          scale * h *
          (upsi[0] * uphi[0] + upsi[1] * uphi[1] + upsi[2] * uphi[2]);

      stiff->addPointwiseMassDVSens(pt, &rho_alpha, dvSens, dvLen);
    }
  }
}

/*
  Evaluate the derivative of the inner product of the stiffness or mass
  matrix with respect to the state variables or load path variables. This is
  much more efficient than computing the derivative of the stiffness/mass
  matrix, then computing the product for each state variable.

  input:
  matType:     the matrix type (e.g. MASS_MATRIX)
  res:         the derivative inner product w.r.t. the state variables
  psi:         the left inner-product vector
  phi:         the right inner-product vector
  Xpts:        the nodal locations
  vars:        the state variable values

  output:
  dvSens:      vector of the design sensitivity
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getMatSVSensInnerProduct(
    ElementMatrixType matType, TacsScalar res[], const TacsScalar psi[],
    const TacsScalar phi[], const TacsScalar Xpts[], const TacsScalar vars[]) {
  if (matType == GEOMETRIC_STIFFNESS_MATRIX) {
    memset(res, 0, NUM_VARIABLES * sizeof(TacsScalar));

    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

    // Get the number of quadrature points
    int numGauss = getNumGaussPts();

    for (int n = 0; n < numGauss; n++) {
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getGaussWtsPts(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb, Nc);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[9];
      TacsScalar h = FElibrary::jacobian3d(Xa, J);
      h = h * weight;

      // Compute the strain derived from load path
      TacsScalar strain[NUM_STRESSES];
      evalStrain(strain, J, Na, Nb, Nc, vars);

      // The derivative of the stress with respect to the strain
      TacsScalar B[NUM_STRESSES * NUM_VARIABLES];
      // Get the derivative of the strain with respect to the nodal
      // displacements
      getBmat(B, J, Na, Nb, Nc, vars);

      // Compute dN/dx*phi
      // Compute the product of psi^{T}*G^{T} and G*phi
      TacsScalar gpsi[9], gphi[9];
      memset(gpsi, 0, 9 * sizeof(TacsScalar));
      memset(gphi, 0, 9 * sizeof(TacsScalar));
      for (int j = 0; j < NUM_NODES; j++) {
        TacsScalar Dx = Na[j] * J[0] + Nb[j] * J[3] + Nc[j] * J[6];
        TacsScalar Dy = Na[j] * J[1] + Nb[j] * J[4] + Nc[j] * J[7];
        TacsScalar Dz = Na[j] * J[2] + Nb[j] * J[5] + Nc[j] * J[8];

        gpsi[0] += Dx * psi[3 * j];
        gphi[0] += Dx * phi[3 * j];
        gpsi[1] += Dx * psi[3 * j + 1];
        gphi[1] += Dx * phi[3 * j + 1];
        gpsi[2] += Dx * psi[3 * j + 2];
        gphi[2] += Dx * phi[3 * j + 2];

        gpsi[3] += Dy * psi[3 * j];
        gphi[3] += Dy * phi[3 * j];
        gpsi[4] += Dy * psi[3 * j + 1];
        gphi[4] += Dy * phi[3 * j + 1];
        gpsi[5] += Dy * psi[3 * j + 2];
        gphi[5] += Dy * phi[3 * j + 2];

        gpsi[6] += Dz * psi[3 * j];
        gphi[6] += Dz * phi[3 * j];
        gpsi[7] += Dz * psi[3 * j + 1];
        gphi[7] += Dz * phi[3 * j + 1];
        gpsi[8] += Dz * psi[3 * j + 2];
        gphi[8] += Dz * phi[3 * j + 2];
      }

      TacsScalar sumN[NUM_STRESSES];
      memset(sumN, 0, NUM_STRESSES * sizeof(TacsScalar));
      TacsScalar *gs = gpsi, *gh = gphi;
      for (int j = 0; j < 3; j++) {
        sumN[j] = gs[0] * gh[0] + gs[1] * gh[1] + gs[2] * gh[2];
        gs += 3;
        gh += 3;
      }
      for (int j = 0; j < 3; j++) {
        sumN[3] += gpsi[3 + j] * gphi[6 + j] + gpsi[6 + j] * gphi[3 + j];
      }
      for (int j = 0; j < 3; j++) {
        sumN[4] += gpsi[j] * gphi[6 + j] + gpsi[6 + j] * gphi[j];
      }
      for (int j = 0; j < 3; j++) {
        sumN[5] += gpsi[j] * gphi[3 + j] + gpsi[3 + j] * gphi[j];
      }

      // Get D*B
      TacsScalar DB[NUM_STRESSES * NUM_VARIABLES];
      memset(DB, 0, NUM_STRESSES * NUM_VARIABLES * sizeof(TacsScalar));
      for (int j = 0; j < NUM_VARIABLES; j++) {
        stiff->calculateStress(pt, &B[NUM_STRESSES * j], &DB[NUM_STRESSES * j]);
      }

      TacsScalar *db = DB;
      for (int j = 0; j < NUM_VARIABLES; j++) {
        res[j] += h * (db[0] * sumN[0] + db[1] * sumN[1] + db[2] * sumN[2] +
                       db[3] * sumN[3] + db[4] * sumN[4] + db[5] * sumN[5]);
        db += NUM_STRESSES;
      }
    }
  }
}

/*
  Get the element matrix of the specified type (e.g. mass matrix)
  from the element.

  output:
  mat:         the element matrix of the specified type

  input:
  matType:     the matrix type (e.g. MASS_MATRIX)
  scaleFactor: scale factor such that mat = scaleFactor*M
  vars:        the element variables
  Xpts:        the nodal coordinates in R^{3}
  matOr:       the matrix orientation either NORMAL or TRANSPOSE
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getMatType(ElementMatrixType matType,
                                          TacsScalar mat[],
                                          const TacsScalar Xpts[],
                                          const TacsScalar vars[]) {
  memset(mat, 0, NUM_VARIABLES * NUM_VARIABLES * sizeof(TacsScalar));

  // The mass matrix
  if (matType == MASS_MATRIX) {
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

    // Get the number of quadrature points
    int numGauss = getNumGaussPts();

    for (int n = 0; n < numGauss; n++) {
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getGaussWtsPts(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb, Nc);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[9];
      TacsScalar h = FElibrary::jacobian3d(Xa, J);
      h = h * weight;

      // Get the pointwise mass
      TacsScalar ptmass[3];
      stiff->getPointwiseMass(pt, ptmass);

      // Fill-in the upper-portion of the matrix
      for (int j = 0; j < NUM_NODES; j++) {
        for (int i = 0; i <= j; i++) {
          TacsScalar d = h * ptmass[0] * N[i] * N[j];

          mat[3 * i + 3 * j * NUM_VARIABLES] += d;
          mat[3 * i + 1 + (3 * j + 1) * NUM_VARIABLES] += d;
          mat[3 * i + 2 + (3 * j + 2) * NUM_VARIABLES] += d;
        }
      }
    }

    // Apply symmetry to the matrix
    for (int j = 0; j < NUM_VARIABLES; j++) {
      for (int i = 0; i < j; i++) {
        mat[j + i * NUM_VARIABLES] = mat[i + j * NUM_VARIABLES];
      }
    }
  } else if (matType == GEOMETRIC_STIFFNESS_MATRIX) {
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

    // Get the number of quadrature points
    int numGauss = getNumGaussPts();

    for (int n = 0; n < numGauss; n++) {
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getGaussWtsPts(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb, Nc);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[9];
      TacsScalar h = FElibrary::jacobian3d(Xa, J);
      h = h * weight;

      // Compute the strain derived from load path
      TacsScalar strain[NUM_STRESSES];
      evalStrain(strain, J, Na, Nb, Nc, vars);

      // Compute the corresponding stress
      TacsScalar stress[NUM_STRESSES];
      stiff->calculateStress(pt, strain, stress);

      for (int j = 0; j < NUM_NODES; j++) {
        TacsScalar Dxj = Na[j] * J[0] + Nb[j] * J[3] + Nc[j] * J[6];
        TacsScalar Dyj = Na[j] * J[1] + Nb[j] * J[4] + Nc[j] * J[7];
        TacsScalar Dzj = Na[j] * J[2] + Nb[j] * J[5] + Nc[j] * J[8];

        for (int i = 0; i < NUM_NODES; i++) {
          TacsScalar Dxi = Na[i] * J[0] + Nb[i] * J[3] + Nc[i] * J[6];
          TacsScalar Dyi = Na[i] * J[1] + Nb[i] * J[4] + Nc[i] * J[7];
          TacsScalar Dzi = Na[i] * J[2] + Nb[i] * J[5] + Nc[i] * J[8];

          // Add the contributions to the stiffness matrix
          TacsScalar scale =
              h * (stress[0] * Dxi * Dxj + stress[1] * Dyi * Dyj +
                   stress[2] * Dzi * Dzj + stress[3] * (Dyi * Dzj + Dyj * Dzi) +
                   stress[4] * (Dxi * Dzj + Dxj * Dzi) +
                   stress[5] * (Dxi * Dyj + Dxj * Dyi));

          mat[3 * i + 3 * j * NUM_VARIABLES] += scale;
          mat[3 * i + 1 + (3 * j + 1) * NUM_VARIABLES] += scale;
          mat[3 * i + 2 + (3 * j + 2) * NUM_VARIABLES] += scale;
        }
      }
    }
  } else if (matType == STIFFNESS_MATRIX) {
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

    // The derivative of the stress with respect to the strain
    TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

    // Get the number of quadrature points
    int numGauss = getNumGaussPts();

    for (int n = 0; n < numGauss; n++) {
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getGaussWtsPts(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb, Nc);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[9];
      TacsScalar h = FElibrary::jacobian3d(Xa, J);
      h = h * weight;

      // Compute the strain
      TacsScalar strain[NUM_STRESSES];
      evalStrain(strain, J, Na, Nb, Nc, vars);

      // Compute the corresponding stress
      TacsScalar stress[NUM_STRESSES];
      stiff->calculateStress(pt, strain, stress);

      // Get the derivative of the strain with respect to the nodal
      // displacements
      getBmat(B, J, Na, Nb, Nc, vars);

      // Fill-in the upper portion of the matrix
      TacsScalar *bj = B;
      for (int j = 0; j < NUM_VARIABLES; j++) {
        // Compute the stress at the given point
        TacsScalar bs[NUM_STRESSES];
        stiff->calculateStress(pt, bj, bs);

        TacsScalar *bi = B;
        for (int i = 0; i <= j; i++) {
          mat[i + j * NUM_VARIABLES] +=
              h * (bi[0] * bs[0] + bi[1] * bs[1] + bi[2] * bs[2] +
                   bi[3] * bs[3] + bi[4] * bs[4] + bi[5] * bs[5]);
          bi += NUM_STRESSES;
        }
        bj += NUM_STRESSES;
      }
    }

    // Apply symmetry to the matrix
    for (int j = 0; j < NUM_VARIABLES; j++) {
      for (int i = 0; i < j; i++) {
        mat[j + i * NUM_VARIABLES] = mat[i + j * NUM_VARIABLES];
      }
    }
  }
}

/*
  Evaluate the determinant of the Jacobian for numerical integration

  returns: the determinant of the Jacobian

  input:
  pt:    the parametric point within the element
  Xpts:  the element nodes
*/
template <int NUM_NODES>
TacsScalar TACS3DElement<NUM_NODES>::getDetJacobian(const double pt[],
                                                    const TacsScalar Xpts[]) {
  // Compute the element shape functions
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of the shape functions w.r.t. the
  // parametric locations
  TacsScalar X[3], Xa[9];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

  return FElibrary::jacobian3d(Xa);
}

/*
  Evaluate the derivative of the determinant of the Jacobian with
  respect to the element nodal locations

  output:
  hXptSens: the derivative of the determinant w.r.t. the nodal locations

  returns:  the determinant of the Jacobian

  input:
  pt:    the parametric point within the element
  Xpts:  the element nodes
*/
template <int NUM_NODES>
TacsScalar TACS3DElement<NUM_NODES>::getDetJacobianXptSens(
    TacsScalar *hXptSens, const double pt[], const TacsScalar Xpts[]) {
  // Compute the element shape functions
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of the shape functions w.r.t. the
  // parametric locations
  TacsScalar X[3], Xa[9];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

  // Evaluate the determinant of the Jacobian
  TacsScalar J[9];
  TacsScalar h = FElibrary::jacobian3d(Xa, J);

  for (int i = 0; i < NUM_NODES; i++) {
    for (int k = 0; k < 3; k++) {
      TacsScalar XaSens[9];
      XaSens[0] = XaSens[1] = XaSens[2] = 0.0;
      XaSens[3] = XaSens[4] = XaSens[5] = 0.0;
      XaSens[6] = XaSens[7] = XaSens[8] = 0.0;
      XaSens[3 * k] = Na[i];
      XaSens[3 * k + 1] = Nb[i];
      XaSens[3 * k + 2] = Nc[i];

      FElibrary::jacobian3dSens(Xa, XaSens, &hXptSens[0]);
      hXptSens++;
    }
  }

  return h;
}

/*
  Evaluate the strain at the specified point using the provided set of
  variables

  output:
  strain:   the strain evaluate at the specific parametric point

  input:
  vars:     the element variable values
  Xpts:     the element nodal locations
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getStrain(TacsScalar strain[], const double pt[],
                                         const TacsScalar Xpts[],
                                         const TacsScalar vars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

  // Compute the element shape functions
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of X with respect to the coordinate directions
  TacsScalar X[3], Xa[9];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

  // Compute the determinant of Xa and the transformation
  TacsScalar J[9];
  FElibrary::jacobian3d(Xa, J);

  // Compute the strain
  evalStrain(strain, J, Na, Nb, Nc, vars);
}

/*
  Compute the derivative of the point-wise strain multiplied by a
  specified vector with respect to the element variables and add the
  result to an array. This can be used to evaluate the derivative of a
  function of interest with respect to the element variables.

  output:
  sens:        the output array - same length as the number of elem variables

  input:
  pt:          parametric point used to evaluate the derivative [-1, 1]^{3}
  scaleFactor: scale ther result by this scalar
  strainSens:  the sensitivity of each strain component
  vars:        the element variables
  Xpts:        the element nodal locations
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::addStrainSVSens(TacsScalar sens[],
                                               const double pt[],
                                               const TacsScalar scale,
                                               const TacsScalar strainSens[],
                                               const TacsScalar Xpts[],
                                               const TacsScalar vars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

  // Compute the element shape functions
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of X with respect to the coordinate
  // directions
  TacsScalar X[3], Xa[9];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

  // Compute the determinant of Xa and the transformation
  TacsScalar J[9];
  FElibrary::jacobian3d(Xa, J);

  // Get the derivative of the strain with respect to the nodal
  // displacements
  getBmat(B, J, Na, Nb, Nc, vars);

  TacsScalar *b = B;
  for (int i = 0; i < NUM_VARIABLES; i++) {
    sens[i] += scale * (strainSens[0] * b[0] + strainSens[1] * b[1] +
                        strainSens[2] * b[2] + strainSens[3] * b[3] +
                        strainSens[4] * b[4] + strainSens[5] * b[5]);
    b += NUM_STRESSES;
  }
}

/*
  Compute the strain and the derivative of the strain with respect to
  the nodal locations.

  output:
  strain:        the strain evaluate at the pamametric point pt
  strainXptSens: the derivative of the straint w.r.t. the nodal locations

  input:
  pt:        the parametric point within the element
  vars:      the element variables
  Xpts:      the nodal locations
*/
template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::addStrainXptSens(TacsScalar strainXptSens[],
                                                const double pt[],
                                                const TacsScalar scale,
                                                const TacsScalar strainSens[],
                                                const TacsScalar Xpts[],
                                                const TacsScalar vars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];

  // Compute the element shape functions
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of X with respect to the coordinate
  // directions
  TacsScalar X[3], Xa[9];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

  // Compute the determinant of Xa and the transformation
  TacsScalar J[9];
  FElibrary::jacobian3d(Xa, J);

  // Compute the derivative of the strain w.r.t. nocal coordinates
  addStrainXptSens(strainXptSens, scale, strainSens, J, Xa, Na, Nb, Nc, vars);
}

#endif  // TACS_3D_ELEMENT_H
