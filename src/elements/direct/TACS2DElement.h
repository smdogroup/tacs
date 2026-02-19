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

#ifndef TACS_2D_ELEMENT_H
#define TACS_2D_ELEMENT_H

/*
  The following file contains the general definition of a
  three-dimensional element that can be used in TACS.
*/

#include "FElibrary.h"
#include "PlaneStressStiffness.h"
#include "TACSElement.h"

/*
  The following class defines a generic two-dimensional element
  without defining the shape functions or quadrature scheme. This
  class can be used to implement multiple 2D elements.  This does not
  significantly impact the computational performance since the cost of
  the element computations is consumed in the inner product of the
  B-matrix with the constitutive matrix.
*/
template <int NUM_NODES>
class TACS2DElement : public TACSElement {
 public:
  // Define some constants for this element type
  static const int NUM_DISPS = 2;
  static const int NUM_STRESSES = 3;
  static const int NUM_EXTRAS = 4;
  static const int NUM_VARIABLES = 2 * NUM_NODES;

  TACS2DElement(PlaneStressStiffness *_stiff, ElementBehaviorType type,
                int _component);
  ~TACS2DElement();

  // Retrieve the shape functions
  // ----------------------------
  virtual void getShapeFunctions(const double pt[], double N[], double Na[],
                                 double Nb[]) = 0;

  // Compute the position vector and its derivative
  // ----------------------------------------------
  void planeJacobian(TacsScalar X[], TacsScalar Xa[], const double N[],
                     const double Na[], const double Nb[],
                     const TacsScalar Xpts[]);

  // Compute the displacement gradient
  // ---------------------------------
  void getDisplacement(TacsScalar U[], const double N[],
                       const TacsScalar vars[]);
  void getDisplGradient(TacsScalar Ud[], const TacsScalar J[],
                        const double Na[], const double Nb[],
                        const TacsScalar vars[]);
  void getDisplGradientSens(TacsScalar Ud[], TacsScalar UdSens[],
                            const TacsScalar J[], const TacsScalar JSens[],
                            const double Na[], const double Nb[],
                            const TacsScalar vars[]);

  // Evaluate the strain in the element
  // ---------------------------------
  void evalStrain(TacsScalar strain[], const TacsScalar J[], const double Na[],
                  const double Nb[], const TacsScalar vars[]);

  // Compute the derivative of the strain with respect to vars
  // ---------------------------------------------------------
  void getBmat(TacsScalar B[], const TacsScalar J[], const double Na[],
               const double Nb[], const TacsScalar vars[]);

  // Add the second derivatives of the strain times the stress to the
  // upper portion of the matrix
  // ----------------------------------------------------------------
  void addGeoStiffness(TacsScalar kmat[], TacsScalar h,
                       const TacsScalar stress[], const TacsScalar J[],
                       const double Na[], const double Nb[]);

  // Compute the derivative of the strain with respect to the nodal coordinates
  // --------------------------------------------------------------------------
  void addStrainXptSens(TacsScalar sens[], TacsScalar scale,
                        const TacsScalar strainSens[], const TacsScalar J[],
                        const TacsScalar Xa[], const double Na[],
                        const double Nb[], const TacsScalar vars[]);

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

  // Retrieve a specific time-independent matrix from the element
  // ------------------------------------------------------------
  void getMatType(ElementMatrixType matType, TacsScalar mat[],
                  const TacsScalar Xpts[], const TacsScalar vars[]);

  // Functions for evaluating global functionals of interest
  // -------------------------------------------------------
  TACSConstitutive *getConstitutive() { return stiff; }

  // Evaluate the determinant of the Jacobian and its derivative
  // -----------------------------------------------------------
  TacsScalar getDetJacobian(const double *pt, const TacsScalar Xpts[]);
  TacsScalar getDetJacobianXptSens(TacsScalar *sh, const double *pt,
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
  PlaneStressStiffness *stiff;

 private:
  static const char *dispNames[NUM_DISPS];
  static const char *stressNames[NUM_STRESSES];
  static const char *strainNames[NUM_STRESSES];
  static const char *extraNames[NUM_EXTRAS];
};

/*
  The constructor for the 2D element matrix
*/
template <int NUM_NODES>
TACS2DElement<NUM_NODES>::TACS2DElement(PlaneStressStiffness *_stiff,
                                        ElementBehaviorType type, int component)
    : TACSElement(component) {
  strain_type = type;
  stiff = _stiff;
  stiff->incref();
}

template <int NUM_NODES>
TACS2DElement<NUM_NODES>::~TACS2DElement() {
  stiff->decref();
}

/*
  Provide the names for the different components of the displacements
  and stresses
*/
template <int NUM_NODES>
const char *TACS2DElement<NUM_NODES>::dispNames[] = {"u", "v"};

template <int NUM_NODES>
const char *TACS2DElement<NUM_NODES>::stressNames[] = {"sxx", "syy", "sxy"};

template <int NUM_NODES>
const char *TACS2DElement<NUM_NODES>::strainNames[] = {"exx", "eyy", "exy"};

template <int NUM_NODES>
const char *TACS2DElement<NUM_NODES>::extraNames[] = {"lambda", "buckling",
                                                      "dv1", "dv2"};

/*
  Get the names of the displacements/stress etc.
*/
template <int NUM_NODES>
const char *TACS2DElement<NUM_NODES>::displacementName(int i) {
  if (i >= 0 && i < NUM_DISPS) {
    return dispNames[i];
  }
  return NULL;
}

template <int NUM_NODES>
const char *TACS2DElement<NUM_NODES>::stressName(int i) {
  if (i >= 0 && i < NUM_STRESSES) {
    return stressNames[i];
  }
  return NULL;
}

template <int NUM_NODES>
const char *TACS2DElement<NUM_NODES>::strainName(int i) {
  if (i >= 0 && i < NUM_STRESSES) {
    return strainNames[i];
  }
  return NULL;
}

template <int NUM_NODES>
const char *TACS2DElement<NUM_NODES>::extraName(int i) {
  if (i >= 0 && i < NUM_EXTRAS) {
    return extraNames[i];
  }
  return NULL;
}

/*
  Retrieve information about the number of displacements/stress etc.
*/
template <int NUM_NODES>
int TACS2DElement<NUM_NODES>::numDisplacements() {
  return NUM_DISPS;
}

template <int NUM_NODES>
int TACS2DElement<NUM_NODES>::numStresses() {
  return NUM_STRESSES;
}

template <int NUM_NODES>
int TACS2DElement<NUM_NODES>::numNodes() {
  return NUM_NODES;
}

template <int NUM_NODES>
int TACS2DElement<NUM_NODES>::numVariables() {
  return NUM_VARIABLES;
}

template <int NUM_NODES>
int TACS2DElement<NUM_NODES>::numExtras() {
  return NUM_EXTRAS;
}

template <int NUM_NODES>
ElementType TACS2DElement<NUM_NODES>::getElementType() {
  return TACS_PLANE_STRESS;
}

/*
  Set the design variable values
*/
template <int NUM_NODES>
void TACS2DElement<NUM_NODES>::setDesignVars(const TacsScalar dvs[],
                                             int numDVs) {
  stiff->setDesignVars(dvs, numDVs);
}

/*
  Retrive the design variable numbers
*/
template <int NUM_NODES>
void TACS2DElement<NUM_NODES>::getDesignVars(TacsScalar dvs[], int numDVs) {
  stiff->getDesignVars(dvs, numDVs);
}

/*
  Set the design variable lower/upper bounds in the provided arrays
*/
template <int NUM_NODES>
void TACS2DElement<NUM_NODES>::getDesignVarRange(TacsScalar lowerBound[],
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
  Na, Nb:      the derivative of the shape functions
  Xpts:        the nodal locations
*/
template <int NUM_NODES>
void TACS2DElement<NUM_NODES>::planeJacobian(TacsScalar X[], TacsScalar Xa[],
                                             const double N[],
                                             const double Na[],
                                             const double Nb[],
                                             const TacsScalar Xpts[]) {
  X[0] = X[1] = 0.0;
  Xa[0] = Xa[1] = Xa[2] = Xa[3] = 0.0;

  for (int i = 0; i < NUM_NODES; i++) {
    X[0] += Xpts[0] * N[0];
    X[1] += Xpts[1] * N[0];

    Xa[0] += Xpts[0] * Na[0];
    Xa[1] += Xpts[0] * Nb[0];

    Xa[2] += Xpts[1] * Na[0];
    Xa[3] += Xpts[1] * Nb[0];

    N++;
    Na++;
    Nb++;
    Xpts += 3;
  }
}

/*
  Compute the displacement given the provided values of the shape functions
*/
template <int NUM_NODES>
void TACS2DElement<NUM_NODES>::getDisplacement(TacsScalar U[], const double N[],
                                               const TacsScalar vars[]) {
  // zero the displacement values
  U[0] = U[1] = 0.0;

  for (int i = 0; i < NUM_NODES; i++) {
    U[0] += N[0] * vars[0];
    U[1] += N[0] * vars[1];

    N++;
    vars += 2;
  }
}

/*
  Compute the displacement gradient using the provided basis functions.

  The displacement gradient is computed as follows:
  Ud = Ua*J =
  [ Ua[0]  Ua[1] ][ J[0]  J[1] ]
  [ Ua[2]  Ua[3] ][ J[2]  J[3] ]

  output:
  Ud:    the displacement gradient

  input:
  J:          the transformation J = [X,a]^{-1}
  Na, Nb:     the derivatives of the shape functions
  vars:       the element variables
*/
template <int NUM_NODES>
void TACS2DElement<NUM_NODES>::getDisplGradient(TacsScalar Ud[],
                                                const TacsScalar J[],
                                                const double Na[],
                                                const double Nb[],
                                                const TacsScalar vars[]) {
  TacsScalar Ua[4];
  Ua[0] = Ua[1] = Ua[2] = Ua[3] = 0.0;

  // Compute the derivative of the u,v displacements with
  // respect to the parametric locations within the element
  for (int i = 0; i < NUM_NODES; i++) {
    Ua[0] += vars[0] * Na[0];
    Ua[1] += vars[0] * Nb[0];

    Ua[2] += vars[1] * Na[0];
    Ua[3] += vars[1] * Nb[0];

    Na++;
    Nb++;
    vars += 2;
  }
  // Compute the displacement gradient: Ud = Ua*J
  Ud[0] = Ua[0] * J[0] + Ua[1] * J[2];
  Ud[2] = Ua[2] * J[0] + Ua[3] * J[2];

  Ud[1] = Ua[0] * J[1] + Ua[1] * J[3];
  Ud[3] = Ua[2] * J[1] + Ua[3] * J[3];
}

/*
  Compute the displacement gradient and the derivative of the
  displacement gradient using the provided basis functions.

  The displacement gradient is computed as follows:
  Ud = Ua*J =
  [ Ua[0]  Ua[1] ][ J[0]  J[1] ]
  [ Ua[2]  Ua[3] ][ J[2]  J[3] ]

  output:
  Ud:       the displacement gradient
  UdSens:   the derivative of the displacement gradient

  input:
  J:          the transformation J = [X,a]^{-1}
  JSens:      the derivative of the inverse Jacobian
  Na, Nb:     the derivatives of the shape functions
  vars:       the element variables
*/
template <int NUM_NODES>
void TACS2DElement<NUM_NODES>::getDisplGradientSens(
    TacsScalar Ud[], TacsScalar UdSens[], const TacsScalar J[],
    const TacsScalar JSens[], const double Na[], const double Nb[],
    const TacsScalar vars[]) {
  TacsScalar Ua[4];
  Ua[0] = Ua[1] = Ua[2] = Ua[3] = 0.0;

  // Compute the derivative of the u,v displacements with
  // respect to the parametric locations within the element
  for (int i = 0; i < NUM_NODES; i++) {
    Ua[0] += vars[0] * Na[0];
    Ua[1] += vars[0] * Nb[0];

    Ua[2] += vars[1] * Na[0];
    Ua[3] += vars[1] * Nb[0];

    Na++;
    Nb++;
    vars += 2;
  }

  // Compute the displacement gradient: Ud = Ua*J
  Ud[0] = Ua[0] * J[0] + Ua[1] * J[2];  // u,x
  Ud[2] = Ua[2] * J[0] + Ua[3] * J[2];  // v,x

  Ud[1] = Ua[0] * J[1] + Ua[1] * J[3];  // u,y
  Ud[3] = Ua[2] * J[1] + Ua[3] * J[3];  // v,y

  // Compute the derivative of the displacement gradient
  UdSens[0] = Ua[0] * JSens[0] + Ua[1] * JSens[2];
  UdSens[2] = Ua[2] * JSens[0] + Ua[3] * JSens[2];

  UdSens[1] = Ua[0] * JSens[1] + Ua[1] * JSens[3];
  UdSens[3] = Ua[2] * JSens[1] + Ua[3] * JSens[3];
}

/*
  Compute the strain using the specified transformation, shape
  functions and variables

  output:
  strain:   the strain

  input:
  J:          the Jacobian of the transformation
  Na, Nb:     the derivatives of the basis functions
  vars:       the variables
*/
template <int NUM_NODES>
void TACS2DElement<NUM_NODES>::evalStrain(TacsScalar strain[],
                                          const TacsScalar J[],
                                          const double Na[], const double Nb[],
                                          const TacsScalar vars[]) {
  // Compute the displacement gradient
  TacsScalar Ud[4];
  getDisplGradient(Ud, J, Na, Nb, vars);

  // Compute the strain using either linear or nonlinear expression
  if (strain_type == LINEAR) {
    strain[0] = Ud[0];
    strain[1] = Ud[3];
    strain[2] = Ud[1] + Ud[2];
  } else {
    strain[0] = Ud[0] + 0.5 * (Ud[0] * Ud[0] + Ud[2] * Ud[2]);
    strain[1] = Ud[3] + 0.5 * (Ud[1] * Ud[1] + Ud[3] * Ud[3]);
    strain[2] = Ud[1] + Ud[2] + (Ud[0] * Ud[1] + Ud[2] * Ud[3]);
  }
}

/*
  Compute the derivative of the strain with respect to vars
*/
template <int NUM_NODES>
void TACS2DElement<NUM_NODES>::getBmat(TacsScalar B[], const TacsScalar J[],
                                       const double Na[], const double Nb[],
                                       const TacsScalar vars[]) {
  if (strain_type == LINEAR) {
    // If this is a linear element, then things are relative easy to
    // deal with - we just compute B alternatively by row
    for (int i = 0; i < NUM_NODES; i++) {
      TacsScalar Dx = Na[0] * J[0] + Nb[0] * J[2];
      TacsScalar Dy = Na[0] * J[1] + Nb[0] * J[3];

      B[0] = Dx;
      B[1] = 0.0;
      B[2] = Dy;
      B += 3;

      B[0] = 0.0;
      B[1] = Dy;
      B[2] = Dx;
      B += 3;

      Na++;
      Nb++;
    }
  } else {
    // Compute the displacement gradient: Ud = Ua*J
    TacsScalar Ud[4];
    getDisplGradient(Ud, J, Na, Nb, vars);

    // Compute the derivative of the strain with respect to
    // the nodal displacements
    for (int i = 0; i < NUM_NODES; i++) {
      TacsScalar Dx = Na[0] * J[0] + Nb[0] * J[2];
      TacsScalar Dy = Na[0] * J[1] + Nb[0] * J[3];

      B[0] = Dx + Ud[0] * Dx;
      B[1] = Ud[1] * Dy;
      B[2] = Dy + Dx * Ud[1] + Ud[0] * Dy;
      B += 3;

      B[0] = Ud[2] * Dx;
      B[1] = Dy + Ud[3] * Dy;
      B[2] = Dx + Dx * Ud[3] + Ud[2] * Dy;
      B += 3;

      Na++;
      Nb++;
    }
  }
}

/*
  Add the second derivatives of the strain times the stress to the
  upper portion of the matrix
*/
template <int NUM_NODES>
void TACS2DElement<NUM_NODES>::addGeoStiffness(TacsScalar mat[], TacsScalar h,
                                               const TacsScalar stress[],
                                               const TacsScalar J[],
                                               const double Na[],
                                               const double Nb[]) {
  if (!(strain_type == LINEAR)) {
    for (int j = 0; j < NUM_NODES; j++) {
      TacsScalar Dxj = Na[j] * J[0] + Nb[j] * J[2];
      TacsScalar Dyj = Na[j] * J[1] + Nb[j] * J[3];

      for (int i = 0; i <= j; i++) {
        TacsScalar Dxi = Na[i] * J[0] + Nb[i] * J[2];
        TacsScalar Dyi = Na[i] * J[1] + Nb[i] * J[3];

        // Add the contributions to the stiffness matrix
        TacsScalar scale = h * (stress[0] * Dxi * Dxj + stress[1] * Dyi * Dyj +
                                stress[2] * (Dxi * Dyj + Dxj * Dyi));
        mat[2 * i + 2 * j * NUM_VARIABLES] += scale;
        mat[2 * i + 1 + (2 * j + 1) * NUM_VARIABLES] += scale;
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
  [ J[0]  J[1] ][ Na[0]  Nb[0] ][ J[0]  J[1] ]
  [ J[2]  J[3] ][     0      0 ][ J[2]  J[3] ]

  (where the middle term is the derivative d(Xa)/dx). These
  derivatives take a special form as follows:

  For the x-derivative:
  dJ/dXa =
  [ J[0]  J[1] ][ d1  d2 ]
  [ J[2]  J[3] ][  0   0 ]
  =
  [ d1*J[0]  d2*J[0] ]
  [ d1*J[2]  d2*J[2] ]

  where:
  d1 = Na[0]*J[0] + Nb[0]*J[2]
  d2 = Na[0]*J[1] + Nb[0]*J[3]
*/
template <int NUM_NODES>
void TACS2DElement<NUM_NODES>::addStrainXptSens(
    TacsScalar sens[], TacsScalar scale, const TacsScalar strainSens[],
    const TacsScalar J[], const TacsScalar Xa[], const double Na[],
    const double Nb[], const TacsScalar vars[]) {
  // Compute the derivative of the u,v,w displacements with
  // respect to the parametric locations within the element
  TacsScalar Ua[4];
  Ua[0] = Ua[1] = Ua[2] = Ua[3] = 0.0;

  const double *na = Na, *nb = Nb;
  for (int i = 0; i < NUM_NODES; i++) {
    Ua[0] += vars[0] * na[0];
    Ua[1] += vars[0] * nb[0];

    Ua[2] += vars[1] * na[0];
    Ua[3] += vars[1] * nb[0];

    na++;
    nb++;
    vars += 2;
  }

  // Compute the scaled strain sensitivity
  TacsScalar eSens[3];
  eSens[0] = scale * strainSens[0];
  eSens[1] = scale * strainSens[1];
  eSens[2] = scale * strainSens[2];

  if (strain_type == LINEAR) {
    for (int i = 0; i < NUM_NODES; i++) {
      // JSens = -J*d(Xa)/dx*J
      TacsScalar d1 = -(Na[0] * J[0] + Nb[0] * J[2]);
      TacsScalar d2 = -(Na[0] * J[1] + Nb[0] * J[3]);

      // The derivative of the inverse transformation matrix
      // in each coordinate direction
      TacsScalar Ja[4], Jb[4];
      Ja[0] = d1 * J[0];
      Ja[1] = d2 * J[0];
      Ja[2] = d1 * J[2];
      Ja[3] = d2 * J[2];

      Jb[0] = d1 * J[1];
      Jb[1] = d2 * J[1];
      Jb[2] = d1 * J[3];
      Jb[3] = d2 * J[3];

      // Compute the derivative of the displacement gradient
      // with respect to each of the three coordinate directions
      // displacement gradient: Ud = Ua*J
      TacsScalar Uda[4];
      Uda[0] = Ua[0] * Ja[0] + Ua[1] * Ja[2];
      Uda[2] = Ua[2] * Ja[0] + Ua[3] * Ja[2];

      Uda[1] = Ua[0] * Ja[1] + Ua[1] * Ja[3];
      Uda[3] = Ua[2] * Ja[1] + Ua[3] * Ja[3];

      TacsScalar Udb[4];
      Udb[0] = Ua[0] * Jb[0] + Ua[1] * Jb[2];
      Udb[2] = Ua[2] * Jb[0] + Ua[3] * Jb[2];

      Udb[1] = Ua[0] * Jb[1] + Ua[1] * Jb[3];
      Udb[3] = Ua[2] * Jb[1] + Ua[3] * Jb[3];

      sens[0] += (Uda[0] * eSens[0] + Uda[3] * eSens[1] +
                  (Uda[1] + Uda[2]) * eSens[2]);
      sens++;

      sens[0] += (Udb[0] * eSens[0] + Udb[3] * eSens[1] +
                  (Udb[1] + Udb[2]) * eSens[2]);
      sens += 2;

      Na++;
      Nb++;
    }
  } else {
    // Compute the displacement gradient: Ud = Ua*J
    TacsScalar Ud[4];
    Ud[0] = Ua[0] * J[0] + Ua[1] * J[2];
    Ud[2] = Ua[2] * J[0] + Ua[3] * J[2];

    Ud[1] = Ua[0] * J[1] + Ua[1] * J[3];
    Ud[3] = Ua[2] * J[1] + Ua[3] * J[3];

    for (int i = 0; i < NUM_NODES; i++) {
      // JSens = -J*d(Xa)/dx*J
      TacsScalar d1 = -(Na[0] * J[0] + Nb[0] * J[2]);
      TacsScalar d2 = -(Na[0] * J[1] + Nb[0] * J[3]);

      // The derivative of the inverse transformation matrix
      // in each coordinate direction
      TacsScalar Ja[4], Jb[4];
      Ja[0] = d1 * J[0];
      Ja[1] = d2 * J[0];
      Ja[2] = d1 * J[2];
      Ja[3] = d2 * J[2];

      Jb[0] = d1 * J[1];
      Jb[1] = d2 * J[1];
      Jb[2] = d1 * J[3];
      Jb[3] = d2 * J[3];

      // Compute the derivative of the displacement gradient
      // with respect to each of the three coordinate directions
      // displacement gradient: Ud = Ua*J
      TacsScalar Uda[4];
      Uda[0] = Ua[0] * Ja[0] + Ua[1] * Ja[2];
      Uda[2] = Ua[2] * Ja[0] + Ua[3] * Ja[2];

      Uda[1] = Ua[0] * Ja[1] + Ua[1] * Ja[3];
      Uda[3] = Ua[2] * Ja[1] + Ua[3] * Ja[3];

      TacsScalar Udb[4];
      Udb[0] = Ua[0] * Jb[0] + Ua[1] * Jb[2];
      Udb[2] = Ua[2] * Jb[0] + Ua[3] * Jb[2];

      Udb[1] = Ua[0] * Jb[1] + Ua[1] * Jb[3];
      Udb[3] = Ua[2] * Jb[1] + Ua[3] * Jb[3];

      // Compute the derivative of the strain
      sens[0] += ((Uda[0] + (Ud[0] * Uda[0] + Ud[2] * Uda[2])) * eSens[0] +
                  (Uda[3] + (Ud[1] * Uda[1] + Ud[3] * Uda[3])) * eSens[1] +
                  (Uda[1] + Uda[2] +
                   (Uda[0] * Ud[1] + Uda[2] * Ud[3] + Ud[0] * Uda[1] +
                    Ud[2] * Uda[3])) *
                      eSens[2]);
      sens++;

      sens[0] += ((Udb[0] + (Ud[0] * Udb[0] + Ud[2] * Udb[2])) * eSens[0] +
                  (Udb[3] + (Ud[1] * Udb[1] + Ud[3] * Udb[3])) * eSens[1] +
                  (Udb[1] + Udb[2] +
                   (Udb[0] * Ud[1] + Udb[2] * Ud[3] + Ud[0] * Udb[1] +
                    Ud[2] * Udb[3])) *
                      eSens[2]);
      sens += 2;

      Na++;
      Nb++;
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
void TACS2DElement<NUM_NODES>::computeEnergies(double time, TacsScalar *_Te,
                                               TacsScalar *_Pe,
                                               const TacsScalar Xpts[],
                                               const TacsScalar vars[],
                                               const TacsScalar dvars[]) {
  // Compute the kinetic and potential energy
  TacsScalar Te = 0.0, Pe = 0.0;

  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES];

  // Get the number of quadrature points
  int numGauss = getNumGaussPts();

  for (int n = 0; n < numGauss; n++) {
    // Retrieve the quadrature points and weight
    double pt[3];
    double weight = getGaussWtsPts(n, pt);

    // Compute the element shape functions
    getShapeFunctions(pt, N, Na, Nb);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[9];
    planeJacobian(X, Xa, N, Na, Nb, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[4];
    TacsScalar h = FElibrary::jacobian2d(Xa, J);
    h = h * weight;

    // Compute the strain
    TacsScalar strain[NUM_STRESSES];
    evalStrain(strain, J, Na, Nb, vars);

    // Compute the corresponding stress
    TacsScalar stress[NUM_STRESSES];
    stiff->calculateStress(pt, strain, stress);

    // Compute the contribution from the potential energy
    Pe +=
        0.5 * h *
        (stress[0] * strain[0] + stress[1] * strain[1] + stress[2] * strain[2]);

    // Get value of the mass/area at this point
    TacsScalar mass;
    stiff->getPointwiseMass(pt, &mass);

    // Compute the contribution from the kinetic energy
    TacsScalar dUdt[2];
    getDisplacement(dUdt, N, dvars);
    Te += 0.5 * h * mass * (dUdt[0] * dUdt[0] + dUdt[1] * dUdt[1]);
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
void TACS2DElement<NUM_NODES>::addResidual(double time, TacsScalar *res,
                                           const TacsScalar Xpts[],
                                           const TacsScalar vars[],
                                           const TacsScalar dvars[],
                                           const TacsScalar ddvars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

  // Get the number of quadrature points
  int numGauss = getNumGaussPts();

  for (int n = 0; n < numGauss; n++) {
    // Retrieve the quadrature points and weight
    double pt[3];
    double weight = getGaussWtsPts(n, pt);

    // Compute the element shape functions
    getShapeFunctions(pt, N, Na, Nb);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[9];
    planeJacobian(X, Xa, N, Na, Nb, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[4];
    TacsScalar h = FElibrary::jacobian2d(Xa, J);
    h = h * weight;

    // Compute the strain
    TacsScalar strain[NUM_STRESSES];
    evalStrain(strain, J, Na, Nb, vars);

    // Compute the corresponding stress
    TacsScalar stress[NUM_STRESSES];
    stiff->calculateStress(pt, strain, stress);

    // Get the derivative of the strain with respect to the nodal
    // displacements
    getBmat(B, J, Na, Nb, vars);

    TacsScalar *b = B;
    for (int i = 0; i < NUM_VARIABLES; i++) {
      res[i] += h * (b[0] * stress[0] + b[1] * stress[1] + b[2] * stress[2]);
      b += NUM_STRESSES;
    }

    // Get value of the mass/area at this point
    TacsScalar mass;
    stiff->getPointwiseMass(pt, &mass);

    // Add the contribution from the inertial terms
    TacsScalar d2Udt2[2];
    getDisplacement(d2Udt2, N, ddvars);
    for (int i = 0; i < NUM_NODES; i++) {
      res[2 * i] += h * mass * N[i] * d2Udt2[0];
      res[2 * i + 1] += h * mass * N[i] * d2Udt2[1];
    }
  }
}

/*
  Add the Jacobian of the governing equations - the exact Jacobian of
  the residual expression.

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
void TACS2DElement<NUM_NODES>::addJacobian(
    double time, TacsScalar mat[], double alpha, double beta, double gamma,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

  // Get the number of quadrature points
  int numGauss = getNumGaussPts();

  for (int n = 0; n < numGauss; n++) {
    // Retrieve the quadrature points and weight
    double pt[3];
    double weight = getGaussWtsPts(n, pt);

    // Compute the element shape functions
    getShapeFunctions(pt, N, Na, Nb);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[4];
    planeJacobian(X, Xa, N, Na, Nb, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[4];
    TacsScalar h = FElibrary::jacobian2d(Xa, J);
    h = h * weight;

    if (alpha != 0.0) {
      // Compute the strain
      TacsScalar strain[NUM_STRESSES];
      evalStrain(strain, J, Na, Nb, vars);

      // Compute the corresponding stress
      TacsScalar stress[NUM_STRESSES];
      stiff->calculateStress(pt, strain, stress);

      // Add the stress times the second derivative of the strain
      addGeoStiffness(mat, alpha * h, stress, J, Na, Nb);

      // Get the derivative of the strain with respect to the nodal
      // displacements
      getBmat(B, J, Na, Nb, vars);

      // Fill-in the upper portion of the matrix
      TacsScalar *bj = B;
      for (int j = 0; j < NUM_VARIABLES; j++) {
        // Compute the stress at the given point
        TacsScalar bs[NUM_STRESSES];
        stiff->calculateStress(pt, bj, bs);

        TacsScalar *bi = B;
        for (int i = 0; i <= j; i++) {
          mat[i + j * NUM_VARIABLES] +=
              alpha * h * (bi[0] * bs[0] + bi[1] * bs[1] + bi[2] * bs[2]);
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
          mat[2 * i + 2 * j * NUM_VARIABLES] += scale * N[i] * N[j];
          mat[2 * i + 1 + (2 * j + 1) * NUM_VARIABLES] += scale * N[i] * N[j];
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
  mat:     the element tangent stiffness matrix
  res:     the element residual

  input:
  psi:     the element adjoint variables
  vars:    the element variables
  Xpts:    the element nodal locations in R^{3}
*/
template <int NUM_NODES>
void TACS2DElement<NUM_NODES>::addAdjResProduct(
    double time, double scale, TacsScalar dvSens[], int dvLen,
    const TacsScalar psi[], const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

  // Get the number of quadrature points
  int numGauss = getNumGaussPts();

  for (int n = 0; n < numGauss; n++) {
    // Retrieve the quadrature points and weights
    double pt[3];
    double weight = getGaussWtsPts(n, pt);

    // Compute the element shape functions
    getShapeFunctions(pt, N, Na, Nb);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[4];
    planeJacobian(X, Xa, N, Na, Nb, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[4];
    TacsScalar h = FElibrary::jacobian2d(Xa, J);
    h = h * weight;

    // Compute the strain
    TacsScalar strain[NUM_STRESSES];
    evalStrain(strain, J, Na, Nb, vars);

    // Get the derivative of the strain with respect to the nodal
    // displacements
    getBmat(B, J, Na, Nb, vars);

    // Compute the product of psi^{T}*B^{T}
    TacsScalar bpsi[NUM_STRESSES];
    memset(bpsi, 0, NUM_STRESSES * sizeof(TacsScalar));

    TacsScalar *b = B;
    const TacsScalar *ps = psi;
    for (int i = 0; i < NUM_VARIABLES; i++) {
      bpsi[0] += ps[0] * b[0];
      bpsi[1] += ps[0] * b[1];
      bpsi[2] += ps[0] * b[2];
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
  Add the product of the adjoint vector times the derivative of the
  residuals multiplied by a scalar to the given derivative vector.

  output:
  fXptSens:  derivative of the product of psi^{T}*R w.r.t. node locations

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
void TACS2DElement<NUM_NODES>::addAdjResXptProduct(
    double time, double scale, TacsScalar fXptSens[], const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[]) {
  /*
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  // Get the number of quadrature points
  int numGauss = getNumGaussPts();

  for ( int n = 0; n < numGauss; n++ ){
    // Retrieve the quadrature points and weights
    double pt[3];
    double weight = getGaussWtsPts(n, pt);

    // Compute the element shape functions
    getShapeFunctions(pt, N, Na, Nb);

    // Compute the derivative of X with respect to the
    // coordinate directions
    TacsScalar X[3], Xa[4];
    planeJacobian(X, Xa, N, Na, Nb, Xpts);

    // Compute the determinant of Xa and the transformation
    TacsScalar J[4];
    TacsScalar h = FElibrary::jacobian2d(Xa, J);
    h = h*weight;

    // Compute the derivative of the determinant w.r.t. node locations
    TacsScalar hXptSens[3*NUM_NODES];
    TacsScalar h = getDetJacobianXptSens(hXptSens, pt, XptSens);

    // Compute the strain
    TacsScalar strain[NUM_STRESSES];
    evalStrain(strain, J, Na, Nb, vars);

    // Get the derivative of the strain with respect to the nodal
    // displacements
    getBmat(B, J, Na, Nb, vars);
  }
  */
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
void TACS2DElement<NUM_NODES>::getMatType(ElementMatrixType matType,
                                          TacsScalar mat[],
                                          const TacsScalar Xpts[],
                                          const TacsScalar vars[]) {
  memset(mat, 0, NUM_VARIABLES * NUM_VARIABLES * sizeof(TacsScalar));

  // The mass matrix
  if (matType == MASS_MATRIX) {
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES];

    // Get the number of quadrature points
    int numGauss = getNumGaussPts();

    for (int n = 0; n < numGauss; n++) {
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getGaussWtsPts(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[4];
      planeJacobian(X, Xa, N, Na, Nb, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[4];
      TacsScalar h = FElibrary::jacobian2d(Xa, J);
      h = h * weight;

      // Get the pointwise mass
      TacsScalar ptmass[3];
      stiff->getPointwiseMass(pt, ptmass);

      // Fill-in the upper-portion of the matrix
      for (int j = 0; j < NUM_NODES; j++) {
        for (int i = 0; i <= j; i++) {
          TacsScalar d = h * ptmass[0] * N[i] * N[j];
          mat[2 * i + 2 * j * NUM_VARIABLES] += d;
          mat[2 * i + 1 + (2 * j + 1) * NUM_VARIABLES] += d;
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
    double Na[NUM_NODES], Nb[NUM_NODES];

    // Get the number of quadrature points
    int numGauss = getNumGaussPts();

    for (int n = 0; n < numGauss; n++) {
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getGaussWtsPts(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      planeJacobian(X, Xa, N, Na, Nb, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[4];
      TacsScalar h = FElibrary::jacobian2d(Xa, J);
      h = h * weight;

      // Compute the strain derived from load path
      TacsScalar strain[NUM_STRESSES];
      evalStrain(strain, J, Na, Nb, vars);

      // Compute the corresponding stress
      TacsScalar stress[NUM_STRESSES];
      stiff->calculateStress(pt, strain, stress);

      for (int j = 0; j < NUM_NODES; j++) {
        TacsScalar Dxj = Na[j] * J[0] + Nb[j] * J[2];
        TacsScalar Dyj = Na[j] * J[1] + Nb[j] * J[3];

        for (int i = 0; i < NUM_NODES; i++) {
          TacsScalar Dxi = Na[i] * J[0] + Nb[i] * J[2];
          TacsScalar Dyi = Na[i] * J[1] + Nb[i] * J[3];

          // Add the contributions to the stiffness matrix
          TacsScalar scale =
              h * (stress[0] * Dxi * Dxj + stress[1] * Dyi * Dyj +
                   stress[2] * (Dxi * Dyj + Dxj * Dyi));

          mat[2 * i + 2 * j * NUM_VARIABLES] += scale;
          mat[2 * i + 1 + (2 * j + 1) * NUM_VARIABLES] += scale;
        }
      }
    }
  } else if (matType == STIFFNESS_MATRIX) {
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES];

    // The derivative of the stress with respect to the strain
    TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

    // Get the number of quadrature points
    int numGauss = getNumGaussPts();

    for (int n = 0; n < numGauss; n++) {
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getGaussWtsPts(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      planeJacobian(X, Xa, N, Na, Nb, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[4];
      TacsScalar h = FElibrary::jacobian2d(Xa, J);
      h = h * weight;

      // Compute the strain
      TacsScalar strain[NUM_STRESSES];
      evalStrain(strain, J, Na, Nb, vars);

      // Compute the corresponding stress
      TacsScalar stress[NUM_STRESSES];
      stiff->calculateStress(pt, strain, stress);

      // Get the derivative of the strain with respect to the nodal
      // displacements
      getBmat(B, J, Na, Nb, vars);

      // Fill-in the upper portion of the matrix
      TacsScalar *bj = B;
      for (int j = 0; j < NUM_VARIABLES; j++) {
        // Compute the stress at the given point
        TacsScalar bs[NUM_STRESSES];
        stiff->calculateStress(pt, bj, bs);

        TacsScalar *bi = B;
        for (int i = 0; i <= j; i++) {
          mat[i + j * NUM_VARIABLES] +=
              h * (bi[0] * bs[0] + bi[1] * bs[1] + bi[2] * bs[2]);
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
void TACS2DElement<NUM_NODES>::addMatDVSensInnerProduct(
    ElementMatrixType matType, double scale, TacsScalar dvSens[], int dvLen,
    const TacsScalar psi[], const TacsScalar phi[], const TacsScalar Xpts[],
    const TacsScalar vars[]) {
  if (matType == STIFFNESS_MATRIX) {
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES];

    // The derivative of the stress with respect to the strain
    TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

    // Get the number of quadrature points
    int numGauss = getNumGaussPts();

    for (int n = 0; n < numGauss; n++) {
      // Retrieve the quadrature points and weights
      double pt[2];
      double weight = getGaussWtsPts(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      planeJacobian(X, Xa, N, Na, Nb, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[4];
      TacsScalar h = FElibrary::jacobian2d(Xa, J);
      h = h * weight;

      // Get the derivative of the strain with respect to the nodal
      // displacements
      getBmat(B, J, Na, Nb, vars);

      // Compute the product of psi^{T}*B^{T}
      TacsScalar bpsi[NUM_STRESSES], bphi[NUM_STRESSES];
      memset(bpsi, 0, NUM_STRESSES * sizeof(TacsScalar));
      memset(bphi, 0, NUM_STRESSES * sizeof(TacsScalar));

      TacsScalar *b = B;
      const TacsScalar *ps = psi, *ph = phi;
      for (int i = 0; i < NUM_VARIABLES; i++) {
        bpsi[0] += ps[0] * b[0];
        bpsi[1] += ps[0] * b[1];
        bpsi[2] += ps[0] * b[2];

        bphi[0] += ph[0] * b[0];
        bphi[1] += ph[0] * b[1];
        bphi[2] += ph[0] * b[2];

        b += NUM_STRESSES;
        ps++;
        ph++;
      }

      // Add the result to the design variable vector
      stiff->addStressDVSens(pt, bphi, scale * h, bpsi, dvSens, dvLen);
    }
  } else if (matType == GEOMETRIC_STIFFNESS_MATRIX) {
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES];

    // Get the number of quadrature points
    int numGauss = getNumGaussPts();

    for (int n = 0; n < numGauss; n++) {
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getGaussWtsPts(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      planeJacobian(X, Xa, N, Na, Nb, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[4];
      TacsScalar h = FElibrary::jacobian2d(Xa, J);
      h = h * weight;

      // Compute the strain derived from load path
      TacsScalar strain[NUM_STRESSES];
      evalStrain(strain, J, Na, Nb, vars);

      // Compute dN/dx*phi
      // Compute the product of psi^{T}*G^{T} and G*phi
      TacsScalar gpsi[4], gphi[4];
      memset(gpsi, 0, 4 * sizeof(TacsScalar));
      memset(gphi, 0, 4 * sizeof(TacsScalar));
      for (int j = 0; j < NUM_NODES; j++) {
        TacsScalar Dx = Na[j] * J[0] + Nb[j] * J[2];
        TacsScalar Dy = Na[j] * J[1] + Nb[j] * J[3];

        gpsi[0] += Dx * psi[2 * j];
        gphi[0] += Dx * phi[2 * j];
        gpsi[1] += Dx * psi[2 * j + 1];
        gphi[1] += Dx * phi[2 * j + 1];

        gpsi[2] += Dy * psi[2 * j];
        gphi[2] += Dy * phi[2 * j];
        gpsi[3] += Dy * psi[2 * j + 1];
        gphi[3] += Dy * phi[2 * j + 1];
      }

      TacsScalar sumN[NUM_STRESSES];
      memset(sumN, 0, NUM_STRESSES * sizeof(TacsScalar));

      TacsScalar *gs = gpsi, *gh = gphi;
      for (int j = 0; j < 2; j++) {
        sumN[j] = gs[0] * gh[0] + gs[1] * gh[1];
        gs += 2;
        gh += 2;
      }
      for (int j = 0; j < 2; j++) {
        sumN[2] += gpsi[j] * gphi[2 + j] + gpsi[2 + j] * gphi[j];
      }

      // Add the result to the design variable vector
      stiff->addStressDVSens(pt, sumN, scale * h, strain, dvSens, dvLen);
    }
  } else if (matType == MASS_MATRIX) {
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES];

    // Get the number of quadrature points
    int numGauss = getNumGaussPts();

    for (int n = 0; n < numGauss; n++) {
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getGaussWtsPts(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      planeJacobian(X, Xa, N, Na, Nb, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[4];
      TacsScalar h = FElibrary::jacobian2d(Xa, J);
      h = h * weight;

      // Compute the nodal accelerations at the quadrature point
      TacsScalar upsi[2], uphi[2];
      upsi[0] = upsi[1] = 0.0;
      uphi[0] = uphi[1] = 0.0;

      double *ns = N;
      const TacsScalar *ps = psi, *ph = phi;
      for (int i = 0; i < NUM_NODES; i++) {
        upsi[0] += ns[0] * ps[0];
        upsi[1] += ns[0] * ps[1];

        uphi[0] += ns[0] * ph[0];
        uphi[1] += ns[0] * ph[1];

        ps += 2;
        ph += 2;
        ns++;
      }

      // Add the result to the design variable vector
      TacsScalar rho_alpha =
          scale * h * (upsi[0] * uphi[0] + upsi[1] * uphi[1]);

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
void TACS2DElement<NUM_NODES>::getMatSVSensInnerProduct(
    ElementMatrixType matType, TacsScalar res[], const TacsScalar psi[],
    const TacsScalar phi[], const TacsScalar Xpts[], const TacsScalar vars[]) {
  if (matType == GEOMETRIC_STIFFNESS_MATRIX) {
    // The shape functions associated with the element
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES];

    // Get the number of quadrature points
    int numGauss = getNumGaussPts();

    for (int n = 0; n < numGauss; n++) {
      // Retrieve the quadrature points and weight
      double pt[3];
      double weight = getGaussWtsPts(n, pt);

      // Compute the element shape functions
      getShapeFunctions(pt, N, Na, Nb);

      // Compute the derivative of X with respect to the
      // coordinate directions
      TacsScalar X[3], Xa[9];
      planeJacobian(X, Xa, N, Na, Nb, Xpts);

      // Compute the determinant of Xa and the transformation
      TacsScalar J[4];
      TacsScalar h = FElibrary::jacobian2d(Xa, J);
      h = h * weight;

      // Compute the strain derived from load path
      TacsScalar strain[NUM_STRESSES];
      evalStrain(strain, J, Na, Nb, vars);

      // The derivative of the stress with respect to the strain
      TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

      // Get the derivative of the strain with respect to the nodal
      // displacements
      getBmat(B, J, Na, Nb, vars);

      // Compute dN/dx*phi
      // Compute the product of psi^{T}*G^{T} and G*phi
      TacsScalar gpsi[4], gphi[4];
      memset(gpsi, 0, 4 * sizeof(TacsScalar));
      memset(gphi, 0, 4 * sizeof(TacsScalar));
      for (int j = 0; j < NUM_NODES; j++) {
        TacsScalar Dx = Na[j] * J[0] + Nb[j] * J[2];
        TacsScalar Dy = Na[j] * J[1] + Nb[j] * J[3];

        gpsi[0] += Dx * psi[2 * j];
        gphi[0] += Dx * phi[2 * j];
        gpsi[1] += Dx * psi[2 * j + 1];
        gphi[1] += Dx * phi[2 * j + 1];

        gpsi[2] += Dy * psi[2 * j];
        gphi[2] += Dy * phi[2 * j];
        gpsi[3] += Dy * psi[2 * j + 1];
        gphi[3] += Dy * phi[2 * j + 1];
      }

      TacsScalar sumN[NUM_STRESSES];
      memset(sumN, 0, NUM_STRESSES * sizeof(TacsScalar));
      TacsScalar *gs = gpsi, *gh = gphi;
      for (int j = 0; j < 2; j++) {
        sumN[j] = gs[0] * gh[0] + gs[1] * gh[1];
        gs += 2;
        gh += 2;
      }
      for (int j = 0; j < 2; j++) {
        sumN[2] += gpsi[j] * gphi[2 + j] + gpsi[2 + j] * gphi[j];
      }

      // Get D*sumN
      TacsScalar s[NUM_STRESSES];
      stiff->calculateStress(pt, sumN, s);

      TacsScalar *b = B;
      for (int j = 0; j < NUM_VARIABLES; j++) {
        res[j] += h * (b[0] * s[0] + b[1] * s[1] + b[2] * s[2]);
        b += NUM_STRESSES;
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
TacsScalar TACS2DElement<NUM_NODES>::getDetJacobian(const double pt[],
                                                    const TacsScalar Xpts[]) {
  // Compute the element shape functions
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES];
  getShapeFunctions(pt, N, Na, Nb);

  // Compute the derivative of the shape functions w.r.t. the
  // parametric locations
  TacsScalar X[3], Xa[4];
  planeJacobian(X, Xa, N, Na, Nb, Xpts);

  return FElibrary::jacobian2d(Xa);
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
TacsScalar TACS2DElement<NUM_NODES>::getDetJacobianXptSens(
    TacsScalar *hXptSens, const double *pt, const TacsScalar Xpts[]) {
  // Compute the element shape functions
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES];
  getShapeFunctions(pt, N, Na, Nb);

  // Compute the derivative of the shape functions w.r.t. the
  // parametric locations
  TacsScalar X[3], Xa[4];
  planeJacobian(X, Xa, N, Na, Nb, Xpts);

  // Evaluate the determinant of the Jacobian
  TacsScalar J[4];
  TacsScalar h = FElibrary::jacobian2d(Xa, J);

  for (int i = 0; i < NUM_NODES; i++) {
    for (int k = 0; k < 2; k++) {
      TacsScalar XaSens[4];
      XaSens[0] = XaSens[1] = XaSens[2] = XaSens[3] = 0.0;
      XaSens[2 * k] = Na[i];
      XaSens[2 * k + 1] = Nb[i];

      FElibrary::jacobian2dSens(Xa, XaSens, &hXptSens[0]);
      hXptSens++;
    }

    // Zero the z-component of the derivative
    hXptSens[0] = 0.0;
    hXptSens++;
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
void TACS2DElement<NUM_NODES>::getStrain(TacsScalar strain[], const double pt[],
                                         const TacsScalar Xpts[],
                                         const TacsScalar vars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES];

  // Compute the element shape functions
  getShapeFunctions(pt, N, Na, Nb);

  // Compute the derivative of X with respect to the coordinate directions
  TacsScalar X[3], Xa[4];
  planeJacobian(X, Xa, N, Na, Nb, Xpts);

  // Compute the determinant of Xa and the transformation
  TacsScalar J[4];
  FElibrary::jacobian2d(Xa, J);

  // Compute the strain
  evalStrain(strain, J, Na, Nb, vars);
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
  strainSens:  the sensitivity of each straint component
  vars:        the element variables
  Xpts:        the element nodal locations
*/
template <int NUM_NODES>
void TACS2DElement<NUM_NODES>::addStrainSVSens(TacsScalar strainSVSens[],
                                               const double pt[],
                                               const TacsScalar scale,
                                               const TacsScalar strainSens[],
                                               const TacsScalar Xpts[],
                                               const TacsScalar vars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

  // Compute the element shape functions
  getShapeFunctions(pt, N, Na, Nb);

  // Compute the derivative of X with respect to the coordinate
  // directions
  TacsScalar X[3], Xa[4];
  planeJacobian(X, Xa, N, Na, Nb, Xpts);

  // Compute the determinant of Xa and the transformation
  TacsScalar J[4];
  FElibrary::jacobian2d(Xa, J);

  // Get the derivative of the strain with respect to the nodal
  // displacements
  getBmat(B, J, Na, Nb, vars);

  TacsScalar *b = B;
  for (int i = 0; i < NUM_VARIABLES; i++) {
    strainSVSens[i] += scale * (strainSens[0] * b[0] + strainSens[1] * b[1] +
                                strainSens[2] * b[2]);
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
void TACS2DElement<NUM_NODES>::addStrainXptSens(TacsScalar strainXptSens[],
                                                const double pt[],
                                                const TacsScalar scale,
                                                const TacsScalar strainSens[],
                                                const TacsScalar Xpts[],
                                                const TacsScalar vars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES];

  // Compute the element shape functions
  getShapeFunctions(pt, N, Na, Nb);

  // Compute the derivative of X with respect to the coordinate
  // directions
  TacsScalar X[3], Xa[4];
  planeJacobian(X, Xa, N, Na, Nb, Xpts);

  // Compute the determinant of Xa and the transformation
  TacsScalar J[4];
  FElibrary::jacobian2d(Xa, J);

  // Compute the derivative of the strain w.r.t. nocal coordinates
  addStrainXptSens(strainXptSens, scale, strainSens, J, Xa, Na, Nb, vars);
}

#endif
