#ifndef TACS_2D_COUPLED_THERMO_ELEMENT_H
#define TACS_2D_COUPLED_THERMO_ELEMENT_H

/*
  The following file contains the general definition of a
  two-dimensional element using a coupled heat transfer problem
  with thermoelastic analysis

  Copyright (c) 2017 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

#include "CoupledThermoPlaneStressStiffness.h"
#include "FElibrary.h"
#include "ThermoElements.h"

/*
  The following class defines a generic two-dimensional element
  using a coupled heat transfer problem with thermoelasticity
  without defining the shape functions or quadrature scheme. This
  class can be used to implement multiple 2D elements.  This does not
  significantly impact the computational performance since the cost of
  the element computations is consumed in the inner product of the
  B-matrix with the constitutive matrix.
*/
template <int NUM_NODES>
class TACS2DCoupledThermoElement : public ThermoQuad {
 public:
  // Define some constants for this element type
  static const int NUM_DISPS = 3;  // u, v, dT
  static const int NUM_STRESSES = 3;
  static const int NUM_EXTRAS = 4;
  static const int NUM_VARIABLES = 3 * NUM_NODES;  // u, v, dT

  TACS2DCoupledThermoElement(CoupledThermoPlaneStressStiffness *_stiff,
                             ElementBehaviorType type, int _component);
  ~TACS2DCoupledThermoElement();

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
  void getTemperature(TacsScalar T[], const double N[],
                      const TacsScalar vars[]);
  // Evaluate the strain in the element
  // ----------------------------------
  void evalStrain(TacsScalar strain[], const TacsScalar J[], const double Na[],
                  const double Nb[], const TacsScalar vars[]);
  void evalBT(TacsScalar strain[], const TacsScalar J[], const double Na[],
              const double Nb[], const TacsScalar vars[]);

  // Compute the derivative of the strain with respect to vars
  // ---------------------------------------------------------
  void getBmat(TacsScalar B[], const TacsScalar J[], const double Na[],
               const double Nb[], const TacsScalar vars[]);
  void getBmatTemp(TacsScalar B[], const TacsScalar J[], const double Na[],
                   const double Nb[], const TacsScalar vars[]);

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

  // Retrieve a specific time-independent matrix from the element
  // ------------------------------------------------------------
  void getMatType(ElementMatrixType matType, TacsScalar mat[],
                  const TacsScalar Xpts[], const TacsScalar vars[]);

  // Functions for evaluating global functionals of interest
  // -------------------------------------------------------
  CoupledThermoPlaneStressStiffness *getConstitutive() { return stiff; }

  // Evaluate the determinant of the Jacobian and its derivative
  // -----------------------------------------------------------
  TacsScalar getDetJacobian(const double *pt, const TacsScalar Xpts[]);
  TacsScalar getDetJacobianXptSens(TacsScalar *sh, const double *pt,
                                   const TacsScalar Xpts[]);

  // Compute the point-wise strain and its derivative
  // ------------------------------------------------
  void getStrain(TacsScalar strain[], const double pt[],
                 const TacsScalar Xpts[], const TacsScalar vars[]);
  void getEffStrain(TacsScalar strain[], const double pt[],
                    const TacsScalar Xpts[], const TacsScalar vars[]);
  void addStrainXptSens(TacsScalar strainXptSens[], const double pt[],
                        const TacsScalar scale, const TacsScalar strainSens[],
                        const TacsScalar Xpts[], const TacsScalar vars[]);

  void addStrainSVSens(TacsScalar strainSVSens[], const double pt[],
                       const TacsScalar scale, const TacsScalar strainSens[],
                       const TacsScalar Xpts[], const TacsScalar vars[]);
  void addEffStrainSVSens(TacsScalar strainSVSens[], const double pt[],
                          const TacsScalar scale, const TacsScalar strainSens[],
                          const TacsScalar Xpts[], const TacsScalar vars[],
                          int vars_j = 0);
  void getBT(TacsScalar strain[], const double pt[], const TacsScalar Xpts[],
             const TacsScalar vars[]);
  void addBTSVSens(TacsScalar strainSVSens[], const double pt[],
                   const TacsScalar scale, const TacsScalar strainSens[],
                   const TacsScalar Xpts[], const TacsScalar vars[]);

 protected:
  ElementBehaviorType strain_type;
  CoupledThermoPlaneStressStiffness *stiff;

 private:
  static const char *dispNames[NUM_DISPS];
  static const char *stressNames[NUM_STRESSES];
  static const char *strainNames[NUM_STRESSES];
  static const char *extraNames[NUM_EXTRAS];
  int conduction, convection, radiation;
};
/*
  The constructor for the 2D element matrix with thermal loading
*/
template <int NUM_NODES>
TACS2DCoupledThermoElement<NUM_NODES>::TACS2DCoupledThermoElement(
    CoupledThermoPlaneStressStiffness *_stiff, ElementBehaviorType type,
    int component)
    : ThermoQuad(component) {
  strain_type = type;
  stiff = _stiff;
  stiff->incref();
  conduction = 1;
  convection = 0;
  radiation = 0;
}

template <int NUM_NODES>
TACS2DCoupledThermoElement<NUM_NODES>::~TACS2DCoupledThermoElement() {
  stiff->decref();
}

/*
  Provide the names for the different components of the displacements
  and stresses
*/
template <int NUM_NODES>
const char *TACS2DCoupledThermoElement<NUM_NODES>::dispNames[] = {"u", "v",
                                                                  "dT"};

template <int NUM_NODES>
const char *TACS2DCoupledThermoElement<NUM_NODES>::stressNames[] = {
    "sxx", "syy", "sxy"};

template <int NUM_NODES>
const char *TACS2DCoupledThermoElement<NUM_NODES>::strainNames[] = {
    "exx", "eyy", "exy"};

template <int NUM_NODES>
const char *TACS2DCoupledThermoElement<NUM_NODES>::extraNames[] = {
    "lambda", "buckling", "dv1", "dv2"};

/*
  Get the names of the displacements/stress etc.
*/
template <int NUM_NODES>
const char *TACS2DCoupledThermoElement<NUM_NODES>::displacementName(int i) {
  if (i >= 0 && i < NUM_DISPS) {
    return dispNames[i];
  }
  return NULL;
}
template <int NUM_NODES>
const char *TACS2DCoupledThermoElement<NUM_NODES>::stressName(int i) {
  if (i >= 0 && i < NUM_STRESSES) {
    return stressNames[i];
  }
  return NULL;
}

template <int NUM_NODES>
const char *TACS2DCoupledThermoElement<NUM_NODES>::strainName(int i) {
  if (i >= 0 && i < NUM_STRESSES) {
    return strainNames[i];
  }
  return NULL;
}

template <int NUM_NODES>
const char *TACS2DCoupledThermoElement<NUM_NODES>::extraName(int i) {
  if (i >= 0 && i < NUM_EXTRAS) {
    return extraNames[i];
  }
  return NULL;
}

/*
  Retrieve information about the number of displacements/stress etc.
*/
template <int NUM_NODES>
int TACS2DCoupledThermoElement<NUM_NODES>::numDisplacements() {
  return NUM_DISPS;
}

template <int NUM_NODES>
int TACS2DCoupledThermoElement<NUM_NODES>::numStresses() {
  return NUM_STRESSES;
}

template <int NUM_NODES>
int TACS2DCoupledThermoElement<NUM_NODES>::numNodes() {
  return NUM_NODES;
}

template <int NUM_NODES>
int TACS2DCoupledThermoElement<NUM_NODES>::numVariables() {
  return NUM_VARIABLES;
}

template <int NUM_NODES>
int TACS2DCoupledThermoElement<NUM_NODES>::numExtras() {
  return NUM_EXTRAS;
}

template <int NUM_NODES>
ElementType TACS2DCoupledThermoElement<NUM_NODES>::getElementType() {
  return TACS_PLANE_STRESS;
}

/*
  Set the design variable values
*/
template <int NUM_NODES>
void TACS2DCoupledThermoElement<NUM_NODES>::setDesignVars(
    const TacsScalar dvs[], int numDVs) {
  stiff->setDesignVars(dvs, numDVs);
}
/*
  Retrive the design variable numbers
*/
template <int NUM_NODES>
void TACS2DCoupledThermoElement<NUM_NODES>::getDesignVars(TacsScalar dvs[],
                                                          int numDVs) {
  stiff->getDesignVars(dvs, numDVs);
}
/*
  Set the design variable lower/upper bounds in the provided arrays
*/
template <int NUM_NODES>
void TACS2DCoupledThermoElement<NUM_NODES>::getDesignVarRange(
    TacsScalar lowerBound[], TacsScalar upperBound[], int numDVs) {
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
void TACS2DCoupledThermoElement<NUM_NODES>::planeJacobian(
    TacsScalar X[], TacsScalar Xa[], const double N[], const double Na[],
    const double Nb[], const TacsScalar Xpts[]) {
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
void TACS2DCoupledThermoElement<NUM_NODES>::getDisplacement(
    TacsScalar U[], const double N[], const TacsScalar vars[]) {
  // Zero the displacement values
  U[0] = U[1] = 0.0;

  for (int i = 0; i < NUM_NODES; i++) {
    U[0] += N[0] * vars[0];
    U[1] += N[0] * vars[1];

    N++;
    vars += 3;
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
void TACS2DCoupledThermoElement<NUM_NODES>::getDisplGradient(
    TacsScalar Ud[], const TacsScalar J[], const double Na[], const double Nb[],
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
    vars += 3;
  }
  // Compute the displacement gradient: Ud = Ua*J
  Ud[0] = Ua[0] * J[0] + Ua[1] * J[2];
  Ud[2] = Ua[2] * J[0] + Ua[3] * J[2];

  Ud[1] = Ua[0] * J[1] + Ua[1] * J[3];
  Ud[3] = Ua[2] * J[1] + Ua[3] * J[3];
}
/*
  Compute the term B*T where B is the derivative of the shape functions wrt x,y
  and B =
  [ J[0]  J[2] ] [ Na[0] Na[1] Na[2] Na[3] ]
  [ J[1]  J[3] ] [ Nb[0] Nb[1] Nb[2] Nb[3] ]
  output:
  strain:     the [2x1] vector B*T

  input:
  J:          the transformation J = [X,a]^{-1}
  Na, Nb:     the derivatives of the shape functions
  vars:       the element variables
*/
template <int NUM_NODES>
void TACS2DCoupledThermoElement<NUM_NODES>::evalBT(TacsScalar strain[],
                                                   const TacsScalar J[],
                                                   const double Na[],
                                                   const double Nb[],
                                                   const TacsScalar vars[]) {
  TacsScalar Ua[2];
  Ua[0] = Ua[1] = 0.0;

  // Compute the derivative of the u,v displacements with
  // respect to the parametric locations within the element
  for (int i = 0; i < NUM_NODES; i++) {
    Ua[0] += vars[2] * Na[0];
    Ua[1] += vars[2] * Nb[0];

    Na++;
    Nb++;
    vars += 3;
  }
  strain[0] = Ua[0] * J[0] + Ua[1] * J[2];
  strain[1] = Ua[0] * J[1] + Ua[1] * J[3];
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
void TACS2DCoupledThermoElement<NUM_NODES>::getDisplGradientSens(
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
    vars += 3;
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
  Compute the temperature given the provided values of the shape functions
*/
template <int NUM_NODES>
void TACS2DCoupledThermoElement<NUM_NODES>::getTemperature(
    TacsScalar T[], const double N[], const TacsScalar vars[]) {
  // Zero the temperature values
  T[0] = 0.0;
  for (int i = 0; i < NUM_NODES; i++) {
    T[0] += N[0] * vars[2];

    N++;
    vars += 3;
  }
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
void TACS2DCoupledThermoElement<NUM_NODES>::evalStrain(
    TacsScalar strain[], const TacsScalar J[], const double Na[],
    const double Nb[], const TacsScalar vars[]) {
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
void TACS2DCoupledThermoElement<NUM_NODES>::getBmat(TacsScalar B[],
                                                    const TacsScalar J[],
                                                    const double Na[],
                                                    const double Nb[],
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
  Compute the derivative of the shape functions [2xNUM_NODES]
*/
template <int NUM_NODES>
void TACS2DCoupledThermoElement<NUM_NODES>::getBmatTemp(
    TacsScalar B[], const TacsScalar J[], const double Na[], const double Nb[],
    const TacsScalar vars[]) {
  if (strain_type == LINEAR) {
    for (int i = 0; i < NUM_NODES; i++) {
      TacsScalar Dx = Na[0] * J[0] + Nb[0] * J[2];
      TacsScalar Dy = Na[0] * J[1] + Nb[0] * J[3];

      B[0] = Dx;
      B[1] = Dy;
      B += 2;
      Na++;
      Nb++;
    }
  }
}

//-----------------------------------------------------------------------------------
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
void TACS2DCoupledThermoElement<NUM_NODES>::addStrainXptSens(
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
    vars += 3;
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
void TACS2DCoupledThermoElement<NUM_NODES>::computeEnergies(
    double time, TacsScalar *_Te, TacsScalar *_Pe, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[]) {
  // Compute the kinetic and potential energy
  TacsScalar Te = 0.0, Pe = 0.0;

  /* // The shape functions associated with the element */
  /* double N[NUM_NODES]; */
  /* double Na[NUM_NODES], Nb[NUM_NODES]; */

  /* // Get the number of quadrature points */
  /* int numGauss = getNumGaussPts(); */

  /* for ( int n = 0; n < numGauss; n++ ){ */
  /*   // Retrieve the quadrature points and weight */
  /*   double pt[3]; */
  /*   double weight = getGaussWtsPts(n, pt); */

  /*   // Compute the element shape functions */
  /*   getShapeFunctions(pt, N, Na, Nb); */

  /*   // Compute the derivative of X with respect to the */
  /*   // coordinate directions */
  /*   TacsScalar X[3], Xa[9]; */
  /*   planeJacobian(X, Xa, N, Na, Nb, Xpts); */

  /*   // Compute the determinant of Xa and the transformation */
  /*   TacsScalar J[4]; */
  /*   TacsScalar h = FElibrary::jacobian2d(Xa, J); */
  /*   h = h*weight; */

  /*   // Compute the strain B*u */
  /*   TacsScalar strain[NUM_STRESSES]; */
  /*   evalStrain(strain, J, Na, Nb, vars); */
  /*   // Get alpha */
  /*   TacsScalar aT = stiff->getThermalAlpha(); */
  /*   // ------------------------------------------------------------ */
  /*   // Compute the product of alpha and the change in temperature T */
  /*   TacsScalar T = 0.0; */
  /*   for (int i = 0; i < NUM_NODES; i++){ */
  /*     T += N[i]*vars[3*i+2]; */
  /*   } */
  /*   aT *= T; */

  /*   // Compute the corresponding stress D*(B*u-e_{th}) */
  /*   TacsScalar stress[NUM_STRESSES], estrain[NUM_STRESSES]; */
  /*   estrain[0] = strain[0]-aT; */
  /*   estrain[1] = strain[1]-aT; */
  /*   estrain[2] = strain[2]; */
  /*   // ------------------------------------------------------------ */
  /*   stiff->calculateStress(pt, estrain, stress); */
  /*   // Compute the contribution from the potential energy */
  /*   // (B*u-e_{th})^{T}*D*(B*u-e_{th}) */
  /*   // Pe += 0.5*h*(stress[0]*strain[0] + stress[1]*strain[1] +
   * stress[2]*strain[2]); */
  /*   Pe += 0.5*h*(stress[0]*estrain[0] + stress[1]*estrain[1] +
   * stress[2]*estrain[2]); */
  /*   //
   * -------------------------------------------------------------------------------
   */
  /*   // Add thermal energy due to heat conduction */
  /*   strain[0] = strain[1] = 0.0; */
  /*   // Compute the term B*T */
  /*   evalBT(strain, J, Na,Nb, vars); */
  /*   stress[0] = stress[1] = 0.0; */
  /*   // Compute the term D*B*T */
  /*   stiff->calculateConduction(pt, strain,stress); */
  /*   // Add the term T^{T}*B^{T}*D*B*T*0.5 */
  /*   Pe += 0.5*h*(stress[0]*strain[0] + stress[1]*strain[1] +
   * stress[2]*strain[2]);     */
  /*   //
   * ------------------------------------------------------------------------------
   */

  /*   // Get value of the mass/area at this point */
  /*   TacsScalar mass; */
  /*   stiff->getPointwiseMass(pt, &mass); */

  /*   // Compute the contribution from the kinetic energy */
  /*   TacsScalar dUdt[2]; */
  /*   getDisplacement(dUdt, N, dvars); */
  /*   Te += 0.5*h*mass*(dUdt[0]*dUdt[0] + dUdt[1]*dUdt[1]); */
  /* } */

  // Set the output values
  *_Te = Te;
  *_Pe = Pe;
}
// --------------------------------------------------------------------------
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
void TACS2DCoupledThermoElement<NUM_NODES>::addResidual(
    double time, TacsScalar *res, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES * 2 * NUM_NODES];

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
    // ---------------------------------------------------
    // Compute the product D*B*u
    TacsScalar stress[NUM_STRESSES];
    stiff->calculateStress(pt, strain, stress);
    // -----------------------------------------------------
    // Get the derivative of the strain with respect to the nodal
    // displacements
    getBmat(B, J, Na, Nb, vars);
    TacsScalar *b = B;
    // Input into corresponding residual position
    for (int i = 0; i < NUM_NODES; i++) {
      // Contribution to u
      res[3 * i] +=
          h * (b[0] * stress[0] + b[1] * stress[1] + b[2] * stress[2]);
      b += NUM_STRESSES;
      // Contribution to v
      res[3 * i + 1] +=
          h * (b[0] * stress[0] + b[1] * stress[1] + b[2] * stress[2]);
      b += NUM_STRESSES;
    }

    b = NULL;
    // Get value of the mass/area at this point
    TacsScalar mass;
    stiff->getPointwiseMass(pt, &mass);

    // Add the contribution from the inertial terms
    TacsScalar d2Udt2[2];
    getDisplacement(d2Udt2, N, ddvars);
    for (int i = 0; i < NUM_NODES; i++) {
      res[3 * i] += h * mass * N[i] * d2Udt2[0];
      res[3 * i + 1] += h * mass * N[i] * d2Udt2[1];
    }

    if (conduction) {
      // Add heat transfer component to the residual
      // ----------------------------------------------------
      // Vector phi [3x1]
      TacsScalar phi[] = {1.0, 1.0, 0.0};
      stress[0] = stress[1] = stress[2] = 0.0;
      // Compute D*alpha*phi [3x1]
      stiff->calculateThermal(pt, phi, stress);
      // Get the derivative of the strain with respect to the nodal
      // displacements [NUM_STRESSESx2*NUM_NODES]
      memset(B, 0.0, NUM_STRESSES * 2 * NUM_NODES * sizeof(TacsScalar));
      getBmat(B, J, Na, Nb, vars);
      b = NULL;
      b = B;
      TacsScalar q[2 * NUM_NODES];
      memset(q, 0.0, 2 * NUM_NODES * sizeof(TacsScalar));
      // Compute B^{T}*D*alpha*phi [2xNUM_NODESx1]
      for (int i = 0; i < 2 * NUM_NODES; i++) {
        q[i] += (b[0] * stress[0] + b[1] * stress[1]);
        b += NUM_STRESSES;
      }
      // Get the temperature N^{T}*T
      TacsScalar T = 0.0;
      for (int i = 0; i < NUM_NODES; i++) {
        T += N[i] * vars[3 * i + 2];
      }
      // Input into the displacement components of the residual
      for (int i = 0; i < NUM_NODES; i++) {
        res[3 * i] -= q[2 * i] * h * T;
        res[3 * i + 1] -= q[2 * i + 1] * h * T;
      }
      // Evaluate the product B*dT [2x1]
      strain[0] = strain[1] = 0.0;
      evalBT(strain, J, Na, Nb, vars);
      // Get the [2xNUM_NODES] B matrix
      getBmatTemp(B, J, Na, Nb, vars);
      // Evaluate the product D*B*dT [2x1]
      stress[0] = stress[1] = 0.0;
      stiff->calculateConduction(pt, strain, stress);
      b = NULL;
      b = B;
      // Evaluate the product B^{T}*D*B*T to input into the
      // temperature components
      for (int i = 0; i < NUM_NODES; i++) {
        res[3 * i + 2] += h * (b[0] * stress[0] + b[1] * stress[1]);
        b += 2;
      }
    }  // end if conduction
  }    // end for int n = 0; n < numGauss
}
/*
  Add the Jacobian of the governing equations - the exact Jacobian of
  the residual expression.

  output:
  mat:     the element Jacobian of size 12x12

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
void TACS2DCoupledThermoElement<NUM_NODES>::addJacobian(
    double time, TacsScalar mat[], double alpha, double beta, double gamma,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES * 2 * NUM_NODES];

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
      // Compute the strain which is the product B*u
      TacsScalar strain[NUM_STRESSES];
      evalStrain(strain, J, Na, Nb, vars);

      // Compute the product D*B*u
      TacsScalar stress[NUM_STRESSES];
      stiff->calculateStress(pt, strain, stress);

      // Get the derivative of the strain with respect to the nodal
      // displacements
      getBmat(B, J, Na, Nb, vars);
      // Fill-in the upper portion of the matrix
      TacsScalar *bj = B;
      for (int j = 0; j < NUM_NODES; j++) {
        for (int jj = 0; jj < 2; jj++) {
          // Compute the product of D*B at the given point
          TacsScalar bs[NUM_STRESSES];
          stiff->calculateStress(pt, bj, bs);
          bj += NUM_STRESSES;
          for (int i = 0; i < NUM_NODES; i++) {
            for (int ii = 0; ii < 2; ii++) {
              const TacsScalar *bi = &B[3 * (2 * i + ii)];
              mat[NUM_VARIABLES * (3 * j + jj) + (3 * i + ii)] +=
                  alpha * h * (bi[0] * bs[0] + bi[1] * bs[1] + bi[2] * bs[2]);
            }  // end for int ii = 0; ii < 2
          }    // end for int i = 0; i < NUM_NODES
        }      // end for int ii = 0; ii < 2
      }        // end for int j = 0; j < NUM_NODES
    }          // end if alpha != 0.0
    // ----------------------------------------------------------------
    if (gamma != 0.0) {
      // Get value of the mass/area at this point
      TacsScalar mass;
      stiff->getPointwiseMass(pt, &mass);
      // Add the contributions from the stiffness matrix
      TacsScalar scale = gamma * h * mass;
      for (int j = 0; j < NUM_NODES; j++) {
        for (int i = 0; i < NUM_NODES; i++) {
          mat[(3 * j) * NUM_VARIABLES + 3 * i] += scale * N[i] * N[j];
          mat[(3 * j + 1) * NUM_VARIABLES + 3 * i + 1] += scale * N[i] * N[j];
        }
      }
    }  // end if (gamma != 0.0)
    // ---------------------------------------------------------------
    // Add contribution of heat transfer to the Jacobian
    // Add heat conduction to the 'bottom right' of the Jacobian [H]
    // [NUM_NODESxNUM_NODES]
    if (conduction) {
      // Compute the nonsymmetric part of the Jacobian i.e. matrix L
      memset(B, 0.0, NUM_STRESSES * 2 * NUM_NODES * sizeof(TacsScalar));
      TacsScalar s[NUM_STRESSES], BDs[2 * NUM_NODES];
      memset(BDs, 0.0, 2 * NUM_NODES * sizeof(TacsScalar));
      // Get the coefficient of thermal expansion
      TacsScalar phi[] = {1.0, 1.0, 0.0};
      // Get B [NUM_STRESSESx2*NUM_NODES];
      getBmat(B, J, Na, Nb, vars);
      TacsScalar *bj = B;
      // Compute the vector D*phi*aT [NUM_STRESSESx1]
      stiff->calculateThermal(pt, phi, s);
      // Compute the vector B^{T}*D*phi*aT [2xNUM_NODESx1]
      for (int i = 0; i < 2 * NUM_NODES; i++) {
        BDs[i] += (bj[0] * s[0] + bj[1] * s[1]);
        bj += NUM_STRESSES;
      }
      for (int j = 0; j < NUM_NODES; j++) {
        for (int jj = 0; jj < 2; jj++) {
          for (int i = 0; i < NUM_NODES; i++) {
            for (int ii = 2; ii < 3; ii++) {
              mat[NUM_VARIABLES * (3 * j + jj) + (3 * i + ii)] -=
                  h * alpha * BDs[2 * j + jj] * N[i];
            }
          }
        }
      }
      // Get the B matrix [2xNUM_NODES] associated with the heat transfer
      // problem
      memset(B, 0.0, 2 * NUM_NODES * sizeof(TacsScalar));
      getBmatTemp(B, J, Na, Nb, vars);
      bj = B;
      for (int j = 0; j < NUM_NODES; j++) {
        for (int jj = 2; jj < 3; jj++) {
          // Compute the product of Dtemp*Btemp [2xNUM_NODES] at the given point
          TacsScalar bs[2];  // column of B
          stiff->calculateConduction(pt, bj, bs);
          bj += 2;
          for (int i = 0; i < NUM_NODES; i++) {
            for (int ii = 2; ii < 3; ii++) {
              const TacsScalar *bi = &B[2 * i];
              mat[NUM_VARIABLES * (3 * j + jj) + (3 * i + ii)] +=
                  alpha * h * (bi[0] * bs[0] + bi[1] * bs[1]);
            }  // end for int ii = 2; ii < 3
          }    // end for int i = 0; i < NUM_NODES
        }      // end for int jj = 2; jj < 3
      }
    }  // end if conduction
  }    // end for int n = 0; n < numGauss
}
/*
  Add the product of the adjoint vector times the derivative of the
  residuals multiplied by a scalar to the given derivative vector.

  output:
  dvSens:  the array of the product of adjoint vector times the
           derivative of the residual

  input:
  psi:     the element adjoint variables
  vars:    the element variables
  Xpts:    the element nodal locations in R^{3}
*/
template <int NUM_NODES>
void TACS2DCoupledThermoElement<NUM_NODES>::addAdjResProduct(
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
    // -------------------------------------------------
    // Split psi into its displacement and temperature components and
    // reassemble at the end of the function
    TacsScalar psi_u[2 * NUM_NODES], psi_t[NUM_NODES];
    int iu = 0, it = 0;
    for (int i = 0; i < NUM_NODES; i++) {
      psi_u[iu] = 1.0 * psi[3 * i];
      psi_u[iu + 1] = 1.0 * psi[3 * i + 1];
      psi_t[it] = 1.0 * psi[3 * i + 2];
      iu += 2;
      it++;
    }
    // Compute the product B*u [NUM_STRESSESx1]
    TacsScalar strain[NUM_STRESSES];
    evalStrain(strain, J, Na, Nb, vars);

    // Get the derivative of the strain with respect to the nodal
    // displacements [NUM_STRESSESx2*NUM_NODES]
    getBmat(B, J, Na, Nb, vars);

    // Compute the product of psi_u^{T}*B^{T}
    TacsScalar bpsi[NUM_STRESSES];
    memset(bpsi, 0, NUM_STRESSES * sizeof(TacsScalar));

    TacsScalar *b = B;
    TacsScalar *ps = psi_u;
    for (int i = 0; i < 2 * NUM_NODES; i++) {
      bpsi[0] += ps[0] * b[0];
      bpsi[1] += ps[0] * b[1];
      bpsi[2] += ps[0] * b[2];
      b += NUM_STRESSES;
      ps++;
    }
    // Compute the term alpha*psi_u^{T}*B^{T}*dC/dx*B*u (contribution
    // of structural model)
    stiff->addStressDVSens(pt, strain, scale * h, bpsi, dvSens, dvLen);
    // Add the contribution of the coupling terms dL/dx
    TacsScalar T = 0.0;
    // Compute the term N^{T}*T
    for (int i = 0; i < NUM_NODES; i++) {
      T += N[i] * vars[3 * i + 2];
    }
    // Compute the vector -{1,1,0}*N^{T}*T
    TacsScalar phi[] = {-T, -T, 0.0};
    // Add the term alpha*psi_u^{T}*B^{T}*dC/dx*a*{1,1,0}*N^{T}*T
    stiff->addThermalDVSens(pt, phi, scale * h, bpsi, dvSens, dvLen);
    // Add the contribution of the heat transfer terms
    strain[0] = strain[1] = 0.0;
    // Compute the term B*T
    evalBT(strain, J, Na, Nb, vars);
    // Get the B matrix [NUM_DIRECTIONxNUM_NODES] associated with the heat
    // transfer problem
    getBmatTemp(B, J, Na, Nb, vars);
    // Initialize the term psi_t^{T}*B^{T}
    memset(bpsi, 0, NUM_STRESSES * sizeof(TacsScalar));
    b = NULL;
    b = B;
    ps = NULL;
    ps = psi_t;
    for (int i = 0; i < NUM_NODES; i++) {
      bpsi[0] += ps[0] * b[0];
      bpsi[1] += ps[0] * b[1];
      b += 2;
      ps++;
    }
    // Compute the term alpha*psi_t^{T}*B^{T}*dC/dx*B*T (contribution
    // of heat transfer model)
    stiff->addConductionDVSens(pt, strain, scale * h, bpsi, dvSens, dvLen);
  }  // end for int n = 0; n < numGauss
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
void TACS2DCoupledThermoElement<NUM_NODES>::addAdjResXptProduct(
    double time, double scale, TacsScalar fXptSens[], const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[]) {
  /* // The shape functions associated with the element */
  /* double N[NUM_NODES]; */
  /* double Na[NUM_NODES], Nb[NUM_NODES]; */

  /* // The derivative of the stress with respect to the strain */
  /* TacsScalar B[NUM_STRESSES*NUM_VARIABLES]; */

  /* // Get the number of quadrature points */
  /* int numGauss = getNumGaussPts(); */

  /* for ( int n = 0; n < numGauss; n++ ){ */
  /*   // Retrieve the quadrature points and weights */
  /*   double pt[3]; */
  /*   double weight = getGaussWtsPts(n, pt); */

  /*   // Compute the element shape functions */
  /*   getShapeFunctions(pt, N, Na, Nb); */

  /*   // Compute the derivative of X with respect to the */
  /*   // coordinate directions */
  /*   TacsScalar X[3], Xa[4]; */
  /*   planeJacobian(X, Xa, N, Na, Nb, Xpts); */

  /*   // Compute the determinant of Xa and the transformation */
  /*   TacsScalar J[4]; */
  /*   TacsScalar h = FElibrary::jacobian2d(Xa, J); */
  /*   h = h*weight; */

  /*   // Compute the derivative of the determinant w.r.t. node locations */
  /*   TacsScalar hXptSens[3*NUM_NODES]; */
  /*   TacsScalar h = getDetJacobianXptSens(hXptSens, pt, XptSens); */

  /*   // Compute the strain */
  /*   TacsScalar strain[NUM_STRESSES]; */
  /*   evalStrain(strain, J, Na, Nb, vars); */

  /*   // Get the derivative of the strain with respect to the nodal */
  /*   // displacements */
  /*   getBmat(B, J, Na, Nb, vars); */
  /* } */
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
void TACS2DCoupledThermoElement<NUM_NODES>::getMatType(
    ElementMatrixType matType, TacsScalar mat[], const TacsScalar Xpts[],
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
          mat[(3 * i) * NUM_VARIABLES + 3 * j] += d;
          mat[(3 * i + 1) * NUM_VARIABLES + 3 * j + 1] += d;
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
}

/*
  Evaluate the determinant of the Jacobian for numerical integration

  returns: the determinant of the Jacobian

  input:
  pt:    the parametric point within the element
  Xpts:  the element nodes
*/
template <int NUM_NODES>
TacsScalar TACS2DCoupledThermoElement<NUM_NODES>::getDetJacobian(
    const double pt[], const TacsScalar Xpts[]) {
  // Compute the element shape functions
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES];
  getShapeFunctions(pt, N, Na, Nb);

  // Compute the derivative of the shape functions w.r.t. the
  // parametric locations
  TacsScalar X[3], Xa[9];
  planeJacobian(X, Xa, N, Na, Nb, Xpts);

  TacsScalar J[4];
  return FElibrary::jacobian2d(Xa, J);
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
TacsScalar TACS2DCoupledThermoElement<NUM_NODES>::getDetJacobianXptSens(
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
void TACS2DCoupledThermoElement<NUM_NODES>::getStrain(TacsScalar strain[],
                                                      const double pt[],
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

  // Compute the strain B*u
  evalStrain(strain, J, Na, Nb, vars);
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
void TACS2DCoupledThermoElement<NUM_NODES>::getEffStrain(
    TacsScalar strain[], const double pt[], const TacsScalar Xpts[],
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
void TACS2DCoupledThermoElement<NUM_NODES>::addStrainSVSens(
    TacsScalar strainSVSens[], const double pt[], const TacsScalar scale,
    const TacsScalar strainSens[], const TacsScalar Xpts[],
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
  for (int i = 0; i < NUM_NODES; i++) {
    strainSVSens[3 * i] +=
        scale *
        (strainSens[0] * b[0] + strainSens[1] * b[1] + strainSens[2] * b[2]);
    b += NUM_STRESSES;
    strainSVSens[3 * i + 1] +=
        scale *
        (strainSens[0] * b[0] + strainSens[1] * b[1] + strainSens[2] * b[2]);
    b += NUM_STRESSES;
  }
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
void TACS2DCoupledThermoElement<NUM_NODES>::addEffStrainSVSens(
    TacsScalar strainSVSens[], const double pt[], const TacsScalar scale,
    const TacsScalar strainSens[], const TacsScalar Xpts[],
    const TacsScalar vars[], int vars_j) {
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

  // --------------------------------------------------
  TacsScalar *b = B;
  for (int i = 0; i < NUM_NODES; i++) {
    // Derivative wrt u
    strainSVSens[3 * i] +=
        scale *
        (strainSens[0] * b[0] + strainSens[1] * b[1] + strainSens[2] * b[2]);
    b += NUM_STRESSES;
    // Derivative wrt v
    strainSVSens[3 * i + 1] +=
        scale *
        (strainSens[0] * b[0] + strainSens[1] * b[1] + strainSens[2] * b[2]);
    b += NUM_STRESSES;
    // Derivative wrt dT
    TacsScalar aT = stiff->getEffThermalAlpha(vars_j);

    aT *= N[i];
    strainSVSens[3 * i + 2] -=
        scale * (strainSens[0] * aT + strainSens[1] * aT);
  }
  // ----------------------------------------------------
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
void TACS2DCoupledThermoElement<NUM_NODES>::addStrainXptSens(
    TacsScalar strainXptSens[], const double pt[], const TacsScalar scale,
    const TacsScalar strainSens[], const TacsScalar Xpts[],
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
/*
  Evaluate the derivative of dT at the specified point using the provided set of
  variables

  output:
  strain:   the strain evaluate at the specific parametric point

  input:
  vars:     the element variable values
  Xpts:     the element nodal locations
*/
template <int NUM_NODES>
void TACS2DCoupledThermoElement<NUM_NODES>::getBT(TacsScalar strain[],
                                                  const double pt[],
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

  // Compute the strain B*dT
  evalBT(strain, J, Na, Nb, vars);
}

/*
  Compute the derivative of the point-wise BT multiplied by a
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
void TACS2DCoupledThermoElement<NUM_NODES>::addBTSVSens(
    TacsScalar strainSVSens[], const double pt[], const TacsScalar scale,
    const TacsScalar strainSens[], const TacsScalar Xpts[],
    const TacsScalar vars[]) {
  // The shape functions associated with the element
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[2 * NUM_NODES];

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
  getBmatTemp(B, J, Na, Nb, vars);
  TacsScalar *b = B;
  for (int i = 0; i < NUM_NODES; i++) {
    strainSVSens[3 * i + 2] +=
        scale * (strainSens[0] * b[0] + strainSens[1] * b[1]);
    b += 2;
  }
}

#endif  // TACS_2D_COUPLED_THERMO_ELEMENT_H
