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

#ifndef TACS_CONSTITUTIVE_H
#define TACS_CONSTITUTIVE_H

#include "TACSObject.h"

/**
  This class defines basic linear constitutive behaviour for some
  elements in TACS. This includes linear

  The object defines mass, stiffness, failure and buckling behavior.
  These properties do not all need to be defined to implement the
  object.  The failure and buckling properties are not needed unless
  the object will be used by TACSFunction objects that evaluate
  failure or buckling criterion.

  The TACSConstitutive object also implements methods to visualize the
  failure and buckling envelopes. In these methods, two stress
  coordinate axes are defined and these are used to compute a 2D
  representation of the failure or buckling envelope.
*/
class TACSConstitutive : public TACSObject {
 public:
  TACSConstitutive() {}
  virtual ~TACSConstitutive() {}

  /**
    Set the object name
  */
  const char* getObjectName();

  /**
    Return the number of stress and strain components
  */
  virtual int getNumStresses() = 0;

  /**
    Get the number of design variables per "design node"
  */
  virtual int getDesignVarsPerNode() { return 1; }

  /**
    Retrieve the global design variable numbers

    Note when the dvNums argument is NULL, then the result is a query
    on the number of design variables and the array is not set.

    @param dvLen The length of the array dvNums
    @param dvNums An array of the design variable numbers
    @return The number of design variable numbers defined
  */
  virtual int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]) {
    return 0;
  }

  /**
    Set the design variables from the design vector

    @param elemIndex The local element index
    @param dvLen The length of the design array
    @param dvs The design variable values
    @return The number of design variable numbers defined
  */
  virtual int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]) {
    return 0;
  }

  /**
    Get the design variables values

    @param elemIndex The local element index
    @param dvLen The length of the design array
    @param dvs The design variable values
    @return The number of design variable numbers defined
  */
  virtual int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]) {
    return 0;
  }

  /**
    Get the lower and upper bounds for the design variable values

    @param elemIndex The local element index
    @param dvLen The length of the design array
    @param lowerBound The design variable lower bounds
    @param lowerBound The design variable upper bounds
    @return The number of design variable numbers defined
  */
  virtual int getDesignVarRange(int elemIndex, int dvLen,
                                TacsScalar lowerBound[],
                                TacsScalar upperBound[]) {
    return 0;
  }

  /**
    Evaluate the mass per unit length, area or volume for the element

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @return The density
  */
  virtual TacsScalar evalDensity(int elemIndex, const double pt[],
                                 const TacsScalar X[]) = 0;

  /**
    Evaluate the mass per unit length, area or volume for the element

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @param vars The state variable array
    @return The density
  */
  virtual TacsScalar evalDensity(int elemIndex, const double pt[],
                                 const TacsScalar X[],
                                 const TacsScalar vars[]) {
    return evalDensity(elemIndex, pt, X);
  }

  /**
    Add the derivative of the pointwise mass times the given scalar

    @param elemIndex The local element index
    @param scale Scale factor for the derivative
    @param pt The parametric location
    @param X The point location
    @param dvLen The length of the sensitivity array
    @param dfdx The sensitivity array
  */
  virtual void addDensityDVSens(int elemIndex, TacsScalar scale,
                                const double pt[], const TacsScalar X[],
                                int dvLen, TacsScalar dfdx[]) {}

  /**
    Add the derivative of the pointwise mass times the given scalar

    @param elemIndex The local element index
    @param scale Scale factor for the derivative
    @param pt The parametric location
    @param X The point location
    @param dvLen The length of the sensitivity array
    @param dfdx The sensitivity array
    @param vars The state variable array
  */
  virtual void addDensityDVSens(int elemIndex, TacsScalar scale,
                                const double pt[], const TacsScalar X[],
                                int dvLen, TacsScalar dfdx[],
                                const TacsScalar vars[]) {
    return addDensityDVSens(elemIndex, scale, pt, X, dvLen, dfdx);
  }

  /**
    Add the derivative of the pointwise mass with respect to state variables

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @param dfdu The state variable sensitivity array
  */
  virtual void addDensitySVSens(int elemIndex, const double pt[],
                                const TacsScalar X[], TacsScalar dfdu[]) {}

  /**
    Add the derivative of the pointwise mass with respect to state variables

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @param dfdu The state variable sensitivity array
    @param vars The state variable array
  */
  virtual void addDensitySVSens(int elemIndex, const double pt[],
                                const TacsScalar X[], TacsScalar dfdu[],
                                const TacsScalar vars[]) {
    return addDensitySVSens(elemIndex, pt, X, dfdu);
  }

  /**
    Evaluate the mass per unit length, area or volume for the element for
    the mass matrix evaluation.

    For some topology optimization problems the density computed for the mass
    matrix computation is different than the element density.

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @return The density
  */
  virtual TacsScalar evalMassMatrixDensity(int elemIndex, const double pt[],
                                           const TacsScalar X[]) {
    return evalDensity(elemIndex, pt, X);
  }

  /**
    Add the derivative of the pointwise mass times the given scalar

    @param elemIndex The local element index
    @param scale Scale factor for the derivative
    @param pt The parametric location
    @param X The point location
    @param dvLen The length of the sensitivity array
    @param dfdx The sensitivity array
  */
  virtual void addMassMatrixDensityDVSens(int elemIndex, TacsScalar scale,
                                          const double pt[],
                                          const TacsScalar X[], int dvLen,
                                          TacsScalar dfdx[]) {
    addDensityDVSens(elemIndex, scale, pt, X, dvLen, dfdx);
  }

  /**
    Evaluate the specific heat (heat capacity per unit mass)

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @return The specific heat of the material
  */
  virtual TacsScalar evalSpecificHeat(int elemIndex, const double pt[],
                                      const TacsScalar X[]) = 0;

  /**
    Evaluate the specific heat (heat capacity per unit mass)

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @param vars The state variable array
    @return The specific heat of the material
  */
  virtual TacsScalar evalSpecificHeat(int elemIndex, const double pt[],
                                      const TacsScalar X[],
                                      const TacsScalar vars[]) {
    return evalSpecificHeat(elemIndex, pt, X);
  }

  /**
    Add the derivative of the pointwise mass times the given scalar

    @param elemIndex The local element index
    @param scale Scale factor for the derivative
    @param pt The parametric location
    @param X The point location
    @param dvLen The length of the sensitivity array
    @param dfdx The sensitivity array
  */
  virtual void addSpecificHeatDVSens(int elemIndex, TacsScalar scale,
                                     const double pt[], const TacsScalar X[],
                                     int dvLen, TacsScalar dfdx[]) {}

  /**
    Add the derivative of the pointwise mass times the given scalar

    @param elemIndex The local element index
    @param scale Scale factor for the derivative
    @param pt The parametric location
    @param X The point location
    @param dvLen The length of the sensitivity array
    @param dfdx The sensitivity array
    @param vars The state variable array
  */
  virtual void addSpecificHeatDVSens(int elemIndex, TacsScalar scale,
                                     const double pt[], const TacsScalar X[],
                                     int dvLen, TacsScalar dfdx[],
                                     const TacsScalar vars[]) {
    return addSpecificHeatDVSens(elemIndex, scale, pt, X, dvLen, dfdx);
  }

  /**
    Add the specific heat derivative with respect to state variables

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @param dvLen The length of the sensitivity array
    @param dfdu The state variable sensitivity array
  */
  virtual void addSpecificHeatSVSens(int elemIndex, const double pt[],
                                     const TacsScalar X[], TacsScalar dfdu[]) {}

  /**
    Add the specific heat derivative with respect to state variables

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @param dvLen The length of the sensitivity array
    @param dfdu The state variable sensitivity array
    @param vars The state variable array
  */
  virtual void addSpecificHeatSVSens(int elemIndex, const double pt[],
                                     const TacsScalar X[], TacsScalar dfdu[],
                                     const TacsScalar vars[]) {
    return addSpecificHeatSVSens(elemIndex, pt, X, dfdu);
  }

  /**
    Return the stress as a function of the strain at the Gauss point

    @param elemIndex The local element index
    @param pt The parametric point within the element
    @param X The physical point location
    @param strain The strain evaluated at the point
    @param stress The components of the stress
  */
  virtual void evalStress(int elemIndex, const double pt[],
                          const TacsScalar X[], const TacsScalar strain[],
                          TacsScalar stress[]) = 0;

  /**
    Compute the tangent stiffness matrix

    Note that the format of the entries of the C matrix are constitutive
    implementation dependent. The matrix is not always of the same size.

    @param elemIndex The local element index
    @param pt The parametric point within the element
    @param X The physical point location
    @param C The components of the tangent stiffness
  */
  virtual void evalTangentStiffness(int elemIndex, const double pt[],
                                    const TacsScalar X[], TacsScalar C[]) = 0;

  /**
    Add the derivative of the stress times an input vector

    @param elemIndex The local element index
    @param scale A scalar factor
    @param pt The parametric point within the element
    @param X The physical point location
    @param psi The adjoint components of the strain
    @param dvLen The length of the design vector array
    @param dfdx The sensitivity vector
  */
  virtual void addStressDVSens(int elemIndex, TacsScalar scale,
                               const double pt[], const TacsScalar X[],
                               const TacsScalar strain[],
                               const TacsScalar psi[], int dvLen,
                               TacsScalar dfdx[]) {}

  /**
    Evaluate the tangent stiffness used for the geometric stiffness
    matrix computations.

    Note that by default this uses the original tangent stiffness
    computation, but some topology optimization techniques use a
    different definition.

    @param elemIndex The local element index
    @param pt The parametric point within the element
    @param X The physical point location
    @param C The components of the geometric tangent stiffness
  */
  virtual void evalGeometricTangentStiffness(int elemIndex, const double pt[],
                                             const TacsScalar X[],
                                             TacsScalar C[]) {
    evalTangentStiffness(elemIndex, pt, X, C);
  }

  /**
    Add the derivative of the geometric tangent stress matrix

    @param elemIndex The local element index
    @param pt The parametric point within the element
    @param X The physical point location
    @param scale A scalar factor
    @param psi The adjoint components of the strain
    @param dvLen The length of the design vector array
    @param dfdx The sensitivity vector
  */
  virtual void addGeometricTangentStressDVSens(int elemIndex, TacsScalar scale,
                                               const double pt[],
                                               const TacsScalar X[],
                                               const TacsScalar strain[],
                                               const TacsScalar psi[],
                                               int dvLen, TacsScalar dfdx[]) {
    addStressDVSens(elemIndex, scale, pt, X, strain, psi, dvLen, dfdx);
  }

  /**
    Evaluate the thermal strain at a point

    @param elemIndex The local element index
    @param pt The parametric point within the element
    @param X The physical point location
    @param theta The local change in temperature
    @param strain The components of the thermal strain
  */
  virtual void evalThermalStrain(int elemIndex, const double pt[],
                                 const TacsScalar X[], TacsScalar theta,
                                 TacsScalar strain[]) {
    memset(strain, 0, getNumStresses() * sizeof(TacsScalar));
  }

  /**
    Add the contributions from the derivatives to the thermal strain

    @param elemIndex The local element index
    @param pt The parametric point within the element
    @param X The physical point location
    @param theta The local change in temperature
    @param psi Multiplier vector (same size as the strain)
    @param dvLen The length of the sensitivity array
    @param dfdx The sensitivity
  */
  virtual void addThermalStrainDVSens(int elemIndex, const double pt[],
                                      const TacsScalar X[], TacsScalar theta,
                                      const TacsScalar psi[], int dvLen,
                                      TacsScalar dfdx[]) {}

  /**
    Given the thermal gradient, compute the heat flux

    @param elemIndex The local element index
    @param pt The parametric point within the element
    @param X The physical point location
  */
  virtual void evalHeatFlux(int elemIndex, const double pt[],
                            const TacsScalar X[], const TacsScalar grad[],
                            TacsScalar flux[]) {}

  /**
    Given the thermal gradient, compute the heat flux

    @param elemIndex The local element index
    @param pt The parametric point within the element
    @param X The physical point location
    @param vars The state variable array
  */
  virtual void evalHeatFlux(int elemIndex, const double pt[],
                            const TacsScalar X[], const TacsScalar grad[],
                            TacsScalar flux[], const TacsScalar vars[]) {
    return evalHeatFlux(elemIndex, pt, X, grad, flux);
  }

  /**
    Compute the coefficients of the heat flux

    @param elemIndex The local element index
    @param pt The parametric point within the element
    @param X The physical point location
    @param C The tangent heat flux matrix
  */
  virtual void evalTangentHeatFlux(int elemIndex, const double pt[],
                                   const TacsScalar X[], TacsScalar C[]) {}

  /**
    Compute the coefficients of the heat flux

    @param elemIndex The local element index
    @param pt The parametric point within the element
    @param X The physical point location
    @param C The tangent heat flux matrix
    @param vars The state variable array
  */
  virtual void evalTangentHeatFlux(int elemIndex, const double pt[],
                                   const TacsScalar X[], TacsScalar C[],
                                   const TacsScalar vars[]) {
    return evalTangentHeatFlux(elemIndex, pt, X, C);
  }

  /**
    Add the derivative of the heat flux to the sensitivity array

    @param elemIndex The local element index
    @param scale A scalar factor
    @param pt The parametric point within the element
    @param X The physical point location
    @param psi Multiplier vector (same size as the strain)
    @param dvLen The length of the sensitivity array
    @param dfdx The sensitivity
  */
  virtual void addHeatFluxDVSens(int elemIndex, TacsScalar scale,
                                 const double pt[], const TacsScalar X[],
                                 const TacsScalar grad[],
                                 const TacsScalar psi[], int dvLen,
                                 TacsScalar dfdx[]) {}

  /**
    Add the derivative of the heat flux to the sensitivity array

    @param elemIndex The local element index
    @param scale A scalar factor
    @param pt The parametric point within the element
    @param X The physical point location
    @param psi Multiplier vector (same size as the strain)
    @param dvLen The length of the sensitivity array
    @param dfdx The sensitivity
    @param vars The state variable array
  */
  virtual void addHeatFluxDVSens(int elemIndex, TacsScalar scale,
                                 const double pt[], const TacsScalar X[],
                                 const TacsScalar grad[],
                                 const TacsScalar psi[], int dvLen,
                                 TacsScalar dfdx[], const TacsScalar vars[]) {
    return addHeatFluxDVSens(elemIndex, scale, pt, X, grad, psi, dvLen, dfdx);
  }

  /**
    Evaluate the failure index at a quadrature point

    @param elemIndex The local element index
    @param pt The parametric point
    @param X The physical node location
    @param strain The strain value
    @return The failure index value
  */
  virtual TacsScalar evalFailure(int elemIndex, const double pt[],
                                 const TacsScalar X[],
                                 const TacsScalar strain[]) {
    return 0.0;
  }

  /**
    Evaluate the failure index at a quadrature point

    @param elemIndex The local element index
    @param pt The parametric point
    @param X The physical node location
    @param strain The strain value
    @param sens The derivative of the failure index w.r.t. the strain
    @return The failure index value
  */
  virtual TacsScalar evalFailureStrainSens(int elemIndex, const double pt[],
                                           const TacsScalar X[],
                                           const TacsScalar strain[],
                                           TacsScalar sens[]) {
    memset(sens, 0, getNumStresses() * sizeof(TacsScalar));
    return 0.0;
  }

  /**
    Add the derivative of the failure w.r.t. design variables

    @param elemIndex The local element index
    @param scale Scale The derivative of the failure index w.r.t. the strain
    @param pt The parametric point
    @param X The physical node location
    @param strain The strain value
    @param dvLen The length of the design vector
    @param dfdx The sensitivity contribution
  */
  virtual void addFailureDVSens(int elemIndex, TacsScalar scale,
                                const double pt[], const TacsScalar X[],
                                const TacsScalar strain[], int dvLen,
                                TacsScalar dfdx[]) {}

  /**
    Evaluate a design field (if defined) at the given quadrature point

    @param elemIndex The local element index
    @param pt The parametric point
    @param X The physical node location
    @param index The design field index
    @return The value of the design field
  */
  virtual TacsScalar evalDesignFieldValue(int elemIndex, const double pt[],
                                          const TacsScalar X[], int index) {
    return 0.0;
  }

  /**
    Compute a two-dimensional representation of the failure envelope
  */
  void getFailureEnvelope(int npts, int elemIndex, const double pt[],
                          const TacsScalar X[], const TacsScalar x_stress[],
                          const TacsScalar y_stress[], TacsScalar x_vals[],
                          TacsScalar y_vals[]);

 private:
  static const char* constName;
};

#endif  // TACS_CONSTITUTIVE_H
