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
  TACSConstitutive(){}
  virtual ~TACSConstitutive(){}

  /**
    Set the object name
  */
  const char* getObjectName();

  /**
    Return the number of stress and strain components
  */
  virtual int getNumStresses() = 0;

  /**
    Retrieve the global design variable numbers

    Note when the dvNums argument is NULL, then the result is a query
    on the number of design variables and the array is not set.

    @param dvLen The length of the array dvNums
    @param dvNums An array of the design variable numbers
    @return The number of design variable numbers defined
  */
  virtual int getDesignVarNums( int elemIndex, int dvLen, int dvNums[] ){
    return 0;
  }

  /**
    Set the design variables from the design vector

    @param elemIndex The local element index
    @param dvLen The length of the design array
    @param dvs The design variable values
  */
  virtual void setDesignVars( int elemIndex,
                              int dvLen, const TacsScalar dvs[] ){}

  /**
    Get the design variables values

    @param elemIndex The local element index
    @param dvLen The length of the design array
    @param dvs The design variable values
  */
  virtual void getDesignVars( int elemIndex,
                              int dvLen, TacsScalar dvs[] ){}

  /**
    Get the lower and upper bounds for the design variable values

    @param elemIndex The local element index
    @param dvLen The length of the design array
    @param lowerBound The design variable lower bounds
    @param lowerBound The design variable upper bounds
  */
  virtual void getDesignVarRange( int elemIndex, int dvLen,
                                  TacsScalar lowerBound[],
                                  TacsScalar upperBound[] ){}

  /**
    Return the stress as a function of the strain at the Gauss point

    @param elemIndex The local element index
    @param pt The parametric point within the element
    @param X The physical point location
    @param strain The strain evaluated at the point
    @param stress The components of the stress
  */
  virtual void evalStress( int elemIndex,
                           const double pt[],
                           const TacsScalar X[],
                           const TacsScalar strain[],
                           TacsScalar stress[] ) = 0;

  /**
    Compute the tangent stiffness matrix

    Note that the format of the entries of the C matrix are constitutive
    implementation dependent. The matrix is not always of the same size.

    @param elemIndex The local element index
    @param pt The parametric point within the element
    @param X The physical point location
    @param C The components of the tangent stiffness
  */
  virtual void evalTangentStiffness( int elemIndex,
                                     const double pt[],
                                     const TacsScalar X[],
                                     TacsScalar C[] ) = 0;

  /**
    Add the derivative of the stress times an input vector to dvSens

    @param elemIndex The local element index
    @param pt The parametric point within the element
    @param X The physical point location
    @param scale The scalar factor
    @param psi The adjoint components of the strain
    @param dvLen The length of the design vector array
    @param dvSens The sensitivity vector
  */
  virtual void addStressDVSens( int elemIndex,
                                const double pt[],
                                const TacsScalar X[],
                                const TacsScalar strain[],
                                TacsScalar scale,
                                const TacsScalar psi[],
                                int dvLen, TacsScalar dvSens[] ){}

  /**
    Evaluate the thermal strain at a point

    @param elemIndex The local element index
    @param pt The parametric point within the element
    @param X The physical point location
    @param strain The components of the thermal strain
  */
  virtual void evalThermalStrain( int elemIndex,
                                  const double pt[],
                                  const TacsScalar X[],
                                  TacsScalar strain[] ){
    memset(strain, 0, getNumStresses()*sizeof(TacsScalar));
  }

  /**
    Evaluate the mass per unit length, area or volume for the element

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @return The density
  */
  virtual TacsScalar evalDensity( int elemIndex,
                                  const double pt[],
                                  const TacsScalar X[] ) = 0;

  /**
    Add the derivative of the pointwise mass times the given scalar

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @param scale Scale factor for the derivative
    @param dvLen the length of the sensitivity array
    @param dvSens The sensitivity array
  */
  virtual void addDensityDVSens( int elemIndex,
                                 const double pt[],
                                 const TacsScalar X[],
                                 const TacsScalar scale,
                                 int dvLen, TacsScalar dvSens[] ){}

  /**
    Evaluate the failure index at a quadrature point

    @param pt The parametric point
    @param X The physical node location
    @param strain the strain value
    @return The failure index value
  */
  virtual TacsScalar failure( int elemIndex,
                              const double pt[],
                              const TacsScalar X[],
                              const TacsScalar strain[] ){
    return 0.0;
  }

  /**
    Evaluate the failure index at a quadrature point

    @param pt The parametric point
    @param X The physical node location
    @param strain the strain value
    @param sens The derivative of the failure index w.r.t. the strain
    @return The failure index value
  */
  virtual TacsScalar failureStrainSens( int elemIndex,
                                        const double pt[],
                                        const TacsScalar X[],
                                        const TacsScalar strain[],
                                        TacsScalar sens[] ){
    memset(sens, 0, getNumStresses()*sizeof(TacsScalar));
    return 0.0;
  }

  /**
    Add the derivative of the failure w.r.t. design variables

    @param pt The parametric point
    @param X The physical node location
    @param strain the strain value
    @param scale Scale The derivative of the failure index w.r.t. the strain
    @param dvLen The length of the design vector
    @param dvSens The sensitivity contribution
  */
  virtual void addFailureDVSens( int elemIndex,
                                 const double pt[],
                                 const TacsScalar X[],
                                 const TacsScalar strain[],
                                 TacsScalar scale,
                                 int dvLen, TacsScalar dvSens[] ){}

  /**
    Write out a two-dimensional representation of the failure envelope
  */
  void writeFailureEnvelope( const char *file_name, int npts,
                             int elemIndex,
                             const double pt[],
                             const TacsScalar X[],
                             const TacsScalar x_stress[],
                             const TacsScalar y_stress[] );

 private:
  static const char *constName;
};

#endif // TACS_CONSTITUTIVE_H
