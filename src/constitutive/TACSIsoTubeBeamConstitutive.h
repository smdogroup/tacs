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

#ifndef TACS_ISO_TUBE_BEAM_CONSTITUTIVE_H
#define TACS_ISO_TUBE_BEAM_CONSTITUTIVE_H

/*
  Base class for the Timoshenko beam constitutive object
*/

#include "TACSBeamConstitutive.h"
#include "TACSMaterialProperties.h"

class TACSIsoTubeBeamConstitutive : public TACSBeamConstitutive {
 public:
  TACSIsoTubeBeamConstitutive( TACSMaterialProperties *properties,
                               TacsScalar inner_diameter,
                               TacsScalar wall_thickness,
                               int inner_diameter_dv );
  ~TACSIsoTubeBeamConstitutive(){}

  int getDesignVarNums( int elemIndex, int dvLen, int dvNums[] );
  int setDesignVars( int elemIndex,
                     int dvLen, const TacsScalar dvs[] );
  int getDesignVars( int elemIndex,
                     int dvLen, TacsScalar dvs[] );
  int getDesignVarRange( int elemIndex, int dvLen,
                         TacsScalar lowerBound[],
                         TacsScalar upperBound[] );
  void evalMassMoments( int elemIndex,
                        const double pt[],
                        const TacsScalar X[],
                        TacsScalar moments[] );
  void addMassMomentsDVSens( int elemIndex,
                             const double pt[],
                             const TacsScalar X[],
                             const TacsScalar scale[],
                             int dvLen, TacsScalar dfdx[] );
  TacsScalar evalDensity( int elemIndex,
                          const double pt[],
                          const TacsScalar X[] );
  void addDensityDVSens( int elemIndex,
                         TacsScalar scale,
                         const double pt[],
                         const TacsScalar X[],
                         int dvLen, TacsScalar dfdx[] );
  void evalStress( int elemIndex,
                   const double pt[],
                   const TacsScalar X[],
                   const TacsScalar strain[],
                   TacsScalar stress[] );
  void evalTangentStiffness( int elemIndex,
                             const double pt[],
                             const TacsScalar X[],
                             TacsScalar C[] );
  void addStressDVSens( int elemIndex,
                        TacsScalar scale,
                        const double pt[],
                        const TacsScalar X[],
                        const TacsScalar strain[],
                        const TacsScalar psi[],
                        int dvLen, TacsScalar dfdx[] );
  TacsScalar evalFailure( int elemIndex,
                          const double pt[],
                          const TacsScalar X[],
                          const TacsScalar strain[] );
  TacsScalar evalFailureStrainSens( int elemIndex,
                                    const double pt[],
                                    const TacsScalar X[],
                                    const TacsScalar strain[],
                                    TacsScalar sens[] );
  void addFailureDVSens( int elemIndex,
                         TacsScalar scale,
                         const double pt[],
                         const TacsScalar X[],
                         const TacsScalar strain[],
                         int dvLen, TacsScalar dfdx[] );
  TacsScalar evalDesignFieldValue( int elemIndex,
                                   const double pt[],
                                   const TacsScalar X[],
                                   int index );
 private:
  TACSMaterialProperties *props;

  TacsScalar inner, wall;
  int innerDV, wallDV;
  TacsScalar lowerInner, upperInner;
  TacsScalar lowerWall, upperWall;
};

#endif // TACS_ISO_TUBE_BEAM_CONSTITUTIVE_H
