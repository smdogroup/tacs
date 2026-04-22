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

#ifndef TACS_COMPOSITE_TUBE_BEAM_CONSTITUTIVE_H
#define TACS_COMPOSITE_TUBE_BEAM_CONSTITUTIVE_H

#include "TACSBeamConstitutive.h"

/*
  Timoshenko beam constitutive object for a hollow circular composite tube.

  Uses Classical Lamination Theory (CLT) to compute smeared effective axial
  modulus E_eff and shear modulus G_eff from ply-level properties and layup
  angles.  All plies are assumed equal thickness.  The layup must be symmetric
  and balanced (B=0, A16=A26=0) for the smearing to be valid.

  DV convention (mirrors TACSIsoTubeBeamConstitutive):
    inner (d)  — inner diameter [m], DV index dNum
    wall  (tw) — diametric wall thickness = OD - ID [m], DV index twNum
  Outer diameter: D_o = inner + wall.

  Failure modes aggregated by a KS smooth-max (ks_weight = 100):
    1. Fiber-direction compression at the outer fibre: -E11*eps_1 / X_c
    2. Fiber-direction tension  at the outer fibre:  +E11*eps_1 / X_t
       (skipped when X_t <= 0)
    3. Euler column buckling: Nx / Pcr, where Pcr = pi^2 * E_eff * Ia / (Kb*Lb)^2
       (skipped when buckleLengthFactor == 0, i.e. Kb == 0)

  Analytic adjoints are provided for all strain and DV sensitivities.
*/
class TACSCompositeTubeBeamConstitutive : public TACSBeamConstitutive {
 public:
  /*
    Constructor.

    Ply properties: E11, E22, G12, nu12 — orthotropic single-ply properties
    rho  — ply density (kg/m^3)
    X_c  — longitudinal compressive strength (Pa), must be > 0
    X_t  — longitudinal tensile strength (Pa); pass 0 to skip tensile check
    layup_angles — array of n_plies ply angles in radians
    n_plies      — number of plies (all equal thickness)

    Geometry DVs:
    inner_init, wall_init — initial inner diameter and diametric wall thickness [m]
    inner_dv, wall_dv     — global DV indices (-1 = fixed)
    inner_lb/ub, wall_lb/ub — DV bounds

    Buckling:
    buckle_length        — member length Lb [m]
    buckle_length_factor — effective length factor Kb (0 = disable buckling check)
  */
  TACSCompositeTubeBeamConstitutive(
      TacsScalar E11, TacsScalar E22, TacsScalar G12, TacsScalar nu12,
      TacsScalar rho, TacsScalar X_c, TacsScalar X_t,
      const TacsScalar *layup_angles, int n_plies,
      TacsScalar inner_init, TacsScalar wall_init,
      int inner_dv, int wall_dv,
      TacsScalar inner_lb, TacsScalar inner_ub,
      TacsScalar wall_lb, TacsScalar wall_ub,
      TacsScalar buckle_length = 1.0,
      TacsScalar buckle_length_factor = 0.0,
      int x_dv = -1,
      TacsScalar p_penalty = 3.0,
      TacsScalar eps_m = 1e-9);
  ~TACSCompositeTubeBeamConstitutive();

  // Retrieve the global design variable numbers
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);

  // Set the element design variable from the design vector
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);

  // Get the element design variables values
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);

  // Get the lower and upper bounds for the design variable values
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]);

  // Evaluate the material density
  TacsScalar evalDensity(int elemIndex, const double pt[],
                         const TacsScalar X[]);

  // Add the derivative of the density
  void addDensityDVSens(int elemIndex, TacsScalar scale, const double pt[],
                        const TacsScalar X[], int dvLen, TacsScalar dfdx[]);

  // Evaluate the mass moments
  void evalMassMoments(int elemIndex, const double pt[], const TacsScalar X[],
                       TacsScalar moments[]);

  // Add the sensitivity of the mass moments
  void addMassMomentsDVSens(int elemIndex, const double pt[],
                            const TacsScalar X[], const TacsScalar scale[],
                            int dvLen, TacsScalar dfdx[]);

  // Evaluate the specific heat (returns 0 — not defined for composites here)
  TacsScalar evalSpecificHeat(int elemIndex, const double pt[],
                              const TacsScalar X[]);

  // Evaluate the stress resultants
  void evalStress(int elemIndex, const double pt[], const TacsScalar X[],
                  const TacsScalar strain[], TacsScalar stress[]);

  // Evaluate the tangent stiffness matrix
  void evalTangentStiffness(int elemIndex, const double pt[],
                            const TacsScalar X[], TacsScalar C[]);

  // Add the contribution of stress to the DV sensitivity
  void addStressDVSens(int elemIndex, TacsScalar scale, const double pt[],
                       const TacsScalar X[], const TacsScalar strain[],
                       const TacsScalar psi[], int dvLen, TacsScalar dfdx[]);

  // Calculate the point-wise failure criteria
  TacsScalar evalFailure(int elemIndex, const double pt[], const TacsScalar X[],
                         const TacsScalar e[]);

  // Evaluate the derivative of the failure criteria w.r.t. the strain
  TacsScalar evalFailureStrainSens(int elemIndex, const double pt[],
                                   const TacsScalar X[], const TacsScalar e[],
                                   TacsScalar sens[]);

  // Add the derivative of the failure criteria w.r.t. the design variables
  void addFailureDVSens(int elemIndex, TacsScalar scale, const double pt[],
                        const TacsScalar X[], const TacsScalar strain[],
                        int dvLen, TacsScalar dfdx[]);

  // The name of the constitutive object
  const char *getObjectName();

  // Retrieve the design variable for plotting purposes
  TacsScalar evalDesignFieldValue(int elemIndex, const double pt[],
                                  const TacsScalar X[], int index);

 private:
  // CLT-smeared effective properties (precomputed at construction, const)
  TacsScalar E_eff;   // axial modulus [Pa]
  TacsScalar G_eff;   // shear modulus [Pa]
  TacsScalar nu_eff;  // effective Poisson ratio

  // Ply-level strength properties
  TacsScalar E11;   // fibre-direction modulus [Pa] (used in failure criterion)
  TacsScalar rho;   // ply density [kg/m^3]
  TacsScalar X_c;   // fibre-direction compressive strength [Pa]
  TacsScalar X_t;   // fibre-direction tensile strength [Pa] (0 = disabled)

  // Geometry DVs (same convention as TACSIsoTubeBeamConstitutive)
  TacsScalar inner, wall;  // inner diameter, diametric wall thickness [m]
  int innerDV, wallDV;
  TacsScalar innerLb, innerUb;
  TacsScalar wallLb, wallUb;

  // Euler buckling parameters
  TacsScalar buckleLength, buckleLengthFactor;

  TacsScalar ks_weight;

  // SIMP topology variables
  int xDV;
  TacsScalar x_val, p_penalty, eps_m;

  static const char *constName;
};

#endif  // TACS_COMPOSITE_TUBE_BEAM_CONSTITUTIVE_H
