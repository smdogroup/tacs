/*
=============================================================================

=============================================================================
@File    :   TACSBladeStiffenedShellConstitutive.cpp
@Date    :   2023/04/06
@Author  :   Alasdair Christison Gray
@Description : Constitutive model for a blade-stiffened shell. Based on the
bladeFSDT model from previous versions of TACS developed by Graeme Kennedy.
*/

// =============================================================================
// Standard Library Includes
// =============================================================================

// =============================================================================
// Extension Includes
// =============================================================================
#include "TACSBladeStiffenedShellConstitutive.h"

#include "TACSMaterialProperties.h"
#include "TACSShellConstitutive.h"

// ==============================================================================
// Constructor
// ==============================================================================

TACSBladeStiffenedShellConstitutive::TACSBladeStiffenedShellConstitutive(
    TACSOrthotropicPly* _panelPly, TACSOrthotropicPly* _stiffenerPly,
    TacsScalar _kcorr, TacsScalar _panelLength, int _panelLengthNum,
    TacsScalar _stiffenerPitch, int _stiffenerPitchNum,
    TacsScalar _stiffenerHeight, int _stiffenerHeightNum,
    TacsScalar _stiffenerThick, int _stiffenerThickNum, TacsScalar _panelThick,
    int _panelThickNum, int _numPanelPlies, TacsScalar _panelPlyAngles[],
    int _panelPlyFracNums[], int _numStiffenerPlies,
    TacsScalar _stiffenerPlyAngles[], int _stiffenerPlyFracNums[],
    TacsScalar _flangeFraction) {
  this->panelPly = _panelPly;
  this->panelPly->incref();

  this->stiffenerPly = _stiffenerPly;
  this->stiffenerPly->incref();

  this->kcorr = _kcorr;

  this->numDesignVars = 0;

  // --- Panel length values ---
  this->panelLength = _panelLength;
  this->panelLengthNum = _panelLengthNum;
  this->panelLengthLocalNum = _panelLengthNum;
  if (_panelLengthNum >= 0) {
    this->panelLengthLocalNum = this->numDesignVars;
    this->numDesignVars++;
  }
  this->panelLengthLowerBound = 0.000;
  this->panelLengthUpperBound = 1e20;

  // --- Stiffener pitch values ---
  this->stiffenerPitch = _stiffenerPitch;
  this->stiffenerPitchNum = _stiffenerPitchNum;
  this->stiffenerPitchLocalNum = _stiffenerPitchNum;
  if (_stiffenerPitchNum >= 0) {
    this->stiffenerPitchLocalNum = this->numDesignVars;
    this->numDesignVars++;
  }
  this->stiffenerPitchLowerBound = 1e-3;
  this->stiffenerPitchUpperBound = 1e20;

  // --- Stiffener height values ---
  this->stiffenerHeight = _stiffenerHeight;
  this->stiffenerHeightNum = _stiffenerHeightNum;
  this->stiffenerHeightLocalNum = _stiffenerHeightNum;
  if (_stiffenerHeightNum >= 0) {
    this->stiffenerHeightLocalNum = this->numDesignVars;
    this->numDesignVars++;
  }
  this->stiffenerHeightLowerBound = 1e-3;
  this->stiffenerHeightUpperBound = 1e20;

  // --- Stiffener thickness values ---
  this->stiffenerThick = _stiffenerThick;
  this->stiffenerThickNum = _stiffenerThickNum;
  this->stiffenerThickLocalNum = _stiffenerThickNum;
  if (_stiffenerThickNum >= 0) {
    this->stiffenerThickLocalNum = this->numDesignVars;
    this->numDesignVars++;
  }
  this->stiffenerThickLowerBound = 1e-4;
  this->stiffenerThickUpperBound = 1e20;

  // --- Panel thickness values ---
  this->panelThick = _panelThick;
  this->panelThickNum = _panelThickNum;
  this->panelThickLocalNum = _panelThickNum;
  if (_panelThickNum >= 0) {
    this->panelThickLocalNum = this->numDesignVars;
    this->numDesignVars++;
  }
  this->panelThickLowerBound = 1e-4;
  this->panelThickUpperBound = 1e20;

  // --- Panel ply values ---
  this->numPanelPlies = _numPanelPlies;
  this->panelPlyAngles = new TacsScalar[_numPanelPlies];
  this->panelPlyFracs = new TacsScalar[_numPanelPlies];
  this->panelPlyFracNums = new int[_numPanelPlies];
  this->panelPlyFracLocalNums = new int[_numPanelPlies];
  this->panelPlyFracLowerBounds = new TacsScalar[_numPanelPlies];
  this->panelPlyFracUpperBounds = new TacsScalar[_numPanelPlies];
  for (int ii = 0; ii < _numPanelPlies; ii++) {
    this->panelPlyAngles[ii] = _panelPlyAngles[ii];
    this->panelPlyFracNums[ii] = _panelPlyFracNums[ii];
    this->panelPlyFracLocalNums[ii] = _panelPlyFracNums[ii];
    if (_panelPlyFracNums[ii] >= 0) {
      this->panelPlyFracLocalNums[ii] = this->numDesignVars;
      this->numDesignVars++;
    }
    this->panelPlyFracs[ii] = 1.0 / _numPanelPlies;
    this->panelPlyFracLowerBounds[ii] = 0.1;
    this->panelPlyFracUpperBounds[ii] = 0.9;
  }

  // --- Stiffener ply values ---
  this->numStiffenerPlies = _numStiffenerPlies;
  this->stiffenerPlyAngles = new TacsScalar[_numStiffenerPlies];
  this->stiffenerPlyFracs = new TacsScalar[_numStiffenerPlies];
  this->stiffenerPlyFracNums = new int[_numStiffenerPlies];
  this->stiffenerPlyFracLocalNums = new int[_numStiffenerPlies];
  this->stiffenerPlyFracLowerBounds = new TacsScalar[_numStiffenerPlies];
  this->stiffenerPlyFracUpperBounds = new TacsScalar[_numStiffenerPlies];
  for (int ii = 0; ii < _numStiffenerPlies; ii++) {
    this->stiffenerPlyAngles[ii] = _stiffenerPlyAngles[ii];
    this->stiffenerPlyFracNums[ii] = _stiffenerPlyFracNums[ii];
    this->stiffenerPlyFracLocalNums[ii] = _stiffenerPlyFracNums[ii];
    if (_stiffenerPlyFracNums[ii] >= 0) {
      this->stiffenerPlyFracLocalNums[ii] = this->numDesignVars;
      this->numDesignVars++;
    }
    this->stiffenerPlyFracs[ii] = 1.0 / _numStiffenerPlies;
    this->stiffenerPlyFracLowerBounds[ii] = 0.1;
    this->stiffenerPlyFracUpperBounds[ii] = 0.9;
  }

  // --- Stiffener flange fraction ---
  this->flangeFraction = _flangeFraction;
}

// ==============================================================================
// Destructor
// ==============================================================================
TACSBladeStiffenedShellConstitutive::~TACSBladeStiffenedShellConstitutive() {
  this->panelPly->decref();
  this->stiffenerPly->decref();

  delete[] this->panelPlyAngles;
  this->panelPlyAngles = nullptr;

  delete[] this->panelPlyFracs;
  this->panelPlyFracs = nullptr;

  delete[] this->panelPlyFracNums;
  this->panelPlyFracNums = nullptr;

  delete[] this->panelPlyFracLocalNums;
  this->panelPlyFracLocalNums = nullptr;

  delete[] this->stiffenerPlyAngles;
  this->stiffenerPlyAngles = nullptr;

  delete[] this->stiffenerPlyFracs;
  this->stiffenerPlyFracs = nullptr;

  delete[] this->stiffenerPlyFracNums;
  this->stiffenerPlyFracNums = nullptr;

  delete[] this->stiffenerPlyFracLocalNums;
  this->stiffenerPlyFracLocalNums = nullptr;

  delete[] this->panelPlyFracLowerBounds;
  this->panelPlyFracLowerBounds = nullptr;

  delete[] this->panelPlyFracUpperBounds;
  this->panelPlyFracUpperBounds = nullptr;

  delete[] this->stiffenerPlyFracLowerBounds;
  this->stiffenerPlyFracLowerBounds = nullptr;

  delete[] this->stiffenerPlyFracUpperBounds;
  this->stiffenerPlyFracUpperBounds = nullptr;
}

// ==============================================================================
// Set non-default values
// ==============================================================================

void TACSBladeStiffenedShellConstitutive::setStiffenerPitchBounds(
    TacsScalar lowerBound, TacsScalar upperBound) {
  this->stiffenerPitchLowerBound = lowerBound;
  this->stiffenerPitchUpperBound = upperBound;
}

void TACSBladeStiffenedShellConstitutive::setStiffenerHeightBounds(
    TacsScalar lowerBound, TacsScalar upperBound) {
  this->stiffenerHeightLowerBound = lowerBound;
  this->stiffenerHeightUpperBound = upperBound;
}

void TACSBladeStiffenedShellConstitutive::setStiffenerThicknessBounds(
    TacsScalar lowerBound, TacsScalar upperBound) {
  this->stiffenerThickLowerBound = lowerBound;
  this->stiffenerThickUpperBound = upperBound;
}

void TACSBladeStiffenedShellConstitutive::setPanelThicknessBounds(
    TacsScalar lowerBound, TacsScalar upperBound) {
  this->panelThickLowerBound = lowerBound;
  this->panelThickUpperBound = upperBound;
}

void TACSBladeStiffenedShellConstitutive::setStiffenerPlyFractionBounds(
    TacsScalar lowerBounds[], TacsScalar upperBounds[]) {
  for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
    this->stiffenerPlyFracLowerBounds[ii] = lowerBounds[ii];
    this->stiffenerPlyFracUpperBounds[ii] = upperBounds[ii];
  }
}

void TACSBladeStiffenedShellConstitutive::setPanelPlyFractionBounds(
    TacsScalar lowerBounds[], TacsScalar upperBounds[]) {
  for (int ii = 0; ii < this->numPanelPlies; ii++) {
    this->panelPlyFracLowerBounds[ii] = lowerBounds[ii];
    this->panelPlyFracUpperBounds[ii] = upperBounds[ii];
  }
}

void TACSBladeStiffenedShellConstitutive::setStiffenerPlyFractions(
    TacsScalar plyFractions[]) {
  for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
    this->stiffenerPlyFracs[ii] = plyFractions[ii];
  }
}

void TACSBladeStiffenedShellConstitutive::setPanelPlyFractions(
    TacsScalar plyFractions[]) {
  for (int ii = 0; ii < this->numPanelPlies; ii++) {
    this->panelPlyFracs[ii] = plyFractions[ii];
  }
}

// ==============================================================================
// Setting/getting design variable information
// ==============================================================================
// Retrieve the global design variable numbers
int TACSBladeStiffenedShellConstitutive::getDesignVarNums(int elemIndex,
                                                          int dvLen,
                                                          int dvNums[]) {
  if (dvNums && dvLen >= this->numDesignVars) {
    int localIndex = 0;
    if (this->panelLengthNum >= 0) {
      dvNums[localIndex] = panelLengthNum;
      localIndex++;
    }
    if (this->stiffenerPitchNum >= 0) {
      dvNums[localIndex] = stiffenerPitchNum;
      localIndex++;
    }
    if (this->stiffenerHeightNum >= 0) {
      dvNums[localIndex] = stiffenerHeightNum;
      localIndex++;
    }
    if (this->stiffenerThickNum >= 0) {
      dvNums[localIndex] = stiffenerThickNum;
      localIndex++;
    }
    if (this->panelThickNum >= 0) {
      dvNums[localIndex] = panelThickNum;
      localIndex++;
    }
    for (int ii = 0; ii < this->numPanelPlies; ii++) {
      if (this->panelPlyFracNums[ii] >= 0) {
        dvNums[localIndex] = panelPlyFracNums[ii];
        localIndex++;
      }
    }
    for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
      if (this->stiffenerPlyFracNums[ii] >= 0) {
        dvNums[localIndex] = stiffenerPlyFracNums[ii];
        localIndex++;
      }
    }
  }
  return numDesignVars;
}

// Set the element design variable from the design vector
int TACSBladeStiffenedShellConstitutive::setDesignVars(int elemIndex, int dvLen,
                                                       const TacsScalar dvs[]) {
  if (dvLen >= this->numDesignVars) {
    int localIndex = 0;
    if (this->panelLengthNum >= 0) {
      this->panelLength = dvs[localIndex];
      localIndex++;
    }
    if (this->stiffenerPitchNum >= 0) {
      this->stiffenerPitch = dvs[localIndex];
      localIndex++;
    }
    if (this->stiffenerHeightNum >= 0) {
      this->stiffenerHeight = dvs[localIndex];
      localIndex++;
    }
    if (this->stiffenerThickNum >= 0) {
      this->stiffenerThick = dvs[localIndex];
      localIndex++;
    }
    if (this->panelThickNum >= 0) {
      this->panelThick = dvs[localIndex];
      localIndex++;
    }
    for (int ii = 0; ii < this->numPanelPlies; ii++) {
      if (this->panelPlyFracNums[ii] >= 0) {
        this->panelPlyFracs[ii] = dvs[localIndex];
        localIndex++;
      }
    }
    for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
      if (this->stiffenerPlyFracNums[ii] >= 0) {
        this->stiffenerPlyFracs[ii] = dvs[localIndex];
        localIndex++;
      }
    }
  }
  return this->numDesignVars;
}

// Get the element design variables values
int TACSBladeStiffenedShellConstitutive::getDesignVars(int elemIndex, int dvLen,
                                                       TacsScalar dvs[]) {
  if (dvLen >= this->numDesignVars) {
    int localIndex = 0;
    if (this->panelLengthNum >= 0) {
      dvs[localIndex] = this->panelLength;
      localIndex++;
    }
    if (this->stiffenerPitchNum >= 0) {
      dvs[localIndex] = this->stiffenerPitch;
      localIndex++;
    }
    if (this->stiffenerHeightNum >= 0) {
      dvs[localIndex] = this->stiffenerHeight;
      localIndex++;
    }
    if (this->stiffenerThickNum >= 0) {
      dvs[localIndex] = this->stiffenerThick;
      localIndex++;
    }
    if (this->panelThickNum >= 0) {
      dvs[localIndex] = this->panelThick;
      localIndex++;
    }
    for (int ii = 0; ii < this->numPanelPlies; ii++) {
      if (this->panelPlyFracNums[ii] >= 0) {
        dvs[localIndex] = this->panelPlyFracs[ii];
        localIndex++;
      }
    }
    for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
      if (this->stiffenerPlyFracNums[ii] >= 0) {
        dvs[localIndex] = this->stiffenerPlyFracs[ii];
        localIndex++;
      }
    }
  }
  return this->numDesignVars;
}

// Get the lower and upper bounds for the design variable values
int TACSBladeStiffenedShellConstitutive::getDesignVarRange(int elemIndex,
                                                           int dvLen,
                                                           TacsScalar lb[],
                                                           TacsScalar ub[]) {
  if (dvLen >= this->numDesignVars) {
    int localIndex = 0;
    if (this->panelLengthNum >= 0) {
      lb[localIndex] = this->panelLengthLowerBound;
      ub[localIndex] = this->panelLengthUpperBound;
      localIndex++;
    }
    if (this->stiffenerPitchNum >= 0) {
      lb[localIndex] = this->stiffenerPitchLowerBound;
      ub[localIndex] = this->stiffenerPitchUpperBound;
      localIndex++;
    }
    if (this->stiffenerHeightNum >= 0) {
      lb[localIndex] = this->stiffenerHeightLowerBound;
      ub[localIndex] = this->stiffenerHeightUpperBound;
      localIndex++;
    }
    if (this->stiffenerThickNum >= 0) {
      lb[localIndex] = this->stiffenerThickLowerBound;
      ub[localIndex] = this->stiffenerThickUpperBound;
      localIndex++;
    }
    if (this->panelThickNum >= 0) {
      lb[localIndex] = this->panelThickLowerBound;
      ub[localIndex] = this->panelThickUpperBound;
      localIndex++;
    }
    for (int ii = 0; ii < this->numPanelPlies; ii++) {
      if (this->panelPlyFracNums[ii] >= 0) {
        lb[localIndex] = this->panelPlyFracLowerBounds[ii];
        ub[localIndex] = this->panelPlyFracUpperBounds[ii];
        localIndex++;
      }
    }
    for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
      if (this->stiffenerPlyFracNums[ii] >= 0) {
        lb[localIndex] = this->stiffenerPlyFracLowerBounds[ii];
        ub[localIndex] = this->stiffenerPlyFracUpperBounds[ii];
        localIndex++;
      }
    }
  }
  return this->numDesignVars;
}

// ==============================================================================
// Evaluate mass properties
// ==============================================================================
// Evaluate the mass per unit area
TacsScalar TACSBladeStiffenedShellConstitutive::evalDensity(
    int elemIndex, const double pt[], const TacsScalar X[]) {
  // Density due to the panel = thickness * rho
  TacsScalar density = this->panelPly->getDensity() * this->panelThick;
  // Density due to the stiffeners = rho * A / pitch
  density += this->stiffenerPly->getDensity() * this->computeStiffenerArea() /
             this->stiffenerPitch;
  return density;
}

// Add the derivative of the density w.r.t. the design variables
void TACSBladeStiffenedShellConstitutive::addDensityDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    int dvLen, TacsScalar dfdx[]) {
  TacsScalar panelDensity = this->panelPly->getDensity();
  TacsScalar stiffenerDensity = this->stiffenerPly->getDensity();
  TacsScalar stiffenerArea = this->computeStiffenerArea();
  TacsScalar dAdt, dAdh;
  this->computeStiffenerAreaSens(dAdt, dAdh);

  if (this->stiffenerPitchLocalNum >= 0) {
    TacsScalar sp = this->stiffenerPitch;
    dfdx[this->stiffenerPitchLocalNum] +=
        -scale * stiffenerDensity * stiffenerArea / (sp * sp);
  }
  if (this->stiffenerHeightLocalNum >= 0) {
    dfdx[this->stiffenerHeightLocalNum] +=
        scale * stiffenerDensity * dAdh / this->stiffenerPitch;
  }
  if (this->stiffenerThickLocalNum >= 0) {
    dfdx[this->stiffenerThickLocalNum] +=
        scale * stiffenerDensity * dAdt / this->stiffenerPitch;
  }
  if (this->panelThickLocalNum >= 0) {
    dfdx[this->panelThickLocalNum] += scale * panelDensity;
  }
}

// Evaluate the mass moments
void TACSBladeStiffenedShellConstitutive::evalMassMoments(
    int elemIndex, const double pt[], const TacsScalar X[],
    TacsScalar moments[]) {
  TacsScalar sPitchInv = 1.0 / this->stiffenerPitch;
  TacsScalar sHeight = this->stiffenerHeight;
  TacsScalar sThick = this->stiffenerThick;
  TacsScalar pThick = this->panelThick;
  TacsScalar kf = this->flangeFraction;
  TacsScalar panelDensity = this->panelPly->getDensity();
  TacsScalar stiffenerDensity = this->stiffenerPly->getDensity();
  TacsScalar stiffenerArea = this->computeStiffenerArea();
  TacsScalar stiffenerOffset =
      this->computeStiffenerCentroidHeight() + 0.5 * pThick;
  TacsScalar stiffenerMOI = this->computeStiffenerMOI();

  moments[0] =
      panelDensity * pThick + stiffenerDensity * stiffenerArea * sPitchInv;

  // First moment of area is non-zero because the stiffener makes the section
  // asymmteric
  moments[1] = -stiffenerDensity * sHeight * sThick *
               (kf * pThick + kf * sThick + pThick + sHeight + 2 * sThick) *
               0.5 * sPitchInv;

  // Panel contribution to second moment of area
  moments[2] = 0.5 * panelDensity * pThick * pThick * pThick / 12.0;

  // Add stiffener MOI about it's own centroid + contribution from parallel axis
  // theorm
  moments[2] += (stiffenerMOI + 0.5 * stiffenerDensity * stiffenerArea *
                                    stiffenerOffset * stiffenerOffset) *
                sPitchInv;
}

/**
  Add the derivative of the pointwise mass times the given scalar

  @param elemIndex The local element index
  @param pt The parametric location
  @param X The point location
  @param scale Scale factor for the moments
  @param dvLen the length of the sensitivity array
  @param dfdx The sensitivity array
*/
void TACSBladeStiffenedShellConstitutive::addMassMomentsDVSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar scale[], int dvLen, TacsScalar dfdx[]) {
  TacsScalar sPitchInv = 1.0 / this->stiffenerPitch;
  TacsScalar sHeight = this->stiffenerHeight;
  TacsScalar sThick = this->stiffenerThick;
  TacsScalar pThick = this->panelThick;
  TacsScalar kf = this->flangeFraction;
  TacsScalar panelDensity = this->panelPly->getDensity();
  TacsScalar stiffenerDensity = this->stiffenerPly->getDensity();
  TacsScalar stiffenerArea = this->computeStiffenerArea();
  TacsScalar stiffenerOffset =
      this->computeStiffenerCentroidHeight() + 0.5 * pThick;
  TacsScalar stiffenerMOI = this->computeStiffenerMOI();

  TacsScalar dzdt, dzdh;
  this->computeStiffenerCentroidHeightSens(dzdt, dzdh);

  TacsScalar dAdt, dAdh;
  this->computeStiffenerAreaSens(dAdt, dAdh);

  TacsScalar dMOIdt, dMOIdh;
  this->computeStiffenerMOISens(dMOIdt, dMOIdh);

  // --- Stiffener pitch sensitivity ---
  if (this->stiffenerPitchLocalNum >= 0) {
    int ii = this->stiffenerPitchLocalNum;
    TacsScalar sPitchInv2 = sPitchInv * sPitchInv;
    // Density contribution
    dfdx[ii] -= scale[0] * stiffenerDensity * stiffenerArea * sPitchInv2;
    // First moment of area contribution
    dfdx[ii] += scale[1] * stiffenerDensity * sHeight * sThick *
                (kf * pThick + kf * sThick + pThick + sHeight + 2 * sThick) *
                0.5 * sPitchInv2;
    // Second moment of area contribution
    dfdx[ii] -= scale[2] *
                (stiffenerMOI + 0.5 * stiffenerDensity * stiffenerArea *
                                    stiffenerOffset * stiffenerOffset) *
                sPitchInv2;
  }

  // --- Stiffener height sensitivity ---
  if (this->stiffenerHeightLocalNum >= 0) {
    int ii = this->stiffenerHeightLocalNum;
    // --- Density contribution ---
    dfdx[ii] += scale[0] * stiffenerDensity * dAdh * sPitchInv;
    // --- First moment of area contribution ---
    dfdx[ii] +=
        scale[1] *
        (sThick * stiffenerDensity *
         (-kf * pThick - kf * sThick - pThick - 2 * sHeight - 2 * sThick) *
         0.5);
    // --- Second moment of area contribution ---
    // d/dh(MOI + 0.5*rho*A*z^2) = d/dh(MOI) + 0.5*rho*(dAdh*z^2 + 2*A*z*dzdh)
    dfdx[ii] +=
        scale[2] *
        (dMOIdh + 0.5 * stiffenerDensity * stiffenerOffset *
                      (stiffenerOffset * dAdh + 2.0 * stiffenerArea * dzdh));
  }

  // --- Stiffener thickness sensitivity ---
  if (this->stiffenerThickLocalNum >= 0) {
    int ii = this->stiffenerThickLocalNum;
    // --- Density contribution ---
    dfdx[ii] += scale[0] * stiffenerDensity * dAdt * sPitchInv;
    // --- First moment of area contribution ---
    dfdx[ii] +=
        scale[1] *
        (sHeight * stiffenerDensity *
         (-kf * pThick - 2 * kf * sThick - pThick - sHeight - 4 * sThick) / 2);
    // --- Second moment of area contribution ---
    // d/dt(MOI + 0.5*rho*A*z^2) = d/dt(MOI) + 0.5*rho*(dAdt*z^2 + 2*A*z*dzdt)
    dfdx[ii] +=
        scale[2] *
        (dMOIdt + 0.5 * stiffenerDensity * stiffenerOffset *
                      (stiffenerOffset * dAdt + 2.0 * stiffenerArea * dzdt));
  }

  // --- Panel thickness sensitivity ---
  if (this->panelThickLocalNum >= 0) {
    int ii = this->panelThickLocalNum;
    // Density contribution
    dfdx[ii] += scale[0] * panelDensity;
    // Second moment of area contribution
    dfdx[ii] += scale[2] * panelDensity * 0.25 * pThick * pThick;
  }
}

// ==============================================================================
// Helper functions for computing stiffener cross-section properties
// ==============================================================================
TacsScalar TACSBladeStiffenedShellConstitutive::computeStiffenerArea() {
  return (1 + this->flangeFraction) * this->stiffenerHeight *
         this->stiffenerThick;
}

void TACSBladeStiffenedShellConstitutive::computeStiffenerAreaSens(
    TacsScalar& dAdt, TacsScalar& dAdh) {
  dAdh = (1 + this->flangeFraction) * this->stiffenerThick;
  dAdt = (1 + this->flangeFraction) * this->stiffenerHeight;
}

TacsScalar
TACSBladeStiffenedShellConstitutive::computeStiffenerCentroidHeight() {
  return ((1 + this->flangeFraction) * this->stiffenerThick +
          0.5 * this->stiffenerHeight) /
         (1 + this->flangeFraction);
}

void TACSBladeStiffenedShellConstitutive::computeStiffenerCentroidHeightSens(
    TacsScalar& dzdt, TacsScalar& dzdh) {
  dzdh = 0.5 * (1 + this->flangeFraction);
  dzdt = (1 + 0.5 * this->flangeFraction) / (1 + this->flangeFraction);
}

TacsScalar TACSBladeStiffenedShellConstitutive::computeStiffenerIzz() {
  TacsScalar sh = this->stiffenerHeight;
  TacsScalar st = this->stiffenerThick;
  TacsScalar kf = this->flangeFraction;
  TacsScalar sh2 = sh * sh;
  TacsScalar st2 = st * st;
  TacsScalar kf2 = kf * kf;
  return sh * st *
         (kf2 * st2 + 4 * kf * sh2 + 6 * kf * sh * st + 4 * kf * st2 + sh2) /
         (12 * (kf + 1));
}

void TACSBladeStiffenedShellConstitutive::computeStiffenerIzzSens(
    TacsScalar& dIdt, TacsScalar& dIdh) {
  TacsScalar sh = this->stiffenerHeight;
  TacsScalar st = this->stiffenerThick;
  TacsScalar kf = this->flangeFraction;
  TacsScalar sh2 = sh * sh;
  TacsScalar st2 = st * st;
  TacsScalar kf2 = kf * kf;

  dIdt =
      sh *
      (3 * kf2 * st2 + 4 * kf * sh2 + 12 * kf * sh * st + 12 * kf * st2 + sh2) /
      (12 * (kf + 1));

  dIdh =
      st *
      (kf2 * st2 + 12 * kf * sh2 + 12 * kf * sh * st + 4 * kf * st2 + 3 * sh2) /
      (12 * (kf + 1));
}

TacsScalar TACSBladeStiffenedShellConstitutive::computeStiffenerJxx() {
  TacsScalar sh = this->stiffenerHeight;
  TacsScalar st = this->stiffenerThick;
  TacsScalar kf = this->flangeFraction;
  return sh * (st * st * st) * (kf + 1) / 3;
}

void TACSBladeStiffenedShellConstitutive::computeStiffenerJxxSens(
    TacsScalar& dJdt, TacsScalar& dJdh) {
  TacsScalar sh = this->stiffenerHeight;
  TacsScalar st = this->stiffenerThick;
  TacsScalar kf = this->flangeFraction;

  dJdt = sh * (st * st) * (kf + 1);
  dJdh = (st * st * st) * (kf + 1) / 3;
}

TacsScalar TACSBladeStiffenedShellConstitutive::computeStiffenerMOI() {
  TacsScalar rho = this->stiffenerPly->getDensity();
  TacsScalar sh = this->stiffenerHeight;
  TacsScalar st = this->stiffenerThick;
  TacsScalar kf = this->flangeFraction;
  TacsScalar A1 = sh * st;        // Area of the stiffener web
  TacsScalar A2 = kf * A1;        // Area of the stiffener flange
  TacsScalar z1 = st + 0.5 * sh;  // Centroid of the stiffener web
  TacsScalar z2 = 0.5 * st;       // Centroid of the stiffener flange
  TacsScalar zc =
      this->computeStiffenerCentroidHeight();  // Centroid of the whole
                                               // stiffener section
  // Offsets of each area from the centroid of the whole stiffener section
  TacsScalar dz1 = z1 - zc;
  TacsScalar dz2 = z2 - zc;

  TacsScalar MOI1 = rho * A1 * sh * sh /
                    12.0;  // MOI of the stiffener web about it's own centroid
  TacsScalar MOI2 = rho * A2 * st * st / 12.0;  // MOI of the stiffener flange
                                                // about it's own centroid

  // Parallel axis theorem to get the MOI of the stiffener web about its
  // centroid
  return MOI1 + MOI2 + 0.5 * rho * (A1 * dz1 * dz1 + A2 * dz2 * dz2);
}

void TACSBladeStiffenedShellConstitutive::computeStiffenerMOISens(
    TacsScalar& dMOIdt, TacsScalar& dMOIdh) {
  TacsScalar rho = this->stiffenerPly->getDensity();
  TacsScalar sh = this->stiffenerHeight;
  TacsScalar st = this->stiffenerThick;
  TacsScalar kf = this->flangeFraction;
  TacsScalar sh2 = sh * sh;
  TacsScalar st2 = st * st;
  TacsScalar kf2 = kf * kf;

  dMOIdt =
      sh * rho *
      (3 * kf2 * st2 + 4 * kf * sh2 + 12 * kf * sh * st + 12 * kf * st2 + sh2) /
      (12 * (kf + 1));

  dMOIdh =
      st * rho *
      (kf2 * st2 + 12 * kf * sh2 + 12 * kf * sh * st + 4 * kf * st2 + 3 * sh2) /
      (12 * (kf + 1));
}
