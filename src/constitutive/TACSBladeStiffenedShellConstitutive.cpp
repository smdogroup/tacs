/*
=============================================================================
Blade-Stiffened Shell Constitutive Model
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

const char* TACSBladeStiffenedShellConstitutive::constName =
    "TACSBladeStiffenedShellConstitutive";

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
    TacsScalar _panelPlyFracs[], int _panelPlyFracNums[],
    int _numStiffenerPlies, TacsScalar _stiffenerPlyAngles[],
    TacsScalar _stiffenerPlyFracs[], int _stiffenerPlyFracNums[],
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
    this->panelPlyFracs[ii] = _panelPlyFracs[ii];
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
    this->stiffenerPlyFracs[ii] = _stiffenerPlyFracs[ii];
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

  // Arrays for storing failure values, need values for each ply angle at the
  // top and bottom of the panel and at the tip of the stiffener
  this->panelPlyFailValues = new TacsScalar[2 * _numPanelPlies];
  this->stiffenerPlyFailValues = new TacsScalar[_numStiffenerPlies];
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

  delete[] this->panelPlyFailValues;
  this->panelPlyFailValues = nullptr;

  delete[] this->stiffenerPlyFailValues;
  this->stiffenerPlyFailValues = nullptr;
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
               (kf * pThick + kf * sThick + pThick + sHeight + 2.0 * sThick) *
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
                (kf * pThick + kf * sThick + pThick + sHeight + 2.0 * sThick) *
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
         (-kf * pThick - kf * sThick - pThick - 2.0 * sHeight - 2.0 * sThick) *
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
         (-kf * pThick - 2.0 * kf * sThick - pThick - sHeight - 4.0 * sThick) /
         2.0);

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
// Evaluate thermal properties
// ==============================================================================

// Evaluate the specific heat
TacsScalar TACSBladeStiffenedShellConstitutive::evalSpecificHeat(
    int elemIndex, const double pt[], const TacsScalar X[]) {
  return 0.0;
}

// ==============================================================================
// Compute stress/strain/stiffness
// ==============================================================================

// Evaluate the stress
void TACSBladeStiffenedShellConstitutive::evalStress(int elemIndex,
                                                     const double pt[],
                                                     const TacsScalar X[],
                                                     const TacsScalar e[],
                                                     TacsScalar s[]) {
  // Compute the panel stresses
  this->computePanelStress(e, s);

  // Compute the stiffener beam stresses then  transform them back to shell
  // stresses
  TacsScalar stiffenerStrain[TACSBeamConstitutive::NUM_STRESSES],
      stiffenerStress[TACSBeamConstitutive::NUM_STRESSES];
  this->transformStrain(e, stiffenerStrain);
  this->computeStiffenerStress(stiffenerStrain, stiffenerStress);
  this->addStiffenerStress(stiffenerStress, s);
}

// Add the derivative of the product of stress with a vector psi to dfdx
void TACSBladeStiffenedShellConstitutive::addStressDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], const TacsScalar psi[], int dvLen,
    TacsScalar dfdx[]) {
  // TODO: Implement this
}

// Evaluate the tangent stiffness
void TACSBladeStiffenedShellConstitutive::evalTangentStiffness(
    int elemIndex, const double pt[], const TacsScalar X[], TacsScalar C[]) {
  this->computeStiffness(C);
}

// ==============================================================================
// Compute failure criteria
// ==============================================================================
// Calculate the point-wise failure criteria
TacsScalar TACSBladeStiffenedShellConstitutive::evalFailure(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar e[]) {
  TacsScalar fail[2];
  fail[0] = this->computePanelFailure(e);

  TacsScalar stiffenerStrain[TACSBeamConstitutive::NUM_STRESSES];
  this->transformStrain(e, stiffenerStrain);
  fail[1] = this->computeStiffenerFailure(stiffenerStrain);

  // TODO: Add the buckling calculation here

  return ksAggregation(fail, 2, this->ksWeight);
}

// Evaluate the derivative of the failure criteria w.r.t. the strain
TacsScalar TACSBladeStiffenedShellConstitutive::evalFailureStrainSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar e[], TacsScalar sens[]) {
  // TODO: Implement this
  return 0.0;
}

// Add the derivative of the failure criteria w.r.t. the design variables
void TACSBladeStiffenedShellConstitutive::addFailureDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], int dvLen, TacsScalar dfdx[]) {
  // TODO: Implement this
}

// ==============================================================================
// Compute output quantities
// ==============================================================================
// Retrieve the design variable for plotting purposes
TacsScalar TACSBladeStiffenedShellConstitutive::evalDesignFieldValue(
    int elemIndex, const double pt[], const TacsScalar X[], int index) {
  // TODO: Implement this
  return 0.0;
}

// ==============================================================================
// Helper functions for computing the stiffness matrices
// ==============================================================================

// Compute the stiffness matrix
void TACSBladeStiffenedShellConstitutive::computeStiffness(TacsScalar C[]) {
  TacsScalar* A = &C[0];
  TacsScalar* B = &C[6];
  TacsScalar* D = &C[12];
  TacsScalar* As = &C[18];
  TacsScalar* drill = &C[21];

  // --- Zero out the C matrix ---
  memset(C, 0, this->NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));

  // Add the panel's stiffness contributions
  this->computePanelStiffness(C);

  // Compute the stiffener's beam stiffness matrix then transform it to a shell
  // stiffness matrix and add it
  TacsScalar Cstiff[TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
  this->computeStiffenerStiffness(Cstiff);
  this->addStiffenerStiffness(C, Cstiff);

  // --- Compute drilling stiffness ---
  drill[0] = 0.5 * DRILLING_REGULARIZATION * (As[0] + As[2]);
}

void TACSBladeStiffenedShellConstitutive::computeSmearedStiffness(
    TACSOrthotropicPly* ply, const int numPlies, const TacsScalar plyAngles[],
    const TacsScalar plyFractions[], TacsScalar Q[], TacsScalar ABar[]) {
  // Temporary arrays for the matrices of each ply
  TacsScalar QTemp[6];
  TacsScalar ABarTemp[3];

  // Zero out the matrices
  for (int entry = 0; entry < 6; entry++) {
    Q[entry] = 0.0;
  }
  for (int entry = 0; entry < 3; entry++) {
    ABar[entry] = 0.0;
  }

  // Add the weighted contribution of this ply to the whole laminate
  for (int plyNum = 0; plyNum < numPlies; plyNum++) {
    ply->calculateQbar(plyAngles[plyNum], QTemp);
    ply->calculateAbar(plyAngles[plyNum], ABarTemp);
    for (int entry = 0; entry < 6; entry++) {
      Q[entry] += plyFractions[plyNum] * QTemp[entry];
    }
    for (int entry = 0; entry < 3; entry++) {
      ABar[entry] += plyFractions[plyNum] * ABarTemp[entry];
    }
  }
}

// ==============================================================================
// Helper functions for transforming strains/stresses/stiffnesses between the
// panel and stiffener
// ==============================================================================

// Given the shell mid-plane strains, compute the equivalent beam strains at the
// stiffener centroid
void TACSBladeStiffenedShellConstitutive::transformStrain(
    const TacsScalar panelStrain[], TacsScalar stiffenerStrain[]) {
  // Compute the offset of the stiffener centroid from the shell mid-plane
  TacsScalar z =
      this->computeStiffenerCentroidHeight() - 0.5 * this->panelThick;

  // Axial strain (contains contribution from panel bending)
  stiffenerStrain[0] = panelStrain[0] + z * panelStrain[3];
  // Torsional strain (rotation around the 1 axis)
  stiffenerStrain[1] = -0.5 * panelStrain[5];
  // Vertical bending strain (rotation around 2 axis)
  stiffenerStrain[2] = panelStrain[3];
  // Horizontal bending strain (rotation around 3 axis)
  stiffenerStrain[3] = 0.0;
  // Vertical shear strain
  stiffenerStrain[4] = panelStrain[7];
  // Horizontal shear strain (contains contribution from panel twisting)
  stiffenerStrain[5] = 0.5 * (panelStrain[2] + z * panelStrain[5]);
}

// Add the contribution of the stiffener stress to the panel stress
void TACSBladeStiffenedShellConstitutive::addStiffenerStress(
    const TacsScalar stiffenerStress[], TacsScalar panelStress[]) {
  TacsScalar pInv = 1.0 / this->stiffenerPitch;
  // Compute the offset of the stiffener centroid from the shell mid-plane
  TacsScalar z =
      this->computeStiffenerCentroidHeight() - 0.5 * this->panelThick;

  panelStress[0] += stiffenerStress[0] * pInv;        // N11 = F1 / P
  panelStress[2] += 0.5 * stiffenerStress[4] * pInv;  // N12 = 1 / 2 * V12 / P
  panelStress[3] += (z * stiffenerStress[0] + stiffenerStress[2]) *
                    pInv;  // M11 = (z*F1 + M2) / P
  panelStress[5] += 0.5 * (-stiffenerStress[1] + z * stiffenerStress[4]) *
                    pInv;                       // M12 = 1/2 (- M1 + V12*z)/P
  panelStress[7] += stiffenerStress[4] * pInv;  // Q13 = V13 / P
}

// Add the contribution of the stiffener stiffness to the panel stiffness
void TACSBladeStiffenedShellConstitutive::addStiffenerStiffness(
    const TacsScalar stiffenerStiffness[], TacsScalar panelStiffness[]) {
  TacsScalar pInv = 1.0 / this->stiffenerPitch;
  // Compute the offset of the stiffener centroid from the shell mid-plane
  TacsScalar z =
      this->computeStiffenerCentroidHeight() - 0.5 * this->panelThick;

  // Some shorthand for the entries of the stiffness matrix
  TacsScalar* A = &panelStiffness[0];
  TacsScalar* B = &panelStiffness[6];
  TacsScalar* D = &panelStiffness[12];
  TacsScalar* As = &panelStiffness[18];
  const TacsScalar* Cs = stiffenerStiffness;

  // A:
  A[0] = pInv * (Cs[0]);
  A[2] = pInv * (Cs[5] / 2.0);
  A[5] = pInv * (Cs[20] / 4.0);

  // B:
  B[0] = pInv * (z * Cs[0] + Cs[2]);
  B[2] = pInv * (z * Cs[5] / 2.0 - Cs[1] / 2.0);
  B[5] = pInv * (z * Cs[20] / 4.0 - Cs[10] / 4.0);

  // D:
  D[0] = pInv * (z * (z * Cs[0] + Cs[2]) + z * Cs[2] + Cs[11]);
  D[2] =
      pInv * (z * (z * Cs[5] + Cs[14]) / 2.0 - z * Cs[1] / 2.0 - Cs[7] / 2.0);
  D[5] =
      pInv * (z * (z * Cs[20] - Cs[10]) / 4.0 - z * Cs[10] / 4.0 + Cs[6] / 4.0);

  // As:
  As[2] = pInv * (Cs[4, 4]);
}

// ==============================================================================
// Helper functions for computing the panel stress/stiffness/failure
// ==============================================================================
// In future, these methods should be replaced by calls to another shell
// constitutive model

// Compute the panel stress given the panel strains
void TACSBladeStiffenedShellConstitutive::computePanelStress(
    const TacsScalar strain[], TacsScalar stress[]) {
  TacsScalar C[this->NUM_TANGENT_STIFFNESS_ENTRIES];
  this->computePanelStiffness(C);

  TacsScalar* A = &C[0];
  TacsScalar* B = &C[6];
  TacsScalar* D = &C[12];
  TacsScalar* As = &C[18];
  TacsScalar drill = C[21];

  this->computeStress(A, B, D, As, drill, strain, stress);
}

void TACSBladeStiffenedShellConstitutive::computePanelStiffness(
    TacsScalar C[]) {
  TacsScalar* A = &C[0];
  TacsScalar* B = &C[6];
  TacsScalar* D = &C[12];
  TacsScalar* As = &C[18];

  // --- Zero out the C matrix ---
  memset(C, 0, this->NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));

  // Compute the smeared laminate properties
  TacsScalar QPanel[this->NUM_Q_ENTRIES], ABarPanel[this->NUM_ABAR_ENTRIES];

  this->computeSmearedStiffness(this->panelPly, this->numPanelPlies,
                                this->panelPlyAngles, this->panelPlyFracs,
                                QPanel, ABarPanel);

  // Add the panel's contributions to the A and D matrices
  TacsScalar t = this->panelThick;
  TacsScalar DFactor = t * t * t / 12.0;

  for (int ii = 0; ii < 6; ii++) {
    A[ii] += t * QPanel[ii];
    D[ii] += DFactor * QPanel[ii];
  }

  // Add the pane;'s contribution to the transverse shear matrix
  for (int ii = 0; ii < 3; ii++) {
    As[ii] += t * ABarPanel[ii];
  }
}

// Compute the failure criteria in the panel
TacsScalar TACSBladeStiffenedShellConstitutive::computePanelFailure(
    const TacsScalar strain[]) {
  TacsScalar t = this->panelThick;
  TacsScalar plyStrain[3];
  TacsScalar* fail = this->panelPlyFailValues;
  int numPly = this->numPanelPlies;

  // Compute the strain state at the top of the panel
  plyStrain[0] = strain[0] + 0.5 * t * strain[3];
  plyStrain[1] = strain[1] + 0.5 * t * strain[4];
  plyStrain[2] = strain[2] + 0.5 * t * strain[5];

  // Compute the failure criteria for each ply angle at this strain state
  for (int ii = 0; ii < numPly; ii++) {
    fail[ii] = this->panelPly->failure(this->panelPlyAngles[ii], plyStrain);
  }

  // Now repeat for the bottom of the panel
  plyStrain[0] = strain[0] - 0.5 * t * strain[3];
  plyStrain[1] = strain[1] - 0.5 * t * strain[4];
  plyStrain[2] = strain[2] - 0.5 * t * strain[5];

  for (int ii = 0; ii < numPly; ii++) {
    fail[numPly + ii] =
        this->panelPly->failure(this->panelPlyAngles[ii], plyStrain);
  }

  // Return the aggregated failure value
  return ksAggregation(fail, 2 * numPly, this->ksWeight);
}

// ==============================================================================
// Helper functions for computing the stiffner's strain/stress/stiffness
// ==============================================================================
// In future, these methods should be replaced by calls to a beam constitutive
// model

// Compute the beam stresses in the stiffener
void TACSBladeStiffenedShellConstitutive::computeStiffenerStress(
    const TacsScalar stiffenerStrain[], TacsScalar stiffenerStress[]) {
  int n = TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES;
  TacsScalar C[n];
  memset(C, 0, n * sizeof(TacsScalar));
  this->computeStiffenerStiffness(C);
  TACSBeamConstitutive::computeStress(stiffenerStrain, C, stiffenerStress);
}

// Compute the stiffener's beam stiffness matrix
void TACSBladeStiffenedShellConstitutive::computeStiffenerStiffness(
    TacsScalar C[]) {
  // --- Zero out the C matrix ---
  memset(
      C, 0,
      TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));

  TacsScalar A = this->computeStiffenerArea();
  TacsScalar Izz = this->computeStiffenerIzz();
  TacsScalar J = this->computeStiffenerJxx();

  // Compute the smeared laminate properties
  TacsScalar QPanel[this->NUM_Q_ENTRIES], ABarPanel[this->NUM_ABAR_ENTRIES];

  this->computeSmearedStiffness(this->panelPly, this->numPanelPlies,
                                this->panelPlyAngles, this->panelPlyFracs,
                                QPanel, ABarPanel);

  // Compute the effective elastic and shear moduli
  TacsScalar Q11 = QPanel[0];
  TacsScalar Q12 = QPanel[1];
  TacsScalar Q22 = QPanel[3];
  TacsScalar Q66 = QPanel[5];
  TacsScalar E = Q11 - Q12 * Q12 / Q22;
  TacsScalar G = Q66;

  // Populate the matrix
  C[0] = E * A;                 // C[0, 0]
  C[6] = G * J;                 // C[1, 1]
  C[11] = E * Izz;              // C[2, 2]
  C[18] = this->kcorr * G * A;  // C[4, 4]
  C[20] = this->kcorr * G * A;  // C[5, 5]
}

// Compute the failure criteria for the stiffener
TacsScalar TACSBladeStiffenedShellConstitutive::computeStiffenerFailure(
    const TacsScalar stiffenerStrain[]) {
  TACSOrthotropicPly* ply = this->stiffenerPly;

  // Compute the strain state at the tip of the stiffener
  TacsScalar zTipOffset = -(this->stiffenerHeight + this->stiffenerThick) -
                          this->computeStiffenerCentroidHeight();
  TacsScalar plyStrain[3];
  memset(plyStrain, 0, 3 * sizeof(TacsScalar));
  plyStrain[0] = stiffenerStrain[0] + zTipOffset * stiffenerStrain[2];

  // Compute the failure criteria at this strain state for each ply angle
  for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
    this->stiffenerPlyFailValues[ii] =
        ply->failure(this->stiffenerPlyAngles[ii], plyStrain);
  }

  // Returned the aggregated value over all plies
  return ksAggregation(this->stiffenerPlyFailValues, this->numStiffenerPlies,
                       this->ksWeight);
}

// ==============================================================================
// Helper functions for computing stiffener cross-section properties
// ==============================================================================
TacsScalar TACSBladeStiffenedShellConstitutive::computeStiffenerArea() {
  return (1.0 + this->flangeFraction) * this->stiffenerHeight *
         this->stiffenerThick;
}

void TACSBladeStiffenedShellConstitutive::computeStiffenerAreaSens(
    TacsScalar& dAdt, TacsScalar& dAdh) {
  dAdh = (1.0 + this->flangeFraction) * this->stiffenerThick;
  dAdt = (1.0 + this->flangeFraction) * this->stiffenerHeight;
}

TacsScalar
TACSBladeStiffenedShellConstitutive::computeStiffenerCentroidHeight() {
  return -((1.0 + this->flangeFraction) * this->stiffenerThick +
           0.5 * this->stiffenerHeight) /
         (1.0 + this->flangeFraction);
}

void TACSBladeStiffenedShellConstitutive::computeStiffenerCentroidHeightSens(
    TacsScalar& dzdt, TacsScalar& dzdh) {
  dzdh = -0.5 * (1.0 + this->flangeFraction);
  dzdt = -(1.0 + 0.5 * this->flangeFraction) / (1.0 + this->flangeFraction);
}

TacsScalar TACSBladeStiffenedShellConstitutive::computeStiffenerIzz() {
  TacsScalar sh = this->stiffenerHeight;
  TacsScalar st = this->stiffenerThick;
  TacsScalar kf = this->flangeFraction;
  TacsScalar sh2 = sh * sh;
  TacsScalar st2 = st * st;
  TacsScalar kf2 = kf * kf;
  return sh * st *
         (kf2 * st2 + 4.0 * kf * sh2 + 6.0 * kf * sh * st + 4.0 * kf * st2 +
          sh2) /
         (12.0 * (kf + 1.0));
}

void TACSBladeStiffenedShellConstitutive::computeStiffenerIzzSens(
    TacsScalar& dIdt, TacsScalar& dIdh) {
  TacsScalar sh = this->stiffenerHeight;
  TacsScalar st = this->stiffenerThick;
  TacsScalar kf = this->flangeFraction;
  TacsScalar sh2 = sh * sh;
  TacsScalar st2 = st * st;
  TacsScalar kf2 = kf * kf;

  dIdt = sh *
         (3.0 * kf2 * st2 + 4.0 * kf * sh2 + 12.0 * kf * sh * st +
          12.0 * kf * st2 + sh2) /
         (12.0 * (kf + 1.0));

  dIdh = st *
         (kf2 * st2 + 12.0 * kf * sh2 + 12.0 * kf * sh * st + 4.0 * kf * st2 +
          3.0 * sh2) /
         (12.0 * (kf + 1.0));
}

TacsScalar TACSBladeStiffenedShellConstitutive::computeStiffenerJxx() {
  TacsScalar sh = this->stiffenerHeight;
  TacsScalar st = this->stiffenerThick;
  TacsScalar kf = this->flangeFraction;
  return sh * (st * st * st) * (kf + 1.0) / 3.0;
}

void TACSBladeStiffenedShellConstitutive::computeStiffenerJxxSens(
    TacsScalar& dJdt, TacsScalar& dJdh) {
  TacsScalar sh = this->stiffenerHeight;
  TacsScalar st = this->stiffenerThick;
  TacsScalar kf = this->flangeFraction;

  dJdt = sh * (st * st) * (kf + 1.0);
  dJdh = (st * st * st) * (kf + 1.0) / 3.0;
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

  dMOIdt = sh * rho *
           (3.0 * kf2 * st2 + 4.0 * kf * sh2 + 12.0 * kf * sh * st +
            12.0 * kf * st2 + sh2) /
           (12.0 * (kf + 1.0));

  dMOIdh = st * rho *
           (kf2 * st2 + 12.0 * kf * sh2 + 12.0 * kf * sh * st + 4.0 * kf * st2 +
            3.0 * sh2) /
           (12.0 * (kf + 1.0));
}
