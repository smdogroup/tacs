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
    TacsScalar _stiffenerPitch, int _stiffenerPitchNum, TacsScalar _panelThick,
    int _panelThickNum, int _numPanelPlies, TacsScalar _panelPlyAngles[],
    TacsScalar _panelPlyFracs[], int _panelPlyFracNums[],
    TacsScalar _stiffenerHeight, int _stiffenerHeightNum,
    TacsScalar _stiffenerThick, int _stiffenerThickNum, int _numStiffenerPlies,
    TacsScalar _stiffenerPlyAngles[], TacsScalar _stiffenerPlyFracs[],
    int _stiffenerPlyFracNums[], TacsScalar _flangeFraction) {
  this->panelPly = _panelPly;
  this->panelPly->incref();

  this->stiffenerPly = _stiffenerPly;
  this->stiffenerPly->incref();

  this->kcorr = _kcorr;

  this->numDesignVars = 0;

  // --- General DVs ---
  // --- Panel length values ---
  this->panelLength = _panelLength;
  this->panelLengthNum = _panelLengthNum;
  this->panelLengthLocalNum = 0;
  if (_panelLengthNum >= 0) {
    this->panelLengthLocalNum = this->numDesignVars;
    this->numDesignVars++;
  }
  this->panelLengthLowerBound = 0.000;
  this->panelLengthUpperBound = 1e20;

  // --- Stiffener pitch values ---
  this->stiffenerPitch = _stiffenerPitch;
  this->stiffenerPitchNum = _stiffenerPitchNum;
  this->stiffenerPitchLocalNum = 0;
  if (_stiffenerPitchNum >= 0) {
    this->stiffenerPitchLocalNum = this->numDesignVars;
    this->numDesignVars++;
  }
  this->stiffenerPitchLowerBound = 1e-3;
  this->stiffenerPitchUpperBound = 1e20;

  // --- Panel DVs ---
  // --- Panel thickness values ---
  this->panelDVStartNum = this->numDesignVars;
  this->panelThick = _panelThick;
  this->panelThickNum = _panelThickNum;
  this->panelThickLocalNum = 0;
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
    this->panelPlyFracLocalNums[ii] = 0;
    if (_panelPlyFracNums[ii] >= 0) {
      this->panelPlyFracLocalNums[ii] = this->numDesignVars;
      this->numDesignVars++;
    }
    this->panelPlyFracs[ii] = 1.0 / _numPanelPlies;
    this->panelPlyFracLowerBounds[ii] = 0.1;
    this->panelPlyFracUpperBounds[ii] = 0.9;
  }

  // --- Stiffener DVs ---
  // --- Stiffener height values ---
  this->stiffenerDVStartNum = this->numDesignVars;
  this->stiffenerHeight = _stiffenerHeight;
  this->stiffenerHeightNum = _stiffenerHeightNum;
  this->stiffenerHeightLocalNum = 0;
  if (_stiffenerHeightNum >= 0) {
    this->stiffenerHeightLocalNum = this->numDesignVars;
    this->numDesignVars++;
  }
  this->stiffenerHeightLowerBound = 1e-3;
  this->stiffenerHeightUpperBound = 1e20;

  // --- Stiffener thickness values ---
  this->stiffenerThick = _stiffenerThick;
  this->stiffenerThickNum = _stiffenerThickNum;
  this->stiffenerThickLocalNum = 0;
  if (_stiffenerThickNum >= 0) {
    this->stiffenerThickLocalNum = this->numDesignVars;
    this->numDesignVars++;
  }
  this->stiffenerThickLowerBound = 1e-4;
  this->stiffenerThickUpperBound = 1e20;

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
    this->stiffenerPlyFracLocalNums[ii] = 0;
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

  // --- Work arrays, these are created to avoid needing to allocate memory
  // inside compute methods like evalFailure and evalFailureStrainSens ---

  // Arrays for storing failure values, need values for each ply angle at the
  // top and bottom of the panel and at the tip of the stiffener
  this->panelPlyFailValues = new TacsScalar[2 * _numPanelPlies];
  this->stiffenerPlyFailValues = new TacsScalar[_numStiffenerPlies];

  // Arrays for storing failure strain sensitivities
  this->panelPlyFailStrainSens = new TacsScalar*[2 * _numPanelPlies];
  this->stiffenerPlyFailStrainSens = new TacsScalar*[_numStiffenerPlies];
  for (int ii = 0; ii < 2 * _numPanelPlies; ii++) {
    this->panelPlyFailStrainSens[ii] = new TacsScalar[this->NUM_STRESSES];
  }
  for (int ii = 0; ii < _numStiffenerPlies; ii++) {
    this->stiffenerPlyFailStrainSens[ii] =
        new TacsScalar[TACSBeamConstitutive::NUM_STRESSES];
  }

  // Arrays for storing failure DV sensitivities
  this->panelPlyFailDVSens = new TacsScalar[2 * this->numPanelPlies];
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

  for (int ii = 0; ii < 2 * this->numPanelPlies; ii++) {
    delete[] this->panelPlyFailStrainSens[ii];
    this->panelPlyFailStrainSens[ii] = nullptr;
  }
  delete[] this->panelPlyFailStrainSens;
  this->panelPlyFailStrainSens = nullptr;

  for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
    delete[] this->stiffenerPlyFailStrainSens[ii];
    this->stiffenerPlyFailStrainSens[ii] = nullptr;
  }
  delete[] this->stiffenerPlyFailStrainSens;
  this->stiffenerPlyFailStrainSens = nullptr;

  delete[] this->panelPlyFailDVSens;
  this->panelPlyFailDVSens = nullptr;
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
    if (this->panelLengthNum >= 0) {
      dvNums[this->panelLengthLocalNum] = panelLengthNum;
    }
    if (this->stiffenerPitchNum >= 0) {
      dvNums[this->stiffenerPitchLocalNum] = stiffenerPitchNum;
    }
    if (this->panelThickNum >= 0) {
      dvNums[this->panelThickLocalNum] = panelThickNum;
    }
    for (int ii = 0; ii < this->numPanelPlies; ii++) {
      if (this->panelPlyFracNums[ii] >= 0) {
        dvNums[this->panelPlyFracLocalNums[ii]] = panelPlyFracNums[ii];
      }
    }
    if (this->stiffenerHeightNum >= 0) {
      dvNums[this->stiffenerHeightLocalNum] = stiffenerHeightNum;
    }
    if (this->stiffenerThickNum >= 0) {
      dvNums[this->stiffenerThickLocalNum] = stiffenerThickNum;
    }
    for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
      if (this->stiffenerPlyFracNums[ii] >= 0) {
        dvNums[this->stiffenerPlyFracLocalNums[ii]] = stiffenerPlyFracNums[ii];
      }
    }
  }
  return numDesignVars;
}

// Set the element design variable from the design vector
int TACSBladeStiffenedShellConstitutive::setDesignVars(int elemIndex, int dvLen,
                                                       const TacsScalar dvs[]) {
  if (dvLen >= this->numDesignVars) {
    if (this->panelLengthNum >= 0) {
      this->panelLength = dvs[this->panelLengthLocalNum];
    }
    if (this->stiffenerPitchNum >= 0) {
      this->stiffenerPitch = dvs[this->stiffenerPitchLocalNum];
    }
    if (this->panelThickNum >= 0) {
      this->panelThick = dvs[this->panelThickLocalNum];
    }
    for (int ii = 0; ii < this->numPanelPlies; ii++) {
      if (this->panelPlyFracNums[ii] >= 0) {
        this->panelPlyFracs[ii] = dvs[this->panelPlyFracLocalNums[ii]];
      }
    }
    if (this->stiffenerHeightNum >= 0) {
      this->stiffenerHeight = dvs[this->stiffenerHeightLocalNum];
    }
    if (this->stiffenerThickNum >= 0) {
      this->stiffenerThick = dvs[this->stiffenerThickLocalNum];
    }
    for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
      if (this->stiffenerPlyFracNums[ii] >= 0) {
        this->stiffenerPlyFracs[ii] = dvs[this->stiffenerPlyFracLocalNums[ii]];
      }
    }
  }
  return this->numDesignVars;
}

// Get the element design variables values
int TACSBladeStiffenedShellConstitutive::getDesignVars(int elemIndex, int dvLen,
                                                       TacsScalar dvs[]) {
  if (dvLen >= this->numDesignVars) {
    if (this->panelLengthNum >= 0) {
      dvs[this->panelLengthLocalNum] = this->panelLength;
    }
    if (this->stiffenerPitchNum >= 0) {
      dvs[this->stiffenerPitchLocalNum] = this->stiffenerPitch;
    }
    if (this->panelThickNum >= 0) {
      dvs[this->panelThickLocalNum] = this->panelThick;
    }
    for (int ii = 0; ii < this->numPanelPlies; ii++) {
      if (this->panelPlyFracNums[ii] >= 0) {
        dvs[this->panelPlyFracLocalNums[ii]] = this->panelPlyFracs[ii];
      }
    }
    if (this->stiffenerHeightNum >= 0) {
      dvs[this->stiffenerHeightLocalNum] = this->stiffenerHeight;
    }
    if (this->stiffenerThickNum >= 0) {
      dvs[this->stiffenerThickLocalNum] = this->stiffenerThick;
    }
    for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
      if (this->stiffenerPlyFracNums[ii] >= 0) {
        dvs[this->stiffenerPlyFracLocalNums[ii]] = this->stiffenerPlyFracs[ii];
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
    if (this->panelLengthNum >= 0) {
      lb[this->panelLengthLocalNum] = this->panelLengthLowerBound;
      ub[this->panelLengthLocalNum] = this->panelLengthUpperBound;
    }
    if (this->stiffenerPitchNum >= 0) {
      lb[this->stiffenerPitchLocalNum] = this->stiffenerPitchLowerBound;
      ub[this->stiffenerPitchLocalNum] = this->stiffenerPitchUpperBound;
    }
    if (this->panelThickNum >= 0) {
      lb[this->panelThickLocalNum] = this->panelThickLowerBound;
      ub[this->panelThickLocalNum] = this->panelThickUpperBound;
    }
    for (int ii = 0; ii < this->numPanelPlies; ii++) {
      if (this->panelPlyFracNums[ii] >= 0) {
        lb[this->panelPlyFracLocalNums[ii]] = this->panelPlyFracLowerBounds[ii];
        ub[this->panelPlyFracLocalNums[ii]] = this->panelPlyFracUpperBounds[ii];
      }
    }
    if (this->stiffenerHeightNum >= 0) {
      lb[this->stiffenerHeightLocalNum] = this->stiffenerHeightLowerBound;
      ub[this->stiffenerHeightLocalNum] = this->stiffenerHeightUpperBound;
    }
    if (this->stiffenerThickNum >= 0) {
      lb[this->stiffenerThickLocalNum] = this->stiffenerThickLowerBound;
      ub[this->stiffenerThickLocalNum] = this->stiffenerThickUpperBound;
    }
    for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
      if (this->stiffenerPlyFracNums[ii] >= 0) {
        lb[this->stiffenerPlyFracLocalNums[ii]] =
            this->stiffenerPlyFracLowerBounds[ii];
        ub[this->stiffenerPlyFracLocalNums[ii]] =
            this->stiffenerPlyFracUpperBounds[ii];
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

  // Compute the stiffener beam stresses then transform them back to shell
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
  TacsScalar fail[this->NUM_FAILURES];
  return this->computeFailureValues(e, fail);
}

// Compute the failure values for each failure mode of the stiffened panel
TacsScalar TACSBladeStiffenedShellConstitutive::computeFailureValues(
    const TacsScalar e[], TacsScalar fail[]) {
  fail[0] = this->computePanelFailure(e);

  TacsScalar stiffenerStrain[TACSBeamConstitutive::NUM_STRESSES];
  this->transformStrain(e, stiffenerStrain);
  fail[1] = this->computeStiffenerFailure(stiffenerStrain);

  // TODO: Add the buckling calculation here

  return ksAggregation(fail, this->NUM_FAILURES, this->ksWeight);
}

// Evaluate the derivative of the failure criteria w.r.t. the strain
TacsScalar TACSBladeStiffenedShellConstitutive::evalFailureStrainSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar e[], TacsScalar sens[]) {
  TacsScalar fails[2], dKSdf[2];
  // First compute the sensitivity of the panel failure value
  TacsScalar panelFailSens[this->NUM_STRESSES];
  fails[0] = this->evalPanelFailureStrainSens(e, panelFailSens);

  // And now for the stiffener failure value, first in terms of the beam
  // strains, and then transformed back to shell strains
  TacsScalar stiffenerStrain[TACSBeamConstitutive::NUM_STRESSES],
      stiffenerStrainSens[TACSBeamConstitutive::NUM_STRESSES],
      stiffenerFailSens[this->NUM_STRESSES];
  this->transformStrain(e, stiffenerStrain);
  fails[1] = this->evalStiffenerFailureStrainSens(stiffenerStrain,
                                                  stiffenerStrainSens);
  this->transformStrainSens(stiffenerStrainSens, stiffenerFailSens);

  // Compute the sensitivity of the aggregate failure value to the panel and
  // stiffener failure values
  TacsScalar fail = ksAggregationSens(fails, 2, this->ksWeight, dKSdf);

  // Compute the total sensitivity
  for (int i = 0; i < this->NUM_STRESSES; i++) {
    sens[i] = dKSdf[0] * panelFailSens[i] + dKSdf[1] * stiffenerFailSens[i];
  }

  return fail;
}

// Add the derivative of the failure criteria w.r.t. the design variables
void TACSBladeStiffenedShellConstitutive::addFailureDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], int dvLen, TacsScalar dfdx[]) {
  // TODO: Implement this
  const int numDV = this->numDesignVars;

  // First, compute the panel and stiffener failure values and then compute the
  // sensitivity of the aggregate failure value w.r.t. them
  TacsScalar fails[this->NUM_FAILURES], dKSdf[this->NUM_FAILURES];
  TacsScalar fail = this->computeFailureValues(strain, fails);
  ksAggregationSens(fails, this->NUM_FAILURES, this->ksWeight, dKSdf);

  // Sensitivity of the panel failure value to it's DVs
  this->addPanelFailureDVSens(strain, scale * dKSdf[0], dfdx);

  // Next, add the direct sensitivity of the stiffener failure value w.r.t DVs
  // Sensitivity of the panel failure value to it's DVs
  // TacsScalar stiffenerStrain[TACSBeamConstitutive::NUM_STRESSES];
  // this->transformStrain(strain, stiffenerStrain);
  // this->addStiffenerFailureDVSens(stiffenerStrain, dKSdf[1],
  //                                 this->numStiffenerDV,
  //                                 &dfdx[this->stiffenerDVStartNum]);

  // Finally, add the sensitivity of the stiffener failure value w.r.t. the DVs
  // due to the dependence of the stiffener strains on the DVs
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

void TACSBladeStiffenedShellConstitutive::transformStrainSens(
    const TacsScalar stiffenerStrainSens[], TacsScalar panelStrainSens[]) {
  memset(panelStrainSens, 0,
         TACSBeamConstitutive::NUM_STRESSES * sizeof(TacsScalar));
  TacsScalar z =
      this->computeStiffenerCentroidHeight() - 0.5 * this->panelThick;

  panelStrainSens[0] = stiffenerStrainSens[0];
  panelStrainSens[2] = 0.5 * stiffenerStrainSens[5];
  panelStrainSens[3] = stiffenerStrainSens[2] + z * stiffenerStrainSens[0];
  panelStrainSens[5] =
      0.5 * (z * stiffenerStrainSens[5] - stiffenerStrainSens[1]);
  panelStrainSens[7] = stiffenerStrainSens[4];
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

// Compute the derivative of the failure criteria in the panel w.r.t the DVs
TacsScalar TACSBladeStiffenedShellConstitutive::evalPanelFailureStrainSens(
    const TacsScalar strain[], TacsScalar sens[]) {
  TACSOrthotropicPly* ply = this->panelPly;
  const int numPlies = this->numPanelPlies;
  const int numStrain = TACSBeamConstitutive::NUM_STRESSES;
  TacsScalar** dFaildStrain = this->panelPlyFailStrainSens;
  TacsScalar* fails = this->panelPlyFailValues;
  const TacsScalar* angles = this->panelPlyAngles;
  const TacsScalar t = this->panelThick;

  // Compute the strain state at the top of the panel
  TacsScalar plyStrain[3];
  plyStrain[0] = strain[0] + 0.5 * t * strain[3];
  plyStrain[1] = strain[1] + 0.5 * t * strain[4];
  plyStrain[2] = strain[2] + 0.5 * t * strain[5];

  // Compute the failure criteria at this strain state for each ply angle and
  // the sensitivity of the failure criteria w.r.t the strains
  for (int ii = 0; ii < numPlies; ii++) {
    TacsScalar plyFailStrainSens[3];
    fails[ii] =
        ply->failureStrainSens(angles[ii], plyStrain, plyFailStrainSens);
    // Convert the sensitivity w.r.t the tip strain to the sensitivity w.r.t
    // beam strains
    memset(dFaildStrain[ii], 0, numStrain * sizeof(TacsScalar));
    dFaildStrain[ii][0] = plyFailStrainSens[0];
    dFaildStrain[ii][1] = plyFailStrainSens[1];
    dFaildStrain[ii][2] = plyFailStrainSens[2];
    dFaildStrain[ii][3] = 0.5 * t * plyFailStrainSens[0];
    dFaildStrain[ii][4] = 0.5 * t * plyFailStrainSens[1];
    dFaildStrain[ii][5] = 0.5 * t * plyFailStrainSens[2];
  }

  // Now repeat for the bottom of the panel
  plyStrain[0] = strain[0] - 0.5 * t * strain[3];
  plyStrain[1] = strain[1] - 0.5 * t * strain[4];
  plyStrain[2] = strain[2] - 0.5 * t * strain[5];

  for (int ii = 0; ii < numPlies; ii++) {
    TacsScalar plyFailStrainSens[3];
    fails[numPlies + ii] =
        ply->failureStrainSens(angles[ii], plyStrain, plyFailStrainSens);

    memset(dFaildStrain[numPlies + ii], 0, numStrain * sizeof(TacsScalar));
    dFaildStrain[numPlies + ii][0] = plyFailStrainSens[0];
    dFaildStrain[numPlies + ii][1] = plyFailStrainSens[1];
    dFaildStrain[numPlies + ii][2] = plyFailStrainSens[2];
    dFaildStrain[numPlies + ii][3] = -0.5 * t * plyFailStrainSens[0];
    dFaildStrain[numPlies + ii][4] = -0.5 * t * plyFailStrainSens[1];
    dFaildStrain[numPlies + ii][5] = -0.5 * t * plyFailStrainSens[2];
  }

  TacsScalar fail = ksAggregationSensProduct(
      fails, 2 * numPlies, numStrain, this->ksWeight, dFaildStrain, sens);

  return fail;
}

// Add the derivative of the panel's failure w.r.t it's DVs
void TACSBladeStiffenedShellConstitutive::addPanelFailureDVSens(
    const TacsScalar strain[], const TacsScalar scale, TacsScalar dfdx[]) {
  // In order to add values directly to the dfdx array, we need to first compute
  // the failure values for each ply, then compute the sensitivity of the
  // aggregated failure value with respect to each ply's value. Then we can
  // compute the sensitivity of each ply failure value with respect to the panel
  // thickness and add the weighted sensitivity to the dfdx array.

  // TODO: This is currently giving the wrong results

  if (this->panelThickNum >= 0) {
    TACSOrthotropicPly* ply = this->panelPly;
    const int numPlies = this->numPanelPlies;
    const int numStrain = TACSBeamConstitutive::NUM_STRESSES;
    TacsScalar* dKSdFail = this->panelPlyFailDVSens;
    TacsScalar* fails = this->panelPlyFailValues;
    const TacsScalar* angles = this->panelPlyAngles;
    const TacsScalar t = this->panelThick;
    const int dvInd = this->panelThickLocalNum;

    // Compute the strain state at the top of the panel
    TacsScalar plyStrain[3];
    plyStrain[0] = strain[0] + 0.5 * t * strain[3];
    plyStrain[1] = strain[1] + 0.5 * t * strain[4];
    plyStrain[2] = strain[2] + 0.5 * t * strain[5];

    // Compute the failure criteria at this strain state for each ply angle
    for (int ii = 0; ii < numPlies; ii++) {
      fails[ii] = ply->failure(angles[ii], plyStrain);
    }

    // Now repeat for the bottom of the panel
    plyStrain[0] = strain[0] - 0.5 * t * strain[3];
    plyStrain[1] = strain[1] - 0.5 * t * strain[4];
    plyStrain[2] = strain[2] - 0.5 * t * strain[5];

    for (int ii = 0; ii < numPlies; ii++) {
      fails[numPlies + ii] = ply->failure(angles[ii], plyStrain);
    }

    // Compute the sensitivity of the aggregated failure w.r.t each ply failure
    ksAggregationSens(fails, 2 * numPlies, this->ksWeight, dKSdFail);

    // Compute the the sensitivity of the failure criteria w.r.t the strains
    plyStrain[0] = strain[0] + 0.5 * t * strain[3];
    plyStrain[1] = strain[1] + 0.5 * t * strain[4];
    plyStrain[2] = strain[2] + 0.5 * t * strain[5];
    for (int ii = 0; ii < numPlies; ii++) {
      TacsScalar plyFailStrainSens[3];
      ply->failureStrainSens(angles[ii], plyStrain, plyFailStrainSens);
      // Convert the sensitivity w.r.t the strain to the sensitivity w.r.t
      // panel thickness
      dfdx[dvInd] +=
          scale * dKSdFail[ii] * 0.5 *
          (strain[3] * plyFailStrainSens[0] + strain[4] * plyFailStrainSens[1] +
           strain[5] * plyFailStrainSens[2]);
    }

    // Now repeat for the bottom of the panel
    plyStrain[0] = strain[0] - 0.5 * t * strain[3];
    plyStrain[1] = strain[1] - 0.5 * t * strain[4];
    plyStrain[2] = strain[2] - 0.5 * t * strain[5];

    for (int ii = 0; ii < numPlies; ii++) {
      TacsScalar plyFailStrainSens[3];
      ply->failureStrainSens(angles[ii], plyStrain, plyFailStrainSens);

      dfdx[dvInd] -=
          scale * dKSdFail[numPlies + ii] * 0.5 *
          (strain[3] * plyFailStrainSens[0] + strain[4] * plyFailStrainSens[1] +
           strain[5] * plyFailStrainSens[2]);
    }
  }
}

// ==============================================================================
// Helper functions for computing the stiffner's strain/stress/stiffness
// ==============================================================================
// In future, these methods should be replaced by calls to a beam
// constitutive model

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
  TacsScalar tipStrain[3];
  memset(tipStrain, 0, 3 * sizeof(TacsScalar));
  tipStrain[0] = stiffenerStrain[0] + zTipOffset * stiffenerStrain[2];

  // Compute the failure criteria at this strain state for each ply angle
  for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
    this->stiffenerPlyFailValues[ii] =
        ply->failure(this->stiffenerPlyAngles[ii], tipStrain);
  }

  // Returned the aggregated value over all plies
  return ksAggregation(this->stiffenerPlyFailValues, this->numStiffenerPlies,
                       this->ksWeight);
}

TacsScalar TACSBladeStiffenedShellConstitutive::evalStiffenerFailureStrainSens(
    const TacsScalar strain[], TacsScalar sens[]) {
  TACSOrthotropicPly* ply = this->stiffenerPly;
  const int numPlies = this->numStiffenerPlies;
  const int numStrain = TACSBeamConstitutive::NUM_STRESSES;
  TacsScalar** dFaildStrain = this->stiffenerPlyFailStrainSens;
  TacsScalar* fails = this->stiffenerPlyFailValues;
  const TacsScalar* angles = this->stiffenerPlyAngles;

  // Compute the strain state at the tip of the stiffener
  TacsScalar zTipOffset = -(this->stiffenerHeight + this->stiffenerThick) -
                          this->computeStiffenerCentroidHeight();
  TacsScalar tipStrain[3];
  memset(tipStrain, 0, 3 * sizeof(TacsScalar));
  tipStrain[0] = strain[0] + zTipOffset * strain[2];

  // Compute the failure criteria at this strain state for each ply angle and
  // the sensitivity of the failure criteria w.r.t the strains
  for (int ii = 0; ii < numPlies; ii++) {
    TacsScalar plyFailStrainSens[3];
    fails[ii] =
        ply->failureStrainSens(angles[ii], tipStrain, plyFailStrainSens);
    // Convert the sensitivity w.r.t the tip strain to the sensitivity w.r.t
    // beam strains
    memset(dFaildStrain[ii], 0, numStrain * sizeof(TacsScalar));
    dFaildStrain[ii][0] = plyFailStrainSens[0];
    dFaildStrain[ii][2] = zTipOffset * plyFailStrainSens[0];
  }

  TacsScalar fail = ksAggregationSensProduct(
      fails, numPlies, numStrain, this->ksWeight, dFaildStrain, sens);

  return fail;
}

void TACSBladeStiffenedShellConstitutive::addStiffenerFailureDVSens(
    const TacsScalar strain[], const TacsScalar scale, const int dvLen,
    TacsScalar dfdx[]) {
  // TODO: Implement this
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
