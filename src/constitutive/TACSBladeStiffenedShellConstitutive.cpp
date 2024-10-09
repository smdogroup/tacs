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

// Explicit definition of static constexpr member, for some reason if this is
// not included, TACS will throw a linker error complaining that
// DUMMY_FAIL_VALUE is an undefined symbol when compiled in complex mode.
constexpr TacsScalar TACSBladeStiffenedShellConstitutive::DUMMY_FAIL_VALUE;

const char* const TACSBladeStiffenedShellConstitutive::constName =
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

  this->numDesignVars = this->numGeneralDV = this->numPanelDV =
      this->numStiffenerDV = 0;

  // --- General DVs ---
  // --- Panel length values ---
  this->panelLength = _panelLength;
  this->panelLengthNum = _panelLengthNum;
  this->panelLengthLocalNum = -1;
  if (_panelLengthNum >= 0) {
    this->panelLengthLocalNum = this->numDesignVars;
    this->numDesignVars++;
    this->numGeneralDV++;
  }
  this->panelLengthLowerBound = 0.000;
  this->panelLengthUpperBound = 1e20;

  // --- Stiffener pitch values ---
  this->stiffenerPitch = _stiffenerPitch;
  this->stiffenerPitchNum = _stiffenerPitchNum;
  this->stiffenerPitchLocalNum = -1;
  if (_stiffenerPitchNum >= 0) {
    this->stiffenerPitchLocalNum = this->numDesignVars;
    this->numDesignVars++;
    this->numGeneralDV++;
  }
  this->stiffenerPitchLowerBound = 1e-3;
  this->stiffenerPitchUpperBound = 1e20;

  // --- Panel DVs ---
  // --- Panel thickness values ---
  this->panelDVStartNum = this->numDesignVars;
  this->panelThick = _panelThick;
  this->panelThickNum = _panelThickNum;
  this->panelThickLocalNum = -1;
  if (_panelThickNum >= 0) {
    this->panelThickLocalNum = this->numDesignVars;
    this->numDesignVars++;
    this->numPanelDV++;
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
    this->panelPlyFracLocalNums[ii] = -1;
    if (_panelPlyFracNums[ii] >= 0) {
      this->panelPlyFracLocalNums[ii] = this->numDesignVars;
      this->numDesignVars++;
      this->numPanelDV++;
    }
    this->panelPlyFracLowerBounds[ii] = 0.0;
    this->panelPlyFracUpperBounds[ii] = 1.0;
  }

  // --- Stiffener DVs ---
  // --- Stiffener height values ---
  this->stiffenerDVStartNum = this->numDesignVars;
  this->stiffenerHeight = _stiffenerHeight;
  this->stiffenerHeightNum = _stiffenerHeightNum;
  this->stiffenerHeightLocalNum = -1;
  if (_stiffenerHeightNum >= 0) {
    this->stiffenerHeightLocalNum = this->numDesignVars;
    this->numDesignVars++;
    this->numStiffenerDV++;
  }
  this->stiffenerHeightLowerBound = 1e-3;
  this->stiffenerHeightUpperBound = 1e20;

  // --- Stiffener thickness values ---
  this->stiffenerThick = _stiffenerThick;
  this->stiffenerThickNum = _stiffenerThickNum;
  this->stiffenerThickLocalNum = -1;
  if (_stiffenerThickNum >= 0) {
    this->stiffenerThickLocalNum = this->numDesignVars;
    this->numDesignVars++;
    this->numStiffenerDV++;
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
    this->stiffenerPlyFracLocalNums[ii] = -1;
    if (_stiffenerPlyFracNums[ii] >= 0) {
      this->stiffenerPlyFracLocalNums[ii] = this->numDesignVars;
      this->numDesignVars++;
      this->numStiffenerDV++;
    }
    this->stiffenerPlyFracLowerBounds[ii] = 0.0;
    this->stiffenerPlyFracUpperBounds[ii] = 1.0;
  }

  // --- Stiffener flange fraction ---
  this->flangeFraction = _flangeFraction;

  // --- panel and stiffener ply stiffness matrices ---
  // Since the ply angles are fixed, the Q and Abar matrices for each ply in the
  // panel and stiffener will remain constant, so we can pre-compute them here.
  this->panelQMats = new TacsScalar[_numPanelPlies * this->NUM_Q_ENTRIES];
  this->panelAbarMats = new TacsScalar[_numPanelPlies * this->NUM_ABAR_ENTRIES];

  TACSOrthotropicPly* ply = this->panelPly;
  TacsScalar* angles = this->panelPlyAngles;
  for (int plyNum = 0; plyNum < _numPanelPlies; plyNum++) {
    ply->calculateQbar(angles[plyNum], &panelQMats[plyNum * NUM_Q_ENTRIES]);
    ply->calculateAbar(angles[plyNum],
                       &panelAbarMats[plyNum * NUM_ABAR_ENTRIES]);
  }

  this->stiffenerQMats =
      new TacsScalar[_numStiffenerPlies * this->NUM_Q_ENTRIES];
  this->stiffenerAbarMats =
      new TacsScalar[_numStiffenerPlies * this->NUM_ABAR_ENTRIES];

  ply = this->stiffenerPly;
  angles = this->stiffenerPlyAngles;
  for (int plyNum = 0; plyNum < _numStiffenerPlies; plyNum++) {
    ply->calculateQbar(angles[plyNum], &stiffenerQMats[plyNum * NUM_Q_ENTRIES]);
    ply->calculateAbar(angles[plyNum],
                       &stiffenerAbarMats[plyNum * NUM_ABAR_ENTRIES]);
  }

  // --- Work arrays, these are created to avoid needing to allocate memory
  // inside compute methods like evalFailure and evalFailureStrainSens ---

  // Arrays for storing failure values, need values for each ply angle at the
  // top and bottom of the panel and at the tip of the stiffener
  this->panelPlyFailValues = new TacsScalar[2 * _numPanelPlies];
  this->stiffenerPlyFailValues = new TacsScalar[2 * _numStiffenerPlies];

  // Arrays for storing failure strain sensitivities
  this->panelPlyFailStrainSens = new TacsScalar*[2 * _numPanelPlies];
  this->stiffenerPlyFailStrainSens = new TacsScalar*[2 * _numStiffenerPlies];
  for (int ii = 0; ii < 2 * _numPanelPlies; ii++) {
    this->panelPlyFailStrainSens[ii] = new TacsScalar[this->NUM_STRESSES];
  }
  for (int ii = 0; ii < 2 * _numStiffenerPlies; ii++) {
    this->stiffenerPlyFailStrainSens[ii] =
        new TacsScalar[TACSBeamConstitutive::NUM_STRESSES];
  }

  // Arrays for storing ply failure sensitivities
  this->panelPlyFailSens = new TacsScalar[2 * this->numPanelPlies];
  this->stiffenerPlyFailSens = new TacsScalar[2 * this->numStiffenerPlies];
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

  delete[] this->panelPlyFailSens;
  this->panelPlyFailSens = nullptr;

  delete[] this->stiffenerPlyFailSens;
  this->stiffenerPlyFailSens = nullptr;

  delete[] this->panelQMats;
  this->panelQMats = nullptr;
  delete[] this->panelAbarMats;
  this->panelAbarMats = nullptr;

  delete[] this->stiffenerQMats;
  this->stiffenerQMats = nullptr;
  delete[] this->stiffenerAbarMats;
  this->stiffenerAbarMats = nullptr;
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
      -this->computeStiffenerCentroidHeight() - 0.5 * pThick;
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
      -this->computeStiffenerCentroidHeight() - 0.5 * pThick;
  TacsScalar stiffenerMOI = this->computeStiffenerMOI();

  TacsScalar dzdt, dzdh;
  this->computeStiffenerCentroidHeightSens(dzdt, dzdh);
  dzdt *= -1;
  dzdh *= -1;

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
        scale[1] * sPitchInv *
        (sThick * stiffenerDensity *
         (-kf * pThick - kf * sThick - pThick - 2.0 * sHeight - 2.0 * sThick) *
         0.5);

    // --- Second moment of area contribution ---
    // d/dh(MOI + 0.5*rho*A*z^2) = d/dh(MOI) + 0.5*rho*(dAdh*z^2 + 2*A*z*dzdh)
    dfdx[ii] +=
        scale[2] * sPitchInv *
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
        scale[1] * sPitchInv *
        (sHeight * stiffenerDensity *
         (-kf * pThick - 2.0 * kf * sThick - pThick - sHeight - 4.0 * sThick) /
         2.0);

    // --- Second moment of area contribution ---
    // d/dt(MOI + 0.5*rho*A*z^2) = d/dt(MOI) + 0.5*rho*(dAdt*z^2 + 2*A*z*dzdt)
    dfdx[ii] +=
        scale[2] * sPitchInv *
        (dMOIdt + 0.5 * stiffenerDensity * stiffenerOffset *
                      (stiffenerOffset * dAdt + 2.0 * stiffenerArea * dzdt));
  }

  // --- Panel thickness sensitivity ---
  if (this->panelThickLocalNum >= 0) {
    int ii = this->panelThickLocalNum;
    // Density contribution
    dfdx[ii] += scale[0] * panelDensity;
    // First moment of area contribution
    dfdx[ii] -= scale[1] * (stiffenerDensity * sHeight * sThick * (kf + 1.0) *
                            0.5 * sPitchInv);
    // Second moment of area contributions
    dfdx[ii] += scale[2] * panelDensity * 0.125 * pThick * pThick;
    dfdx[ii] -= scale[2] * sPitchInv * 0.5 * stiffenerDensity * 2.0 *
                stiffenerArea * stiffenerOffset * 0.5;
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
  // NOTE: The commented out code below is a theoretically nicer way to compute
  // the stress, but it makes it difficult to include the stiffener's
  // constribution to the drilling stress. So I've just settled for computing
  // the full stiffness matrix and multiplying by the strain.
  // ========================================================
  // // Compute the panel stresses
  // memset(s, 0, this->NUM_STRESSES * sizeof(TacsScalar));
  // this->computePanelStress(e, s);

  // // Compute the stiffener beam stresses then transform them back to shell
  // // stresses
  // TacsScalar stiffenerStrain[TACSBeamConstitutive::NUM_STRESSES],
  //     stiffenerStress[TACSBeamConstitutive::NUM_STRESSES];
  // this->transformStrain(e, stiffenerStrain);

  // this->computeStiffenerStress(stiffenerStrain, stiffenerStress);
  // this->addStiffenerStress(stiffenerStress, s);

  // Just compute the stiffness matrix and multiply by the strain
  TacsScalar C[NUM_TANGENT_STIFFNESS_ENTRIES];
  this->computeStiffness(C);
  TacsScalar* A = &C[0];
  TacsScalar* B = &C[6];
  TacsScalar* D = &C[12];
  TacsScalar* As = &C[18];
  TacsScalar drill = C[21];
  this->computeStress(A, B, D, As, drill, e, s);
}

// Add the derivative of the product of stress with a vector psi to dfdx
void TACSBladeStiffenedShellConstitutive::addStressDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], const TacsScalar psi[], int dvLen,
    TacsScalar dfdx[]) {
  TacsScalar pInv = 1.0 / this->stiffenerPitch;
  TacsScalar stiffScale = pInv * scale;

  TacsScalar A, dAdt, dAdh, dIzzdt, dIzzdh, dJdt, dJdh, dzdt, dzdh, E, G;
  A = this->computeStiffenerArea();
  this->computeStiffenerAreaSens(dAdt, dAdh);
  this->computeStiffenerIzz();
  this->computeStiffenerIzzSens(dIzzdt, dIzzdh);
  this->computeStiffenerJxxSens(dJdt, dJdh);
  this->computeEffectiveModulii(this->numStiffenerPlies, this->stiffenerQMats,
                                this->stiffenerPlyFracs, &E, &G);
  this->computeStiffenerCentroidHeightSens(dzdt, dzdh);
  dzdt *= -1;
  dzdh *= -1;

  // Sensitivity of the panel stress values to it's DVs (this has been proven
  // correct)
  this->addPanelStressDVSens(scale, strain, psi, &dfdx[this->panelDVStartNum]);

  // Transform the psi vector the same way we do for strains
  TacsScalar stiffenerPsi[TACSBeamConstitutive::NUM_STRESSES];
  this->transformStrain(psi, stiffenerPsi);

  TacsScalar stiffenerStrain[TACSBeamConstitutive::NUM_STRESSES],
      stiffenerStress[TACSBeamConstitutive::NUM_STRESSES];
  this->transformStrain(strain, stiffenerStrain);
  this->computeStiffenerStress(stiffenerStrain, stiffenerStress);

  // Add scale * 1/pitch * stiffenerStress^T * d/dx(Te) * stiffenerPsi
  this->addStrainTransformProductDVsens(stiffenerStress, psi, stiffScale, dfdx);

  // Add scale * 1/pitch * Cstiffener * stiffenerPsi * d/dx(Te) * strain
  TacsScalar psiStress[TACSBeamConstitutive::NUM_STRESSES];
  this->computeStiffenerStress(stiffenerPsi, psiStress);
  this->addStrainTransformProductDVsens(psiStress, strain, stiffScale, dfdx);

  // Add scale*1/pitch * stiffenerPsi^T * d/dx(stiffenerStress)
  this->addStiffenerStressDVSens(stiffScale, stiffenerStrain, stiffenerPsi,
                                 &dfdx[this->stiffenerDVStartNum]);

  // Add the direct dependence on the stiffener pitch (this has been proven
  // correct)
  if (this->stiffenerPitchLocalNum >= 0) {
    int index = this->stiffenerPitchLocalNum;
    TacsScalar panelStress[this->NUM_STRESSES];
    TacsScalar stress[this->NUM_STRESSES];

    this->computePanelStress(strain, panelStress);
    this->evalStress(elemIndex, pt, X, strain, stress);

    for (int jj = 0; jj < this->NUM_STRESSES; jj++) {
      dfdx[index] -= stiffScale * (stress[jj] - panelStress[jj]) * psi[jj];
    }
  }

  // Add terms related to the stiffener's contribution to the drilling stress
  // (the code below here has been proven correct)

  if (this->stiffenerThickNum >= 0) {
    dfdx[this->stiffenerThickLocalNum] +=
        scale * 0.5 * pInv * DRILLING_REGULARIZATION * this->kcorr * G * dAdt *
        psi[8] * strain[8];
  }
  if (this->stiffenerHeightNum >= 0) {
    dfdx[this->stiffenerHeightLocalNum] +=
        scale * 0.5 * pInv * DRILLING_REGULARIZATION * this->kcorr * G * dAdh *
        psi[8] * strain[8];
  }

  for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
    if (this->stiffenerPlyFracNums[ii] >= 0) {
      TacsScalar* Q = &(this->stiffenerQMats[ii * NUM_Q_ENTRIES]);
      TacsScalar dGdf = Q[5];
      dfdx[this->stiffenerPlyFracLocalNums[ii]] +=
          psi[8] * scale * 0.5 * pInv * DRILLING_REGULARIZATION * this->kcorr *
          A * dGdf * strain[8];
    }
  }
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
  TacsScalar fails[this->NUM_FAILURES];
  TacsScalar aggFail = this->computeFailureValues(e, fails);
  return aggFail;
}

// Compute the failure values for each failure mode of the stiffened panel
TacsScalar TACSBladeStiffenedShellConstitutive::computeFailureValues(
    const TacsScalar e[], TacsScalar fails[]) {
  // Initialize the failure values to some very large and negative that won't
  // contribute to the KS aggregate
  for (int ii = 0; ii < this->NUM_FAILURES; ii++) {
    fails[ii] = DUMMY_FAIL_VALUE;
  }
  TacsScalar stiffenerStrain[TACSBeamConstitutive::NUM_STRESSES];
  this->transformStrain(e, stiffenerStrain);

  // --- Panel material failure ---
  if (this->includePanelMaterialFailure) {
    fails[0] = this->computePanelFailure(e);
  }

  // --- Stiffener material failure ---
  if (this->includeStiffenerMaterialFailure) {
    fails[1] = this->computeStiffenerFailure(stiffenerStrain);
  }

  // --- Local panel buckling ---
  if (this->includeLocalBuckling) {
    fails[2] = this->evalLocalPanelBuckling(e);
  }

  // --- Global buckling ---
  if (this->includeGlobalBuckling) {
    fails[3] = this->evalGlobalPanelBuckling(e);
  }

  // --- Stiffener column buckling ---
  if (this->includeStiffenerColumnBuckling) {
    fails[4] = this->evalStiffenerColumnBuckling(stiffenerStrain);
  }

  // --- Stiffener crippling ---
  if (this->includeStiffenerCrippling) {
    fails[5] = this->evalStiffenerCrippling(stiffenerStrain);
  }

  TacsScalar ksFail = ksAggregation(fails, this->NUM_FAILURES, this->ksWeight);

  return ksFail;
}

// Evaluate the derivative of the failure criteria w.r.t. the strain
TacsScalar TACSBladeStiffenedShellConstitutive::evalFailureStrainSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar e[], TacsScalar sens[]) {
  memset(sens, 0, this->NUM_STRESSES * sizeof(TacsScalar));

  TacsScalar fails[this->NUM_FAILURES], dKSdf[this->NUM_FAILURES];

  // Initialize the failure values to some very large and negative that won't
  // contribute to the KS aggregate
  for (int ii = 0; ii < this->NUM_FAILURES; ii++) {
    fails[ii] = DUMMY_FAIL_VALUE;
  }

  // --- Material failure ---
  TacsScalar panelFailSens[this->NUM_STRESSES];
  memset(panelFailSens, 0, this->NUM_STRESSES * sizeof(TacsScalar));
  if (this->includePanelMaterialFailure) {
    // First compute the sensitivity of the panel failure value
    fails[0] = this->evalPanelFailureStrainSens(e, panelFailSens);
  }

  TacsScalar stiffenerStrainSens[TACSBeamConstitutive::NUM_STRESSES],
      stiffenerMatFailSens[this->NUM_STRESSES];
  memset(stiffenerStrainSens, 0,
         TACSBeamConstitutive::NUM_STRESSES * sizeof(TacsScalar));
  memset(stiffenerMatFailSens, 0, this->NUM_STRESSES * sizeof(TacsScalar));

  TacsScalar stiffenerStrain[TACSBeamConstitutive::NUM_STRESSES];
  this->transformStrain(e, stiffenerStrain);
  if (this->includeStiffenerMaterialFailure) {
    // And now for the stiffener failure value, first in terms of the beam
    // strains, and then transformed back to shell strains
    fails[1] = this->evalStiffenerFailureStrainSens(stiffenerStrain,
                                                    stiffenerStrainSens);
    this->transformStrainSens(stiffenerStrainSens, stiffenerMatFailSens);
  }

  // --- Local panel buckling ---
  TacsScalar localBucklingSens[this->NUM_STRESSES];
  memset(localBucklingSens, 0, this->NUM_STRESSES * sizeof(TacsScalar));
  if (this->includeLocalBuckling) {
    fails[2] = this->evalLocalPanelBucklingStrainSens(e, localBucklingSens);
  }

  // --- Global buckling ---
  TacsScalar globalBucklingSens[this->NUM_STRESSES];
  memset(globalBucklingSens, 0, this->NUM_STRESSES * sizeof(TacsScalar));
  if (this->includeGlobalBuckling) {
    fails[3] = this->evalGlobalPanelBucklingStrainSens(e, globalBucklingSens);
  }

  // --- Stiffener column buckling ---
  TacsScalar stiffenerBucklingSens[this->NUM_STRESSES];
  memset(stiffenerBucklingSens, 0, this->NUM_STRESSES * sizeof(TacsScalar));
  if (this->includeStiffenerColumnBuckling) {
    fails[4] = evalStiffenerColumnBucklingStrainSens(stiffenerStrain,
                                                     stiffenerStrainSens);
    this->transformStrainSens(stiffenerStrainSens, stiffenerBucklingSens);
  }

  // --- Stiffener crippling ---
  TacsScalar stiffenerCripplingSens[this->NUM_STRESSES];
  memset(stiffenerCripplingSens, 0, this->NUM_STRESSES * sizeof(TacsScalar));
  if (this->includeStiffenerCrippling) {
    fails[5] =
        evalStiffenerCripplingStrainSens(stiffenerStrain, stiffenerStrainSens);
    this->transformStrainSens(stiffenerStrainSens, stiffenerCripplingSens);
  }

  // Compute the sensitivity of the aggregate failure value to the individual
  // failure mode values
  TacsScalar fail =
      ksAggregationSens(fails, this->NUM_FAILURES, this->ksWeight, dKSdf);

  // Sum up the strain sens for each failure mode
  memset(sens, 0, this->NUM_STRESSES * sizeof(TacsScalar));
  for (int ii = 0; ii < this->NUM_STRESSES; ii++) {
    sens[ii] =
        dKSdf[0] * panelFailSens[ii] + dKSdf[1] * stiffenerMatFailSens[ii] +
        dKSdf[2] * localBucklingSens[ii] + dKSdf[3] * globalBucklingSens[ii] +
        dKSdf[4] * stiffenerBucklingSens[ii] +
        dKSdf[5] * stiffenerCripplingSens[ii];
  }

  return fail;
}

// Add the derivative of the failure criteria w.r.t. the design variables
void TACSBladeStiffenedShellConstitutive::addFailureDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], int dvLen, TacsScalar dfdx[]) {
  // TODO: The implementation of this function is not great, it uses a mix of
  // forward and backward differentiation and ends up recomputing a lot of
  // stuff. It should be rewritten to use only forward or only backward
  // differentiation in future.

  // Compute the failure values and then compute the
  // sensitivity of the aggregate failure value w.r.t. them
  TacsScalar fails[this->NUM_FAILURES], dKSdf[this->NUM_FAILURES];
  this->computeFailureValues(strain, fails);
  ksAggregationSens(fails, this->NUM_FAILURES, this->ksWeight, dKSdf);

  TacsScalar stiffenerStrain[TACSBeamConstitutive::NUM_STRESSES];
  this->transformStrain(strain, stiffenerStrain);

  // Sensitivity of the panel failure value to it's DVs
  if (this->includePanelMaterialFailure) {
    this->addPanelFailureDVSens(strain, scale * dKSdf[0],
                                &dfdx[this->panelDVStartNum]);
  }

  if (this->includeStiffenerMaterialFailure) {
    // Add the direct sensitivity of the stiffener failure value w.r.t DVs
    // Sensitivity of the panel failure value to it's DVs
    this->addStiffenerFailureDVSens(stiffenerStrain, scale * dKSdf[1],
                                    &dfdx[this->stiffenerDVStartNum]);

    // Add the sensitivity of the stiffener failure value w.r.t. the DVs
    // due to the dependence of the stiffener strains on the DVs
    TacsScalar stiffenerFailStrainSens[TACSBeamConstitutive::NUM_STRESSES];
    this->evalStiffenerFailureStrainSens(stiffenerStrain,
                                         stiffenerFailStrainSens);
    this->addStrainTransformProductDVsens(stiffenerFailStrainSens, strain,
                                          scale * dKSdf[1], dfdx);
  }

  // --- Local panel buckling ---
  if (this->includeLocalBuckling) {
    this->addLocalPanelBucklingDVSens(elemIndex, scale * dKSdf[2], pt, X,
                                      strain, dvLen, dfdx);
  }

  // --- Global buckling sensitivity ---
  if (this->includeGlobalBuckling) {
    this->addGlobalPanelBucklingDVSens(elemIndex, scale * dKSdf[3], pt, X,
                                       strain, dvLen, dfdx);
  }

  // --- Stiffener column buckling ---
  if (this->includeStiffenerColumnBuckling) {
    TacsScalar stiffenerStress[TACSBeamConstitutive::NUM_STRESSES];
    this->computeStiffenerStress(stiffenerStrain, stiffenerStress);
    const TacsScalar columnBucklingLoad =
        this->computeStiffenerColumnBucklingLoad();
    const TacsScalar stiffenerAxialLoad = -stiffenerStress[0];
    this->addStiffenerColumnBucklingDVSens(dKSdf[4] * scale, strain,
                                           stiffenerStrain, stiffenerAxialLoad,
                                           columnBucklingLoad, dfdx);
  }

  // --- Stiffener crippling ---
  if (this->includeStiffenerCrippling) {
    // Direct dependence of the stiffener crippling on the DVs
    this->addStiffenerCripplingDVSens(scale * dKSdf[5], stiffenerStrain,
                                      &dfdx[this->stiffenerDVStartNum]);
    // Sensitivity of the stiffener crippling due to effect of DVs on the
    // stiffener strains
    TacsScalar stiffenerCripplingStrainSens[TACSBeamConstitutive::NUM_STRESSES];
    this->evalStiffenerCripplingStrainSens(stiffenerStrain,
                                           stiffenerCripplingStrainSens);
    this->addStrainTransformProductDVsens(stiffenerCripplingStrainSens, strain,
                                          scale * dKSdf[5], dfdx);
  }
}

// ==============================================================================
// Compute output quantities
// ==============================================================================
// Retrieve the design variable for plotting purposes
TacsScalar TACSBladeStiffenedShellConstitutive::evalDesignFieldValue(
    int elemIndex, const double pt[], const TacsScalar X[], int index) {
  switch (index) {
    case 0:
      return this->computeEffectiveThickness();
    case 1:
      return this->computeEffectiveBendingThickness();
    case 2:
      return this->panelLength;
    case 3:
      return this->stiffenerPitch;
    case 4:
      return this->panelThick;
    case 5:
      return this->stiffenerHeight;
    case 6:
      return this->stiffenerThick;
  }
  return 0.0;
}

TacsScalar TACSBladeStiffenedShellConstitutive::evalFailureFieldValue(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], int failIndex) {
  if (failIndex == 0) {
    return this->evalFailure(elemIndex, pt, X, strain);
  } else if (failIndex >= 1 && failIndex <= this->NUM_FAILURES) {
    TacsScalar fails[this->NUM_FAILURES];
    computeFailureValues(strain, fails);
    TacsScalar failVal = fails[failIndex - 1];
    if (TacsRealPart(failVal) == TacsRealPart(this->DUMMY_FAIL_VALUE)) {
      return 0.0;
    } else {
      return failVal;
    }
  } else {
    return 0.0;
  }
}

TacsScalar
TACSBladeStiffenedShellConstitutive::computeEffectiveBendingThickness() {
  TacsScalar IStiff = this->computeStiffenerIzz();
  TacsScalar zStiff = -this->computeStiffenerCentroidHeight();
  TacsScalar t = this->panelThick;
  TacsScalar AStiff = this->computeStiffenerArea();
  TacsScalar Ieff = t * t * t / 12.0 +
                    (IStiff + AStiff * zStiff * zStiff) / this->stiffenerPitch;
  return cbrt(12.0 * TacsRealPart(Ieff));
}

// ==============================================================================
// Helper functions for computing the stiffness matrices
// ==============================================================================

// Compute the stiffness matrix
void TACSBladeStiffenedShellConstitutive::computeStiffness(TacsScalar C[]) {
  // --- Zero out the C matrix ---
  memset(C, 0, this->NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));

  // Add the panel's stiffness contributions
  this->computePanelStiffness(C);

  // Compute the stiffener's beam stiffness matrix ,then transform it to a shell
  // stiffness matrix and add it
  TacsScalar Cstiff[TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
  this->computeStiffenerStiffness(Cstiff);
  this->addStiffenerStiffness(Cstiff, C);
}

void TACSBladeStiffenedShellConstitutive::computeSmearedStiffness(
    const int numPlies, const TacsScalar* const QMats,
    const TacsScalar* const AbarMats, const TacsScalar plyFractions[],
    TacsScalar Q[], TacsScalar ABar[]) {
  // Zero out the Q and ABar matrices
  memset(Q, 0, this->NUM_Q_ENTRIES * sizeof(TacsScalar));
  memset(ABar, 0, this->NUM_ABAR_ENTRIES * sizeof(TacsScalar));

  // Q = sum_i (f_i * Q_i)
  for (int plyNum = 0; plyNum < numPlies; plyNum++) {
    for (int entry = 0; entry < this->NUM_Q_ENTRIES; entry++) {
      Q[entry] += plyFractions[plyNum] * QMats[NUM_Q_ENTRIES * plyNum + entry];
    }
    for (int entry = 0; entry < this->NUM_ABAR_ENTRIES; entry++) {
      ABar[entry] +=
          plyFractions[plyNum] * AbarMats[NUM_ABAR_ENTRIES * plyNum + entry];
    }
  }
}

void TACSBladeStiffenedShellConstitutive::computeSmearedStiffness(
    const int numPlies, const TacsScalar* const QMats,
    const TacsScalar plyFractions[], TacsScalar Q[]) {
  // Zero out the Q and ABar matrices
  memset(Q, 0, this->NUM_Q_ENTRIES * sizeof(TacsScalar));

  // Q = sum_i (f_i * Q_i)
  for (int plyNum = 0; plyNum < numPlies; plyNum++) {
    for (int entry = 0; entry < this->NUM_Q_ENTRIES; entry++) {
      Q[entry] += plyFractions[plyNum] * QMats[NUM_Q_ENTRIES * plyNum + entry];
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
      -this->computeStiffenerCentroidHeight() - 0.5 * this->panelThick;

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
  TacsScalar z =
      -this->computeStiffenerCentroidHeight() - 0.5 * this->panelThick;

  panelStrainSens[0] = stiffenerStrainSens[0];
  panelStrainSens[1] = 0.0;
  panelStrainSens[2] = 0.5 * stiffenerStrainSens[5];
  panelStrainSens[3] = stiffenerStrainSens[2] + z * stiffenerStrainSens[0];
  panelStrainSens[4] = 0.0;
  panelStrainSens[5] =
      0.5 * (z * stiffenerStrainSens[5] - stiffenerStrainSens[1]);
  panelStrainSens[6] = 0.0;
  panelStrainSens[7] = stiffenerStrainSens[4];
  panelStrainSens[8] = 0.0;
}

// Add the contribution of the stiffener stress to the panel stress
void TACSBladeStiffenedShellConstitutive::addStiffenerStress(
    const TacsScalar stiffenerStress[], TacsScalar panelStress[]) {
  TacsScalar pInv = 1.0 / this->stiffenerPitch;
  // Compute the offset of the stiffener centroid from the shell mid-plane
  TacsScalar z =
      -this->computeStiffenerCentroidHeight() - 0.5 * this->panelThick;

  panelStress[0] += pInv * stiffenerStress[0];
  panelStress[2] += pInv * 0.5 * stiffenerStress[5];
  panelStress[3] += pInv * stiffenerStress[2] + z * stiffenerStress[0];
  panelStress[5] += pInv * 0.5 * (z * stiffenerStress[5] - stiffenerStress[1]);
  panelStress[7] += pInv * stiffenerStress[4];
}

// Add the contribution of the stiffener stiffness to the panel stiffness
void TACSBladeStiffenedShellConstitutive::addStiffenerStiffness(
    const TacsScalar stiffenerStiffness[], TacsScalar panelStiffness[]) {
  TacsScalar pInv = 1.0 / this->stiffenerPitch;
  // Compute the offset of the stiffener centroid from the shell mid-plane
  TacsScalar z =
      -this->computeStiffenerCentroidHeight() - 0.5 * this->panelThick;

  // Some shorthand for the entries of the stiffness matrix
  TacsScalar* A = &(panelStiffness[0]);
  TacsScalar* B = &(panelStiffness[6]);
  TacsScalar* D = &(panelStiffness[12]);
  TacsScalar* As = &(panelStiffness[18]);
  TacsScalar* drill = &(panelStiffness[21]);

  // A:
  A[0] += pInv * (stiffenerStiffness[0]);
  A[2] += pInv * (stiffenerStiffness[5] / 2.0);
  A[5] += pInv * (stiffenerStiffness[20] / 4.0);

  // B:
  B[0] += pInv * (z * stiffenerStiffness[0] + stiffenerStiffness[2]);
  B[2] +=
      pInv * (z * stiffenerStiffness[5] / 2.0 - stiffenerStiffness[1] / 2.0);
  B[5] +=
      pInv * (z * stiffenerStiffness[20] / 4.0 - stiffenerStiffness[10] / 4.0);

  // D:
  D[0] += pInv * (z * (z * stiffenerStiffness[0] + stiffenerStiffness[2]) +
                  z * stiffenerStiffness[2] + stiffenerStiffness[11]);
  D[2] +=
      pInv * (z * (z * stiffenerStiffness[5] + stiffenerStiffness[14]) / 2.0 -
              z * stiffenerStiffness[1] / 2.0 - stiffenerStiffness[7] / 2.0);
  D[5] +=
      pInv * (z * (z * stiffenerStiffness[20] - stiffenerStiffness[10]) / 4.0 -
              z * stiffenerStiffness[10] / 4.0 + stiffenerStiffness[6] / 4.0);

  // As:
  As[2] += pInv * (stiffenerStiffness[18]);

  drill[0] += pInv * 0.5 * (stiffenerStiffness[18]) * DRILLING_REGULARIZATION;
}

void TACSBladeStiffenedShellConstitutive::addStrainTransformProductDVsens(
    const TacsScalar lhs[], const TacsScalar rhs[], const TacsScalar scale,
    TacsScalar dfdx[]) {
  // First compute the sensitivity of the stiffener centroid height w.r.t the
  // design variables (panel thickness, stiffener height and stiffener
  // thickness) zc = -panelThick/2 + computeStiffenerCentroidHeight()
  TacsScalar dzdtp, dzdhs, dzdts;
  dzdtp = -0.5;
  this->computeStiffenerCentroidHeightSens(dzdts, dzdhs);
  dzdts *= -1;
  dzdhs *= -1;

  // The sensitivities of the transformation matrix w.r.t the offset are:
  // dTe[0,3]/dz = 1
  // dTe[5,5]/dz = 1/2

  // Therefore:
  // df/dfx[i] = dTe[0,3]/dz * lhs[0] *rhs[3] * dz/dx[i] + dTe[5,5]/dz * lhs[5]
  // * rhs[5] * dz/dx[i]

  if (this->panelThickNum >= 0) {
    dfdx[this->panelThickLocalNum] +=
        scale * (lhs[0] * rhs[3] + lhs[5] * 0.5 * rhs[5]) * dzdtp;
  }
  if (this->stiffenerHeightNum >= 0) {
    dfdx[this->stiffenerHeightLocalNum] +=
        scale * (lhs[0] * rhs[3] + lhs[5] * 0.5 * rhs[5]) * dzdhs;
  }
  if (this->stiffenerThickNum >= 0) {
    dfdx[this->stiffenerThickLocalNum] +=
        scale * (lhs[0] * rhs[3] + lhs[5] * 0.5 * rhs[5]) * dzdts;
  }
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

// Add the derivative of the product of panel stresses with a vector psi to dfdx
void TACSBladeStiffenedShellConstitutive::addPanelStressDVSens(
    const TacsScalar scale, const TacsScalar strain[], const TacsScalar psi[],
    TacsScalar dfdx[]) {
  // The stress calculation is:
  // s = C * e
  // [ s[0:3] ] = [ A 0 0  0     ][ e[0:3] ]
  // [ s[3:6] ] = [ 0 D 0  0     ][ e[3:6] ]
  // [ s[6:8] ] = [ 0 0 As 0     ][ e[6:8] ]
  // [ s[8]   ] = [ 0 0 0  drill ][ e[8]   ]

  // Where: A = t * sum_i (plyFrac[i] * Q(theta[i]))
  //        D = t^3/12 * sum_i (plyFrac[i] * Q(theta[i]))
  //        As = t * sum_i (plyFrac[i] * Abar(theta[i])) * kcorr
  //        drill = DRILLING_REGULARIZATION * t/2 * (As[0,0] + As[1,1]))

  // Therefore, the derivative w.r.t the panel thickness is:
  // d/dt (psi^T * s) = psi^T * (C * e) =
  // psi[0:3]^T * [sum_i (plyFrac[i] * Q(theta[i]))] * e[0:3] +
  // psi[3:6]^T * [sum_i (plyFrac[i] * Q(theta[i]))] * e[3:6] +
  // psi[6:8]^T * [sum_i (plyFrac[i] * Abar(theta[i]))] * e[6:8] +
  // psi[8] * DRILLING_REGULARIZATION * t/2 * (As[0,0] + As[1,1])) * e[8]

  // And the derivative w.r.t the ply fractions is:
  // d/dplyFrac[i] (psi^T * s) = psi^T * (C * e) =
  // t * psi[0:3]^T * Q(theta[i]) * e[0:3] +
  // t^3/12 * psi[3:6]^T * Q(theta[i]) * e[3:6] +
  // t * psi[6:8]^T * Abar(theta[i]) * e[6:8] *kcorr +
  // t * psi[8] * DRILLING_REGULARIZATION * 1/2 * (As[0,0] + As[1,1])) *
  // e[8]*kcorr

  // --- Panel thickness sensitivity ---
  if (this->panelThickNum >= 0) {
    int index = this->panelThickLocalNum - this->panelDVStartNum;
    TacsScalar t24 = this->panelThick * this->panelThick / 4.0;
    TacsScalar AMatProd, DMatProd, AsMatProd, drillProd;

    TacsScalar QPanel[this->NUM_Q_ENTRIES];
    TacsScalar AbarPanel[this->NUM_ABAR_ENTRIES];
    this->computeSmearedStiffness(this->numPanelPlies, this->panelQMats,
                                  this->panelAbarMats, this->panelPlyFracs,
                                  QPanel, AbarPanel);

    AMatProd = psi[0] * (QPanel[0] * strain[0] + QPanel[1] * strain[1] +
                         QPanel[2] * strain[2]) +
               psi[1] * (QPanel[1] * strain[0] + QPanel[3] * strain[1] +
                         QPanel[4] * strain[2]) +
               psi[2] * (QPanel[2] * strain[0] + QPanel[4] * strain[1] +
                         QPanel[5] * strain[2]);

    DMatProd = psi[3] * (QPanel[0] * strain[3] + QPanel[1] * strain[4] +
                         QPanel[2] * strain[5]) +
               psi[4] * (QPanel[1] * strain[3] + QPanel[3] * strain[4] +
                         QPanel[4] * strain[5]) +
               psi[5] * (QPanel[2] * strain[3] + QPanel[4] * strain[4] +
                         QPanel[5] * strain[5]);

    AsMatProd =
        this->kcorr *
        (psi[6] * (AbarPanel[0] * strain[6] + AbarPanel[1] * strain[7]) +
         psi[7] * (AbarPanel[1] * strain[6] + AbarPanel[2] * strain[7]));

    drillProd = this->kcorr * psi[8] * DRILLING_REGULARIZATION *
                (0.5 * (AbarPanel[0] + AbarPanel[2])) * strain[8];
    dfdx[index] += scale * (AMatProd + t24 * DMatProd + AsMatProd + drillProd);
  }

  // --- Ply fraction sensitivity ---
  TacsScalar t = this->panelThick;
  TacsScalar t3 = t * t * t / 12.0;
  for (int ii = 0; ii < this->numPanelPlies; ii++) {
    if (this->panelPlyFracNums[ii] >= 0) {
      int index = this->panelPlyFracLocalNums[ii] - this->panelDVStartNum;
      TacsScalar* Q = &(this->panelQMats[NUM_Q_ENTRIES * ii]);
      TacsScalar* Abar = &(this->panelAbarMats[NUM_ABAR_ENTRIES * ii]);

      dfdx[index] +=
          scale * t *
          (psi[0] * (Q[0] * strain[0] + Q[1] * strain[1] + Q[2] * strain[2]) +
           psi[1] * (Q[1] * strain[0] + Q[3] * strain[1] + Q[4] * strain[2]) +
           psi[2] * (Q[2] * strain[0] + Q[4] * strain[1] + Q[5] * strain[2]));

      dfdx[index] +=
          scale * t3 *
          (psi[3] * (Q[0] * strain[3] + Q[1] * strain[4] + Q[2] * strain[5]) +
           psi[4] * (Q[1] * strain[3] + Q[3] * strain[4] + Q[4] * strain[5]) +
           psi[5] * (Q[2] * strain[3] + Q[4] * strain[4] + Q[5] * strain[5]));

      dfdx[index] += scale * t * this->kcorr *
                     (psi[6] * (Abar[0] * strain[6] + Abar[1] * strain[7]) +
                      psi[7] * (Abar[1] * strain[6] + Abar[2] * strain[7]));

      dfdx[index] += this->kcorr * scale * t * psi[8] *
                     DRILLING_REGULARIZATION * (0.5 * (Abar[0] + Abar[2])) *
                     strain[8];
    }
  }
}

void TACSBladeStiffenedShellConstitutive::computePanelStiffness(
    TacsScalar C[]) {
  TacsScalar* A = &C[0];
  TacsScalar* D = &C[12];
  TacsScalar* As = &C[18];

  // --- Zero out the C matrix ---
  memset(C, 0, this->NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));

  // Compute the smeared laminate properties
  TacsScalar QPanel[this->NUM_Q_ENTRIES], ABarPanel[this->NUM_ABAR_ENTRIES];

  this->computeSmearedStiffness(this->numPanelPlies, this->panelQMats,
                                this->panelAbarMats, this->panelPlyFracs,
                                QPanel, ABarPanel);

  // Add the panel's contributions to the A and D matrices
  TacsScalar t = this->panelThick;
  TacsScalar DFactor = t * t * t / 12.0;

  for (int ii = 0; ii < NUM_Q_ENTRIES; ii++) {
    A[ii] += t * QPanel[ii];
    D[ii] += DFactor * QPanel[ii];
  }

  // Add the pane;'s contribution to the transverse shear matrix
  for (int ii = 0; ii < NUM_ABAR_ENTRIES; ii++) {
    As[ii] += t * ABarPanel[ii] * this->kcorr;
  }

  // Add the drill stiffness
  C[21] = DRILLING_REGULARIZATION * 0.5 * (As[0] + As[2]);
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
  const int numStrain = TACSShellConstitutive::NUM_STRESSES;
  TacsScalar** dFaildStrain = this->panelPlyFailStrainSens;
  TacsScalar* fails = this->panelPlyFailValues;
  const TacsScalar* angles = this->panelPlyAngles;
  const TacsScalar t = this->panelThick;

  memset(sens, 0, numStrain * sizeof(TacsScalar));

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

  if (this->panelThickNum >= 0) {
    TACSOrthotropicPly* ply = this->panelPly;
    const int numPlies = this->numPanelPlies;
    TacsScalar* dKSdFail = this->panelPlyFailSens;
    TacsScalar* fails = this->panelPlyFailValues;
    const TacsScalar* angles = this->panelPlyAngles;
    const TacsScalar t = this->panelThick;

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
      dfdx[0] +=
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

      dfdx[0] -=
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
  TACSBeamConstitutive::computeStress(C, stiffenerStrain, stiffenerStress);
}

void TACSBladeStiffenedShellConstitutive::addStiffenerStressDVSens(
    const TacsScalar scale, const TacsScalar strain[], const TacsScalar psi[],
    TacsScalar dfdx[]) {
  TacsScalar psiStrain[TACSBeamConstitutive::NUM_STRESSES];
  for (int ii = 0; ii < TACSBeamConstitutive::NUM_STRESSES; ii++) {
    psiStrain[ii] = psi[ii] * strain[ii];
  }

  // Compute the cross-section properties and their derivatives
  TacsScalar A = this->computeStiffenerArea();
  TacsScalar Izz = this->computeStiffenerIzz();
  TacsScalar J = this->computeStiffenerJxx();
  TacsScalar dAdt, dAdh, dIzzdt, dIzzdh, dJdt, dJdh;
  this->computeStiffenerAreaSens(dAdt, dAdh);
  this->computeStiffenerIzzSens(dIzzdt, dIzzdh);
  this->computeStiffenerJxxSens(dJdt, dJdh);

  // Compute the beam laminate properties
  TacsScalar E, G, K;
  this->computeEffectiveModulii(this->numStiffenerPlies, this->stiffenerQMats,
                                this->stiffenerPlyFracs, &E, &G);
  K = this->kcorr;

  // --- Stiffener height sensitivities ---
  if (this->stiffenerHeightNum >= 0) {
    int index = this->stiffenerHeightLocalNum - this->stiffenerDVStartNum;
    dfdx[index] +=
        scale *
        ((E * psiStrain[0] + K * G * (psiStrain[4] + psiStrain[5])) * dAdh +
         (E * psiStrain[2]) * dIzzdh + (G * psiStrain[1]) * dJdh);
  }

  // --- Stiffener thickness sensitivities ---
  if (this->stiffenerThickNum >= 0) {
    int index = this->stiffenerThickLocalNum - this->stiffenerDVStartNum;
    dfdx[index] +=
        scale *
        ((E * psiStrain[0] + K * G * (psiStrain[4] + psiStrain[5])) * dAdt +
         (E * psiStrain[2]) * dIzzdt + (G * psiStrain[1]) * dJdt);
  }

  // --- Ply fraction sensitivities ---
  for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
    if (this->stiffenerPlyFracNums[ii] >= 0) {
      int index =
          this->stiffenerPlyFracLocalNums[ii] - this->stiffenerDVStartNum;

      TacsScalar* Q = &(this->stiffenerQMats[ii * NUM_Q_ENTRIES]);

      TacsScalar dEdx = Q[0] - (Q[1] * Q[1]) / Q[3];
      TacsScalar dGdx = Q[5];
      dfdx[index] += scale * (dEdx * (psiStrain[0] * A + psiStrain[2] * Izz) +
                              dGdx * (psiStrain[1] * J + psiStrain[4] * A * K +
                                      psiStrain[5] * A * K));
    }
  }
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

  TacsScalar E, G;
  this->computeEffectiveModulii(this->numStiffenerPlies, this->stiffenerQMats,
                                this->stiffenerPlyFracs, &E, &G);

  // Populate the matrix
  C[0] = E * A;                 // C[0, 0]
  C[6] = G * J;                 // C[1, 1]
  C[11] = E * Izz;              // C[2, 2]
  C[18] = this->kcorr * G * A;  // C[4, 4]
  C[20] = this->kcorr * G * A;  // C[5, 5]
}

void TACSBladeStiffenedShellConstitutive::computeEffectiveModulii(
    const int numPlies, const TacsScalar QMats[], const TacsScalar plyFracs[],
    TacsScalar* const E, TacsScalar* const G) {
  *E = 0.;
  *G = 0.;
  for (int plyNum = 0; plyNum < numPlies; plyNum++) {
    const TacsScalar* const Q = &(QMats[plyNum * NUM_Q_ENTRIES]);
    *E += plyFracs[plyNum] * (Q[0] - Q[1] * Q[1] / Q[3]);
    *G += plyFracs[plyNum] * Q[5];
  }

  // In theory the code below should produce exactly the same results as the
  // code above, but for some reason (probably related to floating point
  // arithmetic), it produces results that don't quite match complex step
  // derivatives w.r.t the ply fractions
  // TacsScalar Q[NUM_Q_ENTRIES];  //, ABar[this->NUM_ABAR_ENTRIES];
  // for (int ii = 0; ii < NUM_Q_ENTRIES; ii++) {
  //   Q[ii] = 0.0;
  // }
  // for (int plyNum = 0; plyNum < numPlies; plyNum++) {
  //   const TacsScalar* Qply = &(QMats[plyNum * NUM_Q_ENTRIES]);
  //   for (int ii = 0; ii < NUM_Q_ENTRIES; ii++) {
  //     Q[ii] += plyFracs[plyNum] * Qply[ii];
  //   }
  // }

  // // Compute the effective elastic and shear moduli
  // TacsScalar Q11 = Q[0];
  // TacsScalar Q12 = Q[1];
  // TacsScalar Q22 = Q[3];
  // TacsScalar Q66 = Q[5];
  // *E = Q11 - Q12 * Q12 / Q22;
  // *G = Q66;
}

// Compute the failure criteria for the stiffener
TacsScalar TACSBladeStiffenedShellConstitutive::computeStiffenerFailure(
    const TacsScalar stiffenerStrain[]) {
  TACSOrthotropicPly* ply = this->stiffenerPly;

  // Compute the strain state at the tip of the stiffener
  TacsScalar zCentroid = -this->computeStiffenerCentroidHeight();
  TacsScalar zTipOffset =
      -(this->stiffenerHeight + this->stiffenerThick) - zCentroid;
  TacsScalar plyStrain[3];
  memset(plyStrain, 0, 3 * sizeof(TacsScalar));
  plyStrain[0] = stiffenerStrain[0] + zTipOffset * stiffenerStrain[2];

  // Compute the failure criteria at this strain state for each ply angle
  for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
    this->stiffenerPlyFailValues[ii] =
        ply->failure(this->stiffenerPlyAngles[ii], plyStrain);
  }

  // Now do the same thing for the bottom of the stiffener
  zTipOffset = -zCentroid;
  plyStrain[0] = stiffenerStrain[0] + zTipOffset * stiffenerStrain[2];

  // Compute the failure criteria at this strain state for each ply angle
  for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
    this->stiffenerPlyFailValues[this->numStiffenerPlies + ii] =
        ply->failure(this->stiffenerPlyAngles[ii], plyStrain);
  }

  // Returned the aggregated value over all plies
  return ksAggregation(this->stiffenerPlyFailValues,
                       2 * this->numStiffenerPlies, this->ksWeight);
}

TacsScalar TACSBladeStiffenedShellConstitutive::evalStiffenerFailureStrainSens(
    const TacsScalar strain[], TacsScalar sens[]) {
  TACSOrthotropicPly* ply = this->stiffenerPly;
  const int numPlies = this->numStiffenerPlies;
  const int numStrain = TACSBeamConstitutive::NUM_STRESSES;
  TacsScalar** dFaildStrain = this->stiffenerPlyFailStrainSens;
  TacsScalar* fails = this->stiffenerPlyFailValues;
  const TacsScalar* angles = this->stiffenerPlyAngles;

  memset(sens, 0, numStrain * sizeof(TacsScalar));

  // Compute the strain state at the tip of the stiffener
  TacsScalar zCentroid = -this->computeStiffenerCentroidHeight();
  TacsScalar zTipOffset =
      -(this->stiffenerHeight + this->stiffenerThick) - zCentroid;
  TacsScalar tipStrain[3];
  memset(tipStrain, 0, 3 * sizeof(TacsScalar));
  tipStrain[0] = strain[0] + zTipOffset * strain[2];
  TacsScalar bottomStrain[3];
  memset(bottomStrain, 0, 3 * sizeof(TacsScalar));
  bottomStrain[0] = strain[0] - zCentroid * strain[2];

  // Compute the failure criteria at this strain state for each ply angle and
  // the sensitivity of the failure criteria w.r.t the strains
  for (int ii = 0; ii < numPlies; ii++) {
    TacsScalar plyFailStrainSens[3];
    fails[ii] =
        ply->failureStrainSens(angles[ii], tipStrain, plyFailStrainSens);
    // Convert the sensitivity w.r.t the tip strains to the sensitivity
    // w.r.t beam strains
    memset(dFaildStrain[ii], 0, numStrain * sizeof(TacsScalar));
    dFaildStrain[ii][0] = plyFailStrainSens[0];
    dFaildStrain[ii][2] = zTipOffset * plyFailStrainSens[0];

    // Now do the same for the bottom of the stiffener
    fails[numPlies + ii] =
        ply->failureStrainSens(angles[ii], bottomStrain, plyFailStrainSens);
    memset(dFaildStrain[numPlies + ii], 0, numStrain * sizeof(TacsScalar));
    dFaildStrain[numPlies + ii][0] = plyFailStrainSens[0];
    dFaildStrain[numPlies + ii][2] = -zCentroid * plyFailStrainSens[0];
  }

  TacsScalar fail = ksAggregationSensProduct(
      fails, 2 * numPlies, numStrain, this->ksWeight, dFaildStrain, sens);

  return fail;
}

void TACSBladeStiffenedShellConstitutive::addStiffenerFailureDVSens(
    const TacsScalar strain[], const TacsScalar scale, TacsScalar dfdx[]) {
  TACSOrthotropicPly* ply = this->stiffenerPly;
  TacsScalar* fails = this->stiffenerPlyFailValues;
  const TacsScalar* angles = this->stiffenerPlyAngles;
  TacsScalar* dKSdFail = this->stiffenerPlyFailSens;
  const int hNum = this->stiffenerHeightLocalNum - this->stiffenerDVStartNum;
  const int tNum = this->stiffenerThickLocalNum - this->stiffenerDVStartNum;

  TacsScalar zCentroid = -this->computeStiffenerCentroidHeight();
  TacsScalar zTipOffset =
      -(this->stiffenerHeight + this->stiffenerThick) - zCentroid;
  TacsScalar dZdh, dZdt;
  this->computeStiffenerCentroidHeightSens(dZdt, dZdh);
  dZdt *= -1;
  dZdh *= -1;
  TacsScalar dTipStraindh, dTipStraindt;
  dTipStraindh = -dZdh - 1.0;
  dTipStraindt = -dZdt - 1.0;

  // Compute the strain state at the tip of the stiffener
  TacsScalar tipStrain[3];
  memset(tipStrain, 0, 3 * sizeof(TacsScalar));
  tipStrain[0] = strain[0] + zTipOffset * strain[2];
  dTipStraindh *= strain[2];
  dTipStraindt *= strain[2];

  // Compute the failure criteria at this strain state for each ply angle
  for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
    fails[ii] = ply->failure(angles[ii], tipStrain);
  }

  // Now do the same thing for the bottom of the stiffener
  TacsScalar bottomStrain[3];
  memset(bottomStrain, 0, 3 * sizeof(TacsScalar));
  zTipOffset = -zCentroid;
  bottomStrain[0] = strain[0] + zTipOffset * strain[2];

  for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
    fails[this->numStiffenerPlies + ii] =
        ply->failure(this->stiffenerPlyAngles[ii], bottomStrain);
  }

  // Compute the sensitivity of the KS aggregation w.r.t the failure values
  ksAggregationSens(fails, 2 * this->numStiffenerPlies, this->ksWeight,
                    dKSdFail);

  // Now go back through each ply, compute the strain sensitivity of it's
  // failure, then convert it to a DV sensitivity and add it to the dfdx array
  for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
    TacsScalar tipFailStrainSens[3], bottomFailStrainSens[3];
    ply->failureStrainSens(angles[ii], tipStrain, tipFailStrainSens);
    ply->failureStrainSens(angles[ii], bottomStrain, bottomFailStrainSens);
    if (hNum >= 0) {
      dfdx[hNum] += scale * dKSdFail[ii] * dTipStraindh * tipFailStrainSens[0];
      dfdx[hNum] += scale * dKSdFail[this->numStiffenerPlies + ii] * -dZdh *
                    strain[2] * bottomFailStrainSens[0];
    }
    if (tNum >= 0) {
      dfdx[tNum] += scale * dKSdFail[ii] * dTipStraindt * tipFailStrainSens[0];
      dfdx[tNum] += scale * dKSdFail[this->numStiffenerPlies + ii] * -dZdt *
                    strain[2] * bottomFailStrainSens[0];
    }
  }
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
  return ((1.0 + 0.5 * this->flangeFraction) * this->stiffenerThick +
          0.5 * this->stiffenerHeight) /
         (1.0 + this->flangeFraction);
}

void TACSBladeStiffenedShellConstitutive::computeStiffenerCentroidHeightSens(
    TacsScalar& dzdt, TacsScalar& dzdh) {
  TacsScalar kf = this->flangeFraction;
  dzdh = 0.5 / (1.0 + kf);
  dzdt = (1.0 + 0.5 * kf) / (1.0 + kf);
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
  TacsScalar A1 = sh * st;           // Area of the stiffener web
  TacsScalar A2 = kf * A1;           // Area of the stiffener flange
  TacsScalar z1 = -(st + 0.5 * sh);  // Centroid of the stiffener web
  TacsScalar z2 = -(0.5 * st);       // Centroid of the stiffener flange
  TacsScalar zc =
      -this->computeStiffenerCentroidHeight();  // Centroid of the whole
                                                // stiffener section

  // Offsets of each area from the centroid of the whole stiffener section
  TacsScalar dz1 = z1 - zc;
  TacsScalar dz2 = z2 - zc;

  TacsScalar MOI1 = rho * A1 * sh * sh /
                    12.0;  // MOI of the stiffener web about it's own centroid
  TacsScalar MOI2 = rho * A2 * st * st / 12.0;  // MOI of the stiffener flange
                                                // about it's own centroid

  // Parallel axis theorem to get the MOI of the whole stiffener about its
  // centroid
  TacsScalar MOI = MOI1 + MOI2 + rho * (A1 * dz1 * dz1 + A2 * dz2 * dz2);
  return MOI;
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

  dMOIdt = 0.25 * sh * rho *
           (1.0 * kf2 * st2 + 4.0 / 3.0 * kf * sh2 + 4.0 * kf * sh * st +
            4.0 * kf * st2 + sh2 / 3.0) /
           (kf + 1.0);

  dMOIdh = 0.25 * st * rho *
           (kf2 / 3.0 * st2 + 4.0 * kf * sh2 + 4.0 * kf * sh * st +
            4.0 / 3.0 * kf * st2 + 1.0 * sh2) /
           (kf + 1.0);
}

// ==============================================================================
// Buckling functions
// ==============================================================================

TacsScalar TACSBladeStiffenedShellConstitutive::evalGlobalPanelBuckling(
    const TacsScalar e[]) {
  TacsScalar stress[TACSShellConstitutive::NUM_STRESSES];
  TacsScalar D1, D2, D3;
  this->computeCriticalGlobalBucklingStiffness(&D1, &D2, &D3);
  const TacsScalar L = this->panelLength;

  const TacsScalar N1Crit = computeCriticalGlobalAxialLoad(D1, L);
  const TacsScalar N12Crit = this->computeCriticalShearLoad(D1, D2, D3, L);

  this->evalStress(0, NULL, NULL, e, stress);
  return this->bucklingEnvelope(-stress[0], N1Crit, stress[2], N12Crit);
}

void TACSBladeStiffenedShellConstitutive::
    computeCriticalGlobalBucklingStiffness(TacsScalar* const D1,
                                           TacsScalar* const D2,
                                           TacsScalar* const D3) {
  // Compute the Q matrices for the panel and stiffener
  TacsScalar QPanel[this->NUM_Q_ENTRIES], QStiffener[this->NUM_Q_ENTRIES];
  this->computeSmearedStiffness(this->numPanelPlies, this->panelQMats,
                                this->panelPlyFracs, QPanel);
  this->computeSmearedStiffness(this->numStiffenerPlies, this->stiffenerQMats,
                                this->stiffenerPlyFracs, QStiffener);

  TacsScalar E1p, E1s, Ap, As, Ip, Is, Jp, Js, zs, tp, ps, tp3;
  tp = this->panelThick;
  tp3 = tp * tp * tp;
  ps = this->stiffenerPitch;
  Ap = tp * ps;
  As = this->computeStiffenerArea();
  Ip = ps * tp3 / 12.0;
  Is = this->computeStiffenerIzz();
  zs = -this->computeStiffenerCentroidHeight();
  Js = this->computeStiffenerJxx();
  Jp = ps * tp3 / 12.0;

  // Compute the effective modulii of the panel and stiffener
  E1p = QPanel[0] - QPanel[1] * QPanel[1] / QPanel[3];
  E1s = QStiffener[0] - QStiffener[1] * QStiffener[1] / QStiffener[3];

  // --- 1-Direction bending stiffness ---
  // Compute the modulus weighted centroid of the panel and stiffener
  TacsScalar zn = E1s * As * (-0.5 * tp + zs) / (E1s * As + E1p * Ap);

  *D1 = (E1p * (Ip + Ap * zn * zn) + E1s * (Is + As * (zn - zs) * (zn - zs))) /
        ps;

  // --- 2-direction bending stiffness ---
  // TacsScalar D2Panel = (tp3 * QPanel[3]) / 12.0;
  // TacsScalar D2Stiff = ((tp + 0.5 * ts) * (tp + 0.5 * ts) * (tp + 0.5 * ts)
  // *
  //                       (QPanel[3] + QStiffener[3])) /
  //                      6.0;
  // *D2 = ps / ((ps - hs * kf) / D2Panel + (hs * kf) / D2Stiff);
  // NOTE: I am ignoring the contribution of the stiffener flange to the
  // 2-direction bending stiffness so that the stiffness used for the buckling
  // calculation is consistent with the stiffness matrix. Not sure whether
  // this is a good idea or not.
  *D2 = tp3 / 12.0 * QPanel[3];

  // --- Twisting stiffness ---
  // Compute the shear modulus weighted centroid of the panel and stiffener
  TacsScalar zg = 0.25 * QStiffener[5] * As * (-0.5 * tp + zs) /
                  (QStiffener[5] * As + QPanel[5] * Ap);
  *D3 = (QPanel[5] * (Jp + Ap * zg * zg) +
         0.25 * QStiffener[5] * (Js + As * (zg - zs) * (zg - zs))) /
        ps;
}

TacsScalar
TACSBladeStiffenedShellConstitutive::evalGlobalPanelBucklingStrainSens(
    const TacsScalar e[], TacsScalar sens[]) {
  TacsScalar stiffness[NUM_TANGENT_STIFFNESS_ENTRIES], stress[NUM_STRESSES];
  this->computeStiffness(stiffness);
  const TacsScalar *A, *B, *D, *As;
  TacsScalar drill;
  this->extractTangentStiffness(stiffness, &A, &B, &D, &As, &drill);
  this->computeStress(A, B, D, As, drill, e, stress);
  TacsScalar N1GlobalSens, N1CritGlobalSens, N12GlobalSens, N12CritGlobalSens;
  TacsScalar D1, D2, D3;
  this->computeCriticalGlobalBucklingStiffness(&D1, &D2, &D3);
  const TacsScalar L = this->panelLength;
  TacsScalar N1CritGlobal = computeCriticalGlobalAxialLoad(D1, L);
  TacsScalar N12CritGlobal = this->computeCriticalShearLoad(D1, D2, D3, L);

  const TacsScalar strengthRatio = this->bucklingEnvelopeSens(
      -stress[0], N1CritGlobal, stress[2], N12CritGlobal, &N1GlobalSens,
      &N1CritGlobalSens, &N12GlobalSens, &N12CritGlobalSens);

  sens[0] += N1GlobalSens * -A[0] + N12GlobalSens * A[2];
  sens[1] += N1GlobalSens * -A[1] + N12GlobalSens * A[4];
  sens[2] += N1GlobalSens * -A[2] + N12GlobalSens * A[5];
  sens[3] += N1GlobalSens * -B[0] + N12GlobalSens * B[2];
  sens[4] += N1GlobalSens * -B[1] + N12GlobalSens * B[4];
  sens[5] += N1GlobalSens * -B[2] + N12GlobalSens * B[5];

  return strengthRatio;
}

void TACSBladeStiffenedShellConstitutive::addGlobalPanelBucklingDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], int dvLen, TacsScalar dfdx[]) {
  TacsScalar stress[NUM_STRESSES];
  this->evalStress(0, NULL, NULL, strain, stress);
  TacsScalar dfdN1Global, dfdN12Global, dfdN1CritGlobal, dfdN12CritGlobal;
  TacsScalar D1, D2, D3;
  TacsScalar N1Crit, N12Crit;
  this->computeCriticalGlobalBucklingStiffness(&D1, &D2, &D3);
  const TacsScalar L = this->panelLength;
  N1Crit = computeCriticalGlobalAxialLoad(D1, L);
  N12Crit = this->computeCriticalShearLoad(D1, D2, D3, L);

  this->bucklingEnvelopeSens(-stress[0], N1Crit, stress[2], N12Crit,
                             &dfdN1Global, &dfdN1CritGlobal, &dfdN12Global,
                             &dfdN12CritGlobal);

  // Add the sensitivity of the buckling failure criteria due to the
  // dependence of the applied loads on the DVs
  TacsScalar dfdStress[NUM_STRESSES];
  memset(dfdStress, 0, NUM_STRESSES * sizeof(TacsScalar));
  dfdStress[0] = -dfdN1Global;
  dfdStress[2] = dfdN12Global;
  this->addStressDVSens(elemIndex, scale, pt, X, strain, dfdStress, dvLen,
                        dfdx);

  // Propogate the sensitivity of the buckling failure criteria w.r.t the
  // critical loads back to the DVs
  TacsScalar dfdD1, dfdD2, dfdD3, dfdPanelLength;
  this->computeCriticalShearLoadSens(D1, D2, D3, L, &dfdD1, &dfdD2, &dfdD3,
                                     &dfdPanelLength);
  dfdD1 *= dfdN12CritGlobal;
  dfdD2 *= dfdN12CritGlobal;
  dfdD3 *= dfdN12CritGlobal;
  dfdPanelLength *= dfdN12CritGlobal;
  const TacsScalar N1CritGlobal = M_PI * M_PI * D1 / (L * L);
  dfdD1 += dfdN1CritGlobal * N1CritGlobal / D1;
  dfdPanelLength += dfdN1CritGlobal * -2.0 * N1CritGlobal / L;

  if (this->panelLengthNum >= 0) {
    dfdx[this->panelLengthLocalNum] += scale * dfdPanelLength;
  }

  TacsScalar globalBucklingspSens, globalBucklingtpSens, globalBucklinghsSens,
      globalBucklingtsSens, globalBucklingQstiffSens[NUM_Q_ENTRIES],
      globalBucklingQpanelSens[NUM_Q_ENTRIES];
  this->computeCriticalGlobalBucklingStiffnessSens(
      dfdD1, dfdD2, dfdD3, &globalBucklingspSens, &globalBucklingtpSens,
      &globalBucklinghsSens, &globalBucklingtsSens, globalBucklingQpanelSens,
      globalBucklingQstiffSens);

  // --- Panel thickness sensitivity ---
  if (this->panelThickNum >= 0) {
    const int dvNum = this->panelThickLocalNum;
    dfdx[dvNum] += scale * globalBucklingtpSens;
  }

  // --- Stiffener height contribution ---
  if (this->stiffenerHeightNum >= 0) {
    const int dvNum = this->stiffenerHeightLocalNum;
    // --- Global buckling contribution ---
    dfdx[dvNum] += scale * globalBucklinghsSens;
  }

  // --- Stiffener thickness contribution ---
  if (this->stiffenerThickNum >= 0) {
    const int dvNum = this->stiffenerThickLocalNum;
    // --- Global buckling contribution ---
    dfdx[dvNum] += scale * globalBucklingtsSens;
  }

  // --- Panel Ply fraction sensitivities ---
  for (int ii = 0; ii < this->numPanelPlies; ii++) {
    if (this->panelPlyFracNums[ii] >= 0) {
      const int dvNum = this->panelPlyFracLocalNums[ii];
      const TacsScalar* const Q = &(this->panelQMats[ii * NUM_Q_ENTRIES]);

      // --- Global buckling contribution ---
      dfdx[dvNum] += scale * (globalBucklingQpanelSens[0] * Q[0] +
                              globalBucklingQpanelSens[1] * Q[1] +
                              globalBucklingQpanelSens[2] * Q[2] +
                              globalBucklingQpanelSens[3] * Q[3] +
                              globalBucklingQpanelSens[4] * Q[4] +
                              globalBucklingQpanelSens[5] * Q[5]);
    }
  }

  // --- Stiffener ply fraction contributions ---
  for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
    if (this->stiffenerPlyFracNums[ii] >= 0) {
      const int dvNum = this->stiffenerPlyFracLocalNums[ii];
      const TacsScalar* const Q = &(this->stiffenerQMats[ii * NUM_Q_ENTRIES]);

      // --- Global buckling contribution ---
      dfdx[dvNum] += scale * (globalBucklingQstiffSens[0] * Q[0] +
                              globalBucklingQstiffSens[1] * Q[1] +
                              globalBucklingQstiffSens[2] * Q[2] +
                              globalBucklingQstiffSens[3] * Q[3] +
                              globalBucklingQstiffSens[4] * Q[4] +
                              globalBucklingQstiffSens[5] * Q[5]);
    }
  }

  // --- Stiffener pitch sensitivity ---
  if (this->stiffenerPitchNum >= 0) {
    dfdx[this->stiffenerPitchLocalNum] += scale * globalBucklingspSens;
  }
}

void TACSBladeStiffenedShellConstitutive::
    computeCriticalGlobalBucklingStiffnessSens(
        const TacsScalar dfdD1, const TacsScalar dfdD2, const TacsScalar dfdD3,
        TacsScalar* const psSens, TacsScalar* const tpSens,
        TacsScalar* const hsSens, TacsScalar* const tsSens,
        TacsScalar QpanelSens[], TacsScalar QstiffSens[]) {
  // Zero the sensitivities
  *tpSens = 0.0;
  *psSens = 0.0;
  *hsSens = 0.0;
  *tsSens = 0.0;
  memset(QstiffSens, 0, NUM_Q_ENTRIES * sizeof(TacsScalar));
  memset(QpanelSens, 0, NUM_Q_ENTRIES * sizeof(TacsScalar));

  // Compute the Q matrices for the panel and stiffener
  TacsScalar QPanel[this->NUM_Q_ENTRIES], QStiffener[this->NUM_Q_ENTRIES];
  this->computeSmearedStiffness(this->numPanelPlies, this->panelQMats,
                                this->panelPlyFracs, QPanel);
  this->computeSmearedStiffness(this->numStiffenerPlies, this->stiffenerQMats,
                                this->stiffenerPlyFracs, QStiffener);

  TacsScalar E1p, E1s, Ap, As, Ip, Is, Jp, Js, zs, tp, ps, tp3;
  tp = this->panelThick;
  tp3 = tp * tp * tp;
  ps = this->stiffenerPitch;
  Ap = tp * ps;
  As = this->computeStiffenerArea();
  Ip = ps * tp3 / 12.0;
  Is = this->computeStiffenerIzz();
  zs = -this->computeStiffenerCentroidHeight();
  Js = this->computeStiffenerJxx();
  Jp = ps * tp3 / 12.0;

  // Compute the effective modulii of the panel and stiffener
  E1p = QPanel[0] - QPanel[1] * QPanel[1] / QPanel[3];
  E1s = QStiffener[0] - QStiffener[1] * QStiffener[1] / QStiffener[3];

  // --- 1-Direction bending stiffness ---

  // Compute the modulus weighted centroid of the panel and stiffener
  TacsScalar zn = E1s * As * (-0.5 * tp + zs) / (E1s * As + E1p * Ap);
  TacsScalar zn2 = zn * zn;
  TacsScalar zns2 = (zn - zs) * (zn - zs);
  TacsScalar pInv = 1.0 / ps;

  // D1 = (E1p * (Ip + Ap * zn * zn) + E1s * (Is + As * (zn - zs) * (zn -
  // zs))) / ps;

  TacsScalar E1pSens = dfdD1 * ((Ap * zn2 + Ip) / ps);
  TacsScalar IpSens = dfdD1 * (E1p / ps);
  TacsScalar ApSens = dfdD1 * (E1p * zn2 / ps);
  TacsScalar ZnSens =
      dfdD1 * (2.0 * (Ap * E1p * zn + As * E1s * (zn - zs)) / ps);
  TacsScalar E1sSens = dfdD1 * ((As * zns2 + Is) / ps);
  TacsScalar IsSens = dfdD1 * (E1s / ps);
  TacsScalar AsSens = dfdD1 * (E1s * zns2 / ps);
  TacsScalar ZsSens = dfdD1 * (2.0 * As * E1s * (-zn + zs) / ps);
  *psSens +=
      dfdD1 * ((-E1p * (Ap * zn2 + Ip) - E1s * (As * zns2 + Is)) * pInv * pInv);

  // --- 2-direction bending stiffness ---
  // D2 = tp3 / 12.0 * QPanel[3];
  *tpSens = dfdD2 * tp * tp / 4.0 * QPanel[3];
  QpanelSens[3] = dfdD2 * tp3 / 12.0;

  // --- Twisting stiffness ---
  // D3 = (Gp * (Jp + Ap * zg**2) + 1 / 4 * Gs * (Js + As * (zg - zs) ** 2)) /
  // ps Compute the shear modulus weighted centroid of the panel and stiffener
  TacsScalar Gp = QPanel[5];
  TacsScalar Gs = QStiffener[5];
  TacsScalar zg = 0.25 * Gs * As * (-0.5 * tp + zs) / (Gs * As + Gp * Ap);
  TacsScalar zg2 = zg * zg;
  TacsScalar zgs2 = (zg - zs) * (zg - zs);

  TacsScalar GpSens = dfdD3 * ((Ap * zg2 + Jp) / ps);
  TacsScalar JpSens = dfdD3 * (Gp / ps);
  TacsScalar zgSens =
      dfdD3 * ((2.0 * Ap * Gp * zg + 0.5 * As * Gs * (zg - zs)) / ps);
  TacsScalar GsSens = dfdD3 * (0.25 * (As * zgs2 + Js) / ps);
  TacsScalar JsSens = dfdD3 * (0.25 * Gs / ps);
  ApSens += dfdD3 * (Gp * zg2 / ps);
  AsSens += dfdD3 * (0.25 * Gs * zgs2 / ps);
  ZsSens += dfdD3 * (0.5 * As * Gs * (-zg + zs) / ps);
  *psSens += dfdD3 * ((-Gp * (Ap * zg2 + Jp) - 0.25 * Gs * (As * zgs2 + Js)) *
                      pInv * pInv);

  // Now we need to propogate the sensitivities w.r.t intermediate variables
  // back to the inputs we want :
  // zn = E1s * As * (-0.5 * tp + zs) / (E1s * As + E1p * Ap);
  TacsScalar AE2 = (Ap * E1p + As * E1s) * (Ap * E1p + As * E1s);
  E1sSens += ZnSens * (-Ap * As * E1p * (0.5 * tp - zs) / AE2);
  AsSens += ZnSens * (-Ap * E1p * E1s * (0.5 * tp - zs) / AE2);
  *tpSens += ZnSens * (-0.5 * As * E1s / (Ap * E1p + As * E1s));
  ZsSens += ZnSens * (As * E1s / (Ap * E1p + As * E1s));
  E1pSens += ZnSens * (Ap * As * E1s * (0.5 * tp - zs) / AE2);
  ApSens += ZnSens * (As * E1p * E1s * (0.5 * tp - zs) / AE2);

  // zg = 0.25 * QStiffener[5] * As * (-0.5 * tp + zs) / (QStiffener[5] * As +
  // QPanel[5] * Ap);
  TacsScalar AG2 = (Ap * Gp + As * Gs) * (Ap * Gp + As * Gs);
  GsSens += zgSens * (-0.25 * Ap * As * Gp * (0.5 * tp - zs) / AG2);
  AsSens += zgSens * (-0.25 * Ap * Gp * Gs * (0.5 * tp - zs) / AG2);
  *tpSens += zgSens * (-0.125 * As * Gs / (Ap * Gp + As * Gs));
  ZsSens += zgSens * (0.25 * As * Gs / (Ap * Gp + As * Gs));
  GpSens += zgSens * (0.25 * Ap * As * Gs * (0.5 * tp - zs) / AG2);
  ApSens += zgSens * (0.25 * As * Gp * Gs * (0.5 * tp - zs) / AG2);

  // E1pSens
  // E1p = QPanel[0] - QPanel[1] * QPanel[1] / QPanel[3];
  QpanelSens[0] += E1pSens;
  QpanelSens[1] += E1pSens * (-2.0 * QPanel[1] / QPanel[3]);
  QpanelSens[3] +=
      E1pSens * ((QPanel[1] * QPanel[1]) / (QPanel[3] * QPanel[3]));

  // E1sSens
  // E1s = QStiff[0] - QStiff[1] * QStiff[1] / QStiff[3];
  QstiffSens[0] += E1sSens;
  QstiffSens[1] += E1sSens * (-2.0 * QStiffener[1] / QStiffener[3]);
  QstiffSens[3] += E1sSens * ((QStiffener[1] * QStiffener[1]) /
                              (QStiffener[3] * QStiffener[3]));

  // GpSens
  // Gp = QPanel[5]
  QpanelSens[5] += GpSens;

  // GsSens
  // Gs = QStiffener[5]
  QstiffSens[5] += GsSens;

  // IpSens
  // Ip = ps * tp3 / 12.0;
  *tpSens += IpSens * ps * (tp * tp / 4.0);
  *psSens += IpSens * (tp3 / 12.0);

  // JpSens
  // Jp = ps * tp3 / 12.0
  *psSens += JpSens * tp3 / 12.0;
  *tpSens += JpSens * ps * tp * tp / 4.0;

  // ApSens
  // Ap = tp * ps
  *psSens += ApSens * tp;
  *tpSens += ApSens * ps;

  // AsSens
  TacsScalar dAdt, dAdh;
  this->computeStiffenerAreaSens(dAdt, dAdh);
  *tsSens += AsSens * dAdt;
  *hsSens += AsSens * dAdh;

  // Js
  TacsScalar dJdt, dJdh;
  this->computeStiffenerJxxSens(dJdt, dJdh);
  *tsSens += JsSens * dJdt;
  *hsSens += JsSens * dJdh;

  // IsSens
  TacsScalar dIdt, dIdh;
  this->computeStiffenerIzzSens(dIdt, dIdh);
  *tsSens += IsSens * dIdt;
  *hsSens += IsSens * dIdh;

  // ZsSens
  TacsScalar dZdt, dZdh;
  this->computeStiffenerCentroidHeightSens(dZdt, dZdh);
  dZdt *= -1;
  dZdh *= -1;
  *tsSens += ZsSens * dZdt;
  *hsSens += ZsSens * dZdh;
}

void TACSBladeStiffenedShellConstitutive::testGlobalBucklingStiffnessSens() {
#ifdef TACS_USE_COMPLEX
  const double eps = 1e-200;
#else
  const double eps = 1e-6;
#endif
  // Get the number of DVs
  int num_dvs = this->getDesignVars(0, 0, NULL);
  TacsScalar* DV0 = new TacsScalar[num_dvs];
  TacsScalar* DVPert = new TacsScalar[num_dvs];
  this->getDesignVars(0, num_dvs, DV0);
  for (int i = 0; i < num_dvs; i++) {
    DVPert[i] = DV0[i];
  }

  TacsScalar D1, D2, D3, tpSens[3], psSens[3], tsSens[3], hsSens[3],
      QpSens[NUM_Q_ENTRIES], QsSens[NUM_Q_ENTRIES];
  this->computeCriticalGlobalBucklingStiffness(&D1, &D2, &D3);
  this->computeCriticalGlobalBucklingStiffnessSens(1.0, 0.0, 0.0, &psSens[0],
                                                   &tpSens[0], &hsSens[0],
                                                   &tsSens[0], QpSens, QsSens);
  this->computeCriticalGlobalBucklingStiffnessSens(0.0, 1.0, 0.0, &psSens[1],
                                                   &tpSens[1], &hsSens[1],
                                                   &tsSens[1], QpSens, QsSens);
  this->computeCriticalGlobalBucklingStiffnessSens(0.0, 0.0, 1.0, &psSens[2],
                                                   &tpSens[2], &hsSens[2],
                                                   &tsSens[2], QpSens, QsSens);

  TacsScalar tpSensRelError[3], psSensRelError[3], tsSensRelError[3],
      hsSensRelError[3];

  // ==============================================================================
  // stiffener pitch sensititivity
  // ==============================================================================
  printf("\nStiffener pitch sensitivity: \n");
  printf(" Analytic: D1sens = % 011.7e, D2sens = % 011.7e, D3sens = % 011.7e\n",
         TacsRealPart(psSens[0]), TacsRealPart(psSens[1]),
         TacsRealPart(psSens[2]));
  TacsScalar D1Pert, D2Pert, D3Pert;
#ifdef TACS_USE_COMPLEX
  DVPert[this->stiffenerPitchLocalNum] += TacsScalar(0.0, eps);
#else
  DVPert[this->stiffenerPitchLocalNum] += eps;
#endif
  this->setDesignVars(0, num_dvs, DVPert);
  this->computeCriticalGlobalBucklingStiffness(&D1Pert, &D2Pert, &D3Pert);

#ifdef TACS_USE_COMPLEX
  TacsScalar psSensFD = TacsImagPart(D1Pert) / eps;
#else
  TacsScalar psSensFD = (D1Pert - D1) / eps;
#endif
  if (psSens[0] == TacsScalar(0.)) {
    psSensRelError[0] = 0.;
  } else {
    psSensRelError[0] = (psSensFD - psSens[0]) / psSens[0];
  }
  printf("       FD: D1sens = % 011.7e, ", TacsRealPart(psSensFD));

#ifdef TACS_USE_COMPLEX
  psSensFD = TacsImagPart(D2Pert) / eps;
#else
  psSensFD = (D2Pert - D2) / eps;
#endif
  if (psSens[1] == TacsScalar(0.)) {
    psSensRelError[1] = 0.;
  } else {
    psSensRelError[1] = (psSensFD - psSens[1]) / psSens[1];
  }
  printf("D2sens = % 011.7e, ", TacsRealPart(psSensFD));

#ifdef TACS_USE_COMPLEX
  psSensFD = TacsImagPart(D3Pert) * 1e200;
#else
  psSensFD = (D3Pert - D3) / eps;
#endif
  if (psSens[2] == TacsScalar(0.)) {
    psSensRelError[2] = 0.;
  } else {
    psSensRelError[2] = (psSensFD - psSens[2]) / psSens[2];
  }
  printf("D3sens = % 011.7e\n", TacsRealPart(psSensFD));
  printf("Rel Error: D1sens = % 011.7e, D2sens = % 011.7e, D3sens = % 011.7e\n",
         TacsRealPart(psSensRelError[0]), TacsRealPart(psSensRelError[1]),
         TacsRealPart(psSensRelError[2]));

  DVPert[this->stiffenerPitchLocalNum] = DV0[this->stiffenerPitchLocalNum];
  this->setDesignVars(0, num_dvs, DVPert);

  // ==============================================================================
  // Panel thickness sensitivity
  // ==============================================================================
  printf("\nPanel thickness sensitivity: \n");
  printf(" Analytic: D1sens = % 011.7e, D2sens = % 011.7e, D3sens = % 011.7e\n",
         TacsRealPart(tpSens[0]), TacsRealPart(tpSens[1]),
         TacsRealPart(tpSens[2]));
#ifdef TACS_USE_COMPLEX
  DVPert[this->panelThickLocalNum] += TacsScalar(0.0, 1e-200);
#else
  DVPert[this->panelThickLocalNum] += eps;
#endif
  this->setDesignVars(0, num_dvs, DVPert);
  this->computeCriticalGlobalBucklingStiffness(&D1Pert, &D2Pert, &D3Pert);

#ifdef TACS_USE_COMPLEX
  TacsScalar tpSensFD = TacsImagPart(D1Pert) * 1e200;
#else
  TacsScalar tpSensFD = (D1Pert - D1) / eps;
#endif
  if (tpSens[0] == TacsScalar(0.)) {
    tpSensRelError[0] = 0.;
  } else {
    tpSensRelError[0] = (tpSensFD - tpSens[0]) / tpSens[0];
  }
  printf("       FD: D1sens = % 011.7e, ", TacsRealPart(tpSensFD));

#ifdef TACS_USE_COMPLEX
  tpSensFD = TacsImagPart(D2Pert) * 1e200;
#else
  tpSensFD = (D2Pert - D2) / eps;
#endif
  if (tpSens[1] == TacsScalar(0.)) {
    tpSensRelError[1] = 0.;
  } else {
    tpSensRelError[1] = (tpSensFD - tpSens[1]) / tpSens[1];
  }
  printf("D2sens = % 011.7e, ", TacsRealPart(tpSensFD));

#ifdef TACS_USE_COMPLEX
  tpSensFD = TacsImagPart(D3Pert) * 1e200;
#else
  tpSensFD = (D3Pert - D3) / eps;
#endif
  if (tpSens[2] == TacsScalar(0.)) {
    tpSensRelError[2] = 0.;
  } else {
    tpSensRelError[2] = (tpSensFD - tpSens[2]) / tpSens[2];
  }
  printf("D3sens = % 011.7e\n", TacsRealPart(tpSensFD));
  printf("Rel Error: D1sens = % 011.7e, D2sens = % 011.7e, D3sens = % 011.7e\n",
         TacsRealPart(tpSensRelError[0]), TacsRealPart(tpSensRelError[1]),
         TacsRealPart(tpSensRelError[2]));

  DVPert[this->panelThickLocalNum] = DV0[this->panelThickLocalNum];
  this->setDesignVars(0, num_dvs, DVPert);

  // ==============================================================================
  // Stiffener Height sensitivity
  // ==============================================================================
  printf("\nStiffener Height sensitivity: \n");
  printf(" Analytic: D1sens = % 011.7e, D2sens = % 011.7e, D3sens = % 011.7e\n",
         TacsRealPart(hsSens[0]), TacsRealPart(hsSens[1]),
         TacsRealPart(hsSens[2]));
#ifdef TACS_USE_COMPLEX
  DVPert[this->stiffenerHeightLocalNum] += TacsScalar(0.0, 1e-200);
#else
  DVPert[this->stiffenerHeightLocalNum] += eps;
#endif
  this->setDesignVars(0, num_dvs, DVPert);
  this->computeCriticalGlobalBucklingStiffness(&D1Pert, &D2Pert, &D3Pert);

#ifdef TACS_USE_COMPLEX
  TacsScalar hsSensFD = TacsImagPart(D1Pert) * 1e200;
#else
  TacsScalar hsSensFD = (D1Pert - D1) / eps;
#endif
  if (hsSens[0] == TacsScalar(0.)) {
    hsSensRelError[0] = 0.;
  } else {
    hsSensRelError[0] = (hsSensFD - hsSens[0]) / hsSens[0];
  }
  printf("       FD: D1sens = % 011.7e, ", TacsRealPart(hsSensFD));

#ifdef TACS_USE_COMPLEX
  hsSensFD = TacsImagPart(D2Pert) * 1e200;
#else
  hsSensFD = (D2Pert - D2) / eps;
#endif
  if (hsSens[1] == TacsScalar(0.)) {
    hsSensRelError[1] = 0.;
  } else {
    hsSensRelError[1] = (hsSensFD - hsSens[1]) / hsSens[1];
  }
  printf("D2sens = % 011.7e, ", TacsRealPart(hsSensFD));

#ifdef TACS_USE_COMPLEX
  hsSensFD = TacsImagPart(D3Pert) * 1e200;
#else
  hsSensFD = (D3Pert - D3) / eps;
#endif
  if (hsSens[2] == TacsScalar(0.)) {
    hsSensRelError[2] = 0.;
  } else {
    hsSensRelError[2] = (hsSensFD - hsSens[2]) / hsSens[2];
  }
  printf("D3sens = % 011.7e\n", TacsRealPart(hsSensFD));
  printf("Rel Error: D1sens = % 011.7e, D2sens = % 011.7e, D3sens = % 011.7e\n",
         TacsRealPart(hsSensRelError[0]), TacsRealPart(hsSensRelError[1]),
         TacsRealPart(hsSensRelError[2]));

  DVPert[this->stiffenerHeightLocalNum] = DV0[this->stiffenerHeightLocalNum];
  this->setDesignVars(0, num_dvs, DVPert);

  // ==============================================================================
  // Stiffener thickness sensitivity
  // ==============================================================================
  printf("\nStiffener thickness sensitivity: \n");
  printf(" Analytic: D1sens = % 011.7e, D2sens = % 011.7e, D3sens = % 011.7e\n",
         TacsRealPart(tsSens[0]), TacsRealPart(tsSens[1]),
         TacsRealPart(tsSens[2]));
#ifdef TACS_USE_COMPLEX
  DVPert[this->stiffenerThickLocalNum] += TacsScalar(0.0, 1e-200);
#else
  DVPert[this->stiffenerThickLocalNum] += eps;
#endif
  this->setDesignVars(0, num_dvs, DVPert);
  this->computeCriticalGlobalBucklingStiffness(&D1Pert, &D2Pert, &D3Pert);

#ifdef TACS_USE_COMPLEX
  TacsScalar tsSensFD = TacsImagPart(D1Pert) * 1e200;
#else
  TacsScalar tsSensFD = (D1Pert - D1) / eps;
#endif
  if (tsSens[0] == TacsScalar(0.)) {
    tsSensRelError[0] = 0.;
  } else {
    tsSensRelError[0] = (tsSensFD - tsSens[0]) / tsSens[0];
  }
  printf("       FD: D1sens = % 011.7e, ", TacsRealPart(tsSensFD));

#ifdef TACS_USE_COMPLEX
  tsSensFD = TacsImagPart(D2Pert) * 1e200;
#else
  tsSensFD = (D2Pert - D2) / eps;
#endif
  if (tsSens[1] == TacsScalar(0.)) {
    tsSensRelError[1] = 0.;
  } else {
    tsSensRelError[1] = (tsSensFD - tsSens[1]) / tsSens[1];
  }
  printf("D2sens = % 011.7e, ", TacsRealPart(tsSensFD));

#ifdef TACS_USE_COMPLEX
  tsSensFD = TacsImagPart(D3Pert) * 1e200;
#else
  tsSensFD = (D3Pert - D3) / eps;
#endif
  if (tsSens[2] == TacsScalar(0.)) {
    tsSensRelError[2] = 0.;
  } else {
    tsSensRelError[2] = (tsSensFD - tsSens[2]) / tsSens[2];
  }
  printf("D3sens = % 011.7e\n", TacsRealPart(tsSensFD));
  printf("Rel Error: D1sens = % 011.7e, D2sens = % 011.7e, D3sens = % 011.7e\n",
         TacsRealPart(tsSensRelError[0]), TacsRealPart(tsSensRelError[1]),
         TacsRealPart(tsSensRelError[2]));

  DVPert[this->stiffenerThickLocalNum] = DV0[this->stiffenerThickLocalNum];
  this->setDesignVars(0, num_dvs, DVPert);

  delete[] DV0;
  delete[] DVPert;
}

TacsScalar TACSBladeStiffenedShellConstitutive::evalLocalPanelBuckling(
    const TacsScalar e[]) {
  // Compute panel stiffness matrix and loads
  TacsScalar panelStiffness[NUM_TANGENT_STIFFNESS_ENTRIES],
      stress[NUM_STRESSES];
  this->computePanelStiffness(panelStiffness);
  const TacsScalar *A, *D;
  this->extractTangentStiffness(panelStiffness, &A, NULL, &D, NULL, NULL);
  this->computePanelStress(e, stress);

  // Compute the critical local loads
  const TacsScalar D11 = D[0], D12 = D[1], D22 = D[3], D66 = D[5],
                   L = this->stiffenerPitch;
  const TacsScalar N1Crit =
      this->computeCriticalLocalAxialLoad(D11, D22, D12, D66, L);
  const TacsScalar N12Crit =
      this->computeCriticalShearLoad(D11, D22, D12 + 2.0 * D66, L);

  // Compute the buckling criteria
  return this->bucklingEnvelope(-stress[0], N1Crit, stress[2], N12Crit);
}

TacsScalar TACSBladeStiffenedShellConstitutive::computeCriticalShearLoad(
    TacsScalar D1, TacsScalar D2, TacsScalar D3, TacsScalar L) {
  constexpr double ks = 50.0;
  const TacsScalar xi = sqrt(D1 * D2) / D3;

  // NOTE: sqrt(D3 * D1 * xi) = (D1^3 * D2)^0.25
  const TacsScalar N12_crit_1 =
      (4.0 / (L * L)) * sqrt(D3 * D1 * xi) * (8.125 + 5.045 / xi);
  const TacsScalar N12_crit_2 =
      (4.0 / (L * L)) * sqrt(D1 * D3) * (11.7 + 0.532 * xi + 0.938 * xi * xi);

  TacsScalar N12_min = 0.0;
  if (TacsRealPart(N12_crit_1) < TacsRealPart(N12_crit_2)) {
    N12_min = N12_crit_1;
  } else {
    N12_min = N12_crit_2;
  }

  const TacsScalar N12_diff = fabs(N12_crit_1 - N12_crit_2);
  const TacsScalar N12_crit = N12_min - log(1.0 + exp(-ks * N12_diff)) / ks;

  return N12_crit;
}

TacsScalar
TACSBladeStiffenedShellConstitutive::evalLocalPanelBucklingStrainSens(
    const TacsScalar e[], TacsScalar sens[]) {
  // Compute panel stiffness matrix and loads
  TacsScalar panelStiffness[NUM_TANGENT_STIFFNESS_ENTRIES],
      panelStress[NUM_STRESSES];
  this->computePanelStiffness(panelStiffness);
  const TacsScalar *APanel, *DPanel;
  this->extractTangentStiffness(panelStiffness, &APanel, NULL, &DPanel, NULL,
                                NULL);
  this->computePanelStress(e, panelStress);

  // Compute the critical local loads (no need to compute their
  // sensitivities because they're not dependent on the strain))
  TacsScalar N1CritLocal, N12CritLocal;
  TacsScalar D11 = DPanel[0], D12 = DPanel[1], D22 = DPanel[3], D66 = DPanel[5],
             L = this->stiffenerPitch;
  N1CritLocal = this->computeCriticalLocalAxialLoad(D11, D22, D12, D66, L);
  N12CritLocal = this->computeCriticalShearLoad(D11, D22, D12 + 2.0 * D66, L);

  // Compute the buckling criteria and it's sensitivities
  TacsScalar N1LocalSens, N12LocalSens, N1CritLocalSens, N12CritLocalSens;
  const TacsScalar strengthRatio = this->bucklingEnvelopeSens(
      -panelStress[0], N1CritLocal, panelStress[2], N12CritLocal, &N1LocalSens,
      &N1CritLocalSens, &N12LocalSens, &N12CritLocalSens);

  sens[0] = N1LocalSens * -APanel[0] + N12LocalSens * APanel[2];
  sens[1] = N1LocalSens * -APanel[1] + N12LocalSens * APanel[4];
  sens[2] = N1LocalSens * -APanel[2] + N12LocalSens * APanel[5];
  for (int ii = 3; ii < NUM_STRESSES; ii++) {
    sens[ii] = 0.0;
  }

  return strengthRatio;
}

void TACSBladeStiffenedShellConstitutive::addLocalPanelBucklingDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], int dvLen, TacsScalar dfdx[]) {
  // Compute panel stiffness matrix and loads
  TacsScalar panelStiffness[NUM_TANGENT_STIFFNESS_ENTRIES],
      panelStress[NUM_STRESSES];
  const TacsScalar t = this->panelThick;
  this->computePanelStiffness(panelStiffness);
  const TacsScalar *A, *D;
  this->extractTangentStiffness(panelStiffness, &A, NULL, &D, NULL, NULL);
  this->computePanelStress(strain, panelStress);

  // Compute the critical local loads and their sensitivities
  TacsScalar N1Crit, N12Crit;
  const TacsScalar D11 = D[0], D12 = D[1], D22 = D[3], D66 = D[5],
                   L = this->stiffenerPitch;

  // Create arrays for the sensitivities of the critical loads:
  // [dN/dD11, dNdD22, dNdD12, dNdD66, dNdL]
  TacsScalar N1CritSens[5], N12CritSens[5];
  N1Crit = this->computeCriticalLocalAxialLoadSens(
      D11, D22, D12, D66, L, &N1CritSens[0], &N1CritSens[1], &N1CritSens[2],
      &N1CritSens[3], &N1CritSens[4]);
  N12Crit = this->computeCriticalShearLoadSens(
      D11, D22, D12 + 2.0 * D66, L, &N12CritSens[0], &N12CritSens[1],
      &N12CritSens[2], &N12CritSens[4]);

  // N12CritSens[2] is currently dN12Crit/d(D12 + 2D66)
  N12CritSens[3] = 2.0 * N12CritSens[2];

  // Compute the buckling criteria and it's sensitivities to the applied and
  // critical loads
  TacsScalar dfdN1Local, dfdN12Local, dfdN1CritLocal, dfdN12CritLocal;
  this->bucklingEnvelopeSens(-panelStress[0], N1Crit, panelStress[2], N12Crit,
                             &dfdN1Local, &dfdN1CritLocal, &dfdN12Local,
                             &dfdN12CritLocal);
  // Convert sensitivity w.r.t applied loads into sensitivity w.r.t DVs
  TacsScalar dfdPanelStress[NUM_STRESSES];
  memset(dfdPanelStress, 0, NUM_STRESSES * sizeof(TacsScalar));
  dfdPanelStress[0] = -dfdN1Local;
  dfdPanelStress[2] = dfdN12Local;
  this->addPanelStressDVSens(scale, strain, dfdPanelStress,
                             &dfdx[this->panelDVStartNum]);

  // Convert the sensitivities of  the critical loads w.r.t the D matrix
  // entries to sensitivities of the buckling failure criteria w.r.t the D
  // matrix entries
  TacsScalar localBucklingDMatSens[4];
  for (int ii = 0; ii < 4; ii++) {
    localBucklingDMatSens[ii] =
        N1CritSens[ii] * dfdN1CritLocal + N12CritSens[ii] * dfdN12CritLocal;
  }

  // --- Panel thickness sensitivity ---
  if (this->panelThickNum >= 0) {
    const int dvNum = this->panelThickLocalNum;
    // --- Local buckling contribution ---
    const TacsScalar dMatdt = 0.25 * t * t;  // d/dt(t^3/12) = t^2/4
    for (int ii = 0; ii < this->numPanelPlies; ii++) {
      TacsScalar* Q = &(this->panelQMats[ii * NUM_Q_ENTRIES]);
      dfdx[dvNum] +=
          scale * dMatdt * this->panelPlyFracs[ii] *
          (localBucklingDMatSens[0] * Q[0] + localBucklingDMatSens[1] * Q[3] +
           localBucklingDMatSens[2] * Q[1] + localBucklingDMatSens[3] * Q[5]);
    }
  }

  // --- Panel Ply fraction sensitivities ---
  for (int ii = 0; ii < this->numPanelPlies; ii++) {
    if (this->panelPlyFracNums[ii] >= 0) {
      const int dvNum = this->panelPlyFracLocalNums[ii];
      // --- Local buckling contribution ---
      const TacsScalar DFactor = t * t * t / 12.0;
      const TacsScalar* const Q = &(this->panelQMats[ii * NUM_Q_ENTRIES]);
      // Do df/dx += df/dD11 * dD11/dx + df/dD22 * dD22/dx + df/dD12 * dD12/dx
      // + df/dD66 * dD66/dx
      dfdx[dvNum] +=
          scale * DFactor *
          (localBucklingDMatSens[0] * Q[0] + localBucklingDMatSens[1] * Q[3] +
           localBucklingDMatSens[2] * Q[1] + localBucklingDMatSens[3] * Q[5]);
    }
  }

  // --- Stiffener pitch sensitivity ---
  if (this->stiffenerPitchNum >= 0) {
    dfdx[this->stiffenerPitchLocalNum] +=
        scale *
        (N12CritSens[4] * dfdN12CritLocal + N1CritSens[4] * dfdN1CritLocal);
  }
}

TacsScalar
TACSBladeStiffenedShellConstitutive::computeCriticalLocalAxialLoadSens(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar D12,
    const TacsScalar D66, const TacsScalar L, TacsScalar* const D11Sens,
    TacsScalar* const D22Sens, TacsScalar* const D12Sens,
    TacsScalar* const D66Sens, TacsScalar* const LSens) {
  const double pi2 = M_PI * M_PI;
  const TacsScalar L2 = L * L;
  const TacsScalar root = sqrt(D11 * D22);

  *D11Sens = pi2 * root / (D11 * L2);
  *D12Sens = 2.0 * pi2 / L2;
  *D22Sens = pi2 * root / (D22 * L2);
  *D66Sens = 4.0 * pi2 / L2;
  *LSens = 4.0 * pi2 * (-D12 - 2.0 * D66 - root) / (L2 * L);

  return 2.0 * pi2 / L2 * (root + D12 + 2.0 * D66);
}

TacsScalar TACSBladeStiffenedShellConstitutive::computeCriticalShearLoadSens(
    const TacsScalar D1, const TacsScalar D2, const TacsScalar D3,
    const TacsScalar L, TacsScalar* const D1Sens, TacsScalar* const D2Sens,
    TacsScalar* const D3Sens, TacsScalar* const LSens) {
  const TacsScalar L2 = L * L;
  const TacsScalar L3 = L2 * L;
  const TacsScalar D32 = D3 * D3;
  const TacsScalar root = sqrt(D1 * D2);
  const TacsScalar xi = root / D3;

  TacsScalar N12CritVals[2];
  N12CritVals[0] = -(4.0 / L2) * sqrt(D3 * D1 * xi) * (8.125 + 5.045 / xi);
  N12CritVals[1] =
      -(4.0 / L2) * sqrt(D1 * D3) * (11.7 + 0.532 * xi + 0.938 * xi * xi);

  TacsScalar N12CritSens[2];
  const TacsScalar N12Crit =
      -ksAggregationSens(N12CritVals, 2, 50., N12CritSens);

  const TacsScalar dN12Crit1_dD1 =
      sqrt(D1 * root) * (5.045 * D3 + 24.375 * root) / (D1 * L2 * root);
  const TacsScalar dN12Crit2_dD1 =
      sqrt(D1 * D3) * (5.628 * D1 * D2 + 23.4 * D32 + 2.128 * D3 * root) /
      (D1 * D32 * L2);

  *D1Sens = dN12Crit1_dD1 * N12CritSens[0] + dN12Crit2_dD1 * N12CritSens[1];

  const TacsScalar dN12Crit1_dD2 =
      sqrt(D1 * root) * (-5.045 * D3 + 8.125 * root) / (D2 * L2 * root);
  const TacsScalar dN12Crit2_dD2 =
      sqrt(D1 * D3) * (3.752 * D1 * D2 + 1.064 * D3 * root) / (D2 * D32 * L2);

  *D2Sens = dN12Crit1_dD2 * N12CritSens[0] + dN12Crit2_dD2 * N12CritSens[1];

  const TacsScalar dN12Crit1_dD3 = 20.18 * sqrt(D1 * root) / (L2 * root);
  const TacsScalar dN12Crit2_dD3 =
      sqrt(D1 * D3) * (-5.628 * D1 * D2 + 23.4 * D32 - 1.064 * D3 * root) /
      (D32 * D3 * L2);

  *D3Sens = dN12Crit1_dD3 * N12CritSens[0] + dN12Crit2_dD3 * N12CritSens[1];

  const TacsScalar dN12Crit1_dL =
      sqrt(D1 * root) * (-40.36 * D3 - 65.0 * root) / (L3 * root);

  const TacsScalar dN12Crit2_dL =
      sqrt(D1 * D3) * (-7.504 * D1 * D2 - 93.6 * D32 - 4.256 * D3 * root) /
      (D32 * L3);

  *LSens = dN12Crit1_dL * N12CritSens[0] + dN12Crit2_dL * N12CritSens[1];

  return N12Crit;
}

// Test that the derivatives computed by computeCriticalShearLoadSens are
// correct by testing them against finite-difference values
bool TACSBladeStiffenedShellConstitutive::testCriticalShearLoadSens(
    const TacsScalar D1, const TacsScalar D2, const TacsScalar D3,
    const TacsScalar L) {
  TacsScalar D1Sens, D2Sens, D3Sens, LSens;
  TacsScalar N12Crit = computeCriticalShearLoadSens(D1, D2, D3, L, &D1Sens,
                                                    &D2Sens, &D3Sens, &LSens);

  TacsScalar eps = 1e-6;
  const double tol = 1e-5;
  TacsScalar D1SensFD, D2SensFD, D3SensFD, LSensFD;
  TacsScalar D1SensRelError, D2SensRelError, D3SensRelError, LSensRelError;
  TacsScalar N12CritFD = computeCriticalShearLoad(D1 + eps, D2, D3, L);
  D1SensFD = (N12CritFD - N12Crit) / eps;
  D1SensRelError = (D1Sens - D1SensFD) / (D1SensFD);

  N12CritFD = computeCriticalShearLoad(D1, D2 + eps, D3, L);
  D2SensFD = (N12CritFD - N12Crit) / eps;
  D2SensRelError = (D2Sens - D2SensFD) / (D2SensFD);

  N12CritFD = computeCriticalShearLoad(D1, D2, D3 + eps, L);
  D3SensFD = (N12CritFD - N12Crit) / eps;
  D3SensRelError = (D3Sens - D3SensFD) / (D3SensFD);

  N12CritFD = computeCriticalShearLoad(D1, D2, D3, L + eps);
  LSensFD = (N12CritFD - N12Crit) / eps;
  LSensRelError = (LSens - LSensFD) / (LSensFD);

  printf("testCriticalShearLoadSens results:\n");
  printf("----------------------------------\n");
  printf("D1 = % 011.7e\n", TacsRealPart(D1));
  printf("D2 = % 011.7e\n", TacsRealPart(D2));
  printf("D3 = % 011.7e\n", TacsRealPart(D3));
  printf("L = % 011.7e\n", TacsRealPart(L));
  printf("N12Crit = % 011.7e\n", TacsRealPart(N12Crit));
  printf(
      "D1Sens = % 011.7e | D1SensFD = % 011.7e | D1SensRelError = % 011.7e\n",
      TacsRealPart(D1Sens), TacsRealPart(D1SensFD),
      TacsRealPart(D1SensRelError));
  printf(
      "D2Sens = % 011.7e | D2SensFD = % 011.7e | D2SensRelError = % 011.7e\n",
      TacsRealPart(D2Sens), TacsRealPart(D2SensFD),
      TacsRealPart(D2SensRelError));
  printf(
      "D3Sens = % 011.7e | D3SensFD = % 011.7e | D3SensRelError = % 011.7e\n",
      TacsRealPart(D3Sens), TacsRealPart(D3SensFD),
      TacsRealPart(D3SensRelError));
  printf(
      " LSens = % 011.7e |  LSensFD = % 011.7e |  LSensRelError = % 011.7e\n",
      TacsRealPart(LSens), TacsRealPart(LSensFD), TacsRealPart(LSensRelError));
  printf(
      "----------------------------------------------------------------------"
      "--"
      "\n");

  return fabs(TacsRealPart(D1SensRelError)) < tol &&
         fabs(TacsRealPart(D2SensRelError)) < tol &&
         fabs(TacsRealPart(D3SensRelError)) < tol &&
         fabs(TacsRealPart(LSensRelError)) < tol;
}

TacsScalar TACSBladeStiffenedShellConstitutive::bucklingEnvelope(
    const TacsScalar N1, const TacsScalar N1Crit, const TacsScalar N12,
    const TacsScalar N12Crit) {
  const TacsScalar N1Term = N1 / N1Crit;
  const TacsScalar N12Term = N12 / N12Crit;
  const TacsScalar root = sqrt(N1Term * N1Term + 4.0 * N12Term * N12Term);
  return 0.5 * (N1Term + root);
}

// Compute the sensitivity of the buckling failure criterion w.r.t the loads
// and critical loads
TacsScalar TACSBladeStiffenedShellConstitutive::bucklingEnvelopeSens(
    const TacsScalar N1, const TacsScalar N1Crit, const TacsScalar N12,
    const TacsScalar N12Crit, TacsScalar* const dfdN1,
    TacsScalar* const dfdN1Crit, TacsScalar* const dfdN12,
    TacsScalar* const dfdN12Crit) {
  const TacsScalar N1Term = N1 / N1Crit;
  const TacsScalar N12Term = N12 / N12Crit;
  TacsScalar root = sqrt(N1Term * N1Term + 4.0 * N12Term * N12Term);
  if (TacsRealPart(root) == 0.0) {
    root = 1e-13;
  }
  if (dfdN1 != NULL) {
    // *dfdN1 = 1.0 / N1Crit;
    *dfdN1 = N1 / (2.0 * N1Crit * N1Crit * root) + 1.0 / (2.0 * N1Crit);
  }
  if (dfdN12 != NULL) {
    *dfdN12 = 2.0 * N12 / (N12Crit * N12Crit * root);
  }
  if (dfdN1Crit != NULL) {
    *dfdN1Crit = -N1 * N1 / (2.0 * N1Crit * N1Crit * N1Crit * root) -
                 N1 / (2.0 * N1Crit * N1Crit);
  }
  if (dfdN12Crit != NULL) {
    *dfdN12Crit = -2.0 * N12 * N12 / (N12Crit * N12Crit * N12Crit * root);
  }
  return bucklingEnvelope(N1, N1Crit, N12, N12Crit);
}

bool TACSBladeStiffenedShellConstitutive::testBucklingEnvelopeSens(
    const TacsScalar N1, const TacsScalar N1Crit, const TacsScalar N12,
    const TacsScalar N12Crit) {
  const double tol = 1e-5;
  TacsScalar dfdN1, dfdN1Crit, dfdN12, dfdN12Crit;
  double f = TacsRealPart(bucklingEnvelopeSens(
      N1, N1Crit, N12, N12Crit, &dfdN1, &dfdN1Crit, &dfdN12, &dfdN12Crit));

  double eps = 1e-6;
  TacsScalar N1p = N1 + eps;
  TacsScalar N1Critp = N1Crit + eps;
  TacsScalar N12p = N12 + eps;
  TacsScalar N12Critp = N12Crit + eps;

  double fp = TacsRealPart(bucklingEnvelope(N1p, N1Crit, N12, N12Crit));
  double dfdN1FD = TacsRealPart((fp - f) / eps);
  double dfdN1RelError = TacsRealPart((dfdN1 - dfdN1FD) / dfdN1FD);

  fp = TacsRealPart(bucklingEnvelope(N1, N1Critp, N12, N12Crit));
  double dfdN1CritFD = TacsRealPart((fp - f) / eps);
  double dfdN1CritRelError =
      TacsRealPart((dfdN1Crit - dfdN1CritFD) / dfdN1CritFD);

  fp = TacsRealPart(bucklingEnvelope(N1, N1Crit, N12p, N12Crit));
  double dfdN12FD = TacsRealPart((fp - f) / eps);
  double dfdN12RelError = TacsRealPart((dfdN12 - dfdN12FD) / dfdN12FD);

  fp = TacsRealPart(bucklingEnvelope(N1, N1Crit, N12, N12Critp));
  double dfdN12CritFD = TacsRealPart((fp - f) / eps);
  double dfdN12CritRelError =
      TacsRealPart((dfdN12Crit - dfdN12CritFD) / dfdN12CritFD);

  printf("testBucklingEnvelopeSens results:");
  printf(
      "N1 = % 011.7e, N1Crit = % 011.7e, N12 = % 011.7e, N12Crit = % "
      "011.7e\n",
      TacsRealPart(N1), TacsRealPart(N1Crit), TacsRealPart(N12),
      TacsRealPart(N12Crit));
  printf("f: % 011.7e\n", f);
  printf("dfdN1: % 011.7e, dfdN1p: % 011.7e, rel error: % 011.7e\n",
         TacsRealPart(dfdN1), TacsRealPart(dfdN1FD),
         TacsRealPart(dfdN1RelError));
  printf("dfdN1Crit: % 011.7e, dfdN1Critp: % 011.7e, rel error: % 011.7e\n",
         TacsRealPart(dfdN1Crit), TacsRealPart(dfdN1CritFD),
         TacsRealPart(dfdN1CritRelError));
  printf("dfdN12: % 011.7e, dfdN12p: % 011.7e, rel error: % 011.7e\n",
         TacsRealPart(dfdN12), TacsRealPart(dfdN12FD),
         TacsRealPart(dfdN1CritRelError));
  printf("dfdN12Crit: % 011.7e, dfdN12Critp: % 011.7e, rel error: % 011.7e\n",
         TacsRealPart(dfdN12Crit), TacsRealPart(dfdN12CritFD),
         TacsRealPart(dfdN12CritRelError));

  return (fabs(dfdN1RelError) < tol && fabs(dfdN1CritRelError) < tol &&
          fabs(dfdN12RelError) < tol && fabs(dfdN12CritRelError) < tol);
}

TacsScalar TACSBladeStiffenedShellConstitutive::evalStiffenerColumnBuckling(
    const TacsScalar stiffenerStrain[]) {
  TacsScalar stiffenerStress[TACSBeamConstitutive::NUM_STRESSES];
  this->computeStiffenerStress(stiffenerStrain, stiffenerStress);
  // The first component of the stiffener stress is the axial force
  TacsScalar fCrit = this->computeStiffenerColumnBucklingLoad();
  return -stiffenerStress[0] / fCrit;
}

TacsScalar
TACSBladeStiffenedShellConstitutive::computeStiffenerColumnBucklingLoad() {
  TacsScalar E, G;
  this->computeEffectiveModulii(this->numStiffenerPlies, this->stiffenerQMats,
                                this->stiffenerPlyFracs, &E, &G);
  const TacsScalar Izz = this->computeStiffenerIzz();
  const TacsScalar fCrit =
      this->computeColumnBucklingLoad(E, Izz, this->panelLength);
  return fCrit;
}

void TACSBladeStiffenedShellConstitutive::addStiffenerColumnBucklingLoadDVSens(
    const TacsScalar scale, TacsScalar dfdx[]) {
  // First compute the sensitivity of the critical load to the E, I and L
  // values of the stiffener
  TacsScalar E, G;
  this->computeEffectiveModulii(this->numStiffenerPlies, this->stiffenerQMats,
                                this->stiffenerPlyFracs, &E, &G);
  const TacsScalar Izz = this->computeStiffenerIzz();
  TacsScalar dFdE, dFdI, dFdL;
  this->computeColumnBucklingLoadSens(E, Izz, this->panelLength, dFdE, dFdI,
                                      dFdL);

  // Now convert those sensitivities into sensitivities with respect to the
  // DVs Beam length contributions (only relevant DV is panel length)
  if (this->panelLengthNum >= 0) {
    dfdx[this->panelLengthLocalNum] += scale * dFdL;
  }

  // Moment of inertia contributions (only relevant DVs are stiffener height
  // and thickness)
  TacsScalar dIzzdt, dIzzdh;
  this->computeStiffenerIzzSens(dIzzdt, dIzzdh);
  if (this->stiffenerHeightNum >= 0) {
    dfdx[this->stiffenerHeightLocalNum] += scale * dFdI * dIzzdh;
  }
  if (this->stiffenerThickNum >= 0) {
    dfdx[this->stiffenerThickLocalNum] += scale * dFdI * dIzzdt;
  }

  // Young's modulus contributions (only relevant DVs are stiffener ply
  // fractions)
  for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
    if (this->stiffenerPlyFracNums[ii] >= 0) {
      const int index = this->stiffenerPlyFracLocalNums[ii];

      const TacsScalar* const Q = &(this->stiffenerQMats[ii * NUM_Q_ENTRIES]);

      const TacsScalar dEdx = Q[0] - (Q[1] * Q[1]) / Q[3];
      dfdx[index] += scale * dFdE * dEdx;
    }
  }
}

TacsScalar TACSBladeStiffenedShellConstitutive::computeColumnBucklingLoad(
    const TacsScalar E, const TacsScalar I, const TacsScalar L) {
  return M_PI * M_PI * E * I / (L * L);
}

TacsScalar TACSBladeStiffenedShellConstitutive::computeColumnBucklingLoadSens(
    const TacsScalar E, const TacsScalar I, const TacsScalar L,
    TacsScalar& dFdE, TacsScalar& dFdI, TacsScalar& dFdL) {
  const double pi2 = M_PI * M_PI;
  const TacsScalar L2inv = 1.0 / (L * L);
  const TacsScalar F = pi2 * E * I * L2inv;
  dFdE = pi2 * I * L2inv;
  dFdI = pi2 * E * L2inv;
  dFdL = -2.0 * F / L;
  return F;
}

TacsScalar
TACSBladeStiffenedShellConstitutive::evalStiffenerColumnBucklingStrainSens(
    const TacsScalar stiffenerStrain[], TacsScalar sens[]) {
  const int numStiff = TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES;
  TacsScalar C[numStiff];
  this->computeStiffenerStiffness(C);
  const int numStress = TACSBeamConstitutive::NUM_STRESSES;
  TacsScalar stiffenerStress[numStress];
  TACSBeamConstitutive::computeStress(C, stiffenerStrain, stiffenerStress);

  const TacsScalar stiffenerAxialLoad = -stiffenerStress[0];
  const TacsScalar fCrit = this->computeStiffenerColumnBucklingLoad();
  memset(sens, 0, numStress * sizeof(TacsScalar));

  for (int ii = 0; ii < numStress; ii++) {
    sens[ii] = -C[ii] / fCrit;
  }
  return stiffenerAxialLoad / fCrit;
}

void TACSBladeStiffenedShellConstitutive::addStiffenerColumnBucklingDVSens(
    const TacsScalar scale, const TacsScalar shellStrain[],
    const TacsScalar stiffenerStrain[], const TacsScalar stiffenerAxialLoad,
    const TacsScalar fCrit, TacsScalar dfdx[]) {
  // F = stiffenerAxialLoad / fCrit
  // Using the Quotient rule:
  // dF/dx =  1/fCrit * dS11/dx - stiffenerAxialLoad/fCrit^2 * dfCrit/dx

  // First we add the 1/fCrit * dS11/dx term
  // We need to consider both the direct dependence of the stiffener stress on
  // the DVs, and also the dependence of the stiffener stress on the stiffener
  // centroid strains, which depend on the DVS:
  // psi^T dS/dx = psi^T pS/px + psi^T pS/pe * dstiffenerStrain/dx
  // = psi^T pS/px + C * psi * d/dx(Te) * shellStrain
  TacsScalar psi[TACSBeamConstitutive::NUM_STRESSES];
  memset(psi, 0, TACSBeamConstitutive::NUM_STRESSES * sizeof(TacsScalar));
  psi[0] = 1.0;
  this->addStiffenerStressDVSens(-scale / fCrit, stiffenerStrain, psi,
                                 &dfdx[this->stiffenerDVStartNum]);

  TacsScalar psiStress[TACSBeamConstitutive::NUM_STRESSES];
  this->computeStiffenerStress(psi, psiStress);
  this->addStrainTransformProductDVsens(psiStress, shellStrain, -scale / fCrit,
                                        dfdx);

  // Now we add the - stiffenerAxialLoad/fCrit^2 * dfCrit/dx term
  this->addStiffenerColumnBucklingLoadDVSens(
      -scale * stiffenerAxialLoad / (fCrit * fCrit), dfdx);
}

void TACSBladeStiffenedShellConstitutive::computeStiffenerCripplingValues(
    const TacsScalar stiffenerStrain[], TacsScalar plyFailValues[]) {
  // We use the semi-empirical formula for one-edge-free crippling described
  // in the DoD composite materials handbook Volume 3F (a.k.a MIL-HDBK-17-3F)
  // (available at
  // http://everyspec.com/MIL-HDBK/MIL-HDBK-0001-0099/MIL_HDBK_17_3F_216/)
  //
  // The formula is: F_crippling / F_ult = 1.63*(b/t)^-0.717
  // Where:
  // F_crippling is the crippling load
  // F_ult is the ultimate compressive load
  // b is the flange width
  // t is the flange thickness
  //
  // Here we compute the ultimate compressive load using the same criteria
  // used to compute the stiffener material failure value but using only the
  // axial strain component. This gives: SR_Crippling = SR_comp/1.63 *
  // (b/t)^0.717 Where SR_Crippling is the crippling strength ratio SR_comp is
  // the strength ratio computed using only the compressive component of the
  // strain

  const int numPlies = this->numStiffenerPlies;
  TacsScalar* const plyAngles = this->stiffenerPlyAngles;
  memset(plyFailValues, 0, 2 * numPlies * sizeof(TacsScalar));

  // We compute the crippling criteria for both the web and
  // the flange of the stiffener.
  TACSOrthotropicPly* ply = this->stiffenerPly;
  TacsScalar zCentroid = -this->computeStiffenerCentroidHeight();
  TacsScalar zFlange = -0.5 * this->stiffenerThick - zCentroid;

  TacsScalar flangeCrippleFactor = computeCripplingFactor(
      0.5 * this->flangeFraction * this->stiffenerHeight, this->stiffenerThick);

  TacsScalar plyStrain[3];
  memset(plyStrain, 0, 3 * sizeof(TacsScalar));
  plyStrain[0] = stiffenerStrain[0] + zFlange * stiffenerStrain[2];

  if (TacsRealPart(plyStrain[0]) < 0.0) {
    for (int ii = 0; ii < numPlies; ii++) {
      plyFailValues[ii] =
          ply->failure(plyAngles[ii], plyStrain) * flangeCrippleFactor;
    }
  }

  // For the stiffener web, we will use the method described by Kassapoglou to
  // account for the variation in tension/compression over the height of the
  // web. (See section 8.5.3 of Kassapoglou's "Design and Analysis of
  // Composite Structures with Applications to Aerospace Structures").
  const TacsScalar zWebTop =
      -(this->stiffenerHeight + this->stiffenerThick) - zCentroid;
  const TacsScalar zWebBottom = -this->stiffenerThick - zCentroid;
  const TacsScalar axStrainTop =
      stiffenerStrain[0] + zWebTop * stiffenerStrain[2];
  const TacsScalar axStrainBottom =
      stiffenerStrain[0] + zWebBottom * stiffenerStrain[2];

  // There are 3 cases we have to consider:
  // 1. Both strains are tensile
  // 2. Both strains are compressive
  // 3. One strain is tensile and the other is compressive
  //
  // In case 1 the strength ratio is simply zero as there can be no crippling.
  // In case 2 we compute the crippling strength ratio as usual using the
  // average of the two strains and the full height of the stiffener as the b
  // value
  // In cases 3  we compute the crippling strength ratio considering that
  // only the portion of the web that is in compression can cripplie, we
  // therefore use the length of that region as b and the average strain in
  // that region as the strain value.

  // In theory, the point at which the strain changes sign can be treated as a
  // simple support, therefore we should either treat the web as a
  // one-edge-free or no-edge-free flange depending whether it is the top or
  // bottom of the web that is in compression. However, this would lead to a
  // discontinuity in the strength ratio, so we use the one-edge-free formula
  // for both cases.
  if (TacsRealPart(axStrainTop) < 0.0 || TacsRealPart(axStrainBottom) < 0.0) {
    TacsScalar compressiveLength = 0.0, averageStrain = 0.0,
               webCrippleFactor = 0.0;
    if (TacsRealPart(axStrainTop) < 0 && TacsRealPart(axStrainBottom) < 0) {
      // --- Case 2 ---
      compressiveLength = this->stiffenerHeight;
      averageStrain = 0.5 * (axStrainTop + axStrainBottom);
    } else if (TacsRealPart(axStrainTop) < 0 &&
               TacsRealPart(axStrainBottom) >= 0) {
      // --- Case 3 top in compression ---
      compressiveLength =
          this->stiffenerHeight * axStrainTop / (axStrainTop - axStrainBottom);
      averageStrain = 0.5 * axStrainTop;
    } else if (TacsRealPart(axStrainTop) >= 0 &&
               TacsRealPart(axStrainBottom) < 0) {
      // --- Case 3 bottom in compression ---
      compressiveLength = this->stiffenerHeight * axStrainBottom /
                          (axStrainBottom - axStrainTop);
      averageStrain = 0.5 * axStrainBottom;
    }

    webCrippleFactor =
        computeCripplingFactor(compressiveLength, this->stiffenerThick);

    plyStrain[0] = averageStrain;
    for (int ii = 0; ii < numPlies; ii++) {
      plyFailValues[ii + numPlies] =
          ply->failure(plyAngles[ii], plyStrain) * webCrippleFactor;
    }
  }
}

TacsScalar TACSBladeStiffenedShellConstitutive::evalStiffenerCrippling(
    const TacsScalar stiffenerStrain[]) {
  const int numPlies = this->numStiffenerPlies;
  this->computeStiffenerCripplingValues(stiffenerStrain,
                                        this->stiffenerPlyFailValues);
  TacsScalar fail =
      ksAggregation(this->stiffenerPlyFailValues, 2 * numPlies, this->ksWeight);

  return fail;
}

TacsScalar
TACSBladeStiffenedShellConstitutive::evalStiffenerCripplingStrainSens(
    const TacsScalar stiffenerStrain[], TacsScalar sens[]) {
  const int numPlies = this->numStiffenerPlies;
  const int numStrain = TACSBeamConstitutive::NUM_STRESSES;
  TacsScalar* const plyFailValues = this->stiffenerPlyFailValues;
  TacsScalar** const dFaildStrain = this->stiffenerPlyFailStrainSens;
  TacsScalar* const plyAngles = this->stiffenerPlyAngles;
  memset(plyFailValues, 0, 2 * numPlies * sizeof(TacsScalar));
  memset(sens, 0, numStrain * sizeof(TacsScalar));

  for (int ii = 0; ii < 2 * numPlies; ii++) {
    memset(dFaildStrain[ii], 0, numStrain * sizeof(TacsScalar));
  }

  // We compute the crippling criteria for both the web and
  // the flange of the stiffener.

  // --- Flange crippling ---
  TACSOrthotropicPly* ply = this->stiffenerPly;
  TacsScalar zCentroid = -this->computeStiffenerCentroidHeight();
  TacsScalar zFlange = -0.5 * this->stiffenerThick - zCentroid;

  TacsScalar flangeCrippleFactor = computeCripplingFactor(
      0.5 * this->flangeFraction * this->stiffenerHeight, this->stiffenerThick);

  TacsScalar plyStrain[3];
  memset(plyStrain, 0, 3 * sizeof(TacsScalar));
  plyStrain[0] = stiffenerStrain[0] + zFlange * stiffenerStrain[2];

  if (TacsRealPart(plyStrain[0]) < 0.0) {
    for (int ii = 0; ii < numPlies; ii++) {
      TacsScalar plyFailStrainSens[3];
      plyFailValues[ii] =
          ply->failureStrainSens(plyAngles[ii], plyStrain, plyFailStrainSens) *
          flangeCrippleFactor;
      // Convert the sensitivity w.r.t the ply strains to the sensitivity
      // w.r.t beam strains
      dFaildStrain[ii][0] = flangeCrippleFactor * plyFailStrainSens[0];
      dFaildStrain[ii][2] =
          flangeCrippleFactor * zFlange * plyFailStrainSens[0];
    }
  }

  // --- Web crippling ---
  TacsScalar zWebTop =
      -(this->stiffenerHeight + this->stiffenerThick) - zCentroid;
  TacsScalar zWebBottom = -this->stiffenerThick - zCentroid;
  TacsScalar axStrainTop = stiffenerStrain[0] + zWebTop * stiffenerStrain[2];
  TacsScalar axStrainBottom =
      stiffenerStrain[0] + zWebBottom * stiffenerStrain[2];

  if (TacsRealPart(axStrainTop) < 0.0 || TacsRealPart(axStrainBottom) < 0.0) {
    TacsScalar compressiveLength = 0.0, averageStrain = 0.0,
               webCrippleFactor = 0.0;
    TacsScalar averageStrainSens[2] = {0.0, 0.0};
    TacsScalar lengthStrainSens[2] = {0.0, 0.0};
    if (TacsRealPart(axStrainTop) < 0 && TacsRealPart(axStrainBottom) < 0) {
      // --- Case 2 ---
      compressiveLength = this->stiffenerHeight;
      averageStrain = 0.5 * (axStrainTop + axStrainBottom);
      averageStrainSens[0] = 0.5;
      averageStrainSens[1] = 0.5;
    } else if (TacsRealPart(axStrainTop) < 0 &&
               TacsRealPart(axStrainBottom) >= 0) {
      // --- Case 3 top in compression ---
      const TacsScalar strainDiff = axStrainTop - axStrainBottom;
      compressiveLength = this->stiffenerHeight * axStrainTop / strainDiff;
      lengthStrainSens[0] =
          -this->stiffenerHeight * axStrainBottom / (strainDiff * strainDiff);
      lengthStrainSens[1] =
          this->stiffenerHeight * axStrainTop / (strainDiff * strainDiff);

      averageStrain = 0.5 * axStrainTop;
      averageStrainSens[0] = 0.5;
    } else if (TacsRealPart(axStrainTop) >= 0 &&
               TacsRealPart(axStrainBottom) < 0) {
      // --- Case 3 bottom in compression ---
      const TacsScalar strainDiff = axStrainBottom - axStrainTop;
      compressiveLength = this->stiffenerHeight * axStrainBottom / strainDiff;
      lengthStrainSens[0] =
          this->stiffenerHeight * axStrainBottom / (strainDiff * strainDiff);
      lengthStrainSens[1] =
          -this->stiffenerHeight * axStrainTop / (strainDiff * strainDiff);

      averageStrain = 0.5 * axStrainBottom;
      averageStrainSens[1] = 0.5;
    }

    TacsScalar dFactordL, dFactordt;
    webCrippleFactor = computeCripplingFactorSens(
        compressiveLength, this->stiffenerThick, dFactordL, dFactordt);

    plyStrain[0] = averageStrain;
    for (int ii = 0; ii < numPlies; ii++) {
      TacsScalar plyFailStrainSens[3];
      const TacsScalar strengthRatio =
          ply->failureStrainSens(plyAngles[ii], plyStrain, plyFailStrainSens);
      plyFailValues[ii + numPlies] = strengthRatio * webCrippleFactor;

      // Convert sensitivity w.r.t the average strain to sensitivity w.r.t
      // the top and bottom strains
      TacsScalar topStrainSens =
          strengthRatio * dFactordL * lengthStrainSens[0] +
          webCrippleFactor * plyFailStrainSens[0] * averageStrainSens[0];

      TacsScalar bottomStrainSens =
          strengthRatio * dFactordL * lengthStrainSens[1] +
          webCrippleFactor * plyFailStrainSens[0] * averageStrainSens[1];

      // Convert the sensitivity w.r.t the top and bottom strains to the
      // sensitivity w.r.t beam strains
      dFaildStrain[ii + numPlies][0] = topStrainSens + bottomStrainSens;
      dFaildStrain[ii + numPlies][2] =
          zWebTop * topStrainSens + zWebBottom * bottomStrainSens;
    }
  }

  TacsScalar fail =
      ksAggregationSensProduct(plyFailValues, 2 * numPlies, numStrain,
                               this->ksWeight, dFaildStrain, sens);

  return fail;
}

void TACSBladeStiffenedShellConstitutive::addStiffenerCripplingDVSens(
    const TacsScalar scale, const TacsScalar stiffenerStrain[],
    TacsScalar dfdx[]) {
  TACSOrthotropicPly* ply = this->stiffenerPly;
  TacsScalar* plyFailValues = this->stiffenerPlyFailValues;
  const TacsScalar* plyAngles = this->stiffenerPlyAngles;
  TacsScalar* dKSdFail = this->stiffenerPlyFailSens;
  const int hNum = this->stiffenerHeightLocalNum - this->stiffenerDVStartNum;
  const int tNum = this->stiffenerThickLocalNum - this->stiffenerDVStartNum;
  const int numPlies = this->numStiffenerPlies;
  const bool computeThicknessSens = (this->stiffenerThickNum >= 0);
  const bool computeHeightSens = (this->stiffenerHeightNum >= 0);

  this->computeStiffenerCripplingValues(stiffenerStrain, plyFailValues);

  ksAggregationSens(plyFailValues, 2 * numPlies, this->ksWeight, dKSdFail);

  // --- Flange crippling ---
  TacsScalar dzCentroiddh, dzCentroiddt;
  TacsScalar zCentroid = -this->computeStiffenerCentroidHeight();
  this->computeStiffenerCentroidHeightSens(dzCentroiddt, dzCentroiddh);
  dzCentroiddt *= -1;
  dzCentroiddh *= -1;
  TacsScalar zFlange = -0.5 * this->stiffenerThick - zCentroid;
  const TacsScalar dzFlangedt = -0.5 - dzCentroiddt;
  const TacsScalar dzFlangedh = -dzCentroiddh;

  TacsScalar dFactordh, dFactordt;
  TacsScalar flangeCrippleFactor = computeCripplingFactorSens(
      0.5 * this->flangeFraction * this->stiffenerHeight, this->stiffenerThick,
      dFactordh, dFactordt);
  dFactordh *= 0.5 * this->flangeFraction;

  TacsScalar plyStrain[3];
  memset(plyStrain, 0, 3 * sizeof(TacsScalar));
  plyStrain[0] = stiffenerStrain[0] + zFlange * stiffenerStrain[2];
  if (TacsRealPart(plyStrain[0]) < 0.0) {
    for (int ii = 0; ii < numPlies; ii++) {
      TacsScalar plyFailStrainSens[3];
      const TacsScalar strengthRatio =
          ply->failureStrainSens(plyAngles[ii], plyStrain, plyFailStrainSens);
      if (computeThicknessSens) {
        dfdx[tNum] += dKSdFail[ii] * scale *
                      (strengthRatio * dFactordt +
                       plyFailStrainSens[0] * dzFlangedt * stiffenerStrain[2] *
                           flangeCrippleFactor);
      }
      if (computeHeightSens) {
        dfdx[hNum] += dKSdFail[ii] * scale *
                      (strengthRatio * dFactordh +
                       plyFailStrainSens[0] * dzFlangedh * stiffenerStrain[2] *
                           flangeCrippleFactor);
      }
    }
  }

  // --- Flange crippling ---
  const TacsScalar zWebTop =
      -(this->stiffenerHeight + this->stiffenerThick) - zCentroid;
  const TacsScalar dzWebTopdh = -1.0 - dzCentroiddh;
  const TacsScalar dzWebTopdt = -1.0 - dzCentroiddt;

  const TacsScalar zWebBottom = -this->stiffenerThick - zCentroid;
  const TacsScalar dzWebBottomdt = -1.0 - dzCentroiddt;
  const TacsScalar dzWebBottomdh = -dzCentroiddh;

  const TacsScalar axStrainTop =
      stiffenerStrain[0] + zWebTop * stiffenerStrain[2];
  const TacsScalar dStrainTopdh = stiffenerStrain[2] * dzWebTopdh;
  const TacsScalar dStrainTopdt = stiffenerStrain[2] * dzWebTopdt;

  const TacsScalar axStrainBottom =
      stiffenerStrain[0] + zWebBottom * stiffenerStrain[2];
  const TacsScalar dStrainBottomdh = stiffenerStrain[2] * dzWebBottomdh;
  const TacsScalar dStrainBottomdt = stiffenerStrain[2] * dzWebBottomdt;

  if (TacsRealPart(axStrainTop) < 0.0 || TacsRealPart(axStrainBottom) < 0.0) {
    TacsScalar compressiveLength = 0.0, averageStrain = 0.0,
               webCrippleFactor = 0.0;
    TacsScalar dAverageStraindh = 0.0;
    TacsScalar dAverageStraindt = 0.0;
    TacsScalar dLdh = 0.0;
    TacsScalar dLdt = 0.0;
    if (TacsRealPart(axStrainTop) < 0 && TacsRealPart(axStrainBottom) < 0) {
      // --- Case 2 ---
      compressiveLength = this->stiffenerHeight;
      dLdh = 1.0;

      averageStrain = 0.5 * (axStrainTop + axStrainBottom);
      dAverageStraindh = 0.5 * (dStrainTopdh + dStrainBottomdh);
      dAverageStraindt = 0.5 * (dStrainTopdt + dStrainBottomdt);
    } else if (TacsRealPart(axStrainTop) < 0 &&
               TacsRealPart(axStrainBottom) >= 0) {
      // --- Case 3 top in compression ---
      const TacsScalar strainDiff = axStrainTop - axStrainBottom;
      compressiveLength = this->stiffenerHeight * axStrainTop / strainDiff;

      const TacsScalar dLdstrainTop =
          -this->stiffenerHeight * axStrainBottom / (strainDiff * strainDiff);
      const TacsScalar dLdstrainBottom =
          this->stiffenerHeight * axStrainTop / (strainDiff * strainDiff);

      // dLdx = pLpx + pLpet * detdx + pLpeb * debdx
      dLdh = axStrainTop / strainDiff + dLdstrainTop * dStrainTopdh +
             dLdstrainBottom * dStrainBottomdh;
      dLdt = dLdstrainTop * dStrainTopdt + dLdstrainBottom * dStrainBottomdt;

      averageStrain = 0.5 * axStrainTop;
      dAverageStraindh = 0.5 * dStrainTopdh;
      dAverageStraindt = 0.5 * dStrainTopdt;
    } else if (TacsRealPart(axStrainTop) >= 0 &&
               TacsRealPart(axStrainBottom) < 0) {
      // --- Case 3 bottom in compression ---
      const TacsScalar strainDiff = axStrainBottom - axStrainTop;

      compressiveLength = this->stiffenerHeight * axStrainBottom / strainDiff;

      const TacsScalar dLdstrainTop =
          this->stiffenerHeight * axStrainBottom / (strainDiff * strainDiff);
      const TacsScalar dLdstrainBottom =
          -this->stiffenerHeight * axStrainTop / (strainDiff * strainDiff);

      // dLdx = pLpx + pLpet * detdx + pLpeb * debdx
      dLdh = axStrainBottom / strainDiff + dLdstrainTop * dStrainTopdh +
             dLdstrainBottom * dStrainBottomdh;
      dLdt = dLdstrainTop * dStrainTopdt + dLdstrainBottom * dStrainBottomdt;

      averageStrain = 0.5 * axStrainBottom;
      dAverageStraindh = 0.5 * dStrainBottomdh;
      dAverageStraindt = 0.5 * dStrainBottomdt;
    }

    TacsScalar dFactordL, dFactordt;
    webCrippleFactor = computeCripplingFactorSens(
        compressiveLength, this->stiffenerThick, dFactordL, dFactordt);
    const TacsScalar dFactordh = dFactordL * dLdh;
    dFactordt += dFactordL * dLdt;

    plyStrain[0] = averageStrain;
    for (int ii = 0; ii < numPlies; ii++) {
      TacsScalar plyFailStrainSens[3];
      const TacsScalar strengthRatio =
          ply->failureStrainSens(plyAngles[ii], plyStrain, plyFailStrainSens);
      if (computeThicknessSens) {
        const TacsScalar dFaildt =
            strengthRatio * dFactordt +
            webCrippleFactor * plyFailStrainSens[0] * dAverageStraindt;

        dfdx[tNum] += dKSdFail[ii + numPlies] * scale * dFaildt;
      }
      if (computeHeightSens) {
        const TacsScalar dFaildh =
            strengthRatio * dFactordh +
            webCrippleFactor * plyFailStrainSens[0] * dAverageStraindh;

        dfdx[hNum] += dKSdFail[ii + numPlies] * scale * dFaildh;
      }
    }
  }
}

TacsScalar TACSBladeStiffenedShellConstitutive::computeCripplingFactorSens(
    const TacsScalar b, const TacsScalar t, TacsScalar& dkdb,
    TacsScalar& dkdt) {
  const TacsScalar factor = pow(b / t, 0.717) / 1.63;
  dkdb = 0.717 * factor / b;
  dkdt = -0.717 * factor / t;
  return factor;
}

// ==============================================================================
// Utility functions
// ==============================================================================

void printStiffnessMatrix(const TacsScalar* const C) {
  const TacsScalar* A = &C[0];
  const TacsScalar* B = &C[6];
  const TacsScalar* D = &C[12];
  const TacsScalar* As = &C[18];
  TacsScalar drill = C[21];

  printf("[\n");
  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      TacsRealPart(A[0]), TacsRealPart(A[1]), TacsRealPart(A[2]),
      TacsRealPart(B[0]), TacsRealPart(B[1]), TacsRealPart(B[2]), 0., 0., 0.);
  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      TacsRealPart(A[1]), TacsRealPart(A[3]), TacsRealPart(A[4]),
      TacsRealPart(B[1]), TacsRealPart(B[3]), TacsRealPart(B[4]), 0., 0., 0.);
  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      TacsRealPart(A[2]), TacsRealPart(A[4]), TacsRealPart(A[5]),
      TacsRealPart(B[2]), TacsRealPart(B[4]), TacsRealPart(B[5]), 0., 0., 0.);

  printf(
      "--------------------------------------------------------------------"
      "----"
      "--------------------------------------------------------\n");

  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      TacsRealPart(B[0]), TacsRealPart(B[1]), TacsRealPart(B[2]),
      TacsRealPart(D[0]), TacsRealPart(D[1]), TacsRealPart(D[2]), 0., 0., 0.);
  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      TacsRealPart(B[1]), TacsRealPart(B[3]), TacsRealPart(B[4]),
      TacsRealPart(D[1]), TacsRealPart(D[3]), TacsRealPart(D[4]), 0., 0., 0.);
  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      TacsRealPart(B[2]), TacsRealPart(B[4]), TacsRealPart(B[5]),
      TacsRealPart(D[2]), TacsRealPart(D[4]), TacsRealPart(D[5]), 0., 0., 0.);

  printf(
      "--------------------------------------------------------------------"
      "----"
      "--------------------------------------------------------\n");

  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      0., 0., 0., 0., 0., 0., TacsRealPart(As[0]), TacsRealPart(As[1]), 0.);
  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      0., 0., 0., 0., 0., 0., TacsRealPart(As[1]), TacsRealPart(As[2]), 0.);

  printf(
      "--------------------------------------------------------------------"
      "----"
      "--------------------------------------------------------\n");

  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      0., 0., 0., 0., 0., 0., 0., 0., TacsRealPart(drill));
  printf("]\n");
}
