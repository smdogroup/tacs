/*
=============================================================================
Blade-Stiffened Shell Constitutive Model
=============================================================================
@File    :   TACSBladeStiffenedShellConstitutive.h
@Date    :   2023/04/05
@Author  :   Alasdair Christison Gray
@Description : Constitutive model for a blade-stiffened shell. Based on the
bladeFSDT model from previous versions of TACS developed by Graeme Kennedy.
*/

#pragma once

// =============================================================================
// Standard Library Includes
// =============================================================================

// =============================================================================
// Extension Includes
// =============================================================================
#include "TACSMaterialProperties.h"
#include "TACSShellConstitutive.h"

// =============================================================================
// Global constant definitions
// =============================================================================

// =============================================================================
// Class Declaration
// =============================================================================

/*
    This constitutive class models a shell stiffened with T-shaped stiffeners.
    The stiffeners are not explicitly modelled. Instead, their stiffness is
"smeared" across the shell.

                        |                                    |         ^
                        |                                    |         |
                        |                                    |         |
                        |                                    |       height
                        | <--------------pitch-------------> |         |
Stiffener               |                                    |         |
thickness  ---> -----------------                    ----------------- v
    ---------------------------------------------------------------------- <---
Panel thickness

The panel and stiffener are modelled as laminates, both of which are
parameterised by an arbitrary number of ply angles and fractions. The laminate
properties are then computed using a smeared stiffness approach, which assumes
that the various ply angles are distributed evenly throughout the thicknes of
the panel. This is a valid assumption for thick laminates.

WARNING: If you use ply fractions as design variables, it is currently up to you
to add constraints to ensure that the sum of the fractions is equal to 1.0. This
is not done automatically (yet).
*/
class TACSBladeStiffenedShellConstitutive : public TACSShellConstitutive {
 public:
  /**
   * @brief Construct a new TACSBladeStiffenedShellConstitutive object
   *
   * @param _panelPly Orthotropic ply object for the panel
   * @param _stiffenerPly Orthotropic ply object for the stiffener
   * @param _kcorr Shear correction factor
   * @param _panelLength Panel length value
   * @param _panelLengthNum Panel length design variable number
   * @param _stiffenerPitch Stiffener pitch value
   * @param _stiffenerPitchNum Stiffener pitch design variable number
   * @param _stiffenerHeight Stiffener height value
   * @param _stiffenerHeightNum Stiffener height design variable number
   * @param _stiffenerThick Stiffener thickness value
   * @param _stiffenerThickNum Stiffener thickness design variable number
   * @param _panelThick Panel thickness value
   * @param _panelThickNum Panel thickness design variable number
   * @param _numPanelPlies Number of ply angles in the panel laminate
   * @param _panelPlyAngles Panel ply angles
   * @param _panelPlyFracNums Panel ply fraction design variable numbers
   * @param _numStiffenerPlies Number of ply angles in the stiffener laminate
   * @param _stiffenerPlyAngles Stiffener ply angles
   * @param _stiffenerPlyFracNums Stiffener ply fraction design variable numbers
   * @param _flangeFraction Stiffener base width as a fraction of the stiffener
   * height
   */
  TACSBladeStiffenedShellConstitutive(
      TACSOrthotropicPly* _panelPly, TACSOrthotropicPly* _stiffenerPly,
      TacsScalar _kcorr, TacsScalar _panelLength, int _panelLengthNum,
      TacsScalar _stiffenerPitch, int _stiffenerPitchNum,
      TacsScalar _stiffenerHeight, int _stiffenerHeightNum,
      TacsScalar _stiffenerThick, int _stiffenerThickNum,
      TacsScalar _panelThick, int _panelThickNum, int _numPanelPlies,
      TacsScalar _panelPlyAngles[], int _panelPlyFracNums[],
      int _numStiffenerPlies, TacsScalar _stiffenerPlyAngles[],
      int _stiffenerPlyFracNums[], TacsScalar _flangeFraction = 1.0);

  ~TACSBladeStiffenedShellConstitutive();

  // ==============================================================================
  // Set non-default values
  // ==============================================================================
  /**
   * @brief Set the Stiffener Pitch DV Bounds
   *
   * @param lowerBound Lower bound
   * @param upperBound Upper bound
   */
  void setStiffenerPitchBounds(TacsScalar lowerBound, TacsScalar upperBound);

  /**
   * @brief Set the Stiffener Height DV Bounds
   *
   * @param lowerBound Lower bound
   * @param upperBound Upper bound
   */
  void setStiffenerHeightBounds(TacsScalar lowerBound, TacsScalar upperBound);

  /**
   * @brief Set the Stiffener Thickness DV Bounds
   *
   * @param lowerBound Lower bound
   * @param upperBound Upper bound
   */
  void setStiffenerThicknessBounds(TacsScalar lowerBound,
                                   TacsScalar upperBound);

  /**
   * @brief Set the Panel Pitch DV Bounds
   *
   * @param lowerBound Lower bound
   * @param upperBound Upper bound
   */
  void setPanelThicknessBounds(TacsScalar lowerBound, TacsScalar upperBound);

  /**
   * @brief Set the Stiffener Ply Fraction DV Bounds
   *
   * @param lowerBounds Lower bounds for each ply fraction
   * @param upperBounds Upper bounds for each ply fraction
   */
  void setStiffenerPlyFractionBounds(TacsScalar lowerBounds[],
                                     TacsScalar upperBounds[]);

  /**
   * @brief Set the Panel Ply Fraction DV Bounds
   *
   * @param lowerBounds Lower bounds for each ply fraction
   * @param upperBounds Upper bounds for each ply fraction
   */
  void setPanelPlyFractionBounds(TacsScalar lowerBounds[],
                                 TacsScalar upperBounds[]);

  /**
   * @brief Set the Stiffener Ply Fractions
   *
   * @param plyFractions Ply fractions
   */
  void setStiffenerPlyFractions(TacsScalar plyFractions[]);

  /**
   * @brief Set the Panel Ply Fractions
   *
   * @param plyFractions Ply fractions
   */
  void setPanelPlyFractions(TacsScalar plyFractions[]);

  /**
   * @brief Set the weight used for aggregating the failure criteria through the
   * shell thickness
   *
   * @param _ksWeight
   */
  inline void setKSWeight(double _ksWeight) { ksWeight = _ksWeight; }

  // ==============================================================================
  // Setting/getting design variable information
  // ==============================================================================
  // Retrieve the global design variable numbers
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);

  // Set the element design variable from the design vector
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);

  // Get the element design variables values
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);

  // Get the lower and upper bounds for the design variable values
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]);

  // ==============================================================================
  // Evaluate mass properties
  // ==============================================================================
  // Evaluate the mass per unit area
  TacsScalar evalDensity(int elemIndex, const double pt[],
                         const TacsScalar X[]);

  // Add the derivative of the density w.r.t. the design variables
  void addDensityDVSens(int elemIndex, TacsScalar scale, const double pt[],
                        const TacsScalar X[], int dvLen, TacsScalar dfdx[]);

  // Evaluate the mass moments
  void evalMassMoments(int elemIndex, const double pt[], const TacsScalar X[],
                       TacsScalar moments[]);

  // Add the sensitivity of the mass moments
  void addMassMomentsDVSens(int elemIndex, const double pt[],
                            const TacsScalar X[], const TacsScalar scale[],
                            int dvLen, TacsScalar dfdx[]);

 protected:
  // ==============================================================================
  // Helper functions for computing stiffener cross-section properties
  // ==============================================================================
  /**
   * @brief Compute the stiffener area,
   *
   * @return TacsScalar Stiffener area
   */
  inline TacsScalar computeStiffenerArea();

  /**
   * @brief Compute the derivative of the stiffener's cross-sectional area with
   * respect to the stiffener thickness and height
   *
   * @param dAdt Derivative of the stiffener's cross-sectional area with respect
   * to the stiffener thickness
   * @param dAdh Derivative of the stiffener's cross-sectional area with respect
   * to the stiffener height
   */
  inline void computeStiffenerAreaSens(TacsScalar& dAdt, TacsScalar& dAdh);

  /**
   * @brief Compute the stiffener centroid height, this is the height relative
   * to the bottom surface of the stiffener, NOT the mid-plane of the panel.
   *
   * @return TacsScalar The stiffener centroid height
   */
  inline TacsScalar computeStiffenerCentroidHeight();

  /**
   * @brief Compute the derivative of the stiffener centroid height with respect
   * to the stiffener thickness and height
   *
   * @param dzdt Derivative of the stiffener centroid height with respect to the
   * stiffener thickness
   * @param dzdh Derivative of the stiffener centroid height with respect to the
   * stiffener height
   */
  inline void computeStiffenerCentroidHeightSens(TacsScalar& dzdt,
                                                 TacsScalar& dzdh);

  /**
   * @brief Compute the second moment of area about the stiffener centroid
   *
   * @return TacsScalar Second moment of area about the stiffener centroid
   */
  inline TacsScalar computeStiffenerIzz();

  /**
   * @brief Compute the derivative of the stiffener's second moment of area
   * about the stiffener centroid with respect to the stiffener thickness and
   * height
   *
   * @param dIdt Derivative of Izz with respect to the stiffener thickness
   * @param dIdh Derivative of Izz with respect to the stiffener height
   */
  inline void computeStiffenerIzzSens(TacsScalar& dIdt, TacsScalar& dIdh);

  /**
   * @brief Compute the torsion constant of the stiffener cross section
   *
   * The equation used is taken from page 9 of "TORSIONAL SECTION PROPERTIES OF
   * STEEL SHAPES" by the Canadian institute of steel construction
   * (http://dir.cisc-icca.ca/files/technical/techdocs/updates/torsionprop.pdf).
   * The use of this equation is obviously a little suspect for a composite
   * stiffener, but it is the best we have for now.
   *
   * @return TacsScalar The stiffener torsion constant
   */
  inline TacsScalar computeStiffenerJxx();

  /**
   * @brief Compute the derivative of the stiffener's torsion constant with
   * respect to the stiffener thickness and height
   *
   * @param dJdt Derivative of Jxx with respect to the stiffener thickness
   * @param dJdh Derivative of Jxx with respect to the stiffener height
   */
  inline void computeStiffenerJxxSens(TacsScalar& dJdt, TacsScalar& dJdh);

  /**
   * @brief Compute the moment of inertia of the stiffener cross section about
   * it's centroid, NOT the shell mid-plane
   *
   * @return TacsScalar The stiffener moment of inertia (Iyy) about it's
   * centroid
   */
  inline TacsScalar computeStiffenerMOI();

  /**
   * @brief Compute the derivative of the stiffener's moment of inertia about
   * its centroid with respect to the stiffener thickness and height
   *
   * @param dMOIdt Derivative of Iyy with respect to the stiffener thickness
   * @param dMOIdh Derivative of Iyy with respect to the stiffener height
   */
  inline void computeStiffenerMOISens(TacsScalar& dMOIdt, TacsScalar& dMOIdh);

  // ==============================================================================
  // Attributes
  // ==============================================================================

  // --- General attributes ---
  int numPanelPlies;       ///< Number of plies in the panel laminate
  int numStiffenerPlies;   ///< Number of plies in the stiffener laminate
  double ksWeight = 80.0;  ///< Failure criteria KS aggregation weight
  int numDesignVars;       ///< Number of design variables

  // --- Material properties ---
  TACSOrthotropicPly* panelPly;      ///< Orthotropic ply for the panel
  TACSOrthotropicPly* stiffenerPly;  ///< Orthotropic ply for the stiffener
  TacsScalar kcorr;                  ///< Shear correction factor

  // --- Design variable values ---
  TacsScalar panelLength;         ///< Panel length
  TacsScalar stiffenerPitch;      ///< Stiffener pitch
  TacsScalar stiffenerHeight;     ///< Stiffener height
  TacsScalar stiffenerThick;      ///< Stiffener thickness
  TacsScalar panelThick;          ///< Panel thickness
  TacsScalar* panelPlyFracs;      ///< Panel ply Fracs
  TacsScalar* stiffenerPlyFracs;  ///< Stiffener ply Fracs

  // --- Fixed parameters ---
  TacsScalar* panelPlyAngles;      ///< Panel ply angles
  TacsScalar* stiffenerPlyAngles;  ///< Stiffener ply angles
  TacsScalar flangeFraction;  ///< Flange fraction (stiffener base width as a
                              ///< fraction of stiffener height)

  // --- Design variable bounds ---
  TacsScalar panelLengthLowerBound = 0.0;
  TacsScalar panelLengthUpperBound = 1e20;
  TacsScalar stiffenerPitchLowerBound = 1e-3;
  TacsScalar stiffenerPitchUpperBound = 1e20;
  TacsScalar stiffenerHeightLowerBound = 1e-3;
  TacsScalar stiffenerHeightUpperBound = 1e20;
  TacsScalar stiffenerThickLowerBound = 1e-4;
  TacsScalar stiffenerThickUpperBound = 1e20;
  TacsScalar panelThickLowerBound = 1e-4;
  TacsScalar panelThickUpperBound = 1e20;
  TacsScalar* panelPlyFracLowerBounds;
  TacsScalar* panelPlyFracUpperBounds;
  TacsScalar* stiffenerPlyFracLowerBounds;
  TacsScalar* stiffenerPlyFracUpperBounds;

  // --- Design variable numbers ---
  int panelLengthNum;         ///< Panel length DV number
  int stiffenerPitchNum;      ///< Stiffener pitch DV number
  int stiffenerHeightNum;     ///< Stiffener height DV number
  int stiffenerThickNum;      ///< Stiffener thickness DV number
  int panelThickNum;          ///< Panel thickness DV number
  int* panelPlyFracNums;      ///< Panel ply fraction DV numbers
  int* stiffenerPlyFracNums;  ///< Stiffener ply fraction DV numbers

  // --- Local design variable numbers, these are useful when returning
  // sensitivity values ---
  int panelLengthLocalNum;         ///< Panel length local DV number
  int stiffenerPitchLocalNum;      ///< Stiffener pitch local DV number
  int stiffenerHeightLocalNum;     ///< Stiffener height local DV number
  int stiffenerThickLocalNum;      ///< Stiffener thickness local DV number
  int panelThickLocalNum;          ///< Panel thickness local DV number
  int* panelPlyFracLocalNums;      ///< Panel ply fraction local DV numbers
  int* stiffenerPlyFracLocalNums;  ///< Stiffener ply fraction local DV numbers

  static const int PANEL_LENGTH_INDEX = 0;
  static const int STIFFENER_PITCH_INDEX = 1;
  static const int STIFFENER_HEIGHT_INDEX = 2;
  static const int STIFFENER_THICKNESS_INDEX = 3;
  static const int PANEL_THICKNESS_INDEX = 4;
};
