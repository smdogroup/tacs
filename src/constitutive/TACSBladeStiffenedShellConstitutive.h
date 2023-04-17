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
#include "TACSBeamConstitutive.h"
#include "TACSMaterialProperties.h"
#include "TACSShellConstitutive.h"
#include "TacsUtilities.h"

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
   * @param _panelThick Panel thickness value
   * @param _panelThickNum Panel thickness design variable number
   * @param _numPanelPlies Number of ply angles in the panel laminate
   * @param _panelPlyAngles Panel ply angles
   * @param _panelPlyFracs Panel ply fractions
   * @param _panelPlyFracNums Panel ply fraction design variable numbers
   * @param _stiffenerHeight Stiffener height value
   * @param _stiffenerHeightNum Stiffener height design variable number
   * @param _stiffenerThick Stiffener thickness value
   * @param _stiffenerThickNum Stiffener thickness design variable number
   * @param _numStiffenerPlies Number of ply angles in the stiffener laminate
   * @param _stiffenerPlyAngles Stiffener ply angles
   * @param _stiffenerPlyFracs Stiffener ply fractions
   * @param _stiffenerPlyFracNums Stiffener ply fraction design variable numbers
   * @param _flangeFraction Stiffener base width as a fraction of the stiffener
   * height
   */
  TACSBladeStiffenedShellConstitutive(
      TACSOrthotropicPly* _panelPly, TACSOrthotropicPly* _stiffenerPly,
      TacsScalar _kcorr, TacsScalar _panelLength, int _panelLengthNum,
      TacsScalar _stiffenerPitch, int _stiffenerPitchNum,
      TacsScalar _panelThick, int _panelThickNum, int _numPanelPlies,
      TacsScalar _panelPlyAngles[], TacsScalar _panelPlyFracs[],
      int _panelPlyFracNums[], TacsScalar _stiffenerHeight,
      int _stiffenerHeightNum, TacsScalar _stiffenerThick,
      int _stiffenerThickNum, int _numStiffenerPlies,
      TacsScalar _stiffenerPlyAngles[], TacsScalar _stiffenerPlyFracs[],
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

  // ==============================================================================
  // Evaluate thermal properties
  // ==============================================================================
  // Evaluate the specific heat, Not implemented for this class
  TacsScalar evalSpecificHeat(int elemIndex, const double pt[],
                              const TacsScalar X[]);

  // ==============================================================================
  // Compute stress/strain/stiffness
  // ==============================================================================
  // Evaluate the stress
  void evalStress(int elemIndex, const double pt[], const TacsScalar X[],
                  const TacsScalar e[], TacsScalar s[]);

  // Add the derivative of the product of stress with a vector psi to dfdx
  void addStressDVSens(int elemIndex, TacsScalar scale, const double pt[],
                       const TacsScalar X[], const TacsScalar strain[],
                       const TacsScalar psi[], int dvLen, TacsScalar dfdx[]);

  // Evaluate the tangent stiffness
  void evalTangentStiffness(int elemIndex, const double pt[],
                            const TacsScalar X[], TacsScalar C[]);

  // ==============================================================================
  // Compute failure criteria
  // ==============================================================================
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

  // ==============================================================================
  // Compute output quantities
  // ==============================================================================
  // Get the object name
  const char* getObjectName() { return constName; }

  // Retrieve the design variable for plotting purposes
  TacsScalar evalDesignFieldValue(int elemIndex, const double pt[],
                                  const TacsScalar X[], int index);

 protected:
  /**
   * @brief Compute the stiffness matrix of the stiffened shell
   *
   * @param C Array to store the stiffness matrix in (will be zeroed out within
   * this function)
   */
  void computeStiffness(TacsScalar C[]);

  /**
   * @brief Compute the Q and ABar Matrices for a laminate based on the Q and
   * ABar matrices for each ply
   *
   * @param numPlies Number of plies in the laminate
   * @param plyFractions Ply fractions for each ply angle
   * @param Q The Q matrix for the laminate, stored as a flattened 6 entry array
   * ([Q11, Q12, Q16, Q22, Q26, Q66])
   * @param ABar The ABar matrix for the laminate, stored as a flattened 3 entry
   * array ([Q44, Q45, Q55])
   */
  inline void computeSmearedStiffness(const int numPlies,
                                      const TacsScalar* const QMats,
                                      const TacsScalar* const AbarMats,
                                      const TacsScalar plyFractions[],
                                      TacsScalar Q[], TacsScalar ABar[]);

  /**
   * @brief Compute the failure values for each failure mode of the stiffened
   * panel
   *
   * @param e Shell strains
   * @param fail Array to store the failure values in
   * @return TacsScalar The aggregated failure value
   */
  TacsScalar computeFailureValues(const TacsScalar e[], TacsScalar fail[]);

  // ==============================================================================
  // Helper functions for transforming strains/stresses/stiffnesses between the
  // panel and stiffener
  // ==============================================================================
  /**
   * @brief Given the shell mid-plane strains, compute the equivalent beam
   * strains at the stiffener centroid
   *
   * The transformation can be written as [eb] = [Te] [es], or in full:
   * [e11b]   [1, 0, 0,   zc, 0, 0,    0, 0]  [e11s]
   * [xi1b]   [0, 0, 0,   0,  0, -1/2, 0, 0]  [e22s]
   * [xi2b] = [0, 0, 0,   1,  0, 0,    0, 0]  [y12s]
   * [xi3b]   [0, 0, 0,   0,  0, 0,    0, 0]  [k11s]
   * [y13b]   [0, 0, 0,   0,  0, 0,    0, 1]  [k22s]
   * [y12b]   [0, 0, 1/2, 0,  0, zc/2, 0, 0]  [k12s]
   *                                          [y23s]
   *                                          [y13s]
   *
   * Where:
   * - eb is the beam strain vector
   * - es is the shell strain vector
   * - zc is the distance from the shell mid-plane to the stiffener centroid
   *
   * @param panelStrain The shell mid-plane strains [e11, e22, y12, k11, k22,
   * k12, y23, y13]
   * @param stiffenerStrain The stiffener centroid strains [e11, xi1, xi2, xi3,
   * y13, y12]
   */
  inline void transformStrain(const TacsScalar panelStrain[],
                              TacsScalar stiffenerStrain[]);

  /**
   * @brief Transform a sensitivity w.r.t the stiffener centroid strains to a
   * sensitivity w.r.t the shell mid-plane strains
   *
   * This transformation can be written as [df/des] = [df/deb] * [deb/des] =
   * [df/deb] * [Te], or in full:
   * [df/de11s]   [df/de11b]   [1, 0, 0,   zc, 0, 0,    0, 0]
   * [df/de22s]   [df/dxi1b]   [0, 0, 0,   0,  0, -1/2, 0, 0]
   * [df/dy12s] = [df/dxi2b] * [0, 0, 0,   1,  0, 0,    0, 0]
   * [df/dk11s]   [df/dxi3b]   [0, 0, 0,   0,  0, 0,    0, 0]
   * [df/dk22s]   [df/dy13b]   [0, 0, 0,   0,  0, 0,    0, 1]
   * [df/dk12s]   [df/dy12b]   [0, 0, 1/2, 0,  0, zc/2, 0, 0]
   * [df/dy23s]
   * [df/dy13s]
   *
   * Where:
   * - df/des is the sensitivity w.r.t the shell mid-plane strains
   * - df/deb is the sensitivity w.r.t the stiffener centroid strains
   * - deb/des = Te is the strain transformation matrix
   *
   * @param stiffenerStrainSens Array containing sensitivity of an output w.r.t
   * beam strains
   * @param panelStrainSens Array to store the sensitivity of the output w.r.t
   * shell mid-plane strains
   */
  inline void transformStrainSens(const TacsScalar stiffenerStrainSens[],
                                  TacsScalar panelStrainSens[]);

  /**
   * @brief Add the DV sensitivity of the product of the strain
   * transformation matrix with two vectors
   *
   * This can be written as dfdx += scale * [lhs] * d/dx([Te]) * [rhs]
   *
   * A common use for this function is to convert the derivative of a quantity
   * w.r.t the stiffener centroid strain to a derivative w.r.t the design
   * variables, since the stiffener centroid strains depend on the offset of the
   * stiffener centroid from the shell mid-plane, which depends on the design
   * variables.
   *
   * This operation can be written as [df/dx] += scale * [df/deBeam] *
   * d/dx([Te]) * [eShell]
   *
   * @param lhs The left-hand-side vector, should be as long as a beam strain
   * array
   * @param rhs The right hand side vector, should be as long as a shell strain
   * array
   * @param scale The scaling factor
   * @param dfdx The sensitivities of the output w.r.t the design variables to
   * be added to
   */
  inline void addStrainTransformProductDVsens(const TacsScalar lhs[],
                                              const TacsScalar rhs[],
                                              const TacsScalar scale,
                                              TacsScalar dfdx[]);

  /**
   * @brief Add the contribution of the stiffener stress to the panel stress
   *
   * The transformation of stiffener to panel stress consists of two steps:
   * 1. Transform the stiffener forces to the panel mid-plane, this is the
   * transpose of the strain transformation
   * 2. Add the transformed stiffener forces to the panel mid-plane forces,
   * divided by the stiffener pitch to account for the stiffener spacing
   *
   * This can be written as [ss] += 1/pitch * [Te]^{T} [sb]
   * where:
   * - [Te] is the strain transformation matrix
   * - [ss] is the shell mid-plane stress vector
   * - [sb] is the stiffener beam stress vector
   *
   * @param stiffenerStress Stiffener beam stresses [F1, M1, M2, M3, V13, V12]
   * @param panelStress The shell mid-plane stresses to be added to [N11, N22,
   * N12, M11, M22, M12, Q23, Q13]
   */
  inline void addStiffenerStress(const TacsScalar stiffenerStress[],
                                 TacsScalar panelStress[]);

  /**
   * @brief Add the contribution of the stiffener stiffness to the panel
   * stiffness
   *
   * This operation can be written as: [Cs] += 1 / pitch * [Te]^{T} * [Cb] *
   * [Te] Where:
   * - [Te] is the strain transformation matrix
   * - [Cs] is the shell stiffness matrix
   * - [Cb] is the stiffener beam stiffness matrix
   *
   * @param stiffenerStiffness The stiffener stiffness matrix, computed by
   * `computeStiffenerStiffness`
   * @param panelStiffness The panel stiffness matrix to be added to
   */
  inline void addStiffenerStiffness(const TacsScalar stiffenerStiffness[],
                                    TacsScalar panelStiffness[]);

  // ==============================================================================
  // Helper functions for computing the panel stress/stiffness/failure
  // ==============================================================================
  // In future, these methods should be replaced by calls to another shell
  // constitutive model

  /**
   * @brief Compute the stiffness matrix of the panel (without stiffeners)
   *
   * @param C Array to store the stiffness matrix in (will be zeroed out within
   * this function)
   */
  void computePanelStiffness(TacsScalar C[]);

  /**
   * @brief Compute the panel stress given the panel strains
   *
   * @param strain Shell strains [e11, e22, y12, k11, k22, k12, y23, y13]
   * @param stress Shell stresses [N11, N22, N12, M11, M22, M12, Q23, Q13]
   */
  void computePanelStress(const TacsScalar strain[], TacsScalar stress[]);

  /**
   * @brief Add the derivative of the product of panel stresses with a vector
   * psi to dfdx
   *
   * @param scale Scale factor to apply to the derivatives
   * @param strain SHell mid-plane strains
   * @param psi Array multiplying stresses
   * @param dfdx Design variable sensitivity array to add to
   */
  void addPanelStressDVSens(const TacsScalar scale, const TacsScalar strain[],
                            const TacsScalar psi[], TacsScalar dfdx[]);

  /**
   * @brief Compute the failure criteria in the panel
   *
   * @param strain Shell strains [e11, e22, y12, k11, k22, k12, y23, y13]
   * @return TacsScalar The aggregated failure value for the panel
   */
  TacsScalar computePanelFailure(const TacsScalar strain[]);

  /**
   * @brief Compute the sensitivity of the panel failure w.r.t the panel strains
   *
   * @param strain Shell strains [e11, e22, y12, k11, k22, k12, y23, y13]
   * @param sens Array to store the sensitivity of the failure w.r.t the strains
   * @return TacsScalar The aggregated failure value for the panel
   */
  TacsScalar evalPanelFailureStrainSens(const TacsScalar strain[],
                                        TacsScalar sens[]);

  /**
   * @brief Add the derivative of the panel's failure w.r.t it's DVs
   *
   * @param strain Shell strains [e11, e22, y12, k11, k22, k12, y23, y13]
   * @param scale The scale factor to apply to the derivatives
   * @param dfdx The array to add the derivative to
   */
  void addPanelFailureDVSens(const TacsScalar strain[], const TacsScalar scale,
                             TacsScalar dfdx[]);

  // ==============================================================================
  // Helper functions for computing the stiffner's strain/stress/stiffness
  // ==============================================================================
  // In future, these methods should be replaced by calls to a beam constitutive
  // model

  /**
   * @brief Compute the stiffener stress given the stiffener strains
   *
   * @param stiffenerStrain The stiffener centroid strains  [e11, xi1, xi2, xi3,
   * y13, y12]
   * @param stiffenerStress The stiffener stresses/forces [N11, M1, M2, M3, V12,
   * V13]
   */
  inline void computeStiffenerStress(const TacsScalar stiffenerStrain[],
                                     TacsScalar stiffenerStress[]);

  /**
   * @brief Add the derivative of the product of stiffener stresses with a
   * vector w.r.t the design variables
   *
   * @param scale Scale factor to apply to the derivatives
   * @param strain Stiffener centroid beam strains
   * @param psi Array multiplying stresses
   * @param dfdx Design variable sensitivity array to add to
   */
  void addStiffenerStressDVSens(const TacsScalar scale,
                                const TacsScalar strain[],
                                const TacsScalar psi[], TacsScalar dfdx[]);

  /**
   * @brief Compute the stiffener's stiffness matrix
   *
   * @param C The stiffener local stiffness matrix (The derivative of the beam
   * stresses w.r.t the beam strains)
   */
  inline void computeStiffenerStiffness(TacsScalar C[]);

  inline void computeStiffenerModuli(TacsScalar& E, TacsScalar& G);

  /**
   * @brief Compute the failure criterion for the stiffener in a given strain
   * state
   *
   * @param stiffenerStrain The stiffener centroid strains  [e11, xi1, xi2, xi3,
   * y13, y12]
   * @return TacsScalar The failure criterion for the stiffener
   */
  TacsScalar computeStiffenerFailure(const TacsScalar stiffenerStrain[]);

  /**
   * @brief Compute the sensitivity of the stiffener failure w.r.t the stiffener
   * strains
   *
   * @param strain Beam strains [e11, xi1, xi2, xi3, y13, y12]
   * @param sens Array to store the sensitivity of the failure w.r.t the strains
   * @return TacsScalar The aggregated failure value for the stiffener
   */
  TacsScalar evalStiffenerFailureStrainSens(const TacsScalar strain[],
                                            TacsScalar sens[]);

  /**
   * @brief Add the derivative of the stiffener's failure w.r.t it's DVs
   *
   * @param strain Beam strains [e11, e22, y12, k11, k22, k12, y23, y13]
   * @param scale The scale factor to apply to the derivatives
   * @param dvLen The number of stiffener design variables
   * @param dfdx The array to add the derivative to
   */
  void addStiffenerFailureDVSens(const TacsScalar strain[],
                                 const TacsScalar scale, TacsScalar dfdx[]);

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
  int numGeneralDV;        ///< Number of general DVs
  int numPanelDV;          ///< Number of panel DVs
  int numStiffenerDV;      ///< Number of stiffener DVs

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
  TacsScalar* panelQMats;     ///< Panel Q-matrices (they are constant since the
                              ///< ply angles are fixed)
  TacsScalar* stiffenerQMats;  ///< Panel Q-matrices (they are constant since
                               ///< the ply angles are fixed)
  TacsScalar* panelAbarMats;
  TacsScalar* stiffenerAbarMats;

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
  int panelThickLocalNum;          ///< Panel thickness local DV number
  int* panelPlyFracLocalNums;      ///< Panel ply fraction local DV numbers
  int stiffenerPitchLocalNum;      ///< Stiffener pitch local DV number
  int stiffenerHeightLocalNum;     ///< Stiffener height local DV number
  int stiffenerThickLocalNum;      ///< Stiffener thickness local DV number
  int* stiffenerPlyFracLocalNums;  ///< Stiffener ply fraction local DV numbers
  int panelDVStartNum;             ///< Panel DV start number
  int stiffenerDVStartNum;         ///< Stiffener DV start number

  // --- Arrays for storing failure values and their
  // sensitivities ---
  TacsScalar* panelPlyFailValues;
  TacsScalar* stiffenerPlyFailValues;
  TacsScalar** panelPlyFailStrainSens;
  TacsScalar** stiffenerPlyFailStrainSens;
  TacsScalar* panelPlyFailSens;
  TacsScalar* stiffenerPlyFailSens;

  static const char* constName;        ///< Constitutive model name
  static const int NUM_Q_ENTRIES = 6;  ///< Number of entries in the Q matrix
  static const int NUM_ABAR_ENTRIES =
      3;                              ///< Number of entries in the ABar matrix
  static const int NUM_FAILURES = 2;  ///< Number of failure modes
};
