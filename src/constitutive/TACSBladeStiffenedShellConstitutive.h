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

void printStiffnessMatrix(const TacsScalar* const C);

// =============================================================================
// Class Declaration
// =============================================================================

/**
 * @brief Constitutive model for a blade-stiffened shell
 *
 * This constitutive class models a shell stiffened with T-shaped stiffeners.
 * The stiffeners are not explicitly modelled. Instead, their stiffness is
 * "smeared" across the shell.
 *
 *                         |                      |         ^
 *                         |                      |         |
 *                         |                      |         |
 *                         |                      |       height
 *                         | <-------pitch------> |         |
 * Stiffener               |                      |         v
 * thickness ---> -----------------      -----------------
 *     --------------------------------------------------- <--- Panel thickness
 *
 *     ----------> 2-direction
 *     |
 *     |
 *     |
 *     |
 *     v
 *     3-direction (shell normal)
 *
 * The stiffener is modelled as a T-shaped beam as shown above, and is modelled
 * as pointing in the opposite of the shell normal. The thickness of the
 * stiffener web and flange are the same. The ratio of the stiffener flange
 * width to the stiffener height is controlled by the flangeFraction parameter,
 * which defaults to 1.0 (i.e. the flange width is equal to the stiffener
 * height). This value can be set when creating the constitutive object and is
 * then fixed.
 *
 * The panel and stiffener are modelled as laminates, both of which are
 * parameterised by an arbitrary number of ply angles and fractions. The
 * laminate properties are then computed using a smeared stiffness approach,
 * which assumes that the various ply angles are distributed evenly throughout
 * the thickness of the panel. This is a valid assumption for thick laminates.
 *
 * The design variables (in order) for this constitutive model are:
 * - Panel length
 * - Stiffener pitch
 * - Panel thickness
 * - Panel ply fractions
 * - Stiffener height
 * - Stiffener thickness
 * - Stiffener ply fractions
 *
 * The failure criterion returned by this model combines numerous possible
 * failure modes into a single value. A material failure criterion (which one
 * depends on the Orthotropic ply objects you pass to this class) is computed at
 * the upper and lower surface of the panel and at the tip of the stiffener,
 * this calculation is performed for every ply angle present in the panel and
 * stiffener laminate. Additionally buckling criteria are computed for combined
 * shear and axial buckling for both a global buckling mode (i.e the entire
 * panel buckles) and a local buckling mode (i.e. the panel buckles between a
 * pair of stiffeners). These buckling failure values are aggregated along with
 * the material failure values into a single failure value using KS aggregation.
 * The smoothness and conservatism of this aggregation can be controlled using
 * the `setKSWeight` method.
 *
 * The panel length design variables do not directly affect the stiffness or
 * stresses computed by the model, their only role is to allow the computation
 * of the critical buckling loads used in the buckling failure criteria. When
 * using this constitutive model, you must also add a set of
 * PanelLengthConstraints to force the panel length design variables to be equal
 * to the true physical panel length. Despite seeming needlessly complex, this
 * approach is necessary because this constitutive model has no information
 * about the geometry of your finite element model and so cannot compute the
 * panel length itself.
 *
 * WARNING: If you use ply fractions as design variables, it is currently up to
 * you to add DVConstraints to ensure that the sum of the fractions is equal
 * to 1.0. This is not done automatically (yet).
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
  inline void setKSWeight(double _ksWeight) { this->ksWeight = _ksWeight; }

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

  /**
   * @brief Compute the effective tensile thickness of the stiffened shell
   *
   * This is the thickness of an unstiffened shell with the same tensile
   * stiffness and density as the stiffened shell
   *
   * @return TacsScalar The effective thickness
   */
  inline TacsScalar computeEffectiveThickness() {
    return this->panelThick +
           this->computeStiffenerArea() / this->stiffenerPitch;
  }

  /**
   * @brief Compute the effective bending thickness of the stiffened shell
   *
   * This is the thickness of an unstiffened shell with the same bending
   * stiffness as the stiffened shell. Note this calculation does not currently
   * account for any difference in the modulii of the stiffener and panel.
   *
   * @return TacsScalar The effective thickness
   */
  TacsScalar computeEffectiveBendingThickness();

  // ==============================================================================
  // Tests
  // ==============================================================================

  /**
   * @brief Test `computeCriticalGlobalBucklingStiffnessSens` against finite
   * difference/complex-step
   *
   */
  void testGlobalBucklingStiffnessSens();

  /**
   * @brief Test `computeCriticalShearLoadSens` against finite difference
   */
  static bool testCriticalShearLoadSens(const TacsScalar D1,
                                        const TacsScalar D2,
                                        const TacsScalar D3,
                                        const TacsScalar L);

  /**
   * @brief Test `bucklingEnvelopeSens` against finite difference
   */
  static bool testBucklingEnvelopeSens(const TacsScalar N1,
                                       const TacsScalar N1Crit,
                                       const TacsScalar N12,
                                       const TacsScalar N12Crit);

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
   * @param QMats The Q matrices for each ply, stored as a flattened 6*numPlies
   * entry array
   * @param AbarMats The Abar matrices for each ply, stored as a flattened
   * 3*numPlies entry array
   * @param plyFractions Ply fractions for each ply angle
   * @param Q The Q matrix for the laminate, stored as a flattened 6 entry array
   * ([Q11, Q12, Q16, Q22, Q26, Q66])
   * @param ABar The ABar matrix for the laminate, stored as a flattened 3 entry
   * array ([Q44, Q45, Q55])
   */
  void computeSmearedStiffness(const int numPlies,
                               const TacsScalar* const QMats,
                               const TacsScalar* const AbarMats,
                               const TacsScalar plyFractions[], TacsScalar Q[],
                               TacsScalar ABar[]);

  /**
   * @brief Compute the Q matrix for a laminate based on the Q matrices for each
   * ply
   *
   * @param numPlies Number of plies in the laminate
   * @param QMats The Q matrices for each ply, stored as a flattened 6*numPlies
   * entry array
   * @param plyFractions Ply fractions for each ply angle
   * @param Q The Q matrix for the laminate, stored as a flattened 6 entry array
   * ([Q11, Q12, Q16, Q22, Q26, Q66])
   */
  void computeSmearedStiffness(const int numPlies,
                               const TacsScalar* const QMats,
                               const TacsScalar plyFractions[], TacsScalar Q[]);

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
  void transformStrain(const TacsScalar panelStrain[],
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
  void transformStrainSens(const TacsScalar stiffenerStrainSens[],
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
  void addStrainTransformProductDVsens(const TacsScalar lhs[],
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
  void addStiffenerStress(const TacsScalar stiffenerStress[],
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
  void addStiffenerStiffness(const TacsScalar stiffenerStiffness[],
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
  void computeStiffenerStress(const TacsScalar stiffenerStrain[],
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
  void computeStiffenerStiffness(TacsScalar C[]);

  static void computeEffectiveModulii(const int numPlies,
                                      const TacsScalar QMats[],
                                      const TacsScalar plyFracs[],
                                      TacsScalar* E, TacsScalar* G);

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
  TacsScalar computeStiffenerArea();

  /**
   * @brief Compute the derivative of the stiffener's cross-sectional area with
   * respect to the stiffener thickness and height
   *
   * @param dAdt Derivative of the stiffener's cross-sectional area with respect
   * to the stiffener thickness
   * @param dAdh Derivative of the stiffener's cross-sectional area with respect
   * to the stiffener height
   */
  void computeStiffenerAreaSens(TacsScalar& dAdt, TacsScalar& dAdh);

  /**
   * @brief Compute the stiffener centroid height, this is the height relative
   * to the bottom surface of the stiffener, NOT the mid-plane of the panel.
   *
   * @return TacsScalar The stiffener centroid height
   */
  TacsScalar computeStiffenerCentroidHeight();

  /**
   * @brief Compute the derivative of the stiffener centroid height with respect
   * to the stiffener thickness and height
   *
   * @param dzdt Derivative of the stiffener centroid height with respect to the
   * stiffener thickness
   * @param dzdh Derivative of the stiffener centroid height with respect to the
   * stiffener height
   */
  void computeStiffenerCentroidHeightSens(TacsScalar& dzdt, TacsScalar& dzdh);

  /**
   * @brief Compute the second moment of area about the stiffener centroid
   *
   * @return TacsScalar Second moment of area about the stiffener centroid
   */
  TacsScalar computeStiffenerIzz();

  /**
   * @brief Compute the derivative of the stiffener's second moment of area
   * about the stiffener centroid with respect to the stiffener thickness and
   * height
   *
   * @param dIdt Derivative of Izz with respect to the stiffener thickness
   * @param dIdh Derivative of Izz with respect to the stiffener height
   */
  void computeStiffenerIzzSens(TacsScalar& dIdt, TacsScalar& dIdh);
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
  TacsScalar computeStiffenerJxx();

  /**
   * @brief Compute the derivative of the stiffener's torsion constant with
   * respect to the stiffener thickness and height
   *
   * @param dJdt Derivative of Jxx with respect to the stiffener thickness
   * @param dJdh Derivative of Jxx with respect to the stiffener height
   */
  void computeStiffenerJxxSens(TacsScalar& dJdt, TacsScalar& dJdh);

  /**
   * @brief Compute the moment of inertia of the stiffener cross section about
   * it's centroid, NOT the shell mid-plane
   *
   * @return TacsScalar The stiffener moment of inertia (Iyy) about it's
   * centroid
   */
  TacsScalar computeStiffenerMOI();

  /**
   * @brief Compute the derivative of the stiffener's moment of inertia about
   * its centroid with respect to the stiffener thickness and height
   *
   * @param dMOIdt Derivative of Iyy with respect to the stiffener thickness
   * @param dMOIdh Derivative of Iyy with respect to the stiffener height
   */
  void computeStiffenerMOISens(TacsScalar& dMOIdt, TacsScalar& dMOIdh);

  // ==============================================================================
  // Buckling functions
  // ==============================================================================

  /**
   * @brief Compute the panel + stiffener stiffness values used to compute the
   * global buckling loads
   *
   * @output D1 1-direction bending stiffness
   * @output D2 2-direction bending stiffness
   * @output D3 Twisting stiffness
   */
  void computeCriticalGlobalBucklingStiffness(TacsScalar* D1, TacsScalar* D2,
                                              TacsScalar* D3);

  /**
   * @brief Compute the sensitivities of the panel + stiffener stiffness values
   * used to compute the global buckling loads
   *
   * @param dfdD1 Sensitivity of the output to the 1-direction bending stiffness
   * @param dfdD2 Sensitivity of the output to the 2-direction bending stiffness
   * @param dfdD3 Sensitivity of the output to the twisting stiffness
   * @output spSens Sensitivity of the output w.r.t the stiffener pitch
   * @output tpSens Sensitivity of the output w.r.t the panel thickness
   * @output hsSens Sensitivity of the output w.r.t the stiffener height
   * @output tsSens Sensitivity of the output w.r.t the stiffener thickness
   * @output QstiffSens Sensitivity of the output w.r.t the stiffener stiffener
   * Q matrix entries
   * @output QpanelSens Sensitivity of the output w.r.t the stiffener panel Q
   * matrix entries
   */
  void computeCriticalGlobalBucklingStiffnessSens(
      const TacsScalar dfdD1, const TacsScalar dfdD2, const TacsScalar dfdD3,
      TacsScalar* spSens, TacsScalar* tpSens, TacsScalar* hsSens,
      TacsScalar* tsSens, TacsScalar QstiffSens[], TacsScalar QpanelSens[]);

  /**
   * @brief Compute the critical axial load for the global buckling of the
   * stiffened panel
   *
   * @param D1 1-direction bending stiffness, computed by
   * `computeCriticalGlobalBucklingStiffness`
   * @param L Panel length
   * @return TacsScalar The critical load
   */
  static inline TacsScalar computeCriticalGlobalAxialLoad(const TacsScalar D1,
                                                          const TacsScalar L) {
    return M_PI * M_PI * D1 / (L * L);
  }

  /**
   * @brief Compute the critical axial load for local buckling of the panel
   *
   * @param D11 D11 stiffness
   * @param D22 D22 stiffness
   * @param D12 D12 stiffness
   * @param D66 D66 stiffness
   * @param L Length (in this case, the stiffener pitch)
   * @return TacsScalar The critical load
   */
  static inline TacsScalar computeCriticalLocalAxialLoad(const TacsScalar D11,
                                                         const TacsScalar D22,
                                                         const TacsScalar D12,
                                                         const TacsScalar D66,
                                                         const TacsScalar L) {
    return 2.0 * M_PI * M_PI / (L * L) * (sqrt(D11 * D22) + D12 + 2.0 * D66);
  }

  /**
   * @brief Compute the sensitivity of the critical axial load for local
   * buckling of the panel
   *
   * @param D11 D11 stiffness
   * @param D22 D22 stiffness
   * @param D12 D12 stiffness
   * @param D66 D66 stiffness
   * @param L Length (in this case, the stiffener pitch)
   * @param D11Sens Sensitivity w.r.t the D11 stiffness
   * @param D22Sens Sensitivity w.r.t the D22 stiffness
   * @param D12Sens Sensitivity w.r.t the D12 stiffness
   * @param D66Sens Sensitivity w.r.t the D66 stiffness
   * @param LSens Sensitivity w.r.t the Length
   * @return TacsScalar The critical load
   */
  static TacsScalar computeCriticalLocalAxialLoadSens(
      const TacsScalar D11, const TacsScalar D22, const TacsScalar D12,
      const TacsScalar D66, const TacsScalar L, TacsScalar* D11Sens,
      TacsScalar* D22Sens, TacsScalar* D12Sens, TacsScalar* D66Sens,
      TacsScalar* LSens);

  /**
   * @brief Compute the critical shear load for either local or global buckling
   *
   * The shear buckling loads are calculated based on an infinite-plate
   * solution simply supported along the panel sides. For the local skin
   * buckling calculations we use the infinite plate solution along edges
   * supported by the blade stiffeners, therefore the panel width is
   * equal to the stiffener pitch (sp). For the global-level
   * calculations, we use the panel length equal to the rib pitch, thus
   * the panel width is equal to the local panel length (Lx).
   *
   * In Stroud and Arganoff, the following formula are suggested for the
   * calculation of the buckling loads:
   *
   * xi = sqrt(D1*D2)/D3
   *
   * if xi > 1.0:
   *   Nxy,crit = (4.0/Ly^2)*(D2*D1^3)^(1/4)*(8.125 + 5.05/xi)
   * else:
   *   Nxy,crit = (4.0/Ly^2)*sqrt(D1*D3)*(11.7 + 0.532*xi + 0.938*xi^2)
   *
   * Note that if xi = 1, then D3 = sqrt(D1*D2) and so
   * (4.0/Ly^2)*sqrt(D1*D3) = (4.0/Ly^2)*(D1*D1^3)^(1/4)
   *
   * However, when xi = 1, the two terms inside the brackets are not
   * equal. As a result, there is a discontinuity between the two
   * expressions. To avoid this, we adjust the first formula in a
   * conservative manner. As shown in Lekhnitskii, the limit for xi ->
   * infty is 8.125, while for xi = 1 is 13.17 and for xi = 0, 11.71.
   *
   * The formula in Stroud and Arganoff do not meet these end conditions.
   * We adjust the formula as follows:
   *
   * if xi > 1.0:
   *   Nxy,crit = (4.0/Ly^2)*(D2*D1^3)^(1/4)*(8.125 + 5.045/xi)
   * else:
   *   Nxy,crit = (4.0/Ly^2)*sqrt(D1*D3)*(11.7 + 0.532*xi + 0.938*xi^2)
   *
   * input:
   * @param D1 the longitudinal bending stiffness
   * @param D2 the transverse bending stiffness
   * @param D3 the shear bending stiffness
   * @param L  the side-length of the panel
   *
   * returns:
   * @return N12Crit  the approximate critical buckling load
   */
  static TacsScalar computeCriticalShearLoad(TacsScalar D1, TacsScalar D2,
                                             TacsScalar D3, TacsScalar L);

  /**
   * @brief Compute the sensitivity of the critical shear buckling load
   *
   * @param D1 the longitudinal bending stiffness
   * @param D2 the transverse bending stiffness
   * @param D3 the shear bending stiffness
   * @param L the side-length of the panel
   * @param sD1 Sensitivity of the critical load w.r.t the longitudinal bending
   * stiffness
   * @param sD2 Sensitivity of the critical load w.r.t the transverse bending
   * stiffness
   * @param sD3 Sensitivity of the critical load w.r.t the shear bending
   * stiffness
   * @param sL Sensitivity of the critical load w.r.t the side-length of the
   * panel
   * @return TacsScalar The critical shear buckling load
   */
  static TacsScalar computeCriticalShearLoadSens(
      const TacsScalar D1, const TacsScalar D2, const TacsScalar D3,
      const TacsScalar L, TacsScalar* D1Sens, TacsScalar* D2Sens,
      TacsScalar* D3Sens, TacsScalar* LSens);

  /**
   * @brief Compute the buckling failure criterion
   *
   * The failure criterion is: f = N1/N1Crit + (N12/N12Crit)^2
   *
   * @param N1 Axial load
   * @param N1Crit Critical axial load
   * @param N12 Shear load
   * @param N12Crit Critical shear load
   * @return TacsScalar The failure criterion
   */
  static TacsScalar bucklingEnvelope(const TacsScalar N1,
                                     const TacsScalar N1Crit,
                                     const TacsScalar N12,
                                     const TacsScalar N12Crit);

  /**
   * @brief Compute the sensitivity of the buckling failure criterion w.r.t the
   * loads and critical loads
   *
   * @param N1 Axial load
   * @param N1Crit Critical axial load
   * @param N12 Shear load
   * @param N12Crit Critical shear load
   * @param N1Sens Sensitivity of the failure criterion w.r.t the axial load
   * @param N1CritSens Sensitivity of the failure criterion w.r.t the critical
   * axial load
   * @param N12Sens Sensitivity of the failure criterion w.r.t the shear load
   * @param N12CritSens Sensitivity of the failure criterion w.r.t the critical
   * shear load
   */
  static TacsScalar bucklingEnvelopeSens(
      const TacsScalar N1, const TacsScalar N1Crit, const TacsScalar N12,
      const TacsScalar N12Crit, TacsScalar* N1Sens, TacsScalar* N1CritSens,
      TacsScalar* N12Sens, TacsScalar* N12CritSens);

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
  static const int NUM_FAILURES = 4;  ///< Number of failure modes
};
