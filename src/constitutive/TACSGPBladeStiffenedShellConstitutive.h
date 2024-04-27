/*
====================================================================================
Blade-Stiffened Shell Constitutive Model using Machine Learning Buckling Constraints
====================================================================================
@File    :   TACSGPBladeStiffenedShellConstutive.h
@Date    :   2024/04/24
@Author  :   Sean Phillip Engelstad, Alasdair Christian Gray
@Description : Constitutive model for a blade-stiffened shell. Based on the FSDT blade models
adopted by Alasdair Christian Gray in the TACSBladeStiffenedShellConstutitutive.h class and the
original one by Graeme Kennedy. Gaussian Processes for Machine Learning are used for the buckling
constraints of the stiffened panels.
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
#include "TACSBladeStiffenedShellConstitutive.h"

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
 * - Stiffener pitch
 * - Panel thickness
 * - Panel ply fractions
 * - Stiffener height
 * - Stiffener thickness
 * - Stiffener ply fractions
 * NOTE : the panelLength design variable is removed from the TACSBladeStiffenedShellConstitutive.h model
 *
 * The failure criterion returned by this model combines numerous possible
 * failure modes into a single value. A material failure criterion (which one
 * depends on the Orthotropic ply objects you pass to this class) is computed at
 * the upper and lower surface of the panel and at the tip of the stiffener,
 * this calculation is performed for every ply angle present in the panel and
 * stiffener laminate. Additionally buckling criteria are computed for combined
 * shear and axial buckling for both a global buckling mode (i.e the entire
 * panel buckles) and a local buckling mode (i.e. the panel buckles between a
 * pair of stiffeners) and a stiffener crippling mode. These buckling failure 
 * values are aggregated along with the material failure values into a single 
 * failure value using KS aggregation. The smoothness and conservatism of 
 * this aggregation can be controlled using the `setKSWeight` method.
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
class TACSGPBladeStiffenedShellConstitutive : public TACSBladeStiffenedShellConstitutive {
 public:
  /**
   * @brief Construct a new TACSBladeStiffenedShellConstitutive object
   *
   * @param _panelPly Orthotropic ply object for the panel
   * @param _stiffenerPly Orthotropic ply object for the stiffener
   * @param _kcorr Shear correction factor
   * @param _panelLength Panel length value
   * @param _panelWidth Panel width value
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
   * @param useGPs whether to use GPs (Machine Learning) or closed-form for buckling constraints
   * height
   */
  TACSGPBladeStiffenedShellConstitutive(
      TACSOrthotropicPly* _panelPly, TACSOrthotropicPly* _stiffenerPly,
      TacsScalar _kcorr, TacsScalar _panelLength, TACSScalar _panelWidth,
      TacsScalar _stiffenerPitch, int _stiffenerPitchNum,
      TacsScalar _panelThick, int _panelThickNum, int _numPanelPlies,
      TacsScalar _panelPlyAngles[], TacsScalar _panelPlyFracs[],
      int _panelPlyFracNums[], TacsScalar _stiffenerHeight,
      int _stiffenerHeightNum, TacsScalar _stiffenerThick,
      int _stiffenerThickNum, int _numStiffenerPlies,
      TacsScalar _stiffenerPlyAngles[], TacsScalar _stiffenerPlyFracs[],
      int _stiffenerPlyFracNums[], TacsScalar _flangeFraction = 1.0, bool useGPs = false);

  ~TACSGPBladeStiffenedShellConstitutive();

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
  static bool testCriticalShearLoadSens();

 protected:

  // ==============================================================================
  // Override Failure constraint and sensitivities
  // ==============================================================================

    // Compute the failure values for each failure mode of the stiffened panel
    TacsScalar computeFailureValues(const TacsScalar e[], TacsScalar fails[]);

    // Evaluate the derivative of the failure criteria w.r.t. the strain
    TacsScalar evalFailureStrainSens(int elemIndex, const double pt[], const TacsScalar X[],
        const TacsScalar e[], TacsScalar sens[]);

    // Add the derivative of the failure criteria w.r.t. the design variables
    void addFailureDVSens(int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
        const TacsScalar strain[], int dvLen, TacsScalar dfdx[]);

  // ==============================================================================
  // Buckling functions
  // ==============================================================================

  /**
   * @brief Compute the non-dimensional affine aspect ratio rho_0 = a/b * (D22/D11)**0.25
   * 
   * @param D11 D11 stiffness
   * @param D22 D22 stiffness
   * @param a panel length
   * @param b panel width
   * 
   */
  static inline TacsScalar computeAffineAspectRatio(const TacsScalar D11,
                                                    const TacsScalar D22,
                                                    const TacsScalar a,
                                                    const TacsScalar b) {
    return a / b * pow(D22 / D11, 0.25);
  }

  /**
   * @brief Compute the derivatives non-dimensional affine aspect ratio rho_0 = a/b * (D22/D11)**0.25
   * 
   * @param D11 D11 stiffness
   * @param D22 D22 stiffness
   * @param a panel length
   * @param b panel width
   * @param D11sens Sensitivity w.r.t. the D11 stiffness
   * @param D22sens Sensitivity w.r.t. the D22 stiffness
   * @param asens Sensitivity w.r.t. the panel length a
   * @param bsens Sensitivity w.r.t. the panel width b
   * 
   */
  static TacsScalar computeAffineAspectRatioSens(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar a,
    const TacsScalar b, TacsScalar* D11sens, TacsScalar* D22sens,
    TacsScalar* asens, TacsScalar* bsens);

  /**
   * @brief Compute the non-dimensional generalized rigidity xi = (D12 + 2 D66) / sqrt(D11 * D22)
   * 
   * @param D11 D11 stiffness
   * @param D22 D22 stiffness
   * @param D12 D12 stiffness
   * @param D66 D66 stiffness
   * 
   */
  static inline TacsScalar computeGeneralizedRigidity(const TacsScalar D11,
                                                         const TacsScalar D22,
                                                         const TacsScalar D12,
                                                         const TacsScalar D66) {
    return (D12 + 2.0 * D66) / sqrt(D11 * D22);
  }

  /**
   * @brief Compute the derivatives of the non-dimensional generalized rigidity xi = (D12 + 2 D66) / sqrt(D11 * D22)
   * 
   * @param D11 D11 stiffness
   * @param D22 D22 stiffness
   * @param D12 D12 stiffness
   * @param D66 D66 stiffness
   * @param D11sens Sensitivity w.r.t. the D11 stiffness
   * @param D22sens Sensitivity w.r.t. the D22 stiffness
   * @param D12sens Sensitivity w.r.t. the D12 stiffness
   * @param D66sens Sensitivity w.r.t. the D66 stiffness
   * 
   */
  static TacsScalar computeGeneralizedRigiditySens(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar D12,
    const TacsScalar D66, TacsScalar* D11sens, TacsScalar* D22sens,
    TacsScalar* D12sens, TacsScalar* D66sens);

  /**
   * @brief Compute the non-dimensional generalized Poisson's ratio eps = 1/xi * D12 / sqrt(D11*D22) = (D12 + 2D66) / D12
   * 
   * @param D12 D12 stiffness
   * @param D66 D66 stiffness
   * 
   */
  static inline TacsScalar computeGeneralizedPoissonsRatio(const TacsScalar D12,
                                                         const TacsScalar D66) {
    return (D12 + 2.0 * D66) / D12;
  }

  /**
   * @brief Compute the derivatives of the non-dimensional generalized Poisson's ratio eps = 1/xi * D12 / sqrt(D11*D22) = (D12 + 2D66) / D12
   * 
   * @param D12 D12 stiffness
   * @param D66 D66 stiffness
   * @param D12sens Sensitivity w.r.t. the D12 stiffness
   * @param D66sens Sensitivity w.r.t. the D66 stiffness
   * 
   */
  static TacsScalar computeGeneralizedPoissonsRatioSens(
    const TacsScalar D12, const TacsScalar D66, 
    TacsScalar* D12sens, TacsScalar* D66sens);

  /**
   * @brief Compute the non-dimensional stiffener area ratio delta = E1s * As / (E1p * s_p * h)
   *    based on the DVs => stiffenerPitch, panelThick, stiffenerHeight, stiffenerThick
   * 
   */
  TacsScalar computeStiffenerAreaRatio(); 

  /**
   * @brief Compute the sensitivities of the non-dimensional stiffener area ratio delta = E1s * As / (E1p * s_p * h)
   * 
   * @param sthickSens stiffener thickness sens
   * @param sheightSens stiffener height sens
   * @param spitchSens stiffener pitch sens
   * @param pthickSens panel thickness sens
   * 
   */
  TacsScalar computeStiffenerAreaRatioSens(
    TacsScalar* sthickSens, TacsScalar* sheightSens, TacsScalar* spitchSens, TacsScalar* pthickSens);

  /**
   * @brief Compute the non-dimensional stiffener-to-panel stiffness ratio gamma = E1s * Is / (sp * D11)
   *    based on the DVs => stiffenerPitch, panelThick, stiffenerHeight, stiffenerThick
   * 
   * @param D11 the D11 stiffness of the plate
   */
  TacsScalar computeStiffenerStiffnessRatio(const TacsScalar D11); 

  /**
   * @brief Compute the sensitivities of the non-dimensional  stiffener-to-panel stiffness ratio gamma = E1s * Is / (sp * D11)
   * 
   * @param D11 the D11 stiffness of the plate
   * @param D11sens sensitivity w.r.t. the D11 stiffness
   * @param sthickSens stiffener thickness sens
   * @param sheightSens stiffener height sens
   * @param spitchSens stiffener pitch sens
   * 
   */
  TacsScalar computeStiffenerStiffnessRatioSens(
    TacsScalar D11, TacsScalar* D11sens, TacsScalar* sthickSens, 
    TacsScalar* sheightSens, TacsScalar* spitchSens);

  /**
   * @brief Compute the critical axial load for the global buckling of the
   * stiffened panel
   * @param D11 D11 stiffness of the plate
   * @param D22 D22 stiffness of the plate
   * @param b panel width
   * @param delta stiffener-to-(local panel) area ratio
   * @param rho_0 affine aspect ratio
   * @param xi generalized rigidity
   * @param gamma stiffener-to-panel 11-bending stiffness ratio
   *
   * @return TacsScalar The critical axial load for the global buckling mode
   */
  TacsScalar computeCriticalGlobalAxialLoad(const TacsScalar D11,
                                            const TacsScalar D22,
                                            const TacsScalar b,
                                            const TacsScalar delta,
                                            const TacsScalar rho_0,
                                            const TacsScalar xi,
                                            const TacsScalar gamma);

  /**
   * @brief Compute the sensitivities w.r.t. the critical axial load 
   * for the global buckling of the stiffened panel
   * @param D11 D11 stiffness of the plate
   * @param D22 D22 stiffness of the plate
   * @param b panel width
   * @param delta stiffener-to-(local panel) area ratio
   * @param rho_0 affine aspect ratio
   * @param xi generalized rigidity
   * @param gamma stiffener-to-panel 11-bending stiffness ratio
   * @param D11sens Sensitivity w.r.t. the D11 stiffness
   * @param D22sens Sensitivity w.r.t.  D22 stiffness
   * @param bsens Sensitivity w.r.t.  panel width
   * @param deltasens Sensitivity w.r.t.  stiffener-to-(local panel) area ratio
   * @param rho_0sens Sensitivity w.r.t.  affine aspect ratio
   * @param xisens Sensitivity w.r.t.  generalized rigidity
   * @param gammasens Sensitivity w.r.t. stiffener-to-panel 11-bending stiffness ratio
   *
   * @return TacsScalar The critical axial load for the global buckling mode
   */
  TacsScalar computeCriticalGlobalAxialLoadSens(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar b,
    const TacsScalar delta, const TacsScalar rho_0, const TacsScalar xi, 
    const TacsScalar gamma, TacsScalar* D11sens, TacsScalar* D22sens, 
    TacsScalar* bsens, TacsScalar* deltasens, TacsScalar* rho_0sens, 
    TacsScalar* xisens, TacsScalar* gammasens);

  /**
   * @brief Compute the critical axial load for the local buckling mode of the
   * stiffened panel
   * @param D11 D11 stiffness of the plate
   * @param D22 D22 stiffness of the plate
   * @param rho_0 affine aspect ratio
   * @param xi generalized rigidity
   *
   * @return TacsScalar The critical axial load for the local buckling mode
   */
  TacsScalar computeCriticalLocalAxialLoad(const TacsScalar D11,
                                                    const TacsScalar D22,
                                                    const TacsScalar rho_0,
                                                    const TacsScalar xi);

  /**
   * @brief Compute the sensitivities w.r.t. the critical axial load 
   * for the global buckling of the stiffened panel
   * @param D11 D11 stiffness of the plate
   * @param D22 D22 stiffness of the plate
   * @param rho_0 affine aspect ratio
   * @param xi generalized rigidity
   * @param D11sens Sensitivity w.r.t. the D11 stiffness
   * @param D22sens Sensitivity w.r.t.  D22 stiffness
   * @param rho_0sens Sensitivity w.r.t.  affine aspect ratio
   * @param xisens Sensitivity w.r.t.  generalized rigidity
   * @param spitchsens Sensitivity w.r.t. the stiffener pitch
   *
   * @return TacsScalar The critical axial load for the global buckling mode
   */
  TacsScalar computeCriticalLocalAxialLoadSens(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar rho_0, 
    const TacsScalar xi, TacsScalar* D11sens, TacsScalar* D22sens, 
    TacsScalar* rho_0sens, TacsScalar* xisens, TacsScalar* spitchsens);

  /**
   * @brief Compute the critical shear load for global buckling
   *
   *
   * input:
   * @param D11 the D11 stiffness of the plate
   * @param D22 the D22 stiffness of the plate
   * @param b the panel width
   * @param xi the generalized rigidity
   * @param gamma the stiffener-to-(local panel) bending stiffness ratio
   *
   * returns:
   * @return N12Crit  the approximate critical buckling load
   */
  TacsScalar computeCriticalGlobalShearLoad(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar b,
    const TacsScalar xi, const TacsScalar gamma,);

  /**
   * @brief Compute the sensitivity of the critical shear load for global buckling
   *
   * @param D11 the D11 stiffness of the plate
   * @param D22 the D22 stiffness of the plate
   * @param b the panel width
   * @param xi the generalized rigidity
   * @param gamma the stiffener-to-(local panel) bending stiffness ratio
   * @param D11sens sensitivity w.r.t. the D11 stiffness of the plate
   * @param D22sens sensitivity w.r.t. the D22 stiffness of the plate
   * @param bsens sensitivity w.r.t. the panel width
   * @param xisens sensitivity w.r.t. the generalized rigidity
   * @param gammasens sensitivity w.r.t. the stiffener-to-(local panel) bending stiffness ratio
   * @return TacsScalar The critical shear buckling load
   */
  TacsScalar computeCriticalGlobalShearLoadSens( 
    const TacsScalar D11, const TacsScalar D22, const TacsScalar b,
    const TacsScalar xi, const TacsScalar gamma,
    TacsScalar* D11sens, TacsScalar* D22sens, TacsScalar* bsens,
    TacsScalar* xisens, TacsScalar* gammasens);

  /**
   * @brief Compute the critical shear load for local buckling
   *
   *
   * input:
   * @param D11 the D11 stiffness of the plate
   * @param D22 the D22 stiffness of the plate
   * @param xi the generalized rigidity
   *
   * returns:
   * @return N12Crit  the approximate critical buckling load
   */
  TacsScalar computeCriticalLocalShearLoad(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar xi);

  /**
   * @brief Compute the sensitivity of the critical shear buckling load (local mode)
   *
   * @param D11 the D11 stiffness of the plate
   * @param D22 the D22 stiffness of the plate
   * @param b the panel width
   * @param xi the generalized rigidity
   * @param gamma the stiffener-to-(local panel) bending stiffness ratio
   * @param D11sens sensitivity w.r.t. the D11 stiffness of the plate
   * @param D22sens sensitivity w.r.t. the D22 stiffness of the plate
   * @param xisens sensitivity w.r.t. the generalized rigidity
   * @param spitchsens sensitivity w.r.t. the stiffener pitch
   * @return TacsScalar The critical shear buckling load
   */
  TacsScalar computeCriticalLocalShearLoadSens( 
    const TacsScalar D11, const TacsScalar D22, const TacsScalar xi, 
    TacsScalar* D11sens, TacsScalar* D22sens, TacsScalar* xisens,
    TacsScalar* spitchsens);

  /**
   * @brief Compute the non-dimensional parameters lam1bar, lam2bar using Newton iteration for the critical shear load closed-form solution.
   *
   * @param xi the non-dimensional generalized rigidity of the plate
   * @param gamma the stiffener-to-plate stiffness ratio
   * @param lam1bar return value of lam1bar
   * @param lam2bar return value of lam2bar which is always positive since it shows up as squared
   */
  static void nondimShearParams(const TacsScalar xi, const TacsScalar gamma, TacsScalar* lam1bar, TacsScalar* lam2bar);
  
    /**
   * @brief Compute the derivatives of lam2sq  for the critical shear load closed-form solution.
   *
   * @param xi the non-dimensional generalized rigidity of the plate
   * @param gamma the stiffener-to-plate stiffness ratio
   * @param lam1bar return value of lam1bar
   * @param lam2bar return value of lam2bar which is always positive since it shows up as squared
   * @param dl1xi dlam1_bar/dxi sens
   * @param dl1gamma dlam1_bar/dgamma sens
   * @param dl2xi dlam2_bar/dxi sens
   * @param dl2gamma dlam2_bar/dgamma sens
   */
  static void nondimShearParamsSens(const TacsScalar xi, const TacsScalar gamma,
  TacsScalar* lam1bar, TacsScalar* lam2bar,
  TacsScalar* dl1xi, TacsScalar* dl1gamma, TacsScalar* dl2xi, TacsScalar* dl2gamma);

  static TacsScalar lam2Constraint(const TacsScalar lam2sq, const TacsScalar xi, const TacsScalar gamma);
  static TacsScalar lam2ConstraintDeriv(const TacsScalar lam2sq, const TacsScalar xi, const TacsScalar gamma);

    /**
   * @brief Compute the critical stiffener crippling load
   *
   *
   * input:
   * @param D11 the D11 stiffness of the stiffener
   * @param D22 the D22 stiffness of the stiffener
   * @param xi the generalized rigidity
   * @param genPoisss the generalized Poisson's ratio of the stiffener
   *
   * returns:
   * @return N11crit  the approximate critical crippling in-plane load of the stiffener
   */
  TacsScalar computeStiffenerCripplingLoad(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar xi, const TacsScalar genPoiss);

  /**
   * @brief Compute the sensitivity of the critical stiffener crippling load
   *
   * @param D11 the D11 stiffness of the stiffener
   * @param D22 the D22 stiffness of the stiffener
   * @param xi the generalized rigidity
   * @param genPoisss the generalized Poisson's ratio of the stiffener
   * @param D11sens Sensitivity of the the D11 stiffness of the stiffener
   * @param D22sens Sensitivity of the the D22 stiffness of the stiffener
   * @param xisens Sensitivity of the the generalized rigidity
   * @param genPoiss_sens Sensitivity of the the generalized Poisson's ratio of the stiffener
   * @param sheightsens Sensitivity w.r.t. the stiffener height DV
   * @return N11crit  the approximate critical crippling in-plane load of the stiffener
   */
  TacsScalar computeStiffenerCripplingLoadSens( 
    const TacsScalar D11, const TacsScalar D22, const TacsScalar xi, 
    const TacsScalar genPoiss, TacsScalar* D11sens, TacsScalar* D22sens,
    TacsScalar* xisens, TacsScalar* genPoiss_sens, TacsScalar* sheightsens);

  // ==============================================================================
  // Attributes
  // ==============================================================================

  // --- Design variable values ---
  TacsScalar panelWidth;          ///< Panel width
  bool useGPs;
  int NUM_CF_MODES = 50; // number of modes used in closed-form method
};
