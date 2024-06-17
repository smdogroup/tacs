/*
====================================================================================
Blade-Stiffened Shell Constitutive Model using Machine Learning Buckling
Constraints
====================================================================================
@File    :   TACSGPBladeStiffenedShellConstutive.h
@Date    :   2024/04/24
@Author  :   Sean Phillip Engelstad, Alasdair Christian Gray
@Description : Constitutive model for a blade-stiffened shell. Based on the FSDT
blade models adopted by Alasdair Christian Gray in the
TACSBladeStiffenedShellConstutitutive.h class and the original one by Graeme
Kennedy. Gaussian Processes for Machine Learning are used for the buckling
constraints of the stiffened panels.
*/

#pragma once

// =============================================================================
// Standard Library Includes
// =============================================================================

// =============================================================================
// Extension Includes
// =============================================================================
#include "TACSBladeStiffenedShellConstitutive.h"
#include "TACSGaussianProcessModel.h"
#include "TACSPanelGPs.h"
#include "TacsUtilities.h"

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
 * - Panel Length
 * - Stiffener pitch
 * - Panel thickness
 * - Panel ply fractions
 * - Stiffener height
 * - Stiffener thickness
 * - Stiffener ply fractions
 * - Panel Width
 * NOTE : the panelLength design variable is removed from the
 * TACSBladeStiffenedShellConstitutive.h model
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
class TACSGPBladeStiffenedShellConstitutive
    : public TACSBladeStiffenedShellConstitutive {
 public:
  /**
   * @brief Construct a new TACSBladeStiffenedShellConstitutive object
   *
   * @param panelPly Orthotropic ply object for the panel
   * @param stiffenerPly Orthotropic ply object for the stiffener
   * @param kcorr Shear correction factor
   * @param panelLength Panel length value
   * @param panelLengthNum Panel Length design variable number
   * @param stiffenerPitch Stiffener pitch value
   * @param stiffenerPitchNum Stiffener pitch design variable number
   * @param panelThick Panel thickness value
   * @param panelThickNum Panel thickness design variable number
   * @param numPanelPlies Number of ply angles in the panel laminate
   * @param panelPlyAngles Panel ply angles
   * @param panelPlyFracs Panel ply fractions
   * @param panelPlyFracNums Panel ply fraction design variable numbers
   * @param stiffenerHeight Stiffener height value
   * @param stiffenerHeightNum Stiffener height design variable number
   * @param stiffenerThick Stiffener thickness value
   * @param stiffenerThickNum Stiffener thickness design variable number
   * @param numStiffenerPlies Number of ply angles in the stiffener laminate
   * @param stiffenerPlyAngles Stiffener ply angles
   * @param stiffenerPlyFracs Stiffener ply fractions
   * @param stiffenerPlyFracNums Stiffener ply fraction design variable numbers
   * @param panelWidth Panel Width value
   * @param panelWidthNum Panel Width design variable number
   * @param flangeFraction Stiffener base width as a fraction of the stiffener
   * @param panelGPs PanelGP object (one for each TACS component) that contains
   * the GP objects or None if using CF
   */
  TACSGPBladeStiffenedShellConstitutive(
      TACSOrthotropicPly* panelPly, TACSOrthotropicPly* stiffenerPly,
      TacsScalar kcorr, TacsScalar panelLength, int panelLengthNum,
      TacsScalar stiffenerPitch, int stiffenerPitchNum, TacsScalar panelThick,
      int panelThickNum, int numPanelPlies, TacsScalar panelPlyAngles[],
      TacsScalar panelPlyFracs[], int panelPlyFracNums[],
      TacsScalar stiffenerHeight, int stiffenerHeightNum,
      TacsScalar stiffenerThick, int stiffenerThickNum, int numStiffenerPlies,
      TacsScalar stiffenerPlyAngles[], TacsScalar stiffenerPlyFracs[],
      int stiffenerPlyFracNums[], TacsScalar panelWidth, int panelWidthNum,
      TacsScalar flangeFraction = 1.0, TACSPanelGPs* panelGPs = nullptr);

  ~TACSGPBladeStiffenedShellConstitutive();

  // ==============================================================================
  // Tests
  // ==============================================================================

  /**
   *
   * @brief Test the nondimensional parameter computations against finite
   * difference.
   */
  TacsScalar testNondimensionalParameters(TacsScalar epsilon, int printLevel);

  /**
   *
   * @brief Test the axial critical loads.
   */
  TacsScalar testAxialCriticalLoads(TacsScalar epsilon, int printLevel);

  /**
   *
   * @brief Test the shear critical loads.
   */
  TacsScalar testShearCriticalLoads(TacsScalar epsilon, int printLevel);

  /**
   *
   * @brief Test the crippling critical loads.
   */
  TacsScalar testStiffenerCripplingLoad(TacsScalar epsilon, int printLevel);

  /**
   *
   * @brief Test all GP tests
   */
  TacsScalar testAllTests(TacsScalar epsilon, int printLevel);

  // ==============================================================================
  // Getter and setters
  // ==============================================================================

  // get the three Gaussian Process model pointers
  TACSAxialGaussianProcessModel* getAxialGP() {
    if (this->panelGPs) {
      return this->panelGPs->getAxialGP();
    } else {
      return nullptr;
    }
  }
  TACSShearGaussianProcessModel* getShearGP() {
    if (this->panelGPs) {
      return this->panelGPs->getShearGP();
    } else {
      return nullptr;
    }
  }
  TACSCripplingGaussianProcessModel* getCripplingGP() {
    if (this->panelGPs) {
      return this->panelGPs->getCripplingGP();
    } else {
      return nullptr;
    }
  }

  // Retrieve the global design variable numbers
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);

  // Set the element design variable from the design vector
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);

  // Get the element design variables values
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);

  // Get the lower and upper bounds for the design variable values
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]);

  // Retrieve the design variable for plotting purposes
  TacsScalar evalDesignFieldValue(int elemIndex, const double pt[],
                                  const TacsScalar X[], int index);

  // set the KS weight for the failure constraints and the GP models (if GP
  // models are active)
  void setKSWeight(double ksWeight);

  // set the closed-form shear mode [1 - regular, 2 - analytic surogate]
  void setCFShearMode(int newMode) { CFshearMode = newMode; }

  // set the DV write out modes for f5 files [0 - regular DVs, 1 - ND, 2 - fail
  // indexes]
  void setWriteDVMode(int newMode) { writeDVmode = newMode; }

  // ==============================================================================
  // Verification purpose public routines
  // ==============================================================================

  /**
   * @brief Compute the non-dimensional critical axial load for the global
   * buckling of the stiffened panel (for output purposes only)
   * @param rho_0 affine aspect ratio
   * @param xi generalized rigidity
   * @param gamma stiffener-to-panel 11-bending stiffness ratio
   * @param zeta the transverse shear stiffness parameter
   *
   * @return TacsScalar The nondimensional critical axial load for the global
   * buckling mode
   */
  TacsScalar nondimCriticalGlobalAxialLoad(const TacsScalar rho_0,
                                           const TacsScalar xi,
                                           const TacsScalar gamma,
                                           const TacsScalar zeta);

  /**
   * @brief Compute the non-dimensional critical axial load for the local
   * buckling of the stiffened panel (for output purposes only)
   * @param rho_0 affine aspect ratio
   * @param xi generalized rigidity
   * @param zeta the transverse shear stiffness parameter
   *
   * @return TacsScalar The nondimensional critical axial load for the local
   * buckling mode
   */
  TacsScalar nondimCriticalLocalAxialLoad(const TacsScalar rho_0,
                                          const TacsScalar xi,
                                          const TacsScalar zeta);

  /**
   * @brief Compute the non-dimensional critical shear load for the global
   * buckling of the stiffened panel (for output purposes only)
   * @param rho_0 affine aspect ratio
   * @param xi generalized rigidity
   * @param gamma stiffener-to-panel 11-bending stiffness ratio
   * @param zeta the transverse shear stiffness parameter
   *
   * @return TacsScalar The nondimensional critical shear load for the global
   * buckling mode
   */
  TacsScalar nondimCriticalGlobalShearLoad(const TacsScalar rho_0,
                                           const TacsScalar xi,
                                           const TacsScalar gamma,
                                           const TacsScalar zeta);

  /**
   * @brief Compute the non-dimensional critical shear load for the local
   * buckling of the stiffened panel (for output purposes only)
   * @param rho_0 affine aspect ratio
   * @param xi generalized rigidity
   * @param zeta the transverse shear stiffness parameter
   *
   * @return TacsScalar The nondimensional critical shear load for the local
   * buckling mode
   */
  TacsScalar nondimCriticalLocalShearLoad(const TacsScalar rho_0,
                                          const TacsScalar xi,
                                          const TacsScalar zeta);

  /**
   * @brief Compute the non-dimensional critical shear load for the local
   * buckling of the stiffened panel (for output purposes only)
   * @param rho_0 affine aspect ratio
   * @param xi generalized rigidity
   * @param genPoiss generalized poisson's ratio
   * @param zeta the transverse shear stiffness parameter
   *
   * @return TacsScalar The nondimensional critical shear load for the local
   * buckling mode
   */
  TacsScalar nondimStiffenerCripplingLoad(const TacsScalar rho_0,
                                          const TacsScalar xi,
                                          const TacsScalar genPoiss,
                                          const TacsScalar zeta);

 protected:
  // ==============================================================================
  // Override Failure constraint and sensitivities
  // ==============================================================================

  // override from superclass so that evalFailure uses the computeFailureValues
  // from this subroutine
  TacsScalar evalFailure(int elemIndex, const double pt[], const TacsScalar X[],
                         const TacsScalar e[]);

  // Compute the failure values for each failure mode of the stiffened panel
  TacsScalar computeFailureValues(const TacsScalar e[], TacsScalar fails[]);

  // Evaluate the derivative of the failure criteria w.r.t. the strain
  TacsScalar evalFailureStrainSens(int elemIndex, const double pt[],
                                   const TacsScalar X[], const TacsScalar e[],
                                   TacsScalar sens[]);

  // Add the derivative of the failure criteria w.r.t. the design variables
  void addFailureDVSens(int elemIndex, TacsScalar scale, const double pt[],
                        const TacsScalar X[], const TacsScalar strain[],
                        int dvLen, TacsScalar dfdx[]);

  // ==============================================================================
  // Stiffener crippling helper functions
  // ==============================================================================

  TacsScalar computeStiffenerInPlaneLoad(const TacsScalar stiffenerStrain[],
                                         TacsScalar* A11s);

  TacsScalar computeStiffenerInPlaneLoadSens(const TacsScalar scale,
                                             const TacsScalar panelStrain[],
                                             const TacsScalar stiffenerStrain[],
                                             const TacsScalar dN11_stiff,
                                             TacsScalar dfdx[]);

  void computeStiffenerCripplingStiffness(TacsScalar C[]);

  // ==============================================================================
  // Buckling functions
  // ==============================================================================

  /**
   * @brief Compute the non-dimensional affine aspect ratio rho_0 = a/b *
   * (D22/D11)**0.25
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
   * @brief Compute the derivatives non-dimensional affine aspect ratio rho_0 =
   * a/b * (D22/D11)**0.25
   *
   * @param rho0sens backpropagated sensitivity for rho0 non-dimensional
   * parameter
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
  TacsScalar computeAffineAspectRatioSens(
      const TacsScalar rho0sens, const TacsScalar D11, const TacsScalar D22,
      const TacsScalar a, const TacsScalar b, TacsScalar* D11sens,
      TacsScalar* D22sens, TacsScalar* asens, TacsScalar* bsens);

  /**
   *
   * @brief Test the affine aspect ratio sensitivity.
   */
  TacsScalar testAffineAspectRatio(const TacsScalar epsilon, int printLevel);

  /**
   * @brief Compute the non-dimensional generalized rigidity xi = (D12 + 2 D66)
   * / sqrt(D11 * D22)
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
   * @brief Compute the derivatives of the non-dimensional generalized rigidity
   * xi = (D12 + 2 D66) / sqrt(D11 * D22)
   *
   * @param xisens backpropagated sensitivity w.r.t. xi
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
  TacsScalar computeGeneralizedRigiditySens(
      const TacsScalar xisens, const TacsScalar D11, const TacsScalar D22,
      const TacsScalar D12, const TacsScalar D66, TacsScalar* D11sens,
      TacsScalar* D22sens, TacsScalar* D12sens, TacsScalar* D66sens);

  /**
   *
   * @brief Test the generalized rigidity sensitivity.
   */
  TacsScalar testGeneralizedRigidity(const TacsScalar epsilon, int printLevel);

  /**
   * @brief Compute the non-dimensional generalized Poisson's ratio eps = 1/xi *
   * D12 / sqrt(D11*D22) = (D12 + 2D66) / D12
   *
   * @param D12 D12 stiffness
   * @param D66 D66 stiffness
   *
   */
  static inline TacsScalar computeGeneralizedPoissonsRatio(
      const TacsScalar D12, const TacsScalar D66) {
    return D12 / (D12 + 2.0 * D66);
  }

  /**
   * @brief Compute the derivatives of the non-dimensional generalized Poisson's
   * ratio eps = 1/xi * D12 / sqrt(D11*D22) = (D12 + 2D66) / D12
   *
   * @param epssens backpropagated sensitivity for the output gen eps
   * @param D12 D12 stiffness
   * @param D66 D66 stiffness
   * @param D12sens Sensitivity w.r.t. the D12 stiffness
   * @param D66sens Sensitivity w.r.t. the D66 stiffness
   *
   */
  TacsScalar computeGeneralizedPoissonsRatioSens(const TacsScalar epssens,
                                                 const TacsScalar D12,
                                                 const TacsScalar D66,
                                                 TacsScalar* D12sens,
                                                 TacsScalar* D66sens);

  /**
   *
   * @brief Test the generalized poisson's ratio sensitivity.
   */
  TacsScalar testGeneralizedPoissonsRatio(const TacsScalar epsilon,
                                          int printLevel);

  /**
   * @brief Compute the non-dimensional stiffener area ratio delta = E1s * As /
   * (E1p * s_p * h) based on the DVs => stiffenerPitch, panelThick,
   * stiffenerHeight, stiffenerThick
   *
   */
  TacsScalar computeStiffenerAreaRatio();

  /**
   * @brief Compute the sensitivities of the non-dimensional stiffener area
   * ratio delta = E1s * As / (E1p * s_p * h)
   *
   * @param deltasens backpropagated sensitivity w.r.t the stiffener area ratio
   * delta
   * @param sthickSens stiffener thickness sens
   * @param sheightSens stiffener height sens
   * @param spitchSens stiffener pitch sens
   * @param pthickSens panel thickness sens
   *
   */
  TacsScalar computeStiffenerAreaRatioSens(const TacsScalar deltasens,
                                           TacsScalar* sthickSens,
                                           TacsScalar* sheightSens,
                                           TacsScalar* spitchSens,
                                           TacsScalar* pthickSens);

  /**
   *
   * @brief Test the stiffener area ratio sensitivity aka delta parameter
   */
  TacsScalar testStiffenerAreaRatio(const TacsScalar epsilon, int printLevel);

  /**
   * @brief Compute the non-dimensional stiffener-to-panel stiffness ratio gamma
   * = E1s * Is / (sp * D11) based on the DVs => stiffenerPitch, panelThick,
   * stiffenerHeight, stiffenerThick
   *
   * @param D11 the D11 stiffness of the plate
   */
  TacsScalar computeStiffenerStiffnessRatio(const TacsScalar D11);

  /**
   *
   * @brief Return the bending stiffness of an individual stiffener
   *
   **/
  TacsScalar computeStiffenerBendingStiffness();

  /**
   *
   * @brief Returns the 1x2 Jacobian of the stiffener bending stiffness
   *computation in 2 separate scalars
   *
   **/
  TacsScalar computeStiffenerBendingStiffnessSens(TacsScalar& sthickSens,
                                                  TacsScalar& sheightSens);

  /**
   * @brief Compute the sensitivities of the non-dimensional  stiffener-to-panel
   * stiffness ratio gamma = E1s * Is / (sp * D11)
   *
   * @param gammasens backpropagated derivative through gamma output
   * @param D11 the D11 stiffness of the plate
   * @param D11sens sensitivity w.r.t. the D11 stiffness
   * @param sthickSens stiffener thickness sens
   * @param sheightSens stiffener height sens
   * @param spitchSens stiffener pitch sens
   *
   */
  TacsScalar computeStiffenerStiffnessRatioSens(
      const TacsScalar gammasens, const TacsScalar D11, TacsScalar* D11sens,
      TacsScalar* sthickSens, TacsScalar* sheightSens, TacsScalar* spitchSens);

  /**
   *
   * @brief Test the stiffener stiffness ratio sensitivity aka gamma parameter
   */
  TacsScalar testStiffenerStiffnessRatio(const TacsScalar epsilon,
                                         int printLevel);

  /**
   * @brief Compute the non-dimensional transverse shear parameter
   * zeta = A11 / A66 * (h/b)**2
   *
   * @param A66 the A66 stiffness of the plate
   * @param A11 the A11 stiffness of the plate
   * @param b the panel widtdh
   * @param h the panel height
   *
   */
  TacsScalar computeTransverseShearParameter(TacsScalar A66, TacsScalar A11,
                                             TacsScalar b, TacsScalar h);

  /**
   * @brief Compute the sensitivities of the non-dimensional transverse shear
   * parameter zeta = A11 / A66 * (h/b)**2
   *
   * @param zetasens backpropagated sensitivity for output zeta
   * @param A66 the A66 stiffness of the plate
   * @param A11 the A11 stiffness of the plate
   * @param b the panel widtdh
   * @param h the panel height
   * @param A66sens sensitivity w.r.t the A66 stiffness of the plate
   * @param A11sens sensitivity w.r.t the A11 stiffness of the plate
   * @param bsens sensitivity w.r.t the panel widtdh
   * @param hsens sensitivity w.r.t the panel height
   *
   */
  TacsScalar computeTransverseShearParameterSens(
      const TacsScalar zetasens, const TacsScalar A66, const TacsScalar A11,
      const TacsScalar b, const TacsScalar h, TacsScalar* A66sens,
      TacsScalar* A11sens, TacsScalar* bsens, TacsScalar* hsens);

  /**
   *
   * @brief Test the transverse shear parameter sensitivity aka zeta parameter
   */
  TacsScalar testTransverseShearParameter(const TacsScalar epsilon,
                                          int printLevel);

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
   * @param zeta the transverse shear stiffness parameter
   *
   * @return TacsScalar The critical axial load for the global buckling mode
   */
  TacsScalar computeCriticalGlobalAxialLoad(
      const TacsScalar D11, const TacsScalar D22, const TacsScalar b,
      const TacsScalar delta, const TacsScalar rho_0, const TacsScalar xi,
      const TacsScalar gamma, const TacsScalar zeta);

  /**
   * @brief Compute the sensitivities w.r.t. the critical axial load
   * for the global buckling of the stiffened panel
   * @param N1sens backpropagated derivative through N1crit computation of
   * output
   * @param D11 D11 stiffness of the plate
   * @param D22 D22 stiffness of the plate
   * @param b panel width
   * @param delta stiffener-to-(local panel) area ratio
   * @param rho_0 affine aspect ratio
   * @param xi generalized rigidity
   * @param gamma stiffener-to-panel 11-bending stiffness ratio
   * @param zeta the transverse shear stiffness parameter
   * @param D11sens Sensitivity w.r.t. the D11 stiffness
   * @param D22sens Sensitivity w.r.t.  D22 stiffness
   * @param bsens Sensitivity w.r.t.  panel width
   * @param deltasens Sensitivity w.r.t.  stiffener-to-(local panel) area ratio
   * @param rho_0sens Sensitivity w.r.t.  affine aspect ratio
   * @param xisens Sensitivity w.r.t.  generalized rigidity
   * @param gammasens Sensitivity w.r.t. stiffener-to-panel 11-bending stiffness
   * @param zetasens Sensitivity w.r.t.the transverse shear stiffness parameter
   * ratio
   *
   * @return TacsScalar The critical axial load for the global buckling mode
   */
  TacsScalar computeCriticalGlobalAxialLoadSens(
      const TacsScalar N1sens, const TacsScalar D11, const TacsScalar D22,
      const TacsScalar b, const TacsScalar delta, const TacsScalar rho_0,
      const TacsScalar xi, const TacsScalar gamma, const TacsScalar zeta,
      TacsScalar* D11sens, TacsScalar* D22sens, TacsScalar* bsens,
      TacsScalar* deltasens, TacsScalar* rho_0sens, TacsScalar* xisens,
      TacsScalar* gammasens, TacsScalar* zetasens);

  /**
   *
   * @brief Test the critical global axial load function
   */
  TacsScalar testCriticalGlobalAxialLoad(const TacsScalar epsilon,
                                         int printLevel);

  /**
   * @brief Compute the critical axial load for the local buckling mode of the
   * stiffened panel
   * @param D11 D11 stiffness of the plate
   * @param D22 D22 stiffness of the plate
   * @param rho_0 affine aspect ratio
   * @param xi generalized rigidity
   * @param zeta Transverse shear stiffness parameter
   *
   * @return TacsScalar The critical axial load for the local buckling mode
   */
  TacsScalar computeCriticalLocalAxialLoad(const TacsScalar D11,
                                           const TacsScalar D22,
                                           const TacsScalar rho_0,
                                           const TacsScalar xi,
                                           const TacsScalar zeta);

  /**
   * @brief Compute the sensitivities w.r.t. the critical axial load
   * for the global buckling of the stiffened panel
   * @param N1sens backpropagated derivative through N1crit computation of local
   * mode
   * @param D11 D11 stiffness of the plate
   * @param D22 D22 stiffness of the plate
   * @param rho_0 affine aspect ratio
   * @param xi generalized rigidity
   * @param zeta Transverse shear stiffness parameter
   * @param D11sens Sensitivity w.r.t. the D11 stiffness
   * @param D22sens Sensitivity w.r.t.  D22 stiffness
   * @param rho_0sens Sensitivity w.r.t.  affine aspect ratio
   * @param xisens Sensitivity w.r.t.  generalized rigidity
   * @param spitchsens Sensitivity w.r.t. the stiffener pitch
   * @param zetasens Sensitivity w.r.t. the transverse shear stiffness parameter
   *
   * @return TacsScalar The critical axial load for the global buckling mode
   */
  TacsScalar computeCriticalLocalAxialLoadSens(
      const TacsScalar N1sens, const TacsScalar D11, const TacsScalar D22,
      const TacsScalar rho_0, const TacsScalar xi, const TacsScalar zeta,
      TacsScalar* D11sens, TacsScalar* D22sens, TacsScalar* rho_0sens,
      TacsScalar* xisens, TacsScalar* spitchsens, TacsScalar* zetasens);

  /**
   *
   * @brief Test the critical local axial load function
   */
  TacsScalar testCriticalLocalAxialLoad(const TacsScalar epsilon,
                                        int printLevel);

  /**
   * @brief Compute the critical shear load for global buckling
   *
   *
   * input:
   * @param D11 the D11 stiffness of the plate
   * @param D22 the D22 stiffness of the plate
   * @param b the panel width
   * @param xi the generalized rigidity
   * @param rho_0 the affine aspect ratio
   * @param gamma the stiffener-to-(local panel) bending stiffness ratio
   * @param zeta the transverse shear stiffness parameter
   *
   * returns:
   * @return N12Crit  the approximate critical buckling load
   */
  TacsScalar computeCriticalGlobalShearLoad(
      const TacsScalar D11, const TacsScalar D22, const TacsScalar b,
      const TacsScalar xi, const TacsScalar rho_0, const TacsScalar gamma,
      const TacsScalar zeta);

  /**
   * @brief Compute the sensitivity of the critical shear load for global
   * buckling
   *
   * @param N12sens the backpropagated sensitivity of the output N12crit
   * @param D11 the D11 stiffness of the plate
   * @param D22 the D22 stiffness of the plate
   * @param b the panel width
   * @param xi the generalized rigidity
   * @param rho_0 the generalized affine aspect ratio
   * @param gamma the stiffener-to-(local panel) bending stiffness ratio
   * @param zeta the transverse shear stiffness parameter
   * @param D11sens sensitivity w.r.t. the D11 stiffness of the plate
   * @param D22sens sensitivity w.r.t. the D22 stiffness of the plate
   * @param bsens sensitivity w.r.t. the panel width
   * @param xisens sensitivity w.r.t. the generalized rigidity
   * @param rho_0sens sensitivity w.r.t. the affine aspect ratio
   * @param gammasens sensitivity w.r.t. the stiffener-to-(local panel) bending
   * @param zetasens sensitivity w.r.t. the transverse shear stiffness parameter
   * stiffness ratio
   * @return TacsScalar The critical shear buckling load
   */
  TacsScalar computeCriticalGlobalShearLoadSens(
      const TacsScalar N12sens, const TacsScalar D11, const TacsScalar D22,
      const TacsScalar b, const TacsScalar xi, const TacsScalar rho_0,
      const TacsScalar gamma, const TacsScalar zeta, TacsScalar* D11sens,
      TacsScalar* D22sens, TacsScalar* bsens, TacsScalar* xisens,
      TacsScalar* rho_0sens, TacsScalar* gammasens, TacsScalar* zetasens);

  /**
   *
   * @brief Test the critical global shear load function
   */
  TacsScalar testCriticalGlobalShearLoad(const TacsScalar epsilon,
                                         int printLevel);

  /**
   * @brief Compute the critical shear load for local buckling
   *
   *
   * input:
   * @param D11 the D11 stiffness of the plate
   * @param D22 the D22 stiffness of the plate
   * @param xi the generalized rigidity
   * @param rho_0 the affine aspect ratio
   * @param zeta the transverse shear stiffness parameter
   *
   * returns:
   * @return N12Crit  the approximate critical buckling load
   */
  TacsScalar computeCriticalLocalShearLoad(const TacsScalar D11,
                                           const TacsScalar D22,
                                           const TacsScalar xi,
                                           const TacsScalar rho_0,
                                           const TacsScalar zeta);

  /**
   * @brief Compute the sensitivity of the critical shear buckling load (local
   * mode)
   *
   * @param N12sens backpropagated sensitivity w.r.t the output N12crit
   * computation
   * @param D11 the D11 stiffness of the plate
   * @param D22 the D22 stiffness of the plate
   * @param b the panel width
   * @param xi the generalized rigidity
   * @param gamma the stiffener-to-(local panel) bending stiffness ratio
   * @param D11sens sensitivity w.r.t. the D11 stiffness of the plate
   * @param D22sens sensitivity w.r.t. the D22 stiffness of the plate
   * @param spitchsens sensitivity w.r.t. the stiffener pitch
   * @param xisens sensitivity w.r.t. the generalized rigidity
   * @param rho_0sens sensitivity w.r.t. the affine aspect ratio
   * @param zeta sensitivity w.r.t. the transverse shear stiffness parameter
   * @return TacsScalar The critical shear buckling load
   */
  TacsScalar computeCriticalLocalShearLoadSens(
      const TacsScalar N12sens, const TacsScalar D11, const TacsScalar D22,
      const TacsScalar xi, const TacsScalar rho_0, const TacsScalar zeta,
      TacsScalar* D11sens, TacsScalar* D22sens, TacsScalar* spitchsens,
      TacsScalar* xisens, TacsScalar* rho_0sens, TacsScalar* zetasens);

  /**
   *
   * @brief Test the critical local shear load function
   */
  TacsScalar testCriticalLocalShearLoad(const TacsScalar epsilon,
                                        int printLevel);

  /**
   * @brief Compute the non-dimensional parameters lam1bar, lam2bar using Newton
   * iteration for the critical shear load closed-form solution.
   *
   * @param xi the non-dimensional generalized rigidity of the plate
   * @param gamma the stiffener-to-plate stiffness ratio
   * @param lam1bar return value of lam1bar
   * @param lam2bar return value of lam2bar which is always positive since it
   * shows up as squared
   */
  void nondimShearParams(const TacsScalar xi, const TacsScalar gamma,
                         TacsScalar* lam1bar, TacsScalar* lam2bar);

  /**
   * @brief Compute the derivatives of lam2sq  for the critical shear load
   * closed-form solution.
   *
   * @param xi the non-dimensional generalized rigidity of the plate
   * @param gamma the stiffener-to-plate stiffness ratio
   * @param lam1bar return value of lam1bar
   * @param lam2bar return value of lam2bar which is always positive since it
   * shows up as squared
   * @param dl1xi dlam1_bar/dxi sens
   * @param dl1gamma dlam1_bar/dgamma sens
   * @param dl2xi dlam2_bar/dxi sens
   * @param dl2gamma dlam2_bar/dgamma sens
   */
  void nondimShearParamsSens(const TacsScalar xi, const TacsScalar gamma,
                             TacsScalar* lam1bar, TacsScalar* lam2bar,
                             TacsScalar* dl1xi, TacsScalar* dl1gamma,
                             TacsScalar* dl2xi, TacsScalar* dl2gamma);

  /**
   *
   * @brief Test the nondimShearParams function
   */
  TacsScalar testNondimShearParams(const TacsScalar epsilon, int printLevel);

  TacsScalar lam2Constraint(const TacsScalar lam2sq, const TacsScalar xi,
                            const TacsScalar gamma);
  TacsScalar lam2ConstraintDeriv(const TacsScalar lam2sq, const TacsScalar xi,
                                 const TacsScalar gamma);

  /**
   *
   * @brief Test the lam2Constraint function
   */
  TacsScalar testLam2Constraint(const TacsScalar epsilon, int printLevel);

  /**
   * @brief Compute the critical stiffener crippling load
   *
   *
   * input:
   * @param D11 the D11 stiffness of the stiffener
   * @param D22 the D22 stiffness of the stiffener
   * @param xi the generalized rigidity
   * @param rho_0 the affine aspect ratio
   * @param genPoisss the generalized Poisson's ratio of the stiffener
   * @param zeta the transverse shear stiffness parameter
   *
   * returns:
   * @return N11crit  the approximate critical crippling in-plane load of the
   * stiffener
   */
  TacsScalar computeStiffenerCripplingLoad(
      const TacsScalar D11, const TacsScalar D22, const TacsScalar xi,
      const TacsScalar rho_0, const TacsScalar genPoiss, const TacsScalar zeta);

  /**
   * @brief Compute the sensitivity of the critical stiffener crippling load
   *
   * @param N1sens backpropagated sensitivity w.r.t. the stiffener crit load
   * N1crit
   * @param D11 the D11 stiffness of the stiffener
   * @param D22 the D22 stiffness of the stiffener
   * @param xi the generalized rigidity
   * @param rho_0 the affine aspect ratio
   * @param genPoisss the generalized Poisson's ratio of the stiffener
   * @param zeta the transverse shear stiffness parameter
   * @param D11sens Sensitivity of the the D11 stiffness of the stiffener
   * @param D22sens Sensitivity of the the D22 stiffness of the stiffener
   * @param sheightsens Sensitivity w.r.t. the stiffener height DV
   * @param xisens Sensitivity of the the generalized rigidity
   * @param rho_0sens Sensitivity w.r.t the affine aspect ratio
   * @param genPoiss_sens Sensitivity of the the generalized Poisson's ratio of
   * the stiffener
   * @param zetasens Sensitivity w.r.t the transverse shear stiffness parameter
   * @return N11crit  the approximate critical crippling in-plane load of the
   * stiffener
   */
  TacsScalar computeStiffenerCripplingLoadSens(
      const TacsScalar N1sens, const TacsScalar D11, const TacsScalar D22,
      const TacsScalar xi, const TacsScalar rho_0, const TacsScalar genPoiss,
      const TacsScalar zeta, TacsScalar* D11sens, TacsScalar* D22sens,
      TacsScalar* sheightsens, TacsScalar* xisens, TacsScalar* rho_0sens,
      TacsScalar* genPoiss_sens, TacsScalar* zetasens);

  // ==============================================================================
  // Attributes
  // ==============================================================================

  // Machine learning Gaussian Process models stored in panel GP class
  TACSPanelGPs* panelGPs;

  // --- Design variable values ---
  TacsScalar panelWidth;  ///< Panel width
  int NUM_CF_MODES = 50;  // number of modes used in closed-form method

  TacsScalar panelWidthLowerBound = 0.0;
  TacsScalar panelWidthUpperBound = 1e20;

  // --- Design variable numbers ---
  int panelWidthNum;       ///< Panel width DV number
  int panelWidthLocalNum;  ///< Panel width local DV number

  // overwrite number of failure modes
  // panel stress, stiffener stress, global buckling, local buckling, stiffener
  // crippling
  static const int NUM_FAILURES = 5;  ///< Number of failure modes

  // different mode for computing the shear closed form solution
  int CFshearMode = 2;

  // debugging modes
  // should all be false if not debugging
  int writeDVmode = 0;  // 0 - normal DVs, 1 - NDparams, 2 - failure indexes

 private:
  // private so that subclass constName for GP buckling constraints doesn't
  // conflict with superclass
  static const char* constName;  ///< Constitutive model name
};
