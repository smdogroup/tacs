/*
====================================================================================
Blade-Stiffened Shell Constitutive Model using Machine Learning Buckling
Constraints
====================================================================================
@File    :   TACSGPBladeStiffenedShellConstutive.h
@Date    :   2024/04/24
@Author  :   Sean Phillip Engelstad
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
 * this class inherits from TACSBladeStiffenedShellConstitutive.cpp, and
 * overwrites some buckling predictions using machine learning or more accurate
 * closed-form solutions for global and local buckling.
 *
 * The TACSPanelGPs object is input into this class for each panel (see
 * smdogroup/ml_buckling repo on examples from the python side), which is useful
 * in saving and reloading buckling predictions for each element in the same
 * panel / TACS component (to speedup computation). To use Machine Learning
 * buckling predictions with a Gaussian Process model for global and local
 * buckling, the user should create TACSAxialGaussianProcessModel and
 * TACSShearGaussianProcessModel objects in the TACSPanelGPs. If the
 * TACSAxialGaussianProcessModel or TACSShearGaussianProcessModel objects are
 * empty or the TACSPanelGPs is empty, then the axial or shear buckling loads
 * will use closed-form solution buckling predictions instead. Note that either
 * the Axial GP or the shear GP can each be provided or not. Data for the
 * trained GP models is located on the smdogroup/ml_buckling github repo.
 *
 * The main differences in this class and its superclass
 * TACSBladeStiffenedShellConstitutive.cpp for the constructor are the new
 * inputs for panelWidth and panelWidthNum, an additional flag for
 * CPTstiffenerCrippling, and the PanelGPs object which holds machine learning
 * objects.
 *
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
   * @param CPTstiffenerCrippling whether to use CPT solution (true) or to use
   * experimental crippling solution (false)
   * @param panelGPs PanelGP object (one for each TACS component) that contains
   * the GP objects or None if using CF
   */
  TACSGPBladeStiffenedShellConstitutive(
      TACSOrthotropicPly *panelPly, TACSOrthotropicPly *stiffenerPly,
      TacsScalar kcorr, TacsScalar panelLength, int panelLengthNum,
      TacsScalar stiffenerPitch, int stiffenerPitchNum, TacsScalar panelThick,
      int panelThickNum, int numPanelPlies, TacsScalar panelPlyAngles[],
      TacsScalar panelPlyFracs[], int panelPlyFracNums[],
      TacsScalar stiffenerHeight, int stiffenerHeightNum,
      TacsScalar stiffenerThick, int stiffenerThickNum, int numStiffenerPlies,
      TacsScalar stiffenerPlyAngles[], TacsScalar stiffenerPlyFracs[],
      int stiffenerPlyFracNums[], TacsScalar panelWidth, int panelWidthNum,
      TacsScalar flangeFraction = 1.0, bool CPTstiffenerCrippling = false,
      TACSPanelGPs *panelGPs = nullptr);

  ~TACSGPBladeStiffenedShellConstitutive();

  // ==============================================================================
  // Tests
  // ==============================================================================

  /**
   *
   * @brief Test the nondimensional parameter computations against finite
   * difference. Tests all non-dimensional parameter subtests inside this one
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error among all non-dimensional parameter
   * subtests
   */
  TacsScalar testNondimensionalParameters(TacsScalar epsilon, int printLevel);

  /**
   *
   * @brief Test the axial critical loads => includes global axial and local
   * axial buckling subtests
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of global axial and local axial subtests
   */
  TacsScalar testAxialCriticalLoads(TacsScalar epsilon, int printLevel);

  /**
   *
   * @brief Test the shear critical loads => includes local shear and global
   * shear subtests
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testShearCriticalLoads(TacsScalar epsilon, int printLevel);

  /**
   *
   * @brief Test the crippling critical loads
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testStiffenerCripplingLoad(TacsScalar epsilon, int printLevel);

  /**
   *
   * @brief Test all GP tests
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error among all the subtest derivative tests
   */
  TacsScalar testAllTests(TacsScalar epsilon, int printLevel);

  // ==============================================================================
  // Getter and setters
  // ==============================================================================

  // get the three Gaussian Process model pointers
  TACSBucklingGaussianProcessModel *getAxialGP() {
    if (this->panelGPs) {
      return this->panelGPs->getAxialGP();
    } else {
      return nullptr;
    }
  }
  TACSBucklingGaussianProcessModel *getShearGP() {
    if (this->panelGPs) {
      return this->panelGPs->getShearGP();
    } else {
      return nullptr;
    }
  }
  TACSBucklingGaussianProcessModel *getCripplingGP() {
    if (this->panelGPs) {
      return this->panelGPs->getCripplingGP();
    } else {
      return nullptr;
    }
  }

  // Retrieve the global design variable numbers
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]) override;

  // Set the element design variable from the design vector
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]) override;

  // Get the element design variables values
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]) override;

  // Get the lower and upper bounds for the design variable values
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]) override;

  // Retrieve the design variable for plotting purposes
  TacsScalar evalDesignFieldValue(int elemIndex, const double pt[],
                                  const TacsScalar X[], int index) override;

  // set the KS weight for the failure constraints and the GP models (if GP
  // models are active)
  void setKSWeight(double ksWeight);

  // set the DV write out modes for f5 files [0 - regular DVs, 1 - ND, 2 - fail
  // indexes]
  void setWriteDVMode(int newMode) { writeDVmode = newMode; }

  // choose whether to use CPT (classical plate theory) analytical solution
  // (True) or DOD experimental buckling solution (False)
  void setCPTstiffenerCrippling(bool _CPTstiffenerCrippling) {
    CPTstiffenerCrippling = _CPTstiffenerCrippling;
  }

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
  // Override buckling load methods from super class
  // ==============================================================================

  /**
   * @brief Compute the strength ratio for the local buckling of the panel skin
   * between stiffeners
   *
   * @param e Shell strains
   * @return TacsScalar Strength ratio
   */
  TacsScalar evalLocalPanelBuckling(const TacsScalar e[]) override;

  /**
   * @brief Compute the strength ratio for the global buckling of the panel
   *
   * @param e Shell strains
   * @return TacsScalar Strength ratio
   */
  TacsScalar evalGlobalPanelBuckling(const TacsScalar e[]) override;

  /**
   * @brief Compute the strength ratio with respect to stiffener crippling
   *
   * Uses methods described in section 8.5 of "Design and Analysis of Composite
   * Structures with Application to Aerospace Structures, 2nd Edition" by
   * Christos Kassapoglou for Ali's solution.
   *
   * Optional ML model included from Sean's work as well.
   *
   * @param stiffenerStrain Stiffener centroid beam strains
   * @return TacsScalar Strength ratio
   */
  TacsScalar evalStiffenerCrippling(
      const TacsScalar stiffenerStrain[]) override;

  /**
   * @brief Compute the sensitivity of the local panel buckling strength ratio
   *
   * @param e Shell strains
   * @param sens Sensitivity of the output w.r.t the shell strains
   * @return TacsScalar Strength Ratio
   */
  TacsScalar evalLocalPanelBucklingStrainSens(const TacsScalar e[],
                                              TacsScalar sens[]) override;

  /**
   * @brief Compute the sensitivity of the global buckling strength ratio w.r.t
   * the shell strains
   *
   * @param e Shell strains
   * @param sens Sensitivity of the output w.r.t the shell strains
   * @return TacsScalar Strength Ratio
   */
  TacsScalar evalGlobalPanelBucklingStrainSens(const TacsScalar e[],
                                               TacsScalar sens[]) override;

  /**
   * @brief Compute the sensitivity of the stiffener crippling strength ratio
   * w.r.t the stiffener strains
   *
   * @param stiffenerStrain stiffener shell strains strains
   * @param sens Sensitivity of the output w.r.t the shell strains
   * @return TacsScalar Strength Ratio
   */
  TacsScalar evalStiffenerCripplingStrainSens(
      const TacsScalar stiffenerStrain[], TacsScalar sens[]) override;

  /**
   * @brief Add the derivative of the local panel buckling strength ratio w.r.t
   * the design variables

    @param elemIndex The local element index (not used)
    @param scale Value by which to scale the derivatives
    @param pt The parametric point (not used)
    @param X The physical node location (not used)
    @param strain The shell strains
    @param dvLen The length of the design vector (not used)
    @param dfdx The DV sensitivity array to add to
   */
  void addLocalPanelBucklingDVSens(int elemIndex, TacsScalar scale,
                                   const double pt[], const TacsScalar X[],
                                   const TacsScalar strain[], int dvLen,
                                   TacsScalar dfdx[]) override;

  /**
   * @brief Add the derivative of the global panel buckling strength ratio w.r.t
   * the design variables

    @param elemIndex The local element index (not used)
    @param scale Value by which to scale the derivatives
    @param pt The parametric point (not used)
    @param X The physical node location (not used)
    @param strain The shell strains
    @param dvLen The length of the design vector (not used)
    @param dfdx The DV sensitivity array to add to
   */
  void addGlobalPanelBucklingDVSens(int elemIndex, TacsScalar scale,
                                    const double pt[], const TacsScalar X[],
                                    const TacsScalar strain[], int dvLen,
                                    TacsScalar dfdx[]) override;

  /**
   * @brief Compute the sensitivity of the stiffener crippling strength ratio
   * w.r.t the design variables
   *
   * @param scale the derivative scalar for the stiffener crippling failure
   * index coming in
   * @param stiffenerStrain the shell strains in the stiffener
   * @param dfdx the output DV derivatives sensitivity
   */
  void addStiffenerCripplingDVSens(const TacsScalar scale,
                                   const TacsScalar stiffenerStrain[],
                                   TacsScalar dfdx[]) override;

  // ==============================================================================
  // Stiffener crippling helper functions
  // ==============================================================================

  /**
   * @brief compute the in plane load in the stiffener N11,stiff
   *
   * @param stiffenerStrain the midplane strains in the stiffener
   * @param A11s the A11 material constant for the stiffener
   * @return the stiffener in plane load N11,stiff
   */
  TacsScalar computeStiffenerInPlaneLoad(const TacsScalar stiffenerStrain[],
                                         TacsScalar *A11s);

  /**
   * @brief compute the DV sensitivities of the stiffener in plane load
   * N11,stiff
   *
   * @param scale the derivative df/dN11,stiff of the output in plane load
   * @param stiffenerStrain the midplane strains eps11 in the stiffener shell
   * @param dN11_stiff the jacobian dstiff_fail_index / dN11,stiff
   * @param dfdx the final updated DV derivatives df/dx
   */
  TacsScalar computeStiffenerInPlaneLoadSens(const TacsScalar scale,
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
      const TacsScalar a, const TacsScalar b, TacsScalar *D11sens,
      TacsScalar *D22sens, TacsScalar *asens, TacsScalar *bsens);

  /**
   *
   * @brief Test the affine aspect ratio sensitivity.
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testAffineAspectRatio(const TacsScalar epsilon, int printLevel);

  /**
   * @brief Compute the non-dimensional laminate isotropy xi = (D12 + 2 D66)
   * / sqrt(D11 * D22)
   *
   * @param D11 D11 stiffness
   * @param D22 D22 stiffness
   * @param D12 D12 stiffness
   * @param D66 D66 stiffness
   *
   */
  static inline TacsScalar computeLaminateIsotropy(const TacsScalar D11,
                                                   const TacsScalar D22,
                                                   const TacsScalar D12,
                                                   const TacsScalar D66) {
    return (D12 + 2.0 * D66) / sqrt(D11 * D22);
  }

  /**
   * @brief Compute the derivatives of the non-dimensional laminate isotropy
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
  TacsScalar computeLaminateIsotropySens(
      const TacsScalar xisens, const TacsScalar D11, const TacsScalar D22,
      const TacsScalar D12, const TacsScalar D66, TacsScalar *D11sens,
      TacsScalar *D22sens, TacsScalar *D12sens, TacsScalar *D66sens);

  /**
   *
   * @brief Test the generalized rigidity sensitivity.
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testLaminateIsotropy(const TacsScalar epsilon, int printLevel);

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
                                                 TacsScalar *D12sens,
                                                 TacsScalar *D66sens);

  /**
   *
   * @brief Test the generalized poisson's ratio sensitivity.
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testGeneralizedPoissonsRatio(const TacsScalar epsilon,
                                          int printLevel);

  /**
   * @brief Compute the non-dimensional stiffener area ratio delta = E1s * As /
   * (E1p * s_p * h) based on the DVs => stiffenerPitch, panelThick,
   * stiffenerHeight, stiffenerThick
   *
   * @param E1p effective E11p modulus
   * @param E1s effective E11s modulus
   */
  TacsScalar computeStiffenerAreaRatio(const TacsScalar E1p,
                                       const TacsScalar E1s);

  /**
   * @brief Compute the sensitivities of the non-dimensional stiffener area
   * ratio delta = E1s * As / (E1p * s_p * h)
   *
   * @param deltasens backpropagated sensitivity w.r.t the stiffener area ratio
   * delta
   * @param E1p effective E11p modulus
   * @param E1s effective E11s modulus
   * @param sthickSens stiffener thickness sens
   * @param sheightSens stiffener height sens
   * @param spitchSens stiffener pitch sens
   * @param pthickSens panel thickness sens
   * @param E1pSens effective E11p modulus sens
   * @param E1sSens effective E11s modulus sens
   *
   */
  TacsScalar computeStiffenerAreaRatioSens(
      const TacsScalar deltasens, const TacsScalar E1p, const TacsScalar E1s,
      TacsScalar *sthickSens, TacsScalar *sheightSens, TacsScalar *spitchSens,
      TacsScalar *pthickSens, TacsScalar *E1psens, TacsScalar *E1ssens);

  /**
   *
   * @brief Test the stiffener area ratio sensitivity aka delta parameter
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testStiffenerAreaRatio(const TacsScalar epsilon, int printLevel);

  /**
   * @brief Compute the non-dimensional stiffener-to-panel stiffness ratio gamma
   * = E1s * Is / (sp * D11) based on the DVs => stiffenerPitch, panelThick,
   * stiffenerHeight, stiffenerThick
   *
   * @param D11 the D11 stiffness of the plate
   * @param E1s the E1s effective stiffener modulus
   * @param zn the overall centroid
   */
  TacsScalar computeStiffenerStiffnessRatio(const TacsScalar D11,
                                            const TacsScalar E1s,
                                            const TacsScalar zn);

  /**
   * @brief Compute the sensitivities of the non-dimensional  stiffener-to-panel
   * stiffness ratio gamma = E1s * Is / (sp * D11)
   *
   * @param gammaSens backpropagated derivative through gamma output
   * @param D11 the D11 stiffness of the plate
   * @param E1s the E1s effective stiffener modulus
   * @param zn the overall centroid
   * @param D11Sens sensitivity w.r.t. the D11 stiffness
   * @param sthickSens stiffener thickness sens
   * @param sheightSens stiffener height sens
   * @param spitchSens stiffener pitch sens
   * @param E1sSens E1s effective stiffener modulus sens
   * @param znSens overall centroid sens
   *
   */
  TacsScalar computeStiffenerStiffnessRatioSens(
      const TacsScalar gammaSens, const TacsScalar D11, const TacsScalar E1s,
      const TacsScalar zn, TacsScalar *D11Sens, TacsScalar *sthickSens,
      TacsScalar *sheightSens, TacsScalar *spitchSens, TacsScalar *E1sSens,
      TacsScalar *znSens);

  /**
   *
   * @brief Test the stiffener stiffness ratio sensitivity aka gamma parameter
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
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
      const TacsScalar b, const TacsScalar h, TacsScalar *A66sens,
      TacsScalar *A11sens, TacsScalar *bsens, TacsScalar *hsens);

  /**
   *
   * @brief Test the transverse shear parameter sensitivity aka zeta parameter
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
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

  /*
   * computes the overall centroid zc of the panel and stiffener cross-section
   * in the 23-plane
   *
   * @param
   */
  TacsScalar computeOverallCentroid(const TacsScalar E1p, const TacsScalar E1s);

  void computeOverallCentroidSens(const TacsScalar znSens, const TacsScalar E1p,
                                  const TacsScalar E1s, TacsScalar *sthickSens,
                                  TacsScalar *sheightSens,
                                  TacsScalar *pthickSens,
                                  TacsScalar *spitchSens, TacsScalar *E1pSens,
                                  TacsScalar *E1sSens);

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
      TacsScalar *D11sens, TacsScalar *D22sens, TacsScalar *bsens,
      TacsScalar *deltasens, TacsScalar *rho_0sens, TacsScalar *xisens,
      TacsScalar *gammasens, TacsScalar *zetasens);

  /**
   *
   * @brief Test the critical global axial load function
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testCriticalGlobalAxialLoad(const TacsScalar epsilon,
                                         int printLevel);

  /**
   *
   * @brief Test the overall centroid method
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testOverallCentroid(const TacsScalar epsilon, int printLevel);

  /**
   *
   * @brief Test the computation of D11p with overall centroid
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testPanelGlobalBucklingStiffness(const TacsScalar epsilon,
                                              int printLevel);

  /**
   *
   * @brief Test the other util tests
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testOtherTests(const TacsScalar epsilon, int printLevel);

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
   * @param spitchsens Sensitivity w.r.t. the stiffener pitch
   * @param rho_0sens Sensitivity w.r.t.  affine aspect ratio
   * @param xisens Sensitivity w.r.t.  generalized rigidity
   * @param zetasens Sensitivity w.r.t. the transverse shear stiffness parameter
   *
   * @return TacsScalar The critical axial load for the global buckling mode
   */
  TacsScalar computeCriticalLocalAxialLoadSens(
      const TacsScalar N1sens, const TacsScalar D11, const TacsScalar D22,
      const TacsScalar rho_0, const TacsScalar xi, const TacsScalar zeta,
      TacsScalar *D11sens, TacsScalar *D22sens, TacsScalar *spitchsens,
      TacsScalar *rho_0sens, TacsScalar *xisens, TacsScalar *zetasens);

  /**
   *
   * @brief Test the critical local axial load function
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
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
   * @param rho_0 the affine aspect ratio
   * @param xi the generalized rigidity
   * @param gamma the stiffener-to-(local panel) bending stiffness ratio
   * @param zeta the transverse shear stiffness parameter
   *
   * returns:
   * @return N12Crit  the approximate critical buckling load
   */
  TacsScalar computeCriticalGlobalShearLoad(
      const TacsScalar D11, const TacsScalar D22, const TacsScalar b,
      const TacsScalar rho_0, const TacsScalar xi, const TacsScalar gamma,
      const TacsScalar zeta);

  /**
   * @brief Compute the sensitivity of the critical shear load for global
   * buckling
   *
   * @param N12sens the backpropagated sensitivity of the output N12crit
   * @param D11 the D11 stiffness of the plate
   * @param D22 the D22 stiffness of the plate
   * @param b the panel width
   * @param rho_0 the generalized affine aspect ratio
   * @param xi the generalized rigidity
   * @param gamma the stiffener-to-(local panel) bending stiffness ratio
   * @param zeta the transverse shear stiffness parameter
   * @param D11sens sensitivity w.r.t. the D11 stiffness of the plate
   * @param D22sens sensitivity w.r.t. the D22 stiffness of the plate
   * @param bsens sensitivity w.r.t. the panel width
   * @param rho_0sens sensitivity w.r.t. the affine aspect ratio
   * @param xisens sensitivity w.r.t. the generalized rigidity
   * @param gammasens sensitivity w.r.t. the stiffener-to-(local panel) bending
   * @param zetasens sensitivity w.r.t. the transverse shear stiffness parameter
   * stiffness ratio
   * @return TacsScalar The critical shear buckling load
   */
  TacsScalar computeCriticalGlobalShearLoadSens(
      const TacsScalar N12sens, const TacsScalar D11, const TacsScalar D22,
      const TacsScalar b, const TacsScalar rho_0, const TacsScalar xi,
      const TacsScalar gamma, const TacsScalar zeta, TacsScalar *D11sens,
      TacsScalar *D22sens, TacsScalar *bsens, TacsScalar *rho_0sens,
      TacsScalar *xisens, TacsScalar *gammasens, TacsScalar *zetasens);

  /**
   *
   * @brief Test the critical global shear load function
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testCriticalGlobalShearLoad(const TacsScalar epsilon,
                                         int printLevel);

  /**
   *
   * @brief Test the critical global shear load function at low rho0
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testCriticalGlobalShearLoad_LowAR(const TacsScalar epsilon,
                                               int printLevel);

  /**
   * @brief Compute the critical shear load for local buckling
   *
   *
   * input:
   * @param D11 the D11 stiffness of the plate
   * @param D22 the D22 stiffness of the plate
   * @param rho_0 the affine aspect ratio
   * @param xi the generalized rigidity
   * @param zeta the transverse shear stiffness parameter
   *
   * returns:
   * @return N12Crit  the approximate critical buckling load
   */
  TacsScalar computeCriticalLocalShearLoad(const TacsScalar D11,
                                           const TacsScalar D22,
                                           const TacsScalar rho_0,
                                           const TacsScalar xi,
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
   * @param rho0 the affine aspect ratio
   * @param xi the generalized rigidity
   * @param gamma the stiffener-to-(local panel) bending stiffness ratio
   * @param D11sens sensitivity w.r.t. the D11 stiffness of the plate
   * @param D22sens sensitivity w.r.t. the D22 stiffness of the plate
   * @param spitchsens sensitivity w.r.t. the stiffener pitch
   * @param rho_0sens sensitivity w.r.t. the affine aspect ratio
   * @param xisens sensitivity w.r.t. the generalized rigidity
   * @param zeta sensitivity w.r.t. the transverse shear stiffness parameter
   * @return TacsScalar The critical shear buckling load
   */
  TacsScalar computeCriticalLocalShearLoadSens(
      const TacsScalar N12sens, const TacsScalar D11, const TacsScalar D22,
      const TacsScalar rho_0, const TacsScalar xi, const TacsScalar zeta,
      TacsScalar *D11sens, TacsScalar *D22sens, TacsScalar *spitchsens,
      TacsScalar *rho_0sens, TacsScalar *xisens, TacsScalar *zetasens);

  /**
   *
   * @brief Test the critical local shear load function at high AR rho0
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testCriticalLocalShearLoad(const TacsScalar epsilon,
                                        int printLevel);

  /**
   *
   * @brief Test the critical local shear load function at low AR rho0
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testCriticalLocalShearLoad_LowAR(const TacsScalar epsilon,
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
                         TacsScalar *lam1bar, TacsScalar *lam2bar);

  /**
   * @brief Compute the panel global buckling constant D11 adjusted for the
   * overall centroid of the panel & stiffener
   *
   * @param E1p the effective modulus of the panel Q11 - Q12^2 / Q22
   * @param zn the overall centroid in the 23-plane or 1-plane of the panel &
   * stiffeners
   * @return the overall centroid-adjusted D11 constant
   */
  void computePanelGlobalBucklingStiffness(const TacsScalar E1p,
                                           const TacsScalar zn, TacsScalar *D1);

  /**
   * @brief compute sensitivities of the overall centroid-adjusted D11 material
   * constant
   *
   * @param D1Sens the derivative at the D11 level df/dD11, input for
   * backpropagation
   * @param E1p the effective modulus of the panel Q11 - Q12^2 / Q22
   * @param zn the overall centroid in the 23-plane or 1-plane of the panel &
   * stiffeners
   * @param spitchsens the output sensitivity of stiffener pitch df/dspitch
   * @param pthicksens the output sensitivity of panel thickness df/dpthick
   * @param E1pSens the output sensitivity of effective modulus of the panel
   * df/dE1p
   * @param znSens the output sensitivity of the overall centroid, df/dzn
   */
  void computePanelGlobalBucklingStiffnessSens(
      const TacsScalar D1Sens, const TacsScalar E1p, const TacsScalar zn,
      TacsScalar *spitchSens, TacsScalar *pthickSens, TacsScalar *E1pSens,
      TacsScalar *znSens);

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
                             TacsScalar *lam1bar, TacsScalar *lam2bar,
                             TacsScalar *dl1xi, TacsScalar *dl1gamma,
                             TacsScalar *dl2xi, TacsScalar *dl2gamma);

  /**
   *
   * @brief Test the nondimShearParams function used in the shear buckling loads
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
   */
  TacsScalar testNondimShearParams(const TacsScalar epsilon, int printLevel);

  /**
   * @brief the constraint for shear buckling closed-form solution R(lam1, lam2)
   *  which is reformulated as constraint on just lam2_bar R(lam2_bar)
   *
   * @param lam2sq input lam2_bar^2 non-dimensional mode shape parameter for the
   * shear mode
   * @param xi the non-dimensional laminate isotropy
   * @param gamma the non-dimensional stiffener stiffness ratio
   * @return the residual R(lam2_bar^2)
   */
  TacsScalar lam2Constraint(const TacsScalar lam2sq, const TacsScalar xi,
                            const TacsScalar gamma);

  /**
   * @brief the derivative of the shear buckling closed-form solution constraint
   * reformulated as R(lam2_bar) = 0
   *
   * @param lam2sq input lam2_bar^2 non-dimensional mode shape parameter for the
   * shear mode
   * @param xi the non-dimensional laminate isotropy
   * @param gamma the non-dimensional stiffener stiffness ratio
   * @return returns the Jacobian of the constraint, dR/dlam2_bar
   */
  TacsScalar lam2ConstraintDeriv(const TacsScalar lam2sq, const TacsScalar xi,
                                 const TacsScalar gamma);

  /**
   *
   * @brief Test the lam2Constraint function used in the shear buckling
   * closed-form solution
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error of the derivative test
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
      const TacsScalar zeta, TacsScalar *D11sens, TacsScalar *D22sens,
      TacsScalar *sheightsens, TacsScalar *xisens, TacsScalar *rho_0sens,
      TacsScalar *genPoiss_sens, TacsScalar *zetasens);

  // ==============================================================================
  // Attributes
  // ==============================================================================

  // Machine learning Gaussian Process models stored in panel GP class
  TACSPanelGPs *panelGPs;

  // --- Design variable values ---
  TacsScalar panelWidth;  ///< Panel width
  static const int NUM_CF_MODES =
      50;  // number of modes used in closed-form method

  TacsScalar panelWidthLowerBound = 0.0;
  TacsScalar panelWidthUpperBound = 1e20;

  // --- Design variable numbers ---
  int panelWidthNum;       ///< Panel width DV number
  int panelWidthLocalNum;  ///< Panel width local DV number

  // stiffener crippling prediction
  bool CPTstiffenerCrippling = false;

  // debugging modes
  int writeDVmode = 0;  // 0 - normal DVs, 1 - NDparams

  // pointers for Xtest vectors used in GP computation
  TacsScalar *XtestAxial, *XtestShear, *XtestCrippling;
  TacsScalar *XtestAxialSens, *XtestShearSens, *XtestCripplingSens;

 private:
  // private so that subclass constName for GP buckling constraints doesn't
  // conflict with superclass
  static const char *constName;  ///< Constitutive model name
};
