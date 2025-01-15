/*
=====================================================================================================
Blade-Stiffened Shell Constitutive Model using Gaussian Process Machine Learning
Buckling Constraints
=====================================================================================================
@File    :   TACSGPBladeStiffenedShellConstutive.cpp
@Date    :   2024/04/24
@Author  :   Sean Phillip Engelstad
@Description : Constitutive model for a blade-stiffened shell. Based on the FSDT
blade models adopted by Alasdair Christison Gray in the
TACSBladeStiffenedShellConstutitutive.h class and the original one by Graeme
Kennedy. Gaussian Processes for Machine Learning are used for the buckling
constraints of the stiffened panels.
*/

// =============================================================================
// Standard Library Includes
// =============================================================================

// =============================================================================
// Extension Includes
// =============================================================================
#include "TACSGPBladeStiffenedShellConstitutive.h"

#include "TACSMaterialProperties.h"
#include "TACSShellConstitutive.h"

const char* TACSGPBladeStiffenedShellConstitutive::constName =
    "TACSGPBladeStiffenedShellConstitutive";

// ==============================================================================
// Constructor
// ==============================================================================

TACSGPBladeStiffenedShellConstitutive::TACSGPBladeStiffenedShellConstitutive(
    TACSOrthotropicPly* panelPly, TACSOrthotropicPly* stiffenerPly,
    TacsScalar kcorr, TacsScalar panelLength, int panelLengthNum,
    TacsScalar stiffenerPitch, int stiffenerPitchNum, TacsScalar panelThick,
    int panelThickNum, int numPanelPlies, TacsScalar panelPlyAngles[],
    TacsScalar panelPlyFracs[], int panelPlyFracNums[],
    TacsScalar stiffenerHeight, int stiffenerHeightNum,
    TacsScalar stiffenerThick, int stiffenerThickNum, int numStiffenerPlies,
    TacsScalar stiffenerPlyAngles[], TacsScalar stiffenerPlyFracs[],
    int stiffenerPlyFracNums[], TacsScalar panelWidth, int panelWidthNum,
    TacsScalar flangeFraction, bool CPTstiffenerCrippling,
    TACSPanelGPs* panelGPs)
    : TACSBladeStiffenedShellConstitutive(
          panelPly, stiffenerPly, kcorr, panelLength, panelLengthNum,
          stiffenerPitch, stiffenerPitchNum, panelThick, panelThickNum,
          numPanelPlies, panelPlyAngles, panelPlyFracs, panelPlyFracNums,
          stiffenerHeight, stiffenerHeightNum, stiffenerThick,
          stiffenerThickNum, numStiffenerPlies, stiffenerPlyAngles,
          stiffenerPlyFracs, stiffenerPlyFracNums, flangeFraction) {
  // DVs section, only one new DV - panelWidth
  // --- Panel width values ---
  this->panelWidth = panelWidth;
  this->panelWidthNum = panelWidthNum;
  this->panelWidthLocalNum = -1;
  if (panelWidthNum >= 0) {
    this->panelWidthLocalNum = this->numDesignVars;
    this->numDesignVars++;
    this->numGeneralDV++;
  }
  this->panelWidthLowerBound = 0.000;
  this->panelWidthUpperBound = 1e20;

  // set whether to use CPT stiffener crippling from Sean's paper vs.
  //    experimental DOD manuscript crippling solution from Ali's superclass
  //    (default false since DOD probably more accurate)
  this->CPTstiffenerCrippling = CPTstiffenerCrippling;

  // set Gaussian process models in
  this->panelGPs = panelGPs;
  if (this->panelGPs) {
    this->panelGPs->incref();
  }

  // allocate Xtest temporary vectors for each GP
  if (this->getAxialGP()) {
    XtestAxial = new TacsScalar[this->getAxialGP()->getNparam()];
    XtestAxialSens = new TacsScalar[this->getAxialGP()->getNparam()];
  }
  if (this->getShearGP()) {
    XtestShear = new TacsScalar[this->getShearGP()->getNparam()];
    XtestShearSens = new TacsScalar[this->getShearGP()->getNparam()];
  }
  if (this->getCripplingGP()) {
    XtestCrippling = new TacsScalar[this->getCripplingGP()->getNparam()];
    XtestCripplingSens = new TacsScalar[this->getCripplingGP()->getNparam()];
  }

  // default value of ksWeight
  this->setKSWeight(100.0);
}

// ==============================================================================
// Destructor
// ==============================================================================
TACSGPBladeStiffenedShellConstitutive::
    ~TACSGPBladeStiffenedShellConstitutive() {
  // Don't call the base class destructor C++ does that automatically (and will
  // seg fault if you call it here.)

  if (this->panelGPs) {
    this->panelGPs->decref();
    this->panelGPs = nullptr;
  }

  if (this->getAxialGP()) {
    delete[] XtestAxial;
    XtestAxial = nullptr;

    delete[] XtestAxialSens;
    XtestAxialSens = nullptr;
  }

  if (this->getShearGP()) {
    delete[] XtestShear;
    XtestShear = nullptr;

    delete[] XtestShearSens;
    XtestShearSens = nullptr;
  }

  if (this->getCripplingGP()) {
    delete[] XtestCrippling;
    XtestCrippling = nullptr;

    delete[] XtestCripplingSens;
    XtestCripplingSens = nullptr;
  }
}

// ==============================================================================
// Override Failure constraint and sensitivities
// ==============================================================================

TacsScalar TACSGPBladeStiffenedShellConstitutive::evalLocalPanelBuckling(
    const TacsScalar e[]) {
  // this routine computes N11,cr for the local panel section with size a x s_p
  // (in between stiffeners)

  // compute panel stiffness matrix
  TacsScalar panelStiffness[NUM_TANGENT_STIFFNESS_ENTRIES],
      panelStress[NUM_STRESSES];
  this->computePanelStiffness(panelStiffness);
  const TacsScalar *Ap, *Dp;
  this->extractTangentStiffness(panelStiffness, &Ap, NULL, &Dp, NULL, NULL);
  this->computePanelStress(e, panelStress);

  // extract panel stiffnesses and dimensions
  TacsScalar D11Local = Dp[0];
  TacsScalar D12p = Dp[1], D66p = Dp[5], D22p = Dp[3];
  TacsScalar A11p = Ap[0], A66p = Ap[5];
  TacsScalar a = this->panelLength;
  TacsScalar b = this->panelWidth;
  TacsScalar s_p = this->stiffenerPitch;

  // compute non-dimensional parameters for the local panel
  TacsScalar rho0Local = computeAffineAspectRatio(
      D11Local, D22p, a, s_p);  // local panel is a x s_p
  TacsScalar xiLocal = computeLaminateIsotropy(D11Local, D22p, D12p, D66p);
  TacsScalar zetaPanel =
      computeTransverseShearParameter(A66p, A11p, b, this->panelThick);

  // compute the pure axial and pure shear buckling loads
  TacsScalar N1CritLocal = computeCriticalLocalAxialLoad(
      D11Local, D22p, rho0Local, xiLocal, zetaPanel);
  TacsScalar N12CritLocal = computeCriticalLocalShearLoad(
      D11Local, D22p, rho0Local, xiLocal, zetaPanel);

  // compute the combined loading buckling failure index
  return this->bucklingEnvelope(-panelStress[0], N1CritLocal, panelStress[2],
                                N12CritLocal);
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::evalGlobalPanelBuckling(
    const TacsScalar e[]) {
  // this routine computes N11,cr for the global panel with the stiffeners
  // applied

  // compute panel stiffness matrix
  TacsScalar panelStiffness[NUM_TANGENT_STIFFNESS_ENTRIES],
      panelStress[NUM_STRESSES];
  this->computePanelStiffness(panelStiffness);
  const TacsScalar *Ap, *Dp;
  this->extractTangentStiffness(panelStiffness, &Ap, NULL, &Dp, NULL, NULL);
  this->computePanelStress(e, panelStress);

  // compute effective moduli, overall centroid
  TacsScalar E1s, E1p, _;
  this->computeEffectiveModulii(this->numPanelPlies, this->panelQMats,
                                this->panelPlyFracs, &E1p, &_);
  this->computeEffectiveModulii(this->numStiffenerPlies, this->stiffenerQMats,
                                this->stiffenerPlyFracs, &E1s, &_);
  TacsScalar zn = this->computeOverallCentroid(E1p, E1s);

  // get the global buckling D11 with overall stiffener + panel centroid
  TacsScalar D11Global;
  computePanelGlobalBucklingStiffness(E1p, zn, &D11Global);

  // extract panel stiffnesses and dimensions
  TacsScalar D12p = Dp[1], D66p = Dp[5], D22p = Dp[3];
  TacsScalar A11p = Ap[0], A66p = Ap[5];
  TacsScalar a = this->panelLength;
  TacsScalar b = this->panelWidth;
  TacsScalar s_p = this->stiffenerPitch;

  // compute non-dimensional parameters for the global panel
  TacsScalar delta = computeStiffenerAreaRatio(E1p, E1s);
  TacsScalar rho0Global =
      computeAffineAspectRatio(D11Global, D22p, a, b);  // global panel is a x b
  TacsScalar xiGlobal = computeLaminateIsotropy(D11Global, D22p, D12p, D66p);
  TacsScalar gamma = computeStiffenerStiffnessRatio(D11Global, E1s, zn);
  TacsScalar zetaPanel =
      computeTransverseShearParameter(A66p, A11p, b, this->panelThick);

  // compute the pure axial and pure shear buckling loads
  TacsScalar N1CritGlobal = computeCriticalGlobalAxialLoad(
      D11Global, D22p, b, delta, rho0Global, xiGlobal, gamma, zetaPanel);
  TacsScalar N12CritGlobal = computeCriticalGlobalShearLoad(
      D11Global, D22p, b, rho0Global, xiGlobal, gamma, zetaPanel);

  // compute the combined loading buckling failure index
  return this->bucklingEnvelope(-panelStress[0], N1CritGlobal, panelStress[2],
                                N12CritGlobal);
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::evalStiffenerCrippling(
    const TacsScalar stiffenerStrain[]) {
  if (CPTstiffenerCrippling) {  // use predictions for Sean's paper
    // compute D matrix of the stiffener (treating it like a panel for
    // crippling)
    TacsScalar stiffenerCripplingStiffness[NUM_TANGENT_STIFFNESS_ENTRIES];
    const TacsScalar *As_crippling, *Ds_crippling;
    this->computeStiffenerCripplingStiffness(stiffenerCripplingStiffness);
    this->extractTangentStiffness(stiffenerCripplingStiffness, &As_crippling,
                                  NULL, &Ds_crippling, NULL, NULL);

    // get stiffener material props and dimensions
    TacsScalar A11s = As_crippling[0], A66s = As_crippling[5];
    TacsScalar D11s = Ds_crippling[0], D12s = Ds_crippling[1];
    TacsScalar D22s = Ds_crippling[3], D66s = Ds_crippling[5];
    TacsScalar bStiff = this->stiffenerHeight;
    TacsScalar hStiff = this->stiffenerThick;
    TacsScalar a = this->panelLength;

    // compute stiffener non-dimensional parameters
    TacsScalar rho0Stiff = computeAffineAspectRatio(D11s, D22s, a, bStiff);
    TacsScalar xiStiff = computeLaminateIsotropy(D11s, D22s, D12s, D66s);
    TacsScalar genPoiss = computeGeneralizedPoissonsRatio(D12s, D66s);
    TacsScalar zetaStiff =
        computeTransverseShearParameter(A66s, A11s, bStiff, hStiff);

    // Compute stiffener in plane load and crippling failure index
    TacsScalar A11s_beam;
    TacsScalar N1stiff =
        computeStiffenerInPlaneLoad(stiffenerStrain, &A11s_beam);
    TacsScalar N1CritCrippling = computeStiffenerCripplingLoad(
        D11s, D22s, xiStiff, rho0Stiff, genPoiss, zetaStiff);

    // return the stiffener failure index
    return N1stiff / N1CritCrippling;

  } else {  // CPTstiffenerCrippling = false
    // the DOD experimental stiffener crippling solution from Ali's superclass
    return TACSBladeStiffenedShellConstitutive::evalStiffenerCrippling(
        stiffenerStrain);
  }
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::evalLocalPanelBucklingStrainSens(
    const TacsScalar e[], TacsScalar sens[]) {
  // this routine computes N11,cr for the local panel section with size a x s_p
  // (in between stiffeners)

  // compute panel stiffness matrix
  TacsScalar panelStiffness[NUM_TANGENT_STIFFNESS_ENTRIES],
      panelStress[NUM_STRESSES];
  this->computePanelStiffness(panelStiffness);
  const TacsScalar *Ap, *Dp;
  this->extractTangentStiffness(panelStiffness, &Ap, NULL, &Dp, NULL, NULL);
  this->computePanelStress(e, panelStress);

  // extract panel stiffnesses and dimensions
  TacsScalar D11Local = Dp[0];
  TacsScalar D12p = Dp[1], D66p = Dp[5], D22p = Dp[3];
  TacsScalar A11p = Ap[0], A66p = Ap[5];
  TacsScalar a = this->panelLength;
  TacsScalar b = this->panelWidth;
  TacsScalar s_p = this->stiffenerPitch;

  // compute non-dimensional parameters for the local panel
  TacsScalar rho0Local = computeAffineAspectRatio(
      D11Local, D22p, a, s_p);  // local panel is a x s_p
  TacsScalar xiLocal = computeLaminateIsotropy(D11Local, D22p, D12p, D66p);
  TacsScalar zetaPanel =
      computeTransverseShearParameter(A66p, A11p, b, this->panelThick);

  // compute the pure axial and pure shear buckling loads
  TacsScalar N1CritLocal = computeCriticalLocalAxialLoad(
      D11Local, D22p, rho0Local, xiLocal, zetaPanel);
  TacsScalar N12CritLocal = computeCriticalLocalShearLoad(
      D11Local, D22p, rho0Local, xiLocal, zetaPanel);

  // Compute the buckling criteria and it's sensitivities
  TacsScalar N1LocalSens, N12LocalSens, N1CritLocalSens, N12CritLocalSens;
  const TacsScalar strengthRatio = this->bucklingEnvelopeSens(
      -panelStress[0], N1CritLocal, panelStress[2], N12CritLocal, &N1LocalSens,
      &N1CritLocalSens, &N12LocalSens, &N12CritLocalSens);

  sens[0] = N1LocalSens * -Ap[0] + N12LocalSens * Ap[2];
  sens[1] = N1LocalSens * -Ap[1] + N12LocalSens * Ap[4];
  sens[2] = N1LocalSens * -Ap[2] + N12LocalSens * Ap[5];
  for (int ii = 3; ii < NUM_STRESSES; ii++) {
    sens[ii] = 0.0;
  }

  return strengthRatio;
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::evalGlobalPanelBucklingStrainSens(
    const TacsScalar e[], TacsScalar sens[]) {
  // this routine computes N11,cr for the global panel with the stiffeners
  // applied

  // compute panel stiffness matrix
  TacsScalar panelStiffness[NUM_TANGENT_STIFFNESS_ENTRIES],
      panelStress[NUM_STRESSES];
  this->computePanelStiffness(panelStiffness);
  const TacsScalar *Ap, *Dp;
  this->extractTangentStiffness(panelStiffness, &Ap, NULL, &Dp, NULL, NULL);
  this->computePanelStress(e, panelStress);

  // compute effective moduli, overall centroid
  TacsScalar E1s, E1p, _;
  this->computeEffectiveModulii(this->numPanelPlies, this->panelQMats,
                                this->panelPlyFracs, &E1p, &_);
  this->computeEffectiveModulii(this->numStiffenerPlies, this->stiffenerQMats,
                                this->stiffenerPlyFracs, &E1s, &_);
  TacsScalar zn = this->computeOverallCentroid(E1p, E1s);

  // get the global buckling D11 with overall stiffener + panel centroid
  TacsScalar D11Global;
  computePanelGlobalBucklingStiffness(E1p, zn, &D11Global);

  // extract panel stiffnesses and dimensions
  TacsScalar D12p = Dp[1], D66p = Dp[5], D22p = Dp[3];
  TacsScalar A11p = Ap[0], A66p = Ap[5];
  TacsScalar a = this->panelLength;
  TacsScalar b = this->panelWidth;
  TacsScalar s_p = this->stiffenerPitch;

  // compute non-dimensional parameters for the global panel
  TacsScalar delta = computeStiffenerAreaRatio(E1p, E1s);
  TacsScalar rho0Global =
      computeAffineAspectRatio(D11Global, D22p, a, b);  // global panel is a x b
  TacsScalar xiGlobal = computeLaminateIsotropy(D11Global, D22p, D12p, D66p);
  TacsScalar gamma = computeStiffenerStiffnessRatio(D11Global, E1s, zn);
  TacsScalar zetaPanel =
      computeTransverseShearParameter(A66p, A11p, b, this->panelThick);

  // compute the pure axial and pure shear buckling loads
  TacsScalar N1CritGlobal = computeCriticalGlobalAxialLoad(
      D11Global, D22p, b, delta, rho0Global, xiGlobal, gamma, zetaPanel);
  TacsScalar N12CritGlobal = computeCriticalGlobalShearLoad(
      D11Global, D22p, b, rho0Global, xiGlobal, gamma, zetaPanel);

  // backprop sensitivities from combined loading to the in-plane loads for
  // strain sens
  TacsScalar N1GlobalSens, N1CritGlobalSens, N12GlobalSens, N12CritGlobalSens;
  const TacsScalar strengthRatio = this->bucklingEnvelopeSens(
      -panelStress[0], N1CritGlobal, panelStress[2], N12CritGlobal,
      &N1GlobalSens, &N1CritGlobalSens, &N12GlobalSens, &N12CritGlobalSens);

  sens[0] += N1GlobalSens * -Ap[0] + N12GlobalSens * Ap[2];
  sens[1] += N1GlobalSens * -Ap[1] + N12GlobalSens * Ap[4];
  sens[2] += N1GlobalSens * -Ap[2] + N12GlobalSens * Ap[5];
  for (int ii = 3; ii < NUM_STRESSES; ii++) {
    sens[ii] = 0.0;
  }

  // compute the combined loading buckling failure index
  return strengthRatio;
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::evalStiffenerCripplingStrainSens(
    const TacsScalar stiffenerStrain[], TacsScalar sens[]) {
  if (CPTstiffenerCrippling) {  // use predictions for Sean's paper
    // compute D matrix of the stiffener (treating it like a panel for
    // crippling)
    TacsScalar stiffenerCripplingStiffness[NUM_TANGENT_STIFFNESS_ENTRIES];
    const TacsScalar *As_crippling, *Ds_crippling;
    this->computeStiffenerCripplingStiffness(stiffenerCripplingStiffness);
    this->extractTangentStiffness(stiffenerCripplingStiffness, &As_crippling,
                                  NULL, &Ds_crippling, NULL, NULL);

    // get stiffener material props and dimensions
    TacsScalar A11s = As_crippling[0], A66s = As_crippling[5];
    TacsScalar D11s = Ds_crippling[0], D12s = Ds_crippling[1];
    TacsScalar D22s = Ds_crippling[3], D66s = Ds_crippling[5];
    TacsScalar bStiff = this->stiffenerHeight;
    TacsScalar hStiff = this->stiffenerThick;
    TacsScalar a = this->panelLength;

    // compute stiffener non-dimensional parameters
    TacsScalar rho0Stiff = computeAffineAspectRatio(D11s, D22s, a, bStiff);
    TacsScalar xiStiff = computeLaminateIsotropy(D11s, D22s, D12s, D66s);
    TacsScalar genPoiss = computeGeneralizedPoissonsRatio(D12s, D66s);
    TacsScalar zetaStiff =
        computeTransverseShearParameter(A66s, A11s, bStiff, hStiff);

    // Compute stiffener in plane load and crippling failure index
    TacsScalar A11s_beam;
    TacsScalar N1stiff =
        computeStiffenerInPlaneLoad(stiffenerStrain, &A11s_beam);
    TacsScalar N1CritCrippling = computeStiffenerCripplingLoad(
        D11s, D22s, xiStiff, rho0Stiff, genPoiss, zetaStiff);

    // applky stiffener strain sens through the in-plane load
    // Sensitivity from stiffener crippling (ignores dKSdf[4] since that is
    // applied in main evalFailureStrainSens here)
    TacsScalar N1stiffSens = 1.0 / N1CritCrippling;
    sens[0] += N1stiffSens * -A11s_beam;
    TacsScalar z = -1.0 * this->computeStiffenerCentroidHeight() -
                   0.5 * this->panelThick;  // stiffener downwards into wing
    sens[3] += N1stiffSens * z * -A11s_beam;

    // return the stiffener failure index
    return N1stiff / N1CritCrippling;

  } else {  // CPTstiffenerCrippling = false
    // the DOD experimental stiffener crippling solution from Ali's superclass
    return TACSBladeStiffenedShellConstitutive::
        evalStiffenerCripplingStrainSens(stiffenerStrain, sens);
  }
}

void TACSGPBladeStiffenedShellConstitutive::addLocalPanelBucklingDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], int dvLen, TacsScalar dfdx[]) {
  // scale is scale * dKSdf[2] here so already backpropped to the combined
  // axial-shear buckling load now we need to backprop through the failure
  // computation back to the DVs

  // compute panel stiffness matrix
  TacsScalar panelStiffness[NUM_TANGENT_STIFFNESS_ENTRIES],
      panelStress[NUM_STRESSES];
  this->computePanelStiffness(panelStiffness);
  const TacsScalar *Ap, *Dp;
  this->extractTangentStiffness(panelStiffness, &Ap, NULL, &Dp, NULL, NULL);
  this->computePanelStress(strain, panelStress);

  // extract panel stiffnesses and dimensions
  TacsScalar D11Local = Dp[0];
  TacsScalar D12p = Dp[1], D66p = Dp[5], D22p = Dp[3];
  TacsScalar A11p = Ap[0], A66p = Ap[5];
  TacsScalar a = this->panelLength;
  TacsScalar b = this->panelWidth;
  TacsScalar s_p = this->stiffenerPitch;

  // compute non-dimensional parameters for the local panel
  TacsScalar rho0Local = computeAffineAspectRatio(
      D11Local, D22p, a, s_p);  // local panel is a x s_p
  TacsScalar xiLocal = computeLaminateIsotropy(D11Local, D22p, D12p, D66p);
  TacsScalar zetaPanel =
      computeTransverseShearParameter(A66p, A11p, b, this->panelThick);

  // compute the pure axial and pure shear buckling loads
  TacsScalar N1CritLocal = computeCriticalLocalAxialLoad(
      D11Local, D22p, rho0Local, xiLocal, zetaPanel);
  TacsScalar N12CritLocal = computeCriticalLocalShearLoad(
      D11Local, D22p, rho0Local, xiLocal, zetaPanel);

  // backprop from combined failure index to the in-plane load sens, critical
  // buckling load sens
  TacsScalar N1LocalSens, N12LocalSens, N1CritLocalSens, N12CritLocalSens;
  const TacsScalar strengthRatio = this->bucklingEnvelopeSens(
      -panelStress[0], N1CritLocal, panelStress[2], N12CritLocal, &N1LocalSens,
      &N1CritLocalSens, &N12LocalSens, &N12CritLocalSens);

  // backprop in-plane load sens (stress sens) to the DVs
  TacsScalar dfdPanelStress[NUM_STRESSES];
  memset(dfdPanelStress, 0, NUM_STRESSES * sizeof(TacsScalar));
  dfdPanelStress[0] = -N1LocalSens;
  dfdPanelStress[2] = N12LocalSens;
  // figure out whether this should be based on
  this->addPanelStressDVSens(scale, strain, dfdPanelStress,
                             &dfdx[this->panelDVStartNum]);

  // backprop critical buckling load sens to the material properties and the
  // non-dimensional parameters
  // --------------------------

  // holder for local DV jacobians (we then multiply them by scale when we add
  // them into the dfdx) DVsens format [0 - panel length, 1 - stiff pitch, 2 -
  // panel thick,
  //                3 - stiff height, 4 - stiff thick, 5 - panel width]
  TacsScalar DVsens[6];
  memset(DVsens, 0, 6 * sizeof(TacsScalar));

  // define the material sensivitiies
  TacsScalar Dpsens[4];  // D11Local,D12,D22,D66
  TacsScalar Apsens[4];  // A11,A12,A22,A66
  memset(Dpsens, 0, 4 * sizeof(TacsScalar));
  memset(Apsens, 0, 4 * sizeof(TacsScalar));

  // define the non-dimensional (ND) sensitivities [rho0Local, xiLocal, zeta]
  TacsScalar NDsens[3];
  memset(NDsens, 0, 3 * sizeof(TacsScalar));

  // backprop from the critical axial + shear local buckling loads to the
  // material, ND params, DVs
  computeCriticalLocalAxialLoadSens(
      N1CritLocalSens, D11Local, D22p, rho0Local, xiLocal, zetaPanel,
      &Dpsens[0], &Dpsens[2], &DVsens[1], &NDsens[0], &NDsens[1], &NDsens[2]);
  computeCriticalLocalShearLoadSens(
      N12CritLocalSens, D11Local, D22p, rho0Local, xiLocal, zetaPanel,
      &Dpsens[0], &Dpsens[2], &DVsens[1], &NDsens[0], &NDsens[1], &NDsens[2]);

  // backprop from ND params to the material and DVs
  // -------------------

  // note we use s_p and spitchSens in place of b and bsens for local panel
  computeAffineAspectRatioSens(NDsens[0], D11Local, D22p, a, s_p, &Dpsens[0],
                               &Dpsens[2], &DVsens[0],
                               &DVsens[1]);  // 0 - backprop from rho0Local

  computeLaminateIsotropySens(NDsens[1], D11Local, D22p, D12p, D66p, &Dpsens[0],
                              &Dpsens[2], &Dpsens[1],
                              &Dpsens[3]);  // 1 - backprop from xiLocal

  // backprop from material sensitivities to DVs
  // --------------------

  // panel thickness term for material part..
  TacsScalar t = this->panelThick;
  if (this->panelThickNum >= 0) {
    int dvNum = this->panelThickLocalNum;
    // backpropagate through the D matrix
    TacsScalar dDfactor_dthick = 0.25 * t * t;  // d/dt(t^3/12) = t^2/4
    for (int ii = 0; ii < this->numPanelPlies; ii++) {
      TacsScalar* Q = &(this->panelQMats[ii * NUM_Q_ENTRIES]);
      dfdx[dvNum] += scale * dDfactor_dthick * this->panelPlyFracs[ii] *
                     (Dpsens[0] * Q[0] + Dpsens[2] * Q[3] + Dpsens[1] * Q[1] +
                      Dpsens[3] * Q[5]);
      dfdx[dvNum] += scale * this->panelPlyFracs[ii] *
                     (Apsens[0] * Q[0] + Apsens[2] * Q[3] + Apsens[1] * Q[1] +
                      Apsens[3] * Q[5]);
    }
  }

  // --- Panel Ply fraction sensitivities for A,D matrices ---
  for (int plyNum = 0; plyNum < this->numPanelPlies; plyNum++) {
    int dvNum = this->panelPlyFracLocalNums[plyNum];
    if (dvNum >= 0) {
      const TacsScalar* Q = &(this->panelQMats[plyNum * NUM_Q_ENTRIES]);
      dfdx[dvNum] += scale * (t * t * t / 12.0) *
                     (Dpsens[0] * Q[0] + Dpsens[2] * Q[3] + Dpsens[1] * Q[1] +
                      Dpsens[3] * Q[5]);
      dfdx[dvNum] += scale * t *
                     (Apsens[0] * Q[0] + Apsens[2] * Q[3] + Apsens[1] * Q[1] +
                      Apsens[3] * Q[5]);
    }
  }

  // combine local DV jacobians with scale to write into dfdx full derivatives
  // --------------------------------

  if (this->panelLengthLocalNum >= 0) {
    dfdx[this->panelLengthLocalNum] += scale * DVsens[0];
  }
  if (this->stiffenerPitchLocalNum >= 0) {
    dfdx[this->stiffenerPitchLocalNum] += scale * DVsens[1];
  }
  if (this->panelThickLocalNum >= 0) {
    dfdx[this->panelThickLocalNum] += scale * DVsens[2];
  }
  if (this->stiffenerHeightLocalNum >= 0) {
    dfdx[this->stiffenerHeightLocalNum] += scale * DVsens[3];
  }
  if (this->stiffenerThickLocalNum >= 0) {
    dfdx[this->stiffenerThickLocalNum] += scale * DVsens[4];
  }
  if (this->panelWidthLocalNum >= 0) {
    dfdx[this->panelWidthLocalNum] += scale * DVsens[5];
  }
}

void TACSGPBladeStiffenedShellConstitutive::addGlobalPanelBucklingDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], int dvLen, TacsScalar dfdx[]) {
  // scale is scale * dKSdf[3] here so already backpropped to the combined
  // axial-shear buckling load now we need to backprop through the failure
  // computation back to the DVs

  // compute panel stiffness matrix
  TacsScalar panelStiffness[NUM_TANGENT_STIFFNESS_ENTRIES],
      panelStress[NUM_STRESSES];
  this->computePanelStiffness(panelStiffness);
  const TacsScalar *Ap, *Dp;
  this->extractTangentStiffness(panelStiffness, &Ap, NULL, &Dp, NULL, NULL);
  this->computePanelStress(strain, panelStress);

  // compute effective moduli, overall centroid
  TacsScalar E1s, E1p, _;
  this->computeEffectiveModulii(this->numPanelPlies, this->panelQMats,
                                this->panelPlyFracs, &E1p, &_);
  this->computeEffectiveModulii(this->numStiffenerPlies, this->stiffenerQMats,
                                this->stiffenerPlyFracs, &E1s, &_);
  TacsScalar zn = this->computeOverallCentroid(E1p, E1s);

  // get the global buckling D11 with overall stiffener + panel centroid
  TacsScalar D11Global;
  computePanelGlobalBucklingStiffness(E1p, zn, &D11Global);

  // extract panel stiffnesses and dimensions
  TacsScalar D12p = Dp[1], D66p = Dp[5], D22p = Dp[3];
  TacsScalar A11p = Ap[0], A66p = Ap[5];
  TacsScalar a = this->panelLength;
  TacsScalar b = this->panelWidth;
  TacsScalar s_p = this->stiffenerPitch;
  TacsScalar h = this->panelThick;

  // compute non-dimensional parameters for the global panel
  TacsScalar delta = computeStiffenerAreaRatio(E1p, E1s);
  TacsScalar rho0Global =
      computeAffineAspectRatio(D11Global, D22p, a, b);  // global panel is a x b
  TacsScalar xiGlobal = computeLaminateIsotropy(D11Global, D22p, D12p, D66p);
  TacsScalar gamma = computeStiffenerStiffnessRatio(D11Global, E1s, zn);
  TacsScalar zetaPanel =
      computeTransverseShearParameter(A66p, A11p, b, this->panelThick);

  // compute the pure axial and pure shear buckling loads
  TacsScalar N1CritGlobal = computeCriticalGlobalAxialLoad(
      D11Global, D22p, b, delta, rho0Global, xiGlobal, gamma, zetaPanel);
  TacsScalar N12CritGlobal = computeCriticalGlobalShearLoad(
      D11Global, D22p, b, rho0Global, xiGlobal, gamma, zetaPanel);

  // backprop sensitivities from combined loading to the in-plane loads for
  // strain sens
  TacsScalar N1GlobalSens, N1CritGlobalSens, N12GlobalSens, N12CritGlobalSens;
  const TacsScalar strengthRatio = this->bucklingEnvelopeSens(
      -panelStress[0], N1CritGlobal, panelStress[2], N12CritGlobal,
      &N1GlobalSens, &N1CritGlobalSens, &N12GlobalSens, &N12CritGlobalSens);

  // backprop in-plane load sens (stress sens) to the DVs
  TacsScalar dfdPanelStress[NUM_STRESSES];
  memset(dfdPanelStress, 0, NUM_STRESSES * sizeof(TacsScalar));
  dfdPanelStress[0] = -N1GlobalSens;
  dfdPanelStress[2] = N12GlobalSens;
  // figure out whether this should be based on
  this->addPanelStressDVSens(scale, strain, dfdPanelStress,
                             &dfdx[this->panelDVStartNum]);

  // backprop critical buckling load sens to the material properties and the
  // non-dimensional parameters
  // --------------------------

  // holder for local DV jacobians (we then multiply them by scale when we add
  // them into the dfdx) DVsens format [0 - panel length, 1 - stiff pitch, 2 -
  // panel thick,
  //                3 - stiff height, 4 - stiff thick, 5 - panel width]
  TacsScalar DVsens[6];
  memset(DVsens, 0, 6 * sizeof(TacsScalar));

  // sensitivities for E1p, E1s, zn which we'll handle later
  TacsScalar E1p_bar = 0.0, E1s_bar = 0.0, zn_bar = 0.0;

  // define the material sensivitiies
  TacsScalar Dpsens[4];  // D11Global,D12,D22,D66
  TacsScalar Apsens[4];  // A11,A12,A22,A66
  memset(Dpsens, 0, 4 * sizeof(TacsScalar));
  memset(Apsens, 0, 4 * sizeof(TacsScalar));

  // define the non-dimensional (ND) sensitivities [rho0Global, xiGlobal, delta,
  // gamma, zeta]
  TacsScalar NDsens[5];
  memset(NDsens, 0, 5 * sizeof(TacsScalar));

  // backprop from the critical axial + shear local buckling loads to the
  // material, ND params, DVs
  computeCriticalGlobalAxialLoadSens(
      N1CritGlobalSens, D11Global, D22p, b, delta, rho0Global, xiGlobal, gamma,
      zetaPanel, &Dpsens[0], &Dpsens[2], &DVsens[5], &NDsens[2], &NDsens[0],
      &NDsens[1], &NDsens[3], &NDsens[4]);
  computeCriticalGlobalShearLoadSens(
      N12CritGlobalSens, D11Global, D22p, b, rho0Global, xiGlobal, gamma,
      zetaPanel, &Dpsens[0], &Dpsens[2], &DVsens[5], &NDsens[0], &NDsens[1],
      &NDsens[3], &NDsens[4]);

  // backprop from ND params to the material and DVs
  // -------------------

  computeAffineAspectRatioSens(NDsens[0], D11Global, D22p, a, b, &Dpsens[0],
                               &Dpsens[2], &DVsens[0],
                               &DVsens[5]);  // 0 - backprop from rho0Global
  computeLaminateIsotropySens(NDsens[1], D11Global, D22p, D12p, D66p,
                              &Dpsens[0], &Dpsens[2], &Dpsens[1],
                              &Dpsens[3]);  // 1 - backprop from xiGlobal
  computeStiffenerAreaRatioSens(NDsens[2], E1p, E1s, &DVsens[4], &DVsens[3],
                                &DVsens[1], &DVsens[2], &E1p_bar,
                                &E1s_bar);  // 2 - backprop from delta
  computeStiffenerStiffnessRatioSens(
      NDsens[3], D11Global, E1s, zn, &Dpsens[0], &DVsens[4], &DVsens[3],
      &DVsens[1], &E1s_bar, &zn_bar);  // 3 - backprop from gamma
  computeTransverseShearParameterSens(
      NDsens[4], A66p, A11p, b, h, &Apsens[3], &Apsens[0], &DVsens[5],
      &DVsens[2]);  // 4 - backprop from zetaPanel

  // backprop global buckling sensitivities
  // ---------------------------------
  // backprop D11p global panel buckling to E1p, E1s and DV sensitivities
  computePanelGlobalBucklingStiffnessSens(Dpsens[0], E1p, zn, &DVsens[1],
                                          &DVsens[2], &E1p_bar, &zn_bar);

  // backpropagate zn_bar (sens) to E1p, E1s sensitivities
  computeOverallCentroidSens(zn_bar, E1p, E1s, &DVsens[4], &DVsens[3],
                             &DVsens[2], &DVsens[1], &E1p_bar, &E1s_bar);

  // now backprop to panel ply fraction sensitivities through E1p_bar
  for (int plyNum = 0; plyNum < this->numPanelPlies; plyNum++) {
    int dvNum = this->panelPlyFracLocalNums[plyNum];
    if (dvNum >= 0) {
      const TacsScalar* Q = &(this->panelQMats[plyNum * NUM_Q_ENTRIES]);
      TacsScalar jac = (Q[0] - Q[1] * Q[1] / Q[3]);  // dE1p / dply_frac[i]
      dfdx[dvNum] += scale * E1p_bar * jac;
    }
  }

  // now backprop to stiffener ply fraction sensitivities through E1s_bar
  for (int plyNum = 0; plyNum < this->numStiffenerPlies; plyNum++) {
    int dvNum = this->stiffenerPlyFracLocalNums[plyNum];
    if (dvNum >= 0) {
      const TacsScalar* Q = &(this->stiffenerQMats[plyNum * NUM_Q_ENTRIES]);
      TacsScalar jac = (Q[0] - Q[1] * Q[1] / Q[3]);  // dE1p / dply_frac[i]
      dfdx[dvNum] += scale * E1s_bar * jac;
    }
  }

  // backprop from material sensitivities to DVs
  // --------------------

  // panel thickness term for material part..
  TacsScalar t = this->panelThick;
  if (this->panelThickNum >= 0) {
    int dvNum = this->panelThickLocalNum;
    // backpropagate through the D matrix
    TacsScalar dDfactor_dthick = 0.25 * t * t;  // d/dt(t^3/12) = t^2/4
    for (int ii = 0; ii < this->numPanelPlies; ii++) {
      TacsScalar* Q = &(this->panelQMats[ii * NUM_Q_ENTRIES]);
      dfdx[dvNum] += scale * dDfactor_dthick * this->panelPlyFracs[ii] *
                     (Dpsens[2] * Q[3] + Dpsens[1] * Q[1] +
                      Dpsens[3] * Q[5]);  // no Dpsens[0] here as global D11
      dfdx[dvNum] += scale * this->panelPlyFracs[ii] *
                     (Apsens[0] * Q[0] + Apsens[2] * Q[3] + Apsens[1] * Q[1] +
                      Apsens[3] * Q[5]);
    }
  }

  // --- Panel Ply fraction sensitivities for A,D matrices ---
  for (int plyNum = 0; plyNum < this->numPanelPlies; plyNum++) {
    int dvNum = this->panelPlyFracLocalNums[plyNum];
    if (dvNum >= 0) {
      const TacsScalar* Q = &(this->panelQMats[plyNum * NUM_Q_ENTRIES]);
      dfdx[dvNum] += scale * (t * t * t / 12.0) *
                     (Dpsens[2] * Q[3] + Dpsens[1] * Q[1] +
                      Dpsens[3] * Q[5]);  // no Dpsens[0] here as global D11
      dfdx[dvNum] += scale * t *
                     (Apsens[0] * Q[0] + Apsens[2] * Q[3] + Apsens[1] * Q[1] +
                      Apsens[3] * Q[5]);
    }
  }

  // combine local DV jacobians with scale to write into dfdx full derivatives
  // --------------------------------

  if (this->panelLengthLocalNum >= 0) {
    dfdx[this->panelLengthLocalNum] += scale * DVsens[0];
  }
  if (this->stiffenerPitchLocalNum >= 0) {
    dfdx[this->stiffenerPitchLocalNum] += scale * DVsens[1];
  }
  if (this->panelThickLocalNum >= 0) {
    dfdx[this->panelThickLocalNum] += scale * DVsens[2];
  }
  if (this->stiffenerHeightLocalNum >= 0) {
    dfdx[this->stiffenerHeightLocalNum] += scale * DVsens[3];
  }
  if (this->stiffenerThickLocalNum >= 0) {
    dfdx[this->stiffenerThickLocalNum] += scale * DVsens[4];
  }
  if (this->panelWidthLocalNum >= 0) {
    dfdx[this->panelWidthLocalNum] += scale * DVsens[5];
  }
}

void TACSGPBladeStiffenedShellConstitutive::addStiffenerCripplingDVSens(
    const TacsScalar scale, const TacsScalar stiffenerStrain[],
    TacsScalar dfdx[]) {
  if (CPTstiffenerCrippling) {  // use predictions for Sean's paper

    // previous section writes directly into dfdx, this section writes into
    // DVsens DVsens format [0 - panel length, 1 - stiff pitch, 2 - panel thick,
    //                3 - stiff height, 4 - stiff thick, 5 - panel width]
    TacsScalar DVsens[6];
    memset(DVsens, 0, 6 * sizeof(TacsScalar));

    // setup new ND parameter array for stiffener computation
    // ND parameter sens [xi, rho0, genPoiss, zeta]
    TacsScalar stiffNDsens[4];
    memset(stiffNDsens, 0, 4 * sizeof(TacsScalar));

    // set initial A,D matrix, nondim parameter and DV sensitivities to
    // backpropagate to.
    TacsScalar Ds_sens[4];  // D11,D12,D22,D66
    TacsScalar As_sens[4];  // A11,A12,A22,A66
    memset(Ds_sens, 0, 4 * sizeof(TacsScalar));
    memset(As_sens, 0, 4 * sizeof(TacsScalar));

    // compute D matrix of the stiffener (treating it like a panel for
    // crippling)
    TacsScalar stiffenerCripplingStiffness[NUM_TANGENT_STIFFNESS_ENTRIES];
    const TacsScalar *As_crippling, *Ds_crippling;
    this->computeStiffenerCripplingStiffness(stiffenerCripplingStiffness);
    this->extractTangentStiffness(stiffenerCripplingStiffness, &As_crippling,
                                  NULL, &Ds_crippling, NULL, NULL);

    // temporarily modify
    TacsScalar A11s = As_crippling[0], A66s = As_crippling[5];
    TacsScalar D11s = Ds_crippling[0], D12s = Ds_crippling[1];
    TacsScalar D22s = Ds_crippling[3], D66s = Ds_crippling[5];
    TacsScalar bStiff = this->stiffenerHeight;
    TacsScalar hStiff = this->stiffenerThick;
    TacsScalar a = this->panelLength;

    // compute stiffener non-dimensional parameters
    TacsScalar rho0Stiff = computeAffineAspectRatio(D11s, D22s, a, bStiff);
    TacsScalar xiStiff = computeLaminateIsotropy(D11s, D22s, D12s, D66s);
    TacsScalar genPoiss = computeGeneralizedPoissonsRatio(D12s, D66s);
    TacsScalar zetaStiff =
        computeTransverseShearParameter(A66s, A11s, bStiff, hStiff);

    // Compute stiffener in plane load and crippling failure index
    TacsScalar A11s_beam;
    TacsScalar N1stiff =
        computeStiffenerInPlaneLoad(stiffenerStrain, &A11s_beam);
    // TacsScalar N1stiff = 1.0;
    TacsScalar N1CritCrippling = computeStiffenerCripplingLoad(
        D11s, D22s, xiStiff, rho0Stiff, genPoiss, zetaStiff);
    TacsScalar fail_index = N1stiff / N1CritCrippling;
    // --- End of computeFailuresValues subsections ---

    // backpropagate from fail index to N1 stiffener in plane load
    TacsScalar stiffN1sens = 1.0 / N1CritCrippling;
    this->computeStiffenerInPlaneLoadSens(scale, stiffenerStrain, stiffN1sens,
                                          &dfdx[this->stiffenerDVStartNum]);

    // backpropagate N1crit sens of stiffener to ND and material sensitivities
    TacsScalar stiffN1critsens = fail_index * -1.0 / N1CritCrippling;
    computeStiffenerCripplingLoadSens(
        stiffN1critsens, D11s, D22s, xiStiff, rho0Stiff, genPoiss, zetaStiff,
        &Ds_sens[0], &Ds_sens[2], &DVsens[3], &stiffNDsens[0], &stiffNDsens[1],
        &stiffNDsens[2], &stiffNDsens[3]);

    // backpropagate ND sensitivities to A, D matrix for stiffener and DV sens
    computeLaminateIsotropySens(stiffNDsens[0], D11s, D22s, D12s, D66s,
                                &Ds_sens[0], &Ds_sens[2], &Ds_sens[1],
                                &Ds_sens[3]);
    computeAffineAspectRatioSens(stiffNDsens[1], D11s, D22s, a, bStiff,
                                 &Ds_sens[0], &Ds_sens[2], &DVsens[0],
                                 &DVsens[3]);
    computeGeneralizedPoissonsRatioSens(stiffNDsens[2], D12s, D66s, &Ds_sens[1],
                                        &Ds_sens[3]);
    computeTransverseShearParameterSens(stiffNDsens[3], A66s, A11s, bStiff,
                                        hStiff, &As_sens[3], &As_sens[0],
                                        &DVsens[3], &DVsens[4]);

    // backpropagate stiffener A,D matrix sensitivities through stiffener DV
    // -----------------------

    TacsScalar t_s = this->stiffenerThick;
    if (this->stiffenerThickLocalNum >= 0) {
      int dvNum = this->stiffenerThickLocalNum;

      TacsScalar dDfactor_dthick = 0.25 * t_s * t_s;  // d/dt(t^3/12) = t^2/4
      for (int ii = 0; ii < this->numStiffenerPlies; ii++) {
        TacsScalar* Q = &(this->stiffenerQMats[ii * NUM_Q_ENTRIES]);
        // backpropagate through the D matrix
        dfdx[dvNum] += scale * dDfactor_dthick * this->stiffenerPlyFracs[ii] *
                       (Ds_sens[0] * Q[0] + Ds_sens[2] * Q[3] +
                        Ds_sens[1] * Q[1] + Ds_sens[3] * Q[5]);
        // backpropagate through the A matrix
        dfdx[dvNum] += scale * this->stiffenerPlyFracs[ii] *
                       (As_sens[0] * Q[0] + As_sens[2] * Q[3] +
                        As_sens[1] * Q[1] + As_sens[3] * Q[5]);
      }
    }

    // backpropagate stiffener A,D matrix sens through stiff ply DVs
    // -----------------------
    for (int plyNum = 0; plyNum < this->numStiffenerPlies; plyNum++) {
      int dvNum = this->stiffenerPlyFracLocalNums[plyNum];
      if (dvNum >= 0) {
        const TacsScalar* Q = &(this->stiffenerQMats[plyNum * NUM_Q_ENTRIES]);
        dfdx[dvNum] += scale * (t_s * t_s * t_s / 12.0) *
                       (Ds_sens[0] * Q[0] + Ds_sens[2] * Q[3] +
                        Ds_sens[1] * Q[1] + Ds_sens[3] * Q[5]);
        dfdx[dvNum] += scale * t_s *
                       (As_sens[0] * Q[0] + As_sens[2] * Q[3] +
                        As_sens[1] * Q[1] + As_sens[3] * Q[5]);
      }
    }

    // 2,3,4 - backpropagate remaining DV sens into dfdx
    // --------------------------------------------------

    // recall DV sens [0 - panel length, 1 - stiff pitch, 2 - panel thick,
    //                 3 - stiff height, 4 - stiff thick, 5 - panel width]

    if (this->panelLengthLocalNum >= 0) {
      dfdx[this->panelLengthLocalNum] += scale * DVsens[0];
    }
    if (this->stiffenerPitchLocalNum >= 0) {
      dfdx[this->stiffenerPitchLocalNum] += scale * DVsens[1];
    }
    if (this->panelThickLocalNum >= 0) {
      dfdx[this->panelThickLocalNum] += scale * DVsens[2];
    }
    if (this->stiffenerHeightLocalNum >= 0) {
      dfdx[this->stiffenerHeightLocalNum] += scale * DVsens[3];
    }
    if (this->stiffenerThickLocalNum >= 0) {
      dfdx[this->stiffenerThickLocalNum] += scale * DVsens[4];
    }
    if (this->panelWidthLocalNum >= 0) {
      dfdx[this->panelWidthLocalNum] += scale * DVsens[5];
    }

  } else {  // CPTstiffenerCrippling = false
    // the DOD experimental stiffener crippling solution from Ali's superclass
    return TACSBladeStiffenedShellConstitutive::addStiffenerCripplingDVSens(
        scale, stiffenerStrain, dfdx);
  }
}

// ==============================================================================
// HELPER ROUTINES
// ==============================================================================
void TACSGPBladeStiffenedShellConstitutive::computePanelGlobalBucklingStiffness(
    const TacsScalar E1p, const TacsScalar zn, TacsScalar* D1) {
  TacsScalar Ap, Ip, tp, sp, tp3;
  tp = this->panelThick;
  tp3 = tp * tp * tp;
  sp = this->stiffenerPitch;
  Ap = tp * sp;
  Ip = sp * tp3 / 12.0;

  // compute panel D11 with overall centroid from panel and stiffener (so don't
  // include stiffener D11 here)
  *D1 = (E1p * (Ip + Ap * zn * zn)) / sp;

  // don't shift centroid of D22 to be consistent with stiffness matrix

  // ignoring this value for now => should maybe add this term back in later and
  // compare results?
  // // --- Twisting stiffness ---
  // // Compute the shear modulus weighted centroid of the panel and stiffener
  // TacsScalar zg = 0.25 * QStiffener[5] * As * (-0.5 * tp + zs) /
  //                 (QStiffener[5] * As + QPanel[5] * Ap);
  // *D3 = (QPanel[5] * (Jp + Ap * zg * zg)) / ps;
}

void TACSGPBladeStiffenedShellConstitutive::
    computePanelGlobalBucklingStiffnessSens(
        const TacsScalar D1Sens, const TacsScalar E1p, const TacsScalar zn,
        TacsScalar* spitchSens, TacsScalar* pthickSens, TacsScalar* E1pSens,
        TacsScalar* znSens) {
  TacsScalar Ap, Ip, tp, sp, tp3;
  tp = this->panelThick;
  tp3 = tp * tp * tp;
  sp = this->stiffenerPitch;
  Ap = tp * sp;
  Ip = sp * tp3 / 12.0;

  TacsScalar D1;
  computePanelGlobalBucklingStiffness(E1p, zn, &D1);

  // now backpropagate derivatives to E1p, zn
  *E1pSens += D1Sens * D1 / E1p;
  *znSens += D1Sens * E1p * Ap * 2.0 * zn / sp;

  // backprop through intermediate vars
  TacsScalar IpSens = D1Sens * E1p / sp;
  *spitchSens += IpSens * tp3 / 12.0;
  *pthickSens += IpSens * Ip * 3.0 / tp;

  TacsScalar ApSens = D1Sens * E1p * zn * zn / sp;
  *spitchSens += ApSens * tp;
  *pthickSens += ApSens * sp;

  // backprop to DVs
  *spitchSens += D1Sens * D1 * -1.0 / sp;
}

// Retrieve the design variable for plotting purposes
TacsScalar TACSGPBladeStiffenedShellConstitutive::evalDesignFieldValue(
    int elemIndex, const double pt[], const TacsScalar X[], int index) {
  if (writeDVmode == 0) {
    switch (index) {
      case 0:
        return this->computeEffectiveThickness();
      case 1:
        // return this->computeEffectiveBendingThickness();
        return this->panelWidth;
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
  } else if (writeDVmode == 1) {
    // view the non-dimensional parameters

    TacsScalar panelStiffness[NUM_TANGENT_STIFFNESS_ENTRIES];
    this->computePanelStiffness(panelStiffness);
    const TacsScalar *Ap, *Dp;
    TacsScalar D11p;  //, D66p;
    this->extractTangentStiffness(panelStiffness, &Ap, NULL, &Dp, NULL, NULL);
    // compute effective moduli, overall centroid
    TacsScalar E1s, E1p, _;
    this->computeEffectiveModulii(this->numPanelPlies, this->panelQMats,
                                  this->panelPlyFracs, &E1p, &_);
    this->computeEffectiveModulii(this->numStiffenerPlies, this->stiffenerQMats,
                                  this->stiffenerPlyFracs, &E1s, &_);
    TacsScalar zn = this->computeOverallCentroid(E1p, E1s);

    computePanelGlobalBucklingStiffness(E1p, zn, &D11p);

    // Compute panel dimensions, material props and non-dimensional parameters
    TacsScalar N1CritGlobal, N12CritGlobal;
    // TacsScalar D11p = Dp[0];
    TacsScalar D12p = Dp[1], D66p = Dp[5], D22p = Dp[3];
    TacsScalar A11p = Ap[0], A66p = Ap[5];
    TacsScalar delta, rho0, xi, gamma, a, b, zeta;
    a = this->panelLength;
    b = this->panelWidth;
    delta = computeStiffenerAreaRatio(E1p, E1s);
    rho0 = computeAffineAspectRatio(D11p, D22p, a, b);
    xi = computeLaminateIsotropy(D11p, D22p, D12p, D66p);
    gamma = computeStiffenerStiffnessRatio(D11p, E1s, zn);
    zeta = computeTransverseShearParameter(A66p, A11p, b, this->panelThick);

    TacsScalar SAR = this->stiffenerHeight / this->stiffenerThick;

    switch (index) {
      case 0:
        return this->stiffenerPitch;
      case 1:
        return SAR;
      case 2:
        return delta;
      case 3:
        return rho0;
      case 4:
        return xi;
      case 5:
        return gamma;
      case 6:
        return zeta;
    }
  }

  return 0.0;
}

void TACSGPBladeStiffenedShellConstitutive::computeStiffenerCripplingStiffness(
    TacsScalar C[]) {
  TacsScalar* A = &C[0];
  TacsScalar* B = &C[6];
  TacsScalar* D = &C[12];
  TacsScalar* As = &C[18];

  // --- Zero out the C matrix ---
  memset(C, 0, this->NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));

  // Compute the smeared laminate properties
  TacsScalar QStiff[this->NUM_Q_ENTRIES], ABarStiff[this->NUM_ABAR_ENTRIES];

  this->computeSmearedStiffness(this->numStiffenerPlies, this->stiffenerQMats,
                                this->stiffenerAbarMats,
                                this->stiffenerPlyFracs, QStiff, ABarStiff);

  // Add the panel's contributions to the A and D matrices
  TacsScalar t = this->stiffenerThick;
  TacsScalar DFactor = t * t * t / 12.0;

  for (int ii = 0; ii < NUM_Q_ENTRIES; ii++) {
    A[ii] += t * QStiff[ii];
    D[ii] += DFactor * QStiff[ii];
  }

  // Add the pane;'s contribution to the transverse shear matrix
  for (int ii = 0; ii < NUM_ABAR_ENTRIES; ii++) {
    As[ii] += t * QStiff[ii] * this->kcorr;
  }

  // Add the drill stiffness
  C[21] = DRILLING_REGULARIZATION * 0.5 * (As[0] + As[2]);
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeStiffenerInPlaneLoad(
    const TacsScalar stiffenerStrain[], TacsScalar* A11s) {
  int n = TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES;
  TacsScalar C[n];
  memset(C, 0, n * sizeof(TacsScalar));
  this->computeStiffenerStiffness(C);

  // divide out the height (width of stiffener when viewed as vertical panel)
  // this way N11 = load/width = stress * thickness is still satisfied
  *A11s = C[0] / this->stiffenerHeight;

  // return compressive strain * the A11 plate constant of the stiffener
  return -stiffenerStrain[0] * (*A11s);
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::computeStiffenerInPlaneLoadSens(
    const TacsScalar scale, const TacsScalar stiffenerStrain[],
    const TacsScalar dN11_stiff, TacsScalar dfdx[]) {
  int n = TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES;
  TacsScalar C[n];
  memset(C, 0, n * sizeof(TacsScalar));
  this->computeStiffenerStiffness(C);

  // get intermediate quantities needed for jacobian dN11_stiff/dx
  TacsScalar E, _;
  this->computeEffectiveModulii(this->numStiffenerPlies, this->stiffenerQMats,
                                this->stiffenerPlyFracs, &E, &_);

  // backpropgate to N11_bar = df/dN11_stiff
  TacsScalar N11_stiff = this->computeStiffenerInPlaneLoad(stiffenerStrain, &_);
  TacsScalar N11_bar = scale * dN11_stiff;

  // derivatives through E the effective modulus
  //    includes stiffener ply fraction DVs only
  for (int plyNum = 0; plyNum < this->numStiffenerPlies; plyNum++) {
    int dvNum =
        this->stiffenerPlyFracLocalNums[plyNum] - this->stiffenerDVStartNum;
    if (dvNum >= 0) {
      const TacsScalar* Q = &(this->stiffenerQMats[plyNum * NUM_Q_ENTRIES]);
      TacsScalar ply_frac_jac =
          (Q[0] - Q[1] * Q[1] / Q[3]);  // dE1p / dply_frac[i]
      TacsScalar E1s_jac = N11_stiff / E;
      dfdx[dvNum] += N11_bar * E1s_jac * ply_frac_jac;
    }
  }

  // derivatives through A11s but not effective modulus part
  //    includes stiffener thick (height cancels out b.c. in-plane load)
  if (this->stiffenerThickLocalNum >= 0) {
    int dvNum = this->stiffenerThickLocalNum - this->stiffenerDVStartNum;
    // use power-series rule since N11_stiff propto sthick^1
    // if f = A * x, then df/dx = A = f / x (useful trick for fewer steps)
    TacsScalar sthick_jac = N11_stiff / this->stiffenerThick;
    dfdx[dvNum] += N11_bar * sthick_jac;
  }

  // no longer consider derivatives directly for the strain here.

  // return usual forward analysis output if desired
  return N11_stiff;
}

// Retrieve the global design variable numbers
int TACSGPBladeStiffenedShellConstitutive::getDesignVarNums(int elemIndex,
                                                            int dvLen,
                                                            int dvNums[]) {
  TACSBladeStiffenedShellConstitutive::getDesignVarNums(elemIndex, dvLen,
                                                        dvNums);
  if (dvNums && dvLen >= this->numDesignVars) {
    if (this->panelWidthNum >= 0) {
      dvNums[this->panelWidthLocalNum] = panelWidthNum;
    }
  }
  return numDesignVars;
}

// Set the element design variable from the design vector
int TACSGPBladeStiffenedShellConstitutive::setDesignVars(
    int elemIndex, int dvLen, const TacsScalar dvs[]) {
  TACSBladeStiffenedShellConstitutive::setDesignVars(elemIndex, dvLen, dvs);
  if (dvLen >= this->numDesignVars) {
    if (this->panelWidthNum >= 0) {
      this->panelWidth = dvs[this->panelWidthLocalNum];
    }
  }

  // NOTE : this is a very important step => here we reset the save on all
  // computed GP models so they recalculate and compute their new values.
  if (this->panelGPs) {
    this->panelGPs->resetSavedData();
  }

  return this->numDesignVars;
}

// Get the element design variables values
int TACSGPBladeStiffenedShellConstitutive::getDesignVars(int elemIndex,
                                                         int dvLen,
                                                         TacsScalar dvs[]) {
  TACSBladeStiffenedShellConstitutive::getDesignVars(elemIndex, dvLen, dvs);
  if (dvLen >= this->numDesignVars) {
    if (this->panelWidthNum >= 0) {
      dvs[this->panelWidthLocalNum] = this->panelWidth;
    }
  }
  return this->numDesignVars;
}

// Get the lower and upper bounds for the design variable values
int TACSGPBladeStiffenedShellConstitutive::getDesignVarRange(int elemIndex,
                                                             int dvLen,
                                                             TacsScalar lb[],
                                                             TacsScalar ub[]) {
  TACSBladeStiffenedShellConstitutive::getDesignVarRange(elemIndex, dvLen, lb,
                                                         ub);
  if (dvLen >= this->numDesignVars) {
    if (this->panelWidthNum >= 0) {
      lb[this->panelWidthLocalNum] = this->panelWidthLowerBound;
      ub[this->panelWidthLocalNum] = this->panelWidthUpperBound;
    }
  }
  return this->numDesignVars;
}

// ==============================================================================
// Buckling functions
// ==============================================================================

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeAffineAspectRatioSens(
    const TacsScalar rho0sens, const TacsScalar D11, const TacsScalar D22,
    const TacsScalar a, const TacsScalar b, TacsScalar* D11sens,
    TacsScalar* D22sens, TacsScalar* asens, TacsScalar* bsens) {
  // compute the derivatives of the affine aspect ratio and return the affine
  // aspect ratio
  TacsScalar rho_0 = computeAffineAspectRatio(D11, D22, a, b);
  // where rho_0 = a/b * (D22/D11)**0.25

  // use power series rules d(x^p) = p * (x^p) / x to cleanly differentiate the
  // expression
  *D11sens += rho0sens * rho_0 * -0.25 / D11;
  *D22sens += rho0sens * rho_0 * 0.25 / D22;
  *asens += rho0sens * rho_0 / a;
  *bsens += rho0sens * rho_0 * -1.0 / b;

  return rho_0;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeLaminateIsotropySens(
    const TacsScalar xisens, const TacsScalar D11, const TacsScalar D22,
    const TacsScalar D12, const TacsScalar D66, TacsScalar* D11sens,
    TacsScalar* D22sens, TacsScalar* D12sens, TacsScalar* D66sens) {
  // compute the derivatives of the laminate isotropy xi
  TacsScalar denominator = sqrt(D11 * D22);
  TacsScalar xi = computeLaminateIsotropy(D11, D22, D12, D66);
  // so that xi = (D12 + 2 * D66) / sqrt(D11*D22)

  // compute the sensitivities
  *D12sens += xisens * 1.0 / denominator;
  *D66sens += xisens * 2.0 / denominator;
  *D11sens += xisens * -0.5 * xi / D11;
  *D22sens += xisens * -0.5 * xi / D22;

  return xi;
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::computeGeneralizedPoissonsRatioSens(
    const TacsScalar epssens, const TacsScalar D12, const TacsScalar D66,
    TacsScalar* D12sens, TacsScalar* D66sens) {
  // compute derivatives of the generalized poisson's ratio
  TacsScalar eps = computeGeneralizedPoissonsRatio(D12, D66);
  // where eps = (D12 + 2 * D66) / D12

  *D12sens += epssens * eps * eps / D12 / D12 * 2.0 * D66;
  *D66sens += epssens * -eps * eps * 2.0 / D12;

  return eps;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeStiffenerAreaRatio(
    const TacsScalar E1p, const TacsScalar E1s) {
  TacsScalar As = this->computeStiffenerArea();
  return E1s * As / (E1p * this->stiffenerPitch * this->panelThick);
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeStiffenerAreaRatioSens(
    const TacsScalar deltasens, const TacsScalar E1p, const TacsScalar E1s,
    TacsScalar* sthickSens, TacsScalar* sheightSens, TacsScalar* spitchSens,
    TacsScalar* pthickSens, TacsScalar* E1psens, TacsScalar* E1ssens) {
  TacsScalar delta = this->computeStiffenerAreaRatio(E1p, E1s);

  *sthickSens += deltasens * delta / this->stiffenerThick;
  *sheightSens += deltasens * delta / this->stiffenerHeight;
  *spitchSens += deltasens * -1.0 * delta / this->stiffenerPitch;
  *pthickSens += deltasens * -1.0 * delta / this->panelThick;

  *E1ssens += deltasens * delta / E1s;
  *E1psens += deltasens * delta * -1.0 / E1p;

  return delta;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeOverallCentroid(
    const TacsScalar E1p, const TacsScalar E1s) {
  TacsScalar Ap, As, zs, tp, ps;
  tp = this->panelThick;
  ps = this->stiffenerPitch;
  Ap = tp * ps;
  As = this->computeStiffenerArea();
  zs = this->computeStiffenerCentroidHeight();

  TacsScalar zn = E1s * As * (-0.5 * tp - zs) / (E1s * As + E1p * Ap);
  return zn;
}

void TACSGPBladeStiffenedShellConstitutive::computeOverallCentroidSens(
    const TacsScalar znSens, const TacsScalar E1p, const TacsScalar E1s,
    TacsScalar* sthickSens, TacsScalar* sheightSens, TacsScalar* pthickSens,
    TacsScalar* spitchSens, TacsScalar* E1pSens, TacsScalar* E1sSens) {
  // compute numerator and denominators
  TacsScalar Ap, As, zs, tp, sp;
  tp = this->panelThick;
  sp = this->stiffenerPitch;
  Ap = tp * sp;
  As = this->computeStiffenerArea();
  zs = this->computeStiffenerCentroidHeight();

  TacsScalar num = E1s * As * (-0.5 * tp - zs);
  TacsScalar den = E1s * As + E1p * Ap;
  TacsScalar zn = num / den;

  // backprop to num, den
  TacsScalar numSens = znSens / den;
  TacsScalar denSens = znSens * zn * -1.0 / den;

  // backprop to any direct DVs or output sens
  *E1sSens += numSens * num / E1s + denSens * As;
  *E1pSens += denSens * Ap;
  *pthickSens += numSens * E1s * As * -0.5;

  // backprop to intermediate quantities
  TacsScalar ApSens = denSens * E1p;
  TacsScalar AsSens = numSens * num / As + denSens * E1s;
  TacsScalar zsSens = numSens * -1.0 * E1s * As;

  // backprop to DVs from intermediate quantities
  *pthickSens += ApSens * sp;
  *spitchSens += ApSens * tp;

  TacsScalar dAsdt, dAsdh;
  computeStiffenerAreaSens(dAsdt, dAsdh);
  *sthickSens += AsSens * dAsdt;
  *sheightSens += AsSens * dAsdh;

  TacsScalar dzsdt, dzsdh;
  computeStiffenerCentroidHeightSens(dzsdt, dzsdh);
  *sthickSens += zsSens * dzsdt;
  *sheightSens += zsSens * dzsdh;
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::computeStiffenerStiffnessRatio(
    const TacsScalar D11, const TacsScalar E1s, const TacsScalar zn) {
  TacsScalar As = this->computeStiffenerArea();
  TacsScalar Is = this->computeStiffenerIzz();
  TacsScalar zs =
      -1.0 * this->computeStiffenerCentroidHeight();  // stiffener down into
                                                      // wing -Z direction

  TacsScalar Is_total = (Is + As * (zn - zs) * (zn - zs));
  return E1s * Is_total / D11 / this->stiffenerPitch;
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::computeStiffenerStiffnessRatioSens(
    const TacsScalar gammaSens, const TacsScalar D11, const TacsScalar E1s,
    const TacsScalar zn, TacsScalar* D11Sens, TacsScalar* sthickSens,
    TacsScalar* sheightSens, TacsScalar* spitchSens, TacsScalar* E1sSens,
    TacsScalar* znSens) {
  // use power series rules and the forward state to differentiate
  TacsScalar gamma = computeStiffenerStiffnessRatio(D11, E1s, zn);

  // backpropagate to DV level or base level without the Is term which is
  // hardest
  *D11Sens += gammaSens * -1.0 * gamma / D11;
  *spitchSens += gammaSens * gamma * -1.0 / this->stiffenerPitch;
  *E1sSens += gammaSens * gamma / E1s;

  // compute total Is
  TacsScalar As = this->computeStiffenerArea();
  TacsScalar Is = this->computeStiffenerIzz();
  TacsScalar zs =
      -1.0 * this->computeStiffenerCentroidHeight();  // stiffener down into
                                                      // wing -Z direction

  TacsScalar Is_total = (Is + As * (zn - zs) * (zn - zs));

  // backpropagation to Is_tot
  TacsScalar IstotSens = gammaSens * gamma / Is_total;

  // backpropagate to intermediate quantities As, Is, zs, zn
  TacsScalar AsSens = IstotSens * (zn - zs) * (zn - zs);
  TacsScalar IsSens = IstotSens;
  TacsScalar zsSens = -1.0 * IstotSens * As * 2.0 * (zs - zn);
  *znSens += IstotSens * As * 2.0 * (zn - zs);

  // backpropgate from intermediate quantities to DVs or material, etc.
  TacsScalar dAdt, dAdh;
  this->computeStiffenerAreaSens(dAdt, dAdh);
  *sthickSens += AsSens * dAdt;
  *sheightSens += AsSens * dAdh;

  TacsScalar dIdt, dIdh;
  this->computeStiffenerIzzSens(dIdt, dIdh);
  *sthickSens += IsSens * dIdt;
  *sheightSens += IsSens * dIdh;

  TacsScalar dZdt, dZdh;
  this->computeStiffenerCentroidHeightSens(dZdt, dZdh);
  *sthickSens += zsSens * dZdt;
  *sheightSens += zsSens * dZdh;

  return gamma;
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::computeTransverseShearParameter(
    TacsScalar A66, TacsScalar A11, TacsScalar b, TacsScalar h) {
  return A11 / A66 * (h / b) * (h / b);
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::computeTransverseShearParameterSens(
    const TacsScalar zetasens, const TacsScalar A66, const TacsScalar A11,
    const TacsScalar b, const TacsScalar h, TacsScalar* A66sens,
    TacsScalar* A11sens, TacsScalar* bsens, TacsScalar* hsens) {
  TacsScalar zeta = computeTransverseShearParameter(A66, A11, b, h);
  TacsScalar dzeta = zetasens * zeta;

  *A66sens += dzeta * -1.0 / A66;
  *A11sens += dzeta / A11;
  *bsens += dzeta * -2.0 / b;
  *hsens += dzeta * 2.0 / h;

  return zeta;
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::computeCriticalGlobalAxialLoad(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar b,
    const TacsScalar delta, const TacsScalar rho_0, const TacsScalar xi,
    const TacsScalar gamma, const TacsScalar zeta) {
  if (this->getAxialGP()) {
    // use Gaussian processes to compute the critical global axial load
    TacsScalar dim_factor =
        M_PI * M_PI * sqrt(D11 * D22) / b / b / (1.0 + delta);
    TacsScalar one = 1.0;
    XtestAxial[0] = log(one + xi);
    XtestAxial[1] = log(rho_0);
    XtestAxial[2] = log(one + gamma);
    XtestAxial[3] = log(one + 1000.0 * zeta);
    TacsScalar nd_factor =
        exp(this->panelGPs->predictMeanTestData(0, XtestAxial));
    return dim_factor * nd_factor;

  } else {
    // use the CPT closed-form solution to compute the critical global axial
    // load
    TacsScalar neg_N11crits[this->NUM_CF_MODES];
    for (int _m1 = 1; _m1 < this->NUM_CF_MODES + 1; _m1++) {
      TacsScalar dim_factor =
          M_PI * M_PI * sqrt(D11 * D22) / b / b / (1.0 + delta);
      TacsScalar m1 = _m1;
      TacsScalar nondim_factor = (1.0 + gamma) * pow(m1 / rho_0, 2.0) +
                                 pow(m1 / rho_0, -2.0) + 2.0 * xi;
      neg_N11crits[_m1 - 1] =
          -1.0 * dim_factor * nondim_factor;  // negated only because we have to
                                              // do KS min aggregate later
    }

    // compute KS aggregation for -N11crit for each mode then negate again
    // (because we want minimum N11crit so maximize negative N11crit)
    TacsScalar neg_N11crit =
        ksAggregation(neg_N11crits, this->NUM_CF_MODES, this->ksWeight);

    TacsScalar dim_factor =
        M_PI * M_PI * sqrt(D11 * D22) / b / b / (1.0 + delta);
    TacsScalar N11crStar = -1.0 * neg_N11crit / dim_factor;
    // printf("N11cr* = %.4f\n", N11crStar);
    return -1.0 * neg_N11crit;
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::nondimCriticalGlobalAxialLoad(
    const TacsScalar rho_0, const TacsScalar xi, const TacsScalar gamma,
    const TacsScalar zeta) {
  if (this->getAxialGP()) {
    // use Gaussian processes to compute the critical global axial load
    TacsScalar dim_factor = 1.0;
    TacsScalar one = 1.0;
    XtestAxial[0] = log(one + xi);
    XtestAxial[1] = log(rho_0);
    XtestAxial[2] = log(one + gamma);
    XtestAxial[3] = log(one + 1000.0 * zeta);

    // don't need to use the saved data here since this routine is only meant to
    // be called on a single constitutive object (not a whole mesh with O(1e4)
    // const objects)
    TacsScalar nd_factor =
        exp(this->getAxialGP()->predictMeanTestData(XtestAxial));
    return dim_factor * nd_factor;

  } else {
    // use the CPT closed-form solution to compute the critical global axial
    // load
    TacsScalar neg_N11crits[this->NUM_CF_MODES];
    for (int _m1 = 1; _m1 < this->NUM_CF_MODES + 1; _m1++) {
      TacsScalar dim_factor = 1.0;
      TacsScalar m1 = _m1;
      TacsScalar nondim_factor = (1.0 + gamma) * pow(m1 / rho_0, 2.0) +
                                 pow(m1 / rho_0, -2.0) + 2.0 * xi;
      neg_N11crits[_m1 - 1] =
          -1.0 * dim_factor * nondim_factor;  // negated only because we have to
                                              // do KS min aggregate later
    }

    // compute KS aggregation for -N11crit for each mode then negate again
    // (because we want minimum N11crit so maximize negative N11crit)
    TacsScalar neg_N11crit =
        ksAggregation(neg_N11crits, this->NUM_CF_MODES, this->ksWeight);
    return -1.0 * neg_N11crit;
  }
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::computeCriticalGlobalAxialLoadSens(
    const TacsScalar N1sens, const TacsScalar D11, const TacsScalar D22,
    const TacsScalar b, const TacsScalar delta, const TacsScalar rho_0,
    const TacsScalar xi, const TacsScalar gamma, const TacsScalar zeta,
    TacsScalar* D11sens, TacsScalar* D22sens, TacsScalar* bsens,
    TacsScalar* deltasens, TacsScalar* rho_0sens, TacsScalar* xisens,
    TacsScalar* gammasens, TacsScalar* zetasens) {
  if (this->getAxialGP()) {
    // use Gaussian processes to compute the critical global axial load
    TacsScalar dim_factor =
        M_PI * M_PI * sqrt(D11 * D22) / b / b / (1.0 + delta);
    TacsScalar one = 1.0;
    XtestAxial[0] = log(one + xi);
    XtestAxial[1] = log(rho_0);
    XtestAxial[2] = log(one + gamma);
    XtestAxial[3] = log(one + 1000.0 * zeta);
    TacsScalar arg = this->panelGPs->predictMeanTestData(0, XtestAxial);
    TacsScalar nondim_factor = exp(arg);
    TacsScalar output = dim_factor * nondim_factor;

    // compute sensitivities backwards propagated out of the GP
    // back to the nondim parameter inputs, (this part differentiates the
    // nondim_factor)
    TacsScalar Ysens = N1sens * output;
    this->panelGPs->predictMeanTestDataSens(0, Ysens, XtestAxial,
                                            XtestAxialSens);
    *xisens +=
        XtestAxialSens[0] / (one + xi);  // chain rule dlog(one + xi)/dxi = 1/xi
    *rho_0sens += XtestAxialSens[1] / rho_0;
    *gammasens += XtestAxialSens[2] / (1.0 + gamma);
    *zetasens += XtestAxialSens[3] / (one + 1000.0 * zeta) * 1000.0;

    // compute the sensivities of inputs in the dimensional constant
    // (this part differentiates the dim factor)
    *D11sens += Ysens * 0.5 / D11;
    *D22sens += Ysens * 0.5 / D22;
    *bsens += Ysens * -2.0 / b;
    *deltasens += Ysens * -1.0 / (1.0 + delta);

    return output;

  } else {
    // use the CPT closed-form solution to compute the critical global axial
    // load forward analysis part here
    TacsScalar neg_N11crits[this->NUM_CF_MODES];
    for (int _m1 = 1; _m1 < this->NUM_CF_MODES + 1; _m1++) {
      TacsScalar dim_factor =
          M_PI * M_PI * sqrt(D11 * D22) / b / b / (1.0 + delta);
      TacsScalar m1 = _m1;
      TacsScalar nondim_factor = (1.0 + gamma) * pow(m1 / rho_0, 2.0) +
                                 pow(m1 / rho_0, -2.0) + 2.0 * xi;
      neg_N11crits[_m1 - 1] =
          -1.0 * dim_factor * nondim_factor;  // negated only because we have to
                                              // do KS min aggregate later
    }

    // compute KS aggregation sensitivity
    TacsScalar neg_N11crits_sens[this->NUM_CF_MODES];
    TacsScalar neg_N11crit = ksAggregationSens(
        neg_N11crits, this->NUM_CF_MODES, this->ksWeight, neg_N11crits_sens);

    // compute sensitivities here
    for (int _m1 = 1; _m1 < this->NUM_CF_MODES + 1; _m1++) {
      // apply output sens
      neg_N11crits_sens[_m1 - 1] *= N1sens;
      // apply -1 last step
      neg_N11crits_sens[_m1 - 1] *= -1.0;

      // convert to double/cmplx type here
      TacsScalar m1 = _m1;

      // forward analysis states
      TacsScalar dim_factor =
          M_PI * M_PI * sqrt(D11 * D22) / b / b / (1.0 + delta);
      TacsScalar nondim_factor = (1.0 + gamma) * pow(m1 / rho_0, 2.0) +
                                 pow(m1 / rho_0, -2.0) + 2.0 * xi;
      neg_N11crits[_m1 - 1] =
          -1.0 * dim_factor * nondim_factor;  // negated only because we have to
                                              // do KS min aggregate later

      // update sensitivities (left factor is dKS/dv_i, right factor is dv_i /
      // dx)
      *D11sens +=
          neg_N11crits_sens[_m1 - 1] * (0.5 * neg_N11crits[_m1 - 1] / D11);
      *D22sens +=
          neg_N11crits_sens[_m1 - 1] * (0.5 * neg_N11crits[_m1 - 1] / D22);
      *bsens += neg_N11crits_sens[_m1 - 1] * (-2.0 * neg_N11crits[_m1 - 1] / b);
      *deltasens += neg_N11crits_sens[_m1 - 1] *
                    (-1.0 * neg_N11crits[_m1 - 1] / (1.0 + delta));
      *rho_0sens += neg_N11crits_sens[_m1 - 1] * -dim_factor *
                    ((1.0 + gamma) * -2.0 * pow(m1 / rho_0, 2.0) / rho_0 +
                     pow(m1 / rho_0, -2.0) * 2.0 / rho_0);
      *xisens += neg_N11crits_sens[_m1 - 1] * -dim_factor * 2.0;
      *gammasens +=
          neg_N11crits_sens[_m1 - 1] * -dim_factor * pow(m1 / rho_0, 2.0);
    }
    return -1.0 * neg_N11crit;
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeCriticalLocalAxialLoad(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar rho_0,
    const TacsScalar xi, const TacsScalar zeta) {
  if (this->getAxialGP()) {
    // use Gaussian processes to compute the critical global axial load
    TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) /
                            this->stiffenerPitch / this->stiffenerPitch;
    TacsScalar one = 1.0;
    XtestAxial[0] = log(one + xi);
    XtestAxial[1] = log(rho_0);
    XtestAxial[2] =
        0.0;  // log(1+gamma) = 0 since gamma=0 for unstiffened panel
    XtestAxial[3] = log(one + 1000.0 * zeta);
    TacsScalar nd_factor =
        exp(this->panelGPs->predictMeanTestData(1, XtestAxial));
    return dim_factor * nd_factor;

  } else {
    // use the CPT closed-form solution to compute the critical global axial
    // load
    TacsScalar neg_N11crits[this->NUM_CF_MODES];
    for (int _m1 = 1; _m1 < this->NUM_CF_MODES + 1; _m1++) {
      TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) /
                              this->stiffenerPitch / this->stiffenerPitch;
      // convert to double/cmplx type here
      TacsScalar m1 = _m1;
      TacsScalar nondim_factor =
          pow(m1 / rho_0, 2.0) + pow(m1 / rho_0, -2.0) + 2.0 * xi;
      neg_N11crits[_m1 - 1] =
          -1.0 * dim_factor * nondim_factor;  // negated only because we have to
                                              // do KS min aggregate later
    }

    // compute KS aggregation for -N11crit for each mode then negate again
    // (because we want minimum N11crit so maximize negative N11crit)
    TacsScalar neg_N11crit =
        ksAggregation(neg_N11crits, this->NUM_CF_MODES, this->ksWeight);
    return -1.0 * neg_N11crit;
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::nondimCriticalLocalAxialLoad(
    const TacsScalar rho_0, const TacsScalar xi, const TacsScalar zeta) {
  if (this->getAxialGP()) {
    // use Gaussian processes to compute the critical global axial load
    TacsScalar dim_factor = 1.0;
    TacsScalar one = 1.0;
    XtestAxial[0] = log(one + xi);
    XtestAxial[1] = log(rho_0);
    XtestAxial[2] =
        0.0;  // log(1+gamma) = 0 since gamma=0 for unstiffened panel
    XtestAxial[3] = log(one + 1000.0 * zeta);
    // don't need to use the saved data here since this routine is only meant to
    // be called on a single constitutive object (not a whole mesh with O(1e4)
    // const objects)
    TacsScalar nd_factor =
        exp(this->getAxialGP()->predictMeanTestData(XtestAxial));
    return dim_factor * nd_factor;

  } else {
    // use the CPT closed-form solution to compute the critical global axial
    // load
    TacsScalar neg_N11crits[this->NUM_CF_MODES];
    for (int _m1 = 1; _m1 < this->NUM_CF_MODES + 1; _m1++) {
      TacsScalar dim_factor = 1.0;
      // convert to double/cmplx type here
      TacsScalar m1 = _m1;
      TacsScalar nondim_factor =
          pow(m1 / rho_0, 2.0) + pow(m1 / rho_0, -2.0) + 2.0 * xi;
      neg_N11crits[_m1 - 1] =
          -1.0 * dim_factor * nondim_factor;  // negated only because we have to
                                              // do KS min aggregate later
    }

    // compute KS aggregation for -N11crit for each mode then negate again
    // (because we want minimum N11crit so maximize negative N11crit)
    TacsScalar neg_N11crit =
        ksAggregation(neg_N11crits, this->NUM_CF_MODES, this->ksWeight);
    return -1.0 * neg_N11crit;
  }
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::computeCriticalLocalAxialLoadSens(
    const TacsScalar N1sens, const TacsScalar D11, const TacsScalar D22,
    const TacsScalar rho_0, const TacsScalar xi, const TacsScalar zeta,
    TacsScalar* D11sens, TacsScalar* D22sens, TacsScalar* spitchsens,
    TacsScalar* rho_0sens, TacsScalar* xisens, TacsScalar* zetasens) {
  if (this->getAxialGP()) {
    // use Gaussian processes to compute the critical global axial load
    TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) /
                            this->stiffenerPitch / this->stiffenerPitch;
    TacsScalar one = 1.0;
    XtestAxial[0] = log(one + xi);
    XtestAxial[1] = log(rho_0);
    XtestAxial[2] =
        0.0;  // log(1+gamma) = 0 since gamma=0 for unstiffened panel
    XtestAxial[3] = log(one + 1000.0 * zeta);

    TacsScalar arg = this->panelGPs->predictMeanTestData(1, XtestAxial);
    TacsScalar nondim_factor = exp(arg);
    TacsScalar output = dim_factor * nondim_factor;

    // backwards propagate sensitivities out of the axialGP model to nondim
    // params, (this part differentiates the nondim_factor)
    TacsScalar Ysens = N1sens * output;
    this->panelGPs->predictMeanTestDataSens(1, Ysens, XtestAxial,
                                            XtestAxialSens);

    *xisens += XtestAxialSens[0] / (one + xi);
    *rho_0sens +=
        XtestAxialSens[1] / rho_0;  // chain rule dlog(rho_0) / drho_0 = 1/rho_0
    *zetasens += XtestAxialSens[3] / (one + 1000.0 * zeta) * 1000.0;

    // backpropagate the dimensional factor terms out to the material and
    // geometric DVs (this part differentiates the dim_factor)
    *D11sens += Ysens * 0.5 / D11;
    *D22sens += Ysens * 0.5 / D22;
    *spitchsens += Ysens * -2.0 / this->stiffenerPitch;

    // return the final forward analysis output
    return output;

  } else {
    // use the CPT closed-form solution to compute the critical global axial
    // load forward analysis part here
    TacsScalar neg_N11crits[this->NUM_CF_MODES];
    for (int _m1 = 1; _m1 < this->NUM_CF_MODES + 1; _m1++) {
      TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) /
                              this->stiffenerPitch / this->stiffenerPitch;
      TacsScalar m1 = _m1;
      TacsScalar nondim_factor =
          pow(m1 / rho_0, 2.0) + pow(m1 / rho_0, -2.0) + 2.0 * xi;
      neg_N11crits[_m1 - 1] =
          -1.0 * dim_factor * nondim_factor;  // negated only because we have to
                                              // do KS min aggregate later
    }

    // compute KS aggregation sensitivity
    TacsScalar neg_N11crits_sens[this->NUM_CF_MODES];
    TacsScalar neg_N11crit = ksAggregationSens(
        neg_N11crits, this->NUM_CF_MODES, this->ksWeight, neg_N11crits_sens);

    // compute sensitivities here
    for (int _m1 = 1; _m1 < this->NUM_CF_MODES + 1; _m1++) {
      TacsScalar m1 = _m1;
      // backpropagate through the output
      neg_N11crits_sens[_m1 - 1] *= N1sens;
      // apply -1 last step
      neg_N11crits_sens[_m1 - 1] *= -1.0;

      // forward analysis states
      TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) /
                              this->stiffenerPitch / this->stiffenerPitch;
      TacsScalar nondim_factor =
          pow(m1 / rho_0, 2.0) + pow(m1 / rho_0, -2.0) + 2.0 * xi;
      neg_N11crits[_m1 - 1] =
          -1.0 * dim_factor * nondim_factor;  // negated only because we have to
                                              // do KS min aggregate later

      // update sensitivities (left factor is dKS/dv_i, right factor is dv_i /
      // dx)
      *D11sens +=
          neg_N11crits_sens[_m1 - 1] * (0.5 * neg_N11crits[_m1 - 1] / D11);
      *D22sens +=
          neg_N11crits_sens[_m1 - 1] * (0.5 * neg_N11crits[_m1 - 1] / D22);
      *spitchsens += neg_N11crits_sens[_m1 - 1] *
                     (-2.0 * neg_N11crits[_m1 - 1] / this->stiffenerPitch);
      *rho_0sens += neg_N11crits_sens[_m1 - 1] * -dim_factor *
                    (-2.0 * pow(m1 / rho_0, 2.0) / rho_0 +
                     pow(m1 / rho_0, -2.0) * 2.0 / rho_0);
      *xisens += neg_N11crits_sens[_m1 - 1] * -dim_factor * 2.0;
      //*zetasens is unchanged
    }
    return -1.0 * neg_N11crit;
  }
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::computeCriticalGlobalShearLoad(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar b,
    const TacsScalar rho_0, const TacsScalar xi, const TacsScalar gamma,
    const TacsScalar zeta) {
  if (this->getShearGP()) {
    // use Gaussian processes to compute the critical global shear load
    TacsScalar dim_factor =
        M_PI * M_PI * pow(D11 * D22 * D22 * D22, 0.25) / b / b;
    TacsScalar one = 1.0;
    XtestShear[0] = log(one + xi);
    XtestShear[1] = log(rho_0);
    XtestShear[2] = log(
        one + gamma);  // log(1+gamma) = 0 since gamma=0 for unstiffened panel
    XtestShear[3] = log(one + 1000.0 * zeta);
    TacsScalar nd_factor =
        exp(this->panelGPs->predictMeanTestData(2, XtestShear));
    TacsScalar output = dim_factor * nd_factor;
    return dim_factor * nd_factor;

  } else {
    // use the CPT closed-form solution to compute the critical global axial
    // load no mode switching in this solution.. (some error at low aspect
    // ratios
    TacsScalar lam1, lam2;  // lam1bar, lam2bar values
    nondimShearParams(xi, gamma, &lam1, &lam2);
    TacsScalar dim_factor =
        M_PI * M_PI * pow(D11 * D22 * D22 * D22, 0.25) / b / b;
    TacsScalar nondim_factor =
        (1.0 + pow(lam1, 4.0) + 6.0 * pow(lam1 * lam2, 2.0) + pow(lam2, 4.0) +
         2.0 * xi * (lam1 * lam1 + lam2 * lam2) + gamma) /
        (2.0 * lam1 * lam1 * lam2);
    // accounts for high and low ARs here
    // smooth max of (1,rho_0^{-2})
    TacsScalar shear_geom[2];
    shear_geom[0] = 1.0;
    shear_geom[1] = 1.0 / rho_0 / rho_0;
    TacsScalar shear_geom_ks = ksAggregation(shear_geom, 2, this->ksWeight);
    nondim_factor *= shear_geom_ks;
    return dim_factor *
           nondim_factor;  // aka N12_crit from CPT closed-form solution
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::nondimCriticalGlobalShearLoad(
    const TacsScalar rho_0, const TacsScalar xi, const TacsScalar gamma,
    const TacsScalar zeta) {
  if (this->getShearGP()) {
    // use Gaussian processes to compute the critical global shear load
    TacsScalar dim_factor = 1.0;
    TacsScalar one = 1.0;
    XtestShear[0] = log(one + xi);
    XtestShear[1] = log(rho_0);
    XtestShear[2] = log(
        one + gamma);  // log(1+gamma) = 0 since gamma=0 for unstiffened panel
    XtestShear[3] = log(one + 1000.0 * zeta);
    // don't need to use the saved data here since this routine is only meant to
    // be called on a single constitutive object (not a whole mesh with O(1e4)
    // const objects)
    TacsScalar nd_factor =
        exp(this->getShearGP()->predictMeanTestData(XtestShear));
    return dim_factor * nd_factor;

  } else {
    // use the CPT closed-form solution to compute the critical global axial
    // load no mode switching in this solution..
    // CPT closed-form solution now accurate for low and high ARs with some area
    // at intermediate..
    TacsScalar lam1, lam2;  // lam1bar, lam2bar values
    nondimShearParams(xi, gamma, &lam1, &lam2);
    TacsScalar dim_factor = 1.0;
    TacsScalar nondim_factor =
        (1.0 + pow(lam1, 4.0) + 6.0 * pow(lam1 * lam2, 2.0) + pow(lam2, 4.0) +
         2.0 * xi * (lam1 * lam1 + lam2 * lam2) + gamma) /
        (2.0 * lam1 * lam1 * lam2);

    // accounts for high and low ARs here
    // smooth max of (1,rho_0^{-2})
    TacsScalar shear_geom[2];
    shear_geom[0] = 1.0;
    shear_geom[1] = 1.0 / rho_0 / rho_0;
    TacsScalar shear_geom_ks = ksAggregation(shear_geom, 2, this->ksWeight);
    nondim_factor *= shear_geom_ks;

    return dim_factor *
           nondim_factor;  // aka N12_crit from CPT closed-form solution
  }
}

void TACSGPBladeStiffenedShellConstitutive::nondimShearParams(
    const TacsScalar xi, const TacsScalar gamma, TacsScalar* lam1bar,
    TacsScalar* lam2bar) {
  // need to iterate over lam2 with the Newton's method
  TacsScalar lam2bar_sq = 0.0;  // starting guess for lambda2_bar^2

  TacsScalar init_resid = lam2Constraint(lam2bar_sq, xi, gamma);
  // Newton iteration for lam2bar squared
  int ct = 0;
  // use rel tolerance of 1e-13
  while (abs(TacsRealPart(lam2Constraint(lam2bar_sq, xi, gamma)) / init_resid) >
             1e-13 &&
         ct < 50) {
    lam2bar_sq -= lam2Constraint(lam2bar_sq, xi, gamma) /
                  lam2ConstraintDeriv(lam2bar_sq, xi, gamma);
    ct += 1;
  }

  // now compute lam1_bar, lam2_bar
  *lam1bar =
      pow(1.0 + 2.0 * lam2bar_sq * xi + lam2bar_sq * lam2bar_sq + gamma, 0.25);
  *lam2bar = pow(lam2bar_sq, 0.5);
}
TacsScalar TACSGPBladeStiffenedShellConstitutive::lam2Constraint(
    const TacsScalar lam2sq, const TacsScalar xi, const TacsScalar gamma) {
  // compute the residual of the combined lam1bar, lam2bar constraint but on
  // lam2bar^2
  // if x = lam2sq, and lam1sq = f(x), lam1^4 = f(x)^2
  // then we have:

  // f(x) = lam1^2(x) where x = lam2^2
  TacsScalar f = pow(1.0 + 2.0 * lam2sq * xi + lam2sq * lam2sq + gamma, 0.5);
  TacsScalar term3 =
      sqrt((3.0 + xi) / 9.0 + 4.0 / 3.0 * f * xi + 4.0 / 3.0 * f * f);

  return lam2sq + f + xi / 3.0 - term3;
}
TacsScalar TACSGPBladeStiffenedShellConstitutive::lam2ConstraintDeriv(
    const TacsScalar lam2sq, const TacsScalar xi, const TacsScalar gamma) {
  // compute the residual derivatives for the lam2bar constraint above w.r.t.
  // the lam2sq input about a lam2sq input
  // if x = lam2sq, and lam1sq = f(x), lam1^4 = f(x)^2
  // then we have:

  TacsScalar deriv = 1.0;
  // f(x) = lam1^2(x) where x = lam2^2
  TacsScalar f = pow(1.0 + 2.0 * lam2sq * xi + lam2sq * lam2sq + gamma, 0.5);
  // f'(x) = (x + xi) / f(x)
  TacsScalar fp = (lam2sq + xi) / f;
  deriv += fp;

  // term3 the sqrt term with many terms inside
  TacsScalar term3 =
      sqrt((3.0 + xi) / 9.0 + 4.0 / 3.0 * f * xi + 4.0 / 3.0 * f * f);

  // simplified chain rule over sqrt long-term expression (i.e. term3)
  deriv -= 2.0 * fp / 3.0 / term3 * (xi + 2.0 * f);
  return deriv;
}

void TACSGPBladeStiffenedShellConstitutive::nondimShearParamsSens(
    const TacsScalar xi, const TacsScalar gamma, TacsScalar* lam1bar,
    TacsScalar* lam2bar, TacsScalar* dl1xi, TacsScalar* dl1gamma,
    TacsScalar* dl2xi, TacsScalar* dl2gamma) {
  // get the lam1, lam2 from Newton's method
  TacsScalar lam1, lam2;
  nondimShearParams(xi, gamma, &lam1, &lam2);

  // also send out the lam1bar, lam2bar again
  *lam1bar = lam1;
  *lam2bar = lam2;

  // differentiate the nonlinear constraints from nondimShearParam subroutine
  // sys eqns [A,B;C,D] * [lam1bar_dot, lam2bar_dot] = [E,F] for each of the two
  // derivatives

  TacsScalar y1 = 1.0 + 2.0 * lam2 * lam2 * xi + pow(lam2, 4.0) + gamma;
  TacsScalar dy1lam2 = 4.0 * lam2 * xi + 4.0 * lam2 * lam2 * lam2;
  TacsScalar dy1xi = 2.0 * lam2 * lam2;
  TacsScalar dy1gamma = 1.0;
  TacsScalar y2 =
      (3.0 + xi) / 9.0 + 4.0 / 3.0 * (lam1 * lam1 * xi + pow(lam1, 4.0));
  TacsScalar dy2lam1 = 4.0 / 3.0 * (2.0 * lam1 * xi + 4.0 * lam1 * lam1 * lam1);
  TacsScalar dy2xi = 1.0 / 9.0 + 4.0 / 3.0 * lam1 * lam1;
  TacsScalar dy2gamma = 0.0;

  // first for the xi sensitivities
  TacsScalar A1, B1, C1, D1, E1, F1;
  A1 = 1.0;
  B1 = -0.25 * lam1 / y1 * dy1lam2;
  E1 = 0.25 * lam1 / y1 * dy1xi;
  C1 = 2.0 * lam1 - 0.5 * pow(y2, -0.5) * dy2lam1;
  D1 = 2.0 * lam2;
  F1 = -1.0 / 3.0 + 0.5 * pow(y2, -0.5) * dy2xi;
  *dl1xi = (D1 * E1 - B1 * F1) / (A1 * D1 - B1 * C1);
  *dl2xi = (A1 * F1 - C1 * E1) / (A1 * D1 - B1 * C1);

  // then for the gamma sensitivities
  TacsScalar A2, B2, C2, D2, E2, F2;
  A2 = A1;
  B2 = B1;
  E2 = 0.25 * lam1 / y1 * dy1gamma;
  C2 = C1;
  D2 = D1;
  F2 = 0.5 * pow(y2, -0.5) * dy2gamma;
  *dl1gamma = (D2 * E2 - B2 * F2) / (A2 * D2 - B2 * C2);
  *dl2gamma = (A2 * F2 - C2 * E2) / (A2 * D2 - B2 * C2);
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::computeCriticalGlobalShearLoadSens(
    const TacsScalar N12sens, const TacsScalar D11, const TacsScalar D22,
    const TacsScalar b, const TacsScalar rho_0, const TacsScalar xi,
    const TacsScalar gamma, const TacsScalar zeta, TacsScalar* D11sens,
    TacsScalar* D22sens, TacsScalar* bsens, TacsScalar* rho_0sens,
    TacsScalar* xisens, TacsScalar* gammasens, TacsScalar* zetasens) {
  if (this->getShearGP()) {
    // use Gaussian processes to compute the critical global shear load
    TacsScalar dim_factor =
        M_PI * M_PI * pow(D11 * D22 * D22 * D22, 0.25) / b / b;
    TacsScalar one = 1.0;
    XtestShear[0] = log(one + xi);
    XtestShear[1] = log(rho_0);
    XtestShear[2] = log(
        one + gamma);  // log(1+gamma) = 0 since gamma=0 for unstiffened panel
    XtestShear[3] = log(one + 1000.0 * zeta);
    TacsScalar arg = this->panelGPs->predictMeanTestData(2, XtestShear);
    TacsScalar nondim_factor = exp(arg);
    TacsScalar output = dim_factor * nondim_factor;

    // backwards propagate sensitivities out of the axialGP model to nondim
    // params, (this part differentiates the nondim_factor)
    TacsScalar Ysens = N12sens * output;
    this->panelGPs->predictMeanTestDataSens(2, Ysens, XtestShear,
                                            XtestShearSens);

    *xisens += XtestShearSens[0] / (one + xi);
    *rho_0sens +=
        XtestShearSens[1] / rho_0;  // chain rule dlog(rho_0) / drho_0 = 1/rho_0
    *gammasens += XtestShearSens[2] / (1.0 + gamma);
    *zetasens += XtestShearSens[3] / (one + 1000.0 * zeta) * 1000.0;

    // backpropagate the dimensional factor terms out to the material and
    // geometric DVs (this part differentiates the dim_factor)
    *D11sens += Ysens * 0.25 / D11;
    *D22sens += Ysens * 0.75 / D22;
    *bsens += Ysens * -2.0 / b;

    // return the final forward analysis output
    return output;

  } else {
    // use the CPT closed-form solution to compute the critical global axial
    // load no mode switching in this solution.. (only accurate for higher
    // aspect ratios => hence the need for machine learning for the actual
    // solution)
    TacsScalar lam1, lam2;  // lam1bar, lam2bar values
    TacsScalar dl1xi, dl2xi, dl1gamma, dl2gamma;

    // compute the derivatives of the nondimensional constraints
    nondimShearParamsSens(xi, gamma, &lam1, &lam2, &dl1xi, &dl1gamma, &dl2xi,
                          &dl2gamma);

    // compute forward analysis states involved in the N12crit load
    TacsScalar dim_factor =
        M_PI * M_PI * pow(D11 * D22 * D22 * D22, 0.25) / b / b;
    TacsScalar num =
        (1.0 + pow(lam1, 4.0) + 6.0 * pow(lam1 * lam2, 2.0) + pow(lam2, 4.0) +
         2.0 * xi * (lam1 * lam1 + lam2 * lam2) + gamma);
    TacsScalar den = 2.0 * lam1 * lam1 * lam2;
    TacsScalar nondim_factor = num / den;

    // accounts for high and low ARs here
    // smooth max of (1,rho_0^{-2})
    TacsScalar shear_geom[2];
    shear_geom[0] = 1.0;
    shear_geom[1] = 1.0 / rho_0 / rho_0;
    TacsScalar shear_geom_ks = ksAggregation(shear_geom, 2, this->ksWeight);

    // final forward output
    TacsScalar N12crit = dim_factor * nondim_factor * shear_geom_ks;

    // sensitivities for the non_dim factor

    TacsScalar dNDlam1 =
        (4.0 * pow(lam1, 3.0) + 12.0 * lam1 * lam2 * lam2 + 4.0 * lam1 * xi) /
            den -
        num * 4.0 * lam1 * lam2 / den / den;
    TacsScalar dNDlam2 =
        (4.0 * pow(lam2, 3.0) + 12.0 * lam2 * lam1 * lam1 + 4.0 * lam2 * xi) /
            den -
        num * 2.0 * lam1 * lam1 / den / den;

    // compute KS aggregation sensitivity
    TacsScalar shear_geom_sens[2];
    ksAggregationSens(shear_geom, 2, this->ksWeight, shear_geom_sens);
    TacsScalar c_rho0_sens = shear_geom_sens[1] * shear_geom[1] * -2.0 / rho_0;

    // compute the overall sensitivities
    *D11sens += N12sens * N12crit * 0.25 / D11;
    *D22sens += N12sens * N12crit * 0.75 / D22;
    *bsens += N12sens * N12crit * -2.0 / b;
    *xisens += N12sens * dim_factor *
               (dNDlam1 * dl1xi + dNDlam2 * dl2xi +
                2.0 * (lam1 * lam1 + lam2 * lam2) / den) *
               shear_geom_ks;
    *gammasens += N12sens * dim_factor *
                  (dNDlam1 * dl1gamma + dNDlam2 * dl2gamma + 1.0 / den) *
                  shear_geom_ks;
    *rho_0sens += N12sens * c_rho0_sens * dim_factor * nondim_factor;
    // *zetasens are unchanged in closed-form

    // return N12crit from closed-form solution
    return N12crit;
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeCriticalLocalShearLoad(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar rho_0,
    const TacsScalar xi, const TacsScalar zeta) {
  if (this->getShearGP()) {
    // use Gaussian processes to compute the critical global shear load
    TacsScalar s_p = this->stiffenerPitch;
    TacsScalar dim_factor =
        M_PI * M_PI * pow(D11 * D22 * D22 * D22, 0.25) / s_p / s_p;
    TacsScalar one = 1.0;
    XtestShear[0] = log(one + xi);
    XtestShear[1] = log(rho_0);
    XtestShear[2] =
        0.0;  // log(1+gamma) = 0 since gamma=0 for unstiffened panel
    XtestShear[3] = log(one + 1000.0 * zeta);
    TacsScalar nd_factor =
        exp(this->panelGPs->predictMeanTestData(3, XtestShear));
    return dim_factor * nd_factor;

  } else {
    // use the CPT closed-form solution to compute the critical global axial
    // load no mode switching in this solution.. (only accurate for higher
    // aspect ratios => hence the need for machine learning for the actual
    // solution)
    TacsScalar lam1, lam2;  // lam1bar, lam2bar values
    nondimShearParams(xi, 0.0, &lam1, &lam2);
    TacsScalar dim_factor = M_PI * M_PI * pow(D11 * D22 * D22 * D22, 0.25) /
                            this->stiffenerPitch / this->stiffenerPitch;
    TacsScalar nondim_factor =
        (1.0 + pow(lam1, 4.0) + 6.0 * pow(lam1 * lam2, 2.0) + pow(lam2, 4.0) +
         2.0 * xi * (lam1 * lam1 + lam2 * lam2)) /
        (2.0 * lam1 * lam1 * lam2);

    // accounts for high and low ARs here
    // smooth max of (1,rho_0^{-2})
    TacsScalar shear_geom[2];
    shear_geom[0] = 1.0;
    shear_geom[1] = 1.0 / rho_0 / rho_0;
    TacsScalar shear_geom_ks = ksAggregation(shear_geom, 2, this->ksWeight);

    return dim_factor * nondim_factor *
           shear_geom_ks;  // aka N12_crit from CPT closed-form solution
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::nondimCriticalLocalShearLoad(
    const TacsScalar rho_0, const TacsScalar xi, const TacsScalar zeta) {
  if (this->getShearGP()) {
    // use Gaussian processes to compute the critical global shear load
    TacsScalar s_p = this->stiffenerPitch;
    TacsScalar dim_factor = 1.0;
    TacsScalar one = 1.0;
    XtestShear[0] = log(one + xi);
    XtestShear[1] = log(rho_0);
    XtestShear[2] =
        0.0;  // log(1+gamma) = 0 since gamma=0 for unstiffened panel
    XtestShear[3] = log(one + 1000.0 * zeta);
    // don't need to use the saved data here since this routine is only meant to
    // be called on a single constitutive object (not a whole mesh with O(1e4)
    // const objects)
    TacsScalar nd_factor =
        exp(this->getShearGP()->predictMeanTestData(XtestShear));
    return dim_factor * nd_factor;

  } else {
    // use the CPT closed-form solution to compute the critical global axial
    // load no mode switching in this solution.. (only accurate for higher
    // aspect ratios => hence the need for machine learning for the actual
    // solution)
    TacsScalar lam1, lam2;  // lam1bar, lam2bar values
    nondimShearParams(xi, 0.0, &lam1, &lam2);
    TacsScalar dim_factor = 1.0;
    TacsScalar nondim_factor =
        (1.0 + pow(lam1, 4.0) + 6.0 * pow(lam1 * lam2, 2.0) + pow(lam2, 4.0) +
         2.0 * xi * (lam1 * lam1 + lam2 * lam2)) /
        (2.0 * lam1 * lam1 * lam2);

    // accounts for high and low ARs here
    // smooth max of (1,rho_0^{-2})
    TacsScalar shear_geom[2];
    shear_geom[0] = 1.0;
    shear_geom[1] = 1.0 / rho_0 / rho_0;
    TacsScalar shear_geom_ks = ksAggregation(shear_geom, 2, this->ksWeight);
    nondim_factor *= shear_geom_ks;

    return dim_factor *
           nondim_factor;  // aka N12_crit from CPT closed-form solution
  }
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::computeCriticalLocalShearLoadSens(
    const TacsScalar N12sens, const TacsScalar D11, const TacsScalar D22,
    const TacsScalar rho_0, const TacsScalar xi, const TacsScalar zeta,
    TacsScalar* D11sens, TacsScalar* D22sens, TacsScalar* spitchsens,
    TacsScalar* rho_0sens, TacsScalar* xisens, TacsScalar* zetasens) {
  if (this->getShearGP()) {
    // use Gaussian processes to compute the critical global shear load
    TacsScalar s_p = this->stiffenerPitch;
    TacsScalar dim_factor =
        M_PI * M_PI * pow(D11 * D22 * D22 * D22, 0.25) / s_p / s_p;
    TacsScalar one = 1.0;
    XtestShear[0] = log(one + xi);
    XtestShear[1] = log(rho_0);
    XtestShear[2] =
        0.0;  // log(1+gamma) = 0 since gamma=0 for unstiffened panel
    XtestShear[3] = log(one + 1000.0 * zeta);
    TacsScalar arg = this->panelGPs->predictMeanTestData(3, XtestShear);
    TacsScalar nondim_factor = exp(arg);
    TacsScalar output = dim_factor * nondim_factor;

    // backwards propagate sensitivities out of the axialGP model to nondim
    // params, (this part differentiates the nondim_factor)
    TacsScalar Ysens = N12sens * output;
    this->panelGPs->predictMeanTestDataSens(3, Ysens, XtestShear,
                                            XtestShearSens);

    *xisens += XtestShearSens[0] / (one + xi);
    *rho_0sens +=
        XtestShearSens[1] / rho_0;  // chain rule dlog(rho_0) / drho_0 = 1/rho_0
    *zetasens += XtestShearSens[3] / (one + 1000.0 * zeta) * 1000.0;

    // backpropagate the dimensional factor terms out to the material and
    // geometric DVs, (this part is differentiating dim_factor)
    *D11sens += Ysens * 0.25 / D11;
    *D22sens += Ysens * 0.75 / D22;
    *spitchsens += Ysens * -2.0 / s_p;

    // return the final forward analysis output
    return output;

  } else {
    // use the CPT closed-form solution to compute the critical global axial
    // load no mode switching in this solution.. (only accurate for higher
    // aspect ratios => hence the need for machine learning for the actual
    // solution)
    TacsScalar lam1, lam2;  // lam1bar, lam2bar values
    TacsScalar dl1xi, dl2xi, _dl1gamma,
        _dl2gamma;  // gamma derivatives are private since unused here (gamma=0
                    // input)

    // compute the derivatives of the nondimensional constraints
    nondimShearParamsSens(xi, 0.0, &lam1, &lam2, &dl1xi, &_dl1gamma, &dl2xi,
                          &_dl2gamma);

    // compute forward analysis states involved in the N12crit load
    TacsScalar dim_factor = M_PI * M_PI * pow(D11 * D22 * D22 * D22, 0.25) /
                            this->stiffenerPitch / this->stiffenerPitch;
    TacsScalar num = (1.0 + pow(lam1, 4.0) + 6.0 * pow(lam1 * lam2, 2.0) +
                      pow(lam2, 4.0) + 2.0 * xi * (lam1 * lam1 + lam2 * lam2));
    TacsScalar den = 2.0 * lam1 * lam1 * lam2;
    TacsScalar nondim_factor = num / den;

    // accounts for high and low ARs here
    // smooth max of (1,rho_0^{-2})
    TacsScalar shear_geom[2];
    shear_geom[0] = 1.0;
    shear_geom[1] = 1.0 / rho_0 / rho_0;
    TacsScalar shear_geom_ks = ksAggregation(shear_geom, 2, this->ksWeight);

    // final forward output
    TacsScalar N12crit = dim_factor * nondim_factor * shear_geom_ks;

    // sensitivities for the non_dim factor
    TacsScalar dNDlam1 =
        (4.0 * pow(lam1, 3.0) + 12.0 * lam1 * lam2 * lam2 + 4.0 * lam1 * xi) /
            den -
        num * 4.0 * lam1 * lam2 / den / den;
    TacsScalar dNDlam2 =
        (4.0 * pow(lam2, 3.0) + 12.0 * lam2 * lam1 * lam1 + 4.0 * lam2 * xi) /
            den -
        num * 2.0 * lam1 * lam1 / den / den;

    // compute KS aggregation sensitivity
    TacsScalar shear_geom_sens[2];
    ksAggregationSens(shear_geom, 2, this->ksWeight, shear_geom_sens);
    TacsScalar c_rho0_sens = shear_geom_sens[1] * shear_geom[1] * -2.0 / rho_0;

    // compute the overall sensitivities
    *D11sens += N12sens * N12crit * 0.25 / D11;
    *D22sens += N12sens * N12crit * 0.75 / D22;
    *spitchsens += N12sens * N12crit * -2.0 / this->stiffenerPitch;
    *xisens += N12sens * dim_factor *
               (dNDlam1 * dl1xi + dNDlam2 * dl2xi +
                2.0 * (lam1 * lam1 + lam2 * lam2) / den) *
               shear_geom_ks;
    *rho_0sens += N12sens * c_rho0_sens * dim_factor * nondim_factor;
    // zetasens unchanged in closed-form

    // return N12crit from closed-form solution
    return N12crit;
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeStiffenerCripplingLoad(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar xi,
    const TacsScalar rho_0, const TacsScalar genPoiss, const TacsScalar zeta) {
  if (this->getCripplingGP()) {
    // use Gaussian processes to compute the critical global axial load
    TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) /
                            this->stiffenerHeight / this->stiffenerHeight;
    TacsScalar one = 1.0;
    XtestCrippling[0] = log(one + xi);
    XtestCrippling[1] = log(rho_0);
    XtestCrippling[2] = log(genPoiss);
    XtestCrippling[3] = log(one + 1000.0 * zeta);
    TacsScalar nd_factor =
        exp(this->panelGPs->predictMeanTestData(4, XtestCrippling));
    return dim_factor * nd_factor;

  } else {
    // use the literature CPT closed-form solution for approximate stiffener
    // crippling, not a function of aspect ratio
    TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) /
                            this->stiffenerHeight / this->stiffenerHeight;
    TacsScalar nondim_factor = (0.476 - 0.56 * (genPoiss - 0.2)) * xi;
    return dim_factor * nondim_factor;
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::nondimStiffenerCripplingLoad(
    const TacsScalar rho_0, const TacsScalar xi, const TacsScalar genPoiss,
    const TacsScalar zeta) {
  if (this->getCripplingGP()) {
    // use Gaussian processes to compute the critical global axial load
    TacsScalar dim_factor = 1.0;
    TacsScalar one = 1.0;
    XtestCrippling[0] = log(one + xi);
    XtestCrippling[1] = log(rho_0);
    XtestCrippling[2] = log(genPoiss);
    XtestCrippling[3] = log(one + 1000.0 * zeta);
    // don't need to use the saved data here since this routine is only meant to
    // be called on a single constitutive object (not a whole mesh with O(1e4)
    // const objects)
    TacsScalar nd_factor =
        exp(this->getCripplingGP()->predictMeanTestData(XtestCrippling));
    return dim_factor * nd_factor;

  } else {
    // use the literature CPT closed-form solution for approximate stiffener
    // crippling, not a function of aspect ratio
    TacsScalar dim_factor = 1.0;
    TacsScalar nondim_factor = (0.476 - 0.56 * (genPoiss - 0.2)) * xi;
    return dim_factor * nondim_factor;
  }
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::computeStiffenerCripplingLoadSens(
    const TacsScalar N1sens, const TacsScalar D11, const TacsScalar D22,
    const TacsScalar xi, const TacsScalar rho_0, const TacsScalar genPoiss,
    const TacsScalar zeta, TacsScalar* D11sens, TacsScalar* D22sens,
    TacsScalar* sheightsens, TacsScalar* xisens, TacsScalar* rho_0sens,
    TacsScalar* genPoiss_sens, TacsScalar* zetasens) {
  if (this->getCripplingGP()) {
    // use Gaussian processes to compute the critical global axial load
    TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) /
                            this->stiffenerHeight / this->stiffenerHeight;
    TacsScalar one = 1.0;
    XtestCrippling[0] = log(one + xi);
    XtestCrippling[1] = log(rho_0);
    XtestCrippling[2] = log(genPoiss);
    XtestCrippling[3] = log(one + 1000.0 * zeta);
    TacsScalar arg = this->panelGPs->predictMeanTestData(4, XtestCrippling);
    TacsScalar nondim_factor = exp(arg);
    TacsScalar output = dim_factor * nondim_factor;

    // backwards propagate sensitivities out of the crippling model to nondim
    // params, (this part differentiates the nondim_factor)
    TacsScalar Ysens = N1sens * output;
    this->panelGPs->predictMeanTestDataSens(4, Ysens, XtestCrippling,
                                            XtestCripplingSens);

    *xisens += XtestCripplingSens[0] / (one + xi);
    *rho_0sens += XtestCripplingSens[1] /
                  rho_0;  // chain rule dlog(rho_0) / drho_0 = 1/rho_0
    *genPoiss_sens += XtestCripplingSens[2] / genPoiss;
    *zetasens += XtestCripplingSens[3] / (one + 1000.0 * zeta) * 1000.0;

    // backpropagate the dimensional factor terms out to the material and
    // geometric DVs, (this part is differentiating the dimensional factor)
    *D11sens += Ysens * 0.5 / D11;
    *D22sens += Ysens * 0.5 / D22;
    *sheightsens += Ysens * -2.0 / this->stiffenerHeight;

    // return the final forward analysis output
    return output;

  } else {
    // use the literature CPT closed-form solution for approximate stiffener
    // crippling, not a function of aspect ratio
    TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) /
                            this->stiffenerHeight / this->stiffenerHeight;
    TacsScalar nondim_factor = (0.476 - 0.56 * (genPoiss - 0.2)) * xi;
    TacsScalar N11crit = dim_factor * nondim_factor;

    // compute the derivatives
    TacsScalar outputsens = N1sens * N11crit;
    *D11sens += outputsens * 0.5 / D11;
    *D22sens += outputsens * 0.5 / D22;
    *xisens += outputsens / xi;
    *genPoiss_sens += outputsens / nondim_factor * -0.56 * xi;
    *sheightsens += outputsens * -2.0 / this->stiffenerHeight;
    // *rho_0sens, *zetasens unchanged here

    // output the critical load here
    return N11crit;
  }
}

// -----------------------------------------------------------
//               DERIVATIVE TESTING SUBROUTINES
// -----------------------------------------------------------
// -----------------------------------------------------------

TacsScalar TACSGPBladeStiffenedShellConstitutive::testAffineAspectRatio(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 4;
  TacsScalar p_input[n_input], x0[n_input], x[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = 10.341;  // D11
  x0[1] = 5.216;   // D22
  x0[2] = 3.124;   // a
  x0[3] = 1.061;   // b

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  f0 = computeAffineAspectRatio(x[0], x[1], x[2], x[3]);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  f2 = computeAffineAspectRatio(x[0], x[1], x[2], x[3]);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  computeAffineAspectRatioSens(p_output, x0[0], x0[1], x0[2], x0[3],
                               &input_sens[0], &input_sens[1], &input_sens[2],
                               &input_sens[3]);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testAffineAspectRatio:\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testLaminateIsotropy(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 4;
  TacsScalar p_input[n_input], x0[n_input], x[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = 10.341;  // D11
  x0[1] = 5.216;   // D22
  x0[2] = 6.132;   // D12
  x0[3] = 2.103;   // D66

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  f0 = computeLaminateIsotropy(x[0], x[1], x[2], x[3]);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  f2 = computeLaminateIsotropy(x[0], x[1], x[2], x[3]);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  computeLaminateIsotropySens(p_output, x0[0], x0[1], x0[2], x0[3],
                              &input_sens[0], &input_sens[1], &input_sens[2],
                              &input_sens[3]);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testLaminateIsotropy:\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testGeneralizedPoissonsRatio(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 2;
  TacsScalar p_input[n_input], x0[n_input], x[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = 10.341;  // D12
  x0[1] = 5.381;   // D66

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  f0 = computeGeneralizedPoissonsRatio(x[0], x[1]);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  f2 = computeGeneralizedPoissonsRatio(x[0], x[1]);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  computeGeneralizedPoissonsRatioSens(p_output, x0[0], x0[1], &input_sens[0],
                                      &input_sens[1]);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testGeneralizedPoissonsRatio:\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testStiffenerAreaRatio(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 6;
  TacsScalar p_input[n_input], x0[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = this->stiffenerThick;
  x0[1] = this->stiffenerHeight;
  x0[2] = this->stiffenerPitch;
  x0[3] = this->panelThick;
  x0[4] = 1.145e8;  // E11p
  x0[5] = 2.17e7;   // E11s

  TacsScalar E1p, E1s;
  E1p = x0[4] * 1.0;
  E1s = x0[5] * 1.0;

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;
  this->stiffenerThick -= p_input[0] * epsilon;
  this->stiffenerHeight -= p_input[1] * epsilon;
  this->stiffenerPitch -= p_input[2] * epsilon;
  this->panelThick -= p_input[3] * epsilon;
  E1p -= p_input[4] * epsilon;
  E1s -= p_input[5] * epsilon;
  f0 = computeStiffenerAreaRatio(E1p, E1s);

  this->stiffenerThick += 2.0 * p_input[0] * epsilon;
  this->stiffenerHeight += 2.0 * p_input[1] * epsilon;
  this->stiffenerPitch += 2.0 * p_input[2] * epsilon;
  this->panelThick += 2.0 * p_input[3] * epsilon;
  E1p += 2.0 * p_input[4] * epsilon;
  E1s += 2.0 * p_input[5] * epsilon;
  f2 = computeStiffenerAreaRatio(E1p, E1s);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // reset the values
  this->stiffenerThick -= p_input[0] * epsilon;
  this->stiffenerHeight -= p_input[1] * epsilon;
  this->stiffenerPitch -= p_input[2] * epsilon;
  this->panelThick -= p_input[3] * epsilon;
  E1p -= p_input[4] * epsilon;
  E1s -= p_input[5] * epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  computeStiffenerAreaRatioSens(p_output, E1p, E1s, &input_sens[0],
                                &input_sens[1], &input_sens[2], &input_sens[3],
                                &input_sens[4], &input_sens[5]);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testStiffenerAreaRatio:\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testStiffenerStiffnessRatio(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 6;
  TacsScalar p_input[n_input], x0[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = 10.2143;  // D11
  x0[1] = this->stiffenerThick;
  x0[2] = this->stiffenerHeight;
  x0[3] = this->stiffenerPitch;
  x0[4] = 1.145e8;                        // E1s
  x0[5] = this->stiffenerHeight * 0.143;  // zn overall centroid

  TacsScalar E1s = x0[4] * 1.0;
  TacsScalar zn = x0[5] * 1.0;

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;
  TacsScalar D11 = x0[0] * 1.0;
  D11 -= p_input[0] * epsilon;
  this->stiffenerThick -= p_input[1] * epsilon;
  this->stiffenerHeight -= p_input[2] * epsilon;
  this->stiffenerPitch -= p_input[3] * epsilon;
  E1s -= p_input[4] * epsilon;
  zn -= p_input[5] * epsilon;
  f0 = computeStiffenerStiffnessRatio(D11, E1s, zn);

  D11 += 2.0 * p_input[0] * epsilon;
  this->stiffenerThick += 2.0 * p_input[1] * epsilon;
  this->stiffenerHeight += 2.0 * p_input[2] * epsilon;
  this->stiffenerPitch += 2.0 * p_input[3] * epsilon;
  E1s += 2.0 * p_input[4] * epsilon;
  zn += 2.0 * p_input[5] * epsilon;
  f2 = computeStiffenerStiffnessRatio(D11, E1s, zn);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // reset the values
  D11 -= p_input[0] * epsilon;
  this->stiffenerThick -= p_input[1] * epsilon;
  this->stiffenerHeight -= p_input[2] * epsilon;
  this->stiffenerPitch -= p_input[3] * epsilon;
  E1s -= p_input[4] * epsilon;
  zn -= p_input[5] * epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  computeStiffenerStiffnessRatioSens(
      p_output, D11, E1s, zn, &input_sens[0], &input_sens[1], &input_sens[2],
      &input_sens[3], &input_sens[4], &input_sens[5]);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testStiffenerStiffnessRatio:\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testTransverseShearParameter(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 4;
  TacsScalar p_input[n_input], x0[n_input], x[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = 100.234;  // A66
  x0[1] = 421.341;  // A11
  x0[2] = 2.134;    // b
  x0[3] = 0.0112;   // h

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  f0 = computeTransverseShearParameter(x[0], x[1], x[2], x[3]);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  f2 = computeTransverseShearParameter(x[0], x[1], x[2], x[3]);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  computeTransverseShearParameterSens(p_output, x0[0], x0[1], x0[2], x0[3],
                                      &input_sens[0], &input_sens[1],
                                      &input_sens[2], &input_sens[3]);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testTransverseShearParameter:\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testNondimensionalParameters(
    TacsScalar epsilon, int printLevel) {
  // run each of the nondim parameter tests and aggregate the max among them
  const int n_tests = 6;
  TacsScalar relErrors[n_tests];

  if (printLevel != 0) {
    printf("\nTACSGPBladeStiffened..testNondimensionalParmeters start::\n");
    printf("--------------------------------------------------------\n\n");
  }

  relErrors[0] = testAffineAspectRatio(epsilon, printLevel);
  relErrors[1] = testLaminateIsotropy(epsilon, printLevel);
  relErrors[2] = testGeneralizedPoissonsRatio(epsilon, printLevel);
  relErrors[3] = testStiffenerAreaRatio(epsilon, printLevel);
  relErrors[4] = testStiffenerStiffnessRatio(epsilon, printLevel);
  relErrors[5] = testTransverseShearParameter(epsilon, printLevel);

  // get max rel error among them
  TacsScalar maxRelError = 0.0;
  for (int i = 0; i < n_tests; i++) {
    if (TacsRealPart(relErrors[i]) > TacsRealPart(maxRelError)) {
      maxRelError = relErrors[i];
    }
  }

  // report the overall test results
  if (printLevel != 0) {
    printf(
        "\nTACSGPBladeStiffened..testNondimensionalParmeters final "
        "results::\n");
    printf("\ttestAffineAspectRatio = %.4e\n", TacsRealPart(relErrors[0]));
    printf("\ttestLaminateIsotropy = %.4e\n", TacsRealPart(relErrors[1]));
    printf("\ttestGeneralizedPoissonsRatio = %.4e\n",
           TacsRealPart(relErrors[2]));
    printf("\ttestStiffenerAreaRatio = %.4e\n", TacsRealPart(relErrors[3]));
    printf("\ttestStiffenerStiffnessRatio = %.4e\n",
           TacsRealPart(relErrors[4]));
    printf("\ttestTransverseShearParameter = %.4e\n",
           TacsRealPart(relErrors[5]));
    printf("\tOverall max rel error = %.4e\n\n", TacsRealPart(maxRelError));
  }

  return maxRelError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testAxialCriticalLoads(
    TacsScalar epsilon, int printLevel) {
  // run each of the nondim parameter tests and aggregate the max among them
  const int n_tests = 2;
  TacsScalar relErrors[n_tests];

  if (printLevel != 0) {
    printf("\nTACSGPBladeStiffened..testAxialCriticalLoads start::\n");
    printf("--------------------------------------------------------\n\n");
  }

  relErrors[0] = testCriticalGlobalAxialLoad(epsilon, printLevel);
  relErrors[1] = testCriticalLocalAxialLoad(epsilon, printLevel);

  // get max rel error among them
  TacsScalar maxRelError = 0.0;
  for (int i = 0; i < n_tests; i++) {
    if (TacsRealPart(relErrors[i]) > TacsRealPart(maxRelError)) {
      maxRelError = relErrors[i];
    }
  }

  // get max rel error among them
  if (printLevel != 0) {
    printf("\nTACSGPBladeStiffened..testAxialCriticalLoads final results::\n");
    printf("\ttestGlobalAxialLoad = %.4e\n", TacsRealPart(relErrors[0]));
    printf("\ttestLocalAxialLoad = %.4e\n", TacsRealPart(relErrors[1]));
    printf("\tOverall max rel error = %.4e\n\n", TacsRealPart(maxRelError));
  }

  return maxRelError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testCriticalGlobalAxialLoad(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 8;
  TacsScalar p_input[n_input], x0[n_input], x[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }

  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = 10.2412;  // D11
  x0[1] = 5.4323;   // D22
  x0[2] = 2.134;    // b
  x0[3] = 0.13432;  // delta
  x0[4] = 2.4545;   // rho0
  x0[5] = 1.24332;  // xi
  x0[6] = 0.2454;   // gamma
  x0[7] = 1e-3;     // zeta

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  this->panelGPs->resetSavedData();
  f0 = computeCriticalGlobalAxialLoad(x[0], x[1], x[2], x[3], x[4], x[5], x[6],
                                      x[7]);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  this->panelGPs->resetSavedData();
  f2 = computeCriticalGlobalAxialLoad(x[0], x[1], x[2], x[3], x[4], x[5], x[6],
                                      x[7]);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i];
  }
  this->panelGPs->resetSavedData();
  computeCriticalGlobalAxialLoadSens(
      p_output, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], &input_sens[0],
      &input_sens[1], &input_sens[2], &input_sens[3], &input_sens[4],
      &input_sens[5], &input_sens[6], &input_sens[7]);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testCriticalGlobalAxialLoad:\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testCriticalLocalAxialLoad(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 6;
  TacsScalar p_input[n_input], x0[n_input], x[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = 10.2412;               // D11
  x0[1] = 5.4323;                // D22
  x0[2] = this->stiffenerPitch;  // s_p
  x0[3] = 2.4545;                // rho0
  x0[4] = 1.24332;               // xi
  x0[5] = 1e-3;                  // zeta

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  this->stiffenerPitch = x0[2] - p_input[2] * epsilon;
  this->panelGPs->resetSavedData();
  f0 = computeCriticalLocalAxialLoad(x[0], x[1], x[3], x[4], x[5]);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  this->stiffenerPitch = x0[2] + p_input[2] * epsilon;
  this->panelGPs->resetSavedData();
  f2 = computeCriticalLocalAxialLoad(x[0], x[1], x[3], x[4], x[5]);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  this->stiffenerPitch = x0[2];
  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i];
  }
  this->panelGPs->resetSavedData();
  computeCriticalLocalAxialLoadSens(
      p_output, x[0], x[1], x[3], x[4], x[5], &input_sens[0], &input_sens[1],
      &input_sens[2], &input_sens[3], &input_sens[4], &input_sens[5]);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testCriticalLocalAxialLoad:\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testOverallCentroid(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 6;
  TacsScalar p_input[n_input], x0[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = this->stiffenerThick;
  x0[1] = this->stiffenerHeight;
  x0[2] = this->panelThick;
  x0[3] = this->stiffenerPitch;
  x0[4] = 1.734e5;  // E11p
  x0[5] = 8.213e4;  // E11s
  TacsScalar E1p = x0[4] * 1.0;
  TacsScalar E1s = x0[5] * 1.0;

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;
  this->stiffenerThick -= p_input[0] * epsilon;
  this->stiffenerHeight -= p_input[1] * epsilon;
  this->panelThick -= p_input[2] * epsilon;
  this->stiffenerPitch -= p_input[3] * epsilon;
  E1p -= p_input[4] * epsilon;
  E1s -= p_input[5] * epsilon;
  f0 = computeOverallCentroid(E1p, E1s);

  this->stiffenerThick += 2.0 * p_input[0] * epsilon;
  this->stiffenerHeight += 2.0 * p_input[1] * epsilon;
  this->panelThick += 2.0 * p_input[2] * epsilon;
  this->stiffenerPitch += 2.0 * p_input[3] * epsilon;
  E1p += 2.0 * p_input[4] * epsilon;
  E1s += 2.0 * p_input[5] * epsilon;
  f2 = computeOverallCentroid(E1p, E1s);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // reset the values
  this->stiffenerThick -= p_input[0] * epsilon;
  this->stiffenerHeight -= p_input[1] * epsilon;
  this->panelThick -= p_input[2] * epsilon;
  this->stiffenerPitch -= p_input[3] * epsilon;
  E1p -= p_input[4] * epsilon;
  E1s -= p_input[5] * epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  computeOverallCentroidSens(p_output, E1p, E1s, &input_sens[0], &input_sens[1],
                             &input_sens[2], &input_sens[3], &input_sens[4],
                             &input_sens[5]);

  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testOverallCentroid:\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::testPanelGlobalBucklingStiffness(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 4;
  TacsScalar p_input[n_input], x0[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = this->stiffenerPitch;
  x0[1] = this->panelThick;
  x0[2] = 1.734e2;  // E11p (lower value to make test easier)
  x0[3] = 0.0413;   // zn
  TacsScalar E1p = x0[2] * 1.0;
  TacsScalar zn = x0[3] * 1.0;

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;
  this->stiffenerPitch -= p_input[0] * epsilon;
  this->panelThick -= p_input[1] * epsilon;
  E1p -= p_input[2] * epsilon;
  zn -= p_input[3] * epsilon;
  computePanelGlobalBucklingStiffness(E1p, zn, &f0);

  this->stiffenerPitch += 2.0 * p_input[0] * epsilon;
  this->panelThick += 2.0 * p_input[1] * epsilon;
  E1p += 2.0 * p_input[2] * epsilon;
  zn += 2.0 * p_input[3] * epsilon;
  computePanelGlobalBucklingStiffness(E1p, zn, &f2);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // reset the values
  this->stiffenerPitch -= p_input[0] * epsilon;
  this->panelThick -= p_input[1] * epsilon;
  E1p -= p_input[2] * epsilon;
  zn -= p_input[3] * epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  computePanelGlobalBucklingStiffnessSens(p_output, E1p, zn, &input_sens[0],
                                          &input_sens[1], &input_sens[2],
                                          &input_sens[3]);

  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testPanelGlobalBucklingStiffness:\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testOtherTests(
    TacsScalar epsilon, int printLevel) {
  // run each of the nondim parameter tests and aggregate the max among them
  const int n_tests = 2;
  TacsScalar relErrors[n_tests];

  if (printLevel != 0) {
    printf("\nTACSGPBladeStiffened..testOtherTests start::\n");
    printf("--------------------------------------------------------\n\n");
  }

  // first two tests are used for the closed-form solution approach
  relErrors[0] = testOverallCentroid(epsilon, printLevel);
  relErrors[1] = testPanelGlobalBucklingStiffness(epsilon, printLevel);

  // get max rel error among them
  TacsScalar maxRelError = 0.0;
  for (int i = 0; i < n_tests; i++) {
    if (TacsRealPart(relErrors[i]) > TacsRealPart(maxRelError)) {
      maxRelError = relErrors[i];
    }
  }

  // get max rel error among them
  if (printLevel != 0) {
    printf("\nTACSGPBladeStiffened..testOtherTests final results::\n");
    printf("\ttestOverallCentroid = %.4e\n", TacsRealPart(relErrors[0]));
    printf("\ttestPanelGlobalBucklingStiffness = %.4e\n",
           TacsRealPart(relErrors[1]));
    printf("\tOverall max rel error = %.4e\n\n", TacsRealPart(maxRelError));
  }

  return maxRelError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testShearCriticalLoads(
    TacsScalar epsilon, int printLevel) {
  // run each of the nondim parameter tests and aggregate the max among them
  const int n_tests = 6;
  TacsScalar relErrors[n_tests];

  if (printLevel != 0) {
    printf("\nTACSGPBladeStiffened..testShearCriticalLoads start::\n");
    printf("--------------------------------------------------------\n\n");
  }

  // first two tests are used for the closed-form solution approach
  relErrors[0] = testLam2Constraint(epsilon, printLevel);
  relErrors[1] = testNondimShearParams(epsilon, printLevel);
  // final crit load tests
  relErrors[2] = testCriticalGlobalShearLoad(epsilon, printLevel);
  relErrors[3] = testCriticalLocalShearLoad(epsilon, printLevel);
  relErrors[4] = testCriticalGlobalShearLoad_LowAR(epsilon, printLevel);
  relErrors[5] = testCriticalLocalShearLoad_LowAR(epsilon, printLevel);

  // get max rel error among them
  TacsScalar maxRelError = 0.0;
  for (int i = 0; i < n_tests; i++) {
    if (TacsRealPart(relErrors[i]) > TacsRealPart(maxRelError)) {
      maxRelError = relErrors[i];
    }
  }

  // get max rel error among them
  if (printLevel != 0) {
    printf("\nTACSGPBladeStiffened..testShearCriticalLoads final results::\n");
    printf("\ttestLam2Constraint = %.4e\n", TacsRealPart(relErrors[0]));
    printf("\ttestNondimShearParams = %.4e\n", TacsRealPart(relErrors[1]));
    printf("\ttestGlobalShearLoad = %.4e\n", TacsRealPart(relErrors[2]));
    printf("\ttestLocalShearLoad = %.4e\n", TacsRealPart(relErrors[3]));
    printf("\ttestGlobalShearLoad_LowAR = %.4e\n", TacsRealPart(relErrors[4]));
    printf("\ttestLocalShearLoad_LowAR = %.4e\n", TacsRealPart(relErrors[5]));
    printf("\tOverall max rel error = %.4e\n\n", TacsRealPart(maxRelError));
  }

  return maxRelError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testLam2Constraint(
    const TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 1;
  TacsScalar p_input[n_input], x0[n_input], x[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = 0.78153;  // lam2Sq
  // x0[1] = 1.24332;   // xi
  // x0[2] = 0.2454;    //gamma
  TacsScalar xi = 1.24332;    // xi
  TacsScalar gamma = 0.2454;  // gamma

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  f0 = lam2Constraint(x[0], xi, gamma);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  f2 = lam2Constraint(x[0], xi, gamma);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i];
  }
  TacsScalar lam2Deriv = lam2ConstraintDeriv(x[0], xi, gamma);
  TacsScalar adjTD = p_output * lam2Deriv * p_input[0];
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testLam2Constraint:\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testNondimShearParams(
    const TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 2;
  TacsScalar p_input[n_input], x0[n_input], x[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  // temporarily fix gamma perturbation
  // p_input[1] = 0.0;

  const int n_output = 2;
  TacsScalar p_output[n_output], f0[n_output], f2[n_output];
  for (int ii = 0; ii < n_output; ii++) {
    p_output[ii] = ((double)rand() / (RAND_MAX));
  }

  // compute initial values
  x0[0] = 1.24332;  // xi
  x0[1] = 0.2454;   // gamma

  // perform central difference over rho_0 function on [D11,D22,a,b]
  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  nondimShearParams(x[0], x[1], &f0[0], &f0[1]);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  nondimShearParams(x[0], x[1], &f2[0], &f2[1]);

  TacsScalar centralDiff = 0.0;
  for (int ii = 0; ii < n_output; ii++) {
    centralDiff += p_output[ii] * (f2[ii] - f0[ii]) / 2.0 / epsilon;
  }

  // now perform the adjoint sensitivity
  TacsScalar _lam[n_input], l1sens[n_input], l2sens[n_input];
  memset(l1sens, 0, n_input * sizeof(TacsScalar));
  memset(l2sens, 0, n_input * sizeof(TacsScalar));
  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i];
  }
  nondimShearParamsSens(x[0], x[1], &_lam[0], &_lam[1], &l1sens[0], &l1sens[1],
                        &l2sens[0], &l2sens[1]);
  TacsScalar adjTD = 0.0;
  for (int i = 0; i < n_input; i++) {
    adjTD += p_output[0] * l1sens[i] * p_input[i];
    adjTD += p_output[1] * l2sens[i] * p_input[i];
  }

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testNondimShearParams:\n");
    printf("\t\tadjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\tcentralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\trel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testCriticalGlobalShearLoad(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 7;

  TacsScalar p_input[n_input], x0[n_input], x[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = 10.2412;  // D11
  x0[1] = 5.4323;   // D22
  x0[2] = 2.134;    // b
  x0[3] = 2.4545;   // rho0
  x0[4] = 1.24332;  // xi
  x0[5] = 0.2454;   // gamma
  x0[6] = 1e-3;     // zeta

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  this->panelGPs->resetSavedData();
  f0 = computeCriticalGlobalShearLoad(x[0], x[1], x[2], x[3], x[4], x[5], x[6]);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  this->panelGPs->resetSavedData();
  f2 = computeCriticalGlobalShearLoad(x[0], x[1], x[2], x[3], x[4], x[5], x[6]);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i];
  }
  this->panelGPs->resetSavedData();
  computeCriticalGlobalShearLoadSens(
      p_output, x[0], x[1], x[2], x[3], x[4], x[5], x[6], &input_sens[0],
      &input_sens[1], &input_sens[2], &input_sens[3], &input_sens[4],
      &input_sens[5], &input_sens[6]);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testCriticalGlobalShearLoad:\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::testCriticalGlobalShearLoad_LowAR(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 7;
  TacsScalar p_input[n_input], x0[n_input], x[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = 10.2412;  // D11
  x0[1] = 5.4323;   // D22
  x0[2] = 2.134;    // b
  x0[3] = 0.2543;   // rho0
  x0[4] = 1.24332;  // xi
  x0[5] = 0.2454;   // gamma
  x0[6] = 1e-3;     // zeta

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  this->panelGPs->resetSavedData();
  f0 = computeCriticalGlobalShearLoad(x[0], x[1], x[2], x[3], x[4], x[5], x[6]);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  this->panelGPs->resetSavedData();
  f2 = computeCriticalGlobalShearLoad(x[0], x[1], x[2], x[3], x[4], x[5], x[6]);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i];
  }
  this->panelGPs->resetSavedData();
  computeCriticalGlobalShearLoadSens(
      p_output, x[0], x[1], x[2], x[3], x[4], x[5], x[6], &input_sens[0],
      &input_sens[1], &input_sens[2], &input_sens[3], &input_sens[4],
      &input_sens[5], &input_sens[6]);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testCriticalGlobalShearLoad:_LowAR\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testCriticalLocalShearLoad(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 6;
  TacsScalar p_input[n_input], x0[n_input], x[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = 10.2412;               // D11
  x0[1] = 5.4323;                // D22
  x0[2] = this->stiffenerPitch;  // s_p
  x0[3] = 1.24332;               // xi
  x0[4] = 2.4545;                // rho0
  x0[5] = 1e-3;                  // zeta

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  this->stiffenerPitch = x0[2] - p_input[2] * epsilon;
  this->panelGPs->resetSavedData();
  f0 = computeCriticalLocalShearLoad(x[0], x[1], x[3], x[4], x[5]);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  this->panelGPs->resetSavedData();
  this->stiffenerPitch = x0[2] + p_input[2] * epsilon;
  f2 = computeCriticalLocalShearLoad(x[0], x[1], x[3], x[4], x[5]);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  this->stiffenerPitch = x0[2];
  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i];
  }
  this->panelGPs->resetSavedData();
  computeCriticalLocalShearLoadSens(
      p_output, x[0], x[1], x[3], x[4], x[5], &input_sens[0], &input_sens[1],
      &input_sens[2], &input_sens[3], &input_sens[4], &input_sens[5]);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testCriticalLocalShearLoad:\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar
TACSGPBladeStiffenedShellConstitutive::testCriticalLocalShearLoad_LowAR(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 6;
  TacsScalar p_input[n_input], x0[n_input], x[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  x0[0] = 10.2412;               // D11
  x0[1] = 5.4323;                // D22
  x0[2] = this->stiffenerPitch;  // s_p
  x0[3] = 1.24332;               // xi
  x0[4] = 0.34251;               // rho0
  x0[5] = 1e-3;                  // zeta

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  this->stiffenerPitch = x0[2] - p_input[2] * epsilon;
  this->panelGPs->resetSavedData();
  f0 = computeCriticalLocalShearLoad(x[0], x[1], x[3], x[4], x[5]);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  this->panelGPs->resetSavedData();
  this->stiffenerPitch = x0[2] + p_input[2] * epsilon;
  f2 = computeCriticalLocalShearLoad(x[0], x[1], x[3], x[4], x[5]);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  this->stiffenerPitch = x0[2];
  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i];
  }
  this->panelGPs->resetSavedData();
  computeCriticalLocalShearLoadSens(
      p_output, x[0], x[1], x[3], x[4], x[5], &input_sens[0], &input_sens[1],
      &input_sens[2], &input_sens[3], &input_sens[4], &input_sens[5]);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\tTACSGPBladeStiffened..testCriticalLocalShearLoad_LowAR:\n");
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testStiffenerCripplingLoad(
    TacsScalar epsilon, int printLevel) {
  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 7;
  TacsScalar p_input[n_input], x0[n_input], x[n_input], input_sens[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  if (printLevel != 0) {
    printf("\nTACSGPBladeStiffened..testStiffenerCripplingLoad start::\n");
    printf("--------------------------------------------------------\n\n");
  }

  // compute initial values
  x0[0] = 10.2412;                // D11
  x0[1] = 5.4323;                 // D22
  x0[2] = this->stiffenerHeight;  // sheight
  x0[3] = 1.24332;                // xi
  x0[4] = 2.4545;                 // rho0
  x0[5] = 0.2454;                 // genPoiss
  x0[6] = 1e-3;                   // zeta

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  this->stiffenerHeight = x0[2] - p_input[2] * epsilon;
  this->panelGPs->resetSavedData();
  f0 = computeStiffenerCripplingLoad(x[0], x[1], x[3], x[4], x[5], x[6]);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  this->stiffenerHeight = x0[2] + p_input[2] * epsilon;
  this->panelGPs->resetSavedData();
  f2 = computeStiffenerCripplingLoad(x[0], x[1], x[3], x[4], x[5], x[6]);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // now perform the adjoint sensitivity
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  this->stiffenerHeight = x0[2];
  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i];
  }
  this->panelGPs->resetSavedData();
  computeStiffenerCripplingLoadSens(
      p_output, x[0], x[1], x[3], x[4], x[5], x[6], &input_sens[0],
      &input_sens[1], &input_sens[2], &input_sens[3], &input_sens[4],
      &input_sens[5], &input_sens[6]);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("TACSGPBladeStiffened..testStiffenerCripplingLoad:\n");
    printf("\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t rel error = %.4e\n", TacsRealPart(relError));
  }

  return relError;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::testAllTests(
    TacsScalar epsilon, int printLevel) {
  // run each of the nondim parameter tests and aggregate the max among them
  const int n_tests = 8;
  TacsScalar relErrors[n_tests];
  memset(relErrors, 0, n_tests * sizeof(TacsScalar));

  if (printLevel != 0) {
    printf("\nTACSGPBladeStiffened..testAllTests start::\n");
    printf("========================================================\n");
    printf("========================================================\n\n");
  }

  relErrors[0] = testNondimensionalParameters(epsilon, printLevel);
  relErrors[1] = testAxialCriticalLoads(epsilon, printLevel);
  relErrors[2] = testShearCriticalLoads(epsilon, printLevel);
  relErrors[3] = testStiffenerCripplingLoad(epsilon, printLevel);
  relErrors[7] = testOtherTests(epsilon, printLevel);

  if (this->getAxialGP()) {
    if (printLevel != 0) {
      printf("\nAxialGP : testAllGPtests start::\n");
      printf("--------------------------------------------------------\n\n");
    }
    relErrors[4] = this->getAxialGP()->testAllGPTests(epsilon, printLevel);
  }
  if (this->getShearGP()) {
    if (printLevel != 0) {
      printf("\nShearGP : testAllGPtests start::\n");
      printf("--------------------------------------------------------\n\n");
    }
    relErrors[5] = this->getShearGP()->testAllGPTests(epsilon, printLevel);
  }
  if (this->getCripplingGP()) {
    if (printLevel != 0) {
      printf("\nCripplingGP : testAllGPtests start::\n");
      printf("--------------------------------------------------------\n\n");
    }
    relErrors[6] = this->getCripplingGP()->testAllGPTests(epsilon, printLevel);
  }

  // get max rel error among them
  TacsScalar maxRelError = 0.0;
  for (int i = 0; i < n_tests; i++) {
    if (TacsRealPart(relErrors[i]) > TacsRealPart(maxRelError)) {
      maxRelError = relErrors[i];
    }
  }

  // get max rel error among them
  if (printLevel != 0) {
    printf("========================================================\n");
    printf("\nTACSGPBladeStiffened..testAllTests full results::\n");
    printf("\ttestNondimensionalParameters = %.4e\n",
           TacsRealPart(relErrors[0]));
    printf("\ttestAxialCriticalLoads = %.4e\n", TacsRealPart(relErrors[1]));
    printf("\ttestShearCriticalLoads = %.4e\n", TacsRealPart(relErrors[2]));
    printf("\ttestStiffenerCripplingLoad = %.4e\n", TacsRealPart(relErrors[3]));
    if (this->getAxialGP()) {
      printf("\ttestAxialGP all tests = %.4e\n", TacsRealPart(relErrors[4]));
    }
    if (this->getShearGP()) {
      printf("\ttestShearGP all tests = %.4e\n", TacsRealPart(relErrors[5]));
    }
    if (this->getCripplingGP()) {
      printf("\ttestCripplingGp all tests = %.4e\n",
             TacsRealPart(relErrors[6]));
    }
    printf("\ttestOtherTests = %.4e\n", TacsRealPart(relErrors[7]));
    printf("\tOverall max rel error = %.4e\n\n", TacsRealPart(maxRelError));
  }

  return maxRelError;
}

void TACSGPBladeStiffenedShellConstitutive::setKSWeight(double ksWeight) {
  // call the superclass method to set the ksWeight for this const class
  TACSBladeStiffenedShellConstitutive::setKSWeight(ksWeight);
}
