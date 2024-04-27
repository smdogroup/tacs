/*
=====================================================================================================
Blade-Stiffened Shell Constitutive Model using Gaussian Process Machine Learning Buckling Constraints
=====================================================================================================
@File    :   TACSGPBladeStiffenedShellConstutive.cpp
@Date    :   2024/04/24
@Author  :   Sean Phillip Engelstad, Alasdair Christian Gray
@Description : Constitutive model for a blade-stiffened shell. Based on the FSDT blade models
adopted by Alasdair Christian Gray in the TACSBladeStiffenedShellConstutitutive.h class and the
original one by Graeme Kennedy. Gaussian Processes for Machine Learning are used for the buckling
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
    TACSOrthotropicPly* _panelPly, TACSOrthotropicPly* _stiffenerPly,
    TacsScalar _kcorr, TacsScalar _panelLength, TacsScalar _panelWidth,
    TacsScalar _stiffenerPitch, int _stiffenerPitchNum, TacsScalar _panelThick,
    int _panelThickNum, int _numPanelPlies, TacsScalar _panelPlyAngles[],
    TacsScalar _panelPlyFracs[], int _panelPlyFracNums[],
    TacsScalar _stiffenerHeight, int _stiffenerHeightNum,
    TacsScalar _stiffenerThick, int _stiffenerThickNum, int _numStiffenerPlies,
    TacsScalar _stiffenerPlyAngles[], TacsScalar _stiffenerPlyFracs[],
    int _stiffenerPlyFracNums[], TacsScalar _flangeFraction) {
  
  this->panelWidth = _panelWidth;
  // call the superclass constructor except with panelLength DV turned off
  TACSBladeStiffenedShellConstitutive(_panelPly, _stiffenerPly, _kcorr, _panelLength, -1, _stiffenerPitch, _stiffenerPitchNum, _panelThick,
  _panelThickNum, _numPanelPlies, _panelPlyAngles, _panelPlyFracs, _panelPlyFracNums, _stiffenerHeight, _stiffenerHeightNum, _stiffenerThick,
  _stiffenerThickNum, _numStiffenerPlies, _stiffenerPlyAngles, _stiffenerPlyFracs, _stiffenerPlyFracNums, _flangeFraction);
  
}

// ==============================================================================
// Destructor
// ==============================================================================
TACSGPBladeStiffenedShellConstitutive::~TACSGPBladeStiffenedShellConstitutive() {
  ~TACSBladeStiffenedShellConstitutive(); // call superclass destructor
}

// ==============================================================================
// Override Failure constraint and sensitivities
// ==============================================================================

// Compute the failure values for each failure mode of the stiffened panel
TacsScalar TACSGPBladeStiffenedShellConstitutive::computeFailureValues(
    const TacsScalar e[], TacsScalar fails[]) {
  // --- Panel material failure ---
  // TODO : need to complete this function
}

// Evaluate the derivative of the failure criteria w.r.t. the strain
TacsScalar TACSGPBladeStiffenedShellConstitutive::evalFailureStrainSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar e[], TacsScalar sens[]) {
      // TODO : need to complete this function
}

// Add the derivative of the failure criteria w.r.t. the design variables
void TACSGPBladeStiffenedShellConstitutive::addFailureDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], int dvLen, TacsScalar dfdx[]) {
      // TODO : need to complete this function
}

// ==============================================================================
// Buckling functions
// ==============================================================================

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeAffineAspectRatioSens(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar a,
    const TacsScalar b, TacsScalar* D11sens, TacsScalar* D22sens,
    TacsScalar* asens, TacsScalar* bsens) {

  // compute the derivatives of the affine aspect ratio and return the affine aspect ratio
  TacsScalar rho_0 = computeAffineAspectRatio(D11,D22,a,b);
  // where rho_0 = a/b * (D22/D11)**0.25

  // use power series rules d(x^p) = p * (x^p) / x to cleanly differentiate the expression
  *asens = rho_0 / a;
  *bsens = -1.0 * rho_0 / b;
  *D11sens = -0.25 * rho_0 / D11;
  *D22sens = 0.25 * rho_0 / D22;

  return rho_0;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeGeneralizedRigiditySens(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar D12,
    const TacsScalar D66, TacsScalar* D11sens, TacsScalar* D22sens,
    TacsScalar* D12sens, TacsScalar* D66sens) {
  
  // compute the derivatives of the generalized rigidity xi
  TacsScalar denominator = sqrt(D11 * D22);
  TacsScalar xi = computeGeneralizedRigidity(D11,D22,D12,D66);
  // so that xi = (D12 + 2 * D66) / sqrt(D11*D22)

  // compute the sensitivities
  *D12sens = 1.0 / denominator;
  *D66sens = 2.0 / denominator;
  *D11sens = -0.5 * xi / D11;
  *D22sens = -0.5 * xi / D22;

  return xi;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeGeneralizedPoissonsRatioSens(
    const TacsScalar D12, const TacsScalar D66, 
    TacsScalar* D12sens, TacsScalar* D66sens) {
  
  // compute derivatives of the generalized poisson's ratio
  TacsScalar eps = computeGeneralizedPoissonsRatio(D12,D66);
  // where eps = (D12 + 2 * D66) / D12

  *D12sens = - 2.0 * D66 / D12 / D12;
  *D66sens = 2.0 / D12;

  return eps;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeStiffenerAreaRatio() {
  // get effective moduli for the panel and stiffener
  TacsScalar E1s, E1p, _;
  // first the effective modulus of the panel/plate
  // need to add derivatives w.r.t. panel plies here then..
  this->computeEffectiveModulii(this->numPanelPlies, this->panelQMats,
                                this->panelPlyFracs, &E1p, &_);

  // then the stiffener
  this->computeEffectiveModulii(this->numStiffenerPlies, this->stiffenerQMats,
                                this->stiffenerPlyFracs, &E1s, &_);

  TacsScalar As = this->computeStiffenerArea();
  return E1s * As / (E1p * this->stiffenerPitch * this->panelThick);
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeStiffenerAreaRatioSens(
    TacsScalar* sthickSens, TacsScalar* sheightSens, TacsScalar* spitchSens, TacsScalar* pthickSens) {
  // get effective moduli for the panel and stiffener
  TacsScalar delta = this->computeStiffenerAreaRatio();

  // TODO : may need to add the derivatives w.r.t. panel ply fractions for E1p computation later..
  *sthickSens = delta / this->stiffenerThick;
  *sheightSens = delta / this->stiffenerHeight;
  *spitchSens = -1.0 * delta / this->stiffenerPitch;
  *pthickSens = -1.0 * delta / this->panelThick;

  return delta;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeStiffenerStiffnessRatio(TacsScalar D11) {
  // get effective moduli for the panel and stiffener
  TacsScalar E1s, E1p, _;
  // first the effective modulus of the panel/plate
  // need to add derivatives w.r.t. panel plies here then..
  this->computeEffectiveModulii(this->numPanelPlies, this->panelQMats,
                                this->panelPlyFracs, &E1p, &_);

  // then the stiffener
  this->computeEffectiveModulii(this->numStiffenerPlies, this->stiffenerQMats,
                                this->stiffenerPlyFracs, &E1s, &_);

  // get the stiffener bending inertia in 11-direction
  // TODO : double check if this bending stiffness calculation is correct..
  TacsScalar Is = this->computeStiffenerIzz();

  return E1s * Is / D11 / this->stiffenerThick;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeStiffenerStiffnessRatioSens(
    TacsScalar D11, TacsScalar* D11sens, TacsScalar* sthickSens, 
    TacsScalar* sheightSens, TacsScalar* spitchSens) {

  // use power series rules and the forward state to differentiate
  TacsScalar gamma = computeStiffenerStiffnessRatio(D11); // this is the stiffener stiffness ratio (final forward state)

  // intermediate states and sensitivities
  TacsScalar Is = this->computeStiffenerIzz();
  TacsScalar dI_dsthick, dI_dsheight;
  this->computeStiffenerIzzSens(dI_dsthick, dI_dsheight);

  // TODO : may need to add stiffener ply fraction derivatives for E1p computation, for now ignoring
  *D11sens = -1.0 * gamma / D11;
  *sthickSens = gamma / Is * dI_dsthick;
  *sheightSens = gamma / Is * dI_dsheight;
  *spitchSens = -1.0 * gamma / this->stiffenerPitch;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeCriticalGlobalAxialLoad(
  const TacsScalar D11, const TacsScalar D22, const TacsScalar b, const TacsScalar delta,
  const TacsScalar rho_0, const TacsScalar xi, const TacsScalar gamma) {
  
  if (this->useGPs) {
    // use Gaussian processes to compute the critical global axial load
    return -1.0; // not written this part yet.
  } else {
    // use the CPT closed-form solution to compute the critical global axial load
    TacsScalar neg_N11crits[this->NUM_CF_MODES];
    for (int m1 = 1; m1 < this->NUM_CF_MODES+1; m1++) {
      TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) / b / b / (1.0 + delta);
      TacsScalar nondim_factor = (1.0 + gamma) * pow(m1/rho_0,2.0) + pow(m1/rho_0,-2.0) + 2.0 * xi;
      neg_N11crits[m1-1] = -1.0 * dim_factor * nondim_factor; // negated only because we have to do KS min aggregate later
    }

    // compute KS aggregation for -N11crit for each mode then negate again (because we want minimum N11crit so maximize negative N11crit)
    TacsScalar neg_N11crit = ksAggregation(neg_N11crits, this->NUM_CF_MODES, this->ksWeight);
    return -1.0 * neg_N11crit;
  }

}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeCriticalGlobalAxialLoadSens(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar b,
    const TacsScalar delta, const TacsScalar rho_0, const TacsScalar xi, 
    const TacsScalar gamma, TacsScalar* D11sens, TacsScalar* D22sens, 
    TacsScalar* bsens, TacsScalar* deltasens, TacsScalar* rho_0sens, 
    TacsScalar* xisens, TacsScalar* gammasens) {
  
  if (this->useGPs) {
    // use Gaussian processes to compute the critical global axial load
    return -1.0; // not written this part yet.
  } else {
    // use the CPT closed-form solution to compute the critical global axial load
    // forward analysis part here
    TacsScalar neg_N11crits[this->NUM_CF_MODES];
    for (int m1 = 1; m1 < this->NUM_CF_MODES+1; m1++) {
      TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) / b / b / (1.0 + delta);
      TacsScalar nondim_factor = (1.0 + gamma) * pow(m1/rho_0,2.0) + pow(m1/rho_0,-2.0) + 2.0 * xi;
      neg_N11crits[m1-1] = -1.0 * dim_factor * nondim_factor; // negated only because we have to do KS min aggregate later
    }

    // compute KS aggregation sensitivity
    TacsScalar neg_N11crits_sens[this->NUM_CF_MODES];
    TacsScalar neg_N11crit = ksAggregationSens(neg_N11crits, this->NUM_CF_MODES, this->ksWeight, neg_N11crits_sens);

    // compute sensitivities here
    *D11sens = *D22sens = *bsens = *deltasens = *rho_0sens = *xisens = *gammasens = 0.0;
    for (int m1 = 1; m1 < this->NUM_CF_MODES+1; m1++) {
      // forward analysis states
      TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) / b / b / (1.0 + delta);
      TacsScalar nondim_factor = (1.0 + gamma) * pow(m1/rho_0,2.0) + pow(m1/rho_0,-2.0) + 2.0 * xi;
      neg_N11crits[m1-1] = -1.0 * dim_factor * nondim_factor; // negated only because we have to do KS min aggregate later

      // update sensitivities (left factor is dKS/dv_i, right factor is dv_i / dx)
      *D11sens += neg_N11crits_sens[m1-1] * (0.5 * neg_N11crits[m1-1] / D11);
      *D22sens += neg_N11crits_sens[m1-1] * (0.5 * neg_N11crits[m1-1] / D22);
      *bsens += neg_N11crits_sens[m1-1] * (-2.0 * neg_N11crits[m1-1] / b);
      *deltasens += neg_N11crits_sens[m1-1] * (-1.0 * neg_N11crits[m1-1] / (1.0 + delta));
      *rho_0sens += neg_N11crits_sens[m1-1] * -dim_factor * ((1.0 + gamma) * -2.0 * pow(m1/rho_0, 2.0) / rho_0 + pow(m1/rho_0, -2.0) * 2.0 / rho_0 );
      *xisens += neg_N11crits_sens[m1-1] * -dim_factor * 2.0;
      *gammasens += neg_N11crits_sens[m1-1] * -dim_factor * pow(m1/rho_0, 2.0);
    }
    return -1.0 * neg_N11crit;
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeCriticalLocalAxialLoad(
  const TacsScalar D11, const TacsScalar D22, const TacsScalar rho_0, 
  const TacsScalar xi) {
  
  if (this->useGPs) {
    // use Gaussian processes to compute the critical global axial load
    return -1.0; // not written this part yet.
  } else {
    // use the CPT closed-form solution to compute the critical global axial load
    TacsScalar neg_N11crits[this->NUM_CF_MODES];
    for (int m1 = 1; m1 < this->NUM_CF_MODES+1; m1++) {
      TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) / this->stiffenerPitch / this->stiffenerPitch;
      TacsScalar nondim_factor = pow(m1/rho_0,2.0) + pow(m1/rho_0,-2.0) + 2.0 * xi;
      neg_N11crits[m1-1] = -1.0 * dim_factor * nondim_factor; // negated only because we have to do KS min aggregate later
    }

    // compute KS aggregation for -N11crit for each mode then negate again (because we want minimum N11crit so maximize negative N11crit)
    TacsScalar neg_N11crit = ksAggregation(neg_N11crits, this->NUM_CF_MODES, this->ksWeight);
    return -1.0 * neg_N11crit;
  }

}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeCriticalLocalAxialLoadSens(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar rho_0, 
    const TacsScalar xi, TacsScalar* D11sens, TacsScalar* D22sens, 
    TacsScalar* rho_0sens, TacsScalar* xisens, TacsScalar* spitchsens) {
  
  if (this->useGPs) {
    // use Gaussian processes to compute the critical global axial load
    return -1.0; // not written this part yet.
  } else {
    // use the CPT closed-form solution to compute the critical global axial load
    // forward analysis part here
    TacsScalar neg_N11crits[this->NUM_CF_MODES];
    for (int m1 = 1; m1 < this->NUM_CF_MODES+1; m1++) {
      TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) / this->stiffenerPitch / this->stiffenerPitch;
      TacsScalar nondim_factor = pow(m1/rho_0,2.0) + pow(m1/rho_0,-2.0) + 2.0 * xi;
      neg_N11crits[m1-1] = -1.0 * dim_factor * nondim_factor; // negated only because we have to do KS min aggregate later
    }

    // compute KS aggregation sensitivity
    TacsScalar neg_N11crits_sens[this->NUM_CF_MODES];
    TacsScalar neg_N11crit = ksAggregationSens(neg_N11crits, this->NUM_CF_MODES, this->ksWeight, neg_N11crits_sens);

    // compute sensitivities here
    *D11sens = *D22sens = *spitchsens = *rho_0sens = *xisens = 0.0;
    for (int m1 = 1; m1 < this->NUM_CF_MODES+1; m1++) {
      // forward analysis states
      TacsScalar dim_factor = M_PI * M_PI * sqrt(D11 * D22) / this->stiffenerPitch / this->stiffenerPitch;
      TacsScalar nondim_factor = pow(m1/rho_0,2.0) + pow(m1/rho_0,-2.0) + 2.0 * xi;
      neg_N11crits[m1-1] = -1.0 * dim_factor * nondim_factor; // negated only because we have to do KS min aggregate later

      // update sensitivities (left factor is dKS/dv_i, right factor is dv_i / dx)
      *D11sens += neg_N11crits_sens[m1-1] * (0.5 * neg_N11crits[m1-1] / D11);
      *D22sens += neg_N11crits_sens[m1-1] * (0.5 * neg_N11crits[m1-1] / D22);
      *rho_0sens += neg_N11crits_sens[m1-1] * -dim_factor * (-2.0 * pow(m1/rho_0, 2.0) / rho_0 + pow(m1/rho_0, -2.0) * 2.0 / rho_0 );
      *xisens += neg_N11crits_sens[m1-1] * -dim_factor * 2.0;
      *spitchsens += neg_N11crits_sens[m1-1] * (-2.0 * neg_N11crits[m1-1] / this->stiffenerPitch);
    }
    return -1.0 * neg_N11crit;
  }
}

  
TacsScalar TACSGPBladeStiffenedShellConstitutive::computeCriticalGlobalShearLoad(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar b,
    const TacsScalar xi, const TacsScalar gamma,) {
  if (this->useGPs) {
    // use Gaussian processes to compute the critical global axial load
    return -1.0; // not written this part yet.
  } else {
    // use the CPT closed-form solution to compute the critical global axial load
    // no mode switching in this solution.. (only accurate for higher aspect ratios => hence the need for machine learning for the actual solution)
    TacsScalar lam1, lam2; // lam1bar, lam2bar values
    nondimShearParams(xi,gamma,&lam1,&lam2);
    TacsScalar dim_factor = M_PI * M_PI * pow(D11 * D22 * D22 * D22, 0.25) / b / b;
    TacsScalar nondim_factor = (1.0 + pow(lam1, 4.0) + 6.0 * pow(lam1 * lam2, 2.0) + pow(lam2, 4.0) + 2.0 * xi) / (2.0 * lam1 * lam1 * lam2);
    return dim_factor * nondim_factor; // aka N12_crit from CPT closed-form solution
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::nondimShearParams(const TacsScalar xi, const TacsScalar gamma, TacsScalar* lam1bar, TacsScalar* lam2bar) {
  // need to iterate over lam2 with the Newton's method
  TacsScalar lam2bar_sq = 0.0; // starting guess for lambda2_bar^2

  // Newton iteration for lam2bar squared
  while (abs(TacsRealPart(lam2Constraint(lam2bar_sq, xi, gamma))) > 1e-10) {
    lam2bar_sq -= lam2Constraint(lam2bar_sq, xi, gamma) / lam2ConstraintDeriv(lam2bar_sq, xi, gamma);
  }

  // now compute lam1_bar, lam2_bar
  *lam1bar = pow(1.0 + 2.0 * lam2bar_sq * xi + lam2bar_sq * lam2bar_sq + gamma, 0.25);
  *lam2bar = pow(lam2bar_sq, 0.5);
}
TacsScalar TACSGPBladeStiffenedShellConstitutive::lam2Constraint(const TacsScalar lam2sq, const TacsScalar xi, const TacsScalar gamma) {
  // compute the residual of the combined lam1bar, lam2bar constraint but on lam2bar^2
  TacsScalar lam1bar = pow(1.0 + 2.0 * lam2sq * xi + lam2sq * lam2sq + gamma, 0.25);
  TacsScalar lam1sq = lam1bar * lam1bar;
  TacsScalar lam14 = lam1sq * lam1sq;
  return lam2sq + lam1sq + xi / 3.0 - sqrt((3.0 + xi) / 9.0 + 4.0 / 3.0 * lam1sq * xi + 4.0/3.0 * lam14);
}
TacsScalar TACSGPBladeStiffenedShellConstitutive::lam2ConstraintDeriv(const TacsScalar lam2sq, const TacsScalar xi, const TacsScalar gamma) {
  // compute the residual derivatives for the lam2bar constraint above w.r.t. the lam2sq input about a lam2sq input
  TacsScalar dfdlam2sq = 1.0;
  TacsScalar temp = 1.0 + 2.0 * lam2sq * xi + lam2sq * lam2sq + gamma;
  TacsScalar lam1 = pow(temp, 0.25);
  TacsScalar lam1sq = lam1 * lam1;
  TacsScalar lam14 = lam1sq * lam1sq;

  TacsScalar term2 = sqrt((3.0 + xi) / 9.0 + 4.0 / 3.0 * lam1sq * xi + 4.0/3.0 * lam14);
  TacsScalar dlam1_dlam2sq = lam1 * 0.25 / temp * (2.0 + 2.0 * lam2sq);
  TacsScalar dfdlam1 = 2 * lam1 - 0.5 / term2 * 4.0 / 3.0 * (2.0 * lam1 * xi + 4.0 * lam1 * lam1sq);
  dfdlam2sq += dfdlam1 * dlam1_dlam2sq;
  return dfdlam2sq;
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::nondimShearParamsSens(const TacsScalar xi, const TacsScalar gamma,
  TacsScalar* lam1bar, TacsScalar* lam2bar,
  TacsScalar* dl1xi, TacsScalar* dl1gamma, TacsScalar* dl2xi, TacsScalar* dl2gamma) {

  // get the lam1, lam2 from Newton's method
  TacsScalar lam1, lam2;
  nondimShearParams(xi,gamma,&lam1,&lam2);

  // also send out the lam1bar, lam2bar again
  *lam1bar = lam1;
  *lam2bar = lam2;

  // differentiate the nonlinear constraints from nondimShearParam subroutine
  // sys eqns [A,B;C,D] * [lam1bar_dot, lam2bar_dot] = [E,F] for each of the two derivatives

  TacsScalar exp1 = 1.0 + 2.0 * lam2 * lam2 * xi + pow(lam2, 4.0) + gamma;
  TacsScalar dexp1lam2 = 4.0 * lam2 * xi + 4.0 * lam2 * lam2 * lam2;
  TacsScalar dexp1xi = 2.0 * lam2 * lam2;
  TacsScalar dexp1gamma = 1.0;
  TacsScalar exp2 = (3.0 + xi) / 9.0 + 4.0/3.0 * (lam1 * lam1 * xi + pow(lam1, 4.0));
  TacsScalar dexp2lam1 = 4.0/3.0 * (2.0 * lam1 * xi + 4.0 * lam1 * lam1 * lam1);
  TacsScalar dexp2xi = 1.0/9.0;
  TacsScalar dexp2gamma = 0.0;

  // first for the xi sensitivities
  TacsScalar A1, B1, C1, D1, E1, F1;
  A1 = 1.0;
  B1 = -0.25 * lam1 / exp1 * dexp1lam2;
  E1 = 0.25 * lam1 / exp1 * dexp1xi;
  C1 = 2.0 * lam1 - 0.5 * pow(exp2, -0.5) * dexp2lam1;
  D1 = 2.0 * lam2;
  F1 = -1.0/3.0 + 0.5 * pow(exp2, -0.5) * dexp2xi;
  *dl1xi = (D1 * E1 - B1 * F1) / (A1 * D1 - B1 * C1);
  *dl2xi = (A1 * F1 - C1 * E1) / (A1 * D1 - B1 * C1);

  // then for the gamma sensitivities
  TacsScalar A2, B2, C2, D2, E2, F2;
  A2 = A1;
  B2 = B1;
  E2 = 0.25 * lam1 / exp1 * dexp1gamma;
  C2 = C1;
  D2 = D1;
  F2 = -1.0/3.0 + 0.5 * pow(exp2, -0.5) * dexp2gamma;
  *dl1gamma = (D2 * E2 - B2 * F2) / (A2 * D2 - B2 * C2);
  *dl2gamma = (A2 * F2 - C2 * E2) / (A2 * D2 - B2 * C2);
  
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeCriticalGlobalShearLoadSens( 
    const TacsScalar D11, const TacsScalar D22, const TacsScalar b,
    const TacsScalar xi, const TacsScalar gamma,
    TacsScalar* D11sens, TacsScalar* D22sens, TacsScalar* bsens,
    TacsScalar* xisens, TacsScalar* gammasens) {
  if (this->useGPs) {
    // use Gaussian processes to compute the critical global axial load
    return -1.0; // not written this part yet.
  } else {
    // use the CPT closed-form solution to compute the critical global axial load
    // no mode switching in this solution.. (only accurate for higher aspect ratios => hence the need for machine learning for the actual solution)
    TacsScalar lam1, lam2; // lam1bar, lam2bar values
    TacsScalar dl1xi, dl2xi, dl1gamma, dl2gamma;

    // compute the derivatives of the nondimensional constraints
    nondimShearParamsSens(xi,gamma,&lam1,&lam2,&dl1xi,&dl1gamma,&dl2xi,&dl2gamma);

    // compute forward analysis states involved in the N12crit load
    TacsScalar dim_factor = M_PI * M_PI * pow(D11 * D22 * D22 * D22, 0.25) / b / b;
    TacsScalar num = 1.0 + pow(lam1, 4.0) + 6.0 * pow(lam1 * lam2, 2.0) + pow(lam2, 4.0) + 2.0 * xi;
    TacsScalar den = 2.0 * lam1 * lam1 * lam2;
    TacsScalar nondim_factor = num / den;
    TacsScalar N12crit = dim_factor * nondim_factor;

    // sensitivities for the non_dim factor

    TacsScalar dNDlam1 = (4.0 * pow(lam1, 3.0) + 12.0 * lam1 * lam2 * lam2) / den - num * 4.0 * lam1 * lam2 / den / den;
    TacsScalar dNDlam2 = (4.0 * pow(lam2, 3.0) + 12.0 * lam2 * lam1 * lam1) / den - num * 2.0 * lam1 * lam1 / den / den;

    // compute the overall sensitivities
    *D11sens = N12crit * 0.25 / D11;
    *D22sens = N12crit * 0.75 / D22;
    *bsens = N12crit * -2.0 / b;
    *xisens = dim_factor * (dNDlam1 * dl1xi + dNDlam2 * dl2xi + 2.0 / den);
    *gammasens = dim_factor * (dNDlam1 * dl1gamma + dNDlam2 * dl2gamma);

    // return N12crit from closed-form solution
    return N12crit;
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeCriticalLocalShearLoad(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar xi) {
  if (this->useGPs) {
    // use Gaussian processes to compute the critical global axial load
    return -1.0; // not written this part yet.
  } else {
    // use the CPT closed-form solution to compute the critical global axial load
    // no mode switching in this solution.. (only accurate for higher aspect ratios => hence the need for machine learning for the actual solution)
    TacsScalar lam1, lam2; // lam1bar, lam2bar values
    nondimShearParams(xi,0.0,&lam1,&lam2);
    TacsScalar dim_factor = M_PI * M_PI * pow(D11 * D22 * D22 * D22, 0.25) / this->stiffenerPitch / this->stiffenerPitch;
    TacsScalar nondim_factor = (1.0 + pow(lam1, 4.0) + 6.0 * pow(lam1 * lam2, 2.0) + pow(lam2, 4.0) + 2.0 * xi) / (2.0 * lam1 * lam1 * lam2);
    return dim_factor * nondim_factor; // aka N12_crit from CPT closed-form solution
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeCriticalLocalShearLoadSens( 
    const TacsScalar D11, const TacsScalar D22, const TacsScalar xi, 
    TacsScalar* D11sens, TacsScalar* D22sens, TacsScalar* xisens,
    TacsScalar* spitchsens) {
  if (this->useGPs) {
    // use Gaussian processes to compute the critical global axial load
    return -1.0; // not written this part yet.
  } else {
    // use the CPT closed-form solution to compute the critical global axial load
    // no mode switching in this solution.. (only accurate for higher aspect ratios => hence the need for machine learning for the actual solution)
    TacsScalar lam1, lam2; // lam1bar, lam2bar values
    TacsScalar dl1xi, dl2xi, _dl1gamma, _dl2gamma; // gamma derivatives are private since unused here (gamma=0 input)

    // compute the derivatives of the nondimensional constraints
    nondimShearParamsSens(xi,0.0,&lam1,&lam2,&dl1xi,&_dl1gamma,&dl2xi,&_dl2gamma);

    // compute forward analysis states involved in the N12crit load
    TacsScalar dim_factor = M_PI * M_PI * pow(D11 * D22 * D22 * D22, 0.25) / this->stiffenerPitch / this->stiffenerPitch;
    TacsScalar num = 1.0 + pow(lam1, 4.0) + 6.0 * pow(lam1 * lam2, 2.0) + pow(lam2, 4.0) + 2.0 * xi;
    TacsScalar den = 2.0 * lam1 * lam1 * lam2;
    TacsScalar nondim_factor = num / den;
    TacsScalar N12crit = dim_factor * nondim_factor;

    // sensitivities for the non_dim factor

    TacsScalar dNDlam1 = (4.0 * pow(lam1, 3.0) + 12.0 * lam1 * lam2 * lam2) / den - num * 4.0 * lam1 * lam2 / den / den;
    TacsScalar dNDlam2 = (4.0 * pow(lam2, 3.0) + 12.0 * lam2 * lam1 * lam1) / den - num * 2.0 * lam1 * lam1 / den / den;

    // compute the overall sensitivities
    *D11sens = N12crit * 0.25 / D11;
    *D22sens = N12crit * 0.75 / D22;
    *spitchsens = N12crit * -2.0 / this->stiffenerPitch;
    *xisens = dim_factor * (dNDlam1 * dl1xi + dNDlam2 * dl2xi + 2.0 / den);

    // return N12crit from closed-form solution
    return N12crit;
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeStiffenerCripplingLoad(
  const TacsScalar D11, const TacsScalar D22, const TacsScalar xi, const TacsScalar genPoiss) {
  if (this->useGPs) {
    // use Gaussian processes to compute the critical global axial load
    return -1.0; // not written this part yet.
  } else {
    // use the literature CPT closed-form solution for approximate stiffener crippling, not a function of aspect ratio
    TacsScalar dim_factor = M_PI * M_PI * sqrt(D11*D22) / this->stiffenerHeight / this->stiffenerHeight;
    TacsScalar nondim_factor = (0.476 - 0.56 * (genPoiss - 0.2)) * xi;
    return dim_factor * nondim_factor;
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeStiffenerCripplingLoadSens( 
  const TacsScalar D11, const TacsScalar D22, const TacsScalar xi, 
  const TacsScalar genPoiss, TacsScalar* D11sens, TacsScalar* D22sens,
  TacsScalar* xisens, TacsScalar* genPoiss_sens, TacsScalar* sheightsens) {
  if (this->useGPs) {
    // use Gaussian processes to compute the critical global axial load
    return -1.0; // not written this part yet.
  } else {
    // use the literature CPT closed-form solution for approximate stiffener crippling, not a function of aspect ratio
    TacsScalar dim_factor = M_PI * M_PI * sqrt(D11*D22) / this->stiffenerHeight / this->stiffenerHeight;
    TacsScalar nondim_factor = (0.476 - 0.56 * (genPoiss - 0.2)) * xi;
    TacsScalar N11crit = dim_factor * nondim_factor;

    // compute the derivatives
    *D11sens = N11crit * 0.5 / D11;
    *D22sens = N11crit * 0.5 / D22;
    *xisens = N11crit / xi;
    *genPoiss_sens = N11crit / nondim_factor * -0.56 * xi;
    *sheightsens = N11crit * -2.0 / this->stiffenerHeight;

    // output the critical load here
    return N11crit;
  }
}
