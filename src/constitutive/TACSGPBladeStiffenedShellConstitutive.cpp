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
    // use the closed-form solution to compute the critical global axial load
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
    // use the closed-form solution to compute the critical global axial load
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
      *rho_0sens += neg_N11crits_sens[m1-1] * dim_factor * ((1.0 + gamma) * -2.0 * pow(m1/rho_0, 2.0) / rho_0 + pow(m1/rho_0, -2.0) * 2.0 / rho_0 );
      *xisens += neg_N11crits_sens[m1-1] * dim_factor * 2.0;
      *gammasens += neg_N11crits_sens[m1-1] * dim_factor * pow(m1/rho_0, 2.0);
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
    // use the closed-form solution to compute the critical global axial load
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
    // use the closed-form solution to compute the critical global axial load
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
      *rho_0sens += neg_N11crits_sens[m1-1] * dim_factor * (-2.0 * pow(m1/rho_0, 2.0) / rho_0 + pow(m1/rho_0, -2.0) * 2.0 / rho_0 );
      *xisens += neg_N11crits_sens[m1-1] * dim_factor * 2.0;
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
    // use the closed-form solution to compute the critical global axial load
    // forward analysis part here
  }
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
    // use the closed-form solution to compute the critical global axial load
    // forward analysis part here
  }
}

TacsScalar TACSGPBladeStiffenedShellConstitutive::computeCriticalLocalShearLoad(
    const TacsScalar D11, const TacsScalar D22, const TacsScalar xi) {
  if (this->useGPs) {
    // use Gaussian processes to compute the critical global axial load
    return -1.0; // not written this part yet.
  } else {
    // use the closed-form solution to compute the critical global axial load
    // forward analysis part here
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
    // use the closed-form solution to compute the critical global axial load
    // forward analysis part here
  }
}