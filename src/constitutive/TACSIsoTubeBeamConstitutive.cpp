#include "TACSIsoTubeBeamConstitutive.h"

TACSIsoTubeBeamConstitutive::TACSIsoTubeBeamConstitutive( TACSMaterialProperties *properties,
                                                          TacsScalar inner_diameter,
                                                          TacsScalar wall_thickness,
                                                          int inner_diameter_dv ){

}

int TACSIsoTubeBeamConstitutive::getDesignVarNums( int elemIndex,
                                                   int dvLen,
                                                   int dvNums[] ){
  return 0;
}


int TACSIsoTubeBeamConstitutive::setDesignVars( int elemIndex,
                                                int dvLen,
                                                const TacsScalar dvs[] ){
  return 0;
}

int TACSIsoTubeBeamConstitutive::getDesignVars( int elemIndex,
                                                int dvLen, TacsScalar dvs[] ){
  return 0;
}

int TACSIsoTubeBeamConstitutive::getDesignVarRange( int elemIndex, int dvLen,
                                                    TacsScalar lowerBound[],
                                                    TacsScalar upperBound[] ){
  return 0;
}

void TACSIsoTubeBeamConstitutive::evalMassMoments( int elemIndex,
                                                   const double pt[],
                                                   const TacsScalar X[],
                                                   TacsScalar moments[] ){
  // TacsScalar A = M_PI * ((inner + wall) * (inner + wall) - inner * inner);
  // TacsScalar Ia = M_PI ;

  // moments[0] = rho * A;
  // moments[1] = rho * Ia;
  // moments[2] = rho * Ia;
  // moments[3] = 0.0;
}

void TACSIsoTubeBeamConstitutive::addMassMomentsDVSens( int elemIndex,
                                                        const double pt[],
                                                        const TacsScalar X[],
                                                        const TacsScalar scale[],
                                                        int dvLen,
                                                        TacsScalar dfdx[] ){
  int index = 0;
  if (innerDV >= 0){
    index++;
  }
  if (wallDV >= 0){

  }
}

TacsScalar TACSIsoTubeBeamConstitutive::evalDensity( int elemIndex,
                                                     const double pt[],
                                                     const TacsScalar X[] ){

}

void TACSIsoTubeBeamConstitutive::addDensityDVSens( int elemIndex,
                                                    TacsScalar scale,
                                                    const double pt[],
                                                    const TacsScalar X[],
                                                    int dvLen, TacsScalar dfdx[] ){
}

void TACSIsoTubeBeamConstitutive::evalStress( int elemIndex,
                                              const double pt[],
                                              const TacsScalar X[],
                                              const TacsScalar strain[],
                                              TacsScalar stress[] ){
}

void TACSIsoTubeBeamConstitutive::evalTangentStiffness( int elemIndex,
                                                        const double pt[],
                                                        const TacsScalar X[],
                                                        TacsScalar C[] ){
}

void TACSIsoTubeBeamConstitutive::addStressDVSens( int elemIndex,
                                                   TacsScalar scale,
                                                   const double pt[],
                                                   const TacsScalar X[],
                                                   const TacsScalar strain[],
                                                   const TacsScalar psi[],
                                                   int dvLen, TacsScalar dfdx[] ){
}

TacsScalar TACSIsoTubeBeamConstitutive::evalFailure( int elemIndex,
                                                     const double pt[],
                                                     const TacsScalar X[],
                                                     const TacsScalar strain[] ){
}

TacsScalar TACSIsoTubeBeamConstitutive::evalFailureStrainSens( int elemIndex,
                                                               const double pt[],
                                                               const TacsScalar X[],
                                                               const TacsScalar strain[],
                                                               TacsScalar sens[] ){

}

void TACSIsoTubeBeamConstitutive::addFailureDVSens( int elemIndex,
                                                    TacsScalar scale,
                                                    const double pt[],
                                                    const TacsScalar X[],
                                                    const TacsScalar strain[],
                                                    int dvLen, TacsScalar dfdx[] ){

}

TacsScalar TACSIsoTubeBeamConstitutive::evalDesignFieldValue( int elemIndex,
                                                              const double pt[],
                                                              const TacsScalar X[],
                                                              int index ){
  if (index == 0){
    return inner;
  }
  else if (index == 1){
    return wall;
  }
  return 0.0;
}
