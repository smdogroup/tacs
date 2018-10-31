#include <stdio.h>
#include "rectangleEBStiffness.h"
#include "YSlibrary.h"

/*
  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

const char * rectangleEBStiffness::constName = "rectangleEBStiffness";

const char * rectangleEBStiffness::constitutiveName() const {
  return constName;
}

rectangleEBStiffness::rectangleEBStiffness( TacsScalar _rho, TacsScalar _E, TacsScalar _G,
                                TacsScalar _yield_stress,
                                TacsScalar _height,
                                TacsScalar _thickness,
                                int _height_num, int _thickness_num,
			                          TacsScalar _ref_dir[3],
				                        TacsScalar z_offset, TacsScalar y_offset,
                                enum EBBeamReferenceDirection _ref_type){

  /*
  The local cross-section coordinate axes are the lowercase letters x, y, and
  z. We reserve X, Y, and Z for the global coordinate axes. The axis of the beam
  element is determined based on the positions of the nodes defining the element
  boundaries. The orientation of the cross-section is determined using the
  _ref_dir and _ref_type inputs.

  If _ref_type == STRONG_AXIS:
  ----------------------------
    - x is along the axis of the beam
    - z is perpendicular to the reference direction and the x-axis (_ref_dir cross x)
    - y is perpendicular to the x-axis and the z-axis (x cross z)

  If _ref_type == WEAK_AXIS (default):
  ------------------------------------
    - x is along the axis of the beam
    - y is perpendicular to the x-axis and the reference direction (x cross _ref_dir)
    - z is perpendicular to the y-axis and the x-axis (y cross x)

        No offset                                     With offset
        ---------                                     -----------
              z
              ^
         _____|______                             ______________________ 1.0
        |     |     |                            |                     |
        |     |     |                            |                     |
              |     |                            |                     |
        |     |     |                            | z                   |
        |     |     |                            | ^                   |
        |     |     |                            | |                   |
        |     |     |                            | |                   |
        |     |     |                            | |                   |
    y <-------o     |                            | |        x  ---     | 0.0
        |           |                            | |            ^      |
        |           |                            | |            |      |
        |           |                            | |           oz      |
        |           |                            | |            |      |
        |           |                            | |            v      |
        |           |                y <-----------o           ---     |
        |           |                            | |<--oy-->|          |
        |           |                            |                     |
        |___________|                            |_____________________|-1.0
                                                 1.0        0.0       -1.0
  */
  height = _height;
  thickness = _thickness;

  thickness_num = _thickness_num;
  height_num = _height_num;

  ref_dir[0] = _ref_dir[0];
  ref_dir[1] = _ref_dir[1];
  ref_dir[2] = _ref_dir[2];
  ref_type = _ref_type;

  num_dvs = 0;
  if (height_num >= 0){
    dv_nums[num_dvs] = height_num;
    num_dvs++;
  }
  if (thickness_num >= 0){
    dv_nums[num_dvs] = thickness_num;
    num_dvs++;
  }
  oy = y_offset;
  oz = z_offset;

  lb_thickness = 0.0;
  ub_thickness = 1.0;

  lb_height = 0.0;
  ub_height = 1.0;

  rho = _rho;
  yield_stress = _yield_stress;
  E = _E;
  G = _G;
  ks_weight = 70.0;
}

void rectangleEBStiffness::setThicknessBounds( TacsScalar lb, TacsScalar ub ){
  lb_thickness = lb;
  ub_thickness = ub;
}

void rectangleEBStiffness::setHeightBounds( TacsScalar lb, TacsScalar ub ){
  lb_height = lb;
  ub_height = ub;
}

int rectangleEBStiffness::ownsDesignVar( const int dvNum ) const {
  return (thickness_num == dvNum || height_num == dvNum);
}

int rectangleEBStiffness::getNumDesignVars() const {
  return num_dvs;
}

int rectangleEBStiffness::getDesignVarNums( int * dvNums,
				       int * dvIndex, int dvLen ) const {
  for ( int k = 0; k < num_dvs; k++ ){
    if (*dvIndex < dvLen){
      dvNums[*dvIndex] = dv_nums[k];
      (*dvIndex)++;
    }
    else {
      return 0;
    }
  }

  return 1;
}

void rectangleEBStiffness::setDesignVars( const TacsScalar dvs[], int dvLen ){
  if (thickness_num >= 0 && thickness_num < dvLen){
    thickness = dvs[thickness_num];
  }
  if (height_num >= 0 && height_num < dvLen){
    height = dvs[height_num];
  }
}

void rectangleEBStiffness::getDesignVars( TacsScalar dvs[], int dvLen ) const {
  if (thickness_num >= 0 && thickness_num < dvLen){
    dvs[thickness_num] = thickness;
  }
  if (height_num >= 0 && height_num < dvLen){
    dvs[height_num] = height;
  }
}

void rectangleEBStiffness::getDesignVarRange( TacsScalar lb[], TacsScalar ub[],
					 int dvLen ) const {
  if (thickness_num >= 0 && thickness_num < dvLen){
    lb[thickness_num] = lb_thickness;
    ub[thickness_num] = ub_thickness;
  }
  if (height_num >= 0 && height_num < dvLen){
    lb[height_num] = lb_height;
    ub[height_num] = ub_height;
  }
}
//  if (diameter_num == dvNum || thickness_num == dvNum){

void rectangleEBStiffness::getStiffness( const double pt[],
                                    TacsScalar Ct[] ){
  memset(Ct, 0, 10*sizeof(TacsScalar));

  TacsScalar A = height*thickness;
  TacsScalar Iz = (1/12.0)*thickness*thickness*thickness*height;
  TacsScalar Iy = (1/12.0)*height*height*height*thickness;
  /*
  Torsion constant from:
  Young, Warren C., and Richard G. Budnyas. Roark’s formulas for stress and
    strain. McGraw-Hill, 2017. (pg. 401)
  */
  TacsScalar AR = thickness/height;
  TacsScalar J = A*A*AR*(1/3.0-0.21*AR*(1-AR*AR*AR*AR/12.0));

  Ct[0] = E*A;
  Ct[1] = -oy*E*A*thickness/2.0;
  Ct[2] = -oz*E*A*height/2.0;
  Ct[4] = E*(Iz+A*oy*oy*thickness*thickness/4.0);
  Ct[5] = oz*oy*E*A*A/4.0;
  Ct[7] = E*(Iy+A*oz*oz*height*height/4.0);
  Ct[9] = G*J;
}

void rectangleEBStiffness::pointwiseMass( const double pt[],
                                     TacsScalar ptmass[] ){
  memset(ptmass, 0, 6*sizeof(TacsScalar));

  TacsScalar A = height*thickness;
  TacsScalar Iz = (1/12.0)*thickness*thickness*thickness*height;
  TacsScalar Iy = (1/12.0)*height*height*height*thickness;

  ptmass[0] = rho*A;
  ptmass[3] = rho*Iz;
  ptmass[5] = rho*Iy;
}

void rectangleEBStiffness::getStiffnessDVSens( int dvNum,
                                          const double pt[],
                                          TacsScalar Ct[] ){
  memset(Ct, 0, 10*sizeof(TacsScalar));

  if (height_num == dvNum){
    TacsScalar A = height*thickness;
    TacsScalar dA = thickness;
    TacsScalar Iz = (1/12.0)*thickness*thickness*thickness*height;
    TacsScalar dIz = (1/12.0)*thickness*thickness*thickness;
    TacsScalar Iy = (1/12.0)*height*height*height*thickness;
    TacsScalar dIy = (1/4.0)*height*height*thickness;
    TacsScalar AR = thickness/height;
    TacsScalar dAR = -thickness/(height*height);
    TacsScalar J = A*A*AR*(1/3.0-0.21*AR*(1-AR*AR*AR*AR/12.0));
    TacsScalar dJ = (2*dA*A*AR+A*A*dAR)*(1/3.0-0.21*AR*(1-AR*AR*AR*AR/12.0))
                         +A*A*AR*(-0.21*(dAR-5*AR*AR*AR*AR*dAR/12.0));

    Ct[0] = E*dA;
    Ct[1] = -oy*E*dA*thickness/2.0;
    Ct[2] = -oz*E*(dA*height + A)/2.0;
    Ct[4] = E*(dIz+dA*oy*oy*thickness*thickness/4.0);
    Ct[5] = oz*oy*E*A*dA/2.0;
    Ct[7] = E*(dIy+oz*oz*(dA*height*height+A*2.0*height)/4.0);
    Ct[9] = G*dJ;
  }
  else if (thickness_num == dvNum){
    TacsScalar A = height*thickness;
    TacsScalar dA = height;
    TacsScalar Iz = (1/12.0)*thickness*thickness*thickness*height;
    TacsScalar dIz = (1/4.0)*thickness*thickness*height;
    TacsScalar Iy = (1/12.0)*height*height*height*thickness;
    TacsScalar dIy = (1/12.0)*height*height*height;
    TacsScalar AR = thickness/height;
    TacsScalar dAR = 1.0/height;
    TacsScalar J = A*A*AR*(1/3.0-0.21*AR*(1-AR*AR*AR*AR/12.0));
    TacsScalar dJ = (2*dA*A*AR+A*A*dAR)*(1/3.0-0.21*AR*(1-AR*AR*AR*AR/12.0))
                         +A*A*AR*(-0.21*(dAR-5*AR*AR*AR*AR*dAR/12.0));

    Ct[0] = E*dA;
    Ct[1] = -oy*E*(dA*thickness+A)/2.0;
    Ct[2] = -oz*E*dA*height/2.0;
    Ct[4] = E*(dIz+oy*oy*(dA*thickness*thickness+A*2.0*thickness)/4.0);
    Ct[5] = oz*oy*E*A*dA/2.0;
    Ct[7] = E*(dIy+dA*oz*oz*height*height/4.0);
    Ct[9] = G*dJ;
  }
}

void rectangleEBStiffness::pointwiseMassDVSens( int dvNum,
                                           const double pt[],
                                           TacsScalar massDVSens[] ){
  if (height_num == dvNum){
    TacsScalar dA = thickness;
    TacsScalar dIz = (1/12.0)*thickness*thickness*thickness;
    TacsScalar dIy = (1/4.0)*height*height*thickness;

    massDVSens[0] = rho*dA;
    massDVSens[3] = rho*dIz;
    massDVSens[5] = rho*dIy;
  }
  else if (thickness_num == dvNum){
    TacsScalar dA = height;
    TacsScalar dIz = (1/4.0)*thickness*thickness*height;
    TacsScalar dIy = (1/12.0)*height*height*height;

    massDVSens[0] = rho*dA;
    massDVSens[3] = rho*dIz;
    massDVSens[5] = rho*dIy;
  }
}

/*
  Calculate the stress at the outer radius
*/

void rectangleEBStiffness::calcFaceStress( TacsScalar y, TacsScalar z,
                                       TacsScalar stress[], const TacsScalar strain[] ) const {
  // sigma_x = z * k1 - y * k2

  // calculate the stress contribution from the linear strain
  stress[0] = E*(strain[0] + z*strain[1] - y*strain[2]);
  stress[1] = stress[2] = stress[3] = TacsScalar(0.0);
  stress[4] = G*y*strain[3]; // sigma_xz
  stress[5] = -G*z*strain[3]; // sigma_xy
}

/*
  The sensitivity of the failure load to the strain must be computed
  through the stresses and the vonMises calculation as well...

  df/d( epsilon_j ) = df/d( sigma_i )*dsigma_i/d ( epsilon_j )

  stressSens is the sensitivity of df/d( sigma_i )
*/

void rectangleEBStiffness::calcFaceStressSens( TacsScalar y, TacsScalar z, TacsScalar sens[],
					const TacsScalar stressSens[],
					const TacsScalar strain[] ) const {


  sens[0] = E*stressSens[0];

  sens[1] = E*z*stressSens[0];

  sens[2] = -E*y*stressSens[0];

  sens[3] = G*y*stressSens[4] - G*z*stressSens[5];
}

/*
  Calculate the stress at the outer radius
*/

void rectangleEBStiffness::calcFaceStressYZSens( TacsScalar dy, TacsScalar dz,
                                       TacsScalar stress[], const TacsScalar strain[] ) const {
  // sigma_x = z * k1 - y * k2

  // calculate the stress contribution from the linear strain
  stress[0] = E*(dz*strain[1] - dy*strain[2]);
  stress[1] = stress[2] = stress[3] = TacsScalar(0.0);
  stress[4] = G*dy*strain[3]; // sigma_xz
  stress[5] = -G*dz*strain[3]; // sigma_xy
}

void rectangleEBStiffness::failure( const double pt[],
                               const TacsScalar strain[],
                               TacsScalar * fail ){
  TacsScalar stress[6];
  TacsScalar temp[4];
  TacsScalar fail_max = 0;
  TacsScalar ks_sum = 0;
  TacsScalar y = thickness/2*(oy-1);
  TacsScalar z = height/2*(oz-1);
  calcFaceStress( y, z, stress, strain);
  temp[0] = VonMisesFailure3D(stress, yield_stress);
  fail_max = temp[0];

  y = thickness/2*(oy+1);
  calcFaceStress( y, z, stress, strain);
  temp[1] = VonMisesFailure3D(stress, yield_stress);

  y = thickness/2*(oy-1);
  z = height/2*(oz+1);
  calcFaceStress( y, z, stress, strain);
  temp[2] = VonMisesFailure3D(stress, yield_stress);
  y = thickness/2*(oy+1);
  calcFaceStress( y, z, stress, strain);
  temp[3] = VonMisesFailure3D(stress, yield_stress);
  for ( int k = 0; k < 4; k++ ){
    if(temp[k]>fail_max){fail_max=temp[k];}
  }
  for ( int k = 0; k < 4; k++ ){
    if(temp[k]>fail_max){fail_max=temp[k];}
  }

  for ( int i = 0; i < 4; i++ ){
       ks_sum += exp(ks_weight*(temp[i] - fail_max));
  }
  *fail = fail_max+log(ks_sum)/ks_weight;
}

void rectangleEBStiffness::failureStrainSens( const double pt[],
                                         const TacsScalar strain[],
					 TacsScalar sens[] ){
  TacsScalar stress[6], stressSens[6];
  TacsScalar temp[4], tempSens[4][4];
  TacsScalar fail_max = 0;
  TacsScalar ks_sum = 0;
  TacsScalar y = thickness/2*(oy-1);
  TacsScalar z = height/2*(oz-1);
       for( int j = 0; j < 4; j++){
            sens[j] = 0;
            }
  calcFaceStress( y, z, stress, strain);
  temp[0] = VonMisesFailure3D(stress, yield_stress);
  VonMisesFailure3DStressSens(stressSens, stress, yield_stress);
  fail_max = temp[0];
  calcFaceStressSens( y, z, tempSens[0], stressSens, strain);

  y = thickness/2*(oy+1);
  calcFaceStress( y, z, stress, strain);
  temp[1] = VonMisesFailure3D(stress, yield_stress);
  VonMisesFailure3DStressSens(stressSens, stress, yield_stress);
  calcFaceStressSens( y, z, tempSens[1], stressSens, strain);

  y = thickness/2*(oy-1);
  z = height/2*(oz+1);
  calcFaceStress( y, z, stress, strain);
  temp[2] = VonMisesFailure3D(stress, yield_stress);
  VonMisesFailure3DStressSens(stressSens, stress, yield_stress);
  calcFaceStressSens( y, z, tempSens[2], stressSens, strain);

  y = thickness/2*(oy+1);
  calcFaceStress( y, z, stress, strain);
  temp[3] = VonMisesFailure3D(stress, yield_stress);
  VonMisesFailure3DStressSens(stressSens, stress, yield_stress);
  calcFaceStressSens( y, z, tempSens[3], stressSens, strain);


  for ( int k = 0; k < 4; k++ ){
    if(temp[k]>fail_max){fail_max=temp[k];}
  }
  for ( int k = 0; k < 4; k++ ){
    if(temp[k]>fail_max){fail_max=temp[k];}
  }

  for ( int i = 0; i < 4; i++ ){
       ks_sum += exp(ks_weight*(temp[i] - fail_max));
  }

  for ( int i = 0; i < 4; i++ ){
       for( int j = 0; j < 4; j++){
            sens[j] += exp(ks_weight*(temp[i] - fail_max))/ks_sum*tempSens[i][j];
            }
  }

}

void rectangleEBStiffness::failureDVSens( int dvNum,
                                     const double pt[],
				     const TacsScalar strain[],
                                     TacsScalar * failSens ){

  TacsScalar stress[6], dstress[6],stressSens[6];
  TacsScalar temp[4], tempSens[4];
  TacsScalar fail_max = 0;
  TacsScalar ks_sum = 0;
  TacsScalar dheight = 0, dthickness = 0;
  if (height_num == dvNum){dheight =1.0;}
  if (thickness_num == dvNum){dthickness=1.0;}

  TacsScalar y = thickness/2*(oy-1);
  TacsScalar z = height/2*(oz-1);
  TacsScalar dy = dthickness/2*(oy-1);
  TacsScalar dz = dheight/2*(oz-1);
  calcFaceStress( y, z, stress, strain);
  temp[0] = VonMisesFailure3D(stress, yield_stress);
  VonMisesFailure3DStressSens(stressSens, stress, yield_stress);
  calcFaceStressYZSens( dy, dz, dstress, strain);
  tempSens[0] = 0;
  for ( int k = 0; k < 6; k++ ){
    tempSens[0]+=dstress[k]*stressSens[k];
  }
  fail_max = temp[0];


  y = thickness/2*(oy+1);
  dy = dthickness/2*(oy+1);
  calcFaceStress( y, z, stress, strain);
  temp[1] = VonMisesFailure3D(stress, yield_stress);
  VonMisesFailure3DStressSens(stressSens, stress, yield_stress);
  calcFaceStressYZSens( dy, dz, dstress, strain);
  tempSens[1] = 0;
  for ( int k = 0; k < 6; k++ ){
    tempSens[1]+=dstress[k]*stressSens[k];
  }

  y = thickness/2*(oy-1);
  z = height/2*(oz+1);
  dy = dthickness/2*(oy-1);
  dz = dheight/2*(oz+1);
  calcFaceStress( y, z, stress, strain);
  temp[2] = VonMisesFailure3D(stress, yield_stress);
  VonMisesFailure3DStressSens(stressSens, stress, yield_stress);
  calcFaceStressYZSens( dy, dz, dstress, strain);
  tempSens[2] = 0;
  for ( int k = 0; k < 6; k++ ){
    tempSens[2]+=dstress[k]*stressSens[k];
  }

  y = thickness/2*(oy+1);
  dy = dthickness/2*(oy+1);
  calcFaceStress( y, z, stress, strain);
  temp[3] = VonMisesFailure3D(stress, yield_stress);
  VonMisesFailure3DStressSens(stressSens, stress, yield_stress);
  calcFaceStressYZSens( dy, dz, dstress, strain);
  tempSens[3] = 0;
  for ( int k = 0; k < 6; k++ ){
    tempSens[3]+=dstress[k]*stressSens[k];
  }

  for ( int k = 0; k < 4; k++ ){
    if(temp[k]>fail_max){fail_max=temp[k];}
  }
  for ( int k = 0; k < 4; k++ ){
    if(temp[k]>fail_max){fail_max=temp[k];}
  }

  for ( int i = 0; i < 4; i++ ){
       ks_sum += exp(ks_weight*(temp[i] - fail_max));
  }
  *failSens = 0;
       for( int j = 0; j < 4; j++){
            *failSens += exp(ks_weight*(temp[j] - fail_max))/ks_sum*tempSens[j];
            }
}

void rectangleEBStiffness::printInfo(){
  double gpt[2] = {0.0, 0.0};


  TacsScalar mass[6];
  pointwiseMass(gpt, mass);

  TacsScalar A = height*thickness;
  TacsScalar Iz = (1/12.0)*height*height*height*thickness;
  TacsScalar Iy = (1/12.0)*thickness*thickness*thickness*height;
  TacsScalar AR = thickness/height;
  TacsScalar J = A*A*AR*(1/3.0-0.21*AR*(1-AR*AR*AR*AR/12.0));

  printf("%s properties (evaluated at centroid)\n", constName);
  printf("beam height:    %15.6f\n", TacsRealPart(height));
  printf("beam thickness: %15.6f\n", TacsRealPart(thickness));
  printf("mass per unit length:  %15.6f\n", TacsRealPart(mass[0]));
  printf("Iz:        %15.6E\n", TacsRealPart(Iz));
  printf("Iy:       %15.6E\n", TacsRealPart(Iy));
  printf("J:         %15.6E\n", TacsRealPart(J));
  printf("A:        %15.6E\n", TacsRealPart(A));
}
