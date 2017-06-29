#include "SolidStiffness.h"

const char * SolidStiffness::constName = "SolidStiffness";

const char * SolidStiffness::constitutiveName(){
  return constName;
}
  
/*
  SolidStiffness member function definitions: Definitions for a simple
  beam element
*/

SolidStiffness::SolidStiffness( TacsScalar _rho, TacsScalar E, 
                                TacsScalar nu ){
  rho = _rho;

  TacsScalar D = E/((1.0 + nu)*(1.0 - 2.0*nu));
  C[0] = C[3] = C[5] = (1.0 - nu)*D;
  C[1] = C[2] = C[4] = nu*D;

  G23 = G13 = G12 = 0.5*E/(1.0 + nu);
}

/*!
  Calculate the stiffness matrix based upon the material properties.

  This material is from Jones (Mechanics of Composite Materials pg. 38)

  The flexibility form in the orthotropic axises is,

  [ e_1 ]   [       1/E1, - nu_21/E2, - nu_31/E3 ][ s_1 ]
  [ e_2 ] = [ - nu_12/E1,       1/E2, - nu_32/E3 ][ s_2 ]
  [ e_3 ]   [ - nu_13/E1, - nu_23/E2,       1/E3 ][ s_3 ]

  The Poisson's ratio is defined as follows,

  nu_ij = - e_j/e_i

  Due to the symmetry of the matrix

  nu_21/E2 = nu_12/E1
  nu_13/E1 = nu_31/E3
  nu_23/E2 = nu_32/E3

  The stiffness matrix is obtained by inverting the matrix above. The 
  determinant of the flexibility form of the constitutive matrix is,

  D = ( (1 - nu_32 * nu_23 )  
  + nu_21 * ( - nu_12 - nu_32 * nu_13 ) 
  - nu_31 * ( nu_12*nu_23 + nu_13 ) )/(E1*E2*E3)

  D = ( 1 
  - nu_32*nu_23 - nu_21*nu_12 - nu_31*nu_13 
  - nu_21*nu_32*nu_13 - nu_31*nu_12*nu_23 )/(E1*E2*E3)

  Since nu_21*nu_32*nu_13 = nu_31*nu_12*nu_23

  D = ( 1 
  - nu_32*nu_23 - nu_21*nu_12 - nu_31*nu_13 
  - 2*nu_21*nu_32*nu_13 )/(E1*E2*E3)

  The inverse may be determined by applying Cramer's rule,

  [ s_1 ]   [ C_11, C_12, C_13 ][ e_1 ]
  [ s_2 ] = [ C_21, C_22, C_23 ][ e_2 ]
  [ s_3 ]   [ C_31, C_32, C_33 ][ e_3 ]
  
  Due to symmetry,

  C_21 = C_12
  C_31 = C_13
  C_32 = C_23

  The diagonal components are,

  C_11 = (1 - nu_32*nu_23)/(E2*E3*D)
  C_22 = (1 - nu_31*nu_13)/(E1*E3*D)
  C_33 = (1 - nu_12*nu_21)/(E1*E2*D)

  The off-diagonal components are,

  C_21 = C_12 = 
  |      1/E1,  1, - nu_31/E3 |
  |- nu_12/E1,  0, - nu_32/E3 |
  |- nu_13/E1,  0,       1/E3 |/D

  C_21 = C_12 = (nu_12 + nu_32*nu_13)/(E1*E3*D)

  ---------------------------------------------
 
  C_31 = C_13 = 
  |       1/E1, - nu_21/E2, 1 |
  | - nu_12/E1,       1/E2, 0 |
  | - nu_13/E1, - nu_23/E2, 0 |/D

  C_31 = C_13 = (nu_13 + nu_12*nu_23)/(E1*E2*D)

  ----------------------------------------------

  C_32 = C_23 = 
  |       1/E1, - nu_21/E2, 0 |
  | - nu_12/E1,       1/E2, 1 |
  | - nu_13/E1, - nu_23/E2, 0 |/D

  C_32 = C_23 = (nu_23 + nu_21*nu_13)/(E1*E2*D)
*/

SolidStiffness::SolidStiffness( TacsScalar _rho, 
				TacsScalar E1, TacsScalar E2, TacsScalar E3, 
				TacsScalar nu_12, TacsScalar nu_13, 
				TacsScalar nu_23, TacsScalar _G23, 
				TacsScalar _G13, TacsScalar _G12 ){
  rho = _rho;

  TacsScalar nu_21 = (E2*nu_12)/E1;
  TacsScalar nu_31 = (E3*nu_13)/E1;
  TacsScalar nu_32 = (E3*nu_23)/E2;

  TacsScalar D = (1.0 - nu_32*nu_23 - nu_21*nu_12 - nu_31*nu_13 
                  - 2.0*nu_21*nu_32*nu_13)/(E1*E2*E3);

  C[0] = (1.0 - nu_32*nu_23)/(E2*E3*D); // C_11
  C[3] = (1.0 - nu_31*nu_13)/(E1*E3*D); // C_22
  C[5] = (1.0 - nu_12*nu_21)/(E1*E2*D); // C_33

  C[1] = (nu_12 + nu_32*nu_13)/(E1*E3*D); // C_12
  C[2] = (nu_13 + nu_12*nu_23)/(E1*E2*D); // C_13
  C[4] = (nu_23 + nu_21*nu_13)/(E1*E2*D); // C_23

  G23 = _G23;
  G13 = _G13;
  G12 = _G12;
}

SolidStiffness::SolidStiffness(){
  C[0] = C[1] = C[2] = C[3] = C[4] = C[5] = 0.0;
  G23 = G13 = G12 = 0.0;
  rho = 0.0;
}

int SolidStiffness::getNumStresses(){ return NUM_STRESSES; }

void SolidStiffness::calculateStress( const double gpt[], 
				      const TacsScalar strain[],
				      TacsScalar stress[] ){
  calcStress(strain, stress);
}
