#ifndef EB_BEAM_STIFFNESS
#define EB_BEAM_STIFFNESS

/*
  The stiffness object for the variety of different beams

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

#include "TACSConstitutive.h"
enum EBBeamReferenceDirection { STRONG_AXIS, WEAK_AXIS };

class EBStiffness : public TACSConstitutive {
 public:
  static const int NUM_STRESSES = 4;

  EBStiffness( TacsScalar rho, TacsScalar E, TacsScalar G,
	       TacsScalar A, TacsScalar Ix, TacsScalar Iy, TacsScalar J,
	       TacsScalar _ref_dir[3],
	       enum EBBeamReferenceDirection _ref_type=WEAK_AXIS );

  virtual ~EBStiffness(){}

  // Calculate the stress
  // --------------------
  int getNumStresses();
  void calculateStress( const double gpt[], const TacsScalar strain[],
			TacsScalar stress[] );
  void calculateStressDVSens( int dvNum,
                              const double gpt[], const TacsScalar strain[],
			      TacsScalar stress[] );

  // Evaluate the stiffness matrix
  // -----------------------------
  virtual void getStiffness( const double pt[], TacsScalar Ct[] );
  virtual void getStiffnessDVSens( int dvNum, const double pt[],
                                   TacsScalar Ct[] );

  // Return the mass moments
  // -----------------------
  int getNumMassMoments(){ return 6; }
  void getPointwiseMass( const double pt[], TacsScalar mass[] );
  void pointwiseMassDVSens( int dvNum, const double gpt[],
                            TacsScalar massDVSens[] );

  const char *constitutiveName();

  int ownsDesignVar( const int dvNum ) const;
  int getNumDesignVars() const;
  int getDesignVarNums( int * dvNums, int * dvIndex, int dvLen ) const;

  // Copy of the input parameters. Necessary for BDF output.
  TacsScalar rho, E, G, A, Ix, Iy, J, ref_dir[3];
  TacsScalar fuck;
  enum EBBeamReferenceDirection ref_type;

 protected:
  EBStiffness();

  TacsScalar mass[6];
  TacsScalar C[10];

  // Calculate the stress resultants
  inline void calcStress( const TacsScalar Ct[], const TacsScalar e[],
                          TacsScalar s[] );
 private:
  static const char * constName;
};

/*
  The constitutive matrix is,

  [ Ct[0] Ct[1] Ct[2] Ct[3] ]
  [ Ct[1] Ct[4] Ct[5] Ct[6] ]
  [ Ct[2] Ct[5] Ct[7] Ct[8] ]
  [ Ct[3] Ct[6] Ct[8] Ct[9] ]
*/
inline void EBStiffness::calcStress( const TacsScalar Ct[],
                                     const TacsScalar e[],
				     TacsScalar s[] ){
  s[0] = Ct[0]*e[0] + Ct[1]*e[1] + Ct[2]*e[2] + Ct[3]*e[3];
  s[1] = Ct[1]*e[0] + Ct[4]*e[1] + Ct[5]*e[2] + Ct[6]*e[3];
  s[2] = Ct[2]*e[0] + Ct[5]*e[1] + Ct[7]*e[2] + Ct[8]*e[3];
  s[3] = Ct[3]*e[0] + Ct[6]*e[1] + Ct[8]*e[2] + Ct[9]*e[3];
}

#endif
