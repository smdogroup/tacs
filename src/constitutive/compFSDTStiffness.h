#ifndef COMP_FSDT_STIFFNESS_H
#define COMP_FSDT_STIFFNESS_H

/*
  Copyright (c) 2011 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include "FSDTStiffness.h"
#include "MaterialProperties.h"

/*
  This class inherits from the FSDTStiffness class. It is not designed
  for gradient-based optimization but is very useful for testing
  purposes or for generating response surface data.

  The stiffness properties are defined by passing in the OrthoPly
  class with given ply angles. This class calculates the stiffness
  based on those inputs. The failure properties are calculated based
  upon the most recent input as well. In this manner, the number of
  plies can vary from call to call - without re-setting the material
  class within the shell elements.
*/
class compFSDTStiffness : public FSDTStiffness {
 public:
  compFSDTStiffness( OrthoPly **_ortho_ply, TacsScalar _kcorr,
		     TacsScalar *_thickness, TacsScalar *_ply_angles, 
		     int _num_plies );
  ~compFSDTStiffness();
 
  const char *constitutiveName();

  // Functions required by FSDTStiffness
  // -----------------------------------
  TacsScalar getStiffness( const double pt[],
                           TacsScalar A[], TacsScalar B[],
                           TacsScalar D[], TacsScalar As[] );
  void getPointwiseMass( const double pt[], TacsScalar mass[] );

  // Compute the failure criteria
  // ----------------------------
  void failure( const double pt[], const TacsScalar strain[], 
                TacsScalar *fail );

 private:
  // Get the strain in a particular lamina -- still in the global ref. axis
  void getLaminaStrain( TacsScalar strain[], 
			const TacsScalar rmStrain[], TacsScalar tp );

  TacsScalar evalFailLoads( const TacsScalar constStrain[],
                            const TacsScalar linStrain[] );

  // The composite data
  // ------------------
  OrthoPly **ortho_ply;
  TacsScalar kcorr;
  TacsScalar *thickness;
  TacsScalar *ply_angles;
  int num_plies;

  static const char *constName;
};

#endif
