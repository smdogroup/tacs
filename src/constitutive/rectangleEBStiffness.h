#ifndef RECTANGLE_EB_BEAM_STIFFNESS_H
#define RECTANGLE_EB_BEAM_STIFFNESS_H

/*
  The stiffness object for the variety of different beams

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

#include "EBStiffness.h"

/*
  This is a stiffness object for the EB (Euler-Bernoulli) beam class.

  These expressions represent a rectangular beam with a specified height
  and thickness.
*/

class rectangleEBStiffness : public EBStiffness {
 public:
  rectangleEBStiffness( TacsScalar _rho, TacsScalar _E, TacsScalar _G,
                                  TacsScalar _yield_stress,
                                  TacsScalar _height,
                                  TacsScalar _thickness,
                                  int _height_num, int _thickness_num,
			                            TacsScalar _ref_dir[3],
				                          TacsScalar z_offset=0, TacsScalar y_offset=0,
                                  enum EBBeamReferenceDirection _ref_type=WEAK_AXIS );
  ~rectangleEBStiffness(){}

  // Set the thickness and height bounds
  // -------------------------------------
  void setThicknessBounds( TacsScalar lb, TacsScalar ub );
  void setHeightBounds( TacsScalar lb, TacsScalar ub );

  // Functions for design variable control
  // -------------------------------------
  const char * constitutiveName() const;

  int ownsDesignVar( const int dvNum ) const;
  int getNumDesignVars() const;
  int getDesignVarNums( int * dvNums, int * numDVs, int dvLen ) const ;
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs ) const;
  void getDesignVarRange( TacsScalar lowerBound[],
			  TacsScalar upperBound[], int numDVs ) const;

  // Functions required by FSDTStiffness
  // -----------------------------------
  void getStiffness( const double gpt[], TacsScalar Ct[] );
  void getStiffnessDVSens( int dvNum, const double gpt[],
                           TacsScalar Ct[] );

  // Get the mass properites
  // -----------------------
  void pointwiseMass( const double gpt[], TacsScalar mass[] );
  void pointwiseMassDVSens( int dvNum, const double gpt[],
                            TacsScalar massDVSens[] );

  // Functions to compute the failure properties
  // -------------------------------------------
  void failure( const double pt[], const TacsScalar strain[],
                TacsScalar * fail );
  void failureStrainSens( const double pt[], const TacsScalar strain[],
                          TacsScalar sens[] );
  void failureDVSens( int dvNum, const double pt[],
                      const TacsScalar strain[], TacsScalar * failSens );

  void printInfo();

  // inline double RealPart(const double& r) {
  // **
  //  *  So the RealPart() statement can be used even with
  //  *  the double version of the code to be complexified.
  //  *  Most useful inside printf statements.
  //  **
  // return r;
  // }
 private:
  // Calculate the stress at the outer radius
  void calcFaceStress( TacsScalar y, TacsScalar z,
                       TacsScalar stress[], const TacsScalar strain[] ) const;
  void calcFaceStressSens( TacsScalar y, TacsScalar z, TacsScalar sens[],
		     const TacsScalar stressSens[], const TacsScalar strain[] ) const;
  void calcFaceStressYZSens( TacsScalar dy, TacsScalar dz,
                             TacsScalar stress[], const TacsScalar strain[] ) const;


  // The material properties
  // -----------------------
  TacsScalar E, G, rho;
  TacsScalar yield_stress;
  TacsScalar height;
  TacsScalar thickness;
  TacsScalar oz, oy;
  TacsScalar ks_weight;

  TacsScalar lb_thickness, ub_thickness;
  TacsScalar lb_height, ub_height;

  // Design variable information
  // ---------------------------
  int thickness_num, height_num;
  int dv_nums[2], num_dvs;
  static const char * constName;
};

#endif
