#ifndef TACS_PLANE_STRESS_BSPLINE_STIFFNESS_H
#define TACS_PLANE_STRESS_BSPLINE_STIFFNESS_H

/*
  Copyright (c) 2016 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/
#include "TACSConstitutive.h"
#include "PlaneStressStiffness.h"
/*
  Constitutive class for plane stress problems using B-spline 
  as filter
*/
class PlaneStressBsplineStiffness : public PlaneStressStiffness {
 public:
  static const int NUM_STRESSES = PlaneStressStiffness::NUM_STRESSES;
  PlaneStressBsplineStiffness( TacsScalar _rho, 
                               TacsScalar _E, TacsScalar _nu,
                               TacsScalar _ys, TacsScalar _epsilon,
                               TacsScalar *_x, double _qpenalty,
                               double _lower,
                               double *_Tu, double *_Tv,
                               int _Lu, int _Lv, 
                               int _pNum, int _order );
  ~PlaneStressBsplineStiffness();
  // Functions for design variable control
  // -------------------------------------
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );
  void getDesignVarRange( TacsScalar lb[], 
			  TacsScalar ub[], int numDVs );
  // Calculate the stress
  // --------------------
  int getNumStresses();
  void calculateStress( const double pt[], 
                        const TacsScalar strain[],
			TacsScalar stress[] );
  void addStressDVSens( const double pt[], const TacsScalar strain[], 
                        TacsScalar alpha, const TacsScalar psi[], 
                        TacsScalar dvSens[], int dvLen );
  // Return the mass moments
  // -----------------------
  void getPointwiseMass( const double pt[], 
                         TacsScalar mass[] );
  void addPointwiseMassDVSens( const double pt[], 
                               const TacsScalar alpha[],
                               TacsScalar dvSens[], int dvLen );
  // Return the failure
  // -------------------
  void failure( const double pt[], 
                const TacsScalar strain[],
                TacsScalar * fail );
  void addFailureDVSens( const double pt[], 
                         const TacsScalar strain[],
                         TacsScalar alpha,
                         TacsScalar dvSens[], int dvLen );
  
  // Extra info about the constitutive class
  // ---------------------------------------
  const char *constitutiveName();
  
  // Misc helper functions
  // ---------------------
  int findPatch();
  void getShapeFunctions( const double pt[], 
                          double N[] );
  void computeIndexList( int **index,
                         int numDVs );
  
 private:
  static const char * constName;
  // Knot vectors and its associated variables
  double *Tu, *Tv;
  int Lu, Lv, pNum, t1, t2;
  int order, *index;
  // The constitutive properties
  double E, nu, ys, rho, epsilon;
  TacsScalar *x, xw;
  TacsScalar D11, D12, G;
  int dvNum;
  // The topology variables
  double q, lowerBound;
};
#endif //TACS_PLANE_STRESS_BSPLINE_STIFFNESS_H
