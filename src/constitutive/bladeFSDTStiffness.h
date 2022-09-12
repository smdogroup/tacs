#ifndef TACS_BLADE_FSDT_STIFFNESS_H
#define TACS_BLADE_FSDT_STIFFNESS_H

/*!
  Copyright (c) 2011 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

#include "FSDTStiffness.h"
#include "MaterialProperties.h"

/*
  A blade-stiffened panel

  A class for calculating the stiffness and failure properties of
  blade-stiffened panels. The panel stiffness characteristics are
  calculated based on smeared though-thicknes properties.

  The following buckling modes are included:

  1. Global axial/shear buckling of the entire panel
  2. Axial and shear buckling of the skin between stiffeners
  3. Axial and shear buckling of the blade-stiffener

  The following failure modes are included

  1. Failure of the skin via Tsai--Wu failure envelope
  - Calculate the strain at the outer-plies:
  \epsilon_{0} +/- t/2*\kappa

  Design variables:

  Ly/Lx = the approximate dimensions of the panel

  sp = stiffener pitch (spacing) = Ly/ny_stiffener (ny_stiffner - the
  number of stiffeners in the transverse direction is not an integer
  although it should be in reality.)

  ts = skin thickness

  sh = stiffener height

  ply_angles = ply angles included in the optimization

  ply_fractions = ply percentages or ply fractions - the relative
  number of plies at a given orientation
*/
class bladeFSDTStiffness : public FSDTStiffness {
 public:
  static const int NUM_PLIES = 4;

  bladeFSDTStiffness(OrthoPly* _ortho_ply, TacsScalar _kcorr, TacsScalar _Lx,
                     int _Lx_num, TacsScalar _sp, int _sp_num, TacsScalar _sh,
                     int _sh_num, TacsScalar _st, int _st_num, TacsScalar _t,
                     int _t_num, int _pf_nums[], int _stiff_pf_nums[]);
  ~bladeFSDTStiffness();

  // Set non-default values
  // ----------------------
  void setStiffenerPitchBounds(TacsScalar _sp_low, TacsScalar _sp_high);
  void setStiffenerHeightBounds(TacsScalar _sh_low, TacsScalar _sh_high);
  void setStiffenerThicknessBounds(TacsScalar _st_low, TacsScalar _st_high);
  void setStiffenerPlyFractionBounds(TacsScalar _pf_low[],
                                     TacsScalar _pf_high[]);
  void setStiffenerPlyFractions(TacsScalar _pf[]);
  void setThicknessBounds(TacsScalar _t_low, TacsScalar _t_high);
  void setPlyFractionBounds(TacsScalar _pf_low[], TacsScalar _pf_high[]);
  void setPlyFractions(TacsScalar _pf[]);

  void setKSWeight(double _ks_weight);

  // Retrieve the critical loads
  // ---------------------------
  void getCriticalGlobalLoads(TacsScalar* Nx, TacsScalar* Nxy);
  void getCriticalPanelLoads(TacsScalar* Nx, TacsScalar* Nxy);
  void getCriticalStiffenerLoads(TacsScalar* Nx, TacsScalar* Nxy);

  // Compute the critical shear loads and the sensitivities
  // ------------------------------------------------------
  TacsScalar calcCriticalShearLoad(TacsScalar D1, TacsScalar D2, TacsScalar D3,
                                   TacsScalar L);
  TacsScalar calcCriticalShearLoadSens(TacsScalar D1, TacsScalar D2,
                                       TacsScalar D3, TacsScalar L,
                                       TacsScalar sD1, TacsScalar sD2,
                                       TacsScalar sD3, TacsScalar sL);

  // Retrieve the local loads
  // ------------------------
  void getPanelLoads(const TacsScalar stress[], TacsScalar* Nx,
                     TacsScalar* Nxy);
  void getStiffenerLoads(const TacsScalar stress[], TacsScalar* Nx,
                         TacsScalar* Nxy);

  // Functions for panel-only optimizations
  // --------------------------------------
  TacsScalar getFailureCriterion(const TacsScalar stress[]);
  TacsScalar getBucklingCriterion(const TacsScalar stress[]);

  // Additional constraint evaluation functions
  // ------------------------------------------
  int isLinear();
  int getNumCon();
  int getConCSRSize();
  int getConRange(int offset, TacsScalar lb[], TacsScalar ub[]);
  int addConCSR(int offset, int rowp[], int cols[]);
  int evalCon(int offset, TacsScalar con[]);
  int evalConDVSens(int offset, TacsScalar Acol[], const int rowp[],
                    const int cols[]);

  // Design variable information set up
  // ----------------------------------
  int ownsDesignVar(const int dvNum) const;
  void setDesignVars(const TacsScalar dvs[], int numDVs);
  void getDesignVars(TacsScalar dvs[], int numDVs);
  void getDesignVarRange(TacsScalar lowerBound[], TacsScalar upperBound[],
                         int numDVs);
  int getLxNum();

  int getNumDesignVars();
  int getDesignVarNums(int dvNums[], int dvLen);

  // Functions required by FSDTStiffness
  // -----------------------------------
  void getPointwiseMass(const double gpt[], TacsScalar mass[]);
  void pointwiseMassDVSens(int dvNum, const double gpt[], TacsScalar mass[]);
  void addPointwiseMassDVSens(const double pt[], const TacsScalar alpha[],
                              TacsScalar dvSens[], int dvLen);
  TacsScalar getDensity();
  // Get the stiffness matrices and their derivatives
  // ------------------------------------------------
  TacsScalar getStiffness(const double gpt[], TacsScalar A[], TacsScalar B[],
                          TacsScalar D[], TacsScalar As[]);
  TacsScalar getStiffnessDVSens(int dvNum, const double gpt[], TacsScalar sA[],
                                TacsScalar sB[], TacsScalar sD[],
                                TacsScalar sAs[]);
  void addStiffnessDVSens(const double pt[], const TacsScalar e[],
                          const TacsScalar psi[], TacsScalar rotPsi,
                          TacsScalar fdvSens[], int dvLen);
  TacsScalar calculateStressDVSens(const double pt[], const TacsScalar strain[],
                                   TacsScalar s_sens[], int dvNum);

  // Calculate the point-wise failure criteria
  // -----------------------------------------
  void failure(const double gpt[], const TacsScalar strain[], TacsScalar* fail);
  void failureStrainSens(const double gpt[], const TacsScalar strain[],
                         TacsScalar sens[]);
  void failureDVSens(int dvNum, const double gpt[], const TacsScalar strain[],
                     TacsScalar* failSens);
  void addFailureDVSens(const double pt[], const TacsScalar strain[],
                        TacsScalar alpha, TacsScalar dvSens[], int dvLen);

  // Calculate the buckling criteria
  // -------------------------------
  void buckling(const TacsScalar strain[], TacsScalar* bval);
  void bucklingStrainSens(const TacsScalar strain[], TacsScalar* bvalSens);
  void bucklingDVSens(int dvNum, const TacsScalar strain[],
                      TacsScalar* bvalSens);
  void addBucklingDVSens(const TacsScalar strain[], TacsScalar alpha,
                         TacsScalar dvSens[], int dvLen);

  const char* constitutiveName() const;

  // Retrieve the design variable for plotting purposes
  // --------------------------------------------------
  TacsScalar getDVOutputValue(int dv_index, const double gpt[]);

  // Write the buckling properties to stdout
  // ---------------------------------------
  void printBucklingInfo();

 private:
  // Update the stiffness information A, B, D, As, k_penalty
  void updateStiffness();

  // Calculate the buckling loads
  void calcBucklingProps();

  // Calculate the derivative of the buckling loads
  void calcBucklingPropsDVSens(int dvNum, TacsScalar* sNx_global_crit,
                               TacsScalar* sNxy_global_crit,
                               TacsScalar* sNx_panel_crit,
                               TacsScalar* sNxy_panel_crit,
                               TacsScalar* sNx_weak_stiff_crit,
                               TacsScalar* sNx_strong_stiff_crit);

  // Calculate the failure criteria
  int calcFail(const TacsScalar strain[], TacsScalar fail[], TacsScalar* _max);
  void calcFailStrainSens(const TacsScalar strain[], const TacsScalar weights[],
                          TacsScalar sens[]);
  TacsScalar calcFailDVSens(int dvNum, const TacsScalar strain[],
                            const TacsScalar weights[]);

  // Calculate the buckling criteria
  int calcBuckling(const TacsScalar strain[], TacsScalar bval[],
                   TacsScalar* _max);
  void calcBucklingStrainSens(const TacsScalar strain[],
                              const TacsScalar weights[], TacsScalar sens[]);
  TacsScalar calcBucklingDVSens(int dvNum, const TacsScalar strain[],
                                const TacsScalar weights[]);

  static const int NUM_BUCKLING = 4;
  static const int MAX_NUM_FAIL = 3 * NUM_PLIES;

  // Stiffener information
  // ---------------------
  // Stiffener pitch variables
  int sp_num;
  TacsScalar sp, sp_low, sp_high;

  // Stiffener height
  int sh_num;
  TacsScalar sh, sh_low, sh_high;

  // Stiffener thickness
  int st_num;
  TacsScalar st, st_low, st_high;

  // Flage fraction - the fraction of the length of flange over the
  // stiffener height
  TacsScalar sf_fraction;

  int num_pf_con;  // Number of ply-fraction constraints: 0, 1 or 2

  // Ply percentage information
  // --------------------------
  int num_stiff_pf_unique;
  int stiff_pf_nums[NUM_PLIES], stiff_pf_unique[NUM_PLIES];
  TacsScalar stiff_pf[NUM_PLIES], stiff_ply_angles[NUM_PLIES],
      stiff_pf_low[NUM_PLIES], stiff_pf_high[NUM_PLIES];

  // Panel information
  // -----------------
  // Panel thickness
  int t_num;
  TacsScalar t, t_low, t_high;

  // Ply percentage information
  int num_pf_unique;
  int pf_nums[NUM_PLIES], pf_unique[NUM_PLIES];
  TacsScalar pf[NUM_PLIES], ply_angles[NUM_PLIES], pf_low[NUM_PLIES],
      pf_high[NUM_PLIES];

  // Panel dimension information
  TacsScalar Lx;
  int Lx_num;

  // Buckling information
  // --------------------
  // Critical buckling loads
  TacsScalar Nx_global_crit, Nxy_global_crit;
  TacsScalar Nx_panel_crit, Nxy_panel_crit;
  TacsScalar Nx_weak_stiff_crit, Nx_strong_stiff_crit;

  // Stiffness information
  // ---------------------
  // Stiffness information for the panel and the blade stiffener
  TacsScalar Abar_panel[3], Qbar_panel[6];
  TacsScalar Abar_stiff[3], Qbar_stiff[6];

  // The stiffness matrices
  TacsScalar A[6], B[6], D[6], As[6];
  TacsScalar k_penalty;

  // Information about the stiffeners
  OrthoPly* ortho_ply;

  TacsScalar kcorr;  // The shear correction factor
  double ks_weight;  // The ks weight for the failure envelope calculations

  static const char* constName;
};

#endif
