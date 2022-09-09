#include "bladeFSDTStiffness.h"

#include "FElibrary.h"
#include "tacslapack.h"

/*
  Copyright (c) 2011 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

/*
  The bladeFSDTStiffness object

  This is the constitutive class for a blade-stiffened panel. The
  contribution of the stiffeners to the panel thickness are handled
  using a smeared-stiffener approach - see Brush and Almroth: Buckling
  of bars, plates and shells.

  This method requires the initial values and the design variable
  numbers for:

  sp: the stiffener spacing
  sh: the stiffener height
  st: the stiffener thickness
  t: the panel thickness
  pf_nums: the ply fraction numbers
  stiff_pf_nums: the ply fraction numbers for the stiffeners

  The failure calculations are based on a failure envelope that
  includes both material failure and buckling constraints. The
  buckling constraints can only be used when a "PanelDimensions"
  object is passed into the object. This is required in order to
  obtain the effective panel dimensions.
*/
bladeFSDTStiffness::bladeFSDTStiffness(OrthoPly* _ortho_ply, TacsScalar _kcorr,
                                       TacsScalar _Lx, int _Lx_num,
                                       TacsScalar _sp, int _sp_num,
                                       TacsScalar _sh, int _sh_num,
                                       TacsScalar _st, int _st_num,
                                       TacsScalar _t, int _t_num,
                                       int _pf_nums[], int _stiff_pf_nums[]) {
  ortho_ply = _ortho_ply;
  ortho_ply->incref();
  kcorr = _kcorr;

  // Record the panel length and design variable number
  Lx = _Lx;
  Lx_num = _Lx_num;

  // Record the stiffener information
  // Note that the default low/high bounds can be changed
  sp = _sp;
  sp_num = _sp_num;
  sp_low = 0.1;   // 10 cm
  sp_high = 0.2;  // 20 cm

  sh = _sh;
  sh_num = _sh_num;
  sh_low = 0.04;   // 4cm
  sh_high = 0.05;  // 5cm

  st = _st;
  st_num = _st_num;
  st_low = 0.002;  // 2 mm
  st_high = 0.01;  // 1 cm

  sf_fraction = 1.0;  // Flange length fraction

  for (int k = 0; k < NUM_PLIES; k++) {
    stiff_ply_angles[k] = -0.5 * M_PI + ((k + 1) * M_PI) / NUM_PLIES;
  }

  for (int k = 0; k < NUM_PLIES; k++) {
    stiff_pf[k] = 1.0 / NUM_PLIES;
    stiff_pf_nums[k] = _stiff_pf_nums[k];
    stiff_pf_low[k] = 0.15;   // 15 %
    stiff_pf_high[k] = 0.70;  // 70 %
  }

  // Record the panel thickness information
  t = _t;
  t_num = _t_num;
  t_low = 0.002;  // 2 mm
  t_high = 0.01;  // 1 cm

  for (int k = 0; k < NUM_PLIES; k++) {
    ply_angles[k] = -0.5 * M_PI + ((k + 1) * M_PI) / NUM_PLIES;
  }

  for (int k = 0; k < NUM_PLIES; k++) {
    pf[k] = 1.0 / NUM_PLIES;
    pf_nums[k] = _pf_nums[k];
    pf_low[k] = 0.15;   // 15 %
    pf_high[k] = 0.70;  // 70 %
  }

  if (sp_num < 0) {
    fprintf(stderr, "bladeFSDTStiffness: sp_num must be non-negative\n");
  }
  if (sh_num < 0) {
    fprintf(stderr, "bladeFSDTStiffness: sh_num must be non-negative\n");
  }
  if (st_num < 0) {
    fprintf(stderr, "bladeFSDTStiffness: st_num must be non-negative\n");
  }
  if (t_num < 0) {
    fprintf(stderr, "bladeFSDTStiffness: t_num must be non-negative\n");
  }

  // Calculate additional information about the design variables
  for (int k = 0; k < NUM_PLIES; k++) {
    pf_unique[k] = pf_nums[k];
    stiff_pf_unique[k] = stiff_pf_nums[k];
  }

  num_pf_unique = FElibrary::uniqueSort(pf_unique, NUM_PLIES);
  num_stiff_pf_unique = FElibrary::uniqueSort(stiff_pf_unique, NUM_PLIES);

  num_pf_con = 0;
  if (num_pf_unique > 0) {
    num_pf_con++;
  }
  if (num_stiff_pf_unique > 0) {
    num_pf_con++;
  }

  // Initialize the stiffness and failure properties
  ks_weight = 80.0;  // Default value of the ks weight parameter
  updateStiffness();
  calcBucklingProps();
}

bladeFSDTStiffness::~bladeFSDTStiffness() { ortho_ply->decref(); }

const char* bladeFSDTStiffness::constName = "bladeFSDTStiffness";

/*
  Set non-default values into the class
*/

void bladeFSDTStiffness::setStiffenerPitchBounds(TacsScalar _sp_low,
                                                 TacsScalar _sp_high) {
  sp_low = _sp_low;
  sp_high = _sp_high;
}

void bladeFSDTStiffness::setStiffenerHeightBounds(TacsScalar _sh_low,
                                                  TacsScalar _sh_high) {
  sh_low = _sh_low;
  sh_high = _sh_high;
}

void bladeFSDTStiffness::setStiffenerThicknessBounds(TacsScalar _st_low,
                                                     TacsScalar _st_high) {
  st_low = _st_low;
  st_high = _st_high;
}

void bladeFSDTStiffness::setStiffenerPlyFractionBounds(TacsScalar _pf_low[],
                                                       TacsScalar _pf_high[]) {
  for (int k = 0; k < NUM_PLIES; k++) {
    pf_low[k] = _pf_low[k];
    pf_high[k] = _pf_high[k];
  }
}

void bladeFSDTStiffness::setStiffenerPlyFractions(TacsScalar _pf[]) {
  for (int k = 0; k < NUM_PLIES; k++) {
    stiff_pf[k] = _pf[k];
  }
}

void bladeFSDTStiffness::setThicknessBounds(TacsScalar _t_low,
                                            TacsScalar _t_high) {
  t_low = _t_low;
  t_high = _t_high;
}

void bladeFSDTStiffness::setPlyFractionBounds(TacsScalar _pf_low[],
                                              TacsScalar _pf_high[]) {
  for (int k = 0; k < NUM_PLIES; k++) {
    stiff_pf_low[k] = _pf_low[k];
    stiff_pf_high[k] = _pf_high[k];
  }
}

void bladeFSDTStiffness::setPlyFractions(TacsScalar _pf[]) {
  for (int k = 0; k < NUM_PLIES; k++) {
    pf[k] = _pf[k];
  }
}

void bladeFSDTStiffness::setKSWeight(double _ks_weight) {
  ks_weight = _ks_weight;
}

/*
  Retrieve the critical panel loads for the given set of design variables
*/
void bladeFSDTStiffness::getCriticalGlobalLoads(TacsScalar* Nx,
                                                TacsScalar* Nxy) {
  *Nx = Nx_global_crit;
  *Nxy = Nxy_global_crit;
}

void bladeFSDTStiffness::getCriticalPanelLoads(TacsScalar* Nx,
                                               TacsScalar* Nxy) {
  *Nx = Nx_panel_crit;
  *Nxy = Nxy_panel_crit;
}

void bladeFSDTStiffness::getCriticalStiffenerLoads(TacsScalar* Nx_weak,
                                                   TacsScalar* Nx_strong) {
  *Nx_weak = Nx_weak_stiff_crit;
  *Nx_strong = Nx_strong_stiff_crit;
}

/*
  Retrieve the failure criteria for panel-only optimization
  purposes (without the use of a structural analysis in TACS).
*/

TacsScalar bladeFSDTStiffness::getFailureCriterion(const TacsScalar stress[8]) {
  TacsScalar strain[8];
  memcpy(strain, stress, 8 * sizeof(TacsScalar));

  double gpt[2] = {0.0, 0.0};
  TacsScalar A[6], B[6], D[6], As[3];
  getStiffness(gpt, A, B, D, As);

  int ipiv[8];
  TacsScalar C[8 * 8];
  memset(C, 0, 64 * sizeof(TacsScalar));

  // Copy the A/B matrices to C
  C[0] = A[0];
  C[8] = A[1];
  C[16] = A[2];
  C[1] = A[1];
  C[9] = A[3];
  C[17] = A[4];
  C[2] = A[2];
  C[10] = A[4];
  C[18] = A[5];

  C[3] = B[0];
  C[11] = B[1];
  C[19] = B[2];
  C[4] = B[1];
  C[12] = B[3];
  C[20] = B[4];
  C[5] = B[2];
  C[13] = B[4];
  C[21] = B[5];

  // Copy the B/D matrices to C
  C[24] = B[0];
  C[32] = B[1];
  C[40] = B[2];
  C[25] = B[1];
  C[33] = B[3];
  C[41] = B[4];
  C[26] = B[2];
  C[34] = B[4];
  C[42] = B[5];

  C[27] = D[0];
  C[35] = D[1];
  C[43] = D[2];
  C[28] = D[1];
  C[36] = D[3];
  C[44] = D[4];
  C[29] = D[2];
  C[37] = D[4];
  C[45] = D[5];

  C[54] = As[0];
  C[62] = As[1];
  C[55] = As[1];
  C[63] = As[2];

  int nstress = 8, one = 1, info = 0;
  LAPACKgesv(&nstress, &one, C, &nstress, ipiv, strain, &nstress, &info);

  TacsScalar fail = 0.0;
  failure(gpt, strain, &fail);
  return fail;
}

TacsScalar bladeFSDTStiffness::getBucklingCriterion(
    const TacsScalar stress[8]) {
  TacsScalar strain[8];
  memcpy(strain, stress, 8 * sizeof(TacsScalar));

  double gpt[2] = {0.0, 0.0};
  TacsScalar A[6], B[6], D[6], As[3];
  getStiffness(gpt, A, B, D, As);

  int ipiv[8];
  TacsScalar C[8 * 8];
  memset(C, 0, 64 * sizeof(TacsScalar));

  // Copy the A/B matrices to C
  C[0] = A[0];
  C[8] = A[1];
  C[16] = A[2];
  C[1] = A[1];
  C[9] = A[3];
  C[17] = A[4];
  C[2] = A[2];
  C[10] = A[4];
  C[18] = A[5];

  C[3] = B[0];
  C[11] = B[1];
  C[19] = B[2];
  C[4] = B[1];
  C[12] = B[3];
  C[20] = B[4];
  C[5] = B[2];
  C[13] = B[4];
  C[21] = B[5];

  // Copy the B/D matrices to C
  C[24] = B[0];
  C[32] = B[1];
  C[40] = B[2];
  C[25] = B[1];
  C[33] = B[3];
  C[41] = B[4];
  C[26] = B[2];
  C[34] = B[4];
  C[42] = B[5];

  C[27] = D[0];
  C[35] = D[1];
  C[43] = D[2];
  C[28] = D[1];
  C[36] = D[3];
  C[44] = D[4];
  C[29] = D[2];
  C[37] = D[4];
  C[45] = D[5];

  C[54] = As[0];
  C[62] = As[1];
  C[55] = As[1];
  C[63] = As[2];

  int nstress = 8, one = 1, info = 0;
  LAPACKgesv(&nstress, &one, C, &nstress, ipiv, strain, &nstress, &info);

  TacsScalar fail = 0.0;
  buckling(&fail, strain);
  return fail;
}

void bladeFSDTStiffness::getPanelLoads(const TacsScalar stress[],
                                       TacsScalar* _Nx, TacsScalar* _Nxy) {
  TacsScalar strain[8];
  memcpy(strain, stress, 8 * sizeof(TacsScalar));

  double gpt[2] = {0.0, 0.0};
  TacsScalar A[6], B[6], D[6], As[3];
  getStiffness(gpt, A, B, D, As);

  int ipiv[8];
  TacsScalar C[8 * 8];
  memset(C, 0, 64 * sizeof(TacsScalar));

  // Copy the A/B matrices to C
  C[0] = A[0];
  C[8] = A[1];
  C[16] = A[2];
  C[1] = A[1];
  C[9] = A[3];
  C[17] = A[4];
  C[2] = A[2];
  C[10] = A[4];
  C[18] = A[5];

  C[3] = B[0];
  C[11] = B[1];
  C[19] = B[2];
  C[4] = B[1];
  C[12] = B[3];
  C[20] = B[4];
  C[5] = B[2];
  C[13] = B[4];
  C[21] = B[5];

  // Copy the B/D matrices to C
  C[24] = B[0];
  C[32] = B[1];
  C[40] = B[2];
  C[25] = B[1];
  C[33] = B[3];
  C[41] = B[4];
  C[26] = B[2];
  C[34] = B[4];
  C[42] = B[5];

  C[27] = D[0];
  C[35] = D[1];
  C[43] = D[2];
  C[28] = D[1];
  C[36] = D[3];
  C[44] = D[4];
  C[29] = D[2];
  C[37] = D[4];
  C[45] = D[5];

  C[54] = As[0];
  C[62] = As[1];
  C[55] = As[1];
  C[63] = As[2];

  int nstress = 8, one = 1, info = 0;
  LAPACKgesv(&nstress, &one, C, &nstress, ipiv, strain, &nstress, &info);

  TacsScalar Nx = (Qbar_panel[0] * strain[0] + Qbar_panel[1] * strain[1]);
  *_Nx = (TacsRealPart(Nx) < 0.0 ? -Nx : TacsScalar(0.0));

  *_Nxy = Qbar_panel[5] * strain[2];
}

void bladeFSDTStiffness::getStiffenerLoads(const TacsScalar stress[],
                                           TacsScalar* _Nx, TacsScalar* _Nxy) {
  TacsScalar strain[8];
  memcpy(strain, stress, 8 * sizeof(TacsScalar));

  double gpt[2] = {0.0, 0.0};
  TacsScalar A[6], B[6], D[6], As[3];
  getStiffness(gpt, A, B, D, As);

  int ipiv[8];
  TacsScalar C[8 * 8];
  memset(C, 0, 64 * sizeof(TacsScalar));

  // Copy the A/B matrices to C
  C[0] = A[0];
  C[8] = A[1];
  C[16] = A[2];
  C[1] = A[1];
  C[9] = A[3];
  C[17] = A[4];
  C[2] = A[2];
  C[10] = A[4];
  C[18] = A[5];

  C[3] = B[0];
  C[11] = B[1];
  C[19] = B[2];
  C[4] = B[1];
  C[12] = B[3];
  C[20] = B[4];
  C[5] = B[2];
  C[13] = B[4];
  C[21] = B[5];

  // Copy the B/D matrices to C
  C[24] = B[0];
  C[32] = B[1];
  C[40] = B[2];
  C[25] = B[1];
  C[33] = B[3];
  C[41] = B[4];
  C[26] = B[2];
  C[34] = B[4];
  C[42] = B[5];

  C[27] = D[0];
  C[35] = D[1];
  C[43] = D[2];
  C[28] = D[1];
  C[36] = D[3];
  C[44] = D[4];
  C[29] = D[2];
  C[37] = D[4];
  C[45] = D[5];

  C[54] = As[0];
  C[62] = As[1];
  C[55] = As[1];
  C[63] = As[2];

  int nstress = 8, one = 1, info = 0;
  LAPACKgesv(&nstress, &one, C, &nstress, ipiv, strain, &nstress, &info);

  // Calculate the stiffener buckling criteria
  TacsScalar Nx = (Qbar_stiff[0] * (strain[0] - 0.5 * sh * strain[3]) +
                   Qbar_stiff[1] * (strain[1] - 0.5 * sh * strain[4]));
  *_Nx = (TacsRealPart(Nx) < 0.0 ? -Nx : TacsScalar(0.0));

  // Qxz: strain[6] = gamma_{xz}
  *_Nxy = Qbar_stiff[5] * strain[6];
}

/*
  Get information about the constraints - here just the ply fraction
  constraints.
*/
int bladeFSDTStiffness::isLinear() { return 1; }

// The ply fractions must add up to one
int bladeFSDTStiffness::getNumCon() { return num_pf_con; }

int bladeFSDTStiffness::getConCSRSize() {
  return num_pf_unique + num_stiff_pf_unique;
}

int bladeFSDTStiffness::getConRange(int offset, TacsScalar lb[],
                                    TacsScalar ub[]) {
  if (num_pf_unique > 0) {
    lb[offset] = 1.0;
    ub[offset] = 1.0;
    offset++;
  }

  if (num_stiff_pf_unique > 0) {
    lb[offset] = 1.0;
    ub[offset] = 1.0;
  }

  return num_pf_con;
}

int bladeFSDTStiffness::addConCSR(int offset, int rowp[], int cols[]) {
  int start = rowp[offset];
  if (num_pf_unique > 0) {
    for (int k = 0; k < num_pf_unique; k++, start++) {
      cols[start] = pf_unique[k];
    }
    offset++;
    rowp[offset] = start;
  }

  if (num_stiff_pf_unique > 0) {
    for (int k = 0; k < num_stiff_pf_unique; k++, start++) {
      cols[start] = stiff_pf_unique[k];
    }
    offset++;
    rowp[offset] = start;
  }

  return num_pf_con;
}

int bladeFSDTStiffness::evalCon(int offset, TacsScalar con[]) {
  if (num_pf_unique > 0) {
    con[offset] = 0.0;
    for (int k = 0; k < NUM_PLIES; k++) {
      con[offset] += pf[k];
    }
    offset++;
  }

  if (num_stiff_pf_unique > 0) {
    con[offset] = 0.0;
    for (int k = 0; k < NUM_PLIES; k++) {
      con[offset] += stiff_pf[k];
    }
  }

  return num_pf_con;
}

int bladeFSDTStiffness::evalConDVSens(int offset, TacsScalar Acol[],
                                      const int rowp[], const int cols[]) {
  int start = rowp[offset];

  if (num_pf_unique > 0) {
    for (int k = 0; k < num_pf_unique; k++, start++) {
      Acol[start] = 0.0;
      for (int j = 0; j < NUM_PLIES; j++) {
        if (pf_unique[k] == pf_nums[j]) {
          Acol[start] += 1.0;
        }
      }
    }
  }

  if (num_stiff_pf_unique > 0) {
    for (int k = 0; k < num_stiff_pf_unique; k++, start++) {
      Acol[start] = 0.0;
      for (int j = 0; j < NUM_PLIES; j++) {
        if (stiff_pf_unique[k] == stiff_pf_nums[j]) {
          Acol[start] += 1.0;
        }
      }
    }
  }

  return num_pf_con;
}

// Design variable information set up
// ----------------------------------
int bladeFSDTStiffness::ownsDesignVar(const int dvNum) const {
  if (Lx_num == dvNum) {
    return 1;
  } else if (sp_num == dvNum) {
    return 1;
  } else if (sh_num == dvNum) {
    return 1;
  } else if (st_num == dvNum) {
    return 1;
  } else if (t_num == dvNum) {
    return 1;
  }

  for (int k = 0; k < NUM_PLIES; k++) {
    if (pf_nums[k] == dvNum || stiff_pf_nums[k] == dvNum) {
      return 1;
    }
  }
  return 0;
}

void bladeFSDTStiffness::setDesignVars(const TacsScalar dvs[], int numDVs) {
  if (Lx_num >= 0 && Lx_num < numDVs) {
    Lx = dvs[Lx_num];
  }
  if (sp_num < numDVs) {
    sp = dvs[sp_num];
  }
  if (sh_num < numDVs) {
    sh = dvs[sh_num];
  }
  if (st_num < numDVs) {
    st = dvs[st_num];
  }
  if (t_num < numDVs) {
    t = dvs[t_num];
  }

  for (int k = 0; k < NUM_PLIES; k++) {
    if (pf_nums[k] >= 0 && pf_nums[k] < numDVs) {
      pf[k] = dvs[pf_nums[k]];
    }
    if (stiff_pf_nums[k] >= 0 && stiff_pf_nums[k] < numDVs) {
      stiff_pf[k] = dvs[stiff_pf_nums[k]];
    }
  }

  updateStiffness();
}

/*
  Set the stiffness properties
*/
void bladeFSDTStiffness::updateStiffness() {
  // Zero the stiffness properties
  for (int j = 0; j < 6; j++) {
    Qbar_panel[j] = Qbar_stiff[j] = 0.0;
  }

  for (int j = 0; j < 3; j++) {
    Abar_panel[j] = Abar_stiff[j] = 0.0;
  }

  // Calculate the constitutive properties of the panel
  for (int k = 0; k < NUM_PLIES; k++) {
    TacsScalar Qp[6], Ap[3];

    // Calculate the stiffness contribution from the panel
    ortho_ply->calculateQbar(Qp, ply_angles[k]);
    ortho_ply->calculateAbar(Ap, ply_angles[k]);

    for (int j = 0; j < 6; j++) {
      Qbar_panel[j] += pf[k] * Qp[j];
    }

    for (int j = 0; j < 3; j++) {
      Abar_panel[j] += pf[k] * Ap[j];
    }

    // Calculate the stiffness contribution from the stiffener
    ortho_ply->calculateQbar(Qp, stiff_ply_angles[k]);
    ortho_ply->calculateAbar(Ap, stiff_ply_angles[k]);

    for (int j = 0; j < 6; j++) {
      Qbar_stiff[j] += stiff_pf[k] * Qp[j];
    }

    for (int j = 0; j < 3; j++) {
      Abar_stiff[j] += stiff_pf[k] * Ap[j];
    }
  }

  A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = 0.0;
  B[0] = B[1] = B[2] = B[3] = B[4] = B[5] = 0.0;
  D[0] = D[1] = D[2] = D[3] = D[4] = D[5] = 0.0;
  As[0] = As[1] = As[2] = 0.0;

  // Add the contributions from the panel
  TacsScalar It = (t * t * t) / 12.0;

  for (int k = 0; k < 6; k++) {
    A[k] = t * Qbar_panel[k];
    B[k] = 0.0;
    D[k] = It * Qbar_panel[k];
  }

  As[0] = t * kcorr * Abar_panel[0];
  As[1] = t * kcorr * Abar_panel[1];
  As[2] = t * kcorr * Abar_panel[2];

  k_penalty = DRILLING_REGULARIZATION * Qbar_panel[5] * t;

  // Add the contributions from the stiffener
  // st == stiffener thickness
  // sh == stiffener height
  // sf_fraction == stiffener flange fraction
  TacsScalar Area = st * sh * (1.0 + sf_fraction);
  TacsScalar es = -0.5 * sh;

  // The moment of inertia of the stiffener about the neutral axis of
  // the plate
  TacsScalar Is = st * (sh * sh * sh / 12.0);
  TacsScalar Js = sh * (st * st * st / 3.0);

  // The effective Young's modulus of the stiffener
  TacsScalar Es =
      (Qbar_stiff[0] - (Qbar_stiff[1] * Qbar_stiff[1]) / Qbar_stiff[3]);
  TacsScalar Gs = Qbar_stiff[5];

  // Add the contributions to the stiffness
  TacsScalar sfact = 1.0 / sp;
  A[0] += (Es * Area) * sfact;
  B[0] += (es * Es * Area) * sfact;
  D[0] += Es * (Area * es * es + Is) * sfact;
  D[5] += Js * Gs * sfact;
  As[0] += kcorr * Gs * sh * st * sfact;

  if (TacsRealPart(B[0] * B[0]) > TacsRealPart(A[0] * D[0])) {
    fprintf(stderr, "bladeFSDTStiffness: A11*D11 - B11^2 < 0\n");
    fprintf(stderr, "t = %8.3f sh = %8.3f st = %8.3f sp = %8.3f\n",
            TacsRealPart(t), TacsRealPart(sh), TacsRealPart(st),
            TacsRealPart(sp));
    fprintf(stderr, "A11 = %10.3e B11 = %10.3e D11 = %10.3e\n",
            TacsRealPart(A[0]), TacsRealPart(B[0]), TacsRealPart(D[0]));
  }

  calcBucklingProps();
}

int bladeFSDTStiffness::getLxNum() { return Lx_num; }

void bladeFSDTStiffness::getDesignVars(TacsScalar dvs[], int numDVs) {
  if (Lx_num >= 0 && Lx_num < numDVs) {
    dvs[Lx_num] = Lx;
  }
  if (sp_num >= 0 && sp_num < numDVs) {
    dvs[sp_num] = sp;
  }
  if (sh_num >= 0 && sh_num < numDVs) {
    dvs[sh_num] = sh;
  }
  if (st_num >= 0 && st_num < numDVs) {
    dvs[st_num] = st;
  }
  if (t_num >= 0 && t_num < numDVs) {
    dvs[t_num] = t;
  }

  for (int k = 0; k < NUM_PLIES; k++) {
    if (pf_nums[k] >= 0 && pf_nums[k] < numDVs) {
      dvs[pf_nums[k]] = pf[k];
    }
    if (stiff_pf_nums[k] >= 0 && stiff_pf_nums[k] < numDVs) {
      dvs[stiff_pf_nums[k]] = stiff_pf[k];
    }
  }
}

void bladeFSDTStiffness::getDesignVarRange(TacsScalar lb[], TacsScalar ub[],
                                           int numDVs) {
  if (Lx_num >= 0 && Lx_num < numDVs) {
    lb[Lx_num] = 0.0;
    ub[Lx_num] = 1e20;
  }
  if (sp_num < numDVs) {
    lb[sp_num] = sp_low;
    ub[sp_num] = sp_high;
  }
  if (sh_num < numDVs) {
    lb[sh_num] = sh_low;
    ub[sh_num] = sh_high;
  }
  if (st_num < numDVs) {
    lb[st_num] = st_low;
    ub[st_num] = st_high;
  }
  if (t_num < numDVs) {
    lb[t_num] = t_low;
    ub[t_num] = t_high;
  }

  for (int k = 0; k < NUM_PLIES; k++) {
    if (pf_nums[k] >= 0 && pf_nums[k] < numDVs) {
      lb[pf_nums[k]] = pf_low[k];
      ub[pf_nums[k]] = pf_high[k];
    }
  }

  for (int k = 0; k < NUM_PLIES; k++) {
    if (stiff_pf_nums[k] >= 0 && stiff_pf_nums[k] < numDVs) {
      lb[stiff_pf_nums[k]] = stiff_pf_low[k];
      ub[stiff_pf_nums[k]] = stiff_pf_high[k];
    }
  }
}

/*
  Calculate the mass per unit area and the moment of inertia per unit
  area of the plate.
*/
void bladeFSDTStiffness::getPointwiseMass(const double gpt[],
                                          TacsScalar mass[]) {
  mass[0] = ortho_ply->getRho() * t;
  mass[1] = 0.0;
  mass[2] = ortho_ply->getRho() * t * t * t / 12.0;

  TacsScalar Area = st * sh * (1.0 + sf_fraction);
  TacsScalar sfact = 1.0 / sp;

  // Add the contribution to the mass
  mass[0] += ortho_ply->getRho() * Area * sfact;
  mass[2] += 0.25 * ortho_ply->getRho() * st * st * Area * sfact;
}

/*
  Calculate the sensitivity of the pointwise mass
*/
void bladeFSDTStiffness::pointwiseMassDVSens(int dvNum, const double gpt[],
                                             TacsScalar massDVSens[]) {
  massDVSens[0] = massDVSens[1] = massDVSens[2] = 0.0;
  if (dvNum >= 0) {
    if (t_num == dvNum) {
      massDVSens[0] = ortho_ply->getRho();
      massDVSens[2] = ortho_ply->getRho() * t * t / 4.0;
    } else if (sp_num == dvNum) {
      TacsScalar Area = st * sh * (1.0 + sf_fraction);
      TacsScalar sfact_sens = -1.0 / (sp * sp);

      massDVSens[0] = ortho_ply->getRho() * Area * sfact_sens;
      massDVSens[2] = 0.25 * ortho_ply->getRho() * st * st * Area * sfact_sens;
    } else if (sh_num == dvNum) {
      TacsScalar Area_sens = st * (1.0 + sf_fraction);
      TacsScalar sfact = 1.0 / sp;

      massDVSens[0] = ortho_ply->getRho() * Area_sens * sfact;
      massDVSens[2] = 0.25 * ortho_ply->getRho() * st * st * Area_sens * sfact;
    } else if (st_num == dvNum) {
      TacsScalar Area_sens = sh * (1.0 + sf_fraction);
      TacsScalar sfact = 1.0 / sp;

      massDVSens[0] += ortho_ply->getRho() * Area_sens * sfact;
      massDVSens[2] =
          (3.0 / 4.0) * ortho_ply->getRho() * st * st * Area_sens * sfact;
    }
  }
}

TacsScalar bladeFSDTStiffness::getDensity() { return ortho_ply->getRho(); }

/*
  Add the derivative of the pointwise mass times alpha to each component
*/
void bladeFSDTStiffness::addPointwiseMassDVSens(const double pt[],
                                                const TacsScalar alpha[],
                                                TacsScalar dvSens[],
                                                int dvLen) {
  TacsScalar rho = ortho_ply->getRho();

  if (t_num >= 0 && t_num < dvLen) {
    dvSens[t_num] += alpha[0] * rho;
  }
  if (sp_num >= 0 && sp_num < dvLen) {
    TacsScalar Area = st * sh * (1.0 + sf_fraction);
    TacsScalar sfact_sens = -1.0 / (sp * sp);

    dvSens[sp_num] += alpha[0] * rho * Area * sfact_sens;
  }
  if (sh_num >= 0 && sh_num < dvLen) {
    TacsScalar Area_sens = st * (1.0 + sf_fraction);
    TacsScalar sfact = 1.0 / sp;

    dvSens[sh_num] += alpha[0] * rho * Area_sens * sfact;
  }
  if (st_num >= 0 && st_num < dvLen) {
    TacsScalar Area_sens = sh * (1.0 + sf_fraction);
    TacsScalar sfact = 1.0 / sp;

    dvSens[st_num] += alpha[0] * rho * Area_sens * sfact;
  }
}

/*
  Calculate the stiffness properties of the panel
*/
TacsScalar bladeFSDTStiffness::getStiffness(const double gpt[], TacsScalar _A[],
                                            TacsScalar _B[], TacsScalar _D[],
                                            TacsScalar _As[]) {
  for (int k = 0; k < 6; k++) {
    _A[k] = A[k];
    _B[k] = B[k];
    _D[k] = D[k];
  }

  for (int k = 0; k < 3; k++) {
    _As[k] = As[k];
  }

  return k_penalty;
}

/*
  Calculate the derivative w.r.t. one of the design variables:

  sp: the stiffener pitch
  sh: the stiffener height
  st: the stiffener thickness
  t: the panel thickness

  Note that sensitivities w.r.t. the ply fractions are handled
  elsewhere.
*/
TacsScalar bladeFSDTStiffness::getStiffnessDVSens(int dvNum, const double gpt[],
                                                  TacsScalar sA[],
                                                  TacsScalar sB[],
                                                  TacsScalar sD[],
                                                  TacsScalar sAs[]) {
  sA[0] = sA[1] = sA[2] = sA[3] = sA[4] = sA[5] = 0.0;
  sB[0] = sB[1] = sB[2] = sB[3] = sB[4] = sB[5] = 0.0;
  sD[0] = sD[1] = sD[2] = sD[3] = sD[4] = sD[5] = 0.0;
  sAs[0] = sAs[1] = sAs[2] = 0.0;

  TacsScalar k_penalty_sens = 0.0;

  if (dvNum >= 0 && ownsDesignVar(dvNum)) {
    if (t_num == dvNum) {
      // Calcualte the senstivity w.r.t. the panel thickness
      TacsScalar It = (t * t) / 4.0;

      for (int k = 0; k < 6; k++) {
        sA[k] = Qbar_panel[k];
        sB[k] = 0.0;
        sD[k] = It * Qbar_panel[k];
      }

      sAs[0] = kcorr * Abar_panel[0];
      sAs[1] = kcorr * Abar_panel[1];
      sAs[2] = kcorr * Abar_panel[2];

      k_penalty_sens = DRILLING_REGULARIZATION * Qbar_panel[5];
    } else if (sp_num == dvNum) {
      // Calculate the derivative of the stiffness w.r.t. stiffener pitch
      TacsScalar Area = st * sh * (1.0 + sf_fraction);
      TacsScalar es = -0.5 * sh;
      TacsScalar Is = st * (sh * sh * sh / 12.0);
      TacsScalar Js = sh * (st * st * st / 3.0);

      // The effective Young's modulus of the stiffener
      TacsScalar Es =
          (Qbar_stiff[0] - (Qbar_stiff[1] * Qbar_stiff[1]) / Qbar_stiff[3]);
      TacsScalar Gs = Qbar_stiff[5];

      // Add the contributions to the stiffness
      TacsScalar sfact_sens = -1.0 / (sp * sp);
      sA[0] += (Es * Area) * sfact_sens;
      sB[0] += (es * Es * Area) * sfact_sens;
      sD[0] += Es * (Area * es * es + Is) * sfact_sens;
      sD[5] += Js * Gs * sfact_sens;
      sAs[0] += kcorr * Gs * sh * st * sfact_sens;
    } else if (sh_num == dvNum) {
      // Calculate the derivative of the stiffness w.r.t. stiffener height
      TacsScalar Area = st * sh * (1.0 + sf_fraction);
      TacsScalar es = -0.5 * sh;

      TacsScalar Area_sens = st * (1.0 + sf_fraction);
      TacsScalar es_sens = -0.5;
      TacsScalar Is_sens = st * sh * sh / 4.0;
      TacsScalar Js_sens = (st * st * st / 3.0);

      // The effective Young's modulus of the stiffener
      TacsScalar Es =
          (Qbar_stiff[0] - (Qbar_stiff[1] * Qbar_stiff[1]) / Qbar_stiff[3]);
      TacsScalar Gs = Qbar_stiff[5];

      // Add the contributions to the stiffness
      TacsScalar sfact = 1.0 / sp;
      sA[0] += (Es * Area_sens) * sfact;
      sB[0] += Es * (es_sens * Area + es * Area_sens) * sfact;
      sD[0] += Es *
               (Area_sens * es * es + 2.0 * Area * es * es_sens + Is_sens) *
               sfact;
      sD[5] += Js_sens * Gs * sfact;
      sAs[0] += kcorr * Gs * st * sfact;
    } else if (st_num == dvNum) {
      // calculate the derivative of the stiffness w.r.t.
      // the stiffener thickness
      TacsScalar es = -0.5 * sh;

      TacsScalar Area_sens = sh * (1.0 + sf_fraction);
      TacsScalar Is_sens = (sh * sh * sh / 12.0);
      TacsScalar Js_sens = sh * st * st;

      // The effective Young's modulus of the stiffener
      TacsScalar Es =
          (Qbar_stiff[0] - (Qbar_stiff[1] * Qbar_stiff[1]) / Qbar_stiff[3]);
      TacsScalar Gs = Qbar_stiff[5];

      // Add the contributions to the stiffness
      TacsScalar sfact = 1.0 / sp;
      sA[0] += (Es * Area_sens) * sfact;
      sB[0] += (es * Es * Area_sens) * sfact;
      sD[0] += Es * (Area_sens * es * es + Is_sens) * sfact;
      sD[5] += Js_sens * Gs * sfact;
      sAs[0] += kcorr * Gs * sh * sfact;
    } else {
      TacsScalar sAbar_panel[3], sQbar_panel[6];
      TacsScalar sAbar_stiff[3], sQbar_stiff[6];
      sAbar_panel[0] = sAbar_panel[1] = sAbar_panel[2] = 0.0;
      sAbar_stiff[0] = sAbar_stiff[1] = sAbar_stiff[2] = 0.0;

      sQbar_panel[0] = sQbar_panel[1] = sQbar_panel[2] = 0.0;
      sQbar_panel[3] = sQbar_panel[4] = sQbar_panel[5] = 0.0;

      sQbar_stiff[0] = sQbar_stiff[1] = sQbar_stiff[2] = 0.0;
      sQbar_stiff[3] = sQbar_stiff[4] = sQbar_stiff[5] = 0.0;

      for (int k = 0; k < NUM_PLIES; k++) {
        if (pf_nums[k] == dvNum) {
          TacsScalar Qp[6], Ap[3];
          ortho_ply->calculateQbar(Qp, ply_angles[k]);
          ortho_ply->calculateAbar(Ap, ply_angles[k]);

          for (int j = 0; j < 6; j++) {
            sQbar_panel[j] += Qp[j];
          }

          for (int j = 0; j < 3; j++) {
            sAbar_panel[j] += Ap[j];
          }
        }
      }

      // Add the contributions from the panel
      TacsScalar It = (t * t * t) / 12.0;

      for (int k = 0; k < 6; k++) {
        sA[k] = t * sQbar_panel[k];
        sB[k] = 0.0;
        sD[k] = It * sQbar_panel[k];
      }

      sAs[0] = t * kcorr * sAbar_panel[0];
      sAs[1] = t * kcorr * sAbar_panel[1];
      sAs[2] = t * kcorr * sAbar_panel[2];

      for (int k = 0; k < NUM_PLIES; k++) {
        if (stiff_pf_nums[k] == dvNum) {
          TacsScalar Qp[6], Ap[3];
          ortho_ply->calculateQbar(Qp, stiff_ply_angles[k]);
          ortho_ply->calculateAbar(Ap, stiff_ply_angles[k]);

          for (int j = 0; j < 6; j++) {
            sQbar_stiff[j] += Qp[j];
          }

          for (int j = 0; j < 3; j++) {
            sAbar_stiff[j] += Ap[j];
          }
        }
      }

      TacsScalar Area = st * sh * (1.0 + sf_fraction);
      TacsScalar es = -0.5 * sh;
      TacsScalar Is = st * (sh * sh * sh / 12.0);
      TacsScalar Js = sh * (st * st * st / 3.0);

      TacsScalar Es_sens =
          (sQbar_stiff[0] -
           2.0 * (Qbar_stiff[1] * sQbar_stiff[1]) / Qbar_stiff[3] +
           sQbar_stiff[3] * (Qbar_stiff[1] * Qbar_stiff[1]) /
               (Qbar_stiff[3] * Qbar_stiff[3]));
      TacsScalar Gs_sens = sQbar_stiff[5];

      TacsScalar sfact = 1.0 / sp;
      sA[0] += (Es_sens * Area) * sfact;
      sB[0] += (es * Es_sens * Area) * sfact;
      sD[0] += Es_sens * (Area * es * es + Is) * sfact;
      sD[5] += Js * Gs_sens * sfact;
      sAs[0] += kcorr * Gs_sens * sh * st * sfact;

      k_penalty_sens = DRILLING_REGULARIZATION * sQbar_panel[5] * t;
    }
  }

  return k_penalty_sens;
}

/*
  Add the derivative of the product of the stress to the design
  derivative array
*/
void bladeFSDTStiffness::addStiffnessDVSens(const double pt[],
                                            const TacsScalar e[],
                                            const TacsScalar psi[],
                                            TacsScalar rotPsi,
                                            TacsScalar fdvSens[], int dvLen) {
  int panel_nums[5] = {t_num, sp_num, sh_num, st_num, Lx_num};
  int dvNum;
  // Store the derivative of the stress values
  TacsScalar s[8], ksens;

  for (int k = 0; k < 5 + 2 * NUM_PLIES; k++) {
    if (k < 5) {
      dvNum = panel_nums[k];
    } else if (5 <= k && k < NUM_PLIES + 5) {
      dvNum = pf_nums[k - 5];
    } else {
      dvNum = stiff_pf_nums[k - (NUM_PLIES + 5)];
    }

    if (dvNum >= 0 && dvNum < dvLen) {
      ksens = calculateStressDVSens(pt, e, s, dvNum);

      // Add the result to the design variable vector
      fdvSens[dvNum] += (s[0] * psi[0] + s[1] * psi[1] + s[2] * psi[2] +
                         s[3] * psi[3] + s[4] * psi[4] + s[5] * psi[5] +
                         s[6] * psi[6] + s[7] * psi[7] + rotPsi * ksens);
    }
  }
}

/*
  Add the derivative of the product of the stress to the design
  derivative array
*/
TacsScalar bladeFSDTStiffness::calculateStressDVSens(const double pt[],
                                                     const TacsScalar strain[],
                                                     TacsScalar s_sens[],
                                                     int dvNum) {
  // Store the derivative of the stress values
  TacsScalar A_sens[6], B_sens[6], D_sens[6], As_sens[3];
  TacsScalar ksens;

  ksens = getStiffnessDVSens(dvNum, pt, A_sens, B_sens, D_sens, As_sens);

  // Compute the in-plane resultants
  s_sens[0] = A_sens[0] * strain[0] + A_sens[1] * strain[1] +
              A_sens[2] * strain[2] + B_sens[0] * strain[3] +
              B_sens[1] * strain[4] + B_sens[2] * strain[5];
  s_sens[1] = A_sens[1] * strain[0] + A_sens[3] * strain[1] +
              A_sens[4] * strain[2] + B_sens[1] * strain[3] +
              B_sens[3] * strain[4] + B_sens[4] * strain[5];
  s_sens[2] = A_sens[2] * strain[0] + A_sens[4] * strain[1] +
              A_sens[5] * strain[2] + B_sens[2] * strain[3] +
              B_sens[4] * strain[4] + B_sens[5] * strain[5];

  // Compute the bending moments
  s_sens[3] = B_sens[0] * strain[0] + B_sens[1] * strain[1] +
              B_sens[2] * strain[2] + D_sens[0] * strain[3] +
              D_sens[1] * strain[4] + D_sens[2] * strain[5];
  s_sens[4] = B_sens[1] * strain[0] + B_sens[3] * strain[1] +
              B_sens[4] * strain[2] + D_sens[1] * strain[3] +
              D_sens[3] * strain[4] + D_sens[4] * strain[5];
  s_sens[5] = B_sens[2] * strain[0] + B_sens[4] * strain[1] +
              B_sens[5] * strain[2] + D_sens[2] * strain[3] +
              D_sens[4] * strain[4] + D_sens[5] * strain[5];

  // Compute the shear resultants
  s_sens[6] = As_sens[0] * strain[6] + As_sens[1] * strain[7];
  s_sens[7] = As_sens[1] * strain[6] + As_sens[2] * strain[7];

  return ksens;
}

/*
  The shear buckling loads are calculated based on an infinite-plate
  solution simply supported along the panel sides. For the local skin
  buckling calculations we use the infinite plate solution along edges
  supported by the blade stiffeners, therefore the panel width is
  equal to the stiffener pitch (sp). For the global-level
  calculations, we use the panel length equal to the rib pitch, thus
  the panel width is equal to the local panel length (Lx).

  In Stroud and Arganoff, the following formula are suggested for the
  calculation of the buckling loads:

  xi = sqrt(D1*D2)/D3

  if xi > 1.0:
    Nxy,crit = (4.0/Ly^2)*(D2*D1^3)^(1/4)*(8.125 + 5.05/xi)
  else:
    Nxy,crit = (4.0/Ly^2)*sqrt(D1*D3)*(11.7 + 0.532*xi + 0.938*xi^2)

  Note that if xi = 1, then D3 = sqrt(D1*D2) and so
  (4.0/Ly^2)*sqrt(D1*D3) = (4.0/Ly^2)*(D1*D1^3)^(1/4)

  However, when xi = 1, the two terms inside the brackets are not
  equal. As a result, there is a discontinuity between the two
  expressions. To avoid this, we adjust the first formula in a
  conservative manner. As shown in Lekhnitskii, the limit for xi ->
  infty is 8.125, while for xi = 1 is 13.17 and for xi = 0, 11.71.

  The formula in Stroud and Arganoff do not meet these end conditions.
  We adjust the formula as follows:

  if xi > 1.0:
    Nxy,crit = (4.0/Ly^2)*(D2*D1^3)^(1/4)*(8.125 + 5.045/xi)
  else:
    Nxy,crit = (4.0/Ly^2)*sqrt(D1*D3)*(11.7 + 0.532*xi + 0.938*xi^2)

  input:
  D1:        the longitudinal bending stiffness
  D2:        the transverse bending stiffness
  D3:        the shear bending stiffness
  L:         the side-length of the panel

  returns:
  Nxy_crit:  the approximate critical buckling load
*/
TacsScalar bladeFSDTStiffness::calcCriticalShearLoad(TacsScalar D1,
                                                     TacsScalar D2,
                                                     TacsScalar D3,
                                                     TacsScalar L) {
  TacsScalar xi = sqrt(D1 * D2) / D3;

  double ks = 50.0;
  TacsScalar Nxy_crit_1 =
      (4.0 / (L * L)) * sqrt(D3 * D1 * xi) * (8.125 + 5.045 / xi);
  TacsScalar Nxy_crit_2 =
      (4.0 / (L * L)) * sqrt(D1 * D3) * (11.7 + 0.532 * xi + 0.938 * xi * xi);

  TacsScalar Nxy_min = 0.0;
  if (TacsRealPart(Nxy_crit_1) < TacsRealPart(Nxy_crit_2)) {
    Nxy_min = Nxy_crit_1;
  } else {
    Nxy_min = Nxy_crit_2;
  }

  TacsScalar Nxy_diff = fabs(Nxy_crit_1 - Nxy_crit_2);
  TacsScalar Nxy_crit = Nxy_min - log(1.0 + exp(-ks * Nxy_diff)) / ks;

  return Nxy_crit;
}

/*
  Calculate the sensitivity of the shear buckling load w.r.t. the
  inputs. This is the differentiated version of the above code.
*/
TacsScalar bladeFSDTStiffness::calcCriticalShearLoadSens(
    TacsScalar D1, TacsScalar D2, TacsScalar D3, TacsScalar L, TacsScalar sD1,
    TacsScalar sD2, TacsScalar sD3, TacsScalar sL) {
  TacsScalar xi = sqrt(D1 * D2) / D3;
  TacsScalar sxi =
      (0.5 * (D1 * sD2 + D2 * sD1) - xi * xi * D3 * sD3) / (xi * D3 * D3);

  double ks = 50.0;
  TacsScalar Nxy_crit_1 =
      (4.0 / (L * L)) * sqrt(D3 * D1 * xi) * (8.125 + 5.045 / xi);
  TacsScalar Nxy_crit_2 =
      (4.0 / (L * L)) * sqrt(D1 * D3) * (11.7 + 0.532 * xi + 0.938 * xi * xi);

  TacsScalar Nxy_min = 0.0;
  if (TacsRealPart(Nxy_crit_1) < TacsRealPart(Nxy_crit_2)) {
    Nxy_min = Nxy_crit_1;
  } else {
    Nxy_min = Nxy_crit_2;
  }

  TacsScalar Nxy_diff = fabs(Nxy_crit_1 - Nxy_crit_2);
  TacsScalar sum = 1.0 + exp(-ks * Nxy_diff);

  // Compute the derivative w.r.t. the first critical load
  TacsScalar y = sqrt(D3 * D1 * xi);
  TacsScalar sy = 0.5 * (sD3 * D1 * xi + sD1 * D3 * xi + sxi * D1 * D3) / y;
  TacsScalar cb = (8.125 + 5.045 / xi);

  TacsScalar Nxy_crit_1_sens =
      (4.0 / (L * L)) * (sy * cb - 5.045 * sxi * y / (xi * xi));
  if (sL != 0.0) {
    Nxy_crit_1_sens -= 8.0 * sL / (L * L * L) * y * cb;
  }

  // Compute the derivative w.r.t. the second critical load
  y = sqrt(D1 * D3);
  sy = 0.5 * (sD1 * D3 + sD3 * D1) / y;
  cb = (11.7 + 0.532 * xi + 0.938 * xi * xi);

  TacsScalar Nxy_crit_2_sens =
      (4.0 / (L * L)) * (sy * cb + y * (0.532 * sxi + 2.0 * 0.938 * xi * sxi));
  if (sL != 0.0) {
    Nxy_crit_2_sens -= 8.0 * sL / (L * L * L) * y * cb;
  }

  // Compute the derivative due to the KS function
  TacsScalar Nxy_crit_sens = 0.0;
  if (TacsRealPart(Nxy_crit_1) < TacsRealPart(Nxy_crit_2)) {
    Nxy_crit_sens =
        (Nxy_crit_1_sens + Nxy_crit_2_sens * exp(-ks * Nxy_diff)) / sum;
  } else {
    Nxy_crit_sens =
        (Nxy_crit_1_sens * exp(-ks * Nxy_diff) + Nxy_crit_2_sens) / sum;
  }

  return Nxy_crit_sens;
}

/*
  Calculate the buckling properties of the panel

  The global buckling loads for the smeared panel:
  Nx_global_crit: The critical axial buckling load
  Nxy_global_crit: The critical shear buckling load

  The buckling loads for the skin in between stiffeners:
  Nx_panel_crit: The critical buckling load of the panel
  Nxy_panel_crit: The critical shear buckling load of the panel

  The buckling load of the stiffener:
  Nx_weak_stiff_crit: The critical lower buckling load of the stiffener
  Nx_strong_stiff_crit: The criticial higher buckling load of the stiffener
*/
void bladeFSDTStiffness::calcBucklingProps() {
  // Calculate the global buckling mode
  // ----------------------------------
  // Compute the neutral axis location:
  TacsScalar Ep =
      Qbar_panel[0] - (Qbar_panel[1] * Qbar_panel[1]) / Qbar_panel[3];
  TacsScalar Es =
      Qbar_stiff[0] - (Qbar_stiff[1] * Qbar_stiff[1]) / Qbar_stiff[3];

  // The modulus-weighted area:
  TacsScalar An = (Ep * sp * t + Es * sh * st * (1.0 + sf_fraction));

  // The modulus-weighted moment of area
  TacsScalar Cn = 0.5 * Es * sh * sh * st;

  // The modulus-weighted centroid
  TacsScalar zn = Cn / An;

  // The bending stiffness
  TacsScalar EI = (zn * zn * (t * sp * Ep + st * sh * sf_fraction * Es) +
                   Es * (st * (sh * sh * sh) / 12.0 +
                         st * sh * (zn - 0.5 * sh) * (zn - 0.5 * sh)));

  // Compute the unit EI stiffness for the buckling
  TacsScalar D1 = EI / sp;
  Nx_global_crit = (M_PI * M_PI * D1) / (Lx * Lx);

  // Compute the required data for the shear loading
  // The transverse bending stiffness
  TacsScalar D12_panel = (t * t * t * Qbar_panel[1]) / 12.0;
  TacsScalar D12_stiff = ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
                          (Qbar_panel[1] + Qbar_stiff[1])) /
                         6.0;

  TacsScalar D2 = sp / ((sp - sh * sf_fraction) / D12_panel +
                        (sh * sf_fraction) / D12_stiff);

  // The twisting stiffness
  TacsScalar D66_panel = (t * t * t * Qbar_panel[5]) / 12.0;
  TacsScalar D66_stiff = ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
                          (Qbar_panel[5] + Qbar_stiff[5])) /
                         6.0;

  TacsScalar D3 = sp / ((sp - sh * sf_fraction) / D66_panel +
                        (sh * sf_fraction) / D66_stiff);
  Nxy_global_crit = calcCriticalShearLoad(D1, D2, D3, Lx);

  // Compute the critical buckling loads for the panel
  // -------------------------------------------------
  TacsScalar I = (t * t * t) / 12.0;
  TacsScalar D11 = I * Qbar_panel[0];
  TacsScalar D12 = I * Qbar_panel[1];
  TacsScalar D22 = I * Qbar_panel[3];
  TacsScalar D66 = I * Qbar_panel[5];

  Nx_panel_crit =
      (2.0 * M_PI * M_PI / (sp * sp)) * (sqrt(D11 * D22) + D12 + 2.0 * D66);
  Nxy_panel_crit = calcCriticalShearLoad(D11, D22, D12 + 2.0 * D66, sp);

  // Compute the critical loads for the stiffener
  // --------------------------------------------
  TacsScalar Is_weak = sh * st * st * st / 12.0;
  TacsScalar Is_strong = st * sh * sh * sh / 12.0;

  Nx_weak_stiff_crit = (M_PI * M_PI * Es * Is_weak) / (Lx * Lx);
  Nx_strong_stiff_crit = (M_PI * M_PI * Es * Is_strong) / (Lx * Lx);
}

/*
  Calculate the derivative of the critical buckling loads w.r.t. the
  design variables

  This function computes the derivative w.r.t:

  t_num: the panel thickness
  sh_num: the stiffener height
  st_num: the stiffener thickness
  sp_num: the stiffener pitch
  Lx_num: the length of the panel

  And... any remaining ply fraction variables through sQbar_panel and
  sQbar_stiff
*/
void bladeFSDTStiffness::calcBucklingPropsDVSens(
    int dvNum, TacsScalar* sNx_global_crit, TacsScalar* sNxy_global_crit,
    TacsScalar* sNx_panel_crit, TacsScalar* sNxy_panel_crit,
    TacsScalar* sNx_weak_stiff_crit, TacsScalar* sNx_strong_stiff_crit) {
  if (dvNum == t_num) {
    // Compute the neutral axis location:
    TacsScalar Ep =
        Qbar_panel[0] - (Qbar_panel[1] * Qbar_panel[1]) / Qbar_panel[3];
    TacsScalar Es =
        Qbar_stiff[0] - (Qbar_stiff[1] * Qbar_stiff[1]) / Qbar_stiff[3];

    // The modulus-weighted area:
    TacsScalar An = (Ep * sp * t + Es * sh * st + Es * sh * st * sf_fraction);
    TacsScalar An_sens = Ep * sp;

    // The modulus-weighted centroid
    TacsScalar Cn = 0.5 * Es * sh * sh * st;

    // The modulus-weighted centroid
    TacsScalar zn = Cn / An;
    TacsScalar zn_sens = -Cn / (An * An) * An_sens;

    // The bending stiffness
    TacsScalar EI = (zn * zn * (t * sp * Ep + st * sh * sf_fraction * Es) +
                     Es * (st * (sh * sh * sh) / 12.0 +
                           st * sh * (zn - 0.5 * sh) * (zn - 0.5 * sh)));
    TacsScalar EI_sens =
        (2.0 * zn * zn_sens * (t * sp * Ep + st * sh * sf_fraction * Es) +
         zn * zn * sp * Ep + Es * (2.0 * st * sh * (zn - 0.5 * sh) * zn_sens));

    // Compute the unit EI stiffness for the buckling
    TacsScalar D1 = EI / sp;
    TacsScalar D1_sens = EI_sens / sp;

    // Compute the sensitivity of the global axial load
    *sNx_global_crit = (M_PI * M_PI * D1_sens) / (Lx * Lx);

    // Compute the required data for the shear loading
    // The transverse bending stiffness
    TacsScalar D12_panel = (t * t * t * Qbar_panel[1]) / 12.0;
    TacsScalar D12_stiff = ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
                            (Qbar_panel[1] + Qbar_stiff[1])) /
                           6.0;

    TacsScalar D2 = sp / ((sp - sh * sf_fraction) / D12_panel +
                          (sh * sf_fraction) / D12_stiff);

    TacsScalar D12_panel_sens = (t * t * Qbar_panel[1]) / 4.0;
    TacsScalar D12_stiff_sens =
        ((t + 0.5 * st) * (t + 0.5 * st) * (Qbar_panel[1] + Qbar_stiff[1])) /
        2.0;

    TacsScalar D2_sens =
        (D2 * D2 / sp) *
        ((sp - sh * sf_fraction) / (D12_panel * D12_panel) * D12_panel_sens +
         (sh * sf_fraction) / (D12_stiff * D12_stiff) * D12_stiff_sens);

    // The twisting stiffness
    TacsScalar D66_panel = (t * t * t * Qbar_panel[5]) / 12.0;
    TacsScalar D66_stiff = ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
                            (Qbar_panel[5] + Qbar_stiff[5])) /
                           6.0;

    TacsScalar D3 = sp / ((sp - sh * sf_fraction) / D66_panel +
                          (sh * sf_fraction) / D66_stiff);

    TacsScalar D66_panel_sens = (t * t * Qbar_panel[5]) / 4.0;
    TacsScalar D66_stiff_sens =
        ((t + 0.5 * st) * (t + 0.5 * st) * (Qbar_panel[5] + Qbar_stiff[5])) /
        2.0;

    TacsScalar D3_sens =
        (D3 * D3 / sp) *
        ((sp - sh * sf_fraction) / (D66_panel * D66_panel) * D66_panel_sens +
         (sh * sf_fraction) / (D66_stiff * D66_stiff) * D66_stiff_sens);

    *sNxy_global_crit = calcCriticalShearLoadSens(D1, D2, D3, Lx, D1_sens,
                                                  D2_sens, D3_sens, 0.0);

    // Compute the remaining critical loads
    TacsScalar I = (t * t * t) / 12.0;
    TacsScalar D11 = I * Qbar_panel[0];
    TacsScalar D12 = I * Qbar_panel[1];
    TacsScalar D22 = I * Qbar_panel[3];
    TacsScalar D66 = I * Qbar_panel[5];

    TacsScalar I_sens = (t * t) / 4.0;
    TacsScalar D11_sens = I_sens * Qbar_panel[0];
    TacsScalar D12_sens = I_sens * Qbar_panel[1];
    TacsScalar D22_sens = I_sens * Qbar_panel[3];
    TacsScalar D66_sens = I_sens * Qbar_panel[5];

    TacsScalar sqrtD11D22 = sqrt(D11 * D22);

    *sNx_panel_crit = (2.0 * M_PI * M_PI / (sp * sp)) *
                      (0.5 * (D11 * D22_sens + D11_sens * D22) / sqrtD11D22 +
                       D12_sens + 2.0 * D66_sens);

    *sNxy_panel_crit =
        calcCriticalShearLoadSens(D11, D22, D12 + 2.0 * D66, sp, D11_sens,
                                  D22_sens, D12_sens + 2.0 * D66_sens, 0.0);

    *sNx_weak_stiff_crit = 0.0;
    *sNx_strong_stiff_crit = 0.0;
  } else if (dvNum == sp_num) {
    // Calculate the derivative w.r.t. the stiffener spacing
    TacsScalar Ep =
        Qbar_panel[0] - (Qbar_panel[1] * Qbar_panel[1]) / Qbar_panel[3];
    TacsScalar Es =
        Qbar_stiff[0] - (Qbar_stiff[1] * Qbar_stiff[1]) / Qbar_stiff[3];

    // The modulus-weighted area:
    TacsScalar An = (Ep * sp * t + Es * sh * st + Es * sh * st * sf_fraction);
    TacsScalar An_sens = Ep * t;

    // The modulus-weighted centroid
    TacsScalar Cn = 0.5 * Es * sh * sh * st;

    // The modulus-weighted centroid
    TacsScalar zn = Cn / An;
    TacsScalar zn_sens = -Cn / (An * An) * An_sens;

    // The bending stiffness
    TacsScalar EI = (zn * zn * (t * sp * Ep + st * sh * sf_fraction * Es) +
                     Es * (st * (sh * sh * sh) / 12.0 +
                           st * sh * (zn - 0.5 * sh) * (zn - 0.5 * sh)));
    TacsScalar EI_sens =
        (2.0 * zn * zn_sens * (t * sp * Ep + st * sh * sf_fraction * Es) +
         Es * (2.0 * st * sh * (zn - 0.5 * sh) * zn_sens) + zn * zn * (t * Ep));

    // Compute the unit EI stiffness for the buckling
    TacsScalar D1 = EI / sp;
    TacsScalar D1_sens = EI_sens / sp - EI / (sp * sp);
    *sNx_global_crit = (M_PI * M_PI * D1_sens) / (Lx * Lx);

    // Compute the required data for the shear loading
    // The transverse bending stiffness
    TacsScalar D12_panel = (t * t * t * Qbar_panel[1]) / 12.0;
    TacsScalar D12_stiff = ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
                            (Qbar_panel[1] + Qbar_stiff[1])) /
                           6.0;

    TacsScalar D2 = sp / ((sp - sh * sf_fraction) / D12_panel +
                          (sh * sf_fraction) / D12_stiff);
    TacsScalar D2_sens = (D2 / sp) - ((D2 * D2) / sp) * (1.0 / D12_panel);

    // The twisting stiffness
    TacsScalar D66_panel = (t * t * t * Qbar_panel[5]) / 12.0;
    TacsScalar D66_stiff = ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
                            (Qbar_panel[5] + Qbar_stiff[5])) /
                           6.0;

    TacsScalar D3 = sp / ((sp - sh * sf_fraction) / D66_panel +
                          (sh * sf_fraction) / D66_stiff);
    TacsScalar D3_sens = (D3 / sp) - ((D3 * D3) / sp) * (1.0 / D66_panel);

    *sNxy_global_crit = calcCriticalShearLoadSens(D1, D2, D3, Lx, D1_sens,
                                                  D2_sens, D3_sens, 0.0);

    // Compute the critical buckling loads for the panel
    TacsScalar I = (t * t * t) / 12.0;
    TacsScalar D11 = I * Qbar_panel[0];
    TacsScalar D12 = I * Qbar_panel[1];
    TacsScalar D22 = I * Qbar_panel[3];
    TacsScalar D66 = I * Qbar_panel[5];

    *sNx_panel_crit = (-4.0 * M_PI * M_PI / (sp * sp * sp)) *
                      (sqrt(D11 * D22) + D12 + 2.0 * D66);
    *sNxy_panel_crit = calcCriticalShearLoadSens(D11, D22, D12 + 2.0 * D66, sp,
                                                 0.0, 0.0, 0.0, 1.0);

    *sNx_weak_stiff_crit = 0.0;
    *sNx_strong_stiff_crit = 0.0;
  } else if (dvNum == sh_num) {
    // Calculate the derivative w.r.t. stiffener height
    TacsScalar Ep =
        Qbar_panel[0] - (Qbar_panel[1] * Qbar_panel[1]) / Qbar_panel[3];
    TacsScalar Es =
        Qbar_stiff[0] - (Qbar_stiff[1] * Qbar_stiff[1]) / Qbar_stiff[3];

    // The modulus-weighted area:
    TacsScalar An = (Ep * sp * t + Es * sh * st * (1.0 + sf_fraction));
    TacsScalar An_sens = Es * st * (1.0 + sf_fraction);

    // The modulus-weighted centroid
    TacsScalar Cn = 0.5 * Es * sh * sh * st;
    TacsScalar Cn_sens = Es * sh * st;

    // The modulus-weighted centroid
    TacsScalar zn = Cn / An;
    TacsScalar zn_sens = (Cn_sens * An - Cn * An_sens) / (An * An);

    // The bending stiffness
    TacsScalar EI = (zn * zn * (t * sp * Ep + st * sh * sf_fraction * Es) +
                     Es * (st * (sh * sh * sh) / 12.0 +
                           st * sh * (zn - 0.5 * sh) * (zn - 0.5 * sh)));
    TacsScalar EI_sens =
        (2.0 * zn * zn_sens * (t * sp * Ep + st * sh * sf_fraction * Es) +
         zn * zn * st * sf_fraction * Es +
         Es * (st * (sh * sh) / 4.0 + st * (zn - 0.5 * sh) * (zn - 0.5 * sh) +
               2.0 * st * sh * (zn - 0.5 * sh) * (zn_sens - 0.5)));

    // Compute the unit EI stiffness for the buckling
    TacsScalar D1 = EI / sp;
    TacsScalar D1_sens = EI_sens / sp;
    *sNx_global_crit = (M_PI * M_PI * D1_sens) / (Lx * Lx);

    // Compute the required data for the shear loading
    // The transverse bending stiffness
    TacsScalar D12_panel = (t * t * t * Qbar_panel[1]) / 12.0;
    TacsScalar D12_stiff = ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
                            (Qbar_panel[1] + Qbar_stiff[1])) /
                           6.0;

    TacsScalar D2 = sp / ((sp - sh * sf_fraction) / D12_panel +
                          (sh * sf_fraction) / D12_stiff);
    TacsScalar D2_sens =
        -((D2 * D2) / sp) * (sf_fraction / D12_stiff - sf_fraction / D12_panel);

    // The twisting stiffness
    TacsScalar D66_panel = (t * t * t * Qbar_panel[5]) / 12.0;
    TacsScalar D66_stiff = ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
                            (Qbar_panel[5] + Qbar_stiff[5])) /
                           6.0;

    TacsScalar D3 = sp / ((sp - sh * sf_fraction) / D66_panel +
                          (sh * sf_fraction) / D66_stiff);
    TacsScalar D3_sens =
        -((D3 * D3) / sp) * (sf_fraction / D66_stiff - sf_fraction / D66_panel);

    *sNxy_global_crit = calcCriticalShearLoadSens(D1, D2, D3, Lx, D1_sens,
                                                  D2_sens, D3_sens, 0.0);

    *sNx_panel_crit = 0.0;
    *sNxy_panel_crit = 0.0;

    TacsScalar Is_weak = st * st * st / 12.0;
    TacsScalar Is_strong = st * sh * sh / 4.0;

    *sNx_weak_stiff_crit = (M_PI * M_PI * Es * Is_weak) / (Lx * Lx);
    *sNx_strong_stiff_crit = (M_PI * M_PI * Es * Is_strong) / (Lx * Lx);
  } else if (dvNum == st_num) {
    // Compute the derivative w.r.t. the stiffener thickness
    TacsScalar Ep =
        Qbar_panel[0] - (Qbar_panel[1] * Qbar_panel[1]) / Qbar_panel[3];
    TacsScalar Es =
        Qbar_stiff[0] - (Qbar_stiff[1] * Qbar_stiff[1]) / Qbar_stiff[3];

    // The modulus-weighted area:
    TacsScalar An = (Ep * sp * t + Es * sh * st * (1.0 + sf_fraction));
    TacsScalar An_sens = Es * sh * (1.0 + sf_fraction);

    // The modulus-weighted centroid
    TacsScalar Cn = 0.5 * Es * sh * sh * st;
    TacsScalar Cn_sens = 0.5 * Es * sh * sh;

    // The modulus-weighted centroid
    TacsScalar zn = Cn / An;
    TacsScalar zn_sens = (Cn_sens * An - Cn * An_sens) / (An * An);

    // The bending stiffness
    TacsScalar EI = (zn * zn * (t * sp * Ep + st * sh * sf_fraction * Es) +
                     Es * (st * (sh * sh * sh) / 12.0 +
                           st * sh * (zn - 0.5 * sh) * (zn - 0.5 * sh)));
    TacsScalar EI_sens =
        (2.0 * zn * zn_sens * (t * sp * Ep + st * sh * sf_fraction * Es) +
         zn * zn * sh * sf_fraction * Es +
         Es * ((sh * sh * sh) / 12.0 + sh * (zn - 0.5 * sh) * (zn - 0.5 * sh) +
               2.0 * st * sh * (zn - 0.5 * sh) * zn_sens));

    // Compute the unit EI stiffness for the buckling
    TacsScalar D1 = EI / sp;
    TacsScalar D1_sens = EI_sens / sp;
    *sNx_global_crit = (M_PI * M_PI * D1_sens) / (Lx * Lx);

    // Compute the required data for the shear loading
    // The transverse bending stiffness
    TacsScalar D12_panel = (t * t * t * Qbar_panel[1]) / 12.0;
    TacsScalar D12_stiff = ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
                            (Qbar_panel[1] + Qbar_stiff[1])) /
                           6.0;
    TacsScalar D12_stiff_sens =
        1.5 *
        ((t + 0.5 * st) * (t + 0.5 * st) * (Qbar_panel[1] + Qbar_stiff[1])) /
        6.0;

    TacsScalar D2 = sp / ((sp - sh * sf_fraction) / D12_panel +
                          (sh * sf_fraction) / D12_stiff);
    TacsScalar D2_sens =
        (D2 * D2) / sp *
        ((sh * sf_fraction) / (D12_stiff * D12_stiff) * D12_stiff_sens);

    // The twisting stiffness
    TacsScalar D66_panel = (t * t * t * Qbar_panel[5]) / 12.0;
    TacsScalar D66_stiff = ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
                            (Qbar_panel[5] + Qbar_stiff[5])) /
                           6.0;
    TacsScalar D66_stiff_sens =
        1.5 *
        ((t + 0.5 * st) * (t + 0.5 * st) * (Qbar_panel[5] + Qbar_stiff[5])) /
        6.0;

    TacsScalar D3 = sp / ((sp - sh * sf_fraction) / D66_panel +
                          (sh * sf_fraction) / D66_stiff);
    TacsScalar D3_sens =
        (D3 * D3) / sp *
        ((sh * sf_fraction) / (D66_stiff * D66_stiff) * D66_stiff_sens);

    *sNxy_global_crit = calcCriticalShearLoadSens(D1, D2, D3, Lx, D1_sens,
                                                  D2_sens, D3_sens, 0.0);

    *sNx_panel_crit = 0.0;
    *sNxy_panel_crit = 0.0;

    TacsScalar Is_weak = sh * st * st / 4.0;
    TacsScalar Is_strong = sh * sh * sh / 12.0;

    *sNx_weak_stiff_crit = (M_PI * M_PI * Es * Is_weak) / (Lx * Lx);
    *sNx_strong_stiff_crit = (M_PI * M_PI * Es * Is_strong) / (Lx * Lx);
  } else if (dvNum == Lx_num) {
    // Compute the neutral axis location:
    TacsScalar Ep =
        Qbar_panel[0] - (Qbar_panel[1] * Qbar_panel[1]) / Qbar_panel[3];
    TacsScalar Es =
        Qbar_stiff[0] - (Qbar_stiff[1] * Qbar_stiff[1]) / Qbar_stiff[3];

    // The modulus-weighted area:
    TacsScalar An = (Ep * sp * t + Es * sh * st * (1.0 + sf_fraction));

    // The modulus-weighted centroid
    TacsScalar Cn = 0.5 * Es * sh * sh * st;

    // The modulus-weighted centroid
    TacsScalar zn = Cn / An;

    // The bending stiffness
    TacsScalar EI = (zn * zn * (t * sp * Ep + st * sh * sf_fraction * Es) +
                     Es * (st * (sh * sh * sh) / 12.0 +
                           st * sh * (zn - 0.5 * sh) * (zn - 0.5 * sh)));

    // Compute the unit EI stiffness for the buckling
    TacsScalar D1 = EI / sp;
    *sNx_global_crit = -2.0 * (M_PI * M_PI * D1) / (Lx * Lx * Lx);

    // Compute the required data for the shear loading
    // The transverse bending stiffness
    TacsScalar D12_panel = (t * t * t * Qbar_panel[1]) / 12.0;
    TacsScalar D12_stiff = ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
                            (Qbar_panel[1] + Qbar_stiff[1])) /
                           6.0;

    TacsScalar D2 = sp / ((sp - sh * sf_fraction) / D12_panel +
                          (sh * sf_fraction) / D12_stiff);

    // The twisting stiffness
    TacsScalar D66_panel = (t * t * t * Qbar_panel[5]) / 12.0;
    TacsScalar D66_stiff = ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
                            (Qbar_panel[5] + Qbar_stiff[5])) /
                           6.0;

    TacsScalar D3 = sp / ((sp - sh * sf_fraction) / D66_panel +
                          (sh * sf_fraction) / D66_stiff);

    *sNxy_global_crit =
        calcCriticalShearLoadSens(D1, D2, D3, Lx, 0.0, 0.0, 0.0, 1.0);

    *sNx_panel_crit = 0.0;
    *sNxy_panel_crit = 0.0;

    TacsScalar Is_weak = sh * st * st * st / 12.0;
    TacsScalar Is_strong = st * sh * sh * sh / 12.0;

    *sNx_weak_stiff_crit = -2.0 * (M_PI * M_PI * Es * Is_weak) / (Lx * Lx * Lx);
    *sNx_strong_stiff_crit =
        -2.0 * (M_PI * M_PI * Es * Is_strong) / (Lx * Lx * Lx);
  } else if (dvNum >= 0) {
    TacsScalar sQbar_panel[6], sQbar_stiff[6];
    sQbar_panel[0] = sQbar_panel[1] = sQbar_panel[2] = 0.0;
    sQbar_panel[3] = sQbar_panel[4] = sQbar_panel[5] = 0.0;

    sQbar_stiff[0] = sQbar_stiff[1] = sQbar_stiff[2] = 0.0;
    sQbar_stiff[3] = sQbar_stiff[4] = sQbar_stiff[5] = 0.0;

    for (int k = 0; k < NUM_PLIES; k++) {
      if (pf_nums[k] == dvNum) {
        TacsScalar Qp[6];
        ortho_ply->calculateQbar(Qp, ply_angles[k]);

        for (int j = 0; j < 6; j++) {
          sQbar_panel[j] += Qp[j];
        }
      }
    }

    for (int k = 0; k < NUM_PLIES; k++) {
      if (stiff_pf_nums[k] == dvNum) {
        TacsScalar Qp[6];
        ortho_ply->calculateQbar(Qp, stiff_ply_angles[k]);

        for (int j = 0; j < 6; j++) {
          sQbar_stiff[j] += Qp[j];
        }
      }
    }

    // Compute the neutral axis location:
    TacsScalar Ep =
        Qbar_panel[0] - (Qbar_panel[1] * Qbar_panel[1]) / Qbar_panel[3];
    TacsScalar Ep_sens =
        (sQbar_panel[0] -
         2.0 * (Qbar_panel[1] * sQbar_panel[1]) / Qbar_panel[3] +
         sQbar_panel[3] * (Qbar_panel[1] * Qbar_panel[1]) /
             (Qbar_panel[3] * Qbar_panel[3]));

    TacsScalar Es =
        Qbar_stiff[0] - (Qbar_stiff[1] * Qbar_stiff[1]) / Qbar_stiff[3];
    TacsScalar Es_sens =
        (sQbar_stiff[0] -
         2.0 * (Qbar_stiff[1] * sQbar_stiff[1]) / Qbar_stiff[3] +
         sQbar_stiff[3] * (Qbar_stiff[1] * Qbar_stiff[1]) /
             (Qbar_stiff[3] * Qbar_stiff[3]));

    // The modulus-weighted area:
    TacsScalar An = (Ep * sp * t + Es * sh * st * (1.0 + sf_fraction));
    TacsScalar An_sens =
        (Ep_sens * sp * t + Es_sens * sh * st * (1.0 + sf_fraction));

    // The modulus-weighted centroid
    TacsScalar Cn = 0.5 * Es * sh * sh * st;
    TacsScalar Cn_sens = 0.5 * Es_sens * sh * sh * st;

    // The modulus-weighted centroid
    TacsScalar zn = Cn / An;
    TacsScalar zn_sens = (Cn_sens * An - Cn * An_sens) / (An * An);

    // The bending stiffness
    TacsScalar EI = (zn * zn * (t * sp * Ep + st * sh * sf_fraction * Es) +
                     Es * (st * (sh * sh * sh) / 12.0 +
                           st * sh * (zn - 0.5 * sh) * (zn - 0.5 * sh)));

    TacsScalar EI_sens =
        (zn * zn * (t * sp * Ep_sens + st * sh * sf_fraction * Es_sens) +
         Es_sens * (st * (sh * sh * sh) / 12.0 +
                    st * sh * (zn - 0.5 * sh) * (zn - 0.5 * sh)) +
         2.0 * zn * zn_sens * (t * sp * Ep + st * sh * sf_fraction * Es) +
         Es * (2.0 * st * sh * (zn - 0.5 * sh) * zn_sens));

    // Compute the unit EI stiffness for the buckling
    TacsScalar D1 = EI / sp;
    TacsScalar D1_sens = EI_sens / sp;
    *sNx_global_crit = (M_PI * M_PI * D1_sens) / (Lx * Lx);

    // Compute the required data for the shear loading
    // The transverse bending stiffness
    TacsScalar D12_panel = (t * t * t * Qbar_panel[1]) / 12.0;
    TacsScalar D12_stiff = ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
                            (Qbar_panel[1] + Qbar_stiff[1])) /
                           6.0;
    TacsScalar D2 = sp / ((sp - sh * sf_fraction) / D12_panel +
                          (sh * sf_fraction) / D12_stiff);

    TacsScalar D12_panel_sens = (t * t * t * sQbar_panel[1]) / 12.0;
    TacsScalar D12_stiff_sens =
        ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
         (sQbar_panel[1] + sQbar_stiff[1])) /
        6.0;
    TacsScalar D2_sens =
        (D2 * D2 / sp) *
        ((sp - sh * sf_fraction) / (D12_panel * D12_panel) * D12_panel_sens +
         (sh * sf_fraction) / (D12_stiff * D12_stiff) * D12_stiff_sens);

    // The twisting stiffness
    TacsScalar D66_panel = (t * t * t * Qbar_panel[5]) / 12.0;
    TacsScalar D66_stiff = ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
                            (Qbar_panel[5] + Qbar_stiff[5])) /
                           6.0;

    TacsScalar D3 = sp / ((sp - sh * sf_fraction) / D66_panel +
                          (sh * sf_fraction) / D66_stiff);

    TacsScalar D66_panel_sens = (t * t * t * sQbar_panel[5]) / 12.0;
    TacsScalar D66_stiff_sens =
        ((t + 0.5 * st) * (t + 0.5 * st) * (t + 0.5 * st) *
         (sQbar_panel[5] + sQbar_stiff[5])) /
        6.0;

    TacsScalar D3_sens =
        (D3 * D3 / sp) *
        ((sp - sh * sf_fraction) / (D66_panel * D66_panel) * D66_panel_sens +
         (sh * sf_fraction) / (D66_stiff * D66_stiff) * D66_stiff_sens);

    *sNxy_global_crit = calcCriticalShearLoadSens(D1, D2, D3, Lx, D1_sens,
                                                  D2_sens, D3_sens, 0.0);

    // Compute the critical buckling loads for the panel
    // -------------------------------------------------
    TacsScalar I = (t * t * t) / 12.0;
    TacsScalar D11 = I * Qbar_panel[0];
    TacsScalar D12 = I * Qbar_panel[1];
    TacsScalar D22 = I * Qbar_panel[3];
    TacsScalar D66 = I * Qbar_panel[5];

    TacsScalar D11_sens = I * sQbar_panel[0];
    TacsScalar D12_sens = I * sQbar_panel[1];
    TacsScalar D22_sens = I * sQbar_panel[3];
    TacsScalar D66_sens = I * sQbar_panel[5];

    TacsScalar sqrtD11D22 = sqrt(D11 * D22);

    *sNx_panel_crit = (2.0 * M_PI * M_PI / (sp * sp)) *
                      (0.5 * (D11 * D22_sens + D11_sens * D22) / sqrtD11D22 +
                       D12_sens + 2.0 * D66_sens);

    *sNxy_panel_crit =
        calcCriticalShearLoadSens(D11, D22, D12 + 2.0 * D66, sp, D11_sens,
                                  D22_sens, D12_sens + 2.0 * D66_sens, 0.0);

    // Compute the critical loads for the stiffener
    // --------------------------------------------
    TacsScalar Is_weak = sh * st * st * st / 12.0;
    TacsScalar Is_strong = st * sh * sh * sh / 12.0;

    *sNx_weak_stiff_crit = (M_PI * M_PI * Es_sens * Is_weak) / (Lx * Lx);
    *sNx_strong_stiff_crit = (M_PI * M_PI * Es_sens * Is_strong) / (Lx * Lx);
  }
}

/*
  Print information to stdout about this stiffener class.

  This is useful for debugging purposes. This prints the design
  variable values and the critical buckling loads for the stiffener
  class.
*/
void bladeFSDTStiffness::printBucklingInfo() {
  double gpt[2] = {0.0, 0.0};
  TacsScalar mass[3];
  getPointwiseMass(gpt, mass);

  printf("%s properties\n", constName);
  printf("Lx panel dimension:  %15.6f\n", TacsRealPart(Lx));
  printf("stiffener pitch:     %15.6f\n", TacsRealPart(sp));
  printf("stiffener height:    %15.6f\n", TacsRealPart(sh));
  printf("stiffener thickness: %15.6f\n", TacsRealPart(st));
  printf("panel thickness:     %15.6f\n", TacsRealPart(t));
  printf("mass per unit area:  %15.6f\n", TacsRealPart(mass[0]));
  printf("global Nx,cr:        %15.6f\n", TacsRealPart(Nx_global_crit));
  printf("global Nxy,cr:       %15.6f\n", TacsRealPart(Nxy_global_crit));
  printf("panel Nx,cr:         %15.6f\n", TacsRealPart(Nx_panel_crit));
  printf("panel Nxy,cr:        %15.6f\n", TacsRealPart(Nxy_panel_crit));
  printf("stiffener Nx,cr      %15.6f\n", TacsRealPart(Nx_weak_stiff_crit));
  printf("stiffener Nx,cr      %15.6f\n", TacsRealPart(Nx_strong_stiff_crit));
}

/*
  Compute all the failure criteria for each ply, for each of the
  top/bottom plies and the plies in the stiffeners.
*/
int bladeFSDTStiffness::calcFail(const TacsScalar strain[], TacsScalar fail[],
                                 TacsScalar* _max) {
  // Keep track of the maximum value of the failure criteria and the
  // total number of load numbers
  TacsScalar max = 0.0;
  int load_num = 0;

  TacsScalar e[3];
  // Compute the strain at the upper panel surface
  e[0] = strain[0] + 0.5 * t * strain[3];
  e[1] = strain[1] + 0.5 * t * strain[4];
  e[2] = strain[2] + 0.5 * t * strain[5];

  for (int k = 0; k < NUM_PLIES; k++) {
    fail[load_num] = ortho_ply->failure(ply_angles[k], e);
    if (TacsRealPart(fail[load_num]) > TacsRealPart(max)) {
      max = fail[load_num];
    }
    load_num++;
  }

  // Compute the strain at the lower panel sufrace
  e[0] = strain[0] - 0.5 * t * strain[3];
  e[1] = strain[1] - 0.5 * t * strain[4];
  e[2] = strain[2] - 0.5 * t * strain[5];

  for (int k = 0; k < NUM_PLIES; k++) {
    fail[load_num] = ortho_ply->failure(ply_angles[k], e);
    if (TacsRealPart(fail[load_num]) > TacsRealPart(max)) {
      max = fail[load_num];
    }
    load_num++;
  }

  // Compute the strain at the furthest fiber of the stiffener
  e[0] = strain[0] - sh * strain[3];
  e[1] = 0.0;
  e[2] = 0.0;

  // Compute the failure criteria for the stiffener
  for (int k = 0; k < NUM_PLIES; k++) {
    fail[load_num] = ortho_ply->failure(stiff_ply_angles[k], e);
    if (TacsRealPart(fail[load_num]) > TacsRealPart(max)) {
      max = fail[load_num];
    }
    load_num++;
  }

  *_max = max;

  return load_num;
}

/*
  Compute the derivative of each of the failure
*/
void bladeFSDTStiffness::calcFailStrainSens(const TacsScalar strain[],
                                            const TacsScalar weights[],
                                            TacsScalar sens[]) {
  sens[0] = sens[1] = sens[2] = sens[3] = 0.0;
  sens[4] = sens[5] = sens[6] = sens[7] = 0.0;

  int load_num = 0;

  for (int k = 0; k < NUM_PLIES; k++) {
    TacsScalar e[3], e_sens[3];
    // Compute the strain at the upper surface
    e[0] = strain[0] + 0.5 * t * strain[3];
    e[1] = strain[1] + 0.5 * t * strain[4];
    e[2] = strain[2] + 0.5 * t * strain[5];

    ortho_ply->failureStrainSens(e_sens, ply_angles[k], e);

    sens[0] += weights[load_num] * e_sens[0];
    sens[1] += weights[load_num] * e_sens[1];
    sens[2] += weights[load_num] * e_sens[2];

    sens[3] += 0.5 * t * weights[load_num] * e_sens[0];
    sens[4] += 0.5 * t * weights[load_num] * e_sens[1];
    sens[5] += 0.5 * t * weights[load_num] * e_sens[2];
    load_num++;
  }

  for (int k = 0; k < NUM_PLIES; k++) {
    TacsScalar e[3], e_sens[3];
    // Compute the strain at the lower sufrace
    e[0] = strain[0] - 0.5 * t * strain[3];
    e[1] = strain[1] - 0.5 * t * strain[4];
    e[2] = strain[2] - 0.5 * t * strain[5];

    ortho_ply->failureStrainSens(e_sens, ply_angles[k], e);

    sens[0] += weights[load_num] * e_sens[0];
    sens[1] += weights[load_num] * e_sens[1];
    sens[2] += weights[load_num] * e_sens[2];

    sens[3] -= 0.5 * t * weights[load_num] * e_sens[0];
    sens[4] -= 0.5 * t * weights[load_num] * e_sens[1];
    sens[5] -= 0.5 * t * weights[load_num] * e_sens[2];
    load_num++;
  }

  // Compute the failure criteria for the stiffener
  for (int k = 0; k < NUM_PLIES; k++) {
    TacsScalar e[3], e_sens[3];
    // Compute the strain at the upper surface
    e[0] = strain[0] - sh * strain[3];
    e[1] = 0.0;
    e[2] = 0.0;

    ortho_ply->failureStrainSens(e_sens, stiff_ply_angles[k], e);

    sens[0] += weights[load_num] * e_sens[0];
    sens[3] -= sh * weights[load_num] * e_sens[0];
    load_num++;
  }
}

TacsScalar bladeFSDTStiffness::calcFailDVSens(int dvNum,
                                              const TacsScalar strain[],
                                              const TacsScalar weights[]) {
  int load_num = 0;
  TacsScalar failSens = 0.0;

  if (dvNum == t_num) {
    for (int k = 0; k < NUM_PLIES; k++) {
      TacsScalar e[3], e_sens[3];
      // Compute the strain at the upper surface
      e[0] = strain[0] + 0.5 * t * strain[3];
      e[1] = strain[1] + 0.5 * t * strain[4];
      e[2] = strain[2] + 0.5 * t * strain[5];

      ortho_ply->failureStrainSens(e_sens, ply_angles[k], e);
      failSens += 0.5 * weights[load_num] *
                  (strain[3] * e_sens[0] + strain[4] * e_sens[1] +
                   strain[5] * e_sens[2]);
      load_num++;
    }

    for (int k = 0; k < NUM_PLIES; k++) {
      TacsScalar e[3], e_sens[3];
      // Compute the strain at the lower sufrace
      e[0] = strain[0] - 0.5 * t * strain[3];
      e[1] = strain[1] - 0.5 * t * strain[4];
      e[2] = strain[2] - 0.5 * t * strain[5];

      ortho_ply->failureStrainSens(e_sens, ply_angles[k], e);
      failSens -= 0.5 * weights[load_num] *
                  (strain[3] * e_sens[0] + strain[4] * e_sens[1] +
                   strain[5] * e_sens[2]);
      load_num++;
    }
  } else if (dvNum == sh_num) {
    load_num = 2 * NUM_PLIES;

    // Compute the failure criteria for the stiffener
    for (int k = 0; k < NUM_PLIES; k++) {
      TacsScalar e[3], e_sens[3];
      // Compute the strain at the upper surface
      e[0] = strain[0] - sh * strain[3];
      e[1] = 0.0;
      e[2] = 0.0;

      ortho_ply->failureStrainSens(e_sens, stiff_ply_angles[k], e);
      failSens -= weights[load_num] * strain[3] * e_sens[0];
      load_num++;
    }
  }

  return failSens;
}

/*
  Functions to determine the failure load in terms of the strain and
  derivatives of the failure load in terms of the strain components.
*/
void bladeFSDTStiffness::failure(const double pt[], const TacsScalar strain[],
                                 TacsScalar* fval) {
  TacsScalar fail[MAX_NUM_FAIL], max;
  int num_fail = calcFail(strain, fail, &max);
  TacsScalar ks_sum = 0.0;
  for (int k = 0; k < num_fail; k++) {
    ks_sum += exp(ks_weight * (fail[k] - max));
  }

  *fval = max + log(ks_sum) / ks_weight;
}

void bladeFSDTStiffness::failureStrainSens(const double pt[],
                                           const TacsScalar strain[],
                                           TacsScalar sens[]) {
  TacsScalar fail[MAX_NUM_FAIL], sens_weights[MAX_NUM_FAIL], max;
  int num_fail = calcFail(strain, fail, &max);
  TacsScalar ks_sum = 0.0;
  for (int k = 0; k < num_fail; k++) {
    ks_sum += exp(ks_weight * (fail[k] - max));
  }

  for (int k = 0; k < num_fail; k++) {
    sens_weights[k] = exp(ks_weight * (fail[k] - max)) / ks_sum;
  }

  calcFailStrainSens(strain, sens_weights, sens);
}

/*
  Functions to determine the derivative of the failure
  load w.r.t. the design variables
*/
void bladeFSDTStiffness::failureDVSens(int dvNum, const double pt[],
                                       const TacsScalar strain[],
                                       TacsScalar* failSens) {
  TacsScalar fail[MAX_NUM_FAIL], sens_weights[MAX_NUM_FAIL], max;
  int num_fail = calcFail(strain, fail, &max);
  TacsScalar ks_sum = 0.0;
  for (int k = 0; k < num_fail; k++) {
    ks_sum += exp(ks_weight * (fail[k] - max));
  }

  for (int k = 0; k < num_fail; k++) {
    sens_weights[k] = exp(ks_weight * (fail[k] - max)) / ks_sum;
  }

  *failSens = calcFailDVSens(dvNum, strain, sens_weights);
}

/*
  Add the derivative of the failure load times alpha
*/
void bladeFSDTStiffness::addFailureDVSens(const double pt[],
                                          const TacsScalar strain[],
                                          TacsScalar alpha, TacsScalar dvSens[],
                                          int dvLen) {
  TacsScalar fail[MAX_NUM_FAIL], sens_weights[MAX_NUM_FAIL], max;
  int num_fail = calcFail(strain, fail, &max);
  TacsScalar ks_sum = 0.0;
  for (int k = 0; k < num_fail; k++) {
    ks_sum += exp(ks_weight * (fail[k] - max));
  }

  // Compute the sensitivities of the arrays
  for (int k = 0; k < num_fail; k++) {
    sens_weights[k] = exp(ks_weight * (fail[k] - max)) / ks_sum;
  }

  // Add the scaled multiple of the derivatives
  if (Lx_num >= 0 && Lx_num < dvLen) {
    dvSens[Lx_num] += alpha * calcFailDVSens(Lx_num, strain, sens_weights);
  }
  if (sp_num >= 0 && sp_num < dvLen) {
    dvSens[sp_num] += alpha * calcFailDVSens(sp_num, strain, sens_weights);
  }
  if (sh_num >= 0 && sh_num < dvLen) {
    dvSens[sh_num] += alpha * calcFailDVSens(sh_num, strain, sens_weights);
  }
  if (st_num >= 0 && st_num < dvLen) {
    dvSens[st_num] += alpha * calcFailDVSens(st_num, strain, sens_weights);
  }
  if (t_num >= 0 && t_num < dvLen) {
    dvSens[t_num] += alpha * calcFailDVSens(t_num, strain, sens_weights);
  }

  // Add the derivatives from the ply fraction variables
  for (int k = 0; k < num_pf_unique; k++) {
    if (pf_unique[k] < dvLen) {
      dvSens[pf_unique[k]] +=
          alpha * calcFailDVSens(pf_unique[k], strain, sens_weights);
    }
  }
  for (int k = 0; k < num_stiff_pf_unique; k++) {
    if (stiff_pf_unique[k] < dvLen) {
      dvSens[stiff_pf_unique[k]] +=
          alpha * calcFailDVSens(stiff_pf_unique[k], strain, sens_weights);
    }
  }
}

/*
  Compute the buckling criteria for every component of the stiffened
  panel.
*/
int bladeFSDTStiffness::calcBuckling(const TacsScalar strain[],
                                     TacsScalar bval[], TacsScalar* _max) {
  TacsScalar max = 0.0;
  int load_num = 0;

  // Compute the buckling contributions to the failure load.
  // These use the interaction criterion:
  // bval = (Nxy/Nxy,crit)^2 + Nx/Nx,crit <= 1
  // -------------------------------------------------------
  // Calculate the global buckling criteria
  TacsScalar stress[8];
  calculateStress(A, B, D, As, strain, stress);

  TacsScalar Nx = -stress[0];
  TacsScalar Nxy = stress[2];

  bval[load_num] = ((Nx / Nx_global_crit) +
                    (Nxy / Nxy_global_crit) * (Nxy / Nxy_global_crit));
  max = bval[load_num];
  load_num++;

  // Calculate the panel buckling criteria
  // This is for the inter-stiffener portion of the panel
  Nx = -t * (Qbar_panel[0] * strain[0] + Qbar_panel[1] * strain[1]);
  Nxy = t * Qbar_panel[5] * strain[2];

  bval[load_num] =
      ((Nx / Nx_panel_crit) + (Nxy / Nxy_panel_crit) * (Nxy / Nxy_panel_crit));
  if (TacsRealPart(bval[load_num]) > TacsRealPart(max)) {
    max = bval[load_num];
  }
  load_num++;

  *_max = max;

  return load_num;
}

void bladeFSDTStiffness::calcBucklingStrainSens(const TacsScalar strain[],
                                                const TacsScalar weights[],
                                                TacsScalar sens[]) {
  sens[0] = sens[1] = sens[2] = sens[3] = 0.0;
  sens[4] = sens[5] = sens[6] = sens[7] = 0.0;

  int load_num = 0;

  // Compute the buckling contributions to the failure load.
  // These use the interaction criterion:
  // (Nxy/Nxy,crit)^2 + Nx/Nx,crit <= 1
  // -------------------------------------------------------

  // Calculate the global buckling criteria
  TacsScalar stress[8], s_sens[8];
  calculateStress(A, B, D, As, strain, stress);

  s_sens[0] = -weights[load_num] / Nx_global_crit;
  s_sens[2] =
      2.0 * weights[load_num] * stress[2] / (Nxy_global_crit * Nxy_global_crit);
  s_sens[1] = s_sens[3] = s_sens[4] = s_sens[5] = s_sens[6] = s_sens[7] = 0.0;

  calculateStress(A, B, D, As, s_sens, stress);

  for (int j = 0; j < 8; j++) {
    sens[j] += stress[j];
  }
  load_num++;

  // Calculate the panel buckling criteria
  // This is for the inter-stiffener portion of the panel
  TacsScalar Nxy = t * Qbar_panel[5] * strain[2];

  sens[0] -= (weights[load_num] * t * Qbar_panel[0]) / Nx_panel_crit;
  sens[1] -= (weights[load_num] * t * Qbar_panel[1]) / Nx_panel_crit;
  sens[2] += 2.0 * (weights[load_num] * t * Qbar_panel[5] * Nxy) /
             (Nxy_panel_crit * Nxy_panel_crit);
  load_num++;
}

TacsScalar bladeFSDTStiffness::calcBucklingDVSens(int dvNum,
                                                  const TacsScalar strain[],
                                                  const TacsScalar weights[]) {
  int load_num = 0;
  TacsScalar bvalSens = 0.0;

  TacsScalar sNx_global_crit, sNxy_global_crit;
  TacsScalar sNx_panel_crit, sNxy_panel_crit;
  TacsScalar sNx_weak_stiff_crit, sNx_strong_stiff_crit;
  calcBucklingPropsDVSens(dvNum, &sNx_global_crit, &sNxy_global_crit,
                          &sNx_panel_crit, &sNxy_panel_crit,
                          &sNx_weak_stiff_crit, &sNx_strong_stiff_crit);

  TacsScalar sQbar_panel[6], sQbar_stiff[6];
  sQbar_panel[0] = sQbar_panel[1] = sQbar_panel[2] = 0.0;
  sQbar_panel[3] = sQbar_panel[4] = sQbar_panel[5] = 0.0;

  sQbar_stiff[0] = sQbar_stiff[1] = sQbar_stiff[2] = 0.0;
  sQbar_stiff[3] = sQbar_stiff[4] = sQbar_stiff[5] = 0.0;

  for (int k = 0; k < NUM_PLIES; k++) {
    if (pf_nums[k] == dvNum) {
      TacsScalar Qp[6];
      ortho_ply->calculateQbar(Qp, ply_angles[k]);

      for (int j = 0; j < 6; j++) {
        sQbar_panel[j] += Qp[j];
      }
    }
  }

  for (int k = 0; k < NUM_PLIES; k++) {
    if (stiff_pf_nums[k] == dvNum) {
      TacsScalar Qp[6];
      ortho_ply->calculateQbar(Qp, stiff_ply_angles[k]);

      for (int j = 0; j < 6; j++) {
        sQbar_stiff[j] += Qp[j];
      }
    }
  }

  // Calculate the global buckling criteria
  TacsScalar stress[8], stress_sens[8];
  calculateStress(A, B, D, As, strain, stress);

  const double gpt[2] = {0.0, 0.0};
  calculateStressDVSens(gpt, strain, stress_sens, dvNum);

  TacsScalar Nx = 0.0, Nxy = 0.0, sNx = 0.0, sNxy = 0.0;
  Nx = -stress[0];
  sNx = -stress_sens[0];

  Nxy = stress[2];
  sNxy = stress_sens[2];

  bvalSens += weights[load_num] *
              ((sNx / Nx_global_crit) -
               (Nx / (Nx_global_crit * Nx_global_crit)) * sNx_global_crit +
               2.0 * (Nxy * sNxy) / (Nxy_global_crit * Nxy_global_crit) -
               2.0 * (Nxy * Nxy) /
                   (Nxy_global_crit * Nxy_global_crit * Nxy_global_crit) *
                   sNxy_global_crit);
  load_num++;

  // Add the panel buckling failure contribution
  Nx = -t * (Qbar_panel[0] * strain[0] + Qbar_panel[1] * strain[1]);
  Nxy = t * Qbar_panel[5] * strain[2];

  if (dvNum == t_num) {
    sNxy = (t * sQbar_panel[5] + Qbar_panel[5]) * strain[2];
  } else {
    sNxy = t * sQbar_panel[5] * strain[2];
  }

  // Compute the axial contribution
  if (dvNum == t_num) {
    sNx = -(t * (sQbar_panel[0] * strain[0] + sQbar_panel[1] * strain[1]) +
            (Qbar_panel[0] * strain[0] + Qbar_panel[1] * strain[1]));
  } else {
    sNx = -t * (sQbar_panel[0] * strain[0] + sQbar_panel[1] * strain[1]);
  }

  bvalSens += weights[load_num] *
              (sNx / Nx_panel_crit -
               (sNx_panel_crit * Nx) / (Nx_panel_crit * Nx_panel_crit) +
               2.0 * (Nxy * sNxy) / (Nxy_panel_crit * Nxy_panel_crit) -
               2.0 * (sNxy_panel_crit * Nxy * Nxy) /
                   (Nxy_panel_crit * Nxy_panel_crit * Nxy_panel_crit));
  load_num++;

  // Calculate the contribution from the stiffeners
  Nx = -st * sh *
       (Qbar_stiff[0] * (strain[0] - 0.5 * sh * strain[3]) +
        Qbar_stiff[1] * (strain[1] - 0.5 * sh * strain[4]));

  sNx = -st * sh *
        (sQbar_stiff[0] * (strain[0] - 0.5 * sh * strain[3]) +
         sQbar_stiff[1] * (strain[1] - 0.5 * sh * strain[4]));

  if (dvNum == sh_num) {
    sNx = -st * (Qbar_stiff[0] * (strain[0] - sh * strain[3]) +
                 Qbar_stiff[1] * (strain[1] - sh * strain[4]));
  } else if (dvNum == st_num) {
    sNx = -sh * (Qbar_stiff[0] * (strain[0] - 0.5 * sh * strain[3]) +
                 Qbar_stiff[1] * (strain[1] - 0.5 * sh * strain[4]));
  }

  return bvalSens;
}

/*
  Functions to determine the critical buckling function in terms of the
  strain.
*/
void bladeFSDTStiffness::buckling(const TacsScalar strain[], TacsScalar* bval) {
  TacsScalar bfuncs[NUM_BUCKLING], max;
  int num_buckling = calcBuckling(strain, bfuncs, &max);
  TacsScalar ks_sum = 0.0;
  for (int k = 0; k < num_buckling; k++) {
    ks_sum += exp(ks_weight * (bfuncs[k] - max));
  }

  *bval = max + log(ks_sum) / ks_weight;
}

void bladeFSDTStiffness::bucklingStrainSens(const TacsScalar strain[],
                                            TacsScalar sens[]) {
  TacsScalar bfuncs[NUM_BUCKLING], sens_weights[NUM_BUCKLING], max;
  int num_buckling = calcBuckling(strain, bfuncs, &max);
  TacsScalar ks_sum = 0.0;
  for (int k = 0; k < num_buckling; k++) {
    ks_sum += exp(ks_weight * (bfuncs[k] - max));
  }

  for (int k = 0; k < num_buckling; k++) {
    sens_weights[k] = exp(ks_weight * (bfuncs[k] - max)) / ks_sum;
  }

  calcBucklingStrainSens(strain, sens_weights, sens);
}

/*
  Functions to determine the derivative of the buckling criteria
  load w.r.t. the design variables
*/
void bladeFSDTStiffness::bucklingDVSens(int dvNum, const TacsScalar strain[],
                                        TacsScalar* failSens) {
  TacsScalar bfuncs[NUM_BUCKLING], sens_weights[NUM_BUCKLING], max;
  int num_buckling = calcBuckling(strain, bfuncs, &max);
  TacsScalar ks_sum = 0.0;
  for (int k = 0; k < num_buckling; k++) {
    ks_sum += exp(ks_weight * (bfuncs[k] - max));
  }

  for (int k = 0; k < num_buckling; k++) {
    sens_weights[k] = exp(ks_weight * (bfuncs[k] - max)) / ks_sum;
  }

  *failSens = calcBucklingDVSens(dvNum, strain, sens_weights);
}

/*
  Add the derivative of the buckling load times alpha to the array
  dvSens
*/
void bladeFSDTStiffness::addBucklingDVSens(const TacsScalar strain[],
                                           TacsScalar alpha,
                                           TacsScalar dvSens[], int dvLen) {
  TacsScalar bfuncs[NUM_BUCKLING], sens_weights[NUM_BUCKLING], max;
  int num_buckling = calcBuckling(strain, bfuncs, &max);
  TacsScalar ks_sum = 0.0;
  for (int k = 0; k < num_buckling; k++) {
    ks_sum += exp(ks_weight * (bfuncs[k] - max));
  }

  for (int k = 0; k < num_buckling; k++) {
    sens_weights[k] = exp(ks_weight * (bfuncs[k] - max)) / ks_sum;
  }

  // Add the scaled multiple of the derivatives
  if (Lx_num >= 0 && Lx_num < dvLen) {
    dvSens[Lx_num] += alpha * calcBucklingDVSens(Lx_num, strain, sens_weights);
  }
  if (sp_num >= 0 && sp_num < dvLen) {
    dvSens[sp_num] += alpha * calcBucklingDVSens(sp_num, strain, sens_weights);
  }
  if (sh_num >= 0 && sh_num < dvLen) {
    dvSens[sh_num] += alpha * calcBucklingDVSens(sh_num, strain, sens_weights);
  }
  if (st_num >= 0 && st_num < dvLen) {
    dvSens[st_num] += alpha * calcBucklingDVSens(st_num, strain, sens_weights);
  }
  if (t_num >= 0 && t_num < dvLen) {
    dvSens[t_num] += alpha * calcBucklingDVSens(t_num, strain, sens_weights);
  }

  // Add the derivatives from the ply fraction variables
  for (int k = 0; k < num_pf_unique; k++) {
    if (pf_unique[k] < dvLen) {
      dvSens[pf_unique[k]] +=
          alpha * calcBucklingDVSens(pf_unique[k], strain, sens_weights);
    }
  }
  for (int k = 0; k < num_stiff_pf_unique; k++) {
    if (stiff_pf_unique[k] < dvLen) {
      dvSens[stiff_pf_unique[k]] +=
          alpha * calcBucklingDVSens(stiff_pf_unique[k], strain, sens_weights);
    }
  }
}

const char* bladeFSDTStiffness::constitutiveName() const { return constName; }

// Retrieve design variable information ...
TacsScalar bladeFSDTStiffness::getDVOutputValue(int dv_index,
                                                const double gpt[]) {
  if (dv_index == 0) {
    return t + sh * st * (1.0 + sf_fraction) / sp;
  } else if (dv_index - 1 < NUM_PLIES) {
    return pf[dv_index - 1];
  }

  return 0.0;
}

/*
  Get the maximum number of design variables ownable by this constitutive class
*/
int bladeFSDTStiffness::getNumDesignVars() { return 5 + 2 * NUM_PLIES; }

int bladeFSDTStiffness::getDesignVarNums(int dvNums[], int dvLen) {
  if (Lx_num >= 0) {
    dvNums[0] = Lx_num;
  } else {
    dvNums[0] = -1;
  }

  if (sp_num >= 0) {
    dvNums[1] = sp_num;
  } else {
    dvNums[1] = -1;
  }

  if (sh_num >= 0) {
    dvNums[2] = sh_num;
  } else {
    dvNums[2] = -1;
  }

  if (st_num >= 0) {
    dvNums[3] = st_num;
  } else {
    dvNums[3] = -1;
  }

  if (t_num >= 0) {
    dvNums[4] = t_num;
  } else {
    dvNums[4] = -1;
  }

  for (int i = 0; i < NUM_PLIES; i++) {
    if (pf_nums[i] >= 0) {
      dvNums[i + 5] = pf_nums[i];
    } else {
      dvNums[i + 5] = -1;
    }

    if (stiff_pf_nums[i] >= 0) {
      dvNums[i + NUM_PLIES + 5] = stiff_pf_nums[i];
    } else {
      dvNums[i + NUM_PLIES + 5] = -1;
    }
  }
  return 1;
}
