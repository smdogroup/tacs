/*
=============================================================================
Testing the BladeStiffenedShell constitutive model
=============================================================================
@File    :   bladeStiffenedShellTest.cpp
@Date    :   2023/04/14
@Author  :   Alasdair Christison Gray
@Description :
*/

// =============================================================================
// Standard Library Includes
// =============================================================================

// =============================================================================
// Extension Includes
// =============================================================================
#define REAL TacsRealPart
#define _USE_MATH_DEFINES
#include <math.h>

#include "TACSBladeStiffenedShellConstitutive.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSMaterialProperties.h"

// =============================================================================
// Global constant definitions
// =============================================================================
const int NUM_STRESS = TACSShellConstitutive::NUM_STRESSES;
const int NUM_STIFF = 22;
const int NUM_Q = 6;
const int NUM_ABAR = 3;
const int NUM_C = TACSShellConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES;
const int NUM_PANEL_PLY = 3;
const int NUM_STIFF_PLY = 2;
const int NUM_DV = 5 + NUM_PANEL_PLY + NUM_STIFF_PLY;
const double DEG2RAD = M_PI / 180.0;
const TacsScalar FD_STEP = 1e-8;
const TacsScalar KDRILL = 10.;
const TacsScalar strain[9] = {1e-3, 1e-3, 1e-3, 1e-3, 1e-3,
                              1e-3, 1e-3, 1e-3, 1e-3};

// =============================================================================
// Function prototypes
// =============================================================================

void printStiffnessMatrix(TacsScalar *C);
void computeStress(TacsScalar C[], const TacsScalar e[], TacsScalar stress[]);

// =============================================================================
// Main
// =============================================================================

int main() {
  // Define isotropic material properties
  TacsScalar rho = 1550.0;
  TacsScalar specific_heat = 921.096;
  TacsScalar E1 = 54e3;
  TacsScalar E2 = 18e3;
  TacsScalar E3 = 18e3;
  TacsScalar nu12 = 0.25;
  TacsScalar nu13 = 0.25;
  TacsScalar nu23 = 0.25;
  TacsScalar G12 = 9e3;
  TacsScalar G13 = 9e3;
  TacsScalar G23 = 5e3;
  TacsScalar T1 = 2410.0;
  TacsScalar C1 = 1040.0;
  TacsScalar T2 = 73.0;
  TacsScalar C2 = 173.0;
  TacsScalar S12 = 71.0;
  TacsScalar cte = 24.0e-6;
  TacsScalar kappa = 230.0;
  TacsScalar kcorr = 5. / 6.;

  // Create material properties and orthotropic ply objects
  TACSMaterialProperties compositeProps =
      TACSMaterialProperties(rho, specific_heat, E1, E2, E3, nu12, nu13, nu23,
                             G12, G13, G23, T1, C1, T2, C2, T2, C2, S12);
  compositeProps.incref();

  TACSOrthotropicPly compositePly = TACSOrthotropicPly(1.0, &compositeProps);
  compositePly.incref();

  compositePly.testFailSens(1e-6, 21.374658736);

  // Define bladeStiffenedShell parameters
  const TacsScalar panelLength = 1.0;
  const int panelLengthNum = 0;
  const TacsScalar stiffenerPitch = 0.25;
  const int stiffenerPitchNum = 1;

  const TacsScalar panelThickness = 1.5e-2;
  const int panelThicknessNum = 2;

  TacsScalar panelPlyAngles[NUM_PANEL_PLY] = {0., 45., 90.};
  for (int ii = 0; ii < NUM_PANEL_PLY; ii++) {
    panelPlyAngles[ii] *= DEG2RAD;
  }
  TacsScalar panelPlyFracs[NUM_PANEL_PLY] = {0.5, 0.2, 0.3};
  int panelPlyFracNums[NUM_PANEL_PLY] = {3, 4, 5};

  const TacsScalar stiffenerHeight = 0.075;
  const int stiffenerHeightNum = 6;
  const TacsScalar stiffenerThickness = 1e-2;
  const int stiffenerThicknessNum = 7;

  TacsScalar stiffenerPlyAngles[NUM_STIFF_PLY] = {0., 60.};
  for (int ii = 0; ii < NUM_STIFF_PLY; ii++) {
    stiffenerPlyAngles[ii] *= DEG2RAD;
  }
  TacsScalar stiffenerPlyFracs[NUM_STIFF_PLY] = {0.7, 0.3};
  int stiffenerPlyFracNums[NUM_STIFF_PLY] = {8, 9};

  TacsScalar DV0[NUM_DV];
  DV0[panelLengthNum] = panelLength;
  DV0[stiffenerPitchNum] = stiffenerPitch;
  DV0[stiffenerHeightNum] = stiffenerHeight;
  DV0[stiffenerThicknessNum] = stiffenerThickness;
  DV0[panelThicknessNum] = panelThickness;

  for (int ii = 0; ii < NUM_PANEL_PLY; ii++) {
    DV0[panelPlyFracNums[ii]] = panelPlyFracs[ii];
  }
  for (int ii = 0; ii < NUM_STIFF_PLY; ii++) {
    DV0[stiffenerPlyFracNums[ii]] = stiffenerPlyFracs[ii];
  }

  const TacsScalar flangeFraction = 0.8;

  // Compute the Q and Abar matrics for the panel and stiffener laminates
  TacsScalar QPanel[NUM_Q], QStiff[NUM_Q], AbarPanel[NUM_ABAR],
      AbarStiff[NUM_ABAR];
  memset(QPanel, 0, NUM_Q * sizeof(TacsScalar));
  memset(QStiff, 0, NUM_Q * sizeof(TacsScalar));
  memset(AbarPanel, 0, NUM_ABAR * sizeof(TacsScalar));
  memset(AbarStiff, 0, NUM_ABAR * sizeof(TacsScalar));

  for (int ii = 0; ii < NUM_PANEL_PLY; ii++) {
    TacsScalar Qtemp[NUM_Q], AbarTemp[NUM_ABAR];
    compositePly.calculateQbar(panelPlyAngles[ii], Qtemp);
    compositePly.calculateAbar(panelPlyAngles[ii], AbarTemp);

    for (int jj = 0; jj < NUM_Q; jj++) {
      QPanel[jj] += panelPlyFracs[ii] * Qtemp[jj];
    }
    for (int jj = 0; jj < NUM_ABAR; jj++) {
      AbarPanel[jj] += panelPlyFracs[ii] * AbarTemp[jj];
    }
  }
  printf("\nQPanel = \n");
  printf("[% 011.7e, % 011.7e, % 011.7e]\n", REAL(QPanel[0]), REAL(QPanel[1]),
         REAL(QPanel[2]));
  printf("[% 011.7e, % 011.7e, % 011.7e]\n", REAL(QPanel[1]), REAL(QPanel[3]),
         REAL(QPanel[4]));
  printf("[% 011.7e, % 011.7e, % 011.7e]\n", REAL(QPanel[2]), REAL(QPanel[4]),
         REAL(QPanel[5]));
  printf("\nAbarPanel = \n");
  printf("[% 011.7e, % 011.7e]\n", REAL(AbarPanel[0]), REAL(AbarPanel[1]));
  printf("[% 011.7e, % 011.7e]\n", REAL(AbarPanel[1]), REAL(AbarPanel[2]));

  TacsScalar EStiff = 0;
  TacsScalar GStiff = 0;
  for (int ii = 0; ii < NUM_STIFF_PLY; ii++) {
    TacsScalar QTemp[NUM_Q], AbarTemp[NUM_ABAR];
    compositePly.calculateQbar(stiffenerPlyAngles[ii], QTemp);
    compositePly.calculateAbar(stiffenerPlyAngles[ii], AbarTemp);

    printf("\nAngle = %f, QTemp = \n", stiffenerPlyAngles[ii]);
    printf("[% 011.7e, % 011.7e, % 011.7e]\n", REAL(QTemp[0]), REAL(QTemp[1]),
           REAL(QTemp[2]));
    printf("[% 011.7e, % 011.7e, % 011.7e]\n", REAL(QTemp[1]), REAL(QTemp[3]),
           REAL(QTemp[4]));
    printf("[% 011.7e, % 011.7e, % 011.7e]\n", REAL(QTemp[2]), REAL(QTemp[4]),
           REAL(QTemp[5]));
    printf("\nAbarTemp = \n");
    printf("[% 011.7e, % 011.7e]\n", REAL(AbarTemp[0]), REAL(AbarTemp[1]));
    printf("[% 011.7e, % 011.7e]\n", REAL(AbarTemp[1]), REAL(AbarTemp[2]));

    for (int jj = 0; jj < NUM_Q; jj++) {
      QStiff[jj] += stiffenerPlyFracs[ii] * QTemp[jj];
    }
    for (int jj = 0; jj < NUM_ABAR; jj++) {
      AbarStiff[jj] += stiffenerPlyFracs[ii] * AbarTemp[jj];
    }
    EStiff +=
        stiffenerPlyFracs[ii] * (QTemp[0] - (QTemp[1] * QTemp[1]) / QTemp[3]);
    GStiff += stiffenerPlyFracs[ii] * QTemp[5];
  }
  printf("\nQStiff = \n");
  printf("[% 011.7e, % 011.7e, % 011.7e]\n", REAL(QStiff[0]), REAL(QStiff[1]),
         REAL(QStiff[2]));
  printf("[% 011.7e, % 011.7e, % 011.7e]\n", REAL(QStiff[1]), REAL(QStiff[3]),
         REAL(QStiff[4]));
  printf("[% 011.7e, % 011.7e, % 011.7e]\n", REAL(QStiff[2]), REAL(QStiff[4]),
         REAL(QStiff[5]));
  printf("\nAbarStiff = \n");
  printf("[%011.7e, % 011.7e]\n", REAL(AbarStiff[0]), REAL(AbarStiff[1]));
  printf("[% 011.7e, %011.7e]\n", REAL(AbarStiff[1]), REAL(AbarStiff[2]));

  printf("EStiff = % 011.7e\nGStiff = % 011.7e\n", REAL(EStiff), REAL(GStiff));

  // ==============================================================================
  // Create the stiffened shell objects
  // ==============================================================================

  // Create bladeStiffenedShell object
  TACSBladeStiffenedShellConstitutive con = TACSBladeStiffenedShellConstitutive(
      &compositePly, &compositePly, kcorr, panelLength, panelLengthNum,
      stiffenerPitch, stiffenerPitchNum, panelThickness, panelThicknessNum,
      NUM_PANEL_PLY, panelPlyAngles, panelPlyFracs, panelPlyFracNums,
      stiffenerHeight, stiffenerHeightNum, stiffenerThickness,
      stiffenerThicknessNum, NUM_STIFF_PLY, stiffenerPlyAngles,
      stiffenerPlyFracs, stiffenerPlyFracNums, flangeFraction);
  con.incref();
  con.testGlobalBucklingStiffnessSens();

  // Create a version with no stiffener so we can check against the isoShell
  TACSBladeStiffenedShellConstitutive conPanelOnly =
      TACSBladeStiffenedShellConstitutive(
          &compositePly, &compositePly, kcorr, panelLength, panelLengthNum,
          stiffenerPitch, stiffenerPitchNum, panelThickness, panelThicknessNum,
          NUM_PANEL_PLY, panelPlyAngles, panelPlyFracs, panelPlyFracNums,
          stiffenerHeight, stiffenerHeightNum, 0., stiffenerThicknessNum,
          NUM_STIFF_PLY, stiffenerPlyAngles, stiffenerPlyFracs,
          stiffenerPlyFracNums, flangeFraction);
  conPanelOnly.incref();

  // Create a version with no panel
  TACSBladeStiffenedShellConstitutive conStiffenerOnly =
      TACSBladeStiffenedShellConstitutive(
          &compositePly, &compositePly, kcorr, panelLength, panelLengthNum,
          stiffenerPitch, stiffenerPitchNum, 0., panelThicknessNum,
          NUM_PANEL_PLY, panelPlyAngles, panelPlyFracs, panelPlyFracNums,
          stiffenerHeight, stiffenerHeightNum, stiffenerThickness,
          stiffenerThicknessNum, NUM_STIFF_PLY, stiffenerPlyAngles,
          stiffenerPlyFracs, stiffenerPlyFracNums, flangeFraction);
  conStiffenerOnly.incref();

  //==================================================
  // Validate the Q and Abar matrices for the panel and stiffener laminates
  //============================
  TacsScalar Q[6], Abar[3];
  con.computeSmearedStiffness(NUM_PANEL_PLY, con.panelQMats, con.panelAbarMats,
                              con.panelPlyFracs, Q, Abar);
  printf("\nError in Q for panel:\n");
  for (int ii = 0; ii < NUM_Q; ii++) {
    printf("Q[%d] rel error = % 011.7e\n", ii,
           REAL((Q[ii] - QPanel[ii]) / QPanel[ii]));
  }
  printf("\nError in Abar for panel:\n");
  for (int ii = 0; ii < NUM_ABAR; ii++) {
    printf("Abar[%d] rel error = % 011.7e\n", ii,
           REAL((Abar[ii] - AbarPanel[ii]) / AbarPanel[ii]));
  }

  con.computeSmearedStiffness(NUM_STIFF_PLY, con.stiffenerQMats,
                              con.stiffenerAbarMats, con.stiffenerPlyFracs, Q,
                              Abar);
  printf("\nError in Q for stiffener:\n");
  for (int ii = 0; ii < NUM_Q; ii++) {
    printf("Q[%d] rel error = % 011.7e\n", ii,
           REAL((Q[ii] - QStiff[ii]) / QStiff[ii]));
  }
  printf("\nError in Abar for stiffener:\n");
  for (int ii = 0; ii < NUM_ABAR; ii++) {
    printf("Abar[%d] rel error = % 011.7e\n", ii,
           REAL((Abar[ii] - AbarStiff[ii]) / AbarStiff[ii]));
  }

  // =============================================================================
  // Validate the effective modulii of the stiffener laminate
  // ==============================================================================
  TacsScalar Econ, Gcon;
  con.computeEffectiveModulii(NUM_STIFF_PLY, con.stiffenerQMats,
                              con.stiffenerPlyFracs, &Econ, &Gcon);
  printf("\nError in computeEffectiveModulii:\n");
  printf("Econ rel error = % 011.7e\nGcon rel error = % 011.7e\n",
         REAL((Econ - EStiff) / EStiff), REAL((Gcon - GStiff) / GStiff));

  //
  // ==============================================================================
  // Validate the stiffener section properties
  // ==============================================================================
  TacsScalar Aref, A1, A2, Iref, I1, I2, Jref, Zref, z1, z2;
  A1 = stiffenerHeight * stiffenerThickness;
  A2 = stiffenerHeight * flangeFraction * stiffenerThickness;
  Aref = A1 + A2;

  printf("\n                      Aref = % 011.7e\n", REAL(Aref));
  printf("con.computeStiffenerArea() = % 011.7e\n",
         REAL(con.computeStiffenerArea()));
  printf("Rel error = % 011.7e\n",
         REAL((con.computeStiffenerArea() - Aref) / Aref));

  z2 = -stiffenerThickness / 2.0;
  z1 = -stiffenerThickness + -stiffenerHeight / 2.0;
  Zref = (A1 * z1 + A2 * z2) / Aref;
  printf("\n                                Zref = % 011.7e\n", REAL(Zref));
  printf("con.computeStiffenerCentroidHeight() = % 011.7e\n",
         REAL(con.computeStiffenerCentroidHeight()));
  printf("Rel error = % 011.7e\n",
         REAL((con.computeStiffenerCentroidHeight() - Zref) / Zref));

  TacsScalar st3 = stiffenerThickness * stiffenerThickness * stiffenerThickness;
  TacsScalar sh3 = stiffenerHeight * stiffenerHeight * stiffenerHeight;

  I1 = stiffenerThickness * sh3 / 12.;
  I2 = flangeFraction * stiffenerHeight * st3 / 12.;
  Iref =
      I1 + I2 + A1 * (z1 - Zref) * (z1 - Zref) + A2 * (z2 - Zref) * (z2 - Zref);
  printf("\n                     Iref = % 011.7e\n", REAL(Iref));
  printf("con.computeStiffenerIzz() = % 011.7e\n",
         REAL(con.computeStiffenerIzz()));
  printf("Rel error = % 011.7e\n",
         REAL((con.computeStiffenerIzz() - Iref) / Iref));

  Jref = (flangeFraction * stiffenerHeight * st3 + stiffenerHeight * st3) / 3.;
  printf("\n                     Jref = % 011.7e\n", REAL(Jref));
  printf("con.computeStiffenerJxx() = % 011.7e\n",
         REAL(con.computeStiffenerJxx()));
  printf("Rel error = % 011.7e\n",
         REAL((con.computeStiffenerJxx() - Jref) / Jref));

  //========================================
  // Validate the sensitivities of the stiffener section properties against FD
  //========================================

  TacsScalar DVPert[NUM_DV];
  for (int ii = 0; ii < NUM_DV; ii++) {
    DVPert[ii] = DV0[ii];
  }

  // --- Area ---
  TacsScalar dAdtRef, dAdhRef, dAdt, dAdh;
  con.computeStiffenerAreaSens(dAdt, dAdh);

  DVPert[stiffenerHeightNum] += FD_STEP;
  con.setDesignVars(0, NUM_DV, DVPert);
  dAdhRef = (con.computeStiffenerArea() - Aref) / FD_STEP;
  DVPert[stiffenerHeightNum] -= FD_STEP;

  DVPert[stiffenerThicknessNum] += FD_STEP;
  con.setDesignVars(0, NUM_DV, DVPert);
  dAdtRef = (con.computeStiffenerArea() - Aref) / FD_STEP;
  DVPert[stiffenerThicknessNum] -= FD_STEP;

  printf("\n\ndAdtRef = % 011.7e\ndAdhRef = % 011.7e\n", REAL(dAdtRef),
         REAL(dAdhRef));
  printf("dAdt = % 011.7e\ndAdh = % 011.7e\n", REAL(dAdt), REAL(dAdh));
  printf("\nError in computeStiffenerAreaSens:\n");
  printf("dAdt rel error = % 011.7e\ndAdh rel error = % 011.7e\n",
         REAL((dAdt - dAdtRef) / dAdtRef), REAL((dAdh - dAdhRef) / dAdhRef));

  // --- Centroid height ---
  TacsScalar dZdtRef, dZdhRef, dZdt, dZdh;
  con.computeStiffenerCentroidHeightSens(dZdt, dZdh);

  DVPert[stiffenerHeightNum] += FD_STEP;
  con.setDesignVars(0, NUM_DV, DVPert);
  dZdhRef = (con.computeStiffenerCentroidHeight() - Zref) / FD_STEP;
  DVPert[stiffenerHeightNum] -= FD_STEP;

  DVPert[stiffenerThicknessNum] += FD_STEP;
  con.setDesignVars(0, NUM_DV, DVPert);
  dZdtRef = (con.computeStiffenerCentroidHeight() - Zref) / FD_STEP;
  DVPert[stiffenerThicknessNum] -= FD_STEP;

  printf("\n\ndZdtRef = % 011.7e\ndZdhRef = % 011.7e\n", REAL(dZdtRef),
         REAL(dZdhRef));
  printf("dZdt = % 011.7e\ndZdh = % 011.7e\n", REAL(dZdt), REAL(dZdh));
  printf("\nError in computeStiffenerCentroidHeightSens:\n");
  printf("dZdt rel error = % 011.7e\ndZdh rel error = % 011.7e\n",
         REAL((dZdt - dZdtRef) / dZdtRef), REAL((dZdh - dZdhRef) / dZdhRef));

  // --- Izz ---
  TacsScalar dIdtRef, dIdhRef, dIdt, dIdh;
  con.computeStiffenerIzzSens(dIdt, dIdh);

  DVPert[stiffenerHeightNum] += FD_STEP;
  con.setDesignVars(0, NUM_DV, DVPert);
  dIdhRef = (con.computeStiffenerIzz() - Iref) / FD_STEP;
  DVPert[stiffenerHeightNum] -= FD_STEP;

  DVPert[stiffenerThicknessNum] += FD_STEP;
  con.setDesignVars(0, NUM_DV, DVPert);
  dIdtRef = (con.computeStiffenerIzz() - Iref) / FD_STEP;
  DVPert[stiffenerThicknessNum] -= FD_STEP;

  printf("\n\ndIdtRef = % 011.7e\ndIdhRef = % 011.7e\n", REAL(dIdtRef),
         REAL(dIdhRef));
  printf("dIdt = % 011.7e\ndIdh = % 011.7e\n", REAL(dIdt), REAL(dIdh));
  printf("\nError in computeStiffenerIzzSens:\n");
  printf("dIdt rel error = % 011.7e\ndIdh rel error = % 011.7e\n",
         REAL((dIdt - dIdtRef) / dIdtRef), REAL((dIdh - dIdhRef) / dIdhRef));

  // --- J ---
  TacsScalar dJdtRef, dJdhRef, dJdt, dJdh;
  con.computeStiffenerJxxSens(dJdt, dJdh);

  DVPert[stiffenerHeightNum] += FD_STEP;
  con.setDesignVars(0, NUM_DV, DVPert);
  dJdhRef = (con.computeStiffenerJxx() - Jref) / FD_STEP;
  DVPert[stiffenerHeightNum] -= FD_STEP;

  DVPert[stiffenerThicknessNum] += FD_STEP;
  con.setDesignVars(0, NUM_DV, DVPert);
  dJdtRef = (con.computeStiffenerJxx() - Jref) / FD_STEP;
  DVPert[stiffenerThicknessNum] -= FD_STEP;

  printf("\n\ndJdtRef = % 011.7e\ndJdhRef = % 011.7e\n", REAL(dJdtRef),
         REAL(dJdhRef));
  printf("dJdt = % 011.7e\ndJdh = % 011.7e\n", REAL(dJdt), REAL(dJdh));
  printf("\nError in computeStiffenerJxxSens:\n");
  printf("dJdt rel error = % 011.7e\ndJdh rel error = % 011.7e\n",
         REAL((dJdt - dJdtRef) / dJdtRef), REAL((dJdh - dJdhRef) / dJdhRef));

  // --- MOI ---
  TacsScalar MOIRef, dMOIdtRef, dMOIdhRef, dMOIdt, dMOIdh;
  con.setDesignVars(0, NUM_DV, DVPert);
  MOIRef = con.computeStiffenerMOI();
  con.computeStiffenerMOISens(dMOIdt, dMOIdh);

  DVPert[stiffenerHeightNum] += FD_STEP;
  con.setDesignVars(0, NUM_DV, DVPert);
  dMOIdhRef = (con.computeStiffenerMOI() - MOIRef) / FD_STEP;
  DVPert[stiffenerHeightNum] -= FD_STEP;

  DVPert[stiffenerThicknessNum] += FD_STEP;
  con.setDesignVars(0, NUM_DV, DVPert);
  dMOIdtRef = (con.computeStiffenerMOI() - MOIRef) / FD_STEP;
  DVPert[stiffenerThicknessNum] -= FD_STEP;
  con.setDesignVars(0, NUM_DV, DVPert);

  printf("\n\ndMOIdtRef = % 011.7e\ndMOIdhRef = % 011.7e\n", REAL(dMOIdtRef),
         REAL(dMOIdhRef));
  printf("dMOIdt = % 011.7e\ndMOIdh = % 011.7e\n", REAL(dMOIdt), REAL(dMOIdh));
  printf("\nError in computeStiffenerMOISens:\n");
  printf("dMOIdt rel error = % 011.7e\ndMOIdh rel error = % 011.7e\n",
         REAL((dMOIdt - dMOIdtRef) / dMOIdtRef),
         REAL((dMOIdh - dMOIdhRef) / dMOIdhRef));

  //==============================================================================
  // Validate the stiffness matrix and stress calculations
  //==============================================================================

  // --- Panel only ---
  TacsScalar CPanelOnly[NUM_C];
  TacsScalar CPanelOnlyRef[NUM_C];
  TacsScalar CPanelOnlyDiff[NUM_C];
  memset(CPanelOnly, 0, sizeof(TacsScalar) * NUM_C);
  memset(CPanelOnlyRef, 0, sizeof(TacsScalar) * NUM_C);
  memset(CPanelOnlyDiff, 0, sizeof(TacsScalar) * NUM_C);

  for (int ii = 0; ii < NUM_Q; ii++) {
    CPanelOnlyRef[ii] = QPanel[ii] * panelThickness;  // A Matrix entries
    CPanelOnlyRef[ii + (NUM_Q * 2)] =
        QPanel[ii] * (panelThickness * panelThickness * panelThickness /
                      12.0);  // D Matrix entries
  }
  for (int ii = 0; ii < NUM_ABAR; ii++) {
    CPanelOnlyRef[ii + (NUM_Q * 3)] =
        AbarPanel[ii] * panelThickness * kcorr;  // As Matrix entries
  }
  CPanelOnlyRef[21] =
      KDRILL * 0.5 * (CPanelOnlyRef[18] + CPanelOnlyRef[20]);  // Drill entry

  conPanelOnly.evalTangentStiffness(0, NULL, NULL, CPanelOnly);

  bool noDiff = true;
  for (int ii = 0; ii < NUM_C; ii++) {
    CPanelOnlyDiff[ii] = CPanelOnly[ii] - CPanelOnlyRef[ii];
    if (fabs((REAL(CPanelOnlyDiff[ii]))) > 1e-10) {
      noDiff = false;
    }
  }
  printf("\n\nPanel Only Stiffness Matrix Ref:\n");
  printStiffnessMatrix(CPanelOnlyRef);
  printf("\nPanel Only Stiffness Matrix:\n");
  printStiffnessMatrix(CPanelOnly);
  printf("\nPanel Only Stiffness Diff:\n");
  printStiffnessMatrix(CPanelOnlyDiff);

  //
  printf("\n\n============================================================\n");
  if (noDiff) {
    printf("Panel Only Stiffness Matrix is correct!\n");
  } else {
    printf("Panel Only Stiffness Matrix is incorrect!\n");
  }
  printf("============================================================\n");

  TacsScalar stressPanelOnly[NUM_STRESS];
  TacsScalar stressPanelOnlyRef[NUM_STRESS];
  TacsScalar stressPanelOnlyDiff[NUM_STRESS];
  memset(stressPanelOnly, 0, sizeof(TacsScalar) * NUM_STRESS);
  memset(stressPanelOnlyRef, 0, sizeof(TacsScalar) * NUM_STRESS);
  memset(stressPanelOnlyDiff, 0, sizeof(TacsScalar) * NUM_STRESS);
  computeStress(CPanelOnlyRef, strain, stressPanelOnlyRef);
  conPanelOnly.evalStress(0, NULL, NULL, strain, stressPanelOnly);

  noDiff = true;
  printf("\n\nPanel Only Stress:\n");
  for (int ii = 0; ii < NUM_STRESS; ii++) {
    stressPanelOnlyDiff[ii] = stressPanelOnly[ii] - stressPanelOnlyRef[ii];
    if (fabs(REAL(stressPanelOnlyDiff[ii])) > 1e-10) {
      noDiff = false;
    }
    printf("stressRef[%d] = % 011.7e, stress[%d] = % 011.7e, diff = % 011.7e\n",
           ii, REAL(stressPanelOnlyRef[ii]), ii, REAL(stressPanelOnly[ii]),
           REAL(stressPanelOnlyDiff[ii]));
  }

  //
  printf("\n\n============================================================\n");
  if (noDiff) {
    printf("Panel Only Stress is correct!\n");
  } else {
    printf("Panel Only Stress is incorrect!\n");
  }
  printf("============================================================\n");

  TacsScalar psi[NUM_STRESS];
  TacsScalar scale = 0.1235;
  for (int ii = 0; ii < NUM_STRESS; ii++) {
    psi[ii] = 1.0;
  }
  TacsScalar stressDVSensPanelOnly[NUM_DV];
  TacsScalar stressDVSensPanelOnlyRef[NUM_DV];
  TacsScalar stressDVSensPanelOnlyDiff[NUM_DV];
  memset(stressDVSensPanelOnly, 0, sizeof(TacsScalar) * NUM_DV);
  memset(stressDVSensPanelOnlyRef, 0, sizeof(TacsScalar) * NUM_DV);
  memset(stressDVSensPanelOnlyDiff, 0, sizeof(TacsScalar) * NUM_DV);

  conPanelOnly.addStressDVSens(0, scale, NULL, NULL, strain, psi, NUM_DV,
                               stressDVSensPanelOnly);

  for (int ii = 0; ii < NUM_DV; ii++) {
    DVPert[ii] += FD_STEP;
    DVPert[stiffenerThicknessNum] = 0.;
    conPanelOnly.setDesignVars(0, NUM_DV, DVPert);

    TacsScalar stressPert[NUM_STRESS];
    memset(stressPert, 0, sizeof(TacsScalar) * NUM_STRESS);
    conPanelOnly.evalStress(0, NULL, NULL, strain, stressPert);

    for (int jj = 0; jj < NUM_STRESS; jj++) {
      TacsScalar dsdx = (stressPert[jj] - stressPanelOnlyRef[jj]) / FD_STEP;
      stressDVSensPanelOnlyRef[ii] += scale * psi[jj] * dsdx;
    }

    DVPert[ii] -= FD_STEP;
    DVPert[stiffenerThicknessNum] = 0.;
    conPanelOnly.setDesignVars(0, NUM_DV, DVPert);
  }

  noDiff = true;
  // printf("\n\nPanel Only Stress DV Sens:\n");
  for (int ii = 0; ii < NUM_DV; ii++) {
    stressDVSensPanelOnlyDiff[ii] =
        stressDVSensPanelOnly[ii] - stressDVSensPanelOnlyRef[ii];
    if (fabs(REAL(stressDVSensPanelOnlyDiff[ii])) > 1e-5) {
      noDiff = false;
    }
    printf(
        "stressDVSensRef[%d] = % 011.7e, stressDVSens[%d] = % 011.7e, diff = % "
        "011.7e\n",
        ii, REAL(stressDVSensPanelOnlyRef[ii]), ii,
        REAL(stressDVSensPanelOnly[ii]), REAL(stressDVSensPanelOnlyDiff[ii]));
  }

  //
  printf("\n\n============================================================\n");
  if (noDiff) {
    printf("Panel Only Stress DV sensitivity is correct!\n");
  } else {
    printf("Panel Only Stress DV sensitivity is incorrect!\n");
  }
  printf("============================================================\n");

  exit(0);

  // --- Stiffener only ---
  TacsScalar CStiffenerOnly[NUM_C];
  TacsScalar CStiffenerOnlyRef[NUM_C];
  TacsScalar CStiffenerOnlyDiff[NUM_C];
  memset(CStiffenerOnly, 0, sizeof(TacsScalar) * NUM_C);
  memset(CStiffenerOnlyRef, 0, sizeof(TacsScalar) * NUM_C);
  memset(CStiffenerOnlyDiff, 0, sizeof(TacsScalar) * NUM_C);

  TacsScalar pInv = 1.0 / stiffenerPitch;

  CStiffenerOnlyRef[0] = pInv * Aref * EStiff;
  CStiffenerOnlyRef[5] = pInv * Aref * GStiff * kcorr / 4.0;
  CStiffenerOnlyRef[6] = pInv * Aref * EStiff * Zref;
  CStiffenerOnlyRef[11] = pInv * kcorr * GStiff * Aref * Zref / 4.0;
  CStiffenerOnlyRef[12] = pInv * EStiff * (Iref + Aref * Zref * Zref);
  CStiffenerOnlyRef[17] =
      pInv * GStiff * (Jref + kcorr * Aref * Zref * Zref) / 4.0;
  CStiffenerOnlyRef[20] = pInv * kcorr * GStiff * Aref;
  CStiffenerOnlyRef[21] = 0.5 * KDRILL * CStiffenerOnlyRef[20];

  conStiffenerOnly.evalTangentStiffness(0, NULL, NULL, CStiffenerOnly);

  noDiff = true;
  for (int ii = 0; ii < NUM_C; ii++) {
    CStiffenerOnlyDiff[ii] = CStiffenerOnly[ii] - CStiffenerOnlyRef[ii];
    if (fabs(REAL(CStiffenerOnlyDiff[ii])) > 1e-10) {
      noDiff = false;
    }
  }
  printf("\n\nStiffener Only Stiffness Matrix Ref:\n");
  printStiffnessMatrix(CStiffenerOnlyRef);
  printf("\nStiffener Only Stiffness Matrix:\n");
  printStiffnessMatrix(CStiffenerOnly);
  printf("\nStiffener Only Stiffness Diff:\n");
  printStiffnessMatrix(CStiffenerOnlyDiff);

  //
  printf("\n\n============================================================\n");
  if (noDiff) {
    printf("Stiffener Only Stiffness Matrix is correct!\n");
  } else {
    printf("Stiffener Only Stiffness Matrix is incorrect!\n");
  }
  printf("============================================================\n");

  TacsScalar stressStiffenerOnly[NUM_STRESS];
  TacsScalar stressStiffenerOnlyRef[NUM_STRESS];
  TacsScalar stressStiffenerOnlyDiff[NUM_STRESS];
  memset(stressStiffenerOnly, 0, sizeof(TacsScalar) * NUM_STRESS);
  memset(stressStiffenerOnlyRef, 0, sizeof(TacsScalar) * NUM_STRESS);
  memset(stressStiffenerOnlyDiff, 0, sizeof(TacsScalar) * NUM_STRESS);
  computeStress(CStiffenerOnlyRef, strain, stressStiffenerOnlyRef);
  conStiffenerOnly.evalStress(0, NULL, NULL, strain, stressStiffenerOnly);

  noDiff = true;
  printf("\n\nStiffener Only Stress:\n");
  for (int ii = 0; ii < NUM_STRESS; ii++) {
    stressStiffenerOnlyDiff[ii] =
        stressStiffenerOnly[ii] - stressStiffenerOnlyRef[ii];
    if (fabs(REAL(stressStiffenerOnlyDiff[ii])) > 1e-10) {
      noDiff = false;
    }
    printf("stressRef[%d] = % 011.7e, stress[%d] = % 011.7e, diff = %011.7e\n",
           ii, REAL(stressStiffenerOnlyRef[ii]), ii,
           REAL(stressStiffenerOnly[ii]), REAL(stressStiffenerOnlyDiff[ii]));
  }

  //
  printf("\n\n============================================================\n");
  if (noDiff) {
    printf("Stiffener Only Stress is correct!\n");
  } else {
    printf("Stiffener Only Stress is incorrect!\n");
  }
  printf("============================================================\n");

  TacsScalar stressDVSensStiffenerOnly[NUM_DV];
  TacsScalar stressDVSensStiffenerOnlyRef[NUM_DV];
  TacsScalar stressDVSensStiffenerOnlyDiff[NUM_DV];
  memset(stressDVSensStiffenerOnly, 0, sizeof(TacsScalar) * NUM_DV);
  memset(stressDVSensStiffenerOnlyRef, 0, sizeof(TacsScalar) * NUM_DV);
  memset(stressDVSensStiffenerOnlyDiff, 0, sizeof(TacsScalar) * NUM_DV);

  conStiffenerOnly.addStressDVSens(0, scale, NULL, NULL, strain, psi, NUM_DV,
                                   stressDVSensStiffenerOnly);

  for (int ii = 0; ii < NUM_DV; ii++) {
    DVPert[ii] = DV0[ii];
  }

  for (int ii = 0; ii < NUM_DV; ii++) {
    DVPert[ii] += FD_STEP;
    if (ii == panelThicknessNum) {
      DVPert[panelThicknessNum] = FD_STEP;
    } else {
      DVPert[panelThicknessNum] = 0.;
    }
    conStiffenerOnly.setDesignVars(0, NUM_DV, DVPert);

    TacsScalar stressPert[NUM_STRESS];
    memset(stressPert, 0, sizeof(TacsScalar) * NUM_STRESS);
    conStiffenerOnly.evalStress(0, NULL, NULL, strain, stressPert);

    for (int jj = 0; jj < NUM_STRESS; jj++) {
      TacsScalar dsdx = (stressPert[jj] - stressStiffenerOnlyRef[jj]) / FD_STEP;
      stressDVSensStiffenerOnlyRef[ii] += scale * psi[jj] * dsdx;
    }

    DVPert[ii] -= FD_STEP;
    DVPert[panelThicknessNum] = 0.;
    conStiffenerOnly.setDesignVars(0, NUM_DV, DVPert);
  }

  noDiff = true;
  printf("\n\nPanel Only Stress DV Sens:\n");
  for (int ii = 0; ii < NUM_DV; ii++) {
    stressDVSensStiffenerOnlyDiff[ii] =
        stressDVSensStiffenerOnly[ii] - stressDVSensStiffenerOnlyRef[ii];
    if (fabs(REAL(stressDVSensStiffenerOnlyDiff[ii])) > 1e-5) {
      noDiff = false;
    }
    printf(
        "stressDVSensRef[%d] = % 011.7e, stressDVSens[%d] = % 011.7e, diff = % "
        "011.7e\n",
        ii, REAL(stressDVSensStiffenerOnlyRef[ii]), ii,
        REAL(stressDVSensStiffenerOnly[ii]),
        REAL(stressDVSensStiffenerOnlyDiff[ii]));
  }

  printf("\n\n============================================================\n");
  if (noDiff) {
    printf("Stiffener Only Stress DV sensitivity is correct!\n");
  } else {
    printf("Stiffener Only Stress DV sensitivity is incorrect!\n");
  }
  printf("============================================================\n");

  //==============================================================================
  // Test the different parts of the stress DV sensitivity
  //==============================================================================
  // Test stiffStress * d/dx(Te) * psi
  TacsScalar stiffStress[6];
  memset(stiffStress, 0, sizeof(TacsScalar) * 6);
  TacsScalar dfdxTest[NUM_DV];
  memset(dfdxTest, 0, sizeof(TacsScalar) * NUM_DV);
  TacsScalar productRef = 0.;
  TacsScalar psiTrans[6];
  memset(psiTrans, 0, sizeof(TacsScalar) * 6);
  conStiffenerOnly.transformStrain(psi, psiTrans);
  conStiffenerOnly.computeStiffenerStress(strain, stiffStress);
  for (int ii = 0; ii < 6; ii++) {
    productRef += stiffStress[ii] * psiTrans[ii];
  }

  conStiffenerOnly.addStrainTransformProductDVsens(stiffStress, psi, scale,
                                                   dfdxTest);
  TacsScalar dfdxTestRef[NUM_DV];
  for (int ii = 0; ii < NUM_DV; ii++) {
    DVPert[ii] += FD_STEP;
    if (ii == panelThicknessNum) {
      DVPert[panelThicknessNum] = FD_STEP;
    } else {
      DVPert[panelThicknessNum] = 0.;
    }
    conStiffenerOnly.setDesignVars(0, NUM_DV, DVPert);

    dfdxTestRef[ii] = 0.;
    memset(psiTrans, 0, sizeof(TacsScalar) * 6);
    conStiffenerOnly.transformStrain(psi, psiTrans);

    TacsScalar product = 0.;
    for (int ii = 0; ii < 6; ii++) {
      product += stiffStress[ii] * psiTrans[ii];
    }
    dfdxTestRef[ii] = scale * (product - productRef) / FD_STEP;

    printf(
        "DV Number %d:    psiTrans = [% 015.11e, % 015.11e, % 015.11e, % "
        "015.11e, % 015.11e, % 015.11e]\n",
        ii, REAL(psiTrans[0]), REAL(psiTrans[1]), REAL(psiTrans[2]),
        REAL(psiTrans[3]), REAL(psiTrans[4]), REAL(psiTrans[5]));

    DVPert[ii] -= FD_STEP;
    DVPert[panelThicknessNum] = 0.;
    conStiffenerOnly.setDesignVars(0, NUM_DV, DVPert);
  }

  printf("\n\nStiffener Only Strain Transform DV Sensitivity Test:\n");
  noDiff = true;
  for (int ii = 0; ii < NUM_DV; ii++) {
    TacsScalar diff = dfdxTest[ii] - dfdxTestRef[ii];
    if (fabs(REAL(diff)) > 1e-5) {
      noDiff = false;
    }
    printf(
        "dfdxTest[%d] = % 011.7e, dfdxTestRef[%d] = % 011.7e, diff = % "
        "011.7e\n",
        ii, REAL(dfdxTest[ii]), ii, REAL(dfdxTestRef[ii]), REAL(diff));
  }

  // --- Panel and Stiffener ---
  TacsScalar C[NUM_C];
  TacsScalar CRef[NUM_C];
  TacsScalar CDiff[NUM_C];
  memset(C, 0, sizeof(TacsScalar) * NUM_C);
  memset(CRef, 0, sizeof(TacsScalar) * NUM_C);
  memset(CDiff, 0, sizeof(TacsScalar) * NUM_C);

  TacsScalar Z = Zref - 0.5 * panelThickness;
  CRef[0] = pInv * Aref * EStiff;
  CRef[5] = pInv * Aref * GStiff * kcorr / 4.0;
  CRef[6] = pInv * Aref * EStiff * Z;
  CRef[11] = pInv * kcorr * GStiff * Aref * Z / 4.0;
  CRef[12] = pInv * EStiff * (Iref + Aref * Z * Z);
  CRef[17] = pInv * GStiff * (Jref + kcorr * Aref * Z * Z) / 4.0;
  CRef[20] = pInv * kcorr * GStiff * Aref;
  CRef[21] = 0.5 * KDRILL * CRef[20];

  con.evalTangentStiffness(0, NULL, NULL, C);
  noDiff = true;
  for (int ii = 0; ii < NUM_C; ii++) {
    CRef[ii] += CPanelOnlyRef[ii];
    CDiff[ii] = C[ii] - CRef[ii];
    if (fabs(REAL(CDiff[ii])) > 1e-10) {
      noDiff = false;
    }
  }

  printf("\n\nPanel and Stiffener Stiffness Matrix Ref:\n");
  printStiffnessMatrix(CRef);
  printf("\nPanel and Stiffener Stiffness Matrix:\n");
  printStiffnessMatrix(C);
  printf("\nPanel and Stiffener Stiffness Diff:\n");
  printStiffnessMatrix(CDiff);

  //
  printf("\n\n============================================================\n");
  if (noDiff) {
    printf("Panel and Stiffener Stiffness Matrix is correct!\n");
  } else {
    printf("Panel and Stiffener Stiffness Matrix is incorrect!\n");
  }
  printf("============================================================\n");

  TacsScalar stress[NUM_STRESS];
  TacsScalar stressRef[NUM_STRESS];
  TacsScalar stressDiff[NUM_STRESS];
  memset(stress, 0, sizeof(TacsScalar) * NUM_STRESS);
  memset(stressRef, 0, sizeof(TacsScalar) * NUM_STRESS);
  memset(stressDiff, 0, sizeof(TacsScalar) * NUM_STRESS);
  computeStress(CRef, strain, stressRef);
  con.evalStress(0, NULL, NULL, strain, stress);

  noDiff = true;
  printf("\n\nPanel and Stiffener Stress:\n");
  for (int ii = 0; ii < NUM_STRESS; ii++) {
    stressDiff[ii] = stress[ii] - stressRef[ii];
    if (fabs(REAL(stressDiff[ii])) > 1e-10) {
      noDiff = false;
    }
    printf("stressRef[%d] = % 011.7e, stress[%d] = % 011.7e, diff = % 011.7e\n",
           ii, REAL(stressRef[ii]), ii, REAL(stress[ii]), REAL(stressDiff[ii]));
  }

  //
  printf("\n\n============================================================\n");
  if (noDiff) {
    printf("Panel and Stiffener Stress is correct!\n");
  } else {
    printf("Panel and Stiffener Stress is incorrect!\n");
  }
  printf("============================================================\n");

  return 0;
}

// =============================================================================
// Function definitions
// =============================================================================

void printStiffnessMatrix(TacsScalar *C) {
  TacsScalar *A = &C[0];
  TacsScalar *B = &C[6];
  TacsScalar *D = &C[12];
  TacsScalar *As = &C[18];
  TacsScalar drill = C[21];

  printf("[\n");
  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      REAL(A[0]), REAL(A[1]), REAL(A[2]), REAL(B[0]), REAL(B[1]), REAL(B[2]),
      0., 0., 0.);
  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      REAL(A[1]), REAL(A[3]), REAL(A[4]), REAL(B[1]), REAL(B[3]), REAL(B[4]),
      0., 0., 0.);
  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      REAL(A[2]), REAL(A[4]), REAL(A[5]), REAL(B[2]), REAL(B[4]), REAL(B[5]),
      0., 0., 0.);

  printf(
      "--------------------------------------------------------------------"
      "----"
      "--------------------------------------------------------\n");

  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      REAL(B[0]), REAL(B[1]), REAL(B[2]), REAL(D[0]), REAL(D[1]), REAL(D[2]),
      0., 0., 0.);
  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      REAL(B[1]), REAL(B[3]), REAL(B[4]), REAL(D[1]), REAL(D[3]), REAL(D[4]),
      0., 0., 0.);
  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      REAL(B[2]), REAL(B[4]), REAL(B[5]), REAL(D[2]), REAL(D[4]), REAL(D[5]),
      0., 0., 0.);

  printf(
      "--------------------------------------------------------------------"
      "----"
      "--------------------------------------------------------\n");

  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      0., 0., 0., 0., 0., 0., REAL(As[0]), REAL(As[1]), 0.);
  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      0., 0., 0., 0., 0., 0., REAL(As[1]), REAL(As[2]), 0.);

  printf(
      "--------------------------------------------------------------------"
      "----"
      "--------------------------------------------------------\n");

  printf(
      "[% 03.5e, % 03.5e, % 03.5e | % 03.5e, % 03.5e, % 03.5e | % 03.5e, % "
      "03.5e | % 03.5e]\n",
      0., 0., 0., 0., 0., 0., 0., 0., REAL(drill));
  printf("]\n");
}

void computeStress(TacsScalar C[], const TacsScalar e[], TacsScalar stress[]) {
  TacsScalar *A = &C[0];
  TacsScalar *B = &C[6];
  TacsScalar *D = &C[12];
  TacsScalar *As = &C[18];
  TacsScalar drill = C[21];

  stress[0] = A[0] * e[0] + A[1] * e[1] + A[2] * e[2] + B[0] * e[3] +
              B[1] * e[4] + B[2] * e[5];
  stress[1] = A[1] * e[0] + A[3] * e[1] + A[4] * e[2] + B[1] * e[3] +
              B[3] * e[4] + B[4] * e[5];
  stress[2] = A[2] * e[0] + A[4] * e[1] + A[5] * e[2] + B[2] * e[3] +
              B[4] * e[4] + B[5] * e[5];

  stress[3] = B[0] * e[0] + B[1] * e[1] + B[2] * e[2] + D[0] * e[3] +
              D[1] * e[4] + D[2] * e[5];
  stress[4] = B[1] * e[0] + B[3] * e[1] + B[4] * e[2] + D[1] * e[3] +
              D[3] * e[4] + D[4] * e[5];
  stress[5] = B[2] * e[0] + B[4] * e[1] + B[5] * e[2] + D[2] * e[3] +
              D[4] * e[4] + D[5] * e[5];

  stress[6] = As[0] * e[6] + As[1] * e[7];
  stress[7] = As[1] * e[6] + As[2] * e[7];

  stress[8] = drill * e[8];
}
