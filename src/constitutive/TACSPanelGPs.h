/*
========================================================================
Gaussian Process Model for TACS Buckling Constraints
========================================================================
@File   :   TACSPanelGPs.h
@Date   :   2024/05/23
@Author :   Sean Phillip Engelstad
@Description : Container for the Axial, Shear, and Crippling GPs for the
GPBladeShellConstitutive class. There is one of these container objects per TACS
component (usually corresponds to each panel in the structure). The container
objects compute the GP test values once for forward and adjoint and then return
them again for each constitutive object within a TACS component (for huge
speedup in runtime).
*/

#pragma once

// =============================================================================
// Extension Includes
// =============================================================================
#include "TACSGaussianProcessModel.h"
#include "TACSObject.h"
#include "TacsUtilities.h"

// =============================================================================
// Class Declaration
// =============================================================================

class TACSPanelGPs : public TACSObject {
 public:
  /*
   * @param TACSAxialGaussianProcessModel an axial gaussian process model (if
   * null closed-form solution is used)
   * @param TACSShearGaussianProcessModel a shear gaussian process model (if
   * null closed-form solution is used)
   * @param TACSCripplingGaussianProcessModel a crippling gaussian process model
   * (if null closed-form solution is used) buckling constraints height
   */
  TACSPanelGPs(TACSAxialGaussianProcessModel* axialGP,
               TACSShearGaussianProcessModel* shearGP,
               TACSCripplingGaussianProcessModel* cripplingGP, bool saveData);
  ~TACSPanelGPs();

  // predict the test data and sens using the GPs
  // prediction indices to help retrieve the saved values.
  // 0 - axial global, 1 - axial local,
  // 2 - shear global, 3 - shear local
  // 4 - crippling
  //    If the values aren't saved yet the computation is performed
  //      and the saved flags updated.
  TacsScalar predictMeanTestData(int predInd, const TacsScalar* Xtest);
  void predictMeanTestDataSens(int predInd, const TacsScalar Ysens,
                               const TacsScalar* Xtest, TacsScalar* Xtestsens);

  // reset the saved forward and adjoint data upon a
  //    setVariables call into the constituive objects
  void resetSavedData();

  // GETTERS AND SETTERS
  // -------------------

  TACSAxialGaussianProcessModel* getAxialGP() { return this->axialGP; }
  TACSShearGaussianProcessModel* getShearGP() { return this->shearGP; }
  TACSCripplingGaussianProcessModel* getCripplingGP() {
    return this->cripplingGP;
  }

 protected:
  const int n_save = 5;
  int n_save_adj = 20;  // 5 * n_params rn
  bool saveData;

  // saved forward data in this class
  TacsScalar* savedYtest;
  bool* savedForward;

  // saved adjoint data in this class
  TacsScalar* savedJacobians;
  bool* savedAdjoint;

  // stored GP model pointers
  TACSAxialGaussianProcessModel* axialGP;
  TACSShearGaussianProcessModel* shearGP;
  TACSCripplingGaussianProcessModel* cripplingGP;
};
