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
  /**
   * TACSPanelGPs is a container object for the AxialGP, ShearGP, and
   * CripplingGP GaussianProcess ML models which saves and restores the buckling
   * predictions of each GP model for faster computation. Each GP model can be
   * null or not. If all GP models are null or not provided, then the
   * TACSGPBladeStiffenedShellConstitutive class uses closed-form buckling
   * predictions instead.
   *
   * @param TACSBucklingGaussianProcessModel an axial gaussian process model (if
   * null closed-form solution is used)
   * @param TACSBucklingGaussianProcessModel a shear gaussian process model (if
   * null closed-form solution is used)
   * @param TACSBucklingGaussianProcessModel a crippling gaussian process model
   * (if null closed-form solution is used) buckling constraints height
   * @param saveData a boolean flag for whether to save and restore data or not.
   */
  /*

   */
  TACSPanelGPs(TACSBucklingGaussianProcessModel *axialGP,
               TACSBucklingGaussianProcessModel *shearGP,
               TACSBucklingGaussianProcessModel *cripplingGP, bool saveData);
  ~TACSPanelGPs();

  /**
   * predict the mean test data Ytest for one test data point
   * using the GP models. The inputs are saved and restored in this method
   * for each kind of buckling prediction and selecting the appropriate GP for
   * each.
   * // 0 - axial global, 1 - axial local,
   * // 2 - shear global, 3 - shear local
   * // 4 - crippling
   * if the buckling predictions aren't saved yet, the computation is performed
   * and the saved data and flags are updated.
   *
   * @param predInd the index (see above) for which buckling prediction to make
   * @param Xtest the test data point, a rank 1-tensor of length 4
   * @return the Ytest mean prediction of the GP
   */
  TacsScalar predictMeanTestData(int predInd, const TacsScalar *Xtest);

  /**
   * derivatives of predictMeanTestData, which also saves and restores the
   * jacobians for faster derivative computations across an entire panel / TACS
   * component.
   *
   *
   * @param predInd the index (see above) for which buckling prediction to make
   * @param Ysens the derivative df/dYtest
   * @param Xtest the test data point, a rank 1-tensor of length 4
   * @return the derivative df/dXtest which we compute by jacobian product
   */
  void predictMeanTestDataSens(int predInd, const TacsScalar Ysens,
                               const TacsScalar *Xtest, TacsScalar *Xtestsens);

  /**
   * clear and reset all the saved data.
   * this also turns off the flags saying we have saved the data
   * and will trigger a new computation the next time the buckling predictions
   * are called.
   *
   * this is important in the constitutive model to ensure that the next time
   * DVs are updated or we have a new forward / adjoint analysis, we make new
   * buckling predictions and reset the saved data.
   */
  void resetSavedData();

  // GETTERS AND SETTERS
  // -------------------

  TACSBucklingGaussianProcessModel *getAxialGP() { return this->axialGP; }
  TACSBucklingGaussianProcessModel *getShearGP() { return this->shearGP; }
  TACSBucklingGaussianProcessModel *getCripplingGP() {
    return this->cripplingGP;
  }

 protected:
  const int n_save = 5;
  int n_save_adj = 20;  // 5 * n_params rn
  bool saveData;

  // saved forward data in this class
  TacsScalar *savedYtest;
  bool *savedForward;

  // saved adjoint data in this class
  TacsScalar *savedJacobians;
  bool *savedAdjoint;

  // stored GP model pointers
  TACSBucklingGaussianProcessModel *axialGP;
  TACSBucklingGaussianProcessModel *shearGP;
  TACSBucklingGaussianProcessModel *cripplingGP;
};
