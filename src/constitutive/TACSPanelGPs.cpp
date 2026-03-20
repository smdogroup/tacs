#include "TACSPanelGPs.h"

TACSPanelGPs::TACSPanelGPs(TACSBucklingGaussianProcessModel *axialGP,
                           TACSBucklingGaussianProcessModel *shearGP,
                           TACSBucklingGaussianProcessModel *cripplingGP,
                           bool saveData) {
  this->axialGP = axialGP;
  if (this->axialGP) {
    this->axialGP->incref();
  }
  this->shearGP = shearGP;
  if (this->shearGP) {
    this->shearGP->incref();
  }
  this->cripplingGP = cripplingGP;
  if (this->cripplingGP) {
    this->cripplingGP->incref();
  }

  // initialize all the saved data to zeros and the right length
  savedForward = new bool[n_save];
  savedYtest = new TacsScalar[n_save];
  savedAdjoint = new bool[n_save];
  savedJacobians = new TacsScalar[n_save_adj];

  this->saveData = saveData;
}

TACSPanelGPs::~TACSPanelGPs() {
  // destructor for the base TACSGaussianProcessModel class
  // destroy the gaussian process model objects if they exist
  if (this->axialGP) {
    // this object is shared, may not want to delete it (unless a local copy is
    // made)
    this->axialGP->decref();
    this->axialGP = nullptr;
  }

  if (this->shearGP) {
    // this object is shared, may not want to delete it (unless a local copy is
    // made)
    this->shearGP->decref();
    this->shearGP = nullptr;
  }

  if (this->cripplingGP) {
    // this object is shared, may not want to delete it (unless a local copy is
    // made)
    this->cripplingGP->decref();
    this->cripplingGP = nullptr;
  }

  // free pointers for saved data
  delete[] savedForward;
  delete[] savedYtest;
  delete[] savedAdjoint;
  delete[] savedJacobians;
}

void TACSPanelGPs::resetSavedData() {
  // goal here is to reset the saved data
  memset(savedYtest, 0.0, n_save * sizeof(TacsScalar));
  memset(savedJacobians, 0.0, n_save_adj * sizeof(TacsScalar));
  for (int i = 0; i < n_save; i++) {
    savedForward[i] = false;
    savedAdjoint[i] = false;
  }
}

TacsScalar TACSPanelGPs::predictMeanTestData(int predInd,
                                             const TacsScalar *Xtest) {
  // assume checking for input calling is in the other class for now
  // otherwise I would return garbage values..
  if (!savedForward[predInd] || !saveData) {
    // axial global or local options
    if (predInd == 0 || predInd == 1) {
      savedYtest[predInd] = this->axialGP->predictMeanTestData(Xtest);
    }
    // shear global or local options
    if (predInd == 2 || predInd == 3) {
      savedYtest[predInd] = this->shearGP->predictMeanTestData(Xtest);
    }
    // crippling GP
    if (predInd == 4) {
      savedYtest[predInd] = this->cripplingGP->predictMeanTestData(Xtest);
    }
    savedForward[predInd] = true;
  }
  return savedYtest[predInd];
}

void TACSPanelGPs::predictMeanTestDataSens(int predInd, const TacsScalar Ysens,
                                           const TacsScalar *Xtest,
                                           TacsScalar *Xtestsens) {
  // assume checking for input calling is in the other class for now
  // otherwise I would return garbage values..

  // just save the jacobians in each cell as the failure input derivatives
  // for backpropagation may not be the same

  // change this part with 4* to not be hardcoded later..
  TacsScalar *localJacobian = &savedJacobians[4 * predInd];

  if (!savedAdjoint[predInd] || !saveData) {
    // axial global or local options
    if (predInd == 0 || predInd == 1) {
      this->axialGP->predictMeanTestDataSens(1.0, Xtest, localJacobian);
    }
    // shear global or local options
    if (predInd == 2 || predInd == 3) {
      this->shearGP->predictMeanTestDataSens(1.0, Xtest, localJacobian);
    }
    // crippling GP
    if (predInd == 4) {
      this->cripplingGP->predictMeanTestDataSens(1.0, Xtest, localJacobian);
    }
    savedAdjoint[predInd] = true;
  }

  // now multiply the Ysens backpropagated derivative by the saved jacobian
  for (int i = 0; i < 4; i++) {
    Xtestsens[i] = Ysens * localJacobian[i];
  }
}
