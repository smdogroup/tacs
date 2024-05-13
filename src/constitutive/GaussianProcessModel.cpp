#include "GaussianProcessModel.h"

GaussianProcessModel::GaussianProcessModel(int n_train, int n_param,
                                           const TacsScalar Xtrain[],
                                           const TacsScalar alpha[]) {
  // constructor for the base GaussianProcessModel class
  this->n_train = n_train;
  this->n_param = n_param;

  // create a new copy of the Xtrain, alpha arrays on this process
  this->Xtrain = new TacsScalar[n_train * n_param];
  for (int ii = 0; ii < n_train * n_param; ii++) {
    this->Xtrain[ii] = Xtrain[ii];
  }

  this->alpha = new TacsScalar[n_train];
  for (int jj = 0; jj < n_train; jj++) {
    this->alpha[jj] = alpha[jj];
  }
}

GaussianProcessModel::~GaussianProcessModel() {
  // destructor for the base GaussianProcessModel class
  delete[] this->Xtrain;
  this->Xtrain = nullptr;

  delete[] this->alpha;
  this->alpha = nullptr;
}

TacsScalar GaussianProcessModel::predictMeanTestData(const TacsScalar* Xtest) {
  // Xtest is an array of size n_param (for one test data point)
  // use the equation mean(Ytest) = cov(Xtest,X_train) @ alpha [this is a dot
  // product] where Ytest is a scalar
  TacsScalar Ytest = 0.0;

  // iterate over each training data point, get the cross-term covariance and
  // add the coefficient alpha for it
  for (int itrain = 0; itrain < n_train; itrain++) {
    const TacsScalar* loc_Xtrain = &Xtrain[n_param * itrain];
    Ytest += kernel(Xtest, loc_Xtrain) * alpha[itrain];
  }
  return Ytest;
}

TacsScalar GaussianProcessModel::predictMeanTestDataSens(
    const TacsScalar* Xtest, TacsScalar* Xtestsens) {
  // Xtest is an array of size n_param (for one test data point)
  // the sensitivity here is on log[nondim-params]
  // use the equation mean(Ytest) = cov(Xtest,X_train) @ alpha [this is a dot
  // product] where Ytest is a scalar
  TacsScalar Ytest = 0.0;
  memset(Xtestsens, 0, n_param * sizeof(TacsScalar));

  // iterate over each training data point, get the cross-term covariance and
  // add the coefficient alpha for it
  for (int itrain = 0; itrain < n_train; itrain++) {
    const TacsScalar* loc_Xtrain = &Xtrain[n_param * itrain];
    Ytest += kernel(Xtest, loc_Xtrain) * alpha[itrain];
    // add the kernel sensitivity for this training point
    kernelSens(Xtest, loc_Xtrain, alpha[itrain], Xtestsens);
  }
  return Ytest;
}

AxialGaussianProcessModel::~AxialGaussianProcessModel() {
  // destructor for the AxialGaussianProcessModel class
  GaussianProcessModel::~GaussianProcessModel();
}

ShearGaussianProcessModel::~ShearGaussianProcessModel() {
  // destructor for the AxialGaussianProcessModel class
  GaussianProcessModel::~GaussianProcessModel();
}

CripplingGaussianProcessModel::~CripplingGaussianProcessModel() {
  // destructor for the AxialGaussianProcessModel class
  GaussianProcessModel::~GaussianProcessModel();
}

void AxialGaussianProcessModel::setDefaultHyperParameters() {
  this->S1 = 1e-1;
  this->S2 = 3e-1;
  this->c = -1;
  this->L1 = 0.2;
  this->S4 = 1.0;
  this->S5 = 1.0;
  this->L2 = 0.3;
  this->alpha1 = 2.0;
  this->L3 = 4.0;
  this->S6 = 1.0;
}

TacsScalar AxialGaussianProcessModel::kernel(const TacsScalar* Xtest,
                                             const TacsScalar* Xtrain) {
  // define the kernel function k(*,*) on training and testing points for one
  // training point the entries are [log(xi), log(rho_0), log(1+gamma),
  // log(zeta)]

  // log(xi) direction 0
  TacsScalar kernel0 = S1 * S1 + S2 * S2 * soft_relu(Xtest[0] - c, alpha1) *
                                     soft_relu(Xtrain[0] - c, alpha1);

  // log(rho0) direction 1
  TacsScalar d1 = Xtest[1] - Xtrain[1];
  TacsScalar fact1 = soft_relu(1 - soft_relu(Xtest[1], 1.0), 10);
  TacsScalar fact2 = soft_relu(1 - soft_relu(Xtrain[1], 1.0), 10);
  TacsScalar rho0term3 = soft_relu(-Xtest[1], 10) * soft_relu(-Xtrain[1], 10);
  TacsScalar kernel1 =
      exp(-0.5 * (d1 * d1 / L1 / L1)) * fact1 * fact2 + S4 + S5 * rho0term3;

  // log(1+gamma) direction 2
  TacsScalar d2 = Xtest[2] - Xtrain[2];
  TacsScalar kernel2 = exp(-0.5 * d2 * d2 / L2 / L2);

  // log(zeta) direction 3
  TacsScalar d3 = Xtest[3] - Xtrain[3];
  TacsScalar kernel3 = S6 * exp(-0.5 * d3 * d3 / L3 / L3);

  return kernel0 * kernel1 * kernel2 + 2.0 * kernel3;
}

void AxialGaussianProcessModel::kernelSens(const TacsScalar* Xtest,
                                           const TacsScalar* Xtrain,
                                           TacsScalar alpha,
                                           TacsScalar* Xtestsens) {
  // add into the Xtestsens (don't reset to zero) for x_test = log[nondim
  // params] vector

  // forward analysis section
  // -----------------------------------
  // log(xi) direction 0
  TacsScalar kernel0 = S1 * S1 + S2 * S2 * soft_relu(Xtest[0] - c, alpha1) *
                                     soft_relu(Xtrain[0] - c, alpha1);

  // log(rho0) direction 1
  TacsScalar d1 = Xtest[1] - Xtrain[1];
  TacsScalar fact1 = soft_relu(1 - soft_relu(Xtest[1], 1.0), 10);
  TacsScalar fact2 = soft_relu(1 - soft_relu(Xtrain[1], 1.0), 10);
  TacsScalar rho0term3 = soft_relu(-Xtest[1], 10) * soft_relu(-Xtrain[1], 10);
  TacsScalar kernel1 =
      exp(-0.5 * (d1 * d1 / L1 / L1)) * fact1 * fact2 + S4 + S5 * rho0term3;

  // log(1+gamma) direction 2
  TacsScalar d2 = Xtest[2] - Xtrain[2];
  TacsScalar kernel2 = exp(-0.5 * d2 * d2 / L2 / L2);

  // log(zeta) direction 3
  TacsScalar d3 = Xtest[3] - Xtrain[3];
  TacsScalar kernel3 = S6 * exp(-0.5 * d3 * d3 / L3 / L3);
  // sensitivity section
  // ---------------------------------------------------
  // log(xi) direction 0
  TacsScalar kernel0sens = S2 * S2 * soft_relu_sens(Xtest[0] - c, alpha1) *
                           soft_relu(Xtrain[0] - c, alpha1);
  Xtestsens[0] += kernel0sens * kernel1 * kernel2;

  // log(rho_0) direction 1
  TacsScalar kernel1sens = 0.0;
  kernel1sens += kernel1 * -d1 / L1 / L1;
  kernel1sens += exp(-0.5 * (d1 * d1 / L1 / L1)) *
                 soft_relu_sens(1 - soft_relu(Xtest[1], 1.0), 10) *
                 -soft_relu_sens(Xtest[1], 1.0) * fact2;
  kernel1sens += S5 * soft_relu_sens(-Xtest[1], 10) * soft_relu(-Xtrain[1], 10);
  Xtestsens[1] += kernel0 * kernel1sens * kernel2;

  // log(1+gamma) direction 2
  TacsScalar kernel2sens = kernel2 * -d2 / L2 / L2;
  Xtestsens[2] += kernel0 * kernel1 * kernel2sens;

  // log(zeta) direction 3
  TacsScalar kernel3sens = kernel3 * -d3 / L3 / L3;
  Xtestsens[3] += kernel3sens * 2.0;
}