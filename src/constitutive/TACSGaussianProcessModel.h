/*
========================================================================
Gaussian Process Model for TACS Buckling Constraints
========================================================================
@File   :   TACSGaussianProcessModel.h
@Date   :   2024/05/10
@Author :   Sean Phillip Engelstad
@Description : Use Gaussian Processes for machine learning techniques to
interpolate learned buckling constraints on a training dataset to new test data
points. This approach is implemented into the
TACSGPBladeStiffenedShellConstitutive class in TACS for more physically accurate
buckling constraints of stiffened panels.
*/

#pragma once

// =============================================================================
// Extension Includes
// =============================================================================
#include "TACSObject.h"
#include "TacsUtilities.h"

// =============================================================================
// Class Declaration
// =============================================================================

class TACSGaussianProcessModel : public TACSObject {
 public:
  TACSGaussianProcessModel(int n_train, int n_param, const TacsScalar Xtrain[],
                           const TacsScalar alpha[], const TacsScalar theta[]);
  ~TACSGaussianProcessModel();

  // predict the test data at a single point using a matrix-vector product
  // this is the mean test data prediction. The offline training beforehand
  // trains the mean surface of the training set to match mu - 3 * sigma, a
  // bounding curve for buckling so that way we don't need to invert a
  // covariance matrix each time. here Xtest is one point prediction with an
  // array of length n_param
  TacsScalar predictMeanTestData(const TacsScalar* Xtest);
  TacsScalar predictMeanTestDataSens(const TacsScalar Ysens,
                                     const TacsScalar* Xtest,
                                     TacsScalar* Xtestsens);

  // TESTING SCRIPTS
  // ---------------
  TacsScalar testAllGPTests(TacsScalar epsilon, int printLevel);
  TacsScalar testPredictMeanTestData(TacsScalar epsilon, int printLevel);
  virtual TacsScalar testKernelSens(TacsScalar epsilon, int printLevel) {
    return 0.0;
  };

  static inline TacsScalar soft_relu(TacsScalar x, TacsScalar rho) {
    TacsScalar one = 1.0;
    return 1.0 / rho * log(one + exp(rho * x));
  };
  static inline TacsScalar soft_relu_sens(TacsScalar x, TacsScalar rho) {
    TacsScalar one = 1.0;
    return exp(rho * x) / (one + exp(rho * x));
  };
  static TacsScalar test_soft_relu(TacsScalar epsilon) {
    TacsScalar x = 1.0,
               rho = 1.0;  // very low rho for smoother function for deriv test
    TacsScalar f0 = soft_relu(x - epsilon, rho);
    TacsScalar f2 = soft_relu(x + epsilon, rho);
    TacsScalar centDiff = (f2 - f0) / 2.0 / epsilon;
    TacsScalar analyDeriv = soft_relu_sens(x, rho);
    TacsScalar relError = (analyDeriv - centDiff) / centDiff;
    relError = abs(TacsRealPart(relError));
    return relError;
  };

  static inline TacsScalar soft_abs(TacsScalar x, TacsScalar rho) {
    return 1.0 / rho * log(exp(-rho * x) + exp(rho * x));
  };
  static inline TacsScalar soft_abs_sens(TacsScalar x, TacsScalar rho) {
    return (exp(rho * x) - exp(-rho * x)) / (exp(-rho * x) + exp(rho * x));
  };
  static TacsScalar test_soft_abs(TacsScalar epsilon) {
    TacsScalar x = 1.0,
               rho = 1.0;  // very low rho for smoother function for deriv test
    TacsScalar f0 = soft_abs(x - epsilon, rho);
    TacsScalar f2 = soft_abs(x + epsilon, rho);
    TacsScalar centDiff = (f2 - f0) / 2.0 / epsilon;
    TacsScalar analyDeriv = soft_abs_sens(x, rho);
    TacsScalar relError = (analyDeriv - centDiff) / centDiff;
    relError = abs(TacsRealPart(relError));
    return relError;
  };

  // GETTERS AND SETTERS
  // -------------------

  int getNtrain() { return n_train; };
  int getNparam() { return n_param; };
  TacsScalar getKS() { return ks; };
  void setKS(TacsScalar ks) { this->ks = ks; };
  void setAlpha(TacsScalar* alpha) {this->alpha = alpha;};
  void setTheta(TacsScalar* theta) {this->theta = theta;};
  void getTrainingData(TacsScalar* Xtrain) {Xtrain = this->Xtrain;};
  void getTheta(TacsScalar* theta) {theta = this->theta;};
  // virtual functions for the kernel definition and its sensitivity
  virtual TacsScalar kernel(const TacsScalar* Xtest, const TacsScalar* Xtrain) {
    return 0.0;
  };

 protected:
  virtual void kernelSens(const TacsScalar ksens, const TacsScalar* Xtest,
                          const TacsScalar* Xtrain, TacsScalar* Xtestsens){};

  int n_train;
  int n_param;
  int n_theta = 14;
  // rank 1-tensor of length [n_param*n_train] [X]
  // if each point of Xtrain has data [rho0, xi, gamma, delta, zeta] with
  // n_Train=5 then the entries are basically [rho01, xi1, gamma1, delta1,
  // zeta1, rho02, xi2, gamma2, delta2, zeta2, ..., zetaN]
  TacsScalar* Xtrain;
  TacsScalar* alpha;
  TacsScalar* theta; // hyperparameters

  TacsScalar ks = 10.0;  // ks setting for smooth kernel functions
};

class TACSAxialGaussianProcessModel : public TACSGaussianProcessModel {
 public:
  TACSAxialGaussianProcessModel(int n_train, const TacsScalar Xtrain[],
                                const TacsScalar alpha[], const TacsScalar theta[])
      : TACSGaussianProcessModel(n_train, N_PARAM, Xtrain, alpha, theta){};
  ~TACSAxialGaussianProcessModel(){};
  TacsScalar testKernelSens(TacsScalar epsilon, int printLevel);

  // here Xtest, Xtrain are each length 5 arrays (N_PARAM) [just uses one train
  // and one test point here] these are overwritten from subclass
  TacsScalar kernel(const TacsScalar* Xtest, const TacsScalar* Xtrain);

 protected:
  
  void kernelSens(const TacsScalar ksens, const TacsScalar* Xtest,
                  const TacsScalar* Xtrain, TacsScalar* Xtestsens);

  // there are 4 parameters [log(xi), log(rho_0), log(1+gamma), log(zeta)] for
  // the axial model
  static const int N_PARAM = 4;
};

class TACSShearGaussianProcessModel : public TACSAxialGaussianProcessModel {
 public:
  TACSShearGaussianProcessModel(int n_train, const TacsScalar Xtrain[],
                                const TacsScalar alpha[], const TacsScalar theta[])
      : TACSAxialGaussianProcessModel(n_train, Xtrain, alpha, theta){};
  ~TACSShearGaussianProcessModel(){};
};

class TACSCripplingGaussianProcessModel : public TACSAxialGaussianProcessModel {
 public:
  TACSCripplingGaussianProcessModel(int n_train, const TacsScalar Xtrain[],
                                    const TacsScalar alpha[], const TacsScalar theta[])
      : TACSAxialGaussianProcessModel(n_train, Xtrain, alpha, theta){};
  ~TACSCripplingGaussianProcessModel(){};
};
