/*
========================================================================
Gaussian Process Model for TACS Buckling Constraints
========================================================================
@File   :   GaussianProcessModel.h
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
#include "TacsUtilities.h"

// =============================================================================
// Class Declaration
// =============================================================================

class GaussianProcessModel {
 public:
  GaussianProcessModel(int n_train, int n_param, const TacsScalar Xtrain[],
                       const TacsScalar alpha[]);
  ~GaussianProcessModel();

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
  TacsScalar testAllGPTests(TacsScalar epsilon);
  TacsScalar testPredictMeanTestData(TacsScalar epsilon);
  virtual TacsScalar testKernelSens(TacsScalar epsilon) { return 0.0; };

  static inline TacsScalar soft_relu(TacsScalar x, TacsScalar rho) {
    return 1.0 / rho * log(1 + exp(rho * x));
  };
  static inline TacsScalar soft_relu_sens(TacsScalar x, TacsScalar rho) {
    return exp(rho * x) / (1 + exp(rho * x));
  }

  // GETTERS AND SETTERS
  // -------------------

  int getNtrain() { return n_train; }
  int getNparam() { return n_param; }

 protected:
  // virtual functions for the kernel definition and its sensitivity
  virtual TacsScalar kernel(const TacsScalar* Xtest, const TacsScalar* Xtrain) {
    return 0.0;
  };
  virtual void kernelSens(const TacsScalar ksens, const TacsScalar* Xtest,
                          const TacsScalar* Xtrain, TacsScalar* Xtestsens){};

  int n_train;
  int n_param;
  // rank 1-tensor of length [n_param*n_train] [X]
  // if each point of Xtrain has data [rho0, xi, gamma, delta, zeta] with
  // n_Train=5 then the entries are basically [rho01, xi1, gamma1, delta1,
  // zeta1, rho02, xi2, gamma2, delta2, zeta2, ..., zetaN]
  TacsScalar* Xtrain;
  TacsScalar* alpha;

 private:
  static const char* GPname;
};

class AxialGaussianProcessModel : public GaussianProcessModel {
 public:
  AxialGaussianProcessModel(int n_train, const TacsScalar Xtrain[],
                            const TacsScalar alpha[])
      : GaussianProcessModel(n_train, N_PARAM, Xtrain, alpha) {
    setDefaultHyperParameters();
  };
  ~AxialGaussianProcessModel(){};
  void setDefaultHyperParameters();

  TacsScalar testKernelSens(TacsScalar epsilon);

 protected:
  // here Xtest, Xtrain are each length 5 arrays (N_PARAM) [just uses one train
  // and one test point here] these are overwritten from subclass
  TacsScalar kernel(const TacsScalar* Xtest, const TacsScalar* Xtrain);
  void kernelSens(const TacsScalar ksens, const TacsScalar* Xtest,
                  const TacsScalar* Xtrain, TacsScalar* Xtestsens);

  // set the default hyperparameters of the model
  TacsScalar S1, S2, c, L1, S4, S5, L2, alpha1, L3, S6;

  // there are 4 parameters [log(xi), log(rho_0), log(1+gamma), log(zeta)] for
  // the axial model
  static const int N_PARAM = 4;

 private:
  static const char* GPname;
};

class ShearGaussianProcessModel : public AxialGaussianProcessModel {
 public:
  ShearGaussianProcessModel(int n_train, const TacsScalar Xtrain[],
                            const TacsScalar alpha[])
      : AxialGaussianProcessModel(n_train, Xtrain, alpha){};
  ~ShearGaussianProcessModel(){};
  // void setdefaultHyperParameters();
  //  protected:
  // set the default hyperparameters of the model
  // for now just use the same routine as the axial one
  // const TacsScalar S1, S2, c, L1, S4, S5, L2, alpha1, L3, S6;
 private:
  static const char* GPname;
};

class CripplingGaussianProcessModel : public AxialGaussianProcessModel {
 public:
  CripplingGaussianProcessModel(int n_train, const TacsScalar Xtrain[],
                                const TacsScalar alpha[])
      : AxialGaussianProcessModel(n_train, Xtrain, alpha){};
  ~CripplingGaussianProcessModel(){};
  // void setdefaultHyperParameters();
  //  protected:
  // set the default hyperparameters of the model
  // for now just use the same routine as the axial one
  // const TacsScalar S1, S2, c, L1, S4, S5, L2, alpha1, L3, S6;
 private:
  static const char* GPname;
};
