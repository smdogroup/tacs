/*
========================================================================
Gaussian Process Model for TACS Buckling Constraints
========================================================================
@File   :   GaussianProcessModel.h
@Date   :   2024/05/10
@Author :   Sean Phillip Engelstad
@Description : Use Gaussian Processes for machine learning techniques to interpolate
learned buckling constraints on a training dataset to new test data points. This
approach is implemented into the TACSGPBladeStiffenedShellConstitutive class in TACS
for more physically accurate buckling constraints of stiffened panels.
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
  GaussianProcessModel(int n_train, int n_param, TacsScalar Xtrain[], TacsScalar alpha[]);
  ~GaussianProcessModel();

  // predict the test data at a single point using a matrix-vector product
  // this is the mean test data prediction. The offline training beforehand trains the mean
  // surface of the training set to match mu - 3 * sigma, a bounding curve for buckling 
  // so that way we don't need to invert a covariance matrix each time.
  // here Xtest is one point prediction with an array of length n_param
  TacsScalar predictMeanTestData(const TacsScalar* Xtest);
  TacsScalar predictMeanTestDataSens(const TacsScalar* Xtest, TacsScalar* Xtestsens);

  static inline TacsScalar soft_relu(TacsScalar x, TacsScalar rho) {
    return 1.0 / rho * log(1 + exp(rho * x));
  };

protected:
    // virtual functions for the kernel definition and its sensitivity
    virtual TacsScalar kernel(const TacsScalar* Xtest, const TacsScalar* Xtrain);
    virtual void kernelSens(const TacsScalar* Xtest, const TacsScalar* Xtrain, TacsScalar alpha, TacsScalar* Xtestsens);

    int n_train;
    int n_param;
    // rank 1-tensor of length [n_param*n_train] [X]
    // if each point of Xtrain has data [rho0, xi, gamma, delta, zeta] with n_Train=5 then
    // the entries are basically [rho01, xi1, gamma1, delta1, zeta1, rho02, xi2, gamma2, delta2, zeta2, ..., zetaN]
    TacsScalar *Xtrain; 
    TacsScalar *alpha;
};

class AxialGaussianProcessModel : GaussianProcessModel {
public:  
  AxialGaussianProcessModel(int n_train, const TacsScalar Xtrain[], const TacsScalar alpha[]) : GaussianProcessModel(n_train, N_PARAM, Xtrain, alpha) {};
  AxialGaussianProcessModel(); // empty constructor for debugging purposes
  ~AxialGaussianProcessModel() : ~GaussianProcessModel() {};
  void setDefaultHyperParameters();
protected:
    // here Xtest, Xtrain are each length 5 arrays (N_PARAM) [just uses one train and one test point here]
    // these are overwritten from subclass
    TacsScalar kernel(const TacsScalar* Xtest, const TacsScalar* Xtrain) {};
    void kernelSens(const TacsScalar* Xtest, const TacsScalar* Xtrain, TacsScalar alpha, TacsScalar* Xtestsens);
    
    // set the default hyperparameters of the model
    const TacsScalar S1, S2, c, L1, S4, S5, L2, alpha1, L3, S6;

    // there are 5 parameters [rho_0, xi, gamma, zeta] for the axial model
    static int N_PARAM = 4;
};

class ShearGaussianProcessModel : AxialGaussianProcessModel {
public:  
  ShearGaussianProcessModel(int n_train, const TacsScalar Xtrain[], const TacsScalar alpha[]) : GaussianProcessModel(n_train, N_PARAM, Xtrain, alpha);
  ShearGaussianProcessModel(); // empty constructor for debugging purposes
  ~ShearGaussianProcessModel();
  void setdefaultHyperParameters();
protected:
    // here Xtest, Xtrain are each length 5 arrays (N_PARAM) [just uses one train and one test point here]
    TacsScalar kernel(const TacsScalar* Xtest, const TacsScalar* Xtrain) {};
    void kernelSens(const TacsScalar* Xtest, const TacsScalar* Xtrain, TacsScalar alpha, TacsScalar* Xtestsens);
    
    // set the default hyperparameters of the model
    const TacsScalar S1, S2, c, L1, S4, S5, L2, alpha1, L3, S6;

    // there are 5 parameters [rho_0, xi, gamma, zeta] for the axial model
    static int N_PARAM = 4;
};

class CripplingGaussianProcessModel : GaussianProcessModel {
public:  
  CripplingGaussianProcessModel(int n_train, const TacsScalar Xtrain[], const TacsScalar alpha[]) : GaussianProcessModel(n_train, N_PARAM, Xtrain, alpha);
  CripplingGaussianProcessModel(); // empty constructor for debugging purposes
  ~CripplingGaussianProcessModel();

  void setdefaultHyperParameters() {
    S1 = 1e-1;
    S2 = 3e-1;
    c = -1;
    L1 = 0.2;
    S4 = 1.0;
    S5 = 1.0;
    L2 = 0.3;
    alpha1 = 2.0;
    L3 = 4.0;
    S6 = 1.0;
  };

protected:
    // here Xtest, Xtrain are each length 5 arrays (N_PARAM) [just uses one train and one test point here]
    TacsScalar kernel(const TacsScalar* Xtest, const TacsScalar* Xtrain) {};
    TacsScalar kernelSens(const TacsScalar* Xtest, const TacsScalar* Xtrain, TacsScalar* Xtestsens);
    
    // set the default hyperparameters of the model
    const TacsScalar S1, S2, c, L1, S4, S5, L2, alpha1, L3, S6;

    // there are 5 parameters [rho_0, xi, gen_eps, zeta] for the axial model
    static int N_PARAM = 3;
};