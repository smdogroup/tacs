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
  TACSGaussianProcessModel(int n_train, int n_param, bool affine,
                           const TacsScalar Xtrain[], const TacsScalar alpha[],
                           const TacsScalar theta[]);
  ~TACSGaussianProcessModel();

  // predict the test data at a single point using a matrix-vector product
  // this is the mean test data prediction. The offline training beforehand
  // trains the mean surface of the training set to match mu - 3 * sigma, a
  // bounding curve for buckling so that way we don't need to invert a
  // covariance matrix each time. here Xtest is one point prediction with an
  // array of length n_param

  /**
   * @brief predict the mean test data Ytest for one test data point
   *
   * @param Xtest the 1-tensor test data inputs [param1, param2, param3, param4]
   * @return the predicted scalar test data Ytest
   */
  TacsScalar predictMeanTestData(const TacsScalar *Xtest);

  /**
   * @brief backpropagate derivatives df/dYtest to df/dXtest for
   * predictMeanTestData
   *
   * @param Ysens the derivative df/dYtest as a scalar
   * @param Xtest the 1-tensor test data inputs [param1, param2, param3, param4]
   * @return the derivative df/dXtest as a 1-tensor
   */
  TacsScalar predictMeanTestDataSens(const TacsScalar Ysens,
                                     const TacsScalar *Xtest,
                                     TacsScalar *Xtestsens);

  // TESTING SCRIPTS
  // ---------------
  /**
   * @brief test all subtests in the GaussianProcess model
   *
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the maximum relative error among all GP subtests
   */
  TacsScalar testAllGPTests(TacsScalar epsilon, int printLevel);

  /**
   * @brief test the backpropagation of the predict mean test data and its sens
   * routine
   *
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the relative error for predictMeanTestData and its sens routine
   */
  TacsScalar testPredictMeanTestData(TacsScalar epsilon, int printLevel);

  /**
   * @brief test the backpropagation of the kernel() method and its sens routine
   *
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the relative error for the kernel() and kernelSens derivatives
   */
  virtual TacsScalar testKernelSens(TacsScalar epsilon, int printLevel) {
    return 0.0;
  };

  /**
   * @brief a differentiable form of the relu function
   *
   * @param x the input for soft_relu(x)
   * @param rho the smoothing parameter rho_KS
   * @return the soft_relu(x)
   */
  static inline TacsScalar soft_relu(TacsScalar x, TacsScalar rho) {
    TacsScalar one = 1.0;
    return 1.0 / rho * log(one + exp(rho * x));
  };

  /**
   * @brief Jacobian dsoft_relu(x)/dx of the soft_relu function
   *
   * @param x the input for soft_relu(x)
   * @param rho the smoothing parameter rho_KS
   * @return the jacobian dsoft_relu(x)/dx
   */
  static inline TacsScalar soft_relu_sens(TacsScalar x, TacsScalar rho) {
    TacsScalar one = 1.0;
    return exp(rho * x) / (one + exp(rho * x));
  };

  /**
   * @brief test the soft_relu jacobian
   *
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the relative error of the soft_relu() and soft_relu_sens() jacobian
   * routines compared to finite diff
   */
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

  /**
   * @brief a differentiable form of the absolute value function
   *
   * @param x the input for soft_abs(x)
   * @param rho the smoothing parameter rho_KS
   * @return the soft_abs(x)
   */
  static inline TacsScalar soft_abs(TacsScalar x, TacsScalar rho) {
    return 1.0 / rho * log(exp(-rho * x) + exp(rho * x));
  };

  /**
   * @brief Jacobian dsoft_abs(x)/dx of the soft_abs function
   *
   * @param x the input for soft_abs(x)
   * @param rho the smoothing parameter rho_KS
   * @return the jacobian dsoft_abs(x)/dx
   */
  static inline TacsScalar soft_abs_sens(TacsScalar x, TacsScalar rho) {
    return (exp(rho * x) - exp(-rho * x)) / (exp(-rho * x) + exp(rho * x));
  };

  /**
   * @brief test the soft_abs jacobian
   *
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the relative error of the soft_abs() and soft_abs_sens() jacobian
   * routines compared to finite diff
   */
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
  void setAlpha(TacsScalar *alpha) { this->alpha = alpha; };
  void setTheta(TacsScalar *theta) { this->theta = theta; };
  void getTrainingData(TacsScalar *Xtrain) { Xtrain = this->Xtrain; };
  void getTheta(TacsScalar *theta) { theta = this->theta; };

  // virtual functions for the kernel definition and its sensitivity

  /**
   * @brief virtual function for the kernel(x,x') for one test point and one
   * training data point each
   *
   * @param Xtest the test data point, rank 1-tensor of length n_param
   * @param Xtrain the training data point, rank 1-tensor of length n_param
   * @return the kernel value k(Xtest,Xtrain) which gives correlation between
   * these two points from our model
   */
  virtual TacsScalar kernel(const TacsScalar *Xtest,
                            const TacsScalar *Xtrain) = 0;

 protected:
  /**
   * @brief backpropagate derivatives of the kernel function to the Xtest input
   * (this is a virtual function here in the base class)
   *
   * @param ksens the derivative df/dkernel
   * @param Xtest the test data point, rank 1-tensor of length n_param
   * @param Xtrain the training data point, rank 1-tensor of length n_param
   * @return the derivatives of the Xtest input df/dXtest through the kernel
   */
  virtual void kernelSens(const TacsScalar ksens, const TacsScalar *Xtest,
                          const TacsScalar *Xtrain, TacsScalar *Xtestsens) = 0;

  int n_train;
  int n_param;
  int n_theta = 6;
  bool affine = false;
  // rank 1-tensor of length [n_param*n_train] [X]
  // if each point of Xtrain has data [rho0, xi, gamma, delta, zeta] with
  // n_Train=5 then the entries are basically [rho01, xi1, gamma1, delta1,
  // zeta1, rho02, xi2, gamma2, delta2, zeta2, ..., zetaN]
  TacsScalar *Xtrain;
  TacsScalar *alpha;
  TacsScalar *theta;  // hyperparameters

  // not using this ks anymore though.. it's a trained hyperparameter, so fixed
  TacsScalar ks;
  // TacsScalar ks = 10.0; // ks setting for smooth kernel functions
};

class TACSBucklingGaussianProcessModel : public TACSGaussianProcessModel {
 public:
  TACSBucklingGaussianProcessModel(int n_train, bool affine,
                                   const TacsScalar Xtrain[],
                                   const TacsScalar alpha[],
                                   const TacsScalar theta[])
      : TACSGaussianProcessModel(n_train, N_PARAM, affine, Xtrain, alpha,
                                 theta) {};
  ~TACSBucklingGaussianProcessModel() {};

  /**
   * @brief test the backpropagation of the kernel() method and its sens routine
   *
   * @param epsilon the step size for a finite difference or complex-step test
   * (multiply by 1j if complex-step)
   * @param printLevel an integer flag, with 0 to not print the test result to
   * terminal and 1 to print to terminal
   * @return the relative error for the kernel() and kernelSens derivatives
   */
  TacsScalar testKernelSens(TacsScalar epsilon, int printLevel) override;

  /**
   * @brief AxialGP kernel(x,x') for one test point and one training data point
   * each
   *
   * @param Xtest the test data point, rank 1-tensor of length 4
   * @param Xtrain the training data point, rank 1-tensor of length 4
   * @return the kernel value k(Xtest,Xtrain) which gives correlation between
   * these two points from our model
   */
  TacsScalar kernel(const TacsScalar *Xtest, const TacsScalar *Xtrain) override;

 protected:
  /**
   * @brief backpropagate derivatives of the kernel function to the Xtest input
   * for AxialGP
   *
   * @param ksens the derivative df/dkernel
   * @param Xtest the test data point, rank 1-tensor of length 4
   * @param Xtrain the training data point, rank 1-tensor of length 4
   * @return the derivatives of the Xtest input df/dXtest through the kernel
   */
  void kernelSens(const TacsScalar ksens, const TacsScalar *Xtest,
                  const TacsScalar *Xtrain, TacsScalar *Xtestsens) override;

  // there are 4 parameters [log(xi), log(rho_0), log(1+gamma), log(zeta)] for
  // the axial model
  static const int N_PARAM = 4;
};