#include "TACSGaussianProcessModel.h"

TACSGaussianProcessModel::TACSGaussianProcessModel(int n_train, int n_param,
                                                   const TacsScalar Xtrain[],
                                                   const TacsScalar alpha[]) {
  // constructor for the base TACSGaussianProcessModel class
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

TACSGaussianProcessModel::~TACSGaussianProcessModel() {
  // destructor for the base TACSGaussianProcessModel class
  delete[] this->Xtrain;
  this->Xtrain = nullptr;

  delete[] this->alpha;
  this->alpha = nullptr;
}

TacsScalar TACSGaussianProcessModel::predictMeanTestData(
    const TacsScalar* Xtest) {
  // Xtest is an array of size n_param (for one test data point)
  // use the equation mean(Ytest) = cov(Xtest,X_train) @ alpha [this is a dot
  // product] where Ytest is a scalar
  TacsScalar Ytest = 0.0;

  // iterate over each training data point, get the cross-term covariance and
  // add the coefficient alpha for it
  for (int itrain = 0; itrain < n_train; itrain++) {
    TacsScalar* loc_Xtrain = &Xtrain[n_param * itrain];
    Ytest += kernel(Xtest, loc_Xtrain) * alpha[itrain];
  }
  return Ytest;
}

TacsScalar TACSGaussianProcessModel::predictMeanTestDataSens(
    const TacsScalar Ysens, const TacsScalar* Xtest, TacsScalar* Xtestsens) {
  // Xtest is an array of size n_param (for one test data point)
  // the sensitivity here is on log[nondim-params]
  // use the equation mean(Ytest) = cov(Xtest,X_train) @ alpha [this is a dot
  // product] where Ytest is a scalar
  TacsScalar Ytest = 0.0;

  // iterate over each training data point, get the cross-term covariance and
  // add the coefficient alpha for it
  memset(Xtestsens, 0, n_param * sizeof(TacsScalar));

  for (int itrain = 0; itrain < n_train; itrain++) {
    TacsScalar* loc_Xtrain = &Xtrain[n_param * itrain];
    Ytest += kernel(Xtest, loc_Xtrain) * alpha[itrain];

    // backwards propagate to the Xtestsens through the kernel computation
    kernelSens(Ysens * alpha[itrain], Xtest, loc_Xtrain, Xtestsens);
  }

  return Ytest;
}

void TACSAxialGaussianProcessModel::setDefaultHyperParameters() {
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

TacsScalar TACSAxialGaussianProcessModel::kernel(const TacsScalar* Xtest,
                                                 const TacsScalar* Xtrain) {
  // define the kernel function k(*,*) on training and testing points for one
  // training point the entries are [log(xi), log(rho_0), log(1+gamma),
  // log(zeta)]

  // log(xi) direction 0
  TacsScalar kernel0 = S1 * S1 + S2 * S2 * soft_relu(Xtest[0] - c, alpha1) *
                                     soft_relu(Xtrain[0] - c, alpha1);

  // log(rho0) direction 1
  TacsScalar d1 = Xtest[1] - Xtrain[1];
  TacsScalar fact1 = soft_relu(1 - soft_abs(Xtest[1], this->ks), this->ks);
  TacsScalar fact2 = soft_relu(1 - soft_abs(Xtrain[1], this->ks), this->ks);
  TacsScalar k1_term1 = exp(-0.5 * (d1 * d1 / L1 / L1)) * fact1 * fact2;
  TacsScalar rho0term3 =
      soft_relu(-Xtest[1], this->ks) * soft_relu(-Xtrain[1], this->ks);
  TacsScalar kernel1 = k1_term1 + S4 + S5 * rho0term3;

  // log(1+gamma) direction 2
  TacsScalar d2 = Xtest[2] - Xtrain[2];
  TacsScalar kernel2 = exp(-0.5 * d2 * d2 / L2 / L2);

  // log(zeta) direction 3
  TacsScalar d3 = Xtest[3] - Xtrain[3];
  TacsScalar kernel3 = S6 * exp(-0.5 * d3 * d3 / L3 / L3);

  return kernel0 * kernel1 * kernel2 + 2.0 * kernel3;
}

void TACSAxialGaussianProcessModel::kernelSens(const TacsScalar ksens,
                                               const TacsScalar* Xtest,
                                               const TacsScalar* Xtrain,
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
  TacsScalar fact1 = soft_relu(1 - soft_abs(Xtest[1], this->ks), this->ks);
  TacsScalar fact2 = soft_relu(1 - soft_abs(Xtrain[1], this->ks), this->ks);
  TacsScalar k1_term1 = exp(-0.5 * (d1 * d1 / L1 / L1)) * fact1 * fact2;
  TacsScalar rho0term3 =
      soft_relu(-Xtest[1], this->ks) * soft_relu(-Xtrain[1], this->ks);
  TacsScalar kernel1 = k1_term1 + S4 + S5 * rho0term3;

  // log(1+gamma) direction 2
  TacsScalar d2 = Xtest[2] - Xtrain[2];
  TacsScalar kernel2 = exp(-0.5 * d2 * d2 / L2 / L2);

  // log(zeta) direction 3
  TacsScalar d3 = Xtest[3] - Xtrain[3];
  TacsScalar kernel3 = S6 * exp(-0.5 * d3 * d3 / L3 / L3);
  // sensitivity section
  // ---------------------------------------------------
  TacsScalar* jacobian = new TacsScalar[4];

  // log(xi) direction 0
  TacsScalar kernel0sens = S2 * S2 * soft_relu_sens(Xtest[0] - c, alpha1) *
                           soft_relu(Xtrain[0] - c, alpha1);
  jacobian[0] = kernel0sens * kernel1 * kernel2;

  // log(rho_0) direction 1
  TacsScalar kernel1sens = 0.0;
  kernel1sens += k1_term1 * -d1 / L1 / L1;
  kernel1sens += k1_term1 / fact1 *
                 soft_relu_sens(1 - soft_abs(Xtest[1], this->ks), this->ks) *
                 -soft_abs_sens(Xtest[1], this->ks);
  kernel1sens += S5 * soft_relu_sens(-Xtest[1], this->ks) * -1.0 *
                 soft_relu(-Xtrain[1], this->ks);
  jacobian[1] = kernel0 * kernel1sens * kernel2;

  // log(1+gamma) direction 2
  TacsScalar kernel2sens = kernel2 * -d2 / L2 / L2;
  jacobian[2] = kernel0 * kernel1 * kernel2sens;

  // log(zeta) direction 3
  TacsScalar kernel3sens = kernel3 * -d3 / L3 / L3;
  jacobian[3] = kernel3sens * 2.0;

  // scale up the Xtestsens by the backpropagated values
  for (int ii = 0; ii < 4; ii++) {
    Xtestsens[ii] += ksens * jacobian[ii];
  }
}

TacsScalar TACSGaussianProcessModel::testAllGPTests(TacsScalar epsilon,
                                                    int printLevel) {
  // run all GP tests
  const int n_tests = 4;
  TacsScalar* relErrors = new TacsScalar[n_tests];

  relErrors[0] = test_soft_relu(epsilon);
  relErrors[1] = test_soft_abs(epsilon);
  relErrors[2] = testPredictMeanTestData(epsilon, printLevel);
  relErrors[3] = testKernelSens(epsilon, printLevel);

  // get max rel error among them
  TacsScalar maxRelError = 0.0;
  for (int i = 0; i < n_tests; i++) {
    if (relErrors[i] > maxRelError) {
      maxRelError = relErrors[i];
    }
  }

  // get max rel error among them
  if (printLevel != 0) {
    printf("\ntestAllGPtests full results::\n");
    printf("\ttest_soft_relu = %.4e\n", relErrors[0]);
    printf("\ttest_soft_abs = %.4e\n", relErrors[1]);
    printf("\ttestPredictMeanTestData = %.4e\n", relErrors[2]);
    printf("\ttestKernelSens = %.4e\n", relErrors[3]);
    printf("\tOverall max rel error = %.4e\n\n", maxRelError);
  }
  return maxRelError;
}

TacsScalar TACSGaussianProcessModel::testPredictMeanTestData(TacsScalar epsilon,
                                                             int printLevel) {
  // test the sensitivities of the kernel computation

  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = this->n_param;
  TacsScalar* p_input = new TacsScalar[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  TacsScalar* x0 = new TacsScalar[n_input];
  for (int i0 = 0; i0 < n_input; i0++) {
    x0[i0] = ((double)rand() / (RAND_MAX));
  }

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;

  TacsScalar* x = new TacsScalar[n_input];
  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  f0 = predictMeanTestData(x);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  f2 = predictMeanTestData(x);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // now perform the adjoint sensitivity
  TacsScalar* input_sens = new TacsScalar[n_input];
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i];
  }
  predictMeanTestDataSens(p_output, x, input_sens);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\t%s..testPredictMeanTestDataSens:\n", typeid(this).name());
    printf("\t\t adjDeriv = %.4e\n", adjTD);
    printf("\t\t centralDiff = %.4e\n", centralDiff);
    printf("\t\t rel error = %.4e\n", relError);
  }
  return relError;
}

TacsScalar TACSAxialGaussianProcessModel::testKernelSens(TacsScalar epsilon,
                                                         int printLevel) {
  // test the sensitivities of the kernel computation

  // perform complex-step or finite difference check (depending on the value of
  // _eps/epsilon) generate random input perturbation and output perturbation
  // test vectors
  const int n_input = 4;
  TacsScalar* p_input = new TacsScalar[n_input];
  for (int ii = 0; ii < n_input; ii++) {
    p_input[ii] = ((double)rand() / (RAND_MAX));
  }
  // temporarily only make p_input[0] active
  // for (int ii = 0; ii < n_input; ii++) {
  //   if (ii == 1) {
  //     p_input[ii] = 0.0;
  //   }
  // }

  TacsScalar p_output = ((double)rand() / (RAND_MAX));

  // compute initial values
  TacsScalar* x0 = new TacsScalar[n_input];
  x0[0] = 0.43243;  // log(xi)
  x0[1] = 1.64243;  // log(rho0)
  x0[2] = 0.12345;  // log(1+gamma)
  x0[3] = 4.13432;  // log(zeta)

  // perform central difference over rho_0 function on [D11,D22,a,b]
  TacsScalar f0, f1, f2;

  TacsScalar* x = new TacsScalar[n_input];
  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] - p_input[i] * epsilon;
  }
  TacsScalar* Xtrain = this->Xtrain;
  f0 = kernel(x, Xtrain);

  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i] + p_input[i] * epsilon;
  }
  f2 = kernel(x, Xtrain);

  TacsScalar centralDiff = p_output * (f2 - f0) / 2.0 / epsilon;

  // now perform the adjoint sensitivity
  TacsScalar* input_sens = new TacsScalar[n_input];
  memset(input_sens, 0, n_input * sizeof(TacsScalar));
  for (int i = 0; i < n_input; i++) {
    x[i] = x0[i];
  }
  kernelSens(p_output, x, Xtrain, input_sens);
  TacsScalar adjTD = 0.0;
  for (int j = 0; j < n_input; j++) {
    adjTD += input_sens[j] * p_input[j];
  }
  adjTD = TacsRealPart(adjTD);

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\t%s..testKernelSens:\n", typeid(this).name());
    printf("\t\t adjDeriv = %.4e\n", adjTD);
    printf("\t\t centralDiff = %.4e\n", centralDiff);
    printf("\t\t rel error = %.4e\n", relError);
  }
  return relError;
}