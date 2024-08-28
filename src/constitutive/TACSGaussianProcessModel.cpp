#include "TACSGaussianProcessModel.h"

TACSGaussianProcessModel::TACSGaussianProcessModel(int n_train, int n_param,
                                                   const TacsScalar Xtrain[],
                                                   const TacsScalar alpha[],
                                                   const TacsScalar theta[]) {
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

  this->theta = new TacsScalar[n_theta];
  for (int kk = 0; kk < n_theta; kk++) {
    this->theta[kk] = theta[kk];
  }

  this->ks = 10.0;

}

TACSGaussianProcessModel::~TACSGaussianProcessModel() {
  // destructor for the base TACSGaussianProcessModel class
  delete[] this->Xtrain;
  this->Xtrain = nullptr;

  delete[] this->alpha;
  this->alpha = nullptr;

  delete[] this->theta;
  this->theta = nullptr;
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
    // TacsScalar mkernel = kernel(Xtest, loc_Xtrain);
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

TacsScalar TACSAxialGaussianProcessModel::kernel(const TacsScalar* Xtest,
                                                 const TacsScalar* Xtrain) {
  // define the kernel function k(*,*) on training and testing points for one
  // training point the entries are [log(1+xi), log(rho_0), log(1+gamma),
  // log(1+10^3*zeta)]

  TacsScalar d0 = Xtest[0] - Xtrain[0];
  TacsScalar d1 = Xtest[1] - Xtrain[1];
  TacsScalar d2 = Xtest[2] - Xtrain[2];
  TacsScalar d3 = Xtest[3] - Xtrain[3];

  TacsScalar bilinear_kernel = this->theta[0] + 
    soft_relu(-Xtest[1], this->ks) * soft_relu(-Xtrain[1], this->ks);
  TacsScalar SE_kernel = this->theta[1] * exp(-0.5 * (pow(d1/theta[2], 2.0) + pow(d2 / theta[3], 2.0)));
  TacsScalar SE_kernel2 = this->theta[4] * exp(-0.5 * (pow(d0/theta[5], 2.0) + pow(d3 / theta[6], 2.0)));
  TacsScalar SE_window = soft_relu(theta[7] - soft_abs(Xtest[1], this->ks), this->ks) * 
                         soft_relu(theta[7] - soft_abs(Xtrain[1], this->ks), this->ks) *
                         soft_relu(theta[8] - Xtest[2], this->ks) * 
                         soft_relu(theta[8] - Xtrain[2], this->ks);
  TacsScalar gamma_kernel = this->theta[9] + this->theta[10] * Xtest[2] * Xtrain[2];
  TacsScalar xi_kernel = this->theta[11] * Xtest[0] * Xtrain[0] + this->theta[12] * pow(Xtest[0] * Xtrain[0], 2.0);
  TacsScalar zeta_kernel = this->theta[13] * Xtest[3] * Xtrain[3] + this->theta[14] * pow(Xtest[3] * Xtrain[3], 2.0);

  return bilinear_kernel * gamma_kernel + SE_kernel * SE_window + SE_kernel2 + xi_kernel + zeta_kernel;
}

void TACSAxialGaussianProcessModel::kernelSens(const TacsScalar ksens,
                                               const TacsScalar* Xtest,
                                               const TacsScalar* Xtrain,
                                               TacsScalar* Xtestsens) {
  // add into the Xtestsens (don't reset to zero) for x_test = log[nondim
  // params] vector

  // forward analysis section
  // -----------------------------------
  // training point the entries are [log(1+xi), log(rho_0), log(1+gamma),
  // log(1+10^3*zeta)]

  TacsScalar d0 = Xtest[0] - Xtrain[0];
  TacsScalar d1 = Xtest[1] - Xtrain[1];
  TacsScalar d2 = Xtest[2] - Xtrain[2];
  TacsScalar d3 = Xtest[3] - Xtrain[3];

  TacsScalar bilinear_kernel = this->theta[0] + 
    soft_relu(-Xtest[1], this->ks) * soft_relu(-Xtrain[1], this->ks);
  TacsScalar SE_kernel = this->theta[1] * exp(-0.5 * (pow(d1/theta[2], 2.0) + pow(d2 / theta[3], 2.0)));
  TacsScalar SE_kernel2 = this->theta[4] * exp(-0.5 * (pow(d0/theta[5], 2.0) + pow(d3 / theta[6], 2.0)));
  TacsScalar SE_window = soft_relu(theta[7] - soft_abs(Xtest[1], this->ks), this->ks) * 
                         soft_relu(theta[7] - soft_abs(Xtrain[1], this->ks), this->ks) *
                         soft_relu(theta[8] - Xtest[3], this->ks) * 
                         soft_relu(theta[8] - Xtrain[3], this->ks);
  TacsScalar gamma_kernel = this->theta[9] + this->theta[10] * Xtest[2] * Xtrain[2];
  TacsScalar xi_kernel = this->theta[11] * Xtest[0] * Xtrain[0] + this->theta[12] * pow(Xtest[0] * Xtrain[0], 2.0);
  TacsScalar zeta_kernel = this->theta[13] * Xtest[3] * Xtrain[3] + this->theta[14] * pow(Xtest[3] * Xtrain[3], 2.0);

  // the main output of the kernel forward analysis
  TacsScalar output = bilinear_kernel * gamma_kernel + SE_kernel * SE_window + SE_kernel2 + xi_kernel + zeta_kernel;

  // sensitivity section
  // ---------------------------------------------------
  // hold derivatives w.r.t. Xtest[0], ..., Xtest[3] of the kernel
  TacsScalar* jacobian = new TacsScalar[4];

  // log(xi) direction 0
  TacsScalar xi_kernel_sens = this->theta[11] * Xtrain[0] + 2.0 * this->theta[12] * Xtest[0] * pow(Xtrain[0], 2.0);
  TacsScalar SE_kernel2_sens = SE_kernel2 * -d0 / pow(theta[5], 2.0);
  jacobian[0] = xi_kernel_sens + SE_kernel2_sens;

  // log(rho_0) direction 1
  TacsScalar BL_kernel_sens = soft_relu_sens(-Xtest[1], this->ks) * -1.0 *
                              soft_relu(-Xtrain[1], this->ks);
  TacsScalar SE_kernel_sens = SE_kernel * -d1 / pow(theta[2], 2.0);
  TacsScalar test_rho0_fact = soft_relu(theta[7] - soft_abs(Xtrain[1], this->ks), this->ks);
  TacsScalar test_rho0_fact_deriv = soft_relu_sens(theta[7] - soft_abs(Xtest[1], this->ks), this->ks) * 
                                    -soft_abs_sens(Xtest[1], this->ks);
  TacsScalar SE_window_sens = SE_window * test_rho0_fact_deriv / test_rho0_fact;
  jacobian[1] = BL_kernel_sens * gamma_kernel + SE_kernel_sens * SE_window + SE_kernel * SE_window_sens;

  // log(1+gamma) direction 2
  TacsScalar gamma_kernel_sens = theta[10] * Xtrain[2];
  TacsScalar SE_kernel_sens2 = SE_kernel * -d2 / pow(theta[3], 2.0);
  TacsScalar SE_window_sens2 = SE_window * soft_relu_sens(theta[8] - Xtest[2], this->ks) /
   soft_relu(theta[8] - Xtest[2], this->ks);

  // TacsScalar gamma_kernel_sens = 0.1 * Xtrain[2];
  jacobian[2] = bilinear_kernel * gamma_kernel_sens + SE_kernel_sens * SE_window + SE_kernel * SE_window_sens2;

  // log(zeta) direction 3
  TacsScalar SE_kernel2_sens2 = SE_kernel2 * -d3 / pow(theta[6], 2.0);
  TacsScalar zeta_kernel_sens = theta[13] * Xtrain[3] + theta[14] * Xtest[3] * pow(Xtest[3], 2.0);
  jacobian[3] = SE_kernel2_sens + zeta_kernel_sens;

  // scale up the Xtestsens by the backpropagated values
  for (int ii = 0; ii < 4; ii++) {
    Xtestsens[ii] += ksens * jacobian[ii];
  }

  // free memory for jacobian pointer
  delete[] jacobian;
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
    if (TacsRealPart(relErrors[i]) > TacsRealPart(maxRelError)) {
      maxRelError = relErrors[i];
    }
  }

  // get max rel error among them
  if (printLevel != 0) {
    printf("\ntestAllGPtests full results::\n");
    printf("\ttest_soft_relu = %.4e\n", TacsRealPart(relErrors[0]));
    printf("\ttest_soft_abs = %.4e\n", TacsRealPart(relErrors[1]));
    printf("\ttestPredictMeanTestData = %.4e\n", TacsRealPart(relErrors[2]));
    printf("\ttestKernelSens = %.4e\n", TacsRealPart(relErrors[3]));
    printf("\tOverall max rel error = %.4e\n\n", TacsRealPart(maxRelError));
  }

  // free memory for relErrors
  delete[] relErrors;

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
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }

  // free pointers
  delete[] p_input;
  delete[] x0;
  delete[] x;
  delete[] input_sens;

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

  // free pointers
  delete[] p_input;
  delete[] x0;
  delete[] x;
  delete[] input_sens;

  // compute relative error
  TacsScalar relError = abs((adjTD - centralDiff) / centralDiff);
  if (printLevel != 0) {
    printf("\t%s..testKernelSens:\n", typeid(this).name());
    printf("\t\t adjDeriv = %.4e\n", TacsRealPart(adjTD));
    printf("\t\t centralDiff = %.4e\n", TacsRealPart(centralDiff));
    printf("\t\t rel error = %.4e\n", TacsRealPart(relError));
  }
  return relError;
}
