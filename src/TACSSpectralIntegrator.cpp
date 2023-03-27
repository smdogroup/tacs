#include "TACSSpectralIntegrator.h"

TACSSpectralVec::TACSSpectralVec(int Nvecs, TACSAssembler *assembler) {
  N = Nvecs;

  vecs = new TACSBVec *[N];
  for (int i = 0; i < N; i++) {
    vecs[i] = assembler->createVec();
    vecs[i]->incref();
  }
}

TACSSpectralVec::TACSSpectralVec(int Nvecs, TACSMat *mat) {
  N = Nvecs;

  vecs = new TACSBVec *[N];
  for (int i = 0; i < N; i++) {
    TACSVec *vec = mat->createVec();
    vec->incref();

    TACSBVec *bvec = dynamic_cast<TACSBVec *>(vec);
    if (bvec) {
      vecs[i] = bvec;
    } else {
      vecs[i] = NULL;
      vec->decref();
    }
  }
}

TACSSpectralVec::~TACSSpectralVec() {
  for (int i = 0; i < N; i++) {
    if (vecs[i]) {
      vecs[i]->decref();
    }
  }
  delete[] vecs;
}

TacsScalar TACSSpectralVec::norm() {
  TacsScalar nrm = 0.0;
  for (int i = 0; i < N; i++) {
    nrm += vecs[i]->dot(vecs[i]);
  }
  return sqrt(nrm);
}

void TACSSpectralVec::scale(TacsScalar alpha) {
  for (int i = 0; i < N; i++) {
    vecs[i]->scale(alpha);
  }
}

TacsScalar TACSSpectralVec::dot(TACSVec *xvec) {
  TACSSpectralVec *x = dynamic_cast<TACSSpectralVec *>(xvec);
  if (x) {
    TacsScalar d = 0.0;
    for (int i = 0; i < N; i++) {
      d += vecs[i]->dot(x->getVec(i));
    }
    return d;
  }

  return 0.0;
}

void TACSSpectralVec::axpy(TacsScalar alpha, TACSVec *xvec) {
  TACSSpectralVec *x = dynamic_cast<TACSSpectralVec *>(xvec);
  if (x) {
    for (int i = 0; i < N; i++) {
      vecs[i]->axpy(alpha, x->getVec(i));
    }
  }
}

void TACSSpectralVec::copyValues(TACSVec *xvec) {
  TACSSpectralVec *x = dynamic_cast<TACSSpectralVec *>(xvec);
  if (x) {
    for (int i = 0; i < N; i++) {
      vecs[i]->copyValues(x->getVec(i));
    }
  }
}

void TACSSpectralVec::axpby(TacsScalar alpha, TacsScalar beta, TACSVec *xvec) {
  TACSSpectralVec *x = dynamic_cast<TACSSpectralVec *>(xvec);
  if (x) {
    for (int i = 0; i < N; i++) {
      vecs[i]->axpby(alpha, beta, x->getVec(i));
    }
  }
}

void TACSSpectralVec::zeroEntries() {
  for (int i = 0; i < N; i++) {
    vecs[i]->zeroEntries();
  }
}

TACSBVec *TACSSpectralVec::getVec(int index) {
  if (index >= 0 && index < N) {
    return vecs[index];
  }
  return NULL;
}

TACSLinearSpectralMat::TACSLinearSpectralMat(int N, const double *D,
                                             TACSMat *Hmat, TACSMat *Cmat)
    : N(N), D(D) {
  C = Cmat;
  H = Hmat;
  C->incref();
  H->incref();
  orient = TACS_MAT_NORMAL;

  TACSVec *vec = H->createVec();
  temp = dynamic_cast<TACSBVec *>(vec);
  if (temp) {
    temp->incref();
  } else {
    vec->decref();
  }
}

TACSLinearSpectralMat::~TACSLinearSpectralMat() {
  C->decref();
  H->decref();
  if (temp) {
    temp->decref();
  }
}

TACSVec *TACSLinearSpectralMat::createVec() {
  return new TACSSpectralVec(N, H);
}

void TACSLinearSpectralMat::mult(TACSVec *xvec, TACSVec *yvec) {
  TACSSpectralVec *x = dynamic_cast<TACSSpectralVec *>(xvec);
  TACSSpectralVec *y = dynamic_cast<TACSSpectralVec *>(yvec);

  if (x && y) {
    // Compute the contributions to the derivative
    for (int i = 0; i < N; i++) {
      // Start the derivative index from the second column of the derivative
      // operator
      temp->zeroEntries();
      const double *Dt = &D[(i + 1) * (N + 1)];
      for (int j = 0; j < N; j++) {
        temp->axpy(Dt[j + 1], x->getVec(j));
      }

      C->mult(temp, y->getVec(i));
    }

    // Add the diagonal contributions
    for (int i = 0; i < N; i++) {
      H->mult(x->getVec(i), temp);
      y->getVec(i)->axpy(1.0, temp);
    }
  }
}

void TACSLinearSpectralMat::getMat(TACSMat **Hmat, TACSMat **Cmat) {
  if (Hmat) {
    *Hmat = H;
  }
  if (Cmat) {
    *Cmat = C;
  }
}

void TACSLinearSpectralMat::setMatrixOrientation(MatrixOrientation matOr) {
  orient = matOr;
}

TACSSpectralIntegrator::TACSSpectralIntegrator(TACSAssembler *_assembler,
                                               double tfinal, int _N) {
  N = _N;
  pts = NULL;
  wts = NULL;
  D = NULL;
  d0 = NULL;

  // Keep track of the assembler object
  assembler = _assembler;
  assembler->incref();

  // Compute the time-spectral part
  tinit = 0.0;
  tfinal = tfinal;
  tfactor = 2.0 / (tfinal - tinit);

  // Initialization of the integration points and weights
  initLGLPointsAndWeights();

  // Initialize the operator
  initOperator();

  // Create the initial conditions
  init = assembler->createVec();
  init->incref();

  // Set the variables
  vars = new TACSSpectralVec(N, assembler);
  vars->incref();
}

TACSSpectralIntegrator::~TACSSpectralIntegrator() {
  assembler->decref();
  init->decref();
  vars->decref();

  if (pts) {
    delete[] pts;
  }
  if (wts) {
    delete[] wts;
  }
  if (D) {
    delete[] D;
  }
  if (d0) {
    delete[] d0;
  }
}

TACSSpectralVec *TACSSpectralIntegrator::createVec() {
  return new TACSSpectralVec(N, assembler);
}

TACSLinearSpectralMat *TACSSpectralIntegrator::createLinearMat() {
  TACSMat *H = assembler->createMat();
  TACSMat *C = assembler->createMat();
  return new TACSLinearSpectralMat(N, D, H, C);
}

void TACSSpectralIntegrator::setInitialConditions(TACSBVec *_init) {
  init->copyValues(_init);
}

void TACSSpectralIntegrator::setVariables(TACSSpectralVec *vec) {
  vars->copyValues(vec);
}

void TACSSpectralIntegrator::assembleRes(TACSSpectralVec *res) {
  TACSBVec *dudt = assembler->createVec();
  dudt->incref();

  for (int i = 0; i < N; i++) {
    const double *Dt = &D[(i + 1) * (N + 1)];

    // Extract the u variables
    TACSBVec *u = vars->getVec(i);

    // Set the values from the initial conditions
    dudt->copyValues(init);
    dudt->scale(Dt[0]);

    for (int j = 0; j < N; j++) {
      dudt->axpy(Dt[j + 1], vars->getVec(j));
    }

    // Set the values of the variables at this point
    assembler->setVariables(u, dudt);

    // Assemble the residual
    assembler->assembleRes(res->getVec(i));
  }

  dudt->decref();
}

void TACSSpectralIntegrator::assembleMat(TACSLinearSpectralMat *mat,
                                         MatrixOrientation matOr) {
  TACSMat *H, *C;
  mat->getMat(&H, &C);
  assembler->assembleJacobian(1.0, 0.0, 0.0, NULL, H, matOr);
  assembler->assembleJacobian(0.0, 1.0, 0.0, NULL, C, matOr);
  mat->setMatrixOrientation(matOr);
}

void TACSSpectralIntegrator::initLGLPointsAndWeights(int max_newton_iters,
                                                     double tol) {
  // Set the initial guess for the points
  if (pts) {
    delete[] pts;
  }
  if (wts) {
    delete[] wts;
  }
  pts = new double[N + 1];
  wts = new double[N + 1];

  // Temporary array to store the
  double *pts0 = new double[N + 1];

  // Set the initial guesses based on the Chebyshev nodes
  for (int k = 0; k < N + 1; k++) {
    wts[k] = 0.0;
    pts[k] = -cos((M_PI * k) / N);
  }

  // Allocate a temporary array to store the polynomials we are going to
  // compute
  double *P = new double[(N + 1) * (N + 1)];

  for (int j = 0; j < N + 1; j++) {
    P[j] = 1.0;
  }

  for (int i = 0; i < max_newton_iters; i++) {
    // Check if we should stop...
    int converged = 1;
    for (int j = 0; j < N + 1; j++) {
      if (fabs(pts[j] - pts0[j]) > tol) {
        converged = 0;
        break;
      }
    }
    if (converged) {
      break;
    }

    // Copy the old values of the points array
    for (int j = 0; j < N + 1; j++) {
      pts0[j] = pts[j];
    }

    for (int j = 0; j < N + 1; j++) {
      P[j] = 1.0;
      P[j + (N + 1)] = pts[j];
    }

    // Use the recursion to solve for the polynomials
    for (int k = 2; k < N + 1; k++) {
      for (int j = 0; j < N + 1; j++) {
        P[j + k * (N + 1)] =
            ((2.0 * k - 1.0) * pts[j] * P[j + (k - 1) * (N + 1)] -
             (k - 1.0) * P[j + (k - 2) * (N + 1)]) /
            (1.0 * k);
      }
    }

    // Now, update the values of points
    for (int j = 0; j < N + 1; j++) {
      pts[j] =
          pts0[j] - (pts[j] * P[j + N * (N + 1)] - P[j + (N - 1) * (N + 1)]) /
                        ((N + 1) * P[j + N * (N + 1)]);
    }
  }

  // Check that the final points are sorted
  for (int j = 0; j < N; j++) {
    if (pts[j + 1] < pts[j]) {
      fprintf(stderr,
              "TACSSpectralIntegrator: LGL points not in ascending order!\n");
    }
  }

  // Compute the weights
  for (int j = 0; j < N + 1; j++) {
    wts[j] = 2.0 / (N * (N + 1) * P[j + N * (N + 1)] * P[j + N * (N + 1)]);
  }

  // Compute the points in time
  tpts = new double[N + 1];
  for (int j = 0; j < N + 1; j++) {
    tpts[j] = tinit + 2.0 * (tfinal - tinit) * (pts[j] + 1.0);
  }

  delete[] pts0;
  delete[] P;
}

/*
  Compute the derivative operator.

  Here D is a matrix that is N x (N + 1).
*/
void TACSSpectralIntegrator::initOperator() {
  // Initialize the operator
  D = new double[N * (N + 1)];

  // Loop over the interpolation points
  for (int p = 0; p < N; p++) {
    double pt = pts[p + 1];
    double *Nx = &D[p * (N + 1)];

    for (int i = 0; i < N + 1; i++) {
      Nx[i] = 0.0;

      // Loop over each point again, except for the current control
      // point, adding the contribution to the shape function
      for (int j = 0; j < N + 1; j++) {
        if (i != j) {
          double d = 1.0 / (pts[i] - pts[j]);

          // Now add up the contribution to the derivative
          for (int k = 0; k < N + 1; k++) {
            if (k != i && k != j) {
              d *= (pt - pts[k]) / (pts[i] - pts[k]);
            }
          }

          // Add the derivative contribution
          Nx[i] += d;
        }
      }
    }
  }

  // Scale the derivative operator
  for (int i = 0; i < N * (N + 1); i++) {
    D[i] *= tfactor;
  }

  // Set the coefficients for the first-order operator
  d0 = new double[N];
  for (int j = 0; j < N; j++) {
    d0[j] = tfactor / (pts[j + 1] - pts[j]);
  }
}