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

/*
  Create a time spectral matrix with the given integrator class
*/
TACSLinearSpectralMat::TACSLinearSpectralMat(TACSSpectralIntegrator *spec) {
  spectral = spec;
  spectral->incref();

  // Set the number of unknown coefficients (exclude ICs)
  N = spectral->getNumLGLNodes() - 1;

  TACSAssembler *assembler = spectral->getAssembler();

  C = assembler->createMat();
  H = assembler->createMat();
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
  spectral->decref();
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
    if (orient == TACS_MAT_NORMAL) {
      for (int i = 0; i < N; i++) {
        int include_ics = 0;  // Do not include the initial conditions
        spectral->computeDeriv(i + 1, x, temp, include_ics);
        C->mult(temp, y->getVec(i));
      }
    } else {
      for (int i = 0; i < N; i++) {
        spectral->computeDerivTranspose(i + 1, x, temp);
        C->mult(temp, y->getVec(i));
      }
    }

    // Add the diagonal contributions
    for (int i = 0; i < N; i++) {
      H->mult(x->getVec(i), temp);
      y->getVec(i)->axpy(1.0, temp);
    }
  }
}

void TACSLinearSpectralMat::getMat(TACSParallelMat **Hmat,
                                   TACSParallelMat **Cmat) {
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

MatrixOrientation TACSLinearSpectralMat::getMatrixOrientation() {
  return orient;
}

int TACSLinearSpectralMat::getFirstOrderCoefficients(const double *d[]) {
  return spectral->getFirstOrderCoefficients(d);
}

/*
  Multigrid preconditioner for a spectral problem in TACS
*/
TACSLinearSpectralMg::TACSLinearSpectralMg(TACSLinearSpectralMat *_mat,
                                           int _nlevels,
                                           TACSAssembler **assembler,
                                           TACSBVecInterp **_interp,
                                           int coarsen_time[]) {
  mat = _mat;
  mat->incref();

  TACSParallelMat *H, *C;
  mat->getMat(&H, &C);

  nlevels = _nlevels;

  int N = 0;

  data = new MgData *[nlevels];
  for (int i = 0; i < nlevels; i++) {
    MgData *fine = NULL;
    TACSBVecInterp *interp = NULL;
    if (i == 0) {
      // Get the first-order coefficients
      const double *d0;
      N = mat->getFirstOrderCoefficients(&d0);

      data[i] = new MgData(fine, N, d0, assembler[i], interp, H, C);
    } else {
      // Check whether to coarsen the time
      if (coarsen_time && coarsen_time[i - 1] && N % 2 == 0) {
        N = N / 2;
      }

      // Pass in NULL for the time coefficients - these are computed internally
      const double *d0 = NULL;

      // Set the finer MgData object and the correct interpolation
      fine = data[i - 1];
      interp = _interp[i - 1];

      // Use a direct method on the coarsest mesh
      int direct = (i == nlevels - 1);
      data[i] =
          new MgData(fine, N, d0, assembler[i], interp, NULL, NULL, direct);
    }
  }
}

TACSLinearSpectralMg::~TACSLinearSpectralMg() {
  mat->decref();

  for (int i = 0; i < nlevels; i++) {
    delete data[i];
  }
  delete[] data;
}

/*
  Factor all the matrix levels
*/
void TACSLinearSpectralMg::factor() {
  for (int i = 0; i < nlevels; i++) {
    data[i]->factor();
    data[i]->setMatrixOrientation(mat->getMatrixOrientation());
  }
}

TACSLinearSpectralMg::TACSComboMat::TACSComboMat(TACSParallelMat *Hmat,
                                                 TacsScalar _alpha,
                                                 TACSParallelMat *Cmat,
                                                 TACSBVec *tmp) {
  alpha = _alpha;
  H = Hmat;
  H->incref();
  C = Cmat;
  C->incref();
  temp = tmp;
  temp->incref();
}

TACSLinearSpectralMg::TACSComboMat::~TACSComboMat() {
  H->decref();
  C->decref();
  temp->decref();
}

TACSVec *TACSLinearSpectralMg::TACSComboMat::createVec() {
  return C->createVec();
}

void TACSLinearSpectralMg::TACSComboMat::mult(TACSVec *xvec, TACSVec *yvec) {
  TACSBVec *x = dynamic_cast<TACSBVec *>(xvec);
  TACSBVec *y = dynamic_cast<TACSBVec *>(yvec);

  if (x && y) {
    C->mult(x, temp);
    H->mult(x, y);
    y->axpy(alpha, temp);
  }
}

/*
  Create the data at the specified mesh level
*/
TACSLinearSpectralMg::MgData::MgData(MgData *_fine, int Nval, const double d[],
                                     TACSAssembler *_assembler,
                                     TACSBVecInterp *_interp,
                                     TACSParallelMat *Hmat,
                                     TACSParallelMat *Cmat, int _direct) {
  fine = _fine;
  N = Nval;
  assembler = _assembler;
  assembler->incref();
  interp = _interp;
  if (interp) {
    interp->incref();
  }

  // Set the default matrix orientations
  orient = TACS_MAT_NORMAL;

  // Compute the non-zero pattern on this level
  if (fine) {
    interp->computeGalerkinNonZeroPattern(fine->H, &H);
    C = dynamic_cast<TACSParallelMat *>(H->createDuplicate());
  } else {
    H = Hmat;
    C = Cmat;
  }

  H->incref();
  C->incref();

  if (fine) {
    const double *d = fine->d0;  // Time coefficients on the finer mesh

    if (N == fine->N) {
      // Copy the coefficients from the finer mesh
      d0 = new double[N];
      for (int i = 0; i < N; i++) {
        d0[i] = d[i];
      }

      w0 = w1 = NULL;
    } else if (N == fine->N / 2) {
      d0 = new double[N];
      for (int i = 0; i < N; i++) {
        d0[i] = (d[2 * i] * d[2 * i + 1]) / (d[2 * i] + d[2 * i + 1]);
      }

      // Compute the interpolation weights between time levels
      w0 = new double[N];
      w1 = new double[N];
      for (int i = 0; i < N; i++) {
        w0[i] = d[2 * i] / (d[2 * i] + d[2 * i + 1]);
        w1[i] = d[2 * i + 1] / (d[2 * i] + d[2 * i + 1]);
      }
    }
  } else {
    // Copy the input coefficients
    d0 = new double[N];
    for (int i = 0; i < N; i++) {
      d0[i] = d[i];
    }

    // Set the interpolation weights to NULL
    w0 = w1 = NULL;
  }

  // Only allocate the solution and right-hand-side vectors on the coarse meshes
  if (fine) {
    x = new TACSSpectralVec(N, assembler);
    x->incref();
    b = new TACSSpectralVec(N, assembler);
    b->incref();
  } else {
    x = b = NULL;
  }
  r = new TACSSpectralVec(N, assembler);
  r->incref();

  // Create a temporary vector on this level
  temp = assembler->createVec();
  temp->incref();

  direct = _direct;
  mats = NULL;
  smoothers = NULL;
  direct_mats = NULL;
  direct_pcs = NULL;

  if (direct) {
    direct_mats = new TACSParallelMat *[N];
    direct_pcs = new TACSBlockCyclicPc *[N];

    for (int i = 0; i < N; i++) {
      direct_mats[i] = dynamic_cast<TACSParallelMat *>(H->createDuplicate());
      direct_mats[i]->incref();
      direct_pcs[i] = new TACSBlockCyclicPc(direct_mats[i]);
      direct_pcs[i]->incref();
    }
  } else {  // Allocate the smoothers
    mats_temp = assembler->createVec();
    mats_temp->incref();

    mats = new TACSComboMat *[N];
    smoothers = new TACSChebyshevSmoother *[N];

    int chebyshev_degree = 4;
    for (int i = 0; i < N; i++) {
      mats[i] = new TACSComboMat(H, d0[i], C, mats_temp);
      mats[i]->incref();

      smoothers[i] = new TACSChebyshevSmoother(mats[i], chebyshev_degree);
      smoothers[i]->incref();
    }
  }
}

TACSLinearSpectralMg::MgData::~MgData() {
  assembler->decref();
  if (interp) {
    interp->decref();
  }
  if (x) {
    x->decref();
  }
  if (b) {
    b->decref();
  }
  r->decref();
  temp->decref();
  H->decref();
  C->decref();

  if (w0 && w1) {
    delete[] w0;
    delete[] w1;
  }

  if (direct) {
    for (int i = 0; i < N; i++) {
      direct_mats[i]->decref();
      direct_pcs[i]->decref();
    }
    delete[] direct_mats;
    delete[] direct_pcs;
  } else {
    mats_temp->decref();
    for (int i = 0; i < N; i++) {
      mats[i]->decref();
      smoothers[i]->decref();
    }
    delete[] mats;
    delete[] smoothers;
  }
}

void TACSLinearSpectralMg::MgData::setMatrixOrientation(
    MatrixOrientation matOr) {
  orient = matOr;
}

/*
  Factor the matrix at the current level
*/
void TACSLinearSpectralMg::MgData::factor() {
  if (fine) {
    interp->computeGalerkin(fine->H, H);
    interp->computeGalerkin(fine->C, C);
  }

  if (direct) {
    for (int i = 0; i < N; i++) {
      direct_mats[i]->copyValues(H);
      direct_mats[i]->axpy(d0[i], C);
      direct_pcs[i]->factor();
    }
  } else {
    for (int i = 0; i < N; i++) {
      smoothers[i]->factor();
    }
  }
}

/*
  Apply the smoother on the current mesh level

  if orient == TACS_MAT_NORMAL:
  [ (H + d0[0] * C) |                                   ][out[0]] = [in[0]]
  [     - d0[1] * C | (H + d0[1] * C) |                 ][out[1]] = [in[1]]
  [                       - d0[2] * C | (H + d0[2] * C) ][out[2]] = [in[2]]

  else:
  [ (H + d0[0] * C) |    - d0[1] * C                    ][out[0]] = [in[0]]
  [                 | (H + d0[1] * C) |    - d0[2] * C  ][out[1]] = [in[1]]
  [                                   | (H + d0[2] * C) ][out[2]] = [in[2]]
*/
void TACSLinearSpectralMg::MgData::applyFactor(TACSSpectralVec *in,
                                               TACSSpectralVec *out) {
  if (orient == TACS_MAT_NORMAL) {
    if (direct) {
      direct_pcs[0]->applyFactor(in->getVec(0), out->getVec(0));

      for (int i = 1; i < N; i++) {
        // Form the right-hand-side = in[i] + d0[i] * C * out[i - 1]
        C->mult(out->getVec(i - 1), temp);
        temp->axpby(1.0, d0[i], in->getVec(i));
        direct_pcs[i]->applyFactor(temp, out->getVec(i));
      }
    } else {
      smoothers[0]->applyFactor(in->getVec(0), out->getVec(0));

      for (int i = 1; i < N; i++) {
        // Form the right-hand-side = in[i] + d0[i] * C * out[i - 1]
        C->mult(out->getVec(i - 1), temp);
        temp->axpby(1.0, d0[i], in->getVec(i));
        smoothers[i]->applyFactor(temp, out->getVec(i));
      }
    }
  } else {
    if (direct) {
      direct_pcs[N - 1]->applyFactor(in->getVec(N - 1), out->getVec(N - 1));

      for (int i = N - 2; i >= 0; i--) {
        // Form the right-hand-side = in[i] + d0[i + 1] * C * out[i + 1]
        C->mult(out->getVec(i + 1), temp);
        temp->axpby(1.0, d0[i + 1], in->getVec(i));
        direct_pcs[i]->applyFactor(temp, out->getVec(i));
      }
    } else {
      smoothers[N - 1]->applyFactor(in->getVec(N - 1), out->getVec(N - 1));

      for (int i = N - 2; i >= 0; i--) {
        // Form the right-hand-side = in[i] + d0[i + 1] * C * out[i + 1]
        C->mult(out->getVec(i + 1), temp);
        temp->axpby(1.0, d0[i + 1], in->getVec(i));
        smoothers[i]->applyFactor(temp, out->getVec(i));
      }
    }
  }
}

/*
  Perform a matrix-vector multiplication

  if orient == TACS_MAT_NORMAL:
  [ (H + d0[0] * C) |                                   ][in[0]] = [out[0]]
  [     - d0[1] * C | (H + d0[1] * C) |                 ][in[1]] = [out[1]]
  [                       - d0[2] * C | (H + d0[2] * C) ][in[2]] = [out[2]]

  else:
  [ (H + d0[0] * C) |    - d0[1] * C                    ][in[0]] = [out[0]]
  [                 | (H + d0[1] * C) |    - d0[2] * C  ][in[1]] = [out[1]]
  [                                   | (H + d0[2] * C) ][in[2]] = [out[2]]
*/
void TACSLinearSpectralMg::MgData::mult(TACSSpectralVec *in,
                                        TACSSpectralVec *out) {
  if (orient == TACS_MAT_NORMAL) {
    for (int i = 0; i < N; i++) {
      TACSBVec *outvec = out->getVec(i);
      TACSBVec *invec = in->getVec(i);

      H->mult(invec, outvec);
      if (i > 0) {
        outvec->axpy(-d0[i], temp);
      }

      C->mult(invec, temp);
      outvec->axpy(d0[i], temp);
    }
  } else {
    for (int i = N - 1; i >= 0; i--) {
      TACSBVec *outvec = out->getVec(i);
      TACSBVec *invec = in->getVec(i);

      H->mult(invec, outvec);
      if (i < N - 1) {
        outvec->axpy(-d0[i + 1], temp);
      }

      C->mult(invec, temp);
      outvec->axpy(d0[i], temp);
    }
  }
}

/*
  Apply a restriction from the input vector to the output vector
*/
void TACSLinearSpectralMg::MgData::restriction(TACSSpectralVec *in,
                                               TACSSpectralVec *out) {
  if (interp) {
    if (w0 && w1) {
      // Compute the restriction using the transpose of the weights
      // [init = 0] = [ 1.0    w0[0]                     ][init = 0]
      // [out[0]  ] = [        w1[0]  1.0    w0[1]       ][ in[0]  ]
      // [out[1]  ] = [                      w1[1]   1.0 ][ in[1]  ]
      //                                                  [ in[2]  ]
      //                                                  [ in[3]  ]

      for (int i = 0; i < N; i++) {
        TACSBVec *t = fine->temp;

        // Compute the scaling factor to normalize the weights
        double scale = w1[i] + 1.0;
        if (i < N - 1) {
          scale += w0[i + 1];
        }
        scale = 1.0 / scale;

        // Compute the weighted input vector
        t->copyValues(in->getVec(2 * i));
        t->scale(scale * w1[i]);

        // Add the contribution from the middle time plane
        t->axpy(scale, in->getVec(2 * i + 1));

        // Add the contribution from the final time plane
        if (i < N - 1) {
          t->axpy(scale * w0[i + 1], in->getVec(2 * (i + 1)));
        }

        TACSBVec *outvec = out->getVec(i);
        interp->multTranspose(t, outvec);

        // Apply boundary conditions at this new level
        outvec->applyBCs(assembler->getBcMap());
      }

    } else {
      // Use a straight injection in time
      for (int i = 0; i < N; i++) {
        TACSBVec *invec = in->getVec(i);
        TACSBVec *outvec = out->getVec(i);
        interp->multTranspose(invec, outvec);

        // Apply boundary conditions at this new level
        outvec->applyBCs(assembler->getBcMap());
      }
    }
  }
}

void TACSLinearSpectralMg::MgData::interpolateAdd(TACSSpectralVec *in,
                                                  TACSSpectralVec *out) {
  if (interp) {
    if (w0 && w1) {
      // Interpolate from one mesh to the next
      // [init = 0] = [ 1.0                  ][init = 0 ]
      // [ out[0] ] = [ w0[0]  w1[0]         ][ in[0]   ]
      // [ out[1] ] = [        1.0           ][ in[1]   ]
      // [ out[2] ] = [        w0[1]  w1[1]  ]
      // [ out[3] ] = [               1.0    ]

      for (int i = 0; i < N; i++) {
        if (i > 0) {
          // Do the direct interpolation
          TACSBVec *invec = in->getVec(i);
          TACSBVec *outvec = out->getVec(2 * i);
          interp->multAdd(invec, outvec, outvec);

          // Apply boundary conditions on the fine mesh
          if (fine) {
            outvec->applyBCs(fine->assembler->getBcMap());
          }

          // Copy the input to the temporary array and scale
          temp->copyValues(in->getVec(i - 1));
          temp->scale(w0[i]);
        } else {
          temp->zeroEntries();
        }

        TACSBVec *outvec = out->getVec(2 * i + 1);
        TACSBVec *invec = in->getVec(i);
        temp->axpy(w1[i], invec);
        interp->multAdd(temp, outvec, outvec);

        // Apply boundary conditions on the fine mesh
        if (fine) {
          outvec->applyBCs(fine->assembler->getBcMap());
        }
      }
    } else {
      // For now, use a straight injection in time
      for (int i = 0; i < N; i++) {
        TACSBVec *invec = in->getVec(i);
        TACSBVec *outvec = out->getVec(i);
        interp->multAdd(invec, outvec, outvec);

        // Apply boundary conditions on the fine mesh
        if (fine) {
          outvec->applyBCs(fine->assembler->getBcMap());
        }
      }
    }
  }
}

void TACSLinearSpectralMg::applyFactor(TACSVec *bvec, TACSVec *xvec) {
  // Set the RHS at the finest level
  data[0]->b = dynamic_cast<TACSSpectralVec *>(bvec);
  data[0]->x = dynamic_cast<TACSSpectralVec *>(xvec);

  if (data[0]->b && data[0]->x) {
    data[0]->x->zeroEntries();
    applyMg(0);
  } else {
    fprintf(stderr,
            "TACSLinearSpectralMg type error: Input/output must be "
            "TACSSpectralVec\n");
  }

  data[0]->b = NULL;
  data[0]->x = NULL;
}

void TACSLinearSpectralMg::applyMg(int level) {
  // If we've made it to the lowest level, apply the direct solver
  // otherwise, perform multigrid on the next-lowest level
  if (level == nlevels - 1) {
    data[level]->applyFactor(data[level]->b, data[level]->x);
    return;
  }

  // Perform iters[level] cycle at the next lowest level
  // for (int k = 0; k < iters[level]; k++) {
  // Pre-smooth at the current level
  data[level]->applyFactor(data[level]->b, data[level]->x);

  // Compute r[level] = b[level] - A*x[level]
  data[level]->mult(data[level]->x, data[level]->r);
  data[level]->r->axpby(1.0, -1.0, data[level]->b);

  // Restrict the residual to the next lowest level
  // to form the RHS at that level
  data[level + 1]->restriction(data[level]->r, data[level + 1]->b);
  data[level + 1]->x->zeroEntries();

  applyMg(level + 1);

  // Interpolate back from the next lowest level
  data[level + 1]->interpolateAdd(data[level + 1]->x, data[level]->x);
  // }

  // Post-Smooth the residual
  data[level]->applyFactor(data[level]->b, data[level]->x);
}

TACSSpectralIntegrator::TACSSpectralIntegrator(TACSAssembler *_assembler,
                                               double _tfinal, int _N) {
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
  tfinal = _tfinal;
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

int TACSSpectralIntegrator::getNumLGLNodes() { return N + 1; }

double TACSSpectralIntegrator::getPointAtLGLNode(int index) {
  if (index >= 0.0 && index <= N) {
    return pts[index];
  }
  return 0.0;
}

double TACSSpectralIntegrator::getTimeAtLGLNode(int index) {
  if (index >= 0.0 && index <= N) {
    return tpts[index];
  }
  return 0.0;
}

double TACSSpectralIntegrator::getWeightAtLGLNode(int index) {
  if (index >= 0.0 && index <= N) {
    return wts[index];
  }
  return 0.0;
}

TACSAssembler *TACSSpectralIntegrator::getAssembler() { return assembler; }

int TACSSpectralIntegrator::getFirstOrderCoefficients(const double *d[]) {
  if (d) {
    *d = d0;
  }
  return N;
}

TACSSpectralVec *TACSSpectralIntegrator::createVec() {
  return new TACSSpectralVec(N, assembler);
}

TACSLinearSpectralMat *TACSSpectralIntegrator::createLinearMat() {
  return new TACSLinearSpectralMat(this);
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
    // Extract the u variable for the i+1 LGL point
    TACSBVec *u = vars->getVec(i);

    // Compute the derivative at the time interval
    computeDeriv(i + 1, vars, dudt);

    // Set the values of the variables at this point
    assembler->setVariables(u, dudt);

    // Assemble the residual
    assembler->assembleRes(res->getVec(i));
  }

  dudt->decref();
}

void TACSSpectralIntegrator::assembleMat(TACSLinearSpectralMat *mat,
                                         MatrixOrientation matOr) {
  TACSParallelMat *H, *C;
  mat->getMat(&H, &C);
  assembler->assembleJacobian(1.0, 0.0, 0.0, NULL, H, matOr);
  assembler->assembleJacobian(0.0, 1.0, 0.0, NULL, C, matOr);
  mat->setMatrixOrientation(matOr);
}

/*
  Use the spectral operator, compute the derivative
*/
void TACSSpectralIntegrator::computeDeriv(int index, TACSSpectralVec *sol,
                                          TACSBVec *dudt,
                                          int include_init_conditions) {
  if (index >= 0 && index < N + 1) {
    const double *Dt = &D[index * (N + 1)];

    // Include the initial conditions
    if (include_init_conditions) {
      dudt->copyValues(init);
      dudt->scale(Dt[0]);
    } else {
      dudt->zeroEntries();
    }

    // Finish adding the contributions from the remainder of the time interval
    for (int j = 0; j < N; j++) {
      dudt->axpy(Dt[j + 1], sol->getVec(j));
    }
  }
}

/*
  Compute the action of the transpose of the spectral derivative operator
*/
void TACSSpectralIntegrator::computeDerivTranspose(int index,
                                                   TACSSpectralVec *sol,
                                                   TACSBVec *dudt) {
  if (index > 0 && index < N + 1) {
    const double *Dt = &D[index + (N + 1)];

    dudt->zeroEntries();
    for (int j = 0; j < N; j++) {
      dudt->axpy(Dt[0], sol->getVec(j));
      Dt += (N + 1);
    }
  }
}

void TACSSpectralIntegrator::computeSolutionAndDeriv(double time,
                                                     TACSSpectralVec *sol,
                                                     TACSBVec *u,
                                                     TACSBVec *dudt) {
  if (!sol) {
    sol = vars;
  }

  if (time >= tinit && time <= tfinal) {
    double *P = new double[N + 1];
    double *Px = new double[N + 1];

    // Compute the parametric point in the interval
    double pt = -1.0 + 2.0 * (time - tinit) / (tfinal - tinit);

    evalInterpolation(pt, P, Px);

    if (u) {
      // Include the initial conditions
      u->copyValues(init);
      u->scale(P[0]);

      for (int j = 0; j < N; j++) {
        u->axpy(P[j + 1], sol->getVec(j));
      }
    }
    if (dudt) {
      // Include the initial conditions
      dudt->copyValues(init);
      dudt->scale(tfactor * Px[0]);

      for (int j = 0; j < N; j++) {
        dudt->axpy(tfactor * Px[j + 1], sol->getVec(j));
      }
    }

    delete[] P;
    delete[] Px;
  }
}

void TACSSpectralIntegrator::evalFunctions(int num_funcs, TACSFunction **funcs,
                                           TacsScalar *fvals) {
  // TODO: integrate functions as the forward problem is integrated
  // Check whether these are two-stage or single-stage functions
  int twoStage = 0;
  for (int n = 0; n < num_funcs; n++) {
    if (funcs[n] && funcs[n]->getStageType() == TACSFunction::TWO_STAGE) {
      twoStage = 1;
      break;
    }
  }

  TACSBVec *dudt = assembler->createVec();
  dudt->incref();

  // Initialize the function if had already not been initialized
  if (twoStage) {
    // First stage
    for (int n = 0; n < num_funcs; n++) {
      if (funcs[n]) {
        funcs[n]->initEvaluation(TACSFunction::INITIALIZE);
      }
    }

    for (int i = 0; i < N + 1; i++) {
      // Get the solution values at the i-th LGL node
      TACSBVec *u = NULL;
      if (i == 0) {
        u = init;
      } else {
        u = vars->getVec(i - 1);
      }
      computeDeriv(i, vars, dudt);

      // Set the simulation time and variables
      assembler->setSimulationTime(tpts[i]);
      assembler->setVariables(u, dudt);

      // Integrate the function
      TacsScalar tcoeff = 0.5 * wts[i] * (tfinal - tinit);
      assembler->integrateFunctions(tcoeff, TACSFunction::INITIALIZE, num_funcs,
                                    funcs);
    }

    for (int n = 0; n < num_funcs; n++) {
      if (funcs[n]) {
        funcs[n]->finalEvaluation(TACSFunction::INITIALIZE);
      }
    }
  }

  // Second stage
  for (int n = 0; n < num_funcs; n++) {
    if (funcs[n]) {
      funcs[n]->initEvaluation(TACSFunction::INTEGRATE);
    }
  }

  for (int i = 0; i < N + 1; i++) {
    TACSBVec *u = NULL;
    if (i == 0) {
      u = init;
    } else {
      u = vars->getVec(i - 1);
    }
    computeDeriv(i, vars, dudt);

    // Set the simulation time and variables
    assembler->setSimulationTime(tpts[i]);
    assembler->setVariables(u, dudt);

    // Integrate the function
    TacsScalar tcoeff = 0.5 * wts[i] * (tfinal - tinit);
    assembler->integrateFunctions(tcoeff, TACSFunction::INTEGRATE, num_funcs,
                                  funcs);
  }

  for (int n = 0; n < num_funcs; n++) {
    if (funcs[n]) {
      funcs[n]->finalEvaluation(TACSFunction::INTEGRATE);
    }
  }

  // Retrieve the function values
  for (int n = 0; n < num_funcs; n++) {
    fvals[n] = 0.0;
    if (funcs[n]) {
      fvals[n] = funcs[n]->getFunctionValue();
    }
  }

  dudt->decref();
}

void TACSSpectralIntegrator::evalSVSens(TACSFunction *func,
                                        TACSSpectralVec *dfdu) {
  TACSBVec *dudt = assembler->createVec();
  dudt->incref();

  for (int i = 1; i < N + 1; i++) {
    // Get the solution values at the i-th LGL node
    TACSBVec *u = NULL;
    if (i == 0) {
      u = init;
    } else {
      u = vars->getVec(i - 1);
    }
    computeDeriv(i, vars, dudt);

    // Set the simulation time and variables
    assembler->setSimulationTime(tpts[i]);
    assembler->setVariables(u, dudt);

    TacsScalar tcoeff = 0.5 * wts[i] * (tfinal - tinit);
    TACSBVec *vec = dfdu->getVec(i - 1);
    assembler->addSVSens(tcoeff, 0.0, 0.0, 1, &func, &vec);
  }

  dudt->decref();
}

void TACSSpectralIntegrator::addDVSens(TACSFunction *func, TACSBVec *dfdx) {
  TACSBVec *dudt = assembler->createVec();
  dudt->incref();

  for (int i = 0; i < N + 1; i++) {
    // Get the solution values at the i-th LGL node
    TACSBVec *u = NULL;
    if (i == 0) {
      u = init;
    } else {
      u = vars->getVec(i - 1);
    }
    computeDeriv(i, vars, dudt);

    // Set the simulation time and variables
    assembler->setSimulationTime(tpts[i]);
    assembler->setVariables(u, dudt);

    TacsScalar tcoeff = 0.5 * wts[i] * (tfinal - tinit);
    assembler->addDVSens(tcoeff, 1, &func, &dfdx);
  }

  dudt->decref();
}

void TACSSpectralIntegrator::addAdjointResProduct(TacsScalar scale,
                                                  TACSSpectralVec *adjoint,
                                                  TACSBVec *dfdx) {
  TACSBVec *dudt = assembler->createVec();
  dudt->incref();

  for (int i = 1; i < N + 1; i++) {
    // Get the solution values at the i-th LGL node
    TACSBVec *u = NULL;
    if (i == 0) {
      u = init;
    } else {
      u = vars->getVec(i - 1);
    }
    computeDeriv(i, vars, dudt);

    // Set the simulation time and variables
    assembler->setSimulationTime(tpts[i]);
    assembler->setVariables(u, dudt);

    TACSBVec *adj = adjoint->getVec(i - 1);
    assembler->addAdjointResProducts(scale, 1, &adj, &dfdx);
  }

  dudt->decref();
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
  for (int i = 0; i < N + 1; i++) {
    pts0[i] = 0.0;
  }

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
    tpts[j] = tinit + 0.5 * (tfinal - tinit) * (pts[j] + 1.0);
  }

  delete[] pts0;
  delete[] P;
}

/*
  Compute the derivative operator.

  Here D is a matrix that is (N + 1) x (N + 1).
*/
void TACSSpectralIntegrator::initOperator() {
  // Initialize the operator
  D = new double[(N + 1) * (N + 1)];
  double *P = new double[N + 1];

  for (int p = 0; p < N + 1; p++) {
    double *Px = &D[p * (N + 1)];
    double pt = pts[p];
    evalInterpolation(pt, P, Px);
  }

  delete[] P;

  // Scale the derivative operator
  for (int i = 0; i < (N + 1) * (N + 1); i++) {
    D[i] *= tfactor;
  }

  // Set the coefficients for the first-order operator
  d0 = new double[N];
  for (int j = 0; j < N; j++) {
    d0[j] = tfactor / (pts[j + 1] - pts[j]);
  }
}

void TACSSpectralIntegrator::evalInterpolation(double pt, double P[],
                                               double Px[]) {
  for (int i = 0; i < N + 1; i++) {
    P[i] = 1.0;
    Px[i] = 0.0;

    // Loop over each point again, except for the current control
    // point, adding the contribution to the shape function
    for (int j = 0; j < N + 1; j++) {
      if (i != j) {
        double d = 1.0 / (pts[i] - pts[j]);
        P[i] *= (pt - pts[j]) * d;

        // Now add up the contribution to the derivative
        for (int k = 0; k < N + 1; k++) {
          if (k != i && k != j) {
            d *= (pt - pts[k]) / (pts[i] - pts[k]);
          }
        }

        // Add the derivative contribution
        Px[i] += d;
      }
    }
  }
}
