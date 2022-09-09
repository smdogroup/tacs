/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "JacobiDavidson.h"

#include "TacsUtilities.h"
#include "tacslapack.h"

TACSJDSimpleOperator::TACSJDSimpleOperator(TACSAssembler *_assembler,
                                           TACSMat *_mat, TACSPc *_pc) {
  assembler = _assembler;
  assembler->incref();
  mat = _mat;
  mat->incref();
  pc = _pc;
  pc->incref();
}

TACSJDSimpleOperator::~TACSJDSimpleOperator() {
  assembler->decref();
  mat->decref();
  pc->decref();
}

// Get the MPI_Comm
MPI_Comm TACSJDSimpleOperator::getMPIComm() { return assembler->getMPIComm(); }

// Create a vector
TACSVec *TACSJDSimpleOperator::createVec() { return assembler->createVec(); }

// Set the eigenvalue estimate (and reset the factorization)
// pc = K in this case
void TACSJDSimpleOperator::setEigenvalueEstimate(double estimate) {
  pc->factor();
}

// Apply the preconditioner
void TACSJDSimpleOperator::applyFactor(TACSVec *x, TACSVec *y) {
  pc->applyFactor(x, y);
  assembler->applyBCs(y);
}

// Apply boundary conditions associated with the matrix
void TACSJDSimpleOperator::applyBCs(TACSVec *x) { assembler->applyBCs(x); }

// Perform a dot-product of two vectors in the operator space e.g.
// output <- x^{T}*y
TacsScalar TACSJDSimpleOperator::dot(TACSVec *x, TACSVec *y) {
  return y->dot(x);
}

// Matrix-vector products with either the A or B operators.
void TACSJDSimpleOperator::multA(TACSVec *x, TACSVec *y) {
  mat->mult(x, y);
  assembler->applyBCs(y);
}

void TACSJDSimpleOperator::multB(TACSVec *x, TACSVec *y) {
  y->copyValues(x);
  assembler->applyBCs(y);
}

/*
  JD operator class for frequency analysis
*/
TACSJDFrequencyOperator::TACSJDFrequencyOperator(TACSAssembler *_assembler,
                                                 TACSMat *_kmat, TACSMat *_mmat,
                                                 TACSMat *_pc_mat,
                                                 TACSPc *_pc) {
  assembler = _assembler;
  assembler->incref();
  kmat = _kmat;
  kmat->incref();
  mmat = _mmat;
  mmat->incref();
  pc_mat = _pc_mat;
  pc_mat->incref();
  pc = _pc;
  pc->incref();
  work = assembler->createVec();
  work->incref();
}

TACSJDFrequencyOperator::~TACSJDFrequencyOperator() {
  assembler->decref();
  kmat->decref();
  mmat->decref();
  pc_mat->decref();
  pc->decref();
  work->decref();
}

// Get the MPI_Comm
MPI_Comm TACSJDFrequencyOperator::getMPIComm() {
  return assembler->getMPIComm();
}

// Create a vector
TACSVec *TACSJDFrequencyOperator::createVec() { return assembler->createVec(); }

// Set the eigenvalue estimate (and reset the factorization)
// pc = K in this case
void TACSJDFrequencyOperator::setEigenvalueEstimate(double estimate) {
  pc_mat->zeroEntries();
  pc_mat->copyValues(kmat);
  if (estimate != 0.0) {
    pc_mat->axpy(-estimate, mmat);
  }
  assembler->applyBCs(pc_mat);
  pc->factor();
}

// Apply the preconditioner
void TACSJDFrequencyOperator::applyFactor(TACSVec *x, TACSVec *y) {
  pc->applyFactor(x, y);
  assembler->applyBCs(y);
}

// Apply boundary conditions associated with the matrix
void TACSJDFrequencyOperator::applyBCs(TACSVec *x) { assembler->applyBCs(x); }

// Perform a dot-product of two vectors in the operator space e.g.
// output <- x^{T}*B*y
TacsScalar TACSJDFrequencyOperator::dot(TACSVec *x, TACSVec *y) {
  mmat->mult(x, work);
  return y->dot(work);
}

// Matrix-vector products with either the A or B operators.
void TACSJDFrequencyOperator::multA(TACSVec *x, TACSVec *y) {
  kmat->mult(x, y);
  assembler->applyBCs(y);
}

void TACSJDFrequencyOperator::multB(TACSVec *x, TACSVec *y) {
  mmat->mult(x, y);
  assembler->applyBCs(y);
}

/*
  Create the Jacobi-Davidson solver

  This allocates all the vectors needed for the Jacobi-Davidson method
  as well as extra data needed for the FGMRES sub-iteraitons and other
  data required for the small LAPACK eigenproblems.

  input:
  oper:              the Jacobi-Davidson operator
  max_eigen_vectors: the number of eigenvectors/values sought
  max_jd_size:       maximum size of the Jacobi-Davidson search space
  max_gmres_size:    maximum number of FGMRES iterations
*/
TACSJacobiDavidson::TACSJacobiDavidson(TACSJacobiDavidsonOperator *_oper,
                                       int _max_eigen_vectors, int _max_jd_size,
                                       int _max_gmres_size) {
  oper = _oper;
  oper->incref();

  // The maximum size of the deflation space
  max_eigen_vectors = _max_eigen_vectors;

  // The maximum size of the Jacobi-Davidson subspace
  max_jd_size = _max_jd_size;

  // The maximum number of gmres iterations
  max_gmres_size = _max_gmres_size;

  // The number of converged eigenvectors
  nconverged = 0;

  // Set the number of vectors to recycle (0 by default)
  num_recycle_vecs = 0;
  recycle_type = JD_NUM_RECYCLE;

  // The eigen tolerance
  eig_rtol = 1e-9;
  eig_atol = 1e-30;
  theta_cutoff = 0.1;

  // The residual tolerances for the GMRES iterations
  rtol = 1e-9;
  atol = 1e-30;

  // Store the matrix
  M = new TacsScalar[(max_jd_size + 1) * (max_jd_size + 1)];
  ritzvecs = new double[(max_jd_size + 1) * (max_jd_size + 1)];
  ritzvals = new double[(max_jd_size + 1)];

  // Allocate space for the vectors
  V = new TACSVec *[max_jd_size + 1];
  Q = new TACSVec *[max_eigen_vectors + 1];
  P = new TACSVec *[max_eigen_vectors + 1];
  eigvals = new TacsScalar[max_eigen_vectors];
  eigerror = new TacsScalar[max_eigen_vectors];
  eigindex = new int[max_eigen_vectors];

  // Allocate the variables
  for (int i = 0; i < max_jd_size + 1; i++) {
    V[i] = oper->createVec();
    V[i]->incref();
  }
  for (int i = 0; i < max_eigen_vectors + 1; i++) {
    Q[i] = oper->createVec();
    P[i] = oper->createVec();
    Q[i]->incref();
    P[i]->incref();
  }

  // Create and populate the pointer into the H columns
  Hptr = new int[max_gmres_size + 1];
  Hptr[0] = 0;
  for (int i = 0; i < max_gmres_size; i++) {
    Hptr[i + 1] = Hptr[i] + i + 2;
  }

  // Allocate the Hessenberg matrix
  H = new TacsScalar[Hptr[max_gmres_size]];
  res = new TacsScalar[max_gmres_size + 1];

  // Allocate space for the unitary matrix Q
  Qsin = new TacsScalar[max_gmres_size];
  Qcos = new TacsScalar[max_gmres_size];

  // Allocate space for the vectors
  W = new TACSVec *[max_gmres_size + 1];
  for (int i = 0; i < max_gmres_size + 1; i++) {
    W[i] = oper->createVec();
    W[i]->incref();
  }

  // Create the work vector
  work = oper->createVec();
  work->incref();

  Z = new TACSVec *[max_gmres_size];
  for (int i = 0; i < max_gmres_size; i++) {
    Z[i] = oper->createVec();
    Z[i]->incref();
  }
}

/*
  Free the data associated with the object
*/
TACSJacobiDavidson::~TACSJacobiDavidson() {
  // Free the data to solve the GMRES problem
  delete[] H;
  delete[] Hptr;
  delete[] res;
  delete[] Qsin;
  delete[] Qcos;

  // Free the Jacobi-Davidson vectors
  for (int i = 0; i < max_jd_size + 1; i++) {
    V[i]->decref();
  }
  delete[] V;
  for (int i = 0; i < max_eigen_vectors + 1; i++) {
    Q[i]->decref();
    P[i]->decref();
  }
  delete[] Q;
  delete[] P;

  // Free the JD matrix
  delete[] M;
  delete[] ritzvecs;
  delete[] ritzvals;
  delete[] eigvals;
  delete[] eigerror;
  delete[] eigindex;

  // Free the GMRES subspace vectors
  W[0]->decref();
  for (int i = 0; i < max_gmres_size; i++) {
    W[i + 1]->decref();
  }
  delete[] W;

  for (int i = 0; i < max_gmres_size; i++) {
    Z[i]->decref();
  }
  delete[] Z;

  work->decref();
}

// Get the MPI_Comm
MPI_Comm TACSJacobiDavidson::getMPIComm() { return oper->getMPIComm(); }

/*
  Get the number of converged eigenvalues
*/
int TACSJacobiDavidson::getNumConvergedEigenvalues() { return nconverged; }

/*!
  Extract the eigenvalue from the analysis
*/
TacsScalar TACSJacobiDavidson::extractEigenvalue(int n, TacsScalar *error) {
  if (error) {
    *error = 0.0;
  }

  if (n >= 0 && n < nconverged) {
    int index = eigindex[n];
    if (error) {
      oper->multA(Q[index], work);
      work->axpy(-eigvals[index], P[index]);
      *error = work->norm();
    }

    return eigvals[index];
  } else if (n >= 0 && n < max_eigen_vectors) {
    int index = n - nconverged;

    // Get the Ritz value
    TacsScalar theta = ritzvals[index];

    // Assemble the predicted Ritz vector
    if (error) {
      // Set the temporary vector space
      TACSVec *qvec = Q[n];
      TACSVec *pvec = P[n];
      qvec->zeroEntries();
      for (int j = 0; j <= max_jd_size; j++) {
        qvec->axpy(ritzvecs[n * (max_jd_size + 1) + j], V[j]);
      }
      oper->applyBCs(qvec);

      // Normalize the predicted eigenvalue
      TacsScalar qnorm = sqrt(oper->dot(qvec, qvec));
      qvec->scale(1.0 / qnorm);

      // Compute the residual: work = A*q - theta*B*q
      // and store it in the work vector
      oper->multA(qvec, work);
      oper->multB(qvec, pvec);
      work->axpy(-theta, pvec);
      oper->applyBCs(work);

      // Set the error
      *error = work->norm();
    }

    return theta;
  }

  return 0.0;
}

/*!
  Extract the eigenvector and eigenvalue from the eigenvalue analysis
*/
TacsScalar TACSJacobiDavidson::extractEigenvector(int n, TACSVec *ans,
                                                  TacsScalar *error) {
  if (error) {
    *error = 0.0;
  }

  if (n >= 0 && n < nconverged) {
    int index = eigindex[n];
    if (error) {
      oper->multA(Q[index], work);
      work->axpy(-eigvals[index], P[index]);
      *error = work->norm();
    }
    if (ans) {
      ans->copyValues(Q[index]);
    }

    return eigvals[index];
  } else if (n >= 0 && n < max_eigen_vectors) {
    int index = n - nconverged;

    // Get the Ritz value
    TacsScalar theta = ritzvals[index];

    // Assemble the predicted Ritz vector
    if (error || ans) {
      // Set which vectors to use
      TACSVec *qvec = Q[n];
      if (ans) {
        qvec = ans;
      }
      TACSVec *pvec = P[n];

      qvec->zeroEntries();
      for (int j = 0; j <= max_jd_size; j++) {
        qvec->axpy(ritzvecs[n * (max_jd_size + 1) + j], V[j]);
      }
      oper->applyBCs(qvec);

      // Normalize the predicted eigenvalue
      TacsScalar qnorm = sqrt(oper->dot(qvec, qvec));
      qvec->scale(1.0 / qnorm);

      if (error) {
        // Compute the residual: work = A*q - theta*B*q
        // and store it in the work vector
        oper->multA(qvec, work);
        oper->multB(qvec, pvec);
        work->axpy(-theta, pvec);
        oper->applyBCs(work);

        // Set the error
        *error = work->norm();
      }
      if (Q[n] != ans) {
        Q[n]->copyValues(ans);
      }
    }

    return theta;
  }

  return 0.0;
}

/*
  Solve the general eigenproblem using the Jacobi-Davidson method.

  This code implements a technique to find the absolute smallest
  eigenvalues. The code maintains a B-orthonormal basis V that
  is used to form a Galerkin approximation of the eigenvalue problem

  V^{T}*A*V*y - theta*V^{T}*B*V*y = M*y - theta*y = 0.

  The updates to the subspace V are obtained by seeking a
  Jacobi-Davidson update that is B-orthogonal to the current
  eigenvector estimate. To obtain multiple eigenvectors the code
  uses a deflation subspace Q which consists of orthogonal
  eigenvectors that have converged. In addition, the current estimate
  of the eigenvector is maintained in the same array of vectors.
  At each iteration, we seek a new vector t that approximately
  satisfies the following equation

  (I - P*Q^{T})*A*(I - Q*P^{T})*t = -r

  where r is orthogonal to the subspace Q and t is orthogonal to the
  subspace t by construction. FGMRES builds the subspaces W and Z.
*/
void TACSJacobiDavidson::solve(KSMPrint *ksm_print, int print_level) {
  // Keep track of the current subspace
  V[0]->setRand(-1.0, 1.0);
  oper->applyBCs(V[0]);

  // The maximum size of the Jacobi--Davidson subspace
  const int m = max_jd_size;
  memset(M, 0, m * m * sizeof(TacsScalar));

  // Allocate the space for the real work vector
  int lwork = 16 * max_jd_size;
  double *rwork = new double[lwork];

  // Store the number of iterations of the inner GMRES algorithm
  int iteration = 0;
  int gmres_iteration = 0;

  // Check if the recycle flag is set
  int kstart = 0;
  if (num_recycle_vecs > 0) {
    // Set the actual number of recycled eigenvectors
    int num_recycle = num_recycle_vecs;
    if (num_recycle > nconverged) {
      num_recycle = nconverged;
    }
    if (num_recycle > 0) {
      if (recycle_type == JD_SUM_TWO) {
        V[0]->setRand(-1.0, 1.0);
        V[1]->zeroEntries();
        for (int i = 0; i < nconverged; i++) {
          V[1]->axpy(1.0, Q[i]);
        }

        for (int k = 0; k < 2; k++) {
          // B-orthogonalize the eigenvectors
          for (int i = 0; i < k; i++) {
            TacsScalar h = oper->dot(V[k], V[i]);
            V[k]->axpy(-h, V[i]);
          }

          // Apply boundary conditions for this vector
          oper->applyBCs(V[k]);

          // Normalize the vector so that it is orthonormal
          TacsScalar vnorm = sqrt(oper->dot(V[k], V[k]));
          V[k]->scale(1.0 / vnorm);

          // Compute work = A*V[k]
          oper->multA(V[k], work);

          // Complete the entries in the symmetric matrix M that is formed by
          // M = V^{T}*A*V
          for (int i = 0; i <= k; i++) {
            M[k * m + i] = V[i]->dot(work);
            M[i * m + k] = M[k * m + i];
          }
        }
        kstart = 2;
      } else if (recycle_type == JD_NUM_RECYCLE) {
        // B-orthogonalize the old eigenvectors with respect to the
        // new matrix for all but the last recycled eigenvector which
        // will be orthogonalized by the first iteration through the
        // solution loop.
        for (int k = 0; k < num_recycle; k++) {
          // Copy the vector from the old eigenvector
          if (k >= 1) {
            V[k]->copyValues(Q[k - 1]);
          }

          // B-orthogonalize the eigenvectors
          for (int i = 0; i < k; i++) {
            TacsScalar h = oper->dot(V[k], V[i]);
            V[k]->axpy(-h, V[i]);
          }

          // Apply boundary conditions for this vector
          oper->applyBCs(V[k]);

          // Normalize the vector so that it is orthonormal
          TacsScalar vnorm = sqrt(oper->dot(V[k], V[k]));
          V[k]->scale(1.0 / vnorm);

          // Compute work = A*V[k]
          oper->multA(V[k], work);

          // Complete the entries in the symmetric matrix M that is formed by
          // M = V^{T}*A*V
          for (int i = 0; i <= k; i++) {
            M[k * m + i] = V[i]->dot(work);
            M[i * m + k] = M[k * m + i];
          }
        }

        // Copy over the last eigenvector
        kstart = num_recycle;
        V[kstart]->copyValues(Q[num_recycle - 1]);
      }
    }
  }

  // Reset the number of converged eigenvectors/eigenvalues
  nconverged = 0;

  if (ksm_print && print_level > 0) {
    char line[256];
    sprintf(line, "%4s %15s %15s %10s\n", "Iter", "JD Residual", "Ritz value",
            "toler");
    ksm_print->print(line);
  }

  for (int k = kstart; k < m; k++) {
    // Orthogonalize against converged eigenvectors as well (if any)...
    for (int i = 0; i < nconverged; i++) {
      TacsScalar h = oper->dot(V[k], Q[i]);
      V[k]->axpy(-h, Q[i]);
    }

    // Orthogonalize V[k] against all other vectors in the current
    // solution subspace. This ensures that the vectors are orthogonal
    // with respect to the operator.
    for (int i = 0; i < k; i++) {
      TacsScalar h = oper->dot(V[k], V[i]);
      V[k]->axpy(-h, V[i]);
    }

    // Apply boundary conditions for this vector
    oper->applyBCs(V[k]);

    // Normalize the vector so that it is orthonormal
    TacsScalar vnorm = sqrt(oper->dot(V[k], V[k]));
    V[k]->scale(1.0 / vnorm);

    // Compute work = A*V[k]
    oper->multA(V[k], work);

    // Complete the entries in the symmetric matrix M that is formed
    // by M = V^{T}*A*V
    for (int i = 0; i <= k; i++) {
      M[k * m + i] = V[i]->dot(work);
      M[i * m + k] = M[k * m + i];
    }

    // Compute the eigenvalues/eigenvectors of the M matrix. Copy over
    // the values from the M matrix into the ritzvecs array.
    for (int j = 0; j <= k; j++) {
      for (int i = 0; i <= k; i++) {
        ritzvecs[i + (k + 1) * j] = TacsRealPart(M[i + m * j]);
      }
    }

    // The input arguments required for LAPACK
    const char *jobz = "V", *uplo = "U";

    // Compute the eigenvalues and eigenvectors of a (small) symmetric
    // matrix using lapack. The eigenvalues are the ritzvalues. The
    // eigenvectors can be used to construct the Ritz vectors (in the
    // code the ritzvecs array contains the vectors of the reduced
    // eigenproblem that are not Ritz vectors themselves but are used
    // to compute them.)
    int info;
    int n = k + 1;
    LAPACKdsyev(jobz, uplo, &n, ritzvecs, &n, ritzvals, rwork, &lwork, &info);

    // The Ritz value that estimates the smallest eigenvalue that
    // has not yet converged
    double theta = ritzvals[0];

    int num_new_eigvals = 0;
    for (int i = 0; i < k + 1 && nconverged < max_eigen_vectors;
         i++, num_new_eigvals++) {
      // Set the Ritz value
      theta = ritzvals[i];

      // Assemble the predicted Ritz vector
      Q[nconverged]->zeroEntries();
      for (int j = 0; j <= k; j++) {
        Q[nconverged]->axpy(ritzvecs[i * (k + 1) + j], V[j]);
      }
      oper->applyBCs(Q[nconverged]);

      // Normalize the predicted eigenvalue
      TacsScalar qnorm = sqrt(oper->dot(Q[nconverged], Q[nconverged]));
      Q[nconverged]->scale(1.0 / qnorm);

      // Compute the residual: work = A*q - theta*B*q
      // and store it in the work vector
      oper->multA(Q[nconverged], work);
      TacsScalar Anorm = work->norm();
      oper->multB(Q[nconverged], P[nconverged]);
      work->axpy(-theta, P[nconverged]);
      oper->applyBCs(work);

      // Compute the norm of the residual: ||A*q - theta*B*q||
      TacsScalar w_norm = work->norm();

      // Compute the norm of the eigenvalue to check if it has converged
      double abs_theta = fabs(TacsRealPart(theta));
      double toler = (eig_atol * (theta_cutoff / (theta_cutoff + abs_theta)) +
                      eig_rtol * (abs_theta / (theta_cutoff + abs_theta)) *
                          TacsRealPart(Anorm));

      if (ksm_print && print_level > 0) {
        char line[256];
        sprintf(line, "%4d %15.5e %15.5e %10.2e\n", iteration,
                TacsRealPart(w_norm), theta, toler);
        ksm_print->print(line);
      }

      if (TacsRealPart(w_norm) <= toler) {
        // Record the Ritz value as the eigenvalue
        eigvals[nconverged] = theta;
        eigerror[nconverged] = w_norm;

        nconverged++;
      } else {
        // The eigenvalue did not converge, repeat the computation
        break;
      }
    }

    // Check if we should quit the loop or continue
    if (nconverged >= max_eigen_vectors) {
      break;
    }

    if (num_new_eigvals > 0) {
      // This eigenvalue has converged, store it in Q[nconverged] and
      // modify the remaining entries of V. Now we need to compute the
      // new starting vector for the next eigenvalue

      // The number of vectors in the new sub-space. This new subspace
      // is constructed from the Ritz problem from the old subspace stored
      // currently in V.
      int new_subspace_vecs = (k + 1) - num_new_eigvals;
      if (new_subspace_vecs > max_gmres_size + 1) {
        new_subspace_vecs = max_gmres_size + 1;
      }

      // Reset the matrix to zero
      memset(M, 0, m * m * sizeof(TacsScalar));

      // Set the vectors that will be used
      for (int i = 0; i < new_subspace_vecs; i++) {
        // Keep only the ritz values that are not yet converged
        int ritz_index = i + num_new_eigvals;

        // Set the ritz value
        M[i + i * m] = ritzvals[ritz_index];

        // Build the new subspace vector from the old V subspace and
        // store it temporarily in W
        W[i]->zeroEntries();
        for (int j = 0; j <= k; j++) {
          W[i]->axpy(ritzvecs[ritz_index * (k + 1) + j], V[j]);
        }

        // Orthogonalize the new starting vector against all other
        // converged eigenvalues
        for (int j = 0; j < nconverged; j++) {
          TacsScalar h = oper->dot(W[i], Q[j]);
          W[i]->axpy(-h, Q[j]);
        }

        // Orthogonalize the new subspace values against themselves
        for (int j = 0; j < i; j++) {
          TacsScalar h = oper->dot(W[i], W[j]);
          W[i]->axpy(-h, W[j]);
        }

        // Apply the boundary conditions
        oper->applyBCs(W[i]);

        // Normalize the vector so that it is orthonormal
        TacsScalar vnorm = sqrt(oper->dot(W[i], W[i]));
        W[i]->scale(1.0 / vnorm);
      }

      // Copy the new subspace vectors to V
      for (int i = 0; i < new_subspace_vecs; i++) {
        V[i]->copyValues(W[i]);
      }

      // Now add one new subspace vector that is randomly generated
      if (new_subspace_vecs < m) {
        int index = new_subspace_vecs;

        // Generate a new random vector to append to the space
        V[index]->setRand(-1.0, 1.0);
        oper->applyBCs(V[index]);

        // B-orthogonalize the eigenvectors against the subspace
        for (int i = 0; i < new_subspace_vecs; i++) {
          TacsScalar h = oper->dot(V[index], V[i]);
          V[index]->axpy(-h, V[i]);
        }

        // Apply the boundary conditions
        oper->applyBCs(V[index]);

        // Compute work = A*V[k]
        oper->multA(V[index], work);

        // Complete the entries in the symmetric matrix M that is formed by
        // M = V^{T}*A*V
        for (int i = 0; i <= index; i++) {
          M[index * m + i] = V[i]->dot(work);
          M[i * m + index] = M[index * m + i];
        }

        // Increment the number of new subspace vectors by one to account
        // for the extra random vector
        new_subspace_vecs++;
      }

      // Reset the iteration counter here to reflect the new size
      // of the V subspace
      k = new_subspace_vecs - 2;

      // Reset the iteration loop and continue
      continue;
    }

    // Now solve the system (K - theta*M)*t = -work
    // Keep track of the number of iterations in GMRES
    int niters = 0;

    // Copy the residual to the first work vector
    W[0]->copyValues(work);
    res[0] = W[0]->norm();
    W[0]->scale(1.0 / res[0]);  // W[0] = b/|| b ||

    // Keep track of the initial norm of the right-hand-side
    double beta = TacsRealPart(res[0]);

    // Using GMRES, solve for the update equation
    for (int i = 0; i < max_gmres_size; i++) {
      // Now compute: work = (I - Q*P^{T})*Z[i]
      // Note that we compute the product in this way since
      // it is numerical more stable and equivalent to the form above
      // since Q^{T}P = 0 so that: (since Q is a B-orthogonal subspace)
      // (I - q1*p1^{T})*(I - q2*p2^{T}) =
      // (I - q1*p1^{T} - q2*p2^{T} - q1*p1^{T}*q2*p2^{T}) =
      // (I - q1*p1^{T} - q2*p2^{T})

      // Apply the preconditioner, Z[i] = (I - Q*P^{T})*M^{-1} W[i]
      oper->applyFactor(W[i], Z[i]);
      for (int j = 0; j <= nconverged; j++) {
        TacsScalar h = Z[i]->dot(P[j]);
        Z[i]->axpy(-h, Q[j]);
      }

      // Form the product (A - theta*B)*Z[i]
      oper->multA(Z[i], W[i + 1]);
      oper->multB(Z[i], work);
      W[i + 1]->axpy(-theta, work);
      oper->applyBCs(W[i + 1]);

      // Finish computing the product W[i+1] = (I - P*Q^{T})*A*Z[i]
      for (int j = 0; j <= nconverged; j++) {
        TacsScalar h = W[i + 1]->dot(Q[j]);
        W[i + 1]->axpy(-h, P[j]);
      }

      // Build the orthogonal basis using MGS
      for (int j = i; j >= 0; j--) {
        H[j + Hptr[i]] = W[i + 1]->dot(W[j]);   // H[j,i] = dot(W[i+1], W[i])
        W[i + 1]->axpy(-H[j + Hptr[i]], W[j]);  // W[i+1] = W[i+1] - H[j,i]*W[j]
      }

      // Complete the basis
      H[i + 1 + Hptr[i]] = W[i + 1]->norm();  // H[i+1,i] = || W[i+1] ||
      W[i + 1]->scale(1.0 /
                      H[i + 1 + Hptr[i]]);  // W[i+1] = W[i+1]/|| W[i+1] ||

      // Apply the existing part of Q to the new components of the
      // Hessenberg matrix
      TacsScalar h1, h2;
      for (int k = 0; k < i; k++) {
        h1 = H[k + Hptr[i]];
        h2 = H[k + 1 + Hptr[i]];
        H[k + Hptr[i]] = h1 * Qcos[k] + h2 * Qsin[k];
        H[k + 1 + Hptr[i]] = -h1 * Qsin[k] + h2 * Qcos[k];
      }

      // Now, compute the rotation for the new column that was just added
      h1 = H[i + Hptr[i]];
      h2 = H[i + 1 + Hptr[i]];
      TacsScalar sq = sqrt(h1 * h1 + h2 * h2);

      // Evaluate the sin/cos of the rotation matrix
      Qcos[i] = h1 / sq;
      Qsin[i] = h2 / sq;
      H[i + Hptr[i]] = h1 * Qcos[i] + h2 * Qsin[i];
      H[i + 1 + Hptr[i]] = -h1 * Qsin[i] + h2 * Qcos[i];

      // Update the residual
      h1 = res[i];
      res[i] = h1 * Qcos[i];
      res[i + 1] = -h1 * Qsin[i];

      if (ksm_print && print_level > 1) {
        ksm_print->printResidual(niters, fabs(TacsRealPart(res[i + 1])));
      }

      niters++;           // Number of iterations for this GMRES solve
      gmres_iteration++;  // Total number of GMRES iterations

      if (fabs(TacsRealPart(res[i + 1])) < atol ||
          fabs(TacsRealPart(res[i + 1])) < rtol * beta) {
        break;
      }
    }

    // Now, compute the solution. The linear combination of the
    // Arnoldi vectors. The matrix H is now upper triangular. First
    // compute the weights for each basis vector.
    for (int i = niters - 1; i >= 0; i--) {
      for (int j = i + 1; j < niters; j++) {
        res[i] = res[i] - H[i + Hptr[j]] * res[j];
      }
      res[i] = res[i] / H[i + Hptr[i]];
    }

    // Compute the next basis vector for the outer Jacobi--Davidson basis
    V[k + 1]->zeroEntries();
    for (int i = 0; i < niters; i++) {
      V[k + 1]->axpy(-res[i], Z[i]);
    }

    // Compute the product to test the error
    // (1 - P*Q^{T})*(A - theta*B)*(1 - Q*P^{T})
    W[0]->copyValues(V[k + 1]);
    for (int j = 0; j <= nconverged; j++) {
      TacsScalar h = W[0]->dot(P[j]);
      W[0]->axpy(-h, Q[j]);
    }

    iteration++;
  }

  // Sort the indices for the converged eigenvalues
  TacsArgSort(nconverged, eigvals, eigindex);

  if (ksm_print) {
    char line[256];
    for (int i = 0; i < nconverged; i++) {
      int index = eigindex[i];
      sprintf(line, "Eigenvalue[%2d]: %25.10e Eig. error[%2d]: %25.10e\n", i,
              TacsRealPart(eigvals[index]), i, TacsRealPart(eigerror[index]));
      ksm_print->print(line);
    }

    sprintf(line, "JD number of outer iterations: %2d\n", iteration);
    ksm_print->print(line);
    sprintf(line, "JD number of inner GMRES iterations: %2d\n",
            gmres_iteration);
    ksm_print->print(line);
  }

  delete[] rwork;
}

/*
  Set the relative and absolute tolerances used for the stopping
  criterion.

  input:
  eig_rtol: the eigenvector tolerance ||(A - theta*B)*eigvec|| <
  eig_rtol*||A*eigvec|| eig_atol: the eigenvector tolerance ||(A -
  theta*B)*eigvec|| < eig_atol rtol: the relative tolerance ||r_k|| <
  rtol*||r_0|| atol: the absolute tolerancne ||r_k|| < atol
*/
void TACSJacobiDavidson::setTolerances(double _eig_rtol, double _eig_atol,
                                       double _rtol, double _atol) {
  eig_rtol = _eig_rtol;
  eig_atol = _eig_atol;
  rtol = _rtol;
  atol = _atol;
}

/*
  Set the paramter that decide whether the convergence check relies more on
  eig_atol or more on eig_rtol.

  input:
  theta_cutoff: between 0.0 and 1.0, usually is a small value e.g. 0.1 or 0.01,
                the smaller it is, the convergence criterion behaves more like
                a step function between eig_atol and eig_rtol*||A*eigvec||
*/
void TACSJacobiDavidson::setThetaCutoff(double _theta_cutoff) {
  theta_cutoff = _theta_cutoff;
}

/*
  Set the number of vectors to recycle if the eigenvectors are converged

  input:
  num_recycle_vecs: number of vectors to recycle
*/
void TACSJacobiDavidson::setRecycle(int _num_recycle_vecs,
                                    JDRecycleType _recycle_type) {
  num_recycle_vecs = _num_recycle_vecs;
  recycle_type = _recycle_type;
}
