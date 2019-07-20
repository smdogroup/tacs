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

#include "tacslapack.h"
#include "JacobiDavidson.h"


TACSJDFrequencyOperator::TACSJDFrequencyOperator( TACSAssembler *_tacs, 
                                                  TACSMat *_kmat,
                                                  TACSMat *_mmat,
                                                  TACSMat *_pc_mat,
                                                  TACSPc *_pc ){
  tacs = _tacs;  tacs->incref();
  kmat = _kmat;  kmat->incref();
  mmat = _mmat;  mmat->incref();
  pc_mat = _pc_mat;  pc_mat->incref();
  pc = _pc;  pc->incref();
  work = tacs->createVec();
  work->incref();
}

TACSJDFrequencyOperator::~TACSJDFrequencyOperator(){
  tacs->decref();
  kmat->decref();
  mmat->decref();
  pc_mat->decref();
  pc->decref();
  work->decref();
}

// Create a vector
TACSVec *TACSJDFrequencyOperator::createVec(){
  return tacs->createVec();
}

// Set the eigenvalue estimate (and reset the factorization)
// pc = K in this case
void TACSJDFrequencyOperator::setEigenvalueEstimate( double estimate ){
  pc_mat->zeroEntries();
  pc_mat->copyValues(kmat);
  if (estimate != 0.0){
    pc_mat->axpy(-estimate, mmat);
  }
  tacs->applyBCs(pc_mat);
  pc->factor();
}

// Apply the preconditioner
void TACSJDFrequencyOperator::applyFactor( TACSVec *x, TACSVec *y ){
  pc->applyFactor(x, y);
  tacs->applyBCs(y);
}

// Apply boundary conditions associated with the matrix
void TACSJDFrequencyOperator::applyBCs( TACSVec *x ){
  tacs->applyBCs(x);
}

// Perform a dot-product of two vectors in the operator space e.g.
// output <- x^{T}*B*y
TacsScalar TACSJDFrequencyOperator::dot( TACSVec *x, TACSVec *y ){
  mmat->mult(x, work);
  return y->dot(work);
}

// Matrix-vector products with either the A or B operators.
void TACSJDFrequencyOperator::multA( TACSVec *x, TACSVec * y ){
  kmat->mult(x, y);
  tacs->applyBCs(y);
}

void TACSJDFrequencyOperator::multB( TACSVec *x, TACSVec * y ){
  mmat->mult(x, y);
  tacs->applyBCs(y);
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
TACSJacobiDavidson::TACSJacobiDavidson( TACSJacobiDavidsonOperator *_oper,
                                        int _max_eigen_vectors,
                                        int _max_jd_size,
                                        int _max_gmres_size ){
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
  recycle = 0;
  recycle_type = JD_NUM_RECYCLE;

  // The eigen tolerance
  eigtol = 1e-9;

  // The residual tolerances for the GMRES iterations
  rtol = 1e-9;
  atol = 1e-30;

  // Store the matrix
  M = new TacsScalar[ (max_jd_size+1)*(max_jd_size+1) ];
  ritzvecs = new double[ (max_jd_size+1)*(max_jd_size+1) ];
  ritzvals = new double[ (max_jd_size+1) ];

  // Allocate space for the vectors
  V = new TACSVec*[ max_jd_size+1 ];
  Q = new TACSVec*[ max_eigen_vectors+1 ];
  P = new TACSVec*[ max_eigen_vectors+1 ];
  eigvals = new TacsScalar[ max_eigen_vectors ];

  // Allocate the variables
  for ( int i = 0; i < max_jd_size+1; i++ ){
    V[i] = oper->createVec();
    V[i]->incref();
  }
  for ( int i = 0; i < max_eigen_vectors+1; i++ ){
    Q[i] = oper->createVec();
    P[i] = oper->createVec();
    Q[i]->incref();
    P[i]->incref();
  }

  // Create and populate the pointer into the H columns
  Hptr = new int[ max_gmres_size+1 ];
  Hptr[0] = 0;
  for ( int i = 0; i < max_gmres_size; i++ ){
    Hptr[i+1] = Hptr[i] + i+2;
  }

  // Allocate the Hessenberg matrix
  H = new TacsScalar[ Hptr[max_gmres_size] ];
  res = new TacsScalar[ max_gmres_size+1 ];

  // Allocate space for the unitary matrix Q
  Qsin = new TacsScalar[ max_gmres_size ];
  Qcos = new TacsScalar[ max_gmres_size ];

  // Allocate space for the vectors
  W = new TACSVec*[ max_gmres_size+1 ];
  for ( int i = 0; i < max_gmres_size+1; i++ ){
    W[i] = oper->createVec();
    W[i]->incref();
  }

  // Create the work vector
  work = oper->createVec();
  work->incref();

  Z = new TACSVec*[ max_gmres_size ];
  for ( int i = 0; i < max_gmres_size; i++ ){
    Z[i] = oper->createVec();
    Z[i]->incref();
  }
}

/*
  Free the data associated with the object
*/
TACSJacobiDavidson::~TACSJacobiDavidson(){
  // Free the data to solve the GMRES problem
  delete [] H;
  delete [] Hptr;
  delete [] res;
  delete [] Qsin;
  delete [] Qcos;

  // Free the Jacobi-Davidson vectors
  for ( int i = 0; i < max_jd_size+1; i++ ){
    V[i]->decref();
  }
  delete [] V;
  for ( int i = 0; i < max_eigen_vectors+1; i++ ){
    Q[i]->decref();
    P[i]->decref();
  }
  delete [] Q;
  delete [] P;

  // Free the JD matrix
  delete [] M;
  delete [] ritzvecs;
  delete [] ritzvals;
  delete [] eigvals;

  // Free the GMRES subspace vectors
  W[0]->decref();
  for ( int i = 0; i < max_gmres_size; i++ ){
    W[i+1]->decref();
  }
  delete [] W;

  for ( int i = 0; i < max_gmres_size; i++ ){
    Z[i]->decref();
  }
  delete [] Z;

  work->decref();
}

/*!
  Extract the eigenvalue from the analysis
*/
TacsScalar TACSJacobiDavidson::extractEigenvalue( int n,
                                                  TacsScalar *error ){
  if (error){
    *error = 0.0;
  }

  if (n >= 0 && n < max_eigen_vectors){
    if (error){
      oper->multA(Q[n], work);
      work->axpy(-eigvals[n], P[n]);
      *error = work->norm();
    }

    return eigvals[n];
  }

  return 0.0;
}

/*!
  Extract the eigenvector and eigenvalue from the eigenvalue analysis
*/
TacsScalar TACSJacobiDavidson::extractEigenvector( int n, TACSVec *ans,
                                                   TacsScalar *error ){
  if (error){
    *error = 0.0;
  }

  if (n >= 0 && n < max_eigen_vectors){
    if (error){
      oper->multA(Q[n], work);
      work->axpy(-eigvals[n], P[n]);
      *error = work->norm();
    }
    if (ans){
      ans->copyValues(Q[n]);
    }

    return eigvals[n];
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
void TACSJacobiDavidson::solve( KSMPrint *ksm_print, 
                                KSMPrint *ksm_file ){
  // Keep track of the current subspace
  V[0]->setRand(-1.0, 1.0);
  oper->applyBCs(V[0]);

  // The maximum size of the Jacobi--Davidson subspace
  const int m = max_jd_size;
  memset(M, 0, m*m*sizeof(TacsScalar));

  // Allocate the space for the real work vector
  int lwork = 16*max_jd_size;
  double *rwork = new double[ lwork ];

  // Store the number of converged eigenvalues
  int iteration = 0;

  // Check if the recycle flag is set
  int kstart = 0;
  if (recycle > 0){
    // Set the actual number of recycled eigenvectors
    int num_recycle = recycle;
    if (num_recycle > nconverged){
      num_recycle = nconverged;
    }
    if (num_recycle > 0){
      if (recycle_type == JD_SUM_TWO){
        V[0]->setRand(-1.0, 1.0);
        V[1]->zeroEntries();
        for (int i = 0; i < nconverged; i++){
          V[1]->axpy(1.0, Q[i]);
        }
        for ( int k = 0; k < 2; k++ ){
          // B-orthogonalize the eigenvectors
          for ( int i = 0; i < k; i++ ){
            TacsScalar h = oper->dot(V[k], V[i]);
            V[k]->axpy(-h, V[i]);
          }
          // Apply boundary conditions for this vector
          oper->applyBCs(V[k]);
          
          // Normalize the vector so that it is orthonormal
          TacsScalar vnorm = sqrt(oper->dot(V[k], V[k]));
          V[k]->scale(1.0/vnorm);
          // Compute work = A*V[k]
          oper->multA(V[k], work);
          
          // Complete the entries in the symmetric matrix M that is formed by 
          // M = V^{T}*A*V
          for ( int i = 0; i <= k; i++ ){
            M[k*m + i] = V[i]->dot(work);
            M[i*m + k] = M[k*m + i];
          }
        }
        kstart = 2;      
      }
      else if (recycle_type == JD_NUM_RECYCLE){ 
        // B-orthogonalize the old eigenvectors with respect to the
        // new matrix for all but the last recycled eigenvector which
        // will be orthogonalized by the first iteration through the
        // solution loop.
        for ( int k = 0; k < num_recycle; k++ ){
          // Copy the vector from the old eigenvector
          if (k >= 1){
            V[k]->copyValues(Q[k-1]);
          }

          // B-orthogonalize the eigenvectors
          for ( int i = 0; i < k; i++ ){
            TacsScalar h = oper->dot(V[k], V[i]);
            V[k]->axpy(-h, V[i]);
          }

          // Apply boundary conditions for this vector
          oper->applyBCs(V[k]);

          // Normalize the vector so that it is orthonormal
          TacsScalar vnorm = sqrt(oper->dot(V[k], V[k]));
          V[k]->scale(1.0/vnorm);

          // Compute work = A*V[k]
          oper->multA(V[k], work);

          // Complete the entries in the symmetric matrix M that is formed by 
          // M = V^{T}*A*V
          for ( int i = 0; i <= k; i++ ){
            M[k*m + i] = V[i]->dot(work);
            M[i*m + k] = M[k*m + i];
          }
        }

        // Copy over the last eigenvector
        kstart = num_recycle;
        V[kstart]->copyValues(Q[num_recycle-1]);
      }
    }
  }

  // Reset the number of converged eigenvectors/eigenvalues
  nconverged = 0;

  for ( int k = kstart; k < m; k++ ){
    // Orthogonalize against converged eigenvectors as well (if any)...
    for ( int i = 0; i < nconverged; i++ ){
      TacsScalar h = oper->dot(V[k], Q[i]);
      V[k]->axpy(-h, Q[i]);
    }

    // Orthogonalize V[k] against all other vectors in the current
    // solution subspace. This ensures that the vectors are orthogonal
    // with respect to the operator.
    for ( int i = 0; i < k; i++ ){
      TacsScalar h = oper->dot(V[k], V[i]);
      V[k]->axpy(-h, V[i]);
    }

    // Apply boundary conditions for this vector
    oper->applyBCs(V[k]);

    // Normalize the vector so that it is orthonormal
    TacsScalar vnorm = sqrt(oper->dot(V[k], V[k]));
    V[k]->scale(1.0/vnorm);

    // Compute work = A*V[k]
    oper->multA(V[k], work);

    // Complete the entries in the symmetric matrix M that is formed
    // by M = V^{T}*A*V
    for ( int i = 0; i <= k; i++ ){
      M[k*m + i] = V[i]->dot(work);
      M[i*m + k] = M[k*m + i];
    }

    // Compute the eigenvalues/eigenvectors of the M matrix. Copy over
    // the values from the M matrix into the ritzvecs array.
    for ( int j = 0; j <= k; j++ ){
      for ( int i = 0; i <= k; i++ ){
        ritzvecs[i + (k+1)*j] = TacsRealPart(M[i + m*j]);
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
    int n = k+1;
    LAPACKdsyev(jobz, uplo, &n, ritzvecs, &n, ritzvals,
                rwork, &lwork, &info);

    // The Ritz value that estimates the smallest eigenvalue that
    // has not yet converged
    double theta = ritzvals[0];

    int num_new_eigvals = 0;
    for ( int i = 0; i < k+1 && nconverged < max_eigen_vectors; 
          i++, num_new_eigvals++ ){
      // Set the Ritz value
      theta = ritzvals[i];

      // Assemble the predicted Ritz vector
      Q[nconverged]->zeroEntries();
      for ( int j = 0; j <= k; j++ ){
        Q[nconverged]->axpy(ritzvecs[i*(k+1) + j], V[j]);
      }
      oper->applyBCs(Q[nconverged]);

      // Normalize the predicted eigenvalue
      TacsScalar qnorm = sqrt(oper->dot(Q[nconverged], Q[nconverged]));
      Q[nconverged]->scale(1.0/qnorm);

      // Compute the residual: work = A*q - theta*B*q 
      // and store it in the work vector
      oper->multA(Q[nconverged], work);
      TacsScalar Anorm = work->norm();
      oper->multB(Q[nconverged], P[nconverged]);
      work->axpy(-theta, P[nconverged]);
      oper->applyBCs(work);
      
      TacsScalar w_norm = work->norm();
      if (ksm_print){
        char line[256];
        sprintf(line, "JD Residual[%2d]: %15.5e  Eigenvalue[%2d]: %20.10e\n", 
                iteration, TacsRealPart(w_norm), nconverged, theta);
        ksm_print->print(line);
      }

      // Compute the norm of the eigenvalue to check if it has converged
      if (TacsRealPart(w_norm) <= TacsRealPart(eigtol*Anorm)){
        // Record the Ritz value as the eigenvalue
        eigvals[nconverged] = theta;

        nconverged++;
      }
      else {
        // The eigenvalue did not converge, repeat the computation
        break;
      }
    }

    // Check if we should quit the loop or continue
    if (nconverged >= max_eigen_vectors){
      break;
    }

    if (num_new_eigvals > 0){
      // This eigenvalue has converged, store it in Q[nconverged] and
      // modify the remaining entries of V. Now we need to compute the
      // new starting vector for the next eigenvalue
      int max_new_vecs = k+1 - num_new_eigvals;
      if (max_new_vecs > max_gmres_size+1){
        max_new_vecs = max_gmres_size+1;
      }

      // Reset the matrix to zero
      memset(M, 0, m*m*sizeof(TacsScalar));

      // Set the vectors that will be used
      for ( int i = 0; i < max_new_vecs; i++ ){
        int ritz_index = i + num_new_eigvals;

        // Set the ritz value
        M[i + i*m] = ritzvals[ritz_index];

        W[i]->zeroEntries();
        for ( int j = 0; j <= k; j++ ){
          W[i]->axpy(ritzvecs[ritz_index*(k+1) + j], V[j]);
        }

        // Orthogonalize the new starting vector against all other
        // converged eigenvalues
        for ( int j = 0; j < nconverged; j++ ){
          TacsScalar h = oper->dot(W[i], Q[j]);
          W[i]->axpy(-h, Q[j]);
        }

        for ( int j = 0; j < i; j++ ){
          TacsScalar h = oper->dot(W[i], W[j]);
          W[i]->axpy(-h, W[j]);
        }
        
        oper->applyBCs(W[i]);

        // Normalize the vector so that it is orthonormal
        TacsScalar vnorm = sqrt(oper->dot(W[i], W[i]));
        W[i]->scale(1.0/vnorm);
      }

      for ( int i = 0; i < max_new_vecs; i++ ){
        V[i]->copyValues(W[i]);
      }
      k = max_new_vecs-2;

      // Reset the iteration loop and continue
      continue;
    }

    // Now solve the system (K - theta*M)*t = -work
    // Keep track of the number of iterations in GMRES
    int niters = 0;

    // Copy the residual to the first work vector
    W[0]->copyValues(work);
    res[0] = W[0]->norm();
    W[0]->scale(1.0/res[0]); // W[0] = b/|| b ||

    // Keep track of the initial norm of the right-hand-side
    double beta = TacsRealPart(res[0]);

    // Using GMRES, solve for the update equation
    for ( int i = 0; i < max_gmres_size; i++ ){
      // Now compute: work = (I - Q*P^{T})*Z[i]
      // Note that we compute the product in this way since
      // it is numerical more stable and equivalent to the form above
      // since Q^{T}P = 0 so that: (since Q is a B-orthogonal subspace)
      // (I - q1*p1^{T})*(I - q2*p2^{T}) = 
      // (I - q1*p1^{T} - q2*p2^{T} - q1*p1^{T}*q2*p2^{T}) =
      // (I - q1*p1^{T} - q2*p2^{T})

      // Apply the preconditioner, Z[i] = (I - Q*P^{T})*M^{-1} W[i]
      oper->applyFactor(W[i], Z[i]);
      for ( int j = 0; j <= nconverged; j++ ){
        TacsScalar h = Z[i]->dot(P[j]);
        Z[i]->axpy(-h, Q[j]);
      }

      // Form the product (A - theta*B)*Z[i]
      oper->multA(Z[i], W[i+1]);
      oper->multB(Z[i], work);
      W[i+1]->axpy(-theta, work);
      oper->applyBCs(W[i+1]);

      // Finish computing the product W[i+1] = (I - P*Q^{T})*A*Z[i]
      for ( int j = 0; j <= nconverged; j++ ){
        TacsScalar h = W[i+1]->dot(Q[j]);
        W[i+1]->axpy(-h, P[j]);
      }

      // Build the orthogonal basis using MGS
      for ( int j = i; j >= 0; j-- ){
        H[j + Hptr[i]] = W[i+1]->dot(W[j]); // H[j,i] = dot(W[i+1], W[i])
        W[i+1]->axpy(-H[j + Hptr[i]], W[j]); // W[i+1] = W[i+1] - H[j,i]*W[j]
      }

      // Complete the basis
      H[i+1 + Hptr[i]] = W[i+1]->norm(); // H[i+1,i] = || W[i+1] ||
      W[i+1]->scale(1.0/H[i+1 + Hptr[i]]); // W[i+1] = W[i+1]/|| W[i+1] ||

      // Apply the existing part of Q to the new components of the
      // Hessenberg matrix
      TacsScalar h1, h2;
      for ( int k = 0; k < i; k++ ){
        h1 = H[k   + Hptr[i]];
        h2 = H[k+1 + Hptr[i]];
        H[k   + Hptr[i]] =  h1*Qcos[k] + h2*Qsin[k];
        H[k+1 + Hptr[i]] = -h1*Qsin[k] + h2*Qcos[k];
      }

      // Now, compute the rotation for the new column that was just added
      h1 = H[i   + Hptr[i]];
      h2 = H[i+1 + Hptr[i]];
      TacsScalar sq = sqrt(h1*h1 + h2*h2);

      // Evaluate the sin/cos of the rotation matrix
      Qcos[i] = h1/sq;
      Qsin[i] = h2/sq;
      H[i   + Hptr[i]] =  h1*Qcos[i] + h2*Qsin[i];
      H[i+1 + Hptr[i]] = -h1*Qsin[i] + h2*Qcos[i];

      // Update the residual
      h1 = res[i];
      res[i]   =   h1*Qcos[i];
      res[i+1] = - h1*Qsin[i];

      // if (ksm_print){
      //   ksm_print->printResidual(niters, fabs(TacsRealPart(res[i+1])));
      // }

      niters++;

      if (fabs(TacsRealPart(res[i+1])) < atol ||
          fabs(TacsRealPart(res[i+1])) < rtol*beta){
        break;
      }
    }

    // Now, compute the solution. The linear combination of the
    // Arnoldi vectors. The matrix H is now upper triangular. First
    // compute the weights for each basis vector.
    for ( int i = niters-1; i >= 0; i-- ){
      for ( int j = i+1; j < niters; j++ ){
        res[i] = res[i] - H[i + Hptr[j]]*res[j];
      }
      res[i] = res[i]/H[i + Hptr[i]];
    }

    // Compute the next basis vector for the outer Jacobi--Davidson basis
    V[k+1]->zeroEntries();
    for ( int i = 0; i < niters; i++ ){
      V[k+1]->axpy(-res[i], Z[i]);
    }
  
    // Compute the product to test the error
    // (1 - P*Q^{T})*(A - theta*B)*(1 - Q*P^{T})
    W[0]->copyValues(V[k+1]);
    for ( int j = 0; j <= nconverged; j++ ){
      TacsScalar h = W[0]->dot(P[j]);
      W[0]->axpy(-h, Q[j]);
    }

    iteration++;
  }

  if (ksm_print){
    for ( int i = 0; i < nconverged; i++ ){
      char line[256];
      sprintf(line, "Eigenvalue[%2d]: %25.10e\n",
              i, TacsRealPart(eigvals[i]));
      ksm_print->print(line);
    }
  }
  // Print the iteration count to file
  if (ksm_file){
    char line[256];
    sprintf(line, "%2d\n", iteration);
    ksm_file->print(line);
  }
  delete [] rwork;
}

/*
  Set the relative and absolute tolerances used for the stopping
  criterion.

  input:
  rtol: the relative tolerance ||r_k|| < rtol*||r_0||
  atol: the absolute tolerancne ||r_k|| < atol
*/
void TACSJacobiDavidson::setTolerances( double _eigtol,
                                        double _rtol, 
                                        double _atol ){
  eigtol = _eigtol;
  rtol = _rtol;
  atol = _atol;
}
/*
  Set the number of vectors to recycle if the eigenvectors are converged 

  input:
  recycle: number of vectors to recycle
*/
void TACSJacobiDavidson::setRecycle( int _recycle, 
                                     JDRecycleType _recycle_type ){ 
  recycle = _recycle;
  recycle_type = _recycle_type;  
}
