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

#include "TACSRigidBody.h"

#include "TACSElementAlgebra.h"
#include "TACSElementQuaternion.h"

/*
  Rigid-body dynamics routines for TACS
*/

/*
  A generic constructor for the rigid body where the user can directly
  supply the information about the geometry
*/
TACSRigidBodyViz::TACSRigidBodyViz(int _npts, int _nelems, TacsScalar *_Xpt,
                                   int _conn[], TACSGibbsVector *vref) {
  init(_npts, _nelems, _Xpt, _conn, vref);
}

/*
  Create a cube for visualization
*/
TACSRigidBodyViz::TACSRigidBodyViz(TacsScalar Lx, TacsScalar Ly,
                                   TacsScalar Lz) {
  // Create an array with the cube
  int index[27];
  int elem_conn[9 * 6];
  TacsScalar X[3 * 26];

  // Set the point locations
  int count = 0;
  for (int kk = 0; kk < 3; kk++) {
    for (int jj = 0; jj < 3; jj++) {
      for (int ii = 0; ii < 3; ii++) {
        if (ii == 1 && jj == 1 && kk == 1) {
          index[ii + 3 * jj + 9 * kk] = -1;
        } else {
          X[3 * count] = (0.5 * ii - 0.5) * Lx;
          X[3 * count + 1] = (0.5 * jj - 0.5) * Ly;
          X[3 * count + 2] = (0.5 * kk - 0.5) * Lz;
          index[ii + 3 * jj + 9 * kk] = count;
          count++;
        }
      }
    }
  }

  // Create the ii-surfaces
  int elem = 0;
  for (int ii = 0; ii < 3; ii += 2) {
    for (int kk = 0; kk < 3; kk++) {
      for (int jj = 0; jj < 3; jj++) {
        elem_conn[9 * elem + jj + 3 * kk] = index[ii + 3 * jj + 9 * kk];
      }
    }
    elem++;
  }

  // Create the jj-surfaces
  for (int jj = 0; jj < 3; jj += 2) {
    for (int kk = 0; kk < 3; kk++) {
      for (int ii = 0; ii < 3; ii++) {
        elem_conn[9 * elem + ii + 3 * kk] = index[ii + 3 * jj + 9 * kk];
      }
    }
    elem++;
  }

  // Create the kk-surfaces
  for (int kk = 0; kk < 3; kk += 2) {
    for (int jj = 0; jj < 3; jj++) {
      for (int ii = 0; ii < 3; ii++) {
        elem_conn[9 * elem + ii + 3 * jj] = index[ii + 3 * jj + 9 * kk];
      }
    }
    elem++;
  }

  // Initialize the mesh
  init(26, 6, X, elem_conn, NULL);
}

/*
  Initialize the rigid body mesh
*/
void TACSRigidBodyViz::init(int _npts, int _nelems, TacsScalar *_Xpt,
                            int _conn[], TACSGibbsVector *vref) {
  // Copy over the inputs
  npts = _npts;
  nelems = _nelems;

  conn = new int[9 * nelems];
  memcpy(conn, _conn, 9 * nelems * sizeof(int));

  Xpts = new TacsScalar[npts * 3];
  if (vref) {
    const TacsScalar *vorig;
    vref->getVector(&vorig);
    for (int i = 0; i < npts; i++) {
      Xpts[i * 3 + 0] = _Xpt[i * 3 + 0] - vorig[0];
      Xpts[i * 3 + 1] = _Xpt[i * 3 + 1] - vorig[1];
      Xpts[i * 3 + 2] = _Xpt[i * 3 + 2] - vorig[2];
    }
  } else {
    memcpy(Xpts, _Xpt, 3 * npts * sizeof(TacsScalar));
  }
}

/*
  Destructor
*/
TACSRigidBodyViz::~TACSRigidBodyViz() {
  if (Xpts) {
    delete[] Xpts;
  }
  if (conn) {
    delete[] conn;
  }
}

/*
  Get the mesh for the rigid body
*/
void TACSRigidBodyViz::getMesh(int *_npts, int *_nelems,
                               const TacsScalar **_Xpts, const int **_conn) {
  if (_npts) {
    *_npts = npts;
  }
  if (_nelems) {
    *_nelems = nelems;
  }
  if (_Xpts) {
    *_Xpts = Xpts;
  }
  if (_conn) {
    *_conn = conn;
  }
}

/*
  Write the relative error for components of a vector for a
  finite-difference or complex-step study to a given file

  input:
  fp:         the output file
  descript:   description of the component
  a:          analytic values
  b:          set of finite-difference/complex-step values
  size:       the number of values
  rel_err:    relative error tolerance
*/
static void writeErrorComponents(FILE *fp, const char *descript, TacsScalar *a,
                                 TacsScalar *fd, int size,
                                 double rel_err = 1e-12) {
  int print_flag = 1;
  for (int i = 0; i < size; i++) {
    double rel = 0.0;
    if (a[i] != 0.0) {
      rel = fabs(TacsRealPart((a[i] - fd[i]) / a[i]));
    } else {
      rel = fabs(TacsRealPart((a[i] - fd[i])));
    }

    if (rel > rel_err || a[i] != a[i] || fd[i] != fd[i]) {
      if (print_flag) {
        fprintf(fp, "%*s[   ] %15s %15s %15s\n", (int)strlen(descript), "Val",
                "Analytic", "Approximate", "Rel. Error");
        print_flag = 0;
      }
      if (a[i] != 0.0) {
        fprintf(fp, "%s[%3d] %15.6e %15.6e %15.4e\n", descript, i,
                TacsRealPart(a[i]), TacsRealPart(fd[i]),
                fabs(TacsRealPart((a[i] - fd[i]) / a[i])));
      } else {
        fprintf(fp, "%s[%3d] %15.6e %15.6e\n", descript, i, TacsRealPart(a[i]),
                TacsRealPart(fd[i]));
      }
    }
  }
}

/*
  A reference coordinate frame

  This code generates a reference frame from three vectors. The
  vectors define a primary direction and a secondary direction which
  are used to three orthonormal/right-handed vectors. These position
  vectors must be specified with reference to a global reference
  frame.

  input:
  r0:    the base point of the coordinate frame
  r1:    (r1 - r0)/||r1 - r0|| forms the first direction
  r2:    direction perpendicular to (r1 - r0) used as second direction
*/
TACSRefFrame::TACSRefFrame(TACSGibbsVector *_r0, TACSGibbsVector *_r1,
                           TACSGibbsVector *_r2) {
  r0 = _r0;
  r1 = _r1;
  r2 = _r2;
  r0->incref();
  r1->incref();
  r2->incref();
  initialize();
}

/*
  Free the allocated memory
*/
TACSRefFrame::~TACSRefFrame() {
  r0->decref();
  r1->decref();
  r2->decref();
}

/*
  Retrieve the rotation matrix from the global to local coordinates
*/
void TACSRefFrame::getRotation(const TacsScalar **_C) { *_C = C; }

/*
  Update the coordinate frame to reflect any changes to the vectors
  that form the basis of the initial transformation matrix.

  This should be called once during the initial condition
  calculation. This will enable the modification of reference frames
  based solely on moving points.
*/
void TACSRefFrame::initialize() {
  // Get the counts of the design variables
  int elemIndex = 0;
  r0offset = 0;
  r1offset = r0->getDesignVarNums(elemIndex, 0, NULL);
  r2offset = r1offset + r1->getDesignVarNums(elemIndex, 0, NULL);

  // Compute the initial position based on v1
  const TacsScalar *x0, *x1, *x2;
  r0->getVector(&x0);
  r1->getVector(&x1);
  r2->getVector(&x2);

  // Compute the first row of the transformation matrix
  TacsScalar d1[3];
  d1[0] = x1[0] - x0[0];
  d1[1] = x1[1] - x0[1];
  d1[2] = x1[2] - x0[2];

  // Compute the first row of the transformation matrix
  TacsScalar C1norm = sqrt(vec3Dot(d1, d1));
  TacsScalar invC1norm = 1.0 / C1norm;
  C[0] = invC1norm * d1[0];
  C[1] = invC1norm * d1[1];
  C[2] = invC1norm * d1[2];

  // Compute the second vector for the off-axis direction
  TacsScalar d2[3];
  d2[0] = x2[0] - x0[0];
  d2[1] = x2[1] - x0[1];
  d2[2] = x2[2] - x0[2];

  // Orthogonalize the first and second vectors
  TacsScalar dot = vec3Dot(&C[0], d2);
  TacsScalar s2[3];
  s2[0] = d2[0] - dot * C[0];
  s2[1] = d2[1] - dot * C[1];
  s2[2] = d2[2] - dot * C[2];

  // Normalize the vector
  TacsScalar C2norm = sqrt(vec3Dot(s2, s2));
  TacsScalar invC2norm = 1.0 / C2norm;
  C[3] = invC2norm * s2[0];
  C[4] = invC2norm * s2[1];
  C[5] = invC2norm * s2[2];

  // Compute the third basis vector (stored as the third row in
  // rotation matrix). We use compute handed dextral set of basis
  // vectors.
  crossProduct(1.0, &C[0], &C[3], &C[6]);

  //-----------------------------------------------------------------//
  // Code to compute the derivative
  //-----------------------------------------------------------------//

  // The derivative of C1 w.r.t. d1 = (x1 - x0)
  vec3NormDeriv(C1norm, d1, dC1d1);

  // Compute the derivative of C2 w.r.t. d2
  TacsScalar dC2ds2[9];
  vec3NormDeriv(C2norm, s2, dC2ds2);

  TacsScalar ds2d2[9];
  ds2d2[0] = 1.0 - C[0] * C[0];
  ds2d2[1] = -C[0] * C[1];
  ds2d2[2] = -C[0] * C[2];

  ds2d2[3] = -C[1] * C[0];
  ds2d2[4] = 1.0 - C[1] * C[1];
  ds2d2[5] = -C[1] * C[2];

  ds2d2[6] = -C[2] * C[0];
  ds2d2[7] = -C[2] * C[1];
  ds2d2[8] = 1.0 - C[2] * C[2];

  // Compute dC2d2 = dC2ds2*ds2d2
  mat3x3MatMult(dC2ds2, ds2d2, dC2d2);

  // Compute the derivative dC2d1 = dC2ds2*ds2dC1*dC1d1
  // First find the derivative of s2 w.r.t. C1
  TacsScalar ds2dC1[9];
  ds2dC1[0] = -(dot + C[0] * d2[0]);
  ds2dC1[1] = -C[0] * d2[1];
  ds2dC1[2] = -C[0] * d2[2];

  ds2dC1[3] = -C[1] * d2[0];
  ds2dC1[4] = -(dot + C[1] * d2[1]);
  ds2dC1[5] = -C[1] * d2[2];

  ds2dC1[6] = -C[2] * d2[0];
  ds2dC1[7] = -C[2] * d2[1];
  ds2dC1[8] = -(dot + C[2] * d2[2]);

  // Compute the product to find the derivative of C2 w.r.t. d1
  TacsScalar tmp[9];
  mat3x3MatMult(dC2ds2, ds2dC1, tmp);
  mat3x3MatMult(tmp, dC1d1, dC2d1);
}

/*
  Get the design variables numbers from the relevant vectors
*/
int TACSRefFrame::getDesignVarNums(int elemIndex, int dvLen, int dvNums[]) {
  int count = 0;
  count += r0->getDesignVarNums(elemIndex, dvLen, &dvNums[count]);
  count += r1->getDesignVarNums(elemIndex, dvLen, &dvNums[count]);
  count += r2->getDesignVarNums(elemIndex, dvLen, &dvNums[count]);

  return count;
}

/*
  Set the design variables into the three points which define the
  reference frame and recompute the rotation matrix with new set of
  design varaibles

  input:
  dvs:     the design variable values
  numDVs:  the number of design vars/array length
*/
int TACSRefFrame::setDesignVars(int elemIndex, int dvLen,
                                const TacsScalar dvs[]) {
  // Set the design varibles
  int count = 0;
  count += r0->setDesignVars(elemIndex, dvLen, &dvs[count]);
  count += r1->setDesignVars(elemIndex, dvLen, &dvs[count]);
  count += r2->setDesignVars(elemIndex, dvLen, &dvs[count]);

  // Recompute the rotation matrix and other variables
  initialize();

  return count;
}

/*
  Get the design variable values from the points that define the
  reference frame.

  input:
  numDVs:  the number of design vars/array length

  output:
  dvs:     the design variable values
*/
int TACSRefFrame::getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]) {
  // Set the design varibles
  int count = 0;
  count += r0->getDesignVars(elemIndex, dvLen, &dvs[count]);
  count += r1->getDesignVars(elemIndex, dvLen, &dvs[count]);
  count += r2->getDesignVars(elemIndex, dvLen, &dvs[count]);

  return count;
}

/*
  Add the derivative of the inner product of the rotation matrix with
  two input vectors to an array.

  The code computes the following:

  fdvSens = d(psi^{T}*C*phi)/dx

  The code adds the result to the design variable locations of the r0,
  r1, and r2 points. The derivatives of the rotation matrix with
  respect to the points are computed when the point locations are
  updated. These derivatives are stored in dC1d1, dC2d1, and dC2d2.
  These represent the derivatives of the first and second rows of the
  rotation matrix with respect to the difference between (r1 - r0) and
  (r2 - r0), respectively.

  input:
  dvLen:   the number of design variables in the array
  psi:      the pre-multiplying vector
  phi:      the post-multiplying vector

  output:
  fdvSens:  the derivatives are added to this array
*/
void TACSRefFrame::addRotationAdjResProduct(const TacsScalar psi[],
                                            const TacsScalar phi[], int dvLen,
                                            TacsScalar fdvSens[]) {
  // Temporary vector
  TacsScalar tmp[3];

  // Components that store the contributions to each point
  TacsScalar dx0[3], dx1[3], dx2[3];

  // Compute the contribution from the derivative of the first row
  mat3x3MultTrans(dC1d1, phi, dx1);
  vec3Scale(psi[0], dx1);

  // Compute the contributions from the derviative of the second row
  mat3x3MultTrans(dC2d1, phi, tmp);
  vec3Axpy(psi[1], tmp, dx1);

  mat3x3MultTrans(dC2d2, phi, dx2);
  vec3Scale(psi[1], dx2);

  // Compute the contribution from the third row such that c3 = c1 x c2
  crossProduct(-psi[2], &C[0], phi, tmp);
  mat3x3MultTransAdd(dC2d1, tmp, dx1);
  mat3x3MultTransAdd(dC2d2, tmp, dx2);

  crossProduct(psi[2], &C[3], phi, tmp);
  mat3x3MultTransAdd(dC1d1, tmp, dx1);

  // Add the derivative to the sensitivity vector
  r1->addPointAdjResProduct(1.0, dx1, dvLen, &fdvSens[r1offset]);
  r2->addPointAdjResProduct(1.0, dx2, dvLen, &fdvSens[r2offset]);

  // Add the result from the base point
  dx0[0] = -(dx1[0] + dx2[0]);
  dx0[1] = -(dx1[1] + dx2[1]);
  dx0[2] = -(dx1[2] + dx2[2]);
  r0->addPointAdjResProduct(1.0, dx0, dvLen, &fdvSens[r0offset]);
}

/*
  Test the implementation of the derivative of the rotation matrix

  input:
  numDVs: the number of design variables
  dh:     the step size for the central-difference formula
*/
void TACSRefFrame::testRotation(double dh) {
  // Allocate design variable arrays
  int elemIndex = 0;
  int numDVs = getDesignVarNums(0, 0, NULL);
  TacsScalar *x = new TacsScalar[numDVs];
  TacsScalar *product = new TacsScalar[numDVs];

  // Set random values into the psi/phi vectors
  TacsScalar psi[3], phi[3];
  for (int k = 0; k < 3; k++) {
    psi[k] = -0.5 + (1.0 * rand()) / RAND_MAX;
    phi[k] = -0.5 + (1.0 * rand()) / RAND_MAX;
  }

  // Get the variables and make sure that they are consistent
  getDesignVars(elemIndex, numDVs, x);
  setDesignVars(elemIndex, numDVs, x);

  // Compute the derivative and store it in fdvSens
  memset(product, 0, numDVs * sizeof(TacsScalar));
  addRotationAdjResProduct(psi, phi, numDVs, product);

  // Set the vectors to zero
  TacsScalar fd[9];
  memset(fd, 0, sizeof(fd));

  // Set the vectors
  for (int i = 0; i < numDVs; i++) {
    TacsScalar t[3];
    TacsScalar xtmp = x[i];

#ifdef TACS_USE_COMPLEX
    // Evaluate the matrix at x + j*dh
    x[i] = xtmp + TacsScalar(0.0, dh);
    setDesignVars(elemIndex, numDVs, x);
    mat3x3Mult(C, phi, t);
    fd[i] = TacsImagPart(vec3Dot(psi, t)) / dh;
#else
    // Evaluate C at (x + dh)
    x[i] = xtmp + dh;
    setDesignVars(elemIndex, numDVs, x);
    mat3x3Mult(C, phi, t);
    TacsScalar f1 = vec3Dot(psi, t);

    // Evaluate C at (x - dh)
    x[i] = xtmp - dh;
    setDesignVars(elemIndex, numDVs, x);
    mat3x3Mult(C, phi, t);
    TacsScalar f2 = vec3Dot(psi, t);

    fd[i] = 0.5 * (f1 - f2) / dh;
#endif  // TACS_USE_COMPLEX
    x[i] = xtmp;
  }

  // Set the design variable values back to their original values
  setDesignVars(elemIndex, numDVs, x);

  writeErrorComponents(stdout, "TACSRefFrame", product, fd, 9);

  // Free the allocated memory
  delete[] x;
  delete[] product;
}

/*
  The constructor for the rigid body dynamics

  input:
  mass:       the mass of the body
  c:          the first moment of inertia
  J:          the symmetric second moment of inertia
  g:          the acceleration due to gravity in the global frame
*/
TACSRigidBody::TACSRigidBody(TACSRefFrame *_CRef, const TacsScalar _mass,
                             const TacsScalar _cRef[], const TacsScalar _JRef[],
                             TACSGibbsVector *_rInit, TACSGibbsVector *_vInit,
                             TACSGibbsVector *_omegaInit,
                             TACSGibbsVector *_gvec) {
  // Copy over the property values
  mass = _mass;

  // Copy over the inertial properties in the Ref reference frame
  for (int k = 0; k < 3; k++) {
    cRef[k] = _cRef[k];
  }
  for (int k = 0; k < 6; k++) {
    JRef[k] = _JRef[k];
  }

  // Copy over the reference frame
  CRef = _CRef;

  // Copy over the initial vectors. Note that these vectors
  // are in the global reference frame.
  gvec = _gvec;
  rInit = _rInit;
  vInit = _vInit;
  omegaInit = _omegaInit;

  // Increment the reference counts for things
  CRef->incref();
  gvec->incref();
  rInit->incref();
  vInit->incref();
  omegaInit->incref();

  // Initialize the design variable numbers for the inertial properties
  massDV = -1;
  cDV[0] = cDV[1] = cDV[2] = -1;
  JDV[0] = JDV[1] = JDV[2] = JDV[3] = JDV[4] = JDV[5] = -1;

  viz = NULL;

  // Update the inertial properties
  updateInertialProperties();
}

/*
  Decrease the reference counts to everything
*/
TACSRigidBody::~TACSRigidBody() {
  CRef->decref();
  gvec->decref();
  rInit->decref();
  vInit->decref();
  omegaInit->decref();
  if (viz) {
    viz->decref();
  }
}

// Set the element name
const char *TACSRigidBody::elem_name = "TACSRigidBody";

int TACSRigidBody::getVarsPerNode() { return 8; }

int TACSRigidBody::getNumNodes() { return 1; }

ElementLayout TACSRigidBody::getLayoutType() { return TACS_POINT_ELEMENT; }

/*
  Set the design variable numbers associated with the inertial
  properties of the body

  input:
  massDV:   design variable number for the mass
  cDV:      design variable numbers for the first moment of mass
  JDV:      design variable numbers for the second moments of mass
*/
void TACSRigidBody::setDesignVarNums(int _massDV, const int _cDV[],
                                     const int _JDV[]) {
  massDV = _massDV;
  if (_cDV) {
    cDV[0] = _cDV[0];
    cDV[1] = _cDV[1];
    cDV[2] = _cDV[2];
  }
  if (_JDV) {
    JDV[0] = _JDV[0];
    JDV[1] = _JDV[1];
    JDV[2] = _JDV[2];
    JDV[3] = _JDV[3];
    JDV[4] = _JDV[4];
    JDV[5] = _JDV[5];
  }
}

/*
  Set the design variable values
*/
int TACSRigidBody::setDesignVars(int elemIndex, int dvLen,
                                 const TacsScalar dvs[]) {
  // Set the mass design variable
  int count = 0;
  if (massDV >= 0 && massDV < dvLen) {
    mass = dvs[count];
    count++;
  }

  // Set the moment of mass variable
  for (int k = 0; k < 3; k++) {
    if (cDV[k] >= 0 && cDV[k] < dvLen) {
      cRef[k] = dvs[count];
      count++;
    }
  }

  // Set the second moment of mass variables
  for (int k = 0; k < 6; k++) {
    if (JDV[k] >= 0 && JDV[k] < dvLen) {
      JRef[k] = dvs[count];
      count++;
    }
  }

  // Set the reference frame design variables
  count += CRef->setDesignVars(elemIndex, dvLen, &dvs[count]);

  // Set the design variable values for the initial vectors
  count += gvec->setDesignVars(elemIndex, dvLen, &dvs[count]);
  count += rInit->setDesignVars(elemIndex, dvLen, &dvs[count]);
  count += vInit->setDesignVars(elemIndex, dvLen, &dvs[count]);
  count += omegaInit->setDesignVars(elemIndex, dvLen, &dvs[count]);

  // Update the inertial properties based on the design variable
  // values
  updateInertialProperties();

  return count;
}

/*
  Retrieve the design variable numbers
*/
int TACSRigidBody::getDesignVarNums(int elemIndex, int dvLen, int *dvNums) {
  // Get the mass design variable
  int count = 0;
  if (massDV >= 0 && massDV < dvLen) {
    dvNums[count] = massDV;
    count++;
  }

  // Get the moment of mass variables
  for (int k = 0; k < 3; k++) {
    if (cDV[k] >= 0 && cDV[k] < dvLen) {
      dvNums[count] = cDV[k];
      count++;
    }
  }

  // Get the second moment of mass variables
  for (int k = 0; k < 6; k++) {
    if (JDV[k] >= 0 && JDV[k] < dvLen) {
      dvNums[count] = JDV[k];
      count++;
    }
  }

  // Get the reference frame design variables
  count += CRef->getDesignVarNums(elemIndex, dvLen, &dvNums[count]);

  // Get the design variable values for the initial vectors
  count += gvec->getDesignVarNums(elemIndex, dvLen, &dvNums[count]);
  count += rInit->getDesignVarNums(elemIndex, dvLen, &dvNums[count]);
  count += vInit->getDesignVarNums(elemIndex, dvLen, &dvNums[count]);
  count += omegaInit->getDesignVarNums(elemIndex, dvLen, &dvNums[count]);

  return count;
}

/*
  Retrieve the design variable values
*/
int TACSRigidBody::getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]) {
  // Get the mass design variable
  int count = 0;
  if (massDV >= 0 && massDV < dvLen) {
    dvs[count] = mass;
    count++;
  }

  // Get the moment of mass variables
  for (int k = 0; k < 3; k++) {
    if (cDV[k] >= 0 && cDV[k] < dvLen) {
      dvs[count] = cRef[k];
      count++;
    }
  }

  // Get the second moment of mass variables
  for (int k = 0; k < 6; k++) {
    if (JDV[k] >= 0 && JDV[k] < dvLen) {
      dvs[count] = JRef[k];
      count++;
    }
  }

  // Get the reference frame design variables
  count += CRef->getDesignVars(elemIndex, dvLen, &dvs[count]);

  // Get the design variable values for the initial vectors
  count += gvec->getDesignVars(elemIndex, dvLen, &dvs[count]);
  count += rInit->getDesignVars(elemIndex, dvLen, &dvs[count]);
  count += vInit->getDesignVars(elemIndex, dvLen, &dvs[count]);
  count += omegaInit->getDesignVars(elemIndex, dvLen, &dvs[count]);

  return count;
}

/*
  Retrieve the design variable range
*/
int TACSRigidBody::getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                                     TacsScalar ub[]) {
  return 0;
}

/*
  Set the inertial properties in the global frame based on the
  inertial properties in the reference frame
*/
void TACSRigidBody::updateInertialProperties() {
  const TacsScalar *C;
  CRef->getRotation(&C);

  // Convert the first moment of inertial from the local to the
  // inertial reference frame
  mat3x3MultTrans(C, cRef, c);

  // Copy the J values to a row
  TacsScalar Jtmp[9], CJtmp[9];

  // Compute CJtmp = C^{T}*JRef;
  Jtmp[0] = JRef[0];
  Jtmp[1] = JRef[1];
  Jtmp[2] = JRef[2];
  mat3x3MultTrans(C, Jtmp, &CJtmp[0]);

  Jtmp[0] = JRef[1];
  Jtmp[1] = JRef[3];
  Jtmp[2] = JRef[4];
  mat3x3MultTrans(C, Jtmp, &CJtmp[3]);

  Jtmp[0] = JRef[2];
  Jtmp[1] = JRef[4];
  Jtmp[2] = JRef[5];
  mat3x3MultTrans(C, Jtmp, &CJtmp[6]);

  // Compute Jtmp = C^{T}*[C^{T}*J]^{T} = C^{T}*[CJtmp]^{T}. Note that
  // the matrix CJtmp is stored in column-major order so this
  // multiplication is in fact C^{T}*[CJtmp]^{T}
  mat3x3TransMatMult(C, CJtmp, Jtmp);

  // Copy the symmetric values from the computation
  J[0] = Jtmp[0];
  J[1] = Jtmp[1];
  J[2] = Jtmp[2];

  J[3] = Jtmp[4];
  J[4] = Jtmp[5];

  J[5] = Jtmp[8];
}

/*
  Retrieve the initial values if the kinematic and dynamic variables

  output:
  qkin:  the kinematic variables
  qdyn:  the dynamic variable values
*/
void TACSRigidBody::getInitConditions(int elemIndex, const TacsScalar X[],
                                      TacsScalar vars[], TacsScalar dvars[],
                                      TacsScalar ddvars[]) {
  // Set everything to zero first
  memset(vars, 0, 8 * sizeof(TacsScalar));
  memset(dvars, 0, 8 * sizeof(TacsScalar));

  // Set eta = 1.0
  vars[3] = 1.0;

  // Get the initial position
  const TacsScalar *r;
  rInit->getVector(&r);

  // Get the initial velocity
  const TacsScalar *v;
  vInit->getVector(&v);

  // Get the initial angular velocity
  const TacsScalar *w;
  omegaInit->getVector(&w);

  // Set the initial velocity
  TacsScalar t[3];
  crossProduct(1.0, w, r, t);
  dvars[0] = v[0] + t[0];
  dvars[1] = v[1] + t[1];
  dvars[2] = v[2] + t[2];

  // Set the time-derivative of the quaternions
  dvars[4] = 0.5 * w[0];
  dvars[5] = 0.5 * w[1];
  dvars[6] = 0.5 * w[2];

  // Set the accelerations omega^{x}*omega^{x}*r
  crossProduct(1.0, w, t, &ddvars[0]);
}

/*
  Retrieve the position of the rigid body
*/
TACSGibbsVector *TACSRigidBody::getInitPosition() { return rInit; }

/*
  Compute the kinematic and potential energies of the rigid body

  The kinetic energy is given by:

  T = 0.5*m*dot{r}^{T}*dot{r} + 0.5*omega^{T}*J*omega
  .   + dot{r}^{T}*C^{T}*omega^{x}*c

  where r0 and v0 are expressed in the global reference frame, and C
  is the rotation matrix from the global to body-fixed frames. The
  potential energy is due to the force of gravity

  U = -m*g^{T}*r0 - g^{T}*C^{T}*c

  input:
  time:   the simulation time
  X:      the reference node location
  vars:   the state vector
  dvars:  the time derivative of the states

  output:
  Te:   the kinematic energy
  Pe:   the potential energy
*/
void TACSRigidBody::computeEnergies(int elemIndex, double time,
                                    const TacsScalar X[],
                                    const TacsScalar vars[],
                                    const TacsScalar dvars[], TacsScalar *Te,
                                    TacsScalar *Pe) {
  // Get the acceleration due to gravity
  const TacsScalar *g;
  gvec->getVector(&g);

  // Set the location
  const TacsScalar *r0 = &vars[0];
  const TacsScalar *v0 = &dvars[0];

  // Set the pointers to the Euler parameters
  TacsScalar eta = vars[3];
  TacsScalar deta = dvars[3];
  const TacsScalar *eps = &vars[4];
  const TacsScalar *deps = &dvars[4];

  // Compute the rotation matrix
  TacsScalar C[9];
  computeRotationMat(eta, eps, C);

  // Compute the angular velocity from the Euler parameters
  TacsScalar omega[3];
  computeSRateProduct(eta, eps, deta, deps, omega);

  // Add the contribution from the linear motion
  *Te = 0.5 * mass * vec3Dot(v0, v0);

  // Add the contribution from the angular velocity
  TacsScalar tmp[3];
  mat3x3SymmMult(J, omega, tmp);
  *Te += 0.5 * vec3Dot(omega, tmp);

  // Transform the velocity from the inertial to body-fixed frame
  TacsScalar v[3];
  mat3x3Mult(C, v0, v);

  // Add the coupled contribution from the angular velocity/rotation
  crossProduct(1.0, omega, c, tmp);
  *Te += vec3Dot(v, tmp);

  // Compute the potential energy
  mat3x3Mult(C, g, tmp);
  *Pe = -(mass * vec3Dot(r0, g) + vec3Dot(c, tmp));
}

/*
  Compute the residual of the governing equations for the rigid body
  motion.

  The equations of motion for the rigid body are derived using
  Lagrange's equations of motion with constraints. Within each body,
  the constraint enforcing the norm on the quaternion is imposed
  directly. The motion equations are written as follows:

  m*ddot{r} + C^{T}*(omega^{x}*omega^{x}*c - c^{x}*dot{omega}) - mass*g = 0

  S^{T}*dot{y} + 2*dot{S}^{T}*y
  + D(dot{r})^{T}*c^{x}*omega - D(g)^{T}*c - A^{T}*lambda = 0

  where y = J*omega + c^{x}*C*r and D(v) = d(C*v)/dq

  Note that the matrix S relates the quaternion rates to the angular
  acceleration such that omega = S*dot{q}. The matrix S is given by:

  S(q) = 2[ -eps | (eta*I - eps^{x}) ]

  Note that the matrix S has the property that dot{S} = S(dot{q}).
  The transpose of S is given as follows:

  S^{T} = [      2*eps^{T}      ]
  .       [ 2*(eta*I + eps^{x}) ]

  input:
  time:    the simulation time
  X:       the nodal locations
  vars:    the variables
  dvars:   the first time derivative of the variables
  ddvars:  the second time derivative of the variables

  output:
  res:     the residual of the governing equations
*/
void TACSRigidBody::addResidual(int elemIndex, double time,
                                const TacsScalar *Xpts, const TacsScalar *vars,
                                const TacsScalar *dvars,
                                const TacsScalar *ddvars, TacsScalar *res) {
  // Get the acceleration due to gravity
  const TacsScalar *g;
  gvec->getVector(&g);

  // Set the location and its time derivatives
  const TacsScalar *v0 = &dvars[0];
  const TacsScalar *a0 = &ddvars[0];

  // Set the pointers to the Euler parameters and all their time
  // derivatives
  TacsScalar eta = vars[3];
  TacsScalar deta = dvars[3];
  TacsScalar ddeta = ddvars[3];
  const TacsScalar *eps = &vars[4];
  const TacsScalar *deps = &dvars[4];
  const TacsScalar *ddeps = &ddvars[4];

  // Compute the rotation matrix
  TacsScalar C[9], dotC[9];
  computeRotationMat(eta, eps, C);
  computeRotationMatDeriv(eta, eps, deta, deps, dotC);

  // Compute the angular velocity and acceleration from the Euler
  // parameters
  TacsScalar omega[3], domega[3];
  computeSRateProduct(eta, eps, deta, deps, omega);
  computeSRateProduct(eta, eps, ddeta, ddeps, domega);

  // Compute C^{T}*c^{x}*domega
  TacsScalar t1[3], t2[3];
  crossProduct(1.0, c, domega, t2);
  mat3x3MultTrans(C, t2, t1);

  // Add dot{C}^{T}*c^{x}*omega
  crossProduct(1.0, c, omega, t2);
  mat3x3MultTransAdd(dotC, t2, t1);

  // Complete the governing equations for the translational
  // degrees of freedom
  res[0] += mass * (a0[0] - g[0]) - t1[0];
  res[1] += mass * (a0[1] - g[1]) - t1[1];
  res[2] += mass * (a0[2] - g[2]) - t1[2];

  // Compute the residual of the governing equations for
  // the rotational terms
  // ---------------------------------------------------
  // Compute t1 = C*ddot{r} + dot{C}*dot{r}
  mat3x3Mult(C, a0, t1);
  mat3x3MultAdd(dotC, v0, t1);

  // Compute t2 = c^{x}*(C*ddot{r} + dot{C}*dot{r})
  crossProduct(1.0, c, t1, t2);

  // Add t2 += J*domega
  mat3x3SymmMultAdd(J, domega, t2);

  // Add res += S^{T}*t2
  addSRateTransProduct(1.0, eta, eps, t2, &res[3], &res[4]);

  // Compute t2 = J*omega + c^{x}*C*dot{r}
  mat3x3Mult(C, v0, t1);
  crossProduct(1.0, c, t1, t2);
  mat3x3SymmMultAdd(J, omega, t2);

  // Add res[3:] += 2.0*dot{S}^{T}*t2
  addSRateTransProduct(2.0, deta, deps, t2, &res[3], &res[4]);

  // Add res += D(dot(r))^{T}*c^{x}*omega
  crossProduct(1.0, c, omega, t1);
  addDMatTransProduct(1.0, v0, t1, eta, eps, &res[3], &res[4]);

  // Add res -= D(g)^{T}*c
  addDMatTransProduct(-1.0, g, c, eta, eps, &res[3], &res[4]);

  // Add the Lagrange multiplier term
  res[3] += 2.0 * eta * vars[7];
  res[4] += 2.0 * eps[0] * vars[7];
  res[5] += 2.0 * eps[1] * vars[7];
  res[6] += 2.0 * eps[2] * vars[7];

  // Compute the quaternion constraint
  res[7] += eta * eta + vec3Dot(eps, eps) - 1.0;
}

/*
  Compute the Jacobian of the governing equations
*/
void TACSRigidBody::addJacobian(int elemIndex, double time, TacsScalar alpha,
                                TacsScalar beta, TacsScalar gamma,
                                const TacsScalar *Xpts, const TacsScalar *vars,
                                const TacsScalar *dvars,
                                const TacsScalar *ddvars, TacsScalar *res,
                                TacsScalar *mat) {
  addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);

  // Get the acceleration due to gravity
  const TacsScalar *g;
  gvec->getVector(&g);

  // Set the location and its time derivatives
  const TacsScalar *v0 = &dvars[0];
  const TacsScalar *a0 = &ddvars[0];

  // Set the pointers to the Euler parameters and all their time
  // derivatives
  TacsScalar eta = vars[3];
  const TacsScalar *eps = &vars[4];
  TacsScalar deta = dvars[3];
  const TacsScalar *deps = &dvars[4];
  TacsScalar ddeta = ddvars[3];
  const TacsScalar *ddeps = &ddvars[4];

  // Compute the rotation matrix
  TacsScalar C[9], dotC[9];
  computeRotationMat(eta, eps, C);
  computeRotationMatDeriv(eta, eps, deta, deps, dotC);

  // Compute the rotation rate matrices
  TacsScalar S[12], Sdot[12], Sddot[12];
  computeSRateMat(eta, eps, S);
  computeSRateMat(deta, deps, Sdot);
  computeSRateMat(ddeta, ddeps, Sddot);

  // Compute the angular velocity and acceleration from the Euler
  // parameters
  TacsScalar omega[3], domega[3];
  computeSRateProduct(eta, eps, deta, deps, omega);
  computeSRateProduct(eta, eps, ddeta, ddeps, domega);

  // Add the components from the position governing equation
  // -------------------------------------------------------
  mat[0] += gamma * mass;
  mat[9] += gamma * mass;
  mat[18] += gamma * mass;

  // Temporary 3x4 matrices
  TacsScalar A[12], B[12];

  // Add the term C^{T}*c^{x}*S
  setMatSkew(-gamma, c, A);
  mat3x3TransMatMult(C, A, B);
  addBlock3x3x4Product(B, S, &mat[3], 8);

  // Add the term C^{T}*c^{x}*Sddot
  setMatSkew(alpha, c, A);
  mat3x3TransMatMult(C, A, B);
  addBlock3x3x4Product(B, Sddot, &mat[3], 8);

  // Add the term -E(c^{x}*domega)
  // t1 = c^{x}*domega
  TacsScalar t1[3], t2[3];
  crossProduct(1.0, c, domega, t1);
  addBlockEMat(-alpha, eta, eps, t1, &mat[3], 8);

  // Add the term -dot{C}^{T}*c^{x}*S
  setMatSkew(-beta, c, A);
  mat3x3TransMatMult(dotC, A, B);
  addBlock3x3x4Product(B, S, &mat[3], 8);

  // Add the term dot{C}^{T}*c^{x}*Sdot
  setMatSkew(alpha, c, A);
  mat3x3TransMatMult(dotC, A, B);
  addBlock3x3x4Product(B, Sdot, &mat[3], 8);

  // Add the terms from the derivative of the dot{C} matrix
  crossProduct(1.0, c, omega, t1);
  addBlockEMat(-alpha, deta, deps, t1, &mat[3], 8);
  addBlockEMat(-beta, eta, eps, t1, &mat[3], 8);

  // Add the terms from the governing equations for the quaternions
  // --------------------------------------------------------------

  // Add the term S^{T}*(gamma*c^{x}*C + beta*c^{x}*dot{C})
  setMatSkew(gamma, c, A);
  mat3x3MatMult(A, C, B);
  addBlock4x3x3Product(S, B, &mat[24], 8);

  setMatSkew(beta, c, A);
  mat3x3MatMult(A, dotC, B);
  addBlock4x3x3Product(S, B, &mat[24], 8);

  // Add the term  2*beta*dot{S}^{T}*c^{x}*C
  setMatSkew(2.0 * beta, c, A);
  mat3x3MatMult(A, C, B);
  addBlock4x3x3Product(Sdot, B, &mat[24], 8);

  // Add E^{T}(c^{x}*omega)
  crossProduct(beta, c, omega, t1);
  addBlockEMatTrans(1.0, eta, eps, t1, &mat[24], 8);

  // Add terms to the derivative of the quaternion governing
  // equations with respect to the quaternions

  // Add S^{T}*J*S to the matrix
  matSymmMat3x4Mult(J, S, A);
  addBlock3x4Product(gamma, S, A, &mat[27], 8);

  // Add dot{S}^{T}*J*S
  addBlock3x4Product(2.0 * beta, Sdot, A, &mat[27], 8);

  // Add -alpha*S^{T}*J*ddot{S}
  addBlock3x4Product(-alpha, A, Sddot, &mat[27], 8);

  // Add the term -2*alpha*dot{S}^{T}*J*dot{S}
  matSymmMat3x4Mult(J, Sdot, A);
  addBlock3x4Product(-2.0 * alpha, Sdot, A, &mat[27], 8);

  // Compute S^{T}*c^{x}*(C*v0)^{x}*S
  TacsScalar v[3];
  mat3x3Mult(C, v0, v);

  // Compute t1 = J*omega + c^{x}*(C*dot{r})
  mat3x3SymmMult(J, omega, t1);
  crossProductAdd(1.0, c, v, t1);
  addSRateMatTransDeriv(2.0 * beta, t1, &mat[27], 8);

  // Add the term D(dot{r})^{T}*c^{x}*S
  setMatSkew(beta, c, A);
  matMat3x4Mult(A, S, B);
  computeDMat(eta, eps, v0, A);
  addBlock3x4Product(1.0, A, B, &mat[27], 8);

  // Compute t1 = C*ddot{r} + dot{C}*dot{r}
  mat3x3Mult(C, a0, t1);
  mat3x3MultAdd(dotC, v0, t1);

  // Compute t2 = c^{x}*(C*ddot{r} + dot{C}*dot{r})
  crossProduct(1.0, c, t1, t2);

  // Add t2 += J*domega
  mat3x3SymmMultAdd(J, domega, t2);

  // Add the derivative TS(t1) to the Jacobian matrix
  addSRateMatTransDeriv(alpha, t2, &mat[27], 8);

  // Compute the derivatives of C*v0 and C*a0 w.r.t. q
  TacsScalar Dv0[12], Da0[12];
  computeDMat(eta, eps, v0, Dv0);
  computeDMat(eta, eps, a0, Da0);

  // Add the term alpha*S^{T}*c^{x}*D(a0)
  setMatSkew(1.0, c, A);
  matMat3x4Mult(A, Da0, B);
  addBlock3x4Product(alpha, S, B, &mat[27], 8);

  // Compute the term alpha*S^{T}*c^{x}*d(dot{C}*v0)/dq
  TacsScalar dotDv0[12];
  computeDMat(deta, deps, v0, dotDv0);
  setMatSkew(1.0, c, A);
  matMat3x4Mult(A, dotDv0, B);
  addBlock3x4Product(alpha, S, B, &mat[27], 8);

  // Compute the term beta*S^{T}*c^{x}*d(dot{C}*v0)/ddot{q}
  setMatSkew(1.0, c, A);
  matMat3x4Mult(A, Dv0, B);
  addBlock3x4Product(beta, S, B, &mat[27], 8);

  // Add the term 2*dot{S}^{T}*c^{x}*D(dot{r})
  setMatSkew(1.0, c, A);
  matMat3x4Mult(A, Dv0, B);
  addBlock3x4Product(2.0 * alpha, Sdot, B, &mat[27], 8);

  // Add the derivative: d(D(v0)^{T}*t1)/dq with t1 = c^{x}*omega
  crossProduct(1.0, c, omega, t1);
  addBlockDMatTransDeriv(alpha, v0, t1, &mat[27], 8);

  // Add the derivative -alpha*D(v0)^{T}*c^{x}*dot{S}
  setMatSkew(1.0, c, A);
  matMat3x4Mult(A, Sdot, B);
  addBlock3x4Product(-alpha, Dv0, B, &mat[27], 8);

  // Add the derivative: -d(D(g)^{T}*c)/dq
  addBlockDMatTransDeriv(-alpha, g, c, &mat[27], 8);

  // Add the terms from the Lagrange multipliers
  mat[31] += 2.0 * alpha * eta;
  mat[39] += 2.0 * alpha * eps[0];
  mat[47] += 2.0 * alpha * eps[1];
  mat[55] += 2.0 * alpha * eps[2];

  mat[59] += 2.0 * alpha * eta;
  mat[60] += 2.0 * alpha * eps[0];
  mat[61] += 2.0 * alpha * eps[1];
  mat[62] += 2.0 * alpha * eps[2];

  mat[27] += 2.0 * alpha * vars[7];
  mat[36] += 2.0 * alpha * vars[7];
  mat[45] += 2.0 * alpha * vars[7];
  mat[54] += 2.0 * alpha * vars[7];
}

/*
  The following function tests the consistency of the implementation
  of the residuals and the energy expressions, relying on Lagrange's
  equations.

  This function uses finite-differences to compute the derivatives
  within Lagrange's equations and compares the result with the
  residuals of the EOM computed using the residual routine.

  Lagrange's equations of motion are given as follows:

  d/dt(dL/d(dot{q})^{T}) - dL/dq^{T} = 0

  This can be evaluated using finite-differencing as follows:

  dL/dqi(q, dq) .= (L(q, dq + h*ei) - L(q, dq - h*ei))/h

  d(f(q, dq))/dt .=
  (f(q + dt*dq, dq + dt*ddq) - f(q - dt*dq, dq - dt*ddq))/dt
*/
void TACSRigidBody::testResidual(double dh) {
  double time = 0.0;

  // Set the position vector
  TacsScalar X[3] = {1.1, 2.2, 3.3};

  // Set the variable values
  TacsScalar vars[8], dvars[8], ddvars[8];

  // Compute the variable values
  for (int i = 0; i < 7; i++) {
    vars[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    dvars[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    ddvars[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Set the Largange multipliers to zero
  vars[7] = dvars[7] = ddvars[7] = 0.0;

  // Normalize the rotation variables
  TacsScalar e = 1.0 / sqrt(vars[3] * vars[3] + vec3Dot(&vars[4], &vars[4]));
  for (int i = 0; i < 4; i++) {
    vars[3 + i] *= e;
  }

  // Normalize the time derivatives so that they lie within
  // the null space of the quaternion constraint
  dvars[3] = -vec3Dot(&vars[4], &dvars[4]) / vars[3];
  ddvars[3] = -((dvars[3] * dvars[3] + vec3Dot(&dvars[4], &dvars[4])) +
                vec3Dot(&vars[4], &ddvars[4])) /
              vars[3];

  // Temporary vectors containing the derivative
  TacsScalar fd[8], res1[8], res2[8];

  // The values of the variables at the perturbed locations
  TacsScalar q[8], dq[8];

  // Compute the values of the variables at (t + dt)
  for (int i = 0; i < 8; i++) {
    q[i] = vars[i] + dh * dvars[i];
    dq[i] = dvars[i] + dh * ddvars[i];
  }

  // Evaluate the derivative w.r.t. dot{q}
  for (int i = 0; i < 8; i++) {
    // Evaluate the finite-difference for component i
    TacsScalar T1, P1, T2, P2;
    TacsScalar dqtmp = dq[i];
    dq[i] = dqtmp + dh;
    computeEnergies(0, time, X, q, dq, &T1, &P1);

    dq[i] = dqtmp - dh;
    computeEnergies(0, time, X, q, dq, &T2, &P2);

    // Compute and store the approximation
    res1[i] = 0.5 * ((T1 - P1) - (T2 - P2)) / dh;
    dq[i] = dqtmp;
  }

  // Compute the values of the variables at (t - dt)
  for (int i = 0; i < 8; i++) {
    q[i] = vars[i] - dh * dvars[i];
    dq[i] = dvars[i] - dh * ddvars[i];
  }

  // Evaluate the derivative w.r.t. dot{q}
  for (int i = 0; i < 8; i++) {
    // Evaluate the finite-difference for component i
    TacsScalar T1, P1, T2, P2;
    TacsScalar dqtmp = dq[i];
    dq[i] = dqtmp + dh;
    computeEnergies(0, time, X, q, dq, &T1, &P1);

    dq[i] = dqtmp - dh;
    computeEnergies(0, time, X, q, dq, &T2, &P2);

    // Compute and store the approximation
    res2[i] = 0.5 * ((T1 - P1) - (T2 - P2)) / dh;
    dq[i] = dqtmp;
  }

  // Evaluate the finite-difference for the first term in Largrange's
  // equations of motion
  for (int i = 0; i < 8; i++) {
    fd[i] = 0.5 * (res1[i] - res2[i]) / dh;
  }

  // Reset the values of q and dq at time t
  for (int i = 0; i < 8; i++) {
    q[i] = vars[i];
    dq[i] = dvars[i];
  }

  // Compute the contribution from dL/dq^{T}
  for (int i = 0; i < 8; i++) {
    // Evaluate the finite-difference for component i
    TacsScalar T1, P1, T2, P2;
    TacsScalar qtmp = q[i];
    q[i] = qtmp + dh;
    computeEnergies(0, time, X, q, dq, &T1, &P1);

    q[i] = qtmp - dh;
    computeEnergies(0, time, X, q, dq, &T2, &P2);

    // Compute and store the approximation
    res1[i] = 0.5 * ((T1 - P1) - (T2 - P2)) / dh;
    q[i] = qtmp;
  }

  // Add the result to the finite-difference result
  for (int i = 0; i < 8; i++) {
    fd[i] -= res1[i];
  }

  // Evaluate the residual using the code
  memset(res1, 0, 8 * sizeof(TacsScalar));
  addResidual(0, time, X, vars, dvars, ddvars, res1);

  // Write out the error components
  writeErrorComponents(stdout, "Res error", res1, fd, 8);
}

/*
  The following function tests the consistency between the
  implementation of the residuals and the implementation of the system
  Jacobian.

  input:
  dh:      the finite-difference step size
  alpha:   coefficient for the variables
  beta:    coefficient for the time derivative variables
  gamma:   coefficient for the second time derivative variables
*/
void TACSRigidBody::testJacobian(double dh, TacsScalar alpha, TacsScalar beta,
                                 TacsScalar gamma) {
  double time = 0.0;

  // Set the position vector
  TacsScalar X[3] = {0.0, 0.0, 0.0};

  // Set the variable values
  TacsScalar vars[8], dvars[8], ddvars[8];

  // Compute the variable values
  for (int i = 0; i < 8; i++) {
    vars[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    dvars[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    ddvars[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // The computed Jacobian of the element matrix
  TacsScalar mat[64];

  // The finite-difference result
  TacsScalar fd[8], res[8];

  // The perturb direction to test
  TacsScalar perb[8];

  // Temporary variables and their time derivatives
  TacsScalar q[8], dq[8], ddq[8];

  // Set random perburbed values
  for (int i = 0; i < 8; i++) {
    perb[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  for (int ii = 0; ii < 8; ii++) {
    memset(perb, 0, 8 * sizeof(TacsScalar));
    perb[ii] = 1.0;

#ifdef TACS_USE_COMPLEX
    // Set the values for the first evaluation
    TacsScalar a = TacsScalar(0.0, dh * TacsRealPart(alpha));
    TacsScalar b = TacsScalar(0.0, dh * TacsRealPart(beta));
    TacsScalar c = TacsScalar(0.0, dh * TacsRealPart(gamma));
    for (int i = 0; i < 8; i++) {
      q[i] = vars[i] + a * perb[i];
      dq[i] = dvars[i] + b * perb[i];
      ddq[i] = ddvars[i] + c * perb[i];
    }

    // Get the residual at vars + alpha*perb, ... etc.
    memset(fd, 0, 8 * sizeof(TacsScalar));
    addResidual(0, time, X, q, dq, ddq, fd);

    // Form the finite-difference matrix-vector approximation
    for (int i = 0; i < 8; i++) {
      fd[i] = TacsImagPart(fd[i]) / dh;
    }
#else
    // Set the values for the first evaluation
    for (int i = 0; i < 8; i++) {
      q[i] = vars[i] + dh * alpha * perb[i];
      dq[i] = dvars[i] + dh * beta * perb[i];
      ddq[i] = ddvars[i] + dh * gamma * perb[i];
    }

    // Get the residual at vars + alpha*perb, ... etc.
    memset(fd, 0, 8 * sizeof(TacsScalar));
    addResidual(0, time, X, q, dq, ddq, fd);

    // Set the values for the first evaluation
    for (int i = 0; i < 8; i++) {
      q[i] = vars[i] - dh * alpha * perb[i];
      dq[i] = dvars[i] - dh * beta * perb[i];
      ddq[i] = ddvars[i] - dh * gamma * perb[i];
    }

    // Get the residual at vars + alpha*perb, ... etc.
    memset(res, 0, 8 * sizeof(TacsScalar));
    addResidual(0, time, X, q, dq, ddq, res);

    // Form the finite-difference matrix-vector approximation
    for (int i = 0; i < 8; i++) {
      fd[i] = 0.5 * (fd[i] - res[i]) / dh;
    }
#endif  // TACS_USE_COMPLEX

    // Get the Jacobian computed by the element
    memset(mat, 0, 64 * sizeof(TacsScalar));
    addJacobian(0, time, alpha, beta, gamma, X, vars, dvars, ddvars, res, mat);

    // Compute the product: res = J*perb
    // Recall that the Jacobian matrix is stored in row-major order
    memset(res, 0, 8 * sizeof(TacsScalar));
    for (int i = 0; i < 8; i++) {
      for (int j = 0; j < 8; j++) {
        res[i] += mat[8 * i + j] * perb[j];
      }
    }

    // Print out the results to stdout
    char outname[128];
    sprintf(outname, "Jacobian col %d", ii);
    writeErrorComponents(stdout, outname, res, fd, 8);
  }
}

/*
  Retrieve the data associated with the element
*/
void TACSRigidBody::getOutputData(int elemIndex, ElementType etype,
                                  int write_flag, const TacsScalar Xpts[],
                                  const TacsScalar vars[],
                                  const TacsScalar dvars[],
                                  const TacsScalar ddvars[], int ld_data,
                                  TacsScalar *data) {
  // Get the initial vector location
  const TacsScalar *rinit;
  rInit->getVector(&rinit);

  // Write out the displacements at each node
  const TacsScalar *u0 = &vars[0];
  TacsScalar eta = vars[3];
  const TacsScalar *eps = &vars[4];

  // Compute the rotation matrix
  TacsScalar C[9];
  computeRotationMat(eta, eps, C);

  // Get the reference rotation
  const TacsScalar *Cr;
  CRef->getRotation(&Cr);

  // Get the nodal locations for the body
  int npts, nelems;
  const TacsScalar *X;
  viz->getMesh(&npts, &nelems, &X, NULL);

  // Set the locations/variables for all the points from the
  // body-fixed reference frame
  for (int i = 0; i < npts; i++) {
    int index = 0;

    // Compute the initial base-points for each node
    const TacsScalar *x = &X[3 * i];

    if (write_flag & TACS_OUTPUT_NODES) {
      // Write out the nodal locations
      for (int k = 0; k < 3; k++) {
        data[index + k] = TacsRealPart(rinit[k] + x[k]);
      }
      index += 3;
    }
    if (write_flag & TACS_OUTPUT_DISPLACEMENTS) {
      // Compute the new point location
      TacsScalar xr[3], xpt[3];
      mat3x3MultTrans(Cr, x, xr);
      mat3x3MultTrans(C, xr, xpt);

      // xinit = xinit + Cr^{T}*x
      // xfinal = xinit + x + u0 + C^{T}*x
      for (int k = 0; k < 3; k++) {
        data[index + k] = TacsRealPart(u0[k] + xpt[k] - x[k]);
      }
      index += 3;

      // Add the eta variable
      data[index] = TacsRealPart(eta);
      index++;

      // Add the epsilon quaternion components
      for (int k = 0; k < 3; k++) {
        data[index + k] = TacsRealPart(eps[k]);
      }
      index += 3;

      // Add the Lagrange multiplier
      data[index] = TacsRealPart(vars[7]);
    }
    data += ld_data;
  }
}

/*
  Sets the visualization information for the rigid body
*/
void TACSRigidBody::setVisualization(TACSRigidBodyViz *_viz) {
  if (viz) {
    viz->decref();
  }
  viz = _viz;
  viz->incref();
}
