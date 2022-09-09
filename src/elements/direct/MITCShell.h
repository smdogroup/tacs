/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_MITC_SHELL_H
#define TACS_MITC_SHELL_H

/*
  MITC (Mixed interpolation of tensorial components)-based shell
  element.
*/

#include "FElibrary.h"
#include "FSDTStiffness.h"
#include "LargeRotUtils.h"
#include "ShellUtils.h"
#include "TACSElement.h"
#include "TACSShell.h"
#include "TensorToolbox.h"

// Include all the functions from the namespace shellutils
using namespace shellutils;
using namespace largerot;

/*
  A Shell element for general linear and geometrically nonlinear
  analysis.

  This element employs a mixed interpolation of tensorial (strain)
  components (MITC) method to avoid shear locking problems. The flag
  MITCShellType signals whether to use the full nonlinear terms in the
  strain expressions. However, the element unknowns are in-plane
  displacements and small-angle rotations. This limits the element to
  relative small rotations. Stability and post-buckling can be explored
  as long as the rotations remain moderate.

  In the MITC approach, the strain compoennts susceptible to locking
  are interpolated with appropriate low-order polynomials. In plance
  of interpolating the shear and in-plane components directly, the
  tensorial strains are interpolated and transformed back to the local
  coordinate strains at the interpolation points (also called tying
  points). This means that the evaluation of the residuals, stiffness
  matrix and design-type derivatives are more computationally expensive
  than an equivalent displacement-based shell. However, avoiding shear
  locking is a key advantage and these elements should be used whenever
  possible!
*/
template <int order, int tying_order = order>
class MITCShell : public TACSShell {
 public:
  MITCShell(FSDTStiffness *_stiff, ElementBehaviorType type = LINEAR,
            int _componentNum = 0, int _use_lobatto_quadrature = 0);
  ~MITCShell();

  // How many nodes/variables per element
  // ------------------------------------
  int numNodes();
  int numVariables();

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies(double time, TacsScalar *_Te, TacsScalar *_Pe,
                       const TacsScalar Xpts[], const TacsScalar vars[],
                       const TacsScalar dvars[]);

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual(double time, TacsScalar res[], const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[]);

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian(double time, TacsScalar J[], double alpha, double beta,
                   double gamma, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[]);

  // Add the product of the adjoint with the derivative of the design variables
  // --------------------------------------------------------------------------
  void addAdjResProduct(double time, double scale, TacsScalar dvSens[],
                        int dvLen, const TacsScalar psi[],
                        const TacsScalar Xpts[], const TacsScalar vars[],
                        const TacsScalar dvars[], const TacsScalar ddvars[]);

  // Add the product of the adjoint with the derivative of the design variables
  // --------------------------------------------------------------------------
  void addAdjResXptProduct(double time, double scale, TacsScalar XptSens[],
                           const TacsScalar psi[], const TacsScalar Xpts[],
                           const TacsScalar vars[], const TacsScalar dvars[],
                           const TacsScalar ddvars[]);

  // Retrieve a specific time-independent matrix from the element
  // ------------------------------------------------------------
  void getMatType(ElementMatrixType matType, TacsScalar mat[],
                  const TacsScalar Xpts[], const TacsScalar vars[]);

  // Compute the derivative of the inner product w.r.t. design variables
  // -------------------------------------------------------------------
  void addMatDVSensInnerProduct(ElementMatrixType matType, double scale,
                                TacsScalar dvSens[], int dvLen,
                                const TacsScalar psi[], const TacsScalar phi[],
                                const TacsScalar Xpts[],
                                const TacsScalar vars[]);

  // Compute the derivative of the inner product w.r.t. vars[]
  // ---------------------------------------------------------
  void getMatSVSensInnerProduct(ElementMatrixType matType, TacsScalar res[],
                                const TacsScalar psi[], const TacsScalar phi[],
                                const TacsScalar Xpts[],
                                const TacsScalar vars[]);

  // Member functions for evaluating global functions of interest
  // ------------------------------------------------------------
  TACSConstitutive *getConstitutive() { return stiff; }

  // Get the number of Gauss quadrature points
  // -----------------------------------------
  int getNumGaussPts();

  // Get the quadrature points and weights
  // -------------------------------------
  double getGaussWtsPts(const int num, double *pt);

  // Get the shape functions from the element
  // ----------------------------------------
  void getShapeFunctions(const double pt[], double N[]);

  // Return the determinant of the Jacobian at this point
  // ----------------------------------------------------
  TacsScalar getDetJacobian(const double *pt, const TacsScalar Xpts[]);

  // Return the determinant of the Jacobian and its sensitivity at this point
  // ------------------------------------------------------------------------
  TacsScalar getDetJacobianXptSens(TacsScalar *hXptSens, const double *pt,
                                   const TacsScalar Xpts[]);

  // This function returns the strain evaluated at pt
  // ------------------------------------------------
  void getStrain(TacsScalar strain[], const double pt[],
                 const TacsScalar Xpts[], const TacsScalar vars[]);

  // This function returns the sensitivity of the strain w.r.t. Xpts
  // ---------------------------------------------------------------
  void addStrainXptSens(TacsScalar fXptSens[], const double pt[],
                        const TacsScalar scale, const TacsScalar strainSens[],
                        const TacsScalar Xpts[], const TacsScalar vars[]);

  // This function adds the sensitivity of the strain to the state variables
  // -----------------------------------------------------------------------
  void addStrainSVSens(TacsScalar strainSVSens[], const double pt[],
                       const TacsScalar scale, const TacsScalar strainSens[],
                       const TacsScalar Xpts[], const TacsScalar vars[]);

  // Function used for localizing the error to nodes with PU-weights
  // ---------------------------------------------------------------
  void addLocalizedError(double time, TacsScalar err[],
                         const TacsScalar adjoint[], const TacsScalar Xpts[],
                         const TacsScalar vars[]);

  // Functions for post-processing
  // -----------------------------
  void addOutputCount(int *nelems, int *nnodes, int *ncsr);
  void getOutputData(unsigned int out_type, double *data, int ld_data,
                     const TacsScalar Xpts[], const TacsScalar vars[]);
  void getOutputConnectivity(int *con, int node);

  // Set the number of nodes/number of variables per element
  // -------------------------------------------------------
  static const int NUM_NODES = order * order;
  static const int NUM_VARIABLES = 6 * order * order;

 private:
  // Get the partition of unity constraint
  void getPartUnityShapeFunctions(const double pt[], double N[], double Na[],
                                  double Nb[]);

  static const int NUM_G11 = (tying_order - 1) * tying_order;
  static const int NUM_G22 = (tying_order - 1) * tying_order;
  static const int NUM_G12 = (tying_order - 1) * (tying_order - 1);
  static const int NUM_G13 = (tying_order - 1) * tying_order;
  static const int NUM_G23 = (tying_order - 1) * tying_order;

  inline static TacsScalar strain_product(const TacsScalar a[],
                                          const TacsScalar b[]) {
    return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3] +
            a[4] * b[4] + a[5] * b[5] + a[6] * b[6] + a[7] * b[7]);
  }

  // The type of strain expressions to use
  ElementBehaviorType type;
  // The quadrature scheme for integrating the residual/stiffness
  // matrix -- this is Gauss quadrature
  int numGauss;
  const double *gaussWts, *gaussPts;

  // The knot locations
  const double *knots;   // "tying_order" Gauss points
  const double *pknots;  // "tying_order"-1 Gauss points
};

const double MITCShellFirstOrderKnots[2] = {-1.0, 1.0};

template <int order, int tying_order>
MITCShell<order, tying_order>::MITCShell(FSDTStiffness *_stiff,
                                         ElementBehaviorType _type,
                                         int _componentNum,
                                         int use_lobatto_quadrature)
    : TACSShell(_stiff, _componentNum) {
  type = _type;

  if (use_lobatto_quadrature) {
    if (order == 2) {
      numGauss = 3;
      gaussPts = FElibrary::lobattoPts3;
      gaussWts = FElibrary::lobattoWts3;
    } else if (order == 3) {
      numGauss = 4;
      gaussPts = FElibrary::lobattoPts4;
      gaussWts = FElibrary::lobattoWts4;
    } else if (order == 4) {
      numGauss = 5;
      gaussPts = FElibrary::lobattoPts5;
      gaussWts = FElibrary::lobattoWts5;
    } else {
      unsigned int gaussOrder = order;
      numGauss = FElibrary::getGaussPtsWts(gaussOrder, &gaussPts, &gaussWts);
    }
  } else {
    unsigned int gaussOrder = order;
    numGauss = FElibrary::getGaussPtsWts(gaussOrder, &gaussPts, &gaussWts);
  }
  // Get the knot points - the order and order-1-th Gauss points
  if (tying_order == 2) {
    knots = MITCShellFirstOrderKnots;
  } else {
    FElibrary::getGaussPtsWts(tying_order, &knots, NULL);
  }
  FElibrary::getGaussPtsWts(tying_order - 1, &pknots, NULL);
}

template <int order, int tying_order>
MITCShell<order, tying_order>::~MITCShell() {}

template <int order, int tying_order>
int MITCShell<order, tying_order>::numNodes() {
  return NUM_NODES;
}

template <int order, int tying_order>
int MITCShell<order, tying_order>::numVariables() {
  return NUM_VARIABLES;
}

/*
  Compute the kinetic energy and strain energy (potential energy)
  contribution from this element. These are assigned (not added)
  to the output values of Pe and Te.

  input:
  time:    the simulation time
  Xpts:    the nodal locations
  vars:    the element variables
  dvars:   the time derivative of the element variables

  output:
  Te:      the kinetic energy within the body
  Pe:      the potential energy of the body
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::computeEnergies(
    double time, TacsScalar *_Te, TacsScalar *_Pe, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[]) {
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9];

  TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial shear components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the displacement strain at the tying points
  TacsScalar b11[3 * NUM_NODES * NUM_G11], b22[3 * NUM_NODES * NUM_G22];
  TacsScalar b12[3 * NUM_NODES * NUM_G12];
  TacsScalar b13[NUM_VARIABLES * NUM_G13], b23[NUM_VARIABLES * NUM_G23];

  // The stress and strain variables
  TacsScalar strain[NUM_STRESSES], stress[NUM_STRESSES];

  // The derivative of rthe rotation penalty term
  TacsScalar drot[NUM_VARIABLES];

  // The element kinetic and potential energies
  TacsScalar Pe = 0.0, Te = 0.0;

  // Compute the values of the tensorial strain at the tying points
  if (type == LARGE_ROTATION) {
    compute_lr_tying_bmat<order, tying_order>(g11, g22, g12, g23, g13, b11, b22,
                                              b12, b23, b13, knots, pknots,
                                              vars, Xpts);
  } else {
    compute_tying_bmat<order, tying_order>((type == LINEAR), g11, g22, g12, g23,
                                           g13, b11, b22, b12, b23, b13, knots,
                                           pknots, vars, Xpts);
  }

  for (int m = 0; m < numGauss; m++) {
    for (int n = 0; n < numGauss; n++) {
      // Set the quadrature point
      double pt[2];
      pt[0] = gaussPts[n];
      pt[1] = gaussPts[m];

      // Evaluate the stiffness at the parametric point within the
      // element
      TacsScalar At[6], Bt[6], Dt[6], Ats[3];
      TacsScalar kpenalty = stiff->getStiffness(pt, At, Bt, Dt, Ats);

      // Calculate the shape functions and the Jacobian/Hessian of the
      // shell position at the quadrature point
      shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, pt, Xpts);
      compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

      // Compute the transformation from the global coordinates to
      // local shell coordinates
      TacsScalar h = 0.0;
      if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
        h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, Xd,
                              Xdd);
      } else {
        const TacsScalar *axis = stiff->getRefAxis();
        h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi, normal_eta,
                                      axis, Xd, Xdd);
      }
      h = gaussWts[n] * gaussWts[m] * h;

      // Compute the strain and rotation at the qudrature point
      TacsScalar rot = 0.0;
      if (type == LINEAR) {
        linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal, normal_xi,
                           normal_eta);
      } else if (type == NONLINEAR) {
        nonlinear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal,
                              normal_xi, normal_eta);
      } else {
        // Rotation matrix data
        TacsScalar C[9], Ct[27], Ctt[54];

        // Compute the rotation matrices
        TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
        TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
        TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
        compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
        compute_2nd_rate_matrix(Ctt, c1, s1, c2, s2, c3, s3);

        // Evaluate the in-plane rotation term
        rot =
            compute_inplane_penalty(drot, NUM_NODES, Xd, Ud, C, Ct, N, Na, Nb);

        // Calculate the deformation at the current point...
        large_rot_bend_strain(strain, U, Ud, C, Ct, t, tx, ztx, normal,
                              normal_xi, normal_eta);
      }

      // Evaluate the strain interpolation at this point
      tying_interpolation<tying_order>(pt, N11, N22, N12, knots, pknots);
      add_tying_strain<tying_order>(strain, tx, g11, g22, g12, g23, g13, N11,
                                    N22, N12);

      // Compute the stress at the current Gauss point
      stiff->calculateStress(At, Bt, Dt, Ats, strain, stress);

      Pe += 0.5 * h * (strain_product(stress, strain) + kpenalty * rot * rot);
    }
  }

  *_Te = Te;
  *_Pe = Pe;
}

/*
  Compute the residuals of the governing equations of motion.
  Get the element residuals corresponding to the strain energy
  contributions (e.g. not the work terms)

  output:
  res:     the element residual

  input:
  vars:    the element variables
  Xpts:    the element nodal locations in R^{3}
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::addResidual(double time, TacsScalar res[],
                                                const TacsScalar Xpts[],
                                                const TacsScalar vars[],
                                                const TacsScalar dvars[],
                                                const TacsScalar ddvars[]) {
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9];

  TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial shear components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the displacement strain
  TacsScalar b11[3 * NUM_NODES * NUM_G11], b22[3 * NUM_NODES * NUM_G22];
  TacsScalar b12[3 * NUM_NODES * NUM_G12];
  TacsScalar b13[NUM_VARIABLES * NUM_G13], b23[NUM_VARIABLES * NUM_G23];

  // The strain an rotation matrices
  TacsScalar strain[NUM_STRESSES], stress[NUM_STRESSES];
  TacsScalar drot[NUM_VARIABLES];
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

  if (type == LARGE_ROTATION) {
    compute_lr_tying_bmat<order, tying_order>(g11, g22, g12, g23, g13, b11, b22,
                                              b12, b23, b13, knots, pknots,
                                              vars, Xpts);
  } else {
    compute_tying_bmat<order, tying_order>((type == LINEAR), g11, g22, g12, g23,
                                           g13, b11, b22, b12, b23, b13, knots,
                                           pknots, vars, Xpts);
  }

  for (int m = 0; m < numGauss; m++) {
    for (int n = 0; n < numGauss; n++) {
      // Set the quadrature point
      double pt[2];
      pt[0] = gaussPts[n];
      pt[1] = gaussPts[m];

      // Evaluate the stiffness at the parametric point within the
      // element
      TacsScalar At[6], Bt[6], Dt[6], Ats[3];
      TacsScalar kpenalty = stiff->getStiffness(pt, At, Bt, Dt, Ats);

      // Calculate the shape functions and the Jacobian/Hessian of the
      // shell position at the quadrature point
      shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, pt, Xpts);
      compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

      // Compute the transformation from the global coordinates to
      // local shell coordinates
      TacsScalar h = 0.0;
      if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
        h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, Xd,
                              Xdd);
      } else {
        const TacsScalar *axis = stiff->getRefAxis();
        h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi, normal_eta,
                                      axis, Xd, Xdd);
      }
      h = gaussWts[n] * gaussWts[m] * h;

      // Compute the strain and rotation at the qudrature point
      TacsScalar rot = 0.0;
      if (type == LINEAR) {
        linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal, normal_xi,
                           normal_eta);
        linear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, t, tx, ztx, normal,
                         normal_xi, normal_eta);
      } else if (type == NONLINEAR) {
        nonlinear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal,
                              normal_xi, normal_eta);
        nonlinear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, U, Ud, t, tx, ztx,
                            normal, normal_xi, normal_eta);
      } else {
        // Rotation matrix data
        TacsScalar C[9], Ct[27], Ctt[54];

        // Compute the rotation matrices
        TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
        TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
        TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
        compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
        compute_2nd_rate_matrix(Ctt, c1, s1, c2, s2, c3, s3);

        // Evaluate the in-plane rotation term
        rot =
            compute_inplane_penalty(drot, NUM_NODES, Xd, Ud, C, Ct, N, Na, Nb);

        // Calculate the deformation at the current point...
        large_rot_bend_strain(strain, U, Ud, C, Ct, t, tx, ztx, normal,
                              normal_xi, normal_eta);
        large_rot_bend_bmat(B, NUM_NODES, N, Na, Nb, U, Ud, C, Ct, Ctt, t, tx,
                            ztx, normal, normal_xi, normal_eta);
      }

      // Evaluate the strain interpolation at this point
      tying_interpolation<tying_order>(pt, N11, N22, N12, knots, pknots);
      add_tying_strain<tying_order>(strain, tx, g11, g22, g12, g23, g13, N11,
                                    N22, N12);
      add_tying_bmat<tying_order>(B, NUM_NODES, tx, b11, b22, b12, b23, b13,
                                  N11, N22, N12);

      // Compute the stress at the current Gauss point
      stiff->calculateStress(At, Bt, Dt, Ats, strain, stress);

      // Get the pointwise mass at the quadrature point
      TacsScalar mass[2];
      stiff->getPointwiseMass(pt, mass);

      // Compute the accelerations at the current point due to the
      // ddvars input. Store the accelerations in the variable U[].
      TacsScalar Uddot[NUM_DISPS];
      compute_shell_U(NUM_NODES, Uddot, ddvars, N);

      TacsScalar *r = res;
      const TacsScalar *b = B;
      const TacsScalar *br = drot;
      for (int i = 0; i < NUM_NODES; i++) {
        for (int ii = 0; ii < NUM_DISPS; ii++) {
          r[ii] += h * (strain_product(b, stress) + kpenalty * rot * br[0]);
          b += NUM_STRESSES;
          br++;
        }

        // Add the inertial terms from the accelerations
        r[0] += h * N[i] * mass[0] * Uddot[0];
        r[1] += h * N[i] * mass[0] * Uddot[1];
        r[2] += h * N[i] * mass[0] * Uddot[2];

        // Add the inertial terms from the rotations
        TacsScalar d =
            normal[0] * Uddot[3] + normal[1] * Uddot[4] + normal[2] * Uddot[5];
        r[3] += h * N[i] * mass[1] * (Uddot[3] - normal[0] * d);
        r[4] += h * N[i] * mass[1] * (Uddot[4] - normal[1] * d);
        r[5] += h * N[i] * mass[1] * (Uddot[5] - normal[2] * d);

        r += NUM_DISPS;
      }
    }
  }
}

/*
  Add the element tangent stiffness matrix - the exact Jacobian of the
  residual expressions.

  output:
  mat:     the element tangent stiffness matrix
  res:     the element residual

  input:
  vars:    the element variables
  Xpts:    the element nodal locations in R^{3}
  matOr:   the matrix orientation (NORMAL or TRANSPOSE)
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::addJacobian(
    double time, TacsScalar J[], double alpha, double beta, double gamma,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[]) {
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9];

  TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial shear components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the displacement strain
  TacsScalar b11[3 * NUM_NODES * NUM_G11], b22[3 * NUM_NODES * NUM_G22];
  TacsScalar b12[3 * NUM_NODES * NUM_G12];
  TacsScalar b13[NUM_VARIABLES * NUM_G13], b23[NUM_VARIABLES * NUM_G23];

  // The second derivatives of the tying strain
  TacsScalar n13[12 * order * order * (order * order + 1) * NUM_G13];
  TacsScalar n23[12 * order * order * (order * order + 1) * NUM_G23];

  // The stress and strain information
  TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];
  TacsScalar BStress[NUM_STRESSES];
  TacsScalar drot[NUM_VARIABLES];

  // Evaluate the strain and derivative of the strain at the
  // tying points within the element
  if (type == LARGE_ROTATION) {
    compute_lr_tying_nmat<order, tying_order>(g11, g22, g12, g23, g13, b11, b22,
                                              b12, b23, b13, n23, n13, knots,
                                              pknots, vars, Xpts);
  } else {
    compute_tying_bmat<order, tying_order>((type == LINEAR), g11, g22, g12, g23,
                                           g13, b11, b22, b12, b23, b13, knots,
                                           pknots, vars, Xpts);
  }

  for (int m = 0; m < numGauss; m++) {
    for (int n = 0; n < numGauss; n++) {
      // Set the quadrature point
      double pt[2];
      pt[0] = gaussPts[n];
      pt[1] = gaussPts[m];

      // Evaluate the stiffness at the parametric point within the
      // element
      TacsScalar At[6], Bt[6], Dt[6], Ats[3];
      TacsScalar kpenalty = stiff->getStiffness(pt, At, Bt, Dt, Ats);

      // Compute the shape functions and evaluate the surface derivatives
      // at the quadrature point
      shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, pt, Xpts);
      compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

      // Compute the transformation matrix from the global coordinate
      // system to the shell-aligned frame
      TacsScalar h = 0.0;
      if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
        h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, Xd,
                              Xdd);
      } else {
        const TacsScalar *axis = stiff->getRefAxis();
        h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi, normal_eta,
                                      axis, Xd, Xdd);
      }

      // Scale the determinant of the Jacobian transformation by the
      // quadrature weight at this point
      h = gaussWts[n] * gaussWts[m] * h;

      // Store the difference between the rotation variable
      // and the in-plane rotation
      TacsScalar rot = 0.0;

      // Rotation matrix data
      TacsScalar C[9], Ct[27], Ctt[54], Cttt[63];

      // Compute the strain at the current point
      if (type == LINEAR) {
        linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal, normal_xi,
                           normal_eta);
        linear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, t, tx, ztx, normal,
                         normal_xi, normal_eta);
      } else if (type == NONLINEAR) {
        nonlinear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal,
                              normal_xi, normal_eta);
        nonlinear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, U, Ud, t, tx, ztx,
                            normal, normal_xi, normal_eta);
      } else {
        // Compute the rotation matrices
        TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
        TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
        TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
        compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
        compute_2nd_rate_matrix(Ctt, c1, s1, c2, s2, c3, s3);
        compute_3rd_rate_matrix(Cttt, c1, s1, c2, s2, c3, s3);

        // Evaluate the in-plane rotation term
        rot =
            compute_inplane_penalty(drot, NUM_NODES, Xd, Ud, C, Ct, N, Na, Nb);

        // Calculate the strain/bmat at the current point
        large_rot_bend_strain(strain, U, Ud, C, Ct, t, tx, ztx, normal,
                              normal_xi, normal_eta);
        large_rot_bend_bmat(B, NUM_NODES, N, Na, Nb, U, Ud, C, Ct, Ctt, t, tx,
                            ztx, normal, normal_xi, normal_eta);
      }

      // Evaluate the strain interpolation at this point
      tying_interpolation<tying_order>(pt, N11, N22, N12, knots, pknots);

      // Add the interpolated strain and the interpolated b-matrix to the
      // point-wise strain and strain-derivative (B)
      add_tying_strain<tying_order>(strain, tx, g11, g22, g12, g23, g13, N11,
                                    N22, N12);
      add_tying_bmat<tying_order>(B, NUM_NODES, tx, b11, b22, b12, b23, b13,
                                  N11, N22, N12);

      stiff->calculateStress(At, Bt, Dt, Ats, strain, stress);

      // Get the pointwise mass at the quadrature point
      TacsScalar mass[2];
      stiff->getPointwiseMass(pt, mass);

      // Scale the determinant of the Jacobian matrix by the alpha
      // scaling factor
      TacsScalar ha = h * alpha;
      TacsScalar hg = h * gamma;

      for (int i = 0; i < NUM_NODES; i++) {
        for (int ii = 0; ii < NUM_DISPS; ii++) {
          int row = ii + NUM_DISPS * i;

          // Calculate the stress associated with B
          stiff->calculateStress(At, Bt, Dt, Ats, &B[row * NUM_STRESSES],
                                 BStress);

          for (int j = 0; j <= i; j++) {
            int end = (j == i ? ii : NUM_DISPS - 1);
            for (int jj = 0; jj <= end; jj++) {
              int col = jj + NUM_DISPS * j;

              // The regular element matrix
              J[col + row * NUM_VARIABLES] +=
                  ha * (strain_product(BStress, &B[col * NUM_STRESSES]) +
                        kpenalty * drot[row] * drot[col]);
            }
          }
        }
      }

      if (type == NONLINEAR) {
        nonlinear_bend_stress_bmat(J, NUM_NODES, ha, stress, N, Na, Nb, t, tx,
                                   ztx, normal, normal_xi, normal_eta);
        add_nonlinear_tying_stress_nmat<order, tying_order>(
            J, ha, stress, tx, N11, N22, N12, knots, pknots, Xpts);
      } else if (type == LARGE_ROTATION) {
        // Add the second-derivative contributions from the bending strain
        add_large_rot_bend_stress_bmat(J, NUM_NODES, ha, stress, N, Na, Nb, U,
                                       Ud, C, Ct, Ctt, Cttt, t, tx, ztx, normal,
                                       normal_xi, normal_eta);

        // Add the contributions to the second derivative of the tying strain
        add_lr_tying_stress_nmat<order, tying_order>(
            J, ha, stress, n13, n23, tx, N11, N22, N12, knots, pknots);

        // Add the second derivative of the in-plane penalty
        add_inplane_penalty(J, NUM_NODES, ha * kpenalty * rot, Xd, Ud, Ct, Ctt,
                            N, Na, Nb);
      }

      // Add the kinetic energy terms from the displacement
      TacsScalar Am = hg * mass[0];
      for (int i = 0; i < NUM_NODES; i++) {
        for (int j = 0; j < NUM_NODES; j++) {
          for (int ii = 0; ii < 3; ii++) {
            int row = ii + NUM_DISPS * i;
            int col = ii + NUM_DISPS * j;
            J[col + row * NUM_VARIABLES] += Am * N[i] * N[j];
          }
        }
      }

      // Add the contribution due to the rotational terms
      TacsScalar Dm = hg * mass[1];
      TacsScalar D[9];
      D[0] = Dm * (1.0 - normal[0] * normal[0]);
      D[1] = -Dm * normal[0] * normal[1];
      D[2] = -Dm * normal[0] * normal[2];

      D[3] = -Dm * normal[1] * normal[0];
      D[4] = Dm * (1.0 - normal[1] * normal[1]);
      D[5] = -Dm * normal[1] * normal[2];

      D[6] = -Dm * normal[2] * normal[0];
      D[7] = -Dm * normal[2] * normal[1];
      D[8] = Dm * (1.0 - normal[2] * normal[2]);

      // Add the values to the matrix
      for (int i = 0; i < NUM_NODES; i++) {
        for (int j = 0; j < NUM_NODES; j++) {
          for (int ii = 0; ii < 3; ii++) {
            int row = 3 + ii + NUM_DISPS * i;
            for (int jj = 0; jj < 3; jj++) {
              int col = 3 + jj + NUM_DISPS * j;
              J[col + row * NUM_VARIABLES] += D[ii + 3 * jj] * N[i] * N[j];
            }
          }
        }
      }
    }
  }

  // Copy over the matrix
  // Take the lower triangle and copy to the upper triangle
  for (int row = 0; row < NUM_VARIABLES; row++) {
    for (int col = row + 1; col < NUM_VARIABLES; col++) {
      J[col + row * NUM_VARIABLES] = J[row + col * NUM_VARIABLES];
    }
  }
}

/*
  Evaluate the element matrix of a specified type

  output:
  mat:         the element matrix of the specified type

  input:
  matType:     the matrix type (e.g. MASS_MATRIX)
  vars:        the element variables
  Xpts:        the nodal coordinates in R^{3}
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::getMatType(ElementMatrixType matType,
                                               TacsScalar *mat,
                                               const TacsScalar Xpts[],
                                               const TacsScalar vars[]) {
  memset(mat, 0, NUM_VARIABLES * NUM_VARIABLES * sizeof(TacsScalar));

  if (matType == STIFFNESS_MATRIX) {
    // Geometric data
    TacsScalar X[3], Xd[9], Xdd[9];
    TacsScalar normal[3], normal_xi[3], normal_eta[3];

    // Transformation and the transformation derivative w.r.t. zeta
    TacsScalar t[9], tx[9], ztx[9];

    TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

    // Interpolations for the shear components
    double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

    // The interpolated tensorial shear components
    TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
    TacsScalar g13[NUM_G13], g23[NUM_G23];

    // The derivatives of the displacement strain
    TacsScalar b11[3 * NUM_NODES * NUM_G11], b22[3 * NUM_NODES * NUM_G22];
    TacsScalar b12[3 * NUM_NODES * NUM_G12];
    TacsScalar b13[NUM_VARIABLES * NUM_G13], b23[NUM_VARIABLES * NUM_G23];

    // The second derivatives of the tying strain
    TacsScalar n13[12 * order * order * (order * order + 1) * NUM_G13];
    TacsScalar n23[12 * order * order * (order * order + 1) * NUM_G23];

    // The stress and strain information
    TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];
    TacsScalar B[NUM_STRESSES * NUM_VARIABLES];
    TacsScalar BStress[NUM_STRESSES];
    TacsScalar drot[NUM_VARIABLES];

    // Evaluate the strain and derivative of the strain at the
    // tying points within the element
    if (type == LARGE_ROTATION) {
      compute_lr_tying_bmat<order, tying_order>(g11, g22, g12, g23, g13, b11,
                                                b22, b12, b23, b13, knots,
                                                pknots, vars, Xpts);
    } else {
      compute_tying_bmat<order, tying_order>((type == LINEAR), g11, g22, g12,
                                             g23, g13, b11, b22, b12, b23, b13,
                                             knots, pknots, vars, Xpts);
    }

    for (int m = 0; m < numGauss; m++) {
      for (int n = 0; n < numGauss; n++) {
        // Set the quadrature point
        double pt[2];
        pt[0] = gaussPts[n];
        pt[1] = gaussPts[m];

        // Evaluate the stiffness at the parametric point within the
        // element
        TacsScalar At[6], Bt[6], Dt[6], Ats[3];
        TacsScalar kpenalty = stiff->getStiffness(pt, At, Bt, Dt, Ats);

        // Compute the shape functions and evaluate the surface derivatives
        // at the quadrature point
        shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, pt, Xpts);
        compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

        // Compute the transformation matrix from the global coordinate
        // system to the shell-aligned frame
        TacsScalar h = 0.0;
        if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
          h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, Xd,
                                Xdd);
        } else {
          const TacsScalar *axis = stiff->getRefAxis();
          h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi,
                                        normal_eta, axis, Xd, Xdd);
        }

        // Scale the determinant of the Jacobian transformation by the
        // quadrature weight at this point
        h = gaussWts[n] * gaussWts[m] * h;

        // Store the difference between the rotation variable
        // and the in-plane rotation
        TacsScalar rot = 0.0;

        // Rotation matrix data
        TacsScalar C[9], Ct[27], Ctt[54], Cttt[63];

        // Compute the strain at the current point
        if (type == LINEAR) {
          linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal, normal_xi,
                             normal_eta);
          linear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, t, tx, ztx, normal,
                           normal_xi, normal_eta);
        } else if (type == NONLINEAR) {
          nonlinear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal,
                                normal_xi, normal_eta);
          nonlinear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, U, Ud, t, tx, ztx,
                              normal, normal_xi, normal_eta);
        } else {
          // Compute the rotation matrices
          TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
          TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
          TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
          compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
          compute_2nd_rate_matrix(Ctt, c1, s1, c2, s2, c3, s3);
          compute_3rd_rate_matrix(Cttt, c1, s1, c2, s2, c3, s3);

          // Evaluate the in-plane rotation term
          rot = compute_inplane_penalty(drot, NUM_NODES, Xd, Ud, C, Ct, N, Na,
                                        Nb);

          // Calculate the strain/bmat at the current point
          large_rot_bend_strain(strain, U, Ud, C, Ct, t, tx, ztx, normal,
                                normal_xi, normal_eta);
          large_rot_bend_bmat(B, NUM_NODES, N, Na, Nb, U, Ud, C, Ct, Ctt, t, tx,
                              ztx, normal, normal_xi, normal_eta);
        }

        // Evaluate the strain interpolation at this point
        tying_interpolation<tying_order>(pt, N11, N22, N12, knots, pknots);

        // Add the interpolated strain and the interpolated b-matrix to the
        // point-wise strain and strain-derivative (B)
        add_tying_strain<tying_order>(strain, tx, g11, g22, g12, g23, g13, N11,
                                      N22, N12);
        add_tying_bmat<tying_order>(B, NUM_NODES, tx, b11, b22, b12, b23, b13,
                                    N11, N22, N12);

        stiff->calculateStress(At, Bt, Dt, Ats, strain, stress);

        // Get the pointwise mass at the quadrature point
        TacsScalar mass[2];
        stiff->getPointwiseMass(pt, mass);

        for (int i = 0; i < NUM_NODES; i++) {
          for (int ii = 0; ii < NUM_DISPS; ii++) {
            int row = ii + NUM_DISPS * i;

            // Calculate the stress associated with B
            stiff->calculateStress(At, Bt, Dt, Ats, &B[row * NUM_STRESSES],
                                   BStress);

            for (int j = 0; j <= i; j++) {
              int end = (j == i ? ii : NUM_DISPS - 1);
              for (int jj = 0; jj <= end; jj++) {
                int col = jj + NUM_DISPS * j;

                // The regular element matrix
                mat[col + row * NUM_VARIABLES] +=
                    h * (strain_product(BStress, &B[col * NUM_STRESSES]) +
                         kpenalty * drot[row] * drot[col]);
              }
            }
          }
        }

        if (type == NONLINEAR) {
          nonlinear_bend_stress_bmat(mat, NUM_NODES, h, stress, N, Na, Nb, t,
                                     tx, ztx, normal, normal_xi, normal_eta);
          add_nonlinear_tying_stress_nmat<order, tying_order>(
              mat, h, stress, tx, N11, N22, N12, knots, pknots, Xpts);
        } else if (type == LARGE_ROTATION) {
          // Add the second-derivative contributions from the bending strain
          add_large_rot_bend_stress_bmat(mat, NUM_NODES, h, stress, N, Na, Nb,
                                         U, Ud, C, Ct, Ctt, Cttt, t, tx, ztx,
                                         normal, normal_xi, normal_eta);

          // Add the contributions to the second derivative of the tying strain
          add_lr_tying_stress_nmat<order, tying_order>(
              mat, h, stress, n13, n23, tx, N11, N22, N12, knots, pknots);

          // Add the second derivative of the in-plane penalty
          add_inplane_penalty(mat, NUM_NODES, h * kpenalty * rot, Xd, Ud, Ct,
                              Ctt, N, Na, Nb);
        }

        // Copy over the matrix
        // Take the lower triangle and copy to the upper triangle
        for (int row = 0; row < NUM_VARIABLES; row++) {
          for (int col = row + 1; col < NUM_VARIABLES; col++) {
            mat[col + row * NUM_VARIABLES] = mat[row + col * NUM_VARIABLES];
          }
        }
      }
    }
  } else if (matType == MASS_MATRIX) {
    for (int m = 0; m < numGauss; m++) {
      for (int n = 0; n < numGauss; n++) {
        // Set the quadrature point within the element
        double pt[2];
        pt[0] = gaussPts[n];
        pt[1] = gaussPts[m];

        // Evaluate the shape functions and the derivative of the
        // shell surface location along the u/v directions
        TacsScalar X[3], Xd[9], normal[3];
        double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
        shell_jacobian(order, X, Xd, N, Na, Nb, pt, Xpts);

        // Compute the normal vector and set it in the final row
        // of the Xd transformation matrix
        Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
        Tensor::normalize3D(normal);
        for (int i = 0; i < 3; i++) {
          Xd[6 + i] = normal[i];
        }

        // Compute the determinant of the jacobian transformation and
        // scale it by the qudrature weight
        TacsScalar h = FElibrary::jacobian3d(Xd);
        h = gaussWts[n] * gaussWts[m] * h;

        // Get the pointwise mass at the quadrature point
        TacsScalar mass[2];
        stiff->getPointwiseMass(pt, mass);

        // Add the kinetic energy terms from the displacement
        TacsScalar Am = mass[0] * h;
        for (int i = 0; i < NUM_NODES; i++) {
          for (int j = 0; j < NUM_NODES; j++) {
            for (int ii = 0; ii < 3; ii++) {
              int row = ii + NUM_DISPS * i;
              int col = ii + NUM_DISPS * j;
              mat[col + row * NUM_VARIABLES] += Am * N[i] * N[j];
            }
          }
        }

        // Add the contribution due to the rotational terms
        TacsScalar Dm = mass[1] * h;
        TacsScalar D[9];
        D[0] = Dm * (1.0 - normal[0] * normal[0]);
        D[1] = -Dm * normal[0] * normal[1];
        D[2] = -Dm * normal[0] * normal[2];

        D[3] = -Dm * normal[1] * normal[0];
        D[4] = Dm * (1.0 - normal[1] * normal[1]);
        D[5] = -Dm * normal[1] * normal[2];

        D[6] = -Dm * normal[2] * normal[0];
        D[7] = -Dm * normal[2] * normal[1];
        D[8] = Dm * (1.0 - normal[2] * normal[2]);

        // Add the values to the matrix
        for (int i = 0; i < NUM_NODES; i++) {
          for (int j = 0; j < NUM_NODES; j++) {
            for (int ii = 0; ii < 3; ii++) {
              int row = 3 + ii + NUM_DISPS * i;
              for (int jj = 0; jj < 3; jj++) {
                int col = 3 + jj + NUM_DISPS * j;
                mat[col + row * NUM_VARIABLES] += D[ii + 3 * jj] * N[i] * N[j];
              }
            }
          }
        }
      }
    }
  } else if (matType == GEOMETRIC_STIFFNESS_MATRIX && type == LINEAR) {
    // Geometric data
    TacsScalar X[3], Xd[9], Xdd[9];
    TacsScalar normal[3], normal_xi[3], normal_eta[3];

    // Transformation and the transformation derivative w.r.t. zeta
    TacsScalar t[9], tx[9], ztx[9];

    // Displacements/rotations and shape functions
    TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

    // Interpolations for the shear components
    double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

    // The interpolated tensorial shear components
    TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
    TacsScalar g13[NUM_G13], g23[NUM_G23];

    // The stress/strain values
    TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];

    // Compute the tying strain
    const int is_linear = 1;
    compute_tying_strain<order, tying_order>(is_linear, g11, g22, g12, g23, g13,
                                             knots, pknots, vars, Xpts);

    for (int m = 0; m < numGauss; m++) {
      for (int n = 0; n < numGauss; n++) {
        // Set the quadrature point
        double pt[2];
        pt[0] = gaussPts[n];
        pt[1] = gaussPts[m];

        // Evaluate the stiffness at the parametric point within the
        // element
        TacsScalar At[6], Bt[6], Dt[6], Ats[3];
        stiff->getStiffness(pt, At, Bt, Dt, Ats);

        // Calculate the shape functions
        shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, pt, Xpts);
        compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

        // Compute the transformation
        TacsScalar h = 0.0;
        if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
          h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, Xd,
                                Xdd);
        } else {
          const TacsScalar *axis = stiff->getRefAxis();
          h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi,
                                        normal_eta, axis, Xd, Xdd);
        }
        h = gaussWts[n] * gaussWts[m] * h;

        // Calculate the strain/bmat at the current point...
        TacsScalar rot;
        linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal, normal_xi,
                           normal_eta);

        // Evaluate the strain interpolation at this point
        tying_interpolation<tying_order>(pt, N11, N22, N12, knots, pknots);
        add_tying_strain<tying_order>(strain, tx, g11, g22, g12, g23, g13, N11,
                                      N22, N12);

        // Compute the stress at this Gauss point
        stiff->calculateStress(At, Bt, Dt, Ats, strain, stress);

        // Add the geometric stiffness terms
        nonlinear_bend_stress_bmat(mat, NUM_NODES, h, stress, N, Na, Nb, t, tx,
                                   ztx, normal, normal_xi, normal_eta);
        add_nonlinear_tying_stress_nmat<order, tying_order>(
            mat, h, stress, tx, N11, N22, N12, knots, pknots, Xpts);
      }
    }

    // Copy over the matrix
    // Take the lower triangle and copy to the upper triangle
    for (int row = 0; row < NUM_VARIABLES; row++) {
      for (int col = row + 1; col < NUM_VARIABLES; col++) {
        mat[col + row * NUM_VARIABLES] = mat[row + col * NUM_VARIABLES];
      }
    }
  }
}

/*
  Compute the derivative of the product of the residual with the
  adjoint vector psi with respect to the design variables.  The result
  is added to a design-variable array.

  input:
  time:     the simulation time
  scale:    the scaling factor applied to the derivative
  dvLen:    the length of the design variable vector
  psi:      the left-multiplying vector
  phi:      the right-multiplying vector
  vars:     the element state variables
  Xpts:     the nodal locations for this element

  output:
  dvSens:   the result is added to this vector
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::addAdjResProduct(
    double time, double scale, TacsScalar fdvSens[], int dvLen,
    const TacsScalar psi[], const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[]) {
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9];

  TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial shear components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the displacement strain
  TacsScalar b11[3 * NUM_NODES * NUM_G11], b22[3 * NUM_NODES * NUM_G22];
  TacsScalar b12[3 * NUM_NODES * NUM_G12];
  TacsScalar b13[NUM_VARIABLES * NUM_G13], b23[NUM_VARIABLES * NUM_G23];

  TacsScalar strain[NUM_STRESSES];
  TacsScalar drot[NUM_VARIABLES];
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

  // Compute the tying terms in the matrix
  if (type == LARGE_ROTATION) {
    compute_lr_tying_bmat<order, tying_order>(g11, g22, g12, g23, g13, b11, b22,
                                              b12, b23, b13, knots, pknots,
                                              vars, Xpts);
  } else {
    compute_tying_bmat<order, tying_order>((type == LINEAR), g11, g22, g12, g23,
                                           g13, b11, b22, b12, b23, b13, knots,
                                           pknots, vars, Xpts);
  }

  for (int m = 0; m < numGauss; m++) {
    for (int n = 0; n < numGauss; n++) {
      // Set the quadrature point within the element
      double pt[2];
      pt[0] = gaussPts[n];
      pt[1] = gaussPts[m];

      // Compute the shape functions and the derivatives of the
      // position of the shell with respect to the coordinates
      shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, pt, Xpts);
      compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

      // Compute the coordinate transformation
      TacsScalar h = 0.0;
      if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
        h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, Xd,
                              Xdd);
      } else {
        const TacsScalar *axis = stiff->getRefAxis();
        h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi, normal_eta,
                                      axis, Xd, Xdd);
      }

      // Scale the determinant of the Jacobian transformation by the
      // quadrature weight
      h = scale * gaussWts[n] * gaussWts[m] * h;

      // Calculate the deformation at the current point...
      TacsScalar rot = 0.0;
      if (type == LINEAR) {
        linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal, normal_xi,
                           normal_eta);
        linear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, t, tx, ztx, normal,
                         normal_xi, normal_eta);
      } else if (type == NONLINEAR) {
        nonlinear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal,
                              normal_xi, normal_eta);
        nonlinear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, U, Ud, t, tx, ztx,
                            normal, normal_xi, normal_eta);
      } else {
        // Rotation matrix data
        TacsScalar C[9], Ct[27], Ctt[54];

        // Compute the rotation matrices
        TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
        TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
        TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
        compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
        compute_2nd_rate_matrix(Ctt, c1, s1, c2, s2, c3, s3);

        // Evaluate the in-plane rotation term
        rot =
            compute_inplane_penalty(drot, NUM_NODES, Xd, Ud, C, Ct, N, Na, Nb);

        // Calculate the deformation at the current point...
        large_rot_bend_strain(strain, U, Ud, C, Ct, t, tx, ztx, normal,
                              normal_xi, normal_eta);
        large_rot_bend_bmat(B, NUM_NODES, N, Na, Nb, U, Ud, C, Ct, Ctt, t, tx,
                            ztx, normal, normal_xi, normal_eta);
      }

      // Evaluate the strain interpolation at this point
      tying_interpolation<tying_order>(pt, N11, N22, N12, knots, pknots);
      add_tying_strain<tying_order>(strain, tx, g11, g22, g12, g23, g13, N11,
                                    N22, N12);
      add_tying_bmat<tying_order>(B, NUM_NODES, tx, b11, b22, b12, b23, b13,
                                  N11, N22, N12);

      // Compute the product of psi^{T}*B^{T}
      TacsScalar bpsi[NUM_STRESSES], brot = 0.0;
      memset(bpsi, 0, NUM_STRESSES * sizeof(TacsScalar));

      TacsScalar *b = B, *dr = drot;
      const TacsScalar *ps = psi;
      for (int i = 0; i < NUM_VARIABLES; i++) {
        bpsi[0] += ps[0] * b[0];
        bpsi[1] += ps[0] * b[1];
        bpsi[2] += ps[0] * b[2];
        bpsi[3] += ps[0] * b[3];
        bpsi[4] += ps[0] * b[4];
        bpsi[5] += ps[0] * b[5];
        bpsi[6] += ps[0] * b[6];
        bpsi[7] += ps[0] * b[7];

        brot += ps[0] * dr[0];
        b += NUM_STRESSES;
        dr++;
        ps++;
      }

      // Scale the strain and rotational terms by scale*h
      for (int i = 0; i < NUM_STRESSES; i++) {
        bpsi[i] *= h;
      }
      brot *= h;

      // Add the term: scale*psi^{T}*B^{T}*dC/dx*strain to the vector
      // dvSens - Note that this is much more efficient than computing
      // the terms component by component
      stiff->addStiffnessDVSens(pt, strain, bpsi, brot * rot, fdvSens, dvLen);

      // Compute the accelerations at the current point due to the
      // ddvars input. Store the accelerations in the variable U[].
      TacsScalar Uddot[NUM_DISPS];
      compute_shell_U(NUM_NODES, Uddot, ddvars, N);

      // Set the scalar for the mass matrix
      TacsScalar mscale[2] = {0.0, 0.0};
      const TacsScalar *p = psi;
      for (int i = 0; i < NUM_NODES; i++) {
        mscale[0] +=
            h * N[i] * (Uddot[0] * p[0] + Uddot[1] * p[1] + Uddot[2] * p[2]);

        TacsScalar d =
            normal[0] * Uddot[3] + normal[1] * Uddot[4] + normal[2] * Uddot[5];
        mscale[1] += h * N[i] *
                     ((Uddot[3] - normal[0] * d) * p[3] +
                      (Uddot[4] - normal[1] * d) * p[4] +
                      (Uddot[5] - normal[2] * d) * p[5]);
        p += NUM_DISPS;
      }

      // Add the sensitivity term from the mass matrix contribution
      stiff->addPointwiseMassDVSens(pt, mscale, fdvSens, dvLen);
    }
  }
}

/*
  Evaluate the derivative of the element residuals with respect
  to the nodal coordinates e.g res = dR/dXpts

  output:
  res:  the derivative of the residuals w.r.t. the element nodes

  input:
  vars:    the element variables
  Xpts:    the element nodal locations
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::addAdjResXptProduct(
    double time, double scale, TacsScalar fXptSens[], const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[]) {
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];
  TacsScalar dnormal[9 * NUM_NODES], dnormal_xi[9 * NUM_NODES],
      dnormal_eta[9 * NUM_NODES];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9];
  TacsScalar dt[27 * NUM_NODES], dtx[27 * NUM_NODES], dztx[27 * NUM_NODES];

  TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial shear components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the tensorial tying strain components
  TacsScalar dg11[3 * NUM_G11 * NUM_NODES], dg22[3 * NUM_G22 * NUM_NODES],
      dg12[3 * NUM_G12 * NUM_NODES];
  TacsScalar dg13[3 * NUM_G13 * NUM_NODES], dg23[3 * NUM_G23 * NUM_NODES];

  // The derivatives of the displacement strain
  TacsScalar b11[3 * NUM_NODES * NUM_G11], b22[3 * NUM_NODES * NUM_G22];
  TacsScalar b12[3 * NUM_NODES * NUM_G12];
  TacsScalar b13[NUM_VARIABLES * NUM_G13], b23[NUM_VARIABLES * NUM_G23];

  // The stress and strain and their sensitivities
  TacsScalar strain[NUM_STRESSES], stress[NUM_STRESSES];
  TacsScalar dstrain[3 * NUM_STRESSES * NUM_NODES], dstress[NUM_STRESSES];

  // The rotational contributions
  TacsScalar drot[NUM_VARIABLES], srot[NUM_VARIABLES];

  // The derivative of the strain w.r.t. the nodal displacements
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

  // The sensitivity of the determinant of the Jacobian
  TacsScalar dh[3 * NUM_NODES];

  compute_tying_strain_sens<order, tying_order>(
      (type == LINEAR), g11, g22, g12, g23, g13, dg11, dg22, dg12, dg23, dg13,
      knots, pknots, vars, Xpts);
  compute_tying_bmat<order, tying_order>((type == LINEAR), g11, g22, g12, g23,
                                         g13, b11, b22, b12, b23, b13, knots,
                                         pknots, vars, Xpts);

  for (int m = 0; m < numGauss; m++) {
    for (int n = 0; n < numGauss; n++) {
      // Set the quadrature point
      double pt[2];
      pt[0] = gaussPts[n];
      pt[1] = gaussPts[m];

      // Evaluate the stiffness at the parametric point within the
      // element
      TacsScalar At[6], Bt[6], Dt[6], Ats[3];
      TacsScalar k_penalty = stiff->getStiffness(pt, At, Bt, Dt, Ats);

      // Calculate the shape functions and the surface derivatives
      shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, pt, Xpts);

      // Compute the values of U and Ud - the variables and their
      // derivatives w.r.t. the local coordinates
      compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

      // Evaluate the strain interpolation at this point
      tying_interpolation<tying_order>(pt, N11, N22, N12, knots, pknots);

      TacsScalar h = 0.0;
      if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
        h = compute_transform_sens(dh, t, dt, tx, dtx, ztx, dztx, normal,
                                   dnormal, normal_xi, dnormal_xi, normal_eta,
                                   dnormal_eta, Xd, Xdd, Na, Nb, Naa, Nab, Nbb,
                                   NUM_NODES);
      } else {
        const TacsScalar *axis = stiff->getRefAxis();
        h = compute_transform_refaxis_sens(
            dh, t, dt, tx, dtx, ztx, dztx, normal, dnormal, normal_xi,
            dnormal_xi, normal_eta, dnormal_eta, axis, Xd, Xdd, Na, Nb, Naa,
            Nab, Nbb, NUM_NODES);
      }
      h = scale * gaussWts[n] * gaussWts[m] * h;

      // Evaluate the strain and the derivative of the strain w.r.t.
      // the displacements and the nodal coordinates
      TacsScalar rot;
      if (type == LINEAR) {
        linear_bend_strain_sens(strain, dstrain, &rot, srot, U, Ud, t, dt, tx,
                                dtx, ztx, dztx, normal, dnormal, normal_xi,
                                dnormal_xi, normal_eta, dnormal_eta,
                                3 * NUM_NODES);
        linear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, t, tx, ztx, normal,
                         normal_xi, normal_eta);
      } else if (type == NONLINEAR) {
        nonlinear_bend_strain_sens(strain, dstrain, &rot, srot, U, Ud, t, dt,
                                   tx, dtx, ztx, dztx, normal, dnormal,
                                   normal_xi, dnormal_xi, normal_eta,
                                   dnormal_eta, 3 * NUM_NODES);
        nonlinear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, U, Ud, t, tx, ztx,
                            normal, normal_xi, normal_eta);
      }

      // Add the sensitivities of the tying strain points into the
      // strain senstivities
      add_tying_strain_sens<tying_order>(strain, dstrain, tx, dtx, g11, g22,
                                         g12, g23, g13, dg11, dg22, dg12, dg23,
                                         dg13, N11, N22, N12);

      // Add the tying strain to the B-matrix
      add_tying_bmat<tying_order>(B, NUM_NODES, tx, b11, b22, b12, b23, b13,
                                  N11, N22, N12);

      // Calculate the stress at the current point
      stiff->calculateStress(At, Bt, Dt, Ats, strain, stress);

      for (int k = 0; k < 3 * NUM_NODES; k++) {
        dh[k] = scale * gaussWts[n] * gaussWts[m] * dh[k];
        stiff->calculateStress(At, Bt, Dt, Ats, &dstrain[k * NUM_STRESSES],
                               dstress);

        for (int row = 0; row < NUM_VARIABLES; row++) {
          fXptSens[k] +=
              psi[row] *
              (h * (strain_product(&B[NUM_STRESSES * row], dstress) +
                    k_penalty * (srot[k] * drot[row])) +
               dh[k] * (strain_product(&B[NUM_STRESSES * row], stress) +
                        k_penalty * rot * drot[row]));
        }
      }

      if (type == LINEAR) {
        add_linear_bend_bmat_sens(
            fXptSens, psi, NUM_NODES, h, h * k_penalty * rot, stress, N, Na, Nb,
            t, dt, tx, dtx, ztx, dztx, normal, dnormal, normal_xi, dnormal_xi,
            normal_eta, dnormal_eta, 3 * NUM_NODES);
      } else if (type == NONLINEAR) {
        add_nonlinear_bend_bmat_sens(
            fXptSens, psi, NUM_NODES, h, h * k_penalty * rot, stress, N, Na, Nb,
            U, Ud, t, dt, tx, dtx, ztx, dztx, normal, dnormal, normal_xi,
            dnormal_xi, normal_eta, dnormal_eta, 3 * NUM_NODES);
      }

      add_tying_bmat_sens<order, tying_order>((type == LINEAR), fXptSens, psi,
                                              h, stress, tx, dtx, knots, pknots,
                                              vars, Xpts, N11, N22, N12);

      // Get the pointwise mass at the quadrature point
      TacsScalar mass[2];
      stiff->getPointwiseMass(pt, mass);

      // Compute the accelerations at the current point due to the
      // ddvars input. Store the accelerations in the variable U[].
      TacsScalar Uddot[NUM_DISPS];
      compute_shell_U(NUM_NODES, Uddot, ddvars, N);

      // Find the derivative of the inner product w.r.t. the
      // determinant of the Jacobian transformation and the derivative
      // w.r.t. the normal direction
      double alpha = scale * gaussWts[n] * gaussWts[m];
      TacsScalar hd = 0.0;
      TacsScalar nd[3] = {0.0, 0.0, 0.0};
      const TacsScalar *p = psi;
      for (int i = 0; i < NUM_NODES; i++) {
        hd += alpha * mass[0] * N[i] *
              (Uddot[0] * p[0] + Uddot[1] * p[1] + Uddot[2] * p[2]);

        TacsScalar d =
            normal[0] * Uddot[3] + normal[1] * Uddot[4] + normal[2] * Uddot[5];
        hd += alpha * mass[1] * N[i] *
              ((Uddot[3] - normal[0] * d) * p[3] +
               (Uddot[4] - normal[1] * d) * p[4] +
               (Uddot[5] - normal[2] * d) * p[5]);

        // Find the derivative of Psi^{T}*(I - n*n^{T})*U w.r.t n:
        TacsScalar psid =
            (p[3] * normal[0] + p[4] * normal[1] + p[5] * normal[2]);

        nd[0] -= h * mass[1] * N[i] * (p[3] * d + Uddot[3] * psid);
        nd[1] -= h * mass[1] * N[i] * (p[4] * d + Uddot[4] * psid);
        nd[2] -= h * mass[1] * N[i] * (p[5] * d + Uddot[5] * psid);

        p += NUM_DISPS;
      }

      // Add the derivative to the fXptSens array
      for (int i = 0; i < NUM_NODES; i++) {
        for (int k = 0; k < 3; k++) {
          TacsScalar XdSens[9];
          XdSens[0] = XdSens[1] = XdSens[2] = 0.0;
          XdSens[3] = XdSens[4] = XdSens[5] = 0.0;
          XdSens[k] = Na[i];
          XdSens[3 + k] = Nb[i];

          // Compute the derivative of the cross product
          Tensor::crossProduct3DSens(&Xd[6], &XdSens[6], &Xd[0], &Xd[3],
                                     &XdSens[0], &XdSens[3]);

          // Compute the derivative of the normal vector
          TacsScalar snrm;
          Tensor::normalize3DSens(&snrm, &Xd[6], &XdSens[6]);

          TacsScalar hXptSens;
          FElibrary::jacobian3dSens(Xd, XdSens, &hXptSens);

          fXptSens[3 * i + k] +=
              hXptSens * hd +
              (nd[0] * XdSens[6] + nd[1] * XdSens[7] + nd[2] * XdSens[8]);
        }
      }
    }
  }
}

/*
  Compute the derivative of the inner product of the stiffness or mass
  matrices with the given vectors psi and phi with respect to the
  design variables. The result is a vector which is the length of the
  number of design variables.

  input:
  matType:  the type of matrix
  scale:    the scaling factor applied to the derivative
  dvLen:    the length of the design variable vector
  psi:      the left-multiplying vector
  phi:      the right-multiplying vector
  vars:     the element state variables
  Xpts:     the nodal locations for this element

  output:
  dvSens:   the result is added to this vector
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::addMatDVSensInnerProduct(
    ElementMatrixType matType, double scale, TacsScalar fdvSens[], int dvLen,
    const TacsScalar psi[], const TacsScalar phi[], const TacsScalar Xpts[],
    const TacsScalar vars[]) {
  if (matType == STIFFNESS_MATRIX) {
    // Geometric data
    TacsScalar X[3], Xd[9], Xdd[9];
    TacsScalar normal[3], normal_xi[3], normal_eta[3];

    // Transformation and the transformation derivative w.r.t. zeta
    TacsScalar t[9], tx[9], ztx[9];

    TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

    // Interpolations for the shear components
    double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

    // The interpolated tensorial shear components
    TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
    TacsScalar g13[NUM_G13], g23[NUM_G23];

    // The derivatives of the displacement strain
    TacsScalar b11[3 * NUM_NODES * NUM_G11], b22[3 * NUM_NODES * NUM_G22];
    TacsScalar b12[3 * NUM_NODES * NUM_G12];
    TacsScalar b13[NUM_VARIABLES * NUM_G13], b23[NUM_VARIABLES * NUM_G23];

    TacsScalar strain[NUM_STRESSES];
    TacsScalar drot[NUM_VARIABLES];
    TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

    // Compute the tying terms in the matrix
    if (type == LARGE_ROTATION) {
      compute_lr_tying_bmat<order, tying_order>(g11, g22, g12, g23, g13, b11,
                                                b22, b12, b23, b13, knots,
                                                pknots, vars, Xpts);
    } else {
      compute_tying_bmat<order, tying_order>((type == LINEAR), g11, g22, g12,
                                             g23, g13, b11, b22, b12, b23, b13,
                                             knots, pknots, vars, Xpts);
    }

    for (int m = 0; m < numGauss; m++) {
      for (int n = 0; n < numGauss; n++) {
        // Set the quadrature points
        double pt[2];
        pt[0] = gaussPts[n];
        pt[1] = gaussPts[m];

        // Calculate the shape functions
        shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, pt, Xpts);
        compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

        TacsScalar h = 0.0;
        if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
          h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, Xd,
                                Xdd);
        } else {
          const TacsScalar *axis = stiff->getRefAxis();
          h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi,
                                        normal_eta, axis, Xd, Xdd);
        }
        h = gaussWts[n] * gaussWts[m] * h;

        // Calculate the deformation at the current point...
        TacsScalar rot = 0.0;  // Difference between inplane/drilling rotation
        if (type == LINEAR) {
          linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal, normal_xi,
                             normal_eta);
          linear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, t, tx, ztx, normal,
                           normal_xi, normal_eta);
        }

        // Evaluate the strain interpolation at this point
        tying_interpolation<tying_order>(pt, N11, N22, N12, knots, pknots);
        add_tying_strain<tying_order>(strain, tx, g11, g22, g12, g23, g13, N11,
                                      N22, N12);
        add_tying_bmat<tying_order>(B, NUM_NODES, tx, b11, b22, b12, b23, b13,
                                    N11, N22, N12);

        // Compute the product of psi^{T}*B^{T}
        TacsScalar bpsi[NUM_STRESSES], bphi[NUM_STRESSES];
        memset(bpsi, 0, NUM_STRESSES * sizeof(TacsScalar));
        memset(bphi, 0, NUM_STRESSES * sizeof(TacsScalar));

        TacsScalar brpsi = 0.0, brphi = 0.0;
        TacsScalar *b = B, *dr = drot;
        const TacsScalar *ps = psi, *ph = phi;
        for (int i = 0; i < NUM_VARIABLES; i++) {
          bpsi[0] += ps[0] * b[0];
          bpsi[1] += ps[0] * b[1];
          bpsi[2] += ps[0] * b[2];
          bpsi[3] += ps[0] * b[3];
          bpsi[4] += ps[0] * b[4];
          bpsi[5] += ps[0] * b[5];
          bpsi[6] += ps[0] * b[6];
          bpsi[7] += ps[0] * b[7];

          bphi[0] += ph[0] * b[0];
          bphi[1] += ph[0] * b[1];
          bphi[2] += ph[0] * b[2];
          bphi[3] += ph[0] * b[3];
          bphi[4] += ph[0] * b[4];
          bphi[5] += ph[0] * b[5];
          bphi[6] += ph[0] * b[6];
          bphi[7] += ph[0] * b[7];

          brpsi += ps[0] * dr[0];
          brphi += ph[0] * dr[0];

          // Increment the pointers
          b += NUM_STRESSES;
          dr++;
          ps++;
          ph++;
        }

        for (int i = 0; i < NUM_STRESSES; i++) {
          bpsi[i] *= scale * h;
        }

        // Add the term: scale*psi^{T}*B^{T}*dC/dx*strain to the vector
        // dvSens - Note that this is much more efficient than computing
        // the terms component by component
        stiff->addStiffnessDVSens(pt, bpsi, bphi, scale * h * brpsi * brphi,
                                  fdvSens, dvLen);
      }
    }
  } else if (matType == MASS_MATRIX) {
    for (int m = 0; m < numGauss; m++) {
      for (int n = 0; n < numGauss; n++) {
        // Set the quadrature points
        double pt[2];
        pt[0] = gaussPts[n];
        pt[1] = gaussPts[m];

        // Evaluate the shape functions and the derivatives of the
        // position of the shell surface
        TacsScalar X[3], Xd[9], normal[3];
        double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
        shell_jacobian(order, X, Xd, N, Na, Nb, pt, Xpts);

        // Compute the normal to the surface
        Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
        Tensor::normalize3D(normal);
        for (int i = 0; i < 3; i++) {
          Xd[6 + i] = normal[i];
        }

        // Compute the determinant of the shell position
        TacsScalar h = FElibrary::jacobian3d(Xd);
        h = gaussWts[n] * gaussWts[m] * h;

        // Compute the nodal accelerations at the quadrature point
        TacsScalar upsi[6], uphi[6];
        upsi[0] = upsi[1] = upsi[2] = upsi[3] = upsi[4] = upsi[5] = 0.0;
        uphi[0] = uphi[1] = uphi[2] = uphi[3] = uphi[4] = uphi[5] = 0.0;

        double *ns = N;
        const TacsScalar *ps = psi, *ph = phi;
        for (int i = 0; i < NUM_NODES; i++) {
          upsi[0] += ns[0] * ps[0];
          upsi[1] += ns[0] * ps[1];
          upsi[2] += ns[0] * ps[2];
          upsi[3] += ns[0] * ps[3];
          upsi[4] += ns[0] * ps[4];
          upsi[5] += ns[0] * ps[5];

          uphi[0] += ns[0] * ph[0];
          uphi[1] += ns[0] * ph[1];
          uphi[2] += ns[0] * ph[2];
          uphi[3] += ns[0] * ph[3];
          uphi[4] += ns[0] * ph[4];
          uphi[5] += ns[0] * ph[5];
          ps += 6;
          ph += 6;
          ns++;
        }

        // Compute the weights on each component of the mass moments
        TacsScalar rho_scale[2];
        rho_scale[0] =
            scale * h *
            (upsi[0] * uphi[0] + upsi[1] * uphi[1] + upsi[2] * uphi[2]);

        TacsScalar tphi[3], tpsi[3];
        Tensor::crossProduct3D(tphi, &uphi[3], normal);
        Tensor::crossProduct3D(tpsi, &upsi[3], normal);
        rho_scale[1] =
            scale * h *
            (tphi[0] * tpsi[0] + tphi[1] * tpsi[1] + tphi[2] * tpsi[2]);

        // Add the result to the design variable vector
        stiff->addPointwiseMassDVSens(pt, rho_scale, fdvSens, dvLen);
      }
    }
  } else if (matType == GEOMETRIC_STIFFNESS_MATRIX && type == LINEAR) {
    // Geometric data
    TacsScalar X[3], Xd[9], Xdd[9];
    TacsScalar normal[3], normal_xi[3], normal_eta[3];

    // Transformation and the transformation derivative w.r.t. zeta
    TacsScalar t[9], tx[9], ztx[9];

    // Displacements/rotations and shape functions
    TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

    // Interpolations for the shear components
    double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

    // The interpolated tensorial shear components
    TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
    TacsScalar g13[NUM_G13], g23[NUM_G23];

    // The strain values
    TacsScalar strain[NUM_STRESSES];
    TacsScalar bstrain[NUM_STRESSES];

    // Compute the tying strain values
    const int is_linear = 1;
    compute_tying_strain<order, tying_order>(is_linear, g11, g22, g12, g23, g13,
                                             knots, pknots, vars, Xpts);

    for (int m = 0; m < numGauss; m++) {
      for (int n = 0; n < numGauss; n++) {
        double gpt[2];
        gpt[0] = gaussPts[n];
        gpt[1] = gaussPts[m];

        // Calculate the shape functions
        shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, gpt, Xpts);
        compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

        TacsScalar h = 0.0;
        if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
          h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, Xd,
                                Xdd);
        } else {
          const TacsScalar *axis = stiff->getRefAxis();
          h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi,
                                        normal_eta, axis, Xd, Xdd);
        }
        h = gaussWts[n] * gaussWts[m] * h;

        // Calculate the strain/bmat at the current point...
        TacsScalar rot;
        linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal, normal_xi,
                           normal_eta);

        // Evaluate the strain interpolation at this point
        tying_interpolation<tying_order>(gpt, N11, N22, N12, knots, pknots);
        add_tying_strain<tying_order>(strain, tx, g11, g22, g12, g23, g13, N11,
                                      N22, N12);

        // Compute the derivatives of the phi/psi vectors
        TacsScalar Upsi[NUM_DISPS], Udpsi[2 * NUM_DISPS];
        TacsScalar Uphi[NUM_DISPS], Udphi[2 * NUM_DISPS];
        compute_shell_Ud(NUM_NODES, Upsi, Udpsi, psi, N, Na, Nb);
        compute_shell_Ud(NUM_NODES, Uphi, Udphi, phi, N, Na, Nb);

        // Compute the inner product of the second derivatives of the
        // stiffness matrix with the psi and phi vectors
        inner_nonlinear_bend_bmat(bstrain, Upsi, Udpsi, Uphi, Udphi, t, tx, ztx,
                                  normal, normal_xi, normal_eta);
        add_nonlinear_tying_inner_nmat<order, tying_order>(
            bstrain, psi, phi, tx, N11, N22, N12, knots, pknots, Xpts);

        for (int i = 0; i < NUM_STRESSES; i++) {
          strain[i] *= scale * h;
        }

        // Add the term: scale*psi^{T}*B^{T}*dC/dx*strain to the vector
        // dvSens - Note that this is much more efficient than computing
        // the terms component by component
        stiff->addStiffnessDVSens(gpt, strain, bstrain, 0.0, fdvSens, dvLen);
      }
    }
  }
}

/*
  Add the derivative of the inner product of the given matrix type
  with the given psi and phi vectors with respect to the state
  variables.  This only makes sense for nonlinear Jacobian matrices
  such as the geometric stiffness matrix.

  input:
  matType:  the type of matrix
  scale:    the scaling factor applied to the derivative
  psi:      the left-multiplying vector
  phi:      the right-multiplying vector
  vars:     the element state variables
  Xpts:     the nodal locations for this element

  output:
  res:      the derivative of the inner product w.r.t. vars
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::getMatSVSensInnerProduct(
    ElementMatrixType matType, TacsScalar res[], const TacsScalar psi[],
    const TacsScalar phi[], const TacsScalar Xpts[], const TacsScalar vars[]) {
  if (matType == GEOMETRIC_STIFFNESS_MATRIX && type == LINEAR) {
    // Geometric data
    TacsScalar X[3], Xd[9], Xdd[9];
    TacsScalar normal[3], normal_xi[3], normal_eta[3];

    // Transformation and the transformation derivative w.r.t. zeta
    TacsScalar t[9], tx[9], ztx[9];

    // Displacements/rotations and shape functions
    TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

    // Interpolations for the shear components
    double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

    // The interpolated tensorial shear components
    TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
    TacsScalar g13[NUM_G13], g23[NUM_G23];

    // The derivatives of the displacement strain
    TacsScalar b11[3 * NUM_NODES * NUM_G11], b22[3 * NUM_NODES * NUM_G22];
    TacsScalar b12[3 * NUM_NODES * NUM_G12];
    TacsScalar b13[NUM_VARIABLES * NUM_G13], b23[NUM_VARIABLES * NUM_G23];

    // The strain values
    TacsScalar strain[NUM_STRESSES];
    TacsScalar bstrain[NUM_STRESSES], bstress[NUM_STRESSES];

    // The derivative of the strain w.r.t. the displacements
    TacsScalar B[NUM_STRESSES * NUM_VARIABLES];
    TacsScalar drot[NUM_VARIABLES];

    // Compute the tying strain values
    const int is_linear = 1;
    compute_tying_bmat<order, tying_order>(is_linear, g11, g22, g12, g23, g13,
                                           b11, b22, b12, b23, b13, knots,
                                           pknots, vars, Xpts);

    for (int m = 0; m < numGauss; m++) {
      for (int n = 0; n < numGauss; n++) {
        // Set the quadrature point
        double pt[2];
        pt[0] = gaussPts[n];
        pt[1] = gaussPts[m];

        // Evaluate the stiffness properties at the quadrature poitn
        TacsScalar At[6], Bt[6], Dt[6], Ats[3];
        stiff->getStiffness(pt, At, Bt, Dt, Ats);

        // Calculate the shape functions
        shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, pt, Xpts);
        compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

        TacsScalar h = 0.0;
        if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
          h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, Xd,
                                Xdd);
        } else {
          const TacsScalar *axis = stiff->getRefAxis();
          h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi,
                                        normal_eta, axis, Xd, Xdd);
        }
        h = gaussWts[n] * gaussWts[m] * h;

        // Calculate the strain/bmat at the current point...
        TacsScalar rot;
        linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal, normal_xi,
                           normal_eta);
        linear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, t, tx, ztx, normal,
                         normal_xi, normal_eta);

        // Evaluate the strain interpolation at this point
        tying_interpolation<tying_order>(pt, N11, N22, N12, knots, pknots);
        add_tying_strain<tying_order>(strain, tx, g11, g22, g12, g23, g13, N11,
                                      N22, N12);
        add_tying_bmat<tying_order>(B, NUM_NODES, tx, b11, b22, b12, b23, b13,
                                    N11, N22, N12);

        // Compute the derivatives of the phi/psi vectors
        TacsScalar Upsi[NUM_DISPS], Udpsi[2 * NUM_DISPS];
        TacsScalar Uphi[NUM_DISPS], Udphi[2 * NUM_DISPS];
        compute_shell_Ud(NUM_NODES, Upsi, Udpsi, psi, N, Na, Nb);
        compute_shell_Ud(NUM_NODES, Uphi, Udphi, phi, N, Na, Nb);

        // Compute the inner product of the second derivatives of the
        // stiffness matrix with the psi and phi vectors
        inner_nonlinear_bend_bmat(bstrain, Upsi, Udpsi, Uphi, Udphi, t, tx, ztx,
                                  normal, normal_xi, normal_eta);
        add_nonlinear_tying_inner_nmat<order, tying_order>(
            bstrain, psi, phi, tx, N11, N22, N12, knots, pknots, Xpts);

        // Compute the stress at the current Gauss point
        stiff->calculateStress(At, Bt, Dt, Ats, bstrain, bstress);

        for (int i = 0; i < NUM_NODES; i++) {
          for (int ii = 0; ii < NUM_DISPS; ii++) {
            int row = ii + NUM_DISPS * i;
            res[row] += h * strain_product(&B[NUM_STRESSES * row], bstress);
          }
        }
      }
    }
  }
}

/*
  Return the number of points in the specified quadrature scheme
*/
template <int order, int tying_order>
int MITCShell<order, tying_order>::getNumGaussPts() {
  return numGauss * numGauss;
}

/*
  Retrieve the Gauss points and weights for the given Gauss point
  number

  input:
  scheme: the Gauss quadrature scheme
  num:    the Gauss point index for the quadrature point

  returns: the Gauss quadrature weight for the given point

  output:
  pt:   the Gauss point for the given index
*/
template <int order, int tying_order>
double MITCShell<order, tying_order>::getGaussWtsPts(const int num,
                                                     double *pt) {
  int m = (int)(num / numGauss);
  int n = num % numGauss;
  pt[0] = gaussPts[n];
  pt[1] = gaussPts[m];

  return gaussWts[n] * gaussWts[m];
}

/*
  Evaluate the shape functions for this element at the specified point

  output:
  N:  the shape functions values evaluated at the parametric point pt

  input:
  pt: the parametric point within the element
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::getShapeFunctions(const double pt[],
                                                      double N[]) {
  double na[order], nb[order];
  FElibrary::lagrangeSF(na, pt[0], order);
  FElibrary::lagrangeSF(nb, pt[1], order);

  for (int j = 0; j < order; j++) {
    for (int i = 0; i < order; i++) {
      N[0] = na[i] * nb[j];
      N++;
    }
  }
}

/*
  Evaluate the determinant of the Jacobian for numerical integration

  returns: the determinant of the Jacobian

  input:
  pt: the parametric point within the element
  Xpts: the element nodes
*/
template <int order, int tying_order>
TacsScalar MITCShell<order, tying_order>::getDetJacobian(
    const double *pt, const TacsScalar Xpts[]) {
  TacsScalar X[3], Xd[9];
  TacsScalar normal[3];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];

  shell_jacobian(order, X, Xd, N, Na, Nb, pt, Xpts);
  Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
  Tensor::normalize3D(normal);

  for (int i = 0; i < 3; i++) {
    Xd[i + 6] = normal[i];
  }

  return FElibrary::jacobian3d(Xd);
}

/*
  Evaluate the derivative of the determinant of the Jacobian with respect
  to the element nodal locations

  output:
  hXptSens:  the derivative of the determinant w.r.t. the nodal locations

  returns: the determinant of the Jacobian

  input:
  pt: the parametric point within the element
  Xpts: the element nodes
*/
template <int order, int tying_order>
TacsScalar MITCShell<order, tying_order>::getDetJacobianXptSens(
    TacsScalar *hXptSens, const double *pt, const TacsScalar Xpts[]) {
  TacsScalar h = 0.0;
  TacsScalar X[3], Xd[9], XdSens[9];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];

  shell_jacobian(order, X, Xd, N, Na, Nb, pt, Xpts);

  for (int i = 0; i < NUM_NODES; i++) {
    for (int k = 0; k < 3; k++) {
      XdSens[0] = XdSens[1] = XdSens[2] = 0.0;
      XdSens[3] = XdSens[4] = XdSens[5] = 0.0;
      XdSens[k] = Na[i];
      XdSens[3 + k] = Nb[i];

      Tensor::crossProduct3DSens(&Xd[6], &XdSens[6], &Xd[0], &Xd[3], &XdSens[0],
                                 &XdSens[3]);
      TacsScalar snrm;
      Tensor::normalize3DSens(&snrm, &Xd[6], &XdSens[6]);
      h = FElibrary::jacobian3dSens(Xd, XdSens, &hXptSens[0]);
      hXptSens++;
    }
  }

  return h;
}

/*
  Evaluate the strain at the specified point using the provided
  set of variables

  output:
  strain:   the strain evaluate at the specific parametric point

  input:
  vars:     the element variable values
  Xpts:     the element nodal locations
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::getStrain(TacsScalar strain[],
                                              const double pt[],
                                              const TacsScalar Xpts[],
                                              const TacsScalar vars[]) {
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9];

  TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial shear components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the displacement strain at the tying points
  TacsScalar b11[3 * NUM_NODES * NUM_G11], b22[3 * NUM_NODES * NUM_G22];
  TacsScalar b12[3 * NUM_NODES * NUM_G12];
  TacsScalar b13[NUM_VARIABLES * NUM_G13], b23[NUM_VARIABLES * NUM_G23];

  shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, pt, Xpts);
  compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

  if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
    compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, Xd, Xdd);
  } else {
    const TacsScalar *axis = stiff->getRefAxis();
    compute_transform_refaxis(t, tx, ztx, normal, normal_xi, normal_eta, axis,
                              Xd, Xdd);
  }

  if (type == LARGE_ROTATION) {
    compute_lr_tying_bmat<order, tying_order>(g11, g22, g12, g23, g13, b11, b22,
                                              b12, b23, b13, knots, pknots,
                                              vars, Xpts);
  } else {
    compute_tying_strain<order, tying_order>(
        (type == LINEAR), g11, g22, g12, g23, g13, knots, pknots, vars, Xpts);
  }

  TacsScalar rot;
  if (type == LINEAR) {
    linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal, normal_xi,
                       normal_eta);
  } else if (type == NONLINEAR) {
    nonlinear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal, normal_xi,
                          normal_eta);
  } else {
    // Rotation matrix data
    TacsScalar C[9], Ct[27];

    // Compute the rotation matrices
    TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
    TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
    TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
    compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);

    compute_lr_tying_strain<order, tying_order>(g11, g22, g12, g23, g13, knots,
                                                pknots, vars, Xpts);
    large_rot_bend_strain(strain, U, Ud, C, Ct, t, tx, ztx, normal, normal_xi,
                          normal_eta);
  }

  // Evaluate the strain interpolation at this point
  tying_interpolation<tying_order>(pt, N11, N22, N12, knots, pknots);
  add_tying_strain<tying_order>(strain, tx, g11, g22, g12, g23, g13, N11, N22,
                                N12);
}

/*
  Compute the strain and the derivative of the strain with respect to
  the nodal locations.

  output:
  strain:  the strain evaluate at the pamametric point pt
  strainXptSens: the derivative of the strain w.r.t. the nodal locations

  input:
  pt:            the parametric point within the element
  vars:          the element variables
  Xpts:          the nodal locations
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::addStrainXptSens(
    TacsScalar fXptSens[], const double pt[], const TacsScalar scale,
    const TacsScalar strainSens[], const TacsScalar Xpts[],
    const TacsScalar vars[]) {
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];
  TacsScalar dnormal[9 * NUM_NODES], dnormal_xi[9 * NUM_NODES],
      dnormal_eta[9 * NUM_NODES];

  // The derivative of the determinant of the Jacobian and the
  // sensitivity of the drilling rotation
  TacsScalar dh[3 * NUM_NODES], srot[3 * NUM_NODES];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9];
  TacsScalar dt[27 * NUM_NODES], dtx[27 * NUM_NODES], dztx[27 * NUM_NODES];

  TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial strain components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the tensorial tying strain components
  TacsScalar dg11[3 * NUM_G11 * NUM_NODES], dg22[3 * NUM_G22 * NUM_NODES],
      dg12[3 * NUM_G12 * NUM_NODES];
  TacsScalar dg13[3 * NUM_G13 * NUM_NODES], dg23[3 * NUM_G23 * NUM_NODES];

  // The strain and the derivative of the strain w.r.t. nodes
  TacsScalar strain[NUM_STRESSES];
  TacsScalar strainXptSens[3 * NUM_STRESSES * NUM_NODES];

  // Zero out the strain sensitivity
  memset(strainXptSens, 0, 3 * NUM_STRESSES * NUM_NODES * sizeof(TacsScalar));

  // Evaluate the tying interpolation
  tying_interpolation<tying_order>(pt, N11, N22, N12, knots, pknots);

  // Calculate the shape functions
  shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, pt, Xpts);

  // Evaluate the displacements and their derivatives
  compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

  if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
    compute_transform_sens(dh, t, dt, tx, dtx, ztx, dztx, normal, dnormal,
                           normal_xi, dnormal_xi, normal_eta, dnormal_eta, Xd,
                           Xdd, Na, Nb, Naa, Nab, Nbb, NUM_NODES);
  } else {
    const TacsScalar *axis = stiff->getRefAxis();
    compute_transform_refaxis_sens(dh, t, dt, tx, dtx, ztx, dztx, normal,
                                   dnormal, normal_xi, dnormal_xi, normal_eta,
                                   dnormal_eta, axis, Xd, Xdd, Na, Nb, Naa, Nab,
                                   Nbb, NUM_NODES);
  }

  if (type == LARGE_ROTATION) {
    // Rotation matrix data
    TacsScalar C[9], Ct[27];

    // Evaluate the rotation matrices
    TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
    TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
    TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
    compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);

    compute_lr_tying_strain_sens<order, tying_order>(
        g11, g22, g12, g23, g13, dg11, dg22, dg12, dg23, dg13, knots, pknots,
        vars, Xpts);
    large_rot_bend_strain_sens(
        strain, strainXptSens, U, Ud, C, Ct, t, dt, tx, dtx, ztx, dztx, normal,
        dnormal, normal_xi, dnormal_xi, normal_eta, dnormal_eta, 3 * NUM_NODES);
  } else {
    compute_tying_strain_sens<order, tying_order>(
        (type == LINEAR), g11, g22, g12, g23, g13, dg11, dg22, dg12, dg23, dg13,
        knots, pknots, vars, Xpts);

    TacsScalar rot;
    if (type == LINEAR) {
      linear_bend_strain_sens(strain, strainXptSens, &rot, srot, U, Ud, t, dt,
                              tx, dtx, ztx, dztx, normal, dnormal, normal_xi,
                              dnormal_xi, normal_eta, dnormal_eta,
                              3 * NUM_NODES);
    } else {
      nonlinear_bend_strain_sens(strain, strainXptSens, &rot, srot, U, Ud, t,
                                 dt, tx, dtx, ztx, dztx, normal, dnormal,
                                 normal_xi, dnormal_xi, normal_eta, dnormal_eta,
                                 3 * NUM_NODES);
    }
  }

  // Evaluate the strain interpolation at this point
  add_tying_strain_sens<tying_order>(strain, strainXptSens, tx, dtx, g11, g22,
                                     g12, g23, g13, dg11, dg22, dg12, dg23,
                                     dg13, N11, N22, N12);

  // Add the product of the input sensitivity and the derivative of the
  // strain w.r.t. the node locations to the input vector
  TacsScalar *s = strainXptSens;
  for (int i = 0; i < 3 * NUM_NODES; i++) {
    for (int j = 0; j < NUM_STRESSES; j++, s++) {
      fXptSens[i] += scale * s[0] * strainSens[j];
    }
  }
}

/*
  Compute the derivative of the point-wise strain with respect to the
  element variables, multiply the derivative by a strain sensitivity
  vector and add the result, times a scalar multiple, to the outputt
  array. This can be used to evaluate the derivative of functions of
  the strain with respect to the state variables using the chain rule.

  output:
  elementSens: the output array - same length as the number of elem variables

  input:
  pt:          parametric point used to evaluate the derivative [-1, 1]^{2}
  strainSens:  the sensitivity of each strain component
  vars:        the element variables
  Xpts:        the element nodal locations
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::addStrainSVSens(
    TacsScalar strainSVSens[], const double pt[], const TacsScalar scale,
    const TacsScalar strainSens[], const TacsScalar Xpts[],
    const TacsScalar vars[]) {
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];
  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9];

  TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The tensorial components of the strain
  TacsScalar g11[NUM_G11], g22[NUM_G22];
  TacsScalar g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the displacement strain
  TacsScalar b11[3 * NUM_NODES * NUM_G11], b22[3 * NUM_NODES * NUM_G22];
  TacsScalar b12[3 * NUM_NODES * NUM_G12];
  TacsScalar b13[NUM_VARIABLES * NUM_G13], b23[NUM_VARIABLES * NUM_G23];

  TacsScalar dinplane_rot[NUM_VARIABLES];
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

  // Calculate the shape functions
  shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, pt, Xpts);
  compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

  if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
    compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, Xd, Xdd);
  } else {
    const TacsScalar *axis = stiff->getRefAxis();
    compute_transform_refaxis(t, tx, ztx, normal, normal_xi, normal_eta, axis,
                              Xd, Xdd);
  }

  if (type == LARGE_ROTATION) {
    // Rotational matrix data
    TacsScalar C[9], Ct[27], Ctt[54];

    // Compute the rotation matrices
    TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
    TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
    TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
    compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
    compute_2nd_rate_matrix(Ctt, c1, s1, c2, s2, c3, s3);

    compute_lr_tying_bmat<order, tying_order>(g11, g22, g12, g23, g13, b11, b22,
                                              b12, b23, b13, knots, pknots,
                                              vars, Xpts);
    large_rot_bend_bmat(B, NUM_NODES, N, Na, Nb, U, Ud, C, Ct, Ctt, t, tx, ztx,
                        normal, normal_xi, normal_eta);
  } else {
    compute_tying_bmat<order, tying_order>((type == LINEAR), g11, g22, g12, g23,
                                           g13, b11, b22, b12, b23, b13, knots,
                                           pknots, vars, Xpts);

    if (type == LINEAR) {
      linear_bend_bmat(B, dinplane_rot, NUM_NODES, N, Na, Nb, t, tx, ztx,
                       normal, normal_xi, normal_eta);
    } else {
      nonlinear_bend_bmat(B, dinplane_rot, NUM_NODES, N, Na, Nb, U, Ud, t, tx,
                          ztx, normal, normal_xi, normal_eta);
    }
  }

  tying_interpolation<tying_order>(pt, N11, N22, N12, knots, pknots);
  add_tying_bmat<tying_order>(B, NUM_NODES, tx, b11, b22, b12, b23, b13, N11,
                              N22, N12);

  for (int k = 0; k < NUM_VARIABLES; k++) {
    strainSVSens[k] += scale * strain_product(strainSens, &B[k * NUM_STRESSES]);
  }
}

/*
  Get the partition of unity shape functions and their derivatives
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::getPartUnityShapeFunctions(
    const double pt[], double N[], double Na[], double Nb[]) {
  N[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
  N[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
  N[2] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
  N[3] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);

  Na[0] = -0.25 * (1.0 - pt[1]);
  Na[1] = 0.25 * (1.0 - pt[1]);
  Na[2] = -0.25 * (1.0 + pt[1]);
  Na[3] = 0.25 * (1.0 + pt[1]);

  Nb[0] = -0.25 * (1.0 - pt[0]);
  Nb[1] = -0.25 * (1.0 + pt[0]);
  Nb[2] = 0.25 * (1.0 - pt[0]);
  Nb[3] = 0.25 * (1.0 + pt[0]);
}

/*
  Add the localized error term to the input.

  This localization is based on a partition of unity constraint that
  distributes the error to nodes.
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::addLocalizedError(
    double time, TacsScalar err[], const TacsScalar adjoint[],
    const TacsScalar Xpts[], const TacsScalar vars[]) {
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9];

  TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial shear components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the displacement strain
  TacsScalar b11[3 * NUM_NODES * NUM_G11], b22[3 * NUM_NODES * NUM_G22];
  TacsScalar b12[3 * NUM_NODES * NUM_G12];
  TacsScalar b13[NUM_VARIABLES * NUM_G13], b23[NUM_VARIABLES * NUM_G23];

  // The strain an rotation matrices
  TacsScalar strain[NUM_STRESSES], stress[NUM_STRESSES];
  TacsScalar drot[NUM_VARIABLES];
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

  // Determine the knots/pknots from the lower-order element
  const double *low_knots = NULL;
  const double *low_pknots = NULL;
  if (tying_order == 3) {
    low_knots = MITCShellFirstOrderKnots;
  } else {
    FElibrary::getGaussPtsWts(tying_order - 1, &low_knots, NULL);
  }
  FElibrary::getGaussPtsWts(tying_order - 2, &low_pknots, NULL);

  // Evaluate the strain/bmatrix at the tying points
  if (type == LARGE_ROTATION) {
    compute_lr_tying_bmat<order, tying_order>(g11, g22, g12, g23, g13, b11, b22,
                                              b12, b23, b13, knots, pknots,
                                              vars, Xpts);
    compute_lr_tying_strain<order, tying_order - 1>(
        g11, g22, g12, g23, g13, low_knots, low_pknots, vars, Xpts);
  } else {
    // Compute the B-mat contribution from the full residual
    compute_tying_bmat<order, tying_order>((type == LINEAR), g11, g22, g12, g23,
                                           g13, b11, b22, b12, b23, b13, knots,
                                           pknots, vars, Xpts);

    if (order == 3 && tying_order == 3) {
      TacsScalar Xpts2nd[12], vars2nd[24];
      for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
          for (int k = 0; k < 3; k++) {
            Xpts2nd[3 * (i + 2 * j) + k] = Xpts[3 * (2 * i + 6 * j) + k];
          }
          for (int k = 0; k < 6; k++) {
            vars2nd[6 * (i + 2 * j) + k] = vars[6 * (2 * i + 6 * j) + k];
          }
        }
      }

      // Compute the strain contribution from the lower-order element
      compute_tying_strain<order - 1, tying_order - 1>(
          (type == LINEAR), g11, g22, g12, g23, g13, low_knots, low_pknots,
          vars2nd, Xpts2nd);
    } else {
      // Compute the strain contribution from the lower-order element
      compute_tying_strain<order, tying_order - 1>((type == LINEAR), g11, g22,
                                                   g12, g23, g13, low_knots,
                                                   low_pknots, vars, Xpts);
    }
  }

  // Perform the integration over a larger domain
  const double *pts, *wts;
  int npts = FElibrary::getGaussPtsWts(order, &pts, &wts);

  for (int m = 0; m < npts; m++) {
    for (int n = 0; n < npts; n++) {
      // Set the quadrature point
      double pt[2];
      pt[0] = pts[n];
      pt[1] = pts[m];

      // Evaluate the stiffness at the parametric point within the
      // element
      TacsScalar At[6], Bt[6], Dt[6], Ats[3];
      TacsScalar kpenalty = stiff->getStiffness(pt, At, Bt, Dt, Ats);

      // Calculate the shape functions and the Jacobian/Hessian of the
      // shell position at the quadrature point
      shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, pt, Xpts);
      compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

      // Compute the transformation from the global coordinates to
      // local shell coordinates
      TacsScalar h = 0.0;
      if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
        h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, Xd,
                              Xdd);
      } else {
        const TacsScalar *axis = stiff->getRefAxis();
        h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi, normal_eta,
                                      axis, Xd, Xdd);
      }
      h = wts[n] * wts[m] * h;

      // Compute the strain and rotation at the qudrature point
      TacsScalar rot = 0.0;
      if (type == LINEAR) {
        linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal, normal_xi,
                           normal_eta);
        linear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, t, tx, ztx, normal,
                         normal_xi, normal_eta);
      } else if (type == NONLINEAR) {
        nonlinear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal,
                              normal_xi, normal_eta);
        nonlinear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, U, Ud, t, tx, ztx,
                            normal, normal_xi, normal_eta);
      } else {
        // Rotation matrix data
        TacsScalar C[9], Ct[27];

        // Compute the rotation matrices
        TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
        TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
        TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
        compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);

        // Calculate the deformation at the current point...
        large_rot_bend_strain(strain, U, Ud, C, Ct, t, tx, ztx, normal,
                              normal_xi, normal_eta);
      }

      // Evaluate the strain interpolation at this point
      tying_interpolation<tying_order>(pt, N11, N22, N12, knots, pknots);
      add_tying_bmat<tying_order>(B, NUM_NODES, tx, b11, b22, b12, b23, b13,
                                  N11, N22, N12);

      // Evaluate the strain interpolation using the lower-order
      // element tying points
      tying_interpolation<tying_order - 1>(pt, N11, N22, N12, low_knots,
                                           low_pknots);
      add_tying_strain<tying_order - 1>(strain, tx, g11, g22, g12, g23, g13,
                                        N11, N22, N12);

      // Compute the stress at the current Gauss point
      stiff->calculateStress(At, Bt, Dt, Ats, strain, stress);

      const TacsScalar *adj = adjoint;
      const TacsScalar *b = B;
      const TacsScalar *br = drot;

      // Compute the local product of the stress/strain
      TacsScalar product = 0.0;
      for (int i = 0; i < NUM_NODES; i++) {
        for (int ii = 0; ii < NUM_DISPS; ii++) {
          product += adj[ii] * h *
                     (strain_product(b, stress) + kpenalty * rot * br[0]);
          b += NUM_STRESSES;
          br++;
        }
        adj += NUM_DISPS;
      }

      // Evaluate the partition of unity constraint
      double Np[4], Npa[4], Npb[4];
      getPartUnityShapeFunctions(pt, Np, Npa, Npb);

      for (int node = 0; node < 4; node++) {
        // Add the result to the localized error
        err[(node % 2) * (order - 1) + (node / 2) * order * (order - 1)] +=
            Np[node] * product;
      }
    }
  }
}

/*
  Determine the number of nodes and elements for visualization
  generated by the data in this element. Note that higher-order
  elements are broken down into bi-linear elements for visualization.

  output:
  nelems:  the number of visualization elements
  nnodes:  the number of nodes used for visualization
  ncsr:    the number of entries in a CSR-type data structure used
  to store the connectivity
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::addOutputCount(int *nelems, int *nnodes,
                                                   int *ncsr) {
  *nelems += (order - 1) * (order - 1);
  *nnodes += order * order;
  *ncsr += 4 * (order - 1) * (order - 1);
}

/*
  Get the output data from this element and place it in a real
  array for visualization later. The values generated for visualization
  are determined by a bit-wise selection variable 'out_type' which is
  can be used to simultaneously write out different data. Note that this
  is why the bitwise operation & is used below.

  The output may consist of the following:
  - the nodal locations
  - the displacements and rotations
  - the strains or strains within the element
  - extra variables that are used for optimization

  output:
  data:     the data to write to the file (eventually)

  input:
  out_type: the bit-wise variable used to specify what data to generate
  vars:     the element variables
  Xpts:     the element nodal locations
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::getOutputData(unsigned int out_type,
                                                  double *data, int ld_data,
                                                  const TacsScalar Xpts[],
                                                  const TacsScalar vars[]) {
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9];

  TacsScalar U[NUM_DISPS], Ud[2 * NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial shear components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The stress and strain values
  TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];

  if (type == LARGE_ROTATION) {
    compute_lr_tying_strain<order, tying_order>(g11, g22, g12, g23, g13, knots,
                                                pknots, vars, Xpts);
  } else {
    compute_tying_strain<order, tying_order>(
        (type == LINEAR), g11, g22, g12, g23, g13, knots, pknots, vars, Xpts);
  }

  for (int m = 0; m < order; m++) {
    for (int n = 0; n < order; n++) {
      double pt[2];
      pt[0] = -1.0 + 2.0 / (order - 1.0) * n;
      pt[1] = -1.0 + 2.0 / (order - 1.0) * m;

      // Calculate the shape functions
      shell_hessian(order, X, Xd, Xdd, N, Na, Nb, Naa, Nab, Nbb, pt, Xpts);
      compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

      if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
        compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, Xd, Xdd);
      } else {
        const TacsScalar *axis = stiff->getRefAxis();
        compute_transform_refaxis(t, tx, ztx, normal, normal_xi, normal_eta,
                                  axis, Xd, Xdd);
      }

      // Evaluate the strain
      TacsScalar rot;
      if (type == LINEAR) {
        linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal, normal_xi,
                           normal_eta);
      } else if (type == NONLINEAR) {
        nonlinear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, normal,
                              normal_xi, normal_eta);
      } else {
        // Rotational matrix data
        TacsScalar C[9], Ct[27];

        // Compute the rotation matrices
        TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
        TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
        TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
        compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);

        // Evaluate the strain
        large_rot_bend_strain(strain, U, Ud, C, Ct, t, tx, ztx, normal,
                              normal_xi, normal_eta);
      }

      // Evaluate the strain interpolation at this point
      tying_interpolation<tying_order>(pt, N11, N22, N12, knots, pknots);
      add_tying_strain<tying_order>(strain, tx, g11, g22, g12, g23, g13, N11,
                                    N22, N12);

      int index = 0;
      int p = n + m * order;
      if (out_type & TACSElement::OUTPUT_NODES) {
        for (int k = 0; k < 3; k++) {
          data[index + k] = TacsRealPart(Xpts[3 * p + k]);
        }
        index += 3;
      }
      if (out_type & TACSElement::OUTPUT_DISPLACEMENTS) {
        for (int k = 0; k < NUM_DISPS; k++) {
          data[index + k] = TacsRealPart(vars[NUM_DISPS * p + k]);
        }
        index += NUM_DISPS;
      }
      if (out_type & TACSElement::OUTPUT_STRAINS) {
        for (int k = 0; k < NUM_STRESSES; k++) {
          data[index + k] = TacsRealPart(strain[k]);
        }
        index += NUM_STRESSES;
      }
      if (out_type & TACSElement::OUTPUT_STRESSES) {
        // Evaluate the stiffness at the current point
        // and then calculate the stress
        stiff->calculateStress(pt, strain, stress);

        for (int k = 0; k < NUM_STRESSES; k++) {
          data[index + k] = TacsRealPart(stress[k]);
        }
        index += NUM_STRESSES;
      }
      if (out_type & TACSElement::OUTPUT_EXTRAS) {
        // Compute the failure value
        TacsScalar lambda;
        stiff->failure(pt, strain, &lambda);
        data[index] = TacsRealPart(lambda);

        // Compute the buckling constraint value
        TacsScalar bval;
        stiff->buckling(strain, &bval);
        data[index + 1] = TacsRealPart(bval);

        data[index + 2] = TacsRealPart(stiff->getDVOutputValue(0, pt));
        data[index + 3] = TacsRealPart(stiff->getDVOutputValue(1, pt));

        index += NUM_EXTRAS;
      }
      if (out_type & TACSElement::OUTPUT_COORDINATES) {
        // t is the transform from the global coordinates to the
        // local coordinates.
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            data[index + 3 * i + j] = TacsRealPart(t[3 * i + j]);
          }
        }
        index += 9;
      }

      data += ld_data;
    }
  }
}

/*
  Get the element connectivity for visualization purposes. Since each
  element consists of a series of sub-elements used for visualization,
  we also need the connectivity of these visualization elements.

  output:
  con:  the connectivity of the local visualization elements contributed
  by this finite-element

  input:
  node:  the node offset number - so that this connectivity is more or
  less global
*/
template <int order, int tying_order>
void MITCShell<order, tying_order>::getOutputConnectivity(int *con, int node) {
  int p = 0;
  for (int m = 0; m < order - 1; m++) {
    for (int n = 0; n < order - 1; n++) {
      con[4 * p] = node + n + m * order;
      con[4 * p + 1] = node + n + 1 + m * order;
      con[4 * p + 2] = node + n + 1 + (m + 1) * order;
      con[4 * p + 3] = node + n + (m + 1) * order;
      p++;
    }
  }
}

#endif
