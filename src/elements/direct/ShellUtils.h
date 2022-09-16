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

#ifndef TACS_SHELL_UTILS_H
#define TACS_SHELL_UTILS_H

/*
  Common shell related computations.

  To evaluate the strain at the shell surface:

  1. Comute the transformation required at the current point
  2. Evaluate the strain/bmat/product of stress with
  the second derivative of the strain etc
*/

#include "TACSObject.h"
#include "TensorToolbox.h"

/*
  The 'natural' shell transformations. These take the first local
  coordinate axis to be the line formed by the parametrization of the
  element surface. This is always defined for non-degenerate surfaces
*/

TACS_BEGIN_NAMESPACE(shellutils)

/*
  The maximum number of nodes that can be used for the second-derivative
  computations.

  Note that you can increase this value - but more memory will be
  required - even for low-order elements. This value is for cubic
  quadrilateral shell elements. Probably more than sufficient in most
  cases.
*/
static const int MAX_NUM_NODES = 16;

/*
  Compute the transformations required for computing the strain in the
  local shell reference frame.

  input: Xd, Xdd
  output: t, tx, ztx, n, n_xi, n_eta, Xd

  Note that the first 6 entries of Xd must be the derivatives along
  the coordinate directions on input. On output, these are unmodified
  but the last 3 entries of Xd are defined as the normal direction.
*/
TacsScalar compute_transform(TacsScalar t[], TacsScalar tx[], TacsScalar ztx[],
                             TacsScalar n[], TacsScalar n_xi[],
                             TacsScalar n_eta[], TacsScalar Xd[],
                             const TacsScalar Xdd[]);

/*
  Compute the derivative of the transform w.r.t. the sensitivities
  of the input parameters.

  input: num_nodes, Xd, Xdd, Na, Nb, Naa, Nab, Nbb

  output: dh, t, tx, ztx, n, n_xi, n_eta, Xd
  output: the derivatives dt, dtx, dztx, dn, dn_xi, dn_eta

  Note that the number of derivatives produced by this code is
  3*num_nodes.
*/
TacsScalar compute_transform_sens(
    TacsScalar dh[], TacsScalar t[], TacsScalar dt[], TacsScalar tx[],
    TacsScalar dtx[], TacsScalar ztx[], TacsScalar dztx[], TacsScalar n[],
    TacsScalar dn[], TacsScalar n_xi[], TacsScalar dn_xi[], TacsScalar n_eta[],
    TacsScalar dn_eta[], TacsScalar Xd[], const TacsScalar Xdd[],
    const double Na[], const double Nb[], const double Naa[],
    const double Nab[], const double Nbb[], int num_nodes);

/*
  The 'reference axis' shell transformations. Take the first local
  coordinate axis as the projected direction of the 'axis' onto the
  direction. This is not defined when the reference axis is normal
  to the shell surface. This can be problematic, but allows for
  greater over the material properties.
*/

/*
  Compute the transformations required for computing the strain in
  the local shell reference frame. This uses an axis as the first
  direction projected onto the surface directions defined by Xd.

  input: axis, Xd, Xdd
  output: t, tx, ztx, n, n_xi, n_eta, Xd
*/
TacsScalar compute_transform_refaxis(TacsScalar t[], TacsScalar tx[],
                                     TacsScalar ztx[], TacsScalar n[],
                                     TacsScalar n_xi[], TacsScalar n_eta[],
                                     const TacsScalar axis[], TacsScalar Xd[],
                                     const TacsScalar Xdd[]);

/*
  Compute the derivative of the transform w.r.t. the sensitivities
  of the input parameters
*/
TacsScalar compute_transform_refaxis_sens(
    TacsScalar dh[], TacsScalar t[], TacsScalar dt[], TacsScalar tx[],
    TacsScalar dtx[], TacsScalar ztx[], TacsScalar dztx[], TacsScalar n[],
    TacsScalar dn[], TacsScalar n_xi[], TacsScalar dn_xi[], TacsScalar n_eta[],
    TacsScalar dn_eta[], const TacsScalar axis[], TacsScalar Xd[],
    const TacsScalar Xdd[], const double Na[], const double Nb[],
    const double Naa[], const double Nab[], const double Nbb[], int num_nodes);

// Linear implementation
// ---------------------

/*
  Compute the linear strain associated with a shell element

  output:
  strain: the in-plane, bending and shear strain in the local axis
  rotz: the rotation about the unit normal

  input:
  Ux: the displacements at the current point
  Uxd: the derivative of the displacements at the current point
  t: the transformation
  tx: the transformation times the Jacobian
  ztx: the derivative of the transformation times the Jacobian
  n: the unit normal
  n_xi: the derivative of the normal along the xi direction
  n_eta: the derivative of the normal along the eta direction
*/
void linear_strain(TacsScalar strain[], TacsScalar *rotz, const TacsScalar Ux[],
                   const TacsScalar Uxd[], const TacsScalar t[],
                   const TacsScalar tx[], const TacsScalar ztx[],
                   const TacsScalar n[], const TacsScalar n_xi[],
                   const TacsScalar n_eta[]);

/*
  Find the derivative of the strain with respect to the displacements

  output:
  B: the derivative of the strain w.r.t. the nodal displacements
  rotz: the rotation about the unit normal

  input:
  num_point: the number of nodes
  N, Na, Nb: the shape functions and their derivatives
  Ux: the displacements at the current point
  Uxd: the derivative of the displacements at the current point
  t: the transformation
  tx: the transformation times the Jacobian
  ztx: the derivative of the transformation times the Jacobian
  n: the unit normal
  n_xi: the derivative of the normal along the xi direction
  n_eta: the derivative of the normal along the eta direction
*/
void linear_bmat(TacsScalar B[], TacsScalar rotz[], int num_points,
                 const double N[], const double Na[], const double Nb[],
                 const TacsScalar t[], const TacsScalar tx[],
                 const TacsScalar ztx[], const TacsScalar n[],
                 const TacsScalar n_xi[], const TacsScalar n_eta[]);

/*
  Compute the linear strain associated with a shell element

  output:
  B: the derivative of the strain w.r.t. the nodal displacements
  rotz: the rotation about the unit normal

  input:
  num_point: the number of nodes
  N, Na, Nb: the shape functions and their derivatives
  Ux: the displacements at the current point
  Uxd: the derivative of the displacements at the current point
  t: the transformation
  dt: the derivative of the transformation
  tx: the transformation times the Jacobian
  dtx: the deriative of tx
  ztx: the derivative of the transformation times the Jacobian
  dztx: the derivative of the ztx
  n: the unit normal
  dn: the derivative of the normal
  n_xi: the derivative of the normal along the xi direction
  dn_xi: the derivative of n_xi
  n_eta: the derivative of the normal along the eta direction
  dn_eta: the derivative of the n_eta
  num_components: the number of gradient components to compute
*/
void linear_strain_sens(TacsScalar strain[], TacsScalar dstrain[],
                        TacsScalar *rotz, TacsScalar *drotz,
                        const TacsScalar Ux[], const TacsScalar Uxd[],
                        const TacsScalar t[], const TacsScalar dt[],
                        const TacsScalar tx[], const TacsScalar dtx[],
                        const TacsScalar ztx[], const TacsScalar dztx[],
                        const TacsScalar n[], const TacsScalar dn[],
                        const TacsScalar n_xi[], const TacsScalar dn_xi[],
                        const TacsScalar n_eta[], const TacsScalar dn_eta[],
                        int num_components);

/*
  Add the sensitivity of B to the residual

  res = stress_scale * d(B*stress)/dx + rot_scale * d(d rot)/dx

  output:
  res: the product

  input:
  stress_scale: scalar to multiply result by
  rot_scale: scale the rotation by this
  stress: the stress-weighting factors - to multiply by dB/dx
  num_point: the number of nodes
  N, Na, Nb: the shape functions and their derivatives
  Ux: the displacements at the current point
  Uxd: the derivative of the displacements at the current point
  t: the transformation
  dt: the derivative of the transformation
  tx: the transformation times the Jacobian
  dtx: the deriative of tx
  ztx: the derivative of the transformation times the Jacobian
  dztx: the derivative of the ztx
  n: the unit normal
  dn: the derivative of the normal
  n_xi: the derivative of the normal along the xi direction
  dn_xi: the derivative of n_xi
  n_eta: the derivative of the normal along the eta direction
  dn_eta: the derivative of the n_eta
  num_components: the number of gradient components to compute
*/
void add_linear_bmat_sens(
    TacsScalar res[], int num_points, const TacsScalar stress_scale,
    const TacsScalar rot_scale, const TacsScalar stress[], const double N[],
    const double Na[], const double Nb[], const TacsScalar t[],
    const TacsScalar dt[], const TacsScalar tx[], const TacsScalar dtx[],
    const TacsScalar ztx[], const TacsScalar dztx[], const TacsScalar n[],
    const TacsScalar dn[], const TacsScalar n_xi[], const TacsScalar dn_xi[],
    const TacsScalar n_eta[], const TacsScalar dn_eta[], int num_components);

/*
  Compute the B-matrix and the derivative of the Bmatrix with
  respect to a single input perturbation.

  output:
  B: the B-matrix
  dB: dB/dx
  rotz: the derivative of the normal rotation w.r.t. the displacements
  drotz: the derivative of the rotz

  input:
  stress_scale: scalar to multiply result by
  rot_scale: scale the rotation by this
  stress: the stress-weighting factors - to multiply by dB/dx
  num_point: the number of nodes
  N, Na, Nb: the shape functions and their derivatives
  Ux: the displacements at the current point
  Uxd: the derivative of the displacements at the current point
  t: the transformation
  dt: the derivative of the transformation
  tx: the transformation times the Jacobian
  dtx: the deriative of tx
  ztx: the derivative of the transformation times the Jacobian
  dztx: the derivative of the ztx
  n: the unit normal
  dn: the derivative of the normal
  n_xi: the derivative of the normal along the xi direction
  dn_xi: the derivative of n_xi
  n_eta: the derivative of the normal along the eta direction
  dn_eta: the derivative of the n_eta
  num_components: the number of gradient components to compute
*/
void linear_bmat_sens(TacsScalar B[], TacsScalar dB[], TacsScalar rotz[],
                      TacsScalar drotz[], int num_points, const double N[],
                      const double Na[], const double Nb[],
                      const TacsScalar t[], const TacsScalar dt[],
                      const TacsScalar tx[], const TacsScalar dtx[],
                      const TacsScalar ztx[], const TacsScalar dztx[],
                      const TacsScalar n[], const TacsScalar dn[],
                      const TacsScalar n_xi[], const TacsScalar dn_xi[],
                      const TacsScalar n_eta[], const TacsScalar dn_eta[]);

// Nonlinear implementation
// ------------------------

/*
  Compute the nonlinear strain associated with a shell element at a point.
  This includes all nonlinear contributions in the membrane, bending and
  shear terms.
*/
void nonlinear_strain(TacsScalar strain[], TacsScalar *rotz,
                      const TacsScalar Ux[], const TacsScalar Uxd[],
                      const TacsScalar t[], const TacsScalar tx[],
                      const TacsScalar ztx[], const TacsScalar n[],
                      const TacsScalar n_xi[], const TacsScalar n_eta[]);

/*
  Find the derivative of the strain with respect to the displacements.
  Compute the B-matrix - a non-constant matrix since we're dealing with
  a nonlinear strain expression.
*/
void nonlinear_bmat(TacsScalar B[], TacsScalar rotz[], int num_points,
                    const double N[], const double Na[], const double Nb[],
                    const TacsScalar Ux[], const TacsScalar Uxd[],
                    const TacsScalar t[], const TacsScalar tx[],
                    const TacsScalar ztx[], const TacsScalar n[],
                    const TacsScalar n_xi[], const TacsScalar n_eta[]);

/*
  Compute the second derivative of the strain w.r.t. the product of
  the strain and stress.

  This is only strictly required for the nonlinear problems - but may
  be used for linear problems to approximate the geometric stiffness
  matrix.
*/
void nonlinear_stress_bmat(TacsScalar matrix[], int num_points,
                           TacsScalar scale, const TacsScalar stress[],
                           const double N[], const double Na[],
                           const double Nb[], const TacsScalar t[],
                           const TacsScalar tx[], const TacsScalar ztx[],
                           const TacsScalar n[], const TacsScalar n_xi[],
                           const TacsScalar n_eta[]);

/*
  Compute the sensitivity of the strain w.r.t. the nodal locations
  (typically).  In this function the derivatives of the input
  parameters (dt, dtx, dztx, etc.)  are stored in column-major format
  - stacked on top of one another.
*/
void nonlinear_strain_sens(TacsScalar strain[], TacsScalar dstrain[],
                           TacsScalar *rotz, TacsScalar *drotz,
                           const TacsScalar Ux[], const TacsScalar Uxd[],
                           const TacsScalar t[], const TacsScalar dt[],
                           const TacsScalar tx[], const TacsScalar dtx[],
                           const TacsScalar ztx[], const TacsScalar dztx[],
                           const TacsScalar n[], const TacsScalar dn[],
                           const TacsScalar n_xi[], const TacsScalar dn_xi[],
                           const TacsScalar n_eta[], const TacsScalar dn_eta[],
                           int num_components);

/*
  Add the sensitivity of the derivative of the strain w.r.t. the
  displacements (the B-matrix) to res.  Add the derivative:

  res += stress_scale*(stress^{T}B)/dx + rot_scale*d(rot)/dx
*/
void add_nonlinear_bmat_sens(
    TacsScalar res[], int num_points, const TacsScalar stress_scale,
    const TacsScalar rot_scale, const TacsScalar stress[], const double N[],
    const double Na[], const double Nb[], const TacsScalar Ux[],
    const TacsScalar Uxd[], const TacsScalar t[], const TacsScalar dt[],
    const TacsScalar tx[], const TacsScalar dtx[], const TacsScalar ztx[],
    const TacsScalar dztx[], const TacsScalar n[], const TacsScalar dn[],
    const TacsScalar n_xi[], const TacsScalar dn_xi[], const TacsScalar n_eta[],
    const TacsScalar dn_eta[], int num_components);

/*
  Add the sensitivity of the stress-strain product to a portion of the matrix
*/
void nonlinear_stress_bmat_sens(
    TacsScalar matrix[], int num_points, TacsScalar scale, TacsScalar sscale,
    const TacsScalar stress[], const TacsScalar sstress[], const double N[],
    const double Na[], const double Nb[], const TacsScalar t[],
    const TacsScalar dt[], const TacsScalar tx[], const TacsScalar dtx[],
    const TacsScalar ztx[], const TacsScalar dztx[], const TacsScalar n[],
    const TacsScalar dn[], const TacsScalar n_xi[], const TacsScalar dn_xi[],
    const TacsScalar n_eta[], const TacsScalar dn_eta[]);

/*
  For the mixed interpolation of tensorial components (MITC) element,
  the bending strain terms are computed based on the values of the
  displacements directly, while the in-plane and shear terms (the
  components susceptible to shear locking) are computed based on an
  interpolation of the tensorial components of the strain.

  The following functions add the bending-only terms to the strain,
  B-matrix and various design sensitivity terms that are required.
  These are essentially the same routines as above without the
  in-plane components, but are duplicated for efficiency.
*/

/*
  Compute the linear strain associated with a shell element
*/
void linear_bend_strain(TacsScalar strain[], TacsScalar *rotz,
                        const TacsScalar Ux[], const TacsScalar Uxd[],
                        const TacsScalar t[], const TacsScalar tx[],
                        const TacsScalar ztx[], const TacsScalar n[],
                        const TacsScalar n_xi[], const TacsScalar n_eta[]);

/*
  Find the derivative of the strain with respect to the displacements
*/
void linear_bend_bmat(TacsScalar B[], TacsScalar rotz[], int num_points,
                      const double N[], const double Na[], const double Nb[],
                      const TacsScalar t[], const TacsScalar tx[],
                      const TacsScalar ztx[], const TacsScalar n[],
                      const TacsScalar n_xi[], const TacsScalar n_eta[]);

/*
  Compute the linear strain associated with a shell element
*/
void linear_bend_strain_sens(TacsScalar strain[], TacsScalar dstrain[],
                             TacsScalar *rotz, TacsScalar *drotz,
                             const TacsScalar Ux[], const TacsScalar Uxd[],
                             const TacsScalar t[], const TacsScalar dt[],
                             const TacsScalar tx[], const TacsScalar dtx[],
                             const TacsScalar ztx[], const TacsScalar dztx[],
                             const TacsScalar n[], const TacsScalar dn[],
                             const TacsScalar n_xi[], const TacsScalar dn_xi[],
                             const TacsScalar n_eta[],
                             const TacsScalar dn_eta[], int num_components);

/*
  Compute the bending moment strain
*/
void add_linear_bend_bmat_sens(
    TacsScalar fXptSens[], const TacsScalar psi[], int num_points,
    const TacsScalar stress_scale, const TacsScalar rot_scale,
    const TacsScalar stress[], const double N[], const double Na[],
    const double Nb[], const TacsScalar t[], const TacsScalar dt[],
    const TacsScalar tx[], const TacsScalar dtx[], const TacsScalar ztx[],
    const TacsScalar dztx[], const TacsScalar n[], const TacsScalar dn[],
    const TacsScalar n_xi[], const TacsScalar dn_xi[], const TacsScalar n_eta[],
    const TacsScalar dn_eta[], int num_components);

// Nonlinear implementations of the bending-only strain
// ----------------------------------------------------

/*
  Compute the linear strain associated with a shell element
*/
void nonlinear_bend_strain(TacsScalar strain[], TacsScalar *rotz,
                           const TacsScalar Ux[], const TacsScalar Uxd[],
                           const TacsScalar t[], const TacsScalar tx[],
                           const TacsScalar ztx[], const TacsScalar n[],
                           const TacsScalar n_xi[], const TacsScalar n_eta[]);

/*
  Find the derivative of the strain with respect to the displacements
*/
void nonlinear_bend_bmat(TacsScalar B[], TacsScalar rotz[], int num_points,
                         const double N[], const double Na[], const double Nb[],
                         const TacsScalar Ux[], const TacsScalar Uxd[],
                         const TacsScalar t[], const TacsScalar tx[],
                         const TacsScalar ztx[], const TacsScalar n[],
                         const TacsScalar n_xi[], const TacsScalar n_eta[]);

/*
  Compute the second derivative of the strain w.r.t. the product of
  the strain and stress.

  This is only strictly required for the nonlinear problems - but may
  be used for linear problems to approximate the geometric stiffness
  matrix.
*/
void nonlinear_bend_stress_bmat(TacsScalar matrix[], int num_points,
                                TacsScalar scale, const TacsScalar stress[],
                                const double N[], const double Na[],
                                const double Nb[], const TacsScalar t[],
                                const TacsScalar tx[], const TacsScalar ztx[],
                                const TacsScalar n[], const TacsScalar n_xi[],
                                const TacsScalar n_eta[]);

/*
  Compute the second derivatives of the inner product of the strain
  with respect to
*/
void inner_nonlinear_bend_bmat(TacsScalar bstrain[], const TacsScalar Uxpsi[],
                               const TacsScalar Uxdpsi[],
                               const TacsScalar Uxphi[],
                               const TacsScalar Uxdphi[], const TacsScalar t[],
                               const TacsScalar tx[], const TacsScalar ztx[],
                               const TacsScalar n[], const TacsScalar n_xi[],
                               const TacsScalar n_eta[]);
/*
  Compute the sensitivity of the strain
*/
void nonlinear_bend_strain_sens(
    TacsScalar strain[], TacsScalar dstrain[], TacsScalar *rotz,
    TacsScalar *drotz, const TacsScalar Ux[], const TacsScalar Uxd[],
    const TacsScalar t[], const TacsScalar dt[], const TacsScalar tx[],
    const TacsScalar dtx[], const TacsScalar ztx[], const TacsScalar dztx[],
    const TacsScalar n[], const TacsScalar dn[], const TacsScalar n_xi[],
    const TacsScalar dn_xi[], const TacsScalar n_eta[],
    const TacsScalar dn_eta[], int num_components);

/*
  Compute the sensitivity of the nonlinear bmat
*/
void add_nonlinear_bend_bmat_sens(
    TacsScalar fXptSens[], const TacsScalar psi[], int num_points,
    const TacsScalar stress_scale, const TacsScalar rot_scale,
    const TacsScalar stress[], const double N[], const double Na[],
    const double Nb[], const TacsScalar Ux[], const TacsScalar Uxd[],
    const TacsScalar t[], const TacsScalar dt[], const TacsScalar tx[],
    const TacsScalar dtx[], const TacsScalar ztx[], const TacsScalar dztx[],
    const TacsScalar n[], const TacsScalar dn[], const TacsScalar n_xi[],
    const TacsScalar dn_xi[], const TacsScalar n_eta[],
    const TacsScalar dn_eta[], int num_components);

/*
  Compute the displacements at a parametric point within the element
*/
void compute_shell_U(const int num_nodes, TacsScalar U[],
                     const TacsScalar vars[], const double N[]);
/*
  Compute the displacement and the derivative of the displacement
  along the parametric directions

  input:
  num_nodes: the number of nodes
  vars: the x,y,z displacements and rotations for each node
  N: the shape functions
  Na: the derivative of the shape functions along the first
  parametric direction
  Nb: the derivative of the shape functions along the second
  parametric direction

  output:
  U: the displacements and rotations at the parametric point
  Ud: the derivative of the displacements along the parametric directions
*/
void compute_shell_Ud(const int num_nodes, TacsScalar U[], TacsScalar Ud[],
                      const TacsScalar vars[], const double N[],
                      const double Na[], const double Nb[]);

/*
  Compute the position and first derivatives along the shell surface.
  At the same time, compute the shape functions.

  The shell surface is based on the shape functions and nodal
  locations.

  input:
  gpt: the current Gauss point
  Xpts: the nodal locations for the

  output:
  X, Xd: position and first derivative of the surface
  N: the shape functions
  Na, Nb: the first derivatives of the shape functions
*/
void shell_jacobian(const int order, TacsScalar X[], TacsScalar Xd[],
                    double N[], double Na[], double Nb[], const double gpt[],
                    const TacsScalar Xpts[]);

/*
  Compute the postion as well as the first and second derivatives
  along the surface of the shell.

  At the same time compute the shape functions, and their derivatives,
  at the current Gauss point gpt. Here Xpts is a vector of the nodal
  locations.

  input:
  gpt: the current Gauss point
  Xpts: the nodal locations for the

  output:
  X, Xd, Xdd: position, first and second derivatives of the surface
  N: the shape functions
  Na, Nb: the first derivatives of the shape functions
  Naa, Nab, Nbb: the second derivatives of the shape functions
*/
void shell_hessian(const int order, TacsScalar X[], TacsScalar Xd[],
                   TacsScalar Xdd[], double N[], double Na[], double Nb[],
                   double Naa[], double Nab[], double Nbb[], const double gpt[],
                   const TacsScalar Xpts[]);

/*
  Compute only the terms/derivatives required for computing the
  in-plane tensorial strain components.

  input:
  order:     order of the interpolation
  pt:        the parametric point to evaluate the shape functions
  Xpts:      the nodal locations
  vars:      the state varaible values

  output:
  Urot: the rotations at the mid-surface
  Ud: the derivatives of the displacements at the mid-surface
  Xd: the derivatives of the nodal locations along the parametric directions
  N, Na, Nb:  the tensor product shape functions and their derivatives
*/
void compute_tensorial_components(const int order, TacsScalar Urot[],
                                  TacsScalar Ud[], TacsScalar Xd[], double N[],
                                  double Na[], double Nb[], const double pt[],
                                  const TacsScalar Xpts[],
                                  const TacsScalar vars[]);

/*
  Compute only the terms and derivatives required to compute the
  in-plane and out-of-plane shear tensorial strain components for a
  b-spline shell.

  input:
  order:      the order of the b-spline in both u/v <= 8
  tu, tv:     the knot vector in the u/v directions
  intu, intv: the knot intervals in the u/v directions
  pt:         the parametric point to evaluate the shape functions
  Xpts:       the nodal locations
  vars:       the state varaible values

  output:
  Urot: the rotations at the mid-surface
  Ud: the derivatives of the displacements at the mid-surface
  Xd: the derivatives of the nodal locations along the parametric directions
  N, Na, Nb:  the tensor product shape functions and their derivatives
*/
void bspline_tensorial_components(const int order, TacsScalar Urot[],
                                  TacsScalar Ud[], TacsScalar Xd[], double N[],
                                  double Na[], double Nb[], const double pt[],
                                  int intu, int intv, const double tu[],
                                  const double tv[], const TacsScalar Xpts[],
                                  const TacsScalar vars[]);

/*
  Evaluate the tying-interpolation at the appropriate points
*/
template <int tying_order>
void tying_interpolation(const double pt[], double N11[], double N22[],
                         double N12[], const double knots[],
                         const double pknots[]) {
  double na[tying_order], nb[tying_order];
  FElibrary::lagrangeSFKnots(na, pt[0], knots, tying_order);
  FElibrary::lagrangeSFKnots(nb, pt[1], pknots, tying_order - 1);
  for (int j = 0; j < tying_order - 1; j++) {
    for (int i = 0; i < tying_order; i++) {
      N22[i + j * tying_order] = na[i] * nb[j];
    }
  }

  FElibrary::lagrangeSFKnots(na, pt[0], pknots, tying_order - 1);
  for (int j = 0; j < tying_order - 1; j++) {
    for (int i = 0; i < tying_order - 1; i++) {
      N12[i + j * (tying_order - 1)] = na[i] * nb[j];
    }
  }

  FElibrary::lagrangeSFKnots(nb, pt[1], knots, tying_order);
  for (int j = 0; j < tying_order; j++) {
    for (int i = 0; i < tying_order - 1; i++) {
      N11[i + j * (tying_order - 1)] = na[i] * nb[j];
    }
  }
}

/*
  Evaluate the tensorial displacement-based strain at the tying points.

  The tensorial strain is given as follows,

  E_t = U_{x,xi}*X_{,xi}^{T} + X_{,xi}*U_{x,xi}^{T}

  X_{,xi} = Xd^{T}
  U_{x,xi} = [ Ud, r ]

  input:
  knots: the 'order'th (ie 2, 3 or 4) Gauss-points
  pknots: the 'order-1'th Gauss-points
  vars: the displacements and rotations at the nodal locations
  Xpts: the nodal locations for the element

  output:
  g11: the in-plane normal tensorial strain in the 1 direction
  g22: the in-plane normal tensorial strain in the 2 direction
  g12: the in-plane tensorial shear strain

  g23: the out-of-plane tensorial shear strain in the 2-direction
  g13: the out-of-plane tensorial shear strain in the 1-direction
*/
template <int order, int tying_order>
void compute_tying_strain(const int is_linear, TacsScalar g11[],
                          TacsScalar g22[], TacsScalar g12[], TacsScalar g23[],
                          TacsScalar g13[], const double knots[],
                          const double pknots[], const TacsScalar vars[],
                          const TacsScalar Xpts[]) {
  double N[order * order], Na[order * order], Nb[order * order];

  // Evaluate g22 and g23
  if (is_linear) {
    for (int m = 0; m < tying_order - 1; m++) {
      for (int n = 0; n < tying_order; n++) {
        double pt[2];
        pt[0] = knots[n];
        pt[1] = pknots[m];

        TacsScalar Urot[3], Ud[6], Xd[6];
        compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                     vars);

        // Caculate the normal
        TacsScalar normal[3];
        Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
        Tensor::normalize3D(normal);

        // Calculate the rotation
        TacsScalar r[3];
        Tensor::crossProduct3D(r, Urot, normal);

        g22[0] = Xd[3] * Ud[1] + Xd[4] * Ud[3] + Xd[5] * Ud[5];
        g23[0] =
            0.5 * (Xd[3] * r[0] + Xd[4] * r[1] + Xd[5] * r[2] +
                   normal[0] * Ud[1] + normal[1] * Ud[3] + normal[2] * Ud[5]);
        g22++;
        g23++;
      }
    }
  } else {
    for (int m = 0; m < tying_order - 1; m++) {
      for (int n = 0; n < tying_order; n++) {
        double pt[2];
        pt[0] = knots[n];
        pt[1] = pknots[m];

        TacsScalar Urot[3], Ud[6], Xd[6];
        compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                     vars);

        // Caculate the normal
        TacsScalar normal[3];
        Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
        Tensor::normalize3D(normal);

        // Calculate the rotation
        TacsScalar r[3];
        Tensor::crossProduct3D(r, Urot, normal);

        g22[0] = (Xd[3] * Ud[1] + Xd[4] * Ud[3] + Xd[5] * Ud[5] +
                  0.5 * (Ud[1] * Ud[1] + Ud[3] * Ud[3] + Ud[5] * Ud[5]));
        g23[0] =
            0.5 * (Xd[3] * r[0] + Xd[4] * r[1] + Xd[5] * r[2] +
                   normal[0] * Ud[1] + normal[1] * Ud[3] + normal[2] * Ud[5]);
        g22++;
        g23++;
      }
    }
  }

  // Evaluate g12
  if (is_linear) {
    for (int m = 0; m < tying_order - 1; m++) {
      for (int n = 0; n < tying_order - 1; n++) {
        double pt[2];
        pt[0] = pknots[n];
        pt[1] = pknots[m];

        TacsScalar Urot[3], Ud[6], Xd[6];
        compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                     vars);

        g12[0] = 0.5 * (Xd[0] * Ud[1] + Xd[1] * Ud[3] + Xd[2] * Ud[5] +
                        Xd[3] * Ud[0] + Xd[4] * Ud[2] + Xd[5] * Ud[4]);
        g12++;
      }
    }
  } else {
    for (int m = 0; m < tying_order - 1; m++) {
      for (int n = 0; n < tying_order - 1; n++) {
        double pt[2];
        pt[0] = pknots[n];
        pt[1] = pknots[m];

        TacsScalar Urot[3], Ud[6], Xd[6];
        compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                     vars);

        g12[0] = 0.5 * (Xd[0] * Ud[1] + Xd[1] * Ud[3] + Xd[2] * Ud[5] +
                        Xd[3] * Ud[0] + Xd[4] * Ud[2] + Xd[5] * Ud[4] +
                        Ud[0] * Ud[1] + Ud[2] * Ud[3] + Ud[4] * Ud[5]);
        g12++;
      }
    }
  }

  // Evaluate g11 and g13
  if (is_linear) {
    for (int m = 0; m < tying_order; m++) {
      for (int n = 0; n < tying_order - 1; n++) {
        double pt[2];
        pt[0] = pknots[n];
        pt[1] = knots[m];

        TacsScalar Urot[3], Ud[6], Xd[6];
        compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                     vars);

        // Caculate the normal
        TacsScalar normal[3];
        Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
        Tensor::normalize3D(normal);

        // Calculate the rotation
        TacsScalar r[3];
        Tensor::crossProduct3D(r, Urot, normal);

        g11[0] = Xd[0] * Ud[0] + Xd[1] * Ud[2] + Xd[2] * Ud[4];
        g13[0] =
            0.5 * (Xd[0] * r[0] + Xd[1] * r[1] + Xd[2] * r[2] +
                   normal[0] * Ud[0] + normal[1] * Ud[2] + normal[2] * Ud[4]);
        g11++;
        g13++;
      }
    }
  } else {
    for (int m = 0; m < tying_order; m++) {
      for (int n = 0; n < tying_order - 1; n++) {
        double pt[2];
        pt[0] = pknots[n];
        pt[1] = knots[m];

        TacsScalar Urot[3], Ud[6], Xd[6];
        compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                     vars);

        // Caculate the normal
        TacsScalar normal[3];
        Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
        Tensor::normalize3D(normal);

        // Calculate the rotation
        TacsScalar r[3];
        Tensor::crossProduct3D(r, Urot, normal);

        g11[0] = (Xd[0] * Ud[0] + Xd[1] * Ud[2] + Xd[2] * Ud[4] +
                  0.5 * (Ud[0] * Ud[0] + Ud[2] * Ud[2] + Ud[4] * Ud[4]));
        g13[0] =
            0.5 * (Xd[0] * r[0] + Xd[1] * r[1] + Xd[2] * r[2] +
                   normal[0] * Ud[0] + normal[1] * Ud[2] + normal[2] * Ud[4]);
        g11++;
        g13++;
      }
    }
  }
}

/*
  Compute the linear tying strain bmat matrices. This corresponds to
  the derivative of the displacement-based strain at each of the tying
  points.

  input:
  knots: the 'order'th (ie 2, 3 or 4) Gauss-points
  pknots: the 'order-1'th Gauss-points
  Xpts: the nodal locations for the element

  output:
  b11: the derivative of the in-plane normal tensorial strain in the 1
  direction
  b22: the derivative of the in-plane normal tensorial strain in the 2
  direction
  b12: the derivative of the in-plane tensorial shear strain

  b23: the derivative of the out-of-plane tensorial shear strain in
  the 2-direction
  b13: the derivative of the out-of-plane tensorial shear strain in
  the 1-direction
*/
template <int order, int tying_order>
void compute_tying_bmat(const int is_linear, TacsScalar g11[], TacsScalar g22[],
                        TacsScalar g12[], TacsScalar g23[], TacsScalar g13[],
                        TacsScalar b11[], TacsScalar b22[], TacsScalar b12[],
                        TacsScalar b23[], TacsScalar b13[],
                        const double knots[], const double pknots[],
                        const TacsScalar vars[], const TacsScalar Xpts[]) {
  double N[order * order], Na[order * order], Nb[order * order];

  // Evaluate g22 and g23
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order; n++) {
      double pt[2];
      pt[0] = knots[n];
      pt[1] = pknots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      Tensor::normalize3D(normal);

      // Calculate the rotation
      TacsScalar r[3];
      Tensor::crossProduct3D(r, Urot, normal);

      if (is_linear) {
        g22[0] = Xd[3] * Ud[1] + Xd[4] * Ud[3] + Xd[5] * Ud[5];
        g23[0] =
            0.5 * (Xd[3] * r[0] + Xd[4] * r[1] + Xd[5] * r[2] +
                   normal[0] * Ud[1] + normal[1] * Ud[3] + normal[2] * Ud[5]);
        g22++;
        g23++;

        // Add the derivative to the b-matrix
        for (int i = 0; i < order * order; i++) {
          b22[0] = Xd[3] * Nb[i];
          b22[1] = Xd[4] * Nb[i];
          b22[2] = Xd[5] * Nb[i];
          b22 += 3;

          b23[0] = 0.5 * normal[0] * Nb[i];
          b23[1] = 0.5 * normal[1] * Nb[i];
          b23[2] = 0.5 * normal[2] * Nb[i];
          b23[3] = 0.5 * N[i] * (Xd[5] * normal[1] - Xd[4] * normal[2]);
          b23[4] = 0.5 * N[i] * (Xd[3] * normal[2] - Xd[5] * normal[0]);
          b23[5] = 0.5 * N[i] * (Xd[4] * normal[0] - Xd[3] * normal[1]);
          b23 += 6;
        }
      } else {
        g22[0] = (Xd[3] * Ud[1] + Xd[4] * Ud[3] + Xd[5] * Ud[5] +
                  0.5 * (Ud[1] * Ud[1] + Ud[3] * Ud[3] + Ud[5] * Ud[5]));
        g23[0] =
            0.5 * (Xd[3] * r[0] + Xd[4] * r[1] + Xd[5] * r[2] +
                   normal[0] * Ud[1] + normal[1] * Ud[3] + normal[2] * Ud[5]);
        g22++;
        g23++;

        // Add the derivative to the b-matrix
        for (int i = 0; i < order * order; i++) {
          b22[0] = (Xd[3] + Ud[1]) * Nb[i];
          b22[1] = (Xd[4] + Ud[3]) * Nb[i];
          b22[2] = (Xd[5] + Ud[5]) * Nb[i];
          b22 += 3;

          b23[0] = 0.5 * normal[0] * Nb[i];
          b23[1] = 0.5 * normal[1] * Nb[i];
          b23[2] = 0.5 * normal[2] * Nb[i];
          b23[3] = 0.5 * N[i] * (Xd[5] * normal[1] - Xd[4] * normal[2]);
          b23[4] = 0.5 * N[i] * (Xd[3] * normal[2] - Xd[5] * normal[0]);
          b23[5] = 0.5 * N[i] * (Xd[4] * normal[0] - Xd[3] * normal[1]);
          b23 += 6;
        }
      }
    }
  }

  // Evaluate g12
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = pknots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      if (is_linear) {
        g12[0] = 0.5 * (Xd[0] * Ud[1] + Xd[1] * Ud[3] + Xd[2] * Ud[5] +
                        Xd[3] * Ud[0] + Xd[4] * Ud[2] + Xd[5] * Ud[4]);
        g12++;

        for (int i = 0; i < order * order; i++) {
          b12[0] = 0.5 * (Xd[3] * Na[i] + Xd[0] * Nb[i]);
          b12[1] = 0.5 * (Xd[4] * Na[i] + Xd[1] * Nb[i]);
          b12[2] = 0.5 * (Xd[5] * Na[i] + Xd[2] * Nb[i]);
          b12 += 3;
        }
      } else {
        g12[0] = 0.5 * (Xd[0] * Ud[1] + Xd[1] * Ud[3] + Xd[2] * Ud[5] +
                        Xd[3] * Ud[0] + Xd[4] * Ud[2] + Xd[5] * Ud[4] +
                        Ud[0] * Ud[1] + Ud[2] * Ud[3] + Ud[4] * Ud[5]);
        g12++;

        for (int i = 0; i < order * order; i++) {
          b12[0] = 0.5 * ((Xd[3] + Ud[1]) * Na[i] + (Xd[0] + Ud[0]) * Nb[i]);
          b12[1] = 0.5 * ((Xd[4] + Ud[3]) * Na[i] + (Xd[1] + Ud[2]) * Nb[i]);
          b12[2] = 0.5 * ((Xd[5] + Ud[5]) * Na[i] + (Xd[2] + Ud[4]) * Nb[i]);
          b12 += 3;
        }
      }
    }
  }

  // Evaluate g11 and g13
  for (int m = 0; m < tying_order; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = knots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      Tensor::normalize3D(normal);

      // Calculate the rotation
      TacsScalar r[3];
      Tensor::crossProduct3D(r, Urot, normal);

      if (is_linear) {
        g11[0] = Xd[0] * Ud[0] + Xd[1] * Ud[2] + Xd[2] * Ud[4];
        g13[0] =
            0.5 * (Xd[0] * r[0] + Xd[1] * r[1] + Xd[2] * r[2] +
                   normal[0] * Ud[0] + normal[1] * Ud[2] + normal[2] * Ud[4]);
        g11++;
        g13++;

        for (int i = 0; i < order * order; i++) {
          b11[0] = Xd[0] * Na[i];
          b11[1] = Xd[1] * Na[i];
          b11[2] = Xd[2] * Na[i];
          b11 += 3;

          b13[0] = 0.5 * normal[0] * Na[i];
          b13[1] = 0.5 * normal[1] * Na[i];
          b13[2] = 0.5 * normal[2] * Na[i];
          b13[3] = 0.5 * N[i] * (Xd[2] * normal[1] - Xd[1] * normal[2]);
          b13[4] = 0.5 * N[i] * (Xd[0] * normal[2] - Xd[2] * normal[0]);
          b13[5] = 0.5 * N[i] * (Xd[1] * normal[0] - Xd[0] * normal[1]);
          b13 += 6;
        }
      } else {
        g11[0] = (Xd[0] * Ud[0] + Xd[1] * Ud[2] + Xd[2] * Ud[4] +
                  0.5 * (Ud[0] * Ud[0] + Ud[2] * Ud[2] + Ud[4] * Ud[4]));
        g13[0] =
            0.5 * (Xd[0] * r[0] + Xd[1] * r[1] + Xd[2] * r[2] +
                   normal[0] * Ud[0] + normal[1] * Ud[2] + normal[2] * Ud[4]);
        g11++;
        g13++;

        for (int i = 0; i < order * order; i++) {
          b11[0] = (Xd[0] + Ud[0]) * Na[i];
          b11[1] = (Xd[1] + Ud[2]) * Na[i];
          b11[2] = (Xd[2] + Ud[4]) * Na[i];
          b11 += 3;

          b13[0] = 0.5 * normal[0] * Na[i];
          b13[1] = 0.5 * normal[1] * Na[i];
          b13[2] = 0.5 * normal[2] * Na[i];
          b13[3] = 0.5 * N[i] * (Xd[2] * normal[1] - Xd[1] * normal[2]);
          b13[4] = 0.5 * N[i] * (Xd[0] * normal[2] - Xd[2] * normal[0]);
          b13[5] = 0.5 * N[i] * (Xd[1] * normal[0] - Xd[0] * normal[1]);
          b13 += 6;
        }
      }
    }
  }
}

/*
  Evaluate the strain at the tying points and the sensitivity of the
  strain at the tying points w.r.t. the nodal locations.

  input:
  knots: the 'order'th (ie 2, 3 or 4) Gauss-points
  pknots: the 'order-1'th Gauss-points
  Xpts: the nodal locations for the element

  output:
  g11, g22, g12, g23, g13: the in-plane tensorial strains and
  out-of-plane tensorial shear strains at the tying points

  dg11, dg22, dg12, dg23, dg13: the sensitivities of the the in-plane
  tensorial strains and out-of-plane tensorial shear strains at the
  tying points w.r.t. the nodal locataions
*/
template <int order, int tying_order>
void compute_tying_strain_sens(const int is_linear, TacsScalar g11[],
                               TacsScalar g22[], TacsScalar g12[],
                               TacsScalar g23[], TacsScalar g13[],
                               TacsScalar dg11[], TacsScalar dg22[],
                               TacsScalar dg12[], TacsScalar dg23[],
                               TacsScalar dg13[], const double knots[],
                               const double pknots[], const TacsScalar vars[],
                               const TacsScalar Xpts[]) {
  double N[order * order], Na[order * order], Nb[order * order];

  // Evaluate g22 and g23
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order; n++) {
      double pt[2];
      pt[0] = knots[n];
      pt[1] = pknots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      TacsScalar inv_norm = Tensor::normalize3D(normal);
      inv_norm = 1.0 / inv_norm;

      // Calculate the rotation
      TacsScalar r[3];
      Tensor::crossProduct3D(r, Urot, normal);

      if (is_linear) {
        g22[0] = Xd[3] * Ud[1] + Xd[4] * Ud[3] + Xd[5] * Ud[5];
        g23[0] =
            0.5 * (Xd[3] * r[0] + Xd[4] * r[1] + Xd[5] * r[2] +
                   normal[0] * Ud[1] + normal[1] * Ud[3] + normal[2] * Ud[5]);
        g22++;
        g23++;
      } else {
        g22[0] = (Xd[3] * Ud[1] + Xd[4] * Ud[3] + Xd[5] * Ud[5] +
                  0.5 * (Ud[1] * Ud[1] + Ud[3] * Ud[3] + Ud[5] * Ud[5]));
        g23[0] =
            0.5 * (Xd[3] * r[0] + Xd[4] * r[1] + Xd[5] * r[2] +
                   normal[0] * Ud[1] + normal[1] * Ud[3] + normal[2] * Ud[5]);
        g22++;
        g23++;
      }

      for (int i = 0; i < order * order; i++) {
        for (int k = 0; k < 3; k++) {
          TacsScalar normalSens[3];

          if (k == 0) {
            normalSens[0] = 0.0;
            normalSens[1] = Xd[2] * Nb[i] - Na[i] * Xd[5];
            normalSens[2] = Na[i] * Xd[4] - Xd[1] * Nb[i];
          } else if (k == 1) {
            normalSens[0] = Na[i] * Xd[5] - Xd[2] * Nb[i];
            normalSens[1] = 0.0;
            normalSens[2] = Xd[0] * Nb[i] - Na[i] * Xd[3];
          } else {
            normalSens[0] = Xd[1] * Nb[i] - Na[i] * Xd[4];
            normalSens[1] = Na[i] * Xd[3] - Xd[0] * Nb[i];
            normalSens[2] = 0.0;
          }

          TacsScalar normSens =
              (normal[0] * normalSens[0] + normal[1] * normalSens[1] +
               normal[2] * normalSens[2]);

          normalSens[0] = (normalSens[0] - normSens * normal[0]) * inv_norm;
          normalSens[1] = (normalSens[1] - normSens * normal[1]) * inv_norm;
          normalSens[2] = (normalSens[2] - normSens * normal[2]) * inv_norm;

          // Calculate the rotation
          TacsScalar rSens[3];
          Tensor::crossProduct3D(rSens, Urot, normalSens);

          dg22[0] = Ud[2 * k + 1] * Nb[i];
          dg23[0] = 0.5 * (r[k] * Nb[i] + Xd[3] * rSens[0] + Xd[4] * rSens[1] +
                           Xd[5] * rSens[2] + normalSens[0] * Ud[1] +
                           normalSens[1] * Ud[3] + normalSens[2] * Ud[5]);
          dg22++;
          dg23++;
        }
      }
    }
  }

  // Evaluate g12
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = pknots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      if (is_linear) {
        g12[0] = 0.5 * (Xd[0] * Ud[1] + Xd[1] * Ud[3] + Xd[2] * Ud[5] +
                        Xd[3] * Ud[0] + Xd[4] * Ud[2] + Xd[5] * Ud[4]);
        g12++;
      } else {
        g12[0] = 0.5 * (Xd[0] * Ud[1] + Xd[1] * Ud[3] + Xd[2] * Ud[5] +
                        Xd[3] * Ud[0] + Xd[4] * Ud[2] + Xd[5] * Ud[4] +
                        Ud[0] * Ud[1] + Ud[2] * Ud[3] + Ud[4] * Ud[5]);
        g12++;
      }

      for (int i = 0; i < order * order; i++) {
        for (int k = 0; k < 3; k++) {
          dg12[0] = 0.5 * (Nb[i] * Ud[2 * k] + Na[i] * Ud[2 * k + 1]);
          dg12++;
        }
      }
    }
  }

  // Evaluate g11 and g13
  for (int m = 0; m < tying_order; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = knots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      TacsScalar inv_norm = Tensor::normalize3D(normal);
      inv_norm = 1.0 / inv_norm;

      // Calculate the rotation
      TacsScalar r[3];
      Tensor::crossProduct3D(r, Urot, normal);

      if (is_linear) {
        g11[0] = Xd[0] * Ud[0] + Xd[1] * Ud[2] + Xd[2] * Ud[4];
        g13[0] =
            0.5 * (Xd[0] * r[0] + Xd[1] * r[1] + Xd[2] * r[2] +
                   normal[0] * Ud[0] + normal[1] * Ud[2] + normal[2] * Ud[4]);
        g11++;
        g13++;
      } else {
        g11[0] = (Xd[0] * Ud[0] + Xd[1] * Ud[2] + Xd[2] * Ud[4] +
                  0.5 * (Ud[0] * Ud[0] + Ud[2] * Ud[2] + Ud[4] * Ud[4]));
        g13[0] =
            0.5 * (Xd[0] * r[0] + Xd[1] * r[1] + Xd[2] * r[2] +
                   normal[0] * Ud[0] + normal[1] * Ud[2] + normal[2] * Ud[4]);
        g11++;
        g13++;
      }

      for (int i = 0; i < order * order; i++) {
        for (int k = 0; k < 3; k++) {
          TacsScalar normalSens[3];

          if (k == 0) {
            normalSens[0] = 0.0;
            normalSens[1] = Xd[2] * Nb[i] - Na[i] * Xd[5];
            normalSens[2] = Na[i] * Xd[4] - Xd[1] * Nb[i];
          } else if (k == 1) {
            normalSens[0] = Na[i] * Xd[5] - Xd[2] * Nb[i];
            normalSens[1] = 0.0;
            normalSens[2] = Xd[0] * Nb[i] - Na[i] * Xd[3];
          } else {
            normalSens[0] = Xd[1] * Nb[i] - Na[i] * Xd[4];
            normalSens[1] = Na[i] * Xd[3] - Xd[0] * Nb[i];
            normalSens[2] = 0.0;
          }

          TacsScalar normSens =
              (normal[0] * normalSens[0] + normal[1] * normalSens[1] +
               normal[2] * normalSens[2]);

          normalSens[0] = (normalSens[0] - normSens * normal[0]) * inv_norm;
          normalSens[1] = (normalSens[1] - normSens * normal[1]) * inv_norm;
          normalSens[2] = (normalSens[2] - normSens * normal[2]) * inv_norm;

          // Calculate the rotation
          TacsScalar rSens[3];
          Tensor::crossProduct3D(rSens, Urot, normalSens);

          dg11[0] = Ud[2 * k] * Na[i];
          dg13[0] = 0.5 * (Na[i] * r[k] + Xd[0] * rSens[0] + Xd[1] * rSens[1] +
                           Xd[2] * rSens[2] + normalSens[0] * Ud[0] +
                           normalSens[1] * Ud[2] + normalSens[2] * Ud[4]);
          dg11++;
          dg13++;
        }
      }
    }
  }
}

/*
  Add the contribution from the tying strain components to the
  derivative of the residual w.r.t. the nodal coordinates.  This
  function, unlike the others in this file, adds the tying strain
  contributions only at a point - at a specific value of N11, N22 and
  N12. This saves on a significant amount of memory, without
  sacrificing too much in terms of performance.

  inputs:
  scale: scale all contributions by this value
  stress: the values of the stresses at the current point
  tx: the transformation between global and local ref. frames
  dtx: the derivative of tx w.r.t. each coordinate

  knots, pknots: the vectors of Gauss points for the tying strain
  locations

  Xpts: the nodal coordinates
  N11, N22, N12: the shape functions for the tying strain
  interpolation
*/
template <int order, int tying_order>
void add_tying_bmat_sens(const int is_linear, TacsScalar fXptSens[],
                         const TacsScalar psi[], const TacsScalar scale,
                         const TacsScalar stress[], const TacsScalar tx[],
                         const TacsScalar _dtx[], const double knots[],
                         const double pknots[], const TacsScalar vars[],
                         const TacsScalar Xpts[], const double N11[],
                         const double N22[], const double N12[]) {
  TacsScalar B22[3 * order * order], B23[6 * order * order];
  TacsScalar B12[3 * order * order];
  TacsScalar B11[3 * order * order], B13[6 * order * order];

  TacsScalar DB22[9 * order * order * order * order];
  TacsScalar DB23[18 * order * order * order * order];
  TacsScalar DB12[9 * order * order * order * order];
  TacsScalar DB11[9 * order * order * order * order];
  TacsScalar DB13[18 * order * order * order * order];

  // Zero the b-matrices
  memset(B22, 0, 3 * order * order * sizeof(TacsScalar));
  memset(B23, 0, 6 * order * order * sizeof(TacsScalar));
  memset(B12, 0, 3 * order * order * sizeof(TacsScalar));
  memset(B11, 0, 3 * order * order * sizeof(TacsScalar));
  memset(B13, 0, 6 * order * order * sizeof(TacsScalar));

  // Zero the derivatives of the b-matrices
  memset(DB22, 0, 9 * order * order * order * order * sizeof(TacsScalar));
  memset(DB23, 0, 18 * order * order * order * order * sizeof(TacsScalar));
  memset(DB12, 0, 9 * order * order * order * order * sizeof(TacsScalar));
  memset(DB11, 0, 9 * order * order * order * order * sizeof(TacsScalar));
  memset(DB13, 0, 18 * order * order * order * order * sizeof(TacsScalar));

  double N[order * order], Na[order * order], Nb[order * order];

  // Evaluate g22 and g23
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order; n++) {
      double pt[2];
      pt[0] = knots[n];
      pt[1] = pknots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      TacsScalar inv_norm = Tensor::normalize3D(normal);
      inv_norm = 1.0 / inv_norm;

      TacsScalar *b22 = B22;
      TacsScalar *b23 = B23;
      TacsScalar *db22 = DB22;
      TacsScalar *db23 = DB23;

      for (int i = 0; i < order * order; i++) {
        double Ni = N22[0] * N[i];
        double Nbi = N22[0] * Nb[i];

        // Sum up the contributions to the B-matrices
        if (is_linear) {
          b22[0] += Xd[3] * Nbi;
          b22[1] += Xd[4] * Nbi;
          b22[2] += Xd[5] * Nbi;
        } else {
          b22[0] += (Xd[3] + Ud[1]) * Nbi;
          b22[1] += (Xd[4] + Ud[3]) * Nbi;
          b22[2] += (Xd[5] + Ud[5]) * Nbi;
        }
        b22 += 3;

        b23[0] += 0.5 * normal[0] * Nbi;
        b23[1] += 0.5 * normal[1] * Nbi;
        b23[2] += 0.5 * normal[2] * Nbi;
        b23[3] += 0.5 * Ni * (Xd[5] * normal[1] - Xd[4] * normal[2]);
        b23[4] += 0.5 * Ni * (Xd[3] * normal[2] - Xd[5] * normal[0]);
        b23[5] += 0.5 * Ni * (Xd[4] * normal[0] - Xd[3] * normal[1]);
        b23 += 6;

        // Add the contributions to the derivative of the B-matrices
        for (int j = 0; j < order * order; j++) {
          double Naj = Na[j];
          double Nbj = Nb[j];

          for (int k = 0; k < 3; k++) {
            TacsScalar normalSens[3];

            if (k == 0) {
              normalSens[0] = 0.0;
              normalSens[1] = Xd[2] * Nbj - Naj * Xd[5];
              normalSens[2] = Naj * Xd[4] - Xd[1] * Nbj;
            } else if (k == 1) {
              normalSens[0] = Naj * Xd[5] - Xd[2] * Nbj;
              normalSens[1] = 0.0;
              normalSens[2] = Xd[0] * Nbj - Naj * Xd[3];
            } else {
              normalSens[0] = Xd[1] * Nbj - Naj * Xd[4];
              normalSens[1] = Naj * Xd[3] - Xd[0] * Nbj;
              normalSens[2] = 0.0;
            }

            TacsScalar normSens =
                (normal[0] * normalSens[0] + normal[1] * normalSens[1] +
                 normal[2] * normalSens[2]);
            normalSens[0] = (normalSens[0] - normSens * normal[0]) * inv_norm;
            normalSens[1] = (normalSens[1] - normSens * normal[1]) * inv_norm;
            normalSens[2] = (normalSens[2] - normSens * normal[2]) * inv_norm;

            db22[k] += Nbj * Nbi;
            db22 += 3;

            db23[0] += 0.5 * normalSens[0] * Nbi;
            db23[1] += 0.5 * normalSens[1] * Nbi;
            db23[2] += 0.5 * normalSens[2] * Nbi;
            db23[3] +=
                0.5 * Ni * (Xd[5] * normalSens[1] - Xd[4] * normalSens[2]);
            db23[4] +=
                0.5 * Ni * (Xd[3] * normalSens[2] - Xd[5] * normalSens[0]);
            db23[5] +=
                0.5 * Ni * (Xd[4] * normalSens[0] - Xd[3] * normalSens[1]);

            if (k == 0) {
              db23[4] += 0.5 * Ni * Nbj * normal[2];
              db23[5] -= 0.5 * Ni * Nbj * normal[1];
            } else if (k == 1) {
              db23[3] -= 0.5 * Ni * Nbj * normal[2];
              db23[5] += 0.5 * Ni * Nbj * normal[0];
            } else {  // k == 2
              db23[3] += 0.5 * Ni * Nbj * normal[1];
              db23[4] -= 0.5 * Ni * Nbj * normal[0];
            }

            db23 += 6;
          }
        }
      }

      N22++;
    }
  }

  // Evaluate g12
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = pknots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      TacsScalar *b12 = B12;
      TacsScalar *db12 = DB12;
      for (int i = 0; i < order * order; i++) {
        double Nai = N12[0] * Na[i];
        double Nbi = N12[0] * Nb[i];

        if (is_linear) {
          b12[0] += 0.5 * (Xd[3] * Nai + Xd[0] * Nbi);
          b12[1] += 0.5 * (Xd[4] * Nai + Xd[1] * Nbi);
          b12[2] += 0.5 * (Xd[5] * Nai + Xd[2] * Nbi);
        } else {
          b12[0] += 0.5 * ((Xd[3] + Ud[1]) * Nai + (Xd[0] + Ud[0]) * Nbi);
          b12[1] += 0.5 * ((Xd[4] + Ud[3]) * Nai + (Xd[1] + Ud[2]) * Nbi);
          b12[2] += 0.5 * ((Xd[5] + Ud[5]) * Nai + (Xd[2] + Ud[4]) * Nbi);
        }
        b12 += 3;

        for (int j = 0; j < order * order; j++) {
          double Naj = Na[j];
          double Nbj = Nb[j];
          db12[0] += 0.5 * (Nbj * Nai + Naj * Nbi);
          db12 += 3;
          db12[1] += 0.5 * (Nbj * Nai + Naj * Nbi);
          db12 += 3;
          db12[2] += 0.5 * (Nbj * Nai + Naj * Nbi);
          db12 += 3;
        }
      }
      N12++;
    }
  }

  // Evaluate g11 and g13
  for (int m = 0; m < tying_order; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = knots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      TacsScalar inv_norm = Tensor::normalize3D(normal);
      inv_norm = 1.0 / inv_norm;

      // Calculate the rotation
      TacsScalar r[3];
      Tensor::crossProduct3D(r, Urot, normal);

      TacsScalar *b11 = B11;
      TacsScalar *b13 = B13;
      TacsScalar *db11 = DB11;
      TacsScalar *db13 = DB13;

      for (int i = 0; i < order * order; i++) {
        double Ni = N11[0] * N[i];
        double Nai = N11[0] * Na[i];

        if (is_linear) {
          b11[0] += Xd[0] * Nai;
          b11[1] += Xd[1] * Nai;
          b11[2] += Xd[2] * Nai;
        } else {
          b11[0] += (Xd[0] + Ud[0]) * Nai;
          b11[1] += (Xd[1] + Ud[2]) * Nai;
          b11[2] += (Xd[2] + Ud[4]) * Nai;
        }
        b11 += 3;

        b13[0] += 0.5 * normal[0] * Nai;
        b13[1] += 0.5 * normal[1] * Nai;
        b13[2] += 0.5 * normal[2] * Nai;
        b13[3] += 0.5 * Ni * (Xd[2] * normal[1] - Xd[1] * normal[2]);
        b13[4] += 0.5 * Ni * (Xd[0] * normal[2] - Xd[2] * normal[0]);
        b13[5] += 0.5 * Ni * (Xd[1] * normal[0] - Xd[0] * normal[1]);
        b13 += 6;

        for (int j = 0; j < order * order; j++) {
          double Naj = Na[j];
          double Nbj = Nb[j];

          for (int k = 0; k < 3; k++) {
            TacsScalar normalSens[3];

            if (k == 0) {
              normalSens[0] = 0.0;
              normalSens[1] = Xd[2] * Nbj - Naj * Xd[5];
              normalSens[2] = Naj * Xd[4] - Xd[1] * Nbj;
            } else if (k == 1) {
              normalSens[0] = Naj * Xd[5] - Xd[2] * Nbj;
              normalSens[1] = 0.0;
              normalSens[2] = Xd[0] * Nbj - Naj * Xd[3];
            } else {
              normalSens[0] = Xd[1] * Nbj - Naj * Xd[4];
              normalSens[1] = Naj * Xd[3] - Xd[0] * Nbj;
              normalSens[2] = 0.0;
            }

            TacsScalar normSens =
                (normal[0] * normalSens[0] + normal[1] * normalSens[1] +
                 normal[2] * normalSens[2]);

            normalSens[0] = (normalSens[0] - normSens * normal[0]) * inv_norm;
            normalSens[1] = (normalSens[1] - normSens * normal[1]) * inv_norm;
            normalSens[2] = (normalSens[2] - normSens * normal[2]) * inv_norm;

            db11[k] += Naj * Nai;
            db11 += 3;

            db13[0] += 0.5 * normalSens[0] * Nai;
            db13[1] += 0.5 * normalSens[1] * Nai;
            db13[2] += 0.5 * normalSens[2] * Nai;
            db13[3] +=
                0.5 * Ni * (Xd[2] * normalSens[1] - Xd[1] * normalSens[2]);
            db13[4] +=
                0.5 * Ni * (Xd[0] * normalSens[2] - Xd[2] * normalSens[0]);
            db13[5] +=
                0.5 * Ni * (Xd[1] * normalSens[0] - Xd[0] * normalSens[1]);

            if (k == 0) {
              db13[4] += 0.5 * Ni * Naj * normal[2];
              db13[5] -= 0.5 * Ni * Naj * normal[1];
            } else if (k == 1) {
              db13[3] -= 0.5 * Ni * Naj * normal[2];
              db13[5] += 0.5 * Ni * Naj * normal[0];
            } else {  // k == 2
              db13[3] += 0.5 * Ni * Naj * normal[1];
              db13[4] -= 0.5 * Ni * Naj * normal[0];
            }

            db13 += 6;
          }
        }
      }

      N11++;
    }
  }

  // B22, B23, B12, B11 and B13 contain the total contribution to
  // the B-matrix interpolated from the tying points using N11, N12, N22

  // DB22, DB23, DB12, DB11 and DB13 contain the derivatives of the
  // B-matrix interpolated from the tying points using N11, N12, N22

  // Next, transform these contributions and add their values to
  // the derivative of the residuals
  const int nc = 3 * order * order;

  for (int i = 0; i < order * order; i++) {
    for (int ii = 0; ii < 6; ii++) {
      const int row = 6 * i + ii;

      if (ii < 3) {
        const int row3 = 3 * i + ii;

        // Retrieve the values of the B-matrix
        TacsScalar g[6];
        g[0] = B11[row3];
        g[1] = B22[row3];
        g[2] = 0.0;
        g[3] = B23[row];
        g[4] = B13[row];
        g[5] = B12[row3];

        // Complete the first part of the transformation
        // Compute as = tx * G - where G is a matrix formed
        // as follows:
        // G = [ g[0]  g[5]  g[4] ]
        //     [ g[5]  g[1]  g[3] ]
        //     [ g[4]  g[3]  g[2] ]
        // The complete transformation is tx * G * tx^{T}, but
        // only the derivative is required in this computation

        TacsScalar as[9];
        as[0] = tx[0] * g[0] + tx[1] * g[5] + tx[2] * g[4];
        as[3] = tx[3] * g[0] + tx[4] * g[5] + tx[5] * g[4];
        as[6] = tx[6] * g[0] + tx[7] * g[5] + tx[8] * g[4];

        as[1] = tx[0] * g[5] + tx[1] * g[1] + tx[2] * g[3];
        as[4] = tx[3] * g[5] + tx[4] * g[1] + tx[5] * g[3];
        as[7] = tx[6] * g[5] + tx[7] * g[1] + tx[8] * g[3];

        as[2] = tx[0] * g[4] + tx[1] * g[3] + tx[2] * g[2];
        as[5] = tx[3] * g[4] + tx[4] * g[3] + tx[5] * g[2];
        as[8] = tx[6] * g[4] + tx[7] * g[3] + tx[8] * g[2];

        const TacsScalar *dtx = _dtx;
        for (int k = 0; k < nc; k++) {
          // Retrieve the derivative of the B-matrix
          TacsScalar dg[6];
          // 3*i + ii
          // ii + 3*k + 3*nc*i

          dg[0] = DB11[3 * nc * i + 3 * k + ii];  // DB11[nc*row3 + k];
          dg[1] = DB22[3 * nc * i + 3 * k + ii];
          dg[2] = 0.0;
          dg[3] = DB23[6 * nc * i + 6 * k + ii];
          dg[4] = DB13[6 * nc * i + 6 * k + ii];
          dg[5] = DB12[3 * nc * i + 3 * k + ii];

          // Get the derivative of the temporary variable as
          TacsScalar das[9];
          das[0] = (dtx[0] * g[0] + dtx[1] * g[5] + dtx[2] * g[4] +
                    tx[0] * dg[0] + tx[1] * dg[5] + tx[2] * dg[4]);
          das[3] = (dtx[3] * g[0] + dtx[4] * g[5] + dtx[5] * g[4] +
                    tx[3] * dg[0] + tx[4] * dg[5] + tx[5] * dg[4]);
          das[6] = (dtx[6] * g[0] + dtx[7] * g[5] + dtx[8] * g[4] +
                    tx[6] * dg[0] + tx[7] * dg[5] + tx[8] * dg[4]);

          das[1] = (dtx[0] * g[5] + dtx[1] * g[1] + dtx[2] * g[3] +
                    tx[0] * dg[5] + tx[1] * dg[1] + tx[2] * dg[3]);
          das[4] = (dtx[3] * g[5] + dtx[4] * g[1] + dtx[5] * g[3] +
                    tx[3] * dg[5] + tx[4] * dg[1] + tx[5] * dg[3]);
          das[7] = (dtx[6] * g[5] + dtx[7] * g[1] + dtx[8] * g[3] +
                    tx[6] * dg[5] + tx[7] * dg[1] + tx[8] * dg[3]);

          das[2] = (dtx[0] * g[4] + dtx[1] * g[3] + dtx[2] * g[2] +
                    tx[0] * dg[4] + tx[1] * dg[3] + tx[2] * dg[2]);
          das[5] = (dtx[3] * g[4] + dtx[4] * g[3] + dtx[5] * g[2] +
                    tx[3] * dg[4] + tx[4] * dg[3] + tx[5] * dg[2]);
          das[8] = (dtx[6] * g[4] + dtx[7] * g[3] + dtx[8] * g[2] +
                    tx[6] * dg[4] + tx[7] * dg[3] + tx[8] * dg[2]);

          fXptSens[k] +=
              scale * psi[row] *
              (stress[0] * (as[0] * dtx[0] + as[1] * dtx[1] + as[2] * dtx[2] +
                            das[0] * tx[0] + das[1] * tx[1] + das[2] * tx[2]) +
               stress[1] * (as[3] * dtx[3] + as[4] * dtx[4] + as[5] * dtx[5] +
                            das[3] * tx[3] + das[4] * tx[4] + das[5] * tx[5]) +
               2.0 * stress[2] *
                   (as[0] * dtx[3] + as[1] * dtx[4] + as[2] * dtx[5] +
                    das[0] * tx[3] + das[1] * tx[4] + das[2] * tx[5]) +
               2.0 * stress[6] *
                   (as[3] * dtx[6] + as[4] * dtx[7] + as[5] * dtx[8] +
                    das[3] * tx[6] + das[4] * tx[7] + das[5] * tx[8]) +
               2.0 * stress[7] *
                   (as[0] * dtx[6] + as[1] * dtx[7] + as[2] * dtx[8] +
                    das[0] * tx[6] + das[1] * tx[7] + das[2] * tx[8]));

          dtx += 9;
        }
      } else {
        // Retrieve the values of the B-matrix
        TacsScalar g[6];
        g[0] = g[1] = g[2] = 0.0;
        g[3] = B23[row];
        g[4] = B13[row];
        g[5] = 0.0;

        // Complete the first part of the transformation
        // Compute as = tx * G. This computation ignores zero entries
        TacsScalar as[9];
        as[0] = tx[2] * g[4];
        as[3] = tx[5] * g[4];
        as[6] = tx[8] * g[4];

        as[1] = tx[2] * g[3];
        as[4] = tx[5] * g[3];
        as[7] = tx[8] * g[3];

        as[2] = tx[0] * g[4] + tx[1] * g[3];
        as[5] = tx[3] * g[4] + tx[4] * g[3];
        as[8] = tx[6] * g[4] + tx[7] * g[3];

        // Point dtx at the beginning of the derivatives of tx
        const TacsScalar *dtx = _dtx;
        for (int k = 0; k < nc; k++) {
          // Retrieve the derivative of the B-matrix
          TacsScalar dg[6];
          dg[0] = dg[1] = dg[2] = 0.0;
          dg[3] = DB23[6 * nc * i + 6 * k + ii];
          dg[4] = DB13[6 * nc * i + 6 * k + ii];
          dg[5] = 0.0;

          // Get the derivative of the temporary variable as
          TacsScalar das[9];
          das[0] = (dtx[2] * g[4] + tx[2] * dg[4]);
          das[3] = (dtx[5] * g[4] + tx[5] * dg[4]);
          das[6] = (dtx[8] * g[4] + tx[8] * dg[4]);

          das[1] = (dtx[2] * g[3] + tx[2] * dg[3]);
          das[4] = (dtx[5] * g[3] + tx[5] * dg[3]);
          das[7] = (dtx[8] * g[3] + tx[8] * dg[3]);

          das[2] =
              (dtx[0] * g[4] + dtx[1] * g[3] + tx[0] * dg[4] + tx[1] * dg[3]);
          das[5] =
              (dtx[3] * g[4] + dtx[4] * g[3] + tx[3] * dg[4] + tx[4] * dg[3]);
          das[8] =
              (dtx[6] * g[4] + dtx[7] * g[3] + tx[6] * dg[4] + tx[7] * dg[3]);

          // Add the deterivative of tx * G * tx^{T} * stress to res
          fXptSens[k] +=
              scale * psi[row] *
              (stress[0] * (as[0] * dtx[0] + as[1] * dtx[1] + as[2] * dtx[2] +
                            das[0] * tx[0] + das[1] * tx[1] + das[2] * tx[2]) +
               stress[1] * (as[3] * dtx[3] + as[4] * dtx[4] + as[5] * dtx[5] +
                            das[3] * tx[3] + das[4] * tx[4] + das[5] * tx[5]) +
               2.0 * stress[2] *
                   (as[0] * dtx[3] + as[1] * dtx[4] + as[2] * dtx[5] +
                    das[0] * tx[3] + das[1] * tx[4] + das[2] * tx[5]) +
               2.0 * stress[6] *
                   (as[3] * dtx[6] + as[4] * dtx[7] + as[5] * dtx[8] +
                    das[3] * tx[6] + das[4] * tx[7] + das[5] * tx[8]) +
               2.0 * stress[7] *
                   (as[0] * dtx[6] + as[1] * dtx[7] + as[2] * dtx[8] +
                    das[0] * tx[6] + das[1] * tx[7] + das[2] * tx[8]));

          dtx += 9;
        }
      }
    }
  }
}

/*
  Compute the second derivative of the strain w.r.t. the
  displacements/rotations - many of these components are the same.

  input:
  N11, N22, N12: the shape functions associated with the interpolation
  of the strain components

  knots, pknots: the tying interpolation points

  ouput:
  n11, n22, n12: the second derivative of the strain w.r.t. the nodal
  displacements
*/
template <int order, int tying_order>
void nonlinear_tying_nmat(TacsScalar n11[], TacsScalar n22[], TacsScalar n12[],
                          const double N11[], const double N22[],
                          const double N12[], const double knots[],
                          const double pknots[]) {
  double N[order * order], Na[order * order], Nb[order * order];

  // Evaluate g22 and g23
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order; n++) {
      double pt[2];
      pt[0] = knots[n];
      pt[1] = pknots[m];

      FElibrary::biLagrangeSF(N, Na, Nb, pt, order);

      TacsScalar *_n22 = n22;
      for (int i = 0; i < order * order; i++) {
        for (int j = 0; j <= i; j++) {
          _n22[0] += Nb[i] * Nb[j] * N22[0];
          _n22++;
        }
      }
      N22++;
    }
  }

  // Evaluate g12
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = pknots[m];

      FElibrary::biLagrangeSF(N, Na, Nb, pt, order);

      TacsScalar *_n12 = n12;
      for (int i = 0; i < order * order; i++) {
        for (int j = 0; j <= i; j++) {
          _n12[0] += 0.5 * (Na[i] * Nb[j] + Nb[i] * Na[j]) * N12[0];
          _n12++;
        }
      }
      N12++;
    }
  }

  // Evaluate g11 and g13
  for (int m = 0; m < tying_order; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = knots[m];

      FElibrary::biLagrangeSF(N, Na, Nb, pt, order);

      TacsScalar *_n11 = n11;
      for (int i = 0; i < order * order; i++) {
        for (int j = 0; j <= i; j++) {
          _n11[0] += Na[i] * Na[j] * N11[0];
          _n11++;
        }
      }
      N11++;
    }
  }
}

/*
  Add the derivative of the product of stress and the derivative of
  the strain to part of the matrix.  This is used for the nonlinear
  contributions from the in-plane components of the tensorial strain.
*/
template <int order, int tying_order>
void add_nonlinear_tying_stress_nmat(TacsScalar matrix[], TacsScalar scale,
                                     const TacsScalar stress[],
                                     const TacsScalar tx[], const double N11[],
                                     const double N22[], const double N12[],
                                     const double knots[],
                                     const double pknots[],
                                     const TacsScalar Xpts[]) {
  const int nvars = 6 * order * order;
  const int size = (order * order * (order * order + 1)) / 2;

  TacsScalar n11[size], n22[size], n12[size];
  memset(n11, 0, size * sizeof(TacsScalar));
  memset(n22, 0, size * sizeof(TacsScalar));
  memset(n12, 0, size * sizeof(TacsScalar));

  nonlinear_tying_nmat<order, tying_order>(n11, n22, n12, N11, N22, N12, knots,
                                           pknots);

  TacsScalar *_n11 = n11;
  TacsScalar *_n12 = n12;
  TacsScalar *_n22 = n22;

  // Compute the second derivative contributions
  for (int i = 0; i < order * order; i++) {
    const int row = 6 * i;
    for (int j = 0; j <= i; j++) {
      const int col = 6 * j;

      // Evaluate the second derivative of the tensorial strain.
      // Note that this is NOT the engineering strain!
      TacsScalar g[6], s[6];
      g[0] = _n11[0];
      g[1] = _n22[0];
      g[5] = _n12[0];
      g[2] = g[3] = g[4] = 0.0;

      // Transform to the local axis - note that this employs the
      // transformation for stress - since we're using tensorial strain -
      // not engineering strain
      Tensor::transform3DStress(s, g, tx);

      // Evaluate the contribution - note the factor of 2.0 due to the
      // use of tensorial strain
      TacsScalar val =
          scale *
          (stress[0] * s[0] + stress[1] * s[1] +
           2.0 * (stress[2] * s[5] + stress[6] * s[3] + stress[7] * s[4]));

      // Add values to the diagonal of each component
      for (int k = 0; k < 3; k++) {
        matrix[(row + k) * nvars + col + k] += val;
      }

      // Increment the pointers
      _n11++;
      _n22++;
      _n12++;
    }
  }
}

/*
  Add the derivative of the product of stress and the derivative of
  the strain to part of the matrix.  This is used for the nonlinear
  contributions from the in-plane components of the tensorial strain.
*/
template <int order, int tying_order>
void add_nonlinear_tying_inner_nmat(TacsScalar bstrain[],
                                    const TacsScalar psi[],
                                    const TacsScalar phi[],
                                    const TacsScalar tx[], const double N11[],
                                    const double N22[], const double N12[],
                                    const double knots[], const double pknots[],
                                    const TacsScalar Xpts[]) {
  // Evaluate the second derivative of the tensorial strain.
  // Note that this is NOT the engineering strain!
  TacsScalar g[6], e[6];

  // Temporary space for the shape functions
  double N[order * order], Na[order * order], Nb[order * order];

  // Zero the components of the tensorial strain
  memset(g, 0, 6 * sizeof(TacsScalar));

  // Evaluate g11
  for (int m = 0; m < tying_order; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = knots[m];

      FElibrary::biLagrangeSF(N, Na, Nb, pt, order);

      for (int k = 0; k < 3; k++) {
        TacsScalar phia = 0.0, psia = 0.0;
        for (int i = 0; i < order * order; i++) {
          phia += phi[6 * i + k] * Na[i];
          psia += psi[6 * i + k] * Na[i];
        }

        g[0] += phia * psia * N11[0];
      }
      N11++;
    }
  }

  // Evaluate g22
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order; n++) {
      double pt[2];
      pt[0] = knots[n];
      pt[1] = pknots[m];

      FElibrary::biLagrangeSF(N, Na, Nb, pt, order);

      for (int k = 0; k < 3; k++) {
        TacsScalar phib = 0.0, psib = 0.0;
        for (int i = 0; i < order * order; i++) {
          phib += phi[6 * i + k] * Nb[i];
          psib += psi[6 * i + k] * Nb[i];
        }

        g[1] += phib * psib * N22[0];
      }
      N22++;
    }
  }

  // Evaluate g12
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = pknots[m];

      FElibrary::biLagrangeSF(N, Na, Nb, pt, order);

      for (int k = 0; k < 3; k++) {
        TacsScalar phia = 0.0, phib = 0.0;
        TacsScalar psia = 0.0, psib = 0.0;
        for (int i = 0; i < order * order; i++) {
          phia += Na[i] * phi[6 * i + k];
          phib += Nb[i] * phi[6 * i + k];
          psia += Na[i] * psi[6 * i + k];
          psib += Nb[i] * psi[6 * i + k];
        }

        g[5] += 0.5 * (phia * psib + psia * phib) * N12[0];
      }
      N12++;
    }
  }

  // Transform to the local axis - note that this employs the
  // transformation for stress - since we're using tensorial strain -
  // not engineering strain
  Tensor::transform3DStress(e, g, tx);

  bstrain[0] = e[0];
  bstrain[1] = e[1];
  bstrain[2] = 2.0 * e[5];
  bstrain[6] = 2.0 * e[3];
  bstrain[7] = 2.0 * e[4];
}

/*
  Add the tying strain contribution to the strains. This adds the
  result of the interpolation of the displacement-based strain from
  the tying points to the strain. Only the in-plane and shear strain
  components are added.

  input:
  tx: the transformation to the local coordinates
  g11, g22, g12, g23, g13: the components of the tying strain

  N11, N22, N12: the interpolations

  output:
  strain: the strain values are set into the in-plane components (0, 1, 2)
  and the shear strain components (6, 7)
*/
template <int tying_order>
void add_tying_strain(TacsScalar strain[], const TacsScalar tx[],
                      const TacsScalar g11[], const TacsScalar g22[],
                      const TacsScalar g12[], const TacsScalar g23[],
                      const TacsScalar g13[], const double N11[],
                      const double N22[], const double N12[]) {
  TacsScalar g[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  // Use the interpolations
  for (int k = 0; k < tying_order * (tying_order - 1); k++) {
    g[0] += g11[k] * N11[k];
    g[1] += g22[k] * N22[k];

    g[3] += g23[k] * N22[k];
    g[4] += g13[k] * N11[k];
  }

  for (int k = 0; k < (tying_order - 1) * (tying_order - 1); k++) {
    g[5] += g12[k] * N12[k];
  }

  // Transform the tensorial strain g, to the local axis
  TacsScalar s[6];
  Tensor::transform3DStress(s, g, tx);
  strain[0] = s[0];
  strain[1] = s[1];
  strain[2] = 2.0 * s[5];
  strain[6] = 2.0 * s[3];
  strain[7] = 2.0 * s[4];
}

/*
  Add the B matrices (the derivative of the strain w.r.t. the nodal
  displacements and rotations at each tying point) together weighted
  by the tying interpolations: N11, N22, N22. Transform the tensorial
  strain into the local strain.

  input:
  tx: the transformation to the local coordinates
  b11, b22, b12, b23, b13: the tying strain matrices
  N11, N22, N12: the tying strain shape functions

  output:
  B: the derivative of the strain w.r.t. the nodal displacements - this
  only writes in the in-plane and out-of-plane shear components (not
  the bending components.)
*/
template <int tying_order>
void add_tying_bmat(TacsScalar B[], const int num_nodes, const TacsScalar tx[],
                    const TacsScalar b11[], const TacsScalar b22[],
                    const TacsScalar b12[], const TacsScalar b23[],
                    const TacsScalar b13[], const double N11[],
                    const double N22[], const double N12[]) {
  // Variables stored by b11[nvars*(n + m*(order-1)) + row] Points stored
  const int n6pts = 6 * num_nodes;
  const int n3pts = 3 * num_nodes;

  for (int i = 0; i < num_nodes; i++) {
    for (int ii = 0; ii < 6; ii++) {
      TacsScalar g[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

      int row = 6 * i + ii;
      if (ii < 3) {
        int row3 = 3 * i + ii;
        // Use the interpolations
        for (int k = 0; k < tying_order * (tying_order - 1); k++) {
          g[0] += b11[n3pts * k + row3] * N11[k];
          g[1] += b22[n3pts * k + row3] * N22[k];
          g[3] += b23[n6pts * k + row] * N22[k];
          g[4] += b13[n6pts * k + row] * N11[k];
        }

        for (int k = 0; k < (tying_order - 1) * (tying_order - 1); k++) {
          g[5] += b12[n3pts * k + row3] * N12[k];
        }
      } else {
        // Use the interpolations
        for (int k = 0; k < tying_order * (tying_order - 1); k++) {
          g[3] += b23[n6pts * k + row] * N22[k];
          g[4] += b13[n6pts * k + row] * N11[k];
        }
      }

      // Transform the tensorial strain g, to the local axis
      TacsScalar s[6];
      Tensor::transform3DStress(s, g, tx);
      B[0] = s[0];
      B[1] = s[1];
      B[2] = 2.0 * s[5];
      B[6] = 2.0 * s[3];
      B[7] = 2.0 * s[4];
      B += 8;
    }
  }
}

/*
  Add the sensitivity of the tying strain to dstrain

  input:
  tx, dtx: the transformation and the derivative of the transformation

  g11, g22, g12, g23, g13: the tensorial displacement strain values at
  the tying points

  dg11, dg22, dg12, dg23, dg13: the derivatives of the tensorial
  displacement strain values at the tying points

  N11, N22, N12: the shape functions for interpolating the tying strain
*/
template <int tying_order>
void add_tying_strain_sens(TacsScalar strain[], TacsScalar dstrain[],
                           const TacsScalar tx[], const TacsScalar dtx[],
                           const TacsScalar g11[], const TacsScalar g22[],
                           const TacsScalar g12[], const TacsScalar g23[],
                           const TacsScalar g13[], const TacsScalar dg11[],
                           const TacsScalar dg22[], const TacsScalar dg12[],
                           const TacsScalar dg23[], const TacsScalar dg13[],
                           const double N11[], const double N22[],
                           const double N12[]) {
  // Apply the interpolation
  TacsScalar g[6];
  g[0] = g[1] = g[2] = g[3] = g[4] = g[5] = 0.0;

  for (int k = 0; k < tying_order * (tying_order - 1); k++) {
    g[0] += g11[k] * N11[k];
    g[1] += g22[k] * N22[k];
    g[3] += g23[k] * N22[k];
    g[4] += g13[k] * N11[k];
  }

  for (int k = 0; k < (tying_order - 1) * (tying_order - 1); k++) {
    g[5] += g12[k] * N12[k];
  }

  // Now compute the sensitivity of the transformation
  TacsScalar _dg[18 * tying_order * tying_order];

  // The number of derivatives in each of dg11, dg22, dg12 etc.
  const int nc = 3 * tying_order * tying_order;
  memset(_dg, 0, 18 * tying_order * tying_order * sizeof(TacsScalar));

  for (int k = 0; k < tying_order * (tying_order - 1); k++) {
    for (int j = 0; j < nc; j++) {
      _dg[6 * j] += dg11[nc * k + j] * N11[k];
      _dg[6 * j + 1] += dg22[nc * k + j] * N22[k];
      _dg[6 * j + 3] += dg23[nc * k + j] * N22[k];
      _dg[6 * j + 4] += dg13[nc * k + j] * N11[k];
    }
  }

  for (int k = 0; k < (tying_order - 1) * (tying_order - 1); k++) {
    for (int j = 0; j < nc; j++) {
      _dg[6 * j + 5] += dg12[nc * k + j] * N12[k];
    }
  }

  TacsScalar as[9];
  as[0] = tx[0] * g[0] + tx[1] * g[5] + tx[2] * g[4];
  as[3] = tx[3] * g[0] + tx[4] * g[5] + tx[5] * g[4];
  as[6] = tx[6] * g[0] + tx[7] * g[5] + tx[8] * g[4];

  as[1] = tx[0] * g[5] + tx[1] * g[1] + tx[2] * g[3];
  as[4] = tx[3] * g[5] + tx[4] * g[1] + tx[5] * g[3];
  as[7] = tx[6] * g[5] + tx[7] * g[1] + tx[8] * g[3];

  as[2] = tx[0] * g[4] + tx[1] * g[3] + tx[2] * g[2];
  as[5] = tx[3] * g[4] + tx[4] * g[3] + tx[5] * g[2];
  as[8] = tx[6] * g[4] + tx[7] * g[3] + tx[8] * g[2];

  strain[0] = as[0] * tx[0] + as[1] * tx[1] + as[2] * tx[2];
  strain[1] = as[3] * tx[3] + as[4] * tx[4] + as[5] * tx[5];
  strain[2] = 2.0 * (as[0] * tx[3] + as[1] * tx[4] + as[2] * tx[5]);

  strain[6] = 2.0 * (as[3] * tx[6] + as[4] * tx[7] + as[5] * tx[8]);
  strain[7] = 2.0 * (as[0] * tx[6] + as[1] * tx[7] + as[2] * tx[8]);

  TacsScalar *dg = _dg;

  // Now that the derivative dg has been computed, perform the
  // transformation for all components simultaneously.
  for (int k = 0; k < nc; k++) {
    TacsScalar das[9];
    das[0] = (dtx[0] * g[0] + dtx[1] * g[5] + dtx[2] * g[4] + tx[0] * dg[0] +
              tx[1] * dg[5] + tx[2] * dg[4]);
    das[3] = (dtx[3] * g[0] + dtx[4] * g[5] + dtx[5] * g[4] + tx[3] * dg[0] +
              tx[4] * dg[5] + tx[5] * dg[4]);
    das[6] = (dtx[6] * g[0] + dtx[7] * g[5] + dtx[8] * g[4] + tx[6] * dg[0] +
              tx[7] * dg[5] + tx[8] * dg[4]);

    das[1] = (dtx[0] * g[5] + dtx[1] * g[1] + dtx[2] * g[3] + tx[0] * dg[5] +
              tx[1] * dg[1] + tx[2] * dg[3]);
    das[4] = (dtx[3] * g[5] + dtx[4] * g[1] + dtx[5] * g[3] + tx[3] * dg[5] +
              tx[4] * dg[1] + tx[5] * dg[3]);
    das[7] = (dtx[6] * g[5] + dtx[7] * g[1] + dtx[8] * g[3] + tx[6] * dg[5] +
              tx[7] * dg[1] + tx[8] * dg[3]);

    das[2] = (dtx[0] * g[4] + dtx[1] * g[3] + dtx[2] * g[2] + tx[0] * dg[4] +
              tx[1] * dg[3] + tx[2] * dg[2]);
    das[5] = (dtx[3] * g[4] + dtx[4] * g[3] + dtx[5] * g[2] + tx[3] * dg[4] +
              tx[4] * dg[3] + tx[5] * dg[2]);
    das[8] = (dtx[6] * g[4] + dtx[7] * g[3] + dtx[8] * g[2] + tx[6] * dg[4] +
              tx[7] * dg[3] + tx[8] * dg[2]);

    // Add this
    dstrain[0] = (as[0] * dtx[0] + as[1] * dtx[1] + as[2] * dtx[2] +
                  das[0] * tx[0] + das[1] * tx[1] + das[2] * tx[2]);
    dstrain[1] = (as[3] * dtx[3] + as[4] * dtx[4] + as[5] * dtx[5] +
                  das[3] * tx[3] + das[4] * tx[4] + das[5] * tx[5]);
    dstrain[2] = 2.0 * (as[0] * dtx[3] + as[1] * dtx[4] + as[2] * dtx[5] +
                        das[0] * tx[3] + das[1] * tx[4] + das[2] * tx[5]);

    dstrain[6] = 2.0 * (as[3] * dtx[6] + as[4] * dtx[7] + as[5] * dtx[8] +
                        das[3] * tx[6] + das[4] * tx[7] + das[5] * tx[8]);
    dstrain[7] = 2.0 * (as[0] * dtx[6] + as[1] * dtx[7] + as[2] * dtx[8] +
                        das[0] * tx[6] + das[1] * tx[7] + das[2] * tx[8]);

    dstrain += 8;
    dtx += 9;
    dg += 6;
  }
}

template <int order, int tying_order>
void add_tying_bmat_sens(TacsScalar B[], TacsScalar dB[], const TacsScalar tx[],
                         const TacsScalar dtx[], const TacsScalar b11[],
                         const TacsScalar b22[], const TacsScalar b12[],
                         const TacsScalar b23[], const TacsScalar b13[],
                         const TacsScalar db11[], const TacsScalar db22[],
                         const TacsScalar db12[], const TacsScalar db23[],
                         const TacsScalar db13[], const double N11[],
                         const double N22[], const double N12[]) {
  // Variables stored by b11[nvars*(n + m*(order-1)) + row] Points stored
  const int n6pts = 6 * order * order;
  const int n3pts = 3 * order * order;

  for (int i = 0; i < order * order; i++) {
    for (int ii = 0; ii < 6; ii++) {
      TacsScalar g[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      TacsScalar dg[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

      int row = 6 * i + ii;

      if (ii < 3) {
        int row3 = 3 * i + ii;

        // Use the interpolations
        for (int k = 0; k < tying_order * (tying_order - 1); k++) {
          g[0] += b11[n3pts * k + row3] * N11[k];
          g[1] += b22[n3pts * k + row3] * N22[k];
          g[3] += b23[n6pts * k + row] * N22[k];
          g[4] += b13[n6pts * k + row] * N11[k];

          dg[0] += db11[n3pts * k + row3] * N11[k];
          dg[1] += db22[n3pts * k + row3] * N22[k];
          dg[3] += db23[n6pts * k + row] * N22[k];
          dg[4] += db13[n6pts * k + row] * N11[k];
        }

        for (int k = 0; k < (tying_order - 1) * (tying_order - 1); k++) {
          g[5] += b12[n3pts * k + row3] * N12[k];
          dg[5] += db12[n3pts * k + row3] * N12[k];
        }
      } else {
        // Use the interpolations
        for (int k = 0; k < tying_order * (tying_order - 1); k++) {
          g[3] += b23[n6pts * k + row] * N22[k];
          g[4] += b13[n6pts * k + row] * N11[k];
          dg[3] += db23[n6pts * k + row] * N22[k];
          dg[4] += db13[n6pts * k + row] * N11[k];
        }
      }

      // Transform the tensorial strain g, to the local axis
      TacsScalar s[6], ds[6];
      Tensor::transform3DStressSens(s, ds, g, dg, tx, dtx);
      B[0] = s[0];
      B[1] = s[1];
      B[2] = 2.0 * s[5];
      B[6] = 2.0 * s[3];
      B[7] = 2.0 * s[4];

      dB[0] = ds[0];
      dB[1] = ds[1];
      dB[2] = 2.0 * ds[5];
      dB[6] = 2.0 * ds[3];
      dB[7] = 2.0 * ds[4];

      B += 8;
      dB += 8;
    }
  }
}

TACS_END_NAMESPACE

#endif
