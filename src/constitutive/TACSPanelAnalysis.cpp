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

#ifdef TACS_USE_COMPLEX
#warning "TACSPanelAnalysis cannot be used with complex TACS"
#else

#include "TACSPanelAnalysis.h"

#include "FElibrary.h"
#include "tacslapack.h"

/*
  The following code uses a finite-strip panel model to compute the
  critical buckling loads and frequencies of a simply supported
  stiffened panel model. The critical buckling modes can be computed
  for either in-plane compression or shear loading. The panel geometry
  may be sheared into a rhombus in the x-y plane so that the panel
  edges are not perpendicular. For these cases, the stiffeners remain
  straight in the x-y plane. The sheared panels tend to make the
  positive and negative shear buckling modes unequal.

  This code is designed to handle either isotropic or composite cases
  but, the laminates must be balanced such that the in-plane stiffness
  takes the following form:

  A =
  [ A11,  A12,  0   ]
  [ A12,  A22,  0   ]
  [ 0,    0,    A66 ]

  Deviation from this stiffness will lead to less and less accurate
  results.  However, the contributions of the A16, A26 stiffness terms
  vary, depending on their relative magnitudes. Under certain
  conditions, it may be possible to use this analysis even if A16 and
  A26 are not precisely zero.

  We compute the response of the panels based on the following
  kinematics:

  u = sum_{n} u_n(y)*sin(n*pi/L*(x + beta*y))
  v = sum_{n} v_n(y)*sin(n*pi/L*(x + beta*y))
  w = sum_{n} w_n(y)*sin(n*pi/L*(x + beta*y))

  where u_n(y), v_n(y) and w_n(y) are polynomails in the displacement.
  Note that beta is defined as follows:

  beta = - tan(theta)

  such that the rib locations - or locations where simply-supported
  BCs are applied exist along x = - beta*y, and x = - beta*y + Lx

  We use classical lamination theory, rather than first-order shear
  deformation theory. All the constitutive objects used here are FSDT
  stiffness objects, but the shear components are not used. (This is
  for practical reasons, rather than requiring a full new set of
  stiffness objects just for this class.)
*/

/*
  This is the constructor for the TACSPanelAnalysis object.

  The initialization of this object proceeds in a few stages. First,
  the object is initialized with a call to this constructor. This call
  allocates the internal data structures required to store the panel
  information. Next, the geometry, segments and constitutive objects
  are initialized through a series of calls. Finally, once all the
  data has been passed to TACSPanelAnalysis, the call to finialize()
  is made to finalize data internally in TACSPanelAnalysis.

  Note that the way the geometry is specified is a projection onto the
  y-z plane. This projection is then sheared according to the
  parameter theta which controls the angle of shear. The geometry is
  specified as a list of nodes and connectivity, as well as a
  constitutive object for each segment.

  input:
  nnodes:    the number of nodes in the model
  nsegments: the number of segments between nodes
  nbeams:    the number of beams placed at nodes
  nmodes:    the number of modes to use in the compuation
  Lx:        the length of the panel in the x-direction
  theta:     the angle of geometric shearing
*/
TACSPanelAnalysis::TACSPanelAnalysis(int _nnodes, int _nsegments, int _nbeams,
                                     int _nmodes, TacsScalar _Lx,
                                     TacsScalar theta) {
  nnodes = _nnodes;
  nsegments = _nsegments;
  nbeams = _nbeams;
  nmodes = _nmodes;
  if (nmodes < 2) {
    nmodes = 2;
    fprintf(stderr, "TACSPanelAnalysis: Must use a minimum of two modes\n");
  }

  // Set the flag about whether to use the LAPACK eigenvalue solver
  // or not. It is slow b/c it computes the full spectrum, rather than
  // the selected spectrum, which the explicitly restarted Lanczos
  // code does - the default solver.
  use_lapack_eigensolver = 0;

  // Set the default lanczos subspace size - this can be
  // modified by calling setLanczosSubspaceSize - if your eigenvalue
  // problems are not convering, try modifying this, especially if you
  // require many eigenvalues
  lanczos_subspace_size = 100;

  // Set the eigenvalue tolerance
  lanczos_eigen_tol = 1e-12;

  // The slope of the line x = - beta*y, that defines the rib locations
  beta = 0.0;
  if (theta != 0.0) {
    beta = -tan(theta);
  }

  Xpts = new TacsScalar[2 * nnodes];
  panels = new FSDTStiffness *[nsegments];
  nodes = new int[2 * nsegments];

  beams = new EBStiffness *[nbeams];
  bnodes = new int[nbeams];

  segmentType = new int[nsegments];
  memset(segmentType, 0, nsegments * sizeof(int));

  numFailPoints = 0;
  for (int k = 0; k < MAX_NUM_FAIL_POINTS; k++) {
    failSegments[k] = -1;
    failNodes[k] = -1;
    failPointIsSkin[k] = 1;
  }

  // Set the length of the variable
  Lx = _Lx;
  LxDvNum = -1;
  LxLb = LxUb = 0.0;

  nDvGeo = 0;
  XptConst = NULL;
  XptLin = NULL;

  geoDvs = NULL;
  geoDvNums = NULL;
  geoLb = geoUb = NULL;

  for (int k = 0; k < 2 * nsegments; k++) {
    nodes[k] = -1;
  }

  for (int k = 0; k < nbeams; k++) {
    bnodes[k] = -1;
  }

  // Set the default start and end nodes and their boundary
  // conditions
  first_node = 0;
  last_node = nnodes - 1;

  // Note that the boundary conditions are bitwise operators
  first_node_bc = (4 | 8);
  last_node_bc = (4 | 8);

  memset(panels, 0, nsegments * sizeof(FSDTStiffness *));
  memset(beams, 0, nbeams * sizeof(EBStiffness *));
  memset(Xpts, 0, 2 * nnodes * sizeof(TacsScalar));

  // 4*nmodes variables for each node
  nvars = nnodes * (4 * nmodes);

  // Allocate space for the variables, and set their
  // initial values
  vars = new int[nvars];
  for (int k = 0; k < nvars; k++) {
    vars[k] = k;
  }

  // The Gauss points and weights
  numGauss = 2;
  gaussWts = FElibrary::gaussWts2;
  gaussPts = FElibrary::gaussPts2;

  // The number of bands stored in the matrix
  nband = -1;
}

/*
  Destroy the internal data allocated by TACSPanelAnalysis
*/
TACSPanelAnalysis::~TACSPanelAnalysis() {
  for (int k = 0; k < nsegments; k++) {
    if (panels[k]) {
      panels[k]->decref();
    }
  }
  delete[] panels;
  delete[] nodes;

  delete[] segmentType;

  for (int k = 0; k < nbeams; k++) {
    if (beams[k]) {
      beams[k]->decref();
    }
  }
  delete[] beams;
  delete[] bnodes;

  delete[] Xpts;
  delete[] vars;

  // Delete design variable information if it was allocated
  if (XptConst) {
    delete[] XptConst;
    delete[] XptLin;
    delete[] geoDvs;
    delete[] geoDvNums;
  }

  if (geoLb) {
    delete[] geoLb;
    delete[] geoUb;
  }
}

/*
  Set the flag for the eigensolver
*/
void TACSPanelAnalysis::setUseLapackEigensolver(int use_lapack) {
  use_lapack_eigensolver = use_lapack;
}

/*
  Set the size of the Lanczos subspace to use for the eigensolver
*/
void TACSPanelAnalysis::setLanczosSubspaceSize(int subspace_size) {
  lanczos_subspace_size = subspace_size;
}

/*
  Set the nodal locations.

  This code simply copies over the input data

  input:
  Xpts:    the nodal locatiosn
  npoints: the number of points in the y-z plane
*/
void TACSPanelAnalysis::setPoints(TacsScalar *_Xpts, int npoints) {
  memcpy(Xpts, _Xpts, 2 * nnodes * sizeof(TacsScalar));
}

/*
  Set a segment between two nodes

  input:
  seg:      the segment number
  seg_type: the type of segment - stiffener or skin
  stiff:    the constitutive object for this segment
  n1, n2:   the start and end nodes for this segment
*/
void TACSPanelAnalysis::setSegment(int seg, int seg_type, FSDTStiffness *stiff,
                                   int n1, int n2) {
  if (seg >= 0 && seg < nsegments) {
    if (stiff) {
      stiff->incref();
    }
    if (panels[seg]) {
      panels[seg]->decref();
    }
    segmentType[seg] = seg_type;
    panels[seg] = stiff;
    nodes[2 * seg] = n1;
    nodes[2 * seg + 1] = n2;
  }
}

/*
  Set the location of a beam

  input:
  beam:  the beam number
  stiff: the stiffness object associated with the beam
  n:     the node number at which to place the beam
*/
void TACSPanelAnalysis::setBeam(int beam, EBStiffness *stiff, int n) {
  if (beam >= 0 && beam < nbeams) {
    if (stiff) {
      stiff->incref();
    }
    if (beams[beam]) {
      beams[beam]->decref();
    }
    beams[beam] = stiff;
    bnodes[beam] = n;
  }
}

/*
  Set the boundary conditions on the displacements for the first node

  Note that there are 4 displacements and rotations per node: the x,
  y, z displacements as well as the rotation about x.

  input:
  first_node:     the node number to constrain
  first_node_bc:  bitwise integer indicating which variables to constrain
*/
void TACSPanelAnalysis::setFirstNodeBC(int _first_node, int _first_node_bc) {
  if (_first_node >= 0 && _first_node < nnodes) {
    first_node = _first_node;
    first_node_bc = _first_node_bc;
  }
}

/*
  Set the boundary conditions on the displacements for the last node

  Note that there are 4 displacements and rotations per node: the x,
  y, z displacements as well as the rotation about x.

  input:
  last_node:     the node number to constrain
  last_node_bc:  bitwise integer indicating which variables to constrain
*/
void TACSPanelAnalysis::setLastNodeBC(int _last_node, int _last_node_bc) {
  if (_last_node >= 0 && _last_node < nnodes) {
    last_node = _last_node;
    last_node_bc = _last_node_bc;
  }
}

/*
  Set the design variable-dependence of the geometry.

  The geometry must depend linearly on the design variables. This
  can actually produce very general geometry modifications and covers
  most common design modifications.

  The nodal locations are computed as follows:

  Xpts = XptConst + XptLin*x[:nDvGeo]

  where x[:nDvGeo] are the geometric design variables.

  input:
  XptConst:  the constant terms
  XptLin:    the linear terms in the geometry
  geoDvs:    the initial value of the geometric design variables
  dvNums:    the geometric design variable numbers
  nDvGeo:    the number of geometric design variables
*/
void TACSPanelAnalysis::setGeoDesignDependence(TacsScalar *_XptConst,
                                               TacsScalar *_XptLin,
                                               TacsScalar *_geoDvs,
                                               int *_dvNums, int _nDvGeo) {
  if (!XptConst) {
    nDvGeo = _nDvGeo;
    geoDvNums = new int[nDvGeo];
    geoDvs = new TacsScalar[nDvGeo];
    geoLb = new TacsScalar[nDvGeo];
    geoUb = new TacsScalar[nDvGeo];

    XptConst = new TacsScalar[2 * nnodes];
    XptLin = new TacsScalar[2 * nnodes * nDvGeo];

    memcpy(geoDvNums, _dvNums, nDvGeo * sizeof(int));
    memcpy(XptConst, _XptConst, 2 * nnodes * sizeof(TacsScalar));
    memcpy(XptLin, _XptLin, 2 * nnodes * nDvGeo * sizeof(TacsScalar));
    memcpy(geoDvs, _geoDvs, nDvGeo * sizeof(TacsScalar));
    memset(geoLb, 0, nDvGeo * sizeof(TacsScalar));
    memset(geoUb, 0, nDvGeo * sizeof(TacsScalar));
  }
}

/*
  Set the lower and upper bounds on the design variables

  input:
  lb, ub:  the lower/upper bounds on the geometric design variables
*/
void TACSPanelAnalysis::setGeoDVBounds(TacsScalar *lb, TacsScalar *ub) {
  if (!geoLb || !geoUb) {
    geoLb = new TacsScalar[nDvGeo];
    geoUb = new TacsScalar[nDvGeo];
  }

  memcpy(geoLb, lb, nDvGeo * sizeof(TacsScalar));
  memcpy(geoUb, ub, nDvGeo * sizeof(TacsScalar));
}

/*
  Set the design variable number for Lx and the design variable
  bounds

  The length of the panel is typically a design variable of interest.
  This code adds it as a design variable and sets the bounds on Lx as
  well.

  input:
  LxDvNum:  the design variable number
  LxLb:     lower bound on the panel length
  LxUb:     upper bound on the panel length
*/
void TACSPanelAnalysis::setLxDvBounds(int _LxDvNum, TacsScalar _LxLb,
                                      TacsScalar _LxUb) {
  LxDvNum = _LxDvNum;
  LxLb = _LxLb;
  LxUb = _LxUb;
}

/*
  This function must be called after all the segments, beams and
  design variable information has been set, but before any analysis is
  performed.

  This function performs two tasks:

  1. Re-count the number of variables taking the constrained modes
  into account

  2. Count up the number of design variables - both material and
  geometric design variables
*/
void TACSPanelAnalysis::initialize() {
  // Perform a Cuthill McKee reordering of the nodes to minimize
  // the bandwidth of the matrix
  int *work = new int[2 * nnodes];
  int *node_order = &work[0];
  int *node_flag = &work[nnodes];

  for (int k = 0; k < nnodes; k++) {
    node_order[k] = -1;
    node_flag[k] = 0;  // Indicate whether the node is ordered
  }

  // Mark the segment of nodes that
  int first = 0, last = 1;
  node_order[0] = 0;
  node_flag[0] = 1;

  // Note that this is an n^2 algorithm (worst case). If the number
  // of segments becomes too large, this will slow down considerably
  // but we usually have O(50) segments, so this shouldn't be too bad.
  while (last < nnodes) {
    int new_nodes = 0;
    for (int n = 0; n < nsegments; n++) {
      for (int k = first; k < last; k++) {
        // Get the node from the list
        int node = node_order[k];

        // Check if the other node for this segment has not yet
        // been ordered, or is not yet set
        if (nodes[2 * n] == node && !node_flag[nodes[2 * n + 1]]) {
          node_order[last + new_nodes] = nodes[2 * n + 1];
          node_flag[nodes[2 * n + 1]] = 1;
          new_nodes++;
        } else if (nodes[2 * n + 1] == node && !node_flag[nodes[2 * n]]) {
          node_order[last + new_nodes] = nodes[2 * n];
          node_flag[nodes[2 * n]] = 1;
          new_nodes++;
        }
      }
    }

    if (new_nodes == 0) {
      // If no new nodes are added, then we have two disconnected
      // segments of the mesh - we can continue, but this
      // should never happen
      for (int k = 0; k < nnodes; k++) {
        if (!node_flag[k]) {
          node_order[last] = k;
          node_flag[k] = 1;
          first = last;
          last += 1;
          break;
        }
      }
    } else {
      first = last;
      last += new_nodes;
    }
  }

  // Now, order the nodes using the Cuthill McKee ordering
  nvars = 0;
  for (int j = 0; j < nnodes; j++) {
    int node = node_order[j];

    if (node == first_node || node == last_node) {
      // Check the boundary conditions
      int bc = first_node_bc;
      if (node == last_node) {
        bc = last_node_bc;
      }

      // Set the boundary conditions
      for (int j = 0; j < 4; j++) {
        if (bc & (1 << j)) {
          for (int k = 0; k < nmodes; k++) {
            vars[4 * nmodes * node + 4 * k + j] = -1;
          }
        } else {
          for (int k = 0; k < nmodes; k++) {
            vars[4 * nmodes * node + 4 * k + j] = nvars;
            nvars++;
          }
        }
      }
    } else {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < nmodes; k++) {
          vars[4 * nmodes * node + 4 * k + j] = nvars;
          nvars++;
        }
      }
    }
  }

  // Compute the maximum bandwidth of the matrix
  nband = 1;
  for (int k = 0; k < nsegments; k++) {
    int n1 = nodes[2 * k];
    int n2 = nodes[2 * k + 1];

    int bandwidth = node_order[n2] - node_order[n1];
    if (bandwidth < 0) {
      bandwidth *= -1;
    }

    if (bandwidth + 1 > nband) {
      nband = bandwidth + 1;
    }
  }

  // Adjust the bandwidth for the number of nodes
  nband = 4 * nmodes * nband - 1;

  delete[] work;
}

/*
  Set the design variable values in both the constitutive objects as
  well as the geometric design variables

  input:
  dvs:    the design variable values
  numDVs: the length of the design variable array
*/
void TACSPanelAnalysis::setDesignVars(const TacsScalar dvs[], int numDVs) {
  for (int k = 0; k < nsegments; k++) {
    panels[k]->setDesignVars(dvs, numDVs);
  }

  for (int k = 0; k < nbeams; k++) {
    beams[k]->setDesignVars(dvs, numDVs);
  }

  // Set the geometric design variables
  for (int k = 0; k < nDvGeo; k++) {
    if (geoDvNums[k] < numDVs) {
      geoDvs[k] = dvs[geoDvNums[k]];
    }
  }

  if (LxDvNum >= 0 && LxDvNum < numDVs) {
    Lx = dvs[LxDvNum];
  }

  updateGeometry();
}

/*
  Get the design variable values. Retrieve the values of the design
  variables from the internal stiffness objects as well as the
  geometric and panel length variables.

  input:
  numDVs:  the length of the design variable array

  output:
  dvs:     the internal design variable values
*/
void TACSPanelAnalysis::getDesignVars(TacsScalar dvs[], int numDVs) const {
  for (int k = 0; k < nsegments; k++) {
    panels[k]->getDesignVars(dvs, numDVs);
  }

  for (int k = 0; k < nbeams; k++) {
    beams[k]->getDesignVars(dvs, numDVs);
  }

  // Set the geometric design variables
  for (int k = 0; k < nDvGeo; k++) {
    if (geoDvNums[k] < numDVs) {
      dvs[geoDvNums[k]] = geoDvs[k];
    }
  }

  if (LxDvNum >= 0 && LxDvNum < numDVs) {
    dvs[LxDvNum] = Lx;
  }
}

/*
  Get the design variable range

  Retrieve the range of allowable design variable values from the
  stiffness objects, geometric variables and panel length variable.

  input:
  numDVs:  the length of the bound arrays

  output:
  lb, ub:  the lower and upper variable bounds
*/
void TACSPanelAnalysis::getDesignVarRange(TacsScalar lb[], TacsScalar ub[],
                                          int numDVs) const {
  for (int k = 0; k < nsegments; k++) {
    panels[k]->getDesignVarRange(lb, ub, numDVs);
  }

  for (int k = 0; k < nbeams; k++) {
    beams[k]->getDesignVarRange(lb, ub, numDVs);
  }

  for (int k = 0; k < nDvGeo; k++) {
    if (geoDvNums[k] < numDVs) {
      lb[geoDvNums[k]] = geoLb[k];
      ub[geoDvNums[k]] = geoUb[k];
    }
  }

  if (LxDvNum >= 0 && LxDvNum < numDVs) {
    lb[LxDvNum] = LxLb;
    ub[LxDvNum] = LxUb;
  }
}

/*
  Update the geometry of the cross section based on the geometric
  design variables.
*/
void TACSPanelAnalysis::updateGeometry() {
  if (nDvGeo > 0) {
    for (int k = 0; k < 2 * nnodes; k++) {
      Xpts[k] = XptConst[k];
      for (int j = 0; j < nDvGeo; j++) {
        Xpts[k] += XptLin[2 * nnodes * j + k] * geoDvs[j];
      }
    }
  }
}

/*
  Compute the mass per unit area of the panel

  The panel mass per unit area is calculated by computing the mass of
  all segments in the y-z plane and dividing by the length of the
  panel in the y-direction.

  The length of the panel in the y-direction is determined by setting
  start and end nodes that define the outer-extent of the panel (e.g.
  the first and last nodes traveling along the y-direction.)

  returns:
  the panel mass per unit area
*/
TacsScalar TACSPanelAnalysis::computePtMass() {
  TacsScalar ptmass = 0.0;

  for (int k = 0; k < nsegments; k++) {
    int n1 = nodes[2 * k];
    int n2 = nodes[2 * k + 1];

    TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
    TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
    TacsScalar Le = sqrt(c * c + s * s);

    double pt[2] = {0.0, 0.0};
    TacsScalar mass[2];
    panels[k]->getPointwiseMass(pt, mass);
    ptmass += mass[0] * Le;
  }

  TacsScalar cy = Xpts[2 * last_node] - Xpts[2 * first_node];
  TacsScalar sy = Xpts[2 * last_node + 1] - Xpts[2 * first_node + 1];
  TacsScalar Ly = sqrt(cy * cy + sy * sy);
  ptmass /= Ly;

  return ptmass;
}

/*
  Add the derivative of the mass per unit area of the panel with
  respect to the given design variable vector

  This code works by first determining which design variable dvNum
  corresponds to through a call to getDesignVarIndex, then
  computing the derivative w.r.t. that design variable.

  input:
  scale:  scale the design sensitivity by this value

  in/out:
  fdvSens:  the derivative array
  numDVs:   the length of the derivative array
*/
void TACSPanelAnalysis::addPtMassDVSens(TacsScalar scale, TacsScalar fdvSens[],
                                        int numDVs) {
  // Compute the full panel length in the y-direction
  TacsScalar cy = Xpts[2 * last_node] - Xpts[2 * first_node];
  TacsScalar sy = Xpts[2 * last_node + 1] - Xpts[2 * first_node + 1];
  TacsScalar invLy = 1.0 / sqrt(cy * cy + sy * sy);

  TacsScalar ptmass = 0.0;
  for (int k = 0; k < nsegments; k++) {
    int n1 = nodes[2 * k];
    int n2 = nodes[2 * k + 1];

    TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
    TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
    TacsScalar Le = sqrt(c * c + s * s);

    double pt[2] = {0.0, 0.0};
    TacsScalar alpha[2] = {scale * Le * invLy, 0.0};
    panels[k]->addPointwiseMassDVSens(pt, alpha, fdvSens, numDVs);

    // Add the contribution to the pointmass
    TacsScalar mass[2];
    panels[k]->getPointwiseMass(pt, mass);
    ptmass += mass[0] * Le;
  }

  // Scale the total mass/x-direction length by the y-dimension
  // of the panel
  ptmass *= invLy;

  // Now, loop over all of the geometric design variables
  for (int dv = 0; dv < nDvGeo; dv++) {
    // Compute the sensitivity from the y-direction
    TacsScalar scy = (XptLin[2 * nnodes * dv + 2 * last_node] -
                      XptLin[2 * nnodes * dv + 2 * first_node]);
    TacsScalar ssy = (XptLin[2 * nnodes * dv + 2 * last_node + 1] -
                      XptLin[2 * nnodes * dv + 2 * first_node + 1]);
    TacsScalar sLy = (cy * scy + sy * ssy) * invLy;

    // Loop over each segment and compute the derivative
    TacsScalar sptmass = 0.0;
    for (int k = 0; k < nsegments; k++) {
      int n1 = nodes[2 * k];
      int n2 = nodes[2 * k + 1];

      TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
      TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
      TacsScalar sc =
          (XptLin[2 * nnodes * dv + 2 * n2] - XptLin[2 * nnodes * dv + 2 * n1]);
      TacsScalar ss = (XptLin[2 * nnodes * dv + 2 * n2 + 1] -
                       XptLin[2 * nnodes * dv + 2 * n1 + 1]);
      TacsScalar Le = sqrt(c * c + s * s);
      TacsScalar sLe = (c * sc + s * ss) / Le;

      double pt[2] = {0.0, 0.0};
      TacsScalar mass[2];
      panels[k]->getPointwiseMass(pt, mass);
      sptmass += mass[0] * sLe;
    }
    sptmass *= invLy;
    sptmass -= ptmass * invLy * sLy;

    if (geoDvNums[dv] >= 0 && geoDvNums[dv] < numDVs) {
      fdvSens[geoDvNums[dv]] = scale * sptmass;
    }
  }
}

/*
  Set the locations and segments to test for failure

  This code sets the internal data that is used to evaluate a failure
  constraint. The failure constraint is evaluated by aggregating the
  failure criteria at a series of points in the stiffness
  object. These points are offset from the neutral surface in a
  predefined manner and are located at the nodes of the finite-strip
  model.

  input:
  npoints:    the number of points to test < MAX_NUM_FAIL_POINTS
  fail_segs:  an array of the failure segments (index of stiffness objs)
  fail_nodes: an array of node numbers to test the failure
  is_skin:    is the segment a skin segment?
*/
void TACSPanelAnalysis::setFailurePoints(int npoints, int fail_segs[],
                                         int fail_nodes[], int is_skin[]) {
  numFailPoints = npoints;
  if (numFailPoints > MAX_NUM_FAIL_POINTS) {
    numFailPoints = MAX_NUM_FAIL_POINTS;
  }

  memcpy(failSegments, fail_segs, numFailPoints * sizeof(int));
  memcpy(failNodes, fail_nodes, numFailPoints * sizeof(int));
  memcpy(failPointIsSkin, is_skin, numFailPoints * sizeof(int));

  for (int k = 0; k < numFailPoints; k++) {
    if (failSegments[k] < 0 || failSegments[k] >= nsegments) {
      fprintf(stderr,
              "TACSPanelAnalysis::setFailurePoints() \
Segment out of range\n");
      failSegments[k] = 0;
    }
    if (failNodes[k] < 0 || failNodes[k] >= nnodes) {
      fprintf(stderr,
              "TACSPanelAnalysis::setFailurePoints() \
Node out of range\n");
      failNodes[k] = 0;
    }
  }
}

/*
  Compute the failure function over the segments

  The strain in the local segment coordinate system can be obtained
  by the strain in the global coordinate system as follows:

  ex'  = ex
  ey'  = ey * cos(t)**2 + ez * sin(t)**2 + gyz * sin(t)*cos(t)
  gxy' = gxy * cos(t) + gxz * sin(t)

  The local in-plane components of the panel strain tensor are given
  as follows:

  e[0] = strain[0] + z*strain[3];
  e[1] = (strain[1] + z*strain[4])*cos(t)**2
  e[2] = (strain[2] + z*strain[5])*cos(t)

  The bending strain requires an additional multiplication by cos(t)
  because they represent the local through-thickness variation of the
  strain. This local through-thickness direction varies as the angle
  changes. As a result the local bending strain is:

  e[3] = strain[3]*cos(t)
  e[4] = strain[4]*cos(t)**3
  e[5] = strain[5]*cos(t)**2
*/
void TACSPanelAnalysis::failure(const TacsScalar strain[], TacsScalar fail[],
                                int nfail) {
  const double pt[3] = {0.0, 0.0, 0.0};

  for (int k = 0; k < numFailPoints && k < nfail; k++) {
    int seg = failSegments[k];
    if (failPointIsSkin[k]) {
      panels[seg]->failure(pt, strain, &fail[k]);
    } else {
      // Compute the strain at the midpoint of the segment
      TacsScalar e[8];
      TacsScalar z = Xpts[2 * failNodes[k] + 1];

      int n1 = nodes[2 * seg];
      int n2 = nodes[2 * seg + 1];
      TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
      TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
      TacsScalar Le = sqrt(c * c + s * s);
      c = c / Le;
      s = s / Le;

      e[0] = (strain[0] + z * strain[3]);
      e[1] = e[2] = 0.0;

      e[3] = strain[3] * c;
      e[4] = e[5] = e[6] = e[7] = 0.0;

      panels[seg]->failure(pt, e, &fail[k]);
    }
  }
}

/*
  Compute the derivative of the failure function with respect to the
  given design variable number. Add the result times the weight vector
  to the array failDVSens.

  This function can be used to compute the derivative of an aggregated
  failure function w.r.t. the design variables.

  input:
  dvNum:   the design variable number
  strain:  the strain at the point in the constitutive object
  weights: the weight associated with each failure point
  nfail:   the number of failure points = len(weights)

  output:
  failDVSens:  the derivative of the failure function w.r.t. design vars
*/
void TACSPanelAnalysis::addFailureDVSens(const TacsScalar strain[],
                                         const TacsScalar weights[], int nfail,
                                         TacsScalar fdvSens[], int dvLen) {
  const double pt[2] = {0.0, 0.0};

  for (int dv = 0; dv < nDvGeo; dv++) {
    TacsScalar failDVSens = 0.0;
    for (int k = 0; k < numFailPoints && k < nfail; k++) {
      int seg = failSegments[k];
      if (!failPointIsSkin[k]) {
        // Compute the strain at the midpoint of the segment
        TacsScalar e[8];
        TacsScalar z = Xpts[2 * failNodes[k] + 1];
        TacsScalar sz = XptLin[2 * nnodes * dv + 2 * failNodes[k] + 1];

        // Compute the sensitivity of the locations
        int n1 = nodes[2 * seg];
        int n2 = nodes[2 * seg + 1];
        TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
        TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
        TacsScalar sc = (XptLin[2 * nnodes * dv + 2 * n2] -
                         XptLin[2 * nnodes * dv + 2 * n1]);
        TacsScalar ss = (XptLin[2 * nnodes * dv + 2 * n2 + 1] -
                         XptLin[2 * nnodes * dv + 2 * n1 + 1]);
        TacsScalar Le = sqrt(c * c + s * s);
        TacsScalar sLe = (c * sc + s * ss) / Le;
        sc = (sc * Le - sLe * c) / (Le * Le);
        c = c / Le;

        e[0] = (strain[0] + z * strain[3]);
        e[1] = e[2] = 0.0;

        e[3] = strain[3] * c;
        e[4] = e[5] = e[6] = e[7] = 0.0;

        TacsScalar eSens[8];
        panels[seg]->failureStrainSens(pt, e, eSens);
        failDVSens += weights[k] *
                      (eSens[0] * sz * strain[3] + eSens[3] * sc * strain[3]);
      }
    }

    if (geoDvNums[dv] >= 0 && geoDvNums[dv] < dvLen) {
      fdvSens[geoDvNums[dv]] = failDVSens;
    }
  }

  // Loop over all the node locations and add the derivative
  // of the material design variables
  for (int k = 0; k < numFailPoints && k < nfail; k++) {
    int seg = failSegments[k];
    if (failPointIsSkin[k]) {
      panels[seg]->addFailureDVSens(pt, strain, weights[k], fdvSens, dvLen);
    } else {
      // Compute the strain at the midpoint of the segment
      TacsScalar e[8];
      TacsScalar z = Xpts[2 * failNodes[k] + 1];

      int n1 = nodes[2 * seg];
      int n2 = nodes[2 * seg + 1];
      TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
      TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
      TacsScalar Le = sqrt(c * c + s * s);
      c = c / Le;
      s = s / Le;

      e[0] = (strain[0] + z * strain[3]);
      e[1] = e[2] = 0.0;

      e[3] = strain[3] * c;
      e[4] = e[5] = e[6] = e[7] = 0.0;

      panels[seg]->addFailureDVSens(pt, e, weights[k], fdvSens, dvLen);
    }
  }
}

/*
  Compute the derivative of the failure function with respect to the
  input strain. This takes into account the weights from the KS
  aggregation of all the failure points.

  input:
  strain:  the strain at the parametric point in the constitutive obj.
  weights: the weight associated with each failure point
  nfail:   the number of failure points = len(weights)

  output:
  failSens: the derivative of the aggregated failure loads w.r.t. strain
*/
void TACSPanelAnalysis::failureStrainSens(const TacsScalar strain[],
                                          const TacsScalar weights[], int nfail,
                                          TacsScalar failSens[]) {
  const double pt[3] = {0.0, 0.0, 0.0};

  failSens[0] = failSens[1] = failSens[2] = failSens[3] = 0.0;
  failSens[4] = failSens[5] = failSens[6] = failSens[7] = 0.0;

  for (int k = 0; k < numFailPoints && k < nfail; k++) {
    int seg = failSegments[k];
    if (failPointIsSkin[k]) {
      TacsScalar eSens[8];
      panels[seg]->failureStrainSens(pt, strain, eSens);

      failSens[0] += weights[k] * eSens[0];
      failSens[1] += weights[k] * eSens[1];
      failSens[2] += weights[k] * eSens[2];

      failSens[3] += weights[k] * eSens[3];
      failSens[4] += weights[k] * eSens[4];
      failSens[5] += weights[k] * eSens[5];

      failSens[6] += weights[k] * eSens[6];
      failSens[7] += weights[k] * eSens[7];
    } else {
      // Compute the strain at the midpoint of the segment
      TacsScalar e[8];
      int seg = failSegments[k];
      TacsScalar z = Xpts[2 * failNodes[k] + 1];

      int n1 = nodes[2 * seg];
      int n2 = nodes[2 * seg + 1];
      TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
      TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
      TacsScalar Le = sqrt(c * c + s * s);
      c = c / Le;
      s = s / Le;

      e[0] = (strain[0] + z * strain[3]);
      e[1] = e[2] = 0.0;

      e[3] = strain[3] * c;
      e[4] = e[5] = e[6] = e[7] = 0.0;

      TacsScalar eSens[8];
      panels[seg]->failureStrainSens(pt, e, eSens);

      failSens[0] += weights[k] * eSens[0];
      failSens[3] += weights[k] * (z * eSens[0] + c * eSens[3]);
    }
  }
}

/*
  Compute the smeared stiffness of the panel that accounts for both
  the stiffness of the panel itself and the stiffeners.

  The stiffnesses are defined as follows:

  A = 1/Ly*( sum_i A_i*L_i )

  B = 1/Ly*( sum_i 0.5*A_i*(z1 + z2)*L_i^2 + B_i cos(t)*L_i )

  D = 1/Ly*( sum_i 1/3*A_i*(z1^2 + z1*z2 + z2^2)*L_i^3 +
                   B_i*cos(t)*(z1 + z2)*L_i^2 +
                   D_i*cos(t)^2*L_i )

  output:
  A, B, D, As: the in-plane, bending-stretching, bending and
  out-of-plane shear stiffnesses for the smeared stiffness object
*/
void TACSPanelAnalysis::computeStiffness(TacsScalar A[], TacsScalar B[],
                                         TacsScalar D[], TacsScalar As[]) {
  double pt[3] = {0.0, 0.0, 0.0};

  for (int k = 0; k < 6; k++) {
    A[k] = B[k] = D[k] = 0.0;
  }
  As[0] = As[1] = As[2] = 0.0;

  for (int k = 0; k < nsegments; k++) {
    int n1 = nodes[2 * k];
    int n2 = nodes[2 * k + 1];

    TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
    TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
    TacsScalar Le = sqrt(c * c + s * s);
    c = c / Le;

    double kcorr = 5.0 / 6.0;
    TacsScalar At[6], Bt[6], Dt[6], Ast[3];
    ;
    panels[k]->getStiffness(pt, At, Bt, Dt, Ast);

    if (segmentType[k] == STIFFENER_SEGMENT) {
      TacsScalar A11 = (At[0] * At[3] - At[1] * At[1]) / At[3];

      TacsScalar z1 = Xpts[2 * n1 + 1];
      TacsScalar z2 = Xpts[2 * n2 + 1];

      A[0] += Le * A11;
      B[0] += 0.5 * (z1 + z2) * Le * A11;
      D[0] += (z1 * z1 + z1 * z2 + z2 * z2) * Le * A11 / 3.0;

      As[0] += kcorr * At[5] * Le;
    } else {
      for (int i = 0; i < 6; i++) {
        A[i] += At[i] * Le;
        B[i] += Bt[i] * Le;
        D[i] += Dt[i] * Le;
      }

      for (int i = 0; i < 3; i++) {
        As[i] += Ast[i] * Le;
      }
    }
  }

  TacsScalar cy = Xpts[2 * last_node] - Xpts[2 * first_node];
  TacsScalar sy = Xpts[2 * last_node + 1] - Xpts[2 * first_node + 1];
  TacsScalar invLy = 1.0 / sqrt(cy * cy + sy * sy);

  for (int i = 0; i < 6; i++) {
    A[i] *= invLy;
    B[i] *= invLy;
    D[i] *= invLy;
  }

  for (int i = 0; i < 3; i++) {
    As[i] *= invLy;
  }
}

/*
  Compute the derivative of the smeared stiffness of the panel
  w.r.t. the design variable.

  intput
  dvNum: the design variable number

  output:
  A, B, D, As: the derivative of the in-plane, bending-stretching,
  bending and out-of-plane shear stiffnesses for the smeared stiffness
  object w.r.t. the design variables
*/
/*
void TACSPanelAnalysis::computeStiffnessDVSens( int dvNum,
                                                TacsScalar sA[],
                                                TacsScalar sB[],
                                                TacsScalar sD[],
                                                TacsScalar sAs[] ){
  double pt[3] = {0.0, 0.0, 0.0};
  TacsScalar A[6], B[6], D[6], As[3];

  for ( int k = 0; k < 6; k++ ){
    A[k] = B[k] = D[k] = 0.0;
    sA[k] = sB[k] = sD[k] = 0.0;
  }
  As[0] = As[1] = As[2] = 0.0;
  sAs[0] = sAs[1] = sAs[2] = 0.0;

  int queryIndex = getDesignVarIndex(dvNum);

  if (queryIndex >= 0){
    if (designVarTypes[queryIndex] == 1){
      // Determine the index of the geometric design variable
      int dv = 0;
      for ( int j = 0; j < nDvGeo; j++ ){
        if (geoDvNums[j] == dvNum){
          dv = j;
          break;
        }
      }

      for ( int k = 0; k < nsegments; k++ ){
        int n1 = nodes[2*k];
        int n2 = nodes[2*k+1];

        TacsScalar c = (Xpts[2*n2] - Xpts[2*n1]);
        TacsScalar s = (Xpts[2*n2+1] - Xpts[2*n1+1]);
        TacsScalar sc = (XptLin[2*nnodes*dv + 2*n2] -
                         XptLin[2*nnodes*dv + 2*n1]);
        TacsScalar ss = (XptLin[2*nnodes*dv + 2*n2+1] -
                         XptLin[2*nnodes*dv + 2*n1+1]);
        TacsScalar Le = sqrt(c*c + s*s);
        TacsScalar sLe = (c*sc + s*ss)/Le;
        sc = (sc*Le - sLe*c)/(Le*Le);
        c = c/Le;

        double kcorr = 5.0/6.0;
        TacsScalar At[6], Bt[6], Dt[6], Ast[3];
        panels[k]->getStiffness(pt, At, Bt, Dt, Ast);

        if (segmentType[k] == STIFFENER_SEGMENT){
          TacsScalar A11 = (At[0]*At[3] - At[1]*At[1])/At[3];

          TacsScalar z1 = Xpts[2*n1+1];
          TacsScalar z2 = Xpts[2*n2+1];
          TacsScalar sz1 = XptLin[2*nnodes*dv + 2*n1+1];
          TacsScalar sz2 = XptLin[2*nnodes*dv + 2*n2+1];

          A[0] += Le*A11;
          B[0] += 0.5*(z1 + z2)*Le*A11;
          D[0] += (z1*z1 + z1*z2 + z2*z2)*Le*A11/3.0;

          As[0] += kcorr*At[5]*Le;

          sA[0] += sLe*A11;
          sB[0] += 0.5*((sz1 + sz2)*Le + (z1 + z2)*sLe)*A11;
          sD[0] += ((2.0*z1*sz1 + sz1*z2 + z1*sz2 + 2.0*z2*sz2)*Le +
                    (z1*z1 + z1*z2 + z2*z2)*sLe)*A11/3.0;

          sAs[0] += kcorr*At[5]*sLe;
        }
        else {
          for ( int i = 0; i < 6; i++ ){
            A[i] += At[i]*Le;
            B[i] += Bt[i]*Le;
            D[i] += Dt[i]*Le;

            sA[i] += At[i]*sLe;
            sB[i] += Bt[i]*sLe;
            sD[i] += Dt[i]*sLe;
          }

          for ( int i = 0; i < 3; i++ ){
            As[i] += Le*Ast[i];
            sAs[i] += Ast[i]*sLe;
          }
        }
      }

      TacsScalar cy = Xpts[2*last_node] - Xpts[2*first_node];
      TacsScalar sy = Xpts[2*last_node+1] - Xpts[2*first_node+1];
      TacsScalar scy = (XptLin[2*nnodes*dv + 2*last_node] -
                        XptLin[2*nnodes*dv + 2*first_node]);
      TacsScalar ssy = (XptLin[2*nnodes*dv + 2*last_node+1] -
                        XptLin[2*nnodes*dv + 2*first_node+1]);
      TacsScalar invLy = 1.0/sqrt(cy*cy + sy*sy);
      TacsScalar sinvLy = - invLy*invLy*invLy*(scy*cy + ssy*sy);

      for ( int i = 0; i < 6; i++ ){
        sA[i] = invLy*sA[i] + sinvLy*A[i];
        sB[i] = invLy*sB[i] + sinvLy*B[i];
        sD[i] = invLy*sD[i] + sinvLy*D[i];
      }

      for ( int i = 0; i < 3; i++ ){
        sAs[i] = invLy*sAs[i] + sinvLy*As[i];
      }
    }
    else if (designVarTypes[queryIndex] == 0){
      for ( int k = 0; k < nsegments; k++ ){
        int n1 = nodes[2*k];
        int n2 = nodes[2*k+1];

        TacsScalar c = (Xpts[2*n2] - Xpts[2*n1]);
        TacsScalar s = (Xpts[2*n2+1] - Xpts[2*n1+1]);
        TacsScalar Le = sqrt(c*c + s*s);
        c = c/Le;

        double kcorr = 5.0/6.0;
        TacsScalar At[6], Bt[6], Dt[6], Ast[3];
        TacsScalar sAt[6], sBt[6], sDt[6], sAst[3];
        panels[k]->getStiffness(pt, At, Bt, Dt, Ast);
        panels[k]->getStiffnessDVSens(dvNum, pt, sAt, sBt, sDt, sAst);

        if (segmentType[k] == STIFFENER_SEGMENT){
          TacsScalar A11 = (At[0]*At[3] - At[1]*At[1])/At[3];
          TacsScalar sA11 = (At[0]*sAt[3] + sAt[0]*At[3] -
                             2.0*At[1]*sAt[1])/At[3];
          sA11 -= (A11/At[3])*sAt[3];

          TacsScalar z1 = Xpts[2*n1+1];
          TacsScalar z2 = Xpts[2*n2+1];

          sA[0] += Le*sA11;
          sB[0] += 0.5*(z1 + z2)*Le*sA11;
          sD[0] += (z1*z1 + z1*z2 + z2*z2)*Le*sA11/3.0;

          sAs[0] += kcorr*sAt[5]*Le;
        }
        else {
          for ( int i = 0; i < 6; i++ ){
            sA[i] += sAt[i]*Le;
            sB[i] += sBt[i]*Le;
            sD[i] += sDt[i]*Le;
          }

          for ( int i = 0; i < 3; i++ ){
            sAs[i] += sAst[i]*Le;
          }
        }
      }

      TacsScalar cy = Xpts[2*last_node] - Xpts[2*first_node];
      TacsScalar sy = Xpts[2*last_node+1] - Xpts[2*first_node+1];
      TacsScalar invLy = 1.0/sqrt(cy*cy + sy*sy);

      for ( int i = 0; i < 6; i++ ){
        sA[i] *= invLy;
        sB[i] *= invLy;
        sD[i] *= invLy;
      }

      for ( int i = 0; i < 3; i++ ){
        sAs[i] *= invLy;
      }
    }
  }
}
*/

/*
  Compute a few select eigenvalues that solve the eigenvalue problem

  A*x = lambda*B*x

  using an explicitly restarted Lanczos Method. For this problem, B
  must be symmetric positive definite. We then transform the
  eigenproblem using a Cholesky factorization as follows:

  LB*LB^{T} = B   or   UB^{T}*UB = B

  such that the modified eigenproblem is:

  LB^{-1}*(A - s*B)*LB^{-T}*y = lambda*y   or
  UB^{-T}*(A - s*B)*UB^{-1}*y = lambda*y

  where x = LB^{-T}*y   or   x = UB^{-1}*y

  The explicitly restarted Lanczos method provides a fast way to find
  the few largest eigenvalues of a generalized eigenvalue problem
  without computing the full spectrum.

  The Lanczos method builds an orthogonal subspace starting with a
  random initial vector. The subspace utilizes m vectors.  A good
  starting point is to set m = 4*k, or some small multiple of the
  number of desired eigenvectors. Large m leads to a larger subspace
  and since we use full orthogonalization, a more expensive iteration,
  but smaller m usually requires more restarts.

  Note that the size of the work space required is: n*(m+2) + 4*m + m*m

  Note that if lm is set to true, we find the largest magnitude
  eigenvalues: those corresponding to both the largest and smallest
  values.

  input:
  lm:     if true return the largest +ve/-ve eigs, otherwise just +ve
  auplo:  is A stored in the upper or lower half = "U" or "L"
  buplo:  is B stored in the upper or lower half = "U" or "L"
  n:      the size of the square matrices A, B
  ka:     the number of sub or super-diagonals in A
  A:      the matrix A stored in a banded diagonal format
  lda:    the leading dimension of A (>= ka + 1)
  kb:     the number of sub or super-diagonas in B
  B:      the matirx B stored in a banded diagonal format
  ldb:    the leading dimension of B (>= kb + 1)
  k:      the number of requested converged eigenvalues
  m:      the size of the Lanczos factorization
  work:   the work array - at least of size n*(m+2) + 4*m + m*m

  output:
  eigs:    the k (or 2*k if (lm == True)) converged eigenvalues
  eigvecs: the k (or 2*k if (lm == True)) converged eigenvectors
*/
int TACSPanelAnalysis::computeEigenvalues(
    int lm, const char *auplo, const char *buplo, int n, int ka, TacsScalar *A,
    int lda, int kb, TacsScalar *B, int ldb, int k, int m, TacsScalar *work,
    double tol, TacsScalar *eigs, TacsScalar *eigvecs) {
  if (k >= m || (lm && 2 * k >= m)) {
    fprintf(stderr,
            "TACSPanelAnalysis: Error, Lanczos method \
requires a larger subspace\n");
    return -1;
  }

  // Factor the B matrix using a Cholesky factorization
  int info;
  LAPACKpbtrf(buplo, &n, &kb, B, &ldb, &info);
  if (info != 0) {
    fprintf(stderr,
            "TACSPanelAnalysis: Error, Cholesky \
factorizaiton failed: %d\n",
            info);
    return info;
  }

  // Set up the data that we will be accessing
  TacsScalar *V = work;                          // (m+1)*n locations
  TacsScalar *tmp = &work[n * (m + 1)];          // n locations
  TacsScalar *diag = &work[n * (m + 2)];         // m locations
  TacsScalar *subdiag = &work[n * (m + 2) + m];  // m locations

  // Keep track of how many eigenvalues have converged
  int nconv = 0, npconv = 0;

  // Randomly generate an initial vector that will be used
  // as the starting point to generate the Lanczos basis
  TacsScalar *v0 = &V[0];
  for (int i = 0; i < n; i++) {
    v0[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Normalize the initial estimate
  int incx = 1;
  double vnrm = BLASnrm2(&n, v0, &incx);
  vnrm = 1.0 / vnrm;
  BLASscal(&n, &vnrm, v0, &incx);

  // Flag to indicate whether we're converged or not
  int converged = 0;

  // Iterate until we've taken quite a number of iterations
  int iter = 0, max_iter = n / m;
  for (; (iter < max_iter) && (!converged); iter++) {
    // Record the last entry of the Lanczos error
    TacsScalar herr = 0.0;

    // Compute the Lanczos factorization using the full MGS. This
    // ensures that the orthogonal basis is less-sensitive to rounding
    // errors, but does increase the cost of the algorithm.
    int start = nconv + npconv;

    for (int i = start; i < m; i++) {
      // Copy over the most recent vector to temporary storage
      memcpy(tmp, &V[n * i], n * sizeof(TacsScalar));

      // Apply the left factor - note that the third string denotes
      // whether or not to use a unit normal
      if (buplo[0] == 'L') {
        BLAStbsv(buplo, "T", "N", &n, &kb, B, &ldb, tmp, &incx);
      } else {
        BLAStbsv(buplo, "N", "N", &n, &kb, B, &ldb, tmp, &incx);
      }

      // Compute V[i+1] <- A*tmp
      TacsScalar a = 1.0, b = 0.0;
      BLASsbmv(auplo, &n, &ka, &a, A, &lda, tmp, &incx, &b, &V[n * (i + 1)],
               &incx);

      // Apply the right factor V[i+1] <- L^{-1}
      if (buplo[0] == 'L') {
        BLAStbsv(buplo, "N", "N", &n, &kb, B, &ldb, &V[n * (i + 1)], &incx);
      } else {
        BLAStbsv(buplo, "T", "N", &n, &kb, B, &ldb, &V[n * (i + 1)], &incx);
      }

      // Orthogonalize the vector V[i+1] against all previous vectors
      for (int j = i; j >= 0; j--) {
        // Compute the dot product h = dot(V[i+1], V[j])
        TacsScalar h = BLASdot(&n, &V[n * (i + 1)], &incx, &V[n * j], &incx);
        if (i == j) {
          diag[i] = h;
        }

        // Remove the component of V[j] from V[i+1]
        h *= -1.0;
        BLASaxpy(&n, &h, &V[n * j], &incx, &V[n * (i + 1)], &incx);
      }

      if (i == m - 1) {
        herr = BLASnrm2(&n, &V[n * (i + 1)], &incx);
      } else {
        // Normalize the vector
        subdiag[i] = BLASnrm2(&n, &V[n * (i + 1)], &incx);
        vnrm = 1.0 / subdiag[i];
        BLASscal(&n, &vnrm, &V[n * (i + 1)], &incx);
      }
    }

    // Set pointers to the required data
    TacsScalar *Z = &work[n * (m + 2) + 2 * m];          // m*m locations
    TacsScalar *W = &work[n * (m + 2) + 2 * m + m * m];  // 2*m locations

    // Compute full eigenvalue spectrum
    int ldz = m, vinfo;
    LAPACKstev("V", &m, diag, subdiag, Z, &ldz, W, &vinfo);
    if (vinfo > 0) {
      printf(
          "DSTEV failed to converge: %d off-diagonal elements \
did not converge to zero\n",
          vinfo);
    }

    // Check convergence for the desired spectrum
    nconv = 0;
    for (int i = 0; (i < m) && (nconv < k); i++) {
      // Apply the convergence criteria
      if (fabs(herr * Z[(i + 1) * m - 1]) < tol) {
        nconv++;
      } else {
        break;
      }
    }

    // Check the convergence of the other half of the spectrum if
    // requested
    npconv = 0;
    if (lm) {
      for (int i = 0; (i < m) && (npconv < k); i++) {
        // Apply the convergence criteria
        if (fabs(herr * Z[(m - i) * m - 1]) < tol) {
          npconv++;
        } else {
          break;
        }
      }
    }

    // Check to see if we've converged or not
    converged = 0;
    if (lm) {
      converged = ((nconv >= k) && (npconv >= k));
    } else {
      converged = (nconv >= k);
    }

    // If we're on the last iteration, just use what we
    // have and return a false flag
    if (iter == max_iter - 1) {
      nconv = k;
      if (lm) {
        npconv = k;
      }
    }

    // If we want the k largest manitude vectors,
    // adjust the range
    if (lm && npconv > 0) {
      for (int i = 0; i < npconv; i++) {
        // Swap the eigenvectors for the projection
        for (int j = 0; j < m; j++) {
          TacsScalar t = Z[(nconv + i) * m + j];
          Z[(nconv + i) * m + j] = Z[(m - npconv + i) * m + j];
          Z[(m - npconv + i) * m + j] = t;
        }

        // Swap the eigenvalues
        TacsScalar t = diag[nconv + i];
        diag[nconv + i] = diag[m - npconv + i];
        diag[m - npconv + i] = t;
      }
    }

    // Check that enough eigenvalues have converged
    if (converged || iter == max_iter - 1) {
      int nev = nconv + npconv;

      for (int i = 0; i < nev; i++) {
        // Copy over the eigenvalues from the diagonal
        eigs[i] = diag[i];
      }

      // Compute the linear combination of the Lanczos vectors
      // that creates the k-desired eigenvectors
      TacsScalar a = 1.0, b = 0.0;
      BLASgemm("N", "N", &n, &nev, &m, &a, V, &n, Z, &ldz, &b, eigvecs, &n);

      // Compute the transform x = LB^{-T}*y   or   x = UB^{-1}*y
      if (buplo[0] == 'L') {
        for (int i = 0; i < nev; i++) {
          BLAStbsv(buplo, "T", "N", &n, &kb, B, &ldb, &eigvecs[n * i], &incx);
        }
      } else {
        for (int i = 0; i < nev; i++) {
          BLAStbsv(buplo, "N", "N", &n, &kb, B, &ldb, &eigvecs[n * i], &incx);
        }
      }
    } else {
      // Zero the sub-diagonal entries. This makes the first ncov
      // entries an identity
      int nev = nconv + npconv;
      for (int i = 0; i < nev; i++) {
        subdiag[i] = 0.0;
      }

      // If the positive eigenvalues have not converged yet, work on
      // them
      if (nconv >= k && npconv < k) {
        memcpy(&Z[m * nev], &Z[m * (m - npconv - 1)], m * sizeof(TacsScalar));
      }

      // Compute the linear combination of vectors
      int nv = nev + 1;
      TacsScalar a = 1.0, b = 0.0;
      BLASgemm("N", "N", &n, &nv, &m, &a, V, &n, Z, &ldz, &b, eigvecs, &n);

      // Copy over the new, frozen, converged eigenvalues
      memcpy(V, eigvecs, n * nv * sizeof(TacsScalar));

      // Pick a new starting Lanczos vector that is a linear
      // combination of the remaining eigenvalues
      TacsScalar *v0 = &V[n * nev];

      // Normalize the Lanczos vector
      vnrm = BLASnrm2(&n, v0, &incx);
      vnrm = 1.0 / vnrm;
      BLASscal(&n, &vnrm, v0, &incx);
    }
  }

  if (!converged) {
    fprintf(stderr, "TACSPanelAnalysis: Error, eigenvalue solve failure\n");
    return -1;
  }

  return 0;
}

/*
  Compute the solution with a uniform pressure load on the skin
  panels.

  Solve the system of equations:

  K*x = f
*/
int TACSPanelAnalysis::computePressureLoad(TacsScalar p,
                                           const char *file_name) {
  double pt[3] = {0.0, 0.0, 0.0};

  int nentries = (nband + 1) * nvars;
  TacsScalar *K = new TacsScalar[nentries];
  TacsScalar *f = new TacsScalar[nvars];

  memset(K, 0, nentries * sizeof(TacsScalar));
  memset(f, 0, nvars * sizeof(TacsScalar));

  for (int k = 0; k < nsegments; k++) {
    // Compute the stiffness matrix and geometric stiffness matrix
    // corresponding to this segment node
    int n1 = nodes[2 * k];
    int n2 = nodes[2 * k + 1];

    // Add the contribution to the stiffness matrix
    TacsScalar At[6], Bt[6], Dt[6], Ast[3];
    panels[k]->getStiffness(pt, At, Bt, Dt, Ast);
    addStiffMat(K, n1, n2, At, Bt, Dt);

    if (segmentType[k] == SKIN_SEGMENT) {
      addPressureLoad(f, n1, n2, p);
    }
  }

  for (int k = 0; k < nbeams; k++) {
    TacsScalar Ct[10];
    beams[k]->getStiffness(pt, Ct);
    addStiffMatBeam(K, bnodes[k], Ct);
  }

  int info, one = 1;
  int ldk = nband + 1;
  LAPACKdpbsv("U", &nvars, &nband, &one, K, &ldk, f, &nvars, &info);

  if (info != 0) {
    fprintf(stderr,
            "TACSPanelAnalysis: Error in LAPACK Cholesky factorization\
 info = %d\n",
            info);
  }

  if (file_name) {
    printPanelMode(file_name, f, 125);
  }

  delete[] K;
  delete[] f;

  return info;
}

/*
  Compute the first 'nloads' critical buckling loads using a linear
  eigenvalue analysis.

  If a prefix is provided, a file is created for each of the nloads
  buckling modes, otherwise no files are produced.

  input:
  Nx:     the axial force per unit length
  Nxy:    the shear force per unit length
  nloads: the number of loads to compute
  prefix: the directory to put the output file

  output:
  loads:  an array of the first nloads critical buckling loads
*/
int TACSPanelAnalysis::computeBucklingLoads(TacsScalar Nx, TacsScalar Nxy,
                                            TacsScalar loads[], int nloads,
                                            const char *prefix) {
  // Allocate space for the eigenvalues and eigenvectors
  TacsScalar *eigvals = new TacsScalar[nloads];
  TacsScalar *eigvecs = new TacsScalar[nloads * nvars];

  // Allocate space for the segment loads
  TacsScalar *segmentLoads = new TacsScalar[3 * nsegments];
  TacsScalar *beamLoads = new TacsScalar[nbeams];

  // Compute the segment loads
  computeSegmentLoads(Nx, Nxy, segmentLoads, beamLoads);

  // Compute the critical buckling eigenvalues
  int info = computeBucklingLoads(segmentLoads, beamLoads, nloads, eigvals,
                                  eigvecs, 0);

  // Compute the actual buckling load, since we apply a transform
  // to the original buckling problem
  for (int k = 0; k < nloads; k++) {
    loads[k] = -1.0 / eigvals[k];
  }

  // Print out the buckling modes
  if (prefix) {
    int file_len = strlen(prefix) + 81;
    char *file_name = new char[file_len];

    for (int i = 0; i < nloads; i++) {
      sprintf(file_name, "%sbuckling_mode%02d.dat", prefix, i);
      printPanelMode(file_name, &eigvecs[nvars * i], 125);
    }

    delete[] file_name;
  }

  // Delete the allocated memory
  delete[] segmentLoads;
  delete[] beamLoads;
  delete[] eigvals;
  delete[] eigvecs;

  return info;
}

/*
  Compute the maximum magnitude positive and negative critical
  buckling modes.

  Positive and negative buckling loads may be different is there is an
  asymmetry in the problem. For instance, if the panel is skewed at an
  angle along the x-directoin.

  input:
  Nxy:    the shear force per unit length
  nloads: the number of required loads

  output:
  posLoad: the positive critical buckling loads
  negLoad: the negative critical buckling loads
*/
int TACSPanelAnalysis::computeBucklingLoads(TacsScalar Nxy,
                                            TacsScalar posLoads[],
                                            TacsScalar negLoads[], int nloads) {
  // Set the axial load per unit length to zero
  TacsScalar Nx = 0.0;

  // Allocate space for the eigenvalues and eigenvectors
  TacsScalar *eigvals = new TacsScalar[2 * nloads];
  TacsScalar *eigvecs = new TacsScalar[2 * nvars * nloads];

  // Allocate space for the segment loads
  TacsScalar *segmentLoads = new TacsScalar[3 * nsegments];
  TacsScalar *beamLoads = new TacsScalar[nbeams];

  // Compute the segment loads in the segments and beams
  computeSegmentLoads(Nx, Nxy, segmentLoads, beamLoads);

  // Compute the critical buckling loads
  int info = computeBucklingLoads(segmentLoads, beamLoads, nloads, eigvals,
                                  eigvecs, (posLoads && negLoads));

  // Calculate the positive and negative critical loads
  // based on the spectral transformation
  for (int k = 0; k < nloads; k++) {
    posLoads[k] = -1.0 / eigvals[k];
    if (negLoads) {
      negLoads[k] = 1.0 / eigvals[2 * nloads - k - 1];
    }
  }

  // Free the allocated space
  delete[] segmentLoads;
  delete[] beamLoads;
  delete[] eigvals;
  delete[] eigvecs;

  return info;
}

/*
  This code computes the eigenvalues associated with a
  buckling problem.

  The code forms the stiffness matrix K, the geometric stiffness
  matrix G and computes the eigenvalues and eigenvectors of the
  generalized eigenvalue problem:

  G*u = eigs[k]*K*u

  Note that we actually want to compute the loads that are the
  solution of this eigenvalue problem:

  K*u + load[k]*G*u = 0

  these two problems are related by: load[k] = -1.0/eigs[k], but this
  transformation is done outside this function.

  input:
  segmentLoads:  the loads in each segment
  beamLoads:     the loads in the beams
  neigs:         the number of eigenvalues to compute

  output:
  eigvals:       an array of length neigs of the eigenvalues
  eigvecs:       an array of length nvars*neigs of the eigenvectors
*/
int TACSPanelAnalysis::computeBucklingLoads(const TacsScalar segmentLoads[],
                                            const TacsScalar beamLoads[],
                                            int neigs, TacsScalar eigvals[],
                                            TacsScalar eigvecs[],
                                            int two_sided) {
  // The parametric point where we evaluate the stiffness
  double pt[3] = {0.0, 0.0, 0.0};

  // K: the stiffness matrix
  // G: the geometric stiffness matrix
  int nentries = (nband + 1) * nvars;

  TacsScalar *K = new TacsScalar[nentries];
  TacsScalar *G = new TacsScalar[nentries];

  memset(K, 0, nentries * sizeof(TacsScalar));
  memset(G, 0, nentries * sizeof(TacsScalar));

  for (int k = 0; k < nsegments; k++) {
    // Compute the stiffness matrix and geometric stiffness matrix
    // corresponding to this segment node
    int n1 = nodes[2 * k];
    int n2 = nodes[2 * k + 1];

    // Add the contribution to the stiffness matrix
    TacsScalar At[6], Bt[6], Dt[6], Ast[3];
    panels[k]->getStiffness(pt, At, Bt, Dt, Ast);
    addStiffMat(K, n1, n2, At, Bt, Dt);

    // Add the contribution to the geometric stiffness matrix
    addGeoStiffMat(G, &segmentLoads[3 * k], n1, n2);
  }

  for (int k = 0; k < nbeams; k++) {
    TacsScalar Ct[10];
    beams[k]->getStiffness(pt, Ct);
    addStiffMatBeam(K, bnodes[k], Ct);
  }

  int info = 0;
  if (use_lapack_eigensolver) {
    int ldz = nvars;
    int ldk = nband + 1, ldg = nband + 1;

    // Allocate space for the matrices required
    TacsScalar *eigs = new TacsScalar[nvars];
    TacsScalar *Z = new TacsScalar[nvars * nvars];
    TacsScalar *work = new TacsScalar[3 * nvars];

    LAPACKdsbgv("V", "U", &nvars, &nband, &nband, G, &ldg, K, &ldk, eigs, Z,
                &ldz, work, &info);

    // Copy over the eigenvalues and eigenvectors that are required
    memcpy(eigvals, eigs, neigs * sizeof(TacsScalar));
    memcpy(eigvecs, Z, neigs * nvars * sizeof(TacsScalar));

    if (two_sided) {
      // Copy the remaining eigenvalues if required - two sided problem
      memcpy(&eigvals[neigs], &eigs[nvars - neigs], neigs * sizeof(TacsScalar));
      memcpy(&eigvecs[neigs * nvars], &Z[(nvars - neigs) * nvars],
             neigs * nvars * sizeof(TacsScalar));
    }

    if (info != 0) {
      fprintf(stderr,
              "TACSPanelAnalysis: Error in LAPACK eigenvalue solver\
 info = %d\n",
              info);
      if (info > nvars) {
        fprintf(stderr,
                "TACSPanelAnalysis: Stiffness matrix is not \
positive definite\n");
      }
    }

    delete[] eigs;
    delete[] Z;
    delete[] work;
  } else {
    // Solve the eigenvalue problem K*u + load*G*u = 0
    int m = lanczos_subspace_size;
    double tol = lanczos_eigen_tol;
    int lwork = nvars * (m + 2) + m * m + 4 * m;
    TacsScalar *work = new TacsScalar[lwork];

    info = computeEigenvalues(two_sided, "U", "U", nvars, nband, G, nband + 1,
                              nband, K, nband + 1, neigs, m, work, tol, eigvals,
                              eigvecs);
    delete[] work;
  }

  delete[] K;
  delete[] G;

  return info;
}

/*
  Compute the derivative of the critical buckling loads with respect
  to the design variables.

  input:
  Nx:     the axial load per unit length
  Nxy:    the shear load per unit length
  nloads: the number of critical buckling loads

  output:
  loads:       the critical buckling load
  loadsDVSens: the derivative of the critical load w.r.t. design vars
*/
/*
int TACSPanelAnalysis::computeBucklingLoadsDVSens( TacsScalar Nx,
                                                   TacsScalar Nxy,
                                                   TacsScalar loads[],
                                                   TacsScalar loadDVSens[],
                                                   int nloads ){
  return computeBucklingLoadsDVSens(Nx, Nxy, loads,
                                    NULL, loadDVSens, NULL, nloads);
}
*/
/*
  Compute the derivative of both the positive and negative critical
  buckling loads with respect to the design variables.

  input:
  Nxy:    the shear load per unit length
  nloads: the number of critical buckling loads

  output:
  loads:       the critical buckling load
  loadsDVSens: the derivative of the critical load w.r.t. design vars
*/
/*
int TACSPanelAnalysis::computeBucklingLoadsDVSens( TacsScalar Nxy,
                                                   TacsScalar posLoads[],
                                                   TacsScalar negLoads[],
                                                   TacsScalar posLoadDVSens[],
                                                   TacsScalar negLoadDVSens[],
                                                   int nloads ){
  return computeBucklingLoadsDVSens(0.0, Nxy, posLoads, negLoads,
                                    posLoadDVSens, negLoadDVSens, nloads);
}
*/
/*
  Compute the derivative of the buckling loads w.r.t. the design
  variables.

  The original eigenvalue problem is:

  G*Z = K*Z*eigs

  The sensitivity of the eigenvalues can be determined as follows:

  det(Z^{T}*(dG/dx - dK/dx*eigs)*Z + d(eigs)/dx) = 0

  where mu is a diagonal matrix of sensitivities.

  This function is private function. The reason is we don't want to
  request both positive and negative loads for Nx != 0.0. This will
  likely lead to an impossible buckling problem.

  input:
  Nx:   the axial load per unit length
  Nxy:  the shear load per unit length

  output:
  posLoads:      the positive critical buckling loads
  negLoads:      the negative critical buckling loads (possibly NULL)
  posLoadDVSens: the derivative of the positive loads w.r.t. dvs
  negLoadDVSens: the derivative of the negative loads w.r.t. dvs (pos. NULL)
*/
/*
int TACSPanelAnalysis::computeBucklingLoadsDVSens( TacsScalar Nx,
                                                   TacsScalar Nxy,
                                                   TacsScalar posLoads[],
                                                   TacsScalar negLoads[],
                                                   TacsScalar posLoadDVSens[],
                                                   TacsScalar negLoadDVSens[],
                                                   int nloads ){
  double pt[3] = {0.0, 0.0, 0.0};

  // Allocate space for the eigenvalues and eigenvectors
  TacsScalar *eigvals = new TacsScalar[ 2*nloads ];
  TacsScalar *eigvecs = new TacsScalar[ 2*nloads*nvars ];

  // Allocate space for the segment loads
  TacsScalar *segmentLoads = new TacsScalar[ 3*nsegments ];
  TacsScalar *beamLoads = new TacsScalar[ nbeams ];

  // Compute the segment loads
  computeSegmentLoads(Nx, Nxy, segmentLoads, beamLoads);

  // Compute the eigenvalues
  int info = computeBucklingLoads(segmentLoads, beamLoads,
                                  nloads, eigvals, eigvecs,
                                  (posLoads && negLoads));

  // Compute the positive and negative critical loads
  for ( int k = 0; k < nloads; k++ ){
    posLoads[k] = -1.0/eigvals[k];

    if (negLoads){
      negLoads[k] = 1.0/eigvals[2*nloads-k-1];
    }
  }

  // For each design variable compute Z^{T}*dK/dx_{i}* Z
  // and compute Z^{T}*dG/dx*Z
  TacsScalar *segmentLoadsSens = new TacsScalar[3*nsegments];
  TacsScalar *beamLoadsSens = new TacsScalar[nbeams];

  TacsScalar *posKSens = new TacsScalar[nloads];
  TacsScalar *posGSens = new TacsScalar[nloads];

  TacsScalar *negKSens = new TacsScalar[nloads];
  TacsScalar *negGSens = new TacsScalar[nloads];

  for ( int n = 0; n < numDesignVars; n++ ){
    int dvNum = designVarNums[n];

    memset(posKSens, 0, nloads*sizeof(TacsScalar));
    memset(posGSens, 0, nloads*sizeof(TacsScalar));

    memset(negKSens, 0, nloads*sizeof(TacsScalar));
    memset(negGSens, 0, nloads*sizeof(TacsScalar));

    if (designVarTypes[n] == 1){
      // Determine the index of the geometric design variable
      int dv = 0;
      for ( int j = 0; j < nDvGeo; j++ ){
        if (geoDvNums[j] == dvNum){
          dv = j;
          break;
        }
      }

      computeSegmentLoadsGeoSens(Nx, Nxy, dv,
                                 segmentLoadsSens, beamLoadsSens);

      for ( int k = 0; k < nsegments; k++ ){
        // Geometric design variable
        TacsScalar At[6], Bt[6], Dt[6], Ast[3];
        panels[k]->getStiffness(pt, At, Bt, Dt, Ast);

        // Add the positive load contribution
        addStiffMatGeoSens(posKSens, dv, eigvecs, nvars, nloads,
                           nodes[2*k], nodes[2*k+1], At, Bt, Dt);
        addGeoStiffMatGeoSens(posGSens, dv, eigvecs, nvars, nloads,
                              nodes[2*k], nodes[2*k+1],
                              &segmentLoads[3*k], &segmentLoadsSens[3*k]);

        // Add the negative load contribution
        if (negLoadDVSens){
          int offset = nvars*nloads;
          addStiffMatGeoSens(negKSens, dv, &eigvecs[offset], nvars, nloads,
                             nodes[2*k], nodes[2*k+1], At, Bt, Dt);
          addGeoStiffMatGeoSens(negGSens, dv, &eigvecs[offset], nvars, nloads,
                                nodes[2*k], nodes[2*k+1],
                                &segmentLoads[3*k], &segmentLoadsSens[3*k]);
        }
      }
    }
    else if (designVarTypes[n] == 2){
      for ( int k = 0; k < nsegments; k++ ){
        TacsScalar At[6], Bt[6], Dt[6], Ast[3];
        panels[k]->getStiffness(pt, At, Bt, Dt, Ast);

        // Add the positive load contribution
        addStiffMatLxSens(posKSens, eigvecs, nvars, nloads,
                          nodes[2*k], nodes[2*k+1], At, Bt, Dt);
        addGeoStiffMatLxSens(posGSens, eigvecs, nvars, nloads,
                             nodes[2*k], nodes[2*k+1],
                             &segmentLoads[3*k]);

          // Add the negative load contribution
        if (negLoadDVSens){
          int offset = nvars*nloads;
          addStiffMatLxSens(negKSens, &eigvecs[offset], nvars, nloads,
                            nodes[2*k], nodes[2*k+1], At, Bt, Dt);
          addGeoStiffMatLxSens(negGSens, &eigvecs[offset], nvars, nloads,
                               nodes[2*k], nodes[2*k+1],
                               &segmentLoads[3*k]);
        }
      }
    }
    else {
      computeSegmentLoadsDVSens(dvNum, Nx, Nxy,
                                segmentLoadsSens, beamLoadsSens);

      for ( int k = 0; k < nsegments; k++ ){
        // Material design variable
        TacsScalar sAt[6], sBt[6], sDt[6], sAst[3];
        panels[k]->getStiffnessDVSens(dvNum, pt, sAt, sBt, sDt, sAst);

        // Add the positive load contribution
        addStiffMatDVSens(posKSens, eigvecs, nvars, nloads,
                          nodes[2*k], nodes[2*k+1], sAt, sBt, sDt);
        addGeoStiffMatDVSens(posGSens, eigvecs, nvars, nloads,
                             nodes[2*k], nodes[2*k+1],
                             &segmentLoadsSens[3*k]);

          // Add the negative load contribution
        if (negLoadDVSens){
          int offset = nvars*nloads;
          addStiffMatDVSens(negKSens, &eigvecs[offset], nvars, nloads,
                            nodes[2*k], nodes[2*k+1], sAt, sBt, sDt);
          addGeoStiffMatDVSens(negGSens, &eigvecs[offset], nvars, nloads,
                               nodes[2*k], nodes[2*k+1],
                               &segmentLoadsSens[3*k]);
        }
      }
    }

    // Store the current sensitivity into the matrix loadDVSens
    for ( int k = 0; k < nloads; k++ ){
      TacsScalar eig = eigvals[k];
      posLoadDVSens[k*numDesignVars + n] =
        (posGSens[k] - eig*posKSens[k])/(eig*eig);

      // Compute the sensitivity of the negative components
      if (negLoadDVSens){
        eig = eigvals[2*nloads-k-1];
        negLoadDVSens[k*numDesignVars + n] =
          -(negGSens[nloads-k-1] - eig*negKSens[nloads-k-1])/(eig*eig);
      }
    }
  }

  delete [] eigvals;
  delete [] eigvecs;
  delete [] segmentLoads;
  delete [] beamLoads;
  delete [] segmentLoadsSens;
  delete [] beamLoadsSens;
  delete [] posKSens;
  delete [] posGSens;
  delete [] negKSens;
  delete [] negGSens;

  return info;
}
*/
/*
  Set up the stiffness matrices and the mass matrix and compute the
  eigenvalues of the modal decomposition problem.
*/
int TACSPanelAnalysis::computeFrequencies(TacsScalar freq[], int nfreq,
                                          const char *prefix) {
  TacsScalar *eigvecs = new TacsScalar[nvars * nfreq];

  int info = computeFrequencies(nfreq, freq, eigvecs);

  for (int k = 0; k < nfreq; k++) {
    freq[k] = sqrt(freq[k]);
  }

  if (prefix) {
    int file_len = strlen(prefix) + 81;
    char *file_name = new char[file_len];

    for (int i = 0; i < nfreq; i++) {
      sprintf(file_name, "%spanel_mode%02d.dat", prefix, i);
      printPanelMode(file_name, &eigvecs[nvars * i], 125);
    }

    delete[] file_name;
  }

  delete[] eigvecs;

  return info;
}

/*
  Form the stiffness matrix K, the mass matrix M and compute the
  eigenvalues and eigenvectors of the generalized eigenvalue problem

  K*x = lambda*M*x

  The frequencies are given by

  freq = sqrt(lambda)
*/
int TACSPanelAnalysis::computeFrequencies(int neigs, TacsScalar eigvals[],
                                          TacsScalar eigvecs[]) {
  double pt[3] = {0.0, 0.0, 0.0};

  // K: the stiffness matrix
  // M: the mass matrix
  int nentries = (nband + 1) * nvars;

  TacsScalar *K = new TacsScalar[nentries];
  TacsScalar *M = new TacsScalar[nentries];

  memset(K, 0, nentries * sizeof(TacsScalar));
  memset(M, 0, nentries * sizeof(TacsScalar));

  for (int k = 0; k < nsegments; k++) {
    // Compute the stiffness matrix and geometric stiffness matrix
    // corresponding to this segment node
    int n1 = nodes[2 * k];
    int n2 = nodes[2 * k + 1];

    // Add the contribution to the stiffness matrix
    TacsScalar At[6], Bt[6], Dt[6], Ast[3];
    panels[k]->getStiffness(pt, At, Bt, Dt, Ast);
    addStiffMat(K, n1, n2, At, Bt, Dt);

    // Add the contribution to the mass matrix
    double pt[2] = {0.0, 0.0};
    TacsScalar mass[2];
    panels[k]->getPointwiseMass(pt, mass);
    addMassMat(M, n1, n2, mass);
  }

  for (int k = 0; k < nbeams; k++) {
    TacsScalar Ct[10];
    beams[k]->getStiffness(pt, Ct);
    addStiffMatBeam(K, bnodes[k], Ct);

    // Add the contribution to the mass matrix
    double pt[2] = {0.0, 0.0};
    TacsScalar mass[6];
    beams[k]->getPointwiseMass(pt, mass);
    addMassMatBeam(M, bnodes[k], mass);
  }

  int info = 0;
  if (use_lapack_eigensolver) {
    // Set the size of the band matrices
    int ldk = nband + 1, ldm = nband + 1;

    // Allocate space for the solution
    TacsScalar *eigs = new TacsScalar[nvars];
    TacsScalar *Z = new TacsScalar[nvars * nvars];
    TacsScalar *work = new TacsScalar[3 * nvars];

    LAPACKdsbgv("V", "U", &nvars, &nband, &nband, K, &ldk, M, &ldm, eigs, Z,
                &nvars, work, &info);

    // Copy over the eigenvalues and eigenvectors that are required
    memcpy(eigvals, eigs, neigs * sizeof(TacsScalar));
    memcpy(eigvecs, Z, neigs * nvars * sizeof(TacsScalar));

    if (info != 0) {
      fprintf(stderr,
              "TACSPanelAnalysis: Error in LAPACK eigenvalue solver\
 info = %d\n",
              info);
      if (info > nvars) {
        fprintf(stderr,
                "TACSPanelAnalysis: Mass matrix not positive definite\n");
      }
    }

    delete[] eigs;
    delete[] work;
    delete[] Z;
  } else {
    // Solve the eigenvalue problem K*u + eig*M*u = 0
    int m = lanczos_subspace_size;
    double tol = lanczos_eigen_tol;
    int lwork = nvars * (m + 2) + m * m + 4 * m;
    TacsScalar *work = new TacsScalar[lwork];
    int two_sided = 0;

    info = computeEigenvalues(two_sided, "U", "U", nvars, nband, K, nband + 1,
                              nband, M, nband + 1, neigs, m, work, tol, eigvals,
                              eigvecs);
    delete[] work;
  }

  delete[] K;
  delete[] M;

  return info;
}

/*
  Compute the derivative of the frequencies w.r.t. the design
  variables.

  The original eigenvalue problem is:

  K*Z = M*Z*eigs

  The sensitivity of the eigenvalues can be determined using the following
  method:

  det(Z^{T}*(dM/dx - dK/dx*eigs)*Z + mu) = 0

  where mu is a diagonal matrix of sensitivities.
*/
/*
int TACSPanelAnalysis::computeFrequenciesDVSens( TacsScalar freq[],
                                                 TacsScalar freqDVSens[],
                                                 int nfreq ){
  double pt[3] = {0.0, 0.0, 0.0};
  TacsScalar *eigvecs = new TacsScalar[ nfreq*nvars ];

  int info = computeFrequencies(nfreq, freq, eigvecs);
  if (info){
    return info;
  }

  for ( int k = 0; k < nfreq; k++ ){
    freq[k] = sqrt(freq[k]);
  }

  // For each design variable compute Z^{T}*dK/dx_{i}*Z
  // and compute Z^{T}*dM/dx*Z
  TacsScalar *KSens = new TacsScalar[nfreq];
  TacsScalar *MSens = new TacsScalar[nfreq];

  for ( int n = 0; n < numDesignVars; n++ ){
    int dvNum = designVarNums[n];

    memset(KSens, 0, nfreq*sizeof(TacsScalar));
    memset(MSens, 0, nfreq*sizeof(TacsScalar));

    if (designVarTypes[n] == 1){
      // Determine the index of the geometric design variable
      int dv = 0;
      for ( int j = 0; j < nDvGeo; j++ ){
        if (geoDvNums[j] == dvNum){
          dv = j;
          break;
        }
      }

      for ( int k = 0; k < nsegments; k++ ){
        // Geometric design variable
        TacsScalar At[6], Bt[6], Dt[6], Ast[3];
        panels[k]->getStiffness(pt, At, Bt, Dt, Ast);

        double pt[2] = {0.0, 0.0};
        TacsScalar mass[2];
        panels[k]->getPointwiseMass(pt, mass);

        addStiffMatGeoSens(KSens, dv, eigvecs, nvars, nfreq,
                           nodes[2*k], nodes[2*k+1], At, Bt, Dt);
        addMassMatGeoSens(MSens, dv, eigvecs, nvars, nfreq,
                          nodes[2*k], nodes[2*k+1], mass);
      }
    }
    else if (designVarTypes[n] == 2){
      for ( int k = 0; k < nsegments; k++ ){
        TacsScalar At[6], Bt[6], Dt[6], Ast[3];
        panels[k]->getStiffness(pt, At, Bt, Dt, Ast);

        double pt[2] = {0.0, 0.0};
        TacsScalar mass[2];
        panels[k]->getPointwiseMass(pt, mass);

        addStiffMatLxSens(KSens, eigvecs, nvars, nfreq,
                          nodes[2*k], nodes[2*k+1], At, Bt, Dt);
        addMassMatLxSens(MSens, eigvecs, nvars, nfreq,
                         nodes[2*k], nodes[2*k+1], mass);
      }
    }
    else {
      for ( int k = 0; k < nsegments; k++ ){
        // Material design variable
        TacsScalar sAt[6], sBt[6], sDt[6], sAst[3];
        panels[k]->getStiffnessDVSens(dvNum, pt, sAt, sBt, sDt, sAst);

        double pt[2] = {0.0, 0.0};
        TacsScalar mass[2];
        panels[k]->getPointwiseMassDVSens(dvNum, pt, mass);

        addStiffMatDVSens(KSens, eigvecs, nvars, nfreq,
                          nodes[2*k], nodes[2*k+1], sAt, sBt, sDt);
        addMassMatDVSens(MSens, eigvecs, nvars, nfreq,
                         nodes[2*k], nodes[2*k+1], mass);
      }
    }

    // Store the current sensitivity into the matrix loadDVSens
    for ( int k = 0; k < nfreq; k++ ){
      freqDVSens[k*numDesignVars + n] =
        0.5*(KSens[k] - freq[k]*freq[k]*MSens[k])/freq[k];
    }
  }

  delete [] eigvecs;
  delete [] KSens;
  delete [] MSens;

  return info;
}
*/

/*
  Compute the contribution to the stiffness matrix from the given
  segment.
*/
void TACSPanelAnalysis::addStiffMat(TacsScalar mat[], int n1, int n2,
                                    const TacsScalar At[],
                                    const TacsScalar Bt[],
                                    const TacsScalar Dt[]) {
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2 * NUM_NODES], Nahp[2 * NUM_NODES], Naahp[2 * NUM_NODES];

  TacsScalar stress[6];
  TacsScalar Bs[24 * NUM_NODES], Bc[24 * NUM_NODES];
  memset(Bs, 0, 24 * NUM_NODES * sizeof(TacsScalar));
  memset(Bc, 0, 24 * NUM_NODES * sizeof(TacsScalar));

  TacsScalar Ks[16 * NUM_NODES * NUM_NODES];

  TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
  TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);

  TacsScalar Le = sqrt(c * c + s * s);
  c = c / Le;
  s = s / Le;

  // Transformation between coordinate frames
  // local vars = t *global vars
  TacsScalar t[4];
  t[0] = c;
  t[1] = s;
  t[2] = -s;
  t[3] = c;

  // The variables for the n-th and m-th modes
  int n_vars[4 * NUM_NODES], m_vars[4 * NUM_NODES];
  for (int m = 1; m <= nmodes; m++) {
    memset(Ks, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

    TacsScalar lambda_m = (M_PI * m) / Lx;

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
      m_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (m - 1) + k];
    }

    // Calculate the stiffness matrix contributions
    for (int k = 0; k < numGauss; k++) {
      // Evaluate the shape functions
      FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

      TacsScalar h = 2.0 / Le;
      TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;

      computeBsin(Bs, lambda_m, h, N, Na, Nhp, Nahp, Naahp);
      computeBcos(Bc, lambda_m, h, N, Na, Nhp, Nahp, Naahp);

      // Compute the sin and cos-coefficient contributions
      for (int i = 0; i < 4 * NUM_NODES; i++) {
        computeStress(stress, &Bs[6 * i], At, Bt, Dt);
        for (int j = 0; j < 4 * NUM_NODES; j++) {
          Ks[4 * NUM_NODES * i + j] +=
              0.5 * detJ * innerProduct(stress, &Bs[6 * j]);
        }

        computeStress(stress, &Bc[6 * i], At, Bt, Dt);
        for (int j = 0; j < 4 * NUM_NODES; j++) {
          Ks[4 * NUM_NODES * i + j] +=
              0.5 * detJ * innerProduct(stress, &Bc[6 * j]);
        }
      }
    }

    // Transform the matrices to the global reference frame
    transformMat(Ks, t, NUM_NODES, NUM_NODES);
    addValues(mat, nband, 4 * NUM_NODES, m_vars, 4 * NUM_NODES, m_vars, Ks);

    for (int n = 1; n <= nmodes; n++) {
      if (n % 2 != m % 2) {
        // Determine the variables for the m-th mode
        for (int k = 0; k < 4; k++) {
          n_vars[k] = vars[(4 * nmodes) * n1 + 4 * (n - 1) + k];
          n_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (n - 1) + k];
        }

        TacsScalar lambda_n = (M_PI * n) / Lx;
        TacsScalar sc = (2.0 * n) / (M_PI * (n * n - m * m));

        memset(Ks, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

        // Compute sin(n) * cos(m) term
        for (int k = 0; k < numGauss; k++) {
          // Evaluate the shape functions
          FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
          FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

          TacsScalar h = 2.0 / Le;
          TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;

          computeBsin(Bs, lambda_n, h, N, Na, Nhp, Nahp, Naahp);
          computeBcos(Bc, lambda_m, h, N, Na, Nhp, Nahp, Naahp);

          // Compute the sin-coefficient contribution
          for (int i = 0; i < 4 * NUM_NODES; i++) {
            computeStress(stress, &Bs[6 * i], At, Bt, Dt);
            for (int j = 0; j < 4 * NUM_NODES; j++) {
              Ks[4 * NUM_NODES * i + j] +=
                  sc * detJ * innerProduct(stress, &Bc[6 * j]);
            }
          }
        }

        // Transform the result to the global reference frame
        transformMat(Ks, t, NUM_NODES, NUM_NODES);

        // Add values to the matrix
        addValues(mat, nband, 4 * NUM_NODES, n_vars, 4 * NUM_NODES, m_vars, Ks);
        addValuesTranspose(mat, nband, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                           n_vars, Ks);
      }
    }
  }
}

/*
  Compute the contribution to the stiffness matrix from the given
  segment.
*/
void TACSPanelAnalysis::addPressureLoad(TacsScalar F[], int n1, int n2,
                                        TacsScalar p) {
  double Nhp[2 * NUM_NODES], Nahp[2 * NUM_NODES], Naahp[2 * NUM_NODES];
  TacsScalar Fs[4 * NUM_NODES];

  TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
  TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);

  TacsScalar Le = sqrt(c * c + s * s);
  c = c / Le;
  s = s / Le;

  // Transformation between coordinate frames
  // local vars = t * global vars
  TacsScalar t[4];
  t[0] = c;
  t[1] = s;
  t[2] = -s;
  t[3] = c;

  // The variables for the n-th and m-th modes
  int m_vars[4 * NUM_NODES];
  for (int m = 1; m <= nmodes; m += 2) {
    memset(Fs, 0, 4 * NUM_NODES * sizeof(TacsScalar));

    TacsScalar scale = (2.0 * p) / (M_PI * m);

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
      m_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (m - 1) + k];
    }

    // Calculate the stiffness matrix contributions
    for (int k = 0; k < numGauss; k++) {
      // Evaluate the shape functions
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

      TacsScalar h = 2.0 / Le;
      TacsScalar hinv = 1.0 / h;
      TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;

      // Compute the sin and cos-coefficient contributions
      for (int i = 0; i < 2; i++) {
        Fs[4 * i + 2] += scale * detJ * Nhp[2 * i];
        Fs[4 * i + 3] += scale * detJ * hinv * Nhp[2 * i + 1];
      }
    }

    // Transform the vector
    transformVec(Fs, t, NUM_NODES);

    for (int i = 0; i < 4 * NUM_NODES; i++) {
      if (m_vars[i] >= 0) {
        F[m_vars[i]] += Fs[i];
      }
    }
  }
}

/*
  Add the contributions to the stiffness matrix from a longitudinal
  beam
*/
void TACSPanelAnalysis::addStiffMatBeam(TacsScalar mat[], int n1,
                                        const TacsScalar Ct[]) {
  TacsScalar stress[4];
  TacsScalar Bs[16], Bc[16];

  memset(Bs, 0, 16 * sizeof(TacsScalar));
  memset(Bc, 0, 16 * sizeof(TacsScalar));

  TacsScalar Ks[16];

  // The variables for the n-th and m-th modes
  int n_vars[4], m_vars[4];

  for (int m = 1; m <= nmodes; m++) {
    memset(Ks, 0, 16 * sizeof(TacsScalar));

    TacsScalar lambda_m = (M_PI * m) / Lx;

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
    }

    computeBeamBsin(Bs, lambda_m);
    computeBeamBcos(Bc, lambda_m);

    // Compute the sin and cos-coefficient contributions
    for (int i = 0; i < 4; i++) {
      computeBeamStress(stress, &Bs[4 * i], Ct);
      for (int j = 0; j < 4; j++) {
        Ks[4 * i + j] += 0.5 * Lx * innerBeamProduct(stress, &Bs[4 * j]);
      }

      computeBeamStress(stress, &Bc[4 * i], Ct);
      for (int j = 0; j < 4; j++) {
        Ks[4 * i + j] += 0.5 * Lx * innerBeamProduct(stress, &Bc[4 * j]);
      }
    }

    addValues(mat, nband, 4, m_vars, 4, m_vars, Ks);

    for (int n = 1; n <= nmodes; n++) {
      if (n % 2 != m % 2) {
        // Determine the variables for the m-th mode
        for (int k = 0; k < 4; k++) {
          n_vars[k] = vars[(4 * nmodes) * n1 + 4 * (n - 1) + k];
        }

        TacsScalar lambda_n = (M_PI * n) / Lx;
        TacsScalar sc = (2.0 * n) / (M_PI * (n * n - m * m));

        memset(Ks, 0, 16 * sizeof(TacsScalar));

        computeBeamBsin(Bs, lambda_n);
        computeBeamBcos(Bc, lambda_m);

        // Compute the sin-coefficient contribution
        for (int i = 0; i < 4; i++) {
          computeBeamStress(stress, &Bs[4 * i], Ct);
          for (int j = 0; j < 4; j++) {
            Ks[4 * i + j] += sc * Lx * innerBeamProduct(stress, &Bc[4 * j]);
          }
        }

        // Add values to the matrix
        addValues(mat, nband, 4 * NUM_NODES, n_vars, 4 * NUM_NODES, m_vars, Ks);
        addValuesTranspose(mat, nband, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                           n_vars, Ks);
      }
    }
  }
}

/*
  Compute the contribution to the mass matrix from the given segment
*/
void TACSPanelAnalysis::addMassMat(TacsScalar mat[], int n1, int n2,
                                   const TacsScalar mass[]) {
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2 * NUM_NODES], Nahp[2 * NUM_NODES], Naahp[2 * NUM_NODES];

  TacsScalar ks[5 * 4 * NUM_NODES], kc[5 * 4 * NUM_NODES];
  memset(ks, 0, 20 * NUM_NODES * sizeof(TacsScalar));
  memset(kc, 0, 20 * NUM_NODES * sizeof(TacsScalar));

  TacsScalar Ms[16 * NUM_NODES * NUM_NODES];

  TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
  TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);

  TacsScalar Le = sqrt(c * c + s * s);
  c = c / Le;
  s = s / Le;

  // Transformation between coordinate frames
  // local vars = t * global vars
  TacsScalar t[4];
  t[0] = c;
  t[1] = s;
  t[2] = -s;
  t[3] = c;

  // The variables for the n-th and m-th modes
  int n_vars[4 * NUM_NODES], m_vars[4 * NUM_NODES];

  for (int m = 1; m <= nmodes; m++) {
    memset(Ms, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

    TacsScalar lambda_m = (M_PI * m) / Lx;

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
      m_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (m - 1) + k];
    }

    // Calculate the stiffness matrix contributions
    for (int k = 0; k < numGauss; k++) {
      // Evaluate the shape functions
      FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

      TacsScalar h = 2.0 / Le;
      TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;

      computeMsin(ks, lambda_m, h, N, Nhp, Nahp);
      computeMcos(kc, lambda_m, h, N, Nhp, Nahp);

      // Compute the sin and cos-coefficient contributions
      for (int i = 0; i < 4 * NUM_NODES; i++) {
        for (int j = 0; j < 4 * NUM_NODES; j++) {
          Ms[4 * NUM_NODES * i + j] +=
              0.5 * detJ * innerMassProduct(mass, &ks[5 * i], &ks[5 * j]);
        }

        for (int j = 0; j < 4 * NUM_NODES; j++) {
          Ms[4 * NUM_NODES * i + j] +=
              0.5 * detJ * innerMassProduct(mass, &kc[5 * i], &kc[5 * j]);
        }
      }
    }

    // Transform the matrices to the global reference frame
    transformMat(Ms, t, NUM_NODES, NUM_NODES);
    addValues(mat, nband, 4 * NUM_NODES, m_vars, 4 * NUM_NODES, m_vars, Ms);

    for (int n = 1; n <= nmodes; n++) {
      if (n % 2 != m % 2) {
        // Determine the variables for the m-th mode
        for (int k = 0; k < 4; k++) {
          n_vars[k] = vars[(4 * nmodes) * n1 + 4 * (n - 1) + k];
          n_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (n - 1) + k];
        }

        TacsScalar lambda_n = (M_PI * n) / Lx;
        TacsScalar sc = (2.0 * n) / (M_PI * (n * n - m * m));

        memset(Ms, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

        // Compute sin(n) * cos(m) term
        for (int k = 0; k < numGauss; k++) {
          // Evaluate the shape functions
          FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
          FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

          TacsScalar h = 2.0 / Le;
          TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;

          computeMsin(ks, lambda_n, h, N, Nhp, Nahp);
          computeMcos(kc, lambda_m, h, N, Nhp, Nahp);

          // Compute the sine-coefficient contribution
          for (int i = 0; i < 4 * NUM_NODES; i++) {
            for (int j = 0; j < 4 * NUM_NODES; j++) {
              Ms[4 * NUM_NODES * i + j] +=
                  sc * detJ * innerMassProduct(mass, &ks[5 * i], &kc[5 * j]);
            }
          }
        }

        // Transform the result to the global reference frame
        transformMat(Ms, t, NUM_NODES, NUM_NODES);

        // Add values to the matrix
        addValues(mat, nband, 4 * NUM_NODES, n_vars, 4 * NUM_NODES, m_vars, Ms);
        addValuesTranspose(mat, nband, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                           n_vars, Ms);
      }
    }
  }
}

/*
  Add the contributions to the stiffness matrix from a longitudinal
  beam
*/
void TACSPanelAnalysis::addMassMatBeam(TacsScalar mat[], int n1,
                                       const TacsScalar mass[]) {
  TacsScalar ks[20], kc[20];
  memset(ks, 0, 20 * sizeof(TacsScalar));
  memset(kc, 0, 20 * sizeof(TacsScalar));

  TacsScalar Ms[16];

  // The variables for the n-th and m-th modes
  int n_vars[4], m_vars[4];

  for (int m = 1; m <= nmodes; m++) {
    memset(Ms, 0, 16 * sizeof(TacsScalar));

    TacsScalar lambda_m = (M_PI * m) / Lx;

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
    }

    computeBeamMsin(ks, lambda_m);
    computeBeamMcos(kc, lambda_m);

    // Compute the sin and cos-coefficient contributions
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        Ms[4 * i + j] +=
            0.5 * Lx * innerBeamMassProduct(mass, &ks[5 * i], &ks[5 * j]);
      }

      for (int j = 0; j < 4; j++) {
        Ms[4 * i + j] +=
            0.5 * Lx * innerBeamMassProduct(mass, &kc[5 * i], &kc[5 * j]);
      }
    }

    addValues(mat, nband, 4, m_vars, 4, m_vars, Ms);

    for (int n = 1; n <= nmodes; n++) {
      if (n % 2 != m % 2) {
        // Determine the variables for the m-th mode
        for (int k = 0; k < 4; k++) {
          n_vars[k] = vars[(4 * nmodes) * n1 + 4 * (n - 1) + k];
        }

        TacsScalar lambda_n = (M_PI * n) / Lx;
        TacsScalar sc = (2.0 * n) / (M_PI * (n * n - m * m));

        memset(Ms, 0, 16 * sizeof(TacsScalar));

        computeBeamMsin(ks, lambda_n);
        computeBeamMcos(kc, lambda_m);

        // Compute the sin-coefficient contribution
        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            Ms[4 * i + j] +=
                sc * Lx * innerBeamMassProduct(mass, &ks[5 * i], &kc[5 * j]);
          }
        }

        // Add values to the matrix
        addValues(mat, nband, 4 * NUM_NODES, n_vars, 4 * NUM_NODES, m_vars, Ms);
        addValuesTranspose(mat, nband, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                           n_vars, Ms);
      }
    }
  }
}

/*
  Add a segment to the geometric stiffness matrix
*/
void TACSPanelAnalysis::addGeoStiffMat(TacsScalar mat[],
                                       const TacsScalar stress[], int n1,
                                       int n2) {
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2 * NUM_NODES], Nahp[2 * NUM_NODES], Naahp[2 * NUM_NODES];

  TacsScalar Gs[16 * NUM_NODES * NUM_NODES];

  TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
  TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);

  TacsScalar Le = sqrt(c * c + s * s);
  c = c / Le;
  s = s / Le;

  // Transformation between coordinate frames
  // local vars = t * global vars
  TacsScalar t[4];
  t[0] = c;
  t[1] = s;
  t[2] = -s;
  t[3] = c;

  // The variables for the n-th and m-th modes
  int n_vars[4 * NUM_NODES], m_vars[4 * NUM_NODES];

  for (int m = 1; m <= nmodes; m++) {
    memset(Gs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

    TacsScalar lambda_m = (M_PI * m) / Lx;

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
      m_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (m - 1) + k];
    }

    // Calculate the stiffness matrix contributions
    for (int k = 0; k < numGauss; k++) {
      // Evaluate the shape functions
      FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

      TacsScalar h = 2.0 / Le;
      TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;

      addGs(Gs, 0.5 * detJ, stress, lambda_m, h, N, Na, Nhp, Nahp);
    }

    transformMat(Gs, t, NUM_NODES, NUM_NODES);
    addValues(mat, nband, 4 * NUM_NODES, m_vars, 4 * NUM_NODES, m_vars, Gs);

    for (int n = 1; n <= nmodes; n++) {
      if (n % 2 != m % 2) {
        // Determine the variables for the n-th mode
        for (int k = 0; k < 4; k++) {
          n_vars[k] = vars[(4 * nmodes) * n1 + 4 * (n - 1) + k];
          n_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (n - 1) + k];
        }

        TacsScalar lambda_n = (M_PI * n) / Lx;
        TacsScalar sc = (2.0 * n) / (M_PI * (n * n - m * m));

        memset(Gs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

        for (int k = 0; k < numGauss; k++) {
          // Evaluate the shape functions
          FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
          FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

          TacsScalar h = 2.0 / Le;
          TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;

          addGcs(Gs, sc * detJ, stress, lambda_n, lambda_m, h, N, Na, Nhp,
                 Nahp);
        }

        // Transform the result to the global reference frame
        transformMat(Gs, t, NUM_NODES, NUM_NODES);

        // Add values to the matrix
        addValues(mat, nband, 4 * NUM_NODES, n_vars, 4 * NUM_NODES, m_vars, Gs);
        addValuesTranspose(mat, nband, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                           n_vars, Gs);
      }
    }
  }
}

/*
  Compute the product:

  dK = Z^{T} * dK/dx * Z

  where Z are a set of eigenvectors. The matrix dK can then be used
  within a sensitivity analysis problem. Note that the matrix dK is
  computed using a matrix-free method.
*/
void TACSPanelAnalysis::addStiffMatDVSens(
    TacsScalar smat[], const TacsScalar Z[], int ldz, int nz, int n1, int n2,
    const TacsScalar sAt[], const TacsScalar sBt[], const TacsScalar sDt[]) {
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2 * NUM_NODES], Nahp[2 * NUM_NODES], Naahp[2 * NUM_NODES];

  TacsScalar stress[6];
  TacsScalar Bs[24 * NUM_NODES], Bc[24 * NUM_NODES];
  memset(Bs, 0, 24 * NUM_NODES * sizeof(TacsScalar));
  memset(Bc, 0, 24 * NUM_NODES * sizeof(TacsScalar));

  TacsScalar Ks[16 * NUM_NODES * NUM_NODES];

  TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
  TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);

  TacsScalar Le = sqrt(c * c + s * s);
  c = c / Le;
  s = s / Le;

  // Transformation between coordinate frames
  // local vars = t * global vars
  TacsScalar t[4];
  t[0] = c;
  t[1] = s;
  t[2] = -s;
  t[3] = c;

  // The variables for the n-th and m-th modes
  int n_vars[4 * NUM_NODES], m_vars[4 * NUM_NODES];

  for (int m = 1; m <= nmodes; m++) {
    memset(Ks, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

    TacsScalar lambda_m = (M_PI * m) / Lx;

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
      m_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (m - 1) + k];
    }

    // Calculate the stiffness matrix contributions
    for (int k = 0; k < numGauss; k++) {
      // Evaluate the shape functions
      FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

      TacsScalar h = 2.0 / Le;
      TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;

      computeBsin(Bs, lambda_m, h, N, Na, Nhp, Nahp, Naahp);
      computeBcos(Bc, lambda_m, h, N, Na, Nhp, Nahp, Naahp);

      // Compute the sin and cos-coefficient contributions
      for (int i = 0; i < 4 * NUM_NODES; i++) {
        computeStress(stress, &Bs[6 * i], sAt, sBt, sDt);
        for (int j = 0; j < 4 * NUM_NODES; j++) {
          Ks[4 * NUM_NODES * i + j] +=
              0.5 * detJ * innerProduct(stress, &Bs[6 * j]);
        }

        computeStress(stress, &Bc[6 * i], sAt, sBt, sDt);
        for (int j = 0; j < 4 * NUM_NODES; j++) {
          Ks[4 * NUM_NODES * i + j] +=
              0.5 * detJ * innerProduct(stress, &Bc[6 * j]);
        }
      }
    }

    // Transform the matrices to the global reference frame
    transformMat(Ks, t, NUM_NODES, NUM_NODES);
    addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                     m_vars, Ks);

    for (int n = 1; n <= nmodes; n++) {
      if (n % 2 != m % 2) {
        // Determine the variables for the m-th mode
        for (int k = 0; k < 4; k++) {
          n_vars[k] = vars[(4 * nmodes) * n1 + 4 * (n - 1) + k];
          n_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (n - 1) + k];
        }

        TacsScalar lambda_n = (M_PI * n) / Lx;
        TacsScalar sc = (2.0 * n) / (M_PI * (n * n - m * m));

        memset(Ks, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

        // Compute sin(n) * cos(m) term
        for (int k = 0; k < numGauss; k++) {
          // Evaluate the shape functions
          FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
          FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

          TacsScalar h = 2.0 / Le;
          TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;

          computeBsin(Bs, lambda_n, h, N, Na, Nhp, Nahp, Naahp);
          computeBcos(Bc, lambda_m, h, N, Na, Nhp, Nahp, Naahp);

          // Compute the sin-coefficient contribution
          for (int i = 0; i < 4 * NUM_NODES; i++) {
            computeStress(stress, &Bs[6 * i], sAt, sBt, sDt);
            for (int j = 0; j < 4 * NUM_NODES; j++) {
              Ks[4 * NUM_NODES * i + j] +=
                  sc * detJ * innerProduct(stress, &Bc[6 * j]);
            }
          }
        }

        // Transform the result to the global reference frame
        transformMat(Ks, t, NUM_NODES, NUM_NODES);

        // Add values to the matrix
        addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, n_vars, 4 * NUM_NODES,
                         m_vars, Ks);
        addSensMatValuesTranspose(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars,
                                  4 * NUM_NODES, n_vars, Ks);
      }
    }
  }
}

/*
  Add the sensitivity contribution from the mass matrix
*/
void TACSPanelAnalysis::addMassMatDVSens(TacsScalar smat[],
                                         const TacsScalar Z[], int ldz, int nz,
                                         int n1, int n2,
                                         const TacsScalar smass[]) {
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2 * NUM_NODES], Nahp[2 * NUM_NODES], Naahp[2 * NUM_NODES];

  TacsScalar ks[5 * 4 * NUM_NODES], kc[5 * 4 * NUM_NODES];
  memset(ks, 0, 20 * NUM_NODES * sizeof(TacsScalar));
  memset(kc, 0, 20 * NUM_NODES * sizeof(TacsScalar));

  TacsScalar Ms[16 * NUM_NODES * NUM_NODES];

  TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
  TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);

  TacsScalar Le = sqrt(c * c + s * s);
  c = c / Le;
  s = s / Le;

  // Transformation between coordinate frames
  // local vars = t * global vars
  TacsScalar t[4];
  t[0] = c;
  t[1] = s;
  t[2] = -s;
  t[3] = c;

  // The variables for the n-th and m-th modes
  int n_vars[4 * NUM_NODES], m_vars[4 * NUM_NODES];

  for (int m = 1; m <= nmodes; m++) {
    memset(Ms, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

    TacsScalar lambda_m = (M_PI * m) / Lx;

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
      m_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (m - 1) + k];
    }

    // Calculate the stiffness matrix contributions
    for (int k = 0; k < numGauss; k++) {
      // Evaluate the shape functions
      FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

      TacsScalar h = 2.0 / Le;
      TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;

      computeMsin(ks, lambda_m, h, N, Nhp, Nahp);
      computeMcos(kc, lambda_m, h, N, Nhp, Nahp);

      // Compute the sin and cos-coefficient contributions
      for (int i = 0; i < 4 * NUM_NODES; i++) {
        for (int j = 0; j < 4 * NUM_NODES; j++) {
          Ms[4 * NUM_NODES * i + j] +=
              0.5 * detJ * innerMassProduct(smass, &ks[5 * i], &ks[5 * j]);
        }

        for (int j = 0; j < 4 * NUM_NODES; j++) {
          Ms[4 * NUM_NODES * i + j] +=
              0.5 * detJ * innerMassProduct(smass, &kc[5 * i], &kc[5 * j]);
        }
      }
    }

    // Transform the matrices to the global reference frame
    transformMat(Ms, t, NUM_NODES, NUM_NODES);
    addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                     m_vars, Ms);

    for (int n = 1; n <= nmodes; n++) {
      if (n % 2 != m % 2) {
        // Determine the variables for the m-th mode
        for (int k = 0; k < 4; k++) {
          n_vars[k] = vars[(4 * nmodes) * n1 + 4 * (n - 1) + k];
          n_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (n - 1) + k];
        }

        TacsScalar lambda_n = (M_PI * n) / Lx;
        TacsScalar sc = (2.0 * n) / (M_PI * (n * n - m * m));

        memset(Ms, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

        // Compute sin(n) * cos(m) term
        for (int k = 0; k < numGauss; k++) {
          // Evaluate the shape functions
          FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
          FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

          TacsScalar h = 2.0 / Le;
          TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;

          computeMsin(ks, lambda_n, h, N, Nhp, Nahp);
          computeMcos(kc, lambda_m, h, N, Nhp, Nahp);

          // Compute the sine-coefficient contribution
          for (int i = 0; i < 4 * NUM_NODES; i++) {
            for (int j = 0; j < 4 * NUM_NODES; j++) {
              Ms[4 * NUM_NODES * i + j] +=
                  sc * detJ * innerMassProduct(smass, &ks[5 * i], &kc[5 * j]);
            }
          }
        }

        // Transform the result to the global reference frame
        transformMat(Ms, t, NUM_NODES, NUM_NODES);

        // Add values to the matrix
        addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                         n_vars, Ms);
        addSensMatValuesTranspose(smat, Z, ldz, nz, 4 * NUM_NODES, n_vars,
                                  4 * NUM_NODES, m_vars, Ms);
      }
    }
  }
}

/*
  Add the contribution to the sensitivity from a segment to the
  geometric stiffness matrix
*/
void TACSPanelAnalysis::addGeoStiffMatDVSens(TacsScalar smat[],
                                             const TacsScalar Z[], int ldz,
                                             int nz, int n1, int n2,
                                             const TacsScalar sstress[]) {
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2 * NUM_NODES], Nahp[2 * NUM_NODES], Naahp[2 * NUM_NODES];

  TacsScalar Gs[16 * NUM_NODES * NUM_NODES];

  TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
  TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);

  TacsScalar Le = sqrt(c * c + s * s);
  c = c / Le;
  s = s / Le;

  // Transformation between coordinate frames
  // local vars = t * global vars
  TacsScalar t[4];
  t[0] = c;
  t[1] = s;
  t[2] = -s;
  t[3] = c;

  // The variables for the n-th and m-th modes
  int n_vars[4 * NUM_NODES], m_vars[4 * NUM_NODES];

  for (int m = 1; m <= nmodes; m++) {
    memset(Gs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

    TacsScalar lambda_m = (M_PI * m) / Lx;

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
      m_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (m - 1) + k];
    }

    // Calculate the stiffness matrix contributions
    for (int k = 0; k < numGauss; k++) {
      // Evaluate the shape functions
      FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

      TacsScalar h = 2.0 / Le;
      TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;

      addGs(Gs, 0.5 * detJ, sstress, lambda_m, h, N, Na, Nhp, Nahp);
    }

    transformMat(Gs, t, NUM_NODES, NUM_NODES);
    addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                     m_vars, Gs);

    for (int n = 1; n <= nmodes; n++) {
      if (n % 2 != m % 2) {
        // Determine the variables for the n-th mode
        for (int k = 0; k < 4; k++) {
          n_vars[k] = vars[(4 * nmodes) * n1 + 4 * (n - 1) + k];
          n_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (n - 1) + k];
        }

        TacsScalar lambda_n = (M_PI * n) / Lx;
        TacsScalar sc = (2.0 * n) / (M_PI * (n * n - m * m));

        memset(Gs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

        for (int k = 0; k < numGauss; k++) {
          // Evaluate the shape functions
          FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
          FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

          TacsScalar h = 2.0 / Le;
          TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;

          addGcs(Gs, sc * detJ, sstress, lambda_n, lambda_m, h, N, Na, Nhp,
                 Nahp);
        }

        // Transform the result to the global reference frame
        transformMat(Gs, t, NUM_NODES, NUM_NODES);

        // Add values to the matrix
        addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, n_vars, 4 * NUM_NODES,
                         m_vars, Gs);
        addSensMatValuesTranspose(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars,
                                  4 * NUM_NODES, n_vars, Gs);
      }
    }
  }
}

/*
  Add the geometric design variable sensitivity contributions to the
  sensitivity matrix.
*/
void TACSPanelAnalysis::addStiffMatGeoSens(TacsScalar smat[], int dv,
                                           const TacsScalar Z[], int ldz,
                                           int nz, int n1, int n2,
                                           const TacsScalar At[],
                                           const TacsScalar Bt[],
                                           const TacsScalar Dt[]) {
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2 * NUM_NODES], Nahp[2 * NUM_NODES], Naahp[2 * NUM_NODES];

  TacsScalar stress[6], sstress[6];
  TacsScalar Bs[24 * NUM_NODES], Bc[24 * NUM_NODES];
  TacsScalar sBs[24 * NUM_NODES], sBc[24 * NUM_NODES];

  memset(Bs, 0, 24 * NUM_NODES * sizeof(TacsScalar));
  memset(Bc, 0, 24 * NUM_NODES * sizeof(TacsScalar));
  memset(sBs, 0, 24 * NUM_NODES * sizeof(TacsScalar));
  memset(sBc, 0, 24 * NUM_NODES * sizeof(TacsScalar));

  TacsScalar Ks[16 * NUM_NODES * NUM_NODES];
  TacsScalar sKs[16 * NUM_NODES * NUM_NODES];

  TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
  TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
  TacsScalar sc =
      (XptLin[2 * nnodes * dv + 2 * n2] - XptLin[2 * nnodes * dv + 2 * n1]);
  TacsScalar ss = (XptLin[2 * nnodes * dv + 2 * n2 + 1] -
                   XptLin[2 * nnodes * dv + 2 * n1 + 1]);

  TacsScalar Le = sqrt(c * c + s * s);
  TacsScalar sLe = (c * sc + s * ss) / Le;
  sc = (sc * Le - sLe * c) / (Le * Le);
  ss = (ss * Le - sLe * s) / (Le * Le);
  c = c / Le;
  s = s / Le;

  // Transformation between coordinate frames
  // local vars = t * global vars
  TacsScalar t[4], st[4];
  t[0] = c;
  t[1] = s;
  t[2] = -s;
  t[3] = c;
  st[0] = sc;
  st[1] = ss;
  st[2] = -ss;
  st[3] = sc;

  // The variables for the n-th and m-th modes
  int n_vars[4 * NUM_NODES], m_vars[4 * NUM_NODES];

  for (int m = 1; m <= nmodes; m++) {
    memset(Ks, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));
    memset(sKs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

    TacsScalar lambda_m = (M_PI * m) / Lx;

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
      m_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (m - 1) + k];
    }

    // Calculate the stiffness matrix contributions
    for (int k = 0; k < numGauss; k++) {
      // Evaluate the shape functions
      FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

      TacsScalar h = 2.0 / Le;
      TacsScalar sh = -2.0 * sLe / (Le * Le);
      TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;
      TacsScalar sdetJ = 0.5 * gaussWts[k] * sLe * Lx;

      computeBsin(Bs, lambda_m, h, N, Na, Nhp, Nahp, Naahp);
      computeBcos(Bc, lambda_m, h, N, Na, Nhp, Nahp, Naahp);

      computeBsinGeoSens(sBs, lambda_m, 0.0, h, sh, N, Na, Nhp, Nahp, Naahp);
      computeBcosGeoSens(sBc, lambda_m, 0.0, h, sh, N, Na, Nhp, Nahp, Naahp);

      // Compute the sin and cos-coefficient contributions
      for (int i = 0; i < 4 * NUM_NODES; i++) {
        computeStress(stress, &Bs[6 * i], At, Bt, Dt);
        computeStress(sstress, &sBs[6 * i], At, Bt, Dt);
        for (int j = 0; j < 4 * NUM_NODES; j++) {
          TacsScalar p = innerProduct(stress, &Bs[6 * j]);
          Ks[4 * NUM_NODES * i + j] += 0.5 * detJ * p;
          sKs[4 * NUM_NODES * i + j] +=
              0.5 * sdetJ * p + 0.5 * detJ *
                                    (innerProduct(sstress, &Bs[6 * j]) +
                                     innerProduct(stress, &sBs[6 * j]));
        }

        computeStress(stress, &Bc[6 * i], At, Bt, Dt);
        computeStress(sstress, &sBc[6 * i], At, Bt, Dt);
        for (int j = 0; j < 4 * NUM_NODES; j++) {
          TacsScalar p = innerProduct(stress, &Bc[6 * j]);
          Ks[4 * NUM_NODES * i + j] += 0.5 * detJ * p;
          sKs[4 * NUM_NODES * i + j] +=
              0.5 * sdetJ * p + 0.5 * detJ *
                                    (innerProduct(sstress, &Bc[6 * j]) +
                                     innerProduct(stress, &sBc[6 * j]));
        }
      }
    }

    // Transform the matrices to the global reference frame
    transformMatGeoSens(Ks, sKs, t, st, NUM_NODES, NUM_NODES);
    addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                     m_vars, sKs);

    for (int n = 1; n <= nmodes; n++) {
      if (n % 2 != m % 2) {
        // Determine the variables for the m-th mode
        for (int k = 0; k < 4; k++) {
          n_vars[k] = vars[(4 * nmodes) * n1 + 4 * (n - 1) + k];
          n_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (n - 1) + k];
        }

        TacsScalar lambda_n = (M_PI * n) / Lx;
        TacsScalar sc = (2.0 * n) / (M_PI * (n * n - m * m));

        memset(Ks, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));
        memset(sKs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

        // Compute sin(n) * cos(m) term
        for (int k = 0; k < numGauss; k++) {
          // Evaluate the shape functions
          FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
          FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

          TacsScalar h = 2.0 / Le;
          TacsScalar sh = -2.0 * sLe / (Le * Le);
          TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;
          TacsScalar sdetJ = 0.5 * gaussWts[k] * sLe * Lx;

          computeBsin(Bs, lambda_n, h, N, Na, Nhp, Nahp, Naahp);
          computeBcos(Bc, lambda_m, h, N, Na, Nhp, Nahp, Naahp);

          computeBsinGeoSens(sBs, lambda_n, 0.0, h, sh, N, Na, Nhp, Nahp,
                             Naahp);
          computeBcosGeoSens(sBc, lambda_m, 0.0, h, sh, N, Na, Nhp, Nahp,
                             Naahp);

          // Compute the sin-coefficient contribution
          for (int i = 0; i < 4 * NUM_NODES; i++) {
            computeStress(stress, &Bs[6 * i], At, Bt, Dt);
            computeStress(sstress, &sBs[6 * i], At, Bt, Dt);
            for (int j = 0; j < 4 * NUM_NODES; j++) {
              TacsScalar p = innerProduct(stress, &Bc[6 * j]);
              Ks[4 * NUM_NODES * i + j] += sc * detJ * p;
              sKs[4 * NUM_NODES * i + j] +=
                  sc * sdetJ * p + sc * detJ *
                                       (innerProduct(sstress, &Bc[6 * j]) +
                                        innerProduct(stress, &sBc[6 * j]));
            }
          }
        }

        // Transform the result to the global reference frame
        transformMatGeoSens(Ks, sKs, t, st, NUM_NODES, NUM_NODES);

        // Add values to the matrix
        addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, n_vars, 4 * NUM_NODES,
                         m_vars, sKs);
        addSensMatValuesTranspose(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars,
                                  4 * NUM_NODES, n_vars, sKs);
      }
    }
  }
}

/*
  Add the sensitivity contribution from the present segment to mass
  matrix
*/
void TACSPanelAnalysis::addMassMatGeoSens(TacsScalar smat[], int dv,
                                          const TacsScalar Z[], int ldz, int nz,
                                          int n1, int n2,
                                          const TacsScalar mass[]) {
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2 * NUM_NODES], Nahp[2 * NUM_NODES], Naahp[2 * NUM_NODES];

  TacsScalar ks[5 * 4 * NUM_NODES], kc[5 * 4 * NUM_NODES];
  TacsScalar sks[5 * 4 * NUM_NODES], skc[5 * 4 * NUM_NODES];
  memset(ks, 0, 20 * NUM_NODES * sizeof(TacsScalar));
  memset(kc, 0, 20 * NUM_NODES * sizeof(TacsScalar));
  memset(sks, 0, 20 * NUM_NODES * sizeof(TacsScalar));
  memset(skc, 0, 20 * NUM_NODES * sizeof(TacsScalar));

  TacsScalar Ms[16 * NUM_NODES * NUM_NODES];
  TacsScalar sMs[16 * NUM_NODES * NUM_NODES];

  TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
  TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
  TacsScalar sc =
      (XptLin[2 * nnodes * dv + 2 * n2] - XptLin[2 * nnodes * dv + 2 * n1]);
  TacsScalar ss = (XptLin[2 * nnodes * dv + 2 * n2 + 1] -
                   XptLin[2 * nnodes * dv + 2 * n1 + 1]);

  TacsScalar Le = sqrt(c * c + s * s);
  TacsScalar sLe = (c * sc + s * ss) / Le;
  sc = (sc * Le - sLe * c) / (Le * Le);
  ss = (ss * Le - sLe * s) / (Le * Le);
  c = c / Le;
  s = s / Le;

  // Transformation between coordinate frames
  // local vars = t * global vars
  TacsScalar t[4], st[4];
  t[0] = c;
  t[1] = s;
  t[2] = -s;
  t[3] = c;
  st[0] = sc;
  st[1] = ss;
  st[2] = -ss;
  st[3] = sc;

  // The variables for the n-th and m-th modes
  int n_vars[4 * NUM_NODES], m_vars[4 * NUM_NODES];

  for (int m = 1; m <= nmodes; m++) {
    memset(Ms, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));
    memset(sMs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

    TacsScalar lambda_m = (M_PI * m) / Lx;

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
      m_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (m - 1) + k];
    }

    // Calculate the stiffness matrix contributions
    for (int k = 0; k < numGauss; k++) {
      // Evaluate the shape functions
      FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

      TacsScalar h = 2.0 / Le;
      TacsScalar sh = -2.0 * sLe / (Le * Le);
      TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;
      TacsScalar sdetJ = 0.5 * gaussWts[k] * sLe * Lx;

      computeMsin(ks, lambda_m, h, N, Nhp, Nahp);
      computeMcos(kc, lambda_m, h, N, Nhp, Nahp);
      computeMsinGeoSens(sks, lambda_m, 0.0, h, sh, N, Nhp, Nahp);
      computeMcosGeoSens(skc, lambda_m, 0.0, h, sh, N, Nhp, Nahp);

      // Compute the sin and cos-coefficient contributions
      for (int i = 0; i < 4 * NUM_NODES; i++) {
        for (int j = 0; j < 4 * NUM_NODES; j++) {
          TacsScalar p = innerMassProduct(mass, &ks[5 * i], &ks[5 * j]);
          Ms[4 * NUM_NODES * i + j] += 0.5 * detJ * p;
          sMs[4 * NUM_NODES * i + j] +=
              0.5 * sdetJ * p +
              detJ * innerMassProduct(mass, &ks[5 * i], &sks[5 * j]);
        }

        for (int j = 0; j < 4 * NUM_NODES; j++) {
          TacsScalar p = innerMassProduct(mass, &kc[5 * i], &kc[5 * j]);
          Ms[4 * NUM_NODES * i + j] += 0.5 * detJ * p;
          sMs[4 * NUM_NODES * i + j] +=
              0.5 * sdetJ * p +
              detJ * innerMassProduct(mass, &kc[5 * i], &skc[5 * j]);
        }
      }
    }

    // Transform the matrices to the global reference frame
    transformMatGeoSens(Ms, sMs, t, st, NUM_NODES, NUM_NODES);
    addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                     m_vars, sMs);

    for (int n = 1; n <= nmodes; n++) {
      if (n % 2 != m % 2) {
        // Determine the variables for the m-th mode
        for (int k = 0; k < 4; k++) {
          n_vars[k] = vars[(4 * nmodes) * n1 + 4 * (n - 1) + k];
          n_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (n - 1) + k];
        }

        TacsScalar lambda_n = (M_PI * n) / Lx;
        TacsScalar sc = (2.0 * n) / (M_PI * (n * n - m * m));

        memset(Ms, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));
        memset(sMs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

        // Compute sin(n) * cos(m) term
        for (int k = 0; k < numGauss; k++) {
          // Evaluate the shape functions
          FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
          FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

          TacsScalar h = 2.0 / Le;
          TacsScalar sh = -2.0 * sLe / (Le * Le);
          TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;
          TacsScalar sdetJ = 0.5 * gaussWts[k] * sLe * Lx;

          computeMsin(ks, lambda_n, h, N, Nhp, Nahp);
          computeMcos(kc, lambda_m, h, N, Nhp, Nahp);
          computeMsinGeoSens(sks, lambda_n, 0.0, h, sh, N, Nhp, Nahp);
          computeMcosGeoSens(skc, lambda_m, 0.0, h, sh, N, Nhp, Nahp);

          // Compute the sine-coefficient contribution
          for (int i = 0; i < 4 * NUM_NODES; i++) {
            for (int j = 0; j < 4 * NUM_NODES; j++) {
              TacsScalar p = innerMassProduct(mass, &ks[5 * i], &kc[5 * j]);
              Ms[4 * NUM_NODES * i + j] += sc * detJ * p;
              sMs[4 * NUM_NODES * i + j] +=
                  sc * sdetJ * p +
                  sc * detJ *
                      (innerMassProduct(mass, &ks[5 * i], &skc[5 * j]) +
                       innerMassProduct(mass, &sks[5 * i], &kc[5 * j]));
            }
          }
        }

        // Transform the result to the global reference frame
        transformMatGeoSens(Ms, sMs, t, st, NUM_NODES, NUM_NODES);

        // Add values to the matrix
        addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, n_vars, 4 * NUM_NODES,
                         m_vars, sMs);
        addSensMatValuesTranspose(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars,
                                  4 * NUM_NODES, n_vars, sMs);
      }
    }
  }
}

/*
  Add the contribution from a segment to the geometric stiffness matrix
*/
void TACSPanelAnalysis::addGeoStiffMatGeoSens(TacsScalar smat[], int dv,
                                              const TacsScalar Z[], int ldz,
                                              int nz, int n1, int n2,
                                              const TacsScalar stress[],
                                              const TacsScalar sstress[]) {
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2 * NUM_NODES], Nahp[2 * NUM_NODES], Naahp[2 * NUM_NODES];

  TacsScalar Gs[16 * NUM_NODES * NUM_NODES];
  TacsScalar sGs[16 * NUM_NODES * NUM_NODES];

  TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
  TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
  TacsScalar sc =
      (XptLin[2 * nnodes * dv + 2 * n2] - XptLin[2 * nnodes * dv + 2 * n1]);
  TacsScalar ss = (XptLin[2 * nnodes * dv + 2 * n2 + 1] -
                   XptLin[2 * nnodes * dv + 2 * n1 + 1]);

  TacsScalar Le = sqrt(c * c + s * s);
  TacsScalar sLe = (c * sc + s * ss) / Le;
  sc = (sc * Le - sLe * c) / (Le * Le);
  ss = (ss * Le - sLe * s) / (Le * Le);
  c = c / Le;
  s = s / Le;

  // Transformation between coordinate frames
  // local vars = t * global vars
  TacsScalar t[4], st[4];
  t[0] = c;
  t[1] = s;
  t[2] = -s;
  t[3] = c;
  st[0] = sc;
  st[1] = ss;
  st[2] = -ss;
  st[3] = sc;

  // The variables for the n-th and m-th modes
  int n_vars[4 * NUM_NODES], m_vars[4 * NUM_NODES];

  for (int m = 1; m <= nmodes; m++) {
    memset(Gs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));
    memset(sGs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

    TacsScalar lambda_m = (M_PI * m) / Lx;

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
      m_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (m - 1) + k];
    }

    // Calculate the stiffness matrix contributions
    for (int k = 0; k < numGauss; k++) {
      // Evaluate the shape functions
      FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

      TacsScalar h = 2.0 / Le;
      TacsScalar sh = -2.0 * sLe / (Le * Le);
      TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;
      TacsScalar sdetJ = 0.5 * gaussWts[k] * sLe * Lx;

      addGs(Gs, 0.5 * detJ, stress, lambda_m, h, N, Na, Nhp, Nahp);
      addGs(sGs, 0.5 * detJ, sstress, lambda_m, h, N, Na, Nhp, Nahp);
      addGsGeoSens(sGs, 0.5 * detJ, 0.5 * sdetJ, stress, lambda_m, 0.0, h, sh,
                   N, Na, Nhp, Nahp);
    }

    transformMatGeoSens(Gs, sGs, t, st, NUM_NODES, NUM_NODES);
    addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                     m_vars, sGs);

    for (int n = 1; n <= nmodes; n++) {
      if (n % 2 != m % 2) {
        // Determine the variables for the n-th mode
        for (int k = 0; k < 4; k++) {
          n_vars[k] = vars[(4 * nmodes) * n1 + 4 * (n - 1) + k];
          n_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (n - 1) + k];
        }

        TacsScalar lambda_n = (M_PI * n) / Lx;
        TacsScalar sc = (2.0 * n) / (M_PI * (n * n - m * m));

        memset(Gs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));
        memset(sGs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

        // Compute int_{0}^{1} sin(n*pi*x) * cos(m*pi*x) dx term
        for (int k = 0; k < numGauss; k++) {
          // Evaluate the shape functions
          FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
          FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

          TacsScalar h = 2.0 / Le;
          TacsScalar sh = -2.0 * sLe / (Le * Le);
          TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;
          TacsScalar sdetJ = 0.5 * gaussWts[k] * sLe * Lx;

          addGcs(Gs, sc * detJ, stress, lambda_n, lambda_m, h, N, Na, Nhp,
                 Nahp);
          addGcs(sGs, sc * detJ, sstress, lambda_n, lambda_m, h, N, Na, Nhp,
                 Nahp);

          addGcsGeoSens(sGs, sc * detJ, sc * sdetJ, stress, lambda_n, 0.0,
                        lambda_m, 0.0, h, sh, N, Na, Nhp, Nahp);
        }

        // Transform the result to the global reference frame
        transformMatGeoSens(Gs, sGs, t, st, NUM_NODES, NUM_NODES);

        // Add values to the matrix
        addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, n_vars, 4 * NUM_NODES,
                         m_vars, sGs);
        addSensMatValuesTranspose(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars,
                                  4 * NUM_NODES, n_vars, sGs);
      }
    }
  }
}

/*
  Add the sensitivity of the stiffness matrix w.r.t. the length of the
  panel
*/
void TACSPanelAnalysis::addStiffMatLxSens(TacsScalar smat[],
                                          const TacsScalar Z[], int ldz, int nz,
                                          int n1, int n2, const TacsScalar At[],
                                          const TacsScalar Bt[],
                                          const TacsScalar Dt[]) {
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2 * NUM_NODES], Nahp[2 * NUM_NODES], Naahp[2 * NUM_NODES];

  TacsScalar stress[6], sstress[6];
  TacsScalar Bs[24 * NUM_NODES], Bc[24 * NUM_NODES];
  TacsScalar sBs[24 * NUM_NODES], sBc[24 * NUM_NODES];

  memset(Bs, 0, 24 * NUM_NODES * sizeof(TacsScalar));
  memset(Bc, 0, 24 * NUM_NODES * sizeof(TacsScalar));
  memset(sBs, 0, 24 * NUM_NODES * sizeof(TacsScalar));
  memset(sBc, 0, 24 * NUM_NODES * sizeof(TacsScalar));

  TacsScalar sKs[16 * NUM_NODES * NUM_NODES];

  TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
  TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);

  TacsScalar Le = sqrt(c * c + s * s);
  c = c / Le;
  s = s / Le;

  // Transformation between coordinate frames
  // local vars = t * global vars
  TacsScalar t[4];
  t[0] = c;
  t[1] = s;
  t[2] = -s;
  t[3] = c;

  // The variables for the n-th and m-th modes
  int n_vars[4 * NUM_NODES], m_vars[4 * NUM_NODES];

  for (int m = 1; m <= nmodes; m++) {
    memset(sKs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

    TacsScalar lambda_m = (M_PI * m) / Lx;
    TacsScalar slambda_m = -(M_PI * m) / (Lx * Lx);

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
      m_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (m - 1) + k];
    }

    // Calculate the stiffness matrix contributions
    for (int k = 0; k < numGauss; k++) {
      // Evaluate the shape functions
      FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

      TacsScalar h = 2.0 / Le;
      TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;
      TacsScalar sdetJ = 0.5 * gaussWts[k] * Le;

      computeBsin(Bs, lambda_m, h, N, Na, Nhp, Nahp, Naahp);
      computeBcos(Bc, lambda_m, h, N, Na, Nhp, Nahp, Naahp);

      computeBsinGeoSens(sBs, lambda_m, slambda_m, h, 0.0, N, Na, Nhp, Nahp,
                         Naahp);
      computeBcosGeoSens(sBc, lambda_m, slambda_m, h, 0.0, N, Na, Nhp, Nahp,
                         Naahp);

      // Compute the sin and cos-coefficient contributions
      for (int i = 0; i < 4 * NUM_NODES; i++) {
        computeStress(stress, &Bs[6 * i], At, Bt, Dt);
        computeStress(sstress, &sBs[6 * i], At, Bt, Dt);
        for (int j = 0; j < 4 * NUM_NODES; j++) {
          TacsScalar p = innerProduct(stress, &Bs[6 * j]);
          sKs[4 * NUM_NODES * i + j] +=
              0.5 * sdetJ * p + 0.5 * detJ *
                                    (innerProduct(sstress, &Bs[6 * j]) +
                                     innerProduct(stress, &sBs[6 * j]));
        }

        computeStress(stress, &Bc[6 * i], At, Bt, Dt);
        computeStress(sstress, &sBc[6 * i], At, Bt, Dt);
        for (int j = 0; j < 4 * NUM_NODES; j++) {
          TacsScalar p = innerProduct(stress, &Bc[6 * j]);
          sKs[4 * NUM_NODES * i + j] +=
              0.5 * sdetJ * p + 0.5 * detJ *
                                    (innerProduct(sstress, &Bc[6 * j]) +
                                     innerProduct(stress, &sBc[6 * j]));
        }
      }
    }

    // Transform the matrices to the global reference frame
    transformMat(sKs, t, NUM_NODES, NUM_NODES);
    addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                     m_vars, sKs);

    for (int n = 1; n <= nmodes; n++) {
      if (n % 2 != m % 2) {
        // Determine the variables for the m-th mode
        for (int k = 0; k < 4; k++) {
          n_vars[k] = vars[(4 * nmodes) * n1 + 4 * (n - 1) + k];
          n_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (n - 1) + k];
        }

        TacsScalar lambda_n = (M_PI * n) / Lx;
        TacsScalar slambda_n = -(M_PI * n) / (Lx * Lx);
        TacsScalar sc = (2.0 * n) / (M_PI * (n * n - m * m));

        memset(sKs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

        // Compute sin(n) * cos(m) term
        for (int k = 0; k < numGauss; k++) {
          // Evaluate the shape functions
          FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
          FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

          TacsScalar h = 2.0 / Le;
          TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;
          TacsScalar sdetJ = 0.5 * gaussWts[k] * Le;

          computeBsin(Bs, lambda_n, h, N, Na, Nhp, Nahp, Naahp);
          computeBcos(Bc, lambda_m, h, N, Na, Nhp, Nahp, Naahp);

          computeBsinGeoSens(sBs, lambda_n, slambda_n, h, 0.0, N, Na, Nhp, Nahp,
                             Naahp);
          computeBcosGeoSens(sBc, lambda_m, slambda_m, h, 0.0, N, Na, Nhp, Nahp,
                             Naahp);

          // Compute the sin-coefficient contribution
          for (int i = 0; i < 4 * NUM_NODES; i++) {
            computeStress(stress, &Bs[6 * i], At, Bt, Dt);
            computeStress(sstress, &sBs[6 * i], At, Bt, Dt);

            for (int j = 0; j < 4 * NUM_NODES; j++) {
              TacsScalar p = innerProduct(stress, &Bc[6 * j]);
              sKs[4 * NUM_NODES * i + j] +=
                  sc * sdetJ * p + sc * detJ *
                                       (innerProduct(sstress, &Bc[6 * j]) +
                                        innerProduct(stress, &sBc[6 * j]));
            }
          }
        }

        // Transform the result to the global reference frame
        transformMat(sKs, t, NUM_NODES, NUM_NODES);

        // Add values to the matrix
        addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, n_vars, 4 * NUM_NODES,
                         m_vars, sKs);
        addSensMatValuesTranspose(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars,
                                  4 * NUM_NODES, n_vars, sKs);
      }
    }
  }
}

/*
  Add the sensitivity contribution from this segment to the mass
  matrix
*/
void TACSPanelAnalysis::addMassMatLxSens(TacsScalar smat[],
                                         const TacsScalar Z[], int ldz, int nz,
                                         int n1, int n2,
                                         const TacsScalar mass[]) {
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2 * NUM_NODES], Nahp[2 * NUM_NODES], Naahp[2 * NUM_NODES];

  TacsScalar ks[5 * 4 * NUM_NODES], kc[5 * 4 * NUM_NODES];
  TacsScalar sks[5 * 4 * NUM_NODES], skc[5 * 4 * NUM_NODES];

  memset(ks, 0, 20 * NUM_NODES * sizeof(TacsScalar));
  memset(kc, 0, 20 * NUM_NODES * sizeof(TacsScalar));
  memset(sks, 0, 20 * NUM_NODES * sizeof(TacsScalar));
  memset(skc, 0, 20 * NUM_NODES * sizeof(TacsScalar));

  TacsScalar sMs[16 * NUM_NODES * NUM_NODES];

  TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
  TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);

  TacsScalar Le = sqrt(c * c + s * s);
  c = c / Le;
  s = s / Le;

  // Transformation between coordinate frames
  // local vars = t * global vars
  TacsScalar t[4];
  t[0] = c;
  t[1] = s;
  t[2] = -s;
  t[3] = c;

  // The variables for the n-th and m-th modes
  int n_vars[4 * NUM_NODES], m_vars[4 * NUM_NODES];

  for (int m = 1; m <= nmodes; m++) {
    memset(sMs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

    TacsScalar lambda_m = (M_PI * m) / Lx;
    TacsScalar slambda_m = -(M_PI * m) / (Lx * Lx);

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
      m_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (m - 1) + k];
    }

    // Calculate the stiffness matrix contributions
    for (int k = 0; k < numGauss; k++) {
      // Evaluate the shape functions
      FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

      TacsScalar h = 2.0 / Le;
      TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;
      TacsScalar sdetJ = 0.5 * gaussWts[k] * Le;

      computeMsin(ks, lambda_m, h, N, Nhp, Nahp);
      computeMcos(kc, lambda_m, h, N, Nhp, Nahp);
      computeMsinGeoSens(sks, lambda_m, slambda_m, h, 0.0, N, Nhp, Nahp);
      computeMcosGeoSens(skc, lambda_m, slambda_m, h, 0.0, N, Nhp, Nahp);

      // Compute the sin and cos-coefficient contributions
      for (int i = 0; i < 4 * NUM_NODES; i++) {
        for (int j = 0; j < 4 * NUM_NODES; j++) {
          TacsScalar p = innerMassProduct(mass, &ks[5 * i], &ks[5 * j]);
          sMs[4 * NUM_NODES * i + j] +=
              0.5 * sdetJ * p +
              detJ * (innerMassProduct(mass, &ks[5 * i], &sks[5 * j]));
        }

        for (int j = 0; j < 4 * NUM_NODES; j++) {
          TacsScalar p = innerMassProduct(mass, &kc[5 * i], &kc[5 * j]);
          sMs[4 * NUM_NODES * i + j] +=
              0.5 * sdetJ * p +
              detJ * (innerMassProduct(mass, &kc[5 * i], &skc[5 * j]));
        }
      }
    }

    // Transform the matrices to the global reference frame
    transformMat(sMs, t, NUM_NODES, NUM_NODES);
    addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                     m_vars, sMs);

    for (int n = 1; n <= nmodes; n++) {
      if (n % 2 != m % 2) {
        // Determine the variables for the m-th mode
        for (int k = 0; k < 4; k++) {
          n_vars[k] = vars[(4 * nmodes) * n1 + 4 * (n - 1) + k];
          n_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (n - 1) + k];
        }

        TacsScalar lambda_n = (M_PI * n) / Lx;
        TacsScalar slambda_n = -(M_PI * n) / (Lx * Lx);
        TacsScalar sc = (2.0 * n) / (M_PI * (n * n - m * m));

        memset(sMs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

        // Compute sin(n) * cos(m) term
        for (int k = 0; k < numGauss; k++) {
          // Evaluate the shape functions
          FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
          FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

          TacsScalar h = 2.0 / Le;
          TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;
          TacsScalar sdetJ = 0.5 * gaussWts[k] * Le;

          computeMsin(ks, lambda_n, h, N, Nhp, Nahp);
          computeMcos(kc, lambda_m, h, N, Nhp, Nahp);
          computeMsinGeoSens(sks, lambda_n, slambda_n, h, 0.0, N, Nhp, Nahp);
          computeMcosGeoSens(skc, lambda_m, slambda_m, h, 0.0, N, Nhp, Nahp);

          // Compute the sin-coefficient contribution
          for (int i = 0; i < 4 * NUM_NODES; i++) {
            for (int j = 0; j < 4 * NUM_NODES; j++) {
              TacsScalar p = innerMassProduct(mass, &kc[5 * i], &kc[5 * j]);
              sMs[4 * NUM_NODES * i + j] +=
                  sc * sdetJ * p +
                  sc * detJ *
                      (innerMassProduct(mass, &ks[5 * i], &skc[5 * j]) +
                       innerMassProduct(mass, &sks[5 * i], &kc[5 * j]));
            }
          }
        }

        // Transform the result to the global reference frame
        transformMat(sMs, t, NUM_NODES, NUM_NODES);

        // Add values to the matrix
        addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, n_vars, 4 * NUM_NODES,
                         m_vars, sMs);
        addSensMatValuesTranspose(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars,
                                  4 * NUM_NODES, n_vars, sMs);
      }
    }
  }
}

/*
  Add the contribution from a segment to the geometric stiffness matrix
*/
void TACSPanelAnalysis::addGeoStiffMatLxSens(TacsScalar smat[],
                                             const TacsScalar Z[], int ldz,
                                             int nz, int n1, int n2,
                                             const TacsScalar stress[]) {
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2 * NUM_NODES], Nahp[2 * NUM_NODES], Naahp[2 * NUM_NODES];

  TacsScalar Gs[16 * NUM_NODES * NUM_NODES];

  TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
  TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);

  TacsScalar Le = sqrt(c * c + s * s);
  c = c / Le;
  s = s / Le;

  // Transformation between coordinate frames
  // local vars = t * global vars
  TacsScalar t[4];
  t[0] = c;
  t[1] = s;
  t[2] = -s;
  t[3] = c;

  // The variables for the n-th and m-th modes
  int n_vars[4 * NUM_NODES], m_vars[4 * NUM_NODES];

  for (int m = 1; m <= nmodes; m++) {
    memset(Gs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

    TacsScalar lambda_m = (M_PI * m) / Lx;
    TacsScalar slambda_m = -(M_PI * m) / (Lx * Lx);

    // Determine the variables for the m-th mode
    for (int k = 0; k < 4; k++) {
      m_vars[k] = vars[(4 * nmodes) * n1 + 4 * (m - 1) + k];
      m_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (m - 1) + k];
    }

    // Calculate the stiffness matrix contributions
    for (int k = 0; k < numGauss; k++) {
      // Evaluate the shape functions
      FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

      TacsScalar h = 2.0 / Le;
      TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;
      TacsScalar sdetJ = 0.5 * gaussWts[k] * Le;

      addGsGeoSens(Gs, 0.5 * detJ, 0.5 * sdetJ, stress, lambda_m, slambda_m, h,
                   0.0, N, Na, Nhp, Nahp);
    }

    transformMat(Gs, t, NUM_NODES, NUM_NODES);
    addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars, 4 * NUM_NODES,
                     m_vars, Gs);

    for (int n = 1; n <= nmodes; n++) {
      if (n % 2 != m % 2) {
        // Determine the variables for the n-th mode
        for (int k = 0; k < 4; k++) {
          n_vars[k] = vars[(4 * nmodes) * n1 + 4 * (n - 1) + k];
          n_vars[4 + k] = vars[(4 * nmodes) * n2 + 4 * (n - 1) + k];
        }

        TacsScalar lambda_n = (M_PI * n) / Lx;
        TacsScalar slambda_n = -(M_PI * n) / (Lx * Lx);
        TacsScalar sc = (2.0 * n) / (M_PI * (n * n - m * m));

        memset(Gs, 0, 16 * NUM_NODES * NUM_NODES * sizeof(TacsScalar));

        // Compute int_{0}^{1} sin(n*pi*x) * cos(m*pi*x) dx term
        for (int k = 0; k < numGauss; k++) {
          // Evaluate the shape functions
          FElibrary::lagrangeSF(N, Na, gaussPts[k], NUM_NODES);
          FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[k]);

          TacsScalar h = 2.0 / Le;
          TacsScalar detJ = 0.5 * gaussWts[k] * Le * Lx;
          TacsScalar sdetJ = 0.5 * gaussWts[k] * Le;

          addGcsGeoSens(Gs, sc * detJ, sc * sdetJ, stress, lambda_n, slambda_n,
                        lambda_m, slambda_m, h, 0.0, N, Na, Nhp, Nahp);
        }

        // Transform the result to the global reference frame
        transformMat(Gs, t, NUM_NODES, NUM_NODES);

        // Add values to the matrix
        addSensMatValues(smat, Z, ldz, nz, 4 * NUM_NODES, n_vars, 4 * NUM_NODES,
                         m_vars, Gs);
        addSensMatValuesTranspose(smat, Z, ldz, nz, 4 * NUM_NODES, m_vars,
                                  4 * NUM_NODES, n_vars, Gs);
      }
    }
  }
}

/*
  Add values to the matrix:
  mat[rows,cols] += values[:,:]

  mat[]: the matrix stored in the following packed format:
  mat_{i,j} = mat[i + j*(j+1)/2]
*/
void TACSPanelAnalysis::addValues(TacsScalar mat[], int kband, int n,
                                  const int rows[], int m, const int cols[],
                                  const TacsScalar values[]) {
  for (int i = 0; i < n; i++) {
    if (rows[i] >= 0) {
      for (int j = 0; j < m; j++) {
        // Store only the upper diagonal
        if (cols[j] >= 0 && cols[j] >= rows[i]) {
          // Banded storage format
          mat[(kband + rows[i] - cols[j]) + cols[j] * (kband + 1)] +=
              values[m * i + j];
        }
      }
    }
  }
}

/*
  Add the transpose of a matrix of values to the matrix mat:
  mat[rows,cols] += values[:,:]^{T}

  mat[]: the matrix stored in the following packed format:
  mat_{i,j} = mat[i + j*(j+1)/2]
*/
void TACSPanelAnalysis::addValuesTranspose(TacsScalar mat[], int kband, int n,
                                           const int rows[], int m,
                                           const int cols[],
                                           const TacsScalar values[]) {
  for (int i = 0; i < n; i++) {
    if (rows[i] >= 0) {
      for (int j = 0; j < m; j++) {
        // Store only the upper diagonal
        if (cols[j] >= 0 && cols[j] >= rows[i]) {
          // Banded storage format
          mat[(kband + rows[i] - cols[j]) + cols[j] * (kband + 1)] +=
              values[n * j + i];
        }
      }
    }
  }
}

/*
  Add values to the sensitivity matrix

  smat += Z^{T} * values * Z
*/
void TACSPanelAnalysis::addSensMatValues(TacsScalar smat[],
                                         const TacsScalar Z[], int ldz, int nz,
                                         int n, const int rows[], int m,
                                         const int cols[],
                                         const TacsScalar values[]) {
  for (int ii = 0; ii < nz; ii++) {
    for (int i = 0; i < n; i++) {
      if (rows[i] >= 0) {
        TacsScalar Zi = Z[ldz * ii + rows[i]];

        TacsScalar val = 0.0;
        for (int j = 0; j < m; j++) {
          if (cols[j] >= 0) {
            val += values[m * i + j] * Z[ldz * ii + cols[j]];
          }
        }
        smat[ii] += Zi * val;
      }
    }
  }
}

/*
  Add values to the sensitivity matrix

  smat += Z^{T} * values^{T} * Z
*/
void TACSPanelAnalysis::addSensMatValuesTranspose(
    TacsScalar smat[], const TacsScalar Z[], int ldz, int nz, int n,
    const int rows[], int m, const int cols[], const TacsScalar values[]) {
  for (int ii = 0; ii < nz; ii++) {
    for (int i = 0; i < n; i++) {
      if (rows[i] >= 0) {
        TacsScalar Zi = Z[ldz * ii + rows[i]];

        TacsScalar val = 0.0;
        for (int j = 0; j < m; j++) {
          if (cols[j] >= 0) {
            val += values[n * j + i] * Z[ldz * ii + cols[j]];
          }
        }
        smat[ii] += Zi * val;
      }
    }
  }
}

/*
  Transform the vector from the local coordinate frame to the
  global coordinate frame.

  Compute T^{T} * Fs
*/
void TACSPanelAnalysis::transformVec(TacsScalar vec[], const TacsScalar t[],
                                     int nn) {
  for (int i = 0; i < nn; i++) {
    TacsScalar fw = vec[4 * i + 2];
    TacsScalar ft = vec[4 * i + 3];

    vec[4 * i + 2] = t[0] * fw + t[2] * ft;
    vec[4 * i + 3] = t[1] * fw + t[3] * ft;
  }
}

/*
  Transform the given matrix from the local coordinate frame to the
  global reference frame.

  Compute T^{T} * mat * T

  Note that here, mat is stored in a dense format.
*/
void TACSPanelAnalysis::transformMat(TacsScalar mat[], const TacsScalar t[],
                                     int nn, int mm) {
  for (int i = 0; i < nn; i++) {
    for (int j = 0; j < mm; j++) {
      // Record the matrix entries
      TacsScalar m[16];
      for (int k = 0; k < 4; k++) {
        m[k] = mat[4 * mm * (4 * i) + 4 * j + k];
        m[12 + k] = mat[4 * mm * (4 * i + 3) + 4 * j + k];
      }

      // Compute m = T^{T} * mat
      for (int k = 0; k < 4; k++) {
        m[4 + k] = (t[0] * mat[4 * mm * (4 * i + 1) + 4 * j + k] +
                    t[2] * mat[4 * mm * (4 * i + 2) + 4 * j + k]);

        m[8 + k] = (t[1] * mat[4 * mm * (4 * i + 1) + 4 * j + k] +
                    t[3] * mat[4 * mm * (4 * i + 2) + 4 * j + k]);
      }

      // Compute mat = m * T
      // [t[0] t[1]]
      // [t[2] t[3]]
      for (int k = 0; k < 4; k++) {
        mat[4 * mm * (4 * i + k) + 4 * j] = m[4 * k];
        mat[4 * mm * (4 * i + k) + 4 * j + 1] =
            m[4 * k + 1] * t[0] + m[4 * k + 2] * t[2];
        mat[4 * mm * (4 * i + k) + 4 * j + 2] =
            m[4 * k + 1] * t[1] + m[4 * k + 2] * t[3];
        mat[4 * mm * (4 * i + k) + 4 * j + 3] = m[4 * k + 3];
      }
    }
  }
}

/*
  Compute the sensitivity of the matrix transformation given above.

  Compute smat and mat such that:

  smat = sT^{T}*mat*T + T^{T}*smat*T + T^{T}*mat*sT
  mat = T^{T}*mat*T

  Note that here, mat is stored in a full dense format.
*/
void TACSPanelAnalysis::transformMatGeoSens(TacsScalar mat[], TacsScalar smat[],
                                            const TacsScalar t[],
                                            const TacsScalar st[], int nn,
                                            int mm) {
  for (int i = 0; i < nn; i++) {
    for (int j = 0; j < mm; j++) {
      // Record the matrix entries
      TacsScalar m[16], sm[16];
      for (int k = 0; k < 4; k++) {
        m[k] = mat[4 * mm * (4 * i) + 4 * j + k];
        m[12 + k] = mat[4 * mm * (4 * i + 3) + 4 * j + k];
        sm[k] = smat[4 * mm * (4 * i) + 4 * j + k];
        sm[12 + k] = smat[4 * mm * (4 * i + 3) + 4 * j + k];
      }

      // Compute m = T^{T} * mat
      for (int k = 0; k < 4; k++) {
        // compute the matrix entries
        m[4 + k] = (t[0] * mat[4 * mm * (4 * i + 1) + 4 * j + k] +
                    t[2] * mat[4 * mm * (4 * i + 2) + 4 * j + k]);

        m[8 + k] = (t[1] * mat[4 * mm * (4 * i + 1) + 4 * j + k] +
                    t[3] * mat[4 * mm * (4 * i + 2) + 4 * j + k]);

        // compute the sensitivity of m
        sm[4 + k] = (t[0] * smat[4 * mm * (4 * i + 1) + 4 * j + k] +
                     t[2] * smat[4 * mm * (4 * i + 2) + 4 * j + k] +
                     st[0] * mat[4 * mm * (4 * i + 1) + 4 * j + k] +
                     st[2] * mat[4 * mm * (4 * i + 2) + 4 * j + k]);

        sm[8 + k] = (t[1] * smat[4 * mm * (4 * i + 1) + 4 * j + k] +
                     t[3] * smat[4 * mm * (4 * i + 2) + 4 * j + k] +
                     st[1] * mat[4 * mm * (4 * i + 1) + 4 * j + k] +
                     st[3] * mat[4 * mm * (4 * i + 2) + 4 * j + k]);
      }

      // Compute mat = m * T
      // [t[0] t[1]]
      // [t[2] t[3]]
      for (int k = 0; k < 4; k++) {
        mat[4 * mm * (4 * i + k) + 4 * j] = m[4 * k];
        mat[4 * mm * (4 * i + k) + 4 * j + 1] =
            m[4 * k + 1] * t[0] + m[4 * k + 2] * t[2];
        mat[4 * mm * (4 * i + k) + 4 * j + 2] =
            m[4 * k + 1] * t[1] + m[4 * k + 2] * t[3];
        mat[4 * mm * (4 * i + k) + 4 * j + 3] = m[4 * k + 3];

        smat[4 * mm * (4 * i + k) + 4 * j] = sm[4 * k];
        smat[4 * mm * (4 * i + k) + 4 * j + 1] =
            (m[4 * k + 1] * st[0] + m[4 * k + 2] * st[2] +
             sm[4 * k + 1] * t[0] + sm[4 * k + 2] * t[2]);
        smat[4 * mm * (4 * i + k) + 4 * j + 2] =
            (m[4 * k + 1] * st[1] + m[4 * k + 2] * st[3] +
             sm[4 * k + 1] * t[1] + sm[4 * k + 2] * t[3]);
        smat[4 * mm * (4 * i + k) + 4 * j + 3] = sm[4 * k + 3];
      }
    }
  }
}

/*
  Compute the derivatives of the strain coefficient of sin w.r.t. to
  the nodal displacements.

  The displacement field takes the form:

  u(x, y) = sum_{n} u_n(y) sin(lambda*(x + beta*y))

  with similar expressions for v and w. As a result, the strain takes
  the form:

  e = Bs_n(y)*sin(lambda_n*(x + beta*y)) + Bc_n(y)*cos(lambda_n*(x + beta*y)).

  This expression is used in conjunction with the constitutive matrix
  to evaluate the stiffness matrix as follows:

  K =
  int_{A} Bs_n^{T}*Bs_n * sin^2(lambda_n*(x + beta*y)) dA +
  int_{A} Bc_n^{T}*Bc_n * cos^2(lambda_n*(x + beta*y)) dA +
  int_{A} Bs_n^{T}*Bc_m * sin(lambda_n*(x + beta*y))*cos(lambda_m*(x + beta*y))
  dA + int_{A} Bs_m^{T}*Bc_n * sin(lambda_m*(x + beta*y))*cos(lambda_n*(x +
  beta*y)) dA

  Note that the integrals can be evaluated by first integrating along
  x then y as follows:

  int_{y} c(y) int_{x} sin^2(lambda_n*(x + beta*y)) dxdy

  where c(y) depends arbitrarily on y. Note that the limits of
  integration are:

  x = -beta*y        and
  x = -beta*y + Lx

  As a result, we can make the substitution: xi = x + beta*y, and the
  lower and upper limits become:

  xi = 0   and
  xi = Lx

  and the integral becomes:

  int_{y} c(y) int_{xi=0}^{xi=Lx} sin^2(lambda_n*xi) dxi dy

  The inner integral is evalued analytically (in this case it is
  simply 1/2), while the second integral in y is evaluated
  numerically.

  This function computes Bs_n at a given Gauss quadrature point and
  can be used to evaluate the stiffness contribution numerically.

  input:
  lambda:  the scalar coefficient for the n-th mode
  h:       2.0/L where L = the length of the segment in the y-z plane
  N:       the in-plane shape functions
  Na:      the derivative of the in-plane shape functions
  Nhp:     the Hermite polynomial shape functions
  Nahp:    the derivative of the Hermite polynomials
  Naahp:   the second derivative of the Hermite polynomials

  output:
  B:       the derivative of the strain w.r.t. the given modal coefficients
*/
void TACSPanelAnalysis::computeBsin(TacsScalar B[], TacsScalar lambda,
                                    TacsScalar h, const double N[],
                                    const double Na[], const double Nhp[],
                                    const double Nahp[], const double Naahp[]) {
  TacsScalar hinv = 1.0 / h;

  for (int k = 0; k < NUM_NODES; k++) {
    // Compute the derivative w.r.t. the displacement coefficients of sin
    // u0
    B[0] = 0.0;
    B[2] = h * Na[k];
    B += 6;

    // v0
    B[1] = h * Na[k];
    B[2] = 0.0;
    B += 6;

    // w0
    B[3] = lambda * lambda * Nhp[2 * k];
    B[4] = -h * h * Naahp[2 * k] + beta * beta * lambda * lambda * Nhp[2 * k];
    B[5] = 2.0 * beta * lambda * lambda * Nhp[2 * k];
    B += 6;

    // theta
    B[3] = hinv * lambda * lambda * Nhp[2 * k + 1];
    B[4] = -h * Naahp[2 * k + 1] +
           hinv * beta * beta * lambda * lambda * Nhp[2 * k + 1];
    B[5] = 2.0 * hinv * beta * lambda * lambda * Nhp[2 * k + 1];
    B += 6;
  }
}

/*
  Compute the derivative of the cos coefficient of the straint with
  respect to the modal amplitudes.

  This function computes Bc_n at a given Gauss quadrature point and
  can be used to evaluate the stiffness contribution numerically. See
  computeBsin() for a more detailed explanation.

  input:
  lambda:  the scalar coefficient for the n-th mode
  h:       2.0/L where L = the length of the segment in the y-z plane
  N:       the in-plane shape functions
  Na:      the derivative of the in-plane shape functions
  Nhp:     the Hermite polynomial shape functions
  Nahp:    the derivative of the Hermite polynomials
  Naahp:   the second derivative of the Hermite polynomials

  output:
  B:       the derivative of the strain w.r.t. the given modal coefficients
*/
void TACSPanelAnalysis::computeBcos(TacsScalar B[], TacsScalar lambda,
                                    TacsScalar h, const double N[],
                                    const double Na[], const double Nhp[],
                                    const double Nahp[], const double Naahp[]) {
  for (int k = 0; k < NUM_NODES; k++) {
    // Compute the derivative w.r.t. the displacement coefficients of sin
    // u0
    B[0] = lambda * N[k];
    B[2] = beta * lambda * N[k];
    B += 6;

    // v0
    B[1] = beta * lambda * N[k];
    B[2] = lambda * N[k];
    B += 6;

    // w0
    B[3] = 0.0;
    B[4] = -2.0 * beta * lambda * h * Nahp[2 * k];
    B[5] = -2.0 * lambda * h * Nahp[2 * k];
    B += 6;

    // theta
    B[3] = 0.0;
    B[4] = -2.0 * beta * lambda * Nahp[2 * k + 1];
    B[5] = -2.0 * lambda * Nahp[2 * k + 1];
    B += 6;
  }
}

/*
  Compute the derivatives of the strain coefficient of sin w.r.t. to
  the nodal displacements
*/
void TACSPanelAnalysis::computeBsinGeoSens(
    TacsScalar B[], TacsScalar lambda, TacsScalar slambda, TacsScalar h,
    TacsScalar sh, const double N[], const double Na[], const double Nhp[],
    const double Nahp[], const double Naahp[]) {
  TacsScalar hinv = 1.0 / h;

  for (int k = 0; k < NUM_NODES; k++) {
    // Compute the derivative w.r.t. the displacement coefficients of sin
    // u0
    B[2] = sh * Na[k];
    B += 6;

    // v0
    B[1] = sh * Na[k];
    B += 6;

    // w0
    B[3] = 2.0 * slambda * lambda * Nhp[2 * k];
    B[4] = -2.0 * sh * h * Naahp[2 * k] +
           2.0 * beta * beta * lambda * slambda * Nhp[2 * k];
    B[5] = 4.0 * beta * lambda * slambda * Nhp[2 * k];
    B += 6;

    // theta
    B[3] =
        hinv * lambda * Nhp[2 * k + 1] * (2.0 * slambda - lambda * hinv * sh);
    B[4] = -sh * Naahp[2 * k + 1] + hinv * beta * beta * lambda *
                                        Nhp[2 * k + 1] *
                                        (2.0 * slambda - lambda * hinv * sh);
    B[5] = 2.0 * hinv * beta * lambda * Nhp[2 * k + 1] *
           (2.0 * slambda - lambda * hinv * sh);
    B += 6;
  }
}

/*
  Compute the derivatives of the strain coefficient of cosine w.r.t.
  the nodal displacements
*/
void TACSPanelAnalysis::computeBcosGeoSens(
    TacsScalar B[], TacsScalar lambda, TacsScalar slambda, TacsScalar h,
    TacsScalar sh, const double N[], const double Na[], const double Nhp[],
    const double Nahp[], const double Naahp[]) {
  for (int k = 0; k < NUM_NODES; k++) {
    // Compute the derivative w.r.t. the displacement coefficients of sin
    // u0
    B[0] = slambda * N[k];
    B[2] = beta * slambda * N[k];
    B += 6;

    // v0
    B[1] = beta * slambda * N[k];
    B[2] = slambda * N[k];
    B += 6;

    // w0
    B[4] = -2.0 * beta * (lambda * sh + slambda * h) * Nahp[2 * k];
    B[5] = -2.0 * (lambda * sh + slambda * h) * Nahp[2 * k];
    B += 6;

    // theta
    B[4] = -2.0 * beta * slambda * Nahp[2 * k + 1];
    B[5] = -2.0 * slambda * Nahp[2 * k + 1];
    B += 6;
  }
}

// Compute the derivative of the strain for a longitudinal beam
// ------------------------------------------------------------
void TACSPanelAnalysis::computeBeamBsin(TacsScalar B[], TacsScalar lambda) {
  // derivative w.r.t. the sin coefficients
  // u
  B += 4;

  // v
  B[2] = -lambda * lambda;
  B += 4;

  // w
  B[1] = -lambda * lambda;
  B += 4;

  // theta
  B += 4;
}

void TACSPanelAnalysis::computeBeamBcos(TacsScalar B[], TacsScalar lambda) {
  // derivative w.r.t. the sin coefficients
  // u
  B[0] = lambda;
  B += 4;

  // v
  B += 4;

  // w
  B += 4;

  // theta
  B[3] = lambda;
  B += 4;
}

/*
  This function adds the contributions from the squared sin and cos
  terms to the geometric stiffness matrix.

  These terms are derived from the nonlinear components of the strain:

  e_x^{nl}  = 0.5*((du/dx)^2 + (dv/dx)^2 + (dw/dx)^2)
  e_y^{nl}  = 0.5*((du/dy)^2 + (dv/dy)^2 + (dw/dy)^2)
  e_xy^{nl} = ((du/dx)*(du/dy) + (dv/dx)*(dv/dy) + (dw/dx)*(dw/dy))

  The assumed displacement field is expanded as a fourier series in
  the axial direction, and a finite-element type discretization along
  the transverse direction. For instance, the axial displacement is
  assumed as follows:

  u(x, y) = u_n(y)*sin(lambda_n*(x + beta*y))

  where x is the axial coordinate and y is the transverse component.
  The derivatives of u(x, y) with respect to x and y-coordinates are:

  du/dx = lambda_n*u_n(y)*cos(lambda_n*(x + beta*y))
  du/dy = u_n'(y)*sin(lambda_n*(x + beta*y)) +
          beta*lambda_n*u_n(y)*cos(lambda_n*(x + beta*y))

  The geometric stiffness contribution can be evaluated as follows:

  G = d^2/du^2 int_{A} Nx*e_x^{nl} + Ny*e_y^{nl} + Nxy*e_xy^{nl} dA

  This function evaluates the contributions from sin and cos squared
  terms in this G matrix.

  input:
  scale:   scalar to multiply the matrix contributions
  s:       the force per unit length along the panel sides
  lambda:  the lambda scalar multiplier
  h:       2/L where L = the panel segment length
  N:       the shape function in the axial direction
  Na:      the derivative of the shape functions
  Nhp:     the Hermite polynomail shape functions for bending
  Nahp:    the derivative of the Hermite polynomials

  output:
  G:       the geometric stiffness matrix
*/
void TACSPanelAnalysis::addGs(TacsScalar G[], TacsScalar scale,
                              const TacsScalar s[], TacsScalar lambda,
                              TacsScalar h, const double N[], const double Na[],
                              const double Nhp[], const double Nahp[]) {
  TacsScalar a =
      scale * lambda * lambda * (s[0] + 2.0 * beta * s[2] + beta * beta * s[1]);
  TacsScalar b = scale * s[1];
  TacsScalar hinv = 1.0 / h;

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      // u
      G[8 * 4 * i + 4 * j] += a * N[i] * N[j] + b * h * h * Na[i] * Na[j];

      // v
      G[8 * (4 * i + 1) + 4 * j + 1] +=
          a * N[i] * N[j] + b * h * h * Na[i] * Na[j];

      // w
      G[8 * (4 * i + 2) + 4 * j + 2] +=
          (a * Nhp[2 * i] * Nhp[2 * j] + b * h * h * Nahp[2 * i] * Nahp[2 * j]);
      G[8 * (4 * i + 2) + 4 * j + 3] +=
          (a * hinv * Nhp[2 * i] * Nhp[2 * j + 1] +
           b * h * Nahp[2 * i] * Nahp[2 * j + 1]);

      // theta
      G[8 * (4 * i + 3) + 4 * j + 2] +=
          (a * hinv * Nhp[2 * i + 1] * Nhp[2 * j] +
           b * h * Nahp[2 * i + 1] * Nahp[2 * j]);
      G[8 * (4 * i + 3) + 4 * j + 3] +=
          (a * hinv * hinv * Nhp[2 * i + 1] * Nhp[2 * j + 1] +
           b * Nahp[2 * i + 1] * Nahp[2 * j + 1]);
    }
  }
}

/*
  Compute the contributions to the geometric stiffness matrix from
  the terms that contain products: sin(lambda_n*x)*cos(lambda_m*x)

  The geometric stiffness contribution is evaluated as follows:

  G = d^2/du^2 int_{A} Nx*e_x^{nl} + Ny*e_y^{nl} + Nxy*e_xy^{nl} dA

  For instance, the derivative of the axial displacement with respect
  to the x/y directions are:

  du/dx = lambda_n*u_n(y)*cos(lambda_n*(x + beta*y))
  du/dy = u_n'(y)*sin(lambda_n*(x + beta*y)) +
          beta*lambda_n*u_n(y)*cos(lambda_n*(x + beta*y))

  If we set Nx = 0, Nxy = 0, and consider only the axial
  contributions, we arrive at this expression:

  Gu = 0.5*Ny*(du/dy)^2

  This term can be written as follows, where summation is
  implied over the n and m indices:

  Gu = 0.5*Ny*(u_n'(y)*sin(xi_n) + beta*lambda_n*u_n(y)*cos(xi_n))*
              (u_m'(y)*sin(xi_m) + beta*lambda_m*u_m(y)*cos(xi_m))

  The term proportional to the required sin/cos term (sin(lambda_n*(x
  + beta*y))*cos(lambda_m*(x + beta*y))) is:

  Gu = 0.5*Ny*beta*lambda_m*u_n'(y)*sin(xi_n)*cos(xi_m)

  input:
  scale:     scalar to multiply the matrix contributions
  s:         the force per unit length along the panel sides
  lambda_n:  the lambda scalar multiplier
  lambda_m:  the other lambda scalar factor
  h:         2/L where L = the panel segment length
  N:         the shape function in the axial direction
  Na:        the derivative of the shape functions
  Nhp:       the Hermite polynomail shape functions for bending
  Nahp:      the derivative of the Hermite polynomials

  output:
  G:       the geometric stiffness matrix
*/
void TACSPanelAnalysis::addGcs(TacsScalar G[], TacsScalar scale,
                               const TacsScalar s[], TacsScalar lambda_n,
                               TacsScalar lambda_m, TacsScalar h,
                               const double N[], const double Na[],
                               const double Nhp[], const double Nahp[]) {
  // Set the scalar multipliers
  TacsScalar g = scale * lambda_m * (s[2] + 0.5 * beta * s[1]);
  TacsScalar hinv = 1.0 / h;

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      // u
      G[8 * 4 * i + 4 * j] += g * h * Na[i] * N[j];

      // v
      G[8 * (4 * i + 1) + 4 * j + 1] += g * h * Na[i] * N[j];

      // w
      G[8 * (4 * i + 2) + 4 * j + 2] += g * h * Nahp[2 * i] * Nhp[2 * j];
      G[8 * (4 * i + 2) + 4 * j + 3] += g * Nahp[2 * i] * Nhp[2 * j + 1];

      // theta
      G[8 * (4 * i + 3) + 4 * j + 2] += g * Nahp[2 * i + 1] * Nhp[2 * j];
      G[8 * (4 * i + 3) + 4 * j + 3] +=
          g * hinv * Nahp[2 * i + 1] * Nhp[2 * j + 1];
    }
  }
}

/*
  Compute the contribution to the derivative of the geometric
  stiffness matrix w.r.t. to the nodal locations.

  This code computes the contribution from the terms proportional to
  sin(lambda_n*(x + beta*y))^2 and cos(lambda_n*(x + beta*y))^2.

  input:
  scale:   the scalar factor
  sscale:  the derivative of the scalar factor
  s:       the force per unit length
  lambda:  the axial fourier coefficient
  slambda: the derivative of the axial fourier coefficient
  h:       2/L where L = the length of the segment
  sh:      the derivative of 2/L
  N:       the axial shape functions
  Na:      the shape derivative of the shape functions
  Nhp:     the Hermite polynomial
  Nahp:    the shape derivative of the Hermite polynomial

  output:
  G:       the geometric stiffness matrix
*/
void TACSPanelAnalysis::addGsGeoSens(TacsScalar G[], TacsScalar scale,
                                     TacsScalar sscale, const TacsScalar s[],
                                     TacsScalar lambda, TacsScalar slambda,
                                     TacsScalar h, TacsScalar sh,
                                     const double N[], const double Na[],
                                     const double Nhp[], const double Nahp[]) {
  TacsScalar a =
      scale * lambda * lambda * (s[0] + 2.0 * beta * s[2] + beta * beta * s[1]);
  TacsScalar sa = (sscale * lambda * lambda + 2.0 * scale * lambda * slambda) *
                  (s[0] + 2.0 * beta * s[2] + beta * beta * s[1]);
  TacsScalar b = scale * s[1];
  TacsScalar sb = sscale * s[1];

  TacsScalar hinv = 1.0 / h;
  TacsScalar shinv = -sh * hinv * hinv;

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      // u
      G[8 * 4 * i + 4 * j] +=
          sa * N[i] * N[j] + (sb * h * h + 2.0 * b * h * sh) * Na[i] * Na[j];

      // v
      G[8 * (4 * i + 1) + 4 * j + 1] +=
          sa * N[i] * N[j] + (sb * h * h + 2.0 * b * h * sh) * Na[i] * Na[j];

      // w
      G[8 * (4 * i + 2) + 4 * j + 2] +=
          sa * Nhp[2 * i] * Nhp[2 * j] +
          (sb * h * h + 2.0 * b * h * sh) * Nahp[2 * i] * Nahp[2 * j];
      G[8 * (4 * i + 2) + 4 * j + 3] +=
          (sa * hinv + a * shinv) * Nhp[2 * i] * Nhp[2 * j + 1] +
          (sb * h + h * sh) * Nahp[2 * i] * Nahp[2 * j + 1];

      // theta
      G[8 * (4 * i + 3) + 4 * j + 2] +=
          (sa * hinv + a * shinv) * Nhp[2 * i + 1] * Nhp[2 * j] +
          (sb * h + h * sh) * Nahp[2 * i] * Nahp[2 * j + 1];
      G[8 * (4 * i + 3) + 4 * j + 3] +=
          (sa * hinv * hinv + 2.0 * a * shinv * hinv) * Nhp[2 * i + 1] *
              Nhp[2 * j + 1] +
          sb * Nahp[2 * i + 1] * Nahp[2 * j + 1];
    }
  }
}

/*
  Add terms to the geometric stiffness matrix that are proportional to
  sin(lambda_n*x)*cos(lambda_m*x)

  This function adds the contribution to the derivative computed by
  the function addGcs().

  input:
  scale:     scalar to multiply the matrix contributions
  sscale:    derivative of the scalar multiplier for the matrix
  s:         the force per unit length along the panel sides
  lambda_n:  the lambda scalar multiplier
  slambda_n: the derivative of the scalar multiplier
  lambda_m:  the other lambda scalar factor
  slambda_m: the derivative of the other lambda scalar factor
  h:         2/L where L = the panel segment length
  N:         the shape function in the axial direction
  Na:        the derivative of the shape functions
  Nhp:       the Hermite polynomail shape functions for bending
  Nahp:      the derivative of the Hermite polynomials

  output:
  G:       the geometric stiffness matrix
*/
void TACSPanelAnalysis::addGcsGeoSens(TacsScalar G[], TacsScalar scale,
                                      TacsScalar sscale, const TacsScalar s[],
                                      TacsScalar lambda_n, TacsScalar slambda_n,
                                      TacsScalar lambda_m, TacsScalar slambda_m,
                                      TacsScalar h, TacsScalar sh,
                                      const double N[], const double Na[],
                                      const double Nhp[], const double Nahp[]) {
  TacsScalar g = scale * lambda_m * (s[2] + 0.5 * beta * s[1]);
  TacsScalar sg =
      (sscale * lambda_m + scale * slambda_m) * (s[2] + 0.5 * beta * s[1]);

  TacsScalar hinv = 1.0 / h;
  TacsScalar shinv = -sh * hinv * hinv;

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      // u
      G[8 * 4 * i + 4 * j] += (sg * h + g * sh) * Na[i] * N[j];

      // v
      G[8 * (4 * i + 1) + 4 * j + 1] += (sg * h + g * sh) * Na[i] * N[j];

      // w
      G[8 * (4 * i + 2) + 4 * j + 2] +=
          (sg * h + g * sh) * Nahp[2 * i] * Nhp[2 * j];
      G[8 * (4 * i + 2) + 4 * j + 3] += sg * Nahp[2 * i] * Nhp[2 * j + 1];

      // theta
      G[8 * (4 * i + 3) + 4 * j + 2] += sg * Nahp[2 * i + 1] * Nhp[2 * j];
      G[8 * (4 * i + 3) + 4 * j + 3] +=
          (sg * hinv + g * shinv) * Nahp[2 * i + 1] * Nhp[2 * j + 1];
    }
  }
}

/*
  Compute the terms required for the mass matrix that are a function
  of sin
*/
void TACSPanelAnalysis::computeMsin(TacsScalar M[], TacsScalar lambda,
                                    TacsScalar h, const double N[],
                                    const double Nhp[], const double Nahp[]) {
  for (int k = 0; k < NUM_NODES; k++) {
    // u0
    M[0] = N[k];
    M += 5;

    // v0
    M[1] = N[k];
    M += 5;

    // w0
    M[2] = Nhp[2 * k];
    M[4] = h * Nahp[2 * k];
    M += 5;

    // theta
    M[2] = Nhp[2 * k + 1] / h;
    M[4] = Nahp[2 * k + 1];
    M += 5;
  }
}

/*
  Compute the terms required for the mass matrix that are a function
  of cos
*/
void TACSPanelAnalysis::computeMcos(TacsScalar M[], TacsScalar lambda,
                                    TacsScalar h, const double N[],
                                    const double Nhp[], const double Nahp[]) {
  for (int k = 0; k < NUM_NODES; k++) {
    // u0
    M += 5;

    // v0
    M += 5;

    // w0
    M[3] = lambda * Nhp[2 * k];
    M[4] = beta * lambda * Nhp[2 * k];
    M += 5;

    // theta
    M[3] = lambda * Nhp[2 * k + 1] / h;
    M[4] = beta * lambda * Nhp[2 * k] / h;
    M += 5;
  }
}

/*
  Compute the terms required for the mass matrix that are a function
  of sin
*/
void TACSPanelAnalysis::computeMsinGeoSens(TacsScalar M[], TacsScalar lambda,
                                           TacsScalar slambda, TacsScalar h,
                                           TacsScalar sh, const double N[],
                                           const double Nhp[],
                                           const double Nahp[]) {
  for (int k = 0; k < NUM_NODES; k++) {
    // u0
    M += 5;

    // v0
    M += 5;

    // w0
    M[4] = sh * Nahp[2 * k];
    M += 5;

    // theta
    M[2] = -sh * Nhp[2 * k + 1] / (h * h);
    M += 5;
  }
}

/*
  Compute the terms required for the mass matrix that are a function
  of cos
*/
void TACSPanelAnalysis::computeMcosGeoSens(TacsScalar M[], TacsScalar lambda,
                                           TacsScalar slambda, TacsScalar h,
                                           TacsScalar sh, const double N[],
                                           const double Nhp[],
                                           const double Nahp[]) {
  for (int k = 0; k < NUM_NODES; k++) {
    // u0
    M += 5;

    // v0
    M += 5;

    // w0
    M[3] = slambda * Nhp[2 * k];
    M[4] = beta * slambda * Nhp[2 * k];
    M += 5;

    // theta
    M[3] = (slambda * h - lambda * sh) * Nhp[2 * k + 1] / (h * h);
    M[4] = beta * (slambda * h - lambda * sh) * Nhp[2 * k + 1] / (h * h);
    M += 5;
  }
}

// Compute the derivative of the strain for a longitudinal beam
// ------------------------------------------------------------
void TACSPanelAnalysis::computeBeamMsin(TacsScalar M[], TacsScalar lambda) {
  // derivative w.r.t. the sin coefficients
  // u
  M[0] = 1.0;
  M += 5;

  // v
  M[1] = 1.0;
  M += 5;

  // w
  M[2] = 1.0;
  M += 5;

  // theta
  M += 5;
}

void TACSPanelAnalysis::computeBeamMcos(TacsScalar M[], TacsScalar lambda) {
  // derivative w.r.t. the sin coefficients
  // u
  M += 5;

  // v
  M[5] = lambda;
  M += 5;

  // w
  M[4] = lambda;
  M += 5;

  // theta
  M += 5;
}

/*
  Compute the segment and beam loads for the problem

  The axial loads are computed as follows:

  A = sum_i L_i * (A11*A22 - A12^2)/A22

  Nx = Nx*Ly
*/
void TACSPanelAnalysis::computeSegmentLoads(TacsScalar Nx, TacsScalar Nxy,
                                            TacsScalar segmentLoads[],
                                            TacsScalar beamLoads[]) {
  double pt[3] = {0.0, 0.0, 0.0};

  memset(segmentLoads, 0, 3 * nsegments * sizeof(TacsScalar));
  memset(beamLoads, 0, nbeams * sizeof(TacsScalar));

  TacsScalar cy = Xpts[2 * last_node] - Xpts[2 * first_node];
  TacsScalar sy = Xpts[2 * last_node + 1] - Xpts[2 * first_node + 1];
  TacsScalar Ly = sqrt(cy * cy + sy * sy);

  // Sum up the stiffness in the axial direction
  TacsScalar EA = 0.0, G = 0.0;
  for (int k = 0; k < nsegments; k++) {
    int n1 = nodes[2 * k];
    int n2 = nodes[2 * k + 1];

    TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
    TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
    TacsScalar Le = sqrt(c * c + s * s);
    c = c / Le;

    TacsScalar At[6], Bt[6], Dt[6], Ast[3];
    panels[k]->getStiffness(pt, At, Bt, Dt, Ast);
    TacsScalar A11 = (At[0] * At[3] - At[1] * At[1]) / At[3];

    EA += A11 * Le;
    G += At[5] * Le * c * c;
  }

  // Add the contributions from the beam elements
  for (int k = 0; k < nbeams; k++) {
    TacsScalar Ct[10];
    beams[k]->getStiffness(pt, Ct);
    EA += Ct[0];
  }

  // Compute the axial and shear strain
  TacsScalar epx = Nx * Ly / EA;
  TacsScalar gamma = Nxy * Ly / G;

  // Set the stresses
  for (int k = 0; k < nsegments; k++) {
    int n1 = nodes[2 * k];
    int n2 = nodes[2 * k + 1];

    TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
    TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
    TacsScalar Le = sqrt(c * c + s * s);
    c = c / Le;

    TacsScalar At[6], Bt[6], Dt[6], Ast[3];
    panels[k]->getStiffness(pt, At, Bt, Dt, Ast);
    TacsScalar A11 = (At[0] * At[3] - At[1] * At[1]) / At[3];

    segmentLoads[3 * k] = A11 * epx;
    segmentLoads[3 * k + 2] = At[5] * c * c * gamma;
  }

  for (int k = 0; k < nbeams; k++) {
    TacsScalar Ct[10];
    beams[k]->getStiffness(pt, Ct);
    beamLoads[k] = Ct[0] * epx;
  }
}

/*
  Determine the sensitivity of the segment loads to a given design variable
*/
/*
void TACSPanelAnalysis::computeSegmentLoadsDVSens( int dvNum,
                                                   TacsScalar Nx,
                                                   TacsScalar Nxy,
                                                   TacsScalar segmentLoads[],
                                                   TacsScalar beamLoads[] ){
  double pt[3] = {0.0, 0.0, 0.0};

  memset(segmentLoads, 0, 3*nsegments*sizeof(TacsScalar));
  memset(beamLoads, 0, nbeams*sizeof(TacsScalar));

  TacsScalar cy = Xpts[2*last_node] - Xpts[2*first_node];
  TacsScalar sy = Xpts[2*last_node+1] - Xpts[2*first_node+1];
  TacsScalar Ly = sqrt(cy*cy + sy*sy);

  // Sum up the stiffness in the axial direction
  TacsScalar EA = 0.0, G = 0.0;
  TacsScalar sEA = 0.0, sG = 0.0;
  for ( int k = 0; k < nsegments; k++ ){
    int n1 = nodes[2*k];
    int n2 = nodes[2*k+1];

    TacsScalar c = (Xpts[2*n2] - Xpts[2*n1]);
    TacsScalar s = (Xpts[2*n2+1] - Xpts[2*n1+1]);
    TacsScalar Le = sqrt(c*c + s*s);
    c = c/Le;

    TacsScalar At[6], Bt[6], Dt[6], Ast[3];
    TacsScalar sAt[6], sBt[6], sDt[6], sAst[3];
    panels[k]->getStiffness(pt, At, Bt, Dt, Ast);
    panels[k]->getStiffnessDVSens(dvNum, pt, sAt, sBt, sDt, sAst);
    TacsScalar A11 = (At[0]*At[3] - At[1]*At[1])/At[3];
    TacsScalar sA11 = (At[0]*sAt[3] + sAt[0]*At[3] - 2.0*At[1]*sAt[1])/At[3];
    sA11 -= (A11/At[3])*sAt[3];

    EA += A11*Le;
    G += At[5]*Le*c*c;

    sEA += sA11*Le;
    sG += sAt[5]*Le*c*c;
  }

  for ( int k = 0; k < nbeams; k++ ){
    TacsScalar Ct[10], sCt[10];
    beams[k]->getStiffness(pt, Ct);
    beams[k]->getStiffnessDVSens(dvNum, pt, sCt);
    EA += Ct[0];
    sEA += sCt[0];
  }

  TacsScalar epx = Nx*Ly/EA;
  TacsScalar gamma = Nxy*Ly/G;

  TacsScalar sepx = -sEA*Nx*Ly/(EA*EA);
  TacsScalar sgamma = -sG*Nxy*Ly/(G*G);

  // Set the stiffnesses
  for ( int k = 0; k < nsegments; k++ ){
    int n1 = nodes[2*k];
    int n2 = nodes[2*k+1];

    TacsScalar c = (Xpts[2*n2] - Xpts[2*n1]);
    TacsScalar s = (Xpts[2*n2+1] - Xpts[2*n1+1]);
    TacsScalar Le = sqrt(c*c + s*s);
    c = c/Le;

    TacsScalar At[6], Bt[6], Dt[6], Ast[3];
    TacsScalar sAt[6], sBt[6], sDt[6], sAst[3];
    panels[k]->getStiffness(pt, At, Bt, Dt, Ast);
    panels[k]->getStiffnessDVSens(dvNum, pt, sAt, sBt, sDt, sAst);
    TacsScalar A11 = (At[0]*At[3] - At[1]*At[1])/At[3];
    TacsScalar sA11 = (At[0]*sAt[3] + sAt[0]*At[3] - 2.0*At[1]*sAt[1])/At[3];
    sA11 -= (A11/At[3])*sAt[3];

    segmentLoads[3*k] = sA11*epx + A11*sepx;
    segmentLoads[3*k+2] = sAt[5]*c*c*gamma + At[5]*c*c*sgamma;
  }

  for ( int k = 0; k < nbeams; k++ ){
    TacsScalar Ct[10], sCt[10];
    beams[k]->getStiffness(pt, Ct);
    beams[k]->getStiffnessDVSens(dvNum, pt, sCt);
    beamLoads[k] = sCt[0]*epx + Ct[0]*sepx;
  }
}
*/
/*
  Determine the sensitivity of the segment loads to a given geometric
  design variable
*/
void TACSPanelAnalysis::computeSegmentLoadsGeoSens(TacsScalar Nx,
                                                   TacsScalar Nxy, int dv,
                                                   TacsScalar segmentLoads[],
                                                   TacsScalar beamLoads[]) {
  double pt[3] = {0.0, 0.0, 0.0};

  memset(segmentLoads, 0, 3 * nsegments * sizeof(TacsScalar));
  memset(beamLoads, 0, nbeams * sizeof(TacsScalar));

  TacsScalar cy = Xpts[2 * last_node] - Xpts[2 * first_node];
  TacsScalar sy = Xpts[2 * last_node + 1] - Xpts[2 * first_node + 1];
  TacsScalar Ly = sqrt(cy * cy + sy * sy);

  TacsScalar scy = (XptLin[2 * nnodes * dv + 2 * last_node] -
                    XptLin[2 * nnodes * dv + 2 * first_node]);
  TacsScalar ssy = (XptLin[2 * nnodes * dv + 2 * last_node + 1] -
                    XptLin[2 * nnodes * dv + 2 * first_node + 1]);
  TacsScalar sLy = (cy * scy + sy * ssy) / Ly;

  // Sum up the stiffness in the axial direction
  TacsScalar EA = 0.0, G = 0.0;
  TacsScalar sEA = 0.0, sG = 0.0;
  for (int k = 0; k < nsegments; k++) {
    int n1 = nodes[2 * k];
    int n2 = nodes[2 * k + 1];

    TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
    TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
    TacsScalar sc =
        (XptLin[2 * nnodes * dv + 2 * n2] - XptLin[2 * nnodes * dv + 2 * n1]);
    TacsScalar ss = (XptLin[2 * nnodes * dv + 2 * n2 + 1] -
                     XptLin[2 * nnodes * dv + 2 * n1 + 1]);

    TacsScalar Le = sqrt(c * c + s * s);
    TacsScalar sLe = (c * sc + s * ss) / Le;
    sc = (sc * Le - sLe * c) / (Le * Le);
    c = c / Le;

    TacsScalar At[6], Bt[6], Dt[6], Ast[3];
    panels[k]->getStiffness(pt, At, Bt, Dt, Ast);
    TacsScalar A11 = (At[0] * At[3] - At[1] * At[1]) / At[3];

    EA += A11 * Le;
    G += At[5] * Le * c * c;

    sEA += A11 * sLe;
    sG += At[5] * (sLe * c * c + 2.0 * Le * c * sc);
  }

  for (int k = 0; k < nbeams; k++) {
    TacsScalar Ct[10];
    beams[k]->getStiffness(pt, Ct);
    EA += Ct[0];
  }

  TacsScalar gamma = Nxy * Ly / G;

  TacsScalar sepx = Nx * (EA * sLy - sEA * Ly) / (EA * EA);
  TacsScalar sgamma = Nxy * (G * sLy - sG * Ly) / (G * G);

  // Set the stiffnesses
  for (int k = 0; k < nsegments; k++) {
    int n1 = nodes[2 * k];
    int n2 = nodes[2 * k + 1];

    TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
    TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);
    TacsScalar sc =
        (XptLin[2 * nnodes * dv + 2 * n2] - XptLin[2 * nnodes * dv + 2 * n1]);
    TacsScalar ss = (XptLin[2 * nnodes * dv + 2 * n2 + 1] -
                     XptLin[2 * nnodes * dv + 2 * n1 + 1]);

    TacsScalar Le = sqrt(c * c + s * s);
    TacsScalar sLe = (c * sc + s * ss) / Le;
    sc = (sc * Le - sLe * c) / (Le * Le);
    c = c / Le;

    TacsScalar At[6], Bt[6], Dt[6], Ast[3];
    panels[k]->getStiffness(pt, At, Bt, Dt, Ast);
    TacsScalar A11 = (At[0] * At[3] - At[1] * At[1]) / At[3];
    segmentLoads[3 * k] = A11 * sepx;
    segmentLoads[3 * k + 1] = 0.0;
    segmentLoads[3 * k + 2] = At[5] * (c * c * sgamma + 2.0 * sc * c * gamma);
  }
}

/*
  Print the panel mode to a file.

  This code can be used to write out, buckling or natural modes of the
  stiffened panel, or even the displaced panel shape due to a pressure
  load.

  intput:
  file_name:  the name of the output file
  x:          the variables
  nx:         the level of the discretization in the x-direction
*/
void TACSPanelAnalysis::printPanelMode(const char *file_name,
                                       const TacsScalar *x, int nx) {
  // Go through each segment and plot the
  int ns = 5;  // Number of segments per panel to use

  int tnodes = (ns + 1) * nsegments * (nx + 1);
  int telems = ns * nsegments * nx;

  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2 * NUM_NODES], Nahp[2 * NUM_NODES], Naahp[2 * NUM_NODES];

  FILE *fp = fopen(file_name, "w");
  if (fp) {
    TacsScalar *local_vars = new TacsScalar[8 * nmodes];
    TacsScalar *v1 = &local_vars[0];
    TacsScalar *v2 = &local_vars[4 * nmodes];

    fprintf(fp, "Variables = x, y, z, u, v, w\n");
    fprintf(fp, "Zone T = Panel, N=%d, E=%d, ", tnodes, telems);
    fprintf(fp, "DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n");

    // Print out the nodal locations and variables
    for (int k = 0; k < nsegments; k++) {
      int n1 = nodes[2 * k];
      int n2 = nodes[2 * k + 1];

      TacsScalar c = (Xpts[2 * n2] - Xpts[2 * n1]);
      TacsScalar s = (Xpts[2 * n2 + 1] - Xpts[2 * n1 + 1]);

      TacsScalar Le = sqrt(c * c + s * s);
      c = c / Le;
      s = s / Le;

      // Transformation between coordinate frames
      // local vars = t * global vars
      TacsScalar t[4];
      t[0] = c;
      t[1] = s;
      t[2] = -s;
      t[3] = c;

      // Compute the local variables from the set of global variables
      for (int i = 0; i < 4 * nmodes; i++) {
        v1[i] = 0.0;
        v2[i] = 0.0;

        int var = vars[(4 * nmodes) * n1 + i];
        if (var >= 0) {
          v1[i] = x[var];
        }

        var = vars[(4 * nmodes) * n2 + i];
        if (var >= 0) {
          v2[i] = x[var];
        }
      }

      // Transform the local variables to the local coordiante system
      for (int i = 0; i < 2 * nmodes; i++) {
        TacsScalar v = local_vars[4 * i + 1];
        TacsScalar w = local_vars[4 * i + 2];

        local_vars[4 * i + 1] = t[0] * v + t[1] * w;
        local_vars[4 * i + 2] = t[2] * v + t[3] * w;
      }

      TacsScalar h = 2.0 / Le;

      for (int i = 0; i < nx + 1; i++) {
        for (int j = 0; j < ns + 1; j++) {
          TacsScalar x[3];
          double xi = -1.0 + (2.0 * j) / ns;
          x[1] = 0.5 * (Xpts[2 * n1] * (1.0 - xi) + Xpts[2 * n2] * (1.0 + xi));
          x[2] = 0.5 * (Xpts[2 * n1 + 1] * (1.0 - xi) +
                        Xpts[2 * n2 + 1] * (1.0 + xi));
          x[0] = (Lx * i) / nx - beta * x[1];

          FElibrary::lagrangeSF(N, Na, xi, NUM_NODES);
          FElibrary::cubicHP(Nhp, Nahp, Naahp, xi);

          TacsScalar u[3] = {0.0, 0.0, 0.0};

          // Compute the local variable values
          for (int n = 1; n <= nmodes; n++) {
            TacsScalar lambda_n = (M_PI * n) / Lx;

            TacsScalar us[3];
            us[0] = N[0] * v1[4 * (n - 1)] + N[1] * v2[4 * (n - 1)];
            us[1] = N[0] * v1[4 * (n - 1) + 1] + N[1] * v2[4 * (n - 1) + 1];
            us[2] = (Nhp[0] * v1[4 * (n - 1) + 2] +
                     Nhp[1] * v1[4 * (n - 1) + 3] / h +
                     Nhp[2] * v2[4 * (n - 1) + 2] +
                     Nhp[3] * v2[4 * (n - 1) + 3] / h);
            s = sin(lambda_n * (x[0] + beta * x[1]));

            u[0] += s * us[0];
            u[1] += s * us[1];
            u[2] += s * us[2];
          }

          TacsScalar v = u[1];
          TacsScalar w = u[2];
          u[1] = t[0] * v + t[2] * w;
          u[2] = t[1] * v + t[3] * w;

          fprintf(fp, "%e %e %e %e %e %e\n", x[0], x[1], x[2], u[0], u[1],
                  u[2]);
        }
      }
    }

    delete[] local_vars;

    // Print out the connectivity
    int n = 1;
    for (int k = 0; k < nsegments; k++) {
      for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ns; j++) {
          fprintf(fp, "%d %d %d %d\n", n + i * (ns + 1) + j,
                  n + i * (ns + 1) + j + 1, n + (i + 1) * (ns + 1) + j + 1,
                  n + (i + 1) * (ns + 1) + j);
        }
      }

      n += (nx + 1) * (ns + 1);
    }

    fclose(fp);
  }
}

/*
  Code for testing the sensitivity of the various calculations
*/
/*
void TACSPanelAnalysis::testStiffnessDVSens( double dh ){
  int ndvs = 0;

  for ( int k = 0; k < numDesignVars; k++ ){
    if (designVarNums[k] >= ndvs){
      ndvs = designVarNums[k];
    }
  }
  ndvs++;

  TacsScalar * x = new TacsScalar[ ndvs ];
  memset(x, 0, ndvs*sizeof(TacsScalar));
  getDesignVars(x, ndvs);

  // Check the sensitivity of the stiffness calculations
  TacsScalar fA[6], fB[6], fD[6], fAs[3];
  TacsScalar bA[6], bB[6], bD[6], bAs[3];

  TacsScalar sA[6], sB[6], sD[6], sAs[3];

  for ( int k = 0; k < numDesignVars; k++ ){
    int dvNum = designVarNums[k];
    printf("Testing sensitivity w.r.t. design variable %d\n", dvNum);

    setDesignVars(x, ndvs);
    computeStiffnessDVSens(dvNum, sA, sB, sD, sAs);

    TacsScalar xt = x[dvNum];
    x[dvNum] = xt + dh;
    setDesignVars(x, ndvs);
    computeStiffness(fA, fB, fD, fAs);

    x[dvNum] = xt - dh;
    setDesignVars(x, ndvs);
    computeStiffness(bA, bB, bD, bAs);

    x[dvNum] = xt;

    for ( int i = 0; i < 6; i++ ){
      TacsScalar fd = 0.5*(fA[i] - bA[i])/dh;
      if (fd != 0.0 || sA[i] != 0.0){
        TacsScalar rel_err = (fd - sA[i])/fd;
        printf(" sA[%2d] = %15.6e  FD[%2d] = %15.6e Rel Err: %10.2e\n",
               i, sA[i], i, fd, rel_err);
      }
    }

    for ( int i = 0; i < 6; i++ ){
      TacsScalar fd = 0.5*(fB[i] - bB[i])/dh;
      if (fd != 0.0 || sB[i] != 0.0){
        TacsScalar rel_err = (fd - sB[i])/fd;
        printf(" sB[%2d] = %15.6e  FD[%2d] = %15.6e Rel Err: %10.2e\n",
               i, sB[i], i, fd, rel_err);
      }
    }

    for ( int i = 0; i < 6; i++ ){
      TacsScalar fd = 0.5*(fD[i] - bD[i])/dh;
      if (fd != 0.0 || sD[i] != 0.0){
        TacsScalar rel_err = (fd - sD[i])/fd;
        printf(" sD[%2d] = %15.6e  FD[%2d] = %15.6e Rel Err: %10.2e\n",
               i, sD[i], i, fd, rel_err);
      }
    }

    for ( int i = 0; i < 3; i++ ){
      TacsScalar fd = 0.5*(fAs[i] - bAs[i])/dh;
      if (fd != 0.0 || sAs[i] != 0.0){
        TacsScalar rel_err = (fd - sAs[i])/fd;
        printf("sAs[%2d] = %15.6e  FD[%2d] = %15.6e Rel Err: %10.2e\n",
               i, sAs[i], i, fd, rel_err);
      }
    }
  }

  delete [] x;
}
*/
/*
  Test the implementation of the failure criteria with respect to
  the design variables and strains

  input:
  dh:      the finite-difference interval
  strain:  the strain to use as the test
  weights: weights corresponding to each failure point (len(weights) = nfail)
  nfail:   the number of failure points
*/
/*
void TACSPanelAnalysis::testFailureSens( double dh,
                                         const TacsScalar strain[],
                                         const TacsScalar weights[],
                                         int nfail ){
  int ndvs = 0;

  for ( int k = 0; k < numDesignVars; k++ ){
    if (designVarNums[k] >= ndvs){
      ndvs = designVarNums[k];
    }
  }
  ndvs++;

  TacsScalar * x = new TacsScalar[ ndvs ];
  memset(x, 0, ndvs*sizeof(TacsScalar));
  getDesignVars(x, ndvs);

  TacsScalar * fail = new TacsScalar[nfail];
  TacsScalar strainSens[8];

  printf("failureDVSens\n");

  for ( int k = 0; k < numDesignVars; k++ ){
    int dvNum = designVarNums[k];

    setDesignVars(x, ndvs);
    TacsScalar failSens;
    failureDVSens(dvNum, strain, weights, nfail, &failSens);

    TacsScalar xt = x[dvNum];
    x[dvNum] = xt + dh;
    setDesignVars(x, ndvs);
    failure(strain, fail, nfail);
    TacsScalar ffval = 0.0;
    for ( int i = 0; i < nfail; i++ ){
      ffval += weights[i]*fail[i];
    }

    x[dvNum] = xt - dh;
    setDesignVars(x, ndvs);
    failure(strain, fail, nfail);
    TacsScalar bfval = 0.0;
    for ( int i = 0; i < nfail; i++ ){
      bfval += weights[i]*fail[i];
    }

    x[dvNum] = xt;

    TacsScalar fd = 0.5*(ffval - bfval)/dh;
    TacsScalar rel_err = (fd - failSens)/fd;
    printf("Analytic[%2d] = %15.6e  FD[%2d] = %15.6e Rel Err: %10.2e\n",
           k, failSens, k, fd, rel_err);
  }

  failureStrainSens(strain, weights, nfail, strainSens);

  printf("failureStrainSens\n");
  TacsScalar e[8];
  memcpy(e, strain, 8*sizeof(TacsScalar));

  for ( int k = 0; k < 8; k++ ){
    TacsScalar xt = e[k];
    e[k] = xt + dh;
    failure(e, fail, nfail);
    TacsScalar ffval = 0.0;
    for ( int i = 0; i < nfail; i++ ){
      ffval += weights[i]*fail[i];
    }

    e[k] = xt - dh;
    failure(e, fail, nfail);
    TacsScalar bfval = 0.0;
    for ( int i = 0; i < nfail; i++ ){
      bfval += weights[i]*fail[i];
    }

    e[k] = xt;

    TacsScalar fd = 0.5*(ffval - bfval)/dh;
    TacsScalar rel_err = (fd - strainSens[k])/fd;
    printf("Analytic[%2d] = %15.6e  FD[%2d] = %15.6e Rel Err: %10.2e\n",
           k, strainSens[k], k, fd, rel_err);
  }

  delete [] x;
  delete [] fail;
}
*/
/*
  This code tests the implementation of the derivative of the critical
  buckling loads with respect to the design variables using a
  finite-difference method. This code tests whatever design variables
  are set in the problem - constitutive, panel length and geometric
  variables.

  input:
  dh:     the step size
  Nx:     the axial panel load per unit length
  Nxy:    the shear load per unit length
  nloads: the number of critical buckling loads
*/
/*
void TACSPanelAnalysis::testBucklingDVSens( double dh,
                                            TacsScalar Nx,
                                            TacsScalar Nxy,
                                            int nloads ){
  // Find the maximum number of design variables
  int ndvs = 0;
  for ( int k = 0; k < numDesignVars; k++ ){
    if (designVarNums[k] >= ndvs){
      ndvs = designVarNums[k];
    }
  }
  ndvs++;

  TacsScalar * x = new TacsScalar[ ndvs ];
  memset(x, 0, ndvs*sizeof(TacsScalar));
  getDesignVars(x, ndvs);

  TacsScalar * fposLoads = new TacsScalar[nloads];
  TacsScalar * bposLoads = new TacsScalar[nloads];
  TacsScalar * posLoadDVSens = new TacsScalar[nloads*numDesignVars];

  TacsScalar * fnegLoads = new TacsScalar[nloads];
  TacsScalar * bnegLoads = new TacsScalar[nloads];
  TacsScalar * negLoadDVSens = new TacsScalar[nloads*numDesignVars];

  TacsScalar * beamLoads = new TacsScalar[nbeams];
  TacsScalar * fsegLoads = new TacsScalar[3*nsegments];
  TacsScalar * bsegLoads = new TacsScalar[3*nsegments];
  TacsScalar * segLoadDVSens = new TacsScalar[3*nsegments];

  if (Nx == 0.0){
    computeBucklingLoadsDVSens(Nxy, fposLoads, fnegLoads,
                               posLoadDVSens, negLoadDVSens, nloads);
  }
  else {
    computeBucklingLoadsDVSens(Nx, Nxy, fposLoads,
                               posLoadDVSens, nloads);
  }

  printf("Test buckling design variable sensitivities, beta = %5.3f\n",
         beta);

  for ( int k = 0; k < numDesignVars; k++ ){
    int dvNum = designVarNums[k];
    printf("Testing sensitivity w.r.t. design variable %d type = %d\n",
           dvNum, designVarTypes[k]);

    computeSegmentLoads(Nx, Nxy, fsegLoads, beamLoads);
    memset(segLoadDVSens, 0, 3*nsegments*sizeof(TacsScalar));

    if (designVarTypes[k] == 1){
      int dv = 0;
      for ( int j = 0; j < nDvGeo; j++ ){
        if (geoDvNums[j] == dvNum){
          dv = j;
          break;
        }
      }

      computeSegmentLoadsGeoSens(Nx, Nxy, dv,
                                 segLoadDVSens, beamLoads);
    }
    else if (designVarTypes[k] != 2){
      computeSegmentLoadsDVSens(dvNum, Nx, Nxy,
                                segLoadDVSens, beamLoads);
    }

    TacsScalar xt = x[dvNum];
    x[dvNum] = xt + dh;
    setDesignVars(x, ndvs);
    computeSegmentLoads(Nx, Nxy, fsegLoads, beamLoads);
    if (Nx == 0.0){
      computeBucklingLoads(Nxy, fposLoads, fnegLoads, nloads);
    }
    else {
      computeBucklingLoads(Nx, Nxy, fposLoads, nloads);
    }

    x[dvNum] = xt - dh;
    setDesignVars(x, ndvs);
    computeSegmentLoads(Nx, Nxy, bsegLoads, beamLoads);
    if (Nx == 0.0){
      computeBucklingLoads(Nxy, bposLoads, bnegLoads, nloads);
    }
    else {
      computeBucklingLoads(Nx, Nxy, bposLoads, nloads);
    }

    x[dvNum] = xt;
    setDesignVars(x, ndvs);

    for ( int i = 0; i < 3*nsegments; i++ ){
      TacsScalar fd = 0.5*(fsegLoads[i] - bsegLoads[i])/dh;
      if (fd != 0.0 || segLoadDVSens[i] != 0.0){
        TacsScalar rel_err = (fd - segLoadDVSens[i])/fd;
        printf("segmentDVSens[%3d] = %15.6f  FD[%3d] = %15.6f Rel Err:
%10.2e\n", i, segLoadDVSens[i], i, fd, rel_err);
      }
    }

    for ( int i = 0; i < nloads; i++ ){
      TacsScalar fd = 0.5*(fposLoads[i] - bposLoads[i])/dh;
      TacsScalar rel_err = (fd - posLoadDVSens[numDesignVars*i + k])/fd;
      printf("posDVSens[%2d]: %15.6e  FD[%2d]: %15.6e Rel Err: %10.2e\n",
             i, posLoadDVSens[numDesignVars*i + k], i, fd, rel_err);
    }

    if (Nx == 0.0){
      for ( int i = 0; i < nloads; i++ ){
        TacsScalar fd = 0.5*(fnegLoads[i] - bnegLoads[i])/dh;
        TacsScalar rel_err = (fd - negLoadDVSens[numDesignVars*i + k])/fd;
        printf("negDVSens[%2d]: %15.6e  FD[%2d]: %15.6e Rel Err: %10.2e\n",
               i, negLoadDVSens[numDesignVars*i + k], i, fd, rel_err);
      }
    }
  }

  delete [] x;
  delete [] fposLoads;
  delete [] bposLoads;
  delete [] posLoadDVSens;

  delete [] fnegLoads;
  delete [] bnegLoads;
  delete [] negLoadDVSens;
}
*/

#endif  // TACS_USE_COMPLEX
