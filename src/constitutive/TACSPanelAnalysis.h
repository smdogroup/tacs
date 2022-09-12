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

#ifndef TACS_PANEL_ANALYSIS_H
#define TACS_PANEL_ANALYSIS_H

/*
  This code is not compatible with the complex verion of
  TACS. Therefore, when TACS_USE_COMPLEX is defined, we remove the
  TACSPanelAnslysis class from TACS itself.
*/
#ifdef TACS_USE_COMPLEX
#warning "TACSPanelAnalysis cannot be used with complex TACS"
#else

#include "EBStiffness.h"
#include "FSDTStiffness.h"

/*
  This class implements a panel-level buckling analysis for
  determining the buckling loads of the panel and stiffeners.
*/
class TACSPanelAnalysis : public TACSObject {
 public:
  static const int SKIN_SEGMENT = 0;
  static const int STIFFENER_SEGMENT = 1;

  TACSPanelAnalysis(int _nnodes, int _nsegments, int _nbeams, int _nmodes,
                    TacsScalar Lx, TacsScalar theta = 0.0);
  ~TACSPanelAnalysis();

  // Functions for initialization
  // ----------------------------
  void setPoints(TacsScalar *Xpts, int npoints);
  void setSegment(int seg, int seg_type, FSDTStiffness *stiff, int n1, int n2);
  void setBeam(int beam, EBStiffness *stiff, int n);
  void setFirstNodeBC(int _first_node, int _first_node_bc);
  void setLastNodeBC(int _last_node, int _last_node_bc);
  void initialize();

  // Compute the mass per unit area and its derivative
  // -------------------------------------------------
  TacsScalar computePtMass();
  void addPtMassDVSens(TacsScalar scale, TacsScalar fdvSens[], int numDVs);

  // Compute the KS failure function over the segments
  // -------------------------------------------------
  void setFailurePoints(int _npoints, int segments[], int nodes[],
                        int is_skin[]);
  void failure(const TacsScalar strain[], TacsScalar fail[], int nfail);
  void addFailureDVSens(const TacsScalar strain[], const TacsScalar weights[],
                        int nfail, TacsScalar fdvSens[], int dvLen);
  void failureStrainSens(const TacsScalar strain[], const TacsScalar weights[],
                         int nfail, TacsScalar failStrainSens[]);

  // Compute the smeared stiffness and its derivative
  // ------------------------------------------------
  void computeStiffness(TacsScalar A[], TacsScalar B[], TacsScalar D[],
                        TacsScalar As[]);
  // void computeStiffnessDVSens( int dvNum, TacsScalar sA[], TacsScalar sB[],
  //                              TacsScalar sD[], TacsScalar sAs[] );

  // Update the buckling loads or frequencies
  // ----------------------------------------
  int computeFrequencies(TacsScalar freq[], int nfreq,
                         const char *prefix = NULL);
  // Compute just the positive values of the buckling loads
  int computeBucklingLoads(TacsScalar Nx, TacsScalar Nxy, TacsScalar loads[],
                           int nloads, const char *prefix = NULL);
  // Compute the positive and negative loads
  int computeBucklingLoads(TacsScalar Nxy, TacsScalar posLoads[],
                           TacsScalar negLoads[], int nloads);

  // Compute the solution with a uniform surface pressure
  // ----------------------------------------------------
  int computePressureLoad(TacsScalar p, const char *file_name);

  // Compute the design variable sensitivity of the
  // buckling loads or frequencies
  // ----------------------------------------------
  int computeFrequenciesDVSens(TacsScalar freq[], TacsScalar freqDVSens[],
                               int nfreq);
  /*
  int computeBucklingLoadsDVSens( TacsScalar Nx, TacsScalar Nxy,
                                  TacsScalar loads[],
                                  TacsScalar loadDVSens[], int nloads );
  int computeBucklingLoadsDVSens( TacsScalar Nxy,
                                  TacsScalar posLoads[], TacsScalar negLoads[],
                                  TacsScalar posLoadDVSens[],
                                  TacsScalar negLoadDVSens[], int nloads );
  */

  // Set the geometric design variable dependence
  // --------------------------------------------
  void setGeoDesignDependence(TacsScalar *_XptConst, TacsScalar *_XptLin,
                              TacsScalar *_geoDvs, int *_dvNums, int _nDvGeo);
  void setGeoDVBounds(TacsScalar *lb, TacsScalar *ub);
  void setLxDvBounds(int _LxDvNum, TacsScalar _LxLb, TacsScalar _LxUb);

  // Design variable information
  // ---------------------------
  void setDesignVars(const TacsScalar dvs[], int numDVs);
  void getDesignVars(TacsScalar dvs[], int numDVs) const;
  void getDesignVarRange(TacsScalar lb[], TacsScalar ub[], int numDVs) const;

  // Test the implementations of the sensitivities
  // ---------------------------------------------
  /*
  void testStiffnessDVSens( double dh );
  void testFailureSens( double dh,
                        const TacsScalar strain[],
                        const TacsScalar weights[],
                        int nfail );
  void testBucklingDVSens( double dh,
                           TacsScalar Nx, TacsScalar Nxy,
                           int nloads );
  */
  // Set the flag to indicate whether to use LAPACK or not
  // -----------------------------------------------------
  void setUseLapackEigensolver(int use_lapack);
  void setLanczosSubspaceSize(int subspace_size);

  // Retrieve the problem size - for load balancing purposes
  // -------------------------------------------------------
  int getProblemSize() { return nvars; }

 private:
  // Compute the raw eigenvalues and eigenmodes
  // ------------------------------------------
  int computeEigenvalues(int lm, const char *auplo, const char *buplo, int n,
                         int ka, TacsScalar *A, int lda, int kb, TacsScalar *B,
                         int ldb, int k, int m, TacsScalar *work, double tol,
                         TacsScalar *eigs, TacsScalar *Z);
  int computeFrequencies(int neigs, TacsScalar eigvals[], TacsScalar eigvecs[]);
  int computeBucklingLoads(const TacsScalar segmentLoads[],
                           const TacsScalar beamLoads[], int neigs,
                           TacsScalar eigvals[], TacsScalar eigvecs[],
                           int two_sided);

  // Function to compute the derivative of the buckling loads
  // --------------------------------------------------------
  int computeBucklingLoadsDVSens(TacsScalar Nx, TacsScalar Nxy,
                                 TacsScalar posLoads[], TacsScalar negLoads[],
                                 TacsScalar posLoadDVSens[],
                                 TacsScalar negLoadDVSens[], int nloads);

  // Compute the segment loads and their sensitivities
  // -------------------------------------------------
  void computeSegmentLoads(TacsScalar Nx, TacsScalar Nxy,
                           TacsScalar segmentLoads[], TacsScalar beamLoads[]);
  void computeSegmentLoadsDVSens(int dvNum, TacsScalar Nx, TacsScalar Nxy,
                                 TacsScalar segmentLoads[],
                                 TacsScalar beamLoads[]);
  void computeSegmentLoadsGeoSens(TacsScalar Nx, TacsScalar Nxy, int dv,
                                  TacsScalar segmentLoads[],
                                  TacsScalar beamLoads[]);

  // Print the displaced shape of the panel
  // --------------------------------------
  void printPanelMode(const char *file_name, const TacsScalar *x, int nx);

  // Update the Xpt array based on the current values of the geoDvs
  void updateGeometry();

  // The maximum number of nodes in an element
  static const int NUM_NODES = 2;

  // Add a segment to the stiffness matrix
  // -------------------------------------
  void addStiffMat(TacsScalar mat[], int n1, int n2, const TacsScalar At[],
                   const TacsScalar Bt[], const TacsScalar Dt[]);

  void addStiffMatBeam(TacsScalar mat[], int n1, const TacsScalar Ct[]);

  // Add a segment to the mass matrix
  // --------------------------------
  void addMassMat(TacsScalar mat[], int n1, int n2, const TacsScalar mass[]);
  void addMassMatBeam(TacsScalar mat[], int n1, const TacsScalar mass[]);

  // Add a segment to the geometric stiffness matrix
  // -----------------------------------------------
  void addGeoStiffMat(TacsScalar mat[], const TacsScalar s[], int n1, int n2);

  // Add the material-dependent sensitivity from a segment
  // -----------------------------------------------------
  void addMassMatDVSens(TacsScalar smat[], const TacsScalar Z[], int ldz,
                        int nz, int n1, int n2, const TacsScalar mass[]);
  void addStiffMatDVSens(TacsScalar smat[], const TacsScalar Z[], int ldz,
                         int nz, int n1, int n2, const TacsScalar sAt[],
                         const TacsScalar sBt[], const TacsScalar sDt[]);
  void addGeoStiffMatDVSens(TacsScalar smat[], const TacsScalar Z[], int ldz,
                            int nz, int n1, int n2, const TacsScalar sstress[]);

  // Add the geometric-dependent sensitivity from a segment
  // ------------------------------------------------------
  void addMassMatGeoSens(TacsScalar smat[], int dv, const TacsScalar Z[],
                         int ldz, int nz, int n1, int n2,
                         const TacsScalar mass[]);
  void addStiffMatGeoSens(TacsScalar smat[], int dv, const TacsScalar Z[],
                          int ldz, int nz, int n1, int n2,
                          const TacsScalar At[], const TacsScalar Bt[],
                          const TacsScalar Dt[]);
  void addGeoStiffMatGeoSens(TacsScalar smat[], int dv, const TacsScalar Z[],
                             int ldz, int nz, int n1, int n2,
                             const TacsScalar stress[],
                             const TacsScalar sstress[]);

  // Add the sensitivity from the Lx-dependent terms
  // -----------------------------------------------
  void addMassMatLxSens(TacsScalar smat[], const TacsScalar Z[], int ldz,
                        int nz, int n1, int n2, const TacsScalar mass[]);
  void addStiffMatLxSens(TacsScalar smat[], const TacsScalar Z[], int ldz,
                         int nz, int n1, int n2, const TacsScalar At[],
                         const TacsScalar Bt[], const TacsScalar Dt[]);
  void addGeoStiffMatLxSens(TacsScalar smat[], const TacsScalar Z[], int ldz,
                            int nz, int n1, int n2, const TacsScalar stress[]);

  // Add pressure loads to a global force vector
  // -------------------------------------------
  void addPressureLoad(TacsScalar F[], int n1, int n2, TacsScalar p);

  void addValues(TacsScalar mat[], int kband, int n, const int rows[], int m,
                 const int cols[], const TacsScalar values[]);
  void addValuesTranspose(TacsScalar mat[], int kband, int n, const int rows[],
                          int m, const int cols[], const TacsScalar values[]);

  void addSensMatValues(TacsScalar smat[], const TacsScalar Z[], int ldz,
                        int nz, int n, const int rows[], int m,
                        const int cols[], const TacsScalar values[]);
  void addSensMatValuesTranspose(TacsScalar smat[], const TacsScalar Z[],
                                 int ldz, int nz, int n, const int rows[],
                                 int m, const int cols[],
                                 const TacsScalar values[]);

  // Transform a finite-strip element to the global coordinate frame
  // ---------------------------------------------------------------
  void transformVec(TacsScalar vec[], const TacsScalar t[], int nn);
  void transformMat(TacsScalar mat[], const TacsScalar t[], int nn, int mm);
  void transformMatGeoSens(TacsScalar mat[], TacsScalar smat[],
                           const TacsScalar t[], const TacsScalar st[], int nn,
                           int mm);

  // Code to compute the derivative of the strain w.r.t. the nodal displacements
  // ---------------------------------------------------------------------------
  void computeBsin(TacsScalar B[], TacsScalar lambda, TacsScalar h,
                   const double N[], const double Na[], const double Nhp[],
                   const double Nahp[], const double Naahp[]);
  void computeBcos(TacsScalar B[], TacsScalar lambda, TacsScalar h,
                   const double N[], const double Na[], const double Nhp[],
                   const double Nahp[], const double Naahp[]);

  void computeBsinGeoSens(TacsScalar B[], TacsScalar lambda, TacsScalar slambda,
                          TacsScalar h, TacsScalar sh, const double N[],
                          const double Na[], const double Nhp[],
                          const double Nahp[], const double Naahp[]);
  void computeBcosGeoSens(TacsScalar B[], TacsScalar lambda, TacsScalar slambda,
                          TacsScalar h, TacsScalar sh, const double N[],
                          const double Na[], const double Nhp[],
                          const double Nahp[], const double Naahp[]);

  // Code to add the contributions to the geometric stiffness matrix
  // ---------------------------------------------------------------
  void addGs(TacsScalar G[], TacsScalar scale, const TacsScalar s[],
             TacsScalar lambda, TacsScalar h, const double N[],
             const double Na[], const double Nhp[], const double Nahp[]);
  void addGcs(TacsScalar G[], TacsScalar scale, const TacsScalar s[],
              TacsScalar lambda_n, TacsScalar lambda_m, TacsScalar h,
              const double N[], const double Na[], const double Nhp[],
              const double Nahp[]);

  void addGsGeoSens(TacsScalar G[], TacsScalar scale, TacsScalar sscale,
                    const TacsScalar s[], TacsScalar lambda, TacsScalar slambda,
                    TacsScalar h, TacsScalar sh, const double N[],
                    const double Na[], const double Nhp[], const double Nahp[]);
  void addGcsGeoSens(TacsScalar G[], TacsScalar scale, TacsScalar sscale,
                     const TacsScalar s[], TacsScalar lambda_n,
                     TacsScalar slambda_n, TacsScalar lambda_m,
                     TacsScalar slambda_m, TacsScalar h, TacsScalar sh,
                     const double N[], const double Na[], const double Nhp[],
                     const double Nahp[]);

  // Compute inertial terms
  // ----------------------
  void computeMsin(TacsScalar M[], TacsScalar lambda, TacsScalar h,
                   const double N[], const double Nhp[], const double Nahp[]);
  void computeMcos(TacsScalar M[], TacsScalar lambda, TacsScalar h,
                   const double N[], const double Nhp[], const double Nahp[]);

  void computeMsinGeoSens(TacsScalar M[], TacsScalar lambda, TacsScalar slambda,
                          TacsScalar h, TacsScalar sh, const double N[],
                          const double Nhp[], const double Nahp[]);
  void computeMcosGeoSens(TacsScalar M[], TacsScalar lambda, TacsScalar slambda,
                          TacsScalar h, TacsScalar sh, const double N[],
                          const double Nhp[], const double Nahp[]);

  // Compute the derivative of the strain for a longitudinal beam
  // ------------------------------------------------------------
  void computeBeamBsin(TacsScalar B[], TacsScalar lambda);
  void computeBeamBcos(TacsScalar B[], TacsScalar lambda);

  void computeBeamMsin(TacsScalar M[], TacsScalar lambda);
  void computeBeamMcos(TacsScalar M[], TacsScalar lambda);

  // Inline inner-product functions
  // ------------------------------
  TacsScalar innerProduct(const TacsScalar s[], const TacsScalar e[]) {
    return (s[0] * e[0] + s[1] * e[1] + s[2] * e[2] + s[3] * e[3] +
            s[4] * e[4] + s[5] * e[5]);
  }

  TacsScalar innerBeamProduct(const TacsScalar s[], const TacsScalar e[]) {
    return (s[0] * e[0] + s[1] * e[1] + s[2] * e[2] + s[3] * e[3]);
  }

  // Compute the inertial terms for the mass matrix
  // k = d/dt[u, v, w, dw/dx, dv/dx]
  // m = int_{t} [ rho, rho*z, rho*z^2 ] dz
  TacsScalar innerMassProduct(const TacsScalar m[], const TacsScalar k0[],
                              const TacsScalar k1[]) {
    return (m[0] * (k0[0] * k1[0] + k0[1] * k1[1] + k0[2] * k1[2]) +
            m[1] * (k0[3] * k1[3] + k0[4] * k1[4]));
  }

  // Compute the inertial terms for the mass matrix
  // k = d/dt[u, v, w, dw/dx, dw/dy]
  // m = int_{A} [rho, rho*z, rho*y, rho*z^2, rho*y*z, rho*y^2 ] dA
  // where z = strong axis, y = weak axis
  TacsScalar innerBeamMassProduct(const TacsScalar m[], const TacsScalar k0[],
                                  const TacsScalar k1[]) {
    return (m[0] * (k0[0] * k1[0] + k0[1] * k1[1] + k0[2] * k1[2]) -
            m[1] * (k0[0] * k1[3] + k0[3] * k1[0]) -
            m[2] * (k0[1] * k1[4] + k0[4] * k1[1]) + m[3] * k0[3] * k1[3] +
            m[5] * k0[4] * k1[4] - m[4] * (k0[3] * k1[4] + k0[4] * k1[3]));
  }

  // Compute the stress in a segment
  // -------------------------------
  void computeStress(TacsScalar s[], const TacsScalar e[],
                     const TacsScalar At[], const TacsScalar Bt[],
                     const TacsScalar Dt[]) {
    s[0] = At[0] * e[0] + At[1] * e[1] + At[2] * e[2] + Bt[0] * e[3] +
           Bt[1] * e[4] + Bt[2] * e[5];
    s[1] = At[1] * e[0] + At[3] * e[1] + At[4] * e[2] + Bt[1] * e[3] +
           Bt[3] * e[4] + Bt[4] * e[5];
    s[2] = At[2] * e[0] + At[4] * e[1] + At[5] * e[2] + Bt[2] * e[3] +
           Bt[4] * e[4] + Bt[5] * e[5];

    s[3] = Bt[0] * e[0] + Bt[1] * e[1] + Bt[2] * e[2] + Dt[0] * e[3] +
           Dt[1] * e[4] + Dt[2] * e[5];
    s[4] = Bt[1] * e[0] + Bt[3] * e[1] + Bt[4] * e[2] + Dt[1] * e[3] +
           Dt[3] * e[4] + Dt[4] * e[5];
    s[5] = Bt[2] * e[0] + Bt[4] * e[1] + Bt[5] * e[2] + Dt[2] * e[3] +
           Dt[4] * e[4] + Dt[5] * e[5];
  }

  void computeBeamStress(TacsScalar s[], const TacsScalar e[],
                         const TacsScalar Ct[]) {
    s[0] = Ct[0] * e[0] + Ct[1] * e[1] + Ct[2] * e[2] + Ct[3] * e[3];
    s[1] = Ct[1] * e[0] + Ct[4] * e[1] + Ct[5] * e[2] + Ct[6] * e[3];
    s[2] = Ct[2] * e[0] + Ct[5] * e[1] + Ct[7] * e[2] + Ct[8] * e[3];
    s[3] = Ct[3] * e[0] + Ct[6] * e[1] + Ct[8] * e[2] + Ct[9] * e[3];
  }

  // The Gauss points and weights
  int numGauss;
  const double *gaussWts;
  const double *gaussPts;

  // Set a flag to indicate whether to use the LAPACK eigenvalue
  // computational routines, or to use the internal eigenvalue
  // solver (default - it's faster b/c LAPACK computes the full
  // spectrum)
  int use_lapack_eigensolver;

  // Set the default Lanczos subspace size
  int lanczos_subspace_size;
  double lanczos_eigen_tol;

  // Assign the start and end nodes - these are used
  // to compute the mass per unit area of the panel
  int first_node, last_node;
  int first_node_bc, last_node_bc;

  int nDvGeo;          // Number of geometric design variables
  int *geoDvNums;      // Geometric design variable numbers
  TacsScalar *geoDvs;  // The geometric design variable values
  TacsScalar *geoLb, *geoUb;
  TacsScalar *XptConst, *XptLin;

  TacsScalar beta;  // Everything is computed along the diagonal x + beta*y

  int LxDvNum;    // The design variable number associated with the panel
  TacsScalar Lx;  // The length of the panel
  TacsScalar LxLb, LxUb;  // The lower and upper bounds associated with Lx

  TacsScalar *Xpts;  // The positions of the nodes in the panel
  int nnodes;        // The number of nodal locations
  int *nodes;        // The nodes for each segment
  int *bnodes;       // The nodes associated with longitudinal beams

  FSDTStiffness **panels;  // The sub-stiffness properties of each panel
  EBStiffness **beams;     // The beam stiffness objects
  int nsegments;           // The number of panel segments
  int nbeams;              // The number of longitudinal beams
  int nmodes;              // The number of terms in the
  int *segmentType;        // Segment types - skin or stiffener

  // The unknowns associated with each panel segment
  // Note that negative unknowns do not appear in the final matrix
  int *vars;
  int nvars;  // The number of variables
  int nband;  // The number of super-diagonals in the banded matrix

  // The points in the section that are checked for a failure
  static const int MAX_NUM_FAIL_POINTS = 10;
  int numFailPoints;
  int failSegments[MAX_NUM_FAIL_POINTS], failNodes[MAX_NUM_FAIL_POINTS];
  int failPointIsSkin[MAX_NUM_FAIL_POINTS];
};

#endif  // TACS_USE_COMPLEX
#endif  // TACS_PANEL_ANALYSIS_H
