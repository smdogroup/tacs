#ifndef TACS_RBE2_H
#define TACS_RBE2_H

/*
  The RBE2 element.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Copyright (C) 2020 Aerion Technologies Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.
*/

#include "TACSElement.h"

class TACSRBE2 : public TACSElement {
 public:
  TACSRBE2(int _numNodes, int _dof_constrained[], double _C1 = 1e3,
           double _C2 = 1e-3);
  ~TACSRBE2();

  // Info for BDF writer
  // -------------------
  int const* const* getDependentDOFs() { return dof_constrained; }
  int getNumDependentNodes() { return NUM_DEP_NODES; }

  // Get the element properties and names
  // ------------------------------------
  const char* getObjectName();
  const char* displacementName(int i);
  const char* stressName(int i);
  const char* strainName(int i);
  const char* extraName(int i);
  int getVarsPerNode();
  int getNumNodes();
  int numExtras();
  ElementType getElementType();
  int getDesignVarsPerNode() { return 0; }
  int getNumQuadraturePoints() { return 0; }
  double getQuadratureWeight(int n) { return 0.0; }
  double getQuadraturePoint(int n, double pt[]) { return 0.0; }
  int getNumElementFaces() { return 0; }
  int getNumFaceQuadraturePoints(int face) { return 0; }
  double getFaceQuadraturePoint(int face, int n, double pt[],
                                double tangent[]) {
    return 0.0;
  }
  /*int getMultiplierIndex();*/

  // Functions for analysis
  // ----------------------
  void addResidual(int elemIndex, double time, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[], TacsScalar res[]);
  void addJacobian(int elemIndex, double time, TacsScalar alpha,
                   TacsScalar beta, TacsScalar gamma, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[], TacsScalar res[],
                   TacsScalar mat[]);

  // Functions required to determine the derivatives w.r.t. the design variables
  // ---------------------------------------------------------------------------
  void addAdjResXptProduct(int elemIndex, double time, TacsScalar scale,
                           const TacsScalar psi[], const TacsScalar Xpts[],
                           const TacsScalar vars[], const TacsScalar dvars[],
                           const TacsScalar ddvars[], TacsScalar fXptSens[]);

  // Functions for post-processing
  // -----------------------------
  /*void getOutputData( int elemIndex,
                      ElementType etype, int write_flag,
                      const TacsScalar Xpts[],
                      const TacsScalar vars[],
                      const TacsScalar dvars[],
                      const TacsScalar ddvars[],
                      int ld_data, TacsScalar *data );*/

 private:
  void getMaskedMultipliers(TacsScalar maskedLM[], const TacsScalar actualLM[]);

  static const int NUM_DISPS = 6;
  int NUM_NODES;  // Number of nodes (1 indep node + N dep nodes + N dummy
                  // nodes)
  int NUM_DEP_NODES;
  int NUM_VARIABLES;
  static const int NUM_EXTRAS = 6;

  static const char* elemName;
  static const char* dispNames[NUM_DISPS];
  static const char* extraNames[NUM_EXTRAS];

  // Flag which dependent dofs to include
  int** dof_constrained;

  // constraint matrix scaling factor, see ref [2]
  double C1;
  // artificial stiffness scaling factor, see ref [2]
  double C2;
};

#endif
