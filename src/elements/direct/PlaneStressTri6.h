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

#ifndef TACS_PLANE_STRESS_TRI6_H
#define TACS_PLANE_STRESS_TRI6_H

#include "TACS2DElement.h"
#include "TACSElement.h"

class PlaneStressTri6 : public TACS2DElement<6> {
 public:
  PlaneStressTri6(PlaneStressStiffness *_stiff,
                  ElementBehaviorType type = LINEAR, int _componentNum = 0);
  ~PlaneStressTri6();

  // Return the name of this element
  // -------------------------------
  const char *elementName() { return elemName; }

  // Retrieve the shape functions
  // ----------------------------
  void getShapeFunctions(const double pt[], double N[]);
  void getShapeFunctions(const double pt[], double N[], double Na[],
                         double Nb[]);

  // Retrieve the Gauss points/weights
  // ---------------------------------
  int getNumGaussPts();
  double getGaussWtsPts(const int num, double pt[]);

  // Functions for post-processing
  // -----------------------------
  void addOutputCount(int *nelems, int *nnodes, int *ncsr);
  void getOutputData(unsigned int out_type, double *data, int ld_data,
                     const TacsScalar Xpts[], const TacsScalar vars[]);
  void getOutputConnectivity(int *con, int node);

 private:
  static const int NUM_NODES = 6;
  static const char *elemName;
};

#endif  // TACS_PLANE_STRESS_TRI6_H
