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

#ifndef TACS_QUAD_BASIS_H
#define TACS_QUAD_BASIS_H

#include "TACSElementBasis.h"

/**
  Basis class for a linear quad element
*/
class TACSLinearQuadBasis : public TACSElementBasis {
 public:
  ElementLayout getLayoutType();
  void getVisPoint(int n, double pt[]);
  int getNumNodes();
  int getNumParameters();
  int getNumQuadraturePoints();
  double getQuadratureWeight(int n);
  double getQuadraturePoint(int n, double pt[]);
  int getNumElementFaces();
  int getNumFaceQuadraturePoints(int face);
  double getFaceQuadraturePoint(int face, int n, double pt[], double t[]);
  void computeBasis(const double pt[], double N[]);
  void computeBasisGradient(const double pt[], double N[], double Nxi[]);
};

/**
  Basis class for a quadratic quad element
*/
class TACSQuadraticQuadBasis : public TACSElementBasis {
 public:
  ElementLayout getLayoutType();
  void getVisPoint(int n, double pt[]);
  int getNumNodes();
  int getNumParameters();
  int getNumQuadraturePoints();
  double getQuadratureWeight(int n);
  double getQuadraturePoint(int n, double pt[]);
  int getNumElementFaces();
  int getNumFaceQuadraturePoints(int face);
  double getFaceQuadraturePoint(int face, int n, double pt[], double t[]);
  void computeBasis(const double pt[], double N[]);
  void computeBasisGradient(const double pt[], double N[], double Nxi[]);
};

/**
  Basis class for a cubic quad element
*/
class TACSCubicQuadBasis : public TACSElementBasis {
 public:
  ElementLayout getLayoutType();
  void getVisPoint(int n, double pt[]);
  int getNumNodes();
  int getNumParameters();
  int getNumQuadraturePoints();
  double getQuadratureWeight(int n);
  double getQuadraturePoint(int n, double pt[]);
  int getNumElementFaces();
  int getNumFaceQuadraturePoints(int face);
  double getFaceQuadraturePoint(int face, int n, double pt[], double t[]);
  void interpFields(const int n, const double pt[], const int num_fields,
                    const TacsScalar values[], const int incr,
                    TacsScalar field[]);
  void addInterpFieldsTranspose(const int n, const double pt[], const int incr,
                                const TacsScalar field[], const int num_fields,
                                TacsScalar values[]);
  void interpFieldsGrad(const int n, const double pt[], const int num_fields,
                        const TacsScalar values[], TacsScalar grad[]);
  void addInterpFieldsGradTranspose(int n, const double pt[],
                                    const int num_fields,
                                    const TacsScalar grad[],
                                    TacsScalar values[]);
  void computeBasis(const double pt[], double N[]);
  void computeBasisGradient(const double pt[], double N[], double Nxi[]);
};

/**
  Basis class for a quartic quad element
*/
class TACSQuarticQuadBasis : public TACSElementBasis {
 public:
  TACSQuarticQuadBasis();
  ElementLayout getLayoutType();
  void getVisPoint(int n, double pt[]);
  int getNumNodes();
  int getNumParameters();
  int getNumQuadraturePoints();
  double getQuadratureWeight(int n);
  double getQuadraturePoint(int n, double pt[]);
  int getNumElementFaces();
  int getNumFaceQuadraturePoints(int face);
  double getFaceQuadraturePoint(int face, int n, double pt[], double t[]);
  void interpFields(const int n, const double pt[], const int num_fields,
                    const TacsScalar values[], const int incr,
                    TacsScalar field[]);
  void addInterpFieldsTranspose(const int n, const double pt[], const int incr,
                                const TacsScalar field[], const int num_fields,
                                TacsScalar values[]);
  void interpFieldsGrad(const int n, const double pt[], const int num_fields,
                        const TacsScalar values[], TacsScalar grad[]);
  void addInterpFieldsGradTranspose(int n, const double pt[],
                                    const int num_fields,
                                    const TacsScalar grad[],
                                    TacsScalar values[]);
  void computeBasis(const double pt[], double N[]);
  void computeBasisGradient(const double pt[], double N[], double Nxi[]);

 private:
  static const double cosine_pts[5];
  double Nf[25], Nfxi[25];
};

/**
  Basis class for a quintic quad element
*/
class TACSQuinticQuadBasis : public TACSElementBasis {
 public:
  TACSQuinticQuadBasis();
  ElementLayout getLayoutType();
  void getVisPoint(int n, double pt[]);
  int getNumNodes();
  int getNumParameters();
  int getNumQuadraturePoints();
  double getQuadratureWeight(int n);
  double getQuadraturePoint(int n, double pt[]);
  int getNumElementFaces();
  int getNumFaceQuadraturePoints(int face);
  double getFaceQuadraturePoint(int face, int n, double pt[], double t[]);
  void interpFields(const int n, const double pt[], const int num_fields,
                    const TacsScalar values[], const int incr,
                    TacsScalar field[]);
  void addInterpFieldsTranspose(const int n, const double pt[], const int incr,
                                const TacsScalar field[], const int num_fields,
                                TacsScalar values[]);
  void interpFieldsGrad(const int n, const double pt[], const int num_fields,
                        const TacsScalar values[], TacsScalar grad[]);
  void addInterpFieldsGradTranspose(int n, const double pt[],
                                    const int num_fields,
                                    const TacsScalar grad[],
                                    TacsScalar values[]);
  void computeBasis(const double pt[], double N[]);
  void computeBasisGradient(const double pt[], double N[], double Nxi[]);

 private:
  static const double cosine_pts[6];
  double Nf[36], Nfxi[36];
};

#endif  // TACS_QUAD_BASIS_H
