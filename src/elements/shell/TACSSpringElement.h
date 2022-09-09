#ifndef TACS_SPRING_ELEMENT_H
#define TACS_SPRING_ELEMENT_H

/*
  The spring element.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Copyright (C) 2020 Aerion Technologies Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  This element takes as input a spring stiffness  and transform object. This
  object may include the reference direction used to define the local
  y-direction of the spring (x is along spring, z is given by x x y).
  If a reference axis is not defined, the spring axis are assumed tom align
  with the global frame
*/

#include "TACSElement.h"
#include "TACSGeneralSpringConstitutive.h"
#include "TACSSpringElementTransform.h"

class TACSSpringElement : public TACSElement {
 public:
  TACSSpringElement(TACSSpringTransform *_transform,
                    TACSGeneralSpringConstitutive *_springStiff);
  ~TACSSpringElement();

  // Get the element properties and names
  // ------------------------------------
  const char *getObjectName();
  int getVarsPerNode();
  int getNumNodes();
  ElementType getElementType();
  int getDesignVarsPerNode() { return 0; }
  int getNumQuadraturePoints() { return 2; }
  double getQuadratureWeight(int n) { return 0.5; }
  double getQuadraturePoint(int n, double pt[]);
  int getNumElementFaces() { return 0; }
  int getNumFaceQuadraturePoints(int face) { return 0; }
  double getFaceQuadraturePoint(int face, int n, double pt[],
                                double tangent[]) {
    return 0.0;
  }

  // Functions for analysis
  // ----------------------
  void computeEnergies(int elemIndex, double time, const TacsScalar Xpts[],
                       const TacsScalar vars[], const TacsScalar dvars[],
                       TacsScalar *Te, TacsScalar *Pe);

  void addResidual(int elemIndex, double time, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[], TacsScalar res[]);

  void addJacobian(int elemIndex, double time, TacsScalar alpha,
                   TacsScalar beta, TacsScalar gamma, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[], TacsScalar res[],
                   TacsScalar mat[]);

  void getMatType(ElementMatrixType matType, int elemIndex, double time,
                  const TacsScalar Xpts[], const TacsScalar vars[],
                  TacsScalar mat[]);

  // Functions required to determine the derivatives w.r.t. the design variables
  // ---------------------------------------------------------------------------
  void addAdjResXptProduct(int elemIndex, double time, TacsScalar scale,
                           const TacsScalar psi[], const TacsScalar Xpts[],
                           const TacsScalar vars[], const TacsScalar dvars[],
                           const TacsScalar ddvars[], TacsScalar fXptSens[]);

  int evalPointQuantity(int elemIndex, int quantityType, double time, int n,
                        double pt[], const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], TacsScalar *detXd,
                        TacsScalar *quantity);

  void addPointQuantityDVSens(int elemIndex, int quantityType, double time,
                              TacsScalar scale, int n, double pt[],
                              const TacsScalar Xpts[], const TacsScalar vars[],
                              const TacsScalar dvars[],
                              const TacsScalar ddvars[],
                              const TacsScalar dfdq[], int dvLen,
                              TacsScalar dfdx[]) {
    return;
  }

  void addPointQuantitySVSens(int elemIndex, int quantityType, double time,
                              TacsScalar alpha, TacsScalar beta,
                              TacsScalar gamma, int n, double pt[],
                              const TacsScalar Xpts[], const TacsScalar vars[],
                              const TacsScalar dvars[],
                              const TacsScalar ddvars[],
                              const TacsScalar dfdq[], TacsScalar dfdu[]) {
    return;
  }

  void addPointQuantityXptSens(int elemIndex, int quantityType, double time,
                               TacsScalar scale, int n, double pt[],
                               const TacsScalar Xpts[], const TacsScalar vars[],
                               const TacsScalar dvars[],
                               const TacsScalar ddvars[],
                               const TacsScalar dfddetXd,
                               const TacsScalar dfdq[], TacsScalar dfdXpts[]) {
    return;
  }

 private:
  void transformVarsGlobalToLocal(const TacsScalar t[],
                                  const TacsScalar gvars[], TacsScalar vars[]);

  void transformVarsGlobalToLocal(const TacsScalar t[], TacsScalar vars[]);

  void transformResLocalToGlobal(const TacsScalar t[], TacsScalar res[]);

  void transformResLocalToGlobalSens(const TacsScalar t[],
                                     const TacsScalar tSens[],
                                     const TacsScalar resSens[],
                                     TacsScalar res[]);

  void transformStiffnessMat(const TacsScalar t[], TacsScalar mat[]);

  inline TacsScalar strain_product(const TacsScalar a[], const TacsScalar b[]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3] + a[4] * b[4] +
           a[5] * b[5];
  }

  // Function pointers to the strain expressions
  void evalStrain(const TacsScalar vars[], TacsScalar strain[]);
  void evalBmat(TacsScalar B[]);

  static const int NUM_DISPS = 6;
  static const int NUM_NODES = 2;
  static const int NUM_VARIABLES = NUM_NODES * NUM_DISPS;
  static const int NUM_STRESSES = 6;

  static const char *elemName;

  TACSGeneralSpringConstitutive *springStiff;
  TACSSpringTransform *transform;
};

#endif  // TACS_SPRING_ELEMENT_H