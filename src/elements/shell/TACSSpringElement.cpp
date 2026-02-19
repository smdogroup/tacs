#include "TACSSpringElement.h"

/*
  A 6 DOF spring element

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Copyright (C) 2020 Aerion Technologies Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.
*/

/*
  Create the spring. Set the reference axis.
*/
TACSSpringElement::TACSSpringElement(
    TACSSpringTransform *_transform,
    TACSGeneralSpringConstitutive *_springStiff) {
  transform = _transform;
  transform->incref();

  springStiff = _springStiff;
  springStiff->incref();
}

TACSSpringElement::~TACSSpringElement() {
  transform->decref();
  springStiff->decref();
}

/*
  Retrieve information about the names of the element
*/
const char *TACSSpringElement::getObjectName() { return elemName; }

/*
  Retrieve the numbers of displacements, nodes, stress and variables
*/
int TACSSpringElement::getVarsPerNode() { return NUM_DISPS; }
int TACSSpringElement::getNumNodes() { return NUM_NODES; }

ElementType TACSSpringElement::getElementType() { return TACS_SPRING_ELEMENT; }

double TACSSpringElement::getQuadraturePoint(int n, double pt[]) {
  pt[0] = pt[1] = pt[2] = 0.0;
  if (n == 1) {
    pt[0] = 1.0;
  }
  return 0.5;
}

/*
  The element name, variable, stress and strain names.
*/
const char *TACSSpringElement::elemName = "TACSSpringElement";

/*
  Transform the global element variables gvars to the set of local
  element variables vars.
*/
void TACSSpringElement::transformVarsGlobalToLocal(const TacsScalar t[],
                                                   const TacsScalar gvars[],
                                                   TacsScalar vars[]) {
  for (int i = 0; i < 2 * NUM_NODES; i++) {
    vars[0] = t[0] * gvars[0] + t[1] * gvars[1] + t[2] * gvars[2];
    vars[1] = t[3] * gvars[0] + t[4] * gvars[1] + t[5] * gvars[2];
    vars[2] = t[6] * gvars[0] + t[7] * gvars[1] + t[8] * gvars[2];

    vars += 3;
    gvars += 3;
  }
}

/*
  Compute the local variables in-place such that:

  local vars = transform * global vars
*/
void TACSSpringElement::transformVarsGlobalToLocal(const TacsScalar t[],
                                                   TacsScalar vars[]) {
  TacsScalar a[3];
  for (int i = 0; i < 2 * NUM_NODES; i++) {
    a[0] = t[0] * vars[0] + t[1] * vars[1] + t[2] * vars[2];
    a[1] = t[3] * vars[0] + t[4] * vars[1] + t[5] * vars[2];
    a[2] = t[6] * vars[0] + t[7] * vars[1] + t[8] * vars[2];

    vars[0] = a[0];
    vars[1] = a[1];
    vars[2] = a[2];

    vars += 3;
  }
}

/*
  Compute the global residual in-place such that:

  global residual = transform^{T} * local residual
*/
void TACSSpringElement::transformResLocalToGlobal(const TacsScalar t[],
                                                  TacsScalar res[]) {
  TacsScalar a[3];
  for (int i = 0; i < 2 * NUM_NODES; i++) {
    a[0] = t[0] * res[0] + t[3] * res[1] + t[6] * res[2];
    a[1] = t[1] * res[0] + t[4] * res[1] + t[7] * res[2];
    a[2] = t[2] * res[0] + t[5] * res[1] + t[8] * res[2];

    res[0] = a[0];
    res[1] = a[1];
    res[2] = a[2];

    res += 3;
  }
}

/*
  Compute the derivative of the residual in place

  global residual sensitivity = transformSens^{T} local residual +
  transform^{T} * local residual sensitivity
*/
void TACSSpringElement::transformResLocalToGlobalSens(
    const TacsScalar t[], const TacsScalar tSens[], const TacsScalar resSens[],
    TacsScalar res[]) {
  TacsScalar a[3];
  for (int i = 0; i < 2 * NUM_NODES; i++) {
    a[0] = tSens[0] * res[0] + tSens[3] * res[1] + tSens[6] * res[2];
    a[1] = tSens[1] * res[0] + tSens[4] * res[1] + tSens[7] * res[2];
    a[2] = tSens[2] * res[0] + tSens[5] * res[1] + tSens[8] * res[2];

    res[0] = a[0] + t[0] * resSens[0] + t[3] * resSens[1] + t[6] * resSens[2];
    res[1] = a[1] + t[1] * resSens[0] + t[4] * resSens[1] + t[7] * resSens[2];
    res[2] = a[2] + t[2] * resSens[0] + t[5] * resSens[1] + t[8] * resSens[2];

    res += 3;
    resSens += 3;
  }
}

/*
  Compute the global stiffness matrix in-place such that:

  global mat = transform^{T} * mat * transform
*/
void TACSSpringElement::transformStiffnessMat(const TacsScalar t[],
                                              TacsScalar mat[]) {
  for (int i = 0; i < 2 * NUM_NODES; i++) {
    for (int j = 0; j < 2 * NUM_NODES; j++) {
      TacsScalar a[9], b[9];

      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
          a[3 * ii + jj] = mat[(3 * i + ii) * NUM_VARIABLES + 3 * j + jj];
        }
      }

      // Compute b = tramsform^{T} * a
      // b = [t[0]  t[3]  t[6]][a[0]  a[1]  a[2]]
      // .   [t[1]  t[4]  t[7]][a[3]  a[4]  a[5]]
      // .   [t[2]  t[5]  t[8]][a[6]  a[7]  a[8]]

      b[0] = t[0] * a[0] + t[3] * a[3] + t[6] * a[6];
      b[1] = t[0] * a[1] + t[3] * a[4] + t[6] * a[7];
      b[2] = t[0] * a[2] + t[3] * a[5] + t[6] * a[8];

      b[3] = t[1] * a[0] + t[4] * a[3] + t[7] * a[6];
      b[4] = t[1] * a[1] + t[4] * a[4] + t[7] * a[7];
      b[5] = t[1] * a[2] + t[4] * a[5] + t[7] * a[8];

      b[6] = t[2] * a[0] + t[5] * a[3] + t[8] * a[6];
      b[7] = t[2] * a[1] + t[5] * a[4] + t[8] * a[7];
      b[8] = t[2] * a[2] + t[5] * a[5] + t[8] * a[8];

      // Compute a = b * transform
      // a = [b[0]  b[1]  b[2]][t[0]  t[1]  t[2]]
      // .   [b[3]  b[4]  b[5]][t[3]  t[4]  t[5]]
      // .   [b[6]  b[7]  b[8]][t[6]  t[7]  t[8]]

      a[0] = b[0] * t[0] + b[1] * t[3] + b[2] * t[6];
      a[3] = b[3] * t[0] + b[4] * t[3] + b[5] * t[6];
      a[6] = b[6] * t[0] + b[7] * t[3] + b[8] * t[6];

      a[1] = b[0] * t[1] + b[1] * t[4] + b[2] * t[7];
      a[4] = b[3] * t[1] + b[4] * t[4] + b[5] * t[7];
      a[7] = b[6] * t[1] + b[7] * t[4] + b[8] * t[7];

      a[2] = b[0] * t[2] + b[1] * t[5] + b[2] * t[8];
      a[5] = b[3] * t[2] + b[4] * t[5] + b[5] * t[8];
      a[8] = b[6] * t[2] + b[7] * t[5] + b[8] * t[8];

      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
          mat[(3 * i + ii) * NUM_VARIABLES + 3 * j + jj] = a[3 * ii + jj];
        }
      }
    }
  }
}

void TACSSpringElement::computeEnergies(int elemIndex, double time,
                                        const TacsScalar Xpts[],
                                        const TacsScalar vars[],
                                        const TacsScalar dvars[],
                                        TacsScalar *Te, TacsScalar *Pe) {
  *Te = 0.0;
  *Pe = 0.0;

  TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];
  TacsScalar elemVars[NUM_VARIABLES];
  TacsScalar t[9];
  double pt[3] = {0.0, 0.0, 0.0};

  transform->computeTransform(Xpts, t);

  transformVarsGlobalToLocal(t, vars, elemVars);

  evalStrain(elemVars, strain);

  springStiff->evalStress(elemIndex, pt, Xpts, strain, stress);

  for (int i = 0; i < NUM_STRESSES; i++) {
    *Pe += 0.5 * strain[i] * stress[i];
  }
}

/*
  Assemble the element residual associated with the given design
  variables and elements.
*/
void TACSSpringElement::addResidual(int elemIndex, double time,
                                    const TacsScalar Xpts[],
                                    const TacsScalar vars[],
                                    const TacsScalar dvars[],
                                    const TacsScalar ddvars[],
                                    TacsScalar res[]) {
  TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];
  TacsScalar elemVars[NUM_VARIABLES], elemRes[NUM_VARIABLES];
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];
  TacsScalar t[9];
  double pt[3] = {0.0, 0.0, 0.0};

  memset(elemRes, 0, NUM_VARIABLES * sizeof(TacsScalar));

  transform->computeTransform(Xpts, t);

  transformVarsGlobalToLocal(t, vars, elemVars);

  evalStrain(elemVars, strain);
  evalBmat(B);

  springStiff->evalStress(elemIndex, pt, Xpts, strain, stress);

  for (int i = 0; i < NUM_NODES; i++) {
    for (int ii = 0; ii < NUM_DISPS; ii++) {
      int row = ii + i * NUM_DISPS;
      elemRes[row] += strain_product(&B[NUM_STRESSES * row], stress);
    }
  }

  transformResLocalToGlobal(t, elemRes);

  for (int i = 0; i < NUM_NODES; i++) {
    for (int ii = 0; ii < NUM_DISPS; ii++) {
      int row = ii + i * NUM_DISPS;
      res[row] += elemRes[row];
    }
  }
}

/*
  Assemble the stiffness matrix for the spring element.
*/
void TACSSpringElement::addJacobian(int elemIndex, double time,
                                    TacsScalar alpha, TacsScalar beta,
                                    TacsScalar gamma, const TacsScalar Xpts[],
                                    const TacsScalar vars[],
                                    const TacsScalar dvars[],
                                    const TacsScalar ddvars[], TacsScalar res[],
                                    TacsScalar mat[]) {
  TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];
  TacsScalar elemVars[NUM_VARIABLES], elemRes[NUM_VARIABLES];
  TacsScalar elemMat[NUM_VARIABLES * NUM_VARIABLES];
  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];
  TacsScalar BiStress[NUM_STRESSES];
  TacsScalar t[9];
  double pt[3] = {0.0, 0.0, 0.0};
  memset(elemRes, 0, NUM_VARIABLES * sizeof(TacsScalar));
  memset(elemMat, 0, NUM_VARIABLES * NUM_VARIABLES * sizeof(TacsScalar));

  transform->computeTransform(Xpts, t);

  transformVarsGlobalToLocal(t, vars, elemVars);

  evalStrain(elemVars, strain);
  evalBmat(B);

  springStiff->evalStress(elemIndex, pt, Xpts, strain, stress);

  for (int i = 0; i < NUM_NODES; i++) {
    for (int ii = 0; ii < NUM_DISPS; ii++) {
      int row = ii + i * NUM_DISPS;

      elemRes[row] += strain_product(&B[NUM_STRESSES * row], stress);

      springStiff->evalStress(elemIndex, pt, Xpts, &B[NUM_STRESSES * row],
                              BiStress);

      for (int j = 0; j < NUM_NODES; j++) {
        for (int jj = 0; jj < NUM_DISPS; jj++) {
          // generate another B vector
          int col = jj + j * NUM_DISPS;

          // The regular element matrix
          elemMat[col + row * NUM_VARIABLES] +=
              strain_product(BiStress, &B[NUM_STRESSES * col]);
        }
      }
    }
  }

  transformStiffnessMat(t, elemMat);
  transformResLocalToGlobal(t, elemRes);

  for (int row = 0; row < NUM_VARIABLES; row++) {
    if (res) {
      res[row] += elemRes[row];
    }
    for (int col = 0; col < NUM_VARIABLES; col++) {
      mat[col + row * NUM_VARIABLES] +=
          alpha * elemMat[col + row * NUM_VARIABLES];
    }
  }
}

void TACSSpringElement::getMatType(ElementMatrixType matType, int elemIndex,
                                   double time, const TacsScalar Xpts[],
                                   const TacsScalar vars[], TacsScalar mat[]) {
  memset(mat, 0, NUM_VARIABLES * NUM_VARIABLES * sizeof(TacsScalar));
  if (matType == TACS_STIFFNESS_MATRIX) {
    TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];
    TacsScalar elemVars[NUM_VARIABLES];
    TacsScalar B[NUM_STRESSES * NUM_VARIABLES];
    TacsScalar BiStress[NUM_STRESSES];
    TacsScalar t[9];
    double pt[3] = {0.0, 0.0, 0.0};
    memset(mat, 0, NUM_VARIABLES * NUM_VARIABLES * sizeof(TacsScalar));

    transform->computeTransform(Xpts, t);

    transformVarsGlobalToLocal(t, vars, elemVars);

    evalStrain(elemVars, strain);
    evalBmat(B);

    springStiff->evalStress(elemIndex, pt, Xpts, strain, stress);

    for (int i = 0; i < NUM_NODES; i++) {
      for (int ii = 0; ii < NUM_DISPS; ii++) {
        int row = ii + i * NUM_DISPS;

        springStiff->evalStress(elemIndex, pt, Xpts, &B[NUM_STRESSES * row],
                                BiStress);

        for (int j = 0; j < NUM_NODES; j++) {
          for (int jj = 0; jj < NUM_DISPS; jj++) {
            // generate another B vector
            int col = jj + j * NUM_DISPS;

            // The regular element matrix
            mat[col + row * NUM_VARIABLES] +=
                strain_product(BiStress, &B[NUM_STRESSES * col]);
          }
        }
      }
    }

    transformStiffnessMat(t, mat);
  }
}

/*
  Get the derivative of the residual w.r.t. the nodal coordinates.
*/
void TACSSpringElement::addAdjResXptProduct(
    int elemIndex, double time, TacsScalar scale, const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar fXptSens[]) {
  TacsScalar elemRes[NUM_VARIABLES], resSens[NUM_VARIABLES];

  TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];
  TacsScalar stressSens[NUM_STRESSES], strainSens[NUM_STRESSES];

  TacsScalar B[NUM_STRESSES * NUM_VARIABLES];

  TacsScalar elemVars[NUM_VARIABLES], elemVarsSens[NUM_VARIABLES];
  TacsScalar t[9], tSens[9];
  double pt[3] = {0.0, 0.0, 0.0};

  transform->computeTransform(Xpts, t);
  transformVarsGlobalToLocal(t, vars, elemVars);

  evalStrain(elemVars, strain);

  springStiff->evalStress(elemIndex, pt, Xpts, strain, stress);
  evalBmat(B);

  for (int k = 0; k < 3 * NUM_NODES; k++) {
    memset(elemRes, 0, NUM_VARIABLES * sizeof(TacsScalar));
    memset(resSens, 0, NUM_VARIABLES * sizeof(TacsScalar));

    transform->computeTransformXptSens(Xpts, k, t, tSens);

    transformVarsGlobalToLocal(tSens, vars, elemVarsSens);

    evalStrain(elemVarsSens, strainSens);

    springStiff->evalStress(elemIndex, pt, Xpts, strainSens, stressSens);

    for (int i = 0; i < NUM_NODES; i++) {
      for (int ii = 0; ii < NUM_DISPS; ii++) {
        int row = ii + i * NUM_DISPS;

        elemRes[row] += strain_product(&B[NUM_STRESSES * row], stress);

        resSens[row] += strain_product(&B[NUM_STRESSES * row], stressSens);
      }
    }

    transformResLocalToGlobalSens(t, tSens, resSens, elemRes);
    for (int i = 0; i < NUM_VARIABLES; i++) {
      fXptSens[k] += scale * psi[i] * elemRes[i];
    }
  }
}

int TACSSpringElement::evalPointQuantity(
    int elemIndex, int quantityType, double time, int n, double pt[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar *detXd, TacsScalar *quantity) {
  *detXd = 1.0;
  if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    if (quantity) {
      TacsScalar Te = 0.0;
      computeEnergies(elemIndex, time, Xpts, vars, dvars, &Te, quantity);
      *quantity *= 2.0;
    }
    return 1;
  } else if (quantityType == TACS_ELEMENT_DISPLACEMENT) {
    if (quantity) {
      quantity[0] = vars[0] * (1 - pt[0]) + vars[6] * (pt[0]);
      quantity[1] = vars[1] * (1 - pt[0]) + vars[7] * (pt[0]);
      quantity[2] = vars[2] * (1 - pt[0]) + vars[8] * (pt[0]);
    }
    return 3;
  }

  return 0;
}

/*
  The linear strain expressions, for a spring this is
  just the difference between the displacements at each end.
*/

void TACSSpringElement::evalStrain(const TacsScalar vars[],
                                   TacsScalar strain[]) {
  const TacsScalar *ua, *ub;
  ua = &vars[0];
  ub = &vars[NUM_DISPS];

  for (int i = 0; i < NUM_DISPS; i++) {
    strain[i] = ua[i] - ub[i];
  }
}

/*
  Evaluate the derivative of the strain w.r.t. the nodal degrees of
  freedom.
*/
void TACSSpringElement::evalBmat(TacsScalar B[]) {
  memset(B, 0, NUM_STRESSES * NUM_VARIABLES * sizeof(TacsScalar));

  for (int i = 0; i < NUM_STRESSES; i++) {
    B[i] = 1.0;
    B += NUM_STRESSES;
  }

  for (int i = 0; i < NUM_STRESSES; i++) {
    B[i] = -1.0;
    B += NUM_STRESSES;
  }
}
