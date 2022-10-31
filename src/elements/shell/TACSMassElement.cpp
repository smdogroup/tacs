#include "TACSMassElement.h"

#include "TACSMassCentrifugalForce.h"
#include "TACSMassInertialForce.h"

/*
  A 6 DOF point mass element
*/

/*
  Create the MassElement.
*/
TACSMassElement::TACSMassElement(TACSGeneralMassConstitutive *_con) {
  con = _con;
  con->incref();
}

TACSMassElement::~TACSMassElement() { con->decref(); }

/*
  Retrieve information about the names of the element variables
*/
const char *TACSMassElement::getObjectName() { return elemName; }

/*
  Retrieve the numbers of displacements, nodes, stress and variables
*/
int TACSMassElement::getVarsPerNode() { return NUM_DISPS; }

int TACSMassElement::getNumNodes() { return NUM_NODES; }

ElementType TACSMassElement::getElementType() { return TACS_MASS_ELEMENT; }

/*
  The element name, variable, stress and strain names.
*/
const char *TACSMassElement::elemName = "TACSMassElement";

TACSElement *TACSMassElement::createElementInertialForce(
    const TacsScalar inertiaVec[]) {
  return new TACSMassInertialForce(con, inertiaVec);
}

TACSElement *TACSMassElement::createElementCentrifugalForce(
    const TacsScalar omegaVec[], const TacsScalar rotCenter[]) {
  return new TACSMassCentrifugalForce(con, omegaVec, rotCenter);
}

void TACSMassElement::computeEnergies(int elemIndex, double time,
                                      const TacsScalar Xpts[],
                                      const TacsScalar vars[],
                                      const TacsScalar dvars[], TacsScalar *Te,
                                      TacsScalar *Pe) {
  *Pe = 0.0;
  *Te = 0.0;
  TacsScalar f[NUM_DISPS];
  double pt[3] = {0.0, 0.0, 0.0};
  con->evalInertia(elemIndex, pt, Xpts, dvars, f);
  for (int i = 0; i < NUM_DISPS; i++) {
    *Te += 0.5 * dvars[i] * f[i];
  }
}

/*
  Assemble the element residual associated with the given design
  variables and elements.
*/
void TACSMassElement::addResidual(int elemIndex, double time,
                                  const TacsScalar Xpts[],
                                  const TacsScalar vars[],
                                  const TacsScalar dvars[],
                                  const TacsScalar ddvars[], TacsScalar res[]) {
  TacsScalar f[NUM_DISPS];
  double pt[3] = {0.0, 0.0, 0.0};
  con->evalInertia(elemIndex, pt, Xpts, ddvars, f);
  for (int i = 0; i < NUM_DISPS; i++) {
    res[i] += f[i];
  }
}

/*
  Assemble the stiffness matrix for the mass element.
*/
void TACSMassElement::addJacobian(int elemIndex, double time, TacsScalar alpha,
                                  TacsScalar beta, TacsScalar gamma,
                                  const TacsScalar Xpts[],
                                  const TacsScalar vars[],
                                  const TacsScalar dvars[],
                                  const TacsScalar ddvars[], TacsScalar res[],
                                  TacsScalar J[]) {
  double pt[3] = {0.0, 0.0, 0.0};
  for (int j = 0; j < NUM_DISPS; j++) {
    TacsScalar N[NUM_DISPS], f[NUM_DISPS];
    // Shape functions
    memset(N, 0, NUM_DISPS * sizeof(TacsScalar));
    N[j] = 1.0;
    con->evalInertia(elemIndex, pt, Xpts, N, f);
    for (int i = 0; i < NUM_DISPS; i++) {
      J[j + i * NUM_VARIABLES] += gamma * f[i];
      if (res) {
        res[i] += ddvars[j] * f[i];
      }
    }
  }
}

void TACSMassElement::getMatType(ElementMatrixType matType, int elemIndex,
                                 double time, const TacsScalar Xpts[],
                                 const TacsScalar vars[], TacsScalar mat[]) {
  memset(mat, 0, NUM_VARIABLES * NUM_VARIABLES * sizeof(TacsScalar));
  if (matType == TACS_MASS_MATRIX) {
    double pt[3] = {0.0, 0.0, 0.0};
    for (int j = 0; j < NUM_DISPS; j++) {
      TacsScalar N[NUM_DISPS], f[NUM_DISPS];
      // Shape functions
      memset(N, 0, 6 * sizeof(TacsScalar));
      N[j] = 1.0;
      con->evalInertia(elemIndex, pt, Xpts, N, f);
      for (int i = 0; i < NUM_DISPS; i++) {
        mat[j + i * NUM_VARIABLES] = f[i];
      }
    }
  }
}

int TACSMassElement::evalPointQuantity(
    int elemIndex, int quantityType, double time, int n, double pt[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar *detXd, TacsScalar *quantity) {
  if (detXd) {
    *detXd = 1.0;
  }
  if (quantityType == TACS_ELEMENT_DENSITY) {
    if (quantity) {
      *quantity = con->evalDensity(elemIndex, pt, Xpts);
    }
    return 1;
  } else if (quantityType == TACS_ELEMENT_DISPLACEMENT) {
    if (quantity) {
      quantity[0] = vars[0];
      quantity[1] = vars[1];
      quantity[2] = vars[2];
    }
    return 3;
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    if (quantity) {
      TacsScalar mass = con->evalDensity(elemIndex, pt, Xpts);
      quantity[0] = mass * Xpts[0];
      quantity[1] = mass * Xpts[1];
      quantity[2] = mass * Xpts[2];
    }
    return 3;
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    if (quantity) {
      TacsScalar M[21];
      con->evalMassMatrix(elemIndex, pt, Xpts, M);
      quantity[0] =
          M[0] * (Xpts[1] * Xpts[1] + Xpts[2] * Xpts[2]) + M[15];  // Ixx
      quantity[1] = -M[0] * Xpts[0] * Xpts[1] + M[16];             // Ixy
      quantity[2] = -M[0] * Xpts[0] * Xpts[2] + M[17];             // Ixz
      quantity[3] =
          M[0] * (Xpts[0] * Xpts[0] + Xpts[2] * Xpts[2]) + M[18];  // Iyy
      quantity[4] = -M[0] * Xpts[1] * Xpts[2] + M[19];             // Iyz
      quantity[5] =
          M[0] * (Xpts[0] * Xpts[0] + Xpts[1] * Xpts[1]) + M[20];  // Izz
    }
    return 6;
  }

  return 0;
}

void TACSMassElement::addPointQuantityDVSens(
    int elemIndex, int quantityType, double time, TacsScalar scale, int n,
    double pt[], const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[],
    const TacsScalar dfdq[], int dvLen, TacsScalar dfdx[]) {
  if (quantityType == TACS_ELEMENT_DENSITY) {
    con->addDensityDVSens(elemIndex, scale * dfdq[0], pt, Xpts, dvLen, dfdx);
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    TacsScalar dfdmass = 0.0;
    dfdmass += scale * dfdq[0] * Xpts[0];
    dfdmass += scale * dfdq[1] * Xpts[1];
    dfdmass += scale * dfdq[2] * Xpts[2];
    con->addDensityDVSens(elemIndex, dfdmass, pt, Xpts, dvLen, dfdx);
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    TACSElement::addPointQuantityDVSens(elemIndex, quantityType, time, scale, n,
                                        pt, Xpts, vars, dvars, ddvars, dfdq,
                                        dvLen, dfdx);
    return;
  }
}

void TACSMassElement::addPointQuantityXptSens(
    int elemIndex, int quantityType, double time, TacsScalar scale, int n,
    double pt[], const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[],
    const TacsScalar dfddetXd, const TacsScalar dfdq[], TacsScalar dfdXpts[]) {
  if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    TacsScalar mass = con->evalDensity(elemIndex, pt, Xpts);
    dfdXpts[0] += scale * mass * dfdq[0];
    dfdXpts[1] += scale * mass * dfdq[1];
    dfdXpts[2] += scale * mass * dfdq[2];
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    TacsScalar mass = con->evalDensity(elemIndex, pt, Xpts);
    dfdXpts[0] += scale * mass *
                  (2.0 * Xpts[0] * (dfdq[3] + dfdq[5]) - dfdq[1] * Xpts[1] -
                   dfdq[2] * Xpts[2]);
    dfdXpts[1] += scale * mass *
                  (2.0 * Xpts[1] * (dfdq[0] + dfdq[5]) - dfdq[1] * Xpts[0] -
                   dfdq[4] * Xpts[2]);
    dfdXpts[2] += scale * mass *
                  (2.0 * Xpts[2] * (dfdq[0] + dfdq[3]) - dfdq[2] * Xpts[0] -
                   dfdq[4] * Xpts[1]);
  }
}

void TACSMassElement::addPointQuantitySVSens(
    int elemIndex, int quantityType, double time, TacsScalar alpha,
    TacsScalar beta, TacsScalar gamma, int n, double pt[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], const TacsScalar dfdq[], TacsScalar dfdu[]) {
  if (quantityType == TACS_ELEMENT_DISPLACEMENT) {
    dfdu[0] += alpha * dfdq[0];
    dfdu[1] += alpha * dfdq[1];
    dfdu[2] += alpha * dfdq[2];
  }
  return;
}
