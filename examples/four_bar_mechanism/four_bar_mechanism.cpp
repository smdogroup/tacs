#include "MITC3.h"
#include "TACSAssembler.h"
#include "TACSConstitutiveVerification.h"
#include "TACSElementVerification.h"
#include "TACSIntegrator.h"
#include "TACSKSFailure.h"
#include "TACSKinematicConstraints.h"
#include "TACSRigidBody.h"
#include "TACSStructuralMass.h"
#include "TACSTimoshenkoConstitutive.h"

class SquareSection : public TACSTimoshenkoConstitutive {
 public:
  static const double kcorr;

  SquareSection(TacsScalar _density, TacsScalar _E, TacsScalar _G,
                TacsScalar _w, int _wNum, const TacsScalar axis[])
      : TACSTimoshenkoConstitutive(NULL, NULL, axis) {
    density = _density;
    E = _E;
    G = _G;
    w = _w;
    wNum = _wNum;
    if (wNum < 0) {
      wNum = 0;
    }

    computeProperties();
  }

  void computeProperties() {
    // Set the properties based on the width/thickness variables
    TacsScalar A = w * w;
    TacsScalar Iy = w * w * w * w / 12.0;
    TacsScalar Iz = Iy;
    TacsScalar J = Iy + Iz;
    TacsScalar Iyz = 0.0;

    // Set the entries of the stiffness matrix
    memset(C, 0, 36 * sizeof(TacsScalar));
    C[0] = E * A;
    C[7] = G * J;
    C[14] = E * Iy;
    C[21] = E * Iz;
    C[28] = kcorr * G * A;
    C[35] = kcorr * G * A;

    // Set the entries of the density matrix
    rho[0] = density * A;
    rho[1] = density * Iy;
    rho[2] = density * Iz;
    rho[3] = density * Iyz;
  }

  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]) {
    if (dvNums) {
      dvNums[0] = wNum;
    }
    return 1;
  }
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]) {
    w = dvs[0];
    computeProperties();
    return 1;
  }
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]) {
    dvs[0] = w;
    return 1;
  }
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]) {
    lb[0] = 0.0;
    ub[0] = 10.0;
    return 1;
  }
  void addStressDVSens(int elemIndex, TacsScalar scale, const double pt[],
                       const TacsScalar X[], const TacsScalar e[],
                       const TacsScalar psi[], int dvLen, TacsScalar dfdx[]) {
    dfdx[0] +=
        scale *
        (2.0 * w *
             (E * e[0] * psi[0] + kcorr * G * (e[4] * psi[4] + e[5] * psi[5])) +
         (w * w * w / 3.0) *
             (2.0 * G * e[1] * psi[1] + E * (e[2] * psi[2] + e[3] * psi[3])));
  }
  void addMassMomentsDVSens(int elemIndex, TacsScalar scale, const double pt[],
                            const TacsScalar psi[], int dvLen,
                            TacsScalar dfdx[]) {
    dfdx[0] += scale * density *
               (2.0 * w * psi[0] + ((w * w * w) / 3.0) * (psi[1] + psi[2]));
  }

  TacsScalar evalDensity(int elemIndex, const double pt[],
                         const TacsScalar X[]) {
    return density * w * w;
  }
  void addDensityDVSens(int elemIndex, const TacsScalar scale,
                        const double pt[], const TacsScalar X[], int dvLen,
                        TacsScalar dfdx[]) {
    dfdx[0] += 2.0 * scale * density * w;
  }

  TacsScalar evalFailure(int elemIndex, const double pt[], const TacsScalar X[],
                         const TacsScalar e[]) {
    return E * w * w * fabs(e[0]) / 10e3;
  }
  TacsScalar evalFailureStrainSens(int elemIndex, const double pt[],
                                   const TacsScalar X[], const TacsScalar e[],
                                   TacsScalar sens[]) {
    memset(sens, 0, 6 * sizeof(TacsScalar));
    if (TacsRealPart(e[0]) >= 0.0) {
      sens[0] = E * w * w / 10e3;
    } else {
      sens[0] = -E * w * w / 10e3;
    }
    return E * w * w * fabs(e[0]) / 10e3;
  }
  void addFailureDVSens(int elemIndex, TacsScalar scale, const double pt[],
                        const TacsScalar X[], const TacsScalar e[], int dvLen,
                        TacsScalar dfdx[]) {
    dfdx[0] += 2.0 * scale * E * w * fabs(e[0]) / 10e3;
  }

  TacsScalar density, E, G;
  TacsScalar w;
  int wNum;
};

const double SquareSection::kcorr = 5.0 / 6.0;

/*
  Create and return the TACSAssembler object for the four bar
  mechanism as described by Bachau

  B ------------------- C
  |                     |
  |                     |
  |                     |
  A                     D

  Length between A and B = 0.12 m
  Length between B and C = 0.24 m
  Length between C and D = 0.12 m

  A, B and D are revolute joints in the plane perpendicular to the
  plane of the mechanism

  C is a revolute joint in a plane +5 degrees along the DC axis of the
  beam

  Beam properties:

  Young's modulus 207 GPa, nu = 0.3

  Bars 1 and 2 are square and of dimension 16 x 16 mm
  Bar 3 is square and of dimension 8 x 8 mm
*/
TACSAssembler *four_bar_mechanism(int nA, int nB, int nC) {
  // Set the gravity vector
  TACSGibbsVector *gravity = new TACSGibbsVector(0.0, 0.0, -9.81);

  // Set the points b, c and d
  TACSGibbsVector *ptB = new TACSGibbsVector(0.0, 0.12, 0.0);
  TACSGibbsVector *ptC = new TACSGibbsVector(0.24, 0.12, 0.0);
  TACSGibbsVector *ptD = new TACSGibbsVector(0.24, 0.0, 0.0);

  // Create the revolute direction for B and D
  TACSGibbsVector *revDirA = new TACSGibbsVector(0.0, 0.0, 1.0);
  TACSGibbsVector *revDirB = new TACSGibbsVector(0.0, 0.0, 1.0);
  TACSGibbsVector *revDirD = new TACSGibbsVector(0.0, 0.0, 1.0);

  // Create the revolute direction for C
  TacsScalar theta = (5.0 * M_PI / 180.0);
  TACSGibbsVector *revDirC = new TACSGibbsVector(sin(theta), 0.0, cos(theta));

  // Create the revolute constraints
  TacsScalar omega = -0.6;  // rad/seconds
  int fixed_point = 1;
  int not_fixed = 0;

  TACSRevoluteDriver *revDriverA = new TACSRevoluteDriver(revDirA, omega);
  TACSRevoluteConstraint *revB =
      new TACSRevoluteConstraint(not_fixed, ptB, revDirB);
  TACSRevoluteConstraint *revC =
      new TACSRevoluteConstraint(not_fixed, ptC, revDirC);
  TACSRevoluteConstraint *revD =
      new TACSRevoluteConstraint(fixed_point, ptD, revDirD);

  // Set the reference axes for each beam
  TacsScalar axis_A[] = {-1.0, 0.0, 0.0};
  TacsScalar axis_B[] = {0.0, 1.0, 0.0};
  TacsScalar axis_C[] = {1.0, 0.0, 0.0};

  // Set the material properties
  TacsScalar density = 7800.0;
  TacsScalar E = 207e9;
  TacsScalar nu = 0.3;
  TacsScalar G = 0.5 * E / (1.0 + nu);

  TacsScalar wA = 0.016;
  TacsScalar wB = 0.008;
  int wANum = 0, wBNum = 1;

  TACSTimoshenkoConstitutive *stiffA =
      new SquareSection(density, E, G, wA, wANum, axis_A);

  TACSTimoshenkoConstitutive *stiffB =
      new SquareSection(density, E, G, wA, wANum, axis_B);

  TACSTimoshenkoConstitutive *stiffC =
      new SquareSection(density, E, G, wB, wBNum, axis_C);

  // Set up the connectivity
  MITC3 *beamA = new MITC3(stiffA, gravity);
  MITC3 *beamB = new MITC3(stiffB, gravity);
  MITC3 *beamC = new MITC3(stiffC, gravity);

  // Set the number of nodes in the mesh
  int nnodes = (2 * nA + 1) + (2 * nB + 1) + (2 * nC + 1) + 4;

  // Set the number of elements
  int nelems = nA + nB + nC + 4;

  // Create the connectivities
  TacsScalar *X = new TacsScalar[3 * nnodes];
  memset(X, 0, 3 * nnodes * sizeof(TacsScalar));

  int *ptr = new int[nelems + 1];
  int *conn = new int[3 * nelems];
  TACSElement **elems = new TACSElement *[nelems];

  // Set the nodes numbers and locations
  int *nodesA = new int[2 * nA + 1];
  int *nodesB = new int[2 * nB + 1];
  int *nodesC = new int[2 * nC + 1];
  int n = 0;
  for (int i = 0; i < 2 * nA + 1; i++, n++) {
    nodesA[i] = n;
    X[3 * n + 1] = 0.12 * i / (2 * nA);
  }
  for (int i = 0; i < 2 * nB + 1; i++, n++) {
    nodesB[i] = n;
    X[3 * n] = 0.24 * i / (2 * nB);
    X[3 * n + 1] = 0.12;
  }
  for (int i = 0; i < 2 * nC + 1; i++, n++) {
    nodesC[i] = n;
    X[3 * n] = 0.24;
    X[3 * n + 1] = 0.12 * (1.0 - 1.0 * i / (2 * nC));
  }

  // Set the connectivity for the beams
  int elem = 0;
  ptr[0] = 0;
  for (int i = 0; i < nA; i++) {
    conn[ptr[elem]] = nodesA[2 * i];
    conn[ptr[elem] + 1] = nodesA[2 * i + 1];
    conn[ptr[elem] + 2] = nodesA[2 * i + 2];
    elems[elem] = beamA;
    ptr[elem + 1] = ptr[elem] + 3;
    elem++;
  }

  for (int i = 0; i < nB; i++) {
    conn[ptr[elem]] = nodesB[2 * i];
    conn[ptr[elem] + 1] = nodesB[2 * i + 1];
    conn[ptr[elem] + 2] = nodesB[2 * i + 2];
    elems[elem] = beamB;
    ptr[elem + 1] = ptr[elem] + 3;
    elem++;
  }

  for (int i = 0; i < nC; i++) {
    conn[ptr[elem]] = nodesC[2 * i];
    conn[ptr[elem] + 1] = nodesC[2 * i + 1];
    conn[ptr[elem] + 2] = nodesC[2 * i + 2];
    elems[elem] = beamC;
    ptr[elem + 1] = ptr[elem] + 3;
    elem++;
  }

  // Add the connectivities for the constraints
  conn[ptr[elem]] = nodesA[0];
  conn[ptr[elem] + 1] = nnodes - 4;
  elems[elem] = revDriverA;
  ptr[elem + 1] = ptr[elem] + 2;
  elem++;

  conn[ptr[elem]] = nodesA[2 * nA];
  conn[ptr[elem] + 1] = nodesB[0];
  conn[ptr[elem] + 2] = nnodes - 3;
  elems[elem] = revB;
  ptr[elem + 1] = ptr[elem] + 3;
  elem++;

  conn[ptr[elem]] = nodesC[0];
  conn[ptr[elem] + 1] = nodesB[2 * nB];
  conn[ptr[elem] + 2] = nnodes - 2;
  elems[elem] = revC;
  ptr[elem + 1] = ptr[elem] + 3;
  elem++;

  conn[ptr[elem]] = nodesC[2 * nC];
  conn[ptr[elem] + 1] = nnodes - 1;
  elems[elem] = revD;
  ptr[elem + 1] = ptr[elem] + 2;
  elem++;

  delete[] nodesA;
  delete[] nodesB;
  delete[] nodesC;

  // Create the TACSAssembler object
  TACSAssembler *assembler =
      new TACSAssembler(MPI_COMM_WORLD, 8, nnodes, nelems);

  assembler->setElementConnectivity(ptr, conn);
  delete[] conn;
  delete[] ptr;

  assembler->setElements(elems);
  delete[] elems;

  assembler->initialize();

  // Set the node locations
  TACSBVec *Xvec = assembler->createNodeVec();
  Xvec->incref();
  TacsScalar *Xarray;
  Xvec->getArray(&Xarray);
  memcpy(Xarray, X, 3 * nnodes * sizeof(TacsScalar));
  assembler->setNodes(Xvec);
  Xvec->decref();
  delete[] X;

  return assembler;
}

int main(int argc, char *argv[]) {
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Create the finite-element model
  int nA = 4, nB = 8, nC = 4;
  TACSAssembler *assembler = four_bar_mechanism(nA, nB, nC);
  assembler->incref();

  // Set the final time
  double tf = 12.0;

  // The number of total steps (100 per second)
  int num_steps = 1200;

  // Create the integrator class
  TACSIntegrator *integrator =
      new TACSBDFIntegrator(assembler, 0.0, tf, num_steps, 2);
  integrator->incref();

  // Set the integrator options
  integrator->setUseSchurMat(1, TACSAssembler::TACS_AMD_ORDER);
  integrator->setAbsTol(1e-7);
  integrator->setOutputFrequency(10);

  // Integrate the equations of motion forward in time
  integrator->integrate();

  // Create the continuous KS function
  double ksRho = 100.0;
  TACSKSFailure *ksfunc = new TACSKSFailure(assembler, ksRho);
  ksfunc->setKSFailureType(TACSKSFailure::DISCRETE);

  // Set the functions
  TACSFunction *funcs = ksfunc;
  integrator->setFunctions(1, &funcs);

  TacsScalar fval;
  integrator->evalFunctions(&fval);
  printf("Function value: %15.10e\n", TacsRealPart(fval));

  // Evaluate the adjoint
  integrator->integrateAdjoint();

  // Get the gradient
  TACSBVec *dfdx;
  integrator->getGradient(0, &dfdx);

#ifdef TACS_USE_COMPLEX
  integrator->checkGradients(1e-30);
#else
  integrator->checkGradients(1e-6);
#endif  // TACS_USE_COMPLEX

  // Set the output options/locations
  int elem[3];
  elem[0] = nA / 2;
  elem[1] = nA + nB / 2;
  elem[2] = nA + nB + nC / 2;
  double param[][1] = {{-1.0}, {-1.0}, {0.0}};

  // Extra the data to a file
  for (int pt = 0; pt < 3; pt++) {
    char filename[128];
    sprintf(filename, "mid_beam_%d.dat", pt + 1);
    FILE *fp = fopen(filename, "w");

    fprintf(fp, "Variables = t, u0, v0, w0, quantity\n");

    // Write out data from the beams
    TACSBVec *q = NULL;
    for (int k = 0; k < num_steps + 1; k++) {
      TacsScalar X[3 * 3], vars[8 * 3], dvars[8 * 3], ddvars[8 * 3];
      double time = integrator->getStates(k, &q, NULL, NULL);
      assembler->setVariables(q);
      TACSElement *element =
          assembler->getElement(elem[pt], X, vars, dvars, ddvars);

      TacsScalar quantity;
      element->evalPointQuantity(elem[pt], TACS_FAILURE_INDEX, time, 0,
                                 param[pt], X, vars, dvars, ddvars, &quantity);

      fprintf(fp, "%e  %e %e %e  %e\n", time, TacsRealPart(vars[0]),
              TacsRealPart(vars[1]), TacsRealPart(vars[2]),
              TacsRealPart(quantity));
    }
    fclose(fp);
  }

  integrator->decref();
  assembler->decref();

  MPI_Finalize();
  return 0;
}
