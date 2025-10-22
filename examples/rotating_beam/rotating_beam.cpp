#include "KinematicConstraints.h"
#include "MITC3.h"
#include "RigidBody.h"
#include "TACSAssembler.h"
#include "TACSElementAlgebra.h"
#include "TACSIntegrator.h"
#include "TimoshenkoStiffness.h"

int compare(const void *a, const void *b) {
  if (*(TacsScalar *)a < *(TacsScalar *)b) {
    return -1;
  }
  return 1;
}

int main(int argc, char *argv[]) {
  // Initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  const TacsScalar Omega_ref = 109.12;  // rad/s
  double Omega = 1.0;
  int gyroscopic = 1;
  int ne = 20;
  int test_case = 2;
  for (int i = 0; i < argc; i++) {
    if (sscanf(argv[i], "Omega=%lf", &Omega) == 1) {
      printf("Omega = %e\n", Omega * Omega_ref);
    }
    if (sscanf(argv[i], "gyroscopic=%d", &gyroscopic) == 1) {
      printf("gyroscopic = %d\n", gyroscopic);
    }
    if (sscanf(argv[i], "ne=%d", &ne) == 1) {
      printf("ne = %d\n", ne);
    }
  }

  // Set the length of the rotating beam
  TacsScalar chord = 0.121;                 //
  TacsScalar precone = 2.5 * M_PI / 180.0;  // 2.5 degrees
  TacsScalar r = 0.22;
  TacsScalar L = 2.0;

  // Create the Timoshenko stiffness object
  TimoshenkoStiffness *stiff = NULL;

  const TacsScalar angular_rate = Omega * Omega_ref;
  TacsScalar freq_normalization = 1.0;

  if (test_case == 0) {
    // h = sqrt((12*L^2)/beta^2)
    const TacsScalar h = 0.13856406460551018;

    // Set the material properties
    const TacsScalar E = 70e9;      // Pa
    const TacsScalar rho = 2600.0;  // kg/m^3
    const TacsScalar nu = 0.33;     // Poisson's ratio
    const TacsScalar Gmod = 0.5 * E / (1.0 + nu);
    const TacsScalar kappa = 0.5 * (1.0 + nu);  // 5.0/6.0;

    // Set the stiffness properties
    TacsScalar mA = rho * h * h;
    TacsScalar IA = rho * h * h * h * h / 12.0;

    TacsScalar EA = E * h * h;
    TacsScalar GJ = 2.25 * Gmod * h * h * h * h;
    TacsScalar kGAz = kappa * Gmod * h * h;
    TacsScalar EIz = E * h * h * h * h / 12.0;

    // Print out the properties to the screen
    TacsScalar delta = r / L;
    TacsScalar beta = sqrt(mA * L * L / IA);
    TacsScalar eta = EIz / (mA * L * L * L * L);
    TacsScalar T = 1.0 / sqrt(eta);
    TacsScalar gamma = T * Omega;

    printf("delta: %25.15e\n", delta);
    printf("beta:  %25.15e\n", beta);
    printf("eta:   %25.15e\n", eta);
    printf("T:     %25.15e\n", T);
    printf("gamma: %25.15e\n", gamma);
    printf("kG/E:  %25.15e\n", kGAz / EA);

    // Set the reference axes
    TacsScalar axis_A[] = {0.0, 1.0, 0.0};

    stiff = new TimoshenkoStiffness(mA, IA, IA, 0.0, EA, GJ, EIz, EIz, kGAz,
                                    kGAz, axis_A);
  } else if (test_case == 1) {
    // Run uniform HART-II case
    precone = 2.5 * M_PI / 180.0;  // 2.5 degrees
    r = 0.12;                      // cutout
    freq_normalization = 1.0 / Omega_ref;

    // Set the inertial properties
    TacsScalar mA = 0.95;        // kg/m
    TacsScalar m22 = 1.7e-5;     // kg m
    TacsScalar m33 = 7.0128e-4;  // kg m
    TacsScalar m11 = m22 + m33;  // kg m

    // Axial stiffness
    TacsScalar EA = 1.17e7;  // N

    // Torsional stiffness
    TacsScalar GJ = 1.6e2;  // N m^2

    // Bending stiffness
    TacsScalar EI22 = 2.5e2;
    TacsScalar EI33 = 5.2e3;
    TacsScalar EI23 = 0.0;

    // Shear stiffness
    TacsScalar kG22 = 5.85e5;  // N
    TacsScalar kG33 = 5.85e5;  // N

    // Set the reference axes
    TacsScalar axis_A[] = {0.0, 0.0, 1.0};

    stiff = new TimoshenkoStiffness(axis_A, EA, EI22, EI33, 0.0, GJ, kG22, kG33,
                                    0.0, mA, m11, m22, m33, 0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0);
  } else {
    // Run uniform rectangular geometry
    freq_normalization = 1.0 / Omega_ref;

    // Set the inertial properties
    TacsScalar mA = 2.7;         // kg/m
    TacsScalar m22 = 2.25e-5;    // kg m
    TacsScalar m33 = 2.25e-3;    // kg m
    TacsScalar m11 = m22 + m33;  // kg m
    TacsScalar m23 = 0.0;

    // Axial stiffness
    TacsScalar EA = 70.0e6;  // N

    // Torsional stiffness
    TacsScalar GJ = 811.20;  // N m^2

    // Bending stiffness
    TacsScalar EI22 = 583.3333333333333;
    TacsScalar EI33 = 58333.33333333333;
    TacsScalar EI23 = 0.0;

    // Shear stiffness
    TacsScalar kG22 = 21666666.666666668;  // N
    TacsScalar kG33 = 21666666.666666668;  // N

    // Set the reference axes
    TacsScalar axis_A[] = {0.0, 1.0, 0.0};

    stiff = new TimoshenkoStiffness(mA, m22, m33, m23, EA, GJ, EI22, EI33, kG22,
                                    kG33, axis_A);
  }

  TACSGibbsVector *direction = new TACSGibbsVector(0.0, 0.0, 1.0);

  TACSRevoluteDriver *rd = new TACSRevoluteDriver(direction, angular_rate);

  // Create a rigid link
  TACSGibbsVector *zero = new TACSGibbsVector(0.0, 0.0, 0.0);
  TACSGibbsVector *e1 = new TACSGibbsVector(1.0, 0.0, 0.0);
  TACSGibbsVector *e2 = new TACSGibbsVector(0.0, 1.0, 0.0);
  TACSRefFrame *ref = new TACSRefFrame(zero, e1, e2);

  // Create the rigid body
  TacsScalar mass = 1.0;
  TacsScalar c[] = {0.0, 0.0, 0.0};
  TacsScalar J[] = {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
  TACSGibbsVector *rinit = new TACSGibbsVector(0.0, 0.0, 0.0);
  TACSGibbsVector *vinit = new TACSGibbsVector(0.0, 0.0, 0.0);
  TACSGibbsVector *omegainit = new TACSGibbsVector(0.0, 0.0, angular_rate);
  TACSGibbsVector *g = new TACSGibbsVector(0.0, 0.0, 0.0);
  TACSRigidBody *rb =
      new TACSRigidBody(ref, mass, c, J, rinit, vinit, omegainit, g);

  // Create the rigid link
  TACSRigidLink *rl = new TACSRigidLink(rb);

  // Create the element
  TACSGibbsVector *gravity = new TACSGibbsVector(0.0, 0.0, 0.0);
  MITC3 *beam = new MITC3(stiff, gravity, zero, omegainit);
  beam->incref();

  // Set the number of elements in the beam
  int num_elements = 3 + ne;
  int num_nodes = 3 + 2 * ne + 1;

  // Create the input to TACS
  TACSElement *elements[num_elements];
  int *ptr = new int[num_elements + 1];
  int *conn = new int[3 * num_elements];

  // Set the connectivity
  ptr[0] = 0;
  for (int i = 0; i < ne; i++) {
    ptr[i + 1] = 3 * (i + 1);
    conn[ptr[i]] = 2 * i;
    conn[ptr[i] + 1] = 2 * i + 1;
    conn[ptr[i] + 2] = 2 * i + 2;
    elements[i] = beam;
  }

  int nelems = ne;
  int nnodes = 2 * ne + 1;

  // Set the rigid body element
  int rigid_node = nnodes;
  conn[ptr[nelems]] = nnodes;
  ptr[nelems + 1] = ptr[nelems] + 1;
  elements[nelems] = rb;
  nelems++;
  nnodes++;

  // Set the driver element
  conn[ptr[nelems]] = rigid_node;
  conn[ptr[nelems] + 1] = nnodes;
  ptr[nelems + 1] = ptr[nelems] + 2;
  elements[nelems] = rd;
  nelems++;
  nnodes++;

  // Set the rigid link
  conn[ptr[nelems]] = rigid_node;
  conn[ptr[nelems] + 1] = 0;
  conn[ptr[nelems] + 2] = nnodes;
  ptr[nelems + 1] = ptr[nelems] + 3;
  elements[nelems] = rl;
  nelems++;
  nnodes++;

  // Create TACS itself and set the elements
  const int vars_per_node = 8;
  TACSAssembler *tacs =
      new TACSAssembler(comm, vars_per_node, num_nodes, num_elements);

  tacs->setElementConnectivity(conn, ptr);
  tacs->setElements(elements);
  tacs->initialize();

  // Set the node locations
  TacsScalar *X;
  TACSBVec *Xvec = tacs->createNodeVec();
  Xvec->incref();
  Xvec->getArray(&X);
  for (int i = 0; i < 2 * ne + 1; i++) {
    TacsScalar xdist = i * (L - r) / (2 * ne);
    X[3 * i] = r + xdist;
    X[3 * i + 1] = 0.0;
    X[3 * i + 2] = -xdist * sin(precone);
  }

  tacs->setNodes(Xvec);

  // Set the rotational rate
  int steps_per_rotation = 180;
  int num_rotations = 1;
  double angular_freq = angular_rate / (2 * M_PI);
  double tfinal = num_rotations / angular_freq;
  int num_steps = num_rotations * steps_per_rotation;
  int order = 2;

  TACSIntegrator *integrator =
      new TACSBDFIntegrator(tacs, 0.0, tfinal, num_steps, order);
  integrator->incref();

  integrator->setPrintLevel(0);
  integrator->integrate();

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag =
      (TACSElement::OUTPUT_NODES | TACSElement::OUTPUT_DISPLACEMENTS |
       TACSElement::OUTPUT_STRAINS | TACSElement::OUTPUT_STRESSES |
       TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(tacs, TACS_TIMOSHENKO_BEAM, write_flag);
  f5->incref();
  f5->writeToFile("ucrm.f5");

  integrator->setBeamOutput(f5);

  // integrator->writeSolutionToF5();

  TACSBVec *q, *qdot, *qddot;
  integrator->getStates(num_steps - 1, &q, &qdot, &qddot);

  int nfreq = vars_per_node * num_nodes;
  TacsScalar *freq = new TacsScalar[2 * nfreq];
  nfreq =
      integrator->lapackNaturalFrequencies(gyroscopic, q, qdot, qddot, freq);

  printf("number of frequencies = %d\n", nfreq);

  qsort(freq, nfreq, sizeof(TacsScalar), compare);
  printf("Normalized frequencies\n");
  for (int i = 0; i < 12; i++) {
    printf("%25.15e\n", freq[i] * freq_normalization);
  }

  MPI_Finalize();
  return 0;
}
