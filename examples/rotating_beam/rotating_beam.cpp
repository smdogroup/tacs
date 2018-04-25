#include "TACSAssembler.h"
#include "TACSIntegrator.h"
#include "MITC3.h"
#include "KinematicConstraints.h"
#include "TimoshenkoStiffness.h"
#include "TACSElementAlgebra.h"
#include "RigidBody.h"

int compare( const void *a, const void *b ){
  if (*(TacsScalar*)a < *(TacsScalar*)b){
    return -1;
  }
  return 1;
}

int main( int argc, char *argv[] ){
  // Initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  int ne = 10;
  double Omega = 10.0;
  int gyroscopic = 1;

  for ( int i = 0; i < argc; i++ ){
    if (sscanf(argv[i], "Omega=%lf", &Omega) == 1){
      printf("Omega = %e\n", Omega);
    }
    if (sscanf(argv[i], "gyroscopic=%d", &gyroscopic) == 1){
      printf("gyroscopic = %d\n", gyroscopic);
    }
    if (sscanf(argv[i], "ne=%d", &ne) == 1){
      printf("ne = %d\n", ne);
    }
  }

  // Create the motion driver
  const TacsScalar angular_rate = Omega;

  // Set the length of the rotating beam
  const TacsScalar r = 0.4;
  const TacsScalar L = 2.0;

  // h = sqrt((12*L^2)/alpha^2)
  const TacsScalar h = 0.0769800358919501;

  // Set the material properties
  const TacsScalar E = 70e9; // Pa
  const TacsScalar rho = 2600.0; // kg/m^3
  const TacsScalar nu = 0.33; // Poisson's ratio
  const TacsScalar Gmod = 0.5*E/(1.0 + nu);
  const TacsScalar kappa = 5.0/6.0;

  // Set the stiffness properties
  TacsScalar mA = rho*h*h;
  TacsScalar IA = rho*h*h*h*h/12.0;

  TacsScalar EA = E*h*h;
  TacsScalar GJ = 2.25*Gmod*h*h*h*h;
  TacsScalar kGAz = kappa*Gmod*h*h;
  TacsScalar EIz = E*h*h*h*h/12.0;

  // Print out the properties to the screen
  TacsScalar delta = r/L;
  TacsScalar alpha = sqrt(mA*L*L/IA);
  TacsScalar eta = EIz/(mA*L*L*L*L);
  TacsScalar T = 1.0/sqrt(eta);
  TacsScalar gamma = T*Omega;
  printf("delta: %25.15e\n", delta);
  printf("alpha: %25.15e\n", alpha);
  printf("eta:   %25.15e\n", eta);
  printf("T:     %25.15e\n", T);
  printf("gamma: %25.15e\n", gamma);

  // Set the reference axes
  TacsScalar axis_A[] = {0.0, 1.0, 0.0};

  TACSGibbsVector *direction = new TACSGibbsVector(0.0, 0.0, 1.0);

  const int fix_rotations = 1;
  TACSRevoluteDriver *rd = 
    new TACSRevoluteDriver(direction, angular_rate);

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

  // Create the Timoshenko stiffness object
  TimoshenkoStiffness *stiffA =
    new TimoshenkoStiffness(mA, IA, IA, 0.0,
                            EA, GJ, EIz, EIz, kGAz, kGAz,
                            axis_A);

  // Create the element
  TACSGibbsVector *gravity = new TACSGibbsVector(0.0, 0.0, 0.0);
  MITC3 *beam = new MITC3(stiffA, gravity, zero, omegainit);
  beam->incref();

  // Set the number of elements in the beam
  int num_elements = 3 + ne;
  int num_nodes = 3 + 2*ne+1;

  // Create the input to TACS
  TACSElement *elements[num_elements];
  int *ptr = new int[num_elements+1];
  int *conn = new int[3*num_elements];

  // Set the connectivity
  ptr[0] = 0;
  for ( int i = 0; i < ne; i++ ){
    ptr[i+1] = 3*(i+1);
    conn[ptr[i]] = 2*i;
    conn[ptr[i]+1] = 2*i+1;
    conn[ptr[i]+2] = 2*i+2;
    elements[i] = beam;
  }

  int nelems = ne;
  int nnodes = 2*ne+1;

  // Set the rigid body element
  int rigid_node = nnodes;
  conn[ptr[nelems]] = nnodes; 
  ptr[nelems+1] = ptr[nelems]+1;
  elements[nelems] = rb;
  nelems++;
  nnodes++;

  // Set the driver element
  conn[ptr[nelems]] = rigid_node;
  conn[ptr[nelems]+1] = nnodes;
  ptr[nelems+1] = ptr[nelems]+2;
  elements[nelems] = rd;
  nelems++;
  nnodes++;

  // Set the rigid link
  conn[ptr[nelems]] = rigid_node;
  conn[ptr[nelems]+1] = 0;
  conn[ptr[nelems]+2] = nnodes;
  ptr[nelems+1] = ptr[nelems]+3;
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
  for ( int i = 0; i < 2*ne+1; i++ ){
    X[3*i] = r + i*(L - r)/(2*ne);
  }

  tacs->setNodes(Xvec);

  // Set the rotational rate
  int steps_per_rotation = 180;
  int num_rotations = 5;
  double angular_freq = angular_rate/(2*M_PI);
  double tfinal = num_rotations/angular_freq;
  int num_steps = num_rotations*steps_per_rotation + 1;
  int order = 2;

  TACSIntegrator *integrator =
    new TACSBDFIntegrator(tacs, 0.0, tfinal, num_steps, order);
  integrator->incref();

  integrator->setPrintLevel(1);
  integrator->integrate();

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 * f5 = new TACSToFH5(tacs, TACS_TIMOSHENKO_BEAM, write_flag);
  f5->incref();
  f5->writeToFile("ucrm.f5");

  integrator->setBeamOutput(f5);

  integrator->writeSolutionToF5();

  TACSBVec *q, *qdot, *qddot;
  integrator->getStates(num_steps-1, &q, &qdot, &qddot);

  int nfreq = vars_per_node*num_nodes;
  TacsScalar *freq = new TacsScalar[ 2*nfreq ];
  nfreq = integrator->lapackNaturalFrequencies(gyroscopic, 
                                               q, qdot, qddot, freq);

  printf("number of frequencies = %d\n", nfreq);

  qsort(freq, nfreq, sizeof(TacsScalar), compare);
  for ( int i = 0; i < 10; i++ ){
    printf("Normalized frequency[%3d]: %25.15e\n", i, freq[i]*T);
  }

  // // Compute the kinetic energy
  // TacsScalar angular_rate = 5.0;
  // TacsScalar angular_dir[3] = {0.0, 1.0, 0.0};
  // TacsScalar omega[3];
  // omega[0] = angular_rate*angular_dir[0];
  // omega[1] = angular_rate*angular_dir[1];
  // omega[2] = angular_rate*angular_dir[2];

  // // Test the straight beam for bending/torsion/extension 
  // // relationships
  // TacsScalar X[] = {0.0, 0.0, 0.0, 
  //                   0.5*L, 0.0, 0.0,
  //                   1.0*L, 0.0, 0.0};
  // TacsScalar vars[24], dvars[24], ddvars[24];

  // for ( int i = 0; i < 3; i++ ){
  //   // Set the displacement
  //   TacsScalar u[3] = {0.0, 0.0, 0.0};
  //   TacsScalar q[4] = {1.0, 0.0, 0.0, 0.0};

  //   // Set the velocity/time derivative components
  //   TacsScalar v[3], qdot[4];
  //   crossProduct(1.0, omega, &X[3*i], v);
  //   qdot[0] = 0.0;
  //   qdot[1] = 0.5*angular_rate*angular_dir[0];
  //   qdot[2] = 0.5*angular_rate*angular_dir[1];
  //   qdot[3] = 0.5*angular_rate*angular_dir[2];

  //   // Set the acceleration
  //   TacsScalar a[3], qddot[4];
  //   crossProduct(1.0, omega, v, a);
  //   qddot[0] = -0.25*angular_rate*angular_rate;
  //   qddot[1] = 0.0;
  //   qddot[2] = 0.0;
  //   qddot[3] = 0.0;
    
  //   // Set the displacement variables
  //   vars[8*i] = u[0];
  //   vars[8*i+1] = u[1];
  //   vars[8*i+2] = u[2];
  //   vars[8*i+3] = q[0];
  //   vars[8*i+4] = q[1];
  //   vars[8*i+5] = q[2];
  //   vars[8*i+6] = q[3];
  //   vars[8*i+7] = 0.0;

  //   // Set the velocity
  //   dvars[8*i] = v[0];
  //   dvars[8*i+1] = v[1];
  //   dvars[8*i+2] = v[2];
  //   dvars[8*i+3] = qdot[0];
  //   dvars[8*i+4] = qdot[1];
  //   dvars[8*i+5] = qdot[2];
  //   dvars[8*i+6] = qdot[3];
  //   dvars[8*i+7] = 0.0;

  //   // Set the acceleration
  //   ddvars[8*i] = a[0];
  //   ddvars[8*i+1] = a[1];
  //   ddvars[8*i+2] = a[2];
  //   ddvars[8*i+3] = qddot[0];
  //   ddvars[8*i+4] = qddot[1];
  //   ddvars[8*i+5] = qddot[2];
  //   ddvars[8*i+6] = qddot[3];
  //   ddvars[8*i+7] = 0.0;
  // }

  // double time = 0.0;

  // // Compute the kinetic energy
  // TacsScalar Te, Pe;
  // beam->computeEnergies(time, &Te, &Pe, X, vars, dvars);

  // printf("Kinetic energy: %25.15e\n", Te);
  // printf("Kinetic energy: %25.15e\n",
  //        0.5*angular_rate*angular_rate*mA*L*L*L/3.0);

  // int multipliers[3] = {7, 15, 23};
  // beam->setStepSize(1e-5);
  // beam->setPrintLevel(2);
  // beam->testResidual(time, X, vars, dvars, ddvars, multipliers, 3);

  MPI_Finalize();
  return 0;
}