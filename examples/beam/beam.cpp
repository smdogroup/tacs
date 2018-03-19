#include "TimoshenkoStiffness.h"
#include "MITC3.h"
#include "TACSAssembler.h"
#include "RigidBody.h"
#include "TACSIntegrator.h"
#include "KinematicConstraints.h"

/*
  Create the one bar
*/
TACSAssembler *one_bar_beam( int nA ){
    // Set the gravity vector
  TACSGibbsVector *gravity = new TACSGibbsVector(0.0, 0.0, -9.81);

  // Set the points b, c and d
  TACSGibbsVector *ptA = new TACSGibbsVector(0.0, 0.0, 0.0);
  TACSGibbsVector *ptB = new TACSGibbsVector(1.0, 0.0, 0.0);

  // Create the revolute direction for B and D
  TACSGibbsVector *revDirA = new TACSGibbsVector(0.0, 0.0, 1.0);
  TACSGibbsVector *revDirB = new TACSGibbsVector(0.0, 0.0, 1.0);

  // Create the revolute constraints
  int fixed_point = 1;
  TACSRevoluteConstraint *revA = 
    new TACSRevoluteConstraint(fixed_point, ptA, revDirA);
  TACSRevoluteConstraint *revB = 
    new TACSRevoluteConstraint(fixed_point, ptB, revDirB);

  // Create the stiffness objects for each element
  TacsScalar mA = 1.997; // kg/m
  TacsScalar IA = 42.60e-6; // kg*m

  TacsScalar EA_A = 52.99e6;
  TacsScalar GJ_A = 733.5;
  TacsScalar kGAz_A = 16.88e6; 
  TacsScalar EIz_A = 1131.0;

  // Set the reference axes for each beam
  TacsScalar axis_A[] = {0.0, 1.0, 0.0};
  
  // Create the Timoshenko stiffness object
  TimoshenkoStiffness *stiffA =
    new TimoshenkoStiffness(mA, IA, IA, 0.0,
                            EA_A, GJ_A, EIz_A, EIz_A, kGAz_A, kGAz_A,
                            axis_A);

  // Set up the connectivity
  MITC3 *beamA = new MITC3(stiffA, gravity);

  // Set the number of nodes in the mesh
  int nnodes = 2*nA + 3;

  // Set the number of elements
  int nelems = nA + 2;

  // Create the connectivities
  TacsScalar *X = new TacsScalar[ 3*nnodes ];
  memset(X, 0, 3*nnodes*sizeof(TacsScalar));

  int *ptr = new int[ nelems+1 ];
  int *conn = new int[ 3*nelems ];
  TACSElement **elems = new TACSElement*[ nelems ];

  // Set the nodes numbers and locations
  int *nodesA = new int[ 2*nA+1 ];
  int n = 0;
  for ( int i = 0; i < 2*nA+1; i++, n++ ){
    nodesA[i] = n;
    X[3*n+1] = 0.12*i/(2*nA); 
  }

  // Set the connectivity for the beams
  int elem = 0;
  ptr[0] = 0;
  for ( int i = 0; i < nA; i++ ){
    conn[ptr[elem]] = nodesA[2*i];
    conn[ptr[elem]+1] = nodesA[2*i+1];
    conn[ptr[elem]+2] = nodesA[2*i+2];
    elems[elem] = beamA;
    ptr[elem+1] = ptr[elem] + 3;
    elem++;
  }

  // Add the connectivities for the constraints
  conn[ptr[elem]] = nodesA[0];
  conn[ptr[elem]+1] = nnodes-2;
  elems[elem] = revA;
  ptr[elem+1] = ptr[elem] + 2;
  elem++;

  conn[ptr[elem]] = nodesA[2*nA];
  conn[ptr[elem]+1] = nnodes-1;
  elems[elem] = revB;
  ptr[elem+1] = ptr[elem] + 2;
  elem++;

  delete [] nodesA;

  // Create the TACSAssembler object
  TACSAssembler *tacs = new TACSAssembler(MPI_COMM_WORLD, 8, nnodes, nelems);

  tacs->setElementConnectivity(conn, ptr);
  delete [] conn;
  delete [] ptr;

  tacs->setElements(elems);
  delete [] elems;

  tacs->initialize();

  // Set the node locations
  TACSBVec *Xvec = tacs->createNodeVec();
  Xvec->incref();
  TacsScalar *Xarray;
  Xvec->getArray(&Xarray);
  memcpy(Xarray, X, 3*nnodes*sizeof(TacsScalar));
  tacs->setNodes(Xvec);
  Xvec->decref();
  delete [] X;

  return tacs;
}

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
TACSAssembler *four_bar_mechanism( int nA, int nB, int nC ){
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
  TacsScalar theta = (5.0*M_PI/180.0);
  TACSGibbsVector *revDirC = new TACSGibbsVector(sin(theta), 0.0, cos(theta));

  // Create the revolute constraints
  TacsScalar omega = -0.6; // rad/seconds
  int fixed_point = 1;
  int not_fixed = 0;
  
  TACSRevoluteDriver *revDriverA =
    new TACSRevoluteDriver(revDirA, omega);
  TACSRevoluteConstraint *revB = 
    new TACSRevoluteConstraint(not_fixed, ptB, revDirB);
  TACSRevoluteConstraint *revC = 
    new TACSRevoluteConstraint(not_fixed, ptC, revDirC);
  TACSRevoluteConstraint *revD = 
    new TACSRevoluteConstraint(fixed_point, ptD, revDirD);

  // Create the stiffness objects for each element
  TacsScalar mA = 1.997; // kg/m
  TacsScalar IA = 42.60e-6; // kg*m

  TacsScalar EA_A = 52.99e6;
  TacsScalar GJ_A = 733.5;
  TacsScalar kGAz_A = 16.88e6; 
  TacsScalar EIz_A = 1131.0;

  // The properties of the second beam
  TacsScalar mB = 0.4992; // kg*m^2/m
  TacsScalar IB = 2.662e-6; // kg*m^2/m

  TacsScalar EA_B = 13.25e6;
  TacsScalar GJ_B = 45.84;
  TacsScalar kGAz_B = 4.220e6;
  TacsScalar EIz_B = 70.66;

  // Set the reference axes for each beam
  TacsScalar axis_A[] = {-1.0, 0.0, 0.0};
  TacsScalar axis_B[] = {0.0, 1.0, 0.0};
  TacsScalar axis_C[] = {1.0, 0.0, 0.0};
  
  // Create the Timoshenko stiffness object
  TimoshenkoStiffness *stiffA =
    new TimoshenkoStiffness(mA, IA, IA, 0.0,
                            EA_A, GJ_A, EIz_A, EIz_A, kGAz_A, kGAz_A,
                            axis_A);

  TimoshenkoStiffness *stiffB =
    new TimoshenkoStiffness(mA, IA, IA, 0.0,
                            EA_A, GJ_A, EIz_A, EIz_A, kGAz_A, kGAz_A,
                            axis_B);

  TimoshenkoStiffness *stiffC =
    new TimoshenkoStiffness(mB, IB, IB, 0.0,
                            EA_B, GJ_B, EIz_B, EIz_B, kGAz_B, kGAz_B,
                            axis_C);

  // Set up the connectivity
  MITC3 *beamA = new MITC3(stiffA, gravity);
  MITC3 *beamB = new MITC3(stiffB, gravity);
  MITC3 *beamC = new MITC3(stiffC, gravity);

  // Set the number of nodes in the mesh
  int nnodes = (2*nA+1) + (2*nB+1) + (2*nC+1) + 4;

  // Set the number of elements
  int nelems = nA + nB + nC + 4;

  // Create the connectivities
  TacsScalar *X = new TacsScalar[ 3*nnodes ];
  memset(X, 0, 3*nnodes*sizeof(TacsScalar));

  int *ptr = new int[ nelems+1 ];
  int *conn = new int[ 3*nelems ];
  TACSElement **elems = new TACSElement*[ nelems ];

  // Set the nodes numbers and locations
  int *nodesA = new int[ 2*nA+1 ];
  int *nodesB = new int[ 2*nB+1 ];
  int *nodesC = new int[ 2*nC+1 ];
  int n = 0;
  for ( int i = 0; i < 2*nA+1; i++, n++ ){
    nodesA[i] = n;
    X[3*n+1] = 0.12*i/(2*nA); 
  }
  for ( int i = 0; i < 2*nB+1; i++, n++ ){
    nodesB[i] = n;
    X[3*n] = 0.24*i/(2*nB);
    X[3*n+1] = 0.12; 
  }
  for ( int i = 0; i < 2*nC+1; i++, n++ ){
    nodesC[i] = n;
    X[3*n] = 0.24;
    X[3*n+1] = 0.12*(1.0 - 1.0*i/(2*nC));
  }

  // Set the connectivity for the beams
  int elem = 0;
  ptr[0] = 0;
  for ( int i = 0; i < nA; i++ ){
    conn[ptr[elem]] = nodesA[2*i];
    conn[ptr[elem]+1] = nodesA[2*i+1];
    conn[ptr[elem]+2] = nodesA[2*i+2];
    elems[elem] = beamA;
    ptr[elem+1] = ptr[elem] + 3;
    elem++;
  }

  for ( int i = 0; i < nB; i++ ){
    conn[ptr[elem]] = nodesB[2*i];
    conn[ptr[elem]+1] = nodesB[2*i+1];
    conn[ptr[elem]+2] = nodesB[2*i+2];
    elems[elem] = beamB;
    ptr[elem+1] = ptr[elem] + 3;
    elem++;
  }

  for ( int i = 0; i < nC; i++ ){
    conn[ptr[elem]] = nodesC[2*i];
    conn[ptr[elem]+1] = nodesC[2*i+1];
    conn[ptr[elem]+2] = nodesC[2*i+2];
    elems[elem] = beamC;
    ptr[elem+1] = ptr[elem] + 3;
    elem++;
  }

  // Add the connectivities for the constraints
  conn[ptr[elem]] = nodesA[0];
  conn[ptr[elem]+1] = nnodes-4;
  elems[elem] = revDriverA;
  ptr[elem+1] = ptr[elem] + 2;
  elem++;

  conn[ptr[elem]] = nodesA[2*nA];
  conn[ptr[elem]+1] = nodesB[0];
  conn[ptr[elem]+2] = nnodes-3;
  elems[elem] = revB;
  ptr[elem+1] = ptr[elem] + 3;
  elem++;

  conn[ptr[elem]] = nodesC[0];
  conn[ptr[elem]+1] = nodesB[2*nB];
  conn[ptr[elem]+2] = nnodes-2;
  elems[elem] = revC;
  ptr[elem+1] = ptr[elem] + 3;
  elem++;

  conn[ptr[elem]] = nodesC[2*nC];
  conn[ptr[elem]+1] = nnodes-1;
  elems[elem] = revD;
  ptr[elem+1] = ptr[elem] + 2;
  elem++;

  delete [] nodesA;
  delete [] nodesB;
  delete [] nodesC;

  // Create the TACSAssembler object
  TACSAssembler *tacs = new TACSAssembler(MPI_COMM_WORLD, 8, nnodes, nelems);

  tacs->setElementConnectivity(conn, ptr);
  delete [] conn;
  delete [] ptr;

  tacs->setElements(elems);
  delete [] elems;

  tacs->initialize();

  // Set the node locations
  TACSBVec *Xvec = tacs->createNodeVec();
  Xvec->incref();
  TacsScalar *Xarray;
  Xvec->getArray(&Xarray);
  memcpy(Xarray, X, 3*nnodes*sizeof(TacsScalar));
  tacs->setNodes(Xvec);
  Xvec->decref();
  delete [] X;

  return tacs;
}

/*
  Test the element implementation
*/
void test_beam_element(){
    // Set the reference axis
  TacsScalar axis[] = {0.0, 1.0, 0.0};

  // Set the gravity vector
  TACSGibbsVector *gravity = new TACSGibbsVector(0.0, 0.0, -9.81);
  
  // Set the element properties  
  TacsScalar rhoA = 1.5;
  TacsScalar rhoIy = 0.15;
  TacsScalar rhoIz = 0.15;
  TacsScalar rhoIyz = 0.0;

  TacsScalar EA = 1e4;
  TacsScalar GJ = 1.50e4;
  TacsScalar EIy = 2.4e4;
  TacsScalar EIz = 3.24e4;
  TacsScalar kGAy = 2.5e3;
  TacsScalar kGAz = 5.2e3;
  
  // Create the Timoshenko stiffness object
  TimoshenkoStiffness *stiff =
    new TimoshenkoStiffness(rhoA, rhoIy, rhoIz, rhoIyz,
                            EA, GJ, EIy, EIz, kGAy, kGAz,
                            axis);
  stiff->incref();

  // Create the MITC3 element
  MITC3 *beam = new MITC3(stiff, gravity);
  beam->incref();

  int test_element = 1;
  if (test_element){
    TacsScalar X[] = {0.0, 0.0, 0.0,
                      0.375, 0.5, 0.1,
                      1.0, 1.0, 0.2};
    beam->testStrain(X);
   
    int multipliers[3] = {7, 15, 23};
    TacsScalar vars[24], dvars[24], ddvars[24];
    for ( int i = 0; i < 24; i++ ){
      vars[i] = -1.0 + 2.0*rand()/RAND_MAX;
      dvars[i] = -1.0 + 2.0*rand()/RAND_MAX;
      ddvars[i] = -1.0 + 2.0*rand()/RAND_MAX;
    }
    
    beam->setStepSize(5e-6);
    beam->setPrintLevel(2);
    beam->testResidual(0.0, X, vars, dvars, ddvars, multipliers, 3);
    beam->testJacobian(0.0, X, vars, dvars, ddvars);
  }

  int test_average = 1;
  if (test_average){
    TacsScalar X[] = {0.0, 0.0, 0.0,
                      0.0, 3.0, 0.1,
                      1.0, 3.0, 0.2,
                      2.0, 4.0, 0.3,
                      3.0, 3.0, 1.0};
    
    int nmultipliers = 6;
    int multipliers[] = {32, 33, 34, 35, 36, 37};
    TacsScalar vars[40], dvars[40], ddvars[40];
    for ( int i = 0; i < 40; i++ ){
      vars[i] = -1.0 + 2.0*rand()/RAND_MAX;
      dvars[i] = -1.0 + 2.0*rand()/RAND_MAX;
      ddvars[i] = -1.0 + 2.0*rand()/RAND_MAX;
    }
        
    // Construct the frame of reference
    TACSGibbsVector *rAInitVec = new TACSGibbsVector(5.2, 5.3, 5.4); 
    TACSGibbsVector *rA1Vec = new TACSGibbsVector(5.2+1.0, 5.3, 5.4);
    TACSGibbsVector *rA2Vec = new TACSGibbsVector(5.2, 5.3+1.0, 5.4);
    TACSRefFrame *refFrame = new TACSRefFrame(rAInitVec, rA1Vec, rA2Vec);

    // Define the inertial properties
    const TacsScalar mA    = 6.0;
    const TacsScalar cA[3] = {20.0, 14.0, 42.0};
    const TacsScalar JA[6] = {1.0, 0.8, -0.7,
                              2.0, 1.4,
                              3.0};
    // Construct a rigid body
    TACSRigidBody *bodyA = new TACSRigidBody(refFrame,
                                             mA, cA, JA,
                                             rAInitVec, rAInitVec, rAInitVec,
                                             gravity);

    // Test the revolute constraint
    TACSGibbsVector *point = new TACSGibbsVector(0.5, 1.0, -2.5);

    int moment_flag = 7;
    TACSAverageConstraint *avg = 
      new TACSAverageConstraint(bodyA, point, refFrame, moment_flag);

    avg->setStepSize(5e-6);
    avg->setPrintLevel(2);
    avg->testResidual(0.0, X, vars, dvars, ddvars, 
                      multipliers, nmultipliers);
    avg->testJacobian(0.0, X, vars, dvars, ddvars);
  }
}

int main( int argc, char *argv[] ){
  // Initialize MPI
  MPI_Init(&argc, &argv);

  int test_beam = 0;
  for ( int k = 0; k < argc; k++ ){
    if (strcmp(argv[k], "test_beam") == 0){
      test_beam = 1;
    }
  }

  if (test_beam){
    test_beam_element();

    // Create the stiffness objects for each element
    TacsScalar mA = 1.997; // kg/m
    TacsScalar IA = 42.60e-6; // kg*m

    TacsScalar EA_A = 52.99e6;
    TacsScalar GJ_A = 733.5;
    TacsScalar kGAz_A = 16.88e6; 
    TacsScalar EIz_A = 1131.0;

    // Set the reference axes
    TacsScalar axis_A[] = {0.0, 1.0, 0.0};
  
    // Create the Timoshenko stiffness object
    TimoshenkoStiffness *stiffA =
      new TimoshenkoStiffness(mA, IA, IA, 0.0,
                              EA_A, GJ_A, EIz_A, EIz_A, kGAz_A, kGAz_A,
                              axis_A);

    // Create the element
    TACSGibbsVector *gravity = new TACSGibbsVector(0.0, 0.0, -9.81);
    MITC3 *beam = new MITC3(stiffA, gravity);
    beam->incref();

    // Test the straight beam for bending/torsion/extension 
    // relationships
    TacsScalar X[] = {-1.5, 0., 0., 0.0, 0., 0., 2.5, 0., 0.};
    TacsScalar vars[24];

    // Set the distribution of displacements/quaternions
    TacsScalar eps = 1e-3;
    TacsScalar u[] = {0., 0., 0.};
    TacsScalar ux[] = {0., 0., 0.};
    TacsScalar theta0[] = {0., 0., 0.};
    TacsScalar thetax[] = {1., 0., 0.};

    memset(vars, 0, 24*sizeof(TacsScalar));
    for ( int k = 0; k < 3; k++ ){
      vars[8*k] = u[0] + eps*X[3*k]*ux[0];
      vars[8*k+1] = u[1] + eps*X[3*k]*ux[1];
      vars[8*k+2] = u[2] + eps*X[3*k]*ux[2];

      vars[8*k+3] = 1.0;
      vars[8*k+4] = 0.5*eps*(theta0[0] + X[3*k]*thetax[0]);
      vars[8*k+5] = 0.5*eps*(theta0[1] + X[3*k]*thetax[1]);
      vars[8*k+6] = 0.5*eps*(theta0[2] + X[3*k]*thetax[2]);
    }

    printf("%5s %15s %15s %15s %15s %15s %15s\n", 
           "pt", "ex0", "et0", "ey1", "ez2", "exy", "exz");

    for ( int i = 1; i < 11; i += 2 ){
      double pt = -1.0 + 0.2*i;

      TacsScalar e[6];
      beam->getStrain(e, &pt, X, vars);
      printf("%5.2f %15.4e %15.4e %15.4e %15.4e %15.4e %15.4e\n",
             pt, e[0], e[1], e[2], e[3], e[4], e[5]);
    }

    beam->decref();
  }
  else {
    // Create the finite-element model
    int nA = 4, nB = 8, nC = 4;
    TACSAssembler *tacs = four_bar_mechanism(nA, nB, nC);
    tacs->incref();

    // Set the final time
    double tf = 12.0;

    // The number of total steps (100 per second)
    int num_steps = 1200;

    // Create the integrator class
    TACSIntegrator *integrator = 
      new TACSBDFIntegrator(tacs, 0.0, tf, num_steps, 2);
    integrator->incref();

    // Set the integrator options
    integrator->setUseFEMat(1, TACSAssembler::TACS_AMD_ORDER);
    integrator->setAbsTol(1e-7);
    integrator->setOutputFrequency(10);
    integrator->setPrintLevel(-1);
    integrator->integrate();

    // Set the output options/locations
    int elem[3];
    // elem[0] = nA/2;
    // elem[1] = nA + nB/2;
    // elem[2] = nA + nB + nC/2;
    // double param[][1] = {{-1.0}, {-1.0}, {-1.0}}; 
    elem[0] = nA/2;
    elem[1] = nA + nB/2;
    elem[2] = nA + nB + nC/2;
    double param[][1] = {{-1.0}, {-1.0}, {0.0}}; 

    // Extra the data to a file
    for ( int pt = 0; pt < 3; pt++ ){
      char filename[128];
      sprintf(filename, "mid_beam_%d.dat", pt+1);
      FILE *fp = fopen(filename, "w");

      fprintf(fp, "Variables = t, u0, v0, w0, sx0, st0, sy1, sz1, sxy0, sxz0\n");

      // Write out data from the beams
      TACSBVec *q = NULL;
      for ( int k = 0; k < num_steps+1; k++ ){
        TacsScalar X[3*3], vars[8*3];
        double time = integrator->getStates(k, &q, NULL, NULL);
        tacs->setVariables(q);
        TACSElement *element = tacs->getElement(elem[pt], X, vars);

        TacsScalar e[6], s[6];
        element->getStrain(e, param[pt], X, vars);
        TACSConstitutive *con = element->getConstitutive();
        con->calculateStress(param[pt], e, s);

        fprintf(fp, "%e  %e %e %e  %e %e %e  %e %e %e\n",
                time, vars[0], vars[1], vars[2], 
                s[0], s[1], s[2], s[3], s[4], s[5]);
      }
      fclose(fp);
    }

    integrator->decref();
    tacs->decref();
  }

  MPI_Finalize();
  return 0;
}
