#include "Compliance.h"
#include "InducedFailure.h"
#include "KSFailure.h"
#include "MITC9.h"
#include "MITCShell.h"
#include "RigidBody.h"
#include "StructuralMass.h"
#include "TACSCreator.h"
#include "TACSIntegrator.h"
#include "TACSShellTraction.h"
#include "isoFSDTStiffness.h"

/*
  Code for testing adjoints with plate example. Use command line
  arguments as necessary.

  BDF1 BDF2 BDF3    : for BDF integrators
  DIRK2 DIRK3 DIRK4 : for DIRK integrators
  ABM1-6            : for ABM integrators
  NBG               : for Newmark integrator

  test_gradient : to do a complex step verification of the adjoint gradient
  test_element  : to test the element implementation
  num_funcs     : 1 to 3 for the adjoint
  num_threads   : number of threads
  write_solution: write solution to f5 frequency
  print_level: 0, 1, 2
*/
int main(int argc, char **argv) {
  // Intialize MPI and declare communicator
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Default values for important parameters
  int scale = 2;
  int num_funcs = 12;
  int num_threads = 1;
  int test_gradient = 1;
  enum IntegratorType type = DIRK2;

  // Parse command line arguments
  for (int i = 0; i < argc; i++) {
    // Determine the problem size
    if (sscanf(argv[i], "scale=%d", &scale) == 1) {
      if (scale < 0) {
        scale = 1;
      }
      if (scale > 100) {
        scale = 100;
      }
      printf("Problem scale = %d\n", scale);
    }
    // Determine the number of functions for adjoint
    if (sscanf(argv[i], "num_funcs=%d", &num_funcs) == 1) {
      if (num_funcs < 0) {
        num_funcs = 1;
      }
      if (rank == 0) {
        printf("Number of functions : %d\n", num_funcs);
      }
    }

    // Determine the number of threads
    if (sscanf(argv[i], "num_threads=%d", &num_threads) == 1) {
      if (num_threads < 0) {
        num_threads = 1;
      }
      if (num_threads > 24) {
        num_threads = 24;
      }
      if (rank == 0) {
        printf("Number of threads : %d\n", num_threads);
      }
    }

    // Determine whether or not to test gradients with complex step
    if (strcmp("test_gradient", argv[i]) == 0) {
      test_gradient = 1;
      if (rank == 0) {
        printf("Enabled complex-step verification of gradients...\n");
      }
    }

    // Determine the type of integration scheme to use
    // Backward Difference Formulae
    if (strcmp("BDF1", argv[i]) == 0) {
      type = BDF1;
    } else if (strcmp("BDF2", argv[i]) == 0) {
      type = BDF2;
    } else if (strcmp("BDF3", argv[i]) == 0) {
      type = BDF3;

      // Adams-Bashforth-Moulton
    } else if (strcmp("ABM1", argv[i]) == 0) {
      type = ABM1;
    } else if (strcmp("ABM2", argv[i]) == 0) {
      type = ABM2;
    } else if (strcmp("ABM3", argv[i]) == 0) {
      type = ABM3;
    } else if (strcmp("ABM4", argv[i]) == 0) {
      type = ABM4;
    } else if (strcmp("ABM5", argv[i]) == 0) {
      type = ABM5;
    } else if (strcmp("ABM6", argv[i]) == 0) {
      type = ABM6;

      // Diagonally Implicit Runge Kutta
    } else if (strcmp("DIRK2", argv[i]) == 0) {
      type = DIRK2;
    } else if (strcmp("DIRK3", argv[i]) == 0) {
      type = DIRK3;
    } else if (strcmp("DIRK4", argv[i]) == 0) {
      type = DIRK4;

      // Newmark-Beta-Gamma method
    } else if (strcmp("NBGE", argv[i]) == 0) {
      type = NBGE;
    } else if (strcmp("NBG2", argv[i]) == 0) {
      type = NBG2;
    } else if (strcmp("NBG3", argv[i]) == 0) {
      type = NBG3;
    }
  }

  int vars_per_node = 8;
  TACSCreator *creator = new TACSCreator(comm, vars_per_node);
  creator->incref();

  if (rank == 0) {
    int nl_elems = scale * 5;  // number of elements along the length
    int nw_elems = scale * 1;  // number of element along a side of the box

    // Set the number of rigid elements/constraints: 1 rigid body
    // class and 2 constraint classes for each end of the flexible
    // bar.
    int nrigid = 0;
    int num_elements = 4 * nl_elems * nw_elems + nrigid;

    // Set the total number of nodes
    int num_nodes = (2 * nl_elems + 1) * (8 * nw_elems) + nrigid;

    // Set the node locations
    TacsScalar *Xpts = new TacsScalar[3 * num_nodes];

    // Allocate the connectivity
    int conn_size = 9 * num_elements;
    int *elem_ptr = new int[num_elements + 1];
    int *elem_conn = new int[conn_size];
    int *elem_id_nums = new int[num_elements];
    memset(elem_id_nums, 0, num_elements * sizeof(int));

    // Set up the connectivity
    int *ptr = elem_ptr;
    int *conn = elem_conn;
    for (int j = 0; j < 4 * nw_elems; j++) {
      for (int i = 0; i < nl_elems; i++) {
        ptr[0] = conn - elem_conn;
        ptr++;

        for (int jj = 0; jj < 3; jj++) {
          for (int ii = 0; ii < 3; ii++) {
            if (2 * j + jj == 8 * nw_elems) {
              conn[ii + 3 * jj] = 2 * i + ii;
            } else {
              conn[ii + 3 * jj] =
                  2 * i + ii + (2 * nl_elems + 1) * (2 * j + jj);
            }
          }
        }
        conn += 9;
      }
    }
    ptr[0] = conn - elem_conn;

    // Set the length of the flexible bar
    double L = 0.25;
    double Lw = 0.05;

    // Set the node locations
    int nw = 2 * nw_elems;
    int nlen = 2 * nl_elems + 1;
    TacsScalar *x = Xpts;
    for (int j = 0; j < 4 * nw; j++) {
      for (int i = 0; i < nlen; i++) {
        x[0] = (L * i) / (nlen - 1);
        if (j < nw) {
          x[1] = -0.5 * Lw;
          x[2] = -0.5 * Lw + (Lw * j) / nw;
        } else if (j < 2 * nw) {
          x[1] = -0.5 * Lw + (Lw * (j - nw)) / nw;
          x[2] = 0.5 * Lw;
        } else if (j < 3 * nw) {
          x[1] = 0.5 * Lw;
          x[2] = 0.5 * Lw - (Lw * (j - 2 * nw)) / nw;
        } else {
          x[1] = 0.5 * Lw - (Lw * (j - 3 * nw)) / nw;
          x[2] = -0.5 * Lw;
        }
        x += 3;
      }
    }

    int bc = 0;
    int bc_ptr[] = {0, 3};
    int bc_vars[] = {0, 1, 2};
    creator->setBoundaryConditions(1, &bc, bc_ptr, bc_vars);

    // Set the connectivity
    creator->setGlobalConnectivity(num_nodes, num_elements, elem_ptr, elem_conn,
                                   elem_id_nums);

    // Set the nodal locations
    creator->setNodes(Xpts);

    // Free the data
    delete[] Xpts;
    delete[] elem_ptr;
    delete[] elem_conn;
    delete[] elem_id_nums;
  }

  // Create the objects associated with the rigid-bodies
  // ---------------------------------------------------
  // The acceleration due to gravity in global frame of reference
  TACSGibbsVector *gravVec = new TACSGibbsVector(0.0, 0.0, -9.81);

  // Create the constitutive objects
  TacsScalar rho = 2500.0;
  TacsScalar E = 70e9;
  TacsScalar nu = 0.3;
  TacsScalar kcorr = 5.0 / 6.0;
  TacsScalar yield_stress = 464.0e6;
  TacsScalar thickness = 0.01;

  // Create the stiffness object
  isoFSDTStiffness *stiff =
      new isoFSDTStiffness(rho, E, nu, kcorr, yield_stress, thickness, 0);

  // Create the element class
  TACSElement *elem = new MITC9(stiff, gravVec);

  // This call must occur on all processor
  creator->setElements(&elem, 1);

  // Set the reordering type
  creator->setReorderingType(TACSAssembler::RCM_ORDER,
                             TACSAssembler::DIRECT_SCHUR);

  // Create the TACS object
  TACSAssembler *tacs = creator->createTACS();
  tacs->incref();

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag =
      (TACSElement::OUTPUT_NODES | TACSElement::OUTPUT_DISPLACEMENTS |
       TACSElement::OUTPUT_STRAINS | TACSElement::OUTPUT_STRESSES |
       TACSElement::OUTPUT_EXTRAS);

  TACSToFH5 *f5 = new TACSToFH5(tacs, TACS_SHELL, write_flag);
  f5->incref();

  /*-----------------------------------------------------------------*/
  /*------------------ Time Integration and Adjoint Solve -----------*/
  /*-----------------------------------------------------------------*/

  int num_dvs = 1;

  // Create functions of interest
  TACSFunction *func[num_funcs];
  if (num_funcs == 1) {
    func[0] = new TACSCompliance(tacs);
  } else if (num_funcs == 2) {
    func[0] = new TACSKSFailure(tacs, 100.0);

    // Set the induced norm failure types
    TACSInducedFailure *ifunc = new TACSInducedFailure(tacs, 20.0);
    // ifunc->setInducedType(TACSInducedFailure::EXPONENTIAL);
    func[1] = ifunc;

    // func[1] = new TACSCompliance(tacs);
  } else if (num_funcs == 3) {
    func[0] = new TACSKSFailure(tacs, 100.0);
    func[1] = new TACSCompliance(tacs);
    func[2] = new TACSStructuralMass(tacs);
  } else if (num_funcs == 12) {
    // Place functions into the func list
    func[0] = new TACSStructuralMass(tacs);
    func[1] = new TACSCompliance(tacs);

    // Set the discrete and continuous KS functions
    TACSKSFailure *ksfunc = new TACSKSFailure(tacs, 20.0);
    ksfunc->setKSFailureType(TACSKSFailure::DISCRETE);
    func[2] = ksfunc;

    ksfunc = new TACSKSFailure(tacs, 20.0);
    ksfunc->setKSFailureType(TACSKSFailure::CONTINUOUS);
    func[3] = ksfunc;

    // Set the induced norm failure types
    TACSInducedFailure *ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::EXPONENTIAL);
    func[4] = ifunc;

    ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::DISCRETE_EXPONENTIAL);
    func[5] = ifunc;

    ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::DISCRETE_EXPONENTIAL_SQUARED);
    func[6] = ifunc;

    ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::EXPONENTIAL_SQUARED);
    func[7] = ifunc;

    ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::POWER);
    func[8] = ifunc;

    ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::DISCRETE_POWER);
    func[9] = ifunc;

    ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::POWER_SQUARED);
    func[10] = ifunc;

    ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::DISCRETE_POWER_SQUARED);
    func[11] = ifunc;
  }

  for (int i = 0; i < num_funcs; i++) {
    func[i]->incref();
  }

  TacsScalar *funcVals = new TacsScalar[num_funcs];     // adjoint
  TacsScalar *funcValsTmp = new TacsScalar[num_funcs];  // CSD
  TacsScalar *funcVals1 = new TacsScalar[num_funcs];    // forward/reverse

  TacsScalar *dfdx = new TacsScalar[num_funcs * num_dvs];     // adjoint
  TacsScalar *dfdx1 = new TacsScalar[num_funcs * num_dvs];    // CSD
  TacsScalar *dfdxTmp = new TacsScalar[num_funcs * num_dvs];  // forward/reverse

  // Create an array of design variables
  TacsScalar *x = new TacsScalar[num_dvs];
  x[0] = 0.01;

  // Set paramters for time marching
  double tinit = 0.0;
  double tfinal = 1.0e-3;
  for (int k = 0; k < argc; k++) {
    if (sscanf(argv[k], "tfinal=%lf", &tfinal) == 1) {
    }
  }
  int num_steps_per_sec = 10000;

  TACSIntegrator *obj =
      TACSIntegrator::getInstance(tacs, tinit, tfinal, num_steps_per_sec, type);
  obj->incref();

  // Set optional parameters
  obj->setOrderingType(TACSAssembler::RCM_ORDER);
  obj->setRelTol(1.0e-8);
  obj->setAbsTol(1.0e-6);
  obj->setMaxNewtonIters(25);
  obj->setPrintLevel(2);
  obj->setJacAssemblyFreq(1);
  obj->setOutputFrequency(0);
  obj->setShellOutput(0);

  // Set functions of interest for adjoint solve
  obj->setFunction(func, num_funcs);

  // Get the adjoint gradient
  double t0 = MPI_Wtime();
  obj->getFuncGrad(num_dvs, x, funcVals, dfdx);
  t0 = MPI_Wtime() - t0;

  // Print a summary of time taken
  if (rank == 0) {
    obj->printWallTime(t0, 2);
  }

  // Test the adjoint derivatives if sought
  if (test_gradient) {
    // The maximum number of gradient components to test
    // using finite-difference/complex-step

    // Scan any remaining arguments that may be required
#ifdef TACS_USE_COMPLEX
    double dh = 1.0e-30;
#else
    double dh = 1.0e-8;
#endif
    for (int k = 0; k < argc; k++) {
      if (sscanf(argv[k], "dh=%lf", &dh) == 1) {
      }
    }

    // Complex step verification
    obj->getFDFuncGrad(num_dvs, x, funcValsTmp, dfdxTmp, dh);

    // Print out the finite-difference interval
    if (rank == 0) {
#ifdef TACS_USE_COMPLEX
      printf("Complex-step interval: %le\n", dh);
#else
      printf("Finite-difference interval: %le\n", dh);
#endif
    }

    if (rank == 0) {
      printf("Structural sensitivities\n");
      for (int j = 0; j < num_funcs; j++) {
        printf("Sensitivities for function %s %25.15e\n",
               func[j]->functionName(), TacsRealPart(funcValsTmp[j]));
        printf("%25s %25s %25s %25s\n", "Adjoint", "FD/CS", "Abs. error",
               "Rel. error");
        for (int k = 0; k < num_dvs; k++) {
          printf("%25.15e %25.15e %25.15e %25.15e\n",
                 TacsRealPart(dfdx[k + j * num_dvs]),
                 TacsRealPart(dfdxTmp[k + j * num_dvs]),
                 TacsRealPart(dfdx[k + j * num_dvs]) -
                     TacsRealPart(dfdxTmp[k + j * num_dvs]),
                 (TacsRealPart(dfdx[k + j * num_dvs]) -
                  TacsRealPart(dfdxTmp[k + j * num_dvs])) /
                     TacsRealPart(dfdxTmp[k + j * num_dvs]));
        }
      }
    }
  }

  /*
  // Delete objects
  bodyA->decref();
  bodyB->decref();
  conA->decref();
  conB->decref();

  tacs->decref();
  f5->decref();
  obj->decref();
  */

  creator->decref();
  tacs->decref();
  obj->decref();

  for (int i = 0; i < num_funcs; i++) {
    func[i]->decref();
  }

  delete[] x;
  delete[] funcVals;
  delete[] funcVals1;
  delete[] funcValsTmp;

  MPI_Finalize();

  return (0);
}

/*
// Define the zero vector
TACSGibbsVector *zero = new TACSGibbsVector(0.0, 0.0, 0.0);

// Construct the frame of reference
TACSGibbsVector *rA0Vec = new TACSGibbsVector(0.0, 0.0, 0.0);
TACSGibbsVector *rA1Vec = new TACSGibbsVector(1.0, 0.0, 0.0);
TACSGibbsVector *rA2Vec = new TACSGibbsVector(0.0, 1.0, 0.0);
TACSRefFrame *refFrameA = new TACSRefFrame(rA0Vec, rA1Vec, rA2Vec);

// Define the inertial properties
const TacsScalar mA    = 1.0;
const TacsScalar cA[3] = {0.0, 0.0, 0.0};
const TacsScalar JA[6] = {1.0/3.0, 0.0, 0.0,
1.0/3.0, 0.0,
1.0/3.0};

// Define dynamics properties
TACSGibbsVector *rAInitVec = new TACSGibbsVector(0.0, 2.5, 0.0);

// Create visualization
TACSRigidBodyViz *vizA = new TACSRigidBodyViz(0.5, 5.0, 0.5);

// Construct a rigid body
TACSRigidBody *bodyA = new  TACSRigidBody(refFrameA,
mA, cA, JA,
rAInitVec, zero, zero, gravVec);
bodyA->setVisualization(vizA);
bodyA->incref();
*/
