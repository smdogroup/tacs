#include "TACSCreator.h"
#include "TACSIntegrator.h"
#include "RigidBody.h"
#include "MITC9.h"
#include "MITCShell.h"
#include "TACSShellTraction.h"
#include "isoFSDTStiffness.h"
#include "KSFailure.h"
#include "StructuralMass.h"
#include "Compliance.h"

/*
  Set up a double pendulum with a flexible attachment
*/
int main( int argc, char *argv[] ){
  // Initialize MPI
  MPI_Init(&argc, &argv);  
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank; 
  MPI_Comm_rank(comm, &rank); 

  // Default values for important parameters
  int scale = 2;
  int num_funcs = 1;
  int num_threads = 1;
  int test_gradient = 0;
  enum IntegratorType type = BDF2;

  // Parse command line arguments
  for ( int i = 0; i < argc; i++ ){
    // Determine the problem size
    if (sscanf(argv[i], "scale=%d", &scale) == 1){
      if (scale < 0){ scale = 1; }
      if (scale > 10){ scale = 10; }
      if (rank ==0){ printf("Problem scale = %d\n", scale); }
    }
    // Determine the number of functions for adjoint
    if (sscanf(argv[i], "num_funcs=%d", &num_funcs) == 1){
      if (num_funcs < 0){ num_funcs = 1; }
      if (num_funcs > 3){ num_funcs = 3; }
      if (rank == 0){ printf("Number of functions : %d\n", num_funcs); }
    }

    // Determine the number of functions for adjoint
    if (sscanf(argv[i], "num_threads=%d", &num_threads) == 1){
      if (num_threads < 0){ num_threads = 1; }
      if (num_threads > 24){ num_threads = 24; }
      if (rank ==0){ printf("Number of threads : %d\n", num_threads); }
    }

    // Determine whether or not to test gradients with complex step
    if (strcmp("test_gradient", argv[i]) == 0){
      test_gradient = 1;
      if (rank ==0){ printf("Enabled complex-step verification of gradients...\n"); }
    }

    // Determine the type of integration scheme to use
    // Backward Difference Formulae
    if (strcmp("BDF1", argv[i]) == 0){
      type = BDF1;
    } else if (strcmp("BDF2", argv[i]) == 0){
      type = BDF2;
    } else if (strcmp("BDF3", argv[i]) == 0){
      type = BDF3;

      // Adams-Bashforth-Moulton
    } else if (strcmp("ABM1", argv[i]) == 0){
      type = ABM1;
    } else if (strcmp("ABM2", argv[i]) == 0){
      type = ABM2;
    } else if (strcmp("ABM3", argv[i]) == 0){
      type = ABM3;
    } else if (strcmp("ABM4", argv[i]) == 0){
      type = ABM4;
    } else if (strcmp("ABM5", argv[i]) == 0){
      type = ABM5;
    } else if (strcmp("ABM6", argv[i]) == 0){
      type = ABM6;

      // Diagonally Implicit Runge Kutta
    } else if (strcmp("DIRK2", argv[i]) == 0){
      type = DIRK2;
    } else if (strcmp("DIRK3", argv[i]) == 0){
      type = DIRK3;
    } else if (strcmp("DIRK4", argv[i]) == 0){
      type = DIRK4;

      // Newmark-Beta-Gamma method
    } else if (strcmp("NBG", argv[i]) == 0){
      type = NBG;
    }
  }

  int vars_per_node = 8;
  TACSCreator *creator = new TACSCreator(comm, vars_per_node);
  creator->incref();

  if (rank == 0){
    int nl_elems = scale*5; // number of elements along the length
    int nw_elems = scale*1; // number of element along a side of the box

    // Set the number of rigid elements/constraints: 1 rigid body
    // class and 2 constraint classes for each end of the flexible
    // bar.
    int nrigid = 0;
    int num_elements = 4*nl_elems*nw_elems + nrigid;
    
    // Set the total number of nodes 
    int num_nodes = (2*nl_elems+1)*(8*nw_elems) + nrigid;

    // Set the node locations
    TacsScalar *Xpts = new TacsScalar[ 3*num_nodes ];
    
    // Allocate the connectivity
    int conn_size = 9*num_elements;
    int *elem_ptr = new int[ num_elements+1 ];
    int *elem_conn = new int[ conn_size ];
    int *elem_id_nums = new int[ num_elements ];
    memset(elem_id_nums, 0, num_elements*sizeof(int));

    // Set up the connectivity
    int *ptr = elem_ptr;
    int *conn = elem_conn;
    for ( int j = 0; j < 4*nw_elems; j++ ){
      for ( int i = 0; i < nl_elems; i++ ){
        ptr[0] = conn - elem_conn;
        ptr++;

        for ( int jj = 0; jj < 3; jj++ ){
          for ( int ii = 0; ii < 3; ii++ ){          
            if (2*j + jj == 8*nw_elems){
              conn[ii + 3*jj] = 2*i + ii;
            }
            else {
              conn[ii + 3*jj] = 2*i + ii + (2*nl_elems+1)*(2*j + jj);
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
    int nw = 2*nw_elems;
    int nlen = 2*nl_elems+1;
    TacsScalar *x = Xpts;
    for ( int j = 0; j < 4*nw; j++ ){
      for ( int i = 0; i < nlen; i++ ){
        x[0] = (L*i)/(nlen-1);
        if (j < nw){
          x[1] = -0.5*Lw;
          x[2] = -0.5*Lw + (Lw*j)/nw;
        }
        else if (j < 2*nw){
          x[1] = -0.5*Lw + (Lw*(j-nw))/nw;
          x[2] = 0.5*Lw;
        }
        else if (j < 3*nw){
          x[1] = 0.5*Lw;
          x[2] = 0.5*Lw - (Lw*(j-2*nw))/nw;
        }
        else {
          x[1] = 0.5*Lw - (Lw*(j-3*nw))/nw;
          x[2] = -0.5*Lw;
        }
        x += 3;
      }
    }

    int bc = 0;
    int bc_ptr[] = {0, 3};
    int bc_vars[] = {0, 1, 2};
    creator->setBoundaryConditions(1, &bc, bc_ptr, bc_vars);

    // Set the connectivity
    creator->setGlobalConnectivity(num_nodes, num_elements,
  				   elem_ptr, elem_conn,
  				   elem_id_nums);
    
    // Set the nodal locations
    creator->setNodes(Xpts);

    // Free the data
    delete [] Xpts;
    delete [] elem_ptr;
    delete [] elem_conn;
    delete [] elem_id_nums;
  } 

  // Create the objects associated with the rigid-bodies
  // ---------------------------------------------------
  // The acceleration due to gravity in global frame of reference
  TACSGibbsVector *gravVec = new TACSGibbsVector(0.0, 0.0, -9.8);

  // Create the constitutive objects
  TacsScalar rho = 2500.0;
  TacsScalar E = 70e9;
  TacsScalar nu = 0.3;
  TacsScalar kcorr = 5.0/6.0;
  TacsScalar yield_stress = 464.0e6;
  TacsScalar thickness = 0.05;
  
  // Create the stiffness object
  isoFSDTStiffness *stiff = new isoFSDTStiffness(rho, E, nu, kcorr, 
                                                 yield_stress, thickness, 0);
  
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
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(tacs, SHELL, write_flag);

  // Set the number of threads
  tacs->setNumThreads(num_threads);

  // Set up the parameters for adjoint solve
  int num_dvs = tacs->getNumComponents();
  TACSFunction *func[num_funcs];

  if (num_funcs == 1){
    func[0] = new TACSKSFailure(tacs, 100.0); func[0]->incref();
  }
  else if (num_funcs == 2){
    func[0] = new TACSKSFailure(tacs, 100.0); func[0]->incref();
    func[1] = new TACSCompliance(tacs); func[1]->incref();
  } 
  else if (num_funcs == 3){
    func[0] = new TACSKSFailure(tacs, 100.0); func[0]->incref();
    func[1] = new TACSCompliance(tacs); func[1]->incref();
    func[2] = new TACSStructuralMass(tacs); func[2]->incref();
  }
  
  TacsScalar *funcVals     = new TacsScalar[num_funcs]; // adjoint
  TacsScalar *funcValsTmp  = new TacsScalar[num_funcs]; // CSD
  TacsScalar *funcVals1    = new TacsScalar[num_funcs]; // forward/reverse

  TacsScalar *dfdx    = new TacsScalar[num_funcs*num_dvs]; // adjoint
  TacsScalar *dfdx1   = new TacsScalar[num_funcs*num_dvs]; // CSD
  TacsScalar *dfdxTmp = new TacsScalar[num_funcs*num_dvs]; // forward/reverse

  // Create an array of design variables
  TacsScalar *x = new TacsScalar[ num_dvs ]; 
  x[0] = 0.05; 

  // Set up the parameters for the TACSIntegrator
  double tinit = 0.0;
  double tfinal = 0.1;
  int steps_per_second = 1000; 
  TACSIntegrator *obj = TACSIntegrator::getInstance(tacs, tinit, tfinal, 
                                                    steps_per_second, 
                                                    type);
  obj->incref();
  
  // Set optional parameters
  obj->setRelTol(1.0e-8);
  obj->setAbsTol(1.0e-12);
  obj->setMaxNewtonIters(50);
  obj->setPrintLevel(1);
  obj->setJacAssemblyFreq(1);
  obj->setUseLapack(0);
  
  // Set functions of interest
  obj->setFunction(func, num_funcs);
  
  // Adjoint gradient
  double t0 = MPI_Wtime();
  obj->getFuncGrad(num_dvs, x, funcVals, dfdx);
  t0 = MPI_Wtime() - t0;
    
  // Print a summary of time taken
  if (rank == 0){
    obj->printWallTime(t0, 2);

    // Print the adjoint derivative values
    for( int j = 0; j < num_funcs; j++) {
      printf("[%d] Adj NEW  func: %d fval: %15.8e dfdx:", rank, j, RealPart(funcVals[j]));
      for ( int n = 0; n < num_dvs; n++) {
        printf(" %15.8e ",  RealPart(dfdx[n+j*num_dvs]));
      }
      printf("\n");
    }
    printf("\n");
  }

  // Test the adjoint derivatives if sought
  if (test_gradient){

    printf("Finding Complex-step gradient...\n");

    // Complex step verification
    obj->getFDFuncGrad(num_dvs, x, funcValsTmp, dfdxTmp, 1.0e-16);

    if ( rank == 0) { 
      // Print complex step derivative values
      for( int j = 0; j < num_funcs; j++) {
        printf("[%d] CSD      func: %d fval: %15.8e dfdx:", rank, j, RealPart(funcValsTmp[j]));
        for ( int n = 0; n < num_dvs; n++) {
          printf(" %15.8e ",  RealPart(dfdxTmp[n+j*num_dvs]));
        }
        printf("\n");
      }
      printf("\n");

      // Print the differences between complex step and adjoint derivtives
      for ( int j = 0; j < num_funcs; j++ ) {
        printf("[%d] Error Adj NEW  func: %d ferror: %15.8e dfdx error:", rank, j, RealPart(funcValsTmp[j])-RealPart(funcVals[j]) );
        for ( int n = 0; n < num_dvs; n++ ) {
          printf(" %15.8e ",  RealPart(dfdxTmp[j*num_dvs+n]) -  RealPart(dfdx[j*num_dvs+n]) );
        }
        printf("\n");
      }
      printf("\n");
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

  if (num_funcs == 1){
    func[0]->decref();
  }
  else if (num_funcs == 2){
    func[0]->decref();
    func[1]->decref();
  } 
  else if (num_funcs == 3){
    func[0]->decref();
    func[1]->decref();
    func[2]->decref();
  }

  delete [] x;
  delete [] funcVals;
  delete [] funcVals1;
  delete [] funcValsTmp;

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
