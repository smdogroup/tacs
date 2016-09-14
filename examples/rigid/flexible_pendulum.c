#include "TACSCreator.h"
#include "TACSIntegrator.h"
#include "RigidBody.h"
#include "MITC9.h"
#include "isoFSDTStiffness.h"

/*
  Set up a double pendulum with a flexible attachment
*/
int main( int argc, char *argv[] ){
  // Initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  // Initialize the TACSCreator object on all processors
  int rank; 
  MPI_Comm_rank(comm, &rank); 
  int vars_per_node = 8;
  TACSCreator *creator = new TACSCreator(comm, vars_per_node);
  creator->incref();

  if (rank == 0){
    int nl_elems = 5; // number of elements along the length
    int nw_elems = 1; // number of element along a side of the box

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

  // Create the constitutive objects
  TacsScalar rho = 2500.0;
  TacsScalar E = 10e6;
  TacsScalar nu = 0.3;
  TacsScalar kcorr = 5.0/6.0;
  TacsScalar yield_stress = 464.0e6;
  TacsScalar thickness = 0.005;
  
  // Create the stiffness object
  isoFSDTStiffness *stiff = new isoFSDTStiffness(rho, E, nu, kcorr, 
                                                 yield_stress, thickness);
  
  // Create the element class
  TACSElement *elem = new MITC9(stiff, gravVec);

  // This call must occur on all processor
  creator->setElements(&elem, 1);

  // Create the TACS object
  TACSAssembler *tacs = creator->createTACS();
  tacs->incref();

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS);
  TACSToFH5 *f5 = new TACSToFH5(tacs, SHELL, write_flag);
  
  // Set up the parameters for the TACSIntegrator
  double tinit = 0.0;
  double tfinal = 2.0;
  int steps_per_second = 250; 
  int num_stages = 2;
  int max_bdf_order = 2;
  TACSBDFIntegrator *bdf = new TACSBDFIntegrator(tacs, tinit, tfinal,
                                                 steps_per_second, 
                                                 max_bdf_order);
  bdf->incref();
  
  // Set optional parameters
  bdf->setRelTol(1.0e-8);
  bdf->setAbsTol(1.0e-12);
  bdf->setMaxNewtonIters(50);
  bdf->setPrintLevel(1);
  bdf->setJacAssemblyFreq(1);
  bdf->setUseLapack(1);
  bdf->configureOutput(f5, 1, "flexible-pendulum/pendulum_%04d.f5");

  // Integrate and write solution to file
  bdf->integrate();
  bdf->writeSolution("solutionBDF.dat");

  /*
  // Delete objects
  bodyA->decref();
  bodyB->decref();
  conA->decref();
  conB->decref();

  tacs->decref();
  f5->decref();
  bdf->decref();
  */

  creator->decref();
  tacs->decref();
  bdf->decref();

  MPI_Finalize();

  return (0);
}
