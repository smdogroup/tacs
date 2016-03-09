#include "PlaneStressTri6.h"
#include "TACSCreator.h"
#include "TACSToFH5.h"

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);

  // Create a regular mesh of triangular elements
  int nx = 10, ny = 10;
  int num_elements = 2*nx*ny;
  int num_nodes = (2*nx + 1)*(2*ny + 1);

  // Set the lengths along the x and y directions
  double Lx = 10.0, Ly = 10.0;
  
  // Create the nodes
  TacsScalar *Xpts = new TacsScalar[ 3*num_nodes ];
  for ( int j = 0; j < 2*ny+1; j++ ){
    for ( int i = 0; i < 2*nx+1; i++ ){
      Xpts[3*(i + j*(2*nx+1))] = 0.5*Lx*i/nx;
      Xpts[3*(i + j*(2*nx+1))+1] = 0.5*Ly*j/ny;
      Xpts[3*(i + j*(2*nx+1))+2] = 0.0;
    }
  }

  // Set up the element connectivity arrays
  int *elem_node_conn = new int[ 6*num_elements ];
  int *elem_node_ptr = new int[ num_elements+1 ];

  int n = 0;
  int *conn = elem_node_conn;
  elem_node_ptr[0] = 0;
  for ( int j = 0; j < ny; j++ ){
    for ( int i = 0; i < nx; i++ ){
      conn[0] = 2*i + (2*nx+1)*(2*j);
      conn[1] = 2*i+2 + (2*nx+1)*(2*j);
      conn[2] = 2*i + (2*nx+1)*(2*j+2);
      conn[3] = 2*i+1 + (2*nx+1)*(2*j);
      conn[4] = 2*i+1 + (2*nx+1)*(2*j+1);
      conn[5] = 2*i + (2*nx+1)*(2*j+1);
      conn += 6;
      elem_node_ptr[n+1] = elem_node_ptr[n] + 6;
      n++;

      conn[0] = 2*i+2 + (2*nx+1)*(2*j+2);
      conn[1] = 2*i + (2*nx+1)*(2*j+2);
      conn[2] = 2*i+2 + (2*nx+1)*(2*j);
      conn[3] = 2*i+1 + (2*nx+1)*(2*j+2);
      conn[4] = 2*i+1 + (2*nx+1)*(2*j+1);
      conn[5] = 2*i+2 + (2*nx+1)*(2*j+1);
      conn += 6;
      elem_node_ptr[n+1] = elem_node_ptr[n] + 6;
      n++;
    }
  }

  // Set the identity numbers
  int *elem_id_nums = new int[ num_elements ];
  memset(elem_id_nums, 0, num_elements*sizeof(int));

  // Set the boundary conditions
  int num_bcs = 2*ny+1;
  int *bc_nodes = new int[ num_bcs ];
  for ( int j = 0; j < 2*ny+1; j++ ){
    bc_nodes[j] = j*(2*nx+1);
  }

  // Allocate the TACS creator
  TACSCreator *creator = new TACSCreator(MPI_COMM_WORLD, 2);
  creator->incref();
  
  // Set the connectivity
  creator->setGlobalConnectivity(num_nodes, num_elements,
				 elem_node_ptr, elem_node_conn,
				 elem_id_nums);

  // Set the boundary conditions
  creator->setBoundaryConditions(num_bcs, bc_nodes, NULL, NULL);
  
  // Set the nodal locations
  creator->setNodes(Xpts);

  // Create the stiffness object
  PlaneStressStiffness *stiff = new PlaneStressStiffness(1.0, 70.0, 0.3);
  stiff->incref();

  // Create the plane stress 
  TACSElement *elem = new PlaneStressTri6(stiff);
  elem->incref();

  creator->setElements(&elem, 1);

  // Create the TACSAssembler object
  TACSAssembler *tacs = creator->createTACS();
  tacs->incref();

  // Create the preconditioner
  BVec *res = tacs->createVec();
  BVec *ans = tacs->createVec();
  FEMat *mat = tacs->createFEMat();
  int lev = 4500;
  double fill = 10.0;
  int reorder_schur = 1;
  PcScMat *pc = new PcScMat(mat, lev, fill, reorder_schur);
  pc->incref();

  // Assemble and factor the stiffness/Jacobian matrix
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  tacs->assembleJacobian(res, mat, alpha, beta, gamma);			 
  pc->factor();

  res->set(1.0);
  res->applyBCs();
  pc->applyFactor(res, ans);
  tacs->setVariables(ans);

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 * f5 = new TACSToFH5(tacs, PLANE_STRESS, write_flag);
  f5->incref();
  f5->writeToFile("triangle.f5");

  f5->decref();
  stiff->decref();
  elem->decref();
  tacs->decref();
  creator->decref();

  MPI_Finalize();
  return (0);
}
