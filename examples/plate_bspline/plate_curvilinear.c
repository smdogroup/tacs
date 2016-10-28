#include "PlaneStressBspline.h"
#include "PlaneStressBsplineAll.h"
#include "PlaneStressBsplineStiffness.h"
#include "TACSCreator.h"
#include "TACSToFH5.h"
// Convert patch mesh to control point mesh
void controlPointMesh(int *num_nodes_x, int *num_nodes_y,
                      int *num_nodes,
                      int *num_patches_x, int *num_patches_y,
                      int **ptr, int **conn, int order){
  *num_nodes_x = *num_patches_x+(order-1);
  *num_nodes_y = *num_patches_y+(order-1);
  *num_nodes = (*num_nodes_x)*(*num_nodes_y);
 
  int num_patches = (*num_patches_x)*(*num_patches_y);
  int *_ptr = new int[num_patches+1];
  int *_conn = new int[num_patches*order*order];
  for (int i = 0; i < num_patches+1; i++){
    _ptr[i] = i*order*order;
  }
  int node = 0;
  int k = 0;
  if (order == 3){
    for (int j = 0; j < *num_patches_y; j++){
      for (int i = 0; i < *num_patches_x; i++){
        if (i == 0){
          node = j*(*num_nodes_x);
        }
        _conn[9*k] = node;
        _conn[9*k+1] = node+1;
        _conn[9*k+2] = node+2;
      
        _conn[9*k+3] = node+*num_nodes_x;
        _conn[9*k+4] = node+*num_nodes_x+1;
        _conn[9*k+5] = node+*num_nodes_x+2;
      
        _conn[9*k+6] = node+*num_nodes_x*2;
        _conn[9*k+7] = node+*num_nodes_x*2+1;
        _conn[9*k+8] = node+*num_nodes_x*2+2;
        node++;
        k++;
      }
    }
  }
  else if (order == 4){
    for (int j = 0; j < *num_patches_y; j++){
      for (int i = 0; i < *num_patches_x; i++){
        if (i == 0){
          node = j*(*num_nodes_x);
        }
        _conn[16*k] = node;
        _conn[16*k+1] = node+1;
        _conn[16*k+2] = node+2;
        _conn[16*k+3] = node+3;
        
        _conn[16*k+4] = node+*num_nodes_x;
        _conn[16*k+5] = node+*num_nodes_x+1;
        _conn[16*k+6] = node+*num_nodes_x+2;
        _conn[16*k+7] = node+*num_nodes_x+3;
        
        _conn[16*k+8] = node+*num_nodes_x*2;
        _conn[16*k+9] = node+*num_nodes_x*2+1;
        _conn[16*k+10] = node+*num_nodes_x*2+2;
        _conn[16*k+11] = node+*num_nodes_x*2+3;
        
        _conn[16*k+12] = node+*num_nodes_x*3;
        _conn[16*k+13] = node+*num_nodes_x*3+1;
        _conn[16*k+14] = node+*num_nodes_x*3+2;
        _conn[16*k+15] = node+*num_nodes_x*3+3;
        node++;
        k++;
      }
    }
  }
  *ptr = _ptr;
  *conn = _conn;
}
// Generate the knot vectors for given number of patches in the x and
// y direction
void knotVector(int num_patches_x, int num_patches_y, int order,
                double **Tu, double **Tv){
  int tu = num_patches_x+1+2*(order-1);
  int tv = num_patches_y+1+2*(order-1);

  double *_Tu = new double[tu];
  double *_Tv = new double[tv];
  
  for (int i = 0; i < order-1; i++){
    _Tu[i] = 0.0;
    _Tv[i] = 0.0;
  }
  for (int i = 0; i < num_patches_x; i++){
    _Tu[i+order-1] = 1.0*i;
  }
  for (int i = 0; i < num_patches_y; i++){
    _Tv[i+order-1] = 1.0*i;
  }
  for (int i = num_patches_x; i < tu-(order-1); i++){
    _Tu[i+order-1] = 1.0*num_patches_x;
  } 
  for (int i = num_patches_y; i < tv-(order-1); i++){
    _Tv[i+order-1] = 1.0*num_patches_y;
  } 
  *Tu = _Tu;
  *Tv = _Tv;
}
// Generate the mesh of control points itself based on the four
// external points
void controlXpts(double Xpts[], double **Xpts_c, 
                 int ncp_x, int ncp_y){
  double *_Xpts = new double[3*ncp_x*ncp_y];
  int edge_index[] = {0, 1,
                      0, 2,
                      1, 3, 
                      2, 3};
  double edge_pt[3*(2*ncp_x+2*ncp_y)];
  int ind = 0;
  // Compute the control points on the exterior boundary
  for (int i = 0; i < 4; i++){
   
    double a[] = {Xpts[3*edge_index[2*i]], Xpts[3*edge_index[2*i]+1]};
    double b[] = {Xpts[3*edge_index[2*i+1]],Xpts[3*edge_index[2*i+1]+1]};
    
    double ix = (b[0]-a[0])/(ncp_x-1);
    double iy = (b[1]-a[1])/(ncp_y-1);
    
    if (i == 1 || i == 2){
      for (int j = 0; j < ncp_y; j++){
        edge_pt[3*ind] = a[0]+j*ix;
        edge_pt[3*ind+1] = a[1]+j*iy;
        edge_pt[3*ind+2] = 0.0;
        /* printf("edge: %d, %e %e \n", i, edge_pt[3*ind], edge_pt[3*ind+1]); */
        ind++;      
      }
    }
    else {
      for (int j = 0; j < ncp_x; j++){
        edge_pt[3*ind] = a[0]+j*ix;
        edge_pt[3*ind+1] = a[1]+j*iy;
        edge_pt[3*ind+2] = 0.0;
        /* printf("edge: %d, %e %e \n", i, edge_pt[3*ind], edge_pt[3*ind+1]); */
        ind++;      
      }
    }
  }
 
  // Reset the counter
  ind = 0;
  // Initialize the control points
  for (int i = 0; i < ncp_x; i++){
    _Xpts[3*ind] = edge_pt[3*ind];
    _Xpts[3*ind+1] = edge_pt[3*ind+1];
    _Xpts[3*ind+2] = edge_pt[3*ind+2];
    ind++;
  }
  
  // Compute the internal control points
  for (int i = 0; i < ncp_y-2; i++){
    double a[] = {edge_pt[3*(ncp_x+i+1)],
                  edge_pt[3*(ncp_x+i+1)+1]};
    double b[] = {edge_pt[3*(2*ncp_x+i+1)],
                  edge_pt[3*(2*ncp_x+i+1)+1]};
    /* printf("%d a: %e %e \n", i, a[0], a[1]); */
    /* printf("%d b: %e %e \n", i, b[0], b[1]); */
    double ix = (b[0]-a[0])/(ncp_x-1);
    double iy = (b[1]-a[1])/(ncp_y-1);
    for (int j = 0; j < ncp_x; j++){
      _Xpts[3*ind] = a[0]+j*ix;
      _Xpts[3*ind+1] = a[1]+j*iy;
      _Xpts[3*ind+2] = 0.0;
      ind++;
      /* printf("%d Xpts: %e %e \n", i,  */
      /*        _Xpts[3*ind], _Xpts[3*ind+1]); */
    }
  }
  // Input the last exterior control points
  for (int i = 0; i < ncp_x; i++){
    _Xpts[3*ind] = edge_pt[3*(ncp_x+2*ncp_y+i)+0];
    _Xpts[3*ind+1] = edge_pt[3*(ncp_x+2*ncp_y+i)+1];
    _Xpts[3*ind+2] = 0.0;
    /* printf("%d Xpts: %e %e \n", i,  */
    /*        _Xpts[3*ind], _Xpts[3*ind+1]); */
    ind++;
  }
 
  /* for (int i = 0; i < ind; i++){ */
  /*   printf("%d Xpts: %e %e %e \n", i, _Xpts[3*i], */
  /*          _Xpts[3*i+1], _Xpts[3*i+2]); */
  /* } */
  *Xpts_c = _Xpts;  
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  
  // Get the rank of the processor
  int rank,size; 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  // Allocate the TACS creator
  TACSCreator *creator = new TACSCreator(MPI_COMM_WORLD, 2);
  creator->incref();
  
  // Number of knot intervals without the repeating knots
  int Lu = 50, Lv = 50;
  int ind = 0;
  // Order of the Bspline
  int order = 3;
  for (int k = 0; k < argc; k++){
    int _order, _Lu, _Lv;
    if (sscanf(argv[k], "order=%d", &_order) == 1){
      order = _order;
    }
    if (sscanf(argv[k], "Lu=%d", &_Lu) == 1){
      Lu = _Lu;
    }
    if (sscanf(argv[k], "Lv=%d", &_Lv) == 1){
      Lv = _Lv;
    }
  }
  printf("Order is %d \n", order);
  printf("Nx: %d \n", Lu);
  printf("Ny: %d \n", Lv);

  TACSElement **elem = new TACSElement*[Lu*Lv];
  
  double *Tu, *Tv;
  knotVector(Lu, Lv,order, &Tu, &Tv);
  
  /* for (int i = 0; i < Lu+1+2*(order-1); i++){ */
  /*   printf("Tu[%d]: %e\n", i,Tu[i]); */
  /* } */
  TacsScalar *x = new TacsScalar[Lu*Lv];
  memset(x, 0.5, Lu*Lv*sizeof(TacsScalar));
  if (order == 4){
    for (int j = 0; j < Lv; j++){
      for (int i = 0; i < Lu; i++){
        // Create the stiffness object
        PlaneStressBsplineStiffness *stiff = 
          new PlaneStressBsplineStiffness(2700.0, 70.0e9,
                                          0.3, 280e6, 0.0,
                                          x, 0.0, 1e-3,
                                          Tu, Tv, Lu, Lv,
                                          ind, order);
        stiff->incref();
                                          
        elem[ind] = new PlaneStressBsplineAll<4>(stiff,Tu,Tv,
                                                 Lu, Lv,
                                                 LINEAR,
                                                 0, ind);
        stiff->decref();
        elem[ind]->incref();
        ind++;
      }
    }
  }
  else if (order == 3){
    for (int j = 0; j < Lv; j++){
      for (int i = 0; i < Lu; i++){
        // Create the stiffness object
        /* PlaneStressBsplineStiffness *stiff =  */
        /*   new PlaneStressBsplineStiffness(2700.0, 70.0e9, */
        /*                                   0.3, 280e6, 0.0, */
        /*                                   x, 0.0, 1e-3, */
        /*                                   Tu, Tv, Lu, Lv, */
        /*                                   ind, order); */
        PlaneStressStiffness *stiff = 
          new PlaneStressStiffness(2700.0,70.0e9, 0.3);
        stiff->incref();
        
        /* elem[ind] = new PlaneStressBspline(stiff,Tu,Tv, */
        /*                                    Lu, Lv, */
        /*                                    LINEAR, */
        /*                                    0, ind);  */
        elem[ind] = new PlaneStressBsplineAll<3>(stiff,Tu,Tv,
                                                 Lu, Lv,
                                                 LINEAR,
                                                 0, ind);
        stiff->decref();
        elem[ind]->incref();
        ind++;
      }
    }    
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("rank: %d \n", rank);
  // Create mesh on root processor
  if (rank == 0){
    int num_nodes, num_nodes_x, num_nodes_y;
    int num_elems_x = 1*Lu, num_elems_y = 1*Lv;
    int num_elems = num_elems_x*num_elems_y;
    int *ptr_c, *conn_c;
    
    int *elem_ids = new int[Lu*Lv];
    for (int i = 0; i < Lu*Lv; i++){
      elem_ids[i] = i*1;
    }
    
    controlPointMesh(&num_nodes_x, &num_nodes_y,
                     &num_nodes, &num_elems_x,
                     &num_elems_y,
                     &ptr_c, &conn_c, order);
    /* double Xpts[] = {0.0, 0.0, 0.0,  */
    /*                  2.0, -1.0, 0.0, */
    /*                  1.0, 2.0, 0.0, */
    /*                  4.0, 2.0, 0.0}; */
    double Xpts[] = {0.0,5.0, 0.0,
                     4.0,5.0, 0.0, 
                     1.0, 9.0, 0.0, 
                     5.0, 9.0, 0.0};
    printf("num_nodes: %d, num_elems: %d %d %d\n", num_nodes, num_elems,
           Lu, Lv);
    creator->setGlobalConnectivity(num_nodes, num_elems, 
                                   ptr_c, conn_c,elem_ids);
    
    int num_bcs_c = num_nodes_y*1;
    int bc_nodes_c[num_bcs_c];
    for (int i = 0; i < num_bcs_c; i++){
      bc_nodes_c[i] = i*num_nodes_y;
    }    
    // Set the boundary conditions
    creator->setBoundaryConditions(num_bcs_c, bc_nodes_c);
    if (order == 4){
      double *Xpts_c;
      controlXpts(Xpts, &Xpts_c, Lu+order-1, Lv+order-1);
      /* TacsScalar Xpts_c[] = {0.0, 0.0, 0.0, */
      /*                        0.667, 0.0, 0.0, */
      /*                        1.333, 0.0, 0.0, */
      /*                        2.0, 0.0, 0.0, */
      /*                        2.667, 0.0, 0.0, */
      /*                        3.333, 0.0, 0.0, */
      /*                        4.0, 0.0, 0.0, */
      /*                        0.0, 0.667, 0.0, */
      /*                        0.667, 0.667, 0.0, */
      /*                        1.333, 0.667, 0.0, */
      /*                        2.0, 0.667, 0.0, */
      /*                        2.667, 0.667, 0.0, */
      /*                        3.333, 0.667, 0.0, */
      /*                        4.0, 0.667, 0.0, */
      /*                        0.0, 1.333, 0.0, */
      /*                        0.667, 1.333, 0.0, */
      /*                        1.333, 1.333, 0.0, */
      /*                        2.0, 1.333, 0.0, */
      /*                        2.667, 1.333, 0.0, */
      /*                        3.333, 1.333, 0.0, */
      /*                        4.0, 1.333, 0.0, */
      /*                        0.0, 2.0, 0.0, */
      /*                        0.667, 2.0, 0.0, */
      /*                        1.333, 2.0, 0.0, */
      /*                        2.0, 2.0, 0.0, */
      /*                        2.667, 2.0, 0.0, */
      /*                        3.333, 2.0, 0.0, */
      /*                        4.0, 2.0, 0.0, */
      /*                        0.0, 2.667, 0.0, */
      /*                        0.667, 2.667, 0.0, */
      /*                        1.333, 2.667, 0.0, */
      /*                        2.0, 2.667, 0.0, */
      /*                        2.667, 2.667, 0.0, */
      /*                        3.333, 2.667, 0.0, */
      /*                        4.0, 2.667, 0.0, */
      /*                        0.0, 3.333, 0.0, */
      /*                        0.667, 3.333, 0.0, */
      /*                        1.333, 3.333, 0.0, */
      /*                        2.0, 3.333, 0.0, */
      /*                        2.667, 3.333, 0.0, */
      /*                        3.333, 3.333, 0.0, */
      /*                        4.0, 3.333, 0.0, */
      /*                        0.0, 4.0, 0.0, */
      /*                        0.667, 4.0, 0.0, */
      /*                        1.333, 4.0, 0.0, */
      /*                        2.0, 4.0, 0.0, */
      /*                        2.667, 4.0, 0.0, */
      /*                        3.333, 4.0, 0.0, */
      /*                        4.0, 4.0, 0.0}; */
     
      // Set the nodal locations
      creator->setNodes(Xpts_c);
    }
    if (order == 3){
      double *Xpts_d;
      controlXpts(Xpts, &Xpts_d, Lu+order-1, Lv+order-1);
      
      TacsScalar Xpts_c[] = {0.0, 5.0, 0.0,
                             0.8, 5.0, 0.0,
                             1.6, 5.0, 0.0,
                             2.4, 5.0, 0.0,
                             3.2, 5.0, 0.0,
                             4.0, 5.0, 0.0,
                             0.2, 5.8, 0.0,
                             1.0, 5.8, 0.0,
                             1.8, 5.8, 0.0,
                             2.6, 5.8, 0.0,
                             3.4, 5.8, 0.0,
                             4.2, 5.8, 0.0,
                             0.4, 6.6, 0.0,
                             1.2, 6.6, 0.0,
                             2.0, 6.6, 0.0,
                             2.8, 6.6, 0.0,
                             3.6, 6.6, 0.0,
                             4.4, 6.6, 0.0,
                             0.6, 7.4, 0.0,
                             1.4, 7.4, 0.0,
                             2.2, 7.4, 0.0,
                             3.0, 7.4, 0.0,
                             3.8, 7.4, 0.0,
                             4.6, 7.4, 0.0,
                             0.8, 8.2, 0.0,
                             1.6, 8.2, 0.0,
                             2.4, 8.2, 0.0,
                             3.2, 8.2, 0.0,
                             4.0, 8.2, 0.0,
                             4.8, 8.2, 0.0,
                             1.0, 9.0, 0.0,
                             1.8, 9.0, 0.0,
                             2.6, 9.0, 0.0,
                             3.4, 9.0, 0.0,
                             4.2, 9.0, 0.0,
                             5.0, 9.0, 0.0};
      /* for (int i = 0; i < (Lu+(order-1))*(Lv+(order-1)); i++){ */
      /*   /\* if (Xpts_d[3*i+2] != 0.0){ *\/ */
      /*   /\*   printf("i: %d \n", i); *\/ */
      /*   /\*   Xpts_d[3*i+2] = 0.0; *\/ */
      /*   /\* } *\/ */
      /*   /\* printf("Xpts[%d]: %e %e %e \n", i, Xpts_c[3*i], Xpts_c[3*i+1], *\/ */
      /*   /\*        Xpts_c[3*i+2]); *\/ */
      /*   printf("dXpts[%d]: %e %e %e \n", i, Xpts_d[3*i], Xpts_d[3*i+1], */
      /*          Xpts_d[3*i+2]); */
      /*   /\* printf("%d %d %d \n", Xpts_c[3*i]==Xpts_d[3*i],  *\/ */
      /*   /\*        Xpts_c[3*i+1]==Xpts_d[3*i+1],Xpts_c[3*i+2]==Xpts_d[3*i+2]); *\/ */
      /* } */
      // Set the nodal locations
      creator->setNodes(Xpts_d);
    }
    delete [] elem_ids;
  }
  
  // This call must occur on all processor
  creator->setElements(&elem[0], Lu*Lv);

  // Set the reordering typr
  creator->setReorderingType(TACSAssembler::NATURAL_ORDER,
                             TACSAssembler::APPROXIMATE_SCHUR);

  // Create the TACSAssembler object
  TACSAssembler *tacs = creator->createTACS();
  tacs->incref();
  // Test the element
  /* tacs->testElement(0,2); */
  /* tacs->decref(); */
  /* exit(0); */
  // Test the constitutive class
  /* tacs->testConstitutive(0,2); */
  /* tacs->decref(); */
  /* exit(0); */
  // Create the preconditioner
  TACSBVec *res = tacs->createVec();
  TACSBVec *ans = tacs->createVec();
  FEMat *mat = tacs->createFEMat();
  
  // Increment the reference count to the matrix/vectors
  res->incref();
  ans->incref();
  mat->incref();

  // Allocate the factorization
  int lev = 4500;
  double fill = 10.0;
  int reorder_schur = 1;
  PcScMat *pc = new PcScMat(mat, lev, fill, reorder_schur);
  pc->incref();

  // Assemble and factor the stiffness/Jacobian matrix
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  tacs->assembleJacobian(alpha, beta, gamma, res, mat);
  
  mat->applyBCs();
  pc->factor();
  
  // Number of GMRES iterations
  int gmres_iters = 10;
  // Number of allowed restartions
  int nrestart = 2;
  // Is it a flexible preconditioner
  int is_flexible = 1;
  GMRES *gmres = new GMRES(mat, pc, gmres_iters,
                           nrestart, is_flexible);
  gmres->incref();
  gmres->setTolerances(1e-12, 1e-30);
  gmres->setMonitor(new KSMPrintStdout("GMRES", 0, 1));

  TacsScalar *res_array;
  int res_size = res->getArray(&res_array);
  res_array[res_size-1] = 1.0;
  res->applyBCs();
  gmres->solve(res, ans);
  tacs->setVariables(ans);

  // Constant shear strain
  /* TacsScalar *X, *a; */
  /* TACSBVec *Xvec = tacs->createNodeVec(); */
  /* Xvec->incref(); */
  /* tacs->getNodes(Xvec); */

  /* Xvec->getArray(&X); */
  /* int asize = ans->getArray(&a)/2;  */

  /* for ( int i = 0; i < asize; i++ ){ */
  /*   a[2*i] = X[3*i+1]; */
  /*   a[2*i+1] = X[3*i]; */
  /* } */
  /* tacs->setVariables(ans); */

  // Create TACSToFH5 object for writing to tecplot
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(tacs, PLANE_STRESS, write_flag);
  f5->incref();
  char filename[256];
  sprintf(filename, "plate_curvilinear_order%d_%dx%d.f5",order, Lu, Lv);
  f5->writeToFile(filename);
  
  // Free everything
  f5->decref();
  
  // Decrease the reference count to the linear algebra objects
  gmres->decref();
  pc->decref();
  mat->decref();
  ans->decref();
  res->decref();
  // Decrease the reference count to everything else
  if (elem){
    for (int i = 0; i < Lu*Lv; i++){
      elem[i]->decref();
    }
  }
  creator->decref();  
  MPI_Finalize();
  return 0;
}
