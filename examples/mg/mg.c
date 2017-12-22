#include "TACSAssembler.h"
#include "MITCShell.h"
#include "TACSShellTraction.h"
#include "isoFSDTStiffness.h"
#include "TACSMg.h"
#include "TACSCreator.h"
#include "TACSToFH5.h"
#include "TACSBuckling.h"

/*
  Create the TACSAssembler object and return the associated TACS
  creator object 
*/
void createTACS( MPI_Comm comm, int nx, int ny,
                 TACSAssembler **_tacs, TACSCreator **_creator ){
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Set the number of nodes/elements on this proc
  int varsPerNode = 6;

  // Set up the creator object
  TACSCreator *creator = new TACSCreator(comm, varsPerNode);

  if (rank == 0){
    // Set the number of elements
    int numNodes = (nx+1)*(ny+1);
    int numElements = nx*ny;

    // Allocate the input arrays into the creator object
    int *ids = new int[ numElements ];
    int *ptr = new int[ numElements+1 ];
    int *conn = new int[ 4*numElements ];

    // Use just one element
    memset(ids, 0, numElements*sizeof(int));
    
    ptr[0] = 0;
    for ( int k = 0; k < numElements; k++ ){
      // Back out the i, j coordinates from the corresponding
      // element number
      int i = k % nx;
      int j = k/nx;

      // Set the node connectivity
      conn[4*k] = i + j*(nx+1);
      conn[4*k+1] = i+1 + j*(nx+1);
      conn[4*k+2] = i + (j+1)*(nx+1);
      conn[4*k+3] = i+1 + (j+1)*(nx+1);
      ptr[k+1] = 4*(k+1);
    }
  
    // Set the connectivity
    creator->setGlobalConnectivity(numNodes, numElements,
                                   ptr, conn, ids);
    delete [] conn;
    delete [] ptr;
    delete [] ids;

    // We're over-counting one of the nodes on each edge
    int numBcs = 4*(nx+1);
    int *bcNodes = new int[ numBcs ]; 

    for ( int i = 0; i < (nx+1); i++ ){
      bcNodes[4*i] = i;
      bcNodes[4*i+1] = i + (nx+1)*ny;
      bcNodes[4*i+2] = i*(nx+1);
      bcNodes[4*i+3] = (i+1)*(nx+1)-1;
    }

    // Set the boundary conditions
    creator->setBoundaryConditions(numBcs, bcNodes);
    delete [] bcNodes;

    // Set the node locations
    TacsScalar *Xpts = new TacsScalar[ 3*numNodes ];
    for ( int j = 0; j < ny+1; j++ ){
      for ( int i = 0; i < nx+1; i++ ){
        int node = i + (nx+1)*j;
        Xpts[3*node] = (1.0*i)/nx;
        Xpts[3*node+1] = (1.0*j)/ny;
        Xpts[3*node+2] = 0.0;
      }
    }

    // Set the nodal locations
    creator->setNodes(Xpts);
    delete [] Xpts;
  }

  // Create and set the element
  TacsScalar rho = 2500.0; // Not used
  TacsScalar E = 70e9;
  TacsScalar nu = 0.3;
  TacsScalar kcorr = 5.0/6.0; // The shear correction factor
  TacsScalar yield_stress = 464.0e6;
  TacsScalar thickness = 0.005;
    
  // Create the stiffness object
  isoFSDTStiffness *stiff = new isoFSDTStiffness(rho, E, nu, kcorr, 
                                                 yield_stress, thickness);
  TACSElement *elem = new MITCShell<2>(stiff);
  creator->setElements(&elem, 1);

  creator->setReorderingType(TACSAssembler::MULTICOLOR_ORDER,
                             TACSAssembler::GAUSS_SEIDEL);

  // Create TACS
  TACSAssembler *tacs = creator->createTACS();

  // Get all the element numbers
  int *elemNums;
  int id = 0;
  int numElems = creator->getElementIdNums(&id, 1, &elemNums);
  
  // Create the auxiliary element class
  TACSAuxElements *aux = new TACSAuxElements(tacs->getNumElements());

  // Create a surface traction associated with this element and add
  // it to the auxilary elements. Note that the element number must
  // correspond to the local element number used for this processor.
  TacsScalar tx = 0.0, ty = 0.0, tz = -1e5; 
  TACSShellTraction<2> *trac = new TACSShellTraction<2>(tx, ty, tz);
  for ( int i = 0; i < numElems; i++ ){
    aux->addElement(elemNums[i], trac);
  }
  delete [] elemNums;

  // Set the auxiliary elements into TACSAssembler
  tacs->setAuxElements(aux);

  // Set the pointers
  *_tacs = tacs;
  *_creator = creator;
}

/*
  The following code illustrates the use of the geometric multigrid
  capabilities in TACS using a flat plate example.

  The main work required is the assembly of the inter-grid operators
  themselves. The TACSBVecInterp class is used to simplify the set up
  of these operators which are stored explicitly.
  
  The multigrid operations are performed using the TACSMg object which
  inherits from the TACSPc interface. It can be used directly as a
  preconditioner and needs to be "factored" in the same manner as
  regular preconditioners. 

  TACSMg has a number of member functions which facilitate setting
  state and design variables on all multigrid levels and matrix
  assembly.  Matrix assembly on all levels can be performed by calling
  the assembleJacobian() and the assembleMatType() functions directly
  in the TACSMg interface.
*/
int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);

  // Get the rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Number of different levels
  int nlevels = 3;
  const int max_nlevels = 5;
  TACSAssembler *tacs[max_nlevels];
  TACSCreator *creator[max_nlevels];

  // Set whether to use the direct solver or not
  int use_direct_lanczos = 0;

  // Set the dimension of the largest meshes
  int nx = 128;
  int ny = 128;

  // Get the global size of the mesh from the input
  for ( int k = 0; k < argc; k++ ){
    int xpow;
    if (sscanf(argv[k], "xpow=%d", &xpow) == 1){
      if (xpow < 5){ xpow = 5; }
      if (xpow > 10){ xpow = 10; }
      nx = 1 << xpow;
      ny = 1 << xpow;
    }
    if (sscanf(argv[k], "nlevels=%d", &nlevels) == 1){
      if (nlevels < 2){ nlevels = 2; }
      if (nlevels > max_nlevels){ nlevels = max_nlevels; }
    }
    if (strcmp(argv[k], "use_direct_lanczos") == 0){
      use_direct_lanczos = 1;
    }
  }

  // Create the multigrid object
  double omega = 1.0;
  int sor_iters = 1;
  int sor_symm = 0;
  TACSMg *mg = new TACSMg(comm, nlevels, omega, sor_iters, sor_symm);
  mg->incref();

  // Create the TACS/Creator objects for all levels
  for ( int i = 0; i < nlevels; i++ ){
    int Nx = nx/(1 << i), Ny = ny/(1 << i);
    createTACS(comm, Nx, Ny, &tacs[i], &creator[i]);
    tacs[i]->incref();
    creator[i]->incref();
  }

  // Create the interpolation operators
  for ( int level = 0; level < nlevels-1; level++ ){
    // Allocate the interpolation object
    TACSBVecInterp *interp = 
      new TACSBVecInterp(tacs[level+1], tacs[level]);

    if (rank == 0){
      // Retrieve the node numbers
      const int *nodes, *coarse_nodes; 
      creator[level]->getNodeNums(&nodes);
      creator[level+1]->getNodeNums(&coarse_nodes);

      // Set the weights and variables for 
      // each node in the fine mesh
      int vars[4];
      TacsScalar w[4];

      // Compute the number of nodes on the fine and coarse
      // meshes
      int nx_fine = nx/(1 << level);
      int ny_fine = ny/(1 << level);
      int nx_coarse = nx/(1 << (level+1));

      // In this case, we cheat and set the entire interpolation
      // operator just from the root processor. The operator will be
      // communicated to the appropriate procs during initialization
      // on all procs.
      for ( int j = 0; j < ny_fine+1; j++ ){
        for ( int i = 0; i < nx_fine+1; i++ ){
          int nvars = 0;
          if ((i % 2 == 0) && (j % 2 == 0)){
            w[0] = 1.0;
            vars[0] = coarse_nodes[(i/2) + (j/2)*(nx_coarse+1)];
            nvars = 1;
          }
          else if (i % 2 == 0){
            w[0] = w[1] = 0.5;
            vars[0] = coarse_nodes[(i/2) + (j/2)*(nx_coarse+1)]; 
            vars[1] = coarse_nodes[(i/2) + (j/2+1)*(nx_coarse+1)];
            nvars = 2;
          }
          else if (j % 2 == 0){
            w[0] = w[1] = 0.5;
            vars[0] = coarse_nodes[(i/2) + (j/2)*(nx_coarse+1)]; 
            vars[1] = coarse_nodes[(i/2+1) + (j/2)*(nx_coarse+1)];
            nvars = 2;
          }
          else {
            w[0] = w[1] = w[2] = w[3] = 0.25;
            vars[0] = coarse_nodes[(i/2) + (j/2)*(nx_coarse+1)]; 
            vars[1] = coarse_nodes[(i/2+1) + (j/2)*(nx_coarse+1)];
            vars[2] = coarse_nodes[(i/2) + (j/2+1)*(nx_coarse+1)]; 
            vars[3] = coarse_nodes[(i/2+1) + (j/2+1)*(nx_coarse+1)];
            nvars = 4;
          }

          int node = nodes[i + (nx_fine+1)*j];
          interp->addInterp(node, w, vars, nvars);
        }
      }
    }

    // Initialize the interpolation object. This is a collective
    // call that distributes the interpolation operator.
    interp->initialize();
    
    // Set the multigrid information at this level
    mg->setLevel(level, tacs[level], interp, 1);
  }

  // Set the model at the lowest grid level
  mg->setLevel(nlevels-1, tacs[nlevels-1], NULL);

  // We no longer require any of the creator objects
  for ( int i = 0; i < nlevels; i++ ){
    creator[i]->decref();
  }

  // Create the residual and solution vectors on the finest TACS mesh
  TACSBVec *res = tacs[0]->createVec();  res->incref();
  TACSBVec *ans = tacs[0]->createVec();  ans->incref();

  // Allocate the GMRES solution method
  int gmres_iters = 100;
  int nrestart = 1;
  int is_flexible = 0;
  GMRES *gmres = new GMRES(mg->getMat(0), mg,
                           gmres_iters, nrestart, is_flexible);
  gmres->incref();

  // Set a monitor to check on solution progress
  int freq = 5;
  gmres->setMonitor(new KSMPrintStdout("GMRES", rank, freq));

  // The initial time
  double t0 = MPI_Wtime();

  // Assemble the Jacobian matrix for each level
  mg->assembleJacobian(1.0, 0.0, 0.0, res);

  // "Factor" the preconditioner
  mg->factor();

  mg->setMonitor(new KSMPrintStdout("MG", rank, freq));
  mg->solve(res, ans);

  // Compute the solution using GMRES
  gmres->solve(res, ans);

  t0 = MPI_Wtime() - t0;
  
  // Set the variables into TACS
  ans->scale(-1.0);
  tacs[0]->setVariables(ans);

  tacs[0]->assembleRes(res);
  TacsScalar res_norm = res->norm();
  if (rank == 0){
    printf("||R||: %15.5e\n", TacsRealPart(res_norm));
    printf("Solution time: %e\n", t0);
  }

  // Output for visualization
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(tacs[0], TACS_SHELL, write_flag);
  f5->incref();
  f5->writeToFile("mg_plate.f5");
  
  // Perform the Frequency Analysis
  int max_lanczos = 100;
  int num_eigvals = 20;
  double eigtol = 1e-12;
  TacsScalar sigma = 0.0;

  TACSFrequencyAnalysis *freq_analysis = NULL;

  if (use_direct_lanczos){
    // Create matrix and vectors
    FEMat *kmat = tacs[0]->createFEMat(); // stiffness matrix
    kmat->incref();
    
    FEMat *mmat = tacs[0]->createFEMat();
    mmat->incref();
    
    // Create the preconditioner
    int lev_fill = 10000;
    double fill = 10.0;
    PcScMat *pc = new PcScMat(kmat, lev_fill, fill, 1);
    pc->incref();
    
    GMRES *ksm = new GMRES(kmat, pc, gmres_iters, nrestart, is_flexible);
    ksm->incref();
    ksm->setTolerances(1e-12, 1e-30);

    freq_analysis =
      new TACSFrequencyAnalysis(tacs[0], sigma, mmat, kmat, ksm, 
                                max_lanczos, num_eigvals, eigtol);
    freq_analysis->incref();

    // Decref the pointers that go out of scope
    kmat->decref();
    mmat->decref();
    pc->decref();
  }
  else {
    // Create a mass matrix
    TACSMat *mmat = tacs[0]->createMat();
    mmat->incref();
    
    // Create the frequency analysis object
    freq_analysis = 
      new TACSFrequencyAnalysis(tacs[0], sigma, mmat, 
                                mg->getMat(0), gmres,
                                max_lanczos, num_eigvals, eigtol);
    freq_analysis->incref();
    
    // Decref the poitners that go out of scope
    mmat->decref();
  }
 
  // Solve the frequency analysis problem
  freq_analysis->solve(new KSMPrintStdout("Frequency analysis", rank, freq));

  // Print out the eigenvalues and natural frequencies
  for ( int k = 0; k < num_eigvals; k++ ){
    TacsScalar error;
    TacsScalar eigvalue = 
      freq_analysis->extractEigenvector(k, ans, &error);

    tacs[0]->setVariables(ans);
    char filename[128];
    sprintf(filename, "eigenvector%d.f5", k);
    f5->writeToFile(filename);

    if (rank == 0){
      printf("TACS eigs[%2d]: %15.6f Omega = %15.6f [Hz]\n", 
	     k, TacsRealPart(eigvalue), 
             TacsRealPart(sqrt(eigvalue))/(2.0*3.14159));
    }
  }

  // Free the memory
  freq_analysis->decref();
  f5->decref();
  ans->decref();
  res->decref();
  gmres->decref();
  for ( int i = 0; i < nlevels; i++ ){
    tacs[i]->decref();
  }

  MPI_Finalize();
  return (0);
}
