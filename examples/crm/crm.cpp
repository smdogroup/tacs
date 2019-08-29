#include "TACSMeshLoader.h"
#include "MITCShell.h"
#include "KSFailure.h"
#include "isoFSDTStiffness.h"

int main( int argc, char **argv ){
  // Intialize MPI and declare communicator
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int rank;
  MPI_Comm_rank(comm, &rank);

#ifdef TACS_USE_COMPLEX
  double dh = 1e-30;
#else
  double dh = 1e-6;
#endif // TACS_USE_COMPLEX

  // Scan the arguments to check for a manual step size
  for ( int k = 0; k < argc; k++ ){
    if (sscanf(argv[k], "dh=%lf", &dh) == 1){
      if (rank == 0){
        printf("Using step size dh = %e\n", dh);
      }
    }
  }

  // Write name of BDF file to be load to char array
  const char *filename = "CRM_box_2nd.bdf";

  // Create the mesh loader object and load file
  TACSMeshLoader *mesh = new TACSMeshLoader(comm);
  mesh->incref();
  mesh->scanBDFFile(filename);

  // Get number of components prescribed in BDF file
  int num_components = mesh->getNumComponents();

  // Set properties needed to create stiffness object
  double rho = 2500.0; // density, kg/m^3
  double E = 70e9; // elastic modulus, Pa
  double nu = 0.3; // poisson's ratio
  double kcorr = 5.0/6.0; // shear correction factor
  double ys = 350e6; // yield stress, Pa

  // Loop over components, creating constituitive object for each
  for ( int i = 0; i < num_components; i++ ){
    const char *descriptor = mesh->getElementDescript(i);
    double min_thickness = 0.01;
    double max_thickness = 0.20;
    double thickness = 0.01;
    isoFSDTStiffness *stiff = new isoFSDTStiffness(rho, E, nu, kcorr, ys,
        thickness, i, min_thickness, max_thickness);

    // Initialize element object
    TACSElement *element = NULL;

    // Create element object using constituitive information and
    // type defined in descriptor
    if ( strcmp(descriptor, "CQUAD") == 0 ||
         strcmp(descriptor, "CQUADR") == 0 ||
         strcmp(descriptor, "CQUAD4") == 0) {
      element = new MITCShell<2>(stiff, LINEAR, i);
    }
    mesh->setElement(i, element);
  }

  // Create tacs assembler from mesh loader object
  TACSAssembler *tacs = mesh->createTACS(6);
  tacs->incref();
  mesh->decref();

  // Create the function of interest
  TACSFunction *func = new TACSKSFailure(tacs, 100.0);
  func->incref();

  // Get the deseign variable values
  TacsScalar *x = new TacsScalar[ num_components ];
  memset(x, 0, num_components*sizeof(TacsScalar));

  // Get the design variable values
  tacs->getDesignVars(x, num_components);

  // Create matrix and vectors
  TACSBVec *ans = tacs->createVec(); // displacements and rotations
  TACSBVec *f = tacs->createVec(); // loads
  TACSBVec *res = tacs->createVec(); // The residual
  TACSBVec *adjoint = tacs->createVec();
  FEMat *mat = tacs->createFEMat(); // stiffness matrix

  // Increment reference count to the matrix/vectors
  ans->incref();
  f->incref();
  res->incref();
  adjoint->incref();
  mat->incref();

  // Allocate the factorization
  int lev = 10000;
  double fill = 10.0;
  int reorder_schur = 1;
  PcScMat *pc = new PcScMat(mat, lev, fill, reorder_schur);
  pc->incref();

  // Set all the entries in load vector to specified value
  TacsScalar *force_vals;
  int size = f->getArray(&force_vals);
  for ( int k = 2; k < size; k += 6 ){
    force_vals[k] += 100.0;
  }
  tacs->applyBCs(f);

  // Assemble and factor the stiffness/Jacobian matrix. Factor the
  // Jacobian and solve the linear system for the displacements
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  tacs->assembleJacobian(alpha, beta, gamma, NULL, mat);
  pc->factor(); // LU factorization of stiffness matrix
  pc->applyFactor(f, ans);
  tacs->setVariables(ans);

  // Evaluate the function of interest
  TacsScalar fval;
  tacs->evalFunctions(&func, 1, &fval);

  TACSBVec *dfdXpts = tacs->createNodeVec();
  dfdXpts->incref();

  // Evaluate the total derivative
  TacsScalar *dfdx = new TacsScalar[ num_components ];

  // Solve the adjoint equations
  res->zeroEntries();
  tacs->addSVSens(1.0, 0.0, 0.0, &func, 1, &res);
  pc->applyFactor(res, adjoint);
  tacs->applyBCs(adjoint);

  // Compute the total derivative for the material variables
  memset(dfdx, 0, num_components*sizeof(TacsScalar));
  tacs->addDVSens(1.0, &func, 1, dfdx, num_components);
  tacs->addAdjointResProducts(-1.0, &adjoint, 1, dfdx, num_components);

  // Add the material derivative contributions across all processors
  MPI_Allreduce(MPI_IN_PLACE, dfdx, num_components, TACS_MPI_TYPE,
                MPI_SUM, comm);

  // Compute the total derivative for the adjoint variables
  dfdXpts->zeroEntries();
  tacs->addXptSens(1.0, &func, 1, &dfdXpts);
  tacs->addAdjointResXptSensProducts(-1.0, &adjoint, 1, &dfdXpts);

  // Finish adding everything up the geometry derivatives in parallel
  dfdXpts->beginSetValues(TACS_ADD_VALUES);
  dfdXpts->endSetValues(TACS_ADD_VALUES);

  // Compute the projected derivative
  TacsScalar dfdp = 0.0;
  for ( int i = 0; i < num_components; i++ ){
    dfdp += dfdx[i];
  }

  TacsScalar *xnew = new TacsScalar[ num_components ];
#ifdef TACS_USE_COMPLEX
  for ( int i = 0; i < num_components; i++ ){
    xnew[i] = x[i] + TacsScalar(0.0, dh);
  }
#else
  for ( int i = 0; i < num_components; i++ ){
    xnew[i] = x[i] + dh;
  }
#endif // TACS_USE_COMPLEX

  // Set the perturbed design variable values
  tacs->setDesignVars(xnew, num_components);

  // Re-solve the linear system
  tacs->assembleJacobian(alpha, beta, gamma, NULL, mat);
  pc->factor(); // LU factorization of stiffness matrix
  pc->applyFactor(f, ans);
  tacs->setVariables(ans);

  // Evaluate the function of interest
  TacsScalar fval2;
  tacs->evalFunctions(&func, 1, &fval2);

  if (rank == 0){
    printf("Adjoint:       %15.8e\n", TacsRealPart(dfdp));
#ifdef TACS_USE_COMPLEX
    double fd = TacsImagPart(fval2)/dh;
    printf("Complex step:  %15.8e\n", fd);
#else
    double fd = (fval2 - fval)/dh;
    printf("Finite-diff:   %15.8e\n", fd);
#endif // TACS_USE_COMPLEX
    printf("Rel error:     %15.8e\n", TacsRealPart((fd - dfdp)/dfdp));
  }

  // Reset the design variable values
  tacs->setDesignVars(x, num_components);

  // Create a new vector of points
  TACSBVec *X = tacs->createNodeVec();
  X->incref();
  tacs->getNodes(X);

  // Create a perturbation in the node locations
  TACSBVec *pert = tacs->createNodeVec();
  pert->incref();

  // Perturb the array
  TacsScalar *X_array, *pert_array;
  int X_size = X->getArray(&X_array);
  pert->getArray(&pert_array);
  for ( int i = 0; i < X_size/3; i++ ){
    pert_array[3*i] = X_array[3*i+1];
    pert_array[3*i+1] = X_array[3*i];
    pert_array[3*i+2] = X_array[3*i+2];
  }

  // Compute the projected derivative
  dfdp = dfdXpts->dot(pert);

#ifdef TACS_USE_COMPLEX
  X->axpy(TacsScalar(0.0, dh), pert);
#else
  X->axpy(dh, pert);
#endif

  // Set the new node locations
  tacs->setNodes(X);

  // Solve the equations again
  tacs->assembleJacobian(alpha, beta, gamma, NULL, mat);
  pc->factor(); // LU factorization of stiffness matrix
  pc->applyFactor(f, ans);
  tacs->setVariables(ans);

  // Re-evaluate the function
  tacs->evalFunctions(&func, 1, &fval2);

  if (rank == 0){
    printf("Adjoint:       %15.8e\n", TacsRealPart(dfdp));
#ifdef TACS_USE_COMPLEX
    double fd = TacsImagPart(fval2)/dh;
    printf("Complex step:  %15.8e\n", fd);
#else
    double fd = (fval2 - fval)/dh;
    printf("Finite-diff:   %15.8e\n", fd);
#endif // TACS_USE_COMPLEX
    printf("Rel error:     %15.8e\n", TacsRealPart((fd - dfdp)/dfdp));
  }

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 * f5 = new TACSToFH5(tacs, TACS_SHELL, write_flag);
  f5->incref();
  f5->writeToFile("ucrm.f5");
  f5->decref();

  // Free the KS failure
  func->decref();

  // Free the design variable information
  delete [] dfdx;
  delete [] x;
  delete [] xnew;

  // Decref the matrix/pc objects
  mat->decref();
  pc->decref();

  // Decref the vectors
  ans->decref();
  f->decref();
  res->decref();
  adjoint->decref();

  // Decref the nodal vectors
  dfdXpts->decref();
  pert->decref();
  X->decref();

  // Decref TACS
  tacs->decref();

  MPI_Finalize();

  return 0;
}
