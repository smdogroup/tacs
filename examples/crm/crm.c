#include "TACSMeshLoader.h"
#include "MITCShell.h"
#include "isoFSDTStiffness.h"

int main( int argc, char **argv ){
  // Intialize MPI and declare communicator
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

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
    double thickness = 0.07;
    isoFSDTStiffness *stiff = new isoFSDTStiffness(rho, E, nu, kcorr, ys,
        thickness, i, min_thickness, max_thickness); 

    // Initialize element object
    TACSElement *element = NULL;

    // Create element object using constituitive information and type defined in
    // descriptor
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

#ifdef TACS_USE_COMPLEX
  // Get the deseign variable values
  TacsScalar *x = new TacsScalar[ num_components ];
  memset(x, 0, num_components*sizeof(TacsScalar));

  // Get the design variable values
  tacs->getDesignVars(x, num_components);

  double dh = 1e-30;
  for ( int i = 0; i < num_components; i++ ){
    x[i] = x[i] + TacsScalar(0.0, dh);
  }

  // Set the perturbed design variable values
  tacs->setDesignVars(x, num_components);
#endif // TACS_USE_COMPLEX

  // Create matrix and vectors 
  TACSBVec *ans = tacs->createVec(); // displacements and rotations
  TACSBVec *f = tacs->createVec(); // loads
  FEMat *mat = tacs->createFEMat(); // preconditioner

  // Increment reference count to the matrix/vectors
  ans->incref();
  f->incref();
  mat->incref();

  // Allocate the factorization
  int lev = 10000;
  double fill = 10.0;
  int reorder_schur = 1;
  PcScMat *pc = new PcScMat(mat, lev, fill, reorder_schur); 
  pc->incref();

  // Assemble and factor the stiffness/Jacobian matrix
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  tacs->assembleJacobian(alpha, beta, gamma, NULL, mat);
  mat->applyBCs();
  pc->factor(); // LU factorization of stiffness matrix

  // Set all the entries in load vector to specified value
  f->set(1.0);
  f->applyBCs();

  // Get solution and store in ans
  pc->applyFactor(f, ans);
  tacs->setVariables(ans);

  // Print out the norm of the solution on the root processor
  TacsScalar norm = ans->norm();
  int rank;
  MPI_Comm_rank(comm, &rank);

  if (rank == 0){
#ifdef TACS_USE_COMPLEX
    printf("||ans||: %e + j%e\n", RealPart(norm), ImagPart(norm));
#else
    printf("||ans||: %e\n", ans->norm());
#endif
  }

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 * f5 = new TACSToFH5(tacs, SHELL, write_flag);
  f5->incref();

  // Write the displacements
  f5->writeToFile("ucrm.f5");

  // Delete the viewer
  f5->decref();

  mat->decref();
  pc->decref();
  tacs->decref();

  MPI_Finalize();

  return 0;
}
