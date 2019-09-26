#include "TACSLinearElasticity.h"
#include "TACSQuadBasis.h"
#include "TACSElement2D.h"
#include "TACSToFH5.h"
#include "TACSMeshLoader.h"

/*
  The following test illustrates the use of PlaneStressQuad elements
  and demonstrates how to load a BDF file using the TACSMeshLoader
  object. This example can be run in parallel.
*/
int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);

  // Create the mesh loader object on MPI_COMM_WORLD. The
  // TACSAssembler object will be created on the same comm
  TACSMeshLoader *mesh = new TACSMeshLoader(MPI_COMM_WORLD);
  mesh->incref();

  // Create the isotropic material class
  TacsScalar rho = 2700.0;
  TacsScalar E = 70e3;
  TacsScalar nu = 0.3;
  TacsScalar ys = 270.0;
  TacsScalar cte = 0.0, kappa = 0.0;
  TACSMaterialProperties *props =
    new TACSMaterialProperties(rho, E, nu, ys, cte, kappa);

  // Create the stiffness object
  TACSPlaneStressConstitutive *stiff =
    new TACSPlaneStressConstitutive(props);
  stiff->incref();

  // Create the model class
  TACSLinearElasticity2D model(stiff, TACS_LINEAR_STRAIN);

  // Create the basis
  TACSElementBasis *linear_basis = new TACSLinearQuadBasis();
  TACSElementBasis *quad_basis = new TACSQuadraticQuadBasis();
  TACSElementBasis *cubic_basis = new TACSCubicQuadBasis();

  // Create the element type
  TACSElement2D *linear_element = new TACSElement2D(&model, linear_basis);
  TACSElement2D *quad_element = new TACSElement2D(&model, quad_basis);
  TACSElement2D *cubic_element = new TACSElement2D(&model, cubic_basis);

  // The TACSAssembler object - which should be allocated if the mesh
  // is loaded correctly
  TACSAssembler *assembler = NULL;

  // Try to load the input file as a BDF file through the
  // TACSMeshLoader class
  if (argc > 1){
    const char *filename = argv[1];
    FILE *fp = fopen(filename, "r");
    if (fp){
      fclose(fp);

      // Scan the BDF file
      int fail = mesh->scanBDFFile(filename);

      if (fail){
        fprintf(stderr, "Failed to read in the BDF file\n");
      }
      else {
        // Add the elements to the mesh loader class
        for ( int i = 0; i < mesh->getNumComponents(); i++ ){
          TACSElement *elem = NULL;

          // Get the BDF description of the element
          const char *elem_descript = mesh->getElementDescript(i);
          if (strcmp(elem_descript, "CQUAD4") == 0){
            elem = linear_element;
          }
          else if (strcmp(elem_descript, "CQUAD") == 0 ||
                   strcmp(elem_descript, "CQUAD9") == 0){
            elem = quad_element;
          }
          else if (strcmp(elem_descript, "CQUAD16") == 0){
            elem = cubic_element;
          }

          // Set the element object into the mesh loader class
          if (elem){
            mesh->setElement(i, elem);
          }
        }

        // Now, create the TACSAssembler object
        int vars_per_node = 2;
        assembler = mesh->createTACS(vars_per_node);
        assembler->incref();
      }
    }
    else {
      fprintf(stderr, "File %s does not exist\n", filename);
    }
  }
  else {
    fprintf(stderr, "No BDF file provided\n");
  }

  if (assembler){
    // Create the preconditioner
    TACSBVec *res = assembler->createVec();
    TACSBVec *ans = assembler->createVec();
    TACSSchurMat *mat = assembler->createSchurMat();

    // Increment the reference count to the matrix/vectors
    res->incref();
    ans->incref();
    mat->incref();

    // Allocate the factorization
    int lev = 4500;
    double fill = 10.0;
    int reorder_schur = 1;
    TACSSchurPc *pc = new TACSSchurPc(mat, lev, fill, reorder_schur);
    pc->incref();

    // Allocate the GMRES object
    int gmres_iters = 80;
    int nrestart = 2; // Number of allowed restarts
    int is_flexible = 0; // Is a flexible preconditioner?
    TACSKsm *ksm = new GMRES(mat, pc, gmres_iters,
                             nrestart, is_flexible);
    ksm->incref();

    // Assemble and factor the stiffness/Jacobian matrix
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    assembler->assembleJacobian(alpha, beta, gamma, res, mat);
    pc->factor();

    res->set(1.0);
    assembler->applyBCs(res);
    ksm->solve(res, ans);
    assembler->setVariables(ans);

    // Create an TACSToFH5 object for writing output to files
    ElementType etype = TACS_PLANE_STRESS_ELEMENT;
    int write_flag = (TACS_OUTPUT_CONNECTIVITY |
                      TACS_OUTPUT_NODES |
                      TACS_OUTPUT_DISPLACEMENTS |
                      TACS_OUTPUT_STRAINS |
                      TACS_OUTPUT_STRESSES |
                      TACS_OUTPUT_EXTRAS);
    TACSToFH5 *f5 = new TACSToFH5(assembler, etype, write_flag);
    f5->incref();
    f5->writeToFile("output.f5");

    // Free everything
    f5->decref();

    // Decrease the reference count to the linear algebra objects
    ksm->decref();
    pc->decref();
    mat->decref();
    ans->decref();
    res->decref();
  }

  // Deallocate the objects
  mesh->decref();
  stiff->decref();
  if (assembler){ assembler->decref(); }

  MPI_Finalize();
  return (0);
}
