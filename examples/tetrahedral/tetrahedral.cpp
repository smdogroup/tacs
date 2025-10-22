#include "TACSCreator.h"
#include "TACSElement3D.h"
#include "TACSElementVerification.h"
#include "TACSHeatConduction.h"
#include "TACSKSFailure.h"
#include "TACSKSTemperature.h"
#include "TACSLinearElasticity.h"
#include "TACSMeshLoader.h"
#include "TACSStructuralMass.h"
#include "TACSTetrahedralBasis.h"
#include "TACSThermoelasticity.h"
#include "TACSToFH5.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // Check whether to use elasticity or thoermoelasticity
  int analysis_type = 0;
  int function_type = 0;
  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "conduction") == 0) {
      analysis_type = 1;
    } else if (strcmp(argv[i], "thermoelasticity") == 0) {
      analysis_type = 2;
    } else if (strcmp(argv[i], "failure") == 0) {
      function_type = 1;
    } else if (strcmp(argv[i], "temperature") == 0) {
      function_type = 2;
    }
  }

  // Create the mesh loader object on MPI_COMM_WORLD. The
  // TACSAssembler object will be created on the same comm
  MPI_Comm comm = MPI_COMM_WORLD;
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  TACSMeshLoader *mesh = new TACSMeshLoader(comm);
  mesh->incref();

  // Create the isotropic material class
  TacsScalar rho = 2700.0;
  TacsScalar specific_heat = 921.096;
  TacsScalar E = 70e3;
  TacsScalar nu = 0.3;
  TacsScalar ys = 270.0;
  TacsScalar cte = 24.0e-6;
  TacsScalar kappa = 230.0;
  TACSMaterialProperties *props =
      new TACSMaterialProperties(rho, specific_heat, E, nu, ys, cte, kappa);

  // Create stiffness (need class)
  TACSSolidConstitutive *stiff = new TACSSolidConstitutive(props, 1.0, 0);
  stiff->incref();

  // Create model (need class)
  TACSElementModel *model = NULL;
  if (analysis_type == 1) {
    model = new TACSHeatConduction3D(stiff);
  } else if (analysis_type == 2) {
    model = new TACSLinearThermoelasticity3D(stiff, TACS_LINEAR_STRAIN);
  } else {
    model = new TACSLinearElasticity3D(stiff, TACS_LINEAR_STRAIN);
  }
  int vars_per_node = model->getVarsPerNode();

  // Create basis
  TACSElementBasis *linear_basis = new TACSLinearTetrahedralBasis();
  TACSElementBasis *quad_basis = new TACSQuadraticTetrahedralBasis();

  if (mpi_rank == 0) {
    TacsTestElementBasis(linear_basis);
    TacsTestElementBasis(quad_basis);
  }

  // Create the element type (need 3D element class)
  TACSElement3D *linear_element = new TACSElement3D(model, linear_basis);
  TACSElement3D *quad_element = new TACSElement3D(model, quad_basis);

  // The TACSAssembler object - which should be allocated if the mesh
  // is loaded correctly
  TACSAssembler *assembler = NULL;

  // Try to load the input file as a BDF file through the
  // TACSMeshLoader class
  if (argc > 1) {
    const char *filename = argv[1];
    FILE *fp = fopen(filename, "r");
    if (fp) {
      fclose(fp);

      // Scan the BDF file
      int fail = mesh->scanBDFFile(filename);

      if (fail) {
        fprintf(stderr, "Failed to read in the BDF file\n");
      } else {
        // Add the elements to the mesh loader class
        for (int i = 0; i < mesh->getNumComponents(); i++) {
          TACSElement *elem = NULL;

          // Get the BDF description of the element
          const char *elem_descript = mesh->getElementDescript(i);
          if (strcmp(elem_descript, "CTETRA") == 0 ||
              strcmp(elem_descript, "CTETRA4") == 0) {
            elem = linear_element;
          } else if (strcmp(elem_descript, "CTETRA10") == 0) {
            elem = quad_element;
          }

          // Set the element object into the mesh loader class
          if (elem) {
            mesh->setElement(i, elem);
          }
        }

        // Now, create the TACSAssembler object
        assembler = mesh->createTACS(vars_per_node);
        assembler->incref();
      }
    } else {
      fprintf(stderr, "File %s does not exist\n", filename);
    }
  } else {
    fprintf(stderr, "No BDF file provided\n");
  }

  if (assembler) {
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
    int nrestart = 2;     // Number of allowed restarts
    int is_flexible = 0;  // Is a flexible preconditioner?
    TACSKsm *ksm = new GMRES(mat, pc, gmres_iters, nrestart, is_flexible);
    ksm->incref();

    // Assemble and factor the stiffness/Jacobian matrix
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    assembler->assembleJacobian(alpha, beta, gamma, res, mat);
    pc->factor();

    res->set(1.0);
    assembler->applyBCs(res);
    ksm->solve(res, ans);
    assembler->setVariables(ans);

#ifdef TACS_USE_COMPLEX
    assembler->testElement(0, 2, 1e-30);
#else
    assembler->testElement(0, 2);
#endif

    // The function that we will use: The KS failure function evaluated
    // over all the elements in the mesh
    TACSFunction *func = NULL;
    if (function_type == 1) {
      func = new TACSKSFailure(assembler, 100.0);
    } else if (function_type == 2) {
      func = new TACSKSTemperature(assembler, 100.0);
    } else {
      func = new TACSStructuralMass(assembler);
    }
    func->incref();

    // Evaluate the function
    TacsScalar mass = 0.0;
    assembler->evalFunctions(1, &func, &mass);
    printf("%s: %e\n", func->getObjectName(), TacsRealPart(mass));

#ifdef TACS_USE_COMPLEX
    assembler->testFunction(func, 1e-30);
#else
    assembler->testFunction(func, 1e-6);
#endif  // TACS_USE_COMPLEX

    // Create an TACSToFH5 object for writing output to files
    ElementType etype = TACS_SOLID_ELEMENT;
    int write_flag = (TACS_OUTPUT_CONNECTIVITY | TACS_OUTPUT_NODES |
                      TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                      TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
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
  if (assembler) {
    assembler->decref();
  }

  MPI_Finalize();
  return (0);
}
