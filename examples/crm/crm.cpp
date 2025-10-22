#include "TACSIsoShellConstitutive.h"
#include "TACSKSFailure.h"
#include "TACSMeshLoader.h"
#include "TACSShellElementDefs.h"

int main(int argc, char **argv) {
  // Intialize MPI and declare communicator
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int rank;
  MPI_Comm_rank(comm, &rank);

#ifdef TACS_USE_COMPLEX
  double dh = 1e-30;
#else
  double dh = 1e-6;
#endif  // TACS_USE_COMPLEX

  // Scan the arguments to check for a manual step size
  for (int k = 0; k < argc; k++) {
    if (sscanf(argv[k], "dh=%lf", &dh) == 1) {
      if (rank == 0) {
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
  TacsScalar rho = 2500.0;  // density, kg/m^3
  TacsScalar specific_heat = 921.096;
  TacsScalar E = 70e9;       // elastic modulus, Pa
  TacsScalar nu = 0.3;       // poisson's ratio
  TacsScalar ys = 350e6;     // yield stress, Pa
  TacsScalar cte = 24.0e-6;  // Coefficient of thermal expansion
  TacsScalar kappa = 230.0;  // Thermal conductivity

  TACSMaterialProperties *props =
      new TACSMaterialProperties(rho, specific_heat, E, nu, ys, cte, kappa);

  TACSShellTransform *transform = new TACSShellNaturalTransform();

  // Loop over components, creating constituitive object for each
  for (int i = 0; i < num_components; i++) {
    const char *descriptor = mesh->getElementDescript(i);
    TacsScalar thickness = 0.01;
    int thickness_index = i;
    TacsScalar min_thickness = 0.01;
    TacsScalar max_thickness = 0.20;
    TACSShellConstitutive *con = new TACSIsoShellConstitutive(
        props, thickness, thickness_index, min_thickness, max_thickness);

    // Initialize element object
    TACSElement *shell = TacsCreateShellByName(descriptor, transform, con);

    // Set the shell element
    mesh->setElement(i, shell);
  }

  // Create tacs assembler from mesh loader object
  TACSAssembler *assembler = mesh->createTACS(6);
  assembler->incref();
  mesh->decref();

  // Create the function of interest
  TACSFunction *func = new TACSKSFailure(assembler, 100.0);
  func->incref();

  // Create the design vector
  TACSBVec *x = assembler->createDesignVec();
  x->incref();

  // Get the design variable values
  assembler->getDesignVars(x);

  // Create matrix and vectors
  TACSBVec *ans = assembler->createVec();  // displacements and rotations
  TACSBVec *f = assembler->createVec();    // loads
  TACSBVec *res = assembler->createVec();  // The residual
  TACSBVec *adjoint = assembler->createVec();
  TACSSchurMat *mat = assembler->createSchurMat();  // stiffness matrix

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
  TACSSchurPc *pc = new TACSSchurPc(mat, lev, fill, reorder_schur);
  pc->incref();

  // Set all the entries in load vector to specified value
  TacsScalar *force_vals;
  int size = f->getArray(&force_vals);
  for (int k = 2; k < size; k += 6) {
    force_vals[k] += 100.0;
  }
  assembler->applyBCs(f);

  // Assemble and factor the stiffness/Jacobian matrix. Factor the
  // Jacobian and solve the linear system for the displacements
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  assembler->assembleJacobian(alpha, beta, gamma, NULL, mat);
  pc->factor();  // LU factorization of stiffness matrix
  pc->applyFactor(f, ans);
  assembler->setVariables(ans);

  // Evaluate the function of interest
  TacsScalar fval;
  assembler->evalFunctions(1, &func, &fval);

  TACSBVec *dfdx = assembler->createDesignVec();
  dfdx->incref();

  TACSBVec *dfdXpts = assembler->createNodeVec();
  dfdXpts->incref();

  // Solve the adjoint equations
  res->zeroEntries();
  assembler->addSVSens(1.0, 0.0, 0.0, 1, &func, &res);
  pc->applyFactor(res, adjoint);
  assembler->applyBCs(adjoint);

  // Compute the total derivative for the material variables
  dfdx->zeroEntries();
  assembler->addDVSens(1.0, 1, &func, &dfdx);
  assembler->addAdjointResProducts(-1.0, 1, &adjoint, &dfdx);

  dfdx->beginSetValues(TACS_ADD_VALUES);
  dfdx->endSetValues(TACS_ADD_VALUES);

  // Compute the total derivative for the adjoint variables
  dfdXpts->zeroEntries();
  assembler->addXptSens(1.0, 1, &func, &dfdXpts);
  assembler->addAdjointResXptSensProducts(-1.0, 1, &adjoint, &dfdXpts);

  // Finish adding everything up the geometry derivatives in parallel
  dfdXpts->beginSetValues(TACS_ADD_VALUES);
  dfdXpts->endSetValues(TACS_ADD_VALUES);

  // Create a random vector and set it to a random array
  TACSBVec *p = assembler->createDesignVec();
  p->incref();
  p->setRand(-1.0, 1.0);

  // Compute the projected derivative
  TacsScalar dfdp = dfdx->dot(p);

  // Create a new vector of design variables
  TACSBVec *xnew = assembler->createDesignVec();
  xnew->incref();

#ifdef TACS_USE_COMPLEX
  xnew->copyValues(x);
  xnew->axpy(TacsScalar(0.0, dh), p);
#else
  xnew->copyValues(x);
  xnew->axpy(dh, p);
#endif  // TACS_USE_COMPLEX

  // Set the perturbed design variable values
  assembler->setDesignVars(xnew);

  // Re-solve the linear system
  assembler->assembleJacobian(alpha, beta, gamma, NULL, mat);
  pc->factor();  // LU factorization of stiffness matrix
  pc->applyFactor(f, ans);
  assembler->setVariables(ans);

  // Evaluate the function of interest
  TacsScalar fval2;
  assembler->evalFunctions(1, &func, &fval2);

  if (rank == 0) {
    printf("Adjoint:       %15.8e\n", TacsRealPart(dfdp));
#ifdef TACS_USE_COMPLEX
    double fd = TacsImagPart(fval2) / dh;
    printf("Complex step:  %15.8e\n", fd);
#else
    double fd = (fval2 - fval) / dh;
    printf("Finite-diff:   %15.8e\n", fd);
#endif  // TACS_USE_COMPLEX
    printf("Rel error:     %15.8e\n", TacsRealPart((fd - dfdp) / dfdp));
  }

  // Reset the design variable values
  assembler->setDesignVars(x);

  // Create a new vector of points
  TACSBVec *X = assembler->createNodeVec();
  X->incref();
  assembler->getNodes(X);

  // Create a perturbation in the node locations
  TACSBVec *pert = assembler->createNodeVec();
  pert->incref();
  pert->setRand(-1.0, 1.0);

  // Compute the projected derivative
  dfdp = dfdXpts->dot(pert);

#ifdef TACS_USE_COMPLEX
  X->axpy(TacsScalar(0.0, dh), pert);
#else
  X->axpy(dh, pert);
#endif

  // Set the new node locations
  assembler->setNodes(X);

  // Solve the equations again
  assembler->assembleJacobian(alpha, beta, gamma, NULL, mat);
  pc->factor();  // LU factorization of stiffness matrix
  pc->applyFactor(f, ans);
  assembler->setVariables(ans);

  // Re-evaluate the function
  assembler->evalFunctions(1, &func, &fval2);

  if (rank == 0) {
    printf("Adjoint:       %15.8e\n", TacsRealPart(dfdp));
#ifdef TACS_USE_COMPLEX
    double fd = TacsImagPart(fval2) / dh;
    printf("Complex step:  %15.8e\n", fd);
#else
    double fd = (fval2 - fval) / dh;
    printf("Finite-diff:   %15.8e\n", fd);
#endif  // TACS_USE_COMPLEX
    printf("Rel error:     %15.8e\n", TacsRealPart((fd - dfdp) / dfdp));
  }

  // Create an TACSToFH5 object for writing output to files
  int write_flag = (TACS_OUTPUT_CONNECTIVITY | TACS_OUTPUT_NODES |
                    TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                    TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
  TACSToFH5 *f5 =
      new TACSToFH5(assembler, TACS_BEAM_OR_SHELL_ELEMENT, write_flag);
  f5->incref();
  f5->writeToFile("ucrm.f5");
  f5->decref();

  // Free the KS failure
  func->decref();

  // Free the design variable information
  dfdx->decref();
  x->decref();
  xnew->decref();
  p->decref();

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
  assembler->decref();

  MPI_Finalize();

  return 0;
}
