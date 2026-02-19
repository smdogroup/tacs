#include "JacobiDavidson.h"
#include "TACSBuckling.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSMeshLoader.h"
#include "TACSShellElementDefs.h"

int main(int argc, char **argv) {
  // Intialize MPI and declare communicator
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int use_lanczos = 1;
  int use_tacs_freq = 0;
  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "jd") == 0) {
      use_lanczos = 0;
    }
    if (strcmp(argv[i], "jd_freq") == 0) {
      use_lanczos = 0;
      use_tacs_freq = 1;
    }
  }

  int rank;
  MPI_Comm_rank(comm, &rank);

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

  // Create matrix and vectors
  TACSSchurMat *kmat = assembler->createSchurMat();  // stiffness matrix
  kmat->incref();

  TACSSchurMat *mmat = assembler->createSchurMat();
  mmat->incref();

  // Parameters for the preconditioner
  int lev_fill = 10000;
  double fill = 10.0;
  int reorder_schur = 1;

  if (use_lanczos) {
    TACSSchurPc *pc = new TACSSchurPc(kmat, lev_fill, fill, reorder_schur);
    pc->incref();

    TACSBVec *vec = assembler->createVec();
    vec->incref();

    // Assemble and factor the stiffness/Jacobian matrix. Factor the
    // Jacobian and solve the linear system for the displacements
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    assembler->assembleJacobian(alpha, beta, gamma, NULL, kmat);
    pc->factor();  // LU factorization of stiffness matrix

    int gmres_iters = 15;
    int nrestart = 0;     // Number of allowed restarts
    int is_flexible = 0;  // Is a flexible preconditioner?
    GMRES *ksm = new GMRES(kmat, pc, gmres_iters, nrestart, is_flexible);
    ksm->incref();
    ksm->setTolerances(1e-12, 1e-30);

    // Perform the Frequency Analysis
    int max_lanczos = 60;
    int num_eigvals = 20;
    double eig_tol = 1e-8;
    TacsScalar sigma = 10.0;
    int output_freq = 1;

    KSMPrint *ksm_print = new KSMPrintStdout("KSM", rank, output_freq);
    TACSFrequencyAnalysis *freq_analysis = new TACSFrequencyAnalysis(
        assembler, sigma, mmat, kmat, ksm, max_lanczos, num_eigvals, eig_tol);
    freq_analysis->incref();

    double t1 = MPI_Wtime();
    freq_analysis->solve(ksm_print);
    t1 = MPI_Wtime() - t1;
    printf("Lanczos time: %15.5e\n", t1);

    TacsScalar freq;
    for (int k = 0; k < num_eigvals; k++) {
      TacsScalar error;
      TacsScalar eigvalue = freq_analysis->extractEigenvector(k, vec, &error);
      eigvalue = sqrt(eigvalue);
      freq = TacsRealPart(eigvalue) / (2.0 * 3.14159);

      printf("TACS frequency[%2d]: %15.6f %15.6f\n", k, TacsRealPart(eigvalue),
             TacsRealPart(freq));
    }

    pc->decref();
    vec->decref();
    ksm->decref();
    freq_analysis->decref();
  } else {
    TACSSchurMat *pcmat = assembler->createSchurMat();
    TACSSchurPc *pc = new TACSSchurPc(pcmat, lev_fill, fill, reorder_schur);
    pc->incref();

    int num_eigvals = 20;
    int jd_size = 25;
    int fgmres_size = 5;
    double sigma = 0.0;
    if (use_tacs_freq) {
      double rtol = 1e-3;
      double atol = 1e-16;
      TACSFrequencyAnalysis *freq_analysis = new TACSFrequencyAnalysis(
          assembler, sigma, mmat, kmat, pcmat, pc, jd_size, fgmres_size,
          num_eigvals, rtol, atol);
      freq_analysis->incref();

      double t1 = MPI_Wtime();
      int output_freq = 1;
      KSMPrint *ksm_print = new KSMPrintStdout("KSM", rank, output_freq);
      freq_analysis->solve(ksm_print);
      t1 = MPI_Wtime() - t1;

      printf("Jacobi-Davidson time: %15.5e\n", t1);
      for (int k = 0; k < num_eigvals; k++) {
        TacsScalar eigvalue = freq_analysis->extractEigenvalue(k, NULL);
        eigvalue = sqrt(eigvalue);
        TacsScalar freq = TacsRealPart(eigvalue) / (2.0 * 3.14159);

        printf("TACS frequency[%2d]: %15.6f %15.6f\n", k,
               TacsRealPart(eigvalue), TacsRealPart(freq));
      }
      freq_analysis->decref();
    } else {
      assembler->assembleMatType(TACS_MASS_MATRIX, mmat);
      assembler->assembleMatType(TACS_STIFFNESS_MATRIX, kmat);

      TACSJDFrequencyOperator *oper =
          new TACSJDFrequencyOperator(assembler, kmat, mmat, pcmat, pc);
      oper->setEigenvalueEstimate(0.0);

      TACSJacobiDavidson *jd =
          new TACSJacobiDavidson(oper, num_eigvals, jd_size, fgmres_size);
      jd->incref();

      int output_freq = 1;
      KSMPrint *ksm_print = new KSMPrintStdout("KSM", rank, output_freq);
      double t1 = MPI_Wtime();
      jd->solve(ksm_print);
      t1 = MPI_Wtime() - t1;
      printf("Jacobi-Davidson time: %15.5e\n", t1);

      for (int k = 0; k < num_eigvals; k++) {
        TacsScalar eigvalue = jd->extractEigenvalue(k, NULL);
        eigvalue = sqrt(eigvalue);
        TacsScalar freq = TacsRealPart(eigvalue) / (2.0 * 3.14159);

        printf("TACS frequency[%2d]: %15.6f %15.6f\n", k,
               TacsRealPart(eigvalue), TacsRealPart(freq));
      }

      jd->decref();
    }
    pc->decref();
  }

  // Decrefs
  kmat->decref();
  mmat->decref();
  assembler->decref();

  MPI_Finalize();

  return 0;
}
