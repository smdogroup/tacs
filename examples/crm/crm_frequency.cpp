#include "TACSMeshLoader.h"
#include "MITCShell.h"
#include "isoFSDTStiffness.h"
#include "TACSBuckling.h"
#include "JacobiDavidson.h"

int main( int argc, char **argv ){
  // Intialize MPI and declare communicator
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int use_lanczos = 1;
  int use_tacs_freq = 0;
  for ( int i = 0; i < argc; i++ ){
    if (strcmp(argv[i], "jd") == 0){
      use_lanczos = 0;
    }
    if (strcmp(argv[i], "jd_freq") == 0){
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

  // Get the design variable values
  TacsScalar *x = new TacsScalar[ num_components ];
  memset(x, 0, num_components*sizeof(TacsScalar));

  // Get the design variable values
  tacs->getDesignVars(x, num_components);

  // Set the design variables to a uniform 5mm
  for (int i = 0; i < num_components; i++){
    x[i] = 0.015;
  }
  tacs->setDesignVars(x, num_components);

  int lev_fill = 5000; // ILU(k) fill in
  int fill = 8.0;      // Expected number of non-zero entries

  // Create matrix and vectors
  FEMat *kmat = tacs->createFEMat(); // stiffness matrix
  kmat->incref();

  FEMat *mmat = tacs->createFEMat();
  mmat->incref();

  if (use_lanczos){
    PcScMat *pc = new PcScMat(kmat, lev_fill, fill, 1);
    pc->incref();

    TACSBVec *vec = tacs->createVec();
    vec->incref();

    // Assemble and factor the stiffness/Jacobian matrix. Factor the
    // Jacobian and solve the linear system for the displacements
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    tacs->assembleJacobian(alpha, beta, gamma, NULL, kmat);
    pc->factor(); // LU factorization of stiffness matrix

    int gmres_iters = 15;
    int nrestart = 0; // Number of allowed restarts
    int is_flexible = 0; // Is a flexible preconditioner?
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
            tacs, sigma, mmat, kmat, ksm, max_lanczos, num_eigvals, eig_tol);
    freq_analysis->incref();

    double t1 = MPI_Wtime();
    freq_analysis->solve(ksm_print);
    t1 = MPI_Wtime() - t1;
    printf("Lanczos time: %15.5e\n", t1);

    TacsScalar freq;
    for ( int k = 0; k < num_eigvals; k++ ){
      TacsScalar error;
      TacsScalar eigvalue = freq_analysis->extractEigenvector(k, vec, &error);
      eigvalue = sqrt(eigvalue);
      freq = TacsRealPart(eigvalue)/(2.0*3.14159);

      printf("TACS frequency[%2d]: %15.6f %15.6f\n", 
             k, TacsRealPart(eigvalue), TacsRealPart(freq));
    }

    pc->decref();
    vec->decref();
    ksm->decref();
    freq_analysis->decref();
  }
  else {
    FEMat *pcmat = tacs->createFEMat();
    PcScMat *pc = new PcScMat(pcmat, lev_fill, fill, 1);
    pc->incref();
    
    int num_eigvals = 20;
    int jd_size = 25;
    int fgmres_size = 5;
    double sigma = 0.0;
    if (use_tacs_freq){
      double rtol = 1e-3;
      double atol = 1e-16;
      TACSFrequencyAnalysis *freq_analysis = new TACSFrequencyAnalysis(tacs, sigma, mmat,
                                                                       kmat, pcmat, pc,
                                                                       jd_size,
                                                                       fgmres_size, 
                                                                       num_eigvals, rtol,
                                                                       atol);
      freq_analysis->incref();
      double t1 = MPI_Wtime();
      int output_freq = 1;
      KSMPrint *ksm_print = new KSMPrintStdout("KSM", rank, output_freq);
      freq_analysis->solve(ksm_print);
      t1 = MPI_Wtime() - t1;
      printf("Jacobi-Davidson time: %15.5e\n", t1);
      for ( int k = 0; k < num_eigvals; k++ ){
        TacsScalar eigvalue = freq_analysis->extractEigenvalue(k, NULL);
        eigvalue = sqrt(eigvalue);
        TacsScalar freq = TacsRealPart(eigvalue)/(2.0*3.14159);

        printf("TACS frequency[%2d]: %15.6f %15.6f\n", 
               k, TacsRealPart(eigvalue), TacsRealPart(freq));
      }
      freq_analysis->decref();
    }
    else{
      tacs->assembleMatType(MASS_MATRIX, mmat);
      tacs->assembleMatType(STIFFNESS_MATRIX, kmat);

      TACSJDFrequencyOperator *oper = 
        new TACSJDFrequencyOperator(tacs, kmat, mmat, pcmat, pc);
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

      for ( int k = 0; k < num_eigvals; k++ ){
        TacsScalar eigvalue = jd->extractEigenvalue(k, NULL);
        eigvalue = sqrt(eigvalue);
        TacsScalar freq = TacsRealPart(eigvalue)/(2.0*3.14159);

        printf("TACS frequency[%2d]: %15.6f %15.6f\n", 
               k, TacsRealPart(eigvalue), TacsRealPart(freq));
      }

      jd->decref();
    }
    pc->decref();
  }

  // Decrefs
  kmat->decref();
  mmat->decref();
  tacs->decref();

  MPI_Finalize();

  return 0;
}
