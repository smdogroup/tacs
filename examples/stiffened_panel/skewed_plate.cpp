#include "TACSAssembler.h"
#include "TACSBuckling.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSMeshLoader.h"
#include "TACSShellElementDefs.h"
#include "TACSToFH5.h"
// #include "TACSPanelAnalysis.h"

// /*
//   The TACSPanelAnalysis object cannot be used in complex mode because
//   the underlying implementation uses eigenvalue routines from LAPACK
//   that are not "complexified".
// */

// #ifndef TACS_USE_COMPLEX

// /*
//   Build a model of a skewed plate using a TACSPanelAnalysis object
//   and perform buckling and natural frequency analysis on the model.

//   input:
//   skin:    stiffness object
//   theta:   skew angle
//   a:       the longitudinal plate dimension
//   b:       the transverse plate dimension
//   nloads:  the number of buckling loads
//   nfreq:   the number of natural frequencies of vibration

//   output:
//   loads:   the buckling loads
//   freq:    the natural frequencies of vibration
// */
// void skewed_test( FSDTStiffness *skin, TacsScalar theta,
//                   TacsScalar a, TacsScalar b,
//                   TacsScalar Nx, TacsScalar Nxy,
//                   TacsScalar loads[], int nloads,
//                   TacsScalar freq[], int nfreq ){
//   int nseg = 12;
//   int nnodes = nseg+1;
//   int nbeams = 0;
//   int nmodes = 20;
//   TACSPanelAnalysis *panel = new TACSPanelAnalysis(nnodes, nseg,
//                                                     nbeams, nmodes, a,
//                                                     theta);
//   panel->incref();

//   TacsScalar *Xpts = new TacsScalar[2*nnodes];
//   memset(Xpts, 0, 2*nnodes*sizeof(TacsScalar));

//   for ( int k = 0; k < nnodes; k++ ){
//     Xpts[2*k] = (k*b)/(nnodes-1);
//   }

//   for ( int k = 0; k < nnodes-1; k++ ){
//     panel->setSegment(k, TACSPanelAnalysis::SKIN_SEGMENT,
//                       skin, k, k+1);
//   }

//   panel->setPoints(Xpts, nnodes);

//   int cons[4] = {-1, -1, -1, 1};
//   panel->setFirstNodeBC(0, (1 | 2 | 4));
//   panel->setLastNodeBC(nnodes-1, (1 | 2 | 4));
//   panel->initialize();

//   double tp = MPI_Wtime();
//   panel->computePressureLoad(1.0, "./results/pressure_load.dat");
//   tp = MPI_Wtime() - tp;

//   double tb = MPI_Wtime();
//   panel->computeBucklingLoads(Nx, Nxy,
//                               loads, nloads, "./results/skewed_");
//   tb = MPI_Wtime() - tb;

//   double tf = MPI_Wtime();
//   panel->computeFrequencies(freq, nfreq, "./results/skewed_");
//   tf = MPI_Wtime() - tf;

//   printf("Pressure analysis time:  %f\n", tp);
//   printf("Buckling analysis time:  %f\n", tb);
//   printf("Frequency analysis time: %f\n", tf);

//   delete [] Xpts;
//   panel->decref();
// }

#endif  // TACS_USE_COMPLEX

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // Flag to indicate
  int isotropic_flag = 0;

  // Skew angle for the plate
  double theta = 0.0;

  // Loads applied to the plate
  double Nx = -1.0, Nxy = 0.0;

  // Parse the input arguments
  for (int i = 0; i < argc; i++) {
    if (sscanf(argv[i], "theta=%lf", &theta) == 1) {
      // Convert from degrees to radians
      theta *= (M_PI / 180.0);
      continue;
    }
    if (sscanf(argv[i], "Nx=%lf", &Nx) == 1) {
      continue;
    }
    if (sscanf(argv[i], "Nxy=%lf", &Nxy) == 1) {
      continue;
    }
    if (strcmp(argv[i], "composite") == 0) {
      isotropic_flag = 1;
    }
    if (strcmp(argv[i], "isotropic") == 0) {
      isotropic_flag = 1;
    }
  }

  // Get the processor rank
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // The stiffness object
  FSDTStiffness *stiff = NULL;

  // Equivalent values of the density, thickness, modulus and Poisson
  // ratio
  TacsScalar rho = 0.0, t = 0.0;
  TacsScalar E = 0.0, nu = 0.0;

  if (isotropic_flag) {
    // Set the values of the thickness, modulus and Poisson ratio
    rho = 2750.0e-9;  // kg/mm^3
    t = 1.0;          // mm
    E = 70e3;         // MPa
    nu = 0.3;         // dimensionless

    TacsScalar kcorr = 5.0 / 6.0, yield_stress = 464e6;
    stiff = new isoFSDTStiffness(rho, E, nu, kcorr, yield_stress, t);
    stiff->incref();
  } else {
    TacsScalar tply = 0.125;  // mm

    TacsScalar E1 = 164.0e3;
    TacsScalar E2 = 0.83e3;
    TacsScalar nu12 = 0.34;
    TacsScalar G12 = 2.1e3;
    TacsScalar G13 = 2.1e3;
    TacsScalar G23 = 1.2e3;

    TacsScalar Xt = 2410.0;
    TacsScalar Xc = 1040.0;
    TacsScalar Yt = 73.0;
    TacsScalar Yc = 173.0;
    TacsScalar S = 183.0;

    // Set the equivalent values of the thickness, modulus and Poisson
    // ratio
    rho = 1670.0e-9;  // kg/mm^3
    t = 4.0 * tply;   // mm
    E = E1;           // MPa
    nu = nu12;        // dimensionless

    OrthoPly *ortho_ply =
        new OrthoPly(tply, rho, E1, E2, nu12, G12, G23, G13, Xt, Xc, Yt, Yc, S);

    TacsScalar kcorr = 5.0 / 6.0;
    int num_plies = 4;
    OrthoPly *ortho_plies[4];
    ortho_plies[0] = ortho_plies[1] = ortho_plies[2] = ortho_plies[3] =
        ortho_ply;

    TacsScalar ply_angles[4] = {-M_PI / 4.0, M_PI / 4.0, M_PI / 4.0,
                                -M_PI / 4.0};
    TacsScalar thickness[4] = {0.125, 0.125, 0.125, 0.125};

    // Allocate the FSDTStiffness object
    stiff = new compFSDTStiffness(ortho_plies, kcorr, thickness, ply_angles,
                                  num_plies);
    stiff->incref();
  }

  // Set the reference axis for the beam
  TacsScalar ref_axis[] = {1.0, 0.0, 0.0};
  stiff->setRefAxis(ref_axis);

  // Estimate the buckling loads and geometric parameters for the plate
  double a = 100.0;
  double b = a * cos(theta);
  TacsScalar Dt = (t * t * t * E) / (12.0 * (1.0 - nu * nu));

#ifndef TACS_USE_COMPLEX
  // Extract the stiffness matrices from the constitutive object
  double pt[2] = {0.0, 0.0};
  TacsScalar A[6], B[6], D[6], As[3];
  stiff->getStiffness(pt, A, B, D, As);

  static const int nloads = 5;
  TacsScalar loads[nloads];

  static const int nfreq = 10;
  TacsScalar freq[nfreq];
  skewed_test(stiff, theta, a, b, Nx, Nxy, loads, nloads, freq, nfreq);

  if (rank == 0) {
    stiff->printStiffness();

    printf("Nx    = %8.1f\n", Nx);
    printf("Nxy   = %8.1f\n", Nxy);
    printf("theta = %8.1f\n", theta * 180.0 / M_PI);

    for (int k = 0; k < nloads; k++) {
      printf("load[%2d] = %15.6f  k_crit = %15.6f\n", k, loads[k],
             (loads[k] * a * a) / (M_PI * M_PI * Dt));
    }

    for (int k = 0; k < nfreq; k++) {
      printf("freq[%2d] = %15.6f   Omega = %15.6f\n", k, freq[k],
             freq[k] * a * a * sqrt(rho * t / Dt));
    }
  }
#endif

  TACSElement *skin = new MITCShell<4>(stiff, LINEAR, 0);

  const char *file_name = "skewed_plate.bdf";
  TACSMeshLoader *mesh = new TACSMeshLoader(MPI_COMM_WORLD);
  mesh->incref();
  mesh->scanBDFFile(file_name);
  mesh->setElement(0, skin);

  TACSAssembler *tacs = mesh->createTACS(6);
  tacs->incref();

  // Output for visualization
  int write_flag =
      (TACSElement::OUTPUT_NODES | TACSElement::OUTPUT_DISPLACEMENTS |
       TACSElement::OUTPUT_STRAINS | TACSElement::OUTPUT_STRESSES |
       TACSElement::OUTPUT_EXTRAS);

  TACSToFH5 *f5 = new TACSToFH5(tacs, TACS_SHELL, write_flag);
  f5->incref();

  int lev_fill = 5000;  // ILU(k) fill in
  int fill = 8.0;       // Expected number of non-zero entries

  // These calls compute the symbolic factorization and allocate
  // the space required for the preconditioners
  FEMat *aux_mat = tacs->createFEMat();
  PcScMat *pc = new PcScMat(aux_mat, lev_fill, fill, 1);
  aux_mat->incref();
  pc->incref();

  // Now, set up the solver
  int gmres_iters = 15;
  int nrestart = 0;     // Number of allowed restarts
  int is_flexible = 0;  // Is a flexible preconditioner?

  GMRES *ksm = new GMRES(aux_mat, pc, gmres_iters, nrestart, is_flexible);
  ksm->setTolerances(1e-12, 1e-30);
  ksm->incref();

  // Compute the linearized buckling mode
  int max_lanczos = 60;
  int num_eigvals = 5;
  double eig_tol = 1e-8;
  TacsScalar sigma = -120.0;

  FEMat *gmat = tacs->createFEMat();
  FEMat *kmat = tacs->createFEMat();

  int output_freq = 1;
  KSMPrint *ksm_print = new KSMPrintStdout("KSM", rank, output_freq);

  TACSLinearBuckling *linear_buckling = new TACSLinearBuckling(
      tacs, sigma, gmat, kmat, aux_mat, ksm, max_lanczos, num_eigvals, eig_tol);
  linear_buckling->incref();
  linear_buckling->solve(NULL, ksm_print);
  f5->writeToFile("results/load_path.f5");

  TACSBVec *vec = tacs->createVec();
  vec->incref();

  for (int k = 0; k < num_eigvals; k++) {
    TacsScalar error;
    TacsScalar eigvalue = linear_buckling->extractEigenvector(k, vec, &error);

    TacsScalar N_crit = 0.0;
    if (Nx < -0.1) {
      N_crit = 0.001 * t * E * eigvalue;
    } else {
      N_crit = 0.001 * t * 0.5 * E / (1.0 + nu) * eigvalue;
    }
    if (rank == 0) {
      printf("TACS eigs[%2d]: %15.6f k_crit = %15.6f\n", k,
             TacsRealPart(N_crit),
             TacsRealPart((N_crit * a * a) / (M_PI * M_PI * Dt)));
    }

    // Set the local variables
    tacs->setVariables(vec);
    char file_name[256];
    sprintf(file_name, "results/tacs_buckling_mode%02d.f5", k);
    f5->writeToFile(file_name);
  }

  sigma = 10.0;

  TACSFrequencyAnalysis *freq_analysis = new TACSFrequencyAnalysis(
      tacs, sigma, kmat, aux_mat, ksm, max_lanczos, num_eigvals, eig_tol);
  freq_analysis->incref();
  freq_analysis->solve(ksm_print);

  for (int k = 0; k < num_eigvals; k++) {
    TacsScalar error;
    TacsScalar eigvalue = freq_analysis->extractEigenvector(k, vec, &error);
    eigvalue = sqrt(eigvalue);

    // Set the local variables
    tacs->setVariables(vec);
    char file_name[256];
    sprintf(file_name, "results/tacs_mode%02d.f5", k);
    f5->writeToFile(file_name);

    if (rank == 0) {
      printf("TACS eigs[%2d]: %15.6f Omega = %15.6f\n", k,
             TacsRealPart(eigvalue),
             TacsRealPart(eigvalue * a * a * sqrt(rho * t / Dt)));
    }
  }

  linear_buckling->decref();
  freq_analysis->decref();

  vec->decref();
  mesh->decref();
  ksm->decref();
  pc->decref();
  kmat->decref();
  tacs->decref();
  stiff->decref();

  MPI_Finalize();
  return (0);
}
