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
//   Analyze a stiffened-panel model with a T-stiffener for some combined
//   loading. This uses the TACSPanelAnalysis code which forms a
//   finite-strip model of the panel for buckling analysis.

//   input:
//   skin:       the constitutive object for the skin
//   base:       the constitutive object for the base
//   stiffener:  the constitutive object for the stiffener
//   theta:      the skew angle
//   Nx:         the axial load on the panel
//   Nxy:        the shear load on the panel
//   use_lapack: use (or don't use) LAPACK routines for eigenvalue analysis
// */
// void panel_test( FSDTStiffness *skin, FSDTStiffness *base,
//                  FSDTStiffness *stiffener, TacsScalar theta,
//                  TacsScalar Nx, TacsScalar Nxy, int use_lapack ){
//   // The number of elements/nodes in the finite-strip model
//   int nrepeat = 4; // repeat this many segments
//   int nnodes = 10*nrepeat+1;
//   int nseg = 10*nrepeat;
//   int nbeams = 0;
//   int nmodes = 12;

//   // The baked-in dimensions of the panel in [mm]
//   TacsScalar Lx = 450.0;
//   TacsScalar b = 110.0;
//   TacsScalar hs = 20.0;
//   TacsScalar wb = 35.0;

//   // Create the TACSPanaelAnalysis object
//   TACSPanelAnalysis *panel = new TACSPanelAnalysis(nnodes, nseg,
//                                                    nbeams, nmodes, Lx,
//                                                    theta);
//   panel->incref();

//   // Set the size of the Lanczos subspace for the eigenvalue solver
//   panel->setLanczosSubspaceSize(100);

//   // Allocate an array to store the nodal locations
//   TacsScalar *Xpts = new TacsScalar[2*nnodes];
//   memset(Xpts, 0, 2*nnodes*sizeof(TacsScalar));

//   // Set the initial segment
//   int seg = 0, node = 0;

//   // For each repeating geometry segment (which incorporates a
//   // stiffener) set the nodal locations from the geometry parameters
//   //
//   for ( int k = 0; k < nrepeat; k++ ){
//     // Set the nodal locations
//     Xpts[2*node] = b*k;
//     Xpts[2*(node+1)] = b*k + (1.0/6.0)*(b - wb);
//     Xpts[2*(node+2)] = b*k + (1.0/3.0)*(b - wb);
//     Xpts[2*(node+3)] = b*k + 0.5*(b - wb);
//     Xpts[2*(node+4)] = 0.5*b*(2*k + 1);

//     Xpts[2*(node+5)] = 0.5*b*(2*k + 1);
//     Xpts[2*(node+5)+1] = - hs/2.0;

//     Xpts[2*(node+6)] = 0.5*b*(2*k + 1);
//     Xpts[2*(node+6)+1] = - hs;

//     Xpts[2*(node+7)] = 0.5*b*(2*k + 1) + 0.5*wb;
//     Xpts[2*(node+8)] = 0.5*b*(2*k + 1) + 0.5*wb + (1.0/6.0)*(b - wb);
//     Xpts[2*(node+9)] = 0.5*b*(2*k + 1) + 0.5*wb + (1.0/3.0)*(b - wb);
//     Xpts[2*(node+10)] = 0.5*b*(2*k + 1) + 0.5*wb + 0.5*(b - wb);

//     // Set the connectivity and set the constitutive objects
//     panel->setSegment(seg, TACSPanelAnalysis::SKIN_SEGMENT,
//                       skin, node, node+1);
//     panel->setSegment(seg+1, TACSPanelAnalysis::SKIN_SEGMENT,
//                       skin, node+1, node+2);
//     panel->setSegment(seg+2, TACSPanelAnalysis::SKIN_SEGMENT,
//                       skin, node+2, node+3);
//     panel->setSegment(seg+3, TACSPanelAnalysis::SKIN_SEGMENT,
//                       base, node+3, node+4);

//     panel->setSegment(seg+4, TACSPanelAnalysis::STIFFENER_SEGMENT,
//                       stiffener, node+4, node+5);
//     panel->setSegment(seg+5, TACSPanelAnalysis::STIFFENER_SEGMENT,
//                       stiffener, node+5, node+6);

//     panel->setSegment(seg+6, TACSPanelAnalysis::SKIN_SEGMENT,
//                       base, node+4, node+7);
//     panel->setSegment(seg+7, TACSPanelAnalysis::SKIN_SEGMENT,
//                       skin, node+7, node+8);
//     panel->setSegment(seg+8, TACSPanelAnalysis::SKIN_SEGMENT,
//                       skin, node+8, node+9);
//     panel->setSegment(seg+9, TACSPanelAnalysis::SKIN_SEGMENT,
//                       skin, node+9, node+10);

//     node += 10;
//     seg += 10;
//   }

//   // Set the final node locations
//   Xpts[2*node] = b*nrepeat;

//   // Set the node locations into the panel object
//   panel->setPoints(Xpts, nnodes);

//   // Set which boundary conditions to use within the model
//   int bc = (4 | 8);
//   panel->setFirstNodeBC(0, bc);
//   panel->setLastNodeBC(nnodes-1, bc);
//   panel->initialize();

//   // Set the flag for whether to use LAPACK or not
//   panel->setUseLapackEigensolver(use_lapack);

//   // Compute all of the loads and record the time
//   int nloads = 15;
//   TacsScalar pos_loads[20], neg_loads[20];
//   double t0 = 0.0;

//   if (Nx == 0.0){
//     t0 = MPI_Wtime();
//     panel->computeBucklingLoads(Nx, Nxy,
//                                 pos_loads, nloads, "results/pos_");
//     panel->computeBucklingLoads(Nx, -Nxy,
//                                 neg_loads, nloads, "results/neg_");
//     t0 = MPI_Wtime() - t0;

//     for ( int k = 0; k < nloads; k++ ){
//       printf("pos_load[%2d] = %25.12f\n", k, pos_loads[k]);
//     }
//     for ( int k = 0; k < nloads; k++ ){
//       printf("neg_load[%2d] = %25.12f\n", k, neg_loads[k]);
//     }
//   }
//   else {
//     t0 = MPI_Wtime();
//     panel->computeBucklingLoads(Nx, Nxy,
//                                 pos_loads, nloads, "results/");
//     t0 = MPI_Wtime() - t0;

//     for ( int k = 0; k < nloads; k++ ){
//       printf("load[%2d] = %25.12f\n", k, pos_loads[k]);
//     }
//   }

//   printf("Solution time: %f\n", t0);

//   delete [] Xpts;
//   panel->decref();
// }

// #endif // TACS_USE_COMPLEX

/*
  The following code compares the results of a full buckling
  eigenvalue computation in TACS to the approximate calculations
  performed using the TACSPanelAnalysis class. This calculation uses a
  finite-strip method which restricts the type of panel geometries
  that can be analyzed.

  The useage of this script is the following:

  ./stiffened_panel [options]

  where [options] consist of the following:

  dh:     the step length size
  shear:  use the input file for shear
  lapack: use LAPACK routines for the TACSPanelAnalysis

  Note:

  1. You must generate a .bdf file first by using the other code in
  this directory.
  2. The results will be written to files in the ./results/ directory.
  If no such directory exists, then no results will be written.
*/
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // Depending on the arguments, load in a different
  const char *bdf_file = "axial_stiffened_panel.bdf";
  const char *shear_file = "shear_stiffened_panel.bdf";

  // Interpret the input flags
  double dh = 1e-6;
  int shear_flag = 0;
  int use_lapack = 0;
  for (int k = 0; k < argc; k++) {
    if (strcmp(argv[k], "shear") == 0) {
      bdf_file = shear_file;
      shear_flag = 1;
    }
    if (strcmp(argv[k], "lapack") == 0) {
      use_lapack = 1;
    }
    if (sscanf(argv[k], "fd=%le", &dh) == 1) {
      if (dh > 0.1) {
        dh = 0.1;
      }
      printf("Using difference step size: %le\n", dh);
    }
  }

  // Get processor rank
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // The thickness of the skin, base and stiffener
  TacsScalar tskin = 1.0;
  TacsScalar tbase = tskin;
  TacsScalar tstiff = 1.2 * tskin;

  // The material properties used for the skin/base and stiffener
  TacsScalar rho = 2750.0, E = 70.0e3, nu = 0.3;
  TacsScalar specific_heat = 921.096;
  TacsScalar yield_stress = 464e3;
  TacsScalar cte = 24.0e-6;
  TacsScalar kappa = 230.0;

  TACSMaterialProperties *props = new TACSMaterialProperties(
      rho, specific_heat, E, nu, yield_stress, cte, kappa);

  TacsScalar axis[] = {1.0, 0.0, 0.0};
  TACSShellTransform *transform = new TACSShellRefAxisTransform(axis);

  // Create the stiffness objects for the skin/base and stiffener
  TACSIsoShellConstitutive *stiff_skin =
      new TACSIsoShellConstitutive(props, tskin, 0);
  TACSIsoShellConstitutive *stiff_base =
      new TACSIsoShellConstitutive(props, tbase, 1);
  TACSIsoShellConstitutive *stiff_stiffener =
      new TACSIsoShellConstitutive(props, tstiff, 2);

  // Allocate the elements associated with the skin/stiffener/base
  TACSElement *skin = NULL, *base = NULL, *stiffener = NULL;
  skin = TacsCreateShellByName("TACSQuad16Shell", transform, stiff_skin);
  base = TacsCreateShellByName("TACSQuad16Shell", transform, stiff_base);
  stiffener =
      TacsCreateShellByName("TACSQuad16Shell", transform, stiff_stiffener);

  // Set the loading conditions which depends on whether or not a
  // shear flag is used
  TacsScalar Nx = -1.0, Nxy = 0.0;
  if (shear_flag) {
    Nx = 0.0;
    Nxy = 1.0;
  }

  //   // When using complex mode, we cannot use the tacs panel analysis
  // #ifndef TACS_USE_COMPLEX
  //   TacsScalar theta = -15.0/180.0*M_PI;
  //   if (rank == 0){
  //     printf("theta = %8.1f\n", theta*180.0/M_PI);
  //     panel_test(stiff_skin, stiff_base, stiff_stiffener, theta,
  //                Nx, Nxy, use_lapack);
  //   }
  // #endif // TACS_USE_COMPLEX

  // Load in the .bdf file using the TACS mesh loader
  TACSMeshLoader *mesh = new TACSMeshLoader(MPI_COMM_WORLD);
  mesh->incref();

  // Scan the file
  mesh->scanBDFFile(bdf_file);

  // Set the skin/base/stiffener elements
  mesh->setElement(0, skin);
  mesh->setElement(1, base);
  mesh->setElement(2, stiffener);

  // Create the TACSAssembler object
  int vars_per_node = 6;
  TACSAssembler *assembler = mesh->createTACS(vars_per_node);
  assembler->incref();

  // Output for visualization
  int write_flag = (TACS_OUTPUT_CONNECTIVITY | TACS_OUTPUT_NODES |
                    TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                    TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
  TACSToFH5 *f5 =
      new TACSToFH5(assembler, TACS_BEAM_OR_SHELL_ELEMENT, write_flag);
  f5->incref();

  int lev_fill = 5000;  // ILU(k) fill in
  int fill = 8.0;       // Expected number of non-zero entries

  // These calls compute the symbolic factorization and allocate
  // the space required for the preconditioners
  TACSSchurMat *aux_mat = assembler->createSchurMat();
  TACSSchurPc *pc = new TACSSchurPc(aux_mat, lev_fill, fill, 1);
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
  int max_lanczos = 50;
  int neigvals = 4;
  double eig_tol = 1e-8;
  TacsScalar sigma = 0.15;

  // Allocate matrices for the stiffness/geometric stiffness matrix
  TACSSchurMat *kmat = assembler->createSchurMat();
  TACSSchurMat *gmat = assembler->createSchurMat();

  // Create a print object for writing output to the screen
  int freq = 1;
  KSMPrint *ksm_print = new KSMPrintStdout("KSM", rank, freq);

  // Allocate the linear buckling analysis object
  TACSLinearBuckling *linear_buckling =
      new TACSLinearBuckling(assembler, sigma, gmat, kmat, aux_mat, ksm,
                             max_lanczos, neigvals, eig_tol);
  linear_buckling->incref();
  linear_buckling->solve(NULL, ksm_print);
  f5->writeToFile("results/load_path.f5");
  linear_buckling->checkEigenvector(0);

  // Extract the eigenvalue
  TacsScalar error = 0.0, eigvalue = 0.0;
  eigvalue = linear_buckling->extractEigenvalue(0, &error);

  //   // Compute the derivative of the eigenvalue
  //   TacsScalar *fdvSens = new TacsScalar[ ndvs ];
  //   linear_buckling->evalEigenDVSens(0, fdvSens, ndvs);

  //   // Now, compute the FD approximation
  //   TacsScalar *x = new TacsScalar[ ndvs ];
  //   tacs->getDesignVars(x, ndvs);

  //   // Compute the projected derivative. Set the direction such that it
  //   // lies along the positive gradient component directions.
  //   TacsScalar proj = 0.0;
  //   for ( int i = 0; i < ndvs; i++ ){
  //     proj += fabs(fdvSens[i]);
  //   }

  // #ifdef TACS_USE_COMPLEX
  //   // Use a complex-step perturbation
  //   for ( int i = 0; i < ndvs; i++ ){
  //     if (TacsRealPart(fdvSens[i]) > 0.0){
  //       x[i] = x[i] + TacsScalar(0.0, dh);
  //     }
  //     else {
  //       x[i] = x[i] - TacsScalar(0.0, dh);
  //     }
  //   }
  // #else
  //   // Use finite-difference perturbation
  //   for ( int i = 0; i < ndvs; i++ ){
  //     if (fdvSens[i] > 0.0){
  //       x[i] = x[i] + dh;
  //     }
  //     else {
  //       x[i] = x[i] - dh;
  //     }
  //   }
  // #endif // TACS_USE_COMPLEX

  //   // Set the new design variable values
  //   tacs->setDesignVars(x, ndvs);

  //   // Solve the buckling problem again
  //   linear_buckling->solve(NULL, ksm_print);
  //   TacsScalar eigvalue1 =
  //     linear_buckling->extractEigenvalue(0, &error);

  //   // Evaluate the finite-difference or complex-step approximation
  // #ifdef TACS_USE_COMPLEX
  //   TacsScalar pfd = TacsImagPart(eigvalue1)/dh;
  // #else
  //   TacsScalar pfd = (eigvalue1 - eigvalue)/dh;
  // #endif // TACS_USE_COMPLEX

  //   // Write out the projected error
  //   if (rank == 0){
  //     printf("%15s %15s %15s %15s\n",
  //            "<p, df/dx>", "FD", "Err", "Rel err");
  //     printf("%15.8e %15.8e %15.8e %15.8e\n",
  //            TacsRealPart(proj), TacsRealPart(pfd), TacsRealPart(proj - pfd),
  //            TacsRealPart((proj - pfd)/proj));
  //   }

  //   delete [] fdvSens;
  //   delete [] x;

  //   // Write out the results to a file. Iterate through the eigenvectors
  //   // and print out each one using the TACSToFH5 converter
  //   TACSBVec *vec = tacs->createVec();
  //   vec->incref();

  //   for ( int k = 0; k < neigvals; k++ ){
  //     TacsScalar error = 0.0, eigvalue = 0.0;
  //     eigvalue = linear_buckling->extractEigenvector(k, vec, &error);

  //     double Lx = 450.0;
  //     double b = 110.0;
  //     double hs = 20.0;
  //     double wb = 35.0;

  //     if (shear_flag){
  //       TacsScalar Nxy_crit = G*((tskin*(b - wb) +
  //                                 tbase*wb)/b)*(eigvalue/Lx);

  //       if (rank == 0){
  //         printf("TACS eigs[%2d]: %15.6f\n", k, TacsRealPart(Nxy_crit));
  //       }
  //     }
  //     else {
  //       TacsScalar Nx_crit = E*((tskin*(b - wb) +
  //                                tbase*wb + tstiff*hs)/b)*(eigvalue/Lx);

  //       if (rank == 0){
  //         printf("TACS eigs[%2d]: %15.6f\n", k, TacsRealPart(Nx_crit));
  //       }
  //     }

  //     // Set the local variables
  //     tacs->setVariables(vec);
  //     char file_name[256];
  //     sprintf(file_name, "results/tacs_buckling_mode%02d.f5", k);
  //     f5->writeToFile(file_name);
  //   }

  //   linear_buckling->decref();
  //   vec->decref();
  //   mesh->decref();
  //   ksm->decref();
  //   pc->decref();
  //   tacs->decref();

  MPI_Finalize();
  return (0);
}
