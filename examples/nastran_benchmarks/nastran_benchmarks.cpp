/*
Test case to model a built-up structure (beams and shells) taken from a .bdf 
input file. Stiffened plate will be clamped on all edges. There will be two differnt 
cases: static and dynamic. Both cases will have a single point load in the 
transverse direction placed at the center of the plate. In the dynamic case, the
load will be applied as a step. The static case will use GMRES to solve the system
of equations. The dynamic case will use a choice of either BDF or DIRK for the 
time integration solution.
*/

#include "TACSMeshLoader.h"
#include "MITCShell.h"
#include "EBBeam.h"
#include "isoFSDTStiffness.h"
#include "rectangleEBStiffness.h"
#include "KSM.h"
#include "TACSIntegrator.h"
#include <string>
#include <sstream>

int main( int argc, char **argv) {
  // Intialize MPI and declare communicator
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  double t1, t2;
  t1 = MPI_Wtime();

  // Problem Setup
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // The name of BDF file to be loaded
  const char *filename = "stiff_plate_160x160.bdf"; // Use stiff_plate_NxN.bdf/Plate_NxN.bdf or Plate.bdf for static and dynamic tests, respectively 

  // Create the mesh loader object and load file
  TACSMeshLoader *mesh = new TACSMeshLoader(comm);
  mesh->incref();
  mesh->scanBDFFile(filename);

  // Get number of components prescribed in BDF file
  int num_components = mesh->getNumComponents();

  // Set properties needed to create stiffness object for shells (6061-T6 aluminum)
  double rho_s = 2700.0; // kg/m^3
  double E_s = 68.9e9; // Pa
  double nu_s = 0.33; 
  double kcorr_s = 5.0/6.0;
  double ys_s = 276e6; // Pa
  double thickness = 0.005; // m

  isoFSDTStiffness *stiff_s = new isoFSDTStiffness(rho_s, E_s, nu_s, kcorr_s, ys_s, thickness);

  // Set properties needed to create stiffness object for beams (A36 steel),
  // Rectangular section: 5 cm x 5 cm, offset so it supports bottom side of plate
  TacsScalar rho_b = 7850.0; // kg/m^3
  TacsScalar E_b = 200e9; // Pa
  TacsScalar nu_b = 0.26;
  TacsScalar G_b = E_b/(2*(1+nu_b)); // Pa
  TacsScalar ys_b = 250e6; // Pa
  TacsScalar h = 0.05; // m
  int height_num = 0;
  int thick_num = 0;
  TacsScalar axis_b[] = {0.0, 0.0, 1.0};
  TacsScalar z_offset = -0.0275;
  TacsScalar y_offset = 0;

  rectangleEBStiffness *stiff_b = new rectangleEBStiffness(rho_b, E_b, G_b, ys_b, h, h, height_num, thick_num, axis_b, z_offset, y_offset);
 
  // Loop over components, creating elements for each 
  for (int i = 0; i < num_components; i++) {

  	const char *descriptor = mesh->getElementDescript(i);

  	// Initialize the element object
  	TACSElement *element = NULL;

  	if (strcmp(descriptor, "CQUAD4") == 0) {
      element = new MITCShell<2>(stiff_s, LINEAR, i); 
    }
    else {
    	element = new EBBeam(stiff_b);
    }
    mesh->setElement(i, element);
  }

  // Create TACS assembler from mesh loader object
  TACSAssembler::OrderingType order_type = TACSAssembler::NATURAL_ORDER;
  TACSAssembler::MatrixOrderingType mat_type = TACSAssembler::DIRECT_SCHUR;
  TACSAssembler *tacs = mesh->createTACS(6, order_type, mat_type); 
  tacs->incref();
  mesh->decref();
  printf("The mesh has been created successfully\n");

  t2 = MPI_Wtime();

  // Static Solution
  ///////////////////////////////////////////////////////////////////////////////////////////////

  // Solve the system
  int solve_static = 0; // Switch to solve either static or dynamic case

  if (solve_static) {
    // Solve static response

    // On each partition, search the nodes to find which one has the center node
    TACSVarMap *varMap = tacs->getVarMap();
    int center_node_global = 3280; // Center Nodes for Plate Meshes: 61, 221, 841, 3281, 12961, 51521, 205441
    int node_on_proc = varMap->getOwner(center_node_global);

    // Create the force vector and specify load
    double force = -10.0;
    TACSBVec *f = tacs->createVec(); // loads
    f->incref();
    TacsScalar *force_vals;
    int size_f = f->getArray(&force_vals);
    printf("Size of the solution vector on processor rank [%i] is: %i \n",rank, size_f);

    // Find center node in the correct local force vector and apply load
    if (rank == node_on_proc){
      printf("Center node is on processor of rank [%i] \n",node_on_proc);
      const int *ownerRange = NULL;
      varMap->getOwnerRange(&ownerRange);
      printf("The node range for processor of rank [%i] is: %i - %i\n",node_on_proc,ownerRange[rank],ownerRange[rank+1]);
      int start = ownerRange[rank];
      int center_node_local = center_node_global - start;
      int dof = center_node_local*6+2;
      printf("center_node_local d.o.f.: %i\n", dof);
      force_vals[dof] = force;
      // tacs->applyBCs(f);
    }
    // Serial method to apply force
    // for ( int k = 2; k < size_f; k += 6 ){
    //   if (k==(center_node_global*6)+2){
    //     force_vals[k] += force;
    //   }
    // }

    tacs->applyBCs(f);

    // Create matrix and vectors 
    TACSBVec *ans = tacs->createVec(); // displacements and rotations
    // TACSBVec *res = tacs->createVec(); // The residual
    FEMat *mat = tacs->createFEMat(); // stiffness matrix

    // Increment reference count to the matrix/vectors
    ans->incref();
    // res->incref();
    mat->incref();

    // Allocate the factorization (preconditioner)
    int lev = 10000; //10000
    double fill = 10.0; //10.0
    int reorder_schur = 1;
    PcScMat *pc = new PcScMat(mat, lev, fill, reorder_schur); 
    pc->incref();

    // Create GMRES object
    int gmres_iters = 100, nrestart = 2, isflexible = 0; //gmres_iters = 100
    TACSKsm * ksm = new GMRES(mat, pc, gmres_iters, nrestart, isflexible);
    double rtol = 1e-16, atol = 1e-30;
    ksm->setTolerances(rtol, atol);
    int freq = 1;
    KSMPrint *monitor = new KSMPrintStdout( "GMRES", rank, freq );
    ksm->setMonitor(monitor);
    ksm->incref();

    // Assemble and factor the stiffness/Jacobian matrix. Factor the
    // Jacobian and solve the linear system for the displacements
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    tacs->assembleJacobian(alpha, beta, gamma, NULL, mat); // might change NULL to res
    pc->factor(); // LU factorization of stiffness matrix

    pc->applyFactor(f, ans);
    ksm->solve(f, ans); 
    tacs->setVariables(ans);
    // printf("The static solution has been solved for\n");

    // Create an TACSToFH5 object for writing output to files
    unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                               TACSElement::OUTPUT_DISPLACEMENTS |
                               TACSElement::OUTPUT_STRAINS |
                               TACSElement::OUTPUT_STRESSES |
                               TACSElement::OUTPUT_EXTRAS);
    TACSToFH5 * f5 = new TACSToFH5(tacs, TACS_SHELL, write_flag);
    f5->incref();
    f5->writeToFile("panel_static.f5");
    f5->decref();

    // Decref the matrix/pc objects
    mat->decref();
    pc->decref();
    ksm->decref();

    // Decref the vectors
    ans->decref();
    f->decref();
    // res->decref();

    // Decref TACS
    tacs->decref();
  }

  // Dynamic Solution
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  
  else {
    // Solve dynamic response

    // Create matrices and vectors
    TACSBVec *f = tacs->createVec(); // loads

    // Increment reference count to the matrix/vectors
    f->incref();

    // Set the time integration parameters 
    double ti = 0.0; // Initial time
    double tf = 0.5; // Final time
    double num_steps = 50000; // The number of total steps 

    // Switch to decide time integration scheme
    int use_BDF = 1;

    if (use_BDF) {

      int BDF_order = 2;
      TACSIntegrator *integrator = new TACSBDFIntegrator(tacs, ti, tf, num_steps, BDF_order);
      integrator->incref();

      // Set the integrator options
      integrator->setUseFEMat(1, TACSAssembler::TACS_AMD_ORDER);
      integrator->setAbsTol(1e-7);
      integrator->setOutputFrequency(10);
      integrator->setPrintLevel(2);
      printf("Integrator set up successfully\n");

      // Create file that will take displacements at center node for each time step
      char filename[128];
      sprintf(filename,"center_point_displacements.dat");
      FILE *fp = fopen(filename, "w");
      fprintf(fp, "Variables = t, u0, v0, w0\n");
      
      // Set up displacements and step force parameters
      TACSBVec *q = NULL; // vector to hold displacements
      int force_step = floor(0.2*num_steps); // time step that the step-force is "turned on"
      TacsScalar *force_vals;
      double mag = 10.0; // Magnitude of step force
      int size_f = f->getArray(&force_vals);
      printf("The size of the problem is: %i \n", size_f);

      double tref = t2-t1;

      // Solve
      for (int step = 0; step < num_steps+1; step++){
        if (step <= force_step) {
          // No force applied while step force is "turned off"
          tacs->applyBCs(f);

          integrator->iterate(step, f); 

          // Get center point data and write to a file for plotting later
          TacsScalar X[3*4], vars[6*4]; // X - shells have 4 nodes with 3 nodal locations each; vars - each node has 6 variables defined at it
          double time = integrator->getStates(step, &q, NULL, NULL);
          tacs->setVariables(q);
          TACSElement *element = tacs->getElement(12720, X, vars); // Element number will change depending on mesh size (elem 3160 for 80x80 mesh)
          // Write to file the displacements of each node of the element
          fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                  time, vars[0], vars[1], vars[2],
                        vars[6], vars[7], vars[8],
                        vars[12], vars[13], vars[14],
                        vars[18], vars[19], vars[20]); // Only interested in the u, v, w variables at each node
        }
        else {
          // The step force is now "turned on"
          for ( int k = 2; k < size_f; k += 6 ){
            if (k == (12960*6)+2) { // Node number changes with mesh (19682 for center node of 80x80 mesh, z-force)
              force_vals[k] = -1*mag;
            } 
          }
          tacs->applyBCs(f);

          integrator->iterate(step, f);

          // Get center point data and write to a file for plotting later
          TacsScalar X[3*4], vars[6*4];
          double time = integrator->getStates(step, &q, NULL, NULL);
          tacs->setVariables(q);
          TACSElement *element = tacs->getElement(12720, X, vars);
          fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                  time, vars[0], vars[1], vars[2],
                        vars[6], vars[7], vars[8],
                        vars[12], vars[13], vars[14],
                        vars[18], vars[19], vars[20]);
        }
      
        // // Each time step write a new .f5 files so you can view the full time history response
        // unsigned int write_flag = (TACSElement::OUTPUT_NODES |
        //                            TACSElement::OUTPUT_DISPLACEMENTS |
        //                            TACSElement::OUTPUT_STRAINS |
        //                            TACSElement::OUTPUT_STRESSES |
        //                            TACSElement::OUTPUT_EXTRAS);
        // TACSToFH5 * f5 = new TACSToFH5(tacs, TACS_SHELL, write_flag);
        // f5->incref();
        // char fname[128];
        // sprintf(fname,"panel_dynamic_%d.f5",step); // "stiffened_panel_dynamic_%d.f5" for stiffened case, "panel_dynamic_%d.f5" for no stiffeners
        // f5->writeToFile(fname); //"stiffened_panel.f5"

        // integrator->setShellOutput(f5); 

        // integrator->writeSolutionToF5();
      }
      fclose(fp);
      integrator->printWallTime(tref, 0);
    }

    else { 
      int DIRK_stages = 3;
      TACSIntegrator *integrator = new TACSDIRKIntegrator(tacs, ti, tf, num_steps, DIRK_stages);
      integrator->incref();

      // Set the integrator options
      integrator->setUseFEMat(1, TACSAssembler::TACS_AMD_ORDER);
      integrator->setAbsTol(1e-7);
      integrator->setOutputFrequency(10);
      integrator->setPrintLevel(2);

      // Create file that will take displacements at center node for each time step
      char filename[128];
      sprintf(filename,"center_point_displacements.dat");
      FILE *fp = fopen(filename, "w");
      fprintf(fp, "Variables = t, u0, v0, w0\n");

      // Set up displacements and step force parameters
      TACSBVec *q = NULL; // vector to hold displacements
      int force_step = floor(0.2*num_steps); // time step that the step-force is "turned on"
      TacsScalar *force_vals;
      double mag = 10.0; // Magnitude of step force
      int size_f = f->getArray(&force_vals);
      printf("The size of the problem is: %i \n", size_f);
      
      // Solve
      for (int step = 0; step < num_steps+1; step++){
        if (step <= force_step) {
          // No force applied while step force is "turned off"
          tacs->applyBCs(f);

          integrator->iterate(step, f); 

          // Get center point data and write to a file for plotting later
          TacsScalar X[3*4], vars[6*4]; // X - shells have 4 nodes with 3 nodal locations each; vars - each node has 6 variables defined at it
          double time = integrator->getStates(step, &q, NULL, NULL);
          tacs->setVariables(q);
          TACSElement *element = tacs->getElement(3159, X, vars); // Element number will change depending on mesh size (elem 3160 for 80x80 mesh)
          // Write to file the displacements of each node of the element
          fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                  time, vars[0], vars[1], vars[2],
                        vars[6], vars[7], vars[8],
                        vars[12], vars[13], vars[14],
                        vars[18], vars[19], vars[20]); // Only interested in the u, v, w variables at each node
        }
        else {
          // The step force is now "turned on"
          for ( int k = 2; k < size_f; k += 6 ){
            if (k == 19682) { // Node number changes with mesh (19682 for center node of 80x80 mesh, z-force)
              force_vals[k] = -1*mag;
            } 
          }
          tacs->applyBCs(f);

          integrator->iterate(step, f);

          // Get center point data and write to a file for plotting later
          TacsScalar X[3*4], vars[6*4];
          double time = integrator->getStates(step, &q, NULL, NULL);
          tacs->setVariables(q);
          TACSElement *element = tacs->getElement(3159, X, vars);
          fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                  time, vars[0], vars[1], vars[2],
                        vars[6], vars[7], vars[8],
                        vars[12], vars[13], vars[14],
                        vars[18], vars[19], vars[20]);
        }
      
        // // Each time step write a new .f5 files so you can view the full time history response
        // unsigned int write_flag = (TACSElement::OUTPUT_NODES |
        //                            TACSElement::OUTPUT_DISPLACEMENTS |
        //                            TACSElement::OUTPUT_STRAINS |
        //                            TACSElement::OUTPUT_STRESSES |
        //                            TACSElement::OUTPUT_EXTRAS);
        // TACSToFH5 * f5 = new TACSToFH5(tacs, TACS_SHELL, write_flag);
        // f5->incref();
        // char fname[128];
        // sprintf(fname,"panel_dynamic_%d.f5",step); // "stiffened_panel_dynamic_%d.f5" for stiffened case, "panel_dynamic_%d.f5" for no stiffeners
        // f5->writeToFile(fname); //"stiffened_panel.f5"

        // integrator->setShellOutput(f5); 

        // integrator->writeSolutionToF5();
      }
    fclose(fp);
  }

}

  MPI_Finalize();

  return 0;
}