#include "TACSCreator.h"
#include "TACSElement2D.h"
#include "TACSFH5Loader.h"
#include "TACSLinearElasticity.h"
#include "TACSThermoelasticity.h"
#include "TACSToFH5.h"
#include "TACSTriangularBasis.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // Check whether to use elasticity or thoermoelasticity
  int use_thermoelasticity = 0;
  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "thermoelasticity") == 0) {
      use_thermoelasticity = 1;
    }
  }

  // Get the rank of the processor
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int vars_per_node = 2;
  if (use_thermoelasticity) {
    vars_per_node = 3;
  }

  // Allocate the TACS creator
  TACSCreator *creator = new TACSCreator(MPI_COMM_WORLD, vars_per_node);
  creator->incref();

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

  // Create the stiffness object
  TACSPlaneStressConstitutive *stiff = new TACSPlaneStressConstitutive(props);
  stiff->incref();

  // Create the model class
  TACSElementModel *model = NULL;
  if (use_thermoelasticity) {
    model = new TACSLinearThermoelasticity2D(stiff, TACS_LINEAR_STRAIN);
  } else {
    model = new TACSLinearElasticity2D(stiff, TACS_LINEAR_STRAIN);
  }
  model->incref();

  // Create the basis
  TACSElementBasis *basis = new TACSQuadraticTriangleBasis();

  // Create the element type
  TACSElement2D *elem = new TACSElement2D(model, basis);
  elem->incref();

  // Only set the mesh/boundary conditions etc. on the
  // root processor
  if (rank == 0) {
    // Create a regular mesh of triangular elements
    int nx = 20, ny = 20;
    int num_elements = 10 * nx * ny;

    // Allocate the two halfs of the mesh
    int lena = (4 * nx + 1) * (4 * ny + 1);
    int lenb = (2 * nx + 1) * (2 * ny + 1);
    int *anodes = new int[lena];
    int *bnodes = new int[lenb];

    memset(anodes, 0, lena * sizeof(int));
    for (int i = 0; i < lenb; i++) {
      bnodes[i] = -1;
    }

    int dep_nodes = 0;
    // Set the dependent nodes
    for (int j = 1; j < 4 * ny + 1; j += 2) {
      anodes[4 * nx + (4 * nx + 1) * j] = -dep_nodes - 1;
      dep_nodes++;
    }

    // Number the nodes for the a-block
    int num_nodes = 0;
    for (int j = 0; j < 4 * ny + 1; j++) {
      for (int i = 0; i < 4 * nx + 1; i++) {
        int index = i + (4 * nx + 1) * j;
        if (anodes[index] >= 0) {
          anodes[i + (4 * nx + 1) * j] = num_nodes;
          num_nodes++;
        }
      }
    }

    // Copy over the block indices for the second block
    for (int j = 0; j < 2 * ny + 1; j++) {
      bnodes[(2 * nx + 1) * j] = anodes[(4 * nx + 1) * (2 * j + 1) - 1];
    }

    // Now number all the negative block numbers
    for (int j = 0; j < 2 * ny + 1; j++) {
      for (int i = 0; i < 2 * nx + 1; i++) {
        int index = i + (2 * nx + 1) * j;
        if (bnodes[index] <= 0) {
          bnodes[index] = num_nodes;
          num_nodes++;
        }
      }
    }

    // Now allocate the dependent node data structures
    int *dep_ptr = new int[dep_nodes + 1];
    int *dep_conn = new int[3 * dep_nodes];
    double *dep_weights = new double[3 * dep_nodes];

    // Set the dependent weights
    dep_ptr[0] = 0;
    for (int j = 0; j < 2 * ny; j++) {
      if (j % 2 == 0) {
        dep_weights[3 * j] = 3.0 / 8.0;
        dep_weights[3 * j + 1] = 3.0 / 4.0;
        dep_weights[3 * j + 2] = -1.0 / 8.0;
      } else {
        dep_weights[3 * j] = -1.0 / 8.0;
        dep_weights[3 * j + 1] = 3.0 / 4.0;
        dep_weights[3 * j + 2] = 3.0 / 8.0;
      }

      // Set the dependent node connectivity
      int jj = j / 2;

      dep_conn[3 * j] = anodes[4 * nx + 4 * jj * (4 * nx + 1)];
      dep_conn[3 * j + 1] = anodes[4 * nx + (4 * jj + 2) * (4 * nx + 1)];
      dep_conn[3 * j + 2] = anodes[4 * nx + (4 * jj + 4) * (4 * nx + 1)];

      // Set the dependent nodes
      dep_ptr[j + 1] = dep_ptr[j] + 3;
    }

    // Set the lengths along the x and y directions
    double Lx = 10.0, Ly = 10.0;

    // Create the nodes
    TacsScalar *Xpts = new TacsScalar[3 * num_nodes];

    // Loop over the a-nodes
    for (int j = 0; j < 4 * ny + 1; j++) {
      for (int i = 0; i < 4 * nx + 1; i++) {
        int node = anodes[i + (4 * nx + 1) * j];
        if (node >= 0) {
          Xpts[3 * node] = 0.25 * Lx * i / nx;
          Xpts[3 * node + 1] = 0.25 * Ly * j / ny;
          Xpts[3 * node + 2] = 0.0;
        }
      }
    }

    for (int j = 0; j < 2 * ny + 1; j++) {
      for (int i = 0; i < 2 * nx + 1; i++) {
        int node = bnodes[i + (2 * nx + 1) * j];
        if (node >= 0) {
          Xpts[3 * node] = Lx + 0.5 * Lx * i / nx;
          Xpts[3 * node + 1] = 0.5 * Ly * j / ny;
          Xpts[3 * node + 2] = 0.0;
        }
      }
    }

    // Set up the element connectivity arrays
    int *elem_node_conn = new int[6 * num_elements];
    int *elem_node_ptr = new int[num_elements + 1];

    int n = 0;
    int *conn = elem_node_conn;
    elem_node_ptr[0] = 0;

    // (2*i, 2*j+2) -- (2*i+1, 2*j+2) -- (2*i+2, 2*j+2)
    //     |                        /
    // (2*i, 2*j+1)    (2*i+1, 2*j+1)    (2*i+2, 2*j+1)
    //     |        /                       |
    // (2*i, 2*j) --   (2*i+1, 2*j) ---- (2*i+2, 2*j)

    // Add the elements from the a-block
    int N = 4 * nx + 1;
    for (int j = 0; j < 2 * ny; j++) {
      for (int i = 0; i < 2 * nx; i++) {
        conn[0] = anodes[2 * i + 2 + N * (2 * j)];
        conn[1] = anodes[2 * i + 2 + N * (2 * j + 2)];
        conn[2] = anodes[2 * i + N * (2 * j)];
        conn[3] = anodes[2 * i + 2 + N * (2 * j + 1)];
        conn[4] = anodes[2 * i + 1 + N * (2 * j + 1)];
        conn[5] = anodes[2 * i + 1 + N * (2 * j)];
        conn += 6;
        elem_node_ptr[n + 1] = elem_node_ptr[n] + 6;
        n++;

        conn[0] = anodes[2 * i + N * (2 * j + 2)];
        conn[1] = anodes[2 * i + N * (2 * j)];
        conn[2] = anodes[2 * i + 2 + N * (2 * j + 2)];
        conn[3] = anodes[2 * i + N * (2 * j + 1)];
        conn[4] = anodes[2 * i + 1 + N * (2 * j + 1)];
        conn[5] = anodes[2 * i + 1 + N * (2 * j + 2)];
        conn += 6;
        elem_node_ptr[n + 1] = elem_node_ptr[n] + 6;
        n++;
      }
    }

    // Add the elements from the b-block
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        conn[0] = bnodes[2 * i + 2 + (2 * nx + 1) * (2 * j)];
        conn[1] = bnodes[2 * i + 2 + (2 * nx + 1) * (2 * j + 2)];
        conn[2] = bnodes[2 * i + (2 * nx + 1) * (2 * j)];
        conn[3] = bnodes[2 * i + 2 + (2 * nx + 1) * (2 * j + 1)];
        conn[4] = bnodes[2 * i + 1 + (2 * nx + 1) * (2 * j + 1)];
        conn[5] = bnodes[2 * i + 1 + (2 * nx + 1) * (2 * j)];
        conn += 6;
        elem_node_ptr[n + 1] = elem_node_ptr[n] + 6;
        n++;

        conn[0] = bnodes[2 * i + (2 * nx + 1) * (2 * j + 2)];
        conn[1] = bnodes[2 * i + (2 * nx + 1) * (2 * j)];
        conn[2] = bnodes[2 * i + 2 + (2 * nx + 1) * (2 * j + 2)];
        conn[3] = bnodes[2 * i + (2 * nx + 1) * (2 * j + 1)];
        conn[4] = bnodes[2 * i + 1 + (2 * nx + 1) * (2 * j + 1)];
        conn[5] = bnodes[2 * i + 1 + (2 * nx + 1) * (2 * j + 2)];
        conn += 6;
        elem_node_ptr[n + 1] = elem_node_ptr[n] + 6;
        n++;
      }
    }

    // Set the identity numbers
    int *elem_id_nums = new int[num_elements];
    memset(elem_id_nums, 0, num_elements * sizeof(int));

    // Set the boundary conditions
    int num_bcs = 4 * ny + 1;
    int *bc_nodes = new int[num_bcs];
    for (int j = 0; j < 4 * ny + 1; j++) {
      bc_nodes[j] = anodes[j * (4 * nx + 1)];
    }

    // Set the connectivity
    creator->setGlobalConnectivity(num_nodes, num_elements, elem_node_ptr,
                                   elem_node_conn, elem_id_nums);

    // Set the boundary conditions
    creator->setBoundaryConditions(num_bcs, bc_nodes);

    // Set the nodal locations
    creator->setNodes(Xpts);

    // Set the dependent nodes
    creator->setDependentNodes(dep_nodes, dep_ptr, dep_conn, dep_weights);

    // Free all the allocated data
    delete[] Xpts;
    delete[] elem_node_ptr;
    delete[] elem_node_conn;
    delete[] elem_id_nums;
    delete[] bc_nodes;
    delete[] dep_ptr;
    delete[] dep_conn;
    delete[] dep_weights;
    delete[] anodes;
    delete[] bnodes;
  }

  // This call must occur on all processor
  TACSElement *element = elem;
  creator->setElements(1, &element);

  // Set the reordering type
  creator->setReorderingType(TACSAssembler::TACS_AMD_ORDER,
                             TACSAssembler::APPROXIMATE_SCHUR);

  // Create the TACSAssembler object
  TACSAssembler *assembler = creator->createTACS();
  assembler->incref();

  // Create the preconditioner
  TACSBVec *res = assembler->createVec();
  TACSBVec *ans = assembler->createVec();

  // Increment the reference count to the matrix/vectors
  res->incref();
  ans->incref();

  // Set the parameters for the stiffness/Jacobian matrix
  double alpha = 1.0, beta = 0.0, gamma = 0.0;

  // Allocate, assemble and factor the TACSSchurMat and preconditioner
  TACSSchurMat *schur_mat = assembler->createSchurMat();
  schur_mat->incref();

  int lev = 4500;
  double fill = 10.0;
  int reorder_schur = 1;
  TACSSchurPc *schur_pc = new TACSSchurPc(schur_mat, lev, fill, reorder_schur);
  schur_pc->incref();

  double schur_time = MPI_Wtime();
  assembler->assembleJacobian(alpha, beta, gamma, res, schur_mat);
  schur_pc->factor();
  schur_time = MPI_Wtime() - schur_time;

  if (rank == 0) {
    printf("TACSSchurMat/TACSSchurPc assembly and factorization time: %15.5e\n",
           schur_time);
  }

  // Allocate, assemble and factor the TACSParallelMat using the
  // block-cyclic preconditioner
  TACSParallelMat *par_mat = assembler->createMat();
  par_mat->incref();

  TACSBlockCyclicPc *par_pc = new TACSBlockCyclicPc(par_mat, 4, 1);
  par_pc->incref();

  double par_time = MPI_Wtime();
  assembler->assembleJacobian(alpha, beta, gamma, res, par_mat);
  par_pc->factor();
  par_time = MPI_Wtime() - par_time;

  if (rank == 0) {
    printf(
        "TACSParallelMat/TACSBlockCyclicPc assembly and factorization time: "
        "%15.5e\n",
        par_time);
  }

  // Allocate the GMRES object with the TACSSchurMat matrix
  int gmres_iters = 10;  // Number of GMRES iterations
  int nrestart = 2;      // Number of allowed restarts
  int is_flexible = 1;   // Is a flexible preconditioner?
  GMRES *schur_gmres =
      new GMRES(schur_mat, schur_pc, gmres_iters, nrestart, is_flexible);

  res->set(1.0);
  assembler->applyBCs(res);
  schur_gmres->solve(res, ans);
  assembler->setVariables(ans);

  // Create an TACSToFH5 object for writing output to files
  ElementType etype = TACS_PLANE_STRESS_ELEMENT;
  int write_flag = (TACS_OUTPUT_CONNECTIVITY | TACS_OUTPUT_NODES |
                    TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                    TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(assembler, etype, write_flag);
  f5->incref();
  f5->writeToFile("triangle.f5");
  f5->decref();

  // Decrease the reference count to the linear algebra objects
  schur_pc->decref();
  schur_mat->decref();
  par_pc->decref();
  par_mat->decref();
  ans->decref();
  res->decref();

  // Decrease the reference count to everything else
  stiff->decref();
  elem->decref();
  assembler->decref();
  creator->decref();

  MPI_Finalize();
  return (0);
}
