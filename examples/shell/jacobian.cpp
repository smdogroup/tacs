#include "TACSConstitutiveVerification.h"
#include "TACSElementVerification.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSShellElementDefs.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // Get the rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  int no_drill = 0, no_dynamics = 0, test_models = 0, test_constitutive = 0;
  for (int k = 0; k < argc; k++) {
    if (strcmp(argv[k], "no_drill") == 0) {
      no_drill = 1;
    }
    if (strcmp(argv[k], "no_dynamics") == 0) {
      no_drill = 1;
    }
    if (strcmp(argv[k], "test_models") == 0) {
      test_models = 1;
    }
    if (strcmp(argv[k], "test_constitutive") == 0) {
      test_constitutive = 1;
    }
  }

  int elemIndex = 0;
  TacsScalar rho = 2700.0;
  TacsScalar specific_heat = 921.096;
  TacsScalar E = 70e3;
  TacsScalar nu = 0.3;
  TacsScalar ys = 270.0;
  TacsScalar cte = 24.0e-6;
  TacsScalar kappa = 230.0;

  // Set the density to zero (to remove dynamics)
  if (no_dynamics) {
    printf("Setting the material density to zero\n");
    rho = 0.0;
  }

  TACSMaterialProperties *props =
      new TACSMaterialProperties(rho, specific_heat, E, nu, ys, cte, kappa);

  TacsScalar axis[] = {0.0, 1.0, 1.0};
  TACSShellTransform *transform = new TACSShellRefAxisTransform(axis);

  TacsScalar t = 0.01;
  int t_num = 0;
  TACSShellConstitutive *con = new TACSIsoShellConstitutive(props, t, t_num);

  if (test_constitutive) {
    TacsTestConstitutive(con, elemIndex);
  }

  // Set the drilling regularization to zero
  if (no_drill) {
    printf("Setting the drilling regularization to zero\n");
    con->setDrillingRegularization(0.0);
  }

  // Test the different shell element model expressions
  if (test_models) {
    TacsTestShellModelDerivatives<6, TACSShellQuadBasis<2>,
                                  TACSShellLinearModel>();
    TacsTestShellModelDerivatives<6, TACSShellQuadBasis<2>,
                                  TACSShellNonlinearModel>();
    TacsTestShellModelDerivatives<6, TACSShellQuadBasis<2>,
                                  TACSShellInplaneLinearModel>();
    TacsTestShellModelDerivatives<6, TACSShellQuadBasis<2>,
                                  TACSShellInplaneNonlinearModel>();
  }

  TACSElement *shell = NULL;
  if (argc > 1) {
    for (int k = 0; k < argc; k++) {
      shell = TacsCreateShellByName(argv[k], transform, con);
      if (shell) {
        break;
      }
    }
  }
  if (!shell) {
    shell = new TACSQuad4Shell(transform, con);
  }
  shell->incref();

  int vars_per_node = shell->getVarsPerNode();
  int num_nodes = shell->getNumNodes();
  int num_vars = num_nodes * vars_per_node;
  double time = 0.0;

  // Start and end column to test in the Jacobian matrix
  int start = 0, end = 0;

  if (argc > 1) {
    int temp = 0;
    for (int k = 0; k < argc; k++) {
      if (sscanf(argv[k], "%d", &temp) == 1) {
        if (temp >= 0 && temp <= num_vars) {
          start = 0;
          end = temp;
        }
      }
    }
  }

  TacsScalar *Xpts = new TacsScalar[3 * num_nodes];
  TacsScalar *vars = new TacsScalar[num_vars];
  TacsScalar *dvars = new TacsScalar[num_vars];
  TacsScalar *ddvars = new TacsScalar[num_vars];
  TacsGenerateRandomArray(Xpts, 3 * num_nodes);
  TacsGenerateRandomArray(vars, num_vars);
  TacsGenerateRandomArray(dvars, num_vars);
  TacsGenerateRandomArray(ddvars, num_vars);

  // Zero out the multipliers so the residual test passes
  if (vars_per_node == 8) {
    for (int i = 0; i < num_nodes; i++) {
      vars[vars_per_node * i + 7] = 0.0;
    }
  }

  // Check the residual formulation against Lagrange's equations. Not all
  // elements will pass this test - for instance the thermal shell elements.
  TacsTestElementResidual(shell, elemIndex, time, Xpts, vars, dvars, ddvars);

  // Check the implementation of the element Jacobian matrix against the
  // implementation of the residual
  TacsTestElementJacobian(shell, elemIndex, time, Xpts, vars, dvars, ddvars);

  // Check specific columns of the Jacobian matrix (if specified on the command
  // line)
  for (int col = start; col < end; col++) {
    TacsTestElementJacobian(shell, elemIndex, time, Xpts, vars, dvars, ddvars,
                            col);
  }

  // Test the implementation of the derivative of the adjoint-residual product
  const int dvLen = 1;
  TacsScalar x[dvLen];
  shell->getDesignVars(elemIndex, dvLen, x);
  TacsTestAdjResProduct(shell, elemIndex, time, Xpts, vars, dvars, ddvars,
                        dvLen, x);

  // Test the implementation of the derivative of the adjoint-residual product
  // w.r.t. nodes
  TacsTestAdjResXptProduct(shell, elemIndex, time, Xpts, vars, dvars, ddvars);

  int quantityType = TACS_FAILURE_INDEX;
  TacsTestElementQuantityDVSens(shell, elemIndex, quantityType, time, Xpts,
                                vars, dvars, ddvars);
  TacsTestElementQuantitySVSens(shell, elemIndex, quantityType, time, Xpts,
                                vars, dvars, ddvars);

  delete[] Xpts;
  delete[] vars;
  delete[] dvars;
  delete[] ddvars;
  shell->decref();

  MPI_Finalize();
  return 0;
}
