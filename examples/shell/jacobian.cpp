#include "TACSElementVerification.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSShellElementDefs.h"

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);

  // Get the rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  TacsScalar rho = 0.0; // 2700.0;
  TacsScalar specific_heat = 921.096;
  TacsScalar E = 70e3;
  TacsScalar nu = 0.3;
  TacsScalar ys = 270.0;
  TacsScalar cte = 24.0e-6;
  TacsScalar kappa = 230.0;
  TACSMaterialProperties *props =
    new TACSMaterialProperties(rho, specific_heat, E, nu, ys, cte, kappa);

  TacsScalar axis[] = {0.0, 1.0, 1.0};
  TACSShellTransform *transform = new TACSShellRefAxisTransform(axis);

  TacsScalar t = 0.01;
  int t_num = 0;
  TACSShellConstitutive *con = new TACSIsoShellConstitutive(props, t, t_num);

  // Start and end column to test in the Jacobian matrix
  int start = 0, end = 0;

  TACSElement *shell = NULL;
  if (argc > 1){
    if (strcmp(argv[1], "TACSQuad4ShellModRot") == 0){
      shell = new TACSQuad4ShellModRot(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad9ShellModRot") == 0){
      shell = new TACSQuad9ShellModRot(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad16ShellModRot") == 0){
      shell = new TACSQuad16ShellModRot(transform, con);
    }
    else if (strcmp(argv[1], "TACSTri3ShellModRot") == 0){
      shell = new TACSTri3ShellModRot(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad4Shell") == 0){
      shell = new TACSQuad4Shell(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad9Shell") == 0){
      shell = new TACSQuad9Shell(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad16Shell") == 0){
      shell = new TACSQuad16Shell(transform, con);
    }
    else if (strcmp(argv[1], "TACSTri3Shell") == 0){
      shell = new TACSTri3Shell(transform, con);
    }
    // else if (strcmp(argv[1], "TACSTri6Shell") == 0){
    //   shell = new TACSTri6Shell(transform, con);
    // }
    else if (strcmp(argv[1], "TACSQuad4ThermalShell") == 0){
      shell = new TACSQuad4ThermalShell(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad9ThermalShell") == 0){
      shell = new TACSQuad9ThermalShell(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad16ThermalShell") == 0){
      shell = new TACSQuad16ThermalShell(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad4NonlinearShell") == 0){
      shell = new TACSQuad4NonlinearShell(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad9NonlinearShell") == 0){
      shell = new TACSQuad9NonlinearShell(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad16NonlinearShell") == 0){
      shell = new TACSQuad16NonlinearShell(transform, con);
    }
    else {
      shell = new TACSQuad4Shell(transform, con);
    }
  }
  else {
    shell = new TACSQuad4Shell(transform, con);
  }
  shell->incref();

  int vars_per_node = shell->getVarsPerNode();
  int num_nodes = shell->getNumNodes();
  int num_vars = num_nodes*vars_per_node;
  int elemIndex = 0;
  double time = 0.0;

  if (argc > 2){
    int temp = 0;
    if (sscanf(argv[2], "%d", &temp) == 1){
      if (temp >= 0 && temp <= num_vars){
        start = temp;
        end = num_vars;
      }
    }
  }
  if (argc > 3){
    int temp = 0;
    if (sscanf(argv[3], "%d", &temp) == 1){
      if (temp >= 0){
        end = temp;
        if (end > num_vars){
          end = num_vars;
        }
      }
    }
  }

  TacsScalar *Xpts = new TacsScalar[ 3*num_nodes ];
  TacsScalar *vars = new TacsScalar[ num_vars ];
  TacsScalar *dvars = new TacsScalar[ num_vars ];
  TacsScalar *ddvars = new TacsScalar[ num_vars ];
  TacsGenerateRandomArray(Xpts, 3*num_nodes);
  TacsGenerateRandomArray(vars, num_vars);
  TacsGenerateRandomArray(dvars, num_vars);
  TacsGenerateRandomArray(ddvars, num_vars);

  // Check the residual formulation against Lagrange's equations. Not all elements
  // will pass this test - for instance the thermal shell elements.
  TacsTestElementResidual(shell, elemIndex, time, Xpts, vars, dvars, ddvars);

  // Check the implementation of the element Jacobian matrix against the implementation
  // of the residual
  TacsTestElementJacobian(shell, elemIndex, time, Xpts, vars, dvars, ddvars);

  // Check specific columns of the Jacobian matrix (if specified on the command line)
  for ( int col = start; col < end; col++ ){
    TacsTestElementJacobian(shell, elemIndex, time, Xpts, vars, dvars, ddvars, col);
  }

  delete [] Xpts;
  delete [] vars;
  delete [] dvars;
  delete [] ddvars;
  shell->decref();

  MPI_Finalize();
  return 0;
}
