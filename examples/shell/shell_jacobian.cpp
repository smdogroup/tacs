#include "TACSElementVerification.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSShellElementDefs.h"

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);

  // Get the rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  TacsScalar rho = 2700.0;
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

  TACSElement *shell = NULL;
  if (argc > 1){
    printf("argv[1] = %s", argv[1]);
    if (strcmp(argv[1], "TACSQuad2Shell") == 0){
      shell = new TACSQuad2Shell(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad3Shell") == 0){
      shell = new TACSQuad3Shell(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad4Shell") == 0){
      shell = new TACSQuad4Shell(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad2ThermalShell") == 0){
      shell = new TACSQuad2ThermalShell(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad3ThermalShell") == 0){
      shell = new TACSQuad3ThermalShell(transform, con);
    }
    else if (strcmp(argv[1], "TACSQuad4ThermalShell") == 0){
      shell = new TACSQuad4ThermalShell(transform, con);
    }
    else {
      shell = new TACSQuad2Shell(transform, con);
    }
  }
  else {
    shell = new TACSQuad2Shell(transform, con);
  }
  shell->incref();

  int vars_per_node = shell->getVarsPerNode();
  int num_nodes = shell->getNumNodes();
  int num_vars = num_nodes*vars_per_node;
  int elemIndex = 0;
  double time = 0.0;

  TacsScalar *Xpts = new TacsScalar[ 3*num_nodes ];
  TacsScalar *vars = new TacsScalar[ num_vars ];
  TacsScalar *dvars = new TacsScalar[ num_vars ];
  TacsScalar *ddvars = new TacsScalar[ num_vars ];
  TacsGenerateRandomArray(Xpts, 3*num_nodes);
  TacsGenerateRandomArray(vars, num_vars);
  TacsGenerateRandomArray(dvars, num_vars);
  TacsGenerateRandomArray(ddvars, num_vars);

  TacsTestElementResidual(shell, elemIndex, time, Xpts, vars, dvars, ddvars);
  TacsTestElementJacobian(shell, elemIndex, time, Xpts, vars, dvars, ddvars);

  delete [] Xpts;
  delete [] vars;
  delete [] dvars;
  delete [] ddvars;
  shell->decref();

  MPI_Finalize();
  return 0;
}
