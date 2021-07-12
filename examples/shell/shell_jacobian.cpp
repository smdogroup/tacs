// #include "TACSElementAlgebra.h"
// #include "TACSElementVerification.h"
// #include "TACSShellElementBasis.h"
// #include "TACSShellElement.h"
// #include "TACSElementVerification.h"
// #include "TACSConstitutiveVerification.h"
// #include "TACSElementAlgebra.h"
// #include "TACSIsoShellConstitutive.h"

#include "TACSDirector.h"

// typedef TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
//                          TACSLinearizedRotation, TACSShellLinearModel> TACSQuadLinearShell;

// typedef TACSShellElement<TACSQuadQuadraticQuadrature, TACSShellQuadBasis<3>,
//                          TACSLinearizedRotation, TACSShellLinearModel> TACSQuadQuadraticShell;

// typedef TACSShellElement<TACSTriQuadraticQuadrature, TACSShellTriQuadraticBasis,
//     TACSLinearizedRotation, TACSShellLinearModel> TACSTriQuadraticShell;

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);

  // Get the rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  // TacsScalar rho = 2700.0;
  // TacsScalar specific_heat = 921.096;
  // TacsScalar E = 70e3;
  // TacsScalar nu = 0.3;
  // TacsScalar ys = 270.0;
  // TacsScalar cte = 24.0e-6;
  // TacsScalar kappa = 230.0;
  // TACSMaterialProperties *props =
  //   new TACSMaterialProperties(rho, specific_heat, E, nu, ys, cte, kappa);

  // TacsScalar axis[] = {0.0, 1.0, 0.0};
  // TACSShellTransform *transform = new TACSShellRefAxisTransform(axis);

  // TacsScalar t = 0.01;
  // int t_num = 0;
  // TACSShellConstitutive *con = new TACSIsoShellConstitutive(props, t, t_num);

  // TACSElement *linear_shell = new TACSQuadLinearShell(transform, con);
  // linear_shell->incref();

  // TACSElement *quadratic_shell = new TACSQuadQuadraticShell(transform, con);
  // quadratic_shell->incref();

  // const int OFFSET = 1;
  // const int VARS_PER_NODE = OFFSET + TACSLinearizedRotation::NUM_PARAMETERS;
  const int NUM_NODES = 3;

  // const int NUM_VARS = VARS_PER_NODE*NUM_NODES;
  // int elemIndex = 0;
  // double time = 0.0;
  // TacsScalar Xpts[3*NUM_NODES];
  // TacsScalar vars[NUM_VARS], dvars[NUM_VARS], ddvars[NUM_VARS];

  // // Set the values of the
  // TacsGenerateRandomArray(Xpts, 3*NUM_NODES);
  // TacsGenerateRandomArray(vars, 6*NUM_NODES);
  // TacsGenerateRandomArray(dvars, 6*NUM_NODES);
  // TacsGenerateRandomArray(ddvars, 6*NUM_NODES);

  // TacsTestElementResidual(linear_shell, elemIndex, time, Xpts, vars, dvars, ddvars);
  // TacsTestElementResidual(quadratic_shell, elemIndex, time, Xpts, vars, dvars, ddvars);

  // TacsTestElementJacobian(linear_shell, elemIndex, time, Xpts, vars, dvars, ddvars);
  // TacsTestElementJacobian(quadratic_shell, elemIndex, time, Xpts, vars, dvars, ddvars);

  static const int OFFSET = 1;
  static const int VARS_PER_NODE = OFFSET + TACSQuaternionRotation::NUM_PARAMETERS;

  // double dh = 1e-30;
  // TacsTestDirector<VARS_PER_NODE, OFFSET, NUM_NODES, TACSLinearizedRotation>(dh);
  // TacsTestDirector<VARS_PER_NODE, OFFSET, NUM_NODES, TACSQuadraticRotation>(dh);
  // TacsTestDirector<VARS_PER_NODE, OFFSET, NUM_NODES, TACSQuaternionRotation>(dh);

  double dh_res = 1e-5;
  TacsTestDirectorResidual<VARS_PER_NODE, OFFSET, NUM_NODES, TACSLinearizedRotation>(dh_res);
  TacsTestDirectorResidual<VARS_PER_NODE, OFFSET, NUM_NODES, TACSQuadraticRotation>(dh_res);
  TacsTestDirectorResidual<VARS_PER_NODE, OFFSET, NUM_NODES, TACSQuaternionRotation>(dh_res);


  // TacsTestShellTyingStrain<6, TACSShellQuadLinearBasis, TACSShellLinearModel>();
  // TacsTestShellTyingStrain<6, TACSShellQuadLinearBasis, TACSShellNonlinearModel>();
  // TacsTestShellTyingStrain<6, TACSShellQuadQuadraticBasis, TACSShellLinearModel>();

  // TacsTestShellModelDerivatives<6, TACSShellQuadLinearBasis, TACSShellNonlinearModel>();

  // TacsTestShellUtilities<4, TACSShellQuadQuadraticBasis>(1e-30);
  // TacsTestShellUtilities<4, TACSShellQuadBasis<2>>(1e-6);
  // TacsTestShellUtilities<4, TACSShellQuadBasis<3>>(1e-6);

  // TacsScalar alpha = 1.0, beta = 0.0, gamma = 0.0;
  // double t0;

  // t0 = MPI_Wtime();
  // for ( int i = 0; i < 4*500; i++ ){
  //   linear_shell->addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);
  // }
  // t0 = MPI_Wtime() - t0;
  // printf("2nd order residual Time = %15.10e\n", t0);

  // t0 = MPI_Wtime();
  // for ( int i = 0; i < 4*500; i++ ){
  //   linear_shell->addJacobian(elemIndex, time, alpha, beta, gamma,
  //                             Xpts, vars, dvars, ddvars, res, mat);
  // }
  // t0 = MPI_Wtime() - t0;
  // printf("2nd order jacobian Time = %15.10e\n", t0);

  // t0 = MPI_Wtime();
  // for ( int i = 0; i < 500; i++ ){
  //   quadratic_shell->addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);
  // }
  // t0 = MPI_Wtime() - t0;
  // printf("3rd order residual Time = %15.10e\n", t0);

  // t0 = MPI_Wtime();
  // for ( int i = 0; i < 500; i++ ){
  //   quadratic_shell->addJacobian(elemIndex, time, alpha, beta, gamma,
  //                                Xpts, vars, dvars, ddvars, res, mat);
  // }
  // t0 = MPI_Wtime() - t0;
  // printf("3rd order jacobian Time = %15.10e\n", t0);

  // linear_shell->decref();
  // quadratic_shell->decref();

  MPI_Finalize();
  return 0;
}
