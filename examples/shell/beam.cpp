#include "TACSBeamElement.h"
#include "TACSConstitutiveVerification.h"
#include "TACSElementVerification.h"
#include "TACSIsoTubeBeamConstitutive.h"
#include "TACSShellElementDefs.h"

typedef TACSBeamElement<TACSBeamQuadraticQuadrature, TACSBeamBasis<3>,
                        TACSLinearizedRotation, TACSBeamLinearModel>
    TACSQuadBeam;

typedef TACSBeamElement<TACSBeamQuadraticQuadrature, TACSBeamBasis<3>,
                        TACSQuadraticRotation, TACSBeamLinearModel>
    TACSQuadBeamModRot;

typedef TACSBeamElement<TACSBeamQuadraticQuadrature, TACSBeamBasis<3>,
                        TACSQuaternionRotation, TACSBeamLinearModel>
    TACSQuadBeamQuaternion;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // Get the rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  TacsTestBeamModelDerivatives<6, TACSBeamBasis<3>, TACSBeamLinearModel>();

  TacsScalar axis[] = {0.0, 1.0, 0.0};
  TACSBeamRefAxisTransform *transform = new TACSBeamRefAxisTransform(axis);
  transform->incref();

  TacsScalar rho = 0.0;  // 2700.0;
  TacsScalar specific_heat = 921.096;
  TacsScalar E = 70e3;
  TacsScalar nu = 0.3;
  TacsScalar ys = 270.0;
  TacsScalar cte = 24.0e-6;
  TacsScalar kappa = 230.0;

  TACSMaterialProperties *props =
      new TACSMaterialProperties(rho, specific_heat, E, nu, ys, cte, kappa);

  TacsScalar inner = 0.12;
  TacsScalar wall = 0.05;
  int inner_dv = 0, wall_dv = 1;
  TacsScalar inner_lb = 0.01, wall_lb = 0.01;
  TacsScalar inner_ub = 0.5, wall_ub = 0.5;

  TACSIsoTubeBeamConstitutive *stiff =
      new TACSIsoTubeBeamConstitutive(props, inner, wall, inner_dv, wall_dv,
                                      inner_lb, inner_ub, wall_lb, wall_ub);
  TacsTestConstitutive(stiff, 0);

  TACSElement *beam = new TACSQuadBeam(transform, stiff);
  // TACSElement *beam = new TACSQuadBeamModRot(transform, stiff);
  // TACSElement *beam = new TACSQuadBeamQuaternion(transform, stiff);
  beam->incref();

  int vars_per_node = beam->getVarsPerNode();
  int num_nodes = beam->getNumNodes();
  int num_vars = num_nodes * vars_per_node;
  double time = 0.0;
  int elemIndex = 0;

  // Set the state variables
  TacsScalar *vars = new TacsScalar[num_vars];
  TacsScalar *dvars = new TacsScalar[num_vars];
  TacsScalar *ddvars = new TacsScalar[num_vars];
  TacsGenerateRandomArray(vars, num_vars);
  TacsGenerateRandomArray(dvars, num_vars);
  TacsGenerateRandomArray(ddvars, num_vars);

  // Set the node locations
  TacsScalar *Xpts = new TacsScalar[3 * num_nodes];
  TacsGenerateRandomArray(Xpts, 3 * num_nodes);

  // Zero out the multipliers so the residual test passes
  if (vars_per_node == 8) {
    for (int i = 0; i < num_nodes; i++) {
      vars[vars_per_node * i + 7] = 0.0;
    }
  }

  // Check the residual formulation against Lagrange's equations. Not all
  // elements will pass this test - for instance the thermal beam elements.
  TacsTestElementResidual(beam, elemIndex, time, Xpts, vars, dvars, ddvars);

  // Test the quantity derivative implementation
  TacsTestElementQuantityDVSens(beam, elemIndex, TACS_FAILURE_INDEX, time, Xpts,
                                vars, dvars, ddvars);
  TacsTestElementQuantitySVSens(beam, elemIndex, TACS_FAILURE_INDEX, time, Xpts,
                                vars, dvars, ddvars);
  TacsTestElementQuantityXptSens(beam, elemIndex, TACS_FAILURE_INDEX, time,
                                 Xpts, vars, dvars, ddvars);

  TacsTestAdjResXptProduct(beam, elemIndex, time, Xpts, vars, dvars, ddvars);

  const int dvLen = 3;
  int dvNums[dvLen];
  TacsScalar x[dvLen];
  TacsScalar xelem[dvLen];

  // Get the design variables from the element
  int ndvs = beam->getDesignVarNums(elemIndex, dvLen, dvNums);
  beam->getDesignVars(elemIndex, dvLen, xelem);

  // Set the design variable numbers into the design variable vector
  for (int i = 0; i < ndvs; i++) {
    x[dvNums[i]] = xelem[i];
  }

  TacsTestAdjResProduct(beam, elemIndex, time, Xpts, vars, dvars, ddvars, dvLen,
                        x);

  beam->decref();

  MPI_Finalize();
  return 0;
}
