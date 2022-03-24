#include "TACSElementVerification.h"
#include "TACSBasicBeamConstitutive.h"
#include "TACSBeamElement.h"
#include "TACSShellElementDefs.h"

typedef TACSBeamElement<TACSBeamQuadraticQuadrature, TACSBeamBasis<3>,
                        TACSLinearizedRotation, TACSBeamLinearModel> TACSQuadBeam;

typedef TACSBeamElement<TACSBeamQuadraticQuadrature, TACSBeamBasis<3>,
                        TACSQuadraticRotation, TACSBeamLinearModel> TACSQuadBeamModRot;

typedef TACSBeamElement<TACSBeamQuadraticQuadrature, TACSBeamBasis<3>,
                        TACSQuaternionRotation, TACSBeamLinearModel> TACSQuadBeamQuaternion;

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);

  // Get the rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  TacsTestBeamModelDerivatives<6, TACSBeamBasis<3>, TACSBeamLinearModel>();

  TacsScalar axis[] = {0.1, 0.43, 0.5};
  TACSBeamRefAxisTransform *transform = new TACSBeamRefAxisTransform(axis);
  transform->incref();

  TACSBasicBeamConstitutive *stiff =
    new TACSBasicBeamConstitutive(1.0, 1.0, 1.0, 1.0,
                                  0.23, 0.3251, 1.43,
                                  0.47, 0.71, 1.93);


  // TACSElement *beam = new TACSQuadBeam(transform, stiff);
  // TACSElement *beam = new TACSQuadBeamModRot(transform, stiff);
  TACSElement *beam = new TACSQuadBeamQuaternion(transform, stiff);
  beam->incref();

  int vars_per_node = beam->getVarsPerNode();
  int num_nodes = beam->getNumNodes();
  int num_vars = num_nodes*vars_per_node;
  double time = 0.0;
  int elemIndex = 0;

  TacsScalar *Xpts = new TacsScalar[ 3*num_nodes ];
  TacsScalar *vars = new TacsScalar[ num_vars ];
  TacsScalar *dvars = new TacsScalar[ num_vars ];
  TacsScalar *ddvars = new TacsScalar[ num_vars ];
  TacsGenerateRandomArray(Xpts, 3*num_nodes);
  TacsGenerateRandomArray(vars, num_vars);
  TacsGenerateRandomArray(dvars, num_vars);
  TacsGenerateRandomArray(ddvars, num_vars);

  // Zero out the multipliers so the residual test passes
  if (vars_per_node == 8){
    for ( int i = 0; i < num_nodes; i++ ){
      vars[vars_per_node*i + 7] = 0.0;
    }
  }

  // Check the residual formulation against Lagrange's equations. Not all elements
  // will pass this test - for instance the thermal beam elements.
  TacsTestElementResidual(beam, elemIndex, time, Xpts, vars, dvars, ddvars);

  // Test the quantity derivative implementation
  TacsTestElementQuantityDVSens(beam, elemIndex, TACS_FAILURE_INDEX,
                                time, Xpts, vars, dvars, ddvars);
  TacsTestElementQuantitySVSens(beam, elemIndex, TACS_FAILURE_INDEX,
                                time, Xpts, vars, dvars, ddvars);
  TacsTestElementQuantityXptSens(beam, elemIndex, TACS_FAILURE_INDEX,
                                 time, Xpts, vars, dvars, ddvars);

  int dvLen = 3;
  TacsScalar x[3];
  TacsTestAdjResProduct(beam, elemIndex, time, Xpts, vars, dvars, ddvars, dvLen, x);

  TacsTestAdjResXptProduct(beam, elemIndex, time, Xpts, vars, dvars, ddvars);

  beam->decref();

  MPI_Finalize();
  return 0;
}
