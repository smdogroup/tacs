#include "TACSElementVerification.h"

// Include the models
#include "TACSHeatConduction.h"
#include "TACSLinearElasticity.h"
#include "TACSNeohookean.h"
#include "TACSThermoelasticity.h"

// Include the constitutive classes
#include "TACSIsoShellConstitutive.h"
#include "TACSPlaneStressConstitutive.h"
#include "TACSSolidConstitutive.h"

// Include the basis functions
#include "TACSHexaBasis.h"
#include "TACSQuadBasis.h"
#include "TACSTetrahedralBasis.h"
#include "TACSTriangularBasis.h"

// Include the element classes
#include "TACSElement2D.h"
#include "TACSElement3D.h"

/*
  Apply all the tests to the element
*/
void test_element(TACSElement *element, int elemIndex, double time,
                  const TacsScalar Xpts[], const TacsScalar vars[],
                  const TacsScalar dvars[], const TacsScalar ddvars[]) {
  // Get the design variable numbers
  int dvLen = element->getDesignVarNums(elemIndex, 0, NULL);
  TacsScalar *x = new TacsScalar[dvLen];
  element->getDesignVars(elemIndex, dvLen, x);

  // Test the element
#ifdef TACS_USE_COMPLEX
  double dh = 1e-30;
#else
  double dh = 1e-5;
#endif
  // TacsTestElementResidual(element, elemIndex, time, Xpts, vars, dvars,
  // ddvars, dh);
  TacsTestElementJacobian(element, elemIndex, time, Xpts, vars, dvars, ddvars,
                          -1, dh);
  TacsTestAdjResProduct(element, elemIndex, time, Xpts, vars, dvars, ddvars,
                        dvLen, x, dh);
  TacsTestAdjResXptProduct(element, elemIndex, time, Xpts, vars, dvars, ddvars,
                           dh);
  TacsTestElementMatDVSens(element, TACS_MASS_MATRIX, elemIndex, time, Xpts,
                           vars, dvLen, x, dh);
  TacsTestElementMatDVSens(element, TACS_STIFFNESS_MATRIX, elemIndex, time,
                           Xpts, vars, dvLen, x, dh);
  TacsTestElementMatDVSens(element, TACS_GEOMETRIC_STIFFNESS_MATRIX, elemIndex,
                           time, Xpts, vars, dvLen, x, dh);
  TacsTestElementMatSVSens(element, TACS_GEOMETRIC_STIFFNESS_MATRIX, elemIndex,
                           time, Xpts, vars, dh);

  TACSElementModel *model = element->getElementModel();
  if (model) {
    TacsTestElementModel(model, elemIndex, time, dh);
  }

  delete[] x;
}

/*
  The following code tests the element implementation to see if the
  computation of the residual is consistent with the energy
  formulation, to check if the Jacobian matrix is consistent with the
  residual and to test certain design-dependent code.

  Useage:
  ./profile_elements [fd=value]
*/
int main(int argc, char *argv[]) {
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Test only the input element type - if specified, otherwise
  // test everything
  const char *ename = NULL;
  if (argc > 1) {
    ename = argv[1];
  }

  const int MAX_NODES = 64;
  const int MAX_VARS_PER_NODE = 8;
  const int MAX_VARS = MAX_NODES * MAX_VARS_PER_NODE;

  // Set the// Test the element
#ifdef TACS_USE_COMPLEX
  double dh = 1e-30;
#else
  double dh = 1e-5;
#endif

  // Set the simulation time
  int elemIndex = 0;
  double time = 0.0;

  // Set the variable arrays
  TacsScalar Xpts[3 * MAX_NODES];
  TacsScalar vars[MAX_VARS], dvars[MAX_VARS], ddvars[MAX_VARS];

  // Generate random arrays
  TacsGenerateRandomArray(Xpts, 3 * MAX_NODES, 0.0, 1.0);
  TacsGenerateRandomArray(vars, MAX_VARS);
  TacsGenerateRandomArray(dvars, MAX_VARS);
  TacsGenerateRandomArray(ddvars, MAX_VARS);

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
  props->incref();

  // Create the basis functions for 3D
  const int NUM_3D_BASIS = 5;
  TACSElementBasis *basis3d[NUM_3D_BASIS];
  basis3d[0] = new TACSLinearTetrahedralBasis();
  basis3d[1] = new TACSQuadraticTetrahedralBasis();
  basis3d[2] = new TACSLinearHexaBasis();
  basis3d[3] = new TACSQuadraticHexaBasis();
  basis3d[4] = new TACSCubicHexaBasis();
  for (int i = 0; i < NUM_3D_BASIS; i++) {
    basis3d[i]->incref();
  }

  // Create the basis functions for 2D
  const int NUM_2D_BASIS = 5;
  TACSElementBasis *basis2d[NUM_3D_BASIS];
  basis2d[0] = new TACSLinearTriangleBasis();
  basis2d[1] = new TACSQuadraticTriangleBasis();
  basis2d[2] = new TACSLinearQuadBasis();
  basis2d[3] = new TACSQuadraticQuadBasis();
  basis2d[4] = new TACSCubicQuadBasis();
  for (int i = 0; i < NUM_2D_BASIS; i++) {
    basis2d[i]->incref();
  }

  // Create stiffness (need class)
  TACSSolidConstitutive *con3d = new TACSSolidConstitutive(props, 1.0, 0);
  con3d->incref();

  TACSPlaneStressConstitutive *con2d =
      new TACSPlaneStressConstitutive(props, 1.0, 0);
  con2d->incref();

  TACSShellConstitutive *conShell = new TACSIsoShellConstitutive(props, 1.0, 0);
  conShell->incref();

  // Set the model type
  const int NUM_3D_MODELS = 5;
  TACSElementModel *model3d[NUM_3D_MODELS];
  model3d[0] = new TACSHeatConduction3D(con3d);
  model3d[1] = new TACSLinearElasticity3D(con3d, TACS_LINEAR_STRAIN);
  model3d[2] = new TACSLinearElasticity3D(con3d, TACS_NONLINEAR_STRAIN);
  model3d[3] = new TACSLinearThermoelasticity3D(con3d, TACS_LINEAR_STRAIN);
  model3d[4] = new TACSNeohookean3D(2.34, 5.73);
  for (int i = 0; i < NUM_3D_MODELS; i++) {
    model3d[i]->incref();
  }

  const int NUM_2D_MODELS = 4;
  TACSElementModel *model2d[NUM_2D_MODELS];
  model2d[0] = new TACSHeatConduction2D(con2d);
  model2d[1] = new TACSLinearElasticity2D(con2d, TACS_LINEAR_STRAIN);
  model2d[2] = new TACSLinearElasticity2D(con2d, TACS_NONLINEAR_STRAIN);
  model2d[3] = new TACSLinearThermoelasticity2D(con2d, TACS_LINEAR_STRAIN);
  for (int i = 0; i < NUM_2D_MODELS; i++) {
    model2d[i]->incref();
  }

  for (int j = 0; j < NUM_3D_MODELS; j++) {
    TacsTestElementModel(model3d[j], elemIndex, time, dh);
    for (int i = NUM_3D_BASIS - 1; i < NUM_3D_BASIS; i++) {
      printf("Testing with model %s with basis functions %s\n",
             model3d[j]->getObjectName(), basis3d[i]->getObjectName());
      TACSElement *element = new TACSElement3D(model3d[j], basis3d[i]);
      element->incref();
      test_element(element, elemIndex, time, Xpts, vars, dvars, ddvars);
      element->decref();
    }
  }

  for (int j = 0; j < NUM_2D_MODELS; j++) {
    TacsTestElementModel(model2d[j], elemIndex, time, dh);
    for (int i = 0; i < NUM_2D_BASIS; i++) {
      printf("Testing with model %s with basis functions %s\n",
             model2d[j]->getObjectName(), basis2d[i]->getObjectName());
      TACSElement *element = new TACSElement2D(model2d[j], basis2d[i]);
      element->incref();
      test_element(element, elemIndex, time, Xpts, vars, dvars, ddvars);
      element->decref();
    }
  }

  MPI_Finalize();
  return (0);
}
