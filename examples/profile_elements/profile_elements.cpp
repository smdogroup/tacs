// Include the shell element classes
#include "isoFSDTStiffness.h"
#include "MITCShell.h"
#include "MITC9.h"

// Include the plane stress classes
#include "PlaneStressStiffness.h"
#include "PlaneStressQuad.h"
#include "PlaneStressTri6.h"

// Include the solid/3D element classes
#include "SolidStiffness.h"
#include "Solid.h"

// Include the multibody dynamics code
#include "RigidBody.h"
#include "KinematicConstraints.h"

/*
  Generate a random array of values
*/
static void generate_random_array( TacsScalar *array, int size, 
                                   TacsScalar lower=-1.0, 
                                   TacsScalar upper=1.0 ){
  for ( int i = 0; i < size; i++ ){
    array[i] = (upper - lower)*(rand()/((double)RAND_MAX+1)) + lower;
  }
}

/*
  Apply all the tests to the element
*/
void test_element( TACSElement *element,
                   double time,
                   const TacsScalar Xpts[], 
                   const TacsScalar vars[], 
                   const TacsScalar dvars[], 
                   const TacsScalar ddvars[],
                   int dvLen ){
  TacsScalar *x = new TacsScalar[ dvLen ];

  // Get the design variables from the element
  element->getDesignVars(x, dvLen);

  // Test the element
  TACSElement::setStepSize(1e-5);
#ifdef TACS_USE_COMPLEX
  TACSElement::setStepSize(1e-30);
#endif
  element->testResidual(time, Xpts, vars, dvars, ddvars);
  element->testJacobian(time, Xpts, vars, dvars, ddvars);
  element->testAdjResProduct(x, dvLen, time, Xpts, vars, dvars, ddvars);
  element->testAdjResXptProduct(time, Xpts, vars, dvars, ddvars);
  element->testStrainSVSens(Xpts, vars);
  element->testStrainXptSens(Xpts, vars);
  element->testJacobianXptSens(Xpts);
  element->testMatDVSensInnerProduct(STIFFNESS_MATRIX, x, dvLen, Xpts, vars);
  element->testMatDVSensInnerProduct(MASS_MATRIX, x, dvLen, Xpts, vars);

  delete [] x;
}

/*
  The following code tests the element implementation to see if the
  computation of the residual is consistent with the energy
  formulation, to check if the Jacobian matrix is consistent with the
  residual and to test certain design-dependent code.

  Useage:
  ./profile_elements [fd=value]
*/
int main( int argc, char *argv[] ){
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Test only the input element type - if specified, otherwise
  // test everything
  const char *ename = NULL;
  if (argc > 1){
    ename = argv[1];
  }

  const int MAX_NODES = 27;
  const int MAX_VARS_PER_NODE = 8;
  const int MAX_VARS = MAX_NODES*MAX_VARS_PER_NODE;

  // Set the simulation time
  double time = 0.0;

  // Set the variable arrays
  TacsScalar Xpts[3*MAX_NODES];
  TacsScalar vars[MAX_VARS], dvars[MAX_VARS], ddvars[MAX_VARS];
  
  // Generate random arrays
  generate_random_array(Xpts, 3*MAX_NODES, 0.0, 1.0);
  generate_random_array(vars, MAX_VARS);
  generate_random_array(dvars, MAX_VARS);
  generate_random_array(ddvars, MAX_VARS);

  // Set the tolerances depending on whether we're using complex step
  // or not...
#ifdef TACS_USE_COMPLEX
  TACSElement::setFailTolerances(1e-1, 1e-12);
#else
  TACSElement::setFailTolerances(1e-1, 1e-5);
#endif

  // Set the print level
  TACSElement::setPrintLevel(2);

  // Set up the constitutive relationship
  TacsScalar rho = 2700.0, E = 35e4, nu = 0.3, kcorr = 0.8333;
  TacsScalar ys = 434.0e6, t = 0.01;

  // Set the parameter values from the command line
  for ( int i = 0; i < argc; i++ ){
    double Ec = 0.0, rhoc = 0.0;
    if (sscanf(argv[i], "E=%lf", &Ec)){
      if (Ec >= 0.0){
        E = Ec;
        printf("E = %g\n", TacsRealPart(E));
      }
    }
    if (sscanf(argv[i], "rho=%lf", &rhoc)){
      if (rhoc >= 0.0){
        rho = rhoc;
        printf("rho = %g\n", TacsRealPart(rho));
      }
    }
  }

  int dv_num = 14;
  int num_design_vars = dv_num+1;
  FSDTStiffness *fsdt = new isoFSDTStiffness(rho, E, nu, kcorr, ys, t, 
                                             dv_num);
  fsdt->incref();
  TacsScalar axis[3] = {1.0, -1.0, 0.5};
  fsdt->setRefAxis(axis);

  // Allocate and test all the different types of shell elements
  TACSElement *shell = NULL;
  shell = new MITCShell<2>(fsdt, LINEAR);  shell->incref();
  if (!ename || strcmp(ename, shell->elementName()) == 0){
    test_element(shell, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  shell->decref();

  shell = new MITCShell<3>(fsdt, LINEAR);  shell->incref();
  if (!ename || strcmp(ename, shell->elementName()) == 0){
    test_element(shell, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  shell->decref();

  // Nonlinear elements
  shell = new MITCShell<2>(fsdt, NONLINEAR);  shell->incref();
  if (!ename || strcmp(ename, shell->elementName()) == 0){
    test_element(shell, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  shell->decref();

  shell = new MITCShell<3>(fsdt, NONLINEAR);  shell->incref();
  if (!ename || strcmp(ename, shell->elementName()) == 0){
    test_element(shell, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  shell->decref();

  // Large rotation elements
  shell = new MITCShell<2>(fsdt, LARGE_ROTATION);  shell->incref();
  if (!ename || strcmp(ename, shell->elementName()) == 0){
    test_element(shell, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  shell->decref();

  shell = new MITCShell<3>(fsdt, LARGE_ROTATION);  shell->incref();
  if (!ename || strcmp(ename, shell->elementName()) == 0){
    test_element(shell, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  shell->decref();

  // Normalize the variables for unit quaternions
  for ( int i = 0; i < MAX_NODES; i++ ){
    vars[8*i+7] = 0.0;
    TacsScalar *v = &vars[8*i+3];
    TacsScalar fact = 
      1.0/sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
    for ( int j = 0; j < 4; j++ ){
      v[j] *= fact;
    }
  }
  
  MITC9 *mitc9 = new MITC9(fsdt);
  shell = mitc9;
  shell->incref();
  if (!ename || strcmp(ename, shell->elementName()) == 0){
    test_element(shell, time, Xpts, vars, dvars, ddvars, num_design_vars);
    // mitc9->testStrain(Xpts);
  }
  shell->decref();
  fsdt->decref();

  // Set the variables back to a random array
  generate_random_array(dvars, MAX_VARS);

  // Allocate the plane stress stiffness object
  PlaneStressStiffness *ps = new PlaneStressStiffness(rho, E, nu);
  ps->incref();
  
  TACSElement *elem = NULL;

  elem = new PlaneStressTri6(ps);  elem->incref();
  if (!ename || strcmp(ename, elem->elementName()) == 0){
    test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  elem->decref();

  elem = new PlaneStressQuad<2>(ps);  elem->incref();
  if (!ename || strcmp(ename, elem->elementName()) == 0){
    test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  elem->decref();
  
  elem = new PlaneStressQuad<3>(ps);  elem->incref();
  if (!ename || strcmp(ename, elem->elementName()) == 0){
    test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  elem->decref();

  elem = new PlaneStressTri6(ps, NONLINEAR);  elem->incref();
  if (!ename || strcmp(ename, elem->elementName()) == 0){
    test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  elem->decref();

  elem = new PlaneStressQuad<2>(ps, NONLINEAR);  elem->incref();
  if (!ename || strcmp(ename, elem->elementName()) == 0){
    test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  elem->decref();
    
  elem = new PlaneStressQuad<3>(ps, NONLINEAR);  elem->incref();
  if (!ename || strcmp(ename, elem->elementName()) == 0){
    test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  elem->decref();
  ps->decref();

  // Create the solid stiffness classes
  SolidStiffness *stiff = new SolidStiffness(rho, E, nu, ys);
  stiff->incref();
  
  elem = new Solid<2>(stiff); elem->incref();
  if (!ename || strcmp(ename, elem->elementName()) == 0){
    test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  elem->decref();

  elem = new Solid<3>(stiff); elem->incref();
  if (!ename || strcmp(ename, elem->elementName()) == 0){
    test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  elem->decref();
  
  elem = new Solid<2>(stiff, NONLINEAR); elem->incref();
  if (!ename || strcmp(ename, elem->elementName()) == 0){
    test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  elem->decref();

  elem = new Solid<3>(stiff, NONLINEAR); elem->incref();
  if (!ename || strcmp(ename, elem->elementName()) == 0){
    test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  }
  elem->decref();
  stiff->decref();

  // Test the rigid body code within TACS

  if (!ename || strcmp(ename, "Rigid") == 0){
    // Generate a random arrary of variables conforming to the
    // quaternion constraint
    generate_random_array(vars, MAX_VARS);
    generate_random_array(dvars, MAX_VARS);
    generate_random_array(ddvars, MAX_VARS);
    for ( int i = 0; i < MAX_NODES; i++ ){
      vars[8*i+7] = 0.0;
      TacsScalar *v = &vars[8*i+3];

      // Normalize the quaternion constraints
      TacsScalar fact = 
        1.0/sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
      for ( int j = 0; j < 4; j++ ){
        v[j] *= fact;
      }
    }

    // The acceleration due to gravity in global frame of reference
    TACSGibbsVector *gravVec = new TACSGibbsVector(19.0, 10.0, -9.8);

    // Construct the frame of reference
    TACSGibbsVector *rAInitVec = new TACSGibbsVector(5.2, 5.3, 5.4); 
    TACSGibbsVector *rA1Vec = new TACSGibbsVector(5.2+1.0, 5.3, 5.4);
    TACSGibbsVector *rA2Vec = new TACSGibbsVector(5.2, 5.3+1.0, 5.4);
    TACSRefFrame *refFrameA = new TACSRefFrame(rAInitVec, rA1Vec, rA2Vec);

    // Define the inertial properties
    const TacsScalar mA    = 6.0;
    const TacsScalar cA[3] = {20.0, 14.0, 42.0};
    const TacsScalar JA[6] = {1.0, 0.8, -0.7,
                              2.0, 1.4,
                              3.0};
    // Construct a rigid body
    TACSRigidBody *bodyA = new TACSRigidBody(refFrameA,
                                             mA, cA, JA,
                                             rAInitVec, rAInitVec, rAInitVec,
                                             gravVec);
    bodyA->incref();
    test_element(bodyA, time, Xpts, vars, dvars, ddvars, num_design_vars);

    // Test the motion driver element
    TACSGibbsVector *dir = new TACSGibbsVector(0.0, 0.0, 0.1);
    TacsScalar omega = 0.25; // rad/s
    TACSMotionDriver *mDriver = 
      new TACSMotionDriver(dir, omega);
    mDriver->incref();
    test_element(mDriver, time, Xpts, vars, dvars, ddvars, num_design_vars);

    // Test the revolute constraint
    TACSGibbsVector *point = new TACSGibbsVector(0.5, 1.0, -2.5);
    TACSGibbsVector *eRev = new TACSGibbsVector(1.0, -1.0, 1.0);
    TACSRevoluteConstraint *rev = 
      new TACSRevoluteConstraint(bodyA, point, eRev);
    rev->incref();
    test_element(rev, time, Xpts, vars, dvars, ddvars, num_design_vars);

    // Test the cylindrical constraint
    TACSCylindricalConstraint *cyl = 
      new TACSCylindricalConstraint(bodyA, point, eRev);
    cyl->incref();
    test_element(cyl, time, Xpts, vars, dvars, ddvars, num_design_vars);

    // Test the fixed constraint
    TACSFixedConstraint *fix = new TACSFixedConstraint(bodyA, point);
    fix->incref();
    test_element(fix, time, Xpts, vars, dvars, ddvars, num_design_vars);

    // Test the spherical constraint
    TACSSphericalConstraint *ball = new TACSSphericalConstraint(bodyA, point);
    ball->incref();
    test_element(ball, time, Xpts, vars, dvars, ddvars, num_design_vars);
  
    // Test the rigid link code
    TACSRigidLink *rlink = new TACSRigidLink(bodyA);
    rlink->incref();

    // Test the rigid link
    test_element(rlink, time, Xpts, vars, dvars, ddvars, num_design_vars);

    // Define the inertial properties
    const TacsScalar mB    = 2.0;
    const TacsScalar cB[3] = {2.0, 3.0, 4.0};
    const TacsScalar JB[6] = {2.0, 0.60, 0.7,
                              3.0, 0.80,
                              4.0};
  
    // Define dynamics properties
    TACSGibbsVector *rBInitVec = new TACSGibbsVector(3.2, 3.2, 4.3); 
    TACSGibbsVector *rB1Vec = new TACSGibbsVector(3.2+1.0, 3.2, 4.3);
    TACSGibbsVector *rB2Vec = new TACSGibbsVector(3.2, 3.2+1.0, 4.3);
    TACSRefFrame *refFrameB = new TACSRefFrame(rBInitVec, rB1Vec, rB2Vec);

    // Construct a rigid body
    TACSRigidBody *bodyB = new TACSRigidBody(refFrameB,
                                             mB, cB, JB,
                                             rBInitVec, rBInitVec, rBInitVec,
                                             gravVec);
    bodyB->incref();

    // Test the revolute constraint
    TACSRevoluteConstraint *rev2 = 
      new TACSRevoluteConstraint(bodyA, bodyB, point, eRev);
    rev2->incref();
    test_element(rev2, time, Xpts, vars, dvars, ddvars, num_design_vars);

    // Test the revolute constraint
    TACSCylindricalConstraint *cyl2 =
      new TACSCylindricalConstraint(bodyA, bodyB, point, eRev);
    cyl2->incref();
    test_element(cyl2, time, Xpts, vars, dvars, ddvars, num_design_vars);

    // Test the spherical constraint
    TACSSphericalConstraint *ball2 = 
      new TACSSphericalConstraint(bodyA, bodyB, point);
    ball2->incref();
    test_element(ball2, time, Xpts, vars, dvars, ddvars, num_design_vars);
  
    // Decref everything
    cyl->decref();
    fix->decref();
    rev->decref();
    ball->decref();
    bodyA->decref();
    rev2->decref();
    ball2->decref();
    bodyB->decref();
    rlink->decref();
  }

  MPI_Finalize();
  return (0);
}
