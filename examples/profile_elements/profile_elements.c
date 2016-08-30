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
  element->testResidual(time, Xpts, vars, dvars, ddvars);

#ifdef TACS_USE_COMPLEX
  TACSElement::setStepSize(1e-30);
#endif
  element->testJacobian(time, Xpts, vars, dvars, ddvars, 1.0);
  element->testAdjResProduct(x, dvLen, time, Xpts, vars, dvars, ddvars);
  element->testStrainSVSens(Xpts, vars, dvars, ddvars);
  element->testJacobianXptSens(Xpts);

  delete [] x;
}

/*
  The following code testss the most costly pats of s of the element
  computations: the assembly of the stiffness matrices and residuals
  and the computation of the derivative of the residuals w.r.t. the
  nodes.

  These computations are repeated for a series of elements to mimic a
  constant-size mesh with increasingly element order. These
  computational costs, however, do not include the assembly costs of
  placing the entries in the assembled global stiffness matrix or
  residual.

  Useage:
  ./profile_elements [fd=value]
*/
int main( int argc, char * argv[] ){
  // Initialize MPI
  MPI_Init(&argc, &argv); 

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

  // Set the tolerances depending on whether we're using complex step or not...
#ifdef TACS_USE_COMPLEX
  TACSElement::setFailTolerances(1e-1, 1e-12);
#else
  TACSElement::setFailTolerances(1e-1, 1e-5);
#endif

  // Set the print level
  TACSElement::setPrintLevel(2);

  // Set up the constitutive relationship
  TacsScalar rho = 2700.0, E = 35e8, nu = 0.3, kcorr = 0.8333;
  TacsScalar ys = 434.0e6, t = 0.01;
  int dv_num = 14;
  int num_design_vars = dv_num+1;
  FSDTStiffness *fsdt = new isoFSDTStiffness(rho, E, nu, kcorr, ys, t, 
                                             dv_num);
  fsdt->incref();

  // Allocate and test all the different types of shell elements
  TACSElement *shell = NULL;
  shell = new MITCShell<2>(fsdt, LINEAR);  shell->incref();
  test_element(shell, time, Xpts, vars, dvars, ddvars, num_design_vars);
  shell->decref();

  shell = new MITCShell<3>(fsdt, LINEAR);  shell->incref();
  test_element(shell, time, Xpts, vars, dvars, ddvars, num_design_vars);
  shell->decref();

  // Nonlinear elements
  shell = new MITCShell<2>(fsdt, NONLINEAR);  shell->incref();
  test_element(shell, time, Xpts, vars, dvars, ddvars, num_design_vars);
  shell->decref();

  shell = new MITCShell<3>(fsdt, NONLINEAR);  shell->incref();
  test_element(shell, time, Xpts, vars, dvars, ddvars, num_design_vars);
  shell->decref();

  // Large rotation elements
  shell = new MITCShell<2>(fsdt, LARGE_ROTATION);  shell->incref();
  test_element(shell, time, Xpts, vars, dvars, ddvars, num_design_vars);
  shell->decref();

  shell = new MITCShell<3>(fsdt, LARGE_ROTATION);  shell->incref();
  test_element(shell, time, Xpts, vars, dvars, ddvars, num_design_vars);
  shell->decref();

  // Normalize the variables for unit quaternions
  for ( int i = 0; i < MAX_NODES; i++ ){
    vars[8*i+7] = 0.0;
    TacsScalar *v = &vars[8*i+3];
    TacsScalar fact = 
      1.0/sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
  }

  shell = new MITC9(fsdt); shell->incref();
  test_element(shell, time, Xpts, vars, dvars, ddvars, num_design_vars);
  shell->decref();

  fsdt->decref();

  // Set the variables back to a random array
  generate_random_array(dvars, MAX_VARS);

  // Allocate the plane stress stiffness object
  PlaneStressStiffness *ps = new PlaneStressStiffness(rho, E, nu);
  ps->incref();
  
  TACSElement *elem = NULL;
  elem = new PlaneStressTri6(ps);  elem->incref();
  test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  elem->decref();

  elem = new PlaneStressQuad<2>(ps);  elem->incref();
  test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  elem->decref();
  
  elem = new PlaneStressQuad<3>(ps);  elem->incref();
  test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  elem->decref();

  elem = new PlaneStressTri6(ps, NONLINEAR);  elem->incref();
  test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  elem->decref();

  elem = new PlaneStressQuad<2>(ps, NONLINEAR);  elem->incref();
  test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  elem->decref();
    
  elem = new PlaneStressQuad<3>(ps, NONLINEAR);  elem->incref();
  test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  elem->decref();
  ps->decref();
  
  // Create the solid stiffness classes
  SolidStiffness *stiff = new SolidStiffness(rho, E, nu);
  stiff->incref();
  
  elem = new Solid<2>(stiff); elem->incref();
  test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  elem->decref();

  elem = new Solid<3>(stiff); elem->incref();
  test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  elem->decref();
  
  elem = new Solid<2>(stiff, NONLINEAR); elem->incref();
  test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  elem->decref();

  elem = new Solid<3>(stiff, NONLINEAR); elem->incref();
  test_element(elem, time, Xpts, vars, dvars, ddvars, num_design_vars);
  elem->decref();
  stiff->decref();

  MPI_Finalize();
  return (0);
}
