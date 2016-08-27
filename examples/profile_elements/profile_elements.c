#include "isoFSDTStiffness.h"
#include "MITCShell.h"

/*
  Set the nodes for a 2D Cartesian-product element type of given order
*/
void set_nodes( TacsScalar Xpts[], int order ){
  for ( int j = 0; j < order; j++ ){
    for ( int i = 0; i < order; i++ ){
      Xpts[3*(i + order*j)]   =   1.0*i + 0.2*j - 0.15*(i*i - j*j);
      Xpts[3*(i + order*j)+1] = -0.25*i + 1.0*j + 0.2*(i*i*i + j*j*j);
      Xpts[3*(i + order*j)+2] =   0.1*i - 0.1*j + 0.5*(j*i + i*i + j*j);
    }
  }
}

/*
  Profile the element assembly operations and derivative computations
  for a series of elements with different options. This is designed to
  estimate the cost of the computations without assembly operations.

  The cost of the assembly operations is for 'nelems/(order-1)^2'
  repeated computations. This would be the number of elements required
  for a constant size mesh.

  The profile code also tests the element implementations using either
  finite-difference, if compiled in real mode, or complex-step, if
  compiled with complex variables.

  input:
  con:     the constitutive object
  nelems:  the number of elements in the nominal 'mesh'
  type:    the type of test to run [1 through 6]
  dh:      finite-difference step size
*/
void profile_elements( FSDTStiffness *con, 
		       int nelems, int type, double dh ){
  // Now go through and perform the profiling for the displacement based elements
  TACSElement *shell2 = NULL, *shell3 = NULL, *shell4 = NULL;
  if (type == 1){
    shell2 = new DShell<2>(con, DShell<2>::NONLINEAR);
    shell3 = new DShell<3>(con, DShell<3>::NONLINEAR);
    shell4 = new DShell<4>(con, DShell<4>::NONLINEAR);
    printf("Nonlinear DShell element\n");
  }
  else if (type == 2){
    shell2 = new MITCShell<2>(con, MITCShell<2>::LINEAR);
    shell3 = new MITCShell<3>(con, MITCShell<3>::LINEAR);
    shell4 = new MITCShell<4>(con, MITCShell<4>::LINEAR);
    printf("Linear MITC element\n");
  }
  else if (type == 3){
    shell2 = new MITCShell<2>(con, MITCShell<2>::NONLINEAR);
    shell3 = new MITCShell<3>(con, MITCShell<3>::NONLINEAR);
    shell4 = new MITCShell<4>(con, MITCShell<4>::NONLINEAR);
    printf("Nonlinear MITC element\n");
  }
  else if (type == 4){
    shell2 = new DShell<2>(con, DShell<2>::LARGE_ROTATION);
    shell3 = new DShell<3>(con, DShell<3>::LARGE_ROTATION);
    shell4 = new DShell<4>(con, DShell<4>::LARGE_ROTATION);
    printf("Nonlinear DShell large rotation element\n");
  }
  else if (type == 5){
    shell2 = new MITCShell<2>(con, MITCShell<2>::LARGE_ROTATION);
    shell3 = new MITCShell<3>(con, MITCShell<3>::LARGE_ROTATION);
    shell4 = new MITCShell<4>(con, MITCShell<4>::LARGE_ROTATION);
    printf("Nonlinear MITC large rotation element\n");
  }
  else { // (type == 0){
    shell2 = new DShell<2>(con, DShell<2>::LINEAR);
    shell3 = new DShell<3>(con, DShell<3>::LINEAR);
    shell4 = new DShell<4>(con, DShell<4>::LINEAR);
    printf("Linear DShell element\n");
  }

  shell2->incref();
  shell3->incref();
  shell4->incref();

  // Set the max number of elements
  const int NUM_VARS = 16*6;
  const int NUM_NODES = 16*3;
  TacsScalar vars[NUM_VARS], res[NUM_VARS], mat[NUM_VARS*NUM_VARS];
  TacsScalar Xpts[NUM_NODES];
  TacsScalar resXptSens[NUM_VARS*NUM_NODES];

  // Set random variables
  memset(vars, 0, NUM_VARS*sizeof(TacsScalar));
  for ( int k = 0; k < NUM_VARS; k++ ){
    vars[k] = 1.0*rand()/RAND_MAX;
  }

  set_nodes(Xpts, 2);
  double t2 = MPI_Wtime();
  for ( int k = 0; k < nelems; k++ ){
    shell2->getMat(mat, res, vars, Xpts, NORMAL);
  }
  t2 = MPI_Wtime() - t2;
  
  set_nodes(Xpts, 3);
  double t3 = MPI_Wtime();
  for ( int k = 0; k < nelems/4; k++ ){
    shell3->getMat(mat, res, vars, Xpts, NORMAL);
  }
  t3 = MPI_Wtime() - t3;
  
  set_nodes(Xpts, 4);
  double t4 = MPI_Wtime();
  for ( int k = 0; k < nelems/9; k++ ){
    shell4->getMat(mat, res, vars, Xpts, NORMAL);
  }
  t4 = MPI_Wtime() - t4;
  
  printf("Time for assembly of %d 2nd order %s elements: %f\n", 
	 nelems, shell2->elementName(), t2);
  printf("Time for assembly of %d 3rd order %s elements: %f\n", 
	 nelems/4, shell3->elementName(), t3);
  printf("Time for assembly of %d 4th order %s elements: %f\n", 
	 nelems/9, shell4->elementName(), t4);
  
  set_nodes(Xpts, 2);
  t2 = MPI_Wtime();
  for ( int k = 0; k < nelems; k++ ){
    shell2->getResXptSens(resXptSens, vars, Xpts);
  }
  t2 = MPI_Wtime() - t2;
  
  set_nodes(Xpts, 3);
  t3 = MPI_Wtime();
  for ( int k = 0; k < nelems/4; k++ ){
    shell3->getResXptSens(resXptSens, vars, Xpts);
  }
  t3 = MPI_Wtime() - t3;
  
  set_nodes(Xpts, 4);
  t4 = MPI_Wtime();
  for ( int k = 0; k < nelems/9; k++ ){
    shell4->getResXptSens(resXptSens, vars, Xpts);
  }
  t4 = MPI_Wtime() - t4;
  
  printf("Time for dR/dx of %d 2nd order %s elements: %f\n", 
	 nelems, shell2->elementName(), t2);
  printf("Time for dR/dx of %d 3rd order %s elements: %f\n", 
	 nelems/4, shell3->elementName(), t3);
  printf("Time for dR/dx of %d 4th order %s elements: %f\n", 
	 nelems/9, shell4->elementName(), t4);

  // Test the elements using finite-difference or complex step
  double pt[] = {0.0, 0.0};

  int default_print = 1;
  set_nodes(Xpts, 2);
  TestElement * test2 = new TestElement(shell2, NULL, vars, Xpts);
  test2->incref();
  test2->setPrintLevel(default_print);
  test2->setStepSize(dh);
  if (test2->testResXptSens()){
    test2->setPrintLevel(2);
    test2->testResXptSens();
    test2->setPrintLevel(default_print);
  }
  if (test2->testStrainXptSens(pt)){
    test2->setPrintLevel(2);
    test2->testStrainXptSens(pt);
    test2->setPrintLevel(default_print);
  }
  if (test2->testStrainSVSens(pt)){
    test2->setPrintLevel(2);
    test2->testStrainSVSens(pt);
    test2->setPrintLevel(default_print);
  }
  if (test2->testStiffnessMat()){
    test2->setPrintLevel(2);
    test2->testStiffnessMat();
    test2->setPrintLevel(default_print);
  }
  test2->decref();
  
  set_nodes(Xpts, 3);
  TestElement * test3 = new TestElement(shell3, NULL, vars, Xpts);
  test3->incref();
  test3->setPrintLevel(default_print);
  test3->setStepSize(dh);
  if (test3->testResXptSens()){
    test3->setPrintLevel(2);
    test3->testResXptSens();
    test3->setPrintLevel(default_print);
  }
  if (test3->testStrainXptSens(pt)){
    test3->setPrintLevel(2);
    test3->testStrainXptSens(pt);
    test3->setPrintLevel(default_print);
  }
  if (test3->testStrainSVSens(pt)){
    test3->setPrintLevel(2);
    test3->testStrainSVSens(pt);
    test3->setPrintLevel(default_print);
  }
  if (test3->testStiffnessMat()){
    test3->setPrintLevel(2);
    test3->testStiffnessMat();
    test3->setPrintLevel(default_print);
  }
  test3->decref();
  
  set_nodes(Xpts, 4);
  TestElement * test4 = new TestElement(shell4, NULL, vars, Xpts);
  test4->incref();
  test4->setPrintLevel(default_print);
  test4->setStepSize(dh);
  if (test4->testResXptSens()){
    test4->setPrintLevel(2);
    test4->testResXptSens();
    test4->setPrintLevel(default_print);
  }
  if (test4->testStrainXptSens(pt)){
    test4->setPrintLevel(2);
    test4->testStrainXptSens(pt);
    test4->setPrintLevel(default_print);
  }
  if (test4->testStrainSVSens(pt)){
    test4->setPrintLevel(2);
    test4->testStrainSVSens(pt);
    test4->setPrintLevel(default_print);
  }
  if (test4->testStiffnessMat()){
    test4->setPrintLevel(2);
    test4->testStiffnessMat();
    test4->setPrintLevel(default_print);
  }
  test4->decref();

  shell2->decref();
  shell3->decref();
  shell4->decref();
}

/*
  The following code profiles the most costly pats of s of the element
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
  
  // Set up the constitutive relationship
  TacsScalar rho = 2700.0, E = 35e8, nu = 0.3, kcorr = 0.8333;
  TacsScalar ys = 434.0e6, t = 0.01;
  FSDTStiffness * con = new isoFSDTStiffness(rho, E, nu, kcorr, ys, t);
  con->incref();

  // This number is divisible by 9 and 4 so that the number of
  // elements can be nested for 2nd, 3rd and 4th order meshes
  int nelems = 9*4*10;

  double dh = 1e-6;
  for ( int k = 0; k < argc; k++ ){
    if (sscanf(argv[k], "fd=%le", &dh) == 1){
      printf("Finite-difference interval: %e\n", dh);
    }
  } 

  // Profile all the element types desired
  profile_elements(con, nelems, 0, dh);
  profile_elements(con, nelems, 1, dh);
  profile_elements(con, nelems, 2, dh);
  profile_elements(con, nelems, 3, dh);
  profile_elements(con, nelems, 4, dh);
  profile_elements(con, nelems, 5, dh);

  MPI_Finalize();
  return (0);
}
