#include "TimoshenkoStiffness.h"
#include "MITC3.h"


int main( int argc, char *argv[] ){
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Set the reference axis
  TacsScalar axis[] = {0.0, 1.0, 0.0};
  
  TacsScalar rhoA = 1.5;
  TacsScalar rhoIy = 0.15;
  TacsScalar rhoIz = 0.34;
  TacsScalar rhoIyz = -0.02;

  TacsScalar EA = 1e4;
  TacsScalar GJ = 1.50e4;
  TacsScalar EIy = 2.4e4;
  TacsScalar EIz = 3.24e4;
  TacsScalar kGAy = 2.5e3;
  TacsScalar kGAz = 5.2e3;
  
  // Create the Timoshenko stiffness object
  TimoshenkoStiffness *stiff =
    new TimoshenkoStiffness(rhoA, rhoIy, rhoIz, rhoIyz,
                            EA, GJ, EIy, EIz, kGAy, kGAz,
                            axis);
  stiff->incref();

  // Create the MITC3 element
  MITC3 *beam = new MITC3(stiff);
  beam->incref();

  TacsScalar X[] = {0.0, 0.0, 0.0,
                    0.375, 0.375, 0.1,
                    1.0, 1.0, 0.2};
  beam->testStrain(X);

  TacsScalar vars[24], dvars[24], ddvars[24];
  for ( int i = 0; i < 24; i++ ){
    vars[i] = -1.0 + 2.0*rand()/RAND_MAX;
    dvars[i] = -1.0 + 2.0*rand()/RAND_MAX;
    ddvars[i] = -1.0 + 2.0*rand()/RAND_MAX;
  }

  vars[7] = vars[15] = vars[23] = 0.0;

  beam->setStepSize(1e-5);
  beam->setPrintLevel(2);
  beam->testResidual(0.0, X, vars, dvars, ddvars);
  beam->testJacobian(0.0, X, vars, dvars, ddvars);

  beam->decref();
  stiff->decref();

  MPI_Finalize();
  return 0;
}
