#include "RigidBody.h"
#include "TACSElementAlgebra.h"

int main( int argc, char *argv[] ){
  TacsScalar mass = 4.0;
  TacsScalar c[] = {1.0, 0.35, -0.5};
  TacsScalar J[] = {1.0, -0.25, 0.24,
                    2.0, 0.34,
                    0.75};

  // Zero the Jacobian coefficients
  double alpha = 1.26, beta = 0.35, gamma = 4.34;
  beta = gamma = 0.0;

  TACSRigidBody *body = new TACSRigidBody(mass, c, J);
  body->incref();
  body->testResidual(1e-4);
  body->testJacobian(1e-6, alpha, beta, gamma);
  body->decref();

  return (0);
}
