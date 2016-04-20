#include "RigidBody.h"
#include "TACSElementAlgebra.h"

int main( int argc, char *argv[] ){
  TacsScalar mass = 0.0;
  TacsScalar c[] = {1.0, 0.35, -0.5};
  TacsScalar J[] = {1.0, -0.25, 0.0,
                    2.0, 0.0,
                    0.75};

  TACSRigidBody *body = new TACSRigidBody(mass, c, J);
  body->incref();
  body->testResidual(1e-4);
  body->decref();

  return (0);
}
