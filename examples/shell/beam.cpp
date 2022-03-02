#include "TACSElementVerification.h"
#include "TACSBeamConstitutive.h"
#include "TACSShellElementDefs.h"

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);

  // Get the rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  TacsTestBeamModelDerivatives<6, TACSBeamBasis<3>, TACSBeamLinearModel>();  

  MPI_Finalize();
  return 0;
}
