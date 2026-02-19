#include "TACSBlockCyclicMat.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  TacsInitialize();

  // Set the communicator
  MPI_Comm comm = MPI_COMM_WORLD;

  // Create a square matrix
  int nrows = 50;
  TACSBlockCyclicMat *mat = new TACSBlockCyclicMat(comm, nrows, nrows);
  mat->incref();

  // Set Random values into the matrix
  srand(0);
  mat->setRand();

  // Get the size of the local vector
  int size = mat->getLocalVecSize();

  // Allocate the vectors for the size of matrix
  TacsScalar *x = new TacsScalar[size];
  TacsScalar *y = new TacsScalar[size];
  for (int i = 0; i < size; i++) {
    x[i] = 1.0;
  }

  // Multiply to get the right size
  double tm = MPI_Wtime();
  mat->mult(x, y);
  tm = MPI_Wtime() - tm;
  memcpy(x, y, size * sizeof(TacsScalar));

  // Factor the matrix
  double tf = MPI_Wtime();
  mat->factor();
  tf = MPI_Wtime() - tf;

  // Solve a linear system with the matrix
  double ta = MPI_Wtime();
  mat->applyFactor(y);
  ta = MPI_Wtime() - ta;

  // Compute the error
  TacsScalar error[2] = {0.0, 0.0};
  for (int i = 0; i < size; i++) {
    error[0] += (y[i] - 1.0) * (y[i] - 1.0);
    error[1] += x[i];
  }

  // Create
  TacsScalar err[2];
  MPI_Reduce(error, err, 2, TACS_MPI_TYPE, MPI_SUM, 0, comm);

  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    printf("Multiply time:     %15.10e\n", tm);
    printf("Factor time:       %15.10e\n", tf);
    printf("Apply factor time: %15.10e\n", ta);
    printf("Error:             %15.10e\n", TacsRealPart(sqrt(err[0])));
    printf("e^{T}*A*e:         %15.10e\n", TacsRealPart(err[1]));
  }

  TacsFinalize();
  MPI_Finalize();
  return 0;
}
