#include "TACSAssembler.h"
#include "TACSCreator.h"
#include "TACSElementAlgebra.h"
#include "TACSElementVerification.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSShellElementDefs.h"
#include "TACSShellTraction.h"
#include "TACSToFH5.h"
#include "tacslapack.h"

/*
  Take the difference between what's in the vector and what the exact
  solution is at each point in the mesh.
*/
void setExactSolution(TACSAssembler *assembler, TACSBVec *ans, TacsScalar R,
                      double alpha, double beta, TacsScalar U, TacsScalar V,
                      TacsScalar W, TacsScalar theta, TacsScalar phi) {
  // Get the x,y,z locations of the nodes from TACS
  int nnodes = assembler->getNumNodes();

  TACSBVec *X = assembler->createNodeVec();
  X->incref();
  assembler->getNodes(X);

  // Get the nodal array
  TacsScalar *Xpts;
  X->getArray(&Xpts);

  // Get the finite-element solution from TACS
  TacsScalar *uans;
  ans->getArray(&uans);

  for (int node = 0; node < nnodes; node++) {
    // Compute the u, v, w, theta, phi solution at this point within
    // the shell
    TacsScalar x = Xpts[3 * node];
    TacsScalar y = -R * atan2(TacsRealPart(Xpts[3 * node + 2]),
                              TacsRealPart(Xpts[3 * node + 1]));

    TacsScalar u = U * sin(alpha * y) * cos(beta * x);
    TacsScalar v = V * cos(alpha * y) * sin(beta * x);
    TacsScalar w = W * sin(alpha * y) * sin(beta * x);
    TacsScalar psi_x = theta * sin(alpha * y) * cos(beta * x);
    TacsScalar psi_y = phi * cos(alpha * y) * sin(beta * x);

    TacsScalar n[2], t[2];
    n[0] = cos(y / R);
    n[1] = -sin(y / R);

    t[0] = -sin(y / R);
    t[1] = -cos(y / R);

    // Transform the displacements/rotations into the global ref. frame
    uans[0] = u;
    uans[1] = w * n[0] + v * t[0];
    uans[2] = w * n[1] + v * t[1];
    uans[3] = -psi_y;
    uans[4] = psi_x * t[0];
    uans[5] = psi_x * t[1];

    uans += 6;
  }

  X->decref();
}

/*
  Set the solution to be a linearly varying solution.

  This can be used to check the strain expressions for the cylinder.
*/
void setLinearSolution(TACSAssembler *assembler, TACSBVec *ans, TacsScalar R,
                       TacsScalar coef[]) {
  TACSBVec *X = assembler->createNodeVec();
  X->incref();
  assembler->getNodes(X);

  // Get the x,y,z locations of the nodes from TACS
  int nnodes = assembler->getNumNodes();
  TacsScalar *Xpts;
  X->getArray(&Xpts);

  // Get the finite-element solution from TACS
  TacsScalar *uans;
  ans->getArray(&uans);

  // Set the coefficients from the array
  TacsScalar u0, ux, uy, v0, vx, vy, w0, wx, wy;
  TacsScalar psi_x0, psi_xx, psi_xy;
  TacsScalar psi_y0, psi_yx, psi_yy;

  u0 = coef[0], ux = coef[1], uy = coef[2];
  v0 = coef[3], vx = coef[4], vy = coef[5];
  w0 = coef[6], wx = coef[7], wy = coef[8];
  psi_x0 = coef[9], psi_xx = coef[10], psi_xy = coef[11];
  psi_y0 = coef[12], psi_yx = coef[13], psi_yy = coef[14];

  for (int node = 0; node < nnodes; node++) {
    // Compute the u, v, w, theta, phi solution at this point within
    // the shell
    TacsScalar x = Xpts[3 * node];
    TacsScalar y = -R * atan2(TacsRealPart(Xpts[3 * node + 2]),
                              TacsRealPart(Xpts[3 * node + 1]));

    TacsScalar n[2], t[2];
    n[0] = cos(y / R);
    n[1] = -sin(y / R);
    t[0] = -sin(y / R);
    t[1] = -cos(y / R);

    // Transform the displacements/rotations into the global ref. frame
    uans[0] = u0 + x * ux + y * uy;
    uans[1] = (w0 + x * wx + y * wy) * n[0] + (v0 + x * vx + y * vy) * t[0];
    uans[2] = (w0 + x * wx + y * wy) * n[1] + (v0 + x * vx + y * vy) * t[1];
    uans[3] = -(psi_y0 + psi_yx * x + psi_yy * y);
    uans[4] = (psi_x0 + psi_xx * x + psi_xy * y) * t[0];
    uans[5] = (psi_x0 + psi_xx * x + psi_xy * y) * t[1];

    uans += 6;
  }

  X->decref();
}

/*
  Compute the coefficients of a single term in a Fourier series for a
  specially orthotropic cylinder subjectedt to a sinusoidally varying
  pressure distribution specified as follows:

  p = sin(alpha*y)*cos(beta*x)

  Note that y = r*theta

  The coefficients U, V, W, theta and phi are for the expressions:

  u(x,y) = U*sin(alpha*y)*cos(beta*x)
  v(x,y) = V*cos(alpha*y)*sin(beta*x)
  w(x,y) = W*sin(alpha*y)*sin(beta*x)
  psi_x(x,y) = theta*sin(alpha*y)*cos(beta*x)
  psi_y(x,y) = phi*cos(alpha*y)*sin(beta*x)

  Note that u, v and w correspond to the axial, tangential and normal
  displacements in a reference frame attached to the shell surface -
  they are not the displacements expressed in the global reference
  frame. Likewise, psi_x and psi_y are the rotations of the normal
  along the x-direction and the tangential y-direction.
*/
void computeCoefficients(TacsScalar *U, TacsScalar *V, TacsScalar *W,
                         TacsScalar *theta, TacsScalar *phi, double alpha,
                         double beta, TacsScalar ainv, TacsScalar A11,
                         TacsScalar A12, TacsScalar A22, TacsScalar A33,
                         TacsScalar D11, TacsScalar D12, TacsScalar D22,
                         TacsScalar D33, TacsScalar bA11, TacsScalar bA22,
                         TacsScalar load) {
  TacsScalar A[5 * 5], A2[5 * 5];  // 5 x 5 system of linear equations
  int ipiv[5];
  TacsScalar rhs[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  rhs[2] = -load;

  TacsScalar B[8 * 5];
  memset(B, 0, sizeof(B));
  TacsScalar *b = B;

  // Assign the columns for u
  b[0] = -beta;
  b[2] = alpha;
  b[5] = -alpha * ainv;
  b += 8;

  // Assign columns for v
  b[1] = -alpha;
  b[2] = beta;
  b[4] = alpha * ainv;
  b[6] = -ainv;
  b += 8;

  // Assign the columns for w
  b[1] = ainv;
  b[4] = -ainv * ainv;
  b[6] = alpha;
  b[7] = beta;
  b += 8;

  // Assign the columns for psi_x
  b[3] = -beta;
  b[5] = alpha;
  b[7] = 1.0;
  b += 8;

  // Assign the columns for psi_y
  b[4] = -alpha;
  b[5] = beta;
  b[6] = 1.0;

  for (int j = 0; j < 5; j++) {
    TacsScalar *bj = &B[8 * j];
    for (int i = 0; i < 5; i++) {
      TacsScalar *bi = &B[8 * i];

      A2[i + 5 * j] =
          -((bi[0] * (A11 * bj[0] + A12 * bj[1]) +
             bi[1] * (A12 * bj[0] + A22 * bj[1]) + bi[2] * A33 * bj[2]) +
            (bi[3] * (D11 * bj[3] + D12 * bj[4]) +
             bi[4] * (D12 * bj[3] + D22 * bj[4]) + bi[5] * D33 * bj[5]) +
            bi[6] * bA11 * bj[6] + bi[7] * bA22 * bj[7]);
    }
  }

  // The first equation for u
  A[0] = -(A11 * beta * beta + A33 * alpha * alpha) -
         D33 * ainv * ainv * alpha * alpha;
  A[5] = -(A33 + A12) * alpha * beta;
  A[10] = A12 * beta * ainv;
  A[15] = D33 * ainv * alpha * alpha;
  A[20] = D33 * ainv * alpha * beta;

  // The second equation for v
  A[1] = -(A12 + A33) * alpha * beta;
  A[6] = -(A33 * beta * beta + A22 * alpha * alpha) - ainv * ainv * bA11 -
         D22 * ainv * ainv * alpha * alpha;
  A[11] = (A22 + bA11) * ainv * alpha + D22 * alpha * ainv * ainv * ainv;
  A[16] = D12 * ainv * alpha * beta;
  A[21] = bA11 * ainv + D22 * ainv * alpha * alpha;

  // The third equation for w
  A[2] = A12 * beta * ainv;
  A[7] = (bA11 + A22) * alpha * ainv + D22 * alpha * ainv * ainv * ainv;
  A[12] = -(bA11 * alpha * alpha + bA22 * beta * beta) - A22 * ainv * ainv -
          D22 * ainv * ainv * ainv * ainv;
  A[17] = -bA22 * beta - D12 * beta * ainv * ainv;
  A[22] = -bA11 * alpha - D22 * alpha * ainv * ainv;

  // Fourth equation for theta
  A[3] = D33 * ainv * alpha * alpha;
  A[8] = D12 * ainv * alpha * beta;
  A[13] = -bA22 * beta - D12 * beta * ainv * ainv;
  A[18] = -(D11 * beta * beta + D33 * alpha * alpha) - bA22;
  A[23] = -(D12 + D33) * alpha * beta;

  // Fifth equation for phi
  A[4] = D33 * ainv * alpha * beta;
  A[9] = bA11 * ainv + D22 * ainv * alpha * alpha;
  A[14] = -bA11 * alpha - D22 * alpha * ainv * ainv;
  A[19] = -(D33 + D12) * alpha * beta;
  A[24] = -(D33 * beta * beta + D22 * alpha * alpha) - bA11;

  int print_equations = 1;
  if (print_equations) {
    const char *variables[] = {"U", "V", "W", "theta", "phi"};
    for (int k = 0; k < 5; k++) {
      printf("Equation %d\n", k);
      for (int j = 0; j < 5; j++) {
        printf(" %15.5f %s", TacsRealPart(A[k + 5 * j] - A2[k + 5 * j]),
               variables[j]);
      }
      printf("\n");
    }
  }

  int info = 0;
  int n = 5;
  int nrhs = 1;
  LAPACKgesv(&n, &nrhs, A2, &n, ipiv, rhs, &n, &info);

  // Solve for the coefficients
  *U = rhs[0];
  *V = rhs[1];
  *W = rhs[2];
  *theta = rhs[3];
  *phi = rhs[4];
}

/*
  Create the TACSAssembler object and return the associated TACS
  creator object
*/
void createAssembler(MPI_Comm comm, int order, int nx, int ny,
                     TACSElement *element, TACSAssembler **_assembler,
                     TACSCreator **_creator) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  double load = 1.0;
  double L = 100.0;
  double R = 100.0 / M_PI;
  double defect = 0.1;

  // Set the alpha and beta parameters
  double alpha = 4.0 / R;
  double beta = 3 * M_PI / L;

  // Set the number of nodes/elements on this proc
  int varsPerNode = element->getVarsPerNode();

  // Set up the creator object
  TACSCreator *creator = new TACSCreator(comm, varsPerNode);

  if (rank == 0) {
    // Set the number of elements
    int nnx = (order - 1) * nx + 1;
    int nny = (order - 1) * ny;
    int numNodes = nnx * nny;
    int numElements = nx * ny;

    // Allocate the input arrays into the creator object
    int *ids = new int[numElements];
    int *ptr = new int[numElements + 1];
    int *conn = new int[order * order * numElements];

    // Set the element identifiers to all zero
    memset(ids, 0, numElements * sizeof(int));

    ptr[0] = 0;
    for (int k = 0; k < numElements; k++) {
      // Back out the i, j coordinates from the corresponding
      // element number
      int i = k % nx;
      int j = k / nx;

      // Set the node connectivity
      for (int jj = 0; jj < order; jj++) {
        for (int ii = 0; ii < order; ii++) {
          if (j == ny - 1 && jj == order - 1) {
            conn[order * order * k + ii + order * jj] = ((order - 1) * i + ii);
          } else {
            conn[order * order * k + ii + order * jj] =
                ((order - 1) * i + ii) + ((order - 1) * j + jj) * nnx;
          }
        }
      }

      ptr[k + 1] = order * order * (k + 1);
    }

    // Set the connectivity
    creator->setGlobalConnectivity(numNodes, numElements, ptr, conn, ids);
    delete[] conn;
    delete[] ptr;
    delete[] ids;

    int numBcs = 2 * nny;
    int *bcNodes = new int[numBcs];
    int k = 0;

    for (int j = 0; j < nny; j++) {
      int node = j * nnx;
      bcNodes[k] = node;
      k++;

      node = nnx - 1 + j * nnx;
      bcNodes[k] = node;
      k++;
    }

    // Following copied from TACSCreator
    int *bc_ptr = new int[numBcs + 1];

    // Since the bc_vars array is input as NULL, assume that
    // all the variables at this node are fully restrained.
    int *bc_vars = new int[numBcs * 3 + 1];
    bc_ptr[0] = 0;

    for (int i = 0; i < numBcs; i++) {
      bc_ptr[i + 1] = bc_ptr[i];
      if (i == 0) {
        // Fix u, v, w, and rotx
        for (int j = 0; j < 4; j++) {
          bc_vars[bc_ptr[i + 1]] = j;
          bc_ptr[i + 1]++;
        }
      } else {
        // Fix v, w, and rotx
        for (int j = 0; j < 3; j++) {
          bc_vars[bc_ptr[i + 1]] = j + 1;
          bc_ptr[i + 1]++;
        }
      }
    }

    TacsScalar *bc_vals = new TacsScalar[bc_ptr[numBcs]];
    memset(bc_vals, 0, bc_ptr[numBcs] * sizeof(TacsScalar));

    // Set the boundary conditions
    creator->setBoundaryConditions(numBcs, bcNodes, bc_ptr, bc_vars, bc_vals);

    delete[] bcNodes;
    delete[] bc_ptr;
    delete[] bc_vars;
    delete[] bc_vals;

    // Set the node locations
    TacsScalar *Xpts = new TacsScalar[3 * numNodes];

    for (int j = 0; j < nny; j++) {
      double v = -M_PI + (2.0 * M_PI * j) / nny;
      for (int i = 0; i < nnx; i++) {
        double u = 1.0 * i / (nnx - 1);
        double theta =
            v + 0.25 * M_PI * u + defect * sin(v) * cos(2 * M_PI * u);
        double x = L * (u + defect * cos(v) * sin(2 * M_PI * u));
        // double theta = v;
        // double x = L*u;

        int node = i + j * nnx;
        Xpts[3 * node] = x;
        Xpts[3 * node + 1] = R * cos(theta);
        Xpts[3 * node + 2] = -R * sin(theta);
      }
    }

    // Set the nodal locations
    creator->setNodes(Xpts);
    delete[] Xpts;
  }

  // Set the one element
  creator->setElements(1, &element);

  // Set the reordering type
  creator->setReorderingType(TACSAssembler::MULTICOLOR_ORDER,
                             TACSAssembler::GAUSS_SEIDEL);

  // Create TACS
  TACSAssembler *assembler = creator->createTACS();

  // Set the elements the node vector
  TACSBVec *X = assembler->createNodeVec();
  X->incref();
  assembler->getNodes(X);

  TACSAuxElements *aux = new TACSAuxElements();

  for (int elem = 0; elem < assembler->getNumElements(); elem++) {
    TacsScalar Xelem[3 * 9];
    TacsScalar tr[3 * 9];

    int nnodes;
    const int *nodes;
    assembler->getElement(elem, &nnodes, &nodes);
    X->getValues(nnodes, nodes, Xelem);

    for (int node = 0; node < nnodes; node++) {
      // Compute the pressure at this point in the shell
      double x = TacsRealPart(Xelem[3 * node]);
      double y = -R * atan2(TacsRealPart(Xelem[3 * node + 2]),
                            TacsRealPart(Xelem[3 * node + 1]));

      TacsScalar pval = load * sin(beta * x) * sin(alpha * y);
      TacsScalar ynorm = Xelem[3 * node + 1] / R;
      TacsScalar znorm = Xelem[3 * node + 2] / R;

      tr[3 * node] = 0.0;
      tr[3 * node + 1] = ynorm * pval;
      tr[3 * node + 2] = znorm * pval;
    }

    TACSElement *trac = NULL;
    if (order == 2) {
      trac = new TACSShellTraction<6, TACSQuadLinearQuadrature,
                                   TACSShellQuadBasis<2> >(tr);
    } else if (order == 3) {
      trac = new TACSShellTraction<6, TACSQuadQuadraticQuadrature,
                                   TACSShellQuadBasis<3> >(tr);
    }

    aux->addElement(elem, trac);
  }

  X->decref();

  // Set the auxiliary elements
  assembler->setAuxElements(aux);

  // Set the pointers
  *_assembler = assembler;
  *_creator = creator;
}

/*
  Create the TACSAssembler object and return the associated TACS
  creator object
*/
void createTriAssembler(MPI_Comm comm, int order, int nx, int ny,
                        TACSElement *element, TACSAssembler **_assembler,
                        TACSCreator **_creator) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  double load = 1.0;
  double L = 100.0;
  double R = 100.0 / M_PI;
  double defect = 0.1;

  // Set the alpha and beta parameters
  double alpha = 4.0 / R;
  double beta = 3 * M_PI / L;

  // Set the number of nodes/elements on this proc
  int varsPerNode = element->getVarsPerNode();

  // Set up the creator object
  TACSCreator *creator = new TACSCreator(comm, varsPerNode);

  if (rank == 0) {
    // Set the number of elements
    int nnx = (order - 1) * nx + 1;
    int nny = (order - 1) * ny;
    int numNodes = nnx * nny;
    int numElements = 2 * nx * ny;

    // Allocate the input arrays into the creator object
    int *ids = new int[numElements];
    int *ptr = new int[numElements + 1];
    int *conn = new int[6 * numElements];

    // Set the element identifiers to all zero
    memset(ids, 0, numElements * sizeof(int));

    ptr[0] = 0;
    for (int k = 0; k < numElements / 2; k++) {
      // Back out the i, j coordinates from the corresponding
      // element number
      int i = k % nx;
      int j = k / nx;

      // Set the node connectivity
      int nodes[9];
      for (int jj = 0; jj < order; jj++) {
        for (int ii = 0; ii < order; ii++) {
          if (j == ny - 1 && jj == order - 1) {
            nodes[ii + order * jj] = ((order - 1) * i + ii);
          } else {
            nodes[ii + order * jj] =
                ((order - 1) * i + ii) + ((order - 1) * j + jj) * nnx;
          }
        }
      }

      // Set the connectivity from the nodes
      if (order == 2) {
        int *c = &conn[6 * k];
        c[0] = nodes[0];
        c[1] = nodes[1];
        c[2] = nodes[3];
        ptr[2 * k + 1] = 3 * (2 * k + 1);

        c = &conn[6 * k + 3];
        c[0] = nodes[0];
        c[1] = nodes[3];
        c[2] = nodes[2];
        ptr[2 * k + 2] = 3 * (2 * k + 2);
      } else if (order == 3) {
        int *c = &conn[12 * k];
        c[0] = nodes[0];
        c[1] = nodes[2];
        c[2] = nodes[8];
        c[3] = nodes[1];
        c[4] = nodes[5];
        c[5] = nodes[4];
        ptr[2 * k + 1] = 6 * (2 * k + 1);

        c = &conn[12 * k + 6];
        c[0] = nodes[0];
        c[1] = nodes[8];
        c[2] = nodes[6];
        c[3] = nodes[4];
        c[4] = nodes[7];
        c[5] = nodes[3];
        ptr[2 * k + 2] = 6 * (2 * k + 2);
      }
    }

    // Set the connectivity
    creator->setGlobalConnectivity(numNodes, numElements, ptr, conn, ids);
    delete[] conn;
    delete[] ptr;
    delete[] ids;

    int numBcs = 2 * nny;
    int *bcNodes = new int[numBcs];
    int k = 0;

    for (int j = 0; j < nny; j++) {
      int node = j * nnx;
      bcNodes[k] = node;
      k++;

      node = nnx - 1 + j * nnx;
      bcNodes[k] = node;
      k++;
    }

    // Following copied from TACSCreator
    int *bc_ptr = new int[numBcs + 1];

    // Since the bc_vars array is input as NULL, assume that
    // all the variables at this node are fully restrained.
    int *bc_vars = new int[numBcs * 3 + 1];
    bc_ptr[0] = 0;

    for (int i = 0; i < numBcs; i++) {
      bc_ptr[i + 1] = bc_ptr[i];
      if (i == 0) {
        // Fix u, v, w, and rotx
        for (int j = 0; j < 4; j++) {
          bc_vars[bc_ptr[i + 1]] = j;
          bc_ptr[i + 1]++;
        }
      } else {
        // Fix v, w, and rotx
        for (int j = 0; j < 3; j++) {
          bc_vars[bc_ptr[i + 1]] = j + 1;
          bc_ptr[i + 1]++;
        }
      }
    }

    TacsScalar *bc_vals = new TacsScalar[bc_ptr[numBcs]];
    memset(bc_vals, 0, bc_ptr[numBcs] * sizeof(TacsScalar));

    // Set the boundary conditions
    creator->setBoundaryConditions(numBcs, bcNodes, bc_ptr, bc_vars, bc_vals);

    delete[] bcNodes;
    delete[] bc_ptr;
    delete[] bc_vars;
    delete[] bc_vals;

    // Set the node locations
    TacsScalar *Xpts = new TacsScalar[3 * numNodes];

    for (int j = 0; j < nny; j++) {
      double v = -M_PI + (2.0 * M_PI * j) / nny;
      for (int i = 0; i < nnx; i++) {
        double u = 1.0 * i / (nnx - 1);
        double theta =
            v + 0.25 * M_PI * u + defect * sin(v) * cos(2 * M_PI * u);
        double x = L * (u + defect * cos(v) * sin(2 * M_PI * u));

        int node = i + j * nnx;
        Xpts[3 * node] = x;
        Xpts[3 * node + 1] = R * cos(theta);
        Xpts[3 * node + 2] = -R * sin(theta);
      }
    }

    // Set the nodal locations
    creator->setNodes(Xpts);
    delete[] Xpts;
  }

  // Set the one element
  creator->setElements(1, &element);

  // Set the reordering type
  creator->setReorderingType(TACSAssembler::MULTICOLOR_ORDER,
                             TACSAssembler::GAUSS_SEIDEL);

  // Create TACS
  TACSAssembler *assembler = creator->createTACS();

  // Set the elements the node vector
  TACSBVec *X = assembler->createNodeVec();
  X->incref();
  assembler->getNodes(X);

  TACSAuxElements *aux = new TACSAuxElements();

  for (int elem = 0; elem < assembler->getNumElements(); elem++) {
    TacsScalar Xelem[3 * 9];
    TacsScalar tr[3 * 9];

    int nnodes;
    const int *nodes;
    assembler->getElement(elem, &nnodes, &nodes);
    X->getValues(nnodes, nodes, Xelem);

    for (int node = 0; node < nnodes; node++) {
      // Compute the pressure at this point in the shell
      double x = TacsRealPart(Xelem[3 * node]);
      double y = -R * atan2(TacsRealPart(Xelem[3 * node + 2]),
                            TacsRealPart(Xelem[3 * node + 1]));

      TacsScalar pval = load * sin(beta * x) * sin(alpha * y);
      TacsScalar ynorm = Xelem[3 * node + 1] / R;
      TacsScalar znorm = Xelem[3 * node + 2] / R;

      tr[3 * node] = 0.0;
      tr[3 * node + 1] = ynorm * pval;
      tr[3 * node + 2] = znorm * pval;
    }

    TACSElement *trac = NULL;
    if (order == 2) {
      trac = new TACSShellTraction<6, TACSTriLinearQuadrature,
                                   TACSShellTriLinearBasis>(tr);
    } else if (order == 3) {
      trac = new TACSShellTraction<6, TACSTriQuadraticQuadrature,
                                   TACSShellTriQuadraticBasis>(tr);
    }

    aux->addElement(elem, trac);
  }

  X->decref();

  // Set the auxiliary elements
  assembler->setAuxElements(aux);

  // Set the pointers
  *_assembler = assembler;
  *_creator = creator;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // Get the rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Parameters optionally set from the command line
  int order = 3;
  int nx = 20, ny = 40;
  int mesh_type = 0;

  // Parse the command line arguments
  for (int k = 0; k < argc; k++) {
    if (strcmp(argv[k], "triangle") == 0) {
      mesh_type = 1;
    }
    if (sscanf(argv[k], "nx=%d", &nx) == 1) {
      if (nx < 10) {
        nx = 10;
      }
      if (nx > 200) {
        nx = 200;
      }
    }
    if (sscanf(argv[k], "ny=%d", &ny) == 1) {
      if (ny < 10) {
        ny = 10;
      }
      if (ny > 200) {
        ny = 200;
      }
    }
    if (sscanf(argv[k], "order=%d", &order) == 1) {
      if (order < 2) {
        order = 2;
      }
      if (order > 3) {
        order = 3;
      }
    }
  }

  double load = 1.0;
  double t = 1.0;
  double L = 100.0;
  double R = 100.0 / M_PI;

  // Set the alpha and beta parameters
  double alpha = 4.0 / R;
  double beta = 3 * M_PI / L;

  TacsScalar rho = 2700.0;
  TacsScalar specific_heat = 921.096;
  TacsScalar E = 70e3;
  TacsScalar nu = 0.3;
  TacsScalar ys = 270.0;
  TacsScalar cte = 24.0e-6;
  TacsScalar kappa = 230.0;
  TACSMaterialProperties *props =
      new TACSMaterialProperties(rho, specific_heat, E, nu, ys, cte, kappa);

  TacsScalar axis[] = {1.0, 0.0, 0.0};
  TACSShellTransform *transform = new TACSShellRefAxisTransform(axis);

  TACSShellConstitutive *con = new TACSIsoShellConstitutive(props, t);

  TACSAssembler *assembler = NULL;
  TACSCreator *creator = NULL;
  TACSElement *shell = NULL;
  if (mesh_type == 0) {  // Quadrilateral mesh
    if (order == 2) {
      shell = new TACSQuad4Shell(transform, con);
    } else {  // order == 3
      shell = new TACSQuad9Shell(transform, con);
    }
    shell->incref();
    createAssembler(comm, order, nx, ny, shell, &assembler, &creator);
  } else {
    if (order == 2) {
      shell = new TACSTri3Shell(transform, con);
    } else {         // order == 3
      shell = NULL;  // This will cause a segfault
    }
    shell->incref();
    createTriAssembler(comm, order, nx, ny, shell, &assembler, &creator);
  }

  assembler->incref();
  creator->incref();

  // Free the creator object
  creator->decref();

  // Create matrix and vectors
  TACSBVec *ans = assembler->createVec();  // displacements and rotations
  TACSBVec *res = assembler->createVec();  // The residual
  TACSSchurMat *mat = assembler->createSchurMat();  // stiffness matrix

  // Increment reference count to the matrix/vectors
  ans->incref();
  res->incref();
  mat->incref();

  // Allocate the factorization
  int lev = 10000;
  double fill = 10.0;
  int reorder_schur = 1;
  TACSSchurPc *pc = new TACSSchurPc(mat, lev, fill, reorder_schur);
  pc->incref();

  assembler->assembleJacobian(1.0, 0.0, 0.0, res, mat);
  pc->factor();  // LU factorization of stiffness matrix
  pc->applyFactor(res, ans);

  ans->scale(-1.0);
  assembler->setVariables(ans);

  // Output for visualization
  ElementType etype = TACS_BEAM_OR_SHELL_ELEMENT;
  int write_flag = (TACS_OUTPUT_NODES | TACS_OUTPUT_CONNECTIVITY |
                    TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                    TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(assembler, etype, write_flag);
  f5->incref();
  f5->writeToFile("cylinder_solution.f5");

  // Compute the tangent stiffness
  TacsScalar Cs[TACSShellConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
  int eIndex = 0;
  double pt[] = {0.0, 0.0};
  TacsScalar Xpt[] = {0.0, 0.0, 0.0};
  con->evalTangentStiffness(eIndex, pt, Xpt, Cs);

  // Extract the A/B/D/As matrices
  TacsScalar drill;
  const TacsScalar *A, *B, *D, *As;
  con->extractTangentStiffness(Cs, &A, &B, &D, &As, &drill);

  TacsScalar U, V, W, theta, phi;
  TacsScalar ainv = 1.0 / R;
  computeCoefficients(&U, &V, &W, &theta, &phi, alpha, beta, ainv, A[0], A[1],
                      A[3], A[5], D[0], D[1], D[3], D[5], As[0], As[2], load);

  setExactSolution(assembler, ans, R, alpha, beta, U, V, W, theta, phi);

  // Set the force vector
  assembler->setVariables(ans);
  f5->writeToFile("cylinder_exact_solution.f5");

  shell->decref();
  assembler->decref();

  MPI_Finalize();
}
