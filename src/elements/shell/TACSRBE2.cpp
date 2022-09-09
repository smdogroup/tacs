#include "TACSRBE2.h"

/*
  A Nastran RBE2 rigid body element.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Copyright (C) 2020 Aerion Technologies Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  A rigid body connected to an arbitrary number of grid points.
  The independent degrees-of-freedom are the six components of motion at
  a single grid point. The dependent degrees-of-freedom at the other grid
  points. This element is implemented using a Lagrange multiplier method
  based on RBE2 constraints. Using this method, the entries of the
  stiffness matrix become the coefficients of the constraint matrix:

                 - -     -      -      -    -
                | f |   | 0  B^T |    |  u   |
                | 0 | = | B    0 | .  |lambda|
                 - -     -      -      -    -

  Where B is the constraint matrix (i.e. B . u = 0) and lambda are the
  Lagrange multipliers. The Lagrange multipliers
  represent the reaction forces on the each dependent grid due to the RBE2.

  For reference on this implementation see:

  [1]
  http://www2.me.rochester.edu/courses/ME204/nx_help/index.html#uid:id1322317

*/

/*
  Create the RBE2.
*/
TACSRBE2::TACSRBE2(int _numNodes, int _dof_constrained[], double _C1,
                   double _C2) {
  NUM_NODES = _numNodes;  // Number of nodes (1 indep node + N dep nodes + N
                          // dummy nodes)
  NUM_DEP_NODES = (NUM_NODES - 1) / 2;
  NUM_VARIABLES = NUM_DISPS * NUM_NODES;

  // Identify which dofs should be constrained for dependent nodes
  dof_constrained = new int *[NUM_DEP_NODES];
  for (int j = 0; j < NUM_DEP_NODES; j++) {
    dof_constrained[j] = new int[NUM_DISPS];
    for (int i = 0; i < NUM_DISPS; i++) {
      dof_constrained[j][i] = _dof_constrained[NUM_DISPS * j + i];
    }
  }

  // Default scaling and artificial stiffness parameters
  C1 = _C1;
  C2 = _C2;
}

TACSRBE2::~TACSRBE2() {
  if (dof_constrained) {
    for (int j = 0; j < NUM_DEP_NODES; j++) {
      delete[] dof_constrained[j];
    }
    delete[] dof_constrained;
  }
}

/*
  Retrieve information about the names of the element variables
*/
const char *TACSRBE2::getObjectName() { return elemName; }

const char *TACSRBE2::displacementName(int i) {
  if (i >= 0 && i < NUM_DISPS) {
    return dispNames[i];
  }
  return NULL;
}

const char *TACSRBE2::extraName(int i) {
  if (i >= 0 && i < NUM_DISPS) {
    return extraNames[i];
  }
  return NULL;
}

/*
  Retrieve the numbers of displacements, nodes, stress and variables
*/
int TACSRBE2::getVarsPerNode() { return NUM_DISPS; }

int TACSRBE2::getNumNodes() { return NUM_NODES; }

int TACSRBE2::numExtras() { return NUM_EXTRAS; }

ElementType TACSRBE2::getElementType() { return TACS_RIGID_ELEMENT; }

/*
  Returns the multiplier index
void TACSRBE2::getMultiplierIndex( int *multiplier ){
  *multiplier = NUM_DEP_NODES + 1;
}*/

/*
  The element name, variable, stress and strain names.
*/
const char *TACSRBE2::elemName = "TACSRBE2";

const char *TACSRBE2::dispNames[] = {"u0", "v0", "w0", "rotx", "roty", "rotz"};

const char *TACSRBE2::extraNames[] = {"fx", "fy", "fz", "mx", "my", "mz"};

/*
  Assemble the element residual associated with the given design
  variables and elements.
*/
void TACSRBE2::addResidual(int elemIndex, double time, const TacsScalar Xpts[],
                           const TacsScalar vars[], const TacsScalar dvars[],
                           const TacsScalar ddvars[], TacsScalar res[]) {
  const TacsScalar *Fn, *Mn, *Xn, *X0, *un, *tn, *u0, *t0, *actualLM;
  TacsScalar *maskedLM;

  // memset(res, 0, NUM_VARIABLES*sizeof(TacsScalar));

  // Store the last 6*N variables (Lagrange multipliers) as a force and moment
  int ii = NUM_DISPS * (1 + NUM_DEP_NODES);
  actualLM = &vars[ii + 0];
  maskedLM = new TacsScalar[NUM_DISPS * NUM_DEP_NODES];
  // Get reaction forces/moments at dep nodes from Lagrange multipliers,
  // zero out any terms that have been freed
  getMaskedMultipliers(maskedLM, actualLM);

  Fn = &maskedLM[0];
  Mn = &maskedLM[3];
  Xn = &Xpts[3];  // indep node starting point
  X0 = &Xpts[0];

  for (int n = 0; n < NUM_DEP_NODES; n++) {
    // The first 6 equations enforce equilibrium at the rbe2 indep node
    // caused by the forces on each of the dep nodes
    res[0] -= C1 * Fn[0];  // Fx
    res[1] -= C1 * Fn[1];  // Fy
    res[2] -= C1 * Fn[2];  // Fz

    // Moment equilibrium, we have to account for the couple portion as well as
    // the moment arm about the dep node
    res[3] -= C1 * (-(Xn[2] - X0[2]) * Fn[1] + (Xn[1] - X0[1]) * Fn[2]);  // Mx
    res[4] -= C1 * (+(Xn[2] - X0[2]) * Fn[0] - (Xn[0] - X0[0]) * Fn[2]);  // My
    res[5] -= C1 * (-(Xn[1] - X0[1]) * Fn[0] + (Xn[0] - X0[0]) * Fn[1]);  // Mz

    // Only include coupled moments if rotations are included
    res[3] -= C1 * (Mn[0]);  // Mx
    res[4] -= C1 * (Mn[1]);  // My
    res[5] -= C1 * (Mn[2]);  // Mz

    Fn += NUM_DISPS;
    Mn += NUM_DISPS;
    Xn += 3;
  }

  res += NUM_DISPS;

  // The next 6*NUM_DEP_NODES simply enforce that the forces on dep nodes match
  // whats stored in the Lagrange multipliers
  Fn = &maskedLM[0];
  Mn = &maskedLM[3];
  Xn = &Xpts[3];
  for (int n = 0; n < NUM_DEP_NODES; n++) {
    // Forces
    res[0] += C1 * Fn[0];  // Fx
    res[1] += C1 * Fn[1];  // Fy
    res[2] += C1 * Fn[2];  // Fz

    // Moments
    res[3] += C1 * Mn[0];  // Mx
    res[4] += C1 * Mn[1];  // My
    res[5] += C1 * Mn[2];  // Mz

    res += NUM_DISPS;
    Fn += NUM_DISPS;
    Mn += NUM_DISPS;
    Xn += 3;
  }

  // The last 6*NUM_DEP_NODES equations enforce that the dependent nodes
  // displacements are rigid transformations of the independent one (only add
  // equations for constrianed dofs)
  u0 = &vars[0];  // indep node displacements
  t0 = &vars[3];  // indep node rotations
  Xn = &Xpts[3];
  un = &vars[NUM_DISPS];      // dep nodes displacements
  tn = &vars[NUM_DISPS + 3];  // dep nodes rotations

  for (int n = 0; n < NUM_DEP_NODES; n++) {
    // displacements un = u0 + r x theta
    if (dof_constrained[n][0]) {
      res[0] += C1 * (un[0] - (u0[0] + (Xn[2] - X0[2]) * t0[1] -
                               (Xn[1] - X0[1]) * t0[2]));  // u
    }
    if (dof_constrained[n][1]) {
      res[1] += C1 * (un[1] - (u0[1] - (Xn[2] - X0[2]) * t0[0] +
                               (Xn[0] - X0[0]) * t0[2]));  // v
    }
    if (dof_constrained[n][2]) {
      res[2] += C1 * (un[2] - (u0[2] + (Xn[1] - X0[1]) * t0[0] -
                               (Xn[0] - X0[0]) * t0[1]));  // w
    }

    // rotations
    if (dof_constrained[n][3]) {
      res[3] += C1 * (tn[0] - t0[0]);  // thetax
    }
    if (dof_constrained[n][4]) {
      res[4] += C1 * (tn[1] - t0[1]);  // thetay
    }
    if (dof_constrained[n][5]) {
      res[5] += C1 * (tn[2] - t0[2]);  // thetaz
    }

    // For each unconstrained dof, set the residual to the
    // Lagrange multiplier forcing TACS to zero them
    for (int k = 0; k < NUM_DISPS; k++) {
      if (!dof_constrained[n][k]) {
        res[k] += C1 * actualLM[k];
      }
    }

    res += NUM_DISPS;
    un += NUM_DISPS;
    tn += NUM_DISPS;
    Mn += NUM_DISPS;
    Xn += 3;
    actualLM += NUM_DISPS;
  }

  delete[] maskedLM;
}

/*
  Assemble the stiffness matrix for the RBE2 element.
*/
void TACSRBE2::addJacobian(int elemIndex, double time, TacsScalar alpha,
                           TacsScalar beta, TacsScalar gamma,
                           const TacsScalar Xpts[], const TacsScalar vars[],
                           const TacsScalar dvars[], const TacsScalar ddvars[],
                           TacsScalar res[], TacsScalar J[]) {
  const TacsScalar *Xn, *X0;
  int col, row;
  TacsScalar *mat = new TacsScalar[NUM_VARIABLES * NUM_VARIABLES];
  memset(mat, 0, NUM_VARIABLES * NUM_VARIABLES * sizeof(TacsScalar));

  // Store the last 6*N variables (Lagrange multipliers) as a force and moment
  int ii = NUM_DISPS * (1 + NUM_DEP_NODES);
  Xn = &Xpts[3];  // indep node starting point
  X0 = &Xpts[0];

  // We'll take advantage of the matrix symmetry and fill both halves at once
  for (int n = 0; n < NUM_DEP_NODES; n++) {
    // The first 6 equations enforce equilibrium at the rbe2 indep node
    // caused by the forces on each of the dep nodes
    row = 0;
    col = ii + 0;
    mat[row + col * NUM_VARIABLES] = C1 * -1.0;
    mat[col + row * NUM_VARIABLES] = C1 * -1.0;

    row = 1;
    col = ii + 1;
    mat[row + col * NUM_VARIABLES] = C1 * -1.0;
    mat[col + row * NUM_VARIABLES] = C1 * -1.0;

    row = 2;
    col = ii + 2;
    mat[row + col * NUM_VARIABLES] = C1 * -1.0;
    mat[col + row * NUM_VARIABLES] = C1 * -1.0;

    // Moment equilibrium, we have to account for the couple portion as well as
    // the moment arm about the dep node Mx
    row = 3;
    col = ii + 1;
    mat[row + col * NUM_VARIABLES] = C1 * (Xn[2] - X0[2]);
    mat[col + row * NUM_VARIABLES] = C1 * (Xn[2] - X0[2]);
    col = ii + 2;
    mat[row + col * NUM_VARIABLES] = C1 * -(Xn[1] - X0[1]);
    mat[col + row * NUM_VARIABLES] = C1 * -(Xn[1] - X0[1]);

    // My
    row = 4;
    col = ii + 0;
    mat[row + col * NUM_VARIABLES] = C1 * -(Xn[2] - X0[2]);
    mat[col + row * NUM_VARIABLES] = C1 * -(Xn[2] - X0[2]);
    col = ii + 2;
    mat[row + col * NUM_VARIABLES] = C1 * (Xn[0] - X0[0]);
    mat[col + row * NUM_VARIABLES] = C1 * (Xn[0] - X0[0]);

    // Mz
    row = 5;
    col = ii + 0;
    mat[row + col * NUM_VARIABLES] = C1 * (Xn[1] - X0[1]);
    mat[col + row * NUM_VARIABLES] = C1 * (Xn[1] - X0[1]);
    col = ii + 1;
    mat[row + col * NUM_VARIABLES] = C1 * -(Xn[0] - X0[0]);
    mat[col + row * NUM_VARIABLES] = C1 * -(Xn[0] - X0[0]);

    for (row = 3; row < 6; row++) {
      col = ii + row;
      mat[row + col * NUM_VARIABLES] = C1 * -1.0;
      mat[col + row * NUM_VARIABLES] = C1 * -1.0;
    }

    ii += NUM_DISPS;
    Xn += 3;
  }

  row = NUM_DISPS;
  col = NUM_DISPS * (1 + NUM_DEP_NODES);

  // The next 6*NUM_DEP_NODES simply enforce that the forces on dep nodes match
  // whats stored in the Lagrange multipliers
  for (int n = 0; n < NUM_DEP_NODES; n++) {
    // Forces
    mat[row + col * NUM_VARIABLES] = C1 * 1.0;
    mat[col + row * NUM_VARIABLES] = C1 * 1.0;
    row++;
    col++;

    mat[row + col * NUM_VARIABLES] = C1 * 1.0;
    mat[col + row * NUM_VARIABLES] = C1 * 1.0;
    row++;
    col++;

    mat[row + col * NUM_VARIABLES] = C1 * 1.0;
    mat[col + row * NUM_VARIABLES] = C1 * 1.0;
    row++;
    col++;

    mat[row + col * NUM_VARIABLES] = C1 * 1.0;
    mat[col + row * NUM_VARIABLES] = C1 * 1.0;
    row++;
    col++;

    mat[row + col * NUM_VARIABLES] = C1 * 1.0;
    mat[col + row * NUM_VARIABLES] = C1 * 1.0;
    row++;
    col++;

    mat[row + col * NUM_VARIABLES] = C1 * 1.0;
    mat[col + row * NUM_VARIABLES] = C1 * 1.0;
    row++;
    col++;
  }

  // artificial stiffness for dep node and lagrange multiplier terms
  for (row = 0; row < NUM_NODES * NUM_DISPS; row++) {
    mat[row + row * NUM_VARIABLES] += C2 * 1.0;
  }

  /* NOTE: Lastly, we will loop through each unconstrained dof,
     replace the row and column in the matrix associated with that
     Lagrange multiplier with zero and set its corresponding diagonal
     term to unity. This should give us a consistent matrix with the
     procedure in getRes */
  for (int n = 0, ii = NUM_DISPS * (1 + NUM_DEP_NODES); n < NUM_DEP_NODES;
       n++, ii += NUM_DISPS) {
    for (int k = 0; k < NUM_DISPS; k++) {
      if (!dof_constrained[n][k]) {
        for (int i = 0; i < NUM_VARIABLES; i++) {
          mat[i + (ii + k) * NUM_VARIABLES] = 0.0;
          mat[ii + k + (i)*NUM_VARIABLES] = 0.0;
        }
        mat[ii + k + (ii + k) * NUM_VARIABLES] = C1 * 1.0;
      }
    }
  }

  for (int i = 0; i < NUM_VARIABLES * NUM_VARIABLES; i++) {
    J[i] += alpha * mat[i];
  }

  if (res) {
    addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);
  }

  delete[] mat;
}

/*
  Assemble the element residual associated with the given design
  variables and elements.
*/
void TACSRBE2::addAdjResXptProduct(
    int elemIndex, double time, TacsScalar scale, const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar fXptSens[]) {
  const TacsScalar *Fn, *t0, *actualLM;
  TacsScalar *sXpts, *sXn, *sX0, *maskedLM;

  TacsScalar *residual = new TacsScalar[NUM_NODES * 3 * NUM_VARIABLES];
  memset(residual, 0, NUM_NODES * 3 * NUM_VARIABLES * sizeof(TacsScalar));
  sXpts = new TacsScalar[(NUM_DEP_NODES + 1) * 3];

  // Store the last 6*NUM_DEP_NODES variables (Lagrange multipliers) as a force
  // and moment
  int ii = NUM_DISPS * (NUM_DEP_NODES + 1);
  actualLM = &vars[ii + 0];
  maskedLM = new TacsScalar[NUM_DISPS * NUM_DEP_NODES];
  // Get reaction forces/moments at dep nodes from Lagrange multipliers,
  // zero out any terms that have been freed
  getMaskedMultipliers(maskedLM, actualLM);

  TacsScalar *res = residual;

  for (int k = 0; k < 3 * (NUM_DEP_NODES + 1); k++) {
    memset(sXpts, 0, sizeof(TacsScalar) * (NUM_DEP_NODES + 1) * 3);
    sXpts[k] = 1.0;

    sX0 = &sXpts[0];
    sXn = &sXpts[3];
    Fn = &maskedLM[0];

    for (int n = 0; n < NUM_DEP_NODES; n++) {
      // The first 6 equations enforce equilibrium at the rbe2 indep node
      // caused by the forces on each of the dep nodes

      // Moment equilibrium, we have to account for the couple portion as well
      // as the moment arm about the dep node
      res[3] -=
          C1 * (-(sXn[2] - sX0[2]) * Fn[1] + (sXn[1] - sX0[1]) * Fn[2]);  // Mx
      res[4] -=
          C1 * (+(sXn[2] - sX0[2]) * Fn[0] - (sXn[0] - sX0[0]) * Fn[2]);  // My
      res[5] -=
          C1 * (-(sXn[1] - sX0[1]) * Fn[0] + (sXn[0] - sX0[0]) * Fn[1]);  // Mz

      Fn += NUM_DISPS;
      sXn += 3;
    }

    res += NUM_DISPS;

    // The next 6*NUM_DEP_NODES simply enforce that the forces on dep nodes
    // match whats stored in the Lagrange multipliers
    res += NUM_DISPS * NUM_DEP_NODES;

    // The last 6*NUM_DEP_NODES equations enforce that the dependent nodes
    // displacements are rigid transformations of the independent on
    t0 = &vars[3];  // indep node rotations
    sXn = &sXpts[3];

    for (int n = 0; n < NUM_DEP_NODES; n++) {
      // displacements un = u0 + r x theta
      if (dof_constrained[n][0]) {
        res[0] =
            C1 *
            (-(+(sXn[2] - sX0[2]) * t0[1] - (sXn[1] - sX0[1]) * t0[2]));  // u
      }
      if (dof_constrained[n][1]) {
        res[1] =
            C1 *
            (-(-(sXn[2] - sX0[2]) * t0[0] + (sXn[0] - sX0[0]) * t0[2]));  // v
      }
      if (dof_constrained[n][2]) {
        res[2] =
            C1 *
            (-(+(sXn[1] - sX0[1]) * t0[0] - (sXn[0] - sX0[0]) * t0[1]));  // w
      }

      res += NUM_DISPS;
      sXn += 3;
    }
  }

  // Reset residual array back to beginning
  res = residual;

  // Add sensitivity product to output vector
  for (int n = 0; n < 3 * NUM_NODES; n++) {
    for (int k = 0; k < NUM_VARIABLES; k++) {
      fXptSens[n] += psi[k] * scale * res[k];
    }
    res += NUM_VARIABLES;
  }

  delete[] sXpts;
  delete[] maskedLM;
  delete[] residual;
}

/* This procedure masks the Lagrange multipliers depending on
which dofs have been constrained. If the dof is constrained copy
the component, if not force it to zero. */
void TACSRBE2::getMaskedMultipliers(TacsScalar maskedLM[],
                                    const TacsScalar actualLM[]) {
  // Mask forces
  for (int j = 0; j < NUM_DEP_NODES; j++) {
    for (int i = 0; i < 6; i++) {
      if (dof_constrained[j][i]) {
        maskedLM[i] = actualLM[i];
      } else {
        maskedLM[i] = 0.0;
      }
    }
    maskedLM += NUM_DISPS;
    actualLM += NUM_DISPS;
  }
}

/*
  Retrieve the output data for this element
*/
/*void TACSRBE2::getOutputData( int elemIndex,
                                ElementType etype, int out_type,
                                const TacsScalar Xpts[],
                                const TacsScalar vars[],
                                const TacsScalar dvars[],
                                const TacsScalar ddvars[],
                                int ld_data, TacsScalar *data  ){
  if (etype == TACS_RIGID_ELEMENT){
      // Compute residual to find distributed force on each node
      TacsScalar* res;
      res = new TacsScalar[NUM_VARIABLES];
      memset(res, 0, NUM_VARIABLES*sizeof(TacsScalar));
      addResidual(0, 0.0, Xpts, vars, NULL, NULL, res);
      for ( int n = 0; n < NUM_DEP_NODES+1; n++ ){
        int index = 0;
        if (out_type & TACS_OUTPUT_NODES){
          for ( int k = 0; k < 3; k++ ){
        data[index+k] = TacsRealPart(Xpts[3*n+k]);
          }
          index += 3;
        }
        if (out_type & TACS_OUTPUT_DISPLACEMENTS){
          for ( int k = 0; k < NUM_DISPS; k++ ){
        data[index+k] = TacsRealPart(vars[NUM_DISPS*n+k]);
          }
          index += NUM_DISPS;
        }
        // Output RBE-applied forces on each node
        if (out_type & TACS_OUTPUT_EXTRAS){
          for ( int k = 0; k < NUM_DISPS; k++ ){
              if (n == 0){
                // flip force on indep node to get force applied by RBE
                data[index+k] = TacsRealPart(res[NUM_DISPS*n+k]);
              }
              else{
                data[index+k] = -TacsRealPart(res[NUM_DISPS*n+k]);
              }
          }
          index += NUM_DISPS;
        }

        data += ld_data;
      }
      delete [] res;
  }
}*/
