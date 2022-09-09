#include "TACSRBE3.h"

#include "TACSElementAlgebra.h"
/*
  An RBE3 rigid body element.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Copyright (C) 2020 Aerion Technologies Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  The RBE3 element is a powerful tool for distributing applied
  loads and mass in a model. Unlike the RBAR and RBE2 elements,
  the RBE3 doesn’t add additional stiffness to your structure.
  Forces and moments applied to reference points are distributed
  to a set of independent degrees of freedom based on the RBE3
  geometry and local weight factors. This element is implemented
  using a Lagrange multiplier method based on RBE3 constraints.
  Using this method, the entries of the stiffness matrix become
  the coefficients of the constraint matrix.

                 - -     -      -      -    -
                | f |   | K  B^T |    |  u   |
                | 0 | = | B    0 | .  |lambda|
                 - -     -      -      -    -

  Where B is the constraint matrix (i.e. B . u = 0) and lambda are the
  Lagrange multipliers. The Lagrange multipliers represent the reaction
  forces at the dependent node due to element.

  For reference on this implementation see:

  [1]
  https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/custom/advanced_nonlinear/advanced_nonlinear_tmg.pdf
  [2]
  http://www2.me.rochester.edu/courses/ME204/nx_help/index.html#uid:id1322317

*/

/*
  Create the RBE3.
*/
TACSRBE3::TACSRBE3(int numNodes, int _dep_dof_constrained[], double weights[],
                   int _indep_dof_constrained[], double _C1, double _C2) {
  NUM_NODES =
      numNodes;  // Number of nodes (1 dep node + N indep nodes + 1 dummy node)
  NUM_INDEP_NODES = NUM_NODES - 2;
  NUM_VARIABLES = NUM_DISPS * NUM_NODES;

  /* Setup DOF weights from user input */
  w = new double[NUM_INDEP_NODES];
  memcpy(w, weights, NUM_INDEP_NODES * sizeof(double));
  /* Identify which dofs should be constrained for dependent node */
  memcpy(dep_dof_constrained, _dep_dof_constrained, NUM_DISPS * sizeof(int));
  /* Identify which dofs should be constrained for independent nodes */
  indep_dof_constrained = new int *[NUM_INDEP_NODES];
  for (int j = 0; j < NUM_INDEP_NODES; j++) {
    indep_dof_constrained[j] = new int[NUM_DISPS];
    for (int i = 0; i < NUM_DISPS; i++) {
      indep_dof_constrained[j][i] = _indep_dof_constrained[NUM_DISPS * j + i];
    }
  }

  // Default scaling and artificial stiffness parameters
  C1 = _C1;
  C2 = _C2;
}

TACSRBE3::~TACSRBE3() {
  if (w) {
    delete[] w;
  }
  if (indep_dof_constrained) {
    for (int j = 0; j < NUM_INDEP_NODES; j++) {
      delete[] indep_dof_constrained[j];
    }
    delete[] indep_dof_constrained;
  }
}

/*
  Retrieve information about the names of the element variables
*/
const char *TACSRBE3::getObjectName() { return elemName; }

const char *TACSRBE3::displacementName(int i) {
  if (i >= 0 && i < NUM_DISPS) {
    return dispNames[i];
  }
  return NULL;
}

const char *TACSRBE3::extraName(int i) {
  if (i >= 0 && i < NUM_DISPS) {
    return extraNames[i];
  }
  return NULL;
}

/*
  Retrieve the numbers of displacements, nodes, stress and variables
*/
int TACSRBE3::getVarsPerNode() { return NUM_DISPS; }

int TACSRBE3::getNumNodes() { return NUM_NODES; }

int TACSRBE3::numExtras() { return NUM_EXTRAS; }

ElementType TACSRBE3::getElementType() { return TACS_RIGID_ELEMENT; }

/*
  Returns the multiplier index
*/
/*int TACSRBE3::getMultiplierIndex(){
  return NUM_INDEP_NODES + 1;
}*/

/*
  The element name, variable, stress and strain names.
*/
const char *TACSRBE3::elemName = "TACSRBE3";

/* Tolerance for colinearity test in moment of inertia calculation*/
const double TACSRBE3::SMALL_NUM = 1e-8;

const char *TACSRBE3::dispNames[] = {"u0", "v0", "w0", "rotx", "roty", "rotz"};

const char *TACSRBE3::extraNames[] = {"fx", "fy", "fz", "mx", "my", "mz"};

/*
  Assemble the element residual associated with the given design
  variables and elements.
*/
void TACSRBE3::addResidual(int elemIndex, double time, const TacsScalar Xpts[],
                           const TacsScalar vars[], const TacsScalar dvars[],
                           const TacsScalar ddvars[], TacsScalar residual[]) {
  TacsScalar Xcg[3], Jcg[3][3], W[3], Lc;
  const TacsScalar *X0, *Xn, *un, *F0, *M0, *u0, *t0, *tn;
  TacsScalar *maskedVars;

  TacsScalar *res = new TacsScalar[NUM_VARIABLES];
  memset(res, 0, NUM_VARIABLES * sizeof(TacsScalar));

  // Store the last 6 variables (Lagrange multipliers) as a force and moment
  int ii = NUM_DISPS * (NUM_NODES - 1);
  // Get reaction forces/moments at dep node from Lagrange multipliers,
  // zero out any terms that have been freed
  maskedVars = new TacsScalar[NUM_VARIABLES];
  getMaskedVars(maskedVars, vars);
  F0 = &maskedVars[ii + 0];
  M0 = &maskedVars[ii + 3];

  // Compute the centroid location and inverse moments of inertia tensor
  getCG(Xcg, W, w, Xpts);
  Lc = getMomentsOfInertia(Jcg, w, Xpts, Xcg);

  /* The first 6 equations enforce equilibrium between the
     reaction forces and applied forces at the dependent node.
     This is straightforward. */
  res[0] = F0[0];  // Fx
  res[1] = F0[1];  // Fy
  res[2] = F0[2];  // Fz

  res[3] = M0[0];  // Mx
  res[4] = M0[1];  // My
  res[5] = M0[2];  // Mz

  res += NUM_DISPS;

  /* The next 6*NUM_INDEP_NODES equations distribute the reaction forces
     from the dependent node to the independent nodes. To understand how
     this is done, let us first consider the forces and moments about the
     cg of the independent nodes. If we treat this collection of nodes as
     a rigid body (with each node having a mass of w[n]), then basic rigid
     body dynamic arguments can be used to distribute the loads from the cg.
     For instance, lets start with the distributed force component on the nodes
     due to the force applied on the cg. We know that the force will cause a
     uniform linear acceleration for all points in the body, given by the
     equation below: a = F_cg/W Where, W is the sum of all independent nodal
     weights (masses). The force on each node can be found by multiplying the
     acceleration by the nodal weight (mass): F_n = w[n]/W*F_cg The moment
     distribution can be treated in a similar way, in this case we know that a
     moment applied to the cg of a rigid body produces uniform rotational
     acceleration of all point about the cg: alpha = Icg^-1*M_cg Where Icg is
     the moment of inertia tensor for all nodes about the cg. To get the linear
     acceleration of any point due to the rotational acceleration we use the
     following equation: a_n = alpha x rcn = - rcn x Icg^-1 * M_cg Where rcn is
     the vector pointing from the cg to node n. Finally the force contribution
     can be found by multiplying by the weight (mass): F_n = - w[n] * rcn x
     Icg^-1 * M_cg Putting the force and moment contributions together gives
     full transfer equation: F_n = w[n]/W * F_cg - w[n] * rcn x Icg^-1 * M_cg
     The last step is to write the force and moment terms at the cg in terms of
     the forces we have (at the dependent node): F_n = w[n]/W * F_0 - w[n] * rcn
     x Icg^-1 * (M_0 + rco x F0)

     This formula can be generalized to include point moments using the approach
     defined in ref [1].

     This part is obviously more involved than the dependent node equilibrium...
     */
  X0 = &Xpts[0];  // dependent node location
  // vector pointing from cg to ref node
  TacsScalar rco[3], rnc[3];
  for (int i = 0; i < 3; i++) {
    rco[i] = X0[i] - Xcg[i];
  }
  Xn = &Xpts[3];  // independent node locations

  for (int n = 0; n < NUM_INDEP_NODES; n++) {
    // vector pointing from cg to current indep node
    for (int i = 0; i < 3; i++) {
      rnc[i] = Xn[i] - Xcg[i];
    }
    // Forces
    // Fx
    if (indep_dof_constrained[n][0]) {
      res[0] = -w[n] / W[0] * F0[0];
      res[0] += -w[n] * (Jcg[1][0] * rnc[2] - Jcg[2][0] * rnc[1]) *
                (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
      res[0] += -w[n] * (Jcg[1][1] * rnc[2] - Jcg[2][1] * rnc[1]) *
                (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
      res[0] += -w[n] * (Jcg[1][2] * rnc[2] - Jcg[2][2] * rnc[1]) *
                (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
    }

    // Fy
    if (indep_dof_constrained[n][1]) {
      res[1] = -w[n] / W[1] * F0[1];
      res[1] += -w[n] * (Jcg[2][0] * rnc[0] - Jcg[0][0] * rnc[2]) *
                (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
      res[1] += -w[n] * (Jcg[2][1] * rnc[0] - Jcg[0][1] * rnc[2]) *
                (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
      res[1] += -w[n] * (Jcg[2][2] * rnc[0] - Jcg[0][2] * rnc[2]) *
                (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
    }

    // Fz
    if (indep_dof_constrained[n][2]) {
      res[2] = -w[n] / W[2] * F0[2];
      res[2] += -w[n] * (Jcg[0][0] * rnc[1] - Jcg[1][0] * rnc[0]) *
                (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
      res[2] += -w[n] * (Jcg[0][1] * rnc[1] - Jcg[1][1] * rnc[0]) *
                (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
      res[2] += -w[n] * (Jcg[0][2] * rnc[1] - Jcg[1][2] * rnc[0]) *
                (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
    }

    // Moments
    // Mx
    if (indep_dof_constrained[n][3]) {
      res[3] = -w[n] * Lc * Lc * Jcg[0][0] *
               (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
      res[3] += -w[n] * Lc * Lc * Jcg[0][1] *
                (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
      res[3] += -w[n] * Lc * Lc * Jcg[0][2] *
                (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
    }

    // My
    if (indep_dof_constrained[n][4]) {
      res[4] = -w[n] * Lc * Lc * Jcg[1][0] *
               (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
      res[4] += -w[n] * Lc * Lc * Jcg[1][1] *
                (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
      res[4] += -w[n] * Lc * Lc * Jcg[1][2] *
                (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
    }

    if (indep_dof_constrained[n][5]) {
      res[5] = -w[n] * Lc * Lc * Jcg[2][0] *
               (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
      res[5] += -w[n] * Lc * Lc * Jcg[2][1] *
                (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
      res[5] += -w[n] * Lc * Lc * Jcg[2][2] *
                (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
    }

    res += NUM_DISPS;
    Xn += 3;
  }

  /* The last six equations enforce that the dependent nodes displacements are
     the weighted average of the independent nodes. I don't know that an
     intuitive interpretation of this equation exists, but we can take the
     transpose of the load transfer matrix above (using the Lagrange multiplier
     method) to figure out what the coefficients need to be. Only add equations
     for constrained DOFs */
  // start with dependent node (again, trivial)
  u0 = &vars[0];  // dep node displacements
  t0 = &vars[3];  // dep nodes rotations

  for (int i = 0; i < 3; i++) {
    if (dep_dof_constrained[i]) {
      // constrain displacement
      res[i] = u0[i];
    } else {
      // constrain Lagrange multiplier (to zero)
      res[i] = vars[ii + i];
    }

    if (dep_dof_constrained[i + 3]) {
      // constrain rotation
      res[i + 3] = t0[i];
    } else {
      // constrain Lagrange multiplier (to zero)
      res[i + 3] = vars[ii + 3 + i];
    }
  }

  /* Now the independent node contributions (not so trivial),
     just follow the logic from the previous section */
  Xn = &Xpts[3];
  un = &maskedVars[NUM_DISPS];
  tn = &maskedVars[NUM_DISPS + 3];

  for (int n = 0; n < NUM_INDEP_NODES; n++) {
    // vector pointing from cg to current indep node
    for (int i = 0; i < 3; i++) {
      rnc[i] = Xn[i] - Xcg[i];
    }

    // Displacements
    if (dep_dof_constrained[0]) {
      res[0] +=
          -w[n] / W[0] * un[0] +
          -w[n] * (Jcg[1][1] * rnc[2] - Jcg[2][1] * rnc[1]) * rco[2] * un[0] +
          -w[n] * (Jcg[1][2] * rnc[2] - Jcg[2][2] * rnc[1]) * -rco[1] * un[0] +
          -w[n] * (Jcg[2][1] * rnc[0] - Jcg[0][1] * rnc[2]) * rco[2] * un[1] +
          -w[n] * (Jcg[2][2] * rnc[0] - Jcg[0][2] * rnc[2]) * -rco[1] * un[1] +
          -w[n] * (Jcg[0][1] * rnc[1] - Jcg[1][1] * rnc[0]) * rco[2] * un[2] +
          -w[n] * (Jcg[0][2] * rnc[1] - Jcg[1][2] * rnc[0]) * -rco[1] * un[2] +
          -w[n] * Lc * Lc * Jcg[0][1] * rco[2] * tn[0] +
          -w[n] * Lc * Lc * Jcg[0][2] * -rco[1] * tn[0] +
          -w[n] * Lc * Lc * Jcg[1][1] * rco[2] * tn[1] +
          -w[n] * Lc * Lc * Jcg[1][2] * -rco[1] * tn[1] +
          -w[n] * Lc * Lc * Jcg[2][1] * rco[2] * tn[2] +
          -w[n] * Lc * Lc * Jcg[2][2] * -rco[1] * tn[2];  // u
    }

    if (dep_dof_constrained[1]) {
      res[1] +=
          -w[n] * (Jcg[1][0] * rnc[2] - Jcg[2][0] * rnc[1]) * -rco[2] * un[0] +
          -w[n] * (Jcg[1][2] * rnc[2] - Jcg[2][2] * rnc[1]) * rco[0] * un[0] +
          -w[n] / W[1] * un[1] +
          -w[n] * (Jcg[2][0] * rnc[0] - Jcg[0][0] * rnc[2]) * -rco[2] * un[1] +
          -w[n] * (Jcg[2][2] * rnc[0] - Jcg[0][2] * rnc[2]) * rco[0] * un[1] +
          -w[n] * (Jcg[0][0] * rnc[1] - Jcg[1][0] * rnc[0]) * -rco[2] * un[2] +
          -w[n] * (Jcg[0][2] * rnc[1] - Jcg[1][2] * rnc[0]) * rco[0] * un[2] +
          -w[n] * Lc * Lc * Jcg[0][0] * -rco[2] * tn[0] +
          -w[n] * Lc * Lc * Jcg[0][2] * rco[0] * tn[0] +
          -w[n] * Lc * Lc * Jcg[1][0] * -rco[2] * tn[1] +
          -w[n] * Lc * Lc * Jcg[1][2] * rco[0] * tn[1] +
          -w[n] * Lc * Lc * Jcg[2][0] * -rco[2] * tn[2] +
          -w[n] * Lc * Lc * Jcg[2][2] * rco[0] * tn[2];  // v
    }

    if (dep_dof_constrained[2]) {
      res[2] +=
          -w[n] * (Jcg[1][0] * rnc[2] - Jcg[2][0] * rnc[1]) * rco[1] * un[0] +
          -w[n] * (Jcg[1][1] * rnc[2] - Jcg[2][1] * rnc[1]) * -rco[0] * un[0] +
          -w[n] * (Jcg[2][0] * rnc[0] - Jcg[0][0] * rnc[2]) * rco[1] * un[1] +
          -w[n] * (Jcg[2][1] * rnc[0] - Jcg[0][1] * rnc[2]) * -rco[0] * un[1] +
          -w[n] / W[2] * un[2] +
          -w[n] * (Jcg[0][0] * rnc[1] - Jcg[1][0] * rnc[0]) * rco[1] * un[2] +
          -w[n] * (Jcg[0][1] * rnc[1] - Jcg[1][1] * rnc[0]) * -rco[0] * un[2] +
          -w[n] * Lc * Lc * Jcg[0][0] * rco[1] * tn[0] +
          -w[n] * Lc * Lc * Jcg[0][1] * -rco[0] * tn[0] +
          -w[n] * Lc * Lc * Jcg[1][0] * rco[1] * tn[1] +
          -w[n] * Lc * Lc * Jcg[1][1] * -rco[0] * tn[1] +
          -w[n] * Lc * Lc * Jcg[2][0] * rco[1] * tn[2] +
          -w[n] * Lc * Lc * Jcg[2][1] * -rco[0] * tn[2];  // w
    }
    // Rotations
    if (dep_dof_constrained[3]) {
      res[3] += -w[n] * (Jcg[1][0] * rnc[2] - Jcg[2][0] * rnc[1]) * un[0] +
                -w[n] * (Jcg[2][0] * rnc[0] - Jcg[0][0] * rnc[2]) * un[1] +
                -w[n] * (Jcg[0][0] * rnc[1] - Jcg[1][0] * rnc[0]) * un[2] +
                -w[n] * Lc * Lc * Jcg[0][0] * tn[0] +
                -w[n] * Lc * Lc * Jcg[1][0] * tn[1] +
                -w[n] * Lc * Lc * Jcg[2][0] * tn[2];  // thetax
    }

    if (dep_dof_constrained[4]) {
      res[4] += -w[n] * (Jcg[1][1] * rnc[2] - Jcg[2][1] * rnc[1]) * un[0] +
                -w[n] * (Jcg[2][1] * rnc[0] - Jcg[0][1] * rnc[2]) * un[1] +
                -w[n] * (Jcg[0][1] * rnc[1] - Jcg[1][1] * rnc[0]) * un[2] +
                -w[n] * Lc * Lc * Jcg[0][1] * tn[0] +
                -w[n] * Lc * Lc * Jcg[1][1] * tn[1] +
                -w[n] * Lc * Lc * Jcg[2][1] * tn[2];  // thetay
    }

    if (dep_dof_constrained[5]) {
      res[5] += -w[n] * (Jcg[1][2] * rnc[2] - Jcg[2][2] * rnc[1]) * un[0] +
                -w[n] * (Jcg[2][2] * rnc[0] - Jcg[0][2] * rnc[2]) * un[1] +
                -w[n] * (Jcg[0][2] * rnc[1] - Jcg[1][2] * rnc[0]) * un[2] +
                -w[n] * Lc * Lc * Jcg[0][2] * tn[0] +
                -w[n] * Lc * Lc * Jcg[1][2] * tn[1] +
                -w[n] * Lc * Lc * Jcg[2][2] * tn[2];  // thetaz
    }

    Xn += 3;
    un += NUM_DISPS;
    tn += NUM_DISPS;
  }

  // Step back to beginning of residual vector
  res -= NUM_VARIABLES - NUM_DISPS;
  // Rescale all equations so they are on the same order of magnitude as the
  // global stiffness matrix
  for (int i = 0; i < NUM_VARIABLES; i++) {
    residual[i] += C1 * res[i];
  }

  delete[] res;
  delete[] maskedVars;
}

/*
  Assemble the stiffness matrix for the RBE3 element.
*/
void TACSRBE3::addJacobian(int elemIndex, double time, TacsScalar alpha,
                           TacsScalar beta, TacsScalar gamma,
                           const TacsScalar Xpts[], const TacsScalar vars[],
                           const TacsScalar dvars[], const TacsScalar ddvars[],
                           TacsScalar res[], TacsScalar J[]) {
  TacsScalar Xcg[3], Jcg[3][3], W[3], Lc;
  const TacsScalar *Xn, *X0;
  int col, row;
  TacsScalar *mat = new TacsScalar[NUM_VARIABLES * NUM_VARIABLES];
  memset(mat, 0, NUM_VARIABLES * NUM_VARIABLES * sizeof(TacsScalar));

  // Store the last 6 variables (Lagrange multipliers) as a force and moment
  int ii = NUM_DISPS * (NUM_NODES - 1);

  // Compute the centroid location and moments of inertia
  getCG(Xcg, W, w, Xpts);
  Lc = getMomentsOfInertia(Jcg, w, Xpts, Xcg);

  // NOTE: For now we will populate the matrix as if the problem is fully
  // constrained,
  //       We will wait until the last step to zero out the necessary
  //       components.
  // We'll take advantage of the symmetry of

  /* The first 6 equations enforce equilibrium between the forces at the rbe3
   * dep node */

  // Fx
  row = 0;
  col = ii + 0;
  mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] = C1 * 1.0;

  // Fy
  row++;
  col = ii + 1;
  mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] = C1 * 1.0;

  // Fz
  row++;
  col = ii + 2;
  mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] = C1 * 1.0;

  // Mx
  row++;
  col = ii + 3;
  mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] = C1 * 1.0;

  row++;
  col = ii + 4;
  mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] = C1 * 1.0;

  row++;
  col = ii + 5;
  mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] = C1 * 1.0;

  // The next 6*NUM_INDEP_NODES equations distribute the forces at the dep node
  // among the dependent nodes
  X0 = &Xpts[0];  // dependent node location
  // vector pointing from cg to ref node
  TacsScalar rco[3], rnc[3];
  for (int i = 0; i < 3; i++) {
    rco[i] = X0[i] - Xcg[i];
  }
  Xn = &Xpts[3];  // independent node
  for (int n = 0; n < NUM_INDEP_NODES; n++) {
    // vector pointing from cg to current indep node
    for (int i = 0; i < 3; i++) {
      rnc[i] = Xn[i] - Xcg[i];
    }
    // Forces
    // Fx
    row++;
    col = ii + 0;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] / W[0] +
              -w[n] * (Jcg[1][1] * rnc[2] - Jcg[2][1] * rnc[1]) * rco[2] +
              -w[n] * (Jcg[1][2] * rnc[2] - Jcg[2][2] * rnc[1]) * -rco[1]);
    col = ii + 1;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[1][0] * rnc[2] - Jcg[2][0] * rnc[1]) * -rco[2] +
              -w[n] * (Jcg[1][2] * rnc[2] - Jcg[2][2] * rnc[1]) * rco[0]);
    col = ii + 2;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[1][0] * rnc[2] - Jcg[2][0] * rnc[1]) * rco[1] +
              -w[n] * (Jcg[1][1] * rnc[2] - Jcg[2][1] * rnc[1]) * -rco[0]);
    col = ii + 3;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[1][0] * rnc[2] - Jcg[2][0] * rnc[1]));
    col = ii + 4;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[1][1] * rnc[2] - Jcg[2][1] * rnc[1]));
    col = ii + 5;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[1][2] * rnc[2] - Jcg[2][2] * rnc[1]));

    // Fy
    row++;
    col = ii + 0;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[2][1] * rnc[0] - Jcg[0][1] * rnc[2]) * rco[2] +
              -w[n] * (Jcg[2][2] * rnc[0] - Jcg[0][2] * rnc[2]) * -rco[1]);
    col = ii + 1;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] / W[1] +
              -w[n] * (Jcg[2][0] * rnc[0] - Jcg[0][0] * rnc[2]) * -rco[2] +
              -w[n] * (Jcg[2][2] * rnc[0] - Jcg[0][2] * rnc[2]) * rco[0]);
    col = ii + 2;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[2][0] * rnc[0] - Jcg[0][0] * rnc[2]) * rco[1] +
              -w[n] * (Jcg[2][1] * rnc[0] - Jcg[0][1] * rnc[2]) * -rco[0]);
    col = ii + 3;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[2][0] * rnc[0] - Jcg[0][0] * rnc[2]));
    col = ii + 4;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[2][1] * rnc[0] - Jcg[0][1] * rnc[2]));
    col = ii + 5;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[2][2] * rnc[0] - Jcg[0][2] * rnc[2]));

    // Fz
    row++;
    col = ii + 0;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[0][1] * rnc[1] - Jcg[1][1] * rnc[0]) * rco[2] +
              -w[n] * (Jcg[0][2] * rnc[1] - Jcg[1][2] * rnc[0]) * -rco[1]);
    col = ii + 1;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[0][0] * rnc[1] - Jcg[1][0] * rnc[0]) * -rco[2] +
              -w[n] * (Jcg[0][2] * rnc[1] - Jcg[1][2] * rnc[0]) * rco[0]);
    col = ii + 2;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] / W[2] +
              -w[n] * (Jcg[0][0] * rnc[1] - Jcg[1][0] * rnc[0]) * rco[1] +
              -w[n] * (Jcg[0][1] * rnc[1] - Jcg[1][1] * rnc[0]) * -rco[0]);
    col = ii + 3;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[0][0] * rnc[1] - Jcg[1][0] * rnc[0]));
    col = ii + 4;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[0][1] * rnc[1] - Jcg[1][1] * rnc[0]));
    col = ii + 5;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * (-w[n] * (Jcg[0][2] * rnc[1] - Jcg[1][2] * rnc[0]));

    // Moments
    // Mx
    row++;
    col = ii + 0;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[0][1] * rco[2] +
        C1 * -w[n] * Lc * Lc * Jcg[0][2] * -rco[1];
    col = ii + 1;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[0][0] * -rco[2] +
        C1 * -w[n] * Lc * Lc * Jcg[0][2] * rco[0];
    col = ii + 2;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[0][0] * rco[1] +
        C1 * -w[n] * Lc * Lc * Jcg[0][1] * -rco[0];
    col = ii + 3;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[0][0];
    col = ii + 4;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[0][1];
    col = ii + 5;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[0][2];

    // My
    row++;
    col = ii + 0;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[1][1] * rco[2] +
        C1 * -w[n] * Lc * Lc * Jcg[1][2] * -rco[1];
    col = ii + 1;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[1][0] * -rco[2] +
        C1 * -w[n] * Lc * Lc * Jcg[1][2] * rco[0];
    col = ii + 2;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[1][0] * rco[1] +
        C1 * -w[n] * Lc * Lc * Jcg[1][1] * -rco[0];
    col = ii + 3;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[1][0];
    col = ii + 4;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[1][1];
    col = ii + 5;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[1][2];

    // Mz
    row++;
    col = ii + 0;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[2][1] * rco[2] +
        C1 * -w[n] * Lc * Lc * Jcg[2][2] * -rco[1];
    col = ii + 1;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[2][0] * -rco[2] +
        C1 * -w[n] * Lc * Lc * Jcg[2][2] * rco[0];
    col = ii + 2;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[2][0] * rco[1] +
        C1 * -w[n] * Lc * Lc * Jcg[2][1] * -rco[0];
    col = ii + 3;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[2][0];
    col = ii + 4;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[2][1];
    col = ii + 5;
    mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] =
        C1 * -w[n] * Lc * Lc * Jcg[2][2];

    Xn += 3;
  }

  /* This is not part of the traditional RBE3 formulation,
     but the RBE3 lagrange formulation tends to be singular.
     To prevent this, we add a small artificial stiffness
     to some of the diagonals. */
  for (int k = 0; k < NUM_DISPS; k++) {
    mat[k + k * NUM_VARIABLES] = C2 * 1.0;
    mat[ii + k + (ii + k) * NUM_VARIABLES] = C2 * 1.0;
  }

  /* NOTE: Lastly, we will loop through each unconstrained dof,
     replace the row and column in the matrix associated with that
     Lagrange multiplier with zero and set its corresponding diagonal
     term to unity. This should give us a consistent matrix with the
     procedure in getRes */
  for (int k = 0; k < NUM_DISPS; k++) {
    if (!dep_dof_constrained[k]) {
      for (int i = 0; i < NUM_VARIABLES; i++) {
        mat[i + (ii + k) * NUM_VARIABLES] = 0.0;
        mat[ii + k + (i)*NUM_VARIABLES] = 0.0;
      }
      mat[ii + k + (ii + k) * NUM_VARIABLES] = C1 * 1.0;
    }

    for (int j = 0; j < NUM_INDEP_NODES; j++) {
      if (!indep_dof_constrained[j][k]) {
        for (int i = 0; i < NUM_VARIABLES; i++) {
          int row = (j + 1) * NUM_DISPS + k;
          int col = i;
          mat[row + col * NUM_VARIABLES] = mat[col + row * NUM_VARIABLES] = 0.0;
        }
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
void TACSRBE3::addAdjResXptProduct(
    int elemIndex, double time, TacsScalar scale, const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar fXptSens[]) {
  TacsScalar Xcg[3], sXcg[3], Jcg[3][3], sJcg[3][3], W[3], Lc, sLc;
  const TacsScalar *X0, *Xn, *un, *F0, *M0, *u0, *t0, *tn;
  TacsScalar *sXpts, *sXn, *sX0;
  TacsScalar *maskedVars;

  TacsScalar *residual = new TacsScalar[NUM_NODES * 3 * NUM_VARIABLES];
  memset(residual, 0, NUM_NODES * 3 * NUM_VARIABLES * sizeof(TacsScalar));
  sXpts = new TacsScalar[(NUM_NODES - 1) * 3];

  // Store the last 6 variables (Lagrange multipliers) as a force and moment
  int ii = NUM_DISPS * (NUM_NODES - 1);
  // Get reaction forces/moments at dep node from Lagrange multipliers,
  // zero out any terms that have been freed
  maskedVars = new TacsScalar[NUM_VARIABLES];
  getMaskedVars(maskedVars, vars);
  F0 = &maskedVars[ii + 0];
  M0 = &maskedVars[ii + 3];

  TacsScalar *res = residual;

  for (int k = 0; k < 3 * (NUM_NODES - 1); k++) {
    memset(sXpts, 0, sizeof(TacsScalar) * (NUM_NODES - 1) * 3);
    sXpts[k] = 1.0;
    // Compute the centroid location and moments of inertia
    getCGSens(sXcg, Xcg, W, w, Xpts, k);
    Lc = getMomentsOfInertiaSens(sJcg, Jcg, &sLc, w, Xpts, Xcg, sXcg, k);

    /* The first 6 equations enforce equilibrium between the
       reaction forces and applied forces at the dependent node.
       This is straightforward. */
    res[0] = 0.0;  // Fx
    res[1] = 0.0;  // Fy
    res[2] = 0.0;  // Fz

    res[3] = 0.0;  // Mx
    res[4] = 0.0;  // My
    res[5] = 0.0;  // Mz

    res += NUM_DISPS;

    /* The next 6*NUM_INDEP_NODES equations distribute the reaction forces
       from the dependent node to the independent nodes. According to ref [1]
       for more info the force are distributed using the following equation:

       Fn = w_n/W * Fc + w_n/I * rnc x Mc

       This equation only applies for distributing loads from the cg to the
       independent nodes, this means we have to translate the nodes from the
       dependent node to the cg first, giving:

       Fn = w_n/W * F0 + w_n/I * rnc x (M0 + rco x F0)

       This part is obviously more involved than the dependent node
       equilibrium...
       */
    X0 = &Xpts[0];    // dependent node location
    sX0 = &sXpts[0];  // dependent node location
    // vector pointing from cg to ref node
    TacsScalar rco[3], rnc[3], srco[3], srnc[3];
    for (int i = 0; i < 3; i++) {
      rco[i] = X0[i] - Xcg[i];
      srco[i] = sX0[i] - sXcg[i];
    }
    Xn = &Xpts[3];    // independent node locations
    sXn = &sXpts[3];  // independent node locations

    for (int n = 0; n < NUM_INDEP_NODES; n++) {
      // vector pointing from cg to current indep node
      for (int i = 0; i < 3; i++) {
        rnc[i] = Xn[i] - Xcg[i];
        srnc[i] = sXn[i] - sXcg[i];
      }
      // Forces
      // Fx
      if (indep_dof_constrained[n][0]) {
        res[0] = 0.0;
        res[0] += -w[n] * (sJcg[1][0] * rnc[2] - sJcg[2][0] * rnc[1]) *
                  (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
        res[0] += -w[n] * (sJcg[1][1] * rnc[2] - sJcg[2][1] * rnc[1]) *
                  (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
        res[0] += -w[n] * (sJcg[1][2] * rnc[2] - sJcg[2][2] * rnc[1]) *
                  (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
        res[0] += -w[n] * (Jcg[1][0] * srnc[2] - Jcg[2][0] * srnc[1]) *
                  (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
        res[0] += -w[n] * (Jcg[1][1] * srnc[2] - Jcg[2][1] * srnc[1]) *
                  (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
        res[0] += -w[n] * (Jcg[1][2] * srnc[2] - Jcg[2][2] * srnc[1]) *
                  (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
        res[0] += -w[n] * (Jcg[1][0] * rnc[2] - Jcg[2][0] * rnc[1]) *
                  (srco[1] * F0[2] - srco[2] * F0[1]);
        res[0] += -w[n] * (Jcg[1][1] * rnc[2] - Jcg[2][1] * rnc[1]) *
                  (-srco[0] * F0[2] + srco[2] * F0[0]);
        res[0] += -w[n] * (Jcg[1][2] * rnc[2] - Jcg[2][2] * rnc[1]) *
                  (srco[0] * F0[1] - srco[1] * F0[0]);
      }

      // Fy
      if (indep_dof_constrained[n][1]) {
        res[1] = 0.0;
        res[1] += -w[n] * (sJcg[2][0] * rnc[0] - sJcg[0][0] * rnc[2]) *
                  (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
        res[1] += -w[n] * (sJcg[2][1] * rnc[0] - sJcg[0][1] * rnc[2]) *
                  (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
        res[1] += -w[n] * (sJcg[2][2] * rnc[0] - sJcg[0][2] * rnc[2]) *
                  (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
        res[1] += -w[n] * (Jcg[2][0] * srnc[0] - Jcg[0][0] * srnc[2]) *
                  (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
        res[1] += -w[n] * (Jcg[2][1] * srnc[0] - Jcg[0][1] * srnc[2]) *
                  (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
        res[1] += -w[n] * (Jcg[2][2] * srnc[0] - Jcg[0][2] * srnc[2]) *
                  (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
        res[1] += -w[n] * (Jcg[2][0] * rnc[0] - Jcg[0][0] * rnc[2]) *
                  (srco[1] * F0[2] - srco[2] * F0[1]);
        res[1] += -w[n] * (Jcg[2][1] * rnc[0] - Jcg[0][1] * rnc[2]) *
                  (-srco[0] * F0[2] + srco[2] * F0[0]);
        res[1] += -w[n] * (Jcg[2][2] * rnc[0] - Jcg[0][2] * rnc[2]) *
                  (srco[0] * F0[1] - srco[1] * F0[0]);
      }

      // Fz
      if (indep_dof_constrained[n][2]) {
        res[2] = 0.0;
        res[2] += -w[n] * (sJcg[0][0] * rnc[1] - sJcg[1][0] * rnc[0]) *
                  (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
        res[2] += -w[n] * (sJcg[0][1] * rnc[1] - sJcg[1][1] * rnc[0]) *
                  (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
        res[2] += -w[n] * (sJcg[0][2] * rnc[1] - sJcg[1][2] * rnc[0]) *
                  (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
        res[2] += -w[n] * (Jcg[0][0] * srnc[1] - Jcg[1][0] * srnc[0]) *
                  (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
        res[2] += -w[n] * (Jcg[0][1] * srnc[1] - Jcg[1][1] * srnc[0]) *
                  (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
        res[2] += -w[n] * (Jcg[0][2] * srnc[1] - Jcg[1][2] * srnc[0]) *
                  (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
        res[2] += -w[n] * (Jcg[0][0] * rnc[1] - Jcg[1][0] * rnc[0]) *
                  (srco[1] * F0[2] - srco[2] * F0[1]);
        res[2] += -w[n] * (Jcg[0][1] * rnc[1] - Jcg[1][1] * rnc[0]) *
                  (-srco[0] * F0[2] + srco[2] * F0[0]);
        res[2] += -w[n] * (Jcg[0][2] * rnc[1] - Jcg[1][2] * rnc[0]) *
                  (srco[0] * F0[1] - srco[1] * F0[0]);
      }

      // Moments
      // Mx
      if (indep_dof_constrained[n][3]) {
        res[3] = 2.0 * -w[n] * Lc * sLc * Jcg[0][0] *
                 (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
        res[3] += -w[n] * Lc * Lc * sJcg[0][0] *
                  (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
        res[3] +=
            -w[n] * Lc * Lc * Jcg[0][0] * (srco[1] * F0[2] - srco[2] * F0[1]);
        res[3] += 2.0 * -w[n] * Lc * sLc * Jcg[0][1] *
                  (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
        res[3] += -w[n] * Lc * Lc * sJcg[0][1] *
                  (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
        res[3] +=
            -w[n] * Lc * Lc * Jcg[0][1] * (-srco[0] * F0[2] + srco[2] * F0[0]);
        res[3] += 2.0 * -w[n] * Lc * sLc * Jcg[0][2] *
                  (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
        res[3] += -w[n] * Lc * Lc * sJcg[0][2] *
                  (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
        res[3] +=
            -w[n] * Lc * Lc * Jcg[0][2] * (srco[0] * F0[1] - srco[1] * F0[0]);
      }

      // My
      if (indep_dof_constrained[n][4]) {
        res[4] = 2.0 * -w[n] * Lc * sLc * Jcg[1][0] *
                 (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
        res[4] += -w[n] * Lc * Lc * sJcg[1][0] *
                  (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
        res[4] +=
            -w[n] * Lc * Lc * Jcg[1][0] * (srco[1] * F0[2] - srco[2] * F0[1]);
        res[4] += 2.0 * -w[n] * Lc * sLc * Jcg[1][1] *
                  (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
        res[4] += -w[n] * Lc * Lc * sJcg[1][1] *
                  (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
        res[4] +=
            -w[n] * Lc * Lc * Jcg[1][1] * (-srco[0] * F0[2] + srco[2] * F0[0]);
        res[4] += 2.0 * -w[n] * Lc * sLc * Jcg[1][2] *
                  (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
        res[4] += -w[n] * Lc * Lc * sJcg[1][2] *
                  (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
        res[4] +=
            -w[n] * Lc * Lc * Jcg[1][2] * (srco[0] * F0[1] - srco[1] * F0[0]);
      }

      if (indep_dof_constrained[n][5]) {
        res[5] = 2.0 * -w[n] * Lc * sLc * Jcg[2][0] *
                 (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
        res[5] += -w[n] * Lc * Lc * sJcg[2][0] *
                  (M0[0] + rco[1] * F0[2] - rco[2] * F0[1]);
        res[5] +=
            -w[n] * Lc * Lc * Jcg[2][0] * (srco[1] * F0[2] - srco[2] * F0[1]);
        res[5] += 2.0 * -w[n] * Lc * sLc * Jcg[2][1] *
                  (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
        res[5] += -w[n] * Lc * Lc * sJcg[2][1] *
                  (M0[1] - rco[0] * F0[2] + rco[2] * F0[0]);
        res[5] +=
            -w[n] * Lc * Lc * Jcg[2][1] * (-srco[0] * F0[2] + srco[2] * F0[0]);
        res[5] += 2.0 * -w[n] * Lc * sLc * Jcg[2][2] *
                  (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
        res[5] += -w[n] * Lc * Lc * sJcg[2][2] *
                  (M0[2] + rco[0] * F0[1] - rco[1] * F0[0]);
        res[5] +=
            -w[n] * Lc * Lc * Jcg[2][2] * (srco[0] * F0[1] - srco[1] * F0[0]);
      }

      res += NUM_DISPS;
      Xn += 3;
      sXn += 3;
    }

    /* The last six equations enforce that the dependent nodes displacements are
       the weighted average of the independent nodes. I don't know that an
       intuitive interpretation of this equation exists, but we can take the
       transpose of the load transfer matrix above (using the Lagrange
       multiplier method) to figure out what the coefficients need to be. */
    // start with dependent node (again, trivial)
    u0 = &vars[0];  // dep node displacements
    t0 = &vars[3];  // dep nodes rotations

    // displacements
    res[0] = 0.0;  // u
    res[1] = 0.0;  // v
    res[2] = 0.0;  // w

    // rotations
    res[3] = 0.0;  // thetax
    res[4] = 0.0;  // thetay
    res[5] = 0.0;  // thetaz

    /* Now the independent node contributions (not so trivial),
       just follow the logic from the previous section */
    Xn = &Xpts[3];
    sXn = &sXpts[3];
    un = &maskedVars[NUM_DISPS];
    tn = &maskedVars[NUM_DISPS + 3];

    for (int n = 0; n < NUM_INDEP_NODES; n++) {
      // vector pointing from cg to current indep node
      for (int i = 0; i < 3; i++) {
        rnc[i] = Xn[i] - Xcg[i];
        srnc[i] = sXn[i] - sXcg[i];
      }

      // Displacements
      if (dep_dof_constrained[0]) {
        res[0] += 0.0 +
                  -w[n] * (sJcg[1][1] * rnc[2] - sJcg[2][1] * rnc[1]) * rco[2] *
                      un[0] +
                  -w[n] * (sJcg[1][2] * rnc[2] - sJcg[2][2] * rnc[1]) *
                      -rco[1] * un[0] +
                  -w[n] * (sJcg[2][1] * rnc[0] - sJcg[0][1] * rnc[2]) * rco[2] *
                      un[1] +
                  -w[n] * (sJcg[2][2] * rnc[0] - sJcg[0][2] * rnc[2]) *
                      -rco[1] * un[1] +
                  -w[n] * (sJcg[0][1] * rnc[1] - sJcg[1][1] * rnc[0]) * rco[2] *
                      un[2] +
                  -w[n] * (sJcg[0][2] * rnc[1] - sJcg[1][2] * rnc[0]) *
                      -rco[1] * un[2] +
                  -w[n] * (Jcg[1][1] * srnc[2] - Jcg[2][1] * srnc[1]) * rco[2] *
                      un[0] +
                  -w[n] * (Jcg[1][2] * srnc[2] - Jcg[2][2] * srnc[1]) *
                      -rco[1] * un[0] +
                  -w[n] * (Jcg[2][1] * srnc[0] - Jcg[0][1] * srnc[2]) * rco[2] *
                      un[1] +
                  -w[n] * (Jcg[2][2] * srnc[0] - Jcg[0][2] * srnc[2]) *
                      -rco[1] * un[1] +
                  -w[n] * (Jcg[0][1] * srnc[1] - Jcg[1][1] * srnc[0]) * rco[2] *
                      un[2] +
                  -w[n] * (Jcg[0][2] * srnc[1] - Jcg[1][2] * srnc[0]) *
                      -rco[1] * un[2] +
                  -w[n] * (Jcg[1][1] * rnc[2] - Jcg[2][1] * rnc[1]) * srco[2] *
                      un[0] +
                  -w[n] * (Jcg[1][2] * rnc[2] - Jcg[2][2] * rnc[1]) * -srco[1] *
                      un[0] +
                  -w[n] * (Jcg[2][1] * rnc[0] - Jcg[0][1] * rnc[2]) * srco[2] *
                      un[1] +
                  -w[n] * (Jcg[2][2] * rnc[0] - Jcg[0][2] * rnc[2]) * -srco[1] *
                      un[1] +
                  -w[n] * (Jcg[0][1] * rnc[1] - Jcg[1][1] * rnc[0]) * srco[2] *
                      un[2] +
                  -w[n] * (Jcg[0][2] * rnc[1] - Jcg[1][2] * rnc[0]) * -srco[1] *
                      un[2] +
                  2.0 * -w[n] * Lc * sLc * Jcg[0][1] * rco[2] * tn[0] +
                  -w[n] * Lc * Lc * sJcg[0][1] * rco[2] * tn[0] +
                  -w[n] * Lc * Lc * Jcg[0][1] * srco[2] * tn[0] +
                  2.0 * -w[n] * Lc * sLc * Jcg[0][2] * -rco[1] * tn[0] +
                  -w[n] * Lc * Lc * sJcg[0][2] * -rco[1] * tn[0] +
                  -w[n] * Lc * Lc * Jcg[0][2] * -srco[1] * tn[0] +
                  2.0 * -w[n] * Lc * sLc * Jcg[1][1] * rco[2] * tn[1] +
                  -w[n] * Lc * Lc * sJcg[1][1] * rco[2] * tn[1] +
                  -w[n] * Lc * Lc * Jcg[1][1] * srco[2] * tn[1] +
                  2.0 * -w[n] * Lc * sLc * Jcg[1][2] * -rco[1] * tn[1] +
                  -w[n] * Lc * Lc * sJcg[1][2] * -rco[1] * tn[1] +
                  -w[n] * Lc * Lc * Jcg[1][2] * -srco[1] * tn[1] +
                  2.0 * -w[n] * Lc * sLc * Jcg[2][1] * rco[2] * tn[2] +
                  -w[n] * Lc * Lc * sJcg[2][1] * rco[2] * tn[2] +
                  -w[n] * Lc * Lc * Jcg[2][1] * srco[2] * tn[2] +
                  2.0 * -w[n] * Lc * sLc * Jcg[2][2] * -rco[1] * tn[2] +
                  -w[n] * Lc * Lc * sJcg[2][2] * -rco[1] * tn[2] +
                  -w[n] * Lc * Lc * Jcg[2][2] * -srco[1] * tn[2];  // u
      }

      if (dep_dof_constrained[1]) {
        res[1] += -w[n] * (sJcg[1][0] * rnc[2] - sJcg[2][0] * rnc[1]) *
                      -rco[2] * un[0] +
                  -w[n] * (sJcg[1][2] * rnc[2] - sJcg[2][2] * rnc[1]) * rco[0] *
                      un[0] +
                  0.0 +
                  -w[n] * (sJcg[2][0] * rnc[0] - sJcg[0][0] * rnc[2]) *
                      -rco[2] * un[1] +
                  -w[n] * (sJcg[2][2] * rnc[0] - sJcg[0][2] * rnc[2]) * rco[0] *
                      un[1] +
                  -w[n] * (sJcg[0][0] * rnc[1] - sJcg[1][0] * rnc[0]) *
                      -rco[2] * un[2] +
                  -w[n] * (sJcg[0][2] * rnc[1] - sJcg[1][2] * rnc[0]) * rco[0] *
                      un[2] -
                  w[n] * (Jcg[1][0] * srnc[2] - Jcg[2][0] * srnc[1]) * -rco[2] *
                      un[0] +
                  -w[n] * (Jcg[1][2] * srnc[2] - Jcg[2][2] * srnc[1]) * rco[0] *
                      un[0] +
                  -w[n] * (Jcg[2][0] * srnc[0] - Jcg[0][0] * srnc[2]) *
                      -rco[2] * un[1] +
                  -w[n] * (Jcg[2][2] * srnc[0] - Jcg[0][2] * srnc[2]) * rco[0] *
                      un[1] +
                  -w[n] * (Jcg[0][0] * srnc[1] - Jcg[1][0] * srnc[0]) *
                      -rco[2] * un[2] +
                  -w[n] * (Jcg[0][2] * srnc[1] - Jcg[1][2] * srnc[0]) * rco[0] *
                      un[2] -
                  w[n] * (Jcg[1][0] * rnc[2] - Jcg[2][0] * rnc[1]) * -srco[2] *
                      un[0] +
                  -w[n] * (Jcg[1][2] * rnc[2] - Jcg[2][2] * rnc[1]) * srco[0] *
                      un[0] +
                  -w[n] * (Jcg[2][0] * rnc[0] - Jcg[0][0] * rnc[2]) * -srco[2] *
                      un[1] +
                  -w[n] * (Jcg[2][2] * rnc[0] - Jcg[0][2] * rnc[2]) * srco[0] *
                      un[1] +
                  -w[n] * (Jcg[0][0] * rnc[1] - Jcg[1][0] * rnc[0]) * -srco[2] *
                      un[2] +
                  -w[n] * (Jcg[0][2] * rnc[1] - Jcg[1][2] * rnc[0]) * srco[0] *
                      un[2] +
                  2.0 * -w[n] * Lc * sLc * Jcg[0][0] * -rco[2] * tn[0] +
                  -w[n] * Lc * Lc * sJcg[0][0] * -rco[2] * tn[0] +
                  -w[n] * Lc * Lc * Jcg[0][0] * -srco[2] * tn[0] +
                  2.0 * -w[n] * Lc * sLc * Jcg[0][2] * rco[0] * tn[0] +
                  -w[n] * Lc * Lc * sJcg[0][2] * rco[0] * tn[0] +
                  -w[n] * Lc * Lc * Jcg[0][2] * srco[0] * tn[0] +
                  2.0 * -w[n] * Lc * sLc * Jcg[1][0] * -rco[2] * tn[1] +
                  -w[n] * Lc * Lc * sJcg[1][0] * -rco[2] * tn[1] +
                  -w[n] * Lc * Lc * Jcg[1][0] * -srco[2] * tn[1] +
                  2.0 * -w[n] * Lc * sLc * Jcg[1][2] * rco[0] * tn[1] +
                  -w[n] * Lc * Lc * sJcg[1][2] * rco[0] * tn[1] +
                  -w[n] * Lc * Lc * Jcg[1][2] * srco[0] * tn[1] +
                  2.0 * -w[n] * Lc * sLc * Jcg[2][0] * -rco[2] * tn[2] +
                  -w[n] * Lc * Lc * sJcg[2][0] * -rco[2] * tn[2] +
                  -w[n] * Lc * Lc * Jcg[2][0] * -srco[2] * tn[2] +
                  2.0 * -w[n] * Lc * sLc * Jcg[2][2] * rco[0] * tn[2] +
                  -w[n] * Lc * Lc * sJcg[2][2] * rco[0] * tn[2] +
                  -w[n] * Lc * Lc * Jcg[2][2] * srco[0] * tn[2];  // v
      }

      if (dep_dof_constrained[2]) {
        res[2] +=
            -w[n] * (sJcg[1][0] * rnc[2] - sJcg[2][0] * rnc[1]) * rco[1] *
                un[0] +
            -w[n] * (sJcg[1][1] * rnc[2] - sJcg[2][1] * rnc[1]) * -rco[0] *
                un[0] +
            -w[n] * (sJcg[2][0] * rnc[0] - sJcg[0][0] * rnc[2]) * rco[1] *
                un[1] +
            -w[n] * (sJcg[2][1] * rnc[0] - sJcg[0][1] * rnc[2]) * -rco[0] *
                un[1] +
            0.0 +
            -w[n] * (sJcg[0][0] * rnc[1] - sJcg[1][0] * rnc[0]) * rco[1] *
                un[2] +
            -w[n] * (sJcg[0][1] * rnc[1] - sJcg[1][1] * rnc[0]) * -rco[0] *
                un[2] -
            w[n] * (Jcg[1][0] * srnc[2] - Jcg[2][0] * srnc[1]) * rco[1] *
                un[0] +
            -w[n] * (Jcg[1][1] * srnc[2] - Jcg[2][1] * srnc[1]) * -rco[0] *
                un[0] +
            -w[n] * (Jcg[2][0] * srnc[0] - Jcg[0][0] * srnc[2]) * rco[1] *
                un[1] +
            -w[n] * (Jcg[2][1] * srnc[0] - Jcg[0][1] * srnc[2]) * -rco[0] *
                un[1] +
            -w[n] * (Jcg[0][0] * srnc[1] - Jcg[1][0] * srnc[0]) * rco[1] *
                un[2] +
            -w[n] * (Jcg[0][1] * srnc[1] - Jcg[1][1] * srnc[0]) * -rco[0] *
                un[2] -
            w[n] * (Jcg[1][0] * rnc[2] - Jcg[2][0] * rnc[1]) * srco[1] * un[0] +
            -w[n] * (Jcg[1][1] * rnc[2] - Jcg[2][1] * rnc[1]) * -srco[0] *
                un[0] +
            -w[n] * (Jcg[2][0] * rnc[0] - Jcg[0][0] * rnc[2]) * srco[1] *
                un[1] +
            -w[n] * (Jcg[2][1] * rnc[0] - Jcg[0][1] * rnc[2]) * -srco[0] *
                un[1] +
            -w[n] * (Jcg[0][0] * rnc[1] - Jcg[1][0] * rnc[0]) * srco[1] *
                un[2] +
            -w[n] * (Jcg[0][1] * rnc[1] - Jcg[1][1] * rnc[0]) * -srco[0] *
                un[2] +
            2.0 * -w[n] * Lc * sLc * Jcg[0][0] * rco[1] * tn[0] +
            -w[n] * Lc * Lc * sJcg[0][0] * rco[1] * tn[0] +
            -w[n] * Lc * Lc * Jcg[0][0] * srco[1] * tn[0] +
            2.0 * -w[n] * Lc * sLc * Jcg[0][1] * -rco[0] * tn[0] +
            -w[n] * Lc * Lc * sJcg[0][1] * -rco[0] * tn[0] +
            -w[n] * Lc * Lc * Jcg[0][1] * -srco[0] * tn[0] +
            2.0 * -w[n] * Lc * sLc * Jcg[1][0] * rco[1] * tn[1] +
            -w[n] * Lc * Lc * sJcg[1][0] * rco[1] * tn[1] +
            -w[n] * Lc * Lc * Jcg[1][0] * srco[1] * tn[1] +
            2.0 * -w[n] * Lc * sLc * Jcg[1][1] * -rco[0] * tn[1] +
            -w[n] * Lc * Lc * sJcg[1][1] * -rco[0] * tn[1] +
            -w[n] * Lc * Lc * Jcg[1][1] * -srco[0] * tn[1] +
            2.0 * -w[n] * Lc * sLc * Jcg[2][0] * rco[1] * tn[2] +
            -w[n] * Lc * Lc * sJcg[2][0] * rco[1] * tn[2] +
            -w[n] * Lc * Lc * Jcg[2][0] * srco[1] * tn[2] +
            2.0 * -w[n] * Lc * sLc * Jcg[2][1] * -rco[0] * tn[2] +
            -w[n] * Lc * Lc * sJcg[2][1] * -rco[0] * tn[2] +
            -w[n] * Lc * Lc * Jcg[2][1] * -srco[0] * tn[2];  // w; // w
      }

      // Rotations
      if (dep_dof_constrained[3]) {
        res[3] += -w[n] * (sJcg[1][0] * rnc[2] - sJcg[2][0] * rnc[1]) * un[0] +
                  -w[n] * (sJcg[2][0] * rnc[0] - sJcg[0][0] * rnc[2]) * un[1] +
                  -w[n] * (sJcg[0][0] * rnc[1] - sJcg[1][0] * rnc[0]) * un[2] -
                  w[n] * (Jcg[1][0] * srnc[2] - Jcg[2][0] * srnc[1]) * un[0] +
                  -w[n] * (Jcg[2][0] * srnc[0] - Jcg[0][0] * srnc[2]) * un[1] +
                  -w[n] * (Jcg[0][0] * srnc[1] - Jcg[1][0] * srnc[0]) * un[2] +
                  2.0 * -w[n] * Lc * sLc * Jcg[0][0] * tn[0] +
                  -w[n] * Lc * Lc * sJcg[0][0] * tn[0] +
                  2.0 * -w[n] * Lc * sLc * Jcg[1][0] * tn[1] +
                  -w[n] * Lc * Lc * sJcg[1][0] * tn[1] +
                  2.0 * -w[n] * Lc * sLc * Jcg[2][0] * tn[2] +
                  -w[n] * Lc * Lc * sJcg[2][0] * tn[2];  // thetax
      }

      if (dep_dof_constrained[4]) {
        res[4] += -w[n] * (sJcg[1][1] * rnc[2] - sJcg[2][1] * rnc[1]) * un[0] +
                  -w[n] * (sJcg[2][1] * rnc[0] - sJcg[0][1] * rnc[2]) * un[1] +
                  -w[n] * (sJcg[0][1] * rnc[1] - sJcg[1][1] * rnc[0]) * un[2] -
                  w[n] * (Jcg[1][1] * srnc[2] - Jcg[2][1] * srnc[1]) * un[0] +
                  -w[n] * (Jcg[2][1] * srnc[0] - Jcg[0][1] * srnc[2]) * un[1] +
                  -w[n] * (Jcg[0][1] * srnc[1] - Jcg[1][1] * srnc[0]) * un[2] +
                  2.0 * -w[n] * Lc * sLc * Jcg[0][1] * tn[0] +
                  -w[n] * Lc * Lc * sJcg[0][1] * tn[0] +
                  2.0 * -w[n] * Lc * sLc * Jcg[1][1] * tn[1] +
                  -w[n] * Lc * Lc * sJcg[1][1] * tn[1] +
                  2.0 * -w[n] * Lc * sLc * Jcg[2][1] * tn[2] +
                  -w[n] * Lc * Lc * sJcg[2][1] * tn[2];  // thetay
      }

      if (dep_dof_constrained[5]) {
        res[5] += -w[n] * (sJcg[1][2] * rnc[2] - sJcg[2][2] * rnc[1]) * un[0] +
                  -w[n] * (sJcg[2][2] * rnc[0] - sJcg[0][2] * rnc[2]) * un[1] +
                  -w[n] * (sJcg[0][2] * rnc[1] - sJcg[1][2] * rnc[0]) * un[2] -
                  w[n] * (Jcg[1][2] * srnc[2] - Jcg[2][2] * srnc[1]) * un[0] +
                  -w[n] * (Jcg[2][2] * srnc[0] - Jcg[0][2] * srnc[2]) * un[1] +
                  -w[n] * (Jcg[0][2] * srnc[1] - Jcg[1][2] * srnc[0]) * un[2] +
                  2.0 * -w[n] * Lc * sLc * Jcg[0][2] * tn[0] +
                  -w[n] * Lc * Lc * sJcg[0][2] * tn[0] +
                  2.0 * -w[n] * Lc * sLc * Jcg[1][2] * tn[1] +
                  -w[n] * Lc * Lc * sJcg[1][2] * tn[1] +
                  2.0 * -w[n] * Lc * sLc * Jcg[2][2] * tn[2] +
                  -w[n] * Lc * Lc * sJcg[2][2] * tn[2];  // thetaz
      }

      Xn += 3;
      sXn += 3;
      un += NUM_DISPS;
      tn += NUM_DISPS;
    }

    // Step back to beginning of residual vector
    res -= NUM_VARIABLES - NUM_DISPS;
    // Rescale all equations so they are on the same order of magnitude as the
    // global stiffness matrix
    for (int i = 0; i < NUM_VARIABLES; i++) {
      res[i] *= C1;
    }

    res += NUM_VARIABLES;
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

  delete[] residual;
  delete[] sXpts;
  delete[] maskedVars;
}

/* Find centroid of dependent nodes */
void TACSRBE3::getCG(TacsScalar Xcg[], TacsScalar W[], const double w[],
                     const TacsScalar Xpts[]) {
  Xcg[0] = Xcg[1] = Xcg[2] = 0.0;
  W[0] = W[1] = W[2] = 0.0;
  // Skip the first node, because it's dependent
  Xpts += 3;

  for (int n = 0; n < NUM_INDEP_NODES; n++) {
    if (indep_dof_constrained[n][0]) {
      Xcg[0] += w[n] * Xpts[0];
      W[0] += w[n];
    }
    if (indep_dof_constrained[n][1]) {
      Xcg[1] += w[n] * Xpts[1];
      W[1] += w[n];
    }
    if (indep_dof_constrained[n][2]) {
      Xcg[2] += w[n] * Xpts[2];
      W[2] += w[n];
    }

    Xpts += 3;
  }

  Xcg[0] /= W[0];
  Xcg[1] /= W[1];
  Xcg[2] /= W[2];
}

void TACSRBE3::getCGSens(TacsScalar sXcg[], TacsScalar Xcg[], TacsScalar W[],
                         const double w[], const TacsScalar Xpts[],
                         const int component) {
  sXcg[0] = sXcg[1] = sXcg[2] = 0.0;
  getCG(Xcg, W, w, Xpts);
  // Skip the first node, because it's dependent, don't include dummy node
  if (component > 2 && component < (NUM_NODES - 1) * 3) {
    int k = component % 3;      // find dimension number
    int n = component / 3 - 1;  // find indep node number
    if (indep_dof_constrained[n][k]) {
      sXcg[k] = w[n] / W[k];
    }
  }
}

/* Find inverse moment of inertia of dependent nodes about centroid */
TacsScalar TACSRBE3::getMomentsOfInertia(TacsScalar IcgInv[3][3],
                                         const double w[],
                                         const TacsScalar Xpts[],
                                         const TacsScalar Xcg[]) {
  TacsScalar Icg[9], r[3], temp[9], c[6], Lc;
  Lc = 0.0;
  for (int j = 0; j < 9; j++) {
    Icg[j] = 0.0;
  }
  // Skip the first node, because it's dependent
  Xpts += 3;

  // Compute the average distance of the indep nodes from the cg
  for (int n = 0; n < NUM_INDEP_NODES; n++) {
    r[0] = Xpts[3 * n + 0] - Xcg[0];
    r[1] = Xpts[3 * n + 1] - Xcg[1];
    r[2] = Xpts[3 * n + 2] - Xcg[2];
    Lc += vec3Dot(r, r) / TacsScalar(NUM_INDEP_NODES);
  }

  Lc = sqrt(Lc);
  // If Lc is zero, set it to 1.0
  if (TacsRealPart(Lc) < SMALL_NUM) {
    Lc = 1.0;
  }

  // Compute moment of inertia tensor
  for (int n = 0; n < NUM_INDEP_NODES; n++) {
    // Specify the dof weight for each independent node (see ref [1])
    for (int j = 0; j < 6; j++) {
      if (indep_dof_constrained[n][j]) {
        c[j] = w[n];
      } else {
        c[j] = 0.0;
      }
    }

    r[0] = Xpts[3 * n + 0] - Xcg[0];
    r[1] = Xpts[3 * n + 1] - Xcg[1];
    r[2] = Xpts[3 * n + 2] - Xcg[2];
    Icg[0] += c[2] * r[1] * r[1] + c[1] * r[2] * r[2] + c[3] * Lc * Lc;  // Ixx
    Icg[1] += c[2] * -(r[0] * r[1]);                                     // Ixy
    Icg[2] += c[1] * -(r[0] * r[2]);                                     // Ixz
    Icg[3] = Icg[1];                                                     // Iyx
    Icg[4] += c[2] * r[0] * r[0] + c[0] * r[2] * r[2] + c[4] * Lc * Lc;  // Iyy
    Icg[5] += c[0] * -(r[1] * r[2]);                                     // Iyz
    Icg[6] = Icg[2];                                                     // Izx
    Icg[7] = Icg[5];                                                     // Izy
    Icg[8] += c[1] * r[0] * r[0] + c[0] * r[1] * r[1] + c[5] * Lc * Lc;  // Izz
  }

  // Invert the moment of inertia matrix
  TacsScalar detIcg = inv3x3(Icg, temp);

  // Check if the determinate is zero, if so at least one of the principal
  // moment of inertias is zero ( the points are colinear) and the moment of
  // inertia matrix cannot be inverted. Warn the user that they must specify
  // additional independent DOF's
  if (TacsRealPart(detIcg) < SMALL_NUM) {
    fprintf(stderr,
            "TACSRBE3: Error, points are colinear. "
            "Make sure to specify sufficient independent DOF's to react all "
            "moments.\n");
  }

  // Recast result into 2D array
  for (int row = 0; row < 3; row++) {
    for (int col = 0; col < 3; col++) {
      IcgInv[row][col] = temp[col + row * 3];
    }
  }

  return Lc;
}

/* Find sensitivity of inverse moment of inertia tensor of dependent nodes about
 * centroid */
TacsScalar TACSRBE3::getMomentsOfInertiaSens(
    TacsScalar sIcgInv[3][3], TacsScalar IcgInv[3][3], TacsScalar *sLc,
    const double w[], const TacsScalar Xpts[], const TacsScalar Xcg[],
    const TacsScalar sXcg[], const int component) {
  TacsScalar Icg[9], r[3], sr[3], temp[9], sIcg[9], stemp[9], c[6], Lc;
  Lc = 0.0;
  *sLc = 0.0;
  TacsScalar detIcg, sdetIcg;
  for (int j = 0; j < 9; j++) {
    Icg[j] = sIcg[j] = 0.0;
  }
  // Skip the first node, because it's dependent
  Xpts += 3;

  // Compute the average distance of the indep nodes from the cg
  for (int n = 0; n < NUM_INDEP_NODES; n++) {
    r[0] = Xpts[3 * n + 0] - Xcg[0];
    r[1] = Xpts[3 * n + 1] - Xcg[1];
    r[2] = Xpts[3 * n + 2] - Xcg[2];

    sr[0] = -sXcg[0];
    sr[1] = -sXcg[1];
    sr[2] = -sXcg[2];

    int m = component / 3 - 1;  // find independent node num
    if (n == m) {
      int k = component % 3;  // find dimension number
      sr[k] += 1.0;
    }

    Lc += vec3Dot(r, r) / TacsScalar(NUM_INDEP_NODES);
    *sLc += 2.0 * (r[0] * sr[0] + r[1] * sr[1] + r[2] * sr[2]) /
            TacsScalar(NUM_INDEP_NODES);
  }
  Lc = sqrt(Lc);
  *sLc = *sLc / 2.0 / Lc;
  // If Lc is zero, set it to 1.0
  if (TacsRealPart(Lc) < SMALL_NUM) {
    Lc = 1.0;
    *sLc = 0.0;
  }

  for (int n = 0; n < NUM_INDEP_NODES; n++) {
    // Specify the dof weight for each independent node (see ref [1])
    for (int j = 0; j < 6; j++) {
      if (indep_dof_constrained[n][j]) {
        c[j] = w[n];
      } else {
        c[j] = 0.0;
      }
    }

    r[0] = Xpts[3 * n + 0] - Xcg[0];
    r[1] = Xpts[3 * n + 1] - Xcg[1];
    r[2] = Xpts[3 * n + 2] - Xcg[2];

    sr[0] = -sXcg[0];
    sr[1] = -sXcg[1];
    sr[2] = -sXcg[2];

    int m = component / 3 - 1;  // find independent node num
    if (n == m) {
      int k = component % 3;  // find dimension number
      sr[k] += 1.0;
    }

    Icg[0] += c[2] * r[1] * r[1] + c[1] * r[2] * r[2] + c[3] * Lc * Lc;  // Ixx
    Icg[1] += c[2] * -(r[0] * r[1]);                                     // Ixy
    Icg[2] += c[1] * -(r[0] * r[2]);                                     // Ixz
    Icg[3] = Icg[1];                                                     // Iyx
    Icg[4] += c[2] * r[0] * r[0] + c[0] * r[2] * r[2] + c[4] * Lc * Lc;  // Iyy
    Icg[5] += c[0] * -(r[1] * r[2]);                                     // Iyz
    Icg[6] = Icg[2];                                                     // Izx
    Icg[7] = Icg[5];                                                     // Izy
    Icg[8] += c[1] * r[0] * r[0] + c[0] * r[1] * r[1] + c[5] * Lc * Lc;  // Izz

    sIcg[0] += 2.0 * (c[2] * r[1] * sr[1] + c[1] * r[2] * sr[2] +
                      c[3] * Lc * *sLc);               // Ixx
    sIcg[1] += c[2] * -(sr[0] * r[1] + r[0] * sr[1]);  // Ixy
    sIcg[2] += c[1] * -(sr[0] * r[2] + r[0] * sr[2]);  // Ixz
    sIcg[3] = sIcg[1];                                 // Iyx
    sIcg[4] += 2.0 * (c[2] * r[0] * sr[0] + c[0] * r[2] * sr[2] +
                      c[4] * Lc * *sLc);               // Iyy
    sIcg[5] += c[0] * -(sr[1] * r[2] + r[1] * sr[2]);  // Iyz
    sIcg[6] = sIcg[2];                                 // Izx
    sIcg[7] = sIcg[5];                                 // Izy
    sIcg[8] += 2.0 * (c[1] * r[0] * sr[0] + c[0] * r[1] * sr[1] +
                      c[5] * Lc * *sLc);  // Izz
  }

  // Find sensitivity of inverse matrix
  detIcg = inv3x3Sens(Icg, sIcg, temp, stemp, &sdetIcg);

  // Recast into 2D array
  for (int row = 0; row < 3; row++) {
    for (int col = 0; col < 3; col++) {
      IcgInv[row][col] = temp[col + row * 3];
      sIcgInv[row][col] = stemp[col + row * 3];
    }
  }

  return Lc;
}

/* This procedure masks the vars vector copied in depending on
which dofs have been constrained. If the dof is constrained copy
the component, if not force it to zero. */
void TACSRBE3::getMaskedVars(TacsScalar maskedVars[], const TacsScalar vars[]) {
  // Mask forces
  int i = 0;
  for (; i < NUM_DISPS; i++) {
    maskedVars[i] = vars[i];
  }
  for (int n = 0; n < NUM_INDEP_NODES; n++) {
    for (int j = 0; j < 6; j++, i++) {
      if (indep_dof_constrained[n][j]) {
        maskedVars[i] = vars[i];
      } else {
        maskedVars[i] = 0.0;
      }
    }
  }

  for (int j = 0; j < 6; j++, i++) {
    if (dep_dof_constrained[j]) {
      maskedVars[i] = vars[i];
    } else {
      maskedVars[i] = 0.0;
    }
  }
}

/*
  Retrieve the output data for this element
*/
/*void TACSRBE3::getOutputData( int elemIndex,
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
      for ( int n = 0; n < NUM_NODES-1; n++ ){
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
                data[index+k] = TacsRealPart(res[NUM_DISPS*n+k]);
              }
              else{
                // flip force on indep nodes to get force applied by RBE
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
