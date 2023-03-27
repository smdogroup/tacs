#ifndef TACS_SPECTRAL_INTEGRATOR_H
#define TACS_SPECTRAL_INTEGRATOR_H

#include "TACSAssembler.h"
#include "TACSMg.h"

/*
  The spectral vector contains all the time instances used in the spectral
  expansion, except the initial condition
*/
class TACSSpectralVec : public TACSVec {
 public:
  TACSSpectralVec(int N, TACSAssembler *assembler);
  TACSSpectralVec(int N, TACSMat *mat);
  ~TACSSpectralVec();

  TacsScalar norm();
  void scale(TacsScalar alpha);
  TacsScalar dot(TACSVec *x);
  void axpy(TacsScalar alpha, TACSVec *x);
  void copyValues(TACSVec *x);
  void axpby(TacsScalar alpha, TacsScalar beta, TACSVec *x);
  void zeroEntries();

  TACSBVec *getVec(int index);

 private:
  int N;
  TACSBVec **vecs;
};

/*
  This linear matrix class is based on a constant coefficient in time ODE of the
  form

  C * du/dt + H * u
*/
class TACSLinearSpectralMat : public TACSMat {
 public:
  TACSLinearSpectralMat(int N, const double *D, TACSMat *Hmat, TACSMat *Cmat);
  ~TACSLinearSpectralMat();

  TACSVec *createVec();
  void mult(TACSVec *x, TACSVec *y);

  void getMat(TACSMat **Hmat, TACSMat **Cmat);
  void setMatrixOrientation(MatrixOrientation matOr);

 private:
  int N;
  const double *D;
  TACSMat *H, *C;
  TACSBVec *temp;
  MatrixOrientation orient;
};

class TACSSpectralIntegrator : public TACSObject {
 public:
  TACSSpectralIntegrator(TACSAssembler *_assembler, double tfinal, int N);
  ~TACSSpectralIntegrator();

  TACSSpectralVec *createVec();
  TACSLinearSpectralMat *createLinearMat();
  void setInitialConditions(TACSBVec *init);
  void setVariables(TACSSpectralVec *vec);
  void assembleRes(TACSSpectralVec *res);
  void assembleMat(TACSLinearSpectralMat *mat,
                   MatrixOrientation matOr = TACS_MAT_NORMAL);

 private:
  // Initialize the LGL points and weights
  void initLGLPointsAndWeights(int max_newton_iters = 100, double tol = 1e-15);

  // Initialize the full-order and first-order derivative operators
  void initOperator();

  // The TACSAssembler model
  TACSAssembler *assembler;

  // Initial values
  TACSBVec *init;

  // Spectral values
  TACSSpectralVec *vars;

  // Time values
  double tinit;    // The initial time
  double tfinal;   // The final time
  double tfactor;  // = 2.0 / (tfinal - tinit)

  // Points and weights for the LGL integration
  int N;
  double *pts;   // The N+1 LGL points on the computational domain
  double *wts;   // The N+1 LGL weights on the computational domain
  double *tpts;  // The N+1 points in time corresponding to the LGL points

  // The derivative operator at the quadrature point
  double *D;   // First derivative operator
  double *d0;  // First-order first derivative operator
};

#endif  // TACS_SPECTRAL_INTEGRATOR_H