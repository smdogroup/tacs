#ifndef TACS_SPECTRAL_INTEGRATOR_H
#define TACS_SPECTRAL_INTEGRATOR_H

#include "TACSAssembler.h"
#include "TACSFunction.h"
#include "TACSMg.h"

// Forward declaration of spectral integrator
class TACSSpectralIntegrator;

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
  TACSLinearSpectralMat(TACSSpectralIntegrator *spectral);
  ~TACSLinearSpectralMat();

  TACSVec *createVec();
  void mult(TACSVec *x, TACSVec *y);

  void getMat(TACSParallelMat **Hmat, TACSParallelMat **Cmat);
  void setMatrixOrientation(MatrixOrientation matOr);
  MatrixOrientation getMatrixOrientation();
  int getFirstOrderCoefficients(const double *d[]);

 private:
  TACSSpectralIntegrator *spectral;
  int N;
  TACSParallelMat *H, *C;
  TACSBVec *temp;
  MatrixOrientation orient;
};

class TACSLinearSpectralMg : public TACSPc {
 public:
  TACSLinearSpectralMg(TACSLinearSpectralMat *_mat, int _nlevels,
                       TACSAssembler **assembler, TACSBVecInterp **_interp,
                       int coarsen_time[] = NULL);
  ~TACSLinearSpectralMg();

  void factor();
  void applyFactor(TACSVec *in, TACSVec *out);

 private:
  // Apply multigrid at the given level
  void applyMg(int level);

  // Class that stores a linear combo of matrices
  class TACSComboMat : public TACSMat {
   public:
    TACSComboMat(TACSParallelMat *Hmat, TacsScalar _alpha,
                 TACSParallelMat *Cmat, TACSBVec *tmp);
    ~TACSComboMat();

    TACSVec *createVec();
    void mult(TACSVec *x, TACSVec *y);

   private:
    TacsScalar alpha;
    TACSParallelMat *H, *C;
    TACSBVec *temp;

    friend class MgData;
  };

  // Class to store the multigrid data for each level
  class MgData {
   public:
    MgData(MgData *_fine, int Nval, const double d0[],
           TACSAssembler *_assembler, TACSBVecInterp *_interp,
           TACSParallelMat *Hmat = NULL, TACSParallelMat *Cmat = NULL,
           int _direct = 0);
    ~MgData();

    // Pointer to the next finest level of data
    MgData *fine;

    // The size of the spectral space at this mesh level
    int N;

    // Coefficients for the first-order approx.
    double *d0;

    // Full spectral vectors for multigrid
    TACSSpectralVec *x, *b, *r;

    void factor();  // Factor the matrix
    void applyFactor(TACSSpectralVec *in, TACSSpectralVec *out);
    void mult(TACSSpectralVec *in, TACSSpectralVec *out);
    void restriction(TACSSpectralVec *in, TACSSpectralVec *out);
    void interpolateAdd(TACSSpectralVec *in, TACSSpectralVec *out);

    // Set the matrix orientation
    void setMatrixOrientation(MatrixOrientation matOr);

   private:
    // The matrix orientation
    MatrixOrientation orient;

    // Temporary vector for matrix-vector product operations
    TACSBVec *temp;

    // Smoother at each time step
    TACSBVec *mats_temp;
    TACSComboMat **mats;
    TACSChebyshevSmoother **smoothers;

    // Direct solution data on the coarsest mesh
    int direct;  // Flag indicating if we're using a direct solution strategy
    TACSParallelMat **direct_mats;
    TACSBlockCyclicPc **direct_pcs;

    // Not allocated if no multigrid in time
    double *w0, *w1;  // Weights for the temporal interpolation

    TACSBVecInterp *interp;    // Spatial interpolation
    TACSAssembler *assembler;  // Assembler
    TACSParallelMat *H, *C;    // Matrices at this grid level
  };

  int nlevels;
  MgData **data;
  TACSLinearSpectralMat *mat;
};

class TACSSpectralIntegrator : public TACSObject {
 public:
  TACSSpectralIntegrator(TACSAssembler *_assembler, double tfinal, int N);
  ~TACSSpectralIntegrator();

  int getNumLGLNodes();
  double getPointAtLGLNode(int index);
  double getTimeAtLGLNode(int index);
  double getWeightAtLGLNode(int index);
  TACSAssembler *getAssembler();
  int getFirstOrderCoefficients(const double *d[]);

  TACSSpectralVec *createVec();
  TACSLinearSpectralMat *createLinearMat();
  void setInitialConditions(TACSBVec *init);
  void setVariables(TACSSpectralVec *vec);
  void assembleRes(TACSSpectralVec *res);
  void assembleMat(TACSLinearSpectralMat *mat,
                   MatrixOrientation matOr = TACS_MAT_NORMAL);

  // Evaluate the functions of interest and their derivatives
  void evalFunctions(int num_funcs, TACSFunction **funcs, TacsScalar *fvals);
  void evalSVSens(TACSFunction *func, TACSSpectralVec *dfdu);
  void addDVSens(TACSFunction *func, TACSBVec *dfdx);
  void addAdjointResProduct(TacsScalar scale, TACSSpectralVec *adjoint,
                            TACSBVec *dfdx);

  // Compute the time derivative at a point
  void computeDeriv(int index, TACSSpectralVec *sol, TACSBVec *dudt,
                    int include_init_conditions = 1);
  void computeDerivTranspose(int index, TACSSpectralVec *sol, TACSBVec *dudt);

  // Compute the solution at a point in the time interval
  void computeSolutionAndDeriv(double time, TACSSpectralVec *sol, TACSBVec *u,
                               TACSBVec *dudt = NULL);

 private:
  // Initialize the LGL points and weights
  void initLGLPointsAndWeights(int max_newton_iters = 100, double tol = 1e-15);

  // Initialize the full-order and first-order derivative operators
  void initOperator();

  // Compute the interpolation at a point
  void evalInterpolation(double pt, double P[], double Px[]);

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