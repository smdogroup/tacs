/*
  This example solves the following eigenvalue problem with the JD solver:
    A*u = lambda*u
  where:
    A = K - lambda0*M
  where K, M are stiffness and mass matrix for a finite element system
*/

#include <stdlib.h>

#include "JacobiDavidson.h"
#include "ParOptOptimizer.h"
#include "TACSAssembler.h"
#include "TACSConstitutive.h"
#include "TACSElement3D.h"
#include "TACSHexaBasis.h"
#include "TACSLinearElasticity.h"
#include "TACSMaterialProperties.h"
#include "TACSStructuralMass.h"
#include "TACSToFH5.h"
#include "TMRBoundaryConditions.h"
#include "TMREgads.h"
#include "TMRHelmholtzPUFilter.h"
#include "TMRLagrangeFilter.h"
#include "TMRMatrixFilter.h"
#include "TMRMesh.h"
#include "TMROctConstitutive.h"
#include "TMROctForest.h"
#include "TMRTopoProblem.h"
#include "TMRTopology.h"
#include "TMR_RefinementTools.h"
#include "TMR_TACSTopoCreator.h"
#include "nlopt.hpp"  // Nonlinear optimizer library

class OctCreator : public TMROctConformTACSTopoCreator {
 public:
  OctCreator(TMRBoundaryConditions *_bcs, TMROctForest *_forest,
             TMRStiffnessProperties *_stiff_props);

  TACSElement *createElement(int order, TMROctant *oct, int nweights,
                             const int *index, TMROctForest *filter);

 private:
  TMRBoundaryConditions *bcs;
  TMROctForest *forest;
  TMRStiffnessProperties *stiff_props;
};

class CreatorCallback : public TMREntity {
 public:
  CreatorCallback(TMRBoundaryConditions *_bcs,
                  TMRStiffnessProperties *_stiff_props);

  OctCreator *creator_callback(TMROctForest *forest);

 private:
  TMRBoundaryConditions *bcs;
  TMRStiffnessProperties *stiff_props;
};

/*
  Note: We only consider dim = 3 in this
  specific implementation
*/
class OptFilterWeights {
 public:
  OptFilterWeights(int _diagonal_index, int _npts, int _dim,
                   const TacsScalar _Xpts[], double _H[], int _n_entries_H);

  ~OptFilterWeights();

  void set_alphas(double *w, double *alpha);

  int get_size_A() { return size_A; }

  int get_size_b() { return size_b; }

  int get_ncon() { return ncon; }

  double *get_A() { return A; }

  double *get_b() { return b; }

  void obj_counter_inc() { nfev_obj++; }

  void con_counter_inc() { nfev_con++; }

  int get_nfev_obj() { return nfev_obj; }

  int get_nfev_con() { return nfev_con; }

 private:
  int diagonal_index, dim, npts, n_entries_H;
  double delta;
  const TacsScalar *Xpts;

  int size_A, size_b, ncon;
  double *A, *b, *H;

  int nfev_obj;
  int nfev_con;
};

class Mfilter : public TMRHelmholtzPUFilter {
 public:
  Mfilter(int _N, int _nlevels, TACSAssembler *_assembler[],
          TMROctForest *_filter[], double _r);

  int getInteriorStencil(int diagonal_index, int npts, const TacsScalar Xpts[],
                         double alpha[]);

  int getBoundaryStencil(int diagonal_index, const TacsScalar n[], int npts,
                         const TacsScalar Xpts[], double alpha[]);

 private:
  int mpi_rank;
  double r;
  int dim, vars_per_node;
  OptFilterWeights *opt_data;
};

class MFilterCreator : public TMREntity {
 public:
  MFilterCreator(double _r0_frac, int _N, double _a);

  // TMRLagrangeFilter* create_filter( int nlevels,
  //                                   TACSAssembler *assemblers[],
  //                                   TMROctForest *filters[]);

  Mfilter *create_filter(int nlevels, TACSAssembler *assemblers[],
                         TMROctForest *filters[]);

 private:
  double r0_frac, a;
  int N;
};

class TopoProblemCreator : public TMREntity {
 public:
  TopoProblemCreator(TMROctForest *_forest,
                     CreatorCallback *_creator_callback_obj,
                     MFilterCreator *_filter_creator, int _nlevels,
                     int _use_galerkin);

  ~TopoProblemCreator();

  TMRTopoProblem *createTopoProblem();

 private:
  TMROctForest *forest;
  CreatorCallback *creator_callback_obj;
  MFilterCreator *filter_creator;
  int nlevels;
  int use_galerkin;

  TMROctForest **forests;
  TMROctForest **filters;
  TACSAssembler **assemblers;
  TACSMg *mg;
};

/*
  This class evaluates the structural mass as an objective
*/
class MassObj : public ParOptBase {
 public:
  MassObj(double _m_fixed);

  ~MassObj();

  void evalObj(TMRTopoFilter *_filter, TACSMg *dummy, TacsScalar *val);

  void evalGrad(TMRTopoFilter *dummy1, TACSMg *dummy2, TACSBVec *dfdx);

 private:
  int allocated;
  double m_fixed;
  TMRTopoFilter *filter;
  TACSAssembler *assembler;
  TACSFunction *mass_func;
};

/*
  This class evaluates the frequency constraint
*/
class FrequencyConstr : public ParOptBase {
 public:
  FrequencyConstr(TMROctForest *_forest, double _len0, double _lambda0,
                  int _max_jd_size, int _max_gmres_size, double _ksrho,
                  double _non_design_mass);

  ~FrequencyConstr();

  void evalConstr(TMRTopoFilter *_filter, TACSMg *_mg, int nconstr,
                  TacsScalar *vals);

  void evalConstrGrad(TMRTopoFilter *dummy1, TACSMg *dummy2, int nconstr,
                      TACSBVec **vecs);

  void qn_correction(ParOptVec *x, ParOptScalar z[], ParOptVec *dummy,
                     ParOptVec *s, ParOptVec *y);

 private:
  MPI_Comm comm;
  int rank;

  double len0, lambda0, ksrho, non_design_mass;
  int num_eigenvalues, max_jd_size, max_gmres_size;

  int allocated;

  KSMPrint *ksm_print;

  TMROctForest *forest;
  TMRTopoFilter *filter;
  TACSMg *mg;
  TACSAssembler *assembler;

  TACSBVec *svec, *mvec, *dv;
  TACSParallelMat *mmat, *m0mat, *Amat;

  TacsScalar *eig, *eta;
  TACSBVec **eigv, **deig;

  TACSBVec *rho, *rho_original, *update, *temp;

  TACSJacobiDavidsonOperator *oper;
  TACSJacobiDavidson *jd;

  TacsScalar *mvals;
  TMRPoint *Xpts;

  TACSMat *mgmat;

  ParOptBVecWrap *s_wrap, *y_wrap;
  TACSBVec *s_tacs, *y_tacs;
};

class OutputCallback : public ParOptBase {
 public:
  OutputCallback(TACSAssembler *_assembler, int _iter_offset);

  ~OutputCallback();

  void write_output() {}

 private:
  TACSAssembler *assembler;
  int iter_offset;
  int flag;
  TACSToFH5 *f5;
  TACSBVec *xt;
};

/*
  Create the TMRTopoProblem object and set up the topology optimization
  problem
*/
class ProblemCreator : public ParOptBase {
 public:
  ProblemCreator(TMROctForest *forest, TMRBoundaryConditions *bcs,
                 TMRStiffnessProperties *stiff_props, int nlevels,
                 double lambda0, double ksrho, double vol_frac, double r0_frac,
                 double len0, double density, int qn_correction,
                 double non_design_mass, int max_jd_size, int max_gmres_size);

  ~ProblemCreator();

  TMRTopoProblem *get_problem();

 private:
  TopoProblemCreator *tpc;
  TMRTopoProblem *problem;
  MFilterCreator *filter_creator;
  CreatorCallback *creator_callback_obj;
  MassObj *mass_obj;
  FrequencyConstr *freq_con;
};

/* Evaluate objective and constraint for nlopt_optimize */
double objective(unsigned n, const double *w, double *grad, void *data);

void constraint(unsigned ncon, double *cons, unsigned n, const double *w,
                double *grad, void *data);

/* Helper functions for calling member function via function pointer */

void forwarder_objval(void *massobj, TMRTopoFilter *filter, TACSMg *dummy,
                      TacsScalar *val) {
  static_cast<MassObj *>(massobj)->evalObj(filter, dummy, val);
}

void forwarder_objgrad(void *massobj, TMRTopoFilter *dummy1, TACSMg *dummy2,
                       TACSBVec *dfdx) {
  static_cast<MassObj *>(massobj)->evalGrad(dummy1, dummy2, dfdx);
}

void forwarder_conval(void *freqcon, TMRTopoFilter *filter, TACSMg *mg,
                      int nconstr, TacsScalar *vals) {
  static_cast<FrequencyConstr *>(freqcon)->evalConstr(filter, mg, nconstr,
                                                      vals);
}

void forwarder_congrad(void *freqcon, TMRTopoFilter *dummy1, TACSMg *dummy2,
                       int nconstr, TACSBVec **vecs) {
  static_cast<FrequencyConstr *>(freqcon)->evalConstrGrad(dummy1, dummy2,
                                                          nconstr, vecs);
}

void forwarder_qn(int dummy1, void *freqcon, ParOptVec *x, ParOptScalar z[],
                  ParOptVec *dummy2, ParOptVec *s, ParOptVec *y) {
  static_cast<FrequencyConstr *>(freqcon)->qn_correction(x, z, dummy2, s, y);
}
