/*
  This example solves the following eigenvalue problem with the JD solver:
    A*u = lambda*u
  where:
    A = K - lambda0*M
  where K, M are stiffness and mass matrix for a finite element system
*/

#include "JD.h"

OctCreator::OctCreator(TMRBoundaryConditions *_bcs, TMROctForest *_forest,
                       TMRStiffnessProperties *_stiff_props)
    : TMROctConformTACSTopoCreator(_bcs, 1, _forest, -1,
                                   TMR_GAUSS_LOBATTO_POINTS) {
  bcs = _bcs;
  forest = _forest;
  stiff_props = _stiff_props;
}

TACSElement *OctCreator::createElement(int order, TMROctant *oct, int nweights,
                                       const int *index, TMROctForest *filter) {
  TMROctConstitutive *con = new TMROctConstitutive(stiff_props, forest);
  con->incref();
  TACSLinearElasticity3D *model =
      new TACSLinearElasticity3D(con, TACS_LINEAR_STRAIN);
  model->incref();
  TACSLinearHexaBasis *basis = new TACSLinearHexaBasis();
  basis->incref();
  TACSElement *elem = new TACSElement3D(model, basis);
  elem->incref();
  return elem;
}

CreatorCallback::CreatorCallback(TMRBoundaryConditions *_bcs,
                                 TMRStiffnessProperties *_stiff_props) {
  bcs = _bcs;
  stiff_props = _stiff_props;
}

OctCreator *CreatorCallback::creator_callback(TMROctForest *forest) {
  OctCreator *creator = new OctCreator(bcs, forest, stiff_props);
  return creator;
}

Mfilter::Mfilter(int _N, int _nlevels, TACSAssembler *_assembler[],
                 TMROctForest *_filter[], double _r)
    : TMRHelmholtzPUFilter(_N, _nlevels, _assembler, _filter) {
  // Get rank
  MPI_Comm_rank(_assembler[0]->getMPIComm(), &mpi_rank);

  // We hard-code dim to be 3 and vars per node to be 1
  vars_per_node = 1;
  dim = 3;
  r = _r;
}

int Mfilter::getInteriorStencil(int diagonal_index, int npts,
                                const TacsScalar Xpts[], double alpha[]) {
  // Create a 3-by-3 matrix H in flat(1d) format
  double H[9];
  memset(H, 0, sizeof(double) * 9);
  H[0] = H[4] = H[8] = r * r;  // set diagonal entries to be r^2

  // Create optimization problem data object
  opt_data = new OptFilterWeights(diagonal_index, npts, 3, Xpts, H, 9);

  // Create optimization problem
  int ndim = npts - 1;
  int ncon = opt_data->get_ncon();
  nlopt_opt opt = nlopt_create(NLOPT_LD_SLSQP, ndim);

  // Set objective
  nlopt_set_min_objective(opt, objective, opt_data);

  // Set lower bounds
  double lb[ndim];
  memset(lb, 0, sizeof(double) * ndim);
  nlopt_set_lower_bounds(opt, lb);

  // Add equality constraints
  double cons_tol[ncon];
  for (int i = 0; i < ncon; i++) {
    cons_tol[i] = 1e-6;
  }
  nlopt_add_equality_mconstraint(opt, ncon, constraint, opt_data, cons_tol);

  // Set tolerance
  nlopt_set_ftol_rel(opt, 1e-4);
  nlopt_set_xtol_rel(opt, 1e-4);

  // Set initial point
  double w0[ndim];
  for (int i = 0; i < ndim; i++) {
    w0[i] = 1.0;
  }

  // Optimize
  double minf;
  int flag = nlopt_optimize(opt, w0, &minf);
  if (flag < 0) {
    if (mpi_rank == 0) {
      printf("[Error] Optimization failed for interior stencil!, flag = %d\n",
             flag);
    }
    exit(-1);
  }

  // if (mpi_rank == 0){
  //   printf("[MfilterInterior] Optimization converged, nfev_obj = %d, nfev_con
  //   = %d, flag = %d\n",
  //          opt_data->get_nfev_obj(), opt_data->get_nfev_con(), flag);
  // }

  // Set optimized alphas
  opt_data->set_alphas(w0, alpha);

  // Delete optimization problem data object
  delete opt_data;

  // Delete optimizer instance
  nlopt_destroy(opt);
}

int Mfilter::getBoundaryStencil(int diagonal_index, const TacsScalar n[],
                                int npts, const TacsScalar Xpts[],
                                double alpha[]) {
  // Create a 2-by-2 matrix H in flat format
  double H[4];
  memset(H, 0, sizeof(double) * 4);
  H[0] = H[3] = r * r;

  /* Reduce to a 2d problem */

  int index = 0;
  TacsScalar temp = fabs(n[0]);
  for (int i = 0; i < 3; i++) {
    if (fabs(n[i]) < temp) {
      temp = fabs(n[i]);
      index = i;
    }
  }

  // Compute cross products
  TacsScalar t[3] = {0.0, 0.0, 0.0};
  TacsScalar t1[3] = {0.0, 0.0, 0.0};
  TacsScalar t2[3] = {0.0, 0.0, 0.0};
  t[index] = 1.0;

  // t2 = cross(t, n)
  t2[0] = t[1] * n[2] - t[2] * n[1];
  t2[1] = t[2] * n[0] - t[0] * n[2];
  t2[2] = t[0] * n[1] - t[1] * n[0];

  // t1 = cross(n,t2)
  t1[0] = n[1] * t2[2] - n[2] * t2[1];
  t1[1] = n[2] * t2[0] - n[0] * t2[2];
  t1[2] = n[0] * t2[1] - n[1] * t2[0];

  // Reduce the problem
  TacsScalar Xt[2 * npts];

  // Compute diff = X - X[diag, :]
  TacsScalar Xdiag[3];
  for (int i = 0; i < 3; i++) {
    Xdiag[i] = Xpts[3 * diagonal_index + i];
  }
  TacsScalar diff[3 * npts];
  for (int i = 0; i < npts; i++) {
    for (int j = 0; j < 3; j++) {
      diff[3 * i + j] = Xpts[3 * i + j] - Xdiag[j];
    }
  }

  // Xt[:,0] = dot(diff, t1)
  // Xt[:,1] = dot(diff, t2)
  for (int i = 0; i < npts; i++) {
    Xt[2 * i] = 0.0;
    Xt[2 * i + 1] = 0.0;
    for (int j = 0; j < 3; j++) {
      Xt[2 * i] += diff[3 * i + j] * t1[j];
      Xt[2 * i + 1] += diff[3 * i + j] * t2[j];
    }
  }

  // Create optimization problem data object
  opt_data = new OptFilterWeights(diagonal_index, npts, 2, Xt, H, 4);

  int ndim = npts - 1;
  // Create optimization problem
  int ncon = opt_data->get_ncon();
  nlopt_opt opt = nlopt_create(NLOPT_LD_SLSQP, ndim);

  // Set objective
  nlopt_set_min_objective(opt, objective, opt_data);

  // Set lower bounds
  double lb[ndim];
  memset(lb, 0, sizeof(double) * ndim);
  nlopt_set_lower_bounds(opt, lb);

  // Add equality constraints
  double cons_tol[ncon];
  for (int i = 0; i < ncon; i++) {
    cons_tol[i] = 1e-6;
  }
  nlopt_add_equality_mconstraint(opt, ncon, constraint, opt_data, cons_tol);

  // Set tolerance
  nlopt_set_ftol_rel(opt, 1e-4);
  nlopt_set_xtol_rel(opt, 1e-4);

  // Set initial point
  double w0[ndim];
  for (int i = 0; i < ndim; i++) {
    w0[i] = 1.0;
  }

  // Optimize
  double minf;
  int flag = nlopt_optimize(opt, w0, &minf);
  if (flag < 0) {
    if (mpi_rank == 0) {
      printf("[Error] Optimization failed for boundary stencil!, flag = %d\n",
             flag);
    }
    exit(-1);
  }

  // if (mpi_rank == 0){
  //   printf("[MfilterBoundary] Optimization converged, nfev_obj = %d, nfev_con
  //   = %d, flag = %d\n",
  //          opt_data->get_nfev_obj(), opt_data->get_nfev_con(), flag);
  // }

  // Set optimized alphas
  opt_data->set_alphas(w0, alpha);

  // Delete optimization problem data object
  delete opt_data;

  // Delete optimizer instance
  nlopt_destroy(opt);
}

MFilterCreator::MFilterCreator(double _r0_frac, int _N, double _a) {
  r0_frac = _r0_frac;
  N = _N;
  a = _a;
}

// TMRLagrangeFilter* MFilterCreator::create_filter( int nlevels,
//                                                   TACSAssembler
//                                                   *assemblers[], TMROctForest
//                                                   *filters[]){
//   TMRLagrangeFilter *mfilter = new TMRLagrangeFilter(nlevels, assemblers,
//   filters); mfilter->incref();
//   // mfilter->initialize();
//   return mfilter;
// }

Mfilter *MFilterCreator::create_filter(int nlevels, TACSAssembler *assemblers[],
                                       TMROctForest *filters[]) {
  double r = a * r0_frac;
  Mfilter *mfilter = new Mfilter(N, nlevels, assemblers, filters, r);
  mfilter->incref();
  mfilter->initialize();
  return mfilter;
}

TopoProblemCreator::TopoProblemCreator(TMROctForest *_forest,
                                       CreatorCallback *_creator_callback_obj,
                                       MFilterCreator *_filter_creator,
                                       int _nlevels, int _use_galerkin) {
  forest = _forest;
  creator_callback_obj = _creator_callback_obj;
  filter_creator = _filter_creator;
  nlevels = _nlevels;
  use_galerkin = _use_galerkin;
}

TopoProblemCreator::~TopoProblemCreator() {
  for (int i = 0; i < nlevels; i++) {
    if (assemblers[i]) {
      assemblers[i]->decref();
    }
    if (filters[i]) {
      filters[i]->decref();
    }
    if (forests[i]) {
      forests[i]->decref();
    }
    if (mg) {
      mg->decref();
    }
  }
}

TMRTopoProblem *TopoProblemCreator::createTopoProblem() {
  int repartition = 1;
  int lowest_order = 2;
  TACSAssembler::OrderingType ordering = TACSAssembler::MULTICOLOR_ORDER;

  // Store hierarchy data
  TMROctForest *forests[nlevels];
  TMROctForest *filters[nlevels];
  TACSAssembler *assemblers[nlevels];

  // Balance the forest and repartition across processors
  forest->balance(1);
  if (repartition) {
    forest->repartition();
  }

  // Create the forest object
  OctCreator *creator = creator_callback_obj->creator_callback(forest);
  creator->incref();
  forests[0] = forest;
  filters[0] = forest;
  assemblers[0] = creator->createTACS(forest, ordering);
  assemblers[0]->incref();

  for (int i = 0; i < nlevels - 1; i++) {
    int order = forests[i]->getMeshOrder();
    TMRInterpolationType interp = forests[i]->getInterpType();

    TMROctForest *_forest = NULL;
    if (order > lowest_order) {
      _forest = forests[i]->duplicate();
      _forest->incref();
      order = order - 1;
      _forest->setMeshOrder(order, interp);
    } else {
      _forest = forests[i]->coarsen();
      _forest->incref();
      _forest->setMeshOrder(order, interp);
      _forest->balance();
      if (repartition) {
        _forest->repartition();
      }
    }

    creator->decref();
    creator = creator_callback_obj->creator_callback(_forest);
    forests[i + 1] = _forest;
    filters[i + 1] = _forest;
    assemblers[i + 1] = creator->createTACS(_forest, ordering);
    assemblers[i + 1]->incref();
  }
  creator->decref();

  // Create the multigrid object
  double omega = 1.0;
  int use_coarse_direct_solve = 1;
  int use_chebyshev_smoother = 0;
  TMR_CreateTACSMg(nlevels, assemblers, forests, &mg, omega, use_galerkin,
                   use_coarse_direct_solve, use_chebyshev_smoother);
  mg->incref();

  // Create the TMRTopoFilter object
  Mfilter *filter_obj =
      filter_creator->create_filter(nlevels, assemblers, filters);

  // Create problem
  int gmres_subspace = 50;
  double rtol = 1e-9;
  TMRTopoProblem *problem =
      new TMRTopoProblem(filter_obj, mg, gmres_subspace, rtol);
  return problem;
}

MassObj::MassObj(double _m_fixed) {
  m_fixed = _m_fixed;
  filter = NULL;
  allocated = 0;
}

MassObj::~MassObj() {
  if (allocated) {
    mass_func->decref();
  }
}

void MassObj::evalObj(TMRTopoFilter *_filter, TACSMg *dummy, TacsScalar *val) {
  if (!filter) {
    allocated = 1;

    filter = _filter;
    assembler = filter->getAssembler();
    mass_func = new TACSStructuralMass(assembler);
    mass_func->incref();
  }

  // Evaluate mass
  assembler->evalFunctions(1, &mass_func, val);
  *val = *val / m_fixed;

  return;
}

void MassObj::evalGrad(TMRTopoFilter *dummy1, TACSMg *dummy2, TACSBVec *dfdx) {
  dfdx->zeroEntries();
  assembler->addDVSens(1 / m_fixed, 1, &mass_func, &dfdx);
  return;
}

FrequencyConstr::FrequencyConstr(TMROctForest *_forest, double _len0,
                                 double _lambda0, int _max_jd_size,
                                 int _max_gmres_size, double _ksrho,
                                 double _non_design_mass) {
  allocated = 0;

  forest = _forest;
  len0 = _len0;
  lambda0 = _lambda0;
  num_eigenvalues = 10;
  max_jd_size = _max_jd_size;
  max_gmres_size = _max_gmres_size;
  ksrho = _ksrho;
  non_design_mass = _non_design_mass;

  filter = NULL;
}

FrequencyConstr::~FrequencyConstr() {
  if (allocated) {
    ksm_print->decref();

    svec->decref();
    mvec->decref();
    dv->decref();
    mmat->decref();
    m0mat->decref();
    Amat->decref();

    for (int i = 0; i < num_eigenvalues; i++) {
      eigv[i]->decref();
      deig[i]->decref();
    }

    delete[] eig;
    delete[] eta;
    delete[] eigv;
    delete[] deig;

    rho->decref();
    rho_original->decref();
    update->decref();
    temp->decref();

    oper->decref();
    jd->decref();
  }
}

void FrequencyConstr::evalConstr(TMRTopoFilter *_filter, TACSMg *_mg,
                                 int nconstr, TacsScalar *vals) {
  if (!allocated) {
    allocated = 1;

    // KSMprint object
    ksm_print = new KSMPrintStdout("JD", rank, 1);
    ksm_print->incref();

    // Get objects
    mg = _mg;
    filter = _filter;
    assembler = filter->getAssembler();

    comm = assembler->getMPIComm();
    MPI_Comm_rank(comm, &rank);

    // Allocate space vectors and matrices
    svec = assembler->createDesignVec();
    mvec = assembler->createDesignVec();
    dv = assembler->createDesignVec();
    mmat = assembler->createMat();
    m0mat = assembler->createMat();
    Amat = assembler->createMat();

    svec->incref();
    mvec->incref();
    dv->incref();
    mmat->incref();
    m0mat->incref();
    Amat->incref();

    eig = new TacsScalar[num_eigenvalues];
    eta = new TacsScalar[num_eigenvalues];
    eigv = new TACSBVec *[num_eigenvalues];
    deig = new TACSBVec *[num_eigenvalues];

    for (int i = 0; i < num_eigenvalues; i++) {
      eigv[i] = assembler->createVec();
      eigv[i]->incref();
      deig[i] = assembler->createDesignVec();
      deig[i]->incref();
    }

    // Allocate space for qn correction
    rho = assembler->createDesignVec();
    rho_original = assembler->createDesignVec();
    update = assembler->createDesignVec();
    temp = assembler->createDesignVec();

    rho->incref();
    rho_original->incref();
    update->incref();
    temp->incref();

    // Create operator
    oper = new TACSJDSimpleOperator(assembler, Amat, mg);
    oper->incref();

    // Create the Jacobi-Davidson eigensolver
    jd = new TACSJacobiDavidson(oper, num_eigenvalues, max_jd_size,
                                max_gmres_size);
    jd->incref();
    jd->setTolerances(1e-6, 1e-6, 1e-6, 1e-12);
    jd->setThetaCutoff(0.01);
    jd->setRecycle(num_eigenvalues, JD_NUM_RECYCLE);

    /* Populate the non-design mass matrix */

    int nmval = mvec->getArray(&mvals);
    int n_local_nodes = forest->getPoints(&Xpts);
    int offset = forest->getExtPreOffset();

    // Loop over all owned nodes and set non-design mass
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double lx, ly, lz;
    lx = ly = lz = len0;
    double tol = 1e-6;
    xmin = lx - tol;
    xmax = lx + tol;
    ymin = 0.25 * ly - tol;
    ymax = 0.75 * ly + tol;
    zmin = 0.0 * lz - tol;
    zmax = 0.2 * lz + tol;

    for (int i = offset; i < n_local_nodes; i++) {
      if (xmin < Xpts[i].x < xmax) {
        if (ymin < Xpts[i].y < ymax) {
          if (zmin < Xpts[i].z < zmax) {
            mvals[i - offset] = 1.0;
          }
        }
      }
    }

    // Assemble the m0 matrix
    assembler->getDesignVars(dv);    // Get current dv
    assembler->setDesignVars(mvec);  // set dv to mvec only and assemble mat
    assembler->assembleMatType(TACS_MASS_MATRIX, m0mat);
    m0mat->scale(non_design_mass);
    assembler->setDesignVars(dv);  // Set dv back
  }

  assembler->assembleMatType(TACS_STIFFNESS_MATRIX, Amat);
  mg->getMat(&mgmat);
  mgmat->copyValues(Amat);

  assembler->assembleMatType(TACS_MASS_MATRIX, mmat);
  mmat->axpy(1.0, m0mat);  // --> done assemble mass matrix

  Amat->axpy(-lambda0, mmat);  // --> done assemble A = K - lambda0*M
  assembler->applyBCs(Amat);

  mgmat->axpy(-0.95 * lambda0,
              mmat);  // --> done assemble mgmat = K - 0.95*lambda0*M
  assembler->applyBCs(mgmat);

  mg->assembleGalerkinMat();
  mg->factor();

  // Solve
  int print_level = 1;
  jd->solve(ksm_print, print_level);
  if (jd->getNumConvergedEigenvalues() < num_eigenvalues) {
    if (rank == 0) {
      printf("[Warning] JD failed to converge, starting rerun...\n");
    }
    jd->solve(ksm_print, print_level);
    int nconvd = jd->getNumConvergedEigenvalues();
    if (nconvd < num_eigenvalues) {
      if (rank == 0) {
        printf("[Error] No enough eigenvalues converged! (%d/%d)\n", nconvd,
               num_eigenvalues);
      }
      exit(-1);
    }
  }

  // Extract eigenvalues and eigenvectors
  for (int i = 0; i < num_eigenvalues; i++) {
    eig[i] = jd->extractEigenvector(i, eigv[i], NULL);
  }

  // Set first eigenvector for visualization
  assembler->setVariables(eigv[0]);

  // Compute KS aggregation
  TacsScalar eig_min = eig[0];
  for (int i = 0; i < num_eigenvalues; i++) {
    if (eig[i] < eig_min) {
      eig_min = eig[i];
    }
  }

  TacsScalar beta = 0.0;
  for (int i = 0; i < num_eigenvalues; i++) {
    eta[i] = exp(-ksrho * (eig[i] - eig_min));
    beta += eta[i];
  }

  TacsScalar ks = (eig_min - log(beta) / ksrho);
  for (int i = 0; i < num_eigenvalues; i++) {
    eta[i] /= beta;
  }

  vals[0] = ks;
}

void FrequencyConstr::evalConstrGrad(TMRTopoFilter *dummy1, TACSMg *dummy2,
                                     int nconstr, TACSBVec **vecs) {
  // We only have one constraint
  TACSBVec *dcdrho = vecs[0];

  dcdrho->zeroEntries();

  for (int i = 0; i < num_eigenvalues; i++) {
    TacsScalar coeff = eta[i];
    deig[i]->zeroEntries();
    assembler->addMatDVSensInnerProduct(coeff, TACS_STIFFNESS_MATRIX, eigv[i],
                                        eigv[i], deig[i]);
    assembler->addMatDVSensInnerProduct(-coeff * lambda0, TACS_MASS_MATRIX,
                                        eigv[i], eigv[i], deig[i]);
    deig[i]->beginSetValues(TACS_ADD_VALUES);
    deig[i]->endSetValues(TACS_ADD_VALUES);
    dcdrho->axpy(1.0, deig[i]);
  }
}

void FrequencyConstr::qn_correction(ParOptVec *x, ParOptScalar z[],
                                    ParOptVec *dummy, ParOptVec *s,
                                    ParOptVec *y) {
  s_wrap = dynamic_cast<ParOptBVecWrap *>(s);
  s_tacs = s_wrap->vec;
  y_wrap = dynamic_cast<ParOptBVecWrap *>(y);
  y_tacs = y_wrap->vec;

  TacsScalar h = 1e-8;
  update->zeroEntries();
  assembler->getDesignVars(rho);
  rho_original->copyValues(rho);

  svec->zeroEntries();
  filter->applyFilter(s_tacs, svec);
  rho->axpy(h, svec);
  temp->zeroEntries();

  // set density = rho + h*s
  assembler->setDesignVars(rho);

  for (int i = 0; i < num_eigenvalues; i++) {
    TacsScalar coeff1 = eta[i];
    assembler->addMatDVSensInnerProduct(coeff1, TACS_STIFFNESS_MATRIX, eigv[i],
                                        eigv[i], temp);
    TacsScalar coeff2 = -eta[i] * lambda0;
    assembler->addMatDVSensInnerProduct(coeff2, TACS_MASS_MATRIX, eigv[i],
                                        eigv[i], temp);
  }

  // Set density = rho
  assembler->setDesignVars(rho_original);

  for (int i = 0; i < num_eigenvalues; i++) {
    TacsScalar coeff1 = eta[i];
    assembler->addMatDVSensInnerProduct(-coeff1, TACS_STIFFNESS_MATRIX, eigv[i],
                                        eigv[i], temp);
    TacsScalar coeff2 = -eta[i] * lambda0;
    assembler->addMatDVSensInnerProduct(-coeff2, TACS_MASS_MATRIX, eigv[i],
                                        eigv[i], temp);
  }

  temp->beginSetValues(TACS_ADD_VALUES);
  temp->endSetValues(TACS_ADD_VALUES);

  update->axpy(1.0 / h, temp);

  if (svec->dot(update)) {
    filter->applyTranspose(update, update);
    y_tacs->axpy(z[0], update);
  }
}

OutputCallback::OutputCallback(TACSAssembler *_assembler, int _iter_offset) {
  assembler = _assembler;
  xt = assembler->createDesignVec();
  xt->incref();

  flag = TACS_OUTPUT_CONNECTIVITY | TACS_OUTPUT_NODES |
         TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_EXTRAS;
  f5 = new TACSToFH5(assembler, TACS_SOLID_ELEMENT, flag);
  f5->incref();

  iter_offset = _iter_offset;
}

OutputCallback::~OutputCallback() {
  xt->decref();
  f5->decref();
}

/*
  Create the TMRTopoProblem object and set up the topology optimization
  problem
*/
ProblemCreator::ProblemCreator(TMROctForest *forest, TMRBoundaryConditions *bcs,
                               TMRStiffnessProperties *stiff_props, int nlevels,
                               double lambda0, double ksrho, double vol_frac,
                               double r0_frac, double len0, double density,
                               int qn_correction, double non_design_mass,
                               int max_jd_size, int max_gmres_size) {
  // Create filter object
  int N = 20;
  filter_creator = new MFilterCreator(r0_frac, N, len0);
  creator_callback_obj = new CreatorCallback(bcs, stiff_props);
  filter_creator->incref();
  creator_callback_obj->incref();

  int use_galerkin = 1;

  tpc = new TopoProblemCreator(forest, creator_callback_obj, filter_creator,
                               nlevels, use_galerkin);
  problem = tpc->createTopoProblem();

  // Compute fixed mass
  double AR = 1.0;
  double lx = len0 * AR;
  double ly = len0;
  double lz = len0;
  double vol = lx * ly * lz;
  double m_fixed = vol_frac * vol * density;

  // Add objective callback
  mass_obj = new MassObj(m_fixed);
  mass_obj->incref();

  problem->addObjectiveCallback(mass_obj, forwarder_objval, mass_obj,
                                forwarder_objgrad);

  // Add constraint callback
  freq_con = new FrequencyConstr(forest, len0, lambda0, max_jd_size,
                                 max_gmres_size, ksrho, non_design_mass);
  freq_con->incref();

  problem->addConstraintCallback(1, 1, freq_con, forwarder_conval, freq_con,
                                 forwarder_congrad);

  // Use Quasi-Newton update correction if specified
  if (qn_correction) {
    problem->addQnCorrectionCallback(1, freq_con, forwarder_qn);
  }
}

ProblemCreator::~ProblemCreator() {
  filter_creator->decref();
  creator_callback_obj->decref();
  mass_obj->decref();
  freq_con->decref();
  tpc->decref();
}

TMRTopoProblem *ProblemCreator::get_problem() { return problem; }

OptFilterWeights::OptFilterWeights(int _diagonal_index, int _npts, int _dim,
                                   const TacsScalar _Xpts[], double _H[],
                                   int _n_entries_H) {
  nfev_con = 0;
  nfev_obj = 0;

  diagonal_index = _diagonal_index;
  npts = _npts;
  dim = _dim;    // Dimension of Xpts
  Xpts = _Xpts;  // flat array of length 3*npts
  H = _H;        //  flat array of length _n_entries_H
  n_entries_H = _n_entries_H;

  /*
    Compute the normalization
  */
  TacsScalar Xdiag[dim];
  for (int i = 0; i < dim; i++) {
    Xdiag[i] = Xpts[dim * diagonal_index + i];
  }

  TacsScalar diff[dim * npts];
  for (int i = 0; i < npts; i++) {
    for (int j = 0; j < dim; j++) {
      diff[dim * i + j] = Xpts[dim * i + j] - Xdiag[j];
    }
  }

  TacsScalar max_square_sum = 0.0;
  for (int j = 0; j < dim; j++) {
    max_square_sum += diff[j] * diff[j];
  }
  TacsScalar s;
  for (int i = 0; i < npts; i++) {
    s = 0.0;
    for (int j = 0; j < dim; j++) {
      s += diff[dim * i + j] * diff[dim * i + j];
    }
    if (max_square_sum < s) {
      max_square_sum = s;
    }
  }

  delta = sqrt(max_square_sum);

  /*
    Compute the constraint matrix
  */

  if (dim == 2) {
    ncon = 5;
    size_A = ncon * (npts - 1);
    size_b = ncon;

    // Create A and b
    A = new double[size_A];
    b = new double[size_b];
    memset(A, 0, sizeof(double) * size_A);
    memset(b, 0, sizeof(double) * size_b);

    // Populate b
    b[2] = H[0];        // H[0,0]
    b[3] = H[3];        // H[1,1]
    b[4] = 2.0 * H[1];  // H[0,1]

    // Populate A
    int index = 0;
    for (int i = 0; i < npts; i++) {
      if (i != diagonal_index) {
        TacsScalar dx = (Xpts[2 * i] - Xpts[2 * diagonal_index]) / delta;
        TacsScalar dy =
            (Xpts[2 * i + 1] - Xpts[2 * diagonal_index + 1]) / delta;

        A[0 * (npts - 1) + index] = dx;
        A[1 * (npts - 1) + index] = dy;
        A[2 * (npts - 1) + index] = 0.5 * dx * dx;
        A[3 * (npts - 1) + index] = 0.5 * dy * dy;
        A[4 * (npts - 1) + index] = dx * dy;
        index += 1;
      }
    }
  }

  else if (dim == 3) {
    ncon = 9;
    size_A = ncon * (npts - 1);
    size_b = ncon;

    // Create A matrix and b vector
    A = new double[size_A];
    b = new double[size_b];
    memset(A, 0, sizeof(double) * size_A);
    memset(b, 0, sizeof(double) * size_b);

    // Populate b
    b[3] = H[0];      // H[0, 0]
    b[4] = H[4];      // H[1, 1]
    b[5] = H[8];      // H[2, 2]
    b[6] = 2 * H[5];  // H[1, 2]
    b[7] = 2 * H[2];  // H[0, 2]
    b[8] = 2 * H[1];  // H[0, 1]

    // Populate A
    int index = 0;
    for (int i = 0; i < npts; i++) {
      if (i != diagonal_index) {
        TacsScalar dx = (Xpts[3 * i] - Xpts[3 * diagonal_index]) / delta;
        TacsScalar dy =
            (Xpts[3 * i + 1] - Xpts[3 * diagonal_index + 1]) / delta;
        TacsScalar dz =
            (Xpts[3 * i + 2] - Xpts[3 * diagonal_index + 2]) / delta;

        A[0 * (npts - 1) + index] = dx;
        A[1 * (npts - 1) + index] = dy;
        A[2 * (npts - 1) + index] = dz;

        A[3 * (npts - 1) + index] = 0.5 * dx * dx;
        A[4 * (npts - 1) + index] = 0.5 * dy * dy;
        A[5 * (npts - 1) + index] = 0.5 * dz * dz;

        A[6 * (npts - 1) + index] = dy * dz;
        A[7 * (npts - 1) + index] = dx * dz;
        A[8 * (npts - 1) + index] = dx * dy;
        index += 1;
      }
    }
  }

  else {
    printf("[Error] dim should be only 2 or 3\n");
    exit(-1);
  }
}

OptFilterWeights::~OptFilterWeights() {
  delete[] A;
  delete[] b;
}

void OptFilterWeights::set_alphas(double *w, double *alpha) {
  for (int i = 0; i < npts; i++) {
    alpha[i] = 0.0;
  }
  int index = 0;
  for (int i = 0; i < npts; i++) {
    if (i != diagonal_index) {
      alpha[i] = w[index] / (delta * delta);
      alpha[diagonal_index] += w[index] / (delta * delta);
      index += 1;
    }
  }

  alpha[diagonal_index] += 1.0;
}

double objective(unsigned n, const double *w, double *grad, void *data) {
  // Cast data to be type of OptFilterWeights class
  OptFilterWeights *opt_data = static_cast<OptFilterWeights *>(data);

  // Compute objective value and gradient
  double obj_val = 0.0;
  for (int i = 0; i < n; i++) {
    obj_val += 0.5 * w[i] * w[i];
    if (grad) {
      grad[i] = w[i];
    }
  }

  // Keep track of the counter
  opt_data->obj_counter_inc();

  return obj_val;
}

void constraint(unsigned ncon, double *cons, unsigned n, const double *w,
                double *grad, void *data) {
  // Cast data to be type of OptFilterWeights class
  OptFilterWeights *opt_data = static_cast<OptFilterWeights *>(data);

  // Compute constraint value:
  // c = dot(A, w) - b
  double *A = opt_data->get_A();
  double *b = opt_data->get_b();
  for (int m = 0; m < ncon; m++) {
    cons[m] = -b[m];
    for (int i = 0; i < n; i++) {
      cons[m] += A[m * n + i] * w[i];
    }
  }

  // Compute constraint gradient:
  // dc = A
  if (grad) {
    for (int i = 0; i < ncon * n; i++) {
      grad[i] = A[i];
    }
  }

  // Keep track of the counter
  opt_data->con_counter_inc();
}

int main(int argc, char *argv[]) {
  /* Parameters */

  int check_gradient = 0;

  double qval = 5.0;  // penalization parameter
  int mg_levels = 3;  // Multigrid levels
  int max_iter = 200;
  double lambda0 = 5.0;
  int qn_subspace = 2;
  double ksrho = 50.0;
  double vol_frac = 0.4;
  double r0_frac = 0.05;
  double len0 = 1.0;
  double density = 2600.0;
  int qn_correction = 1;
  double non_design_mass = 5.0;
  int max_jd_size = 100;
  int max_gmres_size = 30;
  const char prefix[] = "results";

  /* MPI-related */

  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  /* Material and stiffness properties */

  TACSMaterialProperties *mat_props =
      new TACSMaterialProperties(2600.0, 921.0, 70e3, 0.3, 100.0, 24e-6, 230.0);
  TMRStiffnessProperties *stiff_props = new TMRStiffnessProperties(
      1, &mat_props, qval, 0.2, 1e-3, TMR_RAMP_PENALTY, qval);
  mat_props->incref();
  stiff_props->incref();

  /* Create geometry, mesh, forest and boundary conditions */

  // Read egads file
  const char egads_name[] = "cantilever_1.0_1.0_1.0.egads";
  TMRModel *geo = TMR_EgadsInterface::TMR_LoadModelFromEGADSFile(egads_name, 0);
  geo->incref();

  // Retrieve faces and volume
  TMRFace **faces;
  TMRVolume **volumes;
  geo->getFaces(NULL, &faces);
  geo->getVolumes(NULL, &volumes);

  // Set name and source
  faces[0]->setName("fixed");
  faces[0]->setSource(volumes[0], faces[1]);

  // Create the mesh object
  TMRMesh *mesh = new TMRMesh(comm, geo);
  mesh->incref();

  // Create TMR options
  TMRMeshOptions *options = new TMRMeshOptions();

  // Mesh the mesh
  double htarget = 1.0;
  mesh->mesh(*options, htarget);

  // Create a topology object with underlying mesh geometry
  TMRModel *model = mesh->createModelFromMesh();
  model->incref();
  TMRTopology *topo = new TMRTopology(comm, model);
  topo->incref();

  // Create the oct forest
  TMROctForest *forest = new TMROctForest(comm);
  forest->incref();

  // Set topology for forest
  forest->setTopology(topo);
  forest->createTrees(mg_levels - 1);

  // Export mesh to vtk
  const char *vtk = "test.vtk";
  forest->writeForestToVTK(vtk);

  /* Set boundary conditions */

  TMRBoundaryConditions *bcs = new TMRBoundaryConditions();
  bcs->incref();
  const int bcs_nums[3] = {0, 1, 2};
  const double bcs_vals[3] = {0.0, 0.0, 0.0};
  bcs->addBoundaryCondition("fixed", 3, bcs_nums, bcs_vals);

  /* Create topology problem instance */
  ProblemCreator *pcr =
      new ProblemCreator(forest, bcs, stiff_props, mg_levels, lambda0, ksrho,
                         vol_frac, r0_frac, len0, density, qn_correction,
                         non_design_mass, max_jd_size, max_gmres_size);
  pcr->incref();
  TMRTopoProblem *problem = pcr->get_problem();
  problem->incref();

  // Set the prefix
  problem->setPrefix(prefix);

  // Initialize the problem
  problem->initialize();

  // Check gradient
  if (check_gradient) {
    // ParOptVec *xt = problem->createDesignVec();
    // xt->incref();
    // xt->set(0.95);
    // problem->setInitDesignVars(xt);
    problem->checkGradients(1e-6);
  }

  // Set up options
  ParOptOptions *paropt_options = new ParOptOptions(comm);
  paropt_options->incref();
  ParOptOptimizer::addDefaultOptions(paropt_options);

  paropt_options->setOption("algorithm", "tr");
  paropt_options->setOption("output_level", 0);
  paropt_options->setOption("norm_type", "l1");
  paropt_options->setOption("tr_init_size", 0.05);
  paropt_options->setOption("tr_min_size", 1e-3);
  paropt_options->setOption("tr_max_size", 1.0);
  paropt_options->setOption("tr_eta", 0.25);
  paropt_options->setOption("tr_infeas_tol", 1e-6);
  paropt_options->setOption("tr_l1_tol", 0.0);
  paropt_options->setOption("tr_linfty_tol", 0.0);
  paropt_options->setOption("tr_adaptive_gamma_update", 0);
  paropt_options->setOption("tr_accept_step_strategy", "filter_method");
  paropt_options->setOption("filter_sufficient_reduction", 1);
  paropt_options->setOption("filter_has_feas_restore_phase", 1);
  paropt_options->setOption("tr_use_soc", 0);
  paropt_options->setOption("tr_max_iterations", max_iter);
  paropt_options->setOption("penalty_gamma", 50.0);
  paropt_options->setOption("qn_subspace_size", qn_subspace);
  paropt_options->setOption("qn_type", "bfgs");
  paropt_options->setOption("qn_diag_type", "yty_over_yts");
  paropt_options->setOption("abs_res_tol", 1e-8);
  paropt_options->setOption("starting_point_strategy", "affine_step");
  paropt_options->setOption("barrier_strategy", "mehrotra_predictor_corrector");
  paropt_options->setOption("tr_steering_barrier_strategy",
                            "mehrotra_predictor_corrector");
  paropt_options->setOption("tr_steering_starting_point_strategy",
                            "affine_step");
  paropt_options->setOption("use_line_search", 0);
  paropt_options->setOption("max_major_iters", 200);

  // Set up optimizer
  ParOptOptimizer *optimizer = new ParOptOptimizer(problem, paropt_options);
  optimizer->incref();

  // Optimize
  optimizer->optimize();

  // Free memory
  mat_props->decref();
  stiff_props->decref();
  geo->decref();
  mesh->decref();
  model->decref();
  topo->decref();
  forest->decref();
  bcs->decref();
  pcr->decref();
  problem->decref();
  paropt_options->decref();
  optimizer->decref();

  delete options;

  // This should be called last
  MPI_Finalize();

  return 0;
}