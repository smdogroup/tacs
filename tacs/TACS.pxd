#  This file is part of TACS: The Toolkit for the Analysis of Composite
#  Structures, a parallel finite-element code for structural and
#  multidisciplinary design optimization.
#
#  Copyright (C) 2014 Georgia Tech Research Corporation
#
#  TACS is licensed under the Apache License, Version 2.0 (the
#  "License"); you may not use this software except in compliance with
#  the License.  You may obtain a copy of the License at
#
#  http://www.apache.org/licenses/LICENSE-2.0

# distutils: language=c++

# For MPI capabilities
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy
from libc.string cimport const_char

# Import the major python version
from cpython.version cimport PY_MAJOR_VERSION

# Import numpy
cimport numpy as np
import numpy as np

# Typdefs required for either real or complex mode
include "TacsTypedefs.pxi"

cdef inline char* convert_to_chars(s):
    if isinstance(s, unicode):
        s = (<unicode>s).encode('utf8')
    return s

cdef inline convert_bytes_to_str(bytes s):
    if PY_MAJOR_VERSION >= 3:
        return s.decode('utf8')
    return s

cdef extern from "TACSElement.h":
    enum ElementType:
        TACS_ELEMENT_NONE
        TACS_POINT_ELEMENT
        TACS_EULER_BEAM
        TACS_TIMOSHENKO_BEAM
        TACS_PLANE_STRESS
        TACS_SHELL
        TACS_SOLID
        TACS_Q3D_ELEMENT
        TACS_RIGID
        TACS_POISSON_2D_ELEMENT
        TACS_POISSON_3D_ELEMENT

    enum ElementMatrixType:
        STIFFNESS_MATRIX
        MASS_MATRIX
        GEOMETRIC_STIFFNESS_MATRIX
        STIFFNESS_PRODUCT_DERIVATIVE

cdef extern from "KSM.h":
    enum MatrixOrientation:
        NORMAL
        TRANSPOSE

cdef extern from "BVecDist.h":
    enum TACSBVecOperation:
        TACS_INSERT_VALUES
        TACS_ADD_VALUES
        TACS_INSERT_NONZERO_VALUES

# Special functions required for converting pointers
cdef extern from "":
    ScMat* _dynamicScMat "dynamic_cast<ScMat*>"(TACSMat*)
    TACSPMat* _dynamicPMat "dynamic_cast<TACSPMat*>"(TACSMat*)
    PcScMat* _dynamicPcScMat "dynamic_cast<PcScMat*>"(TACSPc*)
    TACSMg* _dynamicTACSMg "dynamic_cast<TACSMg*>"(TACSPc*)
    GMRES* _dynamicGMRES "dynamic_cast<GMRES*>"(TACSKsm*)
    void deleteArray "delete []"(void*)

cdef extern from "TACSObject.h":
    cdef cppclass TACSObject:
        void incref()
        void decref()

    cdef MPI_Datatype TACS_MPI_TYPE

cdef extern from "KSM.h":
    cdef cppclass TACSVec(TACSObject):
        TacsScalar norm()
        void scale(TacsScalar alpha)
        TacsScalar dot(TACSVec *x)
        void axpy(TacsScalar alpha, TACSVec *x)
        void copyValues(TACSVec *x)
        void axpby(TacsScalar alpha, TacsScalar beta, TACSVec *x)
        void zeroEntries()
        void setRand(double lower, double upper)

    cdef cppclass TACSMat(TACSObject):
        TACSVec *createVec()
        void zeroEntries()
        void mult(TACSVec *x, TACSVec *y)
        void copyValues(TACSMat *mat)
        void scale(TacsScalar alpha)
        void axpy(TacsScalar alpha, TACSMat *mat)

    cdef cppclass TACSPc(TACSObject):
        void factor()
        void applyFactor(TACSVec *x, TACSVec *y)
        void getMat(TACSMat**)

    cdef cppclass KSMPrint(TACSObject):
        pass

    cdef cppclass KSMPrintStdout(KSMPrint):
        KSMPrintStdout(char *descript, int rank, int freq)

    cdef cppclass TACSKsm(TACSObject):
        TACSVec *createVec()
        void setOperators(TACSMat *_mat, TACSPc *_pc)
        void getOperators(TACSMat **_mat, TACSPc **_pc)
        void solve(TACSVec *b, TACSVec *x, int zero_guess)
        void setTolerances(double _rtol, double _atol)
        void setMonitor(KSMPrint *_monitor)

    cdef cppclass GMRES(TACSKsm):
        GMRES(TACSMat *_mat, TACSPc *_pc, int _m,
              int _nrestart, int _isFlexible )
        void setTimeMonitor()

cdef extern from "BVec.h":
    cdef cppclass TACSBVec(TACSVec):
        TACSBVec(TACSVarMap*, int)
        int getSize(int*)
        int getBlockSize()
        TACSVarMap *getVarMap()
        int getArray(TacsScalar**)
        int readFromFile(const_char*)
        int writeToFile(const_char*)
        int getValues(int, const int*, TacsScalar*)
        void setValues(int, const int*, const TacsScalar*, TACSBVecOperation)
        void beginSetValues(TACSBVecOperation)
        void endSetValues(TACSBVecOperation)
        void beginDistributeValues()
        void endDistributeValues()

cdef extern from "BVecDist.h":
    cdef cppclass TACSVarMap(TACSObject):
        TACSVarMap(MPI_Comm, int)
        MPI_Comm getMPIComm()
        void getOwnerRange(const int**)

    cdef cppclass TACSBVecIndices(TACSObject):
        TACSBVecIndices(int**, int)
        int getIndices(const int **)

    cdef cppclass TACSBVecDistribute(TACSObject):
        TACSBVecIndices *getIndices()

cdef class VarMap:
    cdef TACSVarMap *ptr

cdef inline _init_VarMap(TACSVarMap *ptr):
    vmap = VarMap()
    vmap.ptr = ptr
    vmap.ptr.incref()
    return vmap

cdef class VecIndices:
    cdef TACSBVecIndices *ptr

cdef inline _init_VecIndices(TACSBVecIndices *ptr):
    indices = VecIndices()
    indices.ptr = ptr
    indices.ptr.incref()
    return indices

cdef extern from "BVecInterp.h":
    cdef cppclass TACSBVecInterp(TACSObject):
        TACSBVecInterp(TACSVarMap*, TACSVarMap*, int)
        TACSBVecInterp(TACSAssembler*, TACSAssembler*)
        void mult(TACSBVec*, TACSBVec*)
        void multAdd(TACSBVec*, TACSBVec*, TACSBVec*)
        void multTranspose(TACSBVec*, TACSBVec*)
        void multTransposeAdd(TACSBVec*, TACSBVec*, TACSBVec*)
        void multWeightTranspose(TACSBVec*, TACSBVec*)
        void initialize()

cdef class VecInterp:
    cdef TACSBVecInterp *ptr

cdef inline _init_VecInterp(TACSBVecInterp *ptr):
    interp = VecInterp()
    interp.ptr = ptr
    interp.ptr.incref()
    return interp

cdef extern from "BCSRMat.h":
    cdef cppclass BCSRMat(TACSObject):
        int getBlockSize()
        int getRowDim()
        int getColDim()
        void getDenseColumnMajor(TacsScalar*)

cdef extern from "PMat.h":
    cdef cppclass TACSPMat(TACSMat):
        pass

    cdef cppclass TACSAdditiveSchwarz(TACSPc):
        TACSAdditiveSchwarz(TACSPMat *mat, int levFill, double fill)

    cdef cppclass ApproximateSchur(TACSPc):
        TACSApproximateSchur(TACSPMat *mat, int levFill, double fill,
                             int inner_gmres_iters, double inner_rtol,
                             double inner_atol)

cdef extern from "DistMat.h":
    cdef cppclass TACSDistMat(TACSPMat):
        pass

cdef extern from "ScMat.h":
    cdef cppclass ScMat(TACSMat):
        void getBCSRMat(BCSRMat**, BCSRMat**, BCSRMat**, BCSRMat**)
        TACSBVecDistribute *getLocalMap()
        TACSBVecDistribute *getSchurMap()

    cdef cppclass PcScMat(TACSPc):
        PcScMat(ScMat *mat, int levFill, double fill,
                int reorder_schur_complement)
        void setMonitorFactorFlag(int)
        void setMonitorBackSolveFlag(int)

cdef extern from "FEMat.h":
    cdef cppclass FEMat(ScMat):
        pass

cdef extern from "TACSMg.h":
    cdef cppclass TACSMg(TACSPc):
        TACSMg(MPI_Comm, int, double, int, int)
        void setLevel(int, TACSAssembler*, TACSBVecInterp*, int,
                      TACSMat*, TACSPc*)
        void setVariables(TACSBVec*)
        void assembleJacobian(double, double, double, TACSBVec*,
                              MatrixOrientation)
        void assembleMatType(ElementMatrixType,
                             MatrixOrientation)
        void setMonitor(KSMPrint*)

cdef extern from "TACSElement.h":
    void TACSSetElementFDStepSize"TACSElement::setStepSize"(double)

    cdef cppclass TACSElement(TACSObject):
        int numNodes()
        int numVariables()
        int getComponentNum()
        void setComponentNum(int)
        TACSConstitutive *getConstitutive()

cdef extern from "TACSFunction.h":
    cdef cppclass TACSFunction(TACSObject):
        void setDomain(int*, int)

cdef extern from "TACSConstitutive.h":
    cdef cppclass TACSConstitutive(TACSObject):
        TacsScalar getDVOutputValue(int, const double*)

cdef class Element:
    cdef TACSElement *ptr

cdef inline _init_Element(TACSElement *ptr):
    elem = Element()
    elem.ptr = ptr
    elem.ptr.incref()
    return elem

cdef class Function:
    cdef TACSFunction *ptr

cdef class Constitutive:
    cdef TACSConstitutive *ptr

cdef inline _init_Constitutive(TACSConstitutive *ptr):
    cons = Constitutive()
    cons.ptr = ptr
    cons.ptr.incref()
    return cons

cdef class Vec:
    cdef TACSBVec *ptr

cdef inline _init_Vec(TACSBVec *ptr):
    vec = Vec()
    vec.ptr = ptr
    vec.ptr.incref()
    return vec

cdef class Mat:
    cdef TACSMat *ptr

cdef inline _init_Mat(TACSMat *ptr):
    mat = Mat()
    mat.ptr = ptr
    mat.ptr.incref()
    return mat

cdef class Pc:
    cdef TACSPc *ptr

cdef inline _init_Pc(TACSPc *ptr):
    pc = Pc()
    pc.ptr = ptr
    pc.ptr.incref()
    return pc

cdef class Mg(Pc):
    cdef TACSMg *mg

cdef inline _init_Mg(TACSMg *ptr):
    mg = Mg()
    mg.mg = ptr
    mg.ptr = ptr
    mg.mg.incref()
    return mg

cdef class KSM:
    cdef TACSKsm *ptr

cdef extern from "TACSAuxElements.h":
    cdef cppclass TACSAuxElements(TACSObject):
        TACSAuxElements(int)
        void addElement(int, TACSElement*)

cdef extern from "TACSAssembler.h":
    enum OrderingType"TACSAssembler::OrderingType":
        NATURAL_ORDER"TACSAssembler::NATURAL_ORDER"
        RCM_ORDER"TACSAssembler::RCM_ORDER"
        AMD_ORDER"TACSAssembler::AMD_ORDER"
        ND_ORDER"TACSAssembler::ND_ORDER"
        TACS_AMD_ORDER"TACSAssembler::TACS_AMD_ORDER"
        MULTICOLOR_ORDER"TACSAssembler::MULTICOLOR_ORDER"

    enum MatrixOrderingType"TACSAssembler::MatrixOrderingType":
        ADDITIVE_SCHWARZ"TACSAssembler::ADDITIVE_SCHWARZ"
        APPROXIMATE_SCHUR"TACSAssembler::APPROXIMATE_SCHUR"
        DIRECT_SCHUR"TACSAssembler::DIRECT_SCHUR"
        GAUSS_SEIDEL"TACSAssembler::GAUSS_SEIDEL"

    cdef cppclass TACSAssembler(TACSObject):
        TACSAssembler(MPI_Comm tacs_comm, int varsPerNode,
                      int numOwnedNodes, int numElements,
                      int numDependentNodes)

        # Set the element connectivity
        int setElementConnectivity(int *conn, int *ptr)
        int setElements(TACSElement **elements)
        int setDependentNodes(int *depNodeIndex,
                              int *depNodeToTacs,
                              double *depNodeWeights)

        # Add boundary conditions
        void addBCs(int nnodes, int *nodes,
                    int nbcs, int *vars, TacsScalar *vals)
        void addInitBCs(int nnodes, int *nodes,
                        int nbcs, int *vars, TacsScalar *vals)

        void computeReordering(OrderingType, MatrixOrderingType)

        # Finalize the mesh - no further elements or nodes may be added
        # following this call
        void initialize()

        # Return information about the TACSObject
        int getVarsPerNode()
        int getNumNodes()
        int getNumDependentNodes()
        int getNumOwnedNodes()
        int getNumElements()
        TACSVarMap *getVarMap()

        # Return information about the element
        TACSElement **getElements()
        TACSElement *getElement(int, TacsScalar*, TacsScalar*,
                                TacsScalar*, TacsScalar*)
        TACSElement *getElement(int, const int**, int*)

        # MPI communicator
        MPI_Comm getMPIComm()

        # Set the auxiliary element class
        void setAuxElements(TACSAuxElements*)
        TACSAuxElements *getAuxElements()

        TACSBVec *createNodeVec()
        void setNodes(TACSBVec*)
        void getNodes(TACSBVec*)
        void getDesignVars(TacsScalar*, int)
        void setDesignVars(TacsScalar*, int)
        void getDesignVarRange(TacsScalar*, TacsScalar *, int)

        # Create vectors/matrices
        TACSBVec *createVec()
        TACSDistMat *createMat()
        FEMat *createFEMat(OrderingType)

        # Reorder the vector based on the selected reordering
        void getReordering(int*)
        void reorderVec(TACSBVec*)

        # Set/get the simulation time
        void setSimulationTime(double)
        double getSimulationTime()

        void applyBCs(TACSVec*)
        void applyBCs(TACSMat*)
        void setBCs(TACSVec*)

        # Zero the variables
        void zeroVariables()
        void zeroDotVariables()
        void zeroDDotVariables()

        # Set and retrieve the state variables
        void setVariables(TACSBVec*, TACSBVec*, TACSBVec*)
        void getVariables(TACSBVec*, TACSBVec*, TACSBVec*)
        void copyVariables(TACSBVec*, TACSBVec*, TACSBVec*)

        # Get the initial conditions
        void getInitConditions(TACSBVec*, TACSBVec*, TACSBVec*)

        # Evaluate the kinetic/potential energy
        void evalEnergies(TacsScalar*, TacsScalar*)

        # Assembly routines
        void assembleRes(TACSBVec *residual)
        void assembleJacobian(double alpha, double beta, double gamma,
                              TACSBVec *residual, TACSMat *A,
                              MatrixOrientation matOr)
        void assembleMatType(ElementMatrixType matType,
                             TACSMat *A, MatrixOrientation matOr)
        void addJacobianVecProduct(TacsScalar scale,
                                   double alpha, double beta, double gamma,
                                   TACSBVec *x, TACSBVec *y,
                                   MatrixOrientation matOr)

        # Evaluation routines
        void evalFunctions(TACSFunction **functions, int numFuncs,
                           TacsScalar *funcVals)

        # Derivative evaluation routines
        void addDVSens(double coef, TACSFunction **funcs, int numFuncs,
                       TacsScalar *fdvSens, int numDVs)
        void addSVSens(double alpha, double beta, double gamma,
                       TACSFunction **funcs, int numFuncs,
                       TACSBVec **fuSens)
        void addAdjointResProducts(double scale,
                                   TACSBVec **adjoint, int numAdjoints,
                                   TacsScalar *dvSens, int numDVs)
        void addXptSens(double coef, TACSFunction **funcs, int numFuncs,
                        TACSBVec **fXptSens)
        void addAdjointResXptSensProducts(double scale,
                                          TACSBVec **adjoint, int numAdjoints,
                                          TACSBVec **adjXptSens)

        # Add the derivative of the inner product with a matrix
        void addMatDVSensInnerProduct(double scale,
                                      ElementMatrixType matType,
                                      TACSBVec *psi, TACSBVec *phi,
                                      TacsScalar *dvSens, int numDVs)
        void evalMatSVSensInnerProduct(ElementMatrixType matType,
                                       TACSBVec *psi, TACSBVec *phi,
                                       TACSBVec *res)

        # Test routines
        void testElement(int elemNum, int print_level, double dh,
                         double rtol, double atol)
        void testConstitutive(int elemNum, int print_level)
        void testFunction(TACSFunction *func, int num_dvs,
                          double dh)

        # Set the number of threads
        void setNumThreads(int t)

cdef class Assembler:
    cdef TACSAssembler *ptr

cdef inline _init_Assembler(TACSAssembler *ptr):
    tacs = Assembler()
    tacs.ptr = ptr
    tacs.ptr.incref()
    return tacs

cdef extern from "JacobiDavidson.h":
   # Set the type of recycling scheme
   enum JDRecycleType:
      JD_SUM_TWO
      JD_NUM_RECYCLE

cdef extern from "TACSBuckling.h":
    cdef cppclass TACSFrequencyAnalysis(TACSObject):
        TACSFrequencyAnalysis(TACSAssembler *, TacsScalar,
                              TACSMat*, TACSMat*, TACSKsm*,
                              int, int, double)
        TACSFrequencyAnalysis(TACSAssembler*, TacsScalar,
                              TACSMat*, TACSMat*, TACSMat*,
                              TACSPc*,
                              int, int, int, double, double, double,
                              int, JDRecycleType)
        TACSAssembler* getTACS()
        TacsScalar getSigma()
        void setSigma(TacsScalar)
        void solve(KSMPrint*)
        TacsScalar extractEigenvalue(int, TacsScalar*)
        TacsScalar extractEigenvector(int, TACSBVec*, TacsScalar*)

    cdef cppclass TACSLinearBuckling(TACSObject):
        TACSLinearBuckling( TACSAssembler *,
                            TacsScalar,
                            TACSMat *, TACSMat *,
                            TACSMat *, TACSKsm *,
                            int, int, double)
        TACSAssembler* getTACS()
        TacsScalar getSigma()
        void setSigma(TacsScalar)
        void solve(TACSVec*, KSMPrint*)
        void evalEigenDVSens(int, TacsScalar, int)
        TacsScalar extractEigenvalue(int, TacsScalar*)
        TacsScalar extractEigenvector(int, TACSBVec*, TacsScalar*)

cdef extern from "TACSMeshLoader.h":
    cdef cppclass TACSMeshLoader(TACSObject):
        TACSMeshLoader(MPI_Comm _comm)
        int scanBDFFile(char *file_name)
        int getNumComponents()
        const char *getComponentDescript(int comp_num)
        const char *getElementDescript(int comp_num)
        void setElement(int component_num, TACSElement *_element)
        void setConvertToCoordinate(int flag)
        int getNumNodes()
        int getNumElements()
        TACSAssembler*createTACS(int vars_per_node,
                                  OrderingType order_type,
                                  MatrixOrderingType mat_type)
        void addAuxElement(TACSAuxElements *aux, int comp_num,
	                   TACSElement *_element)

        void addFunctionDomain(TACSFunction * function,
                               int *comp_nums, int num_comps)
        void getConnectivity(int *_num_nodes, int *_num_elements,
                             const int **_elem_node_ptr,
                             const int **_elem_node_conn,
                             const int **elem_compoennts,
                             TacsScalar**_Xpts)
        void getBCs(int *_num_bcs, const int **_bc_nodes,
                    const int **_bc_vars, const int **_bc_ptr,
                    const TacsScalar **_bc_vals)

cdef extern from "TACSCreator.h":
    cdef cppclass TACSCreator(TACSObject):
        TACSCreator(MPI_Comm comm, int _vars_per_node)
        void setGlobalConnectivity(int _num_nodes, int _num_elements,
                                   int *_elem_node_ptr,
                                   int *_elem_node_conn,
                                   int *_elem_id_nums )
        void setBoundaryConditions(int _num_bcs, int *_bc_nodes,
                                   int *_bc_vars, int *_bc_ptr)
        void setDependentNodes(int num_dep_nodes,
                               int *_dep_node_ptr,
                               int *_dep_node_conn,
                               double *_dep_node_weights )
        void setElements(TACSElement **_elements, int _num_elems)
        void setNodes(TacsScalar *_Xpts)
        void setReorderingType(OrderingType _order_type,
                               MatrixOrderingType _mat_type)
        void partitionMesh(int split_size, int *part)
        int getElementPartition(const int **)
        TACSAssembler *createTACS()
        int getNodeNums(const int**)
        void getTacsNodeNums(TACSAssembler*, const int*, int, int**, int*)

cdef extern from "TACSToFH5.h":
    cdef cppclass TACSToFH5(TACSObject):
        TACSToFH5(TACSAssembler *_tacs, ElementType _elem_type, int _out_type)
        void setComponentName(int comp_num, char *group_name)
        void writeToFile(char *filename)

cdef extern from "TACSIntegrator.h":
    # Declare the TACSIntegrator base class
    cdef cppclass TACSIntegrator(TACSObject):
        # Setters for class variables
        void setRelTol(double)
        void setAbsTol(double)
        void setMaxNewtonIters(int)
        void setPrintLevel(int level, const_char *filename)
        void setJacAssemblyFreq(int)
        void setUseLapack(int)
        void setUseFEMat(int,OrderingType)
        void setInitNewtonDeltaFraction(double)
        void setKrylovSubspaceMethod(TACSKsm *_ksm)
        void setTimeInterval(double, double)
        void setFunctions(TACSFunction **funcs, int num_funcs,
                          int num_design_vars,
                          int start_step, int end_step)
        void lapackNaturalFrequencies(int, TACSBVec*, TACSBVec*,
                                      TACSBVec*, TacsScalar*, TacsScalar*)

        # Forward mode functions
        int iterate(int step_num,TACSBVec *forces)
        int integrate()
        void evalFunctions(TacsScalar *fvals)

        # Reverse mode functions
        void iterateAdjoint(int step_num, TACSBVec **adj_rhs)
        void initAdjoint(int step_num)
        void integrateAdjoint()
        void postAdjoint(int step_num)
        void getAdjoint(int step_num, int func_num, TACSBVec **adjoint)
        void getGradient(TacsScalar *_dfdx)
        void getXptGradient(int func_num, TACSBVec **dfdXpt)

        double getStates(int step_num,
                         TACSBVec **q, TACSBVec **qdot, TACSBVec **qddot)

        # Configure output
        void setOutputPrefix(const_char *prefix)
        void setOutputFrequency(int write_freq)
        void setRigidOutput(TACSToFH5 *_rigidf5)
        void setShellOutput(TACSToFH5 *_shellf5)
        void setBeamOutput(TACSToFH5 *_beamf5)
        void setSolidOutput(TACSToFH5 *_solidf5)
        void writeSolution(const_char *filename, int format)
        void writeSolutionToF5(int step_num)

        int getNumTimeSteps()
        void writeRawSolution(const_char *name, int format_flag)

        # Debug adjoint
        void checkGradients(double dh)

    # BDF Implementation of the integrator
    cdef cppclass TACSBDFIntegrator(TACSIntegrator):
        TACSBDFIntegrator(TACSAssembler *tacs,
                          double tinit, double tfinal,
                          double num_steps,
                          int max_bdf_order)

    # DIRK Implementation of the integrator
    cdef cppclass TACSDIRKIntegrator(TACSIntegrator):
        TACSDIRKIntegrator(TACSAssembler *tacs,
                           double tinit, double tfinal,
                           double num_steps,
                           int stages)

    # ABM Implementation of the integrator
    cdef cppclass TACSABMIntegrator(TACSIntegrator):
        TACSABMIntegrator(TACSAssembler *tacs,
                          double tinit, double tfinal,
                          double num_steps,
                          int max_abm_order)

    # NBG Implementation of the integrator
    cdef cppclass TACSNBGIntegrator(TACSIntegrator):
        TACSNBGIntegrator(TACSAssembler *tacs,
                          double tinit, double tfinal,
                          double num_steps,
                          int order)
