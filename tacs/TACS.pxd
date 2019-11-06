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

cdef extern from "TACSElementTypes.h":
    enum:
        TACS_OUTPUT_CONNECTIVITY
        TACS_OUTPUT_NODES
        TACS_OUTPUT_DISPLACEMENTS
        TACS_OUTPUT_STRAINS
        TACS_OUTPUT_STRESSES
        TACS_OUTPUT_EXTRAS

    enum ElementType:
        TACS_ELEMENT_NONE
        TACS_SCALAR_2D_ELEMENT
        TACS_SCALAR_3D_ELEMENT
        TACS_BEAM_OR_SHELL_ELEMENT
        TACS_PLANE_STRESS_ELEMENT
        TACS_SOLID_ELEMENT

    enum ElementLayout:
        TACS_LAYOUT_NONE
        TACS_POINT_ELEMENT
        TACS_LINE_ELEMENT
        TACS_LINE_QUADRATIC_ELEMENT
        TACS_LINE_CUBIC_ELEMENT
        TACS_TRI_ELEMENT
        TACS_TRI_QUADRATIC_ELEMENT
        TACS_TRI_CUBIC_ELEMENT
        TACS_QUAD_ELEMENT
        TACS_QUAD_QUADRATIC_ELEMENT
        TACS_QUAD_CUBIC_ELEMENT
        TACS_QUAD_QUARTIC_ELEMENT
        TACS_QUAD_QUINTIC_ELEMENT
        TACS_TETRA_ELEMENT
        TACS_TETRA_QUADRATIC_ELEMENT
        TACS_TETRA_CUBIC_ELEMENT
        TACS_HEXA_ELEMENT
        TACS_HEXA_QUADRATIC_ELEMENT
        TACS_HEXA_CUBIC_ELEMENT
        TACS_HEXA_QUARTIC_ELEMENT
        TACS_HEXA_QUINTIC_ELEMENT
        TACS_PENTA_ELEMENT
        TACS_PETTA_QUADRATIC_ELEMENT
        TACS_PENTA_CUBIC_ELEMENT

    enum ElementMatrixType:
        TACS_STIFFNESS_MATRIX
        TACS_MASS_MATRIX
        TACS_GEOMETRIC_STIFFNESS_MATRIX
        TACS_STIFFNESS_PRODUCT_DERIVATIVE

cdef extern from "KSM.h":
    enum MatrixOrientation:
        TACS_MAT_NORMAL
        TACS_MAT_TRANSPOSE

# Special functions required for converting pointers
cdef extern from "":
    TACSSchurMat* _dynamicSchurMat "dynamic_cast<TACSSchurMat*>"(TACSMat*)
    TACSSchurPc* _dynamicSchurPc "dynamic_cast<TACSSchurPc*>"(TACSPc*)
    TACSParallelMat* _dynamicParallelMat "dynamic_cast<TACSParallelMat*>"(TACSMat*)
    TACSMg* _dynamicTACSMg "dynamic_cast<TACSMg*>"(TACSPc*)
    GMRES* _dynamicGMRES "dynamic_cast<GMRES*>"(TACSKsm*)
    void deleteArray "delete []"(void*)

cdef extern from "TACSObject.h":
    cdef cppclass TACSObject:
        void incref()
        void decref()
        const char* getObjectName()

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

cdef extern from "TACSBVec.h":
    cdef cppclass TACSBVec(TACSVec):
        TACSBVec(TACSNodeMap*, int)
        int getSize(int*)
        int getBlockSize()
        TACSNodeMap *getNodeMap()
        int getArray(TacsScalar**)
        int readFromFile(const_char*)
        int writeToFile(const_char*)
        int getValues(int, const int*, TacsScalar*)
        void setValues(int, const int*, const TacsScalar*, TACSBVecOperation)
        void beginSetValues(TACSBVecOperation)
        void endSetValues(TACSBVecOperation)
        void beginDistributeValues()
        void endDistributeValues()

cdef extern from "TACSBVecDistribute.h":
    enum TACSBVecOperation:
        TACS_INSERT_VALUES
        TACS_ADD_VALUES
        TACS_INSERT_NONZERO_VALUES

    cdef cppclass TACSNodeMap(TACSObject):
        TACSNodeMap(MPI_Comm, int)
        MPI_Comm getMPIComm()
        void getOwnerRange(const int**)

    cdef cppclass TACSBVecIndices(TACSObject):
        TACSBVecIndices(int**, int)
        int getIndices(const int **)

    cdef cppclass TACSBVecDistribute(TACSObject):
        TACSBVecIndices *getIndices()

cdef class NodeMap:
    cdef TACSNodeMap *ptr

cdef inline _init_NodeMap(TACSNodeMap *ptr):
    vmap = NodeMap()
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

cdef extern from "TACSBVecInterp.h":
    cdef cppclass TACSBVecInterp(TACSObject):
        TACSBVecInterp(TACSNodeMap*, TACSNodeMap*, int)
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

cdef extern from "TACSParallelMat.h":
    cdef cppclass TACSParallelMat(TACSMat):
        pass

    cdef cppclass TACSAdditiveSchwarz(TACSPc):
        TACSAdditiveSchwarz(TACSParallelMat *mat, int levFill, double fill)

    cdef cppclass ApproximateSchur(TACSPc):
        TACSApproximateSchur(TACSParallelMat *mat, int levFill, double fill,
                             int inner_gmres_iters, double inner_rtol,
                             double inner_atol)

cdef extern from "TACSSchurMat.h":
    cdef cppclass TACSSchurMat(TACSMat):
        void getBCSRMat(BCSRMat**, BCSRMat**, BCSRMat**, BCSRMat**)
        TACSBVecDistribute *getLocalMap()
        TACSBVecDistribute *getSchurMap()

    cdef cppclass TACSSchurPc(TACSPc):
        TACSSchurPc(TACSSchurMat *mat, int levFill, double fill,
                    int reorder_schur_complement)
        void setMonitorFactorFlag(int)
        void setMonitorBackSolveFlag(int)

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

cdef extern from "TACSElementBasis.h":
    cdef cppclass TACSElementBasis(TACSObject):
        ElementLayout getLayoutType()
        int getNumNodes()
        int getNumParameters()

cdef class ElementBasis:
    cdef TACSElementBasis *ptr

cdef inline _init_ElementBasis(TACSElementBasis *ptr):
    basis = ElementBasis()
    basis.ptr = ptr
    basis.ptr.incref()
    return basis

cdef extern from "TACSElementModel.h":
    cdef cppclass TACSElementModel(TACSObject):
        int getSpatialDim()
        int getVarsPerNode()

cdef class ElementModel:
    cdef TACSElementModel *ptr

cdef inline _init_ElementModel(TACSElementModel *ptr):
    model = ElementModel()
    model.ptr = ptr
    model.ptr.incref()
    return model

cdef extern from "TACSElement.h":
    cdef cppclass TACSElement(TACSObject):
        void setComponentNum(int)
        int getComponentNum()
        int getVarsPerNode()
        int getNumNodes()
        int getNumVariables()
        int getDesignVarsPerNode()
        int getDesignVarNums(int, int, int*)
        int getDesignVars(int, int, TacsScalar*)
        int getDesignVarRange(int, int, TacsScalar*, TacsScalar*)
        TACSElementBasis* getElementBasis()
        TACSElementModel* getElementModel()

cdef class Element:
    cdef TACSElement *ptr

cdef inline _init_Element(TACSElement *ptr):
    elem = Element()
    elem.ptr = ptr
    elem.ptr.incref()
    return elem

cdef extern from "TACSFunction.h":
    cdef cppclass TACSFunction(TACSObject):
        void setDomain(int, int*)

cdef class Function:
    cdef TACSFunction *ptr

cdef extern from "TACSConstitutive.h":
    cdef cppclass TACSConstitutive(TACSObject):
        # TacsScalar getDVOutputValue(int, const double*)
        int getNumStresses()

cdef class Constitutive:
    cdef TACSConstitutive *ptr

cdef inline _init_Constitutive(TACSConstitutive *ptr):
    cons = Constitutive()
    cons.ptr = ptr
    cons.ptr.incref()
    return cons

cdef extern from "TACSAuxElements.h":
    cdef cppclass TACSAuxElements(TACSObject):
        TACSAuxElements(int)
        void addElement(int, TACSElement*)

cdef extern from "TACSAssembler.h":
    enum OrderingType"TACSAssembler::OrderingType":
        TACS_NATURAL_ORDER"TACSAssembler::NATURAL_ORDER"
        TACS_RCM_ORDER"TACSAssembler::RCM_ORDER"
        TACS_ND_ORDER"TACSAssembler::ND_ORDER"
        TACS_TACS_AMD_ORDER"TACSAssembler::TACS_AMD_ORDER"
        TACS_MULTICOLOR_ORDER"TACSAssembler::MULTICOLOR_ORDER"

    enum MatrixOrderingType"TACSAssembler::MatrixOrderingType":
        TACS_ADDITIVE_SCHWARZ"TACSAssembler::ADDITIVE_SCHWARZ"
        TACS_APPROXIMATE_SCHUR"TACSAssembler::APPROXIMATE_SCHUR"
        TACS_DIRECT_SCHUR"TACSAssembler::DIRECT_SCHUR"
        TACS_GAUSS_SEIDEL"TACSAssembler::GAUSS_SEIDEL"

    cdef cppclass TACSAssembler(TACSObject):
        TACSAssembler(MPI_Comm tacs_comm, int varsPerNode,
                      int numOwnedNodes, int numElements,
                      int numDependentNodes)

        # Set the element connectivity
        int setElementConnectivity(int *ptr, int *conn)
        int setElements(TACSElement **elements)
        int setDependentNodes(int *depNodeIndex,
                              int *depNodeToTacs,
                              double *depNodeWeights)
        void setDesignNodeMap(int _designVarsPerNode,
                              TACSNodeMap *_designVarMap)

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
        TACSNodeMap *getNodeMap()

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

        TACSBVec *createDesignVec()
        void getDesignVars(TACSBVec*)
        void setDesignVars(TACSBVec*)
        void getDesignVarRange(TACSBVec*, TACSBVec*)

        # Create vectors/matrices
        TACSBVec *createVec()
        TACSParallelMat *createMat()
        TACSSchurMat *createSchurMat(OrderingType)

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
        void evalFunctions(int numFuncs, TACSFunction **functions,
                           TacsScalar *funcVals)

        # Derivative evaluation routines
        void addDVSens(double coef, int numFuncs, TACSFunction **funcs,
                       TACSBVec **dfdx)
        void addSVSens(double alpha, double beta, double gamma,
                       int numFuncs, TACSFunction **funcs,
                       TACSBVec **fuSens)
        void addAdjointResProducts(double scale, int numAdjoints,
                                   TACSBVec **adjoint, TACSBVec **dfdx)
        void addXptSens(double coef, int numFuncs, TACSFunction **funcs,
                        TACSBVec **fXptSens)
        void addAdjointResXptSensProducts(double scale, int numAdjoints,
                                          TACSBVec **adjoint,
                                          TACSBVec **adjXptSens)

        # Add the derivative of the inner product with a matrix
        void addMatDVSensInnerProduct(double scale,
                                      ElementMatrixType matType,
                                      TACSBVec *psi, TACSBVec *phi,
                                      TACSBVec *dfdx)
        void evalMatSVSensInnerProduct(ElementMatrixType matType,
                                       TACSBVec *psi, TACSBVec *phi,
                                       TACSBVec *res)

        # Test routines
        void testElement(int elemNum, int print_level, double dh,
                         double rtol, double atol)
        void testFunction(TACSFunction *func, double dh)

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
        int getNumNodes()
        int getNumElements()
        TACSAssembler* createTACS(int vars_per_node,
                                  OrderingType order_type,
                                  MatrixOrderingType mat_type)
        void addAuxElement(TACSAuxElements *aux, int comp_num,
	                         TACSElement *_element)
        void addFunctionDomain(TACSFunction * function,
                               int num_comps, int *comp_nums)
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
        void setElements(int _num_elems, TACSElement **_elements)
        void setNodes(TacsScalar *_Xpts)
        void setReorderingType(OrderingType _order_type,
                               MatrixOrderingType _mat_type)
        void partitionMesh(int split_size, int *part)
        int getElementPartition(const int **)
        TACSAssembler *createTACS()
        int getNodeNums(const int**)
        void getAssemblerNodeNums(TACSAssembler*, int, const int*,
                                  int*, int**)

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
        void setUseSchurMat(int, OrderingType)
        void setInitNewtonDeltaFraction(double)
        void setKrylovSubspaceMethod(TACSKsm *_ksm)
        void setTimeInterval(double, double)
        void setFunctions(int num_funcs, TACSFunction **funcs,
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
        void getGradient(int func_num, TACSBVec **dfdx)
        void getXptGradient(int func_num, TACSBVec **dfdXpt)
        double getStates(int step_num,
                         TACSBVec **q, TACSBVec **qdot, TACSBVec **qddot)

        # Configure output
        void setOutputPrefix(const_char *prefix)
        void setOutputFrequency(int write_freq)
        void setFH5(TACSToFH5 *_f5)
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
